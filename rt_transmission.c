/* This file is part of Exo_Transmit.

    Exo_Transmit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Exo_Transmit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Exo_Transmit.  If not, see <http://www.gnu.org/licenses/>.
*/

/*--------------------- rt_transmission.c ------------------------

Author: Eliza Kempton (kemptone@grinnell.edu)

------------------------------------------------------------------ */

/* Calculate a transmission spectrum with no scattering.  

------------------------------------------------------------------ */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "opac.h"
#include "atmos.h" 
#include "constant.h"
#include "include.h"
#include "nrutil.h"
#include "vars.h"
#include "prototypes.h"

/* --- Global variables ------------------------------------------ */

extern struct Atmos atmos;
extern struct Opac opac;

/* ------- begin ---------------- RT_Transmit -------------------- */


int RT_Transmit()
{
  char **fileArray = getFileArray();
  vars variables = getVars();
  getNTau(&variables, fileArray[0]);
  int NTAU = variables.NTAU;
  double R_PLANET = variables.R_PLANET;
  double R_STAR = variables.R_STAR;
  int NLAMBDA = variables.NLAMBDA;
  int NTEMP = variables.NTEMP;
  int NPRESSURE = variables.NPRESSURE;
  double G = variables.G;

  double **kappa_nu, **tau_nu, **dtau_nu, **tau_tr;
  double *intensity, *flux_st, *flux_pl, *flux_tr;
  double *ds, *theta, *dtheta, R, t_star;
  int i, j, a, b;
  FILE *file;
  
  /*   Allocate memory */
  
  kappa_nu = dmatrix(0, NLAMBDA-1, 0, NTAU-1);
  dtau_nu = dmatrix(0, NLAMBDA-1, 0, NTAU-1);
  tau_nu = dmatrix(0, NLAMBDA-1, 0, NTAU-1);
  tau_tr = dmatrix(0, NLAMBDA-1, 0, NTAU-1);

  intensity = dvector(0, NLAMBDA-1);
  flux_st = dvector(0, NLAMBDA-1);
  flux_pl = dvector(0, NLAMBDA-1);
  flux_tr = dvector(0, NLAMBDA-1);
  ds = dvector(0, NTAU-1);
  theta = dvector(0, NTAU-1);
  dtheta = dvector(0, NTAU-1);

  /*  Populate tau_tr with zeros */
  
  for (i=0; i<NLAMBDA; i++)
    for (j=0; j<NTAU; j++)
      tau_tr[i][j] = 0.;
  
  /*   Fill in kappa_nu, dtau_nu,  and ds */
  
  R = 0.0;

  // Grab the layer ~ 1 mbar --- Teal
  int i_mbar;
  Locate(NTAU, atmos.P, 1e-3, &i_mbar);
  
  for(j=0; j<NTAU; j++){
    
    if(j==0){
      ds[j] = atmos.P[j]* (KBOLTZMANN * atmos.T[j]) /
	(atmos.mu[j] * AMU * atmos.P[j]*G);
    } 
    else{
      ds[j] = (atmos.P[j] - atmos.P[j-1])* (KBOLTZMANN * atmos.T[j]) /
	(atmos.mu[j] * AMU * atmos.P[j]*G);
    }

    R += ds[j];
    for(i=0; i<NLAMBDA; i++){
      dtau_nu[i][j] = atmos.kappa_nu[i][j] * ds[j];
    }
    
    // Write file with opacity at P ~ 1 mbar --- Teal
    if (j == i_mbar) {
        // Write out a file with wavelength-dependent opacities at this layer
        FILE *test_file;
        test_file = fopen("test_opacities.dat", "w");

        fprintf(test_file, "Wavelength (m)\t Kappa_nu (m^2/kg)\n");

        for (i=0; i<NLAMBDA; i++)
            fprintf(test_file, "%e\t%e\n", atmos.lambda[i], atmos.kappa_nu[i][j]);

        fclose(test_file);
    };
  }
  
  /*   Allocate memory and fill in tau_nu */
  
  for(i=0; i<NLAMBDA; i++){
    for(j=0; j<NTAU; j++){
      if(j==0)
	tau_nu[i][j] = dtau_nu[i][j];
      else
	tau_nu[i][j] = tau_nu[i][j-1] + dtau_nu[i][j];
    }
  }
  
  
  /*   Calculate line-of-sight optical depths along NTAU rays.  Uses the
       routines found in geometry.c */
  Tau_LOS(atmos.kappa_nu, tau_tr, ds, NTAU, R_PLANET, NLAMBDA);
  
  Angles(ds, theta, dtheta, NTAU, R_PLANET, R_STAR);
  
  /*   Calculate the incident and emergent intensities.  Print to Output
       file specified in userInput.in  */
  
  R -=  Radius(R_PLANET, ds, NTAU);
  printf("R %f\n", (1.0 -  SQ(R)/SQ(R_STAR)));
  
  file = fopen(fileArray[2], "w");
  
  t_star = 6000.0;  // Arbitrary stellar temperature.  Gets divided out.

  fprintf(file, "Wavelength\tTransit Depth\n");
  fprintf(file, "(meters)\t(percent)\n");
  

  for(i=0; i<NLAMBDA; i++){
    flux_pl[i] = 0.0;
    intensity[i] = Planck(t_star, atmos.lambda[i]);
    
    for(j=0; j<NTAU; j++){
      flux_pl[i] += intensity[i] * exp(-tau_tr[i][j]) *
	cos(theta[j]) * sin(theta[j]) * dtheta[j];
    }
    
    flux_pl[i] *= 2 * PI;
    
    flux_st[i] = PI * intensity[i];
    flux_tr[i] = (1.0 -  SQ(Radius(R_PLANET, ds, NTAU))/SQ(R_STAR)) 
      * flux_st[i] + flux_pl[i];
    
    fprintf(file, "%e\t%e\n", atmos.lambda[i], 100.0*(1.0-flux_tr[i]/flux_st[i]));
  }
  
  fclose(file);
  
  /*   Free memory */
  
  free_dmatrix(kappa_nu, 0, NLAMBDA-1, 0, NTAU-1);
  free_dmatrix(dtau_nu, 0, NLAMBDA-1, 0, NTAU-1);
  free_dmatrix(tau_nu, 0, NLAMBDA-1, 0, NTAU-1);
  free_dmatrix(tau_tr, 0, NLAMBDA-1, 0, NTAU-1);
  free_dvector(intensity, 0, NLAMBDA-1);
  free_dvector(flux_st, 0, NLAMBDA-1);
  free_dvector(flux_pl, 0, NLAMBDA-1);
  free_dvector(flux_tr, 0, NLAMBDA-1);
  free_dvector(ds, 0, NTAU-1);
  free_dvector(theta, 0, NTAU-1);
  free_dvector(dtheta, 0, NTAU-1);
  
  return 0;
}

/* ------- end ------------------ RT_Transmit -------------------- */


