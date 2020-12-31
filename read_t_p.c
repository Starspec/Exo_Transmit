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


/*------------ file ------- read_t_p.c ---------------------------

Author: Eliza Kempton (kemptone@grinnell.edu)

------------------------------------------------------------------ */

/* Reads in the temperature - pressure profile from the file 
   defined in userIinput.in
------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atmos.h"
#include "nrutil.h"
#include "opac.h"
#include "prototypes.h"
#include "vars.h"

/* --- Global variables ------------------------------------------ */

extern struct Atmos atmos;
extern struct Chem chem;
extern struct Opac opac;

/* ------- begin --------------------- ReadTP.c ------------------ */

void ReadTP()
{
  
  /* Get relevant variables */
  char **fileArray = getFileArray();	
  vars variables = getVars();
  getNTau(&variables, fileArray[0]);
  
  /* Rename for convenience */
  int NTAU = variables.NTAU;
  int THRESHOLD = variables.THRESHOLD;
  
  int i;
  char dum[9];
  FILE *file;

  /* Molecular masses */
  /* Masses sourced from http://www.webqc.org/mmcalc.php */
  
  double m_H2 = 2.0158;
  double m_H = 1.0079;
  double m_H2O = 18.0152;
  double m_CH4 = 16.0423;
  double m_CO = 28.010;
  double m_CO2 = 44.010;
  double m_O = 15.9994;
  double m_C = 12.0107;
  double m_N = 14.0067;
  double m_NH3 = 17.031;
  double m_N2 = 28.0134;
  double m_O2 = 31.9988;
  double m_O3 = 47.9982;
  double m_C2H2 = 26.0373;
  double m_C2H4 = 28.0532;
  double m_C2H6 = 30.0690;
  double m_H2CO = 30.0260;
  double m_H2S = 34.0809;
  double m_HCN = 27.0253;
  double m_NO = 30.0061;
  double m_NO2 = 46.00550;
  double m_OCS = 60.0751;
  double m_OH = 17.0073;
  double m_SO2 = 64.0638;
  
  /* Allocate memory for atmos structure */
  
  atmos.P = dvector(0, NTAU-1);
  atmos.T = dvector(0, NTAU-1);
  atmos.mu = dvector(0, NTAU-1);

  atmos.H = dvector(0, NTAU-1);
  atmos.H2O = dvector(0, NTAU-1);
  atmos.OH = dvector(0, NTAU-1);
  atmos.O = dvector(0, NTAU-1);
  atmos.O2 = dvector(0, NTAU-1);
  atmos.O3 = dvector(0, NTAU-1);
  atmos.CO = dvector(0, NTAU-1);
  atmos.CO2 = dvector(0, NTAU-1);
  atmos.C = dvector(0, NTAU-1);
  atmos.CH4 = dvector(0, NTAU-1);
  atmos.C2H2 = dvector(0, NTAU-1);
  atmos.H2CO = dvector(0, NTAU-1);
  atmos.N = dvector(0, NTAU-1);
  atmos.N2 = dvector(0, NTAU-1);
  atmos.NO = dvector(0, NTAU-1);
  atmos.NO2 = dvector(0, NTAU-1);
  atmos.NH3 = dvector(0, NTAU-1);
  atmos.HCN = dvector(0, NTAU-1);
  atmos.H2S = dvector(0, NTAU-1);
  atmos.SO2 = dvector(0, NTAU-1);
  atmos.OCS = dvector(0, NTAU-1);
  atmos.C2H4 = dvector(0, NTAU-1);
  atmos.C2H6 = dvector(0, NTAU-1);
  atmos.H2 = dvector(0, NTAU-1);
  
  
  file = fopen(fileArray[0], "r");						
  if(file == NULL){
    printf("\nread_t_p.c:\nError opening file: No such file or directory\n\n");
    exit(1);
  }
  
  /* Read in T-P profile */

  fscanf(file, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s", dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum);
  for (i=0; i<NTAU; i++){
    fscanf(file, "%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le", &atmos.P[i], &atmos.T[i], &atmos.C[i], &atmos.CH4[i], &atmos.CO[i], &atmos.OCS[i], &atmos.CO2[i], &atmos.C2H2[i], &atmos.C2H4[i], &atmos.C2H6[i], &atmos.H[i], &atmos.HCN[i], &atmos.H2[i], &atmos.H2CO[i], &atmos.H2O[i], &atmos.H2S[i], &atmos.N[i], &atmos.N2[i], &atmos.NO2[i], &atmos.NH3[i], &atmos.NO[i], &atmos.O[i], &atmos.O2[i], &atmos.O3[i], &atmos.OH[i], &atmos.SO2[i]);
  }
  fclose(file);
  
  /* Determine mean molecular weight at each altitude */
  
  for (i=0;i<NTAU;i++){ 			
    
    atmos.mu[i] = 
      atmos.H[i] * m_H +
      atmos.H2O[i] * m_H2O +
      atmos.OH[i] * m_OH +
      atmos.O[i] * m_O +
      atmos.O2[i] * m_O2 +
      atmos.O3[i] * m_O3 +
      atmos.CO[i] * m_CO +
      atmos.CO2[i] * m_CO2 +
      atmos.C[i] * m_C +
      atmos.CH4[i] * m_CH4 +
      atmos.C2H2[i] * m_C2H2 +
      atmos.H2CO[i] * m_H2CO +
      atmos.N[i] * m_N +
      atmos.N2[i] * m_N2 +
      atmos.NO[i] * m_NO +
      atmos.NO2[i] * m_NO2 +
      atmos.NH3[i] * m_NH3 +
      atmos.HCN[i] * m_HCN +
      atmos.H2S[i] * m_H2S +
      atmos.SO2[i] * m_SO2 +
      atmos.OCS[i] * m_OCS +
      atmos.C2H4[i] * m_C2H4 +
      atmos.C2H6[i] * m_C2H6 +
      atmos.H2[i] * m_H2;
  }
  
  /* Cloud layer calculation */

  if(THRESHOLD != 0.0){			 
    double proportion = (log10(THRESHOLD) - log10(atmos.P[NTAU-2]))
      / (log10(atmos.P[NTAU-1]) - log10(atmos.P[NTAU-2]));  
    atmos.P[NTAU-1] = THRESHOLD;
    atmos.T[NTAU-1] = proportion*(atmos.T[NTAU-1] - atmos.T[NTAU-2]) 
      + atmos.T[NTAU-2];
    atmos.mu[NTAU-1] = proportion*(atmos.mu[NTAU-1] - atmos.mu[NTAU-2]) 
      + atmos.mu[NTAU-2];
  }
  
  
  printf("Last line of T-P profile:\n");
  printf("T:\t%e Pa\n", atmos.T[NTAU-1]);
  printf("P:\t%f K\n", atmos.P[NTAU-1]);
  printf("C:\t%e\n", atmos.C[NTAU-1]);
  printf("CH4:\t%e\n", atmos.CH4[NTAU-1]);
  printf("CO:\t%e\n", atmos.CO[NTAU-1]);
  printf("OCS:\t%e\n", atmos.OCS[NTAU-1]);
  printf("CO2:\t%e\n", atmos.CO2[NTAU-1]);
  printf("C2H2:\t%e\n", atmos.C2H2[NTAU-1]);
  printf("C2H4:\t%e\n", atmos.C2H4[NTAU-1]);
  printf("C2H6:\t%e\n", atmos.C2H6[NTAU-1]);
  printf("H:\t%e\n", atmos.H[NTAU-1]);
  printf("HCN:\t%e\n", atmos.HCN[NTAU-1]);
  printf("H2:\t%e\n", atmos.H2[NTAU-1]);
  printf("H2CO:\t%e\n", atmos.H2CO[NTAU-1]);
  printf("H2O:\t%e\n", atmos.H2O[NTAU-1]);
  printf("H2S:\t%e\n", atmos.H2S[NTAU-1]);
  printf("N:\t%e\n", atmos.N[NTAU-1]);
  printf("N2:\t%e\n", atmos.N2[NTAU-1]);
  printf("NO2:\t%e\n", atmos.NO2[NTAU-1]);
  printf("NH3:\t%e\n", atmos.NH3[NTAU-1]);
  printf("NO:\t%e\n", atmos.NO[NTAU-1]);
  printf("O:\t%e\n", atmos.O[NTAU-1]);
  printf("O2:\t%e\n", atmos.O2[NTAU-1]);
  printf("O3:\t%e\n", atmos.O3[NTAU-1]);
  printf("OH:\t%e\n", atmos.OH[NTAU-1]);
  printf("SO2:\t%e\n", atmos.SO2[NTAU-1]);
    
}

/* ------- end ----------------------- ReadTP.c ------------------ */

/* ------- start --------------------- FreeTP.c ------------------ */

void FreeTP(){
  
  /* Frees atmos structure */
  
  vars variables = getVars();
  int NTAU = variables.NTAU;
  int NLAMBDA = variables.NLAMBDA;

  free_dvector(atmos.P, 0, NTAU-1);
  free_dvector(atmos.T, 0, NTAU-1);
  free_dvector(atmos.mu, 0, NTAU-1);
  free_dvector(atmos.lambda, 0, NLAMBDA-1);

}

/* ------- end ----------------------- FreeTP.c ------------------ */
