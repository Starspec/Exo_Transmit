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

/*----------------------- totalopac.c ----------------------------

Author: Eliza Kempton (kemptone@grinnell.edu)

------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "include.h"
#include "constant.h"
#include "atmos.h"
#include "opac.h"
#include "nrutil.h"
#include "vars.h"
#include "prototypes.h"

/* --- Global variables ------------------------------------------ */

extern struct Opac opac;
extern struct Atmos atmos;

struct Chem chem;

struct Opac opacCH4;
struct Opac opacCO2;
struct Opac opacCO;
struct Opac opacH2O;
struct Opac opacNH3;
struct Opac opacO2;
struct Opac opacO3;
struct Opac opacC2H2;
struct Opac opacC2H4;
struct Opac opacC2H6;
struct Opac opacCrH; 
struct Opac opacH2CO;
struct Opac opacH2S; 
struct Opac opacHCl; 
struct Opac opacHCN; 
struct Opac opacHF;
struct Opac opacMgH; 
struct Opac opacN2; 
struct Opac opacNO; 
struct Opac opacNO2;
struct Opac opacOH; 
struct Opac opacOCS;
struct Opac opacPH3;
struct Opac opacSH; 
struct Opac opacSiH; 
struct Opac opacSiO;
struct Opac opacSO2; 
struct Opac opacTiO; 
struct Opac opacVO; 
struct Opac opacK; 
struct Opac opacNa; 

struct Opac opacscat;
struct Opac opacCIA;

/* ---------------------------------------------------------------
 * Computes the total opacity due to all of the atmospheric 
 * constituents.
 * --------------------------------------------------------------- */

/* ------- begin ------------ TotalOpac.c ------------------------ */

void TotalOpac() {

  double **opac_CIA_H2H2, **opac_CIA_H2He, **opac_CIA_H2H, 
    **opac_CIA_H2CH4, **opac_CIA_CH4Ar, **opac_CIA_CH4CH4, 
    **opac_CIA_CO2CO2, **opac_CIA_HeH, **opac_CIA_N2CH4, 
    **opac_CIA_N2H2, **opac_CIA_N2N2, **opac_CIA_O2CO2, 
    **opac_CIA_O2N2, **opac_CIA_O2O2;
  double **kappa_nu;
  int i, j, k, a, b;
  
  char **fileArray = getFileArray(); 	//get file names
  vars variables = getVars(); 		//get planet variables
  int chemSelection[32]; 		//get chemistry selections
  
  getChemSelection(chemSelection); 
  
  int NLAMBDA = variables.NLAMBDA; 	//rename relevant variables
  int NPRESSURE = variables.NPRESSURE;
  int NTEMP = variables.NTEMP;
  int NTAU = variables.NTAU;
  double RAYLEIGH = variables.RAYLEIGH;
  
  /* Molecular masses */
  /* Masses sourced from http://www.webqc.org/mmcalc.php */
  
  double m_H2 = 2.0158;
  double m_H = 1.0079;
  double m_He = 4.002602;
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
  double m_HCl= 36.4609;
  double m_HCN = 27.0253;
  double m_HF = 20.0063;
  double m_MgH = 25.3129;
  double m_NO = 30.0061;
  double m_NO2 = 46.0055;
  double m_OCS = 60.0751;
  double m_OH = 17.0073;
  double m_PH3 = 33.9976;
  double m_SH = 33.0729;
  double m_SiH = 29.0934;
  double m_SiO = 44.0849;
  double m_SO2 = 64.0638;
  double m_TiO = 63.8664;
  double m_VO = 66.9409;
  double m_Na = 22.988977; 
  double m_K = 39.0983; 
	
  FILE *f1;
  
  /* Allocate Memory */
  
  opac_CIA_H2H2 = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_H2He = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_H2H = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_H2CH4 = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_CH4Ar = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_CH4CH4 = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_CO2CO2 = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_HeH = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_N2CH4 = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_N2H2 = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_N2N2 = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_O2CO2 = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_O2N2 = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_O2O2 = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  
  /* Read Chemistry Table */
  
  ReadChemTable();
  printf("ReadChemTable done\n");
  
  /* Allocate for total opacity */

  opac.name = "Total";          //Name it Total
  opac.T = dvector(0, NTEMP-1); //Declare T, P, Plog10, and kappa arrays
  opac.P = dvector(0, NPRESSURE-1);
  opac.Plog10 = dvector(0, NPRESSURE-1);
  opac.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
  //opac.abundance = dvector(0, NTAU-1);
  
  //populate with zeros	
  for (i=0; i<NLAMBDA; i++)
    for (j=0; j<NPRESSURE; j++)
        for (k=0; k<NTEMP; k++)
	    opac.kappa[i][j][k] = 0.;
  
  /* Fill in mean molecular weight (mu) values */
  
  chem.mu = dvector(0, NTAU-1);
  
  for (j=0; j<NTAU; j++) {
      chem.mu[j] = 
	chem.H2[j]* m_H2 
	+ chem.H[j]* m_H 
	+ chem.He[j] * m_He 
	+ chem.H2O[j] * m_H2O 
	+ chem.C[j] * m_C 
	+ chem.CH4[j] * m_CH4 
	+ chem.CO[j] * m_CO 
	+ chem.CO2[j] * m_CO2 
	+ chem.O[j] * m_O 
	+ chem.N[j] * m_N 
	+ chem.NH3[j] * m_NH3 
	+ chem.N2[j] * m_N2 
	+ chem.O2[j] * m_O2 
	+ chem.O3[j] * m_O3
	+ chem.C2H2[j] * m_C2H2 
	+ chem.C2H4[j] * m_C2H4
	+ chem.C2H6[j] * m_C2H6
	+ chem.HCN[j] * m_HCN 
	+ chem.HCl[j] * m_HCl
	+ chem.HF[j] * m_HF 
	+ chem.H2CO[j] * m_H2CO
	+ chem.H2S[j] * m_H2S 
	+ chem.MgH[j] * m_MgH
	+ chem.NO2[j] * m_NO2 
	+ chem.NO[j] * m_NO
	+ chem.OH[j] * m_OH 
	+ chem.PH3[j] * m_PH3
	+ chem.SH[j] * m_SH 
	+ chem.SO2[j] * m_SO2
	+ chem.SiH[j] * m_SiH 
	+ chem.SiO[j] * m_SiO
	+ chem.TiO[j] * m_TiO
	+ chem.VO[j] * m_VO 
	+ chem.OCS[j] * m_OCS 
	+ chem.Na[j] * m_Na 
	+ chem.K[j] * m_K;
  }

  /* This adds in the calculation to atmos.mu */
  atmos.mu = dvector(0, NTAU-1);
  for (j=0; j<NTAU; j++) {
      atmos.mu[j] = chem.mu[j];
  };

  printf("atmos.mu has been populated from chem.mu\n");

  /* Allocate atmos.kappa_nu and kappa_nu, and populate with zeros */
  atmos.kappa_nu = dmatrix(0, NLAMBDA-1, 0, NTAU-1);
  kappa_nu = dmatrix(0, NLAMBDA-1, 0, NTAU-1);
  for (i = 0; i < NLAMBDA; i++) {
      for (j = 0; j < NTAU; j++) {
          atmos.kappa_nu[i][j] = 0.;
          kappa_nu[i][j] = 0.;
      };
  };

  printf("Made it here");
  
  /* Fill in CH4 opacities */
  if(chemSelection[0] == 1){          //If CH4 is selected
    opacCH4.name = "CH4";             //Name it CH4
    opacCH4.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacCH4.P = dvector(0, NPRESSURE-1);
    opacCH4.Plog10 = dvector(0, NPRESSURE-1);
    opacCH4.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacCH4, fileArray[3]);     //Read opacity table for CH4
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacCH4.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacCH4.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacCH4.T[a], opacCH4.P[b], chem.CH4[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacCH4.kappa[i][a][b] * chem.CH4[j];
      }
    };
    printf("Read CH4 Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacCH4);                  //Free CH4 opacity table
  };
  
  //This procedure repeats for all gases!!
  
  if(chemSelection[0] == 1){          //If CO2 is selected
    opacCO2.name = "CO2";             //Name it CO2
    opacCO2.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacCO2.P = dvector(0, NPRESSURE-1);
    opacCO2.Plog10 = dvector(0, NPRESSURE-1);
    opacCO2.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacCO2, fileArray[3]);     //Read opacity table for CO2
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacCO2.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacCO2.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacCO2.T[a], opacCO2.P[b], chem.CO2[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacCO2.kappa[i][a][b] * chem.CO2[j];
      }
    };
    printf("Read CO2 Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacCO2);                  //Free CO2 opacity table
  };
  
  if(chemSelection[0] == 1){          //If CO is selected
    opacCO.name = "CO";             //Name it CO
    opacCO.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacCO.P = dvector(0, NPRESSURE-1);
    opacCO.Plog10 = dvector(0, NPRESSURE-1);
    opacCO.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacCO, fileArray[3]);     //Read opacity table for CO
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacCO.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacCO.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacCO.T[a], opacCO.P[b], chem.CO[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacCO.kappa[i][a][b] * chem.CO[j];
      }
    };
    printf("Read CO Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacCO);                  //Free CO opacity table
  };
  
  if(chemSelection[0] == 1){          //If H2O is selected
    opacH2O.name = "H2O";             //Name it H2O
    opacH2O.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacH2O.P = dvector(0, NPRESSURE-1);
    opacH2O.Plog10 = dvector(0, NPRESSURE-1);
    opacH2O.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacH2O, fileArray[3]);     //Read opacity table for H2O
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacH2O.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacH2O.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacH2O.T[a], opacH2O.P[b], chem.H2O[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacH2O.kappa[i][a][b] * chem.H2O[j];
      }
    };
    printf("Read H2O Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacH2O);                  //Free H2O opacity table
  };
  
  if(chemSelection[0] == 1){          //If NH3 is selected
    opacNH3.name = "NH3";             //Name it NH3
    opacNH3.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacNH3.P = dvector(0, NPRESSURE-1);
    opacNH3.Plog10 = dvector(0, NPRESSURE-1);
    opacNH3.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacNH3, fileArray[3]);     //Read opacity table for NH3
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacNH3.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacNH3.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacNH3.T[a], opacNH3.P[b], chem.NH3[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacNH3.kappa[i][a][b] * chem.NH3[j];
      }
    };
    printf("Read NH3 Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacNH3);                  //Free NH3 opacity table
  };
  
  if(chemSelection[0] == 1){          //If O2 is selected
    opacO2.name = "O2";             //Name it O2
    opacO2.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacO2.P = dvector(0, NPRESSURE-1);
    opacO2.Plog10 = dvector(0, NPRESSURE-1);
    opacO2.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacO2, fileArray[3]);     //Read opacity table for O2
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacO2.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacO2.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacO2.T[a], opacO2.P[b], chem.O2[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacO2.kappa[i][a][b] * chem.O2[j];
      }
    };
    printf("Read O2 Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacO2);                  //Free O2 opacity table
  };
  
  if(chemSelection[0] == 1){          //If O3 is selected
    opacO3.name = "O3";             //Name it O3
    opacO3.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacO3.P = dvector(0, NPRESSURE-1);
    opacO3.Plog10 = dvector(0, NPRESSURE-1);
    opacO3.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacO3, fileArray[3]);     //Read opacity table for O3
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacO3.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacO3.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacO3.T[a], opacO3.P[b], chem.O3[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacO3.kappa[i][a][b] * chem.O3[j];
      }
    };
    printf("Read O3 Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacO3);                  //Free O3 opacity table
  };
  
  if(chemSelection[0] == 1){          //If C2H2 is selected
    opacC2H2.name = "C2H2";             //Name it C2H2
    opacC2H2.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacC2H2.P = dvector(0, NPRESSURE-1);
    opacC2H2.Plog10 = dvector(0, NPRESSURE-1);
    opacC2H2.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacC2H2, fileArray[3]);     //Read opacity table for C2H2
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacC2H2.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacC2H2.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacC2H2.T[a], opacC2H2.P[b], chem.C2H2[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacC2H2.kappa[i][a][b] * chem.C2H2[j];
      }
    };
    printf("Read C2H2 Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacC2H2);                  //Free C2H2 opacity table
  };
  
  if(chemSelection[0] == 1){          //If C2H4 is selected
    opacC2H4.name = "C2H4";             //Name it C2H4
    opacC2H4.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacC2H4.P = dvector(0, NPRESSURE-1);
    opacC2H4.Plog10 = dvector(0, NPRESSURE-1);
    opacC2H4.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacC2H4, fileArray[3]);     //Read opacity table for C2H4
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacC2H4.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacC2H4.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacC2H4.T[a], opacC2H4.P[b], chem.C2H4[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacC2H4.kappa[i][a][b] * chem.C2H4[j];
      }
    };
    printf("Read C2H4 Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacC2H4);                  //Free C2H4 opacity table
  };
  
  if(chemSelection[0] == 1){          //If C2H6 is selected
    opacC2H6.name = "C2H6";             //Name it C2H6
    opacC2H6.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacC2H6.P = dvector(0, NPRESSURE-1);
    opacC2H6.Plog10 = dvector(0, NPRESSURE-1);
    opacC2H6.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacC2H6, fileArray[3]);     //Read opacity table for C2H6
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacC2H6.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacC2H6.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacC2H6.T[a], opacC2H6.P[b], chem.C2H6[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacC2H6.kappa[i][a][b] * chem.C2H6[j];
      }
    };
    printf("Read C2H6 Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacC2H6);                  //Free C2H6 opacity table
  };
  
  if(chemSelection[0] == 1){          //If H2CO is selected
    opacH2CO.name = "H2CO";             //Name it H2CO
    opacH2CO.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacH2CO.P = dvector(0, NPRESSURE-1);
    opacH2CO.Plog10 = dvector(0, NPRESSURE-1);
    opacH2CO.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacH2CO, fileArray[3]);     //Read opacity table for H2CO
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacH2CO.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacH2CO.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacH2CO.T[a], opacH2CO.P[b], chem.H2CO[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacH2CO.kappa[i][a][b] * chem.H2CO[j];
      }
    };
    printf("Read H2CO Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacH2CO);                  //Free H2CO opacity table
  };
  
  if(chemSelection[0] == 1){          //If H2S is selected
    opacH2S.name = "H2S";             //Name it H2S
    opacH2S.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacH2S.P = dvector(0, NPRESSURE-1);
    opacH2S.Plog10 = dvector(0, NPRESSURE-1);
    opacH2S.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacH2S, fileArray[3]);     //Read opacity table for H2S
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacH2S.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacH2S.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacH2S.T[a], opacH2S.P[b], chem.H2S[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacH2S.kappa[i][a][b] * chem.H2S[j];
      }
    };
    printf("Read H2S Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacH2S);                  //Free H2S opacity table
  };
  
  if(chemSelection[0] == 1){          //If HCl is selected
    opacHCl.name = "HCl";             //Name it HCl
    opacHCl.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacHCl.P = dvector(0, NPRESSURE-1);
    opacHCl.Plog10 = dvector(0, NPRESSURE-1);
    opacHCl.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacHCl, fileArray[3]);     //Read opacity table for HCl
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacHCl.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacHCl.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacHCl.T[a], opacHCl.P[b], chem.HCl[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacHCl.kappa[i][a][b] * chem.HCl[j];
      }
    };
    printf("Read HCl Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacHCl);                  //Free HCl opacity table
  };
  
  if(chemSelection[0] == 1){          //If HCN is selected
    opacHCN.name = "HCN";             //Name it HCN
    opacHCN.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacHCN.P = dvector(0, NPRESSURE-1);
    opacHCN.Plog10 = dvector(0, NPRESSURE-1);
    opacHCN.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacHCN, fileArray[3]);     //Read opacity table for HCN
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacHCN.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacHCN.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacHCN.T[a], opacHCN.P[b], chem.HCN[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacHCN.kappa[i][a][b] * chem.HCN[j];
      }
    };
    printf("Read HCN Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacHCN);                  //Free HCN opacity table
  };
  
  if(chemSelection[0] == 1){          //If HF is selected
    opacHF.name = "HF";             //Name it HF
    opacHF.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacHF.P = dvector(0, NPRESSURE-1);
    opacHF.Plog10 = dvector(0, NPRESSURE-1);
    opacHF.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacHF, fileArray[3]);     //Read opacity table for HF
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacHF.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacHF.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacHF.T[a], opacHF.P[b], chem.HF[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacHF.kappa[i][a][b] * chem.HF[j];
      }
    };
    printf("Read HF Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacHF);                  //Free HF opacity table
  };
  
  if(chemSelection[0] == 1){          //If MgH is selected
    opacMgH.name = "MgH";             //Name it MgH
    opacMgH.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacMgH.P = dvector(0, NPRESSURE-1);
    opacMgH.Plog10 = dvector(0, NPRESSURE-1);
    opacMgH.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacMgH, fileArray[3]);     //Read opacity table for MgH
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacMgH.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacMgH.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacMgH.T[a], opacMgH.P[b], chem.MgH[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacMgH.kappa[i][a][b] * chem.MgH[j];
      }
    };
    printf("Read MgH Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacMgH);                  //Free MgH opacity table
  };
  
  if(chemSelection[0] == 1){          //If N2 is selected
    opacN2.name = "N2";             //Name it N2
    opacN2.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacN2.P = dvector(0, NPRESSURE-1);
    opacN2.Plog10 = dvector(0, NPRESSURE-1);
    opacN2.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacN2, fileArray[3]);     //Read opacity table for N2
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacN2.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacN2.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacN2.T[a], opacN2.P[b], chem.N2[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacN2.kappa[i][a][b] * chem.N2[j];
      }
    };
    printf("Read N2 Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacN2);                  //Free N2 opacity table
  };
  
  if(chemSelection[0] == 1){          //If NO is selected
    opacNO.name = "NO";             //Name it NO
    opacNO.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacNO.P = dvector(0, NPRESSURE-1);
    opacNO.Plog10 = dvector(0, NPRESSURE-1);
    opacNO.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacNO, fileArray[3]);     //Read opacity table for NO
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacNO.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacNO.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacNO.T[a], opacNO.P[b], chem.NO[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacNO.kappa[i][a][b] * chem.NO[j];
      }
    };
    printf("Read NO Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacNO);                  //Free NO opacity table
  };
  
  if(chemSelection[0] == 1){          //If NO2 is selected
    opacNO2.name = "NO2";             //Name it NO2
    opacNO2.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacNO2.P = dvector(0, NPRESSURE-1);
    opacNO2.Plog10 = dvector(0, NPRESSURE-1);
    opacNO2.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacNO2, fileArray[3]);     //Read opacity table for NO2
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacNO2.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacNO2.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacNO2.T[a], opacNO2.P[b], chem.NO2[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacNO2.kappa[i][a][b] * chem.NO2[j];
      }
    };
    printf("Read NO2 Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacNO2);                  //Free NO2 opacity table
  };
  
  if(chemSelection[0] == 1){          //If OCS is selected
    opacOCS.name = "OCS";             //Name it OCS
    opacOCS.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacOCS.P = dvector(0, NPRESSURE-1);
    opacOCS.Plog10 = dvector(0, NPRESSURE-1);
    opacOCS.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacOCS, fileArray[3]);     //Read opacity table for OCS
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacOCS.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacOCS.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacOCS.T[a], opacOCS.P[b], chem.OCS[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacOCS.kappa[i][a][b] * chem.OCS[j];
      }
    };
    printf("Read OCS Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacOCS);                  //Free OCS opacity table
  };
  
  if(chemSelection[0] == 1){          //If OH is selected
    opacOH.name = "OH";             //Name it OH
    opacOH.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacOH.P = dvector(0, NPRESSURE-1);
    opacOH.Plog10 = dvector(0, NPRESSURE-1);
    opacOH.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacOH, fileArray[3]);     //Read opacity table for OH
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacOH.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacOH.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacOH.T[a], opacOH.P[b], chem.OH[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacOH.kappa[i][a][b] * chem.OH[j];
      }
    };
    printf("Read OH Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacOH);                  //Free OH opacity table
  };
  
  if(chemSelection[0] == 1){          //If PH3 is selected
    opacPH3.name = "PH3";             //Name it PH3
    opacPH3.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacPH3.P = dvector(0, NPRESSURE-1);
    opacPH3.Plog10 = dvector(0, NPRESSURE-1);
    opacPH3.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacPH3, fileArray[3]);     //Read opacity table for PH3
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacPH3.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacPH3.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacPH3.T[a], opacPH3.P[b], chem.PH3[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacPH3.kappa[i][a][b] * chem.PH3[j];
      }
    };
    printf("Read PH3 Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacPH3);                  //Free PH3 opacity table
  };
  
  if(chemSelection[0] == 1){          //If SH is selected
    opacSH.name = "SH";             //Name it SH
    opacSH.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacSH.P = dvector(0, NPRESSURE-1);
    opacSH.Plog10 = dvector(0, NPRESSURE-1);
    opacSH.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacSH, fileArray[3]);     //Read opacity table for SH
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacSH.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacSH.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacSH.T[a], opacSH.P[b], chem.SH[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacSH.kappa[i][a][b] * chem.SH[j];
      }
    };
    printf("Read SH Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacSH);                  //Free SH opacity table
  };
  
  if(chemSelection[0] == 1){          //If SiH is selected
    opacSiH.name = "SiH";             //Name it SiH
    opacSiH.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacSiH.P = dvector(0, NPRESSURE-1);
    opacSiH.Plog10 = dvector(0, NPRESSURE-1);
    opacSiH.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacSiH, fileArray[3]);     //Read opacity table for SiH
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacSiH.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacSiH.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacSiH.T[a], opacSiH.P[b], chem.SiH[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacSiH.kappa[i][a][b] * chem.SiH[j];
      }
    };
    printf("Read SiH Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacSiH);                  //Free SiH opacity table
  };
  
  if(chemSelection[0] == 1){          //If SiO is selected
    opacSiO.name = "SiO";             //Name it SiO
    opacSiO.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacSiO.P = dvector(0, NPRESSURE-1);
    opacSiO.Plog10 = dvector(0, NPRESSURE-1);
    opacSiO.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacSiO, fileArray[3]);     //Read opacity table for SiO
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacSiO.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacSiO.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacSiO.T[a], opacSiO.P[b], chem.SiO[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacSiO.kappa[i][a][b] * chem.SiO[j];
      }
    };
    printf("Read SiO Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacSiO);                  //Free SiO opacity table
  };
  
  if(chemSelection[0] == 1){          //If SO2 is selected
    opacSO2.name = "SO2";             //Name it SO2
    opacSO2.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacSO2.P = dvector(0, NPRESSURE-1);
    opacSO2.Plog10 = dvector(0, NPRESSURE-1);
    opacSO2.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacSO2, fileArray[3]);     //Read opacity table for SO2
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacSO2.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacSO2.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacSO2.T[a], opacSO2.P[b], chem.SO2[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacSO2.kappa[i][a][b] * chem.SO2[j];
      }
    };
    printf("Read SO2 Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacSO2);                  //Free SO2 opacity table
  };
  
  if(chemSelection[0] == 1){          //If TiO is selected
    opacTiO.name = "TiO";             //Name it TiO
    opacTiO.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacTiO.P = dvector(0, NPRESSURE-1);
    opacTiO.Plog10 = dvector(0, NPRESSURE-1);
    opacTiO.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacTiO, fileArray[3]);     //Read opacity table for TiO
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacTiO.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacTiO.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacTiO.T[a], opacTiO.P[b], chem.TiO[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacTiO.kappa[i][a][b] * chem.TiO[j];
      }
    };
    printf("Read TiO Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacTiO);                  //Free TiO opacity table
  };
  
  if(chemSelection[0] == 1){          //If VO is selected
    opacVO.name = "VO";             //Name it VO
    opacVO.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacVO.P = dvector(0, NPRESSURE-1);
    opacVO.Plog10 = dvector(0, NPRESSURE-1);
    opacVO.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacVO, fileArray[3]);     //Read opacity table for VO
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacVO.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacVO.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacVO.T[a], opacVO.P[b], chem.VO[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacVO.kappa[i][a][b] * chem.VO[j];
      }
    };
    printf("Read VO Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacVO);                  //Free VO opacity table
  };

  /* Atomic opacities */

  if(chemSelection[0] == 1){          //If Na is selected
    opacNa.name = "Na";             //Name it Na
    opacNa.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacNa.P = dvector(0, NPRESSURE-1);
    opacNa.Plog10 = dvector(0, NPRESSURE-1);
    opacNa.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacNa, fileArray[3]);     //Read opacity table for Na
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacNa.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacNa.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacNa.T[a], opacNa.P[b], chem.Na[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacNa.kappa[i][a][b] * chem.Na[j];
      }
    };
    printf("Read Na Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacNa);                  //Free Na opacity table
  };

  if(chemSelection[0] == 1){          //If K is selected
    opacK.name = "K";             //Name it K
    opacK.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacK.P = dvector(0, NPRESSURE-1);
    opacK.Plog10 = dvector(0, NPRESSURE-1);
    opacK.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    
    ReadOpacTable(opacK, fileArray[3]);     //Read opacity table for K
    
    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacK.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacK.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          //if (j == 0) {
          //  printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
          //        i, j, a, b, opacK.T[a], opacK.P[b], chem.K[j]);
          //};

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacK.kappa[i][a][b] * chem.K[j];
      }
    };
    printf("Read K Opacity done\n");	     //Confirmation message
    
    FreeOpacTable(opacK);                  //Free K opacity table
  };
  
  /* Fill in total opacities */
  
  for(k=0; k<NTEMP; k++)
    opac.T[k] = chem.T[k];	      //insert temperatures
  
  for(j=0; j<NPRESSURE; j++){	      //insert pressues
    opac.P[j] = chem.P[j];
    opac.Plog10[j] = log10(chem.P[j]);
  }
  

  /* Fill in collision-induced opacities */
  if(chemSelection[31] == 1){	      //If CIA is turned on ...
    
    /* Allocate collison induced opacities */
    opacCIA.name = "CIA";
    opacCIA.T = dvector(0, NTEMP-1);
    opacCIA.P = dvector(0, NPRESSURE-1);
    opacCIA.Plog10 = dvector(0, NPRESSURE-1);
    opacCIA.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    //opacCIA.abundance = dvector(0, NTAU-1);
    
    /* populate with zeros */	
    for (i=0; i<NLAMBDA; i++) {
      for (j=0; j<NPRESSURE; j++) {
            for (k=0; k<NTEMP; k++) {
	        opacCIA.kappa[i][j][k] = 0.;
            };
      };
    };
    
    /* Read in CIA opacities */
    
    f1 = fopen(fileArray[33], "r");
    if(f1 == NULL){
      printf("\n totalopac.c:\nError opening file: %s -- No such file or directory\n\n", fileArray[33]);
      exit(1);
    }
    
    for(i=0; i<NTEMP; i++){
      fscanf(f1, "%le", &opacCIA.T[i]);
      printf("%1e\n", opacCIA.T[i]);
    }

    for (k=0; k<NTEMP; k++){
      fscanf(f1, "%le", &opacCIA.T[k]);
      for (i=0; i<NLAMBDA; i++){
        fscanf(f1, "%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le",
                &atmos.lambda[i], &opac_CIA_H2H2[k][i],
                &opac_CIA_H2He[k][i], &opac_CIA_H2H[k][i],
                &opac_CIA_H2CH4[k][i], &opac_CIA_CH4Ar[k][i],
                &opac_CIA_CH4CH4[k][i], &opac_CIA_CO2CO2[k][i],
                &opac_CIA_HeH[k][i], &opac_CIA_N2CH4[k][i],
                &opac_CIA_N2H2[k][i], &opac_CIA_N2N2[k][i],
                &opac_CIA_O2CO2[k][i], &opac_CIA_O2N2[k][i],
                &opac_CIA_O2O2[k][i]);
      }
    };
    

    /* Populate atmos.kappa_nu */
    for (i=0; i<NLAMBDA; i++){
      for (j=0; j<NTAU; j++) {
          /* Interpolate from TP grid onto the altitude grid */
          Locate(NTEMP, opacCIA.T, atmos.T[j], &a);
          Locate(NPRESSURE, opacCIA.P, atmos.P[j], &b);


          kappa_nu[i][j] = lint2D(opac.T[a], opac.T[a+1], opac.P[b], opac.P[b+1],
      			      opac.kappa[i][b][a], opac.kappa[i][b][a+1],
      			      opac.kappa[i][b+1][a], opac.kappa[i][b+1][a+1],
      			      atmos.T[j], atmos.P[j]);

          if (j == 0) {
            printf("i = %d, j = %d, a = %d, b = %d, T = %e, P = %e, MR = %e\n",
                  i, j, a, b, opacCIA.T[a], opacCIA.P[b], chem.K[j]);
          };

	  /* Add to overall opac.kappa */
	  atmos.kappa_nu[i][j] += opacCIA.kappa[i][a][b] * chem.K[j];
      }
    };
    
    FreeOpacTable(opacCIA);
  }
 
  /* Rayleigh scattering */
  
  /* (Polarizabilities from the CRC Handbook) */
  
  if(chemSelection[30] == 1){		//If Scattering is activated.. 
    
    /* Declare memory for opacscat structure */
    
    opacscat.name = "Scat";
    opacscat.T = dvector(0, NTEMP-1);
    opacscat.P = dvector(0, NPRESSURE-1);
    opacscat.Plog10 = dvector(0, NPRESSURE-1);
    opacscat.kappa = d3tensor(0, NLAMBDA-1, 0, NTAU-1, 0, 1);
    //opacscat.abundance = dvector(0, NTAU-1);
    
    //populate with zeros	
    for (i=0; i<NLAMBDA; i++)
      for (j=0; j<NPRESSURE; j++)
          for (k=0; k<NTEMP; k++)
	  opacscat.kappa[i][j][k] = 0.;
    
    /* Fill in scattering coefficients */
    for (i=0; i<NLAMBDA; i++) {
      for (j=0; j<NTAU; j++) {
          // This is a temporary bypass of the third loop since this
          // calculation does not directly depend on anything beyond known
          // quantities on the atmosphere grid -- Teal
          k = 0;

	  /* Add Rayleigh scattering polarizability to overall kappa */
	  opacscat.kappa[i][j][k] +=
	    (8.0*PI/3.0) * SQ(0.80e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    chem.H2[j]*chem.P[j] / (KBOLTZMANN * chem.T[j])
	    +
	    (8.0*PI/3.0) * SQ(0.21e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    chem.He[j]*chem.P[j] / (KBOLTZMANN * chem.T[j])
	    +
	    (8.0*PI/3.0) * SQ(1.74e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    chem.N2[j]*chem.P[j] / (KBOLTZMANN * chem.T[j])
	    +
	    (8.0*PI/3.0) * SQ(1.45e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    chem.H2O[j]*chem.P[j] / (KBOLTZMANN * chem.T[j])
	    +
	    (8.0*PI/3.0) * SQ(1.95e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    chem.CO[j]*chem.P[j] / (KBOLTZMANN * chem.T[j])
	    +
	    (8.0*PI/3.0) * SQ(2.91e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    chem.CO2[j]*chem.P[j] / (KBOLTZMANN * chem.T[j])
	    +
	    (8.0*PI/3.0) * SQ(2.26e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    chem.NH3[j]*chem.P[j] / (KBOLTZMANN * chem.T[j])
	    +
	    (8.0*PI/3.0) * SQ(2.59e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    chem.CH4[j]*chem.P[j] / (KBOLTZMANN * chem.T[j])
	    +
	    (8.0*PI/3.0) * SQ(1.58e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    chem.O2[j]*chem.P[j] / (KBOLTZMANN * chem.T[j])
	    +
	    (8.0*PI/3.0) * SQ(3.21e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    chem.O3[j]*chem.P[j] / (KBOLTZMANN * chem.T[j])
	    +
	    (8.0*PI/3.0) * SQ(3.33e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    chem.C2H2[j]*chem.P[j] / (KBOLTZMANN * chem.T[j])
	    +
	    (8.0*PI/3.0) * SQ(4.25e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    chem.C2H4[j]*chem.P[j] / (KBOLTZMANN * chem.T[j])
	    +
	    (8.0*PI/3.0) * SQ(4.47e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    chem.C2H6[j]*chem.P[j] / (KBOLTZMANN * chem.T[j])
	    +
	    (8.0*PI/3.0) * SQ(2.59e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    chem.HCN[j]*chem.P[j] / (KBOLTZMANN * chem.T[j])
	    +
	    (8.0*PI/3.0) * SQ(2.63e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    chem.HCl[j]*chem.P[j] / (KBOLTZMANN * chem.T[j])
	    +
	    (8.0*PI/3.0) * SQ(0.80e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    chem.HF[j]*chem.P[j] / (KBOLTZMANN * chem.T[j])
	    +
	    (8.0*PI/3.0) * SQ(3.78e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    chem.H2S[j]*chem.P[j] / (KBOLTZMANN * chem.T[j])
	    +
	    (8.0*PI/3.0) * SQ(1.70e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    chem.NO[j]*chem.P[j] / (KBOLTZMANN * chem.T[j])
	    +
	    (8.0*PI/3.0) * SQ(3.02e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    chem.NO2[j]*chem.P[j] / (KBOLTZMANN * chem.T[j])
	    +
	    (8.0*PI/3.0) * SQ(4.84e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    chem.PH3[j]*chem.P[j] / (KBOLTZMANN * chem.T[j]);
	  
	  opacscat.kappa[i][j][k] *= 0.; // RAYLEIGH;
          atmos.kappa_nu[i][j] += opacscat.kappa[i][j][k];
	}
      }
    
    FreeOpacTable(opacscat);
  }

  /*  Free memory */

  free_dmatrix(opac_CIA_H2H2, 0, NPRESSURE-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_H2He, 0, NPRESSURE-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_H2H, 0, NPRESSURE-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_H2CH4, 0, NPRESSURE-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_CH4Ar, 0, NPRESSURE-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_CH4CH4, 0, NPRESSURE-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_CO2CO2, 0, NPRESSURE-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_HeH, 0, NPRESSURE-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_N2CH4, 0, NPRESSURE-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_N2H2, 0, NPRESSURE-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_N2N2, 0, NPRESSURE-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_O2CO2, 0, NPRESSURE-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_O2N2, 0, NPRESSURE-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_O2O2, 0, NPRESSURE-1, 0, NLAMBDA-1);
  
}

/* ------- end -------------- TotalOpac.c ------------------------ */

