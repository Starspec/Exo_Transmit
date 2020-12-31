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

struct Opac opacCIA;

/* ---------------------------------------------------------------
 * Computes the total opacity due to all of the atmospheric 
 * constituents.
 * --------------------------------------------------------------- */

/* ------- begin ------------ TotalOpac.c ------------------------ */

void TotalOpac() {

  double CIA_temp;
  double **opac_CIA_H2H2, **opac_CIA_H2He, **opac_CIA_H2H, 
    **opac_CIA_H2CH4, **opac_CIA_CH4Ar, **opac_CIA_CH4CH4, 
    **opac_CIA_CO2CO2, **opac_CIA_HeH, **opac_CIA_N2CH4, 
    **opac_CIA_N2H2, **opac_CIA_N2N2, **opac_CIA_O2CO2, 
    **opac_CIA_O2N2, **opac_CIA_O2O2;
  int i, j, k, ll, a, b;
  
  char **fileArray = getFileArray(); 	//get file names
  vars variables = getVars(); 		//get planet variables
  int chemSelection[32]; 		//get chemistry selections
  
  getChemSelection(chemSelection); 
  
  int NLAMBDA = variables.NLAMBDA; 	//rename relevant variables
  int NPRESSURE = variables.NPRESSURE;
  int NTEMP = variables.NTEMP;
  int NTAU = variables.NTAU;
  double RAYLEIGH = variables.RAYLEIGH;
  
  FILE *f1;
  
  /* Allocate Memory */

  atmos.kappa = dmatrix(0, NLAMBDA-1, 0, NTAU-1);
  
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
  
  //populate with zeros	
  for (i=0; i<NLAMBDA; i++)
    for (ll=0; ll<NTAU; ll++)
      atmos.kappa[i][ll] = 0.;
  
  /* Fill in CH4 opacities */
  if(chemSelection[0] == 1){          //If CH4 is selected
    opacCH4.name = "CH4";             //Name it CH4
    opacCH4.T = dvector(0, NTEMP-1);  //Declare T, P, Plog10, and kappa arrays
    opacCH4.P = dvector(0, NPRESSURE-1);
    opacCH4.Plog10 = dvector(0, NPRESSURE-1);
    opacCH4.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    /* Declare abundance */
    opacCH4.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
      
    ReadOpacTable(opacCH4, fileArray[3]);     //Read opacity table for CH4
    
    printf("Read CH4 Opacity done\n");	     //Confirmation message

    /* Interpolate and add to overall atmos.kappa */
    for(ll=0; ll<NTAU; ll++){
      Locate(NTEMP, opacCH4.T, atmos.T[ll], &a);
      Locate(NPRESSURE, opacCH4.P, atmos.P[ll], &b);
      for(i=0; i<NLAMBDA; i++){
	atmos.kappa[i][ll] += lint2D(
		opacCH4.T[a], opacCH4.T[a+1], opacCH4.P[b], opacCH4.P[b+1],
		opacCH4.kappa[i][b][a], opacCH4.kappa[i][b][a+1],
		opacCH4.kappa[i][b+1][a], opacCH4.kappa[i][b+1][a+1],
		atmos.T[ll], atmos.P[ll])
				     * atmos.CH4[ll];
      }
    }
    
    FreeOpacTable(opacCH4);                  //Free CH4 opacity table
  }
  
  //This procedure repeats for all gases!!
  
  /* Fill in CO2 opacities */
  if(chemSelection[1] == 1){
    opacCO2.name = "CO2";
    opacCO2.T = dvector(0, NTEMP-1);
    opacCO2.P = dvector(0, NPRESSURE-1);
    opacCO2.Plog10 = dvector(0, NPRESSURE-1);
    opacCO2.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacCO2.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacCO2, fileArray[8]);
    
    printf("Read CO2 Opacity done\n");

    /* Interpolate and add to overall atmos.kappa */
    for(ll=0; ll<NTAU; ll++){
      Locate(NTEMP, opacCO2.T, atmos.T[ll], &a);
      Locate(NPRESSURE, opacCO2.P, atmos.P[ll], &b);
      for(i=0; i<NLAMBDA; i++){
    	atmos.kappa[i][ll] += lint2D(
    		opacCO2.T[a], opacCO2.T[a+1], opacCO2.P[b], opacCO2.P[b+1],
    		opacCO2.kappa[i][b][a], opacCO2.kappa[i][b][a+1],
    		opacCO2.kappa[i][b+1][a], opacCO2.kappa[i][b+1][a+1],
    		atmos.T[ll], atmos.P[ll])
    				     * atmos.CO2[ll];
      }
    }

    FreeOpacTable(opacCO2);
  }
  
  /* Fill in CO opacities */
  if(chemSelection[2] == 1){
    opacCO.name = "CO";
    opacCO.T = dvector(0, NTEMP-1);
    opacCO.P = dvector(0, NPRESSURE-1);
    opacCO.Plog10 = dvector(0, NPRESSURE-1);
    opacCO.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacCO.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacCO, fileArray[7]);
    
    printf("Read CO Opacity done\n");

    /* Interpolate and add to overall atmos.kappa */
    for(ll=0; ll<NTAU; ll++){
      Locate(NTEMP, opacCO.T, atmos.T[ll], &a);
      Locate(NPRESSURE, opacCO.P, atmos.P[ll], &b);
      for(i=0; i<NLAMBDA; i++){
    	atmos.kappa[i][ll] += lint2D(
    		opacCO.T[a], opacCO.T[a+1], opacCO.P[b], opacCO.P[b+1],
    		opacCO.kappa[i][b][a], opacCO.kappa[i][b][a+1],
    		opacCO.kappa[i][b+1][a], opacCO.kappa[i][b+1][a+1],
    		atmos.T[ll], atmos.P[ll])
    				     * atmos.CO[ll];
      }
    }

    FreeOpacTable(opacCO);
  }
  
  /* Fill in H2O opacities */
  if(chemSelection[3] == 1){
    opacH2O.name = "H2O";
    opacH2O.T = dvector(0, NTEMP-1);
    opacH2O.P = dvector(0, NPRESSURE-1);
    opacH2O.Plog10 = dvector(0, NPRESSURE-1);
    opacH2O.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacH2O.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacH2O, fileArray[10]);
    
    printf("Read H2O Opacity done\n");

    /* Interpolate and add to overall atmos.kappa */
    for(ll=0; ll<NTAU; ll++){
      Locate(NTEMP, opacH2O.T, atmos.T[ll], &a);
      Locate(NPRESSURE, opacH2O.P, atmos.P[ll], &b);
      for(i=0; i<NLAMBDA; i++){
    	atmos.kappa[i][ll] += lint2D(
    		opacH2O.T[a], opacH2O.T[a+1], opacH2O.P[b], opacH2O.P[b+1],
    		opacH2O.kappa[i][b][a], opacH2O.kappa[i][b][a+1],
    		opacH2O.kappa[i][b+1][a], opacH2O.kappa[i][b+1][a+1],
    		atmos.T[ll], atmos.P[ll])
    				     * atmos.H2O[ll];
      }
    }

    FreeOpacTable(opacH2O);
  }
  
  /* Fill in NH3 opacities */
  if(chemSelection[4] == 1){
    opacNH3.name = "NH3";
    opacNH3.T = dvector(0, NTEMP-1);
    opacNH3.P = dvector(0, NPRESSURE-1);
    opacNH3.Plog10 = dvector(0, NPRESSURE-1);
    opacNH3.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacNH3.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);  
    ReadOpacTable(opacNH3, fileArray[17]);
    
    printf("Read NH3 Opacity done\n");

    /* Interpolate and add to overall atmos.kappa */
    for(ll=0; ll<NTAU; ll++){
      Locate(NTEMP, opacNH3.T, atmos.T[ll], &a);
      Locate(NPRESSURE, opacNH3.P, atmos.P[ll], &b);
      for(i=0; i<NLAMBDA; i++){
    	atmos.kappa[i][ll] += lint2D(
    		opacNH3.T[a], opacNH3.T[a+1], opacNH3.P[b], opacNH3.P[b+1],
    		opacNH3.kappa[i][b][a], opacNH3.kappa[i][b][a+1],
    		opacNH3.kappa[i][b+1][a], opacNH3.kappa[i][b+1][a+1],
    		atmos.T[ll], atmos.P[ll])
    				     * atmos.NH3[ll];
      }
    }

    FreeOpacTable(opacNH3);
  }
  
  /* Fill in O2 opacities */
  if(chemSelection[5] == 1){
    opacO2.name = "O2";
    opacO2.T = dvector(0, NTEMP-1);
    opacO2.P = dvector(0, NPRESSURE-1);
    opacO2.Plog10 = dvector(0, NPRESSURE-1);
    opacO2.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacO2.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacO2, fileArray[20]);
    
    printf("Read O2 Opacity done\n");

    /* Interpolate and add to overall atmos.kappa */
    for(ll=0; ll<NTAU; ll++){
      Locate(NTEMP, opacO2.T, atmos.T[ll], &a);
      Locate(NPRESSURE, opacO2.P, atmos.P[ll], &b);
      for(i=0; i<NLAMBDA; i++){
    	atmos.kappa[i][ll] += lint2D(
    		opacO2.T[a], opacO2.T[a+1], opacO2.P[b], opacO2.P[b+1],
    		opacO2.kappa[i][b][a], opacO2.kappa[i][b][a+1],
    		opacO2.kappa[i][b+1][a], opacO2.kappa[i][b+1][a+1],
    		atmos.T[ll], atmos.P[ll])
    				     * atmos.O2[ll];
      }
    }

    FreeOpacTable(opacO2);
  }
  
  /* Fill in O3 opacities */
  if(chemSelection[6] == 1){
    opacO3.name = "O3";
    opacO3.T = dvector(0, NTEMP-1);
    opacO3.P = dvector(0, NPRESSURE-1);
    opacO3.Plog10 = dvector(0, NPRESSURE-1);
    opacO3.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacO3.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacO3, fileArray[21]);
    
    printf("Read O3 Opacity done\n");

    /* Interpolate and add to overall atmos.kappa */
    for(ll=0; ll<NTAU; ll++){
      Locate(NTEMP, opacO3.T, atmos.T[ll], &a);
      Locate(NPRESSURE, opacO3.P, atmos.P[ll], &b);
      for(i=0; i<NLAMBDA; i++){
    	atmos.kappa[i][ll] += lint2D(
    		opacO3.T[a], opacO3.T[a+1], opacO3.P[b], opacO3.P[b+1],
    		opacO3.kappa[i][b][a], opacO3.kappa[i][b][a+1],
    		opacO3.kappa[i][b+1][a], opacO3.kappa[i][b+1][a+1],
    		atmos.T[ll], atmos.P[ll])
    				     * atmos.O3[ll];
      }
    }
    
    FreeOpacTable(opacO3);
  }
  
  /* Fill in C2H2 opacities */
  if(chemSelection[7] == 1){
    opacC2H2.name = "C2H2";
    opacC2H2.T = dvector(0, NTEMP-1);
    opacC2H2.P = dvector(0, NPRESSURE-1);
    opacC2H2.Plog10 = dvector(0, NPRESSURE-1);
    opacC2H2.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacC2H2.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacC2H2, fileArray[4]);
    
    printf("Read C2H2 Opacity done\n");

    /* Interpolate and add to overall atmos.kappa */
    for(ll=0; ll<NTAU; ll++){
      Locate(NTEMP, opacC2H2.T, atmos.T[ll], &a);
      Locate(NPRESSURE, opacC2H2.P, atmos.P[ll], &b);
      for(i=0; i<NLAMBDA; i++){
    	atmos.kappa[i][ll] += lint2D(
    		opacC2H2.T[a], opacC2H2.T[a+1], opacC2H2.P[b], opacC2H2.P[b+1],
    		opacC2H2.kappa[i][b][a], opacC2H2.kappa[i][b][a+1],
    		opacC2H2.kappa[i][b+1][a], opacC2H2.kappa[i][b+1][a+1],
    		atmos.T[ll], atmos.P[ll])
    				     * atmos.C2H2[ll];
      }
    }

    FreeOpacTable(opacC2H2);
  }
  
  /* Fill in C2H4 opacities */
  if(chemSelection[8] == 1){
    opacC2H4.name = "C2H4";
    opacC2H4.T = dvector(0, NTEMP-1);
    opacC2H4.P = dvector(0, NPRESSURE-1);
    opacC2H4.Plog10 = dvector(0, NPRESSURE-1);
    opacC2H4.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacC2H4.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacC2H4, fileArray[5]);
    
    printf("Read C2H4 Opacity done\n");

    /* Interpolate and add to overall atmos.kappa */
    for(ll=0; ll<NTAU; ll++){
      Locate(NTEMP, opacC2H4.T, atmos.T[ll], &a);
      Locate(NPRESSURE, opacC2H4.P, atmos.P[ll], &b);
      for(i=0; i<NLAMBDA; i++){
    	atmos.kappa[i][ll] += lint2D(
    		opacC2H4.T[a], opacC2H4.T[a+1], opacC2H4.P[b], opacC2H4.P[b+1],
    		opacC2H4.kappa[i][b][a], opacC2H4.kappa[i][b][a+1],
    		opacC2H4.kappa[i][b+1][a], opacC2H4.kappa[i][b+1][a+1],
    		atmos.T[ll], atmos.P[ll])
    				     * atmos.C2H4[ll];
      }
    }

    FreeOpacTable(opacC2H4);
  }
  
  /* Fill in C2H6 opacities */
  if(chemSelection[9] == 1){
    opacC2H6.name = "C2H6";
    opacC2H6.T = dvector(0, NTEMP-1);
    opacC2H6.P = dvector(0, NPRESSURE-1);
    opacC2H6.Plog10 = dvector(0, NPRESSURE-1);
    opacC2H6.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacC2H6.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacC2H6, fileArray[6]);
    
    printf("Read C2H6 Opacity done\n");

    /* Interpolate and add to overall atmos.kappa */
    for(ll=0; ll<NTAU; ll++){
      Locate(NTEMP, opacC2H6.T, atmos.T[ll], &a);
      Locate(NPRESSURE, opacC2H6.P, atmos.P[ll], &b);
      for(i=0; i<NLAMBDA; i++){
    	atmos.kappa[i][ll] += lint2D(
    		opacC2H6.T[a], opacC2H6.T[a+1], opacC2H6.P[b], opacC2H6.P[b+1],
    		opacC2H6.kappa[i][b][a], opacC2H6.kappa[i][b][a+1],
    		opacC2H6.kappa[i][b+1][a], opacC2H6.kappa[i][b+1][a+1],
    		atmos.T[ll], atmos.P[ll])
    				     * atmos.C2H6[ll];
      }
    }

    FreeOpacTable(opacC2H6);
  }
  
  /* Fill in H2CO opacities */
  if(chemSelection[10] == 1){
    opacH2CO.name = "H2CO";
    opacH2CO.T = dvector(0, NTEMP-1);
    opacH2CO.P = dvector(0, NPRESSURE-1);
    opacH2CO.Plog10 = dvector(0, NPRESSURE-1);
    opacH2CO.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacH2CO.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacH2CO, fileArray[9]);
    
    printf("Read H2CO Opacity done\n");
   
    /* Interpolate and add to overall atmos.kappa */
    for(ll=0; ll<NTAU; ll++){
      Locate(NTEMP, opacH2CO.T, atmos.T[ll], &a);
      Locate(NPRESSURE, opacH2CO.P, atmos.P[ll], &b);
      for(i=0; i<NLAMBDA; i++){
    	atmos.kappa[i][ll] += lint2D(
    		opacH2CO.T[a], opacH2CO.T[a+1], opacH2CO.P[b], opacH2CO.P[b+1],
    		opacH2CO.kappa[i][b][a], opacH2CO.kappa[i][b][a+1],
    		opacH2CO.kappa[i][b+1][a], opacH2CO.kappa[i][b+1][a+1],
    		atmos.T[ll], atmos.P[ll])
    				     * atmos.H2CO[ll];
      }
    }
 
    FreeOpacTable(opacH2CO);
  }
  
  /* Fill in H2S opacities */
  if(chemSelection[11] == 1){
    opacH2S.name = "H2S";
    opacH2S.T = dvector(0, NTEMP-1);
    opacH2S.P = dvector(0, NPRESSURE-1);
    opacH2S.Plog10 = dvector(0, NPRESSURE-1);
    opacH2S.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacH2S.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacH2S, fileArray[11]);
    
    printf("Read H2S Opacity done\n");

    /* Interpolate and add to overall atmos.kappa */
    for(ll=0; ll<NTAU; ll++){
      Locate(NTEMP, opacH2S.T, atmos.T[ll], &a);
      Locate(NPRESSURE, opacH2S.P, atmos.P[ll], &b);
      for(i=0; i<NLAMBDA; i++){
    	atmos.kappa[i][ll] += lint2D(
    		opacH2S.T[a], opacH2S.T[a+1], opacH2S.P[b], opacH2S.P[b+1],
    		opacH2S.kappa[i][b][a], opacH2S.kappa[i][b][a+1],
    		opacH2S.kappa[i][b+1][a], opacH2S.kappa[i][b+1][a+1],
    		atmos.T[ll], atmos.P[ll])
    				     * atmos.H2S[ll];
      }
    }

    FreeOpacTable(opacH2S);
  }
  
  /* Fill in HCl opacities */
  if(chemSelection[12] == 1){
    opacHCl.name = "HCl";
    opacHCl.T = dvector(0, NTEMP-1);
    opacHCl.P = dvector(0, NPRESSURE-1);
    opacHCl.Plog10 = dvector(0, NPRESSURE-1);
    opacHCl.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacHCl.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacHCl, fileArray[13]);
    
    printf("Read HCl Opacity done\n");

    FreeOpacTable(opacHCl);
  }
  
  /* Fill in HCN opacities */
  if(chemSelection[13] == 1){
    opacHCN.name = "HCN";
    opacHCN.T = dvector(0, NTEMP-1);
    opacHCN.P = dvector(0, NPRESSURE-1);
    opacHCN.Plog10 = dvector(0, NPRESSURE-1);
    opacHCN.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacHCN.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacHCN, fileArray[12]);
    
    printf("Read HCN Opacity done\n");

    /* Interpolate and add to overall atmos.kappa */
    for(ll=0; ll<NTAU; ll++){
      Locate(NTEMP, opacHCN.T, atmos.T[ll], &a);
      Locate(NPRESSURE, opacHCN.P, atmos.P[ll], &b);
      for(i=0; i<NLAMBDA; i++){
    	atmos.kappa[i][ll] += lint2D(
    		opacHCN.T[a], opacHCN.T[a+1], opacHCN.P[b], opacHCN.P[b+1],
    		opacHCN.kappa[i][b][a], opacHCN.kappa[i][b][a+1],
    		opacHCN.kappa[i][b+1][a], opacHCN.kappa[i][b+1][a+1],
    		atmos.T[ll], atmos.P[ll])
    				     * atmos.HCN[ll];
      }
    }
    
    FreeOpacTable(opacHCN);
  }
  
  /* Fill in HF opacities */
  if(chemSelection[14] == 1){
    opacHF.name = "HF";
    opacHF.T = dvector(0, NTEMP-1);
    opacHF.P = dvector(0, NPRESSURE-1);
    opacHF.Plog10 = dvector(0, NPRESSURE-1);
    opacHF.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacHF.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacHF, fileArray[14]);
    
    printf("Read HF Opacity done\n");

    FreeOpacTable(opacHF);
  }
  
  /* Fill in MgH opacities */
  if(chemSelection[15] == 1){
    opacMgH.name = "MgH";
    opacMgH.T = dvector(0, NTEMP-1);
    opacMgH.P = dvector(0, NPRESSURE-1);
    opacMgH.Plog10 = dvector(0, NPRESSURE-1);
    opacMgH.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacMgH.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacMgH, fileArray[15]);
    
    printf("Read MgH Opacity done\n");

    FreeOpacTable(opacMgH);
  }
  
  /* Fill in N2 opacities */
  if(chemSelection[16] == 1){
    opacN2.name = "N2";
    opacN2.T = dvector(0, NTEMP-1);
    opacN2.P = dvector(0, NPRESSURE-1);
    opacN2.Plog10 = dvector(0, NPRESSURE-1);
    opacN2.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacN2.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacN2, fileArray[16]);
    
    printf("Read N2 Opacity done\n");

    /* Interpolate and add to overall atmos.kappa */
    for(ll=0; ll<NTAU; ll++){
      Locate(NTEMP, opacN2.T, atmos.T[ll], &a);
      Locate(NPRESSURE, opacN2.P, atmos.P[ll], &b);
      for(i=0; i<NLAMBDA; i++){
    	atmos.kappa[i][ll] += lint2D(
    		opacN2.T[a], opacN2.T[a+1], opacN2.P[b], opacN2.P[b+1],
    		opacN2.kappa[i][b][a], opacN2.kappa[i][b][a+1],
    		opacN2.kappa[i][b+1][a], opacN2.kappa[i][b+1][a+1],
    		atmos.T[ll], atmos.P[ll])
    				     * atmos.N2[ll];
      }
    }

    FreeOpacTable(opacN2);
  }
  
  /* Fill in NO opacities */
  if(chemSelection[17] == 1){
    opacNO.name = "NO";
    opacNO.T = dvector(0, NTEMP-1);
    opacNO.P = dvector(0, NPRESSURE-1);
    opacNO.Plog10 = dvector(0, NPRESSURE-1);
    opacNO.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacNO.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacNO, fileArray[18]);
    
    printf("Read NO Opacity done\n");

    /* Interpolate and add to overall atmos.kappa */
    for(ll=0; ll<NTAU; ll++){
      Locate(NTEMP, opacNO.T, atmos.T[ll], &a);
      Locate(NPRESSURE, opacNO.P, atmos.P[ll], &b);
      for(i=0; i<NLAMBDA; i++){
    	atmos.kappa[i][ll] += lint2D(
    		opacNO.T[a], opacNO.T[a+1], opacNO.P[b], opacNO.P[b+1],
    		opacNO.kappa[i][b][a], opacNO.kappa[i][b][a+1],
    		opacNO.kappa[i][b+1][a], opacNO.kappa[i][b+1][a+1],
    		atmos.T[ll], atmos.P[ll])
    				     * atmos.NO[ll];
      }
    }

    FreeOpacTable(opacNO);
  }
  
  /* Fill in NO2 opacities */
  if(chemSelection[18] == 1){
    opacNO2.name = "NO2";
    opacNO2.T = dvector(0, NTEMP-1);
    opacNO2.P = dvector(0, NPRESSURE-1);
    opacNO2.Plog10 = dvector(0, NPRESSURE-1);
    opacNO2.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacNO2.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacNO2, fileArray[19]);
    
    printf("Read NO2 Opacity done\n");

    /* Interpolate and add to overall atmos.kappa */
    for(ll=0; ll<NTAU; ll++){
      Locate(NTEMP, opacNO2.T, atmos.T[ll], &a);
      Locate(NPRESSURE, opacNO2.P, atmos.P[ll], &b);
      for(i=0; i<NLAMBDA; i++){
    	atmos.kappa[i][ll] += lint2D(
    		opacNO2.T[a], opacNO2.T[a+1], opacNO2.P[b], opacNO2.P[b+1],
    		opacNO2.kappa[i][b][a], opacNO2.kappa[i][b][a+1],
    		opacNO2.kappa[i][b+1][a], opacNO2.kappa[i][b+1][a+1],
    		atmos.T[ll], atmos.P[ll])
    				     * atmos.NO2[ll];
      }
    }

    FreeOpacTable(opacNO2);
  }
  
  /* Fill in OCS opacities */
  if(chemSelection[19] == 1){
    opacOCS.name = "OCS";
    opacOCS.T = dvector(0, NTEMP-1);
    opacOCS.P = dvector(0, NPRESSURE-1);
    opacOCS.Plog10 = dvector(0, NPRESSURE-1);
    opacOCS.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacOCS.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacOCS, fileArray[22]);
    
    printf("Read OCS Opacity done\n");

    /* Interpolate and add to overall atmos.kappa */
    for(ll=0; ll<NTAU; ll++){
      Locate(NTEMP, opacOCS.T, atmos.T[ll], &a);
      Locate(NPRESSURE, opacOCS.P, atmos.P[ll], &b);
      for(i=0; i<NLAMBDA; i++){
    	atmos.kappa[i][ll] += lint2D(
    		opacOCS.T[a], opacOCS.T[a+1], opacOCS.P[b], opacOCS.P[b+1],
    		opacOCS.kappa[i][b][a], opacOCS.kappa[i][b][a+1],
    		opacOCS.kappa[i][b+1][a], opacOCS.kappa[i][b+1][a+1],
    		atmos.T[ll], atmos.P[ll])
    				     * atmos.OCS[ll];
      }
    }

    FreeOpacTable(opacOCS);
  }
  
  /* Fill in OH opacities */
  if(chemSelection[20] == 1){
    opacOH.name = "OH";
    opacOH.T = dvector(0, NTEMP-1);
    opacOH.P = dvector(0, NPRESSURE-1);
    opacOH.Plog10 = dvector(0, NPRESSURE-1);
    opacOH.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacOH.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacOH, fileArray[23]);
    
    printf("Read OH Opacity done\n");

    /* Interpolate and add to overall atmos.kappa */
    for(ll=0; ll<NTAU; ll++){
      Locate(NTEMP, opacOH.T, atmos.T[ll], &a);
      Locate(NPRESSURE, opacOH.P, atmos.P[ll], &b);
      for(i=0; i<NLAMBDA; i++){
    	atmos.kappa[i][ll] += lint2D(
    		opacOH.T[a], opacOH.T[a+1], opacOH.P[b], opacOH.P[b+1],
    		opacOH.kappa[i][b][a], opacOH.kappa[i][b][a+1],
    		opacOH.kappa[i][b+1][a], opacOH.kappa[i][b+1][a+1],
    		atmos.T[ll], atmos.P[ll])
    				     * atmos.OH[ll];
      }
    }

    FreeOpacTable(opacOH);
  }
  
  /* Fill in PH3 opacities */
  if(chemSelection[21] == 1){
    opacPH3.name = "PH3";
    opacPH3.T = dvector(0, NTEMP-1);
    opacPH3.P = dvector(0, NPRESSURE-1);
    opacPH3.Plog10 = dvector(0, NPRESSURE-1);
    opacPH3.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacPH3.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacPH3, fileArray[24]);
    
    printf("Read PH3 Opacity done\n");

    FreeOpacTable(opacPH3);
  }
  
  /* Fill in SH opacities */
  if(chemSelection[22] == 1){
    opacSH.name = "SH";
    opacSH.T = dvector(0, NTEMP-1);
    opacSH.P = dvector(0, NPRESSURE-1);
    opacSH.Plog10 = dvector(0, NPRESSURE-1);
    opacSH.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacSH.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacSH, fileArray[25]);
    
    printf("Read SH Opacity done\n");

    FreeOpacTable(opacSH);
  }
  
  /* Fill in SiH opacities */
  if(chemSelection[23] == 1){
    opacSiH.name = "SiH";
    opacSiH.T = dvector(0, NTEMP-1);
    opacSiH.P = dvector(0, NPRESSURE-1);
    opacSiH.Plog10 = dvector(0, NPRESSURE-1);
    opacSiH.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacSiH.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacSiH, fileArray[27]);
    
    printf("Read SiH Opacity done\n");

    FreeOpacTable(opacSiH); 
  }
  
  /* Fill in SiO opacities */
  if(chemSelection[24] == 1){
    opacSiO.name = "SiO";
    opacSiO.T = dvector(0, NTEMP-1);
    opacSiO.P = dvector(0, NPRESSURE-1);
    opacSiO.Plog10 = dvector(0, NPRESSURE-1);
    opacSiO.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacSiO.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacSiO, fileArray[28]);
    
    printf("Read SiO Opacity done\n");

    FreeOpacTable(opacSiO);
  }
  
  /* Fill in SO2 opacities */
  if(chemSelection[25] == 1){
    opacSO2.name = "SO2";
    opacSO2.T = dvector(0, NTEMP-1);
    opacSO2.P = dvector(0, NPRESSURE-1);
    opacSO2.Plog10 = dvector(0, NPRESSURE-1);
    opacSO2.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacSO2.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacSO2, fileArray[26]);
    
    printf("Read SO2 Opacity done\n");

    /* Interpolate and add to overall atmos.kappa */
    for(ll=0; ll<NTAU; ll++){
      Locate(NTEMP, opacSO2.T, atmos.T[ll], &a);
      Locate(NPRESSURE, opacSO2.P, atmos.P[ll], &b);
      for(i=0; i<NLAMBDA; i++){
    	atmos.kappa[i][ll] += lint2D(
    		opacSO2.T[a], opacSO2.T[a+1], opacSO2.P[b], opacSO2.P[b+1],
    		opacSO2.kappa[i][b][a], opacSO2.kappa[i][b][a+1],
    		opacSO2.kappa[i][b+1][a], opacSO2.kappa[i][b+1][a+1],
    		atmos.T[ll], atmos.P[ll])
    				     * atmos.SO2[ll];
      }
    }

    FreeOpacTable(opacSO2);
  }
  
  /* Fill in TiO opacities */
  if(chemSelection[26] == 1){
    opacTiO.name = "TiO";
    opacTiO.T = dvector(0, NTEMP-1);
    opacTiO.P = dvector(0, NPRESSURE-1);
    opacTiO.Plog10 = dvector(0, NPRESSURE-1);
    opacTiO.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacTiO.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacTiO, fileArray[29]);
    
    printf("Read TiO Opacity done\n");

    FreeOpacTable(opacTiO);
  }
  
  /* Fill in VO opacities */
  if(chemSelection[27] == 1){
    opacVO.name = "VO";
    opacVO.T = dvector(0, NTEMP-1);
    opacVO.P = dvector(0, NPRESSURE-1);
    opacVO.Plog10 = dvector(0, NPRESSURE-1);
    opacVO.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacVO.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacVO, fileArray[30]);
    
    printf("Read VO Opacity done\n");
    
    FreeOpacTable(opacVO);
  }

  /* Atomic opacities */

  /* Fill in Na opacities */
  if(chemSelection[28] == 1){
    opacNa.name = "Na";
    opacNa.T = dvector(0, NTEMP-1);
    opacNa.P = dvector(0, NPRESSURE-1);
    opacNa.Plog10 = dvector(0, NPRESSURE-1);
    opacNa.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacNa.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacNa, fileArray[31]);
    
    printf("Read Na Opacity done\n");
    
    FreeOpacTable(opacNa);
  }

  /* Fill in K opacities */
  if(chemSelection[29] == 1){
    opacK.name = "K";
    opacK.T = dvector(0, NTEMP-1);
    opacK.P = dvector(0, NPRESSURE-1);
    opacK.Plog10 = dvector(0, NPRESSURE-1);
    opacK.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacK.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    ReadOpacTable(opacK, fileArray[32]);
    
    printf("Read K Opacity done\n");
    
    FreeOpacTable(opacK);
  }
  
  /* Fill in collision-induced opacities */
  if(chemSelection[31] == 1){	      //If CIA is turned on ...
    
    /* Allocate collison induced opacities */
    opacCIA.name = "CIA";
    opacCIA.T = dvector(0, NTEMP-1);
    opacCIA.P = dvector(0, NPRESSURE-1);
    opacCIA.Plog10 = dvector(0, NPRESSURE-1);
    opacCIA.kappa = d3tensor(0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);
    opacCIA.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    
    /* populate with zeros */	
    for (i=0; i<NLAMBDA; i++)
      for (j=0; j<NPRESSURE; j++)
	for (k=0; k<NTEMP; k++)
	  opacCIA.kappa[i][j][k] = 0.;
    
    /* Read in CIA opacities */
    
    f1 = fopen(fileArray[33], "r");
    if(f1 == NULL){
      printf("\n totalopac.c:\nError opening file: %s -- No such file or directory\n\n", fileArray[33]);
      exit(1);
    }
    
    for(i=0; i<NTEMP; i++){
      fscanf(f1, "%le", &opacCIA.T[i]);
    }
    for (k=0; k<NTEMP; k++){
      fscanf(f1, "%le", &opacCIA.T[k]);
      for (i=0; i<NLAMBDA; i++){
	fscanf(f1, "%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le", &atmos.lambda[i], &opac_CIA_H2H2[k][i], &opac_CIA_H2He[k][i], &opac_CIA_H2H[k][i], &opac_CIA_H2CH4[k][i], &opac_CIA_CH4Ar[k][i], &opac_CIA_CH4CH4[k][i], &opac_CIA_CO2CO2[k][i], &opac_CIA_HeH[k][i], &opac_CIA_N2CH4[k][i], &opac_CIA_N2H2[k][i], &opac_CIA_N2N2[k][i], &opac_CIA_O2CO2[k][i], &opac_CIA_O2N2[k][i], &opac_CIA_O2O2[k][i]);
      }
    }

    for (i=0; i<NLAMBDA; i++){
      for(ll=0; ll<NTAU; ll++){

	Locate(NTEMP, opacCIA.T, atmos.T[ll], &a);

	CIA_temp = lint(opacCIA.T[a], opac_CIA_H2H2[a][i], 
			opacCIA.T[a+1], opac_CIA_H2H2[a+1][i], atmos.T[ll]); 
	atmos.kappa[i][ll] += CIA_temp *
	  atmos.H2[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]) *
	  atmos.H2[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]);

	CIA_temp = lint(opacCIA.T[a], opac_CIA_H2H[a][i], 
			opacCIA.T[a+1], opac_CIA_H2H[a+1][i], atmos.T[ll]); 
	atmos.kappa[i][ll] += CIA_temp *
	  atmos.H2[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]) *
	  atmos.H[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]);

	CIA_temp = lint(opacCIA.T[a], opac_CIA_H2CH4[a][i], 
			opacCIA.T[a+1], opac_CIA_H2CH4[a+1][i], atmos.T[ll]); 
	atmos.kappa[i][ll] += CIA_temp *
	  atmos.H2[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]) *
	  atmos.CH4[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]);

	/* atmos.kappa[i][ll] += 0.0 * */
	/*   atmos.CH4[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]) *  */
	/*   atmos.Ar[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]); */

	CIA_temp = lint(opacCIA.T[a], opac_CIA_CH4CH4[a][i], 
			opacCIA.T[a+1], opac_CIA_CH4CH4[a+1][i], atmos.T[ll]);
	atmos.kappa[i][ll] += CIA_temp *
	  atmos.CH4[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]) *
	  atmos.CH4[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]);

	CIA_temp = lint(opacCIA.T[a], opac_CIA_CO2CO2[a][i], 
			opacCIA.T[a+1], opac_CIA_CO2CO2[a+1][i], atmos.T[ll]);
	atmos.kappa[i][ll] += CIA_temp *
	  atmos.CO2[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]) *
	  atmos.CO2[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]);

	/* atmos.kappa += 0.0 * */
	/*   atmos.He[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]) *  */
	/*   atmos.H[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]); */

	CIA_temp = lint(opacCIA.T[a], opac_CIA_N2CH4[a][i], 
			opacCIA.T[a+1], opac_CIA_N2CH4[a+1][i], atmos.T[ll]);
	atmos.kappa[i][ll] += CIA_temp *
	  atmos.N2[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]) *
	  atmos.CH4[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]);

	CIA_temp = lint(opacCIA.T[a], opac_CIA_N2H2[a][i], 
			opacCIA.T[a+1], opac_CIA_N2H2[a+1][i], atmos.T[ll]);
	atmos.kappa[i][ll] += CIA_temp *
	  atmos.N2[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]) *
	  atmos.H2[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]);
	
	CIA_temp = lint(opacCIA.T[a], opac_CIA_N2N2[a][i], 
			opacCIA.T[a+1], opac_CIA_N2N2[a+1][i], atmos.T[ll]);
	atmos.kappa[i][ll] += CIA_temp *
	  atmos.N2[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]) *
	  atmos.N2[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]);

	CIA_temp = lint(opacCIA.T[a], opac_CIA_O2CO2[a][i], 
			opacCIA.T[a+1], opac_CIA_O2CO2[a+1][i], atmos.T[ll]);
	atmos.kappa[i][ll] += CIA_temp *
	  atmos.O2[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]) *
	  atmos.CO2[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]);

	CIA_temp = lint(opacCIA.T[a], opac_CIA_O2N2[a][i], 
			opacCIA.T[a+1], opac_CIA_O2N2[a+1][i], atmos.T[ll]);
	atmos.kappa[i][ll] += CIA_temp *
	  atmos.O2[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]) *
	  atmos.N2[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]);
	
	CIA_temp = lint(opacCIA.T[a], opac_CIA_O2O2[a][i], 
			opacCIA.T[a+1], opac_CIA_O2O2[a+1][i], atmos.T[ll]);
	atmos.kappa[i][ll] += CIA_temp *
	  atmos.O2[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]) *
	  atmos.O2[ll] * atmos.P[ll] / (KBOLTZMANN * atmos.T[ll]);
      }
    }

  }
 
  /* Rayleigh scattering */
  
  /* (Polarizabilities from the CRC Handbook) */
  
  if(chemSelection[30] == 1){		//If Scattering is activated.. 
    
    /* Fill in scattering coefficients */
    for (i=0; i<NLAMBDA; i++) {
      for (ll=0; ll<NTAU; ll++) {
	atmos.kappa[i][ll] += RAYLEIGH * (
	    (8.0*PI/3.0) * SQ(0.80e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    atmos.H2[ll]*atmos.P[ll] / (KBOLTZMANN * atmos.T[ll])
	    +
	    (8.0*PI/3.0) * SQ(1.74e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    atmos.N2[ll]*atmos.P[ll] / (KBOLTZMANN * atmos.T[ll])
	    +
	    (8.0*PI/3.0) * SQ(1.45e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    atmos.H2O[ll]*atmos.P[ll] / (KBOLTZMANN * atmos.T[ll])
	    +
	    (8.0*PI/3.0) * SQ(1.95e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    atmos.CO[ll]*atmos.P[ll] / (KBOLTZMANN * atmos.T[ll])
	    +
	    (8.0*PI/3.0) * SQ(2.91e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    atmos.CO2[ll]*atmos.P[ll] / (KBOLTZMANN * atmos.T[ll])
	    +
	    (8.0*PI/3.0) * SQ(2.26e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    atmos.NH3[ll]*atmos.P[ll] / (KBOLTZMANN * atmos.T[ll])
	    +
	    (8.0*PI/3.0) * SQ(2.59e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    atmos.CH4[ll]*atmos.P[ll] / (KBOLTZMANN * atmos.T[ll])
	    +
	    (8.0*PI/3.0) * SQ(1.58e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    atmos.O2[ll]*atmos.P[ll] / (KBOLTZMANN * atmos.T[ll])
	    +
	    (8.0*PI/3.0) * SQ(3.21e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    atmos.O3[ll]*atmos.P[ll] / (KBOLTZMANN * atmos.T[ll])
	    +
	    (8.0*PI/3.0) * SQ(3.33e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    atmos.C2H2[ll]*atmos.P[ll] / (KBOLTZMANN * atmos.T[ll])
	    +
	    (8.0*PI/3.0) * SQ(4.25e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    atmos.C2H4[ll]*atmos.P[ll] / (KBOLTZMANN * atmos.T[ll])
	    +
	    (8.0*PI/3.0) * SQ(4.47e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    atmos.C2H6[ll]*atmos.P[ll] / (KBOLTZMANN * atmos.T[ll])
	    +
	    (8.0*PI/3.0) * SQ(2.59e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    atmos.HCN[ll]*atmos.P[ll] / (KBOLTZMANN * atmos.T[ll])
	    +
	    (8.0*PI/3.0) * SQ(3.78e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    atmos.H2S[ll]*atmos.P[ll] / (KBOLTZMANN * atmos.T[ll])
	    +
	    (8.0*PI/3.0) * SQ(1.70e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    atmos.NO[ll]*atmos.P[ll] / (KBOLTZMANN * atmos.T[ll])
	    +
	    (8.0*PI/3.0) * SQ(3.02e-30) *
	    SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	    atmos.NO2[ll]*atmos.P[ll] / (KBOLTZMANN * atmos.T[ll])
				   );
      }
    }
      
  }

  /*  Free memory */

  free_dmatrix(opac_CIA_H2H2, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_H2He, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_H2H, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_H2CH4, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_CH4Ar, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_CH4CH4, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_CO2CO2, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_HeH, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_N2CH4, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_N2H2, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_N2N2, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_O2CO2, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_O2N2, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_O2O2, 0, NTEMP-1, 0, NLAMBDA-1);

  printf("\n\n***********\nWARNING:  Double check that T-P profile lies within T and P range covered by the opacity data tables!!\n***********\n\n");
  
}

/* ------- end -------------- TotalOpac.c ------------------------ */

