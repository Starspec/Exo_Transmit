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

/*----------------------- readchemtable.c ------------------------

Author: Eliza Kempton (kemptone@grinnell.edu)

------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>

#include "opac.h" 
#include "nrutil.h" 
#include "vars.h"
#include "prototypes.h"
 
/* --- Global variables ------------------------------------------ */

extern struct Chem chem;

/* ---------------------------------------------------------------
 * Read in chemistry files: abundance(pressure, temperature)
 * --------------------------------------------------------------- */

/* ------- begin ------------ ReadChemTable.c -------------------- */

void ReadChemTable() {
  
  int i;//, j, k;
  char dum[8];
  int chemSelection[30];
  
  /* Initialize and obtain variables from other files */
  char **fileArray = getFileArray();
  vars variables = getVars();
  getChemSelection(chemSelection);
  int NTAU = variables.NTAU;
  
  printf("NTAU = %d\n", NTAU);
  
  FILE *f1;
  
  /* Allocate memory for Chem structure */
  
  chem.T = dvector(0, NTAU-1);	
  chem.P = dvector(0, NTAU-1);
  
  chem.total = dvector(0, NTAU-1);
  chem.C = dvector(0, NTAU-1);	
  chem.CH4 = dvector(0, NTAU-1);
  chem.CO = dvector(0, NTAU-1);
  chem.CO2 = dvector(0, NTAU-1);
  chem.C2H2 = dvector(0, NTAU-1);
  chem.C2H4 = dvector(0, NTAU-1);	
  chem.C2H6 = dvector(0, NTAU-1);		
  chem.H = dvector(0, NTAU-1);	
  chem.HCN = dvector(0, NTAU-1);	
  chem.HCl = dvector(0, NTAU-1);	
  chem.HF = dvector(0, NTAU-1);	
  chem.H2 = dvector(0, NTAU-1);	
  chem.H2CO = dvector(0, NTAU-1);	
  chem.H2O = dvector(0, NTAU-1);	
  chem.H2S = dvector(0, NTAU-1);	
  chem.He = dvector(0, NTAU-1);	
  chem.K = dvector(0, NTAU-1);
  chem.MgH = dvector(0, NTAU-1);	
  chem.N = dvector(0, NTAU-1);	
  chem.N2 = dvector(0, NTAU-1);	
  chem.NO2 = dvector(0, NTAU-1);	
  chem.NH3 = dvector(0, NTAU-1);	
  chem.NO = dvector(0, NTAU-1);	
  chem.Na = dvector(0, NTAU-1);
  chem.O = dvector(0, NTAU-1);	
  chem.O2 = dvector(0, NTAU-1);	
  chem.O3 = dvector(0, NTAU-1);	
  chem.OCS = dvector(0, NTAU-1); 	
  chem.OH = dvector(0, NTAU-1);	
  chem.PH3 = dvector(0, NTAU-1);	
  chem.SH = dvector(0, NTAU-1);	
  chem.SO2 = dvector(0, NTAU-1);	
  chem.SiH = dvector(0, NTAU-1);	
  chem.SiO = dvector(0, NTAU-1);	
  chem.TiO = dvector(0, NTAU-1);	
  chem.VO = dvector(0, NTAU-1);	
  
  /* Read in chemistry abundances */	
  
  f1 = fopen(fileArray[1],"r");
  if(f1 == NULL){
    printf("\nreadchemtable.c:\nError opening file: %s \nNo such file or directory.\nMake sure you have the appropriate file path and name specified in userInput.in\n\n", fileArray[1]);
    exit(1);
  }
  
  {
    /* Skip the first line */
    // This must be longer than the first line, so making it unreasonably long
    int _line = 1000;
    char _buffer[_line]; 
    fgets(_buffer, _line, f1);
    printf("Skipped first line with a buffer.\n");
  };
  
  for (i=NTAU-1; i>=0; i--){
    fscanf(f1,"%le", &chem.P[i]);
      if(i == NTAU-1){
      	fscanf(f1,"%le", &chem.T[i]); 
      	fscanf(f1,"%le", &chem.total[i]);
      	fscanf(f1,"%le", &chem.C[i]);
      	fscanf(f1,"%le", &chem.CH4[i]); 

      	/* Checks for errors, but only the first time */
      	errorCheck(chemSelection[0], chem.CH4[i]);
      	fscanf(f1,"%le", &chem.CO[i]);  
      	errorCheck(chemSelection[2], chem.CO[i]);
      	fscanf(f1,"%le", &chem.OCS[i]); 
      	errorCheck(chemSelection[19], chem.OCS[i]);
      	fscanf(f1,"%le", &chem.CO2[i]); 
      	errorCheck(chemSelection[1], chem.CO2[i]);
      	fscanf(f1,"%le", &chem.C2H2[i]); 
      	errorCheck(chemSelection[7], chem.C2H2[i]);
      	fscanf(f1,"%le", &chem.C2H4[i]); 
      	errorCheck(chemSelection[8], chem.C2H4[i]);
      	fscanf(f1,"%le", &chem.C2H6[i]); 
      	errorCheck(chemSelection[8], chem.C2H6[i]);
      	fscanf(f1,"%le", &chem.H[i]);
      	fscanf(f1,"%le", &chem.HCN[i]); 
      	errorCheck(chemSelection[13], chem.HCN[i]);
      	fscanf(f1,"%le", &chem.HCl[i]); 
      	errorCheck(chemSelection[12], chem.HCl[i]);
      	fscanf(f1,"%le", &chem.HF[i]); 
      	errorCheck(chemSelection[14], chem.HF[i]);
      	fscanf(f1,"%le", &chem.H2[i]);
      	fscanf(f1,"%le", &chem.H2CO[i]); 
      	errorCheck(chemSelection[10], chem.H2CO[i]);
      	fscanf(f1,"%le", &chem.H2O[i]); 
      	errorCheck(chemSelection[3], chem.H2O[i]);
      	fscanf(f1,"%le", &chem.H2S[i]); 
      	errorCheck(chemSelection[11], chem.H2S[i]);
      	fscanf(f1,"%le", &chem.He[i]);
      	fscanf(f1,"%le", &chem.K[i]);
      	fscanf(f1,"%le", &chem.MgH[i]); 
      	errorCheck(chemSelection[15], chem.MgH[i]);
      	fscanf(f1,"%le", &chem.N[i]);
      	fscanf(f1,"%le", &chem.N2[i]); 
      	errorCheck(chemSelection[16], chem.N2[i]);
      	fscanf(f1,"%le", &chem.NO2[i]); 
      	errorCheck(chemSelection[18], chem.NO2[i]); 
      	fscanf(f1,"%le", &chem.NH3[i]); 
      	errorCheck(chemSelection[4], chem.NH3[i]);
      	fscanf(f1,"%le", &chem.NO[i]); 
      	errorCheck(chemSelection[17], chem.NO[i]);
      	fscanf(f1,"%le", &chem.Na[i]); 
      	fscanf(f1,"%le", &chem.O[i]);
      	fscanf(f1,"%le", &chem.O2[i]); 
      	errorCheck(chemSelection[5], chem.O2[i]);
      	fscanf(f1,"%le", &chem.O3[i]); 
      	errorCheck(chemSelection[5], chem.O3[i]);
      	fscanf(f1,"%le", &chem.OH[i]); 
      	errorCheck(chemSelection[20], chem.OH[i]);
      	fscanf(f1,"%le", &chem.PH3[i]); 
      	errorCheck(chemSelection[21], chem.PH3[i]); 
      	fscanf(f1,"%le", &chem.SH[i]); 
      	errorCheck(chemSelection[22], chem.SH[i]);
      	fscanf(f1,"%le", &chem.SO2[i]); 
      	errorCheck(chemSelection[25], chem.SO2[i]); 
      	fscanf(f1,"%le", &chem.SiH[i]); 
      	errorCheck(chemSelection[23], chem.SiH[i]); 
      	fscanf(f1,"%le", &chem.SiO[i]); 
      	errorCheck(chemSelection[24], chem.SiO[i]); 
      	fscanf(f1,"%le", &chem.TiO[i]); 
      	errorCheck(chemSelection[26], chem.TiO[i]); 
      	fscanf(f1,"%le", &chem.VO[i]); 
      	errorCheck(chemSelection[27], chem.VO[i]);
            }
            else{
      	fscanf(f1,"%le", &chem.T[i]); 
      	fscanf(f1,"%le", &chem.total[i]);
      	fscanf(f1,"%le", &chem.C[i]);
      	fscanf(f1,"%le", &chem.CH4[i]); 
      	fscanf(f1,"%le", &chem.CO[i]); 
      	fscanf(f1,"%le", &chem.OCS[i]); 
      	fscanf(f1,"%le", &chem.CO2[i]); 
      	fscanf(f1,"%le", &chem.C2H2[i]); 
      	fscanf(f1,"%le", &chem.C2H4[i]); 
      	fscanf(f1,"%le", &chem.C2H6[i]); 
      	fscanf(f1,"%le", &chem.H[i]);
      	fscanf(f1,"%le", &chem.HCN[i]); 
      	fscanf(f1,"%le", &chem.HCl[i]); 
      	fscanf(f1,"%le", &chem.HF[i]); 
      	fscanf(f1,"%le", &chem.H2[i]);
      	fscanf(f1,"%le", &chem.H2CO[i]); 
      	fscanf(f1,"%le", &chem.H2O[i]); 
      	fscanf(f1,"%le", &chem.H2S[i]); 
      	fscanf(f1,"%le", &chem.He[i]);
      	fscanf(f1,"%le", &chem.K[i]);
      	fscanf(f1,"%le", &chem.MgH[i]); 
      	fscanf(f1,"%le", &chem.N[i]);
      	fscanf(f1,"%le", &chem.N2[i]); 
      	fscanf(f1,"%le", &chem.NO2[i]); 
      	fscanf(f1,"%le", &chem.NH3[i]); 
      	fscanf(f1,"%le", &chem.NO[i]); 
      	fscanf(f1,"%le", &chem.Na[i]); 
      	fscanf(f1,"%le", &chem.O[i]);
      	fscanf(f1,"%le", &chem.O2[i]); 
      	fscanf(f1,"%le", &chem.O3[i]); 
      	fscanf(f1,"%le", &chem.OH[i]); 
      	fscanf(f1,"%le", &chem.PH3[i]); 
      	fscanf(f1,"%le", &chem.SH[i]); 
      	fscanf(f1,"%le", &chem.SO2[i]); 
      	fscanf(f1,"%le", &chem.SiH[i]); 
      	fscanf(f1,"%le", &chem.SiO[i]); 
      	fscanf(f1,"%le", &chem.TiO[i]); 
      	fscanf(f1,"%le", &chem.VO[i]); 
            }
          }
        
        fclose(f1);
        
        /* Print out final line of data as a double-check */
        printf("Chemistry: \n");	
        printf("P_0\t%e \n", chem.P[NTAU-1]);
        printf("T_0\t%e \n", chem.T[NTAU-1]);
        printf("total\t%e \n", chem.total[0]);
        printf("C\t%e \n", chem.C[0]);
        printf("CH4\t%e \n", chem.CH4[0]);
        printf("CO\t%e \n", chem.CO[0]);
        printf("CO2\t%e \n", chem.CO2[0]);
        printf("C2H2\t%e \n", chem.C2H2[0]);
        printf("C2H4\t%e \n", chem.C2H4[0]);
        printf("H\t%e \n", chem.H[0]);
        printf("HCN\t%e \n", chem.HCN[0]);
        printf("HCl\t%e \n", chem.HCl[0]);
        printf("HF\t%e \n", chem.HF[0]);
        printf("H2\t%e \n", chem.H2[0]);
        printf("H2CO\t%e \n", chem.H2CO[0]);
        printf("H2O\t%e \n", chem.H2O[0]);
        printf("H2S\t%e \n", chem.H2S[0]);
        printf("He\t%e \n", chem.He[0]);
        printf("MgH\t%e \n", chem.MgH[0]);
        printf("N\t%e \n", chem.N[0]);
        printf("N2\t%e \n", chem.N2[0]);
        printf("NO2\t%e \n", chem.NO2[0]);
        printf("NH3\t%e \n", chem.NH3[0]);
        printf("NO\t%e \n", chem.NO[0]);
        printf("O\t%e \n", chem.O[0]);
        printf("O2\t%e \n", chem.O2[0]);
        printf("OCS\t%e \n", chem.OCS[0]);
        printf("OH\t%e \n", chem.OH[0]);
        printf("PH3\t%e \n", chem.PH3[0]);
        printf("SH\t%e \n", chem.SH[0]);
        printf("SO2\t%e \n", chem.SO2[0]);
        printf("SiH\t%e \n", chem.SiH[0]);
        printf("SiO\t%e \n", chem.SiO[0]);
        printf("TiO\t%e \n", chem.TiO[0]);
        printf("VO\t%e \n", chem.VO[0]);
        printf("O3\t%e \n", chem.O3[0]);
        printf("C2H6\t%e \n", chem.C2H6[0]);
        printf("Na\t%e \n", chem.Na[0]);
        printf("K\t%e \n", chem.K[0]);
  
  return;
}

/* ------- end -------------- ReadChemTable.c -------------------- */

/* ------- start ------------ FreeChemTable.c -------------------- */

void FreeChemTable(){ 

  vars variables = getVars();
  int NPRESSURE = variables.NPRESSURE;
  int NTEMP = variables.NTEMP;
  int NTAU = variables.NTAU;

  free_dvector(chem.T, 0, NTAU-1);	
  free_dvector(chem.P, 0, NTAU-1);
  
  free_dvector(chem.total, 0, NTAU-1);
  free_dvector(chem.C, 0, NTAU-1);	
  free_dvector(chem.CH4, 0, NTAU-1);
  free_dvector(chem.CO, 0, NTAU-1);
  free_dvector(chem.CO2, 0, NTAU-1);
  free_dvector(chem.C2H2, 0, NTAU-1);
  free_dvector(chem.C2H4, 0, NTAU-1);	
  free_dvector(chem.C2H6, 0, NTAU-1);		
  free_dvector(chem.H, 0, NTAU-1);	
  free_dvector(chem.HCN, 0, NTAU-1);	
  free_dvector(chem.HCl, 0, NTAU-1);	
  free_dvector(chem.HF, 0, NTAU-1);	
  free_dvector(chem.H2, 0, NTAU-1);	
  free_dvector(chem.H2CO, 0, NTAU-1);	
  free_dvector(chem.H2O, 0, NTAU-1);	
  free_dvector(chem.H2S, 0, NTAU-1);	
  free_dvector(chem.He, 0, NTAU-1);	
  free_dvector(chem.K, 0, NTAU-1);
  free_dvector(chem.MgH, 0, NTAU-1);	
  free_dvector(chem.N, 0, NTAU-1);	
  free_dvector(chem.N2, 0, NTAU-1);	
  free_dvector(chem.NO2, 0, NTAU-1);	
  free_dvector(chem.NH3, 0, NTAU-1);	
  free_dvector(chem.NO, 0, NTAU-1);	
  free_dvector(chem.Na, 0, NTAU-1);
  free_dvector(chem.O, 0, NTAU-1);	
  free_dvector(chem.O2, 0, NTAU-1);	
  free_dvector(chem.O3, 0, NTAU-1);	
  free_dvector(chem.OCS, 0, NTAU-1); 	
  free_dvector(chem.OH, 0, NTAU-1);	
  free_dvector(chem.PH3, 0, NTAU-1);	
  free_dvector(chem.SH, 0, NTAU-1);	
  free_dvector(chem.SO2, 0, NTAU-1);	
  free_dvector(chem.SiH, 0, NTAU-1);	
  free_dvector(chem.SiO, 0, NTAU-1);	
  free_dvector(chem.TiO, 0, NTAU-1);	
  free_dvector(chem.VO, 0, NTAU-1);	
  free_dvector(chem.mu, 0, NTAU-1);	

}

/* ------- end -------------- FreeChemTable.c -------------------- */

/* ------- start ------------ errorCheck.c ----------------------- */

/* Shows error message if gas not selected is present in chem table. */
void errorCheck(int onOff, double value){
  if(onOff == 0 && value != 0.0){
    printf("Error: Gas not selected for inclusion in Calculation has non-zero Chem table entry\n");
  }
}

/* ------- end -------------- errorCheck.c ----------------------- */
