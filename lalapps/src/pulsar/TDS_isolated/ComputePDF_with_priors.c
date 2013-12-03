/*
*  Copyright (C) 2007 Matt Pitkin,  Rejean Dupuis
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/**
 * \file
 * \ingroup pulsarApps
 * \author Réjean Dupuis
 * \brief
 * Driver code to calculate posterior pdfs for known pulsars
 */

/*********************************************************************************/
/*                    Driver code to calculate posterior pdfs for known pulsars  */
/*                                                                               */
/*			               Réjean Dupuis                             */
/*                                                                               */
/*                       University of Glasgow - last modified 26/03/2004        */
/*********************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <string.h>
#include "FitToPulsar.h"  /* use this local file until up-to-date version commited to lal */

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 1; } else ((void)0)

#define MAXLENGTH 330000

/* PSIMEAN is 125.155-90 = 35.155 degrees = 0.61357rads */
#define PSIMEAN 0.6135704985386066
#define PSISIGMA 0.02015855286053451
/* IOTAMEAN is 62.165 degrees = 1.08498 rads */
#define IOTAMEAN 1.084983929502275
#define IOTASIGMA 0.02321287905152458

int main(int argc, char **argv)
{
  /***** declare variables *****/
  static LALStatus      status;
  FILE *fp1, *fp2, *fp3, *fp4;
  FILE *fp_pdf1, *fp_pdf2, *fp_pdf3, *fp_pdf4;
  FILE *fp_pdf1Phase, *fp_pdf2Phase, *fp_pdf3Phase, *fp_pdf4Phase;
  FILE *fp_pdf1Psi, *fp_pdf2Psi, *fp_pdf3Psi, *fp_pdf4Psi;
  FILE *fp_pdf1CosIota, *fp_pdf2CosIota, *fp_pdf3CosIota, *fp_pdf4CosIota;  
  FILE *fp_joint, *fp_jointPhase, *fp_jointPsi, *fp_jointCosIota;
  FILE *fpmesh, *psrfp;
  char infile1[256], infile2[256], infile3[256], infile4[256];  
  char outfile1[256], outfile2[256], outfile3[256], outfile4[256];
  char outfile1Phase[256], outfile2Phase[256], outfile3Phase[256], outfile4Phase[256];
  char outfile1Psi[256], outfile2Psi[256], outfile3Psi[256], outfile4Psi[256];
  char outfile1CosIota[256], outfile2CosIota[256], outfile3CosIota[256], outfile4CosIota[256];   
  char outfile[256], outfilePhase[256],outfilePsi[256],outfileCosIota[256];   
  char psrinput[256], inmesh[256], txt[32], pulsar_name[64];
  UINT4 iPsi, iPhase, iCosIota, iH0, arg, i;
  INT4 iRA=0, iDEC=0, if0=0, if1=0, if2=0,ifepoch=0, flag =0, irun;
  INT4 num_dect, iGEO = 0, iL1 = 0, iH1 = 0, iH2 = 0;
  REAL4 psr_ra, psr_dec; 
  REAL8 RA, DEC, f0, f1, f2, val, t;
  REAL8 minChi, area, outL1, outH1, outH2, outGEO;
  COMPLEX16 B, var;    
  LIGOTimeGPS fepoch, tgps1[MAXLENGTH],tgps2[MAXLENGTH],tgps3[MAXLENGTH], tgps4[MAXLENGTH];    
  FitInputStudentT  input1, input2, input3, input4;
  CoarseFitInput input1_chi, input2_chi, input3_chi, input4_chi;
  CoarseFitParams params1, params2, params3, params4;
  CoarseFitOutput output1, output2, output3, output4;
  LALDetector 	detectorLHO,detectorLLO, detectorGEO;
  LALSource pulsar;
  PulsarPdfs prob1, prob2, prob3, prob4, prob;


/* get command line arguments and print error if there is wrong number of arguments.  if 
the variable condor is defined, then read in the first argument pulsar_name from stdin */

#ifndef CONDOR 
  if (argc < 8 || argc > 11){
   fprintf(stderr, "1. Name of pulsar (when not using condor)\n"); 
   fprintf(stderr,"2. working directory\n");
   fprintf(stderr, "3. flag=1 for Gaussian with explicit noise estimate provided\n flag=2 for noise analytically marginalised out of Gaussian\n");
   fprintf(stderr, "4. name of mesh file in inputs directory\n");
   fprintf(stderr, "5. irun - 2 for S2, 3 for S3\n");
   fprintf(stderr, "6. num detectors (max 4)\n");
   fprintf(stderr, "7-10. H1, H2, L1, and/or GEO\n");
   fprintf(stderr,"There were argc= %d arguments:\n",argc -1);
   fflush(stderr);
   return 2;
  }
  sprintf(pulsar_name, "%s", argv[1]);
#else
  // read pulsar name from stdin for condor
   scanf("%s", &pulsar_name[0]);
   fprintf(stderr, "pulsar name: %s\n",pulsar_name);
#endif  

  fprintf(stderr, "pulsarname=%s\n", pulsar_name);
   
  flag = atoi(argv[3]);
  irun = atoi(argv[5]);
  num_dect = atoi(argv[6]);
  
  if (num_dect >4)
  {
    fprintf(stderr, "Number of IFOs must be less than 5 (H1, H2, L1, and/or GEO)!\n");
    return(1);
  }
  
  /* check which detectors we are being used for analysis*/
  for (i=0;i<num_dect;i++)
  {
    if (!strcmp(argv[i+7],"GEO")) iGEO = 1;
    else if (!strcmp(argv[i+7],"L1")) iL1 = 1;
    else if (!strcmp(argv[i+7],"H1")) iH1 = 1;
    else if (!strcmp(argv[i+7],"H2")) iH2 = 1;
    else 
    {
      fprintf(stderr, "what is %s? not GEO, L1, H1, or H2!\n",argv[i+7]);
      return(2); 
    }
  }
/************** BEGIN READING PULSAR PARAMETERS ******************/     
/* read input file with pulsar information */
  
  sprintf(psrinput,"%s/inputs/%s", argv[2],pulsar_name);
  psrfp=fopen(psrinput,"r");
  
  while (2==fscanf(psrfp,"%s %lf", &txt[0], &val))
  {
    if( !strcmp(txt,"ra") || !strcmp(txt,"RA")) {
      psr_ra = RA = val;
      iRA = 1;
    }  
    
    else if( !strcmp(txt,"dec") || !strcmp(txt,"DEC")) {
      psr_dec = DEC = val;
      iDEC = 1;
    }
    else if( !strcmp(txt,"f0") || !strcmp(txt,"F0")) {
      f0 = val;
      if0 = 1;
    }
    else if( !strcmp(txt,"f1") || !strcmp(txt,"F1")) {
      f1 = val;
      if1 = 1;
    }
    else if( !strcmp(txt,"f2") || !strcmp(txt,"F2")) {
      f2 = val;
      if2 = 1;
    }
    else if( !strcmp(txt,"fepoch") || !strcmp(txt,"fepoch")) {
      fepoch.gpsSeconds = floor(val*1e-9);
      fepoch.gpsNanoSeconds = (INT8) val - fepoch.gpsSeconds*1e9;
      ifepoch = 1;
    }
  }
  
  if ((iRA+iDEC+ifepoch+if0+if1+if2) != 6) {
    fprintf(stderr, "pulsar input file missing info \nra\txxx\ndec\txxx\nf0\txxx\nf1\txxx\nf2\txxx\nfepoch\txxx\n");
    return(3);
  }

  fclose(psrfp);
  
  /************** END READING PULSAR PARAMETERS ******************/   

  detectorLHO = lalCachedDetectors[LALDetectorIndexLHODIFF];
  detectorLLO = lalCachedDetectors[LALDetectorIndexLLODIFF];
  detectorGEO = lalCachedDetectors[LALDetectorIndexGEO600DIFF];

  /* currently this code assumes that we are using 30 minute stretches of data
     this could be relaxed in the future */
  input1.N = input2.N = input3.N = input4.N = 30;  
  input1_chi.N = input2_chi.N = input3_chi.N = input4_chi.N = 30;
  
  /* set up strings pointing to the file that we want to read as input
     for the gaussian likelihood (flag=1) there is two extra columns in the 
     input files corresponding to the variance of the real and imaginary Bk's. 
     the Bk's used in the gausian likelihood (flag=1) are calculated every 30 minutes
     (average of 30 1-minute Bk).
     for the student-t likelihood (flag=2), the Bk's are from every 60 seconds */
     
  if (flag == 1)
  {
   if (iL1) sprintf(infile1,"%s/dataL1/outfine.%s_L1.S%d_chi_30", argv[2],pulsar_name, irun);   
   if (iH1) sprintf(infile2,"%s/dataH1/outfine.%s_H1.S%d_chi_30", argv[2], pulsar_name, irun);   
   if (iH2) sprintf(infile3,"%s/dataH2/outfine.%s_H2.S%d_chi_30", argv[2], pulsar_name, irun);         
   if (iGEO) sprintf(infile4, "%s/dataGEO/outfine.%s_GEO.S%d_chi_30", argv[2], pulsar_name, irun);
  }
  else if (flag == 2)
  {
    if (iL1) sprintf(infile1,"%s/dataL1/finehet_%s_L1", argv[2], pulsar_name);
    if (iH1) sprintf(infile2,"%s/dataH1/finehet_%s_H1", argv[2], pulsar_name);
    if (iH2) sprintf(infile3,"%s/dataH2/finehet_%s_H2", argv[2], pulsar_name); 
    if (iGEO) sprintf(infile4,"%s/dataG1/finehet_%s_GEO", argv[2],
pulsar_name);  
  }
  else
  {
    fprintf(stderr, "flag should be 1 or 2\n");
    return(2);    
  }
  
  /* open output files */
  if (iL1) fp1 = fopen(infile1, "r"); 
  if (iH1) fp2 = fopen(infile2, "r"); 
  if (iH2) fp3 = fopen(infile3, "r");
  if (iGEO) fp4 = fopen(infile4, "r");

  /* allocate memory for B_k and var */
  if (flag ==1 && iL1 == 1) 
  {
    input1_chi.B = NULL;
    LALZCreateVector( &status, &input1_chi.B, MAXLENGTH);
    TESTSTATUS(&status);

    input1_chi.var = NULL; 
    LALZCreateVector( &status, &input1_chi.var, MAXLENGTH);   
    TESTSTATUS(&status);
  }
  else if (flag == 2 && iL1 == 1)
  {
    input1.B = NULL;
    LALZCreateVector( &status, &input1.B, MAXLENGTH);
    TESTSTATUS(&status); 
  }
  
  if (flag ==1 && iH1 == 1) 
  {
    input2_chi.B = NULL;
    LALZCreateVector( &status, &input2_chi.B, MAXLENGTH);
    TESTSTATUS(&status);
    
    input2_chi.var = NULL; 
    LALZCreateVector( &status, &input2_chi.var, MAXLENGTH); 
    TESTSTATUS(&status);   
  }
  else if (flag == 2 && iH1 == 1)
  {
    input2.B = NULL;
    LALZCreateVector( &status, &input2.B, MAXLENGTH); 
    TESTSTATUS(&status);
  }
  
  if (flag ==1 && iH2 == 1) 
  {
    input3_chi.B = NULL;
    LALZCreateVector( &status, &input3_chi.B, MAXLENGTH);
    TESTSTATUS(&status);
    
    input3_chi.var = NULL; 
    LALZCreateVector( &status, &input3_chi.var, MAXLENGTH);  
    TESTSTATUS(&status);  
  }
  else if (flag == 2 && iH2 == 1)
  {
    input3.B = NULL;
    LALZCreateVector( &status, &input3.B, MAXLENGTH); 
    TESTSTATUS(&status);
  }
  if (flag ==1 && iGEO == 1) 
  {
    input4_chi.B = NULL;
    LALZCreateVector( &status, &input4_chi.B, MAXLENGTH);
    TESTSTATUS(&status);
    
    input4_chi.var = NULL; 
    LALZCreateVector( &status, &input4_chi.var, MAXLENGTH);   
    TESTSTATUS(&status); 
  }
  else if (flag == 2 && iGEO == 1)
  {
    input4.B = NULL;
    LALZCreateVector( &status, &input4.B, MAXLENGTH); 
    TESTSTATUS(&status);
  }
  
       
  params1.detector = detectorLLO;
  params2.detector = detectorLHO;
  params3.detector = detectorLHO;
  params4.detector = detectorGEO;
  
  /* set up RA, DEC, and coordinate system */
  pulsar.equatorialCoords.longitude = psr_ra;	      
  pulsar.equatorialCoords.latitude = psr_dec;	      
  pulsar.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
  
  /* polarization angle will be redefined inside the fitting routine */
  pulsar.orientation = 0.0; 	      
 		      
  params1.pulsarSrc = params2.pulsarSrc = params3.pulsarSrc = params4.pulsarSrc= pulsar;  
 
  /****************** BEGIN READ INPUT DATA (Bk's) ***************************/
 
  /* read data from L1 */
  if (iL1)
  {
    if (flag == 2) /* student-t */
    {
      fscanf(fp1,"%lf\t%lf\t%lf",&t,&B.re, &B.im);  
      i=0;  
      while (!feof(fp1))
      {
        tgps1[i].gpsSeconds = (INT4)floor(t);
        tgps1[i].gpsNanoSeconds = (INT4)floor((fmod(t,1.0)*1.e9));       
        input1.B->data[i].re = B.re;
        input1.B->data[i].im = B.im;    
        fscanf(fp1,"%lf\t%lf\t%lf",&t,&B.re, &B.im);    
        i++;  
      } 
      input1.t = tgps1;
      input1.B->length = i;
    }
    else if (flag == 1) /* chisquare */
    {
      fscanf(fp1,"%lf\t%lf\t%lf\t%lf\t%lf",&t,&B.re, &B.im, &var.re, &var.im);  
      i=0;  
      while (!feof(fp1))
      {
        tgps1[i].gpsSeconds = (INT4)floor(t);
        tgps1[i].gpsNanoSeconds = (INT4)floor((fmod(t,1.0)*1.e9));       
        input1_chi.B->data[i].re = B.re;
        input1_chi.B->data[i].im = B.im; 
        input1_chi.var->data[i].re = var.re;
        input1_chi.var->data[i].im = var.im;   
        fscanf(fp1,"%lf\t%lf\t%lf\t%lf\t%lf",&t,&B.re, &B.im, &var.re, &var.im);          
        i++;  
      } 
      input1_chi.t = tgps1;
      input1_chi.B->length = i;
      input1_chi.var->length = i;
    }
    fclose(fp1);
  }
 
  fprintf(stderr, "I've read in the L1 data.\n");
 /* read data from H1 */
  if (iH1)
  {
    if (flag == 2) /* student-t */
    {
      fscanf(fp2,"%lf\t%lf\t%lf",&t,&B.re, &B.im);  
      i=0;  
      while (!feof(fp2))
      {
        tgps2[i].gpsSeconds = (INT4)floor(t);
        tgps2[i].gpsNanoSeconds = (INT4)floor((fmod(t,1.0)*1.e9));       
        input2.B->data[i].re = B.re;
        input2.B->data[i].im = B.im;    
        fscanf(fp2,"%lf\t%lf\t%lf",&t,&B.re, &B.im);    
        i++;  
      } 
      input2.t = tgps2;
      input2.B->length = i;
    }
    else if (flag == 1) /* chisquare */
    {
      fscanf(fp2,"%lf\t%lf\t%lf\t%lf\t%lf",&t,&B.re, &B.im, &var.re, &var.im);  
      i=0;  
      while (!feof(fp2))
      {
        tgps2[i].gpsSeconds = (INT4)floor(t);
        tgps2[i].gpsNanoSeconds = (INT4)floor((fmod(t,1.0)*1.e9));       
        input2_chi.B->data[i].re = B.re;
        input2_chi.B->data[i].im = B.im; 
        input2_chi.var->data[i].re = var.re;
        input2_chi.var->data[i].im = var.im;   
        fscanf(fp2,"%lf\t%lf\t%lf\t%lf\t%lf",&t,&B.re, &B.im, &var.re, &var.im);    
        i++;  
      } 
      input2_chi.t = tgps2;
      input2_chi.B->length = i;
      input2_chi.var->length = i;
    }
    fclose(fp2);
  }
  
  fprintf(stderr, "I've read in the H1 data.\n");
 /* read data from H2 */
  if (iH2)
  {
    if (flag == 2) /* student-t */
    {
      fscanf(fp3,"%lf\t%lf\t%lf",&t,&B.re, &B.im);  
      i=0;  
      while (!feof(fp3))
      {
        tgps3[i].gpsSeconds = (INT4)floor(t);
        tgps3[i].gpsNanoSeconds = (INT4)floor((fmod(t,1.0)*1.e9));       
        input3.B->data[i].re = B.re;
        input3.B->data[i].im = B.im;    
        fscanf(fp3,"%lf\t%lf\t%lf",&t,&B.re, &B.im);    
        i++;  
      }  
    input3.t = tgps3;
    input3.B->length = i;
  }
    else if (flag == 1) /* chisquare */
    {
      fscanf(fp3,"%lf\t%lf\t%lf\t%lf\t%lf",&t,&B.re, &B.im, &var.re, &var.im);  
      i=0;  
      while (!feof(fp3))
      {
        tgps3[i].gpsSeconds = (INT4)floor(t);
        tgps3[i].gpsNanoSeconds = (INT4)floor((fmod(t,1.0)*1.e9));       
        input3_chi.B->data[i].re = B.re;
        input3_chi.B->data[i].im = B.im; 
        input3_chi.var->data[i].re = var.re;
        input3_chi.var->data[i].im = var.im;   
        fscanf(fp3,"%lf\t%lf\t%lf\t%lf\t%lf",&t,&B.re, &B.im, &var.re, &var.im);    
        i++;  
      } 
      input3_chi.t = tgps3;
      input3_chi.B->length = i;
      input3_chi.var->length = i;
    }
    fclose(fp3);
  }
  
  fprintf(stderr, "I've read in the H2 data.\n");
  /* read data from GEO */
  if (iGEO)
  {
    if (flag == 2) /* student-t */
    {
      fscanf(fp4,"%lf\t%lf\t%lf",&t,&B.re, &B.im);  
      i=0;  
      while (!feof(fp4))
      {
        tgps4[i].gpsSeconds = (INT4)floor(t);
        tgps4[i].gpsNanoSeconds = (INT4)floor((fmod(t,1.0)*1.e9));       
        input4.B->data[i].re = B.re;
        input4.B->data[i].im = B.im;    
        fscanf(fp4,"%lf\t%lf\t%lf",&t,&B.re, &B.im);    
        i++;  
      } 
      input4.t = tgps4;
      input4.B->length = i;
    }
    else if (flag == 1) /* chisquare */
    {
      fscanf(fp4,"%lf\t%lf\t%lf\t%lf\t%lf",&t,&B.re, &B.im, &var.re, &var.im);  
      i=0;  
      while (!feof(fp4))
      {
        tgps4[i].gpsSeconds = (INT4)floor(t);
        tgps4[i].gpsNanoSeconds = (INT4)floor((fmod(t,1.0)*1.e9));       
        input4_chi.B->data[i].re = B.re;
        input4_chi.B->data[i].im = B.im; 
        input4_chi.var->data[i].re = var.re;
        input4_chi.var->data[i].im = var.im;   
        fscanf(fp4,"%lf\t%lf\t%lf\t%lf\t%lf",&t,&B.re, &B.im, &var.re, &var.im);    
        i++;  
      } 
      input4_chi.t = tgps4;
      input4_chi.B->length = i;
      input4_chi.var->length = i;
    }
    fclose(fp4);
  }  
  /****************** END READ INPUT DATA (Bk's) ***************************/
  /******** read mesh input file *****************/
 
  /* construct name of mesh file and open file */
  sprintf(inmesh,"%s/inputs/%s", argv[2],argv[4]);
  fpmesh = fopen(inmesh,"r");
 
  fscanf(fpmesh,"%s\t%lf\t%lf\t%lf",&txt[0],&params1.meshH0[0], &params1.meshH0[1], &params1.meshH0[2]); 
  fprintf(stderr, "%s\t%e\t%e\t%f\t-> %e\n", txt, params1.meshH0[0], params1.meshH0[1],
  params1.meshH0[2],params1.meshH0[0]+params1.meshH0[1]*(params1.meshH0[2]-1.0));
  params2.meshH0[0] = params3.meshH0[0] = params4.meshH0[0] = params1.meshH0[0];
  params2.meshH0[1] = params3.meshH0[1] = params4.meshH0[1] = params1.meshH0[1];
  params2.meshH0[2] = params3.meshH0[2] = params4.meshH0[2] = params1.meshH0[2];

  fscanf(fpmesh,"%s\t%lf\t%lf\t%lf",&txt[0],&params1.meshCosIota[0], &params1.meshCosIota[1], &params1.meshCosIota[2]); 
  fprintf(stderr, "%s\t%e\t%e\t%f\t-> %e\n", txt, params1.meshCosIota[0], params1.meshCosIota[1],
  params1.meshCosIota[2], params1.meshCosIota[0]+ params1.meshCosIota[1]*(params1.meshCosIota[2]-1.0));
  params2.meshCosIota[0] = params3.meshCosIota[0] =  params4.meshCosIota[0] =  params1.meshCosIota[0];
  params2.meshCosIota[1] = params3.meshCosIota[1]  = params4.meshCosIota[1] =  params1.meshCosIota[1];
  params2.meshCosIota[2] = params3.meshCosIota[2] = params4.meshCosIota[2] =   params1.meshCosIota[2];

  fscanf(fpmesh,"%s\t%lf\t%lf\t%lf",&txt[0],&params1.meshPhase[0], &params1.meshPhase[1], &params1.meshPhase[2]);   
  fprintf(stderr,"%s\t%e\t%e\t%f\t-> %e\n", txt, params1.meshPhase[0], params1.meshPhase[1],
  params1.meshPhase[2], params1.meshPhase[0]+params1.meshPhase[1]*(params1.meshPhase[2]-1.0));
  params2.meshPhase[0] = params3.meshPhase[0] = params4.meshPhase[0] = params1.meshPhase[0]; 
  params2.meshPhase[1] = params3.meshPhase[1] = params4.meshPhase[1] = params1.meshPhase[1]; 
  params2.meshPhase[2] = params3.meshPhase[2] = params4.meshPhase[2] = params1.meshPhase[2]; 
 
  fscanf(fpmesh,"%s\t%lf\t%lf\t%lf",&txt[0],&params1.meshPsi[0], &params1.meshPsi[1], &params1.meshPsi[2]); 
  fprintf(stderr,"%s\t%e\t%e\t%f\t-> %e\n", txt, params1.meshPsi[0], params1.meshPsi[1],
  params1.meshPsi[2],params1.meshPsi[0]+params1.meshPsi[1]*(params1.meshPsi[2]-1.0)); 
  params2.meshPsi[0] = params3.meshPsi[0] =  params4.meshPsi[0] =  params1.meshPsi[0];
  params2.meshPsi[1] = params3.meshPsi[1] =  params4.meshPsi[1] =  params1.meshPsi[1];
  params2.meshPsi[2] = params3.meshPsi[2] =  params4.meshPsi[2] =  params1.meshPsi[2];
 
  fclose(fpmesh);

  /* allocate memory for 'chi square' matrix */
  
  if (iL1)
  { 
   output1.mChiSquare = NULL;
   LALDCreateVector( &status, &output1.mChiSquare,params1.meshH0[2]*params1.meshCosIota[2]*params1.meshPhase[2]*params1.meshPsi[2]);
  }
  
  if (iH1)
  {
    output2.mChiSquare = NULL;
    LALDCreateVector( &status, &output2.mChiSquare,params2.meshH0[2]*params2.meshCosIota[2]*params2.meshPhase[2]*params2.meshPsi[2]);
  }
  
  if (iH2)
  {
    output3.mChiSquare = NULL;
    LALDCreateVector( &status, &output3.mChiSquare,params3.meshH0[2]*params3.meshCosIota[2]*params3.meshPhase[2]*params3.meshPsi[2]);
  }
  if (iGEO)
  {    
    output4.mChiSquare = NULL;
    LALDCreateVector( &status, &output4.mChiSquare,params4.meshH0[2]*params4.meshCosIota[2]*params4.meshPhase[2]*params4.meshPsi[2]);
  }
  
  /* allocate memory for pdfs */
  if (iL1)
  {
   prob1.pdf = NULL; LALCreateVector(&status, &prob1.pdf, params1.meshH0[2]);
   prob1.cdf = NULL; LALCreateVector(&status, &prob1.cdf, params1.meshH0[2]);
   prob1.pdfPhase = NULL; LALCreateVector(&status, &prob1.pdfPhase, params1.meshPhase[2]);
   prob1.cdfPhase = NULL; LALCreateVector(&status, &prob1.cdfPhase, params1.meshPhase[2]);
   prob1.pdfPsi = NULL; LALCreateVector(&status, &prob1.pdfPsi, params1.meshPsi[2]);
   prob1.cdfPsi = NULL; LALCreateVector(&status, &prob1.cdfPsi, params1.meshPsi[2]);
   prob1.pdfCosIota = NULL; LALCreateVector(&status, &prob1.pdfCosIota, params1.meshCosIota[2]);
   prob1.cdfCosIota = NULL; LALCreateVector(&status, &prob1.cdfCosIota, params1.meshCosIota[2]);  
  }  
  if (iH1)
  {
   prob2.pdf = NULL; LALCreateVector(&status, &prob2.pdf, params2.meshH0[2]);
   prob2.cdf = NULL; LALCreateVector(&status, &prob2.cdf, params2.meshH0[2]);
   prob2.pdfPhase = NULL; LALCreateVector(&status, &prob2.pdfPhase, params2.meshPhase[2]);
   prob2.cdfPhase = NULL; LALCreateVector(&status, &prob2.cdfPhase, params2.meshPhase[2]);
   prob2.pdfPsi = NULL; LALCreateVector(&status, &prob2.pdfPsi, params2.meshPsi[2]);
   prob2.cdfPsi = NULL; LALCreateVector(&status, &prob2.cdfPsi, params2.meshPsi[2]);
   prob2.pdfCosIota = NULL; LALCreateVector(&status, &prob2.pdfCosIota, params2.meshCosIota[2]);
   prob2.cdfCosIota = NULL; LALCreateVector(&status, &prob2.cdfCosIota, params2.meshCosIota[2]);
  }
  if (iH2)
  { 
   prob3.pdf = NULL; LALCreateVector(&status, &prob3.pdf, params3.meshH0[2]);
   prob3.cdf = NULL; LALCreateVector(&status, &prob3.cdf, params3.meshH0[2]);
   prob3.pdfPhase = NULL; LALCreateVector(&status, &prob3.pdfPhase, params3.meshPhase[2]);
   prob3.cdfPhase = NULL; LALCreateVector(&status, &prob3.cdfPhase, params3.meshPhase[2]);
   prob3.pdfPsi = NULL; LALCreateVector(&status, &prob3.pdfPsi, params3.meshPsi[2]);
   prob3.cdfPsi = NULL; LALCreateVector(&status, &prob3.cdfPsi, params3.meshPsi[2]);
   prob3.pdfCosIota = NULL; LALCreateVector(&status, &prob3.pdfCosIota, params3.meshCosIota[2]);
   prob3.cdfCosIota = NULL; LALCreateVector(&status, &prob3.cdfCosIota, params3.meshCosIota[2]);
  } 
  if (iGEO)
  {
   prob4.pdf = NULL; LALCreateVector(&status, &prob4.pdf, params4.meshH0[2]);
   prob4.cdf = NULL; LALCreateVector(&status, &prob4.cdf, params4.meshH0[2]);
   prob4.pdfPhase = NULL; LALCreateVector(&status, &prob4.pdfPhase, params4.meshPhase[2]);
   prob4.cdfPhase = NULL; LALCreateVector(&status, &prob4.cdfPhase, params4.meshPhase[2]);
   prob4.pdfPsi = NULL; LALCreateVector(&status, &prob4.pdfPsi, params4.meshPsi[2]);
   prob4.cdfPsi = NULL; LALCreateVector(&status, &prob4.cdfPsi, params4.meshPsi[2]);
   prob4.pdfCosIota = NULL; LALCreateVector(&status, &prob4.pdfCosIota, params4.meshCosIota[2]);
   prob4.cdfCosIota = NULL; LALCreateVector(&status, &prob4.cdfCosIota, params4.meshCosIota[2]);
  }  
  
  /* allocate memory for pdfs for joint analysis if there is data from more than one detector */
  if (iH1 + iH2 + iL1 +iGEO > 1) 
  {
   prob.pdf = NULL; LALCreateVector(&status, &prob.pdf, params1.meshH0[2]);
   prob.cdf = NULL; LALCreateVector(&status, &prob.cdf, params1.meshH0[2]);  
   prob.pdfPhase = NULL; LALCreateVector(&status, &prob.pdfPhase, params1.meshPhase[2]);
   prob.cdfPhase = NULL; LALCreateVector(&status, &prob.cdfPhase, params1.meshPhase[2]);  
   prob.pdfPsi = NULL; LALCreateVector(&status, &prob.pdfPsi, params1.meshPsi[2]);
   prob.cdfPsi = NULL; LALCreateVector(&status, &prob.cdfPsi, params1.meshPsi[2]);  
   prob.pdfCosIota = NULL; LALCreateVector(&status, &prob.pdfCosIota, params1.meshCosIota[2]);
   prob.cdfCosIota = NULL; LALCreateVector(&status, &prob.cdfCosIota, params1.meshCosIota[2]);  
  }
        
  minChi = 0.0;
  
  /* calculate chisquare for each IFO and then marginalize over nuissance parameters */
  if (flag == 1)
  {
    if (iL1) {
      fprintf(stderr, "Entering LALCoarseFitToPulsar for L1\n");
      LALCoarseFitToPulsar(&status,&output1, &input1_chi, &params1);
      if(status.statusCode){
        fprintf(stderr,"Unexpectedly got error code %d and message %s\n",
	 status.statusCode, status.statusDescription);
         return 0;}    
      minChi += output1.chiSquare;
      fprintf(stderr, "Entering  LALPulsarMarginalize for L1\n");
      LALPulsarMarginalize(&status, &prob1, &output1, &params1);}
    if (iH1){
      fprintf(stderr, "Entering LALCoarseFitToPulsar for H1\n");
      LALCoarseFitToPulsar(&status,&output2, &input2_chi, &params2);
      if(status.statusCode){
        fprintf(stderr,"Unexpectedly got error code %d and message %s\n",
	 status.statusCode, status.statusDescription);
         return 0;}   
      minChi += output2.chiSquare;
      fprintf(stderr, "Entering  LALPulsarMarginalize for H1\n");
      LALPulsarMarginalize(&status, &prob2, &output2, &params2);}
   
    if (iH2){ 
      fprintf(stderr, "Entering LALCoarseFitToPulsar for H2\n");
      LALCoarseFitToPulsar(&status,&output3, &input3_chi, &params3);
      if(status.statusCode){
        fprintf(stderr,"Unexpectedly got error code %d and message %s\n",
	 status.statusCode, status.statusDescription);
         return 0;}         
      minChi += output3.chiSquare;
      fprintf(stderr, "Entering  LALPulsarMarginalize for H2\n");
      LALPulsarMarginalize(&status, &prob3, &output3, &params3);}
    
    if (iGEO){
      fprintf(stderr, "Entering LALCoarseFitToPulsar for GEO\n");
      LALCoarseFitToPulsar(&status,&output4, &input4_chi, &params4);
      if(status.statusCode){
        fprintf(stderr,"Unexpectedly got error code %d and message %s\n",
	 status.statusCode, status.statusDescription);
         return 0;}   
      minChi += output4.chiSquare;
      fprintf(stderr, "Entering  LALPulsarMarginalize for GEO\n");
      LALPulsarMarginalize(&status, &prob4, &output4, &params4); }  
  }
  else if (flag == 2)
  {
    if (iL1){
      fprintf(stderr, "Entering FitToPulsarStudentT for L1\n");
      LALFitToPulsarStudentT(&status,&output1, &input1, &params1);
      if(status.statusCode){
        fprintf(stderr,"Unexpectedly got error code %d and message %s\n",
	 status.statusCode, status.statusDescription);
         return 0;}   
      minChi += output1.chiSquare;
      fprintf(stderr, "Entering  LALPulsarMarginalize for L1\n");
      LALPulsarMarginalize(&status, &prob1, &output1, &params1);}
    if (iH1){
      fprintf(stderr, "Entering FitToPulsarStudentT for H1\n");
      LALFitToPulsarStudentT(&status,&output2, &input2, &params2);
      if(status.statusCode){
        fprintf(stderr,"Unexpectedly got error code %d and message %s\n",
	 status.statusCode, status.statusDescription);
         return 0;}   
      minChi += output2.chiSquare;
      fprintf(stderr, "Entering  LALPulsarMarginalize for H1\n");
      LALPulsarMarginalize(&status, &prob2, &output2, &params2);}
    if (iH2){
      fprintf(stderr, "Entering FitToPulsarStudentT for H2\n");
      LALFitToPulsarStudentT(&status,&output3, &input3, &params3);
      if(status.statusCode){
        fprintf(stderr,"Unexpectedly got error code %d and message %s\n",
	 status.statusCode, status.statusDescription);
         return 0;}   
      minChi += output3.chiSquare;
      fprintf(stderr, "Entering  LALPulsarMarginalize for H2\n");
      LALPulsarMarginalize(&status, &prob3, &output3, &params3);}
    if (iGEO){
      fprintf(stderr, "Entering FitToPulsarStudentT for GEO\n");
      LALFitToPulsarStudentT(&status,&output4, &input4, &params4);
      if(status.statusCode){
        fprintf(stderr,"Unexpectedly got error code %d and message %s\n",
	 status.statusCode, status.statusDescription);
         return 0;}   
      minChi += output4.chiSquare;
      fprintf(stderr, "Entering  LALPulsarMarginalize for GEO\n");
      LALPulsarMarginalize(&status, &prob4, &output4, &params4);}  
  }
  
  
/************ BEGIN MARGINALIZE FOR JOINT PDF *******************************/  
  /* if there is data from more that one IFO then calculate joint pdfs */
  if (iL1+iH1+iH2+iGEO>1) 
  {
    /* initialize to zero */
    for (iH0 = 0; iH0 < params1.meshH0[2]; iH0++) 
    {
      prob.pdf->data[iH0] = 0.0;
      prob.cdf->data[iH0] = 0.0;
    }
    for (iPhase=0; iPhase < params1.meshPhase[2]; iPhase++)
    {
      prob.pdfPhase->data[iPhase] = 0.0;
      prob.cdfPhase->data[iPhase] = 0.0;
    } 
    for (iPsi=0;iPsi< params1.meshPsi[2]; iPsi++)
    {
      prob.pdfPsi->data[iPsi] = 0.0;
      prob.cdfPsi->data[iPsi] = 0.0;
    }
    for (iCosIota = 0; iCosIota < params1.meshCosIota[2];iCosIota++)
    {
      prob.pdfCosIota->data[iCosIota] = 0.0;
      prob.cdfCosIota->data[iCosIota] = 0.0;
    }
       
     /* marginalize over angles to get p(h0|Bk) */
   for (iPsi = 0; iPsi < params1.meshPsi[2]; iPsi++)
     for (iPhase = 0; iPhase < params1.meshPhase[2]; iPhase++)
       for (iCosIota = 0; iCosIota < params1.meshCosIota[2];iCosIota++)
         for (iH0 = 0; iH0 < params1.meshH0[2]; iH0++){
           REAL8 ci=0., psi=0., prior=0., iota = 0.;
           
           arg = iH0 + params1.meshH0[2]*(iCosIota + 
params1.meshCosIota[2]*(iPhase + params1.meshPhase[2]*iPsi)); 
	   
           /* use Gaussian prior on psi and wrap around at edges */
         psi = params1.meshPsi[0] + (REAL8)iPsi*params1.meshPsi[1];
         if(psi < PSIMEAN - LAL_PI/4.){
          psi += LAL_PI/2.;
         }
         else if(PSIMEAN + LAL_PI/4. < psi){
          psi -= LAL_PI/2.;
         }
         
         prior = exp(-(psi - PSIMEAN)*(psi - PSIMEAN)/(2.*PSISIGMA*PSISIGMA));
         
         /* use Gaussian prior on iota */
         ci = params1.meshCosIota[0] + (REAL8)iCosIota*params1.meshCosIota[1];
         iota = acos(ci);
         
         /* wrap around at 0 and pi */
         if(iota < IOTAMEAN - LAL_PI/2.){
          iota += LAL_PI;
         }
         else if(IOTAMEAN + LAL_PI/2. < iota){
          iota -= LAL_PI;
         }
         
         prior *=
exp(-(iota-IOTAMEAN)*(iota-IOTAMEAN)/(2.*IOTASIGMA*IOTASIGMA));
           
	   if (iL1) outL1 = output1.mChiSquare->data[arg];
	   else outL1 = 0.0;
	   
	   if (iH1) outH1 = output2.mChiSquare->data[arg];
	   else outH1 = 0.0;
	   
	   if (iH2) outH2 = output3.mChiSquare->data[arg];
	   else outH2 = 0.0;
	   
	   if (iGEO) outGEO = output4.mChiSquare->data[arg];
	   else outGEO = 0.0;
	   
	   prob.pdf->data[iH0] +=  exp((minChi - (outL1 + outH1 + outH2 +
outGEO))/2.0) * prior;    
         }
 
   area = 0.0;
   for (iH0 = 0; iH0 < params1.meshH0[2]; iH0++)
     area += prob.pdf->data[iH0]*params1.meshH0[1];
   
   for (iH0 = 0; iH0 < params1.meshH0[2]; iH0++)
      prob.pdf->data[iH0] = prob.pdf->data[iH0]/area;
  
  prob.cdf->data[0] = prob.pdf->data[0]*params1.meshH0[1];
  for (iH0 = 1; iH0 < params1.meshH0[2]; iH0++)  
    prob.cdf->data[iH0] = prob.pdf->data[iH0]*params1.meshH0[1] + prob.cdf->data[iH0-1];
   
   /* marginalize over h0, psi, cosIota to get p(phase|Bk) */
   for (iPsi = 0; iPsi < params1.meshPsi[2]; iPsi++)
     for (iH0 = 0; iH0 < params1.meshH0[2]; iH0++)
       for (iCosIota = 0; iCosIota < params1.meshCosIota[2];iCosIota++)
         for (iPhase = 0; iPhase < params1.meshPhase[2]; iPhase++){
           REAL8 ci=0., psi=0., prior=0., iota = 0.;
           
           arg = iH0 + params1.meshH0[2]*(iCosIota + 
params1.meshCosIota[2]*(iPhase + params1.meshPhase[2]*iPsi)); 
	   
           /* use Gaussian prior on psi and wrap around at edges */
         psi = params1.meshPsi[0] + (REAL8)iPsi*params1.meshPsi[1];
         if(psi < PSIMEAN - LAL_PI/4.){
          psi += LAL_PI/2.;
         }
         else if(PSIMEAN + LAL_PI/4. < psi){
          psi -= LAL_PI/2.;
         }
         
         prior = exp(-(psi - PSIMEAN)*(psi - PSIMEAN)/(2.*PSISIGMA*PSISIGMA));
         
         /* use Gaussian prior on iota */
         ci = params1.meshCosIota[0] + (REAL8)iCosIota*params1.meshCosIota[1];
         iota = acos(ci);
         
         /* wrap around at 0 and pi */
         if(iota < IOTAMEAN - LAL_PI/2.){
          iota += LAL_PI;
         }
         else if(IOTAMEAN + LAL_PI/2. < iota){
          iota -= LAL_PI;
         }
         
         prior *=
exp(-(iota-IOTAMEAN)*(iota-IOTAMEAN)/(2.*IOTASIGMA*IOTASIGMA));

	   if (iL1) outL1 = output1.mChiSquare->data[arg];
	   else outL1 = 0.0;
	   
	   if (iH1) outH1 = output2.mChiSquare->data[arg];
	   else outH1 = 0.0;
	   
	   if (iH2) outH2 = output3.mChiSquare->data[arg];
	   else outH2 = 0.0;
	   
	   if (iGEO) outGEO = output4.mChiSquare->data[arg];
	   else outGEO = 0.0;
	   
	   prob.pdfPhase->data[iPhase] +=  exp((minChi - (outL1 + outH1 + outH2 +
outGEO))/2.0) * prior;    
         }
 
   area = 0.0;
   for (iPhase = 0; iPhase < params1.meshPhase[2]; iPhase++)
     area += prob.pdfPhase->data[iPhase]*params1.meshPhase[1];
   
   for (iPhase = 0; iPhase < params1.meshPhase[2]; iPhase++)
      prob.pdfPhase->data[iPhase] = prob.pdfPhase->data[iPhase]/area;
  
  prob.cdfPhase->data[0] = prob.pdfPhase->data[0]*params1.meshPhase[1];
  for (iPhase = 1; iPhase < params1.meshPhase[2]; iPhase++)  
    prob.cdfPhase->data[iPhase] = prob.pdf->data[iPhase]*params1.meshPhase[1] + prob.cdf->data[iPhase-1];
  
  /* marginalize over phase, h0, and cosIota to get p(Psi|Bk) */
   for (iH0 = 0; iH0 < params1.meshH0[2]; iH0++)
     for (iPhase = 0; iPhase < params1.meshPhase[2]; iPhase++)
       for (iCosIota = 0; iCosIota < params1.meshCosIota[2];iCosIota++)
         for (iPsi = 0; iPsi < params1.meshPsi[2]; iPsi++){
           REAL8 ci=0., psi=0., prior=0., iota = 0.;
           
           arg = iH0 + params1.meshH0[2]*(iCosIota + 
params1.meshCosIota[2]*(iPhase + params1.meshPhase[2]*iPsi)); 
	   
           /* use Gaussian prior on psi and wrap around at edges */
         psi = params1.meshPsi[0] + (REAL8)iPsi*params1.meshPsi[1];
         if(psi < PSIMEAN - LAL_PI/4.){
          psi += LAL_PI/2.;
         }
         else if(PSIMEAN + LAL_PI/4. < psi){
          psi -= LAL_PI/2.;
         }
         
         prior = exp(-(psi - PSIMEAN)*(psi - PSIMEAN)/(2.*PSISIGMA*PSISIGMA));
         
         /* use Gaussian prior on iota */
         ci = params1.meshCosIota[0] + (REAL8)iCosIota*params1.meshCosIota[1];
         iota = acos(ci);
         
         /* wrap around at 0 and pi */
         if(iota < IOTAMEAN - LAL_PI/2.){
          iota += LAL_PI;
         }
         else if(IOTAMEAN + LAL_PI/2. < iota){
          iota -= LAL_PI;
         }
         
         prior *=
exp(-(iota-IOTAMEAN)*(iota-IOTAMEAN)/(2.*IOTASIGMA*IOTASIGMA));

	   if (iL1) outL1 = output1.mChiSquare->data[arg];
	   else outL1 = 0.0;
	   
	   if (iH1) outH1 = output2.mChiSquare->data[arg];
	   else outH1 = 0.0;
	   
	   if (iH2) outH2 = output3.mChiSquare->data[arg];
	   else outH2 = 0.0;
	   
	   if (iGEO) outGEO = output4.mChiSquare->data[arg];
	   else outGEO = 0.0;
	   
	   prob.pdfPsi->data[iPsi] +=  exp((minChi - (outL1 + outH1 + outH2 +
outGEO))/2.0) * prior;    
         }
 
   area = 0.0;
   for (iPsi = 0; iPsi < params1.meshPsi[2]; iPsi++)
     area += prob.pdfPsi->data[iPsi]*params1.meshPsi[1];
   
   for (iPsi = 0; iPsi < params1.meshPsi[2]; iPsi++)
      prob.pdfPsi->data[iPsi] = prob.pdfPsi->data[iPsi]/area;
  
  prob.cdfPsi->data[0] = prob.pdfPsi->data[0]*params1.meshPsi[1];
  for (iPsi = 1; iPsi < params1.meshPsi[2]; iPsi++)  
    prob.cdfPsi->data[iPsi] = prob.pdfPsi->data[iPsi]*params1.meshPsi[1] + prob.cdfPsi->data[iPsi-1];
  
    /* marginalize over phase, psi, and h0 to get p(Cosiota|Bk) */
   for (iPsi = 0; iPsi < params1.meshPsi[2]; iPsi++)
     for (iPhase = 0; iPhase < params1.meshPhase[2]; iPhase++)
       for (iH0 = 0; iH0 < params1.meshH0[2];iH0++)
         for (iCosIota = 0; iCosIota < params1.meshCosIota[2]; iCosIota++){
           REAL8 ci=0., psi=0., prior=0., iota = 0.;
           
           arg = iH0 + params1.meshH0[2]*(iCosIota + 
params1.meshCosIota[2]*(iPhase + params1.meshPhase[2]*iPsi)); 
	   
           /* use Gaussian prior on psi and wrap around at edges */
         psi = params1.meshPsi[0] + (REAL8)iPsi*params1.meshPsi[1];
         if(psi < PSIMEAN - LAL_PI/4.){
          psi += LAL_PI/2.;
         }
         else if(PSIMEAN + LAL_PI/4. < psi){
          psi -= LAL_PI/2.;
         }
         
         prior = exp(-(psi - PSIMEAN)*(psi - PSIMEAN)/(2.*PSISIGMA*PSISIGMA));
         
         /* use Gaussian prior on iota */
         ci = params1.meshCosIota[0] + (REAL8)iCosIota*params1.meshCosIota[1];
         iota = acos(ci);
         
         /* wrap around at 0 and pi */
         if(iota < IOTAMEAN - LAL_PI/2.){
          iota += LAL_PI;
         }
         else if(IOTAMEAN + LAL_PI/2. < iota){
          iota -= LAL_PI;
         }
         
         prior *=
exp(-(iota-IOTAMEAN)*(iota-IOTAMEAN)/(2.*IOTASIGMA*IOTASIGMA));

	   if (iL1) outL1 = output1.mChiSquare->data[arg];
	   else outL1 = 0.0;
	   
	   if (iH1) outH1 = output2.mChiSquare->data[arg];
	   else outH1 = 0.0;
	   
	   if (iH2) outH2 = output3.mChiSquare->data[arg];
	   else outH2 = 0.0;
	   
	   if (iGEO) outGEO = output4.mChiSquare->data[arg];
	   else outGEO = 0.0;
	   
	   prob.pdfCosIota->data[iCosIota] +=  exp((minChi - (outL1 + outH1 + outH2 +
outGEO))/2.0) * prior;    
         }
 
   area = 0.0;
   for (iCosIota = 0; iCosIota < params1.meshCosIota[2]; iCosIota++)
     area += prob.pdfCosIota->data[iCosIota]*params1.meshCosIota[1];
   
   for (iCosIota = 0; iCosIota < params1.meshCosIota[2]; iCosIota++)
      prob.pdfCosIota->data[iCosIota] = prob.pdfCosIota->data[iCosIota]/area;
  
  prob.cdfCosIota->data[0] = prob.pdfCosIota->data[0]*params1.meshCosIota[1];
  for (iCosIota = 1; iCosIota < params1.meshCosIota[2]; iCosIota++)  
    prob.cdfCosIota->data[iCosIota] = prob.pdfCosIota->data[iCosIota]*params1.meshCosIota[1] + prob.cdfCosIota->data[iCosIota-1];
  
  } 
/************ END MARGINALIZE FOR JOINT PDF *******************************/  
       
/********************* BEGIN WRITING OUTPUT FILES *********************************/
  if (iL1)
  {
    if (flag == 1)
    {
      sprintf(outfile1,"%s/%s/pdf.%s_L1", argv[2], pulsar_name, pulsar_name); 
      sprintf(outfile1Phase,"%s/%s/pdfPhase.%s_L1", argv[2], pulsar_name, pulsar_name); 
      sprintf(outfile1Psi,"%s/%s/pdfPsi.%s_L1", argv[2], pulsar_name, pulsar_name); 
      sprintf(outfile1CosIota,"%s/%s/pdfCosIota.%s_L1", argv[2], pulsar_name, pulsar_name); 
    }
    else if (flag == 2)
    {
       sprintf(outfile1,"%s/pdfoutputs/pdf_st.%s_L1_with_priors", argv[2],
pulsar_name); 
       sprintf(outfile1Phase,"%s/pdfoutputs/pdfPhase_st.%s_L1_with_priors",
argv[2], pulsar_name); 
       sprintf(outfile1Psi,"%s/pdfoutputs/pdfPsi_st.%s_L1_with_priors", argv[2],
pulsar_name); 
       sprintf(outfile1CosIota,"%s/pdfoutputs/pdfCosIota_st.%s_L1_with_priors",
argv[2], pulsar_name); 
    }
   
    fp_pdf1 = fopen(outfile1, "w"); 
    fp_pdf1Phase = fopen(outfile1Phase, "w"); 
    fp_pdf1Psi = fopen(outfile1Psi, "w"); 
    fp_pdf1CosIota = fopen(outfile1CosIota, "w"); 
     
    for (i=0;i<params1.meshH0[2];i++)
      fprintf(fp_pdf1,"%e\t%e\t%e\n", params1.meshH0[0]+(float)i*params1.meshH0[1],prob1.pdf->data[i], prob1.cdf->data[i]);
    for (i=0;i<params1.meshPhase[2];i++)
      fprintf(fp_pdf1Phase,"%e\t%e\t%e\n", params1.meshPhase[0]+(float)i*params1.meshPhase[1],prob1.pdfPhase->data[i], prob1.cdfPhase->data[i]);    
    for (i=0;i<params1.meshPsi[2];i++)
      fprintf(fp_pdf1Psi,"%e\t%e\t%e\n", params1.meshPsi[0]+(float)i*params1.meshPsi[1],prob1.pdfPsi->data[i], prob1.cdfPsi->data[i]);
    for (i=0;i<params1.meshCosIota[2];i++)
      fprintf(fp_pdf1CosIota,"%e\t%e\t%e\n", params1.meshCosIota[0]+(float)i*params1.meshCosIota[1],prob1.pdfCosIota->data[i], prob1.cdfCosIota->data[i]);	   	   	       
  
    fclose(fp_pdf1); fclose(fp_pdf1Phase); fclose(fp_pdf1Psi); fclose(fp_pdf1CosIota);
  }
  
  if (iH1)
  {
     if (flag ==1)
     {  
       sprintf(outfile2,"%s/%s/pdf.%s_H1", argv[2], pulsar_name, pulsar_name);
       sprintf(outfile2Phase,"%s/%s/pdfPhase.%s_H1", argv[2], pulsar_name, pulsar_name);
       sprintf(outfile2Psi,"%s/%s/pdfPsi.%s_H1", argv[2], pulsar_name, pulsar_name);
       sprintf(outfile2CosIota,"%s/%s/pdfCosIota.%s_H1", argv[2], pulsar_name, pulsar_name);
     }
     else if (flag == 2)
     {
       sprintf(outfile2,"%s/pdfoutputs/pdf_st.%s_H1_with_priors", argv[2],
pulsar_name);
       sprintf(outfile2Phase,"%s/pdfoutputs/pdfPhase_st.%s_H1_with_priors",
argv[2], pulsar_name);
       sprintf(outfile2Psi,"%s/pdfoutputs/pdfPsi_st.%s_H1_with_priors", argv[2],
pulsar_name);
       sprintf(outfile2CosIota,"%s/pdfoutputs/pdfCosIota_st.%s_H1_with_priors",
argv[2], pulsar_name);
     }    
     fp_pdf2 = fopen(outfile2, "w"); 
     fp_pdf2Phase = fopen(outfile2Phase, "w"); 
     fp_pdf2Psi = fopen(outfile2Psi, "w"); 
     fp_pdf2CosIota = fopen(outfile2CosIota, "w"); 
  
    for (i=0;i<params1.meshH0[2];i++)
      fprintf(fp_pdf2,"%e\t%e\t%e\n", params2.meshH0[0]+(float)i*params2.meshH0[1], prob2.pdf->data[i], prob2.cdf->data[i]);
    for (i=0;i<params1.meshPhase[2];i++)
      fprintf(fp_pdf2Phase,"%e\t%e\t%e\n", params2.meshPhase[0]+(float)i*params2.meshPhase[1], prob2.pdfPhase->data[i], prob2.cdfPhase->data[i]);
    for (i=0;i<params1.meshPsi[2];i++)
      fprintf(fp_pdf2Psi,"%e\t%e\t%e\n", params2.meshPsi[0]+(float)i*params2.meshPsi[1], prob2.pdfPsi->data[i], prob2.cdfPsi->data[i]);
    for (i=0;i<params1.meshCosIota[2];i++)
      fprintf(fp_pdf2CosIota,"%e\t%e\t%e\n", params2.meshCosIota[0]+(float)i*params2.meshCosIota[1], prob2.pdfCosIota->data[i], prob2.cdfCosIota->data[i]);

    fclose(fp_pdf2); fclose(fp_pdf2Phase); fclose(fp_pdf2Psi); fclose(fp_pdf2CosIota); 
  }

  if (iH2)
  { 
     if (flag ==1)
     {  
       sprintf(outfile3,"%s/%s/pdf.%s_H2", argv[2], pulsar_name, pulsar_name);
       sprintf(outfile3Phase,"%s/%s/pdfPhase.%s_H2", argv[2], pulsar_name, pulsar_name);
       sprintf(outfile3Psi,"%s/%s/pdfPsi.%s_H2", argv[2], pulsar_name, pulsar_name);
       sprintf(outfile3CosIota,"%s/%s/pdfCosIota.%s_H2", argv[2], pulsar_name, pulsar_name);
     }
     else if (flag == 2)
     {
       sprintf(outfile3,"%s/pdfoutputs/pdf_st.%s_H2_with_priors", argv[2],
pulsar_name);
       sprintf(outfile3Phase,"%s/pdfoutputs/pdfPhase_st.%s_H2_with_priors",
argv[2], pulsar_name);
       sprintf(outfile3Psi,"%s/pdfoutputs/pdfPsi_st.%s_H2_with_priors", argv[2],
pulsar_name);
       sprintf(outfile3CosIota,"%s/pdfoutputs/pdfCosIota_st.%s_H2_with_priors",
argv[2], pulsar_name);
     }
     
    fp_pdf3 = fopen(outfile3, "w"); 
    fp_pdf3Phase = fopen(outfile3Phase, "w"); 
    fp_pdf3Psi = fopen(outfile3Psi, "w"); 
    fp_pdf3CosIota = fopen(outfile3CosIota, "w"); 
  
    for (i=0;i<params1.meshH0[2];i++)
      fprintf(fp_pdf3,"%e\t%e\t%e\n", params3.meshH0[0]+(float)i*params3.meshH0[1], prob3.pdf->data[i], prob3.cdf->data[i]);
    for (i=0;i<params1.meshPhase[2];i++)
      fprintf(fp_pdf3Phase,"%e\t%e\t%e\n", params3.meshPhase[0]+(float)i*params3.meshPhase[1], prob3.pdfPhase->data[i], prob3.cdfPhase->data[i]);	    
    for (i=0;i<params1.meshPsi[2];i++)
      fprintf(fp_pdf3Psi,"%e\t%e\t%e\n", params3.meshPsi[0]+(float)i*params3.meshPsi[1], prob3.pdfPsi->data[i], prob3.cdfPsi->data[i]);  
    for (i=0;i<params1.meshCosIota[2];i++)
      fprintf(fp_pdf3CosIota,"%e\t%e\t%e\n", params3.meshCosIota[0]+(float)i*params3.meshCosIota[1], prob3.pdfCosIota->data[i], prob3.cdfCosIota->data[i]);  
   
    fclose(fp_pdf3); fclose(fp_pdf3Phase); fclose(fp_pdf3Psi); fclose(fp_pdf3CosIota); 
  }
       
  if (iGEO){ 
       if (flag ==1)
     {  
       sprintf(outfile4,"%s/%s/pdf.%s_GEO", argv[2], pulsar_name, pulsar_name);
       sprintf(outfile4Phase,"%s/%s/pdfPhase.%s_GEO", argv[2], pulsar_name, pulsar_name);
       sprintf(outfile4Psi,"%s/%s/pdfPsi.%s_GEO", argv[2], pulsar_name, pulsar_name);
       sprintf(outfile4CosIota,"%s/%s/pdfCosIota.%s_GEO", argv[2], pulsar_name, pulsar_name);
     }
     else if (flag == 2)
     {
       sprintf(outfile4,"%s/pdfoutputs/pdf_st.%s_GEO_with_priors", argv[2],
pulsar_name);
       sprintf(outfile4Phase,"%s/pdfoutputs/pdfPhase_st.%s_GEO_with_priors",
argv[2], pulsar_name);
       sprintf(outfile4Psi,"%s/pdfoutputs/pdfPsi_st.%s_GEO_with_priors",
argv[2], pulsar_name);
       sprintf(outfile4CosIota,"%s/pdfoutputs/pdfCosIota_st.%s_GEO_with_priors",
argv[2], pulsar_name);
     }
     
    fp_pdf4 = fopen(outfile4, "w"); 
    fp_pdf4Phase = fopen(outfile4Phase, "w"); 
    fp_pdf4Psi = fopen(outfile4Psi, "w"); 
    fp_pdf4CosIota = fopen(outfile4CosIota, "w"); 
  
    for (i=0;i<params1.meshH0[2];i++)
      fprintf(fp_pdf4,"%e\t%e\t%e\n", params4.meshH0[0]+(float)i*params4.meshH0[1], prob4.pdf->data[i], prob4.cdf->data[i]);
    for (i=0;i<params1.meshPhase[2];i++)
      fprintf(fp_pdf4Phase,"%e\t%e\t%e\n", params4.meshPhase[0]+(float)i*params4.meshPhase[1], prob4.pdfPhase->data[i], prob4.cdfPhase->data[i]); 
    for (i=0;i<params1.meshPsi[2];i++)
      fprintf(fp_pdf4Psi,"%e\t%e\t%e\n", params4.meshPsi[0]+(float)i*params4.meshPsi[1], prob4.pdfPsi->data[i], prob4.cdfPsi->data[i]); 
    for (i=0;i<params1.meshCosIota[2];i++)
      fprintf(fp_pdf4CosIota,"%e\t%e\t%e\n", params4.meshCosIota[0]+(float)i*params4.meshCosIota[1], prob4.pdfCosIota->data[i], prob4.cdfCosIota->data[i]);
  
    fclose(fp_pdf4); fclose(fp_pdf4Phase); fclose(fp_pdf4Psi); fclose(fp_pdf4CosIota);  
  }
  
  if (iL1+iH1+iH2+iGEO>1)
  {
       if (flag ==1)
     {  
       sprintf(outfile,"%s/%s/pdf.%s_Joint", argv[2], pulsar_name, pulsar_name);
       sprintf(outfilePhase,"%s/%s/pdfPhase.%s_Joint", argv[2], pulsar_name, pulsar_name);
       sprintf(outfilePsi,"%s/%s/pdfPsi.%s_Joint", argv[2], pulsar_name, pulsar_name);
       sprintf(outfileCosIota,"%s/%s/pdfCosIota.%s_Joint", argv[2], pulsar_name, pulsar_name);
     }
     else if (flag == 2)
     {
       sprintf(outfile,"%s/pdfoutputs/pdf_st.%s_Joint_with_priors", argv[2],
pulsar_name);
       sprintf(outfilePhase,"%s/pdfoutputs/pdfPhase_st.%s_Joint_with_priors",
argv[2], pulsar_name);
       sprintf(outfilePsi,"%s/pdfoutputs/pdfPsi_st.%s_Joint_with_priors",
argv[2], pulsar_name);
      
sprintf(outfileCosIota,"%s/pdfoutputs/pdfCosIota_st.%s_Joint_with_priors",
argv[2], pulsar_name);
     }
     
    fp_joint = fopen(outfile, "w"); 
    fp_jointPhase = fopen(outfilePhase, "w"); 
    fp_jointPsi = fopen(outfilePsi, "w"); 
    fp_jointCosIota = fopen(outfileCosIota, "w"); 
  
    for (i=0;i<params1.meshH0[2];i++)
      fprintf(fp_joint,"%e\t%e\t%e\n", params1.meshH0[0]+(float)i*params1.meshH0[1], prob.pdf->data[i], prob.cdf->data[i]);        
    for (i=0;i<params1.meshPhase[2];i++)
      fprintf(fp_jointPhase,"%e\t%e\t%e\n", params1.meshPhase[0]+(float)i*params1.meshPhase[1], prob.pdfPhase->data[i], prob.cdfPhase->data[i]);      
    for (i=0;i<params1.meshPsi[2];i++)
      fprintf(fp_jointPsi,"%e\t%e\t%e\n", params1.meshPsi[0]+(float)i*params1.meshPsi[1], prob.pdfPsi->data[i], prob.cdfPsi->data[i]);      
    for (i=0;i<params1.meshCosIota[2];i++)
      fprintf(fp_jointCosIota,"%e\t%e\t%e\n", params1.meshCosIota[0]+(float)i*params1.meshCosIota[1], prob.pdfCosIota->data[i], prob.cdfCosIota->data[i]);      
  
    fclose(fp_joint);fclose(fp_jointPhase);fclose(fp_jointPsi);fclose(fp_jointCosIota);
  } 
 /********************* END WRITING OUTPUT FILES*********************************/

/* print out best fit information for each IFO */  
  if (iL1){fprintf(stderr, "BEST FIT FOR L1:\n");
 fprintf(stderr,"h0 = %e\tcosIota = %f\tpsi = %f\tphase = %f\nchisquare = %f\n",
  output1.h0, output1.cosIota, output1.psi, output1.phase, output1.chiSquare);}

  if (iH1){fprintf(stderr, "\nBEST FIT FOR H1:\n");
 fprintf(stderr,"h0 = %e\tcosIota = %f\tpsi = %f\tphase = %f\nchisquare = %f\n",
  output2.h0, output2.cosIota, output2.psi, output2.phase, output2.chiSquare);}

  if (iH2){fprintf(stderr, "\nBEST FIT FOR H2:\n");
 fprintf(stderr,"h0 = %e\tcosIota = %f\tpsi = %f\tphase = %f\nchisquare = %f\n",
  output3.h0, output3.cosIota, output3.psi, output3.phase, output3.chiSquare);} 
  
  if (iGEO){fprintf(stderr, "\nBEST FIT FOR GEO:\n");
 fprintf(stderr,"h0 = %e\tcosIota = %f\tpsi = %f\tphase = %f\nchisquare = %f\n",
  output4.h0, output4.cosIota, output4.psi, output4.phase, output4.chiSquare);}      
 fprintf(stderr,"iL1 = %d\tiH1=%d\tiH2=%d\tiGEO=%d\n", iL1, iH1, iH2, iGEO);
  
/* free allocated memory */  
  if (iL1) {
    if (flag == 1) LALZDestroyVector(&status, &input1_chi.B);
    else if (flag ==2)  LALZDestroyVector(&status, &input1.B);
    
    LALDDestroyVector(&status, &output1.mChiSquare); 
    LALDestroyVector(&status, &prob1.pdf); LALDestroyVector(&status, &prob1.cdf);
    LALDestroyVector(&status, &prob1.pdfPhase); LALDestroyVector(&status, &prob1.cdfPhase);
    LALDestroyVector(&status, &prob1.pdfPsi); LALDestroyVector(&status, &prob1.cdfPsi);
    LALDestroyVector(&status, &prob1.pdfCosIota); LALDestroyVector(&status, &prob1.cdfCosIota);  
  }

  if (iH1){
    if (flag == 1) LALZDestroyVector(&status, &input2_chi.B);
    else if (flag == 2) LALZDestroyVector(&status, &input2.B); 
    LALDDestroyVector(&status, &output2.mChiSquare); 
    LALDestroyVector(&status, &prob2.pdf); LALDestroyVector(&status, &prob2.cdf);
    LALDestroyVector(&status, &prob2.pdfPhase); LALDestroyVector(&status, &prob2.cdfPhase);
    LALDestroyVector(&status, &prob2.pdfPsi); LALDestroyVector(&status, &prob2.cdfPsi);
    LALDestroyVector(&status, &prob2.pdfCosIota); LALDestroyVector(&status, &prob2.cdfCosIota);   
  }

  if (iH2){
    if (flag == 1) LALZDestroyVector(&status, &input3_chi.B); 
    else if (flag == 2) LALZDestroyVector(&status, &input3.B); 
    LALDDestroyVector(&status, &output3.mChiSquare); 
    LALDestroyVector(&status, &prob3.pdf); LALDestroyVector(&status, &prob3.cdf);
    LALDestroyVector(&status, &prob3.pdfPhase); LALDestroyVector(&status, &prob3.cdfPhase);
    LALDestroyVector(&status, &prob3.pdfPsi); LALDestroyVector(&status, &prob3.cdfPsi);
    LALDestroyVector(&status, &prob3.pdfCosIota); LALDestroyVector(&status, &prob3.cdfCosIota); 
  }
 
  if (iGEO){
    if (flag == 1) LALZDestroyVector(&status, &input4_chi.B); 
    else if (flag == 2) LALZDestroyVector(&status, &input4.B); 
    LALDDestroyVector(&status, &output4.mChiSquare); 
    LALDestroyVector(&status, &prob4.pdf); LALDestroyVector(&status, &prob4.cdf);
    LALDestroyVector(&status, &prob4.pdfPhase); LALDestroyVector(&status, &prob4.cdfPhase);
    LALDestroyVector(&status, &prob4.pdfPsi); LALDestroyVector(&status, &prob4.cdfPsi);
    LALDestroyVector(&status, &prob4.pdfCosIota); LALDestroyVector(&status, &prob4.cdfCosIota);    
  }
  
  if (iL1+iH1+iH2+iGEO>1)
  {  
    LALDestroyVector(&status, &prob.pdf); LALDestroyVector(&status, &prob.cdf);
    LALDestroyVector(&status, &prob.pdfPhase); LALDestroyVector(&status, &prob.cdfPhase);
    LALDestroyVector(&status, &prob.pdfPsi); LALDestroyVector(&status, &prob.cdfPsi);
    LALDestroyVector(&status, &prob.pdfCosIota); LALDestroyVector(&status, &prob.cdfCosIota); 
  }
    LALCheckMemoryLeaks();
    return(0);
}  /* main */
 
