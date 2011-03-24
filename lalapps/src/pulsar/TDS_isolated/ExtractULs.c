/*
*  Copyright (C) 2007  Rejean Dupuis
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
 * This code extracts the 95% upper limits (from the results of ComputePFD)
 * and then calculates the corresponding limits on ellipticity given the
 * pulsar frequency and distance.  The output data from ComputePDF are
 * linearly interpolated to get h95. This limit is compared to the spindown
 * based upper limits where applicable.
 */

/*****************************************************************************/
/* This code extracts the 95% upper limits (from the results of ComputePFD)  */     
/* and then calculates the corresponding limits on ellipticity given the     */
/* pulsar frequency and distance.  The output data from ComputePDF are       */
/* linearly interpolated to get h95. This limit is compared to the spindown  */
/* based upper limits where applicable.                                      */
/*                                                                           */
/*			               Réjean Dupuis                         */
/*                                                                           */
/*                       University of Glasgow - last modified 21/04/2004    */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define UL 0.95
#define MX 28

int main(int argc, char **argv)
{
   FILE *fpin=NULL, *fp=NULL, *fpout=NULL;
   char pulsar_name[MX][16], infile[256], psrinput[256];
   int i, num, flag;
   double f0[MX], P0[MX], P1[MX], DIST[MX]; 
   double h=0, pd=0, cp=0, h2,pd2,cp2;
   char outfile[128];
   double H1UL, H2UL, L1UL, JointUL, eccUL, eccSD;
   
  
  if (argv[1]==NULL||argv[2]==NULL)
  {
    printf("1. 1 for chi square, 2 for st\n2. outfile name\n");
    return(1);
  }
  
  flag = atoi(argv[1]);
  sprintf(outfile, "%s", argv[2]);
  printf("writting to %s\n", outfile);  

  /* open file with information on pulsar distance, and spindown, and frequency */
  
  fpin = fopen("pulsar.dist", "r");
    
  i = 0;
  while (5==fscanf(fpin,"%s %lf %lf %lf %lf", &pulsar_name[i][0], &f0[i], &P0[i], &P1[i], &DIST[i]))
    i++; 
    
  fclose(fpin);
 
 /* for each ifo (and joint) read pdfs, get h95, then calculate e95, then write */
  num = i;

  fpout = fopen(outfile, "w");
  fprintf(fpout, "#name freq dist(kpc) P0 P1 ULH1 ULH2 ULL1 JointUL EpsilonUL ratio\n"); 
  for (i=0; i<num;i++)
  {     
    if (flag==1) sprintf(psrinput,"%s/pdf.%s_H1",pulsar_name[i],pulsar_name[i]);
    else if (flag==2) sprintf(psrinput,"%s/pdf_st.%s_H1",pulsar_name[i],pulsar_name[i]);
    else return(1);
    
    fp=fopen(psrinput,"r");
    cp=0.0;
    while (cp<0.95)
    { 
      h2 = h;
      pd2 = pd;
      cp2 = cp;
      fscanf(fp, "%lf %lf %lf", &h, &pd, &cp);    
    }
  
    fclose(fp);
    H1UL = h2 + (UL-cp2)*(h-h2)/(cp-cp2);
      
    cp = 0.0;
    if (flag==1) sprintf(psrinput,"%s/pdf.%s_H2",pulsar_name[i],pulsar_name[i]);
    else if (flag==2) sprintf(psrinput,"%s/pdf_st.%s_H2",pulsar_name[i],pulsar_name[i]);
    else return(1);
    
    fp=fopen(psrinput,"r");
    while (cp<0.95)
    { 
      h2 = h;
      pd2 = pd;
      cp2 = cp;
      fscanf(fp, "%lf %lf %lf", &h, &pd, &cp);    
    }
  
    fclose(fp);
    H2UL = h2 + (UL-cp2)*(h-h2)/(cp-cp2);

    cp = 0.0;
    if (flag==1) sprintf(psrinput,"%s/pdf.%s_L1",pulsar_name[i],pulsar_name[i]);
    else if (flag==2) sprintf(psrinput,"%s/pdf_st.%s_L1",pulsar_name[i],pulsar_name[i]);
    else return(1);
    
    fp=fopen(psrinput,"r");
    while (cp<0.95)
    { 
      h2 = h;
      pd2 = pd;
      cp2 = cp;
      fscanf(fp, "%lf %lf %lf", &h, &pd, &cp);    
    }
  
    fclose(fp);
    L1UL = h2 + (UL-cp2)*(h-h2)/(cp-cp2);

    cp = 0.0;
    if (flag==1) sprintf(psrinput,"%s/pdf.%s_Joint",pulsar_name[i],pulsar_name[i]);
    else if (flag==2) sprintf(psrinput,"%s/pdf_st.%s_Joint",pulsar_name[i],pulsar_name[i]);
    else return(1);
    
    fp=fopen(psrinput,"r");
    while (cp<0.95)
    { 
      h2 = h;
      pd2 = pd;
      cp2 = cp;
      fscanf(fp, "%lf %lf %lf", &h, &pd, &cp);    
    }
  
    fclose(fp);
    JointUL = h2 + (UL-cp2)*(h-h2)/(cp-cp2);
    
   /* if (f0[i]>0) eccUL = JointUL*30352.82952*DIST[i]*1000.0*3.09e16/(f0[i]*f0[i]);   */
    if (f0[i]>0) eccUL = JointUL*30662.10384*DIST[i]*3.08568025e19/(f0[i]*f0[i]);   
    else
    {
      printf("error1\n");
      return(1);
    }
    /* if spindown is positive then calculate ratio */
    
    if (P1[i] > 0)
    {
     eccSD = 0.0057*P0[i]*sqrt(P0[i])*sqrt(P1[i]*1.e15);
     fprintf(fpout, "%s\t%f\t%f\t%e\t%f\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
       pulsar_name[i], f0[i], P0[i], P1[i], DIST[i],H1UL, H2UL, L1UL, JointUL, eccUL, eccUL/eccSD); 
    }
    else {
         fprintf(fpout, "%s\t%f\t%f\t%e\t%f\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%f\n",
       pulsar_name[i], f0[i], P0[i], P1[i], DIST[i],H1UL, H2UL, L1UL, JointUL, eccUL, -1.0);    
    }   
  }
  fclose(fpout);


}
