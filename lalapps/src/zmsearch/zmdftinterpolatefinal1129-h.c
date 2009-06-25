/*
*  Copyright (C) 2007 Jolien Creighton, Saikat Ray-Majumder
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

#include<stdio.h>
#include<math.h>
#include<stdarg.h> 
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/Calibration.h>
#include <lal/AVFactories.h>
#include <lal/LALDatatypes.h>
#include <lal/LALConstants.h>
#include <lal/RealFFT.h>
#include <lal/Interpolate.h>

NRCSID( PLAYC, "$Id$" );

#define PLAYC_ENORM  0
#define PLAYC_ESUB   1
#define PLAYC_EARG   2
#define PLAYC_EVAL   3
#define PLAYC_EFILE  4
#define PLAYC_EINPUT 5
#define PLAYC_EMEM   6

#define PLAYC_MSGENORM  "Normal exit"
#define PLAYC_MSGESUB   "Subroutine failed"
#define PLAYC_MSGEARG   "Error parsing arguments"
#define PLAYC_MSGEVAL   "Input argument out of valid range"
#define PLAYC_MSGEFILE  "Could not open file"
#define PLAYC_MSGEINPUT "Error reading file"
#define PLAYC_MSGEMEM   "Out of memory"

/* Usage format string. */       
#define USAGE "Usage: %s [-o outfile] [-i infile]\n"

/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                           do                                                                             if (lalDebugLevel & LALERROR )                                                 {                                                                                LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"                              "        %s %s\n", (code), *argv, __FILE__,                                      __LINE__, PLAYC, statement ? statement :                                      "", (msg) );                                                 }                                                                              while(0)

INT4 lalDebugLevel = LALMSGLVL3;


void swap4(int n,char *here);
int readZM(const char* fname,  float **time,  float **amplitude); 

int main(int argc, char **argv)
{
  static LALStatus stat;                                                      
  INT4 arg=1;                   /* counters                            */
  
  CHAR *outfile = NULL;           /* name of ascii outfile */
  CHAR *infile = NULL;          /* name of ilwd outfile */
  FILE *fp;
  int i,j,N[80],p,l, k, m;
  float n[80];
  float dt=100/4096.0;            /*time diff in ms.=0.02441ms*/
  float *time, *amp;
  char filename[10], outputfile[30];   
  float t[80][4000],a[80][4000], h[80][4000];  
  float x; 
  float h_intpolat[80][20000], t_intpolat[80][20000];

  float      t_target[5];
  float      h_target[5];
  float      target;
  SInterpolatePar intpar   ={5,t_target,h_target};
  SInterpolateOut intout;





  /*******************************************************************
   * PARSE ARGUMENTS (arg stores the current position)               *
   *******************************************************************/

  if (argc <= 1){
    LALPrintError( USAGE, *argv );
    return 0;
  }

  while ( arg < argc ) {
    /* Parse output file option. */
    if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        outfile = argv[arg++];
      }else{
        ERROR( PLAYC_EARG, PLAYC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return PLAYC_EARG;
          }
    }
    else if ( !strcmp( argv[arg], "-i" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        infile = argv[arg++];
      }else{
        ERROR( PLAYC_EARG, PLAYC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return PLAYC_EARG;
      }                      
    }
    else if ( !strcmp( argv[arg], "-i" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        infile = argv[arg++];
      }else{
        ERROR( PLAYC_EARG, PLAYC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return PLAYC_EARG;
      }
    }
    /* Check for unrecognized options. */
    else if ( argv[arg][0] == '-' ) {
      LALPrintError( USAGE, *argv );
      return PLAYC_EARG;
    }
  } /* End of argument parsing loop. */

  /**********************************************************************************************/
  N[0]=0;
  /* read in Zwerger Muller */
  for(i=0;i<78;i++){

  snprintf(filename,10,"zm-%d.bin",i);

  N[i] = readZM(filename,&time,&amp);
  
  for(j=0;j<N[i];j++){   
  t[i][j] = 1000*(*(time+j));
  a[i][j] = *(amp+j);
     }
 free(time);
 free(amp);
 }
    
  {
      fp = fopen(outfile,"w");
      for(j=0;j<N[60];j++){
      fprintf(fp,"%f %f\n",t[60][j], a[60][j]);
    }  
      fclose(fp); 
  }
  
  /*calculate h+ from A*/
  for(i=0;i<78;i++){
  for(j=0;j<N[i];j++){
    h[i][j] = (0.273)*a[i][j];             /*look into pg.2 of documentation:here angle=90 & R=1cm.*/
  }
  }
  

      /************************ find l *****************************************/

  fp = fopen("no.dat","w");
 for(i=0;i<78;i++){
  n[i] = (t[i][N[i]-1]-t[i][0])/dt;
  fprintf(fp,"%f\n",n[i]);
}
  fclose(fp);
 
 x=0;
   if(n[1]>=n[0])
     x = n[1];
   else if(n[1]<=n[0])
     x=n[0];
   for(i=2;i<78;i++){
     if(n[i]>=x)
       x = n[i];
     else if(n[i]<=x)
       x = x;
   }
 printf("%f\n",x);

 for(p=0;pow(2,p)<x;++p){
 }

 l=pow(2,p);                          /* l is the required arraysize */

 printf("%d\n",l);

 /******************************interpolate******************************************/
 t_intpolat[0][0]=0;
 h_intpolat[0][0]=0; 

 for(i=0;i<78;i++){
   t_intpolat[i][0]=0;
   h_intpolat[i][0]=0;
 }

 for(i=0;i<78;i++){
   t_intpolat[i][0] = t[i][0];  
   for(j=0;j<(l-1);j++){
     t_intpolat[i][j+1] = t_intpolat[i][j]+dt;
     h_intpolat[i][j+1] = 0;
   }
 }
 
 for(i=0;i<78;i++){ 
   for(j=0;t_intpolat[i][j]<=t[i][N[i]-1];j++){
     for(k=0;t[i][k]<=t_intpolat[i][j];k++){

     }
  
    for(m=0;m<4;m++){
     t_target[m]=0;
     h_target[m]=0;
     }

     if(k>=2 && k<=(N[i]-3)){
       for(m=0;m<=4;m++){
	 t_target[m]=t[i][k-2+m];
	 h_target[m]=h[i][k-2+m];
       }
     }
     else if(k==1){
       for(m=0;m<=4;m++){
	 t_target[m]=t[i][k-1+m];
	 h_target[m]=h[i][k-1+m];
       }
     }
     else if(k==0){
       for(m=0;m<=4;m++){
	 t_target[m]=t[i][k+m];
	 h_target[m]=h[i][k+m];
       }
     }
     else if(k==(N[i]-2)){
       for(m=0;m<=4;m++){
	 t_target[m]=t[i][k-3+m];
	 h_target[m]=h[i][k-3+m];
       }
     }
     else if(k==(N[i]-1)){
       for(m=0;m<=4;m++){
	 t_target[m]=t[i][k-4+m];
	 h_target[m]=h[i][k-4+m];
       }
     }

     target = t_intpolat[i][j];
     
     LALSPolynomialInterpolation( &stat, &intout, target, &intpar );
     h_intpolat[i][j]=intout.y;
   }
 }

 for(i=0;i<78;i++){
 snprintf(outputfile,30,"zminpol-%d.dat",i);

 fp = fopen(outputfile,"w");
 
 for(j=0;j<l;j++){
   fprintf(fp,"%f %f\n",t_intpolat[i][j],h_intpolat[i][j]);
 }

 fclose(fp);
 }
 
 return 0;
}


  /* function to read in ZM waveforms */
  int readZM(const char* fname,  float **time,  float **amplitude)
  {
    long n;
    int i;
    FILE *fp;     


    fp = fopen(fname,"r+b");
   
    fread(&i, sizeof(long), 1, fp);
    fread(&n, sizeof(long), 1, fp);
    fread(&i, sizeof(long), 1, fp);
    fread(&i, sizeof(long), 1, fp);
    swap4(1,(char *)&n);
   
    
    (*time) = (float *)malloc( n*sizeof(float));
    fread(*time, sizeof(float), n, fp);
    swap4(n,(char *)(*time));
   
    (*amplitude) = (float *)malloc( n*sizeof(float));
    fread(*amplitude, sizeof(float), n, fp);
    swap4(n,(char *)(*amplitude));     


    fclose(fp);
    
    return n;
} 


    /* byte swap because the data was written on a machine which is big-endian */
     void swap4(int n, char *b1) {
     char *b2,*b3,*b4;
     char temp;
     int i; for (i=0;i<n;i++) {
     b2=b1+1;
     b3=b2+1;
     b4=b3+1; 

     temp=*b1;
     *b1=*b4;
     *b4=temp;
     temp=*b2;
     *b2=*b3;
     *b3=temp;
     b1+=4;
     }
return;
}








