/*
*  Copyright (C) 2007 Xavier Siemens
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

/*********************************************************************************/
/*            Cosmic string burst rate computation code for small loops          */
/*                                                                               */
/*                  Xavier Siemens, Jolien Creighton, Irit Maor                  */
/*                                                                               */
/*                         UWM/Caltech - September 2006                          */
/*********************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <stdarg.h>
#include "gsl/gsl_interp.h"
#include <gsl/gsl_errno.h>
#include <lal/cs_cosmo.h>
#include <lal/cs_lambda_cosmo.h>
#include <lal/LALStdio.h>

#define CUSPS_PER_LOOP 1.0		/* c */
#define LOOP_RAD_POWER 50.0		/* Gamma */

struct CommandLineArgsTag {
  double f;                /* frequency */
  double logGmustart;      
  double logGmuend;        
  int  nGmu;
  double logepsilonstart;
  double logepsilonend;    
  int  nepsilon;
  double n;
  double logpstart;        
  double logpend;          
  int  np;
  char *efficiencyfile;
} CLA;

double *eff, *lnamp, *amp, *Deff, Dlnz;  /* arrays for efficiency, amplitude and efficiency errors */
double *zofA, *dzdA, *dRdz;
double H0 = LAMBDA_H_0;
int Namp;   /* size of amplitude/eff array */


#define PROGRAM_NAME "cs_gamma"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"

/***************************************************************************/
/* Reads the command line */
int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);
int ReadEfficiencyFile(struct CommandLineArgsTag CLA);
int findzofA(double Gmu, double alpha);
int finddRdz(double Gmu, double alpha, double f, double Gamma);

/************************************* MAIN PROGRAM *************************************/

int main( int argc, char *argv[] )
{
	/* parameters of strings */
	double Gamma   = LOOP_RAD_POWER;	/* power radiated by loops */
	double c       = CUSPS_PER_LOOP;	/* cusps per loop */

	double p,epsilon,f,Gmu,alpha;
	double gammaAverage = 0, gammaMin = 0, gammaMax =0;
	char filename[FILENAME_MAX];
	FILE *fp;
	int i,j,k,l;

	/* read the command line */
	if (ReadCommandLine(argc,argv,&CLA)) return 1;
	f=CLA.f;

	/* read efficiency function */
	if (ReadEfficiencyFile(CLA)) return 2;
	
	snprintf( filename, sizeof( filename ), "gamma.dat");
	fp = fopen( filename, "w" );
	fprintf( fp,"%%     p           n           epsilon         Gmu       gammaAverage    gammaMin      gammaMax\n");
	for ( i = 0; i <  CLA.nepsilon; i++ )
	  {
	    epsilon=pow(10.0,CLA.logepsilonstart+i*(CLA.logepsilonend-CLA.logepsilonstart)/(CLA.nepsilon-1));
	    
	    for ( j = 0; j <  CLA.nGmu; j++ )
	      {
		Gmu=pow(10.0,CLA.logGmustart+j*(CLA.logGmuend-CLA.logGmustart)/(CLA.nGmu-1));
		alpha = epsilon * pow( Gamma * Gmu, CLA.n );
			       
		/* find the z's corresponding to those A's */
		if(findzofA(Gmu, alpha)) return 3;
		
		/* compute the rate derivative at those z's */
		if(finddRdz(Gmu, alpha, f, Gamma)) return 4;
			
		for ( k = 0; k <  CLA.np; k++ )
		  {
		    p=pow(10.0, CLA.logpstart+k*(CLA.logpend-CLA.logpstart)/(CLA.np));
		    
 		    fprintf(stdout,"%%Computing effective rate for Gmu=%e, epsilon=%e, p=%e\n ",Gmu, epsilon, p);
		    
		    /* Compute the rate of bursts */
		    gammaAverage=0.0;
		    gammaMin = 0.0;
		    gammaMax = 0.0;
		    for (l = 0; l < Namp-1; l++)
		      {
			Dlnz = (-log(zofA[l+1])+log(zofA[l]));
			gammaAverage += eff[l] * zofA[l] * dRdz[l] * Dlnz;
			gammaMin += fmaxf((eff[l]-Deff[l]),0.0) * zofA[l] * dRdz[l] * Dlnz;
			gammaMax += fminf((eff[l]+Deff[l]),1.0) * zofA[l] * dRdz[l] * Dlnz;
		      }
		    gammaAverage *= c/p;
		    gammaMin *= c/p;
		    gammaMax *= c/p;

		    fprintf( fp,"%e  %e  %e  %e  %e  %e  %e\n", p,CLA.n,epsilon,Gmu,gammaAverage,gammaMin,gammaMax);

		  }
				
	      }
	  }

	fclose( fp );

	free(amp);
	free(lnamp);
	free(eff);
	free(Deff);
	free(zofA);
	free(dzdA);
	free(dRdz);

	return 0;
}

/*******************************************************************************/
int finddRdz(double Gmu, double alpha, double f, double Gamma)
{
  cs_cosmo_functions_t cosmofns;
  int j;

  cosmofns = XLALCSCosmoFunctions( zofA, (size_t) Namp);
  
  for ( j = 0; j < Namp; j++ )
    {


      double theta = pow((1+cosmofns.z[j]) * f * alpha * cosmofns.phit[j] / H0, -1.0/3.0);
      
      if (theta > 1.0)
	{
	  dRdz[j] = 0.0;
	}
      else
	{
	
	  dRdz[j] = 0.5 * H0 * pow(f/H0,-2.0/3.0) * pow(alpha, -5.0/3.0) / (Gamma*Gmu) *
	    pow(cosmofns.phit[j],-14.0/3.0) * cosmofns.phiV[j] * pow(1+cosmofns.z[j],-5.0/3.0);
	}
/*       fprintf(stdout,"%e %e\n", cosmofns.z[j], dRdz[j]); */

    }


  XLALCSCosmoFunctionsFree( cosmofns );

  return 0;
}
/*******************************************************************************/
int findzofA(double Gmu, double alpha)
{
  double lnz_min = log(1e-20), lnz_max = log(1e10), dlnz =0.05;
  size_t numz       = floor( (lnz_max - lnz_min)/dlnz );
  int i,j;
  cs_cosmo_functions_t cosmofns;
  double *fz,*z;
  double a;
  gsl_interp *zofa_interp; 
  gsl_interp_accel *acc_zofa = gsl_interp_accel_alloc(); 

  cosmofns = XLALCSCosmoFunctionsAlloc( exp( lnz_min ), dlnz, numz );

  zofa_interp = gsl_interp_alloc (gsl_interp_linear, cosmofns.n);

  fz   = calloc( cosmofns.n, sizeof( *fz ) ); 
  z   = calloc( cosmofns.n, sizeof( *z ) ); 
  
  /* first compute the function that relates A and z */
  /* invert order; b/c fz is a monotonically decreasing func of z */
  j=0;
  for ( i = cosmofns.n-1 ; i >=  0; i-- )
    {
      z[j]=cosmofns.z[i];
      fz[j] = pow(cosmofns.phit[i],2.0/3.0) * pow(1+z[j],-1.0/3.0) / cosmofns.phiA[i];
      j=j+1;
    }

  gsl_interp_init (zofa_interp, fz, z, cosmofns.n);

  /* now compute the amplitudes (suitably multiplied) that are equal to fz for some z*/
  for ( j = 0; j < Namp; j++ )
    {
      a = amp[j] * pow(H0,-1.0/3.0) * pow(alpha,-2.0/3.0) / Gmu;
      zofA[j] = gsl_interp_eval (zofa_interp, fz, z, a, acc_zofa );
/*       fprintf(stdout,"%e %e %e\n", amp[j],a, zofA[j]); */
    }

  XLALCSCosmoFunctionsFree( cosmofns );
  free(fz);
  free(z);
  gsl_interp_free (zofa_interp);
  gsl_interp_accel_free(acc_zofa);

  return 0;
  
}

/*******************************************************************************/

int ReadEfficiencyFile(struct CommandLineArgsTag CLA)
{

  char line[256];
  
  int i=0;
  FILE *fpEff;

  fpEff=fopen(CLA.efficiencyfile,"r");
  if (fpEff==NULL)
   {
     fprintf(stderr,"Error: File %s doesn't exist!\n",CLA.efficiencyfile);
     return 1;
   }

  /* count the number of lines */
  while(fgets(line,sizeof(line),fpEff))
    {
      if(*line == '#') continue;
      if(*line == '%') continue;
      i=i+1;
    }
  fclose(fpEff);

  fpEff=fopen(CLA.efficiencyfile,"r");
  Namp=i;
  if (Namp>1)
    {
      /* Allocate amplitude, efficiency and efficieny error arrays */
      lnamp  = calloc( Namp, sizeof( *lnamp ) ); 
      eff    = calloc( Namp, sizeof( *eff ) ); 
      Deff   = calloc( Namp, sizeof( *Deff ) ); 
      amp    = calloc( Namp, sizeof( *amp ) ); 
      zofA   = calloc( Namp, sizeof( *zofA ) ); 
      dzdA   = calloc( Namp, sizeof( *dzdA ) ); 
      dRdz   = calloc( Namp, sizeof( *dRdz ) ); 

      /*read them from the file */
      i=0;
      while(fgets(line,sizeof(line),fpEff))
	{
	  if(*line == '#') continue;
	  if(*line == '%') continue;
	  sscanf(line,"%le %le %le",lnamp+i,eff+i,Deff+i);
	  amp[i]= exp(lnamp[i]);
	  i=i+1;
	}
    }
  else
    {
      fprintf(stderr,"Error: File %s does not appear to contain enough data!\n ",CLA.efficiencyfile);
      fprintf(stderr,"Must have at least two lines of data in the file\n");
      return 1;
    }
  fclose(fpEff);

  return 0;
}



/*******************************************************************************/

int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA)
{
  int errflg = 0;
  optarg = NULL;

  struct option long_options[] = {
    {"frequency",                   required_argument, NULL,           'a'},
    {"log-Gmustart",                required_argument, NULL,           'b'},
    {"log-Gmuend",                  required_argument, NULL,           'c'},
    {"nGmu",                        required_argument, NULL,           'd'},
    {"log-epsilonstart",            required_argument, NULL,           'e'},
    {"log-epsilonend",              required_argument, NULL,           'f'},
    {"nepsilon",                    required_argument, NULL,           'g'},
    {"index",                       required_argument, NULL,           'i'},
    {"log-pstart",                  required_argument, NULL,           'j'},
    {"log-pend",                    required_argument, NULL,           'k'},
    {"np",                          required_argument, NULL,           'l'},
    {"efficiency-file",             required_argument, NULL,           'm'},
    {"help",                        no_argument, NULL,                 'h'},
    {0, 0, 0, 0}
  };
  char args[] = "ha:b:c:d:e:f:g:i:j:k:l:";

  CLA->f=               -300;     
  CLA->logGmustart=     -300;      
  CLA->logGmuend=       -300;        
  CLA->nGmu=            -1;
  CLA->logepsilonstart= -300;  
  CLA->logepsilonend=   -300;
  CLA->nepsilon=        -1;
  CLA->n=               -1;
  CLA->logpstart=       -300;
  CLA->logpend=         -300;
  CLA->np=              -1;
  CLA->efficiencyfile=  NULL;
  

  /* Scan through list of command line arguments */
  while ( 1 )
  {
    int option_index = 0; /* getopt_long stores long option here */
    int c;

    c = getopt_long_only( argc, argv, args, long_options, &option_index );
    if ( c == -1 ) /* end of options */
      break;

    switch ( c )
    {

    case 'a':
      /* lowest frequency  */
      CLA->f=atof(optarg);
      break;
    case 'b':
      /* lowest frequency  */
      CLA->logGmustart=atof(optarg);
      break;
    case 'c':
      /* highest frequency */
      CLA->logGmuend=atof(optarg);
      break;
    case 'd':
      /* number of frequencies to do */
      CLA->nGmu=atoi(optarg);
      break;
    case 'e':
      /* lowest frequency  */
      CLA->logepsilonstart=atof(optarg);
      break;
    case 'f':
      /* highest frequency */
      CLA->logepsilonend=atof(optarg);
      break;
    case 'g':
      /* number of frequencies to do */
      CLA->nepsilon=atoi(optarg);
      break;
    case 'i':
      /* number of frequencies to do */
      CLA->n=atof(optarg);
      break;
    case 'j':
      /* highest frequency */
      CLA->logpstart=atof(optarg);
      break;
    case 'k':
      /* number of frequencies to do */
      CLA->logpend=atof(optarg);
      break;
    case 'l':
      /* number of frequencies to do */
      CLA->np=atoi(optarg);
      break;
    case 'm':
      /* number of frequencies to do */
      CLA->efficiencyfile=optarg;
      break;

    case 'h':
      /* print usage/help message */
      fprintf(stdout,"Arguments are :\n");
      fprintf(stdout,"\t--frequency (-a)\t\tFLOAT\t Lowest frequency.\n");
      fprintf(stdout,"\t--log-Gmustart (-b)\t\tFLOAT\t Lowest Gmu.\n");
      fprintf(stdout,"\t--log-Gmuend (-c)\t\tFLOAT\t Largest Gmu.\n");
      fprintf(stdout,"\t--nGmu (-d)\t\t\tINTEGER\t Number of Gmu bins to do.\n");
      fprintf(stdout,"\t--log-epsilonstart (-e)\t\tFLOAT\t Lowest epsilon.\n");
      fprintf(stdout,"\t--log-epsilonend (-f)\t\tFLOAT\t Largest epsilon.\n");
      fprintf(stdout,"\t--nepsilon (-g)\t\t\tINTEGER\t Number of epsilon bins to do.\n");
      fprintf(stdout,"\t--index (-i)\t\t\tFLOAT\t Index for alpha as function of Gmu.\n");
      fprintf(stdout,"\t--log-pstart (-j)\t\tFLOAT\t Lowest p.\n");
      fprintf(stdout,"\t--log-pend (-k)\t\t\tFLOAT\t Largest p.\n");
      fprintf(stdout,"\t--np (-l)\t\t\tINTEGER\t Number of p bins to do.\n");
      fprintf(stdout,"\t--efficiency-file (-m)\t\t\tSTRING\t File with efficiency values and errors.\n");
      fprintf(stdout,"\t--help (-h)\t\t\tFLAG\t Print this message.\n");
      exit(0);
      break;
    default:
      /* unrecognized option */
      errflg++;
      fprintf(stderr,"Unrecognized option argument %c\n",c);
      exit(1);
      break;
    }
    }

  if(CLA->f == -300)
    {
      fprintf(stderr,"No lowest frequency specified.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }
  if(CLA->logGmustart == -300)
    {
      fprintf(stderr,"No lowest Gmu specified.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }
  if(CLA->logGmuend == -300)
    {
      fprintf(stderr,"No highest Gmu specified.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }
  if(CLA->nGmu == -1)
    {
      fprintf(stderr,"No number of Gmu bins specified.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }
  if(CLA->logepsilonstart == -300)
    {
      fprintf(stderr,"No lowest epsilon specified.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }
  if(CLA->logepsilonend == -300)
    {
      fprintf(stderr,"No highest epsilon specified.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }
  if(CLA->nepsilon == -1)
    {
      fprintf(stderr,"No number of epsilon bins specified.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }

  if(CLA->logpstart == -300)
    {
      fprintf(stderr,"No lowest p specified.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }
  if(CLA->logpend == -300)
    {
      fprintf(stderr,"No highest p specified.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }
  if(CLA->np == -1)
    {
      fprintf(stderr,"No number of p bins specified.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }
  if(CLA->n == -1)
    {
      fprintf(stderr,"No index for alpha specified.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }
  if(CLA->efficiencyfile == NULL)
    {
      fprintf(stderr,"No efficiency file specified.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }

  return errflg;
}

