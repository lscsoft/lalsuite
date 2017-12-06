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
/*            Cosmic string burst rate computation code for large loops          */
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
  double logpstart;        
  double logpend;          
  int  np;
  char *efficiencyfile;
} CLA;

double *eff, *lnamp, *amp, *Deff, Dlnz;  /* arrays for efficiency, amplitude and efficiency errors */
double H0 = LAMBDA_H_0, DlnA;
int Namp;   /* size of amplitude/eff array */


#define PROGRAM_NAME "cs_gammaLargeLoops"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"

/****************************************************************************************/
/* Reads the command line */
int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);
int ReadEfficiencyFile(struct CommandLineArgsTag CLA);
double finddRdzdA(double Gmu, double f, double Gamma, double A, double z, double phit, double phiV, double phiA);
double nu(double Gmu, double Gamma, double A, double z, double phit, double phiA);

/************************************* MAIN PROGRAM *************************************/

int main( int argc, char *argv[] )
{
	/* parameters of strings */
	double Gamma   = LOOP_RAD_POWER;	/* power radiated by loops */
	double c       = CUSPS_PER_LOOP;	/* cusps per loop */
	double RAve, RMin, RMax;

	/* limits and resolution for integral over z */
	double lnz_min    = -25.0;
	double lnz_max    =  25.0;
	double dlnz       = 0.1;
	size_t numz       = floor( (lnz_max - lnz_min)/dlnz );

	double p,f,Gmu;
	char filename[FILENAME_MAX];
	FILE *fp;
	int i,j,k;
  
	cs_cosmo_functions_t cosmofns = XLALCSCosmoFunctionsAlloc( exp( lnz_min ), dlnz, numz );

	/* read the command line */
	if (ReadCommandLine(argc,argv,&CLA)) return 1;
	f=CLA.f;

	/* read efficiency function */
	if (ReadEfficiencyFile(CLA)) return 2;
	
	snprintf( filename, sizeof( filename ), "gammaLL.dat");
	fp = fopen( filename, "w" );
	fprintf( fp,"%%     p           Gmu       gammaAverage     gammaMin      gammaMax\n");

	for ( i = 0; i <  CLA.nGmu; i++ )
	  {
	    Gmu = pow(10.0,CLA.logGmustart+i*(CLA.logGmuend-CLA.logGmustart)/(CLA.nGmu-1));
	    	       
	    RAve = 0.0;
	    RMin = 0.0;
	    RMax = 0.0;
	    /* loop over amplitudes */
	    for ( j = 0; j < Namp; j++ )
	      {
		double dRdA = 0.0;
		/* loop over redshifts */
		for ( k = 0; k < (int) cosmofns.n; k++ )
		  {
		    double z = cosmofns.z[k];
		    double phit = cosmofns.phit[k];
		    double phiA = cosmofns.phiA[k];
		    double phiV = cosmofns.phiV[k];
		    double A = amp[j];

		    double dRdzdA = finddRdzdA(Gmu, f, Gamma, A, z, phit, phiV, phiA);

		    dRdA += z*dRdzdA*dlnz;
		  }
		RAve +=      eff[j]      * amp[j] * dRdA * DlnA;
		RMin +=  fmaxf((eff[j]-Deff[j]),0.0) * amp[j] * dRdA * DlnA;
		RMax +=  fminf((eff[j]+Deff[j]),1.0) * amp[j] * dRdA * DlnA;
	      }

	    for ( k = 0; k <  CLA.np; k++ )
	      {
		p=pow(10.0, CLA.logpstart+k*(CLA.logpend-CLA.logpstart)/(CLA.np));    
		fprintf(stdout,"%%Computing effective rate for Gmu=%e, p=%e\n ",Gmu, p);
		fprintf( fp,"%e  %e  %e  %e  %e\n", p, Gmu,RAve*c/p,RMin*c/p,RMax*c/p);
	      }
	  }
		
	fclose( fp );

	free(amp);
	free(lnamp);
	free(eff);
	free(Deff);
        XLALCSCosmoFunctionsFree( cosmofns );

	return 0;
}

/*******************************************************************************/
double finddRdzdA(double Gmu, double f, double Gamma, double A, double z, double phit, double phiV, double phiA)
{

  double l = pow ( A / Gmu / H0 * pow(1+z,1.0/3.0)* phiA, 3.0/2.0);
  double theta = pow(f*(1+z)*l, -1.0/3.0);
  double Delta;
  double dRdzdA;

  if (theta > 1.0)
    return 0.0;

  Delta = 0.25*theta*theta;

  dRdzdA = pow(H0,-3.0) * phiV / (1.0+z) * nu(Gmu, Gamma, A, z, phit, phiA) * Delta;

  return dRdzdA;

}

/*******************************************************************************/
double nu(double Gmu, double Gamma, double A, double z, double phit, double phiA)
{
  double l = pow( A / Gmu / H0 * pow(1.0+z,1.0/3.0)* phiA, 3.0/2.0);
  
  /* Alex's loop distribution */
  double alpha=1e-1;
  double nuR=0.4*15*sqrt(alpha);
  double nuM=0.12*4;

  double t=phit/H0;
  double cuspdist;

  double teq=8.122570474611143e+11;
  double crateR, crateRadStragglers, crateM;

  /* Radiation era loops */
  crateR = 0.0;
  if( (l < alpha*t) && (t < teq) )
    crateR = nuR * pow(t,-1.5) * pow( l + Gamma*Gmu/H0 * phit, -2.5 );
    
  /* Radiation stragglers */
  crateRadStragglers = 0.0; 
  if ( (l <  alpha*teq-Gamma*Gmu*(t-teq) ) && ( t > teq ) )
    crateRadStragglers = nuR*pow(teq, 0.5)/t/t* pow( l + Gamma*Gmu/H0 * phit, -2.5);

  /* matter era loops */
  crateM = 0.0;
  if( (l < alpha*t) && (t > teq) && (l >  alpha*teq-Gamma*Gmu*(t-teq)) )
    crateM = nuM / t / t / ( l + Gamma*Gmu/H0 * phit) / ( l + Gamma*Gmu/H0 * phit);
 
  cuspdist = 3.0/A * (crateR+crateRadStragglers+crateM);

  return cuspdist;
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

      /*read them from the file */
      i=0;
      while(fgets(line,sizeof(line),fpEff))
	{
	  if(*line == '#') continue;
	  if(*line == '%') continue;
	  sscanf(line,"%le %le %le",lnamp+i,eff+i,Deff+i);
	  amp[i]= exp(lnamp[i]);

/* 	  fprintf(stdout,"%e %e\n", amp[i], eff[i]); */

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

  /* Define resolution in amplitude */
  DlnA = lnamp[1]-lnamp[0];

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
    {"log-pstart",                  required_argument, NULL,           'e'},
    {"log-pend",                    required_argument, NULL,           'f'},
    {"np",                          required_argument, NULL,           'g'},
    {"efficiency-file",             required_argument, NULL,           'i'},
    {"help",                        no_argument, NULL,                 'h'},
    {0, 0, 0, 0}
  };
  char args[] = "ha:b:c:d:e:f:g:i:";

  CLA->f=               -300;     
  CLA->logGmustart=     -300;      
  CLA->logGmuend=       -300;        
  CLA->nGmu=            -1;
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
      /* highest frequency */
      CLA->logpstart=atof(optarg);
      break;
    case 'f':
      /* number of frequencies to do */
      CLA->logpend=atof(optarg);
      break;
    case 'g':
      /* number of frequencies to do */
      CLA->np=atoi(optarg);
      break;
    case 'i':
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
  if(CLA->efficiencyfile == NULL)
    {
      fprintf(stderr,"No efficiency file specified.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }

  return errflg;
}

