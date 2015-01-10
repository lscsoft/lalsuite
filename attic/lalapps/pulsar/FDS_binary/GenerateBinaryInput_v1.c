/*
*  Copyright (C) 2007 Chris Messenger
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
 * \author Chris Messenger
 * \brief Program to generate an input file for lalapps_makefakedata
 */

/*********************************************************************************/
/*         Program to generate an input file for lalapps_makefakedata            */
/*                                                                               */
/*			           C. Messenger                                  */
/*                                                                               */
/*                         BIRMINGHAM UNIVERISTY -  2004                         */
/*********************************************************************************/

#include <errno.h>
#include <unistd.h>
#include <sys/types.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <glob.h>
#include <time.h>
#include <getopt.h>
#include <lal/LALDatatypes.h>
#include "ReadSourceFile_v1.h"


REAL8 alpha,delta;
REAL8 sma,period,ecc,argperi;
INT4 tperisec,tperins;
INT4 tsft,start,starttime,reftime;
REAL8 f_min,band,sigma,duration;
CHAR noisedir[256],efiles[56],basename[256],yr[256],ifo[256],stamps[256],outfile[256];
REAL8 phi,psi,cosiota,f0,h0;
REAL8 f1dot,f2dot,f3dot;
REAL8 aplus,across;
BOOLEAN stampsflag=0;
BOOLEAN startflag=0;
BOOLEAN noisedirflag=0;
CHAR sourcefile[256],source[256];
BOOLEAN sourceflag;

extern char *optarg;
extern int optind, opterr, optopt;
int GenTimeStamps(void);
int ReadCommandLine(int argc,char *argv[]); 
int OutputConfigFile(void);

int main(int argc,char *argv[]) 
{

  binarysource sourceparams;
  LIGOTimeGPS *dummyGPS=NULL;

  if (ReadCommandLine(argc,argv)) return 1;
  
  if (sourceflag) {
    if (ReadSource(sourcefile,source,dummyGPS,&sourceparams)) return 2;
    alpha = sourceparams.skypos.ra;
    delta = sourceparams.skypos.dec;
  }

  if ((stampsflag)&&(startflag)) {
    if (GenTimeStamps()) return 2;
  }

  if (OutputConfigFile()) return 3;
  
  return 0;

}

/*******************************************************************************/

int GenTimeStamps() 
{

  FILE *fpstamps;
  INT4 j;
  INT4 nsft;

  nsft=floor(duration/tsft);

  /* opening the output timestamps file */
  fpstamps=fopen(stamps,"w");
  if (fpstamps==NULL) {
    fprintf(stderr,"Unable to open file %s\n",stamps);
    return 1;
  }
  
  /* write out timestamps file */
  for (j=0;j<nsft;j++) {
    fprintf(fpstamps,"%d\t0\n",starttime+(j*tsft));
  }

  fclose(fpstamps);

  return 0;

}

/*******************************************************************************/
   int ReadCommandLine(int argc,char *argv[]) 
{
  INT4 c, errflg = 0;
  char *temp;
  optarg = NULL;
  
  /* Initialize default values */
  sprintf(sourcefile," ");
  sprintf(source," ");
  sourceflag=0;
  alpha=0.0; /* a */
  delta=0.0; /* d */
  sprintf(ifo,"LLO");      /* I */            
  phi=0.0; /* p */
  psi=0.0; /* P */
  cosiota=0.0; /* c */
  h0=1.0; /* H */
  f0=600.0; /* F */
  tsft=60; /* t */
  sprintf(stamps," "); /* T */
  duration=0.0; /* N */
  f_min=540.0; /* f */
  band=100.0; /* b */
  starttime=0; /* S */
  reftime=0; /* r */
  sigma=0.0; /* w */
  
  sprintf(efiles,"./");/* g */  
  sprintf(yr,"00-04"); /* k */
  sprintf(noisedir," "); /* v */
  sprintf(basename," "); /* u */
  sprintf(outfile,"out.cfg"); /* o */

  f1dot=0.0; /* 1 */
  f2dot=0.0; /* 2 */
  f3dot=0.0; /* 3 */
 
  sma=0.0; /* R */
  ecc=0.0; /* e */
  tperisec=0; /* Q */
  tperins=0; /* q */
  period=0.0; /* O */
  argperi=0.0; /* A */

  {
    int option_index = 0;
    static struct option long_options[] = {
      {"sourcefile", required_argument, 0, 'Q'},
      {"source", required_argument, 0, 'q'},
      {"alpha", required_argument, 0, 'a'},
      {"delta", required_argument, 0, 'd'},
      {"ifo", required_argument, 0, 'I'},
      {"phi", required_argument, 0, 'p'},
      {"psi", required_argument, 0, 'P'},
      {"cosiota", required_argument, 0, 'c'},
      {"h0", required_argument, 0, 'H'},
      {"f0", required_argument, 0, 'F'},
      {"tsft", required_argument, 0, 't'},
      {"duration", required_argument, 0, 'N'},
      {"f_min", required_argument, 0, 'f'},
      {"band", required_argument, 0, 'b'},
      {"stamps", required_argument, 0, 's'},
      {"start", required_argument, 0, 'S'},
      {"reftime", required_argument, 0, 'R'},
      {"sigma", required_argument, 0, 'g'},
      {"ephem", required_argument, 0, 'E'},
      {"yr", required_argument, 0, 'y'},
      {"noisedir", required_argument, 0, 'n'},
      {"basename", required_argument, 0, 'm'},
      {"f1dot", required_argument, 0, '1'},
      {"f2dot", required_argument, 0, '2'},
      {"f3dot", required_argument, 0, '3'},
      {"smaxis", required_argument, 0, 'A'},
      {"period", required_argument, 0, 'O'},
      {"tperisec", required_argument, 0, 'X'},
      {"tperinan", required_argument, 0, 'x'},
      {"ecc", required_argument, 0, 'e'},
      {"argperi", required_argument, 0, 'u'},
      {"outfile", required_argument, 0, 'o'},
      {"help", required_argument, 0, 'h'}
    };
  /* Scan through list of command line arguments */
  while (!errflg && ((c = getopt_long (argc, argv,"hQ:q:a:d:I:p:P:c:H:F:t:N:f:b:s:S:R:g:E:y:n:m:1:2:3:A:O:X:x:e:u:o:",long_options, &option_index)))!=-1)
    switch (c) {
    case 'Q':
      temp=optarg;
      sprintf(sourcefile,"%s",temp);
      sourceflag=1;
      break;
    case 'q':
      temp=optarg;
      sprintf(source,"%s",temp);
      break;
    case 'a':
      alpha=atof(optarg);
      break;
    case 'd':
      delta=atof(optarg);
      break;
    case 'I':
      temp=optarg;
      sprintf(ifo,"%s",temp);
      break;
    case 'p':
      phi=atof(optarg);
      break;
    case 'P':
      psi=atof(optarg);
      break;
    case 'c':
      cosiota=atof(optarg);
      break;
    case 'H':
      h0=atof(optarg);
      break;
    case 'F':
      f0=atof(optarg);
      break;
    case 't':
      tsft=atoi(optarg);
      break;
    case 'N':
      duration=atof(optarg);
      break;
    case 'f':
      f_min=atof(optarg);
      break;
    case 'b':
      band=atof(optarg);
      break;
    case 's':
      temp=optarg;
      sprintf(stamps,"%s",temp);
      stampsflag=1;
      break;
    case 'S':
      start=atoi(optarg);
      break;
    case 'R':
      reftime=atoi(optarg);
      break;
    case 'g':
      sigma=atof(optarg);
      break;
    case 'E':
      temp=optarg;
      sprintf(efiles,"%s",temp);
      break;
    case 'y':
      temp=optarg;
      sprintf(yr,"%s",temp);
      break;
    case 'n':
      temp=optarg;
      sprintf(noisedir,"%s",temp);
      noisedirflag=1;
      break;
    case 'm':
      temp=optarg;
      sprintf(basename,"%s",temp);
      break;
    case '1':
      f1dot=atof(optarg);
      break;
    case '2':
      f2dot=atof(optarg);
      break;
    case '3':
      f3dot=atof(optarg);
      break;
    case 'A':
      sma=atof(optarg);
      break;
    case 'O':
      period=atof(optarg);
      break;
    case 'X':
      tperisec=atoi(optarg);
      break;
    case 'x':
      tperins=atoi(optarg);
      break;
    case 'e':
      ecc=atof(optarg);
      break;
    case 'u':
      argperi=atof(optarg);
      break;
    case 'o':
      temp=optarg;
      sprintf(outfile,"%s",temp);
      break;
    case 'h':
      /* print usage/help message */
      fprintf(stdout,"Arguments are:\n");
      fprintf(stdout,"\t--sourcefile STRING\t The name of the source file conotaining source information [DEFAULT= ]\n");
      fprintf(stdout,"\t--source     STRING\t The name of the source [DEFAULT= ]\n");
      fprintf(stdout,"\t--alpha    FLOAT\t Sky position alpha (equatorial coordinates) in radians [DEFAULT= ]\n");
      fprintf(stdout,"\t--delta    FLOAT\t Sky position delta (equatorial coordinates) in radians [DEFAULT= ]\n");
      fprintf(stdout,"\t--ifo      STRING\t Detector (LLO,LHO,GEO,VIRGO,TAMA) [DEFAULT=LLO]\n");
      fprintf(stdout,"\t--phi      FLOAT\t Initial phase at Tref in radians [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--psi      FLOAT\t Polarisation angle in radians [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--cosiota  FLOAT\t Cos(iota) [DEFAULT 0.0]\n");
      fprintf(stdout,"\t--h0       FLOAT\t Gravitational wave amplitude [DEFAULT=1.0]\n");
      fprintf(stdout,"\t--f0       FLOAT\t Gravitational wave frequency in Hz [DEFAULT=600.0]\n");
      fprintf(stdout,"\t--tsft     INTEGER\t Time basefile of SFT's in seconds [DEFAULT=60]\n");
      fprintf(stdout,"\t--duration REAL8\t Number of SFT's to generate [DEFAULT=0]\n");
      fprintf(stdout,"\t--f_min     FLOAT\t Minimum generation frequency in Hz [DEFAULT=40.0] \n");
      fprintf(stdout,"\t--band     FLOAT\t Bandwidth to be generated in Hz [DEFAULT=10.0]\n");
      fprintf(stdout,"\t--stamps   STRING\t Location and name of timestamps file [DEFAULT=NULL]\n");
      fprintf(stdout,"\t--start    INTEGER\t GPS Start time of observation [DEFAULT=0]\n");
      fprintf(stdout,"\t                  \t (If both stamps and start are set then a stamps file is generated)\n");
      fprintf(stdout,"\t--reftime  INTEGER\t SSB time at which pulsar parameters are specified [DEFAULT=0]\n");
      fprintf(stdout,"\t--sigma    FLOAT\t Variance of Guassian noise to be added (if -ve use real data) [DEFAULT=0]\n");
      fprintf(stdout,"\t--ephem    STRING\t Location of ephemeris data [DEFAULT=./]\n");
      fprintf(stdout,"\t--yr       STRING\t Year of ephemeris to obe read [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--noisedir STRING\t Directory containing data into which the sigmal is injected [DEFAULT=NULL]\n");
      fprintf(stdout,"\t--basename STRING\t Location and basename of output SFT's [DEFAULT=NULL]\n");
      fprintf(stdout,"\t--f1dot    FLOAT\t Spin down parameter (df0/dt) [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--f2dot    FLOAT\t Spin down parameter (d^2f0/dt^2) [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--f3dot    FLOAT\t Spin down parameter (d^3f0/dt^3) [DEFAULT=0.0] \n");
      fprintf(stdout,"\t--smaxis   FLOAT\t Projected semi-major axis of orbit in seconds [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--period   FLOAT\t Period of orbit in seconds [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--tperisec INTEGER\t Observed time of periapse passage in seconds [DEFAULT=0] \n");
      fprintf(stdout,"\t--tperinan INTEGER\t Observed time of periapse passage in nanoseconds [DEFAULT=0]\n");
      fprintf(stdout,"\t--ecc      FLOAT\t Orbital eccentricity [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--argperi  FLOAT\t Argument of orbital periapse in radians [DEFAULT=0.0]\n");
      fprintf(stdout,"\t--outfile  STRING\t Name of output configuration file [DEFAULT=NULL]\n");
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
  
  aplus=0.5*h0*(1.0*cosiota*cosiota);
  across=h0*cosiota;
  
  /* update global variable and return */
  return errflg;
}


/*******************************************************************************/


/*******************************************************************************/

int OutputConfigFile() 
{
  FILE *fp;

  /*   %strcpy(filename,inDataFilename); */
  fp=fopen(outfile,"w");
  if (fp==NULL) {
    fprintf(stderr,"Unable to open file %s\n",outfile);
    return 1;
  }

  fprintf(fp,"## Settings generated by lalapps_GenBinaryInput\n");
  fprintf(fp,"##--------------------------------------------------\n");
  fprintf(fp,"\n");
  fprintf(fp,"##--------------------------------------------------\n");
  fprintf(fp,"## REQUIRED user variables\n");
  fprintf(fp,"##--------------------------------------------------\n");
  fprintf(fp,"\n");
  fprintf(fp,"detector\t= %s\t\t\t# Detector: LHO, LLO, VIRGO, GEO, TAMA, CIT, ROME\n",ifo);
  fprintf(fp,"ephemDir\t= %s\t\t\t# directory containing ephemeris files\n",efiles);
  if (!stampsflag) {
    fprintf(fp,"duration\t\t= %f\t\t\t# Duration of requested signal in seconds\n",duration);
  }
  else {
    fprintf(fp,"#duration\t\t= \t\t\t# Duration of requested signal in seconds\n");
  }
  fprintf(fp,"fmin\t\t= %6.12f\t# lowest SFT-frequency in Hz\n",f_min);
  fprintf(fp,"Band\t\t= %6.12f\t# SFT frequency band in Hz\n",band);
  fprintf(fp,"longitude\t= %6.12f\t# source longitude (in_radians)\n",alpha);
  fprintf(fp,"latitude\t= %6.12f\t# source latitude (in radians)\n",delta);
  fprintf(fp,"aPlus\t\t= %6.12f\t# plus-polarization amplitude a_+ (strain)\n",aplus);
  fprintf(fp,"aCross\t\t= %6.12f\t# cross-polarization amplitude a_x (strain)\n",across);
  fprintf(fp,"psi\t\t= %6.12f\t# wave polarization angle Psi\n",psi);
  fprintf(fp,"phi0\t\t= %6.12f\t# initial wave-phase phi0 (at reference-time tRef)\n",phi);
  fprintf(fp,"f0\t\t= %6.12f\t# intrinsic signal frequency f0 (at tRef)\n",f0);
  fprintf(fp,"\n");
  fprintf(fp,"##--------------------------------------------------\n");
  fprintf(fp,"## OPTIONAL user variables, defaults are shown here (if any)\n");
  fprintf(fp,"##--------------------------------------------------\n");
  fprintf(fp,"\n");
  fprintf(fp,"Tsft\t\t= %d\t\t\t# length of SFTs in seconds\n",tsft);
  fprintf(fp,"ephemYear\t= %s\t\t\t# ephemeris years to be read\n",yr);
  fprintf(fp,"f1dot\t\t= %6.12f\t# spindowns: d/dt f0\n",f1dot);
  fprintf(fp,"f2dot\t\t= %6.12f\t# d^2/dt^2 f0\n",f2dot);
  fprintf(fp,"f3dot\t\t= %6.12f\t# d^3/dt^3 f0\n",f3dot);
  fprintf(fp,"\n");
  fprintf(fp,"## --- exactly ONE of the following two has to be specified\n");
  if (stampsflag) {
    fprintf(fp,"#startTime\t=\t\t# GPS start time of (contiguous) output time-series\n");
    fprintf(fp,"timestampsFile\t\t= %s\t\t# file containing timestamps to produce (nSFT) SFTs\n",stamps);
  }
  else {
    fprintf(fp,"startTime\t= %d\t\t\t# GPS start time of (contiguous) output time-series\n",start);
    fprintf(fp,"#timestampsFile\t= \t\t# file containing timestamps to produce (nSFT) SFTs\n");
  }
  fprintf(fp,"##--------------------------------------------------\n");
  fprintf(fp,"\n");
  if (reftime<=0) 
    {
      fprintf(fp,"#refTime\t=\t # reference-time tRef in SSB of pulsar-parameters\n");
    }
  else 
    {
      fprintf(fp,"refTime\t= %d\t# reference-time tRef in SSB of pulsar-parameters\n",reftime);
    }
  fprintf(fp,"\n");
  fprintf(fp,"## --- maximally ONE of the following two can be specified\n");
  if (noisedirflag)
    {
      fprintf(fp,"#noiseSigma\t= \t# variance for Gaussian noise to be added to output\n");
      fprintf(fp,"noiseSFTs\t= %s\t# Directory with real noise SFTs to be added\n",noisedir);
    }
  else
    {
      fprintf(fp,"noiseSigma\t= %6.12f\t# variance for Gaussian noise to be added to output\n",sigma);
      fprintf(fp,"#noiseSFTs\t= \t# Directory with real noise SFTs to be added\n");
    }
   fprintf(fp,"\n");
   fprintf(fp,"outSFTbname\t= %s\t# path+basename for final SFTs, e.g. 'test_new/SFT'\n",basename);
   if (sma>0.0) 
     {
       fprintf(fp,"\n");
       fprintf(fp,"##------------------------------------------------------------\n");
       fprintf(fp,"## ORBITAL parameters for binary neutron stars [OPTIONAL]\n");
       fprintf(fp,"## NOTE: you need to specify either NONE or ALL of these (except if SemiMajorAxis==0)\n");
       fprintf(fp,"##------------------------------------------------------------\n");
       fprintf(fp,"\n");
       fprintf(fp,"orbitSemiMajorAxis\t= %6.12f\t# Projected orbital semi-major axis a in seconds (i.e. a*sin(i)/c)\n",sma);
       fprintf(fp,"orbitEccentricity\t= %6.12f\t# Orbital eccentricity\n",ecc);
       fprintf(fp,"orbitTperiSSBsec\t= %d\t\t# 'observed' (SSB) time of periapsis passage. Seconds.\n",tperisec);
       fprintf(fp,"orbitTperiSSBns\t\t= %d\t\t# 'observed' (SSB) time of periapsis passage. Nanoseconds.\n",tperins);
       fprintf(fp,"orbitPeriod\t\t= %6.12f\t# Orbital period (seconds)\n",period);
       fprintf(fp,"orbitArgPeriapse\t= %6.12f\t# Argument of periapsis (radians)\n",argperi); 
     }

  fclose(fp);
  return 0;
  
}
/*******************************************************************************/

