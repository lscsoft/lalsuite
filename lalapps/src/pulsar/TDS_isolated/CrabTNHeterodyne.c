/*
*  Copyright (C) 2007 Matt Pitkin
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
 * \author Matt Pitkin
 * \brief
 * This driver code performs the timing noise heterodyne for 
 * the Crab pulsar on completely heterodyned data.
 * No calibration or noise estimation is needed or performed
 */

/* This driver code performs the timing noise heterodyne for 
   the Crab pulsar on completely heterodyned data.
   No calibration or noise estimation is needed or performed */
 
/* Matt Pitkin CrabTNHeterodyne.c v0.1 */

/* Compile with e.g.
 * gcc CrabTNHeterodyne.c HeterodyneCrabPulsar.c -o CrabTNHeterodyne -L${LAL_PREFIX}/lib -I${LAL_PREFIX}/include -lm -llalsupport -llal -llalpulsar -llalframe -llalxml -llalmetaio -llalsimulation -llalburst -llalinspiral -llalinference -llalstochastic -std=c99
 */

/* headers */
#include <stdio.h>
#include <string.h>
#include <math.h>

/* Crab heterodyne specific headers */
#include "HeterodyneCrabPulsar.h"


#define MAXLENGTH 10000010
/* max num of lines in crab ephem file (ok for next 61 years file contains 266 lines as of 15 Jan
2004) */
#define NUM 1000

#define EARTHFILE "earth00-19-DE200.dat.gz"
#define SUNFILE "sun00-19-DE200.dat.gz"
#define TIMEFILE "tdb_2000-2019.dat.gz"

int main(int argc, char *argv[]){
  static LALStatus status;
  FILE *fpin=NULL, *fpout=NULL;
  CHAR inputFile[256];
  CHAR outputFile[256];
  
  FILE *phifp=NULL;
  
  CHAR psrInput[256]; /* pulsar input file containing params f0, f1 etc. */
  
  UINT4 i=0, j; /* counter */
  REAL8Vector *times=NULL;
  COMPLEX16Vector *B=NULL;
  
  /* vars for Crab */
  GetCrabEphemerisInput input;
  CrabSpindownParamsInput crabEphemerisData;
  CrabSpindownParamsOutput crabOutput;
  ParamsForHeterodyne hetParams;
  TNHeterodyneInput TNInput;
  TNHeterodyneOutput TNOutput;
  
  /* vars for solar system ephemeris */
  EarthState earth;
  EphemerisData *edat=NULL;
  TimeCorrectionData *tdat=NULL;
  BarycenterInput baryinput;
  CHAR detName[5];
  LALDetector det;
  CHAR *lalpath = NULL, *lalpulsarpath = NULL;
  CHAR *efile=NULL, *sfile=NULL, *tfile=NULL;
    
  /* pulsar params */
  BinaryPulsarParams pulsarParams;
  
  LIGOTimeGPS dataEpoch;
  
  /* check command line inputs */
  if(argc!=6){
    fprintf(stderr, "Wrong number of input params:\n\tparfile inputfile\
 outputfile detector ephemerisfile\n");
    return 0;
  }
  
  sprintf(psrInput, "%s", argv[1]);
  sprintf(inputFile, "%s", argv[2]);
  sprintf(outputFile, "%s", argv[3]);
  sprintf(detName, "%s", argv[4]);
  input.filename = argv[5]; /* ephemeris file */
   
  /* read in Crab pulsar parameters used in heterodyning */
  XLALReadTEMPOParFile(&pulsarParams, psrInput);

  /* allocate memory for crab ephemeris */
  crabEphemerisData.f1 = NULL;
  crabEphemerisData.f1 = XLALCreateREAL8Vector(NUM);
  
  crabEphemerisData.f0 = NULL;
  crabEphemerisData.f0 = XLALCreateREAL8Vector(NUM);
  
  crabEphemerisData.tArr = NULL;
  crabEphemerisData.tArr = XLALCreateREAL8Vector(NUM);
  
  crabOutput.tArr = NULL;
  crabOutput.tArr = XLALCreateREAL8Vector(NUM);
      
  crabOutput.f0 = NULL;
  crabOutput.f0 = XLALCreateREAL8Vector(NUM);
  
  crabOutput.f1 = NULL;
  crabOutput.f1 = XLALCreateREAL8Vector(NUM);
  
  crabOutput.f2 = NULL;
  crabOutput.f2 = XLALCreateREAL8Vector(NUM);
  
  crabOutput.f3 = NULL;
  crabOutput.f3 = XLALCreateREAL8Vector(NUM);
  
  crabOutput.f4 = NULL;
  crabOutput.f4 = XLALCreateREAL8Vector(NUM);

  /* read Crab ephemeris data */
  LALGetCrabEphemeris( &status, &crabEphemerisData, &input );
  
  /* compute frequency derivatives */
  LALComputeFreqDerivatives ( &status, &crabOutput, &crabEphemerisData );
  
  /* allocate mem for input data */
  times = XLALCreateREAL8Vector(MAXLENGTH);
  B = XLALCreateCOMPLEX16Vector(MAXLENGTH);
  
  /* read in input data */
  fpin = fopen(inputFile, "r");
  
  i=0;
  while(!feof(fpin)){
    if(i >= MAXLENGTH){
      fprintf(stderr, "MAXLENGTH at %d is not big enough.\n", MAXLENGTH);
      return 1;
    }
    
    REAL8 tmpr, tmpi;
    
    if( fscanf(fpin, "%lf%lf%lf", &times->data[i], &tmpr, &tmpi) != 3 ){
      fprintf(stderr, "Error... problem reading line %d of the heterodyne file!\n", i+1);
      exit(1);
    }
    B->data[i] = tmpr + I*tmpi;
    i++;
  }
  
  fclose(fpin);
  
  /* set input params for timing noise heterodyne (2*pulsar spin freq) */
  TNInput.f0 = 2.0*pulsarParams.f0;
  TNInput.f1 = 2.0*pulsarParams.f1;
  TNInput.f2 = 2.0*pulsarParams.f2;
  TNInput.t0 = pulsarParams.pepoch;
  
  fprintf(stderr, "f0 = %lf, f1 = %le, f2 = %le, epoch = %lf.\n", TNInput.f0,
TNInput.f1, TNInput.f2, TNInput.t0); 
  
  fpout = fopen(outputFile, "w");
  phifp = fopen("DPhase.txt", "w");
  
  /* set SSB ephemeris files */
  lalpath = getenv("LALPULSAR_PREFIX");
  lalpulsarpath = XLALStringDuplicate( lalpath );
  lalpulsarpath = XLALStringAppend(lalpulsarpath, "/share/lalpulsar/");
  efile = XLALStringDuplicate(lalpulsarpath);
  sfile = XLALStringDuplicate(lalpulsarpath);
  tfile = XLALStringDuplicate(lalpulsarpath);
  
  efile = XLALStringAppend(efile, EARTHFILE);
  sfile = XLALStringAppend(sfile, SUNFILE);
  tfile = XLALStringAppend(tfile, TIMEFILE);
  
  edat = XLALInitBarycenter(efile, sfile);
  tdat = XLALInitTimeCorrections(tfile);
  
  det = *XLALGetSiteInfo( detName );

  baryinput.dInv = 0.;
  baryinput.site.location[0] = det.location[0]/LAL_C_SI;
  baryinput.site.location[1] = det.location[1]/LAL_C_SI;
  baryinput.site.location[2] = det.location[2]/LAL_C_SI;

  /* perform timing noise heterodyne */
  for(j=0;j<i-1;j++){
    REAL8 dtpos = times->data[j] - TNInput.t0;

    dataEpoch.gpsSeconds = (INT8)floor(times->data[j]);
    dataEpoch.gpsNanoSeconds = 0;
    
    TNInput.epoch.gpsSeconds = dataEpoch.gpsSeconds;
    TNInput.epoch.gpsNanoSeconds = 0;
    
    TNInput.Vh = B->data[j];
    
    /* set values for timing noise removal */
    LALSetSpindownParams( &status, &hetParams, &crabOutput, dataEpoch );

    baryinput.tgps.gpsSeconds = TNInput.epoch.gpsSeconds;
    baryinput.tgps.gpsNanoSeconds = TNInput.epoch.gpsNanoSeconds;

    /* set up RA, DEC, and distance variables for LALBarycenter*/
    baryinput.delta = pulsarParams.dec + dtpos*pulsarParams.pmdec;
    baryinput.alpha = pulsarParams.ra +
      dtpos*pulsarParams.pmra/cos(baryinput.delta);
    
    if(j==0)
    {
      REAL8 dt;
      dt = times->data[j]-hetParams.epoch;
      fprintf(stderr, "f0 = %f,\nf1 = %e,\nf2 = %e,\nf3 = %e,\nf4 = %e.\n", 
      hetParams.f0, hetParams.f1, hetParams.f2, hetParams.f3, hetParams.f4);
      fprintf(stderr, "dt = %f, \n(1/120)*f4*dt^5 = %e.\n", dt, 
      (1.0/120.0)*hetParams.f4*dt*dt*dt*dt*dt);
    }

    XLALBarycenterEarthNew(&earth, &baryinput.tgps, edat, tdat, TIMECORRECTION_TDB);

    /* perform TN heterodyne */
    LALTimingNoiseHeterodyne( &status, &TNOutput, &TNInput, &hetParams,
      baryinput, earth );
    
    fprintf(fpout, "%lf\t%e\t%e\n", times->data[j], creal(TNOutput.Vh),
            cimag(TNOutput.Vh));

    fprintf(phifp, "%f\t%f\t%f\n", TNOutput.phi0, TNOutput.phi1, TNOutput.Dphase);
  }

  fclose(fpout);
  fclose(phifp);

  /* destroy vectors */
  XLALDestroyCOMPLEX16Vector( B );
  XLALDestroyREAL8Vector( times );
  XLALDestroyREAL8Vector(crabEphemerisData.f1);
  XLALDestroyREAL8Vector(crabEphemerisData.f0);
  XLALDestroyREAL8Vector(crabEphemerisData.tArr);
  XLALDestroyREAL8Vector(crabOutput.tArr);
  XLALDestroyREAL8Vector(crabOutput.f0);
  XLALDestroyREAL8Vector(crabOutput.f1);
  XLALDestroyREAL8Vector(crabOutput.f2);
  XLALDestroyREAL8Vector(crabOutput.f3);
  XLALDestroyREAL8Vector(crabOutput.f4);

  XLALDestroyEphemerisData(edat);
  XLALDestroyTimeCorrectionData(tdat);

  return 0;
}
