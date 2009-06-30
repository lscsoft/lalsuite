/*
*  Copyright (C) 2007 Jolien Creighton
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

/* <lalVerbatim file="fctTestCV">
 * Author: Edlund, Jeffrey A.
 * $Id$
 * </lalVerbatim>
 */

/* LALfctTest: This program shows the basic use of the LALfct package.
 *
 * In the future it will probably also have the standard error code tests that
 * the rest of the LAL packages have.
 *
 * Usage: LALfctTest accepts a number of commandline parameters.
 *  Valid Options are:
 *    -t   Activate timing (for performance testing)
 *    -r   read input data from input.dat
 *    -m   Locate the maximum
 *    -DX  dataSizeMax = X
 *    -dX  dataSizeMin = X
 *    -nX  numRuns = X
 *    -NX  numOfDims = X
 *    -zX  noise = X
 *
 */

#include <config.h>
#if defined LAL_FFTW3_ENABLED
/* fftw3 not yet supported */
int main( void ) { return 77; }
#else /* fftw2 implementation */

#include <lal/LALfct.h>  /* Include any required headers */
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <errno.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <sys/times.h>
#include <sys/time.h>
#include "timingfunctions.h"
#include <lal/AVFactories.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALConstants.h>


/* Define RCS ID string */
NRCSID(LALFCTTESTC,"$Id$");

/* Define local constants and macros */

/* I got these from a LAL package.  They are needed for the TestStatus function.*/
#define CODES_(x) #x
#define CODES(x) CODES_(x)

/* Declare and set global debuglevel */
/*   LAL 0.5 will change this variable be sure to change it then. */
int lalDebugLevel = 0; /*LALALLDBG; */

/* Declare local (static) functions (definitions can go here or at the
   end of the file) */
LALFCTREAL phaseFunction1(LALFCTREAL x);
LALFCTREAL phaseFunction2(LALFCTREAL x);
LALFCTREAL phaseFunction3(LALFCTREAL x);
LALFCTREAL mygasdev(void);
void generateData(LALFCTCOMPVector *inputDataVector, UINT2 numOfDims, \
		  LALFCTREAL noise);
void readData(LALFCTCOMPVector *data);
static void TestStatus( LALStatus *status, const char *ignored, int exitcode );
void fctTestCheckOutput( LALfctCalcOutput *fctCalcOutput);
void fctTestFindMax( LALfctCalcOutput *fctCalcOutput);

/* These are for getopt */
extern char *optarg;
extern int optind, opterr, optopt;

/* Define main function */
int main( int argc, char **argv ) {
  static LALStatus status; /* status initialization will complain if status is
			     not set to zero. */
  LALFCTCOMPVector *inputDataVector=0;
  LALfctInitParams fctInitParams;
  LALfctAddPhaseFuncParams fctAddPhaseFuncParams;
  LALfctCalcParams fctCalcParams;
  LALfctGenRowIndexParams fctGenRowIndexParams;
  LALfctCalcOutput fctCalcOutput;
  LALFCTPlan *fctPlan=0;

  UINT4 dataSize;
  UINT4 I;
  /* REMOVE  INT2 check; */
  UINT4 numRuns = 10;
  UINT2 numOfDims = 2;

  FILE *OUTPUT = 0;
  FILE *OUTPUTAVG = 0;

  REAL8 avgWallTime = 0;
  REAL8 wallDeviation = 0;
  REAL8 wallTimeTotal = 0;
  REAL8 wallTimeSqTotal = 0;
  REAL8 wallTempTime = 0;

  REAL8 avgUserTime = 0;
  REAL8 userDeviation = 0;
  REAL8 userTimeTotal = 0;
  REAL8 userTimeSqTotal = 0;
  REAL8 userTempTime = 0;

  LALFCTREAL noise = 0;

  BOOLEAN timeFlag = 0;
  BOOLEAN readFlag = 0;
  BOOLEAN maxFlag = 0;

  UINT4 dataSizeMin = pow(2,10);
  UINT4 dataSizeMax = pow(2,17);

  int index;
  int c;

  /* These have to be set to 0 or the CreateVector functions will complain */
  fctCalcOutput.rowIndex = 0;
  fctCalcOutput.outputData = 0;

  /* Check the commandline parameters */
  opterr = 0;
  while ((c = getopt (argc, argv, "trmD:d:n:N:z:")) != -1) {
    switch (c) {
      case 't':
	timeFlag = 1;
	printf("Timing\n");
	break;
      case 'r':
	readFlag = 1;
	printf("readflag = %d\n", readFlag);
	break;
      case 'm':
	maxFlag = 1;
	break;
      case 'D':
        dataSizeMax = atoi(optarg);
	printf("dataSizeMax = %d\n", dataSizeMax);
        break;
      case 'd':
        dataSizeMin = atoi(optarg);
	printf("dataSizeMin = %d\n", dataSizeMin);
        break;
      case 'n':
	numRuns = atoi(optarg);
	printf("numRuns = %d\n", numRuns);
	break;
      case 'N':
	numOfDims = atoi(optarg);
	printf("numOfDims = %d\n", numOfDims);
	break;
      case 'z':
	noise = atof(optarg);
	printf("noise = %f\n", noise);
	break;
      case '?':
	fprintf (stderr, "Unknown option `-%c'.\n", optopt);
	printf("Valid Options are:\n");
	printf(" -t   Activate timing (for performance testing)\n");
	printf(" -r   read input data from input.dat\n");
	printf(" -m   Locate the maximum\n");
	printf(" -DX  dataSizeMax = X\n");
	printf(" -dX  dataSizeMin = X\n");
	printf(" -nX  numRuns = X\n");
	printf(" -NX  numOfDims = X\n");
	printf(" -zX  noise = X\n");
	exit(1);
	break;
      default:
        exit(1);
    }
  }
  for (index = optind; index < argc; index++) {
    printf ("Non-option argument %s\n", argv[index]);
  }

  /* If we're not timing then only do one run at the minimum data size */
  if (!timeFlag) {
    dataSizeMax = dataSizeMin;
    numRuns = 1;
  }

  /* This loops over the dataSize */
  for ( dataSize = dataSizeMin; dataSize <= dataSizeMax; dataSize += 1024 ) {
    printf("dataSize = %d: ", dataSize);
    fflush(stdout);

    /* Initialize the fct_plan */
    fctInitParams.numOfDataPoints = dataSize;
    fctInitParams.numOfDims = numOfDims;
    LALfctInitialize( &status, &fctPlan, &fctInitParams );
    TestStatus( &status, CODES( 0 ), 1 );

    /* Add the proper phase functions */
    switch (numOfDims) {
      case 2:
	fctAddPhaseFuncParams.dim = 1;
	fctAddPhaseFuncParams.lengthOfDim = 128;
	fctAddPhaseFuncParams.phaseFuncForDim =  &phaseFunction1;
	LALfctAddPhaseFunc(&status, fctPlan, &fctAddPhaseFuncParams);
	TestStatus( &status, CODES( 0 ), 1 );
	break;
      case 3:
	fctAddPhaseFuncParams.dim = 1;
	fctAddPhaseFuncParams.lengthOfDim = 128;
	fctAddPhaseFuncParams.phaseFuncForDim =  &phaseFunction1;
	LALfctAddPhaseFunc(&status, fctPlan, &fctAddPhaseFuncParams);
	TestStatus( &status, CODES( 0 ), 1 );

	fctAddPhaseFuncParams.dim = 2;
	fctAddPhaseFuncParams.lengthOfDim = 64;
	fctAddPhaseFuncParams.phaseFuncForDim =  &phaseFunction2;
	LALfctAddPhaseFunc(&status, fctPlan, &fctAddPhaseFuncParams);
	TestStatus( &status, CODES( 0 ), 1 );
	break;
      case 4:
	fctAddPhaseFuncParams.dim = 1;
	fctAddPhaseFuncParams.lengthOfDim = 32;
	fctAddPhaseFuncParams.phaseFuncForDim =  &phaseFunction1;
	LALfctAddPhaseFunc(&status, fctPlan, &fctAddPhaseFuncParams);
	TestStatus( &status, CODES( 0 ), 1 );

	fctAddPhaseFuncParams.dim = 2;
	fctAddPhaseFuncParams.lengthOfDim = 32;
	fctAddPhaseFuncParams.phaseFuncForDim =  &phaseFunction2;
	LALfctAddPhaseFunc(&status, fctPlan, &fctAddPhaseFuncParams);
	TestStatus( &status, CODES( 0 ), 1 );

	fctAddPhaseFuncParams.dim = 3;
	fctAddPhaseFuncParams.lengthOfDim = 32;
	fctAddPhaseFuncParams.phaseFuncForDim =  &phaseFunction3;
	LALfctAddPhaseFunc(&status, fctPlan, &fctAddPhaseFuncParams);
	TestStatus( &status, CODES( 0 ), 1 );
	break;
      default:
	printf("Sorry this test program doesn't have phase functions");
	printf("for %d dimensions.\n", numOfDims);
	exit(1);
    }

    /* Create the proper size inputDataVector */
    LALFCTCOMPCreateVector ( &status, &inputDataVector, dataSize );
    TestStatus( &status, CODES( 0 ), 1 );

    /* Fill the inputDataVector with either data from a file or generated data. */
    if ((dataSize == 1024) && (readFlag)) {
      readData(inputDataVector);
    } else {
      generateData(inputDataVector, numOfDims, noise);
    }

    /* Fill in fctCalcParams */
    fctCalcParams.fctPlan = fctPlan;
    fctCalcParams.fctGenRowIndexParams = &fctGenRowIndexParams;
    fctCalcParams.fctGenRowIndexFunc = 0;

    fctGenRowIndexParams.fctPlan = fctPlan;
    fctGenRowIndexParams.goToEndOfRows = 1;
    fctGenRowIndexParams.createIndex = 1;
    fctGenRowIndexParams.skipRows = 0;
    fctGenRowIndexParams.numOfRows = 0; /* This is ignored */

    /* Reset the timing varibles */
    if (timeFlag) {
      wallTimeTotal = 0;
      wallTimeSqTotal = 0;
      userTimeTotal = 0;
      userTimeSqTotal = 0;
    }

    /* Loop around numRuns times to get an average run time. */
    for(I = 0; I < numRuns; I++) {

      /* Free the fctCalcOutput components if needed. */
      if (fctCalcOutput.outputData != 0) {
	LALFCTCOMPDestroyArray(&status, &(fctCalcOutput.outputData));
	fctCalcOutput.outputData = 0;
      }
      if (fctCalcOutput.rowIndex != 0) {
	LALU4DestroyArray(&status, &(fctCalcOutput.rowIndex));
	fctCalcOutput.rowIndex = 0;
      }

      /* Start timing */
      if (timeFlag) {
	WallTimeVal();
	UserTimeVal();
      }

      /* Run LALfctCalc */
      LALfctCalc(&status, &fctCalcOutput, inputDataVector, &fctCalcParams);

      /* Read the time from the timers */
      if (timeFlag) {
	userTempTime = UserTimeVal();
	wallTempTime = WallTimeVal();
      }

      /* Test the status that was returned from LALfctCalc*/

      TestStatus( &status, CODES( 0 ), 1 );

      /* Add the time from this run to the totals */
      if (timeFlag) {
	userTimeTotal += userTempTime;
	userTimeSqTotal += userTempTime * userTempTime;
	wallTimeTotal += wallTempTime;
	wallTimeSqTotal += wallTempTime * wallTempTime;
      }
    } /* for numRuns */

    if (timeFlag) {

      /* Calculate the Averages and the Deviations */
      avgWallTime = wallTimeTotal / (double) numRuns;
      wallDeviation = ((wallTimeSqTotal / (double)  numRuns) - \
		       avgWallTime * avgWallTime) *(numRuns/(numRuns-1));
      avgUserTime = userTimeTotal / (double) numRuns;
      userDeviation = ((userTimeSqTotal / (double)  numRuns) - \
		       avgUserTime * avgUserTime) * (numRuns/(numRuns-1));

      /* Print them out to the screen. */
      printf("AvgWallTime = %g, WallDeviation = %g,", \
	     avgWallTime, wallDeviation);
      printf("AvgUserTime = %g, UserDeviation = %g\n", \
	     avgUserTime, userDeviation);

      /* Print them out to a file. */
      OUTPUTAVG = fopen ("./2dtimesAvg.dat", "a");
      fprintf(OUTPUTAVG, "%d, %d, %g, %g, %g, %g\n", \
	      dataSize, numRuns, avgWallTime, wallDeviation, \
	      avgUserTime, userDeviation);
      fclose(OUTPUTAVG);
    }

    if ((!timeFlag) || (maxFlag)) {

      /* Move the Cursor to the next line. */
      printf("\n");

      /* Locate the Maximum in the Output. */
      fctTestFindMax( &fctCalcOutput );

      /* Check the Output to make sure that the power in all the rows
	 is the same. */
      fctTestCheckOutput( &fctCalcOutput );

    }

    if (!timeFlag) {

      /* Output the data to a file. */
      if (numOfDims == 2) {
	OUTPUT = fopen ("./output.dat", "w");
	fwrite(fctCalcOutput.outputData->data, sizeof(REAL4), 2 * \
	       fctCalcOutput.outputData->dimLength->data[0] * \
	       fctCalcOutput.outputData->dimLength->data[1], OUTPUT);
	fclose(OUTPUT);
      }

    }

    /* free the components of fctCalcOutput */
    LALFCTCOMPDestroyArray(&status, &(fctCalcOutput.outputData));
    fctCalcOutput.outputData = 0;
    LALU4DestroyArray(&status, &(fctCalcOutput.rowIndex));
    fctCalcOutput.rowIndex = 0;

    /* free the inputData Vector */
    LALFCTCOMPDestroyVector(&status, &inputDataVector);

    /*    printf("numOfDim = %d\n", fctPlan.numOfDims);*/
    /* Destroy the fctPlan */
    LALfctDestroyPlan(&status, &fctPlan);
    TestStatus( &status, CODES( 0 ), 1 );

    /*Check for Memory leaks.  This is only effective if you set debuglevel
      to something like LALALLDBG. (some value that causes memory checks at
      least.) */
    LALCheckMemoryLeaks();

  } /* for dataSize */

  return 0;
}

LALFCTREAL phaseFunction1(LALFCTREAL x) {
  return(-x*x);
  /*return(pow(x, -5.0/3.0));*/
}

LALFCTREAL phaseFunction2(LALFCTREAL x) {
  return(x*x*x);
  /*return(1/x);*/
}

LALFCTREAL phaseFunction3(LALFCTREAL x) {
  return(x*x*x*x);
  /*return(pow(x, -2.0/3.0));*/
}


void generateData(LALFCTCOMPVector *inputDataVector, UINT2 numOfDims, \
		  LALFCTREAL noise) {
  /* This function fills the inputDataVector with a generated chirp.
     using the noise parameter, you can add gaussian noise. */
  LALFCTREAL x;
  LALFCTREAL phase = 0;
  UINT4 i;

  for( i = 0; i < inputDataVector->length; ++i ){
    x = i / ( (double) inputDataVector->length );
    switch ( numOfDims ) {
    case 2:
      phase = 2*LAL_PI * x * 30. + 2*LAL_PI*20*phaseFunction1(x);
      break;
    case 3:
      phase = 2*LAL_PI * x * 30. + 2*LAL_PI*20*phaseFunction1(x) + \
	2*LAL_PI*11*phaseFunction2(x);
      break;
    case 4:
      phase = 2*LAL_PI * x * 30. + 2*LAL_PI*20*phaseFunction1(x) + \
	2*LAL_PI*11*phaseFunction2(x) + 2*LAL_PI*5*phaseFunction3(x);
      break;
    default:
      printf("Sorry this test program doesn't have phase functions");
      printf("for %d dimensions.\n", numOfDims);
      exit(1);
      break;
    }
    inputDataVector->data[i].re = cos(phase) + noise*mygasdev();
    inputDataVector->data[i].im = sin(phase) + noise*mygasdev();
  }
}

void readData(LALFCTCOMPVector *inputDataVector) {
  /* This function fills the inputDataVector with data from input.dat */
  FILE *INPUT;
  int Check;

  INPUT = fopen("./input.dat", "r");
  if (INPUT == NULL) {
    printf("Error couldn't read Data File.\n");
    exit(1);
  }
  Check = fread(inputDataVector->data, sizeof(LALFCTREAL), \
		2 * inputDataVector->length, INPUT);
  if ( Check != (int)(2 * inputDataVector->length) ) {
    printf("Error couldn't read Data File.  Read Check = %d\n", Check);
    exit(1);
  }
}

LALFCTREAL mygasdev(void) {
  /* This function generates gaussian noise. */
  static INT2 iset = 0;
  static LALFCTREAL gset;
  LALFCTREAL fac, r, v1, v2;
  srand(time(0));
  /*float ran1(long*);*/
  /*printf("Entered mygasdev\n");*/
  if (iset==0) {
    do {
      v1 = 2.0 * ((LALFCTREAL) rand())/(RAND_MAX+1.0) - 1.0;
      v2 = 2.0 * ((LALFCTREAL) rand())/(RAND_MAX+1.0) - 1.0;
      r = v1 * v1 + v2 * v2;
    } while (r >= 1.0);

    fac = sqrt(-2.0 * log(r)/r);
    gset = v1 * fac;
    iset = 1;
    return (v2 * fac);
  } else {
    iset = 0;
    return gset;
  }
}



void fctTestCheckOutput( LALfctCalcOutput *fctCalcOutput){
  /* This function calculates the power in each row and keeps track of the
     minimum and maximum power.  It then prints out the minimum, maximum, and
     the difference.  It's just for debugging and sanity checks. */

  UINT4 row;
  UINT4 location;
  UINT4 elementCursor;
  LALFCTREAL power;
  LALFCTREAL powermax = 0;
  LALFCTREAL powermin = 0;

  for (row = 0; row < fctCalcOutput->outputData->dimLength->data[0]; row++) {
    location =  row * (fctCalcOutput->outputData->dimLength->data[1]);
    power = 0;
    for (elementCursor = 0; \
	   elementCursor < fctCalcOutput->outputData->dimLength->data[1]; \
	   elementCursor++) {
      power += fctCalcOutput->outputData->data[location+elementCursor].im * \
	fctCalcOutput->outputData->data[location+elementCursor].im + \
	fctCalcOutput->outputData->data[location+elementCursor].re * \
	fctCalcOutput->outputData->data[location+elementCursor].re;
    }
    /* printf("row = %d, power = %f\n", row, power);*/
    if (power > powermax) {
      powermax = power;
    }
    if (power < powermin) {
      powermin = power;
    } else if (powermin < 1) {
      powermin = power;
    }
  }
  printf("powermax = %f, powermin = %f, Difference = %g\n", powermax, powermin, \
	 powermax-powermin);
}


void fctTestFindMax( LALfctCalcOutput *fctCalcOutput){
  /* This function finds the maximum in the output and prints it's location. */
  UINT2 J;
  UINT4 row;
  UINT4 location;
  UINT4 elementCursor;
  UINT4 maxRow = 0;
  UINT4 maxRowElement = 0;

  LALFCTREAL temp;
  LALFCTREAL max = 0;

  for (row = 0; row < fctCalcOutput->outputData->dimLength->data[0]; row++) {
    location =  row * (fctCalcOutput->outputData->dimLength->data[1]);
    for (elementCursor = 0; \
	   elementCursor < fctCalcOutput->outputData->dimLength->data[1]; \
	   elementCursor++) {
      temp = fctCalcOutput->outputData->data[location+elementCursor].im * \
	fctCalcOutput->outputData->data[location+elementCursor].im + \
	fctCalcOutput->outputData->data[location+elementCursor].re * \
	fctCalcOutput->outputData->data[location+elementCursor].re;
      if (temp > max) {
	max = temp;
	maxRow = row;
	maxRowElement = elementCursor;
      }
    }
  }
  printf("maxRow = %d, maxRowElement = %d: ", maxRow, maxRowElement);

  location =  maxRow * fctCalcOutput->rowIndex->dimLength->data[1];
  printf("(%d", maxRowElement);
  for (J = 0; J < fctCalcOutput->rowIndex->dimLength->data[1]; J++) {
    printf(", %d", fctCalcOutput->rowIndex->data[location+J]);
  }
  printf("): max = %f\n", max);

}


static void TestStatus( LALStatus *status, const char *ignored, int exitcode ) {
  /* This function tests the status structure. I took most of this from
     one of the LAL fft test programs. */
  char  str[64];
  char *tok;

  if ( strncpy( str, ignored, sizeof( str ) ) ) {
    if ( (tok = strtok( str, " " ) ) ) {
      do {
	if ( status->statusCode == atoi( tok ) ) {
	  return;
	}
      }
      while ( ( tok = strtok( NULL, " " ) ) );
    } else {
      if ( status->statusCode == atoi( tok ) ) {
	return;
      }
    }
  }

  REPORTSTATUS( status );
  fprintf( stderr, "\nExiting to system with code %d\n", exitcode );
  exit( exitcode );
}

#endif /* fftw2 implementation */
