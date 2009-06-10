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

/* <lalVerbatim file="LALfctSampleCV">
 * Author: Edlund, Jeffrey A.
 * $Id$
 * </lalVerbatim>
 */

/* <lalLaTeX>
   \section{Sample Program \texttt{LALfctSample.c}}
   \label{ss:LALfctSample.c}

   This program shows the basic use of the LALfct package.
 It also allows the user to specify the amount of memory that the output
 array should use. It then calculates the output of the FCT in chunks that
 will fit into memory and appends each chunk to the output file.

 This document goes step by step through the code of LALfctSample.c as an
 example of how to use the LALfct library. It may seem a little redundant,
 but hopefully it will make it easier for someone to use the library.
 </lalLaTeX> */

/* <lalLaTeX>
   \subsubsection*{Command line Usage}
   </lalLaTeX> <lalVerbatim> */
 /* Usage: LALfctTest accepts these commandline parameters.
   Valid Options are:
      -d X  dataSize = X  (Where X is an integer)
      -m F  maxSizeOfOutputArrayInMegs = F (Where F is a float)
            Use this to deal with memory constraints.
      -i inputFileName
      -o outputFileName

  An index of the Rows output will also be dumped out to inputFileName.index

  Example Commandline:
  LALfctSample -d 1024 -m 128 -i input.dat -o output.dat */
/* </lalVerbatim> */

/* <lalLaTeX>
   \subsubsection*{Required Headers}
   The following headers are probably the minimum required to use LALfct.
   But your needs may vary.
   </lalLaTeX> */

#include <config.h>
#if defined LAL_FFTW3_ENABLED
/* fftw3 not yet supported */
int main( void ) { return 77; }
#else /* fftw3 implementation */

/*<lalVerbatim> */
#include <lal/LALfct.h>  /* Include any required headers */
#include <stdio.h>
#include <lal/AVFactories.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALfct.h>

#include <unistd.h> /* This is for getopt */
  /* </lalVerbatim> */

/* <lalLaTeX>
   \subsubsection*{Other definitions and declarations}
   LAL standards require the following command to define the RCS ID string.
   </lalLaTeX>
   <lalVerbatim> */
NRCSID(LALFCTSAMPLEC,"$Id$");
/* </lalVerbatim> */


/* <lalLaTeX>
   I got these CODES defines from a LAL package.  They are needed for the TestStatus
   function. Hopefully they'll have a standard way of checking the return status of
   functions soon.
   </lalLaTeX> <lalVerbatim>*/
#define CODES_(x) #x
#define CODES(x) CODES_(x)
/* </lalVerbatim> */

/* <lalLaTeX>
   Declare and set global laldebuglevel. Set to 0 for no debugging info.  Set
   to LALALLDBG for full debugging.
   </lalLaTeX> <lalVerbatim>*/
int lalDebugLevel = 0;
/* </lalVerbatim> */

/* <lalLaTeX>
   Local (static) function declarations (definitions are after the main function.):
   </lalLaTeX> <lalVerbatim>*/
LALFCTREAL phaseFunction1(LALFCTREAL x);
void readData(LALFCTCOMPVector *inputDataVector, char *inputFileName);
static void TestStatus( LALStatus *status, const char *ignored, int exitcode );
static INT4 CheckStatus( LALStatus *status);
/* </lalVerbatim> */

/*  <lalLaTeX>
    These are for getopt
    </lalLaTeX> <lalVerbatim> */
extern char *optarg;
extern int optind, opterr, optopt;
int indx;
int c;
/* </lalVerbatim> */

/* <lalLaTeX>
   \subsubsection*{Main Function}
   </lalLaTeX> <lalVerbatim>  */
int main( int argc, char **argv ) {
  /* </lalVerbatim> */

  /* <lalLaTeX>
     Declare and zero the LALStatus structure that will be passed to the
     the LAL functions.
     </lalLaTeX> <lalVerbatim> */
  static LALStatus status; /* status initialization will complain if status is
			     not set to zero. */
  /* </lalVerbatim> <lalLaTeX>
     Declare the LALfct Structures that will be passed to the LALfct functions.
     </lalLaTeX> <lalVerbatim> */
  LALFCTCOMPVector *inputDataVector=0;
  LALfctInitParams fctInitParams;
  LALfctAddPhaseFuncParams fctAddPhaseFuncParams;
  LALfctCalcParams fctCalcParams;
  LALfctGenRowIndexParams fctGenRowIndexParams;
  LALfctGenRowIndexOutput fctGenRowIndexOutput;
  LALfctCalcOutput fctCalcOutput;
  LALFCTPlan *fctPlan=0;
  /* </lalVerbatim> */

  /* <lalLaTeX>
     The following variables will be used to divide the FCT calculation
     into reasonably sized chunks.
     </lalLaTeX>
     <lalVerbatim> */
  LALFCTREAL maxSizeOfOutputArrayInMegs = 0;
  LALFCTREAL megsPerRow = 0;
  UINT4 numOfIndices = 0;
  UINT4 numOfRowsPerIndex = 0;
  /* </lalVerbatim>
     <lalLaTeX>
     Define input and output variables:
     </lalLaTeX> <lalVerbatim> */
  char inputFileNameDefault[] = "./input.dat";
  char outputFileNameDefault[] = "./output.dat";
  char outputIndexFileNameDefault[] = "./output.dat.index";
  char *inputFileName = inputFileNameDefault;
  char *outputFileName = outputFileNameDefault;
  char *outputIndexFileName = outputIndexFileNameDefault;
  FILE *OUTPUT = NULL;
  BOOLEAN outputIndexFileNameMalloced = 0;
  UINT2 outputIndexFileNameSize = 0;
  /* </lalVerbatim> <lalLaTeX>
     Define variables that will determine the length of the data
     that will be processed and the number of phase functions that
     will be used.
     </lalLaTeX> <lalVerbatim> */
  UINT4 dataSize=0;
  UINT2 numOfDims = 2;
  /* </lalVerbatim> */

  /* <lalLaTeX>
     These have to be set to 0 or the CreateVector functions will complain.
     </lalLaTeX> <lalVerbatim>*/
  fctCalcOutput.rowIndex = 0;
  fctCalcOutput.outputData = 0;
  /* </lalVerbatim> */

  /* <lalLaTeX>
     I'm going to skip over the code that calls getopt to check the
     commandline parameters.  Take a look at the code if you're
     interested.
     </lalLaTeX>*/
  /* Check the commandline parameters */
  printf("Reading commandline options:\n");
  opterr = 0;
  while ((c = getopt (argc, argv, "d:m:i:o:")) != -1) {
    switch (c) {
      case 'd':
        dataSize = atoi(optarg);
	printf("  dataSize = %d\n", dataSize);
        break;
      case 'm':
	maxSizeOfOutputArrayInMegs = atof(optarg);
	printf("  maxSizeOfOutputArrayInMegs = %g\n", \
	       maxSizeOfOutputArrayInMegs);
	break;
      case 'i':
	inputFileName = optarg;
	printf("  inputFileName = %s\n", inputFileName);
	break;
      case 'o':
	outputFileName = optarg;
	printf("  outputFileName = %s\n", outputFileName);
	if (outputIndexFileNameMalloced) {
	  LALFree(outputIndexFileName);
	}
	outputIndexFileNameSize = strlen(outputFileName) + strlen(".index");
	outputIndexFileName = LALCalloc(outputIndexFileNameSize + 1, \
					sizeof(char));
	strncpy(outputIndexFileName, outputFileName, outputIndexFileNameSize);
	strncat(outputIndexFileName, ".index", 6);

	break;
      case '?':
	fprintf (stderr, "Unknown option `-%c'.\n", optopt);
	printf("Usage:\n");
	printf(" -d X  dataSize = X  (Where X is an integer)\n");
	printf(" -m F  maxSizeOfOutputArrayInMegs = F (Where F is a float)\n");
	printf(" -i inputFileName\n");
	printf(" -o outputFileName\n");

	exit(1);
	break;
      default:
        exit(1);
    }
  }
  for (indx = optind; indx < argc; indx++) {
    printf ("Non-option argument %s\n", argv[indx]);
  }

  printf("Done reading commandline options.\n\n");

  /* <lalLaTeX>
     dataSize must be specified on the command line.
     </lalLaTeX> <lalVerbatim> */
  if (dataSize == 0) {
    printf("Please specify the dataSize on the Commandline");
    printf("using the -dX option (where X is an integer.)\n");
    exit(1);
  }
  /* </lalVerbatim> */
  printf("\n");


  /* <lalLaTeX>
     Initialize the fctPlan.
  </lalLaTeX> <lalVerbatim> */
  fctInitParams.numOfDataPoints = dataSize; /* Number of points in input data. */
  fctInitParams.numOfDims = numOfDims; /* Number of phase functions that will be used. */
  LALfctInitialize( &status, &fctPlan, &fctInitParams );
  TestStatus( &status, CODES( 0 ), 1 );
  /* </lalVerbatim> */

  /*  <lalLaTeX>
      Add the phase function to the fctPlan.
      </lalLaTeX> <lalVerbatim> */
  fctAddPhaseFuncParams.dim = 1;  /* The dimension in the N dimensional
				     space that this phase function effects. */
  fctAddPhaseFuncParams.lengthOfDim = 128; /* This determines the range
					      that you want to use to
					      resolve the dimension. */
  fctAddPhaseFuncParams.phaseFuncForDim =  &phaseFunction1; /* A reference to
							       the phase
							       function. */
  LALfctAddPhaseFunc(&status, fctPlan, &fctAddPhaseFuncParams);
  TestStatus( &status, CODES( 0 ), 1 );
  /* </lalVerbatim> */

  /* <lalLaTeX>
     Create the proper size inputDataVector
     </lalLaTeX> <lalVerbatim> */
  LALFCTCOMPCreateVector ( &status, &inputDataVector, dataSize );
  TestStatus( &status, CODES( 0 ), 1 );
  /* </lalVerbatim> */

  /* <lalLaTeX>
     Fill the inputDataVector with data from an input file.
  </lalLaTeX> <lalVerbatim>*/
  readData(inputDataVector, inputFileName);
  /* </lalVerbatim> */

  /* <lalLaTeX>
     Fill in fctCalcParams
     </lalLaTeX> <lalVerbatim>*/
  fctCalcParams.fctPlan = fctPlan;
  fctCalcParams.fctGenRowIndexParams = &fctGenRowIndexParams;
  fctCalcParams.fctGenRowIndexFunc = 0; /* Setting this to zero tells
					   LALfctCalc to use the internal
					   LALfctGenRowIndex */
  /* </lalVerbatim> */


  /* <lalLaTeX>
     With some applications of the FCT, the entire output will not fit
     into the memory of the machine you're running on.  The following
     section of code will divide the problem into chunks and process
     each chunk of the problem individually.
     </lalLaTeX> */

  /* <lalLaTeX>
     Initialize the fctGenRowIndexParams structure.
     </lalLaTeX> <lalVerbatim> */
  fctGenRowIndexParams.fctPlan = fctPlan;
  fctGenRowIndexParams.goToEndOfRows = 1; /* Keep going to the end */
  fctGenRowIndexParams.createIndex = 0; /* Don't create the index yet. */
  fctGenRowIndexParams.skipRows = 0;
  fctGenRowIndexParams.numOfRows = 0; /* This is ignored  because
					 goToEndOfRows is set.*/
  /* </lalVerbatim>  */

  /* <lalLaTeX> Let's figure of the number of chunks that we need to
     calculate. </lalLaTeX> <lalVerbatim> */
  if (maxSizeOfOutputArrayInMegs != 0) {
    megsPerRow = ((LALFCTREAL) dataSize*sizeof(LALFCTCOMP)) / 1048576.0;
    /* </lalVerbatim> <lalLaTeX> Let's ask LalfctGenRowIndex to return the
       number of rows that we will be calculating. This is returned in
       fctGenRowIndexOutput.numOfRows.  Note that no index is created because
       fctGenRowIndexParams.createIndex is set to zero.
       </lalLaTeX> <lalVerbatim> */
    LALfctGenRowIndex(&status, &fctGenRowIndexOutput, &fctGenRowIndexParams);
    TestStatus( &status, CODES( 0 ), 1 );
    /* </lalVerbatim> <lalLaTeX>
       Using fctGenRowIndexOutput.numOfRows, lets calculate the number of
       Indices (or chunks) needed.
    </lalLaTeX> <lalVerbatim> */
    numOfRowsPerIndex = maxSizeOfOutputArrayInMegs / megsPerRow;
    numOfIndices =  fctGenRowIndexOutput.numOfRows / numOfRowsPerIndex;
    if (numOfIndices * numOfRowsPerIndex < fctGenRowIndexOutput.numOfRows) {
      numOfIndices++;
    }
    /* </lalVerbatim> <lalLaTeX>
       If maxSizeOfOutputArrayInMegs was not specified on the command line
       then just do the whole thing in one chunk.
       </lalLaTeX> <lalVerbatim> */
  } else {
    numOfIndices = 1;
    numOfRowsPerIndex = 0; /* This will be ignored, but I set it anyway. */
  }
  /* </lalVerbatim> <lalLaTeX>
     Okay, now lets Calculate the chunks.
     </lalLaTeX> <lalVerbatim> */
  printf("Calculating the FCT:\n");
  fctGenRowIndexParams.createIndex = 1; /* Create a row index this time. */
  fctGenRowIndexParams.numOfRows = numOfRowsPerIndex; /* Calculate only
							 numOfRowsPerIndex
							 number of rows each time. */
  /* loop through the chunks */
  for (indx = 1; indx <= (int)numOfIndices; indx++) {

    if (indx == (int)numOfIndices) {
      fctGenRowIndexParams.goToEndOfRows = 1; /* If it's the last one then
						 ignore the numOfRows and
						 just go to the end of the
						 space. */
    } else {
      fctGenRowIndexParams.goToEndOfRows = 0;
    }

    /* Skip the rows that we've already calculated. */
    fctGenRowIndexParams.skipRows = numOfRowsPerIndex * (indx - 1);

    /* Run LALfctCalc on the chunk */
    LALfctCalc(&status, &fctCalcOutput, inputDataVector, &fctCalcParams);
    TestStatus( &status, CODES( 0 ), 1 );

    /* Store the output somewhere */
    if (indx == 1) {
      OUTPUT = fopen (outputFileName, "w");
    } else {
      OUTPUT = fopen (outputFileName, "a");
    }
    fwrite(fctCalcOutput.outputData->data, sizeof(REAL4), 2 * \
	   fctCalcOutput.outputData->dimLength->data[0] * \
	   fctCalcOutput.outputData->dimLength->data[1], OUTPUT);
    fclose(OUTPUT);

    /* Store the index to the output somewhere. */
    if (indx == 1) {
      OUTPUT = fopen (outputIndexFileName, "w");
    } else {
      OUTPUT = fopen (outputIndexFileName, "a");
    }
    fwrite(fctCalcOutput.rowIndex->data, sizeof(UINT2),
	   fctCalcOutput.rowIndex->dimLength->data[0] * \
	   fctCalcOutput.rowIndex->dimLength->data[1], OUTPUT);
    fclose(OUTPUT);


    printf("  Calculated chunk %d of %d.\n", indx, numOfIndices);

    /* Destroy the output Data Array */
    LALFCTCOMPDestroyArray(&status, &(fctCalcOutput.outputData));
    if (CheckStatus(&status)) {
      LALU4DestroyArray(&status, &(fctCalcOutput.rowIndex));
      TestStatus( &status, CODES( 0 ), 1 );
    }
    fctCalcOutput.outputData = 0; /* make sure that it's zero. */


    /* Destroy the output index Array */
    LALU4DestroyArray(&status, &(fctCalcOutput.rowIndex));
    if (CheckStatus(&status)) {
      TestStatus( &status, CODES( 0 ), 1 );
    }
    fctCalcOutput.rowIndex = 0; /* Make sure that it's zero. */

  }
  /*</lalVerbatim> <lalLaTeX>
    We're done calculating the FCT.  So Let's clean things up and exit.
    </lalLaTeX> <lalVerbatim> */

  /* free the inputData Vector */
  LALFCTCOMPDestroyVector(&status, &inputDataVector);

  /* Destroy the fctPlan */
  LALfctDestroyPlan(&status, &fctPlan);
  TestStatus( &status, CODES( 0 ), 1 );

  if (outputIndexFileNameMalloced) {
    LALFree(outputIndexFileName);
  }

  /*Check for Memory leaks.  This is only effective if you set debuglevel
    to something like LALALLDBG. (or some other value that causes memory
    checks at least.) */
  LALCheckMemoryLeaks();

  return 0;
}
/*</lalVerbatim> */

/*<lalLaTeX>
  \subsubsection*{Miscellaneous Functions}

  This is a sample phase function.
  </lalLaTeX> <lalVerbatim> */
LALFCTREAL phaseFunction1(LALFCTREAL x) {
  return(x*x);
}
/*</lalVerbatim> <lalLaTeX>

This is the function that reads the data from the input data file
</lalLaTeX> <lalVerbatim> */
void readData(LALFCTCOMPVector *inputDataVector, char *inputFileName) {
  /* This function fills the inputDataVector with data from input.dat */
  FILE *INPUT;
  int Check;

  INPUT = fopen( inputFileName, "r" );
  if (INPUT == NULL) {
    printf("Error couldn't read Data File. inputFileName = %s\n", \
	   inputFileName);
    exit(1);
  }
  Check = fread(inputDataVector->data, sizeof(LALFCTREAL), \
		2 * inputDataVector->length, INPUT);
  if ( Check != (int)(2 * inputDataVector->length) ) {
    printf("Error couldn't read Data File.  Read Check = %d\n", Check);
    exit(1);
  }
}
/* </lalVerbatim> <lalLaTeX>
   The following are two functions that I found in another LAL package
   to test the status returned from LAL functions.
   </lalLaTeX> <lalVerbatim> */

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

static INT4 CheckStatus( LALStatus *status) {
  /* This function reports nonzero values of status. */
  if ( status->statusCode != 0) {
    REPORTSTATUS( status );
    return(1);
  }
  return(0);
}
/*</lalVerbatim>*/
#endif /* fftw2 implementation */
