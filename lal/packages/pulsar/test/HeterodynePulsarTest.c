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

/******************************** <lalVerbatim file="HeterodynePulsarTestCV">
Author: Dupuis, R. J.
$Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{HeterodynePulsarTest.c}}

This program demonstrates the usage of the functions \texttt{LALCoarseHeterodyne()}
and \texttt{LALFineHeterodyneToPulsar()}.

\subsubsection*{Usage}
\begin{verbatim}
HeterodynePulsarTest
\end{verbatim}

\subsubsection*{Description}

This test program heterodynes, averages, and resamples an artificial signal using the functions
\texttt{LALCoarseHeterodyne()} and \texttt{LALFineHeterodyneToPulsar()}.

\subsubsection*{Exit codes}
\input{HeterodynePulsarTestCE}

\subsubsection*{Uses}
\begin{verbatim}
LALCoarseHeterodyne()
LALFineHeterodyneToPulsar()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{HeterodynePulsarTestCV}}
******************************************************* </lalLaTeX> */

/******* INCLUDE STANDARD LIBRARY HEADERS; ************/
/* note LALStdLib.h already includes stdio.h and stdarg.h */
#include <math.h>

/******* INCLUDE ANY LDAS LIBRARY HEADERS ************/

/******* INCLUDE ANY LAL HEADERS ************/
#include <lal/LALStdlib.h>
#include <lal/HeterodynePulsar.h>
#include <lal/AVFactories.h>
#include <lal/Random.h>
#include <lal/IIRFilter.h>
#include <lal/ZPGFilter.h>
#include <lal/LALInitBarycenter.h>

/******* DEFINE RCS ID STRING ************/

NRCSID( HETERODYNEPULSARTESTC, "$Id$" );

/******* DEFINE LOCAL CONSTANTS AND MACROS ************/
#define SQRT3_2 0.8660254037844386467637231707529361L

/***************************** <lalErrTable file="HeterodynePulsarTestCE"> */
#define HETERODYNEPULSARTESTC_ENOM 0
#define HETERODYNEPULSARTESTC_ECHK 1
#define HETERODYNEPULSARTESTC_EFLS 2

#define HETERODYNEPULSARTESTC_MSGENOM "Nominal exit"
#define HETERODYNEPULSARTESTC_MSGECHK "Error checking failed to catch bad data"
#define HETERODYNEPULSARTESTC_MSGEFLS "Incorrect answer for valid data"
/***************************** </lalErrTable> */

/* might also wish to define parameters and expected results for test
   cases here, for example */
#define HETERODYNEPULSARTEST_FH 0
#define HETERODYNEPULSARTEST_LENGTH 1048576
#define HETERODYNEPULSARTEST_DT 6.1035e-5  /* Jan 1, 2000, 00:00:00 */
#define HETERODYNEPULSARTEST_A0 100
#define HETERODYNEPULSARTEST_T0 630720013  /* Jan 1, 2000, 00:00:00 */
#define HETERODYNEPULSARTEST_BOXM 128
#define HETERODYNEPULSARTEST_IIRM 32
#define HETERODYNEPULSARTEST_FINEIIRM 16
#define HETERODYNEPULSARTEST_FC1 256
#define HETERODYNEPULSARTEST_FC2 4
#define HETERODYNEPULSARTEST_FINEFC 0.25
#define HETERODYNEPULSARTEST_F0 0    /* pulsar frequency */
#define HETERODYNEPULSARTEST_F1 0   /* pulsar spindown */
#define HETERODYNEPULSARTEST_F2 0  /* 2nd time derivative of pulsar frequency */
/******* DECLARE AND SET GLOBAL lalDebugLevel ************/

int lalDebugLevel = 0;
int main(void)
{
  UINT4 		i;
  REAL4 		dt = HETERODYNEPULSARTEST_DT;
  static LALStatus      status;
  static RandomParams 	*randomParams;
  static REAL4Vector  	*noise;
  INT4 			seed = 0;
  REAL4 		fh = HETERODYNEPULSARTEST_FH;
  CoarseHeterodyneOutput      coarseOutput;
  CoarseHeterodyneInput       coarseInput;
  CoarseHeterodyneParams      coarseParams;
  REAL4 		wc;
  COMPLEX8ZPGFilter *zpgFilter = NULL;
  REAL4IIRFilter *iirFilterRe = NULL; /* IIR filter */
  REAL4IIRFilter *iirFilterIm = NULL; /* IIR filter */
  UINT4 npoints = HETERODYNEPULSARTEST_LENGTH / HETERODYNEPULSARTEST_BOXM*HETERODYNEPULSARTEST_IIRM;
  REAL4 phase;
  FineHeterodyneInput fineInput;
  FineHeterodyneOutput fineOutput;
  FineHeterodyneParams fineParams;
  EphemerisData *edat = NULL;
  UINT4 itmp;
  char earth[] = "earth00.dat";
  char sun[] = "sun00.dat";


  /******* ALLOCATE MEMORY *************/

  coarseInput.V.data = NULL;
  LALCreateVector( &status, &coarseInput.V.data, HETERODYNEPULSARTEST_LENGTH );

  coarseOutput.Vh.data = NULL;
  LALCCreateVector( &status, &coarseOutput.Vh.data, npoints);

  fineInput.Vh.data = NULL;
  LALCCreateVector( &status, &fineInput.Vh.data, npoints);

  fineInput.varh.data = NULL;
  LALCCreateVector( &status, &fineInput.varh.data, npoints);

  fineOutput.B.data = NULL;
  LALZCreateVector( &status, &fineOutput.B.data, npoints /  HETERODYNEPULSARTEST_FINEIIRM);

  fineOutput.var.data = NULL;
  LALZCreateVector( &status, &fineOutput.var.data, npoints / HETERODYNEPULSARTEST_FINEIIRM);

  noise = NULL;
  LALCreateVector( &status, &noise, HETERODYNEPULSARTEST_LENGTH );
  LALCreateRandomParams( &status, &randomParams, seed);

  LALNormalDeviates( &status, noise, randomParams );

  if (fh == 0)
  {
    phase = 0;
  }
  else
  {
    phase = 2.0*LAL_PI*HETERODYNEPULSARTEST_T0/fh;
   }

   /******** GENERATE FAKE INPUT **********/
  for (i = 0;i < HETERODYNEPULSARTEST_LENGTH; i++)
    coarseInput.V.data->data[i] = 100*noise->data[i];

  coarseInput.V.data->length = HETERODYNEPULSARTEST_LENGTH;
  coarseInput.f0 =fh;

  coarseInput.V.epoch.gpsSeconds = HETERODYNEPULSARTEST_T0;
  coarseInput.V.epoch.gpsNanoSeconds = 0;

  coarseInput.V.deltaT = dt;

  /*******  TEST RESPONSE OF LALCoarseHeterodyne TO VALID DATA  ************/

  /* Create the first IIR filter. */

  wc = tan(LAL_PI * dt * HETERODYNEPULSARTEST_FC1);

  /* First create ZPG filter used to define IIR filter. */

  LALCreateCOMPLEX8ZPGFilter( &status, &zpgFilter, 0, 3 );
  zpgFilter->poles->data[0].re = wc*SQRT3_2;
  zpgFilter->poles->data[0].im = wc*0.5;
  zpgFilter->poles->data[1].re = 0.0;
  zpgFilter->poles->data[1].im = wc;
  zpgFilter->poles->data[2].re = -wc*SQRT3_2;
  zpgFilter->poles->data[2].im = wc*0.5;
  zpgFilter->gain.re = 0.0;
  zpgFilter->gain.im = wc*wc*wc;
  LALWToZCOMPLEX8ZPGFilter( &status, zpgFilter );
  LALCreateREAL4IIRFilter( &status, &iirFilterRe, zpgFilter );
  LALCreateREAL4IIRFilter( &status, &iirFilterIm, zpgFilter );
  LALDestroyCOMPLEX8ZPGFilter( &status, &zpgFilter );

  coarseParams.iirFilter1Re = iirFilterRe;
  coarseParams.iirFilter1Im = iirFilterIm;

  /* Create the second IIR filter */

  wc = tan(LAL_PI * dt* HETERODYNEPULSARTEST_BOXM * HETERODYNEPULSARTEST_FC2);

  LALCreateCOMPLEX8ZPGFilter( &status, &zpgFilter, 0, 3 );
  zpgFilter->poles->data[0].re = wc*SQRT3_2;
  zpgFilter->poles->data[0].im = wc*0.5;
  zpgFilter->poles->data[1].re = 0.0;
  zpgFilter->poles->data[1].im = wc;
  zpgFilter->poles->data[2].re = -wc*SQRT3_2;
  zpgFilter->poles->data[2].im = wc*0.5;
  zpgFilter->gain.re = 0.0;
  zpgFilter->gain.im = wc*wc*wc;
  LALWToZCOMPLEX8ZPGFilter( &status, zpgFilter );
  LALCreateREAL4IIRFilter( &status, &iirFilterRe, zpgFilter );
  LALCreateREAL4IIRFilter( &status, &iirFilterIm, zpgFilter );
  LALDestroyCOMPLEX8ZPGFilter( &status, &zpgFilter );

  coarseParams.iirFilter2Re = iirFilterRe;
  coarseParams.iirFilter2Im = iirFilterIm;

  coarseParams.boxM = HETERODYNEPULSARTEST_BOXM;
  coarseParams.iirM = HETERODYNEPULSARTEST_IIRM;

  coarseParams.stats = 2;

  LALCoarseHeterodyne( &status, &coarseOutput, &coarseInput, &coarseParams );

  if(status.statusCode)
  {
    printf("Unexpectedly got error code %d and message %s\n",
	    status.statusCode, status.statusDescription);
    return HETERODYNEPULSARTESTC_EFLS;
  }

    if(coarseOutput.Vh.data->length != HETERODYNEPULSARTEST_LENGTH / (HETERODYNEPULSARTEST_BOXM * HETERODYNEPULSARTEST_IIRM))
  {
    printf("Got incorrect length of output vector %d when expecting %d in LALCoarseHeterodyne\n",
	    coarseOutput.Vh.data->length, HETERODYNEPULSARTEST_LENGTH / (HETERODYNEPULSARTEST_BOXM * HETERODYNEPULSARTEST_IIRM));
    return  HETERODYNEPULSARTESTC_EFLS;
  }

  /*******  TEST RESPONSE OF LALCoarseHeterodyne TO INVALID DATA  ************/
/* Test that all the error conditions are correctly detected by the function */

#ifndef LAL_NDEBUG
if ( ! lalNoDebug ) {

 LALCoarseHeterodyne(&status, NULL, &coarseInput, &coarseParams);

  if (status.statusCode != HETERODYNEPULSARH_ENULLOUTPUT
       || strcmp(status.statusDescription, HETERODYNEPULSARH_MSGENULLOUTPUT))
  {
    printf( "Got error code %d and message %s\n",
	    status.statusCode, status.statusDescription);
    printf( "Expected error code %d and message %s\n",
	    HETERODYNEPULSARH_ENULLOUTPUT, HETERODYNEPULSARH_MSGENULLOUTPUT);
    return HETERODYNEPULSARTESTC_ECHK;
  }

 LALCoarseHeterodyne(&status, &coarseOutput, NULL, &coarseParams);

  if (status.statusCode != HETERODYNEPULSARH_ENULLINPUT
       || strcmp(status.statusDescription, HETERODYNEPULSARH_MSGENULLINPUT))
  {
    printf( "Got error code %d and message %s\n",
	    status.statusCode, status.statusDescription);
    printf( "Expected error code %d and message %s\n",
	    HETERODYNEPULSARH_ENULLINPUT, HETERODYNEPULSARH_MSGENULLINPUT);
    return HETERODYNEPULSARTESTC_ECHK;
  }

 LALCoarseHeterodyne(&status, &coarseOutput, &coarseInput, NULL);

  if (status.statusCode != HETERODYNEPULSARH_ENULLPARAMS
       || strcmp(status.statusDescription, HETERODYNEPULSARH_MSGENULLPARAMS))
  {
    printf( "Got error code %d and message %s\n",
	    status.statusCode, status.statusDescription);
    printf( "Expected error code %d and message %s\n",
	    HETERODYNEPULSARH_ENULLPARAMS, HETERODYNEPULSARH_MSGENULLPARAMS);
    return HETERODYNEPULSARTESTC_ECHK;
  }

 itmp = coarseParams.iirM;
 coarseParams.iirM -= 1;

 LALCoarseHeterodyne(&status, &coarseOutput, &coarseInput, &coarseParams);

  if (status.statusCode != HETERODYNEPULSARH_ERFACTOR
       || strcmp(status.statusDescription, HETERODYNEPULSARH_MSGERFACTOR))
  {
    printf( "Got error code %d and message %s\n",
	    status.statusCode, status.statusDescription);
    printf( "Expected error code %d and message %s\n",
	    HETERODYNEPULSARH_ERFACTOR, HETERODYNEPULSARH_MSGERFACTOR);
    return HETERODYNEPULSARTESTC_ECHK;
  }
  coarseParams.iirM = itmp;

  coarseInput.f0 = -1.0;

  LALCoarseHeterodyne(&status, &coarseOutput, &coarseInput, &coarseParams);

  if (status.statusCode != HETERODYNEPULSARH_EINVALIDF0
       || strcmp(status.statusDescription, HETERODYNEPULSARH_MSGEINVALIDF0))
  {
    printf( "Got error code %d and message %s\n",
	    status.statusCode, status.statusDescription);
    printf( "Expected error code %d and message %s\n",
	    HETERODYNEPULSARH_EINVALIDF0, HETERODYNEPULSARH_MSGEINVALIDF0);
    return HETERODYNEPULSARTESTC_ECHK;
  }

} /* if ( ! lalNoDebug ) */
#endif /* LAL_NDEBUG */

  /*******  TEST RESPONSE OF LALFineHeterodyneToPulsar TO VALID DATA  ************/

  /* set up detector */
  fineParams.detector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];

  /* set up pointer to ephemeris data */

  edat = (EphemerisData *)LALMalloc(sizeof(EphemerisData));

  (*edat).ephiles.earthEphemeris = earth;
  (*edat).ephiles.sunEphemeris = sun;

  LALInitBarycenter(&status, edat);

  fineParams.edat = edat;

  /******** GENERATE FAKE INPUT **********/
  for (i = 0;i < HETERODYNEPULSARTEST_LENGTH / (HETERODYNEPULSARTEST_BOXM*HETERODYNEPULSARTEST_IIRM); i++)
  {
    fineInput.Vh.data->data[i].re = 100. + 10.0*noise->data[i];
    fineInput.Vh.data->data[i].im = 100. + 10.0*noise->data[HETERODYNEPULSARTEST_LENGTH-i];
    fineInput.varh.data->data[i].re = 100. + 10.0*noise->data[HETERODYNEPULSARTEST_LENGTH-i];
    fineInput.varh.data->data[i].im = 100. + 10.0*noise->data[i];
  }

  fineInput.Vh.data->length = HETERODYNEPULSARTEST_LENGTH / (HETERODYNEPULSARTEST_BOXM*HETERODYNEPULSARTEST_IIRM);
  fineInput.varh.data->length = HETERODYNEPULSARTEST_LENGTH / (HETERODYNEPULSARTEST_BOXM*HETERODYNEPULSARTEST_IIRM);

  fineInput.f0 = HETERODYNEPULSARTEST_F0;
  fineInput.f1 = HETERODYNEPULSARTEST_F1;
  fineInput.f2 = HETERODYNEPULSARTEST_F2;
  fineInput.source.longitude = 0.0;
  fineInput.source.latitude = 0.0;
  fineInput.fEpochGPS = 230720013.0;
  fineInput.pmRA = 0.0;
  fineInput.pmDEC = 0.0;
  fineInput.posEpochGPS = 230720013.0;

  /* isolated pulsar - note that the function LALFineHeterodyneToPulsar doesn't accept binaries yet */
  fineInput.model = 0;

  fineParams.M = HETERODYNEPULSARTEST_FINEIIRM;

  fineParams.iirFlag = 0;

  /* Create the time-domain filter. */

  wc = tan(LAL_PI * dt * HETERODYNEPULSARTEST_BOXM *HETERODYNEPULSARTEST_IIRM* HETERODYNEPULSARTEST_FINEFC);

 /* First create ZPG filter used to define IIR filter. */

  LALCreateCOMPLEX8ZPGFilter( &status, &zpgFilter, 0, 3 );
  zpgFilter->poles->data[0].re = wc*SQRT3_2;
  zpgFilter->poles->data[0].im = wc*0.5;
  zpgFilter->poles->data[1].re = 0.0;
  zpgFilter->poles->data[1].im = wc;
  zpgFilter->poles->data[2].re = -wc*SQRT3_2;
  zpgFilter->poles->data[2].im = wc*0.5;
  zpgFilter->gain.re = 0.0;
  zpgFilter->gain.im = wc*wc*wc;
  LALWToZCOMPLEX8ZPGFilter( &status, zpgFilter );

  LALDestroyREAL4IIRFilter( &status, &iirFilterRe );
  LALDestroyREAL4IIRFilter( &status, &iirFilterIm );
  iirFilterRe = NULL;
  iirFilterIm = NULL;

  LALCreateREAL4IIRFilter( &status, &iirFilterRe, zpgFilter );
  LALCreateREAL4IIRFilter( &status, &iirFilterIm, zpgFilter );
  LALDestroyCOMPLEX8ZPGFilter( &status, &zpgFilter );

  fineParams.iirFilterRe = iirFilterRe;
  fineParams.iirFilterIm = iirFilterIm;



  LALFineHeterodyneToPulsar(&status, &fineOutput, &fineInput, &fineParams);

  if(status.statusCode)
  {
    printf("Unexpectedly got error code %d and message %s\n",
	    status.statusCode, status.statusDescription);
    return HETERODYNEPULSARTESTC_EFLS;
  }

    if(fineOutput.B.data->length != HETERODYNEPULSARTEST_LENGTH /
       (HETERODYNEPULSARTEST_BOXM*HETERODYNEPULSARTEST_IIRM*HETERODYNEPULSARTEST_FINEIIRM) )
  {
    printf("Got incorrect length of output vector %d when expecting %d\n",
	    fineOutput.B.data->length, HETERODYNEPULSARTEST_LENGTH /
       (HETERODYNEPULSARTEST_BOXM*HETERODYNEPULSARTEST_IIRM*HETERODYNEPULSARTEST_FINEIIRM));
    return  HETERODYNEPULSARTESTC_EFLS;
  }

 /*******  TEST RESPONSE OF LALFineHeterodyneToPulsar TO INVALID DATA  ************/
#ifndef LAL_NDEBUG
if ( ! lalNoDebug ) {

 LALFineHeterodyneToPulsar(&status, NULL, &fineInput, &fineParams);

  if (status.statusCode != HETERODYNEPULSARH_ENULLOUTPUT
       || strcmp(status.statusDescription, HETERODYNEPULSARH_MSGENULLOUTPUT))
  {
    printf( "Got error code %d and message %s\n",
	    status.statusCode, status.statusDescription);
    printf( "Expected error code %d and message %s\n",
	    HETERODYNEPULSARH_ENULLOUTPUT, HETERODYNEPULSARH_MSGENULLOUTPUT);
    return HETERODYNEPULSARTESTC_ECHK;
  }

 LALFineHeterodyneToPulsar(&status, &fineOutput, NULL, &fineParams);

  if (status.statusCode != HETERODYNEPULSARH_ENULLINPUT
       || strcmp(status.statusDescription, HETERODYNEPULSARH_MSGENULLINPUT))
  {
    printf( "Got error code %d and message %s\n",
	    status.statusCode, status.statusDescription);
    printf( "Expected error code %d and message %s\n",
	    HETERODYNEPULSARH_ENULLINPUT, HETERODYNEPULSARH_MSGENULLINPUT);
    return HETERODYNEPULSARTESTC_ECHK;
  }

 LALFineHeterodyneToPulsar(&status, &fineOutput, &fineInput, NULL);

  if (status.statusCode != HETERODYNEPULSARH_ENULLPARAMS
       || strcmp(status.statusDescription, HETERODYNEPULSARH_MSGENULLPARAMS))
  {
    printf( "Got error code %d and message %s\n",
	    status.statusCode, status.statusDescription);
    printf( "Expected error code %d and message %s\n",
	    HETERODYNEPULSARH_ENULLPARAMS, HETERODYNEPULSARH_MSGENULLPARAMS);
    return HETERODYNEPULSARTESTC_ECHK;
  }

 itmp = fineParams.M;
 fineParams.M -= 1;

 LALFineHeterodyneToPulsar(&status, &fineOutput, &fineInput, &fineParams);

  if (status.statusCode != HETERODYNEPULSARH_ERFACTOR
       || strcmp(status.statusDescription, HETERODYNEPULSARH_MSGERFACTOR))
  {
    printf( "Got error code %d and message %s\n",
	    status.statusCode, status.statusDescription);
    printf( "Expected error code %d and message %s\n",
	    HETERODYNEPULSARH_ERFACTOR, HETERODYNEPULSARH_MSGERFACTOR);
    return HETERODYNEPULSARTESTC_ECHK;
  }
 fineParams.M = itmp;

 fineInput.varh.data->length -= 1;

 LALFineHeterodyneToPulsar(&status, &fineOutput, &fineInput, &fineParams);

  if (status.statusCode != HETERODYNEPULSARH_ELENGTH
       || strcmp(status.statusDescription, HETERODYNEPULSARH_MSGELENGTH))
  {
    printf( "Got error code %d and message %s\n",
	    status.statusCode, status.statusDescription);
    printf( "Expected error code %d and message %s\n",
	    HETERODYNEPULSARH_ELENGTH, HETERODYNEPULSARH_MSGELENGTH);
    return HETERODYNEPULSARTESTC_ECHK;
  }

} /* if ( ! lalNoDebug ) */
#endif /* LAL_NDEBUG */

 /*******  CLEAN UP  ************/
 LALDestroyVector(&status, &coarseInput.V.data);
 LALCDestroyVector(&status, &coarseOutput.Vh.data);
 LALCDestroyVector(&status, &fineInput.Vh.data);
 LALCDestroyVector(&status, &fineInput.varh.data);

 LALZDestroyVector(&status, &fineOutput.B.data);
 LALZDestroyVector(&status, &fineOutput.var.data);
 LALDestroyVector(&status, &noise);
 LALDestroyRandomParams(&status, &randomParams);
 LALDestroyREAL4IIRFilter( &status, &iirFilterRe );
 LALDestroyREAL4IIRFilter( &status, &iirFilterIm );

 LALFree(edat->ephemE);
 LALFree(edat->ephemS);
 LALFree(edat);

 LALCheckMemoryLeaks();

 return HETERODYNEPULSARTESTC_ENOM;
}













