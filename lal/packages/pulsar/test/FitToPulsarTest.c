/******************************** <lalVerbatim file="FitToPulsarTestCV">
Author: Dupuis, R.J.
$Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{FitToPulsarTest.c}}

This test program demonstrates the correct usage of the functions
\texttt{LALCoarseFitToPulsar} and \texttt{LALFineFitToPulsar}.

\subsubsection*{Usage}
\begin{verbatim}
FitToPulsarTest
\end{verbatim}

\subsubsection*{Description}
To be completed.
\subsubsection*{Exit codes}
\input{FitToPulsarTestCE}

\subsubsection*{Uses}
\begin{verbatim}
LALCoarseFitToPulsar()
LALFineFitToPulsar()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{FitToPulsarCV}}
******************************************************* </lalLaTeX> */

/******* INCLUDE STANDARD LIBRARY HEADERS; ************/
/* note LALStdLib.h already includes stdio.h and stdarg.h */

/******* INCLUDE ANY LDAS LIBRARY HEADERS ************/

/******* INCLUDE ANY LAL HEADERS ************/
#include <lal/FitToPulsar.h>
#include <lal/Random.h>

/******* DEFINE RCS ID STRING ************/

NRCSID( FitToPulsarTestC, "$Id$" );

/******* DEFINE LOCAL CONSTANTS AND MACROS ************/

/***************************** <lalErrTable file="FitToPulsarTestCE"> */
#define FITTOPULSARTESTC_ENOM 0
#define FITTOPULSARTESTC_ECHK 1
#define FITTOPULSARTESTC_EFLS 2

#define FITTOPULSARTESTC_MSGENOM "Nominal exit"
#define FITTOPULSARTESTC_MSGECHK "Error checking failed to catch bad data"
#define FITTOPULSARTESTC_MSGEFLS "Incorrect answer for valid data"
/***************************** </lalErrTable> */

/* might also wish to define parameters and expected results for test
   cases here, for example */

#define FITTOPULSARTEST_LENGTH 24*10
#define FITTOPULSARTEST_T0 630720013  /* Jan 1, 2000, 00:00:00 */

/* int lalDebugLevel = LALMSGLVL1; */
int lalDebugLevel = 0;
int main(void)
{
  static LALStatus     	status;
  LALDetector 		detector;
  LALSource 		pulsar;
  CoarseFitOutput  	output;
  CoarseFitInput       	input;
  CoarseFitParams      	params;
  LIGOTimeGPS 		tgps[FITTOPULSARTEST_LENGTH];
  REAL4			cosIota;
  REAL4			phase;
  REAL4 		psi;
  REAL4 		h0;
  REAL4 		cos2phase, sin2phase;
  static RandomParams 	*randomParams;
  static REAL4Vector  	*noise;
  INT4 			seed = 0;
  LALDetAndSource 	detAndSource;
  LALDetAMResponseSeries 	pResponseSeries = {NULL,NULL,NULL};
  REAL4TimeSeries    		Fp, Fc, Fs;
  LALTimeIntervalAndNSample 	time_info;
  UINT4				i;
  
  /* Allocate memory */
  input.B = NULL;
  input.var = NULL;
  
  LALZCreateVector( &status, &input.B, FITTOPULSARTEST_LENGTH);
  LALZCreateVector( &status, &input.var, FITTOPULSARTEST_LENGTH);
    
  noise = NULL;
  LALCreateVector( &status, &noise, FITTOPULSARTEST_LENGTH);
  LALCreateRandomParams( &status, &randomParams, seed);
  
  Fp.data = NULL;
  Fc.data = NULL;
  Fs.data = NULL;
  
  pResponseSeries.pPlus   = &(Fp);
  pResponseSeries.pCross  = &(Fc);
  pResponseSeries.pScalar = &(Fs);

  LALSCreateVector(&status, &(pResponseSeries.pPlus->data), 1);
  LALSCreateVector(&status, &(pResponseSeries.pCross->data), 1);
  LALSCreateVector(&status, &(pResponseSeries.pScalar->data), 1);
  
  input.t = tgps;
  /******** GENERATE FAKE INPUT **********/
  
  time_info.epoch.gpsSeconds     = FITTOPULSARTEST_T0;
  time_info.epoch.gpsNanoSeconds = 0;
  time_info.deltaT               = 60;
  time_info.nSample              = FITTOPULSARTEST_LENGTH;
 
  cosIota = 0.5;
  psi = 0.1;
  phase = 0.4;
  h0 = 5.0;
  
  cos2phase = cos(2.0*phase);
  sin2phase = sin(2.0*phase); 
   
  detector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];     /* use GEO 600 detector for tests */
  pulsar.equatorialCoords.longitude = 1.4653;  			/* right ascention of pulsar */
  pulsar.equatorialCoords.latitude = -1.2095;			/* declination of pulsar */
  pulsar.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;	/* coordinate system */
  pulsar.orientation = psi; 					/* polarization angle */
  strcpy(pulsar.name, "fakepulsar");			 	/* name of pulsar */
  
  detAndSource.pDetector = &detector;
  detAndSource.pSource = &pulsar;
  
  LALNormalDeviates( &status, noise, randomParams ); 
  
  LALComputeDetAMResponseSeries(&status, &pResponseSeries, &detAndSource, &time_info);
  
  for (i = 0;i < FITTOPULSARTEST_LENGTH; i++)
  {
    input.t[i].gpsSeconds = FITTOPULSARTEST_T0 + 60*i;
    input.t[i].gpsNanoSeconds = 0;
        
    input.B->data[i].re =  pResponseSeries.pPlus->data->data[i]*h0*(1.0 + cosIota*cosIota)*cos2phase
                         + 2.0*pResponseSeries.pCross->data->data[i]*h0*cosIota*sin2phase;
    input.B->data[i].im =  pResponseSeries.pPlus->data->data[i]*h0*(1.0 + cosIota*cosIota)*sin2phase
       		         - 2.0*pResponseSeries.pCross->data->data[i]*h0*cosIota*cos2phase;		     
    input.var->data[i].re = noise->data[FITTOPULSARTEST_LENGTH-i-1]*noise->data[FITTOPULSARTEST_LENGTH-i-1];
    input.var->data[i].im = noise->data[i]*noise->data[i];
  }

  input.B->length = FITTOPULSARTEST_LENGTH;
  input.var->length = FITTOPULSARTEST_LENGTH;  
  
  /*******  TEST RESPONSE TO VALID DATA  ************/
/* Test that valid data generate the correct answers */
  params.detector = detector;
  params.pulsarSrc = pulsar;
  
  params.meshH0[0] = 4.0;
  params.meshH0[1] = 0.2;
  params.meshH0[2] = 10;
  
  params.meshCosIota[0] = 0.0;
  params.meshCosIota[1] =  0.1;
  params.meshCosIota[2] =  10;

  params.meshPhase[0] = 0.0;
  params.meshPhase[1] = 0.1;
  params.meshPhase[2] = 10;
  
  params.meshPsi[0] =  0.0;
  params.meshPsi[1] =  0.1;
  params.meshPsi[2] =  10;
  
  output.mChiSquare = NULL;
  LALDCreateVector( &status, &output.mChiSquare,params.meshH0[2]*params.meshCosIota[2]*params.meshPhase[2]*params.meshPsi[2]);
  
  
  LALCoarseFitToPulsar(&status,&output, &input, &params);
 
  if(status.statusCode) 
  {
    printf("Unexpectedly got error code %d and message %s\n",
	    status.statusCode, status.statusDescription);
    return FITTOPULSARTESTC_EFLS;
  }

  
  if(output.phase > phase + params.meshPhase[2] || output.phase < phase - params.meshPhase[2])
  {
    printf("Got incorrect phase %f when expecting %f \n",
	    output.phase, phase);
    return FITTOPULSARTESTC_EFLS;
  }

  if(output.cosIota > cosIota + params.meshCosIota[2] || output.cosIota < cosIota - params.meshCosIota[2])
  {
    printf("Got incorrect cosIota %f when expecting %f \n",
	    output.cosIota, cosIota);
    return FITTOPULSARTESTC_EFLS;
  }
  
    if(output.psi > psi + params.meshPsi[2] || output.psi < psi - params.meshPsi[2])
  {
    printf("Got incorrect psi %f when expecting %f \n",
	    output.psi, psi);
    return FITTOPULSARTESTC_EFLS;
  }

  /*******  TEST RESPONSE OF LALCoarseFitToPulsar TO INVALID DATA  ************/

#ifndef LAL_NDEBUG
if ( ! lalNoDebug ) {
 
 /* Test that all the error conditions are correctly detected by the function */
 
 LALCoarseFitToPulsar(&status, NULL, &input, &params);
 
  if (status.statusCode != FITTOPULSARH_ENULLOUTPUT
       || strcmp(status.statusDescription, FITTOPULSARH_MSGENULLOUTPUT)) 
  {
    printf( "Got error code %d and message %s\n",
	    status.statusCode, status.statusDescription);
    printf( "Expected error code %d and message %s\n",
	    FITTOPULSARH_ENULLOUTPUT, FITTOPULSARH_MSGENULLOUTPUT);
    return FITTOPULSARTESTC_ECHK;
  }
 
 LALCoarseFitToPulsar(&status, &output, NULL, &params);
 
  if (status.statusCode != FITTOPULSARH_ENULLINPUT
       || strcmp(status.statusDescription, FITTOPULSARH_MSGENULLINPUT)) 
  {
    printf( "Got error code %d and message %s\n",
	    status.statusCode, status.statusDescription);
    printf( "Expected error code %d and message %s\n",
	    FITTOPULSARH_ENULLINPUT, FITTOPULSARH_MSGENULLINPUT);
    return FITTOPULSARTESTC_ECHK;
  }
 
 LALCoarseFitToPulsar(&status, &output, &input, NULL);
 
  if (status.statusCode != FITTOPULSARH_ENULLPARAMS
       || strcmp(status.statusDescription, FITTOPULSARH_MSGENULLPARAMS)) 
  {
    printf( "Got error code %d and message %s\n",
	    status.statusCode, status.statusDescription);
    printf( "Expected error code %d and message %s\n",
	    FITTOPULSARH_ENULLPARAMS, FITTOPULSARH_MSGENULLPARAMS);
    return FITTOPULSARTESTC_ECHK;
  }
 
  /* try having two input vectors of different length */
  input.var->length = 1;
  LALCoarseFitToPulsar(&status, &output, &input, &params);
 
  if (status.statusCode != FITTOPULSARH_EVECSIZE
       || strcmp(status.statusDescription, FITTOPULSARH_MSGEVECSIZE)) 
  {
    printf( "Got error code %d and message %s\n",
	    status.statusCode, status.statusDescription);
    printf( "Expected error code %d and message %s\n",
	    FITTOPULSARH_EVECSIZE, FITTOPULSARH_MSGEVECSIZE);
    return FITTOPULSARTESTC_ECHK;
  }
 
  input.var->length = FITTOPULSARTEST_LENGTH;

} /* if ( ! lalNoDebug ) */
#endif /* LAL_NDEBUG */

  /*******  CLEAN UP  ************/
  
  LALZDestroyVector(&status, &input.B);
  LALZDestroyVector(&status, &input.var); 
  LALDestroyVector(&status, &noise);  
  LALDDestroyVector(&status, &output.mChiSquare);
  LALDestroyRandomParams(&status, &randomParams);
  LALSDestroyVector(&status, &(pResponseSeries.pPlus->data));
  LALSDestroyVector(&status, &(pResponseSeries.pCross->data));
  LALSDestroyVector(&status, &(pResponseSeries.pScalar->data));  
  
  LALCheckMemoryLeaks();

  return FITTOPULSARTESTC_ENOM;
}
