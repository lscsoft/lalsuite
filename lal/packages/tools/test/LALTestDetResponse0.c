/* <lalVerbatim file="LALTestDetResponse0CV">
   $Id$
   </lalVerbatim
*/

/*
<lalLaTeX>

\subsection{Program {\texttt{LALTestDetResponse0.c}}
\label{ss:LALTestDetResponse0.c}

\subsubsection*{Usage}

\subsubsection*{Description}

Performs zeroth-order test of \texttt{LALComputeDetAMResponse()} and
\texttt{LALComputeDetAMResponseSeries()}. 

\subsubsection*{Exit codes}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALComputeDetAMResponse()
\end{verbatim}

\subsubsection*{Notes}

</lalLaTeX> 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <lal/LALConfig.h>

#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALStdlib.h>
#include <lal/DetectorSite.h>
#include <lal/TimeDelay.h>
#include <lal/DetResponse.h>

NRCSID( LALTESTDETRESPONSE0C, "$Id$" );

int lalDebugLevel = 0;

static void REAL4VectorSubtraction(REAL4Vector *pA, REAL4Vector *pB,
                                   REAL4Vector *pAminusB);
static REAL4 REAL4VectorRMS(REAL4Vector *pVector);

int main(int argc, char *argv[])
{
  static LALStatus  status;
  LALSource         pulsar;
  LALFrDetector     frdet;    /* Framelib detector info */
  LALDetector       detector;
  LIGOTimeGPS       gps;
  LALDetAndSource   det_and_pulsar;
  LALDetAMResponse  am_response;

  LALDetAMResponseSeries    am_response_series = {NULL,NULL,NULL};
  REAL4TimeSeries           plus_series, cross_series, scalar_series;
  /* REAL4Vector               diffVector; */
  LALTimeIntervalAndNSample time_info;

  UINT4 i;


  
  if (argc == 2)
    lalDebugLevel = atoi(argv[1]);

  /*
   * Set up a source at (RA=0, Dec=0, orientation=0)
   */
  strcpy(pulsar.name, "TEST PULSAR");
  pulsar.equatorialCoords.longitude = 0.;  /* RA */
  pulsar.equatorialCoords.latitude  = 0.;  /* Dec */
  pulsar.equatorialCoords.system    = COORDINATESYSTEM_EQUATORIAL;
  pulsar.orientation                = 0.;  /* orientation */

  /*
   * As per John Whelan's suggestion, directly create a LALDetector 
   * structure rather than starting with a LALFrDetector structure.
   */
  detector.location[0] = LAL_AWGS84_SI;
  detector.location[1] = 0.;
  detector.location[2] = 0.;
  detector.response[0][0] = 0.;
  detector.response[1][1] = 0.5;
  detector.response[2][2] = -0.5;
  detector.response[0][1] = detector.response[1][0] = 0.;
  detector.response[0][2] = detector.response[2][0] = 0.;
  detector.response[1][2] = detector.response[2][1] = 0.;
  detector.type = LALDETECTORTYPE_ABSENT;

  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr,
              "LALTestDetResponse0: LALCreateDetector failed, line %i, %s\n",
              __LINE__, LALTESTDETRESPONSE0C);
      REPORTSTATUS(&status);
      return status.statusCode;
    }
  /* REPORTSTATUS(&status); */

  /*
   * Expect the location vector to be (R, 0, 0): R = radius of Earth
   *                                                 at Equator
   */
  if (lalDebugLevel > 0)
    printf("Det #1 location: (%7.4e, %7.4e, %7.4e)\n",
           detector.location[0], detector.location[1],
           detector.location[2]);

  
  /*
   * Set a GPS time that's close to 0h GMST1. (Found this by trial and
   * error.) 
   */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;


  /*
   * Stuff Detector and Source into structure
   */
  det_and_pulsar.pDetector = &detector;
  det_and_pulsar.pSource   = &pulsar;

  /*
   * Compute instantaneous AM response
   */
  LALComputeDetAMResponse(&status, &am_response, &det_and_pulsar, &gps);

  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr,
              "LALTestDetResponse0: error in LALComputeDetAMResponse, line %i, %s\n",
              __LINE__, LALTESTDETRESPONSE0C);
      REPORTSTATUS(&status);
      return status.statusCode;
    }
  
  /*
   * Expect plus = 1., cross = 0., scalar = undef.
   */
  if ((REAL4)fabs((double)(am_response.plus - (REAL4)1.)) > LAL_REAL4_EPS)
    {
      fprintf(stderr,
              "LALTestDetResponse0: ERROR: expect plus=1, got plus=%g, cross=%g\n",
              am_response.plus, am_response.cross);
      return 1;
    }
  
  if ((REAL4)fabs((double)(am_response.cross - (REAL4)0.)) > LAL_REAL4_EPS)
    {
      fprintf(stderr,
              "LALTestDetResponse0: ERROR: expect cross=0, got plus=%g, cross=%g\n",
              am_response.plus, am_response.cross);
      return 2;
    }


  /*
   * Compute a time series AM response
   */
  if (lalDebugLevel > 0)
    {
      printf("Starting vector test\n");
    }

  plus_series.data = NULL;
  cross_series.data = NULL;
  scalar_series.data = NULL;
  
  am_response_series.pPlus   = &(plus_series);
  am_response_series.pCross  = &(cross_series);
  am_response_series.pScalar = &(scalar_series);

  LALSCreateVector(&status, &(am_response_series.pPlus->data), 1);
  LALSCreateVector(&status, &(am_response_series.pCross->data), 1);
  LALSCreateVector(&status, &(am_response_series.pScalar->data), 1);

  if (lalDebugLevel > 0)
    {
      printf("am_response_series.pPlus->data->length = %d\n",
             am_response_series.pPlus->data->length);
      printf("am_response_series.pCros->data->length = %d\n",
             am_response_series.pCross->data->length);
      printf("am_response_series.pScalar->data->length = %d\n",
             am_response_series.pScalar->data->length);
    }
    
  time_info.epoch.gpsSeconds     = 61094;
  time_info.epoch.gpsNanoSeconds = 640000000;
  time_info.deltaT               = 20.5;
  time_info.nSample              = 13;

  LALComputeDetAMResponseSeries(&status,
                                &am_response_series,
                                &det_and_pulsar,
                                &time_info);

  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr,
              "LALTestDetResponse0: error in LALComputeDetAMResponseSeries, line %i, %s\n",
              __LINE__, LALTESTDETRESPONSE0C);
      return status.statusCode;
    }
  
  if (lalDebugLevel > 0)
    {
      printf("Done computing AM response vectors\n");

      printf("am_response_series.pPlus->data->length = %d\n",
             am_response_series.pPlus->data->length);
      printf("am_response_series.pCross->data->length = %d\n",
             am_response_series.pCross->data->length);
      printf("am_response_series.pScalar->data->length = %d\n",
             am_response_series.pScalar->data->length);
    }

  if (lalDebugLevel >= 8)
    {
      printf("plus: (");
      for (i = 0; i < time_info.nSample; ++i)
        {
          printf("%1.6e, ", am_response_series.pPlus->data->data[i]);
        }
      printf(")\n");

      printf("cross: (");
      for (i = 0; i < time_info.nSample; ++i)
        {
          printf("%1.6e, ", am_response_series.pCross->data->data[i]);
        }
      printf(")\n");

      printf("scalar: (");
      for (i = 0; i < time_info.nSample; ++i)
        {
          printf("%1.6e, ", am_response_series.pScalar->data->data[i]);
        }
      printf(")\n");


      /* print out quadrature sum of plus- and cross-response */
      printf("sqrt(PLUS^2 + CROSS^2): (");
      for (i = 0; i < time_info.nSample; ++i)
        {
          printf("%1.6e, ", sqrt(am_response_series.pPlus->data->data[i] *
                                 am_response_series.pPlus->data->data[i] +
                                 am_response_series.pCross->data->data[i] *
                                 am_response_series.pCross->data->data[i]));
        }
      printf(")\n");
    }


  LALSDestroyVector(&status, &(am_response_series.pPlus->data));
  LALSDestroyVector(&status, &(am_response_series.pCross->data));
  LALSDestroyVector(&status, &(am_response_series.pScalar->data));

  LALCheckMemoryLeaks();

  return 0;
}


/*
 * subtracts two REAL4Vectors; user must do all allocation beforehand
 */
static void REAL4VectorSubtraction(REAL4Vector *pA, REAL4Vector *pB,
                                   REAL4Vector *pAminusB)
{
  UINT4 i;

  /* Check for compatible dimensions */
  if ((pA->length != pB->length) || (pAminusB->length != pA->length) ||
      (pAminusB->length != pB->length)) 
    {
      fprintf(stderr, "VectorSubtraction: ERROR: incompatible dimensions\n");
      exit(13);
    }

  for (i = 0; i < pA->length; ++i)
    {
      pAminusB->data[i] = pA->data[i] - pB->data[i];
    }

    return;
}

static REAL4 REAL4VectorRMS(REAL4Vector *pVector)
{
  UINT4 i;
  REAL4 result = 0.;

  for (i = 0; i < pVector->length; ++i)
    {
      result += pVector->data[i] * pVector->data[i];
    }

  result /= pVector->length;

  return sqrt(result);
}
