/*----------------------------------------------------------------------- 
 * 
 * File Name: CoherentInspiralInput.c
 *
 * Author: Seader, S. E.  Brown, D.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <regex.h>
#include <time.h>
#include <math.h>

#include <lal/LALRCSID.h>
#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/FrameStream.h>
#include <lal/DataBuffer.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/FindChirpBCV.h>
#include <lal/FindChirpBCVSpin.h>
#include <lal/FindChirpChisq.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/DetectorSite.h>
#include <lal/Random.h>
#include <lal/LALInspiral.h>
#include <lal/CoherentInspiral.h>
#include <lal/LALStatusMacros.h>

NRCSID( COHERENTINSPIRALINPUTC, "$Id$");

/* JC: Please don't use rint ... it is not standard C89 */
/* double rint(double x); */
#define rint(x) floor((x)+0.5)

static void
LALFindChirpCreateCoherentInput( 
     LALStatus                  *status,
     COMPLEX8TimeSeries         **coherentInputData,
     COMPLEX8TimeSeries         *input,
     SnglInspiralTable          *tmplt,
     REAL4                      coherentSegmentLength,
     INT4                       corruptedDataLength 
     )
{
  COMPLEX8TimeSeries      *cohInputData = NULL;
  LIGOTimeGPS              end_time;
  UINT8                    tmpltID = 0;
  UINT8                    eventID = 0;
  INT4                     numPoints = 0;
  REAL4                    cohSegLength = 0.0;
  INT4                     inputEpochSeconds = 0;
  INT4                     inputEpochNanoSeconds = 0;
  REAL8                    deltaT = 0.0;
  INT4                     crupDataLength = 0;
  INT4                     cohSegEnd = 0;
  INT4                     cohSegStart = 0;
  INT4                     nonCorruptEnd = 0;
  INT4                     nonCorruptStart = 0;

  INITSTATUS( status, "LALFindChirpCreateCoherentInput",
	      COHERENTINSPIRALINPUTC );
  ATTATCHSTATUSPTR( status );

  /* Put series of ASSERT statements here to check validity of input params */


  /* Get necessary info from input structures */
  cohInputData = *coherentInputData = (COMPLEX8TimeSeries *)
    LALCalloc(1, sizeof(COMPLEX8TimeSeries) );
  end_time = tmplt->end_time;
  eventID  = tmplt->event_id->id;
  numPoints = input->data->length;
  cohSegLength = coherentSegmentLength;
  inputEpochSeconds = input->epoch.gpsSeconds;
  inputEpochNanoSeconds = input->epoch.gpsNanoSeconds;
  deltaT = input->deltaT;
  crupDataLength = corruptedDataLength;

  /* Now determine if the event lies in corrupted data */

  nonCorruptStart = crupDataLength;
  nonCorruptEnd = numPoints - crupDataLength;
  cohSegStart = (INT4) rint( ((end_time.gpsSeconds + end_time.gpsNanoSeconds*1.0e-9) - (inputEpochSeconds + inputEpochNanoSeconds*1.0e-9)-cohSegLength)/deltaT - 1.0);
  cohSegEnd = cohSegStart + 2* (INT4)cohSegLength/deltaT + 1;

  if( cohSegEnd < nonCorruptEnd && cohSegStart >= nonCorruptStart )
    {
      /* The coherent chunk is outside of the corrupted data */

      INT4 fullCohSegLength = 0;
      fullCohSegLength = 2*cohSegLength/deltaT;

      LALCCreateVector(status->statusPtr, &(cohInputData->data), fullCohSegLength);
      memcpy(cohInputData->data->data, &(input->data->data[cohSegStart]), 2*cohSegLength/deltaT);
    }
  else if( (cohSegStart >= nonCorruptEnd) || (cohSegEnd <= nonCorruptStart) )
    {
      /* The coherent chunk is entirely within corrupted data */
      /* The cohInputData will return a null value on exit of function */ 
    }
  else if( cohSegEnd >= nonCorruptEnd )
    {
      /* The coherent chunk is in corrupted data at the end of the segment */

      INT4 overlap          = 0;
      INT4 fullCohSegLength = 0;

      overlap = cohSegEnd - nonCorruptEnd;
      fullCohSegLength = 2*cohSegLength/deltaT - overlap;

      LALCCreateVector( status->statusPtr, &(cohInputData->data), fullCohSegLength );
      memcpy( cohInputData->data->data, &(input->data->data[cohSegStart]),fullCohSegLength);
    }
  else if( cohSegStart <= nonCorruptStart )
    {
      /* The coherent chunk is in corrupted data at the start of the segment */

      INT4 overlap          = 0;
      INT4 fullCohSegLength = 0;

      overlap = nonCorruptStart - cohSegStart;
      fullCohSegLength = 2*cohSegLength/deltaT - overlap;

      LALCCreateVector( status->statusPtr, &(cohInputData->data), fullCohSegLength );
      memcpy( cohInputData->data->data, &(input->data->data[cohSegStart + overlap]), fullCohSegLength);
    }
}


static void
TestStatus (LALStatus *status, const char *expectedCodes, int exitCode);

static void
ClearStatus (LALStatus *status);

int verbose = 1;

int main( int argc, char *argv[] )
{
  /* Test LALFindChirpCreateCoherentInput() */
  static LALStatus       status;
  COMPLEX8TimeSeries    *cohInputData = NULL;
  COMPLEX8TimeSeries    *input = NULL;
  SnglInspiralTable     *tmplt = NULL;
  INT4                   corruptedDataLength = 0;
  REAL4                  cohSegLength = 2.0;
  UINT4                  numPoints = 524288;

  FrStream     *frStream = NULL;
  FrChanIn      frChan;

  if( !(tmplt = (SnglInspiralTable *) LALCalloc( 1, sizeof(SnglInspiralTable)) ) )
    {
      fprintf( stderr, "could not allocate memory for tmplt\n");
      exit( 1 );
    }
  tmplt->event_id = (EventIDColumn *) LALCalloc( 1, sizeof(EventIDColumn));
  tmplt->end_time.gpsSeconds = 751976330;
  tmplt->end_time.gpsNanoSeconds = 234567999;
  tmplt->event_id->id = 1;
  corruptedDataLength = numPoints/4;

  fprintf(stdout,"end_time(s) = %d\nend_time(ns) = %d\nevent_id->id = %d\n",tmplt->end_time.gpsSeconds,tmplt->end_time.gpsNanoSeconds,tmplt->event_id->id);
  fprintf(stdout,"cohSegLength = %f\n",cohSegLength);

  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  input = (COMPLEX8TimeSeries *) LALCalloc( 1, sizeof(COMPLEX8TimeSeries));
  LALCCreateVector(&status, &(input->data), numPoints);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  LALFrOpen(&status,&frStream,"/home/sseader/LIGO/INSPIRAL/S3/playground/figures6/","L1-INSPIRAL-751976300-512.gwf");
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  if(!frStream)
    {
      fprintf(stderr,"The file %s does not exist - exiting...\n",
	      "L1-INSPIRAL-751976300-512.gwf");
      exit( 1 );
    }
  frChan.name = "L1:LSC-AS_Q_CData_0";
  LALFrGetCOMPLEX8TimeSeries( &status, input, &frChan, frStream);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  LALFrClose( &status, &frStream);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  fprintf(stdout,"calling LALFindChirpCreateCoherentInput...\n");
  LALFindChirpCreateCoherentInput( &status, &cohInputData, input, tmplt, cohSegLength, corruptedDataLength );
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  fprintf(stdout,"same data passed out: %f %f\n",cohInputData->data->data[0].re,cohInputData->data->data[0].im);

  goto cleanexit;

 cleanexit:
  LALCDestroyVector( &status, &(input->data) );
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  LALFree( input );

  exit( 0 );
}

/*
 * TestStatus ()
 *
 * Routine to check that the status code status->statusCode agrees with one of
 * the codes specified in the space-delimited string ignored; if not,
 * exit to the system with code exitcode.
 *
 */
static void
TestStatus (
    LALStatus  *status, 
    const char *ignored, 
    int         exitcode
           )
{
  char  str[64];
  char *tok;

  if (verbose)
  {
    REPORTSTATUS (status);
  }

  if (strncpy (str, ignored, sizeof (str)))
  {
    if ((tok = strtok (str, " ")))
    {
      do
      {
        if (status->statusCode == atoi (tok))
        {
          return;
        }
      }
      while ((tok = strtok (NULL, " ")));
    }
    else
    {
      if (status->statusCode == atoi (tok))
      {
        return;
      }
    }
  }

  fprintf (stderr, "\nExiting to system with code %d\n", exitcode);
  exit (exitcode);
}


/*
 *
 * ClearStatus ()
 *
 * Recursively applies DETATCHSTATUSPTR() to status structure to destroy
 * linked list of statuses.
 *
 */
void
ClearStatus (
    LALStatus   *status
            )
{
  if (status->statusPtr)
  {
    ClearStatus      (status->statusPtr);
    DETATCHSTATUSPTR (status);
  }
}
