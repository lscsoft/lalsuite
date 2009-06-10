/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton, Sean Seader
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

/*-----------------------------------------------------------------------
 *
 * File Name: CoherentInspiralInputTest.c
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
#include <getopt.h>

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

#include <lal/LALRCSID.h>
NRCSID (COHERENTINSPIRALINPUTTESTC,"$Id$");

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
  tmplt->end_time.gpsSeconds = 751976331;
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
  LALCheckMemoryLeaks();

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
