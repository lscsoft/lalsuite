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


static void
LALFindChirpCreateCoherentInput( 
     LALStatus                  *status,
     COMPLEX8TimeSeries         **coherentInputData,
     COMPLEX8TimeSeries         *input,
     SnglInspiralTable          *tmplt,
     REAL4                      corruptedDataLength 
     )
{
  LIGOTimeGPS end_time;
  UINT8 tmpltID = 0;
  UINT8 eventID = 0;


  /* Put series of ASSERT statements here to check validity of input params */

  end_time = tmplt->end_time;
  tmpltID  = tmplt->tmpltID;
  eventID  = tmplt->event_id->id;

  fprintf(stdout, "end_time.gpsSeconds=%d\n end_time.gpsNanoSeconds=%d\n tmpltID=%u\n eventID=%u\n",end_time.gpsSeconds,end_time.gpsNanoSeconds,tmpltID,eventID);

}


int main( int argc, char *argv[] )
{
  /* Test LALFindChirpCreateCoherentInput() */
  LALStatus             status;
  COMPLEX8TimeSeries    *cohInputData = NULL;
  COMPLEX8TimeSeries    *input = NULL;
  SnglInspiralTable     *tmplt = NULL;
  REAL4                 corruptedDataLength = 0;

  FrStream     *frStream = NULL;
  FrChanIn      frChan;


  tmplt->tmpltID = 1;
  tmplt->end_time.gpsSeconds = 712444333;
  tmplt->end_time.gpsNanoSeconds = 542377;
  tmplt->event_id->id = 1;

  corruptedDataLength = 2;

  LALFrOpen(&status,&frStream,NULL,"L1-INSPIRAL-751976300-512.gwf");
  if(!frStream)
    {
      fprintf(stderr,"The file %s does not exist - exiting...\n",
	      "L1-INSPIRAL-751976300-512.gwf");
      exit( 1 );
    }
  frChan.name = "L1:LSC-AS_Q_CData_0";
  LALFrGetCOMPLEX8TimeSeries( &status, input, &frChan, frStream);
  LALFrClose( &status, &frStream);

  LALFindChirpCreateCoherentInput( &status, &cohInputData, input, tmplt, corruptedDataLength );

  exit( 0 );
}

