/*
 * readframes.c 
 *
 * Tania Regimbau <regimbau@obs-nice.fr>  
 *
 *
 * $Id$
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>

#include <FrameL.h>

#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/DetectorSite.h>
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALRCSID.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/PrintFTSeries.h>
#include <lal/PrintVector.h>
#include <lal/StreamInput.h>
#include <lal/Units.h>
#include <lal/Window.h>
#include <lal/VectorOps.h>
#include <lalapps.h>

#include "stochastic_virgo.h"

NRCSID (READFRAMESC, "$Id$");
RCSID ("$Id$");

/* cvs info */
#define CVS_ID "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "lalapps_readframes"


CHAR frameCache [200]= "test.cache";
CHAR channel[LALNameLength]= "H1_PSD";
UINT8 startTime = 815412745;
INT4 segmentDuration = 60;
INT4 length=7041;
REAL8 deltaF=0.25;


INT4 main()
 {

  /* variable declarations */

  /* status pointer */
  LALStatus status;
  
  /* frame variables */
  FrCache *frCache;
  FrStream *frStream;
  FrChanIn frChanIn;
 
 /* input data segment */
  LIGOTimeGPS gpsStartTime;
  
  /* data structures for PSDs */

  REAL8FrequencySeries PSD;
  LALUnit PSDUnits = {1,{0,0,1,0,0,2,0},{0,0,0,0,0,0,0}};

  /* error handler */
  status.statusPtr = NULL;

  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "1" );

  
   /* initialize gps time structure */  
  gpsStartTime.gpsSeconds = startTime;
  gpsStartTime.gpsNanoSeconds = 0.;  

  /* set metadata fields for PSDs */
  strncpy(PSD.name, "PSD", LALNameLength);
  PSD.sampleUnits = PSDUnits;
  PSD.epoch  = gpsStartTime;
  PSD.deltaF = deltaF;

  /* allocate memory for PSDs */
  PSD.data = NULL;
  LAL_CALL(LALDCreateVector(&status,&(PSD.data),length),
	    &status );
  memset(PSD.data->data, 0, 
         PSD.data->length * sizeof(*PSD.data->data));

  /* set channels */
  frChanIn.name = channel;
  frChanIn.type = ADCDataChannel;
  
  /* open frame cache */
  fprintf(stdout, "Opening frame cache...\n");
  frCache=NULL;
  frStream=NULL;
  LAL_CALL( LALFrCacheImport( &status, &frCache, frameCache), &status );
  LAL_CALL( LALFrCacheOpen( &status, &frStream, frCache), &status );

   /* set the mode of the frame stream to fail on gaps or time errors */
  frStream->mode = LAL_FR_VERBOSE_MODE;
                
  /* read data */	
  fprintf(stdout, "Reading in channel \"%s\"...\n", frChanIn.name); 
  PSD.data = NULL;      
  LAL_CALL(LALFrSeek( &status, &(gpsStartTime), frStream), &status );
  LAL_CALL(LALFrGetREAL8FrequencySeries(&status,&PSD,&frChanIn,frStream),&status );
               
  /* print */
  fprintf(stdout, "Printing to ascii file...\n"); 
  LALDPrintFrequencySeries(&PSD, "PSD.dat");     

   /* cleanup */
  fprintf(stdout, "close frames and exit\n");
  LAL_CALL( LALFrClose( &status, &frStream), &status );  
  LAL_CALL(LALDDestroyVector(&status, &(PSD.data)),&status );
 }
