/*
*  Copyright (C) 2007 Tania Regimbau
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

/*
 * readframes.c 
 *
 * Tania Regimbau <regimbau@obs-nice.fr>  
 *
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>

#include <lal/LALFrameL.h>

#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/DetectorSite.h>
#include <lal/LALCache.h>
#include <lal/LALFrStream.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
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
  LALCache *frCache;
  LALFrStream *frStream;
  FrChanIn frChanIn;
 
 /* input data segment */
  LIGOTimeGPS gpsStartTime;
  
  /* data structures for PSDs */

  REAL8FrequencySeries PSD;
  LALUnit PSDUnits = {1,{0,0,1,0,0,2,0},{0,0,0,0,0,0,0}};

  /* error handler */
  status.statusPtr = NULL;

  lal_errhandler = LAL_ERR_EXIT;
  
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
  frCache = XLALCacheImport( frameCache );
  LAL_CALL( LALFrCacheOpen( &status, &frStream, frCache), &status );

   /* set the mode of the frame stream to fail on gaps or time errors */
  frStream->mode = LAL_FR_STREAM_VERBOSE_MODE;
                
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
