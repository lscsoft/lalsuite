/*
 * frames.c - produces downsampled frames for the stochastic pipeline
 *
 *
 * Tania Regimbau <Tania.Regimbau@astro.cf.ac.uk>
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

#include <lalapps.h>

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
#include <lal/ResampleTimeSeries.h>
#include <lal/StreamInput.h>
#include <lal/Units.h>
#include <lal/Window.h>
#include <lal/IIRFilter.h>
#include <lal/BandPassTimeSeries.h>

NRCSID (FRAMESC, "$Id$");
RCSID ("$Id$");

/* cvs info */
#define CVS_ID "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "lalapps_frames"

/* variables for getopt options parsing */
char *optarg;
int optind;


/* sampling parameters */
INT4 sampleRate = 16384;
INT4 resampleRate = 1024;
REAL8 deltaF = 0.25;

/* data parameters */
LIGOTimeGPS gpsStartTime;
UINT8 startTime;
UINT8 stopTime;
INT4 segmentDuration = 64;
INT4 segmentShift = 60;
CHAR frameCache1[100] = "cachefiles/H-730793097.cache";
CHAR frameCache2[100] = "cachefiles/L-730793097.cache";
CHAR channel1[LALNameLength]= "H1:LSC-AS_Q";
CHAR channel2[LALNameLength]= "L1:LSC-AS_Q";
CHAR ifo1[LALNameLength] = "H1";
CHAR ifo2[LALNameLength] = "L1";
INT4 site1 = 0;
INT4 site2 = 1;


/* output file */
CHAR outputFilePath[200] = "/usr1/tregimba/";

INT4 main(INT4 argc, CHAR *argv[])
 {
  /* variable declarations */

  /* status pointer */
  LALStatus status;

  /* frame parameters */
  FrOutPar opar1, opar2;

  /* output file */
  FILE *out;
  CHAR outputFilename[LALNameLength];
  CHAR frameName[LALNameLength];
  /* counters */
  INT4 i,j, segLoop;


  /* input data segment */
  INT4 numSegments;
  INT4 segmentLength;
  INT4 segmentShift;
  INT4 padData = 0;
  INT4 buffer;
  ReadDataPairParams streamParams;
  StreamPair streamPair;
  REAL4TimeSeries segment1,segment2;


  /* error handler */
  status.statusPtr = NULL;

  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "1" );

  /* parse command line options */
  parseOptions(argc, argv);


  /* read parameters into input parameter file */
 
  fscanf(stdin,"%d\n",&startTime, &startTime);
  fscanf(stdin,"%d\n",&stopTime, &stopTime);
  fscanf(stdin,"%s\n%s\n",&frameCache1,&frameCache2);
   

  /* get number of segments */
  buffer = (segmentDuration - segmentShift) / 2;
  numSegments = (INT4)(stopTime - startTime + 2 * buffer) / segmentShift };
  

  /* set length for data segments */
  segmentLength = segmentDuration * resampleRate;

  /* set metadata fields for data segments */
     
  strncpy(segment1.name, "segment1", LALNameLength);
  strncpy(segment2.name, "segment2", LALNameLength);
  segment1.sampleUnits = segment2.sampleUnits = lalADCCountUnit;
  segment1.epoch = segment2.epoch = gpsStartTime;
  segment1.deltaT = segment2.deltaT = 1./(REAL8)resampleRate;
  segment1.f0 = segment2.f0 = 0;

 

  segment1.data = segment2.data = NULL;
  LAL_CALL( LALSCreateVector( &status, &(segment1.data), segmentLength), 
            &status );
  LAL_CALL( LALSCreateVector( &status, &(segment2.data), segmentLength), 
            &status );
  memset( segment1.data->data, 0,
          segment1.data->length * sizeof(*segment1.data->data));
  memset( segment2.data->data, 0,
          segment2.data->length * sizeof(*segment2.data->data));

  /* set segment input parameters */
  streamParams.duration = segmentDuration;
  streamParams.frameCache1 = frameCache1;
  streamParams.frameCache2 = frameCache2;
  streamParams.ifo1 = ifo1;
  streamParams.ifo2 = ifo2;
  streamParams.channel1 = channel1;
  streamParams.channel2 = channel2;
  streamParams.startTime = startTime;
  streamParams.buffer = padData;
  streamParams.sampleRate = sampleRate;
  streamParams.resampleRate = resampleRate;

  /* set stream data structures */
  streamPair.stream1 = &segment1;
  streamPair.stream2 = &segment2;

  
   
   /** loop over segments **/
        
   lal_errhandler = LAL_ERR_RTRN;

   for (segLoop = 0; segLoop < numSegments; segLoop++)
    {
     /* define segment epoch */
     gpsStartTime.gpsSeconds = startTime + (segLoop * segmentShift);
     segment1.epoch = segment2.epoch = gpsStartTime;
 

     /* read data and downsample */
     streamParams.startTime = gpsStartTime.gpsSeconds - buffer;
     LAL_CALL(readDataPair(&status, &streamPair, &streamParams), &status);
       
     /* skip segment if data not found or corrupted with 0 values */           
     if ((status.statusCode !=0)||
          (segmentTemp1.data==NULL)||(segmentTemp2.data==NULL))
       {
	clear_status(&status);
        if (segLoop < (numSegments - 1)) continue; 
	else break;   
       }
       
      /* set frame parameters */ 
      LALSnprintf( frameName, LALNameLength, 
                     "-%d-%d.gwf",(INT4)startTime,segmentDuration);
      opar1 = { ifo1, frameName, ADCDataChannel, 1, 0, 0 };
      opar2 = { ifo2, frameName, ADCDataChannel, 1, 0, 0 };  
  
      /* write to frames */
      LAL_CALL(LALFrWriteINT4TimeSeries( &status,&(segment1.data),&opar1 ),&status);
      LAL_CALL(LALFrWriteINT4TimeSeries( &status,&(segment2.data),&opar2 ),&status);
     
      
     
    }
       
   lal_errhandler = LAL_ERR_EXIT;

   /* cleanup */

   LAL_CALL( LALDestroyVector(&status, &(segment1.data)), &status );
   LAL_CALL( LALDestroyVector(&status, &(segment2.data)), &status );
   

  return 0;
 }


/* parse command line options */
void parseOptions(INT4 argc, CHAR *argv[])
 {
  int c = -1;

  while(1)
   {
    static struct option long_options[] =
     {
      /* options that set a flag */
      
      /* options that don't set a flag */
      {"help", no_argument, 0, 'h'},
      {"gps-start-time", required_argument, 0, 't'},
      {"gps-end-time", required_argument, 0, 'T'},
      {"segment-duration", required_argument, 0, 'l'},
      {"sample-rate", required_argument, 0, 'A'},
      {"resample-rate", required_argument, 0, 'a'},
      {"ifo-one", required_argument, 0, 'i'},
      {"ifo-two", required_argument, 0, 'I'},
      {"frame-cache-one", required_argument, 0, 'd'},
      {"frame-cache-two", required_argument, 0, 'D'},
      {"debug-level", required_argument, 0, 'z'},
      {"version", no_argument, 0, 'V'},
      {0, 0, 0, 0}
     };

    /* getopt_long stores the option here */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long(argc, argv, 
                    "ht:T:l:A:a:i:I:d:D:z:V",
 		    long_options, &option_index);

    if (c == -1)
     {
      /* end of options, break loop */
      break;
     }

    switch(c)
     {
      case 0:
             /* If this option set a flag, do nothing else now. */
             if (long_options[option_index].flag != 0)
              break;
             printf ("option %s", long_options[option_index].name);
             if (optarg)
              printf (" with arg %s", optarg);
             printf ("\n");
             break;

      case 'h':
               /* HELP!!! */
               displayUsage(0);
               break;

      case 't':
               /* start time */
	       startTime = atoi(optarg);
	       break;

      case 'T':
	       /* stop time */
	       stopTime = atoi(optarg);
	       break;

      case 'l':
	       /* duration */
	       segmentDuration = atoi(optarg);
	       break;
      case 'A':
               /* sample rate */
               sampleRate = atoi(optarg);
               break;

      case 'a':
	       /* resampling */
	       resampleRate = atoi(optarg);
	       break;

      case 'i':
	       /* ifo for first stream */
	       strncpy(ifo1, optarg, LALNameLength);

	       /* set site and channel */
	       if (strncmp(ifo1, "H1", 2) == 0)
		{
		 site1 = 0;
		 strncpy(channel1, "H1:LSC-AS_Q", LALNameLength);
		}
	       else if (strncmp(ifo1, "H2", 2) == 0)
		{
		 site1 = 0;
		 strncpy(channel1, "H2:LSC-AS_Q", LALNameLength);
		}
	       else if (strncmp(ifo1, "L1", 2) == 0)
		{
		 site1 = 1;
		 strncpy(channel1, "L1:LSC-AS_Q", LALNameLength);
	        }
	       else
		{
		 fprintf(stderr, "First IFO not recognised...\n");
		 exit(1);
		}

	       break;

       case 'I':
		/* ifo for second stream */
		strncpy(ifo2, optarg, LALNameLength);

		/* set site and channel */
		if (strncmp(ifo2, "H1", 2) == 0)
		 {
		  site2 = 0;
		  strncpy(channel2, "H1:LSC-AS_Q", LALNameLength);
		 }
		else if (strncmp(ifo2, "H2", 2) == 0)
		 {
		  site2 = 0;
		  strncpy(channel2, "H2:LSC-AS_Q", LALNameLength);
		 }
		else if (strncmp(ifo2, "L1", 2) == 0)
		 {
		  site2 = 1;
		  strncpy(channel2, "L1:LSC-AS_Q", LALNameLength);
		 }
		 else
		  {
		   fprintf(stderr, "Second IFO not recognised...\n");
		   exit(1);
		  }

		 break;

	case 'd':
         	 /* data cache one */
                 strncpy(frameCache1, optarg, LALNameLength);
        	 break;

        case 'D':
                /* data cache two */
                strncpy(frameCache2, optarg, LALNameLength);
                break;
   
        case 'z':
		/* set debug level */
		set_debug_level( optarg );
		break;

	case 'V':
		/* display version info and exit */
		fprintf(stdout, "Standalone SGWB Search Engine\n" CVS_ID "\n");
		exit(0);
	        break;

		default:
		displayUsage(1);
       }
    }

   if (optind < argc)
    {
     displayUsage(1);
    }

  return;
}

/* display program usage */
void displayUsage(INT4 exitcode)
 {
  fprintf(stderr, "Usage: pipeline [options]\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, " -h                    print this message\n");
  fprintf(stderr, " -V                    display version\n");
  fprintf(stderr, " --verbose             verbose mode\n");
  fprintf(stderr, " -z                    set lalDebugLevel\n");
  fprintf(stderr, " -t                    GPS start time\n");
  fprintf(stderr, " -T                    GPS stop time\n");
  fprintf(stderr, " -l                    segment duration\n");
  fprintf(stderr, " -A                    sample rate\n");
  fprintf(stderr, " -a                    resample rate\n");  
  fprintf(stderr, " -i                    ifo for first stream\n");
  fprintf(stderr, " -I                    ifo for second stream\n");
  fprintf(stderr, " -d                    cache file for first stream\n");
  fprintf(stderr, " -D                    cache file for second stream\n");
       
  exit(exitcode);
}



/* function to read data in frames */
void readDataPair(LALStatus *status,
		  StreamPair *streamPair,
		  ReadDataPairParams *params)
 {
  /* counters */
  INT4 i;

  /* variables */
  FrCache *frCache1 = NULL;
  FrStream *frStream1 = NULL;
  FrCache *frCache2 = NULL;
  FrStream *frStream2 = NULL;
  FrChanIn frChanIn1, frChanIn2;
  REAL4TimeSeries dataStream1, dataStream2;
  ResampleTSParams resampleParams;
  LIGOTimeGPS bufferStartTime;
  UINT8 startTime;
  INT4 buffer;
  INT4 resampleRate, sampleRate;

  /* read parameters */
  startTime = params->startTime;
  buffer = params->buffer;
  resampleRate = params->resampleRate;
  sampleRate = params->sampleRate;

  /* initialise status pointer */
  INITSTATUS( status, "readDataPair", STOCHASTICC );
  ATTATCHSTATUSPTR( status );

  /* buffer start time */
  bufferStartTime.gpsSeconds = startTime - buffer;
  bufferStartTime.gpsNanoSeconds = 0;

  /* set channels */
  frChanIn1.name = params->channel1;
  frChanIn2.name = params->channel2;
  frChanIn2.type = ADCDataChannel;
  frChanIn1.type = ADCDataChannel;

  /* initial data structures */
  dataStream1.epoch =  dataStream2.epoch = bufferStartTime;



  /* allocate memory */
  dataStream1.data = dataStream2.data = NULL;
  LALSCreateVector( status->statusPtr, &(dataStream1.data),
                    sampleRate * (params->duration + (2 * buffer)));
  CHECKSTATUSPTR (status);	
  LALSCreateVector( status->statusPtr, &(dataStream2.data), 
                    sampleRate * (params->duration + (2 * buffer)));
  CHECKSTATUSPTR (status);	
  memset( dataStream1.data->data, 0, 
          dataStream1.data->length * sizeof(*dataStream1.data->data));
  memset( dataStream2.data->data, 0, 
          dataStream2.data->length * sizeof(*dataStream2.data->data));

 

  /* open first frame cache */
  LALFrCacheImport(status->statusPtr, &frCache1, params->frameCache1);
  CHECKSTATUSPTR (status);
  LALFrCacheOpen(status->statusPtr, &frStream1, frCache1);
  CHECKSTATUSPTR (status);
       

  /* read first channel */
  LALFrSeek(status->statusPtr, &(bufferStartTime), frStream1);
  CHECKSTATUSPTR (status);
  LALFrGetREAL4TimeSeries(status->statusPtr,
                          &dataStream1,&frChanIn1, frStream1);
  CHECKSTATUSPTR (status);
  if (strcmp(params->frameCache1, params->frameCache2) == 0)
   {
  

  /* read in second channel */
  LALFrSeek(status->statusPtr, &(bufferStartTime), frStream1);
  CHECKSTATUSPTR (status);
	    
  LALFrGetREAL4TimeSeries(status->statusPtr, &dataStream2, 
                          &frChanIn2, frStream1);
  CHECKSTATUSPTR (status);
		


  /* close frame cache */
  LALFrClose(status->statusPtr, &frStream1);
  CHECKSTATUSPTR (status);
		
   }
  else
   {
    

    /* close first frame cache */
    LALFrClose(status->statusPtr, &frStream1);
    CHECKSTATUSPTR (status);
    

    /* open second frame cache and read in second channel */
    LALFrCacheImport(status->statusPtr, &frCache2, params->frameCache2);
    CHECKSTATUSPTR (status);
    LALFrCacheOpen(status->statusPtr, &frStream2, frCache2);
    CHECKSTATUSPTR (status);
    

    /* read in second channel */
    LALFrSeek(status->statusPtr, &(bufferStartTime), frStream2);
    CHECKSTATUSPTR (status);		
    LALFrGetREAL4TimeSeries(status->statusPtr, &dataStream2,
                            &frChanIn2, frStream2);
    CHECKSTATUSPTR (status);	
    

    /* close second frame stream */
    LALFrClose(status->statusPtr, &frStream2);
    CHECKSTATUSPTR (status);
		
   }

  /* resample */
  if (resampleRate != sampleRate)
   {
  
   /* set resample parameters */
   resampleParams.deltaT = 1.0 / (REAL8)resampleRate;
   resampleParams.filterType = defaultButterworth;

   /* resample */
   LALResampleREAL4TimeSeries(status->statusPtr, &dataStream1,&resampleParams);
   CHECKSTATUSPTR (status);
   LALResampleREAL4TimeSeries(status->statusPtr, &dataStream2,&resampleParams);
   CHECKSTATUSPTR (status);
		
  }

 /* build output */
 strncpy(streamPair->stream1->name,dataStream1.name, LALNameLength);
 strncpy(streamPair->stream2->name,dataStream2.name, LALNameLength);
 streamPair->stream1->epoch.gpsSeconds = startTime;
 streamPair->stream2->epoch.gpsSeconds = startTime;
 streamPair->stream1->epoch.gpsNanoSeconds = 0;
 streamPair->stream2->epoch.gpsNanoSeconds = 0;
 streamPair->stream1->deltaT = 1./(REAL8)resampleRate;
 streamPair->stream2->deltaT = 1./(REAL8)resampleRate;
 streamPair->stream1->f0 = streamPair->stream2->f0 = 0;
 streamPair->stream1->sampleUnits = dataStream1.sampleUnits;
 streamPair->stream2->sampleUnits = dataStream2.sampleUnits;

 /* remove buffer, and hence corruption due to resampling */
 for (i = 0; i < params->duration * resampleRate; i++)
  {
   streamPair->stream1->data->data[i] = 
     dataStream1.data->data[i + (resampleRate * buffer)];
   streamPair->stream2->data->data[i] = 
     dataStream2.data->data[i + (resampleRate * buffer)];
  }

 /* clean up */
 LALSDestroyVector(status->statusPtr, &(dataStream1.data));
 CHECKSTATUSPTR (status);
 LALSDestroyVector(status->statusPtr, &(dataStream2.data));
 CHECKSTATUSPTR (status);
	
 /* return status */
 DETATCHSTATUSPTR( status );
 RETURN( status );
}


