/*
 * coherence.c - SGWB Standalone Analysis Pipeline
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
#include <lal/Calibration.h>
#include <lal/ComplexFFT.h>
#include <lal/CoarseGrainFrequencySeries.h>
#include <lal/Date.h>
#include <lal/DetectorSite.h>
#include <lal/FrameCache.h>
#include <lal/FrameCalibration.h>
#include <lal/FrameStream.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALRCSID.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/PrintFTSeries.h>
#include <lal/PrintVector.h>
#include <lal/Random.h>
#include <lal/RealFFT.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/StreamInput.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Units.h>
#include <lal/Window.h>
#include <lal/IIRFilter.h>
#include <lal/BandPassTimeSeries.h>
#include "stochastic.h"

NRCSID (STOCHASTICC, "$Id$");
RCSID ("$Id$");

/* cvs info */
#define CVS_ID "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "lalapps_coherence"

/* variables for getopt options parsing */
char *optarg;
int optind;

/* flags for getopt_long */

static int overlap_hann_flag = 0;
static int verbose_flag = 0;
static int condor_flag = 0;

/* parameters for the coherence analysis  */

/* sampling parameters */
INT4 sampleRate = 16384;
INT4 resampleRate = 1024;
REAL8 deltaF = 0.25;
INT4 numFFT = 10.;
/* data parameters */
LIGOTimeGPS gpsStartTime, gpsCalibTime;
UINT8 startTime = 730793098;
UINT8 stopTime = 730793138;
CHAR frameCache1[100] = "cachefiles/H-730793097.cache";
CHAR frameCache2[100] = "cachefiles/L-730793097.cache";
CHAR channel1[LALNameLength]= "H1:LSC-AS_Q";
CHAR channel2[LALNameLength]= "L1:LSC-AS_Q";
CHAR ifo1[LALNameLength] = "H1";
CHAR ifo2[LALNameLength] = "L1";

INT4 main(INT4 argc, CHAR *argv[])
 {
  /* variable declarations */

  /* status pointer */
  LALStatus status;

  /* counters */
  INT4 i,j, segLoop;

  /* input data segment */
  INT4 numSegments;
  INT4 segmentDuration;
  INT4 segmentLength;
  INT4 segmentShift;
  INT4 padData;
  ReadDataPairParams streamParams;
  StreamPair streamPair;
  REAL4TimeSeries segment1, segment2;

    
  /* data structures for PSDs */
  INT4 overlapPSDLength;
  INT4 psdLength;
  INT4 windowPSDLength;
  LALWindowParams winparPSD;
  AverageSpectrumParams specparPSD;
  REAL4FrequencySeries psd1,psd2;
  COMPLEX8FrequencySeries csd, coh, cohTot;
  LALUnit psdUnits = {0,{0,0,1,0,0,0,2},{0,0,0,0,0,0,0}};

  /* error handler */
  status.statusPtr = NULL;

  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "1" );

  /* parse command line options */
  parseOptions(argc, argv);


  /* read parameters into input parameter file */
  if (condor_flag == 1)
   { 
     fscanf(stdin,"%d\n",&startTime, &startTime);
     fscanf(stdin,"%d\n",&stopTime, &stopTime);
     fscanf(stdin,"%s\n%s\n",&frameCache1,&frameCache2);
   }
       
  if (verbose_flag)
   {fprintf(stdout, "Calculating number of segments...\n");}

  /* get number of segments */
  segmentDuration = numFFT * (1. / deltaF);
  numSegments = (INT4)((stopTime - startTime) / segmentDuration );
  segmentShift = segmentDuration / 2;
  if (overlap_hann_flag)
   {
    numSegments = 2 * numSegments - 1;
    segmentShift = segmentDuration / 2;
   }

  /* set length for data segments */
  segmentLength = segmentDuration * resampleRate;
  if (sampleRate == resampleRate)
    { padData = 0;}
  else {padData = 1;}     
  strncpy(segment1.name, "segment1", LALNameLength);
  strncpy(segment2.name, "segment2", LALNameLength);
  segment1.sampleUnits = segment2.sampleUnits = lalADCCountUnit;
  segment1.epoch = segment2.epoch = gpsStartTime;
  segment1.deltaT = segment2.deltaT = 1./(REAL8)resampleRate;
  segment1.f0 = segment2.f0 = 0;

  if (verbose_flag)
   {fprintf(stdout, "Allocating memory for data segments...\n");}

  /* allocate memory for data segments */
  

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


  /* set PSD window length */
  windowPSDLength = (UINT4)(resampleRate / deltaF);

  /* set parameters for PSD estimation */
  overlapPSDLength = windowPSDLength / 2;
  psdLength = (windowPSDLength / 2) + 1;

  /* set metadata fields for PSDs */
  strncpy(psd1.name, "psd1", LALNameLength);
  strncpy(psd2.name, "psd2", LALNameLength);
  psd1.sampleUnits = psd2.sampleUnits = psdUnits;
  psd1.deltaF = psd2.deltaF = deltaF;
  psd1.f0 = psd2.f0 = 0;

  /* set metadata fields for CSD */
  strncpy(csd.name, "csd", LALNameLength);
  csd.sampleUnits = psdUnits;
  csd.deltaF  = deltaF;
  csd.f0 = 0;

  /* set metadata fields for coherence */
  strncpy(coh.name, "coh", LALNameLength);
  strncpy(cohTot.name, "cohTot", LALNameLength);
  coh.deltaF = cohTot.deltaF = deltaF;
  coh.f0 = cohTot.f0 = 0;


  if (verbose_flag)
   { fprintf(stdout, "Allocating memory for PSDs...\n");}

  /* allocate memory for PSDs */

  psd1.data = psd2.data = NULL;
  LAL_CALL( LALCreateVector(&status, &(psd1.data), psdLength),
	    &status );
  LAL_CALL( LALCreateVector(&status, &(psd2.data), psdLength), 
	    &status );
  memset( psd1.data->data, 0, 
          psd1.data->length * sizeof(*psd1.data->data));
  memset( psd2.data->data, 0, 
          psd2.data->length * sizeof(*psd2.data->data));

  
  /* allocate memory for CSD  */ 
     
  csd.data = NULL;
  LAL_CALL( LALCCreateVector(&status, &(csd.data), psdLength),
	    &status );
  memset( csd.data->data, 0, 
          csd.data->length * sizeof(*csd.data->data));

  /* allocate memory for coherence */

  coh.data = cohTot.data = NULL;
  LAL_CALL( LALCCreateVector(&status, &(coh.data), psdLength),
	    &status );
  LAL_CALL( LALCCreateVector(&status, &(cohTot.data), psdLength), 
	    &status );
  memset( coh.data->data, 0, 
          coh.data->length * sizeof(*coh.data->data));
  memset( cohTot.data->data, 0, 
          cohTot.data->length * sizeof(*cohTot.data->data));
 
  

  /* set window parameters for PSD estimation */
  winparPSD.length = windowPSDLength;
  winparPSD.type = Hann;

  /* set parameters for PSD estimation */
  specparPSD.method = useMean;
  specparPSD.overlap = overlapPSDLength;
  specparPSD.plan = NULL;
  specparPSD.window = NULL;

  if (verbose_flag)
   {fprintf(stdout, "Creating FFT plan for PSD estimation...\n");}

  /* create fft plan */
  LAL_CALL ( LALCreateForwardRealFFTPlan(&status, &specparPSD.plan,
             windowPSDLength, 0), &status );

  if (verbose_flag)
   {fprintf(stdout, "Creating window for PSD estimation...\n");}

  /* create window for PSD estimation */
  LAL_CALL( LALCreateREAL4Window(&status, &specparPSD.window, &winparPSD), 
            &status );

 
   
   /** loop over segments **/
   

   if (verbose_flag)
    { fprintf(stdout, "Looping over %d segments...\n", numSegments);}
        
   lal_errhandler = LAL_ERR_RTRN;

   for (segLoop = 0; segLoop < numSegments; segLoop++)
    {
     /* define segment epoch */
     gpsStartTime.gpsSeconds = startTime + (segLoop * segmentShift);
     segment1.epoch = segment2.epoch = gpsStartTime;
    
     if (verbose_flag)
      {
       fprintf( stdout, "Performing search on segment %d of %d...\n", 
                segLoop + 1, numSegments);
      }

     /* read data and downsample */
     if (verbose_flag)
       { fprintf(stdout, "Reading data...\n");}

     /* read data */
     streamParams.startTime = gpsStartTime.gpsSeconds;
     LAL_CALL(readDataPair(&status, &streamPair, &streamParams), &status);
       
     /* skip segment if data not found or corrupted with 0 values */           
     if ((status.statusCode !=0)||
         (segment1.data==NULL)||(segment2.data==NULL))
      {
       clear_status(&status);
       if (segLoop < (numSegments - 1)) continue; 
       else break;   
      }
               
     /* save */
     if (verbose_flag)
      {
       LALSPrintTimeSeries(&segment1, "segment1.dat");
       LALSPrintTimeSeries(&segment2, "segment2.dat");
      }
       
     if (verbose_flag)
      { fprintf(stdout, "Estimating PSDs...\n");}

     /* compute CSD and  PSDs */

     LAL_CALL( LALCOMPLEX8AverageSpectrum(&status, &csd, &segment1,&segment2, 
               &specparPSD), &status );
     LAL_CALL( LALREAL4AverageSpectrum(&status, &psd1, &segment1, 
               &specparPSD), &status );
     LAL_CALL( LALREAL4AverageSpectrum(&status, &psd2, &segment2, 
               &specparPSD), &status );


    /* output the results */
    if (verbose_flag)
     {
      LALCPrintFrequencySeries(&csd, "csd.dat");
      LALSPrintFrequencySeries(&psd1, "psd1.dat");
      LALSPrintFrequencySeries(&psd2, "psd2.dat");
     }
 
    /* compute coherence */
    for (i = 0; i < psdLength; i ++)
     {
      coh.data->data[i].re = csd.data->data[i].re / sqrt(psd1.data->data[i] 
                             * psd2.data->data[i]);  
      coh.data->data[i].im = csd.data->data[i].im / sqrt(psd1.data->data[i] 
                             * psd2.data->data[i]);
  
     }
 
   /* sum over segments */
   for (i = 0; i < psdLength; i ++)
    {
     cohTot.data->data[i].re = cohTot.data->data[i].re + coh.data->data[i].re;
     cohTot.data->data[i].im = cohTot.data->data[i].im + coh.data->data[i].im; 
    }
	 
   }

   /* output the results */
   fprintf(stdout,"%d\n", j);
   for ( i = 0; i < psdLength; i++ )
    {
     fprintf(stdout,"%e\t%e\t%e\n", 
      i * cohTot.deltaF,cohTot.data->data[i].re,cohTot.data->data[i].im );
     }         

       
   lal_errhandler = LAL_ERR_EXIT;

   /* cleanup */

   LAL_CALL( LALDestroyRealFFTPlan(&status, &(specparPSD.plan)), &status );
   LAL_CALL( LALDestroyVector(&status, &(segment1.data)), &status );
   LAL_CALL( LALDestroyVector(&status, &(segment2.data)), &status );
   LAL_CALL( LALDestroyVector(&status, &(psd1.data)), &status );
   LAL_CALL( LALDestroyVector(&status, &(psd2.data)), &status );
   LAL_CALL( LALCDestroyVector(&status, &(csd.data)), &status );
   LAL_CALL( LALCDestroyVector(&status, &(coh.data)), &status );
   LAL_CALL( LALCDestroyVector(&status, &(cohTot.data)), &status );



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
      {"overlap-hann", no_argument, &overlap_hann_flag, 1},
      {"condor", no_argument, &condor_flag,1},
      {"verbose", no_argument, &verbose_flag, 1},
      /* options that don't set a flag */
      {"help", no_argument, 0, 'h'},
      {"gps-start-time", required_argument, 0, 't'},
      {"gps-end-time", required_argument, 0, 'T'},
      {"sample-rate", required_argument, 0, 'A'},
      {"resample-rate", required_argument, 0, 'a'},
      {"deltaF", required_argument, 0, 'f'},
      {"number-FFT", required_argument, 0, 'n'},
      {"ifo-one", required_argument, 0, 'i'},
      {"ifo-two", required_argument, 0, 'I'},
      {"debug-level", required_argument, 0, 'z'},
      {"version", no_argument, 0, 'V'},
      {0, 0, 0, 0}
     };

    /* getopt_long stores the option here */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long(argc, argv, 
                    "ht:T:l:A:a:f:n:i:I:z:V",
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

      case 'A':
               /* sample rate */
               sampleRate = atoi(optarg);
               break;

      case 'a':
	       /* resampling */
	       resampleRate = atoi(optarg);
	       break;

      case 'f':
	       /* frequency resolution */
	       deltaF = atof(optarg);
	       break;                            

      case 'n':
               /* number of FFT */
               deltaF = atoi(optarg);
               break;
                          
      case 'i':
	       /* ifo for first stream */
	       strncpy(ifo1, optarg, LALNameLength);

	       /* set site and channel */
	       if (strncmp(ifo1, "H1", 2) == 0)
		{
		 strncpy(channel1, "H1:LSC-AS_Q", LALNameLength);
		}
	       else if (strncmp(ifo1, "H2", 2) == 0)
		{
		 strncpy(channel1, "H2:LSC-AS_Q", LALNameLength);
		}
	       else if (strncmp(ifo1, "L1", 2) == 0)
		{
		
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

		/* set channel */
		if (strncmp(ifo2, "H1", 2) == 0)
		 {
		  strncpy(channel2, "H1:LSC-AS_Q", LALNameLength);
		 }
		else if (strncmp(ifo2, "H2", 2) == 0)
		 {
		  
		  strncpy(channel2, "H2:LSC-AS_Q", LALNameLength);
		 }
		else if (strncmp(ifo2, "L1", 2) == 0)
		 {
		 
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
  fprintf(stderr, " -A                    sample rate\n");
  fprintf(stderr, " -a                    resample rate\n");
  fprintf(stderr, " -f                    frequency resolution\n");
  fprintf(stderr, " -n                    number of FFTs\n");      
  fprintf(stderr, " --overlap-hann        use overlap window\n");             
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
  INITSTATUS( status, "readDataPair",COHERENCEC );
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

  if (verbose_flag)
   { fprintf(stdout, "Allocating memory for raw data streams...\n");}

  /* allocate memory */
  dataStream1.data = dataStream2.data = NULL;
  LALSCreateVector( status->statusPtr, &(dataStream1.data),
                    sampleRate * (params->duration + (2 * buffer)));
	
  LALSCreateVector( status->statusPtr, &(dataStream2.data), 
                    sampleRate * (params->duration + (2 * buffer)));
	
  memset( dataStream1.data->data, 0, 
          dataStream1.data->length * sizeof(*dataStream1.data->data));
  memset( dataStream2.data->data, 0, 
          dataStream2.data->length * sizeof(*dataStream2.data->data));

  if (verbose_flag)
   { fprintf(stdout, "Opening first frame cache...\n");}

  /* open first frame cache */
  LALFrCacheImport(status->statusPtr, &frCache1, params->frameCache1);
  LALFrCacheOpen(status->statusPtr, &frStream1, frCache1);
	
  if (verbose_flag)
   { fprintf(stdout, "Reading in channel \"%s\"...\n", frChanIn1.name);}

  /* read first channel */
  LALFrSeek(status->statusPtr, &(bufferStartTime), frStream1);
  LALFrGetREAL4TimeSeries(status->statusPtr,
                          &dataStream1,&frChanIn1, frStream1);

  if (strcmp(params->frameCache1, params->frameCache2) == 0)
   {
    if (verbose_flag)
     { fprintf(stdout, "Reading in channel \"%s\" from same cache...\n",
               frChanIn2.name);
     }

  /* read in second channel */
  LALFrSeek(status->statusPtr, &(bufferStartTime), frStream1);
	    
  LALFrGetREAL4TimeSeries(status->statusPtr, &dataStream2, 
                          &frChanIn2, frStream1);
		
  if (verbose_flag)
   { fprintf(stdout, "Closing frame cache...\n");}

  /* close frame cache */
  LALFrClose(status->statusPtr, &frStream1);
		
   }
  else
   {
    if (verbose_flag)
     { fprintf(stdout, "Closing first frame cache...\n");}

    /* close first frame cache */
    LALFrClose(status->statusPtr, &frStream1);
    if (verbose_flag)
     { fprintf(stdout, "Opening second frame cache...\n");}

    /* open second frame cache and read in second channel */
    LALFrCacheImport(status->statusPtr, &frCache2, params->frameCache2);
    LALFrCacheOpen(status->statusPtr, &frStream2, frCache2);
    if (verbose_flag)
     { fprintf(stdout, "Reading in channel \"%s\"...\n", frChanIn2.name);}

    /* read in second channel */
    LALFrSeek(status->statusPtr, &(bufferStartTime), frStream2);
		
    LALFrGetREAL4TimeSeries(status->statusPtr, &dataStream2,
                            &frChanIn2, frStream2);
		
    if (verbose_flag)
     { fprintf(stdout, "Closing second frame cache...\n");}

    /* close second frame stream */
    LALFrClose(status->statusPtr, &frStream2);
		
   }

  /* resample */
  if (resampleRate != sampleRate)
   {
    if (verbose_flag)
     { fprintf(stdout, "Resampling to %d Hz...\n", resampleRate);}

   /* set resample parameters */
   resampleParams.deltaT = 1.0 / (REAL8)resampleRate;
   resampleParams.filterType = defaultButterworth;

   /* resample */
   LALResampleREAL4TimeSeries(status->statusPtr, &dataStream1,&resampleParams);
   LALResampleREAL4TimeSeries(status->statusPtr, &dataStream2,&resampleParams);
		
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
 LALSDestroyVector(status->statusPtr, &(dataStream2.data));
	
 /* return status */
 DETATCHSTATUSPTR( status );
 RETURN( status );
}

