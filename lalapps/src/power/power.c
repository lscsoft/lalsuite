#include <stdio.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/AVFactories.h>
#include <lal/LALDatatypes.h>
#include <lal/LALConstants.h>
#include <lal/FrameStream.h>
#include <lal/EPData.h>
#include <lal/EPSearch.h>
#include <lal/BurstSearch.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXML.h>
#include <lal/PrintFTSeries.h>
#include <lal/Random.h>
#include <lal/IIRFilter.h>
#include <lalapps.h>

typedef struct
{
   INT4    argc;
   CHAR**  argv;
}LALInitSearchParams;

NRCSID( POWERC, "power $Id$");
RCSID( "power $Id$");

#include <config.h>
#ifndef HAVE_LIBLALFRAME
int main( void )
{
  fputs( "Disabled: LALApps compiled with non-frame-enabled LAL\n", stderr );
  return 77;
}
#else

#define POWERC_NARGS 19 

#define POWERC_ENORM  0
#define POWERC_ESUB   1
#define POWERC_EARG   2
#define POWERC_EVAL   3
#define POWERC_EFILE  4
#define POWERC_EINPUT 5
#define POWERC_EMEM   6

#define POWERC_MSGENORM  "Normal exit"
#define POWERC_MSGESUB   "Subroutine failed"
#define POWERC_MSGEARG   "Error parsing arguments"
#define POWERC_MSGEVAL   "Input argument out of valid range"
#define POWERC_MSGEFILE  "Could not open file"
#define POWERC_MSGEINPUT "Error reading file"
#define POWERC_MSGEMEM   "Out of memory"

#define POWERC_ARGVNPTS     "--npts"
#define POWERC_ARGNPTS      1
#define POWERC_ARGVNSEG     "--nseg"
#define POWERC_ARGNSEG      2
#define POWERC_ARGVOLAP     "--olap"
#define POWERC_ARGOLAP      3
#define POWERC_ARGVOLAPFCTR "--olapfctr"
#define POWERC_ARGOLAPFCTR  4
#define POWERC_ARGVMINFBIN  "--minfbin"
#define POWERC_ARGMINFBIN   5
#define POWERC_ARGVMINTBIN  "--mintbin"
#define POWERC_ARGMINTBIN   6
#define POWERC_ARGVFLOW     "--flow"
#define POWERC_ARGFLOW      7
#define POWERC_ARGVDELF     "--delf"
#define POWERC_ARGDELF      8
#define POWERC_ARGVLNGTH    "--lngth"
#define POWERC_ARGLNGTH     9
#define POWERC_ARGVNSIGMA   "--nsigma"
#define POWERC_ARGNSIGMA    10
#define POWERC_ARGVALPHDEF  "--alphdef"
#define POWERC_ARGALPHDEF   11
#define POWERC_ARGVSEGDCLE  "--segdcle"
#define POWERC_ARGSEGDCLE   12
#define POWERC_ARGVALPHTH   "--threshold"
#define POWERC_ARGALPHTH    13
#define POWERC_ARGVE2MSTR   "--etomstr"
#define POWERC_ARGE2MSTR    14
#define POWERC_ARGVCHANNEL  "--channel"
#define POWERC_ARGCHANNEL   15
#define POWERC_ARGVSIMTYPE  "--simtype"
#define POWERC_ARGSIMTYPE   16
#define POWERC_ARGVSPECTYPE "--spectype"
#define POWERC_ARGSPECTYPE  17
#define POWERC_ARGVWINTYPE  "--window"
#define POWERC_ARGWINTYPE   18
#define POWERC_ARGVHELP  "--help"

#define TRUE       1
#define FALSE      0

/* Usage format string. */
#define USAGE "Usage: %s --npts npoints --nseg nsegments --olap overlap \
    --olapfctr olapfactor --minfbin nfbin --mintbin ntbin --flow flow \
    --delf df --lngth bandwidth --nsigma sigma --alphdef alpha \
    --segdcle nsegs --threshold threshold --etomstr nevents \
    --framedir framedir --channel channel --simtype simtype \
    --spectype spectype --window window --epoch sec nsec --numpts npts \
    [--outfile filename] [--framecache filename] [--printSpectrum] [--help]\n"

#define RESPONSE_REGEX "CAL-RESPONSE"
#define CAV_GAIN_REGEX "CAL-CAV_GAIN"
#define CAV_FAC_REGEX "CAL-CAV_FAC"
#define OLOOP_FAC_REGEX "CAL-OLOOP_FAC"

/******** <lalVerbatim file="TFTilesToBurstEventsCP"> ********/
void
LALPrintBurstEvent (
               FILE             *fpout,
               SnglBurstTable   *burstEvent
               )
/******** </lalVerbatim> ********/
{
    fprintf(fpout,"%i %09i %f %f %f %f %f %e\n",
            burstEvent->start_time.gpsSeconds,
            burstEvent->start_time.gpsNanoSeconds,
            burstEvent->duration,
            burstEvent->central_freq,
            burstEvent->bandwidth,
            burstEvent->amplitude,
            burstEvent->snr,
            burstEvent->confidence
           );
    fflush(fpout);
}

/* Fill an array with white noise */
static void makeWhiteNoise(
        LALStatus             *status,
        REAL4TimeSeries       *series, 
        REAL4                  noiseAmpl,
        INT4                   seed
        )
{
    long   length  = series->data->length;
    long   index   = 0;
    static RandomParams *rparams = NULL;

    INITSTATUS (status, "makeWhiteNoise", POWERC);
    ATTATCHSTATUSPTR (status);

    /* generate Gaussian white noise with unit variance */
    LALCreateRandomParams(status->statusPtr, &rparams, seed);
    CHECKSTATUSPTR (status);
    LALNormalDeviates(status->statusPtr, series->data, rparams);
    CHECKSTATUSPTR (status);
    LALDestroyRandomParams(status->statusPtr, &rparams);
    CHECKSTATUSPTR (status);

    /* apply gain factor to noise in time series */
    for (index = 0; index < length; index++)
    {
        double sample = series->data->data[index];
        sample *= noiseAmpl;
        series->data->data[index] = sample;
    }

    DETATCHSTATUSPTR (status);
    RETURN (status);
}




int main( int argc, char *argv[])
{
    static LALStatus      stat;
    INT4                  inarg         = 1;
    void                 *searchParams  = NULL;
    LALInitSearchParams   initSearchParams;
    EPSearchParams       *params        = NULL;
    FrStream             *stream        = NULL;
    FrCache              *frameCache    = NULL;
    FrChanIn              channelIn;
    INT4                  numPoints     = 4096;
    CHAR                 *dirname       = NULL;
    CHAR                 *cachefile     = NULL;
    CHAR                 *outputfile    = "tmp.out";
    REAL4TimeSeries       series;
    LIGOTimeGPS           epoch;
    BOOLEAN               epochSet      = FALSE;
    BOOLEAN               printSpectrum = FALSE;
    SnglBurstTable       *burstEvent    = NULL;
    LIGOLwXMLStream     xmlStream;
    MetadataTable        myTable;

    /* new application instance variables added 7/17/02 MSW */
    BOOLEAN               whiteNoise    = FALSE;   /* insertion of Gaussian white noise */
    BOOLEAN               sineBurst     = FALSE;   /* insertion of shaped sine burst */
    REAL4                 sineFreq      = 100.0;   /* nominal frequency of sine burst */
    REAL4                 sineOffset    = 1.0;     /* sin-burst center in time series */
    REAL4                 sineAmpl      = 1.0;     /* peak amplitude of sine burst */
    REAL4                 sineWidth     = 1.0;     /* width (in sigmas) of sine burst */
    INT4                  index         = 0;       /* global loop index */
    BOOLEAN               applyFilter   = FALSE;   /* enable IIR filter */
    REAL4                 noiseAmpl     = 1.0;     /* gain factor for white noise */
    INT4                  seed = 1;  /* TODO set non-zero for debugging 6/11/02 MSW */
    FILE                 *fpout         = NULL;    /* output file */

    initSearchParams.argv = (CHAR **) LALMalloc( (size_t)(POWERC_NARGS * sizeof(CHAR*)) );

    /*******************************************************************
    * PARSE ARGUMENTS (inarg stores the current position)               *
    *******************************************************************/

    if (argc <= 1){
        LALPrintError( USAGE, *argv );
        return 0;
    }
    
    initSearchParams.argc = 1;
    initSearchParams.argv[0] = "-filterparams";

    while ( inarg < argc ) {
        /* Parse output file option. */
        if ( !strcmp( argv[inarg], POWERC_ARGVNPTS ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                initSearchParams.argv[POWERC_ARGNPTS] = argv[inarg++];
                initSearchParams.argc++;
            }else{
                LALPrintError( USAGE, *argv );
                return POWERC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], POWERC_ARGVHELP ) ) {
                LALPrintError( USAGE, *argv );
                return POWERC_EARG;
        }
        else if ( !strcmp( argv[inarg], POWERC_ARGVNSEG ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                initSearchParams.argv[POWERC_ARGNSEG] = argv[inarg++];
                initSearchParams.argc++;
            }else{
                LALPrintError( USAGE, *argv );
                return POWERC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], POWERC_ARGVOLAP ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                initSearchParams.argv[POWERC_ARGOLAP] = argv[inarg++];
                initSearchParams.argc++;
            }else{
                LALPrintError( USAGE, *argv );
                return POWERC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], POWERC_ARGVOLAPFCTR ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                initSearchParams.argv[POWERC_ARGOLAPFCTR] = argv[inarg++];
                initSearchParams.argc++;
            }else{
                LALPrintError( USAGE, *argv );
                return POWERC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], POWERC_ARGVMINFBIN ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                initSearchParams.argv[POWERC_ARGMINFBIN] = argv[inarg++];
                initSearchParams.argc++;
            }else{
                LALPrintError( USAGE, *argv );
                return POWERC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], POWERC_ARGVMINTBIN ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                initSearchParams.argv[POWERC_ARGMINTBIN] = argv[inarg++];
                initSearchParams.argc++;
            }else{
                LALPrintError( USAGE, *argv );
                return POWERC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], POWERC_ARGVFLOW ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                initSearchParams.argv[POWERC_ARGFLOW] = argv[inarg++];
                initSearchParams.argc++;
            }else{
                LALPrintError( USAGE, *argv );
                return POWERC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], POWERC_ARGVDELF ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                initSearchParams.argv[POWERC_ARGDELF] = argv[inarg++];
                initSearchParams.argc++;
            }else{
                LALPrintError( USAGE, *argv );
                return POWERC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], POWERC_ARGVLNGTH ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                initSearchParams.argv[POWERC_ARGLNGTH] = argv[inarg++];
                initSearchParams.argc++;
            }else{
                LALPrintError( USAGE, *argv );
                return POWERC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], POWERC_ARGVNSIGMA ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                initSearchParams.argv[POWERC_ARGNSIGMA] = argv[inarg++];
                initSearchParams.argc++;
            }else{
                LALPrintError( USAGE, *argv );
                return POWERC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], POWERC_ARGVALPHDEF ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                initSearchParams.argv[POWERC_ARGALPHDEF] = argv[inarg++];
                initSearchParams.argc++;
            }else{
                LALPrintError( USAGE, *argv );
                return POWERC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], POWERC_ARGVSEGDCLE ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                initSearchParams.argv[POWERC_ARGSEGDCLE] = argv[inarg++];
                initSearchParams.argc++;
            }else{
                LALPrintError( USAGE, *argv );
                return POWERC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], POWERC_ARGVALPHTH ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                initSearchParams.argv[POWERC_ARGALPHTH] = argv[inarg++];
                initSearchParams.argc++;
            }else{
                LALPrintError( USAGE, *argv );
                return POWERC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], POWERC_ARGVE2MSTR ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                initSearchParams.argv[POWERC_ARGE2MSTR] = argv[inarg++];
                initSearchParams.argc++;
            }else{
                LALPrintError( USAGE, *argv );
                return POWERC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], POWERC_ARGVCHANNEL ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                initSearchParams.argv[POWERC_ARGCHANNEL] = argv[inarg++];
                initSearchParams.argc++;
                channelIn.name = initSearchParams.argv[POWERC_ARGCHANNEL];
                channelIn.type = ADCDataChannel;
            }else{
                LALPrintError( USAGE, *argv );
                return POWERC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], POWERC_ARGVSPECTYPE ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                initSearchParams.argv[POWERC_ARGSPECTYPE] = argv[inarg++];
                initSearchParams.argc++;
            }else{
                LALPrintError( USAGE, *argv );
                return POWERC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], POWERC_ARGVSIMTYPE ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                initSearchParams.argv[POWERC_ARGSIMTYPE] = argv[inarg++];
                initSearchParams.argc++;
            }else{
                LALPrintError( USAGE, *argv );
                return POWERC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], POWERC_ARGVWINTYPE ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                initSearchParams.argv[POWERC_ARGWINTYPE] = argv[inarg++];
                initSearchParams.argc++;
            }else{
                LALPrintError( USAGE, *argv );
                return POWERC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], "--epoch" ) ) {
            if ( argc > inarg + 2 ) {
                inarg++;
                epoch.gpsSeconds = atoi( argv[inarg++] );
                epoch.gpsNanoSeconds = atoi( argv[inarg++] );
                epochSet = TRUE;
            }else{
                LALPrintError( USAGE, *argv );
                return POWERC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], "--printSpectrum" ) ) {
            inarg++;
            printSpectrum = TRUE;
        }
        else if ( !strcmp( argv[inarg], "--numpts" ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                numPoints = atoi( argv[inarg++] );
            }else{
                LALPrintError( USAGE, *argv );
                return POWERC_EARG;
            }
        }
        else if ( !strcmp( argv[inarg], "--framedir" ) )
        {
          /* enable loading frame data from disk */
          if ( argc > inarg + 1 ) {
              inarg++;
              dirname = argv[inarg++];
          }else{
              LALPrintError( USAGE, *argv );
              return POWERC_EARG;
          }
        }
        else if ( !strcmp( argv[inarg], "--framecache" ) )
        {
          /* enable loading frame data from disk using a cache file */
          if ( argc > inarg + 1 ) {
              inarg++;
              cachefile = argv[inarg++];
          }else{
              LALPrintError( USAGE, *argv );
              return POWERC_EARG;
          }
        }
        else if ( !strcmp( argv[inarg], "--outfile" ) )
        {
          /* enable loading frame data from disk using a cache file */
          if ( argc > inarg + 1 ) {
              inarg++;
              outputfile = argv[inarg++];
          }else{
              LALPrintError( USAGE, *argv );
              return POWERC_EARG;
          }
        }
        else if ( !strcmp( argv[inarg], "--noise" ) )
        {
          /* enable insertion of Gaussian white noise */
          /* up to 1 numeric parameter may be entered */
          whiteNoise = TRUE;
          inarg++;
          if (argv[inarg] && (argv[inarg][0] != '-'))
          {
              noiseAmpl = atof(argv[inarg++]);
              if (argv[inarg] && (argv[inarg][0] != '-'))
              {
                  seed = atoi(argv[inarg++]);
              }
          }
        }
        else if ( !strcmp( argv[inarg], "--sine" ) )
        {
          /* enable insertion of shaped sine burst */
          /* up to 4 numeric parameters may be entered */
          sineBurst = TRUE;
          inarg++;
          if (argv[inarg] && (argv[inarg][0] != '-'))
          {
            sineFreq = atof(argv[inarg++]);
            if (argv[inarg] && (argv[inarg][0] != '-'))
            {
              sineOffset = atof(argv[inarg++]);
              if (argv[inarg] && (argv[inarg][0] != '-'))
              {
                sineAmpl = atof(argv[inarg++]);
                if (argv[inarg] && (argv[inarg][0] != '-'))
                {
                  sineWidth = atof(argv[inarg++]);
                }
              }
            }
          }
        }
        else if ( !strcmp( argv[inarg], "--filter" ) )
        {
          /* enable IIR filter */
          /* no numeric parameters are expected */
          applyFilter = TRUE;
          inarg++;
        }
        /* Check for unrecognized options. */
        else if ( argv[inarg][0] == '-' ) {
            LALPrintError( "Unrecognized options\n" );
            return POWERC_EARG;
        }
        /* Default case for unknown arguments */
        else
        {
          LALPrintError( "Unknown arguments\n" );
          return POWERC_EARG;
        }
    } /* End of argument parsing loop. */

    /*******************************************************************
    * INITIALIZE EVERYTHING                                            *
    *******************************************************************/

    /* Which error handler to use */
    lal_errhandler = LAL_ERR_EXIT;
    set_debug_level( "3" );

    /* Make sure all pointers are NULL */
    xmlStream.fp=NULL;

    /* create and initialize the time series vector */
    series.data = NULL;
    LAL_CALL( LALCreateVector( &stat, &series.data, numPoints), &stat);
    for (index = 0; index < numPoints; index++) {
        series.data->data[index] = 0.0;
    }  
    
    if ( epochSet )
    {
      series.epoch.gpsSeconds     = epoch.gpsSeconds;
      series.epoch.gpsNanoSeconds = epoch.gpsNanoSeconds;
    }
    strcpy(series.name, "Data");
    series.deltaT = 4.8828e-4;
    series.f0 = 0.0;
    series.sampleUnits.powerOfTen = 0;
    for (index = 0; index < LALNumUnits; index++)
    {
      series.sampleUnits.unitNumerator[index] = 0;
      series.sampleUnits.unitDenominatorMinusOne[index] = 0;
    }
    
    LAL_CALL( EPInitSearch( &stat, &searchParams, initSearchParams.argv, 
                initSearchParams.argc), &stat);
    params = (EPSearchParams *) searchParams;
    params->printSpectrum = printSpectrum;
    
    /*******************************************************************
    * GET AND CONDITION THE DATA                                       *
    *******************************************************************/
    
    /* only try to load frame if name is specified */
    if (dirname || cachefile)
    {
        if(dirname){
            /* Open frame stream */
            LAL_CALL( LALFrOpen( &stat, &stream, dirname, "*.gwf" ), &stat);
        }
        else if (cachefile){
            /* Open frame cache */
            LAL_CALL( LALFrCacheImport( &stat, &frameCache, cachefile ), &stat);
            LAL_CALL( LALFrCacheOpen( &stat, &stream, frameCache ), &stat);
            LAL_CALL( LALDestroyFrCache( &stat, &frameCache ), &stat );
        }
        /*
         * Determine information about the channel and seek to the
         * right place in the fram files 
         */
        LAL_CALL( LALFrGetREAL4TimeSeries( &stat, &series, &channelIn, stream), &stat);
        if ( epochSet )
        {
            series.epoch.gpsSeconds     = epoch.gpsSeconds;
            series.epoch.gpsNanoSeconds = epoch.gpsNanoSeconds;
        }
        LAL_CALL( LALFrSeek(&stat, &(series.epoch), stream), &stat);

        /* get the data */
        LAL_CALL( LALFrGetREAL4TimeSeries( &stat, &series, &channelIn, stream), &stat);

        /* close the frame stream */
        LAL_CALL( LALFrClose( &stat, &stream ), &stat);
    }

    /* populate time series with white noise if specified */
    if (whiteNoise)
    {
        makeWhiteNoise(&stat, &series, noiseAmpl, seed);
    }

    /* insert shaped sine burst into time series if specified */
    if (sineBurst)
    {
      double deltaT  = series.deltaT;
      long   length  = series.data->length;
      long   index   = 0;
      long   offset  = floor (sineOffset / deltaT);
      double sigmaSq = sineWidth * sineWidth / deltaT / deltaT;
      double incr    = 2.0 * 3.141593 * sineFreq * deltaT;
 
      for (index = 0; index < length; index++)
      {
        double sample = series.data->data[index];
        sample += exp(-1.0 / 2.0 * (index - offset) * (index - offset) / sigmaSq)
               * sineAmpl * sin(index * incr);
        series.data->data[index] = sample;
      }
    }
    
    /* apply IIR filter to input time series */
    if (applyFilter)
    {
      /* IIR filter coefficients obtained from MatLab */
      /* numerator coefficients are 'direct', denominator coefficients are 'recursive' */
      /* NOTE: sign change required for all recursive coefficients from MatLab */
      /* coefficients shown are for 4th order band reject at 1 kHz */
#define FILTER_ORDER 4
#define DEGRADE   0.125
      REAL4 mydirectCoef[FILTER_ORDER + 1] =
        {
          0.6640677019464 + DEGRADE * -1.0,
         -2.4578369166290 + DEGRADE *  2.9719655202650,
          3.6023622225730 + DEGRADE * -3.4861344923810,
         -2.4578369166290 + DEGRADE *  1.9437083129930,
          0.6640677019464 + DEGRADE * -0.4443631340847
        };
      REAL4 myrecursCoef[FILTER_ORDER + 1] =
        {
         -1.0000000000000,       /* this coefficient must = -1 */
          2.9719655202650,
         -3.4861344923810,
          1.9437083129930,
         -0.4443631340847
        };
      REAL4 myhistory[FILTER_ORDER] =
        {0.0, 0.0, 0.0, 0.0}; /* history length equals filter order */

      /* allocate local data structures */
      REAL4IIRFilter myfilter;
      REAL4Vector mydirectVector;
      REAL4Vector myrecursVector;
      REAL4Vector myhistoryVector;

      /* populate filter structs */
      myfilter.name = "my_filter";
      mydirectVector.length  = FILTER_ORDER + 1;  /* same as array sizes, above */
      myrecursVector.length  = FILTER_ORDER + 1;
      myhistoryVector.length = FILTER_ORDER;
      mydirectVector.data    = mydirectCoef;
      myrecursVector.data    = myrecursCoef;
      myhistoryVector.data   = myhistory;
      myfilter.deltaT        = series.deltaT;
      myfilter.directCoef    = &mydirectVector;
      myfilter.recursCoef    = &myrecursVector;
      myfilter.history       = &myhistoryVector;

      /* apply instant-by-instant filter to fill pipeline */
      for (index = 0; index < numPoints; index++)
        {LALSIIRFilter(series.data->data[index], &myfilter);}  /* discard filter output */

      /* apply IIR in place to modify input vector */
      LALIIRFilterREAL4Vector(&stat, series.data, &myfilter);
    }
    
    /* write diagnostic info to disk */
    if ( printSpectrum ){
        LALPrintTimeSeries( &series, "./timeseries.dat" );
    }

    /* Finally call condition data */
    LAL_CALL( EPConditionData( &stat, &series, searchParams), &stat);

    
    /*******************************************************************
    * DO THE SEARCH                                                    *
    *******************************************************************/
    LAL_CALL( LALOpenLIGOLwXMLFile(&stat, &xmlStream, outputfile), &stat);
    LAL_CALL( LALBeginLIGOLwXMLTable (&stat, &xmlStream, sngl_burst_table), &stat);
    while ( params->currentSegment < params->initParams->numSegments )
    {
        UINT4       tmpDutyCycle=0;
        UINT4       dumCurrentSeg=params->currentSegment;

        /* count the segments to analyze */
        for ( tmpDutyCycle=0 ; tmpDutyCycle < params->initParams->segDutyCycle && 
                dumCurrentSeg < params->initParams->numSegments ; tmpDutyCycle++ )
        {
            dumCurrentSeg++;
        }
                
        /* tell operator how we are doing */
        fprintf(stderr,"Analyzing segments %i -- %i\n", params->currentSegment,
                params->currentSegment + tmpDutyCycle - 1);

        /* This is the main engine of the excess power method */ 
        LAL_CALL( EPSearch (&stat, params, &burstEvent, tmpDutyCycle), &stat);

        /* Write the results to the burst table */
        myTable.snglBurstTable=burstEvent;
        LAL_CALL( LALWriteLIGOLwXMLTable (&stat, &xmlStream, myTable, 
                    sngl_burst_table), &stat);

        /* free the temporary memory containing the events */
        while (burstEvent)
        {
            SnglBurstTable *thisEvent;
            thisEvent = burstEvent;
            burstEvent = burstEvent->next;
            LALFree( thisEvent );
        }

        /* increment to the next segment number to be analyzed */
        params->currentSegment += tmpDutyCycle;
    }
    LAL_CALL( LALEndLIGOLwXMLTable (&stat, &xmlStream), &stat);
    LAL_CALL( LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);


    /*******************************************************************
    * FINALIZE EVERYTHING                                            *
    *******************************************************************/
    LAL_CALL( LALSDestroyVector( &stat, &(series.data) ), &stat);
    LAL_CALL( EPFinalizeSearch( &stat, &searchParams), &stat);
    LALCheckMemoryLeaks();

    return 0;
}

#endif
