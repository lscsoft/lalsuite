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
#include <lal/PrintFTSeries.h>
#include <lal/Random.h>
#include <lal/IIRFilter.h>
#include "LALWrapperInterface.h"

INT4 lalDebugLevel = LALMSGLVL3;

#include <config.h>
#ifndef HAVE_LIBLALFRAME
int main( void )
{
  fputs( "Disabled: LALApps compiled with non-frame-enabled LAL\n", stderr );
  return 77;
}
#else

#define POWERC_NARGS  17

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

#define TRUE       1
#define FALSE      0

/* Usage format string. */
#define USAGE "Usage: %s \n"

/******** <lalVerbatim file="TFTilesToBurstEventsCP"> ********/
void
LALPrintBurstEvent (
               BurstEvent                           *burstEvent
               )
/******** </lalVerbatim> ********/
{
    if ( burstEvent->confidence < 0.5 ){
        fprintf(stdout,"%i %09i %f %f %f %f %f %e\n",
                burstEvent->startTime,
                burstEvent->startTimeNS,
                burstEvent->duration,
                burstEvent->centralFrequency,
                burstEvent->bandwidth,
                burstEvent->amplitude,
                burstEvent->excessPower,
                burstEvent->confidence
               );
    }
}

int main( int argc, char *argv[])
{
    static LALStatus      stat;
    INT4                  inarg        = 1;
    void                 *searchParams = NULL;
    LALInitSearchParams   initSearchParams;
    EPSearchParams       *params       = NULL;
    FrStream             *stream       = NULL;
    FrChanIn              channelIn;
    INT4                  numPoints    = 4096;
    CHAR                 *dirname      = NULL;
    REAL4TimeSeries       series;
    LIGOTimeGPS           epoch;
    BOOLEAN               epochSet     = FALSE;
    LALSearchInput        inout;
    BurstEvent           *burstEvent   = NULL;
    LALMPIParams         *mpiParams    = NULL;
    
    /* new application instance variables added 7/17/02 MSW */
    BOOLEAN               whiteNoise   = FALSE;   /* enable insertion of Gaussian white noise */
    BOOLEAN               sineBurst    = FALSE;   /* enable insertion of shaped sine burst */
    REAL4                 sineFreq     = 100.0;   /* nominal frequency of sine burst */
    REAL4                 sineOffset   = 1.0;     /* position of sine burst center in time series */
    REAL4                 sineAmpl     = 1.0;     /* peak amplitude of sine burst */
    REAL4                 sineWidth    = 1.0;     /* width (in sigmas) of sine burst */
    INT4                  index        = 0;       /* global loop index */
    BOOLEAN               applyFilter  = FALSE;   /* enable IIR filter */
    EPMethod              searchMethod = useMean; /* control search method */
    REAL4                 noiseAmpl    = 1.0;     /* gain factor for white noise */

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
        else if ( !strcmp( argv[inarg], "--noise" ) )
        {
          /* enable insertion of Gaussian white noise */
          /* up to 1 numeric parameter may be entered */
          whiteNoise = TRUE;
          inarg++;
          if (argv[inarg] && (argv[inarg][0] != '-'))
          {
            noiseAmpl = atof(argv[inarg++]);
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
        /* enum to control search method added 7/19/02 MSW */
        else if ( !strcmp( argv[inarg], "--mean" ) )
        {
          searchMethod = useMean;
          inarg++;
        }
        else if ( !strcmp( argv[inarg], "--median" ) )
        {
          searchMethod = useMedian;
          inarg++;
        }
        else if ( !strcmp( argv[inarg], "--unity" ) )
        {
          searchMethod = useUnity;
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
    initSearchParams.nodeClass = 0;
    initSearchParams.startTime = 0;
    initSearchParams.dataDuration = 0;  
    initSearchParams.realtimeRatio = 0;
    series.data = NULL;
    LALCreateVector( &stat, &series.data, numPoints);
    for (index = 0; index < numPoints; index++)
      {series.data->data[index] = 0.0;}  /* clear out vector */
    if ( epochSet )
    {
      series.epoch.gpsSeconds     = epoch.gpsSeconds;
      series.epoch.gpsNanoSeconds = epoch.gpsNanoSeconds;
    }
    strcpy(series.name, "what");
    series.deltaT = 61.03515625e-6;
    series.f0 = 0.0;
    series.sampleUnits.powerOfTen = 0;
    for (index = 0; index < LALNumUnits; index++)
    {
      series.sampleUnits.unitNumerator[index] = 0;
      series.sampleUnits.unitDenominatorMinusOne[index] = 0;
    }
    
    LALInitSearch( &stat, &searchParams, &initSearchParams);
    params = (EPSearchParams *) searchParams;
    
    /* set method in search parameters */
    params->initParams->method = searchMethod;

    /*******************************************************************
    * GET AND CONDITION THE DATA                                       *
    *******************************************************************/
    
    /* only try to load frame if name is specified */
    if (dirname)
    {
      /* Open frame stream */
      LALFrOpen( &stat, &stream, dirname, "*.gwf" );

      /* Determine information about the channel */
      LALFrGetREAL4TimeSeries( &stat, &series, &channelIn, stream);
      LALFrSeek(&stat, &(series.epoch), stream);

      /* get the data */
      LALFrGetREAL4TimeSeries( &stat, &series, &channelIn, stream);

      /* close the frame stream */
      LALFrClose( &stat, &stream );
    }

    /* populate time series with white noise if specified */
    if (whiteNoise)
    {
      long   length  = series.data->length;
      long   index   = 0;
      static RandomParams *rparams = NULL;
      INT4 seed = 1;  /* TODO set non-zero for debugging 6/11/02 MSW */

      /* generate Gaussian white noise with unit variance */
      LALCreateRandomParams(&stat, &rparams, seed);
      LALNormalDeviates(&stat, series.data, rparams);
      LALDestroyRandomParams(&stat, &rparams);
      
      /* apply gain factor to noise in time series */
      for (index = 0; index < length; index++)
      {
        double sample = series.data->data[index];
        sample *= noiseAmpl;
        series.data->data[index] = sample;
      }
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
    LALPrintTimeSeries( &series, "./timeseries.dat" );
    
    /* set up the data input for ConditionData */
    inout.numberSequences = 3;
    inout.sequences = (multiDimData *) LALMalloc( (size_t) 
            (inout.numberSequences * sizeof(multiDimData) ) );
    strcpy(inout.sequences[0].name, "channel");
    inout.sequences[0].data.real4 = series.data->data;
    inout.sequences[0].range.dTime.startSec = series.epoch.gpsSeconds;
    inout.sequences[0].range.dTime.startNan = series.epoch.gpsNanoSeconds;
    inout.sequences[0].numberDimensions = 1;
    inout.sequences[0].dimensions = (UINT4 *)  LALMalloc( (size_t) 
            (inout.sequences[0].numberDimensions * sizeof(UINT4) ) );
    inout.sequences[0].dimensions[0] = series.data->length;
    inout.sequences[0].range.dTime.timeStepSize = series.deltaT;

    /* set up the spectrum input for ConditionData */
    inout.sequences[1].range.dFreq.startFreq = 0.0;
    strcpy(inout.sequences[1].name, "spectrum");
    inout.sequences[1].range.dFreq.numberSamples
             = (params->initParams->numPoints/2 + 1);
    inout.sequences[1].range.dFreq.freqStepSize = 1.0 / 
        (series.deltaT * params->initParams->numPoints);
    inout.sequences[1].range.dTime.startSec = series.epoch.gpsSeconds;
    inout.sequences[1].range.dTime.startNan = series.epoch.gpsNanoSeconds;
    inout.sequences[1].data.real4 = (REAL4 *) LALCalloc(
            inout.sequences[1].range.dFreq.numberSamples , sizeof(REAL4) );

    /* set up the response input for ConditionData */
    strcpy(inout.sequences[2].name, "response");
    inout.sequences[2].range.dFreq.numberSamples
             = (params->initParams->numPoints/2 + 1);
    inout.sequences[2].range.dFreq.freqStepSize = 1.0 / 
        (series.deltaT * params->initParams->numPoints);
    inout.sequences[2].range.dTime.startSec = series.epoch.gpsSeconds;
    inout.sequences[2].range.dTime.startNan = series.epoch.gpsNanoSeconds;
    inout.sequences[2].data.complex8 = (COMPLEX8 *) LALCalloc(
            inout.sequences[2].range.dFreq.numberSamples , sizeof(COMPLEX8) );

    /* Finally call condition data */
    LALConditionData( &stat, &inout, searchParams, mpiParams);


    /*******************************************************************
    * DO THE SEARCH                                                    *
    *******************************************************************/
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
        fprintf(stderr,"Analyzing segments %i -- %i", params->currentSegment,
          params->currentSegment + tmpDutyCycle - 1);
        switch(searchMethod)
        {
          case useMean:
            fprintf(stderr, ", using mean method.\n");
            break;
            
          case useMedian:
            fprintf(stderr, ", using median method.\n");
            break;
            
          case useUnity:
            fprintf(stderr, ", using unity method.\n");
            break;
            
          default:
            fprintf(stderr, ", using unknown method.\n");
            break;
        }

        /* This is the main engine of the excess power method */ 
        EPSearch (&stat, params, &burstEvent, tmpDutyCycle);

        /* free the temporary memory containing the events */
        while (burstEvent)
        {
            BurstEvent *thisEvent;
            thisEvent = burstEvent;
            burstEvent = burstEvent->nextEvent;
            LALPrintBurstEvent(thisEvent);
            LALFree( thisEvent );
        }

        /* increment to the next segment number to be analyzed */
        params->currentSegment += tmpDutyCycle;
    }


    /*******************************************************************
    * FINALIZE EVERYTHING                                            *
    *******************************************************************/
    LALFree(inout.sequences[0].dimensions);
    LALFree(inout.sequences[1].data.real4);
    LALFree(inout.sequences[2].data.complex8);
    LALFree(inout.sequences);
    LALSDestroyVector( &stat, &(series.data) );
    LALFinalizeSearch( &stat, &searchParams);
    LALFree(initSearchParams.argv);
    LALCheckMemoryLeaks();

    return 0;
}

#endif
