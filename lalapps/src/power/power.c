#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/AVFactories.h>
#include <lal/LALDatatypes.h>
#include <lal/LALConstants.h>
#include <lal/FrameStream.h>
#include <lal/EPSearch.h>
#include <lal/BurstSearch.h>
#include "LALWrapperInterface.h"

INT4 lalDebugLevel = LALMSGLVL3;

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
    INT4                  inarg = 1;
    void                 *searchParams = NULL;
    LALInitSearchParams   initSearchParams;
    EPSearchParams       *params = NULL;
    FrStream             *stream = NULL;
    FrChanIn              channelIn;
    INT4                  numPoints=4096;
    CHAR                 *dirname = NULL;
    REAL4TimeSeries       series;
    LIGOTimeGPS           epoch;
    BOOLEAN               epochSet = FALSE;
    LALSearchInput        inout;
    BurstEvent           *burstEvent = NULL;
    LALMPIParams          *mpiParams = NULL;

    initSearchParams.argv = (CHAR **) LALMalloc( (size_t)(POWERC_NARGS * sizeof(CHAR*)) );

    /*******************************************************************
    * PARSE ARGUMENTS (arg stores the current position)               *
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
        else if ( !strcmp( argv[inarg], "--framedir" ) ) {
            if ( argc > inarg + 1 ) {
                inarg++;
                dirname = argv[inarg++];
            }else{
                LALPrintError( USAGE, *argv );
                return POWERC_EARG;
            }
        }
        /* Check for unrecognized options. */
        else if ( argv[inarg][0] == '-' ) {
            LALPrintError( "Unrecognized options\n" );
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

    LALInitSearch( &stat, &searchParams, &initSearchParams);
    params = (EPSearchParams *) searchParams;

   
    /*******************************************************************
    * GET AND CONDITION THE DATA                                       *
    *******************************************************************/
    
    /* Open frame stream */
    LALFrOpen( &stat, &stream, dirname, "*.gwf" );

    /* Determine information about the channel */
    series.data = NULL;
    LALFrGetREAL4TimeSeries( &stat, &series, &channelIn, stream);
    if ( epochSet ){
        series.epoch.gpsSeconds = epoch.gpsSeconds;
        series.epoch.gpsNanoSeconds = epoch.gpsNanoSeconds;
    }
    LALFrSeek(&stat, &(series.epoch), stream);

    /* allocate time series */
    LALCreateVector( &stat, &series.data, numPoints);

    /* get the data */
    LALFrGetREAL4TimeSeries( &stat, &series, &channelIn, stream);

    /* close the frame stream */
    LALFrClose( &stat, &stream );

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
        fprintf(stderr,"Analyzing segments %i -- %i\n", params->currentSegment,
                 params->currentSegment + tmpDutyCycle);

        /* This is the main engine of the excess power method */ 
        EPSearch (&stat, params, &burstEvent, tmpDutyCycle);

        /* free the temporary memory containin the events */
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
