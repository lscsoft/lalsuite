/*----------------------------------------------------------------------- 
 * 
 * File Name: event_utils.c
 *
 * Author: Brady, P. R.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */
#include <lal/LALRCSID.h>


NRCSID (EVENT_UTILSC, "$Id$");


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/Thresholds.h>
#include <lal/Date.h>
#include "metaio.h"
#include "donald.h"
#include "event_utils.h"

float fabsf(float x);
/*******************************************************************
 *
 * This file provides the following functions:
 *      int buildVetoTimes( vetoParams *thisentry)
 *      int resolveVetoTimes( timeWindow **vwindows, vetoParams *thisentry)
 *      int buildEventList( candEvent **eventhead, timeWindow *vwindows, 
 *                      candParams candidates, int injectflag, int maxflag)
 *      int build2DHistogram(candEvent *eventhead, const char *outputfile,
 *                       int **histogram, int numbins, float minsnr, 
 *                       float maxchisq)
 *
 * TODO:   All functions should be brough up to date with LAL data
 * types.
 *      
 *
 *
 *******************************************************************/


void
buildSnglInspiralIndex(
        LALStatus             *status,
        const MetaioParseEnv   triggerEnv,
        SnglInspiralIndex        *params
        )
{

  INITSTATUS (status, "buildSnglInspiralIndex", EVENT_UTILSC);
  ATTATCHSTATUSPTR (status);

  params->ifoIndex = MetaioFindColumn( triggerEnv, "IFO");
  params->searchIndex = MetaioFindColumn( triggerEnv, "search");
  params->channelIndex = MetaioFindColumn( triggerEnv, "channel");
  params->end_timeIndex = MetaioFindColumn( triggerEnv, "end_time");
  params->end_time_nsIndex = MetaioFindColumn( triggerEnv, "end_time_ns");
  params->impulse_timeIndex = MetaioFindColumn( triggerEnv, "impulse_time");
  params->impulse_time_nsIndex = MetaioFindColumn( triggerEnv, "impulse_time_ns");
  params->template_durationIndex = MetaioFindColumn( triggerEnv, "template_duration");
  params->event_durationIndex = MetaioFindColumn( triggerEnv, "event_duration");
  params->amplitudeIndex = MetaioFindColumn( triggerEnv, "amplitude");
  params->eff_distanceIndex = MetaioFindColumn( triggerEnv, "eff_distance");
  params->coa_phaseIndex = MetaioFindColumn( triggerEnv, "coa_phase");
  params->mass1Index = MetaioFindColumn( triggerEnv, "mass1");
  params->mass2Index = MetaioFindColumn( triggerEnv, "mass2");
  params->mchirpIndex = MetaioFindColumn( triggerEnv, "mchirp");
  params->etaIndex = MetaioFindColumn( triggerEnv, "eta");
  params->tau0Index = MetaioFindColumn( triggerEnv, "tau0");
  params->tau2Index = MetaioFindColumn( triggerEnv, "tau2");
  params->tau3Index = MetaioFindColumn( triggerEnv, "tau3");
  params->tau4Index = MetaioFindColumn( triggerEnv, "tau4");
  params->tau5Index = MetaioFindColumn( triggerEnv, "tau5");
  params->ttotalIndex = MetaioFindColumn( triggerEnv, "ttotal");
  params->snrIndex = MetaioFindColumn( triggerEnv, "snr");
  params->chisqIndex = MetaioFindColumn( triggerEnv, "chisq");
  params->chisq_dofIndex = MetaioFindColumn( triggerEnv, "chisq_dof");
  params->sigmasqIndex = MetaioFindColumn( triggerEnv, "sigmasq");

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


void
getSnglInspiralEvent(
        LALStatus                *status,
        const MetaioParseEnv      triggerEnv,
        SnglInspiralTable        *inspiralEvent,
        SnglInspiralIndex        *params
        )
{

  INITSTATUS (status, "getSnglInspiralEvent", EVENT_UTILSC);
  ATTATCHSTATUSPTR (status);
  
  strcpy(inspiralEvent->ifo,
          triggerEnv->ligo_lw.table.elt[params->ifoIndex].data.lstring.data);
  strcpy(inspiralEvent->search,
          triggerEnv->ligo_lw.table.elt[params->searchIndex].data.lstring.data);
  strcpy(inspiralEvent->channel,
          triggerEnv->ligo_lw.table.elt[params->channelIndex].data.lstring.data);
  inspiralEvent->end_time.gpsSeconds = 
      triggerEnv->ligo_lw.table.elt[params->end_timeIndex].data.int_4s;  
  inspiralEvent->end_time.gpsNanoSeconds = 
      triggerEnv->ligo_lw.table.elt[params->end_time_nsIndex].data.int_4s;  
  inspiralEvent->impulse_time.gpsSeconds = 
      triggerEnv->ligo_lw.table.elt[params->impulse_timeIndex].data.int_4s;  
  inspiralEvent->impulse_time.gpsNanoSeconds = 
      triggerEnv->ligo_lw.table.elt[params->impulse_time_nsIndex].data.int_4s;  
  inspiralEvent->template_duration = 
      triggerEnv->ligo_lw.table.elt[params->template_durationIndex].data.real_8;  
  inspiralEvent->event_duration = 
      triggerEnv->ligo_lw.table.elt[params->event_durationIndex].data.real_8;  
  inspiralEvent->amplitude = 
      triggerEnv->ligo_lw.table.elt[params->amplitudeIndex].data.real_4;  
  inspiralEvent->eff_distance = 
      triggerEnv->ligo_lw.table.elt[params->eff_distanceIndex].data.real_4;  
  inspiralEvent->coa_phase = 
      triggerEnv->ligo_lw.table.elt[params->coa_phaseIndex].data.real_4;  
  inspiralEvent->mass1 = 
      triggerEnv->ligo_lw.table.elt[params->mass1Index].data.real_4;  
  inspiralEvent->mass2 = 
      triggerEnv->ligo_lw.table.elt[params->mass2Index].data.real_4;  
  inspiralEvent->mchirp = 
      triggerEnv->ligo_lw.table.elt[params->mchirpIndex].data.real_4;  
  inspiralEvent->eta = 
      triggerEnv->ligo_lw.table.elt[params->etaIndex].data.real_4;  
  inspiralEvent->tau0 = 
      triggerEnv->ligo_lw.table.elt[params->tau0Index].data.real_4;  
  inspiralEvent->tau2 = 
      triggerEnv->ligo_lw.table.elt[params->tau2Index].data.real_4;  
  inspiralEvent->tau3 = 
      triggerEnv->ligo_lw.table.elt[params->tau3Index].data.real_4;  
  inspiralEvent->tau4 = 
      triggerEnv->ligo_lw.table.elt[params->tau4Index].data.real_4;  
  inspiralEvent->tau5 = 
      triggerEnv->ligo_lw.table.elt[params->tau5Index].data.real_4;  
  inspiralEvent->ttotal = 
      triggerEnv->ligo_lw.table.elt[params->ttotalIndex].data.real_4;  
  inspiralEvent->snr = 
      triggerEnv->ligo_lw.table.elt[params->snrIndex].data.real_4;  
  inspiralEvent->chisq = 
      triggerEnv->ligo_lw.table.elt[params->chisqIndex].data.real_4;  
  inspiralEvent->chisq_dof = 
      triggerEnv->ligo_lw.table.elt[params->chisq_dofIndex].data.real_4;  
  inspiralEvent->sigmasq = 
      triggerEnv->ligo_lw.table.elt[params->sigmasqIndex].data.real_8;  
  
  DETATCHSTATUSPTR (status);
  RETURN (status);
}

/****************************************************************************
 * 
 * TEMPORARY FUNCTION TO READ IN INSPIRAL TRIGGER FILE UNTIL DUNCAN's
 * CODE GETS INTO LAL
 * 
 ***************************************************************************/
void readInspiralTriggers( 
        LALStatus *status, 
        SnglInspiralTable **eventList, 
        const CHAR *fname
        )
{
    INT4 retVal=0;
    struct MetaioParseEnvironment triggerParseEnv;
    const MetaioParseEnv triggerEnv = &triggerParseEnv;

    SnglInspiralIndex         tableIndex;
    SnglInspiralTable         inspiralEvent; 
    SnglInspiralTable        *currentEvent=NULL, *prevEvent=NULL;

    INITSTATUS (status, "readInspiralTriggers", EVENT_UTILSC);
    ATTATCHSTATUSPTR (status);

    /* open xml file at inspiral table */
    if ( (retVal = MetaioOpenTable( triggerEnv, fname, "sngl_inspiral")) !=0 ){
        fprintf(stderr, "Error opening injection file %s\n", fname );
        MetaioAbort( triggerEnv ); 
        ABORT( status, EVENTUTILSH_EFILE, EVENTUTILSH_MSGEFILE );
    }

    /* Locate the relevant columns in the inspiral table */
    buildSnglInspiralIndex(status->statusPtr, triggerEnv, &tableIndex);
    CHECKSTATUSPTR(status);

    /* Loop over the triggers */
    while (1) {

        /* get the next row */
        retVal = MetaioGetRow(triggerEnv);
        if ( retVal == -1 ) {
            printf( "Error while getting row from file %s\n",fname);
            MetaioAbort( triggerEnv ); 
            ABORT( status, EVENTUTILSH_EFILE, EVENTUTILSH_MSGEFILE );
        } else if ( retVal == 0 ) {
            /*-- Reached end of file --*/
            break;
        }

        /* get the inspiral event */
        getSnglInspiralEvent(status->statusPtr, triggerEnv, &inspiralEvent, &tableIndex);
        CHECKSTATUSPTR(status);
        
        /* allocate memory for the inspiral event */
        if ( (*eventList) == NULL )
        {
            currentEvent=(*eventList)=
                (SnglInspiralTable *) LALMalloc( sizeof(SnglInspiralTable) );
        }
        else 
        {
            currentEvent =
                (SnglInspiralTable *) LALMalloc( sizeof(SnglInspiralTable) );
        }

        /* copy event into linked list element */
        memcpy(currentEvent, &inspiralEvent, sizeof(SnglInspiralTable) );

        /* point to the next event */
        currentEvent->next = NULL;
        if (prevEvent != NULL) prevEvent->next = currentEvent;
        prevEvent = currentEvent;
        currentEvent = currentEvent->next;

    }

    MetaioAbort(triggerEnv);

    DETATCHSTATUSPTR (status);
    RETURN (status);

}



void
buildSimInspiralIndex(
        LALStatus             *status,
        const MetaioParseEnv   triggerEnv,
        SimInspiralIndex        *params
        )
{

  INITSTATUS (status, "buildSnglInspiralIndex", EVENT_UTILSC);
  ATTATCHSTATUSPTR (status);

  params->geocent_end_timeIndex = MetaioFindColumn( triggerEnv, "geocent_end_time");
  params->geocent_end_time_nsIndex = MetaioFindColumn( triggerEnv, "geocent_end_time_ns");
  params->end_time_gmstIndex = MetaioFindColumn( triggerEnv, "end_time_gmst");
  params->sourceIndex = MetaioFindColumn( triggerEnv, "source");
  params->mtotalIndex = MetaioFindColumn( triggerEnv, "mtotal");
  params->etaIndex = MetaioFindColumn( triggerEnv, "eta");
  params->distanceIndex = MetaioFindColumn( triggerEnv, "distance");
  params->longitudeIndex = MetaioFindColumn( triggerEnv, "longitude");
  params->latitudeIndex = MetaioFindColumn( triggerEnv, "latitude");
  params->inclinationIndex = MetaioFindColumn( triggerEnv, "inclination");
  params->coa_phaseIndex = MetaioFindColumn( triggerEnv, "coa_phase");
  params->polarizationIndex = MetaioFindColumn( triggerEnv, "polarization");
  params->eff_dist_hIndex = MetaioFindColumn( triggerEnv, "eff_dist_h");
  params->eff_dist_lIndex = MetaioFindColumn( triggerEnv, "eff_dist_l");
  params->simulation_idIndex = MetaioFindColumn( triggerEnv, "simulation_id");

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


void
getSimInspiralVars(
        LALStatus                *status,
        const MetaioParseEnv      triggerEnv,
        SimInspiralTable        *inspiralEvent,
        SimInspiralIndex        *params
        )
{

  INITSTATUS (status, "getSimInspiralVars", EVENT_UTILSC);
  ATTATCHSTATUSPTR (status);
  
  inspiralEvent->geocent_end_time.gpsSeconds = 
      triggerEnv->ligo_lw.table.elt[params->geocent_end_timeIndex].data.int_4s;  
  inspiralEvent->geocent_end_time.gpsNanoSeconds = 
      triggerEnv->ligo_lw.table.elt[params->geocent_end_time_nsIndex].data.int_4s;  
  inspiralEvent->end_time_gmst = 
      triggerEnv->ligo_lw.table.elt[params->end_time_gmstIndex].data.real_4;  
  strcpy(inspiralEvent->source,
          triggerEnv->ligo_lw.table.elt[params->sourceIndex].data.lstring.data);
  inspiralEvent->mtotal = 
      triggerEnv->ligo_lw.table.elt[params->mtotalIndex].data.real_4;  
  inspiralEvent->eta = 
      triggerEnv->ligo_lw.table.elt[params->etaIndex].data.real_4;  
  inspiralEvent->distance = 
      triggerEnv->ligo_lw.table.elt[params->distanceIndex].data.real_4;  
  inspiralEvent->longitude = 
      triggerEnv->ligo_lw.table.elt[params->longitudeIndex].data.real_4;  
  inspiralEvent->latitude = 
      triggerEnv->ligo_lw.table.elt[params->latitudeIndex].data.real_4;  
  inspiralEvent->inclination = 
      triggerEnv->ligo_lw.table.elt[params->inclinationIndex].data.real_4;  
  inspiralEvent->coa_phase = 
      triggerEnv->ligo_lw.table.elt[params->coa_phaseIndex].data.real_4;  
  inspiralEvent->polarization = 
      triggerEnv->ligo_lw.table.elt[params->polarizationIndex].data.real_4;  
  inspiralEvent->eff_dist_h = 
      triggerEnv->ligo_lw.table.elt[params->eff_dist_hIndex].data.real_4;  
  inspiralEvent->eff_dist_l = 
      triggerEnv->ligo_lw.table.elt[params->eff_dist_lIndex].data.real_4;  
  
  DETATCHSTATUSPTR (status);
  RETURN (status);
}

void
buildSearchSummaryIndex(
        LALStatus             *status,
        const MetaioParseEnv   triggerEnv,
        SearchSummaryIndex        *params
        )
{

  INITSTATUS (status, "buildSearchSummaryIndex", EVENT_UTILSC);
  ATTATCHSTATUSPTR (status);

  params->process_idIndex = MetaioFindColumn( triggerEnv, "process_id");
  params->shared_objectIndex = MetaioFindColumn( triggerEnv, "shared_object");
  params->lalwrapper_cvs_tagIndex = MetaioFindColumn( triggerEnv, "lalwrapper_cvs_tag");
  params->lal_cvs_tagIndex = MetaioFindColumn( triggerEnv, "lal_cvs_tag");
  params->commentIndex = MetaioFindColumn( triggerEnv, "comment");
  params->in_start_timeIndex = MetaioFindColumn( triggerEnv, "in_start_time");
  params->in_start_time_nsIndex = MetaioFindColumn( triggerEnv, "in_start_time_ns");
  params->in_end_timeIndex = MetaioFindColumn( triggerEnv, "in_end_time");
  params->in_end_time_nsIndex = MetaioFindColumn( triggerEnv, "in_end_time_ns");
  params->out_start_timeIndex = MetaioFindColumn( triggerEnv, "out_start_time");
  params->out_start_time_nsIndex = MetaioFindColumn( triggerEnv, "out_start_time_ns");
  params->out_end_timeIndex = MetaioFindColumn( triggerEnv, "out_end_time");
  params->out_end_time_nsIndex = MetaioFindColumn( triggerEnv, "out_end_time_ns");
  params->neventsIndex = MetaioFindColumn( triggerEnv, "nevents");
  params->nnodesIndex = MetaioFindColumn( triggerEnv, "nnodes");

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


void
getSearchSummaryTable(
        LALStatus                *status,
        const MetaioParseEnv      triggerEnv,
        SearchSummaryTable        *table,
        SearchSummaryIndex        *params
        )
{
    INT4 retVal=0;

    INITSTATUS (status, "getSearchSummaryTable", EVENT_UTILSC);
    ATTATCHSTATUSPTR (status);

    retVal = MetaioGetRow(triggerEnv);
    if ( retVal <= 0 ) {
        ABORT( status, EVENTUTILSH_EGETRO, EVENTUTILSH_MSGEGETRO);
    } 

    table->in_start_time.gpsSeconds = 
        triggerEnv->ligo_lw.table.elt[params->in_start_timeIndex].data.int_4s;  
    table->in_start_time.gpsNanoSeconds = 
        triggerEnv->ligo_lw.table.elt[params->in_start_time_nsIndex].data.int_4s;  
    table->in_end_time.gpsSeconds = 
        triggerEnv->ligo_lw.table.elt[params->in_end_timeIndex].data.int_4s;  
    table->in_end_time.gpsNanoSeconds = 
        triggerEnv->ligo_lw.table.elt[params->in_end_time_nsIndex].data.int_4s;  
    table->out_start_time.gpsSeconds = 
        triggerEnv->ligo_lw.table.elt[params->out_start_timeIndex].data.int_4s;  
    table->out_start_time.gpsNanoSeconds = 
        triggerEnv->ligo_lw.table.elt[params->out_start_time_nsIndex].data.int_4s;  
    table->out_end_time.gpsSeconds = 
        triggerEnv->ligo_lw.table.elt[params->out_end_timeIndex].data.int_4s;  
    table->out_end_time.gpsNanoSeconds = 
        triggerEnv->ligo_lw.table.elt[params->out_end_time_nsIndex].data.int_4s;  
    table->nevents = 
        triggerEnv->ligo_lw.table.elt[params->neventsIndex].data.int_4s;  
    table->nnodes = 
        triggerEnv->ligo_lw.table.elt[params->nnodesIndex].data.int_4s;  

    DETATCHSTATUSPTR (status);
    RETURN (status);
}




void
LALCompareSnglInspiral(
        LALStatus                *status,
        SnglInspiralTable        *aPtr,
        SnglInspiralTable        *bPtr,
        SnglInspiralAccuracy     *params
    )
{
    INT8 ta, tb;
    REAL4 dm1, dm2;
    REAL4 sigmaRatio;

    INITSTATUS (status, "LALCompareSnglInspiral", EVENT_UTILSC);
    ATTATCHSTATUSPTR (status);

    params->match=0;

    LALGPStoINT8( status->statusPtr, &ta, &(aPtr->end_time) );
    LALGPStoINT8( status->statusPtr, &tb, &(bPtr->end_time) );

    if( labs(ta-tb) < params->dtime ){

        dm1 = fabsf( aPtr->mass1 - bPtr->mass1 );
        dm2 = fabsf( aPtr->mass2 - bPtr->mass2 );

        if ( dm1 < params->dm && dm2 < params->dm ){

            sigmaRatio = sqrt(bPtr->sigmasq / aPtr->sigmasq);

            if ( sigmaRatio * aPtr->snr -  bPtr->snr < params->dRhoPlus &&
                    bPtr->snr - sigmaRatio * aPtr->snr < params->dRhoMinus ){
                params->match = 1;
            }
        }
    }

    DETATCHSTATUSPTR (status);
    RETURN (status);
}





/*******************************************************************
* Function builds a set of times that are vetoed according to the
* rules in vetoParams.
*******************************************************************/
int buildVetoTimes( vetoParams *thisentry)
{
    int status,iVetoS=-1,iVetoNS=-1,iVetoSNR=-1,first=1,timeBase=0,epoch=1;
    long stime,etime;
    timeWindow *thiswindow=NULL, *prevwindow=NULL;
    double tVtemp,SNR;
    struct MetaioParseEnvironment vetoParseEnv;
    const MetaioParseEnv vetoEnv = &vetoParseEnv;
    FILE *fpin;

    /* Allocate the head of the veto time list */
    if ((*thisentry).vwindows == NULL){
        thiswindow = (*thisentry).vwindows = (timeWindow *)malloc(sizeof(timeWindow));
        (*thiswindow).next_time_window = NULL;
        (*thiswindow).start_time=0.0;
        (*thiswindow).end_time=0.0;
        (*thiswindow).end_time=0.0;
        (*thiswindow).snr=0.0;
        (*thiswindow).ratio=(*thisentry).ratio;
    }

    epoch=strcmp((*thisentry).table_column,"epoch");

    if ( epoch != 0 ){
        /* Open the veto file */
        if ( (status = MetaioOpen( vetoEnv, (*thisentry).filename)) !=0 ){
            fprintf(stderr, "Error opening veto file %s\n", (*thisentry).filename );
            MetaioAbort( vetoEnv ); 
            return 2;
        }

        /* Locate the relevant columns */
        iVetoS   = MetaioFindColumn( vetoEnv, "start_time" );
        iVetoNS  = MetaioFindColumn( vetoEnv, "start_time_ns" );
        if ( iVetoS < 0 || iVetoNS < 0 ) {
            /* try to use end_time */
            iVetoS  = MetaioFindColumn( vetoEnv, "end_time" );
            iVetoNS = MetaioFindColumn( vetoEnv, "end_time_ns" );
            if ( iVetoS < 0 || iVetoNS < 0 ) {
                fprintf(stderr, "Veto file %s does not contain start_time or end_time\n",
                        (*thisentry).filename );
                MetaioAbort( vetoEnv ); 
                return 2;
            }
        }
        iVetoSNR = MetaioFindColumn( vetoEnv, (*thisentry).table_column);

        /* Read in the veto data and construct a list of veto times */
        while (1) {

            status = MetaioGetRow(vetoEnv);
            if ( status == -1 ) {
                printf( "Error while getting row from file %s\n", (*thisentry).filename );
                MetaioAbort( vetoEnv ); 
                return 6;
            } else if ( status == 0 ) {
                /*-- Reached end of file --*/
                break;
            }
            if ( strstr( (*thisentry).filename, "glitchmon" ) ){
                SNR = vetoEnv->ligo_lw.table.elt[iVetoSNR].data.real_8;
            }
            else{
                SNR = vetoEnv->ligo_lw.table.elt[iVetoSNR].data.real_4;
            }
            if ( SNR > (*thisentry).threshold ){

                tVtemp = (double) ( vetoEnv->ligo_lw.table.elt[iVetoS].data.int_4s - timeBase )
                    + 1.0e-9 * (double) vetoEnv->ligo_lw.table.elt[iVetoNS].data.int_4s;

                /* Must treat the first veto event different */
                if ( first ){
                    (*thiswindow).start_time = tVtemp - (*thisentry).minusdtime;
                    (*thiswindow).end_time = tVtemp + (*thisentry).plusdtime;
                    (*thiswindow).snr = SNR;
                    (*thiswindow).ratio=(*thisentry).ratio;
                    first = 0;
                }
                /* If this event is within last veto window,  update veto times */
                else if ( tVtemp <= ((*thiswindow).end_time + (*thisentry).minusdtime + 4.0) ){
                    (*thiswindow).end_time = tVtemp + (*thisentry).plusdtime;
                    if (SNR > (*thiswindow).snr){
                        (*thiswindow).snr=SNR;
                    }
                }
                /* Otherwise allocate next node and update veto window */
                else
                {
                    prevwindow=thiswindow;
                    (*thiswindow).next_time_window = (timeWindow *)malloc(sizeof(timeWindow));
                    thiswindow = (*thiswindow).next_time_window;
                    (*thiswindow).start_time = tVtemp - (*thisentry).minusdtime;
                    (*thiswindow).end_time = tVtemp + (*thisentry).plusdtime;
                    (*thiswindow).snr=SNR;
                    (*thiswindow).ratio=(*thisentry).ratio;
                    (*thiswindow).next_time_window = NULL;
                }

            }

        }

        MetaioAbort(vetoEnv);
    }
    else
    {

        /* Open the veto file */
        if ( ! (fpin = fopen((*thisentry).filename,"r") ) ){
            fprintf(stderr, "Error opening veto file %s\n", (*thisentry).filename );
            return 2;
        }

        while ( fscanf(fpin,"%i\t%i\n",&stime,&etime) != EOF ){

            /* Must treat the first veto event different */
            if ( first ){
                (*thiswindow).start_time = (double)stime;
                (*thiswindow).end_time = (double)etime;
                (*thiswindow).snr = 0.0;
                (*thiswindow).ratio = 1.0;
                first = 0;
            }
            /* Otherwise allocate next node and update veto window */
            else
            {
                prevwindow=thiswindow;
                (*thiswindow).next_time_window = (timeWindow *)malloc(sizeof(timeWindow));
                thiswindow = (*thiswindow).next_time_window;
                (*thiswindow).start_time = (double)stime;
                (*thiswindow).end_time = (double)etime;
                (*thiswindow).snr=0.0;
                (*thiswindow).ratio=1.0;
                (*thiswindow).next_time_window = NULL;
            }
        }

        /* close the veto file */
        fclose(fpin);

    }

    return 0;
}


/*******************************************************************
* Function takes a list of vetoes and resolves all veto times into 
* a single list
*******************************************************************/
int resolveVetoTimes( timeWindow **vwindows, vetoParams *thisentry)
{
    timeWindow *thiswindow=NULL,*prevwindow=NULL,*tempwindow=NULL;
    int updated_window = 0;

    /**************************************************************
     * Allocate memory for the first veto window
     *************************************************************/
    if (thisentry != NULL && (*thisentry).vwindows != NULL ){
        tempwindow = (*thisentry).vwindows;
        thiswindow = (*vwindows) = (timeWindow *)malloc( sizeof(timeWindow));
        (*thiswindow).next_time_window = NULL;
        (*thiswindow).start_time = (*tempwindow).start_time;
        (*thiswindow).end_time = (*tempwindow).end_time;
    }

    /**************************************************************
     * Loop over the vetoes melding into one list
     *************************************************************/
    while ( thisentry != NULL ){

        /* Loop over the time windows for each veto */
        tempwindow = (*thisentry).vwindows;
        while ( tempwindow != NULL ){
                
            double tempSTime=(*tempwindow).start_time;
            double tempETime=(*tempwindow).end_time;
            double tempSNR=(*tempwindow).snr;
            double tempRatio=(*tempwindow).ratio;

            /* Loop over concatenated list of vetoes */
            thiswindow = (*vwindows);
            updated_window = 0;
            while ( thiswindow != NULL){
                
                double thisSTime=(*thiswindow).start_time;
                double thisETime=(*thiswindow).end_time;
                double thisSNR=(*thiswindow).snr;
                double thisRatio=(*thiswindow).ratio;

                if ( tempSTime >= thisSTime && tempSTime <= thisETime){
                    (*thiswindow).end_time = tempETime > thisETime ? tempETime : thisETime;
                    (*thiswindow).snr = tempSNR > thisSNR ? tempSNR : thisSNR;
                    if (thisRatio != tempRatio){
                        fprintf(stderr,"Warning: different ratio factor\n");
                        fprintf(stderr,"Warning: Using one for latest window\n");
                        (*thiswindow).ratio = thisRatio;
                    }
                    updated_window=1;
                    break;
                }
                prevwindow = thiswindow;
                thiswindow = (*thiswindow).next_time_window;
                
            }
            if ( updated_window == 0 ){
                
                (*prevwindow).next_time_window = (timeWindow *)malloc( sizeof(timeWindow));
                thiswindow = (*prevwindow).next_time_window;
                (*thiswindow).start_time = tempSTime;
                (*thiswindow).end_time = tempETime;
                (*thiswindow).next_time_window = NULL;
            }

            tempwindow = (*tempwindow).next_time_window;
            
        }
        thisentry = (*thisentry).next_veto;

    }

    return 0;
}



/*******************************************************************
* Function builds an event list given veto information and 
* rules in candParams.
*******************************************************************/
int buildEventList( candEvent **eventhead, timeWindow *vwindows, candParams candidates,
        int injectflag, int maxflag, float calfudge)
{
    int status,iVetoS,iVetoNS,iVetoSNR,iCandCHISQ,iCandEDIST,iCandMCHIRP;
    int iCandMASS1,iCandMASS2;
    int pass,first=1,timeBase=0, crossingflag=0;
    timeWindow *thiswindow=NULL;
    double tVtemp,snrVtemp,chiVtemp,edistVtemp,lastVtemp,mchirpVtemp;
    double mass1Vtemp,mass2Vtemp;
    struct MetaioParseEnvironment candParseEnv;
    const MetaioParseEnv candEnv = &candParseEnv;
    vetoParams *thisentry=NULL;
    candEvent  *thisCEvent=NULL, *prevCEvent=NULL;
    float distfudge = 1.0, snrfudge = 1.0;

    if (injectflag == INJECTIONS){
        if ( (status = MetaioOpen( candEnv, candidates.injectfile)) !=0 ){
            fprintf(stderr, "Error opening injection file %s\n", candidates.injectfile );
            MetaioAbort( candEnv ); 
            exit(2);
        }
        snrfudge = calfudge;
        distfudge = 1.0;
    }else{
        if ( (status = MetaioOpen( candEnv, candidates.triggerfile)) !=0 ){
            fprintf(stderr, "Error opening trigger file %s\n", candidates.triggerfile );
            MetaioAbort( candEnv ); 
            exit(2);
        }
        snrfudge = 1.0;
        distfudge = calfudge;
    }


    /* Locate the relevant columns */
    iVetoS   = MetaioFindColumn( candEnv, "start_time" );
    iVetoNS  = MetaioFindColumn( candEnv, "start_time_ns" );
    if ( iVetoS < 0 || iVetoNS < 0 ) {
        /* try to use end_time */
        iVetoS  = MetaioFindColumn( candEnv, "end_time" );
        iVetoNS = MetaioFindColumn( candEnv, "end_time_ns" );
        if ( iVetoS < 0 || iVetoNS < 0 ) {
            fprintf(stderr, "File does not contain start_time or end_time\n");
            MetaioAbort( candEnv ); 
            return 2;
        }
    }
    iVetoSNR = MetaioFindColumn( candEnv, "SNR");
    iCandCHISQ = MetaioFindColumn( candEnv, "CHISQ");
    iCandEDIST = MetaioFindColumn( candEnv, "EFF_DISTANCE");
    iCandMCHIRP = MetaioFindColumn( candEnv, "MCHIRP");
    iCandMASS1 = MetaioFindColumn( candEnv, "MASS1");
    iCandMASS2 = MetaioFindColumn( candEnv, "MASS2");

    /* Read in the veto data and construct a list of veto times */
    first = 1;
    while (1) {

        /* assume candidate is an event unless proven otherwise */
        pass=1;
        status = MetaioGetRow(candEnv);
        if ( status == -1 ) {
            printf( "Error while getting row from injection or trigger file\n");
            MetaioAbort( candEnv ); 
            return 6;
        } else if ( status == 0 ) {
            /*-- Reached end of file --*/
            break;
        }

        snrVtemp = snrfudge * candEnv->ligo_lw.table.elt[iVetoSNR].data.real_4;  
        chiVtemp = candEnv->ligo_lw.table.elt[iCandCHISQ].data.real_4;  
        edistVtemp = distfudge * candEnv->ligo_lw.table.elt[iCandEDIST].data.real_4;  
        mchirpVtemp = candEnv->ligo_lw.table.elt[iCandMCHIRP].data.real_4;  
        mass1Vtemp = candEnv->ligo_lw.table.elt[iCandMASS1].data.real_4;  
        mass2Vtemp = candEnv->ligo_lw.table.elt[iCandMASS2].data.real_4;  
        /* Store the time of the candidate event */
        tVtemp = (double) ( candEnv->ligo_lw.table.elt[iVetoS].data.int_4s - timeBase )
            + 1.0e-9 * (double) candEnv->ligo_lw.table.elt[iVetoNS].data.int_4s;

        /* Require the SNR > threshold and CHISQ < threshold */
        if ( snrVtemp > candidates.snr_threshold
                && chiVtemp < candidates.chi_threshold ){

            /* Loop over the time windows for each veto */
            thiswindow = vwindows;
            while ( thiswindow != NULL ){
                /* Is the candidate in one of the vetoed list of times? */
                if ( tVtemp > (*thiswindow).start_time &&
                        tVtemp < (*thiswindow).end_time ){
                    /* this candidate is vetoed */
                        pass = 0;
                        break;
                }
                thiswindow = (*thiswindow).next_time_window;
            }

            if ( pass ){
                /* Must treat the first veto event difficult */
                if ( first ){
                    thisCEvent = (*eventhead) = (candEvent *)malloc(sizeof(candEvent)); 
                    (*thisCEvent).time = tVtemp;
                    (*thisCEvent).snr = snrVtemp;
                    (*thisCEvent).chisq = chiVtemp;
                    (*thisCEvent).eff_distance = edistVtemp;
                    (*thisCEvent).mchirp = mchirpVtemp;
                    (*thisCEvent).mass1 = mass1Vtemp;
                    (*thisCEvent).mass2 = mass2Vtemp;
                    (*thisCEvent).significance = 0;
                    (*thisCEvent).candidate = 0;
                    (*thisCEvent).coincident = 0;
                    first = 0;
                }
                /* If this event is within last veto window,  update veto times */
                else if ( tVtemp <= (lastVtemp + candidates.dtime) && maxflag){
                    if ( (*thisCEvent).snr < snrVtemp 
                            && (*thisCEvent).chisq > chiVtemp
                       ){
                        (*thisCEvent).time = tVtemp;
                        (*thisCEvent).snr = snrVtemp;
                        (*thisCEvent).chisq = chiVtemp;
                        (*thisCEvent).eff_distance = edistVtemp;
                        (*thisCEvent).mchirp = mchirpVtemp;
                        (*thisCEvent).mass1 = mass1Vtemp;
                        (*thisCEvent).mass2 = mass2Vtemp;
                        (*thisCEvent).significance = 0;
                        (*thisCEvent).candidate = 0;
                        (*thisCEvent).coincident = 0;
                    }
                }
                /* Otherwise allocate next node and update veto window */
                else
                {
                    (*thisCEvent).next_event = (candEvent *)malloc(sizeof(candEvent)); 
                    thisCEvent = (*thisCEvent).next_event;
                    (*thisCEvent).time = tVtemp;
                    (*thisCEvent).snr = snrVtemp;
                    (*thisCEvent).chisq = chiVtemp;
                    (*thisCEvent).eff_distance = edistVtemp;
                    (*thisCEvent).mchirp = mchirpVtemp;
                    (*thisCEvent).mass1 = mass1Vtemp;
                    (*thisCEvent).mass2 = mass2Vtemp;
                    (*thisCEvent).significance = 0;
                    (*thisCEvent).candidate = 0;
                    (*thisCEvent).coincident = 0;
                    (*thisCEvent).next_event = NULL;
                }
                lastVtemp = tVtemp; 
            }
        }
    }

    MetaioAbort(candEnv);

    return 0;
}


/*******************************************************************
 * Build an array giving data quality information for use in
 * coincidence studies
 *
 * TODO:  This should take a time series,  that way the sampling rate
 * would be carried around properly with it. 
 *******************************************************************/
int buildDataQaulity(int **coincident_times, snglIFO *ifo, int numIFO,
        double *dummyStart, double *dummyEnd)
{
    int ifoMask[1024],dummyMask=1,maskMax=0;
    int i,j,numPts;
    timeWindow *thiswindow=NULL;

    (*dummyEnd)   = ifo[0].end_time;
    (*dummyStart) = ifo[0].start_time;
    for ( i=0 ; i<numIFO ; i++ ){
        ifoMask[i] = dummyMask;
        dummyMask *= 2;
        maskMax += ifoMask[i];
        if ((*dummyStart) > ifo[i].start_time )
            (*dummyStart) = ifo[i].start_time;
        if ( (*dummyEnd) < ifo[i].end_time )
            (*dummyEnd) = ifo[i].end_time;
    }

    /* Sampling rate is 10 Hz */
    numPts = 10 * (int)((*dummyEnd)-(*dummyStart));
    (*coincident_times) = (int *)calloc( numPts, sizeof(int));

    for ( i=0 ; i<numIFO ; i++ ){
        dummyMask = ifoMask[i];
        thiswindow = ifo[i].awindows;
        while (thiswindow != NULL){
            int jmin,jmax,j;
            jmin = (int)( 10 * ( (*thiswindow).start_time-(*dummyStart) ));
            jmax = (int)( 10 * ( (*thiswindow).end_time-(*dummyStart) ));
            for(j=jmin; j<jmax ; j++) 
                (*coincident_times)[j]+=dummyMask;
            thiswindow =  (*thiswindow).next_time_window;
        }
        thiswindow = ifo[i].vwindows;
        while (thiswindow != NULL){
            int jmin,jmax,j;
            jmin = (int)( 10 * ( (*thiswindow).start_time-(*dummyStart) ));
            jmax = (int)( 10 * ( (*thiswindow).end_time-(*dummyStart) ));
            if (jmin < 0) {
                fprintf(stderr,"Warning: analysis and veto times out of range\n");
                jmin=0;
            }
            if (jmax < 0) {
                fprintf(stderr,"Warning: analysis and veto times out of range\n");
                jmax=0;
            }
            if (jmin > numPts) {
                fprintf(stderr,"Warning: analysis and veto times out of range\n");
                jmin=numPts;
            }
            if (jmax > numPts) {
                fprintf(stderr,"Warning: analysis and veto times out of range\n");
                jmax=numPts;
            }
            for(j=jmin; j<jmax ; j++){ 
                if ( (*coincident_times)[j] >= dummyMask ){
                    (*coincident_times)[j] -= dummyMask;
                }
            }
            thiswindow =  (*thiswindow).next_time_window;
        }
    }
    return 0;
}

/*******************************************************************
* Function builds multiple inspiral event lists using coincidence.  It
* is not general purpose,  but tuned to S1 analysis.  This needs to be
* fixed.  TODO
*******************************************************************/
int cpySnglToMultiInspiral(multiInspiral *thisMEvent, candEvent *myevent, int ifo)
{

    thisMEvent->snr[ifo]=myevent->snr;
    thisMEvent->time[ifo]=myevent->time;
    thisMEvent->effDistance[ifo]=myevent->eff_distance;
    thisMEvent->mchirp[ifo]=myevent->mchirp;

    return 0;
}

int buildMultiInspiralEvents(multiInspiral **multInspEv, int *coincident_times,
        snglIFO *ifo, int numIFO, int injectflag, double dummyStart, float delm,
        float distance, float coincidence_window)
{
    int i, numEvents=0, first=1,dummyMask;
    candEvent *thisCEvent=NULL, *myevent=NULL;
    multiInspiral *thisMEvent;
    FILE *fpout;

    fprintf(stderr,"distance = %f, delm = %f, coinc = %f\n", distance, delm, coincidence_window);

    /* deteremine coincidences according to our rule */
    for ( i=1 ; i<numIFO ; i++ ){
        myevent = (injectflag) ? ifo[0].Ieventhead : ifo[0].eventhead;
        while ( myevent != NULL ){
            if ( (*myevent).eff_distance < distance &&
                    coincident_times[(int)( 10 * ((*myevent).time - dummyStart))] == 3 ){
                thisCEvent = (injectflag) ? ifo[i].Ieventhead : ifo[i].eventhead;
                while ( thisCEvent != NULL && (*myevent).significance == 0){
                    float chirpfrac=( (*myevent).mchirp - (*thisCEvent).mchirp )/
                        (*thisCEvent).mchirp;
                    float m1frac=( (*myevent).mass1 - (*thisCEvent).mass1 )/
                        (*thisCEvent).mass1;
                    float m2frac=( (*myevent).mass2 - (*thisCEvent).mass2 )/
                        (*thisCEvent).mass2;
                    if ( chirpfrac*chirpfrac < delm ){ 
                    /* if ( m1frac*m1frac < delm && m2frac*m2frac < delm ){ */
                        if ( (*myevent).time > (*thisCEvent).time - coincidence_window &&
                                (*myevent).time < (*thisCEvent).time + coincidence_window){
                            numEvents++;
                            thisMEvent = (multiInspiral *)malloc(sizeof(multiInspiral));
                            thisMEvent->next_event = NULL;
                            if (first){
                                (*multInspEv) = thisMEvent; first=0;
                            }
                            cpySnglToMultiInspiral(thisMEvent,myevent,0);
                            cpySnglToMultiInspiral(thisMEvent,thisCEvent,i);
                            thisMEvent=thisMEvent->next_event;
                            (*thisCEvent).significance = (*myevent).significance = 3;
                            (*thisCEvent).coincident = (*myevent).coincident = 1;
                        }
                    }    
                    thisCEvent = (*thisCEvent).next_event;
                }
            } 
            else if ( coincident_times[(int)(10 * ((*myevent).time - dummyStart))] == 3 ){
                (*myevent).significance = 3;
            }
            myevent = (*myevent).next_event;
        }
    }

    /* Identify singles in each of the insterferometers */
    dummyMask=1;
    for ( i=0 ; i<numIFO ; i++ ){
        thisCEvent = (injectflag) ? ifo[i].Ieventhead : ifo[i].eventhead;
        while ( thisCEvent != NULL ) {
            if ( coincident_times[(int)(10 * ((*thisCEvent).time - dummyStart))] == dummyMask ){
                (*thisCEvent).significance = dummyMask;
            }
            thisCEvent = (*thisCEvent).next_event;
        }
        dummyMask *= 2;
    }

    return numEvents;
}

/*******************************************************************
 * Print out the events which have been flagged as significant
 ******************************************************************/
int printInspiralEvents(FILE *fpout, snglIFO *ifo, int significance, int injectflag) 
{
    double tmpTime;
    float  time_interval;
    int    i;
    candEvent *dumEvent=NULL,*myevent=NULL;

    myevent = (injectflag) ? ifo->Ieventhead : ifo->eventhead;
    while ( myevent != NULL && (*myevent).significance != significance ) {
        myevent = (*myevent).next_event;
    }

    if ( myevent == NULL )  return 0;

    time_interval=ifo->candidates.dtime;
    dumEvent = myevent;
    tmpTime = (*dumEvent).time;
    myevent = (*myevent).next_event;
    i=1;
    fprintf(fpout,"# Time  SNR  CHISQ  EFF_DIST  MCHIRP  COINC\n");
    while ( myevent != NULL ){
        if ( (*myevent).significance == significance ){
            if ( (*myevent).time <= ((*dumEvent).time + time_interval) ){
                if ( (*myevent).snr > (*dumEvent).snr ){
                    dumEvent = myevent;
                }
            } else {
                fprintf(fpout,"%i %f %f %f %f %f %i %i %f %f\n",i, (*dumEvent).time,
                        (*dumEvent).snr,(*dumEvent).chisq,
                        (*dumEvent).eff_distance,(*dumEvent).mchirp,
                        (*dumEvent).coincident, (*dumEvent).significance,
                        (*dumEvent).mass1, (*dumEvent).mass2);
                (*dumEvent).candidate = 1;
                dumEvent = myevent;
                i++;
            }
            tmpTime = (*myevent).time;
        }
        myevent = (*myevent).next_event;
    }
    fprintf(fpout,"%i %f %f %f %f %f %i %i %f %f\n",i, (*dumEvent).time,
            (*dumEvent).snr,(*dumEvent).chisq,
            (*dumEvent).eff_distance,(*dumEvent).mchirp,
            (*dumEvent).coincident, (*dumEvent).significance,
                        (*dumEvent).mass1, (*dumEvent).mass2);
    (*dumEvent).candidate = 1;

    return 0;
}




/*******************************************************************
* Function builds a 2 dimensional histogram based on the event list
* passed into it.
*******************************************************************/
int build2DHistogram(candEvent *eventhead, const char *outputfile,
        int **histogram, int numbins, float minsnr, float maxchisq)
{
    FILE *fpout;
    float *snrbin, *chisqbin, snrdiv, chisqdiv;
    float maxsnr=minsnr+10.0, minchisq=0.0;
    int i,j,k;
    candEvent *thisCEvent=NULL;

    fpout = fopen(outputfile,"w");

    /* set up the bin boundaries for the histogram */
    snrbin     = calloc( numbins + 1, sizeof(float) );
    chisqbin   = calloc( numbins + 1, sizeof(float) );
    snrdiv     = ( maxsnr - minsnr ) / (float) numbins;
    chisqdiv   = ( maxchisq - minchisq ) / (float) numbins;
    for ( i = 0; i <= numbins; ++i )
    {
        snrbin[i] = minsnr + (float) i * snrdiv;
        chisqbin[i] = minchisq + (float) i * chisqdiv;
        /*fprintf( stderr, "snrbin[%d] = %f\tchisqbin[%d] = %f\n", 
                i, snrbin[i], i, chisqbin[i] ); */
    }

    for ( j = 0; j < numbins; ++j )
    {
        for ( thisCEvent = eventhead; thisCEvent; thisCEvent = thisCEvent->next_event )
        {
            if ( thisCEvent->chisq > chisqbin[j + 1] ) continue;
            for ( k = 0; k < numbins; ++k )
            {
                if ( thisCEvent->snr >= snrbin[k] ) ++histogram[k][j];
            }
        }
    }

    fprintf(fpout,"# minsnr = %f, snrdiv = %f\n",minsnr,snrdiv);
    fprintf(fpout,"# minchisq = %f, chisqdiv = %f\n",minchisq,chisqdiv);
    for ( j = 0; j < numbins; ++j )
    {
        for ( k = 0; k < numbins; ++k )
        {
            fprintf( fpout, "%d\t", histogram[k][j] );
        }
        fprintf( fpout, "\n" );
    }

    fclose(fpout);
    free(snrbin);
    free(chisqbin);
    return 0;
}

/************************************************************************
 * Compute the event rate limit using the loudest event method
 *
 * NOTE:  the histograms are indexed by histogram[snr][chisq]
 *
 ***********************************************************************/
int computeUL(const char *outputfile, int **triggerHistogram,
        int **injectHistogram, int numbins, float minsnr, float maxchisq,
        int ninject, float time_analyzed)
{
    FILE *fpout;
    float *snrbin, *chisqbin, snrdiv, chisqdiv, snrloud;
    float maxsnr=minsnr+10.0, minchisq=0.0, efficiency;
    double mu;
    int i,j,k,snrloudk=numbins;
    Chi2ThresholdIn  thresholdIn;
    static LALStatus        status;

    /* set up the bin boundaries for the histogram */
    snrbin     = calloc( numbins + 1, sizeof(float) );
    chisqbin   = calloc( numbins + 1, sizeof(float) );
    snrdiv     = ( maxsnr - minsnr ) / (float) numbins;
    chisqdiv   = ( maxchisq - minchisq ) / (float) numbins;
    
    for ( i = 0; i <= numbins; ++i )
    {
        snrbin[i] = minsnr + (float) i * snrdiv;
        chisqbin[i] = minchisq + (float) i * chisqdiv;
        /*fprintf( stderr, "snrbin[%d] = %f\tchisqbin[%d] = %f\n", 
                i, snrbin[i], i, chisqbin[i] ); */
    }


    /*****************************************************************
     * compute the upper limit using the loudest event
     ****************************************************************/
    fpout = fopen(outputfile,"w");
    for ( j = 0; j < numbins; ++j )
    {
        for ( k = 0; k < numbins; ++k )
        {
            if ( triggerHistogram[k][j] == 0 ){
                snrloud = snrbin[k];
                snrloudk = k;
                break;
            }
        }
        efficiency = ((float)injectHistogram[snrloudk][j]+1e-10)/((float)ninject);
        fprintf(fpout,"%f %f %f %f\n", chisqbin[j],
                3.89/(efficiency*time_analyzed), efficiency, snrloud);
    }
    fclose(fpout);

    /*****************************************************************
     * compute the upper limit using the counting method
     ****************************************************************/
    fpout = fopen("counting.rates","w");
    j = numbins-1;
    
    for ( k = 0; k < numbins; ++k )
    {
        thresholdIn.falseAlarm=0.1;
        thresholdIn.dof=2.0*(triggerHistogram[k][j]+1);
        LALChi2Threshold(&status, &mu, &thresholdIn);
        mu /= 2.0;
        efficiency = ((float)injectHistogram[k][j]+1e-10)/((float)ninject);
        fprintf(fpout,"%f %f %f\n", snrbin[k],
                ((float)mu)/(efficiency*time_analyzed), chisqbin[j]);
    }
    fclose(fpout);
    free(snrbin);
    free(chisqbin);
    return 0;
}


/***********************************************************************
 * 
 * Function returns the appropriate data quality bit for interferometer
 * 
 ***********************************************************************/
void LALDQBit(char *ifo, int *dqbit){

    if ( ! strcmp(ifo,"L1") ){
        *dqbit = 10;
        return;
    }
    else if ( ! strcmp(ifo,"H1") ){
        *dqbit = 8;
        return;
    }
    else {
        *dqbit = -1;
        return;
    }
}


static int EventCompare( const void *t1, const void *t2 )
{
  candEvent * const *tiles1 = t1;
  candEvent * const *tiles2 = t2;
  if ( (*tiles1)->time > (*tiles2)->time )
    return 1;
  if ( (*tiles1)->time < (*tiles2)->time )
    return -1;
  return 0;
}

/******** <lalVerbatim file="SortTFTilingCP"> ********/
void
LALSortTriggers (
	      LALStatus         *status,
              snglIFO           *ifo, 
              int               numIFO
	      )
/******** </lalVerbatim> ********************************/
{
  INT4               eventCount;
  INT4               numEvents;
  INT4               i,inject;
  candEvent             *thisEvent;
  candEvent             **events;

  INITSTATUS (status, "LALSortTFTiling", EVENT_UTILSC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not NULL */

  /* make sure excess power has already been computed */
  for (i=0 ; i<numIFO ; i++){
      for (inject=0 ; inject<2; inject++){

          /* compute number of tiles */
          if (inject){
              thisEvent = ifo[i].Ieventhead;
          } else {
              thisEvent = ifo[i].eventhead;
          }
          eventCount=0;
          while (thisEvent != NULL)
          {
              eventCount++;
              thisEvent = thisEvent->next_event;
          }
          numEvents = eventCount;

          /* 
           *
           *  Make an array of pointers to be used to sort the tiles.
           *
           */

          /* allocate memory for array of pointers to tiles */
          events = NULL;
          fprintf(stderr,"numEvents=%i\n",numEvents);
          events = (candEvent **) LALMalloc (numEvents * sizeof(candEvent *));

          /*  Make sure that the allocation was succesful */
          if ( !(events) ){
              ABORT (status, EVENTUTILSH_ENULLP, EVENTUTILSH_MSGENULLP);
          }

          /* copy out pointers into array */
          eventCount=0;
          if (inject){
              thisEvent = ifo[i].Ieventhead;
          } else {
              thisEvent = ifo[i].eventhead;
          }
          while (thisEvent != NULL)
          {
              eventCount++;
              *(events + eventCount-1) = thisEvent;
              thisEvent = thisEvent->next_event;
          }

          qsort( events, numEvents, sizeof( candEvent * ), EventCompare );

          /* copy sorted array back into linked list */
          { 
              candEvent **currentEvent = NULL;
              if (inject){
                  thisEvent = ifo[i].Ieventhead;
              } else {
                  thisEvent = ifo[i].eventhead;
              }
              currentEvent = &(thisEvent);
              

              eventCount=0;
              while (eventCount < numEvents)
              {
                  *currentEvent = *(events + eventCount);
                  eventCount++;
                  currentEvent = &((*currentEvent)->next_event);
              }

              /* correctly terminate the linked list */
              *currentEvent = NULL;
          }

          LALFree (events);
      }
  }

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}



static int EventCompareTime( const void *t1, const void *t2 )
{
  static LALStatus   stat;
  SnglInspiralTable * const *event1 = t1;
  SnglInspiralTable * const *event2 = t2;
  INT8 time1, time2;

  
  LALGPStoINT8(&stat, &time1, &((*event1)->end_time));
  LALGPStoINT8(&stat, &time2, &((*event2)->end_time));
  if ( time1 > time2 )
    return 1;
  if ( time1 < time2 )
    return -1;
  return 0;
}

/******** <lalVerbatim file="SortTFTilingCP"> ********/
void
LALSortSnglInspiralTable (
	      LALStatus         *status,
              SnglInspiralTable *inspiralEvent,
              INT4               numEvents
	      )
/******** </lalVerbatim> ********************************/
{
  INT4                   eventCount;
  INT4                   i,inject;
  SnglInspiralTable     *thisEvent=NULL;
  SnglInspiralTable    **events=NULL;

  INITSTATUS (status, "LALSortSnglInspiralTable", EVENT_UTILSC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not NULL */


  /* 
   *
   *  Make an array of pointers to be used to sort the tiles.
   *
   */

  /* allocate memory for array of pointers to tiles */
  events = NULL;
  events = (SnglInspiralTable **) 
      LALMalloc (numEvents * sizeof(SnglInspiralTable *));

  /*  Make sure that the allocation was succesful */
  if ( !(events) ){
      ABORT (status, EVENTUTILSH_ENULLP, EVENTUTILSH_MSGENULLP);
  }

  /* copy out pointers into array */
  eventCount=0;
  thisEvent = inspiralEvent;
  while (thisEvent != NULL)
  {
      eventCount++;
      *(events + eventCount-1) = thisEvent;
      thisEvent = thisEvent->next;
  }

  qsort( events, numEvents, sizeof( SnglInspiralTable * ), EventCompareTime );

  /* copy sorted array back into linked list */
  { 
      SnglInspiralTable **currentEvent = NULL;

      thisEvent = inspiralEvent;
      currentEvent = &(thisEvent);

      eventCount=0;
      while (eventCount < numEvents)
      {
          *currentEvent = *(events + eventCount);
          eventCount++;
          currentEvent = &((*currentEvent)->next);
      }

      /* correctly terminate the linked list */
      *currentEvent = NULL;
  }

  LALFree (events);

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}



void
LALClusterSnglInspiralTable (
	      LALStatus         *status,
              SnglInspiralTable *inspiralEvent,
              INT4              dtime
	      )
{
  SnglInspiralTable     *thisEvent=NULL,*prevEvent=NULL;

  INITSTATUS (status, "LALClusterSnglInspiralTable", EVENT_UTILSC);
  ATTATCHSTATUSPTR (status);

  thisEvent = inspiralEvent->next;
  prevEvent = inspiralEvent;
  while (thisEvent != NULL)
  {
      INT8 currTime, prevTime;

      /* compute the time in nanosec for each event trigger */
      LALGPStoINT8(status->statusPtr, &currTime, &(thisEvent->end_time));
      CHECKSTATUSPTR(status);
      LALGPStoINT8(status->statusPtr, &prevTime, &(prevEvent->end_time));
      CHECKSTATUSPTR(status);

      /* find events within the cluster window */
      if ( (currTime - prevTime) < 1000000LL * dtime){
          /* displace previous event in cluster ...... */
          if(thisEvent->snr > prevEvent->snr &&
                  thisEvent->chisq < prevEvent->chisq){
              memcpy( prevEvent, thisEvent, sizeof(SnglInspiralTable));
          }
          /* otherwise just dump this event from cluster */
              prevEvent->next = thisEvent->next;
              LALFree(thisEvent);
              thisEvent = prevEvent->next;
      }
      else {
          /* otherwise we keep this unique event trigger */
          prevEvent = thisEvent;
          thisEvent = thisEvent->next;
      }
  }

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


