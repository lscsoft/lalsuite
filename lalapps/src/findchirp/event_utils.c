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


NRCSID (SORTTFTILINGC, "$Id$");


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/Thresholds.h>
#include "metaio.h"
#include "donald.h"
#include "event_utils.h"


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
        int injectflag, int maxflag)
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

    if (injectflag == INJECTIONS){
        if ( (status = MetaioOpen( candEnv, candidates.injectfile)) !=0 ){
            fprintf(stderr, "Error opening injection file %s\n", candidates.injectfile );
            MetaioAbort( candEnv ); 
            exit(2);
        }
    }else{
        if ( (status = MetaioOpen( candEnv, candidates.triggerfile)) !=0 ){
            fprintf(stderr, "Error opening trigger file %s\n", candidates.triggerfile );
            MetaioAbort( candEnv ); 
            exit(2);
        }
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

        snrVtemp = candEnv->ligo_lw.table.elt[iVetoSNR].data.real_4;  
        chiVtemp = candEnv->ligo_lw.table.elt[iCandCHISQ].data.real_4;  
        edistVtemp = candEnv->ligo_lw.table.elt[iCandEDIST].data.real_4;  
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
        snglIFO *ifo, int numIFO, int injectflag, double dummyStart)
{
    float coincidence_window = COINWINDOW;
    float delm = DELM;
    float distance = H1L1DISTANCE;
    int i, numEvents=0, first=1,dummyMask;
    candEvent *thisCEvent=NULL, *myevent=NULL;
    multiInspiral *thisMEvent;
    FILE *fpout;

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
                    if ( chirpfrac*chirpfrac < delm ){
                    /* if ((*myevent).mchirp > (*thisCEvent).mchirp - delm &&
                            (*myevent).mchirp < (*thisCEvent).mchirp + delm){ */
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

  INITSTATUS (status, "LALSortTFTiling", SORTTFTILINGC);
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

