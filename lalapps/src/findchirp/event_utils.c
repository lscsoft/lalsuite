#include "trump.h"

/*******************************************************************
* Function builds a set of times that are vetoed according to the
* rules in vetoParams.
*******************************************************************/
int buildVetoTimes( vetoParams *thisentry)
{
    int status,iVetoS,iVetoNS,iVetoSNR,first=1,timeBase=0;
    timeWindow *thiswindow=NULL, *prevwindow=NULL;
    double tVtemp;
    struct MetaioParseEnvironment vetoParseEnv;
    const MetaioParseEnv vetoEnv = &vetoParseEnv;

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

    /* Allocate the head of the veto time list */
    if ((*thisentry).vwindows == NULL){
        thiswindow = (*thisentry).vwindows = (timeWindow *)malloc(sizeof(timeWindow));
        (*thiswindow).next_time_window = NULL;
        (*thiswindow).start_time=0.0;
        (*thiswindow).end_time=0.0;
    }

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
        if (vetoEnv->ligo_lw.table.elt[iVetoSNR].data.real_4 > (*thisentry).threshold ){
           
            tVtemp = (double) ( vetoEnv->ligo_lw.table.elt[iVetoS].data.int_4s - timeBase )
                + 1.0e-9 * (double) vetoEnv->ligo_lw.table.elt[iVetoNS].data.int_4s;

            /* Must treat the first veto event difficult */
            if ( first ){
                (*thiswindow).start_time = tVtemp - (*thisentry).dtime;
                (*thiswindow).end_time = tVtemp + (*thisentry).dtime;
                first = 0;
            }
            /* If this event is within last veto window,  update veto times */
            else if ( tVtemp <= ((*thiswindow).end_time + (*thisentry).dtime) ){
                (*thiswindow).end_time = tVtemp + (*thisentry).dtime;
            }
            /* Otherwise allocate next node and update veto window */
            else
            {
                prevwindow=thiswindow;
                (*thiswindow).next_time_window = (timeWindow *)malloc(sizeof(timeWindow));
                thiswindow = (*thiswindow).next_time_window;
                (*thiswindow).start_time = tVtemp - (*thisentry).dtime;
                (*thiswindow).end_time = tVtemp + (*thisentry).dtime;
                (*thiswindow).next_time_window = NULL;
            }
            
        }
        
    }

    MetaioAbort(vetoEnv);
    
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

            /* Loop over concatenated list of vetoes */
            thiswindow = (*vwindows);
            updated_window = 0;
            while ( thiswindow != NULL){
                
                double thisSTime=(*thiswindow).start_time;
                double thisETime=(*thiswindow).end_time;
                if ( tempSTime >= thisSTime && tempSTime <= thisETime){
                    (*thiswindow).end_time = tempETime > thisETime ? tempETime : thisETime;
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
        int injectflag)
{
    int status,iVetoS,iVetoNS,iVetoSNR,iCandCHISQ,pass,first=1,timeBase=0;
    timeWindow *thiswindow=NULL;
    double tVtemp,snrVtemp,chiVtemp;
    struct MetaioParseEnvironment candParseEnv;
    const MetaioParseEnv candEnv = &candParseEnv;
    vetoParams *thisentry=NULL;
    candEvent  *thisCEvent=NULL, *prevCEvent=NULL;

    if (injectflag == INJECTIONS){
        if ( (status = MetaioOpen( candEnv, candidates.injectfile)) !=0 ){
            fprintf(stderr, "Error opening veto file %s\n", candidates.injectfile );
            MetaioAbort( candEnv ); 
            return 2;
        }
    }else{
        if ( (status = MetaioOpen( candEnv, candidates.triggerfile)) !=0 ){
            fprintf(stderr, "Error opening veto file %s\n", candidates.triggerfile );
            MetaioAbort( candEnv ); 
            return 2;
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
        /* Require the SNR > threshold and CHISQ < threshold */
        if ( snrVtemp > candidates.snr_threshold
                && chiVtemp < candidates.chi_threshold ){

            /* Store the time of the candidate event */
            tVtemp = (double) ( candEnv->ligo_lw.table.elt[iVetoS].data.int_4s - timeBase )
                + 1.0e-9 * (double) candEnv->ligo_lw.table.elt[iVetoNS].data.int_4s;

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
                    first = 0;
                }
                /* If this event is within last veto window,  update veto times */
                else if ( tVtemp <= ((*thisCEvent).time + candidates.dtime) ){
                    if ( (*thisCEvent).snr < snrVtemp && (*thisCEvent).chisq > chiVtemp){
                        (*thisCEvent).time = tVtemp;
                        (*thisCEvent).snr = snrVtemp;
                        (*thisCEvent).chisq = chiVtemp;
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
                    (*thisCEvent).next_event = NULL;
                }

            }

        }

    }

    MetaioAbort(candEnv);

    return 0;
}

/*******************************************************************
* Function builds a 2 dimensional histogram based on the event list
* passed into it.
*******************************************************************/
int build2DHistogram(candEvent *eventhead, const char *outputfile,
        int **histogram, int numbins){
    
    FILE *fpout;
    float *snrbin, *chisqbin, snrdiv, chisqdiv;
    float maxsnr=17.0, minsnr=7.0, maxchisq=10.0, minchisq=0.0;
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

    for ( j = 0; j < numbins; ++j )
    {
        for ( k = 0; k < numbins; ++k )
        {
            fprintf( fpout, "%d\t", histogram[k][j] );
        }
        fprintf( fpout, "\n" );
    }

    fclose(fpout);

    return 0;
}

/************************************************************************
 * Compute the event rate limit using the loudest event method
 ***********************************************************************/


