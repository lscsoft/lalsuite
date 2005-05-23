#include <lal/LALStdlib.h>
#include <lal/LALNoiseModels.h>
#include <lal/LALConstants.h>

NRCSID (LALDISCOVERINSPIRALEVENTSC, "$Id$");

static void findValidInspiralEvents (REAL4Vector, REAL4Vector, DiscoverInspiralEventsIn, INT4 *, INT4 *); 
static void findMax ( REAL4 *vec, INT4 *ok, INT4 n, REAL4 *mVal, INT4 *mPos);
static INT4 computeSum (INT4 *ivec, INT4 n);

void  LALDiscoverInspiralEvents ( 
        LALStatus                    *status,
        INT4                         *nEvents,
        DiscoverInspiralEventsList   **eventlist,
        DiscoverInspiralEventsIn     *findeventsin )
{ 
    INT4                               i, j, nBegin, nEnd;
    REAL8                              dt, df, distanceNorm, eSec;
    REAL4Vector                        filter1, filter2, output1, output2, correlation; 

    InspiralChisqInput                 chisqin;
    InspiralChisqOutput                *chisqout=NULL;

    INT4                               nTempEvents;
    INT4                               *tempEventsIdx=NULL;

    InspiralWaveOverlapIn              overlapin;
    InspiralWaveOverlapOut             overlapout;

    INITSTATUS (status, "LALDiscoverInspiralEvents", LALDISCOVERINSPIRALEVENTSC);
    ATTATCHSTATUSPTR(status);

    ASSERT (findeventsin->psd.data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
    ASSERT (findeventsin->signal.data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);

    output1.length = output2.length = findeventsin->signal.length;
    filter1.length = filter2.length = findeventsin->signal.length;

    ASSERT (findeventsin->nEnd >= 0,  status, LALNOISEMODELSH_ECHOICE, LALNOISEMODELSH_MSGECHOICE);
    ASSERT (findeventsin->nEnd <= (INT4)output1.length,  status, LALNOISEMODELSH_ECHOICE, LALNOISEMODELSH_MSGECHOICE);
    ASSERT (findeventsin->nBegin >= 0,  status, LALNOISEMODELSH_ECHOICE, LALNOISEMODELSH_MSGECHOICE);
    ASSERT (findeventsin->nBegin <= (INT4)output1.length,  status, LALNOISEMODELSH_ECHOICE, LALNOISEMODELSH_MSGECHOICE);
    
    /* Check for reasonable values of chisqBins - 0 <= chisqBins <= 20 */
    ASSERT (findeventsin->chisqBins >= 0,  status, LALNOISEMODELSH_ECHOICE, LALNOISEMODELSH_MSGECHOICE);
    ASSERT (findeventsin->chisqBins <= 20,  status, LALNOISEMODELSH_ECHOICE, LALNOISEMODELSH_MSGECHOICE);

    overlapin.nBegin = findeventsin->nBegin;
    overlapin.nEnd   = findeventsin->nEnd;
    overlapin.signal = findeventsin->signal;
    overlapin.psd    = findeventsin->psd;
    overlapin.param  = findeventsin->param;
    overlapin.fwdp   = findeventsin->fwdp;
    overlapin.revp   = findeventsin->revp;
    overlapin.ifExtOutput = 1;

    correlation.length = overlapin.signal.length;;
    correlation.data   = (REAL4*) LALMalloc(sizeof(REAL4)*correlation.length);
    output1.length     = overlapin.signal.length;;
    output1.data       = (REAL4*) LALMalloc(sizeof(REAL4)*output1.length);
    output2.length     = overlapin.signal.length;;
    output2.data       = (REAL4*) LALMalloc(sizeof(REAL4)*output2.length);
    filter1.length     = overlapin.signal.length;;
    filter1.data       = (REAL4*) LALMalloc(sizeof(REAL4)*filter1.length);
    output2.length     = overlapin.signal.length;;
    filter2.data       = (REAL4*) LALMalloc(sizeof(REAL4)*filter2.length);

    overlapout.xcorr1  = &output1;
    overlapout.xcorr2  = &output2;
    overlapout.filter1 = &filter1;
    overlapout.filter2 = &filter2;

    /* Actually carry out the correlation here ... */
    LALInspiralWaveOverlap(status->statusPtr, &correlation, &overlapout, &overlapin);
    CHECKSTATUSPTR (status);

    /************************************************************************************/
    /*--------------- Are we displaying Correlation statistics ? -----------------------*/
    /************************************************************************************/
    if (findeventsin->displayCorrelationStats) {

        /* These are required in case it is required to calculate the statistics */
        REAL4Vector              statsBuffer;
        StatsREAL4VectorOut      statsout1;
        StatsREAL4VectorOut      statsout2;
        StatsREAL4VectorOut      statsout3;

        statsBuffer.length = findeventsin->signal.length - (findeventsin->nEnd + findeventsin->nBegin);
        statsBuffer.data   = NULL;
        statsBuffer.data   = (REAL4*) LALMalloc(sizeof(REAL4)*statsBuffer.length);

        for (i=findeventsin->nBegin; i<(findeventsin->signal.length - findeventsin->nEnd); i++) 
              statsBuffer.data[i-findeventsin->nBegin] = output1.data[i];
        LALStatsREAL4Vector(status->statusPtr, &statsout1, &statsBuffer);
        CHECKSTATUSPTR(status);

        for (i=findeventsin->nBegin; i<(findeventsin->signal.length - findeventsin->nEnd); i++) 
              statsBuffer.data[i-findeventsin->nBegin] = output2.data[i];
        LALStatsREAL4Vector(status->statusPtr, &statsout2, &statsBuffer);
        CHECKSTATUSPTR(status);

        for (i=findeventsin->nBegin; i<(findeventsin->signal.length - findeventsin->nEnd); i++) 
              statsBuffer.data[i-findeventsin->nBegin] = ( (output1.data[i]*output1.data[i]) +
                      (output2.data[i]*output2.data[i]) );
        LALStatsREAL4Vector(status->statusPtr, &statsout3, &statsBuffer);
        CHECKSTATUSPTR(status);

        fprintf(stderr, "Quad 1     mean=%.5e std=%.5e min=%.5e max=%.5e\n",
                statsout1.mean, statsout1.stddev, statsout1.min, statsout1.max);   
        fprintf(stderr, "Quad 2     mean=%.5e std=%.5e min=%.5e max=%.5e\n",
                statsout2.mean, statsout2.stddev, statsout2.min, statsout2.max);   
        fprintf(stderr, "Rho square mean=%.5e std=%.5e min=%.5e max=.5%e\n",
                statsout3.mean, statsout3.stddev, statsout3.min, statsout3.max);

        if (statsBuffer.data != NULL) LALFree (statsBuffer.data);
    } 
    /*-------------------------------------------------------------------------------------*/

    /*************************************************/
    /*--------- Are we displaying Correlation ? -----*/
    /*************************************************/
    if (findeventsin->displayCorrelation) {
        FILE *fptr=NULL;
        fptr = fopen("correlation.out","w");

        dt     = 1.L/findeventsin->param.tSampling;
        for (i=0; i<findeventsin->signal.length; i++) {
            fprintf (fptr, "%.8e \t %.8e \t  %.8e \t  %.8e\n",
                    dt*(REAL4)(i), output1.data[i], output2.data[2],
                    ((output1.data[i]*output1.data[i]) + (output2.data[i]*output2.data[i]))
                    );
        }
        fclose(fptr);
        fprintf (stderr, "Written correlation to [1] correlation.out\n");
    }
    /*----------------------------------------------*/

    /* We will now calculate the distanceNorm. It is the effDistance for snr of 1.0 */
    dt     = 1.L/findeventsin->param.tSampling;
    df     = 1.L/(output1.length * dt);
    nBegin = findeventsin->nBegin;
    nEnd   = findeventsin->signal.length - findeventsin->nEnd;

    LALEstimateEffectiveDistance (status->statusPtr, findeventsin->param, df, &findeventsin->psd, 1.0, &distanceNorm);
    CHECKSTATUSPTR (status);

    /* Remember to scale things by dynRangeScalingFactor */
    distanceNorm *= sqrt(findeventsin->dynRangeScalingFac);

    /* We now pick loudest events separated by +/- one coalescence time */
    nTempEvents = 0;
    tempEventsIdx = (INT4 *) LALCalloc (1, sizeof(INT4)*(nEnd-nBegin+1));
    findValidInspiralEvents (output1, output2, *findeventsin, &nTempEvents, tempEventsIdx) ;

    /* If some events are returned, calculate chisq and return results in eventlist */
    if (nTempEvents) {

        /* Calculate the chisq in one shot for ALL time lags if chisqBins is
         * greaer than zero. If it is equal to zero, it means chisq need not
         * be calculated. 
         * */
        if (findeventsin->chisqBins) {
            /* Calculate chisquare */
            chisqin.chisqBins    = findeventsin->chisqBins;
            chisqin.flso         = 1.L/(pow(6.L,1.5L)*(findeventsin->param.totalMass*LAL_MTSUN_SI)*LAL_PI);
            chisqin.findEventsIn = findeventsin; 
            chisqin.filter1      = filter1; 
            chisqin.filter2      = filter2; 
            chisqin.rho1         = output1;
            chisqin.rho2         = output2;

            chisqout = (InspiralChisqOutput *) LALCalloc (1, sizeof(InspiralChisqOutput));
            chisqout->chisqZERO    = NULL; chisqout->chisqZERO    = (REAL8 *) LALCalloc (1, sizeof(REAL8) * output1.length);
            chisqout->chisqPIbyTWO = NULL; chisqout->chisqPIbyTWO = (REAL8 *) LALCalloc (1, sizeof(REAL8) * output1.length);
            chisqout->chisq        = NULL; chisqout->chisq        = (REAL8 *) LALCalloc (1, sizeof(REAL8) * output1.length);

            LALEvaluateInspiralChisqTest (status->statusPtr, chisqout, &chisqin); 
            CHECKSTATUSPTR (status); 
        } /* If chisqBins > 0 */
        
        /* This is where we populate the event list */
        for (i=0; i<nTempEvents; i++) {

            j = tempEventsIdx[i];

            if (!(*eventlist = (DiscoverInspiralEventsList*) LALRealloc(*eventlist, sizeof(DiscoverInspiralEventsList)*((*nEvents)+1))))
            {
                ABORT(status, LALNOISEMODELSH_EMEM, LALNOISEMODELSH_MSGEMEM);
            }

            (*eventlist)[*nEvents].snr             = sqrt( output1.data[j]*output1.data[j] + output2.data[j]*output2.data[j] );
            (*eventlist)[*nEvents].cmax1           = output1.data[j];
            (*eventlist)[*nEvents].cmax2           = output2.data[j];
            (*eventlist)[*nEvents].t0              = findeventsin->param.t0;
            (*eventlist)[*nEvents].t3              = findeventsin->param.t3;
            (*eventlist)[*nEvents].m1              = findeventsin->param.mass1;
            (*eventlist)[*nEvents].m2              = findeventsin->param.mass2;
            (*eventlist)[*nEvents].phase           = atan2(output2.data[j], output1.data[j]);
            (*eventlist)[*nEvents].effDistance     = distanceNorm / ((*eventlist)[*nEvents].snr);
            (*eventlist)[*nEvents].effD8           = distanceNorm / 8.0;

            (*eventlist)[*nEvents].amplitude       = 4.*(findeventsin->param.eta)
                    *(findeventsin->param.totalMass*LAL_MTSUN_SI/((*eventlist)[*nEvents].effD8))
                    *pow(LAL_PI*findeventsin->param.totalMass*LAL_MTSUN_SI*100.,2.L/3.L);;

            (*eventlist)[*nEvents].bin             = j;

            eSec = (double) j / findeventsin->param.tSampling;
            (*eventlist)[*nEvents].impulseTime     = findeventsin->currentGPSTime + (int) eSec;
            (*eventlist)[*nEvents].impulseTimeNS   = (int) (1.e9 * (eSec - (int) eSec));

            eSec += findeventsin->param.tC;
            (*eventlist)[*nEvents].endTime         = findeventsin->currentGPSTime + (int) eSec;
            (*eventlist)[*nEvents].endTimeNS       = (int) (1.e9 * (eSec - (int) eSec));

            (*eventlist)[*nEvents].sigmasq         = 0.0;

            /* If chisq was calculated we fill it accordingly. If chisq was
             * not calculated, we return zero
             * */
            if (findeventsin->chisqBins) {
                (*eventlist)[*nEvents].chisq           = chisqout->chisq[j];
                (*eventlist)[*nEvents].chisq1          = chisqout->chisqZERO[j];
                (*eventlist)[*nEvents].chisq2          = chisqout->chisqPIbyTWO[j];
            }
            else {
                (*eventlist)[*nEvents].chisq           = 0.0;
                (*eventlist)[*nEvents].chisq1          = 0.0;
                (*eventlist)[*nEvents].chisq2          = 0.0;
            }
            (*nEvents)++;
        } /* End of populating the event list */
       
        /* If we calculated chisq, we should free the temporary vectors */
        if (findeventsin->chisqBins) {
            LALFree (chisqout->chisqZERO);
            LALFree (chisqout->chisqPIbyTWO);
            LALFree (chisqout->chisq);
            LALFree (chisqout);
        }
    }

    if (tempEventsIdx != NULL) LALFree (tempEventsIdx);
    if (correlation.data != NULL) LALFree (correlation.data);
    if (output1.data != NULL) LALFree (output1.data);
    if (output2.data != NULL) LALFree (output2.data);
    if (filter1.data != NULL) LALFree (filter1.data);
    if (filter2.data != NULL) LALFree (filter2.data);

    DETATCHSTATUSPTR(status);
    RETURN(status);
}


/***************************************************
  THE FOLLOWING SUBROUTINE CHOSES TC SEPARATED EVENTS 

  Given a template, a valid inspiral event is one that 
  is separated from the others by at least one chirp 
  time. Input is the zero and pi/2 phase lag correlations. 
  Output is the number of valid events and a list of 
  their position (index).
 ****************************************************/

static void findValidInspiralEvents ( REAL4Vector c1, 
        REAL4Vector c2, 
        DiscoverInspiralEventsIn in, 
        INT4 *nEvents, 
        INT4 *eventIdx ) 
{
    INT4      i, events, i1, i2, j, n, unclassified;
    REAL4     x, y, z;
    REAL4     mVal;
    INT4      mPos, win;

    REAL4     *rho    = NULL;
    INT4      *Ones   = NULL, *bin = NULL;

    i1 = in.nBegin;
    i2 = in.signal.length - in.nEnd;

    n = 0;
    for (i=i1; i<=i2; i++) {
        x = c1.data[i];
        y = c2.data[i];
        z = sqrt(x*x + y*y);
        if (z >= in.Threshold)  n++;      
    }

    if (n>0) {

        rho  = (REAL4 *)LALCalloc(1,sizeof(REAL4)*n);
        bin  = (INT4 *)LALCalloc(1,sizeof(INT4)*n);
        Ones = (INT4 *)LALCalloc(1,sizeof(INT4)*n);

        j = 0;
        for (i=i1; i<=i2; i++) {

            x = c1.data[i];
            y = c2.data[i];
            z = sqrt(x*x + y*y);

            if (z >= in.Threshold)  {
                rho[j]  = z;
                bin[j]  = i;
                Ones[j] = 1;
                j++;      
            }
        } 

        /* We rely on the fact that the bin[] is naturally sorted */
        win    = (INT4)(in.param.tC * in.param.tSampling);
        events = 0;
        unclassified = computeSum(Ones,n);
        (*nEvents)   = 0;

        while (unclassified > 0) {

            findMax (rho, Ones, n, &mVal, &mPos);
            events ++;

            /*fprintf (stderr, "Found event %d at %d %f\n", events, bin[mPos], mVal);
	     */
            eventIdx[(*nEvents)] = bin[mPos];
            (*nEvents) ++;

            for (i=0; i<n; i++) { 
                if ( (bin[i] >= bin[mPos] - win) &&
                        (bin[i] <= bin[mPos] + win) ) Ones[i] = 0;     
            }
            unclassified = computeSum(Ones,n);
        }

        if ( rho != NULL ) LALFree (rho);
        if ( bin != NULL ) LALFree (bin);
        if ( Ones != NULL ) LALFree (Ones);
    }
}

/* ------
   Given a real vector vec and a int vector ok (containing binary values 0 and 1) 
   this routine returns the max (vec) over those indices where ok is non zero.
   The position of the maximum is also returned. 
   ------*/ 

static void findMax (REAL4 *vec, INT4 *ok, INT4 n, REAL4 *mVal, INT4 *mPos)
{
    REAL4 max=0.0;
    INT4  i, pos=0;

    for (i=0; i<n; i++) {
        if (ok[i] && vec[i]>max ){
            max = vec[i];
            pos = i;
        }
    }

    (*mVal) = max;
    (*mPos) = pos;
}

/* ------
   Returns the sum of values an integer vector 
   ------*/ 
static INT4 computeSum (INT4 *ivec, INT4 n)
{
    INT4 i;
    INT4 s;

    s = 0;
    for (i=0; i<n; i++) s += ivec[i];

    return s;
}

