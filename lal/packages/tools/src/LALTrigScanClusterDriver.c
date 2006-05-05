/*----------------------------------------------------------------------------
 *
 * File Name: LALTrigScanClusterDriver.c
 *
 * Author: Sengupta, Anand. S. and Gupchup, Jayant A.
 *
 * $Id$
 *
 *---------------------------------------------------------------------------*/

#if 0
<lalVerbatim file="LALTrigScanClusterDriverCV">
Author: Sengupta, Anand. S. and Gupchup, Jayant A.
$Id$"
</lalVerbatim>
#endif

#include <LALTrigScanCluster.h>

NRCSID (LALTRIGSCANCLUSTERDRIVERC,
  "$Id$");

/* Local function proto-type which returns the step length along         */
/* x and y directions (in this case tau0 and tau3 directions) and the    */ 
/* angle alpha between x-axis and the major eigendirection based         */
/* on the metric on the x-y space.                                       */
static void LALTrigScanFindStepT0T3 (
        LALStatus               *status,
        REAL8                   *x, 
        REAL8                   *y,
        trigScanClusterIn       *condenseIn, 
        REAL8                   *deltaX, 
        REAL8                   *deltaY,
        REAL8                   *alpha
        );

/*---------------------------------------------------------------------------
 * This wrapper function is called from outside world. It is responsible for
 * populating the masterList correctly and  calling the underlying driver 
 * function which actually carries out the trigScan clustering                  
 *--------------------------------------------------------------------------*/
void LALClusterSnglInspiralOverTemplatesAndEndTime ( 
        LALStatus              *status,
        SnglInspiralTable      **eventHead,
        trigScanClusterIn      *condenseIn
        )
{ 
    INT4                nclusters = 0;
    trigScanClusterOut  *condenseOut=NULL;
    SnglInspiralTable   *clusteredList = NULL;

    INITSTATUS (status, "LALClusterSnglInspiralOverTemplatesAndEndTime", 
            LALTRIGSCANCLUSTERDRIVERC);
    ATTATCHSTATUSPTR(status);

    ASSERT ((*eventHead), status, 
            LALTRIGSCANCLUSTERH_ENULL, LALTRIGSCANCLUSTERH_MSGENULL);
    ASSERT (condenseIn, status, 
            LALTRIGSCANCLUSTERH_ENULL, LALTRIGSCANCLUSTERH_MSGENULL);
    ASSERT (condenseIn->coarseShf.data, status, 
            LALTRIGSCANCLUSTERH_ENULL, LALTRIGSCANCLUSTERH_MSGENULL);
    ASSERT (condenseIn->n > 0, status, 
            LALTRIGSCANCLUSTERH_ENULL, LALTRIGSCANCLUSTERH_MSGENULL);
    ASSERT (condenseIn->scanMethod > 0 && condenseIn->scanMethod <=2, 
            status, LALTRIGSCANCLUSTERH_ENULL, LALTRIGSCANCLUSTERH_MSGENULL);

    if (condenseIn->vrbflag)
          fprintf(stderr, "\nInput contains %d points \n", condenseIn->n);

    condenseIn->masterList = (trigScanInputPoint *)
            LALCalloc (1, condenseIn->n * sizeof(trigScanInputPoint));

    /* Copy the relevant fields from the inspiralEventList */
    {
        SnglInspiralTable   *thisEvent = NULL;
        INT4 l=0;

        for ( thisEvent = (*eventHead); 
                thisEvent; thisEvent = thisEvent->next )
        {
            if (condenseIn->scanMethod == T0T3Tc)
            {
                condenseIn->masterList[l].x = thisEvent->tau0;
                condenseIn->masterList[l].y = thisEvent->tau3;
            }
            condenseIn->masterList[l].tc_sec = 
                    (thisEvent->end_time).gpsSeconds;
            condenseIn->masterList[l].tc_ns  = 
                    (thisEvent->end_time).gpsNanoSeconds;
            condenseIn->masterList[l].rho    = thisEvent->snr;
            condenseIn->masterList[l].effD   = thisEvent->eff_distance;
            condenseIn->masterList[l].chisq  = 0.0;
            condenseIn->masterList[l].isValidEvent = 1;

            ++l;
        }
    }

    /*-- Call the clustering driver --*/
    LALTrigScanClusterDriver( status->statusPtr, 
            condenseIn, &condenseOut, &nclusters );
    CHECKSTATUSPTR( status );

    if ( nclusters ) /* some clusters were found */
    {
        /* Clustering is now over. Create a new list containing the clustered
         * elements only. The index is present in condenseOut->master_idx
         */
        clusteredList = NULL;
        clusteredList = 
                XLALTrimSnglInspiralTable( eventHead, condenseOut, nclusters );

        /* Now del the inspiralEventList and replace it with clusteredList */
        if ( clusteredList )
        {
            XLALDeleteSnglInspiralTable(eventHead);
            if (!(*eventHead) && nclusters >0 )
            {
                (*eventHead) = clusteredList;
            }
        }
    }
    else /* No clusters were found by our algorithm */
    {
        /* In this case delete the original list. */
        XLALDeleteSnglInspiralTable( eventHead );
    }

    /* Free the masterlist */
    if ( condenseIn->masterList )
          LALFree( condenseIn->masterList );

    /* Normal Exit */
    DETATCHSTATUSPTR(status);
    RETURN(status);
}

/*------------------------------------------------------------------------*/
/* Function which clusters points in three-D (t0,t3 and tc) using the key */
/* TrigScanCluster algorithm                                              */
/*------------------------------------------------------------------------*/
void LALTrigScanClusterDriver (
        LALStatus            *status,
        trigScanClusterIn    *condenseIn,
        trigScanClusterOut   **condenseOut,
        INT4                 *nclusters
        )
{
    INT4    i, n;
    REAL8   dt_a, *nx, *ny, *nstepX, *nstepY, *nalpha;

    /* For clustering */
    trigScanInputPoint     *masterList=NULL, *list=NULL;
    ExpandClusterInput   expandClusterIn;  

    INITSTATUS (status, 
            "LALTrigScanClusterDriver.c", LALTRIGSCANCLUSTERDRIVERC);
    ATTATCHSTATUSPTR(status);

    ASSERT (condenseIn, 
            status, LALTRIGSCANCLUSTERH_ENULL, LALTRIGSCANCLUSTERH_MSGENULL);
    ASSERT (condenseIn->masterList, 
            status, LALTRIGSCANCLUSTERH_ENULL, LALTRIGSCANCLUSTERH_MSGENULL);
    ASSERT (condenseIn->n > 0, status, 
            LALTRIGSCANCLUSTERH_ECHOICE, LALTRIGSCANCLUSTERH_MSGECHOICE);
    ASSERT (condenseIn->appendStragglers >= 0 && 
            condenseIn->appendStragglers <= 1, status, 
            LALTRIGSCANCLUSTERH_ECHOICE, LALTRIGSCANCLUSTERH_MSGECHOICE);
    ASSERT (condenseIn->coarseShf.data->data, 
            status, LALTRIGSCANCLUSTERH_ENULL, LALTRIGSCANCLUSTERH_MSGENULL);
    ASSERT ((condenseIn->coarseShf.data->length > 1) && 
            (condenseIn->coarseShf.deltaF > 0.0), status, 
            LALTRIGSCANCLUSTERH_ECHOICE, LALTRIGSCANCLUSTERH_MSGECHOICE);

    n            = condenseIn->n;
    masterList   = condenseIn->masterList;
    dt_a         = 0.5/ (condenseIn->coarseShf.data->length - 1);
    dt_a        /= condenseIn->coarseShf.deltaF;

    /*-- Create space to store the putative events into temp list --*/
    nx = (REAL8 *) LALMalloc( n * sizeof(REAL8));
    ny = (REAL8 *) LALMalloc( n * sizeof(REAL8));

    /*-- Copy the temp list of x and y co-ords and init clusterID --*/
    for (i = 0; i < n; i++)
    {
        nx[i] = masterList[i].x;
        ny[i] = masterList[i].y;

        masterList[i].clusterID = TRIGSCAN_UNCLASSIFIED; 
    }

    /*-- The stepSize to be calculated at all the points in x-y --*/
    nstepX = (REAL8 *) LALCalloc( 1, n * sizeof(REAL8));
    nstepY = (REAL8 *) LALCalloc( 1, n * sizeof(REAL8));
    nalpha = (REAL8 *) LALCalloc( 1, n * sizeof(REAL8));

    /*-- Calculate the step size now by calling the appropriate subroutine --*/
    switch (condenseIn->scanMethod) 
    { 
        case  T0T3Tc:
            if (condenseIn->vrbflag)
            {
                fprintf (stderr, "\t Setting massChoice t03. Will call "
                        "LALTrigScanFindStepT0T3 to calculate step size\n");
            }
            condenseIn->massChoice = t03;
            LALTrigScanFindStepT0T3 (status->statusPtr, nx , ny, condenseIn, 
                    nstepX, nstepY, nalpha);
            CHECKSTATUSPTR (status);
            break;

        case Psi0Psi3Tc:
        case trigScanNone:
        default:
            /* Print a warning message */
            LALWarning( status->statusPtr, 
                    "trigScanCluster is not available for "
                    "this choice of scanMethod." );
            CHECKSTATUSPTR (status);
            /* Reset nclusters to 0 */
            (*nclusters) = 0;
            /* Return control */
            DETATCHSTATUSPTR(status);
            RETURN(status);
            break;
    }

    condenseIn->a     = nstepX;
    condenseIn->b     = nstepY;
    condenseIn->theta = nalpha;

    /*-----------------------*/
    /*-- CLUSTERING BEGINS --*/
    /*-----------------------*/

    expandClusterIn.nInputPoints  = n;
    expandClusterIn.currClusterID = 1;

    /*-- For each point in the masterList --*/
    for (i = 0; i < n; i++)
    {
        /* if the point is'nt a part of a cluster */
        if (masterList[i].clusterID == TRIGSCAN_UNCLASSIFIED)
        {
            /* create a list of size 1 */
            list = (trigScanInputPoint *) 
                    LALMalloc (sizeof(trigScanInputPoint));

            /* assume that this seed point is a cluster */
            masterList[i].clusterID = expandClusterIn.currClusterID;

            /* assign the values of the seed-point to the list */
            list[0].x            =  masterList[i].x;
            list[0].y            =  masterList[i].y;
            list[0].tc_sec       =  masterList[i].tc_sec;
            list[0].tc_ns        =  masterList[i].tc_ns;
            list[0].rho          =  masterList[i].rho;
            list[0].isValidEvent =  masterList[i].isValidEvent;
            list[0].clusterID    =  masterList[i].clusterID;

            expandClusterIn.epsTc  = condenseIn->bin_time * dt_a;  
            expandClusterIn.epsX   = sqrt(condenseIn->sf_area) * nstepX[i];  
            expandClusterIn.epsY   = sqrt(condenseIn->sf_area) * nstepY[i];  
            expandClusterIn.alpha  = nalpha[i];

            /*-- Try to expand this seed and see if it agglomerates more --*/
            /*-- more points around it. If it does not then assign it to --*/
            /*-- noise. Otherwise increment the current clusterID for    --*/
            /*-- the next cluster.                                       --*/

            if (!(XLALTrigSCanExpandCluster ( list, masterList, 
                            expandClusterIn))
                    )
            {
                /* the seed point did not agglomerate into a cluster */
                masterList[i].clusterID = TRIGSCAN_NOISE;
            }
            else
            {
                /* a cluster has been found .. eureka !*/
                expandClusterIn.currClusterID++;
            }
        }
    }

    /*---------------------*/
    /*-- CLUSTERING ENDS --*/
    /*---------------------*/

    /*--- Output ----------*/
    (*nclusters) = expandClusterIn.currClusterID-1;  

    /* Worry about output only if there are clusters */
    if (*nclusters)
    {
        /*-- Allocate memory for output --*/
        if ( !(*condenseOut = (trigScanClusterOut*) 
                    LALRealloc(*condenseOut, 
                        sizeof(trigScanClusterOut)*(*nclusters)))
           )
        {
            fprintf (stderr, "LALRealloc error in condenseout. \n"); 
            abort ();
        }

        LALTrigScanClusterMakeOutput ( status->statusPtr, 
                condenseIn, condenseOut, (*nclusters));
        CHECKSTATUSPTR (status);

    } /* If nclusters */

    if ( condenseIn->appendStragglers )
    {
        /* Append the unclustered points */
        LALTrigScanAppendIsolatedTriggers( status->statusPtr, 
                condenseIn, condenseOut, nclusters);
        CHECKSTATUSPTR (status);
    }

    /* Free-up memory */
    if (nx) LALFree (nx);
    if (ny) LALFree (ny);
    if (nstepX) LALFree (nstepX);
    if (nstepY) LALFree (nstepY);
    if (nalpha) LALFree (nalpha);

    /* Normal exit */
    DETATCHSTATUSPTR(status);
    RETURN(status);
}





/*-------------------------------------------------------------------------*/
/* Local function that calculates the step size along tau0/tau3 directions */
/* The convention we follow is x[] = tau0[], y[] = tau3[].                 */ 
/* For each of these points we calculate deltaX, deltaY and alpha          */
/*-------------------------------------------------------------------------*/
static void LALTrigScanFindStepT0T3 (
        LALStatus               *status,
        REAL8                   *x, 
        REAL8                   *y,
        trigScanClusterIn       *condenseIn, 
        REAL8                   *deltaX, 
        REAL8                   *deltaY,
        REAL8                   *alpha
        )
{
    INT4                N;
    REAL8               dx, dy;     
    InspiralMetric      *metric   = NULL;
    InspiralTemplate    *tempPars = NULL;
    InspiralMomentsEtc  *moments  = NULL;

    N = condenseIn->n;

    INITSTATUS (status, "LALTrigScanFindStepT0T3", LALTRIGSCANCLUSTERDRIVERC);
    ATTATCHSTATUSPTR(status);

    ASSERT (x && y && deltaX && deltaY && alpha, 
            status, LALTRIGSCANCLUSTERH_ENULL, LALTRIGSCANCLUSTERH_MSGENULL);

    tempPars = (InspiralTemplate *) LALMalloc (sizeof (InspiralTemplate));
    tempPars->nStartPad   = 0;
    tempPars->nEndPad     = 0;
    tempPars->startPhase  = 0.0;
    tempPars->t0          = x[0];
    tempPars->t3          = y[0];
    tempPars->fLower      = condenseIn->fLower;
    tempPars->fCutoff     = condenseIn->fUpper;
    tempPars->tSampling   = condenseIn->tSampling;
    tempPars->massChoice  = condenseIn->massChoice;
    tempPars->order       = condenseIn->order;
    tempPars->approximant = condenseIn->approximant;

    LALInspiralParameterCalc (status->statusPtr, tempPars);
    CHECKSTATUSPTR (status);

    /* Get the noise moments */
    moments = (InspiralMomentsEtc *) 
            LALCalloc (1, sizeof (InspiralMomentsEtc));

    LALGetInspiralMoments (status->statusPtr, 
            moments, &condenseIn->coarseShf, tempPars);
    CHECKSTATUSPTR (status);

    /* Now get the metric at that point */
    {
        int kk;

        for (kk=0; kk<N; kk++) 
        {
            tempPars->t0       = x[kk];
            tempPars->t3       = y[kk];
            LALInspiralParameterCalc (status->statusPtr, tempPars);
            CHECKSTATUSPTR (status);

            metric = (InspiralMetric *) 
                    LALCalloc (1, sizeof (InspiralMetric));

            LALInspiralComputeMetric (status->statusPtr, 
                    metric, tempPars, moments);
            CHECKSTATUSPTR (status);

            /* Now calculate the step size in X and Y direction. 
             * Remember that these are approximate.  However, this 
             * is the best we can do. */
            dx = sqrt(2.0*(1. - condenseIn->mmCoarse)/metric->g00);
            dy = sqrt(2.0*(1. - condenseIn->mmCoarse)/metric->g11);

            deltaX[kk] = dx;
            deltaY[kk] = dy;
            alpha[kk]  = metric->theta;
        }
    }

    if (tempPars)   LALFree (tempPars);
    if (moments)    LALFree (moments);
    if (metric)     LALFree (metric);

    /* Normal exit */
    DETATCHSTATUSPTR(status);
    RETURN(status);
}

