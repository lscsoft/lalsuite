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
$Id$
</lalVerbatim>
#endif

#include <lal/LALTrigScanCluster.h>

NRCSID (LALTRIGSCANCLUSTERDRIVERC,
  "$Id$");

/* Local function prototype. This function calculates the 3d metric in the
 * space of tc,tau0,tau3 for every trigger. Note that y direction is \tau_0
 * and z direction corresponds to \tau_3. Also, the metric is scaled down by
 * S(1-MM) where S is a safety factor and MM is the minimal match of the
 * template bank used in the search.
 */
static void LALTrigScan3DMetricCoeff_TcT0T3  (
        LALStatus               *status,
        REAL8                   *y, 
        REAL8                   *z,
        trigScanClusterIn       *condenseIn 
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

    /* Allocate memory for the masterList */
    condenseIn->masterList = (trigScanInputPoint *)
            LALCalloc (1, condenseIn->n * sizeof(trigScanInputPoint));

    /* Copy the relevant fields from the inspiralEventList to masterList */
    {
        SnglInspiralTable   *thisEvent = NULL;
        INT4 k=0;

        for ( thisEvent = (*eventHead); 
                thisEvent; thisEvent = thisEvent->next )
        {
            if (condenseIn->scanMethod == T0T3Tc)
            {
                condenseIn->masterList[k].y = thisEvent->tau0;
                condenseIn->masterList[k].z = thisEvent->tau3;
            }
            condenseIn->masterList[k].tc_sec = 
                    (thisEvent->end_time).gpsSeconds;
            condenseIn->masterList[k].tc_ns  = 
                    (thisEvent->end_time).gpsNanoSeconds;
            condenseIn->masterList[k].rho    = thisEvent->snr;

            ++k;
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

    /* Free the invGamma gsl_matrix in the masterList */
    {
        int k;

        for (k=0; k<condenseIn->n; k++) {
              gsl_matrix_free ( condenseIn->masterList[k].invGamma );
              condenseIn->masterList[k].invGamma = NULL;
        }
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
    REAL8   *ny, *nz;

    /* For clustering */
    trigScanInputPoint   *masterList=NULL;
    INT4                 *list = NULL;
    INT4                 currClusterID = 1;

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

    /*-- Create space to store the putative events into temp list --*/
    ny = (REAL8 *) LALMalloc( n * sizeof(REAL8));
    nz = (REAL8 *) LALMalloc( n * sizeof(REAL8));

    /*-- Copy the temp list of y and z co-ords and init clusterID --*/
    for (i = 0; i < n; i++)
    {
        ny[i] = masterList[i].y;
        nz[i] = masterList[i].z;

        masterList[i].clusterID = TRIGSCAN_UNCLASSIFIED; 
    }

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

            LALTrigScan3DMetricCoeff_TcT0T3 (status->statusPtr, ny , nz, condenseIn);
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

    /*-----------------------*/
    /*-- CLUSTERING BEGINS --*/
    /*-----------------------*/

    currClusterID = 1;

    /*-- For each point in the masterList --*/
    for (i = 0; i < n; i++)
    {
        /* if the point is'nt a part of a cluster */
        if (masterList[i].clusterID == TRIGSCAN_UNCLASSIFIED)
        {
            /* create a list of size 1 */
            list = (INT4 *) 
                    LALMalloc (sizeof(INT4));

            /* assume that this seed point is a cluster */
            masterList[i].clusterID = currClusterID;

            /* assign the values of the seed-point to the list */
            list[0]  =  i;

            /*-- Try to expand this seed and see if it agglomerates more --*/
            /*-- more points around it. If it does not then assign it to --*/
            /*-- noise. Otherwise increment the current clusterID for    --*/
            /*-- the next cluster.                                       --*/

            if (!(XLALTrigScanExpandCluster ( 
                            list, masterList, n, currClusterID 
                          )))
           {
                /* the seed point did not agglomerate into a cluster */
                masterList[i].clusterID = TRIGSCAN_NOISE;
            }
            else
            {
                /* a cluster has been found .. eureka !*/
                currClusterID++;
            }
        }
    }

    /*---------------------*/
    /*-- CLUSTERING ENDS --*/
    /*---------------------*/

    /*--- Output ----------*/
    (*nclusters) = currClusterID-1;  

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

    } 

    /* Append the unclustered points if required */
    if ( condenseIn->appendStragglers )
    {
        LALTrigScanAppendIsolatedTriggers( status->statusPtr, 
                condenseIn, condenseOut, nclusters);
        CHECKSTATUSPTR (status);
    }

    /* Free-up memory */
    if (ny) LALFree (ny);
    if (nz) LALFree (nz);

    /* Normal exit */
    DETATCHSTATUSPTR(status);
    RETURN(status);
}

/*-------------------------------------------------------------------------*/
/* Local function that calculates the 3d metric                            */
/* The convention we follow is x[] = tc, y[] = tau0[], z[] = tau3[].       */
/* Note that the metric is scaled down by S(1-MM)                          */
/*-------------------------------------------------------------------------*/
static void LALTrigScan3DMetricCoeff_TcT0T3  (
        LALStatus               *status,
        REAL8                   *y, 
        REAL8                   *z,
        trigScanClusterIn       *condenseIn 
        )
{
    INT4                N;
    InspiralMetric      *metric   = NULL;
    InspiralTemplate    *tempPars = NULL;
    InspiralMomentsEtc  *moments  = NULL;
    trigScanInputPoint  *masterList=NULL;

    masterList = condenseIn->masterList;
    N = condenseIn->n;

    INITSTATUS (status, "LALTrigScan3DMetricCoeff_TcT0T3", LALTRIGSCANCLUSTERDRIVERC);
    ATTATCHSTATUSPTR(status);

    ASSERT (y && z && masterList, 
            status, LALTRIGSCANCLUSTERH_ENULL, LALTRIGSCANCLUSTERH_MSGENULL);

    tempPars = (InspiralTemplate *) LALMalloc (sizeof (InspiralTemplate));
    tempPars->nStartPad   = 0;
    tempPars->nEndPad     = 0;
    tempPars->startPhase  = 0.0;
    tempPars->t0          = y[0];
    tempPars->t3          = z[0];
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
            tempPars->t0       = y[kk];
            tempPars->t3       = z[kk];
            LALInspiralParameterCalc (status->statusPtr, tempPars);
            CHECKSTATUSPTR (status);

            metric = (InspiralMetric *) 
                    LALCalloc (1, sizeof (InspiralMetric));

            LALInspiralComputeMetric (status->statusPtr, 
                    metric, tempPars, moments);
            CHECKSTATUSPTR (status);

            /* Copy out the 3d metric co-effs */
            masterList[kk].Gamma[0] = metric->Gamma[0]/(condenseIn->sf_volume*(1.0L - condenseIn->mmCoarse));
            masterList[kk].Gamma[1] = metric->Gamma[1]/(condenseIn->sf_volume*(1.0L - condenseIn->mmCoarse));
            masterList[kk].Gamma[2] = metric->Gamma[2]/(condenseIn->sf_volume*(1.0L - condenseIn->mmCoarse));
            masterList[kk].Gamma[3] = metric->Gamma[3]/(condenseIn->sf_volume*(1.0L - condenseIn->mmCoarse));
            masterList[kk].Gamma[4] = metric->Gamma[4]/(condenseIn->sf_volume*(1.0L - condenseIn->mmCoarse));
            masterList[kk].Gamma[5] = metric->Gamma[5]/(condenseIn->sf_volume*(1.0L - condenseIn->mmCoarse));

            /* This is the part where we calculate the inverse of the Gamma
             * matrix and store it in the element of the masterList structure.
             */
            {
                int               i, j, k, s1;
                gsl_permutation   *p1 = NULL;
                gsl_matrix        *GG = NULL;

                GG = gsl_matrix_calloc (3,3);

                k = 0;
                for (i=0; i<3; i++) {
                      for(j=i; j<3; j++) {
                          gsl_matrix_set (GG, i, j, masterList[kk].Gamma[k]);
                          gsl_matrix_set (GG, j, i, gsl_matrix_get (GG, i, j));
                          k = k+1;
                      }
                }
                            
                p1 = gsl_permutation_alloc( 3 );
                gsl_linalg_LU_decomp( GG, p1, &s1 ); 
                
                masterList[kk].invGamma = NULL;
                masterList[kk].invGamma = gsl_matrix_alloc( 3, 3 );
                gsl_linalg_LU_invert( GG, p1, masterList[kk].invGamma );

                gsl_matrix_free (GG);
                gsl_permutation_free (p1);
            } 
        }
    }

    if (tempPars)   LALFree (tempPars);
    if (moments)    LALFree (moments);
    if (metric)     LALFree (metric);

    /* Normal exit */
    DETATCHSTATUSPTR(status);
    RETURN(status);
}

