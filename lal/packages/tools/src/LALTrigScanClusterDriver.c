/*
*  Copyright (C) 2007 Anand Sengupta, Craig Robinson
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

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
        const SnglInspiralTable **eventHead,
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

            condenseIn->masterList[k].clusterID = TRIGSCAN_UNCLASSIFIED;

            ++k;
        }
    }

    /*-- Call the appropriate subroutine to calculate metric and find clusters --*/
    switch (condenseIn->scanMethod)
    {
        case  T0T3Tc:
            if (condenseIn->vrbflag)
            {
                fprintf (stderr, "\t Setting massChoice t03. Will call "
                        "LALTrigScanFindStepT0T3 to calculate step size\n");
            }
            condenseIn->massChoice = t03;

            LALTrigScan3DMetricCoeff_TcT0T3 (status->statusPtr,
                    (const SnglInspiralTable **)(eventHead), condenseIn);
            CHECKSTATUSPTR (status);

            if ( condenseIn->vrbflag )
            {
                fprintf (stderr, "Max tc footprint = %e\n",
                        condenseIn->maxTcFootPrint);
            }

            /*-- Call the clustering driver --*/
            LALTrigScanClusterDriver( status->statusPtr,
                    condenseIn, &condenseOut, &nclusters );
            CHECKSTATUSPTR( status );

            /* Free the invGamma gsl_matrix in the masterList */
            {
                int k;

                for (k=0; k<condenseIn->n; k++) {
                    gsl_matrix_free ( condenseIn->masterList[k].invGamma );
                    condenseIn->masterList[k].invGamma = NULL;
                }
            }

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
            nclusters = 0;

            break;
    }

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

    if ( condenseOut )
          LALFree( condenseOut );

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

    n            = condenseIn->n;
    masterList   = condenseIn->masterList;

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
                            list, condenseIn, n, currClusterID, condenseOut
					)))
            {
                /* the seed point did not agglomerate into a cluster */
                /* So re-assign it as noise                          */
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

    /*--- How many clusters ---*/
    (*nclusters) = currClusterID-1;

    /* Append the stragglers if required */
    if ( condenseIn->appendStragglers )
    {
        LALTrigScanAppendIsolatedTriggers( status->statusPtr, condenseIn,
                condenseOut, nclusters);
        CHECKSTATUSPTR (status);
    }

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
        const SnglInspiralTable **eventHead,
        trigScanClusterIn       *condenseIn
        )
{
    INT4                N;
    trigScanInputPoint  *masterList=NULL;

    masterList = condenseIn->masterList;
    N = condenseIn->n;

    INITSTATUS (status, "LALTrigScan3DMetricCoeff_TcT0T3", LALTRIGSCANCLUSTERDRIVERC);
    ATTATCHSTATUSPTR(status);

    ASSERT (masterList, status, LALTRIGSCANCLUSTERH_ENULL, LALTRIGSCANCLUSTERH_MSGENULL);
    ASSERT (fabs((*eventHead)->Gamma[0]) > 0.0L, status, LALTRIGSCANCLUSTERH_ECHOICE, LALTRIGSCANCLUSTERH_MSGECHOICE);
    ASSERT (fabs((*eventHead)->Gamma[1]) > 0.0L, status, LALTRIGSCANCLUSTERH_ECHOICE, LALTRIGSCANCLUSTERH_MSGECHOICE);
    ASSERT (fabs((*eventHead)->Gamma[2]) > 0.0L, status, LALTRIGSCANCLUSTERH_ECHOICE, LALTRIGSCANCLUSTERH_MSGECHOICE);
    ASSERT (fabs((*eventHead)->Gamma[3]) > 0.0L, status, LALTRIGSCANCLUSTERH_ECHOICE, LALTRIGSCANCLUSTERH_MSGECHOICE);
    ASSERT (fabs((*eventHead)->Gamma[4]) > 0.0L, status, LALTRIGSCANCLUSTERH_ECHOICE, LALTRIGSCANCLUSTERH_MSGECHOICE);
    ASSERT (fabs((*eventHead)->Gamma[5]) > 0.0L, status, LALTRIGSCANCLUSTERH_ECHOICE, LALTRIGSCANCLUSTERH_MSGECHOICE);

    {
        const SnglInspiralTable *thisEvent = NULL;
        INT4  kk;
        REAL8 a11, a23, a22, a33, a12, a13, tcFootPrint, denom;

        if (condenseIn->vrbflag)
              fprintf (stderr, "Input trigger list seems to have the metric ... using it \n");

        /* assume that the maxFootPrint in the tc dimensions to be 0.0 */
        condenseIn->maxTcFootPrint = 0.0L;

        /* Now get the metric at that point */
        for (kk=0, thisEvent = (*eventHead);
                kk<N && thisEvent; kk++, thisEvent = thisEvent->next  )
        {
            /* Copy out the 3d metric co-effs */
            masterList[kk].Gamma[0] = thisEvent->Gamma[0]/condenseIn->ts_scaling;
            masterList[kk].Gamma[1] = thisEvent->Gamma[1]/condenseIn->ts_scaling;
            masterList[kk].Gamma[2] = thisEvent->Gamma[2]/condenseIn->ts_scaling;
            masterList[kk].Gamma[3] = thisEvent->Gamma[3]/condenseIn->ts_scaling;
            masterList[kk].Gamma[4] = thisEvent->Gamma[4]/condenseIn->ts_scaling;
            masterList[kk].Gamma[5] = thisEvent->Gamma[5]/condenseIn->ts_scaling;

            /* Now figure out the largest tc and tau0 footPrints of the triggers */
            a11 = masterList[kk].Gamma[0] ;
            a12 = masterList[kk].Gamma[1] ;
            a13 = masterList[kk].Gamma[2] ;
            a22 = masterList[kk].Gamma[3] ;
            a23 = masterList[kk].Gamma[4] ;
            a33 = masterList[kk].Gamma[5] ;

            tcFootPrint = (a23 * a23 - a22 * a33) * a22;
            denom = (a12*a23 - a22*a13) * (a12*a23 - a22*a13)
                    - (a23*a23 - a22*a33) * (a12*a12 - a22*a11);
            tcFootPrint = sqrt( tcFootPrint / denom );

            if ( tcFootPrint  > condenseIn->maxTcFootPrint )
                  condenseIn->maxTcFootPrint = tcFootPrint ;
        }

    }

    /* Fill-out the Gamma matrix */
    {
        int               kk, i, j, k, s1;
        gsl_permutation   *p1 = NULL;
        gsl_matrix        *GG = NULL;

        for (kk=0; kk<N; kk++)
        {

            p1 = NULL;
            GG = NULL;

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

            if (GG) gsl_matrix_free (GG);
            if (p1) gsl_permutation_free (p1);
        }
    }

    /* Normal exit */
    DETATCHSTATUSPTR(status);
    RETURN(status);
}

