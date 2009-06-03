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
 * File Name: LALTrigScanCluster.c
 *
 * Author: Sengupta, Anand. S. and Gupchup, Jayant A.
 *
 * $Id$
 *
 *---------------------------------------------------------------------------*/

#if 0
<lalVerbatim file="LALTrigScanClusterCV">
Author: Sengupta, Anand. S. and Gupchup, Jayant A.
$Id$
</lalVerbatim>
#endif

#include <lal/LALTrigScanCluster.h>

NRCSID (LALTRIGSCANCLUSTERC,
        "$Id$");

/* ------------------------------------------------------------------*/
/* This function takes in a seed point and 'expands' that seed point */
/* to include other points from a masterList which are part of the   */
/* same cluster. If the input seed cannot be 'expand'ed then it      */
/* returns false. Otherwise it returns true.                         */
/* ------------------------------------------------------------------*/
trigScanValidEvent XLALTrigScanExpandCluster (
        INT4                  *list,
        trigScanClusterIn     *condenseIn,
        INT4                  nPoints,
        INT4                  currClusterID,
        trigScanClusterOut    **condenseOut
        )
{
    static LALStatus        status;
    trigScanValidEvent      flag = trigScanFalse;
    INT4                    seed;
    trigScanEpsSearchInput  epsSearchIn; /* Data structure given as input to*/
                                         /* getEpsNeighbourhood fn.         */

    INT4   pointer=0, size=1; /* when pointer < size, all elements */
                              /* have been accessed                */

    /* Create the epsSearchIn data structure */
    epsSearchIn.masterList     = condenseIn->masterList;
    epsSearchIn.nInputPoints   = nPoints;
    epsSearchIn.clusterID      = currClusterID;
    epsSearchIn.maxTcFootPrint = condenseIn->maxTcFootPrint;
    epsSearchIn.minLoopIdx     = list[0];

    while (pointer < size) {

        /* Pick the index of the first point in this list as the seed */
        seed = list[pointer];

        /* call function which returns the points which are inside the */
        /* eps-Neighbourhood. Allow the list to grow as update size as */
        /* more seeds are added to the list                            */
        XLALTrigScanGetEpsNeighbourhood (seed, &list, &size, &epsSearchIn);

        /* if valid seed then insert and update size by */
        /* number of points retrieved                   */

        pointer++;
    }

    /* set flag to true if (size > 1) indicating that a valid cluster
     * has been discovered
     * */
    if (size > 1)
    {


        LALTrigScanStoreThisCluster (
                &status, (const INT4 *)( list ),
                (const trigScanClusterIn  *)( condenseIn ),
                (const INT4)( size ), (const INT4)( currClusterID ),
                condenseOut
                );

        flag = trigScanTrue;

    }

    /* deallocate the list before returning */
    LALFree (list);
    list = NULL;

    return flag;
}

/* ------------------------------------------------------------------------*/
/* The following function takes in a seed point and a list of other points */
/* in the master list. Each unclassified point in the masterList is a      */
/* possible candidate for 'epsilon neighbour' of this seed point.          */
/* What it does is the following :                                         */
/*  For each unclassified point in the masterList                          */
/*      Calculate if it is an eps-neighbor of seed point.                  */
/*  If true, that unclassified point in the masterList and the seed        */
/*  become part of the same cluster. Further, the size of the temporary    */
/*  list of epsilon neighbors is incremented by 1 - and this brand-new     */
/*  member added to it. The 'size' variable is the length of this list.    */
/* ------------------------------------------------------------------------*/
void XLALTrigScanGetEpsNeighbourhood (
        INT4                    seed,
        INT4                    **list,
        INT4                    *size,
        trigScanEpsSearchInput  *epsSearchIn
        )
{

    static const char *func = "XLALTrigScanGetEpsNeighbourhood";

    INT4              i;
    REAL8             distance;
    REAL8             q1[3], q2[3]; /* Position vectors */
    gsl_vector_view   vq1, vq2;     /* vector views from position vectors */
    fContactWorkSpace *workSpace;   /* reqd for ellipsoid overlap */

    /* Init the workSpace required for checking ellipsoid overlaps */
    workSpace = XLALInitFContactWorkSpace( 3, NULL, NULL, gsl_min_fminimizer_brent, 1.0e-2 );

    /* Set the position vector (q1) of the seed point */
    q1[0] = epsSearchIn->masterList[seed].tc_sec
            + 1.e-9*(epsSearchIn->masterList[seed].tc_ns);
    q1[1] = epsSearchIn->masterList[seed].y;
    q1[2] = epsSearchIn->masterList[seed].z;

    /* create a vector view from the q1 array */
    vq1 = gsl_vector_view_array(q1,3);

    /* Set the shape matrix of the seed point */
    workSpace->invQ1  = epsSearchIn->masterList[seed].invGamma;

    for (i = epsSearchIn->minLoopIdx; i < epsSearchIn->nInputPoints; i++)
    {
        /* check if the point has been classified */
        if (epsSearchIn->masterList[i].clusterID == TRIGSCAN_UNCLASSIFIED)
        {
            /* Set the position vector (q2) of the i-th point */
            q2[0] = epsSearchIn->masterList[i].tc_sec
                    + 1.e-9*(epsSearchIn->masterList[i].tc_ns);
            q2[1] = epsSearchIn->masterList[i].y;
            q2[2] = epsSearchIn->masterList[i].z;

            /* Notice that the triggers are actually time ordered. This means
             * So, if the endTime difference between trigger and seed exceeds
             * some value (say 6 * maxTcFootPrint), quit the loop
             */
            if ( ( (q2[0] - q1[0]) >= 6.0 * epsSearchIn->maxTcFootPrint ) )
            {
                /* Note that the order q2 - q1 is VERY important in the above
                 * test. This tests for triggers that are further ahead in time
                 * than the seed are not too far away. If they are, there is no
                 * need to continue the loop any further.
                 fprintf (stderr, "Breaking out as q2 (%d) is too far from q1 (%d) index = %d \n", i, seed, i);
                 */
                break;
            }


            /* Worry about overlaps only if the triggers times differ by a
             * fixed amount. If the trigger times differ by more than this, it
             * is assumed that the overlap will fail anyway
             */
            if ( (fabs(q2[0] - q1[0]) <= 6.0 * epsSearchIn->maxTcFootPrint ) )
            {
                /* create a vector view from the q2 array */
                vq2 = gsl_vector_view_array(q2,3);

                /* Set the shape matrix of the i-th point */
                workSpace->invQ2    = epsSearchIn->masterList[i].invGamma;

                /* Figure out if the above ellipsoids overlap */
                distance = XLALCheckOverlapOfEllipsoids (&(vq1.vector), &(vq2.vector), workSpace);
                if ( XLAL_IS_REAL8_FAIL_NAN( distance ) )
                {
                  XLALFreeFContactWorkSpace( workSpace );
                  XLAL_ERROR_VOID( func, XLAL_EFUNC );
                }


                if ( distance <= 1.)
                {
                    /* set the clusterID to the currClusterID */
                    epsSearchIn->masterList[i].clusterID = epsSearchIn->clusterID;

                    /* increment the size variable and hence realloc the list */
                    (*size)++;

                    if ( !(*list = (INT4*)
                                LALRealloc(*list,
                                    sizeof(INT4)*(*size)))
                       )
                    {
                        fprintf (stderr, "LALRealloc error. Aborting at %d\n",
                                *size);
                        abort ();
                    }

                    /* add the shortlisted point to the list at position size - 1*/
                    (*list)[*size - 1]  = i;

                } /* If the two triggers overlap */

            } /* If two triggers are within n times maxTcFootPrint only then */

        } /*if unclassified trigger */

    } /* Loop over triggers */


    /* De-allocate workSpace allocated earlier */
    /* This workSpace is used to check overlap */
    /* of ambiguity ellipsoids                 */
    XLALFreeFContactWorkSpace( workSpace );
}

/*----------------------------------------------------------------------
 * Once a cluster has been discovered we need to store it in the
 * trigScanClusterOut structure. This particular function does that job. The
 * inputs are the list of triggers in the newly discovered cluster and its
 * size, the corresponding clusterID.
 ---------------------------------------------------------------------*/
void LALTrigScanStoreThisCluster (
        LALStatus                *status,
        const INT4               *list,
        const trigScanClusterIn  *condenseIn,
        const INT4               size,
        const INT4               currClusterID,
        trigScanClusterOut       **condenseOut
        )
{
    REAL4Vector *unSortedSnr = NULL;
    INT4Vector  *heapSortIndex = NULL;
    INT4        i, mid;

    INITSTATUS (status, "LALTrigScanStoreThisEvent", LALTRIGSCANCLUSTERC);
    ATTATCHSTATUSPTR(status);

    /* We can now figure out the maximum SNR here */
    unSortedSnr    = XLALCreateREAL4Vector( size );
    heapSortIndex  = XLALCreateINT4Vector( size );

    for (i=0; i<size; i++)
    {
        mid = list[i];
        unSortedSnr->data[i] = condenseIn->masterList[mid].rho;
    }

    LALSHeapIndex(status->statusPtr, heapSortIndex, unSortedSnr);
    CHECKSTATUSPTR( status );

    mid = list[heapSortIndex->data[size-1]];

    /*-- Allocate memory for output --*/
    if ( !(*condenseOut = (trigScanClusterOut*)
                LALRealloc(*condenseOut,
                    sizeof(trigScanClusterOut)*(currClusterID)))
       )
    {
        fprintf (stderr, "LALRealloc error in condenseout. \n");
        abort ();
    }

    (*condenseOut)[currClusterID-1].y          = condenseIn->masterList[mid].y;
    (*condenseOut)[currClusterID-1].z          = condenseIn->masterList[mid].z;
    (*condenseOut)[currClusterID-1].tc_sec     = condenseIn->masterList[mid].tc_sec;
    (*condenseOut)[currClusterID-1].tc_ns      = condenseIn->masterList[mid].tc_ns;
    (*condenseOut)[currClusterID-1].rho        = condenseIn->masterList[mid].rho;
    (*condenseOut)[currClusterID-1].master_idx = mid;
    (*condenseOut)[currClusterID-1].nelements  = size;

    /* Delete the Vectors used for heap sort */
    XLALDestroyREAL4Vector( unSortedSnr );
    XLALDestroyINT4Vector( heapSortIndex );

    /* Normal exit */
    DETATCHSTATUSPTR(status);
    RETURN(status);
}

/* ---------------------------------------------------------------------
 * This function appends the stragglers to the condenseOut list
 ----------------------------------------------------------------------*/
void LALTrigScanAppendIsolatedTriggers (
        LALStatus               *status,
        trigScanClusterIn       *condenseIn,
        trigScanClusterOut      **condenseOut,
        INT4                    *nclusters
        )
{
    INT4                 i, j, n, ni, n1;
    REAL8                xi, xj, dxij;
    trigScanInputPoint   *masterList=NULL;
    REAL8                *xx=NULL, *vv=NULL;
    INT4                 *mid = NULL;

    INITSTATUS (status,
            "LALTrigScanAppendIsolatedTriggers", LALTRIGSCANCLUSTERC);
    ATTATCHSTATUSPTR(status);

    n            = condenseIn->n;
    masterList   = condenseIn->masterList;

    /* AT first figure out how many stragglers */
    for (i=1, n1=0; i<=n; i++) {
        if (masterList[i-1].clusterID < 1 ) {
            n1++;
        }
    }

    /* Now allocate memory */
    mid  = LALCalloc (1, (n1+1)*sizeof(INT4));
    xx   = LALCalloc (1, (n1+1)*sizeof(REAL8));
    vv   = LALCalloc (1, (n1+1)*sizeof(REAL8));

    ASSERT (condenseOut,
            status, LALTRIGSCANCLUSTERH_ENULL, LALTRIGSCANCLUSTERH_MSGENULL);
    ASSERT (masterList,
            status, LALTRIGSCANCLUSTERH_ENULL, LALTRIGSCANCLUSTERH_MSGENULL);
    ASSERT (n > 0, status,
            LALTRIGSCANCLUSTERH_ECHOICE, LALTRIGSCANCLUSTERH_MSGECHOICE);
    ASSERT (*nclusters >= 0, status,
            LALTRIGSCANCLUSTERH_ECHOICE, LALTRIGSCANCLUSTERH_MSGECHOICE);
    ASSERT (xx && vv && mid,
            status, LALTRIGSCANCLUSTERH_ENULL, LALTRIGSCANCLUSTERH_MSGENULL);

#if 0
    if ( condenseIn->vrbflag )
          fprintf (stderr, "--------- BEGINNING TO APPEND ----------\n");
#endif
    for (i=1, n1=0; i<=n; i++) {
        if (masterList[i-1].clusterID < 1 ) {
            n1++;
            xx[n1]  = masterList[i-1].tc_sec + masterList[i-1].tc_ns/1.0e9;
            vv[n1]  = masterList[i-1].rho;
            mid[n1] = i-1;

#if 0
            if ( condenseIn->vrbflag )
                  fprintf (stderr, "%4d   %d  %d   %d   %d   %1.12e\n",
                          n1, mid[n1],
                          (masterList[i-1].clusterID),
                          (int)(masterList[i-1].tc_sec),
                          (int)(masterList[i-1].tc_ns),
                          masterList[i-1].rho);
#endif
        }
    }


    i = 1;
    while (i <= n1)
    {
        xi = xx[i];
        ni = i + 1;
        for (j = i+1; j<=n1; j++)
        {
            xj = xx[j];
            dxij = fabs( xi - xj );

            /* We can now the following simplification - since the stragglers
             * are also time ordered - if (xj-xi) exceeds condenseIn->bin_time
             * then it is meaningless to continue this loop - so ABORT. Note
             * that the test has to be done on (xj-xi) (order is important
             * here).
             */
            if ( (xj-xi) > condenseIn->bin_time )
            {
                break;
            }

            if ( dxij <= condenseIn->bin_time )
            {
                ni = ni + 1;
                if ( vv[j] > vv[i] )
                {
                    i = j;
                }
            }
        }

        /* We are now ready to append the i-th element in this list. The
         * master index has to be de-referenced to k first in order to do
         * this correctly.
         */
        {
            INT4 k = mid[i];

            /*-- Allocate memory for output --*/
            if ( !(*condenseOut = (trigScanClusterOut*)
                        LALRealloc(*condenseOut,
                            sizeof(trigScanClusterOut)*(*nclusters+1)))
               )
            {
                fprintf (stderr, "LALRealloc error in condenseout. \n");
                abort ();
            }

            /* Copy the elements to the newly created memory */
            (*condenseOut)[(*nclusters)].y          = masterList[k].y;
            (*condenseOut)[(*nclusters)].z          = masterList[k].z;
            (*condenseOut)[(*nclusters)].tc_sec     = masterList[k].tc_sec;
            (*condenseOut)[(*nclusters)].tc_ns      = masterList[k].tc_ns;
            (*condenseOut)[(*nclusters)].rho        = masterList[k].rho;
            (*condenseOut)[(*nclusters)].master_idx = k;
            (*condenseOut)[(*nclusters)].nelements  = 1;

            /* increment nclusters as we have added a new one */
            (*nclusters) ++;

#if 0
            /* After adding this element, print it to stderr and stdout */
            if (condenseIn->vrbflag)
            {
                fprintf (stderr, "Added cluster %3d after %3d (%3d members) max snr index "
                        "%3d %9d %9d %e\n",
                        (*nclusters), (*nclusters)-1, (*condenseOut)[(*nclusters)-1].nelements,
                        (*condenseOut)[(*nclusters)-1].master_idx,
                        (INT4)((*condenseOut)[(*nclusters)-1].tc_sec),
                        (INT4)((*condenseOut)[(*nclusters)-1].tc_ns),
                        (*condenseOut)[(*nclusters)-1].rho);
            }
#endif
        }

        i = ni;
    }

#if 0
    if ( condenseIn->vrbflag )
          fprintf (stderr, "--------- DONE APPENDING ----------\n");
#endif

    /* Free up memory */
    if (xx)  LALFree (xx);
    if (vv)  LALFree (vv);
    if (mid) LALFree (mid);

    /* Normal exit */
    DETATCHSTATUSPTR(status);
    RETURN(status);
}

/* ---------------------------------------------------------------------
 * This general purpose function can be used to delete all the elements
 * of a SnglInspiralTable.
 ----------------------------------------------------------------------*/
INT4 XLALDeleteSnglInspiralTable (
        SnglInspiralTable **eventHead
        )
{
    SnglInspiralTable    *thisEvent=NULL;

    /* Delete the original masterList */
    while ( (*eventHead) )
    {
        thisEvent = (*eventHead);
        (*eventHead) = (*eventHead)->next;
        XLALFreeSnglInspiral ( &thisEvent );
    }

    return (0);
}

/* ---------------------------------------------------------------------
 * This function populates some (not all) members of the trigScanClusterIn
 * structure using the information available in other structures. This is
 * typically called from inspiral.c
 ----------------------------------------------------------------------*/
INT4 XLALPopulateTrigScanInput (
        trigScanClusterIn     **condenseIn,
        FindChirpDataParams   *fcDataParams,
        FindChirpTmpltParams  *fcTmpltParams,
        FindChirpFilterParams *fcFilterParams,
        InspiralTemplate      *bankHead
        )
{
    static const char *func = "XLALPopulateTrigScanInput";
    UINT4        numPoints, ki;
    REAL8        sampleRate, deltaF, deltaT;

    if ( !fcFilterParams || !fcTmpltParams || !fcDataParams || !bankHead )
        XLAL_ERROR( func, XLAL_ENOMEM );

    (*condenseIn) = (trigScanClusterIn *)
            LALCalloc (1, sizeof(trigScanClusterIn));

    if ( !condenseIn )
        XLAL_ERROR( func, XLAL_ENOMEM );

    sampleRate = 1.0L/fcTmpltParams->deltaT;
    numPoints  = 2*(fcDataParams->wtildeVec->length - 1);
    deltaT     = fcTmpltParams->deltaT;
    deltaF     = sampleRate / (REAL8)(numPoints);

    (*condenseIn)->rho_th1     = sqrt(fcFilterParams->rhosqThresh);
    (*condenseIn)->chisq_th1   = fcFilterParams->chisqThresh;

    return 0;
}

/* ---------------------------------------------------------------------
 * This function takes a SnglInspiralTable and returns another SnglInspiral
 * table consisting of ncluster elements contained in clusterOut[i]->masterIdx
 * where i runs from 0 to ncluster-1.
 ----------------------------------------------------------------------*/
SnglInspiralTable *
XLALTrimSnglInspiralTable (
        SnglInspiralTable   **eventHead,
        trigScanClusterOut  *clusterOut,
        INT4                nclusters
        )
{
    static const char    *func = "XLALTrimSnglInspiralTable";
    static LALStatus     status;
    SnglInspiralTable    *output     = NULL;
    SnglInspiralTable    *tempList   = NULL;
    SnglInspiralTable    *thisEvent  = NULL;
    REAL4Vector          *clusterMasterIndex = NULL;
    INT4Vector           *heapSortIndex = NULL;
    INT4                 l, j, nOrgTrigs;

    /* if there are no events, then no-op */
    if ( ! *eventHead )
          return (0);

    /* We will assume that the clusters returned are NOT sorted in time */
    clusterMasterIndex = XLALCreateREAL4Vector( nclusters );
    heapSortIndex      = XLALCreateINT4Vector( nclusters );
    for (j=0; j<nclusters; j++)
          clusterMasterIndex->data[j] = (REAL4)(clusterOut[j].master_idx);

    LALSHeapIndex(&status,heapSortIndex,clusterMasterIndex);
    if(status.statusCode){
        LALPrintError("%s: %s\n", LALTRIGSCANCLUSTERC, "Error in Heap Sort");
        REPORTSTATUS(&status);
    }


    nOrgTrigs = XLALCountSnglInspiral(*eventHead);
    /* Copy the clustered events to tempList */
    for ( l=0, j=0; (l<nOrgTrigs) && (j<nclusters); l++ )
    {
        /* Point thisEvent to the correct event */
        if (l<1) /* if first event */
        {
            thisEvent = (*eventHead);
        }
        else /* subsequent event */
        {
            thisEvent = thisEvent->next;
        }

        if (l == clusterOut[heapSortIndex->data[j]].master_idx)
        {
            j = j+1;

            if ( !tempList )
            {
                output = tempList = (SnglInspiralTable *)
                        LALCalloc ( 1, sizeof(SnglInspiralTable) );
            }
            else
            {
                tempList = tempList->next = (SnglInspiralTable *)
                        LALCalloc ( 1, sizeof(SnglInspiralTable) );
            }

            /* If memory allocation failed we should free up everything */
            if ( !tempList )
            {
                XLALDeleteSnglInspiralTable (&output);
                XLAL_ERROR_NULL(func, XLAL_ENOMEM);
            }

            memcpy (tempList, thisEvent, sizeof(SnglInspiralTable));
            tempList->next = NULL;

            /* This is the place to append cluster specific information
             * to the sngl_inspiral table. At the moment we want to know
             * how big was the cluster. Since there is no element of the
             * sngl_inspiral table structure dedicated to store this
             * information, we will use the element alpha. However this is
             * not a good practice.
             */
            tempList->alpha = (REAL4)(clusterOut[heapSortIndex->data[j-1]].nelements);
        }
    }

    /* Delete the Vectors used for heap sort */
    XLALDestroyREAL4Vector( clusterMasterIndex );
    XLALDestroyINT4Vector( heapSortIndex );

    return (output);
}


