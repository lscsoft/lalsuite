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
        trigScanInputPoint    *masterList,
        INT4                  nPoints,
        INT4                  currClusterID
        )
{
    trigScanValidEvent      flag = trigScanFalse;
    INT4                    seed;
    trigScanEpsSearchInput  epsSearchIn; /* Data structure given as input to*/
                                         /* getEpsNeighbourhood fn.         */ 

    INT4   pointer=0, size=1; /* when pointer < size, all elements */
                              /* have been accessed                */

    /* Create the epsSearchIn data structure */
    epsSearchIn.masterList   = masterList;
    epsSearchIn.nInputPoints = nPoints;
    epsSearchIn.clusterID    = currClusterID;

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
    if (size > 1) flag = trigScanTrue;

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

    for (i = 0; i < epsSearchIn->nInputPoints; i++)
    {
        /* check if the point has been classified */
        if (epsSearchIn->masterList[i].clusterID == TRIGSCAN_UNCLASSIFIED)
        {
            /* Set the position vector (q2) of the i-th point */
            q2[0] = epsSearchIn->masterList[i].tc_sec 
                    + 1.e-9*(epsSearchIn->masterList[i].tc_ns);
            q2[1] = epsSearchIn->masterList[i].y;
            q2[2] = epsSearchIn->masterList[i].z;

            /* create a vector view from the q2 array */
            vq2 = gsl_vector_view_array(q2,3);

            /* Set the shape matrix of the i-th point */
            workSpace->invQ2    = epsSearchIn->masterList[i].invGamma;

            /* Figure out if the above ellipsoids overlap */
            distance = XLALCheckOverlapOfEllipsoids (&(vq1.vector), &(vq2.vector), workSpace);

            if (distance > 0 && distance <= 1.) 
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
            }
        }
    }

    /* De-allocate workSpace allocated earlier */
    /* This workSpace is used to check overlap */
    /* of ambiguity ellipsoids                 */
    XLALFreeFContactWorkSpace( workSpace );
}

/* --------------------------------------------------------------------- 
 * This function is used to fillout the trigScanClusterOut structure after
 * the clustering has been done. 
 ---------------------------------------------------------------------*/ 
void LALTrigScanClusterMakeOutput (
        LALStatus               *status,
        trigScanClusterIn       *condenseIn, 
        trigScanClusterOut      **condenseOut,
        INT4                    nclusters
        )
{ 
    INT4          i, j, n, cSize;
    INT4          *t_idx = NULL, maxResultIdx;
    REAL8         maxResult;
    REAL8Vector   *t_rho = NULL ;

    trigScanInputPoint     *masterList=NULL;

    INITSTATUS (status, "LALTrigScanClusterMakeOutput", LALTRIGSCANCLUSTERC);
    ATTATCHSTATUSPTR(status);

    n            = condenseIn->n;
    masterList   = condenseIn->masterList;

    ASSERT (condenseOut, 
            status, LALTRIGSCANCLUSTERH_ENULL, LALTRIGSCANCLUSTERH_MSGENULL);
    ASSERT (masterList, 
            status, LALTRIGSCANCLUSTERH_ENULL, LALTRIGSCANCLUSTERH_MSGENULL);
    ASSERT (n > 0, status, 
            LALTRIGSCANCLUSTERH_ECHOICE, LALTRIGSCANCLUSTERH_MSGECHOICE);
    ASSERT (nclusters > 0, status, 
            LALTRIGSCANCLUSTERH_ECHOICE, LALTRIGSCANCLUSTERH_MSGECHOICE);

    t_idx  = (INT4 *) LALMalloc (n * sizeof(INT4));
    LALDCreateVector (status->statusPtr, &t_rho, n);
    CHECKSTATUSPTR (status);

    /*-- Loop over the clusters --*/
    for (i = 1; i <= nclusters; i++) {

        cSize = 0;

        /* loop over all the elements of masterList */
        for (j = 0; j < n; j++) {
            if (masterList[j].clusterID == i) {
                t_rho->data [cSize] = masterList[j].rho;
                t_idx  [cSize] = j;

#if 0
                if (condenseIn->vrbflag) 
                {
                    fprintf (stdout, "%d %e %e %9.f %9.f %.8e\n", 
                            i, masterList[j].y, masterList[j].z, 
                            masterList[j].tc_sec, masterList[j].tc_ns,
                            masterList[j].rho);
                }
#endif

                cSize ++;
            }
        }
        t_rho->length = cSize;

        /* Find the loudest event in this cluster */
        LALDMax ( status->statusPtr, &maxResult, t_rho, &maxResultIdx);
        CHECKSTATUSPTR (status);

        (*condenseOut)[i-1].y          = masterList[t_idx[maxResultIdx]].y;
        (*condenseOut)[i-1].z          = masterList[t_idx[maxResultIdx]].z;
        (*condenseOut)[i-1].tc_sec     = masterList[t_idx[maxResultIdx]].tc_sec;
        (*condenseOut)[i-1].tc_ns      = masterList[t_idx[maxResultIdx]].tc_ns;
        (*condenseOut)[i-1].rho        = masterList[t_idx[maxResultIdx]].rho;
        (*condenseOut)[i-1].master_idx = t_idx[maxResultIdx];
        (*condenseOut)[i-1].nelements  = cSize;

#if 0
        if (condenseIn->vrbflag)
        {
            fprintf (stderr, 
                    "Found cluster %3d of %3d (%3d members) max snr index %3d %9d %9d %e\n", 
                    i, nclusters, cSize, 
                    (*condenseOut)[i-1].master_idx,
                    (INT4)((*condenseOut)[i-1].tc_sec),
                    (INT4)((*condenseOut)[i-1].tc_ns),
                    (*condenseOut)[i-1].rho);
        }
#endif

    } /* Loop over clusters */

    /* Free memory */
    if (t_idx) LALFree (t_idx);
    if (t_rho) 
    {
        LALDDestroyVector (status->statusPtr, &t_rho);
        CHECKSTATUSPTR (status);
    }

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
    REAL8                *xx=NULL, *vv=NULL, *nn=NULL;
    INT4                 *mid = NULL;

    INITSTATUS (status, 
            "LALTrigScanAppendIsolatedTriggers", LALTRIGSCANCLUSTERC);
    ATTATCHSTATUSPTR(status);

    n            = condenseIn->n;
    masterList   = condenseIn->masterList;

    mid  = LALCalloc (1, (n+1)*sizeof(INT4));
    xx   = LALCalloc (1, (n+1)*sizeof(REAL8));
    vv   = LALCalloc (1, (n+1)*sizeof(REAL8));
    nn   = LALCalloc (1, (n+1)*sizeof(REAL8));

    ASSERT (condenseOut, 
            status, LALTRIGSCANCLUSTERH_ENULL, LALTRIGSCANCLUSTERH_MSGENULL);
    ASSERT (masterList, 
            status, LALTRIGSCANCLUSTERH_ENULL, LALTRIGSCANCLUSTERH_MSGENULL);
    ASSERT (n > 0, status, 
            LALTRIGSCANCLUSTERH_ECHOICE, LALTRIGSCANCLUSTERH_MSGECHOICE);
    ASSERT (*nclusters >= 0, status, 
            LALTRIGSCANCLUSTERH_ECHOICE, LALTRIGSCANCLUSTERH_MSGECHOICE);
    ASSERT (xx && vv && nn && mid,
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
            nn[n1]  = (REAL8)(masterList[i-1].clusterID);
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
            dxij = condenseIn->tSampling * fabs( xi - xj );
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
         * this correctly
         */
        {
            INT4 k;

            k = mid[i];

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
                fprintf (stdout, "%4d %e %e %9.f %9.f %.8e\n", 
                        (*nclusters), masterList[k].y, masterList[k].z, 
                        masterList[k].tc_sec, masterList[k].tc_ns,
                        masterList[k].rho);
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
    if (nn)  LALFree (nn);
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

    (*condenseIn)->fLower      = fcDataParams->fLow;
    (*condenseIn)->fUpper      = sampleRate/2.0L - deltaF;
    (*condenseIn)->tSampling   = sampleRate;
    (*condenseIn)->mmCoarse    = bankHead->minMatch;
    (*condenseIn)->rho_th1     = sqrt(fcFilterParams->rhosqThresh);
    (*condenseIn)->chisq_th1   = fcFilterParams->chisqThresh;
    (*condenseIn)->order       = twoPN;
    (*condenseIn)->approximant = fcFilterParams->approximant;

    (*condenseIn)->coarseShf.f0     = 0.0L;
    (*condenseIn)->coarseShf.deltaF = deltaF;
    (*condenseIn)->coarseShf.data   = 
            XLALCreateREAL8Vector( (numPoints/2 + 1) );
    
    if ( !(*condenseIn)->coarseShf.data )
          XLAL_ERROR( func, XLAL_ENOMEM );

    /* Populate the shf vector from the wtilde vector */
    for (ki=0; ki<fcDataParams->wtildeVec->length ; ki++) 
    {
        if (fcDataParams->wtildeVec->data[ki].re) 
        {
            (*condenseIn)->coarseShf.data ->data[ki] = (REAL8) 
                    (1./fcDataParams->wtildeVec->data[ki].re);
        }
        else 
        {
            /* Note that we can safely set shf to be zero as this is    */
            /* correctly handled in the LALGetInspiralMoments function  */
            (*condenseIn)->coarseShf.data->data[ki] = 0.0;
        }
    }

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


