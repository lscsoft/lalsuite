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
$Id$"
</lalVerbatim>
#endif

#include <LALTrigScanCluster.h>

NRCSID (LALTRIGSCANCLUSTERC, 
        "$Id$");

/* ------------------------------------------------------------------*/
/* This function takes in a seed point and 'expands' that seed point */
/* to include other points from a masterList which are part of the   */
/* same cluster. If the input seed cannot be 'expand'ed then it      */
/* returns false. Otherwise it returns true.                         */
/* ------------------------------------------------------------------*/
trigScanValidEvent XLALTrigSCanExpandCluster (
        trigScanInputPoint    *list, 
        trigScanInputPoint    *masterList,
        ExpandClusterInput  expandClusterIn
        )
{
    trigScanValidEvent      flag = false;
    trigScanInputPoint      seed;
    trigScanEpsSearchInput  epsSearchIn; /* Data structure given as input to*/
                                         /* getEpsNeighbourhood fn.         */ 

    INT4   pointer=0, size=1; /* when pointer < size, all elements */
                              /* have been accessed                */

    /* Create the epsSearchIn data structure */
    epsSearchIn.masterList   = masterList;
    epsSearchIn.nInputPoints = expandClusterIn.nInputPoints;
    epsSearchIn.epsX         = expandClusterIn.epsX;
    epsSearchIn.epsY         = expandClusterIn.epsY;
    epsSearchIn.epsTc        = expandClusterIn.epsTc;
    epsSearchIn.alpha        = expandClusterIn.alpha;
    epsSearchIn.minPoints    = TRIGSCAN_CLUSTER_MIN_PTS;
    epsSearchIn.clusterID    = expandClusterIn.currClusterID;

    while (pointer < size) { 

        /* Pick the first point in this list as the seed */
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
    if (size > 1) flag = true;

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
        trigScanInputPoint      seed, 
        trigScanInputPoint      **list,
        INT4                  *size,
        trigScanEpsSearchInput  *epsSearchIn
        )
{
    INT4     i;
    REAL8    epsX  = epsSearchIn->epsX;
    REAL8    epsY  = epsSearchIn->epsY;
    REAL8    epsTc = epsSearchIn->epsTc;
    REAL8    alpha = epsSearchIn->alpha;
    REAL8    distance;
    REAL8    dx, dy, dtc;

    for (i = 0; i < epsSearchIn->nInputPoints; i++)
    {
        /* check if the point has been classified */
        if (epsSearchIn->masterList[i].clusterID == TRIGSCAN_UNCLASSIFIED)
        {
            dx   = epsSearchIn->masterList[i].x - seed.x;
            dy   = epsSearchIn->masterList[i].y - seed.y;
            dtc  = epsSearchIn->masterList[i].tc_sec - seed.tc_sec;
            dtc += 1.e-9*(epsSearchIn->masterList[i].tc_ns - seed.tc_ns);

            distance = XLALTrigScanGetDistance (dx, dy, dtc, epsX, epsY, 
                    alpha, epsTc);

            if (distance > 0 && distance <= 1.) 
            {
                /* set the clusterID to the currClusterID */
                epsSearchIn->masterList[i].clusterID = epsSearchIn->clusterID; 

                /* increment the size variable and hence realloc the list */
                (*size)++;

                if ( !(*list = (trigScanInputPoint*) 
                            LALRealloc(*list, 
                                sizeof(trigScanInputPoint)*(*size)))
                   )
                {
                    fprintf (stderr, "LALRealloc error. Aborting at %d\n", 
                            *size);
                    abort ();
                }

                /* add the shortlisted point to the list at position size - 1*/
                (*list)[*size - 1].x  = epsSearchIn->masterList[i].x;
                (*list)[*size - 1].y  = epsSearchIn->masterList[i].y;
                (*list)[*size - 1].tc_sec = epsSearchIn->masterList[i].tc_sec;
                (*list)[*size - 1].tc_ns  = epsSearchIn->masterList[i].tc_ns;	
                (*list)[*size - 1].rho    = epsSearchIn->masterList[i].rho;
                (*list)[*size - 1].isValidEvent =
                        epsSearchIn->masterList[i].isValidEvent;
                (*list)[*size - 1].clusterID = epsSearchIn->clusterID;
            }
        }
    }
}

/* --------------------------------------------------------------------- 
 * returns a value >=0 given the offsets dx, dy, dz and 
 * a, b, c and angle alpha for an ellipsoid. If the value is <=1.0 it
 * implies that the point is inside the ellipsoid. 
 ---------------------------------------------------------------------*/ 
REAL8 XLALTrigScanGetDistance (
        REAL8 dx,  REAL8 dy, 
        REAL8 dtc, REAL8 a, 
        REAL8 b,   REAL8 alpha, 
        REAL8 c)
{
    REAL8 xpp, ypp;

    xpp =  dx*cos(alpha) + dy*sin(alpha);
    ypp = -dx*sin(alpha) + dy*cos(alpha);

    return (pow(xpp/a, 2.) + pow (ypp/b, 2.) + pow (dtc/c, 2.));
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
                
                if (condenseIn->vrbflag) 
                {
                    fprintf (stdout, "%d %e %e %9.f %9.f %.8e %e %e %e %e\n", 
                            i, masterList[j].x, masterList[j].y, 
                            masterList[j].tc_sec, masterList[j].tc_ns,
                            masterList[j].rho, condenseIn->a[j], 
                            condenseIn->b[j], condenseIn->theta[j], masterList[j].effD);
                }

                cSize ++;
            }
        }
        t_rho->length = cSize;

        /* Find the loudest event in this cluster */
        LALDMax ( status->statusPtr, &maxResult, t_rho, &maxResultIdx);
        CHECKSTATUSPTR (status);

        (*condenseOut)[i-1].x          = masterList[t_idx[maxResultIdx]].x;
        (*condenseOut)[i-1].y          = masterList[t_idx[maxResultIdx]].y;
        (*condenseOut)[i-1].tc_sec     = masterList[t_idx[maxResultIdx]].tc_sec;
        (*condenseOut)[i-1].tc_ns      = masterList[t_idx[maxResultIdx]].tc_ns;
        (*condenseOut)[i-1].rho        = masterList[t_idx[maxResultIdx]].rho;
        (*condenseOut)[i-1].master_idx = t_idx[maxResultIdx];
        (*condenseOut)[i-1].nelements  = cSize;

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
    ASSERT (nclusters >= 0, status, 
            LALTRIGSCANCLUSTERH_ECHOICE, LALTRIGSCANCLUSTERH_MSGECHOICE);
    ASSERT (xx && vv && nn && mid,
            status, LALTRIGSCANCLUSTERH_ENULL, LALTRIGSCANCLUSTERH_MSGENULL);


    if ( condenseIn->vrbflag )
          fprintf (stderr, "--------- BEGINNING TO APPEND ----------\n");
    for (i=1, n1=0; i<=n; i++) {
        if (masterList[i-1].clusterID < 1 ) {
            n1++;
            xx[n1]  = masterList[i-1].tc_sec + masterList[i-1].tc_ns/1.0e9;
            vv[n1]  = masterList[i-1].rho;
            nn[n1]  = (REAL8)(masterList[i-1].clusterID);
            mid[n1] = i-1; 

            if ( condenseIn->vrbflag )
                  fprintf (stderr, "%4d   %d  %d   %d   %d   %1.12e\n", 
                          n1, mid[n1],
                          (masterList[i-1].clusterID), 
                          (int)(masterList[i-1].tc_sec),
                          (int)(masterList[i-1].tc_ns),
                          masterList[i-1].rho);
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
            (*condenseOut)[(*nclusters)].x          = masterList[k].x;
            (*condenseOut)[(*nclusters)].y          = masterList[k].y;
            (*condenseOut)[(*nclusters)].tc_sec     = masterList[k].tc_sec;
            (*condenseOut)[(*nclusters)].tc_ns      = masterList[k].tc_ns;
            (*condenseOut)[(*nclusters)].rho        = masterList[k].rho;
            (*condenseOut)[(*nclusters)].master_idx = k;
            (*condenseOut)[(*nclusters)].nelements  = 1;

            /* increment nclusters as we have added a new one */
            (*nclusters) ++;

            /* After adding this element, print it to stderr and stdout */
            if (condenseIn->vrbflag) 
            {
                fprintf (stdout, "%4d %e %e %9.f %9.f %.8e %e %e %e %e\n", 
                        (*nclusters), masterList[k].x, masterList[k].y, 
                        masterList[k].tc_sec, masterList[k].tc_ns,
                        masterList[k].rho, condenseIn->a[k], 
                        condenseIn->b[k], condenseIn->theta[k], masterList[k].effD);
                fprintf (stderr, "Added cluster %3d after %3d (%3d members) max snr index "
                        "%3d %9d %9d %e\n", 
                        (*nclusters), (*nclusters)-1, (*condenseOut)[(*nclusters)-1].nelements, 
                        (*condenseOut)[(*nclusters)-1].master_idx, 
                        (INT4)((*condenseOut)[(*nclusters)-1].tc_sec), 
                        (INT4)((*condenseOut)[(*nclusters)-1].tc_ns), 
                        (*condenseOut)[(*nclusters)-1].rho); 
            }
        }

        i = ni;
    }

    if ( condenseIn->vrbflag )
          fprintf (stderr, "--------- DONE APPENDING ----------\n");

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
    UINT4        numPoints, ki;
    REAL8        sampleRate, deltaF, deltaT;

    (*condenseIn) = (trigScanClusterIn *)
            LALCalloc (1, sizeof(trigScanClusterIn));

    sampleRate = 1.0L/fcTmpltParams->deltaT;
    numPoints  = fcTmpltParams->xfacVec->length;
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
        }
    }

    /* Delete the Vectors used for heap sort */
    XLALDestroyREAL4Vector( clusterMasterIndex );
    XLALDestroyINT4Vector( heapSortIndex );

    return (output);
}


