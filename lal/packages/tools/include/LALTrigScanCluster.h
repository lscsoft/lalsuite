
/*-------------------------------------------------------------------------
 *
 * File Name: LALTrigScanCluster.h
 *
 * Author: Sengupta, Anand. S. and Gupchup, Jayant A.
 *
 * Revision: $Id$
 *
 *----------------------------------------------------------------------- */

#if 0
<lalVerbatim file="LALTrigScanClusterHV">
Author: Sengupta, Anand. S. and Gupchup, Jayant A.
$Id$
</lalVerbatim>
#endif

#ifndef _LALTRIGSCANCLUSTER_H
#define _LALTRIGSCANCLUSTER_H

#include    <math.h>
#include    <lal/LALStdlib.h>
#include    <lal/LALInspiralBank.h>
#include    <lal/LALNoiseModels.h>
#include    <lal/LALInspiral.h>
#include    <lal/LALConstants.h>
#include    <lal/LALStdlib.h>
#include    <lal/Matrix.h>
#include    <lal/LIGOMetadataTables.h>
#include    <lal/LIGOMetadataUtils.h>
#include    <lal/AVFactories.h>
#include    <lal/FindChirp.h>
#include    <lal/Sort.h>
#include    <lal/EllipsoidOverlapTools.h>

NRCSID( LALTRIGSCANCLUSTERH, 
        "$Id$");

/* Cluster classification: 
 *      -1 for unclassified, 
 *      0 for noise and 
 *      an integer > 0 for a valid cluster 
 */

#define TRIGSCAN_UNCLASSIFIED     (-1)
#define TRIGSCAN_NOISE            (0)
#define TRIGSCAN_CLUSTER_MIN_PTS  (2)

/* Error messages */
#define LALTRIGSCANCLUSTERH_ENULL          1
#define LALTRIGSCANCLUSTERH_ECHOICE        2
#define LALTRIGSCANCLUSTERH_ESIZE          4

#define LALTRIGSCANCLUSTERH_MSGENULL      "Uexpected NULL pointer"
#define LALTRIGSCANCLUSTERH_MSGECHOICE    "Invalid input parameter"
#define LALTRIGSCANCLUSTERH_MSGESIZE      "Invalid input size"

typedef enum {
    trigScanNone,
    T0T3Tc,
    Psi0Psi3Tc
}
trigScanType;

typedef enum {
    trigScanFalse,
    trigScanTrue
}
trigScanValidEvent;

typedef struct tagTrigScanInputPoint
{ 
    REAL8       y, z, rho, effD, chisq;
    REAL8       tc_sec, tc_ns;  
    REAL8       figOfMerit;
    INT4        isValidEvent;
    INT4        clusterID;
    REAL8       Gamma[6];
    gsl_matrix  *invGamma;
}
trigScanInputPoint;

typedef struct tagTrigScanClusterIn
{ 
    INT2                  vrbflag, appendStragglers;
    REAL8                 fLower, fUpper, tSampling;
    REAL8                 mmCoarse, rho_th1, chisq_th1;
    REAL8FrequencySeries  coarseShf;
    Order                 order;
    Approximant           approximant;
    InputMasses           massChoice;
    REAL8                 bin_time, sf_volume;
    trigScanInputPoint    *masterList;
    trigScanType          scanMethod;
    INT4                  n;
}
trigScanClusterIn;

typedef struct tagTrigScanClusterOut
{
    REAL8  y, z, tc_sec, tc_ns;
    REAL8  rho, effD;
    INT4   master_idx;
    INT4   cluster_id, nelements;
}
trigScanClusterOut;

typedef struct tagTrigScanEpsSearchIn
{
    trigScanInputPoint   *masterList;
    INT4                 nInputPoints;
    INT4                 minPoints;
    INT4                 clusterID;
    INT4                 seedID;
}
trigScanEpsSearchInput;

typedef struct tagExpandClusterIn
{
    INT4   nInputPoints;  
    INT4   currClusterID;
    INT4   currSeedID;
}
ExpandClusterInput;

/*--- Function prototypes ---*/
void LALTrigScanClusterDriver (
        LALStatus           *status,
        trigScanClusterIn   *clusterIn, 
        trigScanClusterOut  **clusterOut, 
        INT4                *nclusters
        );

/*--- Core functions which carry out clustering ---*/
trigScanValidEvent XLALTrigScanExpandCluster (
        trigScanInputPoint    *list, 
        trigScanInputPoint    *masterList,
        ExpandClusterInput    expandClusterIn
        );

void XLALTrigScanGetEpsNeighbourhood (
        trigScanInputPoint      seed, 
        trigScanInputPoint      **list,
        INT4                    *size,
        trigScanEpsSearchInput  *epsSearchIn
        );

void LALTrigScanClusterMakeOutput (
        LALStatus               *status,
        trigScanClusterIn       *condenseIn, 
        trigScanClusterOut      **condenseOut,
        INT4                    nclusters
        );

void LALTrigScanAppendIsolatedTriggers (
        LALStatus               *status,
        trigScanClusterIn       *condenseIn, 
        trigScanClusterOut      **condenseOut,
        INT4                    *nclusters
        );

INT4 XLALDeleteSnglInspiralTable (
        SnglInspiralTable **eventHead
        );

SnglInspiralTable *  
XLALTrimSnglInspiralTable (
        SnglInspiralTable   **inspiralEventList,
        trigScanClusterOut  *clusterOut,
        INT4                nclusters
        );

INT4 XLALPopulateTrigScanInput (
        trigScanClusterIn     **condenseIn,
        FindChirpDataParams   *fcDataParams,
        FindChirpTmpltParams  *fcTmpltParams,
        FindChirpFilterParams *fcFilterParams,
        InspiralTemplate      *bankHead
        );

void LALClusterSnglInspiralOverTemplatesAndEndTime ( 
        LALStatus              *status,
        SnglInspiralTable      **eventHead,
        trigScanClusterIn      *condenseIn
        );

#endif /* _LALTRIGSCANCLUSTER_H */

