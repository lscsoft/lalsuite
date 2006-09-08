
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
 *     -1 for unclassified, 
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
    REAL8       y, z, tc_sec, tc_ns;  
    REAL8       rho;
    INT4        clusterID;
    REAL8       Gamma[6];
    gsl_matrix  *invGamma;
}
trigScanInputPoint;

typedef struct tagTrigScanClusterIn
{ 
    INT2                  vrbflag, appendStragglers;
    REAL8                 mmCoarse, rho_th1, chisq_th1;
    InputMasses           massChoice;
    REAL8                 bin_time, sf_volume;
    trigScanInputPoint    *masterList;
    trigScanType          scanMethod;
    INT4                  n;
    REAL8                 maxTcFootPrint;
}
trigScanClusterIn;

typedef struct tagTrigScanClusterOut
{
    REAL8  y, z, tc_sec, tc_ns;
    REAL8  rho;
    INT4   master_idx;
    INT4   cluster_id, nelements;
}
trigScanClusterOut;

typedef struct tagTrigScanEpsSearchIn
{
    trigScanInputPoint   *masterList;
    INT4                 nInputPoints;
    INT4                 clusterID;
    REAL8                maxTcFootPrint;
    INT4                 minLoopIdx;
}
trigScanEpsSearchInput;

/*--- Function prototypes ---*/
void LALTrigScanClusterDriver (
        LALStatus           *status,
        trigScanClusterIn   *clusterIn, 
        trigScanClusterOut  **clusterOut, 
        INT4                *nclusters
        );

/*--- Core functions which carry out clustering ---*/
trigScanValidEvent XLALTrigScanExpandCluster (
        INT4                  *list, 
        trigScanClusterIn     *condenseIn,
        INT4                  nPoints,
        INT4                  currClusterID,
        trigScanClusterOut    **condenseOut
        );

void XLALTrigScanGetEpsNeighbourhood (
        INT4                    seed, 
        INT4                    **list,
        INT4                    *size,
        trigScanEpsSearchInput  *epsSearchIn
        );

void LALTrigScanStoreThisCluster (
        LALStatus                *status,
        const INT4               *list,
        const trigScanClusterIn  *condenseIn,
        const INT4               size,
        const INT4               currClusterID,
        trigScanClusterOut       **condenseOut
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

