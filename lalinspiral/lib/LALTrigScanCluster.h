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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/


/*-------------------------------------------------------------------------
 *
 * File Name: LALTrigScanCluster.h
 *
 * Author: Sengupta, Anand. S., Gupchup, Jayant A. and Robinson, C. A. K.
 *
 *----------------------------------------------------------------------- */

#ifndef _LALTRIGSCANCLUSTER_H
#define _LALTRIGSCANCLUSTER_H

#include    <lal/LALStdlib.h>
#include    <lal/LIGOMetadataTables.h>
#include    <lal/LIGOMetadataUtils.h>
#include    <lal/Date.h>

#include    <lal/EllipsoidOverlapTools.h>
#include    <lal/TrigScanEThincaCommon.h>

/**
 * \defgroup LALTrigScanCluster_h Header LALTrigScanCluster.h
 * \ingroup lalinspiral_UNCLASSIFIED
 * \author Sengupta, Anand. S., Gupchup, Jayant A. and Robinson, C. A. K.
 * \brief NONE
 */
/** @{ */

/** UNDOCUMENTED */
typedef enum
tagtrigScanType
{
    trigScanNone,	/**< UNDOCUMENTED */
    T0T3Tc,		/**< UNDOCUMENTED */
    Psi0Psi3Tc,		/**< UNDOCUMENTED */
    NUM_TRIGSCAN_TYPE	/**< UNDOCUMENTED */
}
trigScanType;

/** UNDOCUMENTED */
typedef struct
tagTrigScanCluster
{
  INT4                      nelements;	/**< UNDOCUMENTED */
  TriggerErrorList          *element;	/**< UNDOCUMENTED */
  struct tagTrigScanCluster *next;	/**< UNDOCUMENTED */
} TrigScanCluster;

/** UNDOCUMENTED */
typedef enum
tagTrigScanStatus
{
  TRIGSCAN_SUCCESS,	/**< UNDOCUMENTED */
  TRIGSCAN_ERROR,	/**< UNDOCUMENTED */
  TRIGSCAN_NUM_STATUS	/**< UNDOCUMENTED */
} TrigScanStatus;	/**< UNDOCUMENTED */


int XLALTrigScanClusterTriggers( SnglInspiralTable **table,
                                 trigScanType      method,
                                 REAL8             scaleFactor,
                                 INT4              appendStragglers );

TrigScanCluster * XLALTrigScanCreateCluster( TriggerErrorList **errorListHead,
                                             REAL8            tcMax );

int XLALTrigScanRemoveStragglers( TrigScanCluster **clusters );

int XLALTrigScanKeepLoudestTrigger( TrigScanCluster *cluster );

int XLALTrigScanReLinkLists( TrigScanCluster *clusterHead );

void XLALTrigScanDestroyCluster( TrigScanCluster *cluster,
                                TrigScanStatus   status
                              );

/** @} */ /* end:LALTrigScanCluster_h */

#endif /* _LALTRIGSCANCLUSTER_H */
