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


/*-------------------------------------------------------------------------
 *
 * File Name: LALTrigScanCluster.h
 *
 * Author: Sengupta, Anand. S., Gupchup, Jayant A. and Robinson, C. A. K.
 *
 * Revision: $Id: LALTrigScanCluster.h,v 1.10 2007/06/08 14:41:57 bema Exp $
 *
 *----------------------------------------------------------------------- */

#if 0
<lalVerbatim file="LALTrigScanClusterHV">
Author: Sengupta, Anand. S., Gupchup, Jayant A. and Robinson, C. A. K.
$Id: LALTrigScanCluster.h,v 1.10 2007/06/08 14:41:57 bema Exp $
</lalVerbatim>
#endif

#ifndef _LALTRIGSCANCLUSTER_H
#define _LALTRIGSCANCLUSTER_H

#include    <lal/LALStdlib.h>
#include    <lal/LIGOMetadataTables.h>
#include    <lal/LIGOMetadataUtils.h>
#include    <lal/Date.h>

#include    <lal/EllipsoidOverlapTools.h>
#include    <lal/CoincInspiralEllipsoid.h>

NRCSID( LALTRIGSCANCLUSTERH,
        "$Id$");


typedef enum {
    trigScanNone,
    T0T3Tc,
    Psi0Psi3Tc,
    NUM_TRIGSCAN_TYPE
}
trigScanType;

typedef struct
tagTrigScanCluster
{
  INT4                      nelements;
  TriggerErrorList          *element;
  struct tagTrigScanCluster *next;
} TrigScanCluster;

typedef enum
{
  TRIGSCAN_SUCCESS,
  TRIGSCAN_ERROR,
  TRIGSCAN_NUM_STATUS
} TrigScanStatus;


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
#endif /* _LALTRIGSCANCLUSTER_H */

