/*
*  Copyright (C) 2007 Craig Robinson
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

/*-----------------------------------------------------------------------
 *
 * File Name: TrigScanEThincaCommon.h
 *
 * Author: Robinson, C. A. K.
 *
 *-----------------------------------------------------------------------
 */

#ifndef _TRIGSCANETHINCACOMMON_H
#define _TRIGSCANETHINCACOMMON_H

#include <lal/LALAtomicDatatypes.h>
#include <lal/LIGOMetadataTables.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \defgroup TrigScanEThincaCommon_h Header TrigScanEThincaCommon.h
 * \ingroup lalinspiral_UNCLASSIFIED
 * \author Robinson, C. A. K.
 *
 * \brief Provides helper functions common to TrigScan and E-thinca.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/TrigScanEThincaCommon.h>
 * \endcode
 *
 * This header provides functions used for creating and destroying the
 * linked lists used in TrigScan and E-thinca.
 *
 */
/** @{ */

/**
 * The \c TriggerErrorList is a linked list used within e-thinca. It
 * contains pointers to the \c SnglInspiralTable for a given trigger,
 * and its associated error matrix and position vector.
 */
typedef struct tagTriggerErrorList
{
  SnglInspiralTable          *trigger;
  gsl_matrix                 *err_matrix;
  gsl_vector                 *position;
  struct tagTriggerErrorList *next;
}
TriggerErrorList;

TriggerErrorList * XLALCreateTriggerErrorList( SnglInspiralTable *tableHead,
                                               REAL8             scaleFactor,
                                               REAL8             *tcMax );

void XLALDestroyTriggerErrorList( TriggerErrorList *errorListHead );

REAL8 XLALSnglInspiralTimeError(const SnglInspiralTable *table, REAL8 eMatch);

/** @} */ /* end:TrigScanEThincaCommon.h */

#ifdef  __cplusplus
}
#endif

#endif /* _TRIGSCANETHINCACOMMON_H */
