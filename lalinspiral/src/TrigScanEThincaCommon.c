/*
*  Copyright (C) 2007  Craig Robinson
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
 * File Name: TrigScanEThincaCommon.c
 *
 * Author: Robinson, C. A. K.
 *
 *-----------------------------------------------------------------------
 */

#include <lal/TrigScanEThincaCommon.h>

/**
 * \author Robinson, C. A. K.
 * \file
 *
 * \brief Provides helper functions used in TrigScan and E-thinca.
 *
 * \heading{Description}
 * The function <tt>XLALCreateTriggerErrorList()</tt> creates a linked list of
 * structures pointing to the trigger and their associated position vector and
 * shape matrix. If required, the maximum difference in tC associated with the
 * triggers in the list of \c ::SnglInspiralTable's will be passed back in
 * \c tcMax.
 *
 * The function <tt>XLALDestroyTriggerErrorList()</tt> frees all memory associated
 * with the \c ::TriggerErrorList, with the exception of the wrapped
 * ::SnglInspiralTable's, which will normally still be required after
 * TrigScan and E-thinca have completed.
 *
 */

TriggerErrorList * XLALCreateTriggerErrorList( SnglInspiralTable *tableHead,
                                               REAL8             scaleFactor,
                                               REAL8             *tcMax )

{
  REAL8 timeError = 0.0;

  TriggerErrorList *errorListHead = NULL;
  TriggerErrorList *thisErrorList = NULL;

  SnglInspiralTable *currentTrigger = NULL;

#ifndef LAL_NDEBUG
  if ( !tableHead )
    XLAL_ERROR_NULL( XLAL_EFAULT );

  if ( scaleFactor <= 0 )
    XLAL_ERROR_NULL( XLAL_EINVAL );
#endif

  /* Loop through triggers and assign each of them an error ellipsoid */
  for (currentTrigger = tableHead; currentTrigger;
      currentTrigger = currentTrigger->next)
  {
    REAL8 thisTimeError;

    if (!errorListHead)
    {
      errorListHead = (TriggerErrorList *) LALCalloc(1, sizeof(TriggerErrorList));
      thisErrorList = errorListHead;
    }
    else
    {
      thisErrorList->next = (TriggerErrorList *) LALCalloc(1, sizeof(TriggerErrorList));
      thisErrorList = thisErrorList->next;
    }
    if ( !thisErrorList )
    {
      XLALDestroyTriggerErrorList( errorListHead );
      XLAL_ERROR_NULL( XLAL_ENOMEM );
    }

    thisErrorList->trigger    = currentTrigger;
    thisErrorList->err_matrix = XLALGetErrorMatrixFromSnglInspiral( currentTrigger,
                                  scaleFactor );
    if ( !thisErrorList->err_matrix )
    {
      XLALDestroyTriggerErrorList( errorListHead );
      XLAL_ERROR_NULL( XLAL_EFUNC );
    }

    thisErrorList->position   = XLALGetPositionFromSnglInspiral( currentTrigger );
    if ( !thisErrorList->position )
    {
      XLALDestroyTriggerErrorList( errorListHead );
      XLAL_ERROR_NULL( XLAL_EFUNC );
    }
    thisTimeError = XLALSnglInspiralTimeError(currentTrigger, scaleFactor );
    if (thisTimeError > timeError)
    {
      timeError = thisTimeError;
    }
  }

  /* Set the max timing error if required and return the list */
  if ( tcMax ) *tcMax = timeError;
  return errorListHead;
}


void XLALDestroyTriggerErrorList( TriggerErrorList *errorListHead )

{

  TriggerErrorList *thisErrorList;

  thisErrorList = errorListHead;
  while (thisErrorList)
  {
    errorListHead = thisErrorList->next;
    if ( thisErrorList->err_matrix )
      gsl_matrix_free( thisErrorList->err_matrix );

    if ( thisErrorList->position )
      gsl_vector_free( thisErrorList->position );

    LALFree( thisErrorList );
    thisErrorList = errorListHead;
  }
}


/**
 * Using the waveform metric components, translate an "e-thinca" treshold
 * into a \f$\Delta t\f$ error interval.
 */


REAL8 XLALSnglInspiralTimeError(const SnglInspiralTable *table, REAL8 eMatch)
{
  REAL8 a11 = table->Gamma[0] / eMatch;
  REAL8 a12 = table->Gamma[1] / eMatch;
  REAL8 a13 = table->Gamma[2] / eMatch;
  REAL8 a22 = table->Gamma[3] / eMatch;
  REAL8 a23 = table->Gamma[4] / eMatch;
  REAL8 a33 = table->Gamma[5] / eMatch;
  REAL8 x;
  REAL8 denom;

  x = (a23 * a23 - a22 * a33) * a22;
  denom = (a12*a23 - a22*a13) * (a12*a23 - a22*a13)
              - (a23*a23 - a22*a33) * (a12*a12 - a22*a11);

  if (denom == 0)
    XLAL_ERROR_REAL8(XLAL_EFPDIV0);
  if ((x < 0) ^ (denom < 0))
    XLAL_ERROR_REAL8(XLAL_EFPINVAL);

  return sqrt( x / denom );
}
