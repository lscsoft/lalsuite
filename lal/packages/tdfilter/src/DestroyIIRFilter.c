/*
*  Copyright (C) 2007 Jolien Creighton, Teviet Creighton
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

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/IIRFilter.h>


/**
 * \addtogroup DestroyIIRFilter_c
 * \author Creighton, T. D.
 *
 * \brief Destroys IIR filter objects.
 *
 * ### Description ###
 *
 * These functions destroy an object <tt>**input</tt> of type
 * \c REAL4IIRFilter or \c REAL8IIRFilter, and set
 * <tt>*input</tt> to \c NULL.
 *
 */
/*@{*/

/** \see See \ref DestroyIIRFilter_c for documentation */
void XLALDestroyREAL4IIRFilter( REAL4IIRFilter *filter )
{
  if ( filter )
  {
    XLALDestroyREAL4Vector( filter->directCoef );
    XLALDestroyREAL4Vector( filter->recursCoef );
    XLALDestroyREAL4Vector( filter->history );
    LALFree( filter );
  }
  return;
}

/** \see See \ref DestroyIIRFilter_c for documentation */
void XLALDestroyREAL8IIRFilter( REAL8IIRFilter *filter )
{
  if ( filter )
  {
    XLALDestroyREAL8Vector( filter->directCoef );
    XLALDestroyREAL8Vector( filter->recursCoef );
    XLALDestroyREAL8Vector( filter->history );
    LALFree( filter );
  }
  return;
}

/**
 * Deprecated.
 * \deprecated Use XLALDestroyREAL4IIRFilter() instead
 */
void
LALDestroyREAL4IIRFilter( LALStatus      *stat,
			  REAL4IIRFilter **input )
{
  INITSTATUS(stat);

  /* Make sure handle is non-null, and points to a non-null pointer.
     (The routine LALSDestroyVector will check that the data fields are
     non-null.) */
  ASSERT(input,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(*input,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);

  /* Free the filter, then point the handle to NULL. */
  XLALDestroyREAL4IIRFilter( *input );
  *input=NULL;

  /* Normal exit */
  RETURN(stat);
}


/**
 * Deprecated.
 * \deprecated Use XLALDestroyREAL8IIRFilter() instead
 */
void
LALDestroyREAL8IIRFilter( LALStatus      *stat,
			  REAL8IIRFilter **input )
{
  INITSTATUS(stat);

  /* Make sure handle is non-null, and points to a non-null pointer.
     (The routine LALDDestroyVector will check that the data fields are
     non-null.) */
  ASSERT(input,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(*input,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);

  /* Free the filter, then point the handle to NULL. */
  XLALDestroyREAL8IIRFilter( *input );

  /* Normal exit */
  RETURN(stat);
}

/*@}*/
