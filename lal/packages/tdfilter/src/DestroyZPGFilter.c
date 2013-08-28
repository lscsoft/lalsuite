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
#include <lal/ZPGFilter.h>


/**
 * \addtogroup DestroyZPGFilter_c
 * \author Creighton, T. D.
 *
 * \brief Destroys ZPG filter objects.
 *
 * ### Description ###
 *
 * These functions destroy an object <tt>**output</tt> of type
 * \c COMPLEX8ZPGFilter or \c COMPLEX16ZPGFilter, and set
 * <tt>*output</tt> to \c NULL.
 *
 */
/*@{*/

/** \see See \ref DestroyZPGFilter_c for documentation */
void XLALDestroyCOMPLEX8ZPGFilter( COMPLEX8ZPGFilter *filter )
{
  if ( filter )
  {
    XLALDestroyCOMPLEX8Vector( filter->zeros );
    XLALDestroyCOMPLEX8Vector( filter->poles );
    LALFree( filter );
  }
  return;
}

/** \see See \ref DestroyZPGFilter_c for documentation */
void XLALDestroyCOMPLEX16ZPGFilter( COMPLEX16ZPGFilter *filter )
{
  if ( filter )
  {
    XLALDestroyCOMPLEX16Vector( filter->zeros );
    XLALDestroyCOMPLEX16Vector( filter->poles );
    LALFree( filter );
  }
  return;
}

/** \see See \ref DestroyZPGFilter_c for documentation */
void
LALDestroyCOMPLEX8ZPGFilter( LALStatus         *stat,
			     COMPLEX8ZPGFilter **input )
{
  INITSTATUS(stat);

  /* Make sure handle is non-null, and points to a non-null
     pointer. */
  ASSERT(input,stat,ZPGFILTERH_ENUL,ZPGFILTERH_MSGENUL);
  ASSERT(*input,stat,ZPGFILTERH_ENUL,ZPGFILTERH_MSGENUL);

  XLALDestroyCOMPLEX8ZPGFilter( *input );
  *input = NULL;

  /* Normal exit */
  RETURN(stat);
}

/** \see See \ref DestroyZPGFilter_c for documentation */
void
LALDestroyCOMPLEX16ZPGFilter( LALStatus          *stat,
			      COMPLEX16ZPGFilter **input )
{
  INITSTATUS(stat);

  /* Make sure handle is non-null, and points to a non-null
     pointer. */
  ASSERT(input,stat,ZPGFILTERH_ENUL,ZPGFILTERH_MSGENUL);
  ASSERT(*input,stat,ZPGFILTERH_ENUL,ZPGFILTERH_MSGENUL);

  XLALDestroyCOMPLEX16ZPGFilter( *input );
  *input = NULL;

  /* Normal exit */
  RETURN(stat);
}
/*@}*/
