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
 * \addtogroup CreateZPGFilter_c
 * \author Creighton, T. D.
 *
 * \brief Creates ZPG filter objects.
 *
 * ### Description ###
 *
 * These functions create an object <tt>**output</tt>, of type
 * \c COMPLEX8ZPGFilter or \c COMPLEX16ZPGFilter, having
 * \c numZeros zeros and \c numPoles poles.  The values of those
 * zeros and poles are not set by these routines (in general they will
 * start out as garbage).  The handle passed into the functions must be a
 * valid handle (ie \c output\f$\neq\f$\c NULL), but must not
 * point to an existing object (ie <tt>*output</tt>=\c NULL).
 *
 */
/*@{*/

/** \see See \ref CreateZPGFilter_c for documentation */
COMPLEX8ZPGFilter *XLALCreateCOMPLEX8ZPGFilter( INT4 numZeros, INT4 numPoles )
{
  COMPLEX8ZPGFilter *output;
  if ( numZeros < 0 || numPoles < 0 )
    XLAL_ERROR_NULL( XLAL_EINVAL );
  output = LALCalloc( 1, sizeof(*output) );
  if ( ! output )
    XLAL_ERROR_NULL( XLAL_ENOMEM );

  /* Allocate the data fields.  If the number of poles or zeros is 0,
     the corresponding field(s) should remain null. */
  if ( numZeros > 0 )
    if ( ! ( output->zeros = XLALCreateCOMPLEX8Vector( numZeros ) ) )
    {
      XLALDestroyCOMPLEX8ZPGFilter( output );
      XLAL_ERROR_NULL( XLAL_EFUNC );
    }
  if ( numPoles > 0 )
    if ( ! ( output->poles = XLALCreateCOMPLEX8Vector( numPoles ) ) )
    {
      XLALDestroyCOMPLEX8ZPGFilter( output );
      XLAL_ERROR_NULL( XLAL_EFUNC );
    }

  return output;
}

/** \see See \ref CreateZPGFilter_c for documentation */
COMPLEX16ZPGFilter *XLALCreateCOMPLEX16ZPGFilter( INT4 numZeros, INT4 numPoles )
{
  COMPLEX16ZPGFilter *output;
  if ( numZeros < 0 || numPoles < 0 )
    XLAL_ERROR_NULL( XLAL_EINVAL );
  output = LALCalloc( 1, sizeof(*output) );
  if ( ! output )
    XLAL_ERROR_NULL( XLAL_ENOMEM );

  /* Allocate the data fields.  If the number of poles or zeros is 0,
     the corresponding field(s) should remain null. */
  if ( numZeros > 0 )
    if ( ! ( output->zeros = XLALCreateCOMPLEX16Vector( numZeros ) ) )
    {
      XLALDestroyCOMPLEX16ZPGFilter( output );
      XLAL_ERROR_NULL( XLAL_EFUNC );
    }
  if ( numPoles > 0 )
    if ( ! ( output->poles = XLALCreateCOMPLEX16Vector( numPoles ) ) )
    {
      XLALDestroyCOMPLEX16ZPGFilter( output );
      XLAL_ERROR_NULL( XLAL_EFUNC );
    }

  return output;
}

/** \see See \ref CreateZPGFilter_c for documentation */
void
LALCreateCOMPLEX8ZPGFilter( LALStatus         *stat,
			    COMPLEX8ZPGFilter **output,
			    INT4              numZeros,
			    INT4              numPoles )
{
  INITSTATUS(stat);

  /* Make sure that the output handle exists, but points to a null
     pointer. */
  ASSERT(output,stat,ZPGFILTERH_ENUL,ZPGFILTERH_MSGENUL);
  ASSERT(!*output,stat,ZPGFILTERH_EOUT,ZPGFILTERH_MSGEOUT);

  /* Make sure that numZeros and numPoles are non-negative. */
  ASSERT(numZeros>=0,stat,ZPGFILTERH_EBAD,ZPGFILTERH_MSGEBAD);
  ASSERT(numPoles>=0,stat,ZPGFILTERH_EBAD,ZPGFILTERH_MSGEBAD);

  /* Create the output structure. */
  *output = XLALCreateCOMPLEX8ZPGFilter( numZeros, numPoles );
  if ( ! *output )
  {
    ABORT(stat,ZPGFILTERH_EMEM,ZPGFILTERH_MSGEMEM);
  }

  /* Normal exit */
  RETURN(stat);
}


/** \see See \ref CreateZPGFilter_c for documentation */
void
LALCreateCOMPLEX16ZPGFilter( LALStatus          *stat,
			     COMPLEX16ZPGFilter **output,
			     INT4               numZeros,
			     INT4               numPoles )
{
  INITSTATUS(stat);

  /* Make sure that the output handle exists, but points to a null
     pointer. */
  ASSERT(output,stat,ZPGFILTERH_ENUL,ZPGFILTERH_MSGENUL);
  ASSERT(!*output,stat,ZPGFILTERH_EOUT,ZPGFILTERH_MSGEOUT);

  /* Make sure that numZeros and numPoles are non-negative. */
  ASSERT(numZeros>=0,stat,ZPGFILTERH_EBAD,ZPGFILTERH_MSGEBAD);
  ASSERT(numPoles>=0,stat,ZPGFILTERH_EBAD,ZPGFILTERH_MSGEBAD);

  /* Create the output structure. */
  *output = XLALCreateCOMPLEX16ZPGFilter( numZeros, numPoles );
  if ( ! *output )
  {
    ABORT(stat,ZPGFILTERH_EMEM,ZPGFILTERH_MSGEMEM);
  }

  /* Normal exit */
  RETURN(stat);
}
/*@}*/
