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
#include <lal/IIRFilter.h>

/**
 * \addtogroup IIRFilter_c
 * \author Creighton, T. D.
 *
 * \brief Computes an instant-by-instant IIR filter response.
 *
 * ### Description ###
 *
 * These functions pass a time-domain datum to an object <tt>*filter</tt>
 * of type \c REAL4IIRFilter or \c REAL8IIRFilter, and return the
 * filter response.  This is done using the auxiliary data series
 * formalism described in \ref IIRFilter_h.
 *
 * There are two pairs of routines in this module.  The functions
 * <tt>LALIIRFilterREAL4()</tt> and <tt>LALIIRFilterREAL8()</tt> conform to the LAL
 * standard, with status handling and error trapping; the input datum is
 * passed in as \c input and the response is returned in
 * <tt>*output</tt>.  The functions <tt>LALSIIRFilter()</tt> and
 * <tt>LALDIIRFilter()</tt> are non-standard lightweight routines, which may
 * be more suitable for multiple callings from the inner loops of
 * programs; they have no status handling or error trapping.  The input
 * datum is passed in by the variable \c x, and the response is
 * returned through the function's return statement.
 *
 */
/*@{*/

/**
 * WARNING: THIS FUNCTION IS OBSOLETE.
 * \deprecated
 */
void
LALIIRFilterREAL4( LALStatus      *stat,
		   REAL4          *output,
		   REAL4          input,
		   REAL4IIRFilter *filter )
{
  INT4 j;      /* Index for filter coefficients. */
  INT4 jmax;   /* Number of filter coefficients. */
  REAL4 *coef; /* Values of filter coefficients. */
  REAL4 *hist; /* Values of filter history. */

  INITSTATUS(stat);

  /* Check all the passed parameters for null pointers. */
  ASSERT(output,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->directCoef,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->recursCoef,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->history,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->directCoef->data,stat,IIRFILTERH_ENUL,
	 IIRFILTERH_MSGENUL);
  ASSERT(filter->recursCoef->data,stat,IIRFILTERH_ENUL,
	 IIRFILTERH_MSGENUL);
  ASSERT(filter->history->data,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);

  /* Compute the auxiliary datum. */
  jmax=filter->recursCoef->length;
  coef=filter->recursCoef->data+1;
  hist=filter->history->data;
  for(j=1;j<jmax;j++)
    input+=(*(coef++))*(*(hist++));
  hist-=(jmax-1);

  /* Compute the filter output. */
  jmax=filter->directCoef->length;
  coef=filter->directCoef->data;
  *output=(*(coef++))*input;
  for(j=1;j<jmax;j++)
    *output+=(*(coef++))*(*(hist++));
  hist-=(jmax-1);

  /* Updata the filter history. */
  jmax=filter->history->length-1;
  hist+=jmax;
  for(j=jmax;j>0;j--,hist--)
    *hist=hist[-1];
  *hist=input;

  /* Normal exit */
  RETURN(stat);
}


/**
 * WARNING: THIS FUNCTION IS OBSOLETE.
 * \deprecated
 */
void
LALIIRFilterREAL8( LALStatus      *stat,
		   REAL8          *output,
		   REAL8          input,
		   REAL8IIRFilter *filter )
{
  INITSTATUS(stat);

  /* Check all the passed parameters for null pointers. */
  ASSERT(output,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->directCoef,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->recursCoef,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->history,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->directCoef->data,stat,IIRFILTERH_ENUL,
	 IIRFILTERH_MSGENUL);
  ASSERT(filter->recursCoef->data,stat,IIRFILTERH_ENUL,
	 IIRFILTERH_MSGENUL);
  ASSERT(filter->history->data,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);

  *output=XLALIIRFilterREAL8(input,filter);

  /* Normal exit */
  RETURN(stat);
}


/**
 * WARNING: THIS FUNCTION IS OBSOLETE.
 * \deprecated
 */
REAL4
LALSIIRFilter( REAL4 x, REAL4IIRFilter *filter )
{
  INT4 j;      /* Index for filter coefficients. */
  INT4 jmax;   /* Number of filter coefficients. */
  REAL4 *coef; /* Values of filter coefficients. */
  REAL4 *hist; /* Values of filter history. */
  REAL4 w;     /* Auxiliary datum. */
  REAL4 y;     /* Output datum. */

  /* Compute the auxiliary datum. */
  jmax=filter->recursCoef->length;
  coef=filter->recursCoef->data+1;
  hist=filter->history->data;
  w=x;
  for(j=1;j<jmax;j++)
    w+=(*(coef++))*(*(hist++));
  hist-=(jmax-1);

  /* Compute the filter output. */
  jmax=filter->directCoef->length;
  coef=filter->directCoef->data;
  y=(*(coef++))*w;
  for(j=1;j<jmax;j++)
    y+=(*(coef++))*(*(hist++));
  hist-=(jmax-1);

  /* Updata the filter history. */
  jmax=filter->history->length-1;
  hist+=jmax;
  for(j=jmax;j>0;j--,hist--)
    *hist=hist[-1];
  *hist=w;

  /* Normal exit */
  return y;
}


/** \see See \ref IIRFilter_c for documentation */
REAL8 XLALIIRFilterREAL8( REAL8 x, REAL8IIRFilter *filter )
{
  INT4 j;      /* Index for filter coefficients. */
  INT4 jmax;   /* Number of filter coefficients. */
  REAL8 *coef; /* Values of filter coefficients. */
  REAL8 *hist; /* Values of filter history. */
  REAL8 w;     /* Auxiliary datum. */
  REAL8 y;     /* Output datum. */

  /* Compute the auxiliary datum. */
  jmax=filter->recursCoef->length;
  coef=filter->recursCoef->data+1;
  hist=filter->history->data;
  w=x;
  for(j=1;j<jmax;j++)
    w+=(*(coef++))*(*(hist++));
  hist-=(jmax-1);

  /* Compute the filter output. */
  jmax=filter->directCoef->length;
  coef=filter->directCoef->data;
  y=(*(coef++))*w;
  for(j=1;j<jmax;j++)
    y+=(*(coef++))*(*(hist++));
  hist-=(jmax-1);

  /* Updata the filter history. */
  jmax=filter->history->length-1;
  hist+=jmax;
  for(j=jmax;j>0;j--,hist--)
    *hist=hist[-1];
  *hist=w;

  /* Normal exit */
  return y;
}

/** \see See \ref IIRFilter_c for documentation */
REAL4 XLALIIRFilterREAL4( REAL4 x, REAL8IIRFilter *filter )
{
  return (REAL4)XLALIIRFilterREAL8((REAL8)x,filter);
}
/*@}*/
