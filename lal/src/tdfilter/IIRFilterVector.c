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
 * \addtogroup IIRFilterVector_c
 * \author Creighton, T. D.
 *
 * \brief Applies an IIR filter to a data stream.
 *
 * ### Description ###
 *
 * These functions apply a generic time-domain filter given by an object
 * <tt>*filter</tt> of type \c REAL4IIRFilter or \c REAL8IIRFilter
 * to a list <tt>*vector</tt> of data representing a time series.  This is
 * done in place using the auxiliary data series formalism described in
 * \ref IIRFilter.h, so as to accomodate potentially large data series.
 * To filter a piece of a larger dataset, the calling routine may pass a
 * vector structure whose data pointer and length fields specify a subset
 * of a larger vector.
 *
 * The routine <tt>LALDIIRFilterREAL4Vector()</tt> applies a
 * double-precision filter to single-precision data.  It makes a single
 * pass through the data, continuously updating the filter history at
 * each step rather than storing the auxiliary array in-place.  This
 * reduces roundoff error by keeping \e all intermediate results to
 * double-precision.
 *
 * ### Algorithm ###
 *
 * The implementation of <tt>LALDIIRFilterREAL4Vector()</tt> not only has
 * lower truncation errors than <tt>LALIIRFilterREAL4Vector()</tt>, but
 * also appears to be more computationally efficient, for reasons I have
 * not yet determined; see the documentation for \ref IIRFilterTest.c.
 * These combine to suggest that <tt>LALDIIRFilterREAL4Vector()</tt> is the
 * better overall algorithm for filtering \c REAL4Vectors.
 *
 */
/*@{*/

#undef COMPLEX_DATA
#undef SINGLE_PRECISION

#define COMPLEX_DATA
#define SINGLE_PRECISION
#include "IIRFilterVector_source.c"
#undef SINGLE_PRECISION
#include "IIRFilterVector_source.c"
#undef COMPLEX_DATA
#define SINGLE_PRECISION
#include "IIRFilterVector_source.c"
#undef SINGLE_PRECISION
#include "IIRFilterVector_source.c"

/**
 * WARNING: THIS FUNCTION IS OBSOLETE.
 * \deprecated
 */
void
LALIIRFilterREAL4Vector( LALStatus      *stat,
			 REAL4Vector    *vector,
			 REAL4IIRFilter *filter )
{
  INT4 i;            /* Loop counter for data vector. */
  INT4 j;            /* Index for filter coeficients. */
  INT4 k;            /* Index for filter history. */
  INT4 length;       /* Length of vector. */
  REAL4 *data;       /* Vector data. */
  REAL8 datum;       /* Temporary working variable. */
  INT4 directOrder;  /* Number of direct filter coefficients. */
  INT4 recursOrder;  /* Number of recursive filter coefficients. */
  INT4 numHist;      /* The number of history data. */
  REAL4 *directCoef; /* Direct filter coefficients. */
  REAL4 *recursCoef; /* Recursive filter coefficients. */
  REAL4 *history;    /* Filter history. */
  REAL4 *temp=NULL;  /* Temporary storage for the filter history. */

  INITSTATUS(stat);

  /* Make sure all the structures have been initialized. */
  ASSERT(vector,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(vector->data,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->directCoef,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->recursCoef,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->history,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->directCoef->data,stat,IIRFILTERH_ENUL,
      IIRFILTERH_MSGENUL);
  ASSERT(filter->recursCoef->data,stat,IIRFILTERH_ENUL,
      IIRFILTERH_MSGENUL);
  ASSERT(filter->history->data,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  length=vector->length;
  data=vector->data;
  directOrder=filter->directCoef->length;
  recursOrder=filter->recursCoef->length;
  numHist=filter->history->length;
  directCoef=filter->directCoef->data;
  recursCoef=filter->recursCoef->data;
  history=filter->history->data;
  temp=(REAL4 *)LALMalloc(numHist*sizeof(REAL4));
  if ( !temp ) {
    ABORT(stat,IIRFILTERH_EMEM,IIRFILTERH_MSGEMEM);
  }

  /* Compute the auxiliary data series. */
  for(i=0;(i<recursOrder)&&(i<length);i++,data++){
    datum=*data;
    for(j=1;j<=i;j++)
      datum+=data[-j]*recursCoef[j];
    for(k=0;j<recursOrder;j++,k++)
      datum+=history[k]*recursCoef[j];
    *data=datum;
  }
  for(;i<length;i++,data++){
    datum=*data;
    for(j=1;j<recursOrder;j++)
      datum+=data[-j]*recursCoef[j];
    *data=datum;
  }
  data--;

  /* Store the last few auxiliary data to the temporary history. */
  for(k=numHist-1;k>=length;k--)
    temp[k]=history[k-length];
  for(;k>=0;k--)
    temp[k]=data[-k];

  /* Compute the output data series. */
  for(;i>directOrder;i--,data--){
    datum=*data*directCoef[0];
    for(j=1;j<directOrder;j++)
      datum+=data[-j]*directCoef[j];
    *data=datum;
  }
  for(;i>0;i--,data--){
    datum=*data*directCoef[0];
    for(j=1;j<i;j++)
      datum+=data[-j]*directCoef[j];
    for(k=0;j<directOrder;j++,k++)
      datum+=history[k]*directCoef[j];
    *data=datum;
  }

  /* Update the filter history from the temporary history. */
  for(k=0;k<numHist;k++)
    history[k]=temp[k];
  LALFree(temp);

  /* Normal exit */
  RETURN(stat);
}


/**
 * WARNING: THIS FUNCTION IS OBSOLETE.
 * \deprecated
 */
void
LALIIRFilterREAL8Vector( LALStatus      *stat,
    REAL8Vector    *vector,
    REAL8IIRFilter *filter )
{
  INITSTATUS(stat);

  /* Make sure all the structures have been initialized. */
  ASSERT(vector,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(vector->data,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->directCoef,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->recursCoef,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->history,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->directCoef->data,stat,IIRFILTERH_ENUL,
      IIRFILTERH_MSGENUL);
  ASSERT(filter->recursCoef->data,stat,IIRFILTERH_ENUL,
      IIRFILTERH_MSGENUL);
  ASSERT(filter->history->data,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);

  if ( XLALIIRFilterREAL8Vector( vector, filter ) < 0 )
  {
    XLALClearErrno();
    ABORT(stat,IIRFILTERH_EMEM,IIRFILTERH_MSGEMEM);
  }

  /* Normal exit */
  RETURN(stat);
}


/**
 * WARNING: THIS FUNCTION IS OBSOLETE.
 * \deprecated
 */
void
LALDIIRFilterREAL4Vector( LALStatus      *stat,
			  REAL4Vector    *vector,
			  REAL8IIRFilter *filter )
{
  INITSTATUS(stat);

  /* Make sure all the structures have been initialized. */
  ASSERT(vector,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(vector->data,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->directCoef,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->recursCoef,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->history,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->directCoef->data,stat,IIRFILTERH_ENUL,
	 IIRFILTERH_MSGENUL);
  ASSERT(filter->recursCoef->data,stat,IIRFILTERH_ENUL,
	 IIRFILTERH_MSGENUL);
  ASSERT(filter->history->data,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);

  if ( XLALIIRFilterREAL4Vector( vector, filter ) < 0 )
  {
    XLALClearErrno();
    ABORT(stat,IIRFILTERH_EMEM,IIRFILTERH_MSGEMEM);
  }

  /* Normal exit */
  RETURN(stat);
}
/*@}*/
