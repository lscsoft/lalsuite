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

#ifndef _IIRFILTER_H
#define _IIRFILTER_H

#include <lal/LALStdlib.h>
#include <lal/ZPGFilter.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
   \addtogroup IIRFilter_h
   \author Creighton, T. D.

   \brief Provides routines to make and apply IIR filters.

   \heading{Synopsis}
   \code
   #include <lal/IIRFilter.h>
   \endcode

This header covers routines that create, destroy, and apply
generic time-domain filters, given by objects of type
<tt>\<datatype\>IIRFilter</tt>, where <tt>\<datatype\></tt> is either
\c REAL4 or \c REAL8.

An IIR (Infinite Impulse Response) filter is a generalized linear
causal time-domain filter, in which the filter output \f$y_n=y(t_n)\f$ at
any sampled time \f$t_n=t_0+n\Delta t\f$ is a linear combination of the
input \f$x\f$ \e and output \f$y\f$ at previous sampled times:
\f[
y_n = \sum_{k=0}^M c_k x_{n-k} + \sum_{l=1}^N d_l y_{n-l} \; .
\f]
The coefficients \f$c_k\f$ are called the direct filter coefficients, and
the coefficients \f$d_l\f$ are the recursive filter coefficients.  The
filter order is the larger of \f$M\f$ or \f$N\f$, and determines how far back
in time the filter must look to determine its next output.  However,
the recursive nature of the filter means that the output can depend on
input arbitrarily far in the past; hence the name "infinite impulse
response".  Nonetheless, for a well-designed, stable filter, the
actual filter response to an impulse should diminish rapidly beyond
some characteristic timescale.

Note that nonrecursive FIR (Finite Impulse Response) filters are
considered a subset of IIR filters, having \f$N=0\f$.

For practical implementation, it is convenient to express the bilinear
equation above as two linear equations involving an auxiliary sequence
\f$w\f$:
\f[
w_n = x_n + \sum_{l=1}^N d_l w_{n-l} \; ,
\f]
\f[
y_n = \sum_{k=0}^M c_k w_{n-k} \; .
\f]
The equivalence of this to the first expression is not obvious, but
can be proven by mathematical induction.  The advantage of the
auxiliary variable representation is twofold.  First, when one is
feeding data point by point to the filter, the filter needs only
"remember" the previous \f$M\f$ or \f$N\f$ (whichever is larger) values of
\f$w\f$, rather than remembering the previous \f$M\f$ values of \f$x\f$ \e and
the previous \f$N\f$ values of \f$y\f$.  Second, when filtering a large stored
data vector, the filter response can be computed in place: one first
runs forward through the vector replacing \f$x\f$ with \f$w\f$, and then
backward replacing \f$w\f$ with \f$y\f$.

Although the IIR filters in these routines are explicitly real, one
can consider formally their complex response.  A sinusoidal input can
thus be written as \f$x_n=X\exp(2\pi ifn\Delta t)=Xz^n\f$, where \f$X\f$ is a
complex amplitude and \f$z=\exp(2\pi if\Delta t)\f$ is a complex
parametrization of the frequency.  By linearity, the output must also
be sinusoidal: \f$y_m=Y\exp(2\pi ifm\Delta t)=Yz^m\f$.  Putting these into
the bilinear equation, one can easily compute the filter's complex
transfer function:
\f[
T(z) = \frac{Y}{X} = \frac{\sum_{k=0}^M c_k z^{-k}}
                      {1 - \sum_{l=1}^N d_l z^{-l}}
\f]
This can be readily converted to and from the "zeros, poles, gain"
representation of a filter, which expresses \f$T(z)\f$ as a factored
rational function of \f$z\f$.

It should also be noted that, in the routines covered by this header,
I have adopted the convention of including a redundant recursive
coefficient \f$d_0\f$, in order to make the indexing more intuitive.  For
formal correctness \f$d_0\f$ should be set to \f$-1\f$, although the filtering
routines never actually use this coefficient.

*/
/*@{*/

/**
@{
\defgroup CreateIIRFilter_c 	Module CreateIIRFilter.c
\defgroup DestroyIIRFilter_c 	Module DestroyIIRFilter.c
\defgroup IIRFilter_c 		Module IIRFilter.c
\defgroup IIRFilterVector_c 	Module IIRFilterVector.c
\defgroup IIRFilterVectorR_c 	Module IIRFilterVectorR.c
@}
*/

/** \name Error Codes */
/*@{*/
#define IIRFILTERH_ENUL  1	/**< Unexpected null pointer in arguments */
#define IIRFILTERH_EOUT  2	/**< Output handle points to a non-null pointer */
#define IIRFILTERH_EMEM  3	/**< Memory allocation error */
#define IIRFILTERH_EPAIR 4	/**< Input has unpaired nonreal poles or zeros */
/*@}*/

/** \cond DONT_DOXYGEN */
#define IIRFILTERH_MSGENUL  "Unexpected null pointer in arguments"
#define IIRFILTERH_MSGEOUT  "Output handle points to a non-null pointer"
#define IIRFILTERH_MSGEMEM  "Memory allocation error"
#define IIRFILTERH_MSGEPAIR "Input has unpaired nonreal poles or zeros"
/** \endcond */

/** This structure stores the direct and recursive REAL4 filter coefficients, as
 * well as the history of the auxiliary sequence \f$w\f$.
 * The length of the history vector gives the order of the filter
 */
#ifdef SWIG /* SWIG interface directives */
SWIGLAL(IMMUTABLE_MEMBERS(tagREAL4IIRFilter, name));
#endif /* SWIG */
typedef struct tagREAL4IIRFilter{
  const CHAR *name;        /**< User assigned name. */
  REAL8 deltaT;            /**< Sampling time interval of the filter; If \f$\leq0\f$, it will be ignored (ie it will be taken from the data stream) */
  REAL4Vector *directCoef; /**< The direct filter coefficients. */
  REAL4Vector *recursCoef; /**< The recursive filter coefficients. */
  REAL4Vector *history;    /**< The previous values of w. */
} REAL4IIRFilter;

/** This structure stores the direct and recursive REAL8 filter coefficients, as
 * well as the history of the auxiliary sequence \f$w\f$.
 * The length of the history vector gives the order of the filter
 */
#ifdef SWIG /* SWIG interface directives */
SWIGLAL(IMMUTABLE_MEMBERS(tagREAL8IIRFilter, name));
#endif /* SWIG */
typedef struct tagREAL8IIRFilter{
  const CHAR *name;        /**< User assigned name. */
  REAL8 deltaT;            /**< Sampling time interval of the filter; If \f$\leq0\f$, it will be ignored (ie it will be taken from the data stream). */
  REAL8Vector *directCoef; /**< The direct filter coefficients. */
  REAL8Vector *recursCoef; /**< The recursive filter coefficients. */
  REAL8Vector *history;    /**< The previous values of w. */
} REAL8IIRFilter;

/*@}*/

/* Function prototypes. */
REAL4IIRFilter *XLALCreateREAL4IIRFilter( COMPLEX8ZPGFilter *input );
REAL8IIRFilter *XLALCreateREAL8IIRFilter( COMPLEX16ZPGFilter *input );
void XLALDestroyREAL4IIRFilter( REAL4IIRFilter *filter );
void XLALDestroyREAL8IIRFilter( REAL8IIRFilter *filter );

int XLALIIRFilterREAL4Vector( REAL4Vector *vector, REAL8IIRFilter *filter );
int XLALIIRFilterREAL8Vector( REAL8Vector *vector, REAL8IIRFilter *filter );
int XLALIIRFilterReverseREAL4Vector( REAL4Vector *vector, REAL8IIRFilter *filter );
int XLALIIRFilterReverseREAL8Vector( REAL8Vector *vector, REAL8IIRFilter *filter );

REAL4 XLALIIRFilterREAL4( REAL4 x, REAL8IIRFilter *filter );
REAL8 XLALIIRFilterREAL8( REAL8 x, REAL8IIRFilter *filter );
/* WARNING: THIS FUNCTION IS OBSOLETE */
REAL4 LALSIIRFilter( REAL4 x, REAL4IIRFilter *filter );
/* REAL8 LALDIIRFilter( REAL8 x, REAL8IIRFilter *filter ); */
#define LALDIIRFilter(x,f) XLALIIRFilterREAL8(x,f)



/* ----- CreateIIRFilter.c ---------- */
void
LALCreateREAL4IIRFilter( LALStatus         *status,
			 REAL4IIRFilter    **output,
			 COMPLEX8ZPGFilter *input );

void
LALCreateREAL8IIRFilter( LALStatus          *status,
			 REAL8IIRFilter     **output,
			 COMPLEX16ZPGFilter *input );

/* ----- DestroyIIRFilter.c ---------- */
void
LALDestroyREAL4IIRFilter( LALStatus      *status,
			  REAL4IIRFilter **input );

void
LALDestroyREAL8IIRFilter( LALStatus      *status,
			  REAL8IIRFilter **input );

/* ----- IIRFilter.c ---------- */
void
LALIIRFilterREAL4( LALStatus      *status,
		   REAL4          *output,
		   REAL4          input,
		   REAL4IIRFilter *filter );

void
LALIIRFilterREAL8( LALStatus      *status,
		   REAL8          *output,
		   REAL8          input,
		   REAL8IIRFilter *filter );

/* ----- IIRFilterVector.c ---------- */
void
LALIIRFilterREAL4Vector( LALStatus      *status,
			 REAL4Vector    *vector,
			 REAL4IIRFilter *filter );

void
LALIIRFilterREAL8Vector( LALStatus      *status,
			 REAL8Vector    *vector,
			 REAL8IIRFilter *filter );

void
LALDIIRFilterREAL4Vector( LALStatus      *status,
			  REAL4Vector    *vector,
			  REAL8IIRFilter *filter );

/* ----- IIRFilterVectorR.c ---------- */
void
LALIIRFilterREAL4VectorR( LALStatus      *status,
			  REAL4Vector    *vector,
			  REAL4IIRFilter *filter );

void
LALIIRFilterREAL8VectorR( LALStatus      *status,
			  REAL8Vector    *vector,
			  REAL8IIRFilter *filter );

void
LALDIIRFilterREAL4VectorR( LALStatus      *status,
			   REAL4Vector    *vector,
			   REAL8IIRFilter *filter );


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _IIRFILTER_H */
