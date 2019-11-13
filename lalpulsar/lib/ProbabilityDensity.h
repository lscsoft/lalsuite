/*
 * Copyright (C) 2010 Reinhard Prix
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

/*********************************************************************************/
/**
 * \author R. Prix
 * \file
 * \brief
 * Header file containing the exported API for the ProbabilityDensity Module.
 *
 */

#ifndef _PROBABILITY_DENSITY_H
#define _PROBABILITY_DENSITY_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/* ---------- System includes ---------- */
#include <math.h>

/* gsl-includes */
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* LAL-includes */
#include <lal/LALDatatypes.h>

/* ---------- exported API defines ---------- */

/* ---------- exported API types ---------- */

/**
 * Encode a pdf(x) as a discretized function probDens[i] = pdf( x[i] ) with user-specified bins xBin[i].
 * We store the actual probability *density values, such that prob(x in xBin[i]) = probDens[i] * xBin[i].
 *
 * NOTE: Allows for some special encodings for simplicity and efficiency:
 * - one x0 known with certainty:  pdf(x) = delta(x-x0): ==> xTics = {x0},         probDens=NULL
 * - uniform pdf over [xMin,xMax]: pdf(x) = const.       ==> xTics = {xMin, xMax}, probDens=NULL
 *
 * NOTE2: the optional field 'sampling' is used to buffer gsl-precomputed internals to allow using
 * gsl_ran_discrete() to efficiently draw samples from that distribution (cost ~ O(1)).
 * This field is only for internal use, and will be created automatically by XLALDrawFromPDF1D()
 * if it's not initialized yet.
 *
 */
struct tagpdf1D_t
{
  REAL8Vector *xTics;		/**< N+1-dim vector of ordered x 'tics', i.e. bin-boundaries {x[0], x[1], x[2], ... x[N]} */
  REAL8Vector *probDens;	/**< N-dim vector of binned probability densities probDens[i] = prob( x in [ x[i],x[i+i] )/xBin[i]  */
  BOOLEAN isNormalized;		/**< true if the prob is normalized, ie 1 = int P(x) dx ~ sum_i probDens[i] xBin[i] */
  gsl_ran_discrete_t *sampling;	/**< internal: buffer preprocessed sampling distribution for drawing samples using gsl_ran_discrete() */
};


/**
 * Probability-density function object type.
 * This could be made into an *opaque* type, such that only methods defined in ProbabilityDensity.c
 * can operate on the internals of such objects. Everyone would only be able to pass them around.
 * This requires lots more methods though to be useful, so for now this move is postponed.
 */
typedef struct tagpdf1D_t pdf1D_t;


/* empty struct initializers */


/* ---------- exported API prototypes ---------- */
pdf1D_t *XLALCreateSingularPDF1D ( REAL8 x0 );
pdf1D_t *XLALCreateUniformPDF1D ( REAL8 xMin, REAL8 xMax );
pdf1D_t *XLALCreateDiscretePDF1D ( REAL8 xMin, REAL8 xMax, UINT4 numBins );

REAL8 XLALDrawFromPDF1D ( pdf1D_t *pdf, const gsl_rng *rng );
int XLALCheckValidPDF1D ( const pdf1D_t *pdf );
int XLALNormalizePDF1D ( pdf1D_t *pdf );

REAL8 XLALFindModeOfPDF1D ( const pdf1D_t *pdf );

int XLALOutputPDF1D_to_fp ( FILE* fp, const pdf1D_t *pdf, const char *name );


void XLALDestroyPDF1D ( pdf1D_t *pdf );


#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
