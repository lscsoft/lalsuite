/*
 *  Copyright (C) 2019 John T. Whelan
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

/**
 * \author John T. Whelan
 * \date 2019
 * \file
 * \brief Functions for tiling the orbital phase and period space for
 * Sco X-1 and other binaries
 */


/*---------- INCLUDES ----------*/

#include <config.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_roots.h>

#include <lal/TascPorbTiling.h>
#include <lal/LALConstants.h>
#include <lal/LogPrintf.h>

#include <lal/GSLHelpers.h>

/*---------- DEFINES ----------*/

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

// Square of a number
#define SQR(x)       ((x)*(x))

// Real part of square root of a number, i.e. zero if x < ~0
#define RE_SQRT(x)   sqrt(GSL_MAX(DBL_EPSILON, (x)))

typedef struct {
  size_t tasc_dim;
  int norb;
  double T0p;
  double P0;
  double sigTp;
  double pmsigT;
  double sigP;
  double ksq;
  BOOLEAN sheared;
} PorbEllipticalBoundInfo;

static double PorbEllipticalBound(
  const void *data,
  const size_t dim UNUSED,
  const gsl_matrix *cache UNUSED,
  const gsl_vector *point
)
{

  // Get bounds info
  const PorbEllipticalBoundInfo *info = ( const PorbEllipticalBoundInfo * )data;

  // Get current value of Tasc
  const double Tasc = gsl_vector_get( point, info->tasc_dim );

  double Tscaled = ( Tasc - info->T0p ) / info->sigTp;
  double Pscaled = info->pmsigT * RE_SQRT( info->ksq - SQR( Tscaled ) ) / info->sigTp;
  if ( !( info->sheared ) ) {
    Pscaled += info->norb * info->sigP * Tscaled / info->sigTp;
  }
  // fprintf(stdout,"n=%d\n",info->norb);
  // fprintf(stdout,"sigP=%g\n",info->sigP);
  // fprintf(stdout,"sigTp=%g\n",info->sigTp);
  // fprintf(stdout,"pmsigT=%g\n",info->pmsigT);
  // fprintf(stdout,"Tscaled=%g\n",Tscaled);
  // fprintf(stdout,"Pscaled=%g\n",Pscaled);
  // fprintf(stdout,"Tasc = %.0f + %g\n",info->T0p,(Tasc-info->T0p));
  // fprintf(stdout,"+/-sigT=%g\n", info->pmsigT);
  // fprintf(stdout,"n*sigP*Tscaled/sigTp=%g\n",
  //      info->norb * info->sigP * Tscaled / info->sigTp);
  // fprintf(stdout,"+/-sigT*sqrt(disc)/sigTp=%g\n",
  //      info->pmsigT * RE_SQRT( info->ksq - SQR(Tscaled) ) / info->sigTp);
  // fprintf(stdout,"n*sigP/sigTp=%g\n",info->norb * info->sigP / info->sigTp);
  // fprintf(stdout,"+/-sigT/sigTp=%g\n",info->pmsigT / info->sigTp);
  // fprintf(stdout,"Porb bound = %g + %g\n",info->P0,info->sigP*Pscaled);

  // Return bounds on Porb
  return info->P0 + info->sigP * Pscaled;

}

int XLALSetLatticeTilingPorbEllipticalBound(
  LatticeTiling *tiling,
  const size_t tasc_dimension,
  const size_t porb_dimension,
  const double P0,
  const double sigP,
  const double T0,
  const double sigT,
  const int norb,
  const double nsigma,
  const BOOLEAN useShearedPeriod
)
{
  // Check input
  XLAL_CHECK( tiling != NULL, XLAL_EFAULT );
  XLAL_CHECK( tasc_dimension < porb_dimension, XLAL_EINVAL );
  XLAL_CHECK( nsigma >= 0.0, XLAL_EINVAL );
  XLAL_CHECK( P0 > 0.0, XLAL_EINVAL );
  XLAL_CHECK( sigP >= 0.0, XLAL_EINVAL );
  XLAL_CHECK( T0 > 0.0, XLAL_EINVAL );
  XLAL_CHECK( sigT >= 0.0, XLAL_EINVAL );

  // Set the parameter-space bound
  PorbEllipticalBoundInfo XLAL_INIT_DECL( info_lower );
  PorbEllipticalBoundInfo XLAL_INIT_DECL( info_upper );
  info_lower.tasc_dim = info_upper.tasc_dim = tasc_dimension;
  info_lower.norb = info_upper.norb = norb;
  info_lower.T0p = info_upper.T0p = T0 + norb * P0;
  info_lower.P0 = info_upper.P0 = P0;
  info_lower.sigTp = info_upper.sigTp = sqrt( SQR( sigT ) + SQR( norb ) * SQR( sigP ) );
  info_lower.pmsigT = -sigT;
  info_upper.pmsigT = sigT;
  info_lower.sigP = info_upper.sigP = sigP;
  info_lower.ksq = info_upper.ksq = SQR( nsigma );
  info_lower.sheared = info_upper.sheared = useShearedPeriod;

  XLAL_CHECK( XLALSetLatticeTilingBound( tiling, porb_dimension, PorbEllipticalBound, sizeof( info_lower ), &info_lower, &info_upper ) == XLAL_SUCCESS, XLAL_EFAILED );

  return XLAL_SUCCESS;

}
