//
// Copyright (C) 2019--2023 Benjamin Grace
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307 USA
//

#ifndef _PIECEWISEMODEL_H
#define _PIECEWISEMODEL_H

#include <stdio.h>
#include <math.h>
#include <lal/LatticeTiling.h>

#ifdef __cplusplus
extern "C" {
#endif

///
/// Sets the bounds for the piecewise model
///
int XLALSetLatticeTilingPiecewiseBounds(
  LatticeTiling *tiling,        ///< Lattice tiling
  const size_t s,               ///< Number of frequency/spindown parameters per knot
  const double fmin,            ///< Minimum initial frequency
  const double fmax,            ///< Maximum initial frequency
  const double nmin,            ///< Minimum braking index
  const double nmax,            ///< Maximum braking index
  const double kmin,            ///< Minimum k value
  const double kmax,            ///< Maximum k value
  const gsl_vector *knots,      ///< List of knots
  const gsl_vector *bboxpad,    ///< Vector containing fractional bounding box padding
  const gsl_vector_int* intpad  ///< Vector containing number of integer points to use for padding
);

#ifdef __cplusplus
}
#endif

#endif //_PIECEWISEMODEL_H
