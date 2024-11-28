/*
 *  Copyright (C) 2024 Leo Singer
 *
 *  Provide a workaround for removal of gsl_bspline_eval_nonzero in GSL 2.8.
 *  The function gsl_bspline_eval_nonzero was replaced with gsl_bspline_basis.
 *  FIXME: remove this file once we require GSL >= 2.8.
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

#ifndef _GSL_BPSLINE_BASIS_H
#define _GSL_BPSLINE_BASIS_H

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_version.h>

#if GSL_MAJOR_VERSION < 2 || (GSL_MAJOR_VERSION == 2 && GSL_MINOR_VERSION < 8)
static int gsl_bspline_basis(const double x, gsl_vector *Bk, size_t *istart, gsl_bspline_workspace *w) {
    size_t iend;
    return gsl_bspline_eval_nonzero(x, Bk, istart, &iend, w);
}
#endif

#endif
