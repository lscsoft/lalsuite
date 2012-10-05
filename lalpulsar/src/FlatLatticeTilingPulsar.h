/*
 *  Copyright (C) 2007, 2008 Karl Wette
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

/**
 * \author Karl Wette
 * \file
 * \brief Pulsar-specific routines for FlatLatticeTiling
 */

#ifndef _FLATLATTICETILINGPULSAR_H
#define _FLATLATTICETILINGPULSAR_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <lal/LALDatatypes.h>
#include <lal/FlatLatticeTiling.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Functions
 */
int XLALSetFlatLatticeTilingAgeBrakingIndexBounds(FlatLatticeTiling*, double, double, double, double, double);

#ifdef __cplusplus
}
#endif

#endif
