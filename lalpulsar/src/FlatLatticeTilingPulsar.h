//
// Copyright (C) 2007, 2008, 2012 Karl Wette
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA  02111-1307  USA
//

///
/// \addtogroup FlatLatticeTilingPulsar_h
/// \author Karl Wette
/// \brief Lattice-based template generation for continuous wave parameter spaces
///

/// @{

#ifndef _FLATLATTICETILINGPULSAR_H
#define _FLATLATTICETILINGPULSAR_H

#include <lal/LALStdlib.h>
#include <lal/FlatLatticeTiling.h>

#ifdef __cplusplus
extern "C" {
#endif

///
/// Types of tiling for XLALSetFlatLatticeSuperSkyNZBound()
///
typedef enum {
  FLSSNZ_LOWER = 0,	///< Tile over the lower hemisphere
  FLSSNZ_PLANE,		///< Tile over the azimuthal plane
  FLSSNZ_UPPER,		///< Tile over the upper hemisphere
  FLSSNZ_SPHERE,	///< Tile over the whole sphere
  FLSSNZ_LAST
} FLSSNZ;

///
/// Set a singular bound(s) on the Z dimension of the sky position in super-sky coordinates
///
int XLALSetFlatLatticeSuperSkyNZBound(
  FlatLatticeTiling* tiling,		///< [in] Tiling state
  const size_t nz_dimension,		///< [in] Sky position Z dimension
  const FLSSNZ type			///< [in] Tiling type
  );

///
/// Set a constant frequency/spindown parameter space bound, given by the
/// minimum and maximum of the two supplied bounds, on the flat lattice tiling
///
int XLALSetFlatLatticeFnDotConstantBound(
  FlatLatticeTiling* tiling,	///< [in] Tiling state
  const size_t nx_dimension,	///< [in] Sky position X dimension
  const double offset[3],	///< [in] Sky position offset vector
  const size_t dimension,	///< [in] Dimension on which bound applies
  const double bound1,		///< [in] First bound on dimension
  const double bound2		///< [in] Second bound on dimension
  );

///
/// Set a first spindown bound derived from spindown age and braking indices
///
int XLALSetFlatLatticeF1DotAgeBrakingBound(
  FlatLatticeTiling* tiling,		///< [in] Tiling state
  const size_t freq_dimension,		///< [in] Frequency dimension
  const size_t f1dot_dimension,		///< [in] First spindown dimension
  const double age,			///< [in] Spindown age
  const double min_braking,		///< [in] Minimum braking index
  const double max_braking		///< [in] Maximum braking index
  );

///
/// Set a second spindown bound derived from braking indices
///
int XLALSetFlatLatticeF2DotBrakingBound(
  FlatLatticeTiling* tiling,		///< [in] Tiling state
  const size_t freq_dimension,		///< [in] Frequency dimension
  const size_t f1dot_dimension,		///< [in] First spindown dimension
  const size_t f2dot_dimension,		///< [in] Second spindown dimension
  const double min_braking,		///< [in] Minimum braking index
  const double max_braking		///< [in] Maximum braking index
  );

#ifdef __cplusplus
}
#endif

#endif

/// @}
