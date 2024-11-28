//
// Copyright (C) 2019 John T. Whelan
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,x5
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
//  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//  MA  02110-1301  USA
//

#ifndef _TASCPORBTILING_H
#define _TASCPORBTILING_H


/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/*---------- INCLUDES ----------*/
#include <lal/LALDatatypes.h>
#include <lal/LatticeTiling.h>

/*---------- DEFINES ----------*/


///
/// Set an orbital period as a function of time of ascension
///
int XLALSetLatticeTilingPorbEllipticalBound(
  LatticeTiling *tiling,        ///< [in] Lattice tiling
  const size_t tasc_dimension,  ///< [in] Time of ascension dimension
  const size_t porb_dimension,  ///< [in] Orbital period dimension
  const double P0,              ///< [in] Most likely orbital period
  const double sigP,            ///< [in] One-sigma uncertainty on orbital period
  const double T0,              ///< [in] Most likely time of ascension (uncorrelated with orbital period
  const double sigT,            ///< [in] One-sigma uncertainty on time of ascension
  const int norb,            ///< [in] Number of orbits between time of ascention estimate and search region
  const double nsigma,           ///< [in] Radius in sigma of circular search region in parameter space scaled by uncertainties
  const BOOLEAN useShearedPeriod  ///< [in] Whether to use sheared Porb coordinate so the centerline of the search ellipse is horizontal
);

#ifdef __cplusplus
}
#endif

#endif // _TASCPORBTILING_H

// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
