//
// Copyright (C) 2015 Karl Wette
// Copyright (C) 2014 Reinhard Prix
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
// Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA  02110-1301  USA
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <lal/ComputeFstat.h>
#include <lal/Factorial.h>
#include <lal/SinCosLUT.h>

///
/// \file ComputeFstat_DemodHL_Generic.c
/// \ingroup ComputeFstat_Demod_c
/// \brief Old (pre-Akos) LALDemod hotloop variant (unrestricted Dterms)
///
/// \snippet ComputeFstat_DemodHL_Generic.i hotloop
///

#define FUNC XLALComputeFaFb_Generic
#define HOTLOOP_SOURCE "ComputeFstat_DemodHL_Generic.i"
#include "ComputeFstat_Demod_ComputeFaFb.c"
