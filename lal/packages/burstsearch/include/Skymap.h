/*
 * $Id$
 *
 * Copyright (C) 2008-9 Antony Searle
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

#ifndef SKYMAP_H
#define SKYMAP_H

#ifdef __cplusplus
extern "C" {
#endif

#include <lal/LALDetectors.h>
#include <lal/LALRCSID.h>
NRCSID(EPSEARCHH, "$Id$");

// Manipulate log-represented values without overflowing

double XLALSkymapLogSumExp(double a, double b);
double XLALSkymapLogDifferenceExp(double a, double b);
double XLALSkymapLogTotalExp(double* begin, double* end);

// Lightweight coordinate transformations

void XLALSkymapCartesianFromSpherical(double a[3], double b[2]);
void XLALSkymapSphericalFromCartesian(double a[2], double b[3]);

// Cubic interpolation

void XLALSkymap2InterpolationWeights(double t, double* w);
double XLALSkymap2Interpolate(double* x, double* w);

// Largest network of interest (needed to allocate storage)

#define XLALSKYMAP2_N 5

// Use detector names from lal/packages/tools/include/LALDetectors.h :
//	LAL_TAMA_300_DETECTOR	=	0
//	LAL_VIRGO_DETECTOR      =	1
//	LAL_GEO_600_DETECTOR	=	2
//	LAL_LHO_2K_DETECTOR     =	3
//	LAL_LHO_4K_DETECTOR     =	4
//	LAL_LLO_4K_DETECTOR     =	5

// Struct to store basic properties of the analysis: the sample rate, and
// the detectors involved (described by their LAL ID numbers)

typedef struct
{
    int sampleFrequency;
    int n;
    LALDetector site[XLALSKYMAP2_N];
} XLALSkymap2PlanType;

void XLALSkymap2PlanConstruct(
    int sampleFrequency,
    int n,
    int* detectors,
    XLALSkymap2PlanType* plan
    );

// Struct to store reuseable pre-computed quantities for a specific
// direction, set of detectors, and sample rate

typedef struct
{
    double f[XLALSKYMAP2_N][2];
    int delay[XLALSKYMAP2_N];
	double weight[XLALSKYMAP2_N][4];
} XLALSkymap2DirectionPropertiesType;

void XLALSkymap2DirectionPropertiesConstruct(
    XLALSkymap2PlanType* plan,
    double* direction,
    XLALSkymap2DirectionPropertiesType* properties
    );

// Struct to store reuseable pre-computed kernel for a specific direction,
// power spectra, and sample rate

typedef struct
{
    double k[XLALSKYMAP2_N][XLALSKYMAP2_N];
    double logNormalization;
} XLALSkymap2KernelType;

void XLALSkymap2KernelConstruct(
    XLALSkymap2PlanType* plan,
    XLALSkymap2DirectionPropertiesType* properties,
    double* wSw,
    XLALSkymap2KernelType* kernel
    );

// Compute the Bayesian marginalization integral for the specified system

void XLALSkymap2Apply(
    XLALSkymap2PlanType* plan,
    XLALSkymap2DirectionPropertiesType* properties,
    XLALSkymap2KernelType* kernel,
    double** xSw,
    int tau,
    double* logPosterior
    );

#ifdef __cplusplus
}
#endif

#endif // SKYMAP_H

