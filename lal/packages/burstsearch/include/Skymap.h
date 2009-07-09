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

/*
 * Stores pre-computed properties of each sky direction
 */
typedef struct tagLALSkymapPixelType
{
    double area;
    double direction[3];
    double f[3][2];
} XLALSkymapPixelType;

/*
 * Stores pre-computed information for generating a sky map
 */
typedef struct tagXLALSkymapPlanType
{
    int sampleFrequency;
    LALDetector site[3];
    double siteNormal[3];

    /*  physical delays between the detectors,
        rounded to the nearest integer sample,
        range over [-hl, +hl] and [-hv, +hv] */
    int hl, hv, lv;

    XLALSkymapPixelType* pixel;
    int pixelCount;

    /* raster size */
    int m;
    int n;
} XLALSkymapPlanType;

/*
 * Construct an analysis plan for a given frequency
 */
XLALSkymapPlanType* XLALSkymapConstructPlanMN(int sampleFrequency, int m, int n);

/*
 * Destroy an analysis plan
 */
void XLALSkymapDestroyPlan(XLALSkymapPlanType* plan);

/*
 * Produce a skymap according to a plan
 */
int XLALSkymapSignalHypothesisWithLimits(XLALSkymapPlanType* plan, double* p, double sigma, double w[3], int begin[3], int end[3], double** x, int *counts, int *modes, int delay_limits[8]);

/* deprecated legacy interface */
int XLALSkymapSignalHypothesis(XLALSkymapPlanType* plan, double* p, double sigma, double w[3], int begin[3], int end[3], double** x, int *counts, int *modes);
/* deprecated legacy interface */
int XLALSkymapEllipticalHypothesis(XLALSkymapPlanType* plan, double* p, double sigma, double w[3], int begin[3], int end[3], double** x, int* bests);
/* deprecated legacy interface*/
int XLALSkymapAnalyzeElliptical(double* p, XLALSkymapPlanType* plan, double sigma, double w[3], int n, double** x);

/*
 * Investigate the glitch hypothesis for model selection
 */

int XLALSkymapGlitchHypothesis(XLALSkymapPlanType* plan, double* p, double sigma, double w[3], int begin[3], int end[3], double** x);

/*
 * Manipulate log-represented values together, preventing overflow
 */
double XLALSkymapLogSumExp(double a, double b);
double XLALSkymapLogDifferenceExp(double a, double b);
double XLALSkymapLogTotalExp(double* begin, double* end);

/*
 * Add skymaps together, taking account of their log representation and valid
 * pixels
 */
void XLALSkymapSum(XLALSkymapPlanType* plan, double* a, const double* b, const double* c);

/*
 * Lightweight coordinate transformations
 */
void XLALSkymapCartesianFromSpherical(double a[3], double b[2]);
void XLALSkymapSphericalFromCartesian(double a[2], double b[3]);
void XLALSkymapDelaysFromDirection(XLALSkymapPlanType* plan, int delays[3], double direction[3]);
int XLALSkymapIndexFromDirection(XLALSkymapPlanType* plan, double direction[3]);

/*
 * Find the most plausible direction (the mode of the distribution) and return
 * its theta and phi coordinates
 */
void XLALSkymapModeThetaPhi(XLALSkymapPlanType* plan, double* p, double thetaphi[2]);

/*
 * Render the skymap from the internal format to a variety of map projections
 */
int XLALSkymapRender(double* q, XLALSkymapPlanType* plan, double* p);

///////////////////////////////////////////////////////////////////////////

// VERSION 2


typedef double XLALSkymap2SphericalPolarType[2];

//  Stores pre-computed properties of each sky direction

typedef struct
{
    double f[3][2];
    int delay[3];
} XLALSkymap2DirectionPropertiesType;

// Stores pre-computed kernel information for each sky direction

typedef struct
{
    double k[3][3];
    double logNormalization;
} XLALSkymap2KernelType;

// Stores pre-computed information for generating a sky map

typedef struct
{
    int sampleFrequency;
    LALDetector site[3];
} XLALSkymap2PlanType;

// Initialize sample frequency and HLV site information

void XLALSkymap2PlanConstruct(
    int sampleFrequency,
    XLALSkymap2PlanType* plan
    );

// Turn directions into sample delays and antenna patterns

void XLALSkymap2DirectionPropertiesConstruct(
    XLALSkymap2PlanType* plan,
    XLALSkymap2SphericalPolarType* direction,
    XLALSkymap2DirectionPropertiesType* properties
    );

// Turn antenna patterns and waveform normalization into kernel matrix

void XLALSkymap2KernelConstruct(
    XLALSkymap2DirectionPropertiesType* properties,
    double wSw[3],
    XLALSkymap2KernelType* kernel
    );

// Apply the kernel to some data

void XLALSkymap2Apply(
    XLALSkymap2DirectionPropertiesType* properties,
    XLALSkymap2KernelType* kernel,
    double* xSw[3],
    int tau,
    double* logPosterior
    );




#ifdef __cplusplus
}
#endif

#endif /* SKYMAP_H */

