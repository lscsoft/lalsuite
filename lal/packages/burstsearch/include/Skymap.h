/*
 * $Id$
 *
 * Copyright (C) 2008 Antony Searle
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
    int hl, hv; /* lv */

    XLALSkymapPixelType* pixel;
    int pixelCount;
} XLALSkymapPlanType;

/* 
 * Construct an analysis plan for a given frequency
 */
XLALSkymapPlanType* XLALSkymapConstructPlan(int sampleFrequency);

/* 
 * Destroy an analysis plan
 */
void XLALSkymapDestroyPlan(XLALSkymapPlanType* plan);

/* 
 * Produce a skymap according to a plan
 */
int XLALSkymapEllipticalHypothesis(XLALSkymapPlanType* plan, double* p, double sigma, double w[3], int begin[3], int end[3], double** x, int* bests);
int XLALSkymapAnalyzeElliptical(double* p, XLALSkymapPlanType* plan, double sigma, double w[3], int n, double** x);

/*
 * Investigate the glitch hypothesis for model selection
 */

int XLALSkymapGlitchHypothesis(XLALSkymapPlanType* plan, double* p, double sigma, double w[3], int begin[3], int end[3], double** x);

/*
 * Add skymaps together, taking account of their log representation and valid
 * pixels
 */
void XLALSkymapSum(XLALSkymapPlanType* plan, double* a, const double* b, const double* c);

/* 
 * Find the most plausible direction (the mode of the distribution) and return
 * its theta and phi coordinates 
 */
void XLALSkymapModeThetaPhi(XLALSkymapPlanType* plan, double* p, double thetaphi[2]);

/*
 * Render the skymap from the internal format to a variety of map projections
 */
int XLALSkymapRenderEqualArea(int m, int n, double* q, XLALSkymapPlanType* plan, double* p);

int XLALSkymapRenderEquirectangular(int m, int n, double* q, XLALSkymapPlanType* plan, double* p);

int XLALSkymapRenderMollweide(int m, int n, double* q, XLALSkymapPlanType* plan, double* p);

#ifdef __cplusplus
}
#endif

#endif /* SKYMAP_H */

