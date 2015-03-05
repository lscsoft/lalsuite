/* Copyright (C) 2012 Walter Del Pozzo, Tjonnie Li
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
 
#ifndef LALCOSMOLOGYCALCULATOR_H
#define LALCOSMOLOGYCALCULATOR_H

#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include <math.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_integration.h>

typedef struct tagLALCosmologicalParameters
{
    double h;   /* normalised hubble constant h = H0/100km/Mpc/s */
    double om;  /* matter energy density */
    double ol;  /* cosmological constant density */
    double ok;  /* curvature */
    double w0;  /* 0th order dark energy equation of state parameter */
    double w1;  /* 1st order dark energy equation of state parameter */
    double w2;  /* 2nd order dark energy equation of state parameter */
}   LALCosmologicalParameters;

typedef struct tagLALCosmologicalRateParameters
{
    double r0;  /* local coalescence rate in units of Mpc^-3 yr^-1 */ 
    double W;           /* phenomenological parameter for the analytical fit to the SFR rate */
    double Q;           /* phenomenological parameter for the analytical fit to the SFR rate */
    double R;           /* phenomenological parameter for the analytical fit to the SFR rate */
} LALCosmologicalRateParameters;

typedef struct tagLALCosmologicalParametersAndRate
{
    LALCosmologicalParameters *omega;
    LALCosmologicalRateParameters *rate;
} LALCosmologicalParametersAndRate;

double XLALLuminosityDistance(
            LALCosmologicalParameters *omega, 
            double z);

double XLALAngularDistance(
            LALCosmologicalParameters *omega, 
            double z);

double XLALComovingLOSDistance(
            LALCosmologicalParameters *omega, 
            double z);

double XLALComovingTransverseDistance(
            LALCosmologicalParameters *omega, 
            double z);            

double XLALHubbleDistance(
            LALCosmologicalParameters *omega
            );
            
double XLALHubbleParameter(double z,
            void *omega
            );
            
double XLALIntegrateHubbleParameter(
            LALCosmologicalParameters *omega, 
            double z);

double XLALUniformComovingVolumeDensity(
            double z,
            void *omega);

double XLALUniformComovingVolumeDistribution(
            LALCosmologicalParameters *omega, 
            double z,
            double zmax);

double XLALIntegrateComovingVolumeDensity(LALCosmologicalParameters *omega, double z);
            
LALCosmologicalParameters *XLALCreateCosmologicalParameters(double h, double om, double ol, double w0, double w1, double w2);

void XLALDestroyCosmologicalParameters(LALCosmologicalParameters *omega);

double XLALGetHubbleConstant(LALCosmologicalParameters *omega);
double XLALGetOmegaMatter(LALCosmologicalParameters *omega);
double XLALGetOmegaLambda(LALCosmologicalParameters *omega);
double XLALGetOmegaK(LALCosmologicalParameters *omega);
double XLALGetW0(LALCosmologicalParameters *omega);
double XLALGetW1(LALCosmologicalParameters *omega);
double XLALGetW2(LALCosmologicalParameters *omega);
void XLALSetCosmologicalParametersDefaultValue(LALCosmologicalParameters *omega);


LALCosmologicalRateParameters *XLALCreateCosmologicalRateParameters(double r0, double W, double Q, double R);
void XLALDestroyCosmologicalRateParameters(LALCosmologicalRateParameters *rate);
double XLALGetLocalRate(LALCosmologicalRateParameters *rate);
double XLALStarFormationDensity(double z, void *rate);

double XLALRateWeightedUniformComovingVolumeDensity(double z, void *params);
double XLALIntegrateRateWeightedComovingVolumeDensity(LALCosmologicalParametersAndRate *p, double z);
double XLALRateWeightedComovingVolumeDistribution(LALCosmologicalParametersAndRate *p, double z, double zmax);

LALCosmologicalParametersAndRate *XLALCreateCosmologicalParametersAndRate(void);
void XLALDestroyCosmologicalParametersAndRate(LALCosmologicalParametersAndRate *p);
void XLALSetCosmologicalRateParametersDefaultValue(LALCosmologicalRateParameters *params);
#endif

