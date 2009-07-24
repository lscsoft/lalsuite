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

// include guards deliberately omitted

#ifdef SKYMAP_N

#define CAT2(A, B) A##B
#define CAT(A, B) CAT2(A, B)

#ifdef __cplusplus
extern "C" {
#endif

//  Stores pre-computed properties of each sky direction

typedef struct
{
    double f[SKYMAP_N][2];
    int delay[SKYMAP_N];
} CAT(XLALSkymap2DirectionPropertiesType, SKYMAP_N);

// Stores pre-computed kernel information for each sky direction

typedef struct
{
    double k[SKYMAP_N][SKYMAP_N];
    double logNormalization;
} CAT(XLALSkymap2KernelType, SKYMAP_N);

// Stores pre-computed information for generating a sky map

typedef struct
{
    int sampleFrequency;
    LALDetector site[SKYMAP_N];
} CAT(XLALSkymap2PlanType, SKYMAP_N);

// Initialize sample frequency and HLV site information

// Use detector names from lal/packages/tools/include/LALDetectors.h :
//	LAL_TAMA_300_DETECTOR	=	0,
//	LAL_VIRGO_DETECTOR      =	1,
//	LAL_GEO_600_DETECTOR	=	2,
//	LAL_LHO_2K_DETECTOR     =	3,
//	LAL_LHO_4K_DETECTOR     =	4,
//	LAL_LLO_4K_DETECTOR     =	5,

void CAT(XLALSkymap2PlanConstruct, SKYMAP_N)(
    int sampleFrequency,
    int siteNumbers[SKYMAP_N],
    CAT(XLALSkymap2PlanType, SKYMAP_N)* plan
    );

// Turn directions into sample delays and antenna patterns

void CAT(XLALSkymap2DirectionPropertiesConstruct, SKYMAP_N)(
    CAT(XLALSkymap2PlanType, SKYMAP_N)* plan,
    XLALSkymap2SphericalPolarType* direction,
    CAT(XLALSkymap2DirectionPropertiesType, SKYMAP_N)* properties
    );

// Turn antenna patterns and waveform normalization into kernel matrix

void CAT(XLALSkymap2KernelConstruct, SKYMAP_N)(
    CAT(XLALSkymap2DirectionPropertiesType, SKYMAP_N)* properties,
    double wSw[SKYMAP_N],
    CAT(XLALSkymap2KernelType, SKYMAP_N)* kernel
    );

// Apply the kernel to some data

void CAT(XLALSkymap2Apply, SKYMAP_N)(
    CAT(XLALSkymap2DirectionPropertiesType, SKYMAP_N)* properties,
    CAT(XLALSkymap2KernelType, SKYMAP_N)* kernel,
    double* xSw[SKYMAP_N],
    int tau,
    double* logPosterior
    );


#ifdef __cplusplus
}
#endif

#undef CAT2
#undef CAT

#endif // SKYMAP_N

