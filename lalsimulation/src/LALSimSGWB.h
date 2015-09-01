/*
*  Copyright (C) 2011 Jolien Creighton
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

#ifndef _LALSIMSGWB_H
#define _LALSIMSGWB_H

#include <stddef.h>
#include <gsl/gsl_rng.h>
#include <lal/LALDatatypes.h>
#include <lal/LALDetectors.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * @defgroup LALSimSGWB_h Header LALSimSGWB.h
 * @ingroup lalsimulation_stochastic
 * @author Jolien Creighton
 * @brief Routines for simulating a stochastic gravitational-wave background.
 *
 * @{
 * @defgroup LALSimSGWB_c     Module LALSimSGWB.c
 * @defgroup LALSimSGWBORF_c  Module LALSimSGWBORF.c
 * @}
 */

/*
 * OVERLAP REDUCTION FUNCTION ROUTINE
 * in module LALSimSGWBORF.c
 */

double XLALSimSGWBOverlapReductionFunction(double f, const LALDetector *detector1, const LALDetector *detector2);


/*
 * ROUTINES TO GENERATE SGWB SPECTRA
 * in module LALSimSGWB.c
 */

REAL8FrequencySeries *XLALSimSGWBOmegaGWFlatSpectrum(double Omega0, double flow, double deltaF, size_t length);
REAL8FrequencySeries *XLALSimSGWBOmegaGWPowerLawSpectrum(double Omegaref, double alpha, double fref, double flow, double deltaF, size_t length);


/*
 * SGWB GENERATION ROUTINES
 * in module LALSimSGWB.c
 */

int XLALSimSGWB(REAL8TimeSeries **h, const LALDetector *detectors, size_t numDetectors, size_t stride, const REAL8FrequencySeries *OmegaGW, double H0, gsl_rng *rng);
int XLALSimSGWBFlatSpectrum(REAL8TimeSeries **h, const LALDetector *detectors, size_t numDetectors, size_t stride, double Omega0, double flow, double H0, gsl_rng *rng);
int XLALSimSGWBPowerLawSpectrum(REAL8TimeSeries **h, const LALDetector *detectors, size_t numDetectors, size_t stride, double Omegaref, double alpha, double fref, double flow, double H0, gsl_rng *rng);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMSGWB_H */
