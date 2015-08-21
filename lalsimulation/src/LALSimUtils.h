/*
 * Copyright (C) 2015 J. Creighton
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
 * @addtogroup LALSimUtils_h Header LALSimUtils.h
 * @ingroup lalsimulation_general
 * @brief Miscellaneous routines.
 */

#ifndef _LALSIMUTILS_H
#define _LALSIMUTILS_H

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

#include <lal/LALDatatypes.h>

/** @{ */

/**
 * @brief Ratio of horizon distance to sense-monitor range.
 *
 * This factor is used in XLALMeasureStandardSirenSenseMonitorRange().
 * 
 * `sensemon_range = horizon_dist / LAL_HORIZON_DISTANCE_OVER_SENSEMON_RANGE`
 *
 * The factor can be computed using Monte Carlo methods; its value has been
 * found to be 2.264778 +- 0.000002.  This constant keeps it to only 5
 * decimal places however.
 *
 * @sa Appendix D of
 * Bruce Allen, Warren G. Anderson, Patrick R. Brady, Duncan A. Brown, and
 * Jolien D. E. Creighton, "FINDCHIRP: An algorithm for detection of
 * gravitational waves from inspiraling compact binaries", Phys. Rev. D @b 85,
 * 122006 (2012) http://dx.doi.org/10.1103/PhysRevD.85.122006
 */
#define LAL_HORIZON_DISTANCE_OVER_SENSEMON_RANGE 2.26478

double XLALMeasureStandardSirenSenseMonitorRange(REAL8FrequencySeries *psd, double f_min, double f_max);
double XLALMeasureStandardSirenHorizonDistance(REAL8FrequencySeries *psd, double f_min, double f_max);
double XLALMeasureStandardSirenSNR(REAL8FrequencySeries *psd, double f_min, double f_max);
double XLALMeasureSNRFD(COMPLEX16FrequencySeries *htilde, REAL8FrequencySeries *psd, double f_min, double f_max);
double XLALMeasureSNR(REAL8TimeSeries *h, REAL8FrequencySeries *psd, double f_min, double f_max);

/** @} */

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif
