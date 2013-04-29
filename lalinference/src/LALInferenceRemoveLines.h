/*
 * LALInferenceRemoveLines.h    Utility functions for identifying lines 
 * in IFO data to be removed in LALInference
 *
 * Copyright (C) 2013 Michael Coughlin
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
*/

#ifndef LALInferenceRemoveLines_h
#define LALInferenceRemoveLines_h

/**
 * \defgroup LALInferenceRemoveLines_h Header LALInferenceRemoveLines.h
 * \ingroup pkg_LALInference
 * \brief Utility functions for identifying lines
 * in IFO data to be removed in LALInference
 */
/*@{*/

/** \brief Determine correlated frequency bands using cross correlation.
 * This function reads command line arguments and returns a REAL8 array. 
 * \author Michael Coughlin
 */

int LALInferenceXCorrBands(
    REAL8FrequencySeries        *spectrum,
    const REAL8TimeSeries       *tseries,
    UINT4                        seglen,
    UINT4                        stride,
    const REAL8Window           *window,
    const REAL8FFTPlan          *plan,
    REAL8                       *pvalues,
    char*                        filename
    );

/** \brief Determine non-Gaussian frequency bins using a chi-squared test.
 * This function reads command line arguments and returns a REAL8 array.
 * \author Michael Coughlin
*/

int LALInferenceRemoveLinesChiSquared(
    REAL8FrequencySeries        *spectrum,
    const REAL8TimeSeries       *tseries,
    UINT4                        seglen,
    UINT4                        stride,
    const REAL8Window           *window,
    const REAL8FFTPlan          *plan,
    REAL8                       *pvalues
    );

/** \brief Determine non-Gaussian frequency bins using a K-S test.
 * This function reads command line arguments and returns a REAL8 array.
 * \author Michael Coughlin
*/

int LALInferenceRemoveLinesKS(
    REAL8FrequencySeries        *spectrum,
    const REAL8TimeSeries       *tseries,
    UINT4                        seglen,
    UINT4                        stride,
    const REAL8Window           *window,
    const REAL8FFTPlan          *plan,
    REAL8                       *pvalues
    );

/** \brief Calculate PSD by fitting bins to lines.
 * This function reads command line arguments and returns a modified
 * REAL8FrequencySeries.
 * \author Michael Coughlin
*/

int LALInferenceAverageSpectrumBinFit(
    REAL8FrequencySeries        *spectrum,
    const REAL8TimeSeries       *tseries,
    UINT4                        seglen,
    UINT4                        stride,
    const REAL8Window           *window,
    const REAL8FFTPlan          *plan,
    char*                        filename,
    LIGOTimeGPS                 GPStime
    );

/** \brief Determine large amplitude frequency bins using power law fit.
 * This function reads command line arguments and returns a REAL8 array.
 * \author Michael Coughlin
 */

int LALInferenceRemoveLinesPowerLaw(
    REAL8FrequencySeries        *spectrum,
    const REAL8TimeSeries       *tseries,
    UINT4                        seglen,
    UINT4                        stride,
    const REAL8Window           *window,
    const REAL8FFTPlan          *plan,
    REAL8                       *pvalues
    );

/*@}*/

#endif
