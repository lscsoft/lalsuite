/*
*  Copyright (C) 2007 Jolien Creighton, Kaice T. Reilly, Tania Regimbau, John Whelan
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
 * \addtogroup StochasticCrossCorrelation_h
 *
 * \author UTB Relativity Group; contact whelan@phys.utb.edu (original by S. Drasco)
 *
 * \brief Provides prototype and error code information for the \ref StochasticCrossCorrelation_c.
 *
 * Provides prototype and error code information for the modules needed
 * to calculate the standard optimally-filtered cross-correlation
 * statistic for stochastic background searches, given a pair of data
 * segments, along with appropriate representations of the detector
 * transfer function and the (uncalibrated) power spectral density of the
 * noise in each detector.  The relationship among these modules is
 * illustrated in Fig.\figref{stochastic_CrossCorrFlowchart}.
 *
 * \image html stochastic_CrossCorrFlowchart.png "Fig. [stochastic_CrossCorrFlowchart]: Relationship between the modules in StochasticCrossCorrelation.h"
 * \image latex stochastic_CrossCorrFlowchart.pdf "Relationship between the modules in StochasticCrossCorrelation.h"
 *
 * Figure \figref{stochastic_CrossCorrFlowchart} illustrates the relationship among the modules
 * dependent on \ref StochasticCrossCorrelation_h, which are used to calculate the cross-correlation
 * statistic \f$Y\f$ and its theoretical variance per unit time \f$\sigma^2/T\f$ from (uncalibrated)
 * stretches of data \f$h_1(t)\f$, \f$h_2(t)\f$, from two detectors, using metadata on the power
 * spectral densities \f$P_1(f)\f$, \f$P_2(f)\f$ and transfer functions \f$\tilde{R}_1(f)\f$,
 * \f$\tilde{R}_2(f)\f$ for each detector.
 *
 * \c CrossCorr represents the \ref StochasticCrossCorrelation_c (containing the functions
 * <tt>LALStochasticCrossCorrelationStatistic()</tt>,  <tt>LALStochasticHeterodynedCrossCorrelationStatistic()</tt>,
 * and <tt>LALStochasticCrossCorrelationSpectrum()</tt>),
 *
 * \c ZeroPadAndFFT represents the \ref ZeroPadAndFFT_c (containing the functions
 * <tt>LALSZeroPadAndFFT()</tt> and <tt>LALCZeroPadAndFFT()</tt>)
 *
 * \c OptimalFilter represents the \ref StochasticOptimalFilter_c (containing the function
 * <tt>LALStochasticOptimalFilter()</tt>)
 *
 * \c Normalization represents the \ref StochasticOptimalFilterNormalization_c
 * (containing the function <tt>LALStochasticOptimalFilterNormalization()</tt>)
 *
 * \c InverseNoise represents the \ref StochasticInverseNoise_c  (containing the function
 * <tt>LALStochasticInverseNoise()</tt>)
 * \c OmegaGW represents the \ref StochasticOmegaGW_c (containing the function <tt>LALStochasticOmegaGW()</tt>)
 *
 * \c Overlap represents the \ref OverlapReductionFunction_c (containing the function <tt>OverlapReductionFunction()</tt>)
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/StochasticCrossCorrelation.h>
 * \endcode
 *
 * @{
 */

#ifndef _STOCHASTICCROSSCORRELATION_H
#define _STOCHASTICCROSSCORRELATION_H

#include <lal/LALStdlib.h>
#include <lal/DetectorSite.h>
#include <lal/Units.h>
#include <lal/TimeFreqFFT.h>

#ifdef  __cplusplus
extern "C" {
#endif

/**\name Error Codes */ /*@{*/
#define STOCHASTICCROSSCORRELATIONH_ENULLPTR        1	/**< Null pointer */
#define STOCHASTICCROSSCORRELATIONH_ESAMEPTR        2	/**< Input and Output pointers the same */
#define STOCHASTICCROSSCORRELATIONH_EZEROLEN        3	/**< Zero length for data member of series */
#define STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF   4	/**< Negative or zero frequency spacing */
#define STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAT   5	/**< Negative or zero time spacing */
#define STOCHASTICCROSSCORRELATIONH_ENEGFMIN        6	/**< Negative start frequency */
#define STOCHASTICCROSSCORRELATIONH_EMMTIME         7	/**< Mismatch in epochs */
#define STOCHASTICCROSSCORRELATIONH_EMMHETERO       8	/**< Mismatch in heterodyning frequencies */
#define STOCHASTICCROSSCORRELATIONH_EMMFMIN         9	/**< Mismatch in start frequencies */
#define STOCHASTICCROSSCORRELATIONH_EMMDELTAF      10	/**< Mismatch in frequency spacings */
#define STOCHASTICCROSSCORRELATIONH_EMMLEN         11	/**< Mismatch in sequence lengths */
#define STOCHASTICCROSSCORRELATIONH_EOORFREF       12	/**< Out of range reference frequency */
#define STOCHASTICCROSSCORRELATIONH_ENONPOSOMEGA   13	/**< Negative stochastic background strength */
#define STOCHASTICCROSSCORRELATIONH_ENONSYMDIJ     14	/**< Non-symmetric response tensor */
#define STOCHASTICCROSSCORRELATIONH_ENONZEROHETERO 15	/**< Non-zero heterodyning frequency specified for real time series */
#define STOCHASTICCROSSCORRELATIONH_EWRONGUNITS    16	/**< Inconsistent input units */
#define STOCHASTICCROSSCORRELATIONH_ENONPOSWIN     17	/**< Zero or negative total for window functions */
#define STOCHASTICCROSSCORRELATIONH_EMEMORY        18	/**< Memory error */
#define STOCHASTICCROSSCORRELATIONH_ENOTYETHETERO 255	/**< Non-zero heterodyning frequency not yet implemented */
/*@}*/
/** \cond DONT_DOXYGEN */
#define STOCHASTICCROSSCORRELATIONH_MSGENULLPTR    "Null pointer"
#define STOCHASTICCROSSCORRELATIONH_MSGESAMEPTR    "Input and Output pointers the same"
#define STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN    "Zero length for data member of series"
#define STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF "Negative or zero frequency spacing"
#define STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAT "Negative or zero time spacing"
#define STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN "Negative start frequency"
#define STOCHASTICCROSSCORRELATIONH_MSGEMMTIME     "Mismatch in epochs"
#define STOCHASTICCROSSCORRELATIONH_MSGEMMHETERO   "Mismatch in heterodyning frequencies"
#define STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN     "Mismatch in start frequencies"
#define STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF   "Mismatch in frequency spacings"
#define STOCHASTICCROSSCORRELATIONH_MSGEMMLEN      "Mismatch in sequence lengths"
#define STOCHASTICCROSSCORRELATIONH_MSGEOORFREF    "Out of range reference frequency"
#define STOCHASTICCROSSCORRELATIONH_MSGENONPOSOMEGA "Negative stochastic background strength"
#define STOCHASTICCROSSCORRELATIONH_MSGENONSYMDIJ   "Non-symmetric response tensor"
#define STOCHASTICCROSSCORRELATIONH_MSGENONZEROHETERO "Non-zero heterodyning frequency specified for real time series"
#define STOCHASTICCROSSCORRELATIONH_MSGEWRONGUNITS "Inconsistent input units"
#define STOCHASTICCROSSCORRELATIONH_MSGENONPOSWIN  "Zero or negative total for window functions"
#define STOCHASTICCROSSCORRELATIONH_MSGEMEMORY     "Memory error"
#define STOCHASTICCROSSCORRELATIONH_MSGENOTYETHETERO   "Non-zero heterodyning frequency not yet implemented"
/** \endcond */

  /*************************************************************
   *                                                           *
   *       Structures and prototypes associated with           *
   *             StochasticCrossCorrelation.c                  *
   *                                                           *
   *************************************************************/

/**
 * Represents a dimensionful number as a 4-byte float with an associated
 * units structure, which is the output of
 * <tt>LALStochasticCrossCorrelationStatistic()</tt>.  The fields are:
 *
 * <dl>
 * <dt><tt>REAL4 value</tt></dt><dd>
 * The numerical value.</dd>
 *
 * <dt><tt>LALUnit units</tt></dt><dd>
 * The units.</dd>
 * </dl>
 */
typedef struct tagREAL4WithUnits {
  REAL4     value;
  LALUnit   units;
} REAL4WithUnits;

/**
 * Represents a dimensionful number as a single-precision (8-byte) complex
 * number with an associated
 * units structure, which is the output of
 * <tt>LALStochasticHeterodynedCrossCorrelationStatistic()</tt>.  The fields are:
 *
 * <dl>
 * <dt><tt>COMPLEX8 value</tt></dt><dd>
 * The numerical value.</dd>
 *
 * <dt><tt>LALUnit units</tt></dt><dd>
 * The units.</dd>
 * </dl>
 */
typedef struct tagCOMPLEX8WithUnits {
  COMPLEX8  value;
  LALUnit   units;
} COMPLEX8WithUnits;

/**
 * Contains the input data needed by
 * <tt>LALStochasticCrossCorrelationStatistic()</tt>
 * to calculate the value of the standard optimally-filtered
 * cross-correlation statistic.  The fields are:
 *
 * <dl>
 * <dt><tt>COMPLEX8FrequencySeries  *hBarTildeOne</tt></dt><dd>
 * Fourier transform of the first zero-padded data stream.</dd>
 *
 * <dt><tt>COMPLEX8FrequencySeries  *hBarTildeTwo</tt></dt><dd>
 * Fourier transform of the second zero-padded data stream.</dd>
 *
 * <dt><tt>COMPLEX8FrequencySeries  *optimalFilter</tt></dt><dd>
 * Optimal filter function in the frequency domain.</dd>
 * </dl>
 */
typedef struct tagStochasticCrossCorrelationInput {
  COMPLEX8FrequencySeries  *hBarTildeOne;
  COMPLEX8FrequencySeries  *hBarTildeTwo;
  COMPLEX8FrequencySeries  *optimalFilter;
} StochasticCrossCorrelationInput;

typedef struct tagStochasticCrossCorrelationStrainInput {
  COMPLEX8FrequencySeries  *hBarTildeOne;
  COMPLEX8FrequencySeries  *hBarTildeTwo;
  REAL4FrequencySeries  *optimalFilter;
} StochasticCrossCorrelationStrainInput;


typedef struct tagStochasticCrossCorrelationCalInput {
  COMPLEX8FrequencySeries  *hBarTildeOne;
  COMPLEX8FrequencySeries  *hBarTildeTwo;
  REAL4FrequencySeries     *optimalFilter;
  COMPLEX8FrequencySeries  *responseFunctionOne;
  COMPLEX8FrequencySeries  *responseFunctionTwo;
} StochasticCrossCorrelationCalInput;




void
LALStochasticCrossCorrelationStatistic(
            LALStatus                              *status,
            REAL4WithUnits                         *output,
            const StochasticCrossCorrelationInput  *input,
            BOOLEAN                                 epochsMatch);

void
LALStochasticHeterodynedCrossCorrelationStatistic(
            LALStatus                              *status,
            COMPLEX8WithUnits                      *output,
            const StochasticCrossCorrelationInput  *input,
            BOOLEAN                                 epochsMatch);

void
LALStochasticCrossCorrelationSpectrum(
            LALStatus                              *status,
            COMPLEX8FrequencySeries                *output,
            const StochasticCrossCorrelationInput  *input,
            BOOLEAN                                 epochsMatch);

void
LALStochasticCrossCorrelationStatisticStrain(
            LALStatus                              *status,
            REAL4WithUnits                         *output,
            const StochasticCrossCorrelationStrainInput  *input,
            BOOLEAN                                 epochsMatch);

void
LALStochasticHeterodynedCrossCorrelationStatisticStrain(
            LALStatus                              *status,
            COMPLEX8WithUnits                      *output,
            const StochasticCrossCorrelationStrainInput  *input,
            BOOLEAN                                 epochsMatch);

void
LALStochasticCrossCorrelationSpectrumStrain(
            LALStatus                              *status,
            COMPLEX8FrequencySeries                *output,
            const StochasticCrossCorrelationStrainInput  *input,
            BOOLEAN                                 epochsMatch);

void
LALStochasticCrossCorrelationStatisticCal(
            LALStatus                                 *status,
            REAL4WithUnits                            *output,
            const StochasticCrossCorrelationCalInput  *input,
            BOOLEAN                                    epochsMatch);

void
LALStochasticHeterodynedCrossCorrelationStatisticCal(
            LALStatus                                 *status,
            COMPLEX8WithUnits                         *output,
            const StochasticCrossCorrelationCalInput     *input,
            BOOLEAN                                    epochsMatch);

void
LALStochasticCrossCorrelationSpectrumCal(
            LALStatus                                 *status,
            COMPLEX8FrequencySeries                   *output,
            const StochasticCrossCorrelationCalInput  *input,
            BOOLEAN                                   epochsMatch);



  /*************************************************************
   *                                                           *
   * Structures and prototypes associated with ZeroPadAndFFT.c *
   *                                                           *
   *************************************************************/

/**
 * Contains the parameters of <tt>LALSZeroPadAndFFT()</tt>.
 * The fields are:
 *
 * <dl>
 * <dt><tt>RealFFTPlan *fftPlan</tt></dt><dd>
 * The FFT plan to be used by FFTW</dd>
 * <dt><tt>REAL4Vector *window</tt></dt><dd>
 * The window which is to be applied to the data</dd>
 * <dt><tt>UINT4 length</tt></dt><dd>
 * The length of the data after zero-padding</dd>
 * </dl>
 */
typedef struct tagSZeroPadAndFFTParameters {
  RealFFTPlan           *fftPlan;
  REAL4Window           *window;
  UINT4                  length;
} SZeroPadAndFFTParameters;

/**
 * Contains the parameters of <tt>LALCZeroPadAndFFT()</tt>.
 * The fields are:
 *
 * <dl>
 * <dt><tt>ComplexFFTPlan *fftPlan</tt></dt><dd>
 * The FFT plan to be used by FFTW</dd>
 * <dt><tt>REAL4Vector *window</tt></dt><dd>
 * The window which is to be applied to the data</dd>
 * <dt><tt>UINT4 length</tt></dt><dd>
 * The length of the data after zero-padding</dd>
 * </dl>
 *
 */
typedef struct tagCZeroPadAndFFTParameters {
  ComplexFFTPlan        *fftPlan;
  REAL4Window           *window;
  UINT4                  length;
} CZeroPadAndFFTParameters;

void
LALSZeroPadAndFFT(LALStatus                *status,
                  COMPLEX8FrequencySeries  *output,
                  const REAL4TimeSeries    *input,
                  SZeroPadAndFFTParameters *parameters);

void
LALCZeroPadAndFFT(LALStatus                *status,
                  COMPLEX8FrequencySeries  *output,
                  const COMPLEX8TimeSeries *input,
                  CZeroPadAndFFTParameters *parameters);



  /*************************************************************
   *                                                           *
   *   Structures and prototypes associated with               *
   *               StochasticOptimalFilter.c                   *
   *                                                           *
   *************************************************************/

/**
 * Contains the inputs of <tt>LALStochasticOptimalFilter()</tt>.
 * The fields are:
 *
 * <dl>
 * <dt><tt>REAL4FrequencySeries *overlapReductionFunction</tt></dt><dd>
 * The overlap reduction function \f$\gamma(f)\f$ describing the pair of detector
 * sites.</dd>
 * <dt><tt>REAL4FrequencySeries *omegaGW</tt></dt><dd> The spectrum
 * \f$\Omega_{\mathrm{GW}}(f)\f$ of the stochastic gravitational-wave
 * background.</dd>
 * <dt><tt>COMPLEX8FrequencySeries *halfCalibratedInverseNoisePSD1</tt></dt><dd>
 * The reciprocal
 * \f$1/P_1^{\mathrm{HC}}(f)
 * =1/(\tilde{R_1}(f)P_1(f))
 * =\tilde{R_1}(f)^* / P_1(f)\f$ of the
 * half-calibrated noise power spectral density for the first detector.</dd>
 * <dt><tt>COMPLEX8FrequencySeries *halfCalibratedInverseNoisePSD2</tt></dt><dd>
 * The reciprocal
 * \f$1/P_2^{\mathrm{HC}}(f)
 * =1/(\tilde{R_2}(f)P_2(f))
 * =\tilde{R_2}(f)^* / P_2(f)\f$ of the
 * half-calibrated noise power spectral density for the second detector.</dd>
 * </dl>
 *
 */
typedef struct tagStochasticOptimalFilterInput {
  REAL4FrequencySeries     *overlapReductionFunction;
  REAL4FrequencySeries     *omegaGW;
  COMPLEX8FrequencySeries  *halfCalibratedInverseNoisePSD1;
  COMPLEX8FrequencySeries  *halfCalibratedInverseNoisePSD2;
} StochasticOptimalFilterInput;

typedef struct tagStochasticOptimalFilterCalInput {
  REAL4FrequencySeries     *overlapReductionFunction;
  REAL4FrequencySeries     *omegaGW;
  REAL4FrequencySeries     *calibratedInverseNoisePSD1;
  REAL4FrequencySeries     *calibratedInverseNoisePSD2;
} StochasticOptimalFilterCalInput;



void
LALStochasticOptimalFilter(
            LALStatus                                *status,
            COMPLEX8FrequencySeries                  *optimalFilter,
            const StochasticOptimalFilterInput       *input,
            const REAL4WithUnits                     *lambda);

void
LALStochasticOptimalFilterCal(
            LALStatus                                *status,
            REAL4FrequencySeries                     *optimalFilter,
            const StochasticOptimalFilterCalInput    *input,
            const REAL4WithUnits                     *lambda);



  /*************************************************************
   *                                                           *
   *        Structures and prototypes associated with          *
   *         StochasticOptimalFilterNormalization.c            *
   *                                                           *
   *************************************************************/

/**
 * Contains the outputs of <tt>LALStochasticOptimalFilterNormalization()</tt>.
 * The fields are:
 *
 * <dl>
 * <dt><tt>REAL4WithUnits *normalization</tt></dt><dd>
 * The normalization parameter \f$\lambda\f$.</dd>
 * <dt><tt>REAL4WithUnits *variance</tt></dt><dd>
 * The variance per unit time \f$\sigma^2/T\f$ of the cross-correlation statistic.</dd>
 * </dl>
 *
 */
typedef struct tagStochasticOptimalFilterNormalizationOutput {
  REAL4WithUnits           *normalization;
  REAL4WithUnits           *variance;
} StochasticOptimalFilterNormalizationOutput;

/**
 * Contains the inputs of <tt>LALStochasticOptimalFilterNormalization()</tt>.
 * The fields are:
 *
 * <dl>
 * <dt><tt>REAL4FrequencySeries *overlapReductionFunction</tt></dt><dd>
 * The overlap reduction function \f$\gamma(f)\f$ describing the pair of detector
 * sites.</dd>
 * <dt><tt>REAL4FrequencySeries *omegaGW</tt></dt><dd> The spectrum
 * \f$\Omega_{\mathrm{GW}}(f)\f$ of the stochastic gravitational-wave
 * background.</dd>
 * <dt><tt>REAL4FrequencySeries *inverseNoisePSD1</tt></dt><dd>
 * The reciprocal
 * \f$1/P_1(f)=|\tilde{R_1}(f)|^2/P_1(f)\f$ of the
 * ununcalibrated noise power spectral density for the first detector.</dd>
 * <dt><tt>REAL4FrequencySeries *inverseNoisePSD2</tt></dt><dd>
 * The reciprocal
 * \f$1/P_2(f)=|\tilde{R_2}(f)|^2/P_2(f)\f$ of the
 * ununcalibrated noise power spectral density for the second detector.</dd>
 * </dl>
 */
typedef struct tagStochasticOptimalFilterNormalizationInput {
  REAL4FrequencySeries     *overlapReductionFunction;
  REAL4FrequencySeries     *omegaGW;
  REAL4FrequencySeries     *inverseNoisePSD1;
  REAL4FrequencySeries     *inverseNoisePSD2;
} StochasticOptimalFilterNormalizationInput;

/**
 * Contains the parameters of <tt>LALStochasticOptimalFilterNormalization()</tt>.
 * The fields are:
 *
 * <dl>
 * <dt><tt>REAL8 fRef</tt></dt><dd>
 * The reference frequency used in defining the normalization.</dd>
 * <dt><tt>BOOLEAN heterodyned</tt></dt><dd>
 * Indicates whether the filter is to be used on heterodyned data or not.</dd>
 * <dt><tt>REAL4Vector window1</tt></dt><dd>
 * The windowing function with which the first data stream was windowed</dd>
 * <dt><tt>REAL4Vector window2</tt></dt><dd>
 * The windowing function with which the second data stream was windowed</dd>
 * </dl>
 */
typedef struct tagStochasticOptimalFilterNormalizationParameters {
  REAL8               fRef;
  BOOLEAN             heterodyned;
  REAL4Vector        *window1;
  REAL4Vector        *window2;
} StochasticOptimalFilterNormalizationParameters;



void
LALStochasticOptimalFilterNormalization(
            LALStatus                                            *status,
            StochasticOptimalFilterNormalizationOutput           *output,
            const StochasticOptimalFilterNormalizationInput      *input,
            const StochasticOptimalFilterNormalizationParameters *parameters);



  /*************************************************************
   *                                                           *
   * Structures and prototypes associated with StochasticInverseNoise.c  *
   *                                                           *
   *************************************************************/

/**
 * Contains the outputs of <tt>LALStochasticInverseNoise()</tt>.
 * The fields are:
 *
 * <dl>
 * <dt><tt>REAL4FrequencySeries *calibratedInverseNoisePSD</tt></dt><dd>
 * The reciprocal
 * \f$1/P^{\mathrm{C}}(f)=|\tilde{R}(f)|^2/P(f)\f$ of the
 * ununcalibrated noise power spectral density.</dd>
 *
 * <dt><tt>COMPLEX8FrequencySeries *halfCalibratedInverseNoisePSD</tt></dt><dd>
 * The reciprocal \\
 * \f$1/P^{\mathrm{HC}}(f)=\tilde{R}(f)^* / P(f)\f$
 * of the half-calibrated noise power spectral density.</dd>
 * </dl>
 */
typedef struct tagStochasticInverseNoiseOutput {
  REAL4FrequencySeries     *calibratedInverseNoisePSD;
  COMPLEX8FrequencySeries  *halfCalibratedInverseNoisePSD;
} StochasticInverseNoiseOutput;

typedef struct tagStochasticInverseNoiseCalOutput {
  REAL4FrequencySeries     *calibratedInverseNoisePSD;
} StochasticInverseNoiseCalOutput;

/**
 * Contains the inputs to <tt>LALStochasticInverseNoise()</tt>.
 * The fields are:
 *
 * <dl>
 * <dt><tt>REAL4FrequencySeries *unCalibratedNoisePSD</tt></dt><dd>
 * The power spectral density \f$P(f)\f$ of the noise
 * contribution to the detector output.</dd>
 *
 * <dt><tt>COMPLEX8FrequencySeries *responseFunction</tt></dt><dd>
 * The frequency-domain reponse function \f$\tilde{R}(f)\f$.</dd>
 * </dl>
 */
typedef struct tagStochasticInverseNoiseInput {
  REAL4FrequencySeries     *unCalibratedNoisePSD ;
  COMPLEX8FrequencySeries  *responseFunction;
} StochasticInverseNoiseInput;



void
LALStochasticInverseNoise(
            LALStatus                             *status,
            StochasticInverseNoiseOutput          *output,
            const StochasticInverseNoiseInput     *input);

void
LALStochasticInverseNoiseCal(
            LALStatus                             *status,
            StochasticInverseNoiseCalOutput       *output,
            const StochasticInverseNoiseInput     *input);



  /*************************************************************
   *                                                           *
   *    Structures and prototypes associated with StochasticOmegaGW.c    *
   *                                                           *
   *************************************************************/

/**
 * Contains the parameters used by <tt>LALStochasticOmegaGW()</tt> to define a
 * power law: \f$\Omega_{\mathrm{GW}}(f)
 * = \Omega_{\mathrm{R}}(f/f_{\mathrm{R}})^\alpha\f$.
 * The fields are:
 *
 * <dl>
 * <dt><tt>REAL4 alpha</tt></dt><dd> The power-law exponent.</dd>
 *
 * <dt><tt>REAL8 fRef</tt></dt><dd> The reference frequency \f$f_{\mathrm{R}}\f$ used to define the normalization.</dd>
 *
 * <dt><tt>REAL4 omegaRef</tt></dt><dd> The amplitude
 * \f$\Omega_{\mathrm{R}}</tt>
 * =\Omega_{\mathrm{GW}}(f_{\mathrm{R}})\f$
 * at reference frequency.</dd>
 *
 * <dt><tt>UINT4 length</tt></dt><dd>
 * The number of points in the output frequency series.</dd>
 *
 * <dt><tt>REAL8 f0</tt></dt><dd>
 * The start frequency of the output frequency series.</dd>
 *
 * <dt><tt>REAL8 deltaF</tt></dt><dd>
 * The frequency spacing of the output frequency series.</dd>
 * </dl>
 */
typedef struct tagStochasticOmegaGWParameters {
  REAL4     alpha;    /**< exponent in power law: omegaGW(f) = f^alpha */
  UINT4     length;   /**< length of vector containing omegaGW(f) values */
  REAL8     f0;       /**< start frequency */
  REAL8     deltaF;   /**< frequency spacing */
  REAL8     fRef;     /**< reference normalization frequency */
  REAL4     omegaRef; /**< refenence omega coefficent for normalization */
}
StochasticOmegaGWParameters;



void
LALStochasticOmegaGW (
            LALStatus                          *status,
            REAL4FrequencySeries               *output,
            const StochasticOmegaGWParameters  *parameters);



  /*************************************************************
   *                                                           *
   *      Structures and prototypes associated with            *
   *            OverlapReductionFunction.c                     *
   *                                                           *
   *************************************************************/

/**
 * Contains the parameters used by
 * <tt>LALOverlapReductionFunction()</tt> to determine the format of its
 * output for the overlap reduction function.  The fields are:
 *
 * <dl>
 * <dt><tt>UINT4 length</tt></dt><dd>
 * The number of points in the output frequency series.</dd>
 *
 * <dt><tt>REAL8 f0</tt></dt><dd>
 * The start frequency of the output frequency series.</dd>
 *
 * <dt><tt>REAL8 deltaF</tt></dt><dd>
 * The frequency spacing of the output frequency series.</dd>
 * </dl>
 */
typedef struct tagOverlapReductionFunctionParameters {
  UINT4     length;   /**< length of vector containing overlap red function */
  REAL8     f0;       /**< start frequency */
  REAL8     deltaF;   /**< frequency spacing for overlap reduction function */
}
OverlapReductionFunctionParameters;

/**
 * Holds structures defining the location and orientation of a
 * pair of gravitational wave detectors.  This is the input to
 * <tt>LALOverlapReductionFunction()</tt>.  The fields are:
 *
 * <dl>
 * <dt><tt>LALDetector detectorOne</tt></dt><dd>
 * The first interferometer.</dd>
 *
 * <dt><tt>LALDetector detectorTwo</tt></dt><dd>
 * The second interferometer.</dd>
 * </dl>
 */
typedef struct tagLALDetectorPair {
  LALDetector    detectorOne;
  LALDetector    detectorTwo;
}
LALDetectorPair;



void
LALOverlapReductionFunction(
                   LALStatus                                  *status,
                   REAL4FrequencySeries                       *output,
                   const LALDetectorPair                      *detectors,
                   const OverlapReductionFunctionParameters   *parameters);

/** @} */

#ifdef  __cplusplus
}
#endif /* C++ protection */

#endif /* _STOCHASTICCROSSCORRELATION_H */
