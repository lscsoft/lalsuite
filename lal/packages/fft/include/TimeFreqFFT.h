/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Kipp Cannon, Patrick Brady, Tania Regimbau
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

#ifndef _TIMEFREQFFT_H
#define _TIMEFREQFFT_H

/* remove SWIG interface directives */
#if !defined(SWIG) && !defined(SWIGLAL_STRUCT)
#define SWIGLAL_STRUCT(...)
#endif

#include <lal/LALDatatypes.h>
#include <lal/ComplexFFT.h>
#include <lal/RealFFT.h>
#include <lal/Window.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

NRCSID( TIMEFREQFFTH, "$Id$" );

/**
 * \addtogroup TimeFreqFFT_h
 *
 * \brief Performs real-to-complex, complex-to-real FFTs and average power spectrum estimation.
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/TimeFreqFFT.h>
 * \endcode
 *
 * Perform time-to-frequency and frequency-to-time fast Fourier
 * transforms. Also provides a function to compute mean and median power
 * spectra with user specified windowning.
 *
 * The definition of the Fourier transform is \f$ e^{2 \pi i f t} \f$
 * and the inline equation version of that is
 * \f[
 *  \tilde{h}_k = \sum
 * \f]
 *
 */
/*@{*/

/**\name Error Codes */ /*@{*/
#define TIMEFREQFFTH_ENULL 1	/**< Null pointer */
#define TIMEFREQFFTH_ESIZE 2	/**< Invalid size */
#define TIMEFREQFFTH_ERATE 4	/**< Invalid rate */
#define TIMEFREQFFTH_ESIGN 4	/**< Incorrect plan sign */
#define TIMEFREQFFTH_EALLOC 16	/**< Pointer has already been allocated, should be null */
#define TIMEFREQFFTH_EPOSARG 32	/**< Argument must be positive */
#define TIMEFREQFFTH_EMALLOC 64	/**< Malloc failure */
#define TIMEFREQFFTH_EINCOMP 128/**< Incompatible arguments */
#define TIMEFREQFFTH_ENNUL 256	/**< Non-null pointer */
#define TIMEFREQFFTH_EZSEG 512	/**< Segment length is zero */
#define TIMEFREQFFTH_EZOVR 1024	/**< Overlap length is zero */
#define TIMEFREQFFTH_EMISM 2048	/**< Mismatch beteen segment, overlap and data length */
#define TIMEFREQFFTH_EUAVG 4096	/**< Unknown average power spectum method */
/*@}*/

/** \cond DONT_DOXYGEN */
#define TIMEFREQFFTH_MSGENULL "Null pointer"
#define TIMEFREQFFTH_MSGESIZE "Invalid size"
#define TIMEFREQFFTH_MSGERATE "Invalid rate"
#define TIMEFREQFFTH_MSGESIGN "Incorrect plan sign"
#define TIMEFREQFFTH_MSGEALLOC "Pointer has already been allocated, should be null"
#define TIMEFREQFFTH_MSGEPOSARG "Argument must be positive"
#define TIMEFREQFFTH_MSGEMALLOC "Malloc failure"
#define TIMEFREQFFTH_MSGEINCOMP "Incompatible arguments"
#define TIMEFREQFFTH_MSGENNUL "Non-null pointer"
#define TIMEFREQFFTH_MSGEZSEG "Segment length is zero"
#define TIMEFREQFFTH_MSGEZOVR "Overlap length is zero"
#define TIMEFREQFFTH_MSGEMISM "Mismatch beteen segment, overlap and data length"
#define TIMEFREQFFTH_MSGEUAVG "Unknown average power spectum method"
/** \endcond */


/** This type determines the method the type of average that will be used to
 * compute the power sperum estimate by the LALREAL4AverageSpectrum()
 * function. The function computes a series of (possibly overlapping) power
 * spectra and computes the average using one of the following methods.
 */
typedef enum
{
  useUnity,	/**< A constant PSD of value unity will be returned independent of the input data given; this is used for testing purposes */
  useMean,	/**< The arithmetic mean of the individual power spectra computed will be used to compute the output power spectrum */
  useMedian,	/**< The median value of the individual power spectra computed will be used to compute the output power spectrum */
  NumberAvgSpecMethods	/**< gives the number of defined methods */
}
AvgSpecMethod;

/** This structure controls the behaviour of the LALREAL4AverageSpectrum() function.
 */
typedef struct
tagAverageSpectrumParams
{
  SWIGLAL_STRUCT(AverageSpectrumParams);
  REAL4Window          *window;	/**< The windowing function to use when computing the individual power spectra from the input time series.
                                 * The input time series is broken into smaller time series to compute power spectra
                                 * for the estimate. The length of these time series is determined by the \c length parameter of the window vector.
                                 */
  UINT4                 overlap; /**< The overlap between sucessive time series used to compute the power spectra. */
  AvgSpecMethod         method;	/**< The method of computing the average as describe under ::AvgSpecMethod */
  RealFFTPlan          *plan;	/**< The FFT plan to be used in the computation of the power spectrum. */
}
AverageSpectrumParams;

/** UNDOCUMENTED */
typedef struct
tagLALPSDRegressor
{
  SWIGLAL_STRUCT(LALPSDRegressor);
  unsigned average_samples;
  unsigned median_samples;
  unsigned n_samples;
  REAL8Sequence **history;
  REAL8FrequencySeries *mean_square;
}
LALPSDRegressor;

/*
 *
 * XLAL Functions
 *
 */
int XLALREAL4TimeFreqFFT(
    COMPLEX8FrequencySeries *freq,
    const REAL4TimeSeries   *tser,
    const REAL4FFTPlan      *plan
    );

int XLALREAL4FreqTimeFFT(
    REAL4TimeSeries               *tser,
    const COMPLEX8FrequencySeries *freq,
    const REAL4FFTPlan            *plan
    );

int XLALREAL8TimeFreqFFT(
    COMPLEX16FrequencySeries *freq,
    const REAL8TimeSeries    *tser,
    const REAL8FFTPlan       *plan
    );

int XLALREAL8FreqTimeFFT(
    REAL8TimeSeries                *tser,
    const COMPLEX16FrequencySeries *freq,
    const REAL8FFTPlan             *plan
    );

int XLALCOMPLEX8TimeFreqFFT(
    COMPLEX8FrequencySeries  *freq,
    const COMPLEX8TimeSeries *tser,
    const COMPLEX8FFTPlan    *plan
    );

int XLALCOMPLEX8FreqTimeFFT(
    COMPLEX8TimeSeries            *tser,
    const COMPLEX8FrequencySeries *freq,
    const COMPLEX8FFTPlan         *plan
    );

int XLALCOMPLEX16TimeFreqFFT(
    COMPLEX16FrequencySeries  *freq,
    const COMPLEX16TimeSeries *tser,
    const COMPLEX16FFTPlan    *plan
    );

int XLALCOMPLEX16FreqTimeFFT(
    COMPLEX16TimeSeries            *tser,
    const COMPLEX16FrequencySeries *freq,
    const COMPLEX16FFTPlan         *plan
    );

int XLALREAL4ModifiedPeriodogram(
    REAL4FrequencySeries        *periodogram,
    const REAL4TimeSeries       *tseries,
    const REAL4Window           *window,
    const REAL4FFTPlan          *plan
    );

int XLALREAL8ModifiedPeriodogram(
    REAL8FrequencySeries        *periodogram,
    const REAL8TimeSeries       *tseries,
    const REAL8Window           *window,
    const REAL8FFTPlan          *plan
    );

int XLALREAL4AverageSpectrumWelch(
    REAL4FrequencySeries        *spectrum,
    const REAL4TimeSeries       *tseries,
    UINT4                        seglen,
    UINT4                        stride,
    const REAL4Window           *window,
    const REAL4FFTPlan          *plan
    );

int XLALREAL8AverageSpectrumWelch(
    REAL8FrequencySeries        *spectrum,
    const REAL8TimeSeries       *tseries,
    UINT4                        seglen,
    UINT4                        stride,
    const REAL8Window           *window,
    const REAL8FFTPlan          *plan
    );

REAL8 XLALMedianBias( UINT4 nn );

REAL8 XLALLogMedianBiasGeometric( UINT4 nn );

int XLALREAL4AverageSpectrumMedian(
    REAL4FrequencySeries        *spectrum,
    const REAL4TimeSeries       *tseries,
    UINT4                        seglen,
    UINT4                        stride,
    const REAL4Window           *window,
    const REAL4FFTPlan          *plan
    );

int XLALREAL8AverageSpectrumMedian(
    REAL8FrequencySeries        *spectrum,
    const REAL8TimeSeries       *tseries,
    UINT4                        seglen,
    UINT4                        stride,
    const REAL8Window           *window,
    const REAL8FFTPlan          *plan
    );

int XLALREAL4AverageSpectrumMedianMean(
    REAL4FrequencySeries        *spectrum,
    const REAL4TimeSeries       *tseries,
    UINT4                        seglen,
    UINT4                        stride,
    const REAL4Window           *window,
    const REAL4FFTPlan          *plan
    );

int XLALREAL8AverageSpectrumMedianMean(
    REAL8FrequencySeries        *spectrum,
    const REAL8TimeSeries       *tseries,
    UINT4                        seglen,
    UINT4                        stride,
    const REAL8Window           *window,
    const REAL8FFTPlan          *plan
    );

int XLALREAL4SpectrumInvertTruncate(
    REAL4FrequencySeries        *spectrum,
    REAL4                        lowfreq,
    UINT4                        seglen,
    UINT4                        trunclen,
    REAL4FFTPlan                *fwdplan,
    REAL4FFTPlan                *revplan
    );

int XLALREAL8SpectrumInvertTruncate(
    REAL8FrequencySeries        *spectrum,
    REAL8                        lowfreq,
    UINT4                        seglen,
    UINT4                        trunclen,
    REAL8FFTPlan                *fwdplan,
    REAL8FFTPlan                *revplan
    );

REAL4TimeSeries *XLALRespFilt(
    REAL4TimeSeries             *strain,
    COMPLEX8FrequencySeries     *transfer
    );

REAL4TimeSeries *XLALREAL4Convolution(
    REAL4TimeSeries             *strain,
    REAL4TimeSeries             *transfer
    );

COMPLEX8FrequencySeries *XLALWhitenCOMPLEX8FrequencySeries(
    COMPLEX8FrequencySeries     *fseries,
    const REAL4FrequencySeries  *psd
    );

COMPLEX16FrequencySeries *XLALWhitenCOMPLEX16FrequencySeries(
    COMPLEX16FrequencySeries    *fseries,
    const REAL8FrequencySeries  *psd
);

void
LALTimeFreqRealFFT(
    LALStatus               *status,
    COMPLEX8FrequencySeries *fser,
    REAL4TimeSeries         *tser,
    RealFFTPlan             *plan
    );

void
LALFreqTimeRealFFT(
    LALStatus               *status,
    REAL4TimeSeries         *tser,
    COMPLEX8FrequencySeries *fser,
    RealFFTPlan             *plan
    );

void
LALREAL4AverageSpectrum (
    LALStatus                   *status,
    REAL4FrequencySeries        *fSeries,
    REAL4TimeSeries             *tSeries,
    AverageSpectrumParams       *params
    );
void
LALCOMPLEX8AverageSpectrum (
    LALStatus                   *status,
    COMPLEX8FrequencySeries     *fSeries,
    REAL4TimeSeries             *tSeries0,
    REAL4TimeSeries             *tSeries1,
    AverageSpectrumParams       *params
    );
void
LALTimeFreqComplexFFT(
    LALStatus               *status,
    COMPLEX8FrequencySeries *fser,
    COMPLEX8TimeSeries      *tser,
    ComplexFFTPlan          *plan
    );

void
LALFreqTimeComplexFFT(
    LALStatus               *status,
    COMPLEX8TimeSeries      *tser,
    COMPLEX8FrequencySeries *fser,
    ComplexFFTPlan          *plan
    );

LALPSDRegressor *
XLALPSDRegressorNew(
    unsigned average_samples,
    unsigned median_samples
);

void
XLALPSDRegressorFree(
    LALPSDRegressor *r
);

void
XLALPSDRegressorReset(
    LALPSDRegressor *r
);

int XLALPSDRegressorSetAverageSamples(
    LALPSDRegressor *r,
    unsigned average_samples
);


unsigned XLALPSDRegressorGetAverageSamples(
    const LALPSDRegressor *r
);

int XLALPSDRegressorSetMedianSamples(
    LALPSDRegressor *r,
    unsigned median_samples
);

unsigned XLALPSDRegressorGetMedianSamples(
    const LALPSDRegressor *r
);

int
XLALPSDRegressorAdd(
    LALPSDRegressor *r,
    const COMPLEX16FrequencySeries *sample
);

REAL8FrequencySeries *
XLALPSDRegressorGetPSD(
    const LALPSDRegressor *r
);

int
XLALPSDRegressorSetPSD(
    LALPSDRegressor *r,
    const REAL8FrequencySeries *psd,
    unsigned weight
);


/*@}*/

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _TIMEFREQFFT_H */
