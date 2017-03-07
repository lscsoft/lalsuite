/*
*  Copyright (C) 2013 Chris Messenger
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

#ifndef _SEMICOHERENT_H
#define _SEMICOHERENT_H

#include <lal/LALDatatypes.h>
#include <lal/Date.h>

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <sys/stat.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>        /* needed for the gsl interpolation */
#include <gsl/gsl_spline.h>        /* needed for the gsl interpolation */
#include <gsl/gsl_rng.h>           /* for random number generation */
#include <gsl/gsl_randist.h>       /* for random number generation */
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sf_log.h>        /* for log computation */
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/LALDatatypes.h>
#include <lal/Units.h>
#include <lal/SFTutils.h>
#include <lal/SFTfileIO.h>
#include <lal/ComplexFFT.h>
#include <lal/UserInput.h>
#include <lal/LogPrintf.h>
#include <lalapps.h>

/* includes */
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <lal/TimeSeries.h>
#include <lal/LALDatatypes.h>
#include <lal/Units.h>
#include <lal/SFTutils.h>
#include <lal/SFTfileIO.h>
#include <lal/ComplexFFT.h>
#include <lal/UserInput.h>
#include <lal/LogPrintf.h>
#include <lalapps.h>
#include <lal/BandPassTimeSeries.h>

#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/LALConstants.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/RealFFT.h>

#include <lal/TimeSeries.h>
#include <lal/GeneratePulsarSignal.h>
#include <lal/TransientCW_utils.h>
#include <lal/LALString.h>
#include <lal/StringVector.h>
#include <lal/NormalizeSFTRngMed.h>

#include <lal/SinCosLUT.h>

#include <gsl/gsl_rng.h>           /* for random number generation */
#include <gsl/gsl_randist.h>       /* for random number generation */


/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

  /***********************************************************************************************/
  /* some global constants */

  /* #define STRINGLENGTH 256              /\* the length of general string *\/ */
#define MINBAND 0.1                    /* the minimum amnount of bandwidth to read in */
/* #define SAMPFREQ 8192                    /\* the orginal sampling frequency *\/ */
#define LONGSTRINGLENGTH 1024         /* the length of general string */
#define NFREQMAX 4                    /* the max dimensionality of the frequency derivitive grid */
#define NBINMAX 4                        /* the number of binary parameter dimensions */
#define NBINS 4                       /* the number of bins to add to each side of the fft for safety */
#define WINGS_FACTOR 2                /* the safety factor in reading extra frequency from SFTs */
#define MAXNTSERIES 1073741824        /* the max number of samples in a timeseries at one time */

  /***********************************************************************************************/
  /* structure definitions */

  /** A single parameter dimensions boundaries
   */
  typedef struct {
    REAL8 min;                        /**< the parameter space minimum */
    REAL8 max;                        /**< the parameter space maximium */
    REAL8 span;                       /**< the parameter space span */
    CHAR name[LALNameLength];         /**< string containing the name of the dimension */
  } REAL8Dimension;

  /** A vector of parameter space boundary information
   */
  typedef struct {
    REAL8Dimension *data;             /**< the boundaries, span, etc for a single dimension */
    UINT4 ndim;                       /**< the number of dimensions */
  } REAL8Space;

  /** Stores the gridding parameters for a single dimension
   */
  typedef struct {
    REAL8 min;                        /**< the starting points of the grid */
    REAL8 delta;                      /**< the grid spacings */
    REAL8 oneoverdelta;               /**< the inverse of the spacing */
    UINT4 length;                     /**< the number of templates in each dimension */
    CHAR name[LALNameLength];         /**< string containing the name of the dimension */
  } Grid;

  /** Stores the current location in a hyper-cubic parameter space
   */
  typedef struct {
    REAL8 *x;                         /**< the location in parameter space */
    INT4 *idx;                        /**< the index of each dimension for this template */
    UINT4 ndim;                       /**< the dimension of the parameter space */
    UINT4 currentidx;                 /**< the current index value of the template */
  } Template;

  /** Stores the gridding parameters for a hypercubic grid of templates
   */
  typedef struct {
    Grid *grid;                       /**< stores the parameters defining a single dimension */
    UINT4 ndim;                       /**< the number of dimensions */
    UINT4 *prod;                      /**< internal variable used to store the size of sub-dimensions */
    UINT4 max;                        /**< the maximum (total) number of templates */
    REAL8 mismatch;                   /**< the mismatch */
    INT4 Nr;
  } GridParameters;

   /** Stores the gridding parameters for a hypercubic grid of templates
   */
  typedef struct {
    GridParameters **segment;             /**< stores the parameters defining a single dimension */
    UINT4 length;                        /**< the number of dimensions */
  } GridParametersVector;

  /** Stores the gridding parameters for a hypercubic grid of templates
   */
  typedef struct {
    COMPLEX8TimeSeries **data;             /**< stores the parameters defining a single dimension */
    UINT4 length;                        /**< the number of dimensions */
  } COMPLEX8TimeSeriesArray;

  /** Stores the gridding parameters for a hypercubic grid of templates
   */
  typedef struct {
    REAL4Vector *data;             /**< stores the parameters defining a single dimension */
    LIGOTimeGPS epoch;                        /**< the number of dimensions */
  } REAL4DemodulatedPower;

  /** Stores the gridding parameters for a hypercubic grid of templates
   */
  typedef struct {
    REAL4DemodulatedPower **segment;             /**< stores the parameters defining a single dimension */
    UINT4 length;                        /**< the number of dimensions */
  } REAL4DemodulatedPowerVector;

  /** Stores the gridding parameters for a hypercubic grid of templates
   */
  typedef struct {
    REAL8Space *space;             /**< stores the parameters defining a single dimension */
    REAL8 span;                        /**< the number of dimensions */
    LIGOTimeGPS epoch;
    REAL8 tseg;
  } ParameterSpace;

  /** Stores an array of REAL4 vectors
   */
  typedef struct {
    REAL4Vector **data;                  /**< stores REAL4 Vectors */
    UINT4 length;                        /**< the number of vectors */
  } REAL4VectorArray;

  /** Parameters for BinaryToSFT function
   */
  typedef struct {
    REAL8 tsamp;                      /**< the sampling time of the data */
    INT4 tsft;                        /**< the length of the SFTs */
    LIGOTimeGPS tstart;               /**< the first sft time stamp */
    REAL8 freq;                       /**< the starting frequency */
    REAL8 freqband;                   /**< the band width */
    REAL8 highpassf;                  /**< the high pass filter frequency */
    REAL8 amp_inj;                    /**< if set we inject a fake signal with this fractional amplitude */
    LIGOTimeGPS tref;
    REAL8 f_inj;
    REAL8 asini_inj;
    LIGOTimeGPS tasc_inj;
    REAL8 P_inj;
    REAL8 phi_inj;
    gsl_rng *r;
  } BinaryToSFTparams;

  int XLALReadSFTs(SFTVector**,CHAR *,REAL8,REAL8,INT4,INT4,INT4);
  int XLALComputeFreqGridParamsVector(GridParametersVector**,REAL8Space*,SFTVector*,REAL8);
  int XLALComputeFreqGridParams(GridParameters **freqgridparams,REAL8Space *pspace, REAL8 tmid,REAL8 tsft, REAL8 mu);
  int XLALSFTVectorToCOMPLEX8TimeSeriesArray(COMPLEX8TimeSeriesArray **dstimevec, SFTVector *sftvec);
  int XLALSFTToCOMPLEX8TimeSeries(COMPLEX8TimeSeries **ts, COMPLEX8FrequencySeries *sft,COMPLEX8FFTPlan **plan);
  int XLALCOMPLEX8TimeSeriesToCOMPLEX8FrequencySeries(COMPLEX8FrequencySeries **fs,const COMPLEX8TimeSeries *ts,GridParameters **gridparams);
  int XLALCOMPLEX8TimeSeriesArrayToDemodPowerVector(UINT8 *num_coarse,REAL4DemodulatedPowerVector**,COMPLEX8TimeSeriesArray*,GridParametersVector*,FILE *fp);
  int XLALApplyPhaseCorrection(COMPLEX8TimeSeries **outts, COMPLEX8TimeSeries *ints, Template *fn);
  /* int XLALNormaliseSFTs(SFTVector **sftvec,UINT4 blocksize); */
  int XLALComputeBinaryGridParams(GridParameters **binarygridparams,REAL8Space *space,REAL8 T,REAL8 DT,REAL8 mu,REAL8 coverage);
  int XLALGetNextTemplate(Template **,GridParameters *, ParameterSpace *,UNUSED void *);
  int XLALGetNextRandomBinaryTemplate(Template **,GridParameters *, ParameterSpace *,void *);
  int XLALComputeBinaryFreqDerivitives(Template *fdots,Template *bintemp,REAL8 tmid);
  int XLALFreeParameterSpace(ParameterSpace *pspace);
  int XLALFreeREAL4DemodulatedPowerVector(REAL4DemodulatedPowerVector *power);
  int XLALReplaceSFTVectornoise(SFTVector **sftvec,REAL8 background,INT4 seed);
  int XLALInitgslrand(gsl_rng **gslrnd,INT8 seed);
  int XLALCopySFT (SFTtype *dest, const SFTtype *src);
  int XLALAppendSFT2Vector (SFTVector *vect,const SFTtype *sft);
  int XLALBinaryToSFTVector(SFTVector **SFTvect,CHAR *filename,LIGOTimeGPS *fileStart,BinaryToSFTparams *par,INT8Vector **np, REAL8Vector **R);
  int XLALNormalizeSFTVectMean ( SFTVector *sftVect, REAL8Vector *means, INT4 flag);
  int XLALNormalizeSFTMean ( SFTtype  *sft, REAL8 *mean, INT4 flag);
  int XLALNormalizeSFTVectMedian ( SFTVector *sftVect, REAL8Vector *means, INT4 flag);
  int XLALNormalizeSFTMedian ( SFTtype  *sft, REAL8 *mean, INT4 flag);
#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
