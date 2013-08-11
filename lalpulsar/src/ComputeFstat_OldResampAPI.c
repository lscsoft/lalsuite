/*
 * Copyright (C) 2009 Chris Messenger, Reinhard Prix, Pinkesh Patel, Xavier Siemens, Holger Pletsch
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
 * \author Chris Messenger, Reinhard Prix, Pinkesh Patel, Xavier Siemens, Holger Pletsch
 * \ingroup pulsarCoherent
 * \file
 * \brief
 * Functions to calculate the so-called multi-detector F-statistic for a given
 * frequency band in parameter space using a resampling technique.
 *
 */

/*---------- INCLUDES ----------*/
#define __USE_ISOC99 1
#include <math.h>

#include <lal/ExtrapolatePulsarSpins.h>
#include <lal/FindRoot.h>
#include <lal/CWFastMath.h>

/* GSL includes */
#include <lal/LALGSL.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>
#include <lal/LFTandTSutils.h>
#include <lal/ComplexAM.h>
#include <lal/TimeSeries.h>

#include "ComputeFstat_OldResampAPI.h"

/*---------- local DEFINES ----------*/
#define TRUE (1==1)
#define FALSE (1==0)


#define LD_SMALL4       (2.0e-4)		/**< "small" number for REAL4*/
#define OOTWOPI         (1.0 / LAL_TWOPI)	/**< 1/2pi */

#define TWOPI_FLOAT     6.28318530717958f  	/**< single-precision 2*pi */
#define OOTWOPI_FLOAT   (1.0f / TWOPI_FLOAT)	/**< single-precision 1 / (2pi) */

#define EA_ACC          1E-9                    /* the timing accuracy of LALGetBinaryTimes in seconds */

/*----- Macros ----- */
/** convert GPS-time to REAL8 */
#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

#define MYSIGN(x) ( ((x) < 0) ? (-1.0):(+1.0) )
#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )
#define SQ(x) ((x)*(x))

/** Simple Euklidean scalar product for two 3-dim vectors in cartesian coords */
#define SCALAR(u,v) ((u)[0]*(v)[0] + (u)[1]*(v)[1] + (u)[2]*(v)[2])

#define NhalfPosDC(N) ((UINT4)(ceil ( ((N)/2.0 - 1e-6 ))))	/* round up */
#define NhalfNeg(N) ((UINT4)( (N) - NhalfPosDC(N) ))		/* round down (making sure N+ + N- = (N-1) */

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- Global variables ----------*/
#define NUM_FACT 7
static const REAL8 inv_fact[NUM_FACT] = { 1.0, 1.0, (1.0/2.0), (1.0/6.0), (1.0/24.0), (1.0/120.0), (1.0/720.0) };

/* empty initializers  */
static const LALStatus empty_status;
static const AMCoeffs empty_AMCoeffs;

/*---------- internal prototypes ----------*/
int finite(double x);

/*==================== FUNCTION DEFINITIONS ====================*/


/**
 * Function to compute a vector of Fstatistic values for a number of frequency bins.
 * The output, i.e. fstatVector must be properly allocated
 * before this function is called.  The values of the start frequency, the step size
 * in the frequency and the number of frequency values for which the Fstatistic is
 * to be calculated are read from fstatVector.  The other parameters are not checked and
 * they must be correctly set outside this function.
 */
void ComputeFStatFreqBand_RS ( LALStatus *status,				/**< pointer to LALStatus structure */
			       REAL4FrequencySeries *fstatVector, 		/**< [out] Vector of Fstat values */
			       const PulsarDopplerParams *doppler,		/**< parameter-space point to compute F for */
			       MultiSFTVector *multiSFTs, 		        /**< normalized (by DOUBLE-sided Sn!) data-SFTs of all IFOs */
			       const MultiNoiseWeights *multiWeights,	        /**< noise-weights of all SFTs */
			       ComputeFParams *params		                /**< addition computational params */
			       )
{
  UINT4 numDetectors;
  ComputeFBuffer_RS *cfBuffer = NULL;
  MultiDetectorStateSeries *multiDetStates = NULL;
  MultiSSBtimes *multiSSB = NULL;
  MultiAMCoeffs *multiAMcoef = NULL;
  MultiCOMPLEX8TimeSeries *multiTimeseries = NULL;
  REAL8 Ad, Bd, Cd, Dd_inv, AdX, BdX, CdX, DdX_inv;
  SkyPosition skypos;
  UINT4 i,j,k;
  MultiCOMPLEX8TimeSeries *multiFa_resampled = NULL;
  MultiCOMPLEX8TimeSeries *multiFb_resampled = NULL;
  COMPLEX8Vector *Faf_resampled = NULL;
  COMPLEX8Vector *Fbf_resampled = NULL;
  UINT4 numSamples, kmax;
  REAL8 f0_shifted;
  REAL8 f0_shifted_single;
  REAL8 dt;
  REAL8 df_out;
  ComplexFFTPlan *pfwd = NULL;  /* this will store the FFT plan */
  COMPLEX8Vector *outa = NULL;  /* this will contain the FFT output of Fa for this detector */
  COMPLEX8Vector *outb = NULL;  /* this will contain the FFT output of Fb for this detector */
  COMPLEX8Vector *outaSingle = NULL; /* this will contain Faf_resampled for a single IFO */
  COMPLEX8Vector *outbSingle = NULL; /* this will contain Fbf_resampled for a single IFO */

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* check that the input data and parameters structures don't point to NULL */
  ASSERT ( multiSFTs, status, COMPUTEFSTATRSC_ENULL, COMPUTEFSTATRSC_MSGENULL );
  ASSERT ( doppler, status, COMPUTEFSTATRSC_ENULL, COMPUTEFSTATRSC_MSGENULL );
  ASSERT ( params, status, COMPUTEFSTATRSC_ENULL, COMPUTEFSTATRSC_MSGENULL );

  cfBuffer = params->buffer;                      /* set local pointer to the buffer location */
  numDetectors = multiSFTs->length;               /* set the number of detectors to the number of sets of SFTs */
  // unused: SFTtype * firstSFT = &(multiSFTs->data[0]->data[0]);      /* use data from the first SFT from the first detector to set other params */

  /* check that the pre-allocated output vector doesn't point to NULL */
  ASSERT ( fstatVector, status, COMPUTEFSTATRSC_ENULL, COMPUTEFSTATRSC_MSGENULL );
  ASSERT ( fstatVector->data, status, COMPUTEFSTATRSC_ENULL, COMPUTEFSTATRSC_MSGENULL );
  ASSERT ( fstatVector->data->data, status, COMPUTEFSTATRSC_ENULL, COMPUTEFSTATRSC_MSGENULL );
  ASSERT ( fstatVector->data->length > 0, status, COMPUTEFSTATRSC_EINPUT, COMPUTEFSTATRSC_MSGEINPUT );

  df_out = fstatVector->deltaF;                   /* the user defined frequency resolution */

  /* check that the multidetector noise weights have the same length as the multiSFTs */
  if ( multiWeights ) {
    ASSERT ( multiWeights->length == numDetectors , status, COMPUTEFSTATRSC_EINPUT, COMPUTEFSTATRSC_MSGEINPUT );
  }

  /* first, if there is no buffer allocate space for one */
  /* IMPORTANT - use Calloc here so that all pointers within the structure are NULL */
  if ( cfBuffer == NULL ) {
    if ( (cfBuffer = (ComputeFBuffer_RS*)XLALCalloc(1,sizeof(ComputeFBuffer_RS))) == NULL ) {
      XLALPrintError("\nXLALMalloc() failed with error = %d\n\n", xlalErrno );
      ABORT ( status, COMPUTEFSTATRSC_EXLAL, COMPUTEFSTATRSC_MSGEXLAL );
    }
  }

  /************************************************************************************************************/
  /* Dealing with the input SFT -> timeseries conversion and whether it has already been done and is buffered */

  /* generate bandpassed and downsampled timeseries for each detector                         */
  /* we only ever do this once for a given dataset so we read it from the buffer if it exists */
  /* in future implementations we will pass this directly to the function instead of SFTs     */

  /* check if there is an not existing timeseries and if the start time in the buffer does not match the start time of the SFTs */
  if ( !cfBuffer->multiTimeseries || ( XLALGPSCmp(&cfBuffer->segstart,&multiSFTs->data[0]->data[0].epoch) != 0) ) {

    XLALPrintInfo ("*** New segment : recomputing timeseries and detstates\n");
    if ( !cfBuffer->multiTimeseries) XLALPrintInfo("timeseries pointer was null\n");
    if ( XLALGPSCmp(&cfBuffer->segstart,&multiSFTs->data[0]->data[0].epoch) != 0)
      XLALPrintInfo("segstart changed from %d to %d\n",cfBuffer->segstart.gpsSeconds,multiSFTs->data[0]->data[0].epoch.gpsSeconds);

    /* if there was no existing timeseries we need to recompute the timeseries from the SFTs */
    /* generate multiple coincident timeseries - one for each detector spanning start -> end */
    /* we need each timeseries to span the exact same amount of time and to start at the same time */
    /* because for the multi-detector Fstat we need frequency bins to be coincident */
    /* The memory allocated here is freed when the buffer is cleared in the calling program */

    /* generate a new timeseries from the input SFTs */
    if ( ( multiTimeseries = XLALMultiSFTVectorToCOMPLEX8TimeSeries(multiSFTs)) == NULL ) {
      XLALPrintError("\nXLALMultiSFTVectorToCOMPLEX8TimeSeries() failed with error = %d\n\n", xlalErrno );
      ABORT ( status, COMPUTEFSTATRSC_EXLAL, COMPUTEFSTATRSC_MSGEXLAL );
    }

    /* recompute the multidetector states for the possibly time shifted SFTs */
    /* the function XLALMultiSFTVectorToCOMPLEX8TimeSeries may have shifted the SFT start times around */
    /* and since these times will be used later on for the resampling we also need to recompute the */
    /* MultiDetectorStates because their timestamps are used later on to compute SSB times which we */
    /* need to be accurate at the midpoints of the SFTs.  Understand ? */

    /* recompute the multiDetStates for the new SFT start times */
    TRY ( LALGetMultiDetectorStates ( status->statusPtr, &multiDetStates, multiSFTs, params->edat ), status );

    /* set all other segment dependent quantity pointers to NULL */
    /* this will basically mean that we will have to compute all sky dependent quantities again */
    XLALEmptyComputeFBuffer_RS( cfBuffer );

    /* buffer the multitimeseries, detstates and the current start time of the input data */
    cfBuffer->multiTimeseries = multiTimeseries;
    cfBuffer->multiDetStates = multiDetStates;
    cfBuffer->segstart.gpsSeconds = multiSFTs->data[0]->data[0].epoch.gpsSeconds;
    cfBuffer->segstart.gpsNanoSeconds = multiSFTs->data[0]->data[0].epoch.gpsNanoSeconds;

  }  /* if (!cfBuffer->multiTimeseries || (buffered-start != SFT-start) ) */

  /* End of the SFT -> timeseries buffering checks                                                            */
  /************************************************************************************************************/

  /************************************************************************************************************/
  /* Dealing with sky position dependent quantities and buffering them                                        */

  /* if the sky position has changed or if any of the sky position dependent quantities are not buffered
     i.e the multiDetstates, the multiAMcoefficients, the multiSSB times and the resampled multiTimeSeries Fa and Fb,
     then we need to recompute these and buffer them */
  if ( (cfBuffer->Alpha != doppler->Alpha )                                 /* and alpha hasn't changed */
       || ( cfBuffer->Delta != doppler->Delta )                             /* and delta hasn't changed */
       || ( cfBuffer->multiAMcoef == NULL )                                 /* and we have a buffered multiAMcoefficents */
       || ( cfBuffer->multiSSB == NULL )                                    /* and we have buffered multiSSB times */
       || ( cfBuffer->multiFa_resampled == NULL )                           /* and we have buffered multiFa_resampled  */
       || ( cfBuffer->multiFb_resampled == NULL )                           /* and we have multiFb_resampled */
       )
    {

      XLALPrintInfo("*** New sky position : recomputing SSB times, AM coefficients and Fa and Fb\n");

      /* temporary timeseries used to store the unbarycentred antenna weighted Fa and Fb timeseries */
      MultiCOMPLEX8TimeSeries *multiFa = NULL;
      MultiCOMPLEX8TimeSeries *multiFb = NULL;

      XLALPrintInfo("freeing all of the sky dependent buffered quantities\n");
      /* free all of the non-null quantities we're about to regenerate and point them to NULL */
      /* if ( cfBuffer->multiAMcoef ) XLALDestroyMultiAMCoeffs( cfBuffer->multiAMcoef ); */
/*       if ( cfBuffer->multiSSB ) XLALDestroyMultiSSBtimes( cfBuffer->multiSSB ); */
/*       if ( cfBuffer->multiFa_resampled) XLALDestroyMultiCOMPLEX8TimeSeries( cfBuffer->multiFa_resampled ); */
/*       if ( cfBuffer->multiFb_resampled) XLALDestroyMultiCOMPLEX8TimeSeries( cfBuffer->multiFb_resampled ); */
/*       cfBuffer->multiAMcoef = NULL; */
/*       cfBuffer->multiSSB = NULL;  */
/*       cfBuffer->multiFa_resampled = NULL;  */
/*       cfBuffer->multiFb_resampled = NULL;  */

      /* compute the SSB times corresponding to the midpoints of each SFT for the current sky position for all detectors */
      skypos.system = COORDINATESYSTEM_EQUATORIAL;
      skypos.longitude = doppler->Alpha;
      skypos.latitude  = doppler->Delta;
      if ( (multiSSB = XLALGetMultiSSBtimes ( cfBuffer->multiDetStates, skypos, doppler->refTime, params->SSBprec )) == NULL )
        {
          XLALPrintError("XLALGetMultiSSBtimes() failed with error = %d\n\n", xlalErrno );
          ABORT ( status, COMPUTEFSTATC_EXLAL, COMPUTEFSTATC_MSGEXLAL );
        }

      /* compute the AM parameters for each detector */
      LALGetMultiAMCoeffs ( status->statusPtr, &multiAMcoef, cfBuffer->multiDetStates, skypos );
      BEGINFAIL ( status ) {
	XLALDestroyMultiSSBtimes ( multiSSB );
      } ENDFAIL (status);

      /* noise-weight Antenna-patterns and compute A,B,C */
      if ( XLALWeightMultiAMCoeffs ( multiAMcoef, multiWeights ) != XLAL_SUCCESS ) {
	XLALPrintError("\nXLALWeightMultiAMCoeffs() failed with error = %d\n\n", xlalErrno );
	ABORT ( status, COMPUTEFSTATRSC_EXLAL, COMPUTEFSTATRSC_MSGEXLAL );
      }

      /* Generate a(t) and b(t) weighted heterodyned downsampled timeseries */
      if ( XLALAntennaWeightMultiCOMPLEX8TimeSeries ( &multiFa, &multiFb, cfBuffer->multiTimeseries,
						      multiAMcoef, multiSFTs) != XLAL_SUCCESS ) {
	XLALPrintError("\nXLALAntennaWeightMultiCOMPLEX8TimeSeries() failed with error = %d\n\n", xlalErrno );
	ABORT ( status, COMPUTEFSTATRSC_EXLAL, COMPUTEFSTATRSC_MSGEXLAL );
      }

      /* Perform barycentric resampling on the multi-detector timeseries */
      if ( XLALBarycentricResampleMultiCOMPLEX8TimeSeries ( &multiFa_resampled, &multiFb_resampled,
							    multiFa, multiFb, multiSSB, multiSFTs,
							    df_out) != XLAL_SUCCESS ) {
	XLALPrintError("\nXLALBarycentricResampleMultiCOMPLEX8TimeSeries() failed with error = %d\n\n", xlalErrno );
	ABORT ( status, COMPUTEFSTATRSC_EXLAL, COMPUTEFSTATRSC_MSGEXLAL );
      }

      /* free multiFa and MultiFb - we won't need them again since we're storing the resampled versions */
      XLALDestroyMultiCOMPLEX8TimeSeries ( multiFa );
      XLALDestroyMultiCOMPLEX8TimeSeries ( multiFb );

      /* buffer all new sky position dependent values - after clearing them */
      cfBuffer->Alpha = doppler->Alpha;
      cfBuffer->Delta = doppler->Delta;

      XLALDestroyMultiSSBtimes( cfBuffer->multiSSB );
      cfBuffer->multiSSB = multiSSB;

      XLALDestroyMultiAMCoeffs( cfBuffer->multiAMcoef );
      cfBuffer->multiAMcoef = multiAMcoef;

      XLALDestroyMultiCOMPLEX8TimeSeries( cfBuffer->multiFa_resampled );
      XLALDestroyMultiCOMPLEX8TimeSeries( cfBuffer->multiFb_resampled );
      cfBuffer->multiFa_resampled = multiFa_resampled;
      cfBuffer->multiFb_resampled = multiFb_resampled;

    } /* could not reuse previously buffered quantities */

  /* End of the sky position dependent quantity buffering                                                     */
  /************************************************************************************************************/

  /* store AM coefficient integrals in local variables */
  if ( cfBuffer->multiAMcoef ) {
    Ad = cfBuffer->multiAMcoef->Mmunu.Ad;
    Bd = cfBuffer->multiAMcoef->Mmunu.Bd;
    Cd = cfBuffer->multiAMcoef->Mmunu.Cd;
    Dd_inv = 1.0 / cfBuffer->multiAMcoef->Mmunu.Dd;
  }
  else {
    XLALPrintError ( "Programming error: 'multiAMcoef' not available!\n");
    ABORT ( status, COMPUTEFSTATRSC_ENULL, COMPUTEFSTATRSC_MSGENULL );
  }

  // *copy* complete resampled multi-complex8 timeseries so we can apply spindown-corrections to it
  MultiCOMPLEX8TimeSeries *multiFa_spin, *multiFb_spin;
  if ( ( multiFa_spin = XLALDuplicateMultiCOMPLEX8TimeSeries ( cfBuffer->multiFa_resampled ) ) == NULL )
    ABORT ( status, COMPUTEFSTATRSC_EXLAL, COMPUTEFSTATRSC_MSGEXLAL );
  if ( ( multiFb_spin = XLALDuplicateMultiCOMPLEX8TimeSeries ( cfBuffer->multiFb_resampled ) ) == NULL )
    ABORT ( status, COMPUTEFSTATRSC_EXLAL, COMPUTEFSTATRSC_MSGEXLAL );


  /* compute the fractional bin offset between the user requested initial frequency */
  /* and the closest output frequency bin */
  REAL8 diff = cfBuffer->multiTimeseries->data[0]->f0 - doppler->fkdot[0]; /* the difference between the new timeseries heterodyne frequency and the user requested lowest frequency */
  INT4 bins = (INT4)round( diff / fstatVector->deltaF );           /* the rounded number of output frequency bins difference */
  REAL8 shift = diff - fstatVector->deltaF * bins;                       /* the fractional bin frequency offset */

  /* shift the timeseries by a fraction of a frequency bin so that user requested frequency is exactly resolved */
  if (shift != 0.0) {
    if ( XLALFrequencyShiftMultiCOMPLEX8TimeSeries ( &multiFa_spin, shift ) != XLAL_SUCCESS ) {
      XLALPrintError("\nXLALMultiSFTVectorToCOMPLEX8TimeSeries() failed with error = %d\n\n", xlalErrno );
      ABORT ( status, COMPUTEFSTATRSC_EXLAL, COMPUTEFSTATRSC_MSGEXLAL );
    }
    if ( XLALFrequencyShiftMultiCOMPLEX8TimeSeries ( &multiFb_spin, shift ) != XLAL_SUCCESS ) {
      XLALPrintError("\nXLALMultiSFTVectorToCOMPLEX8TimeSeries() failed with error = %d\n\n", xlalErrno );
      ABORT ( status, COMPUTEFSTATRSC_EXLAL, COMPUTEFSTATRSC_MSGEXLAL );
    }
  }

  /* apply spin derivitive correction to resampled timeseries */
  /* this function only applies a correction if there are any non-zero spin derivitives */
  if ( XLALSpinDownCorrectionMultiFaFb ( &multiFa_spin, &multiFb_spin, doppler ) != XLAL_SUCCESS ) {
    XLALPrintError("\nXLALSpinDownCorrectionMultiFaFb() failed with error = %d\n\n", xlalErrno );
    ABORT ( status, COMPUTEFSTATRSC_EXLAL, COMPUTEFSTATRSC_MSGEXLAL );
  }

  /**********************************************************************************/
  /* we now compute the FFTs of the resampled functions Fa and Fb for each detector */
  /* and combine them into the multi-detector F-statistic */

  /* we use the first detector Fa time series to obtain the number of time samples and the sampling time */
  /* these should be the same for all Fa and Fb timeseries */
  numSamples = multiFa_spin->data[0]->data->length;
  dt = multiFa_spin->data[0]->deltaT;

  /* allocate memory for Fa(f) and Fb(f) and individual detector FFT outputs */
  if ( (Faf_resampled = XLALCreateCOMPLEX8Vector(numSamples)) == NULL )
    ABORT ( status, xlalErrno, "XLALCreateCOMPLEX8Vector() failed!\n");
  if ( (Fbf_resampled = XLALCreateCOMPLEX8Vector(numSamples)) == NULL )
    ABORT ( status, xlalErrno, "XLALCreateCOMPLEX8Vector() failed!\n");
  if ( (outa = XLALCreateCOMPLEX8Vector(numSamples)) == NULL )
    ABORT ( status, xlalErrno, "XLALCreateCOMPLEX8Vector() failed!\n");
  if ( (outb = XLALCreateCOMPLEX8Vector(numSamples)) == NULL )
    ABORT ( status, xlalErrno, "XLALCOMPLEX8VectorFFT() failed!\n");
  if ( params->returnSingleF ) {
    if ( (outaSingle = XLALCreateCOMPLEX8Vector(numSamples)) == NULL )
      ABORT ( status, xlalErrno, "XLALCreateCOMPLEX8Vector() failed!\n");
    if ( (outbSingle = XLALCreateCOMPLEX8Vector(numSamples)) == NULL )
      ABORT ( status, xlalErrno, "XLALCreateCOMPLEX8Vector() failed!\n");
  }

  /* initialise output vectors to zero since it will be added to */
  memset(Faf_resampled->data,0,numSamples*sizeof(COMPLEX8));
  memset(Fbf_resampled->data,0,numSamples*sizeof(COMPLEX8));

   /* make forwards FFT plan - this will be re-used for each detector */
  if ( (pfwd = XLALCreateCOMPLEX8FFTPlan(numSamples,1,0)) == NULL )
    ABORT ( status, xlalErrno, "XLALCreateCOMPLEX8FFTPlan() failed!\n");

  /* loop over detectors */
  for (i=0;i<numDetectors;i++) {

    COMPLEX8Vector *ina = multiFa_spin->data[i]->data; /* we point the input to the current detector Fa timeseries */
    COMPLEX8Vector *inb = multiFb_spin->data[i]->data; /* we point the input to the current detector Fb timeseries */

    /* initialise output vectors to zero for safety */
    memset(outa->data,0,numSamples*sizeof(COMPLEX8));
    memset(outb->data,0,numSamples*sizeof(COMPLEX8));
    if ( params->returnSingleF ) {
      memset(outaSingle->data,0,numSamples*sizeof(COMPLEX8));
      memset(outbSingle->data,0,numSamples*sizeof(COMPLEX8));
    }

    /* Fourier transform the resampled Fa(t) and Fb(t) */
    if (XLALCOMPLEX8VectorFFT(outa,ina,pfwd)!= XLAL_SUCCESS)
      ABORT ( status, xlalErrno, "XLALCOMPLEX8VectorFFT() failed!\n");
    if (XLALCOMPLEX8VectorFFT(outb,inb,pfwd)!= XLAL_SUCCESS)
      ABORT ( status, xlalErrno, "XLALCOMPLEX8VectorFFT() failed!\n");

    /*  add to summed Faf and Fbf and normalise by dt */
    for (j=0;j<numSamples;j++) {
      Faf_resampled->data[j] += (outa->data[j] * ((REAL4) dt));
      Fbf_resampled->data[j] += (outb->data[j] * ((REAL4) dt));
    }

    /* compute single-IFO F-stats, if requested */
    if ( params->returnSingleF )
      {
       if (params->buffer == NULL ) {
         AdX = multiAMcoef->data[i]->A;
         BdX = multiAMcoef->data[i]->B;
         CdX = multiAMcoef->data[i]->C;
         DdX_inv = 1.0 / multiAMcoef->data[i]->D;
       } else {
         AdX = cfBuffer->multiAMcoef->data[i]->A;
         BdX = cfBuffer->multiAMcoef->data[i]->B;
         CdX = cfBuffer->multiAMcoef->data[i]->C;
         DdX_inv = 1.0 / cfBuffer->multiAMcoef->data[i]->D;
       }

       /* normalize by dt */
       for (UINT4 l=0; l < numSamples; l++) {
         outaSingle->data[l] = (outa->data[l] * ((REAL4) dt));
         outbSingle->data[l] = (outb->data[l] * ((REAL4) dt));
       }

       /* the complex FFT output is shifted such that the heterodyne frequency is at DC */
       /* we need to shift the negative frequencies to before the positive ones */
       if ( XLALFFTShiftCOMPLEX8Vector(&outaSingle) != XLAL_SUCCESS )
         ABORT ( status, xlalErrno, "XLALCOMPLEX8VectorFFT() failed!\n");
       if ( XLALFFTShiftCOMPLEX8Vector(&outbSingle) != XLAL_SUCCESS )
         ABORT ( status, xlalErrno, "XLALCOMPLEX8VectorFFT() failed!\n");

       /* define new initial frequency of the frequency domain representations of Fa and Fb */
       /* before the shift the zero bin was the heterodyne frequency */
       /* now we've shifted it by N - NhalfPosDC(N) bins */
       f0_shifted_single = multiFa_spin->data[i]->f0 - NhalfNeg(numSamples) * df_out;

       /* define number of bins offset from the internal start frequency bin to the user requested bin */
       UINT4 offset_single = floor(0.5 + (doppler->fkdot[0] - f0_shifted_single)/fstatVector->deltaF);

       /* compute final single-IFO F-stat */
       UINT4 numFreqBins = (fstatVector->data->length)/(numDetectors + 1);
       for (UINT4 m = 0; m < numFreqBins; m++) {
         UINT4 idy = m + offset_single;
         fstatVector->data->data[((i+1)*numFreqBins) + m] = DdX_inv * (  BdX * (SQ(crealf(outaSingle->data[idy])) + SQ(cimagf(outaSingle->data[idy])) )
                                  + AdX * ( SQ(crealf(outbSingle->data[idy])) + SQ(cimagf(outbSingle->data[idy])) )
                                  - 2.0 * CdX *( crealf(outaSingle->data[idy]) * crealf(outbSingle->data[idy]) + cimagf(outaSingle->data[idy]) * cimagf(outbSingle->data[idy]) )
                                   );
       } /* end loop over samples */
    } /* if returnSingleF */

  } /* end loop over detectors */

/* the complex FFT output is shifted such that the heterodyne frequency is at DC */
  /* we need to shift the negative frequencies to before the positive ones */
  if ( XLALFFTShiftCOMPLEX8Vector(&Faf_resampled) != XLAL_SUCCESS )
    ABORT ( status, xlalErrno, "XLALCOMPLEX8VectorFFT() failed!\n");
  if ( XLALFFTShiftCOMPLEX8Vector(&Fbf_resampled) != XLAL_SUCCESS )
    ABORT ( status, xlalErrno, "XLALCOMPLEX8VectorFFT() failed!\n");

  /* define new initial frequency of the frequency domain representations of Fa and Fb */
  /* before the shift the zero bin was the heterodyne frequency */
  /* now we've shifted it by N - NhalfPosDC(N) bins */
  f0_shifted = multiFa_spin->data[0]->f0 - NhalfNeg(numSamples) * df_out;

  /* loop over requested output frequencies and construct F *NOT* 2F */
  {

    /* define number of bins offset from the internal start frequency bin to the user requested bin */
    UINT4 offset = floor(0.5 + (doppler->fkdot[0] - f0_shifted)/fstatVector->deltaF);
    if ( params->returnSingleF ) {
      kmax = (fstatVector->data->length)/(numDetectors + 1);
    } else {
      kmax = fstatVector->data->length;
    }
    for (k=0; k < kmax; k++) {

      UINT4 idx = k + offset;
      /* ----- compute final Fstatistic-value ----- */

      /* NOTE: the data MUST be normalized by the DOUBLE-SIDED PSD (using LALNormalizeMultiSFTVect),
       * therefore there is a factor of 2 difference with respect to the equations in JKS, which
       * where based on the single-sided PSD.
       */
      fstatVector->data->data[k] = Dd_inv * (
                   Bd * (SQ(crealf(Faf_resampled->data[idx])) + SQ(cimagf(Faf_resampled->data[idx])))
					       + Ad * (SQ(crealf(Fbf_resampled->data[idx])) + SQ(cimagf(Fbf_resampled->data[idx])))
					       - 2.0 * Cd *( crealf(Faf_resampled->data[idx]) * crealf(Fbf_resampled->data[idx]) +
							                 cimagf(Faf_resampled->data[idx]) * cimagf(Fbf_resampled->data[idx]) )
					       );
    }
  }

  /* free memory not stored in the buffer */
  XLALDestroyCOMPLEX8Vector( Faf_resampled );
  XLALDestroyCOMPLEX8Vector( Fbf_resampled );
  XLALDestroyCOMPLEX8Vector( outa );
  XLALDestroyCOMPLEX8Vector( outb );
  if ( params->returnSingleF ) {
    XLALDestroyCOMPLEX8Vector ( outaSingle );
    XLALDestroyCOMPLEX8Vector ( outbSingle );
  }
  XLALDestroyCOMPLEX8FFTPlan ( pfwd );

  XLALDestroyMultiCOMPLEX8TimeSeries ( multiFa_spin );
  XLALDestroyMultiCOMPLEX8TimeSeries ( multiFb_spin );

  /* IMPORTANT - point the input buffer pointer to the buffered data */
  params->buffer = cfBuffer;

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* ComputeFStatFreqBand_RS() */


/**
 * Destruction of a ComputeFBuffer *contents*,
 * i.e. the multiSSB and multiAMcoeff, while the
 * buffer-container is not freed (which is why it's passed
 * by value and not by reference...)
 */
void
XLALEmptyComputeFBuffer_RS ( ComputeFBuffer_RS *buffer)
{

  if ( buffer->multiSSB ) XLALDestroyMultiSSBtimes( buffer->multiSSB );
  buffer->multiSSB = NULL;
  if ( buffer->multiBinary ) XLALDestroyMultiSSBtimes( buffer->multiBinary );
  buffer->multiBinary = NULL;
  if ( buffer->multiAMcoef) XLALDestroyMultiAMCoeffs( buffer->multiAMcoef );
  buffer->multiAMcoef = NULL;
  if ( buffer->multiCmplxAMcoef) XLALDestroyMultiCmplxAMCoeffs( buffer->multiCmplxAMcoef );
  buffer->multiCmplxAMcoef = NULL;
  if ( buffer->multiTimeseries) XLALDestroyMultiCOMPLEX8TimeSeries( buffer->multiTimeseries );
  buffer->multiTimeseries = NULL;
  if ( buffer->multiFa_resampled) XLALDestroyMultiCOMPLEX8TimeSeries( buffer->multiFa_resampled );
  buffer->multiFa_resampled = NULL;
  if ( buffer->multiFb_resampled) XLALDestroyMultiCOMPLEX8TimeSeries( buffer->multiFb_resampled );
  buffer->multiFb_resampled = NULL;
  if ( buffer->multiDetStates) XLALDestroyMultiDetectorStateSeries( buffer->multiDetStates);
  buffer->multiDetStates = NULL;
  /* if ( buffer ) XLALFree(buffer); */

  return;

} /* XLALEmptyComputeFBuffer_RS() */
