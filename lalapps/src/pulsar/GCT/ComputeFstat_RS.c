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
#include "ComputeFstat_RS.h"
#include "../FDS_isolated/Fstat_v3.h"
#include <lal/ComplexAM.h>
#include <lal/TimeSeries.h>

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
static LALUnit empty_LALUnit;

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
 * Turn the given multiSFTvector into multiple long COMPLEX8TimeSeries, properly dealing with gaps.
 * Memory allocation for the output MultiCOMPLEX8TimeSeries is done within this function.
 *
 * NOTE : We enforce that each detectors timeseries has <b>equal</b> start times and time spans.
 * Also, the input MultiSFTs get <b>modified</b> in place.
 */
MultiCOMPLEX8TimeSeries *XLALMultiSFTVectorToCOMPLEX8TimeSeries (
                     MultiSFTVector *multisfts  /**< [in/out] multi SFT vector, gets modified! */					       )
{
  UINT4 i;
  MultiCOMPLEX8TimeSeries *out = NULL;		/* long time-series corresponding to full set of SFTs */
  LIGOTimeGPS start,end;

  /* check sanity of input */
  if ( !multisfts || (multisfts->length == 0) ) {
    XLALPrintError ("%s: empty multi SFT input!\n", __func__ );
    XLAL_ERROR_NULL (XLAL_EINVAL);
  }
  for (i=0;i<multisfts->length;i++)
    {
      if ( !multisfts->data[i] || (multisfts->data[i]->length == 0) ) {
	XLALPrintError ("%s: empty multiSFT->data[%d] input!\n", __func__,i );
	XLAL_ERROR_NULL (XLAL_EINVAL);
      }
    }

  /* determine the start and end times of the multiSFT observation */
  if ( XLALEarliestMultiSFTsample( &start, multisfts) != XLAL_SUCCESS ) {
    XLALPrintError("%s: Failed to run XLALEarliestMultiSFTsample()\n", __func__ );
    XLAL_ERROR_NULL (XLAL_EFAULT );
  }
  if ( XLALLatestMultiSFTsample( &end, multisfts) != XLAL_SUCCESS ) {
    XLALPrintError("%s: Failed to run XLALLatestMultiSFTsample()\n", __func__ );
    XLAL_ERROR_NULL (XLAL_EFAULT );
  }
  /*
  printf("earliest sample at %d %d\n",start.gpsSeconds,start.gpsNanoSeconds);
  printf("latest sample at %d %d\n",end.gpsSeconds,end.gpsNanoSeconds);
  */

  /* check that earliest is before latest */
  if ( (XLALGPSDiff ( &end, &start ) ) < 0 )
    {
      XLALPrintError ("%s: start time after end time!\n", __func__ );
      XLAL_ERROR_NULL (XLAL_EINVAL);
    }

  /* allocate memory for the output structure */
  if ((out = (MultiCOMPLEX8TimeSeries*)XLALMalloc(sizeof(MultiCOMPLEX8TimeSeries))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", __func__, sizeof(MultiCOMPLEX8TimeSeries));
    XLAL_ERROR_NULL (XLAL_ENOMEM);
  }
  out->length = multisfts->length;
  if ((out->data = XLALMalloc(out->length*sizeof(COMPLEX8TimeSeries*))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", __func__, out->length*sizeof(COMPLEX8TimeSeries*));
    XLAL_ERROR_NULL (XLAL_ENOMEM);
  }

  /* loop over detectors */
  for (i=0;i<multisfts->length;i++) {

    /* call XLALSFTVectorToCOMPLEX8TimeSeries for each detector */
    if ((out->data[i] = XLALSFTVectorToCOMPLEX8TimeSeries(multisfts->data[i],&start,&end)) == NULL) {
      XLALPrintError ("%s: Failed to run XLALSFTVectorToCOMPLEX8TimeSeries()\n", __func__);
      XLAL_ERROR_NULL (XLAL_EFAULT );
    }

  }

  return out;

} /* XLALMultiSFTVectorToCOMPLEX8TimeSeries() */



/**
 * Finds the earliest timestamp in a multi-SFT data structure
 *
 */
int XLALEarliestMultiSFTsample ( LIGOTimeGPS *out,              /**< [out] earliest GPS time */
				 MultiSFTVector *multisfts      /**< [in] multi SFT vector */
				 )
{
  UINT4 i,j;

  /* check sanity of input */
  if ( !multisfts || (multisfts->length == 0) )
    {
      XLALPrintError ("%s: empty multiSFT input!\n", __func__ );
      XLAL_ERROR (XLAL_EINVAL);
    }
  for (i=0;i<multisfts->length;i++)
    {
      if ( !multisfts->data[i] || (multisfts->data[i]->length == 0) )
	{
	  XLALPrintError ("%s: empty multiSFT->data[%d] input!\n", __func__,i );
	  XLAL_ERROR (XLAL_EINVAL);
	}
    }

  /* initialise the earliest sample value */
  out->gpsSeconds = multisfts->data[0]->data[0].epoch.gpsSeconds;
  out->gpsNanoSeconds = multisfts->data[0]->data[0].epoch.gpsNanoSeconds;

  /* loop over detectors */
  for (i=0;i<multisfts->length;i++) {

    /* loop over all SFTs to determine the earliest SFT epoch */
    for (j=0;j<multisfts->data[i]->length;j++) {

      /* compare current SFT epoch with current earliest */
      if ( (XLALGPSCmp(out,&multisfts->data[i]->data[j].epoch) == 1 ) ) {
	out->gpsSeconds = multisfts->data[i]->data[j].epoch.gpsSeconds;
	out->gpsNanoSeconds = multisfts->data[i]->data[j].epoch.gpsNanoSeconds;
      }

    }

  }

  /* success */
  return XLAL_SUCCESS;

} /* XLALEarliestMultiSFTsample() */

/**
 * Find the time of the end of the latest SFT in a multi-SFT data structure
 */
int XLALLatestMultiSFTsample ( LIGOTimeGPS *out,              /**< [out] latest GPS time */
			       MultiSFTVector *multisfts      /**< [in] multi SFT vector */
			       )
{
  UINT4 i,j;
  SFTtype *firstSFT;
  REAL8 Tsft;

  /* check sanity of input */
  if ( !multisfts || (multisfts->length == 0) )
    {
      XLALPrintError ("%s: empty multiSFT input!\n", __func__ );
      XLAL_ERROR (XLAL_EINVAL);
    }
  for (i=0;i<multisfts->length;i++)
    {
      if ( !multisfts->data[i] || (multisfts->data[i]->length == 0) )
	{
	  XLALPrintError ("%s: empty multiSFT->data[%d] input!\n", __func__,i );
	  XLAL_ERROR (XLAL_EINVAL);
	}
    }

  /* define some useful quantities */
  firstSFT = (multisfts->data[0]->data);        /* a pointer to the first SFT of the first detector */
  Tsft = 1.0 / firstSFT->deltaF;                /* the length of the SFTs in seconds assuming 1/T freq resolution */

  /* initialise the latest sample value */
  out->gpsSeconds = firstSFT->epoch.gpsSeconds;
  out->gpsNanoSeconds = firstSFT->epoch.gpsNanoSeconds;

  /* loop over detectors */
  for (i=0;i<multisfts->length;i++) {

    /* loop over all SFTs to determine the earliest SFT midpoint of the input data in the SSB frame */
    for (j=0;j<multisfts->data[i]->length;j++) {

      /* compare current SFT epoch with current earliest */
      if ( (XLALGPSCmp(out,&multisfts->data[i]->data[j].epoch) == -1 ) ) {
	out->gpsSeconds = multisfts->data[i]->data[j].epoch.gpsSeconds;
	out->gpsNanoSeconds = multisfts->data[i]->data[j].epoch.gpsNanoSeconds;
      }
    }

  }

  /* add length of SFT to the result so that we output the end of the SFT */
  if ( XLALGPSAdd(out,Tsft) == NULL )
    {
      XLALPrintError ("%s: NULL pointer returned from XLALGPSAdd()!\n", __func__ );
      XLAL_ERROR (XLAL_EFAULT);
    }

  /* success */
  return XLAL_SUCCESS;

} /* XLALLatestMultiSFTsample() */


/**
 * Computed the weighted timeseries Fa(t) = x(t).a(t) and Fb(t) = x(t).b(t) for a multi-detector timeseries
 */
int XLALAntennaWeightCOMPLEX8TimeSeries (
            COMPLEX8TimeSeries **Faoft,                   /**< [out] the timeseries weighted by a(t) */
					  COMPLEX8TimeSeries **Fboft,                   /**< [out] the timeseries weighted by b(t) */
					  const COMPLEX8TimeSeries *timeseries,         /**< [in] the input timeseries */
					  const AMCoeffs *AMcoef,                       /**< [in] the AM coefficients */
					  const SFTVector *sfts                         /**< [in] the SFT data */
					  )
{
  UINT4 j,k;

  /* do sanity checks */


  /* local copies */
  REAL8 start = GPS2REAL8(timeseries->epoch);
  REAL8 fHet = timeseries->f0;
  REAL8 deltaT = timeseries->deltaT;
  UINT4 numTimeSamples = timeseries->data->length;
  REAL8 dfSFT = sfts->data[0].deltaF;
  REAL8 Tsft = 1.0 / dfSFT;
  UINT4 nbins = (UINT4)floor(0.5 + Tsft/deltaT);

  /* create empty timeseries structures for Fa(t) and Fb(t) */
  if ( ((*Faoft) = XLALCreateCOMPLEX8TimeSeries ( sfts->data[0].name, &(timeseries->epoch), fHet, deltaT, &empty_LALUnit, numTimeSamples )) == NULL ) {
    XLALPrintError ("%s: XLALCreateCOMPLEX8TimeSeries() for %d timesteps failed! errno = %d!\n", __func__, numTimeSamples, xlalErrno );
    XLAL_ERROR (XLAL_ENOMEM);
  }
  if ( ((*Fboft) = XLALCreateCOMPLEX8TimeSeries ( sfts->data[0].name,&(timeseries->epoch) , fHet, deltaT, &empty_LALUnit, numTimeSamples )) == NULL ) {
    XLALPrintError ("%s: XLALCreateCOMPLEX8TimeSeries() for %d timesteps failed! errno = %d!\n", __func__, numTimeSamples, xlalErrno );
    XLAL_ERROR (XLAL_ENOMEM);
  }
  memset ( (*Faoft)->data->data, 0, numTimeSamples * sizeof(*(*Faoft)->data->data)); 	/* set all time-samples to zero (in case there are gaps) */
  memset ( (*Fboft)->data->data, 0, numTimeSamples * sizeof(*(*Fboft)->data->data)); 	/* set all time-samples to zero (in case there are gaps) */


  /* loop over SFTs */
  for (j=0;j<sfts->length;j++) {

    REAL8 t = GPS2REAL8(sfts->data[j].epoch);                         /* the GPS time at the start of the SFT */
    UINT4 start_index = (UINT4)floor(0.5 + (t - start)/deltaT);                     /* index of timesample corresponding to the start of the SFT */
    REAL8 a = (REAL8)AMcoef->a->data[j];                              /* value of the antenna pattern a(t) at the MID-POINT of the SFT */
    REAL8 b = (REAL8)AMcoef->b->data[j];                              /* value of the antenna pattern b(t) at the MID-POINT of the SFT */

    /* loop over samples from this SFT */
    for (k=0;k<nbins;k++) {

      UINT4 time_index = start_index + k;

      /* weight the complex timeseries by the antenna patterns */
      (*Faoft)->data->data[time_index] = (((REAL4) a) * timeseries->data->data[time_index]);
      (*Fboft)->data->data[time_index] = (((REAL4) b) * timeseries->data->data[time_index]);

      }

    }

  /* success */
  return XLAL_SUCCESS;

}

/**
 * Computed the weighted timeseries Fa(t) = x(t).a(t) and Fb(t) = x(t).b(t) for a multi-detector timeseries
 */
int XLALAntennaWeightMultiCOMPLEX8TimeSeries (
                 MultiCOMPLEX8TimeSeries **Faoft,                        /**< [out] the timeseries weighted by a(t) */
					       MultiCOMPLEX8TimeSeries **Fboft,                        /**< [out] the timeseries weighted by b(t) */
					       const MultiCOMPLEX8TimeSeries *multiTimeseries,         /**< [in] the input multi-detector timeseries */
					       const MultiAMCoeffs *multiAMcoef,                       /**< [in] the multi-detector AM coefficients */
					       const MultiSFTVector *multisfts                         /**< [in] the multi-detector SFT data */
					       )
{
  UINT4 i;

  /* do sanity checks */
  if ( !multiTimeseries || (multiTimeseries->length == 0) ) {
    XLALPrintError ("%s: empty multiTimeseries input!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }
  if ( multiTimeseries->length != multisfts->length ) {
    XLALPrintError ("%s: incorrect length of multiTimeseries or multisfts input!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }

  /* allocate memory for the output structure */
  if ( ((*Faoft) = XLALMalloc(sizeof(MultiCOMPLEX8TimeSeries))) == NULL ) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", __func__, sizeof(MultiCOMPLEX8TimeSeries));
    XLAL_ERROR (XLAL_ENOMEM);
  }
  if ( ((*Fboft) = XLALMalloc(sizeof(MultiCOMPLEX8TimeSeries))) == NULL ) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", __func__, sizeof(MultiCOMPLEX8TimeSeries));
    XLAL_ERROR (XLAL_ENOMEM);
  }

  (*Faoft)->length = multisfts->length;
  (*Fboft)->length = multisfts->length;

  if (((*Faoft)->data = XLALMalloc((multisfts->length)*sizeof(COMPLEX8TimeSeries*))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", __func__, (multisfts->length)*sizeof(COMPLEX8TimeSeries*));
    XLAL_ERROR (XLAL_ENOMEM);
  }
  if (((*Fboft)->data = XLALMalloc((multisfts->length)*sizeof(COMPLEX8TimeSeries*))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", __func__, (multisfts->length)*sizeof(COMPLEX8TimeSeries*));
    XLAL_ERROR (XLAL_ENOMEM);
  }


  /* loop over detectors */
  for (i=0;i<multisfts->length;i++) {

    /* point to current detector params */
    COMPLEX8TimeSeries *timeseries = multiTimeseries->data[i];
    AMCoeffs *AMcoef = multiAMcoef->data[i];
    SFTVector *SFTs = multisfts->data[i];

    if ( XLALAntennaWeightCOMPLEX8TimeSeries(&((*Faoft)->data[i]),&((*Fboft)->data[i]),timeseries,AMcoef,SFTs) != XLAL_SUCCESS ) {
      XLALPrintError("\nXLALAntennaWeightMultiCOMPLEX8TimeSeries() failed with error = %d\n\n", xlalErrno );
      XLAL_ERROR (XLAL_EFAULT);
    }

  }

  /* success */
  return XLAL_SUCCESS;

}

/**
 * Performs barycentric resampling on a multi-detector timeseries
 */
int XLALBarycentricResampleMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries **Faoft_RS,                         /**< [out] the resampled timeseries Fa(t_SSB) */
						     MultiCOMPLEX8TimeSeries **Fboft_RS,                         /**< [out] the resampled timeseries Fb(t_SSB) */
						     const MultiCOMPLEX8TimeSeries *Faoft,                       /**< [in] the detector frame timeseries Fa(t) */
						     const MultiCOMPLEX8TimeSeries *Fboft,                       /**< [in] the detector frame timeseries Fb(t) */
						     const MultiSSBtimes *multiSSB,                              /**< [in] the multi-detector SSB times data */
						     const MultiSFTVector *multiSFTs,                            /**< [in] the multi-detector SFT data */
						     const REAL8 deltaF                                          /**< [in] the user defined output frequency resolution */
						     )
{
  UINT4 i;
  LIGOTimeGPS earliest,latest;
  UINT4 numTimeSamples;
  UINT4 numDetectors;
  REAL8 Tspan;
  REAL8 deltaT;
  REAL8 Teff;
  REAL8 fHet;
  REAL8 Tsft;

  /* do sanity checks on input */
  if ( !Faoft || (Faoft->length == 0) ) {
    XLALPrintError ("%s: empty Faoft input!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }
  numDetectors = Faoft->length;

  if ( !Fboft || (Fboft->length != numDetectors ) ) {
    XLALPrintError ("%s: empty or invalid length of Fboft input!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }

  /* check that SFTs exist */
  if ( !multiSFTs || (multiSFTs->length != numDetectors) ) {
    XLALPrintError ("%s: empty or invalid length of SFTs input!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }

  /* check if the SSB times vector has the correct number of entries - one for each SFT */
  if ( !multiSSB || (multiSSB->length != numDetectors ) ) {
    XLALPrintError ("%s: empty or incorrect length of multiSSB input!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }

  /* define the length of an SFT (assuming 1/T resolution) */
  Tsft = 1.0 / multiSFTs->data[0]->data[0].deltaF;

  /* find earliest SSB time */
  if ( (XLALEarliestMultiSSBtime (&earliest,multiSSB,Tsft)) != XLAL_SUCCESS ) {
    XLALPrintError("\nXLALEarliestMultiSSBtime() failed with error = %d\n\n", xlalErrno );
    XLAL_ERROR (XLAL_EFAULT);
  }

  /* find latest SSB time */
  if ( (XLALLatestMultiSSBtime (&latest,multiSSB,Tsft)) != XLAL_SUCCESS ) {
    XLALPrintError("\nXLALLatestMultiSSBtime() failed with error = %d\n\n", xlalErrno );
    XLAL_ERROR (XLAL_EFAULT);
  }

  /* determine resampled timeseries parameters */
  Tsft = 1.0 / multiSFTs->data[0]->data[0].deltaF;         /* define the length of an SFT (assuming 1/T resolution) */
  Tspan = XLALGPSDiff(&latest,&earliest);                  /* the time span from the earliest sample to the end of the latest sample over *all* detectors */
  deltaT = Faoft->data[0]->deltaT;                         /* the sample rate of the downsampled detector frame timeseries */
  Teff = 1.0/deltaF;                                       /* the effective observation time based on the requested frequency resolution (for zero padding) */
  fHet = Faoft->data[0]->f0;                               /* the input timeseries heterodyne frequency */

  /* if the requested frequency resolution gives an effective observation time less than the actual data span then we exit */
  if (Tspan > Teff) {
    XLALPrintError ("%s: requested frequency resolution too coarse for resampling!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }

  /* redefine sample rate and compute number of samples in the new timeseries */
  numTimeSamples = (UINT4)ceil(Teff/deltaT);      /* we use ceil() so that we artificially widen the band rather than reduce it */
  deltaT = Teff/numTimeSamples;

  /* allocate memory for the output resampled timeseries for Fa and Fb */
  if (((*Faoft_RS) = XLALMalloc(sizeof(MultiCOMPLEX8TimeSeries))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", __func__, sizeof(MultiCOMPLEX8TimeSeries));
    XLAL_ERROR (XLAL_ENOMEM);
  }
  if (((*Fboft_RS) = XLALMalloc(sizeof(MultiCOMPLEX8TimeSeries))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", __func__, sizeof(MultiCOMPLEX8TimeSeries));
    XLAL_ERROR (XLAL_ENOMEM);
  }
  (*Faoft_RS)->length = multiSSB->length;
  (*Fboft_RS)->length = multiSSB->length;

  if (((*Faoft_RS)->data = XLALMalloc(multiSSB->length*sizeof(COMPLEX8TimeSeries*))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", __func__, multiSSB->length*sizeof(COMPLEX8TimeSeries*));
    XLAL_ERROR (XLAL_ENOMEM);
  }
  if (((*Fboft_RS)->data = XLALMalloc(multiSSB->length*sizeof(COMPLEX8TimeSeries*))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", __func__, multiSSB->length*sizeof(COMPLEX8TimeSeries*));
    XLAL_ERROR (XLAL_ENOMEM);
  }

  /* loop over detectors */
  for (i=0;i<multiSSB->length;i++) {

    /* shorthand pointers */
    SSBtimes *SSB = multiSSB->data[i];
    SFTVector *SFTs = multiSFTs->data[i];
    COMPLEX8TimeSeries *Fa = Faoft->data[i];
    COMPLEX8TimeSeries *Fb = Fboft->data[i];

    /* create empty timeseries structures for the resampled Fa(t) and Fb(t) */
    if ( ((*Faoft_RS)->data[i] = XLALCreateCOMPLEX8TimeSeries ( Faoft->data[i]->name, &earliest, fHet, deltaT, &empty_LALUnit, numTimeSamples )) == NULL ) {
      XLALPrintError ("%s: XLALCreateCOMPLEX8TimeSeries() for %d timesteps failed! errno = %d!\n", __func__, numTimeSamples, xlalErrno );
      XLAL_ERROR (XLAL_ENOMEM);
    }
    if ( ((*Fboft_RS)->data[i] = XLALCreateCOMPLEX8TimeSeries ( Fboft->data[i]->name, &earliest , fHet, deltaT, &empty_LALUnit, numTimeSamples )) == NULL ) {
      XLALPrintError ("%s: XLALCreateCOMPLEX8TimeSeries() for %d timesteps failed! errno = %d!\n", __func__, numTimeSamples, xlalErrno );
      XLAL_ERROR (XLAL_ENOMEM);
    }
    memset ( (*Faoft_RS)->data[i]->data->data, 0, numTimeSamples * sizeof(*(*Faoft_RS)->data[i]->data->data)); 	/* set all time-samples to zero (in case there are gaps) */
    memset ( (*Fboft_RS)->data[i]->data->data, 0, numTimeSamples * sizeof(*(*Fboft_RS)->data[i]->data->data)); 	/* set all time-samples to zero (in case there are gaps) */

    /* perform resampling on current detector timeseries */
    if ( XLALBarycentricResampleCOMPLEX8TimeSeries(&((*Faoft_RS)->data[i]),&((*Fboft_RS)->data[i]),Fa,Fb,SSB,SFTs) != XLAL_SUCCESS ) {
      XLALPrintError("\nXLALBarycentricResampleCOMPLEX8TimeSeries() failed with error = %d\n\n", xlalErrno );
      XLAL_ERROR (XLAL_EFAULT);
    }

  }

  /* success */
  return XLAL_SUCCESS;

}

/**
 * Performs barycentric resampling on a COMPLEX8TimeSeries
 * We expect that the output timeseries has already been allocated correctly.
 *
 */
int XLALBarycentricResampleCOMPLEX8TimeSeries ( COMPLEX8TimeSeries **Faoft_RS,                         /**< [out] the resampled timeseries Fa(t_SSB) */
						COMPLEX8TimeSeries **Fboft_RS,                         /**< [out] the resampled timeseries Fb(t_SSB) */
						const COMPLEX8TimeSeries *Faoft,                       /**< [in] the input detector frame timeseries Fa(t) */
						const COMPLEX8TimeSeries *Fboft,                       /**< [in] the input detector frame timeseries Fb(t) */
						const SSBtimes *SSB,                                   /**< [in] the SSB times at the midpoints of each SFT */
						const SFTVector *SFTs                                  /**< [in] the SFT data */
						)
{
  UINT4 i,j,k;
  REAL8Vector *detectortimes = NULL;                /* used to store a vector of *non-uniform* time values in the detector frame (used for interpolation) */
  UINT4 FAFB_LENGTH = 4;                          /* the total number of elements in an instance of Fa and Fb */
  REAL8Vector* FaFb[FAFB_LENGTH];                 /* used to store Fa and Fb real and imaginary parts in 4 seperate vectors */
  REAL8Vector *t_DET = NULL;                        /* used to store a vector of *uniform* time values in the detector frame (used for interpolation) */
  gsl_spline* spline_FaFb[FAFB_LENGTH];           /* used to store spline coefficients for Fa and Fb real and imaginary parts in 4 seperate vectors */
  UINT4 numTimeSamples_DET;
  UINT4 numSFTs;
  REAL8 Tsft;
  REAL8 refTime;
  REAL8 start_DET;
  REAL8 fHet;
  REAL8 deltaT_DET;
  REAL8 start_SSB;
  REAL8 deltaT_SSB;

  /* do sanity checks on input */
  if ( !Faoft || (Faoft->data->length == 0) ) {
    XLALPrintError ("%s: empty Faoft input!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }
  numTimeSamples_DET = Faoft->data->length;         /* set the number of time samples in the detector frame as that defined in Fa */

  /* check if the Fa and Fb timeseries have equal lengths */
  if ( !Fboft || (Fboft->data->length != numTimeSamples_DET ) ) {
    XLALPrintError ("%s: empty or incorrect length of Fboft input!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }

  /* check that the Fa and Fb epochs are equal */
  if ( (XLALGPSCmp(&Faoft->epoch,&Faoft->epoch) != 0 ) ) {
    XLALPrintError ("%s: Faoft and Fboft epochs do not match!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }

  /* check that the output vectors are not NULL */
  if ( !(*Faoft_RS) || ((*Faoft_RS)->data->length == 0) || !(*Fboft_RS) || ((*Fboft_RS)->data->length == 0) ) {
    XLALPrintError ("%s: empty output vectors Faoft_RS and/or Fboft_RS!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }

  /* check that SFTs exist */
  if ( !SFTs || (SFTs->length == 0) ) {
    XLALPrintError ("%s: empty SFTs input!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }
  numSFTs = SFTs->length;                           /* set the number of SFTs */

  /* check if the SSB times vector has the correct number of entries - one for each SFT */
  if ( !SSB || (SSB->DeltaT->length != numSFTs ) || ( SSB->Tdot->length != numSFTs ) ) {
    XLALPrintError ("%s: empty or incorrect length of SSB input!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }

  /* define some useful shorthands */
  Tsft = 1.0 / SFTs->data[0].deltaF;                  /* define the length of an SFT (assuming 1/T resolution) */
  refTime = GPS2REAL8(SSB->refTime);                   /* set the reftime at which the doppler parameters are defined */
  start_DET = GPS2REAL8(Faoft->epoch);                 /* set the start time of the timeseries at the detector */
  fHet = Faoft->f0;                                    /* set the heterodyne frequency of the input timeseries */
  deltaT_DET = Faoft->deltaT;                          /* set the sampling time at the detector */
  start_SSB = GPS2REAL8((*Faoft_RS)->epoch);           /* set the start time of the resampled output SSB timeseries */
  deltaT_SSB = (*Faoft_RS)->deltaT;                    /* set the sampling time of the resampled output SSB timeseries */
  REAL8 end_DET = start_DET + (numTimeSamples_DET-1) * deltaT_DET;	/* time of last sample in detector timeseries */

  /* allocate memory for the uniformly sampled detector time samples (Fa and Fb real and imaginary) */
  for (i=0;i<FAFB_LENGTH;i++) {
    if ( (FaFb[i] = XLALCreateREAL8Vector(numTimeSamples_DET)) == NULL ) {
      XLALPrintError ("%s: XLALCreateREAL8Vector() failed!  errno = %d!\n", __func__, xlalErrno );
      XLAL_ERROR (XLAL_ENOMEM);
    }
  }
  /* allocate memory for the *uniform* detector time vector required for interpolation */
  if ( (t_DET = XLALCreateREAL8Vector(numTimeSamples_DET)) == NULL ) {
    XLALPrintError ("%s: XLALCreateREAL8Vector() failed!  errno = %d!\n", __func__, xlalErrno );
    XLAL_ERROR (XLAL_ENOMEM);
  }

  /* place the timeseries into REAL8Vectors for gsl to be able to interpolate them */
  /* this is annoying because the data is currently stored in a COMPLEX8 format so we can't just point to it */
  for (j=0;j<numTimeSamples_DET;j++) {
    t_DET->data[j] = start_DET + j*deltaT_DET;              /* fill in the uniform detector time vector */
    FaFb[0]->data[j] = crealf(Faoft->data->data[j]);
    FaFb[1]->data[j] = cimagf(Faoft->data->data[j]);
    FaFb[2]->data[j] = crealf(Fboft->data->data[j]);
    FaFb[3]->data[j] = cimagf(Fboft->data->data[j]);
  }

  /* initialise the gsl spline interpolation for each of the 4 timeseries */
  for (i=0;i<FAFB_LENGTH;i++) {
    if ( (XLALGSLInitInterpolateREAL8Vector(&(spline_FaFb[i]),t_DET,FaFb[i])) != XLAL_SUCCESS ) {
      XLALPrintError("\nXLALGSLInitInterpolateREAL8Vector() failed with error = %d\n\n", xlalErrno );
      XLAL_ERROR (XLAL_EFAULT);
    }
  }

  /* loop over SFTs to compute the detector frame time samples corresponding to uniformly sampled SSB time samples */
  for (j=0;j<numSFTs;j++) {

    REAL8Vector *out_FaFb[4];                                                              /* output vectors for the real and imaginary parts of the resampled Fa and Fb */

    /* define some useful shorthands */
    REAL8 Tdot = SSB->Tdot->data[j];                                                       /* the instantaneous time derivitive dt_SSB/dt_DET at the MID-POINT of the SFT */
    REAL8 SFTmid_SSB = refTime + SSB->DeltaT->data[j];                                     /* MID-POINT time of the SFT at the SSB */
    REAL8 SFTstart_SSB = SFTmid_SSB - 0.5*Tsft*Tdot;                                       /* START time of the SFT at the SSB */
    REAL8 SFTend_SSB = SFTmid_SSB + 0.5*Tsft*Tdot;                                         /* END time of the SFT at the SSB */
    REAL8 SFTstart_DET = GPS2REAL8(SFTs->data[j].epoch);                                   /* START time of the SFT at the detector */
    REAL8 SFTmid_DET = SFTstart_DET + 0.5*Tsft;                                            /* MID-POINT time of the SFT at the detector */

    /* define some indices */
    UINT4 idx_start_SSB = floor(0.5 + (SFTstart_SSB - start_SSB)/deltaT_SSB);              /* the index of the resampled timeseries corresponding to the start of the SFT */
    UINT4 idx_end_SSB = floor(0.5 + (SFTend_SSB - start_SSB)/deltaT_SSB);                  /* the index of the resampled timeseries corresponding to the end of the SFT */
    UINT4 numSamples_SSB = idx_end_SSB - idx_start_SSB + 1;                                /* the number of samples in the SSB for this SFT */

    /* allocate memory for the *non-uniform* detector time samples for this SFT */
    /* have to allocate it inside the loop because it may have different lengths for each SFT */
    if ( (detectortimes = XLALCreateREAL8Vector(numSamples_SSB)) == NULL ) {
      XLALPrintError ("%s: XLALCreateREAL8Vector() failed!  errno = %d!\n", __func__, xlalErrno );
      XLAL_ERROR (XLAL_ENOMEM);
    }

    /* for each time sample in the SSB frame for this SFT we estimate the detector time. */
    /* We use a linear approximation expanding around the midpoint of an SFT where  */
    /* t_DET = SFTmid_DET + (t_SSB - SFTmid_SSB)*dt_DET/dt_SSB */
    for (k=0;k<numSamples_SSB;k++) {

      REAL8 t_SSB = start_SSB + (k+idx_start_SSB)*deltaT_SSB;                                  /* the SSB time of the current resampled time sample */
      detectortimes->data[k] = SFTmid_DET + (t_SSB - SFTmid_SSB)/Tdot;                         /* the approximated DET time of the current resampled time sample */

      /*
       * NOTE: we need to be careful that none of the times falls outside
       * of the range of detector timesamples, in order to avoid problems in the interpolation
       * therefore we truncate the detector-times to fully fall within the detector timeseries span
       */
      if ( detectortimes->data[k] > end_DET )
        {
          detectortimes->data[k] = end_DET;
          XLALPrintWarning ("%s: time-sample jSFT=%d, kSample=%d at t=%f to interpolate is *after* detector-timeseries, nudged back to end (end=%f)\n",
                            __func__, j, k, detectortimes->data[k], end_DET );
        }
      if ( detectortimes->data[k] < start_DET )
        {
          detectortimes->data[k] = start_DET;
          XLALPrintWarning ("%s: time-sample jSFT=%d, kSample=%d at t=%f to interpolate is *before* detector-timeseries, nudged to beginning (start=%f)\n",
                            __func__, j, k, detectortimes->data[k], start_DET );
        }

    } /* for k < numSamples_SSB */

    /* interpolate on the non-uniformly sampled detector time vector for this SFT for re and im parts of Fa and Fb */
    /* this function allocates memory for the output vectors */
    for (i=0;i<FAFB_LENGTH;i++) {
      if ( XLALGSLInterpolateREAL8Vector(&(out_FaFb[i]),detectortimes,spline_FaFb[i]) != XLAL_SUCCESS ) {
	XLALPrintError("\nXLALInterpolateMultiCOMPLEX8TimeSeries() failed with error = %d\n\n", xlalErrno );
	XLAL_ERROR (XLAL_EFAULT);
      }
    }

    /* place these interpolated timeseries into the output */
    /* and apply correction due to non-zero heterodyne frequency of input */
    for (k=0;k<numSamples_SSB;k++) {

      UINT4 idx = k + idx_start_SSB;                                                                     /* the full resampled timeseries index */
      REAL8 tDiff = start_SSB + idx*deltaT_SSB - detectortimes->data[k];                                 /* the difference between t_SSB and t_DET */
      REAL8 cycles = fmod ( fHet*tDiff, 1 );                                                            /* the accumulated heterodyne cycles */
      REAL4 cosphase,sinphase;                                                                           /* the real and imaginary parts of the phase correction */

      /* use a look-up-table for speed to compute real and imaginary phase */
      sin_cos_2PI_LUT ( &sinphase, &cosphase, -cycles );

      /* printf("j = %d t = %6.12f tb = %6.12f tDiff = %6.12f\n",j,detectortimes->data[k],start_SSB + idx*deltaT_SSB,tDiff); */
      (*Faoft_RS)->data->data[idx] = crectf( out_FaFb[0]->data[k]*cosphase - out_FaFb[1]->data[k]*sinphase, out_FaFb[1]->data[k]*cosphase + out_FaFb[0]->data[k]*sinphase );
      (*Fboft_RS)->data->data[idx] = crectf( out_FaFb[2]->data[k]*cosphase - out_FaFb[3]->data[k]*sinphase, out_FaFb[3]->data[k]*cosphase + out_FaFb[2]->data[k]*sinphase );

    }

    /* free memory used for this SFT */
    for (i=0;i<FAFB_LENGTH;i++) XLALDestroyREAL8Vector(out_FaFb[i]);
    XLALDestroyREAL8Vector(detectortimes);

  } /* end loop over SFTs */

  /* free memory */
  for (i=0;i<4;i++) {
    XLALDestroyREAL8Vector(FaFb[i]);
    gsl_spline_free(spline_FaFb[i]);
  }
  XLALDestroyREAL8Vector(t_DET);

  /* success */
  return XLAL_SUCCESS;

}

/**
 * Find the earliest timestamp in a multi-SSB data structure
 */
int XLALEarliestMultiSSBtime ( LIGOTimeGPS *out,              /**< output earliest GPS time */
			       const MultiSSBtimes *multiSSB,      /**< input multi SSB SFT-midpoint timestamps */
			       const REAL8 Tsft                    /**< the length of an SFT */
			       )
{
  UINT4 i,j;
  LIGOTimeGPS t;
  REAL8 delta;

  /* check sanity of input */
  if ( !multiSSB || (multiSSB->length == 0) )
    {
      XLALPrintError ("%s: empty multiSSBtimes input!\n", __func__ );
      XLAL_ERROR (XLAL_EINVAL);
    }
  if ( !multiSSB->data[0] || (multiSSB->data[0]->DeltaT->length == 0) )
    {
      XLALPrintError ("%s: empty multiSSBtimes input!\n", __func__ );
      XLAL_ERROR (XLAL_EINVAL);
    }


  /* initialise the earliest and latest sample value */
  out->gpsSeconds = multiSSB->data[0]->refTime.gpsSeconds;
  out->gpsNanoSeconds = multiSSB->data[0]->refTime.gpsNanoSeconds;
  delta = multiSSB->data[0]->DeltaT->data[0] - 0.5*Tsft*multiSSB->data[0]->Tdot->data[0];
  if ( (XLALGPSAdd(out,delta)) == NULL) {
    XLALPrintError ("%s: XLALGPSAdd() failed!  errno = %d!\n", __func__, xlalErrno );
    XLAL_ERROR (XLAL_ENOMEM);
  }
  /* loop over detectors */
  for (i=0;i<multiSSB->length;i++) {

    /* loop over all SSB times and find the earliest SSB SFT start time */
    for (j=0;j<multiSSB->data[i]->DeltaT->length;j++) {

      /* reset the reference time */
      t.gpsSeconds = multiSSB->data[i]->refTime.gpsSeconds;
      t.gpsNanoSeconds = multiSSB->data[i]->refTime.gpsNanoSeconds;

      /* compute SSB time - we approximate the SFT start time in the SSB as t_mid_SSB - 0.5*Tsft*dt_SSB/dt_det */
      delta = multiSSB->data[i]->DeltaT->data[j] - 0.5*Tsft*multiSSB->data[i]->Tdot->data[j];
      if ( (XLALGPSAdd(&t,delta)) == NULL) {
	XLALPrintError ("%s: XLALGPSAdd() failed!  errno = %d!\n", __func__, xlalErrno );
	XLAL_ERROR (XLAL_ENOMEM);
      }

      /* compare it to the existing earliest */
      if ( (XLALGPSCmp(out,&t) == 1 ) ) {
	out->gpsSeconds = t.gpsSeconds;
	out->gpsNanoSeconds = t.gpsNanoSeconds;
      }

    }

  }

  /* success */
  return XLAL_SUCCESS;

} /* XLALEarliestMultiSSBtime() */

/**
 * Find the latest timestamp in a multi-SSB data structure
 */
int XLALLatestMultiSSBtime ( LIGOTimeGPS *out,                   /**< output latest GPS time */
			     const MultiSSBtimes *multiSSB,      /**< input multi SSB SFT-midpoint timestamps */
			     const REAL8 Tsft                    /**< the length of an SFT */
			     )
{
  UINT4 i,j;
  LIGOTimeGPS t;
  REAL8 delta;

  /* check sanity of input */
  if ( !multiSSB || (multiSSB->length == 0) )
    {
      XLALPrintError ("%s: empty multiSSBtimes input!\n", __func__ );
      XLAL_ERROR (XLAL_EINVAL);
    }
  if ( !multiSSB->data[0] || (multiSSB->data[0]->DeltaT->length == 0) )
    {
      XLALPrintError ("%s: empty multiSSBtimes input!\n", __func__ );
      XLAL_ERROR (XLAL_EINVAL);
    }


  /* initialise the earliest and latest sample value */
  out->gpsSeconds = multiSSB->data[0]->refTime.gpsSeconds;
  out->gpsNanoSeconds = multiSSB->data[0]->refTime.gpsNanoSeconds;
  delta = multiSSB->data[0]->DeltaT->data[0] + 0.5*Tsft*multiSSB->data[0]->Tdot->data[0];
  if ( (XLALGPSAdd(out,delta)) == NULL) {
    XLALPrintError ("%s: XLALGPSAdd() failed!  errno = %d!\n", __func__, xlalErrno );
    XLAL_ERROR (XLAL_ENOMEM);
  }
  /* loop over detectors */
  for (i=0;i<multiSSB->length;i++) {

    /* loop over all SSB times and find the latest SSB SFT start time */
    for (j=0;j<multiSSB->data[i]->DeltaT->length;j++) {

      /* reset the reference time */
      t.gpsSeconds = multiSSB->data[i]->refTime.gpsSeconds;
      t.gpsNanoSeconds = multiSSB->data[i]->refTime.gpsNanoSeconds;

      /* compute SSB time - we approximate the SFT end time in the SSB as t_mid_SSB + 0.5*Tsft*dt_SSB/dt_det */
      delta = multiSSB->data[i]->DeltaT->data[j] + 0.5*Tsft*multiSSB->data[i]->Tdot->data[j];
      if ( (XLALGPSAdd(&t,delta)) == NULL) {
	XLALPrintError ("%s: XLALGPSAdd() failed!  errno = %d!\n", __func__, xlalErrno );
	XLAL_ERROR (XLAL_ENOMEM);
      }

      /* compare it to the existing earliest */
      if ( (XLALGPSCmp(out,&t) == -1 ) ) {
	out->gpsSeconds = t.gpsSeconds;
	out->gpsNanoSeconds = t.gpsNanoSeconds;
      }

    }

  }

  /* success */
  return XLAL_SUCCESS;

} /* XLALLatestMultiSSBtime() */

/**
 * Find the latest timestamp in a multi-SSB data structure
 */
int XLALGSLInterpolateREAL8Vector( REAL8Vector **yi,
				   REAL8Vector *xi,
				   gsl_spline *spline
				   )
 {
   /* check input */

   UINT4 i;
   UINT4 numSamples = xi->length;
   gsl_interp_accel *acc = gsl_interp_accel_alloc();

   /* allocate memory for output vector */
   if ( ((*yi) = XLALCreateREAL8Vector(numSamples)) == NULL ) {
     XLALPrintError ("%s: XLALCreateREAL8Vector() failed!  errno = %d!\n", __func__, xlalErrno );
     XLAL_ERROR (XLAL_ENOMEM);
   }

   /* perform inerpolation */
   for (i=0;i<numSamples;i++) {
     (*yi)->data[i] = gsl_spline_eval(spline,xi->data[i],acc);
   }

   /* free memory */
   gsl_interp_accel_free(acc);

   return XLAL_SUCCESS;

}


/**
 * Find the latest timestamp in a multi-SSB data structure
 */
int XLALGSLInitInterpolateREAL8Vector( gsl_spline **spline,
				       REAL8Vector *x,
				       REAL8Vector *y
				       )
 {
   /* check input */

   UINT4 numSamples_in = x->length;
   REAL8 *xtemp = x->data;
   REAL8 *ytemp = y->data;

   /* compute spline interpolation coefficients */
   if ( ((*spline) = gsl_spline_alloc(gsl_interp_cspline,numSamples_in)) == NULL ) {
     XLALPrintError ("%s: XLALInitGSLInterpolateREAL8Vector() failed!  errno = %d!\n", __func__, xlalErrno );
     XLAL_ERROR (XLAL_ENOMEM);
   }
   gsl_spline_init((*spline),xtemp,ytemp,numSamples_in);

   return XLAL_SUCCESS;

}

/**
 * Shifts an FFT output vector such that the Niquest frequency is the central bin
 */
int XLALFFTShiftCOMPLEX8Vector(COMPLEX8Vector **x)
{
  UINT4 N = (*x)->length;
  UINT4 NQ = NhalfPosDC(N);
  UINT4 NminusNQ = N - NQ;
 /*  printf("NQ = %d NminusNQ = %d N = %d\n",NQ,NminusNQ,N); */

  /* allocate temp memory */
  COMPLEX8 *temp = XLALMalloc(NQ*sizeof(COMPLEX8));

  /* copy the bins 0 -> NQ - 1 to temp memory */
  if ( memcpy(temp,(*x)->data,NQ*sizeof(COMPLEX8)) == NULL ) {
    XLALPrintError ("%s: memcpy() failed!  errno = %d!\n", __func__, xlalErrno );
    XLAL_ERROR (XLAL_ENOMEM);
  }

  /* copy the bins NQ -> N - 1 to the start */
  if (memcpy((*x)->data,&((*x)->data[NQ]),NminusNQ*sizeof(COMPLEX8)) == NULL ) {
    XLALPrintError ("%s: memcpy() failed!  errno = %d!\n", __func__, xlalErrno );
    XLAL_ERROR (XLAL_ENOMEM);
  }

  /* copy the temp bins to the end */
  if (memcpy(&((*x)->data[NminusNQ]),temp,NQ*sizeof(COMPLEX8)) == NULL ) {
    XLALPrintError ("%s: memcpy() failed!  errno = %d!\n", __func__, xlalErrno );
    XLAL_ERROR (XLAL_ENOMEM);
  }

  /* free temp memory */
  XLALFree(temp);

  return XLAL_SUCCESS;

}

/**
 * Multi-detector wrapper for XLALFrequencyShiftCOMPLEX8TimeSeries
 * NOTE: this <b>modifies</b> the MultiCOMPLEX8Timeseries in place
 */
int
XLALFrequencyShiftMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries **x,	/**< [in/out] timeseries to time-shift */
					    const REAL8 shift )	                /**< [in] freq-shift in Hz */
{
  UINT4 i;

  if ( !(*x) || (*x)->length == 0 )
    {
      XLALPrintError ("%s: empty input COMPLEX8timeSeries!\n", __func__ );
      XLAL_ERROR (XLAL_EINVAL);
    }

  /* loop over detectors */
  for ( i=0; i < (*x)->length; i++)
    {

      /* shift the frequency of each detector's data */
      if ( XLALFrequencyShiftCOMPLEX8TimeSeries ( &((*x)->data[i]), shift) != XLAL_SUCCESS ) {
	XLALPrintError ("%s: memcpy() failed!  errno = %d!\n", __func__, xlalErrno );
	XLAL_ERROR (XLAL_EFAULT);
      }

    }

  return XLAL_SUCCESS;

} /* XLALFrequencyShiftMultiCOMPLEX8TimeSeries() */

/**
 * Freq-shift the given COMPLEX8Timeseries by an amount of 'shift' Hz,
 * using the time-domain expression y(t) = x(t) * e^(-i 2pi df t),
 * which shifts x(f) into y(f) = x(f - df)
 *
 * NOTE: this <b>modifies</b> the COMPLEX8TimeSeries in place
 */
int
XLALFrequencyShiftCOMPLEX8TimeSeries ( COMPLEX8TimeSeries **x,	        /**< [in/out] timeseries to time-shift */
				       const REAL8 shift )	        /**< [in] freq-shift in Hz */
{
  UINT4 k;
  REAL8 deltat;

  if ( !(*x) || !(*x)->data )
    {
      XLALPrintError ("%s: empty input COMPLEX8TimeSeries!\n", __func__ );
      XLAL_ERROR (XLAL_EINVAL);
    }

  /* get timeseries epoch */
  deltat = (*x)->deltaT;

  /* loop over COMPLEX8TimeSeries elements */
  for ( k=0; k < (*x)->data->length; k++)
    {
      REAL8 tk = k * deltat;	/* time of k-th bin */
      REAL8 shiftCycles = shift * tk;
      REAL4 fact_re, fact_im;			/* complex phase-shift factor e^(-2pi f tau) */
      REAL4 yRe, yIm;

      /* use a sin/cos look-up-table for speed */
      sin_cos_2PI_LUT ( &fact_im, &fact_re, shiftCycles );

      /* apply the phase shift */
      yRe = fact_re * crealf((*x)->data->data[k]) - fact_im * cimagf((*x)->data->data[k]);
      yIm = fact_re * cimagf((*x)->data->data[k]) + fact_im * crealf((*x)->data->data[k]);
      (*x)->data->data[k] = crectf( yRe, yIm );

    } /* for k < numBins */

  /* adjust timeseries heterodyne frequency to the shift */
  (*x)->f0 -= shift;

  return XLAL_SUCCESS;

} /* XLALFrequencyShiftCOMPLEX8TimeSeries() */

/**
 * Apply a spin-down correction to the Fa and Fb complex timeseries
 * using the time-domain expression y(t) = x(t) * e^(-i 2pi sum f_k * (t-tref)^(k+1)),
 *
 * NOTE: this <b>modifies</b> the COMPLEX8TimeSeries Fa and Fb in place
 */
int
XLALSpinDownCorrectionMultiFaFb ( MultiCOMPLEX8TimeSeries **Fa,	                /**< [in/out] timeseries to time-shift */
				  MultiCOMPLEX8TimeSeries **Fb,	                /**< [in/out] timeseries to time-shift */
				  const PulsarDopplerParams *doppler		/**< parameter-space point to correct for */
				  )
{
  INT4 nspins = PULSAR_MAX_SPINS - 1;
  LIGOTimeGPS *epoch;
  UINT4 numSamples,numDetectors;
  REAL8 deltaref,deltaT;
  UINT4 i,j,k;

  /* sanity checks */
  if ( !(*Fa) || !(*Fa)->data || !(*Fb) || !(*Fb)->data )
    {
      XLALPrintError ("%s: empty input MultiCOMPLEX8TimeSeries!\n", __func__ );
      XLAL_ERROR (XLAL_EINVAL);
    }

  /* check if Fa and Fb have the same timeseries parameters */
  epoch = &(*Fa)->data[0]->epoch;
  numSamples = (*Fa)->data[0]->data->length;
  numDetectors = (*Fa)->length;
  deltaT = (*Fa)->data[0]->deltaT;
  if ( ( (*Fa)->length != numDetectors ) || ( (*Fb)->length != numDetectors ) )
    {
      XLALPrintError ("%s: Different numbers of detectors within the Fa and Fb MultiCOMPLEX8TimeSeries!\n", __func__ );
      XLAL_ERROR (XLAL_EINVAL);
    }
  for (i=0;i<(*Fa)->length;i++)
    {
      if ( ( XLALGPSCmp(epoch,&((*Fa)->data[i]->epoch)) != 0 ) || ( XLALGPSCmp(epoch,&((*Fb)->data[i]->epoch)) != 0 ) )
	{
	  XLALPrintError ("%s: Different start times for within Fa and Fb MultiCOMPLEX8TimeSeries!\n", __func__ );
	  XLAL_ERROR (XLAL_EINVAL);
	}
      if ( ( (*Fa)->data[i]->data->length != numSamples ) || ( (*Fb)->data[i]->data->length != numSamples ) )
	{
	  XLALPrintError ("%s: Different MultiCOMPLEX8TimeSeries lengths within Fa and Fb!\n", __func__ );
	  XLAL_ERROR (XLAL_EINVAL);
	}
      if ( ( (*Fa)->data[i]->deltaT != deltaT ) || ( (*Fb)->data[i]->deltaT != deltaT ) )
	{
	  XLALPrintError ("%s: Different MultiCOMPLEX8TimeSeries deltaT within Fa and Fb!\n", __func__ );
	  XLAL_ERROR (XLAL_EINVAL);
	}
    }

  /* determine number of spin down's and check if sensible */
  while (doppler->fkdot[nspins]==0.0) nspins--;
  if ( ( nspins < 0 ) || ( nspins > PULSAR_MAX_SPINS - 1 ) ) {
    XLALPrintError ("%s: Invalid number of spin derivitives, nspins = %d!\n", __func__,nspins );
    XLAL_ERROR (XLAL_EINVAL);
  }
 /*  printf("number of spin downs = %d\n",nspins); */

  /* compute the time difference between timeseries epoch and reference time */
  deltaref = XLALGPSGetREAL8(epoch) - XLALGPSGetREAL8(&(doppler->refTime));
 /*  printf("deltaref = %6.12f\n",deltaref); */
/*   printf("f0 = %6.12f\n",doppler->fkdot[0]); */
/*   printf("f1 = %6.12e\n",doppler->fkdot[1]); */

  /* apply spin derivitive correction to resampled timeseries */
  /* loop over spin derivitives (nspins = 1 means first derivitive, = 2 means second derivitive etc.. ) */
  for (j=1;j<=(UINT4)nspins;j++) {

    /* loop over time samples  and compute the spin down phase correction */
    for (k=0;k<numSamples;k++) {

      /* compute fractional number of cycles the spin-derivitive has added since the reftime */
      REAL8 cycles = fmod ( inv_fact[j+1]*doppler->fkdot[j]*pow(deltaref + k*deltaT,(REAL8)(j+1)), 1);
      REAL4 cosphase, sinphase;

      /* use look-up-table for speed to compute real and imaginary phase */
      sin_cos_2PI_LUT (&sinphase, &cosphase, -cycles );

      /* loop over detectors */
      for (i=0;i<numDetectors;i++) {

	/* apply phase correction to Fa and Fb */
	REAL8 Fare = crealf((*Fa)->data[i]->data->data[k])*cosphase - cimagf((*Fa)->data[i]->data->data[k])*sinphase;
	REAL8 Faim = cimagf((*Fa)->data[i]->data->data[k])*cosphase + crealf((*Fa)->data[i]->data->data[k])*sinphase;
	REAL8 Fbre = crealf((*Fb)->data[i]->data->data[k])*cosphase - cimagf((*Fb)->data[i]->data->data[k])*sinphase;
	REAL8 Fbim = cimagf((*Fb)->data[i]->data->data[k])*cosphase + crealf((*Fb)->data[i]->data->data[k])*sinphase;

	(*Fa)->data[i]->data->data[k] = crectf( Fare, Faim );
	(*Fb)->data[i]->data->data[k] = crectf( Fbre, Fbim );

      } /* (i<numDetectors) */

    } /* (k<numSamples) */

  } /* (j<nspins) */

  return XLAL_SUCCESS;

 } /* XLALSpinDownCorrectionMultiFaFb */



/* ===== Object creation/destruction functions ===== */

/**
 * Destroy a MultiCOMPLEX8TimeSeries structure.
 * Note, this is "NULL-robust" in the sense that it will not crash
 * on NULL-entries anywhere in this struct, so it can be used
 * for failure-cleanup even on incomplete structs
 */
void
XLALDestroyMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries *multiTimes )
{
  UINT4 X;
  COMPLEX8TimeSeries *tmp;

  if ( ! multiTimes ) {
    return;
  }
  if ( multiTimes->data ) {

    for ( X=0; X < multiTimes->length; X ++ ) {

      if ( (tmp = multiTimes->data[X]) != NULL ) {
	XLALDestroyCOMPLEX8TimeSeries(tmp);
      }

    }
    LALFree ( multiTimes->data );
  }
  LALFree ( multiTimes );

  return;

} /* XLALDestroyMultiCOMPLEX8TimeSeries() */

/**
 * Duplicates a MultiCOMPLEX8TimeSeries structure.
 * Allocates memory and copies contents.
 */
MultiCOMPLEX8TimeSeries *
XLALDuplicateMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries *multiTimes )
{

  if ( ! multiTimes )
    XLAL_ERROR_NULL ( XLAL_EINVAL, "Invalid NULL input for multi-timeseries 'multiTimes'\n" );

  if ( multiTimes->length == 0 || multiTimes->data == NULL )
    XLAL_ERROR_NULL ( XLAL_EINVAL, "Invalid empty input timeseries, 0 detectors or data==NULL\n");


  UINT4 numDet = multiTimes->length;

  // ----- prepare memory for multicomplex8timeseries container
  MultiCOMPLEX8TimeSeries *out;
  if ( ( out = XLALCalloc ( 1, sizeof(*out) ) ) == NULL )
    XLAL_ERROR_NULL ( XLAL_ENOMEM, "Failed to calloc ( 1, %d )\n", sizeof(*out) );
  out->length = numDet;
  if ( ( out->data = XLALCalloc ( numDet, sizeof(*out->data) ) ) == NULL )
    XLAL_ERROR_NULL ( XLAL_ENOMEM, "Failed to calloc ( %d, %d )\n", numDet, sizeof(*out->data) );

  // ----- copy each of the numDet complex8timeseries contents
  for ( UINT4 X = 0; X < numDet; X ++ )
    {
      if ( (out->data[X] = XLALDuplicateCOMPLEX8TimeSeries ( multiTimes->data[X] )) == NULL )
        XLAL_ERROR_NULL ( XLAL_EFUNC );
    }

  return out;

} /* XLALDuplicateMultiCOMPLEX8TimeSeries() */

/**
 * Duplicates a COMPLEX8TimeSeries structure.
 * Allocates memory and copies contents.
 */
COMPLEX8TimeSeries *
XLALDuplicateCOMPLEX8TimeSeries ( COMPLEX8TimeSeries *times )
{
  if ( !times )
    XLAL_ERROR_NULL ( XLAL_EINVAL, "Invalid NULL input for timeseries 'times'\n" );

  if ( times->data == NULL || times->data->length == 0 || times->data->data == NULL )
    XLAL_ERROR_NULL ( XLAL_EINVAL, "Invalid empty input timeseries, 0 bins or data==NULL\n");

  COMPLEX8TimeSeries *out;
  if ( (out = XLALCalloc ( 1, sizeof(*out) ) ) == NULL )
    XLAL_ERROR_NULL ( XLAL_ENOMEM, "Failed to calloc ( 1, %d )\n", sizeof(*out) );

  // copy header info [including data-pointer, will be reset]
  memcpy ( out, times, sizeof(*times) );

  UINT4 numBins = times->data->length;
  if ( ( out->data = XLALCreateCOMPLEX8Vector ( numBins )) == NULL )
    XLAL_ERROR_NULL ( XLAL_EFUNC, "XLALCreateCOMPLEX8Vector( %d ) failed.\n", numBins );

  // copy contents of COMPLEX8 vector
  memcpy ( out->data->data, times->data->data, numBins * sizeof(*times->data->data) );

  return out;

} /* XLALDuplicateCOMPLEX8TimeSeries() */

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
