/*
 * Copyright (C) 2009 Chris Messenger, Pinkesh Patel
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

/** \author R. Prix, J. T. Whelan
 * \ingroup pulsarCoherent
 * \file
 * \brief
 * Functions to calculate the so-called F-statistic for a given point in parameter-space,
 * following the equations in \ref JKS98.
 *
 * This code is partly a descendant of an earlier implementation found in
 * LALDemod.[ch] by Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve Berukoff, Xavier Siemens, Bruce Allen
 * ComputSky.[ch] by Jolien Creighton, Reinhard Prix, Steve Berukoff
 * LALComputeAM.[ch] by Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve Berukoff, Xavier Siemens
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
#include "../../FDS_isolated/Fstat_v3.h"
#include <lal/ComplexAM.h>
#include <lal/TimeSeries.h>

NRCSID( COMPUTEFSTATRSC, "$Id$");

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

/* empty initializers  */
static const LALStatus empty_status;
static const AMCoeffs empty_AMCoeffs;

const SSBtimes empty_SSBtimes;
const MultiSSBtimes empty_MultiSSBtimes;
const AntennaPatternMatrix empty_AntennaPatternMatrix;
const MultiAMCoeffs empty_MultiAMCoeffs;
const Fcomponents empty_Fcomponents;
const ComputeFParams empty_ComputeFParams;
const ComputeFBuffer empty_ComputeFBuffer;
const EarthState empty_EarthState;

/*---------- internal prototypes ----------*/
int finite(double x);

/*==================== FUNCTION DEFINITIONS ====================*/


/** Function to compute a vector of Fstatistic values for a number of frequency bins.
    The output, i.e. fstatVector must be properly allocated
    before this function is called.  The values of the start frequency, the step size
    in the frequency and the number of frequency values for which the Fstatistic is
    to be calculated are read from fstatVector.  The other parameters are not checked and
    they must be correctly set outside this function.
*/ 
void ComputeFStatFreqBand_RS ( LALStatus *status,
			       REAL4FrequencySeries *fstatVector, 		/**< [out] Vector of Fstat values */
			       const PulsarDopplerParams *doppler,		/**< parameter-space point to compute F for */
			       MultiSFTVector *multiSFTs, 		/**< normalized (by DOUBLE-sided Sn!) data-SFTs of all IFOs */
			       const MultiNoiseWeights *multiWeights,	/**< noise-weights of all SFTs */
			       MultiDetectorStateSeries *multiDetStates,/**< 'trajectories' of the different IFOs */
			       ComputeFParams *params		/**< addition computational params */
			       )
{
/*   static const CHAR *fn = "ComputeFStatFreqBand_RS()"; */

  UINT4 numDetectors; 
  ComputeFBuffer *cfBuffer = NULL;
  MultiSSBtimes *multiSSB = NULL;
  MultiAMCoeffs *multiAMcoef = NULL;
  MultiCOMPLEX8TimeSeries *multiTimeseries = NULL;
  REAL8 Ad, Bd, Cd, Dd_inv;
  SkyPosition skypos;
  UINT4 i,j,k;
  MultiCOMPLEX8TimeSeries *multiFa_resampled = NULL;
  MultiCOMPLEX8TimeSeries *multiFb_resampled = NULL;
  COMPLEX8Vector *Faf_resampled = NULL;
  COMPLEX8Vector *Fbf_resampled = NULL;
  UINT4 numSamples;
  REAL8 Tsft;
  REAL8 f0_sft;
  REAL8 dt;
  SFTtype *firstSFT;
  REAL8 df_out;
  REAL8 f0_out;
  ComplexFFTPlan *pfwd = NULL;                                                /* this will store the FFT plan */
  COMPLEX8Vector *outa = NULL;                                                /* this will contain the FFT output of Fa for this detector */
  COMPLEX8Vector *outb = NULL;                                                /* this will contain the FFT output of Fb for this detector */

  INITSTATUS( status, "ComputeFStatFreqBand_RS", COMPUTEFSTATRSC );
  ATTATCHSTATUSPTR (status);

  /* check that the input data and parameters structures don't point to NULL */
  ASSERT ( multiSFTs, status, COMPUTEFSTATRSC_ENULL, COMPUTEFSTATRSC_MSGENULL );
  ASSERT ( doppler, status, COMPUTEFSTATRSC_ENULL, COMPUTEFSTATRSC_MSGENULL );
  ASSERT ( multiDetStates, status, COMPUTEFSTATRSC_ENULL, COMPUTEFSTATRSC_MSGENULL );
  ASSERT ( params, status, COMPUTEFSTATRSC_ENULL, COMPUTEFSTATRSC_MSGENULL );

  cfBuffer = params->buffer;                      /* set local pointer to the buffer location */
  numDetectors = multiSFTs->length;               /* set the number of detectors to the number of sets of SFTs */
  firstSFT = &(multiSFTs->data[0]->data[0]);      /* use data from the first SFT from the first detector to set other params */
  Tsft = 1.0 / firstSFT->deltaF;                  /* define the length of an SFT (assuming 1/T resolution) */
  f0_sft = firstSFT->f0;                          /* define the frequency of the first bin in the SFT */  

  /* check that the pre-allocated output vector doesn't point to NULL */ 
  ASSERT ( multiDetStates->length == numDetectors, status, COMPUTEFSTATRSC_EINPUT, COMPUTEFSTATRSC_MSGEINPUT );
  ASSERT ( fstatVector, status, COMPUTEFSTATRSC_ENULL, COMPUTEFSTATRSC_MSGENULL );
  ASSERT ( fstatVector->data, status, COMPUTEFSTATRSC_ENULL, COMPUTEFSTATRSC_MSGENULL );
  ASSERT ( fstatVector->data->data, status, COMPUTEFSTATRSC_ENULL, COMPUTEFSTATRSC_MSGENULL );
  ASSERT ( fstatVector->data->length > 0, status, COMPUTEFSTATRSC_EINPUT, COMPUTEFSTATRSC_MSGEINPUT );
 
  df_out = fstatVector->deltaF;                   /* the user defined frequency resolution */
  f0_out = fstatVector->f0;                       /* the user defined first output frequency bin */    

  /* check that the multidetector states have the same length as the multiSFTs */
  ASSERT ( multiDetStates->length == numDetectors, status, COMPUTEFSTATRSC_EINPUT, COMPUTEFSTATRSC_MSGEINPUT );
  if ( multiWeights ) {
    ASSERT ( multiWeights->length == numDetectors , status, COMPUTEFSTATRSC_EINPUT, COMPUTEFSTATRSC_MSGEINPUT );
  }

  /* generate bandpassed and downsampled timeseries for each detector                         */
  /* we only ever do this once for a given dataset so we read it from the buffer if it exists */
  /* in future implementations we will pass this directly to the function instead of SFTs     */
  if (cfBuffer) {
    if (cfBuffer->multiTimeseries) multiTimeseries = cfBuffer->multiTimeseries;
    else {
      LALPrintError("\nComputeFstat buffer non-null but multiTimeseries is NULL, strange!\n\n");
      ABORT ( status, COMPUTEFSTATRSC_ENULL, COMPUTEFSTATRSC_MSGENULL );
    }
  }
  else {
    
    /* generate multiple coincident timeseries - one for each detector spanning start -> end */ 
    /* we need each timeseries to span the exact same amount of time and to start at the same time */
    /* because for the multi-detector Fstat we need frequency bins to be coincident */
    if ( (multiTimeseries = XLALMultiSFTVectorToCOMPLEX8TimeSeries_CHRIS(multiSFTs)) == NULL ) {
      LALPrintError("\nXLALMultiSFTVectorToCOMPLEX8TimeSeries() failed with error = %d\n\n", xlalErrno );
      ABORT ( status, COMPUTEFSTATRSC_EXLAL, COMPUTEFSTATRSC_MSGEXLAL );
    }
  
    /* compute the fractional bin offset between the user requested initial frequency */
    /* and the closest output frequency bin */
    {
      REAL8 diff = multiTimeseries->data[0]->f0 - fstatVector->f0;         /* the difference between the new timeseries heterodyne frequency and the user requested frequency */
      UINT4 bins = (UINT4)floor(0.5 + diff/fstatVector->deltaF);           /* the rounded number of output frequency bins difference */ 
      REAL8 shift = diff - fstatVector->deltaF*bins;                       /* the fractional bin frequency offset */
      printf("Shift = %6.12f Hz\n",shift);
      printf("time series heterodyne freq = %6.12f\n",multiTimeseries->data[0]->f0);
      printf("requested first frequency = %6.12f\n",fstatVector->f0);
      
      /* shift the timeseries by a fraction of a frequency bin so that user requested frequency is exactly resolved */ 
      if (shift != 0.0) {
	if ( (XLALFrequencyShiftMultiCOMPLEX8TimeSeries(&multiTimeseries,shift)) != XLAL_SUCCESS ) {
	  LALPrintError("\nXLALMultiSFTVectorToCOMPLEX8TimeSeries() failed with error = %d\n\n", xlalErrno );
	  ABORT ( status, COMPUTEFSTATRSC_EXLAL, COMPUTEFSTATRSC_MSGEXLAL );
	}
      }
    }
    
    /* recompute the multidetector states for the possibly time shifted SFTs */
    /* the function XLALMultiSFTVectorToCOMPLEX8TimeSeries may have shifted the SFT start times around */
    /* and since these times will be used later on for the resampling we also need to recompute the */
    /* MultiDetectorStates because their timestamps are used later on to compute SSB times which we */
    /* need to be accurate at the midpoints of the SFTs.  Understand ? */
    XLALDestroyMultiDetectorStateSeries(multiDetStates);                                          /* destroy the input multiDetStates */
    multiDetStates = NULL;                                                                        /* and set to NULL */
    LALGetMultiDetectorStates ( status->statusPtr, &multiDetStates, multiSFTs, params->edat );    /* recompute the multiDetStates for the new SFT start times */

    /* point the buffer to the newly generated timeseries and set all other buffer pointers to NULL */
    /* if we're in this if statement then we will always need to allocate mem for the buffer pointers */
    if ( (cfBuffer = XLALMalloc(sizeof(ComputeFBuffer))) == NULL ) {
      LALPrintError("\nXLALMalloc() failed with error = %d\n\n", xlalErrno );
      ABORT ( status, COMPUTEFSTATRSC_EXLAL, COMPUTEFSTATRSC_MSGEXLAL );
    }			  
    cfBuffer->multiTimeseries = multiTimeseries;
    cfBuffer->multiSSB = NULL;
    cfBuffer->multiBinary = NULL;
    cfBuffer->multiAMcoef = NULL;
    cfBuffer->multiCmplxAMcoef = NULL; 
    cfBuffer->multiFa_resampled = NULL;                
    cfBuffer->multiFb_resampled = NULL;         

  }  /* if (cfBuffer->multiTimeseries) */
 
  /* if the sky position has changed or if any of the sky position dependent quantities are not buffered */
  /* i.e the multiDetstates, the multiAMcoefficients, the multiSSB times and the resampled multiTimeSeries Fa and Fb */
  /* then we need to recompute these and buffer them */
  if ( cfBuffer                                                     /* if we have a buffer */
       && ( cfBuffer->Alpha == doppler->Alpha )                     /* and alpha hasn't changed */
       && ( cfBuffer->Delta == doppler->Delta )                     /* and delta hasn't changed */
       && ( cfBuffer->multiDetStates == multiDetStates )            /* and the buffered multiDetStates hasn't changed */
       && ( cfBuffer->multiAMcoef )                                 /* and we have a buffered multiAMcoefficents */
       && ( cfBuffer->multiSSB )                                    /* and we have buffered multiSSB times */ 
       && ( cfBuffer->multiFa_resampled )                           /* and we have buffered multiFa_resampled  */
       && ( cfBuffer->multiFb_resampled )                           /* and we have multiFb_resampled */
       )
    { /* if already buffered then we reuse all sky dependent buffered values */
         
      multiSSB = cfBuffer->multiSSB;                                /* re-use the multiSSB times for the SFT midpoints */
      multiAMcoef = cfBuffer->multiAMcoef;                          /* re-use AM coefficients */
      multiFa_resampled = cfBuffer->multiFa_resampled;              /* re-use resampled multidetector Fa timeseries */
      multiFb_resampled = cfBuffer->multiFb_resampled;              /* re-use resampled multidetector Fb timeseries */
      
    } 
  /* otherwise we need to (re)compute all of the sky position dependent quantities */
  else {
  
    MultiCOMPLEX8TimeSeries *multiFa = NULL;
    MultiCOMPLEX8TimeSeries *multiFb = NULL;
    
    /* compute the SSB times corresponding to the midpoints of each SFT for the current sky position for all detectors */
    skypos.system =   COORDINATESYSTEM_EQUATORIAL;
    skypos.longitude = doppler->Alpha;
    skypos.latitude  = doppler->Delta;
    TRY ( LALGetMultiSSBtimes ( status->statusPtr, &multiSSB, multiDetStates, skypos, doppler->refTime, params->SSBprec ), status );
  
    /* compute the AM parameters for each detector */
    LALGetMultiAMCoeffs ( status->statusPtr, &multiAMcoef, multiDetStates, skypos );
    BEGINFAIL ( status ) {
      XLALDestroyMultiSSBtimes ( multiSSB );
    } ENDFAIL (status);
    
    /* noise-weight Antenna-patterns and compute A,B,C */
    if ( XLALWeighMultiAMCoeffs ( multiAMcoef, multiWeights ) != XLAL_SUCCESS ) {
      LALPrintError("\nXLALWeighMultiAMCoeffs() failed with error = %d\n\n", xlalErrno );
      ABORT ( status, COMPUTEFSTATRSC_EXLAL, COMPUTEFSTATRSC_MSGEXLAL );
    }
 
    /* Generate a(t) and b(t) weighted heterodyned downsampled timeseries */
    if ( XLALAntennaWeightMultiCOMPLEX8TimeSeries(&multiFa,&multiFb,multiTimeseries,multiAMcoef,multiSFTs) != XLAL_SUCCESS ) {
      LALPrintError("\nXLALAntennaWeightMultiCOMPLEX8TimeSeries() failed with error = %d\n\n", xlalErrno );
      ABORT ( status, COMPUTEFSTATRSC_EXLAL, COMPUTEFSTATRSC_MSGEXLAL );
    }
    printf("done XLALAntennaWeightMultiCOMPLEX8TimeSeries()\n");
            
    /* Perform barycentric resampling on the multi-detector timeseries */
    if ( XLALBarycentricResampleMultiCOMPLEX8TimeSeries(&multiFa_resampled,&multiFb_resampled,multiFa,multiFb,multiSSB,multiSFTs,Tsft,df_out) != XLAL_SUCCESS ) {
      LALPrintError("\nXLALBarycentricResampleMultiCOMPLEX8TimeSeries() failed with error = %d\n\n", xlalErrno );
      ABORT ( status, COMPUTEFSTATRSC_EXLAL, COMPUTEFSTATRSC_MSGEXLAL );
    }
    printf("done XLALBarycentricResampleMultiCOMPLEX8TimeSeries()\n");

    /* now buffer all of this for the next call to this function in case it uses the same sky position */
    /* first we destroy all of the previous buffered data */
    if ( cfBuffer->multiSSB ) XLALDestroyMultiSSBtimes ( cfBuffer->multiSSB );
    if ( cfBuffer->multiAMcoef ) XLALDestroyMultiAMCoeffs ( cfBuffer->multiAMcoef );
    if ( cfBuffer->multiFa_resampled ) XLALDestroyMultiCOMPLEX8TimeSeries ( cfBuffer->multiFa_resampled );
    if ( cfBuffer->multiFb_resampled ) XLALDestroyMultiCOMPLEX8TimeSeries ( cfBuffer->multiFb_resampled );
    /* then we buffer the new data */
    cfBuffer->multiSSB = multiSSB;
    cfBuffer->Alpha = doppler->Alpha;
    cfBuffer->Delta = doppler->Delta;
    cfBuffer->multiDetStates = multiDetStates;
    cfBuffer->multiAMcoef = multiAMcoef;
    cfBuffer->multiFa_resampled = multiFa_resampled;
    cfBuffer->multiFb_resampled = multiFb_resampled;
    
  } /* could not reuse previously buffered quantities */

  /* store AM coefficient integrals in local variables */
  if ( multiAMcoef ) {
    Ad = multiAMcoef->Mmunu.Ad;
    Bd = multiAMcoef->Mmunu.Bd;
    Cd = multiAMcoef->Mmunu.Cd;
    Dd_inv = 1.0 / multiAMcoef->Mmunu.Dd;
  }
  else {
    LALPrintError ( "Programming error: 'multiAMcoef' not available!\n");
    ABORT ( status, COMPUTEFSTATRSC_ENULL, COMPUTEFSTATRSC_MSGENULL );
  }
    
  /* apply spin derivitive correction to resampled timeseries */
  /* this function only applies a correction if there are any non-zero spin derivitives */
  if ( XLALSpinDownCorrectionMultiFaFb(&multiFa_resampled,&multiFb_resampled,doppler) != XLAL_SUCCESS ) {
    LALPrintError("\nXLALSpinDownCorrectionMultiFaFb() failed with error = %d\n\n", xlalErrno );
    ABORT ( status, COMPUTEFSTATRSC_EXLAL, COMPUTEFSTATRSC_MSGEXLAL );
  }

  /* we now compute the FFTs of the resampled functions Fa and Fb for each detector */
  /* and combine them into the F-statistic */

  /* we use the first detector Fa time series to obtain the number of time samples and the sampling time */
  /* these should be the same for all Fa and Fb timeseries */
  numSamples = multiFa_resampled->data[0]->data->length;
  dt = multiFa_resampled->data[0]->deltaT;
  
  /* allocate memory for Fa(f) and Fb(f) and individual detector FFT outputs */
  if ( (Faf_resampled = XLALCreateCOMPLEX8Vector(numSamples)) == NULL )
    ABORT ( status, xlalErrno, "XLALCreateCOMPLEX8Vector() failed!\n");
  if ( (Fbf_resampled = XLALCreateCOMPLEX8Vector(numSamples)) == NULL ) 
    ABORT ( status, xlalErrno, "XLALCreateCOMPLEX8Vector() failed!\n");
  if ( (outa = XLALCreateCOMPLEX8Vector(numSamples)) == NULL ) 
    ABORT ( status, xlalErrno, "XLALCreateCOMPLEX8Vector() failed!\n");
  if ( (outb = XLALCreateCOMPLEX8Vector(numSamples)) == NULL ) 
    ABORT ( status, xlalErrno, "XLALCOMPLEX8VectorFFT() failed!\n");
 
  /* initialise output vectors to zero since it will be added to */
  memset(Faf_resampled->data,0,numSamples*sizeof(COMPLEX8));
  memset(Fbf_resampled->data,0,numSamples*sizeof(COMPLEX8));
 
   /* make forwards FFT plan - this will be re-used for each detector */
  if ( (pfwd = XLALCreateCOMPLEX8FFTPlan(numSamples,1,0)) == NULL ) 
    ABORT ( status, xlalErrno, "XLALCreateCOMPLEX8FFTPlan() failed!\n");
 
  /* loop over detectors */
  for (i=0;i<numDetectors;i++) {
  
    COMPLEX8Vector *ina = (COMPLEX8Vector*)multiFa_resampled->data[i]->data;    /* we point the input to the current detector Fa timeseries */
    COMPLEX8Vector *inb = (COMPLEX8Vector*)multiFb_resampled->data[i]->data;    /* we point the input to the current detector Fb timeseries */

    /* printf("start = %d %d reftime = %d %d\n",multiFa_resampled->data[0]->epoch.gpsSeconds,multiFa_resampled->data[0]->epoch.gpsNanoSeconds,doppler->refTime.gpsSeconds,doppler->refTime.gpsNanoSeconds); */
/*     printf("fdot = %6.12e\n",doppler->fkdot[1]); */
/*     printf("SSBstart = %d %d\n",multiFa_resampled->data[0]->epoch.gpsSeconds,multiFa_resampled->data[0]->epoch.gpsNanoSeconds); */
/*     printf("SSBrefTime = %d %d\n",doppler->refTime.gpsSeconds,doppler->refTime.gpsNanoSeconds); */
/*     printf("numTimeSamples = %d\n",numSamples); */

    /* initialise output vectors to zero for safety */
    memset(outa->data,0,numSamples*sizeof(COMPLEX8));
    memset(outb->data,0,numSamples*sizeof(COMPLEX8));
   
    /* Fourier transform the resampled Fa(t) and Fb(t) */
    if (XLALCOMPLEX8VectorFFT(outa,ina,pfwd)!= XLAL_SUCCESS)
      ABORT ( status, xlalErrno, "XLALCOMPLEX8VectorFFT() failed!\n");
    if (XLALCOMPLEX8VectorFFT(outb,inb,pfwd)!= XLAL_SUCCESS)
      ABORT ( status, xlalErrno, "XLALCOMPLEX8VectorFFT() failed!\n");
  
    /*  add to summed Faf and Fbf */
    for (j=0;j<numSamples;j++) {
      Faf_resampled->data[j].re += outa->data[j].re*dt;
      Faf_resampled->data[j].im += outa->data[j].im*dt;
      Fbf_resampled->data[j].re += outb->data[j].re*dt;
      Fbf_resampled->data[j].im += outb->data[j].im*dt;
    }
  
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
  f0_sft = multiFa_resampled->data[0]->f0 - NhalfNeg(numSamples)*df_out;
  
  /* loop over frequency and construct 2F */
  {
    UINT4 offset = floor(0.5 + (fstatVector->f0 - f0_sft)/fstatVector->deltaF);
    for (k=0;k<fstatVector->data->length;k++) {
      
      UINT4 idx = k + offset;
      /* ----- compute final Fstatistic-value ----- */
      
      /* NOTE: the data MUST be normalized by the DOUBLE-SIDED PSD (using LALNormalizeMultiSFTVect),
       * therefore there is a factor of 2 difference with respect to the equations in JKS, which
       * where based on the single-sided PSD.
       */
      fstatVector->data->data[k] = Dd_inv * (  Bd * (SQ(Faf_resampled->data[idx].re) + SQ(Faf_resampled->data[idx].im) )
					       + Ad * ( SQ(Fbf_resampled->data[idx].re) + SQ(Fbf_resampled->data[idx].im) )
					       - 2.0 * Cd *( Faf_resampled->data[idx].re * Fbf_resampled->data[idx].re + 
							     Faf_resampled->data[idx].im * Fbf_resampled->data[idx].im )
					       );
      
    }
  } 

  /* free memory if no buffer was available */
  if ( !cfBuffer ) {
    XLALDestroyMultiSSBtimes ( multiSSB );
    XLALDestroyMultiAMCoeffs ( multiAMcoef );
  } /* if !cfBuffer */
  
  
  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* ComputeFStatFreqBand_RS() */


/** Turn the given multiSFTvector into multiple long COMPLEX8TimeSeries, properly dealing with gaps.
 * Memeory allocation for the output MultiCOMPLEX8TimeSeries is done within this function.
 *
 * NOTE : We enforce that each detectors timeseries has <b>equal</b> start times and time spans.
 * Also, the input MultiSFTs get <b>modified</b> in place.
 */
MultiCOMPLEX8TimeSeries *XLALMultiSFTVectorToCOMPLEX8TimeSeries_CHRIS ( MultiSFTVector *multisfts  /**< [in/out] multi SFT vector, gets modified! */							
									)	
{
  static const CHAR *fn = "XLALMultiSFTVectorToCOMPLEX8TimeSeries_CHRIS()";

  UINT4 i;
  MultiCOMPLEX8TimeSeries *out = NULL;		/* long time-series corresponding to full set of SFTs */
  LIGOTimeGPS start,end;

  /* check sanity of input */
  if (!multisfts || (multisfts->length == 0)) {
    XLALPrintError ("%s: empty multi SFT input!\n", fn );
    XLAL_ERROR_NULL (fn, XLAL_EINVAL);
  }
  for (i=0;i<multisfts->length;i++)
    {
      if (!multisfts->data[i] || (multisfts->data[i]->length == 0)) {
	XLALPrintError ("%s: empty multiSFT->data[%d] input!\n", fn,i );
	XLAL_ERROR_NULL (fn, XLAL_EINVAL);
      }
    }

  /* determine the start and end times of the multiSFT observation */
  if ( XLALEarliestMultiSFTsample(&start,multisfts) != XLAL_SUCCESS ) {
    XLALPrintError("%s: Failed to run XLALEarliestMultiSFTsample()\n", fn );
    XLAL_ERROR_NULL (fn, XLAL_EFAULT );
  }
  if ( XLALLatestMultiSFTsample(&end,multisfts) != XLAL_SUCCESS ) {
    XLALPrintError("%s: Failed to run XLALLatestMultiSFTsample()\n", fn );
    XLAL_ERROR_NULL (fn, XLAL_EFAULT );
  }
  printf("earliest sample at %d %d\n",start.gpsSeconds,start.gpsNanoSeconds);
  printf("latest sample at %d %d\n",end.gpsSeconds,end.gpsNanoSeconds);

  /* check that earliest is before latest */
  if ( (XLALGPSDiff ( &end, &start ) ) < 0 ) 
    {
      XLALPrintError ("%s: start time after end time!\n", fn );
      XLAL_ERROR_NULL (fn, XLAL_EINVAL);
    }

  /* allocate memory for the output structure */
  if ((out = XLALMalloc(sizeof(MultiCOMPLEX8TimeSeries))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", fn, sizeof(MultiCOMPLEX8TimeSeries));
    XLAL_ERROR_NULL (fn, XLAL_ENOMEM);
  }
  out->length = multisfts->length;
  if ((out->data = XLALMalloc(out->length*sizeof(COMPLEX8TimeSeries*))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", fn, out->length*sizeof(COMPLEX8TimeSeries*));
    XLAL_ERROR_NULL (fn, XLAL_ENOMEM);
  }

  /* loop over detectors */
  for (i=0;i<multisfts->length;i++) {
    
    /* call XLALSFTVectorToCOMPLEX8TimeSeries for each detector */
    if ((out->data[i] = XLALSFTVectorToCOMPLEX8TimeSeries_CHRIS(multisfts->data[i],&start,&end)) == NULL) {
      XLALPrintError ("%s: Failed to run XLALSFTVectorToCOMPLEX8TimeSeries()\n", fn);
      XLAL_ERROR_NULL (fn, XLAL_EFAULT );
    }
    
  }

  return out;
  
} /* XLALMultiSFTVectorToCOMPLEX8TimeSeries() */

/** Turn the given SFTvector into one long time-series, properly dealing with gaps.
 *
 * NOTE: this function <b>modifies</b> the input SFTs in the process!
 * If you need to reuse the SFTvector afterwards, you need to copy it before
 * passing it into this function.  Also, if start and end are NULL then the 
 * start and end of the SFT vector is used.
 *
 */
COMPLEX8TimeSeries *XLALSFTVectorToCOMPLEX8TimeSeries_CHRIS ( SFTVector *sfts,             /**< [in/out] SFT vector, gets modified! */
							      const LIGOTimeGPS *start_in,    /**< [in] start time */
							      const LIGOTimeGPS *end_in       /**< [in] input end time */
							      )
{
  static const CHAR *fn = "XLALSFTVectorToCOMPLEX8TimeSeries_CHRIS()";

  COMPLEX8FFTPlan *SFTplan;

  COMPLEX8TimeSeries *lTS = NULL;		/* long time-series corresponding to full set of SFTs */
  COMPLEX8TimeSeries *sTS = NULL; 		/* short time-series corresponding to a single SFT */

  REAL8 fHet;				/* heterodyning frequency */
  LIGOTimeGPS epoch = {0,0};
  LIGOTimeGPS start;
  LIGOTimeGPS end;
  
  /* constant quantities for all SFTs */
  SFTtype *firstSFT, *lastSFT;
  UINT4 numBinsSFT;
  REAL8 dfSFT;
  REAL8 Tsft;
  REAL8 deltaT;
  UINT4 numSFTs;
  UINT4 n;
  REAL8 Tspan;
  UINT4 numTimeSamples;
  UINT4 NnegSFT;
  REAL8 f0SFT;
  REAL8 SFTFreqBand;

  /* check sanity of input */
  if ( !sfts || (sfts->length == 0) )
    {
      XLALPrintError ("%s: empty SFT input!\n", fn );
      XLAL_ERROR_NULL (fn, XLAL_EINVAL);
    }

  /* define some useful shorthands */
  numSFTs = sfts->length;
  firstSFT = &(sfts->data[0]);
  lastSFT = &(sfts->data[numSFTs-1]);
  numBinsSFT = firstSFT->data->length;
  dfSFT = firstSFT->deltaF;
  Tsft = 1.0 / dfSFT;
  SFTFreqBand = numBinsSFT * dfSFT;
  deltaT = 1.0 / SFTFreqBand;	/* we'll put DC into the middle of [f0, f0+Band], so sampling at fSamp=Band is sufficient */
  f0SFT = firstSFT->f0;

  /* if the start and end input pointers are NOT NULL then determine start and time-span of the final long time-series */
  if (start_in && end_in) 
    {
      start.gpsSeconds = start_in->gpsSeconds;
      start.gpsNanoSeconds = start_in->gpsNanoSeconds;
      end.gpsSeconds = end_in->gpsSeconds;
      end.gpsNanoSeconds = end_in->gpsNanoSeconds;

      /* do sanity checks */
      if ( (XLALGPSDiff ( &end, &firstSFT->epoch ) ) < 0 ) 
	{
	  XLALPrintError ("%s: end time before first SFT!\n", fn );
	  XLAL_ERROR_NULL (fn, XLAL_EINVAL);
	}
      if ( (XLALGPSDiff ( &start, &sfts->data[numSFTs-1].epoch) ) > Tsft ) 
	{
	  XLALPrintError ("%s: start time after end of data!\n", fn );
	  XLAL_ERROR_NULL (fn, XLAL_EINVAL);
	}
    }
  else {   /* otherwise we use the start and end of the sft vector */
    start.gpsSeconds = firstSFT->epoch.gpsSeconds;
    start.gpsNanoSeconds = firstSFT->epoch.gpsNanoSeconds;
    end.gpsSeconds = lastSFT->epoch.gpsSeconds;
    end.gpsNanoSeconds = lastSFT->epoch.gpsNanoSeconds;
    if ( XLALGPSAdd(&end,Tsft) == NULL )
    {
      XLALPrintError ("%s: NULL pointer returned from XLALGPSAdd()!\n", fn );
      XLAL_ERROR_NULL (fn, XLAL_EFAULT);
    }
    
  }
 
  /* determine output time span */
  if ( (Tspan = XLALGPSDiff ( &end, &start ) ) < 0 ) 
    {
      XLALPrintError ("%s: start time after end time!\n", fn );
      XLAL_ERROR_NULL (fn, XLAL_EINVAL);
    }
  
  numTimeSamples = (UINT4)floor(Tspan / deltaT + 0.5);	/* round */
  printf("new Tspan = %6.12f\n",Tspan);
  printf("new end time = %6.12f\n",GPS2REAL8(start)+Tspan);

  /* determine the heterodyning frequency */
  /* fHet = DC of our internal DFTs */
  NnegSFT = NhalfNeg ( numBinsSFT );
  fHet = f0SFT + 1.0 * NnegSFT * dfSFT;
  printf("fHet = %6.12f\n",fHet);

  /* ----- Prepare invFFT of SFTs: compute plan for FFTW */
  if ( (SFTplan = XLALCreateReverseCOMPLEX8FFTPlan( numBinsSFT, 0 )) == NULL )
    {
      XLALPrintError ( "%s: XLALCreateReverseCOMPLEX8FFTPlan(%d, ESTIMATE ) failed! errno = %d!\n", fn, numBinsSFT, xlalErrno );
      goto failed;
    }

  /* ----- Prepare short time-series holding ONE invFFT of a single SFT */
  if ( (sTS = XLALCreateCOMPLEX8TimeSeries ( "short timeseries", &epoch, 0, deltaT, &empty_LALUnit, numBinsSFT )) == NULL )
    {
      XLALPrintError ( "%s: XLALCreateCOMPLEX8TimeSeries() for %d timesteps failed! errno = %d!\n", fn, numBinsSFT, xlalErrno );
      goto failed;
    }

  /* ----- prepare long TimeSeries container ---------- */
  if ( (lTS = XLALCreateCOMPLEX8TimeSeries ( firstSFT->name, &start, fHet, deltaT, &empty_LALUnit, numTimeSamples )) == NULL )
    {
      XLALPrintError ("%s: XLALCreateCOMPLEX8TimeSeries() for %d timesteps failed! errno = %d!\n", fn, numTimeSamples, xlalErrno );
      goto failed;
    }
  memset ( lTS->data->data, 0, numTimeSamples * sizeof(*lTS->data->data)); 	/* set all time-samples to zero (in case there are gaps) */


  /* ---------- loop over all SFTs and inverse-FFT them ---------- */
  for ( n = 0; n < numSFTs; n ++ )
    {
      SFTtype *thisSFT = &(sfts->data[n]);
      UINT4 bin0_n;
      REAL8 offset_n, nudge_n;
      UINT4 copyLen, binsLeft;

      REAL8 offset0, offsetEff, hetCycles;
      COMPLEX8 hetCorrection;

      /* find bin in long timeseries corresponding to starttime of *this* SFT */
      offset_n = XLALGPSDiff ( &thisSFT->epoch, &start );
      bin0_n = (UINT4) ( offset_n / deltaT + 0.5 );	/* round to closest bin */

      nudge_n = bin0_n * deltaT - offset_n;		/* rounding error */
      nudge_n = 1e-9 * (floor)(nudge_n * 1e9 + 0.5);	/* round to closest nanosecond */
      {
	REAL8 t0 = XLALGPSGetREAL8 ( &start );
	XLALPrintInfo ("n = %d: t0_n = %f, sft_tn =(%d,%d), bin-offset = %g s, corresponding to %g timesteps\n",
		n, t0 + bin0_n * deltaT, sfts->data[n].epoch.gpsSeconds,  sfts->data[n].epoch.gpsNanoSeconds, nudge_n, nudge_n/deltaT );
      }
      
     /*  printf("nudge_n = %6.12f\n",nudge_n); */
      /* nudge SFT into integer timestep bin if necessary */
      if ( nudge_n != 0 )
	{
	  if  ( XLALTimeShiftSFT ( thisSFT, nudge_n ) != XLAL_SUCCESS )
	    {
	      XLALPrintError ( "%s: XLALTimeShiftSFT(sft-%d, %g) failed! errno = %d!\n", fn, n, nudge_n, xlalErrno );
	      goto failed;
	    }
	}

      /* determine heterodyning phase-correction for this SFT */
      offset0 = XLALGPSDiff ( &thisSFT->epoch, &start );

      /* fHet * Tsft is an integer, because fHet is a frequency-bin of the input SFTs, so we only need the remainder offset_t0 % Tsft */
      offsetEff = fmod ( offset0, Tsft );
      offsetEff = 1e-9 * (floor)( offsetEff * 1e9 + 0.5 );	/* round to closest integer multiple of nanoseconds */
      hetCycles = fmod ( fHet * offsetEff, 1);			/* required heterodyning phase-correction for this SFT */
    
      sin_cos_2PI_LUT (&hetCorrection.im, &hetCorrection.re, -hetCycles );
     
      /* Note: we also bundle the overall normalization of 'df' into the het-correction.
       * This ensures that the resulting timeseries will have the correct normalization, according to
       * x_l = invFT[sft]_l = df * sum_{k=0}^{N-1} xt_k * e^(i 2pi k l / N )
       * where x_l is the l-th timestamp, and xt_k is the k-th frequency bin of the SFT.
       * See the LAL-conventions on FFTs:  http://www.ligo.caltech.edu/docs/T/T010095-00.pdf
       * (the FFTw convention does not contain the factor of 'df', which is why we need to
       * apply it ourselves)
       *
       */
      hetCorrection.re *= dfSFT;
      hetCorrection.im *= dfSFT;

      /* FIXME: check how time-critical this step is, using proper profiling! */
      if ( XLALMultiplySFTbyCOMPLEX8 ( thisSFT, hetCorrection ) != XLAL_SUCCESS )
	{
	  XLALPrintError ( "%s: XLALMultiplySFTbyCOMPLEX8(sft-%d) failed! errno = %d!\n", fn, n, xlalErrno );
	  goto failed;
	}

      XLALPrintInfo ("SFT n = %d: (tn - t0) = %g s EQUIV %g s, hetCycles = %g, ==> fact = %g + i %g\n",
                     n, offset0, offsetEff, hetCycles, hetCorrection.re, hetCorrection.im );

      /* FIXME: check if required */
      if ( XLALReorderSFTtoFFTW (thisSFT->data) != XLAL_SUCCESS )
	{
	  XLALPrintError ( "%s: XLALReorderSFTtoFFTW() failed! errno = %d!\n", fn, xlalErrno );
	  goto failed;
	}

      if ( XLALCOMPLEX8VectorFFT( sTS->data, thisSFT->data, SFTplan ) != XLAL_SUCCESS )
	{
	  XLALPrintError ( "%s: XLALCOMPLEX8VectorFFT() failed! errno = %d!\n", fn, xlalErrno );
	  goto failed;
	}

      /* copy short (shifted) timeseries into correct location within long timeseries */
      binsLeft = numTimeSamples - bin0_n - 1;
      copyLen = MYMIN ( numBinsSFT, binsLeft );		/* make sure not to write past the end of the long TS */
      memcpy ( &lTS->data->data[bin0_n], sTS->data->data, copyLen * sizeof(lTS->data->data[0]) );

    } /* for n < numSFTs */

  goto success;

 failed:
  /* cleanup memory */
  XLALDestroyCOMPLEX8TimeSeries ( sTS );	/* short invSFT timeseries */
  XLALDestroyCOMPLEX8FFTPlan ( SFTplan );
  XLALDestroyCOMPLEX8TimeSeries ( lTS );

  XLAL_ERROR_NULL ( fn, XLAL_EFUNC );

 success:
  /* cleanup memory */
  XLALDestroyCOMPLEX8TimeSeries ( sTS );	/* short invSFT timeseries */
  XLALDestroyCOMPLEX8FFTPlan ( SFTplan );

  return lTS;

} /* XLALSFTVectorToCOMPLEX8TimeSeries() */

/** Finds the earliest timestamp in a multi-SFT data structure 
 *
*/
int XLALEarliestMultiSFTsample ( LIGOTimeGPS *out,              /**< [out] earliest GPS time */
				 MultiSFTVector *multisfts      /**< [in] multi SFT vector */
				 )
{
  static const CHAR *fn = "XLALEarliestMultiSFTsample()";

  UINT4 i,j;

  /* check sanity of input */
  if ( !multisfts || (multisfts->length == 0) )
    {
      XLALPrintError ("%s: empty multiSFT input!\n", fn );
      XLAL_ERROR (fn, XLAL_EINVAL);
    }
  for (i=0;i<multisfts->length;i++)
    {
      if ( !multisfts->data[i] || (multisfts->data[i]->length == 0) )
	{
	  XLALPrintError ("%s: empty multiSFT->data[%d] input!\n", fn,i );
	  XLAL_ERROR (fn, XLAL_EINVAL);
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

/** Find the time of the end of the latest SFT in a multi-SFT data structure 
 *
*/
int XLALLatestMultiSFTsample ( LIGOTimeGPS *out,              /**< [out] latest GPS time */
			       MultiSFTVector *multisfts      /**< [in] multi SFT vector */
			       )
{
  static const CHAR *fn = "XLALLatestMultiSFTsample()";

  UINT4 i,j;
  SFTtype *firstSFT;
  REAL8 Tsft;

  /* check sanity of input */
  if ( !multisfts || (multisfts->length == 0) )
    {
      XLALPrintError ("%s: empty multiSFT input!\n", fn );
      XLAL_ERROR (fn, XLAL_EINVAL);
    }
  for (i=0;i<multisfts->length;i++)
    {
      if ( !multisfts->data[i] || (multisfts->data[i]->length == 0) )
	{
	  XLALPrintError ("%s: empty multiSFT->data[%d] input!\n", fn,i );
	  XLAL_ERROR (fn, XLAL_EINVAL);
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
      XLALPrintError ("%s: NULL pointer returned from XLALGPSAdd()!\n", fn );
      XLAL_ERROR (fn, XLAL_EFAULT);
    }

  /* success */
  return XLAL_SUCCESS;

} /* XLALLatestMultiSFTsample() */

/** Computed the weighted timeseries Fa(t) = x(t).a(t) and Fb(t) = x(t).b(t) for a multi-detector timeseries
 *
*/
int XLALAntennaWeightCOMPLEX8TimeSeries ( COMPLEX8TimeSeries **Faoft,                         /**< [out] the timeseries weighted by a(t) */
					  COMPLEX8TimeSeries **Fboft,                         /**< [out] the timeseries weighted by b(t) */
					  const COMPLEX8TimeSeries *timeseries,         /**< [in] the input timeseries */
					  const AMCoeffs *AMcoef,                       /**< [in] the AM coefficients */
					  const SFTVector *sfts                         /**< [in] the SFT data */
					  )
{
  static const CHAR *fn = "XLALAntennaWeightCOMPLEX8TimeSeries()";
  UINT4 j,k;

  /* do sanity checks */
  REAL8 start = GPS2REAL8(timeseries->epoch);
  REAL8 fHet = timeseries->f0;
  REAL8 deltaT = timeseries->deltaT;
  UINT4 numTimeSamples = timeseries->data->length;
  REAL8 dfSFT = sfts->data[0].deltaF;
  REAL8 Tsft = 1.0 / dfSFT;
  UINT4 nbins = (UINT4)floor(0.5 + Tsft/deltaT);
 
  /* create empty timeseries structures for Fa(t) and Fb(t) */
  if ( ((*Faoft) = XLALCreateCOMPLEX8TimeSeries ( sfts->data[0].name, &(timeseries->epoch), fHet, deltaT, &empty_LALUnit, numTimeSamples )) == NULL ) {
    XLALPrintError ("%s: XLALCreateCOMPLEX8TimeSeries() for %d timesteps failed! errno = %d!\n", fn, numTimeSamples, xlalErrno );
    XLAL_ERROR (fn, XLAL_ENOMEM);
  }
  if ( ((*Fboft) = XLALCreateCOMPLEX8TimeSeries ( sfts->data[0].name,&(timeseries->epoch) , fHet, deltaT, &empty_LALUnit, numTimeSamples )) == NULL ) {
    XLALPrintError ("%s: XLALCreateCOMPLEX8TimeSeries() for %d timesteps failed! errno = %d!\n", fn, numTimeSamples, xlalErrno );
    XLAL_ERROR (fn, XLAL_ENOMEM);
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
      (*Faoft)->data->data[time_index].re = a*timeseries->data->data[time_index].re;
      (*Faoft)->data->data[time_index].im = a*timeseries->data->data[time_index].im;
      (*Fboft)->data->data[time_index].re = b*timeseries->data->data[time_index].re;
      (*Fboft)->data->data[time_index].im = b*timeseries->data->data[time_index].im;

      }
      
    }
 
  /* success */
  return XLAL_SUCCESS;

}

/** Computed the weighted timeseries Fa(t) = x(t).a(t) and Fb(t) = x(t).b(t) for a multi-detector timeseries
 *
*/
int XLALAntennaWeightMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries **Faoft,                         /**< [out] the timeseries weighted by a(t) */
					       MultiCOMPLEX8TimeSeries **Fboft,                         /**< [out] the timeseries weighted by b(t) */
					       const MultiCOMPLEX8TimeSeries *multiTimeseries,         /**< [in] the input multi-detector timeseries */
					       const MultiAMCoeffs *multiAMcoef,                       /**< [in] the multi-detector AM coefficients */
					       const MultiSFTVector *multisfts                         /**< [in] the multi-detector SFT data */
					       )
{
  static const CHAR *fn = "XLALAntennaWeightMultiCOMPLEX8TimeSeries()";
 
  UINT4 i;
  
  /* do sanity checks */
 
  
  /* allocate memory for the output structure */
  if (((*Faoft) = XLALMalloc(sizeof(MultiCOMPLEX8TimeSeries))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", fn, sizeof(MultiCOMPLEX8TimeSeries));
    XLAL_ERROR (fn, XLAL_ENOMEM);
  }
  if (((*Fboft) = XLALMalloc(sizeof(MultiCOMPLEX8TimeSeries))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", fn, sizeof(MultiCOMPLEX8TimeSeries));
    XLAL_ERROR (fn, XLAL_ENOMEM);
  }
  (*Faoft)->length = multisfts->length;
  (*Fboft)->length = multisfts->length;
  
  if (((*Faoft)->data = XLALMalloc(multisfts->length*sizeof(COMPLEX8TimeSeries*))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", fn, multisfts->length*sizeof(COMPLEX8TimeSeries*));
    XLAL_ERROR (fn, XLAL_ENOMEM);
  }
  if (((*Fboft)->data = XLALMalloc(multisfts->length*sizeof(COMPLEX8TimeSeries*))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", fn, multisfts->length*sizeof(COMPLEX8TimeSeries*));
    XLAL_ERROR (fn, XLAL_ENOMEM);
  }

  /* loop over detectors */
  for (i=0;i<multisfts->length;i++) {

    /* point to current detector params */
    COMPLEX8TimeSeries *timeseries = multiTimeseries->data[i];
    AMCoeffs *AMcoef = multiAMcoef->data[i];
    SFTVector *SFTs = multisfts->data[i];

    if ( XLALAntennaWeightCOMPLEX8TimeSeries(&((*Faoft)->data[i]),&((*Fboft)->data[i]),timeseries,AMcoef,SFTs) != XLAL_SUCCESS ) {
      LALPrintError("\nXLALAntennaWeightMultiCOMPLEX8TimeSeries() failed with error = %d\n\n", xlalErrno );
      XLAL_ERROR (fn, XLAL_EFAULT);
    }
   
   
  }

  /* success */
  return XLAL_SUCCESS;
  
}

/** Performs barycentric resampling on a multi-detector timeseries 
 *
*/
int XLALBarycentricResampleMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries **Faoft_RS,                         /**< [out] the timeseries weighted by a(t) */
						     MultiCOMPLEX8TimeSeries **Fboft_RS,                         /**< [out] the timeseries weighted by b(t) */
						     const MultiCOMPLEX8TimeSeries *Faoft,                   /**< [in] the input multi-detector timeseries */
						     const MultiCOMPLEX8TimeSeries *Fboft,                       /**< [in] the multi-detector AM coefficients */
						     const MultiSSBtimes *multiSSB,                         /**< [in] the multi-detector SFT data */
						     const MultiSFTVector *multiSFTs,                         /**< [in] the multi-detector SFT data */
						     const REAL8 Tsft,                                             /**< [in] the length of an SFT */
						     const REAL8 deltaF                                             /**< [in] the user defined frequency resolution */ 
						     )
{
  static const CHAR *fn = "XLALBarycentricResampleWeightMultiCOMPLEX8TimeSeries()";

  UINT4 i;
  LIGOTimeGPS earliest,latest;
  UINT4 numTimeSamples;
  REAL8 Tspan;
  REAL8 deltaT;
  REAL8 Teff;
  REAL8 fHet;
  
  /* do sanity checks */
  if ( !Faoft || (Faoft->length == 0) ) {
    XLALPrintError ("%s: empty Faoft input!\n", fn );
    XLAL_ERROR (fn, XLAL_EINVAL);
  }
  /* do sanity checks */
  if ( !Fboft || (Fboft->length == 0) ) {
    XLALPrintError ("%s: empty Fboft input!\n", fn );
    XLAL_ERROR (fn, XLAL_EINVAL);
  }

  /* find earliest SSB time */
  if ( (XLALEarliestMultiSSBtime (&earliest,multiSSB,Tsft)) != XLAL_SUCCESS ) {
    LALPrintError("\nXLALEarliestMultiSSBtime() failed with error = %d\n\n", xlalErrno );
    XLAL_ERROR (fn, XLAL_EFAULT);
  }

  /* find latest SSB time */
  if ( (XLALLatestMultiSSBtime (&latest,multiSSB,Tsft)) != XLAL_SUCCESS ) {
    LALPrintError("\nXLALLatestMultiSSBtime() failed with error = %d\n\n", xlalErrno );
    XLAL_ERROR (fn, XLAL_EFAULT);
  }

  /* determine resampled timeseries parameters */
  Tspan = XLALGPSDiff(&latest,&earliest);     /* the time span from the earliest sample to the end of the latest sample over *all* detectors */
  deltaT = Faoft->data[0]->deltaT;      /* the sample rate of the downsampled detector frame timeseries */
  Teff = 1.0/deltaF;                          /* the effective observation time based on the requested frequency resolution (for zero padding) */
  fHet = Faoft->data[0]->f0;                  /* the input timeseries heterodyne frequency */

  /* if the requested frequency resolution gives an effective observation time less than the actual data span then we exit */
  if (Tspan > Teff) {
    XLALPrintError ("%s: requested frequency resolution too coarse for resampling!\n", fn );
    XLAL_ERROR (fn, XLAL_EINVAL);
  }

  /* redefine sample rate and compute number of samples in the new timeseries */
  numTimeSamples = (UINT4)ceil(Teff/deltaT);      /* we use ceil() so that we artificially widen the band rather than reduce it */
  deltaT = Teff/numTimeSamples;

  /* allocate memory for the output structure */
  if (((*Faoft_RS) = XLALMalloc(sizeof(MultiCOMPLEX8TimeSeries))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", fn, sizeof(MultiCOMPLEX8TimeSeries));
    XLAL_ERROR (fn, XLAL_ENOMEM);
  }
  if (((*Fboft_RS) = XLALMalloc(sizeof(MultiCOMPLEX8TimeSeries))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", fn, sizeof(MultiCOMPLEX8TimeSeries));
    XLAL_ERROR (fn, XLAL_ENOMEM);
  }
  (*Faoft_RS)->length = multiSSB->length;
  (*Fboft_RS)->length = multiSSB->length;
 
  if (((*Faoft_RS)->data = XLALMalloc(multiSSB->length*sizeof(COMPLEX8TimeSeries*))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", fn, multiSSB->length*sizeof(COMPLEX8TimeSeries*));
    XLAL_ERROR (fn, XLAL_ENOMEM);
  }
  if (((*Fboft_RS)->data = XLALMalloc(multiSSB->length*sizeof(COMPLEX8TimeSeries*))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", fn, multiSSB->length*sizeof(COMPLEX8TimeSeries*));
    XLAL_ERROR (fn, XLAL_ENOMEM);
  }
  
  /* loop over detectors */
  for (i=0;i<multiSSB->length;i++) {

    SSBtimes *SSB = multiSSB->data[i];
    SFTVector *SFTs = multiSFTs->data[i];
    COMPLEX8TimeSeries *Fa = Faoft->data[i];
    COMPLEX8TimeSeries *Fb = Fboft->data[i];

    /* create empty timeseries structures for the resampled Fa(t) and Fb(t) */
    if ( ((*Faoft_RS)->data[i] = XLALCreateCOMPLEX8TimeSeries ( Faoft->data[i]->name, &earliest, fHet, deltaT, &empty_LALUnit, numTimeSamples )) == NULL ) {
      XLALPrintError ("%s: XLALCreateCOMPLEX8TimeSeries() for %d timesteps failed! errno = %d!\n", fn, numTimeSamples, xlalErrno );
      XLAL_ERROR (fn, XLAL_ENOMEM);
    }
    if ( ((*Fboft_RS)->data[i] = XLALCreateCOMPLEX8TimeSeries ( Fboft->data[i]->name, &earliest , fHet, deltaT, &empty_LALUnit, numTimeSamples )) == NULL ) {
      XLALPrintError ("%s: XLALCreateCOMPLEX8TimeSeries() for %d timesteps failed! errno = %d!\n", fn, numTimeSamples, xlalErrno );
      XLAL_ERROR (fn, XLAL_ENOMEM);
    }
    memset ( (*Faoft_RS)->data[i]->data->data, 0, numTimeSamples * sizeof(*(*Faoft_RS)->data[i]->data->data)); 	/* set all time-samples to zero (in case there are gaps) */
    memset ( (*Fboft_RS)->data[i]->data->data, 0, numTimeSamples * sizeof(*(*Fboft_RS)->data[i]->data->data)); 	/* set all time-samples to zero (in case there are gaps) */

    /* perform resampling on current detector timeseries */ 
    if ( XLALBarycentricResampleCOMPLEX8TimeSeries(&((*Faoft_RS)->data[i]),&((*Fboft_RS)->data[i]),Fa,Fb,SSB,SFTs,Tsft) != XLAL_SUCCESS ) {
      LALPrintError("\nXLALAntennaWeightMultiCOMPLEX8TimeSeries() failed with error = %d\n\n", xlalErrno );
      XLAL_ERROR (fn, XLAL_EFAULT);
    }
   
  }
 
  /* success */
  return XLAL_SUCCESS;
  
}

/**  Performs barycentric resampling on a timeseries 
 *
 * We expect that the output timeseries has already been allocated correctly
 *
*/
int XLALBarycentricResampleCOMPLEX8TimeSeries ( COMPLEX8TimeSeries **Faoft_RS,                         /**< [out] the timeseries weighted by a(t) */
						COMPLEX8TimeSeries **Fboft_RS,                         /**< [out] the timeseries weighted by b(t) */
						const COMPLEX8TimeSeries *Faoft,                      /**< [in] the input timeseries */
						const COMPLEX8TimeSeries *Fboft,                      /**< [in] the AM coefficients */
						const SSBtimes *SSB,                                   /**< [in] the SFT data */
						const SFTVector *SFTs,
						const REAL8 Tsft
						)
{
  static const CHAR *fn = "XLALBarycentricResampleCOMPLEX8TimeSeries()";

  UINT4 i,j,k;
  REAL8Vector *detectortimes = NULL;
  REAL8Vector* FaFb[4];          /* used to store Fa and Fb real and imaginary parts in 4 seperate vectors */
  REAL8Vector *t_DET = NULL;
  gsl_spline* spline_FaFb[4];           /* used to store spline coefficients for Fa and Fb real and imaginary parts in 4 seperate vectors */

  /* do sanity checks */
  printf("sanity\n");
  
  /* define some useful shorthands */
  REAL8 refTime = GPS2REAL8(SSB->refTime);
  REAL8 start_DET = GPS2REAL8(Faoft->epoch);
  REAL8 fHet = Faoft->f0;
  REAL8 deltaT_DET = Faoft->deltaT;
  UINT4 numTimeSamples_DET = Faoft->data->length;
  REAL8 start_SSB = GPS2REAL8((*Faoft_RS)->epoch);
  REAL8 deltaT_SSB = (*Faoft_RS)->deltaT;
  
  /* allocate memory for the uniformly sampled detector time samples for this SFT */
  for (i=0;i<4;i++) {
    if ( (FaFb[i] = XLALCreateREAL8Vector(numTimeSamples_DET)) == NULL ) {
      XLALPrintError ("%s: XLALCreateREAL8Vector() failed!  errno = %d!\n", fn, xlalErrno );
      XLAL_ERROR (fn, XLAL_ENOMEM);
    }   
  }
  if ( (t_DET = XLALCreateREAL8Vector(numTimeSamples_DET)) == NULL ) { 
    XLALPrintError ("%s: XLALCreateREAL8Vector() failed!  errno = %d!\n", fn, xlalErrno ); 
    XLAL_ERROR (fn, XLAL_ENOMEM); 
  }    
  
  /* place the timeseries into REAL8Vectors for gsl to be able to interpolate them */
  /* this is annoying because the data is currently stored in a COMPLEX8 format */
  for (j=0;j<numTimeSamples_DET;j++) {
    t_DET->data[j] = start_DET + j*deltaT_DET;
    FaFb[0]->data[j] = Faoft->data->data[j].re;
    FaFb[1]->data[j] = Faoft->data->data[j].im;
    FaFb[2]->data[j] = Fboft->data->data[j].re;
    FaFb[3]->data[j] = Fboft->data->data[j].im;
  }
  
  /* initialise the gsl spline interpolation for each of the 4 timeseries */
  for (i=0;i<4;i++) {
    if ( (XLALGSLInitInterpolateREAL8Vector(&(spline_FaFb[i]),t_DET,FaFb[i])) != XLAL_SUCCESS ) {
      LALPrintError("\nXLALGSLInitInterpolateREAL8Vector() failed with error = %d\n\n", xlalErrno );
      XLAL_ERROR (fn, XLAL_EFAULT);
    }
  }
	
 /*  { */
/*     /\* TESTING *\/ */
/*     FILE *fp = NULL; */
/*     fp = fopen("/home/chmess/test/out/pre_interp.txt","w"); */
/*     for (i=0;i<numTimeSamples_DET;i++) fprintf(fp,"%6.12f %6.12f %6.12f %6.12f %6.12f\n", */
/* 					       t_DET->data[i], */
/* 					       FaFb[0]->data[i],FaFb[1]->data[i], */
/* 					       FaFb[2]->data[i],FaFb[3]->data[i]); */
/*     fclose(fp); */
/*   } */
			  
  /* loop over SFTs to compute the detector frame time samples corresponding to the uniformly sampled SSB time samples */
  for (j=0;j<SSB->DeltaT->length;j++) { 
    
    REAL8Vector *out_FaFb[4];                                                              /* output vectors for the real and imaginary parts of Fa and Fb */
    REAL8 Tdot = SSB->Tdot->data[j];                                                       /* the instantaneous time derivitive dt_SSB/dt_DET at the MID-POINT of the SFT */
    REAL8 SFTmid_SSB = refTime + SSB->DeltaT->data[j];                                     /* MID-POINT time of the SFT at the SSB */
    REAL8 SFTstart_SSB = SFTmid_SSB - 0.5*Tsft*Tdot;                                       /* START time of the SFT at the SSB */
    REAL8 SFTend_SSB = SFTmid_SSB + 0.5*Tsft*Tdot;                                         /* END time of the SFT at the SSB */
    REAL8 SFTstart_DET = GPS2REAL8(SFTs->data[j].epoch);                                   /* START time of the SFT at the detector */
    REAL8 SFTmid_DET = SFTstart_DET + 0.5*Tsft;                                            /* MID-POINT time of the SFT at the detector */
  /*   REAL8 Tsft_SSB = Tsft*Tdot;    */                                                         /* the approximated span of the SFT in the SSB -> tSFT_SSB = tSFT_DET*(dt_SSB/dt_DET) */

    UINT4 idx_start_SSB = floor(0.5 + (SFTstart_SSB - start_SSB)/deltaT_SSB);                  /* the index of the resampled timeseries corresponding to the start of the SFT */
    UINT4 idx_end_SSB = floor(0.5 + (SFTend_SSB - start_SSB)/deltaT_SSB);                      /* the index of the resampled timeseries corresponding to the end of the SFT */
    UINT4 numSamples_SSB = idx_end_SSB - idx_start_SSB + 1;                                    /* the number of samples in the SSB for this SFT */
    
  /*   printf("SFT epoch = %d %d\n", SFTs->data[j].epoch.gpsSeconds,SFTs->data[j].epoch.gpsNanoSeconds); */
/*     printf("SFTmid_DET = %6.12f SFTmid_SSB = %6.12f\n",SFTmid_DET,SFTmid_SSB); */
/*     printf("SFTstart_DET = %6.12f SFTstart_SSB = %6.12f\n",SFTstart_DET,SFTstart_SSB); */
/*     printf("idx_start_SSB = %d idx_end_SSB = %d\n",idx_start_SSB,idx_end_SSB); */
/*     printf("Tsft_SSB = %6.12f\n",Tsft_SSB); */
/*     printf("Tdot = %6.12f\n",Tdot); */
    
    /* allocate memory for the detector time samples for this SFT */
    if ( (detectortimes = XLALCreateREAL8Vector(numSamples_SSB)) == NULL ) {
      XLALPrintError ("%s: XLALCreateREAL8Vector() failed!  errno = %d!\n", fn, xlalErrno );
      XLAL_ERROR (fn, XLAL_ENOMEM);
    }   

    /* for each time sample in the SSB frame for this SFT we estimate the detector time */
    for (k=0;k<numSamples_SSB;k++) {
      
      REAL8 t_SSB = start_SSB + (k+idx_start_SSB)*deltaT_SSB;                                  /* the SSB time of the current resampled time sample */
      detectortimes->data[k] = SFTmid_DET + (t_SSB - SFTmid_SSB)/Tdot;                         /* the approximated DET time of the current resampled time sample */
 
    }
 
    /* interpolate on the non-uniformly sampled detector time vector for this SFT for re and im parts of Fa and Fb */
    for (i=0;i<4;i++) {
      if ( XLALGSLInterpolateREAL8Vector(&(out_FaFb[i]),detectortimes,spline_FaFb[i]) != XLAL_SUCCESS ) {
	LALPrintError("\nXLALInterpolateMultiCOMPLEX8TimeSeries() failed with error = %d\n\n", xlalErrno );
	XLAL_ERROR (fn, XLAL_EFAULT);
      }
    }
       
    /* { */
/*       /\* TESTING *\/ */
/*       FILE *fp = NULL; */
/*       fp = fopen("/home/chmess/test/out/post_interp.txt","a"); */
/*       for (k=0;k<numSamples_SSB;k++) fprintf(fp,"%6.12f %6.12f %6.12f %6.12f %6.12f\n", */
/* 						 detectortimes->data[k], */
/* 						 out_FaFb[0]->data[k],out_FaFb[1]->data[k], */
/* 						 out_FaFb[2]->data[k],out_FaFb[3]->data[k]); */
/*       fclose(fp); */
/*     } */
  

    /* place these interpolated timeseries into the output */
    /* and apply correction due to non-zero heterodyne frequency of input */
    for (k=0;k<numSamples_SSB;k++) {

      UINT4 idx = k + idx_start_SSB;                                                                     /* the full resampled timeseries index */
      REAL8 tDiff = start_SSB + idx*deltaT_SSB - detectortimes->data[k];                                 /* the difference between t_SSB and t_DET */
      REAL8 phase = -LAL_TWOPI*fHet*tDiff;                                                               /* the accumulated heterodyne phase */            
      REAL8 cosphase = cos(phase);                                                                       /* the real part of the phase correction */
      REAL8 sinphase = sin(phase);                                                                       /* the imaginary part of the phase correction */
      /* printf("j = %d t = %6.12f tb = %6.12f tDiff = %6.12f\n",j,detectortimes->data[k],start_SSB + idx*deltaT_SSB,tDiff); */
      (*Faoft_RS)->data->data[idx].re = out_FaFb[0]->data[k]*cosphase - out_FaFb[1]->data[k]*sinphase;
      (*Faoft_RS)->data->data[idx].im = out_FaFb[1]->data[k]*cosphase + out_FaFb[0]->data[k]*sinphase;
      (*Fboft_RS)->data->data[idx].re = out_FaFb[2]->data[k]*cosphase - out_FaFb[3]->data[k]*sinphase;
      (*Fboft_RS)->data->data[idx].im = out_FaFb[3]->data[k]*cosphase + out_FaFb[2]->data[k]*sinphase;

    }

    /* free memory used for this SFT */
    for (i=0;i<4;i++) XLALDestroyREAL8Vector(out_FaFb[i]);  
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

/** Find the earliest timestamp in a multi-SSB data structure 
 *
*/
int XLALEarliestMultiSSBtime ( LIGOTimeGPS *out,              /**< output earliest GPS time */
			       const MultiSSBtimes *multiSSB,      /**< input multi SSB SFT-midpoint timestamps */
			       const REAL8 Tsft                    /**< the length of an SFT */ 
			       )
{
  static const CHAR *fn = "XLALEarliestMultiSSBtime()";

  UINT4 i,j;
  LIGOTimeGPS t;
  REAL8 delta;

  /* check sanity of input */
  if ( !multiSSB || (multiSSB->length == 0) )
    {
      XLALPrintError ("%s: empty multiSSBtimes input!\n", fn );
      XLAL_ERROR (fn, XLAL_EINVAL);
    }
  if ( !multiSSB->data[0] || (multiSSB->data[0]->DeltaT->length == 0) )
    {
      XLALPrintError ("%s: empty multiSSBtimes input!\n", fn );
      XLAL_ERROR (fn, XLAL_EINVAL);
    }
  
 
  /* initialise the earliest and latest sample value */
  out->gpsSeconds = multiSSB->data[0]->refTime.gpsSeconds;
  out->gpsNanoSeconds = multiSSB->data[0]->refTime.gpsNanoSeconds;
  delta = multiSSB->data[0]->DeltaT->data[0] - 0.5*Tsft*multiSSB->data[0]->Tdot->data[0];
  if ( (XLALGPSAdd(out,delta)) == NULL) {
    XLALPrintError ("%s: XLALGPSAdd() failed!  errno = %d!\n", fn, xlalErrno );
    XLAL_ERROR (fn, XLAL_ENOMEM);
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
	XLALPrintError ("%s: XLALGPSAdd() failed!  errno = %d!\n", fn, xlalErrno );
	XLAL_ERROR (fn, XLAL_ENOMEM);
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

/** Find the latest timestamp in a multi-SSB data structure 
 *
*/ 
int XLALLatestMultiSSBtime ( LIGOTimeGPS *out,                   /**< output latest GPS time */
			     const MultiSSBtimes *multiSSB,      /**< input multi SSB SFT-midpoint timestamps */
			     const REAL8 Tsft                    /**< the length of an SFT */ 
			     )
{
  static const CHAR *fn = "XLALLatestMultiSSBtime()";

  UINT4 i,j;
  LIGOTimeGPS t;
  REAL8 delta;

  /* check sanity of input */
  if ( !multiSSB || (multiSSB->length == 0) )
    {
      XLALPrintError ("%s: empty multiSSBtimes input!\n", fn );
      XLAL_ERROR (fn, XLAL_EINVAL);
    }
  if ( !multiSSB->data[0] || (multiSSB->data[0]->DeltaT->length == 0) )
    {
      XLALPrintError ("%s: empty multiSSBtimes input!\n", fn );
      XLAL_ERROR (fn, XLAL_EINVAL);
    }
  
 
  /* initialise the earliest and latest sample value */
  out->gpsSeconds = multiSSB->data[0]->refTime.gpsSeconds;
  out->gpsNanoSeconds = multiSSB->data[0]->refTime.gpsNanoSeconds;
  delta = multiSSB->data[0]->DeltaT->data[0] + 0.5*Tsft*multiSSB->data[0]->Tdot->data[0];
  if ( (XLALGPSAdd(out,delta)) == NULL) {
    XLALPrintError ("%s: XLALGPSAdd() failed!  errno = %d!\n", fn, xlalErrno );
    XLAL_ERROR (fn, XLAL_ENOMEM);
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
	XLALPrintError ("%s: XLALGPSAdd() failed!  errno = %d!\n", fn, xlalErrno );
	XLAL_ERROR (fn, XLAL_ENOMEM);
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

/** Find the latest timestamp in a multi-SSB data structure 
 *
*/ 
int XLALGSLInterpolateREAL8Vector( REAL8Vector **yi, 
				   REAL8Vector *xi, 
				   gsl_spline *spline
				   )
 {
   static const CHAR *fn = "XLALGSLInterpolateREAL8Vector()";
     
   /* check input */

   UINT4 i;
   UINT4 numSamples = xi->length; 
   gsl_interp_accel *acc = gsl_interp_accel_alloc();

   /* allocate memory for output vector */
   if ( ((*yi) = XLALCreateREAL8Vector(numSamples)) == NULL ) {
     XLALPrintError ("%s: XLALCreateREAL8Vector() failed!  errno = %d!\n", fn, xlalErrno );
     XLAL_ERROR (fn, XLAL_ENOMEM);
   }

   /* perform inerpolation */
   for (i=0;i<numSamples;i++) {
     (*yi)->data[i] = gsl_spline_eval(spline,xi->data[i],acc);
   }
   
   /* free memory */
   gsl_interp_accel_free(acc);

   return XLAL_SUCCESS;
 
}


/** Find the latest timestamp in a multi-SSB data structure 
 *
*/ 
int XLALGSLInitInterpolateREAL8Vector( gsl_spline **spline, 
				       REAL8Vector *x, 
				       REAL8Vector *y
				       )
 {
   static const CHAR *fn = "XLALGSLInitInterpolateREAL8Vector()";
     
   /* check input */

   UINT4 numSamples_in = x->length; 
   REAL8 *xtemp = x->data;
   REAL8 *ytemp = y->data;

   /* compute spline interpolation coefficients */
   if ( ((*spline) = gsl_spline_alloc(gsl_interp_cspline,numSamples_in)) == NULL ) {
     XLALPrintError ("%s: XLALInitGSLInterpolateREAL8Vector() failed!  errno = %d!\n", fn, xlalErrno );
     XLAL_ERROR (fn, XLAL_ENOMEM);
   }
   gsl_spline_init((*spline),xtemp,ytemp,numSamples_in);

   return XLAL_SUCCESS;
 
}

/** Shifts an FFT output vector such that the Niquest frequency is the central bin 
 *
*/ 
int XLALFFTShiftCOMPLEX8Vector(COMPLEX8Vector **x)
{
  static const CHAR *fn = "XLALFFTShiftCOMPLEX8Vector()";

  UINT4 N = (*x)->length;
  UINT4 NQ = NhalfPosDC(N);
  UINT4 NminusNQ = N - NQ;
 /*  printf("NQ = %d NminusNQ = %d N = %d\n",NQ,NminusNQ,N); */
       
  /* allocate temp memory */
  COMPLEX8 *temp = XLALMalloc(NQ*sizeof(COMPLEX8));
  
  /* copy the bins 0 -> NQ - 1 to temp memory */
  if ( memcpy(temp,(*x)->data,NQ*sizeof(COMPLEX8)) == NULL ) {
    XLALPrintError ("%s: memcpy() failed!  errno = %d!\n", fn, xlalErrno );
    XLAL_ERROR (fn, XLAL_ENOMEM);
  }
  
  /* copy the bins NQ -> N - 1 to the start */
  if (memcpy((*x)->data,&((*x)->data[NQ]),NminusNQ*sizeof(COMPLEX8)) == NULL ) {
    XLALPrintError ("%s: memcpy() failed!  errno = %d!\n", fn, xlalErrno );
    XLAL_ERROR (fn, XLAL_ENOMEM);
  }
  
  /* copy the temp bins to the end */
  if (memcpy(&((*x)->data[NminusNQ]),temp,NQ*sizeof(COMPLEX8)) == NULL ) {
    XLALPrintError ("%s: memcpy() failed!  errno = %d!\n", fn, xlalErrno );
    XLAL_ERROR (fn, XLAL_ENOMEM);
  }
  
  /* free temp memory */
  XLALFree(temp);
  
  return XLAL_SUCCESS;
  
}

/** Multi-detector wrapper for XLALFrequencyShiftCOMPLEX8TimeSeries
 *
 * NOTE: this <b>modifies</b> the MultiCOMPLEX8Timeseries in place
 */
int
XLALFrequencyShiftMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries **x,	/**< [in/out] timeseries to time-shift */
					    REAL8 shift )	        /**< freq-shift in Hz */
{
  const CHAR *fn = "XLALFrequencyShiftMultiCOMPLEX8TimeSeries()";
  UINT4 i;
  
  if ( !(*x) )
    {
      XLALPrintError ("%s: empty input COMPLEX8timeSeries!\n", fn );
      XLAL_ERROR (fn, XLAL_EINVAL);
    }

  /* loop over detectors */
  for ( i=0; i < (*x)->length; i++)
    {
      
      if ( XLALFrequencyShiftCOMPLEX8TimeSeries ( &((*x)->data[i]), shift) != XLAL_SUCCESS ) {
	XLALPrintError ("%s: memcpy() failed!  errno = %d!\n", fn, xlalErrno );
	XLAL_ERROR (fn, XLAL_EFAULT);
      }
      
    }
  
  return XLAL_SUCCESS;
  
} /* XLALFrequencyShiftMultiCOMPLEX8TimeSeries() */

/** Freq-shift the given COMPLEX8Timeseries by an amount of 'shift' Hz,
 * using the time-domain expression y(t) = x(t) * e^(-i 2pi df t),
 * which shifts x(f) into y(f) = x(f - df)
 *
 * NOTE: this <b>modifies</b> the COMPLEX8TimeSeries in place
 */
int
XLALFrequencyShiftCOMPLEX8TimeSeries ( COMPLEX8TimeSeries **x,	/**< [in/out] timeseries to time-shift */
				       REAL8 shift )	        /**< freq-shift in Hz */
{
  const CHAR *fn = "XLALFrequencyShiftCOMPLEX8TimeSeries()";
  UINT4 k;
  REAL8 deltat;
  
  if ( !(*x) || !(*x)->data )
    {
      XLALPrintError ("%s: empty input COMPLEX8TimeSeries!\n", fn );
      XLAL_ERROR (fn, XLAL_EINVAL);
    }
   
  /* get timeseries epoch */
  deltat = (*x)->deltaT;

  for ( k=0; k < (*x)->data->length; k++)
    {
      REAL8 tk = k * deltat;	/* time of k-th bin */
      REAL8 shiftCycles = shift * tk;
      REAL4 fact_re, fact_im;			/* complex phase-shift factor e^(-2pi f tau) */
      REAL4 yRe, yIm;
      /* printf("tk = %6.12f shiftCycles = %6.12f\n",tk,shiftCycles); */

      sin_cos_2PI_LUT ( &fact_im, &fact_re, shiftCycles );
      
      yRe = fact_re * (*x)->data->data[k].re - fact_im * (*x)->data->data[k].im;
      yIm = fact_re * (*x)->data->data[k].im + fact_im * (*x)->data->data[k].re;
      /* printf("data = %f %f -> %f %f\n",(*x)->data->data[k].re,(*x)->data->data[k].im,yRe,yIm); */

      (*x)->data->data[k].re = yRe;
      (*x)->data->data[k].im = yIm;
      
    } /* for k < numBins */

  /* adjust timeseries heterodyne frequency to the shift */
  (*x)->f0 -= shift;
/*   printf("new f0 = %6.12f\n",(*x)->f0); */

  return XLAL_SUCCESS;

} /* XLALFrequencyShiftCOMPLEX8TimeSeries() */

/** Apply a spin-down correction to the Fa and Fb complex timeseries
 * using the time-domain expression y(t) = x(t) * e^(-i 2pi sum f_k * (t-tref)^(k+1)),
 *
 * NOTE: this <b>modifies</b> the COMPLEX8TimeSeries Fa and Fb in place
 */
int
XLALSpinDownCorrectionMultiFaFb ( MultiCOMPLEX8TimeSeries **Fa,	/**< [in/out] timeseries to time-shift */
				  MultiCOMPLEX8TimeSeries **Fb,	/**< [in/out] timeseries to time-shift */
				  const PulsarDopplerParams *doppler		/**< parameter-space point to correct for */
				  )
{
  const CHAR *fn = "XLALSpinDownCorrectionMultiFaFb()";
  UINT4 nspins = PULSAR_MAX_SPINS - 1;
  LIGOTimeGPS *epoch;
  UINT4 N,numDetectors;
  REAL8 deltaref,deltaT;
  UINT4 i,j,k;

  /* sanity checks */
  if ( !(*Fa) || !(*Fa)->data || !(*Fb) || !(*Fb)->data ) 
    {
      XLALPrintError ("%s: empty input MultiCOMPLEX8TimeSeries!\n", fn );
      XLAL_ERROR (fn, XLAL_EINVAL);
    }
  
  /* check if Fa and Fb have the same timeseries parameters */
  epoch = &(*Fa)->data[0]->epoch;
  N = (*Fa)->data[0]->data->length;
  numDetectors = (*Fa)->length;
  deltaT = (*Fa)->data[0]->deltaT;
  if ( ( (*Fa)->length != numDetectors ) || ( (*Fb)->length != numDetectors ) ) 
    {
      XLALPrintError ("%s: Different numbers of detectors within the Fa and Fb MultiCOMPLEX8TimeSeries!\n", fn );
      XLAL_ERROR (fn, XLAL_EINVAL);
    }
  for (i=0;i<(*Fa)->length;i++) 
    {
      if ( ( XLALGPSCmp(epoch,&((*Fa)->data[i]->epoch)) != 0 ) || ( XLALGPSCmp(epoch,&((*Fb)->data[i]->epoch)) != 0 ) ) 
	{
	  XLALPrintError ("%s: Different start times for within Fa and Fb MultiCOMPLEX8TimeSeries!\n", fn );
	  XLAL_ERROR (fn, XLAL_EINVAL);
	}
      if ( ( (*Fa)->data[i]->data->length != N ) || ( (*Fb)->data[i]->data->length != N ) ) 
	{
	  XLALPrintError ("%s: Different MultiCOMPLEX8TimeSeries lengths within Fa and Fb!\n", fn );
	  XLAL_ERROR (fn, XLAL_EINVAL);
	}
      if ( ( (*Fa)->data[i]->deltaT != deltaT ) || ( (*Fb)->data[i]->deltaT != deltaT ) ) 
	{
	  XLALPrintError ("%s: Different MultiCOMPLEX8TimeSeries deltaT within Fa and Fb!\n", fn );
	  XLAL_ERROR (fn, XLAL_EINVAL);
	}
    }

  /* determine number of spin down's */
  while (doppler->fkdot[nspins]==0.0) nspins--;
 /*  printf("number of spin downs = %d\n",nspins); */
  
  /* compute the time difference between timeseries epoch and reference time */
  deltaref = XLALGPSGetREAL8(epoch) - XLALGPSGetREAL8(&(doppler->refTime));
 /*  printf("deltaref = %6.12f\n",deltaref); */
/*   printf("f0 = %6.12f\n",doppler->fkdot[0]); */
/*   printf("f1 = %6.12e\n",doppler->fkdot[1]); */

  /* apply spin derivitive correction to resampled timeseries */
  /* loop over spin derivitives (nspins = 1 means first derivitive, = 2 means second derivitive etc.. ) */
  for (j=1;j<=nspins;j++) {

    /* loop over time samples  and compute the spin down phase correction */
    for (k=0;k<N;k++) {
      
      REAL8 cycles = fmod ( 0.5*doppler->fkdot[j]*pow(deltaref + k*deltaT,(REAL8)(j+1)), 1);
      REAL4 cosphase, sinphase;

      sin_cos_2PI_LUT (&sinphase, &cosphase, -cycles );

      /* loop over detectors */
      for (i=0;i<numDetectors;i++) {
	
	REAL8 Fare = (*Fa)->data[i]->data->data[k].re*cosphase - (*Fa)->data[i]->data->data[k].im*sinphase;
	REAL8 Faim = (*Fa)->data[i]->data->data[k].im*cosphase + (*Fa)->data[i]->data->data[k].re*sinphase;
	REAL8 Fbre = (*Fb)->data[i]->data->data[k].re*cosphase - (*Fb)->data[i]->data->data[k].im*sinphase;
	REAL8 Fbim = (*Fb)->data[i]->data->data[k].im*cosphase + (*Fb)->data[i]->data->data[k].re*sinphase;
	
	(*Fa)->data[i]->data->data[k].re = Fare;
	(*Fa)->data[i]->data->data[k].im = Faim;
	(*Fb)->data[i]->data->data[k].re = Fbre;
	(*Fb)->data[i]->data->data[k].im = Fbim;

      }
    
    }
  
  }

  return XLAL_SUCCESS;

 } /* XLALSpinDownCorrectionMultiFaFb */

/* ===== Object creation/destruction functions ===== */

/** Destroy a MultiCOMPLEX8TimeSeries structure.
 * Note, this is "NULL-robust" in the sense that it will not crash
 * on NULL-entries anywhere in this struct, so it can be used
 * for failure-cleanup even on incomplete structs
 */
void
XLALDestroyMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries *multiTimes )
{
  UINT4 X;
  COMPLEX8TimeSeries *tmp;

  if ( ! multiTimes )
    return;

  if ( multiTimes->data )
    {
      for ( X=0; X < multiTimes->length; X ++ )
	{
	  if ( (tmp = multiTimes->data[X]) != NULL )
	    {
	      XLALDestroyCOMPLEX8TimeSeries(tmp);
	    } /* if multiTimes->data[X] */
	} /* for X < numDetectors */
      LALFree ( multiTimes->data );
    }
  LALFree ( multiTimes );
  
  return;
  
} /* XLALDestroyMultiCOMPLEX8TimeSeries() */
