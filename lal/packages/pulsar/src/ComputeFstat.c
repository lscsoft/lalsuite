/*
 * Copyright (C) 2006 John T. Whelan
 * Copyright (C) 2005, 2006 Reinhard Prix
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
 * \file 
 * \brief
 * Functions to calculate the so-called F-statistic for a given point in parameter-space, 
 * following the equations in \ref JKS98.
 *                                                                          
 */

/*---------- INCLUDES ----------*/
#define __USE_ISOC99 1
#include <math.h>

#include <lal/ExtrapolatePulsarSpins.h>

/* GSL includes */
#include <lal/LALGSL.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>


#include <lal/AVFactories.h>
#include "ComputeFstat.h"

NRCSID( COMPUTEFSTATC, "$Id$");

/*---------- local DEFINES ----------*/
#define TRUE (1==1)
#define FALSE (1==0)


#define LD_SMALL4       (1.0e-6)		/**< "small" number for REAL4*/
#define OOTWOPI         (1.0 / LAL_TWOPI)	/**< 1/2pi */

#define TWOPI_FLOAT     6.28318530717958f  	/**< single-precision 2*pi */
#define OOTWOPI_FLOAT   (1.0f / TWOPI_FLOAT)	/**< single-precision 1 / (2pi) */ 


/*----- Macros ----- */
/** convert GPS-time to REAL8 */
#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

#define MYSIGN(x) ( ((x) < 0) ? (-1.0):(+1.0) )

/** Simple Euklidean scalar product for two 3-dim vectors in cartesian coords */
#define SCALAR(u,v) ((u)[0]*(v)[0] + (u)[1]*(v)[1] + (u)[2]*(v)[2])

#define SQ(x) ( (x) * (x) )

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- Global variables ----------*/
#define NUM_FACT 6
static const REAL8 inv_fact[NUM_FACT] = { 1.0, 1.0, (1.0/2.0), (1.0/6.0), (1.0/24.0), (1.0/120.0) };

/* empty initializers  */
static const LALStatus empty_status;

const SSBtimes empty_SSBtimes;
const MultiSSBtimes empty_MultiSSBtimes;
const AntennaPatternMatrix empty_AntennaPatternMatrix;
const MultiAMCoeffs empty_MultiAMCoeffs;
const Fcomponents empty_Fcomponents;
const ComputeFParams empty_ComputeFParams;
const ComputeFBuffer empty_ComputeFBuffer;


/*---------- internal prototypes ----------*/
int finite(double x);

static int printGSLmatrix4 ( FILE *fp, const CHAR *prefix, const gsl_matrix *gij );

/*==================== FUNCTION DEFINITIONS ====================*/


/** Function to compute a vector of Fstatistic values for a number of frequency bins.
    This function is simply a wrapper for ComputeFstat() which is called repeatedly for
    every frequency value.  The output, i.e. fstatVector must be properly allocated
    before this function is called.  The values of the start frequency, the step size
    in the frequency and the number of frequency values for which the Fstatistic is 
    to be calculated are read from fstatVector.  The other parameters are not checked and 
    they must be correctly set outside this function. 
*/
void ComputeFStatFreqBand ( LALStatus *status, 
			    REAL8FrequencySeries *fstatVector, /**< [out] Vector of Fstat values */
			    const PulsarDopplerParams *doppler,/**< parameter-space point to compute F for */
			    const MultiSFTVector *multiSFTs, /**< normalized (by DOUBLE-sided Sn!) data-SFTs of all IFOs */
			    const MultiNoiseWeights *multiWeights,	/**< noise-weights of all SFTs */
			    const MultiDetectorStateSeries *multiDetStates,/**< 'trajectories' of the different IFOs */
			    const ComputeFParams *params	/**< addition computational params */
			    )
{

  UINT4 numDetectors, numBins, k;	
  REAL8 deltaF;
  Fcomponents Fstat;
  PulsarDopplerParams thisPoint;
  ComputeFBuffer cfBuffer = empty_ComputeFBuffer;

  INITSTATUS( status, "ComputeFStatFreqBand", COMPUTEFSTATC );
  ATTATCHSTATUSPTR (status);

  ASSERT ( multiSFTs, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( doppler, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( multiDetStates, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( params, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );

  numDetectors = multiSFTs->length;
  ASSERT ( multiDetStates->length == numDetectors, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );
  ASSERT ( fstatVector, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( fstatVector->data, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( fstatVector->data->data, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( fstatVector->data->length > 0, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );

  /* copy values from 'doppler' to local variable 'thisPoint' */
  thisPoint = *doppler;

  numBins = fstatVector->data->length;
  deltaF = fstatVector->deltaF;

  /* loop over frequency values and fill up values in fstatVector */
  for ( k = 0; k < numBins; k++) {
 
    TRY (ComputeFStat ( status->statusPtr, &Fstat, &thisPoint, multiSFTs, multiWeights, 
			multiDetStates, params, &cfBuffer ), status);

    fstatVector->data->data[k] = Fstat.F;
      
    thisPoint.fkdot[0] += deltaF;
  }

  XLALEmptyComputeFBuffer ( cfBuffer );

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* ComputeFStatFreqBand() */




/** Function to compute (multi-IFO) F-statistic for given parameter-space point ::psPoint,
 *  normalized SFT-data (normalized by <em>double-sided</em> PSD Sn), noise-weights
 *  and detector state-series 
 *
 * NOTE: for better efficiency some quantities that need to be recomputed only for different 
 * sky-positions are buffered in \a cfBuffer if given. 
 * - In order to 'empty' this buffer (at the end) use XLALEmptyComputeFBuffer()
 * - You CAN pass NULL for the \a cfBuffer if you don't want to use buffering (slower).
 *
 * NOTE2: there's a spaceholder for binary-pulsar parameters in \a psPoint, but this 
 * it not implemented yet.
 *
 */
void
ComputeFStat ( LALStatus *status, 
	       Fcomponents *Fstat,                 /**< [out] Fstatistic + Fa, Fb */
	       const PulsarDopplerParams *doppler, /**< parameter-space point to compute F for */
	       const MultiSFTVector *multiSFTs,    /**< normalized (by DOUBLE-sided Sn!) data-SFTs of all IFOs */
	       const MultiNoiseWeights *multiWeights,/**< noise-weights of all SFTs */
	       const MultiDetectorStateSeries *multiDetStates,/**< 'trajectories' of the different IFOs */
	       const ComputeFParams *params,       /**< addition computational params */
	       ComputeFBuffer *cfBuffer            /**< CF-internal buffering structure */
	       )
{
  Fcomponents retF = empty_Fcomponents;
  UINT4 X, numDetectors;	
  MultiSSBtimes *multiSSB = NULL;
  MultiAMCoeffs *multiAMcoef = NULL;
  REAL8 Ad, Bd, Cd, Dd_inv;

  INITSTATUS( status, "ComputeFStat", COMPUTEFSTATC );
  ATTATCHSTATUSPTR (status);

  /* check input */
  ASSERT ( Fstat, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( multiSFTs, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( doppler, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( multiDetStates, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( params, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );

  numDetectors = multiSFTs->length;
  ASSERT ( multiDetStates->length == numDetectors, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );
  if ( multiWeights ) {
    ASSERT ( multiWeights->length == numDetectors , status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );
  }

  if ( doppler->orbit ) {
    LALPrintError ("\nSorry, binary-pulsar search not yet implemented in LALComputeFStat()\n\n");
    ABORT ( status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );
  }

  /* check if that skyposition SSB+AMcoef were already buffered */
  if ( cfBuffer 
       && ( cfBuffer->multiDetStates == multiDetStates )
       && ( cfBuffer->Alpha == doppler->Alpha )
       && ( cfBuffer->Delta == doppler->Delta ) 
       && cfBuffer->multiSSB
       && cfBuffer->multiAMcoef )
    { /* yes ==> reuse */
      multiSSB = cfBuffer->multiSSB;
      multiAMcoef = cfBuffer -> multiAMcoef;
    }
  else 
    {
      SkyPosition skypos;
      skypos.system =   COORDINATESYSTEM_EQUATORIAL;
      skypos.longitude = doppler->Alpha;
      skypos.latitude  = doppler->Delta;
      /* compute new AM-coefficients and SSB-times */
      TRY ( LALGetMultiSSBtimes ( status->statusPtr, &multiSSB, multiDetStates, skypos, doppler->refTime, params->SSBprec ), status );

      LALGetMultiAMCoeffs ( status->statusPtr, &multiAMcoef, multiDetStates, skypos );
      BEGINFAIL ( status ) {
	XLALDestroyMultiSSBtimes ( multiSSB );
      } ENDFAIL (status);

      /* noise-weigh Antenna-patterns and compute A,B,C */
      if ( XLALWeighMultiAMCoeffs ( multiAMcoef, multiWeights ) != XLAL_SUCCESS ) {
	LALPrintError("\nXLALWeighMultiAMCoeffs() failed with error = %d\n\n", xlalErrno );
	ABORT ( status, COMPUTEFSTATC_EXLAL, COMPUTEFSTATC_MSGEXLAL );
      }

      /* store these in buffer if available */
      if ( cfBuffer )
	{
	  XLALEmptyComputeFBuffer ( *cfBuffer );
	  cfBuffer->multiSSB = multiSSB;
	  cfBuffer->multiAMcoef = multiAMcoef;
	  cfBuffer->Alpha = doppler->Alpha;
	  cfBuffer->Delta = doppler->Delta;
	  cfBuffer->multiDetStates = multiDetStates ;
	} /* if cfBuffer */

    } /* if no buffer, different skypos or different detStates */

  Ad = multiAMcoef->Mmunu.Ad;
  Bd = multiAMcoef->Mmunu.Bd;
  Cd = multiAMcoef->Mmunu.Cd;
  Dd_inv = 1.0 / (Ad * Bd - Cd * Cd );

  /* ----- loop over detectors and compute all detector-specific quantities ----- */
  for ( X=0; X < numDetectors; X ++)
    {
      Fcomponents FcX = empty_Fcomponents;	/* for detector-specific FaX, FbX */
 		  
      if ( XLALComputeFaFb (&FcX, multiSFTs->data[X], doppler->fkdot, multiSSB->data[X], multiAMcoef->data[X], params->Dterms) != 0)
	{
	  LALPrintError ("\nXALNewLALDemod() failed\n");
	  ABORT ( status, COMPUTEFSTATC_EXLAL, COMPUTEFSTATC_MSGEXLAL );
	}

#ifndef LAL_NDEBUG
      if ( !finite(FcX.Fa.re) || !finite(FcX.Fa.im) || !finite(FcX.Fb.re) || !finite(FcX.Fb.im) ) {
	LALPrintError("XLALComputeFaFb() return non-finite: Fa=(%f,%f), Fb=(%f,%f)\n", 
		      FcX.Fa.re, FcX.Fa.im, FcX.Fb.re, FcX.Fb.im );
	ABORT (status,  COMPUTEFSTATC_EIEEE,  COMPUTEFSTATC_MSGEIEEE);
      }
#endif
 		 
      /* Fa = sum_X Fa_X */
      retF.Fa.re += FcX.Fa.re;
      retF.Fa.im += FcX.Fa.im;

      /* Fb = sum_X Fb_X */ 		  
      retF.Fb.re += FcX.Fb.re;
      retF.Fb.im += FcX.Fb.im;
  		  
    } /* for  X < numDetectors */
 
  /* ----- compute final Fstatistic-value ----- */

  /* NOTE: the data MUST be normalized by the DOUBLE-SIDED PSD (using LALNormalizeMultiSFTVect),
   * therefore there is a factor of 2 difference with respect to the equations in JKS, which 
   * where based on the single-sided PSD.
   */ 
 		       
  retF.F = Dd_inv * (  Bd * (SQ(retF.Fa.re) + SQ(retF.Fa.im) ) 
                     + Ad * ( SQ(retF.Fb.re) + SQ(retF.Fb.im) )
                     - 2.0 * Cd *( retF.Fa.re * retF.Fb.re + retF.Fa.im * retF.Fb.im )  
		   );

  (*Fstat) = retF;

  /* free memory if no buffer was available */
  if ( !cfBuffer )
    {
      XLALDestroyMultiSSBtimes ( multiSSB );
      XLALDestroyMultiAMCoeffs ( multiAMcoef );
    } /* if !cfBuffer */

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* ComputeFStat() */


/** Revamped version of LALDemod() (based on TestLALDemod() in CFS).
 * Compute JKS's Fa and Fb, which are ingredients for calculating the F-statistic.
 */
int
XLALComputeFaFb ( Fcomponents *FaFb,
		  const SFTVector *sfts,
		  const PulsarSpins fkdot,
		  const SSBtimes *tSSB,
		  const AMCoeffs *amcoe,
		  UINT4 Dterms) 
{ 
  UINT4 alpha;                 	/* loop index over SFTs */
  UINT4 spdnOrder;		/* maximal spindown-orders */
  UINT4 numSFTs;		/* number of SFTs (M in the Notes) */
  COMPLEX16 Fa, Fb;
  REAL8 f;			/* !! MUST be REAL8, or precision breaks down !! */
  REAL8 Tsft; 			/* length of SFTs in seconds */
  INT4 freqIndex0;		/* index of first frequency-bin in SFTs */
  INT4 freqIndex1;		/* index of last frequency-bin in SFTs */

  REAL4 *a_al, *b_al;		/* pointer to alpha-arrays over a and b */
  REAL8 *DeltaT_al, *Tdot_al;	/* pointer to alpha-arrays of SSB-timings */
  SFTtype *SFT_al;		/* SFT alpha  */

  /* ----- check validity of input */
#ifndef LAL_NDEBUG
  if ( !FaFb ) {
    LALPrintError ("\nOutput-pointer is NULL !\n\n");
    XLAL_ERROR ( "XLALComputeFaFb", XLAL_EINVAL);
  }

  if ( !sfts || !sfts->data ) {
    LALPrintError ("\nInput SFTs are NULL!\n\n");
    XLAL_ERROR ( "XLALComputeFaFb", XLAL_EINVAL);
  }
  
  if ( !tSSB || !tSSB->DeltaT || !tSSB->Tdot || !amcoe || !amcoe->a || !amcoe->b )
    {
      LALPrintError ("\nIllegal NULL in input !\n\n");
      XLAL_ERROR ( "XLALComputeFaFb", XLAL_EINVAL);
    }

  if ( PULSAR_MAX_SPINS > NUM_FACT )
    {
      LALPrintError ("\nInverse factorials table only up to order s=%d, can't handle %d spin-order\n\n",
		     NUM_FACT, PULSAR_MAX_SPINS - 1 );
      XLAL_ERROR ( "XLALComputeFaFb", XLAL_EINVAL);
    }
#endif

  /* ----- prepare convenience variables */
  numSFTs = sfts->length;
  Tsft = 1.0 / sfts->data[0].deltaF;

  freqIndex0 = (UINT4) ( sfts->data[0].f0 / sfts->data[0].deltaF + 0.5); /* lowest freqency-index */
  freqIndex1 = freqIndex0 + sfts->data[0].data->length;

  /* find highest non-zero spindown-entry */
  for ( spdnOrder = PULSAR_MAX_SPINS - 1;  spdnOrder > 0 ; spdnOrder --  )
    if ( fkdot[spdnOrder] )
      break;

  f = fkdot[0];

  Fa.re = 0.0f;
  Fa.im = 0.0f;
  Fb.re = 0.0f;
  Fb.im = 0.0f;

  a_al = amcoe->a->data;	/* point to beginning of alpha-arrays */
  b_al = amcoe->b->data;
  DeltaT_al = tSSB->DeltaT->data;
  Tdot_al = tSSB->Tdot->data;
  SFT_al = sfts->data;

  /* Loop over all SFTs  */
  for ( alpha = 0; alpha < numSFTs; alpha++ )
    {
      REAL4 a_alpha, b_alpha;

      INT4 kstar;		/* central frequency-bin k* = round(xhat_alpha) */
      INT4 k0, k1;

      COMPLEX8 *Xalpha = SFT_al->data->data; /* pointer to current SFT-data */
      COMPLEX8 *Xalpha_l; 	/* pointer to frequency-bin k in current SFT */
      REAL4 s_alpha, c_alpha;	/* sin(2pi kappa_alpha) and (cos(2pi kappa_alpha)-1) */
      REAL4 realQ, imagQ;	/* Re and Im of Q = e^{-i 2 pi lambda_alpha} */
      REAL4 realXP, imagXP;	/* Re/Im of sum_k X_ak * P_ak */
      REAL4 realQXP, imagQXP;	/* Re/Im of Q_alpha R_alpha */

      REAL8 lambda_alpha, kappa_alpha;
      REAL8 remainder;

      /* ----- calculate kappa_alpha and lambda_alpha */
      {
	UINT4 s; 		/* loop-index over spindown-order */
	REAL8 phi_alpha, Dphi_alpha, DT_al;
	REAL8 Tas;	/* temporary variable to calculate (DeltaT_alpha)^s */
	
	/* init for s=0 */
	phi_alpha = 0.0;
	Dphi_alpha = 0.0;
	DT_al = (*DeltaT_al);
	Tas = 1.0;		/* DeltaT_alpha ^ 0 */

	for (s=0; s <= spdnOrder; s++)
	  {
	    REAL8 fsdot = fkdot[s];
	    Dphi_alpha += fsdot * Tas * inv_fact[s]; 	/* here: DT^s/s! */
	    Tas *= DT_al;				/* now: DT^(s+1) */
	    phi_alpha += fsdot * Tas * inv_fact[s+1];
	  } /* for s <= spdnOrder */

	/* Step 3: apply global factors to complete Dphi_alpha */
	Dphi_alpha *= Tsft * (*Tdot_al);		/* guaranteed > 0 ! */

	lambda_alpha = phi_alpha - 0.5 * Dphi_alpha;
	
	/* real- and imaginary part of e^{-i 2 pi lambda_alpha } */
	if ( sin_cos_2PI_LUT ( &imagQ, &realQ, - lambda_alpha ) ) {
	  XLAL_ERROR ( "XLALComputeFaFb", XLAL_EFUNC);
	}

	kstar = (INT4) (Dphi_alpha + 0.5);	/* k* = round(Dphi_alpha) for positive Dphi */
	remainder = (Dphi_alpha - kstar);	/* rem(Dphi_alpha) */
	kappa_alpha = remainder + Dterms;

	k0 = kstar - Dterms;	
	k1 = k0 + 2 * Dterms;
	if ( (k0 < freqIndex0) || (k1 > freqIndex1) ) 
	  {
	    LALPrintError ("Required frequency-bins [%d, %d] not covered by SFT-interval [%d, %d]\n\n",
			   k0, k1, freqIndex0, freqIndex1 );
	    XLAL_ERROR("XLALComputeFaFb", XLAL_EDOM);
	  }

      } /* compute kappa_alpha, lambda_alpha */

      /* NOTE: sin[ 2pi (Dphi_alpha - k) ] = sin [ 2pi Dphi_alpha ], therefore
       * the trig-functions need to be calculated only once!
       * We choose the value sin[ 2pi(Dphi_alpha - kstar) ] because it is the 
       * closest to zero and will pose no numerical difficulties !
       */
      sin_cos_2PI_LUT ( &s_alpha, &c_alpha, remainder );
      c_alpha -= 1.0f; 

      /* ---------- calculate the (truncated to Dterms) sum over k ---------- */

      /* ---------- ATTENTION: this the "hot-loop", which will be 
       * executed many millions of times, so anything in here 
       * has a HUGE impact on the whole performance of the code.
       * 
       * DON'T touch *anything* in here unless you really know 
       * what you're doing !!
       *------------------------------------------------------------
       */

      Xalpha_l = Xalpha + k0 - freqIndex0;  /* first frequency-bin in sum */

      realXP = 0;
      imagXP = 0;

      /* if no danger of denominator -> 0 */
      if ( ( remainder > LD_SMALL4 ) || (remainder < -LD_SMALL4) )	
	{ 
	  /* improved hotloop algorithm by Fekete Akos: 
	   * take out repeated divisions into a single common denominator,
	   * plus use extra cleverness to compute the nominator efficiently...
	   */
	  COMPLEX8 Xal = *Xalpha_l;
	  REAL4 Sn = Xal.re;
	  REAL4 Tn = Xal.im;
	  REAL4 pn = kappa_alpha;
	  REAL4 qn = pn;
	  REAL4 U_alpha, V_alpha;
	  
	  /* recursion with 2*Dterms steps */
	  UINT4 l;
	  for ( l = 1; l <= 2*Dterms; l ++ )
	    {
	      Xalpha_l ++;
	      Xal = *Xalpha_l;
	      
	      pn = pn - 1.0f;			/* p_(n+1) */
	      Sn = pn * Sn + qn * Xal.re;	/* S_(n+1) */
	      Tn = pn * Tn + qn * Xal.im;	/* T_(n+1) */
	      qn *= pn;				/* q_(n+1) */
	    } /* for l <= 2*Dterms */
	  
	  U_alpha = Sn / qn;
	  V_alpha = Tn / qn;

#ifndef LAL_NDEBUG
	  if ( !finite(U_alpha) || !finite(V_alpha) || !finite(pn) || !finite(qn) || !finite(Sn) || !finite(Tn) ) {
	    XLAL_ERROR ("XLALComputeFaFb()", COMPUTEFSTATC_EIEEE);
	  }
#endif

	  realXP = s_alpha * U_alpha - c_alpha * V_alpha;
	  imagXP = c_alpha * U_alpha + s_alpha * V_alpha;
	  
	} /* if |remainder| > LD_SMALL4 */
      else
	{ /* otherwise: lim_{rem->0}P_alpha,k  = 2pi delta_{k,kstar} */
	  realXP = TWOPI_FLOAT * Xalpha_l[Dterms].re;
	  imagXP = TWOPI_FLOAT * Xalpha_l[Dterms].im;
	} /* if |remainder| <= LD_SMALL4 */
      
      realQXP = realQ * realXP - imagQ * imagXP;
      imagQXP = realQ * imagXP + imagQ * realXP;
      
      /* we're done: ==> combine these into Fa and Fb */
      a_alpha = (*a_al);
      b_alpha = (*b_al);

      Fa.re += a_alpha * realQXP;
      Fa.im += a_alpha * imagQXP;
      
      Fb.re += b_alpha * realQXP;
      Fb.im += b_alpha * imagQXP;


      /* advance pointers over alpha */
      a_al ++;
      b_al ++;
      DeltaT_al ++;
      Tdot_al ++;
      SFT_al ++;

    } /* for alpha < numSFTs */
      
  /* return result */
  FaFb->Fa.re = Fa.re * OOTWOPI;
  FaFb->Fa.im = Fa.im * OOTWOPI;
  FaFb->Fb.re = Fb.re * OOTWOPI;
  FaFb->Fb.im = Fb.im * OOTWOPI;

  return XLAL_SUCCESS;

} /* XLALComputeFaFb() */

/** Compute the 'amplitude coefficients' \f$a(t), b(t)\f$ as defined in 
 * \ref JKS98 for a series of timestamps.
 * 
 * The input consists of the DetectorState-timeseries, which contains
 * the detector-info and the LMST's corresponding to the different times.
 * 
 * In order to allow re-using the output-structure AMCoeffs for subsequent
 * calls, we require the REAL4Vectors a and b to be allocated already and 
 * to have the same length as the DetectoStates-timeseries.
 *
 * \note This is an alternative implementation to LALComputeAM() with 
 * the aim to be both simpler and faster.
 * The difference being that we don't implicitly re-derive the final expression
 * here but simply try to implement the final expressions (12), (13) in \ref JKS98
 * in the most economical way possible.
 */ 
void
LALGetAMCoeffs(LALStatus *status,
	       AMCoeffs *coeffs,			/**< [out] amplitude-coeffs {a(t_i), b(t_i)} */
	       const DetectorStateSeries *DetectorStates,	/**< timeseries of detector states */
	       SkyPosition skypos			/**< {alpha,delta} of the source */
	       )
{
  REAL4 ah1, ah2, ah3, ah4, ah5;
  REAL4 a1, a2, a3, a4, a5;
  
  REAL4 bh1, bh2, bh3, bh4;
  REAL4 b1, b2, b3, b4;

  REAL4 delta, alpha;
  REAL4 sin1delta, cos1delta, sin2delta, cos2delta;

  REAL4 gamma, lambda;
  REAL4 norm;
  UINT4 i, numSteps;

  INITSTATUS (status, "LALGetAMCoeffs", COMPUTEFSTATC);

  /*---------- check input ---------- */
  ASSERT ( DetectorStates, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);

  numSteps = DetectorStates->length;

  /* require the coeffients-vectors to be allocated and consistent with timestamps */
  ASSERT ( coeffs, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( coeffs->a && coeffs->b, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( (coeffs->a->length == numSteps) && (coeffs->b->length == numSteps), status,
	   COMPUTEFSTATC_EINPUT,  COMPUTEFSTATC_MSGEINPUT);

  /* require sky-pos to be in equatorial coordinates */
  ASSERT ( skypos.system == COORDINATESYSTEM_EQUATORIAL, status, 
	   SKYCOORDINATESH_ESYS, SKYCOORDINATESH_MSGESYS );

  /*---------- detector paramters: lambda, L, gamma */
  {
    /* FIXME: put into DetectorStateSeries */
    /* orientation of detector arms */
    REAL8 xAzi = DetectorStates->detector.frDetector.xArmAzimuthRadians;
    REAL8 yAzi = DetectorStates->detector.frDetector.yArmAzimuthRadians;

    /* get detector orientation gamma */
    gamma = LAL_PI_2 - 0.5 * (xAzi + yAzi);
    /* get detector position latitude (lambda) */
    lambda = DetectorStates->detector.frDetector.vertexLatitudeRadians;
    /*
    printf ("IFO = %s: sin(zeta) = %f\n", DetectorStates->detector.frDetector.name, sin( xAzi - yAzi ) );
    */
  }

  /*---------- coefficient ahN, bhN dependent ONLY on detector-position  ---------- */
  /* FIXME: put these coefficients into DetectorStateSeries */
  {
    REAL4 sin2gamma, cos2gamma;
    REAL4 sin1lambda, cos1lambda;
    REAL4 sin2lambda, cos2lambda;

    sin_cos_LUT (&sin2gamma, &cos2gamma, 2.0f * gamma );
    sin_cos_LUT (&sin1lambda, &cos1lambda, lambda );

    sin2lambda = 2.0f * sin1lambda * cos1lambda;
    cos2lambda = cos1lambda * cos1lambda - sin1lambda * sin1lambda;

    /* coefficients for a(t) */
    ah1 = 0.0625f * sin2gamma * (3.0f - cos2lambda);	/* 1/16 = 0.0625 */
    ah2 = - 0.25f * cos2gamma * sin1lambda;
    ah3 =   0.25f * sin2gamma * sin2lambda;
    ah4 =  -0.5f  * cos2gamma * cos1lambda;
    ah5 =  0.75f  * sin2gamma * cos1lambda * cos1lambda;

    /* coefficients for b(t) */
    bh1 =           cos2gamma * sin1lambda;
    bh2 =   0.25f * sin2gamma * (3.0f - cos2lambda);
    bh3 =           cos2gamma * cos1lambda;
    bh4 =   0.5f  * sin2gamma * sin2lambda;
  }

  /*---------- coefficients aN, bN dependent ONLY on {ahN, bhN} and source-latitude delta */
  alpha = skypos.longitude;
  delta = skypos.latitude;

  sin_cos_LUT (&sin1delta, &cos1delta, delta );
  sin2delta = 2.0f * sin1delta * cos1delta;
  cos2delta = cos1delta * cos1delta - sin1delta * sin1delta;

  /* coefficients for a(t) */
  a1 = ah1 * ( 3.0f - cos2delta );
  a2 = ah2 * ( 3.0f - cos2delta );
  a3 = ah3 * sin2delta;
  a4 = ah4 * sin2delta;
  a5 = ah5 * cos1delta * cos1delta;

  /* coefficients for b(t) */
  b1 = bh1 * sin1delta;
  b2 = bh2 * sin1delta;
  b3 = bh3 * cos1delta;
  b4 = bh4 * cos1delta;


  /*---------- Compute the a(t_i) and b(t_i) ---------- */
  coeffs->A = 0;
  coeffs->B = 0;
  coeffs->C = 0;
  coeffs->D = 0;
  for ( i=0; i < numSteps; i++ )
    {
      REAL4 ah;
      REAL4 cos1ah, sin1ah, cos2ah, sin2ah;
      REAL4 ai, bi;

      ah = alpha - DetectorStates->data[i].LMST;

      sin_cos_LUT ( &sin1ah, &cos1ah, ah );
      sin2ah = 2.0f * sin1ah * cos1ah;
      cos2ah = cos1ah * cos1ah - sin1ah * sin1ah;

      ai = a1 * cos2ah + a2 * sin2ah + a3 * cos1ah + a4 * sin1ah + a5;
      bi = b1 * cos2ah + b2 * sin2ah + b3 * cos1ah + b4 * sin1ah;
      coeffs->a->data[i] = ai;
      coeffs->b->data[i] = bi;

      /* sum A, B, C on the fly */
      coeffs->A += ai * ai;
      coeffs->B += bi * bi;
      coeffs->C += ai * bi;

    } /* for i < numSteps */

  /* finish calculation of A,B,C, D */
  norm = 2.0f / numSteps;
  coeffs->A *= norm;
  coeffs->B *= norm;
  coeffs->C *= norm;

  coeffs->D = coeffs->A * coeffs->B - coeffs->C * coeffs->C;

  RETURN(status);

} /* LALGetAMCoeffs() */

/** Compute the 'amplitude coefficients' \f$a(t)\sin\zeta\f$,
 * \f$b(t)\sin\zeta\f$ as defined in \ref JKS98 for a series of
 * timestamps.
 * 
 * The input consists of the DetectorState-timeseries, which contains
 * the detector-info and the LMST's corresponding to the different times.
 * 
 * In order to allow re-using the output-structure AMCoeffs for subsequent
 * calls, we require the REAL4Vectors a and b to be allocated already and 
 * to have the same length as the DetectoStates-timeseries.
 *
 * \note This is an alternative implementation to both LALComputeAM()
 * and LALGetAMCoeffs(), which uses the geometrical definition of
 * \f$a\sin\zeta\f$ and \f$b\sin\zeta\f$ as detector response
 * coefficients in a preferred polarization basis.  (It is thereby
 * more general than the JKS expressions and could be used e.g., with
 * the response tensor of a bar detector with no further modification
 * needed.)  It is not currently in use.
 */ 
void
LALNewGetAMCoeffs(LALStatus *status,
	       AMCoeffs *coeffs,			/**< [out] amplitude-coeffs {a(t_i), b(t_i)} */
	       const DetectorStateSeries *DetectorStates,	/**< timeseries of detector states */
	       SkyPosition skypos			/**< {alpha,delta} of the source */
	       )
{
  REAL4 delta, alpha;
  REAL4 sin1delta, cos1delta;
  REAL4 sin1alpha, cos1alpha;

  REAL4 xi1, xi2;
  REAL4 eta1, eta2, eta3;
  REAL4 norm;
  UINT4 i, numSteps;

  INITSTATUS (status, "LALNewGetAMCoeffs", COMPUTEFSTATC);

  /*---------- check input ---------- */
  ASSERT ( DetectorStates, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);

  numSteps = DetectorStates->length;

  /* require the coeffients-vectors to be allocated and consistent with timestamps */
  ASSERT ( coeffs, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( coeffs->a && coeffs->b, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( (coeffs->a->length == numSteps) && (coeffs->b->length == numSteps), status,
	   COMPUTEFSTATC_EINPUT,  COMPUTEFSTATC_MSGEINPUT);

  /* require sky-pos to be in equatorial coordinates */
  ASSERT ( skypos.system == COORDINATESYSTEM_EQUATORIAL, status, 
	   SKYCOORDINATESH_ESYS, SKYCOORDINATESH_MSGESYS );

  /*---------- We write components of xi and eta vectors in SSB-fixed coords */
  alpha = skypos.longitude;
  delta = skypos.latitude;

  sin_cos_LUT (&sin1delta, &cos1delta, delta );
  sin_cos_LUT (&sin1alpha, &cos1alpha, alpha );
  xi1 = - sin1alpha;
  xi2 =  cos1alpha;
  eta1 = sin1delta * cos1alpha;
  eta2 = sin1delta * sin1alpha;
  eta3 = - cos1delta;

  /*---------- Compute the a(t_i) and b(t_i) ---------- */
  coeffs->A = 0;
  coeffs->B = 0;
  coeffs->C = 0;
  coeffs->D = 0;
  for ( i=0; i < numSteps; i++ )
    {
      REAL4 ai, bi;

      DetectorTensor *d = &(DetectorStates->data[i].detT);
      
      ai =    d->d11 * ( xi1 * xi1 - eta1 * eta1 )
	+ 2 * d->d12 * ( xi1*xi2 - eta1*eta2 )
	- 2 * d->d13 *             eta1 * eta3
	+     d->d22 * ( xi2*xi2 - eta2*eta2 )
	- 2 * d->d23 *             eta2 * eta3
	-     d->d33 *             eta3*eta3;

      bi =    d->d11 * 2 * xi1 * eta1
	+ 2 * d->d12 *   ( xi1 * eta2 + xi2 * eta1 )
	+ 2 * d->d13 *     xi1 * eta3
	+     d->d22 * 2 * xi2 * eta2
	+ 2 * d->d23 *     xi2 * eta3;

      /*
      printf("xi = (%f,%f)\n",xi1,xi2);
      printf("eta = (%f,%f,%f)\n",eta1,eta2,eta3);
      printf("d = (%f %f %f\n",d->d11,d->d12,d->d13);
      printf("     %f %f %f\n",d->d12,d->d22,d->d23);
      printf("     %f %f %f)\n",d->d13,d->d23,d->d33);
      */

      coeffs->a->data[i] = ai;
      coeffs->b->data[i] = bi;

      /* sum A, B, C on the fly */
      coeffs->A += ai * ai;
      coeffs->B += bi * bi;
      coeffs->C += ai * bi;

    } /* for i < numSteps */

  /* finish calculation of A,B,C, D */
  norm = 2.0f / numSteps;
  coeffs->A *= norm;
  coeffs->B *= norm;
  coeffs->C *= norm;

  coeffs->D = coeffs->A * coeffs->B - coeffs->C * coeffs->C;

  RETURN(status);

} /* LALNewGetAMCoeffs() */

/** For a given DetectorStateSeries, calculate the time-differences
 *  \f$\Delta T_\alpha\equiv T(t_\alpha) - T_0\f$, and their
 *  derivatives \f$Tdot_\alpha \equiv d T / d t (t_\alpha)\f$.
 * 
 *  \note The return-vectors \a DeltaT and \a Tdot must be allocated already
 *  and have the same length as the input time-series \a DetStates.
 *
 */
void
LALGetSSBtimes (LALStatus *status, 
		SSBtimes *tSSB,			/**< [out] DeltaT_alpha = T(t_alpha) - T_0; and Tdot(t_alpha) */
		const DetectorStateSeries *DetectorStates,/**< [in] detector-states at timestamps t_i */
		SkyPosition pos,		/**< source sky-location */
		LIGOTimeGPS refTime,		/**< SSB reference-time T_0 of pulsar-parameters */
		SSBprecision precision		/**< relativistic or Newtonian SSB transformation? */
		)
{
  UINT4 numSteps, i;
  REAL8 vn[3];		/* unit-vector pointing to source in Cart. coord. */
  REAL8 alpha, delta;	/* source position */
  REAL8 refTimeREAL8;

  INITSTATUS( status, "LALGetSSBtimes", COMPUTEFSTATC);
  ATTATCHSTATUSPTR (status);

  ASSERT (DetectorStates, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  numSteps = DetectorStates->length;		/* number of timestamps */

  ASSERT (tSSB, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (tSSB->DeltaT, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (tSSB->Tdot, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);

  ASSERT (tSSB->DeltaT->length == numSteps, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
  ASSERT (tSSB->Tdot->length == numSteps, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);

  ASSERT (precision < SSBPREC_LAST, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
  ASSERT ( pos.system == COORDINATESYSTEM_EQUATORIAL, status,
	   COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
  

  /* convenience variables */
  alpha = pos.longitude;
  delta = pos.latitude;
  refTimeREAL8 = GPS2REAL8(refTime);

  /*----- now calculate the SSB transformation in the precision required */
  switch (precision)
    {
    case SSBPREC_NEWTONIAN:	/* use simple vr.vn to calculate time-delay */

      /*----- get the cartesian source unit-vector */
      vn[0] = cos(alpha) * cos(delta);
      vn[1] = sin(alpha) * cos(delta);
      vn[2] = sin(delta);

      for (i=0; i < numSteps; i++ )
	{
	  LIGOTimeGPS *ti = &(DetectorStates->data[i].tGPS);
	  /* DeltaT_alpha */
	  tSSB->DeltaT->data[i]  = GPS2REAL8 ( (*ti) );
	  tSSB->DeltaT->data[i] += SCALAR(vn, DetectorStates->data[i].rDetector);
	  tSSB->DeltaT->data[i] -= refTimeREAL8;

	  /* Tdot_alpha */
	  tSSB->Tdot->data[i] = 1.0 + SCALAR(vn, DetectorStates->data[i].vDetector);
	  
	} /* for i < numSteps */

      break;

    case SSBPREC_RELATIVISTIC:	/* use LALBarycenter() to get SSB-times and derivative */
      for (i=0; i < numSteps; i++ )
	{
	  BarycenterInput baryinput = empty_BarycenterInput;
	  EmissionTime emit;
	  DetectorState *state = &(DetectorStates->data[i]);

	  baryinput.tgps = state->tGPS;
	  baryinput.site = DetectorStates->detector;
	  /* ARGHHH!!! */
	  baryinput.site.location[0] /= LAL_C_SI;
	  baryinput.site.location[1] /= LAL_C_SI;
	  baryinput.site.location[2] /= LAL_C_SI;

	  baryinput.alpha = alpha;
	  baryinput.delta = delta;
	  baryinput.dInv = 0;

	  TRY ( LALBarycenter(status->statusPtr, &emit, &baryinput, &(state->earthState)), status);

	  tSSB->DeltaT->data[i] = GPS2REAL8 ( emit.te ) - refTimeREAL8;
	  tSSB->Tdot->data[i] = emit.tDot;

	} /* for i < numSteps */

      break;
    default:
      LALPrintError ("\n?? Something went wrong.. this should never be called!\n\n");
      ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
      break;
    } /* switch precision */

  /* finally: store the reference-time used into the output-structure */
  tSSB->refTime = refTime;

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALGetSSBtimes() */

/** Multi-IFO version of LALGetSSBtimes(). 
 * Get all SSB-timings for all input detector-series.
 *
 * NOTE: contrary to LALGetSSBtimes(), this functions *allocates* the output-vector,
 * use XLALDestroyMultiSSBtimes() to free this.
 */
void
LALGetMultiSSBtimes (LALStatus *status, 
		     MultiSSBtimes **multiSSB,		/**< [out] SSB-timings for all input detector-state series */
		     const MultiDetectorStateSeries *multiDetStates, /**< [in] detector-states at timestamps t_i */
		     SkyPosition skypos,		/**< source sky-position [in equatorial coords!] */
		     LIGOTimeGPS refTime,		/**< SSB reference-time T_0 for SSB-timing */
		     SSBprecision precision		/**< use relativistic or Newtonian SSB timing?  */
		     )
{
  UINT4 X, numDetectors;
  MultiSSBtimes *ret = NULL;

  INITSTATUS( status, "LALGetMultiSSBtimes", COMPUTEFSTATC);
  ATTATCHSTATUSPTR (status);

  /* check input */
  ASSERT (multiDetStates, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (multiDetStates->length, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (multiSSB, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( *multiSSB == NULL, status,COMPUTEFSTATC_ENONULL, COMPUTEFSTATC_MSGENONULL);
  ASSERT ( skypos.system == COORDINATESYSTEM_EQUATORIAL, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );

  numDetectors = multiDetStates->length;

  if ( ( ret = LALCalloc( 1, sizeof( *ret ) )) == NULL ) {
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);    
  }
  ret->length = numDetectors;
  if ( ( ret->data = LALCalloc ( numDetectors, sizeof ( *ret->data ) )) == NULL ) {
    LALFree ( ret );
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }

  for ( X=0; X < numDetectors; X ++ )
    {
      SSBtimes *SSBtimesX = NULL;
      UINT4 numStepsX = multiDetStates->data[X]->length;

      ret->data[X] = LALCalloc ( 1, sizeof ( *(ret->data[X]) ) );
      SSBtimesX = ret->data[X];
      SSBtimesX->DeltaT = XLALCreateREAL8Vector ( numStepsX );
      if ( (SSBtimesX->Tdot = XLALCreateREAL8Vector ( numStepsX )) == NULL ) {
	LALPrintError ("\nOut of memory!\n\n");
	goto failed;
      }

      LALGetSSBtimes (status->statusPtr, SSBtimesX, multiDetStates->data[X], skypos, refTime, precision );
      if ( status->statusPtr->statusCode ) 
	{
	  LALPrintError ( "\nCall to LALGetSSBtimes() has failed ... \n\n");
	  goto failed;
	}
 
    } /* for X < numDet */

  goto success;

 failed:
  /* free all memory allocated so far */
  XLALDestroyMultiSSBtimes ( ret );
  ABORT ( status, -1, "LALGetMultiSSBtimes failed" );

 success:
  (*multiSSB) = ret;

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALGetMultiSSBtimes() */

/** Multi-IFO version of LALGetAMCoeffs(). 
 * Get all antenna-pattern coefficients for all input detector-series.
 *
 * NOTE: contrary to LALGetAMCoeffs(), this functions *allocates* the output-vector,
 * use XLALDestroyMultiAMCoeffs() to free this.
 */
void
LALGetMultiAMCoeffs (LALStatus *status, 
		     MultiAMCoeffs **multiAMcoef,	/**< [out] AM-coefficients for all input detector-state series */
		     const MultiDetectorStateSeries *multiDetStates, /**< [in] detector-states at timestamps t_i */
		     SkyPosition skypos			/**< source sky-position [in equatorial coords!] */
		     )
{
  UINT4 X, numDetectors;
  MultiAMCoeffs *ret = NULL;

  INITSTATUS( status, "LALGetMultiAMCoeffs", COMPUTEFSTATC);
  ATTATCHSTATUSPTR (status);

  /* check input */
  ASSERT (multiDetStates, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (multiDetStates->length, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (multiAMcoef, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( *multiAMcoef == NULL, status,COMPUTEFSTATC_ENONULL, COMPUTEFSTATC_MSGENONULL);
  ASSERT ( skypos.system == COORDINATESYSTEM_EQUATORIAL, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );

  numDetectors = multiDetStates->length;

  if ( ( ret = LALCalloc( 1, sizeof( *ret ) )) == NULL ) {
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);    
  }
  ret->length = numDetectors;
  if ( ( ret->data = LALCalloc ( numDetectors, sizeof ( *ret->data ) )) == NULL ) {
    LALFree ( ret );
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }

  for ( X=0; X < numDetectors; X ++ )
    {
      AMCoeffs *amcoeX = NULL;
      UINT4 numStepsX = multiDetStates->data[X]->length;

      ret->data[X] = LALCalloc ( 1, sizeof ( *(ret->data[X]) ) );
      amcoeX = ret->data[X];
      amcoeX->a = XLALCreateREAL4Vector ( numStepsX );
      if ( (amcoeX->b = XLALCreateREAL4Vector ( numStepsX )) == NULL ) {
	LALPrintError ("\nOut of memory!\n\n");
	goto failed;
      }

      /* LALGetAMCoeffs (status->statusPtr, amcoeX, multiDetStates->data[X], skypos ); */
      LALNewGetAMCoeffs (status->statusPtr, amcoeX, multiDetStates->data[X], skypos );
      if ( status->statusPtr->statusCode ) 
	{
	  LALPrintError ( "\nCall to LALNewGetAMCoeffs() has failed ... \n\n");
	  goto failed;
	}
 
    } /* for X < numDetectors */

  goto success;

 failed:
  /* free all memory allocated so far */
  XLALDestroyMultiAMCoeffs ( ret );
  ABORT ( status, -1, "LALGetMultiAMCoeffs() failed" );

 success:
  (*multiAMcoef) = ret;

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALGetMultiAMCoeffs() */



/* ===== Object creation/destruction functions ===== */

/** Destroy a MultiSSBtimes structure. 
 * Note, this is "NULL-robust" in the sense that it will not crash 
 * on NULL-entries anywhere in this struct, so it can be used
 * for failure-cleanup even on incomplete structs 
 */
void
XLALDestroyMultiSSBtimes ( MultiSSBtimes *multiSSB )
{
  UINT4 X;
  SSBtimes *tmp;

  if ( ! multiSSB )
    return;

  if ( multiSSB->data )
    {
      for ( X=0; X < multiSSB->length; X ++ ) 
	{
	  if ( (tmp = multiSSB->data[X]) != NULL )
	    {
	      if ( tmp->DeltaT )
		XLALDestroyREAL8Vector ( tmp->DeltaT );
	      if ( tmp->Tdot )
		XLALDestroyREAL8Vector ( tmp->Tdot );
	      LALFree ( tmp );
	    } /* if multiSSB->data[X] */
	} /* for X < numDetectors */
      LALFree ( multiSSB->data );
    }
  LALFree ( multiSSB );

  return;

} /* XLALDestroyMultiSSBtimes() */

/** Destroy a MultiAMCoeffs structure. 
 * Note, this is "NULL-robust" in the sense that it will not crash 
 * on NULL-entries anywhere in this struct, so it can be used
 * for failure-cleanup even on incomplete structs 
 */
void
XLALDestroyMultiAMCoeffs ( MultiAMCoeffs *multiAMcoef )
{
  UINT4 X;
  AMCoeffs *tmp;

  if ( ! multiAMcoef )
    return;

  if ( multiAMcoef->data )
    {
      for ( X=0; X < multiAMcoef->length; X ++ ) 
	{
	  if ( (tmp = multiAMcoef->data[X]) != NULL )
	    {
	      if ( tmp->a )
		XLALDestroyREAL4Vector ( tmp->a );
	      if ( tmp->b )
		XLALDestroyREAL4Vector ( tmp->b );
	      LALFree ( tmp );
	    } /* if multiAMcoef->data[X] */
	} /* for X < numDetectors */
      LALFree ( multiAMcoef->data );
    }
  LALFree ( multiAMcoef );

  return;

} /* XLALDestroyMultiAMCoeffs() */


/** Destruction of a ComputeFBuffer *contents*, 
 * i.e. the multiSSB and multiAMcoeff, while the 
 * buffer-container is not freed (which is why it's passed
 * by value and not by reference...) */
void
XLALEmptyComputeFBuffer ( ComputeFBuffer cfb )
{
  XLALDestroyMultiSSBtimes ( cfb.multiSSB );
  cfb.multiSSB = NULL;
  XLALDestroyMultiAMCoeffs ( cfb.multiAMcoef );
  cfb.multiAMcoef = NULL;

  return;
} /* XLALDestroyComputeFBuffer() */


/**< Multiply AM-coeffs \f$a_{X\alpha}, b_{X\alpha}\f$ by weights \f$\sqrt(w_{X\alpha})\f$ and 
 * compute the resulting \f$A_d, B_d, C_d\f$ by simply *SUMMING* them, i.e.
 * \f$A_d \equiv \sum_{X,\alpha} \sqrt{w_{X\alpha} a_{X\alpha}^2\f$ etc.
 * 
 * NOTE: this function modifies the AMCoeffs *in place* !
 * NOTE2: if the weights = NULL, we assume unit-weights.
 */
int
XLALWeighMultiAMCoeffs (  MultiAMCoeffs *multiAMcoef, const MultiNoiseWeights *multiWeights )
{
  UINT4 numDetectors, X;
  REAL8 Ad, Bd, Cd, Sinv_Tsft;
  UINT4 alpha;

  if ( !multiAMcoef )
    XLAL_ERROR( "XLALWeighMultiAMCoeffs", XLAL_EINVAL );

  numDetectors = multiAMcoef->length;

  if ( multiWeights && ( multiWeights->length != numDetectors ) )
    {
      LALPrintError("\nmultiWeights must have same length as mulitAMcoef!\n\n");
      XLAL_ERROR( "XLALWeighMultiAMCoeffs", XLAL_EINVAL );
    }
  
  /* noise-weight Antenna-patterns and compute A,B,C */
  Ad = Bd = Cd = 0;

  if ( multiWeights  )
    {
      for ( X=0; X < numDetectors; X ++)
	{
	  AMCoeffs *amcoeX = multiAMcoef->data[X];
	  UINT4 numSteps = amcoeX->a->length;
	  
	  REAL8Vector *weightsX = multiWeights->data[X];;
	  if ( weightsX->length != numSteps ) 
	    {
	      LALPrintError("\nmultiWeights must have same length as mulitAMcoef!\n\n");
	      XLAL_ERROR( "XLALWeighMultiAMCoeffs", XLAL_EINVAL );
	    }
	  
	  for(alpha = 0; alpha < numSteps; alpha++)
	    {
	      REAL8 Sqwi = sqrt ( weightsX->data[alpha] );
	      REAL8 ahat = Sqwi * amcoeX->a->data[alpha] ;
	      REAL8 bhat = Sqwi * amcoeX->b->data[alpha] ;
	      
	      /* *replace* original a(t), b(t) by noise-weighed version! */
	      amcoeX->a->data[alpha] = ahat;
	      amcoeX->b->data[alpha] = bhat;
	      
	      /* sum A, B, C on the fly */
	      Ad += ahat * ahat;
	      Bd += bhat * bhat;
	      Cd += ahat * bhat;
	    } /* for alpha < numSFTsX */
	} /* for X < numDetectors */
      Sinv_Tsft = multiWeights->Sinv_Tsft;
    }
  else /* if no noise-weights: simply add to get A,B,C */
    {
      for ( X=0; X < numDetectors; X ++)
	{
	  AMCoeffs *amcoeX = multiAMcoef->data[X];
	  UINT4 numSteps = amcoeX->a->length;
	  
	  for(alpha = 0; alpha < numSteps; alpha++)
	    {
	      REAL8 ahat = amcoeX->a->data[alpha] ;
	    REAL8 bhat = amcoeX->b->data[alpha] ;
	    
	    /* sum A, B, C on the fly */
	    Ad += ahat * ahat;
	    Bd += bhat * bhat;
	    Cd += ahat * bhat;
	    } /* for alpha < numSFTsX */
	} /* for X < numDetectors */

    } /* if multiWeights == NULL */

  multiAMcoef->Mmunu.Ad = Ad;
  multiAMcoef->Mmunu.Bd = Bd;
  multiAMcoef->Mmunu.Cd = Cd;
  multiAMcoef->Mmunu.Sinv_Tsft = Sinv_Tsft;

  return XLAL_SUCCESS;

} /* XLALWeighMultiAMCoefs() */
  

/* ===== General internal helper functions ===== */

/** Calculate sin(x) and cos(x) to roughly 1e-7 precision using 
 * a lookup-table and Taylor-expansion.
 *
 * NOTE: this function will fail for arguments larger than
 * |x| > INT4_MAX = 2147483647 ~ 2e9 !!!
 *
 * return = 0: OK, nonzero=ERROR
 */
int
sin_cos_LUT (REAL4 *sinx, REAL4 *cosx, REAL8 x)
{
  return sin_cos_2PI_LUT ( sinx, cosx, x * OOTWOPI );
} /* sin_cos_LUT() */

#define LUT_RES         64      /* resolution of lookup-table */
#define LUT_RES_F	(1.0 * LUT_RES)
#define OO_LUT_RES	(1.0 / LUT_RES)

#define X_TO_IND	(1.0 * LUT_RES * OOTWOPI )
#define IND_TO_X	(LAL_TWOPI * OO_LUT_RES)
int
sin_cos_2PI_LUT (REAL4 *sin2pix, REAL4 *cos2pix, REAL8 x)
{
  REAL8 xt;
  INT4 i0;
  REAL8 d, d2;
  REAL8 ts, tc;
  REAL8 dummy;

  static BOOLEAN firstCall = TRUE;
  static REAL4 sinVal[LUT_RES+1], cosVal[LUT_RES+1];

  /* the first time we get called, we set up the lookup-table */
  if ( firstCall )
    {
      UINT4 k;
      for (k=0; k <= LUT_RES; k++)
        {
          sinVal[k] = sin( LAL_TWOPI * k * OO_LUT_RES );
          cosVal[k] = cos( LAL_TWOPI * k * OO_LUT_RES );
        }
      firstCall = FALSE;
    }

  /* we only need the fractional part of 'x', which is number of cylces,
   * this was previously done using
   *   xt = x - (INT4)x;
   * which is numerically unsafe for x > LAL_INT4_MAX ~ 2e9
   * for saftey we therefore rather use modf(), even if that 
   * will be somewhat slower... 
   */
  xt = modf(x, &dummy);/* xt in (-1, 1) */

  if ( xt < 0.0 )
    xt += 1.0;			/* xt in [0, 1 ) */
#ifndef LAL_NDEBUG
  if ( xt < 0.0 || xt > 1.0 )
    {
      LALPrintError("\nFailed numerica in sin_cos_2PI_LUT(): xt = %f not in [0,1)\n\n", xt );
      return XLAL_FAILURE;
    }
#endif

  i0 = (INT4)( xt * LUT_RES_F + 0.5 );	/* i0 in [0, LUT_RES ] */
  d = d2 = LAL_TWOPI * (xt - OO_LUT_RES * i0);
  d2 *= 0.5 * d;

  ts = sinVal[i0];
  tc = cosVal[i0];
   
  /* use Taylor-expansions for sin/cos around LUT-points */
  (*sin2pix) = ts + d * tc - d2 * ts;
  (*cos2pix) = tc - d * ts - d2 * tc;

  return XLAL_SUCCESS;
} /* sin_cos_2PI_LUT() */



/** Parameter-estimation: based on large parts on Yousuke's notes and implemention (in CFSv1),
 * extended for error-estimation.
 * This implementation follows closely the derivations found in 
 * http://www.lsc-group.phys.uwm.edu/cgi-bin/enote.pl?nb=puls5knownpulsardemod&action=view&page=12
 */
void
LALEstimatePulsarAmplitudeParams (LALStatus * status,
				  PulsarCandidate *pulsarParams,  	/**< [out] estimated params {h0,cosi,phi0,psi} plus error-estimates */
				  const Fcomponents *Fstat,	 	/**<  Fstat-components Fa, Fb */
				  const LIGOTimeGPS *FstatRefTime,	/**<  reference-time for the phase of Fa, Fb */
				  const AntennaPatternMatrix *Mmunu 	/**<  antenna-pattern A,B,C and normalization S_inv*Tsft */
				  )
{
  REAL8 A1h, A2h, A3h, A4h;
  REAL8 Ad, Bd, Cd, Dd;
  REAL8 normAmu;
  REAL8 A1check, A2check, A3check, A4check;

  REAL8 Asq, Da, disc;
  REAL8 aPlus, aCross;
  REAL8 Ap2, Ac2;
  REAL8 beta;
  REAL8 phi0, psi;
  REAL8 b1, b2, b3;
  REAL8 h0, cosi;
  
  REAL8 cosphi0, sinphi0, cos2psi, sin2psi;

  REAL8 tolerance = LAL_REAL4_EPS;

  gsl_vector *x_mu, *A_Mu;
  gsl_matrix *M_Mu_Nu;
  gsl_matrix *Jh_Mu_nu;
  gsl_permutation *permh;
  gsl_matrix *tmp, *tmp2;
  int signum;

  INITSTATUS (status, "LALEstimatePulsarAmplitudeParams", COMPUTEFSTATC );
  ATTATCHSTATUSPTR (status);

  ASSERT ( pulsarParams, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( Fstat, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( FstatRefTime, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( Mmunu, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );

  Ad = Mmunu->Ad;
  Bd = Mmunu->Bd;
  Cd = Mmunu->Cd;
  Dd = Ad * Bd - Cd * Cd;

  normAmu = 2.0 / sqrt(Mmunu->Sinv_Tsft);	/* generally *very* small!! */

  /* ----- GSL memory allocation ----- */
  if ( ( x_mu = gsl_vector_calloc (4) ) == NULL ) {
    ABORT ( status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM );
  }
  if ( ( A_Mu = gsl_vector_calloc (4) ) == NULL ) {
    ABORT ( status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM );
  }
  if ( ( M_Mu_Nu = gsl_matrix_calloc (4, 4) ) == NULL ) {
    ABORT ( status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM );
  }
  if ( ( Jh_Mu_nu = gsl_matrix_calloc (4, 4) ) == NULL ) {
    ABORT ( status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM );
  }

  if ( ( permh = gsl_permutation_calloc ( 4 )) == NULL ) {
    ABORT ( status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM );
  }
  if ( ( tmp = gsl_matrix_calloc (4, 4) ) == NULL ) {
    ABORT ( status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM );
  }
  if ( ( tmp2 = gsl_matrix_calloc (4, 4) ) == NULL ) {
    ABORT ( status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM );
  }

  /* ----- fill vector x_mu */
  gsl_vector_set (x_mu, 0,   Fstat->Fa.re );	/* x_1 */
  gsl_vector_set (x_mu, 1,   Fstat->Fb.re ); 	/* x_2 */
  gsl_vector_set (x_mu, 2, - Fstat->Fa.im );	/* x_3 */
  gsl_vector_set (x_mu, 3, - Fstat->Fb.im );	/* x_4 */

  /* ----- fill matrix M^{mu,nu} [symmetric: use UPPER HALF ONLY!!]*/
  gsl_matrix_set (M_Mu_Nu, 0, 0,   Bd / Dd );
  gsl_matrix_set (M_Mu_Nu, 1, 1,   Ad / Dd );
  gsl_matrix_set (M_Mu_Nu, 0, 1, - Cd / Dd );  

  gsl_matrix_set (M_Mu_Nu, 2, 2,   Bd / Dd );
  gsl_matrix_set (M_Mu_Nu, 3, 3,   Ad / Dd );
  gsl_matrix_set (M_Mu_Nu, 2, 3, - Cd / Dd );  

  /* get (un-normalized) MLE's for amplitudes A^mu  = M^{mu,nu} x_nu */

  /* GSL-doc: int gsl_blas_dsymv (CBLAS_UPLO_t Uplo, double alpha, const gsl_matrix * A, 
   *                              const gsl_vector * x, double beta, gsl_vector * y )
   * 
   * compute the matrix-vector product and sum: y = alpha A x + beta y 
   * for the symmetric matrix A. Since the matrix A is symmetric only its 
   * upper half or lower half need to be stored. When Uplo is CblasUpper 
   * then the upper triangle and diagonal of A are used, and when Uplo 
   * is CblasLower then the lower triangle and diagonal of A are used. 
   */
  TRYGSL(gsl_blas_dsymv (CblasUpper, 1.0, M_Mu_Nu, x_mu, 0.0, A_Mu), status);

  A1h = gsl_vector_get ( A_Mu, 0 );
  A2h = gsl_vector_get ( A_Mu, 1 );
  A3h = gsl_vector_get ( A_Mu, 2 );
  A4h = gsl_vector_get ( A_Mu, 3 );

  /* LogPrintf (LOG_DEBUG, "norm= %g; A1 = %g, A2 = %g, A3 = %g, A4 = %g\n",  normAmu, A1h, A2h, A3h, A4h ); 
   */
  Asq = SQ(A1h) + SQ(A2h) + SQ(A3h) + SQ(A4h);
  Da = A1h * A4h - A2h * A3h;
  disc = sqrt ( SQ(Asq) - 4.0 * SQ(Da) );

  Ap2  = 0.5 * ( Asq + disc );
  aPlus = sqrt(Ap2);		/* not yet normalized */

  Ac2 = 0.5 * ( Asq - disc );
  aCross = sqrt( Ac2 );	
  aCross *= MYSIGN ( Da ); 	/* not yet normalized */

  beta = aCross / aPlus;
  
  b1 =   A4h - beta * A1h;
  b2 =   A3h + beta * A2h;
  b3 = - A1h + beta * A4h ;

  psi  = 0.5 * atan ( b1 /  b2 );	/* in [-pi/4,pi/4] (gauge used also by TDS) */
  phi0 =       atan ( b2 / b3 );	/* in [-pi/2,pi/2] */

  /* Fix remaining sign-ambiguity by checking sign of reconstructed A1 */
  A1check = aPlus * cos(phi0) * cos(2.0*psi) - aCross * sin(phi0) * sin(2*psi);
  if ( A1check * A1h <  0 )
    phi0 += LAL_PI;

  cosphi0 = cos(phi0);
  sinphi0 = sin(phi0);
  cos2psi = cos(2*psi);
  sin2psi = sin(2*psi);

  /* check numerical consistency of estimated Amu and reconstructed */
  A1check =   aPlus * cosphi0 * cos2psi - aCross * sinphi0 * sin2psi;  
  A2check =   aPlus * cosphi0 * sin2psi + aCross * sinphi0 * cos2psi;  
  A3check = - aPlus * sinphi0 * cos2psi - aCross * cosphi0 * sin2psi;  
  A4check = - aPlus * sinphi0 * sin2psi + aCross * cosphi0 * cos2psi;  

  /* LogPrintf (LOG_DEBUG, "reconstructed:    A1 = %g, A2 = %g, A3 = %g, A4 = %g\n", 
     A1check, A2check, A3check, A4check ); 
  */

  if ( ( fabs( (A1check - A1h)/A1h ) > tolerance ) ||
       ( fabs( (A2check - A2h)/A2h ) > tolerance ) ||
       ( fabs( (A3check - A3h)/A3h ) > tolerance ) ||
       ( fabs( (A4check - A4h)/A4h ) > tolerance ) )
    {
      if ( lalDebugLevel )
	LALPrintError ( "WARNING LALEstimatePulsarAmplitudeParams(): Difference between estimated and reconstructed Amu exceeds tolerance of %g\n",  
			tolerance );
    }

  /* translate A_{+,x} into {h_0, cosi} */
  h0 = aPlus + sqrt ( disc );  /* not yet normalized ! */
  cosi = aCross / h0;


  /* ========== Estimate the errors ========== */
  
  /* ----- compute derivatives \partial A^\mu / \partial B^\nu, where
   * we consider the output-variables B^\nu = (h0, cosi, phi0, psi) 
   * where aPlus = 0.5 * h0 * (1 + cosi^2)  and aCross = h0 * cosi
   */
  { /* Ahat^mu is defined as A^mu with the replacements: A_+ --> A_x, and A_x --> h0 */
    REAL8 A1hat =   aCross * cosphi0 * cos2psi - h0 * sinphi0 * sin2psi;  
    REAL8 A2hat =   aCross * cosphi0 * sin2psi + h0 * sinphi0 * cos2psi;  
    REAL8 A3hat = - aCross * sinphi0 * cos2psi - h0 * cosphi0 * sin2psi;  
    REAL8 A4hat = - aCross * sinphi0 * sin2psi + h0 * cosphi0 * cos2psi;  

    /* ----- A1 =   aPlus * cosphi0 * cos2psi - aCross * sinphi0 * sin2psi; ----- */
    gsl_matrix_set (Jh_Mu_nu, 0, 0,   A1h / h0 );	/* dA1/h0 */
    gsl_matrix_set (Jh_Mu_nu, 0, 1,   A1hat ); 		/* dA1/dcosi */
    gsl_matrix_set (Jh_Mu_nu, 0, 2,   A3h );		/* dA1/dphi0 */
    gsl_matrix_set (Jh_Mu_nu, 0, 3, - 2.0 * A2h );	/* dA1/dpsi */
  
    /* ----- A2 =   aPlus * cosphi0 * sin2psi + aCross * sinphi0 * cos2psi; ----- */ 
    gsl_matrix_set (Jh_Mu_nu, 1, 0,   A2h / h0 );	/* dA2/h0 */
    gsl_matrix_set (Jh_Mu_nu, 1, 1,   A2hat ); 		/* dA2/dcosi */
    gsl_matrix_set (Jh_Mu_nu, 1, 2,   A4h );		/* dA2/dphi0 */
    gsl_matrix_set (Jh_Mu_nu, 1, 3,   2.0 * A1h );	/* dA2/dpsi */

    /* ----- A3 = - aPlus * sinphi0 * cos2psi - aCross * cosphi0 * sin2psi; ----- */ 
    gsl_matrix_set (Jh_Mu_nu, 2, 0,   A3h / h0 );	/* dA3/h0 */ 
    gsl_matrix_set (Jh_Mu_nu, 2, 1,   A3hat ); 		/* dA3/dcosi */
    gsl_matrix_set (Jh_Mu_nu, 2, 2, - A1h );		/* dA3/dphi0 */  
    gsl_matrix_set (Jh_Mu_nu, 2, 3, - 2.0 * A4h );	/* dA3/dpsi */

    /* ----- A4 = - aPlus * sinphi0 * sin2psi + aCross * cosphi0 * cos2psi; ----- */ 
    gsl_matrix_set (Jh_Mu_nu, 3, 0,   A4h / h0 );	/* dA4/h0 */
    gsl_matrix_set (Jh_Mu_nu, 3, 1,   A4hat ); 		/* dA4/dcosi */
    gsl_matrix_set (Jh_Mu_nu, 3, 2, - A2h );		/* dA4/dphi0 */
    gsl_matrix_set (Jh_Mu_nu, 3, 3,   2.0 * A3h );	/* dA4/dpsi */
  }

  /* ----- compute inverse matrices Jh^{-1} by LU-decomposition ----- */
  TRYGSL( gsl_linalg_LU_decomp (Jh_Mu_nu, permh, &signum ), status);

  /* inverse matrix */
  TRYGSL(gsl_linalg_LU_invert (Jh_Mu_nu, permh, tmp ), status);
  gsl_matrix_memcpy ( Jh_Mu_nu, tmp );

  /* ----- compute Jh^-1 . Minv . (Jh^-1)^T ----- */
  
  /* GSL-doc: gsl_blas_dgemm (CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, double alpha, 
   *                          const gsl_matrix *A, const gsl_matrix *B, double beta, gsl_matrix *C)
   * These functions compute the matrix-matrix product and sum 
   * C = \alpha op(A) op(B) + \beta C 
   * where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans 
   * and similarly for the parameter TransB.
   */

  /* first tmp = Minv . (Jh^-1)^T */
  TRYGSL( gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, M_Mu_Nu, Jh_Mu_nu, 0.0, tmp ), status);
  /* then J^-1 . tmp , store result in tmp2 */
  TRYGSL( gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Jh_Mu_nu, tmp, 0.0, tmp2 ), status);
  gsl_matrix_memcpy ( Jh_Mu_nu, tmp2 );
  
  /* ===== debug-output resulting matrices ===== */
  /* if ( lalDebugLevel ) */
  /* printGSLmatrix4 ( stdout, "var(dBh^mu, dBh^nu) = \n", Jh_Mu_nu ); */

  /* propagate initial-phase from Fstat-reference-time to refTime of Doppler-params */
  TRY ( LALExtrapolatePulsarPhase (status->statusPtr, &phi0, pulsarParams->Doppler.fkdot, pulsarParams->Doppler.refTime, phi0, (*FstatRefTime) ),
	status );
  
  if ( phi0 < 0 )	      /* make sure phi0 in [0, 2*pi] */
    phi0 += LAL_TWOPI;
  phi0 = fmod ( phi0, LAL_TWOPI );

  /* fill candidate-struct with the obtained signal-parameters and error-estimations */
  pulsarParams->Amp.h0     = normAmu * h0;
  pulsarParams->Amp.cosi   = cosi;
  pulsarParams->Amp.phi0   = phi0;
  pulsarParams->Amp.psi    = psi;

  /* read out principal estimation-errors from diagonal elements of inverse Fisher-matrix*/
  pulsarParams->dAmp.h0     = normAmu * sqrt( gsl_matrix_get (Jh_Mu_nu, 0, 0 ) );
  pulsarParams->dAmp.cosi   = sqrt( gsl_matrix_get (Jh_Mu_nu, 1, 1 ) );
  pulsarParams->dAmp.phi0   = sqrt( gsl_matrix_get (Jh_Mu_nu, 2, 2 ) );
  pulsarParams->dAmp.psi    = sqrt( gsl_matrix_get (Jh_Mu_nu, 3, 3 ) );

  /* ----- free GSL memory ----- */
  gsl_vector_free ( x_mu );
  gsl_vector_free ( A_Mu );
  gsl_matrix_free ( M_Mu_Nu );
  gsl_matrix_free ( Jh_Mu_nu );
  gsl_permutation_free ( permh );
  gsl_matrix_free ( tmp );
  gsl_matrix_free ( tmp2 );

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LALEstimatePulsarAmplitudeParams() */


/** Output 4x4 gsl_matrix in octave-format.
 * return -1 on error, 0 if OK.
 */
static int
printGSLmatrix4 ( FILE *fp, const CHAR *prefix, const gsl_matrix *gij )
{
  if ( !gij )
    return -1;
  if ( !fp )
    return -1;
  if ( (gij->size1 != 4) || (gij->size2 != 4 ) )
    return -1;
  if ( !prefix ) 
    return -1;

  fprintf (fp, prefix );
  fprintf (fp, " [ %.9f, %.9f, %.9f, %.9f;\n",
	   gsl_matrix_get ( gij, 0, 0 ), gsl_matrix_get ( gij, 0, 1 ), 
	   gsl_matrix_get ( gij, 0, 2 ), gsl_matrix_get ( gij, 0, 3 ) );
  fprintf (fp, "   %.9f, %.9f, %.9f, %.9f;\n",
	   gsl_matrix_get ( gij, 1, 0 ), gsl_matrix_get ( gij, 1, 1 ),
	   gsl_matrix_get ( gij, 1, 2 ), gsl_matrix_get ( gij, 1, 3 ) );
  fprintf (fp, "   %.9f, %.9f, %.9f, %.9f;\n",
	   gsl_matrix_get ( gij, 2, 0 ), gsl_matrix_get ( gij, 2, 1 ),
	   gsl_matrix_get ( gij, 2, 2 ), gsl_matrix_get ( gij, 2, 3 ) );
  fprintf (fp, "   %.9f, %.9f, %.9f, %.9f ];\n",
	   gsl_matrix_get ( gij, 3, 0 ), gsl_matrix_get ( gij, 3, 1 ),
	   gsl_matrix_get ( gij, 3, 2 ), gsl_matrix_get ( gij, 3, 3 ) );
  
  return 0;

} /* printGSLmatrix4() */
