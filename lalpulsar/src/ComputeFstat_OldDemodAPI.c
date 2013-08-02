/*
 * Copyright (C) 2007 Chris Messenger
 * Copyright (C) 2006 John T. Whelan, Badri Krishnan
 * Copyright (C) 2005, 2006, 2007, 2010 Reinhard Prix
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

/*
 * Functions to calculate the so-called F-statistic for a given point in parameter-space,
 * following the equations in \cite JKS98.
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
#include <gsl/gsl_roots.h>

#include <lal/AVFactories.h>
#include "ComputeFstat_OldDemodAPI.h"
#include "ComplexAM.h"

/*---------- local DEFINES ----------*/
#define LD_SMALL4       (2.0e-4)		/**< "small" number for REAL4*/
#define OOTWOPI         (1.0 / LAL_TWOPI)	/**< 1/2pi */

#define TWOPI_FLOAT     6.28318530717958f  	/**< single-precision 2*pi */
#define OOTWOPI_FLOAT   (1.0f / TWOPI_FLOAT)	/**< single-precision 1 / (2pi) */

#define EA_ACC          1E-9                    /* the timing accuracy of LALGetBinaryTimes in seconds */

/*----- Macros ----- */
#define INIT_MEM(x) memset(&(x), 0, sizeof((x)))

#define MYSIGN(x) ( ((x) < 0) ? (-1.0):(+1.0) )

/** Simple Euklidean scalar product for two 3-dim vectors in cartesian coords */
#define SCALAR(u,v) ((u)[0]*(v)[0] + (u)[1]*(v)[1] + (u)[2]*(v)[2])

#define SQ(x) ( (x) * (x) )

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- Global variables ----------*/
#define NUM_FACT 7
static const REAL8 inv_fact[NUM_FACT] = { 1.0, 1.0, (1.0/2.0), (1.0/6.0), (1.0/24.0), (1.0/120.0), (1.0/720.0) };

/* empty initializers  */
const SSBtimes empty_SSBtimes;
const MultiSSBtimes empty_MultiSSBtimes;

const Fcomponents empty_Fcomponents;
const ComputeFParams empty_ComputeFParams;
const ComputeFBuffer empty_ComputeFBuffer;

/*---------- internal prototypes ----------*/

static double gsl_E_solver ( double E, void *p );

struct E_solver_params {
  double p, q, r, fracOrb;
};

static double gsl_E_solver ( REAL8 E, void *par )
{
  struct E_solver_params *params = (struct E_solver_params*) par;
  double p = params->p;
  double q = params->q;
  double r = params->r;
  double x0 = params->fracOrb * LAL_TWOPI;

  /* this is the function relating the observed time since periapse in the SSB to the true eccentric anomoly E */
  double diff = - x0 + ( E + p * sin(E) + q * ( cos(E) - 1.0 ) + r );

  return diff;
} // gsl_E_solver()

/*==================== FUNCTION DEFINITIONS ====================*/


/**
 * Function to compute a vector of Fstatistic values for a number of frequency bins.
 * This function is simply a wrapper for ComputeFstat() which is called repeatedly for
 * every frequency value.  The output, i.e. fstatVector must be properly allocated
 * before this function is called.  The values of the start frequency, the step size
 * in the frequency and the number of frequency values for which the Fstatistic is
 * to be calculated are read from fstatVector.  The other parameters are not checked and
 * they must be correctly set outside this function.
 */
void ComputeFStatFreqBand ( LALStatus *status,				/**< pointer to LALStatus structure */
			    REAL4FrequencySeries *fstatVector, 		/**< [out] Vector of Fstat values */
			    const PulsarDopplerParams *doppler,		/**< parameter-space point to compute F for */
			    const MultiSFTVector *multiSFTs, 		/**< normalized (by DOUBLE-sided Sn!) data-SFTs of all IFOs */
			    const MultiNoiseWeights *multiWeights,	/**< noise-weights of all SFTs */
			    const MultiDetectorStateSeries *multiDetStates,/**< 'trajectories' of the different IFOs */
			    const ComputeFParams *params		/**< addition computational params */
			    )
{
  UINT4 numBins, k;
  REAL8 deltaF, fStart;
  Fcomponents Fstat;
  PulsarDopplerParams thisPoint;
  ComputeFBuffer cfBuffer = empty_ComputeFBuffer;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( multiSFTs, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( doppler, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( multiDetStates, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( params, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );

  ASSERT ( multiDetStates->length == multiSFTs->length, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );
  ASSERT ( fstatVector, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( fstatVector->data, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( fstatVector->data->data, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( fstatVector->data->length > 0, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );

  if ( params->returnAtoms )
    {
      XLALPrintError ("%s: using the option 'returnAtoms' is not supported in this function!\n", __func__ );
      ABORT (status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );
    }

/**
 * something to improve/cleanup -- the start frequency is available both
 * from the fstatvector and from the input doppler point -- they could be inconsistent
 * or the user of this function could misunderstand
 */

  /* a check that the f0 values from thisPoint and fstatVector are
     at least close to each other -- this is only meant to catch
     stupid errors but not subtle ones */
  ASSERT ( fabs(fstatVector->f0 - doppler->fkdot[0]) < fstatVector->deltaF,
	   status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );


  /* copy values from 'doppler' to local variable 'thisPoint' */
  thisPoint = *doppler;

  numBins = fstatVector->data->length;
  deltaF = fstatVector->deltaF;
  fStart = thisPoint.fkdot[0];

  /* loop over frequency values and fill up values in fstatVector */
  for ( k = 0; k < numBins; k++) {

    thisPoint.fkdot[0] = fStart + k*deltaF;

    TRY (ComputeFStat ( status->statusPtr, &Fstat, &thisPoint, multiSFTs, multiWeights,
			multiDetStates, params, &cfBuffer ), status);

    fstatVector->data->data[k] = Fstat.F;

  }

  XLALEmptyComputeFBuffer ( &cfBuffer );

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* ComputeFStatFreqBand() */




/**
 * Function to compute (multi-IFO) F-statistic for given parameter-space point \a doppler,
 * normalized SFT-data (normalized by <em>double-sided</em> PSD Sn), noise-weights
 * and detector state-series
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
ComputeFStat ( LALStatus *status,				/**< pointer to LALStatus structure */
	       Fcomponents *Fstat,                 		/**< [out] Fstatistic + Fa, Fb */
	       const PulsarDopplerParams *doppler, 		/**< parameter-space point to compute F for */
	       const MultiSFTVector *multiSFTs,    		/**< normalized (by DOUBLE-sided Sn!) data-SFTs of all IFOs */
	       const MultiNoiseWeights *multiWeights,		/**< noise-weights of all SFTs */
	       const MultiDetectorStateSeries *multiDetStates,	/**< 'trajectories' of the different IFOs */
	       const ComputeFParams *params,       		/**< addition computational params */
	       ComputeFBuffer *cfBuffer            		/**< CF-internal buffering structure */
	       )
{
  Fcomponents retF = empty_Fcomponents;
  UINT4 X, numDetectors;
  MultiSSBtimes *multiSSB = NULL;
  MultiSSBtimes *multiBinary = NULL;
  const MultiSSBtimes *multiSSBTotal = NULL;
  MultiAMCoeffs *multiAMcoef = NULL;
  MultiCmplxAMCoeffs *multiCmplxAMcoef = NULL;
  REAL8 Ad, Bd, Cd, Dd_inv, Ed;
  SkyPosition skypos;

  INITSTATUS(status);
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

  /* ----- prepare return of 'FstatAtoms' if requested */
  if ( params->returnAtoms )
    {
      if ( (retF.multiFstatAtoms = LALMalloc ( sizeof(*retF.multiFstatAtoms) )) == NULL ){
	ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
      }
      retF.multiFstatAtoms->length = numDetectors;
      if ( (retF.multiFstatAtoms->data = LALMalloc ( numDetectors * sizeof(*retF.multiFstatAtoms->data) )) == NULL ) {
	LALFree ( retF.multiFstatAtoms );
	ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
      }

    } /* if returnAtoms */


  /* ----- check if that skyposition SSB+AMcoef were already buffered */
  if ( cfBuffer
       && ( cfBuffer->multiDetStates == multiDetStates )
       && ( cfBuffer->Alpha == doppler->Alpha )
       && ( cfBuffer->Delta == doppler->Delta )
       && cfBuffer->multiSSB )
    { /* yes ==> reuse */
      multiSSB = cfBuffer->multiSSB;

      /* re-use (LWL) AM coefficients whenever available */
      if ( cfBuffer->multiAMcoef )
	multiAMcoef = cfBuffer->multiAMcoef;

      /* re-use RAA AM coefficients *only* if bufferedRAA is TRUE !*/
      if ( params->bufferedRAA && cfBuffer->multiCmplxAMcoef  )
	multiCmplxAMcoef = cfBuffer->multiCmplxAMcoef;

    } /* if have buffered stuff to reuse */
  else
    {
      skypos.system =   COORDINATESYSTEM_EQUATORIAL;
      skypos.longitude = doppler->Alpha;
      skypos.latitude  = doppler->Delta;
      if ( (multiSSB = XLALGetMultiSSBtimes ( multiDetStates, skypos, doppler->refTime, params->SSBprec )) == NULL )
        {
          XLALPrintError("XLALGetMultiSSBtimes() failed with error = %d\n\n", xlalErrno );
          ABORT ( status, COMPUTEFSTATC_EXLAL, COMPUTEFSTATC_MSGEXLAL );
        }
      if ( cfBuffer )
	{
	  XLALDestroyMultiSSBtimes ( cfBuffer->multiSSB );
	  cfBuffer->multiSSB = multiSSB;
	  cfBuffer->Alpha = doppler->Alpha;
	  cfBuffer->Delta = doppler->Delta;
	  cfBuffer->multiDetStates = multiDetStates ;
	} /* buffer new SSB times */

    } /* could not reuse previously buffered quantites */

    /* new orbital parameter corrections if not already buffered */
  if ( doppler->orbit )
    {
      /* compute binary time corrections to the SSB time delays and SSB time derivitive */
      if ( XLALAddMultiBinaryTimes ( &multiBinary, multiSSB, doppler->orbit ) != XLAL_SUCCESS )
        {
          XLALPrintError("XLALAddMultiBinaryTimes() failed with xlalErrno = %d\n\n", xlalErrno );
          ABORTXLAL ( status );
        }
      multiSSBTotal = multiBinary;
    }
  else
    multiSSBTotal = multiSSB;

  /* special treatment of AM coefficients */
  if ( params->useRAA && !multiCmplxAMcoef )
    {
      /* compute new RAA AM-coefficients */
      LALGetMultiCmplxAMCoeffs ( status->statusPtr, &multiCmplxAMcoef, multiDetStates, *doppler );
      BEGINFAIL ( status ) {
	XLALDestroyMultiSSBtimes ( multiSSB );
      } ENDFAIL (status);

      /* noise-weight Antenna-patterns and compute A,B,C */
      if ( XLALWeightMultiCmplxAMCoeffs ( multiCmplxAMcoef, multiWeights ) != XLAL_SUCCESS ) {
	XLALPrintError("\nXLALWeightMultiCmplxAMCoeffs() failed with error = %d\n\n", xlalErrno );
	ABORT ( status, COMPUTEFSTATC_EXLAL, COMPUTEFSTATC_MSGEXLAL );
      }

      /* store in buffer if available */
      if ( cfBuffer )
	{
	  XLALDestroyMultiCmplxAMCoeffs ( cfBuffer->multiCmplxAMcoef );
	  cfBuffer->multiCmplxAMcoef = multiCmplxAMcoef;
	}

    } /* if RAA AM coefficients need to be computed */

  if ( !params->useRAA && !multiAMcoef )
    {
      /* compute new AM-coefficients */
      LALGetMultiAMCoeffs ( status->statusPtr, &multiAMcoef, multiDetStates, skypos );
      BEGINFAIL ( status ) {
	XLALDestroyMultiSSBtimes ( multiSSB );
      } ENDFAIL (status);

      /* noise-weight Antenna-patterns and compute A,B,C */
      if ( XLALWeightMultiAMCoeffs ( multiAMcoef, multiWeights ) != XLAL_SUCCESS ) {
	XLALPrintError("\nXLALWeightMultiAMCoeffs() failed with error = %d\n\n", xlalErrno );
	ABORT ( status, COMPUTEFSTATC_EXLAL, COMPUTEFSTATC_MSGEXLAL );
      }

      /* store these in buffer if available */
      if ( cfBuffer )
	{
	  XLALDestroyMultiAMCoeffs ( cfBuffer->multiAMcoef );
	  cfBuffer->multiAMcoef = multiAMcoef;
	} /* if cfBuffer */

    } /* if LWL AM coefficient need to be computed */

  if ( multiAMcoef )
    {
      Ad = multiAMcoef->Mmunu.Ad;
      Bd = multiAMcoef->Mmunu.Bd;
      Cd = multiAMcoef->Mmunu.Cd;
      Dd_inv = 1.0 / multiAMcoef->Mmunu.Dd;
      Ed = 0;
    }
  else if ( multiCmplxAMcoef )
    {
      Ad = multiCmplxAMcoef->Mmunu.Ad;
      Bd = multiCmplxAMcoef->Mmunu.Bd;
      Cd = multiCmplxAMcoef->Mmunu.Cd;
      Ed = multiCmplxAMcoef->Mmunu.Ed;
      Dd_inv = 1.0 / multiCmplxAMcoef->Mmunu.Dd;
    }
  else
    {
      XLALPrintError ( "Programming error: neither 'multiAMcoef' nor 'multiCmplxAMcoef' are available!\n");
      ABORT ( status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
    }

  /* if requested, prepare for returning single-IFO F-stat vector */
  if ( params->returnSingleF )
    {
      retF.numDetectors = numDetectors;
      if ( numDetectors > PULSAR_MAX_DETECTORS ) {
        XLALPrintError ("%s: numDetectors = %d exceeds currently allowed upper limit of detectors (%d) for returnSingleF=TRUE\n", __func__, numDetectors, PULSAR_MAX_DETECTORS );
        ABORT ( status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT );
      }
    }

  /* ----- loop over detectors and compute all detector-specific quantities ----- */
  for ( X=0; X < numDetectors; X ++)
    {
      Fcomponents FcX = empty_Fcomponents;	/* for detector-specific FaX, FbX */

      if ( params->useRAA )
	{
	  if ( XLALComputeFaFbCmplx (&FcX, multiSFTs->data[X], doppler->fkdot, multiSSBTotal->data[X], multiCmplxAMcoef->data[X], params) != 0)
	    {
	      XLALPrintError ("\nXALComputeFaFbCmplx() failed\n");
	      ABORT ( status, COMPUTEFSTATC_EXLAL, COMPUTEFSTATC_MSGEXLAL );
	    }
	}
      else if ( params->upsampling > 1)
	{
	  if ( XLALComputeFaFbXavie (&FcX, multiSFTs->data[X], doppler->fkdot, multiSSBTotal->data[X], multiAMcoef->data[X], params) != 0)
	    {
	      XLALPrintError ("\nXALComputeFaFbXavie() failed\n");
	      ABORT ( status, COMPUTEFSTATC_EXLAL, COMPUTEFSTATC_MSGEXLAL );
	    }
	}
      else
	{
	  if ( XLALComputeFaFb (&FcX, multiSFTs->data[X], doppler->fkdot, multiSSBTotal->data[X], multiAMcoef->data[X], params) != 0)
	    {
	      XLALPrintError ("\nXALComputeFaFb() failed\n");
	      ABORT ( status, COMPUTEFSTATC_EXLAL, COMPUTEFSTATC_MSGEXLAL );
	    }
	  if ( params->returnAtoms )
	    {
	      retF.multiFstatAtoms->data[X] = FcX.multiFstatAtoms->data[0];	/* copy pointer to IFO-specific Fstat-atoms 'contents' */
	      /* free 'container', but not *contents*, which have been linked above */
	      LALFree ( FcX.multiFstatAtoms->data );
	      LALFree ( FcX.multiFstatAtoms );
	    }
	}

#ifndef LAL_NDEBUG
      if ( !isfinite(creal(FcX.Fa)) || !isfinite(cimag(FcX.Fa)) || !isfinite(creal(FcX.Fb)) || !isfinite(cimag(FcX.Fb)) ) {
	XLALPrintError("XLALComputeFaFb() returned non-finite: Fa=(%f,%f), Fb=(%f,%f)\n",
		      creal(FcX.Fa), cimag(FcX.Fa), creal(FcX.Fb), cimag(FcX.Fb) );
	ABORT (status,  COMPUTEFSTATC_EIEEE,  COMPUTEFSTATC_MSGEIEEE);
      }
#endif

      /* compute single-IFO F-stats, if requested */
      if ( params->returnSingleF )
        {
         REAL8 AdX = multiAMcoef->data[X]->A;
         REAL8 BdX = multiAMcoef->data[X]->B;
         REAL8 CdX = multiAMcoef->data[X]->C;
         REAL8 DdX_inv = 1.0 / multiAMcoef->data[X]->D;

	 /* compute final single-IFO F-stat */
	 retF.FX[X] = DdX_inv * (  BdX * (SQ(creal(FcX.Fa)) + SQ(cimag(FcX.Fa)) )
	                           + AdX * ( SQ(creal(FcX.Fb)) + SQ(cimag(FcX.Fb)) )
		                   - 2.0 * CdX *( creal(FcX.Fa) * creal(FcX.Fb) + cimag(FcX.Fa) * cimag(FcX.Fb) )
		                   );
        } /* if returnSingleF */

      /* Fa = sum_X Fa_X */
      retF.Fa += FcX.Fa;

      /* Fb = sum_X Fb_X */
      retF.Fb += FcX.Fb;

    } /* for  X < numDetectors */

  /* ----- compute final Fstatistic-value ----- */

  /* NOTE: the data MUST be normalized by the DOUBLE-SIDED PSD (using LALNormalizeMultiSFTVect),
   * therefore there is a factor of 2 difference with respect to the equations in JKS, which
   * where based on the single-sided PSD.
   */
  retF.F = Dd_inv * (  Bd * (SQ(creal(retF.Fa)) + SQ(cimag(retF.Fa)) )
		       + Ad * ( SQ(creal(retF.Fb)) + SQ(cimag(retF.Fb)) )
		       - 2.0 * Cd *( creal(retF.Fa) * creal(retF.Fb) + cimag(retF.Fa) * cimag(retF.Fb) )
		       );

  if ( Ed != 0 ) /* extra term in RAA case */
    retF.F += - 2.0 * Dd_inv * Ed *( - creal(retF.Fa) * cimag(retF.Fb) + cimag(retF.Fa) * creal(retF.Fb) ); /* -2 E Im(Fa Fb^* ) / D */

  /* set correct F-stat reference time (taken from template 'doppler') [relevant only for phase of {Fa,Fb}] */
  retF.refTime = doppler->refTime;

  /* free memory if no buffer was available */
  if ( !cfBuffer )
    {
      XLALDestroyMultiSSBtimes ( multiSSB );
      XLALDestroyMultiAMCoeffs ( multiAMcoef );
      XLALDestroyMultiCmplxAMCoeffs ( multiCmplxAMcoef );
    } /* if !cfBuffer */

  /* this always needs to be free'ed, as it's no longer buffered */
  XLALDestroyMultiSSBtimes ( multiBinary );

  /* return final Fstat result */
  (*Fstat) = retF;

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* ComputeFStat() */


/**
 * Revamped version of LALDemod() (based on TestLALDemod() in CFS).
 * Compute JKS's Fa and Fb, which are ingredients for calculating the F-statistic.
 */
int
XLALComputeFaFb ( Fcomponents *FaFb,		      	/**< [out] Fa,Fb (and possibly atoms) returned */
		  const SFTVector *sfts,		/**< [in] input SFTs */
		  const PulsarSpins fkdot,		/**< [in] frequency and derivatives fkdot = d^kf/dt^k */
		  const SSBtimes *tSSB,			/**< [in] SSB timing series for particular sky-direction */
		  const AMCoeffs *amcoe,		/**< [in] antenna-pattern coefficients for this sky-direction */
		  const ComputeFParams *params )       	/**< addition computational params */
{
  UINT4 alpha;                 	/* loop index over SFTs */
  UINT4 spdnOrder;		/* maximal spindown-orders */
  UINT4 numSFTs;		/* number of SFTs (M in the Notes) */
  COMPLEX16 Fa, Fb;
  REAL8 Tsft; 			/* length of SFTs in seconds */
  INT4 freqIndex0;		/* index of first frequency-bin in SFTs */
  INT4 freqIndex1;		/* index of last frequency-bin in SFTs */

  REAL4 *a_al, *b_al;		/* pointer to alpha-arrays over a and b */
  REAL8 *DeltaT_al, *Tdot_al;	/* pointer to alpha-arrays of SSB-timings */
  SFTtype *SFT_al;		/* SFT alpha  */
  UINT4 Dterms = params->Dterms;

  REAL8 norm = OOTWOPI;

  /* ----- check validity of input */
#ifndef LAL_NDEBUG
  if ( !FaFb ) {
    XLALPrintError ("\nOutput-pointer is NULL !\n\n");
    XLAL_ERROR ( XLAL_EINVAL);
  }

  if ( !sfts || !sfts->data ) {
    XLALPrintError ("\nInput SFTs are NULL!\n\n");
    XLAL_ERROR ( XLAL_EINVAL);
  }

  if ( !tSSB || !tSSB->DeltaT || !tSSB->Tdot || !amcoe || !amcoe->a || !amcoe->b || !params)
    {
      XLALPrintError ("\nIllegal NULL in input !\n\n");
      XLAL_ERROR ( XLAL_EINVAL);
    }

  if ( PULSAR_MAX_SPINS > NUM_FACT )
    {
      XLALPrintError ("\nInverse factorials table only up to order s=%d, can't handle %d spin-order\n\n",
		     NUM_FACT, PULSAR_MAX_SPINS - 1 );
      XLAL_ERROR ( XLAL_EINVAL);
    }
#endif

  if ( params->upsampling > 1 ) {
    fprintf (stderr, "\n===== WARNING: XLALComputeFaFb() should not be used with upsampled-SFTs!\n");
    XLAL_ERROR ( XLAL_EINVAL);
  }

  /* ----- prepare convenience variables */
  numSFTs = sfts->length;
  Tsft = 1.0 / sfts->data[0].deltaF;
  {
    REAL8 dFreq = sfts->data[0].deltaF;
    freqIndex0 = (UINT4) ( sfts->data[0].f0 / dFreq + 0.5); /* lowest freqency-index */
    freqIndex1 = freqIndex0 + sfts->data[0].data->length;
  }

  /* ----- prepare return of 'FstatAtoms' if requested */
  if ( params->returnAtoms )
    {
      if ( (FaFb->multiFstatAtoms = LALMalloc ( sizeof(*FaFb->multiFstatAtoms) )) == NULL ){
	XLAL_ERROR ( XLAL_ENOMEM );
      }
      FaFb->multiFstatAtoms->length = 1;	/* in this function: single-detector only */
      if ( (FaFb->multiFstatAtoms->data = LALMalloc ( 1 * sizeof( *FaFb->multiFstatAtoms->data) )) == NULL ){
	LALFree (FaFb->multiFstatAtoms);
	XLAL_ERROR ( XLAL_ENOMEM );
      }
      if ( (FaFb->multiFstatAtoms->data[0] = XLALCreateFstatAtomVector ( numSFTs )) == NULL ) {
	LALFree ( FaFb->multiFstatAtoms->data );
	LALFree ( FaFb->multiFstatAtoms );
	XLAL_ERROR( XLAL_ENOMEM );
      }

      FaFb->multiFstatAtoms->data[0]->TAtom = Tsft;	/* time-baseline of returned atoms is Tsft */

    } /* if returnAtoms */

  /* ----- find highest non-zero spindown-entry */
  for ( spdnOrder = PULSAR_MAX_SPINS - 1;  spdnOrder > 0 ; spdnOrder --  )
    if ( fkdot[spdnOrder] )
      break;

  Fa = 0.0;
  Fb = 0.0;

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
      REAL4 s_alpha=0, c_alpha=0;/* sin(2pi kappa_alpha) and (cos(2pi kappa_alpha)-1) */
      REAL4 realQ, imagQ;	/* Re and Im of Q = e^{-i 2 pi lambda_alpha} */
      REAL4 realXP, imagXP;	/* Re/Im of sum_k X_ak * P_ak */
      REAL4 realQXP, imagQXP;	/* Re/Im of Q_alpha R_alpha */

      REAL8 lambda_alpha, kappa_max, kappa_star;
      COMPLEX8 Fa_alpha, Fb_alpha;

      /* ----- calculate kappa_max and lambda_alpha */
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
	  XLAL_ERROR ( XLAL_EFUNC);
	}

	kstar = (INT4) (Dphi_alpha);	/* k* = floor(Dphi_alpha) for positive Dphi */
	kappa_star = Dphi_alpha - 1.0 * kstar;	/* remainder of Dphi_alpha: >= 0 ! */
	kappa_max = kappa_star + 1.0 * Dterms - 1.0;

	/* ----- check that required frequency-bins are found in the SFTs ----- */
	k0 = kstar - Dterms + 1;
	k1 = k0 + 2 * Dterms - 1;
	if ( (k0 < freqIndex0) || (k1 > freqIndex1) )
	  {
	    XLALPrintError ("Required frequency-bins [%d, %d] not covered by SFT-interval [%d, %d]\n\n",
			   k0, k1, freqIndex0, freqIndex1 );
	    XLAL_ERROR(XLAL_EDOM);
	  }

      } /* compute kappa_star, lambda_alpha */

      /* NOTE: sin[ 2pi (Dphi_alpha - k) ] = sin [ 2pi Dphi_alpha ], therefore
       * the trig-functions need to be calculated only once!
       * We choose the value sin[ 2pi(Dphi_alpha - kstar) ] because it is the
       * closest to zero and will pose no numerical difficulties !
       */
      sin_cos_2PI_LUT ( &s_alpha, &c_alpha, kappa_star );
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
      if ( ( kappa_star > LD_SMALL4 ) && (kappa_star < 1.0 - LD_SMALL4) )
	{
	  /* improved hotloop algorithm by Fekete Akos:
	   * take out repeated divisions into a single common denominator,
	   * plus use extra cleverness to compute the nominator efficiently...
	   */
	  REAL4 Sn = crealf(*Xalpha_l);
	  REAL4 Tn = cimagf(*Xalpha_l);
	  REAL4 pn = kappa_max;
	  REAL4 qn = pn;
	  REAL4 U_alpha, V_alpha;

	  /* recursion with 2*Dterms steps */
	  UINT4 l;
	  for ( l = 1; l < 2*Dterms; l ++ )
	    {
	      Xalpha_l ++;

	      pn = pn - 1.0f; 			/* p_(n+1) */
	      Sn = pn * Sn + qn * crealf(*Xalpha_l);	/* S_(n+1) */
	      Tn = pn * Tn + qn * cimagf(*Xalpha_l);	/* T_(n+1) */
	      qn *= pn;				/* q_(n+1) */
	    } /* for l <= 2*Dterms */

	  U_alpha = Sn / qn;
	  V_alpha = Tn / qn;

#ifndef LAL_NDEBUG
          if ( !isfinite(U_alpha) || !isfinite(V_alpha) || !isfinite(pn) || !isfinite(qn) || !isfinite(Sn) || !isfinite(Tn) ) {
            XLALPrintError("XLALComputeFaFb() returned non-finite: U_alpha=%f, V_alpha=%f, pn=%f, qn=%f, Sn=%f, Tn=%f\n",
                           U_alpha, V_alpha, pn, qn, Sn, Tn);
            XLAL_ERROR (XLAL_EFPINVAL);
          }
#endif

	  realXP = s_alpha * U_alpha - c_alpha * V_alpha;
	  imagXP = c_alpha * U_alpha + s_alpha * V_alpha;

	} /* if |remainder| > LD_SMALL4 */
      else
	{ /* otherwise: lim_{rem->0}P_alpha,k  = 2pi delta_{k,kstar} */
	  UINT4 ind0;
  	  if ( kappa_star <= LD_SMALL4 ) ind0 = Dterms - 1;
  	  else ind0 = Dterms;
	  realXP = TWOPI_FLOAT * crealf(Xalpha_l[ind0]);
	  imagXP = TWOPI_FLOAT * cimagf(Xalpha_l[ind0]);
	} /* if |remainder| <= LD_SMALL4 */

      realQXP = realQ * realXP - imagQ * imagXP;
      imagQXP = realQ * imagXP + imagQ * realXP;

      /* we're done: ==> combine these into Fa and Fb */
      a_alpha = (*a_al);
      b_alpha = (*b_al);

      Fa_alpha = crectf( a_alpha * realQXP, a_alpha * imagQXP );
      Fa += crect( crealf(Fa_alpha), cimagf(Fa_alpha) );

      Fb_alpha = crectf( b_alpha * realQXP, b_alpha * imagQXP );
      Fb += crect( crealf(Fb_alpha), cimagf(Fb_alpha) );

      /* store per-SFT F-stat 'atoms' for transient-CW search */
      if ( params->returnAtoms )
	{
	  FaFb->multiFstatAtoms->data[0]->data[alpha].timestamp = (UINT4)XLALGPSGetREAL8( &SFT_al->epoch );
	  FaFb->multiFstatAtoms->data[0]->data[alpha].a2_alpha   = a_alpha * a_alpha;
	  FaFb->multiFstatAtoms->data[0]->data[alpha].b2_alpha   = b_alpha * b_alpha;
	  FaFb->multiFstatAtoms->data[0]->data[alpha].ab_alpha   = a_alpha * b_alpha;
	  FaFb->multiFstatAtoms->data[0]->data[alpha].Fa_alpha = (((REAL4) norm) * Fa_alpha);
	  FaFb->multiFstatAtoms->data[0]->data[alpha].Fb_alpha = (((REAL4) norm) * Fb_alpha);
	}

      /* advance pointers over alpha */
      a_al ++;
      b_al ++;
      DeltaT_al ++;
      Tdot_al ++;
      SFT_al ++;

    } /* for alpha < numSFTs */

  /* return result */
  FaFb->Fa = (((REAL8) norm) * Fa);
  FaFb->Fb = (((REAL8) norm) * Fb);

  return XLAL_SUCCESS;

} /* XLALComputeFaFb() */


/**
 * Revamped version of XLALComputeFaFb() for the case where a and b
 * are complex.
 * Compute JKS's Fa and Fb, which are ingredients for
 * calculating the F-statistic.
 */
int
XLALComputeFaFbCmplx ( Fcomponents *FaFb,		/**< [out] Fa,Fb (and possibly atoms) returned */
		  const SFTVector *sfts,               	/**< [in] input SFTs */
		  const PulsarSpins fkdot,             	/**< [in] frequency and derivatives fkdot = d^kf/dt^k */
		  const SSBtimes *tSSB,                	/**< [in] SSB timing series for particular sky-direction */
		  const CmplxAMCoeffs *amcoe,          	/**< [in] antenna-pattern coefficients for this sky-direction */
		  const ComputeFParams *params)      	/**< addition computational params */
{
  UINT4 alpha;                 	/* loop index over SFTs */
  UINT4 spdnOrder;		/* maximal spindown-orders */
  UINT4 numSFTs;		/* number of SFTs (M in the Notes) */
  COMPLEX16 Fa, Fb;
  REAL8 Tsft; 			/* length of SFTs in seconds */
  INT4 freqIndex0;		/* index of first frequency-bin in SFTs */
  INT4 freqIndex1;		/* index of last frequency-bin in SFTs */

  COMPLEX8 *a_al, *b_al;	/* pointer to alpha-arrays over a and b */
  REAL8 *DeltaT_al, *Tdot_al;	/* pointer to alpha-arrays of SSB-timings */
  SFTtype *SFT_al;		/* SFT alpha  */
  UINT4 Dterms = params->Dterms;

  REAL8 norm = OOTWOPI;

  /* ----- check validity of input */
#ifndef LAL_NDEBUG
  if ( !FaFb ) {
    XLALPrintError ("\nOutput-pointer is NULL !\n\n");
    XLAL_ERROR ( XLAL_EINVAL);
  }

  if ( !sfts || !sfts->data ) {
    XLALPrintError ("\nInput SFTs are NULL!\n\n");
    XLAL_ERROR ( XLAL_EINVAL);
  }

  if ( !tSSB || !tSSB->DeltaT || !tSSB->Tdot || !amcoe || !amcoe->a || !amcoe->b || !params)
    {
      XLALPrintError ("\nIllegal NULL in input !\n\n");
      XLAL_ERROR ( XLAL_EINVAL);
    }

  if ( PULSAR_MAX_SPINS > NUM_FACT )
    {
      XLALPrintError ("\nInverse factorials table only up to order s=%d, can't handle %d spin-order\n\n",
		     NUM_FACT, PULSAR_MAX_SPINS - 1 );
      XLAL_ERROR ( XLAL_EINVAL);
    }
  if ( params->upsampling > 1 ) {
    fprintf (stderr, "\n===== WARNING: XLALComputeFaFbCmplx() should not be used with upsampled-SFTs!\n");
    XLAL_ERROR ( XLAL_EINVAL);
  }

  if ( params->returnAtoms )
    {
      XLALPrintError ("%s: using the option 'returnAtoms' is not supported in this function!\n", __func__ );
      XLAL_ERROR ( XLAL_EINVAL);
    }
#endif

  /* ----- prepare convenience variables */
  numSFTs = sfts->length;
  Tsft = 1.0 / sfts->data[0].deltaF;
  {
    REAL8 dFreq = sfts->data[0].deltaF;
    freqIndex0 = (UINT4) ( sfts->data[0].f0 / dFreq + 0.5); /* lowest freqency-index */
    freqIndex1 = freqIndex0 + sfts->data[0].data->length;
  }

  /* find highest non-zero spindown-entry */
  for ( spdnOrder = PULSAR_MAX_SPINS - 1;  spdnOrder > 0 ; spdnOrder --  )
    if ( fkdot[spdnOrder] )
      break;

  Fa = 0.0;
  Fb = 0.0;

  a_al = amcoe->a->data;	/* point to beginning of alpha-arrays */
  b_al = amcoe->b->data;
  DeltaT_al = tSSB->DeltaT->data;
  Tdot_al = tSSB->Tdot->data;
  SFT_al = sfts->data;

  /* Loop over all SFTs  */
  for ( alpha = 0; alpha < numSFTs; alpha++ )
    {
      COMPLEX8 a_alpha, b_alpha;

      INT4 kstar;		/* central frequency-bin k* = round(xhat_alpha) */
      INT4 k0, k1;

      COMPLEX8 *Xalpha = SFT_al->data->data; /* pointer to current SFT-data */
      COMPLEX8 *Xalpha_l; 	/* pointer to frequency-bin k in current SFT */
      REAL4 s_alpha=0, c_alpha=0;/* sin(2pi kappa_alpha) and (cos(2pi kappa_alpha)-1) */
      REAL4 realQ, imagQ;	/* Re and Im of Q = e^{-i 2 pi lambda_alpha} */
      REAL4 realXP, imagXP;	/* Re/Im of sum_k X_ak * P_ak */
      REAL4 realQXP, imagQXP;	/* Re/Im of Q_alpha R_alpha */

      REAL8 lambda_alpha, kappa_max, kappa_star;

      /* ----- calculate kappa_max and lambda_alpha */
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
	  XLAL_ERROR ( XLAL_EFUNC);
	}

	kstar = (INT4) (Dphi_alpha);	/* k* = floor(Dphi_alpha) for positive Dphi */
	kappa_star = Dphi_alpha - 1.0 * kstar;	/* remainder of Dphi_alpha: >= 0 ! */
	kappa_max = kappa_star + 1.0 * Dterms - 1.0;

	/* ----- check that required frequency-bins are found in the SFTs ----- */
	k0 = kstar - Dterms + 1;
	k1 = k0 + 2 * Dterms - 1;
	if ( (k0 < freqIndex0) || (k1 > freqIndex1) )
	  {
	    XLALPrintError ("Required frequency-bins [%d, %d] not covered by SFT-interval [%d, %d]\n\n",
			   k0, k1, freqIndex0, freqIndex1 );
	    XLAL_ERROR(XLAL_EDOM);
	  }

      } /* compute kappa_star, lambda_alpha */

      /* NOTE: sin[ 2pi (Dphi_alpha - k) ] = sin [ 2pi Dphi_alpha ], therefore
       * the trig-functions need to be calculated only once!
       * We choose the value sin[ 2pi(Dphi_alpha - kstar) ] because it is the
       * closest to zero and will pose no numerical difficulties !
       */
      sin_cos_2PI_LUT ( &s_alpha, &c_alpha, kappa_star );
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
      if ( ( kappa_star > LD_SMALL4 ) && (kappa_star < 1.0 - LD_SMALL4) )
	{
	  /* improved hotloop algorithm by Fekete Akos:
	   * take out repeated divisions into a single common denominator,
	   * plus use extra cleverness to compute the nominator efficiently...
	   */
	  REAL4 Sn = crealf(*Xalpha_l);
	  REAL4 Tn = cimagf(*Xalpha_l);
	  REAL4 pn = kappa_max;
	  REAL4 qn = pn;
	  REAL4 U_alpha, V_alpha;

	  /* recursion with 2*Dterms steps */
	  UINT4 l;
	  for ( l = 1; l < 2*Dterms; l ++ )
	    {
	      Xalpha_l ++;

	      pn = pn - 1.0f; 			/* p_(n+1) */
	      Sn = pn * Sn + qn * crealf(*Xalpha_l);	/* S_(n+1) */
	      Tn = pn * Tn + qn * cimagf(*Xalpha_l);	/* T_(n+1) */
	      qn *= pn;				/* q_(n+1) */
	    } /* for l <= 2*Dterms */

	  U_alpha = Sn / qn;
	  V_alpha = Tn / qn;

#ifndef LAL_NDEBUG
          if ( !isfinite(U_alpha) || !isfinite(V_alpha) || !isfinite(pn) || !isfinite(qn) || !isfinite(Sn) || !isfinite(Tn) ) {
            XLALPrintError("XLALComputeFaFbCmplx() returned non-finite: U_alpha=%f, V_alpha=%f, pn=%f, qn=%f, Sn=%f, Tn=%f\n",
                           U_alpha, V_alpha, pn, qn, Sn, Tn);
            XLAL_ERROR (XLAL_EFPINVAL);
          }
#endif

	  realXP = s_alpha * U_alpha - c_alpha * V_alpha;
	  imagXP = c_alpha * U_alpha + s_alpha * V_alpha;

	} /* if |remainder| > LD_SMALL4 */
      else
	{ /* otherwise: lim_{rem->0}P_alpha,k  = 2pi delta_{k,kstar} */
	  UINT4 ind0;
  	  if ( kappa_star <= LD_SMALL4 ) ind0 = Dterms - 1;
  	  else ind0 = Dterms;
	  realXP = TWOPI_FLOAT * crealf(Xalpha_l[ind0]);
	  imagXP = TWOPI_FLOAT * cimagf(Xalpha_l[ind0]);
	} /* if |remainder| <= LD_SMALL4 */

      realQXP = realQ * realXP - imagQ * imagXP;
      imagQXP = realQ * imagXP + imagQ * realXP;

      /* we're done: ==> combine these into Fa and Fb */
      a_alpha = (*a_al);
      b_alpha = (*b_al);

      /* Fa contains complex conjugate of a */
      Fa += crect( crealf(a_alpha) * realQXP + cimagf(a_alpha) * imagQXP, crealf(a_alpha) * imagQXP - cimagf(a_alpha) * realQXP );

      /* Fb contains complex conjugate of b */
      Fb += crect( crealf(b_alpha) * realQXP + cimagf(b_alpha) * imagQXP, crealf(b_alpha) * imagQXP - cimagf(b_alpha) * realQXP );

      /* advance pointers over alpha */
      a_al ++;
      b_al ++;
      DeltaT_al ++;
      Tdot_al ++;
      SFT_al ++;

    } /* for alpha < numSFTs */

  /* return result */
  FaFb->Fa = (((REAL8) norm) * Fa);
  FaFb->Fb = (((REAL8) norm) * Fb);

  return XLAL_SUCCESS;

} /* XLALComputeFaFbCmplx() */


/**
 * Modified version of ComputeFaFb() based on Xavies trick:
 * need sufficiently oversampled SFTs and uses ZERO Dterms.
 * Compute JKS's Fa and Fb, which are ingredients for calculating the F-statistic.
 */
int
XLALComputeFaFbXavie ( Fcomponents *FaFb,		/**< [out] Fa,Fb (and possibly atoms) returned */
		       const SFTVector *sfts,          	/**< [in] input SFTs */
		       const PulsarSpins fkdot,        	/**< [in] frequency and derivatives fkdot = d^kf/dt^k */
		       const SSBtimes *tSSB,           	/**< [in] SSB timing series for particular sky-direction */
		       const AMCoeffs *amcoe,          	/**< [in] antenna-pattern coefficients for this sky-direction */
		       const ComputeFParams *params    	/**< additional computational params */
                       )
{
  UINT4 alpha;                 	/* loop index over SFTs */
  UINT4 spdnOrder;		/* maximal spindown-orders */
  UINT4 numSFTs;		/* number of SFTs (M in the Notes) */
  COMPLEX16 Fa, Fb;
  REAL8 Tsft; 			/* length of SFTs in seconds */
  INT4 freqIndex0;		/* index of first frequency-bin in SFTs */
  INT4 freqIndex1;		/* index of last frequency-bin in SFTs */

  REAL4 *a_al, *b_al;		/* pointer to alpha-arrays over a and b */
  REAL8 *DeltaT_al, *Tdot_al;	/* pointer to alpha-arrays of SSB-timings */
  SFTtype *SFT_al;		/* SFT alpha  */

  REAL4 Upsampling;

  /* ----- check validity of input */
#ifndef LAL_NDEBUG
  if ( !FaFb ) {
    XLALPrintError ("\nOutput-pointer is NULL !\n\n");
    XLAL_ERROR ( XLAL_EINVAL);
  }

  if ( !sfts || !sfts->data ) {
    XLALPrintError ("\nInput SFTs are NULL!\n\n");
    XLAL_ERROR ( XLAL_EINVAL);
  }

  if ( !tSSB || !tSSB->DeltaT || !tSSB->Tdot || !amcoe || !amcoe->a || !amcoe->b || !params)
    {
      XLALPrintError ("\nIllegal NULL in input !\n\n");
      XLAL_ERROR ( XLAL_EINVAL);
    }

  if ( PULSAR_MAX_SPINS > NUM_FACT )
    {
      XLALPrintError ("\nInverse factorials table only up to order s=%d, can't handle %d spin-order\n\n",
		     NUM_FACT, PULSAR_MAX_SPINS - 1 );
      XLAL_ERROR ( XLAL_EINVAL);
    }
  if ( params->returnAtoms )
    {
      XLALPrintError ("%s: using the option 'returnAtoms' is not supported in this function!\n", __func__ );
      XLAL_ERROR ( XLAL_EINVAL);
    }
#endif

  /* ----- prepare convenience variables */
  Upsampling = (REAL4) params->upsampling;

  numSFTs = sfts->length;
  Tsft = 1.0 / sfts->data[0].deltaF;
  {
    REAL8 dFreq = sfts->data[0].deltaF;
    freqIndex0 = (UINT4) ( sfts->data[0].f0 / dFreq + 0.5); /* lowest freqency-index */
    freqIndex0 *= Upsampling;
    freqIndex1 = freqIndex0 + sfts->data[0].data->length;
  }

  /* find highest non-zero spindown-entry */
  for ( spdnOrder = PULSAR_MAX_SPINS - 1;  spdnOrder > 0 ; spdnOrder --  )
    if ( fkdot[spdnOrder] )
      break;

  Fa = 0.0;
  Fb = 0.0;

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

      COMPLEX8 *Xalpha = SFT_al->data->data; /* pointer to current SFT-data */
      COMPLEX8 Xalpha_l; 	/* frequency-bin k in current SFT */
      REAL4 realQ, imagQ;	/* Re and Im of Q = e^{-i 2 pi lambda_alpha} */
      REAL4 realQXP, imagQXP;	/* Re/Im of Q_alpha R_alpha */

      REAL8 lambda_alpha;	/* !NOTE!: this MUST be REAL8!!! otherwise you lose the signal! */

      /* ----- calculate kappa_max and lambda_alpha */
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
	  XLAL_ERROR (XLAL_EFUNC);
	}

	kstar = (INT4) (Dphi_alpha * Upsampling + 0.5f - freqIndex0);	/* k* = round(Dphi_alpha*chi) for positive Dphi */

	/* ----- check that required frequency-bins are found in the SFTs ----- */
	if ( (kstar < 0) || (kstar > freqIndex1 - freqIndex0) )
	  {
	    XLALPrintError ("Required frequency-bin [%d] not covered by SFT-interval [%d, %d]\n\n",
			   freqIndex0 + kstar, freqIndex0, freqIndex1 );
	    XLAL_ERROR(XLAL_EDOM);
	  }

      } /* compute kstar, lambda_alpha */


      /* ---------- calculate the (truncated to ZERO Dterms) sum over k ---------- */

      /* ---------- ATTENTION: this the "hot-loop", which will be
       * executed many millions of times, so anything in here
       * has a HUGE impact on the whole performance of the code.
       *
       * DON'T touch *anything* in here unless you really know
       * what you're doing !!
       *------------------------------------------------------------
       */

      Xalpha_l = Xalpha[kstar];  /* frequency-bin to use */

      /* lim_{kappa_star->0}P_alpha,k  = 2pi delta_{k,kstar} */

      /* combine with e^-i 2pi lambda_alpha */
      realQXP = realQ * crealf(Xalpha_l) - imagQ * cimagf(Xalpha_l);
      imagQXP = realQ * cimagf(Xalpha_l) + imagQ * crealf(Xalpha_l);

      /* we're done: ==> combine these into Fa and Fb */
      a_alpha = (*a_al);
      b_alpha = (*b_al);

      Fa += crect( a_alpha * realQXP, a_alpha * imagQXP );
      Fb += crect( b_alpha * realQXP, b_alpha * imagQXP );

      /* advance pointers over alpha */
      a_al ++;
      b_al ++;
      DeltaT_al ++;
      Tdot_al ++;
      SFT_al ++;

    } /* for alpha < numSFTs */

  /* return result */
  FaFb->Fa = Fa;
  FaFb->Fb = Fb;

  return XLAL_SUCCESS;

} /* XLALComputeFaFbXavie() */

/**
 * For a given OrbitalParams, return a newly-allocated SSBtimes structure
 * containing the input SSBtimes with the *additional* time-differences due to
 * the binary orbital motion *added* to it.
 *
 * NOTE: the output vector 'tSSBOut' can be passed either
 * - unallocated (where it must be (*tSSBOut)==NULL), and it gets allocated here, or
 * - it can also contain a pre-existing vector, which must then be consistent (same vector lenghts)
 * than the input SSB-vectors.
 * This is intended to minimize unnecessary repeated memory allocs+frees on successive binary templates.
 *
 * NOTE2: it is allowed to pass the same in- and output-vectors, i.e. (*tSSBOut) = tSSBIn,
 * in which case the input vector will be safely added to.
 */
int
XLALAddBinaryTimes ( SSBtimes **tSSBOut,			/**< [out] SSB timings tSSBIn with binary offsets added, can be NULL */
                     const SSBtimes *tSSBIn,			/**< [in] SSB timings DeltaT_alpha = T(t_alpha) - T_0; and Tdot(t_alpha) */
                     const BinaryOrbitParams *binaryparams	/**< [in] source binary orbit parameters */
                     )
{
  XLAL_CHECK ( tSSBIn != NULL, XLAL_EINVAL, "Invalid NULL input 'tSSB'\n" );
  XLAL_CHECK ( tSSBIn->DeltaT != NULL, XLAL_EINVAL, "Invalid NULL input 'tSSBIn->DeltaT'\n" );
  XLAL_CHECK ( tSSBIn->Tdot != NULL, XLAL_EINVAL, "Invalid NULL input 'tSSBIn->Tdot'\n" );
  XLAL_CHECK ( binaryparams != NULL, XLAL_EINVAL, "Invalid NULL input 'binaryparams'\n");

  UINT4 numSteps = tSSBIn->DeltaT->length;		/* number of timesteps */
  XLAL_CHECK (tSSBIn->Tdot->length == numSteps, XLAL_EINVAL,
                   "Length tSSBIn->DeltaT = %d, while tSSBIn->Tdot = %d\n", numSteps, tSSBIn->Tdot->length );

  SSBtimes *binaryTimes;
  // ----- prepare output timeseries: either allocate or re-use existing
  if ( (*tSSBOut) == NULL )	// creating new output vector
    {
      XLAL_CHECK ( (binaryTimes = XLALDuplicateSSBtimes ( tSSBIn )) != NULL, XLAL_EFUNC );
    }
  else if ( (*tSSBOut) == tSSBIn )	// input==output vector
    {
      binaryTimes = (*tSSBOut);
    }
  else // input vector given, but not identical to output vector
    {
      binaryTimes = (*tSSBOut);
      // need to do a few more sanity checks
      XLAL_CHECK ( binaryTimes->DeltaT->length == numSteps, XLAL_EINVAL,
                   "Length (*tSSBOut)->DeltaT = %d, while tSSBIn->DeltaT = %d\n", binaryTimes->DeltaT->length, numSteps );
      XLAL_CHECK ( binaryTimes->Tdot->length == numSteps, XLAL_EINVAL,
                        "Length tSSBOut->Tdot = %d, while tSSBIn->Tdot = %d\n", binaryTimes->Tdot->length, numSteps );
      // ... and copy the vector contents from the input SSB vector
      binaryTimes->refTime = tSSBIn->refTime;
      memcpy ( binaryTimes->DeltaT->data, tSSBIn->DeltaT->data, numSteps * sizeof(binaryTimes->DeltaT->data[0]) );
      memcpy ( binaryTimes->Tdot->data,   tSSBIn->Tdot->data,   numSteps * sizeof(binaryTimes->Tdot->data[0]) );
    } // re-using input vector

  /* ----- convenience variables */
  REAL8 Porb = binaryparams->period;		/* binary orbital period */
  REAL8 e = binaryparams->ecc;			/* the eccentricity */
  REAL8 ome = 1.0 - e;				/* 1 minus eccentricity */
  REAL8 asini = binaryparams->asini;		/* the projected orbital semimajor axis */
  REAL8 sinw = sin(binaryparams->argp);		/* the sin and cos of the argument of periapsis */
  REAL8 cosw = cos(binaryparams->argp);

  REAL8 refTimeREAL8 = XLALGPSGetREAL8 ( &tSSBIn->refTime );

  /* compute (global) p, q and r coeeficients for root-finder */
  REAL8 p = (LAL_TWOPI/Porb)*cosw*asini*sqrt(1.0-e*e);
  REAL8 q = (LAL_TWOPI/Porb)*sinw*asini;
  REAL8 r = (LAL_TWOPI/Porb)*sinw*asini*ome;

  /* Calculate the required accuracy for the root finding procedure in the main loop */
  REAL8 acc = LAL_TWOPI*(REAL8)EA_ACC/Porb;   /* EA_ACC is defined above and represents the required timing precision in seconds (roughly) */

  /* loop over the SFTs */
  for ( UINT4 i = 0; i < numSteps; i++ )
    {
      /* define SSB time for the current SFT midpoint */
      REAL8 tSSB_now = refTimeREAL8 + tSSBIn->DeltaT->data[i];

      /* define fractional orbit in SSB frame since periapsis (enforce result 0->1) */
      /* the result of fmod uses the dividend sign hence the second procedure */
      REAL8 temp = fmod((tSSB_now - XLALGPSGetREAL8 ( &binaryparams->tp ) ),Porb)/(REAL8)Porb;
      REAL8 fracorb = temp - floor(temp);	        /* the fraction of orbits completed since current SSB time */

      // ---------- FIXME: can we use GSL for this root-finding?
      REAL8 E;              /* the eccentric anomoly */
      {
        const gsl_root_fsolver_type *T = gsl_root_fsolver_bisection;
        gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
        REAL8 E_lo = 0, E_hi = LAL_TWOPI;
        gsl_function F;
        struct E_solver_params pars = {p, q, r, fracorb};
        F.function = &gsl_E_solver;
        F.params = &pars;

        XLAL_CHECK ( gsl_root_fsolver_set(s, &F, E_lo, E_hi) == 0, XLAL_EFAILED );

        XLALPrintInfo ("%5s [%9s, %9s] %9s %10s %9s\n", "iter", "lower", "upper", "root", "abstol", "err(est)");
        int max_iter = 100;
        int iter = 0;
        int status;
        do
          {
            iter++;
            status = gsl_root_fsolver_iterate(s);
            XLAL_CHECK ( (status == GSL_SUCCESS) || (status == GSL_CONTINUE), XLAL_EFAILED );
            E = gsl_root_fsolver_root(s);
            E_lo = gsl_root_fsolver_x_lower (s);
            E_hi = gsl_root_fsolver_x_upper (s);
            status = gsl_root_test_interval ( E_lo, E_hi, acc, 0 );

            if (status == GSL_SUCCESS) { XLALPrintInfo ("Converged:\n"); }
            XLALPrintInfo ("%5d [%.7f, %.7f] %.7f %+10.7g %10.7g\n", iter, E_lo, E_hi, E, acc, E_hi - E_lo);

          } while ( (status == GSL_CONTINUE) && (iter < max_iter) );

        XLAL_CHECK ( status == GSL_SUCCESS, XLAL_EMAXITER, "Eccentric anomaly: failed to converge within %d iterations\n", max_iter );
        gsl_root_fsolver_free(s);
      }

      /* use our value of E to compute the additional binary time delay */
      binaryTimes->DeltaT->data[i] += - ( asini*sinw*(cos(E)-e) + asini*cosw*sqrt(1.0-e*e)*sin(E) );

      /* combine with Tdot (dtSSB_by_dtdet) -> dtbin_by_dtdet */
      binaryTimes->Tdot->data[i] *= ( (1.0 - e*cos(E))/(1.0 + p*cos(E) - q*sin(E)) );

    } /* for i < numSteps */

  // pass back output SSB timings
  (*tSSBOut) = binaryTimes;

  return XLAL_SUCCESS;

} /* XLALAddBinaryTimes() */

/**
 * Multi-IFO version of XLALAddBinaryTimes().
 * For a given OrbitalParams, return a newly-allocated multiSSBtimes structure
 * containing the input SSBtimes with the *additional* time-differences due to
 * the binary orbital motion *added* to it.
 *
 * NOTE: the output vector 'multiSSBOut' can be passed either
 * - unallocated (where it must be (*multiSSBOut)==NULL), and it gets allocated here, or
 * - it can also contain a pre-existing vector, which must then be consistent (same vector lenghts)
 * than the input SSB-vectors.
 * This is intended to minimize unnecessary repeated memory allocs+frees on successive binary templates.
 *
 * NOTE2: it is allowed to pass the same in- and output-vectors, i.e. (*multiSSBOut) = multiSSBIn,
 * in which case the input vector will be safely added to.
 *
 */
int
XLALAddMultiBinaryTimes ( MultiSSBtimes **multiSSBOut,		/**< [out] output SSB times */
                          const MultiSSBtimes *multiSSBIn,	/**< [in] SSB-timings for all input detector-state series */
                          const BinaryOrbitParams *binaryparams	/**< [in] source binary orbit parameters, NULL = isolated system */
                          )
{
  /* check input */
  XLAL_CHECK ( multiSSBIn != NULL, XLAL_EINVAL, "Invalid NULL input 'multiSSB'\n");
  XLAL_CHECK ( binaryparams != NULL, XLAL_EINVAL, "Invalid NULL input 'binaryparams'\n");

  UINT4 numDetectors = multiSSBIn->length;

  MultiSSBtimes *multiBinaryTimes;
  // ----- prepare output timeseries: either allocate or re-use existing
  if ( (*multiSSBOut) == NULL )	// creating new output vector
    {
      XLAL_CHECK ( (multiBinaryTimes = XLALDuplicateMultiSSBtimes ( multiSSBIn )) != NULL, XLAL_EFUNC );
    }
  else // input vector given
    {
      multiBinaryTimes = (*multiSSBOut);
      XLAL_CHECK ( multiBinaryTimes->length == numDetectors, XLAL_EINVAL,
                   "Inconsistent length (*multiSSBOut)->length = %d, while multiSSBIn->length = %d\n", (*multiSSBOut)->length, numDetectors );
      // we'll leave all other sanity-checks to XLALAddBinaryTimes() calls
    }

  // ----- simply loop over detectors for XLALAddBinaryTimes()
  for ( UINT4 X = 0; X < numDetectors; X ++ )
    {
      int ret = XLALAddBinaryTimes ( &(multiBinaryTimes->data[X]), multiSSBIn->data[X], binaryparams );
      XLAL_CHECK( ret == XLAL_SUCCESS, XLAL_EFUNC, "XLALAddBinaryTimes() failed for X=%d\n", X );
    } /* for X < numDet */

  // pass back result-vector
  (*multiSSBOut) = multiBinaryTimes;

  return XLAL_SUCCESS;

} /* XLALAddMultiBinaryTimes() */


/**
 * Duplicate (ie allocate + copy) an input SSBtimes structure.
 * This can be useful for creating a copy before adding binary-orbital corrections in XLALAddBinaryTimes()
 */
SSBtimes *
XLALDuplicateSSBtimes ( const SSBtimes *tSSB )
{
  XLAL_CHECK_NULL ( tSSB != NULL, XLAL_EINVAL, "Invalid NULL input 'tSSB'\n" );

  UINT4 len;
  SSBtimes *ret;
  ret = XLALCalloc ( 1, len = sizeof (*ret) );
  XLAL_CHECK_NULL ( ret != NULL, XLAL_ENOMEM, "Failed to XLALCalloc ( 1, %d )\n", len );

  ret->refTime = tSSB->refTime;

  if ( tSSB->DeltaT )
    {
      len = tSSB->DeltaT->length;
      ret->DeltaT = XLALCreateREAL8Vector ( len );
      XLAL_CHECK_NULL ( ret->DeltaT != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector(%d) failed\n", len );
      memcpy ( ret->DeltaT->data, tSSB->DeltaT->data, len * sizeof(ret->DeltaT->data[0]) );
    }

  if ( tSSB->Tdot )
    {
      len = tSSB->Tdot->length;
      ret->Tdot = XLALCreateREAL8Vector ( len );
      XLAL_CHECK_NULL ( ret->Tdot != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector(%d) failed\n", len );
      memcpy ( ret->Tdot->data, tSSB->Tdot->data, len * sizeof(ret->Tdot->data[0]) );
    }

  return ret;

} /* XLALDuplicateSSBtimes() */


/**
 * Duplicate (ie allocate + copy) an input MultiSSBtimes structure.
 */
MultiSSBtimes *
XLALDuplicateMultiSSBtimes ( const MultiSSBtimes *multiSSB )
{
  XLAL_CHECK_NULL ( multiSSB != NULL, XLAL_EINVAL, "Invalid NULL input 'multiSSB'\n");

  UINT4 len;
  MultiSSBtimes *ret;
  ret = XLALCalloc ( 1, len=sizeof(*ret) );
  XLAL_CHECK_NULL ( ret != NULL, XLAL_ENOMEM, "Failed to XLALCalloc ( 1, %d )\n", len );

  UINT4 numDetectors = multiSSB->length;
  ret->length = numDetectors;
  ret->data = XLALCalloc ( numDetectors, len = sizeof(ret->data[0]) );
  XLAL_CHECK_NULL ( ret->data != NULL, XLAL_ENOMEM, "Failed to XLALCalloc ( %d, %d )\n", numDetectors, len );

  for ( UINT4 X = 0; X < numDetectors; X ++ )
    {
      ret->data[X] = XLALDuplicateSSBtimes ( multiSSB->data[X] );
      XLAL_CHECK_NULL( ret->data[X] != NULL, XLAL_EFUNC, "XLALDuplicateSSBtimes() failed for detector X=%d\n", X );
    } // for X < numDetectors

  return ret;

} /* XLALDuplicateMultiSSBtimes() */

/**
 * For a given DetectorStateSeries, calculate the time-differences
 * \f$\Delta T_\alpha\equiv T(t_\alpha) - T_0\f$, and their
 * derivatives \f$\dot{T}_\alpha \equiv d T / d t (t_\alpha)\f$.
 *
 * \note The return-vector is allocated here
 *
 */
SSBtimes *
XLALGetSSBtimes ( const DetectorStateSeries *DetectorStates,	/**< [in] detector-states at timestamps t_i */
                  SkyPosition pos,				/**< source sky-location */
                  LIGOTimeGPS refTime,				/**< SSB reference-time T_0 of pulsar-parameters */
                  SSBprecision precision			/**< relativistic or Newtonian SSB transformation? */
                  )
{
  XLAL_CHECK_NULL ( DetectorStates != NULL, XLAL_EINVAL, "Invalid NULL input 'DetectorStates'\n" );
  XLAL_CHECK_NULL ( precision < SSBPREC_LAST, XLAL_EDOM, "Invalid value precision=%d, allowed are [0, %d]\n", precision, SSBPREC_LAST -1 );
  XLAL_CHECK_NULL ( pos.system == COORDINATESYSTEM_EQUATORIAL, XLAL_EDOM, "Only equatorial coordinate system (=%d) allowed, got %d\n", COORDINATESYSTEM_EQUATORIAL, pos.system );

  UINT4 numSteps = DetectorStates->length;		/* number of timestamps */

  // prepare output SSBtimes struct
  int len;
  SSBtimes *ret = XLALCalloc ( 1, len = sizeof(*ret) );
  XLAL_CHECK_NULL ( ret != NULL, XLAL_ENOMEM, "Failed to XLALCalloc(1,%d)\n", len );
  ret->DeltaT = XLALCreateREAL8Vector ( numSteps );
  XLAL_CHECK_NULL ( ret->DeltaT != NULL, XLAL_EFUNC, "ret->DeltaT = XLALCreateREAL8Vector(%d) failed\n", numSteps );
  ret->Tdot = XLALCreateREAL8Vector ( numSteps );
  XLAL_CHECK_NULL ( ret->Tdot != NULL, XLAL_EFUNC, "ret->Tdot = XLALCreateREAL8Vector(%d) failed\n", numSteps );

  /* convenience variables */
  REAL8 alpha = pos.longitude;
  REAL8 delta = pos.latitude;
  REAL8 refTimeREAL8 = XLALGPSGetREAL8 ( &refTime );

  BarycenterInput baryinput = empty_BarycenterInput;

  /*----- now calculate the SSB transformation in the precision required */
  switch (precision)
    {
      REAL8 vn[3];		/* unit-vector pointing to source in Cart. coord. */

    case SSBPREC_NEWTONIAN:	/* use simple vr.vn to calculate time-delay */

      /*----- get the cartesian source unit-vector */
      vn[0] = cos(alpha) * cos(delta);
      vn[1] = sin(alpha) * cos(delta);
      vn[2] = sin(delta);

      for (UINT4 i = 0; i < numSteps; i++ )
	{
	  LIGOTimeGPS *ti = &(DetectorStates->data[i].tGPS);
	  /* DeltaT_alpha */
	  ret->DeltaT->data[i]  = XLALGPSGetREAL8 ( ti );
	  ret->DeltaT->data[i] += SCALAR(vn, DetectorStates->data[i].rDetector);
	  ret->DeltaT->data[i] -= refTimeREAL8;

	  /* Tdot_alpha */
	  ret->Tdot->data[i] = 1.0 + SCALAR(vn, DetectorStates->data[i].vDetector);

	} /* for i < numSteps */

      break;

    case SSBPREC_RELATIVISTIC:	/* use LALBarycenter() to get SSB-times and derivative */

      baryinput.site = DetectorStates->detector;
      baryinput.site.location[0] /= LAL_C_SI;
      baryinput.site.location[1] /= LAL_C_SI;
      baryinput.site.location[2] /= LAL_C_SI;

      baryinput.alpha = alpha;
      baryinput.delta = delta;
      baryinput.dInv = 0;

      for (UINT4 i = 0; i < numSteps; i++ )
	{
	  EmissionTime emit;
	  DetectorState *state = &(DetectorStates->data[i]);

	  baryinput.tgps = state->tGPS;

          if ( XLALBarycenter ( &emit, &baryinput, &(state->earthState) ) != XLAL_SUCCESS )
            XLAL_ERROR_NULL ( XLAL_EFUNC, "XLALBarycenter() failed with xlalErrno = %d\n", xlalErrno );

	  ret->DeltaT->data[i] = XLALGPSGetREAL8 ( &emit.te ) - refTimeREAL8;
	  ret->Tdot->data[i] = emit.tDot;

	} /* for i < numSteps */

      break;

    case SSBPREC_RELATIVISTICOPT:	/* use optimized version XLALBarycenterOpt() */

      baryinput.site = DetectorStates->detector;
      baryinput.site.location[0] /= LAL_C_SI;
      baryinput.site.location[1] /= LAL_C_SI;
      baryinput.site.location[2] /= LAL_C_SI;

      baryinput.alpha = alpha;
      baryinput.delta = delta;
      baryinput.dInv = 0;

      BarycenterBuffer *bBuffer = NULL;
      for ( UINT4 i = 0; i < numSteps; i++ )
        {
          EmissionTime emit;
          DetectorState *state = &(DetectorStates->data[i]);
          baryinput.tgps = state->tGPS;

          if ( XLALBarycenterOpt ( &emit, &baryinput, &(state->earthState), &bBuffer ) != XLAL_SUCCESS )
            XLAL_ERROR_NULL ( XLAL_EFUNC, "XLALBarycenterOpt() failed with xlalErrno = %d\n", xlalErrno );

          ret->DeltaT->data[i] = XLALGPSGetREAL8 ( &emit.te ) - refTimeREAL8;
          ret->Tdot->data[i] = emit.tDot;

        } /* for i < numSteps */
      // free buffer memory
      XLALFree ( bBuffer );

      break;

    default:
      XLAL_ERROR_NULL (XLAL_EFAILED, "\n?? Something went wrong.. this should never be called!\n\n" );
      break;
    } /* switch precision */

  /* finally: store the reference-time used into the output-structure */
  ret->refTime = refTime;

  return ret;

} /* XLALGetSSBtimes() */

/**
 * Multi-IFO version of LALGetSSBtimes().
 * Get all SSB-timings for all input detector-series.
 *
 * NOTE: this functions *allocates* the output-vector,
 * use XLALDestroyMultiSSBtimes() to free this.
 */
MultiSSBtimes *
XLALGetMultiSSBtimes ( const MultiDetectorStateSeries *multiDetStates, /**< [in] detector-states at timestamps t_i */
                       SkyPosition skypos,		/**< source sky-position [in equatorial coords!] */
                       LIGOTimeGPS refTime,		/**< SSB reference-time T_0 for SSB-timing */
                       SSBprecision precision		/**< use relativistic or Newtonian SSB timing?  */
                       )
{
  /* check input */
  XLAL_CHECK_NULL ( multiDetStates != NULL, XLAL_EINVAL, "Invalid NULL input 'multiDetStates'\n");
  XLAL_CHECK_NULL ( multiDetStates->length > 0, XLAL_EINVAL, "Invalid zero-length 'multiDetStates'\n");

  UINT4 numDetectors = multiDetStates->length;

  // prepare return struct
  int len;
  MultiSSBtimes *ret = XLALCalloc ( 1, len = sizeof( *ret ) );
  XLAL_CHECK_NULL ( ret != NULL, XLAL_ENOMEM, "Failed to XLALCalloc(1,%d)\n", len );
  ret->length = numDetectors;
  ret->data = XLALCalloc ( numDetectors, len=sizeof ( *ret->data ) );
  XLAL_CHECK_NULL ( ret->data != NULL, XLAL_ENOMEM, "Failed to XLALCalloc(%d,%d)\n", numDetectors, len );

  // loop over detectors
  for ( UINT4 X = 0; X < numDetectors; X ++ )
    {
      ret->data[X] = XLALGetSSBtimes ( multiDetStates->data[X], skypos, refTime, precision );
      XLAL_CHECK_NULL ( ret->data[X] != NULL, XLAL_EFUNC, "ret->data[%d] = XLALGetSSBtimes() failed with xlalErrno = %d\n", X, xlalErrno );

    } /* for X < numDet */

  return ret;

} /* XLALGetMultiSSBtimes() */

/* ===== Object creation/destruction functions ===== */

/**
 * Destroy a MultiSSBtimes structure.
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


/**
 * Destruction of a ComputeFBuffer *contents*,
 * i.e. the multiSSB and multiAMcoeff, while the
 * buffer-container is not freed (which is why it's passed
 * by value and not by reference...)
 */
void
XLALEmptyComputeFBuffer ( ComputeFBuffer *cfb)
{
  XLALDestroyMultiSSBtimes ( cfb->multiSSB );
  cfb->multiSSB = NULL;
  XLALDestroyMultiSSBtimes ( cfb->multiBinary );
  cfb->multiBinary = NULL;
  XLALDestroyMultiAMCoeffs ( cfb->multiAMcoef );
  cfb->multiAMcoef = NULL;
  XLALDestroyMultiCmplxAMCoeffs ( cfb->multiCmplxAMcoef );
  cfb->multiCmplxAMcoef = NULL;
  return;
} /* XLALDestroyComputeFBuffer() */

/* ===== General internal helper functions ===== */

/**
 * Calculate sin(x) and cos(x) to roughly 1e-7 precision using
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

  static BOOLEAN firstCall = 1;
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
      firstCall = 0;
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
      XLALPrintError("\nFailed numerica in sin_cos_2PI_LUT(): xt = %f not in [0,1)\n\n", xt );
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



/**
 * Parameter-estimation: based on large parts on Yousuke's notes and implemention (in CFSv1),
 * extended for error-estimation.
 */
void
LALEstimatePulsarAmplitudeParams (LALStatus * status,			/**< pointer to LALStatus structure */
				  PulsarCandidate *pulsarParams,  	/**< [out] estimated params {h0,cosi,phi0,psi} plus error-estimates */
				  const Fcomponents *Fstat,	 	/**<  Fstat-components Fa, Fb */
				  const CmplxAntennaPatternMatrix *Mmunu/**<  antenna-pattern A,B,C and normalization S_inv*Tsft */
				  )
{
  REAL8 A1h, A2h, A3h, A4h;
  REAL8 Ad, Bd, Cd, Dd, Ed;
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

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( pulsarParams, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( Fstat, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );
  ASSERT ( Mmunu, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL );

  Ad = Mmunu->Ad;
  Bd = Mmunu->Bd;
  Cd = Mmunu->Cd;
  Ed = Mmunu->Ed;
  Dd = Ad * Bd - Cd * Cd - Ed * Ed;

  normAmu = 2.0 / sqrt(2.0 * Mmunu->Sinv_Tsft);	/* generally *very* small!! */

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
  gsl_vector_set (x_mu, 0,   creal(Fstat->Fa) );	/* x_1 */
  gsl_vector_set (x_mu, 1,   creal(Fstat->Fb) ); 	/* x_2 */
  gsl_vector_set (x_mu, 2, - cimag(Fstat->Fa) );	/* x_3 */
  gsl_vector_set (x_mu, 3, - cimag(Fstat->Fb) );	/* x_4 */

  /* ----- fill matrix M^{mu,nu} [symmetric: use UPPER HALF ONLY!!]*/
  gsl_matrix_set (M_Mu_Nu, 0, 0,   Bd / Dd );
  gsl_matrix_set (M_Mu_Nu, 1, 1,   Ad / Dd );
  gsl_matrix_set (M_Mu_Nu, 0, 1, - Cd / Dd );

  gsl_matrix_set (M_Mu_Nu, 0, 3, - Ed / Dd );
  gsl_matrix_set (M_Mu_Nu, 1, 2,   Ed / Dd );

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

  if ( ( fabs( (A1check - A1h)/A1h ) > tolerance ) ||
       ( fabs( (A2check - A2h)/A2h ) > tolerance ) ||
       ( fabs( (A3check - A3h)/A3h ) > tolerance ) ||
       ( fabs( (A4check - A4h)/A4h ) > tolerance ) )
    {
      if ( lalDebugLevel )
	XLALPrintError ( "WARNING LALEstimatePulsarAmplitudeParams(): Difference between estimated and reconstructed Amu exceeds tolerance of %g\n",
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
  /* propagate initial-phase from Fstat-reference-time to refTime of Doppler-params */
  if (XLALExtrapolatePulsarPhase (&phi0, pulsarParams->Doppler.fkdot, pulsarParams->Doppler.refTime, phi0, Fstat->refTime ) != XLAL_SUCCESS) {
    ABORT (status, COMPUTEFSTATC_EXLAL, COMPUTEFSTATC_MSGEXLAL );
  }

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
  /* return also the full Amplitude-Fisher matrix: invert Jh_Mu_nu */
  TRYGSL( gsl_linalg_LU_decomp (Jh_Mu_nu, permh, &signum ), status);
  TRYGSL(gsl_linalg_LU_invert  (Jh_Mu_nu, permh, tmp ), status);
  pulsarParams->AmpFisherMatrix = tmp;

  /* ----- free GSL memory ----- */
  gsl_vector_free ( x_mu );
  gsl_vector_free ( A_Mu );
  gsl_matrix_free ( M_Mu_Nu );
  gsl_matrix_free ( Jh_Mu_nu );
  gsl_permutation_free ( permh );
  gsl_matrix_free ( tmp2 );

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LALEstimatePulsarAmplitudeParams() */

/**
 * Function to allocate a 'FstatAtomVector' struct of num timestamps, pre-initialized to zero!
 */
FstatAtomVector *
XLALCreateFstatAtomVector ( UINT4 num )
{
  FstatAtomVector *ret;

  if ( (ret = LALCalloc ( 1, sizeof(*ret))) == NULL )
    goto failed;

  ret->length = num;

  if ( num == 0 )	/* allow num=0: return just 'head' with NULL arrays */
    return ret;

  if ( (ret->data = LALCalloc ( num, sizeof( *ret->data) ) ) == NULL )
    goto failed;

  return ret;

 failed:
  XLALDestroyFstatAtomVector ( ret );
  XLAL_ERROR_NULL ( XLAL_ENOMEM );

} /* XLALCreateFstatAtomVector() */

/**
 * Function to destroy an FstatAtomVector
 */
void
XLALDestroyFstatAtomVector ( FstatAtomVector *atoms )
{
  if ( !atoms )
    return;

  if ( atoms->data )   LALFree ( atoms->data );
  LALFree ( atoms );

  return;

} /* XLALDestroyFstatAtomVector() */


/**
 * Function to destroy a multi-FstatAtom struct
 */
void
XLALDestroyMultiFstatAtomVector ( MultiFstatAtomVector *multiFstatAtoms )
{
  UINT4 X;
  if ( !multiFstatAtoms)
    return;

  for ( X=0; X < multiFstatAtoms->length; X++ )
    XLALDestroyFstatAtomVector ( multiFstatAtoms->data[X] );

  LALFree ( multiFstatAtoms->data );
  LALFree ( multiFstatAtoms );

  return;

} /* XLALDestroyMultiFstatAtoms() */


/**
 * Convert amplitude-params from 'physical' coordinates {h0, cosi, psi, phi0} into
 * 'canonical' coordinates A^mu = {A1, A2, A3, A4}. The equations are found in
 * \cite JKS98 or \cite Prix07 Eq.(2).
 *
 */
int
XLALAmplitudeParams2Vect ( PulsarAmplitudeVect A_Mu,		/**< [out] canonical amplitude coordinates A^mu = {A1, A2, A3, A4} */
                           const PulsarAmplitudeParams Amp	/**< [in] 'physical' amplitude params {h0, cosi, psi, phi0} */
                           )
{
  REAL8 aPlus = 0.5 * Amp.h0 * ( 1.0 + SQ(Amp.cosi) );
  REAL8 aCross = Amp.h0 * Amp.cosi;
  REAL8 cos2psi = cos ( 2.0 * Amp.psi );
  REAL8 sin2psi = sin ( 2.0 * Amp.psi );
  REAL8 cosphi0 = cos ( Amp.phi0 );
  REAL8 sinphi0 = sin ( Amp.phi0 );

  A_Mu[0] =  aPlus * cos2psi * cosphi0 - aCross * sin2psi * sinphi0;
  A_Mu[1] =  aPlus * sin2psi * cosphi0 + aCross * cos2psi * sinphi0;
  A_Mu[2] = -aPlus * cos2psi * sinphi0 - aCross * sin2psi * cosphi0;
  A_Mu[3] = -aPlus * sin2psi * sinphi0 + aCross * cos2psi * cosphi0;

  return XLAL_SUCCESS;

} /* XLALAmplitudeParams2Vect() */


/**
 * Compute amplitude params \f$A^{\tilde{\mu}} = \{h_0,cosi,\psi,\phi_0\}\f$ from amplitude-vector \f$A^\mu\f$
 * Adapted from algorithm in LALEstimatePulsarAmplitudeParams().
 */
int
XLALAmplitudeVect2Params ( PulsarAmplitudeParams *Amp,	  /**< [out] output physical amplitude parameters {h0,cosi,psi,phi} */
                           const PulsarAmplitudeVect A_Mu /**< [in] input canonical amplitude vector A^mu = {A1,A2,A3,A4} */
                           )
{
  REAL8 h0Ret, cosiRet, psiRet, phi0Ret;

  REAL8 A1, A2, A3, A4, Asq, Da, disc;
  REAL8 Ap2, Ac2, aPlus, aCross;
  REAL8 beta, b1, b2, b3;

  if ( !A_Mu ) {
    XLALPrintError ( "%s: Invalid NULL input vector A_Mu\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }
  if ( !Amp ) {
    XLALPrintError ("%s: invalid NULL input Amp.\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  A1 = A_Mu[0];
  A2 = A_Mu[1];
  A3 = A_Mu[2];
  A4 = A_Mu[3];

  Asq = SQ(A1) + SQ(A2) + SQ(A3) + SQ(A4);
  Da = A1 * A4 - A2 * A3;

  disc = sqrt ( SQ(Asq) - 4.0 * SQ(Da) );

  Ap2  = 0.5 * ( Asq + disc );
  aPlus = sqrt(Ap2);

  Ac2 = 0.5 * ( Asq - disc );
  aCross = MYSIGN(Da) * sqrt( Ac2 );

  beta = aCross / aPlus;

  b1 =   A4 - beta * A1;
  b2 =   A3 + beta * A2;
  b3 = - A1 + beta * A4 ;

  /* amplitude params in LIGO conventions */
  psiRet  = 0.5 * atan2 ( b1,  b2 );  /* [-pi/2,pi/2] */
  phi0Ret =       atan2 ( b2,  b3 );  /* [-pi, pi] */

  /* Fix remaining sign-ambiguity by checking sign of reconstructed A1 */
  {
    REAL8 A1check = aPlus * cos(phi0Ret) * cos(2.0*psiRet) - aCross * sin(phi0Ret) * sin(2*psiRet);
    if ( A1check * A1 < 0 )
      phi0Ret += LAL_PI;
  }

  h0Ret = aPlus + sqrt ( disc );
  cosiRet = aCross / h0Ret;

  /* make unique by fixing the gauge to be psi in [-pi/4, pi/4], phi0 in [0, 2*pi] */
  while ( psiRet > LAL_PI_4 )
    {
      psiRet  -= LAL_PI_2;
      phi0Ret -= LAL_PI;
    }
  while ( psiRet < - LAL_PI_4 )
    {
      psiRet  += LAL_PI_2;
      phi0Ret += LAL_PI;
    }
  while ( phi0Ret < 0 )
    {
      phi0Ret += LAL_TWOPI;
    }

  while ( phi0Ret > LAL_TWOPI )
    {
      phi0Ret -= LAL_TWOPI;
    }

  /* Return final answer */
  Amp->h0   = h0Ret;
  Amp->cosi = cosiRet;
  Amp->psi  = psiRet;
  Amp->phi0 = phi0Ret;

  return XLAL_SUCCESS;

} /* XLALAmplitudeVect2Params() */
