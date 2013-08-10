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
#include <lal/CWFastMath.h>

/* GSL includes */
#include <lal/LALGSL.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_roots.h>

#include <lal/AVFactories.h>
#include "ComputeFstat.h"
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
const Fcomponents empty_Fcomponents;
const ComputeFParams empty_ComputeFParams;
const ComputeFBuffer empty_ComputeFBuffer;

/*---------- internal prototypes ----------*/

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
	if ( XLALSinCos2PiLUT ( &imagQ, &realQ, - lambda_alpha ) ) {
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
      XLALSinCos2PiLUT ( &s_alpha, &c_alpha, kappa_star );
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
	if ( XLALSinCos2PiLUT ( &imagQ, &realQ, - lambda_alpha ) ) {
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
      XLALSinCos2PiLUT ( &s_alpha, &c_alpha, kappa_star );
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
	if ( XLALSinCos2PiLUT ( &imagQ, &realQ, - lambda_alpha ) ) {
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
