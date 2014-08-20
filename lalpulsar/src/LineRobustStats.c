/*
 *  Copyright (C) 2011-2014 David Keitel
 *  Copyright (C) 2014 Reinhard Prix
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

//---------- INCLUDES ----------
#include <lal/LineRobustStats.h>

//---------- local DEFINES ----------

//----- Macros -----

// ---------- internal types ----------
// ----- module-internal global variables ----------

// opaque internal setup structure for XLALComputeLRstat()
struct tagLRstatSetup {
  UINT4 numDetectors;
  REAL4 Fstar0;
  REAL4 oLGX[PULSAR_MAX_DETECTORS];
  REAL4 C;				// C  = Fstar0 + log(1-pL)
  REAL4 CX[PULSAR_MAX_DETECTORS];	// CX = log( pL rX / Ndet )
};

/*==================== FUNCTION DEFINITIONS ====================*/

/**
 * \f[
 * \newcommand{\Ndet}{N_{\mathrm{det}}}
 * \newcommand{\F}{\mathcal{F}}
 * \newcommand{\Ftho}{\F_*^{(0)}}
 * \newcommand{\oLGX}{o_{\mathrm{LG}}^X}
 * \f]
 * Pre-compute the 'setup' for XLALComputeLRstat() for given prior parameters \f$\Ftho\f$ and \f$\{\oLGX\}\f$.
 */
LRstatSetup *
XLALCreateLRstatSetup ( const UINT4 numDetectors,			//!< [in] number of detectors \f$\Ndet\f$
                        const REAL4 Fstar0,				//!< [in] prior parameter \f$\Ftho\f$
                        const REAL4 oLGX[PULSAR_MAX_DETECTORS]		//!< [in] prior per-detector line odds \f$\{\oLGX\}\f$, if NULL: interpreted as \f$\oLGX=\frac{1}{\Ndet} \forall X\f$
                        )
{
  // check input
  XLAL_CHECK_NULL ( (numDetectors >= 2) && (numDetectors <= PULSAR_MAX_DETECTORS), XLAL_EDOM );
  for ( UINT4 X = 0; X < numDetectors; X ++ ) {
    XLAL_CHECK_NULL ( (oLGX == NULL) || (oLGX[X] >= 0), XLAL_EDOM );
  }

  // create setup struct
  LRstatSetup *setup;
  XLAL_CHECK_NULL ( (setup = XLALCalloc( 1, sizeof(*setup) )) != NULL, XLAL_ENOMEM );
  setup->numDetectors = numDetectors;
  setup->Fstar0       = Fstar0;

  // Compute oLG from Eq.(22)
  REAL4 oLG = 0.0;
  REAL4 oLGX_default = 1.0/numDetectors;  // 'even odds': use oLGX=1/numDetectors for all X ==> oLG=1, rX=1, pL=0.5
  for ( UINT4 X = 0; X < numDetectors; X ++ )
    {
      setup->oLGX[X] = (oLGX != NULL) ? oLGX[X] : oLGX_default;
      oLG += setup->oLGX[X];
    }

  // constant transition scale term
  REAL4 log_1_plus_oLG = log ( 1.0 + oLG ); 	// ln(1+oLG)
  setup->C = Fstar0 - log_1_plus_oLG; 		// Fstar0 + ln(1-pL) = Fstar0 - ln(1+oLG)

  for (UINT4 X = 0; X < numDetectors; X++)
    {
      setup->CX[X] = -LAL_REAL4_MAX; // fallback if oLGX==0 to avoid raising underflow
      if ( setup->oLGX[X] > 0 ) {
        setup->CX[X] = log ( setup->oLGX[X] ) - log_1_plus_oLG; // ln(oLGX) - log(1+oLG)
      }
    } // for X < numDetectors

  return setup;

} // XLALCreateLRstatSetup()


/**
 * \f[
 * \newcommand{\F}{\mathcal{F}}
 * \newcommand{\Ftho}{\F_*^{(0)}}
 * \newcommand{\FpMax}{\F'_{\mathrm{max}}}
 * \newcommand{\Ndet}{N_{\mathrm{det}}}
 * \newcommand{\Signal}{\mathrm{S}}
 * \newcommand{\Gauss}{\mathrm{G}}
 * \newcommand{\Line}{\mathrm{L}}
 * \newcommand{\SGL}{{\Signal\Gauss\Line}}
 * \newcommand{\oLG}{o_{\Line\Gauss}}
 * \newcommand{\oLGX}{\oLG^X}
 * \newcommand{\OSGL}{O_\SGL}
 * \newcommand{\oSGL}{o_\SGL}
 * \newcommand{\data}{\mathrm{data}}
 * \newcommand{\pL}{p_\Line}
 * \newcommand{\logten}{\log_{10}}
 * \f]
 * Compute the "line-robust" statistic 'LRstat' from multi-detector \f$2\F\f$ and single-detector \f$\{2\F^X\}\f$ 'Fstat' values
 * (coherent or semi-coherent sum) and a pre-computed 'setup' from XLALCreateLRstatSetup().
 *
 * The line-robust statistic is defined as \f$\mathrm{LRstat} = \logten B_\SGL\f$, where \f$B_\SGL\f$ is the BayesFactor[Signal-vs-Gauss_OR_Line], i.e.
 * \f$B_\SGL \equiv \frac{P(\data | \text{Signal})}{P(\data | \text{Gauss or Line})}\f$, which is related to the odds ratio
 * via \f$\OSGL = B_\SGL\,\oSGL\f$ in terms of the prior odds \f$\oSGL\f$.
 *
 * Here we use the odds ratio derived in Eq.(36) (and Eq.(55)) of \cite KPPLS2014, from which we can obtain
 * \f{equation}{
 * \ln B_\SGL = \F - \FpMax - \ln\left( e^{C - \FpMax} + \sum_X e^{\F^X + C^X - \FpMax} \right) \,,
 * \f}
 * where \c useLogCorrection controls whether or not to include the (usually small) log-correction of the last term,
 * and we defined
 * \f{equation}{ \FpMax \equiv \max\left[ C, \{ \F^X + C^X \} \right] \f}
 * and \f$C,\{C^X\}\f$ are the prior quantities pre-computed in XLALCreateLRstatSetup(), namely
 * \f{align}{
 * C &\equiv \Ftho + \ln(1-\pL) = \Ftho - \ln(1+\oLG)\,,\\
 * C^X &\equiv \ln\frac{\pL\,r^X}{\Ndet} = \ln \oLGX - \ln(1+\oLG)\,,
 * \f}
 * and we used the following relations from \cite KPPLS2014:
 * \f$ \oLG = \sum_X \oLGX \;\f$ [Eq.(22)],
 * \f$ r^X = \frac{\oLGX\, \Ndet}{\oLG} \;\f$ [Eq.(23)],
 * and \f$\pL = \frac{\oLG}{1 + \oLG} \;\f$ [Eq.(37)].
 *
 * \return NOTE: return is \f$\logten B_\SGL = \ln B_\SGL \, \logten e\f$
 *
 */
REAL4
XLALComputeLRstat ( const REAL4 twoF,				//!< [in] multi-detector Fstat \f$2\F\f$ (coherent or semi-coherent sum(!))
                    const REAL4 twoFX[PULSAR_MAX_DETECTORS],	//!< [in] per-detector Fstats \f$\{2\F^X\}\f$ (coherent or semi-coherent sum(!))
                    const LRstatSetup *setup,			//!< [in] pre-computed setup from XLALCreateLRstatSetup()
                    const BOOLEAN useLogCorrection		//!< [in] include log-term correction [slower] or not [faster, less accurate]
                    )
{
  XLAL_CHECK ( setup != NULL, XLAL_EINVAL );

  REAL4 FpMax = setup->C; // used to keep track of log of maximal denominator sum-term

  // per-detector contributions, including line weights
  REAL4 Xterm[PULSAR_MAX_DETECTORS];
  for ( UINT4 X=0; X < setup->numDetectors; X ++ )
    {
      Xterm[X] = 0.5 * twoFX[X] + setup->CX[X]; 	// FX + CX
      FpMax = fmax ( FpMax, Xterm[X] );
  }

  REAL4 ln_BSGL = 0.5 * twoF - FpMax; // approximate result without log-correction term

  if ( !useLogCorrection ) {
    return ln_BSGL * LAL_LOG10E; // return log10(B_SGL)
  }

  // if useLogCorrection: extraSum = e^(Fstar0+ln(1-pL) - FpMax) + sum_X e^( Xterm[X] - FpMax )
  REAL4 extraSum = exp ( setup->C  - FpMax );

  // ... and add all FX-contributions
  for ( UINT4 X = 0; X < setup->numDetectors; X++ )
    {
      extraSum += exp ( Xterm[X] - FpMax );
    }

  ln_BSGL -= log ( extraSum ); // F - FpMax - ln( ... )

  return ln_BSGL * LAL_LOG10E; // return log10(B_SGL)

} // XLALComputeLRstat()


/**
 * Parse string-vectors (typically input by user) of N per-detector line-to-Gaussian prior ratios
 * \f$ o_{LG}^X = \frac{P(H_L^X)}{P(H_G^X)} \f$ to a flat array of REAL4s.
 */
int
XLALParseLinePriors ( REAL4 oLGX[PULSAR_MAX_DETECTORS],		//!< [out] array of parsed \f$ o_{LG}^X \f$ values
                      const LALStringVector *oLGX_string	//!< [in] string-list of \f$ o_{LG}^X \f$ values
                    )
{
  XLAL_CHECK ( oLGX != NULL, XLAL_EINVAL );
  XLAL_CHECK ( oLGX_string != NULL, XLAL_EINVAL );
  UINT4 numDetectors = oLGX_string->length;
  XLAL_CHECK ( (numDetectors > 0) && (numDetectors <= PULSAR_MAX_DETECTORS), XLAL_EINVAL );

  // parse input string
  for ( UINT4 X = 0; X < numDetectors; X ++ )
    {
      XLAL_CHECK ( XLALParseStringValueToREAL4 ( &oLGX[X], oLGX_string->data[X] ) == XLAL_SUCCESS, XLAL_EFUNC );
    } // for X < numDetectors

  return XLAL_SUCCESS;

} // XLALParseLinePriors()


// ---------- OLD API functions [deprecated!] ------------------------------

/**
 * Deprecated function to compute Line Veto statistics from multi- and single-detector \f$ \mathcal{F} \f$-stats,
 * using outdated \f$ \rho \f$ notation and prior normalization.
 * This is not the log-Bayes-factor!
 * \f$ \mathrm{LV} = \mathcal{F} - \log \left( \frac{\rho_{\mathrm{max,line}}^4}{70} + \sum_X l^X e^{\mathcal{F}^X} \right) \f$
 *
 * \deprecated use XLALComputeLRstat() instead.
 *
 * Also this is just a wrapper for XLALComputeLineVetoArray, which is faster for many LV values at identical priors.
 * This function here just translates REAL4Vectors to fixed REAL4 arrays,
 * and linear priors rhomaxline and lX to logarithmic values.
 */
REAL4 XLALComputeLineVeto ( const REAL4 TwoF,			/**< multi-detector  \f$ \mathcal{F} \f$-stat */
                            const REAL4Vector *TwoFXvec,	/**< vector of single-detector \f$ \mathcal{F} \f$-stats */
                            const REAL8 rhomaxline,		/**< amplitude prior normalization for lines */
                            const REAL8Vector *lXvec,		/**< vector of single-detector prior line odds ratio, default to \f$ l^X=1 \f$ for all \f$ X \f$ if NULL */
                            const BOOLEAN useAllTerms		/**< only use leading term (FALSE) or all terms (TRUE) in log sum exp formula? */
                          )
{

  XLAL_PRINT_DEPRECATION_WARNING("XLALComputeLRstat");

  /* check input parameters and report errors */
  XLAL_CHECK_REAL4 ( TwoF && TwoFXvec && TwoFXvec->data, XLAL_EFAULT, "Empty TwoF or TwoFX pointer as input parameter." );
  UINT4 numDetectors = TwoFXvec->length;
  XLAL_CHECK_REAL4 ( !lXvec || ( lXvec->length == numDetectors ), XLAL_EBADLEN, "Input lX (%d) and TwoFX (%d) vectors have different lengths.", lXvec->length, numDetectors );

  XLAL_CHECK_REAL4 ( rhomaxline >= 0, XLAL_EDOM, "Negative prior range 'rhomaxline' = %g. Must be >= 0!", rhomaxline );
  REAL8 logRhoTerm = 0.0;
  if ( rhomaxline > 0.0 ) {
   logRhoTerm = 4.0 * log(rhomaxline) - log(70.0);
  }
  else { /* if rhomaxline == 0.0, logRhoTerm should become irrelevant in summation */
    logRhoTerm = - LAL_REAL8_MAX;
  }

  REAL8 *loglX = NULL;
  REAL8 loglXtemp[numDetectors];
  if ( lXvec ) {
    for (UINT4 X = 0; X < numDetectors; X++) {
      if ( lXvec->data[X] > 0 ) {
        loglXtemp[X] = log(lXvec->data[X]);
      }
      else if ( lXvec->data[X] == 0 ) { /* if zero prior ratio, approximate log(0)=-inf by -LAL_REAL4_MAX to avoid raising underflow exceptions */
        loglXtemp[X] = - LAL_REAL8_MAX;
      }
      else { /* negative prior ratio is a mistake! */
       XLAL_ERROR_REAL4 ( XLAL_EDOM, "Negative input prior-ratio for detector X=%d: lX[X]=%g\n", X, lXvec->data[X] );
      }
    }  /* for X < numDetectors */
    loglX = loglXtemp;
  } /* if lXvec */

  REAL4 LV = XLALComputeLineVetoArray ( TwoF, numDetectors, TwoFXvec->data, logRhoTerm, loglX, useAllTerms );

  return LV;

} /* XLALComputeLineVeto() */


/**
 * Deprecated function to compute Line Veto statistics from multi- and single-detector \f$ \mathcal{F} \f$-stats,
 * using outdated \f$ \rho \f$ notation and prior normalization.
 * This is not the log-Bayes-factor!
 * \f$ \mathrm{LV} = \mathcal{F} - \log \left( \frac{\rho_{\mathrm{max,line}}^4}{70} + \sum_X l^X e^{\mathcal{F}^X} \right) \f$
 *
 * \deprecated use XLALComputeLRstat() instead.
 *
 * Implemented by log sum exp formula:
 * \f$ \mathrm{LV} = \mathcal{F} - \max(\mathrm{denom. terms}) - \log \left( \sum e^{\mathrm{denom. terms}-\max} \right) \f$.
 *
 * From the analytical derivation, there should be an extra term \f$ + o_{SN} + 4 \log \left( \frac{\rho_{\mathrm{max,line}}}{\rho_{\mathrm{max,sig}}} \right) \f$,
 * but this is irrelevant for toplist sorting, only a normalization which can be replaced arbitrarily.
 *
 * NOTE: priors logRhoTerm, loglX have to be logarithmized already.
 */
REAL4
XLALComputeLineVetoArray ( const REAL4 TwoF,		/**< multi-detector \f$ \mathcal{F} \f$-stat */
                           const UINT4 numDetectors,	/**< number of detectors */
                           const REAL4 *TwoFX,		/**< array of single-detector \f$ \mathcal{F} \f$-stats */
                           const REAL8 logRhoTerm,	/**< extra term coming from prior normalization: \f$ \log \left( \frac{\rho_{\mathrm{max,line}}^4}{70} \right) \f$ */
                           const REAL8 *loglX,		/**< array of logs of single-detector prior line odds ratios, default to \f$ \log(l^X)=\log(1)=0 \f$ for all \f$ X \f$ if NULL */
                           const BOOLEAN useAllTerms	/**< only use leading term (FALSE) or all terms (TRUE) in log sum exp formula? */
                         )
{

  XLAL_PRINT_DEPRECATION_WARNING("XLALComputeLRstat");

  /* check input parameters and report errors */
  XLAL_CHECK_REAL4 ( TwoF && TwoFX, XLAL_EFAULT, "Empty TwoF or TwoFX pointer as input parameter." );

  /* set up temporary variables and structs */
  REAL4 log0  = - LAL_REAL4_MAX;	/* approximates -inf */

  REAL4 maxInSum = log0;           /* keep track of largest summand in denominator, for logsumexp formula below */
  REAL4 FXprior[numDetectors];     /* FXprior equiv log(lX * e^(FX)) = FX + loglX */

  for (UINT4 X = 0; X < numDetectors; X++) {
    FXprior[X] = 0.5 * TwoFX[X];
    if (loglX) { /* if no priors given, just use lX=1 => loglX=0 for all X => do not add anything */
      FXprior[X] += loglX[X];
    }
    /* keep track of maximum value in denominator sum  */
    if ( FXprior[X] > maxInSum ) {
      maxInSum = FXprior[X];
    }
  } /* for X < numDetectors */

  /* special treatment for additional denominator term 'rho^4/70' */
  if ( logRhoTerm > maxInSum ) {
    maxInSum = logRhoTerm;
  }

  REAL4 LV = 0.0;	/* output variable for Line Veto statistics */

  LV = 0.5 * TwoF - maxInSum;	/* dominant term to LV-statistic */

  if ( useAllTerms ) { /* optionally add logsumexp term (possibly negligible in many cases) */

    REAL4 extraSum=0;	/* will be:  e^[-(maxInSum - logRhoTerm)] + sum_X e^[ -(maxInSum - FXprior) ] >= 1 */

    /* need to treat (rho^4/70) term separately */
    extraSum += exp ( logRhoTerm - maxInSum );

    /* now add all FX-contributions */
    for (UINT4 X = 0; X < numDetectors; X++) {
      extraSum += exp ( FXprior[X] - maxInSum );
    }

    LV -= log ( extraSum );

  } /* if useAllTerms */

  return LV;

} /* XLALComputeLineVetoArray() */
