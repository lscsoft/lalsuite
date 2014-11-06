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

// opaque internal setup structure for XLALComputeBSGL()
struct tagBSGLSetup {
  UINT4 numDetectors;
  REAL4 Fstar0;
  REAL4 oLGX[PULSAR_MAX_DETECTORS];
  REAL4 C;				// C  = Fstar0 + log(1-pL)
  REAL4 CX[PULSAR_MAX_DETECTORS];	// CX = log( pL rX / Ndet )
  BOOLEAN useLogCorrection;
};

/*==================== FUNCTION DEFINITIONS ====================*/

/**
 * \f[
 * \newcommand{\Ndet}{N_{\mathrm{det}}}
 * \newcommand{\F}{\mathcal{F}}
 * \newcommand{\Ftho}{\F_*^{(0)}}
 * \newcommand{\oLGX}{o_{\mathrm{LG}}^X}
 * \f]
 * Pre-compute the 'setup' for XLALComputeBSGL() for given prior parameters \f$\Ftho\f$ and \f$\{\oLGX\}\f$.
 */
BSGLSetup *
XLALCreateBSGLSetup ( const UINT4 numDetectors,			//!< [in] number of detectors \f$\Ndet\f$
                      const REAL4 Fstar0,				//!< [in] prior parameter \f$\Ftho\f$
                      const REAL4 oLGX[PULSAR_MAX_DETECTORS],		//!< [in] prior per-detector line odds \f$\{\oLGX\}\f$, if NULL: interpreted as \f$\oLGX=\frac{1}{\Ndet} \forall X\f$
                      const BOOLEAN useLogCorrection			//!< [in] include log-term correction [slower] or not [faster, less accurate]
                      )
{
  // check input
  XLAL_CHECK_NULL ( (numDetectors >= 2) && (numDetectors <= PULSAR_MAX_DETECTORS), XLAL_EDOM, "numDetectors = %d not within [2,%d]\n", numDetectors, PULSAR_MAX_DETECTORS );
  for ( UINT4 X = 0; X < numDetectors; X ++ ) {
    XLAL_CHECK_NULL ( (oLGX == NULL) || (oLGX[X] >= 0), XLAL_EDOM );
  }

  // create setup struct
  BSGLSetup *setup;
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

  setup->useLogCorrection = useLogCorrection;

  return setup;

} // XLALCreateBSGLSetup()

void
XLALDestroyBSGLSetup ( BSGLSetup * setup )
{
  XLALFree ( setup );
  return;
} // XLALDestroyBSGLSetup()

/**
 * \f[
 * \newcommand{\F}{\mathcal{F}}
 * \newcommand{\Ftho}{\F_*^{(0)}}
 * \newcommand{\FpMax}{\F^{\#}_{\mathrm{max}}}
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
 * \newcommand{\BSGL}{B_{\SGL}}
 * \f]
 * Compute the 'line-robust statistic' \f$\logten \BSGL\f$ from multi-detector \f$2\F\f$ and single-detector \f$\{2\F^X\}\f$ 'F-stat values'
 * (coherent or semi-coherent sum) and a pre-computed 'setup' from XLALCreateBSGLSetup().
 *
 * \f$\BSGL\f$ is the BayesFactor[Signal-vs-Gauss_OR_Line], i.e.
 * \f$\BSGL \equiv \frac{P(\data | \text{Signal})}{P(\data | \text{Gauss or Line})}\f$, which is related to the odds ratio
 * via \f$\OSGL = \oSGL\,\BSGL\f$ in terms of the prior odds \f$\oSGL\f$.
 *
 * Here we use the odds ratio derived in Eq.(36) (and Eq.(55)) of \cite KPPLS2014, from which we can obtain
 * \f{equation}{
 * \ln \BSGL = \F - \FpMax - \ln\left( e^{C - \FpMax} + \sum_X e^{\F^X + C^X - \FpMax} \right) \,,
 * \f}
 * where \c setup->useLogCorrection controls whether or not to include the (usually small) log-correction of the last term,
 * and we defined
 * \f{equation}{ \FpMax \equiv \max\left[ C, \{ \F^X + C^X \} \right] \f}
 * and \f$C,\{C^X\}\f$ are the prior quantities pre-computed in XLALCreateBSGLSetup(), namely
 * \f{align}{
 * C &\equiv \Ftho + \ln(1-\pL) = \Ftho - \ln(1+\oLG)\,,\\
 * C^X &\equiv \ln\frac{\pL\,r^X}{\Ndet} = \ln \oLGX - \ln(1+\oLG)\,,
 * \f}
 * and we used the following relations from \cite KPPLS2014:
 * \f$ \oLG = \sum_X \oLGX \;\f$ [Eq.(22)],
 * \f$ r^X = \frac{\oLGX\, \Ndet}{\oLG} \;\f$ [Eq.(23)],
 * and \f$\pL = \frac{\oLG}{1 + \oLG} \;\f$ [Eq.(37)].
 *
 * Here, \f$\FpMax\f$ differs from \f$\F''_{\mathrm{max}}\f$ from [Eq.(41)] by \f$\ln\Ndet\f$, due to summation instead of an average over detectors.
 *
 * \return NOTE: return is \f$\logten B_\SGL = \ln B_\SGL \, \logten e\f$
 *
 */
REAL4
XLALComputeBSGL ( const REAL4 twoF,				//!< [in] multi-detector Fstat \f$2\F\f$ (coherent or semi-coherent sum(!))
                  const REAL4 twoFX[PULSAR_MAX_DETECTORS],	//!< [in] per-detector Fstats \f$\{2\F^X\}\f$ (coherent or semi-coherent sum(!))
                  const BSGLSetup *setup			//!< [in] pre-computed setup from XLALCreateBSGLSetup()
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

  if ( !setup->useLogCorrection ) {
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

} // XLALComputeBSGL()


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
