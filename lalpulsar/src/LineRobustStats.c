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
#include <lal/UserInputParse.h>

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
  REAL4 C;						// constant term C  = Fstar0 + log(1-pL)
  REAL4 CX[PULSAR_MAX_DETECTORS];			// per-detector term CX = log( pL rX / Ndet )
//   REAL4 CXl[PULSAR_MAX_DETECTORS][PULSAR_MAX_SEGMENTS];	// per-detector, per-segment term ClX = (Nseg - 1)*Fstar0 + log(  ptLXl ) FIXME: not implemented in XLALCreateBSGLSetup()
  REAL4 perSegTerm;					// extra term for per-segment transient contributions, (Nseg - 1)*Fstar0
  BOOLEAN useLogCorrection;
};

/*==================== FUNCTION DEFINITIONS ====================*/

/**
 * \f[
 * \newcommand{\Ndet}{N_{\mathrm{det}}}
 * \newcommand{\F}{\mathcal{F}}
 * \newcommand{\scF}{\hat{\F}}
 * \newcommand{\Ftho}{\F_*^{(0)}}
 * \newcommand{\scFtho}{\scF_*^{(0)}}
 * \newcommand{\oLGX}{o_{\mathrm{LG}}^X}
 * \newcommand{\Nseg}{N_{\mathrm{seg}}}
 * \f]
 * Pre-compute the 'setup' for XLALComputeBSGL() for given prior parameters \f$\Ftho\f$ and \f$\{\oLGX\}\f$.
 */
BSGLSetup *
XLALCreateBSGLSetup ( const UINT4 numDetectors,			//!< [in] number of detectors \f$\Ndet\f$
                      const REAL4 Fstar0,			//!< [in] prior transition-scale parameter: pass \f$\Ftho\f$ in coherent case, \f$\scFtho=\Nseg\Ftho\f$ in semi-coherent case
                      const REAL4 oLGX[PULSAR_MAX_DETECTORS],	//!< [in] prior per-detector line odds \f$\{\oLGX\}\f$, if NULL: interpreted as \f$\oLGX=\frac{1}{\Ndet} \forall X\f$
                      const BOOLEAN useLogCorrection,		//!< [in] include log-term correction [slower] or not [faster, less accurate]
                      const UINT4 numSegments			//!< [in] number of segments \f$\Nseg\f$ in a semi-coherent search (=1 for coherent search)
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

  setup->perSegTerm = Fstar0 * ( numSegments - 1 ) / numSegments; // in this function's convention, Fstar0 is actually semi-coherent Fstar0hat, so have to divide by numSegments here!

  // FIXME: implement CXl[X][l] = (numSegments - 1)*Fstar0 + log(  ptLXl[X][l] )

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
 * Here we use the odds ratio derived in Eq.(36) (and Eq.(55)) of \cite KPPLS2014 , from which we can obtain
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
 * and we used the following relations from \cite KPPLS2014 :
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
XLALComputeBSGL ( const REAL4 twoF,				//!< [in] multi-detector F-stat \f$2\F\f$ (coherent or semi-coherent sum(!))
                  const REAL4 twoFX[PULSAR_MAX_DETECTORS],	//!< [in] per-detector F-stats \f$\{2\F^X\}\f$ (coherent or semi-coherent sum(!))
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
 * \f[
 * \newcommand{\coh}[1]{\widetilde{#1}}
 * \newcommand{\sc}[1]{\widehat{#1}}
 * \newcommand{\cohF}{\coh{\F}}
 * \newcommand{\scF}{\sc{\F}}
 * \newcommand{\denomterm}{\sc{D}}
 * \newcommand{\denommax}{\denomterm_{\mathrm{max}}}
 * \newcommand{\Transline}{{\mathrm{t\Line}}}
 * \newcommand{\cppLX}{\sc{p}_{\Line}^X}
 * \newcommand{\cppTLXk}{\coh{p}_{\Transline}^{X\ell}}
 * \newcommand{\cppLTLsc}{\sc{p}_{\Line\Transline}}
 * \f]
 * Compute the common noise-hypothesis denominator for various semi-coherent 'transient-line-robust statistics'
 * as defined in \cite Keitel2015 .
 *
 * NOTE: In the current implementation, only
 * \c setup->useLogCorrection == TRUE
 * is accepted!
 * In this case, the output is simply the maximum of the set of denominator exponents,
 * \f{equation}{
 * \denommax \equiv \max \left\{ \Nseg\Ftho + \ln\left( 1-\cppLTLsc \right) ,\,
 *                               \scF^X + \ln\cppLX ,\,
 *                               \cohF^{X\ell} + (\Nseg-1)\Ftho + \ln\cppTLXk
 *                       \right\} \,,
 * \f}
 * with the simplifying assumption that all
 * \f$\cppTLXk = \cppLTLsc = \cppLX = 1\f$.
 *
 */
REAL4
XLALComputeGLtLDenominator ( const REAL4 twoFX[PULSAR_MAX_DETECTORS],		//!< [in] semi-coherent sums \f$\{2\scF^X\}\f$ of per-detector F-stats
                             const REAL4 maxtwoFXl[PULSAR_MAX_DETECTORS],	//!< [in] maxima \f$\{\max\limits_{\ell}2\{\cohF^{X\ell}\}\}\f$ of per-detector F-stats over segments
                             const BSGLSetup *setup				//!< [in] pre-computed setup from XLALCreateBSGLSetup()
                             )
{

  XLAL_CHECK ( setup != NULL, XLAL_EINVAL );
  XLAL_CHECK ( !setup->useLogCorrection, XLAL_EDOM, "log correction not implemented for GLtL denominator.");

  REAL4 FpMax = setup->C; // used to keep track of log of maximal denominator sum-term

  // per-detector contributions, including line weights
  REAL4 Xterm[PULSAR_MAX_DETECTORS];
  REAL4 Xlterm[PULSAR_MAX_DETECTORS];
  for ( UINT4 X=0; X < setup->numDetectors; X ++ )
    {
      Xterm[X] = 0.5 * twoFX[X] + setup->CX[X]; 	// FX + CX
      FpMax = fmax ( FpMax, Xterm[X] );
//       Xlterm[X] = 0.5 * maxtwoFXl[X] + setup->CXl[X][maxseg[X]]; // FIXME: this would be the general expression
      Xlterm[X] = 0.5 * maxtwoFXl[X] + setup->perSegTerm + setup->CX[X]; // assuming ptLXl = pLX for all segments
      FpMax = fmax ( FpMax, Xlterm[X] );
  }

//   if ( !setup->useLogCorrection ) {
    return FpMax;
//   }

  // if useLogCorrection: extraSum = e^(Fstar0+ln(1-pL) - FpMax) + sum_X e^( Xterm[X] - FpMax )
//   REAL4 extraSum = exp ( setup->C  - FpMax );

  // ... and add all FX, FXl contributions (FIXME: not implemented!)
//   for ( UINT4 X = 0; X < setup->numDetectors; X++ )
//     {
//       extraSum += exp ( Xterm[X] - FpMax );
//     }

//   return = FpMax + log ( extraSum );

}

/**
 * \f[
 * \newcommand{\denomtermset}{\mathcal{\denomterm}}
 * \newcommand{\SGLtL}{{\Signal\Gauss\Line\Transline}}
 * \newcommand{\OSGLtL}{O_\SGLtL}
 * \newcommand{\oSGLtL}{o_\SGLtL}
 * \newcommand{\BSGLtL}{B_{\SGLtL}}
 * \f]
 * Compute the semi-coherent transient-line-robust CW detection statistic \f$\logten \BSGLtL\f$
 * from \cite Keitel2015 .
 *
 * As defined in Eq.(22), \f$\BSGLtL\f$ is the BayesFactor for CW signals vs.
 * either Gaussian noise, a persistent line or a transient line, i.e.
 * \f$\BSGLtL \equiv \frac{P(\data | \text{Signal})}{P(\data | \text{Gauss or Line or transient Line})}\f$.
 * It is related to the odds ratio
 * via \f$\OSGLtL = \oSGLtL\,\BSGLtL\f$ in terms of the prior odds \f$\oSGLtL\f$.
 *
 * This function returns an approximation of \f$\logten B_\SGLtL = \ln B_\SGLtL \, \logten e\f$,
 * with the general expression for \f$\ln B_\SGLtL\f$ given in Eq.(23) of \cite Keitel2015 :
 * \f{equation}{
 * \ln \BSGLtL = \F - \denommax - \ln\left( \sum\limits_{\denomterm \in \denomtermset} e^{\denomterm-\denommax} \right) \,,
 * \f}
 * while the current implementation of XLALComputeGLtLDenominator() uses only
 * the maximum \f$\denommax\f$ of the set of denominator exponents and drops the additional log-correction term.
 *
 * NOTE: This implementation also assumes unit prior weights for the various noise hypotheses.
 *
 */
REAL4
XLALComputeBSGLtL ( const REAL4 twoF,					//!< [in] semi-coherent sum \f$2\scF\f$ of multi-detector F-stats
                    const REAL4 twoFX[PULSAR_MAX_DETECTORS],		//!< [in] semi-coherent sums \f$\{2\scF^X\}\f$ of per-detector F-stats
                    const REAL4 maxtwoFXl[PULSAR_MAX_DETECTORS],	//!< [in] maxima \f$\{\max\limits_{\ell}2\{\cohF^{X\ell}\}\}\f$ of per-detector F-stats over segments
                    const BSGLSetup *setup				//!< [in] pre-computed setup from XLALCreateBSGLSetup()
                    )
{

  REAL4 GLtLDenominator = XLALComputeGLtLDenominator ( twoFX, maxtwoFXl, setup );

  REAL4 ln_BSGLtL = 0.5 * twoF - GLtLDenominator;

  return ln_BSGLtL * LAL_LOG10E; // return log10(B_SGLtL)

} // XLALComputeBSGLtL

/**
 * \f[
 * \newcommand{\Transsig}{{\mathrm{t\Signal}}}
 * \newcommand{\tSGLtL}{{\Transsig\Gauss\Line\Transline}}
 * \newcommand{\BtSGLtL}{B_{\tSGLtL}}
 * \newcommand{\OtSGLtL}{O_\tSGLtL}
 * \newcommand{\otSGLtL}{o_\tSGLtL}
 * \newcommand{\oTSGk}{\coh{o}_{\Transsig/\Gauss}^{\ell}}
 * \f]
 * Compute a semi-coherent transient-line-robust tCW detection statistic \f$\logten \BtSGLtL\f$
 * based on \cite Keitel2015 .
 *
 * Based on the transient-CW hypothesis as given in Eq.(29), \f$\BtSGLtL\f$ is the BayesFactor for tCW signals vs.
 * either Gaussian noise, a persistent line or a transient line, i.e.
 * \f$\BtSGLtL \equiv \frac{P(\data | \text{Signal})}{P(\data | \text{Gauss or Line or transient Line})}\f$.
 * It is related to the odds ratio
 * via \f$\OtSGLtL = \otSGLtL\,\BtSGLtL\f$ in terms of the prior odds \f$\oSGLtL\f$.
 *
 * This function returns an approximation of \f$\logten B_\SGLtL = \ln B_\SGLtL \, \logten e\f$.
 * While the general expression for \f$\ln B_\SGLtL\f$ would be analoguous to Eq.(23) of \cite Keitel2015 ,
 * \f{equation}{
 * \ln \BtSGLtL = \sum\limits_{\ell=1}^{\Nseg} \oTSGk \, e^{\cohF^{\ell} - \Ftho} - \denommax - \ln\left( \sum\limits_{\denomterm \in \denomtermset} e^{\denomterm-\denommax} \right) \,,
 * \f}
 * the current implementation uses only
 * the maximum multi-detector, single-segment F-stat value \f$\max\cohF^{\ell}\f$ for the numerator
 * and, through XLALComputeGLtLDenominator(), the maximum \f$\denommax\f$ of the set of denominator exponents.
 * The  additional log-correction terms are dropped.
 *
 * NOTE: This implementation also assumes unit prior weights for the various signal and noise hypotheses.
 *
 */
REAL4
XLALComputeBtSGLtL ( const REAL4 maxtwoFl,				//!< [in] maximum \f$\max\limits_{\ell}2\cohF^\ell\f$ of multi-detector F-stats over segments
                     const REAL4 twoFX[PULSAR_MAX_DETECTORS],		//!< [in] semi-coherent sums \f$\{2\scF^X\}\f$ of per-detector F-stats
                     const REAL4 maxtwoFXl[PULSAR_MAX_DETECTORS],	//!< [in] maxima \f$\{\max\limits_{\ell}2\{\cohF^{X\ell}\}\}\f$ of per-detector F-stats over segments
                     const BSGLSetup *setup				//!< [in] pre-computed setup from XLALCreateBSGLSetup()
                     )
{
  XLAL_CHECK ( setup != NULL, XLAL_EINVAL );

  REAL4 GLtLDenominator = XLALComputeGLtLDenominator ( twoFX, maxtwoFXl, setup );

  REAL4 ln_BtSGLtL = 0.5 * maxtwoFl + setup->perSegTerm - GLtLDenominator; // assuming equal pS=ptSl priors for all segments

  return ln_BtSGLtL * LAL_LOG10E; // return log10(B_SGLtL)

} // XLALComputeBtSGLtL

/**
 * \f[
 * \newcommand{\StSGLtL}{{\Signal\Transsig\Gauss\Line\Transline}}
 * \newcommand{\OStSGLtL}{O_\StSGLtL}
 * \newcommand{\oStSGLtL}{o_\StSGLtL}
 * \newcommand{\BStSGLtL}{B_{\StSGLtL}}
 * \newcommand{\numerterm}{\sc{E}}
 * \newcommand{\numertermset}{\mathcal{\numerterm}}
 * \newcommand{\numermax}{\numerterm_{\mathrm{max}}}
 * \newcommand{\cppS}{\sc{p}_{\Signal}}
 * \newcommand{\cppTS}{\coh{p}_{\Transsig}}
 * \newcommand{\cppTSk}{\coh{p}_{\Transsig}^{\ell}}
 * \f]
 * Compute the semi-coherent transient-line-robust CW+tCW detection statistic \f$\logten \BStSGLtL\f$
 * from \cite Keitel2015 .
 *
 * As defined in Eq.(34), \f$\BStSGLtL\f$ is the BayesFactor for either persistent CW or transient tCW signals vs.
 * either Gaussian noise, a persistent line or a transient line, i.e.
 * \f$\BStSGLtL \equiv \frac{P(\data | \text{Signal or transient Signal})}{P(\data | \text{Gauss or Line or transient Line})}\f$.
 * It is related to the odds ratio
 * via \f$\OStSGLtL = \oStSGLtL\,\BStSGLtL\f$ in terms of the prior odds \f$\oStSGLtL\f$.
 *
 * This function returns an approximation of \f$\logten B_\StSGLtL = \ln B_\StSGLtL \, \logten e\f$,
 * with the general expression for \f$\ln B_\StSGLtL\f$ given in Eq.(35) of \cite Keitel2015 :
 * \f{equation}{
 * \ln \BStSGLtL =   \numermax + \ln \left( \sum\limits_{\numerterm \in \numertermset} e^{\numerterm-\numermax} \right)
 *                 - \denommax - \ln \left( \sum\limits_{\denomterm \in \denomtermset} e^{\denomterm-\denommax} \right) \,,
 * \f}
 * with the denominator terms defined as in XLALComputeGLtLDenominator() and
 * \f{equation}{
 * \numertermset = \left\{ \scF + \ln\cppS , \;
                           \cohF^{\ell} + (\Nseg-1)\Ftho + \ln\cppTSk
                   \right\} \,,
 * \f}
 * while the current implementation uses only the maximum multi-detector, single-segment F-stat value \f$\max\cohF^{\ell}\f$ for the numerator:
 * \f$\numermax = \max \left\{  \scF, \max\limits_{\ell}\cohF^\ell + \Ftho (\Nseg-1) \right\}\f$
 * and, through XLALComputeGLtLDenominator(), the maximum \f$\denommax\f$ of the set of denominator exponents.
 * dropping the additional log-correction terms.
 *
 * NOTE: This implementation also assumes unit prior weights for the various signal and noise hypotheses.
 *
 */
REAL4
XLALComputeBStSGLtL ( const REAL4 twoF,					//!< [in] semi-coherent sum \f$2\scF\f$ of multi-detector F-stats
                      const REAL4 maxtwoFl,				//!< [in] maximum \f$\max\limits_{\ell}2\cohF^\ell\f$ of multi-detector F-stats over segments
                      const REAL4 twoFX[PULSAR_MAX_DETECTORS],		//!< [in] semi-coherent sums \f$\{2\scF^X\}\f$ of per-detector F-stats
                      const REAL4 maxtwoFXl[PULSAR_MAX_DETECTORS],	//!< [in] maxima \f$\{\max\limits_{\ell}2\{\cohF^{X\ell}\}\}\f$ of per-detector F-stats over segments
                      const BSGLSetup *setup				//!< [in] pre-computed setup from XLALCreateBSGLSetup()
                      )
{
  XLAL_CHECK ( setup != NULL, XLAL_EINVAL );

  REAL4 GLtLDenominator = XLALComputeGLtLDenominator ( twoFX, maxtwoFXl, setup );

  // assuming equal pS=ptSl priors for all segments
  REAL4 multiF = 0.5 * twoF;
  REAL4 tsTerm = 0.5 * maxtwoFl + setup->perSegTerm;
  REAL4 maxNumerator = fmax ( multiF, tsTerm );
  REAL4 minNumerator = fmin ( multiF, tsTerm );

  REAL4 ln_BStSGLtL  = maxNumerator + log ( 1 + exp ( minNumerator - maxNumerator ) ) - GLtLDenominator;

  return ln_BStSGLtL * LAL_LOG10E; // return log10(B_SGL)

} // XLALComputeBStSGLtL

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
      XLAL_CHECK ( XLALParseStringValueAsREAL4 ( &oLGX[X], oLGX_string->data[X] ) == XLAL_SUCCESS, XLAL_EFUNC );
    } // for X < numDetectors

  return XLAL_SUCCESS;

} // XLALParseLinePriors()
