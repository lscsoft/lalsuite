/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Warren Anderson
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

#ifndef _RING_H
#define _RING_H

#include <lal/LALDatatypes.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif


/**
 * \defgroup RingUtils_h Header RingUtils_h
 * \ingroup pkg_ring
 * \author Jolien Creighton
 *
 * \brief Black hole ringdown waveform generation.
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/RingUtils.h>
 * \endcode
 *
 * Routines for generating waveforms for black hole ringdown.
 *
 * The ringdown waveform is an exponentially-damped sinusoid
 * \f{equation}{
 * r(t) = \left\{
 * \begin{array}{ll}
 * e^{-\pi ft/Q}\cos(2\pi ft + \phi_0) & \mbox{for } t\ge0 \\
 * 0 & \mbox{for } t<0
 * \end{array}
 * \right.
 * \f}
 * where \f$f\f$ is the central frequency of the ringdown waveform, \f$Q\f$ is
 * the quality factor, and \f$\phi_0\f$ is the initial phase of the waveform.
 * Note that Ref. [\ref JDECreighton99] adopted the
 * normalization convention \f$q(t)=(2\pi)^{1/2}r(t)\f$.
 *
 * For a black hole ringdown, the gravitational waveform produced, averaged
 * over the various angles, is
 * \f{equation}{
 * h(t) = A_qq(t)
 * \f}
 * where the central frequency and quality of the ringdown are determined from
 * the mass and spin of the black holes.  An analytic approximation
 * yields [\ref EWLeaver85,\ref FEcheverria89]
 * \f{equation}{
 * f \simeq 32\,\textrm{kHz}\times[1-0.63(1-{\hat{a}})^{3/10}](M_\odot/M)
 * \f}
 * and
 * \f{equation}{
 * Q \simeq 2(1-{\hat{a}})^{-9/20}
 * \f}
 * with the black hole mass given by \f$M\f$ and its spin by \f$S={\hat{a}}GM^2/c\f$
 * (where \f$G\f$ is Newton's constant and \f$c\f$ is the speed of light).  The
 * dimensionless spin parameter \f${\hat{a}}\f$ lies between zero (for a
 * Schwarzschild black hole) and unity (for an extreme Kerr black hole).
 * The amplitude of the waveform depends on these quantities as well as the
 * distance \f$r\f$ to the source and the fractional mass loss \f$\epsilon\f$ radiated
 * in gravitational waves [\ref JDECreighton99]:
 * \f{equation}{
 * A_q = 2.415\times10^{-21}Q^{-1/2}[1-0.63(1-{\hat{a}})^{3/10}]^{-1/2}
 * \left(\frac{\textrm{Mpc}}{r}\right)
 * \left(\frac{M}{M_\odot}\right)
 * \left(\frac{\epsilon}{0.01}\right)^{1/2}.
 * \f}
 * Note that this is the amplitude factor for the waveform \f$q(t)\f$, whereas
 * the amplitude factor for \f$r(t)\f$ would be \f$(2\pi)^{1/2}A_q\f$.
 *
 * The mismatch between two nearby templates is given by \f$ds^2\f$, which can be
 * thought of as the line interval for a mismatch-based metric on the \f$(f,Q)\f$
 * parameter space [\ref Owen_96,\ref JDECreighton99]:
 * \f{equation}{
 * ds^2 = \frac{1}{8} \biggl\{ \frac{3+16Q^4}{Q^2(1+4Q^2)^2}\,dQ^2
 * - 2\frac{3+4Q^2}{fQ(1+4Q^2)}\,dQ\,df + \frac{3+8Q^2}{f^2}\,df^2 \biggr\}.
 * \f}
 * When expressed in terms of \f$\log f\f$ rather than \f$f\f$, the metric coefficients
 * depend on \f$Q\f$ alone.  We can exploit this property for the task of template
 * placement.  The method is the following:  First, choose a "surface" of
 * constant \f$Q=Q_{\mathrm{\scriptstyle min}}\f$, and on this surface place
 * templates at intervals in \f$\phi=\log f\f$ of \f$d\phi=d\ell/\surd g_{\phi\phi}\f$
 * for the entire range of \f$\phi\f$.  Here,
 * \f$d\ell=\surd(2ds^2_{\mathrm{\scriptstyle threshold}})\f$.  Then choose the
 * next surface of constant \f$Q\f$ with \f$dQ=d\ell/\surd g_{QQ}\f$ and repeat the
 * placement of templates on this surface.  This can be iterated until the
 * entire range of \f$Q\f$ has been covered; the collection of templates should now
 * cover the entire parameter region with no point in the region being farther
 * than \f$ds^2_{\mathrm{\scriptstyle threshold}}\f$ from the nearest template.
 *
 * \heading{Algorithm}
 *
 * The waveform generation routines use recurrance relations for both the
 * exponentially-decaying envelope and for the co-sinusoid.
 *
 * The template placement algorithm is described above.
 *
 */
/*@{*/

/**\name Error Codes */
/*@{*/
#define RINGH_ENULL 01	/**< Null pointer */
#define RINGH_ENNUL 02	/**< Non-null pointer */
#define RINGH_EALOC 04	/**< Memory allocation error */
/*@}*/

/** \cond DONT_DOXYGEN */
#define RINGH_MSGENULL "Null pointer"
#define RINGH_MSGENNUL "Non-null pointer"
#define RINGH_MSGEALOC "Memory allocation error"
/** \endcond */

/**
 * This structure contains a bank of ringdown waveforms.
 */
typedef struct
tagRingTemplateBank
{
  UINT4              numTmplt;	/**< The number of templates in the bank */
  SnglRingdownTable *tmplt;	/**< Array of ringdown templates */
}
RingTemplateBank;

/**
 * This structure contains the parameters required for generating a ringdown
 * template bank.
 */
typedef struct
tagRingTemplateBankInput
{
  REAL4 minQuality;		/**< The minimum quality factor in the bank */
  REAL4 maxQuality;		/**< The maximum quality factor in the bank */
  REAL4 minFrequency;		/**< The minimum central frequency in the bank (in Hz) */
  REAL4 maxFrequency;		/**< The minimum central frequency in the bank (in Hz) */
  REAL4 maxMismatch;		/**< The maximum mismatch allowed between templates in the bank */
  REAL4 templatePhase;		/**< The phase of the ringdown templates, in radians;
                                 * Zero is a cosine-phase template;
                                 * \f$-\pi/2\f$ is a sine-phase template.
                                 */
  REAL4 templateDistance;	/**< UNDOCUMENTED */
  REAL4 templateEpsilon;	/**< UNDOCUMENTED */
}
RingTemplateBankInput;

/* ---------- Function prototypes ---------- */

REAL4 XLALBlackHoleRingSpin( REAL4 Q );
REAL4 XLALBlackHoleRingMass( REAL4 f, REAL4 Q );
REAL4 XLALBlackHoleRingQuality( REAL4 a );
REAL4 XLALBlackHoleRingFrequency( REAL4 M, REAL4 a );
REAL4 XLALNonSpinBinaryFinalBHSpin( REAL4 eta );
REAL4 XLALNonSpinBinaryFinalBHMass( REAL4 eta, REAL4 mass1, REAL4 mass2 );
REAL4 XLALSpinBinaryFinalBHSpin( REAL4 eta, REAL4 mass1, REAL4 mass2, REAL4 spin1x, REAL4 spin2x,
   REAL4 spin1y, REAL4 spin2y, REAL4 spin1z, REAL4 spin2z );
REAL4 XLALBlackHoleRingAmplitude( REAL4 f, REAL4 Q, REAL4 r, REAL4 epsilon );
REAL4 XLALBlackHoleRingEpsilon( REAL4 f, REAL4 Q, REAL4 r, REAL4 amplitude );
REAL4 XLALBlackHoleRingHRSS( REAL4 f, REAL4 Q, REAL4 amplitude, REAL4 plus, REAL4 cross );
REAL8 XLAL2DRingMetricDistance( REAL8 fa, REAL8 fb, REAL8 Qa, REAL8 Qb );
REAL8 XLAL3DRingMetricDistance( REAL8 fa, REAL8 fb, REAL8 Qa, REAL8 Qb, REAL8 dt );
REAL8 XLAL3DRingTimeMinimum( REAL8 fa, REAL8 fb, REAL8 Qa, REAL8 Qb);
REAL8 XLALRingdownTimeError( const SnglRingdownTable *table,  REAL8 lal_ring_ds_sq );

int XLALComputeRingTemplate(
    REAL4TimeSeries *output, SnglRingdownTable *input );
int XLALComputeBlackHoleRing(
    REAL4TimeSeries *output, SnglRingdownTable *input, REAL4 dynRange );
RingTemplateBank *XLALCreateRingTemplateBank( RingTemplateBankInput *input );
void XLALDestroyRingTemplateBank( RingTemplateBank *bank );


/*@}*/

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _RING_H */
