/**** <lalVerbatim file="RingHV">
 * Author: Jolien Creighton
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 * 
 * \section{Header \texttt{Ring.h}}
 * \label{sec:Ring.h}
 *
 * Black hole ringdown waveform generation.
 *
 * \subsection*{Synopsis}
 * \begin{verbatim}
 * #include <lal/Ring.h>
 * \end{verbatim}
 * 
 * Routines for generating waveforms for black hole ringdown.
 *
 * The ringdown waveform is an exponentially-damped sinusoid
 * \begin{equation}
 *   r(t) = \left\{
 *   \begin{array}{ll}
 *     e^{-\pi ft/Q}\cos(2\pi ft + \phi_0) & \mbox{for $t\ge0$} \\
 *     0 & \mbox{for $t<0$}
 *   \end{array}
 *   \right.
 * \end{equation}
 * where $f$ is the central frequency of the ringdown waveform, $Q$ is
 * the quality factor, and $\phi_0$ is the initial phase of the waveform.
 * Note that Ref.~\cite{JDECreighton} adopted the
 * normalization convention $q(t)=(2\pi)^{1/2}r(t)$.
 *
 * For a black hole ringdown, the gravitational waveform produced, averaged
 * over the various angles, is
 * \begin{equation}
 *   h(t) = A_qq(t)
 * \end{equation}
 * where the central frequency and quality of the ringdown are determined from
 * the mass and spin of the black holes.  An analytic approximation
 * yields~\cite{EWLeaver,FEcheverria}
 * \begin{equation}
 *   f \simeq 32\,\textrm{kHz}\times[1-0.63(1-{\hat{a}})^{3/10}](M_\odot/M)
 * \end{equation}
 * and
 * \begin{equation}
 *   Q \simeq 2(1-{\hat{a}})^{-9/20}
 * \end{equation}
 * with the black hole mass given by $M$ and its spin by $S={\hat{a}}GM^2/c$
 * (where $G$ is Newton's constant and $c$ is the speed of light).  The
 * dimensionless spin parameter ${\hat{a}}$ lies between zero (for a
 * Schwarzschild black hole) and unity (for an extreme Kerr black hole).
 * The amplitude of the waveform depends on these quantities as well as the
 * distance $r$ to the source and the fractional mass loss $\epsilon$ radiated
 * in gravitational waves~\cite{JDECreighton}:
 * \begin{equation}
 *   A_q = 2.415\times10^{-21}Q^{-1/2}[1-0.63(1-{\hat{a}})^{3/10}]^{-1/2}
 *   \left(\frac{\textrm{Mpc}}{r}\right)
 *   \left(\frac{M}{M_\odot}\right)
 *   \left(\frac{\epsilon}{0.01}\right)^{1/2}.
 * \end{equation}
 * Note that this is the amplitude factor for the waveform $q(t)$, whereas
 * the amplitude factor for $r(t)$ would be $(2\pi)^{1/2}A_q$.
 *
 * The mismatch between two nearby templates is given by $ds^2$, which can be
 * thought of as the line interval for a mismatch-based metric on the $(f,Q)$
 * parameter space~\cite{BJOwen,JDECreighton}:
 * \begin{equation}
 *   ds^2 = \frac{1}{8} \biggl\{ \frac{3+16Q^4}{Q^2(1+4Q^2)^2}\,dQ^2
 *   - 2\frac{3+4Q^2}{fQ(1+4Q^2)}\,dQ\,df + \frac{3+8Q^2}{f^2}\,df^2 \biggr\}.
 * \end{equation}
 * When expressed in terms of $\log f$ rather than $f$, the metric coefficients
 * depend on $Q$ alone.  We can exploit this property for the task of template
 * placement.  The method is the following:  First, choose a ``surface'' of
 * constant~$Q=Q_{\mathrm{\scriptstyle min}}$, and on this surface place
 * templates at intervals in $\phi=\log f$ of~$d\phi=d\ell/\surd g_{\phi\phi}$
 * for the entire range of~$\phi$.  Here,
 * $d\ell=\surd(2ds^2_{\mathrm{\scriptstyle threshold}})$.  Then choose the
 * next surface of constant $Q$ with~$dQ=d\ell/\surd g_{QQ}$ and repeat the
 * placement of templates on this surface.  This can be iterated until the
 * entire range of~$Q$ has been covered; the collection of templates should now
 * cover the entire parameter region with no point in the region being farther
 * than~$ds^2_{\mathrm{\scriptstyle threshold}}$ from the nearest template.
 * 
 **** </lalLaTeX> */

#ifndef _RING_H
#define _RING_H

#include <lal/LALDatatypes.h>

NRCSID( RINGH, "$Id$" );

#ifdef  __cplusplus
extern "C" {
#pragma } /** to match the previous brace **/
#endif

/**** <lalLaTeX>
 * \subsection*{Error conditions}
 **** </lalLaTeX> */
/**** <lalErrTable> */
#define RINGH_ENULL 01
#define RINGH_ENNUL 02
#define RINGH_EALOC 04
#define RINGH_MSGENULL "Null pointer"
#define RINGH_MSGENNUL "Non-null pointer"
#define RINGH_MSGEALOC "Memory allocation error"
/**** </lalErrTable> */

/**** <lalLaTeX>
 * 
 * \subsection*{Structures}
 * 
 * \subsubsection*{Type \texttt{RingTemplateInput}}
 * \idx[Type]{RingTemplateInput}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagRingTemplateInput
{
  REAL4 quality;
  REAL4 frequency;
  REAL4 phase;
}
RingTemplateInput;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 *
 * This structure contains the required information for generating a ringdown
 * template $q(t)$.  The fields are:
 * \begin{description}
 * \item[\texttt{quality}] The quality factor $Q$ of the ringdown waveform.
 * \item[\texttt{frequency}] The central frequency of the ringdown waveform
 *     (in Hz).
 * \item[\texttt{phase}] The initial phase of the ringdown in radians.
 *     Zero is a cosine-phase ringdown; $-\pi/2$ is a sine-phase ringdown.
 * \end{description}
 *
 *
 * \subsubsection*{Type \texttt{BlackHoleRingInput}}
 * \idx[Type]{BlackHoleRingInput}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagBlackHoleRingInput
{
  REAL4 solarMasses;
  REAL4 dimensionlessSpin;
  REAL4 percentMassLoss;
  REAL4 distanceMpc;
  REAL4 initialPhase;
}
BlackHoleRingInput;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 *
 * This structure contains the physical parameters for generating the ringdown
 * waveform from a black hole source.  The fields are:
 * \begin{description}
 * \item[\texttt{solarMasses}] The mass $M$ of the black hole (in solar
 *     masses, $M_\odot$).
 * \item[\texttt{dimensionlessSpin}] The dimensionless spin parameter of the
 *     black hole ${\hat{a}}$ where the spin is $S={\hat{a}}GM^2/c$ ($G$ is
 *     Newton's constant and $c$ is the speed of light).
 * \item[\texttt{percentMassLoss}] The fractional mass loss, as a percent of
 *     the initial black hole mass, in ringdown radiation.
 * \item[\texttt{distanceMpc}] The distance of the source in megaparsecs (Mpc).
 * \item[\texttt{initialPhase}] The initial phase of the ringdown in radians.
 *     Zero is a cosine-phase ringdown; $-\pi/2$ is a sine-phase ringdown.
 * \end{description}
 *
 *
 * \subsubsection*{Type \texttt{RingTemplateBank}}
 * \idx[Type]{RingTemplateBank}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagRingTemplateBank
{
  UINT4              numTmplt;
  RingTemplateInput *tmplt;
}
RingTemplateBank;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 *
 * This structure contains a bank of ringdown waveforms.  The fields are:
 * \begin{description}
 * \item[\texttt{numTmplt}] The number of templates in the bank.
 * \item[\texttt{tmplt}] Array of ringdown templates.
 * \end{description}
 *
 * \subsubsection*{Type \texttt{RingTemplateBankInput}}
 * \idx[Type]{RingTemplateBankInput}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagRingTemplateBankInput
{
  REAL4 minQuality;
  REAL4 maxQuality;
  REAL4 minFrequency;
  REAL4 maxFrequency;
  REAL4 maxMismatch;
  REAL4 templatePhase;
}
RingTemplateBankInput;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 *
 * This structure contains the parameters required for generating a ringdown
 * template bank.  The fields are:
 * \begin{description}
 * \item[\texttt{minQuality}] The minimum quality factor in the bank.
 * \item[\texttt{maxQuality}] The maximum quality factor in the bank.
 * \item[\texttt{minFrequency}] The minimum central frequency in the bank
 *     (in Hz).
 * \item[\texttt{maxFrequency}] The minimum central frequency in the bank
 *     (in Hz).
 * \item[\texttt{maxMismatch}] The maximum mismatch allowed between templates 
 *     in the bank.
 * \item[\texttt{templatePhase}] The phase of the ringdown templates, in
 *     radians.  Zero is a cosine-phase template; $-\pi/2$ is a sine-phase
 *     template.
 * \end{description}
 *
 **** </lalLaTeX> */

int XLALComputeRingTemplate( REAL4TimeSeries *output, RingTemplateInput *input );
int XLALComputeBlackHoleRing( REAL4TimeSeries *output, BlackHoleRingInput *input );
RingTemplateBank *XLALCreateRingTemplateBank( RingTemplateBankInput *input );
void XLALDestroyRingTemplateBank( RingTemplateBank *bank );

void
LALComputeRingTemplate(
    LALStatus          *status,
    REAL4TimeSeries    *output,
    RingTemplateInput  *input
    );

void
LALComputeBlackHoleRing(
    LALStatus          *status,
    REAL4TimeSeries    *output,
    BlackHoleRingInput *input
    );

void
LALCreateRingTemplateBank(
    LALStatus              *status,
    RingTemplateBank      **output,
    RingTemplateBankInput  *input
    );

void
LALDestroyRingTemplateBank(
    LALStatus         *status,
    RingTemplateBank **bank
    );

/**** <lalLaTeX>
 * 
 * \vfill{\footnotesize\input{RingHV}}
 * \newpage\input{RingC}
 * \newpage\input{RingTestC}
 * 
 **** </lalLaTeX> */

#ifdef  __cplusplus
#pragma { /** to match the next brace **/
}
#endif

#endif /* _RING_H */
