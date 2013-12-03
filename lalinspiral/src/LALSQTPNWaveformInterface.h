#ifndef LALSQTPNWAVEFORMINTERFACE_H
#define LALSQTPNWAVEFORMINTERFACE_H

#include <lal/Units.h>
#include <lal/LALInspiral.h>
#include <lal/SeqFactories.h>
#include <lal/LALStatusMacros.h>
#include <lal/GenerateInspiral.h>
#include <lal/LALSQTPNWaveform.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \defgroup LALSQTPNWaveformInterface_h Header LALSQTPNWaveformInterface.h
 * \ingroup pkg_CBC_NEW
 *
 * \brief Contains function declarations to integrate the SpinQuadTaylor code into the other parts of the LALSuit.
 * \author László Veréb
 * \date 2010.06.27.
 */
/*@{*/

#define LALSQTPN_MSGPPNPARAMS "the PPNParamsStruct structure is null"
#define LALSQTPN_MSGINSPIRALTEMPLATE "the InspiralTemplate structure is null"
#define LALSQTPN_ZEROLENGTH 3
#define LALSQTPN_MSGZEROLENGTH "the given length is not positive"

int XLALSQTPNWaveformTemplates (REAL4Vector *signalvec1,
		REAL4Vector *signalvec2, InspiralTemplate *params);

// LAL wrapper to XLAL function
void LALSQTPNWaveformTemplates (LALStatus *status, REAL4Vector *signalvec1,
		REAL4Vector *signalvec2, InspiralTemplate *params);

/**
 * The function returns the generated waveform.
 * @param[out]	signalvec	: array containing the waveform \f$(h_+, h_\times)\f$
 * @param[in]	params		: structure containing the inspiral parameters
 */
int XLALSQTPNWaveform (REAL4Vector *signalvec,
		InspiralTemplate *params);

// LAL wrapper to XLAL function
void LALSQTPNWaveform (LALStatus *status, REAL4Vector *signalvec,
		InspiralTemplate *params);

/**
 * The function returns the generated waveform for injection.
 * @param[out]	wave_out	: structure containing the waveform \f$(a_1, a_2, \Phi, \alpha)\f$
 * @param[in]	params		: structure containing the inspiral parameters
 * @param[in]	ppnParams	: parameters for restricted post-Newtonian waveform
 */
int XLALSQTPNWaveformForInjection(CoherentGW *wave_out,
		InspiralTemplate *params, PPNParamStruc *ppnParams);

// LAL wrapper to XLAL function
void LALSQTPNWaveformForInjection(LALStatus *status, CoherentGW *wave_out,
		InspiralTemplate *params, PPNParamStruc *ppnParams);

/**
 * The function allocates memory for the waveform's \f$a_1\f$, \f$a_2\f$,
 * \f$\Phi\f$ and \f$\alpha\f$.
 * @param[out]		wave	: pointer to the allocated waveform
 * @param[in]		length	: the length of the waveform
 */
int XLALSQTPNAllocateCoherentGW(CoherentGW *wave, UINT4 length);

/**
 * The function deallocates memory of the waveform.
 * @param[out]		wave	: pointer to the allocated waveform
 */
void XLALSQTPNDestroyCoherentGW(CoherentGW *wave);

/**
 * The function calculates the parameters from the InspiralTemplate
 * structure. <em>The used parameters are:</em>
 * <ul>
 * <li>masses of the BHs (or NSs) \f$m_i\f$ in \f$M_\odot\f$</li>
 * <li>the spin components \f$\chi_{ij}\f$, the values of \f$\sqrt{\sum_j\chi_{ij}}\f$, are between 0 and 1</li>
 * <li>the quadrupole parameters \f$w_i\in(4,8)\f$ for NSs [1] and \f$w_i=1\f$ for BHs[2] are 1 (default 1)</li>
 * <li>the inclination (angle between the line of sight and Newtonian orbital angular momentum) \f$\iota\f$ in \f$rad\f$
 * <li>the initial frequency \f$f_L\f$ in \f$Hz\f$</li>
 * <li>the distance \f$d\f$ in \f$Mpc\f$</li>
 * <li>the sampling time \f$t_s\f$ in \f$s\f$</li>
 * <li>the PN order, see #LALPNOrder</li>
 * <li>level of accuracy in including spin and quadrupole contributions, see LALSQTPNSpinInteraction </li>
 * </ul><br />
 * <em>The calculated parameters:</em>
 * \f{gather}{
 * \displaystyle M_{in}=m_1+m_2,\quad
 * \mu=\frac{m_1m_2}{M_{in}},\quad
 * \eta=\frac{\mu}{M_{in}},\\
 * \chi_i=\sqrt{\sum_{j}\chi_{ij}^2},\quad
 * \hat{\chi}_{ij}=\dfrac{\chi_{ij}}{\chi_i},\\
 * f_s=t_s^{-1}\\
 * A=\frac{4\cdot\eta M_{in}M_\odot\displaystyle\frac{G}{c^2}}{d\cdot3.0856775807\cdot10^{16}\cdot10^6}
 * \f}
 * and the initial phase \f$\phi=0\f$
 * Assuming that:
 * <ul>
 * <li>masses are positive</li>
 * <li>eta is positive</li>
 * <li>sampling frequency is positive</li>
 * <li>distance is positive</li>
 * </ul>
 * @param[out]	wave		: the used parameters
 * @param[in]		params	: the inspiral parameters
 */
void XLALSQTPNFillParams(LALSQTPNWaveformParams *wave, InspiralTemplate *params);

/*@}*/ /* end:LALSQTPNWaveformInterface_h */

#ifdef __cplusplus
}
#endif

#endif /* LALSQTPNWAVEFOMRINTERFACE_H */
