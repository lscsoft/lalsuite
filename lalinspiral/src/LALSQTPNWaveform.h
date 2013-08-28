#ifndef LALSQTPNWAVEFORM_H
#define LALSQTPNWAVEFORM_H

#include <math.h>

#include <lal/LALStatusMacros.h>
#include <lal/LALInspiral.h>
#include <lal/Units.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/**
 * \defgroup LALSQTPNWaveform_h Header LALSQTPNWaveform.h
 * \ingroup pkg_CBC_NEW
 * \brief Contains the enumerations, structures and functions declarations to create GWforms.
 *
 * \heading{Notations}
 *
 * \f$M=M_{in}M_\odot G/c^3\f$<br />
 * \f$\hat{L_N}=\left(\sin\iota,0,\cos\iota\right)\f$ is the direction of the
 * orbital angular momentum in radiation frame see [1].<br />
 * \f$\omega=\pi f_L\f$<br />
 * \f$\hat{\chi}_i\f$, \f$M_{in}\f$, \f$\iota\f$, \f$f_L\f$ are defined in
 * LALSQTPNFillParams() function.<br />
 * <b>References</b><br />
 * [1] L. E. Kidder, Phys.Rev. D52, 821 (1995)<br />
 * [2] Alessandra Buonanno, Yanbei Chen, and Michele Vallisneri, Phys.Rev. D67 (2003) 104025; Erratum-ibid. D74 (2006) 029904<br />
 * [3] Balázs Mikóczi, Mátyás Vasúth, László Á. Gergely, Phys.Rev. D71 (2005) 124043
 * \author László Veréb
 * \date 2010.05.21.
 */
/*@{*/

/**
 * The macro function returns the square of the argument.
 * Do not use with incrementing or decrementing operators!
 * @param[in] a	 : the number
 * @return the square of the number
 */
#define SQT_SQR(a) ((a)*(a))

typedef struct tagLALSQTPNWave{
	CoherentGW *waveform;
	REAL4Vector *h;
	REAL4Vector *hp;
	REAL4Vector *hc;
	UINT4 length;
} LALSQTPNWave;

/**
 * The structure contains the coefficients for calculating the derivatives
 * of the evolving quantities.
 */
typedef struct tagLALSQTPNCoefficients{
	///@name coefficients for domega
	//@{
	REAL8 domegaGlobal; ///< global coefficient for domega
	REAL8 domega[LAL_PNORDER_PSEUDO_FOUR]; ///< coefficients for domega for every PN order
	REAL8 domegaSO[3]; ///< the spin-orbit coefficients for domega
	REAL8 domegaSS[3]; ///< the spin1-spin2 coefficients for domega
	REAL8 domegaSSself[3]; ///< the spin-selft coefficients for domega
	REAL8 domegaSSselfConst; ///< the constant spin-selft coefficients for domega
	REAL8 domegaQM[3]; ///< the quadropole-monopole coefficients for domega
	REAL8 domegaQMConst; ///< the constant quadropole-monopole coefficients for domega
	REAL8 domegaLN; ///< coefficient for the ln component in domega
	//@}
	///@name coefficients for dchih and dLNh
	REAL8 dchihSO[3]; ///< the spin-orbit coefficients for dchih
	REAL8 dchihSS[3]; ///< the spin1-spin2 coefficientd for dchih
	REAL8 dchihQM[3]; ///< the quadropole-monopole coefficients for dchih
	REAL8 dLNh[3]; ///< coefficients for dLNh
	//@}
	///@name coefficients for MECO
	//@{
	REAL8 meco[8]; ///< coefficients for MECO-test
	REAL8 mecoSO[3]; ///< spin-orbit coefficients for MECO
	REAL8 mecoSS; ///< spin1-spin2 coefficients for MECO
	REAL8 mecoQM; ///< quadropole-monopole coefficients for MECO
	//@}
} LALSQTPNCoefficients;

/**
 * The structure contains the system's and the generator's parameters.
 */
typedef struct tagLALSQTPNWaveformParams{
	///@name mass-Parameters
	//@{
	REAL8 mass[2]; ///< masses of the BHs in \f$M_\odot\f$
	REAL8 totalMass; ///< total mass in \f$M_\odot\f$
	REAL8 chirpMass; ///< chirp mass in \f$M_\odot\f$
	REAL8 mu; ///< reduced mass in \f$M_\odot\f$
	REAL8 eta; ///< symmetric mass ratio
	//@}
	///@name spin-parameters
	//@{
	REAL8 chi[2][3]; ///< components of the normalized spin
	REAL8 chih[2][3]; ///< components of the unity-vectors of the normalized spin
	REAL8 chiAmp[2]; ///< amplitude of the normalized spin
	//@}
	///@name other system-parameters
	//@{
	REAL8 qmParameter[2]; ///< flatness of the BHs or NSs
	REAL8 distance; ///< distance to the source in \f$Mps\f$
	REAL8 inclination; ///< inclination of the system \f$rad\f$
	REAL8 phi; ///< the initial phase (currently not in use)
	//@}
	///@name other parameters
	//@{
	REAL8 signalAmp; ///< the amplitude of the signal
	REAL8 lowerFreq; ///< the detectors sensitivityband's lower border in \f$Hz\f$
	REAL8 finalFreq;	///< the final frequency
	REAL8 samplingFreq; ///< sampling frequency in \f$Hz\f$
	REAL8 samplingTime; ///< sampling time in \f$s\f$
	REAL8 coalescenceTime;	///< the time at the coalescence
	LALPNOrder order; ///< the Post_Newtonian order of the GW generation
	LALSimInspiralInteraction interaction; ///< which spin interaction will be included in the generation
	LALSQTPNCoefficients coeff; ///< coefficients for the deriving the parameters
	//@}
} LALSQTPNWaveformParams;

/**
 * The function fills the #LALSQTPNCoefficients structure with the needed
 * coefficients for calculating the derived dynamic variables with the LALSQTPNDerivator() function.
 *
 * The orders above 2PN are incomplete, so use them if you want to try their
 * effects.
 * @param[in,out]	params	: the LALSQTPN_Generator's parameters
 */
int XLALSQTPNFillCoefficients(LALSQTPNWaveformParams * const params);

/**
 * The function calculates the derived values.
 * The formulae are:
 * \f{center}{
 * \newcommand{\OM}[1]{\left(M\omega\right)^{#1/3}}
 * \newcommand{\derTpM}[1]{\dfrac{d#1}{d\left(t/M\right)}}
 * \newcommand{\SPU}[2]{\left(\hat{#1}\hat{#2}\right)}
 * \newcommand{\VPU}[2]{\left(\hat{#1}\times\hat{#2}\right)}
 * \newcommand{\SP}[2]{\mathbf{\text{UNDEFINDE}[#1][#2]}}
 * \begin{gather}
 * \begin{gathered}\tag{eq:LALSQTPNchih}
 * \derTpM{{\hat{\chi}}_i}={SO}_{\chi i}\OM{5}+\left({SS}_{\chi i}+{QM}_{\chi i}\right)\OM{6},\\
 * {SO}_{\chi i}=\frac{\eta}{2}\left(4+3\frac{m_j}{m_i}\right)\left(\hat{L_N}\times\hat{\chi_i}\right),\\
 * {SS}_{\chi i}=\frac{1}{2}\frac{\chi_jm_j^2}{M^2}\left[\hat{\chi_j}-3\left(\hat{L_N}\hat{\chi_j}\right)\hat{L_N}\right]\times\hat{\chi_i},\\
 * {QM}_{\chi i}=-\frac{3}{2}\eta\chi_iw_i\left(\hat{L_N}\hat{\chi_i}\right)\left(\hat{L_N}\times\hat{\chi_i}\right),
 * \end{gathered}\\[15pt]
 * \derTpM{\hat{L_N}}=\sum_i-\frac{1}{\eta}\frac{\chi_im_i^2}{M^2}\derTpM{\hat{\chi_i}},\\[15pt]
 * \begin{gathered}\tag{eq:LALSQTPNomega}
 * \begin{split}
 * \derTpM{\left(M\omega\right)}&=\frac{96\eta}{5}\OM{11}\bigg[1-\frac{743+924\eta}{336}\OM{2}\bigg.\\&
 * \quad+\left(4\pi+SO_{\omega}\right)\OM{3}+\bigg(\frac{34103}{18144}+\frac{13661}{2016}\eta+\frac{59}{18}\eta^2\bigg.\\&
 * \quad+\bigg.\bigg.SS_{\omega}+SSself_{\omega}+QM_{\omega}\bigg)\OM{4}\bigg],
 * \end{split}\\
 * SO_{\omega}=\sum_{i\ne j}-\frac{1}{12}\frac{\chi_im_i^2}{M^2}\left(113+75\frac{m_j}{m_i}\right)\SPU{L_N}{\chi_i},\\
 * SS_{\omega}=\frac{\eta\chi_1\chi_2}{48}\left[721\SPU{L_N}{\chi_1}\SPU{L_N}{\chi_2}-247\SPU{\chi_1}{\chi_2}\right],\\
 * SSself_{\omega}=\sum_{i}\frac{1}{96}\frac{\chi_im_i}{M^2}\chi_i\left[7-\SPU{L_N}{\chi_i}^2\right],\\
 * QM_{\omega}=\sum_{i}\frac{5}{2}\frac{\chi_im_i^2}{M^2}\chi_iw_i\left[3\SPU{L_N}{\chi_i}^2-1\right],
 * \end{gathered}\\[15pt]
 * \begin{gathered}
 * \begin{split}
 * MECO&=-0.5\eta\frac{2}{3}\OM{-1}+0.5\eta\frac{4}{3}\frac{9+\eta}{12}\OM{1}\\&\quad+
 * SO_{MECO}\OM{2}+\Big(0.5\eta\frac{6}{3}\frac{81-57\eta+\eta^2}{24}\Big.\\&\quad+
 * \Big.SS_{MECO}+QM_{MECO}\Big)\OM{3},
 * \end{split}\\
 * SO_{MECO}=\sum_{i}-\frac{5}{9}\eta\frac{\chi_im_i^2}{M^2}\left(4+3\frac{m_j}{m_i}\right)\SP{\hat{L}_N}{\hat{\chi}_i};\quad
 * QM_{MECO}=2\eta QM_{\omega}\\
 * SS_{MECO}=-\displaystyle\frac{\chi_1m_1^2}{M^2}\frac{\chi_2m_2^2}{M^2}\left[\SP{\hat{\chi}_1}{\hat{\chi}_2}-3\SP{\hat{L}_N}{\hat{\chi}_1}\SP{\hat{L}_N}{\hat{\chi}_2}\right]
 * \end{gathered}\\[15pt]
 * \begin{gathered}
 * \derTpM{\Phi}=M\omega-\derTpM{\alpha}\cos\iota\\\alpha=\arctan\dfrac{\hat{L_N}_y}{\hat{L_N}_x}
 * \end{gathered}
 * \end{gather}
 * \f}
 * Equation (1) is given by (Eq. (2)-(3)) of [2] except of the quadropole-monopole contribution
 * (\f$QM_{\chi_i}\f$) that was added here, equation (2) is given by (Eq. (9)) of [2],
 * equation (3) is given by (Eq. (7)-(9)) of [3],
 * equation (4) is given by (Eq. (11),(13)) of [2] except of the quadropole-monopole
 * contribution (\f$QM_{MECO}\f$) that was added here, equation
 * (5) is given by (Eq. (4.31)) of [1].<br />
 */
/**
 * The equation's constant parts are calculated by the LALSQTPNFillCoefficients() function.
 * This function is used by the LALSQTPNIntegratorFunc() function to evolve the system with a given time.
 * @param[in]	t		: evolution time, not in used
 * @param[in]	values	: the values to be derived
 * @param[out]	dvalues	: the derived values and the last element is the MECO
 * @param[in]	params	: the LALSQTPN_Generator's parameters
 */
int XLALSQTPNDerivator(UNUSED REAL8 t, const REAL8 values[], REAL8 dvalues[],
		void * params);

// LAL wrapper of above XLAL function
int LALSQTPNDerivator(REAL8 t, const REAL8 values[], REAL8 dvalues[],
		void * params);

/**
 * Enumeration to index the dynamic variables in the LALSQTPNGenerator function.
 */
typedef enum {
	LALSQTPN_PHASE,	///< index of the phase
	LALSQTPN_OMEGA,	///< index of the \f$M\omega\f$
	LALSQTPN_LNH_1,	///< index of the \f$\hat{L}_N\f$'s x component
	LALSQTPN_LNH_2,	///< index of the \f$\hat{L}_N\f$'s y component
	LALSQTPN_LNH_3,	///< index of the \f$\hat{L}_N\f$'s z component
	LALSQTPN_CHIH1_1,	///< index of the \f$\hat{\chi}_1\f$'s x component
	LALSQTPN_CHIH1_2,	///< index of the \f$\hat{\chi}_1\f$'s y component
	LALSQTPN_CHIH1_3,	///< index of the \f$\hat{\chi}_1\f$'s z component
	LALSQTPN_CHIH2_1,	///< index of the \f$\hat{\chi}_2\f$'s x component
	LALSQTPN_CHIH2_2,	///< index of the \f$\hat{\chi}_2\f$'s y component
	LALSQTPN_CHIH2_3,	///< index of the \f$\hat{\chi}_2\f$'s z component
	LALSQTPN_MECO,	///< index of the MECO
	LALSQTPN_NUM_OF_VAR	///< number of the dynamic variables
} LALSQTPNGeneratorVariables;

/**
 * The function generates the parameters of the waveform.
 * The formulae are:
 * \f{center}{
 * \begin{gather*}
 * \newcommand{\derTpM}[1]{\dfrac{d#1}{d\left(t/M\right)}}
 * a_1=-0.5A\left(M\omega\right)^{2/3}\left(1+\cos^2\iota\right),\quad
 * a_2=-A\left(M\omega\right)^{2/3}\cos\iota
 * \end{gather*}
 * \f}
 *
 * The waveform-parameters are generated by evolving the system with the LALSQTPNIntegratorFunc() function.
 * The derived values of the dynamic variables [\f$\left(M\omega\right)\f$, \f$\hat{L_N}_x\f$,
 * \f$\hat{L_N}_y\f$, \f$\hat{L_N}_z\f$, \f$\Phi\f$] are calculated by the LALSQTPNDerivator() function.
 * The function returns the amplitude \f$a_1\f$, \f$a_2\f$, polarization
 * shift \f$\alpha\f$, and phase \f$\Phi\f$.
 * The \f$\alpha\f$ is defined with equation (5) in the documentation of the
 * LALSQTPNDerivator() function.
 * @param[out]		waveform	: the generated waveform
 * @param[in]		params		: the input parameters
 */
int XLALSQTPNGenerator(LALSQTPNWave *waveform, LALSQTPNWaveformParams *params);

// LAL wrapper to the XLAL function above
void LALSQTPNGenerator(LALStatus *status, LALSQTPNWave *waveform,
		LALSQTPNWaveformParams *params);

/*@}*/ /* end:LALSQTPNWaveform_h */

#ifdef __cplusplus
}
#endif

#endif /* LALSQTPN_WAVEFORM_H */
