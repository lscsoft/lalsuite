/*
*  Copyright (C) 2007 Duncan Brown, Eirini Messaritaki, Jolien Creighton, Thomas Cokelaer
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

/*-----------------------------------------------------------------------
 *
 * File Name: FindChirpBCV.h
 *
 * Author: Brown, D. A. and Messaritaki, E.
 *
 *-----------------------------------------------------------------------
 */

#ifndef _FINDCHIRPBCVH_H
#define _FINDCHIRPBCVH_H

#include <lal/LALDatatypes.h>
#include <lal/RealFFT.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpChisq.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif


/**
 * \addtogroup FindChirpBCV_h
 * \author Brown, D. A. and Messaritaki, E.
 *
 * \brief Provides structures and functions to condition interferometer data
 * and generate binary inspiral chirps using the BCV detection template
 * family.
 *
 * \section sec_fcbcv_synopsis Synopsis
 * \code
 * #include <lal/FindChirpBCV.h>
 * \endcode
 *
 * \subsection sec_fcbcv_bh Binary Black Holes
 *
 * For the binary black hole inspiral we use the BCV templates.
 *
 * \subsubsection sec_fcbcv_BCV BCV Templates
 *
 * The signal-to-noise ratio (SNR) for a signal \f$s\f$ and a template \f$h\f$ is given by
 * \f{equation}{
 * \rho(h) = \frac{<s,h>}{\sqrt{<h,h>}},
 * \f}
 * with the inner product \f$<s,h>\f$ being defined as
 * \anchor SNRdef \f{equation}{
 * <s,h> = 2 \int_{-\infty}^{\infty} \frac{\tilde{s}^{\ast}(f) \tilde{h}(f)}
 * {S_h(|f|)} df =
 * 4 \Re \int_0^{\infty} \frac{\tilde{s}^{\ast}(f) \tilde{h}(f)}{S_h(f)} df
 * \tag{SNRdef}
 * \f}
 * and \f$S_h(f)\f$ being the one-sided noise power spectral density.
 * The last equality in Eq.\eqref{SNRdef} holds only if \f$\tilde{s}(f)\f$ and
 * \f$\tilde{h}(f)\f$ are the Fourier-transforms of real time-functions.
 *
 * The effective frequency-domain template given by Buonanno, Chen and Vallisneri
 * is
 * \anchor eq_tmplt \f{equation}{
 * \tilde{h}(f) = A(f) e^{i \psi(f)}
 * \tag{eq_tmplt}
 * \f}
 * where
 * \f{equation}{
 * A(f) = f^{-7/6} (1-\alpha f^{2/3}) \theta(f_{cut}-f),
 * \f}
 * \f{equation}{
 * \psi(f) = \phi_0 + 2 \pi f t_0 + f^{-5/3} (\psi_0 + \psi_{1/2} f^{1/3} +
 * \psi_1 f^{2/3} + \psi_{3/2} f + \ldots).
 * \f}
 * In these expressions, \f$t_0\f$ and \f$\phi_0\f$ are the time of arrival and the
 * frequency-domain phase offset respectively,
 * and \f$\theta\f$ is the Heaviside step function.
 * For most inspiral templates approximated by the template\eqref{eq_tmplt},
 * it is sufficient to use the parameters \f$\psi_0\f$ and \f$\psi_{3/2}\f$ and set all
 * other \f$\psi\f$ coefficients equal to 0.
 * So in the following:
 * \f{eqnarray}{
 * \psi(f) &=& \phi_0 + 2 \pi f t_0 + f^{-5/3} (\psi_0 + \psi_{3/2} f)\\
 * &=& \phi_0 + \psi'(f) . \\
 * \f}
 *
 * To simplify the equations, the abbreviation
 * \f{equation}{
 * I_k \equiv 4 \int_0^{f_{cut}} \frac{df}{f^k S_h(f)}
 * \f}
 * is used in the following.
 *
 * Notice that in the code, \f$\psi_0\f$ is \c psi0 and \f$\psi_{3/2}\f$ is
 * \c psi3.
 *
 * \subsubsection NormOfTemplate Normalized Template
 *
 * We begin by normalizing the template \f$\tilde{h}(f)\f$.
 * Specifically, it is assumed that the normalized template is
 * \f{equation}{
 * \hat{h}(f) = N \tilde{h}(f)
 * \f}
 * where \f$N\f$ is a real number. Then:
 * \anchor Normalize \f{eqnarray}{
 * && <\hat{h}, \hat{h}> = 1 \: \Rightarrow  \:
 * 4 \Re \int_0^{\infty} \frac{\hat{h}^{\ast} \hat{h}}{S_h(f)} df = 1 \:
 * \Rightarrow  \\
 * && 4 N^2 \int_0^{\infty} \frac{ [ f^{-7/6} (1-\alpha f^{2/3})]^2 \theta
 * (f_{cut}-f) }{S_h} df = 1 \: \Rightarrow \\
 * && N = \sqrt{ I_{7/3} - 2 \: \alpha \: I_{5/3} + \alpha^2 \: I_1 }
 * \tag{Normalize}
 * \f}
 * So the normalized template is
 * \f{equation}{
 * \hat{h}(f) = \frac{1}{\sqrt{ I_{7/3} - 2 \: \alpha \: I_{5/3} +
 * \alpha^2 \: I_1 }} f^{-7/6} (1-\alpha f^{2/3}) e^{i \phi_0} e^{i \psi'}
 * \theta(f_{cut}-f), \: f>0
 * \f}
 * and \f$\hat{h}(f) = \hat{h}^{\ast}(-f), \: f<0\f$.
 *
 * Next we construct an orthonormal basis \f$\{ \hat{h}_j \}\f$ for the 4-dimensional
 * linear subspace of templates, with \f$\phi_0 \in [0,2\pi)\f$ and
 * \f$\alpha \in (-\infty, \infty)\f$ and all other parameters fixed.
 * Specifically, we want the basis vectors to satisfy
 * \anchor OrthonormBasis \f{equation}{
 * < \hat{h}_i , \hat{h}_j > = \delta_{ij}.
 * \tag{OrthonormBasis}
 * \f}
 * For that we construct two real functions \f$A_1(f)\f$ and \f$A_2(f)\f$, linear
 * combinations of \f$f^{-7/6}\f$ and \f$f^{-1/2}\f$, which are related to the 4 basis
 * vectors via:
 * \f{eqnarray}{
 * \hat{h}_{1,2}(f) &=& A_{1,2}(f) \: e^{i \psi'(f)} \\
 * \hat{h}_{3,4}(f) &=& A_{1,2}(f) \: i \: e^{i \psi'(f)}.
 * \f}
 * Then, Eq.\eqref{OrthonormBasis} becomes:
 * \anchor OrthonormA \f{equation}{
 * 4 \Re \int_0^{\infty} \frac{A_i(f) A_j(f)}{S_h} df = \delta_{ij}.
 * \tag{OrthonormA}
 * \f}
 * So we choose:
 *
 * \f{eqnarray}{
 * A_1(f) &=& a_1 f^{-7/6} \\
 * A_2(f) &=& b_1 f^{-7/6} + b_2 f^{-1/2}.
 * \f}
 * Imposing condition\eqref{OrthonormA} gives:
 * \anchor A1 \anchor A2 \anchor A1A2 \f{eqnarray}{
 * \tag{A1}
 * && 4 \int_0^{\infty} \frac{A_1(f) A_1(f)}{S_h} df = 1 \: \Rightarrow \:
 * a_1 = I_{7/3}^{-1/2} \\
 * \tag{A2}
 * && 4 \int_0^{\infty} \frac{A_2(f) A_2(f)}{S_h} df = 1 \: \Rightarrow \:
 * b_1^2\: I_{7/3} + 2\: b_1\: b_2 \: I_{5/3} +
 * b_2^2 \: I_1 = 1 \\
 * \tag{A1A2}
 * && 4 \int_0^{\infty} \frac{A_1(f) A_2(f)}{S_h} df = 0 \: \Rightarrow \:
 * b_1 = - b_2 \frac{I_{5/3}}{I_{7/3}}
 * \f}
 * Solving Eqs.\eqref{A2} and\eqref{A1A2} we get
 * \f{eqnarray}{
 * b_1 &=& - \frac{I_{5/3}}{I_{7/3}} \Big ( I_1 - \frac{I_{5/3}^2}{I_{7/3}}
 * \Big )^{-1/2} \\
 * b_2 &=& \Big ( I_1 - \frac{I_{5/3}^2}{I_{7/3}} \Big )^{-1/2}.
 * \f}
 *
 * The next step is to write the normalized template in terms of the 4 basis
 * vectors
 * \f{eqnarray}{
 * \hat{h}(f) &=& c_1 \hat{h}_1(f) + c_2 \hat{h}_2(f) + c_3 \hat{h}_3(f) +
 * c_4 \hat{h}_4(f) \: \Rightarrow \\
 * \frac{(f^{-7/6} - \alpha f^{-1/2}) e^{i \phi_0}}{\sqrt{I_{7/3} - 2 \alpha
 * I_{5/3} + \alpha^2 I_1}} &=& (c_1 +i c_3) a_1 f^{-7/6} +
 * (c_2 + i c_4) (b_1 f^{-7/6} + b_2 f^{-1/2} )
 * \f}
 * and matching the terms gives
 * \f{eqnarray}{
 * c_1 &=& \cos\phi_0 \cos\omega \\
 * c_2 &=& \cos\phi_0 \sin\omega \\
 * c_3 &=& \sin\phi_0 \cos\omega \\
 * c_4 &=& \sin\phi_0 \sin\omega
 * \f}
 * if the angle \f$\omega\f$ is defined by
 * \f{equation}{
 * \tan\omega \equiv - \frac{a_1 \alpha}{b_2 + b_1 \alpha}.
 * \f}
 * So
 * \f{equation}{
 * \hat{h}(f) = \cos\phi_0 \cos\omega \:\hat{h}_1(f) +
 * \cos\phi_0 \sin\omega \:\hat{h}_2(f) +
 * \sin\phi_0 \cos\omega \:\hat{h}_3(f) +
 * \sin\phi_0 \sin\omega \:\hat{h}_4(f).
 * \f}
 *
 * \subsubsection SNRMaximization Maximization of the SNR
 *
 * The Fourier-transformed data is \f$s(f)\f$. The SNR is
 * \f{eqnarray}{
 * \rho &=& <s,\hat{h}(f)> \\
 * &=& \cos\phi_0 \cos\omega \: K_1 +
 * \cos\phi_0 \sin\omega \: K_2+
 * \sin\phi_0 \cos\omega \: K_3 +
 * \sin\phi_0 \sin\omega \: K_4
 * \f}
 * where the 4 integrals are defined by
 * \anchor K1 \anchor K2 \anchor K3 \anchor K4 \f{eqnarray}{
 * \tag{K1}
 * K_1 &=& <s,\hat{h}_1>  = \Re \int_0^{f_{cut}} \frac{4 s^{\ast} a_1 f^{-7/6}
 * e^{i \psi'}}{S_h} df \\ \tag{K2}
 * K_2 &=& <s,\hat{h}_2>  = \Re \int_0^{f_{cut}} \frac{4 s^{\ast} (b_1 f^{-7/6}
 * + b_2 f^{-1/2} ) e^{i \psi'}}{S_h} df \\ \tag{K3}
 * K_3 &=& <s,\hat{h}_3>  = \Re \int_0^{f_{cut}} \frac{4 s^{\ast} a_1 f^{-7/6} i
 * e^{i \psi'}}{S_h} df
 * = -\Im \int_0^{f_{cut}} \frac{4 s^{\ast} a_1 f^{-7/6} e^{i \psi'}}{S_h}
 * df  \\ \tag{K4}
 * K_4 &=& <s,\hat{h}_4>  = \Re \int_0^{f_{cut}} \frac{4 s^{\ast} (b_1 f^{-7/6}
 * + b_2 f^{-1/2} ) i e^{i \psi'}}{S_h} df
 * = -\Im \int_0^{f_{cut}} \frac{4 s^{\ast} (b_1 f^{-7/6} + b_2 f^{-1/2})
 * e^{i \psi'}}{S_h} df
 * \f}
 * Now set
 * \anchor OmegaMinusPhi \anchor OmegaPlusPhi \f{eqnarray}{
 * \tag{OmegaMinusPhi}
 * A &=& \omega - \phi_0 \\
 * B &=& \omega + \phi_0 \tag{OmegaPlusPhi}
 * \f}
 * so that the expression for the SNR becomes
 * \anchor MaxSNR \f{eqnarray}{
 * \rho &=& \frac{1}{2} K_1 [\cos(\omega+\phi_0)+\cos(\omega-\phi_0)] +
 * \frac{1}{2} K_2 [\sin(\omega+\phi_0)+\sin(\omega-\phi_0)] +\\
 * && \frac{1}{2} K_3[\sin(\omega+\phi_0)-\sin(\omega-\phi_0)]
 * +\frac{1}{2} K_4 [\cos(\omega-\phi_0)-\cos(\omega+\phi_0)]\Rightarrow\\
 * 2 \rho &=& (K_1+K_4) \cos A + (K_2-K_3) \sin A + (K_1-K_4) \cos B + (K_2+K_3)
 * \sin B.
 * \tag{MaxSNR}
 * \f}
 *
 * To maximize with respect to \f$A\f$ we take the first derivative
 * \f{equation}{
 * \frac{\partial(2 \rho)}{\partial A} = - (K_1+K_4) \sin A
 * +(K_2-K_3) \cos A
 * \f}
 * and set that equal to 0, which gives
 * \anchor TANA \f{equation}{
 * \frac{\partial(2 \rho)}{\partial A} \Big |_{A_0} = 0 \: \Rightarrow
 * \: \tan A_0 =\frac{K_2-K_3}{K_1+K_4}. \\
 * \tag{TANA}
 * \f}
 * Then the sine and cosine of \f$A_0\f$ can be found:
 * \anchor SINA \anchor COSA \f{eqnarray}{
 * \sin A_0 &=& \pm \frac{\tan A_0}{\sqrt{1 + \tan^2 A_0}} = \pm
 * \frac{K_2-K_3}{\sqrt{(K_1+K_4)^2 +(K_2-K_3)^2}},  \tag{SINA} \\
 * \cos A_0 &=& \pm \frac{1}{\sqrt{1 + \tan^2 A_0}} = \pm
 * \frac{K_1+K_4}{\sqrt{(K_1+K_4)^2 +(K_2-K_3)^2}}.
 * \tag{COSA}
 * \f}
 * Notice that for Eq.\eqref{TANA} to be satisfied, the same sign must be kept
 * in Eqs.\eqref{SINA} and\eqref{COSA}.
 * To find the values that correspond to the maximum, we take the second
 * derivative of \f$\rho\f$ with respect to \f$A\f$:
 * \f{equation}{
 * \frac{\partial^2(2 \rho)}{\partial A^2}\Big |_{A_0} < 0 \: \Rightarrow
 * \: \Big [ - (K_1+K_4) \cos A - (K_2-K_3) \sin A \Big ]_{A_0} < 0
 * \f}
 * which is satisfied if the \f$+\f$ sign is considered in Eqs.\eqref{SINA} and
 * \eqref{COSA}.
 *
 * To maximize with respect to \f$B\f$ we take the first derivative
 * \f{equation}{
 * \frac{\partial(2 \rho)}{\partial B}  = - (K_1-K_4) \sin B
 * +(K_2+K_3) \cos B
 * \f}
 * and set that equal to 0, which gives
 * \anchor TANB \f{equation}{
 * \frac{\partial(2 \rho)}{\partial B} \Big |_{B_0} = 0 \: \Rightarrow
 * \: \tan B_0 =\frac{K_2+K_3}{K_1-K_4}. \\
 * \tag{TANB}
 * \f}
 * Then the sine and cosine of \f$B_0\f$ can be found:
 * \anchor SINB \anchor COSB \f{eqnarray}{
 * \sin B_0 &=& \pm \frac{\tan B_0}{\sqrt{1 + \tan^2 B_0}} = \pm
 * \frac{K_2+K_3}{\sqrt{(K_1-K_4)^2 +(K_2+K_3)^2}},  \tag{SINB} \\
 * \cos B_0 &=& \pm \frac{1}{\sqrt{1 + \tan^2 B_0}} = \pm
 * \frac{K_1-K_4}{\sqrt{(K_1-K_4)^2 +(K_2+K_3)^2}}.
 * \tag{COSB}
 * \f}
 * Again, the same sign must be kept in Eqs.\eqref{SINB} and\eqref{COSB}.
 * To find the values that correspond to the maximum, we take the second
 * derivative of \f$\rho\f$ with respect to \f$B\f$:
 * \f{equation}{
 * \frac{\partial^2(2 \rho)}{\partial B^2}\Big |_{B_0} < 0 \: \Rightarrow
 * \: \Big [ - (K_1-K_4) \cos B - (K_2+K_3) \sin B \Big ]_{B_0} < 0
 * \f}
 * which is satisfied if the \f$+\f$ sign is considered in Eqs.\eqref{SINB} and
 * \eqref{COSB}.
 *
 * Substituting the expressions for the sines and cosines of \f$A_0\f$ and
 * \f$B_0\f$
 * into Eq.\eqref{MaxSNR}, the maximum SNR is:
 * \f{eqnarray}{
 * \rho_{max} &=& \frac{1}{2} \sqrt{(K_1+K_4)^2 + (K_2-K_3)^2} +
 * \frac{1}{2} \sqrt{(K_1-K_4)^2 +(K_2+K_3)^2}\\
 * 2\rho_{max}&=& \sqrt{ K_1^2 + K_2^2 + K_3^2 +K_4^2 + 2(K_1 K_4
 * - K_2 K_3)} + \sqrt{  K_1^2 + K_2^2 + K_3^2 +K_4^2 -2(K_1
 * K_4 -  K_2 K_3)}
 * \f}
 * To achieve a simpler form for the SNR, we can use Eqs.\eqref{K1}-\eqref{K4} to
 * combine the integrals \f$K_1\f$, \f$K_2\f$, \f$K_3\f$ and \f$K_4\f$.
 * Specifically:
 * \f{equation}{
 * K_1^2 + K_2^2 + K_3^2 + K_4^2 = \Big | \int_0^{f_{cut}} \frac{4 s^{\ast} a_1
 * f^{-7/6} e^{i \psi'}}{S_h} df \Big |^2 + \Big | \int_0^{f_{cut}}
 * \frac{4 s^{\ast} (b_1 f^{-7/6} + b_2 f^{-1/2}) e^{i \psi'}}{S_h} df
 * \Big |^2
 * \f}
 * and
 * \f{equation}{
 * 2(K_1 K_4- K_2 K_3) =  2 \Im \Big \{ \int_0^{f_{cut}} \frac{4 s^{\ast} a_1
 * f^{-7/6} e^{i \psi'}}{S_h} df \Big (\int_0^{f_{cut}} \frac{4 s^{\ast}
 * (b_1 f^{-7/6} + b_2 f^{-1/2}) e^{i \psi'}}{S_h}df\Big )^{\ast} \Big \}.
 * \f}
 *
 * \subsubsection ChisquaredVeto The \f$\chi^2\f$-veto
 *
 * If we are working with \f$p\f$ bins the maximum SNR for the template must be
 * divided into \f$p\f$ equal parts. In this case, since we have two different
 * amplitude-parts of the template, we have to calculate two sets of bin
 * boundaries.
 * For the first set, the quantity
 * \f[
 * \int_0^{f_{cut}} \frac{4 a_1^2 f^{-7/3}}{S_h(f)} df
 * \f]
 * must be divided into \f$p\f$ equal pieces.
 * For the second set, the quantity
 * \f[
 * \int_0^{f_{cut}} \frac{4 [b_1 f^{-7/6} + b_2 f^{-1/2}]^2}{S_h(f)} df
 * \f]
 * must be divided into \f$p\f$ equal pieces.
 *
 * To check if the total SNR is smoothly distributed over the bins,
 * take:
 * \f{eqnarray}{
 * \chi^2 &=& \: p \sum_{l=1}^p  \Big | \int_{f_l}^{f_{l+1}} \frac{4 a_1
 * f^{-7/6} \tilde{s}^{\ast} e^{i \psi'}}{S_h(f)} df
 * - \frac{1}{p} \int_{0}^{f_{cut}} \frac{4 a_1 f^{-7/6}
 * \tilde{s}^{\ast} e^{i \psi'}}{S_h(f)} df \Big |^2 \\
 * &+&\: p\sum_{l=1}^p\Big | \int_{f_l}^{f_{l+1}} \frac{4[b_1 f^{-7/6}
 * + b_2 f^{-1/2} ]
 * \tilde{s}^{\ast} e^{i \psi'}}{S_h(f)} df-\frac{1}{p}
 * \int_{0}^{\infty} \frac{4[b_1 f^{-7/6} + b_2 f^{-1/2}]
 * \tilde{s}^{\ast} e^{i \psi'}}{S_h(f)} df\Big |^2. \\
 * \f}
 *
 * \subsubsection sec_fcbcv_quants Quantities calculated in the code for the case of BCV templates
 *
 * <ol>
 * <li> The template normalization squared (in \c LALFindChirpBCVTemplate):
 * \f{equation}{
 * \mathtt{tmpltNorm} = d^2 \Big ( \frac{5\mu}{96 M_{\odot}}\Big ) \Big ( \frac{M}        {\pi^2 M_{\odot}} \Big )^{2/3} T_{\odot}^{-1/3}
 * \Big (\frac{2 T_{\odot} c}{1 Mpc} \Big )^2.
 * \f}</li>
 * <li> The exponential \f$e^{i \psi'}\f$ (in \c LALFindChirpBCVTemplate).</li>
 * <li> The BCV Moments (in \c LALFindChirpBCVData): \\
 * \f{eqnarray}{
 * \mathtt{I73} &=& \sum_{k=0}^{N/2} \frac{4 k^{-7/3}}{|dR|^2 S_v(|f_k|)} \\
 * \mathtt{I53} &=& \sum_{k=0}^{N/2} \frac{4 k^{-5/3}}{|dR|^2 S_v(|f_k|)} \\
 * \mathtt{I1} &=& \sum_{k=0}^{N/2} \frac{4 k^{-1}}{|dR|^2 S_v(|f_k|)}. \\
 * \f}
 * These quantities should be multiplied by \f$(\Delta t /N)\f$ (b/c of the FT) and
 * by (\f$\mathtt{tmpltNorm}\f$) (b/c of the template), but that is taken
 * care of in the calculation of the SNR.</li>
 * <li> The BCV normalization factors (in \c LALFindChirpBCVData):
 * \f{eqnarray}{
 * \mathtt{a1} &=& (\mathtt{I73})^{-1/2} \\
 * \mathtt{b2} &=& \Big( \mathtt{I1} - \frac{\mathtt{I53}^2}{\mathtt{I73}}
 * \Big )^{-1/2} \\
 * \mathtt{b1} &=& -\frac{\mathtt{I53}}{\mathtt{I73}}  \mathtt{b2}. \\
 * \f}
 * Again, these should be multiplied by \f$[(\Delta t/N)^{-1/2} \mathtt{tmpltNorm}
 * )^{-1/2}]\f$
 * but that is taken care of later.</li>
 * <li> The two FTs required for the calculation of the SNR (in
 * \c LALFindChirpBCVFilter):
 * \f{eqnarray}{
 * qj &=&\sum_{k=0}^{N/2}e^{2 \pi i j k/N} \frac{(dR \tilde{v}_k^{\ast})
 * \: \mathtt{a1} \: k^{-7/6} e^{i \psi'} }{ |dR|^2 S_v(|f_k|) } \\
 * q^{{}_{BCV}}_j &=& \sum_{k=0}^{N/2}e^{ 2 \pi i j k/N} \frac{(dR \tilde{v}_k
 * ^{\ast})( \mathtt{b1} \:k^{-7/6} + \mathtt{b2}\: k^{-1/2} ) e^{i \psi'}}
 * {|dR|^2 S_v(|f_k|)} \\
 * \f}
 * up to the appropriate normalization factors, namely \f$(\mathtt{tmpltNorm})
 * ^{1/2}( \Delta t/N)\f$</li>
 * <li> The SNR (in \c LALFindChirpBCVFilter):
 * \f{equation}{
 * \rho^2(t_j) = \Big ( \frac{\Delta t}{N} \Big ) \bigg \{ \frac{1}{2}
 * \sqrt{ |q_j|^2 + |q^{{}_{BCV}}_j|^2 + 2 \Im (q_j q^{{}_{BCV}\ast}_j) } +
 * \frac{1}{2}\sqrt{|q_j|^2+|q^{{}_{BCV}}_j|^2 -
 * 2 \Im (q_j q^{{}_{BCV} \ast}_j) }
 * \bigg \}.
 * \f}</li>
 * </ol>
 *
 */
/*@{*/


/**\name Error Codes */
/*@{*/
#define FINDCHIRPBCVH_ENULL 1
#define FINDCHIRPBCVH_ENNUL 2
#define FINDCHIRPBCVH_EALOC 3
#define FINDCHIRPBCVH_ENUMZ 4
#define FINDCHIRPBCVH_ESEGZ 5
#define FINDCHIRPBCVH_EMISM 6
#define FINDCHIRPBCVH_EDELT 7
#define FINDCHIRPBCVH_EFLOW 8
#define FINDCHIRPBCVH_EDYNR 9
#define FINDCHIRPBCVH_EISTN 10
#define FINDCHIRPBCVH_EDIVZ 11
#define FINDCHIRPBCVH_EMAPX 12
#define FINDCHIRPBCVH_EUAPX 13
#define FINDCHIRPBCVH_EZNRM 14
#define FINDCHIRPBCVH_EQLEN 15
#define FINDCHIRPBCVH_ECLUW 16
/*@}*/
/*@}*/

#define FINDCHIRPBCVH_MSGENULL "Null pointer"
#define FINDCHIRPBCVH_MSGENNUL "Non-null pointer"
#define FINDCHIRPBCVH_MSGEALOC "Memory allocation error"
#define FINDCHIRPBCVH_MSGENUMZ "Invalid number of segments"
#define FINDCHIRPBCVH_MSGESEGZ "Invalid number of points in segments"
#define FINDCHIRPBCVH_MSGEMISM "Mismatch between number of points in segments"
#define FINDCHIRPBCVH_MSGEDELT "deltaT is zero or negative"
#define FINDCHIRPBCVH_MSGEFLOW "Low frequency cutoff is negative"
#define FINDCHIRPBCVH_MSGEDYNR "Dynamic range scaling is zero or negative"
#define FINDCHIRPBCVH_MSGEISTN "Truncation of inverse power spectrum is negative"
#define FINDCHIRPBCVH_MSGEDIVZ "Attempting to divide by zero"
#define FINDCHIRPBCVH_MSGEMAPX "Mismatch in waveform approximant"
#define FINDCHIRPBCVH_MSGEUAPX "Unknown approximant"
#define FINDCHIRPBCVH_MSGEZNRM "No non-zero value assigned to one of a1, b1, b2"
#define FINDCHIRPBCVH_MSGEQLEN "params->qVec->length not equal to params->qVecBCV->length"
#define FINDCHIRPBCVH_MSGECLUW "Unacceptable max-over-chirp clustering method for BCV"


void
LALFindChirpBCVData (
    LALStatus                  *status,
    FindChirpSegmentVector     *fcSegVec,
    DataSegmentVector          *dataSegVec,
    FindChirpDataParams        *params
    );

void
LALFindChirpBCVTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *theTmplt,
    FindChirpTmpltParams       *params
    );

void
LALFindChirpBCVChisqVeto (
    LALStatus                  *status,
    REAL4Vector                *chisqVec,
    FindChirpChisqInput        *input,
    FindChirpChisqInput        *inputBCV,
    FindChirpChisqParams       *params
    );

void
LALFindChirpBCVFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    );

void
LALFindChirpBCVCFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    );


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _FINDCHIRPSPH_H */
