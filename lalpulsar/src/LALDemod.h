/*
*  Copyright (C) 2007 Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve Berukoff, Xavier Siemens
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
#ifndef _LALDEMOD_H
#define _LALDEMOD_H

#include <lal/LALDatatypes.h>
#include <lal/LALComputeAM.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \defgroup LALDemod_h  Header LALDemod.h
 * \ingroup pkg_pulsarCoh
 * \author Berukoff, S.J., Papa, M.A.
 *
 * \brief Computes a demodulated transform.
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/LALDemod.h>
 * \endcode
 *
 * The following is a brief synopsis of the demodulation, or 'Coherent Transform', procedure.
 *
 * In order to remove frequency and amplitude modulation of a time series \f$x_a\f$, we need two basic components:
 * <dl>
 * <dt>Frequency modulation information</dt><dd>  This is given through a phase model \f$\Phi\f$.</dd>
 * <dt>Amplitude modulation information</dt><dd>  This is given through two functions \f$\hat{a}\f$ and
 * \f$\hat{b}\f$, which are derived from the beam-pattern functions \f$F_{+}\f$ and
 * \f$F_{\times}\f$.</dd>
 * </dl>
 * Given these, the F statistic in the \f$b^{th}\f$ frequency bin is
 * \f{equation}{
 * \mathcal{F}_{b} = \frac{4}{S_{h}(f_{0})T_{0}} \frac{B|\hat{F_{a}}|^{2}+A|F_{b}|^{2} - 2C \Re(F_{a}F_{b}^{*})}{D}
 * \f}
 *
 * where
 *
 * \anchor eq_e1 \f{eqnarray}{
 * \tag{eq_e1}
 * \hat{F_{\hat{a}}}=\sum_{a=0}^{NM-1} x_{a} \hat{a} e^{-2{\pi}i{\Phi}_{ab}(\vec{\lambda})} \\
 * \hat{F_{\hat{b}}}=\sum_{a=0}^{NM-1} x_{a} \hat{b} e^{-2{\pi}i{\Phi}_{ab}(\vec{\lambda})}
 * \f}
 * \f$T_{0}\f$ is the observation time, \f$S_{h}\f$ is the noise power spectral density, and \f$A\f$,
 * \f$B\f$, \f$C\f$, and \f$D\f$ are constants.
 *
 * In writing the previous equation we have assumed that there is a total of \f$M\cdot N\f$ data
 * samples and \f$0\leq a<MN\f$.  \f$\Phi_{ab}\f$ is the expected phase at time \f$a\f$ for an
 * intrinsic emission frequency \f$\frac{b}{T_{DeFT}}\f$ (where the denominator is the DeFT time
 * baseline). \f$\Phi\f$ depends on \f$\vec\lambda\f$, a vector of parameters that defines the phase
 * model.  Typically these are the source location and the spin-down parameter values of the template
 * source for which one is demodulating.  For simplicity, we will focus only on \f$F_{a}\f$; the
 * analysis for \f$F_{b}\f$ is identical. Let us now suppose that the time series \f$x_a\f$ is composed
 * of \f$M\f$ chunks, each of \f$N\f$ samples.  If we introduce a short-time index \f$0\leq j<N-1\f$
 * and a short time-series index \f$0\leq \alpha <M-1\f$, so that \f$a=N\alpha+j\f$, we can rewrite the
 * above sum as
 * \anchor eq_e2 \f{equation}{
 * \hat{F_{\hat{a}}}({\vec{\lambda}})=\sum_{\alpha=0}^{M-1}\sum_{j=0}^{N-1}x_{\alpha j}a_{\alpha j}e^{-2{\pi}i{\Phi}_{ab}(\vec{\lambda})}
 * \tag{eq_e2}
 * \f}
 * Note that \f$\hat{a}(t)\f$ is a periodic function with period equal to one sidereal day.  Since the
 * sum over \f$N\f$ is on a timescale much shorter than that (say, 1 hour), then \f$\hat{a}(t)\f$ won't
 * change significantly, and thus can be taken outside of that summation, and then is evaluated at the
 * midpoint of each SFT time.  Now, If \f$\tilde{x}_{\alpha k}\f$ is the matrix of FTs formed along the
 * short time index \f$j\f$
 * \anchor eq_e3 \f{equation}{
 * x_{\alpha j}=\frac{1}{N}\sum_{k=0}^{N-1}\tilde{x}_{\alpha k}e^{2\pi{i}\frac{jk}{N}},
 * \tag{eq_e3}
 * \f}
 * making the appropriate substitutions, Eq.\eqref{eq_e1} becomes
 * \anchor eq_e4 \f{equation}{
 * \hat{F_{\hat{a}}}({\vec{\lambda}})=\sum_{\alpha=0}^{M-1}\hat{a}_{\alpha}\sum_{k=0}^{N-1}\tilde{x}_{\alpha k}\left[\frac{1}{N}\sum_{j=0}^{N-1}e^{-2\pi i(\Phi_{\alpha jb}(\vec{\lambda})-\frac{jk}{N})}\right]
 * \tag{eq_e4}
 * \f}
 * We assume that the phase evolution can be described as linear in \f$t\f$ during the time duration
 * \f$T_{SFT}\f$; thus we can Taylor-expand \f$\Phi\f$ around the temporal midpoint of every SFT time
 * data chunk.  For large values of \f$N\f$, the summation over \f$j\f$ in Eq.\eqref{eq_e4} can be
 * expressed in closed form, thus saving computations, and Eq.\eqref{eq_e4} can be rewritten as
 * \anchor DeFT2 \f{equation}{
 * \hat{F_{\hat{a}}}=\sum_{\alpha=0}^{M-1}\hat{a}_{\alpha}e^{i y_\alpha}\sum_{k=0}^{N-1}\tilde{x}_{\alpha\beta} P_{\alpha k}(b,\vec{\lambda}),
 * \tag{DeFT2}
 * \f}
 * with
 * \anchor DeFT_defs \f{eqnarray}{
 * \tag{DeFT_defs}
 * P_{\alpha k}(b,\vec{\lambda})= \frac{\sin{x'}}{x'}-i \frac{1-\cos{x'}}{x'}\\
 * x'=\sum_{s} f_s B_{s\alpha} - k\\
 * y_\alpha=\sum_{s} f_s A_{s\alpha}.
 * \f}
 * In the previous expressions \f$f_s\f$ indicate the spin-down parameters of different orders (labeled
 * by the index \f$s\f$), and \f$A_{s\alpha}\f$ and \f$B_{s\alpha}\f$ are functions that depend on the phase evolution, whose
 * values depend on \f$\alpha\f$ and on \f$\vec\lambda\f$.
 * [Note that when \f$s=0\f$ the values computed are coefficients of the
 * intrinsic frequency and thus must be computed for the value corresponding to the index.]
 * The values of these functions are
 * calculated by the ComputeSky() routine, also in this package.  Incidentally, in the code,
 * these are the values contained in the variable \c skyConst.
 * Note that the function \f$P_{\alpha k}\f$ is peaked around \f$x'=0\f$.  Thus in the summation
 * over \f$k\f$ in Eq.\eqref{DeFT2} one only needs to consider a few values (NTERMS) of \f$k\f$ around
 * \f$k^*\f$ such that \f$x'(k^*)\approx 0\f$.  This approximation again saves
 * computations. Eq.\eqref{DeFT2} can then be rewritten as
 * \anchor DeFT_algo \f{equation}{
 * \tag{DeFT_algo}
 * \hat{F_{\hat{a}}}=\sum_{\alpha=0}^{M-1}\hat{a}_{\alpha}e^{i y_\alpha}\sum_{k=k^*\pm NTERMS} \tilde x_{\alpha\beta}
 * P_{\alpha k}(b,\vec{\lambda}).
 * \f}
 * If \f$NTERMS\f$ is 8 the power loss due to this approximation is less than \f$\sim 5\%\f$.
 *
 * Now, computing \f$\hat{F_{\hat{a}}}\f$ and \f$\hat{F_{\hat{b}}}\f$ can be done in parallel; given
 * the approximations we have made, for each iteration of the \f$\alpha\f$ loop, one computes first
 * \f$P_{\alpha k}\f$ (through the k-loop), multiplies by \f$\tilde{x}_{\alpha k}\f$, and then forms
 * the statistics of\eqref{eq_e1} at the same time.  After all the iterations of the \f$\alpha\f$ loop
 * are complete, that is, when all SFTs have been exhausted, the final statistic is computed.
 *
 * \heading{Types}
 *
 * \heading{Structure \c DemodPar}
 *
 * This structure contains the parameters for the demodulation routine.   The parameters are:
 *
 * <dl>
 *
 * <dt><tt>INT4 spinDwnOrder</tt></dt><dd> Maximum order of spdwn parameter</dd>
 * <dt><tt>REAL8 *skyConst</tt></dt><dd> The array of sky constants.</dd>
 * <dt><tt>REAL8 *spinDwn</tt></dt><dd> The set of template spindown parameters.</dd>
 * <dt><tt>AMCoeffs *amcoe</tt></dt><dd> The values of the function \f$a\f$ and \f$b\f$, plus their scalar products.</dd>
 * <dt><tt>REAL8 f0</tt></dt><dd> The minimum search frequency</dd>
 * <dt><tt>REAL8 df</tt></dt><dd> The search frequency spacing</dd>
 * <dt><tt>INT4 SFTno</tt></dt><dd> The number of SFTs in coherent timescale</dd>
 * <dt><tt>INT4 Dterms</tt></dt><dd> Terms used in the computation of the dirichlet kernel</dd>
 * <dt><tt>INT4 ifMin</tt></dt><dd> The index of the minimum frequency of the SFT frequency band.</dd>
 * <dt><tt>INT4 imax</tt></dt><dd> How many frequencies are serached.</dd>
 * <dt><tt>BOOLEAN returnFaFb</tt></dt><dd> Wether or not to include the values Fa/Fb in the return-structure Fstat.</dd>
 * </dl>
 *
 * \heading{Structure \c LALFstat}
 *
 * This structure contains the results from LALDemod: either
 * only the value of the \f$\mathcal{F}\f$-statistic \f$F\f$, or also the values
 * of the individual "filters" \f$F_a\f$ and \f$F_b\f$, depending on the
 * <tt>DemodPar->returnFaFb</tt>. \\
 * \e NOTE: the memory has to be allocated before calling <tt>LALDemod()</tt>.
 *
 * <dl>
 * <dt><tt>REAL8 *F</tt></dt><dd>  Array of values of the \f$\mathcal{F}\f$ statistic.</dd>
 * <dt><tt>COMPLEX16 *Fa</tt></dt><dd> Results of match filter with \f$a(t)\f$.</dd>
 * <dt><tt>COMPLEX16 *Fb</tt></dt><dd> Results of match filter with \f$b(t)\f$.</dd>
 * </dl>
 */
/*@{*/

/** \name Error Codes */
/*@{*/
#define LALDEMODH_ENULL 		1	/**< Arguments contained an unexpected null pointer */
/*@}*/
/** \cond DONT_DOXYGEN */
#define LALDEMODH_MSGENULL 	"Arguments contained an unexpected null pointer"
/** \endcond */

#define SMALL	0.000000001

/** PARAMETERS */
typedef struct tagDemodPar {
  INT4		spinDwnOrder;	/* Maximum order of spdwn parameter */
  REAL8		*skyConst;	/* Constants computed in ComputeSky.c */
  REAL8		*spinDwn;	/* Spindown parameter set */
  AMCoeffs      *amcoe;         /*Amplitude Modulation coefficients */
  REAL8         f0;            /*Starting Frequency to be demodulated*/
  REAL8         df;            /*Frequency index resolution*/
  INT4          SFTno;          /* No. of SFTs*/
  INT4          Dterms;         /*Terms used in the computation of the dirichlet kernel*/
  INT4          ifmin;          /*smallest frequency index in SFTs */
  INT4          imax;           /*maximum # of values of F to calculate */
  BOOLEAN	returnFaFb;	/* wether or not to include the values Fa/Fb in the return LALFstat */
}DemodPar;


typedef struct tagLALFstat {
  REAL8         *F;            /* Array of value of the F statistic */
  COMPLEX16     *Fa;           /* Results of match filter with a(t) */
  COMPLEX16     *Fb;           /* Results of match filter with b(t) */
} LALFstat;



/*This structure will hold a single FFT*/
typedef struct tagFFT
{
  COMPLEX8FrequencySeries *fft;
} FFT;






/* Function prototypes */

void LALDemod (LALStatus *status, LALFstat *Fstat, FFT **input, DemodPar *params);
void LALDemodFAST (LALStatus *status, LALFstat *Fstat, FFT **input, DemodPar *params);

/*@}*/

#ifdef __cplusplus
}
#endif
#endif /* _LALDEMOD_H */
