/*** <lalVerbatim file="LALDemodHV">
Author: Berukoff, S.J., Papa, M.A.
$Id$
*** </lalVerbatim> */

/* <lalLaTeX>
\section{Header \texttt{LALDemod.h}}
\label{s:LALDemod.h}
Computes a demodulated transform.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALDemod.h>
\end{verbatim}

\noindent The following is a brief synopsis of the demodulation, or 'Coherent Transform', procedure.

In order to remove frequency and amplitude modulation of a time series $x_a$, we need two basic components:
\begin{description}
\item[\bf{Frequency modulation information}]  This is given through a phase model $\Phi$.
\item[\bf{Amplitude modulation information}]  This is given through two functions $\hat{a}$ and $\hat{b}$, which are derived from the beam-pattern functions $F_{+}$ and $F_{\times}$.
\end{description}
Given these, the F statistic in the $b^{th}$ frequency bin is
\begin{equation}
\mathcal{F}_{b} = \frac{4}{S_{h}(f_{0})T_{0}} \frac{B|\hat{F_{a}}|^{2}+A|F_{b}|^{2} - 2C \Re(F_{a}F_{b}^{*})}{D}
\end{equation}

where

\begin{eqnarray}
\label{e1}
\hat{F_{\hat{a}}}=\sum_{a=0}^{\it{NM-1}} x_{a} \hat{a} e^{-2{\pi}i{\Phi}_{ab}(\vec{\lambda})} \\
\hat{F_{\hat{b}}}=\sum_{a=0}^{\it{NM-1}} x_{a} \hat{b} e^{-2{\pi}i{\Phi}_{ab}(\vec{\lambda})}
\end{eqnarray}
$T_{0}$ is the observation time, $S_{h}$ is the noise power spectral density, and $A$, $B$, $C$, and $D$ are constants.

In writing the previous equation we have assumed that there is a total of $M\cdot N$ data samples and $0\leq a<MN$.  $\Phi_{ab}$ is the expected phase at time $a$ for an  
intrinsic emission frequency $b\over T_{DeFT}$ (where the denominator is the DeFT time baseline). $\Phi$ 
depends on $\vec\lambda$, a vector of parameters that defines the phase model.  Typically these are the source location and the spin-down parameter values of the template source for which one is 
demodulating.  For simplicity, we will focus only on $F_{a}$; the analysis for $F_{b}$ is identical.  
Let us now suppose that the time series $x_a$ is composed of $M$ chunks, each of $N$ samples.  If we introduce a short-time index $0\leq j<N-1$ and a short time-series index $0\leq \alpha <M-1$, so that $a=N\alpha+j$, we can rewrite the above sum as
\begin{equation}
\hat{F_{\hat{a}}}({\vec{\lambda}})=\sum_{\alpha=0}^{M-1}\sum_{j=0}^{N-1}x_{\alpha j}a_{\alpha j}e^{-2{\pi}i{\Phi}_{ab}(\vec{\lambda})}
\label{e2}
\end{equation}
\noindent Note that $\hat{a}(t)$ is a periodic function with period equal to one sidereal day.  Since the sum over $N$ is on a timescale much shorter than that (say, 1 hour), then $\hat{a}(t)$ won't change significantly, and thus can be taken outside of that summation, and then is evaluated at the midpoint of each SFT time.  Now, If $\tilde{x}_{\alpha k}$ is the matrix of FTs formed along the short time index $j$
\begin{equation}
x_{\alpha j}=\frac{1}{N}\sum_{k=0}^{N-1}\tilde{x}_{\alpha k}e^{2\pi{i}\frac{jk}{N}},
\label{e3}
\end{equation}
\noindent  making the appropriate substitutions, Eq.\ref{e1} becomes 
\begin{equation}
\hat{F_{\hat{a}}}({\vec{\lambda}})=\sum_{\alpha=0}^{M-1}\hat{a}_{\alpha}\sum_{k=0}^{N-1}\tilde{x}_{\alpha k}\left[\frac{1}{N}\sum_{j=0}^{N-1}e^{-2\pi i(\Phi_{\alpha jb}(\vec{\lambda})-\frac{jk}{N})}\right]
\label{e4}
\end{equation}
We assume that the phase evolution can be described as linear in $t$ during the time duration 
$T_{SFT}$; thus we can Taylor-expand $\Phi$ around the temporal midpoint of every SFT time data chunk.  For large values of $N$,
the summation over $j$ in eq. (\ref{e4}) can be expressed in closed form, thus saving computations, and eq. (\ref{e4}) can be rewritten as
\begin{equation}
\hat{F_{\hat{a}}}=\sum_{\alpha=0}^{M-1}\hat{a}_{\alpha}e^{i y_\alpha}\sum_{k=0}^{N-1}\tilde{x}_{\alpha\beta} P_{\alpha k}(b,\vec{\lambda}),
\label{DeFT2} 
\end{equation}
with
\begin{eqnarray}
\label{DeFT_defs}
P_{\alpha k}(b,\vec{\lambda})={\sin{x'}\over x'}-i{1-\cos{x'}\over x'}\\
x'=\sum_{s} f_s B_{s\alpha} - k\\
y_\alpha=\sum_{s} f_s A_{s\alpha}.
\end{eqnarray}
In the previous expressions $f_s$ indicate the spin-down parameters of different orders (labeled 
by the index $s$\footnote{Note that when $s=0$ the values computed are coefficients of the intrinsic frequency and thus must be computed for the value corresponding to the index $b$.}), and $A_{s\alpha}$ and $B_{s\alpha}$ 
are functions that depend on the phase evolution, whose values depend on $\alpha$ and on $\vec\lambda$.  The values of these functions are calculated by the \verb@ComputeSky()@ routine, also in this package.  Incidentally, in the code, these are the values contained in the variable \verb@skyConst@.  
Note that the function $P_{\alpha k}$ is peaked around $x'=0$.  Thus in the summation 
over $k$ in eq. (\ref{DeFT2}) one only needs to consider a few values (NTERMS) of $k$ around $k^*$ such 
that $x'(k^*)\approx 0$.  This approximation again saves computations. Eq. (\ref{DeFT2}) can then be rewritten as
\begin{equation}
\label{DeFT_algo}
\hat{F_{\hat{a}}}=\sum_{\alpha=0}^{M-1}\hat{a}_{\alpha}e^{i y_\alpha}\sum_{k=k^*\pm NTERMS} \tilde x_{\alpha\beta} 
P_{\alpha k}(b,\vec{\lambda}).
\end{equation} 
If $NTERMS$ is 8 the power loss due to this approximation is less than $\sim 5\%$. 
 
Now, computing $\hat{F_{\hat{a}}}$ and $\hat{F_{\hat{b}}}$ can be done in parallel; given the approximations we have made, for each iteration of the $\alpha$ loop, one computes first $P_{\alpha k}$ (through the k-loop), multiplies by $\tilde{x}_{\alpha k}$, and then forms the statistics of \ref{e1} at the same time.  After all the iterations of the $\alpha$ loop are complete, that is, when all SFTs have been exhausted, the final statistic is computed.

 </lalLaTeX> */


#ifndef _LALDEMOD_H
#define _LALDEMOD_H

#include <lal/LALDatatypes.h>
#include <lal/LALComputeAM.h>

#ifdef __cplusplus
extern "C" {
#endif

NRCSID (LALDEMODH, "$Id$"); 

/********************************************************** <lalLaTeX>
\subsection*{Error codes}
</lalLaTeX>
***************************************************** <lalErrTable> */
#define LALDEMODH_ENULL 		1

#define LALDEMODH_MSGENULL 	"Arguments contained an unexpected null pointer"

/*************************************************** </lalErrTable> */

#define SMALL	0.000000001

/* <lalLaTeX>

\subsection*{Types}

\subsubsection*{Structure \texttt{DemodPar}}
\idx[Type]{DemodPar}

\noindent This structure contains the parameters for the demodulation routine.   The parameters are:

\begin{description}

\item[\texttt{INT4 spinDwnOrder}] Maximum order of spdwn parameter
\item[\texttt{REAL8 *skyConst}] The array of sky constants.
\item[\texttt{REAL8 *spinDwn}] The set of template spindown parameters.
\item[\texttt{AMCoeffs *amcoe}] The values of the function $a$ and $b$, plus their scalar products.
\item[\texttt{REAL8 f0}] The minimum search frequency 
\item[\texttt{REAL8 df}] The search frequency spacing
\item[\texttt{INT4 SFTno}] The number of SFTs in coherent timescale
\item[\texttt{INT4 Dterms}] Terms used in the computation of the dirichlet kernel 
\item[\texttt{INT4 ifMin}] The index of the minimum frequency of the SFT frequency band.
\item[\texttt{INT4 imax}] How many frequencies are serached.
\item[\texttt{BOOLEAN returnFaFb}] Wether or not to include the values Fa/Fb in the return-structure Fstat.
\end{description}


\subsubsection*{Structure \texttt{LALFstat}}
\idx[Type]{LALFstat}

\noindent This structure contains the results from LALDemod: either
only the value of the $\mathcal{F}$-statistic $F$, or also the values
of the individual "filters" $F_a$ and $F_b$, depending on the
\texttt{DemodPar->returnFaFb}. \\
\emph{NOTE:} the memory has to be allocated before calling \texttt{LALDemod()}.

\begin{description}
\item[\texttt{REAL8 *F}]  Array of values of the $\mathcal{F}$ statistic. 
\item[\texttt{COMPLEX16 *Fa}] Results of match filter with $a(t)$.
\item[\texttt{COMPLEX16 *Fb}] Results of match filter with $b(t)$.
\end{description}


</lalLaTeX> */
 

/* PARAMETERS */
typedef struct DemodParTag{
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


typedef struct {
  REAL8         *F;            /* Array of value of the F statistic */
  COMPLEX16     *Fa;           /* Results of match filter with a(t) */
  COMPLEX16     *Fb;           /* Results of match filter with b(t) */
} LALFstat;



/*This structure will hold a single FFT*/
typedef struct FFTTag 
{
  COMPLEX8FrequencySeries *fft;   
} FFT;

/********************************************************** <lalLaTeX>
\vfill{\footnotesize\input{LALDemodHV}}
\newpage\input{LALDemodC}
******************************************************* </lalLaTeX> */

/* Function prototypes */

void LALDemod (LALStatus *stat, LALFstat *Fstat, FFT **input, DemodPar *params);

	
/* <lalLaTeX>
\newpage\input{LALDemodTestC}
</lalLaTeX> */

#ifdef __cplusplus
}
#endif
#endif /* _LALDEMOD_H */
