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
Given these, the long time baseline demodulated Fourier Transform (DeFT) in the $b^{th}$ frequency bin is
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Random.h>
#include <lal/RealFFT.h>
#include <lal/LALConstants.h>
#include <lal/DetResponse.h>
#include <lal/DetectorSite.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GenerateTaylorCW.h>
#include <lal/LALComputeAM.h>
#include <lal/ComputeSky.h>
#include <lal/LALBarycenter.h>

#ifdef __cplusplus
#extern "C" {
#endif

NRCSID (LALDEMODH, "$Id$"); 

/* <lalLaTeX>
\subsection*{Error conditions}
\vspace{0.1in}
\input{LALDemodHErrorTable}
</lalLaTeX> */
/* *****<lalErrTable file="LALDemodHErrorTable">*/
#define LALDEMODH_ENULL 1
#define LALDEMODH_ENNEG 2
#define LALDEMODH_EORDR 4
#define LALDEMODH_ENOFILE 8
#define LALDEMODH_EBADARG 16
#define LALDEMODH_MSGENULL "Null Pointer"
#define LALDEMODH_MSGENNEG "Erroneous Negative Value"
#define LALDEMODH_MSGEORDR "Min Value Larger than Max"
#define LALDEMODH_MSGENOFILE "Non-existent filename"
#define LALDEMODH_MSGEBADARG "Bad Command Line argument"
/********************************************** </lalErrTable> */


#define NTERM_COH_DIV_TWO 4 
#define SMALL	0.000000001
#define SFTOVERLAP 0.05

/* <lalLaTeX>

\subsection*{Structures}

\begin{verbatim}
struct DemodPar
\end{verbatim}
\index{\texttt{DemodPar}}

\noindent This structure contains the parameters for the demodulation routine.   The parameters are:

\begin{description}
\item[\texttt{INT4 if0Max}] The maximum search frequency index in coarse-grained frequency resolution.
\item[\texttt{INT4 if0Min}] The minimum search frequency index in coarse-grained frequency resolution.
\item[\texttt{INT4 mCohSFT}] The number of SFTs in coherent timescale.  In other words, the number of SFTs which make up one DeFT.
\item[\texttt{INT4 mObsCoh}] The number of coherent timescale in observation time.  In other words, this quantifies how many DeFTs will be produced.
\item[\texttt{INT4 ifMin}] The index of the minimum frequency of the SFT frequency band.
\item[\texttt{INT4 iCoh}] The index of the DeFT being computed (0 $\leq$ iCoh \verb@<@ mObsCoh).
\item[\texttt{REAL8 *skyConst}] The array of sky constants.
\item[\texttt{REAL8 *spinDwn}] The set of template spindown parameters.
\item[\texttt{AMCoeffs *amcoe}] The values of the function $a$ and $b$, plus their scalar products.
\end{description}

\begin{verbatim}
struct SkyPos
\end{verbatim}
\index{\texttt{SkyPos}}
\noindent This structure contains the basic parameters of a sky position.  These are:
\begin{description}
\item[\texttt{REAL8 alpha}]  The inclination, in radians.
\item[\texttt{REAL8 delta}]  The declination, in radians.
\item[\texttt{CHAR type}]  Coordinate system used.
\end{description}

\begin{verbatim}
struct Spindown
\end{verbatim}
\index{\texttt{Spindown}}
\noindent This structure contains the values of the spindown parameters.  Included are
\begin{description}
\item[\texttt{INT4 m}]  The maximum order of spindown parameter employed, or the number of spindown parameters per spindown parameter set.
\item[\texttt{REAL8 *spParams}]  An array containing the values.
\end{description}

\begin{verbatim}
struct ParameterSet
\end{verbatim}
\index{\texttt{ParameterSet}}
\noindent This is a structure which contains all the parameters that describe a source.  Included are
\begin{description}
\item[\texttt{SkyPos *skyP}]
\item[\texttt{Spindown *spind}]
\end{description}

\begin{verbatim}
struct FFT
\end{verbatim}
\index{\texttt{FFT}}
\noindent This structure is used to hold all Fourier domain data.  Note that it is different from the LAL definition of a \verb@SequenceOfFrequencySeries@ variable in that the \verb@FFT@ structure allows one to record information specific to each individual FFT; the LAL structure does not allow this functionality at present.  The fields contained within are
\begin{description}
\item[\texttt{COMPLEX8FrequencySeries *fft}]  Holds the data (within \verb@fft->data->data@) and other parameters relevant to each transform.
\item[\texttt{ParameterSet par}]  Holds source parameters that might be relevant to the data in the previous field.  For example, if the \verb@*fft@ structure above contains DeFTs, this field will contain the demodulation parameters.
\end{description}

\begin{verbatim}
struct DeFTPeriodogram
\end{verbatim}
\index{\texttt{DeFTPeriodogram}}
\noindent  This structure contains the statistic we wanted to compute: \ref{e1}, and information regarding its parameters.
\begin{description}
\item[\texttt{REAL8FrequencySeries *fft}]  The periodogram.
\item[\texttt{ParameterSer par}]  The parameter set used in the demodulation.
\end{description}

 </lalLaTeX> */
 

/* PARAMETERS */
typedef struct DemodParTag{
  INT4		if0Max;		/* Index of maximum freq. of search on this proc. */
  INT4		if0Min;		/* Index of minimum freq. of search on this proc. */
  INT4		mCohSFT;	/* This is G_Mcoh_sft */
  INT4 		mObsCoh;	/* This is G_Mobs_sft */
  INT4 		ifMin;		/* Index of minimum freq. of search */
  INT4 		iCoh;	        /* Index of coherent timescale */
  INT4		spinDwnOrder;	/* Maximum order of spdwn parameter */
  REAL8		*skyConst;	/* Constants computed in ComputeSky.c */
  REAL8		*spinDwn;	/* Spindown parameter set */
  AMCoeffs      *amcoe;
}DemodPar;


typedef struct SkyPosTag 
{
  REAL8 alpha;
  REAL8 delta;
  CHAR type;
}SkyPos;

typedef struct SpindownTag 
{
	INT4 m;
	REAL8 *spParams;
}Spindown;

typedef struct ParameterSetTag 
{
	SkyPos *skyP;
	Spindown *spind;
}ParameterSet;

/*This structure will hold a single FFT*/
typedef struct FFTTag 
{
  COMPLEX8FrequencySeries *fft;   
  ParameterSet par;
} FFT;

typedef struct DeFTPeriodogramTag
{
  REAL8FrequencySeries *fft;
  ParameterSet par;
} DeFTPeriodogram;


/* <lalLaTeX>
\input{LALDemodHV}
</lalLaTeX> */

/* <lalLaTeX>
\newpage\input{LALDemodC}
</lalLaTeX> */

void LALDemod (LALStatus 	        *stat, 
			DeFTPeriodogram *xHatCoh, 
			FFT 		**input, 
			DemodPar 	*params);

void times(REAL8 ,
		 INT4, 
		 LIGOTimeGPS *, 
		 INT4 );
	
/* <lalLaTeX>
\newpage\input{LALDemodTestC}
</lalLaTeX> */

#ifdef __cplusplus
}
#endif
#endif /* _LALDEMOD_H */
