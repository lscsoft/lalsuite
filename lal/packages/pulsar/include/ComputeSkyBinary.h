/************************************ <lalVerbatim file="ComputeSkyBinaryHV">
Author: Messenger, C.J., Berukoff, S.J., Papa, M.A. 
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX> 
\section{Header \texttt{ComputeSkyBinary.h}}
\label{s:ComputeSkyBinary.h}
Computes phase coefficients necessary for a correct demodulation for a source in a binary system.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/ComputeSkyBinary.h>
\end{verbatim}

\noindent  The methods employed here follow very closely those used within \verb@ComputeSky()@.  
Note that at present this code simply corrects for the Doppler modulation present in a polynomial 
frequency function for signals from sources in elliptical orbits.  It does not account for general 
relativistic effects.

At the risk of repeating existing documentation, but in the interests of clarity much of the 
following can also be found in the \verb@ComputeSky()@ documentation.  Recall that a demodulated 
Fourier Transform (DeFT) is given by 

\begin{equation}
\hat{x}_b({\vec{\lambda}})=
\sum_{\alpha =0}^{M-1}\sum_{k=0}^{N-1}\tilde{x}_{\alpha k}\left[\frac{1}{N}\sum_{j=0}^{N-1}e^{-2\pi i(\Phi_{\alpha jb}(\vec{\lambda})-\frac{jk}{N})}\right]
\label{demod_FT}
\end{equation}

The index $b$ defines the DeFT frequency bin, the index $\alpha$ loops through
the SFTs that build the DeFT, $k$ runs on all the SFT frequency bins, and $j$
is a time index that runs on each SFT.  As shown in section
\ref{s:LALDemod.h}, the next step in the development of the demodulation
technique involves Taylor expanding the phase model about the temporal
midpoint of each short segment of data, while retaining only first order
terms.  At this point it is neccessary to clearly define some quantities. 
Times as defined at the chosen detector are denoted by $t$, times defined at the 
solar system barycenter (SSB) are denoted by $T$, and the retarded time measured at an inertial
reference point (chosen as the SSB) at a distance from the source are denote by $t^{\prime}$.

The Taylor expansion of $\Phi (t)$ about the temporal midpoint
$t_{\alpha,1/2}$ is

\begin{equation}
\Phi_{\alpha}(t) = \Phi(t_{\alpha,1/2})+\left[t-t_{\alpha,1/2}\right]\frac{d\Phi}{dt}(t_{\alpha,1/2})\label{taylor2} \\
\end{equation}

For each value of $\alpha$, this expression consists of either constant or linear terms in time.  
With the particular time discretization chosen in this code, $t=t_{0}+(N\alpha+j)\ T_{obs}/NM$, we have

\begin{equation}
\label{time}
\left[t-t_{\alpha,1/2}\right]=\frac{\ T_{obs}}{M}\left(\frac{j}{N}-\frac{1}{2}\right)=\mathcal{T}_{s}\left(\frac{j}{N}-\frac{1}{2}\right),
\end{equation}

where $\mathcal{T}_{s}$ is the short time baseline of the $M$ short FTs.  On
the other hand, the phase can also be expressed as a function of SSB time $T$
(i.e. the time at the solar system barycenter).  We will assume the source to
be at rest in this reference frame.  If we now adopt the notation $\Delta
t^{\prime}_{\alpha}\equiv\left[t^{\prime}(T(t_{\alpha,1/2}))-
t^{\prime}(T(t_{0}))\right]$ and $\dot{t^{\prime}}_{\alpha}\equiv dt^{\prime}/dt({\alpha,1/2})$,
the phase terms described in Eq .\ref{taylor2} become (neglecting constants)

\begin{eqnarray}
\Phi(t_{\alpha,1/2})   & = & f_{0}(\Delta t^{\prime}_{\alpha})+\frac{1}{2}f_{1}(\Delta t^{\prime}_{\alpha})^{2}
+\frac{1}{3}f_{2}(\Delta t^{\prime}_{\alpha})^{3}+\frac{1}{4}f_{3}(\Delta t^{\prime}_{\alpha})^{4}+\frac{1}{5}f_{4}(\Delta t^{\prime}_{\alpha})^{5}
+\frac{1}{6}f_{5}(\Delta t^{\prime}_{\alpha})^{6}, \nonumber\label{phi} \\
                                         &   & \\
\frac{d\Phi}{dt}(t_{\alpha,1/2})         & = & \dot{t^{\prime}}_{\alpha}\left(f_{0}+ f_{1}(\Delta t^{\prime}_{\alpha})
+f_{2}(\Delta t^{\prime}_{\alpha})^{2}+f_{3}(\Delta t^{\prime}_{\alpha})^{3}
+f_{4}(\Delta t^{\prime}_{\alpha})^{4}+f_{5}(\Delta t^{\prime}_{\alpha})^{5}\right). \label{dphi}
\end{eqnarray}

Note that the polynomial phase function is expressed as a function of the retarded time $t^{\prime}$ and subsequently
the intrinsic frequency and it`s derivitives ($f_{i}$) are defined at the chosen inertial reference frame (SSB). 

In order to calculate, for each value of $\alpha$, the quantities $\dot{t^{\prime}}_{\alpha}$ and
$\Delta t^{\prime}_{\alpha}$, we must now look at the binary system in more detail.  At present
we reference section \ref{s:GenerateSpinOrbitCW.h} where the definition of all the following orbital variables
can be found.  For a given set of orbital input parameters we obtain the eccentric anomoly $E$ by
numerically solving 

\begin{eqnarray}\label{bintime}
  T(t_{\alpha})-T_{p}&=&\frac{P}{2\pi}\left[E+\left(p\sin{E}+q\left(\cos{E}-1\right)\right)\right].
\end{eqnarray}

where the quantities $p$ and $q$, dependent only on the orbital parameters
of the source system, are given by

\begin{eqnarray}
\label{pq}
p&=&\frac{2\pi}{P}\frac{a}{c}\sin{i}\sqrt{1-e^{2}}\cos{\omega}-e \nonumber \\
q&=&\frac{2\pi}{P}\frac{a}{c}\sin{i}\sin{\omega}.
\end{eqnarray}

$T(t_{\alpha})$ is returned via a call to \verb@LALBarycenter@ and $a\sin{i},P,T_{p},\omega,e$ 
are the projected semi-major axis (projected along the line of sight), the orbital period, the 
time of observed periapse passage as measured 
in the SSB, the argument of periapse, and the orbital eccentricity respectively.  Having defined $E$ 
(where the source is in it`s orbit) at a given detector time $t_{\alpha}$ we can calculate the derivitive 
of the retarded source time $\prime{t}$ with respect to the SSB time $T$.  This is given by 

\begin{equation}\label{binderiv}
 \frac{dt^{\prime}}{dT}=\frac{[1-e\cos{E}]}{\left[1+pcos{E}-q\sin{E}\right]}.
\end{equation}

The quantity $\dot{t^{\prime}}_{\alpha}$ can now be expressed as

\begin{equation}\label{binderivtwo}
\dot{t^{\prime}}_{\alpha}=\frac{dT}{dt}\frac{dt^{\prime}}{dT},
\end{equation}

where $dT/dt$ is returned via a call to \verb@LALBarycenter()@.

We can now rewrite Eq. \ref{taylor2} and by grouping together the terms in $j$ (linear
in $t$) in order to save computations, we have

\begin{equation}
\Phi_{\alpha}(t)=\sum_{s=0}^{n_{spin}}f_{s}A_{s\alpha}+\frac{j}{N}\sum_{s=0}^{n_{spin}}f_{s}B_{s\alpha},
\label{phasecalc}
\end{equation}

where $n_{spin}$ is the maximum order of spindown parameter.

Thus, for a given sky position and set of orbital parameters, the quantities $\dot{t^{\prime}}_{\alpha}$ and
$\Delta t^{\prime}_{\alpha}$ are calculated only once, just as in \verb@ComputeSky()@.  The analytical constants 
defined in Eq \ref{phasecalc} now become

\begin{equation}
A_{s \alpha}=\frac{1}{s+1}\Delta (t^{\prime}_{\alpha})^{s+1}-\frac{1}{2}\mathcal{T}_{SFT}\dot{t^{\prime}}_{\alpha}\Delta (t^{\prime}_{\alpha})^{s}
\end{equation}
\begin{equation}
B_{s \alpha}=\mathcal{T}_{SFT}\dot{t^{\prime}}_{\alpha}\Delta (t^{\prime}_{\alpha})^{s}.
\end{equation}

It is these constants that form the input to the function \verb@LALDemod()@.

</lalLaTeX> */

#ifndef _COMPUTESKYBINARY_H
#define _COMPUTESKYBINARY_H
 
#include <lal/LALStdlib.h>
#include <lal/LALBarycenter.h>
#include <lal/Date.h>
 
#ifdef __cplusplus
extern "C" {
#endif

  NRCSID( COMPUTESKYBINARYH, "$Id$" );
  
  /* Author-defined error codes */
  /* <lalLaTeX>  
     \subsection*{Error conditions}
     \vspace{0.1in}
     \input{ComputeSkyBinaryHErrorTable}
     
     </lalLaTeX> */
  
/* <lalErrTable file="ComputeSkyBinaryHErrorTable"> */
#define COMPUTESKYBINARYH_ENULL 1
#define COMPUTESKYBINARYH_ENNUL 2
#define COMPUTESKYBINARYH_ERANG 3
#define COMPUTESKYBINARYH_ENEGA 4
#define COMPUTESKYBINARYH_MSGENULL "Null Pointer"
#define COMPUTESKYBINARYH_MSGENNUL "Non-Null Pointer"
#define COMPUTESKYBINARYH_MSGERANG "Input parameter out of range"
#define COMPUTESKYBINARYH_MSGENEGA "Bad Negative Value"
/* </lalErrTable>  */

/* <lalLaTeX>
\subsection*{Structures}

\begin{verbatim}
struct CSBParams
\end{verbatim}
\index{\texttt{CSBParams}}

\noindent This structure contains the parameters for the \verb@ComputeSkyBinary()@ routine.  The parameters are:

\begin{description}
\item[\texttt{INT8 spinDwnOrder}] The maximal number of spindown parameters per spindown parameter set.
\item[\texttt{INT8 mObsSFT}] The number of SFTs in the observation time.
\item[\texttt{REAL8 tSFT}] The timescale of one SFT.
\item[\texttt{LIGOTimeGPS *tGPS}] An array containing the GPS times of the first datum from each SFT.
\item[\texttt{REAL8 *skyPos}] The array containing the sky patch coordinates.
\item[\texttt{REAL8 SemimajorAxis}] The projected speed-of-light-normalised semi-major axis of the orbit (in seconds).
\item[\texttt{REAL8 OrbitalPeriod}] The orbital period (in seconds).
\item[\texttt{REAL8 ArgPeriapse}] The argument of the periapse (in radians).
\item[\texttt{REAL8 OrbitalEccentricity}] The orbital eccentricity.
\item[\texttt{LIGOTimeGPS TperiapseSSB}] A time of observed periapse passage as defined in the SSB.
\end{description}

</lalLaTeX> */

/* The following quantity represents the required accuracy for the timing of the binary orbit in seconds */
#define ACC 1e-9

typedef struct
tagCSBParams
{
	INT8		spinDwnOrder;	/* max spindown parameter order */
	INT8		mObsSFT;	/* number of coherent timescales */
	REAL8		tSFT;		/* timescale of SFT */
	LIGOTimeGPS	*tGPS;		/* GPS time of 1st data sample of each SFT */
	REAL8 		*skyPos; 	/* array of sky positions */
        REAL8 		SemiMajorAxis;  /* orbital radius of binary (in sec) */
        REAL8           OrbitalPeriod;         /* Period of binary (in sec) */
        REAL8           OrbitalEccentricity; /* Orbital eccentricy */
        REAL8           ArgPeriapse;    /* Argument of Periapse */
        LIGOTimeGPS     TperiapseSSB;   /* Instance of periapse passage measured in the SSB frame */
	BarycenterInput *baryinput;	
	EmissionTime *emit;
	EarthState *earth;
	EphemerisData *edat;
}
CSBParams;

/* <lalLaTeX>
\newpage\input{ComputeSkyBinaryHV}
\newpage\input{ComputeSkyBinaryC}
</lalLaTeX> */

void ComputeSkyBinary (LALStatus *status, 
			REAL8 		*skyConst, 
			INT8 		iSkyCoh,
			CSBParams 	*params);
		

#ifdef __cplusplus
}
#endif

#endif /* _COMPUTESKYBINARY_H */
