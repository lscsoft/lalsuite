/*
*  Copyright (C) 2007 Jolien Creighton
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
 * File Name: CLR.h
 *
 * Author: Sintes, A. M.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * CLR.h
 *
 * SYNOPSIS
 * #include "CLR.h"
 *
 * DESCRIPTION
 *
 * DIAGNOSTICS
 *
 *-----------------------------------------------------------------------
 */




/************************************ <lalVerbatim file="CLRHV">
Author: Sintes, A. M.
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{CLR.h}}
\label{s:CLR.h}

Provides routines for finding line harmonics,
generating a reference interference signal, and removing
all the interference harmonics.



\subsection*{Synopsis}


\begin{verbatim}
#include "CLR.h"
\end{verbatim}


\noindent The principal of  {\sc clr} is the following:

We assume that the interference has the form
\begin{equation}
y(t)=\sum_n a_n m(t)^n + \left( a_n m(t)^n\right)^* \ ,
\label{e3}
\end{equation}
where $a_n$ are
complex amplitudes and  $m(t)$ is a nearly monochromatic function
 near
a frequency $f_0$.
The idea is to
 use the information in the different harmonics of the interference
to construct a function $M(t)$ that is  as close a replica as possible of
$m(t)$ and then construct a function close to $y(t)$ which
is subtracted from the output  of the system cancelling the interference.
The key is that real gravitational wave signals will not be present
with multiple harmonics and that $M(t)$ is constructed from many
frequency bands with independent noise. Hence, {\sc clr} will little
affect the statistics of the noise in any one band and
any gravitational wave signal masked by the interference can
be recovered without any disturbance.


We assume that the data produced by the system
is just the sum of the interference plus  noise
\begin{equation}
x(t)=y(t)+n(t) \ ,
\label{e6}
\end{equation}
where $y(t)$ is given by Eq.~(\ref{e3}) and the noise $n(t)$ in the
detector  is a zero-mean stationary
stochastic process.
The procedure consists in   defining a set of functions $\tilde z_k(\nu)$
in the frequency domain as
\begin{equation}
\tilde z_k(\nu)\equiv \left\{
\begin{array}{cc}
\tilde x(\nu) & \nu_{ik}<\nu <\nu_{fk}\\
0 & \mbox{elsewhere}\ ,
\end{array}
\right.
\label{e8}
\end{equation}
where  $(\nu_{ik}, \nu_{fk})$ correspond to the upper and lower frequency
limits of the harmonics of the interference
and $k$ denotes the harmonic considered.
These functions are equivalent to
\begin{equation}
\tilde z_k(\nu)= a_k \widetilde{m^k}(\nu) +\tilde n_k(\nu) \ ,
\end{equation}
where $ \tilde n_k(\nu)$ is the noise in the frequency band of the
harmonic considered.  Their inverse  Fourier transforms yield
\begin{equation}
z_k(t)=a_k m(t)^k +n_k(t) \ .
\end{equation}
Since  $m(t)$ is supposed to be a narrow-band function near a frequency $f_0$,
each $z_k(t)$ is a  narrow-band function near $kf_0$.
Then, we  define
\begin{equation}
B_k(t)\equiv \left[ z_k(t)\right]^{1/k}\ ,\label{e10a}
\end{equation}
that can be rewritten as
\begin{equation}
B_k(t)= (a_k)^{1/k}m(t) \beta_k(t) \ , \qquad
\beta_k(t)=\left[ 1+ {n_k(t) \over a_k m(t)^k}\right]^{1/k} \ .
\label{e10}
\end{equation}
All these  functions, $\{B_k(t)\}$, are almost monochromatic around the
fundamental frequency, $f_0$, but they differ basically by a certain
complex amplitude. These factors, $\Gamma_k$, can easily be  calculated,
and  we can construct a set of functions   $\{b_k(t)\}$
\begin{equation}
 b_k(t)=\Gamma_k B_k(t)\ ,
 \end{equation}
such that, they all have the same mean value. Then,   $M(t)$ can be
 constructed  as a function of all $\{b_k(t)\}$
 in such a way
that it has the same mean and minimum variance.
If
$M(t)$ is linear with $\{b_k(t)\}$, the statistically the best is
\begin{equation}
 M(t)=\left(\sum_k {b_k(t) \over {\rm Var}[\beta_k(t)]} \right) {\Big
 { /}}
\left( \sum_k {1 \over {\rm Var}[\beta_k(t)]}\right) \ ,
\end{equation}
where
\begin{equation}
{\rm Var}[\beta_k(t)]= {\langle n_k(t) n_k(t)^*\rangle\over  k^2
\vert a_k m(t)^k\vert^2}+ \mbox{corrections} \ .
\end{equation}
In practice,
we  approximate
\begin{equation}
\vert a_k m(t)^k\vert^2 \approx \vert z_k(t)\vert^2 \ ,
\end{equation}
and we assume  stationary noise. Therefore,
\begin{equation}
\langle n_k(t) n_k(t)^*\rangle= \int_{\nu_{ik}}^{\nu_{fk}} S(\nu) d\nu \ ,
\end{equation}
where $S(\nu)$ is the power spectral density of the noise.

Finally, it only remains to determine the amplitude of the different
harmonics, which can be obtained applying a least square method.


</lalLaTeX> */


/*
 * 2. include-loop protection (see below). Note the naming convention!
 */

#ifndef _CLR_H
#define _CLR_H

/*
 * 3. Includes. This header may include others; if so, they go immediately
 *    after include-loop protection. Includes should appear in the following
 *    order:
 *    a. Standard library includes
 *    b. LDAS includes
 *    c. LAL includes
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>
#include <lal/ComplexFFT.h>

/*
 *  #include "LALRCSID.h"
 *	not needed, it is included in "LALConstants.h"
 */

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (CLRH, "$Id: CLR.h");

/*
 * 5. Macros. But, note that macros are deprecated.
 */

/*
 * 8. Structure, enum, union, etc., typdefs.
 */


/* <lalLaTeX>
\subsection*{Error conditions}
\vspace{0.1in}
\input{CLRHErrorTable}

</lalLaTeX> */



/* <lalErrTable file="CLRHErrorTable"> */

#define CLRH_ENULL 1
#define CLRH_ESIZE 2
#define CLRH_ESZMM 4
#define CLRH_EINT  6
#define CLRH_ESAME 8
#define CLRH_EFREQ 10

#define CLRH_MSGENULL "Null pointer"
#define CLRH_MSGESIZE "Invalid input size"
#define CLRH_MSGESZMM "Size mismatch"
#define CLRH_MSGEINT  "Invalid interval"
#define CLRH_MSGESAME "Input/Output data vectors are the same"
#define CLRH_MSGEFREQ "Invalid frequency"

/* </lalErrTable>  */



/* <lalLaTeX>

\subsection*{Structures}

\begin{verbatim}
struct REAL4TVectorCLR
\end{verbatim}
\idx[Type]{REAL4TVectorCLR}

\noindent This structure stores the time domain data and other
information needed in order to remove all the line interference
harmonics.  The fields are:

\begin{description}
\item[\texttt{UINT4 length}]  The number of elements in data.
\item[\texttt{REAL4  *data}] The (real) time domain data.
\item[\texttt{REAL4  deltaT}] The sample spacing in time (in seconds).
\item[\texttt{REAL4  fLine}] The interference fundamental frequency $f_0$
    (in Hz), e.g., 60 Hz.
\end{description}



\begin{verbatim}
struct REAL4FVectorCLR
\end{verbatim}
\idx[Type]{REAL4FVectorCLR}

\noindent This structure stores the spectrum, $\vert\tilde x(\nu)\vert^2$,
 and other information
needed to find the location of the line harmonics.
  The fields are:

\begin{description}
\item[\texttt{UINT4 length}]  The number of elements in data.
\item[\texttt{REAL4   *data}] The (real)  frequency domain data, e.g.,
                              $\vert\tilde x(\nu)\vert^2$.
\item[\texttt{REAL4   deltaF}] The $\Delta$F offset between samples (in Hz).
\item[\texttt{REAL4   fLine}] The interference fundamental frequency $f_0$
    (in Hz), e.g., 60 Hz.
\end{description}


</lalLaTeX> */


  /* Available CLR Time Vector  */

typedef struct tagREAL4TVectorCLR{
  UINT4   length;  /* number of elements in data */
  REAL4   *data;   /* time domain data */
  REAL4   deltaT;  /* sample spacing in time (in seconds) */
  REAL4   fLine;   /* interference fundamental frequency (e.g., 50 Hz) */
} REAL4TVectorCLR;


  /* Available CLR Frequency Vector  */

typedef struct tagREAL4FVectorCLR{
  UINT4   length;  /* number of elements in data */
  REAL4   *data;   /* frequency domain data */
  REAL4   deltaF;  /* delta F offset between samples (in Hz)*/
  REAL4   fLine;   /* interference fundamental frequency (e.g., 50 Hz) */
} REAL4FVectorCLR;



/* <lalLaTeX>
\vfill{\footnotesize\input{CLRHV}}
</lalLaTeX> */



  /*
   * 9. Functions Declarations (i.e., prototypes).
   */



/* <lalLaTeX>
\newpage\input{HarmonicFinderC}
\newpage\input{RefInterferenceC}
\newpage\input{CleanAllC}
</lalLaTeX> */

void LALRefInterference ( LALStatus          *status,
			  COMPLEX8Vector     *out, /*  M(t), size n */
			  COMPLEX8Vector     *in,  /* x(nu), size n/2+1 */
			  INT4Vector         *par
			  );

void LALCleanAll ( LALStatus        *status,
		   REAL4Vector      *out,  /* clean data */
		   COMPLEX8Vector   *in2,  /* M(t), ref. interference */
		   REAL4TVectorCLR  *in1   /* x(t), data + information */
		 );


void LALHarmonicFinder (
         LALStatus          *status,
         INT4Vector         *out,   /* harmonic index, bin location (3*l) */
         REAL4FVectorCLR    *in2,   /* |x(f)|^2, data + inform.  */
         INT4Vector         *in1    /* the harmonic index (l) */
		 );



/* Test program. */

/* <lalLaTeX>
\newpage\input{CLRTestC}
</lalLaTeX> */

#ifdef  __cplusplus
}
#endif

#endif /* _CLR_H */
