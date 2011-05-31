/*
*  Copyright (C) 2007 Sukanta Bose, Chad Hanna, Darren Woods, Diego Fazi, Drew Keppel, Duncan Brown, Eirini Messaritaki, Gareth Jones, Jolien Creighton, Patrick Brady, Anand Sengupta, Stephen Fairhurst, Craig Robinson , Sean Seader, Thomas Cokelaer
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
 * File Name: FindChirp.h
 *
 * Author: Allen, B., Brown, D. A. and Creighton, J. D. E.
 *
 *-----------------------------------------------------------------------
 */

/**
 * \defgroup FindChirp_h FindChirp_h
 * \ingroup CBC_findchirp
 */

/**
\author Allen, B., Brown, D. A. and Creighton, J. D. E.
\file
\ingroup FindChirp_h

\brief This header provides core prototypes, structures and functions to
filter interferometer data for binary inspiral chirps.

\heading{Synopsis}
\code
#include <lal/FindChirp.h>
\endcode

Each function in findchirp falls into one of four classes:
<ol>
<li> Generate management functions which are independent of the type of
filtering implemented. The prototypes for these functions are provided by
this header file.</li>

<li> Functions for filtering data for time domain and frequency domain
templates with an unknown amplitude and phase. These are the functions
that implement matched filtering for time domain templates (TaylorT1,
TaylorT2, TaylorT2, PadeT1, EOB and GeneratePPN) and matched filtering
for post-Newtonian frequency domain templates (FindChirpSP). The main
filter function <tt>FindChirpFilterSegment()</tt> is prototyped in
this header file. The template generation and data conditioning functions
are prototyped in \ref FindChirpSP.h and \ref FindChirpTD.h.
Full documentation of the filtering algorithm used can be found in the
documentation of the module \ref FindChirpFilter.c.</li>

<li> Functions to filter interferometer data for using the frequency
domain non-spinning black hole detection template family known as BCV.
These functions are protoyped by the header \ref FindChirpBCV.h
which contains documentation of the algorithms used.</li>

<li> Functions to filter interferometer data for using the frequency domain
spinning black hole detection template family known as BCVSpin.  These
functions are protoyped by the header \ref FindChirpBCVSpin.h which
contains documentation of the algorithms used.</li>
</ol>

The goal of all the filtering functions is to determine if the
(calibrated) output of the interferometer \f$s(t)\f$ contains a gravitational wave
\f$h(t)\f$ in the presence of the detector noise \f$n(t)\f$. When the interferometer
is operating properly
\f{equation}{
s(t) = \left\{ \begin{array}{ll}
n(t) + h(t) & \textrm{signal present},\\
n(t) & \textrm{signal absent}.
\end{array}\right.
\f}
The detection of signals of known form in noise is a classic problem of signal
processing [\ref WeinZub62] and can be answered by the construction of a
<em>detection statistic</em> and a test to see if the statistic is above some
pre-assigned threshold. The construction of the various detection
statistics used for each the three types of search are described in the modules
that implement the search.


*/

#ifndef _FINDCHIRPH_H
#define _FINDCHIRPH_H

#include <lal/LALDatatypes.h>
#include <lal/ComplexFFT.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataInspiralUtils.h>
#include <lal/LALInspiral.h>
#include <lal/LALInspiralBank.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/FindChirpDatatypes.h>
#include <lal/FindChirpChisq.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif


NRCSID (FINDCHIRPH, "$Id$");

/**\name Error Codes */ /*@{*/
#define FINDCHIRPH_ENULL 1
#define FINDCHIRPH_ENNUL 2
#define FINDCHIRPH_EALOC 3
#define FINDCHIRPH_ENUMZ 5
#define FINDCHIRPH_ESEGZ 6
#define FINDCHIRPH_ECHIZ 7
#define FINDCHIRPH_EDTZO 8
#define FINDCHIRPH_ETRNC 10
#define FINDCHIRPH_EFLOW 11
#define FINDCHIRPH_EFREE 12
#define FINDCHIRPH_ERHOT 15
#define FINDCHIRPH_ECHIT 16
#define FINDCHIRPH_ECRUP 17
#define FINDCHIRPH_ESMSM 18
#define FINDCHIRPH_EHETR 19
#define FINDCHIRPH_EDFDT 20
#define FINDCHIRPH_EAPRX 21
#define FINDCHIRPH_EUAPX 22
#define FINDCHIRPH_ECHTZ 23
#define FINDCHIRPH_EMASS 24
#define FINDCHIRPH_EWVFM 25
#define FINDCHIRPH_EBCVC 25
#define FINDCHIRPH_EMAPX 26
#define FINDCHIRPH_EPTFW 27
#define FINDCHIRPH_EIGEN 28
#define FINDCHIRPH_EIMRW 29

#define FINDCHIRPH_MSGENULL "Null pointer"
#define FINDCHIRPH_MSGENNUL "Non-null pointer"
#define FINDCHIRPH_MSGEALOC "Memory allocation error"
#define FINDCHIRPH_MSGENUMZ "Invalid number of points in segment"
#define FINDCHIRPH_MSGESEGZ "Invalid number of segments"
#define FINDCHIRPH_MSGECHIZ "Invalid number of chi squared bins"
#define FINDCHIRPH_MSGEDTZO "deltaT is zero or negative"
#define FINDCHIRPH_MSGETRNC "Duration of inverse spectrum in time domain is negative"
#define FINDCHIRPH_MSGEFLOW "Inverse spectrum low frequency cutoff is negative"
#define FINDCHIRPH_MSGEFREE "Error freeing memory"
#define FINDCHIRPH_MSGERHOT "Rhosq threshold is negative"
#define FINDCHIRPH_MSGECHIT "Chisq threshold is negative"
#define FINDCHIRPH_MSGECRUP "Chirp length or invSpecTrunc too long for length of data segment"
#define FINDCHIRPH_MSGESMSM "Size mismatch between vectors"
#define FINDCHIRPH_MSGEHETR "Attempting to simulate heterodyned GW"
#define FINDCHIRPH_MSGEDFDT "Waveform sampling interval is too large"
#define FINDCHIRPH_MSGEAPRX "Incorrect waveform approximant"
#define FINDCHIRPH_MSGEUAPX "Unknown waveform approximant"
#define FINDCHIRPH_MSGECHTZ "Length of chirp is zero or negative"
#define FINDCHIRPH_MSGEMASS "Invalid mass parameters for template generation"
#define FINDCHIRPH_MSGEWVFM "Unknown injection waveform"
#define FINDCHIRPH_MSGEBCVC "BCVC code: thetav not in [-pi, pi]."
#define FINDCHIRPH_MSGEMAPX "Mismatch in waveform approximant"
#define FINDCHIRPH_MSGEPTFW "Error generating PTF waveform"
#define FINDCHIRPH_MSGEIGEN "Error computing eigenvalues"
#define FINDCHIRPH_MSGEIMRW "Error computing IMR waveform"
/*@}*/


/** This structure provides the essential information for the
filter initialisation and memory allocation functions used by findchirp.

<dl>
<dt><tt>UINT4 numSegments</tt></dt><dd> The number of data segments in the input
\c DataSegmentVector and a the \c FindChirpSegmentVector.</dd>

<dt><tt>UINT4 numPoints</tt></dt><dd> The number of discrete data points \f$N\f$ in each
data segment.</dd>

<dt><tt>UINT4 ovrlap</tt></dt><dd> The number of sample points by which each
data segment overlaps.</dd>

<dt><tt>UINT4 numChisqBins</tt></dt><dd> The number of bins \f$p\f$ used to contruct the
\f$\chi^2\f$ veto.</dd>

<dt><tt>BOOLEAN createRhosqVec</tt></dt><dd> Flag that controls whether or not the
filter function should store the output of the matched filter, \f$\rho^2(t)\f$, as
well as the events. Memory is allocated for this vector if the flag is set to
1.</dd>

<dt><tt>BOOLEAN createCVec</tt></dt><dd> Flag that controls whether or not the
filter function should store the complex filter output \f$x(t) + i y(t)\f$ needed
by the coherent inspiral code.  Memory is allocated for this vector if the
flag is set to 1.</dd>

<dt><tt>Approximant approximant</tt></dt><dd> Initialize the findchirp routines
to fiter with templates of type \c approximant. Valid approximants are
TaylorT1, TaylorT2, TaylorT3, PadeT1, EOB, FindChirpSP, BCV and BCVSpin.
</dd>
</dl>

*/
typedef struct
tagFindChirpInitParams
{
  UINT4                         numSegments;
  UINT4                         numPoints;
  UINT4                         ovrlap;
  UINT4                         numChisqBins;
  BOOLEAN                       createRhosqVec;
  BOOLEAN                       createCVec;
  Approximant                   approximant;
  LALPNOrder                    order;
}
FindChirpInitParams;


/** This structure contains the parameters needed to call the data
conditioning functions <tt>FindChirpSPData()</tt>, <tt>FindChirpTDData()</tt>,
<tt>FindChirpBCVData()</tt> or <tt>FindChirpBCVSpinData()</tt>. It should be
initialized by <tt>FindChirpDataInit()</tt> and destroyed by
<tt>FindChirpDataFinalize()</tt>.

<dl>
<dt><tt>REAL4Vector *ampVec</tt></dt><dd> A vector containing the frequency domain
quantity \f$(k/N)^{-7/6}\f$, where \f$k\f$ is the frequency series index and \f$N\f$ is the
number of points in a data segment. NB: for time domain templates, this is set
to unity by the function <tt>FindChirpTDData()</tt>.</dd>

<dt><tt>REAL4Vector *ampVecBCV</tt></dt><dd> A vector containing the frequency domain
quantity \f$(k/N)^{-1/2}\f$, where \f$k\f$ is the frequency series index and \f$N\f$ is the
number of points in a data segment.</dd>

<dt><tt>REAL4Vector *ampVecBCVSpin1</tt></dt><dd> Undocumented spinning BCV
amplitude vector.</dd>

<dt><tt>REAL4Vector *ampVecBCVSpin2</tt></dt><dd> Undocumented spinning BCV
amplitude vector.</dd>

<dt><tt>REAL4Vector *fwdPlan</tt></dt><dd> An FFTW plan used to transform the
time domain interferometer data \f$v(t_j)\f$ into its DFT \f$\tilde{v}_k\f$.</dd>

<dt><tt>REAL4Vector *invPlan</tt></dt><dd> An FFTW plan used to transform the
dimensionless frequency domain interferometer strain \f$\tilde{w}_k\f$ into
the quantity \f$N w(t_j)\f$ to allow time domain trunction of the inverse
power spectrum.</dd>

<dt><tt>REAL4Vector *wVec</tt></dt><dd> A vector used as workspace when truncating
the imverse power spectrum in the time domain.</dd>

<dt><tt>COMPLEX8Vector *wtildeVec</tt></dt><dd> A vector which on exit from
the data conditioning function contains the inverse of the strain one sided
power spectral density, after trunction in the time domain, <em>for the last
data segment conditioned.</em> Typically all the data segments are conditioned
using the same power spectrum, so this quantity is identical for all data
segments. It contains:
\f{equation}{
\tilde{w}_k = {1}/{\ospsd}.
\f}</dd>

<dt><tt>REAL4Vector *tmpltPowerVec</tt></dt><dd> A vector which on exit from
<tt>FindChirpSPData()</tt> or from <tt>FindChirpBCVData()</tt>
contains the quantity
\f{equation}{
\mathtt{tmpltPower[k]} = \frac{f^{-7/3}}{\ospsd}
\f}</dd>

<dt><tt>REAL4Vector *tmpltPowerVecBCV</tt></dt><dd> A vector which on exit from
<tt>FindChirpBCVData()</tt>
contains the quantity
\f{equation}{
\mathtt{tmpltPowerBCV[k]} = \frac{f^{-1}}{\ospsd}
\f}</dd>

<dt><tt>REAL4 fLow</tt></dt><dd> The frequency domain low frequency cutoff
\f$f_\mathrm{low}\f$. All frequency domain data is set to zero below this
frequency.</dd>

<dt><tt>REAL4 dynRange</tt></dt><dd> A dynamic range factor \f$d\f$ which cancels from
the filter output.  This allows quantities to be stored in the range of
\c REAL4 rather than \c REAL8. This must be set to the same value
as \c dynRange in the \c FindChirpTmpltParams. For LIGO data a
value of \f$d = 2^{69}\f$ is appropriate.</dd>

<dt><tt>UINT4 invSpecTrunc</tt></dt><dd> The length to which to truncate the inverse
power spectral density of the data in the time domain. If set to zero, no
truncation is performed.</dd>

<dt><tt>Approximant approximant</tt></dt><dd> Condition the data for templates of
type \c approximant. Valid approximants are TaylorT1, TaylorT2,
TaylorT3, PadeT1, EOB, FindChirpSP, BCV and BCVSpin.
</dd>
</dl>


*/
typedef struct
tagFindChirpDataParams
{
  REAL4Vector                  *ampVec;
  REAL4Vector                  *ampVecBCV;
  REAL8Vector                  *ampVecBCVSpin1;
  REAL8Vector                  *ampVecBCVSpin2;
  RealFFTPlan                  *fwdPlan;
  RealFFTPlan                  *invPlan;
  REAL4Vector                  *wVec;
  COMPLEX8Vector               *wtildeVec;
  REAL4Vector                  *tmpltPowerVec;
  REAL4Vector                  *tmpltPowerVecBCV;
  REAL4                         fLow;
  REAL4                         dynRange;
  UINT4                         invSpecTrunc;
  Approximant                   approximant;
}
FindChirpDataParams;

/** This structure contains the parameters for generation of templates
by the various template generation functions provided in finchirp.

<dl>
<dt><tt>REAL8 deltaT</tt></dt><dd> The sampling interval \f$\Delta t\f$ of the input
data channel.</dd>

<dt><tt>REAL4 fLow</tt></dt><dd> The frequency domain low frequency cutoff
\f$f_\mathrm{low}\f$. All frequency domain data is zero below this frequency.</dd>

<dt><tt>REAL4 dynRange</tt></dt><dd> A dynamic range factor \f$d\f$ which cancels from
the filter output.  This allows quantities to be stored in the range of
\c REAL4 rather than \c REAL8. This must be set to the same value
as \c dynRange in the \c FindChirpDataParams. For LIGO data a
value of \f$d = 2^{69}\f$ is appropriate.</dd>

<dt><tt>REAL4Vector *xfacVec</tt></dt><dd> For frequency domain templates, this is a
vector of length \f$N/2+1\f$ which contains the quantity \f$k^{-1/3}\f$. For time
domain templates, this is a workspace vector of length \f$N\f$ which contains the
time domain template generated by the inspiral package, shifted so that the
end of the template is at the end of the vector. This vector is Fourier
transformed to obtain the quantity findchirp template \f$\tilde{T}_k\f$.</dd>

<dt><tt>REAL4Vector *fwdPlan</tt></dt><dd> For time domain templates, an FFTW plan
used to transform the time domain data stored in \c xfacVec into its DFT
which is stored in the findchirp template.</dd>

<dt><tt>Approximant approximant</tt></dt><dd> Generate templates of type
\c approximant. Valid approximants are TaylorT1, TaylorT2, TaylorT3,
PadeT1, EOB, FindChirpSP, BCV and BCVSpin. For time domain templates the
post-Newtonian order is always two. For stationary phase templates, the
post-Newtonian order is specified by \c order.</dd>

<dt><tt>LALPNOrder order</tt></dt><dd> Specifies the post-Newtonian order of the
templates. Valid pN orders are twoPN, twoPointFivePN, threePN,
threePointFivePN, pseudoFourPN. pseudoFourPN is not the true four PN
correction, but may increase the fitting factor between stationary phase and
numerical relativity waveforms.</dd>

<dt><tt>INT4 reverseChirpBank</tt></dt><dd> Switches a FindChirpSP template bank
to be a reverse chirp template bank if true.</dd>
</dl>

*/
typedef struct
tagFindChirpTmpltParams
{
  REAL8                         deltaT;
  REAL4                         fLow;
  REAL4                         dynRange;
  REAL4Vector                  *xfacVec;
  REAL4VectorSequence          *ACTDVecs;
  REAL4VectorSequence          *PTFQ;
  REAL4Vector                  *PTFphi;
  REAL4Vector                  *PTFomega_2_3;
  REAL4VectorSequence          *PTFe1;
  REAL4VectorSequence          *PTFe2;
  RealFFTPlan                  *fwdPlan;
  Approximant                   approximant;
  LALPNOrder                    order;
  INT4                          reverseChirpBank;
  INT4                          bandPassTmplt;
  InspiralApplyTaper            taperTmplt;
}
FindChirpTmpltParams;

/** This structure contains the possible methods by which
to maximize over a chirp in a data segment.

<dl>
<dt>FindChirpClustering_none</dt><dd>
The decision to do no clustering
of events.</dd>

<dt>FindChirpClustering_tmplt</dt><dd>
Cluster over the length of the data segment.</dd>

<dt>FindChirpClustering_window</dt><dd>
Cluster over a given number of seconds
given by the argument to the flag <tt>--cluster-window</tt>
(required to be less than the length of the data segment).
</dd>
</dl>

*/
typedef enum {
   FindChirpClustering_none,
   FindChirpClustering_tmplt,
   FindChirpClustering_window,
   FindChirpClustering_tmpltwindow
}
FindChirpClustering;


/** This structure provides the parameters for the filter
    output veto.

    <dl>
    <dt><tt>REAL4 rsqvetoWindow</tt></dt><dd> Width of the \f$r^2\f$ veto window in units of
    seconds.</dd>

    <dt><tt>REAL4 rsqvetoThresh</tt></dt><dd> Threshold of the \f$r^2\f$ veto test analogous to the
    \f$r^2\f$ threshold employed in the bns and macho inspiral searches.</dd>
    </dl>
*/
typedef struct
tagFindChirpFilterOutputVetoParams
{
  REAL4          rsqvetoWindow;
  REAL4          rsqvetoThresh;
  REAL4          rsqvetoTimeThresh;
  REAL4          rsqvetoMaxSNR;
  REAL4          rsqvetoCoeff;
  REAL4          rsqvetoPow;
}
FindChirpFilterOutputVetoParams;

/** This structure provides the parameters used by the
    <tt>FindChirpFilterSegment()</tt> function.

<dl>
<dt><tt>REAL8 deltaT</tt></dt><dd> The sampling interval \f$\Delta t\f$.</dd>

<dt><tt>REAL4 rhosqThresh</tt></dt><dd> The signal-to-noise ratio squared threshold
\f$\rho^2_\ast\f$. If the matched filter output exceeds this value, that is
\f$\rho^2(t_j) > \rho^2_\ast\f$, the event processing algorithm is entered and
triggers may be generated (subject to addition vetoes such as the \f$\chi^2\f$
veto). The value of \f$\rho^2_\ast0\f$ must be greater than or equal to zero.</dd>

<dt><tt>REAL4 chisqThresh</tt></dt><dd> The \f$\chi^2\f$ veto threshold on. This
threshold is described in details in the documentation for the \f$\chi^2\f$
veto.</dd>

<dt><tt>REAL4 norm</tt></dt><dd> On exit this contains the normalisation constant
that relates the quantity \f$|q_j|^2\f$ with the signal to noise squared,
\f$\rho^2(t_j)\f$ by
\f{equation}{
\rho^2(t_j) = \textrm{norm} \times \left|q_j\right|^2.
\f}</dd>

<dt><tt>UINT4 maximiseOverChirp</tt></dt><dd> If not zero, use the
maximise over chirp length algorithm to decide which time \f$t_j\f$ should
have an inspiral trigger generated. Otherwise record all points that pass the
\f$\rho^2\f$ and \f$\chi^2\f$ threshold as triggers (this may generate may triggers).</dd>

<dt><tt>Approximant approximant</tt></dt><dd> Filter the data using templates of
type \c approximant. Valid approximants are TaylorT1, TaylorT2,
TaylorT3, PadeT1, EOB, FindChirpSP, BCV and BCVSpin. The value of
\c approximant here must match that in the findchirp data segment and
findchirp template used as input.</dd>

<dt><tt>COMPLEX8Vector *qVec</tt></dt><dd> Pointer to vector of length \f$N\f$ allocated
by <tt>FindChirpFilterInit()</tt> to store the quantity \f$q_j\f$. The pointer
must not be NULL on entry, but the vetor may contain garbage which will be
overwritten with the value of \f$q_j\f$ for the segment filtered on exit.</dd>

<dt><tt>COMPLEX8Vector *qVecBCV</tt></dt><dd> Pointer to the additional vector
required for the BCV templates, allocated by <tt>FindChirpFilterInit()</tt>.</dd>

<dt><tt>COMPLEX8Vector *qVecBCVSpin1</tt></dt><dd> Pointer to the additional vector
required for filtering spinning BCV templates, allocated by
<tt>FindChirpFilterInit()</tt>.</dd>

<dt><tt>COMPLEX8Vector *qVecBCVSpin2</tt></dt><dd> Pointer to the additional vector
required for filtering spinning BCV templates, allocated by
<tt>FindChirpFilterInit()</tt>.</dd>

<dt><tt>COMPLEX8Vector *qtildeVec</tt></dt><dd> Pointer to vector of length \f$N\f$
allocated by <tt>FindChirpFilterInit()</tt> to store the quantity
\f$\tilde{q}_k\f$, given by
\f{equation}{
\tilde{q}_k = \left\{
\begin{array}{ll}
\tilde{F}_k \tilde{T}_k^\ast & \quad 0 < k < \frac{N}{2} \\,
0 & \quad \textrm{otherwise}.
\end{array}
\right.
\f}
The pointer must not be NULL on entry, but the vetor may contain garbage which
will be overwritten with the value of \f$\tilde{q}_k\f$ for the segment filtered
on exit.</dd>

<dt><tt>COMPLEX8Vector *qtildeVecBCV</tt></dt><dd> Pointer to the additional
vector required for filtering BCV templates, allocated by
<tt>FindChirpFilterInit()</tt>.</dd>

<dt><tt>COMPLEX8Vector *qtildeVecBCVSpin1</tt></dt><dd> Pointer to the additional
vector required for filtering spinning BCV templates, allocated by
<tt>FindChirpFilterInit()</tt>.</dd>

<dt><tt>COMPLEX8Vector *qtildeVecBCVSpin2</tt></dt><dd> Pointer to the additional
vector required for filtering spinning BCV templates, allocated by
<tt>FindChirpFilterInit()</tt>.</dd>

<dt><tt>ComplexFFTPlan *invPlan</tt></dt><dd> Pointer to FFTW plan created by
<tt>FindChirpFilterInit()</tt> to transform the quantity \f$\tilde{q}_k\f$ to
\f${q}_j\f$ usimg the inverse DFT. Must not be NULL.</dd>

<dt><tt>REAL4TimeSeries *rhosqVec</tt></dt><dd> Pointer to a time series which
contains a vector of length \f$N\f$. If this is not NULL, the filter output \f$\rho^2(t_j)\f$
is stored in the vector.</dd>

<dt><tt>COMPLEX8Vector *rhosqVec</tt></dt><dd> Pointer to a time series which
contains a vector of length \f$N\f$. If this is not NULL, the complex filter
output \f$\rho(t_j) = x(t_j) + iy(t_j)\f$ is stored in the vector. This
quantity can be used by the coherent filtering code.</dd>

<dt><tt>REAL4Vector *chisqVec</tt></dt><dd> Workspace vector of length \f$N\f$ used to
compute and store \f$\chi^2(t_j)\f$. Must not be NULL if \c numChisqBins is
greater than zero. Contains \f$\chi^2(t_j)\f$ on exit.</dd>

<dt><tt>FindChirpChisqParams *chisqParams</tt></dt><dd> Pointer to parameter
structure for the \f$\chi^2\f$ veto. Must not be NULL if \c numChisqBins is
greater than zero.</dd>

<dt><tt>FindChirpChisqInput *chisqInput</tt></dt><dd> Pointer to input data
structure for the \f$\chi^2\f$ veto.  Must not be NULL if \c numChisqBins is
greater than zero.</dd>

<dt><tt>FindChirpChisqInput *chisqInputBCV</tt></dt><dd> Pointer to input data
structure for the BCV \f$\chi^2\f$ veto. Must not be NULL if the approximant is
BCV and \c numChisqBins is greater than zero.</dd>

<dt><tt>FindChirpFilterOutputVetoParams *filterOutputVetoParams</tt></dt><dd>
Pointer to the parameter structure for the additional signal based veto
function.
</dd>
</dl>
*/
typedef struct
tagFindChirpFilterParams
{
  REAL8                         deltaT;
  REAL4                         clusterWindow;
  REAL4                         rhosqThresh;
  REAL4                         chisqThresh;
  REAL4                         chisqDelta;
  UINT4                         maximiseOverChirp;
  UINT4                         ignoreIndex;
  FindChirpClustering           clusterMethod;
  Approximant                   approximant;
  LALPNOrder                    order;
  COMPLEX8Vector               *qVec;
  COMPLEX8Vector               *qVecBCV;
  COMPLEX8Vector               *qVecBCVSpin1;
  COMPLEX8Vector               *qVecBCVSpin2;
  COMPLEX8Vector               *qtildeVec;
  COMPLEX8Vector               *qtildeVecBCV;
  COMPLEX8Vector               *qtildeVecBCVSpin1;
  COMPLEX8Vector               *qtildeVecBCVSpin2;
  COMPLEX8Vector              **qVecACTD;
  COMPLEX8Vector              **qtildeVecACTD;
  COMPLEX8VectorSequence       *PTFqVec;
  COMPLEX8Vector               *PTFsnrVec;
  REAL4Array                   *PTFA;
  REAL4Array                   *PTFMatrix;
  ComplexFFTPlan               *invPlan;
  REAL4TimeSeries              *rhosqVec;
  COMPLEX8TimeSeries           *cVec;
  REAL4Vector                  *chisqVec;
  FindChirpChisqParams         *chisqParams;
  FindChirpChisqInput          *chisqInput;
  FindChirpChisqInput          *chisqInputBCV;
  FindChirpFilterOutputVetoParams      *filterOutputVetoParams;
}
FindChirpFilterParams;

/*
 *
 * typedefs of input structures used by functions in findchirp
 *
 */


/** This structure groups the input data required for the
    <tt>FindChirpFilterSegment()</tt> function into a single structure.

<dl>
<dt><tt>FindChirpTemplate *fcTmplt</tt></dt><dd> Pointer to the input template
in a form that can be used by <tt>FindChirpFilterSegment()</tt></dd>

<dt><tt>FindChirpSegment *segment</tt></dt><dd> Pointer to the input data segment
in a form that can be used by <tt>FindChirpFilterSegment()</tt></dd>
</dl>

*/

typedef struct
tagFindChirpFilterInput
{
  FindChirpTemplate            *fcTmplt;
  FindChirpSegment             *segment;
}
FindChirpFilterInput;

/** This structure contains data needed for the bank veto.
*/
typedef struct
tagFindChirpBankVetoData
{
  UINT4                   length;
  COMPLEX8Vector        **qtildeVecArray;
  COMPLEX8Vector        **qVecArray;
  FindChirpFilterInput  **fcInputArray;
  COMPLEX8Vector         *ccMat;
  REAL4Vector		 *normMat;
  REAL4Vector		 *spec;
  COMPLEX8Vector         *resp;
  REAL4Vector		 *acorr;
  COMPLEX8Vector	 *workspace;
  REAL4Vector		 *acorrMat;
  REAL4FFTPlan		 *revplan;
  UINT4 		 acorrMatSize;
  UINT4			 autochisqStride;
  REAL4Vector            *timeshift;
  UINT4			 two_sided_auto_chisq;
  UINT4			 time_freq_bank_veto;
}
FindChirpBankVetoData;


/**
<dl>
<dt><tt>Approximant approx</tt></dt><dd> Waveform pproximant to use for injection.</dd>

<dt><tt>Approximant order</tt></dt><dd> Waveform order to use for injection.</dd>

<dt><tt>REAL4 minMass</tt></dt><dd> Minimum mass of injected signals.</dd>

<dt><tt>REAL4 maxMass</tt></dt><dd> Maximum mass of injected signals.</dd>
</dl>

*/
typedef struct
tagFindChirpBankSimParams
{
  Approximant           approx;
  LALPNOrder            order;
  REAL4                 minMass;
  REAL4                 maxMass;
  RandomParams         *randParams;
  INT4                  maxMatch;
  CHAR                 *frameName;
  CHAR                 *frameChan;
  REAL4			f_lower;
}
FindChirpBankSimParams;

/*
 *
 * function prototypes for memory management functions
 *
 */


void
LALFindChirpCreateTmpltNode (
    LALStatus                  *status,
    InspiralTemplate           *thistmplt,
    InspiralTemplateNode      **tmpltNode
    );

void
LALFindChirpDestroyTmpltNode (
    LALStatus                  *status,
    InspiralTemplateNode      **tmpltNode
    );

void
LALInitializeDataSegmentVector (
    LALStatus                  *status,
    DataSegmentVector         **dataSegVec,
    REAL4TimeSeries            *chan,
    REAL4FrequencySeries       *spec,
    COMPLEX8FrequencySeries    *resp,
    FindChirpInitParams        *params
    );

void
LALFinalizeDataSegmentVector (
    LALStatus                  *status,
    DataSegmentVector         **vector
    );

void
LALCreateDataSegmentVector (
    LALStatus                  *status,
    DataSegmentVector         **vector,
    FindChirpInitParams        *params
    );

void
LALDestroyDataSegmentVector (
    LALStatus                  *status,
    DataSegmentVector         **vector
    );

void
LALCreateFindChirpSegmentVector (
    LALStatus                  *status,
    FindChirpSegmentVector    **vector,
    FindChirpInitParams        *params
    );

void
LALDestroyFindChirpSegmentVector (
    LALStatus                  *status,
    FindChirpSegmentVector    **vector
    );


/*
 *
 * function prototypes for initialization, finalization of data functions
 *
 */

void
LALFindChirpDataInit (
    LALStatus                  *status,
    FindChirpDataParams       **output,
    FindChirpInitParams        *params
    );

void
LALFindChirpDataFinalize (
    LALStatus                  *status,
    FindChirpDataParams       **output
    );


/*
 *
 * function prototypes for initialization, finalization of template functions
 *
 */

void
LALFindChirpTemplateInit (
    LALStatus                  *status,
    FindChirpTmpltParams      **output,
    FindChirpInitParams        *params
    );

void
LALFindChirpTemplateFinalize (
    LALStatus                  *status,
    FindChirpTmpltParams      **output
    );



/*
 *
 * function prototypes for initialization, finalization and filter functions
 *
 */

void
LALFindChirpFilterInit (
    LALStatus                  *status,
    FindChirpFilterParams     **output,
    FindChirpInitParams        *params
    );

void
LALFindChirpFilterFinalize (
    LALStatus                  *status,
    FindChirpFilterParams     **output
    );

void
LALCreateFindChirpInput (
    LALStatus                  *status,
    FindChirpFilterInput      **output,
    FindChirpInitParams        *params
    );

void
LALDestroyFindChirpInput (
    LALStatus                  *status,
    FindChirpFilterInput      **output
    );

void
LALFindChirpCreateCoherentInput(
     LALStatus                  *status,
     COMPLEX8TimeSeries         **coherentInputData,
     COMPLEX8TimeSeries         *input,
     SnglInspiralTable          *templt,
     REAL4                      coherentSegmentLength,
     INT4                       corruptedDataLength
     );

void
LALFindChirpFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    );

void
LALFindChirpStoreEvent (
    LALStatus                  *status,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params,
    SnglInspiralTable          *thisEvent,
    COMPLEX8                   *q,
    UINT4                       kmax,
    REAL4                       norm,
    UINT4                       eventStartIdx,
    UINT4                       numChisqBins,
    CHAR                       *searchName
    );

void
LALFindChirpClusterEvents (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params,
    FindChirpBankVetoData      *bankVetoData,
    UINT4                       subBankIndex,
    int                         writeCData
    );


void
LALFindChirpFilterOutputVeto(
    LALStatus                          *status,
    SnglInspiralTable                 **eventList,
    FindChirpFilterInput               *input,
    FindChirpFilterParams              *fcParams
    );


void
LALFindChirpInjectSignals (
    LALStatus                  *status,
    REAL4TimeSeries            *chan,
    SimInspiralTable           *events,
    COMPLEX8FrequencySeries    *resp
    );

void
LALFindChirpInjectIMR (
    LALStatus                     *status,
    REAL4TimeSeries               *chan,
    SimInspiralTable              *events,
    SimRingdownTable              *ringdownevents,
    COMPLEX8FrequencySeries       *resp,
    INT4                           injectSignalType
    );

INT4
XLALFindChirpTagTemplateAndSegment (
    DataSegmentVector       *dataSegVec,
    InspiralTemplate        *tmpltHead,
    SnglInspiralTable       **events,
    CHAR                    *ifo,
    REAL4                   tdFast,
    UINT4                   *analyseThisTmplt
    );


INT4
XLALFindChirpSetAnalyzeSegment (
    DataSegmentVector          *dataSegVec,
    SimInspiralTable           *injections
    );

INT4
XLALFindChirpSetFollowUpSegment (
    DataSegmentVector          *dataSegVec,
    SnglInspiralTable          **events
    );

void
LALFindChirpSetAnalyseTemplate (
    LALStatus                    *status,
    UINT4                        *analyseThisTmplt,
    REAL4                        mmFast,
    REAL8                        deltaF,
    INT4                         sampleRate,
    FindChirpDataParams          *fcDataParams,
    int                          numTmplts,
    InspiralTemplate             *tmpltHead,
    int                          numInjections,
    SimInspiralTable             *injections
    );

UINT4
XLALCmprSgmntTmpltFlags (
    UINT4 numInjections,
    UINT4 TmpltFlag,
    UINT4 SgmntFlag
    );

UINT4
XLALFindChirpBankSimInitialize (
    REAL4FrequencySeries       *spec,
    COMPLEX8FrequencySeries    *resp,
    REAL8                       fLow
    );

SimInspiralTable *
XLALFindChirpBankSimInjectSignal (
    DataSegmentVector          *dataSegVec,
    COMPLEX8FrequencySeries    *resp,
    SimInspiralTable           *injParams,
    FindChirpBankSimParams     *simParams
    );

REAL4
XLALFindChirpBankSimSignalNorm(
    FindChirpDataParams         *fcDataParams,
    FindChirpSegmentVector      *fcSegVec,
    UINT4                        cut
    );

SimInstParamsTable *
XLALFindChirpBankSimMaxMatch (
    SnglInspiralTable         **bestTmplt,
    REAL4                       matchNorm
    );

SimInstParamsTable *
XLALFindChirpBankSimComputeMatch (
    SnglInspiralTable   *inputTmplt,
    REAL4                matchNorm
    );

FindChirpSubBank*
XLALFindChirpCreateSubBanks(
    UINT4                      *maxSubBankSize,
    UINT4                       subBankSize,
    UINT4                       bankSize,
    InspiralTemplate           *bankHead
    );

void
XLALBankVetoCCMat (
    FindChirpBankVetoData 	*bankVetoData,
    REAL4Vector                 *ampVec,
    UINT4                       subBankSize,
    REAL4 			dynRange,
    REAL4 			fLow,
    REAL4 			deltaF,
    REAL4                       deltaT
    );

REAL4
XLALComputeBankVeto( FindChirpBankVetoData *bankVetoData,
                     UINT4 i,
                     UINT4 snrIX,
		     REAL4 deltaT,
                     UINT4 *dof);


void
XLALInitBankVetoData(
    FindChirpBankVetoData *bvdata
    );

void
XLALDestroyBankVetoData (
    FindChirpBankVetoData *bvdata
    );

REAL4
XLALComputeFullChisq(
    FindChirpBankVetoData      *bankVetoData,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params,
    COMPLEX8                   *q,
    UINT4                       i,
    UINT4                       snrIX,
    UINT4                      *dof,
    REAL4                       norm
);

InspiralTemplate *
XLALFindChirpSortTemplates(
  InspiralTemplate *bankHead,
  FindChirpBankVetoData *bvdata,
  UINT4 num,
  UINT4 max_subbank_size
);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _FINDCHIRPH_H */
