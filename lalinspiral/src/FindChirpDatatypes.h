/*
*  Copyright (C) 2007 Chad Hanna, Duncan Brown, Gareth Jones, Jolien Creighton, Patrick Brady, Robert Adam Mercer
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
 * File Name: FindChirpDatatypes.h
 *
 * Author: Brown, D. A.
 *
 *-----------------------------------------------------------------------
 */

/**

\author Brown, D. A.
\file
\ingroup CBC_findchirp

\brief Provides core protypes for the core datatypes using in findchirp.

\heading{Synopsis}
\code
#include <lal/FindChirpDatatypes.h>
\endcode

*/

#ifndef _FINDCHIRPDATATYPESH_H
#define _FINDCHIRPDATATYPESH_H

#include <lal/LALDatatypes.h>
#include <lal/LALInspiral.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif


NRCSID (FINDCHIRPDATATYPESH, "$Id$");


/*
 *
 * typedefs of structures used by the findchirp functions
 *
 */


/** Struture used to contain a binary inspiral standard candle.
\c distance is the distance in Mpc at which an optimally oriented
binary with the mass parameters stored in \c tmplt would produce
the signal-to-noise ratio squared \c rhosq.

<dl>
<dt><tt>CHAR ifo[3]</tt></dt><dd> NULL terminated ifo name.</dd>

<dt><tt>InspiralTemplate tmplt</tt></dt><dd> Binary parameters used to compute the
standard candle.</dd>

<dt><tt>REAL4 rhosq</tt></dt><dd> The signal-to-noise ratio squared \f$\rho^2\f$ of the
candle.</dd>

<dt><tt>REAL4 sigmasq</tt></dt><dd> The variance of the matched filter \f$\sigma^2\f$
for the data used to calculate the standard candle.</dd>

<dt><tt>REAL4 distance</tt></dt><dd> The distance at which an optimally oriented
inspiral with the masses given by \c tmplt would give the
signal-to-noise ratio squared \c rhosq.</dd>
</dl>

*/
typedef struct
tagFindChirpStandardCandle
{
  CHAR                          ifo[3];
  InspiralTemplate              tmplt;
  REAL4                         rhosq;
  REAL4                         sigmasq;
  REAL4                         distance;
}
FindChirpStandardCandle;


typedef struct
tagDataSegment
{
  REAL4TimeSeries         *chan;
  REAL4FrequencySeries    *spec;
  COMPLEX8FrequencySeries *resp;
  INT4                     number;
  UINT4                    analyzeSegment;
}
DataSegment;

/** Structure used to contain an array of ::DataSegment's.
Each \c DataSegment contains an <tt>INT4 number</tt>
used to identify the data segment and pointers to a data channel
(<tt>REAL4TimeSeries *chan</tt>), a power spectral estimate
(<tt>REAL4FrequencySeries *spec</tt>) and a response function
(<tt>COMPLEX8FrequencySeries *resp</tt>).

<dl>
<dt><tt>UINT4 length</tt></dt><dd> Number of \c DataSegment structures in the
vector.</dd>

<dt><tt>DataSegment *data</tt></dt><dd> Pointer to an array of \c DataSegment
structures.</dd>
</dl>

*/
typedef struct
tagDataSegmentVector
{
  UINT4                         length;
  DataSegment                  *data;
}
DataSegmentVector;


/** This structure provides a method of constucting doubly linked
lists of \c InspiralTemplate structures.

<dl>
<dt><tt>struct tagInspiralTemplateNode *next</tt></dt><dd> The next structure in
the linked list.</dd>

<dt><tt>struct tagInspiralTemplateNode *prev</tt></dt><dd> The previous structure in
the linked list.</dd>

<dt><tt>InspiralTemplate *tmpltPtr</tt></dt><dd> A pointer to an
\c InspiralTemplate structure containing the template parameters.</dd>
</dl>

*/
typedef struct
tagInspiralTemplateNode
{
  struct tagInspiralTemplateNode       *next;
  struct tagInspiralTemplateNode       *prev;
  InspiralTemplate                     *tmpltPtr;
}
InspiralTemplateNode;


/** This structure provides a method of constucting linked
lists of \c InspiralTemplateNode structures (as if it were not already
complicated enough).  Actually it is necessary to create a list of
sub banks for the template bank veto so that roughly 20 or so templates
can be filtered and stored in memory at one time.

<dl>
<dt><tt>struct tagInspiralTemplateNodeList *next</tt></dt><dd> The next structure in
the linked list.</dd>

<dt><tt>InspiralTemplateNode *nodePtr</tt></dt><dd> A pointer to an
\c InspiralTemplateNode structure.</dd>
</dl>

*/
typedef struct
tagInspiralTemplateNodeList
{
  struct tagInspiralTemplateNodeList       *next;
  InspiralTemplateNode                     *nodePtr;
}
InspiralTemplateNodeList;

/** This structure provides contains subbanks for the template bank veto.
<dl>
<dt><tt>InspiralTemplate *bankHead</tt></dt><dd> A pointer to an
\c InspiralTemplate structure which is the head of linked list of
templates.</dd>

<dt><tt>struct tagFindChirpSubBank *next</tt></dt><dd> The next structure in
the linked list.</dd>
</dl>
 */
typedef struct
tagFindChirpSubBank
{
  UINT4                         subBankSize;
  InspiralTemplate             *bankHead;
  struct tagFindChirpSubBank   *next;
}
FindChirpSubBank;


/** This structure contains the conditioned input data and its
parameters and is one of the inputs to the <tt>FindChirpFilter()</tt>
function.

<dl>
<dt><tt>COMPLEX8FrequencySeries *data</tt></dt><dd> The conditioned data used as
part of the matched filter correllation. The exact content of this structure
is determined by which data conditioning routine is called (stationary phase,
time domain, BCV or spinning BCV). The data in this structure is denoted
\f$\tilde{F}_k\f$ and the vetor is of length \f$N/2 + 1\f$. For frequency domain
templates (FindChirpSP, BCV and BCVSpin) it contains:
\f{equation}{
\tilde{F}_k = \frac{d\tilde{v}_k  \left(\frac{k}{N}\right)^{-\frac{7}{6}}}
{d^2|R|^2S_v\left(\left|f_k\right|\right)}.
\f}
For time domain templates (GeneratePPN, TaylorT1, TaylorT2, TaylorT3, PadeT1,
EOB) it contains
\f{equation}{
\tilde{F}_k = \frac{d\tilde{v}_k}
{d^2|R|^2S_v\left(\left|f_k\right|\right)}.
\f}</dd>

<dt><tt>COMPLEX8FrequencySeries *dataBCV</tt></dt><dd> Conditioned input data used
only for the BCV templates.  The conditioning performed is as described in the
documentation for the module \ref FindChirpBCVData.c</dd>

<dt><tt>UINT4Vector *chisqBinVec</tt></dt><dd> A vector containing the indices of
the boundaries of the bins of equal power for the \f$\chi^2\f$ veto for this
segment. The vector is of length \f$p+1\f$, where \f$p\f$ is the number of \f$\chi^2\f$
bins. If no \f$\chi^2\f$ veto is performed, this may be NULL.</dd>

<dt><tt>UINT4Vector *chisqBinVecBCV</tt></dt><dd> A vector containing the indices of
the boundaries of the bins of equal power for the second contribution to the
\f$\chi^2\f$ statistic for the BCV templates for this segment.</dd>

<dt><tt>REAL8 deltaT</tt></dt><dd> The time step \f$\Delta t\f$ of the input
data channel used to create the \c FindChirpSegment.</dd>

<dt><tt>REAL4Vector *segNorm</tt></dt><dd> The quantity segment dependent
normalization quantity \f$\mathcal{S}_k\f$. The vector is of length \f$N/2+1\f$.
For stationary phase templates the segment dependent normalization is
\anchor eq_spsegnorm \f{equation}{
\mathcal{S}_k = \sum_{k^\prime=1}^{k}
\frac{\left(\frac{k^\prime}{N}\right)^{-\frac{7}{3}}}{d^2|R|^2S_v\left(\left|f_{k^\prime}\right|\right)} \quad\quad 1 \le k \le N/2
\label{eq_spsegnorm}
\f}
and can be computed once per data segment and re-used for each template.  For
time domain templates, the segment dependent normalization is
\anchor eq_tdsegnorm \f{equation}{
\mathcal{S}_k = \sum_{k^\prime=1}^{k}
\frac{\tilde{h}_{k^\prime}}{d^2|R|^2S_v\left(\left|f_{k^\prime}\right|\right)} \quad\quad 1 \le k \le N/2
\label{eq_tdsegnorm}
\f}
and it must be recomputed for each template \f$\tilde{h}_k\f$.</dd>

<dt><tt>REAL4Vector *tmpltPowerVec</tt></dt><dd> Vector of length \f$N/2+1\f$ containing
the weighted power in the template. For frequency domain templates, this is
the summand in equation\eqref{eq_spsegnorm}
\f{equation}{
\mathtt{templtPowerVec->data[k]} =
\frac{\left(\frac{k}{N}\right)^{-\frac{7}{3}}}{d^2|R|^2S_v\left(\left|f_k\right|\right)}.
\f}
and can be computed once then re-used for all templates.  For time domain
templates, this is the summand in equation\eqref{eq_spsegnorm}
\f{equation}{
\mathtt{templtPowerVec->data[k]} =
\frac{\tilde{h}_{k}}{d^2|R|^2S_v\left(\left|f_k\right|\right)}.
\f}
which must be re-computed for each template \f$\tilde{h}_k\f$.  This quantity is
used in the computation of the \f$\chi^2\f$ bin boundaries and the re-computation of
\f$\mathcal{S}_k\f$ for time domain templates.</dd>

<dt><tt>REAL4 a1</tt></dt><dd> BCV-template normalization parameter.</dd>

<dt><tt>REAL4 b1</tt></dt><dd> BCV-template normalization parameter.</dd>

<dt><tt>REAL4 b2</tt></dt><dd> BCV-template normalization parameter.</dd>

<dt><tt>REAL4Vector *tmpltPowerVecBCV</tt></dt><dd> Additional weighted template
power for BCV templates.</dd>

<dt><tt>REAL4 fLow</tt></dt><dd> The (frequency domain) low frequency cutoff for the
matched filter, \f$f_\mathrm{low}\f$.</dd>

<dt><tt>UINT4 invSpecTrunc</tt></dt><dd> The number of points to which the inverse
power spectrum \f${S\left(\left|f_{k}\right|\right)}\f$ is truncated to in the time domain in order to truncate
the impulse response time of the matched filter.</dd>

<dt><tt>UINT4 number</tt></dt><dd> A unique identification number for the
\c FindChirpDataSegment. This corresponds to the number in
the \c DataSegment from which the conditioned data was computed.</dd>

<dt><tt>INT4 level</tt></dt><dd> A search level, used by the heirarchical search
engine to determine if a particular data segment should be filtered against
a particular template.</dd>
</dl>

*/
typedef struct
tagFindChirpSegment
{
  COMPLEX8FrequencySeries      *data;
  COMPLEX8FrequencySeries      *dataBCV;
  REAL4TimeSeries              *dataPower;
  UINT4Vector                  *chisqBinVec;
  UINT4Vector                  *chisqBinVecBCV;
  REAL8                         deltaT;
  REAL4Vector                  *segNorm;
  REAL4Vector                  *tmpltPowerVec;
  REAL4Vector                  *a1;
  REAL4Vector                  *b1;
  REAL4Vector                  *b2;
  REAL4Vector                  *tmpltPowerVecBCV;
  REAL4                         fLow;
  UINT4                         invSpecTrunc;
  UINT4                         number;
  UINT4                         analyzeSegment;
  INT4                          level;
  Approximant                   approximant;
}
FindChirpSegment;

/** A vector of \c FindChirpSegment structures, defined above.

<dl>
<dt><tt>UINT4 length</tt></dt><dd> Number of \c FindChirpSegment structures in
the vector</dd>

<dt><tt>DataSegment *data</tt></dt><dd> Pointer to an array of
\c FindChirpSegment structures.</dd>
</dl>


*/
typedef struct
tagFindChirpSegmentVector
{
  UINT4                         length;
  FindChirpSegment             *data;
}
FindChirpSegmentVector;

/** This structure contains a frequency domain template used as input
to the <tt>FindChirpFilter()</tt> routine. This may either be a template
generated in the frequency domain or the Fourier transform of template
generated in the time domain.

<dl>
<dt><tt>InspiralTemplate tmplt</tt></dt><dd> The template parameters of this
\c FindChirpTemplate. In addition to the mass parameters the following
fields of \c tmplt should populated by the template generation functions
as the are used by <tt>FindChirpFilterSegment()</tt>:</dd>
<dl>
<dt>approximant</dt><dd> Used to check that the findchirp data segment
and the template have been created for the same type of waveform.</dd>
<dt>tC</dt><dd> The length of the chirp in seconds. Used by the max over
chirp event finding algorithm.</dd>
<dt>fFinal</dt><dd> The highest frequency component of the chirp. Used to
pick the appropriate value of the segment normalization constant
\f$\mathcal{S}_k\f$ for this template.</dd>
</dl>

<dt><tt>COMPLEX8Vector *data</tt></dt><dd> Vector of length \f$N/2+1\f$ containing the
frequency template data \f$\tilde{T}_k\f$. For a template generated in the frequency
domain template (FindChirpSP) this should contain
\f{equation}{
\tilde{T}_k = \exp\left[i\Psi(f_k;M,\eta)\right] \Theta\left(k-k_\mathrm{isco}\right).
\f}
For a template generated in the time domain this should contain the discrete
Fourier transform of the cosine phase chirp
\f{equation}{
\tilde{T}_k = \tilde{h}_{ck} = \mathrm{DFT}\left[ h(t) \right]
\f}
where \f$h(t)\f$ is an inspiral waveform generated by the
<tt>LALInspiralWave()</tt> function if the approximant TaylorT1, TaylorT2,
TaylorT3, PadeT1 or EOB. Alternatively \f$h(t)\f$ can be generated by the
<tt>LALGeneratePPNInspiral()</tt> function if the approximant is GeneratePPN.
Findchirp always uses second order post-Newtonian templates.</dd>

<dt><tt>REAL4 tmpltNorm</tt></dt><dd> The template dependent normalisation constant
\f$\mathcal{T}\f$. For the stationary phase template FindChirpSP this is
\f{equation}{
\mathcal{T}(M,\mu) = \left[
\left(\frac{2dGM_\odot}{(1\,\mathrm{Mpc})c^2}\right)
\left(\frac{5\mu}{96M_\odot}\right)^\frac{1}{2}
\left(\frac{M}{\pi^2M_\odot}\right)^\frac{1}{3}
\left(\frac{GM_\odot}{\Delta tc^3}\right)^{-\frac{1}{6}}
\right]^2
\f}
where \f$d\f$ is the dynamic range parameter \c dynRange.
For time domain templates generated by <tt>LALInspiralWave()</tt> (TaylorT1,
TaylorT2, TaylorT3, PadeT1 and EOB) this is
\f{equation}{
\mathcal{T}(\mu) = \left[ \left(\frac{4dGM_\odot}{(1\,\mathrm{Mpc})c^2}\right)
  \left(\frac{\mu}{M_\odot}\right) \right]^2.
\f}
For time domain templates generated by <tt>LALGeneratePPNInspiral()</tt>
(GeneratePPN) it is
\f{equation}{
\mathcal{T} = \left(\frac{d}{1\,\mathrm{Mpc}}\right)^2.
\f}</dd>

<dt><tt>REAL8 momentI</tt></dt><dd> Undocumented BCV normalization constant.</dd>

<dt><tt>REAL8 momentJ</tt></dt><dd> Undocumented BCV normalization constant.</dd>

<dt><tt>REAL8 momentK</tt></dt><dd> Undocumented BCV normalization constant.</dd>

<dt><tt>REAL8 rootMomentI</tt></dt><dd> Undocumented BCV normalization constant.</dd>

<dt><tt>REAL8 numFactor</tt></dt><dd> Undocumented BCV normalization constant.</dd>

<dt><tt>REAL8 numFactor1</tt></dt><dd> Undocumented BCV normalization constant.</dd>

<dt><tt>REAL8 numFactor2</tt></dt><dd> Undocumented BCV normalization constant.</dd>

<dt><tt>REAL8 numFactor3</tt></dt><dd> Undocumented BCV normalization constant.</dd>

<dt><tt>REAL8 A1BCVSpin</tt></dt><dd> Undocumented spinning BCV template data.</dd>

<dt><tt>REAL8 A2BCVSpin</tt></dt><dd> Undocumented spinning BCV template data.</dd>

<dt><tt>REAL8 A3BCVSpin</tt></dt><dd> Undocumented spinning BCV template data.</dd>
</dl>

*/
typedef struct
tagFindChirpTemplate
{
  InspiralTemplate              tmplt;
  COMPLEX8Vector               *data;
  COMPLEX8VectorSequence       *ACTDtilde;
  COMPLEX8VectorSequence       *PTFQtilde;
  REAL4Array                   *PTFBinverse;
  REAL4Array                   *PTFB;
  REAL4                         tmpltNorm;
  REAL4				norm;
  REAL8                         momentI;
  REAL8                         momentJ;
  REAL8                         momentK;
  REAL8                         rootMomentI;
  REAL8                         numFactor;
  REAL8                         numFactor1;
  REAL8                         numFactor2;
  REAL8                         numFactor3;
  REAL8Vector                  *A1BCVSpin;
  REAL8Vector                  *A2BCVSpin;
  REAL8Vector                  *A3BCVSpin;
}
FindChirpTemplate;

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _FINDCHIRPDATATYPESH_H */
