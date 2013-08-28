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

#ifndef _FINDCHIRPDATATYPESH_H
#define _FINDCHIRPDATATYPESH_H

#include <lal/LALDatatypes.h>
#include <lal/LALInspiral.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * \addtogroup FindChirpDatatypes_h
 * \author Brown, D. A.
 *
 * \brief Provides core protypes for the core datatypes using in findchirp.
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/FindChirpDatatypes.h>
 * \endcode
 *
 */
/*@{*/

/* ---------- typedefs of structures used by the findchirp functions ---------- */

/**
 * Struture used to contain a binary inspiral standard candle.
 * \c distance is the distance in Mpc at which an optimally oriented
 * binary with the mass parameters stored in \c tmplt would produce
 * the signal-to-noise ratio squared \c rhosq.
 */
typedef struct
tagFindChirpStandardCandle
{
  CHAR                          ifo[3];		/**< NULL terminated ifo name */
  InspiralTemplate              tmplt;		/**< Binary parameters used to compute the standard candle */
  REAL4                         rhosq;		/**< The signal-to-noise ratio squared \f$\rho^2\f$ of the candle */
  REAL4                         sigmasq;	/**< The variance of the matched filter \f$\sigma^2\f$ for the data used to calculate the standard candle */
  REAL4                         distance;	/**< The distance at which an optimally oriented inspiral with the masses given by \c tmplt would give the signal-to-noise ratio squared \c rhosq */
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

/**
 * Structure used to contain an array of ::DataSegment's.
 * Each \c DataSegment contains an <tt>INT4 number</tt>
 * used to identify the data segment and pointers to a data channel
 * (<tt>REAL4TimeSeries *chan</tt>), a power spectral estimate
 * (<tt>REAL4FrequencySeries *spec</tt>) and a response function
 * (<tt>COMPLEX8FrequencySeries *resp</tt>).
 */
typedef struct
tagDataSegmentVector
{
  UINT4                         length;		/**< Number of \c DataSegment structures in the vector */
  DataSegment                  *data;		/**< Pointer to an array of \c DataSegment structures */
}
DataSegmentVector;


/**
 * This structure provides a method of constucting doubly linked
 * lists of \c InspiralTemplate structures.
 */
typedef struct
tagInspiralTemplateNode
{
  struct tagInspiralTemplateNode       *next;		/**< The next structure in the linked list */
  struct tagInspiralTemplateNode       *prev;		/**< The previous structure in the linked list */
  InspiralTemplate                     *tmpltPtr;	/**< A pointer to an \c InspiralTemplate structure containing the template parameters */
}
InspiralTemplateNode;


/**
 * This structure provides a method of constucting linked
 * lists of \c InspiralTemplateNode structures (as if it were not already
 * complicated enough).  Actually it is necessary to create a list of
 * sub banks for the template bank veto so that roughly 20 or so templates
 * can be filtered and stored in memory at one time.
 */
typedef struct
tagInspiralTemplateNodeList
{
  struct tagInspiralTemplateNodeList       *next;	/**< The next structure in the linked list */
  InspiralTemplateNode                     *nodePtr;	/**< A pointer to an \c InspiralTemplateNode structure */
}
InspiralTemplateNodeList;

/**
 * This structure provides contains subbanks for the template bank veto.
 */
typedef struct
tagFindChirpSubBank
{
  UINT4                         subBankSize;	/**< UNDOCUMENTED */
  InspiralTemplate             *bankHead;	/**< A pointer to an \c InspiralTemplate structure which is the head of linked list of templates */
  struct tagFindChirpSubBank   *next;		/**< The next structure in the linked list */
}
FindChirpSubBank;


/**
 * This structure contains the conditioned input data and its
 * parameters and is one of the inputs to the <tt>FindChirpFilter()</tt>
 * function.
 */
typedef struct
tagFindChirpSegment
{
  COMPLEX8FrequencySeries      *data;			/**< The conditioned data used as
                                                         * part of the matched filter correllation; The exact content of this structure
                                                         * is determined by which data conditioning routine is called (stationary phase,
                                                         * time domain, BCV or spinning BCV); The data in this structure is denoted
                                                         * \f$\tilde{F}_k\f$ and the vetor is of length \f$N/2 + 1\f$; For frequency domain
                                                         * templates (FindChirpSP, BCV and BCVSpin) it contains:
                                                         * \f{equation}{
                                                         * \tilde{F}_k = \frac{d\tilde{v}_k  \left(\frac{k}{N}\right)^{-\frac{7}{6}}}
                                                         * {d^2|R|^2S_v\left(\left|f_k\right|\right)};
                                                         * \f}
                                                         * For time domain templates (GeneratePPN, TaylorT1, TaylorT2, TaylorT3, PadeT1,
                                                         * EOB) it contains
                                                         * \f{equation}{
                                                         * \tilde{F}_k = \frac{d\tilde{v}_k}
                                                         * {d^2|R|^2S_v\left(\left|f_k\right|\right)}
                                                         * \f}
                                                         */
  COMPLEX8FrequencySeries      *dataBCV;		/**< Conditioned input data used
                                                         * only for the BCV templates;  The conditioning performed is as described in the
                                                         * documentation for \ref LALFindChirpBCVData().
                                                         */
  REAL4TimeSeries              *dataPower;		/**< UNDOCUMENTED */
  UINT4Vector                  *chisqBinVec;		/**< A vector containing the indices of
                                                         * the boundaries of the bins of equal power for the \f$\chi^2\f$ veto for this
                                                         * segment; The vector is of length \f$p+1\f$, where \f$p\f$ is the number of \f$\chi^2\f$
                                                         * bins; If no \f$\chi^2\f$ veto is performed, this may be NULL
                                                         */
  UINT4Vector                  *chisqBinVecBCV;		/**< A vector containing the indices of
                                                         * the boundaries of the bins of equal power for the second contribution to the
                                                         * \f$\chi^2\f$ statistic for the BCV templates for this segment
                                                         */
  REAL8                         deltaT;			/**< The time step \f$\Delta t\f$ of the input data channel used to create the \c FindChirpSegment */
  REAL4Vector                  *segNorm;		/**< The quantity segment dependent
                                                         * normalization quantity \f$\mathcal{S}_k\f$; The vector is of length \f$N/2+1\f$;
                                                         * For stationary phase templates the segment dependent normalization is
                                                         * \anchor eq_spsegnorm \f{equation}{
                                                         * \mathcal{S}_k = \sum_{k^\prime=1}^{k}
                                                         * \frac{\left(\frac{k^\prime}{N}\right)^{-\frac{7}{3}}}{d^2|R|^2S_v\left(\left|f_{k^\prime}\right|\right)} \quad\quad 1 \le k \le N/2
                                                         * \tag{eq_spsegnorm}
                                                         * \f}
                                                         * and can be computed once per data segment and re-used for each template;  For
                                                         * time domain templates, the segment dependent normalization is
                                                         * \anchor eq_tdsegnorm \f{equation}{
                                                         * \mathcal{S}_k = \sum_{k^\prime=1}^{k}
                                                         * \frac{\tilde{h}_{k^\prime}}{d^2|R|^2S_v\left(\left|f_{k^\prime}\right|\right)} \quad\quad 1 \le k \le N/2
                                                         * \tag{eq_tdsegnorm}
                                                         * \f}
                                                         * and it must be recomputed for each template \f$\tilde{h}_k\f$
                                                         */
  REAL4Vector                  *tmpltPowerVec;		/**< %Vector of length \f$N/2+1\f$ containing
                                                         * the weighted power in the template; For frequency domain templates, this is
                                                         * the summand in equation\eqref{eq_spsegnorm}
                                                         * \f{equation}{
                                                         * \mathtt{templtPowerVec->data[k]} =
                                                         * \frac{\left(\frac{k}{N}\right)^{-\frac{7}{3}}}{d^2|R|^2S_v\left(\left|f_k\right|\right)}
                                                         * \f}
                                                         * and can be computed once then re-used for all templates;  For time domain
                                                         * templates, this is the summand in equation\eqref{eq_spsegnorm}
                                                         * \f{equation}{
                                                         * \mathtt{templtPowerVec->data[k]} =
                                                         * \frac{\tilde{h}_{k}}{d^2|R|^2S_v\left(\left|f_k\right|\right)}
                                                         * \f}
                                                         * which must be re-computed for each template \f$\tilde{h}_k\f$;  This quantity is
                                                         * used in the computation of the \f$\chi^2\f$ bin boundaries and the re-computation of
                                                         * \f$\mathcal{S}_k\f$ for time domain templates
                                                         */
  REAL4Vector                  *a1;			/**< BCV-template normalization parameter */
  REAL4Vector                  *b1;			/**< BCV-template normalization parameter */
  REAL4Vector                  *b2;			/**< BCV-template normalization parameter */
  REAL4Vector                  *tmpltPowerVecBCV;	/**< Additional weighted template power for BCV templates */
  REAL4                         fLow;			/**< The (frequency domain) low frequency cutoff for the matched filter, \f$f_\mathrm{low}\f$ */
  UINT4                         invSpecTrunc;		/**< The number of points to which the inverse
                                                         * power spectrum \f${S\left(\left|f_{k}\right|\right)}\f$ is truncated to in the time domain in order to truncate
                                                         * the impulse response time of the matched filter
                                                         */
  UINT4                         number;			/**< A unique identification number for the
                                                         * \c FindChirpDataSegment; This corresponds to the number in
                                                         * the \c DataSegment from which the conditioned data was computed
                                                         */
  UINT4                         analyzeSegment;		/**< UNDOCUMENTED */
  INT4                          level;			/**< A search level, used by the heirarchical search
                                                         * engine to determine if a particular data segment should be filtered against
                                                         * a particular template
                                                         */
  Approximant                   approximant;		/**< UNDOCUMENTED */
}
FindChirpSegment;

/**
 * A vector of \c FindChirpSegment structures, defined above.
 */
typedef struct
tagFindChirpSegmentVector
{
  UINT4                         length;	/**< Number of \c FindChirpSegment structures in the vector */
  FindChirpSegment             *data;	/**< Pointer to an array of \c FindChirpSegment structures */
}
FindChirpSegmentVector;

/**
 * This structure contains a frequency domain template used as input
 * to the <tt>FindChirpFilter()</tt> routine. This may either be a template
 * generated in the frequency domain or the Fourier transform of template
 * generated in the time domain.
 */
typedef struct
tagFindChirpTemplate
{
  InspiralTemplate              tmplt;		/**< The template parameters of this
                                                 * \c FindChirpTemplate; In addition to the mass parameters the following
                                                 * fields of \c tmplt should populated by the template generation functions
                                                 * as the are used by <tt>FindChirpFilterSegment()</tt>:
                                                 * <dl>
                                                 * <dt>approximant</dt><dd> Used to check that the findchirp data segment
                                                 * and the template have been created for the same type of waveform</dd>
                                                 * <dt>tC</dt><dd> The length of the chirp in seconds; Used by the max over
                                                 * chirp event finding algorithm</dd>
                                                 * <dt>fFinal</dt><dd> The highest frequency component of the chirp; Used to
                                                 * pick the appropriate value of the segment normalization constant
                                                 * \f$\mathcal{S}_k\f$ for this template
                                                 * </dl>
                                                 */
  COMPLEX8Vector               *data;		/**< %Vector of length \f$N/2+1\f$ containing the
                                                 * frequency template data \f$\tilde{T}_k\f$; For a template generated in the frequency
                                                 * domain template (FindChirpSP) this should contain
                                                 * \f{equation}{
                                                 * \tilde{T}_k = \exp\left[i\Psi(f_k;M,\eta)\right] \Theta\left(k-k_\mathrm{isco}\right)
                                                 * \f}
                                                 * For a template generated in the time domain this should contain the discrete
                                                 * Fourier transform of the cosine phase chirp
                                                 * \f{equation}{
                                                 * \tilde{T}_k = \tilde{h}_{ck} = \mathrm{DFT}\left[ h(t) \right]
                                                 * \f}
                                                 * where \f$h(t)\f$ is an inspiral waveform generated by the
                                                 * <tt>LALInspiralWave()</tt> function if the approximant TaylorT1, TaylorT2,
                                                 * TaylorT3, PadeT1 or EOB; Alternatively \f$h(t)\f$ can be generated by the
                                                 * <tt>LALGeneratePPNInspiral()</tt> function if the approximant is GeneratePPN;
                                                 * Findchirp always uses second order post-Newtonian templates
                                                 */
  COMPLEX8VectorSequence       *ACTDtilde;	/**< UNDOCUMENTED */
  REAL4VectorSequence          *PTFQ;		/**< UNDOCUMENTED */
  COMPLEX8VectorSequence       *PTFQtilde;	/**< UNDOCUMENTED */
  REAL4Array                   *PTFBinverse;	/**< UNDOCUMENTED */
  REAL4Array                   *PTFB;		/**< UNDOCUMENTED */
  REAL4                         tmpltNorm;	/**< The template dependent normalisation constant
                                                 * \f$\mathcal{T}\f$; For the stationary phase template FindChirpSP this is
                                                 * \f{equation}{
                                                 * \mathcal{T}(M,\mu) = \left[
                                                 * \left(\frac{2dGM_\odot}{(1\,\mathrm{Mpc})c^2}\right)
                                                 * \left(\frac{5\mu}{96M_\odot}\right)^\frac{1}{2}
                                                 * \left(\frac{M}{\pi^2M_\odot}\right)^\frac{1}{3}
                                                 * \left(\frac{GM_\odot}{\Delta tc^3}\right)^{-\frac{1}{6}}
                                                 * \right]^2
                                                 * \f}
                                                 * where \f$d\f$ is the dynamic range parameter \c dynRange;
                                                 * For time domain templates generated by <tt>LALInspiralWave()</tt> (TaylorT1,
                                                 * TaylorT2, TaylorT3, PadeT1 and EOB) this is
                                                 * \f{equation}{
                                                 * \mathcal{T}(\mu) = \left[ \left(\frac{4dGM_\odot}{(1\,\mathrm{Mpc})c^2}\right)
                                                 * \left(\frac{\mu}{M_\odot}\right) \right]^2;
                                                 * \f}
                                                 * For time domain templates generated by <tt>LALGeneratePPNInspiral()</tt>
                                                 * (GeneratePPN) it is
                                                 * \f{equation}{
                                                 * \mathcal{T} = \left(\frac{d}{1\,\mathrm{Mpc}}\right)^2
                                                 * \f}
                                                 */
  REAL4				norm;		/**< UNDOCUMENTED */
  REAL8                         momentI;	/**< Undocumented BCV normalization constant */
  REAL8                         momentJ;	/**< Undocumented BCV normalization constant */
  REAL8                         momentK;	/**< Undocumented BCV normalization constant */
  REAL8                         rootMomentI;	/**< Undocumented BCV normalization constant */
  REAL8                         numFactor;	/**< Undocumented BCV normalization constant */
  REAL8                         numFactor1;	/**< Undocumented BCV normalization constant */
  REAL8                         numFactor2;	/**< Undocumented BCV normalization constant */
  REAL8                         numFactor3;	/**< Undocumented BCV normalization constant */
  REAL8Vector                  *A1BCVSpin;	/**< Undocumented spinning BCV template data */
  REAL8Vector                  *A2BCVSpin;	/**< Undocumented spinning BCV template data */
  REAL8Vector                  *A3BCVSpin;	/**< Undocumented spinning BCV template data */
}
FindChirpTemplate;

/*@}*/ /* end:FindChirpDatatypes_h */

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _FINDCHIRPDATATYPESH_H */
