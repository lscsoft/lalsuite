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
 * \defgroup FindChirpDatatypes_h Header FindChirpDatatypes.h
 * \ingroup lalinspiral_findchirp
 * \author Brown, D. A.
 *
 * \brief Provides core protypes for the core datatypes using in findchirp.
 *
 * ### Synopsis ###
 *
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
