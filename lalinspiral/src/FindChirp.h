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

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif


/**
 * \defgroup FindChirp_h Header FindChirp.h
 * \ingroup lalinspiral_findchirp
 * \author Allen, B., Brown, D. A. and Creighton, J. D. E.
 *
 * \brief This header provides core prototypes, structures and functions to
 * filter interferometer data for binary inspiral chirps.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/FindChirp.h>
 * \endcode
 *
 * Each function in findchirp falls into one of four classes:
 * <ol>
 * <li> Generate management functions which are independent of the type of
 * filtering implemented. The prototypes for these functions are provided by
 * this header file.</li>
 *
 * <li> Functions for filtering data for time domain and frequency domain
 * templates with an unknown amplitude and phase. These are the functions
 * that implement matched filtering for time domain templates (TaylorT1,
 * TaylorT2, TaylorT2, PadeT1, EOB and GeneratePPN) and matched filtering
 * for post-Newtonian frequency domain templates (FindChirpSP). The main
 * filter function <tt>FindChirpFilterSegment()</tt> is prototyped in
 * this header file.
 * Full documentation of the filtering algorithm used can be found in the
 * documentation of the module \ref FindChirpFilter.c.</li>
 * </ol>
 *
 * The goal of all the filtering functions is to determine if the
 * (calibrated) output of the interferometer \f$s(t)\f$ contains a gravitational wave
 * \f$h(t)\f$ in the presence of the detector noise \f$n(t)\f$. When the interferometer
 * is operating properly
 * \f{equation}{
 * s(t) = \left\{ \begin{array}{ll}
 * n(t) + h(t) & \textrm{signal present},\\
 * n(t) & \textrm{signal absent}.
 * \end{array}\right.
 * \f}
 * The detection of signals of known form in noise is a classic problem of signal
 * processing \cite WeinZub62 and can be answered by the construction of a
 * <em>detection statistic</em> and a test to see if the statistic is above some
 * pre-assigned threshold. The construction of the various detection
 * statistics used for each the three types of search are described in the modules
 * that implement the search.
 *
 */
/*@{*/

/**\name Error Codes */
/*@{*/
#define FINDCHIRPH_ENULL 1	/**< Null pointer */
#define FINDCHIRPH_ENNUL 2	/**< Non-null pointer */
#define FINDCHIRPH_EALOC 3	/**< Memory allocation error */
#define FINDCHIRPH_ENUMZ 5	/**< Invalid number of points in segment */
#define FINDCHIRPH_ESEGZ 6	/**< Invalid number of segments */
#define FINDCHIRPH_ECHIZ 7	/**< Invalid number of chi squared bins */
#define FINDCHIRPH_EDTZO 8	/**< deltaT is zero or negative */
#define FINDCHIRPH_ETRNC 10	/**< Duration of inverse spectrum in time domain is negative */
#define FINDCHIRPH_EFLOW 11	/**< Inverse spectrum low frequency cutoff is negative */
#define FINDCHIRPH_EFREE 12	/**< Error freeing memory */
#define FINDCHIRPH_ERHOT 15	/**< Rhosq threshold is negative */
#define FINDCHIRPH_ECHIT 16	/**< Chisq threshold is negative */
#define FINDCHIRPH_ECRUP 17	/**< Chirp length or invSpecTrunc too long for length of data segment */
#define FINDCHIRPH_ESMSM 18	/**< Size mismatch between vectors */
#define FINDCHIRPH_EHETR 19	/**< Attempting to simulate heterodyned GW */
#define FINDCHIRPH_EDFDT 20	/**< Waveform sampling interval is too large */
#define FINDCHIRPH_EAPRX 21	/**< Incorrect waveform approximant */
#define FINDCHIRPH_EUAPX 22	/**< Unknown waveform approximant */
#define FINDCHIRPH_ECHTZ 23	/**< Length of chirp is zero or negative */
#define FINDCHIRPH_EMASS 24	/**< Invalid mass parameters for template generation */
#define FINDCHIRPH_EWVFM 25	/**< Unknown injection waveform */
#define FINDCHIRPH_EBCVC 25	/**< BCVC code: thetav not in [-pi, pi]. */
#define FINDCHIRPH_EMAPX 26	/**< Mismatch in waveform approximant */
#define FINDCHIRPH_EPTFW 27	/**< Error generating PTF waveform */
#define FINDCHIRPH_EIGEN 28	/**< Error computing eigenvalues */
#define FINDCHIRPH_EIMRW 29	/**< Error computing IMR waveform */
#define FINDCHIRPH_EFLOX 30     /**< Error computing variable flower */
/*@}*/

/** \cond DONT_DOXYGEN */
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
#define FINDCHIRPH_MSGEFLOX "Error calculating an appropriate value of f_low"
/** \endcond */

/**
 * This structure contains the parameters for generation of templates
 * by the various template generation functions provided in \ref lalinspiral_findchirp.
 */
typedef struct
tagFindChirpTmpltParams
{
  REAL8                         deltaT;			/**< The sampling interval \f$\Delta t\f$ of the input data channel */
  REAL4                         fLow;			/**< The frequency domain low frequency cutoff \f$f_\mathrm{low}\f$; All frequency domain data is zero below this frequency */
  REAL4                         dynRange;		/**< A dynamic range factor \f$d\f$ which cancels from the filter output; This allows quantities to be stored in the range of
                                                         * \c REAL4 rather than \c REAL8; This must be set to the same value as \c dynRange in the \c FindChirpDataParams; For LIGO data a
                                                         * value of \f$d = 2^{69}\f$ is appropriate
                                                         */
  UINT4                         invSpecTrunc;           /**< The length to which to truncate the inverse power spectral density of the data in the time domain; If set to zero, no
                                                         * truncation is performed
                                                         */

  REAL4Vector                  *xfacVec;		/**< For frequency domain templates, this is a
                                                         * vector of length \f$N/2+1\f$ which contains the quantity \f$k^{-1/3}\f$; For time
                                                         * domain templates, this is a workspace vector of length \f$N\f$ which contains the
                                                         * time domain template generated by the inspiral package, shifted so that the
                                                         * end of the template is at the end of the vector; This vector is Fourier
                                                         * transformed to obtain the quantity findchirp template \f$\tilde{T}_k\f$
                                                         */
  REAL4VectorSequence          *ACTDVecs;		/**< UNDOCUMENTED */
  REAL4Vector                  *PTFphi;			/**< UNDOCUMENTED */
  REAL4Vector                  *PTFomega_2_3;		/**< UNDOCUMENTED */
  REAL4VectorSequence          *PTFe1;			/**< UNDOCUMENTED */
  REAL4VectorSequence          *PTFe2;			/**< UNDOCUMENTED */
  RealFFTPlan                  *fwdPlan;		/**< For time domain templates, an FFTW plan
                                                         * used to transform the time domain data stored in \c xfacVec into its DFT
                                                         * which is stored in the findchirp template
                                                         */
  Approximant                   approximant;		/**< Generate templates of type
                                                         * \c approximant; Valid approximants are #TaylorT1, #TaylorT2, #TaylorT3,
                                                         * #PadeT1, #EOB, #FindChirpSP, #BCV and #BCVSpin; For time domain templates the
                                                         * post-Newtonian order is always two; For stationary phase templates, the
                                                         * post-Newtonian order is specified by \c order
                                                         */
  LALPNOrder                    order;			/**< Specifies the post-Newtonian order of the
                                                         * templates; Valid pN orders are #LAL_PNORDER_TWO, #LAL_PNORDER_TWO_POINT_FIVE, #LAL_PNORDER_THREE,
                                                         * #LAL_PNORDER_THREE_POINT_FIVE, #LAL_PNORDER_PSEUDO_FOUR; The latter is not the true four PN
                                                         * correction, but may increase the fitting factor between stationary phase and
                                                         * numerical relativity waveforms
                                                         */
  INT4                          reverseChirpBank;	/**< Switches a FindChirpSP template bank to be a reverse chirp template bank if true */
  INT4                          bandPassTmplt;		/**< UNDOCUMENTED */
  LALSimInspiralApplyTaper      taperTmplt;		/**< UNDOCUMENTED */
  INT4                          dynamicTmpltFlow;       /**< Use longest template that will fit in pad length */
  REAL4                         maxTempLength;                   /**< This can be used to store the maximum allowed template length, given the pad length and spectrum truncation */
}
FindChirpTmpltParams;

/**
 * This structure contains the possible methods by which
 * to maximize over a chirp in a data segment.
 */
typedef enum
tagFindChirpClustering
{
  FindChirpClustering_none,		/**< The decision to do no clustering of events */
  FindChirpClustering_tmplt,		/**< Cluster over the length of the data segment */
  FindChirpClustering_window,		/**< Cluster over a given number of seconds given by the argument to the flag
                                         * <tt>--cluster-window</tt> (required to be less than the length of the data segment) */
  FindChirpClustering_tmpltwindow	/**< UNDOCUMENTED */
}
FindChirpClustering;


/**
 * This structure provides the parameters for the filter output veto.
 */
typedef struct
tagFindChirpFilterOutputVetoParams
{
  REAL4          rsqvetoWindow;		/**< Width of the \f$r^2\f$ veto window in units of seconds */
  REAL4          rsqvetoThresh;		/**< Threshold of the \f$r^2\f$ veto test analogous to the \f$r^2\f$ threshold employed in the bns and macho inspiral searches */
  REAL4          rsqvetoTimeThresh;	/**< UNDOCUMENTED */
  REAL4          rsqvetoMaxSNR;		/**< UNDOCUMENTED */
  REAL4          rsqvetoCoeff;		/**< UNDOCUMENTED */
  REAL4          rsqvetoPow;		/**< UNDOCUMENTED */
}
FindChirpFilterOutputVetoParams;

/**
 * This structure provides the parameters used by the <tt>FindChirpFilterSegment()</tt> function.
 *
 * \note Original laldoc contained the following documentation for a non-existent struct-member:
 * <tt>REAL4 norm</tt> On exit this contains the normalisation constant
 * that relates the quantity \f$|q_j|^2\f$ with the signal to noise squared,
 * \f$\rho^2(t_j)\f$ by
 * \f{equation}{
 * \rho^2(t_j) = \textrm{norm} \times \left|q_j\right|^2
 * \f}
 */
typedef struct
tagFindChirpFilterParams
{
  REAL8                         deltaT;			/**< The sampling interval \f$\Delta t\f$ */
  REAL4                         clusterWindow;		/**< UNDOCUMENTED */
  REAL4                         rhosqThresh;		/**< The signal-to-noise ratio squared threshold
                                                         * \f$\rho^2_\ast\f$; If the matched filter output exceeds this value, that is
                                                         * \f$\rho^2(t_j) > \rho^2_\ast\f$, the event processing algorithm is entered and
                                                         * triggers may be generated (subject to addition vetoes such as the \f$\chi^2\f$
                                                         * veto); The value of \f$\rho^2_\ast0\f$ must be greater than or equal to zero
                                                         */
  REAL4                         chisqThresh;		/**< The \f$\chi^2\f$ veto threshold on; This threshold is described in details in the documentation for the \f$\chi^2\f$ veto */
  REAL4                         chisqDelta;		/**< UNDOCUMENTED */
  UINT4                         maximiseOverChirp;	/**< If not zero, use the maximise over chirp length algorithm to decide which time \f$t_j\f$ should
                                                         * have an inspiral trigger generated; Otherwise record all points that pass the
                                                         * \f$\rho^2\f$ and \f$\chi^2\f$ threshold as triggers (this may generate may triggers)
                                                         */
  UINT4                         ignoreIndex;		/**< UNDOCUMENTED */
  FindChirpClustering           clusterMethod;		/**< UNDOCUMENTED */
  Approximant                   approximant;		/**< Filter the data using templates of type \c approximant; Valid approximants are #TaylorT1, #TaylorT2,
                                                         * #TaylorT3, #PadeT1, #EOB, #FindChirpSP, #BCV and #BCVSpin; The value of
                                                         * \c approximant here must match that in the findchirp data segment and
                                                         * findchirp template used as input
                                                         */
  LALPNOrder                    order;			/**< UNDOCUMENTED */
  COMPLEX8Vector               *qVec;			/**< Pointer to vector of length \f$N\f$ allocated by <tt>FindChirpFilterInit()</tt> to store the quantity \f$q_j\f$;
                                                         * The pointer must not be NULL on entry, but the vetor may contain garbage which will be overwritten with the value
                                                         * of \f$q_j\f$ for the segment filtered on exit */
  COMPLEX8Vector               *qVecBCV;		/**< Pointer to the additional vector required for the BCV templates, allocated by <tt>FindChirpFilterInit()</tt> */
  COMPLEX8Vector               *qVecBCVSpin1;		/**< Pointer to the additional vector required for filtering spinning BCV templates, allocated by <tt>FindChirpFilterInit()</tt> */
  COMPLEX8Vector               *qVecBCVSpin2;		/**< Pointer to the additional vector required for filtering spinning BCV templates, allocated by <tt>FindChirpFilterInit()</tt> */
  COMPLEX8Vector               *qtildeVec;		/**< Pointer to vector of length \f$N\f$ allocated by <tt>FindChirpFilterInit()</tt> to store the quantity
                                                         * \f$\tilde{q}_k\f$, given by
                                                         * \f{equation}{
                                                         * \tilde{q}_k = \left\{
                                                         * \begin{array}{ll}
                                                         * \tilde{F}_k \tilde{T}_k^\ast & \quad 0 < k < \frac{N}{2} \\,
                                                         * 0 & \quad \textrm{otherwise}
                                                         * \end{array}
                                                         * \right.
                                                         * \f}
                                                         * The pointer must not be NULL on entry, but the vetor may contain garbage which
                                                         * will be overwritten with the value of \f$\tilde{q}_k\f$ for the segment filtered
                                                         * on exit
                                                         */
  COMPLEX8Vector               *qtildeVecBCV;		/**< Pointer to the additional vector required for filtering BCV templates, allocated by <tt>FindChirpFilterInit()</tt> */
  COMPLEX8Vector               *qtildeVecBCVSpin1;	/**< Pointer to the additional vector required for filtering spinning BCV templates, allocated by <tt>FindChirpFilterInit()</tt> */
  COMPLEX8Vector               *qtildeVecBCVSpin2;	/**< Pointer to the additional vector required for filtering spinning BCV templates, allocated by <tt>FindChirpFilterInit()</tt> */
  COMPLEX8Vector              **qVecACTD;		/**< UNDOCUMENTED */
  COMPLEX8Vector              **qtildeVecACTD;		/**< UNDOCUMENTED */
  COMPLEX8VectorSequence       *PTFqVec;		/**< UNDOCUMENTED */
  COMPLEX8Vector               *PTFsnrVec;		/**< UNDOCUMENTED */
  REAL4Array                   *PTFA;			/**< UNDOCUMENTED */
  REAL4Array                   *PTFMatrix;		/**< UNDOCUMENTED */
  ComplexFFTPlan               *invPlan;		/**< Pointer to FFTW plan created by <tt>FindChirpFilterInit()</tt> to transform the quantity \f$\tilde{q}_k\f$ to
                                                         * \f${q}_j\f$ usimg the inverse DFT; Must not be NULL
                                                         */
  REAL4TimeSeries              *rhosqVec;		/**< Pointer to a time series which contains a vector of length \f$N\f$; If this is not NULL, the complex filter
                                                         * output \f$\rho(t_j) = x(t_j) + iy(t_j)\f$ is stored in the vector; This quantity can be used by the coherent filtering code
                                                         */
  COMPLEX8TimeSeries           *cVec;			/**< UNDOCUMENTED */
  REAL4Vector                  *chisqVec;		/**< Workspace vector of length \f$N\f$ used to compute and store \f$\chi^2(t_j)\f$; Must not be NULL if \c numChisqBins is
                                                         * greater than zero; Contains \f$\chi^2(t_j)\f$ on exit
                                                         */
  FindChirpFilterOutputVetoParams *filterOutputVetoParams; /**< Pointer to the parameter structure for the additional signal based veto function */
}
FindChirpFilterParams;

/* ---------- typedefs of input structures used by functions in findchirp ---------- */

/**
 * UNDOCUMENTED
 */
typedef struct
tagFindChirpBankSimParams
{
  Approximant           approx;		/**< Waveform pproximant to use for injection */
  LALPNOrder            order;		/**< Waveform order to use for injection */
  REAL4                 minMass;	/**< Minimum mass of injected signals */
  REAL4                 maxMass;	/**< Maximum mass of injected signals  */
  RandomParams         *randParams;	/**< UNDOCUMENTED */
  INT4                  maxMatch;	/**< UNDOCUMENTED */
  CHAR                 *frameName;	/**< UNDOCUMENTED */
  CHAR                 *frameChan;	/**< UNDOCUMENTED */
  REAL4			f_lower;	/**< UNDOCUMENTED */
}
FindChirpBankSimParams;

/*@}*/


/*
 *
 * function prototypes for initialization, finalization and filter functions
 *
 */

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
LALFindChirpInjectSignals (
    LALStatus                  *status,
    REAL4TimeSeries            *chan,
    SimInspiralTable           *events,
    COMPLEX8FrequencySeries    *resp
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

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _FINDCHIRPH_H */
