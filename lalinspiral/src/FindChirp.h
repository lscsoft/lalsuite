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
#include <lal/LIGOMetadataTables.h>

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
/** @{ */

/**\name Error Codes */
/** @{ */
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
/** @} */

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

/* ---------- typedefs of input structures used by functions in findchirp ---------- */

/** @} */

/*
 *
 * function prototypes for initialization, finalization and filter functions
 *
 */

void
LALFindChirpInjectSignals (
    LALStatus                  *status,
    REAL4TimeSeries            *chan,
    SimInspiralTable           *events,
    COMPLEX8FrequencySeries    *resp
    );

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _FINDCHIRPH_H */
