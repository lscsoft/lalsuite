/*
*  Copyright (C) 2007 Bernd Machenschalk, Duncan Brown
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

#ifndef _GENERATERING_H
#define _GENERATERING_H

#include <lal/LALStdlib.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/SkyCoordinates.h>
#include <lal/LIGOMetadataTables.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * \addtogroup GenerateRing_h
 * \author Goggin, L., and Brown, D. A.
 *
 * \brief Provides routines to generate ringdown waveforms.
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/GenerateRing.h>
 * \endcode
 *
 * Computes the ringdown waveform with specified \f$h_{rss}\f$.
 *
 */
/*@{*/

/** \name Error Codes */
/*@{*/
#define GENERATERINGH_ENUL 1	/**< Unexpected null pointer in arguments */
#define GENERATERINGH_EOUT 2	/**< Output field a, f, phi, or shift already exists */
#define GENERATERINGH_EMEM 3	/**< Out of memory */
#define GENERATERINGH_ETYP 4	/**< Waveform type not implemented */
#define GENERATERINGH_ELEN 5	/**< Waveform length not correctly specified */
/*@}*/

/** \cond DONT_DOXYGEN */
#define GENERATERINGH_MSGENUL "Unexpected null pointer in arguments"
#define GENERATERINGH_MSGEOUT "Output field a, f, phi, or shift already exists"
#define GENERATERINGH_MSGEMEM "Out of memory"
#define GENERATERINGH_MSGETYP "Waveform type not implemented"
#define GENERATERINGH_MSGELEN "Waveform length not correctly specified"
/** \endcond */

/** UNDOCUMENTED */
typedef enum{
  Ringdown
} SimRingType;


/**
 * This structure stores the parameters for constructing a burst gravitational
 * waveform
 */
typedef struct tagRingParamStruc {
  REAL8 deltaT;             /**< requested sampling interval (s) */
  CoordinateSystem system;  /**< coordinate system to assume for simRingdown */
} RingParamStruc;

/* ---------- Function prototypes. ---------- */

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/** \see See \ref GenerateRing_h for documentation */
void
LALGenerateRing(
    LALStatus          *status,
    CoherentGW         *output,
    REAL4TimeSeries    UNUSED *series,
    SimRingdownTable   *simRingdown,
    RingParamStruc     *params
    );

/** \see See \ref GenerateRing_h for documentation */
void
LALRingInjectSignals(
    LALStatus               *status,
    REAL4TimeSeries         *series,
    SimRingdownTable        *injections,
    COMPLEX8FrequencySeries *resp,
    INT4                     calType
    );

/*@}*/

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _GENERATERING_H */
