/*
*  Copyright (C) 2007 Stephen Fairhurst
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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

#ifndef _GENERATEINSPRING_H
#define _GENERATEINSPRING_H

/* includes */
#include <lal/SimulateCoherentGW.h>
#include <lal/LIGOMetadataTables.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

/**
 * \defgroup GenerateInspRing_h Header GenerateInspRing.h
 * \ingroup lalinspiral_inject
 * \author S.Fairhurst
 *
 * \brief Module for pasting a (realistic) ringdown on the end of an inspiral
 *
 */
/** @{ */

CoherentGW *
XLALGenerateInspRing(
    CoherentGW        *waveform,
    SimInspiralTable  *thisEvent,
    SimRingdownTable  *thisRingEvent,
    int                injectSignalType
    );

/** @} */

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif /* _GENERATEINSPRING_H */
