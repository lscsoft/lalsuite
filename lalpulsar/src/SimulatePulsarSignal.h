/*
 * Copyright (C) 2013 Reinhard Prix
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

#ifndef _SIMULATEPULSARSIGNAL_H  /* Double-include protection. */
#define _SIMULATEPULSARSIGNAL_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \addtogroup SimulatePulsarSignal_h
 * \ingroup pkg_pulsarCommon
 * \author Reinhard Prix
 * \date 2004, 2013
 *
 * \brief New pulsar signal-generation routines
 */
/*@{*/

#include <lal/GeneratePulsarSignal.h>

/*---------- Global variables ----------*/

/*---------- Type definitions ----------*/

/* ---------- Function prototypes ---------- */
REAL4TimeSeries *XLALSimulateExactPulsarSignal ( const PulsarSignalParams *params );

/*@}*/

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
