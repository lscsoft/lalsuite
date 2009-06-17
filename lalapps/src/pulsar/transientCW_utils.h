/*
 * Copyright (C) 2009 Reinhard Prix
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

/*********************************************************************************/
/** \author R. Prix
 * \file
 * \brief
 * Some helper functions useful for "transient CWs", mostly applying transient window
 * functions.
 *
 *********************************************************************************/

#ifndef _COMPUTEFSTATISTIC_H
#define _COMPUTEFSTATISTIC_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

#include "config.h"

/* ---------- System includes ---------- */

/* LAL-includes */
#include <lal/AVFactories.h>
#include <lal/LogPrintf.h>
#include <lal/SFTutils.h>

/* ---------- exported API defines ---------- */

/** Struct to define parameters of a 'transient window' to be applied to obtain transient signals */
typedef enum {
  TRANSIENT_NONE = 0,
  TRANSIENT_RECTANGULAR = 1,
  TRANSIENT_EXPONENTIAL,
  TRANSIENT_LAST
} transientWindowType_t;


/* ---------- exported API types ---------- */
typedef struct
{
  transientWindowType_t type;	/**< window-type: none, rectangular, exponential, .... */
  REAL8 t0;			/**< GPS start-time 't0' */
  REAL8 tau;			/**< transient timescale tau in seconds */
} transientWindow_t;




/* ---------- exported API prototypes ---------- */
int XLALApplyTransientWindow ( REAL4TimeSeries *series, transientWindow_t TransientWindowParams );

int XLALApplyTransientWindow2NoiseWeights ( MultiNoiseWeights *multiNoiseWeights,
                                            const MultiLIGOTimeGPSVector *multiTS,
                                            transientWindow_t TransientWindowParams );


#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
