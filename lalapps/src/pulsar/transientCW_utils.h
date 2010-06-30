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

#include <lal/ComputeFstat.h>

/* ---------- exported API defines ---------- */

/** Struct to define parameters of a 'transient window' to be applied to obtain transient signals */
typedef enum {
  TRANSIENT_NONE = 0,
  TRANSIENT_RECTANGULAR = 1,
  TRANSIENT_EXPONENTIAL,
  TRANSIENT_LAST
} transientWindowType_t;


/* ---------- exported API types ---------- */

/** Struct defining one transient window instance */
typedef struct
{
  transientWindowType_t type;	/**< window-type: none, rectangular, exponential, .... */
  REAL8 t0;			/**< GPS start-time 't0' */
  REAL8 tau;			/**< transient timescale tau in seconds */
} transientWindow_t;

/** Struct defining a range of transient windows */
typedef struct
{
  transientWindowType_t type;	/**< window-type: none, rectangular, exponential, .... */
  UINT4 t0_min;			/**< earliest GPS start-time 't0' in seconds */
  UINT4 t0_max;			/**< latest GPS start-time 't0'  in seconds */
  UINT4 tau_min;		/**< shortest transient timescale tau in seconds */
  UINT4 tau_max;		/**< longest transient timescale tau in seconds */
} transientWindowRange_t;

/** Struct holding a transient CW candidate */
typedef struct {
  PulsarDopplerParams doppler;		/**< Doppler params of this 'candidate' */
  REAL8 fullFstat;			/**< 2F obtained in the full search over all SFTs */
  REAL8 maxFstat;			/**< maximal 2F value obtained over transientWindowRange */
  UINT4 t0_maxF;			/**< start-time of max{2F} over transientWindowRange (in GPS seconds)*/
  UINT4 tau_maxF;			/**< duration of max{2F} over transientWindowRange (in seconds) */
  REAL8 logBstat;			/**< log of Bayes-factor, marginalized over transientWindowRange */
} TransientCandidate_t;

/* empty struct initializers */
extern const TransientCandidate_t empty_TransientCandidate;

/* ---------- exported API prototypes ---------- */
int XLALApplyTransientWindow ( REAL4TimeSeries *series, transientWindow_t TransientWindowParams );

int XLALApplyTransientWindow2NoiseWeights ( MultiNoiseWeights *multiNoiseWeights,
                                            const MultiLIGOTimeGPSVector *multiTS,
                                            transientWindow_t TransientWindowParams );
int
write_TransientCandidate_to_fp ( FILE *fp, const TransientCandidate_t *thisTransCand );

/* ---------- Fstat-atoms related functions ----------*/
int XLALoutputMultiFstatAtoms ( FILE *fp, MultiFstatAtomVector *multiAtoms );
CHAR* XLALPulsarDopplerParams2String ( const PulsarDopplerParams *par );
int XLALComputeTransientBstat ( TransientCandidate_t *transientCand, const MultiFstatAtomVector *multiFstatAtoms, transientWindowRange_t windowRange );
FstatAtomVector *XLALmergeMultiFstatAtomsSorted ( const MultiFstatAtomVector *multiAtoms );

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
