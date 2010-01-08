/*
 * Copyright (C) 2007 Reinhard Prix
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

/**
 * \author Reinhard Prix
 * \date 2006
 * \file
 * \brief Header file defining the Doppler lattice-covering
 *
 * $Id$
 */

#ifndef _DOPPLERLATTICECOVERING_H  /* Double-include protection. */
#define _DOPPLERLATTICECOVERING_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/*---------- INCLUDES ----------*/
NRCSID( DOPPLERFLATTICECOVERINGH, "$Id$" );

/*---------- DEFINES ----------*/

/*---------- external types ----------*/
typedef struct {
  DopplerRegion searchRegion;		/**< Doppler-space region to be covered + scanned */
  REAL8 metricMismatch;			/**< for GRID_METRIC and GRID_ISOTROPIC */
  LIGOTimeGPS startTime;		/**< start-time of the observation */
  REAL8 Tspan;				/**< total time spanned by observation */
  const EphemerisData *ephemeris;	/**< ephemeris-data for numerical metrics */
} DopplerLatticeInit;

/** Opaque type holding current state of lattice-scan */
typedef struct tagDopplerLatticeScan DopplerLatticeScan;

/*---------- Global variables ----------*/
/* some empty structs for initializations */

/*---------- external prototypes [API] ----------*/
void InitDopplerLatticeScan(LALStatus *, DopplerLatticeScan **scan, const DopplerLatticeInit *init );
int XLALFreeDopplerLatticeScan ( DopplerLatticeScan **scan );

int XLALgetCurrentLatticeIndex ( gsl_vector_int **lal_index, const DopplerLatticeScan *scan  );
int XLALsetCurrentLatticeIndex ( DopplerLatticeScan *scan, const gsl_vector_int *lal_index );
int XLALgetCurrentDopplerPos ( PulsarDopplerParams *pos, const DopplerLatticeScan *scan, CoordinateSystem skyCoords );
int XLALadvanceLatticeIndex ( DopplerLatticeScan *scan );
REAL8 XLALCountLatticeTemplates ( const DopplerLatticeScan *scan );

/* ----- variout utility functions ----- */

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
