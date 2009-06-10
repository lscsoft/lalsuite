/*
 * Copyright (C) 2006 Reinhard Prix
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

/** \author R. Prix
 * \file
 * \brief
 * Provides empty initializers for pulsar struct data-types.
 *
 */

/*---------- INCLUDES ----------*/
#include "PulsarDataTypes.h"
#include "ComputeFstat.h"


NRCSID( PULSARDATATYPESC, "$Id$");

/*---------- local DEFINES ----------*/

/*----- Macros ----- */

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- Global variables ----------*/

/* empty struct initializers */
const BinaryOrbitParams empty_BinaryOrbitParams;
const PulsarAmplitudeParams empty_PulsarAmplitudeParams;
const PulsarSpinRange empty_PulsarSpinRange;
const PulsarDopplerParams empty_PulsarDopplerParams;
const PulsarParams empty_PulsarParams;
const PulsarCandidate empty_PulsarCandidate;

/* initializers for frequently used non-pulsar types */
const EphemerisData empty_EphemerisData;
const LIGOTimeGPS empty_LIGOTimeGPS = {0,0};
const SkyPosition empty_SkyPosition = {0,0,COORDINATESYSTEM_HORIZON};
const BarycenterInput empty_BarycenterInput = { {0,0},
						{ {0,0,0},
						  {{0,0,0},{0,0,0},{0,0,0}},
						  LALDETECTORTYPE_ABSENT,
						  {{0},{0},0,0,0,0,0,0,0,0,0}
						},0,0,0};
/*---------- internal prototypes ----------*/

/*==================== FUNCTION DEFINITIONS ====================*/
