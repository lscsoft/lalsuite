/*
 * Copyright (C) 2005 Reinhard Prix
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
 * \defgroup flatPulsarMetric
 * \ingroup pulsarCommon
 * \brief Implement a flat approximation to the pulsar-metric, based
 *        on the "linear model I" in \ref JK99.
 *
 */

/**
 * \author Reinhard Prix
 * \date 2005
 * \file 
 * \ingroup flatPulsarMetric
 * \brief Header-file defining the API for the flat pulsar-metric functions.
 *
 * $Id$
 *
 */

#ifndef _FLATPULSARMETRIC_H  /* Double-include protection. */
#define _FLATPULSARMETRIC_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/*---------- INCLUDES ----------*/
#include <lal/LALStdlib.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/DetectorSite.h>


NRCSID( FLATPULSARMETRICH, "$Id$" );


/*---------- DEFINES ----------*/
#define PMETRIC_MIN(x,y) ((x) < (y) ? (x) : (y))
#define PMETRIC_MAX(x,y) ((x) > (y) ? (x) : (y))

/** Translate metrix matrix-indices (a,b) into vector-index l */
#define PMETRIC_INDEX(a,b) (PMETRIC_MIN((a),(b))+PMETRIC_MAX((a),(b))*(PMETRIC_MAX((a),(b)) + 1 ) / 2 )

/*----- Error-codes -----*/
#define FLATPULSARMETRIC_ENULL 		1
#define FLATPULSARMETRIC_ENONULL	2
#define FLATPULSARMETRIC_EMEM		3
#define FLATPULSARMETRIC_EINPUT		4
#define FLATPULSARMETRIC_ELIST		5
#define FLATPULSARMETRIC_EFUNC		6

#define FLATPULSARMETRIC_MSGENULL 	"Arguments contained an unexpected null pointer"
#define FLATPULSARMETRIC_MSGENONULL	"Output pointer is not NULL"
#define FLATPULSARMETRIC_MSGEMEM	"Out of memory"
#define FLATPULSARMETRIC_MSGEINPUT	"Invald input parameter"
#define FLATPULSARMETRIC_MSGELIST	"Error occurred in list-handling ..."
#define FLATPULSARMETRIC_MSGEFUNC	"Sub-routine failed"

/*---------- exported types ----------*/

/*---------- Global variables ----------*/

/*---------- exported prototypes [API] ----------*/
void LALFlatPulsarMetric ( LALStatus *, 
			   REAL8Vector **metric,	
			   LIGOTimeGPS startTime,
			   REAL8 duration, 
			   const LALDetector *site);

#ifdef  __cplusplus
}
#endif  
/* C++ protection. */

#endif  /* Double-include protection. */
