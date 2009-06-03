/*
*  Copyright (C) 2007 Jolien Creighton, Julien Sylvestre
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

/*************************** <lalVerbatim file="CoherentEstimation">
Author: Sylvestre, J.
$Id$
**************************************************** </lalVerbatim> */

#ifndef _COHERENTESTIMATION_H
#define _COHERENTESTIMATION_H

#include <lal/LALStdlib.h>
#include <lal/DetectorSite.h>
#include <lal/SkyCoordinates.h>
#include <lal/IIRFilter.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

NRCSID( COHERENTESTIMATION, "$Id$" );


#define COHERENTESTIMATIONH_EMEM 1
#define COHERENTESTIMATIONH_ENULL 2
#define COHERENTESTIMATIONH_E0DEC 3
#define COHERENTESTIMATIONH_EDST 4
#define COHERENTESTIMATIONH_EICE 5
#define COHERENTESTIMATIONH_ENUM 6
#define COHERENTESTIMATIONH_EUIMP 7

#define COHERENTESTIMATIONH_MSGEMEM "Memory allocation error"
#define COHERENTESTIMATIONH_MSGENULL "NULL pointer"
#define COHERENTESTIMATIONH_MSGE0DEC "invalid DetectorsData structure"
#define COHERENTESTIMATIONH_MSGEDST "input time series don't all have same start time"
#define COHERENTESTIMATIONH_MSGEICE "invalid CoherentEstimation structure"
#define COHERENTESTIMATIONH_MSGENUM "Numerical erorr"
#define COHERENTESTIMATIONH_MSGEUIMP "Implemented only for 3 detectors"

typedef struct tagDetectorsData {

  UINT4 Ndetectors;      /* number of detectors */
  REAL4TimeSeries *data; /* data time series from all detectors */

} DetectorsData;


typedef struct tagCoherentEstimation {

  UINT4 Ndetectors;        /* number of detectors */
  LALDetector *detectors;  /* vector of detectors info */
  REAL8IIRFilter **filters; /* vector of pre-processing filters */

  BOOLEAN preProcessed;    /* set to 0 to for pre-processing */
  UINT2 nPreProcessed;     /* number of times to apply pre-proc filters */

  SkyPosition *position;   /* position of source (equatorial celestial coordinates) */
  REAL8 polAngle;          /* polarization angle: counter-clockwise angle x-axis makes with a line per- pendicular to meridian of source in Westward direction (i.e. North of West), in decimal radians. */

  REAL8 plus2cross; /* ratio |s+|/|sx|, non-zero */

  REAL8 plusDotcross; /* s+ sx / |s+| |sx| */

  REAL8 **CMat; /* correlation matrix */

} CoherentEstimation;



void
LALCoherentEstimation (
		       LALStatus          *status,
		       REAL4TimeSeries *output,
		       CoherentEstimation *params,
		       DetectorsData      *in
	              );

void
LALClearCoherentData (
		      LALStatus     *status,
		      DetectorsData *dat
		      );

void
LALClearCoherentInfo (
		      LALStatus     *status,
		      CoherentEstimation *dat
		      );

#endif
