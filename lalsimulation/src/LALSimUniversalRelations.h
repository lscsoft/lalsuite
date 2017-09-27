/*
 *  Copyright (C) 2017 Andrea Taracchini
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

#ifndef LALSimUniversalRelations_h
#define LALSimUniversalRelations_h

#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>
#include <math.h>

REAL8 XLALSimUniversalRelation( REAL8 x, REAL8 coeffs[] );
REAL8 XLALSimUniversalRelationlambda3TidalVSlambda2Tidal( REAL8 lambda2Tidal );
REAL8 XLALSimUniversalRelationomega02TidalVSlambda2Tidal( REAL8 lambda2Tidal );
REAL8 XLALSimUniversalRelationomega03TidalVSlambda3Tidal( REAL8 lambda3Tidal );

#endif /* LALSimUniversalRelations_h */
