/*
 * Copyright (C) 2007 S.Fairhurst, B. Krishnan, L.Santamaria
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


/** \defgroup SphericalHarmonics
 * \ingroup support
 * \author S.Fairhurst, B. Krishnan, L.Santamaria
 *
 * \brief Library of spin-weighted Spherical Harmonic functions
 *

 *
 */



/* includes */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/NRWaveIO.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

NRCSID( SPHERICALHARMONICSH, "$Id$");

/** Spin weighted Spherical Harmonic  */
INT4 XLALSphHarm ( COMPLEX16 *out, UINT4   L,  INT4 M,  REAL4 theta, REAL4   phi );



#ifdef  __cplusplus
}                /* Close C++ protection */
#endif
