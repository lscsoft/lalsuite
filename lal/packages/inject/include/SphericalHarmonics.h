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

/** Spherical Harmonic for the l=2 mode */
COMPLEX16 SphHarm ( 
    UINT4   L,      /**< value of L */
    INT4    M,      /**< value of M */
    REAL4   theta,  /**< angle with respect to the z axis */
    REAL4   phi     /**< angle with respect to the x axis */)

{
    COMPLEX16  out; /* complex number */
    REAL4      deptheta; /** dependency on theta */

    if (L == 2)
    {
	switch ( M )
	{
	    case -2:
		deptheta = sqrt( 5.0 / ( 64.0 * LAL_PI ) ) * ( 1.0 - cos( theta ))*( 1.0 - cos( theta ));
		out.re = deptheta * cos( -2.0*phi );
		out.im = deptheta * sin( -2.0*phi );
		break;

	    case -1:
		deptheta = sqrt( 5.0 / ( 16.0 * LAL_PI ) ) * sin( theta )*( 1.0 - cos( theta ));
		out.re = deptheta * cos( -phi );
		out.im = deptheta * sin( -phi );
		break;

	    case 0:
		deptheta = sqrt( 15.0 / ( 32.0 * LAL_PI ) ) * sin( theta )*sin( theta );
		out.re = deptheta;
		out.im = deptheta;
		break;

	    case 1:
		deptheta = sqrt( 5.0 / ( 16.0 * LAL_PI ) ) * sin( theta )*( 1.0 + cos( theta ));
		out.re = deptheta * cos( phi );
		out.im = deptheta * sin( phi );
		break;
		
	    case 2:
		deptheta = sqrt( 5.0 / ( 64.0 * LAL_PI ) ) * ( 1.0 + cos( theta ))*( 1.0 + cos( theta ));
		out.re = deptheta * cos( 2.0*phi );
		out.im = deptheta * sin( 2.0*phi );
		break;	   
	    
	    default:
		/* Error message informing that the chosen M is incompatible with L=2*/
		printf("Sorry, the value chosen for m is not compatible with l");
		break;
	}
    }
    else 
    {
	/* Error message informing that L!=2 is not yet implemented*/
	fprintf(stderr, "Sorry, for the moment we haven't implemented anything other than l=2");
    }
    
    return( out );
}

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif
