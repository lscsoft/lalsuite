


#ifndef _SKYANNULUS_H  /* Double-include protection. */
#define _SKYANNULUA_H


#include <lal/DetectorSite.h>
#include <lal/SkyCoordinates.h>
#include <lal/Date.h>
#include <lal/TimeDelay.h>

NRCSID( COMPUTESKYPOSH, "$Id$" );



void
LALComputeSingleAnnulus ( 
		   LALStatus          *status,
		   INT2Vector         *sky,
		   const LALPlaceAndGPS     *placetime1,
		   const LALPlaceAndGPS     *placetime2,
		   REAL4              dtErr,
		   INT2               nx,
		   INT2               ny
		   );

void
LALComputeTripleAnnuli ( 
		       LALStatus          *status,
		       INT2Vector         **sky,
		       LALPlaceAndGPS     *placetime1,
		       LALPlaceAndGPS     *placetime2,
		       LALPlaceAndGPS     *placetime3,
		       REAL4              dtErr,
		       INT2               nx,
		       INT2               ny
		       );

void
LALExistTripleAnnuli( 
		     LALStatus          *status,
		     INT2               *exist,
		     LALPlaceAndGPS     *placetime1,
		     LALPlaceAndGPS     *placetime2,
		     LALPlaceAndGPS     *placetime3,
		     REAL4              dtErr,
		     INT2               nx,
		     INT2               ny
		     );



#endif
