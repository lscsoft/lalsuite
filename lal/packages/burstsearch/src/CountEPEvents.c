/********* <lalVerbatim file="CountEPEventsCV"> ********
Author: Flanagan, E and Brady, P R
$Id$
********** </lalVerbatim> ********/


#include <lal/LALRCSID.h>

NRCSID(COUNTEPEVENTSC, "$Id$");

#include <string.h>
#include <math.h>

#include <lal/BurstSearch.h>
#include <lal/Date.h>
#include <lal/EPSearch.h>
#include <lal/LALConstants.h>
#include <lal/LALErrno.h>
#include <lal/LIGOMetadataTables.h>

/******** <lalVerbatim file="CountEPEventsCP"> ********/
void
LALCountEPEvents(
               LALStatus         *status,
               INT4              *numEvents,
               TFTiling          *tfTiling,
               REAL8              alphaThreshold
)
/******** </lalVerbatim> ********/
{
	TFTile *tile;

	INITSTATUS (status, "LALCountEPEvents", COUNTEPEVENTSC);
	ATTATCHSTATUSPTR (status);

	/* make sure that arguments are not NULL */
	ASSERT(tfTiling, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);
	ASSERT(tfTiling->firstTile, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);
	ASSERT(numEvents, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);
	ASSERT(alphaThreshold > 0.0, status, EXCESSPOWERH_EPOSARG, EXCESSPOWERH_MSGEPOSARG);

	/* check that tiles are sorted */
	ASSERT(tfTiling->tilesSorted, status, EXCESSPOWERH_EORDER, EXCESSPOWERH_MSGEORDER);

	*numEvents= 0;
	for(tile = tfTiling->firstTile; tile && (tile->alpha <= alphaThreshold); tile = tile->nextTile)
		*numEvents++;

	/* normal exit */
	DETATCHSTATUSPTR(status);
	RETURN(status);
}

