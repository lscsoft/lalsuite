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
	*numEvents= 0;

	INITSTATUS (status, "LALCountEPEvents", COUNTEPEVENTSC);
	ATTATCHSTATUSPTR (status);

	/* make sure that arguments are not NULL */
	ASSERT(tfTiling, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);
	ASSERT(tfTiling->firstTile, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);
	ASSERT(numEvents, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);
	ASSERT(alphaThreshold > 0.0, status, EXCESSPOWERH_EPOSARG, EXCESSPOWERH_MSGEPOSARG);

	/* check that tiles are sorted */
	ASSERT(tfTiling->tilesSorted, status, EXCESSPOWERH_EORDER, EXCESSPOWERH_MSGEORDER);

	for(tile = tfTiling->firstTile; tile && (tile->alpha <= alphaThreshold); tile = tile->nextTile)
		*numEvents++;

	/* normal exit */
	DETATCHSTATUSPTR(status);
	RETURN(status);
}


/******** <lalVerbatim file="TFTilesToBurstEventsCP"> ********/
void
LALTFTileToBurstEvent(
	LALStatus       *status,
	SnglBurstTable  *event,
	TFTile          *tile,
	LIGOTimeGPS     *epoch,
	EPSearchParams  *params  
)
/******** </lalVerbatim> ********/
{
	INITSTATUS (status, "LALCountEPEvents", COUNTEPEVENTSC);
	ATTATCHSTATUSPTR(status);

	ASSERT(event, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(tile, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(epoch, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(params, status, LAL_NULL_ERR, LAL_NULL_MSG);

	event->next = NULL;
	strncpy(event->ifo, params->channelName, 2);
	event->ifo[2] = '\0';
	strncpy(event->search, "power", LIGOMETA_SEARCH_MAX);
	strncpy(event->channel, params->channelName, LIGOMETA_CHANNEL_MAX);

	event->start_time = XLALAddFloatToGPS(*epoch, tile->tstart * tile->deltaT);
	event->duration = (tile->tend - tile->tstart + 1) * tile->deltaT;
	event->peak_time = XLALAddFloatToGPS(event->start_time, 0.5 * event->duration);
	event->bandwidth = (tile->fend - tile->fstart + 1) / tile->deltaT;
	event->central_freq = params->tfTilingInput->flow + tile->fstart + (0.5 * event->bandwidth);
	event->amplitude = tile->excessPower;
	event->snr = tile->excessPower;

	/* report log(confidence) to avoid problems with REAL4's */
	if(tile->alpha < 0)
		event->confidence =  LAL_REAL4_MAX;
	else if(tile->alpha == 0)
		event->confidence = -LAL_REAL4_MAX;
	else
		event->confidence =  log(tile->alpha);

	event->event_id = NULL;

	/* normal exit */
	DETATCHSTATUSPTR(status);
	RETURN(status);
}
