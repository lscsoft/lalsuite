/******** <lalVerbatim file="CreateTFPlaneCV"> ********
Author: Flanagan, E. and Cannon, K.
$Id$
********* </lalVerbatim> ********/

#include <lal/LALRCSID.h>

NRCSID(CREATETFPLANEC, "$Id$");

#include <lal/LALMalloc.h>
#include <lal/LALStdlib.h>
#include <lal/TFTransform.h>
#include <lal/Sequence.h>
#include <lal/XLALError.h>


/******** <lalVerbatim file="CreateTFPlaneCP"> ********/
REAL4TimeFrequencyPlane *XLALCreateTFPlane(
	UINT4 timeBins, /* Number of time bins in TF plane */
	REAL8 deltaT,   /* time resolution of the plane */
	UINT4 channels, /* Number of frequency channels in TF plane */
	REAL8 deltaF,   /* frequency resolution of the plane */
	REAL8 flow      /* minimum frequency to search for */
)
/******** </lalVerbatim> ********/
{
	static const char func[] = "XLALCreateTFPlane";
	REAL4TimeFrequencyPlane *plane;
	REAL4Sequence *channel_mean_square;
	REAL4Sequence **channel;
	unsigned i;

	/* Make sure that input parameters are reasonable */
	if((flow < 0) ||
	   (deltaF <= 0.0) ||
	   (deltaT <= 0.0))
		XLAL_ERROR_NULL(func, XLAL_EINVAL);

	/* Allocate memory */
	plane = LALMalloc(sizeof(*plane));
	channel_mean_square = XLALCreateREAL4Sequence(channels);
	channel = LALMalloc(channels * sizeof(*channel));
	if(!plane || !channel_mean_square || !channel) {
		LALFree(plane);
		XLALDestroyREAL4Sequence(channel_mean_square);
		LALFree(channel);
		XLAL_ERROR_NULL(func, XLAL_ENOMEM);
	}
	for(i = 0; i < channels; i++) {
		channel[i] = XLALCreateREAL4Sequence(timeBins);
		if(!channel[i]) {
			while(--i)
				XLALDestroyREAL4Sequence(channel[i]);
			LALFree(plane);
			XLALDestroyREAL4Sequence(channel_mean_square);
			LALFree(channel);
		XLAL_ERROR_NULL(func, XLAL_ENOMEM);
		}
	}

	/* 
	 * Initialize the structure.
	 */

	plane->name[0] = '\0';
	plane->epoch.gpsSeconds = 0;
	plane->epoch.gpsNanoSeconds = 0;
	plane->deltaT = deltaT;
	plane->channels = channels;
	plane->deltaF = deltaF;
	plane->flow = flow;
	plane->channel_mean_square = channel_mean_square;
	plane->channel = channel;

	/* Normal exit */
	return(plane);
}


/******** <lalVerbatim file="DestroyTFPlaneCP"> ********/
void
XLALDestroyTFPlane(
	REAL4TimeFrequencyPlane *plane
)
/******** </lalVerbatim> ********/
{
	unsigned i;

	if(plane) {
		XLALDestroyREAL4Sequence(plane->channel_mean_square);
		for(i = 0; i < plane->channels; i++)
			XLALDestroyREAL4Sequence(plane->channel[i]);
		LALFree(plane->channel);
	}
	LALFree(plane);
}
