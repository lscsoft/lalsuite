/******** <lalVerbatim file="CreateTFPlaneCV"> ********
Author: Flanagan, E. and Cannon, K.
$Id$
********* </lalVerbatim> ********/

#include <lal/LALRCSID.h>

NRCSID(CREATETFPLANEC, "$Id$");

#include <lal/LALMalloc.h>
#include <lal/LALStdlib.h>
#include <lal/TFTransform.h>
#include <lal/XLALError.h>


/******** <lalVerbatim file="CreateTFPlaneCP"> ********/
REAL4TimeFrequencyPlane *XLALCreateTFPlane(
	INT4 timeBins,  /* Number of time bins in TF plane */
	REAL8 deltaT,   /* time resolution of the plane */
	INT4 freqBins,  /* Number of freq bins in TF plane */
	REAL8 deltaF,   /* frequency resolution of the plane */
	REAL8 flow      /* minimum frequency to search for */
)
/******** </lalVerbatim> ********/
{
	static const char *func = "XLALCreateTFPlane";
	REAL4TimeFrequencyPlane *plane;
	REAL4 *data;

	/* Make sure that input parameters are reasonable */
	if((flow < 0) ||
	   (timeBins <= 0) ||
	   (freqBins <= 0) ||
	   (deltaF <= 0.0) ||
	   (deltaT <= 0.0))
		XLAL_ERROR_NULL(func, XLAL_EDATA);

	/* Allocate memory */
	plane = LALMalloc(sizeof(*plane));
	data = LALCalloc(timeBins * freqBins, sizeof(*data));
	if(!plane || !data) {
		LALFree(plane);
		LALFree(data);
		XLAL_ERROR_NULL(func, XLAL_ENOMEM);
	}

	/* 
	 * Fill some of the fields with nominal values.
	 */

	plane->name[0] = '\0';
	plane->epoch.gpsSeconds = 0;
	plane->epoch.gpsNanoSeconds = 0;
	plane->timeBins = timeBins;
	plane->deltaT = deltaT;
	plane->freqBins = freqBins;
	plane->deltaF = deltaF;
	plane->flow = flow;
	plane->data = data;

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
	if(plane)
		LALFree(plane->data);
	LALFree(plane);
}
