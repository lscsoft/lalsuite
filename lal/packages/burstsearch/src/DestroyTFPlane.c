/******** <lalVerbatim file="DestroyTFPlaneCV"> ********
Author: Flanagan, E. and Cannon, K.
$Id$
********* </lalVerbatim> ********/

#include <lal/LALRCSID.h>

NRCSID (DESTROYTFPLANEC, "$Id$");

#include <lal/LALStdlib.h>
#include <lal/TFTransform.h>

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
