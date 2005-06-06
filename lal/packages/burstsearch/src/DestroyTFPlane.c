/******** <lalVerbatim file="DestroyTFPlaneCV"> ********
Author: Flanagan, E
$Id$
********* </lalVerbatim> ********/

#include <lal/LALRCSID.h>

NRCSID (DESTROYTFPLANEC, "$Id$");

#include <lal/LALStdlib.h>
#include <lal/TFTransform.h>

/******** <lalVerbatim file="DestroyTFPlaneCP"> ********/
void
XLALDestroyTFPlane(
	COMPLEX8TimeFrequencyPlane *tfp
)
/******** </lalVerbatim> ********/
{
	if(tfp) {
		LALFree(tfp->params);
		LALFree(tfp->data);
		LALFree(tfp);
	}
}
