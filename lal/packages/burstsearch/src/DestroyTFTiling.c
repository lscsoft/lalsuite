/******** <lalVerbatim file="DestroyTFTilingCV"> ********
Author: Flanagan, E
$Id$
********* </lalVerbatim> ********/
 

#include <lal/LALRCSID.h>

NRCSID (DESTROYTFTILINGC, "$Id$");

#include <lal/ExcessPower.h>
#include <lal/LALStdlib.h>


/******** <lalVerbatim file="DestroyTFTilingCP"> ********/
void
XLALDestroyTFTiling(
	TFTiling *tfTiling
)
/******** </lalVerbatim> ********/
{
	TFTile *tile;

	if(!tfTiling)
		return;

	XLALDestroyTFPlane(tfTiling->tfp);

	while(tfTiling->firstTile) {
		tile = tfTiling->firstTile->nextTile;
		LALFree(tfTiling->firstTile);
		tfTiling->firstTile = tile;
	}

	LALFree(tfTiling);
}
