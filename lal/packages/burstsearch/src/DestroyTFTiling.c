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
	int i;

	if(!tfTiling)
		return;

	for(i = 0; i < tfTiling->numPlanes; i++) {
		XLALDestroyTFPlane(*(tfTiling->tfp + i));
		XLALDestroyComplexDFTParams(*(tfTiling->dftParams + i));
	}

	LALFree(tfTiling->tfp);
	LALFree(tfTiling->dftParams);

	while(tfTiling->firstTile) {
		tile = tfTiling->firstTile->nextTile;
		LALFree(tfTiling->firstTile);
		tfTiling->firstTile = tile;
	}

	LALFree(tfTiling);
}
