/******** <lalVerbatim file="SortTFTilingCV"> ********
Author: Eanna Flanagan, and Cannon, K.
$Id$
********* </lalVerbatim> ********/

#include <stdlib.h>
#include <lal/LALRCSID.h>

NRCSID(SORTTFTILINGC, "$Id$");

#include <lal/ExcessPower.h>
#include <lal/XLALError.h>


static int TileCompareByAlpha(const void *t1, const void *t2)
{
	const TFTile *tile1 = t1;
	const TFTile *tile2 = t2;
	if(tile1->lnalpha > tile2->lnalpha)
		return(1);
	if(tile1->lnalpha < tile2->lnalpha)
		return(-1);
	return(0);
}


/******** <lalVerbatim file="SortTFTilingCP"> ********/
void XLALSortTFTilingByAlpha(
	TFTiling *tiling
)
/******** </lalVerbatim> ********************************/
{
	qsort(tiling->tile, tiling->numtiles, sizeof(*tiling->tile), TileCompareByAlpha);
}
