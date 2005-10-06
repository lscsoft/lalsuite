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
	/*
	 * Sort tiles in order of increasing lnalpha (log of tile
	 * probability).  Since the probability is in [0, 1], lnalpha's are
	 * negative with more negative numbers indicating the less likely
	 * tiles;  therefore the tiles are sorted with the least likely
	 * tiles --- tiles with greatest excess power --- first.
	 */

	qsort(tiling->tile, tiling->numtiles, sizeof(*tiling->tile), TileCompareByAlpha);
}
