/******** <lalVerbatim file="SortTFTilingCV"> ********
Author: Eanna Flanagan
$Id$
********* </lalVerbatim> ********/


#include <lal/LALRCSID.h>

NRCSID(SORTTFTILINGC, "$Id$");

#include <lal/ExcessPower.h>
#include <lal/LALStdlib.h>
#include <lal/XLALError.h>

#define TRUE 1


static int TileCompare(const void *t1, const void *t2)
{
	const TFTile *tile1 = *(TFTile * const *) t1;
	const TFTile *tile2 = *(TFTile * const *) t2;
	if(tile1->alpha > tile2->alpha)
		return(1);
	if(tile1->alpha < tile2->alpha)
		return(-1);
	return(0);
}


/******** <lalVerbatim file="SortTFTilingCP"> ********/
int XLALSortTFTiling(
	TFTiling * tfTiling
)
/******** </lalVerbatim> ********************************/
{
	static const char *func = "XLALSortTFTiling";
	unsigned length;
	unsigned i;
	TFTile *tile;
	TFTile **array, **addpoint;

	/* make sure excess power has already been computed */
	if(!tfTiling->excessPowerComputed)
		XLAL_ERROR(func, XLAL_EDATA);

	/* count tiles */
	length = 0;
	for(tile = tfTiling->firstTile; tile; tile = tile->nextTile)
		length++;

	/* convert linked list to an array */
	array = LALMalloc(length * sizeof(*array));
	if(!array)
		XLAL_ERROR(func, XLAL_EFUNC);
	i = 0;
	for(tile = tfTiling->firstTile; tile; tile = tile->nextTile)
		array[i++] = tile;

	/* sort array */
	qsort(array, length, sizeof(*array), TileCompare);

	/* relink list according to array order */
	addpoint = &tfTiling->firstTile;
	for(i = 0; i < length; i++) {
		*addpoint = array[i];
		addpoint = &(*addpoint)->nextTile;
	}
	*addpoint = NULL;

	LALFree(array);

	/* set flag saying tiles have been sorted */
	tfTiling->tilesSorted = TRUE;

	return(0);
}
