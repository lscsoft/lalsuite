/******** <lalVerbatim file="SortTFTilingCV"> ********
Author: Eanna Flanagan, and Cannon, K.
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
	TFTile **list
)
/******** </lalVerbatim> ********************************/
{
	static const char *func = "XLALSortTFTiling";
	unsigned length;
	unsigned i;
	TFTile *tile;
	TFTile **array, **addpoint;

	/* count tiles */
	length = 0;
	for(tile = *list; tile; tile = tile->nextTile)
		length++;

	/* convert linked list to an array */
	array = LALMalloc(length * sizeof(*array));
	if(!array)
		XLAL_ERROR(func, XLAL_ENOMEM);
	i = 0;
	for(tile = *list; tile; tile = tile->nextTile)
		array[i++] = tile;

	/* sort array */
	qsort(array, length, sizeof(*array), TileCompare);

	/* relink list according to array order */
	addpoint = list;
	for(i = 0; i < length; i++) {
		*addpoint = array[i];
		addpoint = &(*addpoint)->nextTile;
	}
	*addpoint = NULL;

	LALFree(array);

	return(0);
}
