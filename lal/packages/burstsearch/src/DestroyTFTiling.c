/******** <lalVerbatim file="DestroyTFTilingCV"> ********
Author: Flanagan, E. and Cannon, K.
$Id$
********* </lalVerbatim> ********/
 

#include <lal/LALRCSID.h>

NRCSID (DESTROYTFTILINGC, "$Id$");

#include <lal/ExcessPower.h>
#include <lal/LALStdlib.h>


/******** <lalVerbatim file="DestroyTFTilingCP"> ********/
void
XLALDestroyTFTiling(
	TFTile *list
)
/******** </lalVerbatim> ********/
{
	TFTile *next;

	while(list) {
		next = list->nextTile;
		LALFree(list);
		list = next;
	}
}
