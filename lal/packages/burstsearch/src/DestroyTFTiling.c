/******** <lalVerbatim file="DestroyTFTilingCV"> ********
Author: Flanagan, E
$Id$
********* </lalVerbatim> ********/
 

#include <lal/LALRCSID.h>


NRCSID (DESTROYTFTILINGC, "$Id$");


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <lal/ExcessPower.h>
#include <lal/LALConstants.h>
#include <lal/LALErrno.h>
#include <lal/LALStdlib.h>
#include <lal/Random.h>
#include <lal/RealFFT.h>
#include <lal/SeqFactories.h>
#include <lal/Thresholds.h>


#define TRUE 1
#define FALSE 0


extern INT4 lalDebugLevel;


static void 
DestroyTFTile (LALStatus *status, TFTile *tfTile)
{
  /*
   *  this function destroys linked list of tiles
   *  The easiest coding method would be to use a recursive function,
   *  but that would use a lot of memory on the stack, and memory
   *  overflows on stack are hard to trap.  So, code it by copying out the
   *  linked list of pointers into an array of pointers.
   */

  INT4 tileCount;
  INT4 i;
  TFTile *thisTile;
  TFTile **tiles;

  INITSTATUS (status, "DestroyTFTile", DESTROYTFTILINGC);
  ATTATCHSTATUSPTR (status);

  ASSERT(tfTile, status, LAL_NULL_ERR, LAL_NULL_MSG);

  /* count the tiles */
  tileCount=0;
  thisTile = tfTile;
  while (thisTile != NULL)
    {
      tileCount++;
      thisTile = thisTile->nextTile;
    }

  /* allocate memory for array of pointers to tiles */
  tiles = NULL;
  tiles = (TFTile **) LALMalloc (tileCount * sizeof(TFTile *));
  
  /*  Make sure that the allocation was succesful */
  if ( !(tiles) ){
    ABORT (status, LAL_NULL_ERR, LAL_NULL_MSG);
  }

  tileCount=0;
  thisTile = tfTile;

  while (thisTile != NULL)
    {
      tileCount++;
      *(tiles + tileCount-1) = thisTile;
      thisTile = thisTile->nextTile;
    }
  
  for(i=tileCount-1; i>=0; i--)
    {
      thisTile = *(tiles + i);
      LALFree (thisTile);
    }

  LALFree (tiles);

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/******** <lalVerbatim file="DestroyTFTilingCP"> ********/
void
LALDestroyTFTiling (
		 LALStatus                             *status,
		 TFTiling                           **tfTiling
		 )
/******** </lalVerbatim> ********/
{
  INT4                         i;

  INITSTATUS (status, "LALDestroyTFTiling", DESTROYTFTILINGC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not null */
  ASSERT(tfTiling, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(*tfTiling, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT((*tfTiling)->tfp, status, LAL_NULL_ERR, LAL_NULL_MSG); 
  ASSERT ((*tfTiling)->firstTile, status, LAL_NULL_ERR, LAL_NULL_MSG);

  /* make sure that number of TF planes is positive */
  ASSERT((*tfTiling)->numPlanes>0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);

  /* destroy the set of time frequency planes and DFTParams */
  for(i=0;i<(*tfTiling)->numPlanes;i++)
    {
      COMPLEX8TimeFrequencyPlane    **thisPlane;
      ComplexDFTParams              **thisDftParams;

      thisPlane = (*tfTiling)->tfp + i;
      LALDestroyTFPlane( status->statusPtr, thisPlane );
      CHECKSTATUSPTR (status);

      thisDftParams = (*tfTiling)->dftParams + i;
      LALDestroyComplexDFTParams( status->statusPtr, thisDftParams );
      CHECKSTATUSPTR (status);
    }

  /* free the vector of pointers to TF planes */
  LALFree ( (*tfTiling)->tfp );

  /* free the vector of pointers to DFTParams */
  LALFree ( (*tfTiling)->dftParams );

  /* destroy the linked list of TF Tiles */
  DestroyTFTile (status->statusPtr, (*tfTiling)->firstTile);
  CHECKSTATUSPTR (status);

  /* free tfTiling struct itself */
  LALFree ( *tfTiling );           

  *tfTiling = NULL;	    /* make sure we don't point to freed struct */ 

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}
