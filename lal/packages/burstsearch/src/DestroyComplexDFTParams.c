/******** <lalVerbatim file="DestroyComplexDFTParamsCV"> ********
Author: Flanagan, E
$Id$
********* </lalVerbatim> ********/


#include <lal/LALRCSID.h>


NRCSID (DESTROYCOMPLEXDFTPARAMSC, "$Id$");


#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALErrno.h>
#include <lal/SeqFactories.h>
#include <lal/TFTransform.h>

/******** <lalVerbatim file="DestroyComplexDFTParamsCP"> ********/
void
LALDestroyComplexDFTParams (
		         LALStatus                 *status, 
		         ComplexDFTParams       **dftParams
		        )
/******** </lalVerbatim> ********/
{
  INITSTATUS (status, "LALDestroyComplexDFTParams", DESTROYCOMPLEXDFTPARAMSC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not null */
  ASSERT(dftParams, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(*dftParams, status, LAL_NULL_ERR, LAL_NULL_MSG);

  /* make sure that data pointed to is non-null */
  ASSERT((*dftParams)->plan, status, LAL_NULL_ERR, LAL_NULL_MSG); 
  ASSERT((*dftParams)->window, status, LAL_NULL_ERR, LAL_NULL_MSG); 

  /* Ok, now let's free allocated storage */
  LALSDestroyVector (status->statusPtr, &((*dftParams)->window));
  CHECKSTATUSPTR (status);
  LALDestroyComplexFFTPlan (status->statusPtr, &((*dftParams)->plan));
  CHECKSTATUSPTR (status);
  LALFree ( *dftParams );      /* free DFT parameters structure itself */

  *dftParams = NULL;	       /* make sure we don't point to freed struct */

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}
