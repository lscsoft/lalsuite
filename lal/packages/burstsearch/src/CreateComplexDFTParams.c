/******** <lalVerbatim file="CreateComplexDFTParamsCV"> ********
Author: Flanagan, E
$Id$
********* </lalVerbatim> ********/

#include <lal/LALRCSID.h>

NRCSID (CREATECOMPLEXDFTPARAMSC, "$Id$");

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <lal/ComplexFFT.h>
#include <lal/LALErrno.h>
#include <lal/LALStdlib.h>
#include <lal/RealFFT.h>
#include <lal/SeqFactories.h>
#include <lal/TFTransform.h>

/******** <lalVerbatim file="CreateComplexDFTParamsCP"> ********/
void
LALCreateComplexDFTParams (
                     LALStatus          *status,
                     ComplexDFTParams  **dftParams,
                     LALWindowParams    *params,
                     INT2                sign
		     )
/******** </lalVerbatim> ********/
{
  INITSTATUS (status, "LALCreateComplexDFTParams", CREATECOMPLEXDFTPARAMSC);
  ATTATCHSTATUSPTR (status);

  /* 
   * Check return structure: dftParams should point to a valid pointer
   * which should not yet point to anything.
   */
  ASSERT(dftParams, status, LAL_NULL_ERR, LAL_NULL_MSG); 
  ASSERT(*dftParams == NULL, status, LAL_NNULL_ERR, LAL_NNULL_MSG);
  ASSERT(params, status, LAL_NULL_ERR, LAL_NULL_MSG);  
  ASSERT(params->length > 0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);
  ASSERT((sign==1) || (sign==-1), status, LAL_RANGE_ERR, LAL_RANGE_MSG);


  /*  Assign memory for *dftParams   */
  *dftParams = LALMalloc(sizeof(**dftParams));
  ASSERT(*dftParams, status, LAL_NOMEM_ERR, LAL_NOMEM_MSG);

  /* fill in some values */
  (*dftParams)->window = NULL;
  (*dftParams)->plan = NULL;

  /* _Estimate_ the FFT plan */
  if(sign==1)
      LALCreateForwardComplexFFTPlan (status->statusPtr, &((*dftParams)->plan), 
                              params->length, 0);
  else
      LALCreateReverseComplexFFTPlan (status->statusPtr, &((*dftParams)->plan), 
                              params->length, 0);
  CHECKSTATUSPTR (status);

  LALSCreateVector (status->statusPtr, &((*dftParams)->window), params->length);
  CHECKSTATUSPTR (status);

  LALWindow (status->statusPtr, ((*dftParams)->window), params);
  CHECKSTATUSPTR (status);

  (*dftParams)->sumofsquares = params->sumofsquares;
  (*dftParams)->windowType = params->type;
  
  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}
