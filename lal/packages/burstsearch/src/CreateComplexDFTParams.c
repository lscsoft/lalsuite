/******** <lalVerbatim file="CreateComplexDFTParamsCV"> ********
Author: Flanagan, E
$Id$
********* </lalVerbatim> ********/


#include <lal/LALRCSID.h>


NRCSID (CREATECOMPLEXDFTPARAMSC, "$Id$");


#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>
#include <lal/ComplexFFT.h>
#include <lal/TFTransform.h>

/******** <lalVerbatim file="CreateComplexDFTParamsCP"> ********/
void
LALCreateComplexDFTParams ( 
                     LALStatus                         *status, 
                     ComplexDFTParams               **dftParams, 
                     LALWindowParams                   *params,
                     INT2                           sign
		     )
/******** </lalVerbatim> ********/
{


  INITSTATUS (status, "LALCreateComplexDFTParams", CREATECOMPLEXDFTPARAMSC);
  ATTATCHSTATUSPTR (status);

  /* 
   * Check return structure: dftParams should point to a valid pointer
   * which should not yet point to anything.
   *
   */
  ASSERT (dftParams, status, TFTRANSFORMH_ENULLP, TFTRANSFORMH_MSGENULLP); 
  ASSERT (*dftParams == NULL, status, TFTRANSFORMH_EALLOCP, 
          TFTRANSFORMH_MSGEALLOCP);

  ASSERT (params, status, TFTRANSFORMH_ENULLP, TFTRANSFORMH_MSGENULLP);  

  ASSERT (params->length > 0, status, TFTRANSFORMH_EPOSARG, 
          TFTRANSFORMH_MSGEPOSARG);

  ASSERT( (sign==1) || (sign==-1), status, TFTRANSFORMH_EINCOMP,
          TFTRANSFORMH_MSGEINCOMP);


  /*  Assign memory for *dftParams   */
  if ( !( *dftParams = (ComplexDFTParams *) LALMalloc(sizeof(ComplexDFTParams)) ) ){
    ABORT (status, TFTRANSFORMH_EMALLOC, TFTRANSFORMH_MSGEMALLOC);
  }

  /* fill in some values */
  (*dftParams)->window = NULL;
  (*dftParams)->plan = NULL;

  if(sign==1)
    {
      /* _Estimate_ the FFT plan */
      LALCreateForwardComplexFFTPlan (status->statusPtr, &((*dftParams)->plan), 
                              params->length, 0);
    }
  else
    {
      /* _Estimate_ the FFT plan */
      LALCreateReverseComplexFFTPlan (status->statusPtr, &((*dftParams)->plan), 
                              params->length, 0);
    }
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


