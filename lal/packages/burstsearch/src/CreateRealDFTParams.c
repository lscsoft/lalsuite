/******** <lalVerbatim file="CreateRealDFTParamsCV"> ********
Author: Flanagan, E
$Id$
********* </lalVerbatim> ********/

#include <lal/LALRCSID.h>


NRCSID (CREATEREALDFTPARAMSC, "$Id$");


#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>
#include <lal/ComplexFFT.h>
#include <lal/TFTransform.h>

/******** <lalVerbatim file="CreateRealDFTParamsCP"> ********/
void
LALCreateRealDFTParams ( 
                     LALStatus                         *status, 
                     RealDFTParams                  **dftParams, 
                     LALWindowParams                   *params,
                     INT2                           sign
		     )
/******** </lalVerbatim> ********/
{
  INITSTATUS (status, "LALCreateRealDFTParams", CREATEREALDFTPARAMSC);
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

  /*  Assign memory for *dftParams and check allocation */
  if ( !( *dftParams = (RealDFTParams *) LALMalloc(sizeof(RealDFTParams)) ) ){
    ABORT (status, TFTRANSFORMH_EMALLOC, TFTRANSFORMH_MSGEMALLOC);
  }
  
  /* fill in some values */
  (*dftParams)->window = NULL;
  (*dftParams)->plan = NULL;

  if(sign==1)
    {
      /* Estimate the FFT plan */
      LALCreateForwardRealFFTPlan (status->statusPtr, &((*dftParams)->plan), 
                              params->length, 0);
    }
  else
    {
      /* Estimate the FFT plan */
      LALCreateReverseRealFFTPlan (status->statusPtr, &((*dftParams)->plan), 
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



