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
  ASSERT (dftParams, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP); 
  ASSERT (*dftParams == NULL, status, TFTRANSFORM_EALLOCP, 
          TFTRANSFORM_MSGEALLOCP);

  ASSERT (params, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);  

  ASSERT (params->length > 0, status, TFTRANSFORM_EPOSARG, 
          TFTRANSFORM_MSGEPOSARG);

  ASSERT( (sign==1) || (sign==-1), status, TFTRANSFORM_EINCOMP,
          TFTRANSFORM_MSGEINCOMP);


  /*  Assign memory for *dftParams   */
  *dftParams = (ComplexDFTParams *) LALMalloc(sizeof(ComplexDFTParams));

  /*  Make sure that the allocation was succesful */
  ASSERT (*dftParams, status, TFTRANSFORM_EMALLOC, TFTRANSFORM_MSGEMALLOC);


  /* fill in some values */
  (*dftParams)->window = NULL;
  (*dftParams)->plan = NULL;

  if(sign==1)
    {
      LALEstimateFwdComplexFFTPlan (status->statusPtr, &((*dftParams)->plan), 
                              params->length);
    }
  else
    {
      LALEstimateInvComplexFFTPlan (status->statusPtr, &((*dftParams)->plan), 
                              params->length);
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


