/*----------------------------------------------------------------------- 
 * 
 * File Name: CreateTimeFreqParam.c
 * 
 * Author: Chassande-Mottin, E.
 * 
 * Revision: 
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * CreateTimeFreqParam
 * 
 * SYNOPSIS 
 * void CreateTimeFreqParam (Status *, TimeFreqParam **param, CreateTimeFreqIn *in);
 * 
 * DESCRIPTION 
 * Create a TimeFreqParam object. 
 * 
 * DIAGNOSTICS 
 *
 * CALLS
 * LALMalloc
 * 
 * NOTES
 * 
 *-----------------------------------------------------------------------
 */


#include "TimeFreq.h"

NRCSID (CREATETIMEFREQPARAMC, "$Id$");

void CreateTimeFreqParam (Status *status, 
			TimeFreqParam **param,
			CreateTimeFreqIn *in) 
{
  /* Initialize status */
  INITSTATUS (status, "CreateTimeFreqParam", CREATETIMEFREQPARAMC);	
  
  /* Check input structure: report if NULL */
  ASSERT (in != NULL, status, CREATETFP_ENULL, CREATETFP_MSGENULL);
  
  /* Check return structure  */
  ASSERT (param != NULL, status, CREATETFP_ENULL, CREATETFP_MSGENULL);
  ASSERT (*param == NULL, status, CREATETFP_ENNUL, CREATETFP_MSGENNUL);
  
  *param = (TimeFreqParam *) LALMalloc(sizeof(TimeFreqParam));

  /* Check Allocation */
  ASSERT (*param != NULL, status, CREATETFP_EMALL, CREATETFP_MSGEMALL);
  
  (*param)->type=Undefined; /* Undefined until storage allocated */
  (*param)->windowT = NULL; /* NULL data until allocated */ 
  (*param)->windowF = NULL; /* NULL data until allocated */

  switch (in->type) {
  case Spectrogram :

    /* Make sure the window length is a odd number */
    ASSERT (in->wlengthT%2 != 0, status, CREATETFP_EWSIZ, CREATETFP_MSGEWSIZ);
  
    SCreateVector(status,&(*param)->windowT,in->wlengthT);
    (*param)->type = Spectrogram; 

    break;
  case WignerVille :
 
    (*param)->type = WignerVille; 

    break;
  case PSWignerVille :

    /* Make sure the window length is a odd number */
    ASSERT (in->wlengthT%2 != 0, status, CREATETFP_EWSIZ, CREATETFP_MSGEWSIZ);

    /* Make sure the window length is a odd number */
    ASSERT (in->wlengthF%2 != 0, status, CREATETFP_EWSIZ, CREATETFP_MSGEWSIZ);

    SCreateVector(status,&(*param)->windowT,in->wlengthT);
    SCreateVector(status,&(*param)->windowF,in->wlengthF);
    (*param)->type = PSWignerVille; 
    
    break;
  case RSpectrogram :

    /* Make sure the window length is a odd number */
    ASSERT (in->wlengthT%2 != 0, status, CREATETFP_EWSIZ, CREATETFP_MSGEWSIZ);

    SCreateVector(status,&(*param)->windowT,in->wlengthT);
    (*param)->type = RSpectrogram; 

    break;
  default :
    ABORT(status,CREATETFP_ETYPE, CREATETFP_MSGETYPE);
    
  }
  
  /* We be done: Normal exit */

  RETURN (status);
}
