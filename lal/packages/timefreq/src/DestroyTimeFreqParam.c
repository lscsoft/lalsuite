/*----------------------------------------------------------------------- 
 * 
 * File Name: DestroyTimeFreqParam.c
 * 
 * Author: Chassande-Mottin, E.
 * 
 * Revision: $Id: 
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * DestroyTimeFreqParam
 * 
 * SYNOPSIS 
 * void DestroyTimeFreqParam ( Status *,  TimeFreqParam **param );
 * 
 * DESCRIPTION 
 * Returns to system storage allocated by CreateTimeFreqParam
 * 
 * DIAGNOSTICS 
 * param == NULL, *param == NULL, free failure
 *
 * CALLS
 * LALFree
 * 
 * NOTES
 * 
 *-----------------------------------------------------------------------
 */

#include "TimeFreq.h"

NRCSID (DESTROYTIMEFREQPARAMC, "$Id$");

void DestroyTimeFreqParam (Status *status, TimeFreqParam **param)
{
  /*  Initialize status */
  INITSTATUS (status, "DestroyTimeFreqParam", DESTROYTIMEFREQPARAMC);
      
  /* Check param: report if NULL */
  ASSERT (param != NULL, status, DESTROYTFP_ENULL, DESTROYTFP_MSGENULL); 

  /*  Check *param: report if NULL */
  ASSERT (*param != NULL, status, DESTROYTFP_ENULL, DESTROYTFP_MSGENULL);

  switch ((*param)->type) {
  case Spectrogram :
  
    SDestroyVector(status,&(*param)->windowT);
    (*param)->type = Undefined; 

    break;
  case WignerVille :
 
    (*param)->type = Undefined; 

    break;
  case PSWignerVille :

    SDestroyVector(status,&(*param)->windowT);
    SDestroyVector(status,&(*param)->windowF);
    (*param)->type = Undefined; 
    
    break;
  case RSpectrogram :

    SDestroyVector(status,&(*param)->windowT);
    (*param)->type = Undefined; 

    break;
  default :
    ABORT(status,DESTROYTFP_ETYPE, DESTROYTFP_MSGETYPE);
    
  }
  /* Ok, now let's free allocated storage */

  LALFree ( *param );	      /* free param struct itself */
  *param = NULL;	      /* make sure we don't point to freed struct */

  RETURN (status);
}
