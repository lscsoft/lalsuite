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
 * void LALDestroyTimeFreqParam ( LALStatus *,  TimeFreqParam **param );
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

#include <lal/TimeFreq.h>

NRCSID (DESTROYTIMEFREQPARAMC, "$Id$");

void LALDestroyTimeFreqParam (LALStatus *status, TimeFreqParam **param)
{
  /*  Initialize status */
  INITSTATUS (status, "LALDestroyTimeFreqParam", DESTROYTIMEFREQPARAMC);
      
  /* Check param: report if NULL */
  ASSERT (param != NULL, status, DESTROYTFP_ENULL, DESTROYTFP_MSGENULL); 

  /*  Check *param: report if NULL */
  ASSERT (*param != NULL, status, DESTROYTFP_ENULL, DESTROYTFP_MSGENULL);

  switch ((*param)->type) {
  case Spectrogram :
  
    LALSDestroyVector(status,&(*param)->windowT);
    (*param)->type = Undefined; 

    break;
  case WignerVille :
 
    (*param)->type = Undefined; 

    break;
  case PSWignerVille :

    LALSDestroyVector(status,&(*param)->windowT);
    LALSDestroyVector(status,&(*param)->windowF);
    (*param)->type = Undefined; 
    
    break;
  case RSpectrogram :

    LALSDestroyVector(status,&(*param)->windowT);
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
