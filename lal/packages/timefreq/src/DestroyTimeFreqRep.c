/*----------------------------------------------------------------------- 
 * 
 * File Name: DestroyTimeFreqRep.c
 * 
 * Author: Chassande-Mottin, E.
 * 
 * Revision: $Id: 
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * DestroyTimeFreqRep
 * 
 * SYNOPSIS 
 * void LALDestroyTimeFreqRep ( LALStatus *,  TimeFreqRep **tfr );
 * 
 * DESCRIPTION 
 * Returns to system storage allocated by CreateTimeFreqRep
 * 
 * DIAGNOSTICS 
 * tfr == NULL, *tfr == NULL, (*tfr)->map == NULL, free failure
 *
 * CALLS
 * LALFree
 * 
 * NOTES
 * 
 *-----------------------------------------------------------------------
 */

#include "TimeFreq.h"

NRCSID (DESTROYTIMEFREQREPC, "$Id$");

void LALDestroyTimeFreqRep (LALStatus *status, TimeFreqRep **tfr)
{
  /* 
   * Initialize status
   */

  INITSTATUS (status, "LALDestroyTimeFreqRep", DESTROYTIMEFREQREPC);
      
  /* Check tfr: report if NULL */
  ASSERT (tfr != NULL, status, DESTROYTFR_ENULL, DESTROYTFR_MSGENULL); 

  /* Check tfr: report if NULL */
  ASSERT (*tfr != NULL, status, DESTROYTFR_ENULL, DESTROYTFR_MSGENULL);

  /* Check data in tfr: report if NULL */
  ASSERT ((*tfr)->timeInstant  != NULL, status, DESTROYTFR_ENULL, DESTROYTFR_MSGENULL);

  /* Check data in tfr: report if NULL */
  ASSERT ((*tfr)->freqBin  != NULL, status,DESTROYTFR_ENULL, DESTROYTFR_MSGENULL);

  /* Check data in tfr: report if NULL */
  ASSERT ((*tfr)->map != NULL, status, DESTROYTFR_ENULL, DESTROYTFR_MSGENULL);

  /* Ok, now let's free allocated storage */

  LALFree( (*tfr)->timeInstant );
  LALFree( (*tfr)->freqBin );

  {
    INT4 column;
    for (column = 0; column < (*tfr)->tCol; column++)
      LALFree ( (*tfr)->map[column] );
  }

  LALFree ( (*tfr)->map ); /* free allocated data */
  LALFree ( *tfr );	      /* free tfr struct itself */

  *tfr = NULL;		/* make sure we don't point to freed struct */

  RETURN (status);
}
