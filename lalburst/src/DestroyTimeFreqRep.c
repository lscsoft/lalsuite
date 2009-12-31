 /*
  * Copyright (C) 2004, 2005 Cristina V. Torres
  *                          E. Chassande-Mottin
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 2 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with with program; see the file COPYING. If not, write to the
  *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
  *  MA  02111-1307  USA
  */
/*-----------------------------------------------------------------------
 *
 * File Name: DestroyTimeFreqRep.c
 *
 * Maintainer: Torres, C (Univ TX at Browsville)
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

#include <lal/TimeFreq.h>

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
