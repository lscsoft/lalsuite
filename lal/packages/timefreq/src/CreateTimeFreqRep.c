/*----------------------------------------------------------------------- 
 * 
 * File Name: CreateTimeFreqRep.c
 * 
 * Author: Chassande-Mottin, E.
 * 
 * Revision: 
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * CreateTimeFreqRep
 * 
 * SYNOPSIS 
 * void CreateTimeFreqRep (Status *, TimeFreqRep **vseq, CreateTimeFreqIn *in);
 * 
 * DESCRIPTION 
 * Create a TimeFreqRep object. 
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

NRCSID (CREATETIMEFREQREPC, "$Id$");

void CreateTimeFreqRep (Status *status, 
			TimeFreqRep **tfr,
			CreateTimeFreqIn *in) 
{
  /* 
   * Initialize status
   */
  INITSTATUS (status, "CreateTimeFreqRep", CREATETIMEFREQREPC);	
  
  /* Check input structure: report if NULL */
  ASSERT (in != NULL, status, CREATETFR_ENULL, CREATETFR_MSGENULL);
  
  /* Check sequence length: report error if 0 or negative */
   ASSERT (in->fRow > 0, status, CREATETFR_EFROW, CREATETFR_MSGEFROW);

  /* Check sequence length: report error if not a power of 2 */
   /*   ASSERT (in->fRow ??, status, CREATETFR_EFROW, CREATETFR_MSGEFROW); */
  
  /* Check vector length: report error if 0 or negative*/
  ASSERT (in->tCol > 0, status, CREATETFR_ETCOL, CREATETFR_MSGETCOL); 
  
  /* Check return structure */
  ASSERT (tfr != NULL, status, CREATETFR_ENULL, CREATETFR_MSGENULL);
  ASSERT (*tfr == NULL, status, CREATETFR_ENNUL, CREATETFR_MSGENNUL);
    
  *tfr = (TimeFreqRep *) LALMalloc(sizeof(TimeFreqRep));

  /* Check allocation */
  ASSERT (*tfr != NULL, status, CREATETFR_EMALL, CREATETFR_MSGEMALL);
  
  /*   (*tfr)->TimeFreqRepType */
  (*tfr)->type=Undefined; /* Undefined until storage allocated */
  (*tfr)->fRow = 0;	/* length 0 until storage allocated */
  (*tfr)->tCol = 0;	/* length 0 until storage allocated */
  (*tfr)->timeInstant = NULL; /* NULL data until allocated */
  (*tfr)->freqBin = NULL; /* NULL data until allocated */
  (*tfr)->map = NULL;	/* NULL data until allocated */
  
  /* 
   * Allocate storage 
   */

  (*tfr)->timeInstant = (INT4 *) LALMalloc(sizeof(INT4)*in->tCol);

  if ((*tfr)->timeInstant == NULL)
  {
    /* Must free storage pointed to by *tfr */
    LALFree ((void *) *tfr);
    ABORT (status, CREATETFR_EMALL, CREATETFR_MSGEMALL);
    return;
  }

  (*tfr)->freqBin = (REAL4 *) LALMalloc(sizeof(REAL4)*(in->fRow/2+1));

  if ((*tfr)->freqBin == NULL)
  {
    /* Must free storage pointed to by *tfr */
    LALFree ((void *) *tfr);
    ABORT (status, CREATETFR_EMALL, CREATETFR_MSGEMALL);
    return;
  }

  {
    INT4 column;
    (*tfr)->map = (REAL4 **) LALMalloc (in->tCol*sizeof(REAL4*));
    for (column = 0; column < (in->tCol); column++)
      (*tfr)->map[column] = (REAL4 *) LALMalloc ((in->fRow/2+1)*sizeof(REAL4));
  }
  
  if ((*tfr)->map == NULL)
  {
    /* Must free storage pointed to by *tfr */
    LALFree ((void *) *tfr);
    ABORT (status, CREATETFR_EMALL, CREATETFR_MSGEMALL);
    return;
  }
 
  /* Set TimeFreqRepType, fRow and tCol if storage allocated */

  (*tfr)->type = in->type;
  (*tfr)->fRow = in->fRow;	
  (*tfr)->tCol = in->tCol;

  /* We be done: Normal exit */

  RETURN (status);
}
