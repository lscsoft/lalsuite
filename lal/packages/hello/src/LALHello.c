/*----------------------------------------------------------------------- 
 * 
 * File Name: LALHello.c
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _STDIO_H
#include <stdio.h>
#ifndef _STDIO_H
#define _STDIO_H
#endif
#endif

#ifndef _LALSTDLIB_H
#include "LALStdlib.h"
#ifndef _LALSTDLIB_H
#define _LALSTDLIB_H
#endif
#endif

#ifndef _LALHELLO_H
#include "LALHello.h"
#ifndef _LALHELLO_H
#define _LALHELLO_H
#endif
#endif

NRCSID (LALHELLOC, "$Id$");


/* don't want this function to be visible outside this file */
static void
LALPrintMessage (Status *status, const CHAR *message, const CHAR *fileName)
{
  INT4  error;
  INT4  numChar;
  FILE *filePtr;

  INITSTATUS (status, LALHELLOC);

  if (fileName)
  {
    filePtr = LALFopen (fileName, "w");
  }
  else
  {
    filePtr = stdout;
  }

  ASSERT (filePtr, status, LALHELLO_EOPEN, LALHELLO_MSGEOPEN);
    
  numChar = fprintf (filePtr, message);

  ASSERT (numChar >= 0, status, LALHELLO_EWRITE, LALHELLO_MSGEWRITE);

  if (fileName)
  {
    error = LALFclose (filePtr);
    ASSERT (error == 0, status, LALHELLO_ECLOSE, LALHELLO_MSGECLOSE);
  }
  else
  {
    error = fflush (stdout);
    ASSERT (error == 0, status, LALHELLO_EFLUSH, LALHELLO_MSGEFLUSH);
  }

  RETURN (status);
}


void
LALHello (Status *status, const CHAR *fileName)
{
  INITSTATUS (status, LALHELLOC);
  ATTATCHSTATUSPTR (status);

  LALPrintMessage (status->statusPtr, "hello, LSC!\n", fileName);
  CHECKSTATUSPTR (status);
  
  DETATCHSTATUSPTR (status);
  RETURN (status);  
}
