/*----------------------------------------------------------------------- 
 * 
 * File Name: LALError.c
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "LALError.h"

NRCSID (LALERRORC, "$Id$");

int
LALPrintError (const char *fmt, ...)
{
  int n;
  va_list ap;
  va_start (ap, fmt);
  n = vfprintf (stderr, fmt, ap);
  va_end (ap);
  return n;
}


void
LALAbort (const char *msg)
{
  LALPrintError (msg);
  abort ();
}
