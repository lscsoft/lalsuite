/*----------------------------------------------------------------------- 
 * 
 * File Name: LALError.c
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

#ifndef _STDARG_H
#include <stdarg.h>
#ifndef _STDARG_H
#define _STDARG_H
#endif
#endif

#ifndef _LALERROR_H
#include "LALError.h"
#ifndef _LALERROR_H
#define _LALERROR_H
#endif
#endif

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
