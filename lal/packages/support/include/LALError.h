/*----------------------------------------------------------------------- 
 * 
 * File Name: LALError.h
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _LALERROR_H
#define _LALERROR_H

#ifndef _LALRCSID_H
#include "LALRCSID.h"
#ifndef _LALRCSID_H
#define _LALRCSID_H
#endif
#endif

NRCSID (LALERRORH, "$Id$");

int
LALPrintError (const char *fmt, ...);

void
LALAbort (const char *msg);

#endif
