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

#include "LALRCSID.h"

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (LALERRORH, "$Id$");

int
LALPrintError (const char *fmt, ...);

void
LALAbort (const char *msg);


#ifdef  __cplusplus
}
#endif

#endif /* _LALERROR_H */
