/*----------------------------------------------------------------------- 
 * 
 * File Name: LALHelloTest.c
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */


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

NRCSID (MAIN, "$Id$");

int debuglevel = 2;

int
main ()
{
  static Status status;
  LALHello (&status, NULL);
  REPORTSTATUS (&status);
  LALHello (&status, "protected");
  REPORTSTATUS (&status);
  return 0;
}
