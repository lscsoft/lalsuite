/*----------------------------------------------------------------------- 
 * 
 * File Name: tmpltbank.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <config.h>
#ifndef HAVE_GETOPT_H
#include <stdio.h>
int main( void )
{
  fputs( "Disabled: LALApps compiled without getopt.h\n", stderr );
  return 77;
}
#else

#include "ligolwbank.h"
#define BANK_FILE "tmpltbank.xml"

int main ( int argc, char *argv[] )
{
  exit( 0 );
}

#endif

