dnl m4 template for defining test program for different types. 
dnl To create the c-source for the typecode x test program, execute
dnl     m4 -DTYPECODE=x cvt.m4
dnl $Id$
dnl $Log$
dnl Revision 1.1  2000/02/16 01:58:29  jolien
dnl Initial entry
dnl
dnl
ifelse(TYPECODE,`Z',`define(`TYPE',`COMPLEX16')')
ifelse(TYPECODE,`C',`define(`TYPE',`COMPLEX8')')
ifelse(TYPECODE,`D',`define(`TYPE',`REAL8')')
ifelse(TYPECODE,`S',`define(`TYPE',`REAL4')')
ifelse(TYPECODE,`I2',`define(`TYPE',`INT2')')
ifelse(TYPECODE,`I4',`define(`TYPE',`INT4')')
ifelse(TYPECODE,`I8',`define(`TYPE',`INT8')')
ifelse(TYPECODE,`CHAR',`define(`TYPE',`CHAR')')
ifelse(TYPECODE,`',`define(`TYPE',`REAL4')')
define(`VTYPE',`format(`%sVector',TYPE)')
define(`PROGNAME',`format(`%sCreateVector',TYPECODE)')
define(`TESTPROG',`format(`%sTest.c',PROGNAME)')
/*----------------------------------------------------------------------- 
 * 
 * File Name: TESTPROG
 * 
 * Author: Finn, L. S. 
 * 
 * Revision: $Id$ 
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 *   main()
 *
 * SYNOPSIS 
 * 
 * DESCRIPTION 
 *   Test suite for vector factories
 * 
 * DIAGNOSTICS
 *   Writes PASS, FAIL to stdout as tests are passed or failed. 
 *   Returns 0 if all tests passed and -1 if any test failed. 
 * 
 * CALLS
 *   PROGNAME
 * 
 * NOTES
 *   PROGNAME can't test LALMalloc failure without messing with
 *   LALMalloc
 * 
 *   Automatically generated from cvt.tmpl by m4
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _STDLIB_H
#include <stdlib.h>
#ifndef _STDLIB_H
#define _STDLIB_H
#endif
#endif

#ifndef _STDIO_H
#include <stdio.h>
#ifndef _STDIO_H
#define _STDIO_H
#endif
#endif

#ifndef _STRING_H
#include <string.h>
#ifndef _STRING_H
#define _STRING_H
#endif
#endif

#ifndef _LALDATATYPES_H
#include "LALDatatypes.h"
#ifndef _LALDATATYPES_H
#define _LALDATATYPES_H
#endif
#endif

#ifndef _AVFACTORIES_H
#include "AVFactories.h"
#ifndef _AVFACTORIES_H
#define _AVFACTORIES_H
#endif
#endif

NRCSID (MAIN, "$Id$");

int debuglevel = 2;

typedef VTYPE vtype;
int main()
{

  static Status status;		/*  */
  vtype *vector = NULL;		/*  */
  vtype vstor;
  INT4 length = 0;
  int rcode = 0;		/* return code */

  /* Test error returns */

  /* Sequence length */

  rcode++;			/* increment code for each test */
  PROGNAME ( &status, &vector, 0 );
  if (status.statusCode != CREATEVECTOR_ELENGTH) 
    printf("FAIL zero sequence length: wrong status code\n");
  else if (strcmp(CREATEVECTOR_MSGELENGTH, status.statusDescription) != 0)
    printf("FAIL zero sequence length: statusDescription mismatch\n");
  else if (vector != NULL)
    printf("FAIL zero sequence length\n");
  else 
    { 
      rcode--;			/* decrement code if test passed */
      printf("PASS zero sequence length\n");
    }

  /* return structure */

  rcode++;			/* increment code for each test */
  length = 1;
  PROGNAME ( &status, NULL, length );
  if (status.statusCode != CREATEVECTOR_EVPTR) 
    printf("FAIL return struct ptr is NULL: wrong status code\n");
  else if (strcmp(CREATEVECTOR_MSGEVPTR, status.statusDescription) != 0)
    printf("FAIL return struct ptr is NULL: statusDescription mismatch\n");
  else 
    {
      rcode--;			/* decrement code if test passed */
      printf("PASS return struct ptr is NULL\n");
    }

  rcode++;			/* increment code for each test */
  vector = &vstor;
  PROGNAME ( &status, &vector, length );
  if (status.statusCode != CREATEVECTOR_EUPTR) 
    printf("FAIL return struct points to non-NULL: wrong status code\n");
  else if (strcmp(CREATEVECTOR_MSGEUPTR, status.statusDescription) != 0)
    printf("FAIL return struct points to non-NULL: statusDescription mismatch\n");
  else if (vector != &vstor)
    printf("FAIL return struct points to non-NULL: return struct reset\n");
  else 
    {
      rcode--;			/* decrement code if test passed */
      printf("PASS return struct points to non-NULL\n");
    }

  /* Test allocations */

  rcode++;			/* increment code for each test */
  length = 10;			/* allocate this length */
  vector = NULL;		/* make sure we aren't pointing to anything */
  PROGNAME ( &status, &vector, length );
  if (vector == NULL)
    printf("FAIL allocation: vector == NULL\n");
  else if (vector->data == NULL)
    printf("FAIL allocation: vector->data = NULL\n");
  else if (vector->length != length)
    printf("FAIL allocation: vector->length != length\n");
  else if (status.statusCode != 0)
    printf("FAIL allocation: non-zero status code legal allocation\n");
  else
    {
      rcode--;			/* decrement code if test passed */
      printf("PASS allocation\n");
    }

  /* We should have correctly set the status Id field: */

  printf ( "Test results for %s\n", status.Id );
  printf ( "Tested by %s\n", MAIN);

  exit((rcode==0)?0:-1);
}
  
