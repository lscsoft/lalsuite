dnl m4 template for defining test program for different types. 
dnl To create the c-source for the typecode x test program, execute
dnl     m4 -DTYPECODE=x dvt.m4
dnl $Id$
ifelse(TYPECODE,`Z',`define(`TYPE',`COMPLEX16')')
ifelse(TYPECODE,`C',`define(`TYPE',`COMPLEX8')')
ifelse(TYPECODE,`D',`define(`TYPE',`REAL8')')
ifelse(TYPECODE,`S',`define(`TYPE',`REAL4')')
ifelse(TYPECODE,`I4',`define(`TYPE',`INT4')')
ifelse(TYPECODE,`I8',`define(`TYPE',`INT8')')
ifelse(TYPECODE,`I2',`define(`TYPE',`INT2')')
ifelse(TYPECODE,`U4',`define(`TYPE',`UINT4')')
ifelse(TYPECODE,`U8',`define(`TYPE',`UINT8')')
ifelse(TYPECODE,`U2',`define(`TYPE',`UINT2')')
ifelse(TYPECODE,`CHAR',`define(`TYPE',`CHAR')')
ifelse(TYPECODE,`',`define(`TYPE',`REAL4')')
define(`VTYPE',`format(`%sVector',TYPE)')
define(`CPROG',`format(`%sCreateVector',TYPECODE)')
define(`DPROG',`format(`%sDestroyVector',TYPECODE)')
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
 * 
 * CALLS
 *   CPROG
 *   DPROG
 * 
 * NOTES
 *   Automatically generated from dvt.tmpl by m4
 *
 * 
 *-----------------------------------------------------------------------
 */

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
typedef TYPE dtype;

int main()
{

  static Status status;		/*  */
  vtype *vector = NULL;		/*  */
  dtype *data = NULL;
  int rcode = 0;

  /* Announcements */
  printf ( "Testing %s\n" , "DPROG" );

  /* Test error returns */

  /* NULL vector */

  rcode++;	
  DPROG ( &status, NULL );
  if (status.statusCode != DESTROYVECTOR_EVPTR)   
    printf("FAIL vector == NULL: wrong status code\n");
  else if (strcmp(status.statusDescription,DESTROYVECTOR_MSGEVPTR) != 0)   
    printf("FAIL vector == NULL: statusDescription mismatch\n");
  else 
    {
      rcode--;
      printf("PASS vector == NULL\n");
    }

  /* NULL *vector */

  rcode++;
  DPROG ( &status, &vector );
  if (status.statusCode != DESTROYVECTOR_EUPTR)   
    printf("FAIL return struct ptr is NULL: wrong status code\n");
  else if (strcmp(status.statusDescription,DESTROYVECTOR_MSGEUPTR) != 0)   
    printf("FAIL return struct ptr is NULL: statusDescription mismatch\n");
  else 
    {
      rcode--;
      printf("PASS return struct ptr is NULL\n");
    }

  /* vector->data == NULL */

  {
    UINT4 length = 10;
    CPROG ( &status, &vector, length );
  }

  rcode++;
  data = vector->data;		/* save pointer to data */
  vector->data = NULL;		/* make mal-formed vector */
  DPROG ( &status, &vector );
  if (status.statusCode != DESTROYVECTOR_EDPTR)   
    printf("FAIL vector->data == NULL: wrong status code\n");
  else if (strcmp(status.statusDescription,DESTROYVECTOR_MSGEDPTR) != 0)   
    printf("FAIL vector->data == NULL: statusDescription mismatch\n");
  else if (vector == NULL)
    printf("FAIL vector->data == NULL: vector changed\n");
  else 
    { 
      rcode--;
      printf("PASS vector->data == NULL\n");
    }

  /* Test destruction */  

  rcode++;
  vector->data = data;		/* we're correct again */
  DPROG ( &status, &vector );
  if ( vector != NULL )
    printf("FAIL Destruction: vector != NULL\n");
  else if (status.statusCode != 0)
    printf("FAIL Destruction: normal exit with non-zero statusCode");
  else
    {
      rcode--;
      printf("PASS Destruction\n");
    }

  /* We should have correctly set the status Id field: */

  printf ( "Test results for %s\n", status.Id );
  printf ( "Tested by %s\n", MAIN);

  exit(rcode==0?0:-1);
}
