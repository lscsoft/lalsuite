dnl m4 template for defining test program for different types. 
dnl To create the c-source for the typecode x test program, execute
dnl     m4 -DTYPECODE=x dvt.m4
dnl $Id$
ifelse(TYPECODE,`Z',`define(`TYPE',`COMPLEX16')')
ifelse(TYPECODE,`C',`define(`TYPE',`COMPLEX8')')
ifelse(TYPECODE,`D',`define(`TYPE',`REAL8')')
ifelse(TYPECODE,`S',`define(`TYPE',`REAL4')')
ifelse(TYPECODE,`I2',`define(`TYPE',`INT2')')
ifelse(TYPECODE,`I4',`define(`TYPE',`INT4')')
ifelse(TYPECODE,`I8',`define(`TYPE',`INT8')')
ifelse(TYPECODE,`U2',`define(`TYPE',`UINT2')')
ifelse(TYPECODE,`U4',`define(`TYPE',`UINT4')')
ifelse(TYPECODE,`U8',`define(`TYPE',`UINT8')')
ifelse(TYPECODE,`CHAR',`define(`TYPE',`CHAR')')
ifelse(TYPECODE,`',`define(`TYPE',`REAL4')')
define(`STYPE',`format(`%sVectorSequence',TYPE)')
define(`CPROG',`format(`%sCreateVectorSequence',TYPECODE)')
define(`DPROG',`format(`%sDestroyVectorSequence',TYPECODE)')
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
 *   Test suite for VectorSequence factories
 * 
 * DIAGNOSTICS
 * 
 * CALLS
 *   CPROG
 *   DPROG
 * 
 * NOTES
 *   Automatically generated from dvst.tmpl by m4
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

#ifndef _SEQFACTORIES_H
#include "SeqFactories.h"
#ifndef _SEQFACTORIES_H
#define _SEQFACTORIES_H
#endif
#endif

NRCSID (MAIN, "$Id$");

int debuglevel = 2;

typedef STYPE stype;
typedef TYPE dtype;

int main()
{

  static Status status;		/*  */
  stype *seq = NULL;		/*  */
  dtype *data = NULL;
  int rcode = 0;

  /* Announcements */
  printf ( "Testing %s\n" , "DPROG" );

  /* Test error returns */

  /* NULL seq */

  rcode++;	
  DPROG ( &status, NULL );
  if (status.statusCode != DESTROYVECSEQ_EVPTR)   
    printf("FAIL seq == NULL: wrong status code\n");
  else if (strcmp(status.statusDescription,DESTROYSEQ_MSGEVPTR) != 0)   
    printf("FAIL seq == NULL: statusDescription mismatch\n");
  else 
    {
      rcode--;
      printf("PASS seq == NULL\n");
    }

  /* NULL *seq */


  rcode++;
  seq = NULL;
  DPROG ( &status, &seq );
  if (status.statusCode != DESTROYVECSEQ_EUPTR)   
    printf("FAIL return struct ptr is NULL: wrong status code\n");
  else if (strcmp(status.statusDescription,DESTROYVECSEQ_MSGEUPTR) != 0)   
    printf("FAIL return struct ptr is NULL: statusDescription mismatch\n");
  else 
    {
      rcode--;
      printf("PASS return struct ptr is NULL\n");
    }

  /* seq->data == NULL */

  {
    CreateVectorSequenceIn in;
    in.length = 10;
    in.vectorLength = 20;
    (void) CPROG ( &status, &seq, &in );
  }

  rcode++;
  data = seq->data;		/* save pointer to data */
  seq->data = NULL;		/* make mal-formed seq */
  DPROG ( &status, &seq );
  if (status.statusCode != DESTROYVECSEQ_EDPTR)   
    printf("FAIL seq->data == NULL: wrong status code\n");
  else if (strcmp(status.statusDescription,DESTROYVECSEQ_MSGEDPTR) != 0)   
    printf("FAIL seq->data == NULL: statusDescription mismatch\n");
  else if (seq == NULL)
    printf("FAIL seq->data == NULL: seq changed\n");
  else 
    { 
      rcode--;
      printf("PASS seq->data == NULL\n");
    }

  /* Test destruction */  

  rcode++;
  seq->data = data;		/* we're correct again */
  DPROG ( &status, &seq );
  if ( seq != NULL )
    printf("FAIL Destruction: seq != NULL\n");
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
