dnl m4 template for defining test program for different types. 
dnl To create the c-source for the typecode x test program, execute
dnl     m4 -DTYPECODE=x cvt.m4
dnl $Id$
dnl $Log$
dnl Revision 1.2  2000/02/26 21:18:11  jolien
dnl Modified to accomodate new unsigned types
dnl
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
ifelse(TYPECODE,`U2',`define(`TYPE',`UINT2')')
ifelse(TYPECODE,`U4',`define(`TYPE',`UINT4')')
ifelse(TYPECODE,`U8',`define(`TYPE',`UINT8')')
ifelse(TYPECODE,`CHAR',`define(`TYPE',`CHAR')')
ifelse(TYPECODE,`',`define(`TYPE',`REAL4')')
define(`STYPE',`format(`%sVectorSequence',TYPE)')
define(`PROGNAME',`format(`%sCreateVectorSequence',TYPECODE)')
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

#ifndef _SEQFACTORIES_H
#include "SeqFactories.h"
#ifndef _SEQFACTORIES_H
#define _SEQFACTORIES_H
#endif
#endif

NRCSID (MAIN, "$Id$");

int debuglevel = 2;

typedef STYPE stype;
int main()
{

  static Status status;		/*  */
  stype *seq = NULL;		/*  */
  stype sstor;
  CreateVectorSequenceIn in;
  int rcode = 0;		/* return code */

  printf ( "Testing %s\n" , "PROGNAME" );

  /* Test error returns */

  /* Sequence length */

  rcode++;			/* increment code for each test */

  in.length = 0;
  in.vectorLength = 1;
  PROGNAME ( &status, &seq, &in );
  if (status.statusCode != CREATEVECSEQ_ESLENGTH) 
    printf("FAIL zero sequence length: wrong status code\n");
  else if (strcmp(CREATEVECSEQ_MSGESLENGTH, status.statusDescription) != 0)
    printf("FAIL zero sequence length: statusDescription mismatch\n");
  else if (seq != NULL)
    printf("FAIL zero sequence length: seq != NULL\n");
  else 
    { 
      rcode--;			/* decrement code if test passed */
      printf("PASS zero sequence length\n");
    }

  /* vector_sequence length */

  rcode++;			/* increment code for each test */

  in.length = 1;
  in.vectorLength = 0;
  PROGNAME ( &status, &seq, &in );
  if (status.statusCode != CREATEVECSEQ_EVLENGTH) 
    printf("FAIL zero vectorLength: wrong status code\n");
  else if (strcmp(CREATEVECSEQ_MSGEVLENGTH, status.statusDescription) != 0)
    printf("FAIL zero vectorLength: statusDescription mismatch\n");
  else if (seq != NULL)
    printf("FAIL zero vectorLength: seq != NULL\n");
  else 
    { 
      rcode--;			/* decrement code if test passed */
      printf("PASS zero vectorLength\n");
    }

  /* input structure */

  rcode++;			/* increment code for each test */
  PROGNAME ( &status, &seq, NULL );
  if (status.statusCode != CREATEVECSEQ_EINPTR) 
    printf("FAIL input struct ptr is NULL: wrong status code\n");
  else if (strcmp(CREATEVECSEQ_MSGEINPTR, status.statusDescription) != 0)
    printf("FAIL input struct ptr is NULL: statusDescription mismatch\n");
  else 
    {
      rcode--;			/* decrement code if test passed */
      printf("PASS input struct ptr is NULL\n");
    }

  /* return structure */

  rcode++;			/* increment code for each test */
  in.length = 10;
  in.vectorLength = 2;
  PROGNAME ( &status, NULL, &in );
  if (status.statusCode != CREATEVECSEQ_EVPTR) 
    printf("FAIL return struct ptr is NULL: wrong status code\n");
  else if (strcmp(CREATEVECSEQ_MSGEVPTR, status.statusDescription) != 0)
    printf("FAIL return struct ptr is NULL: statusDescription mismatch\n");
  else 
    {
      rcode--;			/* decrement code if test passed */
      printf("PASS return struct ptr is NULL\n");
    }

  rcode++;			/* increment code for each test */
  seq = &sstor;
  PROGNAME ( &status, &seq, &in );
  if (status.statusCode != CREATEVECSEQ_EUPTR) 
    printf("FAIL return struct points to non-NULL: wrong status code\n");
  else if (strcmp(CREATEVECSEQ_MSGEUPTR, status.statusDescription) != 0)
    printf("FAIL return struct points to non-NULL: statusDescription mismatch\n");
  else if (seq != &sstor)
    printf("FAIL return struct points to non-NULL: return struct reset\n");
  else 
    {
      rcode--;			/* decrement code if test passed */
      printf("PASS return struct points to non-NULL\n");
    }

  /* Test allocations */

  rcode++;			/* increment code for each test */
  in.length = 10;		/* allocate this length */
  in.vectorLength = 5;		/* this size vector */
  seq = NULL;		/* make sure we aren't pointing to anything */
  PROGNAME ( &status, &seq, &in );
  if (status.statusCode != status.statusCode)   
    printf("FAIL allocation: status code mismatch\n");
  else if (seq == NULL)
    printf("FAIL allocation: seq == NULL\n");
  else if (seq->data == NULL)
    printf("FAIL allocation: seq->data = NULL\n");
  else if (seq->length != in.length)
    printf("FAIL allocation: seq->length != in->length\n");
  else if (seq->vectorLength != in.vectorLength)
    printf("FAIL allocation: seq->vectorLength != in->vectorLength\n");
  else if (status.statusCode != 0)
    printf("FAIL allocation: non-zero status code on legal allocation\n");
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
  
