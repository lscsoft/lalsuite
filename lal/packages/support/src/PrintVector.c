/*----------------------------------------------------------------------- 
 * 
 * File Name: PrintVector.c
 * 
 * Author: Allen, Bruce ballen@dirac.phys.uwm.edu
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * PrintVector
 * 
 * SYNOPSIS 
 * void PrintVector (REAL4Vector *vector);
 * 
 * DESCRIPTION 
 * Print a REAL4Vector object. 
 * 
 * DIAGNOSTICS 
 * vector != NULL
 *
 * CALLS
 * 
 * NOTES
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _LALSTDLIB_H
#include "LALStdlib.h"
#ifndef _LALSTDLIB_H
#define _LALSTDLIB_H
#endif
#endif

#ifndef _STDIO_H
#include "stdio.h"
#ifndef _STDIO_H
#define _STDIO_H
#endif
#endif

#ifndef _LALDATATYPES_H
#include "LALDatatypes.h"
#ifndef _LALDATATYPES_H
#define _LALDATATYPES_H
#endif
#endif

#ifndef _PRINTVECTOR_H
#include "PrintVector.h"
#ifndef _PRINTVECTOR_H
#define _PRINTVECTOR_H
#endif
#endif

NRCSID (PRINTVECTORC, "$Id$");

void
PrintVector(REAL4Vector *vector)
{
  int i;
  static int fileno=0;
  FILE *fp;
  char fname[256];


  if (vector==NULL) return;

  /* open output file */
  sprintf(fname,"PrintVector.%03d",fileno++);
  fp=fopen(fname,"w");

  for (i=0;i<vector->length;i++)
    fprintf(fp,"%i %f\n",i,vector->data[i]);

  fclose(fp);

  return;
}
