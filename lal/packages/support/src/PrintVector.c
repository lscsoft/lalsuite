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

#include "LALStdlib.h"
#include <stdio.h>
#include "LALDatatypes.h"
#include "PrintVector.h"

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
