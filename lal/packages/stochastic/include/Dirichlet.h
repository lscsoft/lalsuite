/*-----------------------------------------------------------------------
 *
 * File Name: Dirichlet.h
 *
 * Author: UTB Relativity Group
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * Dirichlet.h
 *
 * SYNOPSIS
 * #include <lal/Dirichlet.h>
 *
 * DESCRIPTION
 * Error codes, typedefs, and prototype for function LALDirichlet()
 *
 * DIAGNOSTICS
 *
 *-----------------------------------------------------------------------
 */

#ifndef _DIRICHLET_H
#define _DIRICHLET_H

#include <lal/LALStdlib.h>

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (DIRICHLETH, "$Id$");

#define DIRICHLET_ENULLIP     1
#define DIRICHLET_ENVALUE     2
#define DIRICHLET_ESIZE       4
#define DIRICHLET_EDELTAX     8
#define DIRICHLET_ENULLOP     16
#define DIRICHLET_ESIZEMM     32
#define DIRICHLET_ENULLD      64

#define DIRICHLET_MSGENULLIP  "Pointer to input parameters must be non-null"
#define DIRICHLET_MSGENVALUE  "LALDirichlet parameter N must be > 0"
#define DIRICHLET_MSGESIZE    "Specified length of output vector must be > 0"
#define DIRICHLET_MSGEDELTAX  "Spacing of x values must be > 0"
#define DIRICHLET_MSGENULLOP  "Pointer to ouput vector must be non-null"
#define DIRICHLET_MSGESIZEMM  "Length of output vector must agree with length specified in input parameters"
#define DIRICHLET_MSGENULLD   "Pointer to data member of output vector must be non-null"

typedef struct tagDirichletParameters{
  INT4	 n;       /* LALDirichlet parameter N */
  INT4	 length;  /* specified length of output vector */
  REAL8	 deltaX;  /* spacing of x values */
} DirichletParameters; 

void 
LALDirichlet ( LALStatus*, REAL4Vector*, DirichletParameters* );


#ifdef  __cplusplus
}
#endif

#endif /* _DIRICHLET_H */
