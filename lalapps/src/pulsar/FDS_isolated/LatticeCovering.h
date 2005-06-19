/*
 * Copyright (C) 2005 Reinhard Prix
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the 
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
 *  MA  02111-1307  USA
 */

/**
 * \author Reinhard Prix
 * \date 2005
 * \file 
 * \ingroup moduleLatticeCovering
 * \brief Header-file defining the API for the lattice-covering functions.
 *
 * $Id$
 *
 */

#ifndef _LATTICECOVERING_H  /* Double-include protection. */
#define _LATTICECOVERING_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/*---------- INCLUDES ----------*/
#include <gsl/gsl_matrix.h>

#include <lal/LALStdlib.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>



NRCSID( DOPPLERSCANH, "$Id$" );


/*---------- DEFINES ----------*/

/*----- Error-codes -----<lalLaTeX>
\subsection*{Error codes}</lalLaTeX><lalErrTable> */

#define LATTICECOVERING_ENULL 		1
#define LATTICECOVERING_ENONULL		2
#define LATTICECOVERING_EMEM		3
#define LATTICECOVERING_EINPUT		4
#define LATTICECOVERING_ELIST		5
#define LATTICECOVERING_EFUNC		6

#define LATTICECOVERING_MSGENULL 	"Arguments contained an unexpected null pointer"
#define LATTICECOVERING_MSGENONULL	"Output pointer is not NULL"
#define LATTICECOVERING_MSGEMEM		"Out of memory"
#define LATTICECOVERING_MSGEINPUT	"Invald input parameter"
#define LATTICECOVERING_MSGELIST	"Error occurred in list-handling ..."
#define LATTICECOVERING_MSGEFUNC	"Sub-routine failed"

/*</lalErrTable> */

/*---------- external types ----------*/

/** enum-type for denoting several types of lattice */
typedef enum
{
  LATTICE_TYPE_CUBIC = 0,	/**< standard cubic grid: Zn */
  LATTICE_TYPE_ANSTAR,		/**< An*: optimal covering grid */
  LATTICE_TYPE_LAST
} LatticeType;

/** doubly linked list of INT4-vectors (lattice-vectors) */
typedef struct tagINT4VectorList
{
  INT4Vector entry;
  struct tagINT4VectorList *next;
  struct tagINT4VectorList *prev;
} INT4VectorList;

/** doubly linked list of REAL8-vectors (physical vectors) */
typedef struct tagREAL8VectorList
{
  REAL8Vector entry;
  struct tagREAL8VectorList *next;
  struct tagREAL8VectorList *prev;
} REAL8VectorList;

/*---------- Global variables ----------*/

/*---------- external prototypes [API] ----------*/
void LALLatticeCovering (LALStatus *, REAL8VectorList **covering, REAL8 coveringRadius, 
			 const gsl_matrix *metric, const REAL8Vector *startPoint,
			 BOOLEAN (*isInside)(const REAL8Vector *point) );

void LALLatticeFill (LALStatus *, REAL8VectorList **fillGrid, const gsl_matrix  *generator,
		     const REAL8Vector *startPoint, BOOLEAN (*isInside)(const REAL8Vector *point) );

/* functions for handling lattice's generating matrix */
int XLALFindCoveringGenerator (gsl_matrix **outmatrix, LatticeType type, UINT4 dimension, 
			       REAL8 coveringRadius, const gsl_matrix *gij);
int XLALReduceGenerator2FullRank(gsl_matrix **outmatrix, const gsl_matrix *inmatrix);
int XLALGetLatticeGenerator (gsl_matrix **outmatrix, UINT4 dimension, LatticeType type);

/* some REAL8 list-handling functions that might be useful to users of the above functions */
REAL8VectorList* REAL8VectorListAddEntry (REAL8VectorList *head, const REAL8Vector *entry);
void REAL8VectorListDestroy (REAL8VectorList *head);



#ifdef  __cplusplus
}
#endif  
/* C++ protection. */

#endif  /* Double-include protection. */
