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
 * \defgroup moduleLatticeCovering Lattice Covering
 * \ingroup coreDA
 *
 * \brief Module implementing a practical (approximate) solution to the "covering problem".
 *  This is achieved by producing a lattice-covering with the \f$A_n^*\f$ lattice,
 *  which is the best known covering up to dimension \f$n\le 23\f$
 *  (see \ref CS99 for details).
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

/* we need these until we put them in some better LAL-place (FIXME) */
#include <lal/PtoleMetric.h>



NRCSID( LATTICECOVERINGH, "$Id$" );


/*---------- DEFINES ----------*/

/*----- Error-codes -----*/

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

/*---------- exported types ----------*/

/** enum-type for denoting several types of lattice */
typedef enum
{
  LATTICE_TYPE_ANSTAR = 0,	/**< An*: optimal covering grid */
  LATTICE_TYPE_CUBIC,		/**< standard cubic grid: Zn */
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
/*---------- empty initializers ---------- */
extern INT4VectorList empty_INT4VectorList;
extern REAL8VectorList empty_REAL8VectorList;

/*---------- exported prototypes [API] ----------*/
void LALLatticeCovering (LALStatus *, REAL8VectorList **covering, REAL8 coveringRadius,
			 const gsl_matrix *metric, const REAL8Vector *startPoint,
			 BOOLEAN (*isInside)(const REAL8Vector *point),
			 LatticeType latticeType);

void LALLatticeFill (LALStatus *, REAL8VectorList **fillGrid, const gsl_matrix  *generator,
		     const REAL8Vector *startPoint, BOOLEAN (*isInside)(const REAL8Vector *point) );

/* functions for handling lattice's generating matrix */
int XLALFindCoveringGenerator (gsl_matrix **outmatrix, LatticeType type, REAL8 coveringRadius, const gsl_matrix *gij);
int XLALReduceGenerator2FullRank(gsl_matrix **outmatrix, const gsl_matrix *inmatrix);
int XLALGetLatticeGenerator (gsl_matrix **outmatrix, UINT4 dimension, LatticeType type);

REAL8 XLALMetricScalarProduct (const gsl_vector *vector1, const gsl_vector *vector2,
			       const gsl_matrix *metric);


/* some REAL8 list-handling functions that might be useful to users of the above functions */
REAL8VectorList* XLALREAL8VectorListAddEntry (REAL8VectorList *head, const REAL8Vector *entry);
void XLALREAL8VectorListDestroy (REAL8VectorList *head);

/* convert a symmetric gsl_matrix into a 'LAL-encoded' vector containing only
 * the independent symmetric-matrix elements (see PMETRIC_INDEX) */
REAL8Vector *XLALgsl2LALmetric (const gsl_matrix *gmetric);

/* convert LAL-encoded metric into a symmetric gsl-matrix */
gsl_matrix *XLALmetric2gsl (const REAL8Vector *metric);

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
