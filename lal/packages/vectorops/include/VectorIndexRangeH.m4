/* -*- C -*- */
/******** <lalVerbatim file="VectorIndexRangeHV"> ********
Author: David Chin <dwchin@umich.edu> +1-734-709-9119
$Id$
********* </lalVerbatim> ********/

/**** <lalLaTeX>
 *
 * \section{Header \texttt{VectorIndexRange.h}}
 *
 * Routines to slice up vectors/sequences by
 * specifying starting and ending indices, inclusive.
 * \texttt{VectorIndexHole()}
 * is "complementary" to \texttt{VectorIndexRange()}, in the sense that
 * it returns the two segments of the vector that would remain after
 * \texttt{VectorIndexRange()} acts.
 *
 * \subsection*{Synopsis}
 * \begin{verbatim}
 * #include <lal/VectorIndexRange.h>
 * \end{verbatim}
 *
 * \subsection*{Error conditions}
 * \input{VectorIndexRangeHE}
 *
 * \vfill{\footnotesize\input{VectorIndexRangeHV}}
 * \newpage\input{VectorIndexRangeC}
 *
 ***** </lalLaTeX> */

#include <math.h>
#include <lal/LALRCSID.h>
#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>
#include <lal/AVFactories.h>

#ifndef __LALVECTORINDEXRANGE_H__
#define __LALVECTORINDEXRANGE_H__

#ifdef __cplusplus
extern "C" {
#endif /* C++ */

NRCSID( VECTORINDEXRANGEH, "$Id$" );

/**** <lalErrTable file="VectorIndexRangeHE"> ****/
#define VECTORINDEXRANGEH_EARG 1
#define VECTORINDEXRANGEH_ECHK 2
#define VECTORINDEXRANGEH_EFLS 3
#define VECTORINDEXRANGEH_EUSE 4
#define VECTORINDEXRANGEH_ENULL 5
#define VECTORINDEXRANGEH_EALOC 6
#define VECTORINDEXRANGEH_EFPMS 7
#define VECTORINDEXRANGEH_ENUMZ 8
#define VECTORINDEXRANGEH_ELNTH 9
#define VECTORINDEXRANGEH_ENNUL 10
#define VECTORINDEXRANGEH_MSGEARG "Error parsing command-line arguments"
#define VECTORINDEXRANGEH_MSGECHK "Error checking failed to catch bad data"
#define VECTORINDEXRANGEH_MSGEFLS "Incorrect answer for valid data"
#define VECTORINDEXRANGEH_MSGEUSE "Bad user-entered data"
#define VECTORINDEXRANGEH_MSGENULL "Null Pointer."
#define VECTORINDEXRANGEH_MSGEALOC "Memory Allocation Error"
#define VECTORINDEXRANGEH_MSGEFPMS "Filter Parameter Structure Error"
#define VECTORINDEXRANGEH_MSGENUMZ "Incorrect number of command line arguments"
#define VECTORINDEXRANGEH_MSGELNTH "Vector/Array of Improper Length"
#define VECTORINDEXRANGEH_MSGENNUL "Non-Null Pointer that should be NULL"
/**** </lalErrTable> ****/

/* typedefs */

define(`TYPECODE',`CHAR')
include(`VectorPairTypedefsBase.m4')

define(`TYPECODE',`I2')
include(`VectorPairTypedefsBase.m4')

define(`TYPECODE',`I4')
include(`VectorPairTypedefsBase.m4')

define(`TYPECODE',`I8')
include(`VectorPairTypedefsBase.m4')

define(`TYPECODE',`U2')
include(`VectorPairTypedefsBase.m4')

define(`TYPECODE',`U4')
include(`VectorPairTypedefsBase.m4')

define(`TYPECODE',`U8')
include(`VectorPairTypedefsBase.m4')

define(`TYPECODE',`S')
include(`VectorPairTypedefsBase.m4')

define(`TYPECODE',`D')
include(`VectorPairTypedefsBase.m4')

define(`TYPECODE',`C')
include(`VectorPairTypedefsBase.m4')

define(`TYPECODE',`Z')
include(`VectorPairTypedefsBase.m4')

/* Function protos */

define(`TYPECODE',`CHAR')
include(`VectorIndexRangeBaseH.m4')

define(`TYPECODE',`I2')
include(`VectorIndexRangeBaseH.m4')

define(`TYPECODE',`I4')
include(`VectorIndexRangeBaseH.m4')

define(`TYPECODE',`I8')
include(`VectorIndexRangeBaseH.m4')

define(`TYPECODE',`U2')
include(`VectorIndexRangeBaseH.m4')

define(`TYPECODE',`U4')
include(`VectorIndexRangeBaseH.m4')

define(`TYPECODE',`U8')
include(`VectorIndexRangeBaseH.m4')

define(`TYPECODE',`S')
include(`VectorIndexRangeBaseH.m4')

define(`TYPECODE',`D')
include(`VectorIndexRangeBaseH.m4')

define(`TYPECODE',`C')
include(`VectorIndexRangeBaseH.m4')

define(`TYPECODE',`Z')
include(`VectorIndexRangeBaseH.m4')

#ifdef __cplusplus
}
#endif /* C++ */

#endif /* __LALVECTORINDEXRANGE_H__ */


