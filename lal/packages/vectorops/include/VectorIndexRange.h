/* -*- C -*- */
/******** <lalVerbatim file="VectorIndexRangeCV"> ********
Author: David Chin <dwchin@umich.edu> +1-734-709-9119
$Id$
********* </lalVerbatim> ********/

/**** <lalLaTeX>
 *
 * \section{Header \textttVectorIndexRange.h}}
 *
 * Routines to slice up vectors/sequences by
 * specifying starting and ending indices.
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


/*













*/


void LALCHARVectorIndexRange ( LALStatus *status,
          CHARVector    **result,
          CHARVector     *v,
          const UINT4Vector *indexVector );
         



/*













*/


void LALI2VectorIndexRange ( LALStatus *status,
          INT2Vector    **result,
          INT2Vector     *v,
          const UINT4Vector *indexVector );
         



/*













*/


void LALI4VectorIndexRange ( LALStatus *status,
          INT4Vector    **result,
          INT4Vector     *v,
          const UINT4Vector *indexVector );
         



/*













*/


void LALI8VectorIndexRange ( LALStatus *status,
          INT8Vector    **result,
          INT8Vector     *v,
          const UINT4Vector *indexVector );
         



/*













*/


void LALU2VectorIndexRange ( LALStatus *status,
          UINT2Vector    **result,
          UINT2Vector     *v,
          const UINT4Vector *indexVector );
         



/*













*/


void LALU4VectorIndexRange ( LALStatus *status,
          UINT4Vector    **result,
          UINT4Vector     *v,
          const UINT4Vector *indexVector );
         



/*













*/


void LALU8VectorIndexRange ( LALStatus *status,
          UINT8Vector    **result,
          UINT8Vector     *v,
          const UINT4Vector *indexVector );
         



/*













*/


void LALSVectorIndexRange ( LALStatus *status,
          REAL4Vector    **result,
          REAL4Vector     *v,
          const UINT4Vector *indexVector );
         



/*













*/


void LALDVectorIndexRange ( LALStatus *status,
          REAL8Vector    **result,
          REAL8Vector     *v,
          const UINT4Vector *indexVector );
         



/*













*/


void LALCVectorIndexRange ( LALStatus *status,
          COMPLEX8Vector    **result,
          COMPLEX8Vector     *v,
          const UINT4Vector *indexVector );
         



/*













*/


void LALZVectorIndexRange ( LALStatus *status,
          COMPLEX16Vector    **result,
          COMPLEX16Vector     *v,
          const UINT4Vector *indexVector );
         


#ifdef __cplusplus
}
#endif /* C++ */

#endif /* __LALVECTORINDEXRANGE_H__ */


