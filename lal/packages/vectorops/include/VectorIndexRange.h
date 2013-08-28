/* -*- C -*- */
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>
#include <lal/AVFactories.h>

#ifndef __LALVECTORINDEXRANGE_H__
#define __LALVECTORINDEXRANGE_H__

#ifdef __cplusplus
extern "C" {
#endif /* C++ */

/**
 * \addtogroup VectorIndexRange_h
 * \author David Chin <dwchin@umich.edu> +1-734-709-9119
 *
 * \brief Routines to slice up vectors/sequences by specifying starting and ending indices, inclusive.
 *
 * VectorIndexHole() is "complementary" to VectorIndexRange(), in the sense that
 * it returns the two segments of the vector that would remain after VectorIndexRange() acts.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/VectorIndexRange.h>
 * \endcode
 */
/*@{*/

/**\name Error Codes */
/*@{*/
#define VECTORINDEXRANGEH_EARG 1	/**< Error parsing command-line arguments */
#define VECTORINDEXRANGEH_ECHK 2	/**< Error checking failed to catch bad data */
#define VECTORINDEXRANGEH_EFLS 3	/**< Incorrect answer for valid data */
#define VECTORINDEXRANGEH_EUSE 4	/**< Bad user-entered data */
#define VECTORINDEXRANGEH_ENULL 5	/**< Null Pointer */
#define VECTORINDEXRANGEH_EALOC 6	/**< Memory Allocation Error */
#define VECTORINDEXRANGEH_EFPMS 7	/**< Filter Parameter Structure Error */
#define VECTORINDEXRANGEH_ENUMZ 8	/**< Incorrect number of command line arguments */
#define VECTORINDEXRANGEH_ELNTH 9	/**< Vector/Array of Improper Length */
#define VECTORINDEXRANGEH_ENNUL 10	/**< Non-Null Pointer that should be NULL */
/*@}*/

/** \cond DONT_DOXYGEN */
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
/** \endcond */


/* ---------- typedefs ---------- */

/* CHAR */
typedef struct
tagCHARVectorPair {
  CHARVector **head;
  CHARVector **tail;
} CHARVectorPair;

/* INT2 */
typedef struct
tagINT2VectorPair {
  INT2Vector **head;
  INT2Vector **tail;
} INT2VectorPair;

/* INT4 */
typedef struct
tagINT4VectorPair {
  INT4Vector **head;
  INT4Vector **tail;
} INT4VectorPair;

/* INT8 */
typedef struct
tagINT8VectorPair {
  INT8Vector **head;
  INT8Vector **tail;
} INT8VectorPair;

/* UINT2 */
typedef struct
tagUINT2VectorPair {
  UINT2Vector **head;
  UINT2Vector **tail;
} UINT2VectorPair;

/* UINT4 */
typedef struct
tagUINT4VectorPair {
  UINT4Vector **head;
  UINT4Vector **tail;
} UINT4VectorPair;

/* UINT8 */
typedef struct
tagUINT8VectorPair {
  UINT8Vector **head;
  UINT8Vector **tail;
} UINT8VectorPair;

/* REAL4 */
typedef struct
tagREAL4VectorPair {
  REAL4Vector **head;
  REAL4Vector **tail;
} REAL4VectorPair;

/* REAL8 */
typedef struct
tagREAL8VectorPair {
  REAL8Vector **head;
  REAL8Vector **tail;
} REAL8VectorPair;

/* COMPLEX8 */
typedef struct
tagCOMPLEX8VectorPair {
  COMPLEX8Vector **head;
  COMPLEX8Vector **tail;
} COMPLEX8VectorPair;

/* COMPLEX16 */
typedef struct
tagCOMPLEX16VectorPair {
  COMPLEX16Vector **head;
  COMPLEX16Vector **tail;
} COMPLEX16VectorPair;


/* ---------- Function protos ---------- */

/* CHAR */
void LALCHARVectorIndexRange ( LALStatus *status,
          CHARVector    **result,
          CHARVector     *v,
          const UINT4Vector *indexVector );
void LALCHARVectorIndexHole ( LALStatus  *status,
          CHARVectorPair  *result_pair,
          CHARVector      *v,
          const UINT4Vector *indexVector );

/* INT2 */
void LALI2VectorIndexRange ( LALStatus *status,
          INT2Vector    **result,
          INT2Vector     *v,
          const UINT4Vector *indexVector );
void LALI2VectorIndexHole ( LALStatus  *status,
          INT2VectorPair  *result_pair,
          INT2Vector      *v,
          const UINT4Vector *indexVector );

/* INT4 */
void LALI4VectorIndexRange ( LALStatus *status,
          INT4Vector    **result,
          INT4Vector     *v,
          const UINT4Vector *indexVector );
void LALI4VectorIndexHole ( LALStatus  *status,
          INT4VectorPair  *result_pair,
          INT4Vector      *v,
          const UINT4Vector *indexVector );

/* INT8 */
void LALI8VectorIndexRange ( LALStatus *status,
          INT8Vector    **result,
          INT8Vector     *v,
          const UINT4Vector *indexVector );
void LALI8VectorIndexHole ( LALStatus  *status,
          INT8VectorPair  *result_pair,
          INT8Vector      *v,
          const UINT4Vector *indexVector );

/* UINT2 */
void LALU2VectorIndexRange ( LALStatus *status,
          UINT2Vector    **result,
          UINT2Vector     *v,
          const UINT4Vector *indexVector );
void LALU2VectorIndexHole ( LALStatus  *status,
          UINT2VectorPair  *result_pair,
          UINT2Vector      *v,
          const UINT4Vector *indexVector );

/* UINT4 */
void LALU4VectorIndexRange ( LALStatus *status,
          UINT4Vector    **result,
          UINT4Vector     *v,
          const UINT4Vector *indexVector );
void LALU4VectorIndexHole ( LALStatus  *status,
          UINT4VectorPair  *result_pair,
          UINT4Vector      *v,
          const UINT4Vector *indexVector );

/* UINT8 */
void LALU8VectorIndexRange ( LALStatus *status,
          UINT8Vector    **result,
          UINT8Vector     *v,
          const UINT4Vector *indexVector );
void LALU8VectorIndexHole ( LALStatus  *status,
          UINT8VectorPair  *result_pair,
          UINT8Vector      *v,
          const UINT4Vector *indexVector );

/* REAL4 */
void LALSVectorIndexRange ( LALStatus *status,
          REAL4Vector    **result,
          REAL4Vector     *v,
          const UINT4Vector *indexVector );
void LALSVectorIndexHole ( LALStatus  *status,
          REAL4VectorPair  *result_pair,
          REAL4Vector      *v,
          const UINT4Vector *indexVector );

/* REAL8 */
void LALDVectorIndexRange ( LALStatus *status,
          REAL8Vector    **result,
          REAL8Vector     *v,
          const UINT4Vector *indexVector );
void LALDVectorIndexHole ( LALStatus  *status,
          REAL8VectorPair  *result_pair,
          REAL8Vector      *v,
          const UINT4Vector *indexVector );

/* COMPLEX8 */
void LALCVectorIndexRange ( LALStatus *status,
          COMPLEX8Vector    **result,
          COMPLEX8Vector     *v,
          const UINT4Vector *indexVector );
void LALCVectorIndexHole ( LALStatus  *status,
          COMPLEX8VectorPair  *result_pair,
          COMPLEX8Vector      *v,
          const UINT4Vector *indexVector );

/* COMPLEX16 */
void LALZVectorIndexRange ( LALStatus *status,
          COMPLEX16Vector    **result,
          COMPLEX16Vector     *v,
          const UINT4Vector *indexVector );
void LALZVectorIndexHole ( LALStatus  *status,
          COMPLEX16VectorPair  *result_pair,
          COMPLEX16Vector      *v,
          const UINT4Vector *indexVector );

/*@}*/

#ifdef __cplusplus
}
#endif /* C++ */

#endif /* __LALVECTORINDEXRANGE_H__ */
