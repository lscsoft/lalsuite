/* -*- C -*- */
/**
 * \author David Chin <dwchin@umich.edu> +1-734-709-9119
 * \addtogroup VectorIndexRange_h
 *
 * \section TODOref Header \ref VectorIndexRange.h
 *
 * \brief Routines to slice up vectors/sequences by specifying starting and ending indices, inclusive.
 *
 * VectorIndexHole() is "complementary" to VectorIndexRange(), in the sense that
 * it returns the two segments of the vector that would remain after VectorIndexRange() acts.
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/VectorIndexRange.h>
 * \endcode
 *
 * @{
 */

#include <math.h>
#include <lal/LALRCSID.h>
#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>
#include <lal/AVFactories.h>

#ifndef __LALVECTORINDEXRANGE_H__
#define __LALVECTORINDEXRANGE_H__

/* remove SWIG interface directives */
#if !defined(SWIG) && !defined(SWIGLAL_STRUCT_LALALLOC)
#define SWIGLAL_STRUCT_LALALLOC
#endif

#ifdef __cplusplus
extern "C" {
#endif /* C++ */

/**\name Error Codes */ /*@{*/
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
/** @}*/

/** @}*/

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


NRCSID( VECTORINDEXRANGEH, "$Id$" );

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


