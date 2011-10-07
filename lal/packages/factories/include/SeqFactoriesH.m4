/**
\defgroup SeqFactories_h SeqFactories_h
\ingroup factories

\brief Provides prototype and status code information for use of CreateVectorSequence
and DestroyVectorSequence.

\heading{Synopsis}
\code
#include <lal/SeqFactories.h>
\endcode

\section sec_SF_XLALfcts XLAL Functions

\code
REAL4VectorSequence * XLALCreateVectorSequence(UINT4 length, UINT4 veclen);
void XLALCreateVectorSequence(REAL4VectorSequence *vecseq);

<vecseqtype> * XLALCreate<vecseqtype>(UINT4 length, UINT4 veclen);
void XLALCreate<vecseqtype>(<vecseqtype> *vecseq);
\endcode

Here <tt>\<vecseqtype\></tt> is one of
\c COMPLEX16VectorSequence,
\c COMPLEX8VectorSequence,
\c REAL8VectorSequence,
\c REAL4VectorSequence,
\c INT8VectorSequence,
\c INT4VectorSequence,
\c INT2VectorSequence,
\c UINT8VectorSequence,
\c UINT4VectorSequence,
\c UINT2VectorSequence, or
\c CHARVectorSequence.

\subsection ss_SF_desc Description

The <tt>XLALCreate\<type\>VectorSequence</tt> functions create vector sequences
of type <tt>\<type\></tt>, length and vector length \c veclen.
The function \c XLALCreateVectorSequence is the same as
\c XLALCreateREAL4VectorSequence.

The <tt>XLALDestroy\<type\>VectorSequence</tt> functions deallocate the memory
allocation pointed to by \c vecseq including its contents.  The function
\c XLALDestroyVectorSequence is the same as
\c XLALDestroyREAL4VectorSequence.

\subsection ss_SF_ret Return Values

The create functions return a pointer to the created vector sequence if
successful; upon failure they will return \c NULL and set \c xlalErrno
to one of the following values: \c #XLAL_ENOMEM if memory allocation
failed, or \c #XLAL_EBADLEN if the requested \c length or \c veclen
is zero.

The destroy functions do not have a return value.  They can fail if they are
passed a \c NULL pointer, in which case \c xlalErrno is set to
\c #XLAL_EFAULT, or if the vector sequency passed to the destroy routine
has zero length, vector length, or \c NULL data pointer then
\c xlalErrno is set to \c #XLAL_EINVAL.

*/


#ifndef _SEQFACTORIES_H
#define _SEQFACTORIES_H

/* remove SWIG interface directives */
#if !defined(SWIG) && !defined(SWIGLAL_STRUCT)
#define SWIGLAL_STRUCT(...)
#endif

#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (SEQFACTORIESH, "$Id$");

/**\name Error Codes *//*@{*/
/** \ingroup SeqFactories_h */
#define SEQFACTORIESH_ESLENGTH  1
#define SEQFACTORIESH_EVLENGTH  2
#define SEQFACTORIESH_EALENGTH  4
#define SEQFACTORIESH_EVPTR     8
#define SEQFACTORIESH_EUPTR    16
#define SEQFACTORIESH_EDPTR    32
#define SEQFACTORIESH_EINPTR   64
#define SEQFACTORIESH_EMALLOC 128

#define SEQFACTORIESH_MSGESLENGTH "Illegal sequence length."
#define SEQFACTORIESH_MSGEVLENGTH "Illegal vector length."
#define SEQFACTORIESH_MSGEALENGTH "Illegal array dimension."
#define SEQFACTORIESH_MSGEVPTR    "Null sequence handle."
#define SEQFACTORIESH_MSGEUPTR    "Non-null sequence pointer."
#define SEQFACTORIESH_MSGEDPTR    "Null sequence data."
#define SEQFACTORIESH_MSGEINPTR   "Null input pointer."
#define SEQFACTORIESH_MSGEMALLOC  "Malloc failure."
/*@}*/




/** \ingroup SeqFactories_h
 * \brief This structure stores the input required for creating a vector
 * sequence.  This input includes the length of the sequence (i.e., the number of
 * vectors) and the length of each vector.
 */
typedef struct tagCreateVectorSequenceIn {
  SWIGLAL_STRUCT(CreateVectorSequenceIn);
  UINT4 length; 	/**< The sequence length */
  UINT4 vectorLength;	/**< The length of each vector in the sequence */
} CreateVectorSequenceIn;


/** \ingroup SeqFactories_h
 * \brief This structure stores the input required for creating an array
 * sequence.  This input includes the length of the sequence (i.e., the number of
 * array) and the dimensions of each array index.
 */
typedef struct tagCreateArraySequenceIn {
  SWIGLAL_STRUCT(CreateArraySequenceIn);
  UINT4 length;			/**< The sequence length */
  UINT4Vector *dimLength;	/**< The dimensions of each array index (the same for every array in the sequence) */
} CreateArraySequenceIn;


REAL4VectorSequence * XLALCreateVectorSequence ( UINT4 length, UINT4 veclen );
void XLALDestroyVectorSequence ( REAL4VectorSequence * vecseq );


void LALCreateVectorSequence(LALStatus *, REAL4VectorSequence **,
                             CreateVectorSequenceIn *);
void LALDestroyVectorSequence(LALStatus *, REAL4VectorSequence**);

void LALCreateArraySequence(LALStatus *, REAL4ArraySequence **,
                            CreateArraySequenceIn *);
void LALDestroyArraySequence(LALStatus *, REAL4ArraySequence **);

define(`TYPECODE',`CHAR')
include(`SeqFactoriesBaseH.m4')

define(`TYPECODE',`I2')
include(`SeqFactoriesBaseH.m4')
include(`ArraySeqFactoriesBaseH.m4')

define(`TYPECODE',`I4')
include(`SeqFactoriesBaseH.m4')
include(`ArraySeqFactoriesBaseH.m4')

define(`TYPECODE',`I8')
include(`SeqFactoriesBaseH.m4')
include(`ArraySeqFactoriesBaseH.m4')

define(`TYPECODE',`U2')
include(`SeqFactoriesBaseH.m4')
include(`ArraySeqFactoriesBaseH.m4')

define(`TYPECODE',`U4')
include(`SeqFactoriesBaseH.m4')
include(`ArraySeqFactoriesBaseH.m4')

define(`TYPECODE',`U8')
include(`SeqFactoriesBaseH.m4')
include(`ArraySeqFactoriesBaseH.m4')

define(`TYPECODE',`S')
include(`SeqFactoriesBaseH.m4')
include(`ArraySeqFactoriesBaseH.m4')

define(`TYPECODE',`D')
include(`SeqFactoriesBaseH.m4')
include(`ArraySeqFactoriesBaseH.m4')

define(`TYPECODE',`C')
include(`SeqFactoriesBaseH.m4')
include(`ArraySeqFactoriesBaseH.m4')

define(`TYPECODE',`Z')
include(`SeqFactoriesBaseH.m4')
include(`ArraySeqFactoriesBaseH.m4')

#ifdef  __cplusplus
}
#endif

#endif /* _SEQFACTORIES_H */
