#ifndef _SEQFACTORIES_H
#define _SEQFACTORIES_H

#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>

#ifdef  __cplusplus
extern "C" {
#endif


/**
 * \addtogroup SeqFactories_h
 *
 * \brief Provides prototype and status code information for use of CreateVectorSequence and DestroyVectorSequence.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/SeqFactories.h>
 * \endcode
 *
 * \section sec_SF_XLALfcts XLAL Functions
 *
 * \code
 * REAL4VectorSequence * XLALCreateVectorSequence(UINT4 length, UINT4 veclen);
 * void XLALCreateVectorSequence(REAL4VectorSequence *vecseq);
 *
 * <vecseqtype> * XLALCreate<vecseqtype>(UINT4 length, UINT4 veclen);
 * void XLALCreate<vecseqtype>(<vecseqtype> *vecseq);
 * \endcode
 *
 * Here <tt>\<vecseqtype\></tt> is one of
 * \c COMPLEX16VectorSequence,
 * \c COMPLEX8VectorSequence,
 * \c REAL8VectorSequence,
 * \c REAL4VectorSequence,
 * \c INT8VectorSequence,
 * \c INT4VectorSequence,
 * \c INT2VectorSequence,
 * \c UINT8VectorSequence,
 * \c UINT4VectorSequence,
 * \c UINT2VectorSequence, or
 * \c CHARVectorSequence.
 *
 * \subsection ss_SF_desc Description
 *
 * The <tt>XLALCreate\<type\>VectorSequence</tt> functions create vector sequences
 * of type <tt>\<type\></tt>, length and vector length \c veclen.
 * The function \c XLALCreateVectorSequence is the same as
 * \c XLALCreateREAL4VectorSequence.
 *
 * The <tt>XLALDestroy\<type\>VectorSequence</tt> functions deallocate the memory
 * allocation pointed to by \c vecseq including its contents.  The function
 * \c XLALDestroyVectorSequence is the same as
 * \c XLALDestroyREAL4VectorSequence.
 *
 * \subsection ss_SF_ret Return Values
 *
 * The create functions return a pointer to the created vector sequence if
 * successful; upon failure they will return \c NULL and set \c xlalErrno
 * to one of the following values: \c #XLAL_ENOMEM if memory allocation
 * failed, or \c #XLAL_EBADLEN if the requested \c length or \c veclen
 * is zero.
 *
 * The destroy functions do not have a return value.  They can fail if they are
 * passed a \c NULL pointer, in which case \c xlalErrno is set to
 * \c #XLAL_EFAULT, or if the vector sequency passed to the destroy routine
 * has zero length, vector length, or \c NULL data pointer then
 * \c xlalErrno is set to \c #XLAL_EINVAL.
 *
 */
/*@{*/

/**\name Error Codes */
/*@{*/
#define SEQFACTORIESH_ESLENGTH  1	/**< Illegal sequence length. */
#define SEQFACTORIESH_EVLENGTH  2	/**< Illegal vector length. */
#define SEQFACTORIESH_EALENGTH  4	/**< Illegal array dimension. */
#define SEQFACTORIESH_EVPTR     8	/**< Null sequence handle. */
#define SEQFACTORIESH_EUPTR    16	/**< Non-null sequence pointer. */
#define SEQFACTORIESH_EDPTR    32	/**< Null sequence data. */
#define SEQFACTORIESH_EINPTR   64	/**< Null input pointer. */
#define SEQFACTORIESH_EMALLOC 128	/**< Malloc failure. */
/*@}*/

/** \cond DONT_DOXYGEN */
#define SEQFACTORIESH_MSGESLENGTH "Illegal sequence length."
#define SEQFACTORIESH_MSGEVLENGTH "Illegal vector length."
#define SEQFACTORIESH_MSGEALENGTH "Illegal array dimension."
#define SEQFACTORIESH_MSGEVPTR    "Null sequence handle."
#define SEQFACTORIESH_MSGEUPTR    "Non-null sequence pointer."
#define SEQFACTORIESH_MSGEDPTR    "Null sequence data."
#define SEQFACTORIESH_MSGEINPTR   "Null input pointer."
#define SEQFACTORIESH_MSGEMALLOC  "Malloc failure."
/** \endcond */

/**
 * \brief This structure stores the input required for creating a vector
 * sequence.  This input includes the length of the sequence (i.e., the number of
 * vectors) and the length of each vector.
 */
typedef struct tagCreateVectorSequenceIn {
  UINT4 length; 	/**< The sequence length */
  UINT4 vectorLength;	/**< The length of each vector in the sequence */
} CreateVectorSequenceIn;


/**
 * \brief This structure stores the input required for creating an array
 * sequence.  This input includes the length of the sequence (i.e., the number of
 * array) and the dimensions of each array index.
 */
typedef struct tagCreateArraySequenceIn {
  UINT4 length;			/**< The sequence length */
  UINT4Vector *dimLength;	/**< The dimensions of each array index (the same for every array in the sequence) */
} CreateArraySequenceIn;

/*@}*/
/* ---------- end:SeqFactories_h ---------- */

/**
 * \defgroup ArraySequenceFactories_c Module ArraySequenceFactories.c
 * \ingroup SeqFactories_h
 *
 * \brief Create/destroy \<datatype\>ArraySequence objects.
 *
 * ### Description ###
 *
 * The \c CreateArraySequence family of functions create a \<datatype\>ArraySequence of the appropriate dimensions.
 *
 * The \c DestroyArraySequence family of functions return the storage allocated by the \c CreateArraySequence functions to the system.
 */
/*@{*/

/** \name REAL4 prototypes (default name) */
/*@{*/
void LALCreateArraySequence(LALStatus *, REAL4ArraySequence **,
                            CreateArraySequenceIn *);
void LALDestroyArraySequence(LALStatus *, REAL4ArraySequence **);
/*@}*/

/** \name INT2 prototypes */
/*@{*/
void LALI2CreateArraySequence ( LALStatus *status,
          INT2ArraySequence **arraySequence,
          CreateArraySequenceIn *aSeqParams);

void LALI2DestroyArraySequence ( LALStatus *status,
          INT2ArraySequence **arraySeqence);
/*@}*/


/** \name INT4 prototypes */
/*@{*/
void LALI4CreateArraySequence ( LALStatus *status,
          INT4ArraySequence **arraySequence,
          CreateArraySequenceIn *aSeqParams);

void LALI4DestroyArraySequence ( LALStatus *status,
          INT4ArraySequence **arraySeqence);
/*@}*/

/** \name INT8 prototypes */
/*@{*/
void LALI8CreateArraySequence ( LALStatus *status,
          INT8ArraySequence **arraySequence,
          CreateArraySequenceIn *aSeqParams);

void LALI8DestroyArraySequence ( LALStatus *status,
          INT8ArraySequence **arraySeqence);
/*@}*/

/** \name UINT2 prototypes */
/*@{*/
void LALU2CreateArraySequence ( LALStatus *status,
          UINT2ArraySequence **arraySequence,
          CreateArraySequenceIn *aSeqParams);

void LALU2DestroyArraySequence ( LALStatus *status,
          UINT2ArraySequence **arraySeqence);
/*@}*/

/** \name UINT4 prototypes */
/*@{*/
void LALU4CreateArraySequence ( LALStatus *status,
          UINT4ArraySequence **arraySequence,
          CreateArraySequenceIn *aSeqParams);

void LALU4DestroyArraySequence ( LALStatus *status,
          UINT4ArraySequence **arraySeqence);
/*@}*/

/** \name UINT8 prototypes */
/*@{*/
void LALU8CreateArraySequence ( LALStatus *status,
          UINT8ArraySequence **arraySequence,
          CreateArraySequenceIn *aSeqParams);

void LALU8DestroyArraySequence ( LALStatus *status,
          UINT8ArraySequence **arraySeqence);
/*@}*/

/** \name REAL4 prototypes */
/*@{*/
void LALSCreateArraySequence ( LALStatus *status,
          REAL4ArraySequence **arraySequence,
          CreateArraySequenceIn *aSeqParams);

void LALSDestroyArraySequence ( LALStatus *status,
          REAL4ArraySequence **arraySeqence);
/*@}*/

/** \name REAL8 prototypes */
/*@{*/
void LALDCreateArraySequence ( LALStatus *status,
          REAL8ArraySequence **arraySequence,
          CreateArraySequenceIn *aSeqParams);

void LALDDestroyArraySequence ( LALStatus *status,
          REAL8ArraySequence **arraySeqence);
/*@}*/

/** \name COMPLEX8 prototypes */
/*@{*/
void LALCCreateArraySequence ( LALStatus *status,
          COMPLEX8ArraySequence **arraySequence,
          CreateArraySequenceIn *aSeqParams);

void LALCDestroyArraySequence ( LALStatus *status,
          COMPLEX8ArraySequence **arraySeqence);
/*@}*/

/** \name COMPLEX16 prototypes */
/*@{*/
void LALZCreateArraySequence ( LALStatus *status,
          COMPLEX16ArraySequence **arraySequence,
          CreateArraySequenceIn *aSeqParams);

void LALZDestroyArraySequence ( LALStatus *status,
          COMPLEX16ArraySequence **arraySeqence);
/*@}*/

/*@}*/
/* ---------- end: ArraySequenceFactories_c ---------- */

/**
 * \defgroup VectorSequenceFactories_c Module VectorSequenceFactories.c
 * \ingroup SeqFactories_h
 * \brief Create/destroy \<datatype\>VectorSequence objects.
 *
 * ### Description ###
 *
 * The \c CreateVectorSequence family of functions create a \<datatype\>VectorSequence of the appropriate dimensions.
 *
 * The \c DestroyVectorSequence family of functions return the storage allocated by the \c CreateVectorSequence functions to the system.
 */
/*@{*/

/** \name REAL4 prototypes (default name) */
/*@{*/
REAL4VectorSequence * XLALCreateVectorSequence ( UINT4 length, UINT4 veclen );
#ifndef SWIG   /* exclude from SWIG interface */
void XLALDestroyVectorSequence ( REAL4VectorSequence * vecseq );
#endif   /* SWIG */


void LALCreateVectorSequence(LALStatus *, REAL4VectorSequence **,
                             CreateVectorSequenceIn *);
void LALDestroyVectorSequence(LALStatus *, REAL4VectorSequence**);
/*@}*/

/** \name CHAR prototypes */
/*@{*/
CHARVectorSequence * XLALCreateCHARVectorSequence ( UINT4 length, UINT4 veclen );
void XLALDestroyCHARVectorSequence ( CHARVectorSequence * vecseq );

void LALCHARCreateVectorSequence ( LALStatus *status,
          CHARVectorSequence **vectorSequence,
          CreateVectorSequenceIn *vSeqParams);

void LALCHARDestroyVectorSequence ( LALStatus *status,
          CHARVectorSequence **vectorSequence);
/*@}*/

/** \name INT2 prototypes */
/*@{*/
INT2VectorSequence * XLALCreateINT2VectorSequence ( UINT4 length, UINT4 veclen );
void XLALDestroyINT2VectorSequence ( INT2VectorSequence * vecseq );

void LALI2CreateVectorSequence ( LALStatus *status,
          INT2VectorSequence **vectorSequence,
          CreateVectorSequenceIn *vSeqParams);

void LALI2DestroyVectorSequence ( LALStatus *status,
          INT2VectorSequence **vectorSequence);
/*@}*/

/** \name INT4 prototypes */
/*@{*/
INT4VectorSequence * XLALCreateINT4VectorSequence ( UINT4 length, UINT4 veclen );
void XLALDestroyINT4VectorSequence ( INT4VectorSequence * vecseq );

void LALI4CreateVectorSequence ( LALStatus *status,
          INT4VectorSequence **vectorSequence,
          CreateVectorSequenceIn *vSeqParams);

void LALI4DestroyVectorSequence ( LALStatus *status,
          INT4VectorSequence **vectorSequence);
/*@}*/

/** \name INT8 prototypes */
/*@{*/
INT8VectorSequence * XLALCreateINT8VectorSequence ( UINT4 length, UINT4 veclen );
void XLALDestroyINT8VectorSequence ( INT8VectorSequence * vecseq );

void LALI8CreateVectorSequence ( LALStatus *status,
          INT8VectorSequence **vectorSequence,
          CreateVectorSequenceIn *vSeqParams);

void LALI8DestroyVectorSequence ( LALStatus *status,
          INT8VectorSequence **vectorSequence);
/*@}*/

/** \name UINT2 prototypes */
/*@{*/
UINT2VectorSequence * XLALCreateUINT2VectorSequence ( UINT4 length, UINT4 veclen );
void XLALDestroyUINT2VectorSequence ( UINT2VectorSequence * vecseq );

void LALU2CreateVectorSequence ( LALStatus *status,
          UINT2VectorSequence **vectorSequence,
          CreateVectorSequenceIn *vSeqParams);

void LALU2DestroyVectorSequence ( LALStatus *status,
          UINT2VectorSequence **vectorSequence);
/*@}*/

/** \name UINT4 prototypes */
/*@{*/
UINT4VectorSequence * XLALCreateUINT4VectorSequence ( UINT4 length, UINT4 veclen );
void XLALDestroyUINT4VectorSequence ( UINT4VectorSequence * vecseq );

void LALU4CreateVectorSequence ( LALStatus *status,
          UINT4VectorSequence **vectorSequence,
          CreateVectorSequenceIn *vSeqParams);

void LALU4DestroyVectorSequence ( LALStatus *status,
          UINT4VectorSequence **vectorSequence);
/*@}*/

/** \name UINT8 prototypes */
/*@{*/
UINT8VectorSequence * XLALCreateUINT8VectorSequence ( UINT4 length, UINT4 veclen );
void XLALDestroyUINT8VectorSequence ( UINT8VectorSequence * vecseq );

void LALU8CreateVectorSequence ( LALStatus *status,
          UINT8VectorSequence **vectorSequence,
          CreateVectorSequenceIn *vSeqParams);

void LALU8DestroyVectorSequence ( LALStatus *status,
          UINT8VectorSequence **vectorSequence);
/*@}*/

/** \name REAL4 prototypes */
/*@{*/
REAL4VectorSequence * XLALCreateREAL4VectorSequence ( UINT4 length, UINT4 veclen );
void XLALDestroyREAL4VectorSequence ( REAL4VectorSequence * vecseq );

void LALSCreateVectorSequence ( LALStatus *status,
          REAL4VectorSequence **vectorSequence,
          CreateVectorSequenceIn *vSeqParams);

void LALSDestroyVectorSequence ( LALStatus *status,
          REAL4VectorSequence **vectorSequence);
/*@}*/

/** \name REAL8 prototypes */
/*@{*/
REAL8VectorSequence * XLALCreateREAL8VectorSequence ( UINT4 length, UINT4 veclen );
void XLALDestroyREAL8VectorSequence ( REAL8VectorSequence * vecseq );

void LALDCreateVectorSequence ( LALStatus *status,
          REAL8VectorSequence **vectorSequence,
          CreateVectorSequenceIn *vSeqParams);

void LALDDestroyVectorSequence ( LALStatus *status,
          REAL8VectorSequence **vectorSequence);
/*@}*/

/** \name COMPLEX8 prototypes */
/*@{*/
COMPLEX8VectorSequence * XLALCreateCOMPLEX8VectorSequence ( UINT4 length, UINT4 veclen );
void XLALDestroyCOMPLEX8VectorSequence ( COMPLEX8VectorSequence * vecseq );

void LALCCreateVectorSequence ( LALStatus *status,
          COMPLEX8VectorSequence **vectorSequence,
          CreateVectorSequenceIn *vSeqParams);

void LALCDestroyVectorSequence ( LALStatus *status,
          COMPLEX8VectorSequence **vectorSequence);
/*@}*/

/** \name COMPLEX16 prototypes */
/*@{*/
COMPLEX16VectorSequence * XLALCreateCOMPLEX16VectorSequence ( UINT4 length, UINT4 veclen );
void XLALDestroyCOMPLEX16VectorSequence ( COMPLEX16VectorSequence * vecseq );

void LALZCreateVectorSequence ( LALStatus *status,
          COMPLEX16VectorSequence **vectorSequence,
          CreateVectorSequenceIn *vSeqParams);

void LALZDestroyVectorSequence ( LALStatus *status,
          COMPLEX16VectorSequence **vectorSequence);
/*@}*/

/*@}*/
/* ---------- end: VectorSequenceFactories_c ---------- */

#ifdef  __cplusplus
}
#endif

#endif /* _SEQFACTORIES_H */
