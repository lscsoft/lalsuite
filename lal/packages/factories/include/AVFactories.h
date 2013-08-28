/*-----------------------------------------------------------------------

File Name: AVFactories.h

-------------------------------------------------------------------------*/

#ifndef _AVFACTORIES_H
#define _AVFACTORIES_H

#include <lal/LALDatatypes.h>
#include <stdarg.h>

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \addtogroup AVFactories_h
 *
 * \brief Provides prototype and status code information for use of CreateVector,
 * CreateArray, ResizeVector, ResizeArray, DestroyVector and DestroyArray
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/AVFactories.h>
 * \endcode
 *
 * \section secXLALfcts XLAL Functions
 *
 * \code
 * REAL4Vector * XLALCreateVector(UINT4 length);
 * REAL4Vector * XLALResizeVector(REAL4Vector *vector, UINT4 length);
 * void XLALDestroyVector(REAL4Vector *vector, UINT4 length);
 *
 * <vectype> * XLALCreate<vectype>(UINT4 length );
 * <vectype> * XLALResize<vectype>(<vectype> *vector, UINT4 length );
 * void XLALDestroy<vectype>(<vectype> *vector);
 *
 * REAL4Array * XLALCreateArrayL(UINT4 ndim, ...);
 * REAL4Array * XLALCreateArrayV(UINT4 ndim, UINT4 *dims);
 * REAL4Array * XLALCreateArray(UINT4Vector *dimLength);
 * REAL4Array * XLALResizeArrayL(REAL4Array *array, UINT4 ndim, ...);
 * REAL4Array * XLALResizeArrayV(REAL4Array *array, UINT4 ndim, UINT4 *dims);
 * REAL4Array * XLALResizeArray(REAL4Array *array, UINT4Vector *dimLength);
 * void XLALDestroyArray(REAL4Array *array);
 *
 * <arrtype> * XLALCreate<arrtype>L(UINT4 ndim, ...);
 * <arrtype> * XLALCreate<arrtype>V(UINT4 ndim, UINT4 *dims);
 * <arrtype> * XLALCreate<arrtype>(UINT4Vector *dimLength);
 * <arrtype> * XLALResize<arrtype>L(<arrtype> *array, UINT4 ndim, ...);
 * <arrtype> * XLALResize<arrtype>V(<arrtype> *array, UINT4 ndim, UINT4 *dims);
 * <arrtype> * XLALResize<arrtype>(<arrtype> *array, UINT4Vector *dimLength);
 * void XLALDestroy<arrtype>(<arrtype> *array);
 * \endcode
 *
 * Here <tt>\<vectype\></tt> is one of
 * \c COMPLEX16Vector,
 * \c COMPLEX8Vector,
 * \c REAL8Vector,
 * \c REAL4Vector,
 * \c INT8Vector,
 * \c INT4Vector,
 * \c INT2Vector,
 * \c UINT8Vector,
 * \c UINT4Vector,
 * \c UINT2Vector, or
 * \c CHARVector,
 * and <tt>\<arrtype\></tt> is one of
 * \c COMPLEX16Array,
 * \c COMPLEX8Array,
 * \c REAL8Array,
 * \c REAL4Array,
 * \c INT8Array,
 * \c INT4Array,
 * \c INT2Array,
 * \c UINT8Array,
 * \c UINT4Array, or
 * \c UINT2Array.
 *
 * \subsection ss_AVF_desc Description
 *
 * The <tt>XLALCreate\<type\>%Vector</tt> functions create vectors of the specified
 * \c length number of objects of type <tt>\<type\></tt>.  The function
 * \c XLALCreateVector() is the same as \c XLALCreateREAL4Vector().
 *
 * The <tt>XLALDestroy\<type\>%Vector</tt> functions deallocate the memory allocation
 * pointed to by \c vector including its contents.  The function
 * \c XLALDestroyVector() is the same as \c XLALDestroyREAL4Vector().
 *
 * The <tt>XLALResize\<type\>%Vector</tt> functions resize the supplied vector
 * \c vector to the new size \c length.  If \c vector is \c NULL
 * then this is equivalent to <tt>XLALCreate\<type\>%Vector</tt>.  If \c length
 * is zero then this is equivalent to <tt>XLALDestroy\<type\>%Vector</tt> and the
 * routine returns \c NULL.  Otherwise, the amount of data in the vector
 * is realloced using \c LALRealloc().  The function
 * \c XLALResizeVector() is the same as \c XLALResizeREAL4Vector().
 *
 * The <tt>XLALCreate\<type\>Array</tt>
 * <tt>XLALCreate\<type\>ArrayL</tt>
 * <tt>XLALCreate\<type\>ArrayV</tt>
 * all create an object of type <tt>\<type\>Array</tt>.  They differ in the way
 * that the dimensions of the array are specified.  The function
 * <tt>XLALCreate\<type\>Array</tt> allocates an array with dimensions specified
 * by the \c UINT4Vector \c dimLength which is a vector of dimension
 * lengths for the array.  The function
 * <tt>XLALCreate\<type\>ArrayV</tt> provides these dimensions with two arguments:
 * \c ndim is the number of dimensions and \c dims is an array of
 * \c UINT4 values for the dimensions.  The function
 * <tt>XLALCreate\<type\>ArrayL</tt> also specifies the dimensions as arguments.
 * Here, the first argument, \c ndim, is the number of dimensions, and
 * this is followed by \c ndim arguments that provide the dimensions.
 * Note that for this function, a maximum of 16 dimensions can be provided
 * (that is, \c ndim cannot be more than 16 and there cannot be more than
 * 16 arguments after the first).
 * The \c XLALCreateArray()
 * \c XLALCreateArrayL()
 * \c XLALCreateArrayV()
 * functions are equivalent to the
 * \c XLALCreateREAL4Array()
 * \c XLALCreateREAL4ArrayL()
 * \c XLALCreateREAL4ArrayV()
 * functions respectively.
 *
 * The <tt>XLALDestroy\<type\>Array</tt> functions deallocate the memory allocation
 * pointed to by \c array including its contents.  The function
 * \c XLALDestroyArray() is the same as \c XLALDestroyREAL4Array().
 *
 * The <tt>XLALResize\<type\>Array</tt>
 * <tt>XLALResize\<type\>ArrayL</tt>
 * <tt>XLALResize\<type\>ArrayV</tt>
 * functions resize the provided array \c array.  The arguments after the
 * first are interpreted in the same way as for the
 * <tt>XLALCreate\<type\>Array</tt>
 * <tt>XLALCreate\<type\>ArrayL</tt>
 * <tt>XLALCreate\<type\>ArrayV</tt>
 * functions.  If \c array is \c NULL, the resize functions are equivalent
 * to the corresponding create function.  If \c ndim is zero for
 * <tt>XLALResize\<type\>ArrayL</tt> or
 * <tt>XLALResize\<type\>ArrayV</tt>, or if \c dimLength is \c NULL for
 * <tt>XLALResize\<type\>Array</tt>, then these functions are equivalent to
 * <tt>XLALDestroy\<type\>Array</tt>.
 * The \c XLALResizeArray()
 * \c XLALResizeArrayL()
 * \c XLALResizeArrayV()
 * functions are equivalent to the
 * \c XLALResizeREAL4Array()
 * \c XLALResizeREAL4ArrayL()
 * \c XLALResizeREAL4ArrayV()
 * functions respectively.
 *
 * \subsection ss_AVF_return Return Values
 *
 * If successful, the create and resize functions return a pointer to the same
 * data that was passed to the function.  The resize functions will return
 * a \c NULL pointer if the size of the new object was zero.  Upon failure
 * these routines will return \c NULL and will set \c xlalErrno to
 * one of these values: \c #XLAL_ENOMEM if a memory allocation failed,
 * \c #XLAL_EBADLEN if an invalid length was provided (for example, a
 * zero-size allocation with a create function), \c XLAL_EINVAL if an
 * invalid argument is provided (for example, if the pointer to an
 * array of dimensions is \c NULL).
 *
 * The destroy function does not return any value.  If the function is passed a
 * \c NULL pointer, it will set \c xlalErrno to \c #XLAL_EFAULT.
 * If the function is passed an object that appears corrupted (e.g., a
 * vector with zero length or will a \c NULL data pointer) it will set
 * \c xlalErrno to \c #XLAL_EINVAL.
 *
 */
/*@{ */

/** \name Error Codes */
/*@{*/
#define AVFACTORIESH_ELENGTH 1	/**< Illegal length. */
#define AVFACTORIESH_EVPTR   2	/**< Null vector/array handle. */
#define AVFACTORIESH_EUPTR   4	/**< Non-null vector/array pointer. */
#define AVFACTORIESH_EDPTR   8	/**< Null vector/array data. */
#define AVFACTORIESH_EMALLOC 16	/**< Malloc failure. */
/*@}*/
/*@} */

/** \cond DONT_DOXYGEN */
#define AVFACTORIESH_MSGELENGTH  "Illegal length."
#define AVFACTORIESH_MSGEVPTR    "Null vector/array handle."
#define AVFACTORIESH_MSGEUPTR    "Non-null vector/array pointer."
#define AVFACTORIESH_MSGEDPTR    "Null vector/array data."
#define AVFACTORIESH_MSGEMALLOC  "Malloc failure."
/** \endcond */


/**
 * \defgroup ArrayFactories_c Module ArrayFactories.c
 * \ingroup AVFactories_h
 *
 * \brief Create/destroy \<datatype\>Array objects.
 *
 * \heading{Description}
 *
 * The \c CreateArray family of functions create a \<datatype\>Array of the appropriate dimensions.
 *
 * The \c DestroyArray family of functions return the storage allocated by the \c CreateArray functions to the system.
 *
 */
/*@{*/

/** \name INT2 array prototypes */
/*@{*/
INT2Array * XLALCreateINT2ArrayL ( UINT4, ... );
INT2Array * XLALCreateINT2ArrayV ( UINT4, UINT4 * );
INT2Array * XLALCreateINT2Array ( UINT4Vector * );
INT2Array * XLALResizeINT2ArrayL ( INT2Array *, UINT4, ... );
INT2Array * XLALResizeINT2ArrayV ( INT2Array *, UINT4, UINT4 * );
INT2Array * XLALResizeINT2Array ( INT2Array *, UINT4Vector * );
void XLALDestroyINT2Array ( INT2Array * );
void LALI2CreateArray ( LALStatus *, INT2Array **, UINT4Vector * );
void LALI2ResizeArray ( LALStatus *, INT2Array **, UINT4Vector * );
void LALI2DestroyArray ( LALStatus *, INT2Array ** );
/*@}*/

/** \name INT4 array prototypes */
/*@{*/
INT4Array * XLALCreateINT4ArrayL ( UINT4, ... );
INT4Array * XLALCreateINT4ArrayV ( UINT4, UINT4 * );
INT4Array * XLALCreateINT4Array ( UINT4Vector * );
INT4Array * XLALResizeINT4ArrayL ( INT4Array *, UINT4, ... );
INT4Array * XLALResizeINT4ArrayV ( INT4Array *, UINT4, UINT4 * );
INT4Array * XLALResizeINT4Array ( INT4Array *, UINT4Vector * );
void XLALDestroyINT4Array ( INT4Array * );
void LALI4CreateArray ( LALStatus *, INT4Array **, UINT4Vector * );
void LALI4ResizeArray ( LALStatus *, INT4Array **, UINT4Vector * );
void LALI4DestroyArray ( LALStatus *, INT4Array ** );
/*@}*/

/** \name INT8 array prototypes */
/*@{*/
INT8Array * XLALCreateINT8ArrayL ( UINT4, ... );
INT8Array * XLALCreateINT8ArrayV ( UINT4, UINT4 * );
INT8Array * XLALCreateINT8Array ( UINT4Vector * );
INT8Array * XLALResizeINT8ArrayL ( INT8Array *, UINT4, ... );
INT8Array * XLALResizeINT8ArrayV ( INT8Array *, UINT4, UINT4 * );
INT8Array * XLALResizeINT8Array ( INT8Array *, UINT4Vector * );
void XLALDestroyINT8Array ( INT8Array * );
void LALI8CreateArray ( LALStatus *, INT8Array **, UINT4Vector * );
void LALI8ResizeArray ( LALStatus *, INT8Array **, UINT4Vector * );
void LALI8DestroyArray ( LALStatus *, INT8Array ** );
/*@}*/

/** \name UINT2 array prototypes */
/*@{*/
UINT2Array * XLALCreateUINT2ArrayL ( UINT4, ... );
UINT2Array * XLALCreateUINT2ArrayV ( UINT4, UINT4 * );
UINT2Array * XLALCreateUINT2Array ( UINT4Vector * );
UINT2Array * XLALResizeUINT2ArrayL ( UINT2Array *, UINT4, ... );
UINT2Array * XLALResizeUINT2ArrayV ( UINT2Array *, UINT4, UINT4 * );
UINT2Array * XLALResizeUINT2Array ( UINT2Array *, UINT4Vector * );
void XLALDestroyUINT2Array ( UINT2Array * );
void LALU2CreateArray ( LALStatus *, UINT2Array **, UINT4Vector * );
void LALU2ResizeArray ( LALStatus *, UINT2Array **, UINT4Vector * );
void LALU2DestroyArray ( LALStatus *, UINT2Array ** );
/*@}*/

/** \name UINT4 array prototypes */
/*@{*/
UINT4Array * XLALCreateUINT4ArrayL ( UINT4, ... );
UINT4Array * XLALCreateUINT4ArrayV ( UINT4, UINT4 * );
UINT4Array * XLALCreateUINT4Array ( UINT4Vector * );
UINT4Array * XLALResizeUINT4ArrayL ( UINT4Array *, UINT4, ... );
UINT4Array * XLALResizeUINT4ArrayV ( UINT4Array *, UINT4, UINT4 * );
UINT4Array * XLALResizeUINT4Array ( UINT4Array *, UINT4Vector * );
void XLALDestroyUINT4Array ( UINT4Array * );
void LALU4CreateArray ( LALStatus *, UINT4Array **, UINT4Vector * );
void LALU4ResizeArray ( LALStatus *, UINT4Array **, UINT4Vector * );
void LALU4DestroyArray ( LALStatus *, UINT4Array ** );
/*@}*/

/** \name UINT8 array prototypes */
/*@{*/
UINT8Array * XLALCreateUINT8ArrayL ( UINT4, ... );
UINT8Array * XLALCreateUINT8ArrayV ( UINT4, UINT4 * );
UINT8Array * XLALCreateUINT8Array ( UINT4Vector * );
UINT8Array * XLALResizeUINT8ArrayL ( UINT8Array *, UINT4, ... );
UINT8Array * XLALResizeUINT8ArrayV ( UINT8Array *, UINT4, UINT4 * );
UINT8Array * XLALResizeUINT8Array ( UINT8Array *, UINT4Vector * );
void XLALDestroyUINT8Array ( UINT8Array * );
void LALU8CreateArray ( LALStatus *, UINT8Array **, UINT4Vector * );
void LALU8ResizeArray ( LALStatus *, UINT8Array **, UINT4Vector * );
void LALU8DestroyArray ( LALStatus *, UINT8Array ** );
/*@}*/


/** \name REAL4 array prototypes */
/*@{*/
REAL4Array * XLALCreateREAL4ArrayL ( UINT4, ... );
REAL4Array * XLALCreateREAL4ArrayV ( UINT4, UINT4 * );
REAL4Array * XLALCreateREAL4Array ( UINT4Vector * );
REAL4Array * XLALResizeREAL4ArrayL ( REAL4Array *, UINT4, ... );
REAL4Array * XLALResizeREAL4ArrayV ( REAL4Array *, UINT4, UINT4 * );
REAL4Array * XLALResizeREAL4Array ( REAL4Array *, UINT4Vector * );
void XLALDestroyREAL4Array ( REAL4Array * );
void LALSCreateArray ( LALStatus *, REAL4Array **, UINT4Vector * );
void LALSResizeArray ( LALStatus *, REAL4Array **, UINT4Vector * );
void LALSDestroyArray ( LALStatus *, REAL4Array ** );
/*@}*/

/** \name REAL4 array prototypes (default name) */
/*@{*/
REAL4Array * XLALCreateArrayL ( UINT4, ... );
REAL4Array * XLALCreateArrayV ( UINT4, UINT4 * );
REAL4Array * XLALCreateArray ( UINT4Vector * );
REAL4Array * XLALResizeArrayL ( REAL4Array *, UINT4, ... );
REAL4Array * XLALResizeArrayV ( REAL4Array *, UINT4, UINT4 * );
REAL4Array * XLALResizeArray ( REAL4Array *, UINT4Vector * );
void XLALDestroyArray ( REAL4Array * );
void LALCreateArray ( LALStatus *, REAL4Array **, UINT4Vector * );
void LALResizeArray ( LALStatus *, REAL4Array **, UINT4Vector * );
void LALDestroyArray ( LALStatus *, REAL4Array ** );
/*@}*/

/** \name REAL8 array prototypes */
/*@{*/
REAL8Array * XLALCreateREAL8ArrayL ( UINT4, ... );
REAL8Array * XLALCreateREAL8ArrayV ( UINT4, UINT4 * );
REAL8Array * XLALCreateREAL8Array ( UINT4Vector * );
REAL8Array * XLALResizeREAL8ArrayL ( REAL8Array *, UINT4, ... );
REAL8Array * XLALResizeREAL8ArrayV ( REAL8Array *, UINT4, UINT4 * );
REAL8Array * XLALResizeREAL8Array ( REAL8Array *, UINT4Vector * );
void XLALDestroyREAL8Array ( REAL8Array * );
void LALDCreateArray ( LALStatus *, REAL8Array **, UINT4Vector * );
void LALDResizeArray ( LALStatus *, REAL8Array **, UINT4Vector * );
void LALDDestroyArray ( LALStatus *, REAL8Array ** );
/*@}*/

/** \name COMPLEX8 array prototypes */
/*@{*/
COMPLEX8Array * XLALCreateCOMPLEX8ArrayL ( UINT4, ... );
COMPLEX8Array * XLALCreateCOMPLEX8ArrayV ( UINT4, UINT4 * );
COMPLEX8Array * XLALCreateCOMPLEX8Array ( UINT4Vector * );
COMPLEX8Array * XLALResizeCOMPLEX8ArrayL ( COMPLEX8Array *, UINT4, ... );
COMPLEX8Array * XLALResizeCOMPLEX8ArrayV ( COMPLEX8Array *, UINT4, UINT4 * );
COMPLEX8Array * XLALResizeCOMPLEX8Array ( COMPLEX8Array *, UINT4Vector * );
void XLALDestroyCOMPLEX8Array ( COMPLEX8Array * );
void LALCCreateArray ( LALStatus *, COMPLEX8Array **, UINT4Vector * );
void LALCResizeArray ( LALStatus *, COMPLEX8Array **, UINT4Vector * );
void LALCDestroyArray ( LALStatus *, COMPLEX8Array ** );
/*@}*/

/** \name COMPLEX16 array prototypes */
/*@{*/
COMPLEX16Array * XLALCreateCOMPLEX16ArrayL ( UINT4, ... );
COMPLEX16Array * XLALCreateCOMPLEX16ArrayV ( UINT4, UINT4 * );
COMPLEX16Array * XLALCreateCOMPLEX16Array ( UINT4Vector * );
COMPLEX16Array * XLALResizeCOMPLEX16ArrayL ( COMPLEX16Array *, UINT4, ... );
COMPLEX16Array * XLALResizeCOMPLEX16ArrayV ( COMPLEX16Array *, UINT4, UINT4 * );
COMPLEX16Array * XLALResizeCOMPLEX16Array ( COMPLEX16Array *, UINT4Vector * );
void XLALDestroyCOMPLEX16Array ( COMPLEX16Array * );
void LALZCreateArray ( LALStatus *, COMPLEX16Array **, UINT4Vector * );
void LALZResizeArray ( LALStatus *, COMPLEX16Array **, UINT4Vector * );
void LALZDestroyArray ( LALStatus *, COMPLEX16Array ** );
/*@}*/


/*@}*/
/* ---------- end: ArrayFactories_c ---------- */


/**
 * \defgroup VectorFactories_c Module VectorFactories.c
 * \ingroup AVFactories_h
 *
 * \brief Create/destroy \<datatype\>%Vector objects.
 *
 * \heading{Description}
 *
 * The \c CreateVector family of functions create a \<datatype\>%Vector of the appropriate dimensions.
 *
 * The \c ResizeVector family of functions changes the amount of storage allocated by the \c CreateVector functions.
 *
 * The \c DestroyVector family of functions return the storage allocated by the \c CreateVector functions to the system.
 *
 */
/*@{*/

/** \name CHAR vector prototypes */
/*@{*/
CHARVector * XLALCreateCHARVector ( UINT4 length );
CHARVector * XLALResizeCHARVector ( CHARVector * vector, UINT4 length );
void XLALDestroyCHARVector ( CHARVector * vector );
void LALCHARCreateVector ( LALStatus *, CHARVector **, UINT4 );
void LALCHARResizeVector ( LALStatus *, CHARVector **, UINT4 );
void LALCHARDestroyVector ( LALStatus *, CHARVector ** );
/*@}*/

/** \name INT2 vector prototypes */
/*@{*/
INT2Vector * XLALCreateINT2Vector ( UINT4 length );
INT2Vector * XLALResizeINT2Vector ( INT2Vector * vector, UINT4 length );
void XLALDestroyINT2Vector ( INT2Vector * vector );
void LALI2CreateVector ( LALStatus *, INT2Vector **, UINT4 );
void LALI2ResizeVector ( LALStatus *, INT2Vector **, UINT4 );
void LALI2DestroyVector ( LALStatus *, INT2Vector ** );
/*@}*/

/** \name INT4 vector prototypes */
/*@{*/
INT4Vector * XLALCreateINT4Vector ( UINT4 length );
INT4Vector * XLALResizeINT4Vector ( INT4Vector * vector, UINT4 length );
void XLALDestroyINT4Vector ( INT4Vector * vector );
void LALI4CreateVector ( LALStatus *, INT4Vector **, UINT4 );
void LALI4ResizeVector ( LALStatus *, INT4Vector **, UINT4 );
void LALI4DestroyVector ( LALStatus *, INT4Vector ** );
/*@}*/


/** \name INT8 vector prototypes */
/*@{*/
INT8Vector * XLALCreateINT8Vector ( UINT4 length );
INT8Vector * XLALResizeINT8Vector ( INT8Vector * vector, UINT4 length );
void XLALDestroyINT8Vector ( INT8Vector * vector );
void LALI8CreateVector ( LALStatus *, INT8Vector **, UINT4 );
void LALI8ResizeVector ( LALStatus *, INT8Vector **, UINT4 );
void LALI8DestroyVector ( LALStatus *, INT8Vector ** );
/*@}*/

/** \name UINT2 vector prototypes */
/*@{*/
UINT2Vector * XLALCreateUINT2Vector ( UINT4 length );
UINT2Vector * XLALResizeUINT2Vector ( UINT2Vector * vector, UINT4 length );
void XLALDestroyUINT2Vector ( UINT2Vector * vector );
void LALU2CreateVector ( LALStatus *, UINT2Vector **, UINT4 );
void LALU2ResizeVector ( LALStatus *, UINT2Vector **, UINT4 );
void LALU2DestroyVector ( LALStatus *, UINT2Vector ** );
/*@}*/

/** \name UINT4 vector prototypes */
/*@{*/
UINT4Vector * XLALCreateUINT4Vector ( UINT4 length );
UINT4Vector * XLALResizeUINT4Vector ( UINT4Vector * vector, UINT4 length );
void XLALDestroyUINT4Vector ( UINT4Vector * vector );
void LALU4CreateVector ( LALStatus *, UINT4Vector **, UINT4 );
void LALU4ResizeVector ( LALStatus *, UINT4Vector **, UINT4 );
void LALU4DestroyVector ( LALStatus *, UINT4Vector ** );
/*@}*/

/** \name UINT8 vector prototypes */
/*@{*/
UINT8Vector * XLALCreateUINT8Vector ( UINT4 length );
UINT8Vector * XLALResizeUINT8Vector ( UINT8Vector * vector, UINT4 length );
void XLALDestroyUINT8Vector ( UINT8Vector * vector );
void LALU8CreateVector ( LALStatus *, UINT8Vector **, UINT4 );
void LALU8ResizeVector ( LALStatus *, UINT8Vector **, UINT4 );
void LALU8DestroyVector ( LALStatus *, UINT8Vector ** );
/*@}*/

/** \name REAL4 vector prototypes */
/*@{*/
REAL4Vector * XLALCreateREAL4Vector ( UINT4 length );
REAL4Vector * XLALResizeREAL4Vector ( REAL4Vector * vector, UINT4 length );
void XLALDestroyREAL4Vector ( REAL4Vector * vector );
void LALSCreateVector ( LALStatus *, REAL4Vector **, UINT4 );
void LALSResizeVector ( LALStatus *, REAL4Vector **, UINT4 );
void LALSDestroyVector ( LALStatus *, REAL4Vector ** );
/*@}*/

/** \name REAL4 vector prototypes (default name) */
/*@{*/
REAL4Vector * XLALCreateVector ( UINT4 length );
REAL4Vector * XLALResizeVector ( REAL4Vector * vector, UINT4 length );
void XLALDestroyVector ( REAL4Vector * vector );
void LALCreateVector ( LALStatus *, REAL4Vector **, UINT4 );
void LALResizeVector ( LALStatus *, REAL4Vector **, UINT4 );
void LALDestroyVector ( LALStatus *, REAL4Vector ** );
/*@}*/

/** \name REAL8 vector prototypes */
/*@{*/
REAL8Vector * XLALCreateREAL8Vector ( UINT4 length );
REAL8Vector * XLALResizeREAL8Vector ( REAL8Vector * vector, UINT4 length );
void XLALDestroyREAL8Vector ( REAL8Vector * vector );
void LALDCreateVector ( LALStatus *, REAL8Vector **, UINT4 );
void LALDResizeVector ( LALStatus *, REAL8Vector **, UINT4 );
void LALDDestroyVector ( LALStatus *, REAL8Vector ** );
/*@}*/

/** \name COMPLEX8 vector prototypes */
/*@{*/
COMPLEX8Vector * XLALCreateCOMPLEX8Vector ( UINT4 length );
COMPLEX8Vector * XLALResizeCOMPLEX8Vector ( COMPLEX8Vector * vector, UINT4 length );
void XLALDestroyCOMPLEX8Vector ( COMPLEX8Vector * vector );
void LALCCreateVector ( LALStatus *, COMPLEX8Vector **, UINT4 );
void LALCResizeVector ( LALStatus *, COMPLEX8Vector **, UINT4 );
void LALCDestroyVector ( LALStatus *, COMPLEX8Vector ** );
/*@}*/

/** \name COMPLEX16 vector prototypes */
/*@{*/
COMPLEX16Vector * XLALCreateCOMPLEX16Vector ( UINT4 length );
COMPLEX16Vector * XLALResizeCOMPLEX16Vector ( COMPLEX16Vector * vector, UINT4 length );
void XLALDestroyCOMPLEX16Vector ( COMPLEX16Vector * vector );
void LALZCreateVector ( LALStatus *, COMPLEX16Vector **, UINT4 );
void LALZResizeVector ( LALStatus *, COMPLEX16Vector **, UINT4 );
void LALZDestroyVector ( LALStatus *, COMPLEX16Vector ** );
/*@}*/

/*@}*/

/* ---------- end: VectorFactories_c ---------- */

#ifdef  __cplusplus
}
#endif

#endif /* _AVFACTORIES_H */
