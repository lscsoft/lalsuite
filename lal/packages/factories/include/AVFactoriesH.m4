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

NRCSID (AVFACTORIESH, "$Id$");

/**
\defgroup AVFactories_h AVFactories_h
\ingroup factories

\brief Provides prototype and status code information for use of CreateVector,
CreateArray, ResizeVector, ResizeArray, DestroyVector and DestroyArray

\heading{Synopsis}
\code
#include <lal/AVFactories.h>
\endcode

\section secXLALfcts XLAL Functions

\code
REAL4Vector * XLALCreateVector(UINT4 length);
REAL4Vector * XLALResizeVector(REAL4Vector *vector, UINT4 length);
void XLALDestroyVector(REAL4Vector *vector, UINT4 length);

<vectype> * XLALCreate<vectype>(UINT4 length );
<vectype> * XLALResize<vectype>(<vectype> *vector, UINT4 length );
void XLALDestroy<vectype>(<vectype> *vector);

REAL4Array * XLALCreateArrayL(UINT4 ndim, ...);
REAL4Array * XLALCreateArrayV(UINT4 ndim, UINT4 *dims);
REAL4Array * XLALCreateArray(UINT4Vector *dimLength);
REAL4Array * XLALResizeArrayL(REAL4Array *array, UINT4 ndim, ...);
REAL4Array * XLALResizeArrayV(REAL4Array *array, UINT4 ndim, UINT4 *dims);
REAL4Array * XLALResizeArray(REAL4Array *array, UINT4Vector *dimLength);
void XLALDestroyArray(REAL4Array *array);

<arrtype> * XLALCreate<arrtype>L(UINT4 ndim, ...);
<arrtype> * XLALCreate<arrtype>V(UINT4 ndim, UINT4 *dims);
<arrtype> * XLALCreate<arrtype>(UINT4Vector *dimLength);
<arrtype> * XLALResize<arrtype>L(<arrtype> *array, UINT4 ndim, ...);
<arrtype> * XLALResize<arrtype>V(<arrtype> *array, UINT4 ndim, UINT4 *dims);
<arrtype> * XLALResize<arrtype>(<arrtype> *array, UINT4Vector *dimLength);
void XLALDestroy<arrtype>(<arrtype> *array);
\endcode

Here <tt>\<vectype\></tt> is one of
\c COMPLEX16Vector,
\c COMPLEX8Vector,
\c REAL8Vector,
\c REAL4Vector,
\c INT8Vector,
\c INT4Vector,
\c INT2Vector,
\c UINT8Vector,
\c UINT4Vector,
\c UINT2Vector, or
\c CHARVector,
and <tt>\<arrtype\></tt> is one of
\c COMPLEX16Array,
\c COMPLEX8Array,
\c REAL8Array,
\c REAL4Array,
\c INT8Array,
\c INT4Array,
\c INT2Array,
\c UINT8Array,
\c UINT4Array, or
\c UINT2Array.

\subsection ss_AVF_desc Description

The <tt>XLALCreate\<type\>%Vector</tt> functions create vectors of the specified
\c length number of objects of type <tt>\<type\></tt>.  The function
\c XLALCreateVector() is the same as \c XLALCreateREAL4Vector().

The <tt>XLALDestroy\<type\>%Vector</tt> functions deallocate the memory allocation
pointed to by \c vector including its contents.  The function
\c XLALDestroyVector() is the same as \c XLALDestroyREAL4Vector().

The <tt>XLALResize\<type\>%Vector</tt> functions resize the supplied vector
\c vector to the new size \c length.  If \c vector is \c NULL
then this is equivalent to <tt>XLALCreate\<type\>%Vector</tt>.  If \c length
is zero then this is equivalent to <tt>XLALDestroy\<type\>%Vector</tt> and the
routine returns \c NULL.  Otherwise, the amount of data in the vector
is realloced using \c LALRealloc().  The function
\c XLALResizeVector() is the same as \c XLALResizeREAL4Vector().

The <tt>XLALCreate\<type\>Array</tt>
<tt>XLALCreate\<type\>ArrayL</tt>
<tt>XLALCreate\<type\>ArrayV</tt>
all create an object of type <tt>\<type\>Array</tt>.  They differ in the way
that the dimensions of the array are specified.  The function
<tt>XLALCreate\<type\>Array</tt> allocates an array with dimensions specified
by the \c UINT4Vector \c dimLength which is a vector of dimension
lengths for the array.  The function
<tt>XLALCreate\<type\>ArrayV</tt> provides these dimensions with two arguments:
\c ndim is the number of dimensions and \c dims is an array of
\c UINT4 values for the dimensions.  The function
<tt>XLALCreate\<type\>ArrayL</tt> also specifies the dimensions as arguments.
Here, the first argument, \c ndim, is the number of dimensions, and
this is followed by \c ndim arguments that provide the dimensions.
Note that for this function, a maximum of 16 dimensions can be provided
(that is, \c ndim cannot be more than 16 and there cannot be more than
16 arguments after the first).
The \c XLALCreateArray()
\c XLALCreateArrayL()
\c XLALCreateArrayV()
functions are equivalent to the
\c XLALCreateREAL4Array()
\c XLALCreateREAL4ArrayL()
\c XLALCreateREAL4ArrayV()
functions respectively.

The <tt>XLALDestroy\<type\>Array</tt> functions deallocate the memory allocation
pointed to by \c array including its contents.  The function
\c XLALDestroyArray() is the same as \c XLALDestroyREAL4Array().

The <tt>XLALResize\<type\>Array</tt>
<tt>XLALResize\<type\>ArrayL</tt>
<tt>XLALResize\<type\>ArrayV</tt>
functions resize the provided array \c array.  The arguments after the
first are interpreted in the same way as for the
<tt>XLALCreate\<type\>Array</tt>
<tt>XLALCreate\<type\>ArrayL</tt>
<tt>XLALCreate\<type\>ArrayV</tt>
functions.  If \c array is \c NULL, the resize functions are equivalent
to the corresponding create function.  If \c ndim is zero for
<tt>XLALResize\<type\>ArrayL</tt> or
<tt>XLALResize\<type\>ArrayV</tt>, or if \c dimLength is \c NULL for
<tt>XLALResize\<type\>Array</tt>, then these functions are equivalent to
<tt>XLALDestroy\<type\>Array</tt>.
The \c XLALResizeArray()
\c XLALResizeArrayL()
\c XLALResizeArrayV()
functions are equivalent to the
\c XLALResizeREAL4Array()
\c XLALResizeREAL4ArrayL()
\c XLALResizeREAL4ArrayV()
functions respectively.

\subsection ss_AVF_return Return Values

If successful, the create and resize functions return a pointer to the same
data that was passed to the function.  The resize functions will return
a \c NULL pointer if the size of the new object was zero.  Upon failure
these routines will return \c NULL and will set \c xlalErrno to
one of these values: \c #XLAL_ENOMEM if a memory allocation failed,
\c #XLAL_EBADLEN if an invalid length was provided (for example, a
zero-size allocation with a create function), \c XLAL_EINVAL if an
invalid argument is provided (for example, if the pointer to an
array of dimensions is \c NULL).

The destroy function does not return any value.  If the function is passed a
\c NULL pointer, it will set \c xlalErrno to \c #XLAL_EFAULT.
If the function is passed an object that appears corrupted (e.g., a
vector with zero length or will a \c NULL data pointer) it will set
\c xlalErrno to \c #XLAL_EINVAL.

*/
/*@{ */
/** \name Error Codes *//*@{*/
#define AVFACTORIESH_ELENGTH 1
#define AVFACTORIESH_EVPTR   2
#define AVFACTORIESH_EUPTR   4
#define AVFACTORIESH_EDPTR   8
#define AVFACTORIESH_EMALLOC 16
#define AVFACTORIESH_MSGELENGTH  "Illegal length."
#define AVFACTORIESH_MSGEVPTR    "Null vector/array handle."
#define AVFACTORIESH_MSGEUPTR    "Non-null vector/array pointer."
#define AVFACTORIESH_MSGEDPTR    "Null vector/array data."
#define AVFACTORIESH_MSGEMALLOC  "Malloc failure."
/*@}*/
/*@}*/

dnl There is no CHARArray type
define(`TYPECODE',`CHAR')
include(`VFactoriesBaseH.m4')

define(`TYPECODE',`I2')
include(`VFactoriesBaseH.m4')
include(`AFactoriesBaseH.m4')

define(`TYPECODE',`I4')
include(`VFactoriesBaseH.m4')
include(`AFactoriesBaseH.m4')

define(`TYPECODE',`I8')
include(`VFactoriesBaseH.m4')
include(`AFactoriesBaseH.m4')

define(`TYPECODE',`U2')
include(`VFactoriesBaseH.m4')
include(`AFactoriesBaseH.m4')

define(`TYPECODE',`U4')
include(`VFactoriesBaseH.m4')
include(`AFactoriesBaseH.m4')

define(`TYPECODE',`U8')
include(`VFactoriesBaseH.m4')
include(`AFactoriesBaseH.m4')

define(`TYPECODE',`S')
include(`VFactoriesBaseH.m4')
include(`AFactoriesBaseH.m4')

define(`TYPECODE',`D')
include(`VFactoriesBaseH.m4')
include(`AFactoriesBaseH.m4')

define(`TYPECODE',`C')
include(`VFactoriesBaseH.m4')
include(`AFactoriesBaseH.m4')

define(`TYPECODE',`Z')
include(`VFactoriesBaseH.m4')
include(`AFactoriesBaseH.m4')

define(`TYPECODE',`')
include(`VFactoriesBaseH.m4')
include(`AFactoriesBaseH.m4')

#ifdef  __cplusplus
}
#endif

#endif /* _AVFACTORIES_H */
