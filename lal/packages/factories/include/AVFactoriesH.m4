/*----------------------------------------------------------------------- 

File Name: AVFactories.h

<lalVerbatim file="AVFactoriesHV">
Revision: $Id$
</lalVerbatim>

-------------------------------------------------------------------------*/

/* <lalLaTeX>

\section{Header \texttt{AVFactories.h}}
\label{s:AVFactories.h}

Provides prototype and status code information for use of CreateVector,
CreateArray, ResizeVector, ResizeArray, DestroyVector and DestroyArray

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/AVFactories.h>
\end{verbatim}

</lalLaTeX> */

#ifndef _AVFACTORIES_H
#define _AVFACTORIES_H

#include <lal/LALDatatypes.h>
#include <stdarg.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (AVFACTORIESH, "$Id$");

/* <lalLaTeX>

\subsection*{Error conditions}
\input{AVFactoriesHErrTab}

</lalLaTeX> */

/*
<lalErrTable file="AVFactoriesHErrTab">
*/
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
/*
</lalErrTable>
*/

/* Function prototypes. */
/* <lalLaTeX>
\newpage\input{VectorFactoriesC}
\newpage\input{ArrayFactoriesC}

\newpage
\subsection{XLAL Functions}

\subsubsection*{Synopsis}
\begin{verbatim}
#include <lal/AVFactories.h>

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
\end{verbatim}
\idx{XLALCreateVector}
\idx{XLALResizeVector}
\idx{XLALDestroyVector}
\idx{XLALCreateArrayL}
\idx{XLALCreateArrayV}
\idx{XLALCreateArray}
\idx{XLALResizeArrayL}
\idx{XLALResizeArrayV}
\idx{XLALResizeArray}
\idx{XLALDestroyArray}
\idx{XLALCreate<type>Vector}
\idx{XLALResize<type>Vector}
\idx{XLALDestroy<type>Vector}
\idx{XLALCreate<type>ArrayL}
\idx{XLALCreate<type>ArrayV}
\idx{XLALCreate<type>Array}
\idx{XLALResize<type>ArrayL}
\idx{XLALResize<type>ArrayV}
\idx{XLALResize<type>Array}
\idx{XLALDestroy<type>Array}

Here \verb+<vectype>+ is one of
\verb+COMPLEX16Vector+,
\verb+COMPLEX8Vector+,
\verb+REAL8Vector+,
\verb+REAL4Vector+,
\verb+INT8Vector+,
\verb+INT4Vector+,
\verb+INT2Vector+,
\verb+UINT8Vector+,
\verb+UINT4Vector+,
\verb+UINT2Vector+, or
\verb+CHARVector+,
and \verb+<arrtype>+ is one of
\verb+COMPLEX16Array+,
\verb+COMPLEX8Array+,
\verb+REAL8Array+,
\verb+REAL4Array+,
\verb+INT8Array+,
\verb+INT4Array+,
\verb+INT2Array+,
\verb+UINT8Array+,
\verb+UINT4Array+, or
\verb+UINT2Array+.

\subsubsection*{Description}

The \verb+XLALCreate<type>Vector+ functions create vectors of the specified
\verb+length+ number of objects of type \verb+<type>+.  The function
\verb+XLALCreateVector+ is the same as \verb+XLALCreateREAL4Vector+.

The \verb+XLALDestroy<type>Vector+ functions deallocate the memory allocation
pointed to by \verb+vector+ including its contents.  The function
\verb+XLALDestroyVector+ is the same as \verb+XLALDestroyREAL4Vector+.

The \verb+XLALResize<type>Vector+ functions resize the supplied vector
\verb+vector+ to the new size \verb+length+.  If \verb+vector+ is \verb+NULL+
then this is equivalent to \verb+XLALCreate<type>Vector+.  If \verb+length+
is zero then this is equivalent to \verb+XLALDestroy<type>Vector+ and the
routine returns \verb+NULL+.  Otherwise, the amount of data in the vector
is realloced using \verb+LALRealloc+.  The function
\verb+XLALResizeVector+ is the same as \verb+XLALResizeREAL4Vector+.

The \verb+XLALCreate<type>Array+
\verb+XLALCreate<type>ArrayL+
\verb+XLALCreate<type>ArrayV+
all create an object of type \verb+<type>Array+.  They differ in the way
that the dimensions of the array are specified.  The function
\verb+XLALCreate<type>Array+ allocates an array with dimensions specified
by the \verb+UINT4Vector+ \verb+dimLength+ which is a vector of dimension
lengths for the array.  The function
\verb+XLALCreate<type>ArrayV+ provides these dimensions with two arguments:
\verb+ndim+ is the number of dimensions and \verb+dims+ is an array of
\verb+UINT4+ values for the dimensions.  The function
\verb+XLALCreate<type>ArrayL+ also specifies the dimensions as arguments.
Here, the first argument, \verb+ndim+, is the number of dimensions, and
this is followed by \verb+ndim+ arguments that provide the dimensions.
Note that for this function, a maximum of 16 dimensions can be provided
(that is, \verb+ndim+ cannot be more than 16 and there cannot be more than
16 arguments after the first).
The \verb+XLALCreateArray+
\verb+XLALCreateArrayL+
\verb+XLALCreateArrayV+
functions are equivalent to the
\verb+XLALCreateREAL4Array+
\verb+XLALCreateREAL4ArrayL+
\verb+XLALCreateREAL4ArrayV+
functions respectively.

The \verb+XLALDestroy<type>Array+ functions deallocate the memory allocation
pointed to by \verb+array+ including its contents.  The function
\verb+XLALDestroyArray+ is the same as \verb+XLALDestroyREAL4Array+.

The \verb+XLALResize<type>Array+
\verb+XLALResize<type>ArrayL+
\verb+XLALResize<type>ArrayV+
functions resize the provided array \verb+array+.  The arguments after the
first are interpreted in the same way as for the 
\verb+XLALCreate<type>Array+
\verb+XLALCreate<type>ArrayL+
\verb+XLALCreate<type>ArrayV+
functions.  If \verb+array+ is \verb+NULL+, the resize functions are equivalent
to the corresponding create function.  If \verb+ndim+ is zero for 
\verb+XLALResize<type>ArrayL+ or
\verb+XLALResize<type>ArrayV+, or if \verb+dimLength+ is \verb+NULL+ for
\verb+XLALResize<type>Array+, then these functions are equivalent to
\verb+XLALDestroy<type>Array+.
The \verb+XLALResizeArray+
\verb+XLALResizeArrayL+
\verb+XLALResizeArrayV+
functions are equivalent to the
\verb+XLALResizeREAL4Array+
\verb+XLALResizeREAL4ArrayL+
\verb+XLALResizeREAL4ArrayV+
functions respectively.

\subsubsection*{Return Values}

If successful, the create and resize functions return a pointer to the same
data that was passed to the function.  The resize functions will return
a \verb+NULL+ pointer if the size of the new object was zero.  Upon failure
these routines will return \verb+NULL+ and will set \verb+xlalErrno+ to
one of these values: \verb+XLAL_ENOMEM+ if a memory allocation failed,
\verb+XLAL_EBADLEN+ if an invalid length was provided (for example, a
zero-size allocation with a create function), \verb+XLAL_EINVAL+ if an
invalid argument is provided (for example, if the pointer to an
array of dimensions is \verb+NULL+).

The destroy function does not return any value.  If the function is passed a
\verb+NULL+ pointer, it will set \verb+xlalErrno+ to \verb+XLAL_EFAULT+.
If the function is passed an object that appears corrupted (e.g., a
vector with zero length or will a \verb+NULL+ data pointer) it will set
\verb+xlalErrno+ to \verb+XLAL_EINVAL+.

</lalLaTeX> */

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

/* Test program. */

/* <lalLaTeX>
\newpage\input{VectorFactoriesTestC}
\newpage\input{ArrayFactoriesTestC}
</lalLaTeX> */

#ifdef  __cplusplus
}
#endif

#endif /* _AVFACTORIES_H */
