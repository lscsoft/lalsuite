/*
 *
 * Copyright (C) 2007  Kipp Cannon
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


#ifndef _SEQUENCE_H
#define _SEQUENCE_H


#include <stddef.h>
#include <lal/LALDatatypes.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * \addtogroup SequenceManipulation
 * \author Kipp Cannon <kipp@gravity.phys.uwm.edu>
 *
 * \brief This is a suite of functions for creating, destroying, and manipulating LAL
 * sequences.  For example XLALCreateREAL4Sequence() is available for
 * creating sequences of \c REAL4 data.
 */
/*@{*/


/**
 * \name Creation Functions
 * \heading{Synopsis}
 *
 * \code
 * #include <lal/Sequence.h>
 *
 * XLALCreate<sequencetype>()
 * LALCreate<sequencetype>()
 * \endcode
 *
 * \heading{Description}
 *
 * These functions create LAL sequences.  The parameter \c length
 * specifies the length of the desired sequence.  The return value is a
 * pointer to the newly created sequence or \c NULL on failure.
 */
/*@{*/
COMPLEX8Sequence *XLALCreateCOMPLEX8Sequence ( 	size_t length );
COMPLEX16Sequence *XLALCreateCOMPLEX16Sequence ( size_t length );
REAL4Sequence *XLALCreateREAL4Sequence ( size_t length );
REAL8Sequence *XLALCreateREAL8Sequence ( size_t length );
INT2Sequence *XLALCreateINT2Sequence ( size_t length );
INT4Sequence *XLALCreateINT4Sequence ( size_t length );
INT8Sequence *XLALCreateINT8Sequence ( 	size_t length );
UINT2Sequence *XLALCreateUINT2Sequence ( size_t length );
UINT4Sequence *XLALCreateUINT4Sequence ( size_t length );
UINT8Sequence *XLALCreateUINT8Sequence ( size_t length );
/*@}*/

/**
 * \name Destruction Functions
 *
 * \heading{Synopsis}
 *
 * \code
 * #include <lal/Sequence.h>
 *
 * XLALDestroy<sequencetype>()
 * \endcode
 *
 * \heading{Description}
 *
 * These functions free all memory associated with a LAL sequence.  It is safe
 * to pass \c NULL to these functions.
 *
 */
/*@{*/
void XLALDestroyCOMPLEX8Sequence ( COMPLEX8Sequence *sequence );
void XLALDestroyCOMPLEX16Sequence ( COMPLEX16Sequence *sequence );
void XLALDestroyREAL4Sequence ( REAL4Sequence *sequence );
void XLALDestroyREAL8Sequence ( REAL8Sequence *sequence );
void XLALDestroyINT2Sequence ( INT2Sequence *sequence );
void XLALDestroyINT4Sequence ( INT4Sequence *sequence );
void XLALDestroyINT8Sequence ( INT8Sequence *sequence );
void XLALDestroyUINT2Sequence ( UINT2Sequence *sequence );
void XLALDestroyUINT4Sequence ( UINT4Sequence *sequence );
void XLALDestroyUINT8Sequence ( UINT8Sequence *sequence );
/*@}*/

/**
 * \name Cutting Functions
 * \heading{Synopsis}
 *
 * \code
 * #include <lal/Sequence.h>
 *
 * XLALCut<sequencetype>()
 * XLALCopy<sequencetype>()
 * \endcode
 *
 * \heading{Description}
 *
 * These functions create a new sequence by extracting a section of an
 * existing sequence.
 *
 */
/*@{*/
COMPLEX8Sequence *XLALCutCOMPLEX8Sequence ( COMPLEX8Sequence *sequence, size_t first, size_t length );
COMPLEX16Sequence *XLALCutCOMPLEX16Sequence ( COMPLEX16Sequence *sequence, size_t first, size_t length );
REAL4Sequence *XLALCutREAL4Sequence ( REAL4Sequence *sequence, size_t first, size_t length );
REAL8Sequence *XLALCutREAL8Sequence ( REAL8Sequence *sequence, size_t first, size_t length );
INT2Sequence *XLALCutINT2Sequence ( INT2Sequence *sequence, size_t first, size_t length );
INT4Sequence *XLALCutINT4Sequence ( INT4Sequence *sequence, size_t first, size_t length );
INT8Sequence *XLALCutINT8Sequence ( INT8Sequence *sequence, size_t first, size_t length );
UINT2Sequence *XLALCutUINT2Sequence ( UINT2Sequence *sequence, size_t first, size_t length );
UINT4Sequence *XLALCutUINT4Sequence ( UINT4Sequence *sequence, size_t first, size_t length );
UINT8Sequence *XLALCutUINT8Sequence ( UINT8Sequence *sequence, size_t first, size_t length );

COMPLEX8Sequence *XLALCopyCOMPLEX8Sequence ( COMPLEX8Sequence *sequence );
COMPLEX16Sequence *XLALCopyCOMPLEX16Sequence ( COMPLEX16Sequence *sequence );
REAL4Sequence *XLALCopyREAL4Sequence ( REAL4Sequence *sequence );
REAL8Sequence *XLALCopyREAL8Sequence ( REAL8Sequence *sequence );
INT2Sequence *XLALCopyINT2Sequence ( INT2Sequence *sequence );
INT4Sequence *XLALCopyINT4Sequence ( INT4Sequence *sequence );
INT8Sequence *XLALCopyINT8Sequence ( INT8Sequence *sequence );
UINT2Sequence *XLALCopyUINT2Sequence ( UINT2Sequence *sequence );
UINT4Sequence *XLALCopyUINT4Sequence ( UINT4Sequence *sequence );
UINT8Sequence *XLALCopyUINT8Sequence ( UINT8Sequence *sequence );
/*@}*/

/**
 * \name Shifting Functions
 * \heading{Synopsis}
 *
 * \code
 * #include <lal/Sequence.h>
 *
 * XLALShift<sequencetype>()
 * \endcode
 *
 * \heading{Description}
 *
 * These functions shift the samples in a sequence, with zeros being placed in
 * the space that is freed.
 */
/*@{*/
void XLALShiftCOMPLEX8Sequence ( COMPLEX8Sequence *sequence, int count );
void XLALShiftCOMPLEX16Sequence ( COMPLEX16Sequence *sequence, int count );
void XLALShiftREAL4Sequence ( REAL4Sequence *sequence, int count );
void XLALShiftREAL8Sequence ( REAL8Sequence *sequence, int count );
void XLALShiftINT2Sequence ( INT2Sequence *sequence, int count );
void XLALShiftINT4Sequence ( INT4Sequence *sequence, int count );
void XLALShiftINT8Sequence ( INT8Sequence *sequence, int count );
void XLALShiftUINT2Sequence ( UINT2Sequence *sequence, int count );
void XLALShiftUINT4Sequence ( UINT4Sequence *sequence, int count );
void XLALShiftUINT8Sequence ( UINT8Sequence *sequence, int count );
/*@}*/

/**
 * \name Resizing Functions
 *
 * \heading{Synopsis}
 *
 * \code
 * #include <lal/Sequence.h>
 *
 * XLALResize<sequencetype>()
 * XLALShrink<sequencetype>()
 * \endcode
 *
 * \heading{Description}
 *
 * The resize functions alter the size of an existing sequence.  The sequence
 * is adjusted to have the specified length, and that part of the original
 * sequence starting at sample first is used to fill the new sequence.  If
 * first is negative, then the start of the new sequence is padded by that
 * many samples.  If part of the new sequence does not correspond to some part
 * of the original sequence, then those samples are set to 0.
 *
 * The shrink functions, originally, could only handle the special case in
 * which the new sequence is wholly contained in the original sequence.  Now
 * the shrink functions are wrappers for the resize functions and are only
 * retained for backwards compatibility.
 *
 */
/*@{*/
COMPLEX8Sequence *XLALResizeCOMPLEX8Sequence ( COMPLEX8Sequence *sequence, int first, size_t length );
COMPLEX16Sequence *XLALResizeCOMPLEX16Sequence ( COMPLEX16Sequence *sequence, int first, size_t length );
REAL4Sequence *XLALResizeREAL4Sequence ( REAL4Sequence *sequence, int first, size_t length );
REAL8Sequence *XLALResizeREAL8Sequence ( REAL8Sequence *sequence, int first, size_t length );
INT2Sequence *XLALResizeINT2Sequence ( INT2Sequence *sequence, int first, size_t length );
INT4Sequence *XLALResizeINT4Sequence ( INT4Sequence *sequence, int first, size_t length );
INT8Sequence *XLALResizeINT8Sequence ( INT8Sequence *sequence, int first, size_t length );
UINT2Sequence *XLALResizeUINT2Sequence ( UINT2Sequence *sequence, int first, size_t length );
UINT4Sequence *XLALResizeUINT4Sequence ( UINT4Sequence *sequence, int first, size_t length );
UINT8Sequence *XLALResizeUINT8Sequence ( UINT8Sequence *sequence, int first, size_t length );

COMPLEX8Sequence *XLALShrinkCOMPLEX8Sequence ( COMPLEX8Sequence *sequence, size_t first, size_t length );
COMPLEX16Sequence *XLALShrinkCOMPLEX16Sequence ( COMPLEX16Sequence *sequence, size_t first, size_t length );
REAL4Sequence *XLALShrinkREAL4Sequence ( REAL4Sequence *sequence, size_t first, size_t length );
REAL8Sequence *XLALShrinkREAL8Sequence ( REAL8Sequence *sequence, size_t first, size_t length );
INT2Sequence *XLALShrinkINT2Sequence ( INT2Sequence *sequence, size_t first, size_t length );
INT4Sequence *XLALShrinkINT4Sequence ( INT4Sequence *sequence, size_t first, size_t length );
INT8Sequence *XLALShrinkINT8Sequence ( INT8Sequence *sequence, size_t first, size_t length );
UINT2Sequence *XLALShrinkUINT2Sequence ( UINT2Sequence *sequence, size_t first, size_t length );
UINT4Sequence *XLALShrinkUINT4Sequence ( UINT4Sequence *sequence, size_t first, size_t length );
UINT8Sequence *XLALShrinkUINT8Sequence ( UINT8Sequence *sequence, size_t first, size_t length );
/*@}*/

/**
 * \name Summing Functions
 * \heading{Synopsis}
 *
 * \code
 * #include <lal/Sequence.h>
 *
 * XLAL<datatype>Sum()
 * XLAL<datatype>SumSquares()
 * XLAL<sequencetype>Sum>()
 * XLAL<sequencetype>SumSquares()
 * \endcode
 *
 * \heading{Description}
 *
 * The \c XLAL\<datatype\>Sum() and
 * \c XLAL\<datatype\>SumSquares() functions sum the
 * elements and squares of the elements, respectively, in an array.
 *
 * The \c XLAL\<sequencetype\>Sum() and
 * \c XLAL\<sequencetype\>SumSquares() functions sum the
 * elements and the squares of the elements, respectively in a sequence.
 * Bounds checking is performed.
 *
 * In all cases, the return value is the sum, and these functions cannot fail.
 * In the case of the sequence-related functions, if the sum extends beyond
 * the bounds of the sequence, then the missing values are assumed to be 0.
 *
 * \heading{Bugs}
 *
 * Because the LAL library must conform to the C89 specification, aggregate
 * data types cannot be returned from functions so the COMPLEX8 and COMPLEX16
 * versions of the sum functions (not sum-of-squares functions) are commented
 * out at this time.
 */
/*@{*/
COMPLEX8 XLALCOMPLEX8Sum ( const COMPLEX8 *data, size_t first, size_t count );
REAL4 XLALCOMPLEX8SumSquares ( const COMPLEX8 *data, size_t first, size_t count );
COMPLEX8 XLALCOMPLEX8SequenceSum ( const COMPLEX8Sequence *sequence, size_t first, size_t count );
REAL4 XLALCOMPLEX8SequenceSumSquares ( const COMPLEX8Sequence *sequence, size_t first, size_t count );
COMPLEX16 XLALCOMPLEX16Sum ( const COMPLEX16 *data, size_t first, size_t count );
REAL8 XLALCOMPLEX16SumSquares ( const COMPLEX16 *data, size_t first, size_t count );
COMPLEX16 XLALCOMPLEX16SequenceSum ( const COMPLEX16Sequence *sequence, size_t first, size_t count );
REAL8 XLALCOMPLEX16SequenceSumSquares ( const COMPLEX16Sequence *sequence, size_t first, size_t count );
REAL4 XLALREAL4Sum ( const REAL4 *data, size_t first, size_t count );
REAL4 XLALREAL4SumSquares ( const REAL4 *data, size_t first, size_t count );
REAL4 XLALREAL4SequenceSum ( const REAL4Sequence *sequence, size_t first, size_t count );
REAL4 XLALREAL4SequenceSumSquares ( const REAL4Sequence *sequence, size_t first, size_t count );
REAL8 XLALREAL8Sum ( const REAL8 *data, size_t first, size_t count );
REAL8 XLALREAL8SumSquares ( const REAL8 *data, size_t first, size_t count );
REAL8 XLALREAL8SequenceSum ( const REAL8Sequence *sequence, size_t first, size_t count );
REAL8 XLALREAL8SequenceSumSquares ( const REAL8Sequence *sequence, size_t first, size_t count );
INT2 XLALINT2Sum ( const INT2 *data, size_t first, size_t count );
UINT2 XLALINT2SumSquares ( const INT2 *data, size_t first, size_t count );
INT2 XLALINT2SequenceSum ( const INT2Sequence *sequence, size_t first, size_t count );
UINT2 XLALINT2SequenceSumSquares ( const INT2Sequence *sequence, size_t first, size_t count );
INT4 XLALINT4Sum ( const INT4 *data, size_t first, size_t count );
UINT4 XLALINT4SumSquares ( const INT4 *data, size_t first, size_t count );
INT4 XLALINT4SequenceSum ( const INT4Sequence *sequence, size_t first, size_t count );
UINT4 XLALINT4SequenceSumSquares ( const INT4Sequence *sequence, size_t first, size_t count );
INT8 XLALINT8Sum ( const INT8 *data, size_t first, size_t count );
UINT8 XLALINT8SumSquares ( const INT8 *data, size_t first, size_t count );
INT8 XLALINT8SequenceSum ( const INT8Sequence *sequence, size_t first, size_t count );
UINT8 XLALINT8SequenceSumSquares ( const INT8Sequence *sequence, size_t first, size_t count );
UINT2 XLALUINT2Sum ( const UINT2 *data, size_t first, size_t count );
UINT2 XLALUINT2SumSquares ( const UINT2 *data, size_t first, size_t count );
UINT2 XLALUINT2SequenceSum ( const UINT2Sequence *sequence, size_t first, size_t count );
UINT2 XLALUINT2SequenceSumSquares ( const UINT2Sequence *sequence, size_t first, size_t count );
UINT4 XLALUINT4Sum ( const UINT4 *data, size_t first, size_t count );
UINT4 XLALUINT4SumSquares ( const UINT4 *data, size_t first, size_t count );
UINT4 XLALUINT4SequenceSum ( const UINT4Sequence *sequence, size_t first, size_t count );
UINT4 XLALUINT4SequenceSumSquares ( const UINT4Sequence *sequence, size_t first, size_t count );
UINT8 XLALUINT8Sum ( const UINT8 *data, size_t first, size_t count );
UINT8 XLALUINT8SumSquares ( const UINT8 *data, size_t first, size_t count );
UINT8 XLALUINT8SequenceSum ( const UINT8Sequence *sequence, size_t first, size_t count );
UINT8 XLALUINT8SequenceSumSquares ( const UINT8Sequence *sequence, size_t first, size_t count );
/*@}*/

/**
 * \name Conjugate Functions
 *
 * \heading{Synopsis}
 *
 * \code
 * #include <lal/Sequence.h>
 *
 * XLAL<datatype>Conjugate()
 * \endcode
 *
 * \heading{Description}
 *
 * These functions replace a sequence with its complex conjugate.
 *
 */
/*@{*/
COMPLEX8Sequence *XLALConjugateCOMPLEX8Sequence ( COMPLEX8Sequence *series );
COMPLEX16Sequence *XLALConjugateCOMPLEX16Sequence ( COMPLEX16Sequence *series );
/*@}*/

/*@}*/

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif  /* _SEQUENCE_H */
