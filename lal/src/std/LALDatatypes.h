/*
*  Copyright (C) 2007 Jolien Creighton, Kipp Cannon, Reinhard Prix
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/* ---------- SEE LALDatatypes.dox for doxygen documentation ---------- */

#ifndef _LALDATATYPES_H
#define _LALDATATYPES_H

#include <lal/LALAtomicDatatypes.h>

/** \addtogroup LALDatatypes */ /*@{ */

/* ---------- constants ---------- */

/** Type size constants, see \ref LALDatatypes for more details. */
enum {
    LAL_1_BYTE_TYPE_SIZE = 000,  /**< One byte size 00 = 0 */
    LAL_2_BYTE_TYPE_SIZE = 001,  /**< Two byte size 01 = 1 */
    LAL_4_BYTE_TYPE_SIZE = 002,  /**< Four byte size 010 = 2 */
    LAL_8_BYTE_TYPE_SIZE = 003,  /**< Eight byte size 011 = 3 */
    LAL_16_BYTE_TYPE_SIZE = 004, /**< Sixteen byte size 0100 = 4 */
    LAL_TYPE_SIZE_MASK = 007     /**< Type size mask 0111 = 7 */
};

/** Type flag constants, see \ref LALDatatypes for more details. */
enum {
    LAL_FLTPT_TYPE_FLAG = 010,   /**< Floating-point (vs integer) type 01000 =  8 */
    LAL_CMPLX_TYPE_FLAG = 020,   /**< Complex (vs real) type  010000 = 16 */
    LAL_UNSGN_TYPE_FLAG = 040    /**< Unsigned (vs signed) type 0100000 = 32 */
};

/** Type codes: use these type codes to identify a LAL atomic data type, see \ref LALDatatypes for more details. */
typedef enum {
    LAL_CHAR_TYPE_CODE = LAL_1_BYTE_TYPE_SIZE,  /**< CHAR type code (0) */
    LAL_I2_TYPE_CODE = LAL_2_BYTE_TYPE_SIZE,    /**< INT2 type code (1) */
    LAL_I4_TYPE_CODE = LAL_4_BYTE_TYPE_SIZE,    /**< INT4 type code (2) */
    LAL_I8_TYPE_CODE = LAL_8_BYTE_TYPE_SIZE,    /**< INT8 type code (3) */
    LAL_UCHAR_TYPE_CODE = LAL_1_BYTE_TYPE_SIZE | LAL_UNSGN_TYPE_FLAG,   /**< UCHAR type code (32) */
    LAL_U2_TYPE_CODE = LAL_2_BYTE_TYPE_SIZE | LAL_UNSGN_TYPE_FLAG,      /**< UINT2 type code (33) */
    LAL_U4_TYPE_CODE = LAL_4_BYTE_TYPE_SIZE | LAL_UNSGN_TYPE_FLAG,      /**< UINT4 type code (34) */
    LAL_U8_TYPE_CODE = LAL_8_BYTE_TYPE_SIZE | LAL_UNSGN_TYPE_FLAG,      /**< UINT8 type code (35) */
    LAL_S_TYPE_CODE = LAL_4_BYTE_TYPE_SIZE | LAL_FLTPT_TYPE_FLAG,       /**< REAL4 type code (18) */
    LAL_D_TYPE_CODE = LAL_8_BYTE_TYPE_SIZE | LAL_FLTPT_TYPE_FLAG,       /**< REAL8 type code (19) */
    LAL_C_TYPE_CODE = LAL_8_BYTE_TYPE_SIZE | LAL_CMPLX_TYPE_FLAG | LAL_FLTPT_TYPE_FLAG,         /**< COMPLEX8 type code (27) */
    LAL_Z_TYPE_CODE = LAL_16_BYTE_TYPE_SIZE | LAL_CMPLX_TYPE_FLAG | LAL_FLTPT_TYPE_FLAG         /**< COMPLEX16 type code (28) */
} LALTYPECODE;

#ifdef  __cplusplus
extern "C" {
#elif 0
}       /* so that editors will match preceding brace */
#endif

/* ---------- Vector types ---------- */

/** Vector of type CHAR, see \ref ss_Vector for more details. */
typedef struct tagCHARVector {
#ifdef SWIG     /* SWIG interface directives */
    SWIGLAL(ARRAY_1D(CHARVector, CHAR, data, UINT4, length));
#endif /* SWIG */
    UINT4 length; /**< Number of elements in array. */
    CHAR *data; /**< Pointer to the data array. */
} CHARVector;

/** Vector of type CHAR*, ie 'strings', see \ref ss_Vector for more details.  */
typedef struct tagLALStringVector {
#ifdef SWIG     /* SWIG interface directives */
    SWIGLAL(ARRAY_1D(LALStringVector, CHAR *, data, UINT4, length));
#endif /* SWIG */
    UINT4 length; /**< Number of elements in array. */
    CHAR **data; /**< Pointer to the data array. */
} LALStringVector;

/** Vector of type INT2, see \ref ss_Vector for more details. */
typedef struct tagINT2Vector {
#ifdef SWIG     /* SWIG interface directives */
    SWIGLAL(ARRAY_1D(INT2Vector, INT2, data, UINT4, length));
#endif /* SWIG */
    UINT4 length; /**< Number of elements in array. */
    INT2 *data; /**< Pointer to the data array. */
} INT2Vector;

/** Vector of type UINT2, see \ref ss_Vector for more details. */
typedef struct tagUINT2Vector {
#ifdef SWIG     /* SWIG interface directives */
    SWIGLAL(ARRAY_1D(UINT2Vector, UINT2, data, UINT4, length));
#endif /* SWIG */
    UINT4 length; /**< Number of elements in array. */
    UINT2 *data; /**< Pointer to the data array. */
} UINT2Vector;

/** Vector of type INT4, see \ref ss_Vector for more details. */
typedef struct tagINT4Vector {
#ifdef SWIG     /* SWIG interface directives */
    SWIGLAL(ARRAY_1D(INT4Vector, INT4, data, UINT4, length));
#endif /* SWIG */
    UINT4 length; /**< Number of elements in array. */
    INT4 *data; /**< Pointer to the data array. */
} INT4Vector;

/** Vector of type UINT4, see \ref ss_Vector for more details. */
typedef struct tagUINT4Vector {
#ifdef SWIG     /* SWIG interface directives */
    SWIGLAL(ARRAY_1D(UINT4Vector, UINT4, data, UINT4, length));
#endif /* SWIG */
    UINT4 length; /**< Number of elements in array. */
    UINT4 *data; /**< Pointer to the data array. */
} UINT4Vector;

/** Vector of type INT8, see \ref ss_Vector for more details. */
typedef struct tagINT8Vector {
#ifdef SWIG     /* SWIG interface directives */
    SWIGLAL(ARRAY_1D(INT8Vector, INT8, data, UINT4, length));
#endif /* SWIG */
    UINT4 length; /**< Number of elements in array. */
    INT8 *data; /**< Pointer to the data array. */
} INT8Vector;

/** Vector of type UINT8, see \ref ss_Vector for more details. */
typedef struct tagUINT8Vector {
#ifdef SWIG     /* SWIG interface directives */
    SWIGLAL(ARRAY_1D(UINT8Vector, UINT8, data, UINT4, length));
#endif /* SWIG */
    UINT4 length; /**< Number of elements in array. */
    UINT8 *data; /**< Pointer to the data array. */
} UINT8Vector;

/** Vector of type REAL4, see \ref ss_Vector for more details. */
typedef struct tagREAL4Vector {
#ifdef SWIG     /* SWIG interface directives */
    SWIGLAL(ARRAY_1D(REAL4Vector, REAL4, data, UINT4, length));
#endif /* SWIG */
    UINT4 length; /**< Number of elements in array. */
    REAL4 *data; /**< Pointer to the data array. */
} REAL4Vector;

/** Vector of type REAL8, see \ref ss_Vector for more details. */
typedef struct tagREAL8Vector {
#ifdef SWIG     /* SWIG interface directives */
    SWIGLAL(ARRAY_1D(REAL8Vector, REAL8, data, UINT4, length));
#endif /* SWIG */
    UINT4 length; /**< Number of elements in array. */
    REAL8 *data; /**< Pointer to the data array. */
} REAL8Vector;

/** Vector of type COMPLEX8, see \ref ss_Vector for more details. */
typedef struct tagCOMPLEX8Vector {
#ifdef SWIG     /* SWIG interface directives */
    SWIGLAL(ARRAY_1D(COMPLEX8Vector, COMPLEX8, data, UINT4, length));
#endif /* SWIG */
    UINT4 length; /**< Number of elements in array. */
    COMPLEX8 *data; /**< Pointer to the data array. */
} COMPLEX8Vector;

/** Vector of type COMPLEX16, see \ref ss_Vector for more details. */
typedef struct tagCOMPLEX16Vector {
#ifdef SWIG     /* SWIG interface directives */
    SWIGLAL(ARRAY_1D(COMPLEX16Vector, COMPLEX16, data, UINT4, length));
#endif /* SWIG */
    UINT4 length;    /**< Number of elements in array. */
    COMPLEX16 *data; /**< Pointer to the data array. */
} COMPLEX16Vector;


/* ---------- Array types ---------- */

/** Multidimentional array of INT2, see \ref ss_Array for more details. */
typedef struct tagINT2Array {
    UINT4Vector *dimLength; /**< Vector of array dimensions. */
    INT2 *data;      /**< Pointer to the data array. */
} INT2Array;

/** Multidimentional array of UINT2, see \ref ss_Array for more details. */
typedef struct tagUINT2Array {
    UINT4Vector *dimLength; /**< Vector of array dimensions. */
    UINT2 *data;     /**< Pointer to the data array. */
} UINT2Array;

/** Multidimentional array of INT4, see \ref ss_Array for more details. */
typedef struct tagINT4Array {
    UINT4Vector *dimLength; /**< Vector of array dimensions. */
    INT4 *data;      /**< Pointer to the data array. */
} INT4Array;

/** Multidimentional array of UINT4, see \ref ss_Array for more details. */
typedef struct tagUINT4Array {
    UINT4Vector *dimLength; /**< Vector of array dimensions. */
    UINT4 *data;     /**< Pointer to the data array. */
} UINT4Array;

/** Multidimentional array of INT8, see \ref ss_Array for more details. */
typedef struct tagINT8Array {
    UINT4Vector *dimLength; /**< Vector of array dimensions. */
    INT8 *data;      /**< Pointer to the data array. */
} INT8Array;

/** Multidimentional array of UINT8, see \ref ss_Array for more details. */
typedef struct tagUINT8Array {
    UINT4Vector *dimLength; /**< Vector of array dimensions. */
    UINT8 *data;     /**< Pointer to the data array. */
} UINT8Array;

/** Multidimentional array of REAL4, see \ref ss_Array for more details. */
typedef struct tagREAL4Array {
    UINT4Vector *dimLength; /**< Vector of array dimensions. */
    REAL4 *data;     /**< Pointer to the data array. */
} REAL4Array;

/** Multidimentional array of REAL8, see \ref ss_Array for more details. */
typedef struct tagREAL8Array {
    UINT4Vector *dimLength; /**< Vector of array dimensions. */
    REAL8 *data;     /**< Pointer to the data array. */
} REAL8Array;

/** Multidimentional array of COMPLEX8, see \ref ss_Array for more details. */
typedef struct tagCOMPLEX8Array {
    UINT4Vector *dimLength; /**< Vector of array dimensions. */
    COMPLEX8 *data; /**< Pointer to the data array. */
} COMPLEX8Array;

/** Multidimentional array of COMPLEX16, see \ref ss_Array for more details. */
typedef struct tagCOMPLEX16Array {
    UINT4Vector *dimLength; /**< Vector of array dimensions. */
    COMPLEX16 *data; /**< Pointer to the data array. */
} COMPLEX16Array;


    /* ---------- Sequence types ---------- */

typedef CHARVector CHARSequence;        /**< See \ref ss_Sequence for documentation */
typedef INT2Vector INT2Sequence;        /**< See \ref ss_Sequence for documentation */
typedef UINT2Vector UINT2Sequence;      /**< See \ref ss_Sequence for documentation */
typedef INT4Vector INT4Sequence;        /**< See \ref ss_Sequence for documentation */
typedef UINT4Vector UINT4Sequence;      /**< See \ref ss_Sequence for documentation */
typedef INT8Vector INT8Sequence;        /**< See \ref ss_Sequence for documentation */
typedef UINT8Vector UINT8Sequence;      /**< See \ref ss_Sequence for documentation */
typedef REAL4Vector REAL4Sequence;      /**< See \ref ss_Sequence for documentation */
typedef REAL8Vector REAL8Sequence;      /**< See \ref ss_Sequence for documentation */
typedef COMPLEX8Vector COMPLEX8Sequence; /**< See \ref ss_Sequence for documentation */
typedef COMPLEX16Vector COMPLEX16Sequence; /**< See \ref ss_Sequence for documentation */

    /* ---------- VectorSequence types ---------- */

/** Sequence of CHAR Vectors, see \ref ss_VectorSequence for more details. */
typedef struct tagCHARVectorSequence {
#ifdef SWIG     /* SWIG interface directives */
    SWIGLAL(ARRAY_2D(CHARVectorSequence, CHAR, data, UINT4, length, vectorLength));
#endif /* SWIG */
    UINT4 length; /**< The number \a l of vectors. */
    UINT4 vectorLength; /**< The length \a n of each vector. */
    CHAR *data; /**< Pointer to the data array.  Element \a i of vector \a j is \c data[ \a jn + \a i \c ]. */
} CHARVectorSequence;

/** Sequence of INT2 Vectors, see \ref ss_VectorSequence for more details. */
typedef struct tagINT2VectorSequence {
#ifdef SWIG     /* SWIG interface directives */
    SWIGLAL(ARRAY_2D(INT2VectorSequence, INT2, data, UINT4, length, vectorLength));
#endif /* SWIG */
    UINT4 length; /**< The number \a l of vectors. */
    UINT4 vectorLength; /**< The length \a n of each vector. */
    INT2 *data; /**< Pointer to the data array.  Element \a i of vector \a j is \c data[ \a jn + \a i \c ]. */
} INT2VectorSequence;

/** Sequence of UINT2 Vectors, see \ref ss_VectorSequence for more details. */
typedef struct tagUINT2VectorSequence {
#ifdef SWIG     /* SWIG interface directives */
    SWIGLAL(ARRAY_2D(UINT2VectorSequence, UINT2, data, UINT4, length, vectorLength));
#endif /* SWIG */
    UINT4 length; /**< The number \a l of vectors. */
    UINT4 vectorLength; /**< The length \a n of each vector. */
    UINT2 *data; /**< Pointer to the data array.  Element \a i of vector \a j is \c data[ \a jn + \a i \c ]. */
} UINT2VectorSequence;

/** Sequence of INT4 Vectors, see \ref ss_VectorSequence for more details. */
typedef struct tagINT4VectorSequence {
#ifdef SWIG     /* SWIG interface directives */
    SWIGLAL(ARRAY_2D(INT4VectorSequence, INT4, data, UINT4, length, vectorLength));
#endif /* SWIG */
    UINT4 length; /**< The number \a l of vectors. */
    UINT4 vectorLength; /**< The length \a n of each vector. */
    INT4 *data; /**< Pointer to the data array.  Element \a i of vector \a j is \c data[ \a jn + \a i \c ]. */
} INT4VectorSequence;

/** Sequence of UINT4 Vectors, see \ref ss_VectorSequence for more details. */
typedef struct tagUINT4VectorSequence {
#ifdef SWIG     /* SWIG interface directives */
    SWIGLAL(ARRAY_2D(UINT4VectorSequence, UINT4, data, UINT4, length, vectorLength));
#endif /* SWIG */
    UINT4 length; /**< The number \a l of vectors. */
    UINT4 vectorLength; /**< The length \a n of each vector. */
    UINT4 *data; /**< Pointer to the data array.  Element \a i of vector \a j is \c data[ \a jn + \a i \c ]. */
} UINT4VectorSequence;

/** Sequence of INT8 Vectors, see \ref ss_VectorSequence for more details. */
typedef struct tagINT8VectorSequence {
#ifdef SWIG     /* SWIG interface directives */
    SWIGLAL(ARRAY_2D(INT8VectorSequence, INT8, data, UINT4, length, vectorLength));
#endif /* SWIG */
    UINT4 length; /**< The number \a l of vectors. */
    UINT4 vectorLength; /**< The length \a n of each vector. */
    INT8 *data; /**< Pointer to the data array.  Element \a i of vector \a j is \c data[ \a jn + \a i \c ]. */
} INT8VectorSequence;

/** Sequence of UINT8 Vectors, see \ref ss_VectorSequence for more details. */
typedef struct tagUINT8VectorSequence {
#ifdef SWIG     /* SWIG interface directives */
    SWIGLAL(ARRAY_2D(UINT8VectorSequence, UINT8, data, UINT4, length, vectorLength));
#endif /* SWIG */
    UINT4 length; /**< The number \a l of vectors. */
    UINT4 vectorLength; /**< The length \a n of each vector. */
    UINT8 *data; /**< Pointer to the data array.  Element \a i of vector \a j is \c data[ \a jn + \a i \c ]. */
} UINT8VectorSequence;

/** Sequence of REAL4 Vectors, see \ref ss_VectorSequence for more details. */
typedef struct tagREAL4VectorSequence {
#ifdef SWIG     /* SWIG interface directives */
    SWIGLAL(ARRAY_2D(REAL4VectorSequence, REAL4, data, UINT4, length, vectorLength));
#endif /* SWIG */
    UINT4 length; /**< The number \a l of vectors. */
    UINT4 vectorLength; /**< The length \a n of each vector. */
    REAL4 *data; /**< Pointer to the data array.  Element \a i of vector \a j is \c data[ \a jn + \a i \c ]. */
} REAL4VectorSequence;

/** Sequence of REAL8 Vectors, see \ref ss_VectorSequence for more details. */
typedef struct tagREAL8VectorSequence {
#ifdef SWIG     /* SWIG interface directives */
    SWIGLAL(ARRAY_2D(REAL8VectorSequence, REAL8, data, UINT4, length, vectorLength));
#endif /* SWIG */
    UINT4 length; /**< The number \a l of vectors. */
    UINT4 vectorLength; /**< The length \a n of each vector. */
    REAL8 *data; /**< Pointer to the data array.  Element \a i of vector \a j is \c data[ \a jn + \a i \c ]. */
} REAL8VectorSequence;

/** Sequence of COMPLEX8 Vectors, see \ref ss_VectorSequence for more details. */
typedef struct tagCOMPLEX8VectorSequence {
#ifdef SWIG     /* SWIG interface directives */
    SWIGLAL(ARRAY_2D(COMPLEX8VectorSequence, COMPLEX8, data, UINT4, length, vectorLength));
#endif /* SWIG */
    UINT4 length; /**< The number \a l of vectors. */
    UINT4 vectorLength; /**< The length \a n of each vector. */
    COMPLEX8 *data; /**< Pointer to the data array.  Element \a i of vector \a j is \c data[ \a jn + \a i \c ]. */
} COMPLEX8VectorSequence;

/** Sequence of COMPLEX16 Vectors, see \ref ss_VectorSequence for more details. */
typedef struct tagCOMPLEX16VectorSequence {
#ifdef SWIG     /* SWIG interface directives */
    SWIGLAL(ARRAY_2D(COMPLEX16VectorSequence, COMPLEX16, data, UINT4, length, vectorLength));
#endif /* SWIG */
    UINT4 length;    /**< The number \a l of vectors. */
    UINT4 vectorLength;    /**< The length \a n of each vector. */
    COMPLEX16 *data; /**< Pointer to the data array.  Element \a i of vector \a j is \c data[ \a jn + \a i \c ]. */
} COMPLEX16VectorSequence;

/* ---------- ArraySequence types ---------- */

/** Sequence of INT2 multidimensional arrays, see \ref ss_ArraySequence for more details. */
typedef struct tagINT2ArraySequence {
    UINT4 length;      /**< The number \a l of vectors. */
    UINT4 arrayDim;      /**< The number of data \a N in each array element (this is not the number \a m of indices). */
    UINT4Vector *dimLength; /**< Pointer to a vector of length \a m storing the array dimensions */
    INT2 *data;      /**< Pointer to the data array. */
} INT2ArraySequence;

/** Sequence of UINT2 multidimensional arrays, see \ref ss_ArraySequence for more details. */
typedef struct tagUINT2ArraySequence {
    UINT4 length;      /**< The number \a l of vectors. */
    UINT4 arrayDim;      /**< The number of data \a N in each array element (this is not the number \a m of indices). */
    UINT4Vector *dimLength; /**< Pointer to a vector of length \a m storing the array dimensions */
    UINT2 *data;     /**< Pointer to the data array. */
} UINT2ArraySequence;

/** Sequence of INT4 multidimensional arrays, see \ref ss_ArraySequence for more details. */
typedef struct tagINT4ArraySequence {
    UINT4 length;      /**< The number \a l of vectors. */
    UINT4 arrayDim;      /**< The number of data \a N in each array element (this is not the number \a m of indices). */
    UINT4Vector *dimLength; /**< Pointer to a vector of length \a m storing the array dimensions */
    INT4 *data;      /**< Pointer to the data array. */
} INT4ArraySequence;

/** Sequence of UINT4 multidimensional arrays, see \ref ss_ArraySequence for more details. */
typedef struct tagUINT4ArraySequence {
    UINT4 length;      /**< The number \a l of vectors. */
    UINT4 arrayDim;      /**< The number of data \a N in each array element (this is not the number \a m of indices). */
    UINT4Vector *dimLength; /**< Pointer to a vector of length \a m storing the array dimensions */
    UINT4 *data;     /**< Pointer to the data array. */
} UINT4ArraySequence;

/** Sequence of INT8 multidimensional arrays, see \ref ss_ArraySequence for more details. */
typedef struct tagINT8ArraySequence {
    UINT4 length;      /**< The number \a l of vectors. */
    UINT4 arrayDim;      /**< The number of data \a N in each array element (this is not the number \a m of indices). */
    UINT4Vector *dimLength; /**< Pointer to a vector of length \a m storing the array dimensions */
    INT8 *data;      /**< Pointer to the data array. */
} INT8ArraySequence;

/** Sequence of UINT8 multidimensional arrays, see \ref ss_ArraySequence for more details. */
typedef struct tagUINT8ArraySequence {
    UINT4 length;      /**< The number \a l of vectors. */
    UINT4 arrayDim;      /**< The number of data \a N in each array element (this is not the number \a m of indices). */
    UINT4Vector *dimLength; /**< Pointer to a vector of length \a m storing the array dimensions */
    UINT8 *data;     /**< Pointer to the data array. */
} UINT8ArraySequence;

/** Sequence of REAL4 multidimensional arrays, see \ref ss_ArraySequence for more details. */
typedef struct tagREAL4ArraySequence {
    UINT4 length;      /**< The number \a l of vectors. */
    UINT4 arrayDim;      /**< The number of data \a N in each array element (this is not the number \a m of indices). */
    UINT4Vector *dimLength; /**< Pointer to a vector of length \a m storing the array dimensions */
    REAL4 *data;     /**< Pointer to the data array. */
} REAL4ArraySequence;

/** Sequence of REAL8 multidimensional arrays, see \ref ss_ArraySequence for more details. */
typedef struct tagREAL8ArraySequence {
    UINT4 length;      /**< The number \a l of vectors. */
    UINT4 arrayDim;      /**< The number of data \a N in each array element (this is not the number \a m of indices). */
    UINT4Vector *dimLength; /**< Pointer to a vector of length \a m storing the array dimensions */
    REAL8 *data;     /**< Pointer to the data array. */
} REAL8ArraySequence;

/** Sequence of COMPLEX8 multidimensional arrays, see \ref ss_ArraySequence for more details. */
typedef struct tagCOMPLEX8ArraySequence {
    UINT4 length;      /**< The number \a l of vectors. */
    UINT4 arrayDim;      /**< The number of data \a N in each array element (this is not the number \a m of indices). */
    UINT4Vector *dimLength; /**< Pointer to a vector of length \a m storing the array dimensions */
    COMPLEX8 *data; /**< Pointer to the data array. */
} COMPLEX8ArraySequence;

/** Sequence of COMPLEX16 multidimensional arrays, see \ref ss_ArraySequence for more details. */
typedef struct tagCOMPLEX16ArraySequence {
    UINT4 length;      /**< The number \a l of vectors. */
    UINT4 arrayDim;      /**< The number of data \a N in each array element (this is not the number \a m of indices). */
    UINT4Vector *dimLength; /**< Pointer to a vector of length \a m storing the array dimensions */
    COMPLEX16 *data; /**< Pointer to the data array. */
} COMPLEX16ArraySequence;

    /* ---------- Structured datatypes ---------- */

/** Epoch relative to GPS epoch, see \ref ss_LIGOTimeGPS for more details */
#ifdef SWIG     /* SWIG interface directives */
SWIGLAL(IMMUTABLE_MEMBERS(tagLIGOTimeGPS, gpsSeconds, gpsNanoSeconds));
#endif /* SWIG */
typedef struct tagLIGOTimeGPS {
    INT4 gpsSeconds; /**< Seconds since 0h UTC 6 Jan 1980. */
    INT4 gpsNanoSeconds; /**< Residual nanoseconds. */
} LIGOTimeGPS;

/** Zero-initializer for LIGOTimeGPS structs */
#define LIGOTIMEGPSZERO { 0, 0 }

/**
 * Indices of arrays corresponding to particular units.
 * The ::LALUnit structure has arrays giving the numerators
 * and denominators-minus-one of the powers of various units.
 * These are the indices for the particular units.
 */
enum {
    LALUnitIndexMeter, /**< The meter index. */
    LALUnitIndexKiloGram, /**< The kilogram index. */
    LALUnitIndexSecond, /**< The second index. */
    LALUnitIndexAmpere, /**< The ampere index. */
    LALUnitIndexKelvin, /**< The kelvin index. */
    LALUnitIndexStrain, /**< The strain index. */
    LALUnitIndexADCCount, /**< The ADC counts index. */
    LALNumUnits         /**< The number of units. */
};

/**
 * This structure stores units in the mksA system (plus Kelvin, Strain,
 * and ADC Count).  It also stores an overall power-of-ten scaling factor.
 * Thus, the units are given by
 * \f{equation}{
 * 10^p\times\textrm{m}^{N_0/(1+D_0)}\times\textrm{kg}^{N_1/(1+D_1)}
 * \times\textrm{s}^{N_2/(1+D_2)}\times\textrm{A}^{N_3/(1+D_3)}
 * \times\textrm{K}^{N_4/(1+D_4)}\times\textrm{strain}^{N_5/(1+D_5)}
 * \times\textrm{count}^{N_6/(1+D_6)}
 * \f}
 *
 */
#ifdef SWIG     /* SWIG interface directives */
SWIGLAL(IMMUTABLE_MEMBERS(tagLALUnit, powerOfTen, unitNumerator, unitDenominatorMinusOne));
#endif /* SWIG */
typedef struct tagLALUnit {
    INT2 powerOfTen; /**< Overall power-of-ten scaling is 10^\c powerOfTen. */
    INT2 unitNumerator[LALNumUnits]; /**< Array of unit power numerators. */
    UINT2 unitDenominatorMinusOne[LALNumUnits]; /**< Array of unit power denominators-minus-one. */
} LALUnit;


    /* ---------- TimeSeries types ---------- */

/** Length of name fields of LAL structured data types. */
enum enumLALNameLength { LALNameLength = 64 };

/** Time series of INT2 data, see \ref ss_TimeSeries for more details. */
typedef struct tagINT2TimeSeries {
    CHAR name[LALNameLength];        /**< The name of the time series. */
    LIGOTimeGPS epoch; /**< The start time of the time series. */
    REAL8 deltaT;       /**< The time step between samples of the time series in seconds. */
    REAL8 f0;       /**< The heterodyning frequency, in Hertz (zero if not heterodyned). */
    LALUnit sampleUnits;     /**< The physical units of the quantity being sampled. */
    INT2Sequence *data; /**< The sequence of sampled data. */
} INT2TimeSeries;

/** Time series of UINT2 data, see \ref ss_TimeSeries for more details. */
typedef struct tagUINT2TimeSeries {
    CHAR name[LALNameLength];         /**< The name of the time series. */
    LIGOTimeGPS epoch; /**< The start time of the time series. */
    REAL8 deltaT;        /**< The time step between samples of the time series in seconds. */ REAL8 f0;        /**< The heterodyning frequency, in Hertz (zero if not heterodyned). */
    LALUnit sampleUnits;      /**< The physical units of the quantity being sampled. */
    UINT2Sequence *data; /**< The sequence of sampled data. */
} UINT2TimeSeries;

/** Time series of INT4 data, see \ref ss_TimeSeries for more details. */
typedef struct tagINT4TimeSeries {
    CHAR name[LALNameLength];        /**< The name of the time series. */
    LIGOTimeGPS epoch; /**< The start time of the time series. */
    REAL8 deltaT;       /**< The time step between samples of the time series in seconds. */
    REAL8 f0;       /**< The heterodyning frequency, in Hertz (zero if not heterodyned). */
    LALUnit sampleUnits;     /**< The physical units of the quantity being sampled. */
    INT4Sequence *data; /**< The sequence of sampled data. */
} INT4TimeSeries;

/** Time series of UINT4 data, see \ref ss_TimeSeries for more details. */
typedef struct tagUINT4TimeSeries {
    CHAR name[LALNameLength];         /**< The name of the time series. */
    LIGOTimeGPS epoch; /**< The start time of the time series. */
    REAL8 deltaT;        /**< The time step between samples of the time series in seconds. */
    REAL8 f0;        /**< The heterodyning frequency, in Hertz (zero if not heterodyned). */
    LALUnit sampleUnits;      /**< The physical units of the quantity being sampled. */
    UINT4Sequence *data; /**< The sequence of sampled data. */
} UINT4TimeSeries;

/** Time series of INT8 data, see \ref ss_TimeSeries for more details. */
typedef struct tagINT8TimeSeries {
    CHAR name[LALNameLength];        /**< The name of the time series. */
    LIGOTimeGPS epoch; /**< The start time of the time series. */
    REAL8 deltaT;       /**< The time step between samples of the time series in seconds. */
    REAL8 f0;       /**< The heterodyning frequency, in Hertz (zero if not heterodyned). */
    LALUnit sampleUnits;     /**< The physical units of the quantity being sampled. */
    INT8Sequence *data; /**< The sequence of sampled data. */
} INT8TimeSeries;

/** Time series of UINT8 data, see \ref ss_TimeSeries for more details. */
typedef struct tagUINT8TimeSeries {
    CHAR name[LALNameLength];         /**< The name of the time series. */
    LIGOTimeGPS epoch; /**< The start time of the time series. */
    REAL8 deltaT;        /**< The time step between samples of the time series in seconds. */
    REAL8 f0;        /**< The heterodyning frequency, in Hertz (zero if not heterodyned). */
    LALUnit sampleUnits;      /**< The physical units of the quantity being sampled. */
    UINT8Sequence *data; /**< The sequence of sampled data. */
} UINT8TimeSeries;

/** Time series of REAL4 data, see \ref ss_TimeSeries for more details. */
typedef struct tagREAL4TimeSeries {
    CHAR name[LALNameLength];         /**< The name of the time series. */
    LIGOTimeGPS epoch; /**< The start time of the time series. */
    REAL8 deltaT;        /**< The time step between samples of the time series in seconds. */
    REAL8 f0;        /**< The heterodyning frequency, in Hertz (zero if not heterodyned). */
    LALUnit sampleUnits;      /**< The physical units of the quantity being sampled. */
    REAL4Sequence *data; /**< The sequence of sampled data. */
} REAL4TimeSeries;

/** Time series of REAL8 data, see \ref ss_TimeSeries for more details. */
typedef struct tagREAL8TimeSeries {
    CHAR name[LALNameLength];         /**< The name of the time series. */
    LIGOTimeGPS epoch; /**< The start time of the time series. */
    REAL8 deltaT;        /**< The time step between samples of the time series in seconds. */
    REAL8 f0;        /**< The heterodyning frequency, in Hertz (zero if not heterodyned). */
    LALUnit sampleUnits;      /**< The physical units of the quantity being sampled. */
    REAL8Sequence *data; /**< The sequence of sampled data. */
} REAL8TimeSeries;

/** Time series of COMPLEX8 data, see \ref ss_TimeSeries for more details. */
typedef struct tagCOMPLEX8TimeSeries {
    CHAR name[LALNameLength];            /**< The name of the time series. */
    LIGOTimeGPS epoch;     /**< The start time of the time series. */
    REAL8 deltaT;           /**< The time step between samples of the time series in seconds. */
    REAL8 f0;           /**< The heterodyning frequency, in Hertz (zero if not heterodyned). */
    LALUnit sampleUnits;         /**< The physical units of the quantity being sampled. */
    COMPLEX8Sequence *data; /**< The sequence of sampled data. */
} COMPLEX8TimeSeries;

/** Time series of COMPLEX16 data, see \ref ss_TimeSeries for more details. */
typedef struct tagCOMPLEX16TimeSeries {
    CHAR name[LALNameLength];             /**< The name of the time series. */
    LIGOTimeGPS epoch;      /**< The start time of the time series. */
    REAL8 deltaT;            /**< The time step between samples of the time series in seconds. */
    REAL8 f0;            /**< The heterodyning frequency, in Hertz (zero if not heterodyned). */
    LALUnit sampleUnits;          /**< The physical units of the quantity being sampled. */
    COMPLEX16Sequence *data; /**< The sequence of sampled data. */
} COMPLEX16TimeSeries;


    /* ---------- TimeVectorSeries types ---------- */

/** Time series of INT2 vectors, see \ref ss_TimeVectorSeries for more details. */
typedef struct tagINT2TimeVectorSeries {
    CHAR name[LALNameLength];               /**< The name of the time series of vectors. */
    LIGOTimeGPS epoch;        /**< The start time of the time series of vectors. */
    REAL8 deltaT;              /**< The time step between samples of the time series of vectors in seconds. */
    REAL8 f0;              /**< The heterodyning frequency, in Hertz (zero if not heterodyned). */
    LALUnit sampleUnits;            /**< The physical units of the quantity being sampled. */
    INT2VectorSequence *data; /**< The sequence of sampled data vectors. */
} INT2TimeVectorSeries;

/** Time series of UINT2 vectors, see \ref ss_TimeVectorSeries for more details. */
typedef struct tagUINT2TimeVectorSeries {
    CHAR name[LALNameLength];               /**< The name of the time series of vectors. */
    LIGOTimeGPS epoch;        /**< The start time of the time series of vectors. */
    REAL8 deltaT;              /**< The time step between samples of the time series of vectors in seconds. */
    REAL8 f0;              /**< The heterodyning frequency, in Hertz (zero if not heterodyned). */
    LALUnit sampleUnits;            /**< The physical units of the quantity being sampled. */
    UINT2VectorSequence *data; /**< The sequence of sampled data vectors. */
} UINT2TimeVectorSeries;

/** Time series of INT4 vectors, see \ref ss_TimeVectorSeries for more details. */
typedef struct tagINT4TimeVectorSeries {
    CHAR name[LALNameLength];               /**< The name of the time series of vectors. */
    LIGOTimeGPS epoch;        /**< The start time of the time series of vectors. */
    REAL8 deltaT;              /**< The time step between samples of the time series of vectors in seconds. */
    REAL8 f0;              /**< The heterodyning frequency, in Hertz (zero if not heterodyned). */
    LALUnit sampleUnits;            /**< The physical units of the quantity being sampled. */
    INT4VectorSequence *data; /**< The sequence of sampled data vectors. */
} INT4TimeVectorSeries;

/** Time series of UINT4 vectors, see \ref ss_TimeVectorSeries for more details. */
typedef struct tagUINT4TimeVectorSeries {
    CHAR name[LALNameLength];               /**< The name of the time series of vectors. */
    LIGOTimeGPS epoch;        /**< The start time of the time series of vectors. */
    REAL8 deltaT;              /**< The time step between samples of the time series of vectors in seconds. */
    REAL8 f0;              /**< The heterodyning frequency, in Hertz (zero if not heterodyned). */
    LALUnit sampleUnits;            /**< The physical units of the quantity being sampled. */
    UINT4VectorSequence *data; /**< The sequence of sampled data vectors. */
} UINT4TimeVectorSeries;

/** Time series of INT8 vectors, see \ref ss_TimeVectorSeries for more details. */
typedef struct tagINT8TimeVectorSeries {
    CHAR name[LALNameLength];               /**< The name of the time series of vectors. */
    LIGOTimeGPS epoch;        /**< The start time of the time series of vectors. */
    REAL8 deltaT;              /**< The time step between samples of the time series of vectors in seconds. */
    REAL8 f0;              /**< The heterodyning frequency, in Hertz (zero if not heterodyned). */
    LALUnit sampleUnits;            /**< The physical units of the quantity being sampled. */
    INT8VectorSequence *data; /**< The sequence of sampled data vectors. */
} INT8TimeVectorSeries;

/** Time series of UINT8 vectors, see \ref ss_TimeVectorSeries for more details. */
typedef struct tagUINT8TimeVectorSeries {
    CHAR name[LALNameLength];               /**< The name of the time series of vectors. */
    LIGOTimeGPS epoch;        /**< The start time of the time series of vectors. */
    REAL8 deltaT;              /**< The time step between samples of the time series of vectors in seconds. */
    REAL8 f0;              /**< The heterodyning frequency, in Hertz (zero if not heterodyned). */
    LALUnit sampleUnits;            /**< The physical units of the quantity being sampled. */
    UINT8VectorSequence *data; /**< The sequence of sampled data vectors. */
} UINT8TimeVectorSeries;

/** Time series of REAL4 vectors, see \ref ss_TimeVectorSeries for more details. */
typedef struct tagREAL4TimeVectorSeries {
    CHAR name[LALNameLength];               /**< The name of the time series of vectors. */
    LIGOTimeGPS epoch;        /**< The start time of the time series of vectors. */
    REAL8 deltaT;              /**< The time step between samples of the time series of vectors in seconds. */
    REAL8 f0;              /**< The heterodyning frequency, in Hertz (zero if not heterodyned). */
    LALUnit sampleUnits;            /**< The physical units of the quantity being sampled. */
    REAL4VectorSequence *data; /**< The sequence of sampled data vectors. */
} REAL4TimeVectorSeries;

/** Time series of REAL8 vectors, see \ref ss_TimeVectorSeries for more details. */
typedef struct tagREAL8TimeVectorSeries {
    CHAR name[LALNameLength];               /**< The name of the time series of vectors. */
    LIGOTimeGPS epoch;        /**< The start time of the time series of vectors. */
    REAL8 deltaT;              /**< The time step between samples of the time series of vectors in seconds. */
    REAL8 f0;              /**< The heterodyning frequency, in Hertz (zero if not heterodyned). */
    LALUnit sampleUnits;            /**< The physical units of the quantity being sampled. */
    REAL8VectorSequence *data; /**< The sequence of sampled data vectors. */
} REAL8TimeVectorSeries;

/** Time series of COMPLEX8 vectors, see \ref ss_TimeVectorSeries for more details. */
typedef struct tagCOMPLEX8TimeVectorSeries {
    CHAR name[LALNameLength];                   /**< The name of the time series of vectors. */
    LIGOTimeGPS epoch;            /**< The start time of the time series of vectors. */
    REAL8 deltaT;                  /**< The time step between samples of the time series of vectors in seconds. */
    REAL8 f0;                  /**< The heterodyning frequency, in Hertz (zero if not heterodyned). */
    LALUnit sampleUnits;                /**< The physical units of the quantity being sampled. */
    COMPLEX8VectorSequence *data; /**< The sequence of sampled data vectors. */
} COMPLEX8TimeVectorSeries;

/** Time series of COMPLEX16 vectors, see \ref ss_TimeVectorSeries for more details. */
typedef struct tagCOMPLEX16TimeVectorSeries {
    CHAR name[LALNameLength];                    /**< The name of the time series of vectors. */
    LIGOTimeGPS epoch;             /**< The start time of the time series of vectors. */
    REAL8 deltaT;                   /**< The time step between samples of the time series of vectors in seconds. */
    REAL8 f0;                   /**< The heterodyning frequency, in Hertz (zero if not heterodyned). */
    LALUnit sampleUnits;                 /**< The physical units of the quantity being sampled. */
    COMPLEX16VectorSequence *data; /**< The sequence of sampled data vectors. */
} COMPLEX16TimeVectorSeries;


/* ---------- TimeArraySeries ---------- */

/** See \ref ss_TimeArraySeries for documentation */
typedef struct tagINT2TimeArraySeries {
    CHAR name[LALNameLength];
    LIGOTimeGPS epoch;
    REAL8 deltaT;
    REAL8 f0;
    LALUnit sampleUnits;
    INT2ArraySequence *data;
} INT2TimeArraySeries;

/** See \ref ss_TimeArraySeries for documentation */
typedef struct tagUINT2TimeArraySeries {
    CHAR name[LALNameLength];
    LIGOTimeGPS epoch;
    REAL8 deltaT;
    REAL8 f0;
    LALUnit sampleUnits;
    UINT2ArraySequence *data;
} UINT2TimeArraySeries;

/** See \ref ss_TimeArraySeries for documentation */
typedef struct tagINT4TimeArraySeries {
    CHAR name[LALNameLength];
    LIGOTimeGPS epoch;
    REAL8 deltaT;
    REAL8 f0;
    LALUnit sampleUnits;
    INT4ArraySequence *data;
} INT4TimeArraySeries;

/** See \ref ss_TimeArraySeries for documentation */
typedef struct tagUINT4TimeArraySeries {
    CHAR name[LALNameLength];
    LIGOTimeGPS epoch;
    REAL8 deltaT;
    REAL8 f0;
    LALUnit sampleUnits;
    UINT4ArraySequence *data;
} UINT4TimeArraySeries;

/** See \ref ss_TimeArraySeries for documentation */
typedef struct tagINT8TimeArraySeries {
    CHAR name[LALNameLength];
    LIGOTimeGPS epoch;
    REAL8 deltaT;
    REAL8 f0;
    LALUnit sampleUnits;
    INT8ArraySequence *data;
} INT8TimeArraySeries;

/** See \ref ss_TimeArraySeries for documentation */
typedef struct tagUINT8TimeArraySeries {
    CHAR name[LALNameLength];
    LIGOTimeGPS epoch;
    REAL8 deltaT;
    REAL8 f0;
    LALUnit sampleUnits;
    UINT8ArraySequence *data;
} UINT8TimeArraySeries;

/** See \ref ss_TimeArraySeries for documentation */
typedef struct tagREAL4TimeArraySeries {
    CHAR name[LALNameLength];
    LIGOTimeGPS epoch;
    REAL8 deltaT;
    REAL8 f0;
    LALUnit sampleUnits;
    REAL4ArraySequence *data;
} REAL4TimeArraySeries;

/** See \ref ss_TimeArraySeries for documentation */
typedef struct tagREAL8TimeArraySeries {
    CHAR name[LALNameLength];
    LIGOTimeGPS epoch;
    REAL8 deltaT;
    REAL8 f0;
    LALUnit sampleUnits;
    REAL8ArraySequence *data;
} REAL8TimeArraySeries;

/** See \ref ss_TimeArraySeries for documentation */
typedef struct tagCOMPLEX8TimeArraySeries {
    CHAR name[LALNameLength];
    LIGOTimeGPS epoch;
    REAL8 deltaT;
    REAL8 f0;
    LALUnit sampleUnits;
    COMPLEX8ArraySequence *data;
} COMPLEX8TimeArraySeries;

/** See \ref ss_TimeArraySeries for documentation */
typedef struct tagCOMPLEX16TimeArraySeries {
    CHAR name[LALNameLength];
    LIGOTimeGPS epoch;
    REAL8 deltaT;
    REAL8 f0;
    LALUnit sampleUnits;
    COMPLEX16ArraySequence *data;
} COMPLEX16TimeArraySeries;


/* ---------- FrequencySeries types ---------- */

/** See \ref ss_FrequencySeries for documentation */
typedef struct tagINT2FrequencySeries {
    CHAR name[LALNameLength];
    LIGOTimeGPS epoch;
    REAL8 f0;
    REAL8 deltaF;
    LALUnit sampleUnits;
    INT2Sequence *data;
} INT2FrequencySeries;

/** See \ref ss_FrequencySeries for documentation */
typedef struct tagUINT2FrequencySeries {
    CHAR name[LALNameLength];
    LIGOTimeGPS epoch;
    REAL8 f0;
    REAL8 deltaF;
    LALUnit sampleUnits;
    UINT2Sequence *data;
} UINT2FrequencySeries;

/** See \ref ss_FrequencySeries for documentation */
typedef struct tagINT4FrequencySeries {
    CHAR name[LALNameLength];
    LIGOTimeGPS epoch;
    REAL8 f0;
    REAL8 deltaF;
    LALUnit sampleUnits;
    INT4Sequence *data;
} INT4FrequencySeries;

/** See \ref ss_FrequencySeries for documentation */
typedef struct tagUINT4FrequencySeries {
    CHAR name[LALNameLength];
    LIGOTimeGPS epoch;
    REAL8 f0;
    REAL8 deltaF;
    LALUnit sampleUnits;
    UINT4Sequence *data;
} UINT4FrequencySeries;

/** See \ref ss_FrequencySeries for documentation */
typedef struct tagINT8FrequencySeries {
    CHAR name[LALNameLength];
    LIGOTimeGPS epoch;
    REAL8 f0;
    REAL8 deltaF;
    LALUnit sampleUnits;
    INT8Sequence *data;
} INT8FrequencySeries;

/** See \ref ss_FrequencySeries for documentation */
typedef struct tagUINT8FrequencySeries {
    CHAR name[LALNameLength];
    LIGOTimeGPS epoch;
    REAL8 f0;
    REAL8 deltaF;
    LALUnit sampleUnits;
    UINT8Sequence *data;
} UINT8FrequencySeries;

/** See \ref ss_FrequencySeries for documentation */
typedef struct tagREAL4FrequencySeries {
    CHAR name[LALNameLength];
    LIGOTimeGPS epoch;
    REAL8 f0;
    REAL8 deltaF;
    LALUnit sampleUnits;
    REAL4Sequence *data;
} REAL4FrequencySeries;

/** See \ref ss_FrequencySeries for documentation */
typedef struct tagREAL8FrequencySeries {
    CHAR name[LALNameLength];
    LIGOTimeGPS epoch;
    REAL8 f0;
    REAL8 deltaF;
    LALUnit sampleUnits;
    REAL8Sequence *data;
} REAL8FrequencySeries;

/** See \ref ss_FrequencySeries for documentation */
typedef struct tagCOMPLEX8FrequencySeries {
    CHAR name[LALNameLength];
    LIGOTimeGPS epoch;
    REAL8 f0;
    REAL8 deltaF;
    LALUnit sampleUnits;
    COMPLEX8Sequence *data;
} COMPLEX8FrequencySeries;

/** See \ref ss_FrequencySeries for documentation */
typedef struct tagCOMPLEX16FrequencySeries {
    CHAR name[LALNameLength];
    LIGOTimeGPS epoch;
    REAL8 f0;
    REAL8 deltaF;
    LALUnit sampleUnits;
    COMPLEX16Sequence *data;
} COMPLEX16FrequencySeries;

/* ---------- ZPGFilter types ---------- */

/** See \ref ss_ZPGFilter for details */
typedef struct tagCOMPLEX8ZPGFilter {
    CHAR name[LALNameLength];
    REAL8 deltaT;
    COMPLEX8Vector *zeros;
    COMPLEX8Vector *poles;
    COMPLEX8 gain;
} COMPLEX8ZPGFilter;

/** See \ref ss_ZPGFilter for details */
typedef struct tagCOMPLEX16ZPGFilter {
    CHAR name[LALNameLength];
    REAL8 deltaT;
    COMPLEX16Vector *zeros;
    COMPLEX16Vector *poles;
    COMPLEX16 gain;
} COMPLEX16ZPGFilter;

           /*@} */ /* end of LALDatatypes documentation group */


#ifndef SWIG    /* exclude from SWIG interface */

/**
 * \ingroup LALStatusMacros_h
 * \brief LAL status structure, see \ref ss_LALStatus for more details.
 */
typedef struct tagLALStatus {
    INT4 statusCode;                            /**< A numerical code identifying the type of error, or 0 for nominal status; Negative values are reserved for certain standard error types */
    const CHAR *statusDescription;              /**< An explanatory string corresponding to the numerical status code */
    volatile const CHAR *Id;                    /**< A character string identifying the source file and version number of the function being reported on */
    const CHAR *function;                       /**< The name of the function */
    const CHAR *file;                           /**< The name of the source file containing the function code */
    INT4 line;                                  /**< The line number in the source file where the current \c statusCode was set */
    struct tagLALStatus *statusPtr;             /**< Pointer to the next node in the list; \c NULL if this function is not reporting a subroutine error */
    INT4 level;                                 /**< The nested-function level where any error was reported */
} LALStatus;
#endif /* SWIG */

#if 0
{       /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALDATATYPES_H */
