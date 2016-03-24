//
// Copyright (C) 2016 Karl Wette
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307 USA
//

#ifndef _FITSFILEIO_H
#define _FITSFILEIO_H

#include <stdint.h>
#include <gsl/gsl_matrix.h>
#include <lal/LALStdlib.h>
#include <lal/LALVCSInfoType.h>

#ifdef __cplusplus
extern "C" {
#endif

///
/// \defgroup FITSFileIO_h Header FITSFileIO.h
/// \ingroup lalpulsar_general
/// \author Karl Wette
/// \brief Routines for reading/writing data to/from FITS files
///

/// @{

///
/// \name General Definitions
///
/// @{

///
/// Maximum possible number of FITS array dimensions, and FITS table columns
///
#define FFIO_MAX 99

///
/// Representation of a FITS file
///
typedef struct tagFITSFile FITSFile;

/// @}

///
/// \name Open/Close FITS File
///
/// These functions open a FITS file \p file_name for writing or reading, and close a FITS file
/// represented by \p file. Writing a FITS file required a list of LALVCSInfo from the appropriate
/// LAL library VCS header, e.g. #lalVCSInfoList from \ref LALVCSInfo_h
///
/// @{
void XLALFITSFileClose(FITSFile *file);
FITSFile *XLALFITSFileOpenWrite(const CHAR *file_name, const LALVCSInfo *const vcs_list[]);
FITSFile *XLALFITSFileOpenRead(const CHAR *file_name);
/// @}

///
/// \name Write/Read Keyword-Value Pairs To/From FITS Header Data Unit
///
/// These functions write/read keyword-value pairs (\p keyword, \p value) to/from a FITS Header Data
/// Unit (HDU).  Scalar #BOOLEAN, #INT4, #INT8, #REAL4, #REAL8, #COMPLEX8, and #COMPLEX16 values,
/// strings, string vectors (#LALStringVector) and GPS times (#LIGOTimeGPS) can be written and read.
/// A \p comment string describing the value is required when writing to an HDU, and an arbitrary
/// comment can also be written.
///
/// There are some usage restrictions:
///
/// - When writing to or reading from the primary HDU, these functions can be used as soon as the
///   FITS file is opened.
///
/// - When writing to an extension HDU containing an array or table, these functions can only be
///   used after data has been written to the array or table.
///
/// - When reading from an extension HDU containing an array or table, these functions can only be
///   used after the array or table is opened for reading.
///
/// @{
int XLALFITSHeaderWriteComment(FITSFile *file, const CHAR *comment);
int XLALFITSHeaderWriteBOOLEAN(FITSFile *file, const CHAR *keyword, const BOOLEAN value, const CHAR *comment);
int XLALFITSHeaderReadBOOLEAN(FITSFile *file, const CHAR *keyword, BOOLEAN *value);
int XLALFITSHeaderWriteINT4(FITSFile *file, const CHAR *keyword, const INT4 value, const CHAR *comment);
int XLALFITSHeaderReadINT4(FITSFile *file, const CHAR *keyword, INT4 *value);
int XLALFITSHeaderWriteINT8(FITSFile *file, const CHAR *keyword, const INT8 value, const CHAR *comment);
int XLALFITSHeaderReadINT8(FITSFile *file, const CHAR *keyword, INT8 *value);
int XLALFITSHeaderWriteREAL4(FITSFile *file, const CHAR *keyword, const REAL4 value, const CHAR *comment);
int XLALFITSHeaderReadREAL4(FITSFile *file, const CHAR *keyword, REAL4 *value);
int XLALFITSHeaderWriteREAL8(FITSFile *file, const CHAR *keyword, const REAL8 value, const CHAR *comment);
int XLALFITSHeaderReadREAL8(FITSFile *file, const CHAR *keyword, REAL8 *value);
int XLALFITSHeaderWriteCOMPLEX8(FITSFile *file, const CHAR *keyword, const COMPLEX8 value, const CHAR *comment);
int XLALFITSHeaderReadCOMPLEX8(FITSFile *file, const CHAR *keyword, COMPLEX8 *value);
int XLALFITSHeaderWriteCOMPLEX16(FITSFile *file, const CHAR *keyword, const COMPLEX16 value, const CHAR *comment);
int XLALFITSHeaderReadCOMPLEX16(FITSFile *file, const CHAR *keyword, COMPLEX16 *value);
int XLALFITSHeaderWriteString(FITSFile *file, const CHAR *keyword, const CHAR *value, const CHAR *comment);
int XLALFITSHeaderReadString(FITSFile *file, const CHAR *keyword, CHAR **value);
int XLALFITSHeaderWriteStringVector(FITSFile *file, const CHAR *keyword, const LALStringVector *values, const CHAR *comment);
int XLALFITSHeaderReadStringVector(FITSFile *file, const CHAR *keyword, LALStringVector **values);
int XLALFITSHeaderWriteGPSTime(FITSFile *file, const CHAR *keyword, const LIGOTimeGPS *value, const CHAR *comment);
int XLALFITSHeaderReadGPSTime(FITSFile *file, const CHAR *keyword, LIGOTimeGPS *value);
/// @}

///
/// \name Write/Read Array To/From FITS File
///
/// These function write/read arbitrary-dimensional array to/from a FITS image extension HDU. A call
/// to XLALFITSArrayOpenWrite() or XLALFITSArrayOpenRead() is required first to write/read the
/// dimension count \p ndim and sizes \p dims[]; the functions
/// <tt>XLALFITSArrayOpenWrite<i>N</i>()</tt> and <tt>XLALFITSArrayOpenRead<i>N</i>()</tt> are
/// convenience functions for 1- and 2-dimensional arrays. The functions
/// <tt>XLALFITSArrayWrite<i>TYPE</i>()</tt> and <tt>XLALFITSArrayRead<i>TYPE</i>()</tt> are then
/// used to write scalar #UINT2, #UINT4, #INT2, #INT4, #REAL4, and #REAL8 values to the array
/// element indexed by \p idx. For the convenience the functions XLALFITSArrayWriteGSLMatrix() and
/// XLALFITSArrayReadGSLMatrix() exist for writing/reading elements to/from a <tt>gsl_matrix</tt>.
///
/// @{
int XLALFITSArrayOpenWrite(FITSFile *file, const CHAR *name, const size_t ndim, const size_t dims[], const CHAR *comment);
int XLALFITSArrayOpenRead(FITSFile *file, const CHAR *name, size_t *ndim, size_t dims[]);
int XLALFITSArrayOpenWrite1(FITSFile *file, const CHAR *name, const size_t dim0, const CHAR *comment);
int XLALFITSArrayOpenRead1(FITSFile *file, const CHAR *name, size_t *dim0);
int XLALFITSArrayOpenWrite2(FITSFile *file, const CHAR *name, const size_t dim0, const size_t dim1, const CHAR *comment);
int XLALFITSArrayOpenRead2(FITSFile *file, const CHAR *name, size_t *dim0, size_t *dim1);
int XLALFITSArrayWriteUINT2(FITSFile *file, const size_t idx[], const UINT2 elem);
int XLALFITSArrayReadUINT2(FITSFile *file, const size_t idx[], UINT2 *elem);
int XLALFITSArrayWriteUINT4(FITSFile *file, const size_t idx[], const UINT4 elem);
int XLALFITSArrayReadUINT4(FITSFile *file, const size_t idx[], UINT4 *elem);
int XLALFITSArrayWriteINT2(FITSFile *file, const size_t idx[], const INT2 elem);
int XLALFITSArrayReadINT2(FITSFile *file, const size_t idx[], INT2 *elem);
int XLALFITSArrayWriteINT4(FITSFile *file, const size_t idx[], const INT4 elem);
int XLALFITSArrayReadINT4(FITSFile *file, const size_t idx[], INT4 *elem);
int XLALFITSArrayWriteREAL4(FITSFile *file, const size_t idx[], const REAL4 elem);
int XLALFITSArrayReadREAL4(FITSFile *file, const size_t idx[], REAL4 *elem);
int XLALFITSArrayWriteREAL8(FITSFile *file, const size_t idx[], const REAL8 elem);
int XLALFITSArrayReadREAL8(FITSFile *file, const size_t idx[], REAL8 *elem);
int XLALFITSArrayWriteGSLMatrix(FITSFile *file, const size_t idx[], const gsl_matrix *elems);
int XLALFITSArrayReadGSLMatrix(FITSFile *file, const size_t idx[], gsl_matrix **elems);
/// @}

///
/// \name Write/Read Table To/From FITS File
///
/// These functions write/read arbitrary tables to/from a FITS table extension HDU. A call to
/// XLALFITSTableOpenWrite() or XLALFITSTableOpenRead() is required first; the latter returns the
/// number of rows in the table \p nrows, if needed. The table columns must then be specified using
/// the <tt>XLAL_FITS_TABLE_COLUMN_...()</tt> macros:
///
/// First, \ref XLAL_FITS_TABLE_COLUMN_BEGIN(\p record_type) is called, where \p record_type is the
/// name of a C structure that will store the values of the columns written to or read from the
/// table. For example:
///
/// \code
/// typedef struct { INT4 index; CHAR name[32]; REAL4 x; } Row
/// XLAL_FITS_TABLE_COLUMN_BEGIN(Row);
/// \endcode
///
/// Next, \ref XLAL_FITS_TABLE_COLUMN_ADD(\p file, \p type, \p field) is called, where \p file is
/// the FITS file, \p type is the C type of the field in \p record_type, and \p field is the C
/// structure field. For array fields, including fixed-length strings, \ref
/// XLAL_FITS_TABLE_COLUMN_ADD_ARRAY(\p file, \p type, \p field) should be used instead. An
/// alternative name for the FITS table column can be specified using \ref
/// XLAL_FITS_TABLE_COLUMN_ADD_NAMED(\p file, \p type, \p field, \p col_name). For example:
///
/// \code
/// XLAL_FITS_TABLE_COLUMN_ADD(file, INT4, index);
/// XLAL_FITS_TABLE_COLUMN_ADD_ARRAY(file, CHAR, name);
/// XLAL_FITS_TABLE_COLUMN_ADD_NAMED(file, REAL4, x, "value");
/// \endcode
///
/// Finally, XLALFITSTableWriteRow() or XLALFITSTableReadRow() are called to write/read table rows;
/// the latter returns the number of rows remaining in the table \p rem_nrows, if needed.
///
/// @{
int XLALFITSTableOpenWrite(FITSFile *file, const CHAR *name, const CHAR *comment);
int XLALFITSTableOpenRead(FITSFile *file, const CHAR *name, UINT8 *nrows);
/// \cond DONT_DOXYGEN
int XLALFITSTableColumnAddBOOLEAN(FITSFile *file, const CHAR *col_name, const void *record, const size_t record_size, const BOOLEAN *field, const size_t field_size);
int XLALFITSTableColumnAddINT2(FITSFile *file, const CHAR *col_name, const void *record, const size_t record_size, const INT2 *field, const size_t field_size);
int XLALFITSTableColumnAddINT4(FITSFile *file, const CHAR *col_name, const void *record, const size_t record_size, const INT4 *field, const size_t field_size);
int XLALFITSTableColumnAddREAL4(FITSFile *file, const CHAR *col_name, const void *record, const size_t record_size, const REAL4 *field, const size_t field_size);
int XLALFITSTableColumnAddREAL8(FITSFile *file, const CHAR *col_name, const void *record, const size_t record_size, const REAL8 *field, const size_t field_size);
int XLALFITSTableColumnAddCOMPLEX8(FITSFile *file, const CHAR *col_name, const void *record, const size_t record_size, const COMPLEX8 *field, const size_t field_size);
int XLALFITSTableColumnAddCOMPLEX16(FITSFile *file, const CHAR *col_name, const void *record, const size_t record_size, const COMPLEX16 *field, const size_t field_size);
int XLALFITSTableColumnAddCHAR(FITSFile *file, const CHAR *col_name, const void *record, const size_t record_size, const void *CHAR, const size_t field_size);
int XLALFITSTableColumnAddGPSTime(FITSFile *file, const CHAR *col_name, const void *record, const size_t record_size, const LIGOTimeGPS *field, const size_t field_size);
/// \endcond
#define XLAL_FITS_TABLE_COLUMN_BEGIN(record_type) \
  const record_type _xlal_fits_table_record_
#define XLAL_FITS_TABLE_COLUMN_ADD(file, type, field) \
  XLALFITSTableColumnAdd ## type (file, #field, &_xlal_fits_table_record_, sizeof(_xlal_fits_table_record_), &(_xlal_fits_table_record_.field), sizeof(_xlal_fits_table_record_.field))
#define XLAL_FITS_TABLE_COLUMN_ADD_ARRAY(file, type, field) \
  XLALFITSTableColumnAdd ## type (file, #field, &_xlal_fits_table_record_, sizeof(_xlal_fits_table_record_), &(_xlal_fits_table_record_.field[0]), sizeof(_xlal_fits_table_record_.field))
#define XLAL_FITS_TABLE_COLUMN_ADD_NAMED(file, type, field, col_name) \
  XLALFITSTableColumnAdd ## type (file, col_name, &_xlal_fits_table_record_, sizeof(_xlal_fits_table_record_), &(_xlal_fits_table_record_.field), sizeof(_xlal_fits_table_record_.field))
int XLALFITSTableWriteRow(FITSFile *file, const void *record);
int XLALFITSTableReadRow(FITSFile *file, void *record, UINT8 *rem_nrows);
/// @}

/// @}

#ifdef __cplusplus
}
#endif

#endif // _FITSFILEIO_H
