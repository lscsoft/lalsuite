/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton, Robert Adam Mercer, Xavier Siemens
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

#ifndef _LALFRAMEIO_H
#define _LALFRAMEIO_H

#include <lal/LALDatatypes.h>
#include <lal/LALDetectors.h>
#include <lal/LALFrameU.h>

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

struct tagLALFrFile;
struct tagLALFrameUFrameH;

/**
 * @defgroup LALFrameIO_h Header LALFrameIO.h
 * @ingroup lalframe_general
 *
 * @author Jolien Creighton
 * @brief Provides an intermediate-level interface for working on individual frame-files.
 * @details
 * This provides an intermediate-level LAL interface for working on individual frame files.
 * @sa <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
/** @{ */

/**
 * @defgroup LALFrameIO_h_read Frame Reading Routines
 * @ingroup LALFrameIO_h
 * @brief Intermediate-level routines for reading frame files.
 * @details
 * These routines provide the intermediate-level routines for reading from
 * frame files and converting the contents to LAL data types.
 * @{
 */

/**
 * @brief Incomplete type for a frame file structure.
 * @details
 * This structure contains the information used for reading from a frame file.
 */
typedef struct tagLALFrFile LALFrFile;

/**
 * @name Frame File Open/Close/Checksum Routines
 * @{
 */

/**
 * @brief Close a frame file described by a LALFrFile structure.
 * @note This routine is a no-op if passed a NULL pointer.
 * @param frfile Pointer to the #LALFrFile structure
 * @retval 0 Successfully closed the file and freed memory.
 */
int XLALFrFileClose(LALFrFile * frfile);

/**
 * @brief Open frame file for reading and return a LALFrFile structure.
 * @note Only "file:" protocol is supported in URLs.
 * @param url URL of the frame file to be opened.
 * @return Pointer to a LALFrFile structure that can be used to read the frame
 * file, or NULL if the URL could not be opened.
 */
LALFrFile *XLALFrFileOpenURL(const char *url);

/**
 * @brief Use checksum to determine if a frame file is valid.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @return Logical value indicating if the frame file checksum is correct.
 * @retval 1 The frame file is valid.
 * @retval 0 The frame file has an incorrect checksum or an error occurred.
 */
int XLALFrFileCksumValid(LALFrFile * frfile);

/** @} */

/**
 * @name Routines to Query Frame File Contents
 * @{
 */

/**
 * @brief Query a frame file for the number of frames contained in the file.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @return Then number of frame structures contained.
 * @retval (size_t)(-1) Failure.
 */
size_t XLALFrFileQueryNFrame(const LALFrFile * frfile);

/** 
 * @brief Query a frame file for the start time of a particular frame.
 * @param[out] start Pointer to a #LIGOTimeGPS structure containing the start time.
 * @param[in] frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param[in] pos The index of the frame in the frame file.
 * @returns The pointer to the #LIGOTimeGPS parameter start, with its values
 * set to the start GPS time of the specified frame.
 */
LIGOTimeGPS *XLALFrFileQueryGTime(LIGOTimeGPS * start, const LALFrFile * frfile, size_t pos);

/** 
 * @brief Query a frame file for the duration of a particular frame.
 * @param[in] frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param[in] pos The index of the frame in the frame file.
 * @returns The duration of the frame file in seconds.
 * @retval #LAL_REAL8_FAIL_NAN Failure.
 */
double XLALFrFileQueryDt(const LALFrFile * frfile, size_t pos);

/** 
 * @brief Query a frame file for the data type of a channel in a frame.
 * @param[in] frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param[in] chname String containing the name of the channel.
 * @param[in] pos The index of the frame in the frame file.
 * @returns The #LALTYPECODE value of the data type of the channel.
 * @retval LAL_CHAR_TYPE_CODE Channel is an array of type char.
 * @retval LAL_I2_TYPE_CODE Channel is an array of type int16_t.
 * @retval LAL_I4_TYPE_CODE Channel is an array of type int32_t.
 * @retval LAL_I8_TYPE_CODE Channel is an array of type int64_t.
 * @retval LAL_UCHAR_TYPE_CODE Channel is an array of type unsigned char.
 * @retval LAL_U2_TYPE_CODE Channel is an array of type uint16_t.
 * @retval LAL_U4_TYPE_CODE Channel is an array of type uint32_t.
 * @retval LAL_U8_TYPE_CODE Channel is an array of type uint64_t.
 * @retval LAL_S_TYPE_CODE Channel is an array of type float.
 * @retval LAL_D_TYPE_CODE Channel is an array of type double.
 * @retval LAL_C_TYPE_CODE Channel is an array of type float complex.
 * @retval LAL_Z_TYPE_CODE Channel is an array of type double complex.
 * @retval -1 Failure.
 */
LALTYPECODE XLALFrFileQueryChanType(const LALFrFile * frfile, const char *chname, size_t pos);

/** 
 * @brief Query a frame file for the number of data points in a channel in a frame.
 * @param[in] frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param[in] chname String containing the name of the channel.
 * @param[in] pos The index of the frame in the frame file.
 * @returns The length of the data vector of the channel in the specified frame.
 * @retval (size_t)(-1) Failure.
 */
size_t XLALFrFileQueryChanVectorLength(const LALFrFile * frfile, const char *chname, size_t pos);

/** @} */

/**
 * @name Routines to Read Channel Metadata
 * @{
 */

/**
 * @brief Acquires metadata about a specified channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #INT2TimeSeries containing the correct
 * channel name, sample units, epoch, and sampling interval.  No data is actually
 * read by this routine, and the length of the data vector is set to 0.
 * @retval NULL Failure.
 */
INT2TimeSeries *XLALFrFileReadINT2TimeSeriesMetadata(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Acquires metadata about a specified channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #INT4TimeSeries containing the correct
 * channel name, sample units, epoch, and sampling interval.  No data is actually
 * read by this routine, and the length of the data vector is set to 0.
 * @retval NULL Failure.
 */
INT4TimeSeries *XLALFrFileReadINT4TimeSeriesMetadata(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Acquires metadata about a specified channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #INT8TimeSeries containing the correct
 * channel name, sample units, epoch, and sampling interval.  No data is actually
 * read by this routine, and the length of the data vector is set to 0.
 * @retval NULL Failure.
 */
INT8TimeSeries *XLALFrFileReadINT8TimeSeriesMetadata(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Acquires metadata about a specified channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #UINT2TimeSeries containing the correct
 * channel name, sample units, epoch, and sampling interval.  No data is actually
 * read by this routine, and the length of the data vector is set to 0.
 * @retval NULL Failure.
 */
UINT2TimeSeries *XLALFrFileReadUINT2TimeSeriesMetadata(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Acquires metadata about a specified channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #UINT4TimeSeries containing the correct
 * channel name, sample units, epoch, and sampling interval.  No data is actually
 * read by this routine, and the length of the data vector is set to 0.
 * @retval NULL Failure.
 */
UINT4TimeSeries *XLALFrFileReadUINT4TimeSeriesMetadata(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Acquires metadata about a specified channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #UINT8TimeSeries containing the correct
 * channel name, sample units, epoch, and sampling interval.  No data is actually
 * read by this routine, and the length of the data vector is set to 0.
 * @retval NULL Failure.
 */
UINT8TimeSeries *XLALFrFileReadUINT8TimeSeriesMetadata(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Acquires metadata about a specified channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #REAL4TimeSeries containing the correct
 * channel name, sample units, epoch, and sampling interval.  No data is actually
 * read by this routine, and the length of the data vector is set to 0.
 * @retval NULL Failure.
 */
REAL4TimeSeries *XLALFrFileReadREAL4TimeSeriesMetadata(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Acquires metadata about a specified channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #REAL8TimeSeries containing the correct
 * channel name, sample units, epoch, and sampling interval.  No data is actually
 * read by this routine, and the length of the data vector is set to 0.
 * @retval NULL Failure.
 */
REAL8TimeSeries *XLALFrFileReadREAL8TimeSeriesMetadata(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Acquires metadata about a specified channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #COMPLEX8TimeSeries containing the correct
 * channel name, sample units, epoch, and sampling interval.  No data is actually
 * read by this routine, and the length of the data vector is set to 0.
 * @retval NULL Failure.
 */
COMPLEX8TimeSeries *XLALFrFileReadCOMPLEX8TimeSeriesMetadata(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Acquires metadata about a specified channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #COMPLEX16TimeSeries containing the correct
 * channel name, sample units, epoch, and sampling interval.  No data is actually
 * read by this routine, and the length of the data vector is set to 0.
 * @retval NULL Failure.
 */
COMPLEX16TimeSeries *XLALFrFileReadCOMPLEX16TimeSeriesMetadata(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Acquires metadata about a specified channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #REAL4FrequencySeries containing the correct
 * channel name, sample units, epoch, and sampling interval.  No data is actually
 * read by this routine, and the length of the data vector is set to 0.
 * @retval NULL Failure.
 */
REAL4FrequencySeries *XLALFrFileReadREAL4FrequencySeriesMetadata(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Acquires metadata about a specified channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #REAL8FrequencySeries containing the correct
 * channel name, sample units, epoch, and sampling interval.  No data is actually
 * read by this routine, and the length of the data vector is set to 0.
 * @retval NULL Failure.
 */
REAL8FrequencySeries *XLALFrFileReadREAL8FrequencySeriesMetadata(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Acquires metadata about a specified channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #COMPLEX8FrequencySeries containing the correct
 * channel name, sample units, epoch, and sampling interval.  No data is actually
 * read by this routine, and the length of the data vector is set to 0.
 * @retval NULL Failure.
 */
COMPLEX8FrequencySeries *XLALFrFileReadCOMPLEX8FrequencySeriesMetadata(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Acquires metadata about a specified channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #COMPLEX16FrequencySeries containing the correct
 * channel name, sample units, epoch, and sampling interval.  No data is actually
 * read by this routine, and the length of the data vector is set to 0.
 * @retval NULL Failure.
 */
COMPLEX16FrequencySeries *XLALFrFileReadCOMPLEX16FrequencySeriesMetadata(LALFrFile * frfile, const char *chname, size_t pos);

/** @} */

/**
 * @name Routines to Read Channel Data
 * @{
 */

/**
 * @brief Reads data from a channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #INT2TimeSeries containing the data
 * from the specified channel in the specified frame.
 * @retval NULL Failure.
 */
INT2TimeSeries *XLALFrFileReadINT2TimeSeries(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Reads data from a channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #INT4TimeSeries containing the data
 * from the specified channel in the specified frame.
 * @retval NULL Failure.
 */
INT4TimeSeries *XLALFrFileReadINT4TimeSeries(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Reads data from a channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #INT8TimeSeries containing the data
 * from the specified channel in the specified frame.
 * @retval NULL Failure.
 */
INT8TimeSeries *XLALFrFileReadINT8TimeSeries(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Reads data from a channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #UINT2TimeSeries containing the data
 * from the specified channel in the specified frame.
 * @retval NULL Failure.
 */
UINT2TimeSeries *XLALFrFileReadUINT2TimeSeries(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Reads data from a channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #UINT4TimeSeries containing the data
 * from the specified channel in the specified frame.
 * @retval NULL Failure.
 */
UINT4TimeSeries *XLALFrFileReadUINT4TimeSeries(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Reads data from a channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #UINT8TimeSeries containing the data
 * from the specified channel in the specified frame.
 * @retval NULL Failure.
 */
UINT8TimeSeries *XLALFrFileReadUINT8TimeSeries(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Reads data from a channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #REAL4TimeSeries containing the data
 * from the specified channel in the specified frame.
 * @retval NULL Failure.
 */
REAL4TimeSeries *XLALFrFileReadREAL4TimeSeries(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Reads data from a channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #REAL8TimeSeries containing the data
 * from the specified channel in the specified frame.
 * @retval NULL Failure.
 */
REAL8TimeSeries *XLALFrFileReadREAL8TimeSeries(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Reads data from a channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #COMPLEX8TimeSeries containing the data
 * from the specified channel in the specified frame.
 * @retval NULL Failure.
 */
COMPLEX8TimeSeries *XLALFrFileReadCOMPLEX8TimeSeries(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Reads data from a channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #COMPLEX16TimeSeries containing the data
 * from the specified channel in the specified frame.
 * @retval NULL Failure.
 */
COMPLEX16TimeSeries *XLALFrFileReadCOMPLEX16TimeSeries(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Reads data from a channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #REAL4FrequencySeries containing the data
 * from the specified channel in the specified frame.
 * @retval NULL Failure.
 */
REAL4FrequencySeries *XLALFrFileReadREAL4FrequencySeries(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Reads data from a channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #REAL8FrequencySeries containing the data
 * from the specified channel in the specified frame.
 * @retval NULL Failure.
 */
REAL8FrequencySeries *XLALFrFileReadREAL8FrequencySeries(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Reads data from a channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #COMPLEX8FrequencySeries containing the data
 * from the specified channel in the specified frame.
 * @retval NULL Failure.
 */
COMPLEX8FrequencySeries *XLALFrFileReadCOMPLEX8FrequencySeries(LALFrFile * frfile, const char *chname, size_t pos);

/**
 * @brief Reads data from a channel in a frame.
 * @param frfile Pointer to a #LALFrFile structure associated with a frame file.
 * @param chname String containing the name of the channel.
 * @param pos The index of the frame in the frame file.
 * @returns A pointer to a newly allocated #COMPLEX8FrequencySeries containing the data
 * from the specified channel in the specified frame.
 * @retval NULL Failure.
 */
COMPLEX16FrequencySeries *XLALFrFileReadCOMPLEX16FrequencySeries(LALFrFile * frfile, const char *chname, size_t pos);

/** @} */

/** @} */

/**
 * @defgroup LALFrameIO_h_write Frame Writing Routines
 * @ingroup LALFrameIO_h
 * @brief Intermediate-level routines for writing frame files.
 * @details
 * These routines provide the intermediate-level routines for writing
 * LAL data types to frame files.
 * @{
 */

/**
 * @brief Incomplete type for a frame header structure.
 * @details
 * This structure contains information about an individual frame.  In this
 * interface, these frames are constructed for writing to frame files.
 * @sa Section 4.3.2.3 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
typedef struct tagLALFrameUFrameH LALFrameH;

/**
 * @name Basic Frame Writing Routines
 * @brief These routines can be used to directly output series data to a frame file.
 * @{
 */

/**
 * @brief Creates a frame file holding the data in a #INT2TimeSeries.
 * @details
 * Outputs the data contained in a #INT2TimeSeries to a frame file
 * whose name is derived from the series name.  The name of the file
 * conforms to the convention LIGO-T010150.  If the name begins with
 * @c XmYn:channel_name then this routine adds the appropriate detector
 * structures for prefixes @c Xm and @c Yn to the frame, and the name of
 * the file will be <tt>XmYn-description-\<start\>-\<duration\>.gwf</tt>.
 * The data is written as a FrProcData structure.
 * @param[in] series The series to write to the frame file.
 * @param[in] fnum The frame number of the frame to write.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa LIGO-T010150 [https://dcc.ligo.org/LIGO-T010150/public]
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrWriteINT2TimeSeries(const INT2TimeSeries * series, int fnum);

/**
 * @brief Creates a frame file holding the data in a #INT4TimeSeries.
 * @details
 * Outputs the data contained in a #INT4TimeSeries to a frame file
 * whose name is derived from the series name.  The name of the file
 * conforms to the convention LIGO-T010150.  If the name begins with
 * @c XmYn:description then this routine adds the appropriate detector
 * structures for prefixes @c Xm and @c Yn to the frame, and the name of
 * the file will be <tt>XmYn-description-\<start\>-\<duration\>.gwf</tt>.
 * The data is written as a FrProcData structure.
 * @param[in] series The series to write to the frame file.
 * @param[in] fnum The frame number of the frame to write.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa LIGO-T010150 [https://dcc.ligo.org/LIGO-T010150/public]
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrWriteINT4TimeSeries(const INT4TimeSeries * series, int fnum);

/**
 * @brief Creates a frame file holding the data in a #INT8TimeSeries.
 * @details
 * Outputs the data contained in a #INT8TimeSeries to a frame file
 * whose name is derived from the series name.  The name of the file
 * conforms to the convention LIGO-T010150.  If the name begins with
 * @c XmYn:description then this routine adds the appropriate detector
 * structures for prefixes @c Xm and @c Yn to the frame, and the name of
 * the file will be <tt>XmYn-description-\<start\>-\<duration\>.gwf</tt>.
 * The data is written as a FrProcData structure.
 * @param[in] series The series to write to the frame file.
 * @param[in] fnum The frame number of the frame to write.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa LIGO-T010150 [https://dcc.ligo.org/LIGO-T010150/public]
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrWriteINT8TimeSeries(const INT8TimeSeries * series, int fnum);

/**
 * @brief Creates a frame file holding the data in a #REAL4TimeSeries.
 * @details
 * Outputs the data contained in a #REAL4TimeSeries to a frame file
 * whose name is derived from the series name.  The name of the file
 * conforms to the convention LIGO-T010150.  If the name begins with
 * @c XmYn:description then this routine adds the appropriate detector
 * structures for prefixes @c Xm and @c Yn to the frame, and the name of
 * the file will be <tt>XmYn-description-\<start\>-\<duration\>.gwf</tt>.
 * The data is written as a FrProcData structure.
 * @param[in] series The series to write to the frame file.
 * @param[in] fnum The frame number of the frame to write.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa LIGO-T010150 [https://dcc.ligo.org/LIGO-T010150/public]
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrWriteREAL4TimeSeries(const REAL4TimeSeries * series, int fnum);

/**
 * @brief Creates a frame file holding the data in a #REAL8TimeSeries.
 * @details
 * Outputs the data contained in a #REAL8TimeSeries to a frame file
 * whose name is derived from the series name.  The name of the file
 * conforms to the convention LIGO-T010150.  If the name begins with
 * @c XmYn:description then this routine adds the appropriate detector
 * structures for prefixes @c Xm and @c Yn to the frame, and the name of
 * the file will be <tt>XmYn-description-\<start\>-\<duration\>.gwf</tt>.
 * The data is written as a FrProcData structure.
 * @param[in] series The series to write to the frame file.
 * @param[in] fnum The frame number of the frame to write.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa LIGO-T010150 [https://dcc.ligo.org/LIGO-T010150/public]
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrWriteREAL8TimeSeries(const REAL8TimeSeries * series, int fnum);

/**
 * @brief Creates a frame file holding the data in a #COMPLEX8TimeSeries.
 * @details
 * Outputs the data contained in a #COMPLEX8TimeSeries to a frame file
 * whose name is derived from the series name.  The name of the file
 * conforms to the convention LIGO-T010150.  If the name begins with
 * @c XmYn:description then this routine adds the appropriate detector
 * structures for prefixes @c Xm and @c Yn to the frame, and the name of
 * the file will be <tt>XmYn-description-\<start\>-\<duration\>.gwf</tt>.
 * The data is written as a FrProcData structure.
 * @param[in] series The series to write to the frame file.
 * @param[in] fnum The frame number of the frame to write.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa LIGO-T010150 [https://dcc.ligo.org/LIGO-T010150/public]
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrWriteCOMPLEX8TimeSeries(const COMPLEX8TimeSeries * series, int fnum);

/**
 * @brief Creates a frame file holding the data in a #COMPLEX16TimeSeries.
 * @details
 * Outputs the data contained in a #COMPLEX16TimeSeries to a frame file
 * whose name is derived from the series name.  The name of the file
 * conforms to the convention LIGO-T010150.  If the name begins with
 * @c XmYn:description then this routine adds the appropriate detector
 * structures for prefixes @c Xm and @c Yn to the frame, and the name of
 * the file will be <tt>XmYn-description-\<start\>-\<duration\>.gwf</tt>.
 * The data is written as a FrProcData structure.
 * @param[in] series The series to write to the frame file.
 * @param[in] fnum The frame number of the frame to write.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa LIGO-T010150 [https://dcc.ligo.org/LIGO-T010150/public]
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrWriteCOMPLEX16TimeSeries(const COMPLEX16TimeSeries * series, int fnum);

/**
 * @brief Creates a frame file holding the data in a #REAL4FrequencySeries.
 * @details
 * Outputs the data contained in a #REAL4FrequencySeries to a frame file
 * whose name is derived from the series name.  The name of the file
 * conforms to the convention LIGO-T010150.  If the name begins with
 * @c XmYn:description then this routine adds the appropriate detector
 * structures for prefixes @c Xm and @c Yn to the frame, and the name of
 * the file will be <tt>XmYn-description-\<start\>-\<duration\>.gwf</tt>.
 * The data is written as a FrProcData structure.
 *
 * The subtypes are as follows:
 *
 * Value  | Subtype
 * :----: | :-------------------------
 *   0    | Unknown / User-Defined
 *   1    | Discrete Fourier Transform
 *   2    | Amplitude Spectral Density
 *   3    | Power Spectral Density
 *   4    | Cross Spectral Density
 *   5    | Coherence
 *   6    | Transfer Function
 *
 * @param[in] series The series to write to the frame file.
 * @param[in] fnum The frame number of the frame to write.
 * @param[in] subtype The FrProcData subtype of this data.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa LIGO-T010150 [https://dcc.ligo.org/LIGO-T010150/public]
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrWriteREAL4FrequencySeries(const REAL4FrequencySeries * series, int fnum, int subtype);

/**
 * @brief Creates a frame file holding the data in a #REAL8FrequencySeries.
 * @details
 * Outputs the data contained in a #REAL8FrequencySeries to a frame file
 * whose name is derived from the series name.  The name of the file
 * conforms to the convention LIGO-T010150.  If the name begins with
 * @c XmYn:description then this routine adds the appropriate detector
 * structures for prefixes @c Xm and @c Yn to the frame, and the name of
 * the file will be <tt>XmYn-description-\<start\>-\<duration\>.gwf</tt>.
 * The data is written as a FrProcData structure.
 *
 * The subtypes are as follows:
 *
 * Value  | Subtype
 * :----: | :-------------------------
 *   0    | Unknown / User-Defined
 *   1    | Discrete Fourier Transform
 *   2    | Amplitude Spectral Density
 *   3    | Power Spectral Density
 *   4    | Cross Spectral Density
 *   5    | Coherence
 *   6    | Transfer Function
 *
 * @param[in] series The series to write to the frame file.
 * @param[in] fnum The frame number of the frame to write.
 * @param[in] subtype The FrProcData subtype of this data.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa LIGO-T010150 [https://dcc.ligo.org/LIGO-T010150/public]
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrWriteREAL8FrequencySeries(const REAL8FrequencySeries * series, int fnum, int subtype);

/**
 * @brief Creates a frame file holding the data in a #COMPLEX8FrequencySeries.
 * @details
 * Outputs the data contained in a #COMPLEX8FrequencySeries to a frame file
 * whose name is derived from the series name.  The name of the file
 * conforms to the convention LIGO-T010150.  If the name begins with
 * @c XmYn:description then this routine adds the appropriate detector
 * structures for prefixes @c Xm and @c Yn to the frame, and the name of
 * the file will be <tt>XmYn-description-\<start\>-\<duration\>.gwf</tt>.
 * The data is written as a FrProcData structure.
 *
 * The subtypes are as follows:
 *
 * Value  | Subtype
 * :----: | :-------------------------
 *   0    | Unknown / User-Defined
 *   1    | Discrete Fourier Transform
 *   2    | Amplitude Spectral Density
 *   3    | Power Spectral Density
 *   4    | Cross Spectral Density
 *   5    | Coherence
 *   6    | Transfer Function
 *
 * @param[in] series The series to write to the frame file.
 * @param[in] fnum The frame number of the frame to write.
 * @param[in] subtype The FrProcData subtype of this data.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa LIGO-T010150 [https://dcc.ligo.org/LIGO-T010150/public]
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrWriteCOMPLEX8FrequencySeries(const COMPLEX8FrequencySeries * series, int fnum, int subtype);

/**
 * @brief Creates a frame file holding the data in a #COMPLEX16FrequencySeries.
 * @details
 * Outputs the data contained in a #COMPLEX16FrequencySeries to a frame file
 * whose name is derived from the series name.  The name of the file
 * conforms to the convention LIGO-T010150.  If the name begins with
 * @c XmYn:description then this routine adds the appropriate detector
 * structures for prefixes @c Xm and @c Yn to the frame, and the name of
 * the file will be <tt>XmYn-description-\<start\>-\<duration\>.gwf</tt>.
 * The data is written as a FrProcData structure.
 *
 * The subtypes are as follows:
 *
 * Value  | Subtype
 * :----: | :-------------------------
 *   0    | Unknown / User-Defined
 *   1    | Discrete Fourier Transform
 *   2    | Amplitude Spectral Density
 *   3    | Power Spectral Density
 *   4    | Cross Spectral Density
 *   5    | Coherence
 *   6    | Transfer Function
 *
 * @param[in] series The series to write to the frame file.
 * @param[in] fnum The frame number of the frame to write.
 * @param[in] subtype The FrProcData subtype of this data.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa LIGO-T010150 [https://dcc.ligo.org/LIGO-T010150/public]
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrWriteCOMPLEX16FrequencySeries(const COMPLEX16FrequencySeries * series, int fnum, int subtype);

/* @} */

/**
 * @name Advanced Frame Writing Routines
 * @brief Routines allowing for more control of the data written to a frame file.
 * @details
 * These routines allow the user to add specific components to a #LALFrameH structure,
 * which represents a single frame.  Various channels of and other data types can be
 * added to the #LALFrameH structure, and then the frame can be written to a frame file.
 * @{
 */

/**
 * @brief Frees a frame structure.
 * @note This routine is a no-op if passed a NULL pointer.
 * @param frame Pointer to a #LALFrameH structure.
 */
void XLALFrameFree(LALFrameH * frame);

/**
 * @brief Creates a new frame structure.
 * @details
 * Creates a new frame structure with a specified start time, duration, project
 * name (i.e., "LIGO"), run serial number, frame serial number within the run,
 * and detectors associated with the data that will be added to the frame.  The
 * detectors are specified in a flag field composed of the various detector bits,
 * e.g., <tt>( LAL_LHO_4K_DETECTOR_BIT | LAL_LLO_4K_DETECTOR_BIT )</tt> would
 * attach detector structures to the frame for the LIGO Hanford and Livingston
 * observatories.  See @ref LALDetectors_h.
 * @param[in] epoch Pointer to a #LIGOTimeGPS structure containing the start time of this frame.
 * @param[in] duration Duration of this frame in seconds.
 * @param[in] project String describing the project associated with this frame.
 * @param[in] run Run number associated with this frame.
 * @param[in] frnum Frame number (within the run) associated with this frame.
 * @param[in] detectorFlags Flag field specifying the detectors to be associated with the frame.
 * @returns Pointer to a new frame structure, or NULL if failure.
 * @sa Section 4.3.2.3 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
LALFrameH *XLALFrameNew(const LIGOTimeGPS * epoch, double duration, const char *project, int run, int frnum, int detectorFlags);

/**
 * @brief Adds a history structure to a frame.
 * @remark
 * History records can be added to a frame, e.g., to detail how the data was
 * generated.  For example, the program name, version information, etc., can
 * be added to a frame.  The name of the history record should match the name
 * of the FrProcData channel that it is associated with.
 * @param frame Pointer to a #LALFrameH frame structure to which the history will be added.
 * @param name String containing the name of the history record that will be added.
 * @param comment The history comment that will be added.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa Section 4.3.2.9 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameAddFrHistory(LALFrameH * frame, const char *name, const char *comment);

/**
 * @brief Adds a detector structure to a frame.
 * @param frame Pointer to a #LALFrameH frame structure to which the detector will be added.
 * @param detector Pointer to a #LALFrDetector structure to add to the frame.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa Section 4.3.2.5 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameAddFrDetector(LALFrameH * frame, const LALFrDetector * detector);

/**
 * @brief Adds an #INT2TimeSeries to a frame as a FrAdcData channel.
 * @remark FrAdcData channels contains "raw" interferometer data.
 * @param frame Pointer to a #LALFrameH frame structure to which the series will be added.
 * @param series Pointer to the series to add to the frame.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa Section 4.3.2.4 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameAddINT2TimeSeriesAdcData(LALFrameH * frame, const INT2TimeSeries * series);

/**
 * @brief Adds an #INT4TimeSeries to a frame as a FrAdcData channel.
 * @remark FrAdcData channels contains "raw" interferometer data.
 * @param frame Pointer to a #LALFrameH frame structure to which the series will be added.
 * @param series Pointer to the series to add to the frame.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa Section 4.3.2.4 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameAddINT4TimeSeriesAdcData(LALFrameH * frame, const INT4TimeSeries * series);

/**
 * @brief Adds an #REAL4TimeSeries to a frame as a FrAdcData channel.
 * @remark FrAdcData channels contains "raw" interferometer data.
 * @param frame Pointer to a #LALFrameH frame structure to which the series will be added.
 * @param series Pointer to the series to add to the frame.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa Section 4.3.2.4 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameAddREAL4TimeSeriesAdcData(LALFrameH * frame, const REAL4TimeSeries * series);

/**
 * @brief Adds an #REAL8TimeSeries to a frame as a FrAdcData channel.
 * @remark FrAdcData channels contains "raw" interferometer data.
 * @param frame Pointer to a #LALFrameH frame structure to which the series will be added.
 * @param series Pointer to the series to add to the frame.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa Section 4.3.2.4 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameAddREAL8TimeSeriesAdcData(LALFrameH * frame, const REAL8TimeSeries * series);

/**
 * @brief Adds an #INT2TimeSeries to a frame as a FrSimData channel.
 * @remark FrSimData channels contains simulated interferometer data.
 * @param frame Pointer to a #LALFrameH frame structure to which the series will be added.
 * @param series Pointer to the series to add to the frame.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa Section 4.3.2.14 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameAddINT2TimeSeriesSimData(LALFrameH * frame, const INT2TimeSeries * series);

/**
 * @brief Adds an #INT4TimeSeries to a frame as a FrSimData channel.
 * @remark FrSimData channels contains simulated interferometer data.
 * @param frame Pointer to a #LALFrameH frame structure to which the series will be added.
 * @param series Pointer to the series to add to the frame.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa Section 4.3.2.14 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameAddINT4TimeSeriesSimData(LALFrameH * frame, const INT4TimeSeries * series);

/**
 * @brief Adds an #REAL4TimeSeries to a frame as a FrSimData channel.
 * @remark FrSimData channels contains simulated interferometer data.
 * @param frame Pointer to a #LALFrameH frame structure to which the series will be added.
 * @param series Pointer to the series to add to the frame.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa Section 4.3.2.14 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameAddREAL4TimeSeriesSimData(LALFrameH * frame, const REAL4TimeSeries * series);

/**
 * @brief Adds an #REAL8TimeSeries to a frame as a FrSimData channel.
 * @remark FrSimData channels contains simulated interferometer data.
 * @param frame Pointer to a #LALFrameH frame structure to which the series will be added.
 * @param series Pointer to the series to add to the frame.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa Section 4.3.2.14 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameAddREAL8TimeSeriesSimData(LALFrameH * frame, const REAL8TimeSeries * series);

/**
 * @brief Adds an #INT2TimeSeries to a frame as a FrProcData channel.
 * @remark FrProcData channels contains post-processed data.
 * @param frame Pointer to a #LALFrameH frame structure to which the series will be added.
 * @param series Pointer to the series to add to the frame.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameAddINT2TimeSeriesProcData(LALFrameH * frame, const INT2TimeSeries * series);

/**
 * @brief Adds an #INT4TimeSeries to a frame as a FrProcData channel.
 * @remark FrProcData channels contains post-processed data.
 * @param frame Pointer to a #LALFrameH frame structure to which the series will be added.
 * @param series Pointer to the series to add to the frame.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameAddINT4TimeSeriesProcData(LALFrameH * frame, const INT4TimeSeries * series);

/**
 * @brief Adds an #INT8TimeSeries to a frame as a FrProcData channel.
 * @remark FrProcData channels contains post-processed data.
 * @param frame Pointer to a #LALFrameH frame structure to which the series will be added.
 * @param series Pointer to the series to add to the frame.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameAddINT8TimeSeriesProcData(LALFrameH * frame, const INT8TimeSeries * series);

/**
 * @brief Adds an #UINT2TimeSeries to a frame as a FrProcData channel.
 * @remark FrProcData channels contains post-processed data.
 * @param frame Pointer to a #LALFrameH frame structure to which the series will be added.
 * @param series Pointer to the series to add to the frame.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameAddUINT2TimeSeriesProcData(LALFrameH * frame, const UINT2TimeSeries * series);

/**
 * @brief Adds an #UINT4TimeSeries to a frame as a FrProcData channel.
 * @remark FrProcData channels contains post-processed data.
 * @param frame Pointer to a #LALFrameH frame structure to which the series will be added.
 * @param series Pointer to the series to add to the frame.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameAddUINT4TimeSeriesProcData(LALFrameH * frame, const UINT4TimeSeries * series);

/**
 * @brief Adds an #UINT8TimeSeries to a frame as a FrProcData channel.
 * @remark FrProcData channels contains post-processed data.
 * @param frame Pointer to a #LALFrameH frame structure to which the series will be added.
 * @param series Pointer to the series to add to the frame.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameAddUINT8TimeSeriesProcData(LALFrameH * frame, const UINT8TimeSeries * series);

/**
 * @brief Adds an #REAL4TimeSeries to a frame as a FrProcData channel.
 * @remark FrProcData channels contains post-processed data.
 * @param frame Pointer to a #LALFrameH frame structure to which the series will be added.
 * @param series Pointer to the series to add to the frame.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameAddREAL4TimeSeriesProcData(LALFrameH * frame, const REAL4TimeSeries * series);

/**
 * @brief Adds an #REAL8TimeSeries to a frame as a FrProcData channel.
 * @remark FrProcData channels contains post-processed data.
 * @param frame Pointer to a #LALFrameH frame structure to which the series will be added.
 * @param series Pointer to the series to add to the frame.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameAddREAL8TimeSeriesProcData(LALFrameH * frame, const REAL8TimeSeries * series);

/**
 * @brief Adds an #COMPLEX8TimeSeries to a frame as a FrProcData channel.
 * @remark FrProcData channels contains post-processed data.
 * @param frame Pointer to a #LALFrameH frame structure to which the series will be added.
 * @param series Pointer to the series to add to the frame.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameAddCOMPLEX8TimeSeriesProcData(LALFrameH * frame, const COMPLEX8TimeSeries * series);

/**
 * @brief Adds an #COMPLEX16TimeSeries to a frame as a FrProcData channel.
 * @remark FrProcData channels contains post-processed data.
 * @param frame Pointer to a #LALFrameH frame structure to which the series will be added.
 * @param series Pointer to the series to add to the frame.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameAddCOMPLEX16TimeSeriesProcData(LALFrameH * frame, const COMPLEX16TimeSeries * series);

/**
 * @brief Adds an #REAL4FrequencySeries to a frame as a FrProcData channel.
 * @remark FrProcData channels contains post-processed data.
 *
 * The subtypes are as follows:
 *
 * Value  | Subtype
 * :----: | :-------------------------
 *   0    | Unknown / User-Defined
 *   1    | Discrete Fourier Transform
 *   2    | Amplitude Spectral Density
 *   3    | Power Spectral Density
 *   4    | Cross Spectral Density
 *   5    | Coherence
 *   6    | Transfer Function
 *
 * @param frame Pointer to a #LALFrameH frame structure to which the series will be added.
 * @param series Pointer to the series to add to the frame.
 * @param subtype The FrProcData subtype of this frequency series.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameAddREAL4FrequencySeriesProcData(LALFrameH * frame, const REAL4FrequencySeries * series, int subtype);

/**
 * @brief Adds an #REAL8FrequencySeries to a frame as a FrProcData channel.
 * @remark FrProcData channels contains post-processed data.
 *
 * The subtypes are as follows:
 *
 * Value  | Subtype
 * :----: | :-------------------------
 *   0    | Unknown / User-Defined
 *   1    | Discrete Fourier Transform
 *   2    | Amplitude Spectral Density
 *   3    | Power Spectral Density
 *   4    | Cross Spectral Density
 *   5    | Coherence
 *   6    | Transfer Function
 *
 * @param frame Pointer to a #LALFrameH frame structure to which the series will be added.
 * @param series Pointer to the series to add to the frame.
 * @param subtype The FrProcData subtype of this frequency series.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameAddREAL8FrequencySeriesProcData(LALFrameH * frame, const REAL8FrequencySeries * series, int subtype);

/**
 * @brief Adds an #COMPLEX8FrequencySeries to a frame as a FrProcData channel.
 * @remark FrProcData channels contains post-processed data.
 *
 * The subtypes are as follows:
 *
 * Value  | Subtype
 * :----: | :-------------------------
 *   0    | Unknown / User-Defined
 *   1    | Discrete Fourier Transform
 *   2    | Amplitude Spectral Density
 *   3    | Power Spectral Density
 *   4    | Cross Spectral Density
 *   5    | Coherence
 *   6    | Transfer Function
 *
 * @param frame Pointer to a #LALFrameH frame structure to which the series will be added.
 * @param series Pointer to the series to add to the frame.
 * @param subtype The FrProcData subtype of this frequency series.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameAddCOMPLEX8FrequencySeriesProcData(LALFrameH * frame, const COMPLEX8FrequencySeries * series, int subtype);

/**
 * @brief Adds an #COMPLEX16FrequencySeries to a frame as a FrProcData channel.
 * @remark FrProcData channels contains post-processed data.
 *
 * The subtypes are as follows:
 *
 * Value  | Subtype
 * :----: | :-------------------------
 *   0    | Unknown / User-Defined
 *   1    | Discrete Fourier Transform
 *   2    | Amplitude Spectral Density
 *   3    | Power Spectral Density
 *   4    | Cross Spectral Density
 *   5    | Coherence
 *   6    | Transfer Function
 *
 * @param frame Pointer to a #LALFrameH frame structure to which the series will be added.
 * @param series Pointer to the series to add to the frame.
 * @param subtype The FrProcData subtype of this frequency series.
 * @retval 0 Success.
 * @retval -1 Failure.
 * @sa Section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameAddCOMPLEX16FrequencySeriesProcData(LALFrameH * frame, const COMPLEX16FrequencySeries * series, int subtype);


/**
 * @brief Write a #LALFrameH frame structure to a frame file.
 * @param frame Pointer to the #LALFrameH frame structure to be written.
 * @param fname String with the path name of the frame file to create.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALFrameWrite(LALFrameH * frame, const char *fname);

/** @} */

/** @} */

/** @} */

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
