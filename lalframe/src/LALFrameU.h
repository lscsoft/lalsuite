/*
*  Copyright (C) 2013 Jolien Creighton
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

/**
 * @file
 * @author Jolien Creighton
 * @brief Provides a unified interface to frame libraries.
 * @details
 * This provides a unified interface for reading and writing data files
 * in the "Frame Format for Interferometric Gravitational Wave Detectors".
 * @sa <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 * @sa Virgo Frame Library [http://lappweb.in2p3.fr/virgo/FrameL].
 * @sa FrameCPP Library [http://www.ldas-sw.ligo.caltech.edu/doc/framecpp/html].
 */

#include <stddef.h>
#include <time.h>

#ifndef _LALFRAMEU_H
#define _LALFRAMEU_H

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif
struct tagLALFrameUFrameH;
struct tagLALFrameUFrFile;
struct tagLALFrameUFrTOC;
struct tagLALFrameUFrChan;
struct tagLALFrameUFrDetector;
struct tagLALFrameUFrHistory;

/**
 * @brief Incomplete type for a frame header FrameH structure.
 * @sa Section 4.3.2.3 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
typedef struct tagLALFrameUFrameH LALFrameUFrameH;

/**
 * @brief Incomplete type for a frame file FrFile stream.
 */
typedef struct tagLALFrameUFrFile LALFrameUFrFile;

/**
 * @brief Incomplete type for a table of contents FrTOC structure.
 * @sa Section 4.3.2.19 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
typedef struct tagLALFrameUFrTOC LALFrameUFrTOC;

/**
 * @brief Incomplete type for a data channel FrChan structure.
 * @sa Sections 4.3.2.4, 4.3.2.11, and 4.3.2.14 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
typedef struct tagLALFrameUFrChan LALFrameUFrChan;

/**
 * @brief Incomplete type for a detector data FrDetector structure.
 * @sa Section 4.3.2.5 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
typedef struct tagLALFrameUFrDetector LALFrameUFrDetector;

/**
 * @brief Incomplete type for a history data FrHistory structure.
 * @sa Section 4.3.2.9 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
typedef struct tagLALFrameUFrHistory LALFrameUFrHistory;

/**
 * @brief Compression scheme codes.
 * @details These compression scheme Id codes are from Appendix B of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
enum LALFrameUFrVectCompressionScheme {
    /** @brief Id for uncompressed raw values. */
    LAL_FRAMEU_FR_VECT_COMPRESS_RAW = 0,

    /** @brief Id for gzip compression. */
    LAL_FRAMEU_FR_VECT_COMPRESS_GZIP = 1,

    /** @brief Id for gzip compression of differential values. */
    LAL_FRAMEU_FR_VECT_COMPRESS_DIFF_GZIP = 3,

    /** @brief Id for differentiation and zero suppression for 2-byte word
     * integer types only. */
    LAL_FRAMEU_FR_VECT_COMPRESS_ZERO_SUPPRESS_WORD_2 = 5,

    /** @brief Id for differentiation and zero suppression for 4-byte word
     * (integer or floating-point). */
    LAL_FRAMEU_FR_VECT_COMPRESS_ZERO_SUPPRESS_WORD_4 = 8
};

/**
 * @brief FrVect data type codes.
 * @details These data type Id codes are from Appendix C of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
enum LALFrameUFrVectDataType {
    /** @brief Id for 8-bit signed char data type. */
    LAL_FRAMEU_FR_VECT_C = 0,

    /** @brief Id for 16-bit signed integer data type. */
    LAL_FRAMEU_FR_VECT_2S = 1,

    /** @brief Id for 64-bit double precision floating point data type. */
    LAL_FRAMEU_FR_VECT_8R = 2,

    /** @brief Id for 32-bit single precision floating point data type. */
    LAL_FRAMEU_FR_VECT_4R = 3,

    /** @brief Id for 32-bit signed integer data type. */
    LAL_FRAMEU_FR_VECT_4S = 4,

    /** @brief Id for 64-bit signed integer data type. */
    LAL_FRAMEU_FR_VECT_8S = 5,

    /** @brief Id for 64-bit single precision complex data type. */
    LAL_FRAMEU_FR_VECT_8C = 6,

    /** @brief Id for 128-bit double precision complex data type. */
    LAL_FRAMEU_FR_VECT_16C = 7,

    /** @brief Id for string data type. */
    LAL_FRAMEU_FR_VECT_STRING = 8,

    /** @brief Id for 16-bit unsigned integer data type. */
    LAL_FRAMEU_FR_VECT_2U = 9,

    /** @brief Id for 32-bit unsigned integer data type. */
    LAL_FRAMEU_FR_VECT_4U = 10,

    /** @brief Id for 64-bit unsigned integer data type. */
    LAL_FRAMEU_FR_VECT_8U = 11,

    /** @brief Id for 8-bit unsigned char data type. */
    LAL_FRAMEU_FR_VECT_1U = 12
};

/**
 * @brief FrProcData type codes.
 * @details These type codes are from section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
enum LALFrameUFrProcDataType {
    /** @brief Id for unknown or user-defined data. */
    LAL_FRAMEU_FR_PROC_TYPE_UNKNOWN = 0,

    /** @brief Id for time-series data. */
    LAL_FRAMEU_FR_PROC_TYPE_TIME_SERIES = 1,

    /** @brief Id for frequency-series data. */
    LAL_FRAMEU_FR_PROC_TYPE_FREQUENCY_SERIES = 2,

    /** @brief Id for other one-dimensional series data. */
    LAL_FRAMEU_FR_PROC_TYPE_OTHER_1D_SERIES = 3,

    /** @brief Id for time-frequency data. */
    LAL_FRAMEU_FR_PROC_TYPE_TIME_FREQUENCY = 4,

    /** @brief Id for wavelet data. */
    LAL_FRAMEU_FR_PROC_TYPE_WAVELET = 5,

    /** @brief Id for multi-dimensional data. */
    LAL_FRAMEU_FR_PROC_TYPE_MULTI_DIMENSIONAL = 6
};

/**
 * @brief FrProcData subtype codes for frequency series.
 * @details These type codes are from section 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
enum LALFrameUFrProcDataSubType {
    /** @brief Id for unknown or user-defined data. */
    LAL_FRAMEU_FR_PROC_SUB_TYPE_UNKNOWN = 0,

    /** @brief Id for DFT data. */
    LAL_FRAMEU_FR_PROC_SUB_TYPE_DFT = 1,

    /** @brief Id for amplitude spectral density data. */
    LAL_FRAMEU_FR_PROC_SUB_TYPE_AMPLITUDE_SPECTRAL_DENSITY = 2,

    /** @brief Id for power spectral density data. */
    LAL_FRAMEU_FR_PROC_SUB_TYPE_POWER_SPECTRAL_DENSITY = 3,

    /** @brief Id for cross spectral density data. */
    LAL_FRAMEU_FR_PROC_SUB_TYPE_CROSS_SPECTRAL_DENSITY = 4,

    /** @brief Id for coherence data. */
    LAL_FRAMEU_FR_PROC_SUB_TYPE_COHERENCE = 5,

    /** @brief Id for transfer function data. */
    LAL_FRAMEU_FR_PROC_SUB_TYPE_TRANSFER_FUNCTION = 6
};

/* TODO: add routines:
int XLALFrameUFrFileIGWDVersion(LALFrameUFrFile *stream);
*/


/**
 * @name FrFile Routines
 * @{
 */

/**
 * @brief Close a FrFile stream.
 * @param stream Pointer to the FrFile stream to close.
 */
void XLALFrameUFrFileClose(LALFrameUFrFile * stream);

/**
 * @brief Open a frame file FrFile stream.
 * @param filename Filename of frame file to open.
 * @param mode Access mode: either "r" for read or "w" for write.
 * @return Pointer to a FrFile stream.
 * @retval NULL Failure.
 * @attention The calling routine is responsible for freeing the returned
 * pointer with XLALFrameUFrFileClose().
 */
LALFrameUFrFile *XLALFrameUFrFileOpen(const char *filename, const char *mode);

/**
 * @brief Use checksum to determine if FrFile stream is valid.
 * @param stream Pointer to the FrFile stream to check.
 * @warning Pointers to FrTOC might be left dangling by this function.
 * @return Logical value indicating if the FrFile stream checksum is correct.
 * @retval 1 The FrFile stream is valid.
 * @retval 0 The FrFile stream has an incorrect checksum or an error occurred.
 */
int XLALFrameUFileCksumValid(LALFrameUFrFile * stream);


/** @} */

/**
 * @name FrTOC Routines
 * @{
 */


/**
 * @brief Free a FrTOC structure.
 * @param toc Pointer to the FrTOC structure to free.
 * @sa Section 4.3.2.19 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
void XLALFrameUFrTOCFree(LALFrameUFrTOC * toc);

/** 
 * @brief Read the table of contents FrTOC structure for a FrFile stream.
 * @param stream Pointer to the input FrFile stream from which to read FrTOC.
 * @return Pointer to the FrTOC structure.
 * @retval NULL Failure.
 * @attention The calling routine is responsible for freeing the returned
 * pointer with XLALFrameUFrTOCFree().
 * @warning The pointer returned might be shallow and might be left
 * dangling if the stream is closed or reopened.  The FrTOC structure
 * should be freed before calling XLALFrameUFrFileClose() or
 * XLALFrameUFileCksumValid().
 * @sa Section 4.3.2.19 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
LALFrameUFrTOC *XLALFrameUFrTOCRead(LALFrameUFrFile * stream);


/** 
 * @name FrTOC Query Routines
 * @{
 */

/**
 * @brief Query FrTOC structure for number of FrameH structures contained.
 * @param toc Pointer to the FrTOC structure to be queried.
 * @return Number of FrameH structures contained.
 * @retval (size_t)(-1) Failure.
 * @sa Section 4.3.2.19 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
size_t XLALFrameUFrTOCQueryNFrame(const LALFrameUFrTOC * toc);

/**
 * @brief Query FrTOC structure for start time of a FrameH structure.
 * @param iptr Pointer to the integer number of seconds of the GPS start time.
 * @param toc Pointer to the FrTOC structure to be queried.
 * @param pos Index position of the FrameH structure.
 * @return Fractional part of the GPS start time.
 * @retval LAL_REAL8_FAIL_NAN Failure.
 * @sa Section 4.3.2.19 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
double XLALFrameUFrTOCQueryGTimeModf(double *iptr, const LALFrameUFrTOC * toc,
    size_t pos);

/** 
 * @brief Query FrTOC structure for duration of a FrameH structure.
 * @param toc Pointer to the FrTOC structure to be queried.
 * @param pos Index position of the FrameH structure.
 * @return Duration of the frame in seconds.
 * @retval LAL_REAL8_FAIL_NAN Failure.
 * @sa Section 4.3.2.19 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
double XLALFrameUFrTOCQueryDt(const LALFrameUFrTOC * toc, size_t pos);

/**
 * @brief Query FrTOC structure for number of FrAdcData structures.
 * @param toc Pointer to the FrTOC structure to be queried.
 * @return Number of FrAdcData structures.
 * @retval (size_t)(-1) Failure.
 * @sa Section 4.3.2.19 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
size_t XLALFrameUFrTOCQueryAdcN(const LALFrameUFrTOC * toc);

/**
 * @brief Query FrTOC structure for the name of a FrAdcData structure.
 * @param toc Pointer to the FrTOC structure to be queried.
 * @param adc Index position of the FrAdcData structure.
 * @return Pointer to string with the name of the FrAdcData structure.
 * @retval NULL Failure.
 * @warning The pointer returned is shallow and will be left dangling
 * when the FrTOC structure is freed.  The calling routine should not
 * attempt to free the returned pointer.
 * @sa Section 4.3.2.19 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
const char *XLALFrameUFrTOCQueryAdcName(const LALFrameUFrTOC * toc,
    size_t adc);

/**
 * @brief Query FrTOC structure for number of FrSimData structures.
 * @param toc Pointer to the FrTOC structure to be queried.
 * @return Number of FrSimData structures.
 * @retval (size_t)(-1) Failure.
 * @sa Section 4.3.2.19 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
size_t XLALFrameUFrTOCQuerySimN(const LALFrameUFrTOC * toc);

/**
 * @brief Query FrTOC structure for the name of a FrSimData structure.
 * @param toc Pointer to the FrTOC structure to be queried.
 * @param sim Index position of the FrSimData structure.
 * @return Pointer to string with the name of the FrSimData structure.
 * @retval NULL Failure.
 * @warning The pointer returned is shallow and will be left dangling
 * when the FrTOC structure is freed.  The calling routine should not
 * attempt to free the returned pointer.
 * @sa Section 4.3.2.19 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
const char *XLALFrameUFrTOCQuerySimName(const LALFrameUFrTOC * toc,
    size_t sim);

/**
 * @brief Query FrTOC structure for number of FrProcData structures.
 * @param toc Pointer to the FrTOC structure to be queried.
 * @return Number of FrProcData structures.
 * @retval (size_t)(-1) Failure.
 * @sa Section 4.3.2.19 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
size_t XLALFrameUFrTOCQueryProcN(const LALFrameUFrTOC * toc);

/**
 * @brief Query FrTOC structure for the name of a FrProcData structure.
 * @param toc Pointer to the FrTOC structure to be queried.
 * @param proc Index position of the FrProcData structure.
 * @return Pointer to string with the name of the FrProcData structure.
 * @retval NULL Failure.
 * @warning The pointer returned is shallow and will be left dangling
 * when the FrTOC structure is freed.  The calling routine should not
 * attempt to free the returned pointer.
 * @sa Section 4.3.2.19 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
const char *XLALFrameUFrTOCQueryProcName(const LALFrameUFrTOC * toc,
    size_t proc);

/**
 * @brief Query FrTOC structure for number of FrDetector structures.
 * @return Number of FrDetector structures.
 * @retval (size_t)(-1) Failure.
 * @sa Section 4.3.2.19 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
size_t XLALFrameUFrTOCQueryDetectorN(const LALFrameUFrTOC * toc);

/**
 * @brief Query FrTOC structure for the name of a FrDetector structure.
 * @param toc Pointer to the FrTOC structure to be queried.
 * @param det Index position of the FrDetector structure.
 * @return Pointer to string with the name of the FrDetector structure.
 * @retval NULL Failure.
 * @warning The pointer returned is shallow and will be left dangling
 * when the FrTOC structure is freed.  The calling routine should not
 * attempt to free the returned pointer.
 * @sa Section 4.3.2.19 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
const char *XLALFrameUFrTOCQueryDetectorName(const LALFrameUFrTOC * toc,
    size_t det);


/** @} */
/** @} */


/** 
 * @name FrameH Routines
 * @{
 */

/**
 * @brief Free a FrameH structure.
 * @param frame Pointer to the FrameH structure to free.
 */
void XLALFrameUFrameHFree(LALFrameUFrameH * frame);

/**
 * @brief Allocate memory for a new frame header FrameH structure.
 * @param name Name for this FrameH structure.
 * @param start GPS start time in seconds for this FrameH structure.
 * @param dt Duration in seconds for this FrameH structure.
 * @param frnum Number for this FrameH structure.
 * @return Pointer to a new FrameH structure.
 * @retval NULL Failure.
 * @attention The calling routine is responsible for freeing the returned
 * pointer with XLALFrameUFrameHFree().
 * @sa Section 4.3.2.3 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
LALFrameUFrameH *XLALFrameUFrameHAlloc(const char *name, double start,
    double dt, int frnum);

/**
 * @brief Read a frame header FrameH structure from a FrFile stream.
 * @param stream Pointer to the input FrFile stream.
 * @param pos Index position of the FrameH structure to read.
 * @return Pointer to a frame header FrameH structure.
 * @retval NULL Failure.
 * @attention The calling routine is responsible for freeing the returned
 * pointer with XLALFrameUFrameHFree().
 * @sa Section 4.3.2.3 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
LALFrameUFrameH *XLALFrameUFrameHRead(LALFrameUFrFile * stream, int pos);

/**
 * @brief Write a FrameH structure to an output FrFile stream.
 * @param stream Pointer to the output FrFile stream.
 * @param frame Pointer to the FrameH structure to be written.
 * @retval 0 Success.
 * @retval <0 Failure.
 * @sa Section 4.3.2.3 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameUFrameHWrite(LALFrameUFrFile * stream, LALFrameUFrameH * frame);


/**
 * @name FrameH Add Routines
 * @{
 */

/**
 * @brief Add a FrChan structure to a FrameH structure.
 * @param frame Pointer to the FrameH structure.
 * @param channel Pointer to the FrChan structure to be written.
 * @retval 0 Success.
 * @retval <0 Failure.
 * @sa Section 4.3.2.3 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameUFrameHFrChanAdd(LALFrameUFrameH * frame,
    LALFrameUFrChan * channel);

/**
 * @brief Add a FrDetector structure to a FrameH structure.
 * @param frame Pointer to the FrameH structure.
 * @param detector Pointer to the FrDetector structure to be written.
 * @retval 0 Success.
 * @retval <0 Failure.
 * @sa Section 4.3.2.3 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameUFrameHFrDetectorAdd(LALFrameUFrameH * frame,
    LALFrameUFrDetector * detector);

/**
 * @brief Add a FrHistory structure to a FrameH structure.
 * @param frame Pointer to the FrameH structure.
 * @param history Pointer to the FrHistory structure to be written.
 * @retval 0 Success.
 * @retval <0 Failure.
 * @sa Section 4.3.2.3 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameUFrameHFrHistoryAdd(LALFrameUFrameH * frame,
    LALFrameUFrHistory * history);


/** @} */


/**
 * @name FrameH Query Routines
 * @{
 */

/**
 * @brief Query FrameH structure for its name.
 * @param frame Pointer to the FrameH structure to be queried.
 * @return Pointer to string with the name of the FrameH structure.
 * @retval NULL Failure.
 * @warning The pointer returned is shallow and will be left dangling
 * when the FrameH structure is freed.  The calling routine should not
 * attempt to free the returned pointer.
 * @sa Section 4.3.2.3 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
const char *XLALFrameUFrameHQueryName(const LALFrameUFrameH * frame);

/**
 * @brief Query FrameH structure for the run number.
 * @param frame Pointer to the FrameH structure to be queried.
 * @return The run number associated with the FrameH structure.
 * @retval <0 Failure.
 * @sa Section 4.3.2.3 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameUFrameHQueryRun(const LALFrameUFrameH * frame);

/**
 * @brief Query FrameH structure for the frame number.
 * @param frame Pointer to the FrameH structure to be queried.
 * @return The frame number associated with the FrameH structure.
 * @retval <0 Failure.
 * @sa Section 4.3.2.3 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameUFrameHQueryFrame(const LALFrameUFrameH * frame);

/**
 * @brief Query FrameH structure for the data quality word.
 * @param frame Pointer to the FrameH structure to be queried.
 * @return The 32-bit data quality word expressed as an int type.
 * @retval <0 Failure.
 * @sa Section 4.3.2.3 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameUFrameHQueryDataQuality(const LALFrameUFrameH * frame);

/**
 * @brief Query FrameH structure for the start time.
 * @param iptr Pointer to the integer number of seconds of the GPS start time.
 * @param frame Pointer to the FrameH structure to be queried.
 * @return Fractional part of the GPS start time.
 * @retval LAL_REAL8_FAIL_NAN Failure.
 * @sa Section 4.3.2.3 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
double XLALFrameUFrameHQueryGTimeModf(double *iptr,
    const LALFrameUFrameH * frame);

/**
 * @brief Query FrameH structure for the number of leap seconds.
 * @param frame Pointer to the FrameH structure to be queried.
 * @return The number of leap seconds (TAI-UTC).
 * @retval <0 Failure.
 * @sa Section 4.3.2.3 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameUFrameHQueryULeapS(const LALFrameUFrameH * frame);

/**
 * @brief Query FrameH structure for the duration.
 * @param frame Pointer to the FrameH structure to be queried.
 * @return The duration of the FrameH structure in seconds.
 * @retval LAL_REAL8_FAIL_NAN Failure.
 * @sa Section 4.3.2.3 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
double XLALFrameUFrameHQueryDt(const LALFrameUFrameH * frame);

/** @} */


/**
 * @name FrameH Set Routines
 * @{
 */

/**
 * @brief Set the run number in a FrameH structure.
 * @param frame Pointer to the FrameH structure to be modified.
 * @param run The value of the run number to be set.
 * @retval 0 Success.
 * @retval <0 Failure.
 * @sa Section 4.3.2.3 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameUFrameHSetRun(LALFrameUFrameH * frame, int run);


/** @} */
/** @} */


/** 
 * @name FrChan Routines
 * @{
 */

/**
 * @brief Free a FrChan structure.
 * @param channel Pointer to the FrChan structure to free.
 * @sa Sections 4.3.2.4, 4.3.2.11, and 4.3.2.14 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
void XLALFrameUFrChanFree(LALFrameUFrChan * channel);

/**
 * @brief Read a channel FrChan structure from a FrFile stream.
 * @param stream Pointer to the input FrFile stream.
 * @param name Name of the FrChan structure to be read.
 * @param pos Index position of the FrameH structure containing the FrChan.
 * @return Pointer to a channel FrChan structure.
 * @retval NULL Failure.
 * @attention The calling routine is responsible for freeing the returned
 * pointer with XLALFrameUFrChanFree().
 * @sa Sections 4.3.2.4, 4.3.2.11, and 4.3.2.14 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
LALFrameUFrChan *XLALFrameUFrChanRead(LALFrameUFrFile * stream,
    const char *name, size_t pos);

/**
 * @brief Allocate memory for a new FrChan structure of FrAdcData type.
 * @param name Name for this FrChan structure.
 * @param dtype Data type code given in #LALFrameUFrVectDataType.
 * @param ndata Number of data points contained in this FrChan structure.
 * @return Pointer to a new FrChan structure.
 * @retval NULL Failure.
 * @attention The calling routine is responsible for freeing the returned
 * pointer with XLALFrameUFrChanFree().
 * @sa Sections 4.3.2.4 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
LALFrameUFrChan *XLALFrameUFrAdcChanAlloc(const char *name, int dtype,
    size_t ndata);

/**
 * @brief Allocate memory for a new FrChan structure of FrSimData type.
 * @param name Name for this FrChan structure.
 * @param dtype Data type code given in #LALFrameUFrVectDataType.
 * @param ndata Number of data points contained in this FrChan structure.
 * @return Pointer to a new FrChan structure.
 * @retval NULL Failure.
 * @attention The calling routine is responsible for freeing the returned
 * pointer with XLALFrameUFrChanFree().
 * @sa Sections 4.3.2.11 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
LALFrameUFrChan *XLALFrameUFrSimChanAlloc(const char *name, int dtype,
    size_t ndata);

/**
 * @brief Allocate memory for a new FrChan structure of FrAdcData type.
 * @param name Name for this FrChan structure.
 * @param type FrProcData type code given in #LALFrameUFrProcDataType.
 * @param subtype FrProcData subtype code given in #LALFrameUFrProcDataSubType.
 * @param dtype Data type code given in #LALFrameUFrVectDataType.
 * @param ndata Number of data points contained in this FrChan structure.
 * @return Pointer to a new FrChan structure.
 * @retval NULL Failure.
 * @attention The calling routine is responsible for freeing the returned
 * pointer with XLALFrameUFrChanFree().
 * @sa Sections 4.3.2.14 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
LALFrameUFrChan *XLALFrameUFrProcChanAlloc(const char *name, int type,
    int subtype, int dtype, size_t ndata);


/**
 * @name FrChan Query Routines
 * @{
 */

/**
 * @brief Query FrChan structure for its name.
 * @param channel Pointer to the FrChan structure to be queried.
 * @return Pointer to string with the name of the FrChan structure.
 * @retval NULL Failure.
 * @warning The pointer returned is shallow and will be left dangling
 * when the FrChan structure is freed.  The calling routine should not
 * attempt to free the returned pointer.
 * @sa Sections 4.3.2.4, 4.3.2.11, and 4.3.2.14 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
const char *XLALFrameUFrChanQueryName(const LALFrameUFrChan * channel);

/**
 * @brief Query FrChan structure for time offset for this channel.
 * @param channel Pointer to the FrChan structure to be queried.
 * @return The offset of the first sample relative to the frame start
 * time in seconds.
 * @retval LAL_REAL8_FAIL_NAN Failure.
 * @remark The offset must be positive and smaller than the frame length.
 * Time offsets are always added together: to obtain the time of the
 * first sample, add the time offset to the frame start time as obtained
 * from the routine XLALFrameUFrTOCQueryGTimeModf() or from the routine
 * XLALFrameUFrameHQueryGTimeModf().
 * @sa Sections 4.3.2.4, 4.3.2.11, and 4.3.2.14 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
double XLALFrameUFrChanQueryTimeOffset(const LALFrameUFrChan * channel);

/** @} */


/**
 * @name FrChan Set Routines
 * @{
 */

/**
 * @brief Set the sample rate in a FrChan structure.
 * @param channel Pointer to the FrChan structure to be modified.
 * @param sampleRate The value of the sample rate in Hertz to be set.
 * @retval 0 Success.
 * @retval <0 Failure.
 * @sa Sections 4.3.2.4, 4.3.2.11, and 4.3.2.14 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameUFrChanSetSampleRate(LALFrameUFrChan * channel,
    double sampleRate);

/** @} */
/** @} */


/**
 * @name FrVect Routines
 * @{
 */

/**
 * @brief Allocate memory for a FrVect structure within a FrChan structure.
 * @param channel Pointer to the FrChan structure to be modified.
 * @param dtype Data type code given in #LALFrameUFrVectDataType.
 * @param ndata Number of data points contained the FrVect structure.
 * @retval 0 Success.
 * @retval <0 Failure.
 * @sa Sections 4.3.2.20 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameUFrChanVectorAlloc(LALFrameUFrChan * channel, int dtype,
    size_t ndata);

/**
 * @brief Compress a FrVect structure within a FrChan structure.
 * @param channel Pointer to the FrChan structure to be modified.
 * @param compressLevel Compression scheme given in
 * #LALFrameUFrVectCompressionScheme.
 * @retval 0 Success.
 * @retval <0 Failure.
 * @sa Sections 4.3.2.20 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameUFrChanVectorCompress(LALFrameUFrChan * channel,
    int compressLevel);

/**
 * @brief Expands a FrVect structure within a FrChan structure.
 * @param channel Pointer to the FrChan structure to be modified.
 * @retval 0 Success.
 * @retval <0 Failure.
 * @sa Sections 4.3.2.20 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameUFrChanVectorExpand(LALFrameUFrChan * channel);


/**
 * @name FrVect Query Routines
 * @{
 */

/**
 * @brief Query FrChan structure for the name of its FrVect structure.
 * @param channel Pointer to the FrChan structure to be queried.
 * @return Pointer to string with the name of the FrVect structure.
 * @retval NULL Failure.
 * @warning The pointer returned is shallow and will be left dangling
 * when the FrChan structure is freed.  The calling routine should not
 * attempt to free the returned pointer.
 * @sa Sections 4.3.2.20 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
const char *XLALFrameUFrChanVectorQueryName(const LALFrameUFrChan * channel);

/**
 * @brief Query FrChan structure for the compression scheme of its
 * FrVect structure.
 * @param channel Pointer to the FrChan structure to be queried.
 * @return The compression scheme of the FrVect structure as given in
 * #LALFrameUFrVectCompressionScheme.
 * @retval <0 Failure.
 * @sa Sections 4.3.2.20 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameUFrChanVectorQueryCompress(const LALFrameUFrChan * channel);

/**
 * @brief Query FrChan structure for the data type of its FrVect structure.
 * @param channel Pointer to the FrChan structure to be queried.
 * @return The data type of the FrVect structure as given in
 * #LALFrameUFrVectDataType.
 * @retval <0 Failure.
 * @sa Sections 4.3.2.20 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameUFrChanVectorQueryType(const LALFrameUFrChan * channel);

/**
 * @brief Query FrChan structure for the data pointer in its FrVect structure.
 * @param channel Pointer to the FrChan structure to be queried.
 * @return Pointer to the data in the FrVect structure.
 * @retval NULL Failure.
 * @warning The pointer returned is shallow and will be left dangling
 * when the FrChan structure is freed.  The calling routine should not
 * attempt to free the returned pointer.
 * @sa Sections 4.3.2.20 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
void *XLALFrameUFrChanVectorQueryData(const LALFrameUFrChan * channel);

/**
 * @brief Query FrChan structure for the number of bytes of data in its
 * FrVect structure.
 * @param channel Pointer to the FrChan structure to be queried.
 * @return The number of bytes of data in the FrVect structure.
 * @retval (size_t)(-1) Failure.
 * @sa Sections 4.3.2.20 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
size_t XLALFrameUFrChanVectorQueryNBytes(const LALFrameUFrChan * channel);

/**
 * @brief Query FrChan structure for the number of points of data in its
 * FrVect structure.
 * @param channel Pointer to the FrChan structure to be queried.
 * @return The number of bytes of data in the FrVect structure.
 * @retval (size_t)(-1) Failure.
 * @sa Sections 4.3.2.20 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
size_t XLALFrameUFrChanVectorQueryNData(const LALFrameUFrChan * channel);

/**
 * @brief Query FrChan structure for the number of dimensions of
 * the multi-dimensional data in its FrVect structure.
 * @param channel Pointer to the FrChan structure to be queried.
 * @return The number of dimensions of the data in the FrVect structure.
 * @retval (size_t)(-1) Failure.
 * @sa Sections 4.3.2.20 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
size_t XLALFrameUFrChanVectorQueryNDim(const LALFrameUFrChan * channel);

/**
 * @brief Query FrChan structure for the number of points in the @a dim 
 * dimension of the multi-dimensional data in the FrVect structure.
 * @param channel Pointer to the FrChan structure to be queried.
 * @param dim The dimension of the multi-dimensional data.
 * @return The number of points in the @a dim dimension of the data in the
 * FrVect structure.
 * @retval (size_t)(-1) Failure.
 * @sa Sections 4.3.2.20 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
size_t XLALFrameUFrChanVectorQueryNx(const LALFrameUFrChan * channel,
    size_t dim);

/**
 * @brief Query FrChan structure for the sampling interval in the @a dim 
 * dimension of the multi-dimensional data in the FrVect structure.
 * @param channel Pointer to the FrChan structure to be queried.
 * @param dim The dimension of the multi-dimensional data.
 * @return The sampling interval in the @a dim dimension of the data in the
 * FrVect structure.
 * @retval LAL_REAL8_FAIL_NAN Failure.
 * @sa Sections 4.3.2.20 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
double XLALFrameUFrChanVectorQueryDx(const LALFrameUFrChan * channel,
    size_t dim);

/**
 * @brief Query FrChan structure for the starting value of the @a dim 
 * dimension of the multi-dimensional data in the FrVect structure.
 * @param channel Pointer to the FrChan structure to be queried.
 * @param dim The dimension of the multi-dimensional data.
 * @return The starting value of the @a dim dimension of the data in the
 * FrVect structure.
 * @retval LAL_REAL8_FAIL_NAN Failure.
 * @sa Sections 4.3.2.20 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
double XLALFrameUFrChanVectorQueryStartX(const LALFrameUFrChan * channel,
    size_t dim);

/**
 * @brief Query FrChan structure for the units of the domain of the @a dim 
 * dimension of the multi-dimensional data in the FrVect structure.
 * @param channel Pointer to the FrChan structure to be queried.
 * @param dim The dimension of the multi-dimensional data.
 * @return Pointer to a string containing the units of the domain of the @a dim
 * dimension of the data in the FrVect structure.
 * @retval NULL Failure.
 * @warning The pointer returned is shallow and will be left dangling
 * when the FrChan structure is freed.  The calling routine should not
 * attempt to free the returned pointer.
 * @sa Sections 4.3.2.20 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
const char *XLALFrameUFrChanVectorQueryUnitX(const LALFrameUFrChan * channel,
    size_t dim);

/**
 * @brief Query FrChan structure for sample units of the data in the FrVect
 * structure.
 * @param channel Pointer to the FrChan structure to be queried.
 * @return Pointer to a string containing the sample units of the data in the
 * FrVect structure.
 * @retval NULL Failure.
 * @warning The pointer returned is shallow and will be left dangling
 * when the FrChan structure is freed.  The calling routine should not
 * attempt to free the returned pointer.
 * @sa Sections 4.3.2.20 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
const char *XLALFrameUFrChanVectorQueryUnitY(const LALFrameUFrChan * channel);

/** @} */


/**
 * @name FrVect Set Routines
 * @{
 */

/**
 * @brief Set the name of the FrVect structure contained in a FrChan structure.
 * @param channel Pointer to the FrChan structure to be modified.
 * @param name Pointer to a string with the name for the FrVect structure.
 * @retval 0 Success.
 * @retval <0 Failure.
 * @sa Sections 4.3.2.20 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameUFrChanVectorSetName(LALFrameUFrChan * channel,
    const char *name);

/**
 * @brief Set the sampling interval for the data in the FrVect structure
 * contained in a FrChan structure.
 * @param channel Pointer to the FrChan structure to be modified.
 * @param dx The sampling interval for the data.
 * @retval 0 Success.
 * @retval <0 Failure.
 * @sa Sections 4.3.2.20 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameUFrChanVectorSetDx(LALFrameUFrChan * channel, double dx);

/**
 * @brief Set the starting value for the domain of the data in the FrVect
 * structure contained in a FrChan structure.
 * @param channel Pointer to the FrChan structure to be modified.
 * @param x0 The starting value for the domain of the data in the FrVect
 * structure.
 * @retval 0 Success.
 * @retval <0 Failure.
 * @sa Sections 4.3.2.20 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameUFrChanVectorSetStartX(LALFrameUFrChan * channel, double x0);

/**
 * @brief Set the units of the domain of the data in the FrVect structure
 * contained in a FrChan structure.
 * @param channel Pointer to the FrChan structure to be modified.
 * @param unit Pointer to a string containing the units of the domain of the
 * data in the FrVect structure.
 * @retval 0 Success.
 * @retval <0 Failure.
 * @sa Sections 4.3.2.20 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameUFrChanVectorSetUnitX(LALFrameUFrChan * channel,
    const char *unit);

/**
 * @brief Set the units of the data in the FrVect structure contained in a
 * FrChan structure.
 * @param channel Pointer to the FrChan structure to be modified.
 * @param unit Pointer to a string containing the units of the data in the
 * FrVect structure.
 * @retval 0 Success.
 * @retval <0 Failure.
 * @sa Sections 4.3.2.20 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
int XLALFrameUFrChanVectorSetUnitY(LALFrameUFrChan * channel,
    const char *unit);

/* TODO: a bunch more things to set coming up!!! */

/** @} */
/** @} */


/**
 * @name FrDetector Routines
 * @{
 */

/**
 * @brief Free a FrDetector structure.
 * @param detector Pointer to the FrDetector structure to free.
 * @sa Section 4.3.2.5 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
void XLALFrameUFrDetectorFree(LALFrameUFrDetector * detector);

/**
 * @brief Read a detector FrDetector structure from a FrFile stream.
 * @param stream Pointer to the input FrFile stream.
 * @param name Name of the FrDetector structure to be read.
 * @return Pointer to a detector FrDetector structure.
 * @retval NULL Failure.
 * @attention The calling routine is responsible for freeing the returned
 * pointer with XLALFrameUFrDetectorFree().
 * @sa Section 4.3.2.5 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
LALFrameUFrDetector *XLALFrameUFrDetectorRead(LALFrameUFrFile * stream,
    const char *name);

/**
 * @brief Allocate memory for a new detector FrDetector structure.
 * @param name Name for this FrDetector structure.
 * @param latitude Latitude North of the equator of detector in radians.
 * @param longitude Longitude East of Greenwich meridian of detector in radians.
 * @param elevation Elevation of detector above WGS84 ellipsoid in meters.
 * @param azimuthX Orientation of x-arm in radians East of North.
 * @param azimuthY Orientation of y-arm in radians East of North.
 * @param altitudeX Altitude (pitch) angle of x-arm in radians above horizaon
 * (local tangent to WGS84 ellipsoid).
 * @param altitudeY Altitude (pitch) angle of y-arm in radians above horizaon
 * (local tangent to WGS84 ellipsoid).
 * @param midpointX Distance in meters between detector vertex and the middle
 * of the x-arm.
 * @param midpointY Distance in meters between detector vertex and the middle
 * of the y-arm.
 * @param localTime Local seasonal time minus UTC in seconds.
 * @return Pointer to a new FrDetector structure.
 * @retval NULL Failure.
 * @attention The calling routine is responsible for freeing the returned
 * pointer with XLALFrameUFrDetectorFree().
 * @sa Section 4.3.2.5 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 * @sa <em>Beam Pattern Response Functions and Times of Arrival for Earthbound
 * Interferometer</em>
 * LIGO-T010110 [https://dcc.ligo.org/LIGO-T010110-v1/public]
 */
LALFrameUFrDetector *XLALFrameUFrDetectorAlloc(const char *name,
    const char *prefix, double latitude, double longitude, double elevation,
    double azimuthX, double azimuthY, double altitudeX, double altitudeY,
    double midpointX, double midpointY, int localTime);

/**
 * @name FrDetector Query Routines
 * @{ 
 */

/**
 * @brief Query FrDetector structure for the detector name.
 * @param detector Pointer to the FrDetector structure to be queried.
 * @return Pointer to string with the name of the FrDetector structure.
 * @retval NULL Failure.
 * @warning The pointer returned is shallow and will be left dangling
 * when the FrDetector structure is freed.  The calling routine should not
 * attempt to free the returned pointer.
 * @sa Section 4.3.2.5 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
const char *XLALFrameUFrDetectorQueryName(const LALFrameUFrDetector *
    detector);

/**
 * @brief Query FrDetector structure for the detector prefix.
 * @param detector Pointer to the FrDetector structure to be queried.
 * @return Pointer to string with the prefix of the detector.
 * @retval NULL Failure.
 * @warning The pointer returned is shallow and will be left dangling
 * when the FrDetector structure is freed.  The calling routine should not
 * attempt to free the returned pointer.
 * @sa Section 4.3.2.5 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
const char *XLALFrameUFrDetectorQueryPrefix(const LALFrameUFrDetector *
    detector);

/**
 * @brief Query FrDetector structure for the detector longitude.
 * @param detector Pointer to the FrDetector structure to be queried.
 * @return The detector longitude in radians East of the Greenwich meridian.
 * @retval LAL_REAL8_FAIL_NAN Failure.
 * @sa Section 4.3.2.5 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 * @sa <em>Beam Pattern Response Functions and Times of Arrival for Earthbound
 * Interferometer</em>
 * LIGO-T010110 [https://dcc.ligo.org/LIGO-T010110-v1/public]
 */
double XLALFrameUFrDetectorQueryLongitude(const LALFrameUFrDetector *
    detector);

/**
 * @brief Query FrDetector structure for the detector latitude.
 * @param detector Pointer to the FrDetector structure to be queried.
 * @return The detector latitude in radians North of the equator.
 * @retval LAL_REAL8_FAIL_NAN Failure.
 * @sa Section 4.3.2.5 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 * @sa <em>Beam Pattern Response Functions and Times of Arrival for Earthbound
 * Interferometer</em>
 * LIGO-T010110 [https://dcc.ligo.org/LIGO-T010110-v1/public]
 */
double XLALFrameUFrDetectorQueryLatitude(const LALFrameUFrDetector *
    detector);

/**
 * @brief Query FrDetector structure for the detector elevation.
 * @param detector Pointer to the FrDetector structure to be queried.
 * @return The detector elevation in meters relative to the WGS84 ellipsoid.
 * @retval LAL_REAL8_FAIL_NAN Failure.
 * @sa Section 4.3.2.5 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 * @sa <em>Beam Pattern Response Functions and Times of Arrival for Earthbound
 * Interferometer</em>
 * LIGO-T010110 [https://dcc.ligo.org/LIGO-T010110-v1/public]
 */
double XLALFrameUFrDetectorQueryElevation(const LALFrameUFrDetector *
    detector);

/**
 * @brief Query FrDetector structure for the detector x-arm azimuth.
 * @param detector Pointer to the FrDetector structure to be queried.
 * @return The orientation of the x-arm in radians East of North.
 * @retval LAL_REAL8_FAIL_NAN Failure.
 * @sa Section 4.3.2.5 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 * @sa <em>Beam Pattern Response Functions and Times of Arrival for Earthbound
 * Interferometer</em>
 * LIGO-T010110 [https://dcc.ligo.org/LIGO-T010110-v1/public]
 */
double XLALFrameUFrDetectorQueryArmXAzimuth(const LALFrameUFrDetector *
    detector);

/**
 * @brief Query FrDetector structure for the detector y-arm azimuth.
 * @param detector Pointer to the FrDetector structure to be queried.
 * @return The orientation of the y-arm in radians East of North.
 * @retval LAL_REAL8_FAIL_NAN Failure.
 * @sa Section 4.3.2.5 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 * @sa <em>Beam Pattern Response Functions and Times of Arrival for Earthbound
 * Interferometer</em>
 * LIGO-T010110 [https://dcc.ligo.org/LIGO-T010110-v1/public]
 */
double XLALFrameUFrDetectorQueryArmYAzimuth(const LALFrameUFrDetector *
    detector);

/**
 * @brief Query FrDetector structure for the detector x-arm altitude.
 * @param detector Pointer to the FrDetector structure to be queried.
 * @return The altitude (pitch) of the x-arm in radians above horizontal
 * (local tangent to WGS84 ellipsoid).
 * @retval LAL_REAL8_FAIL_NAN Failure.
 * @sa Section 4.3.2.5 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 * @sa <em>Beam Pattern Response Functions and Times of Arrival for Earthbound
 * Interferometer</em>
 * LIGO-T010110 [https://dcc.ligo.org/LIGO-T010110-v1/public]
 */
double XLALFrameUFrDetectorQueryArmXAltitude(const LALFrameUFrDetector *
    detector);

/**
 * @brief Query FrDetector structure for the detector y-arm altitude.
 * @param detector Pointer to the FrDetector structure to be queried.
 * @return The altitude (pitch) of the y-arm in radians above horizontal
 * (local tangent to WGS84 ellipsoid).
 * @retval LAL_REAL8_FAIL_NAN Failure.
 * @sa Section 4.3.2.5 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 * @sa <em>Beam Pattern Response Functions and Times of Arrival for Earthbound
 * Interferometer</em>
 * LIGO-T010110 [https://dcc.ligo.org/LIGO-T010110-v1/public]
 */
double XLALFrameUFrDetectorQueryArmYAltitude(const LALFrameUFrDetector *
    detector);

/**
 * @brief Query FrDetector structure for the detector x-arm midpoint.
 * @param detector Pointer to the FrDetector structure to be queried.
 * @return The distance in meters between the detector vertex and the
 * midpoint of the x-arm.
 * @retval LAL_REAL8_FAIL_NAN Failure.
 * @sa Section 4.3.2.5 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 * @sa <em>Beam Pattern Response Functions and Times of Arrival for Earthbound
 * Interferometer</em>
 * LIGO-T010110 [https://dcc.ligo.org/LIGO-T010110-v1/public]
 */
double XLALFrameUFrDetectorQueryArmXMidpoint(const LALFrameUFrDetector *
    detector);

/**
 * @brief Query FrDetector structure for the detector y-arm midpoint.
 * @param detector Pointer to the FrDetector structure to be queried.
 * @return The distance in meters between the detector vertex and the
 * midpoint of the y-arm.
 * @retval LAL_REAL8_FAIL_NAN Failure.
 * @sa Section 4.3.2.5 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 * @sa <em>Beam Pattern Response Functions and Times of Arrival for Earthbound
 * Interferometer</em>
 * LIGO-T010110 [https://dcc.ligo.org/LIGO-T010110-v1/public]
 */
double XLALFrameUFrDetectorQueryArmYMidpoint(const LALFrameUFrDetector *
    detector);

/**
 * @brief Query FrDetector structure for the local time offset at the detector.
 * @param detector Pointer to the FrDetector structure to be queried.
 * @return The local seasonal time minus UTC in seconds.
 * @sa Section 4.3.2.5 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 * @sa <em>Beam Pattern Response Functions and Times of Arrival for Earthbound
 * Interferometer</em>
 * LIGO-T010110 [https://dcc.ligo.org/LIGO-T010110-v1/public]
 */
int XLALFrameUFrDetectorQueryLocalTime(const LALFrameUFrDetector * detector);

/** @} */
/** @} */

/**
 * @name FrHistory Routines
 * @{
 */

/**
 * @brief Free a FrHistory structure.
 * @param history Pointer to the FrHistory structure to free.
 * @sa Section 4.3.2.9 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
void XLALFrameUFrHistoryFree(LALFrameUFrHistory * history);

/**
 * @brief Allocate memory for a new detector FrHistory structure.
 * @param name Name for this FrDetector structure.
 * @param gpssec Time stamp in GPS seconds for this FrHistory structure.
 * @param comment History comment for this FrHistory structure.
 * @return Pointer to a new FrHistory structure.
 * @retval NULL Failure.
 * @attention The calling routine is responsible for freeing the returned
 * pointer with XLALFrameUFrHistoryFree().
 * @sa Section 4.3.2.9 of
 * <em>Specification of a Common Data Frame Format for Interferometric
 * Gravitational Wave Detectors (IGWD)</em>
 * LIGO-T970130 [https://dcc.ligo.org/LIGO-T970130-v1/public].
 */
LALFrameUFrHistory *XLALFrameUFrHistoryAlloc(const char *name, double gpssec,
    const char *comment);

/** @} */

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif /* _LALFRAMEU_H */
