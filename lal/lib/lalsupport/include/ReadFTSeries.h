/*
*  Copyright (C) 2007 Jolien Creighton
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

#ifndef _READFTSERIES_H
#define _READFTSERIES_H

#include <lal/LALStdlib.h>

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \addtogroup ReadFTSeries_h
 * \author Torres, C. W.
 *
 * \brief This is a simple utility to Read time and frequency series into a file.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/ReadFTSeries.h>
 * \endcode
 *
 * Provides prototype information for the routines in \ref ReadTimeSeries_c and \ref ReadFrequencySeries_c.
 */
/*@{*/

/**\name Error Codes */ /*@{*/
#define  READFTSERIESH_EFILENOTFOUND       1    /**< Invalid Filename or File Not Found */
#define  READFTSERIESH_EPARSE              2    /**< Error Parsing File */
/*@}*/
/** \cond DONT_DOXYGEN */
#define  READFTSERIESH_MSGEFILENOTFOUND    "Invalid Filename or File Not Found"
#define  READFTSERIESH_MSGEPARSE           "Error Parsing File"
/** \endcond */


#ifndef SWIG /* exclude from SWIG interface */
/** \cond DONT_DOXYGEN */
enum enumLALSupportUnitTextSize {
  LALSupportUnitTextSize = sizeof("10^-32768 m^-32768/32767 "
                                  "kg^-32768/32767 "
                                  "s^-32768/32767 A^-32768/32767 "
                                  "K^-32768/32767 strain^-32768/32767 "
                                  "count^-32768/32767")
};
enum enumMaxLineLength {
  MaxLineLength = LALSupportUnitTextSize + sizeof("Units are ()\n")
};
/** \endcond */
#endif /* SWIG */


/**
 * \defgroup ReadTimeSeries_c Module ReadTimeSeries.c
 * \author Torres, C. V.
 *
 * \brief Each member of this family of functions reads from a file the output of the corresponding \c PrintTimeSeries routine.
 *
 * ### Notes ###
 *
 * These functions perform I/O operations, which are not a part of LAL
 * proper They should only be used for debugging purposes in test
 * functions, not in any production code.
 */
/*@{*/
void LALReadTimeSeries(LALStatus* status,  REAL4TimeSeries *series , const CHAR *filename );
void LALSReadTimeSeries(LALStatus* status,  REAL4TimeSeries *series , const CHAR *filename );
void LALDReadTimeSeries(LALStatus* status,  REAL8TimeSeries *series , const CHAR *filename );
void LALCReadTimeSeries(LALStatus* status,  COMPLEX8TimeSeries *series , const CHAR *filename );
void LALZReadTimeSeries(LALStatus* status,  COMPLEX16TimeSeries *series , const CHAR *filename );
/*@}*/

/**
 * \defgroup ReadFrequencySeries_c Module ReadFrequencySeries.c
 * \author Torres, C. V.
 *
 * \brief Each member of this family of functions reads from a file the output of the corresponding \c PrintFrequencySeries routine.
 *
 * ### Notes ###
 *
 * These functions perform I/O operations, which are not a part of LAL
 * proper. They should only be used for debugging purposes in test
 * functions, not in any production code.
 */
/*@{*/
void LALReadFrequencySeries(LALStatus* status,  REAL4FrequencySeries *series , const CHAR *filename );
void LALSReadFrequencySeries(LALStatus* status,  REAL4FrequencySeries *series , const CHAR *filename );
void LALDReadFrequencySeries(LALStatus* status,  REAL8FrequencySeries *series , const CHAR *filename );
void LALCReadFrequencySeries(LALStatus* status, COMPLEX8FrequencySeries *series , const CHAR *filename );
void LALZReadFrequencySeries(LALStatus* status,  COMPLEX16FrequencySeries *series , const CHAR *filename );
/*@}*/

/*@}*/

#ifdef  __cplusplus
}
#endif

#endif /* _READSERIES_H */
