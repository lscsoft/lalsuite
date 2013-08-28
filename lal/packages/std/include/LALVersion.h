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

#ifndef _LALVERSION_H
#define _LALVERSION_H

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#elif 0
}       /* so that editors will match preceding brace */
#endif

/**
 * \addtogroup LALVersion_h
 *
 * \brief Provides routines for reporting the LAL version.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/LALVersion.h>
 * \endcode
 *
 * \section ss_globals Global variables
 *
 * \code
 * extern const char *lalVersion;
 * extern const int   lalVersionMajor;
 * extern const int   lalVersionMinor;
 * extern const char *lalConfigureArgs;
 * extern const char *lalConfigureDate;
 * \endcode
 *
 * These constant variables are set at compile time and included into the LAL
 * library.  They contain the information about the version of LAL and the
 * configuration information.
 *
 */
/*@{*/

/**\name Error Codes *//*@{ */
#define LALVERSIONH_ENULL 1     /**< Null string pointer */
#define LALVERSIONH_ESIZE 2     /**< Zero string size */
#define LALVERSIONH_ESPRN 4     /**< Error in snprintf */
#define LALVERSIONH_ESHRT 8     /**< String too short */
/*@}*/

#define LALVERSIONH_MSGENULL "Null string pointer."
#define LALVERSIONH_MSGESIZE "Zero string size."
#define LALVERSIONH_MSGESPRN "Error in snprintf."
#define LALVERSIONH_MSGESHRT "String too short."


extern const char *const lalVersion;
extern const int lalVersionMajor;
extern const int lalVersionMinor;
extern const int lalVersionMicro;
extern const int lalVersionDevel;
extern const char *const lalBuildDate;
extern const char *const lalConfigureArgs;
extern const char *const lalConfigureDate;


void LALVersion(LALStatus * status, CHAR * message, UINT4 size,
                INT4 verbose);

/*@}*/

#if 0
{       /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif
#endif /* _LALVERSION_H */
