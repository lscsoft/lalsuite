/*
*  Copyright (C) 2007 Jolien Creighton, Patrick Brady
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

/*-----------------------------------------------------------------------
 *
 * File Name: ReadNoiseSpectrum.h
 *
 * Author: Brady, P. R.
 *
 *
 *-----------------------------------------------------------------------
 */

#ifndef _READNOISESPECTRUMH_H
#define _READNOISESPECTRUMH_H

#include <lal/LALDatatypes.h>
#include <lal/Date.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * \addtogroup ReadNoiseSpectrum_h
 * \author Brady, P. R.
 *
 * \brief Provides function to read in a file containing a possibly unequally sampled
 * noise amplitude spectrum (\f$\textrm{strain}/\sqrt(\textrm{Hz})\f$) and return as
 * a frequency series.
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/ReadNoiseSpectrum.h>
 * \endcode
 */
/*@{*/

/**\name Error Codes */ /*@{*/
#define LALREADNOISESPECTRUMH_ENULL 1   /**< Null pointer */
#define LALREADNOISESPECTRUMH_ENNUL 2   /**< Non-null pointer */
#define LALREADNOISESPECTRUMH_EALOC 3   /**< Memory allocation error */
#define LALREADNOISESPECTRUMH_EOPEN 4   /**< Error opening file */
#define LALREADNOISESPECTRUMH_EFCLO 5   /**< Error closing file */
#define LALREADNOISESPECTRUMH_EPARS 8   /**< Error parsing spectrum file */
/*@}*/
/*@}*/

#define LALREADNOISESPECTRUMH_MSGENULL "Null pointer"
#define LALREADNOISESPECTRUMH_MSGENNUL "Non-null pointer"
#define LALREADNOISESPECTRUMH_MSGEALOC "Memory allocation error"
#define LALREADNOISESPECTRUMH_MSGEOPEN "Error opening file"
#define LALREADNOISESPECTRUMH_MSGEFCLO "Error closing file"
#define LALREADNOISESPECTRUMH_MSGEPARS "Error parsing spectrum file"


#define LALREADNOISESPECTRUM_MAXLINELENGTH 2048

void LALReadNoiseSpectrum(
    LALStatus *status,
    REAL4FrequencySeries *spectrum,
    CHAR *fname
    );


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _READNOISESPECTRUMH_H */
