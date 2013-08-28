/*
 * Copyright (C) 2009 Reinhard Prix
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
 * \defgroup CFSv3 F-statistic (v3)
 * \ingroup pulsarTODO
 * \author Reinhard Prix
 * \brief Module for FFT-based F-statistic calculation
 *
 */

/**
 * \file
 * \ingroup CFSv3
 * \author Reinhard Prix
 * \brief Header file defining the API for the CFSv3 module
 */

#ifndef _FSTAT_V3_H  /* Double-include protection. */
#define _FSTAT_V3_H

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif


/*---------- exported INCLUDES ----------*/
#include <lal/SFTutils.h>
#include <lal/LALDatatypes.h>

/*---------- exported DEFINES ----------*/

/*---------- exported types ----------*/

/*---------- exported Global variables ----------*/
/* empty init-structs for the types defined in here */

/*---------- exported prototypes [API] ----------*/

COMPLEX8TimeSeries *XLALSFTVectorToCOMPLEX8TimeSeries ( SFTVector *sfts,      /**< input SFT vector, gets modified! */
							const LIGOTimeGPS *start_in,    /**< input start time */
							const LIGOTimeGPS *end_in       /**< input end time */
							);

SFTtype *XLALSFTVectorToLFT ( const SFTVector *sfts, REAL8 upsampling );


int XLALReorderFFTWtoSFT (COMPLEX8Vector *X);
int XLALReorderSFTtoFFTW (COMPLEX8Vector *X);
int XLALTimeShiftSFT ( SFTtype *sft, REAL8 shift );
int XLALMultiplySFTbyCOMPLEX8 ( SFTtype *sft, COMPLEX8 factor );


#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif  /* Double-include protection. */



