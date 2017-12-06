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

#ifndef _LALCORRELATION_H
#define _LALCORRELATION_H

#include <lal/LALStdlib.h>
/* ****** INCLUDE ANY OTHER LAL HEADERS needed for header (NOT module) ****/

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \addtogroup LALCorrelation_h
 * \author Yakushin, Igor
 *
 * \brief [One sentence briefly defining scope of the header]
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/LALCorrelation.h>
 * \endcode
 *
 * [Generic documentation on the header; this is the main place to
 * document any stuff not specific to the module]
 *
 */
/*@{*/

/**\name Error Codes */
/*@{*/
#define LALCORRELATIONH_ENULLP        1		/**< Null pointer */
#define LALCORRELATIONH_ESTART        2		/**< Time series do not start simultaneously */
#define LALCORRELATIONH_ESAMPLING     3		/**< Time series are not sampled with the same rate */
/*@}*/

/** \cond DONT_DOXYGEN */
#define LALCORRELATIONH_MSGENULLP     "Null pointer"
#define LALCORRELATIONH_MSGESTART    "Time series do not start simultaneously"
#define LALCORRELATIONH_MSGESAMPLING "Time series are not sampled with the same rate"
/** \endcond */


/* ***** DEFINE OTHER GLOBAL CONSTANTS OR MACROS ************/

/* ***** DEFINE NEW STRUCTURES AND TYPES ************/

/** UNDOCUMENTED */
typedef struct
tagCorrelationParams
{
  REAL4 maxTimeShiftNan;
}
CorrelationParams;

/** UNDOCUMENTED */
typedef struct
tagInputCorrelation
{
  REAL4TimeSeries *one;
  REAL4TimeSeries *two;
}
InputCorrelation;

/** OutputCorrelation structure */
typedef struct
tagOutputCorrelation
{
  REAL4 *timeShiftedCorrelation; /**< cor(x(t-T),y(t)) or cor(x(t),y(t-T)) as a function of T*/

  INT4 maxCorrelationTimeShift;  /**< time step (from shift=0) at which correlation is maximum */
  REAL4 maxCorrelationValue;     /**< maximum value of the correlation for the range of shifts */

  INT4 minCorrelationTimeShift;	/**< UNDOCUMENTED */
  REAL4 minCorrelationValue;	/**< UNDOCUMENTED */

  LIGOTimeGPS start; 		/**< start time for the time series under consideration */
  INT4 length; 			/**< number of time steps in the time series */
  REAL8 deltaT; 		/**< time separation in seconds between different time steps */

  INT4 shift; 			/**< max time step shift between time series */
}
OutputCorrelation;

/* ***** INCLUDE EXTERNAL GLOBAL VARIABLES ************/

void
LALCorrelation( LALStatus                      *status,
		OutputCorrelation              **output,
		const InputCorrelation         *input,
		const CorrelationParams        *params);

/*@}*/

#ifdef  __cplusplus
}
#endif

#endif /* _LDASCAMPMOMENT_H */
