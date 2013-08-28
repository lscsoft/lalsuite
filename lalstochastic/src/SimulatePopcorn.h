/*
*  Copyright (C) 2007 Bruce Allen, Jolien Creighton, Tania Regimbau
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

#ifndef _SIMULATEPOPCOR_H
#define _SIMULATEPOPCOR_H


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/AVFactories.h>
#include <lal/RealFFT.h>
#include <lal/ComplexFFT.h>
#include <lal/Units.h>
#include <lal/Random.h>
#include <lal/DetectorSite.h>


#ifdef __cplusplus
extern "C" {
#endif

/**
 * \addtogroup SimulatePopcorn_h
 * \author Tania Regimbau
 *
 * \brief Provides prototype for simulating whitened time-domain signals in a pair
 * of detectors that arises from low duty cycle astrophysical backgrounds.
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/SimulatePopcorn.h>
 * \endcode
 *
 */
/*@{*/

/** \name Error Codes */
/*@{*/
#define SIMULATEPOPCORNH_ENULLP          1	/**< Null pointer */
#define SIMULATEPOPCORNH_ENONNULLFMIN    2	/**< Non zero start frequency */
#define SIMULATEPOPCORNH_EMMDELTA        3	/**< Mismatch in sequence spacings */
#define SIMULATEPOPCORNH_EMMLEN          4	/**< Mismatch in sequence lengths */
#define SIMULATEPOPCORNH_EBV             5	/**< Bad input or parameter */
/*@}*/

/** \cond DONT_DOXYGEN */
#define SIMULATEPOPCORNH_MSGENULLP         "Null pointer"
#define SIMULATEPOPCORNH_MSGENONNULLFMIN   "Non zero start frequency"
#define SIMULATEPOPCORNH_MSGEMMDELTA       "Mismatch in sequence spacings"
#define SIMULATEPOPCORNH_MSGEMMLEN         "Mismatch in sequence lengths"
#define SIMULATEPOPCORNH_MSGEBV            "Bad input or parameter"
/** \endcond */

/** \name These constants define the cosmological model */
/*@{*/
#define SIMULATEPOPCORN_ho 0.7
#define SIMULATEPOPCORN_OMEGAMATTER 0.3
#define SIMULATEPOPCORN_OMEGAVACUUM 0.7
/*@}*/

/**
 * These are function pointers to functions that model burst waveforms.
 */
typedef void (REAL4LALWform) (REAL4 *output, REAL4 input);


/**
 * This structure contains the input of the simulation.
 */
typedef struct tagSimPopcornInputStruc {
  REAL4LALWform   *inputwform; 	/**< waveform of a single burst*/
  REAL4   inputduration; 	/**< mean duration of a single burst*/
  REAL4   inputlambda; 		/**< mean tims interval between successive bursts*/
  UINT4   inputNdataset; 	/**< number of detector sites 1 for H1/H2, 2 for H/L*/
  INT2   inputsite0; 		/**< first detector code*/
  INT2   inputsite1; 		/**< second detector code*/
  COMPLEX8FrequencySeries *wfilter0; /**< response of the first detector*/
  COMPLEX8FrequencySeries *wfilter1; /**< response of the second detector*/
} SimPopcornInputStruc;

/**
 * This structure contains the parameters of the simulation.
 */
typedef struct tagSimPopcornParamsStruc {
  UINT4   paramsstarttime; 	/**< starting time*/
  UINT4   paramslength; 	/**< length of the time serie in s*/
  UINT4   paramssrate; 		/**< sampling rate of the time serie in Hz*/
  UINT4   paramsseed; 		/**< random generator seed*/
  REAL8   paramsfref; 		/**< reference frequency if normalization, -1 otherwise*/
} SimPopcornParamsStruc;

/**
 * This structure contains the simulated pair time series and \f$\Omega\f$ spectrum
 */
typedef struct tagSimPopcornOutputStruc {
  REAL4TimeSeries   *SimPopcorn0;	/**< UNDOCUMENTED */
  REAL4TimeSeries   *SimPopcorn1;	/**< UNDOCUMENTED */
  REAL4FrequencySeries   *omega0;	/**< UNDOCUMENTED */
  REAL4FrequencySeries   *omega1;	/**< UNDOCUMENTED */
  } SimPopcornOutputStruc;


void
LALSimPopcornTimeSeries (LALStatus *status,SimPopcornOutputStruc *output,
    SimPopcornInputStruc *input,SimPopcornParamsStruc *params);


/*@}*/

#ifdef  __cplusplus
}
#endif
#endif
