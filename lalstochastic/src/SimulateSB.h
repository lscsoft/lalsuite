/*
*  Copyright (C) 2007 Jolien Creighton, Tania Regimbau, Teviet Creighton, John Whelan
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

#ifndef _SIMULATESB_H
#define _SIMULATESB_H


#include <lal/LALStdlib.h>
#include <lal/DetectorSite.h>
#include <lal/Units.h>
#include <lal/RealFFT.h>

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \addtogroup SimulateSB_h
 * \author Sukanta Bose
 *
 * \brief Provides prototype and error code information for the modules needed
 * to simulate a stochastic background signal.
 *
 * Provides prototype and error code information for the modules needed
 * to simulate a stochastic background signal (whitened, if desired) in a pair of
 * detectors, given the appropriate representations of the
 * detector transfer function in each detector.
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/SimulateSB.h>
 * \endcode
 *
 */
/*@{*/

/**\name Error Codes */
/*@{*/
#define SIMULATESBH_ENULLP          1	/**< Null pointer */
#define SIMULATESBH_ENONPOSLEN      2	/**< Negative or zero length for data member of time series */
#define SIMULATESBH_ENONPOSDELTAF   3	/**< Negative or zero frequency spacing */
#define SIMULATESBH_ENONPOSDELTAT   4	/**< Negative or zero time spacing */
#define SIMULATESBH_ENEGFMIN        5	/**< Negative start frequency */
#define SIMULATESBH_EMMTIME         6	/**< Mismatch in epochs */
#define SIMULATESBH_EMMHETERO       7	/**< Mismatch in heterodyning frequencies */
#define SIMULATESBH_EMMFMIN         8	/**< Mismatch in start frequencies */
#define SIMULATESBH_EMMDELTAF       9	/**< Mismatch in frequency spacings */
#define SIMULATESBH_EMMLEN         10	/**< Mismatch in sequence lengths */
#define SIMULATESBH_EOORFREF       11	/**< Out of range reference frequency */
#define SIMULATESBH_ENONPOSOMEGA   12	/**< Negative stochastic background strength */
#define SIMULATESBH_EALOC          13	/**< Memory allocation error */
#define SIMULATESBH_ENONZEROHETERO 14	/**< Non-zero heterodyning frequency specified for real time series */
#define SIMULATESBH_EWRONGUNITS    15	/**< Inconsistent input units */
#define SIMULATESBH_ECOMPTIME      16	/**< Time domain data complex instead of real */
#define SIMULATESBH_ENOTYETHETERO 255	/**< Non-zero heterodyning frequency not yet implemented */
/*@}*/

/** \cond DONT_DOXYGEN */
#define SIMULATESBH_MSGENULLP         "Null pointer"
#define SIMULATESBH_MSGENONPOSLEN     "Negative or zero length for data member of time series"
#define SIMULATESBH_MSGENONPOSDELTAF  "Negative or zero frequency spacing"
#define SIMULATESBH_MSGENONPOSDELTAT  "Negative or zero time spacing"
#define SIMULATESBH_MSGENEGFMIN       "Negative start frequency"
#define SIMULATESBH_MSGEMMTIME        "Mismatch in epochs"
#define SIMULATESBH_MSGEMMHETERO      "Mismatch in heterodyning frequencies"
#define SIMULATESBH_MSGEMMFMIN        "Mismatch in start frequencies"
#define SIMULATESBH_MSGEMMDELTAF      "Mismatch in frequency spacings"
#define SIMULATESBH_MSGEMMLEN         "Mismatch in sequence lengths"
#define SIMULATESBH_MSGEOORFREF       "Out of range reference frequency"
#define SIMULATESBH_MSGENONPOSOMEGA   "Negative stochastic background strength"
#define SIMULATESBH_MSGEALOC         "Memory allocation error"
#define SIMULATESBH_MSGENONZEROHETERO "Non-zero heterodyning frequency specified for real time series"
#define SIMULATESBH_MSGEWRONGUNITS    "Inconsistent input units"
#define SIMULATESBH_MSGECOMPTIME      "Time domain data complex instead of real"
#define SIMULATESBH_MSGENOTYETHETERO  "Non-zero heterodyning frequency not yet implemented"
/** \endcond */

  /*************************************************************
   *                                                           *
   *       Structures and prototypes associated with           *
   *             SimulateSB.c                  *
   *                                                           *
   *************************************************************/

/**
 * Contains the output data produced by
 * <tt>LALSSSimStochBGTimeSeries()</tt>. It comprises of a pair of
 * (real) time-series simulated stochastic background signal in the outputs of
 * a given pair of detectors.
 */
  typedef struct tagSSSimStochBGOutput {
    REAL4TimeSeries    *SSimStochBG1;	/**< Simulated stochastic background signal in the output of the first detector */
    REAL4TimeSeries    *SSimStochBG2;	/**< Simulated stochastic background signal in the output of the second detector */
  } SSSimStochBGOutput;

/**
 * Contains the input data needed by <tt>LALSSSimStochBGTimeSeries()</tt>
 * to calculate the whitened stochastic background signal in the output of
 * a detector.
 */
  typedef struct tagSSSimStochBGInput {
    REAL4FrequencySeries     *omegaGW;		/**< The spectrum \f$\Omega_{\mathrm{GW}}(f)\f$ of the stochastic gravitational-wave background */
    COMPLEX8FrequencySeries  *whiteningFilter1;	/**< The frequency-domain response function \f$\tilde{R}_1(f)\f$ for the first detector */
    COMPLEX8FrequencySeries  *whiteningFilter2;	/**< The frequency-domain response function \f$\tilde{R}_2(f)\f$ for the second detector */
  } SSSimStochBGInput;

  /** UNDOCUMENTED */
  typedef struct tagSSSimStochBGStrainInput {
    REAL4FrequencySeries     *omegaGW;	/**< UNDOCUMENTED */
  } SSSimStochBGStrainInput;


/**
 * Contains the parameters used by <tt>LALSSSimStochBGTimeSeries()</tt>
 * to compute the whitened stochastic background signal in the output of an
 * interferometric detector.
 */
  typedef struct tagSSSimStochBGParams {
    UINT4        length;   			/**< The number of points in the output time series */
    REAL8        deltaT;			/**< The temporal spacing of the output time series */
    INT4         seed;				/**< The random number seed for the stochastic simulation */
    LALDetector  detectorOne;			/**< The site location and orientation information of first detector involved in the stochastic background search */
    LALDetector  detectorTwo;			/**< The site location and orientation information of second detector involved in the stochastic background search */
    LALUnit      SSimStochBGTimeSeries1Unit;	/**< The unit field of the stochastic background, expressed as a Real4 time series, in detector 1 */
    LALUnit      SSimStochBGTimeSeries2Unit;	/**< The unit field of the stochastic background, expressed as a Real4 time series, in detector 2 */
  } SSSimStochBGParams;

  /** UNDOCUMENTED */
  typedef struct tagSSSimStochBGStrainParams {
    UINT4        length1,length2;   		/**< time length of output vector data samples */
    REAL8        deltaT1, deltaT2;   		/**< time spacing */
    INT4         seed;     			/**< for random numbers x, y */
    LALDetector  detectorOne;		  	/**< UNDOCUMENTED */
    LALDetector  detectorTwo;			/**< UNDOCUMENTED */
    LALUnit      SSimStochBGTimeSeries1Unit;	/**< UNDOCUMENTED */
    LALUnit      SSimStochBGTimeSeries2Unit;	/**< UNDOCUMENTED */
  } SSSimStochBGStrainParams;


  void
  LALSSSimStochBGTimeSeries( LALStatus                  *status,
			     SSSimStochBGOutput           *output,
			     SSSimStochBGInput            *input,
			     SSSimStochBGParams           *params );

  void
  LALSSSimStochBGStrainTimeSeries( LALStatus                  *status,
			     SSSimStochBGOutput           *output,
			     SSSimStochBGStrainInput            *input,
			     SSSimStochBGStrainParams           *params );


/*@}*/

#ifdef  __cplusplus
}
#endif

#endif /* _SIMULATESB_H */
