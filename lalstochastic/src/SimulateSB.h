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

/**
\author Sukanta Bose
\file
\ingroup stochastic

\brief Provides prototype and error code information for the modules needed
to simulate a stochastic background signal.

Provides prototype and error code information for the modules needed
to simulate a stochastic background signal (whitened, if desired) in a pair of
detectors, given the appropriate representations of the
detector transfer function in each detector.

\heading{Synopsis}
\code
#include <lal/SimulateSB.h>
\endcode

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

  NRCSID( SIMULATESBH, "$Id$" );

/**\name Error Codes */ /*@{*/
#define SIMULATESBH_ENULLP          1
#define SIMULATESBH_ENONPOSLEN      2
#define SIMULATESBH_ENONPOSDELTAF   3
#define SIMULATESBH_ENONPOSDELTAT   4
#define SIMULATESBH_ENEGFMIN        5
#define SIMULATESBH_EMMTIME         6
#define SIMULATESBH_EMMHETERO       7
#define SIMULATESBH_EMMFMIN         8
#define SIMULATESBH_EMMDELTAF       9
#define SIMULATESBH_EMMLEN         10
#define SIMULATESBH_EOORFREF       11
#define SIMULATESBH_ENONPOSOMEGA   12
#define SIMULATESBH_EALOC          13
#define SIMULATESBH_ENONZEROHETERO 14
#define SIMULATESBH_EWRONGUNITS    15
#define SIMULATESBH_ECOMPTIME      16
#define SIMULATESBH_ENOTYETHETERO 255

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
/*@}*/

  /*************************************************************
   *                                                           *
   *       Structures and prototypes associated with           *
   *             SimulateSB.c                  *
   *                                                           *
   *************************************************************/

/** Contains the output data produced by
<tt>LALSSSimStochBGTimeSeries()</tt>. It comprises of a pair of
(real) time-series simulated stochastic background signal in the outputs of
a given pair of detectors. The fields are:

<dl>
<dt><tt>REAL4TimeSeries *SSimStochBG1</tt></dt><dd>
Simulated stochastic background signal in the output of
the first detector.</dd>
<dt><tt>REAL4TimeSeries *SSimStochBG2</tt></dt><dd>
Simulated stochastic background signal in the output of
the second detector.</dd>
</dl>

*/
  typedef struct tagSSSimStochBGOutput {
    REAL4TimeSeries    *SSimStochBG1;
    REAL4TimeSeries    *SSimStochBG2;
  } SSSimStochBGOutput;

/** Contains the input data needed by
<tt>LALSSSimStochBGTimeSeries()</tt>
to calculate the whitened stochastic background signal in the output of
a detector.
The fields are:

<dl>
<dt><tt>REAL4FrequencySeries *omegaGW</tt></dt><dd> The spectrum
\f$\Omega_{\mathrm{GW}}(f)\f$ of the stochastic gravitational-wave
background.</dd>
<dt><tt>COMPLEX8FrequencySeries *whiteningFilter1</tt></dt><dd>
The frequency-domain response function \f$\tilde{R}_1(f)\f$ for the first detector.</dd>
<dt><tt>COMPLEX8FrequencySeries *whiteningFilter2</tt></dt><dd>
The frequency-domain response function \f$\tilde{R}_2(f)\f$ for the second detector.</dd>
</dl>

*/
  typedef struct tagSSSimStochBGInput {
    REAL4FrequencySeries     *omegaGW;
    COMPLEX8FrequencySeries  *whiteningFilter1;
    COMPLEX8FrequencySeries  *whiteningFilter2;
  } SSSimStochBGInput;

  typedef struct tagSSSimStochBGStrainInput {
    REAL4FrequencySeries     *omegaGW;
  } SSSimStochBGStrainInput;


/** Contains the parameters used by <tt>LALSSSimStochBGTimeSeries()</tt>
to compute the whitened stochastic background signal in the output of an
interferometric detector. The fields are:

<dl>
<dt><tt>UINT4 length</tt></dt><dd>
The number of points in the output time series.</dd>

<dt><tt>REAL8 deltaT</tt></dt><dd>
The temporal spacing of the output time series.</dd>

<dt><tt>INT4 seed</tt></dt><dd>
The random number seed for the stochastic simulation.</dd>

<dt><tt>LALDetector *detector1</tt></dt><dd>
The site location and orientation information of first detector involved in
the stochastic background search.</dd>

<dt><tt>LALDetector *detector2</tt></dt><dd>
The site location and orientation information of second detector involved in
the stochastic background search.</dd>

<dt><tt>LALUnit SSimStochBGTimeSeries1Unit</tt></dt><dd>
The unit field of the stochastic background, expressed as a Real4
time series, in detector 1.</dd>

<dt><tt>LALUnit SSimStochBGTimeSeries2Unit</tt></dt><dd>
The unit field of the stochastic background, expressed as a Real4
time series, in detector 2.
</dd>
</dl>
*/
  typedef struct tagSSSimStochBGParams {
    UINT4        length;   /* time length of output vector data samples */
    REAL8        deltaT;   /* time spacing */
    INT4         seed;     /* for random numbers x, y */
    LALDetector  detectorOne;
    LALDetector  detectorTwo;
    LALUnit      SSimStochBGTimeSeries1Unit;
    LALUnit      SSimStochBGTimeSeries2Unit;
  } SSSimStochBGParams;

  typedef struct tagSSSimStochBGStrainParams {
    UINT4        length1,length2;   /* time length of output vector data samples */
    REAL8        deltaT1, deltaT2;   /* time spacing */
    INT4         seed;     /* for random numbers x, y */
    LALDetector  detectorOne;
    LALDetector  detectorTwo;
    LALUnit      SSimStochBGTimeSeries1Unit;
    LALUnit      SSimStochBGTimeSeries2Unit;
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

#ifdef  __cplusplus
}
#endif

#endif /* _SIMULATESB_H */
