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

#ifndef _LALCALIBRATION_H
#define _LALCALIBRATION_H

#include <lal/LALDatatypes.h>

NRCSID( LALCALIBRATIONH, "$Id$" );

#ifdef __cplusplus
extern "C" {
#endif

/** Structure containing calibration information (reference spectra and
 * factors for updating calibraiton to a particular time). */
typedef struct tagLALCalData
{
  REAL4TimeSeries         *cavityFactors;
  REAL4TimeSeries         *openLoopFactors;
  COMPLEX8FrequencySeries *responseReference;
  COMPLEX8FrequencySeries *cavityGainReference;
  COMPLEX8FrequencySeries *openLoopGainReference;
  COMPLEX8FrequencySeries *actuationReference;
  COMPLEX8FrequencySeries *digitalFilterReference;
}
LALCalData;

/** Destroys a calibration data structure. */
void XLALDestroyCalData( LALCalData *caldata );

/** Averages calibration factors alpha and gamma for a specified period */
int XLALAverageCalibrationFactors(
    REAL8 *cal_alpha,    /**< returned value of alpha factor */
    REAL8 *cal_gamma,    /**< returned value of gamma factor */
    LIGOTimeGPS *epoch,  /**< epoch to begin averaging */
    REAL8 duration,      /**< duration of averaging (0 = use just one point) */
    LALCalData *caldata  /**< calibration reference data */
    );

/** Returns the reference response function updated with given values
 * of alpha and gamma */
COMPLEX8FrequencySeries * XLALUpdateReferenceResponse(
    LALCalData *caldata, /**< calibration reference data */
    REAL8 cal_alpha,     /**< value of the alpha factor */
    REAL8 cal_gamma      /**< value of the gamma factor */
    );

/** Creates a response function at a requested frequency resolution
 * that is valid at a requested epoch. */
COMPLEX8FrequencySeries * XLALCreateCOMPLEX8Response(
    LIGOTimeGPS *epoch,    /**< epoch of response function */
    REAL8        duration, /**< duration of averaging (0 = use one point) */
    REAL8        deltaF,   /**< frequency resolution of requested response */
    UINT4        length,   /**< number of frequency bins in response */
    LALCalData  *caldata   /**< calibration reference data */
    );

/** Creates a response function at a requested frequency resolution
 * that is valid at a requested epoch. */
COMPLEX16FrequencySeries * XLALCreateCOMPLEX16Response(
    LIGOTimeGPS *epoch,    /**< epoch of response function */
    REAL8        duration, /**< duration of averaging (0 = use one point) */
    REAL8        deltaF,   /**< frequency resolution of requested response */
    UINT4        length,   /**< number of frequency bins in response */
    LALCalData  *caldata   /**< calibration reference data */
    );

/** Computes a response function at a requested frequency resolution
 * that is valid at a requested epoch. */
int XLALUpdateResponse(
    COMPLEX8FrequencySeries *response,  /**< response function to return */
    REAL8 duration,                     /**< duration for averaging factors */
    LALCalData *caldata                 /**< calibration reference data */
    );

#ifdef __cplusplus
}
#endif

#endif /* _LALCALIBRATION_H */
