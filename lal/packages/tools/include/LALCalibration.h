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
}
LALCalData;

/** Destroys a calibration data structure. */
void XLALDestroyCalData( LALCalData *caldata );

/** Constructs a response function at a requested frequency resolution
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
