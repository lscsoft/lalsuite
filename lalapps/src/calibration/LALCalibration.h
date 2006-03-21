#ifndef _LALCALIBRATION_H
#define _LALCALIBRATION_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct tagLALCalData
{
  REAL4TimeSeries         *cavityFactors;
  REAL4TimeSeries         *openLoopFactors;
  COMPLEX8FrequencySeries *responseReference;
  COMPLEX8FrequencySeries *cavityGainReference;
  COMPLEX8FrequencySeries *openLoopGainReference;
}
LALCalData;

void XLALDestroyCalData( LALCalData *caldata );
int XLALUpdateResponse( 
    COMPLEX8FrequencySeries *response,  /**< response function to return */
    REAL8 duration,                     /**< duration for averaging factors */
    LALCalData *caldata                 /**< calibration reference data */
    );

#ifdef __cplusplus
}
#endif

#endif /* _LALCALIBRATION_H */
