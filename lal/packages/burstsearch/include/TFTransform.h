/********************************** <lalVerbatim file="TFTransformHV">
Author: Flanagan, E
$Id$
**************************************************** </lalVerbatim> */

#ifndef _TFTRANSFORM_H
#define _TFTRANSFORM_H


#include <lal/LALStdlib.h>
#include <lal/LALDatatypes.h>
#include <lal/Window.h>
#include <lal/RealFFT.h>
#include <lal/ComplexFFT.h>
#include <lal/LALRCSID.h>
#include <lal/TimeFreqFFT.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif


NRCSID (TFTRANSFORMH, "$Id$");


typedef enum {
  verticalPlane,
  /*
   *  Constructed by dividing time domain data into chunks
   *  and FFTing each chunk, thus producing a vertical line (the FFT)
   *  in the TF plane for each chunk in the time domain
   *
   */
  horizontalPlane
  /*
   *  Constructed by dividing frequency domain data into chunks
   *  and FFTing each chunk back to time domain, thus producing a horizontal
   *  line (the FFT) in the TF plane for each chunk in the frequency domain
   *
   */
}
TFPlaneType;

typedef enum {
  useSingleTFPlane, /* tiles are created in a single TF plane */
  useMultipleTFPlane /*tiles are created in multiple TF planes */ 
}
TFPlaneMethod;


typedef struct tagTFPlaneParams
{
  INT4             timeBins;    /* Number of time bins in TF plane    */
  INT4             freqBins;    /* Number of freq bins in TF plane    */
  REAL8            deltaT;      /* time resolution of the plane     */
  REAL8            deltaF;      /* freq. resolution of the plane */
  REAL8            flow;        /* minimum frequency to search for */
  REAL8            fhigh;       /* maximum frequency to search for */
  REAL8            timeDuration;/* length of data to be used to create a TF plane at a time */ 
}
TFPlaneParams;


typedef struct tagComplexDFTParams
{
  WindowType               windowType;
  REAL4Vector              *window;
  REAL4                    sumofsquares;
  ComplexFFTPlan           *plan;
}
ComplexDFTParams;


typedef struct tagCOMPLEX8TimeFrequencyPlane
{
  CHAR                     *name;
  LIGOTimeGPS              epoch;
  CHARVector               *sampleUnits;
  TFPlaneParams            *params;
  TFPlaneType              planeType;
  COMPLEX8                 *data;
  /*
   * data[i*params->freqBins+j] is a complex number
   * corresponding to a time t_i = epoch + i*(deltaT)
   * and a frequency f_j = flow + j / (deltaT)
   */
}
COMPLEX8TimeFrequencyPlane;


typedef struct tagVerticalTFTransformIn
{
  RealDFTParams                *dftParams;
  INT4                         startT;
  UINT4                        windowShift;
}
VerticalTFTransformIn;


typedef struct tagHorizontalTFTransformIn
{
  ComplexDFTParams             *dftParams;
  INT4                         startT;
  UINT4                        windowShift;
}
HorizontalTFTransformIn;




void
LALCreateComplexDFTParams (
                        LALStatus                      *status,
                        ComplexDFTParams            **dftParams,
                        LALWindowParams                *params,
                        INT2                        sign
                        );


void
LALDestroyComplexDFTParams (
                         LALStatus                     *status,
                         ComplexDFTParams           **dftParams
                         );


void
LALComputeFrequencySeries (
                        LALStatus                      *status,
                        COMPLEX8FrequencySeries     *freqSeries,
                        REAL4TimeSeries             *timeSeries,
                        RealDFTParams               *dftParams
                        );


void
LALCreateTFPlane (
               LALStatus                               *status,
               COMPLEX8TimeFrequencyPlane           **tfp,
               TFPlaneParams                        *input
               );


void
LALDestroyTFPlane (
               LALStatus                               *status,
               COMPLEX8TimeFrequencyPlane           **tfp
                );


void
LALTimeSeriesToTFPlane (
                     LALStatus                         *status,
                     COMPLEX8TimeFrequencyPlane     *tfp,
                     REAL4TimeSeries                *timeSeries,
                     VerticalTFTransformIn         *input
                     );


void
LALFreqSeriesToTFPlane (
                     LALStatus                         *status,
                     COMPLEX8TimeFrequencyPlane     *tfp,
                     COMPLEX8FrequencySeries        *freqSeries,
                     HorizontalTFTransformIn        *input
                     );
void
LALModFreqSeriesToTFPlane (
                     LALStatus                         *status,
                     COMPLEX8TimeFrequencyPlane     *tfp,
                     COMPLEX8FrequencySeries        *freqSeries,
                     HorizontalTFTransformIn        *input
                     );
void
LALModModFreqSeriesToTFPlane (
                     LALStatus                         *status,
                     COMPLEX8TimeFrequencyPlane     *tfp,
                     COMPLEX8FrequencySeries        *freqSeries,
                     HorizontalTFTransformIn        *input
                     );


#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif
