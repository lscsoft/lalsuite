/*-----------------------------------------------------------------------
 *
 * File Name: TFTransform.h
 *
 * Author: Eanna Flanagan
 *
 * Revision: $Id$ 
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * TFTransform.h
 *
 * SYNOPSIS
 * #include "TFTransform.h"
 *
 * DESCRIPTION
 * Function to create a time frequency (TF) plane from a segment of data
 * either (i) by stacking Fourier transforms of subsegments of the data
 * in the time domain, or (ii) FFTing the entire data segment and performing
 * inverse FFTs of subsegments of the data in the frequency domain.
 *
 * DIAGNOSTICS
 *
 *-----------------------------------------------------------------------
 */


#ifndef _TFTRANSFORM_H
#define _TFTRANSFORM_H


#include "LALStdlib.h"
#include "LALDatatypes.h"
#include "Window.h"
#include "RealFFT.h"
#include "ComplexFFT.h"
#include "LALRCSID.h"

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif


NRCSID (TFTRANSFORMH, "$Id$");


#define TFTRANSFORM_ENULLP       1
#define TFTRANSFORM_EPOSARG      2
#define TFTRANSFORM_EALLOCP      4
#define TFTRANSFORM_EMALLOC      8
#define TFTRANSFORM_EINCOMP      16





#define TFTRANSFORM_MSGENULLP      "Null pointer"
#define TFTRANSFORM_MSGEPOSARG     "Argument must be positive"
#define TFTRANSFORM_MSGEALLOCP     "Pointer has already been allocated, should be null"
#define TFTRANSFORM_MSGEMALLOC     "Malloc failure"
#define TFTRANSFORM_MSGEINCOMP     "Incompatible arguments"



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


typedef struct tagTFPlaneParams
{
  INT4             timeBins;    /* Number of time bins in TF plane    */
  INT4             freqBins;    /* Number of freq bins in TF plane    */
  REAL8            deltaT;      /* deltaF will always be 1/deltaT     */
  REAL8            flow;        
  /* 
   * frequencies f will lie in range flow <= f <= flow + freqBins * deltaF
   * flow is in Hertz
   */
}
TFPlaneParams;


typedef struct tagRealDFTParams
{
  WindowType               windowType;
  REAL4Vector              *window;
  REAL4                    sumofsquares;
  RealFFTPlan              *plan;
}
RealDFTParams;

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
}
VerticalTFTransformIn;


typedef struct tagHorizontalTFTransformIn
{
  ComplexDFTParams             *dftParams;
  INT4                         startT;
}
HorizontalTFTransformIn;




void
CreateRealDFTParams ( 
                     Status                         *status, 
                     RealDFTParams                  **dftParams, 
                     LALWindowParams                   *params,
                     INT2                           sign
                     );


void
DestroyRealDFTParams (
                      Status                        *status, 
                      RealDFTParams                 **dftParams
                      );


void
CreateComplexDFTParams ( 
                        Status                      *status, 
                        ComplexDFTParams            **dftParams, 
                        LALWindowParams                *params,
                        INT2                        sign
                        );


void
DestroyComplexDFTParams (
                         Status                     *status, 
                         ComplexDFTParams           **dftParams
                         );


void
ComputeFrequencySeries (
                        Status                      *status,
                        COMPLEX8FrequencySeries     *freqSeries,
                        REAL4TimeSeries             *timeSeries,
                        RealDFTParams               *dftParams
                        );


void
CreateTFPlane (
               Status                               *status,
               COMPLEX8TimeFrequencyPlane           **tfp,
               TFPlaneParams                        *input
               );


void
DestroyTFPlane (
               Status                               *status,
               COMPLEX8TimeFrequencyPlane           **tfp
                );


void
TimeSeriesToTFPlane (
                     Status                         *status,
                     COMPLEX8TimeFrequencyPlane     *tfp,
                     REAL4TimeSeries                *timeSeries,
                     VerticalTFTransformIn         *input
                     );


void
FreqSeriesToTFPlane (
                     Status                         *status,
                     COMPLEX8TimeFrequencyPlane     *tfp,
                     COMPLEX8FrequencySeries        *freqSeries,
                     HorizontalTFTransformIn        *input
                     );

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif





