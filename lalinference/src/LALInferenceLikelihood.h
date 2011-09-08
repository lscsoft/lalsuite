/*
 *
 *  LALInferenceLikelihood.c:   Likelihood functions for LALInference codes        
 *  LALInferenceLikelihood.h:   header file
 *
 *  Copyright (C) 2009 Ilya Mandel, Vivien Raymond, Christian Roever, Marc van der Sluys and John Veitch
 *
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
 * \file LALInferenceLikelihood.h
 * \brief Header file for likelihood functions used by LALInference codes
 *
 * LALInferenceLikelihood contains all the necessary routines to compute the likelihood 
 * from a template (computed with LALInferenceTemplate) and the data (initialised with LALInferenceReadData).
 *
 * Likelihood functions follow the basic naming convention: LALInference<type_of>LogLikelihood()
 * Takes as input:
 * - a pointer to a LALInferenceVariable structure containing the parameters to compute the likelihood for,
 * - a pointer to a LALInferenceIFOData structure containing the linked list of interferometer data,
 * - a pointer to the LALInferenceTemplateFunction template function to be used.
 * 
 * Outputs as a REAL8 the natural logarithm value of the likelihood, as defined by:
 *  
 * Likelihood(\vec{x}|\vec{\lambda},M)=\exp(-\tfrac{1}{2}<\vec{x}-\vec{h_M}(\vec{\lambda})|\vec{x}-\vec{h_M}(\vec{\lambda})>)
 * 
 * where: <x|y>=4Re\left ( \int \frac{\tilde{x}\,\tilde{y}^*}{S_f}\, df \right )
 * 
 * Note that the likelihood is reported unnormalised.
 * 
 */



#ifndef LALInferenceLikelihood_h
#define LALInferenceLikelihood_h

#include <lal/LALInference.h>

REAL8 LALInferenceUndecomposedFreqDomainLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData *data, 
                              LALInferenceTemplateFunction *template);

/* For testing purposes, likelihood that returns 0.0 = log(1) every
   time.  Activated with the --zeroLogLike command flag. */
REAL8 LALInferenceZeroLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData *data, LALInferenceTemplateFunction *template);

REAL8 LALInferenceFreqDomainLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData * data,
                              LALInferenceTemplateFunction *template);

REAL8 LALInferenceChiSquareTest(LALInferenceVariables *currentParams, LALInferenceIFOData * data,
                              LALInferenceTemplateFunction *template);
void LALInferenceComputeFreqDomainResponse(LALInferenceVariables *currentParams, LALInferenceIFOData * dataPtr,
                              LALInferenceTemplateFunction *template, COMPLEX16Vector *freqWaveform);

REAL8 LALInferenceTimeDomainLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData * data,
                              LALInferenceTemplateFunction *template);
void LALInferenceComputeTimeDomainResponse(LALInferenceVariables *currentParams, LALInferenceIFOData * dataPtr,
                               LALInferenceTemplateFunction *template, REAL8TimeSeries *timeWaveform);
REAL8 LALInferenceComputeFrequencyDomainOverlap(LALInferenceIFOData * data,
        COMPLEX16Vector * freqData1, COMPLEX16Vector * freqData2);

REAL8 LALInferenceNullLogLikelihood(LALInferenceIFOData *data);
REAL8 LALInferenceTimeDomainNullLogLikelihood(LALInferenceIFOData *data);


REAL8 LALInferenceWhitenedTimeDomainOverlap(const REAL8TimeSeries *whitenedData, const REAL8TimeSeries *data);
void LALInferencePSDToTDW(REAL8TimeSeries *TDW, const REAL8FrequencySeries *PSD, const REAL8FFTPlan *plan,
              const REAL8 fMin, const REAL8 fMax); 

/* UINT4 nextPowerOfTwo(const UINT4 n);
void padREAL8Sequence(REAL8Sequence *padded, const REAL8Sequence *data);
void padWrappedREAL8Sequence(REAL8Sequence *padded, const REAL8Sequence *data); */
/* UINT4 LIGOTimeGPSToNearestIndex(const LIGOTimeGPS *time, const REAL8TimeSeries *series); */
REAL8 LALInferenceIntegrateSeriesProduct(const REAL8TimeSeries *s1, const REAL8TimeSeries *s2);
void LALInferenceConvolveTimeSeries(REAL8TimeSeries *conv, const REAL8TimeSeries *data, const REAL8TimeSeries *response);
/* UINT4 NTDWFromNPSD(const UINT4 NPSD); */
void LALInferenceWrappedTimeSeriesToLinearTimeSeries(REAL8TimeSeries *linear, const REAL8TimeSeries *wrapped);
void LALInferenceLinearTimeSeriesToWrappedTimeSeries(REAL8TimeSeries *wrapped, const REAL8TimeSeries *linear);

REAL8 LALInferenceTimeDomainOverlap(const REAL8TimeSeries *TDW, const REAL8TimeSeries *A, const REAL8TimeSeries *B);

REAL8 LALInferenceFreqDomainStudentTLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData *data,
                                      LALInferenceTemplateFunction *template);

REAL8 LALInferenceAnalyticLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData *data, LALInferenceTemplateFunction *template);

#endif
