//
//  LALInferenceFVectorMultiBanding.h
//  
//
//  Created by John Veitch and Serena Vinciguerra on 24/02/2015.
//
//

#ifndef _LALInferenceFVectorMultiBanding_Flat_h
#define _LALInferenceFVectorMultiBanding_Flat_h

/** Decimate a frequency series */
REAL8Sequence *DecimateREAL8Sequence(REAL8Sequence *old, UINT4 factor);
COMPLEX16Sequence *DecimateCOMPLEX16Sequence(COMPLEX16Sequence *old, UINT4 factor);
/* Populate the bands*/
//void LALInferencePopulateMultiBandData(LALInferenceIFOData *ifodata);
/** Create a list of frequencies */
REAL8Sequence *LALInferenceFrequencySequenceFunction(int NBands, double fmin, double fmax, double deltaF0, double mc);
/** F(t) and T(f) for newtonian waveform */
double LALInferenceTimeFrequencyRelation(double mc, double inPar, const char *type);

#endif
