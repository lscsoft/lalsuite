//
//  LALInferenceFVectorMultiBanding.h
//  
//
//  Created by John Veitch and Serena Vinciguerra on 24/02/2015.
//
//

#ifndef _LALInferenceFVectorMultiBanding_Flat_h
#define _LALInferenceFVectorMultiBanding_Flat_h

/** Create a list of frequencies to use in multiband template generation, between fmin and fmax
 mc is minimum allowable chirp mass (sets freq evolution assumption ) */
REAL8Sequence *LALInferenceMultibandFrequencies(int NBands, double fmin, double fmax, double deltaF0, double mc);

#endif
