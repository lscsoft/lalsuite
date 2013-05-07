/*
 *
 *  LALInferenceInit.h:   Initialisation functions for LALInference codes        
 *
 *  Copyright (C) 2009 Vivien Raymond and John Veitch
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
 * \file LALInferenceInit.h
 * \brief Header file for initialisation functions used by LALInference codes
 *
 */

#ifndef LALInferenceInit_h
#define LALInferenceInit_h

#include <lal/LALInference.h>

/**
Initialise state variables needed for LALInferenceNest or LALInferenceMCMC to run
on a CBC signal. Reads the command line to get user-specified options
*/
LALInferenceVariables *LALInferenceInitCBCVariables(LALInferenceRunState *state);

/**
Initialise the template for a standard CBC signal
*/
void LALInferenceInitCBCTemplate(LALInferenceRunState *runState);

/** Review functions
 * 
 * */

LALInferenceVariables *LALInferenceInitVariablesReviewEvidence(LALInferenceRunState *state);
LALInferenceVariables *LALInferenceInitVariablesReviewEvidence_bimod(LALInferenceRunState *state);
LALInferenceVariables *LALInferenceInitVariablesReviewEvidence_banana(LALInferenceRunState *state);

#endif

