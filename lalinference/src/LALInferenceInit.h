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
 * Register a variable in vars for the model with given name, and a uniform prior.
 * Use the min and max arguments to specify a default range
 * Use startval to give a starting value. This is ignored unless the varytype==LALINFERENCE_PARAM_FIXED,
 * in which case it will be used as the initial value of the variable, unless the user over-rides it on
 * the command line (see below).
 * varytype is a LALInferenceParamVaryType e.g. LALINFERENCE_PARAM_LINEAR, LALINFERENCE_PARAM_CIRCULAR
 *
 * The function will query the state->commandLine to read command line arguments
 * of the form --name-min MINVAL --name-max MAXVAL to allow the user to over-ride
 * the prior limits. The initial value can be given with --name FIXEDVAL
 * If the --fix-name flag is given the variable will be pinned to that value
 * by setting it as LALINFERENCE_PARAM_FIXED
 * Note: The prior is setup in state->priorArgs, but it does not set state->currentParams
 */
void LALInferenceRegisterUniformVariableREAL8(LALInferenceRunState *state, LALInferenceVariables *var, const char name[VARNAME_MAX], REAL8 startval, REAL8 min, REAL8 max, LALInferenceParamVaryType varytype);


/**
 * Initialise state variables needed for LALInferenceNest or LALInferenceMCMC to run
 * on a CBC signal. Reads the command line to get user-specified options
 */
LALInferenceModel *LALInferenceInitCBCModel(LALInferenceRunState *state);

/**
 * Initialise the template for a standard CBC signal
 */
LALInferenceTemplateFunction LALInferenceInitCBCTemplate(LALInferenceRunState *runState);

/**
 Initialise the glitch fitting parameters
 */
void LALInferenceInitGlitchVariables(LALInferenceRunState *runState, LALInferenceVariables *currentParams);

/**
 * Review functions
 *
 * 
 */

LALInferenceModel *LALInferenceInitModelReviewEvidence(LALInferenceRunState *state);
LALInferenceModel *LALInferenceInitModelReviewEvidence_bimod(LALInferenceRunState *state);
LALInferenceModel *LALInferenceInitModelReviewEvidence_banana(LALInferenceRunState *state);

/**
 * Check options consistency 
 **/
void LALInferenceCheckOptionsConsistency(ProcessParamsTable *commandLine);


#endif

