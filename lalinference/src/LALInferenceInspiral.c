/*
 *  LALInferenceInspiral.c:  CBC-specific initialization routines.
 *
 *  Copyright (C) 2015 Ben Farr
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

#include <stdio.h>
#include <lal/Date.h>
#include <lal/GenerateInspiral.h>
#include <lal/LALInference.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/StringInput.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/TimeSeries.h>
#include "LALInferenceMCMCSampler.h"
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALInferenceProposal.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceReadData.h>
#include <lal/LALInferenceInit.h>
#include <lalapps.h>
#include <lal/LALInferenceCalibrationErrors.h>

#include <mpi.h>

LALInferenceRunState *init_runstate(ProcessParamsTable *commandLine);
LALInferenceThreadState *init_threadstate(LALInferenceRunState *run_state);


/* Initialize a bare-bones run-state
 *  calls the "ReadData()" function to gather data & PSD from files,
 *  and initializes other variables accordingly.                     */
LALInferenceRunState *init_runstate(ProcessParamsTable *command_line) {
    LALInferenceRunState *run_state = XLALCalloc(1, sizeof(LALInferenceRunState));

    /* Check that command line is consistent first */
    LALInferenceCheckOptionsConsistency(command_line);
    run_state->commandLine = command_line;

    /* Read data from files or generate fake data */
    run_state->data = LALInferenceReadData(command_line);
    if (run_state->data == NULL)
        return(NULL);

    /* Perform injections if data successful read or created */
    LALInferenceInjectInspiralSignal(run_state->data, command_line);

    /* Apply calibration errors if desired*/
    LALInferenceApplyCalibrationErrors(run_state, command_line);

    /* Initialize parameters structure */
    run_state->algorithmParams = XLALCalloc(1, sizeof(LALInferenceVariables));
    run_state->priorArgs = XLALCalloc(1, sizeof(LALInferenceVariables));
    run_state->proposalArgs = XLALCalloc(1, sizeof(LALInferenceVariables));

    return(run_state);
}

LALInferenceThreadState *init_threadstate(LALInferenceRunState *run_state) {
    LALInferenceThreadState *thread;

    /* Set up CBC model and parameter array */
    thread = XLALCalloc(1, sizeof(LALInferenceThreadState));

    thread->currentPropDensity = -DBL_MAX;
    thread->model = LALInferenceInitCBCModel(run_state);

    /* Setup ROQ */
    LALInferenceSetupROQ(run_state->data, thread->model, command_line);

    thread->currentParams = XLALCalloc(1, sizeof(LALInferenceVariables));

    LALInferenceCopyVariables(thread->model->params, thread->currentParams);

    /* Inherit the proposal from the run state */
    thread->proposal = run_state->proposal;

    /* Explicitly zero out DE, in case it's not used later */
    thread->differentialPoints = NULL;
    thread->differentialPointsLength = 0;
    thread->differentialPointsSize = 0;

    /* If the currentParams are not in the prior, overwrite and pick paramaters from the priors. OVERWRITE EVEN USER CHOICES.
     *     (necessary for complicated prior shapes where LALInferenceCyclicReflectiveBound() is not enough */
    LALInferenceDrawApproxPrior(thread, thread->currentParams, thread->currentParams);
    while (run_state->prior(run_state, thread->currentParams, thread->model) <= -DBL_MAX)
        LALInferenceDrawApproxPrior(thread, thread->currentParams, thread->currentParams);

    /* Make sure that our initial value is within the prior-supported volume. */
    LALInferenceCyclicReflectiveBound(thread->currentParams, priorArgs);

    /* Initialize starting likelihood and prior */
    thread->currentPrior = run_state->prior(run_state, thread->currentParams, thread->model);

    thread->currentLikelihood = run_state->likelihood(run_state, thread->currentParams, thread->model);

    return thread;
}
