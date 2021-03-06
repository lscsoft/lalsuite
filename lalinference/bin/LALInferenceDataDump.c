/*
 *  LALInferenceDataDump.c:  Dump LALInference Data
 *
 *  Copyright (C) 2015 John Veitch
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

#include <stdio.h>
#include <lal/LALInference.h>
#include <lal/LALInferenceInit.h>
#include <lal/LALInferenceReadData.h>
#include <lal/LALInferenceCalibrationErrors.h>

int main(int argc, char *argv[]){
  ProcessParamsTable *procParams = NULL;
  LALInferenceRunState *runState=NULL;

  procParams=LALInferenceParseCommandLine(argc,argv);

  runState = LALInferenceInitRunState(procParams);

  if(runState) {
    LALInferenceInjectInspiralSignal(runState->data,runState->commandLine);

    /* Simulate calibration errors */
    LALInferenceApplyCalibrationErrors(runState->data,runState->commandLine);
  }
}
