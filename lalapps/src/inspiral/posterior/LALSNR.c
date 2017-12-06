/* 
 *  InferenceNest.c:  Nested Sampling using LALInference
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <lalapps.h>

#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>

#include <lal/Date.h>
#include <lal/GenerateInspiral.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/StringInput.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/TimeSeries.h>

#include <lal/LALInspiral.h>
#include <lal/LALSimulation.h>
#include <lal/LALCache.h>
#include <lal/LALFrStream.h>
#include <lal/TimeFreqFFT.h>
#include <lal/LALDetectors.h>
#include <lal/AVFactories.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/VectorOps.h>
#include <lal/Random.h>
#include <lal/LALNoiseModels.h>
#include <lal/XLALError.h>
#include <lal/GenerateInspiral.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/SeqFactories.h>
#include <lal/DetectorSite.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/Inject.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOMetadataInspiralUtils.h>
#include <lal/LIGOMetadataRingdownUtils.h>
#include <lal/LALInspiralBank.h>
#include <lal/FindChirp.h>
#include <lal/GenerateInspiral.h>
#include <lal/NRWaveInject.h>
#include <lal/GenerateInspRing.h>
#include <lal/LALErrno.h>

#include <lal/LALInference.h>
#include <lal/LALInferenceNestedSampler.h>
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceReadData.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALInferenceProposal.h>
#include <lal/LALInferenceInit.h>
#include <lal/LALSimNoise.h>
#include <lal/LALSimBlackHoleRingdownTiger.h>
//#include <LALInferenceRemoveLines.h>


LALInferenceRunState *initialize(ProcessParamsTable *commandLine); 

LALInferenceRunState *initialize(ProcessParamsTable *commandLine)
/* calls the "ReadData()" function to gather data & PSD from files, */
/* and initializes other variables accordingly.                     */
{
	char help[]="\
Initialisation arguments:\n\
(--verbose [N])\tOutput more info. N=1: errors, N=2 (default): warnings, N=3: info \n\n";
	LALInferenceRunState *irs=NULL;
	LALInferenceIFOData *ifoPtr, *ifoListStart;
	ProcessParamsTable *ppt=NULL;
	struct timeval tv;
	FILE *devrandom;
	
	irs = XLALCalloc(1, sizeof(LALInferenceRunState));
	/* read data from files: */
	fprintf(stdout, " readData(): started.\n");
	irs->commandLine=commandLine;
	
	INT4 verbose=0;
	INT4 x=0;
	ppt=LALInferenceGetProcParamVal(commandLine,"--verbose");
	if(ppt) {
	  if(ppt->value){
	    x=atoi(ppt->value);
	    switch(x){
	     case 0:
	       verbose=LALNDEBUG; /* Nothing */
	       break;
	     case 1:
	       verbose=LALMSGLVL1; /* Only errors */
	       break;
	     case 2:
	       verbose=LALMSGLVL2; /* Errors and warnings */
	       break;
	     case 3:
	       verbose=LALMSGLVL3; /* Errors, warnings and info */
	       break;
	     default:
	       verbose=LALMSGLVL2;
	       break;
	   }
	  }
	  else verbose=LALMSGLVL2; /* Errors and warnings */
	}
	
	irs->data = LALInferenceReadData(commandLine);
	/* (this will already initialise each LALIFOData's following elements:  */
  ppt=LALInferenceGetProcParamVal(commandLine,"--help");
  if(ppt)
  {
          fprintf(stdout,"%s",help);
          return(irs);
  }

	/*     fLow, fHigh, detector, timeToFreqFFTPlan, freqToTimeFFTPlan,     */
	/*     window, oneSidedNoisePowerSpectrum, timeDate, freqData         ) */
	fprintf(stdout, " LALInferenceReadData(): finished.\n");
	if (irs->data != NULL) {
		fprintf(stdout, " initialize(): successfully read data.\n");
		
		if (LALInferenceGetProcParamVal(commandLine,"--ringdown")){  
		    fprintf(stdout, " LALInferenceInjectRingdownSignal(): started.\n");
		    LALInferenceInjectRingdownSignal(irs->data, commandLine);
		    fprintf(stdout, " LALInferenceInjectRingdownSignal(): finished.\n");
}		else{
    		fprintf(stdout, " LALInferenceInjectInspiralSignal(): started.\n");
		    LALInferenceInjectInspiralSignal(irs->data,commandLine);
	    	fprintf(stdout, " LALInferenceInjectInspiralSignal(): finished.\n");
		}
		ifoPtr = irs->data;
		ifoListStart = irs->data;
		while (ifoPtr != NULL) {
			/*If two IFOs have the same sampling rate, they should have the same timeModelh*,
			 freqModelh*, and modelParams variables to avoid excess computation 
			 in model waveform generation in the future*/
			LALInferenceIFOData * ifoPtrCompare=ifoListStart;
			int foundIFOwithSameSampleRate=0;
			while(ifoPtrCompare != NULL && ifoPtrCompare!=ifoPtr) {
				if(ifoPtrCompare->timeData->deltaT == ifoPtr->timeData->deltaT){
					ifoPtr->timeModelhPlus=ifoPtrCompare->timeModelhPlus;
					ifoPtr->freqModelhPlus=ifoPtrCompare->freqModelhPlus;
					ifoPtr->timeModelhCross=ifoPtrCompare->timeModelhCross;				
					ifoPtr->freqModelhCross=ifoPtrCompare->freqModelhCross;				
					ifoPtr->modelParams=ifoPtrCompare->modelParams;	
					foundIFOwithSameSampleRate=1;	
					break;
				}
			}
			if(!foundIFOwithSameSampleRate){
				ifoPtr->timeModelhPlus  = XLALCreateREAL8TimeSeries("timeModelhPlus",&(ifoPtr->timeData->epoch),0.0,ifoPtr->timeData->deltaT,&lalDimensionlessUnit,ifoPtr->timeData->data->length);
				ifoPtr->timeModelhCross = XLALCreateREAL8TimeSeries("timeModelhCross",&(ifoPtr->timeData->epoch),0.0,	ifoPtr->timeData->deltaT,&lalDimensionlessUnit,ifoPtr->timeData->data->length);
				ifoPtr->freqModelhPlus = XLALCreateCOMPLEX16FrequencySeries("freqModelhPlus",&(ifoPtr->freqData->epoch),0.0,ifoPtr->freqData->deltaF,&lalDimensionlessUnit,ifoPtr->freqData->data->length);
				ifoPtr->freqModelhCross = XLALCreateCOMPLEX16FrequencySeries("freqModelhCross",&(ifoPtr->freqData->epoch), 0.0, ifoPtr->freqData->deltaF, &lalDimensionlessUnit,ifoPtr->freqData->data->length);
				ifoPtr->modelParams = XLALCalloc(1, sizeof(LALInferenceVariables));
			}
			ifoPtr = ifoPtr->next;
		}
	}
	else
	{
		fprintf(stdout, " initialize(): no data read.\n");
		exit(1);
	}

	return(irs);
}

/*************** MAIN **********************/


int main(int argc, char *argv[]){
        char help[]="\
LALSNR:\n\
Hacked from LALInferenceNest. Everything is removed, \n\
except what is required to output the SNR from xml file.\n\n";

  LALInferenceRunState *state;
  ProcessParamsTable *procParams=NULL;

  /* Read command line and parse */
  procParams=LALInferenceParseCommandLine(argc,argv);  
  
  /* This reads in the data */
  /* performs any injections specified */
  /* and calculates the snr */
  state = initialize(procParams);  
  
  /* Print command line arguments if help requested */
  if(LALInferenceGetProcParamVal(state->commandLine,"--help"))
  {
    fprintf(stdout,"%s",help);
    exit(0);
  }  

  return(0);
}


