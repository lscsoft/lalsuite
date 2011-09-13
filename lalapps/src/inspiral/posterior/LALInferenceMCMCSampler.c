/* 
 *  LALInferenceMCMC.c:  Bayesian Followup, MCMC algorithm.
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
#include <lal/LALInspiral.h>
#include <lal/DetResponse.h>
#include <lal/SeqFactories.h>
#include <lal/Date.h>
#include <lal/VectorOps.h>
#include <lal/TimeFreqFFT.h>
#include <lal/GenerateInspiral.h>
#include <lal/TimeDelay.h>
#include <lalapps.h>
#include <mpi.h>
#include <lal/LALInference.h>
//#include "mpi.h"
#include "LALInferenceMCMCSampler.h"
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALInferenceProposal.h>

#include <LALAppsVCSInfo.h>
#include <lal/LALStdlib.h>

RCSID("$Id$");
#define PROGRAM_NAME "LALInferenceMCMCSampler.c"
#define CVS_ID_STRING "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define CVS_NAME_STRING "$Name$"

int LALwaveformToSPINspiralwaveform(int waveform);

//Test LALAlgorithm
void PTMCMCAlgorithm(struct tagLALInferenceRunState *runState)
{

  const char *USAGE="\
(--appendOutput fname)          Basename of the file to append outputs to.\n";
  
  if(LALInferenceGetProcParamVal(runState->commandLine,"--help"))
	{
		fprintf(stdout,"%s",USAGE);
		return;
	}
  
  
	int i,t,p,lowerRank,upperRank,x; //indexes for for() loops
	REAL8 tempDelta;
	int nChain;
	int MPIrank, MPIsize;
	LALStatus status;
	memset(&status,0,sizeof(status));
	REAL8 dummyR8 = 0.0;
	REAL8 temperature = 1.0;
	INT4 acceptanceCount = 0;
	REAL8 nullLikelihood;
	REAL8 logChainSwap = 0.0;
	int tempIndex;
	//int *tempIndexVec = NULL;
	//int dummyTemp;
	REAL8 *tempLadder = NULL;			//the temperature ladder
	INT4 *acceptanceCountLadder = NULL;	//array of acceptance counts to compute the acceptance ratios.
	double *TcurrentLikelihood = NULL; //the current likelihood for each chain
  INT4 **pdf = NULL;
  REAL8 pdf_count = 0.0;

  //REAL8 *sigmaVec = NULL;
  REAL8Vector *sigmas = NULL;
  //REAL8 *PacceptCountVec = NULL;
  REAL8Vector *PacceptCount = NULL;
  //REAL8 *PproposeCountVec = NULL;
  REAL8Vector *PproposeCount = NULL;
  
  REAL8 *parametersVec = NULL;
  REAL8Vector * parameters = NULL;
  
	//LALVariables* TcurrentParams = malloc(sizeof(LALVariables));	//the current parameters for each chains
	//LALVariables dummyLALVariable;
	INT4 adaptationOn = 0;
  INT4 acceptanceRatioOn = 0;
	INT4 nPar = LALInferenceGetVariableDimensionNonFixed(runState->currentParams);
	INT4 Niter = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "Niter");
	INT4 Nskip = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "Nskip");
	REAL8 tempMax = *(REAL8*) LALInferenceGetVariable(runState->algorithmParams, "tempMax");   //max temperature in the temperature ladder
	UINT4 randomseed = *(UINT4*) LALInferenceGetVariable(runState->algorithmParams,"random_seed");

	UINT4 nIFO=0;
	LALInferenceIFOData *ifodata1=runState->data;
	ProcessParamsTable *ppt;//,*ppt1,*ppt2,*ppt3;	

	ppt=LALInferenceGetProcParamVal(runState->commandLine, "--help");
	if (ppt) {
	  fprintf(stderr, "%s\n", USAGE);
	}

	while(ifodata1){
		nIFO++;
		ifodata1=ifodata1->next;
	}
	

	MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);
	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

	nChain = MPIsize;		//number of parallel chain
	tempIndex = MPIrank;		//set initial temp indices

	tempLadder = malloc(nChain * sizeof(REAL8));			//the temperature ladder
	
    if(MPIrank == 0){
      parametersVec = (REAL8 *)malloc(MPIsize*nPar*sizeof(REAL8)); 
      for (p=0;p<(nChain*nPar);++p){
       parametersVec[p] = 0.0;
     }
   }
  
  parameters = XLALCreateREAL8Vector(nPar);
  
  LALInferenceVariableItem *ptr=runState->currentParams->head;
  p=0;
  while(ptr!=NULL) {
    if (ptr->vary != LALINFERENCE_PARAM_FIXED) {
      parameters->data[p]=*(REAL8 *)ptr->value;
      p++;
    }
    ptr=ptr->next;
  }
  
  
  ppt=LALInferenceGetProcParamVal(runState->commandLine, "--acceptanceRatio");
  if(ppt){
    acceptanceRatioOn = 1;

  }
  ppt=LALInferenceGetProcParamVal(runState->commandLine, "--adapt");
  if (ppt) {
    adaptationOn = 1;

  }

	
	
	if (nChain==1){
    ppt=LALInferenceGetProcParamVal(runState->commandLine,"--tempMax");
    if(ppt){
		  tempLadder[0]=tempMax;
    }else{
      tempLadder[0]=1.0;
		  tempMax=1.0;
    }
	}
	else {
		ppt = LALInferenceGetProcParamVal(runState->commandLine, "--inverseLadder");
		if(ppt){
			tempDelta = (1.0 - 1.0/tempMax)/(REAL8)(nChain-1);
			for (t=0; t<nChain; ++t) {
				tempLadder[t]=1.0/(REAL8)(1.0-t*tempDelta);
			}
		}
    else if(LALInferenceGetProcParamVal(runState->commandLine, "--geomLadder")){
      tempDelta=1.3;
      for (t=0;t<nChain; ++t) {
        tempLadder[t]=pow(tempDelta,t);
      }
    }
		else{
			tempDelta = log(tempMax)/(REAL8)(nChain-1);
			for (t=0; t<nChain; ++t) {
				tempLadder[t]=exp(t*tempDelta);
				}
			}
		
		}
  
  REAL8 s_gamma = 1.0;

  
  pdf=(INT4**)calloc(nPar,sizeof(INT4 *));
  for (p=0;p<nPar;++p){
    pdf[p]=calloc(100,sizeof(INT4));
    for(x=0;x<100;++x){
      pdf[p][x]=0;
    }
  }
  
  char *name = NULL;
  char nameMin[VARNAME_MAX], nameMax[VARNAME_MAX];
  REAL8 priorMin, priorMax, dprior;
  INT4 temperature_test = 0;
  if (LALInferenceGetProcParamVal(runState->commandLine, "--temperatureTest")) temperature_test=1;
  
  
	if (MPIrank == 0) {

		TcurrentLikelihood = (double*) malloc(sizeof(double)*nChain);
		acceptanceCountLadder = (int*) malloc(sizeof(int)*nChain);		//array of acceptance counts to compute the acceptance ratios.

		for (t=0; t<nChain; ++t) {
			acceptanceCountLadder[t] = 0;
			printf("tempLadder[%d]=%f\n",t,tempLadder[t]);
		}
	}
	
	
        if (runState->likelihood==&LALInferenceTimeDomainLogLikelihood) {
          fprintf(stderr, "Computing null likelihood in time domain.\n");
          nullLikelihood = LALInferenceTimeDomainNullLogLikelihood(runState->data);
        } else if (runState->likelihood==&LALInferenceUndecomposedFreqDomainLogLikelihood ||
                   runState->likelihood==&LALInferenceFreqDomainLogLikelihood) {
          nullLikelihood = LALInferenceNullLogLikelihood(runState->data);
        } else if (runState->likelihood==&LALInferenceZeroLogLikelihood) {
          nullLikelihood = 0.0;
        } else if (runState->likelihood==&LALInferenceAnalyticLogLikelihood) {
          nullLikelihood = 0.0;
        } else {
          fprintf(stderr, "Unrecognized log(L) function (in %s, line %d)\n", 
                  __FILE__, __LINE__);
          exit(1);
        }

	// initialize starting likelihood value:
	runState->currentLikelihood = runState->likelihood(runState->currentParams, runState->data, runState->template);
        LALInferenceIFOData *headData = runState->data;
        while (headData != NULL) {
          headData->acceptedloglikelihood = headData->loglikelihood;
          headData = headData->next;
        }
	runState->currentPrior = runState->prior(runState, runState->currentParams);
	
	
	//FILE **chainoutput = (FILE**)calloc(nChain,sizeof(FILE*));
  FILE * chainoutput = NULL;
	//char outfileName[99];

	//char **outfileName = (char**)calloc(nChain,sizeof(char*));
	char *outfileName = NULL;
  
  FILE *stat = NULL;
	char statfilename[256];
  if(MPIrank == 0){
    if (LALInferenceGetProcParamVal(runState->commandLine, "--adaptVerbose") || LALInferenceGetProcParamVal(runState->commandLine, "--acceptanceRatioVerbose")) {
      sprintf(statfilename,"PTMCMC.statistics.%u",randomseed);
      stat = fopen(statfilename, "a");
    }
  }
	//"waveform" and "pnorder" are ints to label the template used. Just to comform to SPINspiral output format. Should be temporary, and replaced by the command line used.
	int waveform = 0;
  if(LALInferenceCheckVariable(runState->currentParams,"LAL_APPROXIMANT")) waveform= LALwaveformToSPINspiralwaveform(*(INT4 *)LALInferenceGetVariable(runState->currentParams,"LAL_APPROXIMANT"));
  double pnorder = 0.0;
	if(LALInferenceCheckVariable(runState->currentParams,"LAL_PNORDER")) pnorder = ((double)(*(INT4 *)LALInferenceGetVariable(runState->currentParams,"LAL_PNORDER")))/2.0;

	char str[999];
	LALInferencePrintCommandLine(runState->commandLine, str);

  REAL8 networkSNR=0.0;
  ifodata1=runState->data;
  while(ifodata1){
    networkSNR+=ifodata1->SNR*ifodata1->SNR;
    ifodata1=ifodata1->next;
  }
  networkSNR=sqrt(networkSNR);
  
  REAL8 SampleRate=4096.0; //default value of the sample rate from LALInferenceReadData()
  
  if(LALInferenceGetProcParamVal(runState->commandLine,"--srate")) SampleRate=atof(LALInferenceGetProcParamVal(runState->commandLine,"--srate")->value);
  
	//for (t=0; t<nChain; ++t) {
		outfileName = (char*)calloc(99,sizeof(char*));
                ppt = LALInferenceGetProcParamVal(runState->commandLine, "--appendOutput");
                if (ppt) {
                  sprintf(outfileName, "%s.%2.2d", ppt->value, MPIrank);
                } else {
                  sprintf(outfileName,"PTMCMC.output.%u.%2.2d",randomseed,MPIrank);
                }
		if (!ppt) { /* Skip header output if we are appending. */
			chainoutput = fopen(outfileName,"w");
			fprintf(chainoutput, "  LALInference version:%s,%s,%s,%s,%s\n", LALAPPS_VCS_ID,LALAPPS_VCS_DATE,LALAPPS_VCS_BRANCH,LALAPPS_VCS_AUTHOR,LALAPPS_VCS_STATUS);
			fprintf(chainoutput,"  %s\n",str);
			fprintf(chainoutput, "%10s  %10s  %6s  %20s  %6s %8s   %6s  %8s  %10s  %12s  %9s  %9s  %8s\n",
					"nIter","Nburn","seed","null likelihood","Ndet","nCorr","nTemps","Tmax","Tchain","Network SNR","Waveform","pN order","Npar");
			fprintf(chainoutput, "%10d  %10d  %u  %20.10lf  %6d %8d   %6d%10d%12.1f%14.6f  %9i  %9.1f  %8i\n",
					Niter,0,randomseed,nullLikelihood,nIFO,0,nChain,(int)tempMax,tempLadder[MPIrank],networkSNR,waveform,(double)pnorder,nPar);
			fprintf(chainoutput, "\n%16s  %16s  %10s  %10s  %10s  %10s  %20s  %15s  %12s  %12s  %12s\n",
					"Detector","SNR","f_low","f_high","before tc","after tc","Sample start (GPS)","Sample length","Sample rate","Sample size","FT size");
			ifodata1=runState->data;
			while(ifodata1){
				fprintf(chainoutput, "%16s  %16.8lf  %10.2lf  %10.2lf  %10.2lf  %10.2lf  %20.8lf  %15.7lf  %.1f  %12d  %12d\n",
							ifodata1->detector->frDetector.name,ifodata1->SNR,ifodata1->fLow,ifodata1->fHigh,atof(LALInferenceGetProcParamVal(runState->commandLine,"--seglen")->value)-2.0,2.00,
							XLALGPSGetREAL8(&(ifodata1->epoch)),atof(LALInferenceGetProcParamVal(runState->commandLine,"--seglen")->value),SampleRate,
							(int)(atof(LALInferenceGetProcParamVal(runState->commandLine,"--seglen")->value)*SampleRate),
							(int)(atof(LALInferenceGetProcParamVal(runState->commandLine,"--seglen")->value)*SampleRate));
				ifodata1=ifodata1->next;
			}
			fprintf(chainoutput, "\n\n%31s\n","");
                        fprintf(chainoutput, "cycle\tlogpost\tlogprior\t");
                        LALInferenceFprintParameterNonFixedHeaders(chainoutput, runState->currentParams);
                        fprintf(chainoutput, "logl\t");
                        LALInferenceIFOData *headIFO = runState->data;
                        while (headIFO != NULL) {
                          fprintf(chainoutput, "logl");
                          fprintf(chainoutput, "%s",headIFO->name);
                          fprintf(chainoutput, "\t");
                          headIFO = headIFO->next;
                        }
			fprintf(chainoutput,"\n");
			fprintf(chainoutput, "%d\t%f\t%f\t", 0,(runState->currentLikelihood - nullLikelihood)+runState->currentPrior, runState->currentPrior);
			LALInferencePrintSampleNonFixed(chainoutput,runState->currentParams);
			fprintf(chainoutput,"%f\t",runState->currentLikelihood - nullLikelihood);
                        headIFO = runState->data;
                        while (headIFO != NULL) {
                          fprintf(chainoutput, "%f\t", headIFO->acceptedloglikelihood - headIFO->nullloglikelihood);
                          headIFO = headIFO->next;
                        }

			fprintf(chainoutput,"\n");
			fclose(chainoutput);
		}	

		chainoutput = fopen(outfileName,"a");

	INT4 parameter=0;

	LALInferenceAddVariable(runState->proposalArgs, "temperature", &temperature,  LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	LALInferenceAddVariable(runState->proposalArgs, "acceptanceCount", &acceptanceCount,  LALINFERENCE_INT4_t, LALINFERENCE_PARAM_LINEAR);

  if (adaptationOn == 1) {
    sigmas = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, SIGMAVECTORNAME)); 
  }
  if (acceptanceRatioOn == 1){
    PacceptCount = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, "PacceptCount"));
    PproposeCount = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, "PproposeCount"));
  }

  if (adaptationOn == 1) {
	LALInferenceAddVariable(runState->proposalArgs, "s_gamma", &s_gamma, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
  }
	LALInferenceAddVariable(runState->algorithmParams, "nChain", &nChain,  LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
	LALInferenceAddVariable(runState->algorithmParams, "nPar", &nPar,  LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
	LALInferenceAddVariable(runState->proposalArgs, "parameter",&parameter, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_LINEAR);
	LALInferenceAddVariable(runState->proposalArgs, "nullLikelihood", &nullLikelihood, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
	
	if (MPIrank == 0) {
		printf(" PTMCMCAlgorithm(); starting parameter values:\n");
		LALInferencePrintVariables(runState->currentParams);
		printf(" MCMC iteration: 0\t");
		printf("%f\t", runState->currentLikelihood - nullLikelihood); 
		printf("0\n");

                /* Do a data dump. */
                ppt = LALInferenceGetProcParamVal(runState->commandLine, "--data-dump");
                if (ppt) {
                  const UINT4 nameLength=256;
                  char filename[nameLength];
                  FILE *out;
                  //LALInferenceIFOData *headData = runState->data;
                  headData = runState->data;
                  UINT4 ui;
                  
                  while (headData != NULL) {

                    snprintf(filename, nameLength, "%s-freqModelhPlus.dat", headData->name);
                    out = fopen(filename, "w");
                    for (ui = 0; ui < headData->freqModelhPlus->data->length; ui++) {
                      REAL8 f = headData->freqModelhPlus->deltaF * ui;
                      COMPLEX16 d = headData->freqModelhPlus->data->data[ui];

                      fprintf(out, "%g %g %g\n", f, d.re, d.im);
                    }
                    fclose(out);

                    snprintf(filename, nameLength, "%s-freqModelhCross.dat", headData->name);
                    out = fopen(filename, "w");
                    for (ui = 0; ui < headData->freqModelhCross->data->length; ui++) {
                      REAL8 f = headData->freqModelhCross->deltaF * ui;
                      COMPLEX16 d = headData->freqModelhCross->data->data[ui];

                      fprintf(out, "%g %g %g\n", f, d.re, d.im);
                    }
                    fclose(out);

                    snprintf(filename, nameLength, "%s-timeModelhPlus.dat", headData->name);
                    out = fopen(filename, "w");
                    for (ui = 0; ui < headData->timeModelhPlus->data->length; ui++) {
                      REAL8 tt = XLALGPSGetREAL8(&(headData->timeModelhPlus->epoch)) + 
                        ui * headData->timeModelhPlus->deltaT;
                      REAL8 d = headData->timeModelhPlus->data->data[ui];

                      fprintf(out, "%.6f %g\n", tt, d);
                    }
                    fclose(out);

                    snprintf(filename, nameLength, "%s-timeModelhCross.dat", headData->name);
                    out = fopen(filename, "w");
                    for (ui = 0; ui < headData->timeModelhCross->data->length; ui++) {
                      REAL8 tt = XLALGPSGetREAL8(&(headData->timeModelhCross->epoch)) + 
                        ui * headData->timeModelhCross->deltaT;
                      REAL8 d = headData->timeModelhCross->data->data[ui];

                      fprintf(out, "%.6f %g\n", tt, d);
                    }
                    fclose(out);

                    snprintf(filename, nameLength, "%s-timeModel.dat", headData->name);
                    out = fopen(filename, "w");
                    for (ui = 0; ui < headData->timeModelhCross->data->length; ui++) {
                      REAL8 tt = XLALGPSGetREAL8(&(headData->timeModelhCross->epoch)) + 
                        headData->timeshift + ui*headData->timeModelhCross->deltaT;
                      REAL8 d = headData->fPlus*headData->timeModelhPlus->data->data[ui] + 
                        headData->fCross*headData->timeModelhCross->data->data[ui];

                      fprintf(out, "%.6f %g\n", tt, d);
                    }
                    fclose(out);

                    headData = headData->next;
                  }
                }

	}
	MPI_Barrier(MPI_COMM_WORLD);
	// iterate:
	for (i=1; i<=Niter; i++) {
    

		LALInferenceSetVariable(runState->proposalArgs, "temperature", &(tempLadder[MPIrank]));  //update temperature of the chain
		LALInferenceSetVariable(runState->proposalArgs, "acceptanceCount", &(acceptanceCount));

    if (adaptationOn == 1 && i == 1000001) {
        adaptationOn = 0;  //turn adaptation off after 10^6 iterations
    }
    if (adaptationOn == 1) {
      s_gamma=10.0*exp(-(1.0/6.0)*log((double)i))-1;
      LALInferenceSetVariable(runState->proposalArgs, "s_gamma", &(s_gamma));
    }    
    
		runState->evolve(runState); //evolve the chain with the parameters TcurrentParams[t] at temperature tempLadder[t]
		acceptanceCount = *(INT4*) LALInferenceGetVariable(runState->proposalArgs, "acceptanceCount");

		if ((i % Nskip) == 0) {
                  fseek(chainoutput, 0L, SEEK_END);
                  fprintf(chainoutput, "%d\t%f\t%f\t", i,(runState->currentLikelihood - nullLikelihood)+runState->currentPrior,runState->currentPrior);
                  LALInferencePrintSampleNonFixed(chainoutput,runState->currentParams);
                  fprintf(chainoutput,"%f\t",runState->currentLikelihood - nullLikelihood);
                  
                  LALInferenceIFOData *headIFO = runState->data;
                  while (headIFO != NULL) {
                    fprintf(chainoutput, "%f\t", headIFO->acceptedloglikelihood - headIFO->nullloglikelihood);
                    headIFO = headIFO->next;
                  }
		
                  fprintf(chainoutput,"\n");
                  fflush(chainoutput);

                  if (adaptationOn == 1) {
                    if (LALInferenceGetProcParamVal(runState->commandLine, "--adaptVerbose") || LALInferenceGetProcParamVal(runState->commandLine, "--acceptanceRatioVerbose")) {
                      if(MPIrank == 0){
            
                        fseek(stat, 0L, SEEK_END);
                        fprintf(stat,"%d\t",i);
                        
                        //printf("MPIrank=%d\tT=%f\ts_gamma=%f\t",MPIrank,*(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "temperature"),s_gamma);
                        if (LALInferenceGetProcParamVal(runState->commandLine, "--adaptVerbose")){
                          fprintf(stat,"%f\t",s_gamma);
                          for (p=0; p<nPar; ++p) {
                            fprintf(stat,"%f\t",sigmas->data[p]);
                          }
                        }
                        if(LALInferenceGetProcParamVal(runState->commandLine, "--acceptanceRatioVerbose")){
                          for (p=0; p<nPar; ++p) {
                            fprintf(stat,"%f\t",PacceptCount->data[p]/( PproposeCount->data[p]==0 ? 1.0 : PproposeCount->data[p] ));
                          }
                        }
                        fprintf(stat,"\n");
                        fflush(stat);
                      }
                    }
                  }
                }
		
    if ((i % Nskip) == 0) {
    ptr=runState->currentParams->head;
    p=0;
    while(ptr!=NULL) {
      if (ptr->vary != LALINFERENCE_PARAM_FIXED) {
        parameters->data[p]=*(REAL8 *)ptr->value;
        
        if(temperature_test==1){
          name = LALInferenceGetVariableName(runState->currentParams, (p+1));
          sprintf(nameMin, "%s_min", name);
          sprintf(nameMax, "%s_max", name);
          priorMin = *((REAL8 *)LALInferenceGetVariable(runState->priorArgs, nameMin));
          priorMax = *((REAL8 *)LALInferenceGetVariable(runState->priorArgs, nameMax));
          dprior = priorMax - priorMin;
          x=(int)(((parameters->data[p] - priorMin)/dprior)*100);
          if(x<0) x=0;
          if(x>99) x=99;
          pdf[p][x]++;
        }
        
        p++;
      }
      ptr=ptr->next;
    }    
  
    if(temperature_test==1){
      for (p=0;p<nPar;++p){
        pdf_count=0;
        for(x=0;x<100;++x){
          if(pdf[p][x]<((double)i)/1000.0) pdf_count++;
        }
        if(pdf_count==0) printf("PDF of parmeter %d is flat at temperature %f, iteration %d\n",p,tempLadder[MPIrank],i);
      }
    }

    
    dprior = priorMax - priorMin;
    

		MPI_Gather(&(runState->currentLikelihood), 1, MPI_DOUBLE, TcurrentLikelihood, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather(&acceptanceCount, 1, MPI_INT, acceptanceCountLadder, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(parameters->data,nPar,MPI_DOUBLE,parametersVec,nPar,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);	

		if (MPIrank == 0) {
			for(lowerRank=0;lowerRank<nChain-1;lowerRank++) { //swap parameters and likelihood between chains
				for(upperRank=lowerRank+1;upperRank<nChain;upperRank++) {
					
					logChainSwap = (1.0/tempLadder[lowerRank]-1.0/tempLadder[upperRank]) * (TcurrentLikelihood[upperRank]-TcurrentLikelihood[lowerRank]);
					
					if ((logChainSwap > 0)
						|| (log(gsl_rng_uniform(runState->GSLrandom)) < logChainSwap )) { //Then swap... 

            for (p=0; p<(nPar); ++p){
              dummyR8=parametersVec[p+nPar*upperRank];
              parametersVec[p+nPar*upperRank]=parametersVec[p+nPar*lowerRank];
              parametersVec[p+nPar*lowerRank]=dummyR8;
            }
            dummyR8 = TcurrentLikelihood[upperRank];
						TcurrentLikelihood[upperRank] = TcurrentLikelihood[lowerRank];
						TcurrentLikelihood[lowerRank] = dummyR8;

					}
				} //for(upperRank=lowerRank+1;upperRank<nChain;upperRank++)
			} //for(lowerRank=0;lowerRank<nChain-1;lowerRank++)
		} //if (MPIrank == 0)
      
    MPI_Scatter(parametersVec,nPar,MPI_DOUBLE,parameters->data,nPar,MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Scatter(TcurrentLikelihood, 1, MPI_DOUBLE, &(runState->currentLikelihood), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatter(acceptanceCountLadder, 1, MPI_INT, &acceptanceCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    ptr=runState->currentParams->head;
    p=0;
    while(ptr!=NULL) {
      if (ptr->vary != LALINFERENCE_PARAM_FIXED) {
        memcpy(ptr->value,&(parameters->data[p]),LALInferenceTypeSize[ptr->type]);
        p++;
      }
      ptr=ptr->next;
    }  
    
    MPI_Barrier(MPI_COMM_WORLD);
    }
	}// for (i=1; i<=Niter; i++)
  
	MPI_Barrier(MPI_COMM_WORLD);

  fclose(chainoutput);

  if(MPIrank == 0){
    if (LALInferenceGetProcParamVal(runState->commandLine, "--adaptVerbose") || LALInferenceGetProcParamVal(runState->commandLine, "--acceptanceRatioVerbose")) {
      fclose(stat);
    }
  }
	
	free(outfileName);
	free(tempLadder);
	
	if (MPIrank == 0) {
		free(TcurrentLikelihood);
		free(acceptanceCountLadder);
	}
}


void PTMCMCOneStep(LALInferenceRunState *runState)
// Metropolis-Hastings sampler.
{
  int MPIrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
  //INT4 p=0;
	REAL8 logPriorCurrent, logPriorProposed;
	REAL8 logLikelihoodCurrent, logLikelihoodProposed;
	LALInferenceVariables proposedParams;
	REAL8 logProposalRatio = 0.0;  // = log(P(backward)/P(forward))
	REAL8 logAcceptanceProbability;
	REAL8 temperature;
	INT4 acceptanceCount;
  INT4 accepted = 0;
  ProcessParamsTable *ppt, *commandLine = runState->commandLine;
	
	// current values:
	logPriorCurrent      = runState->currentPrior;
	logLikelihoodCurrent = runState->currentLikelihood;
  
	temperature = *(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "temperature");
	acceptanceCount = *(INT4*) LALInferenceGetVariable(runState->proposalArgs, "acceptanceCount");
	
	// generate proposal:
	proposedParams.head = NULL;
	proposedParams.dimension = 0;
	runState->proposal(runState, &proposedParams);
	if (LALInferenceCheckVariable(runState->proposalArgs, "logProposalRatio"))
		logProposalRatio = *(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "logProposalRatio");
	
	// compute prior & likelihood:
	logPriorProposed = runState->prior(runState, &proposedParams);
	if (logPriorProposed > -DBL_MAX)
		logLikelihoodProposed = runState->likelihood(&proposedParams, runState->data, runState->template);
	else
		logLikelihoodProposed = -DBL_MAX;
	
	// determine acceptance probability:
	logAcceptanceProbability = (1.0/temperature)*(logLikelihoodProposed - logLikelihoodCurrent) 
	+ (logPriorProposed - logPriorCurrent)
	+ logProposalRatio;
	
	// accept/reject:
	if ((logAcceptanceProbability > 0) 
		|| (log(gsl_rng_uniform(runState->GSLrandom)) < logAcceptanceProbability)) {   //accept
		LALInferenceCopyVariables(&proposedParams, runState->currentParams);
		runState->currentLikelihood = logLikelihoodProposed;
                LALInferenceIFOData *headData = runState->data;
                while (headData != NULL) {
                  headData->acceptedloglikelihood = headData->loglikelihood;
                  headData = headData->next;
                }
		runState->currentPrior = logPriorProposed;
		acceptanceCount++;
    accepted = 1;
		LALInferenceSetVariable(runState->proposalArgs, "acceptanceCount", &acceptanceCount);

	}
  
  INT4 adaptableStep = 0;
  INT4 i = 0;
  ppt = LALInferenceGetProcParamVal(commandLine, "--acceptanceRatio");
  if (ppt) {
    adaptableStep = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "adaptableStep"));
    if (adaptableStep) {
      i = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "proposedArrayNumber"));
      REAL8Vector *PacceptCount = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, "PacceptCount"));
      REAL8Vector *PproposeCount = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, "PproposeCount"));
      PproposeCount->data[i]+=1;
      if(accepted == 1){
        PacceptCount->data[i]+=1;
      }
    }
  }

        /* Adapt if desired. */
        ppt = LALInferenceGetProcParamVal(commandLine, "--adapt");
        if (ppt) {
          adaptableStep = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "adaptableStep"));
          if (adaptableStep) {
      
            i = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "proposedArrayNumber"));
            INT4 varNr = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "proposedVariableNumber"));
            //REAL8 Pacc = (logAcceptanceProbability > 0.0 ? 1.0 : exp(logAcceptanceProbability));
            REAL8 s_gamma = *(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "s_gamma");
            REAL8Vector *sigmas = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, SIGMAVECTORNAME));
            
            REAL8 sigma = sigmas->data[i];
            //REAL8 tau = *((REAL8 *)LALInferenceGetVariable(runState->proposalArgs, "adaptTau"));
            char *name = LALInferenceGetVariableName(&proposedParams, varNr);

            char nameMin[VARNAME_MAX], nameMax[VARNAME_MAX];
            REAL8 priorMin, priorMax, dprior;

            sprintf(nameMin, "%s_min", name);
            sprintf(nameMax, "%s_max", name);

            priorMin = *((REAL8 *)LALInferenceGetVariable(runState->priorArgs, nameMin));
            priorMax = *((REAL8 *)LALInferenceGetVariable(runState->priorArgs, nameMax));

            dprior = priorMax - priorMin;

            //sigma *= (1.0 - 4.0*Pacc - 8.0*tau)/(4.0*Pacc - 8.0*tau - 1.0);
            if(accepted == 1){
              sigma=sigma+s_gamma*(dprior/100.0)*(1.0-0.234);
            }else{
              sigma=sigma-s_gamma*(dprior/100.0)*(0.234);
            }
            
            sigma = (sigma > dprior ? dprior : sigma);
            sigma = (sigma < 0 ? 0 : sigma);
            
            sigmas->data[i] = sigma;

            //PacceptCount->data[i] *= 1.0 - 1.0/tau;
            //PacceptCount->data[i] += Pacc / tau;

            /* Make sure we don't do this again until we take another adaptable step.*/

          }
        }
  adaptableStep = 0;
  LALInferenceSetVariable(runState->proposalArgs, "adaptableStep", &adaptableStep);
	LALInferenceDestroyVariables(&proposedParams);	
}

void VNRPriorOneStep(LALInferenceRunState *runState)
// Von Neumann rejection sampler for the prior !!
{  
  LALInferenceVariables proposedParams;
	proposedParams.head = NULL;
	proposedParams.dimension = 0;
  
  PTMCMCLALInferenceDrawUniformlyFromPrior(runState, &proposedParams);
  runState->currentPrior = runState->prior(runState, &proposedParams);
  LALInferenceCopyVariables(&proposedParams, runState->currentParams);
  LALInferenceDestroyVariables(&proposedParams);
}



int LALwaveformToSPINspiralwaveform(int waveform)
{
	switch (waveform) {
		case 11:
			return 3;
			break;
    case 10:
			return 3;
			break;
		default:
			return 4;
			break;
	}
}


