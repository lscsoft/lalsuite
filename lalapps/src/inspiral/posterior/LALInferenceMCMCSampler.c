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
#include "LALInference.h"
//#include "mpi.h"
#include "LALInferenceMCMCSampler.h"
#include "LALInferencePrior.h"

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
	int i,t,p,lowerRank,upperRank; //indexes for for() loops
	int tempSwapCount=0;
	REAL8 tempDelta;
	int nChain;
	int count = 0;		//temporary counters to monitor the number of swaps between chains
	int MPIrank, MPIsize;
	LALStatus status;
	memset(&status,0,sizeof(status));
	//REAL8 dummyR8 = 0.0;
	REAL8 temperature = 1.0;
	INT4 acceptanceCount = 0;
	REAL8 nullLikelihood;
	REAL8 logChainSwap = 0.0;
	int tempIndex;
	int *tempIndexVec = NULL;
	int dummyTemp;
	REAL8 *tempLadder = NULL;			//the temperature ladder
	INT4 *acceptanceCountLadder = NULL;	//array of acceptance counts to compute the acceptance ratios.
	double *TcurrentLikelihood = NULL; //the current likelihood for each chain
	REAL8 *sigmaVec = NULL;
	//LALVariables* TcurrentParams = malloc(sizeof(LALVariables));	//the current parameters for each chains
	//LALVariables dummyLALVariable;
	
	INT4 nPar = getVariableDimensionNonFixed(runState->currentParams);
	INT4 Niter = *(INT4*) getVariable(runState->algorithmParams, "Niter");
	INT4 Nskip = *(INT4*) getVariable(runState->algorithmParams, "Nskip");
	REAL8 tempMax = *(REAL8*) getVariable(runState->algorithmParams, "tempMax");   //max temperature in the temperature ladder
	UINT4 randomseed = *(UINT4*) getVariable(runState->algorithmParams,"random_seed");
	UINT4 nIFO=0;
	LALIFOData *ifodata1=runState->data;
        const char *USAGE = 
          "[--appendOutput fname\tBasename of the file to append outputs to.]";
	ProcessParamsTable *ppt;//,*ppt1,*ppt2,*ppt3;	

	ppt=getProcParamVal(runState->commandLine, "--help");
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
		sigmaVec = (REAL8 *)malloc(MPIsize*nPar*sizeof(REAL8));//matrix of sigmas per parameter and temperature for adaptatio
		
		for (p=0;p<(nChain*nPar);++p){
			sigmaVec[p]=1.0;
		}				
	}
	//REAL8 *sigma = (REAL8*)calloc(nPar,sizeof(REAL8));
	REAL8 s_gamma = 1.0;
	REAL8Vector *sigma = XLALCreateREAL8Vector(nPar);
		for (p=0; p<nPar; ++p) {
			//sigma->data[p]=1.0;
			sigma->data[p]=0.1;
		}
	
	//REAL8 sigma = 0.1;
  //  gsl_matrix **sigma=calloc(1,sizeof(gsl_matrix *));
	//gsl_matrix * sigma = gsl_matrix_calloc(nChain,nPar);
//	if(NULL==(*sigma=gsl_matrix_alloc(nChain,nPar))) {fprintf(stderr,"Unable to allocate matrix memory\n"); exit(1);}
//	for (i = 0; i < nChain; i++){
//		for (j = 0; j < nPar; j++){
//			gsl_matrix_set (*sigma, i, j, i*j);
//		}
//	}
	
	
	if (nChain==1){
    ppt=getProcParamVal(runState->commandLine,"--tempMax");
    if(ppt){
		  tempLadder[0]=tempMax;
    }else{
      tempLadder[0]=1.0;
		  tempMax=1.0;
    }
	}
	else {
		ppt = getProcParamVal(runState->commandLine, "--inverseLadder");
		if(ppt){
			tempDelta = (1.0 - 1.0/tempMax)/(REAL8)(nChain-1);
			for (t=0; t<nChain; ++t) {
				tempLadder[t]=1.0/(REAL8)(1.0-t*tempDelta);
			}
		}
		else{
			tempDelta = log(tempMax)/(REAL8)(nChain-1);
			for (t=0; t<nChain; ++t) {
				tempLadder[t]=exp(t*tempDelta);
				}
			}
		
		}
	
	if (MPIrank == 0) {
		tempIndexVec = (int*) malloc(sizeof(int)*nChain);	//initialize temp index
		TcurrentLikelihood = (double*) malloc(sizeof(double)*nChain);
		acceptanceCountLadder = (int*) malloc(sizeof(int)*nChain);		//array of acceptance counts to compute the acceptance ratios.

		for (t=0; t<nChain; ++t) {
			tempIndexVec[t] = t;
			acceptanceCountLadder[t] = 0;
			printf("tempLadder[%d]=%f\n",t,tempLadder[t]);
		}
	}
	
	
        if (runState->likelihood==&TimeDomainLogLikelihood) {
          fprintf(stderr, "Computing null likelihood in time domain.\n");
          nullLikelihood = TimeDomainNullLogLikelihood(runState->data);
        } else if (runState->likelihood==&UndecomposedFreqDomainLogLikelihood ||
                   runState->likelihood==&FreqDomainLogLikelihood) {
          nullLikelihood = NullLogLikelihood(runState->data);
        } else if (runState->likelihood==&ZeroLogLikelihood) {
          nullLikelihood = 0.0;
        } else {
          fprintf(stderr, "Unrecognized log(L) function (in %s, line %d)\n", 
                  __FILE__, __LINE__);
          exit(1);
        }

	// initialize starting likelihood value:
	runState->currentLikelihood = runState->likelihood(runState->currentParams, runState->data, runState->template);
        LALIFOData *headData = runState->data;
        while (headData != NULL) {
          headData->acceptedloglikelihood = headData->loglikelihood;
          headData = headData->next;
        }
	runState->currentPrior = runState->prior(runState, runState->currentParams);
	
	
	FILE **chainoutput = (FILE**)calloc(nChain,sizeof(FILE*));
	//char outfileName[99];

	char **outfileName = (char**)calloc(nChain,sizeof(char*));
	
	//"waveform" and "pnorder" are ints to label the template used. Just to comform to SPINspiral output format. Should be temporary, and replaced by the command line used.
	
	int waveform= LALwaveformToSPINspiralwaveform(*(INT4 *)getVariable(runState->currentParams,"LAL_APPROXIMANT"));
	double pnorder = ((double)(*(INT4 *)getVariable(runState->currentParams,"LAL_PNORDER")))/2.0;
	char str[999];
	printCommandLine(runState->commandLine, str);

	for (t=0; t<nChain; ++t) {
		outfileName[t] = (char*)calloc(99,sizeof(char*));
                ppt = getProcParamVal(runState->commandLine, "--appendOutput");
                if (ppt) {
                  sprintf(outfileName[t], "%s.%2.2d", ppt->value, t);
                } else {
                  sprintf(outfileName[t],"PTMCMC.output.%u.%2.2d",randomseed,t);
                }
		if (MPIrank == 0 && !ppt) { /* Skip header output if we are appending. */
			chainoutput[t] = fopen(outfileName[t],"w");
			fprintf(chainoutput[t], "  LALInference version:%s,%s,%s,%s,%s\n", LALAPPS_VCS_ID,LALAPPS_VCS_DATE,LALAPPS_VCS_BRANCH,LALAPPS_VCS_AUTHOR,LALAPPS_VCS_STATUS);
			fprintf(chainoutput[t],"  %s\n",str);
			fprintf(chainoutput[t], "%10s  %10s  %6s  %20s  %6s %8s   %6s  %8s  %10s  %12s  %9s  %9s  %8s\n",
					"nIter","Nburn","seed","null likelihood","Ndet","nCorr","nTemps","Tmax","Tchain","Network SNR","Waveform","pN order","Npar");
			fprintf(chainoutput[t], "%10d  %10d  %u  %20.10lf  %6d %8d   %6d%10d%12.1f%14.6f  %9i  %9.1f  %8i\n",
					Niter,0,randomseed,nullLikelihood,nIFO,0,nChain,(int)tempMax,tempLadder[t],0.0,waveform,(double)pnorder,nPar);
			fprintf(chainoutput[t], "\n%16s  %16s  %10s  %10s  %10s  %10s  %20s  %15s  %12s  %12s  %12s\n",
					"Detector","SNR","f_low","f_high","before tc","after tc","Sample start (GPS)","Sample length","Sample rate","Sample size","FT size");
			ifodata1=runState->data;
			while(ifodata1){
				fprintf(chainoutput[t], "%16s  %16.8lf  %10.2lf  %10.2lf  %10.2lf  %10.2lf  %20.8lf  %15.7lf  %12d  %12d  %12d\n",
							ifodata1->detector->frDetector.name,0.0,ifodata1->fLow,ifodata1->fHigh,atof(getProcParamVal(runState->commandLine,"--seglen")->value)-2.0,2.00,
							XLALGPSGetREAL8(&(ifodata1->epoch)),atof(getProcParamVal(runState->commandLine,"--seglen")->value),atoi(getProcParamVal(runState->commandLine,"--srate")->value),
							(int)atof(getProcParamVal(runState->commandLine,"--seglen")->value)*atoi(getProcParamVal(runState->commandLine,"--srate")->value),
							(int)atof(getProcParamVal(runState->commandLine,"--seglen")->value)*atoi(getProcParamVal(runState->commandLine,"--srate")->value));
				ifodata1=ifodata1->next;
			}
			fprintf(chainoutput[t], "\n\n%31s\n","");
                        fprintf(chainoutput[t], "cycle\tlogpost\tlogprior\t");
                        fprintParameterNonFixedHeaders(chainoutput[t], runState->currentParams);
                        fprintf(chainoutput[t], "logl\t");
                        LALIFOData *headIFO = runState->data;
                        while (headIFO != NULL) {
                          fprintf(chainoutput[t], "logl");
                          fprintf(chainoutput[t], "%s",headIFO->name);
                          fprintf(chainoutput[t], "\t");
                          headIFO = headIFO->next;
                        }
			fprintf(chainoutput[t],"\n");
			fprintf(chainoutput[t], "%d\t%f\t%f\t", 0,(runState->currentLikelihood - nullLikelihood)+runState->currentPrior, runState->currentPrior);
			fprintSampleNonFixed(chainoutput[t],runState->currentParams);
			fprintf(chainoutput[t],"%f\t",runState->currentLikelihood - nullLikelihood);
                        headIFO = runState->data;
                        while (headIFO != NULL) {
                          fprintf(chainoutput[t], "%f\t", headIFO->acceptedloglikelihood - headIFO->nullloglikelihood);
                          headIFO = headIFO->next;
                        }
			//fprintf(chainoutput[t],"%f\t",tempLadder[t]);
			//fprintf(chainoutput[t],"%d\t",MPIrank);
			//fprintf(chainoutput[t],"%d\t",acceptanceCount);
			fprintf(chainoutput[t],"\n");
			fclose(chainoutput[t]);
		}	
		//fprintf(chainoutput[t],"This is temperature chain %d of %d.\n", t, nChain);
		//fclose(chainoutput[t]);
		chainoutput[t] = fopen(outfileName[t],"a");
	} 
	

	/*
	TcurrentParams.head=NULL;
	TcurrentParams.dimension=0;
	copyVariables(runState->currentParams,&(TcurrentParams));  //initiallize all chains
	*/
	
	/*
	dummyLALVariable.head=NULL;
	dummyLALVariable.dimension=0;
	copyVariables(runState->currentParams,&(dummyLALVariable));
	*/
	INT4 parameter=0;
	addVariable(runState->proposalArgs, "temperature", &temperature,  REAL8_t, PARAM_LINEAR);
	addVariable(runState->proposalArgs, "acceptanceCount", &acceptanceCount,  INT4_t, PARAM_LINEAR);
	//addVariable(runState->proposalArgs, "sigma", sigma,  gslMatrix_t, PARAM_LINEAR);
	addVariable(runState->proposalArgs, "sigma", &sigma,  REAL8Vector_t, PARAM_LINEAR);
	addVariable(runState->proposalArgs, "s_gamma", &s_gamma, REAL8_t, PARAM_LINEAR);
	addVariable(runState->algorithmParams, "nChain", &nChain,  INT4_t, PARAM_FIXED);
	addVariable(runState->algorithmParams, "nPar", &nPar,  INT4_t, PARAM_FIXED);
	addVariable(runState->proposalArgs, "parameter",&parameter, INT4_t, PARAM_LINEAR);
	addVariable(runState->proposalArgs, "tempIndex", &tempIndex,  INT4_t, PARAM_LINEAR);
	addVariable(runState->proposalArgs, "nullLikelihood", &nullLikelihood, REAL8_t, PARAM_FIXED);
	
	
	//for(t=0; t<nChain; t++) { TcurrentLikelihood[t] = runState->currentLikelihood; } // initialize the liklelihood for all chains
	
	if (MPIrank == 0) {
		printf(" PTMCMCAlgorithm(); starting parameter values:\n");
		printVariables(runState->currentParams);
		printf(" MCMC iteration: 0\t");
		printf("%f\t", runState->currentLikelihood - nullLikelihood); 
		printf("0\n");

                /* Do a data dump. */
                ppt = getProcParamVal(runState->commandLine, "--data-dump");
                if (ppt) {
                  const UINT4 nameLength=256;
                  char filename[nameLength];
                  FILE *out;
                  //LALIFOData *headData = runState->data;
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
		//printf(" MCMC iteration: %d\t", i+1);
		//copyVariables(&(TcurrentParams),runState->currentParams);
		setVariable(runState->proposalArgs, "temperature", &(tempLadder[tempIndex]));  //update temperature of the chain
		setVariable(runState->proposalArgs, "acceptanceCount", &(acceptanceCount));
		setVariable(runState->proposalArgs, "tempIndex", &(tempIndex));
		//s_gamma=10.0*exp(-(1.0/6.0)*log((double)i));
		//setVariable(runState->proposalArgs, "s_gamma", &(s_gamma));
		//setVariable(runState->proposalArgs, "sigma", sigmaVec[tempIndex]);
		//dummyR8 = runState->currentLikelihood;
		//	if (runState->currentLikelihood != dummyR8) {
		//		printf(" accepted! new parameter values:\n");
		//		printVariables(runState->currentParams);
		//	}
		//copyVariables(runState->currentParams,&(TcurrentParams));
		//TcurrentLikelihood[t] = runState->currentLikelihood; // save the parameters and temperature.
		runState->evolve(runState); //evolve the chain with the parameters TcurrentParams[t] at temperature tempLadder[t]
		acceptanceCount = *(INT4*) getVariable(runState->proposalArgs, "acceptanceCount");
		
                /* if (MPIrank == 0 && (i % Nskip == 0)) { */
                /*   /\* Every 100 steps, print out sigma. *\/ */
                /*   ppt = getProcParamVal(runState->commandLine, "--adapt"); */
                /*   if (ppt) { */
                /*     UINT4 j; */
                /*     REAL8Vector *sigmas = *((REAL8Vector **)getVariable(runState->proposalArgs, SIGMAVECTORNAME)); */
                /*     REAL8Vector *Pacc = *((REAL8Vector **) getVariable(runState->proposalArgs, "adaptPacceptAvg")); */
                /*     fprintf(stderr, "Iteration %10d: sigma = {", i); */
                /*     for (j = 0; j < sigmas->length-1; j++) { */
                /*       fprintf(stderr, "%8.5g, ", sigmas->data[j]); */
                /*     } */
                /*     fprintf(stderr, "%8.5g}\n", sigmas->data[sigmas->length - 1]); */
                /*     fprintf(stderr, "                    Paccept = {"); */
                /*     for (j = 0; j < Pacc->length-1; j++) { */
                /*       fprintf(stderr, "%8.5g, ", Pacc->data[j]); */
                /*     } */
                /*     fprintf(stderr, "%8.5g}\n", Pacc->data[Pacc->length - 1]); */
                /*   } */
                /* } */

		if ((i % Nskip) == 0) {
                  fseek(chainoutput[tempIndex], 0L, SEEK_END);
                  fprintf(chainoutput[tempIndex], "%d\t%f\t%f\t", i,(runState->currentLikelihood - nullLikelihood)+runState->currentPrior,runState->currentPrior);
                  fprintSampleNonFixed(chainoutput[tempIndex],runState->currentParams);
                  fprintf(chainoutput[tempIndex],"%f\t",runState->currentLikelihood - nullLikelihood);
                  
                  LALIFOData *headIFO = runState->data;
                  while (headIFO != NULL) {
                    fprintf(chainoutput[tempIndex], "%f\t", headIFO->acceptedloglikelihood - headIFO->nullloglikelihood);
                    headIFO = headIFO->next;
                  }
		
                  fprintf(chainoutput[tempIndex],"\n");
                  fflush(chainoutput[tempIndex]);
                }
		

		//if (tempIndex == 0) {
		/*	fprintf(stdout, "%8d %12.5lf %9.6lf", i,runState->currentLikelihood - nullLikelihood,1.0);
			
			fprintf(stdout," %9.5f",*(REAL8 *)getVariable(runState->currentParams,"x0"));*/

		/*	fprintf(stdout," %9.5f",*(REAL8 *)getVariable(&(TcurrentParams[t]),"chirpmass"));
			fprintf(stdout," %9.5f",*(REAL8 *)getVariable(&(TcurrentParams[t]),"massratio"));
			fprintf(stdout," %9.5f",*(REAL8 *)getVariable(&(TcurrentParams[t]),"inclination"));
			fprintf(stdout," %9.5f",*(REAL8 *)getVariable(&(TcurrentParams[t]),"phase"));
			fprintf(stdout," %9.5f",*(REAL8 *)getVariable(&(TcurrentParams[t]),"time"));
			fprintf(stdout," %9.5f",*(REAL8 *)getVariable(&(TcurrentParams[t]),"rightascension"));
			fprintf(stdout," %9.5f",*(REAL8 *)getVariable(&(TcurrentParams[t]),"declination"));
			fprintf(stdout," %9.5f",*(REAL8 *)getVariable(&(TcurrentParams[t]),"polarisation"));
			fprintf(stdout," %9.5f",*(REAL8 *)getVariable(&(TcurrentParams[t]),"distance"));*/
			
		//	fprintf(stdout,"\n");
		//}
		MPI_Gather(&(runState->currentLikelihood), 1, MPI_DOUBLE, TcurrentLikelihood, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather(&acceptanceCount, 1, MPI_INT, acceptanceCountLadder, 1, MPI_INT, 0, MPI_COMM_WORLD);
		//MPI_Gather(sigma->data,nPar,MPI_DOUBLE,sigmaVec,nPar,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);	
			
		//printVariables(&(TcurrentParams[0]));
		if (MPIrank == 0) {
			for(lowerRank=0;lowerRank<nChain-1;lowerRank++) { //swap parameters and likelihood between chains
				for(upperRank=lowerRank+1;upperRank<nChain;upperRank++) {
					
					logChainSwap = (1.0/tempLadder[tempIndexVec[lowerRank]]-1.0/tempLadder[tempIndexVec[upperRank]]) * (TcurrentLikelihood[upperRank]-TcurrentLikelihood[lowerRank]);
					
					if ((logChainSwap > 0)
						|| (log(gsl_rng_uniform(runState->GSLrandom)) < logChainSwap )) { //Then swap... 

						/*
						copyVariables(&(TcurrentParams[tempj]),&(dummyLALVariable));
						copyVariables(&(TcurrentParams[tempi]),&(TcurrentParams[tempj]));
						copyVariables(&(dummyLALVariable),&(TcurrentParams[tempi]));
						*/
						
						dummyTemp = tempIndexVec[upperRank];
						tempIndexVec[upperRank] = tempIndexVec[lowerRank];
						tempIndexVec[lowerRank] = dummyTemp;
						++tempSwapCount;
						
						dummyTemp = acceptanceCountLadder[upperRank];
						acceptanceCountLadder[upperRank]=acceptanceCountLadder[lowerRank];
						acceptanceCountLadder[lowerRank] = dummyTemp;
						/*
						dummyR8 = TcurrentLikelihood[tempj];
						TcurrentLikelihood[tempj] = TcurrentLikelihood[tempi];
						TcurrentLikelihood[tempi] = dummyR8;
						count++;
						*/
						//for (p=0; p<(nPar); ++p){
						//	dummyR8=sigmaVec[p+nPar*upperRank];
						//	sigmaVec[p+nPar*upperRank]=sigmaVec[p+nPar*lowerRank];
						//	sigmaVec[p+nPar*lowerRank]=dummyR8;
						//}
					}
				} //upperRank
			} //lowerRank
		} //MPIrank==0
		MPI_Scatter(tempIndexVec, 1, MPI_INT, &tempIndex, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Scatter(acceptanceCountLadder, 1, MPI_INT, &acceptanceCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
		//MPI_Scatter(sigmaVec,nPar,MPI_DOUBLE,sigma->data,nPar,MPI_DOUBLE,0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		//printf("%d\n",count);
		count = 0;
	}// for(i=0; i<100; i++)
	MPI_Barrier(MPI_COMM_WORLD);
	if (MPIrank == 0) printf("Temp swaps %d times\n", tempSwapCount);

	for (t=0; t<nChain; ++t) {
		fclose(chainoutput[t]);
	}
	
	free(chainoutput);

	for (t=0; t<nChain; ++t) {
		free(outfileName[t]);
	}
	
	free(outfileName);

	free(tempLadder);
	

	if (MPIrank == 0) {
		free(tempIndexVec);
		free(TcurrentLikelihood);
		free(acceptanceCountLadder);
	}
}


//Test LALEvolveOneStepFunction
void PTMCMCOneStep(LALInferenceRunState *runState)
// Metropolis-Hastings sampler.
{
	REAL8 logPriorCurrent, logPriorProposed;
	REAL8 logLikelihoodCurrent, logLikelihoodProposed;
	LALVariables proposedParams;
	REAL8 logProposalRatio = 0.0;  // = log(P(backward)/P(forward))
	REAL8 logAcceptanceProbability;
	REAL8 temperature;
	INT4 acceptanceCount;
        ProcessParamsTable *ppt, *commandLine = runState->commandLine;
	
	// current values:
	//logPriorCurrent      = runState->prior(runState, runState->currentParams);
	logPriorCurrent      = runState->currentPrior;
	logLikelihoodCurrent = runState->currentLikelihood;
	temperature = *(REAL8*) getVariable(runState->proposalArgs, "temperature");
	acceptanceCount = *(INT4*) getVariable(runState->proposalArgs, "acceptanceCount");
//	REAL8 nullLikelihood = *(REAL8*) getVariable(runState->proposalArgs, "nullLikelihood");
	
	// generate proposal:
	proposedParams.head = NULL;
	proposedParams.dimension = 0;
	runState->proposal(runState, &proposedParams);
	if (checkVariable(runState->proposalArgs, "logProposalRatio"))
		logProposalRatio = *(REAL8*) getVariable(runState->proposalArgs, "logProposalRatio");
	
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
	
	//fprintf(stdout," %9.5f=(1.0 / %9.5f)*(%9.5f-%9.5f)+(%9.5f-%9.5f)+%9.5f\n",logAcceptanceProbability,temperature,logLikelihoodProposed,logLikelihoodCurrent,logPriorProposed,logPriorCurrent,logProposalRatio);
	//double temp = log(gsl_rng_uniform(runState->GSLrandom));
	// accept/reject:
	if ((logAcceptanceProbability > 0) 
		|| (log(gsl_rng_uniform(runState->GSLrandom)) < logAcceptanceProbability)) {   //accept
		//|| (temp < logAcceptanceProbability)) {   //accept
		//if(logLikelihoodProposed>nullLikelihood){
		copyVariables(&proposedParams, runState->currentParams);
		runState->currentLikelihood = logLikelihoodProposed;
                LALIFOData *headData = runState->data;
                while (headData != NULL) {
                  headData->acceptedloglikelihood = headData->loglikelihood;
                  headData = headData->next;
                }
		runState->currentPrior = logPriorProposed;
		acceptanceCount++;
		setVariable(runState->proposalArgs, "acceptanceCount", &acceptanceCount);
		//}
	}

        /* Adapt if desired. */
        ppt = getProcParamVal(commandLine, "--adapt");
        if (ppt) {
          INT4 adaptableStep = *((INT4 *)getVariable(runState->proposalArgs, "adaptableStep"));
          if (adaptableStep) {
            INT4 i = *((INT4 *)getVariable(runState->proposalArgs, "proposedSigmaNumber"));
            INT4 varNr = *((INT4 *)getVariable(runState->proposalArgs, "proposedVariableNumber"));
            REAL8 Pacc = (logAcceptanceProbability > 0.0 ? 1.0 : exp(logAcceptanceProbability));
            REAL8Vector *sigmaVec = *((REAL8Vector **)getVariable(runState->proposalArgs, SIGMAVECTORNAME));
            REAL8Vector *avgPacc = *((REAL8Vector **)getVariable(runState->proposalArgs, "adaptPacceptAvg"));
            REAL8 sigma = sigmaVec->data[i];
            REAL8 tau = *((REAL8 *)getVariable(runState->proposalArgs, "adaptTau"));
            char *name = getVariableName(&proposedParams, varNr);
            char nameMin[VARNAME_MAX], nameMax[VARNAME_MAX];
            REAL8 priorMin, priorMax, dprior;

            sprintf(nameMin, "%s_min", name);
            sprintf(nameMax, "%s_max", name);

            priorMin = *((REAL8 *)getVariable(runState->priorArgs, nameMin));
            priorMax = *((REAL8 *)getVariable(runState->priorArgs, nameMax));

            dprior = priorMax - priorMin;

            sigma *= (1.0 - 4.0*Pacc - 8.0*tau)/(4.0*Pacc - 8.0*tau - 1.0);
            
            sigma = (sigma > dprior ? dprior : sigma);
            
            sigmaVec->data[i] = sigma;

            avgPacc->data[i] *= 1.0 - 1.0/tau;
            avgPacc->data[i] += Pacc / tau;

            /* Make sure we don't do this again until we take another adaptable step. */
            adaptableStep = 0;
            setVariable(runState->proposalArgs, "adaptableStep", &adaptableStep);
          }
        }
	//fprintf(stdout,"%9.5f < %9.5f\t(%9.5f)\n",temp,logAcceptanceProbability,temperature);
	destroyVariables(&proposedParams);	
}

//Test LALEvolveOneStepFunction
void PTMCMCAdaptationOneStep(LALInferenceRunState *runState)
// Metropolis-Hastings sampler.
{
	INT4 p=0;
	REAL8 logPriorCurrent, logPriorProposed;
	REAL8 logLikelihoodCurrent, logLikelihoodProposed;
	LALVariables proposedParams;
	REAL8 logProposalRatio = 0.0;  // = log(P(backward)/P(forward))
	REAL8 logAcceptanceProbability;
	REAL8 temperature;
	REAL8 *sigma=NULL;
	sigma=(REAL8 *)(*(REAL8Vector **)getVariable(runState->proposalArgs,"sigma"))->data;
	INT4 nPar  = *(INT4*) getVariable(runState->algorithmParams, "nPar");
	REAL8 s_gamma = *(REAL8*) getVariable(runState->proposalArgs, "s_gamma");
	gsl_rng * GSLrandom=runState->GSLrandom;
	
	// current values:
	logPriorCurrent      = runState->prior(runState, runState->currentParams);
	logLikelihoodCurrent = runState->currentLikelihood;
	temperature = *(REAL8*) getVariable(runState->proposalArgs, "temperature");
	//REAL8 nullLikelihood = *(REAL8*) getVariable(runState->proposalArgs, "nullLikelihood");
	
	// generate proposal:
	proposedParams.head = NULL;
	proposedParams.dimension = 0;
	
	
	if(gsl_ran_flat(GSLrandom,0.0,1.0)<=0.1){
	
	PTMCMCLALAdaptationProposal(runState, &proposedParams);
	if (checkVariable(runState->proposalArgs, "logProposalRatio"))
		logProposalRatio = *(REAL8*) getVariable(runState->proposalArgs, "logProposalRatio");
	
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
	
	//fprintf(stdout," %9.5f=(1.0 / %9.5f)*(%9.5f-%9.5f)+(%9.5f-%9.5f)+%9.5f\n",logAcceptanceProbability,temperature,logLikelihoodProposed,logLikelihoodCurrent,logPriorProposed,logPriorCurrent,logProposalRatio);
	//double temp = log(gsl_rng_uniform(runState->GSLrandom));
	// accept/reject:
	if ((logAcceptanceProbability > 0) 
		|| (log(gsl_rng_uniform(runState->GSLrandom)) < logAcceptanceProbability)) {   //accept
		//|| (temp < logAcceptanceProbability)) {   //accept
		//if(logLikelihoodProposed>nullLikelihood){
		copyVariables(&proposedParams, runState->currentParams);
		runState->currentLikelihood = logLikelihoodProposed;
                LALIFOData *headData = runState->data;
                while (headData != NULL) {
                  headData->acceptedloglikelihood = headData->loglikelihood;
                  headData=headData->next;
                }
		for (p=0; p<(nPar); ++p){
			sigma[p]=sigma[p]+s_gamma*(1.0-0.234);
		}
		//}
	}else{
		for (p=0; p<(nPar); ++p){
			sigma[p]=sigma[p]-s_gamma*(0.234);
			if(sigma[p]<0.0){sigma[p]=0.0;}
		}
	}
	
	}else {
		
		for (p=0; p<(nPar); ++p){
		setVariable(runState->proposalArgs, "parameter",&p);
			
		PTMCMCLALAdaptationSingleProposal(runState, &proposedParams);
		if (checkVariable(runState->proposalArgs, "logProposalRatio"))
			logProposalRatio = *(REAL8*) getVariable(runState->proposalArgs, "logProposalRatio");
		
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
		
		//fprintf(stdout," %9.5f=(1.0 / %9.5f)*(%9.5f-%9.5f)+(%9.5f-%9.5f)+%9.5f\n",logAcceptanceProbability,temperature,logLikelihoodProposed,logLikelihoodCurrent,logPriorProposed,logPriorCurrent,logProposalRatio);
		//double temp = log(gsl_rng_uniform(runState->GSLrandom));
		// accept/reject:
		if ((logAcceptanceProbability > 0) 
			|| (log(gsl_rng_uniform(runState->GSLrandom)) < logAcceptanceProbability)) {   //accept
			//|| (temp < logAcceptanceProbability)) {   //accept
			//if(logLikelihoodProposed>nullLikelihood){
			copyVariables(&proposedParams, runState->currentParams);
			runState->currentLikelihood = logLikelihoodProposed;
                        LALIFOData *headData = runState->data;
                        while (headData != NULL) {
                          headData->acceptedloglikelihood = headData->loglikelihood;
                          headData = headData->next;
                        }
							sigma[p]=sigma[p]+s_gamma*(1.0-0.234);
			//}
		}else{
				sigma[p]=sigma[p]-s_gamma*(0.234);
				if(sigma[p]<0.0){sigma[p]=0.0;}
			}
		}
		
	}
	//fprintf(stdout,"%9.5f < %9.5f\t(%9.5f)\n",temp,logAcceptanceProbability,temperature);
	//setVariable(runState->proposalArgs, "sigma", &sigma);
	
	destroyVariables(&proposedParams);	
}



//Test LALPriorFunction
//REAL8 PTUniformLALPrior(LALInferenceRunState *runState, LALVariables *params)
/****************************************/
/* Returns unnormalized (!),            */
/* logarithmic (!) prior density.      	*/
/****************************************/
/* Assumes the following parameters	*/
/* exist (e.g., for TaylorT1):		*/
/* chirpmass, massratio, inclination,	*/
/* phase, time, rightascension,		*/
/* desclination, polarisation, distance.*/
/* Prior is flat if within range	*/
/****************************************/
//{
//	REAL8 mc, eta, iota, phi, tc, ra, dec, psi, dist;	
//	REAL8 logdensity;
//	
//	mc   = *(REAL8*) getVariable(params, "chirpmass");		/* solar masses*/
//	eta  = *(REAL8*) getVariable(params, "massratio");		/* dim-less    */
//	iota = *(REAL8*) getVariable(params, "inclination");		/* radian      */
//	tc   = *(REAL8*) getVariable(params, "time");			/* GPS seconds */
//	phi  = *(REAL8*) getVariable(params, "phase");		/* radian      */
//	ra   = *(REAL8*) getVariable(params, "rightascension");	/* radian      */
//	dec  = *(REAL8*) getVariable(params, "declination");		/* radian      */
//	psi  = *(REAL8*) getVariable(params, "polarisation"); 	/* radian      */
//	dist = *(REAL8*) getVariable(params, "distance");		/* Mpc         */
//	
//	if(mc>2.41 && mc<=9.64 && eta>0.03 && eta<=0.25 && iota>=0.0 && iota<=LAL_PI && phi>=0.0 && phi<=LAL_TWOPI 
//	   && ra>=0.0 && ra<=LAL_TWOPI && dec>=-LAL_PI_2 && dec<=LAL_PI_2 && psi>=0.0 && psi<=LAL_PI && dist>0.0 && dist<=100.0
//	   && tc>=968654557.90 && tc<=968654558.20)	
//		logdensity = 0.0;
//	else
//		logdensity = -DBL_MAX;
//	//TODO: should be properly normalized; pass in range via priorArgs?	
//	
//	return(logdensity);
//}

static void PTMCMCCombinedProposal(LALInferenceRunState *runState, LALVariables *proposedParams,
                                   LALProposalFunction *props[], REAL8 weights[]) {
  REAL8 totalWeight, uniformRand;
  INT4 NProps, i;

  NProps=0;
  while (props[NProps] != NULL) NProps++;

  totalWeight = 0.0;
  for (i = 0; i < NProps; i++) {
    totalWeight += weights[i];
  }

  uniformRand = gsl_rng_uniform(runState->GSLrandom);

  i = 0;
  while (uniformRand > weights[i] / totalWeight) {
    uniformRand -= weights[i]/totalWeight;
    i++;
  }

  (props[i])(runState, proposedParams);

  return;
}

void PTMCMCLALProposal(LALInferenceRunState *runState, LALVariables *proposedParams)
{
	UINT4 nIFO=0;
	LALIFOData *ifo=runState->data;
	REAL8 BLOCKFRAC=0.0, /* Removed block jumps because of
                                roundoff problems in cholesky
                                decomp. */
          SINGLEFRAC=1.0,
          //SKYFRAC=0.0, /* Not symmetric! */
          INCFRAC=0.05,
          PHASEFRAC=0.05,
          SKYLOCSMALLWANDERFRAC=0.0; /* Not symmetric! Was: 0.05; */
        /* No spin rotations, because they are actually not symmetric! */
        REAL8 SPINROTFRAC = 0.0; /* (runState->template == &templateLALSTPN ? 0.05 : 0.0); */
        REAL8 COVEIGENFRAC;
        REAL8 IOTADISTANCEFRAC=0.0; /* Not symmetric! Stop! */
        REAL8 DIFFFULLFRAC;
        REAL8 DIFFPARTIALFRAC;
        ProcessParamsTable *ppt;
  
        ppt=getProcParamVal(runState->commandLine, "--iotaDistance");
        if (ppt) {
          IOTADISTANCEFRAC = atof(ppt->value);
        }
        ppt=getProcParamVal(runState->commandLine, "--covarianceMatrix");
        if (ppt) {
          COVEIGENFRAC = 1.0;
        } else {
          COVEIGENFRAC = 0.0;
        }

        if (getProcParamVal(runState->commandLine, "--differential-evolution")) {
          DIFFFULLFRAC = 1.0;
          DIFFPARTIALFRAC = 1.0 / 4.0;
        } else {
          DIFFFULLFRAC = 0.0;
          DIFFPARTIALFRAC = 0.0;
        }

        nIFO = 0;
        while(ifo){ifo=ifo->next; nIFO++;}
        
        if (nIFO < 2) {
          REAL8 weights[] = {BLOCKFRAC, SINGLEFRAC, INCFRAC, PHASEFRAC, SPINROTFRAC, COVEIGENFRAC, SKYLOCSMALLWANDERFRAC, IOTADISTANCEFRAC, DIFFFULLFRAC, DIFFPARTIALFRAC, DIFFPARTIALFRAC, DIFFPARTIALFRAC, DIFFPARTIALFRAC, DIFFPARTIALFRAC};
          LALProposalFunction *props[] = {&PTMCMCLALBlockCorrelatedProposal,
                                          &PTMCMCLALSingleAdaptProposal,
                                          &PTMCMCLALInferenceInclinationFlip,
                                          &PTMCMCLALInferenceOrbitalPhaseJump,
                                          &PTMCMCLALInferenceRotateSpins,
                                          &PTMCMCLALInferenceCovarianceEigenvectorJump,
                                          &PTMCMCLALInferenceSkyLocWanderJump,
                                          &PTMCMCLALInferenceInclinationDistanceConstAmplitudeJump,
                                          &PTMCMCLALInferenceDifferentialEvolutionFull,
                                          &PTMCMCLALInferenceDifferentialEvolutionMasses,
                                          &PTMCMCLALInferenceDifferentialEvolutionAmp,
                                          &PTMCMCLALInferenceDifferentialEvolutionSpins,
                                          &PTMCMCLALInferenceDifferentialEvolutionSky,
                                          0};
          PTMCMCCombinedProposal(runState, proposedParams, props, weights);
          return;
        } else if (nIFO < 3) {
          /* Removed the rotate sky function from proposal because it's not symmetric. */
          REAL8 weights[] = {BLOCKFRAC, SINGLEFRAC, 0.0 /* SKYFRAC */, INCFRAC, PHASEFRAC, SPINROTFRAC, COVEIGENFRAC, SKYLOCSMALLWANDERFRAC, IOTADISTANCEFRAC, DIFFFULLFRAC, DIFFPARTIALFRAC, DIFFPARTIALFRAC, DIFFPARTIALFRAC, DIFFPARTIALFRAC, DIFFPARTIALFRAC};
          LALProposalFunction *props[] = {&PTMCMCLALBlockCorrelatedProposal,
                                          &PTMCMCLALSingleAdaptProposal,
                                          &PTMCMCLALInferenceRotateSky,
                                          &PTMCMCLALInferenceInclinationFlip,
                                          &PTMCMCLALInferenceOrbitalPhaseJump,
                                          &PTMCMCLALInferenceRotateSpins,
                                          &PTMCMCLALInferenceCovarianceEigenvectorJump,
                                          &PTMCMCLALInferenceSkyLocWanderJump,
                                          &PTMCMCLALInferenceInclinationDistanceConstAmplitudeJump,
                                          &PTMCMCLALInferenceDifferentialEvolutionFull,
                                          &PTMCMCLALInferenceDifferentialEvolutionMasses,
                                          &PTMCMCLALInferenceDifferentialEvolutionAmp,
                                          &PTMCMCLALInferenceDifferentialEvolutionSpins,
                                          &PTMCMCLALInferenceDifferentialEvolutionSky,
                                          0};
          PTMCMCCombinedProposal(runState, proposedParams, props, weights);
        } else {
          /* Removed the rotate sky function because it's not symmetric. */
          REAL8 weights[] = {BLOCKFRAC, SINGLEFRAC, 0.0 /* SKYFRAC */, 0.0 /* SKYFRAC */, INCFRAC, PHASEFRAC, SPINROTFRAC, COVEIGENFRAC, SKYLOCSMALLWANDERFRAC, IOTADISTANCEFRAC, DIFFFULLFRAC, DIFFPARTIALFRAC, DIFFPARTIALFRAC, DIFFPARTIALFRAC, DIFFPARTIALFRAC, DIFFPARTIALFRAC};
          LALProposalFunction *props[] = {&PTMCMCLALBlockCorrelatedProposal,
                                          &PTMCMCLALSingleAdaptProposal,
                                          &PTMCMCLALInferenceRotateSky,
                                          (LALProposalFunction *)(&PTMCMCLALInferenceReflectDetPlane),
                                          &PTMCMCLALInferenceInclinationFlip,
                                          &PTMCMCLALInferenceOrbitalPhaseJump,
                                          &PTMCMCLALInferenceRotateSpins,
                                          &PTMCMCLALInferenceCovarianceEigenvectorJump,
                                          &PTMCMCLALInferenceSkyLocWanderJump,
                                          &PTMCMCLALInferenceInclinationDistanceConstAmplitudeJump,
                                          &PTMCMCLALInferenceDifferentialEvolutionFull,
                                          &PTMCMCLALInferenceDifferentialEvolutionMasses,
                                          &PTMCMCLALInferenceDifferentialEvolutionAmp,
                                          &PTMCMCLALInferenceDifferentialEvolutionSpins,
                                          &PTMCMCLALInferenceDifferentialEvolutionSky,
                                          0};
          PTMCMCCombinedProposal(runState, proposedParams, props, weights);
        }
}

void PTMCMCLALBlockProposal(LALInferenceRunState *runState, LALVariables *proposedParams)
{
	gsl_rng * GSLrandom=runState->GSLrandom;
	LALVariableItem *paraHead=NULL;
	INT4 i;
	copyVariables(runState->currentParams, proposedParams);

        REAL8 T = *(REAL8 *)getVariable(runState->proposalArgs, "temperature");
	
	REAL8 sigma = 0.1*sqrt(T); /* Adapt to temperature. */
	REAL8 big_sigma = 1.0;
	
	if(gsl_ran_ugaussian(GSLrandom) < 1.0e-3) big_sigma = 1.0e1;    //Every 1e3 iterations, take a 10x larger jump in all parameters
	if(gsl_ran_ugaussian(GSLrandom) < 1.0e-4) big_sigma = 1.0e2;    //Every 1e4 iterations, take a 100x larger jump in all parameters
	
	/* loop over all parameters */
	for (paraHead=proposedParams->head,i=0; paraHead; paraHead=paraHead->next)
	{ 
		if(paraHead->vary==PARAM_LINEAR || paraHead->vary==PARAM_CIRCULAR){
			
			if (!strcmp(paraHead->name,"massratio") || !strcmp(paraHead->name,"time") || !strcmp(paraHead->name,"a_spin2") || !strcmp(paraHead->name,"a_spin1")){
				*(REAL8 *)paraHead->value += gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.001;
			}else if (!strcmp(paraHead->name,"polarisation") || !strcmp(paraHead->name,"phase") || !strcmp(paraHead->name,"inclination")){
				*(REAL8 *)paraHead->value += gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.1;
			}else{
				*(REAL8 *)paraHead->value += gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;
			}
			i++;
		}
	}
	
	LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);
	
}

void PTMCMCLALSingleAdaptProposal(LALInferenceRunState *runState, LALVariables *proposedParams) {
  LALVariables *args = runState->proposalArgs;
  ProcessParamsTable *ppt = getProcParamVal(runState->commandLine, "--adapt");
  
  if (!checkVariable(args, SIGMAVECTORNAME) || !ppt) {
    /* We are not adaptive, or for some reason don't have a sigma
       vector---fall back on old proposal. */
    PTMCMCLALSingleProposal(runState, proposedParams);
  } else {
    gsl_rng *rng = runState->GSLrandom;
    LALVariableItem *param = NULL, *dummyParam = NULL;
    REAL8 T = *(REAL8 *)getVariable(args, "temperature");
    REAL8 sqrtT = sqrt(T);
    UINT4 dim;
    UINT4 i;
    UINT4 varNr;
    REAL8Vector *sigmas = *(REAL8Vector **) getVariable(args, SIGMAVECTORNAME);

    copyVariables(runState->currentParams, proposedParams);

    dim = proposedParams->dimension;

    do {
      varNr = 1+gsl_rng_uniform_int(rng, dim);
      param = getItemNr(proposedParams, varNr);
    } while (param->vary == PARAM_FIXED || param->vary == PARAM_OUTPUT);

    for (dummyParam = proposedParams->head, i = 0; dummyParam != NULL; dummyParam = dummyParam->next) {
      if (!strcmp(dummyParam->name, param->name)) {
        /* Found it; i = index into sigma vector. */
        break;
      } else if (dummyParam->vary == PARAM_FIXED || dummyParam->vary == PARAM_OUTPUT) {
        /* Don't increment i, since we're not dealing with a "real" parameter. */
        continue;
      } else {
        i++;
        continue;
      }
    }

    if (param->type != REAL8_t) {
      fprintf(stderr, "Attempting to set non-REAL8 parameter with numerical sigma (in %s, %d)\n",
              __FILE__, __LINE__);
      exit(1);
    } 

    if (i >= sigmas->length) {
      fprintf(stderr, "Attempting to draw single-parameter jump past the end of sigma array.\n(Maybe you used a non-spinning correlation matrix for a spinning run?)\nError in %s, line %d.\n",
              __FILE__, __LINE__);
      exit(1);
    }

    *((REAL8 *)param->value) += gsl_ran_ugaussian(rng)*sigmas->data[i]*sqrtT;

    LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);

    INT4 as = 1;
    setVariable(args, "adaptableStep", &as);

    setVariable(args, "proposedVariableNumber", &varNr);
    
    setVariable(args, "proposedSigmaNumber", &i);
  }
}

void PTMCMCLALSingleProposal(LALInferenceRunState *runState, LALVariables *proposedParams)
{
	gsl_rng * GSLrandom=runState->GSLrandom;
	LALVariableItem *paraHead=NULL;
	copyVariables(runState->currentParams, proposedParams);

        REAL8 T = *(REAL8 *)getVariable(runState->proposalArgs, "temperature");
	
	REAL8 sigma = 0.1*sqrt(T); /* Adapt step to temperature. */
	REAL8 big_sigma = 1.0;
	
	if(gsl_ran_ugaussian(GSLrandom) < 1.0e-3) big_sigma = 1.0e1;    //Every 1e3 iterations, take a 10x larger jump in a parameter
	if(gsl_ran_ugaussian(GSLrandom) < 1.0e-4) big_sigma = 1.0e2;    //Every 1e4 iterations, take a 100x larger jump in a parameter

	do {
          paraHead=getItemNr(proposedParams,1+(int)gsl_rng_uniform_int(GSLrandom,proposedParams->dimension));
	} while ((paraHead->vary==PARAM_FIXED || paraHead->vary==PARAM_OUTPUT));
	//printf("%s\n",paraHead->name);
		
        if (getProcParamVal(runState->commandLine, "--zeroLogLike")) {
          if (!strcmp(paraHead->name, "massratio")) {
            sigma = 0.02;
          } else if (!strcmp(paraHead->name, "chirpmass")) {
            sigma = 1.0;
          } else if (!strcmp(paraHead->name, "time")) {
            sigma = 0.02;
          } else if (!strcmp(paraHead->name, "phase")) {
            sigma = 0.6;
          } else if (!strcmp(paraHead->name, "distance")) {
            sigma = 10.0;
          } else if (!strcmp(paraHead->name, "declination")) {
            sigma = 0.3;
          } else if (!strcmp(paraHead->name, "rightascension")) {
            sigma = 0.6;
          } else if (!strcmp(paraHead->name, "polarisation")) {
            sigma = 0.6;
          } else if (!strcmp(paraHead->name, "inclination")) {
            sigma = 0.3;
          } else if (!strcmp(paraHead->name, "a_spin1")) {
            sigma = 0.1;
          } else if (!strcmp(paraHead->name, "theta_spin1")) {
            sigma = 0.3;
          } else if (!strcmp(paraHead->name, "phi_spin1")) {
            sigma = 0.6;
          } else if (!strcmp(paraHead->name, "a_spin2")) {
            sigma = 0.1;
          } else if (!strcmp(paraHead->name, "theta_spin2")) {
            sigma = 0.3;
          } else if (!strcmp(paraHead->name, "phi_spin2")) {
            sigma = 0.6;
          } else {
            fprintf(stderr, "Could not find parameter %s!", paraHead->name);
            exit(1);
          }
          *(REAL8 *)paraHead->value += gsl_ran_ugaussian(GSLrandom)*sigma;
        } else {
          if (!strcmp(paraHead->name,"massratio") || !strcmp(paraHead->name,"time") || !strcmp(paraHead->name,"a_spin2") || !strcmp(paraHead->name,"a_spin1")){
            *(REAL8 *)paraHead->value += gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.001;
          } else if (!strcmp(paraHead->name,"polarisation") || !strcmp(paraHead->name,"phase") || !strcmp(paraHead->name,"inclination")){
            *(REAL8 *)paraHead->value += gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.1;
          } else {
            *(REAL8 *)paraHead->value += gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;
          }
        }
	LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);
}

void PTMCMCLALBlockCorrelatedProposal(LALInferenceRunState *runState, LALVariables *proposedParams) {
  LALVariables *args = runState->proposalArgs;

  if (!checkVariable(args, COVMATRIXNAME) || !checkVariable(args, UNCORRSAMPNAME)) {
    /* No correlation matrix! */
    PTMCMCLALBlockProposal(runState, proposedParams);
  } else {
    gsl_rng *rng = runState->GSLrandom;
    REAL8 T = *(REAL8 *)getVariable(args, "temperature");
    REAL8 sqrtT = sqrt(T);
    gsl_matrix *covarianceMatrix = *(gsl_matrix **)getVariable(args, COVMATRIXNAME);
    REAL8Vector *uncorrelatedSample = *(REAL8Vector **)getVariable(args, UNCORRSAMPNAME);
    UINT4 i;
    UINT4 N = uncorrelatedSample->length;
    LALVariableItem *param = NULL;

    copyVariables(runState->currentParams, proposedParams);

    if (covarianceMatrix->size1 != N || covarianceMatrix->size2 != N) {
      fprintf(stderr, "ERROR: covariance matrix and sample vector sizes do not agree (in %s, line %d)\n",
              __FILE__, __LINE__);
      exit(1);
    }
    
    for (i = 0; i < N; i++) {
      uncorrelatedSample->data[i] = gsl_ran_ugaussian(rng)*sqrtT; /* Normalized to magnitude sqrt(T) */
    }

    for (i = 0, param = proposedParams->head; param != NULL; param = param->next) {
      if (param->vary != PARAM_FIXED && param->vary != PARAM_OUTPUT) {
        /* Then it's a parameter to set. */
        UINT4 j;
        REAL8 sum;

        for (j = 0, sum = 0.0; j < N; j++) {
          sum += gsl_matrix_get(covarianceMatrix, i, j)*uncorrelatedSample->data[j];
        }

        if (param->type != REAL8_t) {
          fprintf(stderr, "Trying to use covariance matrix to set non-REAL8 parameter (in %s, line %d)\n",
                  __FILE__, __LINE__);
          exit(1);
        } else {
          *((REAL8 *)param->value) += sum;
        }
        
        i++;
      }
    }

    LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);
  }
}

/* Reflect the inclination about the observing plane, iota -> Pi - iota */
void PTMCMCLALInferenceOrbitalPhaseJump(LALInferenceRunState *runState, LALVariables *proposedParams) {
  REAL8 phi;

  copyVariables(runState->currentParams, proposedParams);
  
  phi = *((REAL8 *) getVariable(proposedParams, "phase"));

  phi = fmod(phi+M_PI, 2.0*M_PI);

  setVariable(proposedParams, "phase", &phi);

  /* Probably not needed, but play it safe. */
  LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);
}


//Test LALProposalFunction
//void PTMCMCLALProposaltemp(LALInferenceRunState *runState, LALVariables *proposedParams)
/****************************************/
/* Assumes the following parameters		*/
/* exist (e.g., for TaylorT1):			*/
/* chirpmass, massratio, inclination,	*/
/* phase, time, rightascension,			*/
/* desclination, polarisation, distance.*/
/* Simply picks a new value based on	*/
/* fixed Gaussian;						*/
/* need smarter wall bounces in future.	*/
/****************************************/
//{
//	REAL8 mc, eta, iota, phi, tc, ra, dec, psi, dist;
//	REAL8 a_spin1, a_spin2, theta_spin1, theta_spin2, phi_spin1, phi_spin2;
//	REAL8 mc_proposed, eta_proposed, iota_proposed, phi_proposed, tc_proposed, 
//	ra_proposed, dec_proposed, psi_proposed, dist_proposed;
//	REAL8 a_spin1_proposed, a_spin2_proposed, theta_spin1_proposed, 
//	theta_spin2_proposed, phi_spin1_proposed, phi_spin2_proposed;
//	REAL8 logProposalRatio = 0.0;  // = log(P(backward)/P(forward))
//	gsl_rng * GSLrandom=runState->GSLrandom;
//	LALVariables * currentParams = runState->currentParams;
//	copyVariables(currentParams, proposedParams);
//	
//	REAL8 sigma = 0.1;
//	REAL8 big_sigma = 1.0;
//	
//	if(gsl_ran_ugaussian(GSLrandom) < 1.0e-3) big_sigma = 1.0e1;    //Every 1e3 iterations, take a 10x larger jump in all parameters
//	if(gsl_ran_ugaussian(GSLrandom) < 1.0e-4) big_sigma = 1.0e2;    //Every 1e4 iterations, take a 100x larger jump in all parameters
//	
//	
//	mc   = *(REAL8*) getVariable(currentParams, "chirpmass");		/* solar masses*/
//	eta  = *(REAL8*) getVariable(currentParams, "massratio");		/* dim-less    */
//	iota = *(REAL8*) getVariable(currentParams, "inclination");		/* radian      */
//	tc   = *(REAL8*) getVariable(currentParams, "time");				/* GPS seconds */
//	phi  = *(REAL8*) getVariable(currentParams, "phase");			/* radian      */
//	ra   = *(REAL8*) getVariable(currentParams, "rightascension");	/* radian      */
//	dec  = *(REAL8*) getVariable(currentParams, "declination");		/* radian      */
//	psi  = *(REAL8*) getVariable(currentParams, "polarisation");		/* radian      */
//	dist = *(REAL8*) getVariable(currentParams, "distance");			/* Mpc         */
//	
//	if (checkVariable(currentParams, "a_spin1")){
//		a_spin1 = *(REAL8*) getVariable(currentParams, "a_spin1");
//		a_spin1_proposed = a_spin1 + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.001;
//		setVariable(proposedParams, "a_spin1",      &a_spin1_proposed);
//	}
//	if (checkVariable(currentParams, "theta_spin1")){
//		theta_spin1 = *(REAL8*) getVariable(currentParams, "theta_spin1");
//		theta_spin1_proposed = theta_spin1 + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;
//		setVariable(proposedParams, "theta_spin1",      &theta_spin1_proposed);
//	}
//	if (checkVariable(currentParams, "phi_spin1")){
//		phi_spin1 = *(REAL8*) getVariable(currentParams, "phi_spin1");
//		phi_spin1_proposed = phi_spin1 + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;
//		setVariable(proposedParams, "phi_spin1",      &phi_spin1_proposed);
//	}
//	if (checkVariable(currentParams, "a_spin2")){
//		a_spin2 = *(REAL8*) getVariable(currentParams, "a_spin2");
//		a_spin2_proposed = a_spin2 + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.001;
//		setVariable(proposedParams, "a_spin2",      &a_spin2_proposed);
//	}
//	if (checkVariable(currentParams, "theta_spin2")){
//		theta_spin2 = *(REAL8*) getVariable(currentParams, "theta_spin2");
//		theta_spin2_proposed = theta_spin2 + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;
//		setVariable(proposedParams, "theta_spin2",      &theta_spin2_proposed);
//	}
//	if (checkVariable(currentParams, "phi_spin2")){
//		phi_spin2 = *(REAL8*) getVariable(currentParams, "phi_spin2");
//		phi_spin2_proposed = phi_spin2 + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;
//		setVariable(proposedParams, "phi_spin2",      &phi_spin2_proposed);
//	}
//
//	//mc_proposed   = mc*(1.0+gsl_ran_ugaussian(GSLrandom)*0.01);	/*mc changed by 1% */
//	// (above proposal is not symmetric!)
//	mc_proposed   = mc   + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;	/*mc changed by 0.0001 */
//	//mc_proposed   = mc * exp(gsl_ran_ugaussian(GSLrandom)*sigma*0.01);          /* mc changed by ~0.1% */
//	
//	eta_proposed  = eta  + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.001; /*eta changed by 0.01*/
//	//TODO: if(eta_proposed>0.25) eta_proposed=0.25-(eta_proposed-0.25); etc.
//	iota_proposed = iota + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.5;
//	tc_proposed   = tc   + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.001; /*time changed by 5 ms*/
//	phi_proposed  = phi  + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.1;
//	ra_proposed   = ra   + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;
//	dec_proposed  = dec  + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;
//	psi_proposed  = psi  + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.1;
//	//dist_proposed = dist + gsl_ran_ugaussian(GSLrandom)*0.5;
//	dist_proposed = dist * exp(gsl_ran_ugaussian(GSLrandom)*sigma*0.1); // ~10% change
//	
//	
//	
//	setVariable(proposedParams, "chirpmass",      &mc_proposed);		
//	setVariable(proposedParams, "massratio",      &eta_proposed);
//	setVariable(proposedParams, "inclination",    &iota_proposed);
//	setVariable(proposedParams, "phase",          &phi_proposed);
//	setVariable(proposedParams, "time",           &tc_proposed); 
//	setVariable(proposedParams, "rightascension", &ra_proposed);
//	setVariable(proposedParams, "declination",    &dec_proposed);
//	setVariable(proposedParams, "polarisation",   &psi_proposed);
//	setVariable(proposedParams, "distance",       &dist_proposed);
//	
//	LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);
//	
//	dist_proposed = *(REAL8*) getVariable(proposedParams, "distance");
//	logProposalRatio *= dist_proposed / dist;
//	//logProposalRatio *= mc_proposed / mc;   // (proposal ratio for above "scaled log-normal" proposal)
//	
//	// return ratio of proposal densities (for back & forth jumps) 
//	// in "runState->proposalArgs" vector:
//	if (checkVariable(runState->proposalArgs, "logProposalRatio"))
//		setVariable(runState->proposalArgs, "logProposalRatio", &logProposalRatio);
//	else
//		addVariable(runState->proposalArgs, "logProposalRatio", &logProposalRatio, REAL8_t, PARAM_OUTPUT);
//}



void PTMCMCLALAdaptationProposal(LALInferenceRunState *runState, LALVariables *proposedParams)
/****************************************/
/* Assumes the following parameters		*/
/* exist (e.g., for TaylorT1):			*/
/* chirpmass, massratio, inclination,	*/
/* phase, time, rightascension,			*/
/* desclination, polarisation, distance.*/
/* Simply picks a new value based on	*/
/* fixed Gaussian;						*/
/* need smarter wall bounces in future.	*/
/****************************************/
{
	REAL8 mc, eta, iota, phi, tc, ra, dec, psi, dist;
	REAL8 mc_proposed, eta_proposed, iota_proposed, phi_proposed, tc_proposed, 
	ra_proposed, dec_proposed, psi_proposed, dist_proposed;
	REAL8 logProposalRatio = 1.0;  // = log(P(backward)/P(forward))
	gsl_rng * GSLrandom=runState->GSLrandom;
	LALVariables * currentParams = runState->currentParams;
	//REAL8 sigmat = 0.1;
	//INT4 nPar = getVariableDimensionNonFixed(runState->currentParams);
	//INT4 i,j;
	
	//INT4 nPar  = *(INT4*) getVariable(runState->algorithmParams, "nPar");
	//INT4 nChain  = *(INT4*) getVariable(runState->algorithmParams, "nChain");
	//INT4 tempIndex  = *(INT4*) getVariable(runState->proposalArgs, "tempIndex");
	REAL8 big_sigma = 1.0;
	
	if(gsl_ran_ugaussian(GSLrandom) < 1.0e-3) big_sigma = 1.0e1;    //Every 1e3 iterations, take a 10x larger jump in all parameters
	if(gsl_ran_ugaussian(GSLrandom) < 1.0e-4) big_sigma = 1.0e2;    //Every 1e4 iterations, take a 100x larger jump in all parameters


	REAL8 *sigma=NULL;
	sigma=(REAL8 *)(*(REAL8Vector **)getVariable(runState->proposalArgs,"sigma"))->data;
	
	//gsl_matrix *sigma = *(gsl_matrix **)getVariable(runState->proposalArgs, "sigma");
	
	//printf ("m(%d,%d) = %g\n", 1, 1, gsl_matrix_get (sigma, 1, 1));
	
//	for (i = 0; i < nChain; i++){
//		for (j = 0; j < nPar; j++){
//			printf ("m(%d,%d) = %g\n", i, j, gsl_matrix_get (sigma, i, j));
//		}
//	}
//	printf("%f\n",sigma[0][0]);
//	int t,p;
//	for (t=0; t<5; ++t){
//		for (p=0; p<9; ++p){
		//	printf("sigma[%d][%d]=%f\n",t,p,sigma[t][p]);
//		}
//	}
	
	mc   = *(REAL8*) getVariable(currentParams, "chirpmass");		/* solar masses*/
	eta  = *(REAL8*) getVariable(currentParams, "massratio");		/* dim-less    */
	iota = *(REAL8*) getVariable(currentParams, "inclination");		/* radian      */
	tc   = *(REAL8*) getVariable(currentParams, "time");				/* GPS seconds */
	phi  = *(REAL8*) getVariable(currentParams, "phase");			/* radian      */
	ra   = *(REAL8*) getVariable(currentParams, "rightascension");	/* radian      */
	dec  = *(REAL8*) getVariable(currentParams, "declination");		/* radian      */
	psi  = *(REAL8*) getVariable(currentParams, "polarisation");		/* radian      */
	dist = *(REAL8*) getVariable(currentParams, "distance");			/* Mpc         */
	
	//mc_proposed   = mc*(1.0+gsl_ran_ugaussian(GSLrandom)*0.01);	/*mc changed by 1% */
	// (above proposal is not symmetric!)
	//mc_proposed   = mc   + gsl_ran_ugaussian(GSLrandom)*0.0001;	/*mc changed by 0.0001 */
	mc_proposed   = mc * exp(gsl_ran_gaussian(GSLrandom,sigma[8])*big_sigma*0.001);          /* mc changed by ~0.1% */
	logProposalRatio *= mc_proposed / mc;   // (proposal ratio for above "scaled log-normal" proposal)
	eta_proposed  = eta  + gsl_ran_gaussian(GSLrandom,sigma[7])*big_sigma*0.01; /*eta changed by 0.01*/
	//TODO: if(eta_proposed>0.25) eta_proposed=0.25-(eta_proposed-0.25); etc.
	iota_proposed = iota + gsl_ran_gaussian(GSLrandom,sigma[6])*big_sigma*0.1;
	tc_proposed   = tc   + gsl_ran_gaussian(GSLrandom,sigma[5])*big_sigma*0.005; /*time changed by 5 ms*/
	phi_proposed  = phi  + gsl_ran_gaussian(GSLrandom,sigma[4])*big_sigma*0.5;
	ra_proposed   = ra   + gsl_ran_gaussian(GSLrandom,sigma[3])*big_sigma*0.05;
	dec_proposed  = dec  + gsl_ran_gaussian(GSLrandom,sigma[2])*big_sigma*0.05;
	psi_proposed  = psi  + gsl_ran_gaussian(GSLrandom,sigma[1])*big_sigma*0.1;
	//dist_proposed = dist + gsl_ran_ugaussian(GSLrandom)*0.5;
	dist_proposed = dist * exp(gsl_ran_gaussian(GSLrandom,sigma[0])*big_sigma*0.1); // ~10% change
	logProposalRatio *= dist_proposed / dist;
	
	copyVariables(currentParams, proposedParams);
	setVariable(proposedParams, "chirpmass",      &mc_proposed);		
	setVariable(proposedParams, "massratio",      &eta_proposed);
	setVariable(proposedParams, "inclination",    &iota_proposed);
	setVariable(proposedParams, "phase",          &phi_proposed);
	setVariable(proposedParams, "time",           &tc_proposed); 
	setVariable(proposedParams, "rightascension", &ra_proposed);
	setVariable(proposedParams, "declination",    &dec_proposed);
	setVariable(proposedParams, "polarisation",   &psi_proposed);
	setVariable(proposedParams, "distance",       &dist_proposed);
	
	// return ratio of proposal densities (for back & forth jumps) 
	// in "runState->proposalArgs" vector:
	if (checkVariable(runState->proposalArgs, "logProposalRatio"))
		setVariable(runState->proposalArgs, "logProposalRatio", &logProposalRatio);
	else
		addVariable(runState->proposalArgs, "logProposalRatio", &logProposalRatio, REAL8_t, PARAM_OUTPUT);
}


void PTMCMCLALAdaptationSingleProposal(LALInferenceRunState *runState, LALVariables *proposedParams)
/****************************************/
/* Assumes the following parameters		*/
/* exist (e.g., for TaylorT1):			*/
/* chirpmass, massratio, inclination,	*/
/* phase, time, rightascension,			*/
/* desclination, polarisation, distance.*/
/* Simply picks a new value based on	*/
/* fixed Gaussian;						*/
/* need smarter wall bounces in future.	*/
/****************************************/
{
	REAL8 mc, eta, iota, phi, tc, ra, dec, psi, dist;
	REAL8 mc_proposed, eta_proposed, iota_proposed, phi_proposed, tc_proposed, 
	ra_proposed, dec_proposed, psi_proposed, dist_proposed;
	REAL8 logProposalRatio = 1.0;  // = log(P(backward)/P(forward))
	gsl_rng * GSLrandom=runState->GSLrandom;
	LALVariables * currentParams = runState->currentParams;
	//REAL8 sigmat = 0.1;
	//INT4 nPar = getVariableDimensionNonFixed(runState->currentParams);
	//INT4 i,j;
	
	INT4 p = *(INT4*) getVariable(runState->proposalArgs, "parameter");
	//INT4 nPar  = *(INT4*) getVariable(runState->algorithmParams, "nPar");
	//INT4 nChain  = *(INT4*) getVariable(runState->algorithmParams, "nChain");
	//INT4 tempIndex  = *(INT4*) getVariable(runState->proposalArgs, "tempIndex");
	REAL8 big_sigma = 1.0;
	
	if(gsl_ran_ugaussian(GSLrandom) < 1.0e-3) big_sigma = 1.0e1;    //Every 1e3 iterations, take a 10x larger jump in all parameters
	if(gsl_ran_ugaussian(GSLrandom) < 1.0e-4) big_sigma = 1.0e2;    //Every 1e4 iterations, take a 100x larger jump in all parameters
	
	
	REAL8 *sigma=NULL;
	sigma=(REAL8 *)(*(REAL8Vector **)getVariable(runState->proposalArgs,"sigma"))->data;
	
	//gsl_matrix *sigma = *(gsl_matrix **)getVariable(runState->proposalArgs, "sigma");
	
	//printf ("m(%d,%d) = %g\n", 1, 1, gsl_matrix_get (sigma, 1, 1));
	
	//	for (i = 0; i < nChain; i++){
	//		for (j = 0; j < nPar; j++){
	//			printf ("m(%d,%d) = %g\n", i, j, gsl_matrix_get (sigma, i, j));
	//		}
	//	}
	//	printf("%f\n",sigma[0][0]);
	//	int t,p;
	//	for (t=0; t<5; ++t){
	//		for (p=0; p<9; ++p){
	//	printf("sigma[%d][%d]=%f\n",t,p,sigma[t][p]);
	//		}
	//	}
	
	//mc_proposed   = mc*(1.0+gsl_ran_ugaussian(GSLrandom)*0.01);	/*mc changed by 1% */
	// (above proposal is not symmetric!)
	//mc_proposed   = mc   + gsl_ran_ugaussian(GSLrandom)*0.0001;	/*mc changed by 0.0001 */
		copyVariables(currentParams, proposedParams);
	
	switch ( p ) {
		case 8:
			mc   = *(REAL8*) getVariable(currentParams, "chirpmass");		/* solar masses*/
			mc_proposed   = mc * exp(gsl_ran_gaussian(GSLrandom,sigma[8])*big_sigma*0.001);          /* mc changed by ~0.1% */
			logProposalRatio *= mc_proposed / mc;   // (proposal ratio for above "scaled log-normal" proposal)
			setVariable(proposedParams, "chirpmass",      &mc_proposed);
			break;
		case 7:
			eta  = *(REAL8*) getVariable(currentParams, "massratio");		/* dim-less    */
			eta_proposed  = eta  + gsl_ran_gaussian(GSLrandom,sigma[7])*big_sigma*0.01; /*eta changed by 0.01*/
			//TODO: if(eta_proposed>0.25) eta_proposed=0.25-(eta_proposed-0.25); etc.
			setVariable(proposedParams, "massratio",      &eta_proposed);
			break;
		case 6:
			iota = *(REAL8*) getVariable(currentParams, "inclination");		/* radian      */
			iota_proposed = iota + gsl_ran_gaussian(GSLrandom,sigma[6])*big_sigma*0.1;
			setVariable(proposedParams, "inclination",    &iota_proposed);
			break;
		case 5:
			tc   = *(REAL8*) getVariable(currentParams, "time");				/* GPS seconds */
			tc_proposed   = tc   + gsl_ran_gaussian(GSLrandom,sigma[5])*big_sigma*0.005; /*time changed by 5 ms*/
			setVariable(proposedParams, "time",           &tc_proposed);
			break;
		case 4:
			phi  = *(REAL8*) getVariable(currentParams, "phase");			/* radian      */
			phi_proposed  = phi  + gsl_ran_gaussian(GSLrandom,sigma[4])*big_sigma*0.5;
			setVariable(proposedParams, "phase",          &phi_proposed);
			break;
		case 3:
			ra   = *(REAL8*) getVariable(currentParams, "rightascension");	/* radian      */
			ra_proposed   = ra   + gsl_ran_gaussian(GSLrandom,sigma[3])*big_sigma*0.05;
			setVariable(proposedParams, "rightascension", &ra_proposed);
			break;
		case 2:
			dec  = *(REAL8*) getVariable(currentParams, "declination");		/* radian      */
			dec_proposed  = dec  + gsl_ran_gaussian(GSLrandom,sigma[2])*big_sigma*0.05;
			setVariable(proposedParams, "declination",    &dec_proposed);
			break;
		case 1:
			psi  = *(REAL8*) getVariable(currentParams, "polarisation");		/* radian      */
			psi_proposed  = psi  + gsl_ran_gaussian(GSLrandom,sigma[1])*big_sigma*0.1;
			setVariable(proposedParams, "polarisation",   &psi_proposed);
			break;
		case 0:
			dist = *(REAL8*) getVariable(currentParams, "distance");			/* Mpc         */
			//dist_proposed = dist + gsl_ran_ugaussian(GSLrandom)*0.5;
			dist_proposed = dist * exp(gsl_ran_gaussian(GSLrandom,sigma[0])*big_sigma*0.1); // ~10% change
			logProposalRatio *= dist_proposed / dist;
			setVariable(proposedParams, "distance",       &dist_proposed);
			break;
	}
	
	// return ratio of proposal densities (for back & forth jumps) 
	// in "runState->proposalArgs" vector:
	if (checkVariable(runState->proposalArgs, "logProposalRatio"))
		setVariable(runState->proposalArgs, "logProposalRatio", &logProposalRatio);
	else
		addVariable(runState->proposalArgs, "logProposalRatio", &logProposalRatio, REAL8_t, PARAM_OUTPUT);
}




/*REAL8 GaussianLikelihood(LALVariables *currentParams, LALIFOData * data, LALTemplateFunction *template)
{
	
	double result=0.0;
	double sumsq=0.0;
	//double norm=0.0;
	//int i=0;
	double x[20];
	double xmax=0.0;
	double deltax=0.01;
	
	x[0]=*(REAL8 *)getVariable(currentParams,"x0");
	//for(i=0;i<run.nMCMCpar;i++){
	//	x[i]= par->par[run.parRevID[185+i]];
	//}
	
//	for(i=0;i<run.nMCMCpar;i++){
	//	sumsq+=(x[i]-xmax)*(x[i]-xmax)/(2*deltax);
		//norm+=-0.91893853320468-log(sqrt(deltax));
//	}
	sumsq=(x[0]-xmax)*(x[0]-xmax)/(2*deltax);
    //result=log(100*exp(-sumsq));
	//result=15/(2*deltax)-sumsq;
	result=1.0/(2.0*deltax)-sumsq;
	return result;

}*/

/*REAL8 UnityLikelihood(LALVariables *currentParams, LALIFOData * data, LALTemplateFunction *template)
{
	return 1.0;
}*/



/*REAL8 PTUniformGaussianPrior(LALInferenceRunState *runState, LALVariables *params)
{

	REAL8 x0;	
	REAL8 logdensity;
	
	x0   = *(REAL8*) getVariable(params, "x0");

	if(x0>=-1.0 && x0<=1.0)	
		logdensity = 0.0;
	else
		logdensity = -DBL_MAX;
	//TODO: should be properly normalized; pass in range via priorArgs?	
	
	return(logdensity);
	
}*/

/*void PTMCMCGaussianProposal(LALInferenceRunState *runState, LALVariables *proposedParams)
{
	
	REAL8 x0;
	REAL8 x0_proposed;
	REAL8 logProposalRatio = 0.0;  // = log(P(backward)/P(forward))
	gsl_rng * GSLrandom=runState->GSLrandom;
	LALVariables * currentParams = runState->currentParams;
	REAL8 sigma = 0.1;
	
	x0   = *(REAL8*) getVariable(currentParams, "x0");

	x0_proposed   = x0 + gsl_ran_ugaussian(GSLrandom)*sigma;
	//logProposalRatio *= x0_proposed / x0;   // (proposal ratio for above "scaled log-normal" proposal)

	
	copyVariables(currentParams, proposedParams);
	setVariable(proposedParams, "x0",      &(x0_proposed));		

	
	// return ratio of proposal densities (for back & forth jumps) 
	// in "runState->proposalArgs" vector:
	if (checkVariable(runState->proposalArgs, "logProposalRatio"))
		setVariable(runState->proposalArgs, "logProposalRatio", &logProposalRatio);
	else
		addVariable(runState->proposalArgs, "logProposalRatio", &logProposalRatio, REAL8_t, PARAM_OUTPUT);
	
	
	
}*/



void GetCartesianPos(REAL8 vec[3],REAL8 longitude, REAL8 latitude);
void GetCartesianPos(REAL8 vec[3],REAL8 longitude, REAL8 latitude)
{
	vec[0]=cos(longitude)*cos(latitude);
	vec[1]=sin(longitude)*cos(latitude);
	vec[1]=sin(latitude);
	return;
}

void CartesianToSkyPos(REAL8 pos[3],REAL8 *longitude, REAL8 *latitude);
void CartesianToSkyPos(REAL8 pos[3],REAL8 *longitude, REAL8 *latitude)
{
	REAL8 longi,lat,dist;
	dist=sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
	/*XLALMCMCSetParameter(parameter,"distMpc",dist);*/
	longi=atan2(pos[1]/dist,pos[0]/dist);
	if(longi<0.0) longi=LAL_TWOPI+longi;
	lat=asin(pos[2]/dist);
	*longitude=longi;
	*latitude=lat;	
	return;
}

void crossProduct(REAL8 out[3],REAL8 x[3],REAL8 y[3]);
void crossProduct(REAL8 out[3],REAL8 x[3],REAL8 y[3])
{
	out[0]=x[1]*y[2] - x[2]*y[1];
	out[1]=y[0]*x[2] - x[0]*y[2];
	out[2]=x[0]*y[1] - x[1]*y[0];
	return;
}

void normalise(REAL8 vec[3]);
void normalise(REAL8 vec[3]){
	REAL8 my_abs=0.0;
	my_abs=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
	vec[0]/=my_abs;
	vec[1]/=my_abs;
	vec[2]/=my_abs;
	return;
}


void PTMCMCLALInferenceRotateSky(
						   LALInferenceRunState *state,
						   LALVariables *parameter
						   )
{ /* Function to rotate the current sample around the vector between two random detectors */
	static LALStatus status;
	INT4 IFO1,IFO2;
	REAL4 randnum;
	REAL8 vec[3];
	REAL8 cur[3];
	REAL8 longi,lat;
	REAL8 vec_abs=0.0,theta,c,s;
	UINT4 i,j;
	
	copyVariables(state->currentParams, parameter);
	
	UINT4 nIFO=0;
	LALIFOData *ifodata1=state->data;
	while(ifodata1){
		nIFO++;
		ifodata1=ifodata1->next;
	}
	
	LALIFOData **IFOs=calloc(nIFO,sizeof(LALIFOData *));
	for(i=0,ifodata1=state->data;i<nIFO;i++){
		IFOs[i]=ifodata1;
		ifodata1=ifodata1->next;
	}
	
	
	if(nIFO<2) return;
	if(nIFO==2 && IFOs[0]==IFOs[1]) return;
	
	longi = *(REAL8 *)getVariable(parameter,"rightascension");
	lat = *(REAL8 *)getVariable(parameter,"declination");
	
	/* Convert the RA/dec to geodetic coordinates, as the detectors use these */
	SkyPosition geodetic,equatorial;
	equatorial.longitude=longi;
	equatorial.latitude=lat;
	equatorial.system=COORDINATESYSTEM_EQUATORIAL;
	geodetic.system=COORDINATESYSTEM_GEOGRAPHIC;
	LALEquatorialToGeographic(&status,&geodetic,&equatorial,&(IFOs[0]->epoch));
	longi=geodetic.longitude;
	lat=geodetic.latitude;
	cur[0]=cos(lat)*cos(longi);
	cur[1]=cos(lat)*sin(longi);
	cur[2]=sin(lat);
	
	IFO1 = gsl_rng_uniform_int(state->GSLrandom,nIFO);
	do{ /* Pick random interferometer other than the first one */
		IFO2 = gsl_rng_uniform_int(state->GSLrandom,nIFO);
	}while(IFO2==IFO1 || IFOs[IFO1]->detector==IFOs[IFO2]->detector);
	
	/*	fprintf(stderr,"Rotating around %s-%s vector\n",inputMCMC->ifoID[IFO1],inputMCMC->ifoID[IFO2]);*/
	/* Calc normalised direction vector */
	for(i=0;i<3;i++) vec[i]=IFOs[IFO2]->detector->location[i]-IFOs[IFO1]->detector->location[i];
	for(i=0;i<3;i++) vec_abs+=vec[i]*vec[i];
	vec_abs=sqrt(vec_abs);
	for(i=0;i<3;i++) vec[i]/=vec_abs;
	
	/* Chose random rotation angle */
	randnum=gsl_rng_uniform(state->GSLrandom);
	theta=LAL_TWOPI*randnum;
	c=cos(-theta); s=sin(-theta);
	/* Set up rotation matrix */
	double R[3][3] = {{c+vec[0]*vec[0]*(1.0-c), 
		vec[0]*vec[1]*(1.0-c)-vec[2]*s,
		vec[0]*vec[2]*(1.0-c)+vec[1]*s},
		{vec[1]*vec[0]*(1.0-c)+vec[2]*s,
			c+vec[1]*vec[1]*(1.0-c),
			vec[1]*vec[2]*(1.0-c)-vec[0]*s},
		{vec[2]*vec[0]*(1.0-c)-vec[1]*s,
			vec[2]*vec[1]*(1.0-c)+vec[0]*s,
			c+vec[2]*vec[2]*(1.0-c)}};
	REAL8 new[3]={0.0,0.0,0.0};
	for (i=0; i<3; ++i)
		for (j=0; j<3; ++j)
			new[i] += R[i][j]*cur[j];
	double newlong = atan2(new[1],new[0]);
	if(newlong<0.0) newlong=LAL_TWOPI+newlong;
	
	geodetic.longitude=newlong;
	geodetic.latitude=asin(new[2]);
	/* Convert back into equatorial (sky) coordinates */
	LALGeographicToEquatorial(&status,&equatorial,&geodetic,&(IFOs[0]->epoch));
	newlong=equatorial.longitude;
	double newlat=equatorial.latitude;
	
	/* Compute change in tgeocentre for this change in sky location */
	REAL8 dtold,dtnew,deltat;
	dtold = XLALTimeDelayFromEarthCenter(IFOs[0]->detector->location, longi, lat, &(IFOs[0]->epoch)); /* Compute time delay */
	dtnew = XLALTimeDelayFromEarthCenter(IFOs[0]->detector->location, newlong, newlat, &(IFOs[0]->epoch)); /* Compute time delay */
	deltat=dtold-dtnew; /* deltat is change in arrival time at geocentre */
	deltat+=*(REAL8 *)getVariable(parameter,"time");
	setVariable(parameter,"time",&deltat);	
	setVariable(parameter,"declination",&newlat);
	setVariable(parameter,"rightascension",&newlong);
	/*fprintf(stderr,"Skyrotate: new pos = %lf %lf %lf => %lf %lf\n",new[0],new[1],new[2],newlong,asin(new[2]));*/
	LALInferenceCyclicReflectiveBound(parameter,state->priorArgs);
	free(IFOs);
	return;
}


INT4 PTMCMCLALInferenceReflectDetPlane(
								 LALInferenceRunState *state,
								 LALVariables *parameter
								 )
{ /* Function to reflect a point on the sky about the plane of 3 detectors */
	/* Returns -1 if not possible */
	static LALStatus status;
	UINT4 i;
	int DetCollision=0;
	REAL4 randnum;
	REAL8 longi,lat;
	REAL8 dist;
	REAL8 pos[3];
	REAL8 normal[3];
	REAL8 w1[3]; /* work vectors */
	REAL8 w2[3];
	INT4 IFO1,IFO2,IFO3;
	REAL8 detvec[3];
	
	copyVariables(state->currentParams, parameter);
	
	UINT4 nIFO=0;
	LALIFOData *ifodata1=state->data;
	LALIFOData *ifodata2=NULL;
	while(ifodata1){
		nIFO++;
		ifodata1=ifodata1->next;
	}
	
	LALIFOData **IFOs=calloc(nIFO,sizeof(LALIFOData *));
	if(!IFOs) {
		printf("Unable to allocate memory for %i LALIFOData *s\n",nIFO);
		exit(1);
	}
	for(i=0,ifodata1=state->data;i<nIFO;i++){
		IFOs[i]=ifodata1;
		ifodata1=ifodata1->next;
	}
	
	if(nIFO<3) return(-1) ; /* not enough IFOs to construct a plane */
	for(ifodata1=state->data;ifodata1;ifodata1=ifodata1->next)
		for(ifodata2=ifodata1->next;ifodata2;ifodata2=ifodata2->next)
			if(ifodata1->detector==ifodata2->detector) DetCollision+=1;
	
	if(nIFO-DetCollision<3) return(-1); /* Not enough independent IFOs */
	
	/* Select IFOs to use */
	IFO1=gsl_rng_uniform_int(state->GSLrandom,nIFO);
	do {
		IFO2=gsl_rng_uniform_int(state->GSLrandom,nIFO);
	}while(IFO1==IFO2 || IFOs[IFO1]==IFOs[IFO2]);
	randnum=gsl_rng_uniform(state->GSLrandom);
	do {
		IFO3 = gsl_rng_uniform_int(state->GSLrandom,nIFO);
	}while(IFO3==IFO1
		   || IFO3==IFO2
		   || IFOs[IFO3]==IFOs[IFO1]
		   || IFOs[IFO3]==IFOs[IFO2]);
	/*fprintf(stderr,"Using %s, %s and %s for plane\n",inputMCMC->ifoID[IFO1],inputMCMC->ifoID[IFO2],inputMCMC->ifoID[IFO3]);*/
	
	longi = *(REAL8 *)getVariable(parameter,"rightascension");
	lat = *(REAL8 *)getVariable(parameter,"declination");
	
	double deltalong=0;
	
	/* Convert to earth coordinates */
	SkyPosition geodetic,equatorial;
	equatorial.longitude=longi;
	equatorial.latitude=lat;
	equatorial.system=COORDINATESYSTEM_EQUATORIAL;
	geodetic.system=COORDINATESYSTEM_GEOGRAPHIC;
	LAL_CALL(LALEquatorialToGeographic(&status,&geodetic,&equatorial,&(state->data->epoch)),&status);
	//LALEquatorialToGeographic(&status,&geodetic,&equatorial,&(state->data->epoch));
        deltalong=geodetic.longitude-equatorial.longitude;
	
	/* Add offset to RA to convert to earth-fixed */
	
	/* Calculate cartesian version of earth-fixed sky position */
	GetCartesianPos(pos,geodetic.longitude,lat); /* Get sky position in cartesian coords */
	
	
	/* calculate the unit normal vector of the detector plane */
	for(i=0;i<3;i++){ /* Two vectors in the plane */
		w1[i]=IFOs[IFO2]->detector->location[i] - IFOs[IFO1]->detector->location[i];
		w2[i]=IFOs[IFO3]->detector->location[i] - IFOs[IFO1]->detector->location[i];
		detvec[i]=IFOs[IFO1]->detector->location[i];
	}
	crossProduct(normal,w1,w2);
	normalise(normal);
	normalise(detvec);
	
	/* Calculate the distance between the point and the plane n.(point-IFO1) */
	for(dist=0.0,i=0;i<3;i++) dist+=pow(normal[i]*(pos[i]-detvec[i]),2.0);
	dist=sqrt(dist);
	/* Reflect the point pos across the plane */
	for(i=0;i<3;i++) pos[i]=pos[i]-2.0*dist*normal[i];
	
	REAL8 newLongGeo,newLat;
	CartesianToSkyPos(pos,&newLongGeo,&newLat);
	REAL8 newLongSky=newLongGeo-deltalong;
	
	
	setVariable(parameter,"rightascension",&newLongSky);
	setVariable(parameter,"declination",&newLat);
	
	/* Compute change in tgeocentre for this change in sky location */
	REAL8 dtold,dtnew,deltat;
	dtold = XLALTimeDelayFromEarthCenter(IFOs[0]->detector->location, longi, lat, &(IFOs[0]->epoch)); /* Compute time delay */
	dtnew = XLALTimeDelayFromEarthCenter(IFOs[0]->detector->location, newLongSky, newLat, &(IFOs[0]->epoch)); /* Compute time delay */
	deltat=dtold-dtnew; /* deltat is change in arrival time at geocentre */
	deltat+=*(REAL8 *)getVariable(parameter,"time");
	setVariable(parameter,"time",&deltat);
	
	LALInferenceCyclicReflectiveBound(parameter,state->priorArgs);
	free(IFOs);
	
	return(0);
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

static REAL8 dot(REAL8 v[3], REAL8 w[3]) {
  return v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
}

static REAL8 norm3(REAL8 v[3]) { return sqrt(dot(v,v)); }

static void cross(REAL8 x[3], REAL8 y[3], REAL8 z[3]) {
  z[0] = x[1]*y[2] - x[2]*y[1];
  z[1] = x[2]*y[0] - x[0]*y[2];
  z[2] = x[0]*y[1] - x[1]*y[0];
}

static void rotateVectorAboutVector(REAL8 v[3], REAL8 axis[3], REAL8 theta) {
  REAL8 an = norm3(axis);
  REAL8 x[3], y[3], z[3];
  REAL8 zdotv, xnorm;
  INT4 i;

  for (i = 0; i < 3; i++) {
    z[i] = axis[i]/an; /* zhat = axisHat */
  }

  zdotv = dot(z, v);

  for (i = 0; i < 3; i++) {
    x[i] = v[i] - zdotv*z[i]; /* Remove the z component from v, store in x. */
  }

  xnorm = norm3(x);

  if (xnorm == 0.0) return; /* v is along axis, rotation is done. */

  for (i = 0; i < 3; i++) {
    x[i] /= xnorm;
  }

  cross(z, x, y);  /* y = z \times x*/

  for (i = 0; i < 3; i++) {
    v[i] = zdotv*z[i] + xnorm*(cos(theta)*x[i]+sin(theta)*y[i]);
  }
}

static void thetaPhiToVector(REAL8 norm, REAL8 theta, REAL8 phi, REAL8 v[3]) {
  v[2] = norm*cos(theta);
  v[0] = norm*sin(theta)*cos(phi);
  v[1] = norm*sin(theta)*sin(phi);
}

static void vectorToThetaPhi(REAL8 *nrm, REAL8 *theta, REAL8 *phi, REAL8 v[3]) {
  *nrm = norm3(v);

  *theta = acos(v[2]/(*nrm));

  *phi = atan2(v[1], v[0]);

  if (*phi < 0.0) {
    *phi += 2.0*M_PI; /* We use 0 <= phi < 2*M_PI, while atan2 uses -M_PI < phi <= M_PI. */
  }
}

void PTMCMCLALInferenceRotateSpins(LALInferenceRunState *runState, LALVariables *proposedParams) {
  REAL8 inc, theta1, theta2, phi1, phi2;
  REAL8 L[3], A1[3], A2[3];
  REAL8 rotAngle;
  REAL8 selector;
  REAL8 dummy;

  copyVariables(runState->currentParams, proposedParams);

  inc = *((REAL8 *)getVariable(proposedParams, "inclination"));
  theta1 = *((REAL8 *)getVariable(proposedParams, "theta_spin1"));
  phi1 = *((REAL8 *)getVariable(proposedParams, "phi_spin1"));
  theta2 = *((REAL8 *)getVariable(proposedParams, "theta_spin2"));
  phi2 = *((REAL8 *)getVariable(proposedParams, "phi_spin2"));

  thetaPhiToVector(1.0, inc, 0.0, L); 
  thetaPhiToVector(1.0, theta1, phi1, A1);
  thetaPhiToVector(1.0, theta2, phi2, A2);

  rotAngle = 2.0*M_PI*gsl_rng_uniform(runState->GSLrandom);

  selector = gsl_rng_uniform(runState->GSLrandom);
  
  if (selector < 1.0/3.0) {
    /* Rotate both spins about L */
    rotateVectorAboutVector(A1, L, rotAngle);
    rotateVectorAboutVector(A2, L, rotAngle);
  } else if (selector < 2.0 / 3.0) {
    /* Rotate only A1 */
    rotateVectorAboutVector(A1, L, rotAngle);
  } else {
    /* Rotate only A2 */
    rotateVectorAboutVector(A2, L, rotAngle);
  }

  vectorToThetaPhi(&dummy, &theta1, &phi1, A1);
  vectorToThetaPhi(&dummy, &theta2, &phi2, A2);

  setVariable(proposedParams, "theta_spin1", &theta1);
  setVariable(proposedParams, "phi_spin1", &phi1);
  setVariable(proposedParams, "theta_spin2", &theta2);
  setVariable(proposedParams, "phi_spin2", &phi2);
}

void PTMCMCLALInferenceInclinationFlip(LALInferenceRunState *runState, LALVariables *proposedParams) {
  REAL8 iota;

  copyVariables(runState->currentParams, proposedParams);

  iota = *((REAL8 *) getVariable(proposedParams, "inclination"));
  
  if (runState->template==&templateLALSTPN) {
    /* Handle spins. */
    REAL8 dummyNorm, newIota, newPhi;
    REAL8 theta1, theta2, phi1, phi2;
    REAL8 L[3], a1[3], a2[3], xhat[3] = {1,0,0};

    theta1 = *((REAL8 *) getVariable(proposedParams, "theta_spin1"));
    theta2 = *((REAL8 *) getVariable(proposedParams, "theta_spin2"));

    phi1 = *((REAL8 *) getVariable(proposedParams, "phi_spin1"));
    phi2 = *((REAL8 *) getVariable(proposedParams, "phi_spin2"));

    thetaPhiToVector(1.0, iota, 0.0, L);
    thetaPhiToVector(1.0, theta1, phi1, a1);
    thetaPhiToVector(1.0, theta2, phi2, a2);

    rotateVectorAboutVector(L, xhat, M_PI-2.0*iota);
    rotateVectorAboutVector(a1, xhat, M_PI-2.0*iota);
    rotateVectorAboutVector(a2, xhat, M_PI-2.0*iota);

    vectorToThetaPhi(&dummyNorm, &newIota, &newPhi, L);
    vectorToThetaPhi(&dummyNorm, &theta1, &phi1, a1);
    vectorToThetaPhi(&dummyNorm, &theta2, &phi2, a2);

    if (fabs(newIota + iota - M_PI) > 1e-8 || fabs(newPhi) > 1e-8) {
      fprintf(stderr, "ERROR: inclination swap not implemented properly.\n");
      fprintf(stderr, "ERROR: should have new iota = Pi - iota, phi = 0 instead have\n");
      fprintf(stderr, "ERROR: new iota = %g, old iota = %g, phi = %g\n", newIota, iota, newPhi);
      exit(1);
    }

    setVariable(proposedParams, "phi_spin1", &phi1);
    setVariable(proposedParams, "phi_spin2", &phi2);
    setVariable(proposedParams, "theta_spin1", &theta1);
    setVariable(proposedParams, "theta_spin2", &theta2);
    /* Don't need to set iota because it will happen outside the if statement. */
  }

  iota = M_PI - iota;

  setVariable(proposedParams, "inclination", &iota);

  /* Not really needed (probably), but let's be safe. */
  LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs); 
}

void PTMCMCLALInferenceCovarianceEigenvectorJump(LALInferenceRunState *runState, LALVariables *proposedParams) {
  LALVariables *proposalArgs = runState->proposalArgs;
  gsl_matrix *eigenvectors = *((gsl_matrix **)getVariable(proposalArgs, "covarianceEigenvectors"));
  REAL8Vector *eigenvalues = *((REAL8Vector **)getVariable(proposalArgs, "covarianceEigenvalues"));
  REAL8 temp = *((REAL8 *)getVariable(proposalArgs, "temperature"));
  UINT4 N = eigenvalues->length;
  gsl_rng *rng = runState->GSLrandom;
  UINT4 i = gsl_rng_uniform_int(rng, N);
  REAL8 jumpSize = sqrt(temp*eigenvalues->data[i])*gsl_ran_ugaussian(rng);
  UINT4 j;
  LALVariableItem *proposeIterator;

  copyVariables(runState->currentParams, proposedParams);

  j = 0;
  proposeIterator = proposedParams->head;
  if (proposeIterator == NULL) {
    fprintf(stderr, "Bad proposed params in %s, line %d\n",
            __FILE__, __LINE__);
    exit(1);
  }
  do {
    if (proposeIterator->vary != PARAM_FIXED && proposeIterator->vary != PARAM_OUTPUT) {
      REAL8 tmp = *((REAL8 *)proposeIterator->value);
      REAL8 inc = jumpSize*gsl_matrix_get(eigenvectors, j, i);

      tmp += inc;

      memcpy(proposeIterator->value, &tmp, sizeof(REAL8));

      j++;
    }
  } while ((proposeIterator = proposeIterator->next) != NULL && j < N);

  LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);
}

void PTMCMCLALInferenceSkyLocWanderJump(LALInferenceRunState *runState, LALVariables *proposedParams) {
  gsl_rng *rng = runState->GSLrandom;
  LALVariables *proposalArgs = runState->proposalArgs;
  REAL8 temp = *((REAL8 *)getVariable(proposalArgs, "temperature"));
  REAL8 jumpX = sqrt(temp)*0.01*gsl_ran_ugaussian(rng)/sqrt(2.0);
  REAL8 jumpY = sqrt(temp)*0.01*gsl_ran_ugaussian(rng)/sqrt(2.0);
  REAL8 RA, DEC;

  copyVariables(runState->currentParams, proposedParams);

  RA = *((REAL8 *)getVariable(proposedParams, "rightascension"));
  DEC = *((REAL8 *)getVariable(proposedParams, "declination"));

  RA += jumpX/cos(DEC);
  DEC += jumpY;

  setVariable(proposedParams, "rightascension", &RA);
  setVariable(proposedParams, "declination", &DEC);

  LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);
}

/* Choose a jump in inclination and distance such that the signal
   amplitude in one of the detectors remains constant. */
void PTMCMCLALInferenceInclinationDistanceConstAmplitudeJump(LALInferenceRunState *runState, 
                                                             LALVariables *proposedParams) {
  UINT4 nIFO, iIFO;
  LALIFOData *ifoData;
  REAL8 fPlus, fCross;
  LIGOTimeGPS timeGPS;
  REAL8 ra, dec, psi, t, gmst;
  REAL8 iotaNew, cosIotaNew, iota, cosIota;
  REAL8 dNew, d;
  REAL8 norm;

  copyVariables(runState->currentParams, proposedParams);

  ra = *((REAL8 *)getVariable(proposedParams, "rightascension"));
  dec = *((REAL8 *)getVariable(proposedParams, "declination"));
  psi = *((REAL8 *)getVariable(proposedParams, "polarisation"));
  t = *((REAL8 *)getVariable(proposedParams, "time"));

  XLALGPSSetREAL8(&timeGPS, t);
  gmst = XLALGreenwichMeanSiderealTime(&timeGPS);

  /* Find number of IFO's. */
  nIFO=0;
  ifoData=runState->data;
  do {
    nIFO++;
    ifoData=ifoData->next;
  } while (ifoData != NULL);

  /* Now choose one at random, and get its data. */
  iIFO=gsl_rng_uniform_int(runState->GSLrandom, nIFO);
  ifoData=runState->data;
  while (iIFO > 0) {
    ifoData=ifoData->next;
    iIFO--;
  }

  /* Now we know which IFO we want to keep the magnitude constant in.
     Now compute f+, fx for that IFO. */
  XLALComputeDetAMResponse(&fPlus, &fCross, ifoData->detector->response, ra, dec, psi, gmst);

  iotaNew = M_PI*gsl_rng_uniform(runState->GSLrandom);
  cosIotaNew = cos(iotaNew);

  d = *((REAL8 *)getVariable(proposedParams, "distance"));

  iota = *((REAL8 *)getVariable(proposedParams, "inclination"));
  cosIota = cos(iota);

  norm = fabs((fPlus*(0.5*(1.0 + cosIota*cosIota)) + fCross*cosIota)/d);

  dNew = fabs((fPlus*(0.5*(1.0 + cosIotaNew*cosIotaNew)) + fCross*cosIotaNew) / norm);

  setVariable(proposedParams, "distance", &dNew);
  setVariable(proposedParams, "inclination", &iotaNew);

  LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);
  
}

void PTMCMCLALInferenceDifferentialEvolutionFull(LALInferenceRunState *runState, LALVariables *proposedParams) {
  PTMCMCLALInferenceDifferentialEvolutionNames(runState, proposedParams, NULL);
}

void PTMCMCLALInferenceDifferentialEvolutionNames(LALInferenceRunState *runState, 
                                                  LALVariables *proposedParams,
                                                  const char **names) {
  if (names == NULL) {
    size_t i;
    size_t N = getVariableDimension(runState->currentParams) + 1;
    names = alloca(N*sizeof(char *)); /* Hope we have alloca---saves
                                         having to deallocate after
                                         proposal. */
    for (i = 1; i <= N-1; i++) {
      names[i-1] = getVariableName(runState->currentParams, i);
    }

    names[N-1]=NULL; /* Terminate */
  }


  LALVariables **dePts = runState->differentialPoints;
  size_t nPts = runState->differentialPointsLength;

  if (dePts == NULL || nPts <= 1) {
    fprintf(stderr, "Trying to differentially evolve without enough points (in %s, %d)\n",
            __FILE__, __LINE__);
    exit(1);
  }

  copyVariables(runState->currentParams, proposedParams);

  size_t i,j;

  i = gsl_rng_uniform_int(runState->GSLrandom, nPts);
  do {
    j = gsl_rng_uniform_int(runState->GSLrandom, nPts);
  } while (j == i);

  LALVariables *ptI = dePts[i];
  LALVariables *ptJ = dePts[j];

  for (i = 0; names[i] != NULL; i++) {
    if (!checkVariable(proposedParams, names[i]) || !checkVariable(ptJ, names[i]) || !checkVariable(ptI, names[i])) {
      /* Ignore variable if it's not in each of the params. */
    } else {
      REAL8 x = *((REAL8 *)getVariable(proposedParams, names[i]));
      REAL8 scale = 1.66511*gsl_ran_ugaussian(runState->GSLrandom); 
      /* 1.66511 = number of sigma where Gaussian PDF drops to
         0.25. */
      
      x += scale * (*((REAL8 *) getVariable(ptJ, names[i])));
      x -= scale * (*((REAL8 *) getVariable(ptI, names[i])));
      
      setVariable(proposedParams, names[i], &x);
    }
  }

  LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);
}

void PTMCMCLALInferenceDifferentialEvolutionMasses(LALInferenceRunState *runState, LALVariables *pp) {
  const char *names[] = {"chirpmass", "massratio", NULL};
  PTMCMCLALInferenceDifferentialEvolutionNames(runState, pp, names);
}

void PTMCMCLALInferenceDifferentialEvolutionAmp(LALInferenceRunState *runState, LALVariables *pp) {
  const char *names[] = {"rightascension", "declination", "polarisation", "inclination", "distance", NULL};
  PTMCMCLALInferenceDifferentialEvolutionNames(runState, pp, names);
}

void PTMCMCLALInferenceDifferentialEvolutionSpins(LALInferenceRunState *runState, LALVariables *pp) {
  const char *names[] = {"a_spin1", "a_spin2", "phi_spin1", "phi_spin2", "theta_spin1", "theta_spin2", NULL};
  PTMCMCLALInferenceDifferentialEvolutionNames(runState, pp, names);
}

void PTMCMCLALInferenceDifferentialEvolutionSky(LALInferenceRunState *runState, LALVariables *pp) {
  const char *names[] = {"rightascension", "declination", NULL};
  PTMCMCLALInferenceDifferentialEvolutionNames(runState, pp, names);
}
