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
#include <lal/SeqFactories.h>
#include <lal/Date.h>
#include <lal/VectorOps.h>
#include <lal/TimeFreqFFT.h>
#include <lal/GenerateInspiral.h>
#include <mpi.h>
#include "LALInference.h"
#include "mpi.h"
#include "LALInferenceMCMCMPISampler.h"
#include "LALInferencePrior.h"

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
	REAL8 dummyR8 = 0.0;
	REAL8 temperature = 1.0;
	REAL8 nullLikelihood;
	REAL8 logChainSwap = 0.0;
	int tempIndex;
	int *tempIndexVec = NULL;
	int dummyTemp;
	REAL8 *tempLadder = NULL;			//the temperature ladder
	double *TcurrentLikelihood = NULL; //the current likelihood for each chain
	REAL8 *sigmaVec = NULL;
	//LALVariables* TcurrentParams = malloc(sizeof(LALVariables));	//the current parameters for each chains
	//LALVariables dummyLALVariable;
	
	INT4 nPar = getVariableDimensionNonFixed(runState->currentParams);
	INT4 Niter = *(INT4*) getVariable(runState->algorithmParams, "Niter");
	INT4 Nskip = *(INT4*) getVariable(runState->algorithmParams, "Nskip");
	REAL8 tempMax = *(REAL8*) getVariable(runState->algorithmParams, "tempMax");   //max temperature in the temperature ladder
	INT4 randomseed = *(INT4*) getVariable(runState->algorithmParams,"random_seed");

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
	
	
	if (nChain==1) tempLadder[0]=1.0;
	else {
		tempDelta = log(tempMax)/(REAL8)(nChain-1);
		for (t=0; t<nChain; ++t) tempLadder[t]=exp(t*tempDelta);
		}
	
	if (MPIrank == 0) {
		tempIndexVec = (int*) malloc(sizeof(int)*nChain);	//initialize temp index
		TcurrentLikelihood = (double*) malloc(sizeof(double)*nChain);

		for (t=0; t<nChain; ++t) {
			tempIndexVec[t] = t;
			printf("tempLadder[%d]=%f\n",t,tempLadder[t]);
		}
	}
	
	
	nullLikelihood = NullLogLikelihood(runState->data);
	//nullLikelihood = 0.0;
	// initialize starting likelihood value:
	runState->currentLikelihood = runState->likelihood(runState->currentParams, runState->data, runState->template);
	
	
	
	FILE **chainoutput = (FILE**)calloc(nChain,sizeof(FILE*));
	//char outfileName[99];

	char **outfileName = (char**)calloc(nChain,sizeof(char*));
	
	for (t=0; t<nChain; ++t) {
		outfileName[t] = (char*)calloc(99,sizeof(char*));
		sprintf(outfileName[t],"PTMCMC.output.%d.%2.2d",randomseed,t);
		if (MPIrank == 0) {
			chainoutput[t] = fopen(outfileName[t],"w");
			fprintf(chainoutput[t], "  SPINspiral version:%8.2f\n\n",1.0);
			fprintf(chainoutput[t], "%10s  %10s  %6s  %20s  %6s %8s   %6s  %8s  %10s  %12s  %9s  %9s  %8s\n",
					"nIter","Nburn","seed","null likelihood","Ndet","nCorr","nTemps","Tmax","Tchain","Network SNR","Waveform","pN order","Npar");
			fprintf(chainoutput[t], "%10d  %10d  %6d  %20.10lf  %6d %8d   %6d%10d%12.1f%14.6f  %9i  %9.1f  %8i\n",
					Niter,10,100000,nullLikelihood,1,1,nChain,(int)tempMax,tempLadder[t],50.0,4,2.0,nPar);
			fprintf(chainoutput[t], "\n%16s  %16s  %10s  %10s  %10s  %10s  %20s  %15s  %12s  %12s  %12s\n",
					"Detector","SNR","f_low","f_high","before tc","after tc","Sample start (GPS)","Sample length","Sample rate","Sample size","FT size");
			for(i=0;i<1;i++) {
					fprintf(chainoutput[t], "%16s  %16.8lf  %10.2lf  %10.2lf  %10.2lf  %10.2lf  %20.8lf  %15.7lf  %12d  %12d  %12d\n",
							"Hanford",50.0,40.0,350.0,6.00,1.00,
							864162757.00000,8.00,1024,9152,4577);
			}
			fprintf(chainoutput[t], "\n\n%31s","");
			fprintf(chainoutput[t], " %9i %9i %9i %9i %9i %9i %9i %9i %9i",55,52,33,31,23,41,11,62,61);
			//fprintf(chainoutput[t], " %9i",185);
			fprintf(chainoutput[t],"\n");
			fprintf(chainoutput[t], "%8s %12s %9s","Cycle","log_Post.","Prior");
			fprintf(chainoutput[t], " %9s %9s %9s %9s %9s %9s %9s %9s %9s","iota","psi","dec","R.A.","dist","phi_orb","t_c","eta","Mc");
			//fprintf(chainoutput[t], " %9s","x1");
			fprintf(chainoutput[t],"\n");
			fprintf(chainoutput[t], "%d\t%f\t%f\t", 0,runState->currentLikelihood - nullLikelihood,1.0);
			fprintSampleNonFixed(chainoutput[t],runState->currentParams);
			fprintf(chainoutput[t],"%f\t",tempLadder[t]);
			fprintf(chainoutput[t],"%d\t",MPIrank);
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
	}
	MPI_Barrier(MPI_COMM_WORLD);
	// iterate:
	for (i=1; i<=Niter; i++) {
		//printf(" MCMC iteration: %d\t", i+1);
		//copyVariables(&(TcurrentParams),runState->currentParams);
		setVariable(runState->proposalArgs, "temperature", &(tempLadder[tempIndex]));  //update temperature of the chain
		setVariable(runState->proposalArgs, "tempIndex", &(tempIndex));
		s_gamma=10.0*exp(-(1.0/6.0)*log((double)i));
		setVariable(runState->proposalArgs, "s_gamma", &(s_gamma));
		//setVariable(runState->proposalArgs, "sigma", sigmaVec[tempIndex]);
		//dummyR8 = runState->currentLikelihood;
		//	if (runState->currentLikelihood != dummyR8) {
		//		printf(" accepted! new parameter values:\n");
		//		printVariables(runState->currentParams);
		//	}
		//copyVariables(runState->currentParams,&(TcurrentParams));
		//TcurrentLikelihood[t] = runState->currentLikelihood; // save the parameters and temperature.
		runState->evolve(runState); //evolve the chain with the parameters TcurrentParams[t] at temperature tempLadder[t]
		if ((i % Nskip) == 0){
			//chainoutput[tempIndex] = fopen(outfileName[tempIndex],"a");
			//fprintf(chainoutput[tempIndex], "%8d %12.5lf %9.6lf", i,runState->currentLikelihood - nullLikelihood,1.0);
			fprintf(chainoutput[tempIndex], "%d\t%f\t%f\t", i,runState->currentLikelihood - nullLikelihood,1.0);
			/*fprintf(chainoutput[tempIndex]," %9.5f",*(REAL8 *)getVariable(runState->currentParams,"chirpmass"));
			 fprintf(chainoutput[tempIndex]," %9.5f",*(REAL8 *)getVariable(runState->currentParams,"massratio"));
			 fprintf(chainoutput[tempIndex]," %9.5f",*(REAL8 *)getVariable(runState->currentParams,"time"));
			 fprintf(chainoutput[tempIndex]," %9.5f",*(REAL8 *)getVariable(runState->currentParams,"distance"));
			 fprintf(chainoutput[tempIndex]," %9.5f",*(REAL8 *)getVariable(runState->currentParams,"rightascension"));
			 fprintf(chainoutput[tempIndex]," %9.5f",*(REAL8 *)getVariable(runState->currentParams,"declination"));
			 fprintf(chainoutput[tempIndex]," %9.5f",*(REAL8 *)getVariable(runState->currentParams,"inclination"));
			 fprintf(chainoutput[tempIndex]," %9.5f",*(REAL8 *)getVariable(runState->currentParams,"phase"));
			 fprintf(chainoutput[tempIndex]," %9.5f",*(REAL8 *)getVariable(runState->currentParams,"polarisation"));*/
			//fprintf(chainoutput[tempIndex]," %9.5f",*(REAL8 *)getVariable(runState->currentParams,"x0"));
			fprintSampleNonFixed(chainoutput[tempIndex],runState->currentParams);
			fprintf(chainoutput[tempIndex],"%f\t",tempLadder[tempIndex]);
			fprintf(chainoutput[tempIndex],"%d\t",MPIrank);
			fprintf(chainoutput[tempIndex],"\n");
			fflush(chainoutput[tempIndex]);
			//fclose(chainoutput[tempIndex]);
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
		MPI_Gather(sigma->data,nPar,MPI_DOUBLE,sigmaVec,nPar,MPI_DOUBLE,0,MPI_COMM_WORLD);
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
						/*
						dummyR8 = TcurrentLikelihood[tempj];
						TcurrentLikelihood[tempj] = TcurrentLikelihood[tempi];
						TcurrentLikelihood[tempi] = dummyR8;
						count++;
						*/
						for (p=0; p<(nPar); ++p){
							dummyR8=sigmaVec[p+nPar*upperRank];
							sigmaVec[p+nPar*upperRank]=sigmaVec[p+nPar*lowerRank];
							sigmaVec[p+nPar*lowerRank]=dummyR8;
						}
					}
				} //upperRank
			} //lowerRank
		} //MPIrank==0
		MPI_Scatter(tempIndexVec, 1, MPI_INT, &tempIndex, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Scatter(sigmaVec,nPar,MPI_DOUBLE,sigma->data,nPar,MPI_DOUBLE,0, MPI_COMM_WORLD);
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
	
	// current values:
	//logPriorCurrent      = runState->prior(runState, runState->currentParams);
	logPriorCurrent      = runState->currentPrior;
	logLikelihoodCurrent = runState->currentLikelihood;
	temperature = *(REAL8*) getVariable(runState->proposalArgs, "temperature");
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
		runState->currentPrior = logPriorProposed;
		//}
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



//Test LALProposalFunction
void PTMCMCLALProposaltemp(LALInferenceRunState *runState, LALVariables *proposedParams)
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
	REAL8 a_spin1, a_spin2, theta_spin1, theta_spin2, phi_spin1, phi_spin2;
	REAL8 mc_proposed, eta_proposed, iota_proposed, phi_proposed, tc_proposed, 
	ra_proposed, dec_proposed, psi_proposed, dist_proposed;
	REAL8 a_spin1_proposed, a_spin2_proposed, theta_spin1_proposed, 
	theta_spin2_proposed, phi_spin1_proposed, phi_spin2_proposed;
	REAL8 logProposalRatio = 0.0;  // = log(P(backward)/P(forward))
	gsl_rng * GSLrandom=runState->GSLrandom;
	LALVariables * currentParams = runState->currentParams;
	copyVariables(currentParams, proposedParams);
	
	REAL8 sigma = 0.1;
	REAL8 big_sigma = 1.0;
	
	if(gsl_ran_ugaussian(GSLrandom) < 1.0e-3) big_sigma = 1.0e1;    //Every 1e3 iterations, take a 10x larger jump in all parameters
	if(gsl_ran_ugaussian(GSLrandom) < 1.0e-4) big_sigma = 1.0e2;    //Every 1e4 iterations, take a 100x larger jump in all parameters
	
	
	mc   = *(REAL8*) getVariable(currentParams, "chirpmass");		/* solar masses*/
	eta  = *(REAL8*) getVariable(currentParams, "massratio");		/* dim-less    */
	iota = *(REAL8*) getVariable(currentParams, "inclination");		/* radian      */
	tc   = *(REAL8*) getVariable(currentParams, "time");				/* GPS seconds */
	phi  = *(REAL8*) getVariable(currentParams, "phase");			/* radian      */
	ra   = *(REAL8*) getVariable(currentParams, "rightascension");	/* radian      */
	dec  = *(REAL8*) getVariable(currentParams, "declination");		/* radian      */
	psi  = *(REAL8*) getVariable(currentParams, "polarisation");		/* radian      */
	dist = *(REAL8*) getVariable(currentParams, "distance");			/* Mpc         */
	
	if (checkVariable(currentParams, "a_spin1")){
		a_spin1 = *(REAL8*) getVariable(currentParams, "a_spin1");
		a_spin1_proposed = a_spin1 + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.001;
		setVariable(proposedParams, "a_spin1",      &a_spin1_proposed);
	}
	if (checkVariable(currentParams, "theta_spin1")){
		theta_spin1 = *(REAL8*) getVariable(currentParams, "theta_spin1");
		theta_spin1_proposed = theta_spin1 + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;
		setVariable(proposedParams, "theta_spin1",      &theta_spin1_proposed);
	}
	if (checkVariable(currentParams, "phi_spin1")){
		phi_spin1 = *(REAL8*) getVariable(currentParams, "phi_spin1");
		phi_spin1_proposed = phi_spin1 + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;
		setVariable(proposedParams, "phi_spin1",      &phi_spin1_proposed);
	}
	if (checkVariable(currentParams, "a_spin2")){
		a_spin2 = *(REAL8*) getVariable(currentParams, "a_spin2");
		a_spin2_proposed = a_spin2 + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.001;
		setVariable(proposedParams, "a_spin2",      &a_spin2_proposed);
	}
	if (checkVariable(currentParams, "theta_spin2")){
		theta_spin2 = *(REAL8*) getVariable(currentParams, "theta_spin2");
		theta_spin2_proposed = theta_spin2 + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;
		setVariable(proposedParams, "theta_spin2",      &theta_spin2_proposed);
	}
	if (checkVariable(currentParams, "phi_spin2")){
		phi_spin2 = *(REAL8*) getVariable(currentParams, "phi_spin2");
		phi_spin2_proposed = phi_spin2 + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;
		setVariable(proposedParams, "phi_spin2",      &phi_spin2_proposed);
	}

	//mc_proposed   = mc*(1.0+gsl_ran_ugaussian(GSLrandom)*0.01);	/*mc changed by 1% */
	// (above proposal is not symmetric!)
	mc_proposed   = mc   + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;	/*mc changed by 0.0001 */
	//mc_proposed   = mc * exp(gsl_ran_ugaussian(GSLrandom)*sigma*0.01);          /* mc changed by ~0.1% */
	
	eta_proposed  = eta  + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.001; /*eta changed by 0.01*/
	//TODO: if(eta_proposed>0.25) eta_proposed=0.25-(eta_proposed-0.25); etc.
	iota_proposed = iota + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.5;
	tc_proposed   = tc   + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.001; /*time changed by 5 ms*/
	phi_proposed  = phi  + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.1;
	ra_proposed   = ra   + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;
	dec_proposed  = dec  + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;
	psi_proposed  = psi  + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.1;
	//dist_proposed = dist + gsl_ran_ugaussian(GSLrandom)*0.5;
	dist_proposed = dist * exp(gsl_ran_ugaussian(GSLrandom)*sigma*0.1); // ~10% change
	
	
	
	setVariable(proposedParams, "chirpmass",      &mc_proposed);		
	setVariable(proposedParams, "massratio",      &eta_proposed);
	setVariable(proposedParams, "inclination",    &iota_proposed);
	setVariable(proposedParams, "phase",          &phi_proposed);
	setVariable(proposedParams, "time",           &tc_proposed); 
	setVariable(proposedParams, "rightascension", &ra_proposed);
	setVariable(proposedParams, "declination",    &dec_proposed);
	setVariable(proposedParams, "polarisation",   &psi_proposed);
	setVariable(proposedParams, "distance",       &dist_proposed);
	
	LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);
	
	dist_proposed = *(REAL8*) getVariable(proposedParams, "distance");
	logProposalRatio *= dist_proposed / dist;
	//logProposalRatio *= mc_proposed / mc;   // (proposal ratio for above "scaled log-normal" proposal)
	
	// return ratio of proposal densities (for back & forth jumps) 
	// in "runState->proposalArgs" vector:
	if (checkVariable(runState->proposalArgs, "logProposalRatio"))
		setVariable(runState->proposalArgs, "logProposalRatio", &logProposalRatio);
	else
		addVariable(runState->proposalArgs, "logProposalRatio", &logProposalRatio, REAL8_t, PARAM_OUTPUT);
}



void PTMCMCLALProposal(LALInferenceRunState *runState, LALVariables *proposedParams)
{
	gsl_rng * GSLrandom=runState->GSLrandom;
	LALVariableItem *paraHead=NULL;
	INT4 i;
	copyVariables(runState->currentParams, proposedParams);
	
	REAL8 sigma = 0.1;
	REAL8 big_sigma = 1.0;
	
	if(gsl_ran_ugaussian(GSLrandom) < 1.0e-3) big_sigma = 1.0e1;    //Every 1e3 iterations, take a 10x larger jump in all parameters
	if(gsl_ran_ugaussian(GSLrandom) < 1.0e-4) big_sigma = 1.0e2;    //Every 1e4 iterations, take a 100x larger jump in all parameters

	/* loop over all parameters */
	for (paraHead=proposedParams->head,i=0; paraHead; paraHead=paraHead->next)
	{ 
		if(paraHead->vary==PARAM_LINEAR || paraHead->vary==PARAM_CIRCULAR){
			*(REAL8 *)paraHead->value += gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;
			i++;
		}
	}
	
	LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);

}


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
