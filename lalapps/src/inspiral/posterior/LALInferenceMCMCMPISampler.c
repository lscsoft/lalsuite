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

//Test LALAlgorithm
void PTMCMCAlgorithm(struct tagLALInferenceRunState *runState)
{
	int i,t,tempi,tempj,p;
	int tempSwapCount=0;
	REAL8 tempDelta;
	int nChain;
	int count = 0;		//temporary counters to monitor the number of swaps between chains
	int MPIrank, MPIsize;
	LALStatus status;
	memset(&status,0,sizeof(status));
	//REAL8 dummyR8 = 1.0;
	REAL8 temperature = 1.0;
	REAL8 nullLikelihood;
	REAL8 logChainSwap = 0.0;
	int tempIndex;
	int *tempIndexVec = NULL;
	int dummyTemp;
	REAL8 *tempLadder = NULL;			//the temperature ladder
	double *TcurrentLikelihood = NULL; //the current likelihood for each chain
	//LALVariables* TcurrentParams = malloc(sizeof(LALVariables));	//the current parameters for each chains
	//LALVariables dummyLALVariable;
	
	INT4 nPar = getVariableDimensionNonFixed(runState->currentParams);
	INT4 Niter = *(INT4*) getVariable(runState->algorithmParams, "Niter");
	REAL8 tempMax = *(REAL8*) getVariable(runState->algorithmParams, "tempMax");   //max temperature in the temperature ladder

	MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);
	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

	nChain = MPIsize;		//number of parallel chain
	tempIndex = MPIrank;		//set initial temp indices


	tempLadder = malloc(nChain * sizeof(REAL8));			//the temperature ladder
	
	REAL8 **sigma = (REAL8 **)calloc(nChain,sizeof(REAL8 *));//matrix of sigmas per parameter and temperature for adaptation
	
	for (t=0; t<nChain; ++t) {
		sigma[t] = (REAL8 *)calloc(nPar,sizeof(REAL8));
		for (p=0; p<nPar; ++p) {
			sigma[t][p]=t*p;
		}
	}							 
	//REAL8 sigma = 0.1;
	
	
	if (nChain==1) tempLadder[0]=1.0;
	else {
		tempDelta = log(tempMax)/(REAL8)(nChain-1);
		for (t=0; t<nChain; ++t) tempLadder[t]=exp(t*tempDelta);
		}
	
	if (MPIrank == 0) {
		tempIndexVec = (int*) malloc(sizeof(int)*MPIsize);	//itialize temp index
		TcurrentLikelihood = (double*) malloc(sizeof(double)*nChain);

		for (t=0; t<nChain; ++t) {
			tempIndexVec[t] = t;
			printf("tempLadder[%d]=%f\n",t,tempLadder[t]);
		}
	}
	
	
	
	
	
	
	FILE **chainoutput = (FILE**)calloc(nChain,sizeof(FILE*));
	//char outfileName[99];

	char **outfileName = (char**)calloc(nChain,sizeof(char*));
	
	for (t=0; t<nChain; ++t) {
		outfileName[t] = (char*)calloc(99,sizeof(char*));
		sprintf(outfileName[t],"PTMCMC.output.%2.2d",t);
		chainoutput[t] = fopen(outfileName[t],"w");
		
		fprintf(chainoutput[t], "  SPINspiral version:%8.2f\n\n",1.0);
		fprintf(chainoutput[t], "%10s  %10s  %6s  %20s  %6s %8s   %6s  %8s  %10s  %12s  %9s  %9s  %8s\n",
				"nIter","Nburn","seed","null likelihood","Ndet","nCorr","nTemps","Tmax","Tchain","Network SNR","Waveform","pN order","Npar");
		fprintf(chainoutput[t], "%10d  %10d  %6d  %20.10lf  %6d %8d   %6d%10d%12.1f%14.6f  %9i  %9.1f  %8i\n",
				Niter,10,100000,0.0,1,1,nChain,(int)tempMax,tempLadder[t],50.0,4,2.0,nPar);
		fprintf(chainoutput[t], "\n%16s  %16s  %10s  %10s  %10s  %10s  %20s  %15s  %12s  %12s  %12s\n",
				"Detector","SNR","f_low","f_high","before tc","after tc","Sample start (GPS)","Sample length","Sample rate","Sample size","FT size");
		for(i=0;i<1;i++) {
			fprintf(chainoutput[t], "%16s  %16.8lf  %10.2lf  %10.2lf  %10.2lf  %10.2lf  %20.8lf  %15.7lf  %12d  %12d  %12d\n",
					"Hanford",50.0,40.0,350.0,6.00,1.00,
					864162757.00000,8.00,1024,9152,4577);
		}
		fprintf(chainoutput[t], "\n\n%31s","");
		fprintf(chainoutput[t], " %9i %9i %9i %9i %9i %9i %9i %9i %9i",61,62,11,22,31,32,51,41,52);
		fprintf(chainoutput[t],"\n");
		fprintf(chainoutput[t], "%8s %12s %9s","Cycle","log_Post.","Prior");
		fprintf(chainoutput[t], " %9s %9s %9s %9s %9s %9s %9s %9s %9s","Mc","eta","t_c","log(d)","R.A.","sin(dec)","cos(i)","phi_orb","psi");
		fprintf(chainoutput[t],"\n");
		
		//fprintf(chainoutput[t],"This is temperature chain %d of %d.\n", t, nChain);
		fclose(chainoutput[t]);
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
	
	addVariable(runState->proposalArgs, "temperature", &temperature,  REAL8_t, PARAM_LINEAR);	
	addVariable(runState->proposalArgs, "sigma", sigma,  REAL8_t, PARAM_FIXED);
	nullLikelihood = NullLogLikelihood(runState->data);
	//nullLikelihood = 0.0;
	// initialize starting likelihood value:
	runState->currentLikelihood = runState->likelihood(runState->currentParams, runState->data, runState->template);
	
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
	for (i=0; i<Niter; i++) {
		//printf(" MCMC iteration: %d\t", i+1);
		//copyVariables(&(TcurrentParams),runState->currentParams);
		setVariable(runState->proposalArgs, "temperature", &(tempLadder[tempIndex]));  //update temperature of the chain
		//dummyR8 = runState->currentLikelihood;
		runState->evolve(runState); //evolve the chain with the parameters TcurrentParams[t] at temperature tempLadder[t]
		//	if (runState->currentLikelihood != dummyR8) {
		//		printf(" accepted! new parameter values:\n");
		//		printVariables(runState->currentParams);
		//	}
		//copyVariables(runState->currentParams,&(TcurrentParams));
		//TcurrentLikelihood[t] = runState->currentLikelihood; // save the parameters and temperature.
		chainoutput[tempIndex] = fopen(outfileName[tempIndex],"a");
		fprintf(chainoutput[tempIndex], "%8d %12.5lf %9.6lf", i,runState->currentLikelihood - nullLikelihood,1.0);
		
		fprintf(chainoutput[tempIndex]," %9.5f",*(REAL8 *)getVariable(runState->currentParams,"chirpmass"));
		fprintf(chainoutput[tempIndex]," %9.5f",*(REAL8 *)getVariable(runState->currentParams,"massratio"));
		fprintf(chainoutput[tempIndex]," %9.5f",*(REAL8 *)getVariable(runState->currentParams,"time"));
		fprintf(chainoutput[tempIndex]," %9.5f",*(REAL8 *)getVariable(runState->currentParams,"distance"));
		fprintf(chainoutput[tempIndex]," %9.5f",*(REAL8 *)getVariable(runState->currentParams,"rightascension"));
		fprintf(chainoutput[tempIndex]," %9.5f",*(REAL8 *)getVariable(runState->currentParams,"declination"));
		fprintf(chainoutput[tempIndex]," %9.5f",*(REAL8 *)getVariable(runState->currentParams,"inclination"));
		fprintf(chainoutput[tempIndex]," %9.5f",*(REAL8 *)getVariable(runState->currentParams,"phase"));
		fprintf(chainoutput[tempIndex]," %9.5f",*(REAL8 *)getVariable(runState->currentParams,"polarisation"));
		//fprintf(chainoutput[tempIndex]," %9.5f",*(REAL8 *)getVariable(runState->currentParams,"x0"));
		
		fprintf(chainoutput[tempIndex],"\n");
		//fflush(chainoutput[tempIndex]);
		fclose(chainoutput[tempIndex]);

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
		MPI_Barrier(MPI_COMM_WORLD);	
			
		//printVariables(&(TcurrentParams[0]));
		if (MPIrank == 0) {
			for(tempi=0;tempi<nChain-1;tempi++) { //swap parameters and likelihood between chains
				for(tempj=tempi+1;tempj<nChain;tempj++) {
					
					logChainSwap = (1.0/tempLadder[tempi]-1.0/tempLadder[tempj]) * (TcurrentLikelihood[tempj]-TcurrentLikelihood[tempi]);
					
					if ((logChainSwap > 0)
						|| (log(gsl_rng_uniform(runState->GSLrandom)) < logChainSwap )) { //Then swap... 

						/*
						copyVariables(&(TcurrentParams[tempj]),&(dummyLALVariable));
						copyVariables(&(TcurrentParams[tempi]),&(TcurrentParams[tempj]));
						copyVariables(&(dummyLALVariable),&(TcurrentParams[tempi]));
						*/
						
						dummyTemp = tempIndexVec[tempj];
						tempIndexVec[tempj] = tempIndexVec[tempi];
						tempIndexVec[tempi] = dummyTemp;
						++tempSwapCount;
						/*
						dummyR8 = TcurrentLikelihood[tempj];
						TcurrentLikelihood[tempj] = TcurrentLikelihood[tempi];
						TcurrentLikelihood[tempi] = dummyR8;
						count++;
						*/
					}
				} //tempj
			} //tempi
		} //MPIrank==0
		MPI_Scatter(tempIndexVec, 1, MPI_INT, &tempIndex, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		//printf("%d\n",count);
		count = 0;
	}// for(i=0; i<100; i++)
	MPI_Barrier(MPI_COMM_WORLD);
	if (MPIrank == 0) printf("Temp swaps %d times\n", tempSwapCount);

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
	logPriorCurrent      = runState->prior(runState, runState->currentParams);
	logLikelihoodCurrent = runState->currentLikelihood;
	temperature = *(REAL8*) getVariable(runState->proposalArgs, "temperature");
	
	// generate proposal:
	proposedParams.head = NULL;
	proposedParams.dimension = 0;
	runState->proposal(runState, &proposedParams);
	if (checkVariable(runState->proposalArgs, "logProposalRatio"))
		logProposalRatio = *(REAL8*) getVariable(runState->proposalArgs, "logProposalRatio");
	
	// compute prior & likelihood:
	logPriorProposed = runState->prior(runState, &proposedParams);
	if (logPriorProposed > -HUGE_VAL)
		logLikelihoodProposed = runState->likelihood(&proposedParams, runState->data, runState->template);
	else
		logLikelihoodProposed = -HUGE_VAL;
	
	// determine acceptance probability:
	logAcceptanceProbability = (1.0/temperature)*(logLikelihoodProposed - logLikelihoodCurrent) 
	+ (logPriorProposed - logPriorCurrent)
	+ logProposalRatio;
	
//	fprintf(stdout," %9.5f=(1.0 / %9.5f)*(%9.5f-%9.5f)+(%9.5f-%9.5f)+%9.5f\n",logAcceptanceProbability,temperature,logLikelihoodProposed,logLikelihoodCurrent,logPriorProposed,logPriorCurrent,logProposalRatio);
	
	// accept/reject:
	if ((logAcceptanceProbability > 0) 
		|| (log(gsl_rng_uniform(runState->GSLrandom)) < logAcceptanceProbability)) {   //accept
		copyVariables(&proposedParams, runState->currentParams);
		runState->currentLikelihood = logLikelihoodProposed;
	}
	
	destroyVariables(&proposedParams);	
}

//Test LALPriorFunction
REAL8 PTUniformLALPrior(LALInferenceRunState *runState, LALVariables *params)
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
{
	REAL8 mc, eta, iota, phi, tc, ra, dec, psi, dist;	
	REAL8 logdensity;
	
	mc   = *(REAL8*) getVariable(params, "chirpmass");		/* solar masses*/
	eta  = *(REAL8*) getVariable(params, "massratio");		/* dim-less    */
	iota = *(REAL8*) getVariable(params, "inclination");		/* radian      */
	tc   = *(REAL8*) getVariable(params, "time");			/* GPS seconds */
	phi  = *(REAL8*) getVariable(params, "phase");		/* radian      */
	ra   = *(REAL8*) getVariable(params, "rightascension");	/* radian      */
	dec  = *(REAL8*) getVariable(params, "declination");		/* radian      */
	psi  = *(REAL8*) getVariable(params, "polarisation"); 	/* radian      */
	dist = *(REAL8*) getVariable(params, "distance");		/* Mpc         */
	
	if(mc>0.0 && mc<=30.0 && eta>0.03 && eta<=0.25 && iota>=0.0 && iota<=LAL_PI && phi>=0.0 && phi<=LAL_TWOPI 
	   && ra>=0.0 && ra<=LAL_TWOPI && dec>=-LAL_PI_2 && dec<=LAL_PI_2 && psi>=0.0 && psi<=LAL_TWOPI && dist>0.0 && dist<=100.0)	
		logdensity = 0.0;
	else
		logdensity = -HUGE_VAL;
	//TODO: should be properly normalized; pass in range via priorArgs?	
	
	return(logdensity);
}


//Test LALProposalFunction
void PTMCMCLALProposal(LALInferenceRunState *runState, LALVariables *proposedParams)
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
	REAL8 logProposalRatio = 0.0;  // = log(P(backward)/P(forward))
	gsl_rng * GSLrandom=runState->GSLrandom;
	LALVariables * currentParams = runState->currentParams;
	REAL8 sigma = 0.1;
	
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
	mc_proposed   = mc * exp(gsl_ran_ugaussian(GSLrandom)*sigma*0.001);          /* mc changed by ~0.1% */
	logProposalRatio *= mc_proposed / mc;   // (proposal ratio for above "scaled log-normal" proposal)
	eta_proposed  = eta  + gsl_ran_ugaussian(GSLrandom)*sigma*0.01; /*eta changed by 0.01*/
	//TODO: if(eta_proposed>0.25) eta_proposed=0.25-(eta_proposed-0.25); etc.
	iota_proposed = iota + gsl_ran_ugaussian(GSLrandom)*sigma*0.1;
	tc_proposed   = tc   + gsl_ran_ugaussian(GSLrandom)*sigma*0.005; /*time changed by 5 ms*/
	phi_proposed  = phi  + gsl_ran_ugaussian(GSLrandom)*sigma*0.5;
	ra_proposed   = ra   + gsl_ran_ugaussian(GSLrandom)*sigma*0.05;
	dec_proposed  = dec  + gsl_ran_ugaussian(GSLrandom)*sigma*0.05;
	psi_proposed  = psi  + gsl_ran_ugaussian(GSLrandom)*sigma*0.1;
	//dist_proposed = dist + gsl_ran_ugaussian(GSLrandom)*0.5;
	dist_proposed = dist * exp(gsl_ran_ugaussian(GSLrandom)*sigma*0.1); // ~10% change
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
	REAL8 logProposalRatio = 0.0;  // = log(P(backward)/P(forward))
	gsl_rng * GSLrandom=runState->GSLrandom;
	LALVariables * currentParams = runState->currentParams;
	REAL8 sigmat = 0.1;
	
	REAL8 **sigma = NULL;
	
	sigma = *(REAL8***) getVariable(runState->proposalArgs, "sigma");
	
	
	printf("%f\n",sigma[0][0]);
	int t,p;
	for (t=0; t<5; ++t){
		for (p=0; p<9; ++p){
		//	printf("sigma[%d][%d]=%f\n",t,p,sigma[t][p]);
		}
	}
	
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
	mc_proposed   = mc * exp(gsl_ran_gaussian(GSLrandom,sigmat)*0.001);          /* mc changed by ~0.1% */
	logProposalRatio *= mc_proposed / mc;   // (proposal ratio for above "scaled log-normal" proposal)
	eta_proposed  = eta  + gsl_ran_gaussian(GSLrandom,sigmat)*0.01; /*eta changed by 0.01*/
	//TODO: if(eta_proposed>0.25) eta_proposed=0.25-(eta_proposed-0.25); etc.
	iota_proposed = iota + gsl_ran_gaussian(GSLrandom,sigmat)*0.1;
	tc_proposed   = tc   + gsl_ran_gaussian(GSLrandom,sigmat)*0.005; /*time changed by 5 ms*/
	phi_proposed  = phi  + gsl_ran_gaussian(GSLrandom,sigmat)*0.5;
	ra_proposed   = ra   + gsl_ran_gaussian(GSLrandom,sigmat)*0.05;
	dec_proposed  = dec  + gsl_ran_gaussian(GSLrandom,sigmat)*0.05;
	psi_proposed  = psi  + gsl_ran_gaussian(GSLrandom,sigmat)*0.1;
	//dist_proposed = dist + gsl_ran_ugaussian(GSLrandom)*0.5;
	dist_proposed = dist * exp(gsl_ran_gaussian(GSLrandom,sigmat)*0.1); // ~10% change
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



REAL8 GaussianLikelihood(LALVariables *currentParams, LALIFOData * data, LALTemplateFunction *template)
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

}

REAL8 PTUniformGaussianPrior(LALInferenceRunState *runState, LALVariables *params)
{

	REAL8 x0;	
	REAL8 logdensity;
	
	x0   = *(REAL8*) getVariable(params, "x0");

	if(x0>=-1.0 && x0<=1.0)	
		logdensity = 0.0;
	else
		logdensity = -HUGE_VAL;
	//TODO: should be properly normalized; pass in range via priorArgs?	
	
	return(logdensity);
	
}

void PTMCMCGaussianProposal(LALInferenceRunState *runState, LALVariables *proposedParams)
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
	
	
	
}
