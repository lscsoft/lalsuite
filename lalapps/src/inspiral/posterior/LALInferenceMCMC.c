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
#include "LALInference.h"


//Test LALAlgorithm
void PTMCMCAlgorithm(struct tagLALInferenceRunState *runState)
{
	int i,t,tempi,tempj;
	int nChain = 5;
	int count = 0;
	LALStatus status;
	memset(&status,0,sizeof(status));
	REAL8 dummyR8,temperature = 1.0;
	REAL8 nullLikelihood;
	REAL8 logChainSwap = 0.0;
	REAL8 *tempLadder = malloc(nChain * sizeof(REAL8));
	REAL8 *TcurrentLikelihood = malloc(nChain * sizeof(REAL8));
	LALVariables* TcurrentParams = malloc(nChain * sizeof(LALVariables));
	LALVariables dummyLALVariable;

	for(t=0; t<nChain; t++) {
		TcurrentParams[t].head=NULL;
		TcurrentParams[t].dimension=0;
		copyVariables(runState->currentParams,&(TcurrentParams[t]));
	}
	dummyLALVariable.head=NULL;
	dummyLALVariable.dimension=0;
	copyVariables(runState->currentParams,&(dummyLALVariable));
	
	tempLadder[0]=1.0;
	tempLadder[1]=2.0;
	tempLadder[2]=5.0;
	tempLadder[3]=10.0;
	tempLadder[4]=20.0;
	
	addVariable(runState->proposalArgs, "temperature", &temperature,  REAL8_t);

	nullLikelihood = NullLogLikelihood(runState->data);
	
	printf(" PTMCMCAlgorithm(); starting parameter values:\n");
	printVariables(runState->currentParams);
	// initialize starting likelihood value:
	runState->currentLikelihood = runState->likelihood(runState->currentParams, runState->data, runState->template);
	for(t=0; t<nChain; t++) { TcurrentLikelihood[t] = runState->currentLikelihood; }
	
	// iterate:
	for(i=0; i<100; i++) {
		printf(" MCMC iteration: %d\t", i+1);
		for(t=0; t<nChain; t++) {
			copyVariables(&(TcurrentParams[t]),runState->currentParams);
			setVariable(runState->proposalArgs, "temperature", &(tempLadder[t]));
			//dummyR8 = runState->currentLikelihood;
			runState->evolve(runState);
			//	if (runState->currentLikelihood != dummyR8) {
			//		printf(" accepted! new parameter values:\n");
			//		printVariables(runState->currentParams);
			//	}
			copyVariables(runState->currentParams,&(TcurrentParams[t]));
			TcurrentLikelihood[t] = runState->currentLikelihood;
			printf("%f\t", runState->currentLikelihood - nullLikelihood);
		} //for(t=0; t<nChain; t++)
		
		for(tempi=0;tempi<nChain-1;tempi++) {
			for(tempj=tempi+1;tempj<nChain;tempj++) {
				
				logChainSwap = (1.0/tempLadder[tempi]-1.0/tempLadder[tempj]) * (TcurrentLikelihood[tempj]-TcurrentLikelihood[tempi]);
				
				if ((logChainSwap > 0)
					|| (log(gsl_rng_uniform(runState->GSLrandom)) < logChainSwap )) { //Then swap...

					copyVariables(&(TcurrentParams[tempj]),&(dummyLALVariable));
					copyVariables(&(TcurrentParams[tempi]),&(TcurrentParams[tempj]));
					copyVariables(&(dummyLALVariable),&(TcurrentParams[tempi]));
					
					dummyR8 = TcurrentLikelihood[tempj];
					TcurrentLikelihood[tempj] = TcurrentLikelihood[tempi];
					TcurrentLikelihood[tempi] = dummyR8;
					count++;
				}
			} //tempj
		} //tempi
		
		printf("%d\n",count);
		count = 0;
	}// for(i=0; i<100; i++)	
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
	
	if(eta>0.0 && eta<=0.25 && iota>=0.0 && iota<=LAL_PI && phi>=0.0 && phi<=LAL_TWOPI 
	   && ra>=0.0 && ra<=LAL_TWOPI && dec>=-LAL_PI_2 && dec<=LAL_PI_2 && psi>=0.0 && psi<=LAL_TWOPI)	
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
	mc_proposed   = mc * exp(gsl_ran_ugaussian(GSLrandom)*0.001);          /* mc changed by ~0.1% */
	logProposalRatio *= mc_proposed / mc;   // (proposal ratio for above "scaled log-normal" proposal)
	eta_proposed  = eta  + gsl_ran_ugaussian(GSLrandom)*0.01; /*eta changed by 0.01*/
	//TODO: if(eta_proposed>0.25) eta_proposed=0.25-(eta_proposed-0.25); etc.
	iota_proposed = iota + gsl_ran_ugaussian(GSLrandom)*0.1;
	tc_proposed   = tc   + gsl_ran_ugaussian(GSLrandom)*0.005; /*time changed by 5 ms*/
	phi_proposed  = phi  + gsl_ran_ugaussian(GSLrandom)*0.5;
	ra_proposed   = ra   + gsl_ran_ugaussian(GSLrandom)*0.05;
	dec_proposed  = dec  + gsl_ran_ugaussian(GSLrandom)*0.05;
	psi_proposed  = psi  + gsl_ran_ugaussian(GSLrandom)*0.1;
	//dist_proposed = dist + gsl_ran_ugaussian(GSLrandom)*0.5;
	dist_proposed = dist * exp(gsl_ran_ugaussian(GSLrandom)*0.1); // ~10% change
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
		addVariable(runState->proposalArgs, "logProposalRatio", &logProposalRatio, REAL8_t);
}





