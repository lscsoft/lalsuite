/* Implementation of Nested Sampling for LALInference.
 * (C) John Veitch, 2010
 */

#include "LALInferenceNestedSampler.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <lal/LALStdlib.h>

double logadd(double a,double b){
	if(a>b) return(a+log(1.0+exp(b-a)));
	else return(b+log(1.0+exp(a-b)));
}


REAL8 mean(REAL8 *array,int N){
	REAL8 sum=0.0;
	int i;
	for(i=0;i<N;i++) sum+=array[i];
	return sum/((REAL8) N);
}

REAL8 sample_logt(int Nlive){
	REAL8 t=0.0;
	REAL8 a=0.0;
	while((Nlive--)>1) {a=gsl_rng_uniform(RNG); t = t>a ? t : a;}
	return(log(t));
}


/* NestedSamplingAlgorithm implements the nested sampling algorithm,
 see e.g. Sivia "Data Analysis: A Bayesian Tutorial, 2nd edition.
 REQUIREMENTS:
	Calling routine must have set up runState->livePoints already to
	contain samples from the prior distribution.
	runState->algorithmParams must contain a variable "logLikelihoods"
	which contains a REAL8 array of likelihood values for the live
	points.
 */
void NestedSamplingAlgorithm(LALInferenceRunState *runState)
{
	UINT4 iter=0,i,j,minpos;
	UINT4 Nlive=*(UINT4 *)getVariable(runState->algorithmParams,"Nlive");
	UINT4 Nruns=1;
	REAL8 *logZarray,*oldZarray,*Harray,*logwarray,*Wtarray;
	REAL8 TOLERANCE=0.1;
	REAL8 logZ,logZnew,logLmax=-DBL_MAX,logLtmp;
	LALVariables *temp;
	FILE *fp=NULL;

	/* Operate on parallel runs if requested */
	if(checkVariable(runState->algorithmParams,"Nruns"))
		Nruns = *(UINT4 *) getVariable(runState->algorithmParams,"Nruns");

	if(checkVariable(runState->algorithmParams,"tolerance"))
		TOLERANCE = *(REAL8 *) getVariable(runState->algorithmParams,"tolerance");

	/* Check that necessary parameters are created */
	if(!checkVariable(runState->algorithmParams,"logLmin"))
		addVariable(runState->algorithmParams,"logLmin",-DBL_MAX,REAL8_t,PARAM_OUTPUT);

	if(!checkVariable(runState->algorithmParams,"accept_rate"))
		addVariable(runState->algorithmParams,"accept_rate",0.0,REAL8_t,PARAM_OUTPUT);

	
	/* FIXME: Open output file */
	char *outfile=getProcParamVal(runState->commandLine,"outfile")->value;
	fp=open(outfile,"w");
	if(fpout==NULL) fprintf(stderr,"Unable to open output file %s!\n",outfile);

	/* Set up arrays for parallel runs */
	logZarray = calloc(Nruns,sizeof(REAL8));
	oldZarray = calloc(Nruns,sizeof(REAL8));
	Harray = calloc(Nruns,sizeof(REAL8));
	logwarray = calloc(Nruns,sizeof(REAL8));
	Wtarray = calloc(Nruns,sizeof(REAL8));
	if(logZarray==NULL || Harray==NULL || oldZarray==NULL || logwarray==NULL)
		{fprintf(stderr,"Unable to allocate RAM\n"); exit(-1);}

	logw=log(1.0-exp(-1.0/Nlive));
	for(i=0;i<Nruns;i++)  {logwarray[i]=logw; logZarray[i]=-DBL_MAX; oldZarray[i]=-DBL_MAX; Harray[i]=0.0;}
	i=0;
	/* Find maximum likelihood */
	for(i=0;i<Nlive;i++)
	{
		logLtmp=(REAL8 *)getVariable(runState->algorithmParams,"logLikelihoods")[i];
		logLmax=logLtmp>logLmax? logLtmp : logLmax;
	}

	/* Iterate until termination condition is met */
	do {
		/* Find minimum likelihood sample to replace */
		minpos=0;
		for(i=0;i<Nlive;i++){
			if((REAL8 *)getVariable(runState->algorithmParams,"logLikelihoods")[i]
			   <(REAL8 *)getVariable(runState->algorithmParams,"logLikelihoods")[minpos])
				minpos=i;
		}
		logLmin=(REAL8 *)getVariable(runState->algorithmParams,"logLikelihoods")[minpos];

		/* Update evidence array */
		for(j=0;j<Nruns;j++){
			logZarray[j]=logadd(logZarray[j],Live[minpos]->logLikelihood + logwarray[j]);
			Wtarray[j]=logwarray[j]+Live[minpos]->logLikelihood;
			Harray[j]= exp(Wtarray[j]-logZarray[j])*Live[minpos]->logLikelihood
			+ exp(oldZarray[j]-logZarray[j])*(Harray[j]+oldZarray[j])-logZarray[j];
		}
		logZnew=mean(logZarray,Nruns);
		deltaZ=logZnew-logZ;
		H=mean(Harray,Nruns);
		logZ=logZnew;
		for(j=0;j<Nruns;j++) oldZarray[j]=logZarray[j];

		/* Write out old sample */
		fprintSample(fpout,runState->livePoints[i]);
		fprintf(fpout,"%lf\n",(REAL8 *)getVariable(runState->algorithmParams,"logLikelihoods")[i]);


		/* Generate a new live point */
		do{ /* This loop is here in case it is necessary to find a different sample */
			/* Clone an old live point and evolve it */
			while((j=gsl_rng_uniform_int(runState->GSLrandom,Nlive)==minpos)){};
			copyVariables(runState->livePoints[j],runState->currentParams);
			setVariable(runState->algorithmParams,"logLmin",(void *)&logLmin);
			runState->evolve(runState);
			copyVariables(runState->currentParams,runState->livePoints[minpos]);
			(REAL8 *)getVariable(runState->algorithmParams,"logLikelihoods")[minpos]=runState->currentLikelihood;
		}while(runState->currentLikelihood<=logLmin);

		if (runState->currentLikelihood>logLmax)
			logLmax=runState->currentLikelihood;

		for(j=0;j<Nruns;j++) logwarray[j]+=sample_logt(Nlive);
		logw=mean(logwarray,Nruns);
		if(MCMCinput->verbose) fprintf(stderr,"%i: (%2.1lf%%) accpt: %1.3f H: %3.3lf nats (%3.3lf b) logL:%lf ->%lf logZ: %lf Zratio: %lf db\n",
									   iter,100.0*((REAL8)iter)/(((REAL8) Nlive)*H),*(REAL8 *)getVariable(runState->algorithmParams,"accept_rate"),
									   ,H,H/log(2.0),logLmin,runState->currentLikelihood,logZ,10.0*log10(exp(1.0))*(logZ-*(REAL8 *)getVariable(runState->algorithmParams,"logZnoise")));

		/* Flush output file */
		if(fpout && !(iter%100)) fflush(fpout);
		iter++;
	}
	while(iter<Nlive || logadd(logZ,logLmax-((double iter)/(double)Nlive))-logZ > TOLERANCE)

	/* Sort the remaining points (not essential, just nice)*/
		for(i=0;i<Nlive-1;i++){
			minpos=i;
			logLmin=(REAL8 *)getVariable(runState->algorithmParams,"logLikelihoods")[i];
			for(j=i+1;j<Nlive;j++){
				if((REAL8 *)getVariable(runState->algorithmParams,"logLikelihoods")[j]<logLmin)
					{
						minpos=j;
						logLmin=(REAL8 *)getVariable(runState->algorithmParams,"logLikelihoods")[j];
					}
			}
			temp=runState->livePoints[minpos]; /* Put the minimum remaining point in the current position */
			runState->livePoints[minpos]=runState->livePoints[i];
			runState->livePoints[i]=temp;
		}

	/* final corrections */
	for(i=0;i<Nlive;i++){
		logZ=logadd(logZ,(REAL8 *)getVariable(runState->algorithmParams,"logLikelihoods")[i]+logw);
		for(j=0;j<Nruns;j++){
			logwarray[j]+=sample_logt(Nlive);
			logZarray[j]=logadd(logZarray[j],(REAL8 *)getVariable(runState->algorithmParams,"logLikelihoods")[i]+logwarray[j]);
		}

		fprintSample(fpout,runState->livePoints[i]);
		fprintf(fpout,"%lf\n",(REAL8 *)getVariable(runState->algorithmParams,"logLikelihoods")[i]);
	}

	/* Write out the evidence */

}

/* Evolve nested sampling algorithm by one step, i.e.
 evolve runState->currentParams to a new point with higher
 likelihood than currentLikelihood. Uses the MCMC method.
 */
NestedSamplingOneStep(LALInferenceRunState runState)
{
	LALVariables *newParams;
	UINT4 mcmc_iter=0,accept=0,Naccepted=0;
	UINT4 Nmcmc=*(UINT4 *)getVariable(runState->algorithmParams,"Nmcmc");
	REAL8 logLmin=*(REAL8 *)getVariable(runState->algorithmParams,"logLmin");
	REAL8 logPriorOld,logPriorNew,logLnew;
	/* Make a copy of the parameters passed through currentParams */
	copyVariables(runState->currentParams,newParams);
	/* Evolve the sample until it is accepted */
	logPriorOld=runState->prior(runState,newParams);
	do{
		mcmc_iter++;
		runState->proposal(runState,newParams);
		logPriorNew=runState->prior(runState,newParams);
		/* If rejected, continue to next iteration */
		if(log(gsl_rng_uniform(runState->GSLrandom))>logPriorNew-logPriorOld)
			continue;
		/* Otherwise, check that logL is OK */
		logLnew=runState->likelihood(newParams,runState->data,runState->template);
		if(logLnew>logLmin){
			Naccepted++;
			logPriorOld=logPriorNew;
			copyVariables(newParams,runState->currentParams);
			runState->currentLikelihood=logLnew;
	} while(runState->currentLikelihood<=logLmin || mcmc_iter<Nmcmc);
	destroyVariables(newParams);
	setVariable(runState->algorithmParams,"accept_rate",(REAL8)Naccepted/(REAL8)mcmc_iter);
	return;
}