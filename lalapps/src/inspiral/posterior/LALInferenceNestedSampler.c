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

REAL8 sample_logt(int Nlive,gsl_rng *RNG){
	REAL8 t=0.0;
	REAL8 a=0.0;
	while((Nlive--)>1) {a=gsl_rng_uniform(RNG); t = t>a ? t : a;}
	return(log(t));
}

/* Calculate shortest angular distance between a1 and a2 */
REAL8 ang_dist(REAL8 a1, REAL8 a2){
	double raw = (a2>a1 ? a2-a1 : a1-a2);
	return(raw>LAL_PI ? 2.0*LAL_PI - raw : raw);
}

/* Calculate the variance of a modulo-2pi distribution */
REAL8 ang_var(LALVariables **list,const char *pname, int N){
	int i=0;
	REAL8 ang_mean=0.0;
	REAL8 var=0.0;
	REAL8 ms,mc;
	/* Calc mean */
	for(i=0,ms=0.0,mc=0.0;i<N;i++) {
		ms+=sin(*(REAL8 *)getVariable(list[i],pname));
		mc+=cos(*(REAL8 *)getVariable(list[i],pname));
	}
	ms/=N; mc/=N;
	ang_mean=atan2(ms,mc);
	ang_mean = ang_mean<0? 2.0*LAL_PI + ang_mean : ang_mean;
	/* calc variance */
	for(i=0;i<N;i++) var+=ang_dist(*(REAL8 *)getVariable(list[i],pname),ang_mean)*ang_dist(*(REAL8 *)getVariable(list[i],pname),ang_mean);
	return(var/(REAL8)N);
}

/* estimateCovarianceMatrix reads the list of live points,
 and works out the covariance matrix of the varying parameters
 with varyType==PARAM_LINEAR */
void calcCVM(gsl_matrix *cvm, LALVariables **Live, UINT4 Nlive)
{
	UINT4 i,j,k;
	UINT4 ND=0;
	LALVariableItem *item,*k_item,*j_item;
	REAL8 *means;
	
	/* Find the number of dimensions which vary in the covariance matrix */
	for(item=Live[0]->head;item!=NULL;item=item->next)
		if(item->vary==PARAM_LINEAR || item->vary==PARAM_CIRCULAR) ND++;
	
	/* Set up matrix if necessary */
	if(cvm==NULL)
	{if(NULL==(cvm=gsl_matrix_alloc(ND,ND))) {fprintf(stderr,"Unable to allocate matrix memory\n"); exit(1);}}
	else {
		if(cvm->size1!=cvm->size2 || cvm->size1!=ND)
		{	fprintf(stderr,"ERROR: Matrix wrong size. Something has gone wrong in calcCVM\n");
			exit(1);
		}
	}
	/* clear the matrix */
	for(i=0;i<cvm->size1;i++) for(j=0;j<cvm->size2;j++) gsl_matrix_set(cvm,i,j,0.0);

	/* Find the means */
	if(NULL==(means = malloc((size_t)ND*sizeof(REAL8)))){fprintf(stderr,"Can't allocate RAM"); exit(-1);}
	for(i=0;i<ND;i++) means[i]=0.0;
	for(i=0;i<Nlive;i++){
		for(item=Live[i]->head,j=0;item;item=item->next) {
			if(item->vary==PARAM_LINEAR || item->vary==PARAM_CIRCULAR ) {
				if (item->type==REAL4_t) means[j]+=*(REAL4 *)item->value;
				if (item->type==REAL8_t) means[j]+=*(REAL8 *)item->value;
				j++;
			}
		}
	}
	for(j=0;j<ND;j++) means[j]/=(REAL8)Nlive;
	
	/* Find the (co)-variances */
	for(i=0;i<Nlive;i++){
		k_item = j_item = item = Live[i]->head;
		for( j_item=item,j=0; j_item; j_item=j_item->next ){
			if(j_item->vary!=PARAM_LINEAR || item->vary!=PARAM_CIRCULAR ) continue; /* Skip params that we don't vary */
			for( k_item=item,k=0; k<=j; k_item=k_item->next ){
				if(k_item->vary!=PARAM_LINEAR || k_item->vary!=PARAM_CIRCULAR) continue;
				/* Only set elements for non-circular parameters */
				if(k_item->vary==PARAM_LINEAR && j_item->vary==PARAM_LINEAR)
					gsl_matrix_set(cvm,j,k,gsl_matrix_get(cvm,j,k) +
							   (*(REAL8 *)k_item->value - means[k])*
							   (*(REAL8 *)j_item->value - means[j]));
				k++;
			}
			j++;
		}
	}

	/* Normalise */
	for(i=0;i<ND;i++) for(j=0;j<ND;j++) gsl_matrix_set(cvm,i,j,gsl_matrix_get(cvm,i,j)/((REAL8) Nlive));
	free(means);
	/* Fill in variances for circular parameters */
	for(item=Live[0]->head,j=0;j<ND;j++,item=item->next) {
		if(item->vary==PARAM_CIRCULAR) {
			for(k=0;k<j;k++) gsl_matrix_set(cvm,j,k,0.0);
			gsl_matrix_set(cvm,j,j,ang_var(Live,item->name,Nlive));
			for(k=j+1;k<ND;k++) gsl_matrix_set(cvm,k,j,0.0);
		}
	}
	
	/* the other half */
	for(i=0;i<ND;i++) for(j=0;j<i;j++) gsl_matrix_set(cvm,j,i,gsl_matrix_get(cvm,i,j));
	return;
}


/* NestedSamplingAlgorithm implements the nested sampling algorithm,
 see e.g. Sivia & Skilling "Data Analysis: A Bayesian Tutorial, 2nd edition.
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
	REAL8 logZ,logZnew,logLmin,logLmax=-DBL_MAX,logLtmp,logw,deltaZ,H,logZnoise,dZ=0;
	LALVariables *temp;
	FILE *fpout=NULL;
	gsl_matrix *cvm=NULL;
	REAL8 dblmax=-DBL_MAX;
	REAL8 zero=0.0;
	REAL8 *logLikelihoods=NULL;
	UINT4 verbose=0;
	
	logZnoise=NullLogLikelihood(runState->data);
	addVariable(runState->algorithmParams,"logZnoise",&logZnoise,REAL8_t,PARAM_FIXED);
	logLikelihoods=(REAL8 *)(*(REAL8Vector **)getVariable(runState->algorithmParams,"logLikelihoods"))->data;

	verbose=checkVariable(runState->algorithmParams,"verbose");
	
	/* Operate on parallel runs if requested */
	if(checkVariable(runState->algorithmParams,"Nruns"))
		Nruns = *(UINT4 *) getVariable(runState->algorithmParams,"Nruns");

	if(checkVariable(runState->algorithmParams,"tolerance"))
		TOLERANCE = *(REAL8 *) getVariable(runState->algorithmParams,"tolerance");

	/* Check that necessary parameters are created */
	if(!checkVariable(runState->algorithmParams,"logLmin"))
		addVariable(runState->algorithmParams,"logLmin",&dblmax,REAL8_t,PARAM_OUTPUT);

	if(!checkVariable(runState->algorithmParams,"accept_rate"))
		addVariable(runState->algorithmParams,"accept_rate",&zero,REAL8_t,PARAM_OUTPUT);

	
	/* Open output file */
	ProcessParamsTable *ppt=getProcParamVal(runState->commandLine,"--outfile");
	if(!ppt){
		fprintf(stderr,"Must specify --outfile <filename.dat>\n");
		exit(1);
	}
	char *outfile=ppt->value;
	fpout=fopen(outfile,"w");
	if(fpout==NULL) fprintf(stderr,"Unable to open output file %s!\n",outfile);

	/* Set up arrays for parallel runs */
	logZarray = calloc(Nruns,sizeof(REAL8));
	oldZarray = calloc(Nruns,sizeof(REAL8));
	Harray = calloc(Nruns,sizeof(REAL8));
	logwarray = calloc(Nruns,sizeof(REAL8));
	Wtarray = calloc(Nruns,sizeof(REAL8));
	if(logZarray==NULL || Harray==NULL || oldZarray==NULL || logwarray==NULL || Wtarray==NULL)
		{fprintf(stderr,"Unable to allocate RAM\n"); exit(-1);}

	logw=log(1.0-exp(-1.0/Nlive));
	for(i=0;i<Nruns;i++)  {logwarray[i]=logw; logZarray[i]=-DBL_MAX; oldZarray[i]=-DBL_MAX; Harray[i]=0.0;}
	i=0;
	/* Find maximum likelihood */
	for(i=0;i<Nlive;i++)
	{
		logLtmp=logLikelihoods[i];
		logLmax=logLtmp>logLmax? logLtmp : logLmax;
	}
	/* Add the covariance matrix for proposal distribution */
	calcCVM(cvm,runState->livePoints,Nlive);
	addVariable(runState->proposalArgs,"LiveCVM",&cvm,gslMatrix_t,PARAM_OUTPUT);
	fprintf(stdout,"Starting nested sampling loop!\n");
	/* Iterate until termination condition is met */
	do {
		/* Find minimum likelihood sample to replace */
		minpos=0;
		for(i=1;i<Nlive;i++){
			if(logLikelihoods[i]<logLikelihoods[minpos])
				minpos=i;
		}
		logLmin=logLikelihoods[minpos];

		/* Update evidence array */
		for(j=0;j<Nruns;j++){
			logZarray[j]=logadd(logZarray[j],logLikelihoods[minpos]+ logwarray[j]);
			Wtarray[j]=logwarray[j]+logLikelihoods[minpos];
			Harray[j]= exp(Wtarray[j]-logZarray[j])*logLikelihoods[minpos]
			+ exp(oldZarray[j]-logZarray[j])*(Harray[j]+oldZarray[j])-logZarray[j];
		}
		logZnew=mean(logZarray,Nruns);
		deltaZ=logZnew-logZ;
		H=mean(Harray,Nruns);
		logZ=logZnew;
		for(j=0;j<Nruns;j++) oldZarray[j]=logZarray[j];

		/* Write out old sample */
		fprintSample(fpout,runState->livePoints[minpos]);
		fprintf(fpout,"%lf\n",logLikelihoods[minpos]);


		/* Generate a new live point */
		do{ /* This loop is here in case it is necessary to find a different sample */
			/* Clone an old live point and evolve it */
			while((j=gsl_rng_uniform_int(runState->GSLrandom,Nlive)==minpos)){};
			copyVariables(runState->livePoints[j],runState->currentParams);
			setVariable(runState->algorithmParams,"logLmin",(void *)&logLmin);
			runState->evolve(runState);
			copyVariables(runState->currentParams,runState->livePoints[minpos]);
			logLikelihoods[minpos]=runState->currentLikelihood;
		}while(runState->currentLikelihood<=logLmin);

		if (runState->currentLikelihood>logLmax)
			logLmax=runState->currentLikelihood;

		for(j=0;j<Nruns;j++) logwarray[j]+=sample_logt(Nlive,runState->GSLrandom);
		logw=mean(logwarray,Nruns);
		dZ=logadd(logZ,logLmax-((double) iter)/((double)Nlive))-logZ;
		if(verbose) fprintf(stderr,"%i: (%2.1lf%%) accpt: %1.3f H: %3.3lf nats (%3.3lf b) logL:%lf ->%lf logZ: %lf dZ: %lf Zratio: %lf db\n",
									   iter,100.0*((REAL8)iter)/(((REAL8) Nlive)*H),*(REAL8 *)getVariable(runState->algorithmParams,"accept_rate")
									   ,H,H/log(2.0),logLmin,runState->currentLikelihood,logZ,dZ,10.0*log10(exp(1.0))*(logZ-*(REAL8 *)getVariable(runState->algorithmParams,"logZnoise")));

		/* Flush output file */
		if(fpout && !(iter%100)) fflush(fpout);
		iter++;
		/* Update the covariance matrix */
		if(iter%(Nlive/4)) 	calcCVM(cvm,runState->livePoints,Nlive);
		setVariable(runState->proposalArgs,"LiveCVM",&cvm);
	}
	while(iter<Nlive ||  dZ> TOLERANCE);

	/* Sort the remaining points (not essential, just nice)*/
		for(i=0;i<Nlive-1;i++){
			minpos=i;
			logLmin=logLikelihoods[i];
			for(j=i+1;j<Nlive;j++){
				if(logLikelihoods[j]<logLmin)
					{
						minpos=j;
						logLmin=logLikelihoods[j];
					}
			}
			temp=runState->livePoints[minpos]; /* Put the minimum remaining point in the current position */
			runState->livePoints[minpos]=runState->livePoints[i];
			runState->livePoints[i]=temp;
		}

	/* final corrections */
	for(i=0;i<Nlive;i++){
		logZ=logadd(logZ,logLikelihoods[i]+logw);
		for(j=0;j<Nruns;j++){
			logwarray[j]+=sample_logt(Nlive,runState->GSLrandom);
			logZarray[j]=logadd(logZarray[j],logLikelihoods[i]+logwarray[j]);
		}

		fprintSample(fpout,runState->livePoints[i]);
		fprintf(fpout,"%lf\n",logLikelihoods[i]);
	}

	/* Write out the evidence */
	fclose(fpout);
	char bayesfile[FILENAME_MAX];
	sprintf(bayesfile,"%s_B.txt",outfile);
	fpout=fopen(bayesfile,"w");
	fprintf(fpout,"%lf %lf %lf %lf\n",logZ-logZnoise,logZ,logZnoise,logLmax);
	fclose(fpout);
}

/* Evolve nested sampling algorithm by one step, i.e.
 evolve runState->currentParams to a new point with higher
 likelihood than currentLikelihood. Uses the MCMC method.
 */
void NestedSamplingOneStep(LALInferenceRunState *runState)
{
	LALVariables *newParams=NULL;
	UINT4 mcmc_iter=0,Naccepted=0;
	UINT4 Nmcmc=*(UINT4 *)getVariable(runState->algorithmParams,"Nmcmc");
	REAL8 logLmin=*(REAL8 *)getVariable(runState->algorithmParams,"logLmin");
	REAL8 logPriorOld,logPriorNew,logLnew;
	newParams=calloc(1,sizeof(LALVariables));
	/* Make a copy of the parameters passed through currentParams */
	copyVariables(runState->currentParams,newParams);
	/* Evolve the sample until it is accepted */
	logPriorOld=runState->prior(runState,runState->currentParams);
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
		}
	} while(runState->currentLikelihood<=logLmin || mcmc_iter<Nmcmc);
	destroyVariables(newParams);
	free(newParams);
	REAL8 accept_rate=(REAL8)Naccepted/(REAL8)mcmc_iter;
	setVariable(runState->algorithmParams,"accept_rate",&accept_rate);
	return;
}


void LALInferenceProposalNS(LALInferenceRunState *runState, LALVariables *parameter)
{
	LALInferenceProposalDifferentialEvolution(runState,parameter);
	
	return;	
}

void LALInferenceProposalMultiStudentT(LALInferenceRunState *runState, LALVariables *parameter)
{
	return;
}

void LALInferenceProposalDifferentialEvolution(LALInferenceRunState *runState,
									   LALVariables *parameter)
	{
		LALVariables **Live=runState->livePoints;
		int i=0,j=0,dim=0,same=1;
		REAL4 randnum;
		INT4 Nlive = *(INT4 *)getVariable(runState->algorithmParams,"Nlive");
		LALVariableItem *paraHead=NULL;
		LALVariableItem *paraA=NULL;
		LALVariableItem *paraB=NULL;
		
		dim = parameter->dimension;
		/* Select two other samples A and B*/
		randnum=gsl_rng_uniform(runState->GSLrandom);
		i=(int)(Nlive*randnum);
		/* Draw two different samples from the basket. Will loop back here if the original sample is chosen*/
	drawtwo:
		do {randnum=gsl_rng_uniform(runState->GSLrandom); j=(int)(Nlive*randnum);} while(j==i);
		paraHead=parameter->head;
		paraA=Live[i]->head; paraB=Live[j]->head;
		/* Add the vector B-A */
		same=1;
		for(paraHead=parameter->head,paraA=Live[i]->head,paraB=Live[j]->head;paraHead;paraHead=paraHead->next,paraB=paraB->next,paraA=paraA->next)
		{
			if(paraHead->vary!=PARAM_LINEAR && paraHead->vary!=PARAM_CIRCULAR) continue;
			*(REAL8 *)paraHead->value+=*(REAL8 *)paraB->value;
			*(REAL8 *)paraHead->value-=*(REAL8 *)paraA->value;
			if(*(REAL8 *)paraHead->value!=*(REAL8 *)paraA->value &&
			   *(REAL8 *)paraHead->value!=*(REAL8 *)paraB->value &&
			   *(REAL8 *)paraA->value!=*(REAL8 *)paraB->value) same=0;
		}
		if(same==1) goto drawtwo;
		/* Bring the sample back into bounds */
		LALInferenceCyclicReflectiveBound(parameter,runState->priorArgs);
		return;
	}

void LALInferenceCyclicReflectiveBound(LALVariables *parameter, LALVariables *priorArgs){
/* Apply cyclic and reflective boundaries to parameter to bring it back within
 the prior */
	LALVariableItem *paraHead=NULL;
	REAL8 delta;
	REAL8 min,max;
	for (paraHead=parameter->head;paraHead;paraHead=paraHead->next)
	{
		if(paraHead->vary==PARAM_FIXED || paraHead->vary==PARAM_OUTPUT) continue;
		getMinMaxPrior(priorArgs,paraHead->name, (void *)&min, (void *)&max);
		if(paraHead->vary==PARAM_CIRCULAR) /* For cyclic boundaries */
		{
			delta = max-min;
			while ( *(REAL8 *)paraHead->value > max) 
				*(REAL8 *)paraHead->value -= delta;
			while ( *(REAL8 *)paraHead->value < min) 
				*(REAL8 *)paraHead->value += delta;
		}
		else if(paraHead->vary==PARAM_LINEAR) /* Use reflective boundaries */
		{
			if(max < *(REAL8 *)paraHead->value) *(REAL8 *)paraHead->value-=2.0*(*(REAL8 *)paraHead->value - max);
			if(min > *(REAL8 *)paraHead->value) *(REAL8 *)paraHead->value+=2.0*(min - *(REAL8 *)paraHead->value);
		}
	}	
	return;
}