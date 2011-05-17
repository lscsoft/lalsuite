/* Implementation of Nested Sampling for LALInference.
 * (C) John Veitch, 2010
 */

#include "LALInferenceNestedSampler.h"
#include "LALInferencePrior.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <lal/TimeDelay.h>
#include <lalapps.h>

#include <lal/LALStdlib.h>

RCSID("$Id$");
#define PROGRAM_NAME "LALInferenceNestedSampler.c"
#define CVS_ID_STRING "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define CVS_NAME_STRING "$Name$"

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
void calcCVM(gsl_matrix **cvm, LALVariables **Live, UINT4 Nlive)
{
	UINT4 i,j,k;
	UINT4 ND=0;
	LALVariableItem *item,*k_item,*j_item;
	REAL8 *means;

	/* Find the number of dimensions which vary in the covariance matrix */
	for(item=Live[0]->head;item!=NULL;item=item->next)
		if(item->vary==PARAM_LINEAR || item->vary==PARAM_CIRCULAR) ND++;

	/* Set up matrix if necessary */
	if(*cvm==NULL)
	{if(NULL==(*cvm=gsl_matrix_alloc(ND,ND))) {fprintf(stderr,"Unable to allocate matrix memory\n"); exit(1);}}
	else {
		if((*cvm)->size1!=(*cvm)->size2 || (*cvm)->size1!=ND)
		{	fprintf(stderr,"ERROR: Matrix wrong size. Something has gone wrong in calcCVM\n");
			exit(1);
		}
	}
	/* clear the matrix */
	for(i=0;i<(*cvm)->size1;i++) for(j=0;j<(*cvm)->size2;j++) gsl_matrix_set(*cvm,i,j,0.0);

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
			if(j_item->vary!=PARAM_LINEAR && j_item->vary!=PARAM_CIRCULAR) {
				continue;}

			for( k_item=item, k=0; k<=j; k_item=k_item->next ){
				if(k_item->vary!=PARAM_LINEAR && k_item->vary!=PARAM_CIRCULAR) {
					continue;}

					gsl_matrix_set(*cvm,j,k,gsl_matrix_get(*cvm,j,k) +
							   (*(REAL8 *)k_item->value - means[k])*
							   (*(REAL8 *)j_item->value - means[j]));
					k++;
			}
			j++;
		}
	}

	/* Normalise */
	for(i=0;i<ND;i++) for(j=0;j<ND;j++) gsl_matrix_set(*cvm,i,j,gsl_matrix_get(*cvm,i,j)/((REAL8) Nlive));
	free(means);
	/* Fill in variances for circular parameters */
	for(item=Live[0]->head,j=0;item;item=item->next) {
		if(item->vary!=PARAM_CIRCULAR && item->vary!=PARAM_LINEAR) continue;
		if(item->vary==PARAM_CIRCULAR) {
			for(k=0;k<j;k++) gsl_matrix_set(*cvm,j,k,0.0);
			gsl_matrix_set(*cvm,j,j,ang_var(Live,item->name,Nlive));
			for(k=j+1;k<ND;k++) gsl_matrix_set(*cvm,k,j,0.0);
		}
		j++;
	}

	/* the other half */
	for(i=0;i<ND;i++) for(j=0;j<i;j++) gsl_matrix_set(*cvm,j,i,gsl_matrix_get(*cvm,i,j));
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
	REAL8 TOLERANCE=0.5;
	REAL8 logZ,logZnew,logLmin,logLmax=-DBL_MAX,logLtmp,logw,deltaZ,H,logZnoise,dZ=0;
	LALVariables *temp;
	FILE *fpout=NULL;
	gsl_matrix **cvm=calloc(1,sizeof(gsl_matrix *));
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

	/* Set up the proposal scale factor, for use in the multi-student jump step */
	REAL8 propScale=0.1;
	addVariable(runState->proposalArgs,"proposal_scale",&propScale,REAL8_t,PARAM_FIXED);

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
	addVariable(runState->proposalArgs,"LiveCVM",cvm,gslMatrix_t,PARAM_OUTPUT);
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

		UINT4 itercounter=0;
		/* Generate a new live point */
		do{ /* This loop is here in case it is necessary to find a different sample */
			/* Clone an old live point and evolve it */
			while((j=gsl_rng_uniform_int(runState->GSLrandom,Nlive)==minpos)){};
			copyVariables(runState->livePoints[j],runState->currentParams);
			setVariable(runState->algorithmParams,"logLmin",(void *)&logLmin);
			runState->evolve(runState);
			copyVariables(runState->currentParams,runState->livePoints[minpos]);
			logLikelihoods[minpos]=runState->currentLikelihood;
			itercounter++;
		}while(runState->currentLikelihood<=logLmin || *(REAL8 *)getVariable(runState->algorithmParams,"accept_rate")==0.0);

		if (runState->currentLikelihood>logLmax)
			logLmax=runState->currentLikelihood;

		for(j=0;j<Nruns;j++) logwarray[j]+=sample_logt(Nlive,runState->GSLrandom);
		logw=mean(logwarray,Nruns);
		dZ=logadd(logZ,logLmax-((double) iter)/((double)Nlive))-logZ;
		if(verbose) fprintf(stderr,"%i: (%2.1lf%%) accpt: %1.3f H: %3.3lf nats (%3.3lf b) logL:%lf ->%lf logZ: %lf dZ: %lf Zratio: %lf db\n",
									   iter,100.0*((REAL8)iter)/(((REAL8) Nlive)*H),*(REAL8 *)getVariable(runState->algorithmParams,"accept_rate")/(REAL8)itercounter
									   ,H,H/log(2.0),logLmin,runState->currentLikelihood,logZ,dZ,10.0*log10(exp(1.0))*(logZ-*(REAL8 *)getVariable(runState->algorithmParams,"logZnoise")));

		/* Flush output file */
		if(fpout && !(iter%100)) fflush(fpout);
		iter++;
		/* Update the covariance matrix */
		if(!(iter%(Nlive/4))) 	calcCVM(cvm,runState->livePoints,Nlive);
		setVariable(runState->proposalArgs,"LiveCVM",(void *)cvm);
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
	} while(mcmc_iter<Nmcmc);
	destroyVariables(newParams);
	free(newParams);
	REAL8 accept_rate=(REAL8)Naccepted/(REAL8)mcmc_iter;
	setVariable(runState->algorithmParams,"accept_rate",&accept_rate);
	return;
}


void LALInferenceProposalNS(LALInferenceRunState *runState, LALVariables *parameter)
{

	UINT4 nIFO=0;
	LALIFOData *ifo=runState->data;
	REAL8 randnum;
	REAL8 STUDENTTFRAC=0.8,
	      DIFFEVFRAC=0.1,
	      SKYFRAC=0.1;
	top:
	randnum=gsl_rng_uniform(runState->GSLrandom);
	/* Choose a random type of jump to propose */
	if(randnum<STUDENTTFRAC)
		LALInferenceProposalMultiStudentT(runState, parameter);
	else if(randnum<STUDENTTFRAC+DIFFEVFRAC)
		LALInferenceProposalDifferentialEvolution(runState,parameter);
	else if(randnum<STUDENTTFRAC + DIFFEVFRAC + SKYFRAC){
		/* Check number of detectors */
		while(ifo){ifo=ifo->next; nIFO++;}

		if(nIFO<2) goto top;
		if(nIFO<3)
			LALInferenceRotateSky(runState, parameter);
		else {
			/* Choose to rotate or reflect */
			if(randnum- (STUDENTTFRAC + DIFFEVFRAC)>SKYFRAC/2.0)
				LALInferenceRotateSky(runState, parameter);
			else{
				/* Have to call the diff ev too, in case the reflection happens twice */
				LALInferenceReflectDetPlane(runState, parameter);
				LALInferenceProposalDifferentialEvolution(runState,parameter);
			}
		}
	}
	return;
}

UINT4 LALInferenceCheckPositiveDefinite(
						  gsl_matrix       *matrix,
						  UINT4            dim
						  )
{
	gsl_matrix  *m     = NULL;
	gsl_vector  *eigen = NULL;
	gsl_eigen_symm_workspace *workspace = NULL;
	UINT4 i;

	/* copy input matrix */
	m =  gsl_matrix_alloc( dim,dim );
	gsl_matrix_memcpy( m, matrix);

	/* prepare variables */
	eigen = gsl_vector_alloc ( dim );
	workspace = gsl_eigen_symm_alloc ( dim );

	/* compute the eigen values */
	gsl_eigen_symm ( m,  eigen, workspace );

	/* test the result */
	for (i = 0; i < dim; i++)
    {
		/* printf("diag: %f | eigen[%d]= %f\n", gsl_matrix_get( matrix,i,i), i, eigen->data[i]);*/
		if (eigen->data[i]<0)
		{
			printf("NEGATIVE EIGEN VALUE!!! PANIC\n");
			return 0;
		}
	}

	/* freeing unused stuff */
	gsl_eigen_symm_free( workspace);
	gsl_matrix_free(m);
	gsl_vector_free(eigen);

	return 1;
}

/* Reference: http://www.mail-archive.com/help-gsl@gnu.org/msg00631.html*/

void
XLALMultiNormalDeviates(
						REAL4Vector *vector,
						gsl_matrix *matrix,
						UINT4 dim,
						RandomParams *randParam
						)
{
	static LALStatus status;

	UINT4 i=0;
	gsl_matrix *work=NULL;
	gsl_vector *result = NULL;

	static const char *func = "LALMultiNormalDeviates";

	/* check input arguments */
	if (!vector || !matrix || !randParam)
		XLAL_ERROR_VOID( func, XLAL_EFAULT );

	if (dim<1)
		XLAL_ERROR_VOID( func, XLAL_EINVAL );

	/* copy matrix into workspace */
	work =  gsl_matrix_alloc(dim,dim);
	gsl_matrix_memcpy( work, matrix );

	/* compute the cholesky decomposition */
	gsl_linalg_cholesky_decomp(work);

	/* retrieve the normal distributed random numbers (LAL procedure) */
	LALNormalDeviates( &status, vector, randParam);

	/* store this into a gsl vector */
	result = gsl_vector_alloc ( (int)dim );
	for (i = 0; i < dim; i++)
	{
		gsl_vector_set (result, i, vector->data[i]);
	}

	/* compute the matrix-vector multiplication */
	gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, result);

	/* recopy the results */
	for (i = 0; i < dim; i++)
	{
		vector->data[i]=gsl_vector_get (result, i);
	}

	/* free unused stuff */
	gsl_matrix_free(work);
	gsl_vector_free(result);

}


void
XLALMultiStudentDeviates(
						 REAL4Vector  *vector,
						 gsl_matrix   *matrix,
						 UINT4         dim,
						 UINT4         n,
						 RandomParams *randParam
						 )
{
	static const char *func = "LALMultiStudentDeviates";

	static LALStatus status;

	REAL4Vector *dummy=NULL;
	REAL4 chi=0.0, factor;
	UINT4 i;

	/* check input arguments */
	if (!vector || !matrix || !randParam)
		XLAL_ERROR_VOID( func, XLAL_EFAULT );

	if (dim<1)
		XLAL_ERROR_VOID( func, XLAL_EINVAL );

	if (n<1)
		XLAL_ERROR_VOID( func, XLAL_EINVAL );


	/* first draw from MVN */
	XLALMultiNormalDeviates( vector, matrix, dim, randParam);


	/* then draw from chi-square with n degrees of freedom;
     this is the sum d_i*d_i with d_i drawn from a normal
     distribution. */
	LALSCreateVector( &status, &dummy, n);
	LALNormalDeviates( &status, dummy, randParam);

	/* calculate the chisquare distributed value */
	for (i=0; i<n; i++)
	{
		chi+=dummy->data[i]*dummy->data[i];
	}

	/* destroy the helping vector */
	LALSDestroyVector( &status, &dummy );

	/* now, finally, calculate the distribution value */
	factor=sqrt(n/chi);
	for (i=0; i<dim; i++)
	{
		vector->data[i]*=factor;
	}

}


void LALInferenceProposalMultiStudentT(LALInferenceRunState *runState, LALVariables *parameter)
{
	gsl_matrix *covMat=*(gsl_matrix **)getVariable(runState->proposalArgs,"LiveCVM");

	static LALStatus status;

	LALVariableItem *paraHead=NULL;
	REAL4Vector  *step=NULL;
	gsl_matrix *work=NULL;
	REAL8 aii, aij, ajj;
	INT4 i, j, dim;
	RandomParams *randParam;
	UINT4 randomseed = gsl_rng_get(runState->GSLrandom);


	REAL8 proposal_scale=*(REAL8 *)getVariable(runState->proposalArgs,"proposal_scale");
	randParam=XLALCreateRandomParams(randomseed);

	/* set some values */
	dim=covMat->size1;

	/* draw the mutinormal deviates */
	LALSCreateVector( &status, &step, dim);

	/* copy matrix into workspace and scale it appriopriately */
	work =  gsl_matrix_alloc(dim,dim);

	gsl_matrix_memcpy( work, covMat );
	gsl_matrix_scale( work, proposal_scale);

	/* check if the matrix if positive definite */
	while ( !LALInferenceCheckPositiveDefinite( work, dim) ) {
		printf("WARNING: Matrix not positive definite!\n");
		/* downweight the off-axis elements */
		for (i=0; i<dim; ++i)
		{
			for (j=0; j<i; ++j)
			{
				aij=gsl_matrix_get( work, i, j);
				aii=gsl_matrix_get( work, i, i);
				ajj=gsl_matrix_get( work, j, j);

				if ( fabs(aij) > 0.95* sqrt( aii*ajj ) )
				{
					aij=aij/fabs(aij)*0.95*sqrt( aii*ajj );
				}
				gsl_matrix_set( work, i, j, aij);
				gsl_matrix_set( work, j, i, aij);
				printf(" %f", gsl_matrix_get( work, i, j));
			}
			printf("\n");
		}
		exit(0);
	}

	/* draw multivariate student distribution with n=2 */
	XLALMultiStudentDeviates( step, work, dim, 2, randParam);

	/* loop over all parameters */
	for (paraHead=parameter->head,i=0; paraHead; paraHead=paraHead->next)
	{
		/*  if (inputMCMC->verbose)
		 printf("MCMCJUMP: %10s: value: %8.3f  step: %8.3f newVal: %8.3f\n",
		 paraHead->core->name, paraHead->value, step->data[i] , paraHead->value + step->data[i]);*/
		/* only increment the varying parameters, and only increment the data pointer if it's been done*/
		if((paraHead->vary==PARAM_LINEAR || paraHead->vary==PARAM_CIRCULAR))
		/* && strcmp(paraHead->name,"rightascension") && strcmp(paraHead->name,"declination") && strcmp(paraHead->name,"time") */
		{
			*(REAL8 *)paraHead->value += step->data[i];
			i++;
		}
	}

	LALInferenceCyclicReflectiveBound(parameter,runState->priorArgs);
	/* destroy the vectors */
	LALSDestroyVector(&status, &step);
	gsl_matrix_free(work);

	XLALDestroyRandomParams(randParam);
	/* Check boundary condition */

	return;
}

void LALInferenceProposalDifferentialEvolution(LALInferenceRunState *runState,
									   LALVariables *parameter)
	{
		LALVariables **Live=runState->livePoints;
		int i=0,j=0,dim=0,same=1;
		INT4 Nlive = *(INT4 *)getVariable(runState->algorithmParams,"Nlive");
		LALVariableItem *paraHead=NULL;
		LALVariableItem *paraA=NULL;
		LALVariableItem *paraB=NULL;

		dim = parameter->dimension;
		/* Select two other samples A and B*/
		i=gsl_rng_uniform_int(runState->GSLrandom,Nlive);
		/* Draw two different samples from the basket. Will loop back here if the original sample is chosen*/
	drawtwo:
		do {j=gsl_rng_uniform_int(runState->GSLrandom,Nlive);} while(j==i);
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

void GetCartesianPos(REAL8 vec[3],REAL8 longitude, REAL8 latitude)
{
	vec[0]=cos(longitude)*cos(latitude);
	vec[1]=sin(longitude)*cos(latitude);
	vec[1]=sin(latitude);
	return;
}

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

void LALInferenceRotateSky(
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


INT4 LALInferenceReflectDetPlane(
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
