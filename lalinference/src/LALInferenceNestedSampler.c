/* Implementation of Nested Sampling for LALInference.
 * (C) John Veitch, 2010
 */

#include <lal/LALInferenceNestedSampler.h>
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceProposal.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <lal/TimeDelay.h>
#include <lal/LALInferenceConfig.h>

#include <lal/LALStdlib.h>

#ifdef HAVE_LIBLALXML
#include <lal/LALInferenceXML.h>
#endif

#define PROGRAM_NAME "LALInferenceNestedSampler.c"
#define CVS_ID_STRING "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define CVS_NAME_STRING "$Name$"


 static void LALInferenceProjectSampleOntoEigenvectors(LALInferenceVariables *params, gsl_matrix *eigenvectors, REAL8Vector **projection);
/* Prototypes for private "helper" functions. */
//static void SamplePriorDiscardAcceptance(LALInferenceRunState *runState);
static void crossProduct(REAL8 out[3],REAL8 x[3],REAL8 y[3]);
static void CartesianToSkyPos(REAL8 pos[3],REAL8 *longitude, REAL8 *latitude);
static void GetCartesianPos(REAL8 vec[3],REAL8 longitude, REAL8 latitude);
static double logadd(double a,double b);
static REAL8 mean(REAL8 *array,int N);
static void getMinMaxLivePointValue( LALInferenceVariables **livepoints, 
                                     const CHAR *pname, UINT4 Nlive, 
                                     REAL8 *minval, REAL8 *maxval );

static double logadd(double a,double b){
	if(a>b) return(a+log(1.0+exp(b-a)));
	else return(b+log(1.0+exp(a-b)));
}


static REAL8 mean(REAL8 *array,int N){
	REAL8 sum=0.0;
	int i;
	for(i=0;i<N;i++) sum+=array[i];
	return sum/((REAL8) N);
}

/* Calculate shortest angular distance between a1 and a2 */
REAL8 LALInferenceAngularDistance(REAL8 a1, REAL8 a2){
	double raw = (a2>a1 ? a2-a1 : a1-a2);
	return(raw>LAL_PI ? 2.0*LAL_PI - raw : raw);
}

/** Get the maximum value of a parameter from a set of live points */
static void getMinMaxLivePointValue( LALInferenceVariables **livepoints, 
                                     const CHAR *pname, UINT4 Nlive, 
                                     REAL8 *minval, REAL8 *maxval ){
  REAL8 maxvaltmp = -DBL_MAX, minvaltmp = DBL_MAX;
  UINT4 i = 0;
  
  for ( i = 0; i < Nlive; i++ ){
    REAL8 val = *(REAL8 *)LALInferenceGetVariable( livepoints[i], pname );
    
    if ( val < minvaltmp ) minvaltmp = val;
    if ( val > maxvaltmp ) maxvaltmp = val;
  }
  
  *minval = minvaltmp;
  *maxval = maxvaltmp;
  return;
}
  

 
/* Calculate the variance of a modulo-2pi distribution */
REAL8 LALInferenceAngularVariance(LALInferenceVariables **list,const char *pname, int N){
	int i=0;
	REAL8 ang_mean=0.0;
	REAL8 var=0.0;
	REAL8 ms,mc;
	/* Calc mean */
	for(i=0,ms=0.0,mc=0.0;i<N;i++) {
		ms+=sin(*(REAL8 *)LALInferenceGetVariable(list[i],pname));
		mc+=cos(*(REAL8 *)LALInferenceGetVariable(list[i],pname));
	}
	ms/=N; mc/=N;
	ang_mean=atan2(ms,mc);
	ang_mean = ang_mean<0? 2.0*LAL_PI + ang_mean : ang_mean;
	/* calc variance */
	for(i=0;i<N;i++) var+=LALInferenceAngularDistance(*(REAL8 *)LALInferenceGetVariable(list[i],pname),ang_mean)*LALInferenceAngularDistance(*(REAL8 *)LALInferenceGetVariable(list[i],pname),ang_mean);
	return(var/(REAL8)N);
}

REAL8 LALInferenceNSSample_logt(int Nlive,gsl_rng *RNG){
	REAL8 t=0.0;
	REAL8 a=0.0;
	while((Nlive--)>1) {a=gsl_rng_uniform(RNG); t = t>a ? t : a;}
	return(log(t));
}

/* estimateCovarianceMatrix reads the list of live points,
 and works out the covariance matrix of the varying parameters
 - CIRCULAR parameters are wrapped around before the calculation and must be
 scaled to the range 0 -> 2pi */
void LALInferenceNScalcCVM(gsl_matrix **cvm, LALInferenceVariables **Live, UINT4 Nlive)
{
	UINT4 i,j,k;
	UINT4 ND=0;
	LALInferenceVariableItem *item,*k_item,*j_item;
	REAL8 *means, *ms, *mc;
	
	/* Find the number of dimensions which vary in the covariance matrix */
	for(item=Live[0]->head;item!=NULL;item=item->next)
		if(item->vary==LALINFERENCE_PARAM_LINEAR || item->vary==LALINFERENCE_PARAM_CIRCULAR) ND++;
	
	/* Set up matrix if necessary */
	if(*cvm==NULL)
	{if(NULL==(*cvm=gsl_matrix_alloc(ND,ND))) {fprintf(stderr,"Unable to allocate matrix memory\n"); exit(1);}}
	else {
		if((*cvm)->size1!=(*cvm)->size2 || (*cvm)->size1!=ND)
		{	fprintf(stderr,"ERROR: Matrix wrong size. Something has gone wrong in LALInferenceNScalcCVM\n");
			exit(1);
		}
	}
	/* clear the matrix */
	for(i=0;i<(*cvm)->size1;i++) for(j=0;j<(*cvm)->size2;j++) gsl_matrix_set(*cvm,i,j,0.0);

	/* Find the means */
	if(NULL==(means = malloc((size_t)ND*sizeof(REAL8)))){fprintf(stderr,"Can't allocate RAM"); exit(-1);}
	if(NULL==(ms = malloc((size_t)ND*sizeof(REAL8)))){fprintf(stderr,"Can't allocate RAM"); exit(-1);}
	if(NULL==(mc = malloc((size_t)ND*sizeof(REAL8)))){fprintf(stderr,"Can't allocate RAM"); exit(-1);}
	for(i=0;i<ND;i++){ 
          means[i]=0.0;
          ms[i] = 0.;
          mc[i] = 0.;
        }
	for(i=0;i<Nlive;i++){
                   for(item=Live[i]->head,j=0;item;item=item->next) {
			/*if(item->vary==LALINFERENCE_PARAM_LINEAR || item->vary==LALINFERENCE_PARAM_CIRCULAR ) {
				if (item->type==LALINFERENCE_REAL4_t) means[j]+=*(REAL4 *)item->value;
				if (item->type==LALINFERENCE_REAL8_t) means[j]+=*(REAL8 *)item->value;
				j++;
			}*/
                        if(item->vary==LALINFERENCE_PARAM_LINEAR ) {
                                if (item->type==LALINFERENCE_REAL4_t) means[j]+=*(REAL4 *)item->value;
                                if (item->type==LALINFERENCE_REAL8_t) means[j]+=*(REAL8 *)item->value;
                                j++;
                        }
			else if( item->vary==LALINFERENCE_PARAM_CIRCULAR ){
                               if(item->type==LALINFERENCE_REAL4_t){
                                 ms[j] += sin(*(REAL4 *)item->value);
                                 mc[j] += cos(*(REAL4 *)item->value);
                               }
                               if(item->type==LALINFERENCE_REAL8_t){
                                 ms[j] += sin(*(REAL8 *)item->value);
                                 mc[j] += cos(*(REAL8 *)item->value);
                               }

                               j++;
                        }
		}
	}

        for(item=Live[0]->head,j=0;item;item=item->next){ 
          if( item->vary==LALINFERENCE_PARAM_LINEAR ){
            means[j]/=(REAL8)Nlive;
            j++;
          }
          if( item->vary==LALINFERENCE_PARAM_CIRCULAR ){
            ms[j]/=(REAL8)Nlive;
            mc[j]/=(REAL8)Nlive;
            
            means[j] = atan2(ms[j], mc[j]);
            means[j] = means[j]<0? 2.0*LAL_PI + means[j] : means[j];
            
            j++;
          }
        }
        
        free(ms);
        free(mc);
        
	/* Find the (co)-variances */
	for(i=0;i<Nlive;i++){
		k_item = j_item = item = Live[i]->head;

		for( j_item=item,j=0; j_item; j_item=j_item->next ){
		  REAL8 jval = 0.;	
                  if(j_item->vary!=LALINFERENCE_PARAM_LINEAR && j_item->vary!=LALINFERENCE_PARAM_CIRCULAR) {
				continue;}
			
			if( j_item->vary==LALINFERENCE_PARAM_CIRCULAR )
                          jval = LALInferenceAngularDistance(*(REAL8 *)j_item->value, means[j]);
                        else if( j_item->vary==LALINFERENCE_PARAM_LINEAR )
                          jval = *(REAL8 *)j_item->value - means[j];
			
			for( k_item=item, k=0; k<=j; k_item=k_item->next ){	
                          REAL8 kval = 0.;
                          if(k_item->vary!=LALINFERENCE_PARAM_LINEAR && k_item->vary!=LALINFERENCE_PARAM_CIRCULAR) {
					continue;}
                                        
                                        if( k_item->vary==LALINFERENCE_PARAM_CIRCULAR )
                                          kval = LALInferenceAngularDistance(*(REAL8 *)k_item->value, means[k]);
                                        else if( k_item->vary==LALINFERENCE_PARAM_LINEAR )
                                          kval = *(REAL8 *)k_item->value - means[k];

					gsl_matrix_set(*cvm,j,k,gsl_matrix_get(*cvm,j,k) +
							   kval*jval);
					k++;
			}
			j++;
		}
	}

	/* Normalise */
	for(i=0;i<ND;i++) for(j=0;j<ND;j++) gsl_matrix_set(*cvm,i,j,gsl_matrix_get(*cvm,i,j)/((REAL8) Nlive));
	free(means);
        
	/* Fill in variances for circular parameters */
	/*for(item=Live[0]->head,j=0;item;item=item->next) {
		if(item->vary!=LALINFERENCE_PARAM_CIRCULAR && item->vary!=LALINFERENCE_PARAM_LINEAR) continue;
		if(item->vary==LALINFERENCE_PARAM_CIRCULAR) {
			for(k=0;k<j;k++) gsl_matrix_set(*cvm,j,k,0.0);
			gsl_matrix_set(*cvm,j,j,LALInferenceAngularVariance(Live,item->name,Nlive));
			for(k=j+1;k<ND;k++) gsl_matrix_set(*cvm,k,j,0.0);
		}
		j++;
	}*/
	
	/* the other half */
	for(i=0;i<ND;i++) 
          for(j=0;j<i;j++)
            gsl_matrix_set(*cvm,j,i,gsl_matrix_get(*cvm,i,j));
       
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

void LALInferenceNestedSamplingAlgorithm(LALInferenceRunState *runState)
{
	UINT4 iter=0,i,j,minpos;
	UINT4 Nlive=*(UINT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nlive");
	UINT4 Nruns=100;
	REAL8 *logZarray,*oldZarray,*Harray,*logwarray,*Wtarray,*logtarray,*logt2array;
	REAL8 TOLERANCE=0.1;
	REAL8 logZ,logZnew,logLmin,logLmax=-DBL_MAX,logLtmp,logw,H,logZnoise,dZ=0;//deltaZ - set but not used
	LALInferenceVariables *temp;
	FILE *fpout=NULL;
	gsl_matrix **cvm=calloc(1,sizeof(gsl_matrix *));
	REAL8 dblmax=-DBL_MAX;
	REAL8 zero=0.0;
	REAL8 *logLikelihoods=NULL;
	UINT4 verbose=0;
    REAL8 sloppyfrac;
	UINT4 displayprogress=0;
	LALInferenceVariableItem *param_ptr;
	LALInferenceVariables currentVars;
	memset(&currentVars,0,sizeof(currentVars));
	

	/* Default sample logging functions with and without XML */
#ifdef HAVE_LIBLALXML
  	char *outVOTable=NULL;
	if(!runState->logsample) runState->logsample=LALInferenceLogSampleToArray;
#else
	if(!runState->logsample) runState->logsample=LALInferenceLogSampleToFile;
#endif
	
        if ( !LALInferenceCheckVariable(runState->algorithmParams, "logZnoise" ) ){
          /*if (runState->data->modelDomain == LALINFERENCE_DOMAIN_FREQUENCY )*/
            logZnoise=LALInferenceNullLogLikelihood(runState->data);
	
          LALInferenceAddVariable(runState->algorithmParams,"logZnoise",&logZnoise,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
        }
        else{
          logZnoise = 
            *(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"logZnoise");
        }
        
        logLikelihoods=(REAL8 *)(*(REAL8Vector **)LALInferenceGetVariable(runState->algorithmParams,"logLikelihoods"))->data;

	verbose=LALInferenceCheckVariable(runState->algorithmParams,"verbose");
    displayprogress=verbose;
	
	/* Operate on parallel runs if requested */
	if(LALInferenceCheckVariable(runState->algorithmParams,"Nruns"))
		Nruns = *(UINT4 *) LALInferenceGetVariable(runState->algorithmParams,"Nruns");

	if(LALInferenceCheckVariable(runState->algorithmParams,"tolerance"))
		TOLERANCE = *(REAL8 *) LALInferenceGetVariable(runState->algorithmParams,"tolerance");

	/* Check that necessary parameters are created */
	if(!LALInferenceCheckVariable(runState->algorithmParams,"logLmin"))
		LALInferenceAddVariable(runState->algorithmParams,"logLmin",&dblmax,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);

	if(!LALInferenceCheckVariable(runState->algorithmParams,"accept_rate"))
		LALInferenceAddVariable(runState->algorithmParams,"accept_rate",&zero,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
    if(!LALInferenceCheckVariable(runState->algorithmParams,"sub_accept_rate"))
        LALInferenceAddVariable(runState->algorithmParams,"sub_accept_rate",&zero,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);

	/* Set up the proposal scale factor, for use in the multi-student jump step */
	REAL8 propScale = 0.1;
	LALInferenceAddVariable(runState->proposalArgs,"proposal_scale",&propScale,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);

	/* Open output file */
        ProcessParamsTable *ppt = NULL; 
	ppt=LALInferenceGetProcParamVal(runState->commandLine,"--outfile");
	if(!ppt){
		fprintf(stderr,"Must specify --outfile <filename.dat>\n");
		exit(1);
	}
	char *outfile=ppt->value;
       
    if(LALInferenceGetProcParamVal(runState->commandLine,"--progress"))
        displayprogress=1;

#ifdef HAVE_LIBLALXML
	ppt=LALInferenceGetProcParamVal(runState->commandLine,"--outxml");
	if(!ppt){
		ppt=LALInferenceGetProcParamVal(runState->commandLine,"--outXML");
	}
	if(!ppt){
		fprintf(stderr,"Can specify --outXML <filename.dat> for VOTable output\n");
	}
	else{
		outVOTable=ppt->value;
	}
#endif	
	fpout=fopen(outfile,"w");

	if(fpout==NULL) fprintf(stderr,"Unable to open output file %s!\n",outfile);
	else
	  LALInferenceAddVariable(runState->algorithmParams,"outfile",&fpout,LALINFERENCE_void_ptr_t,LALINFERENCE_PARAM_FIXED);
	
	//fprintf(fpout,"chirpmass\tdistance\tLAL_APPROXIMANT\tLAL_PNORDER\tlogmc\tmassratio\ttime\tphase\tlogdistance\trightascension\tdeclination\tpolarisation\tinclination\ta_spin1\ta_spin2\ttheta_spin1\ttheta_spin2\tphi_spin1\tphi_spin2\t logL\n");	
	/* Set up arrays for parallel runs */
	minpos=0;
	        /*   for(param_ptr=runState->livePoints[minpos]->head;param_ptr;param_ptr=param_ptr->next)
   {
        fprintf(fpout,"%s\t",param_ptr->name);
    }	
	fprintf(fpout,"logL\n");*/
	logZarray = calloc(Nruns,sizeof(REAL8));
	oldZarray = calloc(Nruns,sizeof(REAL8));
	Harray = calloc(Nruns,sizeof(REAL8));
	logwarray = calloc(Nruns,sizeof(REAL8));
	logtarray=calloc(Nruns,sizeof(REAL8));
    logt2array=calloc(Nruns,sizeof(REAL8));
	Wtarray = calloc(Nruns,sizeof(REAL8));
	if(logZarray==NULL || Harray==NULL || oldZarray==NULL || logwarray==NULL || Wtarray==NULL)
		{fprintf(stderr,"Unable to allocate RAM\n"); exit(-1);}

	logw=log(1.0-exp(-1.0/Nlive));
	for(i=0;i<Nruns;i++)  {logwarray[i]=logw; logZarray[i]=-DBL_MAX; oldZarray[i]=-DBL_MAX; Harray[i]=0.0;logtarray[i]=-1.0/Nlive; logt2array[i]=-1.0/Nlive; }
	i=0;
	/* Find maximum likelihood */
	for(i=0;i<Nlive;i++)
	{
		logLtmp=logLikelihoods[i];
		logLmax=logLtmp>logLmax? logLtmp : logLmax;
	}
	/* Add the covariance matrix for proposal distribution */
	LALInferenceNScalcCVM(cvm,runState->livePoints,Nlive);
	runState->differentialPoints=runState->livePoints;
	runState->differentialPointsLength=(size_t) Nlive;
	/* Set up eigenvectors and eigenvalues. */
	UINT4 N=(*cvm)->size1;
	gsl_matrix *covCopy = gsl_matrix_alloc(N,N);
	gsl_matrix *eVectors = gsl_matrix_alloc(N,N);
	gsl_vector *eValues = gsl_vector_alloc(N);
	REAL8Vector *eigenValues = XLALCreateREAL8Vector(N);
	gsl_eigen_symmv_workspace *ws = gsl_eigen_symmv_alloc(N);
	int gsl_status;
	gsl_matrix_memcpy(covCopy, *cvm);
	
	if ((gsl_status = gsl_eigen_symmv(covCopy, eValues, eVectors, ws)) != GSL_SUCCESS) {
	  XLALPrintError("Error in gsl_eigen_symmv (in %s, line %d): %d: %s\n", __FILE__, __LINE__, gsl_status, gsl_strerror(gsl_status));
	  XLAL_ERROR_VOID(XLAL_EFAILED);
	}
	
	for (i = 0; i < N; i++) {
	  eigenValues->data[i] = gsl_vector_get(eValues,i);
	}
	
	LALInferenceAddVariable(runState->proposalArgs, "covarianceEigenvectors", &eVectors, LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_FIXED);
	LALInferenceAddVariable(runState->proposalArgs, "covarianceEigenvalues", &eigenValues, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED);
	
	LALInferenceAddVariable(runState->proposalArgs,"covarianceMatrix",cvm,LALINFERENCE_gslMatrix_t,LALINFERENCE_PARAM_OUTPUT);
      
        /* set up k-D tree if required and not already set */
        if ( ( LALInferenceGetProcParamVal(runState->commandLine,"--kDTree") ||
             LALInferenceGetProcParamVal(runState->commandLine,"--kdtree")) &&
         !LALInferenceCheckVariable( runState->proposalArgs, "kDTree" ) ){
          LALInferenceSetupkDTreeNSLivePoints( runState );
        }
        
	/* Sprinkle points */
	LALInferenceSetVariable(runState->algorithmParams,"logLmin",&dblmax);
	for(i=0;i<Nlive;i++) {
	  runState->currentParams=runState->livePoints[i];
	  runState->evolve(runState);
	  logLikelihoods[i]=runState->likelihood(runState->livePoints[i],runState->data,runState->template);
          LALInferenceAddVariable(runState->livePoints[i],"logw",&logw,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
	  if(XLALPrintProgressBar((double)i/(double)Nlive)) fprintf(stderr,"\n");
	}
	/* re-calculate the k-D tree from the new points if required */
	if ( LALInferenceCheckVariable( runState->proposalArgs, "kDTree" ) ) 
          LALInferenceSetupkDTreeNSLivePoints( runState );

	runState->currentParams=&currentVars;
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
			Wtarray[j]=logwarray[j]+logLikelihoods[minpos]+logadd(0,logt2array[j]*logtarray[j])-log(2.0);
			logZarray[j]=logadd(logZarray[j],Wtarray[j]);
			Harray[j]= exp(Wtarray[j]-logZarray[j])*logLikelihoods[minpos]
			+ exp(oldZarray[j]-logZarray[j])*(Harray[j]+oldZarray[j])-logZarray[j];
		}
		logZnew=mean(logZarray,Nruns);
		//deltaZ=logZnew-logZ; - set but not used
		H=mean(Harray,Nruns);
		logZ=logZnew;
		for(j=0;j<Nruns;j++) oldZarray[j]=logZarray[j];
               
		if(runState->logsample) runState->logsample(runState,runState->livePoints[minpos]);
		
		UINT4 itercounter=0;
		
		/* Generate a new live point */
		do{ /* This loop is here in case it is necessary to find a different sample */
			/* Clone an old live point and evolve it */
                        while((j=gsl_rng_uniform_int(runState->GSLrandom,Nlive))==minpos){};
			LALInferenceCopyVariables(runState->livePoints[j],runState->currentParams);
			runState->currentLikelihood = logLikelihoods[j];
			LALInferenceSetVariable(runState->algorithmParams,"logLmin",(void *)&logLmin);
                        runState->evolve(runState);
                        itercounter++;
		}while( runState->currentLikelihood<=logLmin ||  *(REAL8*)LALInferenceGetVariable(runState->algorithmParams,"accept_rate")==0.0);

                LALInferenceCopyVariables(runState->currentParams,runState->livePoints[minpos]);
                logLikelihoods[minpos]=runState->currentLikelihood;
     
    //            for(param_ptr=runState->livePoints[minpos]->head;param_ptr;param_ptr=param_ptr->next)
   // {
      //  fprintf(fpout,"%s\t",param_ptr->name);
    //}        
	//fprintf(fpout,"chirpmass\tdistance\tLAL_APPROXIMANT\tLAL_PNORDER\tlogmc\tmassratio\ttime\tphase\tlogdistance\trightascension\tdeclination\tpolarisation\tinclination\ta_spin1\ta_spin2\ttheta_spin1\ttheta_spin2\tphi_spin1\tphi_spin2\t logL\n"); 
		if (runState->currentLikelihood>logLmax)
			logLmax=runState->currentLikelihood;
		for(j=0;j<Nruns;j++) {
          REAL8 *tmp=logtarray;
          logtarray=logt2array;
          logt2array=tmp;
		  logtarray[j]=LALInferenceNSSample_logt(Nlive,runState->GSLrandom);
		  logwarray[j]+=logtarray[j];
		}
		logw=mean(logwarray,Nruns);
		LALInferenceAddVariable(runState->livePoints[minpos],"logw",&logw,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
		dZ=logadd(logZ,logLmax-((double) iter)/((double)Nlive))-logZ;
		sloppyfrac=*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"sloppyfraction");
		if(displayprogress) fprintf(stderr,"%i: (%2.1lf%%) accpt: %1.3f Nmcmc: %i sub_accpt: %1.3f slpy: %2.1f%% H: %3.3lf nats (%3.3lf b) logL:%lf ->%lf logZ: %lf dZ: %lf Zratio: %lf db\n",\
		  iter,\
		  100.0*((REAL8)iter)/(((REAL8) Nlive)*H),\
		  *(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"accept_rate")/(REAL8)itercounter,\
		  *(INT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nmcmc"),\
		  *(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"sub_accept_rate"),\
		  100.0*sloppyfrac,\
		  H,\
		  H/LAL_LN2,\
		  logLmin,\
		  runState->currentLikelihood,\
		  logZ,\
		  dZ,\
		  10.0*LAL_LOG10E*( logZ-*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"logZnoise")));
		/* Flush output file */
		if(fpout && !(iter%100)) fflush(fpout);
		iter++;
		/* Update the proposal */
		if(!(iter%(Nlive/4))) {
            /* Update the covariance matrix */
            if ( LALInferenceCheckVariable( runState->proposalArgs,"covarianceMatrix" ) ){
		        LALInferenceNScalcCVM(cvm,runState->livePoints,Nlive);
		        gsl_matrix_memcpy(covCopy, *cvm);
		        if ((gsl_status = gsl_eigen_symmv(covCopy, eValues, eVectors, ws)) != GSL_SUCCESS) {
		            XLALPrintError("Error in gsl_eigen_symmv (in %s, line %d): %d: %s\n", __FILE__, __LINE__, gsl_status, gsl_strerror(gsl_status));
		            XLAL_ERROR_VOID(XLAL_EFAILED);
		        }
		        for (i = 0; i < N; i++) {
		            eigenValues->data[i] = gsl_vector_get(eValues,i);
		        }
		        LALInferenceAddVariable(runState->proposalArgs, "covarianceEigenvectors", &eVectors, LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_FIXED);
		        LALInferenceAddVariable(runState->proposalArgs, "covarianceEigenvalues", &eigenValues, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED);
		        LALInferenceSetVariable(runState->proposalArgs,"covarianceMatrix",(void *)cvm);
            }
                  
            /* update k-d tree */
            if ( LALInferenceCheckVariable( runState->proposalArgs,"kDTree" ) )
                LALInferenceSetupkDTreeNSLivePoints( runState ); 
		
            /* Measure Autocorrelations if asked */
		    if(LALInferenceGetProcParamVal(runState->commandLine,"--auto-chain-length")){
                LALInferenceAddVariable(runState->algorithmParams,"current_iteration",&iter,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_OUTPUT);
                //INT4 meanmax=0.;
                //int Ntries=5;
                /* Calculate ACF of Prior MCMC sampler */
                //INT4 sloppyratio;
                //for (int try=0;try<Ntries;try++){
                //    sloppyratio=*(INT4 *)LALInferenceGetVariable(runState->algorithmParams,"sloppyratio");
                //    LALInferenceVariables *acls=LALInferenceComputeAutoCorrelation(runState, sloppyratio*10, SamplePriorDiscardAcceptance);
                //    INT4 max=1;
                //    for(LALInferenceVariableItem *this=acls->head;this;this=this->next) { if(*(REAL8 *)this->value>max) max=(INT4) *(REAL8 *)this->value;}
                //    LALInferenceDestroyVariables(acls);
                //    free(acls);
                //    meanmax+=max;
                //}
                //meanmax/=Ntries;
                //if(sloppyratio>meanmax)
                 //   LALInferenceSetVariable(runState->algorithmParams,"sloppyratio",&meanmax);
              
                /* Measure ACF of Actual evolution function */
		
		/* Reset the sloppy fraction as we have updated our proposal */
		//LALInferenceSetVariable(runState->algorithmParams,"sloppylogit",&zero);
			    INT4 MAX_MCMC=20000; /* Maximum chain length, set to be higher than expected from a reasonable run */
                INT4 max=4 * *(INT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nmcmc"); /* We will use this to go out 4x last ACL */
                if(max>MAX_MCMC) max=MAX_MCMC;
                LALInferenceVariables *acls=LALInferenceComputeAutoCorrelation(runState, max*4, runState->evolve) ;
                max=10;
                for(LALInferenceVariableItem *this=acls->head;this;this=this->next) { if(*(REAL8 *)this->value>max) max=(INT4) *(REAL8 *)this->value;}
                LALInferenceDestroyVariables(acls);
                free(acls);
                max*=2;
                LALInferenceSetVariable(runState->algorithmParams,"Nmcmc",&max);
		    }
	      }
		
	}
	while( iter <= Nlive ||  dZ> TOLERANCE ); 

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
			logLikelihoods[minpos]=logLikelihoods[i];
			logLikelihoods[i]=logLmin;
		}
	/* final corrections */
        for(i=0;i<Nlive;i++){
		logZ=logadd(logZ,logLikelihoods[i]+logw);
		for(j=0;j<Nruns;j++){
			//logwarray[j]+=LALInferenceNSSample_logt(Nlive,runState->GSLrandom);
			logZarray[j]=logadd(logZarray[j],logLikelihoods[i]+logwarray[j]-log(Nlive));
		}

		if(runState->logsample) runState->logsample(runState,runState->livePoints[i]);

	}

	/* Write out the evidence */
	fclose(fpout);
	char bayesfile[FILENAME_MAX];
	sprintf(bayesfile,"%s_B.txt",outfile);
	fpout=fopen(bayesfile,"w");
	fprintf(fpout,"%lf %lf %lf %lf\n",logZ-logZnoise,logZ,logZnoise,logLmax);
	fclose(fpout);
	double logB=logZ-logZnoise;
	/* Pass output back through algorithmparams */
	LALInferenceAddVariable(runState->algorithmParams,"logZ",(void *)&logZ,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
	LALInferenceAddVariable(runState->algorithmParams,"logB",(void *)&logB,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
	LALInferenceAddVariable(runState->algorithmParams,"logLmax",(void *)&logLmax,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);

#ifdef HAVE_LIBLALXML	
	/* Write out the XML if requested */
    LALInferenceVariables *output_array=NULL;
    UINT4 N_output_array=0;
	if(LALInferenceCheckVariable(runState->algorithmParams,"outputarray")
	  &&LALInferenceCheckVariable(runState->algorithmParams,"N_outputarray") )
	{
	  output_array=*(LALInferenceVariables **)LALInferenceGetVariable(runState->algorithmParams,"outputarray");
	  N_output_array=*(UINT4 *)LALInferenceGetVariable(runState->algorithmParams,"N_outputarray");
	}
	if(output_array && outVOTable && N_output_array>0){
		xmlNodePtr votable=XLALInferenceVariablesArray2VOTTable(output_array, N_output_array, "Nested Samples");
		xmlNewProp(votable, CAST_CONST_XMLCHAR("utype"), CAST_CONST_XMLCHAR("lalinference:results:nestedsamples"));
		
		xmlNodePtr stateResource=XLALInferenceStateVariables2VOTResource(runState, "Run State Configuration");
		
		xmlNodePtr nestResource=XLALCreateVOTResourceNode("lalinference:results","Nested sampling run",votable);
		
		if(stateResource)
			xmlAddChild(nestResource,stateResource);
		

		char *xmlString = XLALCreateVOTStringFromTree ( nestResource );
		
		/* Write to disk */
		fpout=fopen(outVOTable,"w");
		fprintf(fpout,"%s",xmlString);
		fclose(fpout);
		
		
	}
	if(output_array) free(output_array);

#endif
	/* Write out names of parameters */
	FILE *lout=NULL;
	char param_list[FILENAME_MAX];
	sprintf(param_list,"%s_params.txt",outfile);
	lout=fopen(param_list,"w");
	minpos=0;
	LALInferenceSortVariablesByName(runState->livePoints[0]);
	for(param_ptr=runState->livePoints[0]->head;param_ptr;param_ptr=param_ptr->next)
	{
	  fprintf(lout,"%s\t",param_ptr->name);
	}
	fclose(lout);
	free(logtarray); free(logwarray); free(logZarray);
}

/* Calculate the autocorrelation function of the sampler (runState->evolve) for each parameter
 * Evolves the sample starting with the value passed in temp, with a maximum of max_iterations steps.
 Return the ACL for each parameter as a LALInferenceVariables */
LALInferenceVariables *LALInferenceComputeAutoCorrelation(LALInferenceRunState *runState, UINT4 max_iterations, LALInferenceEvolveOneStepFunction *evolve)
{
  ProcessParamsTable *ppt=NULL;
  char chainfilename[128]="";
  char acf_file_name[128]="";
  FILE *chainfile=NULL;
  FILE *acffile=NULL;
  INT4 global_iter;
  UINT4 i,j;
  INT4 nPar=0; // = LALInferenceGetVariableDimensionNonFixed(runState->currentParams);
  REAL8 **data_array=NULL;
  REAL8 **acf_array=NULL;
  LALInferenceVariableItem *this;
  REAL8 tolerance=0.01;
  INT4 thinning=10;
  max_iterations/=thinning;
  for(this=runState->currentParams->head;this;this=this->next) if(this->vary!=LALINFERENCE_PARAM_FIXED && this->vary!=LALINFERENCE_PARAM_OUTPUT && this->type==LALINFERENCE_REAL8_t) nPar++;

  REAL8 ACF,ACL,max=0;
  LALInferenceVariables *acls=calloc(1,sizeof(LALInferenceVariables));

  global_iter=*(INT4 *)LALInferenceGetVariable(runState->algorithmParams,"current_iteration");
  /* Back up the algorithm state and replace with a clean version for logSampletoarray */
  LALInferenceVariables myAlgParams,*oldAlgParams=runState->algorithmParams;
  LALInferenceVariables myCurrentParams,*oldCurrentParams=runState->currentParams;
  memset(&myAlgParams,0,sizeof(LALInferenceVariables));
  memset(&myCurrentParams,0,sizeof(LALInferenceVariables));
  LALInferenceCopyVariables(runState->currentParams,&myCurrentParams);
  LALInferenceCopyVariables(oldAlgParams,&myAlgParams);
  LALInferenceRemoveVariable(&myAlgParams,"outputarray");
  LALInferenceRemoveVariable(&myAlgParams,"N_outputarray");
  LALInferenceRemoveVariable(&myAlgParams,"outfile");
  LALInferenceSetVariable(&myAlgParams,"Nmcmc",&thinning);

  REAL8Vector *projection=NULL;

  gsl_matrix *covarianceMatrix=NULL;
  runState->algorithmParams=&myAlgParams;
  runState->currentParams=&myCurrentParams;
  /* We can record write the MCMC chain to a file too */
  ppt=LALInferenceGetProcParamVal(runState->commandLine,"--acf-chainfile");
  if(ppt){
    sprintf(chainfilename,"%s.%i",ppt->value,global_iter);
    chainfile=fopen(chainfilename,"w");
    LALInferenceAddVariable(&myAlgParams,"outfile",&chainfile,LALINFERENCE_void_ptr_t,LALINFERENCE_PARAM_FIXED);
  }
  ppt=LALInferenceGetProcParamVal(runState->commandLine,"--acf-file");
  if(ppt){
    sprintf(acf_file_name,"%s.%i",ppt->value,global_iter);
    acffile=fopen(acf_file_name,"w");
  }
  LALInferenceVariables **livePoints=runState->livePoints;
  UINT4 Nlive = *(INT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nlive");
  UINT4 BAILOUT=100; /* this should be the same as the bailout in the sampler */
  REAL8 accept=0.0;
  do{ /* Pick a random sample that isn't trapped in some corner*/
	UINT4 idx=gsl_rng_uniform_int(runState->GSLrandom,Nlive);
	runState->currentParams=livePoints[idx];
    i=0;
  	do {
        evolve(runState);
        i++;
        accept=*(REAL8*)LALInferenceGetVariable(runState->algorithmParams,"accept_rate");
        }
        while(accept==0.0 && i<BAILOUT);
  }
  while(0.==accept);
	/* log the first sample*/ 
  LALInferenceLogSampleToArray(runState,runState->currentParams);
  /* Evolve the initial sample (i starts at 1)*/
  for(i=1;i<max_iterations;i++)
  {
   evolve(runState);
   LALInferenceLogSampleToArray(runState,runState->currentParams);
  }
  
  /* Get the location of the sample array */
  LALInferenceVariables *variables_array=*(LALInferenceVariables **)LALInferenceGetVariable(runState->algorithmParams,"outputarray");
  LALInferenceVariables **pointer_array=calloc(max_iterations,sizeof(LALInferenceVariables *));
  for(i=0;i<max_iterations;i++) pointer_array[i]=&(variables_array[i]);
  
   	/* Calculate the eigenvectors of the correlation matrix*/
	LALInferenceNScalcCVM(&covarianceMatrix, pointer_array,max_iterations);

   free(pointer_array);
	/* Set up eigenvectors and eigenvalues. */
	UINT4 N=covarianceMatrix->size1;
	gsl_matrix *covCopy = gsl_matrix_alloc(N,N);
	gsl_matrix *eVectors = gsl_matrix_alloc(N,N);
	gsl_vector *eValues = gsl_vector_alloc(N);

	gsl_eigen_symmv_workspace *ws = gsl_eigen_symmv_alloc(N);
	int gsl_status;
	gsl_matrix_memcpy(covCopy, covarianceMatrix);
	
	if ((gsl_status = gsl_eigen_symmv(covCopy, eValues, eVectors, ws)) != GSL_SUCCESS) {
	  XLALPrintError("Error in gsl_eigen_symmv (in %s, line %d): %d: %s\n", __FILE__, __LINE__, gsl_status, gsl_strerror(gsl_status));
	  return ((LALInferenceVariables *)NULL);
	}
	
	gsl_matrix_free(covarianceMatrix);
	gsl_matrix_free(covCopy);
	gsl_vector_free(eValues);

  /* Convert to a 2D array for ACF calculation */
  data_array=calloc(nPar,sizeof(REAL8 *));
  acf_array=calloc(nPar,sizeof(REAL8 *));
  for (i=0;i<(UINT4)nPar;i++){
    data_array[i]=calloc(max_iterations,sizeof(REAL8));
    acf_array[i]=calloc(max_iterations/2,sizeof(REAL8));
  }

  for (i=0;i<max_iterations;i++){
    j=0;
    LALInferenceProjectSampleOntoEigenvectors( &variables_array[i], eVectors, &projection );
    data_array[j][i]=projection->data[j];
    //for(this=variables_array[i].head;this;this=this->next)
    //{
    //  switch(this->vary){
	//case LALINFERENCE_PARAM_CIRCULAR:
	//case LALINFERENCE_PARAM_LINEAR:
	//{
	 // if(this->type!=LALINFERENCE_REAL8_t) continue;
	 // else {
	 //   data_array[j][i]=*(REAL8 *)this->value;
	    j++;
	 // }
	//}
	//default:
	 // continue;
    //  }
    //}
    LALInferenceDestroyVariables(&variables_array[i]);
  }
  	gsl_matrix_free(eVectors);
	XLALDestroyREAL8Vector(projection);
  free(variables_array);
  this=myCurrentParams.head;
  for(i=0;i<(UINT4)nPar;i++){
   /* Subtract the mean */
   REAL8 this_mean = gsl_stats_mean(&data_array[i][0], 1, max_iterations);
   for(j=0;j<max_iterations;j++) data_array[i][j]-=this_mean;
   ACL=1;
   int startflag=1;
   ACF=1.;
   /* Use GSL to compute the ACF */
   for(UINT4 lag=0;ACF>=tolerance&&lag<max_iterations/2;lag++){
      ACF=(REAL8) gsl_stats_correlation(&data_array[i][0], 1, &data_array[i][lag], 1, max_iterations-lag);
      acf_array[i][lag]=ACF;
      ACL+=2.0*ACF;
      if((ACF<tolerance && startflag) || lag==max_iterations/2-1){
	    startflag=0;
        ACL*=(REAL8)thinning;
	    if(ACL>max) max=ACL;
	    LALInferenceAddVariable(acls,this->name,&ACL,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
        break;
      }
   }   
   do{this=this->next;}while(this && (this->vary==LALINFERENCE_PARAM_FIXED || this->vary==LALINFERENCE_PARAM_OUTPUT || this->type!=LALINFERENCE_REAL8_t));
  }
  if(acffile){
  /* Write out the ACF */
  for(i=0;i<max_iterations/2;i++){
    for(j=0;j<(UINT4)nPar;j++) fprintf(acffile,"%f ",acf_array[j][i]);
    fprintf(acffile,"\n");
  }
  }
/*  
  FILE *aclfile=fopen("acl.dat","a");
  FILE *aclfile_header=fopen("acl_params.txt","w");
  fprintf(aclfile,"%i ",global_iter);
  for(this=acls->head;this;this=this->next) {
    fprintf(aclfile,"%lf ",*(REAL8 *)this->value);
    fprintf(aclfile_header,"%s ",this->name);
  }
  fprintf(aclfile,"\n");
  fprintf(aclfile_header,"\n");
  fclose(aclfile_header);
  fclose(aclfile);
*/  
  /* Clean up */
  for(i=0;i<(UINT4)nPar;i++) {free(data_array[i]); free(acf_array[i]);}
  free(data_array); free(acf_array);
  LALInferenceDestroyVariables(&myAlgParams);
  LALInferenceDestroyVariables(&myCurrentParams);
  runState->currentParams=oldCurrentParams;
  runState->algorithmParams=oldAlgParams;
  if(chainfile) fclose(chainfile);
  if(acffile) fclose(acffile);
  return(acls);
}

//void SamplePriorDiscardAcceptance(LALInferenceRunState *runState)
//{
//    LALInferenceMCMCSamplePrior(runState);
//    return;
//}

/* Perform one MCMC iteration on runState->currentParams. Return 1 if accepted or 0 if not */
UINT4 LALInferenceMCMCSamplePrior(LALInferenceRunState *runState)
{
    LALInferenceVariables tempParams;
    REAL8 logProposalRatio=0.0;
    memset(&tempParams,0,sizeof(tempParams));
    LALInferenceVariables *oldParams=&tempParams;
    UINT4 accepted=0;

    REAL8 logPriorOld=*(REAL8 *)LALInferenceGetVariable(runState->currentParams,"logPrior");
    LALInferenceCopyVariables(runState->currentParams,oldParams);
    runState->proposal(runState,runState->currentParams);
    REAL8 logPriorNew=runState->prior(runState,runState->currentParams);
    if(LALInferenceCheckVariable(runState->proposalArgs,"logProposalRatio"))
       logProposalRatio=*(REAL8 *)LALInferenceGetVariable(runState->proposalArgs,"logProposalRatio");
    if(logPriorNew==-DBL_MAX || isnan(logPriorNew) || log(gsl_rng_uniform(runState->GSLrandom)) > (logPriorNew-logPriorOld) + logProposalRatio) 
    {
        LALInferenceCopyVariables(oldParams,runState->currentParams);
    } 
    else {
        accepted=1;
        LALInferenceSetVariable(runState->currentParams,"logPrior",&logPriorNew);
    }
    LALInferenceDestroyVariables(oldParams);

    return(accepted);
}

/* Sample the prior N times, returns number of acceptances */
UINT4 LALInferenceMCMCSamplePriorNTimes(LALInferenceRunState *runState, UINT4 N)
{
    UINT4 i=0;
    UINT4 Naccepted=0;
    for(i=0;i<N;i++) Naccepted+=LALInferenceMCMCSamplePrior(runState);
    return(Naccepted);
}

void LALInferenceProjectSampleOntoEigenvectors(LALInferenceVariables *params, gsl_matrix *eigenvectors, REAL8Vector **projection)
{
	  LALInferenceVariableItem *proposeIterator = params->head;
	  UINT4 j=0,i=0;
	  UINT4 N=eigenvectors->size1;

     if(!*projection) *projection=XLALCreateREAL8Vector(N);
	  if((*projection)->length==0) *projection=XLALCreateREAL8Vector(N);


	if (proposeIterator == NULL) {
    fprintf(stderr, "Bad proposed params in %s, line %d\n",
            __FILE__, __LINE__);
    exit(1);
  }
  for(i=0;i<N;i++){
  	j=0;
  	proposeIterator = params->head;
  	(*projection)->data[i]=0.0;
		do { 
		    if (proposeIterator->vary != LALINFERENCE_PARAM_FIXED && proposeIterator->vary != LALINFERENCE_PARAM_OUTPUT) {
     			 	(*projection)->data[i]+= *((REAL8 *)proposeIterator->value) * gsl_matrix_get(eigenvectors, j, i);
      			j++;
    		 }
 		}while ((proposeIterator = proposeIterator->next) != NULL && j < N);
  } 
	
}


/* Sample the limited prior distribution using the MCMC method as usual, but
   only check the likelihood bound x fraction of the time. Always returns a fulled checked sample.
   x=LALInferenceGetVariable(runState->algorithmParams,"sloppyfraction")
   */

void LALInferenceNestedSamplingSloppySample(LALInferenceRunState *runState)
{
    LALInferenceVariables oldParams;
    LALInferenceIFOData *data=runState->data;
    REAL8 tmp;
    REAL8 Target=0.3;
    char tmpName[32];
    REAL8 logLold=*(REAL8 *)LALInferenceGetVariable(runState->currentParams,"logL");
    memset(&oldParams,0,sizeof(oldParams));
    LALInferenceCopyVariables(runState->currentParams,&oldParams);
    REAL8 logLmin=*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"logLmin");
    UINT4 Nmcmc=*(UINT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nmcmc");
    REAL8 sloppyfraction=0.;
    REAL8 maxsloppyfraction=((REAL8)Nmcmc-1)/(REAL8)Nmcmc ;
    REAL8 minsloppyfraction=0.;
    if(Nmcmc==1) maxsloppyfraction=minsloppyfraction=0.0;
    if (LALInferenceCheckVariable(runState->algorithmParams,"sloppyfraction"))
      sloppyfraction=*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"sloppyfraction");
    UINT4 mcmc_iter=0,Naccepted=0,sub_accepted=0;
    //INT4 Nlive=*(INT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nlive");
    UINT4 sloppynumber=(UINT4) (sloppyfraction*(REAL8)Nmcmc);
    UINT4 testnumber=Nmcmc-sloppynumber;
    /* +1 for the last iteration which we do check */
    UINT4 subchain_length=(sloppynumber/testnumber) +1;
    //sloppynumber=sloppynumber<1?1:sloppynumber;
    REAL8 logLnew;
    UINT4 sub_iter=0;
    UINT4 tries=0;
    UINT4 BAILOUT=100*testnumber; /* If no acceptance after 100 tries, will exit and the sampler will try a different starting point */
    //REAL8 *logLikelihoods=(REAL8 *)(*(REAL8Vector **)LALInferenceGetVariable(runState->algorithmParams,"logLikelihoods"))->data;
    do{
        /* Draw an independent sample from the prior */
        sub_accepted+=LALInferenceMCMCSamplePriorNTimes(runState,subchain_length);
	    if(sub_accepted==0.) {
            //    INT4 j=gsl_rng_uniform_int(runState->GSLrandom,Nlive);
            //    LALInferenceCopyVariables(runState->livePoints[j],runState->currentParams);
            //    runState->currentLikelihood = logLikelihoods[j];
	        tries++;
	        sub_iter+=subchain_length;
	        mcmc_iter++;
            LALInferenceCopyVariables(&oldParams,runState->currentParams);
            runState->currentLikelihood=logLold;
	        continue;
        }
        tries=0;
        mcmc_iter++;
    	sub_iter+=subchain_length;
        logLnew=runState->likelihood(runState->currentParams,runState->data,runState->template);
        if(logLnew>logLmin || logLmin==-DBL_MAX) /* Accept */
        {
            Naccepted++;
            /* Update information to pass back out */
            LALInferenceAddVariable(runState->currentParams,"logL",(void *)&logLnew,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
            if(LALInferenceCheckVariable(runState->algorithmParams,"logZnoise")){
               tmp=logLnew-*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"logZnoise");
               LALInferenceAddVariable(runState->currentParams,"deltalogL",(void *)&tmp,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
            }
            while(data)
            {
               tmp=data->loglikelihood - data->nullloglikelihood;
               sprintf(tmpName,"deltalogl%s",data->name);
               LALInferenceAddVariable(runState->currentParams,tmpName,&tmp,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
               data=data->next;
            }
            LALInferenceCopyVariables(runState->currentParams,&oldParams);
            logLold=logLnew;
            runState->currentLikelihood=logLnew;
        }
        else /* reject */
        {
            LALInferenceCopyVariables(&oldParams,runState->currentParams);
            runState->currentLikelihood=logLold;
        }
    }while((mcmc_iter<testnumber||logLnew<=logLmin||Naccepted==0)&&(mcmc_iter<BAILOUT));
    REAL8 sub_accept_rate=(REAL8)sub_accepted/(REAL8)sub_iter;
    REAL8 accept_rate=(REAL8)Naccepted/(REAL8)testnumber;
    LALInferenceSetVariable(runState->algorithmParams,"accept_rate",&accept_rate);
    LALInferenceSetVariable(runState->algorithmParams,"sub_accept_rate",&sub_accept_rate);
    if(logLmin!=-DBL_MAX){
        if((REAL8)accept_rate>Target) { sloppyfraction+=1.0/(REAL8)Nmcmc;}
        else { sloppyfraction-=1.0/(REAL8)Nmcmc;}
        if(sloppyfraction>maxsloppyfraction) sloppyfraction=maxsloppyfraction;
	if(sloppyfraction<minsloppyfraction) sloppyfraction=minsloppyfraction;
        //if(sloppylogit<1) sloppy=1;
	
	LALInferenceSetVariable(runState->algorithmParams,"sloppyfraction",&sloppyfraction);
        //LALInferenceSetVariable(runState->algorithmParams,"Nmcmc",&Nmcmc);
    }
    LALInferenceDestroyVariables(&oldParams);
}


/* Evolve nested sampling algorithm by one step, i.e.
 evolve runState->currentParams to a new point with higher
 likelihood than currentLikelihood. Uses the MCMC method.
 */
void LALInferenceNestedSamplingOneStep(LALInferenceRunState *runState)
{
	LALInferenceVariables *newParams=NULL,*lastAcceptedWithLikelihood=NULL;
	LALInferenceIFOData *data=runState->data;
	char tmpName[32];
	UINT4 mcmc_iter=0,Naccepted=0;
	UINT4 Nmcmc=*(UINT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nmcmc");
	REAL8 logLmin=*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"logLmin");
        REAL8 logPriorOld,logPriorNew,logLnew=DBL_MAX,logPriorlastAccepted;
	REAL8 logProposalRatio,tmp=0.0;
	newParams=calloc(1,sizeof(LALInferenceVariables));
	lastAcceptedWithLikelihood=calloc(1,sizeof(LALInferenceVariables));
	INT4 thin_factor=1;
	/* Evolve the sample until it is accepted */
	logPriorOld=runState->prior(runState,runState->currentParams);
	logPriorlastAccepted=logPriorOld;
	LALInferenceCopyVariables(runState->currentParams,lastAcceptedWithLikelihood);
	do{
		mcmc_iter++;
		/* Make a copy of the parameters passed through currentParams */
                LALInferenceCopyVariables(runState->currentParams,newParams);
                runState->proposal(runState,newParams);
		logPriorNew=runState->prior(runState,newParams);
		LALInferenceSetVariable(newParams,"logPrior",&logPriorNew);
		if(LALInferenceCheckVariable(runState->proposalArgs,"logProposalRatio"))
		  logProposalRatio=*(REAL8 *)LALInferenceGetVariable(runState->proposalArgs,"logProposalRatio");
		else logProposalRatio=0.0;
                
		/* If rejected, continue to next iteration */
		if(logPriorNew==-DBL_MAX || isnan(logPriorNew)) continue;
                if(log(gsl_rng_uniform(runState->GSLrandom)) > (logPriorNew-logPriorOld) + logProposalRatio) {
			        continue;
                }
		if(mcmc_iter%thin_factor==(UINT4)(thin_factor-1)){
		  /* Otherwise, check that logL is OK */
		  if(logLmin!=-DBL_MAX ){
		    logLnew=runState->likelihood(newParams,runState->data,runState->template);
		  }
		  if(logLnew > logLmin){
			Naccepted++;
			logPriorOld=logPriorNew;
			LALInferenceCopyVariables(newParams,runState->currentParams);
			logPriorlastAccepted=logPriorNew;
			LALInferenceCopyVariables(newParams,lastAcceptedWithLikelihood);
		  }
		  else{
		   LALInferenceCopyVariables(lastAcceptedWithLikelihood,runState->currentParams);
		   logPriorOld=logPriorlastAccepted;
		  }
		}
		else{
		  logPriorOld=logPriorNew;
		  LALInferenceCopyVariables(newParams,runState->currentParams);
		}
		if(isnan(logPriorNew)){
		  XLALPrintError("Caught NaN prior slipping through sampler\n");
		  XLAL_ERROR_VOID(XLAL_EINVAL);
		}
		if(logPriorNew==-DBL_MAX || logPriorNew==DBL_MAX){
		  XLALPrintError("Caught log prior %lf slipping through sampler\n");
		  XLAL_ERROR_VOID(XLAL_EINVAL);
		}
	} while(mcmc_iter<Nmcmc*thin_factor || (0==Naccepted) );

    logLnew=runState->likelihood(runState->currentParams,runState->data,runState->template);
    if(logLnew<=logLmin) /* This point isn't above threshold, restore last one that was */
    {
        LALInferenceCopyVariables(lastAcceptedWithLikelihood,runState->currentParams);
        logPriorOld=logPriorlastAccepted;
        logLnew=runState->likelihood(runState->currentParams,runState->data,runState->template);
    }

	/* Update information to pass back out */
	LALInferenceAddVariable(runState->currentParams,"logL",(void *)&logLnew,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
	if(LALInferenceCheckVariable(runState->algorithmParams,"logZnoise")){
		tmp=logLnew-*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"logZnoise");
		LALInferenceAddVariable(runState->currentParams,"deltalogL",(void *)&tmp,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
	}
	while(data)
	{
		tmp=data->loglikelihood - data->nullloglikelihood;
		sprintf(tmpName,"deltalogl%s",data->name);
		LALInferenceAddVariable(runState->currentParams,tmpName,&tmp,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
		data=data->next;
	}
	runState->currentLikelihood=logLnew;
	LALInferenceDestroyVariables(lastAcceptedWithLikelihood);
	LALInferenceDestroyVariables(newParams);
	free(newParams);
	free(lastAcceptedWithLikelihood);
	REAL8 accept_rate=(REAL8)Naccepted/((REAL8)mcmc_iter/(REAL8)thin_factor);
	LALInferenceSetVariable(runState->algorithmParams,"accept_rate",&accept_rate);
	return;
}


void LALInferenceProposalNS(LALInferenceRunState *runState, LALInferenceVariables *parameter)
{
	
	UINT4 nIFO=0;
	LALInferenceIFOData *ifo=runState->data;
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
	UINT4 i=0;
	gsl_matrix *work=NULL;
	gsl_vector *result = NULL;
	
	/* check input arguments */
	if (!vector || !matrix || !randParam)
		XLAL_ERROR_VOID( XLAL_EFAULT );
	
	if (dim<1)
		XLAL_ERROR_VOID( XLAL_EINVAL );
	
	/* copy matrix into workspace */
	work =  gsl_matrix_alloc(dim,dim); 
	gsl_matrix_memcpy( work, matrix );
	
	/* compute the cholesky decomposition */
	gsl_linalg_cholesky_decomp(work);
	
	/* retrieve the normal distributed random numbers (LAL procedure) */
	XLALNormalDeviates( vector, randParam );
	
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
	REAL4Vector *dummy=NULL;
	REAL4 chi=0.0, factor;
	UINT4 i;
	
	/* check input arguments */
	if (!vector || !matrix || !randParam)
		XLAL_ERROR_VOID( XLAL_EFAULT );
	
	if (dim<1)
		XLAL_ERROR_VOID( XLAL_EINVAL );
	
	if (n<1)
		XLAL_ERROR_VOID( XLAL_EINVAL );
	
	
	/* first draw from MVN */
        XLALMultiNormalDeviates( vector, matrix, dim, randParam);
       
	/* then draw from chi-square with n degrees of freedom;
     this is the sum d_i*d_i with d_i drawn from a normal 
     distribution. */
        dummy = XLALCreateREAL4Vector( n );
        XLALNormalDeviates( dummy, randParam );

	/* calculate the chisquare distributed value */
	for (i=0; i<n; i++) 
	{
		chi+=dummy->data[i]*dummy->data[i];
	}
	
	/* destroy the helping vector */
        XLALDestroyREAL4Vector( dummy );
	
	/* now, finally, calculate the distribution value */
	factor=sqrt(n/chi);
	for (i=0; i<dim; i++) 
	{
		vector->data[i]*=factor;
	}
	
}


void LALInferenceProposalMultiStudentT(LALInferenceRunState *runState, LALInferenceVariables *parameter)
{
	gsl_matrix *covMat=*(gsl_matrix **)LALInferenceGetVariable(runState->proposalArgs,"covarianceMatrix");
	
	LALInferenceVariableItem *paraHead=NULL;
	REAL4Vector  *step=NULL;
	gsl_matrix *work=NULL; 
	REAL8 aii, aij, ajj;
	INT4 i, j, dim;
	RandomParams *randParam;
	UINT4 randomseed = gsl_rng_get(runState->GSLrandom);
	
	REAL8 proposal_scale=*(REAL8 *)LALInferenceGetVariable(runState->proposalArgs,"proposal_scale");
	randParam=XLALCreateRandomParams(randomseed);
	
	/* set some values */
	dim=covMat->size1;
	
	/* draw the mutinormal deviates */
	step = XLALCreateREAL4Vector(dim);
	
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
		if((paraHead->vary==LALINFERENCE_PARAM_LINEAR || paraHead->vary==LALINFERENCE_PARAM_CIRCULAR) && strcmp(paraHead->name,"rightascension") && strcmp(paraHead->name,"declination") && strcmp(paraHead->name,"time") ){
                  *(REAL8 *)paraHead->value += step->data[i];
			i++;
		}
	}

	/* LALInferenceRotateInitialPhase(parameter); */
        LALInferenceCyclicReflectiveBound(parameter,runState->priorArgs);
        
        /* destroy the vectors */
	XLALDestroyREAL4Vector(step);
	gsl_matrix_free(work);
	
	XLALDestroyRandomParams(randParam);
	/* Check boundary condition */

	return;
}


void LALInferenceProposalDifferentialEvolution(LALInferenceRunState *runState,
									   LALInferenceVariables *parameter)
	{
		LALInferenceVariables **Live=runState->livePoints;
		int i=0,j=0,same=1;//dim=0 - set but not used
		INT4 Nlive = *(INT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nlive");
		LALInferenceVariableItem *paraHead=NULL;
		LALInferenceVariableItem *paraA=NULL;
		LALInferenceVariableItem *paraB=NULL;
		
		//dim = parameter->dimension; - set but not used
		/* Select two other samples A and B*/
		i=gsl_rng_uniform_int(runState->GSLrandom,Nlive);
		/* Draw two different samples from the basket. Will loop back here if the original sample is chosen*/
	drawtwo:
		do {j=gsl_rng_uniform_int(runState->GSLrandom,Nlive);} while(j==i);
		paraHead=parameter->head;
		paraA=Live[i]->head; paraB=Live[j]->head;
		/* Add the vector B-A */
		same=1;
		for(paraHead=parameter->head,paraA=Live[i]->head,paraB=Live[j]->head;paraHead&&paraA&&paraB;paraHead=paraHead->next,paraB=paraB->next,paraA=paraA->next)
		{
			if(paraHead->vary!=LALINFERENCE_PARAM_LINEAR && paraHead->vary!=LALINFERENCE_PARAM_CIRCULAR) continue;
			*(REAL8 *)paraHead->value+=*(REAL8 *)paraB->value;
			*(REAL8 *)paraHead->value-=*(REAL8 *)paraA->value;
			if(*(REAL8 *)paraHead->value!=*(REAL8 *)paraA->value &&
			   *(REAL8 *)paraHead->value!=*(REAL8 *)paraB->value &&
			   *(REAL8 *)paraA->value!=*(REAL8 *)paraB->value) same=0;
		}
		if(same==1) goto drawtwo;
		/* Bring the sample back into bounds */
                /* LALInferenceRotateInitialPhase(parameter); */
		LALInferenceCyclicReflectiveBound(parameter,runState->priorArgs);
		return;
	}
	
static void GetCartesianPos(REAL8 vec[3],REAL8 longitude, REAL8 latitude)
{
	vec[0]=cos(longitude)*cos(latitude);
	vec[1]=sin(longitude)*cos(latitude);
	vec[2]=sin(latitude);
	return;
}

static void CartesianToSkyPos(REAL8 pos[3],REAL8 *longitude, REAL8 *latitude)
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

static void crossProduct(REAL8 out[3],REAL8 x[3],REAL8 y[3])
{
	out[0]=x[1]*y[2] - x[2]*y[1];
	out[1]=y[0]*x[2] - x[0]*y[2];
	out[2]=x[0]*y[1] - x[1]*y[0];
	return;
}

static void normalise(REAL8 vec[3]);
static void normalise(REAL8 vec[3]){
	REAL8 my_abs=0.0;
	my_abs=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
	vec[0]/=my_abs;
	vec[1]/=my_abs;
	vec[2]/=my_abs;
	return;
}

void LALInferenceRotateSky(
					   LALInferenceRunState *state,
					   LALInferenceVariables *parameter
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
	LALInferenceIFOData *ifodata1=state->data;
	while(ifodata1){
		nIFO++;
		ifodata1=ifodata1->next;
	}
	
	LALInferenceIFOData **IFOs=calloc(nIFO,sizeof(LALInferenceIFOData *));
	for(i=0,ifodata1=state->data;i<nIFO;i++){
		IFOs[i]=ifodata1;
		ifodata1=ifodata1->next;
	}
	
	
	if(nIFO<2) return;
	if(nIFO==2 && IFOs[0]==IFOs[1]) return;
	
	longi = *(REAL8 *)LALInferenceGetVariable(parameter,"rightascension");
	lat = *(REAL8 *)LALInferenceGetVariable(parameter,"declination");
	
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
	deltat+=*(REAL8 *)LALInferenceGetVariable(parameter,"time");
	LALInferenceSetVariable(parameter,"time",&deltat);	
	LALInferenceSetVariable(parameter,"declination",&newlat);
	LALInferenceSetVariable(parameter,"rightascension",&newlong);
	/*fprintf(stderr,"Skyrotate: new pos = %lf %lf %lf => %lf %lf\n",new[0],new[1],new[2],newlong,asin(new[2]));*/
	LALInferenceCyclicReflectiveBound(parameter,state->priorArgs);
	free(IFOs);
	return;
}


INT4 LALInferenceReflectDetPlane(
							 LALInferenceRunState *state,
							 LALInferenceVariables *parameter
							 )
{ /* Function to reflect a point on the sky about the plane of 3 detectors */
	/* Returns -1 if not possible */
	static LALStatus status;
	UINT4 i;
	int DetCollision=0;
	//REAL4 randnum; - set but not used
	REAL8 longi,lat;
	REAL8 dist;
	REAL8 pos[3];
	REAL8 normal[3];
	REAL8 w1[3]; /* work vectors */
	REAL8 w2[3];
	INT4 IFO1,IFO2,IFO3;
	REAL8 detvec[3];
	
	UINT4 nIFO=0;
	LALInferenceIFOData *ifodata1=state->data;
	LALInferenceIFOData *ifodata2=NULL;
	while(ifodata1){
		nIFO++;
		ifodata1=ifodata1->next;
	}
	
	LALInferenceIFOData **IFOs=calloc(nIFO,sizeof(LALInferenceIFOData *));
	if(!IFOs) {
		printf("Unable to allocate memory for %i LALInferenceIFOData *s\n",nIFO);
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
	//randnum=gsl_rng_uniform(state->GSLrandom); - set but not used
	do {
		IFO3 = gsl_rng_uniform_int(state->GSLrandom,nIFO);
	}while(IFO3==IFO1
		  || IFO3==IFO2
		  || IFOs[IFO3]==IFOs[IFO1]
		  || IFOs[IFO3]==IFOs[IFO2]);
	/*fprintf(stderr,"Using %s, %s and %s for plane\n",inputMCMC->ifoID[IFO1],inputMCMC->ifoID[IFO2],inputMCMC->ifoID[IFO3]);*/
	
	longi = *(REAL8 *)LALInferenceGetVariable(parameter,"rightascension");
	lat = *(REAL8 *)LALInferenceGetVariable(parameter,"declination");
	
	double deltalong=0;
	
	/* Convert to earth coordinates */
	SkyPosition geodetic,equatorial;
	equatorial.longitude=longi;
	equatorial.latitude=lat;
	equatorial.system=COORDINATESYSTEM_EQUATORIAL;
	geodetic.system=COORDINATESYSTEM_GEOGRAPHIC;
	LALEquatorialToGeographic(&status,&geodetic,&equatorial,&(state->data->epoch));
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
	
    /* Calculate the signed distance between the point and the plane n.(point-IFO1) */
    for(dist=0.0,i=0;i<3;i++) dist+=normal[i]*pos[i];

    /* Reflect the point pos across the plane */
    for(i=0;i<3;i++) pos[i]=pos[i]-2.0*dist*normal[i];
    
	REAL8 newLongGeo,newLat;
	CartesianToSkyPos(pos,&newLongGeo,&newLat);
	REAL8 newLongSky=newLongGeo-deltalong;
	
	
	LALInferenceSetVariable(parameter,"rightascension",&newLongSky);
	LALInferenceSetVariable(parameter,"declination",&newLat);
		
	/* Compute change in tgeocentre for this change in sky location */
	REAL8 dtold,dtnew,deltat;
	dtold = XLALTimeDelayFromEarthCenter(IFOs[0]->detector->location, longi, lat, &(IFOs[0]->epoch)); /* Compute time delay */
	dtnew = XLALTimeDelayFromEarthCenter(IFOs[0]->detector->location, newLongSky, newLat, &(IFOs[0]->epoch)); /* Compute time delay */
	deltat=dtold-dtnew; /* deltat is change in arrival time at geocentre */
	deltat+=*(REAL8 *)LALInferenceGetVariable(parameter,"time");
	LALInferenceSetVariable(parameter,"time",&deltat);
	
	LALInferenceCyclicReflectiveBound(parameter,state->priorArgs);
	free(IFOs);
	
	return(0);
}


void LALInferenceSetupLivePointsArray(LALInferenceRunState *runState){
	/* Set up initial basket of live points, drawn from prior,
	 by copying runState->currentParams to all entries in the array*/
	
	UINT4 Nlive=(UINT4)*(INT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nlive");
	UINT4 i;
	REAL8Vector *logLs;
	REAL8 logPrior=0.0;

	/* Allocate the array */
	/* runState->livePoints=XLALCalloc(Nlive,sizeof(LALVariables *)); */
	runState->livePoints=calloc(Nlive,sizeof(LALInferenceVariables *));
	if(runState->livePoints==NULL)
	{
		fprintf(stderr,"Unable to allocate memory for %i live points\n",Nlive);
		exit(1);
	}
	runState->differentialPoints=runState->livePoints;
	runState->differentialPointsLength=(size_t) Nlive;
	logLs=XLALCreateREAL8Vector(Nlive);

	LALInferenceAddVariable(runState->algorithmParams,"logLikelihoods",&logLs,LALINFERENCE_REAL8Vector_t,LALINFERENCE_PARAM_FIXED);
	fprintf(stdout,"Sprinkling %i live points, may take some time\n",Nlive);
	for(i=0;i<Nlive;i++)
	{
		runState->livePoints[i]=calloc(1,sizeof(LALInferenceVariables));
		
		/* Copy the param structure */
		LALInferenceCopyVariables(runState->currentParams,runState->livePoints[i]);
		
		/* Sprinkle the varying points among prior */
		do{
			LALInferenceDrawFromPrior( runState->livePoints[i], runState->priorArgs, runState->GSLrandom );
			logPrior=runState->prior(runState,runState->livePoints[i]);
		}while(logPrior==-DBL_MAX || isnan(logPrior));
		/* Populate log likelihood */
		logLs->data[i]=runState->likelihood(runState->livePoints[i],runState->data,runState->template);
        LALInferenceAddVariable(runState->livePoints[i],"logL",(void *)&(logLs->data[i]),LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
	LALInferenceAddVariable(runState->livePoints[i],"logPrior",(void*)&logPrior,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
	}
	
}


void LALInferenceSetupkDTreeNSLivePoints( LALInferenceRunState *runState ){
  /* create a k-d tree from the nested sampling live points */
  LALInferenceKDTree *tree;
  REAL8 *low = NULL, *high = NULL; /* upper and lower bounds of tree */
  size_t ndim = 0;
  LALInferenceVariableItem *currentItem;
  UINT4 cnt = 0;
  REAL8 *pt = NULL;
  LALInferenceVariables *template =
    XLALCalloc(1,sizeof(LALInferenceVariables));
  UINT4 Nlive = *(UINT4 *)LALInferenceGetVariable( runState->algorithmParams,
                                                   "Nlive" );
  
  /* if a current tree exists remove it */
  if ( LALInferenceCheckVariable( runState->proposalArgs, "kDTree" ) )
  {
    LALInferenceKDTreeDelete(*(LALInferenceKDTree **)LALInferenceGetVariable(runState->proposalArgs,"kDTree"));
    LALInferenceRemoveVariable( runState->proposalArgs, "kDTree" );
  }
  /* get the upper and lower bounds for each parameter */
  currentItem = runState->currentParams->head;
  while ( currentItem != NULL ) {
    if ( currentItem->vary != LALINFERENCE_PARAM_FIXED &&
         currentItem->vary != LALINFERENCE_PARAM_OUTPUT ) {
      if( LALInferenceCheckMinMaxPrior( runState->priorArgs,
                                        currentItem->name ) ){
        cnt++;
         
        low = XLALRealloc(low, sizeof(REAL8)*cnt);
        high = XLALRealloc(high, sizeof(REAL8)*cnt);
        pt = XLALRealloc(pt, sizeof(REAL8)*cnt);
                
        LALInferenceGetMinMaxPrior( runState->priorArgs, currentItem->name,
                                    &(low[cnt-1]), &(high[cnt-1]) );
      }
      else if( LALInferenceCheckGaussianPrior( runState->priorArgs,
                                               currentItem->name ) ){
        REAL8 mn, stddiv;
        REAL8 livelow, livehigh, difflh;
        
        cnt++;
        
        low = XLALRealloc(low, sizeof(REAL8)*cnt);
        high = XLALRealloc(high, sizeof(REAL8)*cnt);
        pt = XLALRealloc(pt, sizeof(REAL8)*cnt);
        
        LALInferenceGetGaussianPrior( runState->priorArgs, currentItem->name,
                                      &mn, &stddiv );
        
        /* find the maximum and minimum live point values */
        getMinMaxLivePointValue( runState->livePoints, currentItem->name, Nlive,
                                 &livelow, &livehigh );
        difflh = livehigh - livelow;
        
        /* to add a bit of room at either side add on half the difference */
        low[cnt-1] = livelow - difflh/2.;
        high[cnt-1] = livehigh + difflh/2.;
      }
    }
                      
    currentItem = currentItem->next;
  }

  ndim = (size_t)cnt;
          
  /* set up tree */
  tree = LALInferenceKDEmpty( low, high, ndim );
  LALInferenceCopyVariables( runState->currentParams, template );
                    
  /* add points to tree */
  for( cnt = 0; cnt < Nlive; cnt++ ){
    LALInferenceKDVariablesToREAL8( runState->livePoints[cnt], pt, template );
    LALInferenceKDAddPoint( tree, pt );
  }
                  
  /* add tree */
  LALInferenceAddVariable( runState->proposalArgs, "kDTree", &tree,
                           LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_FIXED );
  
  /* if template doesn't exist add it */
  if ( !LALInferenceCheckVariable( runState->proposalArgs,
                                   "kDTreeVariableTemplate" ) ){
    LALInferenceAddVariable( runState->proposalArgs, "kDTreeVariableTemplate",
                             &template, LALINFERENCE_void_ptr_t,
                             LALINFERENCE_PARAM_FIXED );
  }
  
  XLALFree( high );
  XLALFree( low );
  XLALFree( pt );
}
