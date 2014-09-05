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
#include <signal.h>

#ifdef HAVE_LIBLALXML
#include <lal/LALInferenceXML.h>
#endif

#define PROGRAM_NAME "LALInferenceNestedSampler.c"
#define CVS_ID_STRING "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define CVS_NAME_STRING "$Name$"

#define MAX_MCMC 5000 /* Maximum chain length, set to be higher than expected from a reasonable run */
#define ACF_TOLERANCE 0.01 /* Desired maximum correlation of MCMC samples */

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

static int __chainfile_iter;

/**
 * structure holding internal state of the NS integrator
 */
typedef struct tagNSintegralState
{
  UINT4 size,iteration;
  REAL8 *logZarray,
	*oldZarray,
	*Harray,
	*logwarray,
	*logtarray,
	*logt2array;
} NSintegralState;

static struct itimerval checkpoint_timer;

/** Write the current state to a checkpoint file the given filename */
static int ReadNSCheckPoint(CHAR *filename, LALInferenceRunState *runState, NSintegralState *s);
/** Read the given filename to populate a given LALInferenceRunState and NSintegralState */
static int WriteNSCheckPoint(CHAR *filename, LALInferenceRunState *runState, NSintegralState *s);

/* This is checked by the main loop to determine when to checkpoint */
static volatile sig_atomic_t __ns_saveStateFlag = 0;
/* This indicates the main loop should terminate */
static volatile sig_atomic_t __ns_exitFlag = 0;

/* Signal handler for SIGINT, which is produced by condor
 * when evicting a job, or when the user presses Ctrl-C */
static void catch_interrupt(UNUSED int sig, UNUSED siginfo_t *siginfo,UNUSED void *context);
static void catch_interrupt(UNUSED int sig, UNUSED siginfo_t *siginfo,UNUSED void *context)
{
  
 __ns_saveStateFlag=1;
 __ns_exitFlag=1;
}

/* Signal handler for SIGALRM, for periodic checkpointing */
static void catch_alarm(UNUSED int sig, UNUSED siginfo_t *siginfo,UNUSED void *context);
static void catch_alarm(UNUSED int sig, UNUSED siginfo_t *siginfo,UNUSED void *context)
{
  __ns_saveStateFlag=1;
}

static UINT4 UpdateNMCMC(LALInferenceRunState *runState);
/* Prototypes for private "helper" functions. */
//static void SamplePriorDiscardAcceptance(LALInferenceRunState *runState);
static double logadd(double a,double b);
static REAL8 mean(REAL8 *array,int N);
static void getMinMaxLivePointValue( LALInferenceVariables **livepoints, 
                                     const CHAR *pname, UINT4 Nlive, 
                                     REAL8 *minval, REAL8 *maxval );

static double logadd(double a,double b){
	if(a>b) return(a+log(1.0+exp(b-a)));
	else return(b+log(1.0+exp(a-b)));
}

static UINT4 checkForSingleAdapt(LALInferenceRunState *runState);
static void SetupEigenProposals(LALInferenceRunState *runState);

/**
 * Update the internal state of the integrator after receiving the lowest logL
 * value logL
 */
REAL8 incrementEvidenceSamples(LALInferenceRunState *runState, REAL8 logL, NSintegralState *s);
REAL8 incrementEvidenceSamples(LALInferenceRunState *runState, REAL8 logL, NSintegralState *s)
{
  REAL8 Wtarray[s->size];
  UINT4 Nlive=*(UINT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nlive");
  /* Update evidence array */
  for(UINT4 j=0;j<s->size;j++){
    Wtarray[j]=s->logwarray[j]+logL+logadd(0,s->logt2array[j]*s->logtarray[j])-log(2.0);
    s->logZarray[j]=logadd(s->logZarray[j],Wtarray[j]);
    s->Harray[j]= exp(Wtarray[j]-s->logZarray[j])*logL
    + exp(s->oldZarray[j]-s->logZarray[j])*(s->Harray[j]+s->oldZarray[j])-s->logZarray[j];
    REAL8 tmp=s->logtarray[j];
    s->logtarray[j]=s->logt2array[j];
    s->logt2array[j]=tmp;
    s->logtarray[j]=LALInferenceNSSample_logt(Nlive,runState->GSLrandom);
    s->logwarray[j]+=s->logtarray[j];
  }
  s->iteration++;
  return(mean(s->logZarray,s->size));
}

static void printAdaptiveJumpSizes(FILE *file, LALInferenceRunState *runState);
static void printAdaptiveJumpSizes(FILE *file, LALInferenceRunState *runState)
{
    LALInferenceVariableItem *this=runState->currentParams->head;
    REAL8 *val=NULL;
    char tmpname[1000]="";
    fprintf(file,"Adaptive proposal step size:\n");
    while(this)
    {
        sprintf(tmpname,"%s_%s",this->name,ADAPTSUFFIX);
        if(LALInferenceCheckVariable(runState->proposalArgs,tmpname))
        {
            val=(REAL8 *)LALInferenceGetVariable(runState->proposalArgs,tmpname);
            fprintf(file,"%s: %lf\n",this->name,*val);
        }
        this=this->next;
    }

}

static void resetProposalStats(LALInferenceRunState *runState);
static void resetProposalStats(LALInferenceRunState *runState)
{
    LALInferenceProposalStatistics *propStat;
    LALInferenceVariableItem *this;
    this = runState->proposalStats->head;
    while(this){
        propStat = (LALInferenceProposalStatistics *)this->value;
        propStat->accepted = 0;
        propStat->proposed = 0;
        this = this->next;
    } 
}



/* Create Internal arrays for sampling the integral */
static NSintegralState *initNSintegralState(UINT4 Nruns, UINT4 Nlive);
static NSintegralState *initNSintegralState(UINT4 Nruns, UINT4 Nlive)
{
  NSintegralState *s=XLALMalloc(sizeof(NSintegralState));
  s->iteration=0;
  s->size=Nruns;
  s->logZarray = XLALCalloc(Nruns,sizeof(REAL8));
  s->oldZarray = XLALCalloc(Nruns,sizeof(REAL8));
  s->Harray = XLALCalloc(Nruns,sizeof(REAL8));
  s->logwarray = XLALCalloc(Nruns,sizeof(REAL8));
  s->logtarray=XLALCalloc(Nruns,sizeof(REAL8));
  s->logt2array=XLALCalloc(Nruns,sizeof(REAL8));
  REAL8 logw=log(1.0-exp(-1.0/Nlive));

  if(s->logZarray==NULL || s->Harray==NULL || s->oldZarray==NULL || s->logwarray==NULL)
  {fprintf(stderr,"Unable to allocate RAM\n"); exit(-1);}
  for(UINT4 i=0;i<Nruns;i++)  {s->logwarray[i]=logw; s->logZarray[i]=-DBL_MAX;
    s->oldZarray[i]=-DBL_MAX; s->Harray[i]=0.0;s->logtarray[i]=-1.0/Nlive;
    s->logt2array[i]=-1.0/Nlive;
  }
  return s;
}

static REAL8 mean(REAL8 *array,int N){
	REAL8 sum=0.0;
	int i;
	for(i=0;i<N;i++) sum+=array[i];
	return sum/((REAL8) N);
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
  
REAL8 LALInferenceNSSample_logt(int Nlive,gsl_rng *RNG){
	REAL8 t=0.0;
	REAL8 a=0.0;
	while((Nlive--)>1) {a=gsl_rng_uniform(RNG); t = t>a ? t : a;}
	return(log(t));
}

static UINT4 UpdateNMCMC(LALInferenceRunState *runState){
	INT4 max = 0;
	INT4 maxMCMC = MAX_MCMC;
	/* Measure Autocorrelations if the Nmcmc is not over-ridden */
	if(!LALInferenceGetProcParamVal(runState->commandLine,"--Nmcmc") && !LALInferenceGetProcParamVal(runState->commandLine,"--nmcmc")){
        if(LALInferenceCheckVariable(runState->algorithmParams,"maxmcmc"))
            maxMCMC = *(INT4 *)LALInferenceGetVariable(runState->algorithmParams,"maxmcmc");
        if(LALInferenceCheckVariable(runState->algorithmParams,"Nmcmc")) /* if already estimated the length */
            max=4 * *(INT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nmcmc"); /* We will use this to go out 4x last ACL */
        else max=4*maxMCMC; /* otherwise use the MAX_MCMC */
        if(max>4*maxMCMC) max=4*maxMCMC;
	UINT4 Nlive = *(INT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nlive");
	UINT4 rnd=gsl_rng_uniform_int(runState->GSLrandom,Nlive);
	LALInferenceCopyVariables(runState->livePoints[rnd],runState->currentParams);
        LALInferenceVariables *acls=LALInferenceComputeAutoCorrelation(runState, max, LALInferenceNestedSamplingSloppySample) ;
        max=10;
        for(LALInferenceVariableItem *this=acls->head;this;this=this->next) {
            if(LALInferenceCheckVariable(runState->algorithmParams,"verbose"))
                fprintf(stdout,"Autocorrelation length of %s: %i\n",this->name,(INT4) *(REAL8 *)this->value);
            if(*(REAL8 *)this->value>max) {
                max=(INT4) *(REAL8 *)this->value;
            }
        }
        LALInferenceClearVariables(acls);
        XLALFree(acls);
        if(max>maxMCMC){
            fprintf(stderr,"Warning: Estimated chain length %i exceeds maximum %i!\n",max,maxMCMC);
            max=maxMCMC;
        }
        LALInferenceSetVariable(runState->algorithmParams,"Nmcmc",&max);
    }
    return(max);
}

/* estimateCovarianceMatrix reads the list of live points,
 and works out the covariance matrix of the varying parameters
 - CIRCULAR parameters are wrapped around before the calculation and must be
 scaled to the range 0 -> 2pi. Only works with REAL8 variables */
void LALInferenceNScalcCVM(gsl_matrix **cvm, LALInferenceVariables **Live, UINT4 Nlive)
{
	UINT4 i,j,k;
	UINT4 ND=0;
	LALInferenceVariableItem *item,*k_item,*j_item;
	REAL8 *means, *ms, *mc,jval=0.,kval=0.;
	
	/* Find the number of dimensions which vary in the covariance matrix */
	for(item=Live[0]->head;item!=NULL;item=item->next)
		if((item->vary==LALINFERENCE_PARAM_LINEAR || item->vary==LALINFERENCE_PARAM_CIRCULAR) && item->type==LALINFERENCE_REAL8_t) ND++;
	
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
	if(NULL==(means = XLALMalloc((size_t)ND*sizeof(REAL8)))){fprintf(stderr,"Can't allocate RAM"); exit(-1);}
	if(NULL==(ms = XLALMalloc((size_t)ND*sizeof(REAL8)))){fprintf(stderr,"Can't allocate RAM"); exit(-1);}
	if(NULL==(mc = XLALMalloc((size_t)ND*sizeof(REAL8)))){fprintf(stderr,"Can't allocate RAM"); exit(-1);}
	for(i=0;i<ND;i++){ 
          means[i]=0.0;
          ms[i] = 0.;
          mc[i] = 0.;
        }
        for(item=Live[0]->head,j=0;item;item=item->next){
	  if((item->vary!=LALINFERENCE_PARAM_CIRCULAR && item->vary!=LALINFERENCE_PARAM_LINEAR) || item->type!=LALINFERENCE_REAL8_t) continue;
	    for(i=0;i<Nlive;i++){
			void *ptr = LALInferenceGetVariable(Live[i],item->name);
			switch(item->vary)
			{
			  case LALINFERENCE_PARAM_LINEAR:
			  {
                                means[j]+=*(REAL8 *)ptr;
				break;
			  }
			  case LALINFERENCE_PARAM_CIRCULAR:
			  {
                               ms[j] += sin(*(REAL8 *)ptr);
                               mc[j] += cos(*(REAL8 *)ptr);
                               break;
			  }
			  default:
			  {
			    break;
			  }
			}
	  }
	  j++;
	}

        for(item=Live[0]->head,j=0;item;item=item->next){ 
          if( item->vary==LALINFERENCE_PARAM_LINEAR && item->type==LALINFERENCE_REAL8_t ){
            means[j]/=(REAL8)Nlive;
            j++;
          }
          if( item->vary==LALINFERENCE_PARAM_CIRCULAR && item->type==LALINFERENCE_REAL8_t ){
            ms[j]/=(REAL8)Nlive;
            mc[j]/=(REAL8)Nlive;
            
            means[j] = atan2(ms[j], mc[j]);
            means[j] = means[j]<0? 2.0*LAL_PI + means[j] : means[j];
            
            j++;
          }
        }
        
        XLALFree(ms);
        XLALFree(mc);
        
	/* Iterate over the matrix elements one at a time */
	j_item=k_item=Live[0]->head;
	for(j=0,j_item=Live[0]->head;j_item;j_item=j_item->next)
	{
	  /* Skip non-varying parameters */
	  if((j_item->vary!=LALINFERENCE_PARAM_LINEAR && j_item->vary!=LALINFERENCE_PARAM_CIRCULAR )||j_item->type!=LALINFERENCE_REAL8_t ) {continue;}

	  for(k_item=Live[0]->head,k=0;k_item &&k<=j;k_item=k_item->next)
	  {
	    /* Skip non-varying parameters */
	    if((k_item->vary!=LALINFERENCE_PARAM_LINEAR && k_item->vary!=LALINFERENCE_PARAM_CIRCULAR)||k_item->type!=LALINFERENCE_REAL8_t) {continue;}
	    
	    /* Loop over live points updating covariance elements */
	    for(i=0;i<Nlive;i++)
	    {
	      void *jptr=LALInferenceGetVariable(Live[i],j_item->name);
	      void *kptr=LALInferenceGetVariable(Live[i],k_item->name);
	      
	      /* Check for linear or circular */
	      if( j_item->vary==LALINFERENCE_PARAM_CIRCULAR )
		jval = LALInferenceAngularDistance(*(REAL8 *)jptr, means[j]);
	      else if( j_item->vary==LALINFERENCE_PARAM_LINEAR )
		jval = *(REAL8 *)jptr - means[j];
	      
	      if( k_item->vary==LALINFERENCE_PARAM_CIRCULAR )
		kval = LALInferenceAngularDistance(*(REAL8 *)kptr, means[k]);
	      else if( k_item->vary==LALINFERENCE_PARAM_LINEAR )
		kval = *(REAL8 *)kptr - means[k];
	      
	      gsl_matrix_set(*cvm,j,k, gsl_matrix_get(*cvm,j,k) + kval*jval);
	      
	    }
	    
	    k++;
	  } /* end loop over k */
	  j++;
	}/* end loop over j */
	
	/* Normalise */
	for(i=0;i<ND;i++) for(j=0;j<ND;j++) gsl_matrix_set(*cvm,i,j,gsl_matrix_get(*cvm,i,j)/((REAL8) Nlive));
	XLALFree(means);
        
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
  REAL8 *logZarray,*oldZarray,*Harray,*logwarray,*logtarray;
  REAL8 TOLERANCE=0.1;
  REAL8 logZ,logZnew,logLmin,logLmax=-DBL_MAX,logLtmp,logw,H,logZnoise,dZ=0;//deltaZ - set but not used
  LALInferenceVariables *temp;
  FILE *fpout=NULL;
  REAL8 dblmax=-DBL_MAX;
  REAL8 zero=0.0;
  REAL8 *logLikelihoods=NULL;
  UINT4 verbose=0;
  REAL8 sloppyfrac;
  UINT4 displayprogress=0;
  LALInferenceVariables *currentVars=XLALCalloc(1,sizeof(LALInferenceVariables));
  REAL8 kdupdate=0.;
  UINT4 samplePrior=0; //If this flag is set to a positive integer, code will just draw this many samples from the prior
  ProcessParamsTable *ppt=NULL;
  
  /* Default sample logging functions with and without XML */
  #ifdef HAVE_LIBLALXML
  char *outVOTable=NULL;
  if(!runState->logsample) runState->logsample=LALInferenceLogSampleToArray;
  #else
  if(!runState->logsample) runState->logsample=LALInferenceLogSampleToFile;
  #endif
  
  if ( !LALInferenceCheckVariable(runState->algorithmParams, "logZnoise" ) ){
    logZnoise=LALInferenceNullLogLikelihood(runState->data);
    
    LALInferenceAddVariable(runState->algorithmParams,"logZnoise",&logZnoise,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
  }
  else{
    logZnoise = 
      *(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"logZnoise");
  }

  logLikelihoods=(REAL8 *)(*(REAL8Vector **)LALInferenceGetVariable(runState->algorithmParams,"logLikelihoods"))->data;
  
  if((ppt=LALInferenceGetProcParamVal(runState->commandLine,"--sampleprior")))
  {
    samplePrior=atoi(ppt->value);
    fprintf(stdout,"Generating %i samples from the prior\n",samplePrior);
  }

  verbose=LALInferenceCheckVariable(runState->algorithmParams,"verbose");
  displayprogress=verbose;
  
  /* Operate on parallel runs if requested */
  if(LALInferenceCheckVariable(runState->algorithmParams,"Nruns"))
    Nruns = *(UINT4 *) LALInferenceGetVariable(runState->algorithmParams,"Nruns");

  /* Create workspace for arrays */
  NSintegralState *s=NULL;

  if(LALInferenceCheckVariable(runState->algorithmParams,"tolerance"))
    TOLERANCE = *(REAL8 *) LALInferenceGetVariable(runState->algorithmParams,"tolerance");
  
  /* Check that necessary parameters are created */
  if(!LALInferenceCheckVariable(runState->algorithmParams,"logLmin"))
    LALInferenceAddVariable(runState->algorithmParams,"logLmin",&dblmax,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
  
  if(!LALInferenceCheckVariable(runState->algorithmParams,"accept_rate"))
    LALInferenceAddVariable(runState->algorithmParams,"accept_rate",&zero,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
  if(!LALInferenceCheckVariable(runState->algorithmParams,"sub_accept_rate"))
    LALInferenceAddVariable(runState->algorithmParams,"sub_accept_rate",&zero,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
  if(!LALInferenceCheckVariable(runState->algorithmParams,"sloppyfraction"))
    LALInferenceAddVariable(runState->algorithmParams,"sloppyfraction",&zero,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
  
  /* Set up the proposal scale factor, for use in the multi-student jump step */
  REAL8 propScale = 0.1;
  LALInferenceAddVariable(runState->proposalArgs,"proposal_scale",&propScale,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
  
  /* Open output file */
  ppt = NULL; 
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

  minpos=0;
  
  logw=log(1.0-exp(-1.0/Nlive));
  i=0;
  /* sort points for consistent order before creating the matrix */
  for(i=0;i<Nlive;i++) 
    LALInferenceSortVariablesByName(runState->livePoints[i]);

  /* Set up eigenvector proposals */
  SetupEigenProposals(runState);
  
  /* Use the live points as differential evolution points */
  runState->differentialPoints=runState->livePoints;
  runState->differentialPointsLength=(size_t) Nlive;
  
  /* If there is not a proposal counter, put one into the variables, initialized to zero. */
  if (!LALInferenceCheckVariable(runState->proposalArgs, cycleArrayCounterName)) {
    i = 0;
    LALInferenceAddVariable(runState->proposalArgs, cycleArrayCounterName, &i, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_CIRCULAR);
  }
  /* set up k-D tree if required and not already set */
  if ( ( LALInferenceGetProcParamVal(runState->commandLine,"--kDTree") ||
    LALInferenceGetProcParamVal(runState->commandLine,"--kdtree")) &&
    !LALInferenceCheckVariable( runState->proposalArgs, "kDTree" ) )
      LALInferenceSetupkDTreeNSLivePoints( runState );
      
  if(!LALInferenceCheckVariable(runState->algorithmParams,"Nmcmc")){
    INT4 tmp=200;
    LALInferenceAddVariable(runState->algorithmParams,"Nmcmc",&tmp,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_OUTPUT);
  }
  s=initNSintegralState(Nruns,Nlive);
  
  
  /* Check for an interrupted run */
  char resumefilename[FILENAME_MAX];
  sprintf(resumefilename,"%s_resume",outfile);
  char outfilebackup[FILENAME_MAX];
  sprintf(outfilebackup,"%s.bak",outfile);
  int retcode=1;
  if(LALInferenceGetProcParamVal(runState->commandLine,"--resume")){
      retcode=ReadNSCheckPoint(resumefilename,runState,s);
      if(retcode==0){
          for(i=0;i<Nlive;i++) logLikelihoods[i]=*(REAL8 *)LALInferenceGetVariable(runState->livePoints[i],"logL");
          iter=s->iteration;
          /* back up the old file and get ready to append to it from the right place */
          char commandline[1024];
          sprintf(commandline,"cp %s %s",outfile,outfilebackup);
          system(commandline);
          sprintf(commandline,"head -n %i %s > %s",iter+1,outfilebackup,outfile);
          system(commandline);
          fpout=fopen(outfile,"a");
      }
      /* Install a periodic alarm that will trigger a checkpoint */
      int sigretcode=0;
      struct sigaction sa;
      sa.sa_sigaction=catch_alarm;
      sa.sa_flags=SA_SIGINFO;
      sigretcode=sigaction(SIGVTALRM,&sa,NULL);
      if(sigretcode!=0) fprintf(stderr,"WARNING: Cannot establish checkpoint timer!\n");
      /* Condor sends SIGUSR2 to checkpoint and continue */
      sigretcode=sigaction(SIGUSR2,&sa,NULL);
      if(sigretcode!=0) fprintf(stderr,"WARNING: Cannot establish checkpoint on SIGUSR2.\n");
      checkpoint_timer.it_interval.tv_sec=4*3600; /* Default timer 4 hours */
      checkpoint_timer.it_interval.tv_usec=0;
      checkpoint_timer.it_value=checkpoint_timer.it_interval;
      setitimer(ITIMER_VIRTUAL,&checkpoint_timer,NULL);
      /* Install the handler for the condor interrupt signal */
      sa.sa_sigaction=catch_interrupt;
      sigretcode=sigaction(SIGINT,&sa,NULL);
      if(sigretcode!=0) fprintf(stderr,"WARNING: Cannot establish checkpoint on SIGINT.\n");
      /* Condor sends SIGTERM to vanilla universe jobs to evict them */
      sigretcode=sigaction(SIGTERM,&sa,NULL);
      if(sigretcode!=0) fprintf(stderr,"WARNING: Cannot establish checkpoint on SIGTERM.\n");
      /* Condor sends SIGTSTP to standard universe jobs to evict them.
       *I think condor handles this, so didn't add a handler CHECK */
  }
  
  if(retcode!=0)
  {
    if(LALInferenceGetProcParamVal(runState->commandLine,"--resume"))
        fprintf(stdout,"Unable to open resume file %s. Starting anew.\n",resumefilename);
    /* Sprinkle points */
    LALInferenceSetVariable(runState->algorithmParams,"logLmin",&dblmax);
    for(i=0;i<Nlive;i++) {
	runState->currentParams=runState->livePoints[i];
	do
    {
        runState->evolve(runState);
    	logLikelihoods[i]=runState->likelihood(runState->livePoints[i],runState->data,runState->model);
    }while(isnan(logLikelihoods[i]));
	if(XLALPrintProgressBar((double)i/(double)Nlive)) fprintf(stderr,"\n");
    }
    fpout=fopen(outfile,"w");
  }
  if(fpout==NULL) {fprintf(stderr,"Unable to open output file %s!\n",outfile); exit(1);}
  else{
    if(setvbuf(fpout,NULL,_IOFBF,0x100000)) /* Set buffer to 1MB so as to not thrash NFS */
      fprintf(stderr,"Warning: Unable to set output file buffer!");
    LALInferenceAddVariable(runState->algorithmParams,"outfile",&fpout,LALINFERENCE_void_ptr_t,LALINFERENCE_PARAM_FIXED);
  }

  logZarray=s->logZarray;
  logtarray=s->logtarray;
  oldZarray=s->oldZarray;
  logwarray=s->logwarray;
  Harray=s->Harray;
  
  /* Find maximum likelihood and sanity check */
  for(i=0;i<Nlive;i++)
  {
    LALInferenceAddVariable(runState->livePoints[i],"logw",&logw,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
    logLtmp=logLikelihoods[i];
    logLmax=logLtmp>logLmax? logLtmp : logLmax;
    if(isnan(logLikelihoods[i])) {
      fprintf(stderr,"Detected logL[%i]=%lf! Sanity checking...\n",i,logLikelihoods[i]);
      if(LALInferenceSanityCheck(runState))
	exit(1);
    }
  }
  /* sort points for consistent order after attaching delta parameters */
  for(i=0;i<Nlive;i++) 
    LALInferenceSortVariablesByName(runState->livePoints[i]);

  /* Update the covariance matrix for proposal distribution */
  SetupEigenProposals(runState);
//  LALInferenceNScalcCVM(cvm,runState->livePoints,Nlive);

  /* re-calculate the k-D tree from the new points if required */
  if ( LALInferenceCheckVariable( runState->proposalArgs, "kDTree" ) ){
    LALInferenceSetupkDTreeNSLivePoints( runState );
    
    /* get k-d tree update rate (this is how often the tree gets updated
     * as a factor the number of live points - default is 4 */
    if( LALInferenceGetProcParamVal( runState->commandLine, 
      "--kDTreeUpdateFactor") ){
      kdupdate = atof( LALInferenceGetProcParamVal( runState->commandLine,
						    "--kDTreeUpdateFactor")->value );
      }else
	kdupdate = 4.;
  }
  
  /* Reset proposal stats before starting */
  resetProposalStats(runState);
  
  /* Set the number of MCMC points */
  UpdateNMCMC(runState);
  /* Output some information */
  if(verbose){
    LALInferencePrintProposalStatsHeader(stdout,runState->proposalStats);
    LALInferencePrintProposalStats(stdout,runState->proposalStats);
    resetProposalStats(runState);
    printAdaptiveJumpSizes(stdout, runState);
  }
  /* Write out names of parameters */
  FILE *lout=NULL;
  char param_list[FILENAME_MAX];
  sprintf(param_list,"%s_params.txt",outfile);
  lout=fopen(param_list,"w");
  minpos=0;
  LALInferenceFprintParameterHeaders(lout,runState->livePoints[0]);
  fclose(lout);
  runState->currentParams=currentVars;
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
    if(samplePrior) logLmin=-DBL_MAX;

    logZnew=incrementEvidenceSamples(runState, logLikelihoods[minpos], s);
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
  
  if (runState->currentLikelihood>logLmax)
    logLmax=runState->currentLikelihood;

  logw=mean(logwarray,Nruns);
  LALInferenceAddVariable(runState->livePoints[minpos],"logw",&logw,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
  dZ=logadd(logZ,logLmax-((double) iter)/((double)Nlive))-logZ;
  sloppyfrac=*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"sloppyfraction");
  if(displayprogress) fprintf(stderr,"%i: accpt: %1.3f Nmcmc: %i sub_accpt: %1.3f slpy: %2.1f%% H: %3.2lf nats logL:%.3lf ->%.3lf logZ: %.3lf deltalogLmax: %.2lf dZ: %.3lf Zratio: %.3lf \n",\
    iter,\
    *(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"accept_rate")/(REAL8)itercounter,\
    *(INT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nmcmc"),\
    *(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"sub_accept_rate"),\
    100.0*sloppyfrac,\
    H,\
    logLmin,\
    runState->currentLikelihood,\
    logZ,\
    (logLmax - LALInferenceGetREAL8Variable(runState->algorithmParams,"logZnoise")), \
    dZ,\
    ( logZ - LALInferenceGetREAL8Variable(runState->algorithmParams,"logZnoise"))\
  );
  iter++;
  
  /* Save progress */
  if(__ns_saveStateFlag!=0)
    {
      if(__ns_exitFlag) fprintf(stdout,"Saving state to %s.\n",resumefilename);
      WriteNSCheckPoint(resumefilename,runState,s);
      fflush(fpout);
      __ns_saveStateFlag=0;
    }
   /* Have we been told to quit? */
  if(__ns_exitFlag) {
    fclose(fpout);
    exit(0);
  }
  
  /* Update the proposal */
  if(!(iter%(Nlive/4))) {
    /* Update the covariance matrix */
    if ( LALInferenceCheckVariable( runState->proposalArgs,"covarianceMatrix" ) ){
      SetupEigenProposals(runState);
    }
  
    /* Update NMCMC from ACF */
    UpdateNMCMC(runState);
  
    /* Output some information */
    if(verbose){
      LALInferencePrintProposalStatsHeader(stdout,runState->proposalStats);
      LALInferencePrintProposalStats(stdout,runState->proposalStats);
      resetProposalStats(runState);
      printAdaptiveJumpSizes(stdout, runState);
    }
  }
  
  if ( LALInferenceCheckVariable( runState->proposalArgs,"kDTree" )){
    /* update k-d tree */
    if(!(iter%((int)floor((REAL8)Nlive * kdupdate))))
      LALInferenceSetupkDTreeNSLivePoints( runState ); 
  }
  }
  while(samplePrior?((Nlive+iter)<samplePrior):( iter <= Nlive ||  dZ> TOLERANCE)); /* End of NS loop! */
    
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
    if(output_array) XLALFree(output_array);      
    #endif
   
    /* Clean up resume file */
    if(LALInferenceGetProcParamVal(runState->commandLine,"--resume"))
    {
      if(!access(resumefilename,W_OK)) remove(resumefilename);
      if(!access(outfilebackup,W_OK)) remove(outfilebackup);
    }
    
    /* Free memory */
    XLALFree(logtarray); XLALFree(logwarray); XLALFree(logZarray);
}

/* Calculate the autocorrelation function of the sampler (runState->evolve) for each parameter
 * Evolves the sample starting with the value passed in temp, with a maximum of max_iterations steps.
 Return the ACL for each parameter as a LALInferenceVariables */
LALInferenceVariables *LALInferenceComputeAutoCorrelation(LALInferenceRunState *runState, UINT4 max_iterations, LALInferenceEvolveOneStepFunction evolve)
{
  ProcessParamsTable *ppt=NULL;
  char chainfilename[128]="";
  char acf_file_name[128]="";
  FILE *chainfile=NULL;
  FILE *acffile=NULL;
  UINT4 i,j;
  UINT4 nPar=0; // = LALInferenceGetVariableDimensionNonFixed(runState->currentParams);
  REAL8 **data_array=NULL;
  REAL8 **acf_array=NULL;
  LALInferenceVariableItem *this;
  INT4 thinning=10;
  max_iterations/=thinning;
  /* Find the number and names of variables */
  for(this=runState->currentParams->head;this;this=this->next) if(this->vary!=LALINFERENCE_PARAM_FIXED && this->vary!=LALINFERENCE_PARAM_OUTPUT && this->type==LALINFERENCE_REAL8_t) nPar++;
  char **param_names=XLALCalloc(nPar,sizeof(char *));
  for(i=0,this=runState->currentParams->head;this;this=this->next) if(this->vary!=LALINFERENCE_PARAM_FIXED && this->vary!=LALINFERENCE_PARAM_OUTPUT && this->type==LALINFERENCE_REAL8_t) param_names[i++]=this->name;

  REAL8 ACF,ACL,max=0;
  LALInferenceVariables *acls=XLALCalloc(1,sizeof(LALInferenceVariables));

  /* Back up the algorithm state and replace with a clean version for logSampletoarray */
  LALInferenceVariables myAlgParams,*oldAlgParams=runState->algorithmParams;
  LALInferenceVariables myCurrentParams,*oldCurrentParams=runState->currentParams;
  memset(&myAlgParams,0,sizeof(LALInferenceVariables));
  memset(&myCurrentParams,0,sizeof(LALInferenceVariables));
  LALInferenceCopyVariables(oldAlgParams,&myAlgParams);
  if(LALInferenceCheckVariable(&myAlgParams,"outputarray")) LALInferenceRemoveVariable(&myAlgParams,"outputarray");
  if(LALInferenceCheckVariable(&myAlgParams,"N_outputarray")) LALInferenceRemoveVariable(&myAlgParams,"N_outputarray");
  LALInferenceRemoveVariable(&myAlgParams,"outfile");
  LALInferenceRemoveVariable(&myAlgParams,"Nmcmc");
  LALInferenceAddVariable(&myAlgParams,"Nmcmc",&thinning,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_OUTPUT);

  LALInferenceSortVariablesByName(&myCurrentParams);
  runState->algorithmParams=&myAlgParams;
  runState->currentParams=&myCurrentParams;
  LALInferenceVariables **livePoints=runState->livePoints;
  UINT4 Nlive = *(INT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nlive");
  UINT4 BAILOUT=100; /* this should be the same as the bailout in the sampler */
  REAL8 accept=0.0;
  /* We can record write the MCMC chain to a file too */
  ppt=LALInferenceGetProcParamVal(runState->commandLine,"--acf-chainfile");
  if(ppt){
    sprintf(chainfilename,"%s.%i",ppt->value,__chainfile_iter);
    chainfile=fopen(chainfilename,"w");
    LALInferenceCopyVariables(livePoints[0],&myCurrentParams);
    LALInferenceSortVariablesByName(&myCurrentParams);
    for(this=myCurrentParams.head;this;this=this->next) fprintf(chainfile,"%s ",this->name);
    fprintf(chainfile,"\n");
    LALInferenceAddVariable(&myAlgParams,"outfile",&chainfile,LALINFERENCE_void_ptr_t,LALINFERENCE_PARAM_FIXED);
  }
  ppt=LALInferenceGetProcParamVal(runState->commandLine,"--acf-file");
  if(ppt){
    sprintf(acf_file_name,"%s.%i",ppt->value,__chainfile_iter);
    acffile=fopen(acf_file_name,"w");
  }
  __chainfile_iter++;
  do{ /* Pick a random sample that isn't trapped in some corner*/
	UINT4 idx=gsl_rng_uniform_int(runState->GSLrandom,Nlive);
	/* Copy the variable to avoid over-writing one of the live points */
	LALInferenceCopyVariables(livePoints[idx],&myCurrentParams);
	runState->currentParams=&myCurrentParams;
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

  /* Convert to a 2D array for ACF calculation */
  data_array=XLALCalloc(nPar,sizeof(REAL8 *));
  acf_array=XLALCalloc(nPar,sizeof(REAL8 *));
  for (i=0;i<(UINT4)nPar;i++){
    data_array[i]=XLALCalloc(max_iterations,sizeof(REAL8));
    acf_array[i]=XLALCalloc(max_iterations/2,sizeof(REAL8));
  }
  /* Measure autocorrelation in each dimension */
  /* Not ideal, should be measuring something like the det(autocorrelation-crosscorrelation matrix) */
  for (i=0;i<max_iterations;i++){
    for(j=0;j<nPar;j++) data_array[j][i]=*(REAL8 *)LALInferenceGetVariable(&variables_array[i],param_names[j]);
  }
  this=myCurrentParams.head;
  for(i=0;i<(UINT4)nPar;i++){
   /* Subtract the mean */
   REAL8 this_mean = gsl_stats_mean(&data_array[i][0], 1, max_iterations);
   for(j=0;j<max_iterations;j++) data_array[i][j]-=this_mean;
   ACL=1;
   int startflag=1;
   ACF=1.;
   /* Use GSL to compute the ACF */
   for(UINT4 lag=0;ACF>=ACF_TOLERANCE&&lag<max_iterations/2;lag++){
      ACF=(REAL8) gsl_stats_correlation(&data_array[i][0], 1, &data_array[i][lag], 1, max_iterations-lag);
      if(isnan(ACF)) ACF=1.;
      acf_array[i][lag]=ACF;
      ACL+=2.0*ACF;
      if((ACF<ACF_TOLERANCE && startflag) || lag==max_iterations/2-1){
	    startflag=0;
        ACL*=(REAL8)thinning;
	    if(ACL>max) max=ACL;
	    LALInferenceAddVariable(acls,param_names[i],&ACL,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
        break;
      }
   }
   if(LALInferenceCheckVariable(runState->algorithmParams,"verbose")) fprintf(stdout,"%s: mean= %lf, ACL=%lf\n",param_names[i],this_mean,ACL);
   do{this=this->next;}while(this && (this->vary==LALINFERENCE_PARAM_FIXED || this->vary==LALINFERENCE_PARAM_OUTPUT || this->type!=LALINFERENCE_REAL8_t));
  }
  if(acffile){
  /* Write out the ACF */
  for(j=0;j<(UINT4)nPar;j++) fprintf(acffile,"%s ",param_names[j]);
  fprintf(acffile,"\n");
  for(i=0;i<max_iterations/2;i++){
    for(j=0;j<(UINT4)nPar;j++) fprintf(acffile,"%f ",acf_array[j][i]);
    fprintf(acffile,"\n");
  }
  }

  /* Cache the samples */
  if(LALInferenceCheckVariable(oldAlgParams,"proposalcache"))
  {
    INT4 *Ncache=(INT4 *)LALInferenceGetVariable(oldAlgParams,"proposalcachesize");
    LALInferenceVariables **cache_ptr=(LALInferenceVariables **)LALInferenceGetVariable(oldAlgParams,"proposalcache");
    LALInferenceVariables *cache=*cache_ptr;
    UINT4 Nnew=max_iterations/(UINT4)(max/thinning);
    INT4 stride=max/thinning;
    if(LALInferenceCheckVariable(runState->algorithmParams,"verbose")) fprintf(stderr,"Caching %i samples\n",Nnew);

    /* Copy independent samples */
    REAL8 oldLogL=-DBL_MAX;
    for(i=stride,j=*Ncache;j<Nnew+*Ncache&&i<max_iterations;i+=stride,j++)
    {
      REAL8 newlogL=*(REAL8 *)LALInferenceGetVariable(&(variables_array[i]),"logL");
      if(newlogL==oldLogL) {j--; continue;}
      cache=XLALRealloc(cache,(j+1)*sizeof(LALInferenceVariables) );
      if(!cache) fprintf(stderr,"ERROR!!! Could not resize cache to %i!\n",j+1);
      memset(&(cache[j]),0,sizeof(LALInferenceVariables));
      LALInferenceCopyVariables(&(variables_array[i]),&(cache[j]));
      oldLogL=newlogL;
    }
    
    /* Update the state variables */
    *Ncache=j;
    *cache_ptr=cache;
  }
  
  /* Clean up */
  for(i=0;i<(UINT4)nPar;i++) {XLALFree(data_array[i]); XLALFree(acf_array[i]);}
  for (i=0;i<max_iterations;i++){
    LALInferenceClearVariables(&variables_array[i]);
  }
  XLALFree(variables_array);
  XLALFree(data_array); XLALFree(acf_array);
  LALInferenceClearVariables(&myAlgParams);
  LALInferenceClearVariables(&myCurrentParams);
  runState->currentParams=oldCurrentParams;
  runState->algorithmParams=oldAlgParams;
  if(chainfile) fclose(chainfile);
  if(acffile) fclose(acffile);
  XLALFree(param_names);
  return(acls);
}

static UINT4 checkForSingleAdapt(LALInferenceRunState *runState)
{
  UINT4 i = 0;
  LALInferenceProposalFunction *cycle = NULL;
  LALInferenceVariables *propArgs = runState->proposalArgs;
  cycle = *((LALInferenceProposalFunction **)LALInferenceGetVariable(propArgs, cycleArrayName));
  
  /* If there is not a proposal counter, put one into the variables, initialized to zero. */
  if (!LALInferenceCheckVariable(propArgs, cycleArrayCounterName)) {
    i = 0;
    LALInferenceAddVariable(propArgs, cycleArrayCounterName, &i, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_CIRCULAR);
  }
  
  i = *((UINT4 *)LALInferenceGetVariable(propArgs, cycleArrayCounterName));
  
  if (cycle[i]==&LALInferenceSingleAdaptProposal) return 1;
  else return 0;
  
}

/* Perform one MCMC iteration on runState->currentParams. Return 1 if accepted or 0 if not */
UINT4 LALInferenceMCMCSamplePrior(LALInferenceRunState *runState)
{
    UINT4 adaptProp=checkForSingleAdapt(runState);
    UINT4 outOfBounds=0;
    //LALInferenceVariables tempParams;
    REAL8 logProposalRatio=0.0;
    //LALInferenceVariables *oldParams=&tempParams;
    LALInferenceVariables proposedParams;
    memset(&proposedParams,0,sizeof(proposedParams));
    REAL8 logLmin=*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"logLmin");
    REAL8 thislogL;
    UINT4 accepted=0;

    REAL8 logPriorOld=*(REAL8 *)LALInferenceGetVariable(runState->currentParams,"logPrior");
    //LALInferenceCopyVariables(runState->currentParams,oldParams);
    LALInferenceCopyVariables(runState->currentParams,&proposedParams);
    
    if(adaptProp)
    {
      thislogL=runState->likelihood(runState->currentParams,runState->data,runState->model);
      if (logLmin<thislogL) outOfBounds=0;
    }
    
    logProposalRatio = runState->proposal(runState,runState->currentParams,&proposedParams);
    runState->proposal(runState,runState->currentParams,&proposedParams);
    REAL8 logPriorNew=runState->prior(runState, &proposedParams, runState->model);
    if(LALInferenceCheckVariable(runState->proposalArgs,"logProposalRatio"))
       logProposalRatio=*(REAL8 *)LALInferenceGetVariable(runState->proposalArgs,"logProposalRatio");
    if(logPriorNew==-DBL_MAX || isnan(logPriorNew) || log(gsl_rng_uniform(runState->GSLrandom)) > (logPriorNew-logPriorOld) + logProposalRatio) 
    {
	/* Reject - don't need to copy new params back to currentParams */
        /*LALInferenceCopyVariables(oldParams,runState->currentParams); */
    } 
    else {
        accepted=1;
	LALInferenceCopyVariables(&proposedParams,runState->currentParams);
        LALInferenceSetVariable(runState->currentParams,"logPrior",&logPriorNew);
    }
    LALInferenceClearVariables(&proposedParams);
    
    if((!outOfBounds)&&adaptProp)
    {
      thislogL=runState->likelihood(runState->currentParams,runState->data,runState->model);
      if(logLmin<thislogL) LALInferenceUpdateAdaptiveJumps(runState, accepted, 0.35);
      else LALInferenceUpdateAdaptiveJumps(runState, 0, 0.35);

    }
    
    LALInferenceTrackProposalAcceptance(runState, accepted);
    
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

/* Cache wrapper around another sampler */
INT4 LALInferenceNestedSamplingCachedSampler(LALInferenceRunState *runState)
{
  INT4 Naccept=0;
  if(!LALInferenceCheckVariable(runState->algorithmParams,"proposalcache") || !LALInferenceCheckVariable(runState->algorithmParams,"proposalcachesize"))
  {
    fprintf(stderr,"Adding cache variables in the sampler\n");
    /* Add space for the proposal cache */
    LALInferenceVariables *cache=NULL;
    INT4 newNcache=0;
    LALInferenceAddVariable(runState->algorithmParams,"proposalcache",(void *)&cache,LALINFERENCE_void_ptr_t,LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(runState->algorithmParams,"proposalcachesize",&newNcache,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_LINEAR);
  }
  
  INT4 *Ncache=(INT4 *)LALInferenceGetVariable(runState->algorithmParams,"proposalcachesize");
  LALInferenceVariables **cache_ptr=(LALInferenceVariables **)LALInferenceGetVariable(runState->algorithmParams,"proposalcache");
  LALInferenceVariables *cache=*cache_ptr;
  
  if(*Ncache==0 || LALInferenceGetProcParamVal(runState->commandLine,"--no-cache")){
    Naccept = LALInferenceNestedSamplingSloppySample(runState);
    return Naccept;
  }
  REAL8 logL=-DBL_MAX;
  REAL8 logLmin=*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"logLmin");
  /* Draw the last sample from the cache and reduce the size of the cache by one
   until we find one that has a high enough likelihood */
  do {
      LALInferenceVariables *new=&(cache[*Ncache-1]);
      logL=*(REAL8 *)LALInferenceGetVariable(new,"logL");
      if(logL>logLmin){
        runState->currentLikelihood=logL;
        LALInferenceCopyVariables(new,runState->currentParams);
      }
      LALInferenceClearVariables(new);
      new=NULL;
      cache=XLALRealloc(cache,sizeof(LALInferenceVariables) * (*Ncache-1));
      (*Ncache)--;
    } while (logL<=logLmin&&*Ncache>0);
  
  *cache_ptr=cache;
  /* If we didn't get any acceptable samples, call the main sampler */
  if(*Ncache==0 && logL<=logLmin)
  {
    Naccept = LALInferenceNestedSamplingSloppySample(runState);
  }

  return Naccept;
}

/* Sample the limited prior distribution using the MCMC method as usual, but
   only check the likelihood bound x fraction of the time. Always returns a fulled checked sample.
   x=LALInferenceGetVariable(runState->algorithmParams,"sloppyfraction")
   */

INT4 LALInferenceNestedSamplingSloppySample(LALInferenceRunState *runState)
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
    UINT4 sloppynumber=(UINT4) (sloppyfraction*(REAL8)Nmcmc);
    UINT4 testnumber=Nmcmc-sloppynumber;
    /* +1 for the last iteration which we do check */
    UINT4 subchain_length=(sloppynumber/testnumber) +1;
    REAL8 logLnew=0.0;
    UINT4 sub_iter=0;
    UINT4 tries=0;
    UINT4 ifo=0;
    REAL8 counter=1.;
    UINT4 BAILOUT=100*testnumber; /* If no acceptance after 100 tries, will exit and the sampler will try a different starting point */
    do{
        counter=counter-1.;
        subchain_length=0;
        /* Draw an independent sample from the prior */
        do{
            sub_accepted+=LALInferenceMCMCSamplePrior(runState);
            subchain_length++;
            counter+=(1.-sloppyfraction);
        }while(counter<1);
	/* Check that there was at least one accepted point */
	if(sub_accepted==0) {
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
        if(logLmin!=-DBL_MAX) logLnew=runState->likelihood(runState->currentParams,runState->data,runState->model);
        if(logLnew>logLmin || logLmin==-DBL_MAX) /* Accept */
        {
            Naccepted++;
            /* Update information to pass back out */
            LALInferenceAddVariable(runState->currentParams,"logL",(void *)&logLnew,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
            if(LALInferenceCheckVariable(runState->algorithmParams,"logZnoise")){
               tmp=logLnew-*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"logZnoise");
               LALInferenceAddVariable(runState->currentParams,"deltalogL",(void *)&tmp,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
            }
            ifo=0;
            while(data)
            {
               tmp=runState->model->ifo_loglikelihoods[ifo] - data->nullloglikelihood;
               sprintf(tmpName,"deltalogl%s",data->name);
               LALInferenceAddVariable(runState->currentParams,tmpName,&tmp,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
               ifo++;
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
    }while((mcmc_iter<testnumber||runState->currentLikelihood<=logLmin||Naccepted==0)&&(mcmc_iter<BAILOUT));
    /* Make sure likelihood is filled in if it wasn't done during sampling */
    if(logLnew==0.0){
            logLnew=runState->likelihood(runState->currentParams,runState->data,runState->model);
            runState->currentLikelihood=logLnew;
            LALInferenceAddVariable(runState->currentParams,"logL",(void *)&logLnew,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
            if(LALInferenceCheckVariable(runState->algorithmParams,"logZnoise")){
               tmp=logLnew-*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"logZnoise");
               LALInferenceAddVariable(runState->currentParams,"deltalogL",(void *)&tmp,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
            }
            ifo=0;
            while(data)
            {
               tmp=runState->model->ifo_loglikelihoods[ifo] - data->nullloglikelihood;
               sprintf(tmpName,"deltalogl%s",data->name);
               LALInferenceAddVariable(runState->currentParams,tmpName,&tmp,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
               ifo++;
               data=data->next;
            }
    }
    
    /* Compute some statistics for information */
    REAL8 sub_accept_rate=(REAL8)sub_accepted/(REAL8)sub_iter;
    REAL8 accept_rate=(REAL8)Naccepted/(REAL8)testnumber;
    LALInferenceSetVariable(runState->algorithmParams,"accept_rate",&accept_rate);
    LALInferenceSetVariable(runState->algorithmParams,"sub_accept_rate",&sub_accept_rate);
    /* Adapt the sloppy fraction toward target acceptance of outer chain */
    if(logLmin!=-DBL_MAX){
        if((REAL8)accept_rate>Target) { sloppyfraction+=5.0/(REAL8)Nmcmc;}
        else { sloppyfraction-=5.0/(REAL8)Nmcmc;}
        if(sloppyfraction>maxsloppyfraction) sloppyfraction=maxsloppyfraction;
	if(sloppyfraction<minsloppyfraction) sloppyfraction=minsloppyfraction;
	
	LALInferenceSetVariable(runState->algorithmParams,"sloppyfraction",&sloppyfraction);
    }
    /* Cleanup */
    LALInferenceClearVariables(&oldParams);
    return Naccepted;
}


/* Evolve nested sampling algorithm by one step, i.e.
 evolve runState->currentParams to a new point with higher
 likelihood than currentLikelihood. Uses the MCMC method with sloppy sampling.
 */
INT4 LALInferenceNestedSamplingOneStep(LALInferenceRunState *runState)
{
     INT4 Naccept = LALInferenceNestedSamplingCachedSampler(runState);
     return Naccept;
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
	runState->livePoints=XLALCalloc(Nlive,sizeof(LALInferenceVariables *));
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
	LALInferenceVariables *curParsBackup=runState->currentParams;
	for(i=0;i<Nlive;i++)
	{
	   /* Clone the currentParams into LivePoints[i] */
	    runState->livePoints[i]=XLALCalloc(1,sizeof(LALInferenceVariables));
	    /* Copy the param structure */
	    LALInferenceCopyVariables(runState->currentParams,runState->livePoints[i]);
	  
	  /* Sprinkle the varying points among prior */
	  do{
	    LALInferenceDrawFromPrior( runState->livePoints[i], runState->priorArgs, runState->GSLrandom );
	    logPrior=runState->prior(runState,runState->livePoints[i],runState->model);
	  }while(logPrior==-DBL_MAX || isnan(logPrior));
	  /* Populate log likelihood */
	  logLs->data[i]=runState->likelihood(runState->livePoints[i],runState->data,runState->model);
	  LALInferenceAddVariable(runState->livePoints[i],"logL",(void *)&(logLs->data[i]),LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
	  LALInferenceAddVariable(runState->livePoints[i],"logPrior",(void*)&logPrior,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
	}
	runState->currentParams=curParsBackup;
	if(!runState->currentParams) runState->currentParams=XLALCalloc(1,sizeof(LALInferenceVariables));
	
}


void LALInferenceSetupkDTreeNSLivePoints( LALInferenceRunState *runState ){
  /* create a k-d tree from the nested sampling live points */
  LALInferenceKDTree *tree;
  REAL8 *low = NULL, *high = NULL; /* upper and lower bounds of tree */
  size_t ndim = 0;
  LALInferenceVariableItem *currentItem;
  UINT4 cnt = 0;
  REAL8 *pt = NULL;
  LALInferenceVariables *templt =
    XLALCalloc(1,sizeof(LALInferenceVariables));
  UINT4 Nlive = *(UINT4 *)LALInferenceGetVariable( runState->algorithmParams,
                                                   "Nlive" );
 
  /* if a current tree exists remove it */
  if ( LALInferenceCheckVariable( runState->proposalArgs, "kDTree" ) )
  {
    LALInferenceKDTreeDelete( *(LALInferenceKDTree
      **)LALInferenceGetVariable(runState->proposalArgs, "kDTree"));
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
  pt = XLALMalloc(cnt*sizeof(REAL8));
  
  /* set up tree */
  tree = LALInferenceKDEmpty( low, high, ndim );
  LALInferenceCopyVariables( runState->currentParams, templt );
                    
  /* add points to tree */
  for( cnt = 0; cnt < Nlive; cnt++ ){
    LALInferenceKDVariablesToREAL8( runState->livePoints[cnt], pt, templt );
    
    LALInferenceKDAddPoint( tree, pt );
  }
                  
  /* add tree */
  LALInferenceAddVariable( runState->proposalArgs, "kDTree", &tree,
                           LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_FIXED );
  
  /* if template doesn't exist add it */
  if ( !LALInferenceCheckVariable( runState->proposalArgs,
                                   "kDTreeVariableTemplate" ) ){
    LALInferenceAddVariable( runState->proposalArgs, "kDTreeVariableTemplate",
                             &templt, LALINFERENCE_void_ptr_t,
                             LALINFERENCE_PARAM_FIXED );
  }
  
  XLALFree( high );
  XLALFree( low );
  XLALFree( pt );
}

static void SetupEigenProposals(LALInferenceRunState *runState)
{ 
  gsl_matrix *eVectors=NULL;
  gsl_vector *eValues =NULL;
  REAL8Vector *eigenValues=NULL;
  /* Check for existing covariance matrix */
  gsl_matrix **cvm=NULL;
  if(LALInferenceCheckVariable(runState->proposalArgs,"covarianceMatrix"))
    LALInferenceRemoveVariable(runState->proposalArgs,"covarianceMatrix");
//    cvm=(gsl_matrix **)LALInferenceGetVariable(runState->proposalArgs,"covarianceMatrix");
  cvm=XLALCalloc(1,sizeof(gsl_matrix *));
  
  /* Add the covariance matrix for proposal distribution */
  UINT4 *Nlive=LALInferenceGetVariable(runState->algorithmParams,"Nlive");
  /* Sort the variables to ensure consistent order */
  for(UINT4 i=0;i<*Nlive;i++) LALInferenceSortVariablesByName(runState->livePoints[i]);
  LALInferenceNScalcCVM(cvm,runState->livePoints,*Nlive);
  UINT4 N=(*cvm)->size1;


  /* Check for the eigenvectors and values */
  if(LALInferenceCheckVariable(runState->proposalArgs,"covarianceEigenvectors"))
    eVectors=*(gsl_matrix **)LALInferenceGetVariable(runState->proposalArgs,"covarianceEigenvectors");
  else
    eVectors=gsl_matrix_alloc(N,N);
  
  if(LALInferenceCheckVariable(runState->proposalArgs,"covarianceEigenvalues"))
    eigenValues=*(REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs,"covarianceEigenvalues");
  else
    eigenValues=XLALCreateREAL8Vector(N);

  /* Set up eigenvectors and eigenvalues. */
  gsl_matrix *covCopy = gsl_matrix_alloc(N,N);
  eVectors = gsl_matrix_alloc(N,N);
  eValues = gsl_vector_alloc(N);
  gsl_eigen_symmv_workspace *ws = gsl_eigen_symmv_alloc(N);
  int gsl_status;
  gsl_matrix_memcpy(covCopy, *cvm);
  
  if ((gsl_status = gsl_eigen_symmv(covCopy, eValues, eVectors, ws)) != GSL_SUCCESS) {
    XLALPrintError("Error in gsl_eigen_symmv (in %s, line %d): %d: %s\n", __FILE__, __LINE__, gsl_status, gsl_strerror(gsl_status));
    XLAL_ERROR_VOID(XLAL_EFAILED);
  }
  
  for (UINT4 i = 0; i < N; i++) {
    eigenValues->data[i] = gsl_vector_get(eValues,i);
  }
  
  LALInferenceAddVariable(runState->proposalArgs, "covarianceEigenvectors", &eVectors, LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(runState->proposalArgs, "covarianceEigenvalues", &eigenValues, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(runState->proposalArgs,"covarianceMatrix",cvm,LALINFERENCE_gslMatrix_t,LALINFERENCE_PARAM_OUTPUT);
  
  gsl_matrix_free(covCopy);
  gsl_vector_free(eValues);
  gsl_eigen_symmv_free(ws);
  XLALFree(cvm);
}

static int _saveNSintegralState(FILE *fp, NSintegralState *s);
static int _saveNSintegralState(FILE *fp, NSintegralState *s)
{
  UINT4 N=s->size;
  if(1!=fwrite(&N,sizeof(UINT4),1,fp)) return 1;
  if(1!=fwrite(&(s->iteration),sizeof(s->iteration),1,fp)) return 1;
  if(N!=fwrite(s->logZarray,sizeof(REAL8),N,fp)) return 1;
  if(N!=fwrite(s->oldZarray,sizeof(REAL8),N,fp)) return 1;
  if(N!=fwrite(s->Harray,sizeof(REAL8),N,fp)) return 1;
  if(N!=fwrite(s->logwarray,sizeof(REAL8),N,fp)) return 1;
  if(N!=fwrite(s->logtarray,sizeof(REAL8),N,fp)) return 1;
  if(N!=fwrite(s->logt2array,sizeof(REAL8),N,fp)) return 1;
  return 0;
}
static int _loadNSintegralState(FILE *fp, NSintegralState *s);
static int _loadNSintegralState(FILE *fp, NSintegralState *s)
{
  if(1!=fread(& (s->size) , sizeof(UINT4), 1, fp)) return 1;
  UINT4 N=s->size;
  if(1!=fread(&(s->iteration),sizeof(UINT4),1,fp)) return 1;
  if(N!=fread(s->logZarray,sizeof(REAL8),N,fp)) return 1;
  if(N!=fread(s->oldZarray,sizeof(REAL8),N,fp)) return 1;
  if(N!=fread(s->Harray,sizeof(REAL8),N,fp)) return 1;
  if(N!=fread(s->logwarray,sizeof(REAL8),N,fp)) return 1;
  if(N!=fread(s->logtarray,sizeof(REAL8),N,fp)) return 1;
  if(N!=fread(s->logt2array,sizeof(REAL8),N,fp)) return 1;
  return 0;
}

static int WriteNSCheckPoint(CHAR *filename, LALInferenceRunState *runState, NSintegralState *s)
{
  FILE *progfile=fopen(filename,"w");
  if(!progfile)
  {
    fprintf(stderr,"Unable to save resume file %s!\n",filename);
    return 1;
  }
  else
  {
    UINT4 Nlive=*(UINT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nlive");
    int retcode= _saveNSintegralState(progfile,s);
    if(retcode) {
        fprintf(stderr,"Unable to write nested sampling state - will not be able to resume!\n");
        fclose(progfile);
        return 1;
    }
    LALInferenceWriteVariablesArrayBinary(progfile,runState->livePoints, Nlive);
    fclose(progfile);
    return 0;
  }
}

static int ReadNSCheckPoint(CHAR *filename, LALInferenceRunState *runState, NSintegralState *s)
{
  FILE *progfile=fopen(filename,"r");
  if(!progfile)
  {
    fprintf(stderr,"Unable to load resume file %s!\n",filename);
    return 1;
  }
  else
  {
    UINT4 Nlive=*(UINT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nlive");
    int retcode=_loadNSintegralState(progfile,s);
    if(retcode){
        fprintf(stderr,"Unable to read nested sampling state - unable to resume!\n");
        fclose(progfile);
        return 1;
    }
    //if(retcode) return 1;
    LALInferenceReadVariablesArrayBinary(progfile,runState->livePoints,Nlive);
    fclose(progfile);
    return 0;
  }
}
