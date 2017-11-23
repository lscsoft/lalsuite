/* Implementation of Nested Sampling for LALInference.
 * (C) John Veitch, 2010
 */

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <signal.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <lal/TimeDelay.h>
#include <lal/LALInferenceConfig.h>
#include <lal/LALStdlib.h>
#include <lal/LALInferenceNestedSampler.h>
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceProposal.h>
#include <lal/LALInferenceHDF5.h>
#include <lal/LALInferencePriorVolumes.h>

#include "logaddexp.h"

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
  REAL8Vector *logZarray,
	*oldZarray,
	*Harray,
	*logwarray,
	*logtarray,
	*logt2array;
} NSintegralState;

static struct itimerval checkpoint_timer;

/** Utility functions for the resume functionality */
/** Write the current state to a checkpoint file the given filename */
/** Read the given filename to populate a given LALInferenceRunState and NSintegralState */
#ifdef HAVE_HDF5
static int _saveNSintegralStateH5(LALH5File *group, NSintegralState *s);
static int _saveNSintegralStateH5(LALH5File *group, NSintegralState *s)
{
  XLALH5FileAddScalarAttribute(group, "iteration", &(s->iteration), LAL_U4_TYPE_CODE );
  XLALH5FileWriteREAL8Vector(group, "logZarray", s->logZarray);
  XLALH5FileWriteREAL8Vector(group, "oldZarray", s->oldZarray);
  XLALH5FileWriteREAL8Vector(group, "Harray", s->Harray);
  XLALH5FileWriteREAL8Vector(group, "logwarray", s->logwarray);
  XLALH5FileWriteREAL8Vector(group, "logtarray", s->logtarray);
  XLALH5FileWriteREAL8Vector(group, "logt2array", s->logt2array);
  return(0);
}

static int _loadNSintegralStateH5(LALH5File *group, NSintegralState *s);
static int _loadNSintegralStateH5(LALH5File *group, NSintegralState *s)
{
  XLALH5FileQueryScalarAttributeValue(&(s->iteration), group, "iteration");
  s->logZarray = XLALH5FileReadREAL8Vector(group, "logZarray");
  s->oldZarray = XLALH5FileReadREAL8Vector(group, "oldZarray");
  s->Harray = XLALH5FileReadREAL8Vector(group, "Harray");
  s->logwarray = XLALH5FileReadREAL8Vector(group, "logwarray");
  s->logtarray = XLALH5FileReadREAL8Vector(group, "logtarray");
  s->logt2array = XLALH5FileReadREAL8Vector(group, "logt2array");
  s->size = s->logZarray->length;
  return(0);
}


static int ReadNSCheckPointH5(char *filename, LALInferenceRunState *runState, NSintegralState *s);
static int WriteNSCheckPointH5(char *filename, LALInferenceRunState *runState, NSintegralState *s);

static int WriteNSCheckPointH5(char *filename, LALInferenceRunState *runState, NSintegralState *s)
{
  LALH5File *h5file = XLALH5FileOpen(filename,"w");
  if(!h5file)
  {
    fprintf(stderr,"Unable to save resume file %s!\n",filename);
    return(1);
  }
  UINT4 Nlive=*(UINT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nlive");
  LALH5File *group = XLALH5GroupOpen(h5file, "lalinferencenest_checkpoint");
  if(!group) XLAL_ERROR(XLAL_EFAILED,"Unable to read group lalinferencenest_checkpoint\n");
  int retcode = _saveNSintegralStateH5(group,s);
  if(retcode) XLAL_ERROR(XLAL_EFAILED,"Unable to save integral state\n");
  LALInferenceH5VariablesArrayToDataset(group, runState->livePoints, Nlive, "live_points");
  INT4 N_output_array=0;
  if(LALInferenceCheckVariable(runState->algorithmParams,"N_outputarray")) N_output_array=LALInferenceGetINT4Variable(runState->algorithmParams,"N_outputarray");
  if(N_output_array>0)
  {
    LALInferenceVariables **output_array=NULL;
    output_array=*(LALInferenceVariables ***)LALInferenceGetVariable(runState->algorithmParams,"outputarray");
    LALInferenceH5VariablesArrayToDataset(group, output_array, N_output_array, "past_chain");
    XLALH5FileAddScalarAttribute(group, "N_outputarray", &N_output_array, LAL_I4_TYPE_CODE);
  }
  
  XLALH5FileClose(group);
  XLALH5FileClose(h5file);
  return(retcode);
}

static int ReadNSCheckPointH5(char *filename, LALInferenceRunState *runState, NSintegralState *s)
{
  int retcode;
  LALH5File *h5file;
  UINT4 Nlive;
  if( access( filename, F_OK ) == -1 ) return(1);
  XLAL_TRY(h5file = XLALH5FileOpen(filename,"r"),retcode);
  if(retcode!=XLAL_SUCCESS) return(retcode);
  LALH5File *group;
  XLAL_TRY(group = XLALH5GroupOpen(h5file,"lalinferencenest_checkpoint"),retcode);
  if(retcode!=XLAL_SUCCESS) return(retcode);
  UINT4 N_outputarray;
  LALInferenceVariables **outputarray;
  XLALH5FileQueryScalarAttributeValue(&N_outputarray, group, "N_outputarray");
  if(!h5file)
  {
    fprintf(stderr,"Unable to load resume file %s!\n",filename);
    return(1);
  }
  printf("restoring nested sampling integral state\n");
  retcode=_loadNSintegralStateH5(group,s);
  if(retcode){
    fprintf(stderr,"Unable to read nested sampling state - unable to resume!\n");
    XLALH5FileClose(h5file);
    return 1;
  }
  LALH5Dataset *liveGroup = XLALH5DatasetRead(group,"live_points");
  retcode = LALInferenceH5DatasetToVariablesArray(liveGroup , &(runState->livePoints), &Nlive );
  printf("restored %i live points\n",Nlive);
  XLALH5DatasetFree(liveGroup);
  if(N_outputarray>0)
  {
    printf("restoring %i past iterations\n",N_outputarray);
    LALH5Dataset *outputGroup = XLALH5DatasetRead(group, "past_chain");
    retcode |= LALInferenceH5DatasetToVariablesArray(outputGroup, &outputarray, &N_outputarray);
    LALInferenceAddVariable(runState->algorithmParams,"N_outputarray",&N_outputarray,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddVariable(runState->algorithmParams,"outputarray",&outputarray,LALINFERENCE_void_ptr_t,LALINFERENCE_PARAM_OUTPUT);
    XLALH5DatasetFree(outputGroup);
  }
  XLALH5FileClose(group);
  XLALH5FileClose(h5file);
  printf("done restoring\n");
  return(retcode);
  
}

#else
static int ReadNSCheckPoint(CHAR *filename, LALInferenceRunState *runState, NSintegralState *s);
static int WriteNSCheckPoint(CHAR *filename, LALInferenceRunState *runState, NSintegralState *s);

static int _saveNSintegralState(FILE *fp, NSintegralState *s);
static int _saveNSintegralState(FILE *fp, NSintegralState *s)
{
  UINT4 N=s->size;
  if(1!=fwrite(&N,sizeof(UINT4),1,fp)) return 1;
  if(1!=fwrite(&(s->iteration),sizeof(s->iteration),1,fp)) return 1;
  if(N!=fwrite(s->logZarray->data,sizeof(REAL8),N,fp)) return 1;
  if(N!=fwrite(s->oldZarray->data,sizeof(REAL8),N,fp)) return 1;
  if(N!=fwrite(s->Harray->data,sizeof(REAL8),N,fp)) return 1;
  if(N!=fwrite(s->logwarray->data,sizeof(REAL8),N,fp)) return 1;
  if(N!=fwrite(s->logtarray->data,sizeof(REAL8),N,fp)) return 1;
  if(N!=fwrite(s->logt2array->data,sizeof(REAL8),N,fp)) return 1;
  return 0;
}
static int _loadNSintegralState(FILE *fp, NSintegralState *s);
static int _loadNSintegralState(FILE *fp, NSintegralState *s)
{
  if(1!=fread(& (s->size) , sizeof(UINT4), 1, fp)) return 1;
  UINT4 N=s->size;
  s->logZarray = XLALCreateREAL8Vector(N);
  s->oldZarray = XLALCreateREAL8Vector(N);
  s->Harray = XLALCreateREAL8Vector(N);
  s->logwarray = XLALCreateREAL8Vector(N);
  s->logtarray = XLALCreateREAL8Vector(N);
  s->logt2array = XLALCreateREAL8Vector(N);
  if(1!=fread(&(s->iteration),sizeof(UINT4),1,fp)) return 1;
  if(N!=fread(s->logZarray->data,sizeof(REAL8),N,fp)) return 1;
  if(N!=fread(s->oldZarray->data,sizeof(REAL8),N,fp)) return 1;
  if(N!=fread(s->Harray->data,sizeof(REAL8),N,fp)) return 1;
  if(N!=fread(s->logwarray->data,sizeof(REAL8),N,fp)) return 1;
  if(N!=fread(s->logtarray->data,sizeof(REAL8),N,fp)) return 1;
  if(N!=fread(s->logt2array->data,sizeof(REAL8),N,fp)) return 1;
  return 0;
}

static int WriteNSCheckPoint(CHAR *filename, LALInferenceRunState *runState, NSintegralState *s)
{
  FILE *progfile=fopen(filename,"w");
  int errnum;
  if(!progfile)
  {
    fprintf(stderr,"Unable to save resume file %s!\n",filename);
    return 1;
  }
  else
  {
    if(setvbuf(progfile,NULL,_IOFBF,0x100000)) /* Set buffer to 1MB so as to not thrash NFS */
      fprintf(stderr,"Warning: Unable to set resume file buffer!");
    UINT4 Nlive=*(UINT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nlive");
    int retcode= _saveNSintegralState(progfile,s);
    if(retcode) {
      fprintf(stderr,"Unable to write nested sampling state - will not be able to resume!\n");
      fclose(progfile);
      return 1;
    }
    XLAL_TRY(LALInferenceWriteVariablesArrayBinary(progfile,runState->livePoints, Nlive),errnum);
	if(errnum!=XLAL_SUCCESS)
	{
      fprintf(stderr,"Unable to write live points - will not be able to resume!\n");
      fclose(progfile);
      return 1;
	}
    INT4 N_output_array=0;
    if(LALInferenceCheckVariable(runState->algorithmParams,"N_outputarray")) N_output_array=LALInferenceGetINT4Variable(runState->algorithmParams,"N_outputarray");
    fwrite(&N_output_array,sizeof(INT4),1,progfile);
    if(N_output_array!=0 )
    {
      LALInferenceVariables **output_array=NULL;
      output_array=*(LALInferenceVariables ***)LALInferenceGetVariable(runState->algorithmParams,"outputarray");
      XLAL_TRY(LALInferenceWriteVariablesArrayBinary(progfile,output_array, N_output_array),errnum);
	  if(errnum!=XLAL_SUCCESS)
	  {
		      fprintf(stderr,"Unable to write past chain - will not be able to resume!\n");
		      fclose(progfile);
		      return 1;	
	  }
      fprintf(stderr,"Resume --> wrote %d past chain samples\n\n",N_output_array);
    }
    fclose(progfile);
    return 0;
  }
}

static int ReadNSCheckPoint(CHAR *filename, LALInferenceRunState *runState, NSintegralState *s)
{
  int errnum;
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
    XLAL_TRY(LALInferenceReadVariablesArrayBinary(progfile,runState->livePoints,Nlive),errnum);
	if(errnum!=XLAL_SUCCESS)
    {
		fprintf(stderr,"Unable to read live points from resume file!\n");
		fclose(progfile);
		return(1);
	}
    INT4 N_output_array;
    fread(&N_output_array,sizeof(INT4),1,progfile);
    LALInferenceVariables **output_array=NULL;
    if(N_output_array!=0){
      output_array=XLALCalloc(N_output_array,sizeof(LALInferenceVariables *));
      fprintf(stderr,"Resume --> read %d past chain samples\n",N_output_array);
      XLAL_TRY(LALInferenceReadVariablesArrayBinary(progfile,output_array,N_output_array),errnum);
	  if(errnum!=XLAL_SUCCESS)
	  {
		fprintf(stderr,"Unable to read past chain from resume file!\n");
		fclose(progfile);
		return(1);
	  }
      if(LALInferenceCheckVariable(runState->algorithmParams,"N_outputarray")) LALInferenceRemoveVariable(runState->algorithmParams,"N_outputarray");
      if(LALInferenceCheckVariable(runState->algorithmParams,"outputarray")) LALInferenceRemoveVariable(runState->algorithmParams,"outputarray");
      LALInferenceAddVariable(runState->algorithmParams,"N_outputarray",&N_output_array,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_OUTPUT);
      LALInferenceAddVariable(runState->algorithmParams,"outputarray",&output_array,LALINFERENCE_void_ptr_t,LALINFERENCE_PARAM_OUTPUT);
    }
    fclose(progfile);
    return 0;
  }
}

#endif

/** Sync the live points to the differential evolution buffer */
static int syncLivePointsDifferentialPoints(LALInferenceRunState *state, LALInferenceThreadState *thread);

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

static REAL8 LALInferenceNSSample_logt(int Nlive,gsl_rng *RNG);

//static void SamplePriorDiscardAcceptance(LALInferenceRunState *runState);
static REAL8 mean(REAL8 *array,int N);


/** Calculate covariance matrix from a collection of live points */
static void LALInferenceNScalcCVM(gsl_matrix **cvm, LALInferenceVariables **Live, UINT4 Nlive);

/* log( exp(a) - exp(b) ) */
static double logsubexp(double a, double b);
static double logsubexp(double a, double b)
{
		if(b>a)
		{
				fprintf(stderr,"Cannot take log of negative number %lf - %lf = %lf !\n",exp(a),exp(b),exp(a)-exp(b));
				return -INFINITY;
		}
		else
		{
				return a + log1p(-exp(b-a));
		}
}

static void SetupEigenProposals(LALInferenceRunState *runState);

/**
 * Update the internal state of the integrator after receiving the lowest logL
 * value logL
 */
static REAL8 incrementEvidenceSamples(gsl_rng *GSLrandom, UINT4 Nlive, REAL8 logL, NSintegralState *s);
static REAL8 incrementEvidenceSamples(gsl_rng *GSLrandom, UINT4 Nlive, REAL8 logL, NSintegralState *s)
{
  REAL8 Wtarray[s->size];
  /* Update evidence array */
  for(UINT4 j=0;j<s->size;j++){
    s->oldZarray->data[j]=s->logZarray->data[j];
    Wtarray[j]=s->logwarray->data[j]+logL+logsubexp(0,s->logt2array->data[j] + s->logtarray->data[j])-log(2.0);
    s->logZarray->data[j]=logaddexp(s->logZarray->data[j],Wtarray[j]);
    if(isfinite(s->oldZarray->data[j]) && isfinite(s->logZarray->data[j]))
        s->Harray->data[j]= exp(Wtarray[j]-s->logZarray->data[j])*logL
          + exp(s->oldZarray->data[j]-s->logZarray->data[j])*(s->Harray->data[j]+s->oldZarray->data[j])-s->logZarray->data[j];
    REAL8 tmp=s->logtarray->data[j];
    s->logtarray->data[j]=s->logt2array->data[j];
    s->logt2array->data[j]=tmp;
    s->logtarray->data[j]=LALInferenceNSSample_logt(Nlive,GSLrandom);
    s->logwarray->data[j]+=s->logtarray->data[j];
  }
  s->iteration++;
  return(mean(s->logZarray->data,s->logZarray->length));
}

static void printAdaptiveJumpSizes(FILE *file, LALInferenceThreadState *threadState);
static void printAdaptiveJumpSizes(FILE *file, LALInferenceThreadState *threadState)
{
    LALInferenceVariableItem *this=threadState->currentParams->head;
    REAL8 *val=NULL;
    char tmpname[1000]="";
    int first=1;
	while(this)
    {
        sprintf(tmpname,"%s_%s",this->name,ADAPTSUFFIX);
        if(LALInferenceCheckVariable(threadState->proposalArgs,tmpname))
        {
		    if(first) {fprintf(file,"Adaptive proposal step sizes:\n"); first=0;}
            val=(REAL8 *)LALInferenceGetVariable(threadState->proposalArgs,tmpname);
            fprintf(file,"%s: %lf\n",this->name,*val);
        }
        this=this->next;
    }
}


/* Create Internal arrays for sampling the integral */
static NSintegralState *initNSintegralState(UINT4 Nruns, UINT4 Nlive);
static NSintegralState *initNSintegralState(UINT4 Nruns, UINT4 Nlive)
{
  NSintegralState *s=XLALMalloc(sizeof(NSintegralState));
  s->iteration=0;
  s->size=Nruns;
  s->logZarray = XLALCreateREAL8Vector(Nruns);
  s->oldZarray = XLALCreateREAL8Vector(Nruns);
  s->Harray = XLALCreateREAL8Vector(Nruns);
  s->logwarray = XLALCreateREAL8Vector(Nruns);
  s->logtarray=XLALCreateREAL8Vector(Nruns);
  s->logt2array=XLALCreateREAL8Vector(Nruns);
  REAL8 logw=0;

  if(s->logZarray==NULL || s->Harray==NULL || s->oldZarray==NULL || s->logwarray==NULL)
  {fprintf(stderr,"Unable to allocate RAM\n"); exit(-1);}
  for(UINT4 i=0;i<Nruns;i++)  {s->logwarray->data[i]=logw; s->logZarray->data[i]=-INFINITY;
    s->oldZarray->data[i]=-INFINITY; s->Harray->data[i]=0.0;s->logtarray->data[i]=-1.0/Nlive;
    s->logt2array->data[i]=-1.0/Nlive;
  }
  return s;
}

static REAL8 mean(REAL8 *array,int N){
	REAL8 sum=0.0;
	int i;
	for(i=0;i<N;i++) sum+=array[i];
	return sum/((REAL8) N);
}


static REAL8 LALInferenceNSSample_logt(int Nlive,gsl_rng *RNG){
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
        if(LALInferenceCheckVariable(runState->algorithmParams,"Nmcmc")){ /* if already estimated the length */
          if ( LALInferenceGetINT4Variable(runState->algorithmParams,"Nmcmc") != 0 ){
            max=4 * *(INT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nmcmc"); /* We will use this to go out 4x last ACL */
          }
          else { max=4*maxMCMC; } /* otherwise use the MAX_MCMC */
        }
        else max=4*maxMCMC; /* otherwise use the MAX_MCMC */
        if(max>4*maxMCMC) max=4*maxMCMC;
	UINT4 Nlive = *(INT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nlive");
	UINT4 rnd=gsl_rng_uniform_int(runState->GSLrandom,Nlive);
      /* Single threaded here */
	LALInferenceCopyVariables(runState->livePoints[rnd],runState->threads[0]->currentParams);
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
    /* Single threaded here */
    if (LALInferenceGetProcParamVal(runState->commandLine,"--proposal-kde"))
        LALInferenceSetupClusteredKDEProposalFromDEBuffer(runState->threads[0]);
    return(max);
}

/* estimateCovarianceMatrix reads the list of live points,
 and works out the covariance matrix of the varying parameters
 - CIRCULAR parameters are wrapped around before the calculation and must be
 scaled to the range 0 -> 2pi. Only works with REAL8 variables */
static void LALInferenceNScalcCVM(gsl_matrix **cvm, LALInferenceVariables **Live, UINT4 Nlive)
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

void LALInferenceNestedSamplingAlgorithmInit(LALInferenceRunState *runState)
{
  char help[]="\
    ----------------------------------------------\n\
    --- Nested Sampling Algorithm Parameters -----\n\
    ----------------------------------------------\n\
    --Nlive N                        Number of live points to use\n\
    (--Nmcmc M)                      Over-ride auto chain length determination and use <M> MCMC samples\n\
    (--maxmcmc M)                    Use at most M MCMC points when autodetermining the chain (5000)\n\
    (--Nmcmcinitial M)               Use M MCMC points when initially resampling from the prior\n\
                                     (otherwise default is to use maxmcmc)\n\
    (--sloppyratio S)                Number of sub-samples of the prior for every sample from the\n\
                                     limited prior\n\
    (--Nruns R)                      Number of parallel samples from logt to use(1)\n\
    (--tolerance dZ)                 Tolerance of nested sampling algorithm (0.1)\n\
    (--randomseed seed)              Random seed of sampling distribution\n\
    (--prior )                       Set the prior to use (InspiralNormalised,SkyLoc,malmquist)\n\
                                     (default: InspiralNormalised)\n\
    (--sampleprior N)                For Testing: Draw N samples from the prior, will not perform the\n\
                                     nested sampling integral\n\
    (--progress)                     Output some progress information at each iteration\n\
    (--verbose)                      Output more info. N=1: errors, N=2 (default): warnings, N=3: info\n\
    (--resume)                       Allow non-condor checkpointing every 4 hours. If given will check \n\
                                     or OUTFILE_resume and continue if possible\n\
    \n";

  ProcessParamsTable *ppt=NULL;
  /* Print command line arguments if help requested */
  if(runState == NULL || LALInferenceGetProcParamVal(runState->commandLine,"--help"))
  {
    fprintf(stdout,"%s",help);
    return;
  }
  ProcessParamsTable *commandLine=runState->commandLine;

  INT4 verbose=0;
  INT4 x=0;
  ppt=LALInferenceGetProcParamVal(commandLine,"--verbose");
  if(ppt) {
    if(ppt->value[0]){
      x=atoi(ppt->value);
      switch(x){
        case 0:
          verbose=LALNDEBUG; /* Nothing */
          break;
        case 1:
          verbose=LALMSGLVL1; /* Only errors */
          break;
        case 2:
          verbose=LALMSGLVL2; /* Errors and warnings */
          break;
        case 3:
          verbose=LALMSGLVL3; /* Errors, warnings and info */
          break;
        default:
          verbose=LALMSGLVL2;
          break;
      }
    }
    else verbose=LALMSGLVL2; /* Errors and warnings */
    LALInferenceAddVariable(runState->algorithmParams,"verbose", &verbose , LALINFERENCE_INT4_t,
                            LALINFERENCE_PARAM_FIXED);
  }
  INT4 tmpi=0;
  REAL8 tmp=0;

  /* Single thread only */
  LALInferenceThreadState *threadState = runState->threads[0];

  /* Set up the appropriate functions for the nested sampling algorithm */
  runState->algorithm=&LALInferenceNestedSamplingAlgorithm;
  runState->evolve=&LALInferenceNestedSamplingOneStep;

  /* use the ptmcmc proposal to sample prior */
  threadState->proposal=&LALInferenceCyclicProposal;
  REAL8 temp=1.0;
  LALInferenceAddVariable(runState->proposalArgs,"temperature",&temp,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);

  runState->logsample=LALInferenceLogSampleToArray;

  /* Number of live points */
  ppt=LALInferenceGetProcParamVal(commandLine,"--Nlive");
  if (!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--nlive");
  if(ppt)
    tmpi=atoi(ppt->value);
  else {
    fprintf(stderr,"Error, must specify number of live points\n");
    exit(1);
  }
  LALInferenceAddVariable(runState->algorithmParams,"Nlive",&tmpi, LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);

  /* Number of points in MCMC chain */
  ppt=LALInferenceGetProcParamVal(commandLine,"--Nmcmc");
  if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--nmcmc");
  if(ppt){
    tmpi=atoi(ppt->value);
    LALInferenceAddVariable(runState->algorithmParams,"Nmcmc",&tmpi,
                            LALINFERENCE_INT4_t,LALINFERENCE_PARAM_OUTPUT);
    printf("set number of MCMC points, over-riding auto-determination!\n");
  }

  /* Maximum number of points in MCMC chain */
  ppt=LALInferenceGetProcParamVal(commandLine,"--maxmcmc");
  if(ppt){
    tmpi=atoi(ppt->value);
    LALInferenceAddVariable(runState->algorithmParams,"maxmcmc",&tmpi,
                            LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
  }

  /* Set fraction for sloppy sampling */
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--sloppyfraction")))
    tmp=atof(ppt->value);
  else tmp=0.0;
  LALInferenceAddVariable(runState->algorithmParams,"sloppyfraction",&tmp,
                          LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);

  /* Optionally specify number of parallel runs */
  ppt=LALInferenceGetProcParamVal(commandLine,"--Nruns");
  if(ppt) {
    tmpi=atoi(ppt->value);
    LALInferenceAddVariable(runState->algorithmParams,"Nruns",&tmpi,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
  }

  printf("set tolerance.\n");
  /* Tolerance of the Nested sampling integrator */
  ppt=LALInferenceGetProcParamVal(commandLine,"--tolerance");
  if(ppt){
    tmp=strtod(ppt->value,(char **)NULL);
    LALInferenceAddVariable(runState->algorithmParams,"tolerance",&tmp, LALINFERENCE_REAL8_t,
                            LALINFERENCE_PARAM_FIXED);
  }
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
  /* Single thread here */
  LALInferenceThreadState *threadState = runState->threads[0];
  UINT4 HDFOUTPUT=1;
  UINT4 Nlive=*(UINT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nlive");
  UINT4 Nruns=100;
  REAL8 *logZarray,*Harray,*logwarray,*logtarray;
  REAL8 TOLERANCE=0.1;
  REAL8 logZ,logZnew,logLmin,logLmax=-INFINITY,logLtmp,logw,H,logZnoise,dZ=0;
  LALInferenceVariables *temp;
  FILE *fpout=NULL;
  REAL8 neginfty=-INFINITY;
  REAL8 zero=0.0;
  REAL8 *logLikelihoods=NULL;
  UINT4 verbose=0;
  REAL8 sloppyfrac;
  UINT4 displayprogress=0;
  LALInferenceVariables *currentVars=XLALCalloc(1,sizeof(LALInferenceVariables));
  UINT4 samplePrior=0; //If this flag is set to a positive integer, code will just draw this many samples from the prior
  ProcessParamsTable *ppt=NULL;

  if(!runState->logsample) runState->logsample=LALInferenceLogSampleToArray;

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
    LALInferenceAddVariable(runState->algorithmParams,"logLmin",&neginfty,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);

  if(!LALInferenceCheckVariable(runState->algorithmParams,"accept_rate"))
    LALInferenceAddVariable(runState->algorithmParams,"accept_rate",&zero,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
  if(!LALInferenceCheckVariable(runState->algorithmParams,"sub_accept_rate"))
    LALInferenceAddVariable(runState->algorithmParams,"sub_accept_rate",&zero,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
  if(!LALInferenceCheckVariable(runState->algorithmParams,"sloppyfraction"))
    LALInferenceAddVariable(runState->algorithmParams,"sloppyfraction",&zero,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);

  /* Open output file */
  ppt = NULL;
  ppt=LALInferenceGetProcParamVal(runState->commandLine,"--outfile");
  if(!ppt){
      fprintf(stderr,"Must specify --outfile <filename.hdf5>\n");
      exit(1);
  }
  char *outfile=ppt->value;
  /* Check if the output file has hdf5 extension */
  if(strstr(outfile,".h5") || strstr(outfile,".hdf")) HDFOUTPUT=1;
  else HDFOUTPUT=0;
  
#ifndef HAVE_HDF5
  if(HDFOUTPUT)
  {
      fprintf(stderr,"Error: LALSuite was compiled without HDF5 support. Unable to write HDF5 output files\n");
      exit(1);
  }
#endif
  
  double logvolume=0.0;
  if ( LALInferenceCheckVariable( runState->livePoints[0], "chirpmass" ) ){
    /* If a cbc run, calculate the mass-distance volume and store it to file*/ 
    /* Do it before algorithm starts so that we can kill the run and still get this */
    logvolume=log(LALInferenceMassDistancePriorVolume(runState));
  }

  if(LALInferenceGetProcParamVal(runState->commandLine,"--progress"))
    displayprogress=1;

  minpos=0;

  logw=log(1.0-exp(-1.0/Nlive));
  i=0;
  /* sort points for consistent order before creating the matrix */
  for(i=0;i<Nlive;i++)
    LALInferenceSortVariablesByName(runState->livePoints[i]);

  /* Set up eigenvector proposals */
  SetupEigenProposals(runState);

  /* Use the live points as differential evolution points */
  /* Single thread here */
  syncLivePointsDifferentialPoints(runState,threadState);
  threadState->differentialPointsSkip=1;

  if(!LALInferenceCheckVariable(runState->algorithmParams,"Nmcmc")){
    INT4 tmp=MAX_MCMC;
    if(LALInferenceGetProcParamVal(runState->commandLine,"--Nmcmcinitial")){
      tmp=atoi(LALInferenceGetProcParamVal(runState->commandLine,"--Nmcmcinitial")->value);
    }
    LALInferenceAddVariable(runState->algorithmParams,"Nmcmc",&tmp,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_OUTPUT);
  }
  s=initNSintegralState(Nruns,Nlive);

  /* Check for an interrupted run */
  char resumefilename[FILENAME_MAX+10];
  snprintf(resumefilename,sizeof(resumefilename),"%s_resume",outfile);
  int retcode=1;
  if(LALInferenceGetProcParamVal(runState->commandLine,"--resume")){
#ifdef HAVE_HDF5
      retcode=ReadNSCheckPointH5(resumefilename,runState,s);
#else
      retcode=ReadNSCheckPoint(resumefilename,runState,s);
#endif
      if(retcode==0){
          for(i=0;i<Nlive;i++) logLikelihoods[i]=*(REAL8 *)LALInferenceGetVariable(runState->livePoints[i],"logL");
          iter=s->iteration;
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
      checkpoint_timer.it_interval.tv_sec=30*60; /* Default timer 30 mins */
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
    LALInferenceSetVariable(runState->algorithmParams,"logLmin",&neginfty);
    int sprinklewarning=0;
    for(i=0;i<Nlive;i++) {
	    threadState->currentParams=runState->livePoints[i];
	    j=0;
	    do
	    {
		    runState->evolve(runState);
		    logLikelihoods[i]=runState->likelihood(runState->livePoints[i],runState->data,threadState->model);
		    if(!sprinklewarning && j++==100) { fprintf(stderr,"Warning: having difficulty sampling prior, check your prior bounds\n"); sprinklewarning=1;}
	    }while(isnan(logLikelihoods[i]) || isinf(logLikelihoods[i]));
	    if(XLALPrintProgressBar((double)i/(double)Nlive)) fprintf(stderr,"\n");
    }
  }

  logZarray=s->logZarray->data;
  logtarray=s->logtarray->data;
  logwarray=s->logwarray->data;
  Harray=s->Harray->data;

  /* Find maximum likelihood and sanity check */
  for(i=0;i<Nlive;i++)
  {
    LALInferenceAddVariable(runState->livePoints[i],"logw",&logw,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
    logLtmp=logLikelihoods[i];
    logLmax=logLtmp>logLmax? logLtmp : logLmax;
    if(isnan(logLikelihoods[i]) || isinf(logLikelihoods[i])) {
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

  /* Reset proposal stats before starting */
  LALInferenceZeroProposalStats(threadState->cycle);

  /* Set the number of MCMC points */
  UpdateNMCMC(runState);
  /* Output some information */
  if(verbose){
    LALInferencePrintProposalStatsHeader(stdout,threadState->cycle);
    LALInferencePrintProposalStats(stdout,threadState->cycle);
    LALInferenceZeroProposalStats(threadState->cycle);
    printAdaptiveJumpSizes(stdout, threadState);
  }
  /* Write out names of parameters */
  if(!HDFOUTPUT)
  {
      FILE *lout=NULL;
      char param_list[FILENAME_MAX+16];
      snprintf(param_list,sizeof(param_list),"%s_params.txt",outfile);
      lout=fopen(param_list,"w");
      LALInferenceFprintParameterHeaders(lout,runState->livePoints[0]);
      fclose(lout);
  }
  minpos=0;
  threadState->currentParams=currentVars;
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
    if(samplePrior) logLmin=-INFINITY;

    logZnew=incrementEvidenceSamples(runState->GSLrandom, Nlive, logLikelihoods[minpos], s);
    //deltaZ=logZnew-logZ; - set but not used
    H=mean(Harray,Nruns);
    logZ=logZnew;
    if(runState->logsample) runState->logsample(runState->algorithmParams,runState->livePoints[minpos]);
    UINT4 itercounter=0;

    /* Generate a new live point */
    do{ /* This loop is here in case it is necessary to find a different sample */
      /* Clone an old live point and evolve it */
      while((j=gsl_rng_uniform_int(runState->GSLrandom,Nlive))==minpos){};
      LALInferenceCopyVariables(runState->livePoints[j],threadState->currentParams);
      threadState->currentLikelihood = logLikelihoods[j];
      LALInferenceSetVariable(runState->algorithmParams,"logLmin",(void *)&logLmin);
      runState->evolve(runState);
      itercounter++;
    }while( threadState->currentLikelihood<=logLmin ||  *(REAL8*)LALInferenceGetVariable(runState->algorithmParams,"accept_rate")==0.0);

    LALInferenceCopyVariables(threadState->currentParams,runState->livePoints[minpos]);
    logLikelihoods[minpos]=threadState->currentLikelihood;

  if (threadState->currentLikelihood>logLmax)
    logLmax=threadState->currentLikelihood;

  logw=mean(logwarray,Nruns);
  LALInferenceAddVariable(runState->livePoints[minpos],"logw",&logw,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
  dZ=logaddexp(logZ,logLmax-((double) iter)/((double)Nlive))-logZ;
  sloppyfrac=*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"sloppyfraction");
  if(displayprogress) fprintf(stderr,"%i: accpt: %1.3f Nmcmc: %i sub_accpt: %1.3f slpy: %2.1f%% H: %3.2lf nats logL:%.3lf ->%.3lf logZ: %.3lf deltalogLmax: %.2lf dZ: %.3lf Zratio: %.3lf \n",\
    iter,\
    *(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"accept_rate")/(REAL8)itercounter,\
    *(INT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nmcmc"),\
    *(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"sub_accept_rate"),\
    100.0*sloppyfrac,\
    H,\
    logLmin,\
    threadState->currentLikelihood,\
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
#ifdef HAVE_HDF5

      WriteNSCheckPointH5(resumefilename,runState,s);
#else
      WriteNSCheckPoint(resumefilename,runState,s);
#endif
      fflush(fpout);
      __ns_saveStateFlag=0;
    }
   /* Have we been told to quit? */
  if(__ns_exitFlag) {
    exit(0);
  }

  /* Update the proposal */
  if(!(iter%(Nlive/10))) {
    /* Update the covariance matrix */
    //WriteNSCheckPointH5(resumefilename, runState, s);
    if ( LALInferenceCheckVariable( threadState->proposalArgs,"covarianceMatrix" ) ){
      SetupEigenProposals(runState);
    }

    /* Update NMCMC from ACF */
    UpdateNMCMC(runState);

    /* Sync the live points to differential points */
    syncLivePointsDifferentialPoints(runState,threadState);

    /* Output some information */
    if(verbose){
      LALInferencePrintProposalStatsHeader(stdout,threadState->cycle);
      LALInferencePrintProposalStats(stdout,threadState->cycle);
      LALInferenceZeroProposalStats(threadState->cycle);
      printAdaptiveJumpSizes(stdout, threadState);
    }
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
    logZ=incrementEvidenceSamples(runState->GSLrandom, Nlive-i, logLikelihoods[i], s);
    if(runState->logsample) runState->logsample(runState->algorithmParams,runState->livePoints[i]);
  }
  
    LALInferenceVariables **output_array=NULL;
    UINT4 N_output_array=0;
    if(LALInferenceCheckVariable(runState->algorithmParams,"outputarray")
            &&LALInferenceCheckVariable(runState->algorithmParams,"N_outputarray") )
    {
            output_array=*(LALInferenceVariables ***)LALInferenceGetVariable(runState->algorithmParams,"outputarray");
            N_output_array=*(UINT4 *)LALInferenceGetVariable(runState->algorithmParams,"N_outputarray");
    }
    double logB=logZ-logZnoise;
    /* Pass output back through algorithmparams */
    LALInferenceAddVariable(runState->algorithmParams,"logZ",(void *)&logZ,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddVariable(runState->algorithmParams,"logB",(void *)&logB,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddVariable(runState->algorithmParams,"logLmax",(void *)&logLmax,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);

    /* Write out the evidence */
    if(!HDFOUTPUT)
    {
        fpout = fopen(outfile,"w");
        for(i=0;i<N_output_array;i++)
		{
				LALInferencePrintSample(fpout,output_array[i]);
				fprintf(fpout,"\n");
		}
        fclose(fpout);
    
        char bayesfile[FILENAME_MAX+10];
        sprintf(bayesfile,"%s_B.txt",outfile);
        fpout=fopen(bayesfile,"w");
        fprintf(fpout,"%lf %lf %lf %lf\n",logZ-logZnoise,logZ,logZnoise,logLmax);
        fclose(fpout);
    
        
     
    }
    /* Write HDF5 file */
    if(HDFOUTPUT)
    {
      
      LALH5File *h5file=XLALH5FileOpen(outfile, "w");
      // Create group heirarchy 
      char runID[2048];
      if((ppt=LALInferenceGetProcParamVal(runState->commandLine,"--runid")))
        snprintf(runID,sizeof(runID),"%s_%s","lalinference_nest",ppt->value);
      else
        snprintf(runID,sizeof(runID),"lalinference_nest");

      LALH5File *groupPtr = LALInferenceH5CreateGroupStructure(h5file, "lalinference", runID);
      /* Create run identifier group */
      LALInferenceH5VariablesArrayToDataset(groupPtr, output_array, N_output_array, LALInferenceHDF5NestedSamplesDatasetName);
      /* TODO: Write metadata */
      XLALH5FileAddScalarAttribute(groupPtr, "log_evidence", &logZ, LAL_D_TYPE_CODE);
      XLALH5FileAddScalarAttribute(groupPtr, "log_bayes_factor", &logB, LAL_D_TYPE_CODE);
      XLALH5FileAddScalarAttribute(groupPtr, "information_nats", &H, LAL_D_TYPE_CODE);
      XLALH5FileAddScalarAttribute(groupPtr, "log_noise_evidence", &logZnoise, LAL_D_TYPE_CODE );
      XLALH5FileAddScalarAttribute(groupPtr, "log_max_likelihood", &logLmax , LAL_D_TYPE_CODE);
      XLALH5FileAddScalarAttribute(groupPtr, "number_live_points", &Nlive, LAL_U4_TYPE_CODE);
      XLALH5FileAddScalarAttribute(groupPtr, "log_prior_volume", &logvolume, LAL_D_TYPE_CODE);

      XLALH5FileClose(h5file);
    }
  
    if(output_array) {
      for(i=0;i<N_output_array;i++){
        LALInferenceClearVariables(output_array[i]);
        XLALFree(output_array[i]);
      }
      XLALFree(output_array);
    }
  
  /* Clean up resume file */
  if(LALInferenceGetProcParamVal(runState->commandLine,"--resume"))
  {
    if(!access(resumefilename,W_OK)) remove(resumefilename);
  }

  /* Free memory */
  XLALFree(logtarray); XLALFree(logwarray); XLALFree(logZarray);
}

/* Calculate the autocorrelation function of the sampler (runState->evolve) for each parameter
 * Evolves the sample starting with the value passed in temp, with a maximum of max_iterations steps.
 Return the ACL for each parameter as a LALInferenceVariables */
LALInferenceVariables *LALInferenceComputeAutoCorrelation(LALInferenceRunState *runState, UINT4 max_iterations, LALInferenceEvolveOneStepFunction evolve)
{
  /* Single threaded here */
  LALInferenceThreadState *threadState = runState->threads[0];
  ProcessParamsTable *ppt=NULL;
  char chainfilename[2048]="";
  char acf_file_name[2048]="";
  FILE *chainfile=NULL;
  FILE *acffile=NULL;
  UINT4 i,j;
  UINT4 nPar=0; // = LALInferenceGetVariableDimensionNonFixed(runState->currentParams);
  REAL8 **data_array=NULL;
  REAL8 **acf_array=NULL;
  LALInferenceVariableItem *this;
  INT4 thinning=1;
  max_iterations/=thinning;
  /* Find the number and names of variables */
  for(this=threadState->currentParams->head;this;this=this->next) if(this->vary!=LALINFERENCE_PARAM_FIXED && this->vary!=LALINFERENCE_PARAM_OUTPUT && this->type==LALINFERENCE_REAL8_t) nPar++;
  char **param_names=XLALCalloc(nPar,sizeof(char *));
  for(i=0,this=threadState->currentParams->head;this;this=this->next) if(this->vary!=LALINFERENCE_PARAM_FIXED && this->vary!=LALINFERENCE_PARAM_OUTPUT && this->type==LALINFERENCE_REAL8_t) param_names[i++]=this->name;

  REAL8 ACF,ACL,max=0;
  LALInferenceVariables *acls=XLALCalloc(1,sizeof(LALInferenceVariables));

  /* Back up the algorithm state and replace with a clean version for logSampletoarray */
  LALInferenceVariables myAlgParams,*oldAlgParams=runState->algorithmParams;
  LALInferenceVariables myCurrentParams,*oldCurrentParams=threadState->currentParams;
  memset(&myAlgParams,0,sizeof(LALInferenceVariables));
  memset(&myCurrentParams,0,sizeof(LALInferenceVariables));
  LALInferenceCopyVariables(oldAlgParams,&myAlgParams);
  if(LALInferenceCheckVariable(&myAlgParams,"outputarray")) LALInferenceRemoveVariable(&myAlgParams,"outputarray");
  if(LALInferenceCheckVariable(&myAlgParams,"N_outputarray")) LALInferenceRemoveVariable(&myAlgParams,"N_outputarray");
  if(LALInferenceCheckVariable(&myAlgParams,"outfile")) LALInferenceRemoveVariable(&myAlgParams,"outfile");
  LALInferenceRemoveVariable(&myAlgParams,"Nmcmc");
  LALInferenceAddVariable(&myAlgParams,"Nmcmc",&thinning,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_OUTPUT);

  LALInferenceSortVariablesByName(&myCurrentParams);
  runState->algorithmParams=&myAlgParams;
  threadState->currentParams=&myCurrentParams;
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
	threadState->currentParams=&myCurrentParams;
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
  LALInferenceLogSampleToArray(runState->algorithmParams,threadState->currentParams);
  /* Evolve the initial sample (i starts at 1)*/
  for(i=1;i<max_iterations;i++)
  {
   evolve(runState);
   LALInferenceLogSampleToArray(runState->algorithmParams,threadState->currentParams);
  }

  /* Get the location of the sample array */
  LALInferenceVariables **variables_array=*(LALInferenceVariables ***)LALInferenceGetVariable(runState->algorithmParams,"outputarray");

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
    for(j=0;j<nPar;j++) data_array[j][i]=*(REAL8 *)LALInferenceGetVariable(variables_array[i],param_names[j]);
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
    REAL8 oldLogL=-INFINITY;
    for(i=stride,j=*Ncache;j<Nnew+*Ncache&&i<max_iterations;i+=stride,j++)
    {
      REAL8 newlogL=*(REAL8 *)LALInferenceGetVariable(variables_array[i],"logL");
      if(newlogL==oldLogL) {j--; continue;}
      cache=XLALRealloc(cache,(j+1)*sizeof(LALInferenceVariables) );
      if(!cache) fprintf(stderr,"ERROR!!! Could not resize cache to %i!\n",j+1);
      memset(&(cache[j]),0,sizeof(LALInferenceVariables));
      LALInferenceCopyVariables(variables_array[i],&(cache[j]));
      oldLogL=newlogL;
    }

    /* Update the state variables */
    *Ncache=j;
    *cache_ptr=cache;
  }

  /* Clean up */
  for(i=0;i<(UINT4)nPar;i++) {XLALFree(data_array[i]); XLALFree(acf_array[i]);}
  for (i=0;i<max_iterations;i++){
    LALInferenceClearVariables(variables_array[i]);
    XLALFree(variables_array[i]);
  }
  XLALFree(variables_array);
  XLALFree(data_array); XLALFree(acf_array);
  LALInferenceClearVariables(&myAlgParams);
  LALInferenceClearVariables(&myCurrentParams);
  threadState->currentParams=oldCurrentParams;
  runState->algorithmParams=oldAlgParams;
  if(chainfile) fclose(chainfile);
  if(acffile) fclose(acffile);
  XLALFree(param_names);
  return(acls);
}

/* Perform one MCMC iteration on runState->currentParams. Return 1 if accepted or 0 if not */
UINT4 LALInferenceMCMCSamplePrior(LALInferenceRunState *runState)
{
    /* Single threaded here */
    LALInferenceThreadState * threadState=runState->threads[0];
    UINT4 outOfBounds=0;
    UINT4 adaptProp=0;
    //LALInferenceVariables tempParams;
    REAL8 logProposalRatio=0.0;
    //LALInferenceVariables *oldParams=&tempParams;
    LALInferenceVariables proposedParams;
    memset(&proposedParams,0,sizeof(proposedParams));
    REAL8 logLmin=*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"logLmin");
    REAL8 thislogL=-INFINITY;
    UINT4 accepted=0;

    REAL8 logPriorOld=*(REAL8 *)LALInferenceGetVariable(threadState->currentParams,"logPrior");
    if(adaptProp)
    {
          thislogL=runState->likelihood(threadState->currentParams,runState->data,threadState->model);
          if (logLmin<thislogL) outOfBounds=0;
    }
    LALInferenceCopyVariables(threadState->currentParams,&proposedParams);

    logProposalRatio = threadState->proposal(threadState,threadState->currentParams,&proposedParams);
    REAL8 logPriorNew=runState->prior(runState, &proposedParams, threadState->model);
    if(isinf(logPriorNew) || isnan(logPriorNew) || log(gsl_rng_uniform(runState->GSLrandom)) > (logPriorNew-logPriorOld) + logProposalRatio)
    {
	/* Reject - don't need to copy new params back to currentParams */
        /*LALInferenceCopyVariables(oldParams,runState->currentParams); */
    }
    else {
        accepted=1;
        //printf("Accepted line %i\n",__LINE__);
        LALInferenceCopyVariables(&proposedParams,threadState->currentParams);
        LALInferenceSetVariable(threadState->currentParams,"logPrior",&logPriorNew);
    }
    LALInferenceClearVariables(&proposedParams);

    if((!outOfBounds) && adaptProp)
    {
      thislogL=runState->likelihood(threadState->currentParams,runState->data,threadState->model);
      if(logLmin<thislogL) threadState->accepted = accepted;
      else threadState->accepted=accepted=0;
      LALInferenceUpdateAdaptiveJumps(threadState, 0.35);
    }
    //printf("logLnew = %lf, logPriorNew = %lf, logProposalRatio = %lf\n",thislogL,logPriorNew,logProposalRatio);
    threadState->accepted=accepted;
    LALInferenceTrackProposalAcceptance(threadState);
    //printf("Accepted = %i\n",accepted);
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
  /* Single thread here */
  LALInferenceThreadState *threadState = runState->threads[0];
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
  REAL8 logL=-INFINITY;
  REAL8 logLmin=*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"logLmin");

  /* Draw the last sample from the cache and reduce the size of the cache by one
   until we find one that has a high enough likelihood */
  do {
      LALInferenceVariables *new=&(cache[*Ncache-1]);
      logL=*(REAL8 *)LALInferenceGetVariable(new,"logL");
      if(logL>logLmin){
        threadState->currentLikelihood=logL;
        LALInferenceCopyVariables(new,threadState->currentParams);
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
    /* Single thread here */
    LALInferenceThreadState *threadState = runState->threads[0];
    LALInferenceIFOData *data=runState->data;
    REAL8 tmp;
    REAL8 Target=0.3;
    char tmpName[320];
    REAL8 logLold=*(REAL8 *)LALInferenceGetVariable(threadState->currentParams,"logL");
    memset(&oldParams,0,sizeof(oldParams));
    LALInferenceCopyVariables(threadState->currentParams,&oldParams);
    REAL8 logLmin=*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"logLmin");
    UINT4 Nmcmc=*(UINT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nmcmc");
    REAL8 maxsloppyfraction=((REAL8)Nmcmc-1)/(REAL8)Nmcmc ;
    REAL8 sloppyfraction=maxsloppyfraction/2.0;
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
    const char *extra_names[]={"logL","optimal_snr","matched_filter_snr","deltalogL"}; /* Names for parameters to be stripped when sampling prior */
    UINT4 Nnames = 4;
    if ( Nmcmc ){
      do{
        counter=counter-1.;
        subchain_length=0;
        for(UINT4 i=0;i<Nnames;i++)
        {
          if(LALInferenceCheckVariable(threadState->currentParams,extra_names[i]))
            LALInferenceRemoveVariable(threadState->currentParams,extra_names[i]);
        }
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
          LALInferenceCopyVariables(&oldParams,threadState->currentParams);
          threadState->currentLikelihood=logLold;
          continue;
        }
        tries=0;
        mcmc_iter++;
    	sub_iter+=subchain_length;
        if(isfinite(logLmin)) logLnew=runState->likelihood(threadState->currentParams,runState->data,threadState->model);
        if(logLnew>logLmin || isinf(logLmin)) /* Accept */
        {
            Naccepted++;
            /* Update information to pass back out */
            LALInferenceAddVariable(threadState->currentParams,"logL",(void *)&logLnew,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
            if(LALInferenceCheckVariable(runState->algorithmParams,"logZnoise")){
               tmp=logLnew-*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"logZnoise");
               LALInferenceAddVariable(threadState->currentParams,"deltalogL",(void *)&tmp,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
            }
            ifo=0;
            data=runState->data;
            while(data)
            {
               if(!threadState->model->ifo_loglikelihoods) break;
               tmp=threadState->model->ifo_loglikelihoods[ifo] - data->nullloglikelihood;
               sprintf(tmpName,"deltalogl%s",data->name);
               LALInferenceAddVariable(threadState->currentParams,tmpName,&tmp,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
               ifo++;
               data=data->next;
            }
            LALInferenceCopyVariables(threadState->currentParams,&oldParams);
            logLold=logLnew;
            threadState->currentLikelihood=logLnew;
        }
        else /* reject */
        {
            LALInferenceCopyVariables(&oldParams,threadState->currentParams);
            threadState->currentLikelihood=logLold;
        }
      }while((mcmc_iter<testnumber||threadState->currentLikelihood<=logLmin||Naccepted==0)&&(mcmc_iter<BAILOUT));
    }
    /* Make sure likelihood is filled in if it wasn't done during sampling */
    if(logLnew==0.0){
            logLnew=runState->likelihood(threadState->currentParams,runState->data,threadState->model);
            threadState->currentLikelihood=logLnew;
            LALInferenceAddVariable(threadState->currentParams,"logL",(void *)&logLnew,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
            if(LALInferenceCheckVariable(runState->algorithmParams,"logZnoise")){
               tmp=logLnew-*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"logZnoise");
               LALInferenceAddVariable(threadState->currentParams,"deltalogL",(void *)&tmp,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
            }
            ifo=0;
            data=runState->data;
            while(data && threadState->model->ifo_loglikelihoods)
            {
              tmp=threadState->model->ifo_loglikelihoods[ifo] - data->nullloglikelihood;
              sprintf(tmpName,"deltalogl%s",data->name);
              LALInferenceAddVariable(threadState->currentParams,tmpName,&tmp,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
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
    if(isfinite(logLmin)){
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

    /* Single thread here */
    LALInferenceThreadState *threadState = runState->threads[0];

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
	threadState->differentialPoints=runState->livePoints;
	threadState->differentialPointsLength=(size_t) Nlive;
	logLs=XLALCreateREAL8Vector(Nlive);

	LALInferenceAddVariable(runState->algorithmParams,"logLikelihoods",&logLs,LALINFERENCE_REAL8Vector_t,LALINFERENCE_PARAM_FIXED);
	fprintf(stdout,"Sprinkling %i live points, may take some time\n",Nlive);
	LALInferenceVariables *curParsBackup=threadState->currentParams;
	for(i=0;i<Nlive;i++)
	{
	   /* Clone the currentParams into LivePoints[i] */
	    runState->livePoints[i]=XLALCalloc(1,sizeof(LALInferenceVariables));
	    /* Copy the param structure */
	    LALInferenceCopyVariables(threadState->currentParams,runState->livePoints[i]);

	  /* Sprinkle the varying points among prior */
	  do{
	    LALInferenceDrawFromPrior( runState->livePoints[i], runState->priorArgs, runState->GSLrandom );
	    logPrior=runState->prior(runState,runState->livePoints[i],threadState->model);
	  }while(isinf(logPrior) || isnan(logPrior));
	  /* Populate log likelihood */
	  logLs->data[i]=runState->likelihood(runState->livePoints[i],runState->data,threadState->model);
	  LALInferenceAddVariable(runState->livePoints[i],"logL",(void *)&(logLs->data[i]),LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
	  LALInferenceAddVariable(runState->livePoints[i],"logPrior",(void*)&logPrior,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
	}
	threadState->currentParams=curParsBackup;
	if(!threadState->currentParams) threadState->currentParams=XLALCalloc(1,sizeof(LALInferenceVariables));
}


static void SetupEigenProposals(LALInferenceRunState *runState)
{
  /* Single thread here */
  LALInferenceThreadState *threadState = runState->threads[0];
  gsl_matrix *eVectors=NULL;
  gsl_vector *eValues =NULL;
  REAL8Vector *eigenValues=NULL;
  /* Check for existing covariance matrix */
  gsl_matrix **cvm=NULL;
  if(LALInferenceCheckVariable(threadState->proposalArgs,"covarianceMatrix"))
    LALInferenceRemoveVariable(threadState->proposalArgs,"covarianceMatrix");
//    cvm=(gsl_matrix **)LALInferenceGetVariable(threadState->proposalArgs,"covarianceMatrix");
  cvm=XLALCalloc(1,sizeof(gsl_matrix *));

  /* Add the covariance matrix for proposal distribution */
  UINT4 *Nlive=LALInferenceGetVariable(runState->algorithmParams,"Nlive");
  /* Sort the variables to ensure consistent order */
  for(UINT4 i=0;i<*Nlive;i++) LALInferenceSortVariablesByName(runState->livePoints[i]);
  LALInferenceNScalcCVM(cvm,runState->livePoints,*Nlive);
  UINT4 N=(*cvm)->size1;


  /* Check for the eigenvectors and values */
  if(LALInferenceCheckVariable(threadState->proposalArgs,"covarianceEigenvectors"))
    eVectors=*(gsl_matrix **)LALInferenceGetVariable(threadState->proposalArgs,"covarianceEigenvectors");
  else
    eVectors=gsl_matrix_alloc(N,N);

  if(LALInferenceCheckVariable(threadState->proposalArgs,"covarianceEigenvalues"))
    eigenValues=*(REAL8Vector **)LALInferenceGetVariable(threadState->proposalArgs,"covarianceEigenvalues");
  else
    eigenValues=XLALCreateREAL8Vector(N);

  /* Set up eigenvectors and eigenvalues. */
  gsl_matrix *covCopy = gsl_matrix_alloc(N,N);
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

  if(!LALInferenceCheckVariable(threadState->proposalArgs,"covarianceEigenvectors"))
    LALInferenceAddVariable(threadState->proposalArgs, "covarianceEigenvectors", &eVectors, LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_FIXED);
  if(!LALInferenceCheckVariable(threadState->proposalArgs,"covarianceEigenvalues"))
    LALInferenceAddVariable(threadState->proposalArgs, "covarianceEigenvalues", &eigenValues, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(threadState->proposalArgs,"covarianceMatrix",cvm,LALINFERENCE_gslMatrix_t,LALINFERENCE_PARAM_OUTPUT);

  gsl_matrix_free(covCopy);
  gsl_vector_free(eValues);
  gsl_eigen_symmv_free(ws);
  XLALFree(cvm);
}


static int syncLivePointsDifferentialPoints(LALInferenceRunState *state, LALInferenceThreadState *thread)
{
    INT4 N = LALInferenceGetINT4Variable(state->algorithmParams,"Nlive");
    if(!thread->differentialPoints) thread->differentialPoints=XLALCalloc(N,sizeof(LALInferenceVariables *));

    for(INT4 i=0;i<N;i++)
    {
        if(!thread->differentialPoints[i]) thread->differentialPoints[i]=XLALCalloc(1,sizeof(LALInferenceVariables));
        LALInferenceCopyVariables(state->livePoints[i],thread->differentialPoints[i]);
    }
    thread->differentialPointsLength=N;
    return(XLAL_SUCCESS);
}
