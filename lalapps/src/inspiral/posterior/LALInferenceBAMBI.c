/*
 *  InferenceBAMBI.c:  BAMBI with LALInference
 *
 *  Copyright (C) 2009 Ilya Mandel, Vivien Raymond, Christian Roever, Marc van der Sluys, John Veitch, Farhan Feroz, and Philip Graff
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
#include <lal/Date.h>
#include <lal/GenerateInspiral.h>
#include <lal/LALInference.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/StringInput.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/TimeSeries.h>
#include <lalapps.h>
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceReadData.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALInferenceProposal.h>
#include <lal/LALInferenceInit.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#ifdef __INTEL_COMPILER             // if the MultiNest library was compiled with ifor
       #define NESTRUN nested_mp_nestrun_
#elif defined __GNUC__                 // if the MultiNest library was compiled with gfortran
       #define NESTRUN __nested_MOD_nestrun
#else
       #error Do not know how to link to Fortran libraries, check symbol table for your platform (nm libnest3.a | grep nestrun) & edit example_eggbox_C++/eggbox.cc
#endif

LALInferenceRunState *initialize(ProcessParamsTable *commandLine);
void initializeMN(LALInferenceRunState *runState);
void initStudentt(LALInferenceRunState *state);

/******** Defined for BAMBI (start) **********/

#define BAMBI_STRLEN 200

extern float *omicron,tol,thL[3],logLRange;
extern double *maxsigma,logZero;
extern int nNN,nNNerr,totpar,loglcalls,ncheck,myid,nproc;
extern char *root,*networkinputs;
extern bool likenetinit,converged,lastconverged,netres,firstrun,discardpts;
extern int ignoredbambicalls,counter;
extern size_t nlayers,nnodes[10];
extern int doBAMBI,useNN,whitenin,whitenout,resume;

void BAMBIRun(int mmodal, int ceff, int nlive, double tol, double efr, int ndims, int nPar, int nClsPar,
    int maxModes, int updInt, double Ztol, char *root, int seed, int *pWrap, int fb, int resume, int outfile,
    int initMPI, double logZero, int maxiter, void (*LogLike)(double *, int *, int *, double *, void *),
    void (*dumper)(int *, int *, int *, double **, double **, double **, double *, double *, double *, void *),
    int (*bambi)(int *, int *, double **, double *), void *context);

extern void NESTRUN(int *mmodal, int *ceff, int *nlive, double *tol, double *efr, int *ndims, int *nPar, int *nClsPar,
    int *maxModes, int *updInt, double *Ztol, char *root, int *seed, int *pWrap, int *fb, int *resume, int *outfile,
    int *initMPI, double *logZero, int *maxiter, void (*LogLike)(double *, int *, int *, double *, void *),
    void (*dumper)(int *, int *, int *, double **, double **, double **, double *, double *, double *, void *),
    int (*bambi)(int *, int *, double **, double *), void *context);

void LogLike(double *Cube, int *ndim, int *npars, double *lnew, void *context);
int bambi(int *ndata, int *ndim, double **BAMBIData, double *lowlike);
void dumper(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *logZerr, void *context);
void LALInferenceMultiNestAlgorithm(LALInferenceRunState *runState);

/******** Defined for BAMBI (end) **********/

LALInferenceRunState *runStateGlobal;

void BAMBIRun(int mmodal, int ceff, int nlive, double tol2, double efr, int ndims, int nPar, int nClsPar,
int maxModes, int updInt, double Ztol, char *root2, int seed, int *pWrap, int fb, int resume2, int outfile,
int initMPI, double logZero2, int maxiter, void (*LogLike2)(double *, int *, int *, double *, void *),
void (*dumper2)(int *, int *, int *, double **, double **, double **, double *, double *, double *, void *),
int (*bambi2)(int *, int *, double **, double *), void *context)
{
    int i;
    char rootformn[BAMBI_STRLEN];
    strcpy(rootformn, root2);
    for (i = strlen(rootformn); i < BAMBI_STRLEN; i++) rootformn[i] = ' ';

    // Do "nm libbambi.a | grep nestrun" to find the name of the function to put here.
    // Remove one leading underscore for the name.

    NESTRUN(&mmodal, &ceff, &nlive, &tol2, &efr, &ndims, &nPar, &nClsPar, &maxModes, &updInt, &Ztol,
        rootformn, &seed, pWrap, &fb, &resume2, &outfile, &initMPI, &logZero2, &maxiter, LogLike2, dumper2, bambi2, context);
}

void getLogLike(double *Cube, UNUSED int *ndim, UNUSED int *npars, double *lnew, void *context);
void getphysparams(double *Cube, UNUSED int *ndim, UNUSED int *nPar, void *context);
void getallparams(double *Cube, UNUSED int *ndim, UNUSED int *nPar, void *context);

void getLogLike(double *Cube, UNUSED int *ndim, UNUSED int *npars, double *lnew, void *context)
{
    // transform the parameter in the unit hypercube to their physical counterparts according to the prior
    LALInferenceVariables *newParams=NULL;
    newParams=calloc(1,sizeof(LALInferenceVariables));
    /* Make a copy of the parameters passed through currentParams */
    LALInferenceCopyVariables(runStateGlobal->currentParams,newParams);
    int i = runStateGlobal->CubeToPrior(runStateGlobal, newParams, Cube, context);

    // if the parameters violate the prior then set likelihood to log(0);
    if( i == 0 )
    {
        *lnew = -DBL_MAX;
        LALInferenceClearVariables(newParams);
        free(newParams);
        return;
    }

    // calculate the loglike
    *lnew=runStateGlobal->likelihood(newParams, runStateGlobal->data, runStateGlobal->templt);
    *lnew -= (*(REAL8 *)LALInferenceGetVariable(runStateGlobal->algorithmParams, "logZnoise"));
    LALInferenceClearVariables(newParams);
    free(newParams);
}

void dumper(UNUSED int *nSamples, UNUSED int *nlive, UNUSED int *nPar, UNUSED double **physLive,
            UNUSED double **posterior, UNUSED double **paramConstr, UNUSED double *maxLogLike,
            double *logZ, UNUSED double *logZerr, void *context)
{
    char **info = (char **)context;
    char *root2=&info[0][0];
    char *header=&info[1][0];
    char outfile[BAMBI_STRLEN];
    FILE *fileout;

    /* Write evidence to file for use by post-processing */
    strcpy(outfile,root2);
    strcat(outfile,"evidence.dat");
    fileout=fopen(outfile,"w");
    fprintf(fileout,"%g\n",*logZ);
    fclose(fileout);

    /* Write header line to file for use by post-processing */
    if (strcmp(header,"DONOTWRITE")!=0) {
        strcpy(outfile,root2);
        strcat(outfile,"params.txt");
        fileout=fopen(outfile,"w");
        fprintf(fileout,"%s\n",header);
        fclose(fileout);
    }
}

void getphysparams(double *Cube, UNUSED int *ndim, UNUSED int *nPar, void *context)
{
    // CubeToPrior function does this and physical params are in first ndim already
    LALInferenceVariables *newParams=NULL;
    newParams=calloc(1,sizeof(LALInferenceVariables));
    LALInferenceCopyVariables(runStateGlobal->currentParams,newParams);
    runStateGlobal->CubeToPrior(runStateGlobal, newParams, Cube, context);
    free(newParams);

    // Adjust time if necessary
    char **info = (char **)context;
    char *timeID = &info[2][0];
    int id = atoi(timeID);
    if (id >= 0) {
        REAL8 trigtime = *(REAL8 *)LALInferenceGetVariable(runStateGlobal->priorArgs,"trigtime");
        Cube[id] -= trigtime;
    }
}

void getallparams(double *Cube, UNUSED int *ndim, UNUSED int *nPar, void *context)
{
    // CubeToPrior function does this
    LALInferenceVariables *newParams=NULL;
    newParams=calloc(1,sizeof(LALInferenceVariables));
    LALInferenceCopyVariables(runStateGlobal->currentParams,newParams);
    runStateGlobal->CubeToPrior(runStateGlobal, newParams, Cube, context);
    free(newParams);
}

/* MultiNestAlgorithm implements the MultiNest algorithm*/
void LALInferenceMultiNestAlgorithm(LALInferenceRunState *runState)
{
    UINT4 Nlive=*(UINT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nlive");
    REAL8 eff=*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"eff");
    REAL8 MNTol=*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams,"evidencetol");
    UINT4 Ntrain;
    REAL8 logZnoise;
    UINT4 verbose=0,resval=1;

    if (LALInferenceGetProcParamVal(runState->commandLine, "--correlatedGaussianLikelihood")
         || LALInferenceGetProcParamVal(runState->commandLine, "--bimodalGaussianLikelihood")
         || LALInferenceGetProcParamVal(runState->commandLine, "--rosenbrockLikelihood")) {
        logZnoise=0.0;
    } else {
        logZnoise=LALInferenceNullLogLikelihood(runState->data);
    }

    LALInferenceAddVariable(runState->algorithmParams,"logZnoise",&logZnoise,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
    //logLikelihoods=(REAL8 *)(*(REAL8Vector **)LALInferenceGetVariable(runState->algorithmParams,"logLikelihoods"))->data;

    //verbose=LALInferenceCheckVariable(runState->algorithmParams,"verbose");
    if (LALInferenceGetProcParamVal(runState->commandLine, "--progress"))
        verbose=1;

    if (LALInferenceGetProcParamVal(runState->commandLine, "--noresume"))
        resval=0;

    /* output file root */
    ProcessParamsTable *ppt=LALInferenceGetProcParamVal(runState->commandLine,"--outfile");
    if(!ppt){
        fprintf(stderr,"Must specify --outfile <fileroot>\n");
        exit(1);
    }
    char *outfilestr=ppt->value;

    // Do BAMBI?
    doBAMBI = 0;
    if (LALInferenceGetProcParamVal(runState->commandLine, "--BAMBI"))
        doBAMBI=1;

    // Use saved NN?
    useNN = 0;
    if (LALInferenceGetProcParamVal(runState->commandLine, "--useNNslow"))
        useNN = 1;
    else if (LALInferenceGetProcParamVal(runState->commandLine, "--useNNfast"))
        useNN = 2;

    /* NN settings file */
    char *netfilestr=NULL;
    ppt=LALInferenceGetProcParamVal(runState->commandLine,"--NNfile");
    if(doBAMBI || useNN) {
        if (!ppt){
            fprintf(stderr,"Must specify --NNfile <filename.inp>\n");
            exit(1);
        } else
            netfilestr=ppt->value;
    }

    if (LALInferenceCheckVariable(runState->algorithmParams,"Ntrain"))
        Ntrain=*(UINT4 *)LALInferenceGetVariable(runState->algorithmParams,"Ntrain");
    else if (doBAMBI) {
        fprintf(stderr,"Must specify --Ntrain <# points> when running BAMBI, setting it equal to Nlive.\n");
        Ntrain=Nlive;
    } else
        Ntrain=50;

    REAL8 tmin,tmax,tmid;
    if (LALInferenceCheckVariable(runState->priorArgs,"time"))
    {
        LALInferenceGetMinMaxPrior(runState->priorArgs, "time", (void *)&tmin, (void *)&tmax);
        tmid=(tmax+tmin)/2.0;
        LALInferenceAddVariable(runState->priorArgs,"trigtime",&tmid,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
    }

    runStateGlobal = runState;

    // find out the dimensionality of the problem
    int ND = 0;
    LALInferenceVariableItem *item=runState->currentParams->head;
    for(;item;item=item->next)
    {
        if(item->vary==LALINFERENCE_PARAM_LINEAR || item->vary==LALINFERENCE_PARAM_CIRCULAR) ND++;
    }

    if( ND==0 )
    {
        double like = runState->likelihood(runState->currentParams,runState->data,runState->templt);
        like -= (*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams, "logZnoise"));
        fprintf(stdout,"LOG-LIKELIHOOD VALUE RETURNED = %g\n",like);
        double prior = runState->prior(runState,runState->currentParams);
        fprintf(stdout,"LOG-PRIOR VALUE RETURNED = %g\n",prior);
        fprintf(stdout,"LOG-POSTERIOR VALUE RETURNED = %g\n",like+prior);
        return;
    }

    int mmodal = 0;
    int maxModes = 1;
    ppt=LALInferenceGetProcParamVal(runState->commandLine,"--multimodal");
    if(ppt){
        mmodal = 1;
        maxModes = atoi(ppt->value);
    }
    int ceff = 0;
    int nlive = Nlive;
    double efr = eff;
    double mntol = MNTol;
    int ndims = ND;
    int nPar = ndims + 3;
    int nClsPar = fmin(2,ND);
    int updInt = Ntrain;
    double Ztol = -1.e90;
    int pWrap[ndims];
    item=runState->currentParams->head;
    int k = -1;
    for(;item;item=item->next)
    {
        if(item->vary==LALINFERENCE_PARAM_LINEAR || item->vary==LALINFERENCE_PARAM_CIRCULAR)
        {
            k++;
            if(item->vary==LALINFERENCE_PARAM_CIRCULAR)
                pWrap[k] = 1;
            else
                pWrap[k] = 0;
        }
    }
    root=(char *)malloc(BAMBI_STRLEN*sizeof(char));
    networkinputs=(char *)malloc(BAMBI_STRLEN*sizeof(char));
    strcpy(root,outfilestr);
    if (netfilestr!=NULL)
        strcpy(networkinputs,netfilestr);
    int rseed = -1;
    int fb = verbose;
    int bresume = resval;
    int outfile = 1;
    int initMPI = 0;
#ifdef PARALLEL
    initMPI = 1;
#endif
    logZero = -DBL_MAX;
    int maxiter = 0;
    char **info;
    info=(char **)malloc(3*sizeof(char *));
    info[0]=(char *)malloc(BAMBI_STRLEN*sizeof(char));
    info[1]=(char *)malloc(150*sizeof(char));
    info[2]=(char *)malloc(5*sizeof(char));
    strcpy(&info[0][0],outfilestr);
    strcpy(&info[1][0],"DONOTWRITE");
    strcpy(&info[2][0],"-1");
    void *context = (void *)info;

    /* Read injection XML file for parameters if specified */
    /* Used to print injection likelihood and prior */
    ppt=LALInferenceGetProcParamVal(runState->commandLine,"--inj");
    if(ppt){
      SimInspiralTable *injTable=NULL;
      SimInspiralTableFromLIGOLw(&injTable,ppt->value,0,0);
      if(!injTable){
        fprintf(stderr,"Unable to open injection file %s\n",ppt->value);
        exit(1);
      }
      ppt=LALInferenceGetProcParamVal(runState->commandLine,"--event");
      if(ppt){
        int event= atoi(ppt->value);
        int i=0;
        while(i<event) {i++; injTable=injTable->next;} /* select event */
       LALInferenceVariables tempParams;
    //memset(&tempParams,0,sizeof(tempParams));
    LALInferenceInjectionToVariables(injTable,&tempParams);
    LALInferenceVariableItem *node=NULL;
    item=runState->currentParams->head;
    for(;item;item=item->next) {
      node=LALInferenceGetItem(&tempParams,item->name);
      if(node) {
        LALInferenceSetVariable(runState->currentParams,node->name,node->value);
        /*if(strstr(node->name,"LAL")==NULL)
        fprintf(stdout,"Injection variable %s = %g\n",node->name,*(double *)(node->value));
        else
        fprintf(stdout,"Injection variable %s = %d\n",node->name,*(int *)(node->value));*/
      }
    }
    double linj,pinj,lz;
        char finjname[150];
        sprintf(finjname,"%sinjlike.txt",root);
    linj=runState->likelihood(runState->currentParams, runState->data, runState->templt);
        lz = (*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams, "logZnoise"));
    linj -= lz;
    pinj = runState->prior(runState,runState->currentParams);
    FILE *finj=fopen(finjname,"w");
    fprintf(finj,"Log-likelihood value returned = %g\n",linj);
    fprintf(finj,"Log-prior value returned = %g\n",pinj);
    fprintf(finj,"Log-posterior value returned = %g\n",linj+pinj);
    fprintf(finj,"Noise log-evidence value = %g\n",lz);
    fclose(finj);
      }
    }

    BAMBIRun(mmodal, ceff, nlive, mntol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, rseed, pWrap, fb,
    bresume, outfile, initMPI, logZero, maxiter, LogLike, dumper, bambi, context);

    free(info[1]);free(info[0]);free(info);free(root);free(networkinputs);
}


LALInferenceRunState *initialize(ProcessParamsTable *commandLine)
/* calls the "ReadData()" function to gather data & PSD from files, */
/* and initializes other variables accordingly.                     */
{
    char help[]="\
Initialisation arguments:\n\
(--verbose [N])\tOutput more info. N=1: errors, N=2 (default): warnings, N=3: info \n\
(--randomseed seed           Random seed)\n\n";
    LALInferenceRunState *irs=NULL;
    LALInferenceIFOData *ifoPtr, *ifoListStart;
    ProcessParamsTable *ppt=NULL;

    irs = XLALCalloc(1, sizeof(LALInferenceRunState));
    irs->commandLine=commandLine;

    /* Initialise parameters structure */
    irs->algorithmParams=XLALCalloc(1,sizeof(LALInferenceVariables));
    irs->priorArgs=XLALCalloc(1,sizeof(LALInferenceVariables));
    irs->proposalArgs=XLALCalloc(1,sizeof(LALInferenceVariables));

    INT4 verbose=0;
    INT4 x=0;
    ppt=LALInferenceGetProcParamVal(commandLine,"--verbose");
    if(ppt) {
      if(ppt->value){
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
      LALInferenceAddVariable(irs->algorithmParams,"verbose", &verbose , LALINFERENCE_INT4_t,
                  LALINFERENCE_PARAM_FIXED);
    }

    /* read data from files */
    /* (this will already initialise each LALIFOData's following elements:  */
        ppt=LALInferenceGetProcParamVal(commandLine,"--help");
        if(ppt)
        {
                fprintf(stdout,"%s",help);
        irs->data = LALInferenceReadData(commandLine);
                return(irs);
        }
    else
    {
        fprintf(stdout, " readData(): started.\n");
            irs->data = LALInferenceReadData(commandLine);
    }

    /*     fLow, fHigh, detector, timeToFreqFFTPlan, freqToTimeFFTPlan,     */
    /*     window, oneSidedNoisePowerSpectrum, timeDate, freqData         ) */
    fprintf(stdout, " LALInferenceReadData(): finished.\n");
    if (irs->data != NULL) {
        fprintf(stdout, " initialize(): successfully read data.\n");

        fprintf(stdout, " LALInferenceInjectInspiralSignal(): started.\n");
        LALInferenceInjectInspiralSignal(irs->data,commandLine);
        fprintf(stdout, " LALInferenceInjectInspiralSignal(): finished.\n");

        ifoPtr = irs->data;
        ifoListStart = irs->data;
        while (ifoPtr != NULL) {
            /*If two IFOs have the same sampling rate, they should have the same timeModelh*,
             freqModelh*, and modelParams variables to avoid excess computation
             in model waveform generation in the future*/
            LALInferenceIFOData * ifoPtrCompare=ifoListStart;
            int foundIFOwithSameSampleRate=0;
            while(ifoPtrCompare != NULL && ifoPtrCompare!=ifoPtr) {
                if(ifoPtrCompare->timeData->deltaT == ifoPtr->timeData->deltaT){
                    ifoPtr->timeModelhPlus=ifoPtrCompare->timeModelhPlus;
                    ifoPtr->freqModelhPlus=ifoPtrCompare->freqModelhPlus;
                    ifoPtr->timeModelhCross=ifoPtrCompare->timeModelhCross;
                    ifoPtr->freqModelhCross=ifoPtrCompare->freqModelhCross;
                    ifoPtr->modelParams=ifoPtrCompare->modelParams;
                    foundIFOwithSameSampleRate=1;
                    break;
                }
            }
            if(!foundIFOwithSameSampleRate){
                ifoPtr->timeModelhPlus  = XLALCreateREAL8TimeSeries("timeModelhPlus",&(ifoPtr->timeData->epoch),0.0,ifoPtr->timeData->deltaT,&lalDimensionlessUnit,ifoPtr->timeData->data->length);
                ifoPtr->timeModelhCross = XLALCreateREAL8TimeSeries("timeModelhCross",&(ifoPtr->timeData->epoch),0.0,    ifoPtr->timeData->deltaT,&lalDimensionlessUnit,ifoPtr->timeData->data->length);
                ifoPtr->freqModelhPlus = XLALCreateCOMPLEX16FrequencySeries("freqModelhPlus",&(ifoPtr->freqData->epoch),0.0,ifoPtr->freqData->deltaF,&lalDimensionlessUnit,ifoPtr->freqData->data->length);
                ifoPtr->freqModelhCross = XLALCreateCOMPLEX16FrequencySeries("freqModelhCross",&(ifoPtr->freqData->epoch), 0.0, ifoPtr->freqData->deltaF, &lalDimensionlessUnit,ifoPtr->freqData->data->length);
                ifoPtr->modelParams = XLALCalloc(1, sizeof(LALInferenceVariables));
            }
            ifoPtr = ifoPtr->next;
        }
        irs->currentLikelihood=LALInferenceNullLogLikelihood(irs->data);
        printf("Null Log Likelihood: %g\n", irs->currentLikelihood);
    }
    else
    {
        fprintf(stdout, " initialize(): no data read.\n");
        exit(1);
    }

    /* set up GSL random number generator: */
    unsigned long int randomseed;
    struct timeval tv;
    FILE *devrandom;
    gsl_rng_env_setup();
    irs->GSLrandom = gsl_rng_alloc(gsl_rng_mt19937);
    /* (try to) get random seed from command line: */
    ppt = LALInferenceGetProcParamVal(commandLine, "--randomseed");
    if (ppt != NULL)
        randomseed = atoi(ppt->value);
    else { /* otherwise generate "random" random seed: */
        if ((devrandom = fopen("/dev/random","r")) == NULL) {
            gettimeofday(&tv, 0);
            randomseed = tv.tv_sec + tv.tv_usec;
        }
        else {
            if(1!=fread(&randomseed, sizeof(randomseed), 1, devrandom)){
              fprintf(stderr,"Error: Unable to read random seed from /dev/random\n");
              exit(1);
            }
            fclose(devrandom);
        }
    }
    fprintf(stdout, " initialize(): random seed: %lu\n", randomseed);
    gsl_rng_set(irs->GSLrandom, randomseed);

    return(irs);
}

/***** Initialise MultiNest structures *****/
/************************************************/
void initializeMN(LALInferenceRunState *runState)
{
    char help[]="\
               ---------------------------------------------------------------------------------------------------\n\
               --- General Algorithm Parameters ------------------------------------------------------------------\n\
               ---------------------------------------------------------------------------------------------------\n\
               --Nlive N                        Number of live points to use.\n\
               (--Ntrain N)                     Number of training points to use for NN (default=Nlive).\n\
               (--eff e)                        Target efficiency (0.1)\n\
               (--tol tol)                      Tolerance on evidence calculation (0.5)\n\
               (--multimodal maxModes)          Enables multimodal sampling with specified maximum number of modes\n\
                                                  (default is turned off with 1 mode)\
               (--progress)                     Produce progress information.\n\
               (--noresume)                     Do not resume on previous run.\n\
               (--BAMBI)                        Use BAMBI instead of just MultiNest\n\
               (--NNfile filename)              Use specified file for neureal network training inputs.\n\
               (--useNNslow)                    Use previously saved NN - slow method\n\
               (--useNNfast)                    Use previously saved NN - fast method\n\
               \n\
               ---------------------------------------------------------------------------------------------------\n\
               --- Likelihood Functions --------------------------------------------------------------------------\n\
               ---------------------------------------------------------------------------------------------------\n\
               (--zeroLogLike)                  Use flat, null likelihood.\n\
               (--studentTLikelihood)           Use the Student-T Likelihood that marginalizes over noise.\n\
               (--correlatedGaussianLikelihood) Use analytic, correlated Gaussian for Likelihood.\n\
               (--bimodalGaussianLikelihood)    Use analytic, bimodal correlated Gaussian for Likelihood.\n\
               (--rosenbrockLikelihood)         Use analytic, Rosenbrock banana for Likelihood.\n\
               ---------------------------------------------------------------------------------------------------\n\
               --- Prior Functions -------------------------------------------------------------------------------\n\
               ---------------------------------------------------------------------------------------------------\n\
                                                Default prior is currently the S6 prior.\n\
               (--S6Prior)                      Use prior from S6 analysis.\n\
               (--skyLocPrior)                  Use prior from sky localisation project.\n\
               (--AnalyticPrior)                Use prior for analytic likelihood tests.\n\
               \n\
               ---------------------------------------------------------------------------------------------------\n\
               --- Output ----------------------------------------------------------------------------------------\n\
               ---------------------------------------------------------------------------------------------------\n\
               --outfile file                   Root for output files.\n";

    ProcessParamsTable *ppt=NULL;
    ProcessParamsTable *commandLine=runState->commandLine;


    /* Print command line arguments if help requested */
    ppt=LALInferenceGetProcParamVal(commandLine,"--help");
    if(ppt)
    {
        fprintf(stdout,"%s",help);
        return;
    }

    INT4 tmpi=0;
    REAL8 tmpd=0;


    /* Initialise parameters structure */
    runState->algorithmParams=XLALCalloc(1,sizeof(LALInferenceVariables));
    runState->priorArgs=XLALCalloc(1,sizeof(LALInferenceVariables));


    /* Set up the appropriate functions for MultiNest */
    runState->algorithm=&LALInferenceMultiNestAlgorithm;


    /* Set up the prior function */
    if(LALInferenceGetProcParamVal(commandLine,"--skyLocPrior")){
        runState->prior=&LALInferenceInspiralSkyLocPrior;
        runState->CubeToPrior = &LALInferenceInspiralSkyLocCubeToPrior;
    } else if (LALInferenceGetProcParamVal(commandLine, "--S6Prior")) {
        runState->prior=&LALInferenceInspiralPriorNormalised;
        runState->CubeToPrior = &LALInferenceInspiralPriorNormalisedCubeToPrior;
    } else if (LALInferenceGetProcParamVal(commandLine, "--AnalyticPrior")) {
        runState->prior = &LALInferenceAnalyticNullPrior;
        runState->CubeToPrior = &LALInferenceAnalyticCubeToPrior;
    } else {
        runState->prior = &LALInferenceInspiralPriorNormalised;
        runState->CubeToPrior = &LALInferenceInspiralPriorNormalisedCubeToPrior;
    }

    if (LALInferenceGetProcParamVal(commandLine, "--correlatedGaussianLikelihood") ||
        LALInferenceGetProcParamVal(commandLine, "--bimodalGaussianLikelihood") ||
        LALInferenceGetProcParamVal(commandLine, "--rosenbrockLikelihood"))
    {
        runState->prior = &LALInferenceAnalyticNullPrior;
        runState->CubeToPrior = &LALInferenceAnalyticCubeToPrior;
    }


    /* Number of live points */
    //printf("set number of live points.\n");
    ppt=LALInferenceGetProcParamVal(commandLine,"--Nlive");
    if (!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--nlive");
    if(ppt)
        tmpi=atoi(ppt->value);
    else {
        fprintf(stderr,"Error, must specify number of live points\n");
        exit(1);
    }
    LALInferenceAddVariable(runState->algorithmParams,"Nlive",&tmpi, LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
    printf("Setting Nlive = %d\n",tmpi);

    /* Target efficiency */
    ppt=LALInferenceGetProcParamVal(commandLine,"--eff");
    if(ppt)
        tmpd=fabs(atof(ppt->value));
    else {
        tmpd=0.1;
    }
    LALInferenceAddVariable(runState->algorithmParams,"eff",&tmpd, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);

    /* Tolerance of evidence calculation */
    ppt=LALInferenceGetProcParamVal(commandLine,"--tol");
    if(ppt)
	tmpd=fabs(atof(ppt->value));
    else {
	tmpd=0.5;
    }
    LALInferenceAddVariable(runState->algorithmParams,"evidencetol",&tmpd, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);

    /* Number of training points */
    ppt=LALInferenceGetProcParamVal(commandLine,"--Ntrain");
    if(ppt) {
        tmpi=atoi(ppt->value);
        LALInferenceAddVariable(runState->algorithmParams,"Ntrain",&tmpi, LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
    }

    return;

}

/** Initialise student-t extra variables, set likelihood */
void initStudentt(LALInferenceRunState *state)
{
        char help[]="\
Student T Likelihood Arguments:\n\
(--studentTLikelihood)\tUse student-t likelihood function\n";

    ProcessParamsTable *ppt=NULL;
    LALInferenceIFOData *ifo=state->data;

    /* Print command line arguments if help requested */
        if(LALInferenceGetProcParamVal(state->commandLine,"--help"))
        {
                fprintf(stdout,"%s",help);
        while(ifo) {
            fprintf(stdout,"(--dof-%s DoF)\tDegrees of freedom for %s\n",ifo->name,ifo->name);
            ifo=ifo->next;
        }
        return;
        }
    /* Don't do anything unless asked */
    if(!LALInferenceGetProcParamVal(state->commandLine,"--studentTLikelihood")) return;

    /* initialise degrees of freedom parameters for each IFO */
    while(ifo){
        CHAR df_argument_name[128];
        CHAR df_variable_name[64];
        REAL8 dof=10.0; /* Degrees of freedom parameter */

        sprintf(df_argument_name,"--dof-%s",ifo->name);
        if((ppt=LALInferenceGetProcParamVal(state->commandLine,df_argument_name)))
            dof=atof(ppt->value);
            sprintf(df_variable_name,"df_%s",ifo->name);
            LALInferenceAddVariable(state->currentParams,df_variable_name,&dof,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
        fprintf(stdout,"Setting %lf degrees of freedom for %s\n",dof,ifo->name);
        ifo=ifo->next;
    }

    /* Set likelihood to student-t */
    state->likelihood = &LALInferenceFreqDomainStudentTLogLikelihood;

    /* Set the noise model evidence to the student t model value */
    LALInferenceTemplateNullFreqdomain(state->data);
    REAL8 noiseZ=LALInferenceFreqDomainStudentTLogLikelihood(state->currentParams,state->data,&LALInferenceTemplateNullFreqdomain);
    LALInferenceAddVariable(state->algorithmParams,"logZnoise",&noiseZ,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
    fprintf(stdout,"Student-t Noise evidence %lf\n",noiseZ);

    return;
}

/*************** MAIN **********************/


int main(int argc, char *argv[]){
        char help[]="LALInferenceBAMBI:\n\
Bayesian analysis tool using BAMBI algorithm\n\
for CBC analysis. Uses LALInference library for back-end.\n\n\
Arguments for each section follow:\n\n";

    LALInferenceRunState *state;
    ProcessParamsTable *procParams=NULL;

    /* Read command line and parse */
    procParams=LALInferenceParseCommandLine(argc,argv);

    /* Print command line arguments if help requested */
        if(LALInferenceGetProcParamVal(procParams,"--help")) fprintf(stdout,"%s",help);

        /* initialise runstate based on command line */
    /* This includes reading in the data */
    /* And performing any injections specified */
    /* And allocating memory */
    state = initialize(procParams);

    /* Set template function */
    LALInferenceInitCBCTemplate(state);

    /* Set up structures for MultiNest */
    initializeMN(state);

    /* Set up currentParams with variables to be used */
    /* Review task needs special priors */
    if(LALInferenceGetProcParamVal(procParams,"--correlatedGaussianLikelihood"))
        LALInferenceInitVariablesReviewEvidence(state);
    else if(LALInferenceGetProcParamVal(procParams,"--bimodalGaussianLikelihood"))
        LALInferenceInitVariablesReviewEvidence_bimod(state);
    else if(LALInferenceGetProcParamVal(procParams,"--rosenbrockLikelihood"))
        LALInferenceInitVariablesReviewEvidence_banana(state);
    else
        LALInferenceInitCBCVariables(state);

    /* Choose the likelihood */
    LALInferenceInitLikelihood(state);

    /* Exit if help requested */
    if(LALInferenceGetProcParamVal(state->commandLine,"--help")) exit(0);

    /* Call MultiNest algorithm */
    state->algorithm(state);

    /* end */
    return 0;
}
