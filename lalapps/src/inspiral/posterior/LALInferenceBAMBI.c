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
#include <lal/LALInferenceCalibrationErrors.h>

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

void initializeMN(LALInferenceRunState *runState);

/******** Defined for BAMBI (start) **********/

#define BAMBI_STRLEN 1000

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
void setParams(double *Cube, LALInferenceVariables *params, void *context, bool allparams);
void runTestLikelihood(LALInferenceRunState *runState);
void countDimensions(LALInferenceVariables *params, int *nsamp, int *nextra);

void getLogLike(double *Cube, UNUSED int *ndim, UNUSED int *npars, double *lnew, void *context)
{
    // transform the parameter in the unit hypercube to their physical counterparts according to the prior
    LALInferenceVariables *newParams=NULL;
    newParams=calloc(1,sizeof(LALInferenceVariables));
    /* Make a copy of the parameters passed through currentParams */
    LALInferenceCopyVariables(runStateGlobal->threads[0]->currentParams,newParams);
    int i = runStateGlobal->CubeToPrior(runStateGlobal, newParams, runStateGlobal->threads[0]->model, Cube, context);

    // if the parameters violate the prior then set likelihood to log(0);
    if( i == 0 )
    {
        *lnew = -DBL_MAX;
        LALInferenceClearVariables(newParams);
        free(newParams);
        return;
    }

    // calculate the loglike
    *lnew=runStateGlobal->likelihood(newParams, runStateGlobal->data, runStateGlobal->threads[0]->model);
    *lnew -= (*(REAL8 *)LALInferenceGetVariable(runStateGlobal->algorithmParams, "logZnoise"));
    setParams(Cube, newParams, context, true);
    LALInferenceClearVariables(newParams);
    free(newParams);
}

void dumper(UNUSED int *nSamples, UNUSED int *nlive, UNUSED int *nPar, UNUSED double **physLive,
            UNUSED double **posterior, UNUSED double **paramConstr, UNUSED double *maxLogLike,
            double *logZ, double *logZerr, void *context)
{
    char **info = (char **)context;
    char *root2=&info[0][0];
    char *header=&info[1][0];
    char outfile[BAMBI_STRLEN];
    FILE *fileout = NULL;

    /* Write evidence to file for use by post-processing */
    double logZnoise = (*(REAL8 *)LALInferenceGetVariable(runStateGlobal->algorithmParams, "logZnoise"));
    strcpy(outfile,root2);
    strcat(outfile,"evidence.dat");
    fileout=fopen(outfile,"w");
    fprintf(fileout,"%lf\t%lf\t%lf\n",*logZ,*logZerr,logZnoise);
    fclose(fileout);

    /* Write header line to file for use by post-processing */
    if (strcmp(header,"DONOTWRITE")!=0) {
        strcpy(outfile,root2);
        strcat(outfile,"params.txt");
        fileout=fopen(outfile,"w");
        fprintf(fileout,"%s\n",header);
        fclose(fileout);
    }

    /* Prints stats file with template and likelihood evaluation counts */
    /*sprintf(outfile,"%srunstats.txt",root2);
    fileout=fopen(outfile,"w");
    fprintf(fileout,"IFO templates likelihoods\n");
    for(LALInferenceIFOData *p=runStateGlobal->data;p;p=p->next)
        fprintf(fileout,"%s: %u %u\n",p->name,p->templa_counter,p->likeli_counter);
    fclose(fileout);*/
}

void getphysparams(double *Cube, UNUSED int *ndim, UNUSED int *nPar, void *context)
{
    // CubeToPrior function does this and physical params are in first ndim already
    LALInferenceVariables *newParams=NULL;
    newParams=calloc(1,sizeof(LALInferenceVariables));
    LALInferenceCopyVariables(runStateGlobal->threads[0]->currentParams,newParams);
    runStateGlobal->CubeToPrior(runStateGlobal, newParams, runStateGlobal->threads[0]->model, Cube, context);
    setParams(Cube, newParams, context, false);
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
    LALInferenceCopyVariables(runStateGlobal->threads[0]->currentParams,newParams);
    runStateGlobal->CubeToPrior(runStateGlobal, newParams, runStateGlobal->threads[0]->model, Cube, context);
    setParams(Cube, newParams, context, true);
    free(newParams);
}

void setParams(double *Cube, LALInferenceVariables *params, void *context, bool allparams)
{
    char **info = (char **)context;
    char *header = &info[1][0];
    strcpy(header, "");

    UINT4 i = 0, j, k;
    char name[100];
    LALInferenceVariableItem *item=params->head;
    for(;item;item=item->next)
    {
        if (allparams || item->vary == LALINFERENCE_PARAM_LINEAR || item->vary == LALINFERENCE_PARAM_CIRCULAR)
        {
            switch (item->type)
            {
                case LALINFERENCE_INT4_t:
                    Cube[i++] = (double) (*(INT4 *) item->value);
                    strcat(header, LALInferenceTranslateInternalToExternalParamName(item->name));
                    strcat(header, " ");
                    break;
                case LALINFERENCE_INT8_t:
                    Cube[i++] = (double) (*(INT8 *) item->value);
                    strcat(header, LALInferenceTranslateInternalToExternalParamName(item->name));
                    strcat(header, " ");
                    break;
                case LALINFERENCE_UINT4_t:
                    Cube[i++] = (double) (*(UINT4 *) item->value);
                    strcat(header, LALInferenceTranslateInternalToExternalParamName(item->name));
                    strcat(header, " ");
                    break;
                case LALINFERENCE_REAL4_t:
                    Cube[i++] = (double) (*(REAL4 *) item->value);
                    strcat(header, LALInferenceTranslateInternalToExternalParamName(item->name));
                    strcat(header, " ");
                    break;
                case LALINFERENCE_REAL8_t:
                    Cube[i++] = (double) (*(REAL8 *) item->value);
                    strcat(header, LALInferenceTranslateInternalToExternalParamName(item->name));
                    strcat(header, " ");
                    break;
                case LALINFERENCE_COMPLEX8_t: ; // empty statement
                    COMPLEX8 temp1 = *(COMPLEX8 *) item->value;
                    Cube[i++] = (double) creal(temp1);
                    Cube[i++] = (double) cimag(temp1);
                    sprintf(name, "%s_real %s_imag", LALInferenceTranslateInternalToExternalParamName(item->name),
                        LALInferenceTranslateInternalToExternalParamName(item->name));
                    strcat(header, name);
                    strcat(header, " ");
                    break;
                case LALINFERENCE_COMPLEX16_t: ; // empty statement
                    COMPLEX16 temp2 = *(COMPLEX16 *) item->value;
                    Cube[i++] = (double) creal(temp2);
                    Cube[i++] = (double) cimag(temp2);
                    sprintf(name, "%s_real %s_imag", LALInferenceTranslateInternalToExternalParamName(item->name),
                        LALInferenceTranslateInternalToExternalParamName(item->name));
                    strcat(header, name);
                    strcat(header, " ");
                    break;
                case LALINFERENCE_gslMatrix_t: ; // empty statement
                    gsl_matrix *nparams = *((gsl_matrix **)LALInferenceGetVariable(params,"psdscale"));
                    for (j=0; j<(UINT4)nparams->size1; j++) {
                        for (k=0; k<(UINT4)nparams->size2; k++) {
                            Cube[i++] = gsl_matrix_get(nparams, j, k);
                            sprintf(name,"%s_%d_%d",LALInferenceTranslateInternalToExternalParamName(item->name),j,k);
                            strcat(header,name);
                            strcat(header, " ");
                        }
                    }
                    break;
                case LALINFERENCE_REAL8Vector_t: ; // empty statement
                    REAL8Vector *vector1 = *((REAL8Vector **)item->value);
                    for (j=0; j<(UINT4)vector1->length; j++) {
                        Cube[i++] = (double) vector1->data[j];
                        sprintf(name, "%s_%d", LALInferenceTranslateInternalToExternalParamName(item->name), j);
                        strcat(header,name);
                        strcat(header, " ");
                    }
                    break;
                case LALINFERENCE_INT4Vector_t: ; // empty statement
                    INT4Vector *vector2 = *((INT4Vector **)item->value);
                    for (j=0; j<(UINT4)vector2->length; j++) {
                        Cube[i++] = (double) vector2->data[j];
                        sprintf(name, "%s_%d", LALInferenceTranslateInternalToExternalParamName(item->name), j);
                        strcat(header,name);
                        strcat(header, " ");
                    }
                    break;
                case LALINFERENCE_UINT4Vector_t: ; // empty statement
                    UINT4Vector *vector3 = *((UINT4Vector **)item->value);
                    for (j=0; j<(UINT4)vector3->length; j++) {
                        Cube[i++] = (double) vector3->data[j];
                        sprintf(name, "%s_%d", LALInferenceTranslateInternalToExternalParamName(item->name), j);
                        strcat(header,name);
                        strcat(header, " ");
                    }
                    break;
                case LALINFERENCE_COMPLEX16Vector_t: ; // empty statement
                    COMPLEX16Vector *vector4 = *((COMPLEX16Vector **)item->value);
                    for (j=0; j<(UINT4)vector4->length; j++) {
                        Cube[i++] = (double) creal(vector4->data[j]);
                        Cube[i++] = (double) cimag(vector4->data[j]);
                        sprintf(name, "%s_%d_real %s_%d_imag", LALInferenceTranslateInternalToExternalParamName(item->name), j,
                            LALInferenceTranslateInternalToExternalParamName(item->name), j);
                        strcat(header,name);
                        strcat(header, " ");
                    }
                    break;
                default:
                    break;
            }
        }
    }

    // add prior
    Cube[i++] = runStateGlobal->prior(runStateGlobal, params, runStateGlobal->threads[0]->model);
    strcat(header, "logprior ");

    // add logL
    strcat(header,"logl");
}

void runTestLikelihood(LALInferenceRunState *runState)
{
    // create and allocated params struct
    LALInferenceVariables *params=runState->threads[0]->currentParams;

    // count number of sampling dimensions
    int nd = 0, ne = 0;
    countDimensions(params, &nd, &ne);

    // create and allocate Cube[] of 0.5 (middle of prior)
    double *Cube = NULL;
    Cube = (double *)malloc(nd*sizeof(double));
    int i;
    for (i=0; i<nd; i++) Cube[i] = 0.5;

    // create context var
    char **info=(char **)malloc(3*sizeof(char *));
    info[0]=(char *)malloc(5*sizeof(char));
    info[1]=(char *)malloc(5*sizeof(char));
    info[2]=(char *)malloc(5*sizeof(char));
    strcpy(&info[0][0],"");
    strcpy(&info[1][0],"");
    strcpy(&info[2][0],"-1");
    void *context = (void *)info;

    // run CubeToPrior
    runState->CubeToPrior(runState, params, runState->threads[0]->model, Cube, context);

    // do logLikelihood
    runState->likelihood(params, runState->data, runState->threads[0]->model);

    // free context var
    free(info[2]);free(info[1]);free(info[0]);free(info);

    // free temp Cube
    free(Cube);
}

void countDimensions(LALInferenceVariables *params, int *nsamp, int *nextra)
{
    *nsamp = 0;
    *nextra = 0;

    LALInferenceVariableItem *item=params->head;
    for(;item;item=item->next)
    {
        if(item->vary==LALINFERENCE_PARAM_LINEAR || item->vary==LALINFERENCE_PARAM_CIRCULAR)
        {
            if (item->type == LALINFERENCE_gslMatrix_t)
            {
                gsl_matrix *nparams = *((gsl_matrix **)item->value);
                INT4 numdims = nparams->size1 * nparams->size2;
                *nsamp += numdims;
            }
            else if (item->type == LALINFERENCE_REAL8Vector_t)
            {
                REAL8Vector *vector1 = *((REAL8Vector **)item->value);
                *nsamp += vector1->length;
            }
            else
                (*nsamp)++;
        }
        else
        {
            if (item->type == LALINFERENCE_gslMatrix_t)
            {
                gsl_matrix *nparams = *((gsl_matrix **)item->value);
                INT4 numdims = nparams->size1 * nparams->size2;
                *nextra += numdims;
            }
            else if (item->type == LALINFERENCE_REAL8Vector_t)
            {
                REAL8Vector *vector1 = *((REAL8Vector **)item->value);
                *nextra += vector1->length;
            }
            else
                (*nextra)++;
        }
    }
    printf("%d sampled parameters, %d extra parameters\n", *nsamp, *nextra);
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

    INT4 SKY_FRAME=0;
    if(LALInferenceCheckVariable(runState->threads[0]->currentParams,"SKY_FRAME"))
      SKY_FRAME=*(INT4 *)LALInferenceGetVariable(runState->threads[0]->currentParams,"SKY_FRAME");
    REAL8 tmin,tmax,tmid;
    if (SKY_FRAME==1)
    {
        if (LALInferenceCheckVariable(runState->priorArgs,"t0"))
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "t0", (void *)&tmin, (void *)&tmax);
            tmid=(tmax+tmin)/2.0;
            LALInferenceAddVariable(runState->priorArgs,"trigtime",&tmid,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
        }
    }
    else
    {
        if (LALInferenceCheckVariable(runState->priorArgs,"time"))
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "time", (void *)&tmin, (void *)&tmax);
            tmid=(tmax+tmin)/2.0;
            LALInferenceAddVariable(runState->priorArgs,"trigtime",&tmid,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
        }
    }

    runStateGlobal = runState;

    // run a likelihood before counting dimensions
    runTestLikelihood(runState);

    // find out the dimensionality of the problem
    int ND = 0, Nextra = 0;
    countDimensions(runState->threads[0]->currentParams, &ND, &Nextra);

    if( ND==0 )
    {
        double like = runState->likelihood(runState->threads[0]->currentParams,runState->data,runState->threads[0]->model);
        like -= (*(REAL8 *)LALInferenceGetVariable(runState->algorithmParams, "logZnoise"));
        fprintf(stdout,"LOG-LIKELIHOOD VALUE RETURNED = %g\n",like);
        double prior = runState->prior(runState,runState->threads[0]->currentParams,runState->threads[0]->model);
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
    int nPar = ND + Nextra + 1; // one extra for logprior
    printf("ndims = %d, nPar = %d\n", ndims, nPar);
    //if (LALInferenceCheckVariable(runState->threads[0]->currentParams,"f_ref")) nPar++;  // add space for f_ref
    //if (SKY_FRAME==1) nPar += 3;
    int nClsPar = fmin(2,ND);
    int updInt = Ntrain;
    double Ztol = -1.e90;
    int pWrap[ndims];
    LALInferenceVariableItem *item = NULL;
    item=runState->threads[0]->currentParams->head;
    int k = -1;
    for(;item;item=item->next)
    {
        if(item->vary==LALINFERENCE_PARAM_LINEAR || item->vary==LALINFERENCE_PARAM_CIRCULAR)
        {
            if (item->type == LALINFERENCE_gslMatrix_t)
            {
                gsl_matrix *nparams = *((gsl_matrix **)item->value);
                INT4 numdims = nparams->size1 * nparams->size2;
                INT4 kk;
                if (item->vary==LALINFERENCE_PARAM_CIRCULAR)
                {
                    for (kk=0;kk<numdims;kk++)
                    {
                        k++;
                        pWrap[k] = 1;
                    }
                }
                else
                {
                    for (kk=0;kk<numdims;kk++)
                    {
                        k++;
                        pWrap[k] = 0;
                    }
                }
            }
            else
            {
                k++;
                if(item->vary==LALINFERENCE_PARAM_CIRCULAR)
                    pWrap[k] = 1;
                else
                    pWrap[k] = 0;
            }
        }
    }
    root=(char *)malloc(BAMBI_STRLEN*sizeof(char));
    networkinputs=(char *)malloc(BAMBI_STRLEN*sizeof(char));
    strcpy(root,outfilestr);
    if (netfilestr!=NULL)
        strcpy(networkinputs,netfilestr);
    int rseed = -1;
    ppt=LALInferenceGetProcParamVal(runState->commandLine,"--randomseed");
    if(ppt) rseed = atoi(ppt->value) % 30000;
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
    info[1]=(char *)malloc(nPar*15*sizeof(char));
    info[2]=(char *)malloc(5*sizeof(char));
    strcpy(&info[0][0],outfilestr);
    strcpy(&info[1][0],"DONOTWRITE");
    strcpy(&info[2][0],"-1");
    void *context = (void *)info;

    BAMBIRun(mmodal, ceff, nlive, mntol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, rseed, pWrap, fb,
    bresume, outfile, initMPI, logZero, maxiter, LogLike, dumper, bambi, context);

    free(info[2]);free(info[1]);free(info[0]);free(info);
    free(root);free(networkinputs);
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
               (--eff e)                        Target efficiency (0.01)\n\
               (--tol tol)                      Tolerance on evidence calculation (1.0)\n\
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
        tmpd=0.01;
    }
    LALInferenceAddVariable(runState->algorithmParams,"eff",&tmpd, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);

    /* Tolerance of evidence calculation */
    ppt=LALInferenceGetProcParamVal(commandLine,"--tol");
    if(ppt)
	tmpd=fabs(atof(ppt->value));
    else {
	tmpd=1.0;
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
    state = LALInferenceInitRunState(procParams);
    /* Perform injections if data successful read or created */
    LALInferenceInjectInspiralSignal(state->data, state->commandLine);
    
    /* Simulate calibration errors. 
     * NOTE: this must be called after both ReadData and (if relevant) 
     * injectInspiralTD/FD are called! */
    LALInferenceApplyCalibrationErrors(state->data, state->commandLine);

    /* Set up prior */
    LALInferenceInitCBCPrior(state);
    
    /* Set up structures for MultiNest */
    initializeMN(state);

    /* Set up thread structures */
    LALInferenceInitCBCThreads(state,1);
        
    /* Choose the likelihood and set some auxiliary variables */
    LALInferenceInitLikelihood(state);
    

    /* Exit if help requested */
    if(LALInferenceGetProcParamVal(state->commandLine,"--help")) exit(0);

    /* write injection with noise evidence information from algorithm */
    LALInferencePrintInjectionSample(state);
    
    /* Call MultiNest algorithm */
    state->algorithm(state);

    /* end */
    return 0;
}
