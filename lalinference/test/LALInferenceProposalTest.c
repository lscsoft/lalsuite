/* 
 *  LALInferenceProposalTest.c: Testing the jump propsals in LALInferenceProposal.c
 *
 *  Copyright (C) 2011 Will M. Farr, John Veitch
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

#include <math.h>

#include <stdlib.h>

#include <lal/LALInferenceProposal.h>
#include <lal/LALInference.h>
#include <lal/XLALError.h>
#include <lal/LALInferenceReadData.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/Units.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/GenerateInspiral.h>
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceNestedSampler.h>


/* Comparison function for qsorting the arrays later */
static int cmpREAL8p(const void *p1, const void *p2);
REAL8 PriorCDF(const char *name, const REAL8 x, LALInferenceVariables *priorArgs);
REAL8 mchirpminmax(const REAL8 x, const REAL8 min, const REAL8 max);

REAL8 rSquaredCDF(const REAL8 x, const REAL8 min, const REAL8 max);
REAL8 FlatInSine(const REAL8 x,const REAL8 min,const REAL8 max);
REAL8 FlatInCosine(const REAL8 x,const REAL8 min,const REAL8 max);
REAL8 UniformMinMax(const REAL8 x,const REAL8 min,const REAL8 max);
/** Efficient integer power computation. */
REAL8 pow_int(const REAL8, const INT4);
REAL8
pow_int(const REAL8 x, const INT4 n) {
  if (n < 0) {
    return 1.0/pow_int(x, -n);
  } else if (n == 0) {
    return 1.0;
  } else if (n == 1) {
    return x;
  } else if (n % 2 == 0) {
    REAL8 sqrt_x = pow_int(x, n/2);
    return sqrt_x * sqrt_x;
  } else {
    return x*pow_int(x, n-1);
  }
}

/** Cumulative distribution function for KS statistic.  Algorithm from
    Numerical Recipes, Third Edition by Press, Teukolsky, Vetterling
    and Flannery.  Cambridge University Press, 2007. Section
    6.14.12 */
REAL8 PKS(const REAL8);
REAL8
PKS(const REAL8 z) {
  if (z < 0.0) {
    XLALError("PKS", __FILE__, __LINE__, XLAL_FAILURE);
    return -1.0;
  } else if (z < 0.042) {
    return 0.0;
  } else if (z < 1.18) {
    REAL8 x = exp(-1.23370055013616983/(z*z));
    return 2.25675833419102515*sqrt(-log(x))*(x + pow_int(x,9) + pow_int(x, 25) + pow_int(x, 49));
  } else {
    REAL8 x = exp(-2.0*z*z);
    return 1.0 - 2.0*(x - pow_int(x,4) + pow_int(x,9));
  }
}

/** Compliment of PKS(). */
REAL8 QKS(const REAL8);
REAL8
QKS(const REAL8 z) {
  if (z < 0.0) {
    XLALError("QKS", __FILE__, __LINE__, XLAL_FAILURE);
    return -1.0;
  } else if (z == 0.0) {
    return 1.0;
  } else if (z < 1.18) {
    return 1.0 - PKS(z);
  } else {
    REAL8 x = exp(-2.0*z*z);
    return 2.0*(x - pow_int(x, 4) + pow_int(x, 9));
  }
}

/** Computes the p of the KS-statistic comparing the cumulative
    distribution given by the discrete points in \a points with the
    corresponding analytic cumulative distribution values in \a
    cumValues (assumed to be evaluated at the locations in points).
    The input array \a points must be sorted. */
REAL8 KSPValue(const REAL8Vector *, const REAL8Vector *);
REAL8
KSPValue(const REAL8Vector *points, const REAL8Vector *cumValues) {
  UINT4 i;
  REAL8 maxD = 0.0;

  if (points->length != cumValues->length) {
    XLALError("KSPValue", __FILE__, __LINE__, XLAL_FAILURE);
    return -1.0;
  }

  /* Check for sorted points. */
  for (i = 0; i < points->length-1; i++) {
    if (points->data[i+1] < points->data[i]) {
      XLALError("KSPValue", __FILE__, __LINE__, XLAL_FAILURE);
      return -1.0;
    }
  }

  maxD = 0.0;
  for (i = 0; i < points->length; i++) {
    REAL8 D = fabs((REAL8)i/((REAL8) points->length) - cumValues->data[i]);
    maxD = (D > maxD ? D : maxD);
  }

  return QKS((sqrt(points->length) + 0.12 + 0.11/sqrt(points->length))*maxD);
}

LALInferenceRunState *initialize(ProcessParamsTable *commandLine);

LALInferenceRunState *initialize(ProcessParamsTable *commandLine)
/* calls the "ReadData()" function to gather data & PSD from files, */
/* and initializes other variables accordingly.                     */
{
/*	char help[]="\
	Initialisation arguments:\n\
	(--randomseed seed           Random seed for Nested Sampling)\n\n";
 */
	LALInferenceRunState *irs=NULL;
	LALInferenceIFOData *ifoPtr;
	ProcessParamsTable *ppt=NULL;
	unsigned long int randomseed;
	struct timeval tv;
	FILE *devrandom;
	
	irs = calloc(1, sizeof(LALInferenceRunState));
	irs->commandLine=commandLine;

	/* set up GSL random number generator: */
	gsl_rng_env_setup();
	irs->GSLrandom = gsl_rng_alloc(gsl_rng_mt19937);
	/* (try to) get random seed from command line: */
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
	
	/* Add a site for the inclination-distance jump */
	ifoPtr=calloc(1,sizeof(LALInferenceIFOData));
	ifoPtr->detector=calloc(1,sizeof(LALDetector));
	memcpy(ifoPtr->detector,&lalCachedDetectors[LALDetectorIndexLHODIFF],sizeof(LALDetector));
	irs->data=ifoPtr;
	return(irs);
}

/* Setup the variables to control template generation */
/* Includes specification of prior ranges */
/* Stolen from LALInferenceNest.c */
void initVariables(LALInferenceRunState *state);
void initVariables(LALInferenceRunState *state)
{
	LALStatus status;
	LALInferenceVariables *priorArgs=state->priorArgs;
	state->currentParams=XLALCalloc(1,sizeof(LALInferenceVariables));
	LALInferenceVariables *currentParams=state->currentParams;
	ProcessParamsTable *commandLine=state->commandLine;
	REAL8 endtime;
	ProcessParamsTable *ppt=NULL;
	LALPNOrder PhaseOrder=LAL_PNORDER_THREE_POINT_FIVE;
	Approximant approx=TaylorF2;
	REAL8 Dmin=1.0;
	REAL8 Dmax=100.0;
	REAL8 mcMin=1.0;
	REAL8 mcMax=20.5;
//	REAL8 mMin=1.0,mMax=30.0,MTotMax=35.0;
	REAL8 a_spin2_max=1.0, a_spin1_max=1.0;
	REAL8 a_spin2_min=0.0, a_spin1_min=0.0;
	REAL8 phi_spin1_min=0;
	REAL8 phi_spin1_max=2.0*LAL_PI;
	REAL8 theta_spin1_min=0.;
	REAL8 theta_spin1_max=LAL_PI;	
	REAL8 qMin=0.0;
	REAL8 qMax=1.0;
	REAL8 dt=0.1;            /* Width of time prior */
	REAL8 tmpMin,tmpMax,tmpVal;
	REAL8 one=1.0;
	memset(currentParams,0,sizeof(LALInferenceVariables));
	memset(&status,0,sizeof(LALStatus));
	INT4 enable_spin=0;
	INT4 aligned_spin=0;
	gsl_rng *RNG=state->GSLrandom;
	char help[]="\
	Parameter arguments:\n\
	(--qmin eta)\tMinimum eta\n\
	(--qmax eta)\tMaximum eta\n\
	(--dt time)\tWidth of time prior, centred around trigger (0.1s)\n\
	(--trigtime time)\tTrigger time to use\n\
	(--Dmin dist)\tMinimum distance in Mpc (1)\n\
	(--Dmax dist)\tMaximum distance in Mpc (100)\n\
	(--approx ApproximantorderPN)\tSpecify a waveform to use, (default TaylorF2threePointFivePN)\n\
	(--enable-spin)\tEnable spin parameters\n\
	(--aligned-spin)\tUse only aligned spin parameters (uses spins between -1 and 1)\n\
	(--approx ApproximantphaseOrderPN)\tSet approximant (PhenSpin implicitly enables spin)\n\
	(--s1max SPIN)\tMax magnitude of spin (on both bodies!)\n\
	(--s1min SPIN)\tMin magnitude of spin (on both bodies!)\n";
	
	/* Print command line arguments if help requested */
	ppt=LALInferenceGetProcParamVal(commandLine,"--help");
	if(ppt)
	{
		fprintf(stdout,"%s",help);
		return;
	}
	
	/* Over-ride approximant if user specifies */
	ppt=LALInferenceGetProcParamVal(commandLine,"--approx");
	if(ppt){
		if(strstr(ppt->value,"TaylorF2")) approx=TaylorF2;
		else
		{
			approx = XLALGetApproximantFromString(ppt->value);
			if( (int) approx == XLAL_FAILURE)
				ABORTXLAL(&status); 
		}
                PhaseOrder = XLALGetOrderFromString(ppt->value);
                if ( (int) PhaseOrder == XLAL_FAILURE)
                        ABORTXLAL(&status);
	}
	//fprintf(stdout,"Templates will run using Approximant %i, phase order %i\n",approx,PhaseOrder);
	
	/* Over-ride end time if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--trigtime");
	if(ppt){
		endtime=atof(ppt->value);
	}
	
	/* Over-ride time prior if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--dt");
	if(ppt){
		dt=atof(ppt->value);
	}
	
	/* Over-ride Distance min if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--Dmin");
	if(ppt){
		Dmin=(atof(ppt->value));
	}
	
	/* Over-ride Distance max if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--Dmax");
	if(ppt){
		Dmax=(atof(ppt->value));
	}
	ppt=LALInferenceGetProcParamVal(commandLine,"--qmin");
	if(ppt)
		qMin=atof(ppt->value);
	
	ppt=LALInferenceGetProcParamVal(commandLine,"--qmax");
    if(ppt) qMax=atof(ppt->value);
	/* Over-ride Spin prior if specified*/
	
	ppt=LALInferenceGetProcParamVal(commandLine,"--s1max");
	if(ppt){
		a_spin2_max=atof(ppt->value);
		a_spin1_max=atof(ppt->value);
	}
	ppt=LALInferenceGetProcParamVal(commandLine,"--s1min");
	if(ppt){
		a_spin2_min=atof(ppt->value);
		a_spin1_min=atof(ppt->value);
	}
	/* Over-ride component masses */
	//ppt=LALInferenceGetProcParamVal(commandLine,"--compmin");
	//if(ppt)	mMin=atof(ppt->value);
	//fprintf(stderr,"Mmin %f, Mmax %f\n",mMin,mMax);
    //REAL8 zero=0.;
	//LALInferenceAddVariable(priorArgs,"component_min",&zero,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
	//ppt=LALInferenceGetProcParamVal(commandLine,"--compmax");
	//if(ppt)	mMax=atof(ppt->value);
	//LALInferenceAddVariable(priorArgs,"component_max",&mMax,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
    //ppt=LALInferenceGetProcParamVal(commandLine,"--MTotMax");
    //if(ppt)	MTotMax=atof(ppt->value);
    //LALInferenceAddVariable(priorArgs,"MTotMax",&MTotMax,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
	
	endtime=10.;
	
	LALInferenceAddVariable(currentParams, "LAL_APPROXIMANT", &approx,        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
	LALInferenceAddVariable(currentParams, "LAL_PNORDER",     &PhaseOrder,        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
	
	/* Set up the variable parameters */
	tmpVal=(mcMin+gsl_rng_uniform(RNG)*(mcMax-mcMin));
	LALInferenceAddVariable(currentParams,"chirpmass",&tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	LALInferenceAddMinMaxPrior(priorArgs,	"chirpmass",	&mcMin,	&mcMax,		LALINFERENCE_REAL8_t);
	
	tmpVal=gsl_rng_uniform(RNG);
	LALInferenceAddVariable(currentParams, "asym_massratio",       &tmpVal,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	LALInferenceAddMinMaxPrior(priorArgs,	"asym_massratio",	&qMin,	&qMax,	LALINFERENCE_REAL8_t);
	endtime=gsl_rng_uniform(RNG)*dt+endtime-0.5*dt;
	LALInferenceAddVariable(currentParams, "time",            &endtime   ,           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR); 
	tmpMin=endtime-0.5*dt; tmpMax=endtime+0.5*dt;
	LALInferenceAddMinMaxPrior(priorArgs, "time",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);	
	
	tmpVal=gsl_rng_uniform(RNG)*LAL_TWOPI;
	LALInferenceAddVariable(currentParams, "phase",           &tmpVal,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
	tmpMin=0.0; tmpMax=LAL_TWOPI;
	LALInferenceAddMinMaxPrior(priorArgs, "phase",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
	
	tmpVal=Dmin+(Dmax-Dmin)*gsl_rng_uniform(RNG);
	LALInferenceAddVariable(currentParams,"distance", &tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	LALInferenceAddMinMaxPrior(priorArgs, "distance",     &Dmin, &Dmax,   LALINFERENCE_REAL8_t);
	
	tmpVal=gsl_rng_uniform(RNG)*LAL_TWOPI;
	LALInferenceAddVariable(currentParams, "rightascension",  &tmpVal,      LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
	tmpMin=0.0; tmpMax=LAL_TWOPI;
	LALInferenceAddMinMaxPrior(priorArgs, "rightascension",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
	tmpVal=gsl_rng_uniform(RNG)*LAL_PI -LAL_PI/2.;
	LALInferenceAddVariable(currentParams, "declination",     &tmpVal,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	tmpMin=-LAL_PI/2.0; tmpMax=LAL_PI/2.0;
	LALInferenceAddMinMaxPrior(priorArgs, "declination",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
	tmpVal=gsl_rng_uniform(RNG)*LAL_PI;
	LALInferenceAddVariable(currentParams, "polarisation",    &tmpVal,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
	tmpMin=0.0; tmpMax=LAL_PI;
	LALInferenceAddMinMaxPrior(priorArgs, "polarisation",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
	tmpVal=gsl_rng_uniform(RNG)*LAL_PI;
 	LALInferenceAddVariable(currentParams, "inclination",     &tmpVal,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	tmpMin=0.0; tmpMax=LAL_PI;
	LALInferenceAddMinMaxPrior(priorArgs, "inclination",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
	
	/* Additional parameters for spinning waveforms */
	if(LALInferenceGetProcParamVal(commandLine,"--enable-spin")) enable_spin=1;
	
	/* If aligned spins use magnitude in (-1,1) */
	ppt=LALInferenceGetProcParamVal(commandLine,"--aligned-spin");
	if(ppt) {enable_spin=1; aligned_spin=1; a_spin1_min=-1; a_spin2_min=-1;}
	
	if(enable_spin){
		tmpVal=a_spin1_min+(a_spin1_max-a_spin1_min)*gsl_rng_uniform(RNG);
		LALInferenceAddVariable(currentParams, "a_spin1",		&tmpVal,	LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
		LALInferenceAddMinMaxPrior(priorArgs, "a_spin1",     &a_spin1_min, &a_spin1_max,   LALINFERENCE_REAL8_t); 
		
		tmpVal=a_spin2_min+(a_spin2_max-a_spin2_min)*gsl_rng_uniform(RNG);
		LALInferenceAddVariable(currentParams, "a_spin2",		&tmpVal,	LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR); 
		LALInferenceAddMinMaxPrior(priorArgs, "a_spin2",     &a_spin2_min, &a_spin2_max,   LALINFERENCE_REAL8_t); 
		
		
		if(aligned_spin){ /* Set the spin angles to be parallel to orbital */
			tmpVal=0;
			LALInferenceAddVariable(currentParams,"theta_spin1",&tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
			LALInferenceAddVariable(currentParams,"theta_spin2",&tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
			tmpVal=0;
			LALInferenceAddVariable(currentParams,"phi_spin1",&tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
			LALInferenceAddVariable(currentParams,"phi_spin2",&tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
		}
		else{ /* Use full spinning parameters */
			tmpVal=theta_spin1_min+(theta_spin1_max - theta_spin1_min)*gsl_rng_uniform(RNG);
			
			LALInferenceAddVariable(currentParams,"theta_spin1",	&tmpVal,	LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
			LALInferenceAddMinMaxPrior(priorArgs, "theta_spin1",     &theta_spin1_min, &theta_spin1_max,   LALINFERENCE_REAL8_t); 
			
			tmpVal=theta_spin1_min+(theta_spin1_max - theta_spin1_min)*gsl_rng_uniform(RNG);
			LALInferenceAddVariable(currentParams,"theta_spin2",	&tmpVal,	LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
			LALInferenceAddMinMaxPrior(priorArgs, "theta_spin2",     &theta_spin1_min, &theta_spin1_max,   LALINFERENCE_REAL8_t); 
			
			tmpVal=phi_spin1_min+(phi_spin1_max - phi_spin1_min)*gsl_rng_uniform(RNG);
			
			LALInferenceAddVariable(currentParams,"phi_spin1",		&tmpVal,	LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
			LALInferenceAddMinMaxPrior(priorArgs, "phi_spin1",     &phi_spin1_min, &phi_spin1_max,   LALINFERENCE_REAL8_t); 
			
			tmpVal=phi_spin1_min+(phi_spin1_max - phi_spin1_min)*gsl_rng_uniform(RNG);
			LALInferenceAddVariable(currentParams,"phi_spin2",		&tmpVal,	LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
			LALInferenceAddMinMaxPrior(priorArgs, "phi_spin2",     &phi_spin1_min, &phi_spin1_max,   LALINFERENCE_REAL8_t);
		}
	}
    LALInferenceAddVariable(currentParams,"logL",&one,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
    REAL8 prior=state->prior(state,currentParams);
    LALInferenceAddVariable(currentParams,"logPrior",&prior,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
	return;
}

int main(int argc, char *argv[]) {
  
	char help[]="\
	LALInferenceProposalTest:\n\
	Test jump proposal distributions\n\
	--Nprop N\t: Number of jumps to perform (1000)\n\
	--outfile file.dat\t: Optional output file for samples\n\
	--thin N\t:Thin MCMC chain by factor N (default 1000)\
	";
	LALInferenceRunState *state=NULL;
	ProcessParamsTable *procParams=NULL;
	UINT4 Nmcmc=1000,i=1,NvarArray=0,thinfac=1000;
	REAL8 logLmin=-DBL_MAX;
	ProcessParamsTable *ppt=NULL;
	FILE *outfile=NULL;
	char *filename=NULL;
	LALInferenceVariableItem *param=NULL;
	LALInferenceVariables *varArray=NULL;
    /* Number of points to sprinkle in prior to get covariance matrix, differential evolution steps */
    UINT4 NCOV=100;


	/* Read command line and parse */
	procParams=LALInferenceParseCommandLine(argc,argv);
	/* initialise runstate based on command line */
	/* This includes reading in the data */
	/* And performing any injections specified */
	/* And allocating memory */
	ppt=LALInferenceGetProcParamVal(procParams,"--help");
	if(ppt)
	{
		fprintf(stdout,"%s",help);
	}
	if((ppt=LALInferenceGetProcParamVal(procParams,"--Nprop")))
	  Nmcmc=atoi(ppt->value);
	if((ppt=LALInferenceGetProcParamVal(procParams,"--outfile")))
	  filename=ppt->value;
	if((ppt=LALInferenceGetProcParamVal(procParams,"--thin")))
	  thinfac=atoi(ppt->value);
	
	state = initialize(procParams);
	
	/* Set a prior function */
	state->algorithm=NULL;
	state->evolve=NULL; /* Use MCMC for this? */
	state->template=NULL;
    /* Log the samples to an array for later use */
	state->logsample=&LALInferenceLogSampleToArray;
	state->priorArgs=calloc(1,sizeof(LALInferenceVariables));
	state->proposalArgs=calloc(1,sizeof(LALInferenceVariables));
	state->algorithmParams=calloc(1,sizeof(LALInferenceVariables));
	//state->prior=LALInferenceInspiralPriorNormalised;
	state->prior=LALInferenceInspiralPrior;
	state->likelihood=&LALInferenceZeroLogLikelihood;
	state->proposal=&NSWrapMCMCLALProposal;
	state->proposalStats = calloc(1,sizeof(LALInferenceVariables));
	
	/* Set up a sample to evolve */
    LALInferenceVariables **samples=calloc(sizeof(LALInferenceVariables *),NCOV);
    state->livePoints=samples;
    LALInferenceAddVariable(state->algorithmParams,"Nlive",&NCOV,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
    //LALInferenceVariables *backup=state->currentParams;
    initVariables(state);
    for(i=0;i<NCOV;i++){
        samples[i]=calloc(1,sizeof(LALInferenceVariables));
        LALInferenceCopyVariables(state->currentParams,samples[i]);
        LALInferenceDrawFromPrior(samples[i], state->priorArgs, state->GSLrandom );
        /* scatter points for CVM calculation */
        REAL8 prior = state->prior(state,samples[i]);
	LALInferenceSetVariable(samples[i],"logPrior",&prior);
    }
    state->currentParams=samples[0];
    gsl_matrix **cvm=calloc(1,sizeof(gsl_matrix *));
    /* Add the covariance matrix for proposal distribution */
	LALInferenceNScalcCVM(cvm,state->livePoints,NCOV);
	state->differentialPoints=state->livePoints;
	state->differentialPointsLength=(size_t) NCOV;
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
        XLAL_ERROR(XLAL_EFAILED);
	}
	
	for (i = 0; i < N; i++) {
        eigenValues->data[i] = gsl_vector_get(eValues,i);
	}
	
	LALInferenceAddVariable(state->proposalArgs, "covarianceEigenvectors", &eVectors, LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_FIXED);
	LALInferenceAddVariable(state->proposalArgs, "covarianceEigenvalues", &eigenValues, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED);
	
	LALInferenceAddVariable(state->proposalArgs,"covarianceMatrix",cvm,LALINFERENCE_gslMatrix_t,LALINFERENCE_PARAM_OUTPUT);
	
	/* set up k-D tree if required and not already set */
	LALInferenceSetupkDTreeNSLivePoints( state );
	LALInferenceSetupAdaptiveProposals(state);
    
	/* Set up the proposal function requirements */
	LALInferenceAddVariable(state->algorithmParams,"Nmcmc",&thinfac,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
	LALInferenceAddVariable(state->algorithmParams,"logLmin",&logLmin,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);

	/* Open the output file */
	if(filename) outfile=fopen(filename,"w");
	if(!outfile) fprintf(stdout,"No output file specified, internal testing only\n");
	
	/* Burn in */
	LALInferenceNestedSamplingOneStep(state);
	
	/* Evolve with fixed likelihood */
	for(i=0;i<Nmcmc*thinfac;i++){
	  LALInferenceMCMCSamplePrior(state);
	  /* output sample */
	  if(!(i%thinfac)){
	    if(state->logsample) state->logsample(state,state->currentParams);
	    if(outfile) LALInferencePrintSample(outfile,state->currentParams);
	    if(outfile) fprintf(outfile,"\n");
	  } 
	}
    if(outfile) fclose(outfile);
	outfile=fopen("headers.txt","w");
    LALInferenceFprintParameterNonFixedHeaders(outfile,state->currentParams);
    fclose(outfile);
    
    /* Perform K-S test for parameters with analytic distributions */
    varArray = *(LALInferenceVariables**)LALInferenceGetVariable(state->algorithmParams,"outputarray");
    UINT4* pNvarArray = (UINT4 *)LALInferenceGetVariable(state->algorithmParams,"N_outputarray");
    int exists=LALInferenceCheckVariable(state->algorithmParams,"N_outputarray");
    printf("%d\n",exists);    
    NvarArray=*pNvarArray;
    
    if(NvarArray!=Nmcmc) printf("ERROR: Not all iterations were saved\n");

    /* For each parameter */
    for(param=state->currentParams->head; param; param=param->next)
    {
        if(param->type!=LALINFERENCE_REAL8_t || !(param->vary==LALINFERENCE_PARAM_CIRCULAR ||param->vary==LALINFERENCE_PARAM_LINEAR )) continue;
        /* Create sorted parameter vector */
        REAL8Vector *sampvec=XLALCreateREAL8Vector(Nmcmc);
        for(i=0;i<NvarArray;i++)
            sampvec->data[i]=*(REAL8 *)LALInferenceGetVariable(&(varArray[i]),param->name);
        qsort((void *)(sampvec->data),NvarArray,sizeof(REAL8),cmpREAL8p);

        /* Create cumulative distribution */
        REAL8Vector *cumvec=XLALCreateREAL8Vector(Nmcmc);
        for(i=0;i<NvarArray;i++)
            cumvec->data[i]=PriorCDF(param->name,sampvec->data[i],state->priorArgs);

        /* Perform test*/
        REAL8 Pval=KSPValue(sampvec,cumvec);
        printf("%s: P-val = %lf\n",param->name,Pval);
        XLALDestroyREAL8Vector(sampvec);
        XLALDestroyREAL8Vector(cumvec);
    }
   
    FILE *lout=NULL;
    char param_list[FILENAME_MAX];
    sprintf(param_list,"%s_params.txt",filename);
    lout=fopen(param_list,"w");
    LALInferenceSortVariablesByName(state->livePoints[0]);
    for(LALInferenceVariableItem *param_ptr=state->livePoints[0]->head;param_ptr;param_ptr=param_ptr->next)
    {
        fprintf(lout,"%s\t",param_ptr->name);
    }
    fclose(lout);
    
    return(0);
}


static int cmpREAL8p(const void *p1, const void *p2)
{
    REAL8 r1=*(const REAL8 *)p1;
    REAL8 r2=*(const REAL8 *)p2;
    if(r1<r2) return -1;
    if(r1==r2) return 0;
    else return 1;
}

/* Calculatethe CDF at point x for the variable given its name and the prior args */
REAL8 PriorCDF(const char *name, const REAL8 x, LALInferenceVariables *priorArgs)
{
    REAL8 min=0,max=0;
    LALInferenceGetMinMaxPrior(priorArgs,name,&min,&max);
    if(!strcmp(name,"inclination")) return(FlatInCosine(x,min,max));
    if(!strcmp(name,"theta_spin1")) return(FlatInCosine(x,min,max));
    if(!strcmp(name,"theta_spin2")) return(FlatInCosine(x,min,max));
    if(!strcmp(name,"declination")) return(FlatInSine(x,min,max));
    if(!strcmp(name,"distance")) return(rSquaredCDF(x,min,max));
    if(!strcmp(name,"logdistance")) return(rSquaredCDF(exp(x),exp(min),exp(max)));
    if(!strcmp(name,"chirpmass")) return(mchirpminmax(x,min,max));
    else return(UniformMinMax(x,min,max));

}

REAL8 mchirpminmax(const REAL8 x, const REAL8 min, const REAL8 max)
{
    return ( (pow(x,-5./6.) -pow(min,-5./6.) )/(pow(max,-5./6.)-pow(min,-5./6.)));
}

REAL8 UniformMinMax(const REAL8 x,const REAL8 min,const REAL8 max)
{
    return ((x-min)/(max-min));
}

REAL8 FlatInCosine(const REAL8 x,const REAL8 min,const REAL8 max)
{
    return( (cos(min)-cos(x)) / (cos(min)-cos(max)) );
}

REAL8 FlatInSine(const REAL8 x, const REAL8 min, const REAL8 max)
{
    return( (sin(x)-sin(min)) / (sin(max)-sin(min)) );
}

REAL8 rSquaredCDF(const REAL8 x, const REAL8 min, const REAL8 max)
{
    return( (x*x*x - min*min*min)/(max*max*max - min*min*min));
}
