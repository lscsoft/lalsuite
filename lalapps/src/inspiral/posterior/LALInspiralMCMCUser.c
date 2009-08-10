/*  <lalVerbatim file="LALInspiralMCMCUSERCV">
Author: A. Dietz, J. Veitch, C. Roever
$Id: LALInspiralPhase.c,v 1.9 2003/04/14 00:27:22 sathya Exp $
</lalVerbatim>  */


/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralMCMCUser.c} }

The file \texttt{LALInspiralMCMCUser} contains user algorithms for the MCMC computation.
\subsubsection*{Prototypes}
\vspace{0.1in}

\subsubsection*{Description}

The only functions called from outside this package are {\tt XLALMCMCMetro}, as well as functions to set up and initialize the parameter structure.

To set a different chain for the parameter estimation, just set another random seed.

\subsubsection*{Algorithm}

The algorithms used in these functions are explained in detail in [Ref Needed].

\subsubsection*{Uses}


\subsubsection*{Notes}



</lalLaTeX>  */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include "LALInspiralMCMC.h"

#include <lal/Date.h>
#include <lal/Random.h>
#include <lal/AVFactories.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/FindChirpBCV.h>
#include <lal/FindChirpTD.h>
#include <lal/LALError.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>

#include "LALInspiralMCMCUser.h"
#include <fftw3.h>

#define rint(x) floor((x)+0.5)
#define MpcInMeters 3.08568025e22

#define DEBUGMODEL 0
gsl_rng *RNG;
double timewindow;
REAL4Vector *model;
REAL4Vector *Tmodel;
REAL8Sequence **topdown_sum;
REAL8 *normalisations;

/*NRCSID (LALINSPIRALMCMCUSERC, "$Id: LALInspiralPhase.c,v 1.9 2003/04/14 00:27:22 sathya Exp $"); */

double mc2mass1(double mc, double eta)
/* mass 1 (the smaller one) for given mass ratio & chirp mass */
{
 double root = sqrt(0.25-eta);
 double fraction = (0.5+root) / (0.5-root);
 return mc * (pow(1+fraction,0.2) / pow(fraction,0.6));
}


double mc2mass2(double mc, double eta)
/* mass 2 (the greater one) for given mass ratio & chirp mass */
{
 double root = sqrt(0.25-eta);
 double inversefraction = (0.5-root) / (0.5+root);
 return mc * (pow(1+inversefraction,0.2) / pow(inversefraction,0.6));
}

double mc2mt(double mc, double eta)
/* total mass (mt) for given mass ratio & chirp mass */
{
 double root = sqrt(0.25-eta);
 double fraction = (0.5+root) / (0.5-root);
 double inversefraction = (0.5-root) / (0.5+root);
 return mc * ((pow(1+fraction,0.2) / pow(fraction,0.6))
              + (pow(1+inversefraction,0.2) / pow(inversefraction,0.6)));
}

double m2eta(double m1, double m2)
/* component masses to eta */
{
	return(m1*m2/((m1+m2)*(m1+m2)));
}

double m2mc(double m1, double m2)
/* component masses to chirp mass */
{
	return(pow(m2eta(m1,m2),0.6)*(m1+m2));
}

double logJacobianMcEta(double mc, double eta)
/* posterior multiplier for transformed parameters */
/* (jacobian of inverse tranformation)             */
/* (mc & eta  instead of  m1 & m2)                 */
{
 double result;
 double m1, m2, msum, mprod;
 double term1, term2, term3, term4;
 double a,b,c,d;
 /* Factor for the  (m1,m2) --> (mc,eta)  transform: */
 msum  = mc2mt(mc,eta) * LAL_MSUN_SI;
 m1    = mc2mass1(mc, eta) * LAL_MSUN_SI;
 m2    = msum-m1;
 mprod = m1*m2;
 term1 = 0.6*pow(msum,-0.2);
 term2 = 0.2*pow(msum,-1.2)*pow(mprod,0.6);
 term3 = pow(msum,-2.0);
 term4 = 2*mprod*pow(msum,-3.0);
 a = pow(m1,-0.4)*pow(m2,0.6)*term1 - term2;
 b = m2*term3-term4;
 c = pow(m1,0.6)*pow(m2,-0.4)*term1 - term2;
 d = m1*term3-term4;
 result =  -log(fabs(a*d - b*c));
 return result;
}


/*  <lalVerbatim file="MCMCInitTestCP"> */
void MCMCInitTest(
  LALMCMCParameter  *parameter,
  SnglInspiralTable *inspiralTable)
{/* </lalVerbatim> */

  /* create the parameters and set the boundaries */
  parameter->param = NULL;
  parameter->dimension=0;
  XLALMCMCAddParam( parameter, "x",  0.2, -5.0, 5.0, 0);
  XLALMCMCAddParam( parameter, "y",  0.5, -5.0, 5.0, 0);
}


/*  <lalVerbatim file="MCMCLikelihoodTestCP"> */
REAL8 MCMCLikelihoodTest(
    LALMCMCInput      *inputMCMC,
    LALMCMCParameter  *parameter)
{/* </lalVerbatim> */

  double x,y;
  double a,b,d;

  x   = XLALMCMCGetParameter( parameter, "x" );
  y   = XLALMCMCGetParameter( parameter, "y" );

  a = exp(-0.5*x*x/(1.0*1.0));
  b = exp(-0.5*y*y/(2.0*2.0));

  return a*b;
}


/*  <lalVerbatim file="MCMCPriorTestCP"> */
INT4 MCMCPriorTest(
  REAL4            *logPrior,
  LALMCMCInput     *inputMCMC,
  LALMCMCParameter *parameter)
{/* </lalVerbatim> */

  /* always return prior of 1.0 */

  *logPrior = 0.0;
  return 1;
}

int ParamInRange(LALMCMCParameter *parameter)
{
int i;
int inrange=1;
LALMCMCParam *p=parameter->param;
for(i=0;i<parameter->dimension;i++){
	inrange &= (p->value<=p->core->maxVal)&&(p->value>=p->core->minVal);
	p=p->next;
}
if(!inrange) parameter->logPrior = -DBL_MAX;
return inrange;
}

void NestInitInjNINJA(LALMCMCParameter *parameter, void *iT){
REAL8 time,mcmin,mcmax;
SimInspiralTable *injTable = (SimInspiralTable *)iT;
REAL4 mtot,eta,mwindow,localetawin;
parameter->param = NULL;
parameter->dimension = 0;
time = (REAL8) injTable->geocent_end_time.gpsSeconds + (REAL8)injTable->geocent_end_time.gpsNanoSeconds *1.0e-9;
REAL8 mchirp,total;

/*eta = injTable->eta;*/

/*double etamin = eta-0.5*etawindow;
etamin = etamin<0.01?0.01:etamin;*/
double etamin=0.01;
/*double etamax = eta+0.5*etawindow;
etamax = etamax>0.25?0.25:etamax;*/
double etamax=0.25;
localetawin=etamax-etamin;
mcmin=m2mc(25.,25.);
mcmax=m2mc(75.,75.);

do{
	mchirp=mcmin+(mcmax-mcmin)*gsl_rng_uniform(RNG);
	eta=etamin+localetawin*gsl_rng_uniform(RNG);
	total=mc2mt(mchirp,eta);
	}
while(total>475);

/*		parameter structure, name of parameter, initial value of parameter, minimum value parameter, maximum value of parameter, wrapped?) */
XLALMCMCAddParam(parameter,"mchirp",mchirp,mcmin,mcmax,0);
/*XLALMCMCAddParam(parameter,"mtotal",gsl_rng_uniform(RNG)*100.0+50.0,50.0,150.0,0);*/
/*XLALMCMCAddParam(parameter,"mtotal",3.0+27.0*gsl_rng_uniform(RNG),3.0,30.0,0);*/
XLALMCMCAddParam(parameter, "eta", eta , etamin, etamax, 0);
XLALMCMCAddParam(parameter, "time",		(gsl_rng_uniform(RNG)-0.5)*timewindow + time ,time-0.5*timewindow,time+0.5*timewindow,0);
XLALMCMCAddParam(parameter, "phi",		LAL_TWOPI*gsl_rng_uniform(RNG),0.0,LAL_TWOPI,1);
XLALMCMCAddParam(parameter, "distMpc", 499.0*gsl_rng_uniform(RNG)+1.0, 1.0, 500.0, 0);
XLALMCMCAddParam(parameter,"long",LAL_TWOPI*gsl_rng_uniform(RNG),0,LAL_TWOPI,1);
XLALMCMCAddParam(parameter,"lat",LAL_PI*(gsl_rng_uniform(RNG)-0.5),-LAL_PI/2.0,LAL_PI/2.0,0);
XLALMCMCAddParam(parameter,"psi",0.5*LAL_PI*gsl_rng_uniform(RNG),0,LAL_PI/2.0,0);
XLALMCMCAddParam(parameter,"iota",LAL_PI*gsl_rng_uniform(RNG),0,LAL_PI,0);


return;
}

void NestInitInjNINJAHighMass(LALMCMCParameter *parameter, void *iT){
REAL8 time,mcmin,mcmax;
SimInspiralTable *injTable = (SimInspiralTable *)iT;
REAL4 mtot,eta,mwindow,localetawin;
parameter->param = NULL;
parameter->dimension = 0;
time = (REAL8) injTable->geocent_end_time.gpsSeconds + (REAL8)injTable->geocent_end_time.gpsNanoSeconds *1.0e-9;
REAL8 mchirp,total;

/*eta = injTable->eta;*/

/*double etamin = eta-0.5*etawindow;
etamin = etamin<0.01?0.01:etamin;*/
double etamin=0.01;
/*double etamax = eta+0.5*etawindow;
etamax = etamax>0.25?0.25:etamax;*/
double etamax=0.25;
localetawin=etamax-etamin;
mcmin=m2mc(25.,25.);
mcmax=m2mc(237.5,237.5);

do{
	mchirp=mcmin+(mcmax-mcmin)*gsl_rng_uniform(RNG);
	eta=etamin+localetawin*gsl_rng_uniform(RNG);
	total=mc2mt(mchirp,eta);
	}
while(total>475.0);

/*		parameter structure, name of parameter, initial value of parameter, minimum value parameter, maximum value of parameter, wrapped?) */
XLALMCMCAddParam(parameter,"mchirp",mchirp,mcmin,mcmax,0);
/*XLALMCMCAddParam(parameter,"mtotal",gsl_rng_uniform(RNG)*100.0+50.0,50.0,150.0,0);*/
/*XLALMCMCAddParam(parameter,"mtotal",3.0+27.0*gsl_rng_uniform(RNG),3.0,30.0,0);*/
XLALMCMCAddParam(parameter, "eta", eta , etamin, etamax, 0);
XLALMCMCAddParam(parameter, "time",		(gsl_rng_uniform(RNG)-0.5)*timewindow + time ,time-0.5*timewindow,time+0.5*timewindow,0);
XLALMCMCAddParam(parameter, "phi",		LAL_TWOPI*gsl_rng_uniform(RNG),0.0,LAL_TWOPI,1);
XLALMCMCAddParam(parameter, "distMpc", 499.0*gsl_rng_uniform(RNG)+1.0, 1.0, 500.0, 0);
XLALMCMCAddParam(parameter,"long",LAL_TWOPI*gsl_rng_uniform(RNG),0,LAL_TWOPI,1);
XLALMCMCAddParam(parameter,"lat",LAL_PI*(gsl_rng_uniform(RNG)-0.5),-LAL_PI/2.0,LAL_PI/2.0,0);
XLALMCMCAddParam(parameter,"psi",0.5*LAL_PI*gsl_rng_uniform(RNG),0,LAL_PI/2.0,0);
XLALMCMCAddParam(parameter,"iota",LAL_PI*gsl_rng_uniform(RNG),0,LAL_PI,0);


return;
}

REAL8 GRBPrior(LALMCMCInput *inputMCMC,LALMCMCParameter *parameter)
{
  REAL8 mNS,mComp;
  REAL8 mc,eta;
  /* Priors for the GRB component masses */
#define m1min 1.0
#define m1max 3.0
#define m2min 1.0
#define m2max 35.0
  if(XLALMCMCCheckParameter(parameter,"logM")) mc=exp(XLALMCMCGetParameter(parameter,"logM"));
  else mc=XLALMCMCGetParameter(parameter,"mchirp");
  eta=XLALMCMCGetParameter(parameter,"eta");
  parameter->logPrior+=log(fabs(cos(XLALMCMCGetParameter(parameter,"lat"))));
  parameter->logPrior+=log(fabs(sin(XLALMCMCGetParameter(parameter,"iota"))));
  parameter->logPrior+=logJacobianMcEta(mc,eta);
  parameter->logPrior-=2.0*log(XLALMCMCGetParameter(parameter,"distMpc"));
  ParamInRange(parameter);
  /*check GRB component masses */
  mNS=mc2mass1(mc,eta);
  mComp=mc2mass2(mc,eta);
  if(mNS<m1min || mNS>m1max || mComp<m2min || mComp>m2max) parameter->logPrior=-DBL_MAX;
  return(parameter->logPrior);

}

REAL8 NestPriorHighMass(LALMCMCInput *inputMCMC,LALMCMCParameter *parameter)
{
  REAL8 m1,m2,logdl,ampli,a=50,b=21;
  parameter->logPrior=0.0;
  REAL8 mc,eta;
  REAL8 minCompMass = 1.0;
  REAL8 maxCompMass = 100.0;

  /* Check in range */
  if(XLALMCMCCheckParameter(parameter,"logM")) mc=exp(XLALMCMCGetParameter(parameter,"logM"));
  else mc=XLALMCMCGetParameter(parameter,"mchirp");

  eta=XLALMCMCGetParameter(parameter,"eta");
  m1 = mc2mass1(mc,eta);
  m2 = mc2mass2(mc,eta);

  parameter->logPrior+=log(fabs(cos(XLALMCMCGetParameter(parameter,"lat"))));
  parameter->logPrior+=log(fabs(sin(XLALMCMCGetParameter(parameter,"iota"))));
  /*      parameter->logPrior+=logJacobianMcEta(mc,eta);*/
  ParamInRange(parameter);
  if(inputMCMC->approximant==IMRPhenomA && mc2mt(mc,eta)>475.0) parameter->logPrior=-DBL_MAX;
  if(m1<minCompMass || m2<minCompMass) parameter->logPrior=-DBL_MAX;
  if(m1>maxCompMass || m2>maxCompMass) parameter->logPrior=-DBL_MAX;
  return parameter->logPrior;
}

REAL8 NestPrior(LALMCMCInput *inputMCMC,LALMCMCParameter *parameter)
{
	REAL8 m1,m2,logdl,ampli,a=50,b=21;
	parameter->logPrior=0.0;
	REAL8 mc,eta;
	REAL8 minCompMass = 1.0;
	REAL8 maxCompMass = 35.0;
	/* copied from alex's function */
/*	logdl=2.0*XLALMCMCGetParameter(parameter,"distMpc");
	parameter->logPrior+=2.0*logdl;
	m1=XLALMCMCGetParameter(parameter,"mass1");
	m2=XLALMCMCGetParameter(parameter,"mass2");
	ampli = log(sqrt(m1*m2)/(logdl*pow(m1+m2,1.0/6.0)));
    parameter->logPrior+= -log( 1.0+exp((ampli-a)/b) );
*/
/* Check in range */
	if(XLALMCMCCheckParameter(parameter,"logM")) mc=exp(XLALMCMCGetParameter(parameter,"logM"));
	else mc=XLALMCMCGetParameter(parameter,"mchirp");

	eta=XLALMCMCGetParameter(parameter,"eta");
	m1 = mc2mass1(mc,eta);
	m2 = mc2mass2(mc,eta);

	parameter->logPrior+=log(fabs(cos(XLALMCMCGetParameter(parameter,"lat"))));
	parameter->logPrior+=log(fabs(sin(XLALMCMCGetParameter(parameter,"iota"))));
	/*	parameter->logPrior+=logJacobianMcEta(mc,eta);*/
	ParamInRange(parameter);
	if(inputMCMC->approximant==IMRPhenomA && mc2mt(mc,eta)>475.0) parameter->logPrior=-DBL_MAX;
	if(m1<minCompMass || m2<minCompMass) parameter->logPrior=-DBL_MAX;
	if(m1>maxCompMass || m2>maxCompMass) parameter->logPrior=-DBL_MAX;
	return parameter->logPrior;
}

REAL8 MCMCLikelihood1IFO(LALMCMCInput *inputMCMC,LALMCMCParameter *parameter,int idx)
{
/* Function to generate a inspiral wave using LALInspiralWave, then calculate the log
   likelihood for a single interferometer */

/* WARNING! THIS ASSUMES THAT THE DATA IN inputMCMC->stilde IS NOT MULTIPLIED BY THE
   PSD, i.e. IT IS JUST THE WINDOWED FFT OF THE TIME SERIES. We probably want to create a new field
   in the inputMCMC structure to hold this data to avoid confusion */

	static LALStatus status;
	InspiralTemplate template;
	UINT4 Nmodel; /* Length of the model */
	INT4 i;
	int Fdomain;
/*	memset(&status,0,sizeof(LALStatus)); */
	REAL8 chisq=0.0;
	REAL8 real,imag,f,t;
	REAL8 Mchirp;
	REAL8 eta;
	Fdomain = 1;
	memset(&template,0,sizeof(InspiralTemplate));
/* Populate the template */
	if(XLALMCMCCheckParameter(parameter,"logM")) Mchirp=exp(XLALMCMCGetParameter(parameter,"logM"));
	else Mchirp = XLALMCMCGetParameter(parameter,"mchirp");
	eta = XLALMCMCGetParameter(parameter,"eta");
	template.totalMass = mc2mt(Mchirp,eta);
	template.eta = eta;
	template.massChoice = totalMassAndEta;
	template.fLower = inputMCMC->fLow;
	template.distance = XLALMCMCGetParameter(parameter,"distMpc");
	template.order=LAL_PNORDER_TWO;
	template.approximant=inputMCMC->approximant;
	template.tSampling = 1.0/inputMCMC->deltaT;
	template.fCutoff = 0.5/inputMCMC->deltaT -1.0;
	template.nStartPad = 0;
	template.nEndPad =0;
	template.startPhase = XLALMCMCGetParameter(parameter,"phi");
	template.startTime = XLALMCMCGetParameter(parameter,"time");
	template.startTime -= inputMCMC->stilde[idx]->epoch.gpsSeconds + 1e-9*inputMCMC->stilde[idx]->epoch.gpsNanoSeconds;
	template.ieta = 1;
	template.next = NULL;
	template.fine = NULL;

/*fprintf(stderr,"XLALErrNoGlobal=%i\n",*XLALGetErrnoPtr());
XLALClearErrno();*/

/* Is this the correct way to set the end time? */
/*	XLALGPSSetREAL8((&template.end_time),XLALMCMCGetParameter(parameter,"time"));*/
/* Fill the rest of the mass/tc parameters in */
	LALInspiralParameterCalc(&status,&template);

	template.startTime-=template.tC;

	LALInspiralRestrictedAmplitude(&status,&template);

/* Find the number of points needed to store the template - disabled */
	/*LALInspiralWaveLength(&status,&Nmodel,template);*/
/* Set the number of points in the template from the lengh of the data */

	if(Fdomain) Nmodel=inputMCMC->stilde[idx]->data->length*2; /* *2 for real/imag packing format */
	else	Nmodel = inputMCMC->segment[idx]->data->length;

	if(model==NULL)	LALCreateVector(&status,&model,Nmodel);

	LALInspiralWave(&status,model,&template); /* Create the model */

	if(inputMCMC->verbose) fprintf(stderr,"duration: %fs\tfinal freq: %f Hz \n",template.tC,template.fFinal);
/*	FILE *fp=fopen("waveform.txt","w");
	for(i=0;i<model->length;i++) fprintf(fp,"%e %e\n",i/(Nmodel/template.tSampling),model->data[i]);
	fclose(fp);
	fp=fopen("data.txt","w");
	for(i=0;i<inputMCMC->stilde[idx]->data->length;i++) fprintf(fp,"%e %e %e\n",i*inputMCMC->stilde[idx]->deltaF,inputMCMC->stilde[idx]->data->data[i].re,inputMCMC->stilde[idx]->data->data[i].im);
	fclose(fp);
*/
	if(Fdomain){
	/* model is a f-domain vector with complex representation = model[j] + i*model[Nmodel-j] */
	int lowBin = (int)(inputMCMC->fLow / inputMCMC->stilde[idx]->deltaF);
	/* Compute the -logL */
	for(i=lowBin;i<Nmodel/2;i++){
		real=inputMCMC->stilde[idx]->data->data[i].re - (REAL8) model->data[i];
		imag=inputMCMC->stilde[idx]->data->data[i].im - (REAL8) model->data[Nmodel-i];
		chisq+=(real*real + imag*imag)*inputMCMC->invspec[idx]->data->data[i];
/*		printf("%e\t%e\t%e\t%e\n",inputMCMC->invspec[idx]->data->data[i],real,imag,chisq); */
	}
	chisq*=2.0; /* for 2 sigma^2 on denominator */
/*	LALDestroyVector(&status,&model); */ /* Free the model vector */
	/* add the normalisation constant */
	for(i=lowBin;i<inputMCMC->invspec[idx]->data->length-1;i++) chisq+=0.5*log(inputMCMC->invspec[idx]->data->data[i]);
	chisq+=(REAL8)( 0.5 * (inputMCMC->invspec[idx]->data->length-lowBin) * log(2.0*LAL_PI));
	} /* end F domain part */
	else {
	/* T domain part yet to be written */


	}


	return(-chisq);
}

REAL8 MCMCLikelihoodMultiCoherentF(LALMCMCInput *inputMCMC,LALMCMCParameter *parameter)
/* Calculate the likelihood of the signal using multiple interferometer data sets,
in the frequency domain */
{
	REAL8 logL=0.0;
	UINT4 det_i;
	CHAR name[10] = "inspiral";
	REAL8 TimeFromGC; /* Time delay from geocentre */
	static LALStatus status;
	REAL4FFTPlan *likelihoodPlan=NULL;
	REAL8 resp_r,resp_i,ci;
	InspiralTemplate template;
	UINT4 Nmodel; /* Length of the model */
	UINT4 idx;
	INT4 i,NtimeModel;
	LALDetAMResponse det_resp;
	int Fdomain;
	REAL8 chisq=0.0;
	REAL8 real,imag,f,t;
	TofVIn TofVparams;
	FILE *modelout;
	memset(&template,0,sizeof(InspiralTemplate));
/* Populate the template */
	REAL8 ChirpISCOLength;
	REAL8 eta,mtot,mchirp;
	expnFunc expnFunction;
	expnCoeffs ak;
	if(XLALMCMCCheckParameter(parameter,"logM")) mchirp=exp(XLALMCMCGetParameter(parameter,"logM"));
        else mchirp=XLALMCMCGetParameter(parameter,"mchirp");

	eta = XLALMCMCGetParameter(parameter,"eta");
	mtot=mc2mt(mchirp,eta);
	template.totalMass = mtot;
	template.eta = eta;
	template.massChoice = totalMassAndEta;
	template.fLower = inputMCMC->fLow;
	template.distance = XLALMCMCGetParameter(parameter,"distMpc"); /* This must be in Mpc, contrary to the docs */
	template.order=LAL_PNORDER_TWO;
	template.approximant=inputMCMC->approximant;
	template.tSampling = 1.0/inputMCMC->deltaT;
	template.fCutoff = 0.5/inputMCMC->deltaT -1.0;
	template.nStartPad = 0;
	template.nEndPad =0;
	template.startPhase = XLALMCMCGetParameter(parameter,"phi");
	template.startTime = 0.0;
	template.ieta = 1;
	template.next = NULL;
	template.fine = NULL;

	/* Compute frequency-domain waveform in free space */
	LALInspiralParameterCalc(&status,&template);
	LALInspiralRestrictedAmplitude(&status,&template);
	Nmodel = (inputMCMC->stilde[0]->data->length-1)*2; /* *2 for real/imag packing format */

	if(model==NULL)	LALCreateVector(&status,&model,Nmodel); /* Allocate storage for the waveform */

	/* Perform additional setup for time domain model */
	if(template.approximant==TaylorT2 || template.approximant==TaylorT3){
		NtimeModel = inputMCMC->segment[0]->data->length;
		if(Tmodel==NULL) LALCreateVector(&status,&Tmodel,NtimeModel);
		LALInspiralWave(&status,Tmodel,&template);
		float winNorm = sqrt(inputMCMC->window->sumofsquares/inputMCMC->window->data->length);
		float Norm = winNorm * inputMCMC->deltaT;
		for(idx=0;idx<Tmodel->length;idx++) Tmodel->data[idx]*=(REAL4)inputMCMC->window->data->data[idx] * Norm; /* window & normalise */

		if(likelihoodPlan==NULL) {LALCreateForwardREAL4FFTPlan(&status,&likelihoodPlan,(UINT4) NtimeModel,FFTW_PATIENT); fprintf(stderr,"Created FFTW plan\n");}
		LALREAL4VectorFFT(&status,model,Tmodel,likelihoodPlan); /* REAL4VectorFFT doesn't normalise like TimeFreqRealFFT, so we do this above in Norm */


		/*LALDestroyREAL4FFTPlan(&status,&plan);*/
	}
	else{
		if(template.approximant==IMRPhenomA) {
			template.distance*=LAL_PC_SI*1.0e6; /* PhenomA takes distance in metres */
			LALBBHPhenWaveFreqDom(&status,model,&template);
		}
		else LALInspiralWave(&status,model,&template); /* Create the model waveform - StationaryPhaseApprox2 includes deltaF factor */
	}
	memset(&ak,0,sizeof(expnCoeffs));
	memset(&TofVparams,0,sizeof(TofVparams));
	/* Calculate the time of ISCO (v = 6^(-1/2) ) */
	LALInspiralSetup(&status,&ak,&template);
	LALInspiralChooseModel(&status,&expnFunction,&ak,&template);
	TofVparams.coeffs=&ak;
	TofVparams.dEnergy=expnFunction.dEnergy;
	TofVparams.flux=expnFunction.flux;
	TofVparams.v0= ak.v0;
	TofVparams.t0= ak.t0;
	TofVparams.vlso= ak.vlso;
	TofVparams.totalmass=ak.totalmass;
/*	LALInspiralTofV(&status,&ChirpISCOLength,pow(6.0,-0.5),(void *)&TofVparams);*/
	ChirpISCOLength=ak.tn;

	/* This is the time of the start of the wave in the GeoCentre */
	REAL8 TimeShiftToGC=XLALMCMCGetParameter(parameter,"time");

	TimeShiftToGC-=inputMCMC->epoch.gpsSeconds + 1.e-9*inputMCMC->epoch.gpsNanoSeconds;
/*	fprintf(stderr,"time from epoch to end of wave %lf\n",TimeShiftToGC);*/
/*	TimeShiftToGC-=template.tC;*/
	TimeShiftToGC-=ChirpISCOLength;
/*	template.nStartPad = (INT4)(TimeShiftToGC/inputMCMC->deltaT);*/
/*	fprintf(stderr,"ChirpLengthISCO = %lf\n",ChirpISCOLength);*/

	/* Initialise structures for detector response calcs */
	LALSource source; /* The position and polarisation of the binary */
	source.equatorialCoords.longitude = XLALMCMCGetParameter(parameter,"long");
	source.equatorialCoords.latitude = XLALMCMCGetParameter(parameter,"lat");
	source.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
	source.orientation = XLALMCMCGetParameter(parameter,"psi");
/*	source.name = (CHAR *)NULL; */

	ci = cos(XLALMCMCGetParameter(parameter,"iota")); /* cos iota */

	LALPlaceAndGPS det_gps; /* This will hold the detector site and epoch of observation */
	det_gps.p_gps=&(inputMCMC->epoch);
	DetTimeAndASource DTAAS; /* This holds the source and the detector */
	DTAAS.p_source = &(source.equatorialCoords);
	DTAAS.p_det_and_time=&det_gps;

	/* This also holds the source and the detector, LAL has two different structs for this! */
	LALDetAndSource det_source;
	det_source.pSource=&source;

	/* The epoch of observation and the accuracy required ( we don't care about a few leap seconds) */
	LALGPSandAcc GPSandAcc;
	memcpy(&(GPSandAcc.gps),&(inputMCMC->epoch),sizeof(LIGOTimeGPS));
	GPSandAcc.accuracy = LALLEAPSEC_LOOSE; /* Don't need to worry about leap seconds for AM response func */
	for(det_i=0;det_i<inputMCMC->numberDataStreams;det_i++){ /* For each detector */
		#if DEBUGMODEL !=0
			char modelname[100];
			sprintf(modelname,"model_%i.dat",det_i);
			modelout = fopen(modelname,"w");
		#endif
		chisq=0.0;
		/* Compute time delay */
		det_gps.p_detector = (inputMCMC->detector[det_i]); /* Select detector */
		LALTimeDelayFromEarthCenter(&status,&TimeFromGC,&DTAAS); /* Compute time delay */
		REAL8 time_sin;
		REAL8 time_cos;

		/* Compute detector amplitude response */
		det_source.pDetector = (inputMCMC->detector[det_i]); /* select detector */
		LALComputeDetAMResponse(&status,&det_resp,&det_source,&GPSandAcc); /* Compute det_resp */
		det_resp.plus*=0.5*(1.0+ci*ci);
		det_resp.cross*=ci;
		/* Compute the response to the wave in the detector */
		REAL8 deltaF = inputMCMC->stilde[det_i]->deltaF;
		int lowBin = (int)(inputMCMC->fLow / inputMCMC->stilde[det_i]->deltaF);
		int highBin = (int)(template.fFinal / inputMCMC->stilde[det_i]->deltaF);
		if(highBin==0 || highBin>inputMCMC->stilde[det_i]->data->length-1) highBin=inputMCMC->stilde[det_i]->data->length-1;
		
		for(idx=lowBin;idx<=highBin;idx++){
			time_sin = sin(LAL_TWOPI*(TimeFromGC+TimeShiftToGC)*((double) idx)*deltaF);
			time_cos = cos(LAL_TWOPI*(TimeFromGC+TimeShiftToGC)*((double) idx)*deltaF);

/* Version derived 19/08/08 */
			REAL8 hc = (REAL8)model->data[idx]*time_cos + (REAL8)model->data[Nmodel-idx]*time_sin;
			REAL8 hs = (REAL8)model->data[Nmodel-idx]*time_cos - (REAL8)model->data[idx]*time_sin;
			resp_r = det_resp.plus * hc - det_resp.cross * hs;
			resp_i = det_resp.cross * hc + det_resp.plus * hs;

			real=inputMCMC->stilde[det_i]->data->data[idx].re - resp_r/deltaF;
			imag=inputMCMC->stilde[det_i]->data->data[idx].im - resp_i/deltaF;


/* Gaussian version */
/* NOTE: The factor deltaF is to make ratio dimensionless, when using the specific definitions of the vectors
that LAL uses. Please check this whenever any change is made */
			chisq+=(real*real + imag*imag)*inputMCMC->invspec[det_i]->data->data[idx];

/* Student-t version */
/*			chisq+=log(real*real+imag*imag); */
			#if DEBUGMODEL !=0
				fprintf(modelout,"%lf %10.10e %10.10e\n",i*deltaF,resp_r,resp_i);
			#endif
		}
		#if DEBUGMODEL !=0
			fclose(modelout);
		#endif
		if(highBin<inputMCMC->stilde[det_i]->data->length-2 && highBin>lowBin) chisq+=topdown_sum[det_i]->data[highBin+1];
		else if(highBin<=lowBin) chisq+=topdown_sum[det_i]->data[highBin+1];
		chisq*=2.0*deltaF; /* for 2 sigma^2 on denominator, also in student-t version */
		/* add the normalisation constant */

		/*		chisq+=normalisations[det_i]; */ /*Gaussian version*/
		/*chisq+=(Nmodel-lowBin)*log(2);*/ /* student-t version */

		/*chisq+=(REAL8)( 0.5 * (inputMCMC->invspec[det_i]->data->length-lowBin) * log(2.0*LAL_PI));*/
		logL-=chisq;
	}

	/* Add log likelihoods to find global likelihood */
	parameter->logLikelihood=logL;
	return(logL);
}

REAL8 MCMCSTLikelihoodMultiCoherentF(LALMCMCInput *inputMCMC,LALMCMCParameter *parameter)
/* Calculate the likelihood of the signal using multiple interferometer data sets,
in the frequency domain */
{
	REAL8 logL=0.0;
	int det_i;
	CHAR name[10] = "inspiral";
	REAL8 TimeFromGC; /* Time delay from geocentre */
	static LALStatus status;
	REAL8 resp_r,resp_i,ci;
	InspiralTemplate template;
	UINT4 Nmodel; /* Length of the model */
	INT4 i,NtimeModel;
	LALDetAMResponse det_resp;
	REAL4FFTPlan *likelihoodPlan=NULL;
	int Fdomain;
	REAL8 chisq=0.0;
	REAL8 real,imag,f,t;
	TofVIn TofVparams;
	FILE *modelout;
	memset(&template,0,sizeof(InspiralTemplate));
/* Populate the template */
	REAL8 ChirpISCOLength;
	REAL8 eta,mtot,mchirp;
	expnFunc expnFunction;
	expnCoeffs ak;
	if(XLALMCMCCheckParameter(parameter,"logM")) mchirp=exp(XLALMCMCGetParameter(parameter,"logM"));
	else mchirp=XLALMCMCGetParameter(parameter,"mchirp");
	eta = XLALMCMCGetParameter(parameter,"eta");
	mtot=mc2mt(mchirp,eta);
	template.totalMass = mtot;
	template.eta = eta;
	template.massChoice = totalMassAndEta;
	template.fLower = inputMCMC->fLow;
	template.distance = XLALMCMCGetParameter(parameter,"distMpc"); /* This must be in Mpc, contrary to the docs */
	template.order=LAL_PNORDER_TWO;
	template.approximant=inputMCMC->approximant;
	template.tSampling = 1.0/inputMCMC->deltaT;
	template.fCutoff = 0.5/inputMCMC->deltaT -1.0;
	template.nStartPad = 0;
	template.nEndPad =0;
	template.startPhase = XLALMCMCGetParameter(parameter,"phi");
	template.startTime = 0.0;
	template.ieta = 1;
	template.next = NULL;
	template.fine = NULL;

	/* Compute frequency-domain waveform in free space */
	LALInspiralParameterCalc(&status,&template);
	LALInspiralRestrictedAmplitude(&status,&template);
	Nmodel = (inputMCMC->stilde[0]->data->length-1)*2; /* *2 for real/imag packing format */

	if(model==NULL)	LALCreateVector(&status,&model,Nmodel); /* Allocate storage for the waveform */

	/* Perform additional setup for time domain model */
	if(template.approximant==TaylorT2 || template.approximant==TaylorT3){
		NtimeModel = inputMCMC->segment[0]->data->length;
		if(Tmodel==NULL) LALCreateVector(&status,&Tmodel,NtimeModel);
		LALInspiralWave(&status,Tmodel,&template);
		float winNorm = sqrt(inputMCMC->window->sumofsquares/inputMCMC->window->data->length);
		float Norm = winNorm * inputMCMC->deltaT;
		for(i=0;i<Tmodel->length;i++) Tmodel->data[i]*=(REAL4)inputMCMC->window->data->data[i] * Norm; /* window & normalise */

		if(likelihoodPlan==NULL) {LALCreateForwardREAL4FFTPlan(&status,&likelihoodPlan,(UINT4) NtimeModel,FFTW_PATIENT); fprintf(stderr,"Created FFTW plan\n");}
		LALREAL4VectorFFT(&status,model,Tmodel,likelihoodPlan); /* REAL4VectorFFT doesn't normalise like TimeFreqRealFFT, so we do this above in Norm */


		/*LALDestroyREAL4FFTPlan(&status,&plan);*/
	}
	else{
		if(template.approximant==IMRPhenomA) {
			template.distance*=LAL_PC_SI*1.0e6; /* PhenomA takes distance in metres */
			LALBBHPhenWaveFreqDom(&status,model,&template);
		}
		else LALInspiralWave(&status,model,&template); /* Create the model waveform - StationaryPhaseApprox2 includes deltaF factor */
	}
	memset(&ak,0,sizeof(expnCoeffs));
	memset(&TofVparams,0,sizeof(TofVparams));
	/* Calculate the time of ISCO (v = 6^(-1/2) ) */
	LALInspiralSetup(&status,&ak,&template);
	LALInspiralChooseModel(&status,&expnFunction,&ak,&template);
	TofVparams.coeffs=&ak;
	TofVparams.dEnergy=expnFunction.dEnergy;
	TofVparams.flux=expnFunction.flux;
	TofVparams.v0= ak.v0;
	TofVparams.t0= ak.t0;
	TofVparams.vlso= ak.vlso;
	TofVparams.totalmass=ak.totalmass;
/*	LALInspiralTofV(&status,&ChirpISCOLength,pow(6.0,-0.5),(void *)&TofVparams);*/
	ChirpISCOLength=ak.tn;

	/* This is the time of the start of the wave in the GeoCentre */
	REAL8 TimeShiftToGC=XLALMCMCGetParameter(parameter,"time");

	TimeShiftToGC-=inputMCMC->epoch.gpsSeconds + 1.e-9*inputMCMC->epoch.gpsNanoSeconds;
/*	fprintf(stderr,"time from epoch to end of wave %lf\n",TimeShiftToGC);*/
/*	TimeShiftToGC-=template.tC;*/
	TimeShiftToGC-=ChirpISCOLength;
/*	template.nStartPad = (INT4)(TimeShiftToGC/inputMCMC->deltaT);*/
/*	fprintf(stderr,"ChirpLengthISCO = %lf\n",ChirpISCOLength);*/

	/* Initialise structures for detector response calcs */
	LALSource source; /* The position and polarisation of the binary */
	source.equatorialCoords.longitude = XLALMCMCGetParameter(parameter,"long");
	source.equatorialCoords.latitude = XLALMCMCGetParameter(parameter,"lat");
	source.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
	source.orientation = XLALMCMCGetParameter(parameter,"psi");
/*	source.name = (CHAR *)NULL; */

	ci = cos(XLALMCMCGetParameter(parameter,"iota")); /* cos iota */

	LALPlaceAndGPS det_gps; /* This will hold the detector site and epoch of observation */
	det_gps.p_gps=&(inputMCMC->epoch);
	DetTimeAndASource DTAAS; /* This holds the source and the detector */
	DTAAS.p_source = &(source.equatorialCoords);
	DTAAS.p_det_and_time=&det_gps;

	/* This also holds the source and the detector, LAL has two different structs for this! */
	LALDetAndSource det_source;
	det_source.pSource=&source;

	/* The epoch of observation and the accuracy required ( we don't care about a few leap seconds) */
	LALGPSandAcc GPSandAcc;
	memcpy(&(GPSandAcc.gps),&(inputMCMC->epoch),sizeof(LIGOTimeGPS));
	GPSandAcc.accuracy = LALLEAPSEC_LOOSE; /* Don't need to worry about leap seconds for AM response func */
	for(det_i=0;det_i<inputMCMC->numberDataStreams;det_i++){ /* For each detector */
		#if DEBUGMODEL !=0
			char modelname[100];
			sprintf(modelname,"model_%i.dat",det_i);
			modelout = fopen(modelname,"w");
		#endif
		chisq=0.0;
		/* Compute time delay */
		det_gps.p_detector = (inputMCMC->detector[det_i]); /* Select detector */
		LALTimeDelayFromEarthCenter(&status,&TimeFromGC,&DTAAS); /* Compute time delay */
		REAL8 time_sin;
		REAL8 time_cos;

		/* Compute detector amplitude response */
		det_source.pDetector = (inputMCMC->detector[det_i]); /* select detector */
		LALComputeDetAMResponse(&status,&det_resp,&det_source,&GPSandAcc); /* Compute det_resp */
		det_resp.plus*=0.5*(1.0+ci*ci);
		det_resp.cross*=ci;
		/* Compute the response to the wave in the detector */
		REAL8 deltaF = inputMCMC->stilde[det_i]->deltaF;
		int lowBin = (int)(inputMCMC->fLow / inputMCMC->stilde[det_i]->deltaF);
		int highBin = (int)(template.fFinal / inputMCMC->stilde[det_i]->deltaF);

		for(i=lowBin;i<Nmodel/2;i++){
			time_sin = sin(LAL_TWOPI*(TimeFromGC+TimeShiftToGC)*((double) i)*deltaF);
			time_cos = cos(LAL_TWOPI*(TimeFromGC+TimeShiftToGC)*((double) i)*deltaF);

/* Version derived 19/08/08 */
			REAL8 hc = (REAL8)model->data[i]*time_cos + (REAL8)model->data[Nmodel-i]*time_sin;
			REAL8 hs = (REAL8)model->data[Nmodel-i]*time_cos - (REAL8)model->data[i]*time_sin;
			resp_r = det_resp.plus * hc - det_resp.cross * hs;
			resp_i = det_resp.cross * hc + det_resp.plus * hs;

			real=inputMCMC->stilde[det_i]->data->data[i].re - resp_r/deltaF;
			imag=inputMCMC->stilde[det_i]->data->data[i].im - resp_i/deltaF;


/* Gaussian version */
/* NOTE: The factor deltaF is to make ratio dimensionless, when using the specific definitions of the vectors
that LAL uses. Please check this whenever any change is made */
			/* chisq+=(real*real + imag*imag)*inputMCMC->invspec[det_i]->data->data[i]*inputMCMC->invspec[det_i]->deltaF; */

/* Student-t version */
			chisq+=log(deltaF*(real*real+imag*imag));
			#if DEBUGMODEL !=0
				fprintf(modelout,"%lf %10.10e %10.10e\n",i*deltaF,resp_r,resp_i);
			#endif
		}
		#if DEBUGMODEL !=0
			fclose(modelout);
		#endif
		/*		chisq+=topdown_sum[det_i]->data[highBin+1];*/
		chisq*=2.0; /* for 2 sigma^2 on denominator, also in student-t version */
		/* add the normalisation constant */

		/*chisq+=normalisations[det_i];*/  /*Gaussian version*/
		/*		chisq+=(Nmodel-lowBin)*log(2); *//* student-t version */

		/*chisq+=(REAL8)( 0.5 * (inputMCMC->invspec[det_i]->data->length-lowBin) * log(2.0*LAL_PI));*/
		logL-=chisq;
	}

	/* Add log likelihoods to find global likelihood */
	parameter->logLikelihood=logL;
	return(logL);
}



REAL8 MCMCLikelihoodMultiCoherent(LALMCMCInput *inputMCMC,LALMCMCParameter *parameter)
{
	REAL8 logL=0.0;
	int det_i,i;
	static LALStatus status;
	PPNParamStruc PPNparams;
	CoherentGW co_wave; /* output structure */
	DetectorResponse det_resp;
	REAL4TimeSeries *template;
	REAL8TimeSeries *cur_seg;
	COMPLEX8FrequencySeries *Ftemplate;
	REAL4FFTPlan *FFTplan;
	REAL8 mchirp=0;
	REAL8 eta=0;
        if(XLALMCMCCheckParameter(parameter,"logM")) mchirp=exp(XLALMCMCGetParameter(parameter,"logM"));
	else mchirp=XLALMCMCGetParameter(parameter,"mchirp");
	eta = XLALMCMCGetParameter(parameter,"eta");


	/* Populate input struct */
	PPNparams.mTot=(REAL4)mc2mt(mchirp,eta);
	PPNparams.eta=(REAL4)eta;
	PPNparams.d=(REAL4)XLALMCMCGetParameter(parameter,"distMpc")*MpcInMeters;
	PPNparams.phi=(REAL4)XLALMCMCGetParameter(parameter,"phi");
	PPNparams.position.system=COORDINATESYSTEM_EQUATORIAL;
	PPNparams.position.longitude=XLALMCMCGetParameter(parameter,"long");
	PPNparams.position.latitude=XLALMCMCGetParameter(parameter,"lat");
	PPNparams.psi=(REAL4)XLALMCMCGetParameter(parameter,"psi");
	PPNparams.inc=(REAL4)XLALMCMCGetParameter(parameter,"iota");
	PPNparams.deltaT=(REAL8)inputMCMC->deltaT;
	PPNparams.fStartIn=(REAL4)inputMCMC->fLow;
	PPNparams.fStopIn=(REAL4)0.5/inputMCMC->deltaT - 1.0;
	PPNparams.lengthIn=0; /* Let the waveform be as long as necessary */
	PPNparams.ppn=NULL;

	memcpy(&(PPNparams.epoch),&(inputMCMC->epoch),sizeof(LIGOTimeGPS));

	memset(&co_wave,0,sizeof(CoherentGW)); /* Blank the output fields */
	memset(&det_resp,0,sizeof(DetectorResponse)); /* Only need the detector field in this struct, the rest are left blank */

	/* Create coherent waveform */
	LALGeneratePPNInspiral(&status,&co_wave,&PPNparams);


	/* Adjust time */
	REAL8 tC=XLALMCMCGetParameter(parameter,"time");
	tC-=PPNparams.tc; /* tC is now the time the wave reaches the low frequency */
	XLALGPSSetREAL8(&(co_wave.a->epoch),tC); /* which is the starting epoch of the co_wave */
	memcpy(&(co_wave.f->epoch),&(co_wave.a->epoch),sizeof(LIGOTimeGPS));
	memcpy(&(co_wave.phi->epoch),&(co_wave.a->epoch),sizeof(LIGOTimeGPS));

	FFTplan=XLALCreateForwardREAL4FFTPlan(inputMCMC->segment[0]->data->length,0);
	/* Loop over detectors */
	for(det_i=0;det_i<inputMCMC->numberDataStreams;det_i++)
	{
		cur_seg=inputMCMC->segment[det_i];
		/* Allocate template */
		template=(REAL4TimeSeries *)XLALCreateREAL4TimeSeries("template",&(cur_seg->epoch),cur_seg->f0,cur_seg->deltaT,&(cur_seg->sampleUnits),(size_t)cur_seg->data->length);
		/* compute waveform in detector det_i */
		det_resp.site=inputMCMC->detector[det_i];
		LALSimulateCoherentGW(&status,template,&co_wave,&det_resp);

		/* Calculate likelihood in that detector */
/*		char filename[100];
		sprintf(filename,"wave_%s.dat",inputMCMC->ifoID[det_i]);
		FILE *fpout=fopen(filename,"w");
		for(i=0;i<template->data->length;i++) fprintf(fpout,"%e\t%10.10e\n",i*template->deltaT,template->data->data[i]);
		fclose(fpout);
*/
		for(i=0;i<template->data->length;i++) template->data->data[i]*=(REAL4)inputMCMC->window->data->data[i];
		Ftemplate=(COMPLEX8FrequencySeries *)XLALCreateCOMPLEX8FrequencySeries("Ftemplate",
			&(cur_seg->epoch),
			cur_seg->f0,
			1.0/cur_seg->deltaT,
			&(inputMCMC->stilde[det_i]->sampleUnits),
			(size_t)inputMCMC->stilde[det_i]->data->length);

		XLALREAL4TimeFreqFFT(Ftemplate,template,FFTplan);
/*
		sprintf(filename,"waveF_%sdat",inputMCMC->ifoID[det_i]);
		fpout=fopen(filename,"w");
		int N=inputMCMC->stilde[det_i]->data->length;
		for(i=0;i<Ftemplate->data->length;i++) fprintf(fpout,"%e\t%10.10e\t%10.10e\t%10.10e\t%10.10e\n",i*Ftemplate->deltaF,Ftemplate->data->data[i].re,Ftemplate->data->data[i].im,inputMCMC->stilde[det_i]->data->data[i].re,inputMCMC->stilde[det_i]->data->data[i].im);
		fclose(fpout);
*/
		int lowBin = (int)(inputMCMC->fLow / inputMCMC->stilde[det_i]->deltaF);
		/* Compute the -logL */
		REAL8 chisq=0.0;
		for(i=lowBin;i<Ftemplate->data->length;i++){
			REAL8 real=inputMCMC->stilde[det_i]->data->data[i].re - (REAL8) Ftemplate->data->data[i].re;
			REAL8 imag=inputMCMC->stilde[det_i]->data->data[i].im - (REAL8) Ftemplate->data->data[i].im;
			chisq+=(real*real + imag*imag)*inputMCMC->invspec[det_i]->data->data[i];
		}
		chisq*=0.5; /* for 2 sigma^2 on denominator */
		/* add the normalisation constant */
		for(i=lowBin;i<inputMCMC->invspec[det_i]->data->length-1;i++) chisq+=0.5*log(inputMCMC->invspec[det_i]->data->data[i]);
		chisq+=(REAL8)( 0.5 * (inputMCMC->invspec[det_i]->data->length-lowBin) * log(2.0*LAL_PI));

		/* Free template */
		XLALDestroyREAL4TimeSeries(template);
		XLALDestroyCOMPLEX8FrequencySeries(Ftemplate);
		logL-=chisq;
	}
	/* Free co_wave structures */

	XLALDestroyVectorSequence((co_wave.a->data));
	LALFree(co_wave.a);
	XLALDestroyREAL8TimeSeries(co_wave.phi);
	XLALDestroyREAL4TimeSeries(co_wave.f);
	LALDestroyREAL4FFTPlan(&status,&FFTplan);
	parameter->logLikelihood=logL;
	return(logL);
}


REAL8 MCMCLikelihoodMultiIFO(LALMCMCInput *inputMCMC,LALMCMCParameter *parameter)
{
	REAL8 logL=0.0;
	int i;
	/* Loop over detectors/datasets */
	for(i=0;i<inputMCMC->numberDataStreams;i++) logL+= MCMCLikelihood1IFO(inputMCMC,parameter,i);
	parameter->logLikelihood=logL;
	return logL;
}

void MCMCInit1IFO(LALMCMCParameter *parameter, SnglInspiralTable *inspiralTable)
{
REAL4 mass1,mass2;
REAL8 time,length;

mass1=inspiralTable->mass1;
mass2=inspiralTable->mass2;
time = (REAL8)inspiralTable->end_time.gpsSeconds + (REAL8)inspiralTable->end_time.gpsNanoSeconds * 1.0e-9;
length=(REAL8)inspiralTable->template_duration;
time=-length;
	parameter->param=NULL;
	parameter->dimension=0;
	XLALMCMCAddParam(parameter,"mass1", mass1, 0.98*mass1, 1.02*mass1,0);
	XLALMCMCAddParam(parameter,"mass2", mass2, 0.98*mass2, 1.02*mass2,0);
	XLALMCMCAddParam(parameter,"time",time,time-0.01,time+0.01,0);
	XLALMCMCAddParam(parameter,"phi",0.0,0.0,LAL_TWOPI,1);
	XLALMCMCAddParam(parameter,"logdistance",10.0,1.0,90.0,0);

}
