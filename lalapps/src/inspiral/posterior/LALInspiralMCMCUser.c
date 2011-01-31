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
#include <lal/Units.h>
#include <lal/TimeSeries.h>
#include <lal/VectorOps.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/SeqFactories.h>

#include "LALInspiralMCMCUser.h"
#include <fftw3.h>

#define MpcInMeters 3.08568025e22

#define DEBUGMODEL 0
#if DEBUGMODEL
FILE *modelout=NULL;
#endif
gsl_rng *RNG;
double timewindow;
REAL4Vector *model;
REAL4Vector *Tmodel;
REAL8Sequence **topdown_sum;
REAL8 *normalisations;
const LALUnit strainPerCount={0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};

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
  /* FIXME: this function doesn't use the inspiralTable argument, the
   * correct fix is to modify this function not to take inspiralTable as
   * an argument, this simply supresses the unused parameter warning */
  (void)inspiralTable;

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
  double a,b;

  /* FIXME: this function doesn't use the inputMCMC argument, the
   * correct fix is to modify this function not to take inputMCMC as
   * an argument, this simply supresses the unused parameter warning */
  (void)inputMCMC;

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
  (void)inputMCMC;
  (void)parameter;

  /* always return prior of 1.0 */

  *logPrior = 0.0;
  return 1;
}

int ParamInRange(LALMCMCParameter *parameter)
{
UINT4 i;
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
REAL8 trg_time,mcmin,mcmax;
SimInspiralTable *injTable = (SimInspiralTable *)iT;
REAL4 eta,localetawin;
parameter->param = NULL;
parameter->dimension = 0;
trg_time = (REAL8) injTable->geocent_end_time.gpsSeconds + (REAL8)injTable->geocent_end_time.gpsNanoSeconds *1.0e-9;
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
XLALMCMCAddParam(parameter, "time",		(gsl_rng_uniform(RNG)-0.5)*timewindow + trg_time ,trg_time-0.5*timewindow,trg_time+0.5*timewindow,0);
XLALMCMCAddParam(parameter, "phi",		LAL_TWOPI*gsl_rng_uniform(RNG),0.0,LAL_TWOPI,1);
XLALMCMCAddParam(parameter, "distMpc", 499.0*gsl_rng_uniform(RNG)+1.0, 1.0, 500.0, 0);
XLALMCMCAddParam(parameter,"long",LAL_TWOPI*gsl_rng_uniform(RNG),0,LAL_TWOPI,1);
XLALMCMCAddParam(parameter,"lat",LAL_PI*(gsl_rng_uniform(RNG)-0.5),-LAL_PI/2.0,LAL_PI/2.0,0);
XLALMCMCAddParam(parameter,"psi",0.5*LAL_PI*gsl_rng_uniform(RNG),0,LAL_PI/2.0,0);
XLALMCMCAddParam(parameter,"iota",LAL_PI*gsl_rng_uniform(RNG),0,LAL_PI,0);


return;
}

void NestInitInjNINJAHighMass(LALMCMCParameter *parameter, void *iT){
REAL8 trg_time,mcmin,mcmax;
SimInspiralTable *injTable = (SimInspiralTable *)iT;
REAL4 eta,localetawin;
parameter->param = NULL;
parameter->dimension = 0;
trg_time = (REAL8) injTable->geocent_end_time.gpsSeconds + (REAL8)injTable->geocent_end_time.gpsNanoSeconds *1.0e-9;
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
XLALMCMCAddParam(parameter, "time",		(gsl_rng_uniform(RNG)-0.5)*timewindow + trg_time ,trg_time-0.5*timewindow,trg_time+0.5*timewindow,0);
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
  REAL8 mNS,mComp,logmc;
  REAL8 mc,eta;

  /* FIXME: this function doesn't use the inputMCMC argument, the
   * correct fix is to modify this function not to take inputMCMC as
   * an argument, this simply supresses the unused parameter warning */
  (void)inputMCMC;

  /* Priors for the GRB component masses */
#define m1min 1.0
#define m1max 3.0
#define m2min 1.0
#define m2max 35.0
  parameter->logPrior=0.0;
  if(XLALMCMCCheckParameter(parameter,"logM")) mc=exp(XLALMCMCGetParameter(parameter,"logM"));
  else mc=XLALMCMCGetParameter(parameter,"mchirp");
  logmc=log(mc);
  parameter->logPrior+=-(5.0/6.0)*logmc;

  eta=XLALMCMCGetParameter(parameter,"eta");
  parameter->logPrior+=log(fabs(cos(XLALMCMCGetParameter(parameter,"lat"))));
  parameter->logPrior+=log(fabs(sin(XLALMCMCGetParameter(parameter,"iota"))));
  /*parameter->logPrior+=logJacobianMcEta(mc,eta);*/
  parameter->logPrior+=2.0*log(XLALMCMCGetParameter(parameter,"distMpc"));
  ParamInRange(parameter);
  /*check GRB component masses */
  mNS=mc2mass1(mc,eta);
  mComp=mc2mass2(mc,eta);
  if(mNS<m1min || mNS>m1max || mComp<m2min || mComp>m2max) parameter->logPrior=-DBL_MAX;
  return(parameter->logPrior);

}

REAL8 NestPriorHighMass(LALMCMCInput *inputMCMC,LALMCMCParameter *parameter)
{
  REAL8 m1,m2;
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
	parameter->logPrior+=-(5.0/6.0)*log(mc);
  if(XLALMCMCCheckParameter(parameter,"logdist"))
     parameter->logPrior+=3.0*XLALMCMCGetParameter(parameter,"logdist");
  else
     parameter->logPrior+=2.0*log(XLALMCMCGetParameter(parameter,"distMpc"));

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
	REAL8 m1,m2;
	parameter->logPrior=0.0;
	REAL8 mc,eta;
	REAL8 minCompMass = 1.0;
	REAL8 maxCompMass = 34.0;
#define MAX_MTOT 35.0
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
	double logmc=log(mc);
	eta=XLALMCMCGetParameter(parameter,"eta");
	m1 = mc2mass1(mc,eta);
	m2 = mc2mass2(mc,eta);
	/* This term is the sqrt of m-m term in F.I.M, ignoring dependency on f and eta */
	parameter->logPrior+=-(5.0/6.0)*logmc;
	if(XLALMCMCCheckParameter(parameter,"logdist"))
		parameter->logPrior+=3.0*XLALMCMCGetParameter(parameter,"logdist");
	else
		parameter->logPrior+=2.0*log(XLALMCMCGetParameter(parameter,"distMpc"));
	parameter->logPrior+=log(fabs(cos(XLALMCMCGetParameter(parameter,"lat"))));
	parameter->logPrior+=log(fabs(sin(XLALMCMCGetParameter(parameter,"iota"))));
	/*	parameter->logPrior+=logJacobianMcEta(mc,eta);*/
	ParamInRange(parameter);
	if(inputMCMC->approximant==IMRPhenomA && mc2mt(mc,eta)>475.0) parameter->logPrior=-DBL_MAX;
	if(m1<minCompMass || m2<minCompMass) parameter->logPrior=-DBL_MAX;
	if(m1>maxCompMass || m2>maxCompMass) parameter->logPrior=-DBL_MAX;
	if(m1+m2>MAX_MTOT) parameter->logPrior=-DBL_MAX;
	return parameter->logPrior;
}


REAL8 NestPriorPhenSpin(LALMCMCInput *inputMCMC,LALMCMCParameter *parameter)
{

    (void)inputMCMC;
    REAL8 m1,m2;
        parameter->logPrior=0.0;
        REAL8 mc,eta;
        REAL8 minCompMass = 1.0;
        REAL8 maxCompMass = 34.0;
    REAL8 maxMTotal= 68.0;

/* Check in range */
        if(XLALMCMCCheckParameter(parameter,"logMc")) mc=exp(XLALMCMCGetParameter(parameter,"logMc"));
        else mc=XLALMCMCGetParameter(parameter,"mchirp");
        double logmc=log(mc);
        eta=XLALMCMCGetParameter(parameter,"eta");
        m1 = mc2mass1(mc,eta);
        m2 = mc2mass2(mc,eta);
        /* This term is the sqrt of m-m term in F.I.M, ignoring dependency on f and eta */
        parameter->logPrior+=-(5.0/6.0)*logmc;
        if(XLALMCMCCheckParameter(parameter,"logdist"))
                parameter->logPrior+=3.0*XLALMCMCGetParameter(parameter,"logdist");
        else
                parameter->logPrior+=2.0*log(XLALMCMCGetParameter(parameter,"distMpc"));
        parameter->logPrior+=log(fabs(cos(XLALMCMCGetParameter(parameter,"lat"))));
        parameter->logPrior+=log(fabs(sin(XLALMCMCGetParameter(parameter,"iota"))));
        ParamInRange(parameter);

        if(m1<minCompMass || m2<minCompMass) {
          parameter->logPrior=-DBL_MAX;
        }
        if(m1>maxCompMass || m2>maxCompMass) {
          parameter->logPrior=-DBL_MAX;
        }
        if(m1+m2>maxMTotal) {
          parameter->logPrior=-DBL_MAX;
        }
        return parameter->logPrior;
}


REAL8 MCMCLikelihoodMultiCoherentAmpCor(LALMCMCInput *inputMCMC, LALMCMCParameter *parameter){
	/* Calculate the likelihood for an amplitude-corrected waveform */
	/* This template is generated in the time domain */
	REAL8 logL=0.0,chisq=0.0;
	REAL8 mc,eta,end_time,resp_r,resp_i,real,imag;
	UINT4 det_i=0,idx=0;
	UINT4 i=0;
	DetectorResponse det;
	static LALStatus status;
	CoherentGW coherent_gw;
	PPNParamStruc PPNparams;
	LALDetAMResponse det_resp;
	REAL4TimeSeries *h_p_t=NULL,*h_c_t=NULL;
	COMPLEX8FrequencySeries *H_p_t=NULL, *H_c_t=NULL;
	size_t NFD = 0;
	memset(&PPNparams,0,sizeof(PPNparams));
	memset(&coherent_gw,0,sizeof(CoherentGW));
	memset(&status,0,sizeof(LALStatus));
	memset(&det,0,sizeof(DetectorResponse));
	/* Populate the structures */
	if(XLALMCMCCheckParameter(parameter,"logM")) mc=exp(XLALMCMCGetParameter(parameter,"logM"));
	else mc=XLALMCMCGetParameter(parameter,"mchirp");
	eta=XLALMCMCGetParameter(parameter,"eta");
	PPNparams.position.longitude=XLALMCMCGetParameter(parameter,"long");
	PPNparams.position.latitude=XLALMCMCGetParameter(parameter,"lat");
	PPNparams.position.system=COORDINATESYSTEM_EQUATORIAL;
	PPNparams.psi=XLALMCMCGetParameter(parameter,"psi");
	memcpy(&(PPNparams.epoch),&(inputMCMC->epoch),sizeof(LIGOTimeGPS));
	PPNparams.mTot=mc2mt(mc,eta);
	PPNparams.eta=eta;
	if (XLALMCMCCheckParameter(parameter,"logdist")) PPNparams.d=exp(XLALMCMCGetParameter(parameter,"logdist"))*MpcInMeters;
	else PPNparams.d=XLALMCMCGetParameter(parameter,"distMpc")*MpcInMeters;
	PPNparams.inc=XLALMCMCGetParameter(parameter,"iota");
	PPNparams.phi=XLALMCMCGetParameter(parameter,"phi");
	PPNparams.fStartIn=inputMCMC->fLow;
	PPNparams.fStopIn=0.5/inputMCMC->deltaT;
	PPNparams.deltaT=inputMCMC->deltaT;
	PPNparams.ampOrder = inputMCMC->ampOrder;

	/* Call LALGeneratePPNAmpCorInspiral */
	LALGeneratePPNAmpCorInspiral(&status,&coherent_gw,&PPNparams);
	if(status.statusCode)
	{
		REPORTSTATUS(&status);
		chisq=DBL_MAX;
		goto noWaveform;
	}

	/* Set the epoch so that the t_c is correct */
	end_time = XLALMCMCGetParameter(parameter,"time");

	REAL8 adj_epoch = end_time-PPNparams.tc;
	if(coherent_gw.h) XLALGPSSetREAL8(&(coherent_gw.h->epoch),adj_epoch);
	if(coherent_gw.a) XLALGPSSetREAL8(&(coherent_gw.a->epoch),adj_epoch);

	/* Inject h+ and hx into time domain signal of correct length */
	UINT4 NtimeDomain=inputMCMC->segment[det_i]?inputMCMC->segment[det_i]->data->length:2*(inputMCMC->stilde[det_i]->data->length-1);
	h_p_t = XLALCreateREAL4TimeSeries("hplus",&inputMCMC->epoch,
										   inputMCMC->fLow,inputMCMC->deltaT,
										   &lalADCCountUnit,NtimeDomain);
	h_c_t = XLALCreateREAL4TimeSeries("hcross",&inputMCMC->epoch,
										 inputMCMC->fLow,inputMCMC->deltaT,
										 &lalADCCountUnit,NtimeDomain);
	if(!(h_p_t && h_c_t)){
		fprintf(stderr,"Unable to allocate signal buffer\n");
		exit(1);
	}

	/* Separate the + and x parts */
	for(i=0;i< (NtimeDomain<coherent_gw.h->data->length?NtimeDomain:coherent_gw.h->data->length) ;i++){
		h_p_t->data->data[i]=coherent_gw.h->data->data[2*i];
		h_c_t->data->data[i]=coherent_gw.h->data->data[2*i + 1];
	}
	for(;i<NtimeDomain;i++){
		h_p_t->data->data[i]=0.0;
		h_c_t->data->data[i]=0.0;
	}

	/* Get H+ and Hx in the Freq Domain */
	NFD=inputMCMC->stilde[det_i]->data->length;
	H_p_t = XLALCreateCOMPLEX8FrequencySeries("Hplus",&inputMCMC->epoch,0,(REAL8)inputMCMC->deltaF,&lalDimensionlessUnit,(size_t)NFD);
	H_c_t = XLALCreateCOMPLEX8FrequencySeries("Hcross",&inputMCMC->epoch,0,(REAL8)inputMCMC->deltaF,&lalDimensionlessUnit,(size_t)NFD);

	if(!(H_p_t && H_c_t)){
		fprintf(stderr,"Unable to allocate F-domain signal buffer\n");
		exit(1);
	}

	if(inputMCMC->likelihoodPlan==NULL) {LALCreateForwardREAL4FFTPlan(&status,&inputMCMC->likelihoodPlan,
																	  NtimeDomain,
																	  FFTW_PATIENT); fprintf(stderr,"Created FFTW plan\n");}
	LALTimeFreqRealFFT(&status,H_p_t,h_p_t,inputMCMC->likelihoodPlan);
	LALTimeFreqRealFFT(&status,H_c_t,h_c_t,inputMCMC->likelihoodPlan);

#if DEBUGMODEL !=0
	char tmodelname[100];
	sprintf(tmodelname,"tmodel_plus_%i.dat",det_i);
	modelout = fopen(tmodelname,"w");
	for(i=0;i<h_p_t->data->length;i++){
		fprintf(modelout,"%g %g\n",h_p_t->data->data[i],h_c_t->data->data[i]);
	}
	fclose(modelout);
#endif
	XLALDestroyREAL4TimeSeries(h_p_t);
	XLALDestroyREAL4TimeSeries(h_c_t);

	/* The epoch of observation and the accuracy required ( we don't care about a few leap seconds) */
	LALSource source; /* The position and polarisation of the binary */
	source.equatorialCoords.longitude = XLALMCMCGetParameter(parameter,"long");
	source.equatorialCoords.latitude = XLALMCMCGetParameter(parameter,"lat");
	source.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
	source.orientation = XLALMCMCGetParameter(parameter,"psi");

	/* This also holds the source and the detector, LAL has two different structs for this! */
	LALDetAndSource det_source;
	det_source.pSource=&source;

	REAL8 TimeShiftToGC=XLALMCMCGetParameter(parameter,"time");

	TimeShiftToGC-=inputMCMC->epoch.gpsSeconds + 1.e-9*inputMCMC->epoch.gpsNanoSeconds;
	TimeShiftToGC-=PPNparams.tc;


	/* For each IFO */
	for(det_i=0;det_i<inputMCMC->numberDataStreams;det_i++){
		REAL8 TimeFromGC;
		/* Set up the detector */
		det.site=inputMCMC->detector[det_i];
		/* Simulate the response */
#if DEBUGMODEL !=0
		char modelname[100];
		sprintf(modelname,"model_%i.dat",det_i);
		modelout = fopen(modelname,"w");
#endif

		TimeFromGC = XLALTimeDelayFromEarthCenter(inputMCMC->detector[det_i]->location, source.equatorialCoords.longitude, source.equatorialCoords.latitude, &(inputMCMC->epoch)); /* Compute time delay */
		/* Compute detector amplitude response */
		det_source.pDetector = (inputMCMC->detector[det_i]); /* select detector */
		LALComputeDetAMResponse(&status,&det_resp,&det_source,&(inputMCMC->epoch)); /* Compute det_resp */
		/* No need to multiply by cos(iota) as GenerateAmpCorPPNInspiral() takes this into account */
		chisq=0.0;
		/* Calculate the logL */
		REAL8 deltaF = inputMCMC->stilde[det_i]->deltaF;
		UINT4 lowBin = (UINT4)(inputMCMC->fLow / inputMCMC->stilde[det_i]->deltaF);
		UINT4 highBin;
		/* Only compute sum up to maximum frquency of the waveform. PPNparams.fStop is the maximum of the 2*f_orb harmonic */
		
		REAL8 fMultiplier = (inputMCMC->ampOrder + 2.0)/2.0; /* The frequency of the highest harmonic as determined by ampOrder */
		highBin = (UINT4)(PPNparams.fStop * fMultiplier / inputMCMC->stilde[det_i]->deltaF);
		
		if(highBin==0 || highBin>inputMCMC->stilde[det_i]->data->length-1)
			highBin=inputMCMC->stilde[det_i]->data->length-1;  /* AmpCor waveforms don't set the highest frequency of the highest harmonic */
		for(idx=lowBin;idx<=highBin;idx++){
			/* The phase shift angle, determined by the time FROM the geocentre to the Detector, plus the time */
			/* FROM the start of the segment to the t_c at geocentre */
			/* Phase shift is exp(-i*ang), but - signs on sin below take the negative into account */
			/* exp(-i*ang) = cos(ang) - sin(ang) */
			REAL8 ang = 2.0*LAL_PI*(TimeFromGC+TimeShiftToGC)*inputMCMC->stilde[det_i]->deltaF*idx;
			/* Calculate rotated parts of the plus and cross */
			
			/* Negative signs on sins: see comment above for definition of ang */
			REAL4 plus_re,plus_im,cross_re,cross_im;
			plus_re = H_p_t->data->data[idx].re*cos(ang) + H_p_t->data->data[idx].im*sin(ang);
			plus_im = H_p_t->data->data[idx].im*cos(ang) - H_p_t->data->data[idx].re*sin(ang);
			cross_re = H_c_t->data->data[idx].re*cos(ang) + H_c_t->data->data[idx].im*sin(ang);
			cross_im = H_c_t->data->data[idx].im*cos(ang) - H_c_t->data->data[idx].re*sin(ang);
			
			/* Compute total real and imaginary responses */
			resp_r = (REAL8)( plus_re*det_resp.plus + cross_re*det_resp.cross );
			resp_i = (REAL8)( plus_im*det_resp.plus + cross_im*det_resp.cross );
			real=inputMCMC->stilde[det_i]->data->data[idx].re - resp_r;
			imag=inputMCMC->stilde[det_i]->data->data[idx].im - resp_i;
			
			/* Gaussian version */
			chisq+=(real*real + imag*imag)*inputMCMC->invspec[det_i]->data->data[idx];
			
#if DEBUGMODEL !=0
			fprintf(modelout,"%lf %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e\n",idx*deltaF,resp_r,resp_i,H_p_t->data->data[idx].re,H_p_t->data->data[idx].im,H_c_t->data->data[idx].re,H_c_t->data->data[idx].im);
#endif
		}
#if DEBUGMODEL !=0
		fclose(modelout);
#endif
		/* Add on the remaining sum, consulting the lookup table */
		if(highBin<inputMCMC->stilde[det_i]->data->length-2 && highBin>lowBin) chisq+=topdown_sum[det_i]->data[highBin+1];
		else if(highBin<=lowBin) chisq+=topdown_sum[det_i]->data[highBin+1];
		chisq*=2.0*deltaF; /* for 2 sigma^2 on denominator, also in student-t version */
		
		logL-=chisq;
		
	}
	/* Destroy the response series */
	if(coherent_gw.f) XLALDestroyREAL4TimeSeries(coherent_gw.f);
	if(coherent_gw.phi) XLALDestroyREAL8TimeSeries(coherent_gw.phi);
	if(coherent_gw.shift) XLALDestroyREAL4TimeSeries(coherent_gw.shift);
	if(coherent_gw.h) {XLALDestroyREAL4VectorSequence(coherent_gw.h->data); LALFree(coherent_gw.h);}
	if(coherent_gw.a) {XLALDestroyREAL4VectorSequence(coherent_gw.a->data); LALFree(coherent_gw.a);}
	XLALDestroyCOMPLEX8FrequencySeries(H_p_t);
	XLALDestroyCOMPLEX8FrequencySeries(H_c_t);
	
noWaveform:
	/* return logL */
	parameter->logLikelihood=logL;
	return(logL);

}

REAL8 MCMCLikelihoodMultiCoherentF(LALMCMCInput *inputMCMC,LALMCMCParameter *parameter)
/* Calculate the likelihood of the signal using multiple interferometer data sets,
in the frequency domain */
{
	REAL8 logL=0.0;
	UINT4 det_i;
	REAL8 TimeFromGC; /* Time delay from geocentre */
	static LALStatus status;
	REAL8 resp_r,resp_i,ci;
	InspiralTemplate template;
	UINT4 Nmodel; /* Length of the model */
	UINT4 idx;
    CHAR dumpfile[2048];
	LALDetAMResponse det_resp;
	REAL8 chisq=0.0;
	REAL8 real,imag;
	TofVIn TofVparams;
	memset(&template,0,sizeof(InspiralTemplate));
/* Populate the template */
	REAL8 ChirpISCOLength;
	REAL8 eta,mtot,mchirp;
	expnFunc expnFunction;
	expnCoeffs ak;

	if(inputMCMC->numberDataStreams==0){
		parameter->logLikelihood=0.0;
		return 0.0;
	}
	if(XLALMCMCCheckParameter(parameter,"logM")) mchirp=exp(XLALMCMCGetParameter(parameter,"logM"));
        else mchirp=XLALMCMCGetParameter(parameter,"mchirp");

	eta = XLALMCMCGetParameter(parameter,"eta");
	mtot=mc2mt(mchirp,eta);
	template.totalMass = mtot;
	template.eta = eta;
	template.massChoice = totalMassAndEta;
	template.fLower = inputMCMC->fLow;
	if(XLALMCMCCheckParameter(parameter,"distMpc"))
		template.distance = XLALMCMCGetParameter(parameter,"distMpc"); /* This must be in Mpc, contrary to the docs */
	else if(XLALMCMCCheckParameter(parameter,"logdist"))
		template.distance=exp(XLALMCMCGetParameter(parameter,"logdist"));

	template.order=inputMCMC->phaseOrder;
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

	Nmodel = (inputMCMC->stilde[0]->data->length-1)*2; /* *2 for real/imag packing format */

	if(model==NULL)	LALCreateVector(&status,&model,Nmodel); /* Allocate storage for the waveform */

	switch(template.approximant)
	{
		case IMRPhenomFA : IMRPhenomFA_template(&status,&template,parameter,inputMCMC); break;
		case IMRPhenomFB : IMRPhenomFB_template(&status,&template,parameter,inputMCMC); break;
		case EOBNR : EOBNR_template(&status,&template,parameter,inputMCMC); break;
		case TaylorF2 : TaylorF2_template(&status,&template,parameter,inputMCMC); break;
		case TaylorT3 :
		case TaylorT2 :
        case TaylorT4 : TaylorT_template(&status,&template,parameter,inputMCMC); break;
		case SpinTaylor : SpinTaylor_template(&status,&template,parameter,inputMCMC); break;
		case IMRPhenomB: IMRPhenomB_template(&status,&template, parameter, inputMCMC); break;
		default: {fprintf(stderr,"This template is not available. Exiting."); exit(1);}
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

	/* IMRPhenomB calcaultes the chirp length differently */
	if(template.approximant == IMRPhenomB){
		ChirpISCOLength = template.tC;
	}


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

	/* This also holds the source and the detector, LAL has two different structs for this! */
	LALDetAndSource det_source;
	det_source.pSource=&source;

	for(det_i=0;det_i<inputMCMC->numberDataStreams;det_i++){ /* For each detector */

		chisq=0.0;
		/* Compute time delay */
		TimeFromGC = XLALTimeDelayFromEarthCenter(inputMCMC->detector[det_i]->location, source.equatorialCoords.longitude, source.equatorialCoords.latitude, &(inputMCMC->epoch)); /* Compute time delay */
		REAL8 time_sin;
		REAL8 time_cos;

		/* Compute detector amplitude response */
		det_source.pDetector = (inputMCMC->detector[det_i]); /* select detector */
		LALComputeDetAMResponse(&status,&det_resp,&det_source,&inputMCMC->epoch); /* Compute det_resp */
		det_resp.plus*=0.5*(1.0+ci*ci);
		det_resp.cross*=-ci;
		
		/* Compute the response to the wave in the detector */
		REAL8 deltaF = inputMCMC->stilde[det_i]->deltaF;
		UINT4 lowBin = (UINT4)(inputMCMC->fLow / inputMCMC->stilde[det_i]->deltaF);
		UINT4 highBin = (UINT4)(template.fFinal / inputMCMC->stilde[det_i]->deltaF);
		if(highBin==0 || highBin>inputMCMC->stilde[det_i]->data->length-1) highBin=inputMCMC->stilde[det_i]->data->length-1;

		for(idx=lowBin;idx<=highBin;idx++){
			time_sin = sin(LAL_TWOPI*(TimeFromGC+TimeShiftToGC)*((double) idx)*deltaF);
			time_cos = cos(LAL_TWOPI*(TimeFromGC+TimeShiftToGC)*((double) idx)*deltaF);

			/* Version derived 18/08/10 */
			/* This is the time delayed waveforms as it appears at the detector */
			/* data[idx] is real and data[Nmodel-idx] is imaginary part of the waveform at index idx */
			/* H+ = hc + i*hs, and Hx=iH+, ONLY WHERE H+=cos(phi) and Hx=sin(phi) in the time domain (SPA, non-spinning, no-HH) */
			
			/* Model contains h(f)exp(-psi(f)), want h'(f)=h(f)exp(-2pi*i*deltaT)  */
			/* model_re_prime and model_im_prime contain the time delayed part */
			REAL8 model_re_prime = (REAL8)model->data[idx]*time_cos + (REAL8)model->data[Nmodel-idx]*time_sin; /* Plus sign from -i*sin(phi)*i */
			REAL8 model_im_prime = (REAL8)model->data[Nmodel-idx]*time_cos - (REAL8)model->data[idx]*time_sin; /* Minus sign from -i*sin(phi) */

			/* Now, h+=model_prime and hx=i*model_prime */
			/* real(H+ + Hx) = F+(real(model_prime)) + Fx( -imag(model_prime)) : negative sign from multiplication by i*i */
			/* imag(H+ + Hx) = F+(imag(model_prime)) + Fx(real(model_prime)) : No negative sign */

			resp_r = det_resp.plus*model_re_prime - det_resp.cross*model_im_prime;
			resp_i = det_resp.plus*model_im_prime + det_resp.cross*model_re_prime;

			real=inputMCMC->stilde[det_i]->data->data[idx].re - resp_r/deltaF;
			imag=inputMCMC->stilde[det_i]->data->data[idx].im - resp_i/deltaF;
			chisq+=(real*real + imag*imag)*inputMCMC->invspec[det_i]->data->data[idx];

		} /* End loop over frequency */

		/* Dump the template and data if required */
		if(inputMCMC->dumpfile){
			sprintf(dumpfile,"%s.%s",inputMCMC->dumpfile,inputMCMC->ifoID[det_i]);
			FILE *modelout = fopen(dumpfile,"w");
			for(idx=lowBin;idx<inputMCMC->stilde[det_i]->data->length;idx++)
			{
				time_sin = sin(LAL_TWOPI*(TimeFromGC+TimeShiftToGC)*((double) idx)*deltaF);
				time_cos = cos(LAL_TWOPI*(TimeFromGC+TimeShiftToGC)*((double) idx)*deltaF);
				REAL8 model_re_prime = (REAL8)model->data[idx]*time_cos + (REAL8)model->data[Nmodel-idx]*time_sin; /* Plus sign from -i*sin(phi)*i */
				REAL8 model_im_prime = (REAL8)model->data[Nmodel-idx]*time_cos - (REAL8)model->data[idx]*time_sin; /* Minus sign from -i*sin(phi) */
				resp_r = det_resp.plus*model_re_prime - det_resp.cross*model_im_prime;
				resp_i = det_resp.plus*model_im_prime + det_resp.cross*model_re_prime;
				resp_r/=deltaF; resp_i/=deltaF;
				
				fprintf(modelout,"%4.3e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e\n",
						idx*deltaF, inputMCMC->invspec[det_i]->data->data[idx],
						inputMCMC->stilde[det_i]->data->data[idx].re, inputMCMC->stilde[det_i]->data->data[idx].im,
						resp_r, resp_i, 2.0*deltaF*(inputMCMC->stilde[det_i]->data->data[idx].re-resp_r), 2.0*deltaF*(inputMCMC->stilde[det_i]->data->data[idx].im-resp_i));
			}
			fclose(modelout);
		}
		
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

/* Just a likelihood function for PhenSpinTaylorRD_template*/


REAL8 MCMCLikelihoodMultiCoherentF_PhenSpin(LALMCMCInput *inputMCMC,LALMCMCParameter *parameter)
  
/* Calculate the likelihood of the signal using multiple interferometer data 
   sets,in the frequency domain */

{
  static LALStatus status;
  InspiralTemplate template;
  LALDetAMResponse det_resp;
  TofVIn TofVparams;
  REAL8 ChirpISCOLength;
  expnFunc expnFunction;
  expnCoeffs ak;

  REAL8 logL=0.0;
  REAL8 chisq=0.0;

  REAL8 eta,mtot=0,mchirp;
  REAL8 TimeFromGC; /* Time delay from geocentre */
  REAL8 resp_r,resp_i;
  REAL8 ci,real,imag;

  UINT4 NtimeModel;
  UINT4 Nmodel; /* Length of the model */
  UINT4 det_i;
  UINT4 idx;


  eta = XLALMCMCGetParameter(parameter,"eta");

 if (XLALMCMCCheckParameter(parameter,"logMc")) {
    mchirp=exp(XLALMCMCGetParameter(parameter,"logMc"));
    mtot=mchirp/pow(eta,3./5.);
  }
 else {
    if (XLALMCMCCheckParameter(parameter,"mchirp")) {
      mchirp=XLALMCMCGetParameter(parameter,"mchirp");
      mtot=mchirp/pow(eta,3./5.);
    }
    else {
      if (XLALMCMCCheckParameter(parameter,"mtotal")) {
	mtot=XLALMCMCGetParameter(parameter,"mtotal");
      }
    }
  }

  template.totalMass=mtot;
  template.eta=eta;
  template.massChoice=totalMassAndEta;
  template.fLower = inputMCMC->fLow;

  if(XLALMCMCCheckParameter(parameter,"distMpc"))
    template.distance = XLALMCMCGetParameter(parameter,"distMpc"); /* This must be in Mpc, contrary to the docs */
  else if(XLALMCMCCheckParameter(parameter,"logdist"))
    template.distance=exp(XLALMCMCGetParameter(parameter,"logdist"));
  
  template.order=inputMCMC->phaseOrder;
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

  inputMCMC->fLow =40.0;
				
  Nmodel = (inputMCMC->stilde[0]->data->length-1)*2; 
  /* "*2" for real/imag packing format ( see definition of stilde[0] ) */
  NtimeModel=inputMCMC->numPoints;

  if (Nmodel!=NtimeModel) {
    fprintf(stderr,"** LikelihoodPSIRD **: ERROR, T and F-domain templates have unequal length F=%d T=%d\n",Nmodel,NtimeModel);
  }

  UINT4 countplan=0;

  switch(template.approximant){
  case PhenSpinTaylorRD : 
    if (inputMCMC->likelihoodPlan==NULL) {
      countplan++;
      inputMCMC->likelihoodPlan = XLALCreateForwardREAL4FFTPlan((UINT4) NtimeModel,FFTW_ESTIMATE);
      fprintf(stdout,"** Creato plan %d volta **\n",countplan);
    }
    PhenSpinTaylorRD_template(&status,&template,parameter,inputMCMC); 
    break;
  default:
    fprintf(stderr,"Template not available. Exiting test function.\n");
    exit(0);
  }

  memset(&ak,0,sizeof(expnCoeffs));
  memset(&TofVparams,0,sizeof(TofVparams));
  
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
  TimeShiftToGC-=ChirpISCOLength;
  
  /* Initialise structures for detector response calcs */
  LALSource source; /* The position and polarisation of the binary */
  source.equatorialCoords.longitude = XLALMCMCGetParameter(parameter,"long");
  source.equatorialCoords.latitude = XLALMCMCGetParameter(parameter,"lat");
  source.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
  source.orientation = XLALMCMCGetParameter(parameter,"psi");
  /*	source.name = (CHAR *)NULL; */
	
  ci = cos(XLALMCMCGetParameter(parameter,"iota")); /* cos iota */
  
  /* This also holds the source and the detector, LAL has two different structs for this! */
  LALDetAndSource det_source;
  det_source.pSource=&source;

  for(det_i=0;det_i<inputMCMC->numberDataStreams;det_i++){ 
    /* **MAIN LOOP** For each detector */
//    char modelname[100];
//    sprintf(modelname,"%sF_%s.dat","waveform",inputMCMC->ifoID[det_i]);
    /* Compute time delay */
//    modelout = fopen(modelname,"w");
    TimeFromGC = XLALTimeDelayFromEarthCenter(inputMCMC->detector[det_i]->location, source.equatorialCoords.longitude, source.equatorialCoords.latitude, &(inputMCMC->epoch)); /* Compute time delay */
    
    /* Compute detector amplitude response */
    det_source.pDetector = (inputMCMC->detector[det_i]); /* select detector */
    LALComputeDetAMResponse(&status,&det_resp,&det_source,&inputMCMC->epoch); /* Compute det_resp */
    REAL8 deltaF = inputMCMC->stilde[det_i]->deltaF;
    UINT4 lowBin = (UINT4)(inputMCMC->fLow / inputMCMC->stilde[det_i]->deltaF);
    UINT4 highBin = (UINT4)((0.5/inputMCMC->deltaT -1.0) / inputMCMC->stilde[det_i]->deltaF);
    if(highBin==0 || highBin>inputMCMC->stilde[det_i]->data->length-1) highBin=inputMCMC->stilde[det_i]->data->length-1;
    
    /* Plus & cross polarised spinning waveforms are not related by a trivial phase shift: we calculate the comps here (see PhenSpinTaylorRD_template(..) */ 
    REAL8 hplr,hcrr,hpli,hcri; 

//        fprintf(stdout,"MCMC: Nmodel=%d  Lbin=%d  hBin=%d\n",Nmodel,lowBin,highBin);
    
//        printf("Fp=%12.5e  Fc=%12.5e\n",det_resp.plus,det_resp.cross);

    for(idx=lowBin;idx<highBin;idx++){
      
      hplr = (REAL8)inputMCMC->Fwfp->data[idx];
      hcrr = (REAL8)inputMCMC->Fwfc->data[idx];
      hpli = (REAL8)inputMCMC->Fwfp->data[Nmodel-idx];
      hcri = (REAL8)inputMCMC->Fwfc->data[Nmodel-idx];
      
      resp_r = det_resp.plus * hplr + det_resp.cross * hcrr;
      resp_i = det_resp.cross * hcri + det_resp.plus * hpli;
      
      real=inputMCMC->stilde[det_i]->data->data[idx].re - resp_r/deltaF;
      imag=inputMCMC->stilde[det_i]->data->data[idx].im - resp_i/deltaF;
      chisq+=(real*real + imag*imag)*inputMCMC->invspec[det_i]->data->data[idx];
      
//      fprintf(modelout,"%12.6f %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e\n",idx*deltaF,resp_r,resp_i, chisq, det_resp.plus, det_resp.cross,inputMCMC->detector[det_i]->location[0],inputMCMC->detector[det_i]->location[1],inputMCMC->detector[det_i]->location[2], hplr,hcrr);
    }

//    fclose(modelout);

    if(highBin<inputMCMC->stilde[det_i]->data->length-2 && highBin>lowBin) chisq+=topdown_sum[det_i]->data[highBin+1];
    else if(highBin<=lowBin) chisq+=topdown_sum[det_i]->data[highBin+1];
    chisq*=2.0*deltaF;		
    
    logL-=chisq;	
	
  }

  parameter->logLikelihood=logL;
  return(logL);

}


REAL8 MCMCSTLikelihoodMultiCoherentF(LALMCMCInput *inputMCMC,LALMCMCParameter *parameter)
/* Calculate the likelihood of the signal using multiple interferometer data sets,
in the frequency domain */
{
	REAL8 logL=0.0;
	UINT4 det_i;
	REAL8 TimeFromGC; /* Time delay from geocentre */
	static LALStatus status;
	REAL8 resp_r,resp_i,ci;
	InspiralTemplate template;
	UINT4 Nmodel; /* Length of the model */
	UINT4 i,NtimeModel;
	LALDetAMResponse det_resp;
	REAL4FFTPlan *likelihoodPlan=NULL;
	REAL8 chisq=0.0;
	REAL8 real,imag;
	TofVIn TofVparams;
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
	template.order=inputMCMC->phaseOrder;
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

	/* This also holds the source and the detector, LAL has two different structs for this! */
	LALDetAndSource det_source;
	det_source.pSource=&source;

	for(det_i=0;det_i<inputMCMC->numberDataStreams;det_i++){ /* For each detector */
		#if DEBUGMODEL !=0
			char modelname[100];
			sprintf(modelname,"model_%i.dat",det_i);
			modelout = fopen(modelname,"w");
		#endif
		chisq=0.0;
		/* Compute time delay */
		TimeFromGC = XLALTimeDelayFromEarthCenter(inputMCMC->detector[det_i]->location, source.equatorialCoords.longitude, source.equatorialCoords.latitude, &(inputMCMC->epoch)); /* Compute time delay */
		REAL8 time_sin;
		REAL8 time_cos;

		/* Compute detector amplitude response */
		det_source.pDetector = (inputMCMC->detector[det_i]); /* select detector */
		LALComputeDetAMResponse(&status,&det_resp,&det_source,&inputMCMC->epoch); /* Compute det_resp */
		det_resp.plus*=-0.5*(1.0+ci*ci);
		det_resp.cross*=-ci;
		/* Compute the response to the wave in the detector */
		REAL8 deltaF = inputMCMC->stilde[det_i]->deltaF;
		int lowBin = (int)(inputMCMC->fLow / inputMCMC->stilde[det_i]->deltaF);

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


void IMRPhenomFA_template(LALStatus *status,InspiralTemplate *template, LALMCMCParameter *parameter,LALMCMCInput *inputMCMC) {
    (void)parameter;
	/*All components of spin must be set to zero.*/
	template->spin1[0]=0.;
	template->spin1[1]=0.;
	template->spin1[2]=0.;

	template->spin2[0]=0.;
	template->spin2[1]=0.;
	template->spin2[2]=0.;

	/* Fill out parameter structure*/
	LALInspiralParameterCalc(status,template);
	LALInspiralRestrictedAmplitude(status,template);
	/* PhenomA takes distance in metres */
	double distanceMPC = template->distance;

	double distanceSI= LAL_PC_SI*1e6*distanceMPC;
	template->distance = distanceSI/inputMCMC->deltaF;//IMR doesnt multiply by df
	/*Generate IMR waveform in frequency domain*/
	LALBBHPhenWaveFreqDom(status,model,template);

	template->distance = distanceMPC;

}



void IMRPhenomFB_template(LALStatus *status,InspiralTemplate *template, LALMCMCParameter *parameter,LALMCMCInput *inputMCMC) {

	/*'x' and 'y' components of spins must be set to zero.*/
	template->spin1[0]=0.;
	template->spin1[1]=0.;
    template->spin1[2]=0.;

	template->spin2[0]=0.;
	template->spin2[1]=0.;
    template->spin2[2]=0.;

	/*Get aligned spins configuration/magnitude if these are set*/
	if(XLALMCMCCheckParameter(parameter,"spin1z")) {
        template->spin1[2] = XLALMCMCGetParameter(parameter,"spin1z");
    }

    if(XLALMCMCCheckParameter(parameter,"spin2z")) {
        template->spin2[2] = XLALMCMCGetParameter(parameter,"spin2z");
    }

    /* Fill out parameter structure*/
    LALInspiralParameterCalc(status,template);
    LALInspiralRestrictedAmplitude(status,template);

    if(XLALMCMCCheckParameter(parameter,"chiSpin")) {
        double spin1z = 1.0;
        template->spin1[2] = spin1z ;
        double chiSpin = XLALMCMCGetParameter(parameter,"chiSpin");
        double delta = sqrt(1.0- 4.0*template->eta);

        template->spin2[2] = (chiSpin*2.0 - (1+delta)*spin1z)/ (1-delta);
    }

    /* IMRPhenomFB takes distance in metres */
    double distanceMPC = template->distance;
    double distanceSI= LAL_PC_SI*1e6*distanceMPC;
    template->distance = distanceSI/(inputMCMC->deltaF*inputMCMC->deltaF);
	//IMR doesnt normalise by multiplying by df, plus TF2 has a deltaF assumed which is divided out later

    LALBBHPhenWaveFreqDom(status,model,template);
    template->distance = distanceMPC;

}

void TaylorF2_template(LALStatus *status,InspiralTemplate *template, LALMCMCParameter *parameter,LALMCMCInput *inputMCMC) {

    (void)parameter;
    (void)inputMCMC;
	/* Compute frequency-domain waveform in free space */

	LALInspiralParameterCalc(status,template);
	LALInspiralRestrictedAmplitude(status,template);

	LALInspiralWave(status,model,template);

	/*
	FILE* model_output;
	model_output=fopen("output.dat","w");

	fprintf(model_output,"Sampling frequency: %lf\n",template->tSampling);

	fprintf(model_output,"Mass 1: %lf\n",template->mass1);
	fprintf(model_output,"Mass 2: %lf\n",template->mass2);

	for(i=0;i<model->length;i++) {
		fprintf(model_output,"%g\n",model->data[i]);
	}
	fclose(model_output);

	exit(0);
	*/

	return;

}

void SpinTaylor_template(LALStatus *status,InspiralTemplate *template, LALMCMCParameter *parameter,LALMCMCInput *inputMCMC) {

	UINT4 Nmodel,idx,NtimeModel = inputMCMC->numPoints;
    if(Tmodel==NULL)	LALCreateVector(status,&Tmodel,NtimeModel); /* Allocate storage for the waveform */
		/*
	double spin_long1=XLALMCMCGetParameter(parameter,"spin_long1");
	double spin_lat1=XLALMCMCGetParameter(parameter,"spin_lat1");
	double spin_mag1=XLALMCMCGetParameter(parameter,"spin_mag1");

	double spin_long2=XLALMCMCGetParameter(parameter,"spin_long2");
	double spin_lat2=XLALMCMCGetParameter(parameter,"spin_lat2");
	double spin_mag2=XLALMCMCGetParameter(parameter,"spin_mag2");

	template->spin1[0]=cos(spin_long1)*sin(spin_lat1)*spin_mag1;
	template->spin1[1]=sin(spin_long1)*sin(spin_lat1)*spin_mag1;
	template->spin1[2]=cos(spin_lat1)*spin_mag1;

	template->spin2[0]=cos(spin_long2)*sin(spin_lat2)*spin_mag2;
	template->spin2[1]=sin(spin_long2)*sin(spin_lat2)*spin_mag2;
	template->spin2[2]=cos(spin_lat2)*spin_mag2;
*/

	template->spin1[0]=XLALMCMCGetParameter(parameter,"spin1x");
	template->spin1[1]=XLALMCMCGetParameter(parameter,"spin1y");
	template->spin1[2]=XLALMCMCGetParameter(parameter,"spin1z");

 	template->spin2[0]=XLALMCMCGetParameter(parameter,"spin2x");
	template->spin2[1]=XLALMCMCGetParameter(parameter,"spin2y");
	template->spin2[2]=XLALMCMCGetParameter(parameter,"spin2z");

	template->inclination=XLALMCMCGetParameter(parameter,"iota");


	Nmodel = (inputMCMC->stilde[0]->data->length-1)*2; /* *2 for real/imag packing format */
	NtimeModel = 1024*1024;

	if(model==NULL)	LALCreateVector(status,&model,Nmodel); /* Allocate storage for the waveform */

	if(Tmodel==NULL) LALCreateVector(status,&Tmodel,NtimeModel);

	LALInspiralParameterCalc(status,template);

	//SpinTaylor takes distance in metres
	double distance = LAL_PC_SI*1e6*template->distance;
	template->distance=distance;

	LALInspiralWave(status,Tmodel,template);

	float winNorm = sqrt(inputMCMC->window->sumofsquares/inputMCMC->window->data->length);
	float Norm = winNorm * inputMCMC->deltaT;
	for(idx=0;idx<Tmodel->length;idx++) Tmodel->data[idx]*=(REAL4)inputMCMC->window->data->data[idx] * Norm; /* window & normalise */

	if(inputMCMC->likelihoodPlan==NULL) {LALCreateForwardREAL4FFTPlan(status,&inputMCMC->likelihoodPlan,(UINT4) NtimeModel,FFTW_PATIENT);}
	LALREAL4VectorFFT(status,model,Tmodel,inputMCMC->likelihoodPlan); /* REAL4VectorFFT doesn't normalise like TimeFreqRealFFT, so we do this above in Norm */

	return;
}

void TaylorT_template(LALStatus *status,InspiralTemplate *template, LALMCMCParameter *parameter,LALMCMCInput *inputMCMC) {

    (void)parameter;

	UINT4 NtimeModel,idx;

	NtimeModel = inputMCMC->segment[0]->data->length;
	if(Tmodel==NULL) LALCreateVector(status,&Tmodel,NtimeModel);

	LALInspiralParameterCalc(status,template);
	LALInspiralRestrictedAmplitude(status,template);
	LALInspiralWave(status,Tmodel,template);

	float winNorm = sqrt(inputMCMC->window->sumofsquares/inputMCMC->window->data->length);
	float Norm = winNorm * inputMCMC->deltaT;
	for(idx=0;idx<Tmodel->length;idx++) Tmodel->data[idx]*=(REAL4)inputMCMC->window->data->data[idx] * Norm; /* window & normalise */

	if(inputMCMC->likelihoodPlan==NULL) {LALCreateForwardREAL4FFTPlan(status,&inputMCMC->likelihoodPlan,(UINT4) NtimeModel,FFTW_PATIENT);}
	LALREAL4VectorFFT(status,model,Tmodel,inputMCMC->likelihoodPlan); /* REAL4VectorFFT doesn't normalise like TimeFreqRealFFT, so we do this above in Norm */

	return;

}

void IMRPhenomB_template(LALStatus *status, InspiralTemplate *template, LALMCMCParameter *parameter, LALMCMCInput *inputMCMC)
{
	UINT4 NtimeModel, idx;
	NtimeModel = inputMCMC->segment[0]->data->length;
	if(Tmodel ==NULL) LALCreateVector(status, &Tmodel, NtimeModel);
	
        /*'x' and 'y' components of spins must be set to zero.*/
        template->spin1[0]=0.;
        template->spin1[1]=0.;
	template->spin1[2]=0.;

        template->spin2[0]=0.;
        template->spin2[1]=0.;
	template->spin2[2]=0.;

        /*Get aligned spins configuration/magnitude if these are set*/
        if(XLALMCMCCheckParameter(parameter,"spin1z")) {
        template->spin1[2] = XLALMCMCGetParameter(parameter,"spin1z");
    }

    if(XLALMCMCCheckParameter(parameter,"spin2z")) {
        template->spin2[2] = XLALMCMCGetParameter(parameter,"spin2z");
    }

    /* Fill out parameter structure*/
    LALInspiralParameterCalc(status,template);
    LALInspiralRestrictedAmplitude(status,template);

    if(XLALMCMCCheckParameter(parameter,"chiSpin")) {
        double spin1z = 1.0;
        template->spin1[2] = spin1z ;
        double chiSpin = XLALMCMCGetParameter(parameter,"chiSpin");
        double delta = sqrt(1.0- 4.0*template->eta);

        template->spin2[2] = (chiSpin*2.0 - (1+delta)*spin1z)/ (1-delta);
    }

    /* IMRPhenomFB takes distance in metres */
    double distanceMPC = template->distance;
    double distanceSI= LAL_PC_SI*1e6*distanceMPC;
    template->distance = distanceSI/(inputMCMC->deltaF*inputMCMC->deltaF);//IMR doesnt normalise by multiplying by df

    LALInspiralWave(status, Tmodel, template);
    template->distance = distanceMPC;
    
    float winNorm = sqrt(inputMCMC->window->sumofsquares/inputMCMC->window->data->length);
    float Norm = winNorm * inputMCMC->deltaT;
    for(idx=0;idx<Tmodel->length;idx++) Tmodel->data[idx]*=(REAL4)inputMCMC->window->data->data[idx] * Norm; /* window & normalise */
    if(inputMCMC->likelihoodPlan==NULL) {LALCreateForwardREAL4FFTPlan(status,&inputMCMC->likelihoodPlan,(UINT4) NtimeModel,FFTW_PATIENT);}
    LALREAL4VectorFFT(status,model,Tmodel,inputMCMC->likelihoodPlan); /* REAL4VectorFFT doesn't normalise like TimeFreqRealFFT, so we do this above in Norm */
    return;


}


void EOBNR_template(LALStatus *status,InspiralTemplate *template, LALMCMCParameter *parameter,LALMCMCInput *inputMCMC) {

    (void)parameter;

	double qnm223freq;


	UINT4 i,NtimeModel, idx;

	/*Containers for EOBNR data */
	COMPLEX8Vector *modefreqs=NULL;

	 /*Extra EOBNR parameters*/
	template->OmegaS= 0.;
	template->Zeta2=0.;
	template->Theta=1.;

	NtimeModel = (inputMCMC->stilde[0]->data->length);

	if(Tmodel==NULL) {
		LALCreateVector(status,&Tmodel,NtimeModel);
		LALCreateForwardREAL4FFTPlan(status,&inputMCMC->likelihoodPlan,(UINT4)NtimeModel,FFTW_PATIENT);
	}
	/*EOBNR requires template.order to be set to
	* pseudo-fourth order PN approx.*/
	template->order = LAL_PNORDER_PSEUDO_FOUR;


	/*EOBNR may require padding to prevent code from crashing due to
	* longer waveform.?? */
	template->nEndPad = 1.;

	/* Calculate QNM 220 frequency to determine appropriate sampling frequency. Taken from test done in EOBNR code. */
	modefreqs = XLALCreateCOMPLEX8Vector( 3 );
	XLALGenerateQNMFreq( modefreqs, template, 2, 2, 3 );

	qnm223freq = modefreqs->data[0].re / LAL_PI + 50.;

	/*Determine if sampling frequency is suitable for EOBNR template. If qnm223 freq. is greater than nyquist use EOB*/

	if( qnm223freq > template->tSampling ) {
		template->approximant=EOB;

	}

	/*Set cut-off freq. using tSampling*/
	template->fCutoff = template->tSampling/2. - 1.;

	/*Fill out parameter structure*/
	LALInspiralParameterCalc(status,template);

	/* EOBNR function uses SI units for distance so convert (from Mpc) ...*/
	double distance = LAL_PC_SI*1e6*template->distance;
	template->distance = distance;

	printf("qnm %lf samp %lf mass1 %lf mass2 %lf",qnm223freq,template->tSampling,template->mass1,template->mass2);

	/*Generate EOBNR waveform*/
	LALInspiralWave(status,Tmodel,template);


	XLALDestroyCOMPLEX8Vector(modefreqs);

	/*FFT time-domain model to get frequency domain template*/

	float winNorm = sqrt(inputMCMC->window->sumofsquares/inputMCMC->window->data->length);
	float Norm = winNorm * inputMCMC->deltaT;
	for(idx=0;idx<Tmodel->length;idx++) Tmodel->data[idx]*=(REAL4)inputMCMC->window->data->data[idx] * Norm; /* window & normalise */

	LALREAL4VectorFFT(status,model,Tmodel,inputMCMC->likelihoodPlan); /* REAL4VectorFFT doesn't normalise like TimeFreqRealFFT, so we do this above in Norm */


	FILE* model_output;
	model_output=fopen("model_output.dat","w");

	fprintf(model_output,"Sampling frequency: %lf\n",template->tSampling);

	fprintf(model_output,"Mass 1: %lf\n",template->mass1);
	fprintf(model_output,"Mass 2: %lf\n",template->mass2);

	for(i=0;i<model->length;i++) {
		fprintf(model_output,"%g\n",model->data[i]);
	}
	fclose(model_output);

	exit(0);


	return;

}


void PhenSpinTaylorRD_template(LALStatus *status,InspiralTemplate *template, LALMCMCParameter *parameter,LALMCMCInput *inputMCMC){

  UINT4 Nmodel,NtimeModel,idx;
  double aspin, theta, phi;
  NtimeModel =  inputMCMC->numPoints;
  Nmodel = (inputMCMC->stilde[0]->data->length-1)*2;
  REAL4Vector* TmodelPlus;//PhenSpin treats "+" and "x" polarisations differently.
  REAL4Vector* TmodelCross;//PhenSpin treats "+" and "x" polarisations differently.
	
  TmodelPlus=XLALCreateREAL4Vector(NtimeModel); /* Allocate storage for the waveform */
  TmodelCross=XLALCreateREAL4Vector(NtimeModel);/* Allocate storage for the waveform */

  template->spin1[0]=template->spin1[1]=template->spin1[2]=0;
  template->spin2[0]=template->spin2[1]=template->spin2[2]=0;

  if(XLALMCMCCheckParameter(parameter,"Spin1") && XLALMCMCCheckParameter(parameter,"Spin1theta")
     && XLALMCMCCheckParameter(parameter,"Spin1phi"))
  {
     aspin=XLALMCMCGetParameter(parameter,"Spin1");
     theta=XLALMCMCGetParameter(parameter,"Spin1theta");
     phi  =XLALMCMCGetParameter(parameter,"Spin1phi");
     template->spin1[0]=aspin*sin(theta)*cos(phi);
     template->spin1[1]=aspin*sin(theta)*sin(phi);
     template->spin1[2]=aspin*cos(theta);
  }
    if(XLALMCMCCheckParameter(parameter,"Spin2") && XLALMCMCCheckParameter(parameter,"Spin2theta")
         && XLALMCMCCheckParameter(parameter,"Spin2phi"))
  {
      aspin=XLALMCMCGetParameter(parameter,"Spin2");
      theta=XLALMCMCGetParameter(parameter,"Spin2theta");
      phi  =XLALMCMCGetParameter(parameter,"Spin2phi");
      template->spin2[0]=aspin*sin(theta)*cos(phi);
      template->spin2[1]=aspin*sin(theta)*sin(phi);
      template->spin2[2]=aspin*cos(theta);
  }
  template->inclination=XLALMCMCGetParameter(parameter,"iota");
  template->startPhase=XLALMCMCGetParameter(parameter,"phi");
  //template->OmegaS=0.0548;

  LALInspiralParameterCalc(status,template);

  double distanceMPC = template->distance;
  double distanceSI = LAL_PC_SI*1e6*distanceMPC;
  template->distance = distanceSI;
  
  LALPSpinInspiralRDTemplates(status,TmodelPlus,TmodelCross,template);/*PhenSpinTaylorRD generates plus and cross polarizations within this function*/

  template->distance = distanceMPC;

  float winNorm = sqrt(inputMCMC->window->sumofsquares/inputMCMC->window->data->length);
  float Norm = inputMCMC->deltaT/winNorm;
  for(idx=0;idx<TmodelPlus->length;idx++) TmodelPlus->data[idx]*=(REAL4)inputMCMC->window->data->data[idx]*Norm; /* window & normalise */
  for(idx=0;idx<TmodelCross->length;idx++) TmodelCross->data[idx]*=(REAL4)inputMCMC->window->data->data[idx]*Norm; /* window & normalise */  

  memset(inputMCMC->Fwfp->data,0,(inputMCMC->numPoints)*sizeof(REAL4));
  memset(inputMCMC->Fwfc->data,0,(inputMCMC->numPoints)*sizeof(REAL4));
  
  if(inputMCMC->likelihoodPlan==NULL) {
    LALCreateForwardREAL4FFTPlan(status,&inputMCMC->likelihoodPlan,Nmodel,FFTW_PATIENT);
  }
  XLALREAL4VectorFFT(inputMCMC->Fwfp,TmodelPlus,inputMCMC->likelihoodPlan);
  XLALREAL4VectorFFT(inputMCMC->Fwfc,TmodelCross,inputMCMC->likelihoodPlan); /* REAL4VectorFFT doesn't normalise like TimeFreqRealFFT, so we do this above in Norm */



   XLALDestroyREAL4Vector(TmodelPlus);
   XLALDestroyREAL4Vector(TmodelCross);
} /* End of PhenSpinTaylorRD_template */



