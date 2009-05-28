/* Nested sampling algorithm */
/* And support functions */
/* (C) John Veitch 2009 */

#include <lal/Units.h>
#include <lal/LALStdlib.h>
#include "LALInspiralMCMC.h"
#include "LALInspiralMCMCUser.h"
#include <lal/LALError.h>
#include <lal/TimeDelay.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include "nest_calc.h"
#include <float.h>

#define infosafe 1.5

double logadd(double a,double b){
if(a>b) return(a+log(1.0+exp(b-a)));
else return(b+log(1.0+exp(a-b)));
}

void NestInit2PN(LALMCMCParameter *parameter, void *iT){
REAL8 time;
SnglInspiralTable *inspiralTable = (SnglInspiralTable *)iT;
REAL4 mtot,eta,mwindow,etawindow,timewindow;
time = (REAL8)inspiralTable->end_time.gpsSeconds + (REAL8)inspiralTable->end_time.gpsNanoSeconds *1.0e-9;
parameter->param=NULL;
parameter->dimension=0;
mwindow=0.5;
REAL8 m1=inspiralTable->mass1;
REAL8 m2=inspiralTable->mass2;
mtot=inspiralTable->mass1 + inspiralTable->mass2;
eta=inspiralTable->eta;
if (eta==0.0) eta=m1*m2/((m1+m2)*(m1+m2));
double etamin = eta-0.5*etawindow;
etamin = etamin<0.01?0.01:etamin;
double etamax = eta+0.5*etawindow;
etamax = etamax>0.25?0.25:etamax;
etamin=0.05; etamax=0.25;
REAL4 localetawin=etamax-etamin;
fprintf(stderr,"clipped eta prior %lf\n",localetawin);

XLALMCMCAddParam(parameter,"mtotal",2.0+33.0*gsl_rng_uniform(RNG),2.0,35.0,0); /* Low mass */
/*XLALMCMCAddParam(parameter, "mtotal",	mtot*(1.0 + mwindow*(gsl_rng_uniform(RNG)-0.5)),mtot*(1.0-0.5*mwindow),mtot*(1.0+0.5*mwindow),0);*/
/*XLALMCMCAddParam(parameter, "eta", gsl_rng_uniform(RNG)*0.25 , 0, 0.25, 0);*/
XLALMCMCAddParam(parameter, "eta", gsl_rng_uniform(RNG)*localetawin+etamin , etamin, etamax, 0);
/*XLALMCMCAddParam(parameter, "eta",	eta*(1.0 + etawindow*(gsl_rng_uniform(RNG)-0.5)),eta*(1.0-0.5*etawindow),eta*(1.0+0.5*etawindow),0);*/
XLALMCMCAddParam(parameter, "time",		(gsl_rng_uniform(RNG)-0.5)*timewindow + time ,time-0.5*timewindow,time+0.5*timewindow,0);
XLALMCMCAddParam(parameter, "phi",		LAL_TWOPI*gsl_rng_uniform(RNG),0.0,LAL_TWOPI,1);
XLALMCMCAddParam(parameter, "distMpc", 99.0*gsl_rng_uniform(RNG)+1.0, 1.0, 100.0, 0);
XLALMCMCAddParam(parameter,"long",LAL_TWOPI*gsl_rng_uniform(RNG),0,LAL_TWOPI,1);
XLALMCMCAddParam(parameter,"lat",LAL_PI*(gsl_rng_uniform(RNG)-0.5),-LAL_PI*0.5,LAL_PI*0.5,0);
XLALMCMCAddParam(parameter,"psi",0.5*LAL_PI*gsl_rng_uniform(RNG),0,LAL_PI*0.5,0);
XLALMCMCAddParam(parameter,"iota",LAL_PI*gsl_rng_uniform(RNG),0,LAL_PI,0);

return;
}

void Inject2PN(LALMCMCParameter *parameter, LALMCMCInput *inputMCMC, double SNR){
static LALStatus status;
InspiralTemplate template;
REAL8 real,imag,chisq;
UINT4 Nmodel;
INT4 i;
double SNR1,mul_factor=1.0;

memset(&template,0,sizeof(InspiralTemplate));
template.totalMass = XLALMCMCGetParameter(parameter,"mtotal");
template.eta = XLALMCMCGetParameter(parameter,"eta");
template.massChoice = totalMassAndEta;
template.fLower = inputMCMC->fLow;
template.distance = XLALMCMCGetParameter(parameter,"distMpc");
template.order = twoPN;
template.approximant=inputMCMC->approximant;
template.tSampling = 1.0/inputMCMC->deltaT;
	template.fCutoff = 0.5/inputMCMC->deltaT -1.0;
	template.nStartPad = 0;
	template.nEndPad =0;
	template.startPhase = XLALMCMCGetParameter(parameter,"phi");
	template.startTime = XLALMCMCGetParameter(parameter,"time");
	template.startTime -= inputMCMC->stilde[0]->epoch.gpsSeconds + 1e-9*inputMCMC->stilde[0]->epoch.gpsNanoSeconds;
	template.ieta = 1;
	template.next = NULL;
	template.fine = NULL;

	LALInspiralParameterCalc(&status,&template);

	template.startTime-=template.tC;
	
	LALInspiralRestrictedAmplitude(&status,&template);
	Nmodel=inputMCMC->stilde[0]->data->length*2; /* *2 for real/imag packing format */

	LALCreateVector(&status,&model,Nmodel);
	
	LALInspiralWave(&status,model,&template); /* Create the model */

	int lowBin = (int)(inputMCMC->fLow / inputMCMC->stilde[0]->deltaF);
	/* Find the SNR of this wave */
	for(chisq=0.0,i=lowBin;i<Nmodel/2;i++){
		real=model->data[i]; imag=model->data[Nmodel-i];
		chisq+=(real*real + imag*imag)*inputMCMC->invspec[0]->data->data[i];
	}
	SNR1=sqrt(chisq);
	if(SNR>=0.0) mul_factor = SNR/SNR1; /* Multiplicative factor to achieve desired SNR */
	/* Inject the wave */
	for(chisq=0.0,i=lowBin;i<Nmodel/2;i++){
		inputMCMC->stilde[0]->data->data[i].re+=mul_factor*(REAL8) model->data[i];
		inputMCMC->stilde[0]->data->data[i].im+=mul_factor*(REAL8) model->data[Nmodel-i];
		real=model->data[i]; imag=model->data[Nmodel-i];
		chisq+=mul_factor*mul_factor*(real*real + imag*imag)*inputMCMC->invspec[0]->data->data[i];
	}
	fprintf(stderr,"Injected wave, SNR = %e\n",sqrt(chisq));
	model=NULL;
/*	free(model); */
}

REAL8 nestZ(INT4 Nruns, INT4 Nlive, LALMCMCParameter **Live, LALMCMCInput *MCMCinput)
{
	int i=0;
	int j,minpos;
	static LALStatus status;
	REAL4 accept;
	REAL8 *logZarray,*logwarray,*Harray,*oldZarray,*Wtarray;
	REAL8 logw,H=0.0,logLmin,logWt,logZ,logZnew,deltaZ;
	REAL8 MCMCfail=0;
	REAL8 logZnoise=0.0;
	REAL8 logLmax=-DBL_MAX;
	REAL4 rngseed=0;
	FILE *fpout=NULL;
	CHAR outEnd[1000];
	LALMCMCParameter *temp=(LALMCMCParameter *)malloc(sizeof(LALMCMCParameter));
	
	if(!(MCMCinput->randParams)) LALCreateRandomParams(&status,&(MCMCinput->randParams),seed);
	
	MCMCinput->Live=Live;
	MCMCinput->Nlive=(UINT4)Nlive;
	/* Initialise the RNG */
	gsl_rng_env_setup();
	RNG=gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(RNG,seed==0 ? (unsigned long int)time(NULL) : seed);
	
	/* Initialise the optimised tables declared in LALInspiralMCMCUser.h*/

	normalisations = calloc((size_t)MCMCinput->numberDataStreams,sizeof(REAL8));

	for(i=0;i<MCMCinput->numberDataStreams;i++) {
		normalisations[i]=0.0;
		if(MCMCinput->invspec[i]!=NULL){
			for(j=(int)(MCMCinput->fLow/MCMCinput->invspec[i]->deltaF);j<MCMCinput->invspec[i]->data->length-1;j++) normalisations[i]+=0.5*log(MCMCinput->invspec[i]->data->data[j]);
		}
	}
	
	
	if(MCMCinput->injectionTable!=NULL) MCMCinput->funcInit(temp,(void *)MCMCinput->injectionTable);
	else MCMCinput->funcInit(temp,(void *)MCMCinput->inspiralTable);
	
	XLALMCMCSetParameter(temp,"distMpc",DBL_MAX);
	MCMCinput->funcPrior(MCMCinput,temp);
	MCMCinput->funcLikelihood(MCMCinput,temp);
	logZnoise = temp->logLikelihood;
	fprintf(stdout,"Noise evidence: %lf\n",logZnoise);
	fprintf(stderr,"Sprinkling initial points, may take a while");
	/* Set up the parameters for the live points */
	for(i=0;i<Nlive;i++) {
		if(MCMCinput->injectionTable!=NULL) MCMCinput->funcInit(Live[i],(void *)MCMCinput->injectionTable);
		else MCMCinput->funcInit(Live[i],(void *)MCMCinput->inspiralTable);
		MCMCinput->dim=Live[i]->dimension;
		MCMCinput->funcPrior(MCMCinput,Live[i]);
		MCMCinput->funcLikelihood(MCMCinput,Live[i]);
	}
	/* Set up covariance matrix */
	cov_mat = gsl_matrix_alloc(MCMCinput->dim,MCMCinput->dim);


	calcCVM(cov_mat,Live,Nlive);
	
	for(i=0;i<Nlive;i++) {
	  accept=MCMCSampleLimitedPrior(Live[i],temp,MCMCinput,-DBL_MAX,cov_mat,MCMCinput->numberDraw);
	  if(i%50==0)fprintf(stderr,".");
	}
	if(MCMCinput->verbose) fprintf(stderr,"Set up %i live points\n",Nlive);

	/* Find max likelihood currently */
	for(i=1,logLmax=Live[0]->logLikelihood;i<Nlive;i++) logLmax=logLmax>Live[i]->logLikelihood ? logLmax : Live[i]->logLikelihood;

	/* Set up arrays for parallel runs */
	logZarray = calloc(Nruns,sizeof(REAL8));
	oldZarray = calloc(Nruns,sizeof(REAL8));
	Harray = calloc(Nruns,sizeof(REAL8));
	logwarray = calloc(Nruns,sizeof(REAL8));
	Wtarray = calloc(Nruns,sizeof(REAL8));
	if(logZarray==NULL || Harray==NULL || oldZarray==NULL || logwarray==NULL) {fprintf(stderr,"Unable to allocate RAM\n"); exit(-1);}
	
	logw=log(1.0-exp(-1.0/Nlive));
	for(i=0;i<Nruns;i++)  {logwarray[i]=logw; logZarray[i]=-DBL_MAX; oldZarray[i]=-DBL_MAX; Harray[i]=0.0;}
	i=0;
	
	/* open outfile */
	fpout=fopen(outfile,"w");
	if(fpout==NULL) fprintf(stderr,"Unable to open output file %s\n",outfile);
	
	/* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- Nested sampling loop -=-=-=-=--=-=-=-=-==-=-=-=-=-=-= */
/*	while(((REAL8)i)<=((REAL8)Nlive)*infosafe*H || i<3*Nlive) */
	deltaZ=1.0;
	/*	while((REAL8)i<=((REAL8)Nlive)*infosafe*H ? 1 : Nlive*fabs(deltaZ/logZ)>1e-6)*/
	while(((REAL8)i)<=((REAL8)Nlive) || logLmax+logw > logZ-5) /* This termination condition: when remaining prior can't
					     account for more than exp(-5) of the evidence, even
					     if entire support is at Lmax */
	{
	minpos=0;
	/* Find minimum likelihood sample to replace */
	for(j=0;j<Nlive;j++) {if(Live[j]->logLikelihood <= Live[minpos]->logLikelihood) minpos=j;}
	logLmin = Live[minpos]->logLikelihood;
	logWt = logw + Live[minpos]->logLikelihood;
	fprintSample(fpout,Live[minpos]);
	/* update evidence array */
	for(j=0;j<Nruns;j++){
		logZarray[j]=logadd(logZarray[j],Live[minpos]->logLikelihood + logwarray[j]);
		Wtarray[j]=logwarray[j]+Live[minpos]->logLikelihood;
		Harray[j]= exp(Wtarray[j]-logZarray[j])*Live[minpos]->logLikelihood
				 + exp(oldZarray[j]-logZarray[j])*(Harray[j]+oldZarray[j])-logZarray[j];
		}
	logZnew=mean(logZarray,Nruns);
	deltaZ=logZnew-logZ;
	H=mean(Harray,Nruns);
	logZ=logZnew;
	for(j=0;j<Nruns;j++) oldZarray[j]=logZarray[j];
	MCMCfail=0.0;
	/* Update covariance matrix every so often */
	if(!(i%(Nlive/4)))	calcCVM(cov_mat,Live,Nlive);
	/* generate new live point */
	do{
		while((j=gsl_rng_uniform_int(RNG,Nlive))==minpos){};
		XLALMCMCCopyPara(&(Live[minpos]),Live[j]);
		accept = MCMCSampleLimitedPrior(Live[minpos],temp,MCMCinput,logLmin,cov_mat,MCMCinput->numberDraw);
		MCMCinput->funcLikelihood(MCMCinput,Live[minpos]);
		MCMCinput->funcPrior(MCMCinput,Live[minpos]);
		MCMCfail+=1.0;
		}while((Live[minpos]->logLikelihood<=logLmin)||(accept==0.0));
	if(Live[minpos]->logLikelihood > logLmax) logLmax = Live[minpos]->logLikelihood;
	for(j=0;j<Nruns;j++) logwarray[j]+=sample_logt(Nlive);
	logw=mean(logwarray,Nruns);
	if(MCMCinput->verbose) fprintf(stderr,"%i: (%2.1lf%%) accpt: %1.3f H: %3.3lf nats (%3.3lf b) logL:%lf ->%lf logZ: %lf Zratio: %lf db\n",
						i,100.0*((REAL8)i)/(((REAL8) Nlive)*H*infosafe),accept/MCMCfail,H,H/log(2.0),logLmin,Live[minpos]->logLikelihood,logZ,10.0*log10(exp(1.0))*(logZ-logZnoise));
	if(fpout && !(i%50)) fflush(fpout);
	i++;
	}
	
	/* Sort the remaining points (not essential, just nice)*/
	for(i=0;i<Nlive-1;i++){
		minpos=i;
		logLmin=Live[i]->logLikelihood;
		for(j=i+1;j<Nlive;j++){
			if(Live[j]->logLikelihood<logLmin) {minpos=j; logLmin=Live[j]->logLikelihood;}
			}
		temp=Live[minpos]; /* Put the minimum remaining point in the current position */
		Live[minpos]=Live[i];
		Live[i]=temp;
	}

	/* final corrections */
	for(i=0;i<Nlive;i++){
		logZ=logadd(logZ,Live[i]->logLikelihood+logw);
		for(j=0;j<Nruns;j++){
			logwarray[j]+=sample_logt(Nlive);
			logZarray[j]=logadd(logZarray[j],Live[i]->logLikelihood+logwarray[j]-log((double)Nlive));
			}
		fprintSample(fpout,Live[i]);
		}
	double Npoints = MCMCinput->numberDataStreams*MCMCinput->stilde[0]->data->length-(int)(MCMCinput->fLow/MCMCinput->deltaF);
	fprintf(stdout,"MaxL = %lf\nReduced chi squared = %lf\n",Live[Nlive-1]->logLikelihood,-Live[Nlive-1]->logLikelihood/Npoints);
	logZ=mean(logZarray,Nruns);
	fprintf(stdout,"deltaLmax = %lf\n",Live[Nlive-1]->logLikelihood-logZnoise);
	double zscore =( -2.0*Live[Nlive-1]->logLikelihood - Npoints) / sqrt(2.0*Npoints);
	fprintf(stdout,"Z-score = %lf\n",zscore);
	
	close(fpout);
	sprintf(outEnd,"%s_B.txt",outfile);
	fpout=fopen(outEnd,"w");
	fprintf(fpout,"%lf %lf %lf %lf %lf\n",logZ-logZnoise,logZ,logZnoise,Live[Nlive-1]->logLikelihood-logZnoise,zscore);
	fclose(fpout);
	free(Harray); free(logwarray); free(Wtarray); free(oldZarray); free(logZarray);
	fprintf(stdout,"lodds ratio %lf\n",logZ-logZnoise);
	return logZ;

}

REAL8 mean(REAL8 *array,int N){
	REAL8 sum=0.0;
	int i;
	for(i=0;i<N;i++) sum+=array[i];
	return sum/((REAL8) N);
}

REAL4 MCMCSampleLimitedPrior(LALMCMCParameter *sample, LALMCMCParameter *temp, LALMCMCInput *MCMCInput,REAL8 minL,gsl_matrix *covM,INT4 N)
{
/* Sample from prior using MCMC to evolve the existing value of sample subject to the new likelihood being >minL*/
/* Returns acceptance ratio */
#define ROTATEFRAC 0.1
#define REFLECTFRAC 0.1
int i=0;
int a_cnt=0;
int accept=0;
int nreflect=0;
REAL8 logL,logPri,phi;
REAL8 jump_select=0;
REAL8 scale_temp;
REAL8 scalefactor_small=1E-4;
int ret=0;

MCMCInput->funcPrior(MCMCInput,sample);

XLALMCMCCopyPara(&temp,sample);

i=0;
while (i<N || (nreflect==a_cnt && nreflect>0 && nreflect%2==0)){
  i++;
  jump_select = gsl_rng_uniform(RNG);
  if(jump_select<0.1 && MCMCInput->numberDataStreams>1){
	if(MCMCInput->numberDataStreams>1) jump_select = gsl_rng_uniform(RNG);
	else jump_select=0;
	if(jump_select>0.5) {
		ret=0;
		/* Custom proposal for symmetry of the sky sphere */
		if(MCMCInput->numberDataStreams>=3 && jump_select>0.9) ret=XLALMCMCReflectDetPlane(MCMCInput,temp);
		if(ret==0) nreflect++;
		if(ret==-1 || jump_select <=0.9 || MCMCInput->numberDataStreams==2) XLALMCMCRotateSky(MCMCInput,temp);
		/* Custom proposal for mass/eta1 surface of constant 1PN term */
		/* if(gsl_rng_uniform(RNG)<0.1) XLALMCMC1PNMasseta(MCMCInput,temp); */
	}
	else  XLALMCMCDifferentialEvolution(MCMCInput,temp);
  /* Evaluate MH ratio */
  }
  else if( (jump_select=gsl_rng_uniform(RNG))<0.1) XLALMCMCDifferentialEvolution(MCMCInput,temp);
  /*else if(jump_select<0.6) XLALMCMCJumpIntrinsic(MCMCInput,temp,covM);*/
  else XLALMCMCJump(MCMCInput,temp,covM);
  
  
  MCMCInput->funcPrior(MCMCInput,temp);
  if(temp->logPrior!=-DBL_MAX && ( (temp->logPrior - sample->logPrior) > log(gsl_rng_uniform(RNG)) )) {
    /* this would be accepted based on the priors, we can now confirm that its likelihood is above the limit
       with the more expensive calculation */
    MCMCInput->funcLikelihood(MCMCInput,temp);
    if(temp->logLikelihood>minL) accept = 1;
  }
  if(accept==1) {XLALMCMCCopyPara(&sample,temp); a_cnt++; accept=0;}
  else XLALMCMCCopyPara(&temp,sample);
}
return(((REAL4) a_cnt)/((REAL4) i));
}

void calcCVM(gsl_matrix *cvm, LALMCMCParameter **samples,int N)
{ int i,j,k;
int ND=samples[0]->dimension;
REAL8 *means;
LALMCMCParam *p;
LALMCMCParam *jp;
LALMCMCParam *kp;
gsl_matrix *oldcvm;

oldcvm = gsl_matrix_alloc(ND,ND);
gsl_matrix_memcpy (oldcvm, cvm);

/* clear the matrix */
for(i=0;i<cvm->size1;i++) for(j=0;j<cvm->size2;j++) gsl_matrix_set(cvm,i,j,0.0);

/* Find the means */
if(NULL==(means = malloc((size_t)ND*sizeof(REAL8)))){fprintf(stderr,"Can't allocate RAM"); exit(-1);}
for(i=0;i<ND;i++) means[i]=0.0;
for(i=0;i<N;i++){
	p=samples[i]->param;
	for(j=0;j<ND;j++) {if(p->core->wrapping==0) {means[j]+=p->value;} p=p->next;}
	}
for(j=0;j<ND;j++) means[j]/=(REAL8)N;

/* Find the (co)-variances */
for(i=0;i<N;i++){
	kp=jp=p=samples[i]->param;
	for(j=0;j<ND;j++){
		for(k=0,kp=p;k<=j;k++){
			gsl_matrix_set(cvm,j,k,gsl_matrix_get(cvm,j,k) + (kp->value - means[k])*(jp->value - means[j]));
			kp=kp->next;
			}
		jp=jp->next;
		}
	}

/* Normalise */
for(i=0;i<ND;i++) for(j=0;j<ND;j++) gsl_matrix_set(cvm,i,j,gsl_matrix_get(cvm,i,j)/((REAL8) N));
free(means);
/* Fill in variances for angle parameters */
for(p=samples[0]->param,j=0;j<ND;j++,p=p->next) {
	if(p->core->wrapping!=0) {
		for(k=0;k<j;k++) gsl_matrix_set(cvm,j,k,0.0);
		gsl_matrix_set(cvm,j,j,ang_var(samples,p->core->name,N));
		for(k=j+1;k<ND;k++) gsl_matrix_set(cvm,k,j,0.0);
	}
}

/* the other half */
for(i=0;i<ND;i++) for(j=0;j<i;j++) gsl_matrix_set(cvm,j,i,gsl_matrix_get(cvm,i,j));

/*fprintf(stderr,"shrinkage: ");
for(i=0;i<ND;i++) fprintf(stderr,"%lf ",sqrt(gsl_matrix_get(oldcvm,i,i)/gsl_matrix_get(cvm,i,i)));
fprintf(stderr,"\n");
*/
/* Check the 2nd order first moment - indicates non-gaussian structure */
/* double *twopt=calloc(ND,sizeof(double));
double *max=calloc(ND,sizeof(double));
double *min=calloc(ND,sizeof(double));
for(k=0;k<ND;k++) {min[k]=DBL_MAX; max[k]=-DBL_MAX;}
for(i=0;i<N;i++){
	p=samples[i]->param;
	for(k=0;k<ND;k++){
		min[k]=min[k]<p->value?min[k]:p->value;
		max[k]=max[k]>p->value?max[k]:p->value;
		p=p->next;
	}
}
for(i=0;i<N;i++) for(j=i+1;j<N;j++) {
	p=samples[i]->param;
	jp=samples[j]->param;
	for(k=0;k<ND;k++) {
		twopt[k]+=fabs(p->value-jp->value);
		p=p->next; jp=jp->next;}
}


for(k=0;k<ND;k++) twopt[k]/= (1.0/(2.0*sqrt(2.0)))*((double)(N*(N-1)/2))*(max[k]-min[k]);
fprintf(stderr,"structure indicator: ");
for(k=0;k<ND;k++) fprintf(stderr,"%lf ",twopt[k]);
fprintf(stderr,"\n");
free(twopt);
free(max); free(min);
*/
/* Debug output the matrix 
for(i=0;i<ND;i++){
	for(j=0;j<ND;j++){ fprintf(stderr,"%e ",gsl_matrix_get(cvm,i,j));}
	fprintf(stderr,"\n");
	}
*/
	
return;
}

/* Calculate shortest angular distance between a1 and a2 */
REAL8 ang_dist(REAL8 a1, REAL8 a2){
	double raw = (a2>a1 ? a2-a1 : a1-a2);
	return(raw>LAL_PI ? 2.0*LAL_PI - raw : raw); 
}

/* Calculate the variance of a modulo-2pi distribution */
REAL8 ang_var(LALMCMCParameter **list,const char *pname, int N){
int i=0;
REAL8 mean=0.0;
REAL8 var=0.0;
REAL8 ms,mc;
/* Calc mean */
for(i=0,ms=0.0,mc=0.0;i<N;i++) {
	ms+=sin(XLALMCMCGetParameter(list[i],pname));
	mc+=cos(XLALMCMCGetParameter(list[i],pname));
}
ms/=N; mc/=N;
mean=atan2(ms,mc);
mean = mean<0? 2.0*LAL_PI + mean : mean;
/* calc variance */
for(i=0;i<N;i++) var+=ang_dist(XLALMCMCGetParameter(list[i],pname),mean)*ang_dist(XLALMCMCGetParameter(list[i],pname),mean);
return(var/(REAL8)N);
}

void fprintSample(FILE *fp,LALMCMCParameter *sample){
	LALMCMCParam *p=sample->param;
	if(fp==NULL) return;
	while(p!=NULL) {fprintf(fp,"%15.15lf\t",p->value); p=p->next;}
	fprintf(fp,"%lf\n",sample->logLikelihood);
	return;
}

REAL8 sample_logt(int Nlive){
	REAL8 t=0.0;
	REAL8 a=0.0;
	while((Nlive--)>1) {a=gsl_rng_uniform(RNG); t = t>a ? t : a;}
	return(log(t));
}
