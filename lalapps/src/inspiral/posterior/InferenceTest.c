/* 
 *  InferenceTest.c:  Bayesian Followup function testing site
 *
 *  Copyright (C) 2009 Ilya Mandel, Vivien Raymond, Christian Roever, Marc van der Sluys and John Veitch
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

/* example command line: */
/* ./InferenceTest --IFO [H1] --cache [/Users/john/data/triple/H1/frames.cache] --PSDstart 864162143.0 --PSDlength 1000 --srate 1024 --seglen 10 --trigtime 864162943.0  */

#include <stdio.h>
#include <lal/Date.h>
#include "LALInference.h"
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/StringInput.h>


LALVariables variables;
LALVariables variables2;
LALVariables currentParams;
LALIFOData *IfoPtr;
REAL4 number,five;
REAL8 numberR8;
INT4 numberI4;
INT8 numberI8;
COMPLEX8 numberC8;
COMPLEX16 numberC16;
REAL8 likelihood, nulllikelihood;

ProcessParamsTable *ppt, *ptr;
LALInferenceRunState *runstate=NULL;
int i, j, k;

//Test LALProposalFunction
void BasicMCMCLALProposal(LALInferenceRunState *runState, LALVariables *proposedParams);
//Test LALPriorFunction
REAL8 BasicUniformLALPrior(LALInferenceRunState *runState, LALVariables *params);
//Test LALEvolveOneStepFunction
void BasicMCMCOneStep(LALInferenceRunState *runState);
//Test LALAlgorithm
void MCMCAlgorithm (struct tagLALInferenceRunState *runState);

// gsl_rng * InitializeRandomSeed(void);
// unsigned long int random_seed();



int main(int argc, char *argv[]){
  fprintf(stdout," ========== InferenceTest.c ==========\n");

  /* test "LALVariables" stuff: */
  number = 10.0;
  LALStatus status;	
  five=5.0;
  variables.head=NULL;
  variables.dimension=0;
	
  memset(&status,0,sizeof(status));
  addVariable(&variables, "number", &number, REAL4_t);
  numberR8 = 7.0;
  addVariable(&variables, "seven", &numberR8, REAL8_t);
  numberR8 = LAL_PI;
  addVariable(&variables, "pi", &numberR8, REAL8_t);
  numberI4 = 123;
  addVariable(&variables, "small", &numberI4, INT4_t);
  numberI8 = 256*256*256*64;
  addVariable(&variables, "large", &numberI8, INT8_t);
  numberC8.re = 2.0;  numberC8.im = 3.0;
  addVariable(&variables, "complex1", &numberC8, COMPLEX8_t);
  numberC16.re = 1.23;  numberC16.im = -3.45;
  addVariable(&variables, "complex2", &numberC16, COMPLEX16_t);

  number=*(REAL4 *)getVariable(&variables,"number");
  fprintf(stdout,"Got %lf\n",number);
  setVariable(&variables,"number",&five);
  number=*(REAL4 *)getVariable(&variables,"number");
  fprintf(stdout,"Got %lf\n",number);
  fprintf(stdout,"Checkvariable?: %i\n",checkVariable(&variables,"number"));
  printVariables(&variables);
  copyVariables(&variables, &variables2);
  printVariables(&variables2);
  fprintf(stdout,"compareVariables?: %i\n",
          compareVariables(&variables,&variables2));
  numberC16.im = 4.56;
  setVariable(&variables2,"complex2",&numberC16);
  fprintf(stdout,"compareVariables?: %i\n",
          compareVariables(&variables,&variables2));
  numberC16.im = -3.45;
  setVariable(&variables2,"complex2",&numberC16);
  fprintf(stdout,"compareVariables?: %i\n",
          compareVariables(&variables,&variables2));

  removeVariable(&variables,"number");
  fprintf(stdout,"Removed, Checkvariable?: %i\n",checkVariable(&variables,"number"));
  
  fprintf(stdout,"compareVariables?: %i\n",
          compareVariables(&variables,&variables2));
  destroyVariables(&variables);
  destroyVariables(&variables2);
  printVariables(&variables2);

  fprintf(stdout," ----------\n");
  
  /* test "parseCommandLine()" function: */
  ppt = (ProcessParamsTable*) parseCommandLine(argc,argv);
  printf("parsed command line arguments:\n");
  ptr = ppt;
  i=1;
  while (ptr != NULL){
    printf(" (%d)  %s  %s  %s  \"%s\"\n", i, ptr->program, ptr->param, ptr->type, ptr->value);
    ptr = ptr->next;
    ++i;
  }

  fprintf(stdout," ----------\n");

  /* Test the data initialisation &c. */
  runstate = initialize(ppt);
  
  if(runstate->data) {
    fprintf(stdout," data found --> trying some template computations etc.\n");
    
    /* print some information on individual "runstate->data" elements: */
    IfoPtr = runstate->data;  i = 1;
    while (IfoPtr != NULL) {
      if (IfoPtr->timeData)
        fprintf(stdout, " [%d] timeData (\"%s\"): length=%d, deltaT=%f, epoch=%.3f\n", 
                i, IfoPtr->timeData->name, IfoPtr->timeData->data->length, IfoPtr->timeData->deltaT, 
                XLALGPSGetREAL8(&IfoPtr->timeData->epoch));
      if (IfoPtr->freqData)
        fprintf(stdout, "     freqData (\"%s\"): length=%d, deltaF=%f\n", 
                IfoPtr->freqData->name, IfoPtr->freqData->data->length, IfoPtr->freqData->deltaF);
      fprintf(stdout, "     fLow=%.1f Hz,  fHigh=%.1f Hz  (%d freq bins w/in range)\n", 
              IfoPtr->fLow, IfoPtr->fHigh, 
              ((int) (floor(IfoPtr->fHigh / IfoPtr->freqData->deltaF) - ceil(IfoPtr->fLow / IfoPtr->freqData->deltaF)))+1);
      fprintf(stdout, "     detector location: (%.1f, %.1f, %.1f)\n",
              IfoPtr->detector->location[0], IfoPtr->detector->location[1], IfoPtr->detector->location[2]);
      fprintf(stdout, "     detector response matrix:\n");
      for (j=0; j<3; ++j){
        fprintf(stdout, "     ");
        for (k=0; k<3; ++k)
          fprintf(stdout, "%f  ", IfoPtr->detector->response[j][k]);
        fprintf(stdout, "\n");
      }
      IfoPtr = IfoPtr->next;
    }

    REAL4 m1 = 10.0;
    addVariable(runstate->data->modelParams,"m1",&m1,REAL4_t);
    REAL4 m2 = 1.4;
    addVariable(runstate->data->modelParams,"m2",&m2,REAL4_t);
    REAL4 inc = 0.0;
    addVariable(runstate->data->modelParams,"inc",&inc,REAL4_t);
    REAL4 phii = 0.0;
    addVariable(runstate->data->modelParams,"phii",&phii,REAL4_t);
	ProcessParamsTable *procparam=getProcParamVal(ppt,"--trigtime");
	LIGOTimeGPS trigger_time;
	char * chartmp;
	LALStringToGPS(&status,&trigger_time,procparam->value,&chartmp);
	REAL8 tc = XLALGPSGetREAL8(&trigger_time);
	addVariable(runstate->data->modelParams,"time",&tc,REAL8_t);
	
    LALTemplateGeneratePPN(runstate->data);
	  executeFT(runstate->data);
	  
	  FILE *testout=fopen("test_FD.txt","w");
	  for (i=0;i<runstate->data->freqModelhPlus->data->length;i++){
		  fprintf(testout,"%g %g %g %g %g\n",i*runstate->data->freqModelhPlus->deltaF,
				  runstate->data->freqModelhPlus->data->data[i].re,
				  runstate->data->freqModelhPlus->data->data[i].im,
				  runstate->data->freqModelhCross->data->data[i].re,
				  runstate->data->freqModelhCross->data->data[i].im);
	  }
	  fclose(testout);
	  testout=fopen("test_TD.txt","w");
	  for (i=0;i<runstate->data->timeModelhPlus->data->length;i++){
		  fprintf(testout,"%10.10lf %g %g\n",runstate->data->timeData->epoch.gpsSeconds+(1e-9*runstate->data->timeData->epoch.gpsNanoSeconds)+i*runstate->data->timeModelhPlus->deltaT,
				  runstate->data->timeModelhPlus->data->data[i],
				  runstate->data->timeModelhCross->data->data[i]);
	  }
	  fclose(testout);
	  testout=fopen("PSD.txt","w");
	  for (i=0;i<runstate->data->oneSidedNoisePowerSpectrum->data->length;i++){
		  fprintf(testout,"%g %g\n",i*runstate->data->oneSidedNoisePowerSpectrum->deltaF,
				  runstate->data->oneSidedNoisePowerSpectrum->data->data[i]);
	  }
	  fclose(testout);
	  testout=fopen("noise_TD.txt","w");
	  for (i=0;i<runstate->data->timeData->data->length;i++){
		  fprintf(testout,"%10.10lf %g\n",runstate->data->timeData->epoch.gpsSeconds+(1e-9*runstate->data->timeData->epoch.gpsNanoSeconds)+i*runstate->data->timeData->deltaT,
				  runstate->data->timeData->data->data[i]);
	  }
	  fclose(testout);
	  testout=fopen("noise_FD.txt","w");
	  for (i=0;i<runstate->data->freqData->data->length;i++){
	          //fprintf(testout,"%g %g %g %g %g\n",i*runstate->data->freqData->deltaF,
		  fprintf(testout,"%g %g %g\n",i*runstate->data->freqData->deltaF,
				  runstate->data->freqData->data->data[i].re,
			          //runstate->data->freqData->data->data[i].im,
				  //runstate->data->freqData->data->data[i].re,
				  runstate->data->freqData->data->data[i].im);
	  }
	  
	  
    fprintf(stdout," ----------\n");


    //  templateStatPhase() test: 
    fprintf(stdout, " trying out 'templateStatPhase()'...\n");
    REAL8 mc   = 4.0;
    REAL8 eta  = 0.24;
    REAL8 iota = 0.0;
    REAL8 phi  = 2.0;
    REAL8 tcoal   = XLALGPSGetREAL8(&(runstate->data->timeData->epoch)) + (((double)runstate->data->timeData->data->length) * runstate->data->timeData->deltaT) - 1.0;
    printf("TCOAL: %f\n",tcoal);
    destroyVariables(runstate->data->modelParams);
    addVariable(runstate->data->modelParams, "chirpmass",   &mc,    REAL8_t);
    addVariable(runstate->data->modelParams, "massratio",   &eta,   REAL8_t);
    addVariable(runstate->data->modelParams, "inclination", &iota,  REAL8_t);
    addVariable(runstate->data->modelParams, "phase",       &phi,   REAL8_t);
    addVariable(runstate->data->modelParams, "time",        &tcoal, REAL8_t);
    printVariables(runstate->data->modelParams);
    templateStatPhase(runstate->data);
    fprintf(stdout, " ...done.\n");

	  
	  // Parameters for which I am going to compute the likelihood
	  
	  REAL8 ra_current        = 0.0;	/* radian      */
	  REAL8 dec_current       = 0.0;	/* radian      */
	  REAL8 psi_current       = 0.8;	/* radian      */
	  REAL8 distMpc_current   = 10.0;	/* Mpc         */
	  
    addVariable(&currentParams, "chirpmass",       &mc,              REAL8_t);
    addVariable(&currentParams, "massratio",       &eta,             REAL8_t);
    addVariable(&currentParams, "inclination",     &iota,            REAL8_t);
    addVariable(&currentParams, "phase",           &phi,             REAL8_t);
    addVariable(&currentParams, "time",            &tc   ,           REAL8_t); 
    addVariable(&currentParams, "rightascension",  &ra_current,      REAL8_t);
    addVariable(&currentParams, "declination",     &dec_current,     REAL8_t);
    addVariable(&currentParams, "polarisation",    &psi_current,     REAL8_t);
    addVariable(&currentParams, "distance",        &distMpc_current, REAL8_t);
    fprintf(stdout, " trying 'templateLAL' likelihood...\n");
    numberI4 = TaylorT1;
    addVariable(&currentParams, "LAL_APPROXIMANT", &numberI4,        INT4_t);
    numberI4 = LAL_PNORDER_TWO;
    addVariable(&currentParams, "LAL_PNORDER",     &numberI4,        INT4_t);
    likelihood = FreqDomainLogLikelihood(&currentParams, runstate->data, templateLAL);
    fprintf(stdout, " ...done.\n");
    fprintf(stdout," templateLAL log-likelihood %f\n", likelihood);  
    fprintf(stdout," ----------\n");
		
	//Single IFO likelihood test
	fprintf(stdout, "Single IFO likelihood test\n");
	COMPLEX16Vector *freqModel1=XLALCreateCOMPLEX16Vector(runstate->data->freqData->data->length);
	COMPLEX16Vector *freqModel2=XLALCreateCOMPLEX16Vector(runstate->data->freqData->data->length);
	numberI4 = LAL_PNORDER_TWO;
    setVariable(&currentParams, "LAL_PNORDER",     &numberI4);	
	numberI4 = TaylorT1;
    setVariable(&currentParams, "LAL_APPROXIMANT", &numberI4);														  																  
	ComputeFreqDomainResponse(&currentParams, runstate->data, templateLAL, freqModel1);
	freqModel2=runstate->data->freqData->data;
	//ComputeFreqDomainResponse(&currentParams, runstate->data, templateLAL, freqModel2);
	FILE * freqModelFile=fopen("freqModelFile.dat", "w");
	for(i=0; i<runstate->data->freqData->data->length; i++){
		fprintf(freqModelFile, "%g\t %g\t %g\t %g\t %g\t %g\n", 
		((double)i)*1.0/ (((double)runstate->data->timeData->data->length) * runstate->data->timeData->deltaT),
		freqModel1->data[i].re, freqModel1->data[i].im, freqModel2->data[i].re, freqModel2->data[i].im,
		runstate->data->oneSidedNoisePowerSpectrum->data->data[i]);
	}
	fprintf(stdout, "overlap=%g\n", 
		ComputeFrequencyDomainOverlap(runstate->data, freqModel1, freqModel2));
	fprintf(stdout, "<d|d>=%g, <d|h>=%g, <h|h>=%g, <d|h>-1/2<h|h>=%g\n", 
		ComputeFrequencyDomainOverlap(runstate->data, freqModel2, freqModel2),
		ComputeFrequencyDomainOverlap(runstate->data, freqModel1, freqModel2),
		ComputeFrequencyDomainOverlap(runstate->data, freqModel1, freqModel1),
		ComputeFrequencyDomainOverlap(runstate->data, freqModel2, freqModel1)-0.5*ComputeFrequencyDomainOverlap(runstate->data, freqModel1, freqModel1)
		);				
	fprintf(stdout, "likelihood %g\n",
		FreqDomainLogLikelihood(&currentParams, runstate->data, templateLAL));
	fprintf(stdout, "undecomposed likelihood %g \n", 
		UndecomposedFreqDomainLogLikelihood(&currentParams, runstate->data, templateLAL));
	fprintf(stdout, "null likelihood %g decomposed null likelihood %g\n",
		FreqDomainNullLogLikelihood(runstate->data),
		NullLogLikelihood(runstate->data));
    XLALDestroyCOMPLEX16Vector(freqModel1);
	XLALDestroyCOMPLEX16Vector(freqModel2);

    /* NOTE: try out the "forceTimeLocation" flag within the "templateLAL()" function */
    /*       for aligning (time domain) templates.                                    */
    fprintf(stdout," generating templates & writing to files...:\n");
    dumptemplateFreqDomain(&currentParams, runstate->data, templateStatPhase, "test_FTemplate25SP.csv");
    dumptemplateTimeDomain(&currentParams, runstate->data, templateStatPhase, "test_TTemplate25SP.csv");
    dumptemplateFreqDomain(&currentParams, runstate->data, template3525TD, "test_FTemplate3525TD.csv");
    dumptemplateTimeDomain(&currentParams, runstate->data, template3525TD, "test_TTemplate3525TD.csv");

    fprintf(stdout," ----------\n");

    /* These are the LAL templates that (...seem to...) work right now: */
    /* TaylorT1, TaylorT2, TaylorT3, TaylorF2, IMRPhenomA, PadeT1, EOB  */
    numberI4 = LAL_PNORDER_TWO;
    setVariable(&currentParams, "LAL_PNORDER",     &numberI4);
    numberI4 = TaylorF2;
    setVariable(&currentParams, "LAL_APPROXIMANT", &numberI4);
    dumptemplateTimeDomain(&currentParams, runstate->data, templateLAL, "test_TTemplateLAL-TF2.csv");
    numberI4 = TaylorT1;
    setVariable(&currentParams, "LAL_APPROXIMANT", &numberI4);
    dumptemplateTimeDomain(&currentParams, runstate->data, templateLAL, "test_TTemplateLAL-TT1.csv");
    numberI4 = TaylorT2;
    setVariable(&currentParams, "LAL_APPROXIMANT", &numberI4);
    dumptemplateTimeDomain(&currentParams, runstate->data, templateLAL, "test_TTemplateLAL-TT2.csv");
    numberI4 = TaylorT3;
    setVariable(&currentParams, "LAL_APPROXIMANT", &numberI4);
    dumptemplateTimeDomain(&currentParams, runstate->data, templateLAL, "test_TTemplateLAL-TT3.csv");

    numberI4 = IMRPhenomA;
    setVariable(&currentParams, "LAL_APPROXIMANT", &numberI4);
    dumptemplateTimeDomain(&currentParams, runstate->data, templateLAL, "test_TTemplateLAL-Phenom.csv");
	dumptemplateFreqDomain(&currentParams, runstate->data, templateLAL, "test_FTemplateLAL-Phenom.csv");

    numberI4 = PadeT1;
    setVariable(&currentParams, "LAL_APPROXIMANT", &numberI4);
    dumptemplateTimeDomain(&currentParams, runstate->data, templateLAL, "test_TTemplateLAL-PadeT1.csv");

    numberI4 = EOB;
    setVariable(&currentParams, "LAL_APPROXIMANT", &numberI4);
    numberI4 = LAL_PNORDER_PSEUDO_FOUR;
    setVariable(&currentParams, "LAL_PNORDER", &numberI4);
    dumptemplateTimeDomain(&currentParams, runstate->data, templateLAL, "test_TTemplateLAL-EOB.csv");

    numberI4 = BCV;
    setVariable(&currentParams, "LAL_APPROXIMANT", &numberI4);
    numberI4 = LAL_PNORDER_TWO;
    setVariable(&currentParams, "LAL_PNORDER",     &numberI4);
    dumptemplateTimeDomain(&currentParams, runstate->data, templateLAL, "test_TTemplateLAL-BCV.csv");

    numberI4 = EOBNR;
    setVariable(&currentParams, "LAL_APPROXIMANT", &numberI4);
    numberI4 = LAL_PNORDER_PSEUDO_FOUR;
    setVariable(&currentParams, "LAL_PNORDER", &numberI4);
    dumptemplateTimeDomain(&currentParams, runstate->data, templateLAL, "test_TTemplateLAL-EOBNR.csv");

    fprintf(stdout," ----------\n");

    numberR8 = 440;
    addVariable(&currentParams, "frequency", &numberR8, REAL8_t);
    numberR8 = 1e-19;
    addVariable(&currentParams, "amplitude", &numberR8, REAL8_t);
    dumptemplateTimeDomain(&currentParams, runstate->data, templateSinc, "test_TTemplateSinc.csv");

    numberR8 = 0.01;
    addVariable(&currentParams, "sigma", &numberR8, REAL8_t);
    dumptemplateTimeDomain(&currentParams, runstate->data, templateSineGaussian, "test_TTemplateSineGauss.csv");

    numberR8 = 0.01;
    addVariable(&currentParams, "tau", &numberR8, REAL8_t);
    dumptemplateTimeDomain(&currentParams, runstate->data, templateDampedSinusoid, "test_TTemplateDampedSinus.csv");

    destroyVariables(&currentParams);
    fprintf(stdout," ----------\n");
  }

  printf(" ========== main(): finished. ==========\n");
  return 0;
}



//Test LALProposalFunction
void BasicMCMCLALProposal(LALInferenceRunState *runState, LALVariables *proposedParams)
/****************************************/
/* Assumes the following parameters		*/
/* exist (e.g., for TaylorT1):			*/
/* chirpmass, massratio, inclination,	*/
/* phase, time, rightascension,			*/
/* desclination, polarisation, distance.*/
/* Simply picks a new value based on	*/
/* fixed Gaussian if within prior;		*/
/* need smarter wall bounces in future.	*/
/****************************************/
{
	REAL8 mc, eta, iota, phi, tc, ra, dec, psi, dist;
	REAL8 mc_proposed, eta_proposed, iota_proposed, phi_proposed, tc_proposed, ra_proposed, dec_proposed, psi_proposed, dist_proposed;
	gsl_rng * GSLrandom=runState->GSLrandom;
	LALVariables * currentParams = runState->currentParams;	
	//int withinPrior=0;

	destroyVariables(proposedParams);

	mc		= *(REAL8*) getVariable(currentParams, "chirpmass");		/* solar masses*/
	eta		= *(REAL8*) getVariable(currentParams, "massratio");		/* dim-less    */
	iota	= *(REAL8*) getVariable(currentParams, "inclination");		/* radian      */
	tc		= *(REAL8*) getVariable(currentParams, "time");				/* GPS seconds */
	phi		= *(REAL8*) getVariable(currentParams, "phase");			/* radian      */
	ra      = *(REAL8*) getVariable(currentParams, "rightascension");	/* radian      */
	dec     = *(REAL8*) getVariable(currentParams, "declination");		/* radian      */
	psi     = *(REAL8*) getVariable(currentParams, "polarisation");		/* radian      */
	dist	= *(REAL8*) getVariable(currentParams, "distance");			/* Mpc         */

	addVariable(proposedParams, "chirpmass",       &mc,			REAL8_t);
	addVariable(proposedParams, "massratio",       &eta,		REAL8_t);
	addVariable(proposedParams, "inclination",     &iota,		REAL8_t);
	addVariable(proposedParams, "phase",           &phi,		REAL8_t);
	addVariable(proposedParams, "time",            &tc,			REAL8_t); 
	addVariable(proposedParams, "rightascension",  &ra,			REAL8_t);
	addVariable(proposedParams, "declination",     &dec,		REAL8_t);
	addVariable(proposedParams, "polarisation",    &psi,		REAL8_t);
	addVariable(proposedParams, "distance",        &dist,		REAL8_t);  
  
        /* the "withinPrior" restriction would screw up the sampling (not symmetric any more)!! */
	//while(!withinPrior)
	//{
		mc_proposed=mc*(1.0+gsl_ran_ugaussian(GSLrandom)*0.01);	/*mc changed by 1% */
		eta_proposed=eta+gsl_ran_ugaussian(GSLrandom)*0.01; /*eta changed by 0.01*/
		//TODO: if(eta_proposed>0.25) eta_proposed=0.25-(eta_proposed-0.25); etc.
		iota_proposed=iota+gsl_ran_ugaussian(GSLrandom)*0.1;
		tc_proposed=tc+gsl_ran_ugaussian(GSLrandom)*0.005; /*time changed by 5 ms*/
		phi_proposed=phi+gsl_ran_ugaussian(GSLrandom)*0.5;
		ra_proposed=ra+gsl_ran_ugaussian(GSLrandom)*0.05;
		dec_proposed=dec+gsl_ran_ugaussian(GSLrandom)*0.05;
		psi_proposed=psi+gsl_ran_ugaussian(GSLrandom)*0.1;
		dist_proposed=dist*(1.0+gsl_ran_ugaussian(GSLrandom)*0.3);
	
		setVariable(proposedParams, "chirpmass",       &mc_proposed);
		setVariable(proposedParams, "massratio",       &eta_proposed);
		setVariable(proposedParams, "inclination",     &iota_proposed);
		setVariable(proposedParams, "phase",           &phi_proposed);
		setVariable(proposedParams, "time",            &tc_proposed); 
		setVariable(proposedParams, "rightascension",  &ra_proposed);
		setVariable(proposedParams, "declination",     &dec_proposed);
		setVariable(proposedParams, "polarisation",    &psi_proposed);
		setVariable(proposedParams, "distance",        &dist_proposed);
	
	//	withinPrior= runState->prior(runState, proposedParams) > 0;
	//}
}

//Test LALPriorFunction
REAL8 BasicUniformLALPrior(LALInferenceRunState *runState, LALVariables *params)
/****************************************/
/* Assumes the following parameters		*/
/* exist (e.g., for TaylorT1):			*/
/* chirpmass, massratio, inclination,	*/
/* phase, time, rightascension,			*/
/* desclination, polarisation, distance.*/
/* Prior is flat if within range		*/
/****************************************/
{
	REAL8 mc, eta, iota, phi, tc, ra, dec, psi, dist;	
	
	mc		= *(REAL8*) getVariable(params, "chirpmass");		/* solar masses*/
	eta		= *(REAL8*) getVariable(params, "massratio");		/* dim-less    */
	iota	= *(REAL8*) getVariable(params, "inclination");		/* radian      */
	tc		= *(REAL8*) getVariable(params, "time");			/* GPS seconds */
	phi		= *(REAL8*) getVariable(params, "phase");			/* radian      */
	ra      = *(REAL8*) getVariable(params, "rightascension");	/* radian      */
	dec     = *(REAL8*) getVariable(params, "declination");		/* radian      */
	psi     = *(REAL8*) getVariable(params, "polarisation");	/* radian      */
	dist	= *(REAL8*) getVariable(params, "distance");		/* Mpc         */
	
	if(eta>0 && eta<=0.25 && iota>=0 && iota <=2.0*pi && phi>=0 && phi <=2.0*pi && ra>=0 && ra<=24.0 
				&& dec>=-pi/2.0 && dec<=pi/2.0 && psi>=0 && psi<=2.0*pi)	
		return 1.0;
	//TODO: should be properly normalized; pass in range via priorArgs?	
	return 0;
}

//Test LALEvolveOneStepFunction
void BasicMCMCOneStep(LALInferenceRunState *runState)
{
	REAL8 priorCurrent, priorProposed=0;
	LALVariables proposedParams;
	REAL8 likelihoodCurrent, likelihoodProposed;
	REAL8 acceptanceProbability;
	
	priorCurrent=runState->prior(runState, runState->currentParams);

        /* again, the within-prior restriction would screw up the sampler!! */
	//while(priorProposed<=0){	
		runState->proposal(runState, &proposedParams);
		priorProposed=runState->prior(runState, &proposedParams);
	//}

        /* you can still save computation time by doing the following: */
	if (priorProposed > 0) 
          likelihoodProposed = runState->likelihood(&proposedParams, runState->data, runState->template);
        else
          likelihoodProposed = 0;

	likelihoodCurrent=runState->likelihood(runState->currentParams, runState->data, runState->template);
        /* (the above number should have been stored somwehere from the previous iteration) */
	
	acceptanceProbability=likelihoodCurrent*priorCurrent/(likelihoodProposed*priorProposed);
	
	if(acceptanceProbability > gsl_rng_uniform(runState->GSLrandom))	//accept
		copyVariables(&proposedParams, runState->currentParams);

	destroyVariables(&proposedParams);	
}

//Test LALAlgorithm
void MCMCAlgorithm (struct tagLALInferenceRunState *runState)
{
	int i;
	//runState->GSLrandom=InitializeRandomSeed();
	for(i=0; i<100; i++)
	{
		runState->evolve(runState);
		printf("MCMC step %d.  Current parameters are:\n", i);
		printVariables(runState->currentParams);
	}
}


////Random variant algorithm initialization
////see GSL documentation
////http://www.gnu.org/software/gsl/manual/
//gsl_rng * InitializeRandomSeed(void)
//{
//        const gsl_rng_type * T;
//        gsl_rng * rng;
//        gsl_rng_env_setup();
//        T = gsl_rng_default;
//        rng = gsl_rng_alloc (T);
//        //long seed = time(NULL)*getpid();
//        unsigned long seed=random_seed();
//        gsl_rng_set(rng, seed);
//        return rng;
//}
//
//
////http://www.cygwin.com/ml/gsl-discuss/2004-q1/msg00072.html
//unsigned long int random_seed(){
//        unsigned int seed;
//        struct timeval tv;
//        FILE *devrandom;
//
//        if ((devrandom = fopen("/dev/random","r")) == NULL) {
//                gettimeofday(&tv,0);
//                seed = tv.tv_sec + tv.tv_usec;
//        } else {
//                fread(&seed,sizeof(seed),1,devrandom);
//                fclose(devrandom);
//        }
//
//        return seed;
//}
//
// (I moved the above part into the "initialize()" function.  Hope that's OK.  Christian.)
