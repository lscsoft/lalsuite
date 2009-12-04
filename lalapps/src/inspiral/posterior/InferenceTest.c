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
/* 
./InferenceTest --IFO [H1] --cache [/Users/john/data/triple/H1/frames.cache] --PSDstart 864162143.0 --PSDlength 1000 --srate 1024 --seglen 10 --trigtime 864162943.0
*/

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
void ASinOmegaTProposal(LALInferenceRunState *runState, LALVariables *proposedParams);
//Test LALPriorFunction
REAL8 BasicUniformLALPrior(LALInferenceRunState *runState, LALVariables *params);
REAL8 ASinOmegaTPrior(LALInferenceRunState *runState, LALVariables *params);
//Test LALEvolveOneStepFunction
void BasicMCMCOneStep(LALInferenceRunState *runState);
//Test LALAlgorithm
void MCMCAlgorithm (struct tagLALInferenceRunState *runState);
void NelderMeadEval(struct tagLALInferenceRunState *runState,
                    char **names, REAL8 *values, int dim,
                    REAL8 *logprior, REAL8 *loglikelihood);
void NelderMeadAlgorithm(struct tagLALInferenceRunState *runState, LALVariables *subset);

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
		  fprintf(testout,"%10.10lf %g %g\n",runstate->data->timeData->epoch.gpsSeconds
					+(1e-9*runstate->data->timeData->epoch.gpsNanoSeconds)+i*runstate->data->timeModelhPlus->deltaT,
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
		  fprintf(testout,"%10.10lf %g\n",runstate->data->timeData->epoch.gpsSeconds
					+(1e-9*runstate->data->timeData->epoch.gpsNanoSeconds)+i*runstate->data->timeData->deltaT,
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
    REAL8 iota = 0.4;
    REAL8 phi  = 2.0;
    REAL8 tcoal   = XLALGPSGetREAL8(&(runstate->data->timeData->epoch)) + 
		(((double)runstate->data->timeData->data->length) * runstate->data->timeData->deltaT) - 1.0;
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
		ComputeFrequencyDomainOverlap(runstate->data, freqModel2, freqModel1)
			-0.5*ComputeFrequencyDomainOverlap(runstate->data, freqModel1, freqModel1)
		);				
	fprintf(stdout, "likelihood %g\n",
		FreqDomainLogLikelihood(&currentParams, runstate->data, templateLAL));
	fprintf(stdout, "undecomposed likelihood %g \n", 
		UndecomposedFreqDomainLogLikelihood(&currentParams, runstate->data, templateLAL));
	fprintf(stdout, "null likelihood %g decomposed null likelihood %g\n",
		FreqDomainNullLogLikelihood(runstate->data),
		NullLogLikelihood(runstate->data));
    XLALDestroyCOMPLEX16Vector(freqModel1);
//	XLALDestroyCOMPLEX16Vector(freqModel2);
	
	//MCMC basic Sampler test
	fprintf(stdout, "Try MCMC basic Sampler test\n");
	runstate->algorithm=MCMCAlgorithm;
	runstate->evolve=BasicMCMCOneStep;
	runstate->prior=BasicUniformLALPrior;
	runstate->proposal=BasicMCMCLALProposal;
        runstate->proposalArgs = malloc(sizeof(LALVariables));
        runstate->proposalArgs->head=NULL;
        runstate->proposalArgs->dimension=0;
	runstate->likelihood=FreqDomainLogLikelihood;
	//runstate->template=templateLAL;
	runstate->template=templateStatPhase;
	runstate->currentParams=&currentParams;
	MCMCAlgorithm(runstate);
	fprintf(stdout, "End of MCMC basic Sampler test\n");
        NelderMeadAlgorithm(runstate, NULL);

    /* NOTE: try out the "forceTimeLocation" flag within the "templateLAL()" function */
    /*       for aligning (time domain) templates.                                    */
    fprintf(stdout," generating templates & writing to files...:\n");
    dumptemplateFreqDomain(&currentParams, runstate->data, templateStatPhase, "test_FTemplate25SP.csv");
    dumptemplateTimeDomain(&currentParams, runstate->data, templateStatPhase, "test_TTemplate25SP.csv");
    dumptemplateFreqDomain(&currentParams, runstate->data, template3525TD, "test_FTemplate3525TD.csv");
    dumptemplateTimeDomain(&currentParams, runstate->data, template3525TD, "test_TTemplate3525TD.csv");

    fprintf(stdout," ----------\n");
	 
	 double mass1=10.;
	 double mass2=1.4;
    addVariable(&currentParams, "m1",       &mass1,              REAL8_t);
	addVariable(&currentParams, "m2",       &mass2,              REAL8_t);
	  double spin1x = 0.5;
	  double spin1y = 0.1;
	  double spin1z = 0.0;
	  double spin2x = 0.2;
	  double spin2y = 0.0;
	  double spin2z = 0.3;
	  addVariable(&currentParams, "spin1x",       &spin1x,              REAL8_t);	  
	  addVariable(&currentParams, "spin1y",       &spin1y,              REAL8_t);
	  addVariable(&currentParams, "spin1z",       &spin1z,              REAL8_t);
	  addVariable(&currentParams, "spin2x",       &spin2x,              REAL8_t);	  
	  addVariable(&currentParams, "spin2y",       &spin2y,              REAL8_t);	  
	  addVariable(&currentParams, "spin2z",       &spin2z,              REAL8_t);
	  double shift0 = 0.3;
	  addVariable(&currentParams, "shift0",       &shift0,              REAL8_t);
	  double coa_phase = 0.1;
	  addVariable(&currentParams, "coa_phase",    &coa_phase,           REAL8_t);	  
	  double PNorder = 3.5;
	  addVariable(&currentParams, "PNorder",      &PNorder,             REAL8_t);	  
	  dumptemplateTimeDomain(&currentParams, runstate->data, templateLALSTPN, "test_TTemplateLALSTPN.csv");

	  
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
/* fixed Gaussian;						*/
/* need smarter wall bounces in future.	*/
/****************************************/
{
  REAL8 mc, eta, iota, phi, tc, ra, dec, psi, dist;
  REAL8 mc_proposed, eta_proposed, iota_proposed, phi_proposed, tc_proposed, 
        ra_proposed, dec_proposed, psi_proposed, dist_proposed;
  REAL8 logProposalRatio = 0.0;  // = log(P(backward)/P(forward))
  gsl_rng * GSLrandom=runState->GSLrandom;
  LALVariables * currentParams = runState->currentParams;	

  mc   = *(REAL8*) getVariable(currentParams, "chirpmass");		/* solar masses*/
  eta  = *(REAL8*) getVariable(currentParams, "massratio");		/* dim-less    */
  iota = *(REAL8*) getVariable(currentParams, "inclination");		/* radian      */
  tc   = *(REAL8*) getVariable(currentParams, "time");				/* GPS seconds */
  phi  = *(REAL8*) getVariable(currentParams, "phase");			/* radian      */
  ra   = *(REAL8*) getVariable(currentParams, "rightascension");	/* radian      */
  dec  = *(REAL8*) getVariable(currentParams, "declination");		/* radian      */
  psi  = *(REAL8*) getVariable(currentParams, "polarisation");		/* radian      */
  dist = *(REAL8*) getVariable(currentParams, "distance");			/* Mpc         */

  //mc_proposed   = mc*(1.0+gsl_ran_ugaussian(GSLrandom)*0.01);	/*mc changed by 1% */
  // (above proposal is not symmetric!)
  //mc_proposed   = mc   + gsl_ran_ugaussian(GSLrandom)*0.0001;	/*mc changed by 0.0001 */
  mc_proposed   = mc * exp(gsl_ran_ugaussian(GSLrandom)*0.001);          /* mc changed by ~0.1% */
  logProposalRatio *= mc_proposed / mc;   // (proposal ratio for above "scaled log-normal" proposal)
  eta_proposed  = eta  + gsl_ran_ugaussian(GSLrandom)*0.01; /*eta changed by 0.01*/
  //TODO: if(eta_proposed>0.25) eta_proposed=0.25-(eta_proposed-0.25); etc.
  iota_proposed = iota + gsl_ran_ugaussian(GSLrandom)*0.1;
  tc_proposed   = tc   + gsl_ran_ugaussian(GSLrandom)*0.005; /*time changed by 5 ms*/
  phi_proposed  = phi  + gsl_ran_ugaussian(GSLrandom)*0.5;
  ra_proposed   = ra   + gsl_ran_ugaussian(GSLrandom)*0.05;
  dec_proposed  = dec  + gsl_ran_ugaussian(GSLrandom)*0.05;
  psi_proposed  = psi  + gsl_ran_ugaussian(GSLrandom)*0.1;
  //dist_proposed = dist + gsl_ran_ugaussian(GSLrandom)*0.5;
  dist_proposed = dist * exp(gsl_ran_ugaussian(GSLrandom)*0.1); // ~10% change
  logProposalRatio *= dist_proposed / dist;
		
  copyVariables(currentParams, proposedParams);
  setVariable(proposedParams, "chirpmass",      &mc_proposed);		
  setVariable(proposedParams, "massratio",      &eta_proposed);
  setVariable(proposedParams, "inclination",    &iota_proposed);
  setVariable(proposedParams, "phase",          &phi_proposed);
  setVariable(proposedParams, "time",           &tc_proposed); 
  setVariable(proposedParams, "rightascension", &ra_proposed);
  setVariable(proposedParams, "declination",    &dec_proposed);
  setVariable(proposedParams, "polarisation",   &psi_proposed);
  setVariable(proposedParams, "distance",       &dist_proposed);

  // return ratio of proposal densities (for back & forth jumps) 
  // in "runstate->proposalArgs" vector:
  if (checkVariable(runstate->proposalArgs, "logProposalRatio"))
    setVariable(runstate->proposalArgs, "logProposalRatio", &logProposalRatio);
  else
    addVariable(runstate->proposalArgs, "logProposalRatio", &logProposalRatio, REAL8_t);
}



//Test LALPriorFunction
REAL8 BasicUniformLALPrior(LALInferenceRunState *runState, LALVariables *params)
/****************************************/
/* Returns unnormalized (!),            */
/* logarithmic (!) prior density.      	*/
/****************************************/
/* Assumes the following parameters	*/
/* exist (e.g., for TaylorT1):		*/
/* chirpmass, massratio, inclination,	*/
/* phase, time, rightascension,		*/
/* desclination, polarisation, distance.*/
/* Prior is flat if within range	*/
/****************************************/
{
  REAL8 mc, eta, iota, phi, tc, ra, dec, psi, dist;	
  REAL8 logdensity;
  
  mc   = *(REAL8*) getVariable(params, "chirpmass");		/* solar masses*/
  eta  = *(REAL8*) getVariable(params, "massratio");		/* dim-less    */
  iota = *(REAL8*) getVariable(params, "inclination");		/* radian      */
  tc   = *(REAL8*) getVariable(params, "time");			/* GPS seconds */
  phi  = *(REAL8*) getVariable(params, "phase");		/* radian      */
  ra   = *(REAL8*) getVariable(params, "rightascension");	/* radian      */
  dec  = *(REAL8*) getVariable(params, "declination");		/* radian      */
  psi  = *(REAL8*) getVariable(params, "polarisation"); 	/* radian      */
  dist = *(REAL8*) getVariable(params, "distance");		/* Mpc         */

  if(eta>0.0 && eta<=0.25 && iota>=0.0 && iota<=LAL_PI && phi>=0.0 && phi<=LAL_TWOPI 
     && ra>=0.0 && ra<=LAL_TWOPI && dec>=-LAL_PI_2 && dec<=LAL_PI_2 && psi>=0.0 && psi<=LAL_PI)	
    logdensity = 0.0;
  else
    logdensity = -HUGE_VAL;
  //TODO: should be properly normalized; pass in range via priorArgs?	

  return(logdensity);
}



//Test LALProposalFunction
void ASinOmegaTProposal(LALInferenceRunState *runState, LALVariables *proposedParams)
/****************************************/
/* Assumes the following parameters		*/
/* exist:	A, Omega					*/
/* Simply picks a new value based on	*/
/* fixed Gaussian						*/
/****************************************/
{
  REAL8 A, Omega;
  REAL8 A_proposed, Omega_proposed;
  REAL8 logProposalRatio = 0.0;  // = log(P(backward)/P(forward))
  gsl_rng * GSLrandom=runState->GSLrandom;
  LALVariables * currentParams = runState->currentParams;	

  A     = *(REAL8*) getVariable(currentParams, "A");				/* dim-less	   */
  Omega = *(REAL8*) getVariable(currentParams, "Omega");			/* rad/sec     */	

  //A_proposed=A*(1.0+gsl_ran_ugaussian(GSLrandom)*0.1);			/*mc changed by 10% */
  //Omega_proposed=Omega*(1.0+gsl_ran_ugaussian(GSLrandom)*0.01);	/*Omega changed by 0.01*/
  // (above proposals not symmetric!)
  //A_proposed     = A     + gsl_ran_ugaussian(GSLrandom) * 1e-20;   // (insert some sensible number here)
  //Omega_proposed = Omega + gsl_ran_ugaussian(GSLrandom) * 0.01;
  A_proposed     = A     * exp(gsl_ran_ugaussian(GSLrandom)*0.1);   // ~ 10% change
  logProposalRatio *= A_proposed / A;
  Omega_proposed = Omega * exp(gsl_ran_ugaussian(GSLrandom)*0.01);  // ~ 1% change
  logProposalRatio *= Omega_proposed / Omega;
  
  copyVariables(currentParams, proposedParams);
  setVariable(proposedParams, "A",     &A_proposed);		
  setVariable(proposedParams, "Omega", &Omega_proposed);

  // return ratio of proposal densities (for back & forth jumps) 
  // in "runstate->proposalArgs" vector:
  if (checkVariable(runstate->proposalArgs, "logProposalRatio"))
    setVariable(runstate->proposalArgs, "logProposalRatio", &logProposalRatio);
  else
    addVariable(runstate->proposalArgs, "logProposalRatio", &logProposalRatio, REAL8_t);
}



//Test LALPriorFunction
REAL8 ASinOmegaTPrior(LALInferenceRunState *runState, LALVariables *params)
/****************************************/
/* Prior for two-parameter				*/
/* waveform family ASinOmegaT			*/
/* Assumes the following parameters		*/
/* exist:	A, Omega					*/
/* Prior is flat if within range		*/
/****************************************/
{
  REAL8 A, Omega;
  REAL8 logdensity;
  
  A     = *(REAL8*) getVariable(params, "A");				/* dim-less	   */
  Omega = *(REAL8*) getVariable(params, "Omega");			/* rad/sec     */
  
  if ((A>0.0) & (Omega>0))
    logdensity = 0.0;
  else
    logdensity = -HUGE_VAL;

  return logdensity;
}



//Test LALEvolveOneStepFunction
void BasicMCMCOneStep(LALInferenceRunState *runState)
// Metropolis-Hastings sampler.
{
  REAL8 logPriorCurrent, logPriorProposed;
  REAL8 logLikelihoodCurrent, logLikelihoodProposed;
  LALVariables proposedParams;
  REAL8 logProposalRatio = 0.0;  // = log(P(backward)/P(forward))
  REAL8 logAcceptanceProbability;

  // current values:
  logPriorCurrent      = runState->prior(runState, runState->currentParams);
  logLikelihoodCurrent = runState->currentLikelihood;

  // generate proposal:
  proposedParams.head = NULL;
  proposedParams.dimension = 0;
  runState->proposal(runState, &proposedParams);
  if (checkVariable(runstate->proposalArgs, "logProposalRatio"))
    logProposalRatio = *(REAL8*) getVariable(runstate->proposalArgs, "logProposalRatio");

  // compute prior & likelihood:
  logPriorProposed = runState->prior(runState, &proposedParams);
  if (logPriorProposed > -HUGE_VAL)
    logLikelihoodProposed = runState->likelihood(&proposedParams, runState->data, runState->template);
  else
    logLikelihoodProposed = -HUGE_VAL;

  // determine acceptance probability:
  logAcceptanceProbability = (logLikelihoodProposed - logLikelihoodCurrent) 
                             + (logPriorProposed - logPriorCurrent)
                             + logProposalRatio;

  // accept/reject:
  if ((logAcceptanceProbability > 0) 
      || (log(gsl_rng_uniform(runState->GSLrandom)) < logAcceptanceProbability)) {   //accept
    copyVariables(&proposedParams, runState->currentParams);
    runState->currentLikelihood = logLikelihoodProposed;
  }

  destroyVariables(&proposedParams);	
}



//Test LALAlgorithm
void MCMCAlgorithm(struct tagLALInferenceRunState *runState)
{
  int i;
  REAL8 dummyR8;
  
  printf(" MCMCAlgorithm(); starting parameter values:\n");
  printVariables(runState->currentParams);
  // initialize starting likelihood value:
  runState->currentLikelihood = runState->likelihood(runstate->currentParams, runState->data, runState->template);
  // iterate:
  for(i=0; i<100; i++) {
    printf(" MCMC iteration: %d\n", i+1);
    dummyR8 = runState->currentLikelihood;
    runState->evolve(runState);
    if (runState->currentLikelihood != dummyR8) {
      printf(" accepted! new parameter values:\n");
      printVariables(runState->currentParams);
    }
  }
}


void NelderMeadEval(struct tagLALInferenceRunState *runState,
                    char **names, REAL8 *values, int dim,
                    REAL8 *logprior, REAL8 *loglikelihood)
// Auxiliary function for "NelderMeadAlgorithm()" (see below).
// Evaluates Prior & Likelihood for a given (sub-) set of parameters.
//  /!\  Side effect: alters value of "runState->currentParams" !
{
  int i;
  // copy over (subset of) values from "value" argument
  // (other parameter values, if any, remain as they are):
  for (i=0; i<dim; ++i)
    setVariable(runState->currentParams, names[i], &values[i]);
  // evaluate prior & likelihood:
  *logprior = runstate->prior(runstate, runstate->currentParams);
  if (*logprior > -HUGE_VAL)
    *loglikelihood = runState->likelihood(runstate->currentParams, runState->data, runState->template);
  else
    *loglikelihood = -HUGE_VAL;
  runState->currentLikelihood = *loglikelihood;
  // printf(" x");
  return;
}


void NelderMeadAlgorithm(struct tagLALInferenceRunState *runState, LALVariables *subset)
/************************************************************************************/
/*  Nelder-Mead (flexible polyhedron search) algorithm                              */
/*  following D. M. Himmelblau (1972): Applied nonlinear programming. McGraw-Hill.  */
/************************************************************************************/
/* Starting values are generated from the "runState->currentParams" value, by       */
/* using repeated draws from "runState->proposal()" function while ensuring         */
/* non-zero prior density for these.                                                */
/* Depending on the "ML" setting (still hard-coded, see below), it will either      */
/* aim for Maximum-Likelihood (ML) or Maximum-A-Posteriori (MAP) values.            */
/* In future, one should be able to specify the subset of parameters to be          */
/* optimized over (since e.g. one wouldn't want to optimize over PN order, which    */
/* may also be part of the parameters. Or one may want to keep sky location fixed). */
/* By now the algorithm can only handle REAL8 parameters.                           */
/************************************************************************************/
/* TO DO:                                                                           */
/*  - use (named) "subset" argument to determine subset of parameters over which to */
/*    optimize. By now simply all "REAL8" values are used.                          */
/*  - allow to specify "ML" option from outside function.                           */
/*  - allow to supply starting simplex?                                             */
/*  - allow to specify stop criteria from outside.                                  */
/*  - get rid of text output.                                                       */
/*  - somehow allow parameters like phase or rightascension to wrap around          */ 
/*    (i.e. let the simplex move across parameter space bounds)                     */
/************************************************************************************/
{
  int ML = 1; // ML==1 --> Maximum-Likelihood (ML);  ML==0 --> Maximum-A-Posteriori (MAP).
  //REAL8 e = sqrt(LAL_REAL8_EPS); // stop criterion
  REAL8 epsilon = 0.001;  // stop criterion 
  int maxiter = 500;      // also stop criterion

  int i, j;
  LALVariables param;
  LALVariables startval;
  char str[VARNAME_MAX];
  int nmDim;            // dimension of (sub-) space to be optimized over.
  char **nameVec=NULL;  // vector of dimensions' names.
  REAL8 *R8Vec=NULL;
  REAL8 *simplex=NULL;  // (d x (d+1)) - matrix containing simplex vertices.
  REAL8 *val_simplex;   // corresponding function values (likelihood or posterior).
  REAL8 logprior, loglikelihood; // dummy variables.
  REAL8 *centroid, *reflected, *expanded, *contracted; // proposed new vertices...
  REAL8 val_centroid, val_reflected, val_expanded, val_contracted; // ...and corresponding function values
  int iteration;
  int terminate=0;
  int mini, maxi;
  
  printf(" NelderMeadAlgorithm(); current parameter values:\n");
  printVariables(runState->currentParams);
  startval.head=NULL;
  startval.dimension=0;
  copyVariables(runState->currentParams, &startval);

  // initialize "param":
  param.head=NULL;
  param.dimension=0;
  // "subset" specified? If not, simply gather all REAL8 elements of "currentParams" to optimize over:
  if (subset==NULL) {
    if (runstate->currentParams == NULL) 
      die(" ERROR in NelderMeadAlgorithm(): no \"runstate->currentParams\" vector provided.\n");
    i = getVariableDimension(runstate->currentParams);
    if (i==0) 
      die(" ERROR in NelderMeadAlgorithm(): empty \"runstate->currentParams\" vector provided.\n");
    for (j=1; j<=i; ++j) {  // check "currentParams" entries and copy all REAL( values:
      if (getVariableType(runstate->currentParams, j) == REAL8_t){
	strcpy(str, getVariableName(runstate->currentParams, j));
        addVariable(&param, str, getVariable(runstate->currentParams, str), REAL8_t);
      }
    }
  }
  else {
    die(" ERROR in NelderMeadAlgorithm(): \"subset\" feature not yet implemented.\n");
    // TODO: take a (named) "subset" vector of zeroes/ones indicating which variables optimize over and which to keep fixed.
  }

  // Figure out parameter space dimension, names, &c.:
  nmDim   = getVariableDimension(&param);
  R8Vec   = (REAL8*) malloc(sizeof(REAL8) * nmDim);
  nameVec = (char**) malloc(sizeof(char*) * nmDim);
  for (i=0; i<nmDim; ++i) {
    nameVec[i] = (char*) malloc(sizeof(char) * VARNAME_MAX);
    strcpy(nameVec[i], getVariableName(&param, i+1));
    R8Vec[i] = *(REAL8*) getVariable(&param, nameVec[i]);
  }
  simplex       = (REAL8*) malloc(sizeof(REAL8) * nmDim * (nmDim+1));
  val_simplex   = (REAL8*) malloc(sizeof(REAL8) * (nmDim+1));
  centroid   = (REAL8*) malloc(sizeof(REAL8) * nmDim);
  reflected  = (REAL8*) malloc(sizeof(REAL8) * nmDim);
  expanded   = (REAL8*) malloc(sizeof(REAL8) * nmDim);
  contracted = (REAL8*) malloc(sizeof(REAL8) * nmDim);

  // populate simplex;
  // first corner is starting value:
  for (j=0; j<nmDim; ++j)
    simplex[j] = *(REAL8*) getVariable(&param, nameVec[j]);
  NelderMeadEval(runState, nameVec, &simplex[0], nmDim, &logprior, &loglikelihood);
  if (!(loglikelihood>-HUGE_VAL))
    die(" ERROR in NelderMeadAlgorithm(): invalid starting value provided.\n");
  val_simplex[0] = ML ? loglikelihood : logprior+loglikelihood;
  // remaining corners are drawn from "runState->proposal()" function:
  for (i=1; i<(nmDim+1); ++i) {  // (loop over vertices (except 1st))
    logprior = -HUGE_VAL;
    while (!(logprior > -HUGE_VAL)) {
      // draw a proposal & copy over:
      copyVariables(&startval, runState->currentParams);
      runState->proposal(runState, &param);
      for (j=0; j<nmDim; ++j)
        simplex[i*nmDim+j] = *(REAL8*) getVariable(&param, nameVec[j]);
      // compute prior & likelihood:
      NelderMeadEval(runState, nameVec, &simplex[i*nmDim], nmDim, &logprior, &loglikelihood);
      val_simplex[i] = ML ? loglikelihood : logprior+loglikelihood;
    }    
  }
  // determine minimum & maximum in simplex:
  mini = maxi = 0;
  for (i=1; i<(nmDim+1); ++i) {
    if (val_simplex[i] < val_simplex[mini]) mini = i;
    if (val_simplex[i] > val_simplex[maxi]) maxi = i;
  }

  // start actual Nelder-Mead iterations:
  iteration = 0;
  while (!terminate) {
    ++iteration;
    // determine centroid of simplex, excluding the worst (minimal) point:
    for (i=0; i<nmDim; ++i) {      // (loop over parameter dimensions)
      centroid[i] = 0.0;
      for (j=0; j<(nmDim+1); ++j)  // (loop over simplex vertices)
        centroid[i] += (j==mini) ? 0.0 : (1.0/((double)nmDim)) * simplex[j*nmDim+i];
    }
    NelderMeadEval(runState, nameVec, centroid, nmDim, &logprior, &loglikelihood);
    val_centroid = ML ? loglikelihood : logprior+loglikelihood;

    // REFLECT:
    for (i=0; i<nmDim; ++i)
      reflected[i] = centroid[i] + 1.0*(centroid[i] - simplex[mini*nmDim + i]);
    NelderMeadEval(runState, nameVec, reflected, nmDim, &logprior, &loglikelihood);
    val_reflected = ML ? loglikelihood : logprior+loglikelihood;
    if (val_reflected > val_simplex[maxi]) { // reflected better than best so far?
      // EXPAND:
      for (i=0; i<nmDim; ++i)
        expanded[i] = centroid[i] + 2.9*(reflected[i] - centroid[i]);
      NelderMeadEval(runState, nameVec, expanded, nmDim, &logprior, &loglikelihood);
      val_expanded = ML ? loglikelihood : logprior+loglikelihood;
      if (val_expanded > val_simplex[maxi]) { // expanded better than best so far?
        for (i=0; i<nmDim; ++i) // adopt expanded
          simplex[mini*nmDim+i] = expanded[i];
        val_simplex[mini] = val_expanded;
      }
      else {
        for (i=0; i<nmDim; ++i) // adopt reflected
          simplex[mini*nmDim+i] = reflected[i];
        val_simplex[mini] = val_reflected;
      }
    }
    else { // (reflected is worse that best so far)
      // check: reflected better than any of current (except worst)?
      j=0;
      for (i=0; i<(nmDim+1); ++i)
        j += ((i!=mini) && (val_reflected > val_simplex[i]));
      if (j>0) {
        for (i=0; i<nmDim; ++i) // adopt reflected
          simplex[mini*nmDim+i] = reflected[i];
        val_simplex[mini] = val_reflected;
      }
      else { // (reflected is either worst or 2nd worst)
        if (val_reflected > val_simplex[mini]) { // if 2nd worst, adopt
          for (i=0; i<nmDim; ++i) // adopt reflected
            simplex[mini*nmDim+i] = reflected[i];
          val_simplex[mini] = val_reflected;
        }
        // either way: CONTRACT:
        for (i=0; i<nmDim; ++i)
          contracted[i] = centroid[i] + 0.5*(simplex[mini*nmDim+i] - centroid[i]);
        NelderMeadEval(runState, nameVec, contracted, nmDim, &logprior, &loglikelihood);
        val_contracted = ML ? loglikelihood : logprior+loglikelihood;
        if (val_contracted > val_simplex[mini]) { // adopt contracted
          for (i=0; i<nmDim; ++i)
            simplex[mini*nmDim+i] = contracted[i];
          val_simplex[mini] = val_contracted;
        }
        else { // contraction didn't help, REDUCE:
          for (i=0; i<(nmDim+1); ++i)  // loop over vertices
            if (i!=maxi) {
              for (j=0; j<nmDim; ++j)  // loop over parameters
                simplex[i*nmDim+j] += 0.5 * (simplex[maxi*nmDim+j] - simplex[i*nmDim+j]);
              NelderMeadEval(runState, nameVec, &simplex[i*nmDim], nmDim, &logprior, &loglikelihood);
              val_simplex[i] = ML ? loglikelihood : logprior+loglikelihood;
	    }
        }
      }
    }
    // re-determine minimum & maximum:
    mini = maxi = 0;
    for (i=1; i<(nmDim+1); ++i) {
      if (val_simplex[i] < val_simplex[mini]) mini = i;
      if (val_simplex[i] > val_simplex[maxi]) maxi = i;
    }
    printf(" iter=%d,  maxi=%f,  range=%f\n", 
           iteration, val_simplex[maxi], val_simplex[maxi]-val_simplex[mini]);
    // termination condition:
    terminate = ((val_simplex[maxi]-val_simplex[mini]<epsilon) || (iteration>=maxiter));
  }
  // copy optimized value over to "runState->currentParams":
  for (j=0; j<nmDim; ++j)
    setVariable(runState->currentParams, nameVec[j], &simplex[maxi*nmDim+j]);
  runState->currentLikelihood = ML ? val_simplex[maxi] : runState->likelihood(runstate->currentParams, runState->data, runState->template);

  printf(" NelderMeadAlgorithm(); done.\n");
  printVariables(runState->currentParams);

  destroyVariables(&startval);
  destroyVariables(&param);
  free(R8Vec);
  for (i=0; i<nmDim; ++i) free(nameVec[i]);
  free(nameVec);
  free(simplex);
  free(val_simplex);
  free(centroid);
  free(reflected);
  free(expanded);
  free(contracted);
}
