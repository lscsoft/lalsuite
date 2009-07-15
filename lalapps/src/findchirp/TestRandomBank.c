/*
*  Copyright (C) 2007 B.S. Sathyaprakash, Thomas Cokelaer
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

/* ****************************************************************
 * Author: Cokelaer, T., Harry, I. and Sathyaprakash, B. S.
 * ************************************************************** */
#include "TestRandomBank.h"

/* --- version information --- */
NRCSID( BANKEFFICIENCYC, "$Id$");
RCSID(  "$Id$");

#define CVS_ID_STRING_C      "$Id$"
#define CVS_REVISION_C      "$Revision$"


static COMPLEX8Vector *strainSegment = NULL;
static COMPLEX8FrequencySeries       resp;
static int vrbflg = 0;

RandomParams *randParams = NULL;
INT4 randnStartPad = 0;
INT4 ascii2xml=0;

int
main (INT4 argc, CHAR **argv ) 
{

  /* checkpointing */
  CHAR  ckptPath[FILENAME_MAX];           /* input and ckpt file path     */
  CHAR  outputPath[FILENAME_MAX];         /* output data file path        */

  /* --- Variables ---*/
  INT4          ntrials = 0;
  INT4 	        i, thisTemplateIndex;
  Approximant   tempOrder; 

  /* --- input --- */
  UserParametersIn 	userParam;

  /* --- signal related --- */
  REAL4Vector                   signal;
  RandomInspiralSignalIn 	randIn;		

  /* --- template bank related --- */  
  InspiralTemplate      	insptmplt, randtmplt;

  InspiralCoarseBankIn      	coarseBankIn;
  INT4  	                sizeBank = 0;
  INT4                          filter_processed = 0;
   
  /* --- variables needed for the random bank */
  UINT4 randI, randN=1418;
  REAL8 randM1, randM2;
  /* ------------------ */

  MetadataTable                 templateBank;
  SnglInspiralTable            *tmpltHead  = NULL;
  SnglInspiralTable            *tmpltCurrent  = NULL;
  SnglInspiralTable            *bank;
  /* --- filtering related --- */

  REAL4Vector    correlation;
  REAL4Vector    FilterBCVSpin1;
  REAL4Vector    FilterBCVSpin2;
  REAL4Vector    FilterBCVSpin3;
  REAL4Vector    FilterBCV1;
  REAL4Vector    FilterBCV2;
  BEPowerVector  powerVector;
  BEMoments      moments;
  InspiralWaveOverlapIn 	overlapin;		
  InspiralWaveOverlapOut     	overlapout; 

  /* --- results and datamining --- */
  ResultIn              result;
  OverlapOutputIn       OverlapOutputThisTemplate;
  OverlapOutputIn       OverlapOutputBestTemplate;

  /* others */
  RealFFTPlan 		*fwdp = NULL;
  RealFFTPlan 		*revp = NULL;


  LALStatus             status = blank_status;

  FILE                  *Foutput;

  /* --- START main code --- */
  
  /* --- init some histogram structure --- */
  gsl_histogram * histogramNoise = gsl_histogram_alloc (200);  
  gsl_histogram_set_ranges_uniform (histogramNoise, 0.,20.);
    
  /* --- Some initialization --- */ 

  lal_errhandler = LAL_ERR_EXIT;
  lalDebugLevel = 0;

  templateBank.snglInspiralTable = NULL;

  /* --- Initialization of structure --- */
  ParametersInitialization(&coarseBankIn, 
		  	   &randIn,
			   &userParam);

  /* --- Read user parameters --- */
  ParseParameters(&argc,
		  argv,
		  &coarseBankIn,
		  &randIn,
		  &userParam);

  /* --- Check input parameters --- */
  UpdateParams(&coarseBankIn,
	       &randIn,
	       &userParam);
   
  if (ascii2xml==1)
    {
      BEAscii2Xml();
      exit(1);
    }

 /* zero out the checkpoint and output paths */
  memset( ckptPath, 0, FILENAME_MAX * sizeof(CHAR) );
  memset( outputPath, 0, FILENAME_MAX * sizeof(CHAR) );
  snprintf( outputPath, FILENAME_MAX * sizeof(CHAR), "bankefficiency.ckpt");
        

  /* init a random structure using possibly the seed from user input*/
  LAL_CALL(LALCreateRandomParams(&status, &randParams, randIn.useed ), 
	   &status);


  /* 
     This is used if one want to use condor script to keep track of the
     input parameters. It creates a prototype to be read by BEAscii2xml. 
     Useful with condor. See BankEfficiency documentation.
  */
  if (userParam.printPrototype){
    BEPrintProtoXml(coarseBankIn, randIn, userParam);   
    exit(0);
  }
  
  /* --- Estimate the size of the signal --- */
  LAL_CALL(BEGetMaximumSize(&status, randIn, coarseBankIn, &(signal.length)), 
	   &status);
  /* if --num-seconds is provided, we want a specified length of data.
   * Otherwise, we stick to a minimal size given by twice the length of the
   * maximal template length.*/

  if (userParam.numSeconds*randIn.param.tSampling >= signal.length) {
    signal.length = userParam.numSeconds * randIn.param.tSampling;
  }
  else if (userParam.numSeconds != -1) {
    fprintf(stderr, "you asked a length of %d seconds but a (%f %f) system  might be longer than that length...quitting \n", 
	    userParam.numSeconds, randIn.mMin, randIn.mMin);
    exit(0);
  }
  else {
    userParam.numSeconds = (INT4)(signal.length/randIn.param.tSampling);
  }



  /* --- Allocate memory for some vectors --- */
  FilterBCV1.length     = signal.length;
  FilterBCV2.length     = signal.length;
  FilterBCVSpin1.length     = signal.length;
  FilterBCVSpin2.length     = signal.length;
  FilterBCVSpin3.length     = signal.length;

  randIn.psd.length 	= signal.length/2 + 1; 
  correlation.length 	= signal.length;

  signal.data 		= (REAL4*) LALCalloc(1, sizeof(REAL4) * signal.length);
  correlation.data 	= (REAL4*) LALCalloc(1, sizeof(REAL4) * correlation.length);
  randIn.psd.data 	= (REAL8*) LALCalloc(1, sizeof(REAL8) * randIn.psd.length); 
  FilterBCV1.data       = (REAL4*) LALCalloc(1, sizeof(REAL4) * FilterBCV1.length);
  FilterBCV2.data       = (REAL4*) LALCalloc(1, sizeof(REAL4) * FilterBCV2.length);
  FilterBCVSpin1.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) * FilterBCVSpin1.length);
  FilterBCVSpin2.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) * FilterBCVSpin2.length);
  FilterBCVSpin3.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) * FilterBCVSpin3.length);
   

  /* --- create the PSd noise --- */
  if (vrbflg){
    fprintf(stdout, "generating PSD ...");
  }
  LAL_CALL(BECreatePsd(&status, &coarseBankIn, &randIn, userParam), 
	   &status);
  if (vrbflg){
    fprintf(stdout, " ... done \n");
  }

  /* --- create the template bank --- */
 if ( vrbflg )
  {
    fprintf( stdout, "generating template bank parameters... " );
    fflush( stdout );
  }
 {
   /* make sure the pointer to the first template is null */
   templateBank.snglInspiralTable = NULL;
   if (coarseBankIn.approximant == BCVSpin){
     coarseBankIn.mMin = 0;
   }
   /*
   LAL_CALL( LALInspiralBankGeneration( &status, &coarseBankIn, &tmpltHead, &sizeBank),
	     &status );
   */

   randtmplt.fLower = 40.;
   randtmplt.massChoice = m1Andm2;
   randtmplt.order = twoPN;
   randtmplt.fFinal = 1024.;
   randtmplt.beta = 0.;

   bank = (SnglInspiralTable *) LALCalloc(1, sizeof(SnglInspiralTable));
   /*
   if (bank == NULL){
	   ABORT( &status, BANKEFFICIENCY_EMEM, BANKEFFICIENCY_MSGEMEM);
   }
   */
   sizeBank = randN;
   tmpltHead = bank;
   for (randI=0; randI<randN; randI++)
   {
	   fscanf(stdin, "%le %le\n", &randM1, &randM2);
	   randtmplt.mass1 = randM1;
	   randtmplt.mass2 = randM2;
	   LALInspiralParameterCalc( &status, &randtmplt);
      
	   bank = bank->next = (SnglInspiralTable *) LALCalloc( 1, sizeof( SnglInspiralTable ) );
	   /*
	   if (bank == NULL)
	   {
		   ABORT( &status, BANKEFFICIENCY_EMEM, BANKEFFICIENCY_MSGEMEM);
	   }
	   */
	   bank->mass1 = randtmplt.mass1;
	   bank->mass2 = randtmplt.mass2;
	   bank->mchirp = randtmplt.chirpMass;
	   bank->eta = randtmplt.eta;
	   bank->tau0 = randtmplt.t0;
	   bank->tau2 = randtmplt.t2;
	   bank->tau3 = randtmplt.t3;
	   bank->tau4 = randtmplt.t4;
	   bank->tau5 = randtmplt.t5;
	   bank->ttotal = randtmplt.tC;
	   bank->psi0 = randtmplt.psi0;
	   bank->psi3 = randtmplt.psi3;
	   bank->f_final = randtmplt.fFinal;
	   bank->eta = randtmplt.eta;
	   bank->beta = randtmplt.beta;
      
	   /*
	   memcpy (bank->Gamma, coarseList[cnt].metric.Gamma, 10*sizeof(REAL4));
	   */ 
	   fprintf(stderr, "%e %e %e %e\n", bank->mass1, bank->mass2, bank->tau0, bank->tau3);
   }
   bank = tmpltHead->next;
   LALFree(tmpltHead);
   tmpltHead = bank;
	   /*
   exit(0);
	   */
   

   if ( sizeBank )
     {
       templateBank.snglInspiralTable = tmpltHead;
     }
   if (userParam.printBank)
   {
     LALBankPrintXML(templateBank, coarseBankIn, randIn, userParam);
     LALBankPrintAscii(templateBank, sizeBank, coarseBankIn);
   }
   
   if ( vrbflg )
     {
       fprintf( stdout, "done. Got %d templates\n", sizeBank );
       fflush(stdout);
     }
 }



 

  /* --- Estimate the fft's plans --- */
  LAL_CALL(LALCreateForwardRealFFTPlan(&status, &fwdp, signal.length, 0), 
	   &status);
  LAL_CALL(LALCreateReverseRealFFTPlan(&status, &revp, signal.length, 0), 
	   &status);
  
  /* --- The overlap structure --- */
  overlapin.nBegin 	= 0;
  overlapin.nEnd 	= 0;
  overlapin.psd 	= randIn.psd;
  overlapin.fwdp 	= randIn.fwdp = fwdp;
  overlapin.revp 	= revp;
  overlapin.ifExtOutput = 0;  /* new option from anand WaveOverlap */

  /* --- Create Matrix for BCV Maximization once for all --- */
  LAL_CALL( BECreatePowerVector(&status, &powerVector, randIn, signal.length), 
	    &status);
  LALCreateBCVMomentVector(&moments, &coarseBankIn.shf,  randIn.param.tSampling, 
			randIn.param.fLower, signal.length);
  
  /*
   * Main loop is here. We creae random signal and filter 
   * it with the template bank.
   */
  /* while we want to perform one simulation */
  if (vrbflg){
    PrintParameters(coarseBankIn, randIn, userParam);
  }

  while (++ntrials <= userParam.ntrials) 
    {     
      UINT4 currentTemplate = 0;

      if (vrbflg){
	fprintf(stdout,"Simulation number %d/%d\n", ntrials, userParam.ntrials);
      }
      
      /* do some initialisation for the output storage of the best temlpate */
      BEInitOverlapOutputIn(&OverlapOutputBestTemplate);

      if (vrbflg){
	fprintf(stdout,"Init ... done\n");
      }
      
      randIn.param.fCutoff = userParam.signalfFinal; 
      if (vrbflg){
        fprintf(stdout, "Upper frequency cut-off on signal: %e\n", randIn.param.fCutoff);
	fprintf(stdout,"Signal Generation ... ");
	fflush(stdout);
      }
      /* first generate some input signal. GW alone or mixed with gaussian noise. */
      /* if realNoise is requested and simulation type is at least with noise,
       * then we just add real noise to the simulated signal. */
      if (userParam.realNoise && userParam.noiseModel == REALPSD) {
        for (i=0; i<(INT4)signal.length/2; i++){
	  INT4 k = signal.length-i;
	    signal.data[i] = strainSegment->data[i].re; 
	    signal.data[k] = strainSegment->data[i].im; 
        }
      }
      else
	{

	  if (randnStartPad==1)
	    {
	      randIn.param.nStartPad = (int)((float)rand()/(float)RAND_MAX*signal.length/2);
	    }
	  LAL_CALL(BEGenerateInputData(&status,  &signal, &randIn, userParam),
		   &status);

	}

      overlapin.signal 	= signal;
       
      if (vrbflg){
	fprintf(stdout,"... done\n");
	Foutput=  fopen("BE_signal.dat","w");
	for (i=0; i<(INT4)signal.length; i++)
	  fprintf(Foutput, "%e %e\n", (float)i/randIn.param.tSampling,signal.data[i]);
	fclose(Foutput);
	

      }

      /* --- we can process the data through the bank now ---        */
      /* -- first we populate the insptmplt with the random parameter*/      
      insptmplt = randIn.param;      
      filter_processed = sizeBank;

      /* -- and enter in the bank */
      for (tmpltCurrent  = tmpltHead, thisTemplateIndex=0;
	   tmpltCurrent && thisTemplateIndex < sizeBank;
	   tmpltCurrent = tmpltCurrent->next, thisTemplateIndex++)
	{ 
		fprintf(stderr, "%d\n", thisTemplateIndex);
	  /* populate InspiralTemplateList with tmplt */

	  if (vrbflg){	    
	    fprintf(stdout,".");
	    fflush(stdout);
	  }	
	  /* Set dummy values in the future output */
	  BEInitOverlapOutputIn(&OverlapOutputThisTemplate);

	  /* which template family to be used ? */
	  switch(userParam.template)
	    {

	    case BCVSpin:
	      CreateListfromTmplt(&insptmplt, tmpltCurrent);
	      insptmplt.massChoice = psi0Andpsi3;
	      LAL_CALL(LALInspiralParameterCalc( &status,  &(insptmplt) ), &status);

	      if (userParam.faithfulness){
		insptmplt = randIn.param;
		overlapin.param              = randIn.param;
		overlapin.param.approximant  = userParam.template;
		sizeBank                     = 1;  /* we dont need the whole bank for checking */ 
	      }
	      
	      LAL_CALL(LALInspiralOverlapBCVSpin(&status, 
						 &insptmplt,
						 &powerVector,
						 &userParam,
						 &randIn, 
						 &FilterBCVSpin1,
						 &FilterBCVSpin2, 
						 &FilterBCVSpin3, 
						 &overlapin,
						 &OverlapOutputThisTemplate,
						 &correlation, &moments), 
		       &status);
	      
	      OverlapOutputThisTemplate.freq              =  overlapin.param.fFinal;
	      OverlapOutputThisTemplate.templateNumber    = thisTemplateIndex;
	      break;

	    case BCV:	   
	      CreateListfromTmplt(&insptmplt, tmpltCurrent);
	      insptmplt.massChoice = psi0Andpsi3;
	      LAL_CALL(LALInspiralParameterCalc( &status,  &(insptmplt) ), &status);
	      
	      if (userParam.faithfulness){
		insptmplt = randIn.param;
		overlapin.param              = randIn.param;
		overlapin.param.approximant  = userParam.template;
		sizeBank                     = 1;  /* we dont need the whole bank for checking */ 
	      }
	 
	      LAL_CALL(LALInspiralOverlapBCV(&status, 
					     &insptmplt,
					     &powerVector,
					     &userParam,
					     &randIn, 
					     &FilterBCV1,
					     &FilterBCV2, 
					     &overlapin,
					     &OverlapOutputThisTemplate,
					     &correlation, &moments), 
		       &status);

	      OverlapOutputThisTemplate.freq              =  overlapin.param.fFinal;
	      OverlapOutputThisTemplate.templateNumber    = thisTemplateIndex;

	      
	      break;
	    case TaylorT1:
	    case TaylorT2:
	    case TaylorT3:
	    case TaylorF1:
	    case TaylorF2:
	    case EOB:
	    case PadeT1:
	    case PadeF1:
	    case SpinTaylor:
	      CreateListfromTmplt(&insptmplt, tmpltCurrent);
	      insptmplt.massChoice=t03;
	      LAL_CALL(LALInspiralParameterCalc( &status,  &(insptmplt) ), &status);
	      
              overlapin.param 	= insptmplt	   ;
              LAL_CALL(LALInspiralParameterCalc( &status,  &(overlapin.param) ), &status);
              overlapout.max    = -1.; /* let us be sure that is has a value */
              overlapin.param.fCutoff = coarseBankIn.fUpper;
              overlapin.param.fFinal = randIn.param.tSampling/2. - 1;


              
              if (userParam.faithfulness){ /* same parameters but different templates */
               /* randIn.param.fCutoff = coarseBankIn.fUpper;*/
		tempOrder = insptmplt.order;
		insptmplt = randIn.param;
		overlapin.param                    = randIn.param;
                LAL_CALL(LALInspiralParameterCalc( &status,  &(overlapin.param) ), &status);
                overlapin.param.fCutoff = coarseBankIn.fUpper;
                overlapin.param.fFinal = randIn.param.tSampling/2. - 1;
		overlapin.param.approximant        = userParam.template;
		overlapin.param.order              = tempOrder;
		sizeBank                           = 1;  /* we dont need the whole bank for checking */ 
	      }
	        
              if (userParam.template_m1 != -1){
                /* same parameters but given by the user instead of the bank */
		tempOrder = insptmplt.order;
		insptmplt = randIn.param;
                insptmplt.mass1 = userParam.template_m1;
                insptmplt.mass2 = userParam.template_m2;
                LAL_CALL(LALInspiralParameterCalc( &status,  &(insptmplt) ), &status);
		overlapin.param                    = insptmplt;
                LAL_CALL(LALInspiralParameterCalc( &status,  &(overlapin.param) ), &status);
                overlapin.param.fCutoff = randIn.param.tSampling/2. - 1;
                overlapin.param.fFinal = randIn.param.tSampling/2. - 1;
		overlapin.param.approximant        = userParam.template;
		overlapin.param.order              = tempOrder;
		sizeBank                           = 1;  /* we dont need the whole bank for checking */ 
	      }

	      /* if we want to cut integration before the fFinal
		 overlapin.param.fCutoff = 512;
	      */ 
	      if (
		  (userParam.fastSimulation == 1 ) &&
		  		  ((fabs(randIn.param.t0 - insptmplt.t0)> .3) ||
		    (fabs(randIn.param.t3 - insptmplt.t3)> .15) )
		  )
		{		  
		  filter_processed--;
		}
	      else
		{
		  LAL_CALL(LALInspiralWaveOverlap(&status,
						  &correlation,
						  &overlapout,
						  &overlapin), &status);

		  OverlapOutputThisTemplate.rhoMax = overlapout.max;
		  OverlapOutputThisTemplate.templateNumber   = currentTemplate;
		  OverlapOutputThisTemplate.phase    = overlapout.phase;
		  OverlapOutputThisTemplate.rhoBin   = overlapout.bin;
		  OverlapOutputThisTemplate.freq     = overlapin.param.fFinal;
		  OverlapOutputThisTemplate.templateNumber    = thisTemplateIndex;
		  insptmplt.fFinal = overlapin.param.fFinal;
		  OverlapOutputThisTemplate.snrAtCoaTime 
		    =  correlation.data[randIn.param.nStartPad];
	  
		}
	      break;
	    }
	  
	  /* fill histogram of the correlation output. Useful to get a flavour
           * of the output distribution in presence of noise only for
           * instance. */
	  if (userParam.printSNRHisto){
	    for (i=0; i<(INT4)correlation.length; i++){
              /* in the unconstraint case, if alphafCut is applied,
               * correlation =-1 if alphaF > 1. Therefore we do not count it
               * in the correlation histogram*/
	      if (correlation.data[i]!=-1){
                gsl_histogram_increment(histogramNoise, correlation.data[i]);
              }
	    }		    
	  }	

		  
	  /* for the final results using only one event (loudest one) throughtout the bank.*/ 
	  KeepHighestValues(OverlapOutputThisTemplate, 
			    &OverlapOutputBestTemplate);
	  /* We can also keep the best trigger in each template. */


	  
	}/*end of  bank process*/

      
      
      
      /* Then print the maximum overlap over the whole bank and the corresponding 
       * parameter of the templates, injected signal and so on. This is the main results
       * We'll print one line for each simulation (ntrials parameter)*/
      thisTemplateIndex = 0;
      for (tmpltCurrent  = tmpltHead, thisTemplateIndex = 0;
	   tmpltCurrent && thisTemplateIndex<OverlapOutputBestTemplate.templateNumber;
	   tmpltCurrent = tmpltCurrent->next, thisTemplateIndex++)
	{
	  /*nothing to do, we just want to come back to template 
	    in list which gave the best SNR */
	}

      CreateListfromTmplt(&insptmplt, tmpltCurrent);
      GetResult(	   &status, 
			   &insptmplt,
			   randIn.param, 
			   OverlapOutputBestTemplate, 
			   &result, 
			   userParam);
      /* keep track of the number of templates actually used.*/
      result.nfast  = filter_processed;
      result.ntrial = ntrials;


      PrintResults(result, randIn);  
      if (userParam.printResultXml) {
	  BEPrintResultsXml(coarseBankIn,randIn,userParam,result);	 
	}
      
      
      /*Now we might want given the best template, to reproduce results and keep overlp, wavefom of the signal , template ...*/
      /* save signal */
      /* overlap     */
      /* save template and correlation */
      
      if (userParam.template == BCV) {
	if (userParam.printBestOverlap || userParam.printBestTemplate) {
	  userParam.extraFinalPrinting = 1;
	  LAL_CALL(LALInspiralOverlapBCV(&status, 
					 &insptmplt,
					 &powerVector,
					 &userParam,
					 &randIn, 
					 &FilterBCV1,
					 &FilterBCV2, 
					 &overlapin,
					 &OverlapOutputThisTemplate,
					 &correlation, &moments), 
		   &status);	
	}
      }

      if (userParam.template == EOB) {
	if (userParam.printBestOverlap ) {
	  Foutput=  fopen("BE_correlation.dat","w");
	  
	  for (i=0; i<(INT4)correlation.length; i++)
	    {

	      fprintf(Foutput, "%e %e\n", (float)i/randIn.param.tSampling,correlation.data[i]);
	      
	      
	    }
	  fclose(Foutput);
	  
	}
      }
    }  /*end while(trial)*/
  


  if (userParam.printSNRHisto) {
    Foutput=  fopen("BE_histo.dat","w");
    gsl_histogram_fprintf(Foutput, histogramNoise, "%f", "%g");
    fclose(Foutput);
  }


  
  /* --- destroy the plans, correlation and signal --- */

  while ( templateBank.snglInspiralTable )
  {
    tmpltHead = templateBank.snglInspiralTable;
    templateBank.snglInspiralTable = templateBank.snglInspiralTable->next;
    LALFree( tmpltHead );
  }

  LALDestroyRandomParams(&status, &randParams );
  
  LALFree(powerVector.fm5_3.data);
  LALFree(powerVector.fm2_3.data);
  LALFree(powerVector.fm1_2.data);
  LALFree(powerVector.fm7_6.data);

  LALFree(moments.a11.data);
  LALFree(moments.a21.data);
  LALFree(moments.a22.data);
   
  LALFree(FilterBCV2.data);
  LALFree(FilterBCV1.data);
  
  LALFree(randIn.psd.data);
  LALDDestroyVector( &status, &(coarseBankIn.shf.data) );
 
  LALFree(signal.data);
  LALFree(correlation.data);

  
  LALDestroyRealFFTPlan(&status,&fwdp);
  LALDestroyRealFFTPlan(&status,&revp);
  
  gsl_histogram_free(histogramNoise);
   
  LALCheckMemoryLeaks();    
   
   
  return 0;
}



		      
/* ****************************************************************************
 * Some output Results
 * ************************************************************************** */
void
KeepHighestValues(OverlapOutputIn resultThisTemplate,   
		  OverlapOutputIn *resultBestTemplate
		  )
{
  if (resultThisTemplate.rhoMax > resultBestTemplate->rhoMax){    
    resultBestTemplate->rhoMax   = resultThisTemplate.rhoMax;
    resultBestTemplate->phase    = resultThisTemplate.phase;
    resultBestTemplate->alpha    = resultThisTemplate.alpha;
    resultBestTemplate->rhoBin   = resultThisTemplate.rhoBin;
    resultBestTemplate->freq     = resultThisTemplate.freq;
    resultBestTemplate->templateNumber = resultThisTemplate.templateNumber;
    resultBestTemplate->snrAtCoaTime   = resultThisTemplate.snrAtCoaTime;
  }
  
  
  
}


/* ****************************************************************************
 * Some output Results
 * ************************************************************************** */

void 
GetResult(
	  LALStatus                  *status, 
	  InspiralTemplate           *list,
	  InspiralTemplate           injected,
	  OverlapOutputIn 	     bestOverlap, 
	  ResultIn                   *result,
	  UserParametersIn           userParam )
{
  INT4 templateNumber;
  InspiralTemplate trigger;

  
  INITSTATUS (status, "GetResult", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);
     
  templateNumber = bestOverlap.templateNumber;

  trigger = *list;

  if (userParam.template == BCV || userParam.template == BCVSpin){

    /*    LALInspiralParameterCalc( status->statusPtr,  &trigger );
    CHECKSTATUSPTR(status);							
    */
    result->psi0_inject  = injected.psi0;
    result->psi3_inject  = injected.psi3;     
    result->beta_inject  = injected.beta;     
    result->psi0_trigger  = trigger.psi0;
    result->psi3_trigger  = trigger.psi3;     
    result->beta_trigger  = trigger.beta;     
    result->tau0_trigger = 0.;
    result->tau3_trigger = 0.;
    result->tau0_inject  = 0.;
    result->tau3_inject  = 0.; 
  }
  else{
    LALInspiralParameterCalc( status->statusPtr,  &trigger );
    CHECKSTATUSPTR(status);							
    result->psi0_inject   = 0.;
    result->psi3_inject   = 0.;     
    result->beta_inject  = 0.;
    result->psi0_trigger  = 0.;
    result->psi3_trigger  = 0.;     
    result->beta_trigger  = 0.;
    result->tau0_trigger = trigger.t0;
    result->tau3_trigger = trigger.t3;
    result->tau0_inject  = injected.t0;
    result->tau3_inject  = injected.t3; 
  }

  result->mass1_inject      = injected.mass1;
  result->mass2_inject      = injected.mass2;
  result->fend_inject       = injected.fFinal;
  result->fend_trigger      = bestOverlap.freq;

  result->rho_final   = bestOverlap.rhoMax;
  result->alphaF      = bestOverlap.alpha*pow(bestOverlap.freq,2./3.);
  result->bin         = bestOverlap.rhoBin;
  result->phase       = bestOverlap.phase;
  result->layer       = bestOverlap.layer;
  result->phase       = bestOverlap.phase;
  result->snrAtCoaTime = bestOverlap.snrAtCoaTime;
  
    
  DETATCHSTATUSPTR(status);
  RETURN (status);

}


void
PrintResults(   ResultIn result, 
                RandomInspiralSignalIn randIn)
{

  fprintf(stdout, "%e %e %e %e ", 
	  result.tau0_trigger, 
	  result.tau3_trigger,
	  randIn.param.t0,
	  randIn.param.t3);
  
  fprintf(stdout, "%7.2f %7.2f   %e %e ", 
	  result.fend_trigger, 
	  randIn.param.fFinal,
	  randIn.param.mass1,
	  randIn.param.mass2);

  
  fprintf(stdout, " %e %d %d \n",
	  result.rho_final, 
	  result.bin, 
	  randIn.param.nStartPad); 

  fflush(stdout);
}





/* ****************************************************************************
 * Functio to generate and stored the moments in  REAL4vectors. 
 * ************************************************************************* */
void LALCreateBCVMomentVector(BEMoments             *moments,
			      REAL8FrequencySeries  *psd,
			      REAL8                 sampling,
			      REAL8                 fLower,
			      INT4                  length)
{
  REAL8 		m7 = 0.;							/* the moments */
  REAL8 		m5 = 0.;	
  REAL8 		m3 = 0.;
  REAL8 		f;
  UINT4			kMin;
  UINT4		  	kMax;
  UINT4			k;
  InspiralMomentsIn 	in;

  moments->a11.length =  length / 2 ;
  moments->a22.length =  length / 2 ;
  moments->a21.length =  length / 2 ;
  moments->a11.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) * moments->a11.length);
  moments->a22.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) * moments->a22.length);
  moments->a21.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) * moments->a21.length);

  /* inspiral structure */
  in.shf 	= psd; 			/* The spectrum 			*/
  in.xmin 	= fLower;       /* Lower frequency of integral?	is it correct or should we put zero? 	*/
  in.xmax 	= sampling/2.;	/* We dont' need to cut here		*/
  in.norm 	= 1./4.;		/* Some normalization			*/

  in.norm *= sampling * sampling;

  /* The minimum and maximum frequency */
  kMin = (UINT8) floor( (in.xmin - psd->f0) / psd->deltaF );
  kMax = (UINT8) floor( (in.xmax - psd->f0) / psd->deltaF );
  
  /* Skip frequency less than fLower.*/
  for (k = 0; k < kMin ; k++)
    {
      moments->a11.data[k] = 0.;
      moments->a21.data[k] = 0.;
      moments->a22.data[k] = 0.;
    }
  for ( k = kMin; k < kMax; ++k )
    {
      f = psd->f0 + (REAL8) k * psd->deltaF;

      if (psd->data->data[k])
	{
	  m7 +=  pow( f, -(7./3.) ) / psd->data->data[k] * psd->deltaF / in.norm;
	  m5 +=  pow( f, -(5./3) )  / psd->data->data[k] * psd->deltaF / in.norm;
	  m3 +=  pow( f, -(1.) )    / psd->data->data[k] * psd->deltaF / in.norm;
	  moments->a11.data[k] = 1./sqrt(m7);
	  moments->a22.data[k] = 1. / sqrt(m3 - m5*m5/m7);
	  moments->a21.data[k] = -m5 / m7 * moments->a22.data[k];
	}
    }
}


/* ****************************************************************************
 * Create an orthogonal filter. by taking its complex conjuguate
 * ************************************************************************** */
void 
LALGetOrthogonalFilter(REAL4Vector *filter)
{
  UINT4   i;
  UINT4	  n 	 = filter->length;
  UINT4   nby2   = filter->length / 2;
  REAL4	  temp;
  
  for (i = 1; i < nby2; i++)
    {
      temp 		=  filter->data[i];
      filter->data[i] 	= -filter->data[n-i];
      filter->data[n-i] = temp;
    }
}

  

/* *****************************************************************************
 * Create a vector f^{a/b} 
 * ************************************************************************** */
void 
LALCreateVectorFreqPower(REAL4Vector 		*vector,
			 InspiralTemplate 	params,
			 INT4 			a,
			 INT4 			b)
{

  INT4 	        i;
  INT4		n = vector->length;			/* Length of the vector to create 	*/
	  
  REAL8		power = (REAL8)a / (REAL8)b;		/* the frequency power			*/
  REAL8 	f;					/* La frequence				*/
  REAL8		df = params.tSampling/(REAL8)n/2.;	/* sampling frequency			*/

  /* First value equal to zero */
  vector->data[0] = 0.;
  for( i = 1; i < n; i++)
    {
      f = i * df; 
      /* Probably not efficient but this function is 
       * just called once at the beginning */
      vector->data[i] = pow(f, power);
    }
}


/* *****************************************************************************
 * Create the filters for BCV correlation 					
 * ************************************************************************** */
void 
LALCreateBCVFilters(REAL4Vector         *FilterBCV1,
		    REAL4Vector 	*FilterBCV2,
		    BEPowerVector       *powerVector,
		    BEMoments           *moments,	/* moments 				*/
		    UINT4       	kMin,		/* position of the integration		*/
		    UINT4               kMax,
		    REAL4 		psi0,
		    REAL4 		psi3)
{
  UINT4			i;
  UINT4  		n = FilterBCV1->length;
  UINT4	  		nby2 = FilterBCV1->length / 2;
	  
  REAL8 		amplitude, phase;
  REAL8	  		cos_phase;
  REAL8	  		sin_phase;

  REAL4 a11,a21,a22;

  a11 = moments->a11.data[kMax];
  a22 = moments->a22.data[kMax];
  a21 = moments->a21.data[kMax];

  /* Create the templates */	       
  for (i = 0; i< kMin-1; i++)
    {
      FilterBCV1->data[i] = 0;
      FilterBCV2->data[i] = 0;
    }
  /* Should add the upper fendBCV here as well */

  for (i = kMin; i< nby2; i++)
    {
      phase  =   psi0 * powerVector->fm5_3.data[i]	/* The same for the two filters 	*/
	+ psi3 * powerVector->fm2_3.data[i];
      
      cos_phase  = cos(phase);	/* The same for the two filters		*/
      sin_phase  = sin(phase);
      /* Fill the first filter here */
      amplitude  = a11 * powerVector->fm7_6.data[i];
      
      FilterBCV1->data[i]   =  amplitude * cos_phase;      
      FilterBCV1->data[n-i] = -amplitude * sin_phase;

      /* Fill the second one here */
      amplitude =  a21 * powerVector->fm7_6.data[i] + 
	a22 * powerVector->fm1_2.data[i];
      
      FilterBCV2->data[i]   =  amplitude * cos_phase;      
      FilterBCV2->data[n-i] = -amplitude * sin_phase;
    }  
}

void 
LALCreateBCVSpinFilters(REAL4Vector         *FilterBCVSpin1,
			REAL4Vector 	    *FilterBCVSpin2,
			REAL4Vector 	    *FilterBCVSpin3,
			InspiralWaveOverlapIn   *overlapin,
			BEPowerVector       *powerVector,
			BEMoments           *moments,	/* moments 				*/
			UINT4       	     kMin,		/* position of the integration		*/
			UINT4                kMax
			)
{
  UINT4			i;
  UINT4  		n = FilterBCVSpin1->length;
  UINT4	  		nby2 = FilterBCVSpin1->length / 2;
  
  REAL4 a11,a21,a22; 
  REAL8 phase, c_phase, s_phase, c_phase_beta, s_phase_beta;
  REAL8 I=0., J=0., K=0., L=0., M=0., temp;
  REAL4 deltaTto2by3, deltaF;
  REAL4 psi0,psi3,beta;
  REAL4 sampling = 2048;
  beta = overlapin->param.beta;
  psi3 = overlapin->param.psi3;
  psi0 = overlapin->param.psi0;


  fprintf(stderr, "%f %f %f\n", beta, psi0, psi3);


  deltaF = sampling/n;
  deltaTto2by3 = pow(deltaF, -2./3);

  a11 = moments->a11.data[kMax];
  a22 = moments->a22.data[kMax];
  a21 = moments->a21.data[kMax];
  
  /* 
     we must move that loop outside. It is independant of the 
     psi0/psi3 ... but not beta...hum difficult. But bea dmension
     is much less than psi0/psi3. We might have a 2 dimensional 
     set of vector having i as first dimension and beta as second.
  */  
  



  for ( i = kMin; i < kMax; i++ )
    {
      if (overlapin->psd.data[i]==0){
	fprintf(stderr, "error in BCVSpinFilters function\n");
	exit(0);
      }
      temp =  powerVector->fm7_6.data[i] * powerVector->fm7_6.data[i] /overlapin->psd.data[i];
      

      I += temp;
      J += temp *
	cos(beta * powerVector->fm2_3.data[i] * deltaTto2by3);
      K += temp *
	sin(beta * powerVector->fm2_3.data[i] * deltaTto2by3);
      L += temp *
	sin(2 * beta * powerVector->fm2_3.data[i] * deltaTto2by3);
      M += temp * 
	cos(2 * beta * powerVector->fm2_3.data[i] * deltaTto2by3);

    }

  I *= 4*deltaF/sampling/sampling;
  J *= 4*deltaF/sampling/sampling;
  K *= 4*deltaF/sampling/sampling;
  L *= 2*deltaF/sampling/sampling;
  M *= 2*deltaF/sampling/sampling;

  /* Create the templates */	       
  for (i = 0; i< kMin-1; i++)
    {
      FilterBCVSpin1->data[i] = 0;
      FilterBCVSpin2->data[i] = 0;
      FilterBCVSpin3->data[i] = 0;
    }
  /* Should add the upper fendBCV here as well */

  if (beta == 0)
    {
      REAL4 rootI;
      rootI = sqrt(I);
      for (i = kMin; i< kMax; i++)  
	{
	  phase  =   psi0 * powerVector->fm5_3.data[i]	
	    + psi3 * powerVector->fm2_3.data[i];
	  c_phase = cos(phase)  ;
	  s_phase = sin(phase) ;
	  

	  FilterBCVSpin1->data[i]   =  powerVector->fm7_6.data[i] * c_phase / rootI;      
	  FilterBCVSpin1->data[n-i] =  powerVector->fm7_6.data[n-i] * s_phase / rootI;
	  FilterBCVSpin2->data[i]   = 0.;
	  FilterBCVSpin2->data[n-i] = 0.;	  
	  FilterBCVSpin3->data[i]   = 0.;
	  FilterBCVSpin3->data[n-i] = 0.;	  
	}
    }
  else
    {
      REAL8 rootI, denominator, rootDenominator, numerator1, denominator1;
      rootI           = sqrt(I);
      denominator     = I*M  +  0.5*pow(I,2) - pow(J,2);
      rootDenominator = sqrt(denominator);
      numerator1      = (I*L)-(J*K);
      denominator1    =  sqrt( (0.5*pow(I,2)) -(I*M) - pow(K,2)
          -  (pow(numerator1,2)/denominator) );


      for (i = kMin; i< kMax; i++)  
	{
	  phase  =   psi0 * powerVector->fm5_3.data[i]	/* The same for the two filters 	*/
	    + psi3 * powerVector->fm2_3.data[i];
	  c_phase = cos(phase);
	  s_phase = sin(phase);

	  FilterBCVSpin1->data[i]   =  powerVector->fm7_6.data[i] * c_phase / rootI;      
	  FilterBCVSpin1->data[n-i] =  powerVector->fm7_6.data[n-i] * s_phase / rootI;
	  
	  c_phase_beta = cos( beta * powerVector->fm2_3.data[i] * deltaTto2by3);
	  s_phase_beta = sin( beta * powerVector->fm2_3.data[i] * deltaTto2by3);

	  FilterBCVSpin2->data[i]   =  powerVector->fm7_6.data[i]
	    *( c_phase_beta - J/I) * rootI * c_phase / rootDenominator;
	  FilterBCVSpin2->data[n-i]   =  powerVector->fm7_6.data[i]
	    *( c_phase_beta - J/I) *rootI *s_phase / rootDenominator;

	  FilterBCVSpin3->data[i]   =  powerVector->fm7_6.data[i]/denominator1  
	    *(
	      s_phase_beta - (K/I)  - ( numerator1 /denominator * 
					( c_phase_beta - (J/I)))
	      ) * rootI * c_phase;   
	
	  FilterBCVSpin3->data[n-i]   =  powerVector->fm7_6.data[i]/denominator1  
	    *(
	      s_phase_beta - (K/I)  - ( numerator1 /denominator * 
					( c_phase_beta - (J/I)))
	      ) * rootI * s_phase;   
	}
    }
  /*
 {
   REAL8 A1A1=0, A2A2=0, A3A3=0, A1A2=0, A1A3=0,A2A3=0;;   
   fprintf (stdout, "Checking orthonormalisation of amplitude vectors \n");
   
   for (i=kMin; i < kMax; i++)
     {
       
       A1A1 += FilterBCVSpin1->data[i] * FilterBCVSpin1->data[i] /overlapin->psd.data[i]; ;
       A2A2 += FilterBCVSpin2->data[i] * FilterBCVSpin2->data[i]  /overlapin->psd.data[i];;
       A3A3 += FilterBCVSpin3->data[i] * FilterBCVSpin3->data[i]  /overlapin->psd.data[i];;
       A1A2 += FilterBCVSpin1->data[i] * FilterBCVSpin2->data[i] /overlapin->psd.data[i]; ;
       A1A3 += FilterBCVSpin1->data[i] * FilterBCVSpin3->data[i] /overlapin->psd.data[i]; ;
       A2A3 += FilterBCVSpin2->data[i] * FilterBCVSpin3->data[i]  /overlapin->psd.data[i];;
     }
   
   A1A1 *= 4 * deltaF;
   A2A2 *= 4 * deltaF;
   A3A3 *= 4 * deltaF;
   A1A2 *= 4 * deltaF;
   A1A3 *= 4 * deltaF;
   A2A3 *= 4 * deltaF;
   
   fprintf (stdout, "A1hat cross A1hat %e\n", A1A1);
   fprintf (stdout, "A2hat cross A2hat %e\n", A2A2);
   fprintf (stdout, "A3hat cross A3hat %e\n", A3A3);
   fprintf (stdout, "A1hat cross A2hat %e\n", A1A2);
   fprintf (stdout, "A1hat cross A3hat %e\n", A1A3);
   fprintf (stdout, "A2hat cross A3hat %e\n\n", A2A3);
 }
  */
  
}



/*  ****************************************************************************
 *  The Overlap function for BCV templates 
 *  ***************************************************************************/
void
LALWaveOverlapBCV(	     LALStatus               *status,
			     REAL4Vector             *correlation,
			     InspiralWaveOverlapIn   *overlapin,
			     REAL4Vector             *FilterBCV1,
			     REAL4Vector             *FilterBCV2,
			     UserParametersIn         userParam ,
			     OverlapOutputIn          *OverlapOutput, 
			     BEMoments *moments
		  	     )
     /*  </lalVerbatim>  */
{
 
  REAL4   rhoMaxConstraint=0   ,   alphaConstraint=0,   phaseConstraint=0;
  REAL4   rhoMaxUnconstraint=0;
  INT4    rhoBinUnconstraint=0, rhoBinConstraint=0;
  REAL4   alphaUnconstraint=0,   phaseUnconstraint=0;
  REAL4   rhoConstraint=0, rhoUnconstraint=0;

  REAL4 df;

  
  REAL4   BestPhase=0,thetab,thetav=0, alphaMax, a11,a22,a21, 
	  x1_2, x2_2, x3_2, x4_2, 		/* some temporary data 	      		*/
	  V0, V1, V2, 				/* idem 				*/
	  rho		= 0,			/* the SNR 				*/
           alpha, alphaFU, fm23, fp23;
  
  REAL4Vector 
    template, 
    x1, x2, x3, x4,
    phaseV,
    rho1, rho2, rho3, 
    v0, v1, v2;					/* output of the correlations 		*/

  InspiralWaveCorrelateIn	  corrin;/* the correlation input structure 	*/
  UINT4 k,
	  i, 									
	  nBegin, 				/* beginning of valid correlation 	*/
	  nEnd, 				/* end of valid correlation 		*/
	  n = correlation->length; 		/* length of the vectors		*/
  REAL4 phi, phase, cphase, sphase, cphi, sphi;
  FILE *Foutput;

	  
  INITSTATUS (status, "LALWaveOverlapBCV", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);
  
  /* size affectation*/
  x1.length = x2.length = x3.length = x4.length = 
  v1.length = v2.length = v0.length =
    phaseV.length =  rho1.length= rho2.length = rho3.length  =template.length= n;
  
  /* memory allocation*/
  x1.data = (REAL4*) LALMalloc(sizeof(REAL4) * x1.length);
  x2.data = (REAL4*) LALMalloc(sizeof(REAL4) * x2.length);
  x3.data = (REAL4*) LALMalloc(sizeof(REAL4) * x3.length);
  x4.data = (REAL4*) LALMalloc(sizeof(REAL4) * x4.length);

  if (userParam.printBestOverlap && userParam.extraFinalPrinting){
    phaseV.data = (REAL4*) LALMalloc(sizeof(REAL4) * phaseV.length);
    template.data = (REAL4*) LALMalloc(sizeof(REAL4) * template.length);
    rho1.data = (REAL4*) LALMalloc(sizeof(REAL4) * rho1.length);
    rho2.data = (REAL4*) LALMalloc(sizeof(REAL4) * rho2.length);
    rho3.data = (REAL4*) LALMalloc(sizeof(REAL4) * rho3.length);
    v1.data = (REAL4*) LALMalloc(sizeof(REAL4) * v1.length);
    v2.data = (REAL4*) LALMalloc(sizeof(REAL4) * v2.length);
    v0.data = (REAL4*) LALMalloc(sizeof(REAL4) * v0.length);
  }
  /**/
  overlapin->param.nStartPad 	= 0;
  overlapin->param.startPhase 	= 0;
  
  /*  We want to compute <x|h1> and <x|h2>; let's prepare the correlation  
   *  process; let's get the input data namely x and the spectrum namely psd 
   *  */
  corrin.fCutoff        = overlapin->param.fFinal;				/* the integral's cutoff frequency 	*/
  corrin.samplingRate   = overlapin->param.tSampling;				/* the sampling 			*/
  corrin.df             = overlapin->param.tSampling / 
	  		(REAL8) overlapin->signal.length;			/* the resolution in frequency 		*/
  corrin.psd            = overlapin->psd;					/* the spectrum 			*/
  corrin.signal1        = overlapin->signal;					/* the signal x				*/
  corrin.revp           = overlapin->revp;					/* fourier transform flag		*/
  
  /* Filter h1 and its orthonormal h1* */  
  corrin.signal2        = *FilterBCV1;						/* The first template 			*/
  LALInspiralWaveCorrelate(status->statusPtr, &x1, corrin);			/* the correlation			*/
  CHECKSTATUSPTR(status);							
  LALGetOrthogonalFilter(FilterBCV1);					/* get the orthonormalized template     */
  corrin.signal2        = *FilterBCV1;						/* the second template 		*/
  LALInspiralWaveCorrelate(status->statusPtr, &x3, corrin);			/* the correlation 			*/
  CHECKSTATUSPTR(status);
  
  /* Filter h2 and its orthornormal filter */
  corrin.signal2        = *FilterBCV2;						/* The first template			*/
  LALInspiralWaveCorrelate(status->statusPtr, &x2, corrin);			/* the correlation 			*/
  CHECKSTATUSPTR(status);		
  LALGetOrthogonalFilter( FilterBCV2);					/* get the orthonormailzed templates 	*/
  corrin.signal2        = *FilterBCV2;						/* the second template			*/
  LALInspiralWaveCorrelate(status->statusPtr, &x4, corrin);			/* the correlation 			*/
  
  /* nbegin and end for the correlation processus 
   * as well as overlapout affectation */
  nBegin                = overlapin->nBegin;					/* starting bin 			*/
  nEnd              	= FilterBCV1->length - overlapin->nEnd;			/* ending valid bin 			*/


  /* some output for debugging, checking ... */
  if (userParam.printBestOverlap && userParam.extraFinalPrinting)
    {
      Foutput = fopen( "BE_Filter.dat", "w");
      for (i=0; i< x1.length; i++)
	fprintf(Foutput, "%e %e\n", i/corrin.samplingRate,fabs(x1.data[i]));
      fprintf(Foutput,"&\n");	       
      
      for (i=1; i< FilterBCV1->length; i++)
	fprintf(Foutput, "%e %e\n", i/corrin.samplingRate,fabs(x2.data[i]));
      fprintf(Foutput,"&\n");	       
      
      for (i=0; i< FilterBCV1->length; i++)
	fprintf(Foutput, "%e %e\n", i/corrin.samplingRate,fabs(x3.data[i]));
      fprintf(Foutput,"&\n");	       
      
      for (i=0; i< FilterBCV1->length-1; i++)
	fprintf(Foutput, "%e %e\n", i/corrin.samplingRate ,fabs(x4.data[i]));
      fprintf(Foutput,"&\n");	 
      fclose(Foutput);      
    }

  df 	= overlapin->param.tSampling / (REAL8)(n/2.) / 2.;      
  k = floor(overlapin->param.fFinal / df);

  a11 = moments->a11.data[k];
  a22 = moments->a22.data[k];
  a21 = moments->a21.data[k];


  fm23 = pow(corrin.fCutoff,-2./3.);
  fp23 = pow(corrin.fCutoff,2./3.);

  alphaMax = fm23;

  thetab = (-(a11 * alphaMax)/(a22+a21*alphaMax));
  thetab = atan(thetab);
  rho    = 0.;
  
  if(userParam.printBestOverlap && userParam.extraFinalPrinting){
  fprintf(stderr,"theta_b = %e a11=%e a21=%e a22=%e alphaMax = %e fCutoff=%e\n",
	  thetab,  a11, a21, a22, alphaMax, corrin.fCutoff);

  }





  /* Once we have the 4 correlations, we can search for the
   * maximum of correlation and compute the final SNR */
  for (i = nBegin; i < nEnd; i++) 
    {
      x1_2 = x1.data[i] * x1.data[i];
      x2_2 = x2.data[i] * x2.data[i];
      x3_2 = x3.data[i] * x3.data[i];
      x4_2 = x4.data[i] * x4.data[i];      
      
      V0 = x1_2 + x2_2 + x3_2 + x4_2;
      V1 = x1_2 + x3_2 - x2_2 - x4_2;
      V2 = 2*(x1.data[i]*x2.data[i] + x3.data[i]*x4.data[i]);
      
      rhoUnconstraint = sqrt((V0 + sqrt(V1*V1+V2*V2))/2.);

      /* get thetav first in order to constraint alpha_F*/

      thetav =   atan2(V2,V1);
	
      if(userParam.printBestOverlap && userParam.extraFinalPrinting){
	rho1.data[i] = rhoUnconstraint;
	rho2.data[i] = sqrt((V0 + V1)/2.);
	rho3.data[i] = sqrt((V0+V1*cos(2*thetab)+V2*sin(2*thetab))/2.);
	v0.data[i]   = V0;
	v1.data[i]   = V1;
	v2.data[i]   = V2;
        phaseV.data[i] = .5*thetav; /*used to compute alpha*/
      }
             

      if (thetab >= 0){
	if ( 0  <= thetav && thetav <= 2 * thetab){	  
	  rhoConstraint = rhoUnconstraint;
	}
	else if (thetab-LAL_PI <= thetav && thetav < 0) {
	  rhoConstraint = sqrt((V0 + V1)/2.);	  
	}
	else if( (2*thetab  < thetav && thetav<=LAL_PI+1e-4 )
		 || (-LAL_PI-1e-4<=thetav && thetav < -LAL_PI+thetab)){
	  rhoConstraint =sqrt((V0+V1*cos(2*thetab)+V2*sin(2*thetab))/2.);;
	}
	else
	  {
	    fprintf(stderr,"must not enter here  thetav = %e thetab=%e\n ",
		    thetav , thetab);
	    exit(0);
	  }
      }
      else{
	if ( 2*thetab  <= thetav && thetav <= 0){	  
	  rhoConstraint = rhoUnconstraint;
	}
	else if (0 < thetav &&  thetav  <= LAL_PI + thetab ) {
	  rhoConstraint = sqrt((V0 + V1)/2.);	  
	}
	else if( (-LAL_PI-1e-4  <= thetav && thetav < 2*thetab ) || 
		 (LAL_PI +thetab <= thetav && thetav <= LAL_PI+1e-4))
        {
	  rhoConstraint =sqrt((V0+V1*cos(2*thetab)+V2*sin(2*thetab))/2.);;
	}
	else 
	  {
	    fprintf(stderr,"must not enter herethetav = %e thetab=%e %e %e %d\n ",
                thetav , thetab, V1, V2, i);
	    exit(0);
	  }
      }
    
  
      if (userParam.alphaFConstraint == ALPHAFConstraint){
  
	if(userParam.printBestOverlap || userParam.printSNRHisto) {
	  correlation->data[i] = rhoConstraint;
	}
      }
      else 
	{
	  if(userParam.printBestOverlap || userParam.printSNRHisto) {
                alphaFU =   -(a22 * tan(.5*thetav)) 				
        	    / (a11 + a21* tan(.5*thetav)) *fp23;
	        if (alphaFU <=1 )
                  correlation->data[i] = rhoUnconstraint;
                else 
                  correlation->data[i] = -1;                  
	  }
	}
      

      if ( rhoConstraint > rhoMaxConstraint)			
	/* Keep the position of the max only	*/
	{
	  rhoMaxConstraint 	  = rhoConstraint;
	  rhoBinConstraint 	  = i;
	  phaseConstraint   =    .5*thetav;
	  alphaConstraint = -(a22 * tan(phaseConstraint)) 				
	    / (a11 + a21* tan(phaseConstraint));
	}

      alphaFU =   -(a22 * tan(.5*thetav)) 				
	    / (a11 + a21* tan(.5*thetav)) *fp23;


      if ( rhoUnconstraint > rhoMaxUnconstraint && alphaFU<=1)			
	/* Keep the position of the max only	*/
	{
	  rhoMaxUnconstraint 	= rhoUnconstraint  ;
	  rhoBinUnconstraint 	= i;
	  phaseUnconstraint   = .5*thetav ; 
	  alphaUnconstraint = -(a22 * tan(phaseUnconstraint)) 				
	    / (a11 + a21* tan(phaseUnconstraint));
	}
    }
  /* Finally get the alpha value corresponding to the best rho */ 
  OverlapOutput->rhoMax     = rhoMaxConstraint;
  OverlapOutput->rhoBin     = rhoBinConstraint;
  OverlapOutput->alpha      = alphaConstraint;
  OverlapOutput->phase      = phaseConstraint;


  /* The final template to be print ? */
  if (userParam.printBestOverlap  && userParam.extraFinalPrinting){

    /*print overlap, combinaison of the x_i filters*/
    Foutput=fopen("BE_Phase.dat","w");
    for (i=0; i< phaseV.length; i++){
      fprintf(Foutput, "%e\n",  atan2(v2.data[i],v1.data[i]));
    }
    fclose(Foutput);
    
    Foutput=fopen("BE_rho1.dat","w");
    for (i=0; i< phaseV.length; i++){
      fprintf(Foutput, "%e\n", rho1.data[i]);
    }
    fclose(Foutput);
    

  Foutput=fopen("BE_rho2.dat","w");
    for (i=0; i< phaseV.length; i++){
      fprintf(Foutput, "%e\n", rho2.data[i]);
    }
    fclose(Foutput);
    

    Foutput=fopen("BE_rho3.dat","w");
    for (i=0; i< phaseV.length; i++){
      fprintf(Foutput, "%e\n",rho3.data[i]);
    }
    fclose(Foutput);
     
    Foutput=fopen("BE_v0.dat","w");
    for (i=0; i< phaseV.length; i++){
      fprintf(Foutput, "%e\n",v0.data[i]);
    }
    fclose(Foutput);

    Foutput=fopen("BE_v1.dat","w");
    for (i=0; i< phaseV.length; i++){
      fprintf(Foutput, "%e\n", v1.data[i]);
    }
    fclose(Foutput); 

    Foutput=fopen("BE_v2.dat","w");
    for (i=0; i< phaseV.length; i++){
      fprintf(Foutput, "%e\n", v2.data[i]);
    }
    fclose(Foutput);
    
    Foutput=fopen("BE_alpha.dat","w");
    for (i=0; i< correlation->length; i++){
      alpha = -(a22 * tan(phaseV.data[i])) 				/* Compute the final alpha parameter 	*/
	/ (a11 + a21* tan(phaseV.data[i]));
      fprintf(Foutput, "%e \n",alpha*fp23);
    }
    fclose(Foutput);

    /*print overlap, combinaison of the x_i filters*/
    Foutput=fopen("BE_Overlap.dat","w");
    for (i=0; i< correlation->length; i++){
      fprintf(Foutput, "%e\n", fabs(correlation->data[i]));
    }
    fclose(Foutput);
  }


    if (userParam.extraFinalPrinting && userParam.printBestOverlap){
    /* print the final template given phase, alpha ans so on*/
    
    phase 		= BestPhase;
    phi   		= overlapin->param.startPhase; /*how to get it ohterwise ? */
   
    for (i=0; i<(UINT4)correlation->length; i++) 
    {
      cphase = cos(phase);
      cphi   = cos(phi);
      sphase = sqrt(1-cphase*cphase);
      sphi   = sqrt(1-cphi*cphi);
      
      correlation->data[i] = x1.data[i] * cphase * cphi
	+ x2.data[i] * sphase * cphi
	+ x3.data[i] * cphase * sphi
	+ x4.data[i] * sphase * sphi;	
    }

    Foutput=fopen("BE_Correlation.dat", "w");
    for (i=0; i< correlation->length; i++){
      fprintf(Foutput, "%e %e\n", i/corrin.samplingRate,fabs(correlation->data[i]));
    }
    fclose(Foutput);


    LALGetOrthogonalFilter(FilterBCV1);						/* get the orthonormalized template     */
    corrin.signal2        = *FilterBCV1;

    for (i=0; i<(UINT4)correlation->length; i++) 
    { 
      cphase = cos(phaseV.data[i]);
      cphi   = cos(phi);
      sphase = sqrt(1-cphase*cphase);
      sphi   = sqrt(1-cphi*cphi);
      
      template.data[i]=FilterBCV1->data[i]*cphase*cphi; 
      template.data[i]+=corrin.signal2.data[i] * sphase*cphi; 

	
    }

    LALGetOrthogonalFilter(FilterBCV2);						/* get the orthonormalized template     */
    corrin.signal2        = *FilterBCV2;

    for (i=0; i<(UINT4)correlation->length; i++) 
    { 
      cphase = cos(phaseV.data[i]);
      cphi   = cos(phi);
      sphase = sqrt(1-cphase*cphase);
      sphi   = sqrt(1-cphi*cphi);
      
      template.data[i]=FilterBCV2->data[i]*sphase*cphi; 
      template.data[i]+=corrin.signal2.data[i]*sphase*sphi; 

	
    }

    Foutput=fopen("BE_BestTemplate.dat", "w");
    for (i=0; i< correlation->length; i++){
      fprintf(Foutput, "%e %e\n", i/corrin.samplingRate,template.data[i]);
    }
    fclose(Foutput);



  }/*end if printTemplate*/

  /* Free memory */
  LALFree(x1.data);
  LALFree(x2.data);
  LALFree(x3.data);
  LALFree(x4.data);

  if (userParam.extraFinalPrinting && userParam.printBestOverlap){
    LALFree(phaseV.data);
    LALFree(template.data);
    LALFree(rho1.data);
    LALFree(rho2.data);
    LALFree(rho3.data);
    LALFree(v0.data);
    LALFree(v1.data);
    LALFree(v2.data);
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);
}



/*  ****************************************************************************
 *  The Overlap function for BCV templates 
 *  ***************************************************************************/
void
LALWaveOverlapBCVSpin(	     LALStatus               *status,
			     REAL4Vector             *correlation,
			     InspiralWaveOverlapIn   *overlapin,
			     REAL4Vector             *FilterBCVSpin1,
			     REAL4Vector             *FilterBCVSpin2,
			     REAL4Vector             *FilterBCVSpin3,
			     UserParametersIn         userParam ,
			     OverlapOutputIn          *OverlapOutput 
		  	     )
     /*  </lalVerbatim>  */
{
 
  REAL4   rhoMax=0 ;
  INT4    rhoBin=0;
  REAL4   rho=0;

  REAL4   df;

  
  REAL4
    a11,a22,a21, 
    x5_2,x6_2,x1_2, x2_2, x3_2, x4_2, 		/* some temporary data 	      		*/
    alpha, fm23, fp23, phase;
  
  REAL4Vector 
    template, 
    x1, x2, x3, x4,x5,x6;

  InspiralWaveCorrelateIn	  corrin;/* the correlation input structure 	*/
  UINT4 k,
	  i, 									
	  nBegin, 				/* beginning of valid correlation 	*/
	  nEnd, 				/* end of valid correlation 		*/
	  n = correlation->length; 		/* length of the vectors		*/
  FILE *Foutput;

	  
  INITSTATUS (status, "LALWaveOverlapBCV", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);
  
  /* size affectation*/
  x1.length = x2.length = x3.length = x4.length = 
  x5.length = x6.length =  template.length= n;
  
  /* memory allocation*/
  x1.data = (REAL4*) LALMalloc(sizeof(REAL4) * x1.length);
  x2.data = (REAL4*) LALMalloc(sizeof(REAL4) * x2.length);
  x3.data = (REAL4*) LALMalloc(sizeof(REAL4) * x3.length);
  x4.data = (REAL4*) LALMalloc(sizeof(REAL4) * x4.length);
  x5.data = (REAL4*) LALMalloc(sizeof(REAL4) * x5.length);
  x6.data = (REAL4*) LALMalloc(sizeof(REAL4) * x6.length);


  /**/
  overlapin->param.nStartPad 	= 0;
  overlapin->param.startPhase 	= 0;
  
  /*  We want to compute <x|h1> and <x|h2>; let's prepare the correlation  
   *  process; let's get the input data namely x and the spectrum namely psd 
   *  */
  corrin.fCutoff        = overlapin->param.fFinal;				/* the integral's cutoff frequency 	*/
  corrin.samplingRate   = overlapin->param.tSampling;				/* the sampling 			*/
  corrin.df             = overlapin->param.tSampling / 
	  		(REAL8) overlapin->signal.length;			/* the resolution in frequency 		*/
  corrin.psd            = overlapin->psd;					/* the spectrum 			*/
  corrin.signal1        = overlapin->signal;					/* the signal x				*/
  corrin.revp           = overlapin->revp;					/* fourier transform flag		*/
  
  /* Filter h1 and its orthonormal h1* */  
  corrin.signal2        = *FilterBCVSpin1;						/* The first template 			*/
  LALInspiralWaveCorrelate(status->statusPtr, &x1, corrin);			/* the correlation			*/
  CHECKSTATUSPTR(status);							
  LALGetOrthogonalFilter(FilterBCVSpin1);					/* get the orthonormalized template     */
  corrin.signal2        = *FilterBCVSpin1;						/* the second template 		*/
  LALInspiralWaveCorrelate(status->statusPtr, &x2, corrin);			/* the correlation 			*/
  CHECKSTATUSPTR(status);
  
  /* Filter h2 and its orthornormal filter */
  corrin.signal2        = *FilterBCVSpin2;						/* The first template			*/
  LALInspiralWaveCorrelate(status->statusPtr, &x3, corrin);			/* the correlation 			*/
  CHECKSTATUSPTR(status);		
  LALGetOrthogonalFilter( FilterBCVSpin2);					/* get the orthonormailzed templates 	*/
  corrin.signal2        = *FilterBCVSpin2;						/* the second template			*/
  LALInspiralWaveCorrelate(status->statusPtr, &x4, corrin);			/* the correlation 			*/
  
  /* Filter h3 and its orthornormal filter */
  corrin.signal2        = *FilterBCVSpin3;						/* The first template			*/
  LALInspiralWaveCorrelate(status->statusPtr, &x5, corrin);			/* the correlation 			*/
  CHECKSTATUSPTR(status);		
  LALGetOrthogonalFilter( FilterBCVSpin3);					/* get the orthonormailzed templates 	*/
  corrin.signal2        = *FilterBCVSpin3;						/* the second template			*/
  LALInspiralWaveCorrelate(status->statusPtr, &x6, corrin);			/* the correlation 			*/

  /* nbegin and end for the correlation processus 
   * as well as overlapout affectation */
  nBegin                = overlapin->nBegin;					/* starting bin 			*/
  nEnd              	= FilterBCVSpin1->length - overlapin->nEnd;			/* ending valid bin 			*/


  /* some output for debugging, checking ... */

    {
      Foutput = fopen( "BE_Filter.dat", "w");
      for (i=0; i< x1.length; i++)
	fprintf(Foutput, "%e %e\n", i/corrin.samplingRate,fabs(x1.data[i]));
      fprintf(Foutput,"&\n");	       
      
      for (i=1; i< FilterBCVSpin1->length; i++)
	fprintf(Foutput, "%e %e\n", i/corrin.samplingRate,fabs(x2.data[i]));
      fprintf(Foutput,"&\n");	       
      
      for (i=0; i< FilterBCVSpin1->length; i++)
	fprintf(Foutput, "%e %e\n", i/corrin.samplingRate,fabs(x3.data[i]));
      fprintf(Foutput,"&\n");	       
      
      for (i=0; i< FilterBCVSpin1->length-1; i++)
	fprintf(Foutput, "%e %e\n", i/corrin.samplingRate ,fabs(x4.data[i]));
      fprintf(Foutput,"&\n");

      for (i=0; i< FilterBCVSpin1->length-1; i++)
	fprintf(Foutput, "%e %e\n", i/corrin.samplingRate ,fabs(x5.data[i]));
      fprintf(Foutput,"&\n");	 
      
      for (i=0; i< FilterBCVSpin1->length-1; i++)
	fprintf(Foutput, "%e %e\n", i/corrin.samplingRate ,fabs(x6.data[i]));
      fprintf(Foutput,"&\n");	 



    }

  df 	= overlapin->param.tSampling / (REAL8)(n/2.) / 2.;      
  k = floor(overlapin->param.fFinal / df);


  /* Once we have the 4 correlations, we can search for the
   * maximum of correlation and compute the final SNR */
  for (i = nBegin; i < nEnd; i++) 
    {
      x1_2 = x1.data[i] * x1.data[i];
      x2_2 = x2.data[i] * x2.data[i];
      x3_2 = x3.data[i] * x3.data[i];
      x4_2 = x4.data[i] * x4.data[i];      
      x5_2 = x5.data[i] * x5.data[i];      
      x6_2 = x6.data[i] * x6.data[i];      

      rho = sqrt(x1_2+x2_2+x3_2+x4_2+x5_2+x6_2);
      fprintf(Foutput, "%e %e\n", i/corrin.samplingRate ,rho);

 /*      fprintf(stderr, "%e %e %e\n", (float)i/overlapin->param.tSampling, rho, x1_2/rho/rho+x2_2/rho/rho+x3_2/rho/rho+x4_2/rho/rho+x5_2/rho/rho+x6_2/rho/rho );*/
      if ( rho > rhoMax)			
	/* Keep the position of the max only	*/
	{
	  rhoMax 	  = rho;
	  rhoBin 	  = i;
	}


    }      
  fclose(Foutput);      

  /* Finally get the alpha value corresponding to the best rho */ 
  OverlapOutput->rhoMax     = rhoMax;
  OverlapOutput->rhoBin     = rhoBin;
  OverlapOutput->alpha      = alpha;
  OverlapOutput->phase      = phase;

  /* Free memory */
  LALFree(x1.data);
  LALFree(x2.data);
  LALFree(x3.data);
  LALFree(x4.data);
  LALFree(x5.data);
  LALFree(x6.data);


  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/* *****************************************************************************
 *  A draft function to finalized later ; purpose: print the bank 
 *  ****************************************************************************/
void
BEPrintBank(	InspiralCoarseBankIn coarseBankIn,
		InspiralTemplateList **list,
      		UINT4 sizeBank
		) 
{
  UINT4 	i;
  FILE 		*output;

  output = fopen(BANKEFFICIENCY_PRINTBANK_FILEASCII,"w");

  fprintf(output,"#Number of Coarse Bank Templates=%d\n",sizeBank);
  if (coarseBankIn.approximant == BCV)
    {
      fprintf(output, "#psi0Min=%e, psi0Max=%e, psi3Min=%e, psi3Max=%e\n", 
	      coarseBankIn.psi0Min, coarseBankIn.psi0Max, coarseBankIn.psi3Min, coarseBankIn.psi3Max);
      fprintf(output, "#psi0 psi3  totalMass fFinal\n");
      
    }
  else
    {
      fprintf(output, "#mMin=%e, mMax=%e\n", coarseBankIn.mMin,coarseBankIn.mMax);
      fprintf(output, "#tau0, tau3, mass1, mass2\n");
    }

  for ( i = 0; i < sizeBank; i++)
    {
      if (coarseBankIn.approximant == BCV)
	{      		   
      		fprintf(output, "%e %e %e %e\n", 	(*list)[i].params.psi0, 
				  			(*list)[i].params.psi3,
							(*list)[i].params.totalMass,
						       	(*list)[i].params.fFinal);
	} 
      else
	{

      		fprintf(output, "%e %e %e %e\n", 	(*list)[i].params.t0, 
				  			(*list)[i].params.t3, 
							(*list)[i].params.mass1, 
							(*list)[i].params.mass2);
	}
    }	   
  fprintf(output,"&\n");
  fclose(output);
}


#if 0
  if ( numCoarse )
    {
      /* --- we sae one line --- */
      tmplt = templateBank.snglInspiralTable = (SnglInspiralTable *)
	LALCalloc( 1, sizeof(SnglInspiralTable) );
      snprintf( tmplt->ifo, LIGOMETA_IFO_MAX * sizeof(CHAR), ifo );
      snprintf( tmplt->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR), 
		   "BankEfficiency" );
      snprintf( tmplt->channel, LIGOMETA_CHANNEL_MAX * sizeof(CHAR),
		   channelName );
      tmplt->mass1   = (REAL4) coarseList[0].params.mass1;
      tmplt->mass2   = (REAL4) coarseList[0].params.mass2;
      tmplt->mchirp  = (REAL4) coarseList[0].params.chirpMass;
      tmplt->eta     = (REAL4) coarseList[0].params.eta;
      tmplt->tau0    = (REAL4) coarseList[0].params.t0;
      tmplt->tau3    = (REAL4) coarseList[0].params.t3;
      tmplt->tau4    = (REAL4) coarseList[0].params.t4;
      tmplt->tau5    = (REAL4) coarseList[0].params.t5;
      tmplt->ttotal  = (REAL4) coarseList[0].params.tC;
      tmplt->psi0    = (REAL4) coarseList[0].params.psi0;
      tmplt->psi3    = (REAL4) coarseList[0].params.psi3;
      tmplt->f_final = (REAL4) coarseList[0].params.fFinal;
      
      /* --- then the others --- */
      for ( i = 1; i < numCoarse; ++i )
	{
	  tmplt = tmplt->next = (SnglInspiralTable *)
	    LALCalloc( 1, sizeof(SnglInspiralTable) );
	  snprintf( tmplt->ifo, LIGOMETA_IFO_MAX * sizeof(CHAR), ifo );
	  snprintf( tmplt->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR), 
		       "BankEfficiency" );
	  snprintf( tmplt->channel, LIGOMETA_CHANNEL_MAX * sizeof(CHAR),
		       channelName );
	  tmplt->mass1   = (REAL4) coarseList[i].params.mass1;
	  tmplt->mass2   = (REAL4) coarseList[i].params.mass2;
	  tmplt->mchirp  = (REAL4) coarseList[i].params.chirpMass;
	  tmplt->eta     = (REAL4) coarseList[i].params.eta;
	  tmplt->tau0    = (REAL4) coarseList[i].params.t0;
	  tmplt->tau3    = (REAL4) coarseList[i].params.t3;
	  tmplt->tau4    = (REAL4) coarseList[i].params.t4;
	  tmplt->tau5    = (REAL4) coarseList[i].params.t5;
	  tmplt->ttotal  = (REAL4) coarseList[i].params.tC;
	  tmplt->psi0    = (REAL4) coarseList[i].params.psi0;
	  tmplt->psi3    = (REAL4) coarseList[i].params.psi3;
	  tmplt->f_final = (REAL4) coarseList[i].params.fFinal;
	}
    }
  
#endif



void 
BEFillProc(
	       ProcessParamsTable     *this_proc_param,
	       InspiralCoarseBankIn   coarseBankIn,
	       RandomInspiralSignalIn randIn,
	       UserParametersIn           userParam)
{


#define ADD_PROCESS_PARAM( pptype, format, ppname, ppvalue ) \
this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
  calloc( 1, sizeof(ProcessParamsTable) ); \
  snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME ); \
  snprintf( this_proc_param->param,   LIGOMETA_PARAM_MAX,   "%-20s", ppname); \
  snprintf( this_proc_param->type,    LIGOMETA_TYPE_MAX,    "%-10s",   pptype ); \
  snprintf( this_proc_param->value,   LIGOMETA_VALUE_MAX,   format, ppvalue );


#define ADD_2PROCESS_PARAM( pptype, format, ppname, ppvalue1, ppvalue2 ) \
this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
  calloc( 1, sizeof(ProcessParamsTable) ); \
  snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME ); \
  snprintf( this_proc_param->param,   LIGOMETA_PARAM_MAX,   "%-20s", ppname); \
  snprintf( this_proc_param->type,    LIGOMETA_TYPE_MAX,    "%-10s",   pptype ); \
  snprintf( this_proc_param->value,   LIGOMETA_VALUE_MAX,   format, ppvalue1, ppvalue2 );

  ADD_PROCESS_PARAM("float","%f","--bank-alpha",           
      coarseBankIn.alpha);
  ADD_2PROCESS_PARAM("float","%f %f","--bank-fcut-range",
      coarseBankIn.LowGM, coarseBankIn.HighGM);
  ADD_PROCESS_PARAM("float","%d ","--bank-inside-polygon",
      coarseBankIn.insidePolygon);
  ADD_PROCESS_PARAM("float","%f","--bank-ffinal",
      coarseBankIn.fUpper);
  ADD_2PROCESS_PARAM("float","%f %f","--bank-mass-range",
      coarseBankIn.mMin, coarseBankIn.mMax);
  ADD_2PROCESS_PARAM("float","%f %f","--bank-psi0-range",
      coarseBankIn.psi0Min, coarseBankIn.psi0Max);
  ADD_2PROCESS_PARAM("float","%f %f","--bank-psi3-range",
      coarseBankIn.psi3Min, coarseBankIn.psi3Max);
  ADD_PROCESS_PARAM("string","%s","--bank-grid-spacing",
      GetStringFromGridType(coarseBankIn.gridSpacing));
  ADD_PROCESS_PARAM("string", "%s","--channel",
      userParam.chanName);    
  ADD_PROCESS_PARAM("string", "%s","--detector",
      GetStringFromDetector(userParam.detector));
  ADD_PROCESS_PARAM("float","%f","--fl",	
      coarseBankIn.fLower);
  ADD_PROCESS_PARAM("string", "%d","--gps-start-time",
      userParam.startTime);     
  ADD_PROCESS_PARAM("float","%f","--fl-template", 
      coarseBankIn.fLower);
  ADD_PROCESS_PARAM("float","%f","--max-total-mass",
      userParam.maxTotalMass);
  ADD_PROCESS_PARAM("float","%f","--m1",	
      userParam.m1);
  ADD_PROCESS_PARAM("float","%f","--m2",	
      userParam.m2);
  ADD_PROCESS_PARAM("float","%f","--mm",	
      coarseBankIn.mmCoarse);
  ADD_PROCESS_PARAM("int","%d","--ntrial",	
      userParam.ntrials);
  ADD_PROCESS_PARAM("float","%f","--noise-amplitude",
      randIn.NoiseAmp);
  ADD_PROCESS_PARAM("string", "%s","--noise-model",
      GetStringFromNoiseModel(userParam.noiseModel));
  ADD_PROCESS_PARAM("int","%d","--bank-number-fcut",	        
      coarseBankIn.numFcutTemplates);
  ADD_PROCESS_PARAM("int","%d","--num-seconds",
      userParam.numSeconds);
  ADD_PROCESS_PARAM("float","%f","--psi0",	
      userParam.psi0);
  ADD_PROCESS_PARAM("float","%f","--psi3",
      userParam.psi3);
  ADD_PROCESS_PARAM("string", "%s","--run",
      GetStringFromScientificRun(userParam.run));
  ADD_PROCESS_PARAM("float","%f","--sampling",
      coarseBankIn.tSampling);
  ADD_PROCESS_PARAM("string","%s","--simulation-type",
      GetStringFromSimulationType(randIn.type));
  ADD_PROCESS_PARAM("float","%f","--signal-amplitude",
      randIn.SignalAmp);
  ADD_PROCESS_PARAM("float","%f","--signal-alpha",
      randIn.param.alpha);
  ADD_PROCESS_PARAM("float","%f","--signal-ffinal",
      userParam.signalfFinal);
  ADD_PROCESS_PARAM("float","%f","--fl-signal",
      randIn.param.fLower);
  ADD_2PROCESS_PARAM("float","%f %f","--signal-mass-range",
      randIn.mMin, randIn.mMax);
  ADD_2PROCESS_PARAM("float","%f %f","--signal-psi0-range",
      randIn.psi0Min, randIn.psi0Max);
  ADD_2PROCESS_PARAM("float","%f %f","--signal-psi3-range",
      randIn.psi3Min, randIn.psi3Max);
  ADD_PROCESS_PARAM("int","%d","--seed",
      userParam.useed);
  ADD_PROCESS_PARAM("string","%s","--signal",	
      GetStringFromTemplate(userParam.signal));
  ADD_PROCESS_PARAM("int","%d","--signal-order",
      randIn.param.order);
  ADD_PROCESS_PARAM("string","%s","--template",	
      GetStringFromTemplate(userParam.template));
  ADD_PROCESS_PARAM("int","%d","--template-order",
      coarseBankIn.order);
  ADD_PROCESS_PARAM("float","%e","--tau0",	
      userParam.tau0);
  ADD_PROCESS_PARAM("float","%e","--tau3",	
      userParam.tau3);

  if (userParam.printResultXml){
    ADD_PROCESS_PARAM("float", "%s", "--print-result-xml"," ");
  }
  if (!userParam.startPhase){
    ADD_PROCESS_PARAM("float", "%s", "--no-start-phase"," ");
  }
  if (!userParam.alphaFConstraint ){
    ADD_PROCESS_PARAM("float", "%s", "--alpha-constraint"," ");
  } 
  else{
    ADD_PROCESS_PARAM("float", "%s", "--no-alpha-constraint"," ");
  }
  if (userParam.fastSimulation ){
    ADD_PROCESS_PARAM("float", "%s", "--fast-simulation"," ");
  }
  if (userParam.binaryInjection == BHNS){
    ADD_PROCESS_PARAM("float", "%s", "--bhns-injection", " ");
  }
  if (vrbflg){
    ADD_PROCESS_PARAM("float", "%s", "--verbose", " ");
  }

  	  
#undef ADD_PROCESS_PARAM
#undef ADD_2PROCESS_PARAM
}

CHAR* GetStringFromSimulationType(INT4 input)
{
  switch (input){
    case 0:
      return "SignalOnly";
      break;
    case 1:
      return "NoiseOnly";
      break;
    case 2:
      return "NoiseAndSignal";
      break;
  }
  return NULL;
}
CHAR* GetStringFromTemplate(INT4 input)
{
  switch (input){
    case EOB:
      return "EOB";
      break;
    case TaylorT1:
      return "TaylorT1";
      break;
    case TaylorT2:
      return "TaylorT2";
      break;
    case TaylorT3:
      return "TaylorT3";
      break;
    case PadeT1:
      return "PadeT1";
      break;
    case TaylorF2:
      return "TaylorF2";
      break;
    case BCV:
      return "BCV";
      break;
    case SpinTaylor:
      return "SpinTaylor";
      break;
  }
  return NULL;
}

CHAR* GetStringFromGridType(INT4 input)
 {
  switch (input){
  case Square:
    return "Square";
    break;
  case SquareNotOriented:
    return "SquareNotOriented";
    break;
  case HexagonalNotOriented:
    return "HexagonalNotOriented";
    break;
  case Hexagonal:
    return "Hexagonal";
    break;
  }
  return NULL;
}

CHAR* GetStringFromNoiseModel(INT4 input)
 {
  switch (input){
  case LIGOI:
    return "LIGOI";
    break;
  case LIGOA:
    return "LIGOA";
    break;
  case VIRGO:
    return "VIRGO";
    break;
  case GEO:
    return "GEO";
    break;
  case TAMA:    
    return "TAMA"; 
    break;
  case UNITY:    
    return "UNITY";
    break;
  case REALPSD:
    return "REALPSD";
    break;
  } 
  return NULL;
 }

CHAR* GetStringFromDetector(INT4 input)
 {
   
  switch (input){
  case L1: 
    return "L1";
    break;
  case H1: 
    return "H1";
    break;
  case H2: 
    return "H2";
    break;
  case V1: 
    return "V1";
    break;
  case G1: 
    return "G1";
    break;
  }
  return NULL;
 }


CHAR* GetStringFromScientificRun(INT4 input)
 {
  switch (input){
     case S2: 
       return "S2";
       break;
     case S3: 
       return "S3";
       break;
     case S1: 
       return "S1";
       break;
     case S4:
       return "S4";
       break;
     case S5:
       return "S5";
       break;
  }
  return NULL;
 }


/* xml file for the standalone code */
void 
BEPrintResultsXml( InspiralCoarseBankIn         coarseBankIn,
		   RandomInspiralSignalIn       randIn,
		   UserParametersIn             userParam,
		   ResultIn                     trigger
		  )
{
  LALStatus             status = blank_status;


  MetadataTable         templateBank;
  CHAR                  ifo[3];                           /* two character ifo code       */
  LIGOLwXMLStream       xmlStream;
  CHAR                  fname[256];
  LIGOTimeGPS           gpsStartTime 	= { 0, 0 };    /* input data GPS start time    */
  LIGOTimeGPS           gpsEndTime 	= { 0, 0 };      /* input data GPS end time      */
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;
  CHAR                  comment[LIGOMETA_COMMENT_MAX];
  CHAR                  ifoName[MAXIFO][LIGOMETA_IFO_MAX];

  MetadataTable         processParamsTable;
  ProcessParamsTable   *this_proc_param = NULL;

  snprintf( fname, sizeof(fname), 
	       BANKEFFICIENCY_PRINTRESULT_FILEXML ,
	       userParam.detector,
	       gpsStartTime.gpsSeconds,
	       gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );

  if (trigger.ntrial == 1){
    strncpy( ifoName[0], "no", LIGOMETA_IFO_MAX * sizeof(CHAR) );
    strncpy( ifoName[1], "ne", LIGOMETA_IFO_MAX * sizeof(CHAR) );
    memset( ifo, 0, sizeof(ifo) );
    memcpy( ifo, "MC", sizeof(ifo) - 1 );
    
    /* -- we start to fill the xml file here --- */
    memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
    LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlStream, fname), &status );
    
    /* create the process and process params tables */
    templateBank.processTable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
    LAL_CALL( LALGPSTimeNow ( &status, &(templateBank.processTable->start_time),
			      &accuracy ), &status );
    LAL_CALL( populate_process_table( &status, 
				      templateBank.processTable, 
				      PROGRAM_NAME, 
				      CVS_REVISION,
				      CVS_SOURCE,
				      CVS_DATE ), &status );
    this_proc_param = processParamsTable.processParamsTable = 
      (ProcessParamsTable *) calloc( 1, sizeof(ProcessParamsTable) );
    
    BEFillProc(this_proc_param, coarseBankIn, randIn, userParam);
    
    memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );
    
    /* write process table */
    snprintf( templateBank.processTable->ifos, LIGOMETA_IFOS_MAX, "%s%s", 
		 ifoName[0], ifoName[1] );
    LAL_CALL( LALGPSTimeNow ( &status, &(templateBank.processTable->end_time),
			      &accuracy ), &status );
    
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, process_table ), 
	      &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, templateBank, 
				      process_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );
    
    /* write process params table */
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, 
				      process_params_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, processParamsTable, 
				      process_params_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );
    
    /* finally write sngl inspiral table */
    
    PRINT_LIGOLW_XML_BANKEFFICIENCY(xmlStream.fp);
    
    fprintf(xmlStream.fp,BANKEFFICIENCY_PARAMS_ROW,
	    trigger.psi0_trigger,
	    trigger.psi3_trigger,
	    trigger.beta_trigger,
	    randIn.param.psi0, 
	    randIn.param.psi3, 
	    randIn.param.beta, 
	    trigger.tau0_trigger,
	    trigger.tau3_trigger,
	    randIn.param.t0, 
	    randIn.param.t3, 
	    trigger.fend_trigger, 
	    trigger.fend_inject,
	    trigger.mass1_inject,
	    trigger.mass2_inject,
	    randIn.param.startPhase,
	    trigger.rho_final,
	    trigger.snrAtCoaTime,
	    trigger.phase,
	    trigger.alphaF, 
	    trigger.bin,
	    randIn.param.nStartPad,
	    trigger.nfast
	    );
	    
     if (trigger.ntrial == (UINT4)userParam.ntrials){
       PRINT_LIGOLW_XML_TABLE_FOOTER(xmlStream.fp);
       PRINT_LIGOLW_XML_FOOTER(xmlStream.fp);
     }
      else
	{
	  fprintf(xmlStream.fp, ",\n");
	}
     fclose( xmlStream.fp );
     xmlStream.fp = NULL;
  }
  else 
    {
      xmlStream.fp = fopen(fname,"a+");
      
      fprintf(xmlStream.fp,BANKEFFICIENCY_PARAMS_ROW,
	      trigger.psi0_trigger,
	      trigger.psi3_trigger,
	      trigger.beta_trigger,
	      randIn.param.psi0, 
	      randIn.param.psi3, 
	      randIn.param.beta, 
	      trigger.tau0_trigger,
	      trigger.tau3_trigger,
	      randIn.param.t0, 
	      randIn.param.t3, 
	      trigger.fend_trigger, 
	      trigger.fend_inject,
	      trigger.mass1_inject,
	      trigger.mass2_inject,
	      randIn.param.startPhase,
	      trigger.rho_final,
	      trigger.snrAtCoaTime, 
	      trigger.phase,
	      trigger.alphaF, 
	      trigger.bin,
	      randIn.param.nStartPad,
	      trigger.nfast);	

      if (trigger.ntrial == (UINT4)userParam.ntrials){
	PRINT_LIGOLW_XML_TABLE_FOOTER(xmlStream.fp );
	PRINT_LIGOLW_XML_FOOTER(xmlStream.fp );
      }
      else
	{
	  fprintf(xmlStream.fp, ",\n");
	}
      
      
      fclose( xmlStream.fp );
      xmlStream.fp = NULL;
    }  

  /* close the output xml file */

}


/* print prototype for condor code (option --print-result-prototype)*/
void 
BEPrintProtoXml(InspiralCoarseBankIn   coarseBankIn,
		RandomInspiralSignalIn randIn,
		UserParametersIn           userParam
		)
{
  LALStatus             status = blank_status;


  MetadataTable         templateBank;
  CHAR                  ifo[3];                           /* two character ifo code       */
  LIGOLwXMLStream       xmlStream;
  CHAR                  fname[256];
  LIGOTimeGPS           gpsStartTime 	= { 0, 0 };    /* input data GPS start time    */
  LIGOTimeGPS           gpsEndTime 	= { 0, 0 };      /* input data GPS end time      */
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;
  CHAR                  comment[LIGOMETA_COMMENT_MAX];
  CHAR                  ifoName[MAXIFO][LIGOMETA_IFO_MAX];

  MetadataTable         processParamsTable;
  ProcessParamsTable   *this_proc_param = NULL;

  snprintf( fname, sizeof(fname), BANKEFFICIENCY_PRINTPROTO_FILEXML ,
	       ifo, gpsStartTime.gpsSeconds,
	       gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );




    strncpy( ifoName[0], "no", LIGOMETA_IFO_MAX * sizeof(CHAR) );
    strncpy( ifoName[1], "ne", LIGOMETA_IFO_MAX * sizeof(CHAR) );
    memset( ifo, 0, sizeof(ifo) );
    memcpy( ifo, "MC", sizeof(ifo) - 1 );
    
    
    
    /* -- we start to fill the xml file here --- */
    memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
    LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlStream, fname), &status );
    
    
    /* create the process and process params tables */
    templateBank.processTable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
    LAL_CALL( LALGPSTimeNow ( &status, &(templateBank.processTable->start_time),
			      &accuracy ), &status );
    LAL_CALL( populate_process_table( &status, templateBank.processTable, 
				      PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE ), &status );
    this_proc_param = processParamsTable.processParamsTable = 
      (ProcessParamsTable *) calloc( 1, sizeof(ProcessParamsTable) );
    
    BEFillProc(this_proc_param, coarseBankIn, randIn, userParam);
    
    memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );
    
    
    
    /* write process table */
    snprintf( templateBank.processTable->ifos, LIGOMETA_IFOS_MAX, "%s%s", 
		 ifoName[0], ifoName[1] );
    LAL_CALL( LALGPSTimeNow ( &status, &(templateBank.processTable->end_time),
			      &accuracy ), &status );
    
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, process_table ), 
	      &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, templateBank, 
				      process_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );
    
    /* write process params table */
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, 
				      process_params_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, processParamsTable, 
				      process_params_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );
    
    /* finally write sngl inspiral table */
    
     
      fclose( xmlStream.fp );
      xmlStream.fp = NULL;
#undef MAXIFO 


}



void
BEInitOverlapOutputIn(OverlapOutputIn *this){
  this->rhoMax   = -1.;
  this->phase    = -1.;
  this->rhoBin   = 0;
  this->templateNumber = 0;
  this->alpha    = -1.;
  this->layer    = 0;
  this->freq     = -1.;
}



/* Estimate the size of the longest template */
void BEGetMaximumSize(LALStatus  *status, 		      
		      RandomInspiralSignalIn  randIn,
		      InspiralCoarseBankIn coarseBankIn, 
		      UINT4 *length)
{
  /* Use these for the longest template from the bank */
  InspiralTemplate params;
  UINT4 maxTmpltLength = 0;

  INITSTATUS( status, "BEGetMaximumSize", BANKEFFICIENCYC );
  ATTATCHSTATUSPTR( status );
  
  randIn.param.massChoice 	= m1Andm2;
  randIn.param.approximant 	= EOB;
  *length = 0;
 
  LAL_CALL(LALInspiralWaveLength(status->statusPtr, length, randIn.param), 
	   status->statusPtr);
  
  params = randIn.param;
  params.mass1 = params.mass2 = coarseBankIn.mMin;

  LAL_CALL(LALInspiralWaveLength(status->statusPtr, &maxTmpltLength, params),
           status->statusPtr);

  /* Now return the longest of max template or signal */
  if (maxTmpltLength > *length)
    *length = maxTmpltLength;

  DETATCHSTATUSPTR( status );
  RETURN( status );  

}



void BECreatePsd(LALStatus                *status, 
		   InspiralCoarseBankIn   *coarseBankIn, 
		   RandomInspiralSignalIn *randIn,
		   UserParametersIn           userParam)
{
  UINT4 i=0;
  REAL4 df;
  REAL4 dummy; 
  FILE  *Foutput;
  FILE  *Finput;

  INITSTATUS( status, "BECreatePsd", BANKEFFICIENCYC );
  ATTATCHSTATUSPTR( status );


  memset( &(coarseBankIn->shf), 0, sizeof(REAL8FrequencySeries) );
  coarseBankIn->shf.f0 	= 0;
  LAL_CALL(LALDCreateVector( status->statusPtr, &(coarseBankIn->shf.data), randIn->psd.length ),
	   status->statusPtr);
  
  coarseBankIn->shf.deltaF 	= randIn->param.tSampling / ((randIn->psd.length-1)*2);
  
  /* --- Compute Noise Spectral Density --- */
  df = randIn->param.tSampling/(float) ((randIn->psd.length-1)*2);	 

  switch (userParam.noiseModel)						
    {
    case UNITY:
      LAL_CALL(LALNoiseSpectralDensity (status->statusPtr, coarseBankIn->shf.data, &LALLIGOIPsd, df),
	       status->statusPtr);
      for (i = 0; i < coarseBankIn->shf.data->length; i++)
	coarseBankIn->shf.data->data[i] = 1;
      break;
    case LIGOI:  
      LAL_CALL(LALNoiseSpectralDensity (status->statusPtr, coarseBankIn->shf.data, &LALLIGOIPsd, df), 
	       status->statusPtr);
      break;
    case LIGOA:  
      LAL_CALL(LALNoiseSpectralDensity (status->statusPtr, coarseBankIn->shf.data, &LALAdvLIGOPsd, df),
	       status->statusPtr);
      break;
    case VIRGO: 
      LAL_CALL(LALNoiseSpectralDensity (status->statusPtr, coarseBankIn->shf.data, &LALVIRGOPsd, df),
	       status->statusPtr);;
      break;
    case GEO:
      LAL_CALL(LALNoiseSpectralDensity (status->statusPtr, coarseBankIn->shf.data, &LALGEOPsd, df), 
	       status->statusPtr);
      break;
    case TAMA:
      LAL_CALL(LALNoiseSpectralDensity (status->statusPtr, coarseBankIn->shf.data, &LALTAMAPsd, df), 
	       status->statusPtr);
      break;
    case REALPSD:	
      /* Read psd dans le fichier InputPSD.dat */
      LAL_CALL(LALCreateRealPsd (status->statusPtr, coarseBankIn, *randIn , userParam), 
	       status->statusPtr);
      break;
    case READPSD:
     coarseBankIn->shf.data->data[0]=0;
     if ((Finput= fopen(userParam.inputPSD,"r")) == NULL)
       BEPrintError("userParam.inputPSD does not exists\n");
     else
     while(fscanf(Finput, "%f %le\n",&dummy, &coarseBankIn->shf.data->data[i+1]) == 2)  
         i++;

     if ((i+1) != coarseBankIn->shf.data->length){

       fprintf(stderr, "ERROR::number of points read(%d)  in the file (%s) and size of the time vectors (using the sampling(%f) and number of seconds (%f)) do not match.Generate a new file of (%d) points  with an expected df of %f \n",  i,userParam.inputPSD,randIn->param.tSampling, (coarseBankIn->shf.data->length-1)*2/randIn->param.tSampling,  coarseBankIn->shf.data->length-1, df);
       exit(0);
     }

     
     fclose(Finput);
      
     }

  /* --- save the psd in the file psd.dat if requested --- */
  if (userParam.printPsd)
   {
     Foutput= fopen(BANKEFFICIENCY_PRINTPSD_FILE,"w");
     for (i = 1; i < coarseBankIn->shf.data->length; i++)
       fprintf(Foutput, "%f %e\n",(float)i*df, coarseBankIn->shf.data->data[i]);  
     fclose(Foutput);
   }


  /* copy the psd in RandIn.psd */
  for (i = 0; i< coarseBankIn->shf.data->length; i++){
    randIn->psd.data[i] = coarseBankIn->shf.data->data[i];
  }


  

  DETATCHSTATUSPTR( status );
  RETURN( status );  
}



void BEGenerateInputData(LALStatus *status,
			 REAL4Vector * signal,
			 RandomInspiralSignalIn  *randIn,
			 UserParametersIn           userParam)
{
  UINT4 trial ;
  UINT4 success ;
  REAL4 u;

  INITSTATUS( status, "BEGenerateInputData", BANKEFFICIENCYC );
  ATTATCHSTATUSPTR( status );
  
  trial =0 ;
  success =0 ;
  
  /* we might force to compute a non random waveform giving the two
     input parameters m1 and m2 */      
  randIn->param.approximant = userParam.signal; 
  if (userParam.startPhase == 1){
    LALUniformDeviate( status->statusPtr, &u, randParams);
    CHECKSTATUSPTR(status);
    randIn->param.startPhase = u * LAL_PI;
  }
  else{
    randIn->param.startPhase = 0;
  }


  
  if (randIn->type != 1){
    if (userParam.signal == BCV ||userParam.signal == BCVSpin )
      {
	/* user parameters*/
	if (userParam.psi0!=-1 && userParam.psi3!=-1)
	  {
	    randIn->param.massChoice = fixedPsi;
	    fprintf(stderr, "signal = %f %f %f %f\n",randIn->param.beta, randIn->param.psi0,randIn->param.psi3, randIn->param.alpha1 );
	    LALRandomInspiralSignal(status->statusPtr, signal, randIn);
	    CHECKSTATUSPTR(status);
	  }
	else /* random parameters*/
	  {
	    randIn->param.massChoice = fixedPsi;
	    {
	      int valid = 0 ;
	      
	      while (!valid)
		{
		  REAL8 epsilon1 = (float) rand()/(float)RAND_MAX;
		  REAL8 epsilon2 = (float) rand()/(float)RAND_MAX;
		  REAL8 fend = 10000000.;
		  int trial = 0; 
		  
		  randIn->param.psi0 = randIn->psi0Min + epsilon1 * (randIn->psi0Max - randIn->psi0Min);
		  randIn->param.psi3 = randIn->psi3Min + epsilon2 * (randIn->psi3Max - randIn->psi3Min);
		  
		  
		  while (((fend > randIn->param.tSampling/2.) || ( fend < randIn->param.fLower)) && (trial < 10))
		    { 	     
		      
		      REAL8 fLR, fLSO;
		      REAL8 epsilon3 ;	
		      epsilon3 = (float) rand()/(float)RAND_MAX;
		      
		      randIn->param.totalMass = -randIn->param.psi3/(16.L*LAL_PI * LAL_PI * randIn->param.psi0);
		      randIn->param.totalMass *= 2 / LAL_MTSUN_SI;
		      
		      
		      fLR = 1.L/(LAL_PI * pow (3.L,1.5) * randIn->param.totalMass * LAL_MTSUN_SI); 
		      fLSO = 1.L/(LAL_PI * pow (6.L,1.5) * randIn->param.totalMass * LAL_MTSUN_SI);		
		      
		      fend = fLSO + (fLR - fLSO) * epsilon3;
		      fend = fLSO + (fLR - fLSO) * epsilon3;
		      
		      randIn->param.fFinal = fend;		  
		      randIn->param.fCutoff = fend;		  
		      
		      trial++;
		    }
		  if (trial == 10 ) valid = 0;
		  else valid = 1;
		}
	    }
	    
	    LALRandomInspiralSignal(status->statusPtr, signal, randIn);
	    CHECKSTATUSPTR(status);	    	    	    
	}
      }
    else /* EOB , T1 and so on*/
      {
	if (userParam.m1!=-1 && userParam.m2!=-1)
	  {
	    randIn->param.massChoice = fixedMasses;
	    LALRandomInspiralSignal(status->statusPtr, signal, randIn);
	    CHECKSTATUSPTR(status);	  
	  }
	else if (userParam.tau0!=-1 && userParam.tau3!=-1) 
	  {
	    randIn->param.massChoice = fixedTau;
	    LALRandomInspiralSignal(status->statusPtr, signal, randIn);
	    CHECKSTATUSPTR(status);
	  }
	else if (userParam.binaryInjection == BHNS)
	  {
	    randIn->param.massChoice = bhns;
	    LALRandomInspiralSignal(status->statusPtr, signal, randIn);
	    CHECKSTATUSPTR(status);
	  }
	else if (randIn->t0Min != 0)
          {
            randIn->param.massChoice = t03;
            LALRandomInspiralSignal(status->statusPtr, signal, randIn);
            CHECKSTATUSPTR(status);
          }
        else
	  {
	    randIn->param.massChoice = m1Andm2;
	    randIn->param.massChoice = totalMassUAndEta;	    
	    LALRandomInspiralSignal(status->statusPtr, signal, randIn);
	    CHECKSTATUSPTR(status);
	    
	}    
    }
  }


    if (randIn->type == 1)
      {
	randIn->param.massChoice = m1Andm2;	
	LALRandomInspiralSignal(status->statusPtr, signal, randIn);
	CHECKSTATUSPTR(status);
      }    
    

     
  DETATCHSTATUSPTR(status);
  RETURN (status);

}

void 
LALCreateRealPsd(LALStatus *status, 
		 InspiralCoarseBankIn *bankIn,
		 RandomInspiralSignalIn randIn, 
		 UserParametersIn userParam)
{




  enum
    {
      undefined,
      real_4,
      real_8
    } calData = undefined;

  enum { unset, urandom, user } randSeedType = user;    /* sim seed type */
  INT4  randomSeed        = 1;            /* value of sim rand seed       */
  REAL4 gaussVar          = 1.0;         /* variance of gaussian noise   */
  INT4  gaussianNoise     = 0;            /* make input data gaussian   0 non , 1 oui  */
  RandomParams *randParams = NULL;


  /* frame input data */
  FrCache      *frInCache = NULL;
  FrCache      *calCache = NULL;
  FrStream     *frStream = NULL;
  FrChanIn      frChan;
  CHAR tmpChName[LALNameLength];

  COMPLEX8FrequencySeries       injResp;
  COMPLEX8FrequencySeries      *injRespPtr;
  LIGOTimeGPS slideData   = {0,0};        /* slide data for time shifting */

  INT8  durationNS = 0;  
  INT8  gpsEndTimeNS     = 0;             /* input data GPS end time ns   */
  INT8  gpsStartTimeNS   = 0;             /* input data GPS start time ns */

  /* parameters used to generate calibrated power spectrum */
  LIGOTimeGPS gpsStartTime = { 0, 0 };    /* input data GPS start time    */
  LIGOTimeGPS gpsEndTime = { 0, 0 };      /* input data GPS end time      */
  INT4  padData = 8;                      /* saftety margin on input data */
  CHAR  *fqChanName       = NULL;


  /*CHAR *injectionFile = "/home/cokelaer/Work/TestWaveOverlap/HL-INJECTIONS_1-732005208-2048.xml";*/
  CHAR *injectionFile = NULL;
  CHAR  *frInCacheName = NULL;

  INT4  numPoints         = (randIn.psd.length-1)*2;           /* points in a segment          */
  INT4  numSegments       = 15;           /* number of segments           */
  CHAR  ifo[3];                           /* two character ifo code       */

  INT4  inputDataLength =  -1;              /* number of points in input    */
  INT4   resampFiltType   = 0;           /* low pass filter used for res 0 = = ldas*/
  INT4   sampleRate       = randIn.param.tSampling;           /* sample rate of filter data   */
  INT4   highPass         = 1;           /* enable high pass on raw data */
  REAL4  highPassFreq     = bankIn->fLower;            /* high pass frequency          */
  INT4   highPassOrder    = 8;           /* order of the td iir filter   */
  REAL4  highPassAtten    = 0.1;           /* attenuation of the td filter */
  REAL4  fLow             = bankIn->fLower;           /* low frequency cutoff         */
  INT4   specType         = 1;           /* use median (1)or mean(0) psd   or gauss    */

  CHAR  *calCacheName     = NULL;
  INT4   pointCal         = 0;            /* don't average cal over chunk */
  REAL4  dynRangeExponent = 40;            /*onent of dynamic range    */
  REAL4 geoHighPassFreq = -1;             /* GEO high pass frequency      */
  INT4  geoHighPassOrder = -1;            /* GEO high pass filter order   */
  REAL4 geoHighPassAtten = -1;            /* GEO high pass attenuation    */  
  /* raw iut data storage */
  REAL4TimeSeries               chan;
  REAL8TimeSeries               geoChan;
  REAL4FrequencySeries          spec;
  
  /* structures for preconditioning */
  ResampleTSParams              resampleParams;
    
  
  /* counters and other variables */
  UINT4 cut,  j, k;
  const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};
  
  LALUnitPair pair;
  REAL8 respRe, respIm;
  REAL8 shf;
  REAL8 inputLengthNS;
  UINT4 numInputPoints;
  const REAL8 epsilon = 1.0e-8;
  UINT4 resampleChan = 0;
  CalibrationUpdateParams calfacts, inj_calfacts;
  REAL8 dynRange = 0;
  REAL8 inputDeltaT;


  WindowSpectrumIn windowSpectrumParam;
  /*InspiralPipelineIn inspiralPipelineIn;*/

  CHAR         *calGlobPattern;

  INITSTATUS( status, "LALCreatRealPsd", BANKEFFICIENCYC );
  ATTATCHSTATUSPTR( status );

  
  /*  SetInspiralPipelineParams(&inspiralPipelineIn, randIn);*/


  /*  fprintf(stderr,"Generating real PSD\n");*/
  calData = undefined;
  /*calData = real_4;*/

  memset( ifo, 0, sizeof(ifo));

  switch (userParam.detector){
  case L1: 
    memcpy( ifo, "L1", sizeof(ifo) - 1 );
    break;
  case H1: 
    memcpy( ifo, "H1", sizeof(ifo) - 1 );
    break;
  case H2: 
    memcpy( ifo, "H2", sizeof(ifo) - 1 );
    break;
  case V1: 
  case G1: 
    break;
  }


  if ( dynRangeExponent )
  {
    /* compute the dynamic range scaling for the psd computation */
    dynRange = (REAL8) pow( 2.0, dynRangeExponent );
  }
  else
  {
    dynRange = 1.0;
  }
  if ( vrbflg )
    fprintf( stdout, "using dynamic range scaling %e\n", dynRange );




  /** Add on Thomas BEGIN */ 
  
  

  /* S3 */
  /*gpsStartTimeNS += (INT8) (752578908) * 1000000000LL;*/

  /* S2 */
   /*   gpsStartTimeNS += (INT8) (732005200 +padData) * 1000000000LL ;
  gpsStartTimeNS += (INT8) (0);
  */

  gpsStartTimeNS += (INT8) (userParam.startTime) *  1000000000LL;
  fqChanName    = userParam.chanName; 
  calCacheName  = userParam.calCacheName;
  frInCacheName = userParam.frInCacheName;

  inputDataLength = numPoints   * ( numSegments + 1) / 2 ;


  
  
  gpsEndTimeNS = gpsStartTimeNS;
  gpsEndTimeNS += (INT8) ((UINT8)inputDataLength /(UINT8)sampleRate) 
    * 1000000000LL;
  gpsEndTimeNS += (INT8) (0);
 
  LAL_CALL( LALINT8toGPS( status->statusPtr, &(gpsStartTime), &(gpsStartTimeNS) ), 
	    status->statusPtr ); 
  LAL_CALL( LALINT8toGPS( status->statusPtr, &(gpsEndTime), &(gpsEndTimeNS) ), 
	    status->statusPtr );
  /** End Thomas */

  /*
   *
   * read in the input data channel
   *
   */
 

  /* set the time series parameters of the input data and resample params */
  memset( &resampleParams, 0, sizeof(ResampleTSParams) );
  resampleParams.deltaT = 1.0 / (REAL8) sampleRate;
  
  /* set the params of the input data time series */
  memset( &chan, 0, sizeof(REAL4TimeSeries) );
  memset( &geoChan, 0, sizeof(REAL8TimeSeries) );
  chan.epoch = gpsStartTime;
  chan.epoch.gpsSeconds -= padData; /* subtract pad seconds from start */
  /* subtract slide from start */
  chan.epoch.gpsSeconds -= slideData.gpsSeconds; 
  chan.epoch.gpsNanoSeconds -= slideData.gpsNanoSeconds;
  /* copy the start time into the GEO time series */
  geoChan.epoch = chan.epoch;


  
  if ( vrbflg ) 
    fprintf( stdout,  "reading frame file locations from cache file: %s\n", frInCacheName );

  /* read a frame cache from the specified file */
  LAL_CALL( LALFrCacheImport( status->statusPtr, &frInCache, 
			      frInCacheName), status->statusPtr );
  

  /* open the input data frame stream from the frame cache */
  LAL_CALL( LALFrCacheOpen( status->statusPtr, &frStream,frInCache ), status->statusPtr );

  /* set the mode of the frame stream to fail on gaps or time errors */
  frStream->mode = LAL_FR_VERBOSE_MODE;

  /* seek to required epoch and set chan name */
  LAL_CALL( LALFrSeek( status->statusPtr, &(chan.epoch), frStream ), status->statusPtr );
  frChan.name = fqChanName;

  if ( calData == real_8 )
  {
    /* determine the sample rate of the raw data */
    LAL_CALL( LALFrGetREAL8TimeSeries( status->statusPtr, &geoChan, &frChan, frStream ),
        status->statusPtr );

    /* copy the data paramaters from the GEO channel to input data channel */
    snprintf( chan.name, LALNameLength * sizeof(CHAR), "%s", geoChan.name );
    chan.epoch          = geoChan.epoch;
    chan.deltaT         = geoChan.deltaT;
    chan.f0             = geoChan.f0;
    chan.sampleUnits    = geoChan.sampleUnits;
  }
  else
  {
    /* determine the sample rate of the raw data and allocate enough memory */
    LAL_CALL( LALFrGetREAL4TimeSeries( status->statusPtr, &chan, &frChan, frStream ),
        status->statusPtr );
  }


  inputDeltaT = chan.deltaT;



  /* determine if we need to resample the channel */
  if ( vrbflg )
  {
    fprintf( stdout, "resampleParams.deltaT = %e\n", resampleParams.deltaT );
    fprintf( stdout, "chan.deltaT = %e\n", chan.deltaT );
  }
  if ( ! ( fabs( resampleParams.deltaT - chan.deltaT ) < epsilon ) )
  {
    resampleChan = 1;
    if ( vrbflg )
      fprintf( stdout, "input channel will be resampled\n" );

    if ( resampFiltType == 0 )
    {
      resampleParams.filterType = LDASfirLP;
    }
    else if ( resampFiltType == 1 )
    {
      resampleParams.filterType = defaultButterworth;
    }
  }

  /* determine the number of points to get and create storage for the data */
  inputLengthNS = 
    (REAL8) ( gpsEndTimeNS - gpsStartTimeNS + 2000000000LL * padData );
  numInputPoints = (UINT4) floor( inputLengthNS / (chan.deltaT * 1.0e9) + 0.5 );
  if ( calData == real_8 )
    {
    /* create storage for the GEO input data */
      LAL_CALL( LALDCreateVector( status->statusPtr, &(geoChan.data), numInputPoints ), 
		status->statusPtr );
    }
  LAL_CALL( LALSCreateVector( status->statusPtr, &(chan.data), numInputPoints ), 
	    status->statusPtr );

  if ( vrbflg ) fprintf( stdout, "input channel %s has sample interval "
      "(deltaT) = %e\nreading %d points from frame stream\n", fqChanName, 
      chan.deltaT, numInputPoints );

  if ( calData == real_8 )
  {
    /* read in the GEO data here */
    PassBandParamStruc geoHighpassParam;

    /* read the GEO data from the time series into geoChan      */
    /* which already has the correct amount of memory allocated */
    if ( vrbflg ) fprintf( stdout, "reading GEO data from frames... " );

    LAL_CALL( LALFrGetREAL8TimeSeries( status->statusPtr, &geoChan, &frChan, frStream ),
        status->statusPtr);

    if ( vrbflg ) fprintf( stdout, "done\n" );

    /* high pass the GEO data using the parameters specified on the cmd line */
    geoHighpassParam.nMax = geoHighPassOrder;
    geoHighpassParam.f1 = -1.0;
    geoHighpassParam.f2 = (REAL8) geoHighPassFreq;
    geoHighpassParam.a1 = -1.0;
    geoHighpassParam.a2 = (REAL8)(1.0 - geoHighPassAtten);
    if ( vrbflg ) fprintf( stdout, "applying %d order high pass to GEO data: "
        "%3.2f of signal passes at %4.2f Hz\n", 
			   geoHighpassParam.nMax, 
			   geoHighpassParam.a2,
			   geoHighpassParam.f2 );

    LAL_CALL( LALButterworthREAL8TimeSeries( status->statusPtr, &geoChan, 
          &geoHighpassParam ), status->statusPtr );

    /* cast the GEO data to REAL4 in the chan time series       */
    /* which already has the correct amount of memory allocated */
    for ( j = 0 ; j <numInputPoints ; ++j )
    {
      chan.data->data[j] = (REAL4) ( geoChan.data->data[j] * dynRange );
    }

    /* re-copy the data paramaters from the GEO channel to input data channel */
    snprintf( chan.name, LALNameLength * sizeof(CHAR), "%s", geoChan.name );
    chan.epoch          = geoChan.epoch;
    chan.deltaT         = geoChan.deltaT;
    chan.f0             = geoChan.f0;
    chan.sampleUnits    = geoChan.sampleUnits;

    /* free the REAL8 GEO input data */
    LAL_CALL( LALDDestroyVector( status->statusPtr, &(geoChan.data) ), status->statusPtr );
    geoChan.data = NULL;
  }
  else if ( calData == real_4)
    {
      /* read the data channel time series from frames */
      LAL_CALL( LALFrGetREAL4TimeSeries( status->statusPtr, &chan, &frChan, frStream ),
		status->statusPtr );
      
      /* multiply the input data by dynRange */
      for ( j = 0 ; j < numInputPoints ; ++j )
	{
	  chan.data->data[j] *= dynRange;           
	}
    }
  else
    {
      /* read the data channel time series from frames */
      LAL_CALL( LALFrGetREAL4TimeSeries( status->statusPtr, &chan, &frChan, frStream ),
		status->statusPtr );	
    }



  memcpy( &(chan.sampleUnits), &lalADCCountUnit, sizeof(LALUnit) );


  /* close the frame file stream and destroy the cache */
  LAL_CALL( LALFrClose( status->statusPtr, &frStream ), status->statusPtr );


  LAL_CALL( LALDestroyFrCache( status->statusPtr, &frInCache ), status->statusPtr );


  if ( vrbflg ) fprintf( stdout, "read channel %s from frame stream\n"
      "got %d points with deltaT %e\nstarting at GPS time %d sec %d ns\n", 
      chan.name, chan.data->length, chan.deltaT, 
      chan.epoch.gpsSeconds, chan.epoch.gpsNanoSeconds );



  /*
   *
   * create the random seed if it will be needed
   *
   */


  if ( randSeedType != unset )
    {
      
    if ( randSeedType == urandom )
    {
      FILE   *fpRand = NULL;
      INT4    randByte;

      if ( vrbflg ) 
        fprintf( stdout, "obtaining random seed from /dev/urandom: " );

      randomSeed = 0;
      fpRand = fopen( "/dev/urandom", "r" );
      if ( fpRand )
      {
        for ( randByte = 0; randByte < 4 ; ++randByte )
        {
          INT4 tmpSeed = (INT4) fgetc( fpRand );
          randomSeed += tmpSeed << ( randByte * 8 );
        }
        fclose( fpRand );
      }
      else
      {
        perror( "error obtaining random seed from /dev/urandom" );
        exit( 1 );
      }
    }
    else if ( randSeedType == user )
    {
      if ( vrbflg ) 
        fprintf( stdout, "using user specified random seed: " );
    }
    else
    {
      /* should never get here */
      fprintf( stderr, "error obtaining random seed\n" );
      exit( 1 );
    }

    if ( vrbflg ) fprintf( stdout, "%d\n", randomSeed );

    /* create the tmplt bank random parameter structure */
    LAL_CALL( LALCreateRandomParams( status->statusPtr, &randParams, randomSeed ),
        status->statusPtr );
  }

  /* replace the input data with gaussian noise if necessary */
  if ( gaussianNoise )
  {
    if ( vrbflg ) fprintf( stdout, 
        "setting input data to gaussian noise with variance %e... ", gaussVar );
    memset( chan.data->data, 0, chan.data->length * sizeof(REAL4) );
    LAL_CALL( LALNormalDeviates(status->statusPtr, chan.data, randParams ), 
	      status->statusPtr );
    for ( j = 0; j < chan.data->length; ++j )
    {
      chan.data->data[j] *= gaussVar;
    }
    if ( vrbflg ) fprintf( stdout, "done\n" );

  }


  

  /*
   *
   * generate the response function for the requested time
   *
   */
  
  memset( &resp, 0, sizeof(COMPLEX8FrequencySeries) );
  LAL_CALL( LALCCreateVector( status->statusPtr, &(resp.data), numPoints / 2 + 1 ), 
      status->statusPtr );
  
  /* set the parameters of the response to match the data */
  resp.epoch.gpsSeconds = chan.epoch.gpsSeconds + padData;
  resp.epoch.gpsNanoSeconds = chan.epoch.gpsNanoSeconds;
  resp.deltaF = (REAL8) sampleRate / (REAL8) numPoints;
  resp.sampleUnits = strainPerCount;
  strcpy( resp.name, chan.name );

  /* generate the response function for the current time */
  if ( vrbflg ) fprintf( stdout, "generating response at time %d sec %d ns\n",
			 resp.epoch.gpsSeconds, resp.epoch.gpsNanoSeconds );
  
  /* initialize the calfacts */
  memset( &calfacts, 0, sizeof(CalibrationUpdateParams) );
  calfacts.ifo = ifo;

  /* determine length of chunk */
  if ( pointCal )
    {
    calfacts.duration.gpsSeconds = 1;
    calfacts.duration.gpsNanoSeconds = 0;
    }
  else
    {
     durationNS = gpsEndTimeNS - gpsStartTimeNS;
      LAL_CALL( LALINT8toGPS( status->statusPtr, &(calfacts.duration), 
			      &durationNS ), status->statusPtr );
    }

  
  if ( calData )
  {
    /* if we are using calibrated data set the response to unity */

    for( k = 0; k < resp.data->length; ++k )
    {
      resp.data->data[k].re = (REAL4) (1.0 / dynRange);
      resp.data->data[k].im = 0.0;
    }
  }
  else
  {
    calGlobPattern = NULL;        
    if ( vrbflg ) fprintf( stdout, 
          "reading calibration data from cache: %s ....", calCacheName );
    LAL_CALL( LALCreateCalibFrCache( status->statusPtr, &calCache, calCacheName, 
          NULL, calGlobPattern ), status->statusPtr );
  
    
    if ( calGlobPattern ) LALFree( calGlobPattern );


    /* get the response from the frame data */
    LAL_CALL( LALExtractFrameResponse( status->statusPtr, &resp, calCache, 
				       &calfacts), status->statusPtr );
    
    LAL_CALL( LALDestroyFrCache( status->statusPtr, &calCache), status->statusPtr );
    if ( vrbflg ) fprintf( stdout, 
			   "for calibration of data, alpha = %f and alphabeta = %f\n",
			   (REAL4) calfacts.alpha.re, (REAL4) calfacts.alphabeta.re);
  }


  if ( gaussianNoise )
    {
      /* replace the response function with unity if */
      /* we are filtering gaussian noise             */
      if ( vrbflg ) fprintf( stdout, "setting response to unity... " );
      for ( k = 0; k < resp.data->length; ++k )
	{
	  resp.data->data[k].re = 1.0;
	  resp.data->data[k].im = 0;
	}
      if ( vrbflg ) fprintf( stdout, "done\n" );
      
  
    }
  
        injectionFile = NULL;
    


  if ( injectionFile )
  {
    /* get injections within 500 seconds of either end of the segment.   */
    /* a 0.4,0.4 MACHO starting at 30.0 Hz has length 435.374683 seconds */
    /* so this should be plenty of safety. better to waste cpu than miss */
    /* injected signals...                                               */
    INT4 injSafety = 500;
    int  numInjections = 0;
    SimInspiralTable    *injections = NULL;
    SimInspiralTable    *thisInj = NULL;

    /* read in the injection data from XML */
    numInjections = SimInspiralTableFromLIGOLw( &injections, injectionFile,
        gpsStartTime.gpsSeconds - injSafety, 
        gpsEndTime.gpsSeconds + injSafety );

    if ( numInjections < 0 )
    {
      fprintf( stderr, "error: cannot read injection file" );
      exit( 1 );
    }
    else if ( numInjections )
    {

      fprintf(stderr, "resample injection\n");

      /* see if we need a higher resolution response to do the injections */
      if ( resampleChan )	
      {


        /* we need a different resolution of response function for injections */
        UINT4 rateRatio = floor( resampleParams.deltaT / chan.deltaT + 0.5 );
        UINT4 rawNumPoints = rateRatio * numPoints;

        if ( vrbflg ) fprintf( stdout, "rateRatio = %d\nrawNumPoints = %d\n"
            "chan.deltaT = %e\n", rateRatio, rawNumPoints, chan.deltaT );

        memset( &injResp, 0, sizeof(COMPLEX8FrequencySeries) );
        LAL_CALL( LALCCreateVector( status->statusPtr, &(injResp.data), 
              rawNumPoints / 2 + 1 ), status->statusPtr );
        injResp.epoch = resp.epoch;
        injResp.deltaF = 1.0 / ( rawNumPoints * chan.deltaT );
        injResp.sampleUnits = strainPerCount;
        strcpy( injResp.name, chan.name );

	
        if ( calData )
        {
          /* if we are using calibrated data set the response to unity */
          if ( vrbflg ) fprintf( stdout, 
              "setting injection response to inverse dynRange... " );
	  dynRange = (REAL8) pow( 2.0, 40);
          for ( k = 0; k < injResp.data->length; ++k )
          {
            injResp.data->data[k].re = (REAL4)(1.0/dynRange);
            injResp.data->data[k].im = 0.0;
          }
          injRespPtr = &injResp;
        }
        else
        {
          /* generate the response function for the current time */
          if ( vrbflg ) fprintf( stdout, 
              "generating high resolution response at time %d sec %d ns\n"
              "length = %d points, deltaF = %e Hz\n",
              resp.epoch.gpsSeconds, resp.epoch.gpsNanoSeconds,
              injResp.data->length, injResp.deltaF );

          /* initialize the inj_calfacts */
          memset( &inj_calfacts, 0, sizeof(CalibrationUpdateParams) );
          inj_calfacts.ifo = ifo;

	  LAL_CALL( LALINT8toGPS( status->statusPtr, &(inj_calfacts.duration), 
				  &durationNS ), status->statusPtr );
	  
          /* create the lal calibration frame cache */
	  calGlobPattern = NULL;
	  if ( vrbflg ) fprintf( stdout, 
				 "reading calibration data from cache: %s\n", calCacheName );
	  
	  LAL_CALL( LALCreateCalibFrCache( status->statusPtr, &calCache, calCacheName, 
                NULL, calGlobPattern ), status->statusPtr );
          if ( calGlobPattern ) LALFree( calGlobPattern );


	  /* extract the calibration from frames */
          LAL_CALL( LALExtractFrameResponse( status->statusPtr, &injResp, calCache, 
                &inj_calfacts ), status->statusPtr );

	  
	  LAL_CALL( LALDestroyFrCache( status->statusPtr, &calCache), status->statusPtr );
	  

          injRespPtr = &injResp;	  
	}
	if ( gaussianNoise )
	  {
	    /* replace the response function with unity if */
	    /* we are filtering gaussian noise             */
	    if ( vrbflg ) fprintf( stdout, "setting response to unity... " );
	    for ( k = 0; k < injResp.data->length; ++k )
	      {
		injResp.data->data[k].re = 1.0;
		injResp.data->data[k].im = 0;
	      }
	    if ( vrbflg ) fprintf( stdout, "done\n" );
	    
        }

      }
      else
      {
        /* the data is already at the correct sample rate, just do injections */
        injRespPtr = &resp;
        memset( &injResp, 0, sizeof(COMPLEX8FrequencySeries) );
      }

      /* inject the signals, preserving the channel name (Tev mangles it) */


      /* if injectOverhead option, then set chan.name to "ZENITH".  
       * This causes no detector site to be found in the injection code so
       * that the injection is done directly overhead (i.e. with a response 
       * function of F+ = 1; Fx = 0) */
      snprintf( tmpChName, LALNameLength * sizeof(CHAR), "%s", chan.name );
      
      snprintf( chan.name, LALNameLength * sizeof(CHAR), "ZENITH" );


      LAL_CALL( LALFindChirpInjectSignals( status->statusPtr, &chan, injections, 
            injRespPtr ), status->statusPtr );
      snprintf( chan.name,  LALNameLength * sizeof(CHAR), "%s", tmpChName );

      if ( vrbflg ) fprintf( stdout, "injected %d signals from %s into %s\n", 
          numInjections, injectionFile, chan.name );

      while ( injections )
      {
        thisInj = injections;
        injections = injections->next;
        LALFree( thisInj );
      }

      /* write the raw channel data plus injections to the output frame file */

      if ( injResp.data )
        LAL_CALL( LALCDestroyVector( status->statusPtr, &(injResp.data) ), status->statusPtr );
    }
    else
    {
      if ( vrbflg ) fprintf( stdout, "no injections in this chunk\n" );
    }
  }


  /* resample the input data */
  if ( resampleChan )
  {
    if (vrbflg) fprintf( stdout, "resampling input data from %e to %e\n",
        chan.deltaT, resampleParams.deltaT );

    LAL_CALL( LALResampleREAL4TimeSeries( status->statusPtr, &chan, &resampleParams ),
        status->statusPtr );

    if ( vrbflg ) fprintf( stdout, "channel %s resampled:\n"
        "%d points with deltaT %e\nstarting at GPS time %d sec %d ns\n", 
        chan.name, chan.data->length, chan.deltaT, 
        chan.epoch.gpsSeconds, chan.epoch.gpsNanoSeconds );
  }


  /* 
   *
   * high pass the data, removed pad from time series and check length of data
   *
   */
  
  /* iir filter to remove low frequencies from data channel */
  if ( highPass )
  {
    PassBandParamStruc highpassParam;
    highpassParam.nMax = highPassOrder;
    highpassParam.f1 = -1.0;
    highpassParam.f2 = (REAL8) highPassFreq;
    highpassParam.a1 = -1.0;
    highpassParam.a2 = (REAL8)(1.0 - highPassAtten); /* a2 is not attenuation */

    if ( vrbflg ) fprintf( stdout, "applying %d order high pass: "
        "%3.2f of signal passes at %4.2f Hz\n", 
        highpassParam.nMax, highpassParam.a2, highpassParam.f2 );

    LAL_CALL( LALButterworthREAL4TimeSeries( status->statusPtr, &chan, &highpassParam ),
        status->statusPtr );
  }

  /* remove pad from requested data from start and end of time series */
  memmove( chan.data->data, chan.data->data + padData * sampleRate, 
      (chan.data->length - 2 * padData * sampleRate) * sizeof(REAL4) );
  LALRealloc( chan.data->data, 
      (chan.data->length - 2 * padData * sampleRate) * sizeof(REAL4) );
  chan.data->length -= 2 * padData * sampleRate;
  chan.epoch.gpsSeconds += padData;

  if ( vrbflg ) fprintf( stdout, "after removal of %d second padding at "
      "start and end:\ndata channel sample interval (deltaT) = %e\n"
      "data channel length = %d\nstarting at %d sec %d ns\n", 
      padData , chan.deltaT , chan.data->length, 
      chan.epoch.gpsSeconds, chan.epoch.gpsNanoSeconds );


  /*
   *
   * power spectrum estimation
   *
   */


  /* create storage for the power spectral estimate */
  memset( &spec, 0, sizeof(REAL4FrequencySeries) );
  LAL_CALL( LALSCreateVector( status->statusPtr, &(spec.data), numPoints / 2 + 1 ), 
      status->statusPtr );


  windowSpectrumParam.numPoints =  numPoints ;
  windowSpectrumParam.gaussVar =  gaussVar ;
  windowSpectrumParam.inputDeltaT =  inputDeltaT ;
  windowSpectrumParam.specType =  specType ;


  LAL_CALL( LALComputeWindowSpectrum( status->statusPtr, 
				      &windowSpectrumParam,
				      &spec, 
				      &chan ),
	    status->statusPtr );




  /* set low frequency cutoff of power spectrum */
  cut = fLow / spec.deltaF > 1 ?  fLow / spec.deltaF : 1;

  /* compute a calibrated strain power spectrum */
  
  bankIn->shf.epoch = spec.epoch;
  memcpy( bankIn->shf.name, spec.name, LALNameLength * sizeof(CHAR) );
  bankIn->shf.deltaF = spec.deltaF;
  bankIn->shf.f0 = spec.f0;
  bankIn->shf.data = NULL;
  
  pair.unitOne = &(spec.sampleUnits);
  pair.unitTwo = &(resp.sampleUnits);
  LAL_CALL( LALUnitMultiply( status->statusPtr, &(bankIn->shf.sampleUnits), &pair ), 
      status->statusPtr );
  LAL_CALL( LALDCreateVector( status->statusPtr, &(bankIn->shf.data), spec.data->length ),
      status->statusPtr );
  memset( bankIn->shf.data->data, 0, 
      bankIn->shf.data->length * sizeof(COMPLEX8) ); 
  
  shf = spec.data->data[cut] * 
    ( resp.data->data[cut].re * resp.data->data[cut].re +
      resp.data->data[cut].im * resp.data->data[cut].im );

  for ( k = 1; k < cut ; ++k )
  {
    bankIn->shf.data->data[k] = shf;
  }

  for ( k = cut; k < bankIn->shf.data->length; ++k )
  {
    respRe = (REAL8) resp.data->data[k].re;
    respIm = (REAL8) resp.data->data[k].im;
    bankIn->shf.data->data[k] = (REAL8) spec.data->data[k] *
      ( respRe * respRe + respIm * respIm );

    
  }
  

/*  Foutput= fopen("spec.dat","w");

  for (i = 1; i < (int)bankIn->shf.data->length; i++)
    fprintf(Foutput, "%15.10e\n", bankIn->shf.data->data[i]);  
  fclose(Foutput);*/

   LALDPrintFrequencySeries( &(bankIn->shf), "spec.dat" );  
  /*
   Foutput= fopen("uncalibratedpsd.dat","w");

  for (i = 1; i < spec.data->length; i++)
    fprintf(Foutput, "%15.10e\n",spec.data->data[i]);  
  fclose(Foutput); 


  Foutput= fopen("respre.dat","w");

  for (i = 0; i < spec.data->length; i++){
    respRe = (REAL8) resp.data->data[i].re;
    fprintf(Foutput, "%e\n",respRe);  
  }
  fclose(Foutput); 
  
  Foutput= fopen("respim.dat","w");

  for (i = 1; i < spec.data->length; i++){
    respIm = (REAL8) resp.data->data[i].im;
    fprintf(Foutput, "%e\n",respIm);  
  }
  fclose(Foutput); 

  */

  /* ADD on Thomas Create strain data to be process by Overlap or something else */

  {
    REAL4       *dataPtr;
    RealFFTPlan                  *fwdPlan = NULL;


    REAL4Vector    *rawSegment = NULL ; 


    INT4 length;

    length = numPoints/2 + 1;
    strainSegment = NULL;

    LALCreateVector( status->statusPtr, &rawSegment,
		     numPoints);
    CHECKSTATUSPTR( status );
    
    LALCCreateVector( status->statusPtr, &strainSegment,
		      numPoints/2+1);    
    CHECKSTATUSPTR( status );
 

     
    dataPtr = chan.data->data;
    rawSegment->data = dataPtr;
   
    
    LALCreateForwardRealFFTPlan( status->statusPtr, &fwdPlan, 
				 numPoints, 0);
    CHECKSTATUSPTR( status );

    /*     Foutput= fopen("raw0.dat","w");
    for (i = 1; i < rawSegment->length; i++)
      fprintf(Foutput, "%e\n", rawSegment->data[i]);  
    fclose(Foutput); 
    */
   
    
    
    LALForwardRealFFT( status->statusPtr, strainSegment,
		       rawSegment, fwdPlan );
    CHECKSTATUSPTR( status );



    dynRange = (REAL8) pow( 2.0, 0);
    
    
    for ( k = 0; k < strainSegment->length; ++k )
      {
	
	REAL4 p = strainSegment->data[k].re;
	REAL4 q = strainSegment->data[k].im;
	REAL4 x = resp.data->data[k].re * dynRange;
	REAL4 y = resp.data->data[k].im * dynRange;
	
	strainSegment->data[k].re =  p*x - q*y;
	strainSegment->data[k].im =  p*y + q*x;
      }
    
    LALDestroyRealFFTPlan(status->statusPtr,&fwdPlan);
       

    for ( k = 0; k < cut; ++k )
      {
	strainSegment->data[k].re = 0.0;
	strainSegment->data[k].im = 0.0;
      }
      
    dynRange = (REAL8) pow( 2.0, 40);/* pour avoir meme resultat que LAL avant filtrage*/
    /*
    Foutput= fopen("finalstrainre.dat","w");
    
    for (i = 0; i < strainSegment->length; i++)
      fprintf(Foutput, "%15.10e\n", strainSegment->data[i].re*dynRange );  
    fclose(Foutput); 
      
      Foutput= fopen("finalstrainim.dat","w");
      
      for (i = 0; i < strainSegment->length; i++)
	fprintf(Foutput, "%15.10e\n", strainSegment->data[i].im*dynRange );  
	
	fclose(Foutput); 
    */
    
  
  }
  /* END ADD on Thomas*/                                                                                                                        


  /* spec, uncalibrated spec, and raw0.dat identical to 
     inspiral in FindChipSPData. */

  /*if one want to inject BCV directly in the bcv injection*/
#if 0
 { 
   REAL4Vector            signal;
   RandomInspiralSignalIn randIn;
   UINT4 n,k;
  RealFFTPlan 			*fwdp = NULL;


   n 	= (strainSegment->length-1)*2;
   signal.length = n;
   signal.data 		= (REAL4*) LALCalloc(1, sizeof(REAL4) * signal.length);

   LALCreateForwardRealFFTPlan(status->statusPtr, &fwdp, signal.length, 0);
   randIn.mMin = 3; 
   randIn.mMax =5; 
   randIn.MMax = 10;

   randIn.param.tSampling  = sampleRate;
   randIn.param.massChoice = fixedPsi;
   randIn.param.psi0 = 200000;
   randIn.param.psi3 = -500;
   randIn.param.fCutoff = 900;
   randIn.param.alpha = .99;
   randIn.param.approximant = BCV;
   randIn.param.massChoice = fixedPsi;
   randIn.param.fLower  = 40 ; 
   randIn.psd.length 	= strainSegment->length;
   randIn.psd.data 	= (REAL8*) LALMalloc(sizeof(REAL8) * randIn.psd.length); 
   randIn.type = 0;
   randIn.param.order = twoPN;
   randIn.fwdp = fwdp;
   randIn.param.signalAmplitude = 1;
   randIn.param.startPhase = 0;
   randIn.param.nStartPad =  128 *sampleRate;


   fprintf(stderr, "%d %d %d\n",n,bankIn->shf.data->length, strainSegment->length);
   for ( k = 1; k < bankIn->shf.data->length; ++k )
     {    
       randIn.psd.data[k] =  bankIn->shf.data->data[k];
     }
      

      
   LALRandomInspiralSignal(status->statusPtr, &signal, &randIn);
   CHECKSTATUSPTR(status);
   
   
   for ( k = 0; k < n/2; ++k )
     {
       
       strainSegment->data[k].re += signal.data[k] *10;
       strainSegment->data[k].im += signal.data[n - k] *10;
     
       }
      
 }
#endif



  LALCheckMemoryLeaks();    




  DETATCHSTATUSPTR( status );
  RETURN( status );  
}



void LALComputeWindowSpectrum(LALStatus *status, 
			      WindowSpectrumIn *param,
			      REAL4FrequencySeries  *spec,
			      REAL4TimeSeries *chan)
{
  
  LALWindowParams               wpars;
  AverageSpectrumParams         avgSpecParams;
  INT4 k;
  
  INITSTATUS( status, "LALCreatRealPsd", BANKEFFICIENCYC );
  ATTATCHSTATUSPTR( status );
  
  /* compute the windowed power spectrum for the data channel */
  avgSpecParams.window = NULL;
  avgSpecParams.plan   = NULL; 
  
  LAL_CALL( LALCreateForwardRealFFTPlan( status->statusPtr, 
					 &(avgSpecParams.plan), 
					 param->numPoints, 0 ),
	    status->statusPtr );
  
  
  switch ( param->specType )
    {
    case 0:
      avgSpecParams.method = useMean;
      if ( vrbflg ) fprintf( stdout, "computing mean psd" );
      break;
    case 1:
      avgSpecParams.method = useMedian;
      if ( vrbflg ) fprintf( stdout, "computing median psd" );
      break;
    case 2:
      avgSpecParams.method = useUnity;
      if ( vrbflg ) fprintf( stdout, "simulation gaussian noise psd" );
      break;
    }
  
  
  wpars.type = Hann;
  wpars.length = param->numPoints;
  avgSpecParams.overlap = param->numPoints / 2;
  if ( vrbflg ) 
    fprintf( stdout, " with overlap %d\n", avgSpecParams.overlap );
  
  LAL_CALL( LALCreateREAL4Window( status->statusPtr, &(avgSpecParams.window),
				  &wpars ), status->statusPtr );
  LAL_CALL( LALREAL4AverageSpectrum( status->statusPtr, spec, chan, &avgSpecParams ),
	    status->statusPtr );
  XLALDestroyREAL4Window( avgSpecParams.window );
  LAL_CALL( LALDestroyRealFFTPlan( status->statusPtr, &(avgSpecParams.plan) ), status->statusPtr );

  strcpy( spec->name, chan->name );

  if ( param->specType == 2 )
  {
    /* multiply the unit power spectrum to get a gaussian psd */
    REAL4 gaussVarSq = param->gaussVar * param->gaussVar;
    if ( param->inputDeltaT != chan->deltaT )
    {
      /* reduce the variance as we have resampled the data */
      gaussVarSq *= param->inputDeltaT / chan->deltaT;
    }
    for ( k = 0; k < (int)spec->data->length; ++k )
    {
      spec->data->data[k] *= 2.0 * gaussVarSq * (REAL4) chan->deltaT;
    }

    if ( vrbflg ) 
      fprintf( stdout, "set psd to constant value = %e\n", spec->data->data[0] );
  }


  DETATCHSTATUSPTR( status );
  RETURN( status );  

}


void SetInspiralPipelineParam(InspiralPipelineIn *param,
			      RandomInspiralSignalIn randIn)
{

  /* parameters used to generate calibrated power spectrum */
  /*  param->gpsStartTime = {0,0};
    param->gpsEndTime = 0;
  */
  param->padData = 8;                
  param->fqChanName       = NULL;


  param->injectionFile = "/home/cokelaer/Work/TestWaveOverlap/HL-INJECTIONS_1-732005208-2048.xml";
  param->frInCacheName = NULL;

  param->numPoints         = (randIn.psd.length-1)*2;           /* points in a segment          */
  param->numSegments       = 15;           /* number of segments           */
  

  param->inputDataLength =  -1;            
  param->resampFiltType   = 0;           
  param->sampleRate       = randIn.param.tSampling;    
  param->highPass         = 1;           
  param->highPassFreq     = randIn.param.fLower;            
  param->highPassOrder    = 8;           
  param->highPassAtten    = 0.1;         
  param->fLow             = randIn.param.fLower;
  param->specType         = 1;           

  param->calCacheName     = NULL;
  param->pointCal         = 0;            
  param->dynRangeExponent = 40;            
  param->geoHighPassFreq = -1;     
  param->geoHighPassOrder = -1;    
  param->geoHighPassAtten = -1;    
}






void BECreatePowerVector(LALStatus              *status, 
			 BEPowerVector          *powerVector,
			 RandomInspiralSignalIn  randIn, 
			 INT4                    length)
{


  INITSTATUS (status, "BECreatePowerVector", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);

  powerVector->fm2_3.length = length / 2 ;
  powerVector->fm5_3.length = length / 2;       
  powerVector->fm7_6.length = length / 2 ;
  powerVector->fm1_2.length = length / 2;       

  powerVector->fm2_3.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) * powerVector->fm2_3.length);
  powerVector->fm5_3.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) * powerVector->fm5_3.length);
  powerVector->fm7_6.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) * powerVector->fm7_6.length);
  powerVector->fm1_2.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) * powerVector->fm1_2.length);
  
  /* --- Create useful vector once for all --- */
  LAL_CALL(LALCreateVectorFreqPower(&powerVector->fm5_3, randIn.param, -5, 3), 
	   status->statusPtr);
  LAL_CALL(LALCreateVectorFreqPower(&powerVector->fm2_3, randIn.param, -2, 3), 
	   status->statusPtr);
  LAL_CALL(LALCreateVectorFreqPower(&powerVector->fm7_6, randIn.param, -7, 6), 
	   status->statusPtr);
  LAL_CALL(LALCreateVectorFreqPower(&powerVector->fm1_2, randIn.param, -1, 2), 
	   status->statusPtr);

  DETATCHSTATUSPTR( status );
  RETURN( status );  
}



void LALInspiralOverlapBCV(LALStatus *status,
			   InspiralTemplate   *list,
			   BEPowerVector          *powerVector,
			   UserParametersIn           *userParam, 
			   RandomInspiralSignalIn *randIn,
			   REAL4Vector            *FilterBCV1,
			   REAL4Vector            *FilterBCV2,
			   InspiralWaveOverlapIn  *overlapin,
			   OverlapOutputIn        *overlapout,
			   REAL4Vector            *correlation,
			   BEMoments              *moments)
			   

{
  REAL4 fendBCV ;
  REAL4 df;
  INT4 n, kMin, kMax; 

  INITSTATUS (status, "LALInspiralOverlapBCV", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);

  n    = FilterBCV1->length;
  df = randIn->param.tSampling / (REAL8)(n/2) / 2.;      

  kMin  = floor(randIn->param.fLower /df );
  fendBCV = list->fFinal;


  /* if we want to test the BCV metric by injectd BCV signal, we 
   * have to get the same final frequnecy for both template and signal */
  if (userParam->signal == BCV) 
    fendBCV = list->fFinal =  randIn->param.fFinal; 
  
		  
  overlapin->ifExtOutput = 0;
  overlapin->param.fFinal   =  fendBCV;
  overlapin->param.fCutoff  =  fendBCV;               
  overlapin->param 	= *list;	   
 
  overlapout->rhoMax  = -1;
  kMax = floor(fendBCV / df);
  
  /* Extract some values a11,a21,a22 parameter for template creation*/

  
  /* Template creation */
  LAL_CALL(LALCreateBCVFilters(FilterBCV1,
			    FilterBCV2,
			    powerVector,
			    moments,
			    kMin,kMax,     
			    list->psi0,
			    list->psi3
			    ), status->statusPtr);
  
  
  /* The overlap given two filters and the input signal*/
  /*BEInitOverlapOutputIn(&OverlapOutputThisTemplate); */
  /*be sure rhomax = 0 and so on*/
  LAL_CALL(LALWaveOverlapBCV(status->statusPtr, 
			     correlation,
			     overlapin,
			     FilterBCV1, 
			     FilterBCV2, 
			     *userParam, 
			     overlapout, 
			     moments), 
	   status->statusPtr);
  
  
  DETATCHSTATUSPTR( status );
  RETURN( status );  

}

void LALInspiralOverlapBCVSpin(LALStatus *status,
			       InspiralTemplate       *list,
			       BEPowerVector          *powerVector,
			       UserParametersIn       *userParam, 
			       RandomInspiralSignalIn *randIn,
			       REAL4Vector            *FilterBCVSpin1,
			       REAL4Vector            *FilterBCVSpin2,
			       REAL4Vector            *FilterBCVSpin3,
			       InspiralWaveOverlapIn  *overlapin,
			       OverlapOutputIn        *overlapout,
			       REAL4Vector            *correlation,
			       BEMoments              *moments)
  
  
{
  REAL4 fendBCV ;
  REAL4 df;
  INT4 n, kMin, kMax; 

  INITSTATUS (status, "LALInspiralOverlapBCV", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);

  n    = FilterBCVSpin1->length;
  df = randIn->param.tSampling / (REAL8)(n/2) / 2.;      

  kMin  = floor(randIn->param.fLower /df );
  fendBCV = list->fFinal;


  /* if we want to test the BCV metric by injectd BCV signal, we 
   * have to get the same final frequnecy for both template and signal */
  if (userParam->signal == BCVSpin) 
    fendBCV = list->fFinal =  randIn->param.fFinal; 
  
		  
  overlapin->ifExtOutput = 0;
  overlapin->param.fFinal   =  fendBCV;
  overlapin->param.fCutoff  =  fendBCV;               
  overlapin->param 	= *list;	   
 
  overlapout->rhoMax  = -1;
  kMax = floor(fendBCV / df);
  
  /* Extract some values a11,a21,a22 parameter for template creation*/

  
  /* Template creation */
  LAL_CALL(LALCreateBCVSpinFilters(FilterBCVSpin1,
				   FilterBCVSpin2,
				   FilterBCVSpin3,
				   overlapin,
				   powerVector,
				   moments,
				   kMin,kMax     
				   ), status->statusPtr);
  
  
  /* The overlap given two filters and the input signal*/
  /*BEInitOverlapOutputIn(&OverlapOutputThisTemplate); */
  /*be sure rhomax = 0 and so on*/
  LAL_CALL(LALWaveOverlapBCVSpin(status->statusPtr, 
				 correlation,
				 overlapin,
				 FilterBCVSpin1, 
				 FilterBCVSpin2, 
				 FilterBCVSpin3,
				 *userParam, 
				 overlapout),
	   status->statusPtr);
  
  
  DETATCHSTATUSPTR( status );
  RETURN( status );  

}





void
PrintParameters(InspiralCoarseBankIn 	coarse,
		RandomInspiralSignalIn 	rand,
		UserParametersIn    	other)
{

  fprintf(stdout, "coarseBankIn parameters\n"
	          "-----------------------\n");

  fprintf(stdout, "fLower\t\t=\t%f\n"
	  "fUpper\t\t=\t%f\n"
	  "tSampling\t=\t%f\n"
	  "space\t\t=\t%d\n"
	  "mmCoarse\t=\t%f\n" 
	  "mmFine\t\t=\t%f\n" 
	  "iflso\t\t=\t%d\n" 
	  "mMin\t\t=\t%f\n" 
	  "mMax\t\t=\t%f\n" 
	  "MMax\t\t=\t%f\n" 
	  "massRange\t=\t%d\n" 
	  "etamin\t\t=\t%f\n" 
	  "psi0Min\t\t=\t%f\n" 
	  "psi0Max\t\t=\t%f\n" 
	  "psi3Min\t\t=\t%f\n" 
	  "psi3Max\t\t=\t%f\n" 
	  "alpha\t\t=\t%f\n" 
	  "numFcut\t\t=\t%d\n" 
	  "approximant\t=\t%d\n" 
	  "order\t\t=\t%d\n" 
	  "LowGM\t\t=\t%f\n" 
	  "HighGM\t\t=\t%f\n" 
	  "gridType\t=\t%d\n",
	  coarse.fLower,    	
	  coarse.fUpper,    	
	  coarse.tSampling, 	
	  (INT4)(coarse.space),     	
	  coarse.mmCoarse,
	  coarse.mmFine,
	  coarse.iflso,           
	  coarse.mMin,      	
	  coarse.mMax,      	
	  coarse.MMax,      	
	  coarse.massRange,       
	  coarse.etamin, 	        
	  coarse.psi0Min,   	
	  coarse.psi0Max,   	
	  coarse.psi3Min,   	
	  coarse.psi3Max,   	
	  coarse.alpha,     	
	  coarse.numFcutTemplates, 
	  coarse.approximant,     
	  coarse.order,           
	  (REAL8)coarse.LowGM,     	
	  (REAL8)coarse.HighGM,    	
	  (INT4)coarse.gridSpacing);

  
  fprintf(stdout, "bankIn parameters\n"
	  "------------------------\n");

  fprintf(stdout, "useed\t\t=\t%d\n"
	  "type\t\t=\t%d\n"
	  "SignalAmp\t=\t%f\n"
	  "param.order\t=\t%d\n"
	  "param.alpha\t=\t%f\n"
	  "param.ieta\t=\t%d\n"
	  "param.mass1\t=\t%f\n"
	  "param.mass2\t=\t%f\n"
	  "param.fLower\t=\t%f\n"
	  "param.OmegaS\t=\t%f\n"
	  "param.Theta\t=\t%f\n"
	  "mMin\t\t=\t%f\n"
	  "mMax\t\t=\t%f\n"
	  "MMax\t\t=\t%f\n"
	  "etaMin\t\t=\t%f\n"
	  "psi0Min\t\t=\t%f\n"
	  "psi0Max\t\t=\t%f\n"
	  "psi3Min\t\t=\t%f\n"
	  "psi3Max\t\t=\t%f\n"
	  "param.approximant=\t%d\n"
	  "param.tSampling\t=\t%f\n"
	  "param.fCutoff\t=\t%f\n"
	  "param.startTime\t=\t%f\n"
	  "param.startPhase=\t%f\n"
	  "param.nStartPad\t=\t%d\n"
	  "param.signalAmplitude\t=\t%f\n"
	  "param.nEndPad\t=\t%d\n"
	  "NoiseAmp\t=\t%f\n",
	  rand.useed,
	  rand.type,
	  rand.SignalAmp,
	  rand.param.order,
	  rand.param.alpha,
	  rand.param.ieta,
	  rand.param.mass1,
	  rand.param.mass2,
	  rand.param.fLower,
	  rand.param.OmegaS,
	  rand.param.Theta,
	  rand.mMin,
	  rand.mMax,
	  rand.MMax,
	  rand.etaMin,
	  rand.psi0Min,
	  rand.psi0Max,
	  rand.psi3Min,
	  rand.psi3Max,
	  rand.param.approximant,
	  rand.param.tSampling,
	  rand.param.fCutoff,
	  rand.param.startTime,
	  rand.param.startPhase,
	  rand.param.nStartPad,
	  rand.param.signalAmplitude,
	  rand.param.nEndPad,
	  rand.NoiseAmp);

  

  fprintf(stdout, "userParam parameters\n"
	  "------------------------\n");
  
   fprintf(stdout,
	  "alphaFConstraint\t%d\t\n"
	  "extraFinalPrinting\t%d\t\n"
	  "template\t\t%d\t\n"
	  "signalfFinal\t\t%f\t\n"
	  "lalDebug\t\t%d\t\n"
	  "signal\t\t\t%d\t\n"
	  "m1\t\t\t%f\t\n"
	  "m2\t\t\t%f\t\n"
	  "numSeconds\t\t%d\t\n"
	  "psi0\t\t%f\t\n"
	  "psi3\t\t%f\t\n"
	  "tau0\t\t%f\t\n"
	  "tau3\t\t%f\t\n" 
	  "printBestOverlap\t%d\t\n"	
	  "printBestTemplate\t%d\t\n" 	
	  "printSNRHisto\t\t%d\t\n"        
	  "printPsd\t\t%d\t\n"             
	  "printBank\t\t%d\t\n"    	
	  "printResultXml\t\t%d\t\n"    	
	  "printPrototype\t\t%d\t\n"    	
	  "faithfulness\t\t%d\t\n"        	
	  "ntrials\t\t%d\t\n" 	     	
	  "fastSimulation\t\t%d\t\n"       
	  "noiseModel\t\t%d\t\n"           
	  "binaryInjection\t\t%d\t\n"      
	  "maxTotalMass\t\t%f\t\n"         
	  "detector\t\t%d\t\n" 
	  "run\t\t%d\t\n"      
	  "chanName\t\t%s\t\n" 
	  "calCacheName\t\t%s\t\n" 
	  "startTime\t\t%d\t\n"  
	  "numSeconds\t\t%d\t\n" 
	  "L1.chanName\t\t%s\t\n"             
	  "H1.chanName\t\t%s\t\n"             
	  "H2.chanName\t\t%s\t\n"               
	  "L1.dataFile.S3.calCacheName\t%s\t\n"         
	  "L1.dataFile.S3.frInCacheName\t%s\t\n"        
	  "H1.dataFile.S3.calCacheName\t%s\t\n"         
	  "H1.dataFile.S3.frInCacheName\t%s\t\n"        
	  "H2.dataFile.S3.calCacheName\t%s\t\n"         
	  "H2.dataFile.S3.frInCacheName\t%s\t\n"        
	  "L1.dataFile.S2.calCacheName\t%s\t\n"         
	  "L1.dataFile.S2.frInCacheName\t%s\t\n"        
	  "H1.dataFile.S2.calCacheName\t%s\t\n"         
	  "H1.dataFile.S2.frInCacheName\t%s\t\n"        
	  "H2.dataFile.S2.calCacheName\t%s\t\n"         
	  "H2.dataFile.S2.frInCacheName\t%s\t\n",
	  other.alphaFConstraint,
	  other.extraFinalPrinting,
	  other.template,
	  other.signalfFinal,
	  lalDebugLevel,
	  other.signal,
	  other.m1,
	  other.m2,
	  other.numSeconds,
	  other.psi0,
	  other.psi3,
	  other.tau0,
	  other.tau3, 
	  other.printBestOverlap,	
	  other.printBestTemplate, 	
	  other.printSNRHisto,        
	  other.printPsd,             
	  other.printBank,    	
	  other.printResultXml,    	
	  other.printPrototype,    	
	  other.faithfulness,        	
	  other.ntrials, 	     	
	  other.fastSimulation,       
	  other.noiseModel,           
	  other.binaryInjection,      
	  other.maxTotalMass,         
	  other.detector, 
	  other.run,      
	  other.chanName, 
	  other.calCacheName, 
	  other.startTime,  
	  other.numSeconds, 
	  other.L1.chanName,             
	  other.H1.chanName,             
	  other.H2.chanName,               
	  other.L1.dataFile.S3.calCacheName,         
	  other.L1.dataFile.S3.frInCacheName,        
	  other.H1.dataFile.S3.calCacheName,         
	  other.H1.dataFile.S3.frInCacheName,        
	  other.H2.dataFile.S3.calCacheName,         
	  other.H2.dataFile.S3.frInCacheName,        
	  other.L1.dataFile.S2.calCacheName,         
	  other.L1.dataFile.S2.frInCacheName,        
	  other.H1.dataFile.S2.calCacheName,         
	  other.H1.dataFile.S2.frInCacheName,        
	  other.H2.dataFile.S2.calCacheName,         
	  other.H2.dataFile.S2.frInCacheName);
}



void LALBankPrintAscii(MetadataTable          templateBank ,
		       UINT4                  numCoarse,
		       InspiralCoarseBankIn   coarseBankIn )
{
  FILE 		*output;
  SnglInspiralTable            *tmpltCurrent  = NULL;


  output = fopen(BANKEFFICIENCY_PRINTBANK_FILEASCII, "w");

  fprintf(output,"#Number of Coarse Bank Templates=%d\n",numCoarse);
  if (coarseBankIn.approximant == BCV)
    {
      fprintf(output, "#psi0Min=%e, psi0Max=%e, psi3Min=%e, psi3Max=%e\n", 
	      coarseBankIn.psi0Min, coarseBankIn.psi0Max, coarseBankIn.psi3Min, coarseBankIn.psi3Max);
      fprintf(output, "#psi0 psi3  totalMass fFinal\n");
      
    }
  else
    {
      fprintf(output, "#mMin=%e, mMax=%e\n", coarseBankIn.mMin,coarseBankIn.mMax);
      fprintf(output, "#tau0, tau3, mass1, mass2\n");
    }

  for (tmpltCurrent = templateBank.snglInspiralTable;
       tmpltCurrent ;
       tmpltCurrent = tmpltCurrent->next)
    {
      
      if (coarseBankIn.approximant == BCV)
	{      		   
	  fprintf(output, "%e %e %e \n",
		  tmpltCurrent->psi0, 
		  tmpltCurrent->psi3,
		  tmpltCurrent->f_final);
	} 
      else if (coarseBankIn.approximant == BCVSpin )
	{	  
	  fprintf(output, "%e %e %e %e\n", 	
		  tmpltCurrent->psi0, 
		  tmpltCurrent->psi3,
		  tmpltCurrent->beta,
		  tmpltCurrent->f_final);
	}
      else 
	{	  
	  fprintf(output, "%e %e %e %e\n", 	
		  tmpltCurrent->tau0, 
		  tmpltCurrent->tau3, 
		  tmpltCurrent->mass1, 
		  tmpltCurrent->mass2);
	}
    	   
    }
  fprintf(output,"&\n");
  fclose(output);
    
}





void LALBankPrintXML(MetadataTable templateBank ,
		     InspiralCoarseBankIn   coarseBankIn,
		     RandomInspiralSignalIn randIn,
		     UserParametersIn           userParam)

{
  LALStatus             status = blank_status;
  MetadataTable         proctable;
  CHAR                  ifo[3];         /* two character ifo code       */
  LIGOLwXMLStream       xmlStream;
  CHAR                  fname[256];
  LIGOTimeGPS gpsStartTime 	= { 0, 0 }; /* input data GPS start time    */
  LIGOTimeGPS gpsEndTime 	= { 0, 0 }; /* input data GPS end time      */
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;
  CHAR  comment[LIGOMETA_COMMENT_MAX];
  CHAR  ifoName[3][LIGOMETA_IFO_MAX];
  MetadataTable         processParamsTable;
  ProcessParamsTable   *this_proc_param = NULL;


  strncpy( ifoName[0], "no", LIGOMETA_IFO_MAX * sizeof(CHAR) );
  strncpy( ifoName[1], "ne", LIGOMETA_IFO_MAX * sizeof(CHAR) );
  memset( ifo, 0, sizeof(ifo) );
  memcpy( ifo, "MC", sizeof(ifo) - 1 );


  /* --- first we create the filename --- */
  snprintf( fname, sizeof(fname), "TMPLTBANK.xml" ,
	       ifo, gpsStartTime.gpsSeconds,
	       gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  
  /* -- we start to fill the xml file here --- */
  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );

  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlStream, fname), &status );

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
  
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->start_time),
			    &accuracy ), &status );
  LAL_CALL( populate_process_table( &status, proctable.processTable, 
				    PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE ), &status );
  this_proc_param = processParamsTable.processParamsTable = 
    (ProcessParamsTable *) calloc( 1, sizeof(ProcessParamsTable) );

  BEFillProc(this_proc_param, coarseBankIn, randIn, userParam);

  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );



  /* write process table */
  snprintf( proctable.processTable->ifos, LIGOMETA_IFOS_MAX, "%s%s", 
	       ifoName[0], ifoName[1] );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->end_time),
			    &accuracy ), &status );
  
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, process_table ), 
	    &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, proctable, 
				    process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );
  
  /* write process params table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, 
				    process_params_table ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, processParamsTable, 
				    process_params_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );
  
  /* finally write sngl inspiral table */
  if ( templateBank.snglInspiralTable )
    {
      LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, sngl_inspiral_table ), 
		&status );
      LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, templateBank, 
					sngl_inspiral_table ), &status );
      LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );
    }
  
  /* close the output xml file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlStream ), &status );






}

void
 CreateListfromTmplt(InspiralTemplate  *params, 
		     SnglInspiralTable *tmplt)
{
  /* Right now only t03 and psi0psi3 bank exists so we 
     only need those information in principle and nothing else*/
  params->mass1 = tmplt->mass1; 
  params->mass2 = tmplt->mass2; 
  params->fFinal = tmplt->f_final; 
  params->t0 = tmplt->tau0;
  params->t3 = tmplt->tau3;
  params->psi0 = tmplt->psi0;
  params->psi3 = tmplt->psi3;
  params->beta = tmplt->beta;


  /*      tmplt->mchirp  = (REAL4) coarseList[0].params.chirpMass;
      tmplt->eta     = (REAL4) coarseList[0].params.eta;
      tmplt->tau4    = (REAL4) coarseList[0].params.t4;
      tmplt->tau5    = (REAL4) coarseList[0].params.t5;
      tmplt->ttotal  = (REAL4) coarseList[0].params.tC;

  */  


}


/* ****************************************************************************
 * Initialization of the parameters fromn random signal, bank and others 
 * structures
 * ***************************************************************************/
void ParametersInitialization(	InspiralCoarseBankIn 	*coarseBankIn, 
				RandomInspiralSignalIn	*randIn, 
				UserParametersIn	*userParam)
{
  
  InitInspiralCoarseBankIn(coarseBankIn);
  InitRandomInspiralSignalIn(randIn);
  InitUserParametersIn(userParam);
  
}


/* ****************************************************************************
 * CoarseIn initialization
 ***************************************************************************/
void InitInspiralCoarseBankIn(InspiralCoarseBankIn *coarseBankIn)
{

  coarseBankIn->fLower    	= BANKEFFICIENCY_FLOWER;
  coarseBankIn->fUpper    	= -1.;
  coarseBankIn->tSampling 	= BANKEFFICIENCY_TSAMPLING;
  coarseBankIn->space     	= BANKEFFICIENCY_SPACE;
  coarseBankIn->mmCoarse  	= BANKEFFICIENCY_MMCOARSE;
  coarseBankIn->mmFine    	= BANKEFFICIENCY_MMFINE;
  coarseBankIn->iflso           = BANKEFFICIENCY_IFLSO ;
  coarseBankIn->mMin      	= BANKEFFICIENCY_MMIN;
  coarseBankIn->mMax      	= BANKEFFICIENCY_MMAX;
  coarseBankIn->MMax      	= BANKEFFICIENCY_MMAX * 2;
  coarseBankIn->massRange       = MinMaxComponentMass; 
  coarseBankIn->etamin 	        = coarseBankIn->mMin * coarseBankIn->mMax  
			        / pow(coarseBankIn->MMax, 2.);

  coarseBankIn->psi0Min   	= BANKEFFICIENCY_PSI0MIN;
  coarseBankIn->psi0Max   	= BANKEFFICIENCY_PSI0MAX;
  coarseBankIn->psi3Min   	= BANKEFFICIENCY_PSI3MIN;
  coarseBankIn->psi3Max   	= BANKEFFICIENCY_PSI3MAX;
  coarseBankIn->alpha     	= BANKEFFICIENCY_ALPHABANK;
  /*unsigned variable ... */
  coarseBankIn->numFcutTemplates = BANKEFFICIENCY_NFCUT;
  coarseBankIn->approximant     = BANKEFFICIENCY_TEMPLATE;  
  coarseBankIn->order           = BANKEFFICIENCY_ORDER_TEMPLATE;  
  coarseBankIn->LowGM     	= BANKEFFICIENCY_LOWGM;
  coarseBankIn->insidePolygon  	= BANKEFFICIENCY_INSIDEPOLYGON;
  coarseBankIn->HighGM    	= BANKEFFICIENCY_HIGHGM;
  coarseBankIn->gridSpacing     = SquareNotOriented;
}


/* ****************************************************************************
 * Random initialization 
 * ****************************************************************************/
void InitRandomInspiralSignalIn(RandomInspiralSignalIn *randIn)
{
  randIn->useed         	= BANKEFFICIENCY_USEED;				/* seed for random MC simulation 	*/ 
  randIn->type          	= BANKEFFICIENCY_TYPE;				/* type of simulation 			*/
  randIn->SignalAmp     	= BANKEFFICIENCY_SIGNALAMPLITUDE;		/* if x = n +h 				*/
  randIn->param.order   	= BANKEFFICIENCY_ORDER_SIGNAL;   		/* and its order			*/
  randIn->param.alpha		= BANKEFFICIENCY_ALPHASIGNAL; 			/* alpha value of BCV 			*/
  randIn->param.ieta    	= BANKEFFICIENCY_IETA; 				/*					*/
  randIn->param.mass1 		= BANKEFFICIENCY_MMIN;				/* To aalocate memory			*/ 
  randIn->param.mass2  		= BANKEFFICIENCY_MMIN; 				/* idem 				*/
  randIn->param.fLower  	= BANKEFFICIENCY_FLOWER;			/* Lower cutoff freq. of the template	*/
  randIn->param.OmegaS		= 0.;						/* EOB parameter 			*/
  randIn->param.Theta		= 0.;						/* EOB parameter 			*/
  randIn->mMin          	= BANKEFFICIENCY_MMIN;				/* min mass to inject			*/
  randIn->mMax          	= BANKEFFICIENCY_MMAX;				/* max mass to inject 			*/
  randIn->MMax          	= BANKEFFICIENCY_MMAX * 2;			/* total mass max 			*/
  randIn->t0Min                 = 0.;                                           /* min tau0 to inject                   */
  randIn->t0Max                 = 0.;                                           /* max tau0 to inject                   */
  randIn->t0Min                 = 0.;                                           /* min tau3 to inject                   */
  randIn->t0Max                 = 0.;                                           /* max tau3 to inject                   */
  randIn->etaMin                = BANKEFFICIENCY_MMIN
    * (BANKEFFICIENCY_MMAX - BANKEFFICIENCY_MMIN) 
    / (BANKEFFICIENCY_MMAX * 2) / (BANKEFFICIENCY_MMAX * 2);
  randIn->psi0Min  		= BANKEFFICIENCY_PSI0MIN;			/* psi0 range 				*/
  randIn->psi0Max  		= BANKEFFICIENCY_PSI0MAX;	
  randIn->psi3Min  		= BANKEFFICIENCY_PSI3MIN;	
  randIn->psi3Max  		= BANKEFFICIENCY_PSI3MAX;	
  randIn->param.approximant 	= BANKEFFICIENCY_SIGNAL;			/* approximant of h, the signal		*/
  randIn->param.tSampling 	= BANKEFFICIENCY_TSAMPLING;			/* sampling 				*/
  randIn->param.fCutoff		= BANKEFFICIENCY_TSAMPLING/2.-1;		/* sampling 				*/
  randIn->param.startTime       = BANKEFFICIENCY_STARTTIME; 
  randIn->param.startPhase      = BANKEFFICIENCY_STARTPHASE; 
  randIn->param.nStartPad       = BANKEFFICIENCY_NSTARTPAD;
  randIn->param.signalAmplitude = BANKEFFICIENCY_SIGNALAMPLITUDE;
  randIn->param.nEndPad         = BANKEFFICIENCY_NENDPAD;
  randIn->NoiseAmp              = BANKEFFICIENCY_NOISEAMPLITUDE;
}


/* ****************************************************************************
 * Other Param initialization							
 * ***************************************************************************/
void InitUserParametersIn(UserParametersIn *userParam)
{
  userParam->alphaFConstraint   = ALPHAFConstraint;
  userParam->extraFinalPrinting = 0; 
  userParam->template		= BANKEFFICIENCY_TEMPLATE;
  userParam->signalfFinal       =  BANKEFFICIENCY_TSAMPLING/2.-1;
  userParam->signal       	= BANKEFFICIENCY_SIGNAL;
  userParam->m1           	= -1;
  userParam->m2           	= -1;
  userParam->template_m1        = -1;
  userParam->template_m2        = -1;
  userParam->numSeconds         = -1;
  userParam->psi0         	= -1;
  userParam->psi3         	= -1;
  userParam->tau0         	= -1;
  userParam->tau3         	= -1;
  userParam->printBestOverlap 	= BANKEFFICIENCY_PRINTBESTOVERLAP;
  userParam->printBestTemplate 	= BANKEFFICIENCY_PRINTBESTTEMPLATE;
  userParam->printSNRHisto        = BANKEFFICIENCY_PRINTSNRHISTO;
  userParam->printPsd             = BANKEFFICIENCY_PRINTPSD;
  userParam->printBank    	= BANKEFFICIENCY_PRINTBANK;
  userParam->printResultXml    	= BANKEFFICIENCY_PRINTRESULTXML;
  userParam->printPrototype    	= BANKEFFICIENCY_PRINTPROTOTYPE;
  userParam->faithfulness        	= BANKEFFICIENCY_FAITHFULNESS;
  userParam->ntrials 	     	= BANKEFFICIENCY_NTRIALS;
  userParam->fastSimulation       = BANKEFFICIENCY_FASTSIMULATION;
  userParam->noiseModel           = LIGOI;
  userParam->binaryInjection      = NoUserChoice; 
  userParam->maxTotalMass         = -1;
  userParam->startPhase           = 1;

  userParam->detector     = L1;
  userParam->run          = S3;
  userParam->chanName     = NULL;
  userParam->calCacheName = NULL;
  userParam->calCacheName = NULL;
  userParam->startTime    = 751956568;
  userParam->numSeconds   = -1;
  userParam->realNoise    = 0;
  userParam->inputPSD     = NULL;
  
  userParam->L1.chanName             = "L1:LSC-AS_Q";
  userParam->H1.chanName             = "H1:LSC-AS_Q";
  userParam->H2.chanName             = "H2:LSC-AS_Q";
  
  userParam->L1.dataFile.S4.calCacheName         =  "/home/spxcar/inspiralRuns/Calibration/L1-CAL-V03-793126813-796025053.cache";
  userParam->L1.dataFile.S4.frInCacheName        =  "/home/spxcar/inspiralRuns/cacheFiles/L-RDS_R_L3-793131485-795679161.cache";     
  userParam->H1.dataFile.S4.calCacheName         =  "/home/spxcar/inspiralRuns/Calibration/H1-CAL-V03-793126813-796025053.cache";
  userParam->H1.dataFile.S4.frInCacheName        =  "/home/spxcar/inspiralRuns/cacheFiles/H-RDS_R_L3-793139225-795679221.cache";     
  userParam->H2.dataFile.S4.calCacheName         =  "/home/spxcar/inspiralRuns/Calibration/H2-CAL-V03-793126813-796025053.cache";
  userParam->H2.dataFile.S4.frInCacheName        =  "/home/spxcar/inspiralRuns/cacheFiles/H-RDS_R_L3-793139225-795679221.cache";

  userParam->L1.dataFile.S3.calCacheName         =  "/home/spxcar/cache_files/L1-CAL-V02-729273600-734367600.cache";
  userParam->L1.dataFile.S3.frInCacheName        =  "/home/spxcar/cachefiles/inspiralRuns/cacheFiles/CacheFile_L_S3_RDS_R_L3.txt";     
  userParam->H1.dataFile.S3.calCacheName         =  "/home/spxcar/cache_files/H1-CAL-V02-751651244-757699245.cache";
  userParam->H1.dataFile.S3.frInCacheName        =  "/home/spxcar/cachefiles/inspiralRuns/cacheFiles/CacheFile_H_S3_RDS_R_L3.txt";     
  userParam->H2.dataFile.S3.calCacheName         =  "/home/spxcar/cache_files/H2-CAL-V02-751651244-757699245.cache";
  userParam->H2.dataFile.S3.frInCacheName        =  "/home/spxcar/cachefiles/inspiralRuns/cacheFiles/CacheFile_H_S3_RDS_R_L3.txt";     

  userParam->L1.dataFile.S2.calCacheName         =  "/home/spxcar/inspiralRuns/Calibrationcache_files/L1-CAL-V03-729273600-734367600.cache";
  userParam->L1.dataFile.S2.frInCacheName        =  "/home/spxcar/inspiralRuns//cacheFiles/CacheFile_L_S2_RDS_R_L3.txt";     
  userParam->H1.dataFile.S2.calCacheName         =  "/home/spxcar/inspiralRuns/Calibrationcache_files/H1-CAL-V03-729273600-734367600.cache";
  userParam->H1.dataFile.S2.frInCacheName        =  "/home/spxcar/inspiralRuns/cacheFiles/CacheFile_H_S2_RDS_R_L3.txt";     
  userParam->H2.dataFile.S2.calCacheName         =  "/home/spxcar/inspiralRuns/Calibration/cache_files/H2-CAL-V03-731849076-734367576.cache";
  userParam->H2.dataFile.S2.frInCacheName        =  "/home/spxcar/inspiralRuns/cacheFiles/CacheFile_H_S2_RDS_R_L3.txt";      
  userParam->dataCheckPoint                      = 0;
}



/* ***************************************************************************
 *  Function to parse user parameters 
 *  ************************************************************************ */
void 
ParseParameters(	INT4 			*argc, 
			CHAR 			**argv,
			InspiralCoarseBankIn 	*coarseBankIn,
			RandomInspiralSignalIn 	*randIn,
			UserParametersIn    	*userParam)
{
  INT4  i = 1;
  REAL8 tmp1, tmp2;
  
  while(i < *argc)
    {
      if (!strcmp(argv[i],	"--bank-alpha")) {
        BEParseGetDouble(argv,  &i, &(coarseBankIn->alpha)); 
      }
      else if (!strcmp(argv[i],"--bank-fcut-range")) {
	{
	  BEParseGetDouble2(argv,  &i, &tmp1, &tmp2);
	  coarseBankIn->LowGM  = (REAL4)tmp1; 
	  coarseBankIn->HighGM = (REAL4)tmp2;
	}
      }
      else if (!strcmp(argv[i],	"--bank-ffinal")) {
        BEParseGetDouble(argv,  &i, &(coarseBankIn->fUpper)); 
      }    
      else if (!strcmp(argv[i],	"--bank-inside-polygon")) {
        BEParseGetInt(argv,  &i, (INT4*)&(coarseBankIn->insidePolygon)); 
      }

      else if (!strcmp(argv[i], "--bank-grid-spacing")) {
	BEParseGetString(argv, &i);
	if (!strcmp(argv[i],"Square")) {
	  coarseBankIn->gridSpacing = Square;
	}
	else if (!strcmp(argv[i], "Hexagonal")) { 
	  coarseBankIn->gridSpacing = Hexagonal;
	}
	else if (!strcmp(argv[i], "HexagonalNotOriented")) {
	  coarseBankIn->gridSpacing = HexagonalNotOriented;
	}
	else if (!strcmp(argv[i], "SquareNotOriented")) {
	  coarseBankIn->gridSpacing = SquareNotOriented;
	}
	else {
	  coarseBankIn->gridSpacing = None;        	  	
	}
      }        
      else if (!strcmp(argv[i],	"--bank-number-fcut")) {
        BEParseGetInt(argv,  &i, (INT4*)&(coarseBankIn->numFcutTemplates)); 
      }
      else if (!strcmp(argv[i],"--bank-mass-range")) {
	BEParseGetDouble2(argv,  &i, &(coarseBankIn->mMin), &(coarseBankIn->mMax));
      }
      else if (!strcmp(argv[i],	"--bank-psi0-range")) {
        BEParseGetDouble2(argv,  &i, &(coarseBankIn->psi0Min), &(coarseBankIn->psi0Max));
      }
      else if (!strcmp(argv[i],	"--bank-psi3-range")) {
        BEParseGetDouble2(argv,  &i, &(coarseBankIn->psi3Min), &(coarseBankIn->psi3Max));
      }
      else if (!strcmp(argv[i], "--channel")) {
	  BEParseGetString(argv, &i);
	  userParam->chanName = argv[i];
	}
      else if (!strcmp(argv[i],"--debug")) {
        BEParseGetInt(argv,  &i, &(lalDebugLevel)); 
      }
      else if (!strcmp(argv[i], "dataCheckPoint")) {
	userParam->dataCheckPoint = 1;	
      }
      else if (!strcmp(argv[i], "--detector")) {
	BEParseGetString(argv, &i);
	if (!strcmp(argv[i],  "L1")) {
	  userParam->detector = L1;
	}
	else if (!strcmp(argv[i],"H1")) {
	  userParam->detector = H1;
	}
	else if (!strcmp(argv[i],"H2")) {
	  userParam->detector = H2;
	}
	else {
	  userParam->detector  = None;
	}
      }
      else if (!strcmp(argv[i],"--fl-signal")) {
        BEParseGetDouble(argv,  &i, &(randIn->param.fLower));
      }
      else if  (!strcmp(argv[i],"--fl-template")) {
        BEParseGetDouble(argv,  &i, &(coarseBankIn->fLower));
      }
      else if  (!strcmp(argv[i],"--fl")) {
        BEParseGetDouble(argv,  &i, &(coarseBankIn->fLower));
	randIn->param.fLower = coarseBankIn->fLower;
      }  
      else if (!strcmp(argv[i], "--gps-start-time")) {
	BEParseGetInt(argv, &i , &(userParam->startTime));
      }
      else if (!strcmp(argv[i], "--help") || !strcmp(argv[i],"--h")) {
	Help();
      }
      else if (!strcmp(argv[i],"--max-total-mass")) {
	
	BEParseGetDouble(argv,  &i, &tmp1);
	userParam->maxTotalMass = (REAL4)tmp1;	
      }
      else if (!strcmp(argv[i],"--ascii2xml")) {       
	ascii2xml = 1;	
      }
      else if (!strcmp(argv[i], "--m1")){
        BEParseGetDouble(argv,  &i, &(userParam->m1));
      }
      else if (!strcmp(argv[i], "--m2")) {
        BEParseGetDouble(argv,  &i, &(userParam->m2));
      }
      else if (!strcmp(argv[i], "--template-m1")){
        BEParseGetDouble(argv,  &i, &(userParam->template_m1));
      }
      else if (!strcmp(argv[i], "--template-m2")) {
        BEParseGetDouble(argv,  &i, &(userParam->template_m2));
      }
      else if (!strcmp(argv[i],	"--mm")){
        BEParseGetDouble(argv,  &i, &(coarseBankIn->mmCoarse)); 
      }
      else if (!strcmp(argv[i],	"--n") || !strcmp(argv[i], "--ntrial")){
        BEParseGetInt(argv,  &i, &(userParam->ntrials));
      }
      else if (!strcmp(argv[i],	"--noise-amplitude")) {
        BEParseGetDouble(argv,  &i, &(randIn->NoiseAmp)); 
      }
      else if (!strcmp(argv[i],	"--noise-model")){	
	BEParseGetString(argv, &i);
	if (!strcmp(argv[i], "LIGOI")) {
	  userParam->noiseModel = LIGOI;
	}
	else if (!strcmp(argv[i], "LIGOA")) {
	  userParam->noiseModel = LIGOA;
	}
	else if (!strcmp(argv[i], "VIRGO")) {
	  userParam->noiseModel = VIRGO;
	}
	else if (!strcmp(argv[i], "TAMA")) {
	  userParam->noiseModel = TAMA;
	}
	else if (!strcmp(argv[i], "GEO")) {
	  userParam->noiseModel = GEO;
	}
	else if (!strcmp(argv[i], "UNITY")) {
	  userParam->noiseModel = UNITY;
	}
	else if (!strcmp(argv[i], "REALPSD")) {
	  userParam->noiseModel = REALPSD;
	}
	else if (!strcmp(argv[i], "READPSD")) { 
	  userParam->noiseModel = READPSD;
	  userParam->inputPSD = argv[++i];
	}
	else {
	  userParam->noiseModel = None;	  
	}
      }
      else if (!strcmp(argv[i], "--num-seconds")) {
        BEParseGetInt(argv,  &i, &(userParam->numSeconds));
      }
      else if (!strcmp(argv[i], "--psi0")) {
        BEParseGetDouble(argv,  &i, &(userParam->psi0));
      }
      else if (!strcmp(argv[i], "--psi3")) {
        BEParseGetDouble(argv,  &i, &(userParam->psi3));
      }
      else if (!strcmp(argv[i], "--real-noise")) {
	userParam->realNoise = 1;
      }
      else if (!strcmp(argv[i], "--run")) {
	BEParseGetString(argv, &i);	 
	if (!strcmp(argv[i],  "S2")) {
	  userParam->run = S2;
	}
	else if (!strcmp(argv[i],"S3")){
	    userParam->run = S3;	  
	}
	else if (!strcmp(argv[i],"S4")){
	    userParam->run = S4;	  
	}
	else {
	    userParam->run = None;
	}
      } 
      else if (!strcmp(argv[i],"--sampling")) {
        BEParseGetDouble(argv,  &i, &(coarseBankIn->tSampling));
        randIn->param.tSampling = coarseBankIn->tSampling;
      }   
      else if (!strcmp(argv[i],"--seed")){
        BEParseGetInt(argv,  &i, &(randIn->useed));
      }
      else if (!strcmp(argv[i], "--signal")) {
	BEParseGetString(argv, &i);
	
	if (strcmp(argv[i],	  "TaylorT1")	==0) 		userParam->signal = TaylorT1;
	else if (strcmp(argv[i],"TaylorT2")	==0)		userParam->signal = TaylorT2;
	else if (strcmp(argv[i],"TaylorT3")	==0)	    	userParam->signal = TaylorT3;
	else if (strcmp(argv[i],"TaylorF1")	==0)	    	userParam->signal = TaylorF1;
	else if (strcmp(argv[i],"TaylorF2")	==0)	    	userParam->signal = TaylorF2;
	else if (strcmp(argv[i],"PadeT1")	==0)	    	userParam->signal = PadeT1;
	else if (strcmp(argv[i],"PadeF1")	==0)	    	userParam->signal = PadeF1;
	else if (strcmp(argv[i],"EOB")	==0)    	        userParam->signal = EOB;
	else if (strcmp(argv[i],"BCV")	==0)                    userParam->signal = BCV;
	else if (strcmp(argv[i],"BCVSpin")	==0)            userParam->signal = BCVSpin;
	else if (strcmp(argv[i],"SpinTaylorT3")==0)	    	userParam->signal = SpinTaylorT3;
	else  BEPrintError("wrong approximant\n");
	randIn->param.approximant = userParam->signal;		  
      }
      else if (!strcmp(argv[i],	"--signal-alpha")){	     
	BEParseGetDouble(argv,  &i, &(randIn->param.alpha));
      } 
      else if (!strcmp(argv[i],	"--signal-alpha1")){	     
	BEParseGetDouble(argv,  &i, &(randIn->param.alpha1));
      }
      else if (!strcmp(argv[i],	"--signal-alpha2")){	     
	BEParseGetDouble(argv,  &i, &(randIn->param.alpha2));
      }
      else if (!strcmp(argv[i],	"--signal-beta")){	     
	BEParseGetDouble(argv,  &i, &(randIn->param.beta));
      }
      else if (!strcmp(argv[i],	"--signal-amplitude")) {
	BEParseGetDouble(argv,  &i, &(randIn->SignalAmp));
      }
      else if (!strcmp(argv[i],	"--signal-ffinal")) {
	BEParseGetDouble(argv,  &i, &(randIn->param.fCutoff));
	userParam->signalfFinal = randIn->param.fCutoff ;
      }
      else if (!strcmp(argv[i],	"--signal-mass-range")){	
        BEParseGetDouble2(argv,  &i, &(randIn->mMin), &(randIn->mMax));
	randIn->param.mass1 = randIn->param.mass2 = randIn->mMin; /*default values for injection*/
      }   
      else if (!strcmp(argv[i], "--signal-tau0-range")){
        BEParseGetDouble2(argv,  &i, &(randIn->t0Min), &(randIn->t0Max));
      }
      else if (!strcmp(argv[i], "--signal-tau3-range")){
        BEParseGetDouble2(argv,  &i, &(randIn->tnMin), &(randIn->tnMax));
      }
      else if (!strcmp(argv[i],	"--signal-order")) {
	BEParseGetInt(argv,  &i, (INT4*)&(randIn->param.order));
      }
      else if (!strcmp(argv[i],	"--signal-nstartpad")) {
	BEParseGetInt(argv,  &i, (INT4*)&(randIn->param.nStartPad));
      }
      else if (!strcmp(argv[i],	"--signal-psi0-range")) {
        BEParseGetDouble2(argv,  &i, &(randIn->psi0Min), &(randIn->psi0Max));
      }
      else if (!strcmp(argv[i],	"--signal-psi3-range")){
        BEParseGetDouble2(argv,  &i, &(randIn->psi3Min), &(randIn->psi3Max));
      }  
      else if (!strcmp(argv[i],	"--signal-random-nstartpad")) {
	randnStartPad = 1;
      }
      else if (!strcmp(argv[i],	"--simulation-type")) {
	BEParseGetString(argv, &i);
	if (!strcmp(argv[i],"NoiseOnly"))	
	  randIn->type = 1;
	else if (!strcmp(argv[i],"SignalOnly"))
	  randIn->type = 0;
	else if (!strcmp(argv[i],"NoiseAndSignal"))
	  randIn->type = 2;
	else randIn->type = None;
      }
      else if (!strcmp(argv[i], "--no-start-phase")){
	userParam->startPhase = 0;
      }
      else if (!strcmp(argv[i], "--tau0")) {
        BEParseGetDouble(argv,  &i, &(userParam->tau0));
      }
      else if (!strcmp(argv[i], "--tau3")) {
        BEParseGetDouble(argv,  &i, &(userParam->tau3));
      }
      else if (!strcmp(argv[i],"--template")) {
	BEParseGetString(argv, &i);
	
	if (!strcmp(argv[i],	"TaylorT1")	)	userParam->template = TaylorT1;
	else if (!strcmp(argv[i],	"TaylorT2")	)	userParam->template = TaylorT2;
	else if (!strcmp(argv[i],	"TaylorT3")	)	userParam->template = TaylorT3;
	else if (!strcmp(argv[i],	"TaylorF1")	)	userParam->template = TaylorF1;
	else if (!strcmp(argv[i],	"TaylorF2")	)	userParam->template = TaylorF2;
	else if (!strcmp(argv[i],	"PadeT1")	)	userParam->template = PadeT1;
	else if (!strcmp(argv[i],	"PadeF1")	)	userParam->template = PadeF1;
	else if (!strcmp(argv[i],	"EOB")		)	userParam->template = EOB;
	else if (!strcmp(argv[i],	"BCV")		)       userParam->template = BCV;
	else if (!strcmp(argv[i],	"BCVSpin")	)
	  {
	    userParam->template = BCVSpin;
	    coarseBankIn->betaMin = 0.;
	    coarseBankIn->betaMax = 1000.;
	  }
	else if (!strcmp(argv[i],	"SpinTaylorT3")	)	userParam->template = SpinTaylorT3;
	else BEPrintError("wrong approximant\n");
	
	coarseBankIn->approximant = userParam->template;
	if ( coarseBankIn->approximant == BCV || coarseBankIn->approximant == BCVSpin) 	
	  coarseBankIn->space = Psi0Psi3; 
	else 
	  coarseBankIn->space = Tau0Tau3;
      }
      else if (!strcmp(argv[i],	"--template-order")) {
        BEParseGetInt(argv,  &i, (INT4*)(&(coarseBankIn->order)));
      }
      
      /* no user parameters requeste, only flags below*/               
      else if (!strcmp(argv[i],"--alpha-constraint")) {
        userParam->alphaFConstraint       = ALPHAFConstraint;
      }
      else if (!strcmp(argv[i],"--bhns-injection")){
	userParam->binaryInjection        = BHNS;
      }
      else if (!strcmp(argv[i],"--no-alpha-constraint")) {
        userParam->alphaFConstraint	= ALPHAFUnconstraint;
      }
      else if (!strcmp(argv[i],"--print-best-overlap")) {
        userParam->printBestOverlap 	= 1;
      }
      else if (!strcmp(argv[i],"--faithfulness")) {
        userParam->faithfulness    	= 1;
      }
      else if (!strcmp(argv[i],"--check")) {
        userParam->faithfulness    	= 1;
      }
      else if (!strcmp(argv[i],"--real-noise")) {
        randIn->type    	        = RealNoise;
      }
      else if (!strcmp(argv[i],"--print-psd")) {
        userParam->printPsd 		= 1;
      }
      else if (!strcmp(argv[i],"--print-snr-histo")) {
	userParam->printSNRHisto 	        = 1;
      }
      else if (!strcmp(argv[i],"--verbose")) {
        vrbflg 	= 1;
      }
      else if (!strcmp(argv[i],"--version")) {
	fprintf(stderr, "TestRandomBank code"
		"Thomas Cokelaer, Thomas.Cokelaer@astro.cf.ac.uk\n"
		"CVS Version :" CVS_ID_STRING_C "\n"
		"CVS Version :" CVS_ID_STRING "\n"
		"CVS Tag: " CVS_NAME_STRING "\n");
	exit(0);
      }
      else if (!strcmp(argv[i],"--print-bank")) {
        userParam->printBank		= 1;
      }
      else if (!strcmp(argv[i],"--print-bank")) {
	PrintParameters(*coarseBankIn, *randIn,*userParam);
	exit(0);
      }
      else if (!strcmp(argv[i],"--print-result-xml")) {
        userParam->printResultXml        = 1;
      }
      else if (!strcmp(argv[i],"--print-prototype")) {
        userParam->printPrototype        = 1;
      }
      else if (!strcmp(argv[i],"--fast-simulation")) {
        userParam->fastSimulation        = 1;           
      }	   	    
      else {
        fprintf(stderr, "%s option does not exist. Type --help for options\n", argv[i]);
	exit(0);
      }      
      i++;       
    }   
}






void
BEParseGetInt(  CHAR    **argv, 
                INT4    *index,
                INT4    *data)
{  
  CHAR *tmp;
  CHAR msg[2048];
  CHAR *tmp1;
  tmp1=argv[*index+1];

  if (argv[*index+1]!=NULL)
    {
      *data = strtol(argv[*index+1], &tmp  , 10);
      if (*data==0 && tmp1[0] != '0')
	{
	  sprintf(msg, "Expect a float after option %s (got %s)\n ",
		  argv[*index],
		  argv[*index+1]);
	  BEPrintError(msg);
	}
    }
  else  
    {
      sprintf(msg, "Expect a float after option %s (got %s)\n ",
	      argv[*index],
	      argv[*index+1]);
      BEPrintError(msg);
    }
  *index =*index + 1;
}

void
BEParseGetString(  CHAR    **argv, 
		   INT4    *index)		  
{  
  CHAR msg[2048];
  
  if (argv[*index+1]==NULL)
    {
      sprintf(msg, "Expect a string after %s\n ",
	      argv[*index] );
      BEPrintError(msg);
    }
  *index =*index + 1;
}


void
BEParseGetDouble(CHAR    **argv, 
                INT4    *index,
                REAL8   *data)
{  
  CHAR *tmp;
  CHAR msg[2048];
  CHAR *tmp2 ;
  tmp2 = argv[*index+1];

  if (argv[*index+1] != NULL)
    {
      *data = strtod(argv[*index+1], &tmp  );
      if (*data == 0 && tmp2[0]!='0')
	{
	  sprintf(msg, "Expect a float after option %s (got %s)\n ",
		  argv[*index],
		  argv[*index+1]);
	  BEPrintError(msg);
	}
    }
  else  
    {
      sprintf(msg, "Expect a float after option %s (got %s)\n ",
	      argv[*index],
	      argv[*index+1]);
      BEPrintError(msg);
    }
  *index =*index + 1;
}


void
BEParseGetDouble2(CHAR    **argv, 
		  INT4    *index,
		  REAL8   *data1,
		  REAL8   *data2)
{  
  CHAR *tmp; 
  CHAR msg[2048];
  CHAR *tmp2 , *tmp1;

  tmp1 = argv[*index+1];
  tmp2=  argv[*index+2];

  *data1 = 0 ;
  *data2 = 0 ;

  if (argv[*index+1]!=NULL && argv[*index+2]!=NULL)
    {
      *data1 = strtod(argv[*index+1], &tmp  );
      *data2 = strtod(argv[*index+2], &tmp  );
      if ((!(*data1) && tmp1[0]!='0')
	||  (!(*data2) && tmp2[0]!='0'))
	{   
	  sprintf(msg, "Expect 2 floats after option %s (got %s and %s)\n ",
		  argv[*index],
		  argv[*index+1],argv[*index+2]);
	  BEPrintError(msg);
	}
    }
  else  
    {
      sprintf(msg, "Expect 2 floats after option %s (got %s and %s)\n ",
	      argv[*index],
	      argv[*index+1],argv[*index+2]);
      BEPrintError(msg);
    }
  *index = *index +2 ;  
}  
  
    
    /* --- Function to check validity of some parsed parameters (TODO)--- */
void UpdateParams(InspiralCoarseBankIn *coarseBankIn,
		  RandomInspiralSignalIn *randIn,
		  UserParametersIn *userParam)
{
  REAL4  temp;
  CHAR   msg[2048];

  temp =randIn->param.mass1; /*just to avoid boring warning*/

  userParam->useed = randIn->useed;

  if (coarseBankIn->fUpper == -1 ){
    coarseBankIn->fUpper = coarseBankIn->tSampling/2. -1.;
  }
  
  if (coarseBankIn->fUpper >=coarseBankIn->tSampling/2. -1.){
    coarseBankIn->fUpper = coarseBankIn->tSampling/2 -1.;
  }
 
  randIn->param.fCutoff = coarseBankIn->fUpper;
 
  if (randIn->param.fCutoff >=coarseBankIn->tSampling/2. -1.){
      randIn->param.fCutoff = coarseBankIn->tSampling/2 -1.;
  }
  
  if (coarseBankIn->alpha < 0 ){
    sprintf(msg, "--bank-alpha (%f) parameter must be positive in the range [0,1] \n",
	    coarseBankIn->alpha);
    BEPrintError(msg);
    }

  if (coarseBankIn->fUpper <= coarseBankIn->fLower 
      || coarseBankIn->fUpper >= coarseBankIn->tSampling/2){
    sprintf(msg, "--bank-ffinal (%f) paramter must be greater than  bank-fl (%f) and less than sampling/2 %f\n",
	    coarseBankIn->fUpper,
	    coarseBankIn->fLower ,
	    coarseBankIn->tSampling/2); 
    BEPrintError(msg);
  }
  

  if (coarseBankIn->gridSpacing == None  ){
    sprintf(msg, "--bank-grid-spacing (%d) parameter must be < square, hexagonal, hexagonalOriented, squareOriented> \n",
	    coarseBankIn->gridSpacing);
    
    BEPrintError(msg);
  }

  
  if  ((coarseBankIn->mMin >= coarseBankIn->mMax ) || (coarseBankIn->mMin <=0)){
    sprintf(msg, "--bank-mass-range (%f %f) paramter must be sorted and > 0 \n",
	    coarseBankIn->mMin, coarseBankIn->mMax);
    BEPrintError(msg);
  }
  else /* if valid parameters */
    {
      if (userParam->maxTotalMass != -1)
	{
	  if (userParam->maxTotalMass >= (2 * coarseBankIn->mMax))
	    {
	      /*coarseBankIn->MMax =  2 * userParam->maxTotalMass;*/
	      coarseBankIn->etamin = coarseBankIn->mMin * coarseBankIn->mMax 
		/ (coarseBankIn->mMin + coarseBankIn->mMax) 
		/ (coarseBankIn->mMin + coarseBankIn->mMax);
	    }
	  else
	    {
	      coarseBankIn->MMax =  userParam->maxTotalMass;
	      coarseBankIn->etamin = coarseBankIn->mMin * (coarseBankIn->MMax  - coarseBankIn->mMin)
		/ (coarseBankIn->MMax )
		/ (coarseBankIn->MMax );
	    }
	}
      else
	{
	  coarseBankIn->MMax = 2 * coarseBankIn->mMax;
	  coarseBankIn->etamin = coarseBankIn->mMin * (coarseBankIn->MMax  - coarseBankIn->mMin)
	    / (coarseBankIn->MMax )
	    / (coarseBankIn->MMax );
	}
    }
  
  if (coarseBankIn->psi0Min <=0 || coarseBankIn->psi0Min > coarseBankIn->psi0Max) {
    sprintf(msg, "--bank-psi0-range (%f %f) paramter must be sorted and > 0 \n",
	    coarseBankIn->psi0Min, coarseBankIn->psi0Max);
    BEPrintError(msg);
  } 

  if (coarseBankIn->psi3Min >=0 || coarseBankIn->psi3Min > coarseBankIn->psi3Max) {
    sprintf(msg, "--bank-psi0-range (%f %f) paramter must be sorted and >= 0 \n",
	    coarseBankIn->psi3Min, coarseBankIn->psi3Max);
    BEPrintError(msg);
  }


  if (userParam->detector == None)
    {
      sprintf(msg, "--detector (%d) expect option <H1, H2, L1> only",
	      userParam->detector);
      BEPrintError(msg);
    }

  if (coarseBankIn->LowGM > coarseBankIn->HighGM)
    {
      sprintf(msg, "--bank-fcut-range (%f %f) expect sorted , typically 3 and 6",
	      coarseBankIn->LowGM , coarseBankIn->HighGM);
      BEPrintError(msg);
    }
 
  if (coarseBankIn->fLower <10 || randIn->param.fLower < 10)
    {
      sprintf(msg, "--fl or --fl-signal or --fl-template must be >=10 Hz (%f %f)",
	      randIn->param.fLower , coarseBankIn->fLower);
      BEPrintError(msg);
    }

  if (userParam->startTime <= 0 )
    {
      sprintf(msg, "--gps-start-time (%d) invalid. must be > 0 ",
	      userParam->startTime );
      BEPrintError(msg);
    }

  if (userParam->maxTotalMass != -1)
    {
      if (userParam->maxTotalMass < 2*randIn->mMin) 
	{	
	  sprintf(msg, "--max-total-mass (%f) must be > twice the minimla mass (%f) ",
	      userParam->maxTotalMass , randIn->mMin );
	  BEPrintError(msg);
	}
    }

  /* massChoice setting (begin)*/


  if((userParam->template_m1 !=-1) && (userParam->template_m2!=-1))
    {
      if ((userParam->m1 == -1) || (userParam->m2==-1))
	{
	  sprintf(msg, "if --template-m1 and template-m2 are used, you must set --m1 and --m2 \n");
	  BEPrintError(msg);
	}       
    }


  if ((userParam->m1 != -1) && (userParam->m2 != -1)){

    randIn->param.massChoice = fixedMasses; 
    randIn->param.mass1 = userParam->m1;
    randIn->param.mass2 = userParam->m2;
    randIn->param.fCutoff = userParam->signalfFinal;
    /* we need to update the mMin which is used to estimate the length of the vetors*/
    if (userParam->m1  > userParam->m2){
      randIn->mMin = userParam->m2;
      randIn->mMax = userParam->m1+1e-2;
      
    }
    else{
      randIn->mMin = userParam->m1;
      randIn->mMax = userParam->m2+1e-2;
    }
    if (userParam->psi0 != -1 ||userParam->psi3 != -1 || userParam->tau0 != -1 || userParam->tau3 != -1){
      sprintf(msg, "--m1 --m2 --psi0 --psi3 --tau0 --tau3 error. If particular injection is requested,  you must choose either (--m1,--m2) options or (--psi0,--psi3) or (--tau0,--tau3)\n");
      BEPrintError(msg);
    } 
  }
     
  if (userParam->psi0 != -1 && userParam->psi3 != -1){
   randIn->param.massChoice = fixedPsi; 
   randIn->param.psi0 = userParam->psi0;
   randIn->param.psi3 = userParam->psi3;
   randIn->param.fCutoff = userParam->signalfFinal;

   if (userParam->m1 != -1 ||userParam->m2 != -1 || userParam->tau0 != -1 || userParam->tau3 != -1){
     sprintf(msg, "--m1 --m2 --psi0 --psi3 --tau0 --tau3 error. If particular injection is requested,  you must choose either (--m1,--m2) options or (--psi0,--psi3) or (--tau0,--tau3)\n");
     BEPrintError(msg);
   } 
  } 
  
 if (userParam->tau0 != -1 && userParam->tau3 != -1){
   randIn->param.massChoice = fixedTau; 
   randIn->param.t0 = userParam->tau0;
   randIn->param.t3 = userParam->tau3;
   randIn->param.fCutoff = userParam->signalfFinal;
   if (userParam->psi0 != -1 ||userParam->psi3 != -1 || userParam->m1 != -1 || userParam->m2 != -1){
     sprintf(msg, "--m1 --m2 --psi0 --psi3 --tau0 --tau3 error. If particular injection is requested,  you must choose either (--m1,--m2) options or (--psi0,--psi3) or (--tau0,--tau3)\n");
     BEPrintError(msg);  
   } 
  } 

 if  (userParam->binaryInjection == BHNS){	 
   randIn->param.massChoice = bhns;
 }

 /*massChoice setting (end)*/
 
 if (coarseBankIn->mmCoarse <=0 || coarseBankIn->mmCoarse>=1){
   sprintf(msg, "--mm (%f) must be in the range ]0 1[\n", 
	   coarseBankIn->mmCoarse);
     BEPrintError(msg);  
 } 
 if (userParam->noiseModel == None){
   sprintf(msg, "--noise-model must be <LIGOI, LIGOA, VIRGO, GEO, TAMA, REALPSD>\n");
   BEPrintError(msg);  
 } 
  
 if (coarseBankIn->numFcutTemplates <= 0){
   sprintf(msg, "--bank-number-fcut (%d) must be > 0>\n", coarseBankIn->numFcutTemplates);
   BEPrintError(msg);  
 } 
 
 if (userParam->run == None){
   sprintf(msg, "--run must be <S2 or S3>\n");
   BEPrintError(msg);  
 } 

 if  ((randIn->mMin >= randIn->mMax ) || (randIn->mMin <=0)){
    sprintf(msg, "--signal-mass-range (%f %f) paramter must be sorted and > 0 \n",
	    randIn->mMin, randIn->mMax);
    BEPrintError(msg);
  }
 else /* if valid parameters */
   {
     if (userParam->maxTotalMass != -1)
       {  
	 if (userParam->maxTotalMass >= (2 * randIn->mMax))
	   {
	  randIn->MMax = 2* randIn->mMax;
	  randIn->etaMin = randIn->mMin * randIn->mMax 
	    / (randIn->mMin + randIn->mMax)
	    / (randIn->mMin + randIn->mMax);
	   }
	 else
	   {
	     randIn->MMax = userParam->maxTotalMass;
	     randIn->etaMin = randIn->mMin * (randIn->MMax  - randIn->mMin)
	       / (randIn->MMax )
	       / (randIn->MMax );
	   }
       }          
     else 
       {
	 randIn->MMax = 2 * randIn->mMax;
	 randIn->etaMin = randIn->mMin * (randIn->MMax  - randIn->mMin)
	   / (randIn->MMax )
	   / (randIn->MMax );
       }
     randIn->param.fCutoff = userParam->signalfFinal;
   }

  if ((randIn->t0Min != 0 || randIn->t0Max != 0) && 
        ((randIn->t0Min >= randIn->t0Max ) || (randIn->t0Min <=0))){
    sprintf(msg, "--signal-tau0-range (%f %f) paramter must be sorted and > 0 \n",
            randIn->t0Min, randIn->t0Max);
    BEPrintError(msg);
  }

  if ((randIn->tnMin != 0 || randIn->tnMax != 0) && 
        ((randIn->tnMin >= randIn->tnMax ) || (randIn->tnMin <=0))){
    sprintf(msg, "--signal-tau3-range (%f %f) paramter must be sorted and > 0 \n",
            randIn->tnMin, randIn->tnMax);
    BEPrintError(msg);
  }

  switch ( userParam->detector)
    {
    case L1:
      switch (userParam->run)
	{
	case S2:
	  userParam->calCacheName   = userParam->L1.dataFile.S2.calCacheName  ;
	  userParam->frInCacheName  = userParam->L1.dataFile.S2.frInCacheName ;       
	  if (!userParam->chanName)
	    userParam->chanName  =  userParam->L1.chanName;

	  break;
	case S3:
	  userParam->calCacheName   = userParam->L1.dataFile.S3.calCacheName  ;
	  userParam->frInCacheName  = userParam->L1.dataFile.S3.frInCacheName ;
	  if (!userParam->chanName) userParam->chanName  =  userParam->L1.chanName;

	  break;
	case S1:
	case S4:
	  userParam->calCacheName   = userParam->L1.dataFile.S4.calCacheName  ;
	  userParam->frInCacheName  = userParam->L1.dataFile.S4.frInCacheName ;
	  if (!userParam->chanName) userParam->chanName  =  userParam->L1.chanName;
	  break;
	case S5:
	case S6:
	  break;
	}
      break;
    case H1:
      switch (userParam->run)
	{
	case S1:
	case S4:
	  userParam->calCacheName   = userParam->H1.dataFile.S4.calCacheName  ;
	  userParam->frInCacheName  = userParam->H1.dataFile.S4.frInCacheName ;
	 if (!userParam->chanName)  userParam->chanName  =  userParam->H1.chanName;
    	  break;
	case S5:
	case S6:
	  break;
	case S2:
	  userParam->calCacheName   = userParam->H1.dataFile.S2.calCacheName  ;
	  userParam->frInCacheName  = userParam->H1.dataFile.S2.frInCacheName ;
	 if (!userParam->chanName)  userParam->chanName  =  userParam->H1.chanName;
		  
	  break;
	case S3:
	  userParam->calCacheName   = userParam->H1.dataFile.S3.calCacheName  ;
	  userParam->frInCacheName  = userParam->H1.dataFile.S3.frInCacheName        ;
	 if (!userParam->chanName)  userParam->chanName       = userParam->H1.chanName;
	  break;
	}
      break;
     
    case H2:
      switch (userParam->run)
	{
	case S1:
	case S4:
	  userParam->calCacheName   = userParam->H2.dataFile.S4.calCacheName;
	  userParam->frInCacheName  = userParam->H2.dataFile.S4.frInCacheName;
	  if (!userParam->chanName) userParam->chanName  =  userParam->H2.chanName;
	  break;
	case S5:
	case S6:
	  break;
	case S2:
	  userParam->calCacheName   = userParam->H2.dataFile.S2.calCacheName;
	  userParam->frInCacheName  = userParam->H2.dataFile.S2.frInCacheName;
	  if (!userParam->chanName) userParam->chanName  =  userParam->H2.chanName;
	  break;
	case S3:
	  userParam->calCacheName   = userParam->H2.dataFile.S3.calCacheName;
	  userParam->frInCacheName  = userParam->H2.dataFile.S3.frInCacheName;
	  if (!userParam->chanName) userParam->chanName       =  userParam->H2.chanName;
	  break;
	}
      break;
      
    case V1:break;
      
    case G1:break;
    }
  if (userParam->faithfulness ==1 && randIn->type == 1)
      BEPrintError("can not check code if no injection performed. use simulation-type = 0 or 2\n");

  if (coarseBankIn->approximant == (Approximant)(-1))
      BEPrintError("--template,  template approximant must be provided\n");

  if (randIn->param.approximant == (Approximant)(-1))
      BEPrintError("--signal, signal approximant must be provided\n");

  if (userParam->binaryInjection == BHNS){
    if (randIn->mMin >3 || randIn->mMax < 3 ){
      BEPrintError("if you want to inject BHNS systems then adjust the mass-range so that the minimum is less than 3 solar mass and the maximum  is greater than 3 solar mass !! \n");
    }
  }      
}



void
BEPrintError(char *chaine)
{
  fprintf(stderr,"//--->  TestRandomBank code Error <---//\n");
  fprintf(stderr,"----->  %s",chaine);
  fprintf(stderr,"//---> please use \"--help\" option to get help documentation <---//\n");
  exit(0);
}



/* ****************************************************************************
 *  Documenation  on line
 *  **************************************************************************/
void Help()
{
  CHAR msg[2048];


  fprintf(stderr, "[NAME %s]\n ", CVS_NAME_STRING);
  fprintf(stderr, "[VERSION %s]\n ", CVS_ID_STRING);
  fprintf(stderr, "[VERSION %s]\n ", CVS_ID_STRING_C);
  fprintf(stderr, "[DESCRIPTION]\n");
  sprintf(msg, "\t lalapps_TestRandomBank is a standalone code testing the efficiency of\n"
	  "\t a *random* inpiral template bank - an algorithm implemented by I. Harry. \n"
	  "\t By efficiency we mean the match between the bank and an inspiral injection.\n"
	  "\t The code allows to use any design sensitivity curve provided  in lal/noisemodesl\n"
	  "\t or real PSD coming from GW data interferometer. (currently only L1, H1 and H2 \n"
	  "\t are handled but it is straightforwad to implement the GEO case. \n\n");
  fprintf(stderr, msg);
  sprintf(msg,
	  "\t The injections which can be performed uses the inspiral packages and therefore \n"  
	  "\t allows the following approximant [TaylorT1, TaylorT2, TaylorT3, EOB, PadeT1 \n"
	  "\t and SpinTaylorT3. \n\n\tThe bank and filtering uses the noisemodels packages for time- \n"
	  "\t domain approximant and BankEfficiency functions to perform BCV filtering method.\n"
	  "\t Both noisemodesl (for Time-domain) and BankEfficiency (for BCV) and in agreement\n"
	  );
  fprintf(stderr, msg);
  sprintf(msg,
	  "\t with findchirp routines.\n\n"
	  "\t Simulation can be done in absence of noise, in noise only or with an injection in noise.\n"
	  "\t In that code, we stored only the maximum point of the correaltion as our SNR. It is \n"
	  "\t sufficient in the case of overlap studies or when the data are gaussian and \n"
	  "\t injection have large snr.\n"
	  "\t Finally, results are stored in an xml files structures\n\n"
	  "\t SEE lalapps documenation for further descrption and examples.\n");
  fprintf(stderr, "%s\n",msg);

  
  fprintf(stderr, "[SYNOPSIS]\n");
  fprintf(stderr, "\t[--help]\n");
  fprintf(stderr, "\t[--verbose] \t\t\t gives some extra information on screen \n" );
  fprintf(stderr, "\t[--ascii2xml]\t\t read the file Trigger.dat, concatenation of one or several outputs from BankEfficiency and creates a unique XML outputs.\n");


  fprintf(stderr, "\t[--bank-alpha<float>]\t\t set the BCV alpha value in the moments computation\n");
  fprintf(stderr, "\t[--bank-fcut-range<float float>] set the range of BCV fcut (in units of GM) \n");
  fprintf(stderr, "\t[--bank-ffinal<float>]\t\t set the final frequency to be used in the BCV moments computation\n");
  fprintf(stderr, "\t[--bank-grid-spacing <gridSpacing>]\t set the grid type of the BCV bank (Square, SquareNotOriented, Hexagonal, HexagonalNotOriented\t\n");
  
  
  fprintf(stderr, "\t[--bank-number-fcut<integer>]\t set the number of BCV fcut \n");
  fprintf(stderr, "\t[--bank-mass-range<float float>] set the range of mass to be covered by the SPA bank\n");
  fprintf(stderr, "\t[--bank-psi0-range<float float>] set the range of psi0 to be covered by the BCV bank\n");
  fprintf(stderr, "\t[--bank-psi3-range<float float>] set the range of psi3 to be covered by the BCV bank\n");
  fprintf(stderr, "\t[--channel<string>]\t\t set the channel to look at \n");
  fprintf(stderr, "\t[--debug<integer>]\t\t set the debug level (same as in lal)\n");
  fprintf(stderr, "\t[--detector<string>]\t\t set the detector name to look at for real data (H1, H2, L1)\n");
  fprintf(stderr, "\t[--fl-signal<float>]\t\t set the lower cut off frequency of signal to inject\n");
  fprintf(stderr, "\t[--fl-template<float>]\t\t set the lower cut off frequnecy of template \n");
  fprintf(stderr, "\t[--fl<float>]\t\t\t set both template and signal lower cutoff frequency \n");
  fprintf(stderr, "\t[--max-total-mass<float>]\t set maximum total mass to be injected\n");
  fprintf(stderr, "\t[--gps-start-time<integer>]\t set gps start time if real data or psd are requested\n");
  fprintf(stderr, "\t[--m1<float>]\t\t\t force injection first individual mass to be equal to m1. needs to set m2 as well then\n");
  fprintf(stderr, "\t[--m2<float>]\t\t\t force injection second individual mass to be equal to m2. needs to set m1 as well then\n");
  fprintf(stderr, "\t[--mm<float>]\t\t\t set minimal match of the bank\n");
  fprintf(stderr, "\t[--n<float>]\t\t\t set number of trial in the simulation\n");
  fprintf(stderr, "\t[--ntrial<float>]\t\t same as --n\n");
  fprintf(stderr, "\t[--noise-amplitude<float>]\t set noise amplitude when using NoiseAndSignal flag simulation\n");
  fprintf(stderr, "\t[--noise-model<string>]\t\t set noise model curve to be <LIGOI, LIGOA, VIRGO, GEO, TAMA, REALPSD>\n");
  fprintf(stderr, "\t[--num-seconds<integer>]\t set number of seconds of data to look at.\n");
  fprintf(stderr, "\t[--psi0<float>]\t\t\t force injection psi0  value; request to psi3 as well. \n");
  fprintf(stderr, "\t[--psi3<float>]\t\t\t force injection psi3 value; request psi0 as well\n");
  fprintf(stderr, "\t[--run<string>]\t\t\t set run <S2, S3>\n");
  fprintf(stderr, "\t[--sampling<float>]\t\t set sampling frequency.\n");
  fprintf(stderr, "\t[--seed<integer>]\t\t set seed for random generator.\n");
  fprintf(stderr, "\t[--signal<string>]\t\t set signal approximant (TaylorT1, TaylorT2, TaylorT3, TaylorF2, PadeT1, EOB, SpinTaylorT3)\n");
  fprintf(stderr, "\t[--signal-alpha<float>]\t\t set alpha parameter of BCV injection\n");
  fprintf(stderr, "\t[--signal-amplitude<float>]\t set SNR of injection in the case NoiseandSignal simulation\n");
  fprintf(stderr, "\t[--signal-ffinal<float>]\t force final frequency value\n");
  fprintf(stderr, "\t[--signal-mass-range<float float>]\t set range of masses to inject (SPA injection)\n");
  fprintf(stderr, "\t[--signal-tau0-range<float float>]\t set range of tau0 to inject (SPA injection)\n");
  fprintf(stderr, "\t[--signal-tau3-range<float float>]\t set range of tau3 to inject (SPA injection)\n");
  fprintf(stderr, "\t[--signal-order<integer>]\t set PN order of injections \n");
  fprintf(stderr, "\t[--signal-psi0-range<float float>] set range of BCV injection \n");
  fprintf(stderr, "\t[--signal-psi3-range<float float>] set range of BCV injection\n");
  fprintf(stderr, "\t[--simulation-type<string>]\t set type of simulation (SignalOnly, noiseOnly, NoiseAndSignal)\n");
  fprintf(stderr, "\t[--tau0<float>]\t\t\t force injection to have tau0 value \n");
  fprintf(stderr, "\t[--tau3<float>]\t\t\t force injection to have tau3 value\n");
  fprintf(stderr, "\t[--template<string>]\t\tset signal approximant (TaylorT1, TaylorT2, TaylorT3, TaylorF2, PadeT1, EOB, SpinTaylorT3)\n");
  fprintf(stderr, "\t[--template-order<integer>]\t set PN order of template\n");
  fprintf(stderr, "\t[--alpha-constraint]\t\t set BCV code to be constrained \n");
  fprintf(stderr, "\t[--bhns-injection]\t\t set injection to be only bhbs systems\n");
  fprintf(stderr, "\t[--no-alpha-constraint]\t\t set BCV code to be unconstrained\n");
  fprintf(stderr, "\t[--faithfulness]\t check the code. template parameters are equal to injection parameters, size of the bank is therefore unity. It computed the faithfulness instead of effectualness\n");
  fprintf(stderr, "\t[--real-noise]\t\t\t use real data and real PSD.force simulaion type to be Noise Only\n");
  fprintf(stderr, "\t[--no-start-phase]\t\t\t unset random phase which is always set to zero.\n");
  fprintf(stderr, "\t[--print-psd]\t\t\t print the psd in  a file BE_PSD_type_gpstime.dat\n");
  fprintf(stderr, "\t[--print-best-overlap]\t\t print best overlap and other information\n");
  fprintf(stderr, "\t[--print-snr-histo]\t\t print histogram of the correlation output\n");
  fprintf(stderr, "\t[--print-bank]\t\t\t print the bank in ascii and xml format\n");
  fprintf(stderr, "\t[--print-prototype]\t\t print a prototype to be used by condor script\n");
  fprintf(stderr, "\t[--fast-simulation]\t\t perform fast simulation in the case of SPA abnk\n");
  fprintf(stderr, "type --print-default to get the default values of the current version\n");
  exit(0);
}


void 
BEAscii2Xml()
{
  UINT4 line=0;
  UINT8  id=0;
  UINT4 start = 0;
  ResultIn trigger;
  REAL4 tau0, tau3, tau0I, tau3I, psi0, psi3, phaseI, beta;
  FILE *input1, *input2, *bank;
  FILE *output;      
  SnglInspiralTable     *inputData = NULL;
  INT4 numFileTriggers = 0, nStartPad;






  char sbuf[512];
  
  /*  if (argc>1) {
    if (strcmp(argv[1], "--help")==0) {
      BEAscii2XmlHelp();
    }  
    if (strcmp(argv[1], "-h")==0) { 
      BEAscii2XmlHelp();
    }  
    }*/
  /* Main program starts here */
  /* First we open the file containing the ascii results */
  fprintf(stderr,"opening the xml output data file -- %s", BEASCII2XML_OUTPUT);
  output = fopen(BEASCII2XML_OUTPUT,"w");
  fprintf(stderr,"done\n");

  fprintf(stderr,"opening the xml input data file -- %s", BEASCII2XML_INPUT1);
  if  ( (input1  = fopen(BEASCII2XML_INPUT1,"r"))==NULL) {
    fprintf(stderr,"error while opening input file %s\n",BEASCII2XML_INPUT1);
    exit(0);
  }
  fprintf(stderr,"done\n");

  fprintf(stderr,"opening the xml prototype (argument of BankEfficiency code) -- ");
  if  ( (input2  = fopen(BEASCII2XML_INPUT2,"r"))==NULL)
    {
      fprintf(stderr,"error while opening input file %s\n",BEASCII2XML_INPUT2);
      fprintf(stderr,"the xml file will not contains parameters information\n");
      PRINT_LIGOLW_XML_HEADER(output);
      fprintf(stderr,"creating the header file -- done\n");
    }
  else 
    {
      /* read prototype and save in outputfile */
      fprintf(stderr,"parsing the prototype  -- ");
      while(fgets(sbuf,1024,input2) !=NULL)
	fputs(sbuf, output);
      fprintf(stderr," done\n");

      
    }
  

  /* insert the template bank here */
  if  ( (bank  = fopen(BEASCII2XML_BANK,"r"))==NULL)
    {
      fprintf(stderr,"error while opening input file %s\n",BEASCII2XML_BANK);
      fprintf(stderr,"the xml file will not contains the bank table\n");
    }
  else 
    {
      /* read prototype and save in outputfile */
      fprintf(stderr,"parsing the bank  -- ");
      fprintf( stdout, "reading triggers from file: %s\n", BEASCII2XML_BANK );
      numFileTriggers = 
        LALSnglInspiralTableFromLIGOLw( &inputData,BEASCII2XML_BANK , 0, -1 );
      fprintf(stderr," done %d\n", numFileTriggers);      
       myfprintf(output, LIGOLW_XML_SNGL_INSPIRAL );
       while(inputData)
	{
	  /*	  id = inputData->event_id->id;*/
	  
	  fprintf(output, SNGL_INSPIRAL_ROW, 
		  inputData->ifo,
		  inputData->search,
		  inputData->channel,
		  inputData->end_time.gpsSeconds,
		  inputData->end_time.gpsNanoSeconds,
		  inputData->end_time_gmst,
		  inputData->impulse_time.gpsSeconds,
		  inputData->impulse_time.gpsNanoSeconds,
		  inputData->template_duration,
		  inputData->event_duration,
		  inputData->amplitude,
		  inputData->eff_distance,
		  inputData->coa_phase,
		  inputData->mass1,
		  inputData->mass2,
		  inputData->mchirp,
		  inputData->mtotal,
		  inputData->eta,
		  inputData->tau0,
		  inputData->tau2,
		  inputData->tau3,
		  inputData->tau4,
		  inputData->tau5,
		  inputData->ttotal,
		  inputData->psi0,
		  inputData->psi3,
		  inputData->alpha,
		  inputData->alpha1,
		  inputData->alpha2,
		  inputData->alpha3,
		  inputData->alpha4,
		  inputData->alpha5,
		  inputData->alpha6,
		  inputData->beta,
		  inputData->f_final,
		  inputData->snr,
		  inputData->chisq,
		  inputData->chisq_dof,
		  inputData->sigmasq,
		  inputData->rsqveto_duration,
		  id);
	  inputData = inputData->next;
	  fprintf(output, "\n");

	}
       myfprintf(output, LIGOLW_XML_TABLE_FOOTER );

    }



  PRINT_LIGOLW_XML_BANKEFFICIENCY(output);
  fprintf(stderr,"done\n");
  /* read ascii input and save in xml format */
  fprintf(stderr,"reading the ascii file -- and saving xml file");
  do  
    {
      fscanf(input1,BANKEFFICIENCY_PARAMS_ROW_SPACE,
	     &trigger.psi0_trigger,
	     &trigger.psi3_trigger,
	     &trigger.beta_trigger,
	     &psi0, 
	     &psi3,
             &beta,   
	     &tau0, 
	     &tau3, 
	     &tau0I, 
	     &tau3I,
	     &trigger.fend_trigger, 
	     &trigger.fend_inject,
	     &trigger.mass1_inject,
	     &trigger.mass2_inject,
	     &phaseI,
	     &trigger.rho_final,
	     &trigger.snrAtCoaTime, 
	     &trigger.phase,
	     &trigger.alphaF, 
	     &trigger.bin, 
	     &nStartPad, 
	     &trigger.nfast); 

     if (start==0){
	      start+=1;
      }
      else 
      {
	      fprintf(output,",\n");
      }
      fprintf(output, BANKEFFICIENCY_PARAMS_ROW,
	     trigger.psi0_trigger,
	     trigger.psi3_trigger,
	     trigger.beta_trigger,
	     psi0, 
	     psi3,
             beta,   
	     tau0, 
	     tau3, 
	     tau0I, 
	     tau3I,
	     trigger.fend_trigger, 
	     trigger.fend_inject,
	     trigger.mass1_inject,
	     trigger.mass2_inject,
	     phaseI,
	     trigger.rho_final,
	     trigger.snrAtCoaTime, 
	     trigger.phase,
	     trigger.alphaF, 
	     trigger.bin, 
	     nStartPad, 
	     trigger.nfast); 
             line++;
    }
   while(!feof(input1));

    fprintf(stderr,"read %d lines...done\n", line);
  PRINT_LIGOLW_XML_TABLE_FOOTER(output);
  PRINT_LIGOLW_XML_FOOTER(output);

  fclose(output);
    fprintf(stderr,"closing xml file\n");


}


void BEAscii2XmlHelp(void)
{
  fprintf(stderr, "BEAscii2Xml help:\n");
  fprintf(stderr, "=================\n");
  /*  fprintf(stderr, "[PURPOSE]\n\t%s\n \t%s\n \t%s\n \t%s\n", 
		  "That code reads a file in ascii format generated by BankEfficiency",
		  "code and generate the appropriate xml files. It is used if one forget",
		  "to use the appropriate option within BankEfficiency code (--print-result-xml)",
		  "or when the bankEfficiency code has been used within the condor ",
		  "submit file. In that case indeed several ascii files have been ",
		  "generated but no xml files.\n");
  
  fprintf(stderr, "[INPUT/OUTPUT]\n\t%s\n \t%s\n \t%s\n \t%s %s %s %s\n", 
		  "That code do not need any input argument but needs an input file",
		  "in ascii format (the user has to be sure it is in the appropriate",
		  "format) and then it writes an xml file into an output file \n", 
		  "the input and output file name are ", BEASCII2XML_INPUT1, "and", BEASCII2XML_OUTPUT);
  */
  exit(1);
  

}
