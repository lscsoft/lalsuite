/* *****************************************************************************
 * Author: Cokelaer, T. and Sathyaprakash, B. S.
 * ************************************************************************** */
#include "BankEfficiency.h"

/* --- version information --- */
NRCSID( BANKEFFICIENCYC, "$Id$");
RCSID(  "$Id$");


/* --- the main program --- */
size_t nHisto = 200;
double minHisto = 0., maxHisto = 20.;

static REAL4Vector X1;
static REAL4Vector X2;
static REAL4Vector X3;
static REAL4Vector X4;


int
main (INT4 argc, CHAR **argv ) 
{
  /* --- Variables ---*/
  /*comments on variables will come later */
  UINT4		nn;
  UINT4         nby2;
  UINT4		currentTemplateNumber	= 0;
  UINT4		n;
  UINT4		k;
  INT4  	sizeBank = 0;
  UINT4  	kMin = 0;
  INT4  	temporder;
  UINT4         ntrials = 0;
  UINT4 	i;


  
  REAL8		df;
  REAL8		fendBCV	= 0;

  
  REAL4Vector   signal;
  REAL4Vector   correlation;
  REAL4Vector   VectorPhase;
  
  BEPowerVector  powerVector;
  BEMoments      moments; /* store the a11, a21, and a22 moments for any frequency values */

  REAL4Vector   Filter1; 
  REAL4Vector   Filter2;

  OverlapOutputIn OverlapOutputThisTemplate;
  OverlapOutputIn OverlapOutputBestTemplate;

  RandomInspiralSignalIn 	randIn;						/* random signal waveform to inject	*/

  InspiralWaveOverlapIn 	overlapin;					/* structure for Overlap		*/

  InspiralWaveOverlapOut     	overlapout; 

  ResultIn                      result;
  OtherParamIn 			otherIn;					/**/
  
  RealFFTPlan 			*fwdp = NULL;
  RealFFTPlan 			*revp = NULL;
  
  InspiralTemplateList      	*list = NULL;

  InspiralCoarseBankIn      	coarseBankIn; 					/* strcture for the bank of templates	*/

  LALStatus 		        status = blank_status;

  BCVMaximizationMatrix 	matrix;

  FILE                          *Foutput;


  float 			*bankEfficiencyOverlapIn;			
  float 			*bankEfficiencyOverlapInTime;			
  float 			*bankEfficiencyOverlapCIn;			
  float 			*bankEfficiencyOverlapCInTime;			



  
  gsl_histogram * histogramNoise = gsl_histogram_alloc (200);
  
  gsl_histogram_set_ranges_uniform (histogramNoise, 0.,20.);


  /* --- main code start here --- */ 

 
  /* --- Some initialization --- */
  /* Initialization of structure */
  ParametersInitialization(	&coarseBankIn, 
		  		&randIn,
			       	&otherIn);

  /* Read user parameters */
  ParseParameters(		&argc,
		  		argv,
			       	&coarseBankIn,
			       	&randIn,
			       	&otherIn);

  /* Check input parameters*/
  CheckParams(			coarseBankIn,
		 		randIn,
				otherIn);
  
  /* any debugging options ? */
  lalDebugLevel = otherIn.lalDebug ;

  /* This is used if one want to use condor script to keep track of the
     input parameters. It creates a prototype to be read by Ascii2xml */
  if (otherIn.PrintPrototype){
    BEPrintProtoXml(coarseBankIn, randIn, otherIn);    
    exit(0);
  }
  
  /* --- Estimate the size of the signal --- */
  LAL_CALL(BEGetMaximumSize(&status, randIn, &(signal.length)), 
	   &status);
  /* we multiply by 2 just for sanity */
  signal.length *= 2 ;
  signal.data 		= (REAL4*) LALCalloc(1, sizeof(REAL4) * signal.length);


  /* --- now that we have the maximum length of the signal we can go ahead --- */
  /* do we want to increase the size of the segments ? */
  signal.length *= otherIn.lengthFactor;

  /* --- Allocate memory for some vectors --- */
  correlation.length 	= signal.length;
  correlation.data 	= (REAL4*) LALCalloc(1, sizeof(REAL4) * correlation.length);

  randIn.psd.length 	= signal.length/2 + 1;
  randIn.psd.data 	= (REAL8*) LALMalloc(sizeof(REAL8) * randIn.psd.length); 


  /* create the PSd noise */
  LAL_CALL(BECreatePsd(&status, &coarseBankIn, randIn, otherIn), 
	   &status);

  /* copy the psd in RandIn.psd */
  for (i = 0; i< randIn.psd.length; i++){
    randIn.psd.data[i] = coarseBankIn.shf.data->data[i];
  }

  /* Create the bank od templates */
  if (otherIn.inputXMLBank){  
    BEReadXmlBank(&status, otherIn.inputXMLBank, &list, &sizeBank, coarseBankIn );
  }
  else{
    temporder  = coarseBankIn.order;
    coarseBankIn.order = 4; 			/* force to be 2PN otherwise the bank 
						   is not generated. has to be fixed 
						   in the bank itself*/
    LAL_CALL(LALInspiralCreateCoarseBank(&status, &list, &sizeBank, coarseBankIn),
	     &status);
    
    if (sizeBank == 0) {
      fprintf(stderr, "BankEfficiency Error :: bank is empty\n");
      exit(0);
    }
    
    for (i = 0; i < (UINT4)sizeBank; i++) {
      (list[i]).params.order = temporder;
    }
  }

  /* --- do we want to print the bank in asccii ? in xml format ?? --- */
  if (otherIn.PrintBank){
    BEPrintBank(coarseBankIn, &list, sizeBank);		      
  }
  if (otherIn.PrintBankXml){
    BEPrintBankXml(list,  sizeBank, coarseBankIn, randIn, otherIn);	
  }
  

  /* Finally let's  create extra storage */
  bankEfficiencyOverlapIn      = (float*) malloc(sizeof(float*)* sizeBank);              	/* allocate memory to store the bank results	*/
  bankEfficiencyOverlapInTime  = (float*) malloc(sizeof(float*)* sizeBank);              	/* allocate memory to store the bank results 	*/
  bankEfficiencyOverlapCIn     = (float*) malloc(sizeof(float*)* sizeBank);              	/* allocate memory to store the bank results	*/
  bankEfficiencyOverlapCInTime = (float*) malloc(sizeof(float*)* sizeBank);              	/* allocate memory to store the bank results 	*/
  
  
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
  
  list[1].params.startPhase = 0;   				/* why ? */
  
  /* --- Some new memory allocation --- */
  VectorPhase.length      = signal.length / 2;
  VectorPhase.data        = (REAL4*) LALCalloc(1, sizeof(REAL4) * VectorPhase.length);   
  
  moments.a11.length =  signal.length / 2 ;
  moments.a22.length =  signal.length / 2 ;
  moments.a21.length =  signal.length / 2 ;
  moments.a11.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) * moments.a11.length);
  moments.a22.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) * moments.a22.length);
  moments.a21.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) * moments.a21.length);

  powerVector.fm2_3.length = signal.length / 2 ;
  powerVector.fm5_3.length = signal.length / 2;       
  powerVector.fm7_6.length = signal.length / 2 ;
  powerVector.fm1_2.length = signal.length / 2;       
  powerVector.fm2_3.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) * powerVector.fm2_3.length);
  powerVector.fm5_3.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) * powerVector.fm5_3.length);
  powerVector.fm7_6.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) * powerVector.fm7_6.length);
  powerVector.fm1_2.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) * powerVector.fm1_2.length);

   
  Filter1.length          = correlation.length;
  Filter2.length          = correlation.length;
  Filter1.data            = (REAL4*) LALCalloc(1, sizeof(REAL4) * Filter1.length);
  Filter2.data            = (REAL4*) LALCalloc(1, sizeof(REAL4) * Filter2.length);
   
  
  n 	= signal.length / 2;
  df 	= randIn.param.tSampling / (REAL8)n / 2.;      
  
  /* --- Create useful vector once for all --- */
  LAL_CALL(LALCreateVectorFreqPower(&powerVector.fm5_3, randIn.param, -5, 3), 
	   &status);
  LAL_CALL(LALCreateVectorFreqPower(&powerVector.fm2_3, randIn.param, -2, 3), 
	   &status);
  LAL_CALL(LALCreateVectorFreqPower(&powerVector.fm7_6, randIn.param, -7, 6), 
	   &status);
  LAL_CALL(LALCreateVectorFreqPower(&powerVector.fm1_2, randIn.param, -1, 2), 
	   &status);

  /* --- Create Matrix for BCV Maximization once for all --- */
  LALCreateMomentVector(&moments, &coarseBankIn.shf,  &(list[0].params));

  X1.length = correlation.length;
  X2.length = correlation.length;
  X3.length = correlation.length;
  X4.length = correlation.length;


  X1.data = (REAL4*) LALMalloc(sizeof(REAL4) * X1.length);
  X2.data = (REAL4*) LALMalloc(sizeof(REAL4) * X2.length);
  X3.data = (REAL4*) LALMalloc(sizeof(REAL4) * X3.length);
  X4.data = (REAL4*) LALMalloc(sizeof(REAL4) * X4.length);  

  /*
   * Main loop is here. We creae random signal and filter 
   * it with the template bank.
   */
  /* while we want to perform one simulation */
  while (++ntrials <= otherIn.ntrials) 
    {
      /* do some initialisation for the output */
      BEInitOverlapOutputIn(&OverlapOutputBestTemplate);

      /* LALCreateRandomInjection*/
      randIn.param.approximant   = otherIn.signal; /* The waveform parameter for injection */
      /* however, we might need to cahnge some parameters depending on the user parameters */
      if (otherIn.signal == BCV) {
	randIn.param.massChoice = psi0Andpsi3;
	/* ending frequency will be fixed later */
      }
      else
	{ 
	  randIn.param.fCutoff 		= coarseBankIn.fUpper; 			      
          randIn.param.massChoice       = m1Andm2;
	  randIn.param.massChoice       = totalMassUAndEta; /* TODO : an input parameter*/
	}
      /* Let's compute the random parameters of the waveform to inject*/
      LAL_CALL(BEGenerateInputData(&status, &signal, &randIn, otherIn),
	       &status);


      overlapin.signal 	= signal;
      
      /* we might force to use only one template giving the input psi0 and psi3. 
	 Then, we don't need to use the bank*/
      /*if (otherIn.psi0 != -1 && otherIn.psi3 != -1)
	{
	  list[0].params.psi0 = otherIn.psi0;
	  list[0].params.psi3 = otherIn.psi3;
	  sizeBank  = 1;
	}
      */  
      
      nby2 	= Filter1.length / 2;
      nn 	= Filter1.length;
      kMin  	= floor(randIn.param.fLower /df );
      df   	= randIn.param.tSampling / (REAL8)(VectorPhase.length) / 2.;
      
      /* --- we can process the data through the bank now ---        */
      for (currentTemplateNumber = 0; currentTemplateNumber < (UINT4)sizeBank; currentTemplateNumber++) 
	{	
	  /* trick to check the code; template == signal and we use only one template 
	     suppose to be optimal if signal and template are equivalent. Obvisously not 
	     if signal family <> template family -- */
	  if (otherIn.check){
	    list[currentTemplateNumber].params = randIn.param;
	    overlapin.param                    = randIn.param;
	    sizeBank                           = 1;  /* we dont need the whole bank for checking. only one template needed */ 
	  }	     	     
	  else {
	    /* else we really go through the template parameter */
	    overlapin.param 	= list[currentTemplateNumber].params;	   
	    overlapout.max      = -1.;
	  }
	  
	  
	  /* if template==signal family then we can artificially increase code speed
	     by not performing match with different template. Can not be used with 
	     BCV and physical signal though since we do not know the correspondance 
	     between parameters. */
	  if (
	      (otherIn.FastSimulation == 1 ) &&
	      ((fabs(randIn.param.t0 - list[currentTemplateNumber].params.t0)> .1) ||
	       (fabs(randIn.param.t3 - list[currentTemplateNumber].params.t3)> .1) )
	      )
	    {
	      /*nothing dto do except plotting template paremter maybe */
	    }
	  else /*overlap is here below */
	    {	      
	      
	      
	      switch(otherIn.overlapMethod)
		{
		case AlphaMaximization:	   	   	   
		  /* let us fix the ending frequency of the template for LALOverlap function */
		  fendBCV = list[currentTemplateNumber].params.fFinal; 
		  overlapin.param.fFinal   =  fendBCV;
		  overlapin.param.fCutoff  =  fendBCV;
		  
		  /* Extract some values a11,a21,a22 parameter for template creation*/
		  k = floor(fendBCV / df);
		  BEGetMatrixFromVectors(moments, k,  &matrix);
		  
		  /* Template creation */
		  LAL_CALL(LALCreateFilters(&Filter1,
					    &Filter2,
					    &powerVector,
					    matrix,
					    kMin,					       
					    list[currentTemplateNumber].params.psi0,
					    list[currentTemplateNumber].params.psi3
					    ), &status);


		  /* The overlap given two filters and the input signal*/
		  /*BEInitOverlapOutputIn(&OverlapOutputThisTemplate); */
		  /*be sure rhomax = 0 and so on*/
		  LAL_CALL(LALWaveOverlapBCV(&status, 
					      &correlation, 
					      &overlapin,
					      &Filter1,
					      &Filter2, 
					      matrix,
					      otherIn, 
					      &OverlapOutputThisTemplate), 
			   &status);
		  
		  break;
		  /* this is obsolet*/
		case InQuadrature:
		  /* SHOULD be replace by flso of the template in the bank*/
		  
		  
		  fendBCV   = 1./LAL_PI/pow(6, 1.5)/(list[currentTemplateNumber].params.mass1+list[currentTemplateNumber].params.mass2)/LAL_MTSUN_SI;
		  
		  if (coarseBankIn.approximant==EOB && coarseBankIn.order==threePN)
		    fendBCV   = 1./LAL_PI/pow(2.3, 1.5)/(list[currentTemplateNumber].params.mass1+list[currentTemplateNumber].params.mass2)/LAL_MTSUN_SI;
		  
		  if (fendBCV > randIn.param.tSampling/2 ) 
		    fendBCV = randIn.param.tSampling/2. - 1;
		  
		  if (otherIn.check)
		    fendBCV   =  list[currentTemplateNumber].params.fFinal;  
		  
		  
		  /* should simply be 
		     fendBCV   =  list[currentTemplateNumber].params.fFinal;  
		     all the time isnt' it ? 
		  */
		  overlapin.param.fFinal  = fendBCV;
		  overlapin.param.fCutoff = fendBCV;
		  
		  LAL_CALL(LALInspiralWaveOverlap(&status,
						  &correlation,
						  &overlapout,
						  &overlapin), &status);
		  
		  BEFillOverlapOutput(overlapout, 
				      &OverlapOutputThisTemplate);
		  
		  break;	     
		}     /* switch end */

	      
	      /* fill historgranm of the correlation output */
	      if (otherIn.PrintSNRHisto){
		for (i=0; i<correlation.length; i++){
		  gsl_histogram_increment(histogramNoise, correlation.data[i]);
		}		    
	      }	


	      
	      /* --- if overlap is the largest one, then we keep some 
	       * other information . Later we might just use the previous 
	       * vector. Keep everything and use it outside the bank 
	       * process --- we should keep time as well*/	     
	      if (otherIn.alphaFConstraint==ALPHAFConstraint){
		bankEfficiencyOverlapCIn[currentTemplateNumber]     = OverlapOutputThisTemplate.rhoMaxConstraint; 
		bankEfficiencyOverlapCInTime[currentTemplateNumber] = OverlapOutputThisTemplate.rhoBinConstraint; 
		bankEfficiencyOverlapIn[currentTemplateNumber]     = OverlapOutputThisTemplate.rhoMaxUnconstraint; 
		bankEfficiencyOverlapInTime[currentTemplateNumber] = OverlapOutputThisTemplate.rhoBinUnconstraint; 
	 
	      }
	      else {
		bankEfficiencyOverlapIn[currentTemplateNumber]     = OverlapOutputThisTemplate.rhoMaxUnconstraint; 
		bankEfficiencyOverlapInTime[currentTemplateNumber] = OverlapOutputThisTemplate.rhoBinUnconstraint; 
	      }
		

	      
	      /* --- if overlap is the largest one, then we keep some 
	       * other information . Later we might just use the previous 
	       * vector. Keep everything and use it outside the bank 
	       * process --- */	     
	      OverlapOutputThisTemplate.freqConstraint                = fendBCV;
 	      OverlapOutputThisTemplate.freqUnconstraint              = fendBCV;
	      OverlapOutputThisTemplate.templateNumberConstraint      = currentTemplateNumber;
 	      OverlapOutputThisTemplate.templateNumberUnconstraint    = currentTemplateNumber;
	      OverlapOutputThisTemplate.layerConstraint               = list[currentTemplateNumber].nLayer;
	      OverlapOutputThisTemplate.layerUnconstraint             = list[currentTemplateNumber].nLayer;


	      KeepHighestValues(OverlapOutputThisTemplate, 
				&OverlapOutputBestTemplate);
	      
	    }/*end of  cheating code for non BCV code*/  
	  
	  
	}/*end of the bank process*/
      
     
      /* Then print the maximum overlap over the whole bank and the corresponding 
       * parameter of the templates, injected signal and so on. This is the main results
       * We'll print one line for each simulation (ntrials parameter)*/
      GetResult(	   &list, 
			   randIn.param, 
			   OverlapOutputBestTemplate, 
			   &result, 
			   otherIn);
      result.ntrial  = ntrials;
      result.coaTime = randIn.coalescenceTime;
	


         
      /* use the best template and keep the snr at coaTime */
#if 0
      fendBCV = list[OverlapOutputBestTemplate.templateNumberConstraint].params.fFinal;
      overlapin.param.fFinal   =  overlapin.param.fCutoff  =  fendBCV;
      k = floor(fendBCV / df);
      BEGetMatrixFromVectors(moments, k, &matrix);
      LAL_CALL(LALCreateFilters(&Filter1,
				&Filter2,
				&powerVector,
				matrix,
				kMin,					       
				list[OverlapOutputBestTemplate.templateNumberConstraint].params.psi0,
				list[OverlapOutputBestTemplate.templateNumberConstraint].params.psi3
				), &status);
      LAL_CALL(LALWaveOverlapBCV(&status, 
				 &correlation, &overlapin,
				 &Filter1, 
				 &Filter2,
				 matrix, 
				 otherIn , &OverlapOutputThisTemplate), &status);
      LAL_CALL(LALInspiralParameterCalc(&status, &list[OverlapOutputBestTemplate.templateNumberConstraint].params), &status);
      result.snrAtCoaTime = correlation.data[randIn.coalescenceTime];
      result.snrCAtCoaTime = correlation.data[randIn.coalescenceTime];
#endif      
      
      if (!otherIn.quietFlag){
	PrintResults(result);  
      }
      
      
      if (otherIn.PrintResultXml)
	{
	  BEPrintResultsXml(coarseBankIn,randIn,otherIn,result);	 
	}
      
 



      /*Now we might want given the best template, to reproduce results and keep overlp, wavefom of the signal , template ...*/
      /* save signal */
      /* overlap     */
      /* save template and correlation */
      if (otherIn.PrintBestOverlap || otherIn.PrintBestTemplate){
	
	otherIn.extraFinalPrinting = 1;
	

	Foutput=fopen("BE_waveformIn.dat", "w");
	for (i=0; i<signal.length; i++)
	  {
	    fprintf(Foutput, "%e %e\n", i/randIn.param.tSampling , signal.data[i]);
	  }
	fclose(Foutput);
	
	fendBCV = list[OverlapOutputBestTemplate.templateNumberConstraint].params.fFinal;
	overlapin.param.fFinal   =  overlapin.param.fCutoff  =  fendBCV;
	/* Extract some values for template creation*/
	k = floor(fendBCV / df);
	BEGetMatrixFromVectors(moments, k, &matrix);
	/* Template creation*/
	LAL_CALL(LALCreateFilters(&Filter1,
				  &Filter2,
				  &powerVector,
				  matrix,
				  kMin,					       
				  list[OverlapOutputBestTemplate.templateNumberConstraint].params.psi0,
				  list[OverlapOutputBestTemplate.templateNumberConstraint].params.psi3
				  ), &status);
	/* The overlap given two filters and the input signal*/
	LAL_CALL(LALWaveOverlapBCV(&status, 
				   &correlation, &overlapin,
				   &Filter1, 
				   &Filter2,
				   matrix, 
				   otherIn , &OverlapOutputThisTemplate), &status);
	
	LAL_CALL(LALInspiralParameterCalc(&status, &list[OverlapOutputBestTemplate.templateNumberConstraint].params), &status);
	/*	   if (!otherIn.quietFlag){
		   PrintResults(list[jmax].params, randIn.param, overlapoutmax, fMax, list[jmax].nLayer);*/
      }
      
    }  /*end while(trial)*/
  


  if (otherIn.PrintSNRHisto){
    Foutput=  fopen("BE_histo.dat","w");
    gsl_histogram_fprintf(Foutput, histogramNoise, "%f", "%g");
    fclose(Foutput);
  }


  
  /* --- destroy the plans, correlation and signal --- */
  
  LALFree(powerVector.fm5_3.data);
  LALFree(powerVector.fm2_3.data);
  LALFree(powerVector.fm1_2.data);
  LALFree(powerVector.fm7_6.data);

  LALFree(VectorPhase.data);

  LALFree(moments.a11.data);
  LALFree(moments.a21.data);
  LALFree(moments.a22.data);
   
  LALFree(Filter2.data);
  LALFree(Filter1.data);
  
  LALFree(randIn.psd.data);
  LALDDestroyVector( &status, &(coarseBankIn.shf.data) );
  
  LALFree(signal.data);
  LALFree(correlation.data);
  LALFree(list);
  
  LALDestroyRealFFTPlan(&status,&fwdp);
  LALDestroyRealFFTPlan(&status,&revp);
  
  gsl_histogram_free(histogramNoise);
   
  LALCheckMemoryLeaks();    
   
   
  return 0;
}



/* ****************************************************************************
 * Initialization of the parameters fromn random signal, bank and others 
 * structures
 * ***************************************************************************/
void ParametersInitialization(	InspiralCoarseBankIn 	*coarseBankIn, 
				RandomInspiralSignalIn	*randIn, 
				OtherParamIn		*otherIn)
{
  
  InitInspiralCoarseBankIn(coarseBankIn);
  InitRandomInspiralSignalIn(randIn);
  InitOtherParamIn(otherIn);
  
}


/* ****************************************************************************
 * CoarseIn initialization
 ***************************************************************************/
void InitInspiralCoarseBankIn(InspiralCoarseBankIn *coarseBankIn)
{

  coarseBankIn->fLower    	= BANKEFFICIENCY_FLOWER;
  coarseBankIn->fUpper    	= BANKEFFICIENCY_TSAMPLING/2. - 1.;
  coarseBankIn->tSampling 	= BANKEFFICIENCY_TSAMPLING;
  coarseBankIn->space     	= BANKEFFICIENCY_SPACE;
  coarseBankIn->mmCoarse  	= BANKEFFICIENCY_MMCOARSE;
  coarseBankIn->mmFine    	= BANKEFFICIENCY_MMFINE;
  coarseBankIn->iflso       = BANKEFFICIENCY_IFLSO ;
  coarseBankIn->mMin      	= BANKEFFICIENCY_MMIN;
  coarseBankIn->mMax      	= BANKEFFICIENCY_MMAX;
  coarseBankIn->MMax      	= BANKEFFICIENCY_MMAX * 2;
  coarseBankIn->massRange   = MinMaxComponentMass; 
  coarseBankIn->etamin 	= coarseBankIn->mMin * coarseBankIn->mMax  
			/ pow(coarseBankIn->MMax, 2.);

  coarseBankIn->psi0Min   	= BANKEFFICIENCY_PSI0MIN;
  coarseBankIn->psi0Max   	= BANKEFFICIENCY_PSI0MAX;
  coarseBankIn->psi3Min   	= BANKEFFICIENCY_PSI3MIN;
  coarseBankIn->psi3Max   	= BANKEFFICIENCY_PSI3MAX;
  coarseBankIn->alpha     	= BANKEFFICIENCY_ALPHABANK;
  /*unsigned variable ... */
  coarseBankIn->numFcutTemplates = BANKEFFICIENCY_NFCUT;
  coarseBankIn->approximant = BANKEFFICIENCY_TEMPLATE;  
  coarseBankIn->order       = BANKEFFICIENCY_ORDER_TEMPLATE;  
  coarseBankIn->LowGM     	= BANKEFFICIENCY_LOWGM;
  coarseBankIn->HighGM    	= BANKEFFICIENCY_HIGHGM;
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
  randIn->etaMin        = (BANKEFFICIENCY_MMIN * BANKEFFICIENCY_MMAX) 
    /BANKEFFICIENCY_MMAX/ BANKEFFICIENCY_MMAX;
  randIn->psi0Min  		= BANKEFFICIENCY_PSI0MIN;			/* psi0 range 				*/
  randIn->psi0Max  		= BANKEFFICIENCY_PSI0MAX;	
  randIn->psi3Min  		= BANKEFFICIENCY_PSI3MIN;	
  randIn->psi3Max  		= BANKEFFICIENCY_PSI3MAX;	
  randIn->param.approximant 	= BANKEFFICIENCY_SIGNAL;			/* approximant of h, the signal		*/
  randIn->param.tSampling 	= BANKEFFICIENCY_TSAMPLING;			/* sampling 				*/
  randIn->param.fCutoff		= BANKEFFICIENCY_TSAMPLING/2.-1;		/* sampling 				*/
  randIn->param.startTime       = BANKEFFICIENCY_STARTTIME; 
  randIn->param.startPhase      = BANKEFFICIENCY_STARTPHASE; 
  randIn->param.nStartPad       = BANKEFFICIENCY_NSTARTPHASE;
  randIn->param.signalAmplitude = BANKEFFICIENCY_SIGNALAMPLITUDE;
  randIn->param.nEndPad         = BANKEFFICIENCY_NENDPAD;
  randIn->NoiseAmp              = BANKEFFICIENCY_NOISEAMPLITUDE;
}


/* ****************************************************************************
 * Other Param initialization							
 * ***************************************************************************/
void InitOtherParamIn(OtherParamIn *otherIn)
{
  otherIn->alphaFConstraint      = ALPHAFConstraint;
  otherIn->extraFinalPrinting   =0; 
  otherIn->template		= BANKEFFICIENCY_TEMPLATE;
  otherIn->bank         	= -1;
  otherIn->lalDebug         	= 0;
  otherIn->signal       	= BANKEFFICIENCY_SIGNAL;
  otherIn->m1           	= -1;
  otherIn->m2           	= -1;
  otherIn->lengthFactor	        = 1;
  otherIn->psi0         	= -1;
  otherIn->psi3         	= -1;
  otherIn->tau0         	= -1;
  otherIn->tau3         	= -1;
  otherIn->PrintOverlap 	= BANKEFFICIENCY_PRINTOVERLAP;
  otherIn->PrintBestOverlap 	= BANKEFFICIENCY_PRINTBESTOVERLAP;
  otherIn->PrintBestTemplate 	= BANKEFFICIENCY_PRINTBESTTEMPLATE;
  otherIn->PrintSNRHisto = BANKEFFICIENCY_PRINTSNRHISTO;
  otherIn->PrintFilter  	= BANKEFFICIENCY_PRINTFILTER;
  otherIn->PrintPsd             = BANKEFFICIENCY_PRINTPSD;
  otherIn->overlapMethod	= AlphaMaximization;
  otherIn->snrAtCoaTime         = 0;

  otherIn->PrintBank    	= BANKEFFICIENCY_PRINTBANK;
  otherIn->PrintBankXml    	= BANKEFFICIENCY_PRINTBANKXML;
  otherIn->PrintResultXml    	= BANKEFFICIENCY_PRINTRESULTXML;
  otherIn->PrintPrototype    	= BANKEFFICIENCY_PRINTPROTOTYPE;

  otherIn->PrintBankOverlap	= BANKEFFICIENCY_PRINTBANKOVERLAP;
  otherIn->PrintTemplate	= BANKEFFICIENCY_PRINTTEMPLATE; 
  otherIn->check        	= BANKEFFICIENCY_CHECK;
  otherIn->inputPSD          	= NULL;
  otherIn->quietFlag 	     	= BANKEFFICIENCY_QUIETFLAG;
  otherIn->ntrials 	     	= BANKEFFICIENCY_NTRIALS;
  otherIn->FastSimulation       = BANKEFFICIENCY_FASTSIMULATION;
  otherIn->NoiseModel           = LIGOI;
  otherIn->binaryInjection=NoUserChoice; 
  otherIn->inputXMLBank = NULL;
}



/* ***************************************************************************
 *  Function to parse user parameters 
 *  ************************************************************************ */
void 
ParseParameters(	INT4 			*argc, 
			CHAR 			**argv,
			InspiralCoarseBankIn 	*coarseBankIn,
			RandomInspiralSignalIn 	*randIn,
			OtherParamIn    	*otherIn)
{
  INT4 		i = 1;


  while(i < *argc)
    {
      if ( strcmp(argv[i],	"--debug") 	== 0 ) {
	      otherIn->lalDebug = atoi(argv[++i]);
      }
      else if ( strcmp(argv[i],	"--fend-bcv") 	== 0 ) {
		coarseBankIn->LowGM      	= atof(argv[++i]);
       	      	coarseBankIn->HighGM     	= atof(argv[++i]);
	}
      else if ( (strcmp(argv[i],	"--fl") 	== 0 ) 
		|| ( strcmp(argv[i],	"--fLower") 	== 0 ) ){
  		coarseBankIn->fLower = randIn->param.fLower	= atof(argv[++i]);  
      }
      else if ( strcmp(argv[i],	"--sampling") 	== 0 ) {
  		coarseBankIn->tSampling = randIn->param.tSampling = atof(argv[++i]);  
		randIn->param.fCutoff 	= coarseBankIn->tSampling/2. - 1.;
		coarseBankIn->fUpper 	= coarseBankIn->tSampling/2. - 1.;

	}      
      else if ( strcmp(argv[i],	"--mass-range-bank")	== 0 ) 	{	
		coarseBankIn->mMin = atof(argv[++i]);
		coarseBankIn->mMax = atof(argv[++i]); 
		coarseBankIn->MMax = coarseBankIn->mMax * 2.; 
		coarseBankIn->etamin 	= coarseBankIn->mMin *  coarseBankIn->mMax 
		  / pow(coarseBankIn->MMax, 2.);
      }
      else if ( strcmp(argv[i],	"--mass-range")	== 0 ) 	{	
		randIn->mMin = atof(argv[++i]);
		randIn->param.mass1 = randIn->param.mass2 = randIn->mMin;
		randIn->mMax = atof(argv[++i]); 
		randIn->MMax = randIn->mMax * 2.; 
		randIn->etaMin 	= randIn->mMin * randIn->mMax 
			/ pow(randIn->MMax, 2.);
	}
      else if ( strcmp(argv[i], "--m1") 	== 0 ){
	otherIn->m1 = atof(argv[++i]);	      
	randIn->param.massChoice = fixedMasses; 
      }
      else if ( strcmp(argv[i], "--m2") 	== 0 ) { 	
	otherIn->m2 = atof(argv[++i]);
	randIn->param.massChoice = fixedMasses;  	       
      }
      else if ( strcmp(argv[i], "--psi0") 	== 0 ) {
	otherIn->psi0 = atof(argv[++i]);  
	randIn->param.massChoice = fixedPsi; 
	
      }
      else if ( strcmp(argv[i], "--psi3") 	== 0 ) {
	otherIn->psi3 = atof(argv[++i]); 
	randIn->param.massChoice = fixedPsi; 
	
      }
      else if ( strcmp(argv[i], "--tau0") 	== 0 ) {
	otherIn->psi0 = atof(argv[++i]);  
	randIn->param.massChoice = fixedTau; 	
      }
      else if ( strcmp(argv[i], "--tau3") 	== 0 ) {
	otherIn->psi3 = atof(argv[++i]); 
	randIn->param.massChoice = fixedTau; 
      }

      else if (strcmp(argv[i],	"--psi0-range")	==0) {
	coarseBankIn->psi0Min = randIn->psi0Min = atof(argv[++i]);
	coarseBankIn->psi0Max = randIn->psi0Max = atof(argv[++i]);
      }
      else if (strcmp(argv[i],	"--psi3-range")	==0){
		coarseBankIn->psi3Max = randIn->psi0Max = atof(argv[++i]);
		coarseBankIn->psi3Min = randIn->psi3Min = atof(argv[++i]);
	}
      else if ( strcmp(argv[i],	"--mm")		==0){
	     coarseBankIn->mmCoarse = atof(argv[++i]);
      }	  
      else if ( strcmp(argv[i],	"--number-fcut")	==0)  
  	  coarseBankIn->numFcutTemplates = atoi(argv[++i]);	
      else if ( strcmp(argv[i],	"--simulation-type")	==0)	     	
	      randIn->type = atoi(argv[++i]);	
      else if ( strcmp(argv[i],	"--signal-amplitude")	==0)	     
    	      randIn->SignalAmp = atof(argv[++i]);	
      else if ( strcmp(argv[i],	"--noise-amplitude")	==0)	     
    	      randIn->NoiseAmp = atof(argv[++i]);	
      else if ( strcmp(argv[i],	"--alpha-bank")	==0)	     
    	      coarseBankIn->alpha = atof(argv[++i]);      
      else if ( strcmp(argv[i],	"--alpha-signal")==0)	     
    	      randIn->param.alpha = atof(argv[++i]);
      else if ( strcmp(argv[i],	"--fFinal")==0)	     
    	      randIn->param.fCutoff = atof(argv[++i]);
      else if ( strcmp(argv[i],	"--freq-moment-bank")	==0)	     
	coarseBankIn->fUpper = atof(argv[++i]);      
      else if ( strcmp(argv[i],	"--template-order")==0)     
	      coarseBankIn->order = atoi(argv[++i]); 	
      else if ( strcmp(argv[i],	"--signal-order")==0) 
    	      randIn->param.order = atoi(argv[++i]); 	
      else if ( strcmp(argv[i],	"--ntrial")	==0)
	      otherIn->ntrials = atoi(argv[++i]);       	
      else if ( strcmp(argv[i],	"--n")	==0)
	      otherIn->ntrials = atoi(argv[++i]);       	
      else if ( strcmp(argv[i],	"--seed")	==0)
	      randIn->useed = atoi(argv[++i])+1;	
      else if ( strcmp(argv[i],"--lengthFactor") ==0)   
	      otherIn->lengthFactor = atoi(argv[++i]);	
      else if ( strcmp(argv[i],	"--noise-model")==0){
	  i++;

	  if (strcmp(argv[i], "LIGOI")		== 0)	otherIn->NoiseModel = LIGOI;
	  else if (strcmp(argv[i], "LIGOA") 	== 0)  	otherIn->NoiseModel = LIGOA;
	  else if (strcmp(argv[i], "VIRGO") 	== 0)  	otherIn->NoiseModel = VIRGO;
	  else if (strcmp(argv[i], "TAMA") 	== 0)   otherIn->NoiseModel = TAMA;
	  else if (strcmp(argv[i], "GEO") 	== 0)   otherIn->NoiseModel = GEO;
	  else if (strcmp(argv[i], "UNITY") 	== 0)   otherIn->NoiseModel = UNITY;
	  else 
	    {
	      otherIn->NoiseModel = REALPSD;
	      otherIn->filename= argv[i];
	      /*TODO check whether the file exist otherwise quit with error message */
	    }	  
	}
      else if ( strcmp(argv[i],"--alpha-constraint")	==0)     otherIn->alphaFConstraint	= ALPHAFConstraint;
      else if ( strcmp(argv[i],"--no-alpha-constraint")	==0)     otherIn->alphaFConstraint	= ALPHAFUnconstraint;
      else if ( strcmp(argv[i],"--quiet")		==0)     otherIn->quietFlag 		= 1;
      else if ( strcmp(argv[i],"--InQuadrature")	==0)     otherIn->overlapMethod 	= InQuadrature;
      else if ( strcmp(argv[i],"--print-overlap")	==0)  	 otherIn->PrintOverlap 		= 1;
      else if ( strcmp(argv[i],"--print-best-overlap")	==0)  	 otherIn->PrintBestOverlap 	= 1;
      else if ( strcmp(argv[i],"--print-best-template")	==0)  	 otherIn->PrintBestTemplate	= 1;
      else if ( strcmp(argv[i],"--snr-at-coa-time")	==0)  	 otherIn->snrAtCoaTime    	= 1;

      else if ( strcmp(argv[i],"--print-psd")	        ==0)  	 otherIn->PrintPsd 		= 1;
      else if ( strcmp(argv[i],"--PrintFilter")		==0) 	 otherIn->PrintFilter 		= 1; 
      else if ( strcmp(argv[i],"--print-bank-overlap")	==0)  	 otherIn->PrintBankOverlap 	= 1;
      else if ( strcmp(argv[i],"--print-snr-histo") ==0) otherIn->PrintSNRHisto 	= 1;
      else if ( strcmp(argv[i],"--print-template")	==0)  	 otherIn->PrintTemplate 	= 1;
      else if ( strcmp(argv[i],"--print-bank")		==0) 	 otherIn->PrintBank		= 1;
      else if ( strcmp(argv[i],"--print-bank-xml")	==0) 	 otherIn->PrintBankXml	        = 1;
      else if ( strcmp(argv[i],"--print-result-xml")	==0) 	 otherIn->PrintResultXml        = 1;
      else if ( strcmp(argv[i],"--print-prototype")	==0) 	 otherIn->PrintPrototype        = 1;
      else if ( strcmp(argv[i],"--read-bank-xml")       ==0){
	otherIn->inputXMLBank = argv[++i];
      }
      else if ( strcmp(argv[i],"--fast-simulation")      ==0)     otherIn->FastSimulation        = 1;
      else if ( strcmp(argv[i],"--binaryInjection")      ==0) {
	  if ( strcmp(argv[++i],	"BHNS")	==0)   {
	    otherIn->binaryInjection = BHNS;
	    /*should force randIn.MMax to be 3+ maxmass and not 2timesmaxmass. But
	     mass-range option might be given after that option. so we should do that later.
	    */
	  }
	  else {
	    fprintf(stderr, "binaryInjection should be follow by [BHNS]\n");
	    Help(*coarseBankIn, *randIn, *otherIn);
	    exit(0);
	  }
	  
      }


      else if ( strcmp(argv[i],"--check")		==0)     otherIn->check 		= 1;
      else if ( strcmp(argv[i],"--InputPSD")		==0)  	 otherIn->inputPSD 		= argv[++i];
      else if ( strcmp(argv[i],"--template")		==0)
	{
	  if ( strcmp(argv[++i],	"TaylorT1")	==0)	otherIn->template = TaylorT1;
	  else if ( strcmp(argv[i],	"TaylorT2")	==0)	otherIn->template = TaylorT2;
	  else if ( strcmp(argv[i],	"TaylorT3")	==0)	otherIn->template = TaylorT3;
	  else if ( strcmp(argv[i],	"TaylorF1")	==0)	otherIn->template = TaylorF1;
	  else if ( strcmp(argv[i],	"TaylorF2")	==0)	otherIn->template = TaylorF2;
	  else if ( strcmp(argv[i],	"PadeT1")	==0)	otherIn->template = PadeT1;
	  else if ( strcmp(argv[i],	"PadeF1")	==0)	otherIn->template = PadeF1;
	  else if ( strcmp(argv[i],	"EOB")		==0)	otherIn->template = EOB;
	  else if ( strcmp(argv[i],	"BCV")		==0)    otherIn->template = BCV;
	  else if ( strcmp(argv[i],	"SpinTaylorT3")	==0)	otherIn->template = SpinTaylorT3;

	  coarseBankIn->approximant = otherIn->template;
	  if ( coarseBankIn->approximant == BCV ) 	
	    coarseBankIn->space = Psi0Psi3; 
	  else 
	    coarseBankIn->space = Tau0Tau3;
	}	
      else if ( (strcmp(argv[i], "--signal")==0))
	{
	  if (strcmp(argv[++i],	  "TaylorT1")	==0) 		otherIn->signal = TaylorT1;
	  else if (strcmp(argv[i],"TaylorT2")	==0)		otherIn->signal = TaylorT2;
	  else if (strcmp(argv[i],"TaylorT3")	==0)	    	otherIn->signal = TaylorT3;
	  else if (strcmp(argv[i],"TaylorF1")	==0)	    	otherIn->signal = TaylorF1;
	  else if (strcmp(argv[i],"TaylorF2")	==0)	    	otherIn->signal = TaylorF2;
	  else if (strcmp(argv[i],"PadeT1")	==0)	    	otherIn->signal = PadeT1;
	  else if (strcmp(argv[i],"PadeF1")	==0)	    	otherIn->signal = PadeF1;
	  else if (strcmp(argv[i],"EOB")	==0)    	otherIn->signal = EOB;
	  else if (strcmp(argv[i],"BCV")	==0)            otherIn->signal = BCV;
	  else if (strcmp(argv[i],"SpinTaylorT3")==0)	    	otherIn->signal = SpinTaylorT3;
	  randIn->param.approximant = otherIn->signal;	
	  
	}	
      else 
	{
	  Help(*coarseBankIn, *randIn, *otherIn);
	  exit(0);
	}
      
      i++;       
    } 
}




/* --- Function to check validity of some parsed parameters (TODO)--- */
void CheckParams(InspiralCoarseBankIn coarseBankIn,
		 RandomInspiralSignalIn randIn,
		 OtherParamIn otherIn)
{
  REAL4 temp;
  temp =randIn.param.mass1; /*just to avoid boring warning*/

  if (coarseBankIn.psi0Min <=0 ) {
    BEPrintError("Psi0 must be > 0");
    exit(0);
  }
  if (coarseBankIn.psi0Max <=0 ) {
    BEPrintError("psi0 Max must be > 0"); 
	  exit(0);
  }
  if (coarseBankIn.psi3Max >=0 ) {
    BEPrintError("psi3 Max must be < 0"); 
    exit(0);
  }
  if (coarseBankIn.psi3Min >=0 ) {
    BEPrintError("psi3 Min must be < 0"); 
    exit(0);
  }
  if (coarseBankIn.psi0Min > coarseBankIn.psi0Max
      || coarseBankIn.psi3Min > coarseBankIn.psi3Max){
    BEPrintError("psi range should be [psiMin psiMax]; take care of the sign. (i.e. -10 -2000 or 10 2000)\n");
    exit(0);
  }
  if (coarseBankIn.mMin > coarseBankIn.mMax )
    {
      BEPrintError("in option --mass-range, first argument shoukd be < to the second one\n");
    }
  if (coarseBankIn.approximant != BCV && otherIn.overlapMethod != InQuadrature)
    {
      BEPrintError("If the template are not BCV, the overlap \nmethod must be InQuadrature (use the \"--InQuadrature\" option)\n");
      exit(0);
    }
  if (coarseBankIn.tSampling <= 2.*coarseBankIn.fUpper)
    {
      BEPrintError("arg given by  --freq-moment-bank should be < to half the value of the argument given by --sampling \n");
      exit(0);
    }
  if (otherIn.binaryInjection==BHNS){
    if (randIn.mMin >3 || randIn.mMax < 3 ){
      BEPrintError("if you want to inject BHNS systems then adjust the mass-range so that the minimum is less than 3 solar mass and the maximum  is greater than 3solar mass !! \n");
      exit(0);
    }
  }
  
  /* check if we want a well defined injection. If so input paremter should be consistent. */
  if (otherIn.m1 != -1 || otherIn.m2 != -1 ){
     if (otherIn.psi0 != -1 ||otherIn.psi3 != -1 || otherIn.tau0 != -1 || otherIn.tau3 != -1){
       fprintf(stderr, "ERROR::The choice of arbitrary injection is exclusive !! you should choose either (--m1,--m2) options or (--psi0,--psi3) or (--tau0,--tau3)\n");
       exit(0);       
     } 
     if (otherIn.m1 == -1 || otherIn.m2 == -1 ){       
       fprintf(stderr, "ERROR::If you define m1 you must defined m2 as well and vice versa.\n");
       exit(0);       
     }
     
  } 
  if (otherIn.psi0 != -1 || otherIn.psi3 != -1 ){
    if (otherIn.m1 != -1 ||otherIn.m2 != -1 || otherIn.tau0 != -1 || otherIn.tau3 != -1){
      fprintf(stderr, "ERROR::The choice of arbitrary injection is exclusive !! you should choose either (--m1,--m2) options or (--psi0,--psi3) or (--tau0,--tau3)\n");
      exit(0);       
    } 
    if (otherIn.psi0 == -1 || otherIn.psi3 == -1 ){       
      fprintf(stderr, "ERROR::If you define psi3 you must defined psi0 as well and vice versa.\n");
      exit(0);       
    }
    
  } 
  if (otherIn.tau0 != -1 || otherIn.tau3 != -1 ){
    if (otherIn.psi0 != -1 ||otherIn.psi3 != -1 || otherIn.tau0 != -1 || otherIn.tau3 != -1){
      fprintf(stderr, "ERROR::The choice of arbitrary injection is exclusive !! you should choose either (--m1,--m2) options or (--psi0,--psi3) or (--tau0,--tau3)\n");
      exit(0);       
    } 
    if (otherIn.tau0 == -1 || otherIn.tau3 == -1 ){       
      fprintf(stderr, "ERROR::If you define tau0 you must defined tau3 as well and vice versa.\n");
      exit(0);       
    }
  } 
  
}


void
BEPrintError(char *chaine)
{
  fprintf(stderr,"|=============================================================|\n");
  fprintf(stderr,"| BankEfficiency code						\n");
  fprintf(stderr,"| Error while parsing parameters				\n");
  fprintf(stderr,"|=============================================================|\n");
  fprintf(stderr,"| %s \n",chaine);
  fprintf(stderr,"|type BankEfficiency -h to get details\n");
  fprintf(stderr,"|=============================================================|\n");
}



/* ****************************************************************************
 *  Documenation  on line
 *  **************************************************************************/
void Help(	InspiralCoarseBankIn   coarseBankIn,
       		RandomInspiralSignalIn randIn,
		OtherParamIn           otherIn)
{
  INT4 		temp;
	
  fprintf(stderr,"-------------------\n");
  fprintf(stderr,"BankEfficiency Code\n");
  fprintf(stderr,"-------------------\n");
  fprintf(stderr,"\nPURPOSE: Test efficiency of bank of templates using Monte-Carlo simulations\n");

  fprintf(stderr,"\nUSAGE:  [options]\n");
  fprintf(stderr,"There are two categories of options (default values in brackets)\n");
  fprintf(stderr,"First options must be followed by one or two arguments:\n\n");
  fprintf(stderr,"--alpha-bank		: BCV amplitude correction parameter 			(%7.2f)\n", BANKEFFICIENCY_ALPHABANK);
  fprintf(stderr,"--alpha-signal	: BCV amplitude correction parameter 			(%7.2f)\n",BANKEFFICIENCY_ALPHASIGNAL);
  fprintf(stderr,"--binaryInjection[BHNS]: if argument=BHNS then system are Blackhole+neutron star\n");
  fprintf(stderr,"--debug		: debugging option for lal functions 0=no edbug, 1=debug ... (0)\n");
  fprintf(stderr,"--fast-simulation	:fast simulation in the SPA case in a box of 0.1seconds square\n");
  fprintf(stderr,"--fLower			: lower frequency cutoff             		(%7.2f) Hz\n", BANKEFFICIENCY_FLOWER);
  fprintf(stderr,"--fl			: lower frequency cutoff             			(%7.2f) Hz\n", BANKEFFICIENCY_FLOWER);
  fprintf(stderr,"--length-factor	: multiply length signal by some factor(should be of unity order) ( %7.2d)\n", 1);
  fprintf(stderr,"--fend-bcv		: Lower and highest values of bcv frequency cutoff 	(%7.2d, %7.2d)\n", BANKEFFICIENCY_LOWGM, BANKEFFICIENCY_HIGHGM);
  fprintf(stderr,"--mass-range		: minimal mass of component stars    			(%7.2f, %7.2f) Mo\n", BANKEFFICIENCY_MMIN, BANKEFFICIENCY_MMAX);
  fprintf(stderr,"--mm			: minimal match for template bank    			(%7.3f)\n", BANKEFFICIENCY_MMCOARSE);
  fprintf(stderr,"--m1 <> --m2 <>	: Set parameter values of the injected signal \n");
  fprintf(stderr,"--tau0 <> --tau3 <>	: which should be one othe three parameters (maases, psi \n");
  fprintf(stderr,"--psi0 <> --psi3 <>	: tau components). That choice is exclusive. You cant set m1 and then tau0\n");
  fprintf(stderr,"--noise-amplitude	: use it with simulation-type <> 0 			(%7.2f)\n", BANKEFFICIENCY_NOISEAMPLITUDE); 
  fprintf(stderr,"--noise-model		: design sensitivity to use(LIGOI ,LIGOA VIRGO, GEO, TAMA	(%s)\n", BANKEFFICIENCY_NOISEMODEL);
  fprintf(stderr,"--number-fcut		: number of layers in Fcut dimension 			(%7.2d)\n", BANKEFFICIENCY_NFCUT);
  fprintf(stderr,"--ntrial		: number of trials                   			(%7.2d)\n", BANKEFFICIENCY_NTRIALS);
  fprintf(stderr,"--psi0-range		: psi0 range in the bank				(%7.2f, %7.2f) \n", BANKEFFICIENCY_PSI0MIN, BANKEFFICIENCY_PSI0MAX);
  fprintf(stderr,"--psi3-range		: psi3 range in the bank 				(%7.2f, %7.2f) \n", BANKEFFICIENCY_PSI3MAX, BANKEFFICIENCY_PSI3MIN);
  fprintf(stderr,"--template		: PN model for template bank, e.g. BCV 			(%7.2d)\n", BANKEFFICIENCY_TEMPLATE);
  fprintf(stderr,"--sampling		: sampling frequency             			(%7.2f)\n", BANKEFFICIENCY_TSAMPLING);
  fprintf(stderr,"--seed			: seed for random generation         		(%7.2d)\n", BANKEFFICIENCY_USEED);
  fprintf(stderr,"--signal		: model to inject in Monte-Carlo, e.g. EOB            	(%7.2d)\n", BANKEFFICIENCY_SIGNAL);
  fprintf(stderr,"--signal-order		: order of PN model                             (%7.2d)\n", BANKEFFICIENCY_ORDER_SIGNAL);
  fprintf(stderr,"--simulation-type	: type of simulation, 0, 1 or 2      			(%7.2d)\n", BANKEFFICIENCY_TYPE);
  fprintf(stderr,"--template-order	: order of signal to injec                            	(%7.2d)\n", BANKEFFICIENCY_ORDER_TEMPLATE);
  fprintf(stderr,"--signal-amplitude	: amplitude of the signal                             	(%7.2f)\n", BANKEFFICIENCY_SIGNALAMPLITUDE);

  fprintf(stderr,"\n --\nThe second category doesn't take any arguments.\n Those options are mainly used to print results in files\n");
  fprintf(stderr,"Verbose Options \n\n");
  fprintf(stderr,"     --check		: For a given random waveform, we compute overlap with same parameter (should get overlap=1)\n");
  /*fprintf(stderr,"     --FMaximization	: Maximize over Ending frequency 			(%d)\n", BANKEFFICIENCY_FMAXIMIZATION);*/
  fprintf(stderr,"     --FastSimulation		: filters only template around the injection (dt0=.5 and dt3=.25)\n");

  fprintf(stderr,"     --PrintOverlap	: Print Overlap given by the best template 		(%d)\n", BANKEFFICIENCY_PRINTOVERLAP);
  /*  fprintf(stderr,"     --PrintFilter	: Print filters giving the best overlap template 	(%d)\n", BANKEFFICIENCY_PRINTFILTER);
  fprintf(stderr,"     --print-ambiguity-function : print ambiguity function (bank overlap)	(%d)\n", BANKEFFICIENCY_AMBIGUITYFUNCTION);*/
  fprintf(stderr,"     --print-bank-overlap    : Print the overlap given by each template (default=no printing)		\n");

  fprintf(stderr,"     --print-bank	: Print the bank in ascii format (%s)		\n",  BANKEFFICIENCY_PRINTBANK_FILEASCII);
  fprintf(stderr,"     --print-snr-histo  : Print histogram of all the snr output (%d)	\n",  BANKEFFICIENCY_PRINTSNRHISTO);

  fprintf(stderr,"     --print-bank-xml	: Print the bank in xml format (%s)		\n",  BANKEFFICIENCY_PRINTBANK_FILEXML);
  fprintf(stderr,"    --print-result-xml  : Print result in xml format (%s)		\n",  BANKEFFICIENCY_PRINTRESULT_FILEXML);


  fprintf(stderr,"     --print-overlap	: Print the overlap given by the best template	(%s)	(%d)\n", 
		  BANKEFFICIENCY_PRINTOVERLAP_FILE, 
		  BANKEFFICIENCY_PRINTOVERLAP);
  fprintf(stderr,"     --print-psd	: Print the psd in                              	(%s)	(%d)\n", 
		  BANKEFFICIENCY_PRINTPSD_FILE, 
		  BANKEFFICIENCY_PRINTPSD);

  fprintf(stderr,"     --InQuadrature	: Overlap is compute without alpha maximization\n");
  fprintf(stderr,"     --quiet		: if this flag is present, the output is restricted to the minimum\n");
/*  fprintf(stderr,"     --InputPSD	: to give the name of real PSD\n");*/
  /*to avoid boring warning*/
  temp = coarseBankIn.approximant;
  temp =  randIn.useed;
  temp =  otherIn.signal;
}


		      
/* ****************************************************************************
 * Some output Results
 * ************************************************************************** */
void
KeepHighestValues(OverlapOutputIn resultThisTemplate,   
		  OverlapOutputIn *resultBestTemplate
		  )
{
  if (resultThisTemplate.rhoMaxConstraint > resultBestTemplate->rhoMaxConstraint){    
      resultBestTemplate->rhoMaxConstraint   = resultThisTemplate.rhoMaxConstraint;
      resultBestTemplate->phaseConstraint    = resultThisTemplate.phaseConstraint;
      resultBestTemplate->alphaConstraint    = resultThisTemplate.alphaConstraint;
      resultBestTemplate->rhoBinConstraint   = resultThisTemplate.rhoBinConstraint;
      resultBestTemplate->freqConstraint     = resultThisTemplate.freqConstraint;
      resultBestTemplate->layerConstraint    = resultThisTemplate.layerConstraint;
      resultBestTemplate->templateNumberConstraint   = resultThisTemplate.templateNumberConstraint;
      /*we do not need all template params for the time being*/
  }
  
  
  if (resultThisTemplate.rhoMaxUnconstraint > resultBestTemplate->rhoMaxUnconstraint){ 
    resultBestTemplate->rhoMaxUnconstraint   = resultThisTemplate.rhoMaxUnconstraint;
    resultBestTemplate->phaseUnconstraint    = resultThisTemplate.phaseUnconstraint;
    resultBestTemplate->alphaUnconstraint   = resultThisTemplate.alphaUnconstraint;
    resultBestTemplate->rhoBinUnconstraint   = resultThisTemplate.rhoBinUnconstraint;
    resultBestTemplate->freqUnconstraint     = resultThisTemplate.freqUnconstraint;
    resultBestTemplate->layerUnconstraint    = resultThisTemplate.layerUnconstraint;
    resultBestTemplate->templateNumberUnconstraint = resultThisTemplate.templateNumberUnconstraint;
  }
  
}


/* ****************************************************************************
 * Some output Results
 * ************************************************************************** */

void 
GetResult(
	  
	  InspiralTemplateList        **list,
	  InspiralTemplate       	injected,
	  OverlapOutputIn 	        bestOverlap, 
	  ResultIn                      *result,
	  OtherParamIn                  otherIn )
{
  INT4 templateNumber, templateNumberC;
  InspiralTemplate trigger, triggerC;
  LALStatus status = blank_status;
  
  /* INITSTATUS (status, "GetResult", BANKEFFICIENCYC);
     ATTATCHSTATUSPTR(status);
  */
    templateNumberC = bestOverlap.templateNumberConstraint;
    templateNumber = bestOverlap.templateNumberUnconstraint;
    trigger = (*list)[templateNumber].params;
    triggerC = (*list)[templateNumberC].params;
    LALInspiralParameterCalc(&status, &trigger);
    LALInspiralParameterCalc(&status, &triggerC);


  if (otherIn.overlapMethod != InQuadrature){
    result->psi0_trigger = trigger.psi0;
    result->psi3_trigger = trigger.psi3;
    result->psi0_inject  = injected.psi0;
    result->psi3_inject  = injected.psi3;     
    result->psi0_triggerC  = triggerC.psi0;
    result->psi3_triggerC  = triggerC.psi3;     
  }
  else{



    result->psi0_trigger = trigger.t0;
    result->psi3_trigger = trigger.t3;
    result->psi0_inject  = injected.t0;
    result->psi3_inject  = injected.t3; 
    result->psi0_triggerC  = -1;
    result->psi3_triggerC  = -1;     
  }

  result->mass1_inject      = injected.mass1;
  result->mass2_inject      = injected.mass2;
  result->fend_inject       = injected.fFinal;
  result->fend_trigger      = trigger.fFinal;
  result->fend_triggerC      = triggerC.fFinal;
  result->totalMass_trigger = trigger.totalMass;
  result->eta_trigger       = trigger.eta;
  result->totalMass_triggerC = triggerC.totalMass;
  result->eta_triggerC       = triggerC.eta;

  result->rho_finalC   = bestOverlap.rhoMaxConstraint;
  result->alphaC       = bestOverlap.alphaConstraint;
  result->alpha_fC     = bestOverlap.alphaConstraint*pow(triggerC.fFinal,2./3.);
  result->binC         = bestOverlap.rhoBinConstraint;
  result->phaseC       = bestOverlap.phaseConstraint;
  result->layerC       = bestOverlap.layerConstraint;
  
  result->rho_final   = bestOverlap.rhoMaxUnconstraint;
  result->alpha       = bestOverlap.alphaUnconstraint;
  result->alpha_f     = bestOverlap.alphaUnconstraint*pow(trigger.fFinal,2./3.);
  result->bin         = bestOverlap.rhoBinUnconstraint;
  result->phase       = bestOverlap.phaseUnconstraint;
  result->layer       = bestOverlap.layerUnconstraint;
  
    
  /* DETATCHSTATUSPTR(status);
  RETURN (status);
  */
}


void
PrintResults( ResultIn result)
{

  fprintf(stdout, "%e %e %e %e %e %e  %7.2f %7.2f %7.2f  %e %e  %e %e %e %e    %7.5f %e %e %e  %d %d    %7.5f %e %e %e  %d %d %d %e %e\n",
	  result.psi0_trigger, 
	  result.psi3_trigger,
	  result.psi0_triggerC, 
	  result.psi3_triggerC,
	  result.psi0_inject,
	  result.psi3_inject,
	  result.fend_trigger, 
	  result.fend_triggerC, 
	  result.fend_inject,
	  result.totalMass_trigger,
	  result.eta_trigger,
	  result.totalMass_triggerC,
	  result.eta_triggerC,
	  result.mass1_inject,
	  result.mass2_inject,
	  result.rho_final, 
	  result.phase, 
	  result.alpha,
	  result.alpha_f,
	  result.layer,
	  result.bin,
	  result.rho_finalC, 
	  result.phaseC, 
	  result.alphaC,
	  result.alpha_fC,
	  result.layerC,
	  result.binC, 
	  result.coaTime, 
	  result.snrAtCoaTime, 
	  result.snrCAtCoaTime);
  fflush(stdout);

}





/* ****************************************************************************
 * Functio to generate and stored the moments in  REAL4vectors. 
 * ************************************************************************* */
void LALCreateMomentVector(BEMoments             *moments,
			   REAL8FrequencySeries  *psd,
			   InspiralTemplate      *params)
{

  REAL8 		m7 = 0.;							/* the moments */
  REAL8 		m5 = 0.;	
  REAL8 		m3 = 0.;
  REAL8 		f;
	  
  UINT4			kMin;
  UINT4		  	kMax;
  UINT4			k;

  InspiralMomentsIn 	in;
    
  /* inspiral structure */
  in.shf 	= psd;  							/* The spectrum 			*/
  in.xmin 	= params->fLower;						/* Lower frequency of integral?	is it correct or should we put zero? 	*/
  in.xmax 	= params->tSampling/2.;						/* We dont' need to cut here		*/
  in.norm 	= 1./4.;							/* Some normalization			*/


  in.norm *= params->tSampling * params->tSampling;

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
void LALGetOrthogonalFilter2(REAL4Vector *filter)
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
LALCreateFilters(	REAL4Vector 		*Filter1,
		 	REAL4Vector 		*Filter2,
			BEPowerVector           *powerVector,
		 	BCVMaximizationMatrix   matrix,				/* moments 				*/
		 	UINT4 			kMin,				/* position of the integration		*/
		 	double 			psi0,				/* change to real4 or real 8		*/
		 	double 			psi3)				/* idem 				*/
{
  UINT4			i;
  UINT4  		n = Filter1->length;
  UINT4	  		nby2 = Filter1->length / 2;
	  
  REAL8 		amplitude;
  REAL8	  		cos_phase;
  REAL8	  		sin_phase;
  
  /* Create the templates */	       
  for (i = 0; i< kMin-1; i++)
    {
      Filter1->data[i] = 0;
      Filter2->data[i] = 0;
    }
  /* Should add the upper fendBCV here as well */

  for (i = kMin; i< nby2; i++)
    {
      amplitude  =   psi0 * powerVector->fm5_3.data[i]  				/* The same for the two filters 	*/
	+ psi3 * powerVector->fm2_3.data[i];
      
      cos_phase  = cos(amplitude);						/* The same for the two filters		*/
      sin_phase  = sin(amplitude);
      

      /* Fill the first filter here */
      amplitude  = matrix.a11 * powerVector->fm7_6.data[i];
      
      Filter1->data[i]   =  amplitude * cos_phase;      
      Filter1->data[n-i] = -amplitude * sin_phase;

      /* Fill the second one here */
      amplitude =  matrix.a21 * powerVector->fm7_6.data[i] + 
	matrix.a22 * powerVector->fm1_2.data[i];
      
      Filter2->data[i]   =  amplitude * cos_phase;      
      Filter2->data[n-i] = -amplitude * sin_phase;
    }  
}




void
LALWaveOverlapBCV2(	     LALStatus               *status,
			     REAL4Vector             *correlation,
			     InspiralWaveOverlapIn   *overlapin,
			     REAL4Vector             *Filter1,
			     REAL4Vector             *Filter2,
			     BCVMaximizationMatrix    matrix,
			     OtherParamIn             otherIn ,
			     OverlapOutputIn          *OverlapOutput
		  	     )
{
 
  REAL4   rhoMaxConstraint=0   ,   alphaConstraint=0,   phaseConstraint=0;
  REAL4   rhoMaxUnconstraint=0;
  INT4    rhoBinUnconstraint=0, rhoBinConstraint=0;
  REAL4   alphaUnconstraint=0,   phaseUnconstraint=0;
  REAL4   rhoConstraint=0, rhoUnconstraint=0, inter1, inter2;



  
  REAL4  thetab,thetav=0, alphaMax, a11,a22,a21, 
	  V0, V1, V2, 								/* idem 				*/
	  rho		= 0,							/* the SNR 				*/
            w, thetac;
  

  InspiralWaveCorrelateIn 
	  corrin;								/* the correlation input structure 	*/
  UINT4 
	  i, 									
	  nBegin, 								/* beginning of valid correlation 	*/
	  nEnd; 								/* end of valid correlation 		*/
	
  

	  
  INITSTATUS (status, "LALInspiralWaveOverlapBCV2", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);
  

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
  corrin.signal2        = *Filter1;						/* The first template 			*/
  LALInspiralWaveCorrelate(status->statusPtr, &X1, corrin);			/* the correlation			*/
  CHECKSTATUSPTR(status);							
  LALGetOrthogonalFilter2(Filter1);						/* get the orthonormalized template     */
  corrin.signal2        = *Filter1;						/* the second template 		*/
  LALInspiralWaveCorrelate(status->statusPtr, &X3, corrin);			/* the correlation 			*/
  CHECKSTATUSPTR(status);
  
  /* Filter h2 and its orthornormal filter */
  corrin.signal2        = *Filter2;						/* The first template			*/
  LALInspiralWaveCorrelate(status->statusPtr, &X2, corrin);			/* the correlation 			*/
  CHECKSTATUSPTR(status);		
  LALGetOrthogonalFilter2( Filter2);						/* get the orthonormailzed templates 	*/
  corrin.signal2        = *Filter2;						/* the second template			*/
  LALInspiralWaveCorrelate(status->statusPtr, &X4, corrin);			/* the correlation 			*/
  
  /* nbegin and end for the correlation processus 
   * as well as overlapout affectation */
  nBegin                = overlapin->nBegin;					/* starting bin 			*/
  nEnd              	= Filter1->length - overlapin->nEnd;			/* ending valid bin 			*/


  a11 = matrix.a11;
  a22 = matrix.a22;
  a21 = matrix.a21;
  alphaMax = pow(corrin.fCutoff, -2./3.);

  thetab = fabs(-(a11 * alphaMax)/(a22+a21*alphaMax));
  thetab = atan(thetab);
  rho    = 0.;
  
  thetac = atan(-a11/a21);
  if(otherIn.PrintBestOverlap){
  fprintf(stderr, "a11=%e, a21=%e a22=%e %e thetab=%e, thetac=%e\n", 
	  matrix.a11, matrix.a21, matrix.a22,corrin.fCutoff, thetab, thetac);  
  }

  /*fprintf(stderr,"%e %e %e %e %e %e\n", alphaMax , thetab, corrin.fCutoff, a11, a22, a21);*/

  /* Once we have the 4 correlations, we can search for the
   * maximum of correlation and compute the final SNR */
  for (i = nBegin; i < nEnd; i++) 
    {
      
      inter1 = X1.data[i] * X1.data[i] +  X3.data[i] * X3.data[i];
      inter2 = X2.data[i] * X2.data[i] +  X4.data[i] * X4.data[i];
      V0 = inter1 + inter2;
      V1 = inter1 - inter2;
      V2 =  2*(X1.data[i]*X2.data[i] + X3.data[i]*X4.data[i]);
      
      rhoUnconstraint = ((V0 + sqrt(V1*V1+V2*V2))/2.);
      

      /* get thetav first in order to constraint alpha_F*/
      thetav = atan2(V2,V1);
      thetav*=.5;
	
      
      if (otherIn.alphaFConstraint == ALPHAFConstraint){
	
	if (-thetab < thetav && thetav < 0){	  
	  rhoConstraint = rhoUnconstraint;
	  w = thetav;
	}
	else if ((0 <= thetav && thetav < thetac) ){
	  rhoConstraint = ((V0 + V1)/2.);	  
	  w = 0.;
	}
	else if( (thetac <= thetav && thetav<LAL_PI_2+1e-6)||( -LAL_PI_2-1e-6 < thetav && thetav<= -thetab)){
	  rhoConstraint =((V0+V1*cos(2*thetab)+V2*sin(2*thetab))/2.);;
	  w = -thetab;
	}
	else 
	  {
	    fprintf(stderr,"must not enter here thetac = %e thetav = %e thetab=%e\n ",thetac, thetav , thetab);
	    exit(0);
	    w = -1;
	  }
      	  
	 if(otherIn.PrintBestOverlap || otherIn.PrintSNRHisto || otherIn.snrAtCoaTime) {
	   correlation->data[i] = rhoConstraint;
	 }



	if ( rhoConstraint > rhoMaxConstraint)			
	  /* Keep the position of the max only	*/
	  {
	    rhoMaxConstraint 	  = rhoConstraint;
	    rhoBinConstraint 	  = i;
	    phaseConstraint       = w;
	    
	  }          
	/*we save unconstraint results as well*/
	if ( rhoUnconstraint > rhoMaxUnconstraint)							/* Keep the position of the max only	*/
	  {
	    rhoMaxUnconstraint 	= rhoUnconstraint;
	    rhoBinUnconstraint 	= i;
	    phaseUnconstraint   = .5*atan2(V2, V1);
	  }          


      } 
      else{ /*no constraint on alphaF here, we keep only results with respect to unconstraint formulaes*/
	if(otherIn.PrintBestOverlap || otherIn.PrintSNRHisto || otherIn.snrAtCoaTime) {
	  correlation->data[i] = rhoUnconstraint;
	}
	
	if ( rhoUnconstraint > rhoMaxUnconstraint)							/* Keep the position of the max only	*/
	  {
	    rhoMaxUnconstraint 	= rhoUnconstraint;
	    rhoBinUnconstraint 	= i;
	    phaseUnconstraint   = .5*atan2(V2, V1);

	  }          
      }
    }   

  /* Finally get the alpha value corresponding to the best rho */ 


  alphaConstraint = -(matrix.a22 * tan(phaseConstraint)) 				
    / (matrix.a11 + matrix.a21* tan(phaseConstraint));
  
  alphaUnconstraint = -(matrix.a22 * tan(phaseUnconstraint)) 				
    / (matrix.a11 + matrix.a21* tan(phaseUnconstraint));

 
  OverlapOutput->rhoMaxConstraint     = sqrt(rhoMaxConstraint);
  OverlapOutput->rhoBinConstraint     = rhoBinConstraint;
  OverlapOutput->alphaConstraint      = alphaConstraint;
  OverlapOutput->phaseConstraint      = phaseConstraint;
  OverlapOutput->rhoMaxUnconstraint   = sqrt(rhoMaxUnconstraint);
  OverlapOutput->rhoBinUnconstraint   = rhoBinUnconstraint;
  OverlapOutput->alphaUnconstraint    = alphaUnconstraint;
  OverlapOutput->phaseUnconstraint    = phaseUnconstraint;




 
   DETATCHSTATUSPTR(status);
  RETURN(status);
}






/*  ****************************************************************************
 *  The Overlap function for BCV templates 
 *  ***************************************************************************/
void
LALWaveOverlapBCV(	     LALStatus               *status,
			     REAL4Vector             *correlation,
			     InspiralWaveOverlapIn   *overlapin,
			     REAL4Vector             *Filter1,
			     REAL4Vector             *Filter2,
			     BCVMaximizationMatrix    matrix,
			     OtherParamIn             otherIn ,
			     OverlapOutputIn          *OverlapOutput
		  	     )
     /*  </lalVerbatim>  */
{
 
  REAL4   rhoMaxConstraint=0   ,   alphaConstraint=0,   phaseConstraint=0;
  REAL4   rhoMaxUnconstraint=0;
  INT4    rhoBinUnconstraint=0, rhoBinConstraint=0;
  REAL4   alphaUnconstraint=0,   phaseUnconstraint=0;
  REAL4   rhoConstraint=0, rhoUnconstraint=0;



  
  REAL4   BestPhase=0,thetab,thetav=0, alphaMax, a11,a22,a21, 
	  x1_2, x2_2, x3_2, x4_2, 						/* some temporary data 	      		*/
	  V0, V1, V2, 								/* idem 				*/
	  rho		= 0,							/* the SNR 				*/
           alpha, w, thetac;
  
  REAL4Vector 
    template, 
    x1, x2, x3, x4,
    phaseV,
    rho1, rho2, rho3, 
    v0, v1, v2;					/* output of the correlations 		*/

  InspiralWaveCorrelateIn 
	  corrin;								/* the correlation input structure 	*/
  UINT4 
	  i, 									
	  nBegin, 								/* beginning of valid correlation 	*/
	  nEnd, 								/* end of valid correlation 		*/
	  n = correlation->length; 						/* length of the vectors		*/
  REAL4 phi, phase, cphase, sphase, cphi, sphi;
  FILE *Foutput;
	  
  INITSTATUS (status, "LALInspiralWaveOverlapBCV", BANKEFFICIENCYC);
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

  if (otherIn.PrintBestOverlap && otherIn.extraFinalPrinting){
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
  corrin.signal2        = *Filter1;						/* The first template 			*/
  LALInspiralWaveCorrelate(status->statusPtr, &x1, corrin);			/* the correlation			*/
  CHECKSTATUSPTR(status);							
  LALGetOrthogonalFilter2(Filter1);						/* get the orthonormalized template     */
  corrin.signal2        = *Filter1;						/* the second template 		*/
  LALInspiralWaveCorrelate(status->statusPtr, &x3, corrin);			/* the correlation 			*/
  CHECKSTATUSPTR(status);
  
  /* Filter h2 and its orthornormal filter */
  corrin.signal2        = *Filter2;						/* The first template			*/
  LALInspiralWaveCorrelate(status->statusPtr, &x2, corrin);			/* the correlation 			*/
  CHECKSTATUSPTR(status);		
  LALGetOrthogonalFilter2( Filter2);						/* get the orthonormailzed templates 	*/
  corrin.signal2        = *Filter2;						/* the second template			*/
  LALInspiralWaveCorrelate(status->statusPtr, &x4, corrin);			/* the correlation 			*/
  
  /* nbegin and end for the correlation processus 
   * as well as overlapout affectation */
  nBegin                = overlapin->nBegin;					/* starting bin 			*/
  nEnd              	= Filter1->length - overlapin->nEnd;			/* ending valid bin 			*/



  /* some output for debugging, checking ... */
  if (otherIn.PrintBestOverlap && otherIn.extraFinalPrinting)
    {
      Foutput = fopen( "BE_Filter.dat", "w");
      for (i=0; i< x1.length; i++)
	fprintf(Foutput, "%e %e\n", i/corrin.samplingRate,fabs(x1.data[i]));
      fprintf(Foutput,"&\n");	       
      
      for (i=1; i< Filter1->length; i++)
	fprintf(Foutput, "%e %e\n", i/corrin.samplingRate,fabs(x2.data[i]));
      fprintf(Foutput,"&\n");	       
      
      for (i=0; i< Filter1->length; i++)
	fprintf(Foutput, "%e %e\n", i/corrin.samplingRate,fabs(x3.data[i]));
      fprintf(Foutput,"&\n");	       
      
      for (i=0; i< Filter1->length-1; i++)
	fprintf(Foutput, "%e %e\n", i/corrin.samplingRate ,fabs(x4.data[i]));
      fprintf(Foutput,"&\n");	 
      fclose(Foutput);      
    }

  a11 = matrix.a11;
  a22 = matrix.a22;
  a21 = matrix.a21;
  alphaMax = pow(corrin.fCutoff, -2./3.);

  thetab = fabs(-(a11 * alphaMax)/(a22+a21*alphaMax));
  thetab = atan(thetab);
  rho    = 0.;
  
  thetac = atan(-a11/a21);
  if(otherIn.PrintBestOverlap){
  fprintf(stderr, "a11=%e, a21=%e a22=%e %e thetab=%e, thetac=%e\n", 
	  matrix.a11, matrix.a21, matrix.a22,corrin.fCutoff, thetab, thetac);  
  }

  /*fprintf(stderr,"%e %e %e %e %e %e\n", alphaMax , thetab, corrin.fCutoff, a11, a22, a21);*/

  /* Once we have the 4 correlations, we can search for the
   * maximum of correlation and compute the final SNR */
  for (i = nBegin; i < nEnd; i++) 
    {
      x1_2 = x1.data[i] * x1.data[i];
      x2_2 = x2.data[i] * x2.data[i];
      x3_2 = x3.data[i] * x3.data[i];
      x4_2 = x4.data[i] * x4.data[i];      
      
      /*
	inter1 = x1.data[i] * x1.data[i] +  x3.data[i] * x3.data[i];
	inter2 = x2.data[i] * x2.data[i] +  x4.data[i] * x4.data[i];
	V0 = inter1 + inter2;
	V1 = inter1 - inter2;
	V2 =  2*(x1.data[i]*x2.data[i] + x3.data[i]*x4.data[i]);
      */


      V0 = x1_2 + x2_2 + x3_2 + x4_2;
      V1 = x1_2 + x3_2 - x2_2 - x4_2;
      V2 = 2*(x1.data[i]*x2.data[i] + x3.data[i]*x4.data[i]);
      
      rhoUnconstraint = sqrt((V0 + sqrt(V1*V1+V2*V2))/2.);
      

      /* get thetav first in order to constraint alpha_F*/
      thetav = atan2(V2,V1);
      thetav*=.5;
	
      if(otherIn.PrintBestOverlap && otherIn.extraFinalPrinting){
	rho1.data[i] = rhoUnconstraint;
	rho2.data[i] = sqrt((V0 + V1)/2.);
	rho3.data[i] = sqrt((V0+V1*cos(2*thetab)+V2*sin(2*thetab))/2.);
	v0.data[i]   = V0;
	v1.data[i]   = V1;
	v2.data[i]   = V2;
	phaseV.data[i] =thetav;		             
      }
      
      
      if (otherIn.alphaFConstraint == ALPHAFConstraint){
	
	if (-thetab < thetav && thetav < 0){	  
	  rhoConstraint = rhoUnconstraint;
	  w = thetav;
	}
	else if ((0 <= thetav && thetav < thetac) ){
	  rhoConstraint = sqrt((V0 + V1)/2.);	  
	  w = 0.;
	}
	else if( (thetac <= thetav && thetav<LAL_PI_2+1e-6)||( -LAL_PI_2-1e-6 < thetav && thetav<= -thetab)){
	  rhoConstraint =sqrt((V0+V1*cos(2*thetab)+V2*sin(2*thetab))/2.);;
	  w = -thetab;
	}
	else 
	  {
	    fprintf(stderr,"must not enter here thetac = %e thetav = %e thetab=%e\n ",thetac, thetav , thetab);
	    exit(0);
	    w = -1;
	  }
      	  
	 if(otherIn.PrintBestOverlap || otherIn.PrintSNRHisto || otherIn.snrAtCoaTime) {
	   correlation->data[i] = rhoConstraint;
	 }



	if ( rhoConstraint > rhoMaxConstraint)			
	  /* Keep the position of the max only	*/
	  {
	    rhoMaxConstraint 	  = rhoConstraint;
	    rhoBinConstraint 	  = i;
	    phaseConstraint       = w;
	    
	  }          
	/*we save unconstraint results as well*/
	if ( rhoUnconstraint > rhoMaxUnconstraint)							/* Keep the position of the max only	*/
	  {
	    rhoMaxUnconstraint 	= rhoUnconstraint;
	    rhoBinUnconstraint 	= i;
	    phaseUnconstraint   = .5*atan2(V2, V1);
	  }          


      } 
      else{ /*no constraint on alphaF here, we keep only results with respect to unconstraint formulaes*/
	if(otherIn.PrintBestOverlap || otherIn.PrintSNRHisto || otherIn.snrAtCoaTime) {
	  correlation->data[i] = rhoUnconstraint;
	}
	
	if ( rhoUnconstraint > rhoMaxUnconstraint)							/* Keep the position of the max only	*/
	  {
	    rhoMaxUnconstraint 	= rhoUnconstraint;
	    rhoBinUnconstraint 	= i;
	    phaseUnconstraint   = .5*atan2(V2, V1);

	  }          
      }
    }   

  /* Finally get the alpha value corresponding to the best rho */ 


  alphaConstraint = -(matrix.a22 * tan(phaseConstraint)) 				
    / (matrix.a11 + matrix.a21* tan(phaseConstraint));
  
  alphaUnconstraint = -(matrix.a22 * tan(phaseUnconstraint)) 				
    / (matrix.a11 + matrix.a21* tan(phaseUnconstraint));

 
  OverlapOutput->rhoMaxConstraint     = rhoMaxConstraint;
  OverlapOutput->rhoBinConstraint     = rhoBinConstraint;
  OverlapOutput->alphaConstraint      = alphaConstraint;
  OverlapOutput->phaseConstraint      = phaseConstraint;
  OverlapOutput->rhoMaxUnconstraint   = rhoMaxUnconstraint;
  OverlapOutput->rhoBinUnconstraint   = rhoBinUnconstraint;
  OverlapOutput->alphaUnconstraint    = alphaUnconstraint;
  OverlapOutput->phaseUnconstraint    = phaseUnconstraint;



  /* The final template to be print ? */
  if (otherIn.PrintBestOverlap  && otherIn.extraFinalPrinting){

    /*print overlap, combinaison of the x_i filters*/
    Foutput=fopen("BE_Phase.dat","w");
    for (i=0; i< phaseV.length; i++){
      fprintf(Foutput, "%e %e\n", i/corrin.samplingRate,phaseV.data[i]);
    }
    fclose(Foutput);
    

  Foutput=fopen("BE_rho1.dat","w");
    for (i=0; i< phaseV.length; i++){
      fprintf(Foutput, "%e %e\n", i/corrin.samplingRate,rho1.data[i]);
    }
    fclose(Foutput);
    

  Foutput=fopen("BE_rho2.dat","w");
    for (i=0; i< phaseV.length; i++){
      fprintf(Foutput, "%e %e\n", i/corrin.samplingRate,rho2.data[i]);
    }
    fclose(Foutput);
    

    Foutput=fopen("BE_rho3.dat","w");
    for (i=0; i< phaseV.length; i++){
      fprintf(Foutput, "%e %e\n", i/corrin.samplingRate,rho3.data[i]);
    }
    fclose(Foutput);
     
    Foutput=fopen("BE_v0.dat","w");
    for (i=0; i< phaseV.length; i++){
      fprintf(Foutput, "%e %e\n", i/corrin.samplingRate,v0.data[i]);
    }
    fclose(Foutput);

    Foutput=fopen("BE_v1.dat","w");
    for (i=0; i< phaseV.length; i++){
      fprintf(Foutput, "%e %e\n", i/corrin.samplingRate,v1.data[i]);
    }
    fclose(Foutput); 

    Foutput=fopen("BE_v2.dat","w");
    for (i=0; i< phaseV.length; i++){
      fprintf(Foutput, "%e %e\n", i/corrin.samplingRate,v2.data[i]);
    }
    fclose(Foutput);
    
    Foutput=fopen("BE_alpha.dat","w");
    for (i=0; i< correlation->length; i++){
      alpha = -(matrix.a22 * tan(phaseV.data[i])) 				/* Compute the final alpha parameter 	*/
	/ (matrix.a11 + matrix.a21* tan(phaseV.data[i]));
      fprintf(Foutput, "%e %e\n", i/corrin.samplingRate,alpha*pow(corrin.fCutoff,2./3.));
    }
    fclose(Foutput);

    /*print overlap, combinaison of the x_i filters*/
    Foutput=fopen("BE_Overlap.dat","w");
    for (i=0; i< correlation->length; i++){
      fprintf(Foutput, "%e %e\n", i/corrin.samplingRate,fabs(correlation->data[i]));
    }
    fclose(Foutput);
  }


    if (otherIn.extraFinalPrinting && otherIn.PrintBestOverlap){
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


    LALGetOrthogonalFilter2(Filter1);						/* get the orthonormalized template     */
    corrin.signal2        = *Filter1;

    for (i=0; i<(UINT4)correlation->length; i++) 
    { 
      cphase = cos(phaseV.data[i]);
      cphi   = cos(phi);
      sphase = sqrt(1-cphase*cphase);
      sphi   = sqrt(1-cphi*cphi);
      
      template.data[i]=Filter1->data[i]*cphase*cphi; 
      template.data[i]+=corrin.signal2.data[i] * sphase*cphi; 

	
    }

    LALGetOrthogonalFilter2(Filter2);						/* get the orthonormalized template     */
    corrin.signal2        = *Filter2;

    for (i=0; i<(UINT4)correlation->length; i++) 
    { 
      cphase = cos(phaseV.data[i]);
      cphi   = cos(phi);
      sphase = sqrt(1-cphase*cphase);
      sphi   = sqrt(1-cphi*cphi);
      
      template.data[i]=Filter2->data[i]*sphase*cphi; 
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

  if (otherIn.extraFinalPrinting && otherIn.PrintBestOverlap){
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
      fprintf(output, "#psi0 psi3 nLayer totalMass fFinal\n");
      
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
      		fprintf(output, "%e %e %d %e %e\n", 	(*list)[i].params.psi0, 
				  			(*list)[i].params.psi3,
				  			(*list)[i].nLayer, 
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

/* personal function. Should not be used. output hte bank  overlaps in a binary
 * file to be used with DES3.
 * */

void
PrintBankOverlap(
		 InspiralTemplateList 	**list,
		 int 			sizeBank,
		 float 			*overlap,
		 InspiralCoarseBankIn 	coarseBankIn
		)
{


/*to avoid warnings */
int temp, temp3;
float temp2,temp1;

temp = sizeBank;
temp1=*overlap;
temp2 = (*list)[0].metric.g00;
temp3 = coarseBankIn.numFcutTemplates;
}

#if 0
FILE 		*output1;
  FILE 		*output2;
  
  long Npsi0, Npsi3;

  double dx0, dx3;
  double psi0Min = coarseBankIn.psi0Min;
  double psi0Max = coarseBankIn.psi0Max;
  double psi3Min = coarseBankIn.psi3Min;
  double psi3Max = coarseBankIn.psi3Max;
  double psi0, psi3;

  int    numfcut = coarseBankIn.numFcutTemplates;
  int    i,j,l,n;

  float  *a;

  double minimalMatch = coarseBankIn.mmCoarse;
  double theta, myphi, fac;


  output1 = fopen("FF.sr4", "a+");
  output2 = fopen("FF.dim", "w");
   dx0 = sqrt(2.L * (1.L - minimalMatch)/(*list)[0].metric.g00 );
   dx3 = sqrt(2.L * (1.L - minimalMatch)/(*list)[0].metric.g11 );


   fprintf(stderr, "%f %f\n", dx0, dx3);
   
     if ((*list)[0].metric.theta!=0.L)
   {
     myphi = atan2(dx3, dx0);
     theta = fabs((*list)[0].metric.theta);
     if (theta <= myphi) {
       fac = cos(theta);
       dx0  /= fac;
       dx3 *=  fac;
     } 
     else {
       fac = sin(theta);
       dx0 = dx3 / fac;
       dx3 = dx0 * fac;
     }
   }
   
     fprintf(stderr, "%f %f\n", dx0, dx3);
     fprintf(stderr,"Approxi= %d\n", coarseBankIn.approximant); 
     a =( float*)malloc(sizeof(float) ); 

     switch( coarseBankIn.approximant ){
	case BCV: 
       /*some data*/
       
       
       Npsi0 = floor(( psi0Max - psi0Min ) / dx0) + 1 ;
       Npsi3 = floor(-( psi3Min - psi3Max ) / dx3) + 1;

       printf("%f %f %ld %ld\n", dx0, dx3, Npsi0, Npsi3);     
       
/*	Npsi0 = 56; */
	Npsi3 = sizeBank/Npsi0/5;  
	dx0  =  ( psi0Max - psi0Min)/(Npsi0);
	dx3  =  ( psi3Max - psi3Min)/(Npsi3);
	
       printf("%f %f %ld %ld\n", dx0, dx3, Npsi0, Npsi3);     
       fflush(stdout); 
       /* The dimension file */
       fprintf(output2,"%7ld %14.8f %14.8f  %s\n",Npsi3, -psi3Max+dx3, dx3,"psi3");
       fprintf(output2,"%7ld %14.8f %14.8f  %s\n",Npsi0, psi0Min+dx0, dx0,"psi0");
       fprintf(output2,"%7d %14.8f %14.8f  %s\n",numfcut,1.,1.,"layers");
       
       fclose(output2);
       float tab[Npsi0][Npsi3][numfcut];
       for (l = 0; l < numfcut; l++){
	 for (i = 0; i < Npsi0; i++){
	   for (j = 0; j < Npsi3; j++){
	     tab[i][j][l]=0;
	   }
	 }
       }
       n = 0;
       /*fill to zero */
       while (n< sizeBank)
	 {
	   psi0 = (*list)[n].params.psi0;
	   psi3 = (*list)[n].params.psi3;
	   
	   i = (int)(( psi0 - psi0Min ) / dx0);
           j = (int)(-( psi3 - psi3Max) / dx3);;
	 if (i > Npsi0 || (long)j > Npsi3 || (*list)[n].nLayer > (unsigned int)numfcut)
		fprintf(stderr," %d %d %d %f %f \n", n, i, j, psi0, psi3); 
	  
  	  tab[i][j][(*list)[n].nLayer-1] = overlap[n];
	  n++;
	}
/*fprintf(stderr, "ok %d %d %d\n",numfcut, Npsi0, Npsi3);*/
      for (l = 1; l <= numfcut; l++){
       for (i = 0; i < Npsi0; i++){
	for (j = 0; j < Npsi3; j++){
/*        fprintf(stderr,"%ld %ld %ld %lf\n", i , j , l , tab[i][j][l]); */
       if (tab[i][j][l-1]!=0)
	a[0] = tab[i][j][l-1];	
	else a[0] = 0;
/*fprintf(stderr,"#%d %d %d %lf\n",i, j, l-1, a[0]);	*/
	fwrite(a,4,1,output1);
	 }
	}	 	
      }
      
    fclose(output1);
      break;
    case BCVSpin: 
    case SpinTaylorT3: 
    case TaylorT1: 
    case TaylorT2: 
    case TaylorT3: 
    case TaylorF1: 
    case TaylorF2: 
    case PadeT1: 
    case PadeF1: 
    case EOB:
    case SpinTaylor: 
      break;      
  }

}
#endif





/* routine to print the bank in a xml file*/
/* should add the status variable since we call sone lal functions here 
 * this routine is not well commented and won't be in close future.
 * */
/* !! tau2 is used to store the Layer number */
void 
BEPrintBankXml(
	       InspiralTemplateList *coarseList, 
	       UINT4 numCoarse,
	       InspiralCoarseBankIn   coarseBankIn,
	       RandomInspiralSignalIn randIn,
	       OtherParamIn           otherIn)
{
#define MAXIFO 2
  LALStatus             status = blank_status;
  MetadataTable         proctable;

  MetadataTable         templateBank;
  SnglInspiralTable     *tmplt  = NULL;
  CHAR  ifo[3];                           /* two character ifo code       */
  CHAR *channelName = NULL;               /* channel string               */
  UINT4 	i;  
  LIGOLwXMLStream       xmlStream;
  CHAR  fname[256];
  LIGOTimeGPS gpsStartTime 	= { 0, 0 };    /* input data GPS start time    */
  LIGOTimeGPS gpsEndTime 	= { 0, 0 };      /* input data GPS end time      */
  LALLeapSecAccuracy    accuracy = 1;
  CHAR  comment[LIGOMETA_COMMENT_MAX];
  CHAR  ifoName[MAXIFO][LIGOMETA_IFO_MAX];

  MetadataTable         processParamsTable;
  ProcessParamsTable   *this_proc_param = NULL;

  
  strncpy( ifoName[0], "no", LIGOMETA_IFO_MAX * sizeof(CHAR) );
  strncpy( ifoName[1], "ne", LIGOMETA_IFO_MAX * sizeof(CHAR) );
  memset( ifo, 0, sizeof(ifo) );
  memcpy( ifo, "MC", sizeof(ifo) - 1 );


  /* --- first we create the filename --- */
  LALSnprintf( fname, sizeof(fname), BANKEFFICIENCY_PRINTBANK_FILEXML ,
	       ifo, gpsStartTime.gpsSeconds,
	       gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  
  /* --- let's save the bank itself --- */
  templateBank.snglInspiralTable = NULL;

  

  /* --- convert the templates to sngl_inspiral structures 
     in order to  write results in the XML file --- */
  if ( numCoarse )
    {
      /* --- we sae one line --- */
      tmplt = templateBank.snglInspiralTable = (SnglInspiralTable *)
	LALCalloc( 1, sizeof(SnglInspiralTable) );
      LALSnprintf( tmplt->ifo, LIGOMETA_IFO_MAX * sizeof(CHAR), ifo );
      LALSnprintf( tmplt->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR), 
		   "BankEfficiency" );
      LALSnprintf( tmplt->channel, LIGOMETA_CHANNEL_MAX * sizeof(CHAR),
		   channelName );
      tmplt->mass1   = (REAL4) coarseList[0].params.mass1;
      tmplt->mass2   = (REAL4) coarseList[0].params.mass2;
      tmplt->mchirp  = (REAL4) coarseList[0].params.chirpMass;
      tmplt->eta     = (REAL4) coarseList[0].params.eta;
      tmplt->tau0    = (REAL4) coarseList[0].params.t0;
      tmplt->tau2    = (REAL4) coarseList[0].nLayer;
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
	  LALSnprintf( tmplt->ifo, LIGOMETA_IFO_MAX * sizeof(CHAR), ifo );
	  LALSnprintf( tmplt->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR), 
		       "BankEfficiency" );
	  LALSnprintf( tmplt->channel, LIGOMETA_CHANNEL_MAX * sizeof(CHAR),
		       channelName );
	  tmplt->mass1   = (REAL4) coarseList[i].params.mass1;
	  tmplt->mass2   = (REAL4) coarseList[i].params.mass2;
	  tmplt->mchirp  = (REAL4) coarseList[i].params.chirpMass;
	  tmplt->eta     = (REAL4) coarseList[i].params.eta;
	  tmplt->tau0    = (REAL4) coarseList[i].params.t0;
	  tmplt->tau2    = (REAL4) coarseList[i].nLayer;
	  tmplt->tau3    = (REAL4) coarseList[i].params.t3;
	  tmplt->tau4    = (REAL4) coarseList[i].params.t4;
	  tmplt->tau5    = (REAL4) coarseList[i].params.t5;
	  tmplt->ttotal  = (REAL4) coarseList[i].params.tC;
	  tmplt->psi0    = (REAL4) coarseList[i].params.psi0;
	  tmplt->psi3    = (REAL4) coarseList[i].params.psi3;
	  tmplt->f_final = (REAL4) coarseList[i].params.fFinal;
	}
    }
  
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

  BEFillProc(this_proc_param, coarseBankIn, randIn, otherIn);




  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );



  /* write process table */
  LALSnprintf( proctable.processTable->ifos, LIGOMETA_IFOS_MAX, "%s%s", 
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


  /* --- desallocate memory here --- */
  while ( templateBank.snglInspiralTable )
    {
      tmplt = templateBank.snglInspiralTable;
      templateBank.snglInspiralTable = templateBank.snglInspiralTable->next;
      LALFree( tmplt );
    }
  
  /* close the output xml file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlStream ), &status );
}



/* function to get matrix components given 
 * the vectors of the different moments   */
void BEGetMatrixFromVectors(BEMoments moments,
	              UINT4 k,
	              BCVMaximizationMatrix *matrix2fill
	              )
{
  matrix2fill->a11 = moments.a11.data[k];
  matrix2fill->a21 = moments.a21.data[k];
  matrix2fill->a22 = moments.a22.data[k];   
}


void 
BEFillProc(
	       ProcessParamsTable     *this_proc_param,
	       InspiralCoarseBankIn   coarseBankIn,
	       RandomInspiralSignalIn randIn,
	       OtherParamIn           otherIn)
{


#define ADD_PROCESS_PARAM( pptype, format, ppname, ppvalue ) \
this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
  calloc( 1, sizeof(ProcessParamsTable) ); \
  LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME ); \
  LALSnprintf( this_proc_param->param,   LIGOMETA_PARAM_MAX,   "%-20s", ppname); \
  LALSnprintf( this_proc_param->type,    LIGOMETA_TYPE_MAX,    "%-10s",   pptype ); \
  LALSnprintf( this_proc_param->value,   LIGOMETA_VALUE_MAX,   format, ppvalue );

  ADD_PROCESS_PARAM("float",	"%12.5f",	"--alpha-bank",	        coarseBankIn.alpha);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--alpha-signal",	randIn.param.alpha);
  ADD_PROCESS_PARAM("float",	"%d",    	"--debug",	        otherIn.lalDebug);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--fl",		        coarseBankIn.fLower);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--fendbcv",		coarseBankIn.LowGM);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--fendbcv",		coarseBankIn.HighGM);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--fFinal-signal",	randIn.param.fCutoff);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--freq-moment-bank",	coarseBankIn.fUpper);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--sampling",		coarseBankIn.tSampling);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--mass-range",		randIn.mMin);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--mass-range",		randIn.mMax);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--mass-range-bank",	coarseBankIn.mMin);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--mass-range-bank",	coarseBankIn.mMax);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--m1",			otherIn.m1);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--m2",			otherIn.m2);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--mm",			coarseBankIn.mmCoarse);
  ADD_PROCESS_PARAM("int",	"%6d",		"--n",			otherIn.ntrials);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--noise-amplitude",	randIn.NoiseAmp);
  ADD_PROCESS_PARAM("int",	"%6d",		"--number-fcut",	coarseBankIn.numFcutTemplates);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--psi0",		otherIn.psi0);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--psi3",		otherIn.psi3);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--psi0-range",		coarseBankIn.psi0Min);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--psi0-range",		coarseBankIn.psi0Max);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--psi3-range",		coarseBankIn.psi3Min);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--psi3-range",		coarseBankIn.psi3Max);
  ADD_PROCESS_PARAM("int",	"%6d",		"--simulation-type",	randIn.type);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--signal-amplitude",	randIn.SignalAmp);
  ADD_PROCESS_PARAM("int",	"%6d",		"--signal-order",	randIn.param.order);
  ADD_PROCESS_PARAM("int",	"%6d",		"--seed",		randIn.useed);
  switch (otherIn.NoiseModel){
  case LIGOI:
    ADD_PROCESS_PARAM("string",	"%6s",		"--noisemodel",		"LIGOI");
    break;
  case LIGOA:
    ADD_PROCESS_PARAM("string",	"%6s",		"--noisemodel",		"LIGOA");
    break;
  case VIRGO:
    ADD_PROCESS_PARAM("string",	"%6s",		"--noisemodel",		"VIRGO");
    break;
  case GEO:
    ADD_PROCESS_PARAM("string",	"%6s",		"--noisemodel",		"GEO");
    break;
  case TAMA:    
    ADD_PROCESS_PARAM("string",	"%6s",		"--noisemodel",		"TAMA");    
    break; 
  case UNITY:    
    ADD_PROCESS_PARAM("string",	"%6s",		"--noisemodel",		"UNITY");    
    break;
  case REALPSD:
     ADD_PROCESS_PARAM("string", "%6s",		"--noisemodel",		"REALPSD");    
    break;
  }
  ADD_PROCESS_PARAM("int",	"%6d",		"--template",		otherIn.template);
  ADD_PROCESS_PARAM("int",	"%6d",		"--signal",		otherIn.signal);
  ADD_PROCESS_PARAM("int",	"%6d",		"--template-order",	coarseBankIn.order);

  switch (otherIn.overlapMethod){
  case InQuadrature:
    ADD_PROCESS_PARAM("string",	"%s", 	"OverlapMethod", 	"InQuadrature");
    break;
  case AlphaMaximization:
    ADD_PROCESS_PARAM("string",	"%s", 	"OverlapMethod", 	"AlphaMaximization");
    break;
  }
  
	  
#undef ADD_PROCESS_PARAM
}


/* xml file for the standalone code */
void 
BEPrintResultsXml( InspiralCoarseBankIn   coarseBankIn,
		   RandomInspiralSignalIn randIn,
		   OtherParamIn           otherIn,
		   ResultIn trigger
		  )
{
  LALStatus             status = blank_status;


  MetadataTable         templateBank;
  CHAR                  ifo[3];                           /* two character ifo code       */
  LIGOLwXMLStream       xmlStream;
  CHAR                  fname[256];
  LIGOTimeGPS           gpsStartTime 	= { 0, 0 };    /* input data GPS start time    */
  LIGOTimeGPS           gpsEndTime 	= { 0, 0 };      /* input data GPS end time      */
  LALLeapSecAccuracy    accuracy = 1;
  CHAR                  comment[LIGOMETA_COMMENT_MAX];
  CHAR                  ifoName[MAXIFO][LIGOMETA_IFO_MAX];

  MetadataTable         processParamsTable;
  ProcessParamsTable   *this_proc_param = NULL;

  LALSnprintf( fname, sizeof(fname), 
	       BANKEFFICIENCY_PRINTRESULT_FILEXML ,
	       ifo,
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
    
    BEFillProc(this_proc_param, coarseBankIn, randIn, otherIn);
    
    memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );
    
    /* write process table */
    LALSnprintf( templateBank.processTable->ifos, LIGOMETA_IFOS_MAX, "%s%s", 
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
    
    fprintf(xmlStream.fp,BANKEFFICIENCY_PARAMS_ROW"\n",
	    trigger.psi0_trigger,
	    trigger.psi3_trigger,
	    trigger.psi0_triggerC,
	    trigger.psi3_triggerC,
	    trigger.psi0_inject, 
	    trigger.psi3_inject,
	    trigger.fend_trigger, 
	    trigger.fend_triggerC, 
	    trigger.fend_inject,
	    trigger.totalMass_trigger,
	    trigger.eta_trigger,
	    trigger.totalMass_triggerC,
	    trigger.eta_triggerC,
	    trigger.mass1_inject,
	    trigger.mass2_inject,
	    trigger.rho_final,
	    trigger.phase,
	    trigger.alpha,
	    trigger.alpha_f, 
	    trigger.layer,
	    trigger.bin,
	    trigger.rho_finalC,
	    trigger.phaseC,
	    trigger.alphaC,
	    trigger.alpha_fC, 
	    trigger.layerC,
	    trigger.binC,
	    trigger.coaTime,
	    trigger.snrAtCoaTime, 
	    trigger.snrCAtCoaTime
	    );
	    
     if (trigger.ntrial == otherIn.ntrials){
       PRINT_LIGOLW_XML_TABLE_FOOTER(xmlStream.fp);
       PRINT_LIGOLW_XML_FOOTER(xmlStream.fp);
     }
     fclose( xmlStream.fp );
     xmlStream.fp = NULL;
  }
  else 
    {
      xmlStream.fp = fopen(fname,"a+");
      
      if (trigger.ntrial == otherIn.ntrials){
	
	fprintf(xmlStream.fp,BANKEFFICIENCY_PARAMS_ROW,
		trigger.psi0_trigger,
		trigger.psi3_trigger,
		trigger.psi0_triggerC,
		trigger.psi3_triggerC,
		trigger.psi0_inject, 
		trigger.psi3_inject,
		trigger.fend_trigger, 
		trigger.fend_triggerC, 
		trigger.fend_inject,
		trigger.totalMass_trigger,
		trigger.eta_trigger,
		trigger.totalMass_triggerC,
		trigger.eta_triggerC,
		trigger.mass1_inject,
		trigger.mass2_inject,
		trigger.rho_final,
		trigger.phase,
		trigger.alpha,
		trigger.alpha_f, 
		trigger.layer,
		trigger.bin,
		trigger.rho_finalC,
		trigger.phaseC,
		trigger.alphaC,
		trigger.alpha_fC, 
		trigger.layerC,
		trigger.binC, 
		trigger.coaTime,
		trigger.snrAtCoaTime, 
		trigger.snrCAtCoaTime);
	PRINT_LIGOLW_XML_TABLE_FOOTER(xmlStream.fp );
	PRINT_LIGOLW_XML_FOOTER(xmlStream.fp );
      }
      else
	{
	  fprintf(xmlStream.fp,BANKEFFICIENCY_PARAMS_ROW"\n",
		  trigger.psi0_trigger,
		  trigger.psi3_trigger,
		  trigger.psi0_triggerC,
		  trigger.psi3_triggerC,
		  trigger.psi0_inject, 
		  trigger.psi3_inject,
		  trigger.fend_trigger, 
		  trigger.fend_triggerC, 
		  trigger.fend_inject,
		  trigger.totalMass_trigger,
		  trigger.eta_trigger,
		  trigger.totalMass_triggerC,
		  trigger.eta_triggerC,
		  trigger.mass1_inject,
		  trigger.mass2_inject,
		  trigger.rho_final,
		  trigger.phase,
		  trigger.alpha,
		  trigger.alpha_f, 
		  trigger.layer,
		  trigger.bin,
		  trigger.rho_finalC,
		  trigger.phaseC,
		  trigger.alphaC,
		  trigger.alpha_fC, 
		  trigger.layerC,
		  trigger.binC, 
		  trigger.coaTime,
		  trigger.snrAtCoaTime, 
		  trigger.snrCAtCoaTime);
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
		  OtherParamIn           otherIn
		  )
{
  LALStatus             status = blank_status;


  MetadataTable         templateBank;
  CHAR                  ifo[3];                           /* two character ifo code       */
  LIGOLwXMLStream       xmlStream;
  CHAR                  fname[256];
  LIGOTimeGPS           gpsStartTime 	= { 0, 0 };    /* input data GPS start time    */
  LIGOTimeGPS           gpsEndTime 	= { 0, 0 };      /* input data GPS end time      */
  LALLeapSecAccuracy    accuracy = 1;
  CHAR                  comment[LIGOMETA_COMMENT_MAX];
  CHAR                  ifoName[MAXIFO][LIGOMETA_IFO_MAX];

  MetadataTable         processParamsTable;
  ProcessParamsTable   *this_proc_param = NULL;

  LALSnprintf( fname, sizeof(fname), BANKEFFICIENCY_PRINTPROTO_FILEXML ,
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
    
    BEFillProc(this_proc_param, coarseBankIn, randIn, otherIn);
    
    memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );
    
    
    
    /* write process table */
    LALSnprintf( templateBank.processTable->ifos, LIGOMETA_IFOS_MAX, "%s%s", 
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
BEReadXmlBank(  LALStatus  *status, 
		 CHAR *bankFileName, 
		 InspiralTemplateList **list,
		 INT4 *sizeBank,
		InspiralCoarseBankIn coarseIn)
{
  INT4 i;
  InspiralTemplate *tempPars;
  INT4 startTemplate = -1;
  INT4 stopTemplate  = -1;
  INT4                          numTmplts    = 0;
  InspiralTemplate             *bankHead     = NULL;
  InspiralTemplate             *bankCurrent  = NULL;
  InspiralTemplateNode         *tmpltHead    = NULL;
  InspiralTemplateNode         *tmpltCurrent = NULL;
  InspiralMetric metric;


  INITSTATUS( status, 
	      "BEReadXmlBank", BANKEFFICIENCYC );
  ATTATCHSTATUSPTR( status );
    
  /* how many template in the xml file ? */
  numTmplts = InspiralTmpltBankFromLIGOLw( &bankHead, 
					   bankFileName,
					   startTemplate, 
					   stopTemplate );

  if ( numTmplts < 0 )
    {
      fprintf( stderr, "error: unable to read templates from %s\n", 
	       bankFileName );
      exit( 1 );
    }
  else if ( numTmplts == 0 )
    {
      
      fprintf( stderr, "no templates found in template bank file: %s\n"
	       "exiting without searching for events.\n", bankFileName );
       
    }
   else
     fprintf( stderr, "%d templates found in template bank file: %s\n", numTmplts, bankFileName);
  

  
  if ( !
       (tempPars = (InspiralTemplate *)LALCalloc( 1, sizeof(InspiralTemplate) ))
       )
    {
      ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
    }
  
 
  bankCurrent = bankHead; 

  for (i=0; i<numTmplts; i++){
    *sizeBank = i;
    *list = (InspiralTemplateList*) 
      LALRealloc( *list, sizeof(InspiralTemplateList) * (*sizeBank + 1) );
    if ( ! *list )
      {
	LALFree( tempPars );
	ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
      }
    memset( *list + *sizeBank, 0, sizeof(InspiralTemplateList) );
          

    LAL_CALL( LALFindChirpCreateTmpltNode( status->statusPtr, 
					   bankCurrent, 
					   &tmpltCurrent ), status->statusPtr );
    if ( !tmpltHead ) tmpltHead = tmpltCurrent;
  
    
    (*list)[*sizeBank].ID = *sizeBank;
    (*list)[*sizeBank].params = *bankCurrent; 
    /* Extra info which are not in the bank */
    (*list)[*sizeBank].params.fLower = coarseIn.fLower; 
    (*list)[*sizeBank].params.tSampling = coarseIn.tSampling; 
    (*list)[*sizeBank].params.fCutoff = coarseIn.fUpper; 
    (*list)[*sizeBank].params.fFinal = coarseIn.fUpper; 
    (*list)[*sizeBank].params.order = coarseIn.order; 
    (*list)[*sizeBank].params.approximant = coarseIn.approximant; 
    (*list)[*sizeBank].params.massChoice = t03;
    (*list)[*sizeBank].params.distance =  1.;
    (*list)[*sizeBank].params.signalAmplitude= 1.;
    /*  (*list)[*sizeBank].params.nStartPad = 0;
      (*list)[*sizeBank].params.nEndPad = 0;
    */
    LALInspiralParameterCalc( status->statusPtr,  &((*list)[*sizeBank].params ));
    CHECKSTATUSPTR( status );
  
    (*list)[*sizeBank].metric = metric; 

    CHECKSTATUSPTR( status );
    ++(*sizeBank); 


    bankCurrent = bankCurrent->next ;

  }

 for ( i = 0; i < *sizeBank - 1 ; ++i )
  {
    (*list)[i].params.minMatch = (REAL4) coarseIn.mmCoarse;
    (*list)[i].params.next     = &((*list)[i+1].params);
    (*list)[i].params.fine     = NULL;
  }
  (*list)[i].params.minMatch = (REAL4) coarseIn.mmCoarse;
  (*list)[i].params.next = NULL;
  (*list)[i].params.fine = NULL;


  DETATCHSTATUSPTR( status );
  RETURN( status );


}



void
BEInitOverlapOutputIn(OverlapOutputIn *this){
  this->rhoMaxConstraint   = 0.;
  this->phaseConstraint    = 0.;
  this->rhoBinConstraint   = 0;
  this->templateNumberConstraint        = 0;
  this->alphaConstraint   = -1.;
  this->layerConstraint    = -1;
  this->freqConstraint     = -1.;

  this->rhoMaxUnconstraint = 0.;  
  this->phaseUnconstraint  = 0.;  
  this->rhoBinUnconstraint = 0;  
  this->templateNumberUnconstraint      = 0;
  this->alphaUnconstraint = -1.;
  this->layerUnconstraint  = -1;
   this->freqUnconstraint  = -1.;
 /*I do not think we need any other  iniitalisation */
}

void
BEFillOverlapOutput(InspiralWaveOverlapOut overlapout, 
		    OverlapOutputIn *this)
{

  this->rhoMaxConstraint         = overlapout.max;
  this->phaseConstraint          = overlapout.phase;
  this->rhoBinConstraint         = overlapout.bin;
  this->alphaConstraint         = -1.;
  this->freqConstraint           = -1;
  this->layerConstraint          = -1;
  this->templateNumberConstraint = -1;
  
  this->rhoMaxUnconstraint       = overlapout.max;
  this->phaseUnconstraint        = overlapout.phase;
  this->rhoBinUnconstraint       = overlapout.bin;
  this->alphaUnconstraint       = -1.;
  this->freqConstraint           = -1;
  this->layerConstraint          = -1;
  this->templateNumberConstraint = -1;
      
}



/* Estimate the size of the longest template */
void BEGetMaximumSize(LALStatus  *status, 		      
		      RandomInspiralSignalIn  randIn, 
		      UINT4 *length)
{

  INITSTATUS( status, "BEGetMaximumSize", BANKEFFICIENCYC );
  ATTATCHSTATUSPTR( status );
  
  randIn.param.massChoice 	= m1Andm2;
  randIn.param.approximant 	= EOB;
  *length = 0;
 
  LAL_CALL(LALInspiralWaveLength(status->statusPtr, length, randIn.param), 
	   status->statusPtr);

  DETATCHSTATUSPTR( status );
  RETURN( status );  

}



void BECreatePsd(LALStatus *status, 
		   InspiralCoarseBankIn *coarseBankIn, 
		   RandomInspiralSignalIn  randIn,
		   OtherParamIn           otherIn)
{
  UINT4 i,k;
  REAL4 temp;
  REAL4 df; 
  FILE *Finput, *Foutput;

  INITSTATUS( status, "BECreatePsd", BANKEFFICIENCYC );
  ATTATCHSTATUSPTR( status );


  memset( &(coarseBankIn->shf), 0, sizeof(REAL8FrequencySeries) );
  coarseBankIn->shf.f0 	= 0;
  LAL_CALL(LALDCreateVector( status->statusPtr, &(coarseBankIn->shf.data), randIn.psd.length ),
	   status->statusPtr);
  
  coarseBankIn->shf.deltaF 	= randIn.param.tSampling / ((randIn.psd.length-1)*2);
  
  /* --- Compute Noise Spectral Density --- */
  df = randIn.param.tSampling/(float) ((randIn.psd.length-1)*2);	 

  switch (otherIn.NoiseModel)						
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
      Finput = fopen(otherIn.filename,"r");
      for (i = 0; i < coarseBankIn->shf.data->length; i++)
	{
	  /* factor DeltaT since the input data have a sampling frequency of 1./DeltaT */
	  fscanf(Finput, "%f %lf\n", &temp,&coarseBankIn->shf.data->data[i]);
	  for (k = 1; k < (UINT4)( DeltaT * df); k++)
	      fscanf(Finput, "%f %f\n", &temp,&temp);
	}
      fclose(Finput);
      break;
     }

  /* --- save the psd in the file psd.dat if requested --- */
  if (otherIn.PrintPsd)
   {
     Foutput= fopen(BANKEFFICIENCY_PRINTPSD_FILE,"w");
     for (i = 1; i < coarseBankIn->shf.data->length; i++)
       fprintf(Foutput, "%f %e\n",(float)i*df, coarseBankIn->shf.data->data[i]);  
     fclose(Foutput);
   }

 

  DETATCHSTATUSPTR( status );
  RETURN( status );  
}



void BEGenerateInputData(LALStatus *status, 
			 REAL4Vector * signal,
			 RandomInspiralSignalIn  *randIn,
			 OtherParamIn           otherIn)
{
  UINT4  i;
  
  INITSTATUS( status, "BEGenerateInputData", BANKEFFICIENCYC );
  ATTATCHSTATUSPTR( status );
  
  
  for (i = 0; i < signal->length; i++) 
    signal->data[i] = 0.;       
  
  /* we might force to compute a non random waveform giving the two
     input parameters m1 and m2 */      
  
  if (otherIn.m1 == -1 && otherIn.psi0 ==-1 && otherIn.tau0 == -1)
    {   
      if (otherIn.binaryInjection == BHNS){	 
	randIn->param.massChoice = bhns;
	LAL_CALL(LALRandomInspiralSignal(status->statusPtr, signal, randIn),
		 status->statusPtr);
      }
      else
	LAL_CALL(LALRandomInspiralSignal(status->statusPtr, signal, randIn), 
		 status->statusPtr);	  
    }
  else 
    {
      if (otherIn.m1!=-1) {
	randIn->param.massChoice = fixedMasses;
	randIn->param.mass1 = otherIn.m1;
	randIn->param.mass2 = otherIn.m2;
	LAL_CALL(LALRandomInspiralSignal(status->statusPtr, signal, randIn), 
		 status->statusPtr);	  
      }
      if (otherIn.psi0!=-1) {
	randIn->param.massChoice = fixedPsi;
	randIn->param.psi0 = otherIn.psi0;
	randIn->param.psi3 = otherIn.psi3;
	LAL_CALL(LALRandomInspiralSignal(status->statusPtr, signal, randIn), 
		 status->statusPtr);	  
      }
      if (otherIn.tau0!=-1) {
	randIn->param.massChoice = fixedTau;
	randIn->param.t0 = otherIn.tau0;
	randIn->param.t3 = otherIn.tau3;
	LAL_CALL(LALRandomInspiralSignal(status->statusPtr, signal, randIn), 
		 status->statusPtr);	  
      }
    }

  DETATCHSTATUSPTR( status );
  RETURN( status );  

}
