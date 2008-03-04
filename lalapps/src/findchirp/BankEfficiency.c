/*
*  Copyright (C) 2007 Thomas Cokelaer
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

#include "BankEfficiency.h"
/* --- version information --- */
NRCSID( BANKEFFICIENCYC, "$Id$");
RCSID(  "$Id$");

#define CVS_ID_STRING_C      "$Id$"
#define CVS_REVISION_C       "$Revision$"

#define CVS_NAME_STRING_C    "$Name$"
#define CVS_SOURCE_C         "$Source$"
#define CVS_DATE_C           "$Date$"
#define PROGRAM_NAME         "BankEfficiency"



RandomParams *randParams = NULL;
INT4 randnStartPad = 0;
INT4 ascii2xml = 0;



int
main (INT4 argc, CHAR **argv ) 
{
  
  INT4          ntrials = 0; /*number of simulations*/
  INT4 	                 i;
  INT4 	                 j;
  INT4                   thisTemplateIndex;
  Order            tempOrder;  /* temporary phase order */

  /* --- input --- */
  UserParametersIn       userParam;

  /* --- signal related --- */
  REAL4Vector            signal;
  RandomInspiralSignalIn randIn;		

  /* --- template bank related --- */  
  InspiralTemplate       insptmplt;
  InspiralCoarseBankIn   coarseBankIn;
  INT4  	         sizeBank = 0;
  INT4                   filter_processed = 0;
  MetadataTable          templateBank;
  SnglInspiralTable      *tmpltHead = NULL;
  SnglInspiralTable      *tmpltCurrent = NULL;

  /* --- filtering related --- */
  REAL4Vector            correlation;
  REAL4Vector            FilterBCV1;
  REAL4Vector            FilterBCV2;
  BEPowerVector          powerVector;
  BEMoments              moments;
  InspiralWaveOverlapIn  overlapin;		
  InspiralWaveOverlapOut overlapout; 

  /* --- results and data mining --- */
  ResultIn              result;
  OverlapOutputIn       OverlapOutputThisTemplate;
  OverlapOutputIn       OverlapOutputBestTemplate;

  /* --- fft related --- */
  RealFFTPlan 		*fwdp = NULL;
  RealFFTPlan 		*revp = NULL;

  /* --- others ---*/
  LALStatus             status = blank_status;
  FILE                  *Foutput;

  /* --- fast option related --- */
  REAL4                  dt0;
  REAL4                  dt3;
  REAL4                  g00;
  REAL4                  g11;
  REAL4                  g01;
  REAL4                  match;
      
  /* --- eccentricity related to the bank */
  REAL4                  eccentricityTemplate = 0;
  REAL4                  EccentricityMin = 0;
  REAL4                  EccentricityMax = 0.8;
  REAL4                  eccentricityStep;
  REAL4                  eccentricityBins = 8;

  /* --- ambiguity function and statistics --- */
  gsl_histogram         * histogramNoise = gsl_histogram_alloc (200);  
  gsl_matrix            *amb1;
  /* --- START main code --- */
  
    
  /* --- Some initialization --- */ 
  gsl_histogram_set_ranges_uniform (histogramNoise, 0.,20.);
  lal_errhandler = LAL_ERR_EXIT;
  lalDebugLevel = 0;
  templateBank.snglInspiralTable = NULL;

  /* --- Initialization of structure --- */
  ParametersInitialization(&coarseBankIn, &randIn, &userParam);
  
  /* --- Read user parameters --- */
  ParseParameters(&argc, argv, &coarseBankIn, &randIn, &userParam);

  /* --- Check input parameters --- */
  UpdateParams(&coarseBankIn, &randIn, &userParam);
  
  /* init a random structure using possibly the seed from user input*/
  /* this call must be after the userParam have been read because it uses
   * randIn.useed as an input !*/
  LAL_CALL(LALCreateRandomParams(&status, &randParams, randIn.useed ), 
	   &status);

  /* this is a call to the function that converts the ASCII output of
   * BankEfficiency into a standard XML output. Useful when several ASCII
   * output are avalaible and needs to be converted all together into an XML
   * file.
   * */
  if (ascii2xml == 1)
  {
    BEAscii2Xml();
    exit(1);
  }

  /* 
     This is used if one want to use condor script to keep track of the
     input parameters. It creates a prototype to be read by BEAscii2xml. 
     Useful with condor. 
  */
  if (userParam.printPrototype)
  {
    BEPrintProtoXml(coarseBankIn, randIn, userParam);   
    exit(0);
  }
  
  /* --- Estimate the size of the signal --- */
  LAL_CALL(BEGetMaximumSize(&status, randIn, coarseBankIn, userParam, &(signal.length)), 
	   &status);
  
  /* --- Allocate memory for some vectors --- */
  randIn.psd.length 	= signal.length/2 + 1; 
  correlation.length 	= signal.length;
  signal.data 		= (REAL4*) LALCalloc(1, sizeof(REAL4) * signal.length);
  correlation.data 	= (REAL4*) LALCalloc(1, sizeof(REAL4) * correlation.length);
  randIn.psd.data 	= (REAL8*) LALCalloc(1, sizeof(REAL8) * randIn.psd.length); 
  if (userParam.template == BCV)
  {
    FilterBCV1.length     = signal.length;
    FilterBCV2.length     = signal.length;
    FilterBCV1.data       = (REAL4*) LALCalloc(1, sizeof(REAL4) * FilterBCV1.length);
    FilterBCV2.data       = (REAL4*) LALCalloc(1, sizeof(REAL4) * FilterBCV2.length);
  }
  
  /* --- Create the PSD noise --- */
  LAL_CALL(BECreatePsd(&status, &coarseBankIn, &randIn, userParam), 
	   &status);
 
  /* --- Create the template bank --- */
  if ( vrbflg )
  {
    fprintf( stdout, "generating template bank parameters... " ); fflush( stdout );
  }
  /* make sure the pointer to the first template is null */
  templateBank.snglInspiralTable = NULL;
  
  /*Let us replace the template bank with a fine template bank if requested */
  if (userParam.t0FineBin > 0)
  {
    LAL_CALL( LALInspiralBankGeneration2( &status, &coarseBankIn, &tmpltHead, &sizeBank, &userParam),
 	   &status );  
  }
  else
  {
    INT4 temp_order = coarseBankIn.order;
    if (coarseBankIn.order < 4){
     coarseBankIn.order = 4;  
     LAL_CALL( LALInspiralBankGeneration( &status, &coarseBankIn, &tmpltHead, &sizeBank),
 	   &status );
     coarseBankIn.order = temp_order;  
    }
    else
    {
      LAL_CALL( LALInspiralBankGeneration( &status, &coarseBankIn, &tmpltHead, &sizeBank),
 	   &status );
    }

  }
      
  
  
  amb1 = gsl_matrix_alloc(4,sizeBank); /*allocated t0,t3,max*/
  for (i=0;i<4; i++)
    for (j=0;j<sizeBank; j++)
      gsl_matrix_set(amb1, i, j, 0.);      
  
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
  if ((userParam.template == BCV) )
  {
    LAL_CALL( BECreatePowerVector(&status, &powerVector, randIn, signal.length), 
	    &status);
    LALCreateBCVMomentVector(&moments, &coarseBankIn.shf,  randIn.param.tSampling, 
			randIn.param.fLower, signal.length);
  }
  
  /*
   * Main loop is here. We create a random signal and filter 
   * it with the template bank while we want to perform one
    simulation. */

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

      if (randnStartPad==1)
      {
	randIn.param.nStartPad = (int)((float)rand()/(float)RAND_MAX*signal.length/2);
      }
      LAL_CALL(BEGenerateInputData(&status,  &signal, &randIn, userParam),
      	&status);
	

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
      /* -- !! here the approximant is overwritten. Later, we will populate
      the approximant again with userParam.template into overlapin.appoximant */     
      insptmplt = randIn.param;  
      insptmplt.order = coarseBankIn.order;
      filter_processed = 0;
      insptmplt.eccentricity = 0;

      if (userParam.eccentricBank)
        eccentricityStep = (EccentricityMax - EccentricityMin)/eccentricityBins; 
      else
      {
        eccentricityStep = 1;/*larger than max eccentricity so that only EccentricityMin is used, which we set to zero*/
        EccentricityMin = 0;
        EccentricityMax = 1;
      }
      
      for (eccentricityTemplate = EccentricityMin; eccentricityTemplate < EccentricityMax; eccentricityTemplate+=eccentricityStep)
      {
        insptmplt.eccentricity = eccentricityTemplate;
        /* -- and enter in the bank */
        for (tmpltCurrent  = tmpltHead, thisTemplateIndex=0;
	     tmpltCurrent && thisTemplateIndex < sizeBank;
	     tmpltCurrent = tmpltCurrent->next, thisTemplateIndex++)
  	{ 
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

	    case BCV:	   
	      CreateListfromTmplt(&insptmplt, tmpltCurrent);
	      insptmplt.massChoice = psi0Andpsi3;
	      LAL_CALL(LALInspiralParameterCalc( &status,  &(insptmplt) ), &status);
	      
	      if (userParam.faithfulness)
	      {
		insptmplt = randIn.param;
		overlapin.param              = randIn.param;
		overlapin.param.approximant  = userParam.template;
		sizeBank                     = 1;  /* we dont need the whole bank for checking */
                insptmplt.fFinal = randIn.param.fFinal;
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
	    case AmpCorPPN:
	    case TaylorT1: 
            case Eccentricity:
	    case TaylorT2:
	    case TaylorT3:
	    case TaylorF1:
	    case TaylorF2:
	    case EOB:
	    case PadeT1:
	    case PadeF1:
	    case SpinTaylor:
	      CreateListfromTmplt(&insptmplt, tmpltCurrent);
	      insptmplt.massChoice = t03;
	      LAL_CALL(LALInspiralParameterCalc( &status,  &(insptmplt) ), &status);
	      
              overlapin.param = insptmplt;
              if( userParam.template == Eccentricity)
              {
	        insptmplt.massChoice = t03;
                overlapin.param.order = 4; /*insptmplt does not contain the good approximant*/              
              }
              else
              {
                overlapin.param.approximant = userParam.template; /*insptmplt does not contain the good approximant*/              
              }
              LAL_CALL(LALInspiralParameterCalc( &status,  &(overlapin.param) ), &status);
	      if( userParam.template == Eccentricity)
              {
                overlapin.param.order = coarseBankIn.order;               
              }
              overlapout.max    = -1.; /* let us be sure that it has a value */
              overlapin.param.fCutoff = coarseBankIn.fUpper;
              /* we need to set the fLower of the template and the SNR integral */
              overlapin.param.fLower = coarseBankIn.fLower;
              overlapin.param.fFinal = randIn.param.tSampling/2. - 1;
              
              if (userParam.faithfulness)
	      {
	       /* same parameters but different templates */
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
	        		
		
	      /* if we want to cut integration before the fFinal*/
   	      /*overlapin.param.fCutoff = 512;*/
	       

	      dt0 = -(randIn.param.t0 - insptmplt.t0);
	      dt3 = -(randIn.param.t3 - insptmplt.t3);
	      g00 = tmpltCurrent->Gamma[3] - 
	           tmpltCurrent->Gamma[1] * tmpltCurrent->Gamma[1]/tmpltCurrent->Gamma[0];
	      g01 = tmpltCurrent->Gamma[4] - 
	           tmpltCurrent->Gamma[1] * tmpltCurrent->Gamma[2]/tmpltCurrent->Gamma[0];
	      g11 = tmpltCurrent->Gamma[5] - 
	           tmpltCurrent->Gamma[2] * tmpltCurrent->Gamma[2]/tmpltCurrent->Gamma[0];	    
	      

	      match = 1 - (g00*dt0*dt0 + 2*g01*dt0*dt3 + g11*dt3*dt3);
	
	      {
	        if (userParam.fastSimulation == 1 && (match <userParam.eMatch )  )
	        { 		  
/*  		  filter_processed--;*/
                  gsl_matrix_set(amb1,2,thisTemplateIndex, 0); 
                  gsl_matrix_set(amb1,3,thisTemplateIndex, 0); 
                  gsl_matrix_set(amb1,0,thisTemplateIndex, insptmplt.t0); 
                  gsl_matrix_set(amb1,1,thisTemplateIndex, insptmplt.t3); 
	        }		
	        else
	        {
                    LAL_CALL(LALInspiralWaveOverlap(&status,
						  &correlation,
						  &overlapout,
						  &overlapin), &status);
                  /* we compute the averaged ambiguity function a t=ta and the averaged maximizaed ambiguity function over time*/
		  gsl_matrix_set(amb1,2,thisTemplateIndex, gsl_matrix_get(amb1,2,thisTemplateIndex) + overlapout.max); 
                  gsl_matrix_set(amb1,3,thisTemplateIndex, gsl_matrix_get(amb1,3,thisTemplateIndex) + correlation.data[randIn.param.nStartPad]); 
                  gsl_matrix_set(amb1,0,thisTemplateIndex, insptmplt.t0); 
                  gsl_matrix_set(amb1,1,thisTemplateIndex, insptmplt.t3); 

		  OverlapOutputThisTemplate.rhoMax = overlapout.max;
		  OverlapOutputThisTemplate.templateNumber = currentTemplate;
		  OverlapOutputThisTemplate.phase = overlapout.phase;
		  OverlapOutputThisTemplate.rhoBin = overlapout.bin;
		  OverlapOutputThisTemplate.freq = overlapin.param.fFinal;
		  OverlapOutputThisTemplate.templateNumber = thisTemplateIndex;
		  insptmplt.fFinal = overlapin.param.fFinal;
		  OverlapOutputThisTemplate.snrAtCoaTime  =  correlation.data[randIn.param.nStartPad];
  		  filter_processed++;
	  
		}
	      }
		
		
		
		
	      
	      
	      break;
          } /*end of the switch*/
	  
	  /* fill histogram of the correlation output. Useful to get a flavour
           * of the output distribution in presence of noise only for
           * instance. */
	  if (userParam.printSNRHisto)
	  {
	    for (i=0; i<(INT4)correlation.length; i++)
	    {
              /* in the unconstraint case, if alphafCut is applied,
               * correlation =-1 if alphaF > 1. Therefore we do not count it
               * in the correlation histogram*/
	      if (correlation.data[i]!=-1)
	      {
                gsl_histogram_increment(histogramNoise, correlation.data[i]);
              }
	    }		    
	  }	

		  
	  /* for the final results using only one event (loudest one) throughtout the bank.*/ 
	  KeepHighestValues(OverlapOutputThisTemplate, 
			    &OverlapOutputBestTemplate,insptmplt);
			    
	  /* We can also keep the best trigger in each template. */

        }/*end of  bank process*/
      }/*end of eccentricity*/

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


      PrintResults(result, randIn, ntrials, userParam.ntrials);  
      if (userParam.printResultXml) {
	  BEPrintResultsXml(coarseBankIn,randIn,userParam,result, 0);	 
	}
      
      
      /*Now we might want given the best template, to reproduce results and keep overlap, wavefom of the signal , template ...*/
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

      if (userParam.template == EOB ) {
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


  /* we save the ambiguity in a file */
  if (userParam.ambiguity)
  {
    CHAR str[512];
    sprintf(str, "ambiguity_%d.dat", userParam.useed);
    
    Foutput=  fopen(str,"w");
    for (j=0; j<sizeBank; j++)
      fprintf(Foutput, "%e %e %e %e\n", 
          gsl_matrix_get(amb1,0,j),
          gsl_matrix_get(amb1,1,j),
          gsl_matrix_get(amb1,2,j)/result.ntrial,
          gsl_matrix_get(amb1,3,j)/result.ntrial
          );
    fclose(Foutput);
    gsl_matrix_free(amb1);
  }
  
  /* free memory */
  while ( templateBank.snglInspiralTable )
  {
    tmpltHead = templateBank.snglInspiralTable;
    templateBank.snglInspiralTable = templateBank.snglInspiralTable->next;
    LALFree( tmpltHead );
  }
  
  /* --- destroy the plans, correlation and signal --- */
  LALDestroyRandomParams(&status, &randParams );
  if (userParam.template == BCV)
  {
    LALFree(powerVector.fm5_3.data);
    LALFree(powerVector.fm2_3.data);
    LALFree(powerVector.fm1_2.data);
    LALFree(powerVector.fm7_6.data);
    LALFree(moments.a11.data);
    LALFree(moments.a21.data);
    LALFree(moments.a22.data);   
  
    if (userParam.template == BCV)
    {  
      LALFree(FilterBCV2.data);
      LALFree(FilterBCV1.data);
    }
  }
  
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
		  OverlapOutputIn *resultBestTemplate,
         InspiralTemplate insptmplt)
{
  if (resultThisTemplate.rhoMax > resultBestTemplate->rhoMax){    
    resultBestTemplate->rhoMax   = resultThisTemplate.rhoMax;
    resultBestTemplate->phase    = resultThisTemplate.phase;
    resultBestTemplate->alpha    = resultThisTemplate.alpha;
    resultBestTemplate->rhoBin   = resultThisTemplate.rhoBin;
    resultBestTemplate->freq     = resultThisTemplate.freq;
    resultBestTemplate->templateNumber = resultThisTemplate.templateNumber;
    resultBestTemplate->snrAtCoaTime   = resultThisTemplate.snrAtCoaTime;
    resultBestTemplate->eccentricity = insptmplt.eccentricity;
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
	  UserParametersIn           userParam
	  )
{
  INT4 templateNumber;
  InspiralTemplate trigger;

  
  INITSTATUS (status, "GetResult", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);
     
  templateNumber = bestOverlap.templateNumber;

  trigger = *list;

  if (userParam.template == BCV ){

    /*    LALInspiralParameterCalc( status->statusPtr,  &trigger );
    CHECKSTATUSPTR(status);							
    */
    result->psi0_inject  = injected.psi0;
    result->psi3_inject  = injected.psi3;     
    result->psi0_trigger  = trigger.psi0;
    result->psi3_trigger  = trigger.psi3;     
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
    result->psi0_trigger  = 0.;
    result->psi3_trigger  = 0.;     
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

  result->eccentricity = bestOverlap.eccentricity;
  
  DETATCHSTATUSPTR(status);
  RETURN (status);

}


void
PrintResults(   ResultIn result, 
                RandomInspiralSignalIn randIn, INT4 n, INT4 N)
{
  fprintf(stdout, " %e %e %e %e ",  
	  result.psi0_trigger, 
	  result.psi3_trigger,
	  randIn.param.psi0,
	  randIn.param.psi3);
  
  fprintf(stdout, "%e %e %e %e %e %e ", 
	  result.tau0_trigger, 
	  result.tau3_trigger,
	  randIn.param.t0,
	  randIn.param.t3, 
          result.eccentricity, 
          randIn.param.eccentricity);
  
  fprintf(stdout, "%e %e   %e %e %e ", 
	  result.fend_trigger, 
	  randIn.param.fFinal,
	  randIn.param.mass1,
	  randIn.param.mass2,
	  randIn.param.startPhase);

  
  fprintf(stdout, " %e %e  %e   %e %d %d %d %d/%d\n",
	  result.rho_final, 
	  result.snrAtCoaTime,
	  result.phase, 
	  result.alphaF,
	  result.bin , 
	  randIn.param.nStartPad, 
	  result.nfast, n,N);

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
	  Vector0, Vector1, Vector2, 				/* idem 				*/
	  rho		= 0,			/* the SNR 				*/
           alpha, alphaFU, fm23, fp23;
  REAL4Vector 
    template, 
    x1, x2, x3, x4,
    phaseV,
    rho1, rho2, rho3;
  REAL4Vector v0, v1, v2;					/* output of the correlations 		*/

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

  /*if (userParam.printBestOverlap && userParam.extraFinalPrinting){*/
    phaseV.data = (REAL4*) LALMalloc(sizeof(REAL4) * phaseV.length);
    template.data = (REAL4*) LALMalloc(sizeof(REAL4) * template.length);
    rho1.data = (REAL4*) LALMalloc(sizeof(REAL4) * rho1.length);
    rho2.data = (REAL4*) LALMalloc(sizeof(REAL4) * rho2.length);
    rho3.data = (REAL4*) LALMalloc(sizeof(REAL4) * rho3.length);
    v1.data = (REAL4*) LALMalloc(sizeof(REAL4) * v1.length);
    v2.data = (REAL4*) LALMalloc(sizeof(REAL4) * v2.length);
    v0.data = (REAL4*) LALMalloc(sizeof(REAL4) * v0.length);
  /*}*/
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
      
      Vector0 = x1_2 + x2_2 + x3_2 + x4_2;
      Vector1 = x1_2 + x3_2 - x2_2 - x4_2;
      Vector2 = 2*(x1.data[i]*x2.data[i] + x3.data[i]*x4.data[i]);
      
      rhoUnconstraint = sqrt((Vector0 + sqrt(Vector1*Vector1+Vector2*Vector2))/2.);

      /* get thetav first in order to constraint alpha_F*/

      thetav =   atan2(Vector2,Vector1);
	
      if(userParam.printBestOverlap && userParam.extraFinalPrinting){
	rho1.data[i] = rhoUnconstraint;
	rho2.data[i] = sqrt((Vector0 + Vector1)/2.);
	rho3.data[i] = sqrt((Vector0+Vector1*cos(2*thetab)+Vector2*sin(2*thetab))/2.);
	v0.data[i]   = Vector0;
	v1.data[i]   = Vector1;
	v2.data[i]   = Vector2;
        phaseV.data[i] = .5*thetav; /*used to compute alpha*/
      }
             

      if (thetab >= 0){
	if ( 0  <= thetav && thetav <= 2 * thetab){	  
	  rhoConstraint = rhoUnconstraint;
	}
	else if (thetab-LAL_PI <= thetav && thetav < 0) {
	  rhoConstraint = sqrt((Vector0 + Vector1)/2.);	  
	}
	else if( (2*thetab  < thetav && thetav<=LAL_PI+1e-4 )
		 || (-LAL_PI-1e-4<=thetav && thetav < -LAL_PI+thetab)){
	  rhoConstraint =sqrt((Vector0+Vector1*cos(2*thetab)+Vector2*sin(2*thetab))/2.);;
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
	  rhoConstraint = sqrt((Vector0 + Vector1)/2.);	  
	}
	else if( (-LAL_PI-1e-4  <= thetav && thetav < 2*thetab ) || 
		 (LAL_PI +thetab <= thetav && thetav <= LAL_PI+1e-4))
        {
	  rhoConstraint =sqrt((Vector0+Vector1*cos(2*thetab)+Vector2*sin(2*thetab))/2.);;
	}
	else 
	  {
	    fprintf(stderr,"must not enter herethetav = %e thetab=%e %e %e %d\n ",
                thetav , thetab, Vector1, Vector2, i);
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
  if (userParam.printBestOverlap  && userParam.extraFinalPrinting)
  {

    /*print overlap, combinaison of the x_i filters*/
    Foutput=fopen("BE_Phase.dat","w");
    for (i=0; i< phaseV.length; i++)
    {
      fprintf(Foutput, "%e\n",  atan2(v2.data[i],v1.data[i]));
    }
    fclose(Foutput);
    
    Foutput=fopen("BE_rho1.dat","w");
    for (i=0; i< phaseV.length; i++)
    {
      fprintf(Foutput, "%e\n", rho1.data[i]);
    }
    fclose(Foutput);
    

  Foutput=fopen("BE_rho2.dat","w");
    for (i=0; i< phaseV.length; i++)
    {
      fprintf(Foutput, "%e\n", rho2.data[i]);
    }
    fclose(Foutput);
    

    Foutput=fopen("BE_rho3.dat","w");
    for (i=0; i< phaseV.length; i++)
    {
      fprintf(Foutput, "%e\n",rho3.data[i]);
    }
    fclose(Foutput);
     
    Foutput=fopen("BE_v0.dat","w");
    for (i=0; i< phaseV.length; i++)
    {
      fprintf(Foutput, "%e\n",v0.data[i]);
    }
    fclose(Foutput);

    Foutput=fopen("BE_v1.dat","w");
    for (i=0; i< phaseV.length; i++)
    {
      fprintf(Foutput, "%e\n", v1.data[i]);
    }
    fclose(Foutput); 

    Foutput=fopen("BE_v2.dat","w");
    for (i=0; i< phaseV.length; i++)
    {
      fprintf(Foutput, "%e\n", v2.data[i]);
    }
    fclose(Foutput);
    
    Foutput=fopen("BE_alpha.dat","w");
    for (i=0; i< correlation->length; i++)
    {
      alpha = -(a22 * tan(phaseV.data[i])) 				/* Compute the final alpha parameter 	*/
	/ (a11 + a21* tan(phaseV.data[i]));
      fprintf(Foutput, "%e \n",alpha*fp23);
    }
    fclose(Foutput);

    /*print overlap, combinaison of the x_i filters*/
    Foutput=fopen("BE_Overlap.dat","w");
    for (i=0; i< correlation->length; i++)
    {
      fprintf(Foutput, "%e\n", fabs(correlation->data[i]));
    }
    fclose(Foutput);
  }


    if (userParam.extraFinalPrinting && userParam.printBestOverlap)
    {
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
    for (i=0; i< correlation->length; i++)
    {
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
    for (i=0; i< correlation->length; i++)
    {
      fprintf(Foutput, "%e %e\n", i/corrin.samplingRate,template.data[i]);
    }
    fclose(Foutput);



  }/*end if printTemplate*/

  /* Free memory */
  LALFree(x1.data);
  LALFree(x2.data);
  LALFree(x3.data);
  LALFree(x4.data);

  /*if (userParam.extraFinalPrinting && userParam.printBestOverlap)
  {*/
    LALFree(phaseV.data);
    LALFree(template.data);
    LALFree(rho1.data);
    LALFree(rho2.data);
    LALFree(rho3.data);
    LALFree(v0.data);
    LALFree(v1.data);
    LALFree(v2.data);
 /* }*/

  DETATCHSTATUSPTR(status);
  RETURN(status);
}



/*  ****************************************************************************
 *  The Overlap function for BCV templates 
 *  ***************************************************************************/

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
  switch (coarseBankIn.approximant)
  {
	case BCV:
	      	fprintf(output, "#psi0Min=%e, psi0Max=%e, psi3Min=%e, psi3Max=%e\n", 
	      	coarseBankIn.psi0Min, coarseBankIn.psi0Max, coarseBankIn.psi3Min, coarseBankIn.psi3Max);
	      	fprintf(output, "#psi0 psi3  totalMass fFinal\n");
		break;
      default:
      		fprintf(output, "#mMin=%e, mMax=%e\n", coarseBankIn.mMin,coarseBankIn.mMax);
		fprintf(output, "#tau0, tau3, mass1, mass2\n");
		break;
    }

  for ( i = 0; i < sizeBank; i++)
    {
      switch (coarseBankIn.approximant)
	{
	case BCV:	  		   
		fprintf(output, "%e %e %e %e\n",
			(*list)[i].params.psi0, (*list)[i].params.psi3,
			(*list)[i].params.totalMass,(*list)[i].params.fFinal);
		break;
	default:
      		fprintf(output, "%e %e %e %e\n", 	
			(*list)[i].params.t0, (*list)[i].params.t3, 
			(*list)[i].params.mass1, (*list)[i].params.mass2);
		break;
	}
	
    }	   
  fprintf(output,"&\n");
  fclose(output);
}



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
  LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME ); \
  LALSnprintf( this_proc_param->param,   LIGOMETA_PARAM_MAX,   "%-20s", ppname); \
  LALSnprintf( this_proc_param->type,    LIGOMETA_TYPE_MAX,    "%-10s",   pptype ); \
  LALSnprintf( this_proc_param->value,   LIGOMETA_VALUE_MAX,   format, ppvalue );


#define ADD_2PROCESS_PARAM( pptype, format, ppname, ppvalue1, ppvalue2 ) \
this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
  calloc( 1, sizeof(ProcessParamsTable) ); \
  LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME ); \
  LALSnprintf( this_proc_param->param,   LIGOMETA_PARAM_MAX,   "%-20s", ppname); \
  LALSnprintf( this_proc_param->type,    LIGOMETA_TYPE_MAX,    "%-10s",   pptype ); \
  LALSnprintf( this_proc_param->value,   LIGOMETA_VALUE_MAX,   format, ppvalue1, ppvalue2 );

  ADD_PROCESS_PARAM("float","%f","--bank-alpha",           
      coarseBankIn.alpha);
  ADD_2PROCESS_PARAM("float","%f %f","--bank-fcut-range",
      coarseBankIn.LowGM, coarseBankIn.HighGM);
  if (coarseBankIn.numFcutTemplates>0)
  {
    ADD_PROCESS_PARAM("float","%f","--bank-number-fcut",
      coarseBankIn.numFcutTemplates);
  }
  ADD_PROCESS_PARAM("float","%d ","--bank-inside-polygon",
      coarseBankIn.insidePolygon);
  ADD_PROCESS_PARAM("float","%f","--bank-ffinal",
      coarseBankIn.fUpper);
  ADD_2PROCESS_PARAM("float","%f %f","--bank-mass-range",
      coarseBankIn.mMin, coarseBankIn.mMax);
  ADD_PROCESS_PARAM("float","%f","--bank-max-total-mass",
      coarseBankIn.MMax);
  ADD_PROCESS_PARAM("float","%f","--bank-min-total-mass",
      coarseBankIn.MMin);
  ADD_2PROCESS_PARAM("float","%f %f","--bank-psi0-range",
      coarseBankIn.psi0Min, coarseBankIn.psi0Max);
  ADD_2PROCESS_PARAM("float","%f %f","--bank-psi3-range",
      coarseBankIn.psi3Min, coarseBankIn.psi3Max);
  ADD_PROCESS_PARAM("string","%s","--bank-grid-spacing",
      GetStringFromGridType(coarseBankIn.gridSpacing));  
  ADD_2PROCESS_PARAM("float","%f","--eccentricity-range",
      userParam.eccentricityMin, userParam.eccentricityMax);
  ADD_PROCESS_PARAM("float","%f","--fl",	
      coarseBankIn.fLower);
  ADD_PROCESS_PARAM("float","%f","--template-fl", 
      coarseBankIn.fLower);
  ADD_PROCESS_PARAM("float","%f","--signal-max-total-mass",
      userParam.maxTotalMass);
  ADD_PROCESS_PARAM("float","%f","--signal-min-total-mass",
      userParam.minTotalMass);
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
  ADD_PROCESS_PARAM("int","%d","--num-seconds",
      userParam.numSeconds);
  ADD_PROCESS_PARAM("float","%f","--psi0",	
      userParam.psi0);
  ADD_PROCESS_PARAM("float","%f","--psi3",
      userParam.psi3);
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
  ADD_PROCESS_PARAM("float","%f","--signal-fl",
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
  ADD_2PROCESS_PARAM("float","%f %f","--t0-fine-range",
      userParam.t0FineMin, userParam.t0FineMax);
  ADD_2PROCESS_PARAM("float","%f %f","--t3-fine-range",
      userParam.t3FineMin, userParam.t3FineMax);
  ADD_PROCESS_PARAM("float","%f","--t0-fine-bin",
      userParam.t0FineBin);
  ADD_PROCESS_PARAM("float","%f","--t3-fine-bin",
      userParam.t3FineBin);
  ADD_PROCESS_PARAM("float","%e","--e-match",	
      userParam.eMatch);
  if (userParam.printResultXml){
    ADD_PROCESS_PARAM("float", "%s", "--print-xml"," ");
  }
  if (coarseBankIn.computeMoments == 1 ){
    ADD_PROCESS_PARAM("int", "%s", "--compute-moments"," ");
  }
  if (userParam.startPhase==0){
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
  if (userParam.eccentricBank ){
    ADD_PROCESS_PARAM("float", "%s", "--eccentric-bank"," ");
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
  static CHAR this[256]; 

  switch (input){
    case 0:
      LALSnprintf(this, sizeof(this), "SignalOnly");
      break;
    case 1:
      LALSnprintf(this, sizeof(this), "NoiseOnly");
      break;
    case 2:
      LALSnprintf(this, sizeof(this), "NoiseAndSignal");
      break;
  }

  return this;
}
CHAR* GetStringFromTemplate(INT4 input)
{
  static CHAR this[256];

  switch (input){
    case EOB:
      LALSnprintf(this, sizeof(this), "EOB");
      break;
    case AmpCorPPN:
      LALSnprintf(this, sizeof(this),"AmpCorPPN");
      break;
    case TaylorT1:
      LALSnprintf(this, sizeof(this),"TaylorT1");
      break;
    case Eccentricity:
      LALSnprintf(this, sizeof(this),"Eccentricity");
      break;
    case TaylorT2:
      LALSnprintf(this,  sizeof(this),"TaylorT2");
      break;
    case TaylorT3:
      LALSnprintf(this, sizeof(this), "TaylorT3");
      break;
    case PadeT1:
      LALSnprintf(this, sizeof(this), "PadeT1");
      break;
    case TaylorF2:
      LALSnprintf(this, sizeof(this), "TaylorF2");
      break;
    case BCV:
      LALSnprintf(this, sizeof(this), "BCV");
      break;
    case SpinTaylor:
      LALSnprintf(this, sizeof(this), "SpinTaylor");
      break;
  }
  return this;
}

CHAR* GetStringFromGridType(INT4 input)
{
  static CHAR this[256]; 

  switch (input){
  case Square:
    LALSnprintf(this, sizeof(this), "Square");
    break;
  case SquareNotOriented:
    LALSnprintf(this, sizeof(this), "SquareNotOriented");
    break;
  case HexagonalNotOriented:
    LALSnprintf(this,  sizeof(this), "HexagonalNotOriented");
    break;
  case Hexagonal:
    LALSnprintf(this, sizeof(this),"Hexagonal");
    break;
  case HybridHexagonal:
    LALSnprintf(this, sizeof(this),"HybridHexagonal");
    break;
  }
  
  return this;
}

CHAR* GetStringFromNoiseModel(INT4 input)
{
  static CHAR this[256];

  switch (input){
  case LIGOI:
    LALSnprintf(this, sizeof(this),"LIGOI");
    break;
  case LIGOA:
    LALSnprintf(this, sizeof(this),"LIGOA");
    break;
  case VIRGO:
    LALSnprintf(this, sizeof(this),"VIRGO");
    break;
  case EGO:
    LALSnprintf(this, sizeof(this),"EGO");
    break;
  case GEO:
    LALSnprintf(this, sizeof(this),"GEO");
    break;
  case TAMA:    
    LALSnprintf(this, sizeof(this),"TAMA");
    break;
  case UNITY:    
    LALSnprintf(this, sizeof(this),"UNITY");
    break;
  case REALPSD:
    LALSnprintf(this, sizeof(this),"REALPSD");
    break;
  } 
  return this;
}



      
/* xml file for the standalone code */
void 
BEPrintResultsXml( InspiralCoarseBankIn         coarseBankIn,
		   RandomInspiralSignalIn       randIn,
		   UserParametersIn             userParam,
		   ResultIn                     trigger,
                   UINT4                        itbest
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

  LALSnprintf( fname, sizeof(fname), 
	       BANKEFFICIENCY_PRINTRESULT_FILEXML ,
	       "simulation",
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
				      CVS_REVISION_C,
				      CVS_SOURCE_C,
				      CVS_DATE_C ), &status );
    this_proc_param = processParamsTable.processParamsTable = 
      (ProcessParamsTable *) calloc( 1, sizeof(ProcessParamsTable) );
    
    BEFillProc(this_proc_param, coarseBankIn, randIn, userParam);
    
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
    
    PRINT_LIGOLW_XML_BANKEFFICIENCY(xmlStream.fp->fp);
  } 
  else
  {
    xmlStream.fp = XLALFileOpenAppend( fname, 0);
  }
  
  fprintf(xmlStream.fp->fp,BANKEFFICIENCY_PARAMS_ROW,
	trigger.psi0_trigger,
	trigger.psi3_trigger,
	randIn.param.psi0, 
	randIn.param.psi3, 
	trigger.tau0_trigger,
	trigger.tau3_trigger,
	randIn.param.t0, 
	randIn.param.t3, 
        trigger.eccentricity,           
        randIn.param.eccentricity,          
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
	trigger.nfast,
        itbest
	);

  if (trigger.ntrial == (UINT4)userParam.ntrials)
  {
    PRINT_LIGOLW_XML_TABLE_FOOTER(xmlStream.fp->fp);
    PRINT_LIGOLW_XML_FOOTER(xmlStream.fp->fp);
  }
  else
  {
    fprintf(xmlStream.fp->fp, ",\n");
  }
    
     
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlStream ), &status );
      
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
				      PROGRAM_NAME, CVS_REVISION_C, CVS_SOURCE_C, CVS_DATE_C ), &status );
    this_proc_param = processParamsTable.processParamsTable = 
      (ProcessParamsTable *) calloc( 1, sizeof(ProcessParamsTable) );
    
    BEFillProc(this_proc_param, coarseBankIn, randIn, userParam);
    
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
    
     fclose( xmlStream.fp->fp );
      xmlStream.fp->fp = NULL;
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
                      UserParametersIn userParam,
		      UINT4 *length)
{
  /* Use these for the longest template from the bank */
  InspiralTemplate params;
  UINT4 maxTmpltLength = 0;

  INITSTATUS( status, "BEGetMaximumSize", BANKEFFICIENCYC );
  ATTATCHSTATUSPTR( status );
  
  *length = 0;
  randIn.param.massChoice 	= m1Andm2;
  if( randIn.param.approximant 	== Eccentricity)
  {
  params.fLower /=1.;
  }
  else 
  {
    /* why ? because this is probably the longest template (stops at light
     * ring) */
    randIn.param.approximant 	= EOB;
  }
  /*HEre we do not care about the order, so let us set it to a value that wont
   * cause any problems (Newtonian , 1PN not properly implemented eveywhere!*/
  coarseBankIn.order = twoPN;
  randIn.param.order = twoPN;

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
  
  /* if --num-seconds is provided, we want a specified length of data.
   * Otherwise, we stick to a minimal size given by twice the length of the
   * maximal template length.*/

  if (userParam.numSeconds * randIn.param.tSampling >= *length) {
    *length = userParam.numSeconds * randIn.param.tSampling;
  }
  else if (userParam.numSeconds != -1) {
    fprintf(stderr, "you asked for a signal duration of  %d seconds (--num-seconds) but a (%f %f) system  might be longer than that this duration (%f)...quitting \n", 
	    userParam.numSeconds, randIn.mMin, randIn.mMin, (float)(*length)/(float)randIn.param.tSampling);
    exit(0);
  }
  else {
    userParam.numSeconds = (INT4)(*length/randIn.param.tSampling);
  }
    
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
  FILE  *Foutput;

  INITSTATUS( status, "BECreatePsd", BANKEFFICIENCYC );
  ATTATCHSTATUSPTR( status );
  
  
  if (vrbflg){
    fprintf(stdout, "generating PSD ...");
  }

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
    case EGO: 
      LAL_CALL(LALNoiseSpectralDensity (status->statusPtr, coarseBankIn->shf.data, &LALEGOPsd, df),
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
      fprintf(stderr, "to be implemented\n");
      exit(0);
      break;
      

    case READPSD:  
      fprintf(stderr, "to be implemented\n");
      exit(0);
      break;
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


  if (vrbflg){
    fprintf(stdout, " ... done \n");
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
  
  
  trial = 0 ;
  success = 0 ;
  
  /* we might force to compute a non random waveform giving the two
     input parameters m1 and m2 */      
  randIn->param.approximant = userParam.signal; 
  if (userParam.startPhase == 1){
    LALUniformDeviate( status->statusPtr, &u, randParams);
    CHECKSTATUSPTR(status);
    randIn->param.startPhase = u * (2*LAL_PI) - LAL_PI;
  }
  else{
    randIn->param.startPhase = 0.;
  }
 
  
  if (userParam.eccentricityMin >=0 && (userParam.eccentricityMax <= 1)){
    LALUniformDeviate( status->statusPtr, &u, randParams);
    CHECKSTATUSPTR(status);
    randIn->param.eccentricity = userParam.eccentricityMin + u * (userParam.eccentricityMax-userParam.eccentricityMin);
  }
  else{
    randIn->param.eccentricity = 0;
  }




  
  if (randIn->type != 1)
  {
    if (userParam.signal == BCV )
    {
      /* user parameters*/
      if (userParam.psi0!=-1 && userParam.psi3!=-1)
      {
	randIn->param.psi0 = userParam.psi0;
	randIn->param.psi3 = userParam.psi3;
	randIn->param.massChoice = fixedPsi;
	fprintf(stderr, "signal = %f %f %f\n", randIn->param.psi0,randIn->param.psi3, randIn->param.alpha1 );
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
	    int trials = 0; 
		  
	    randIn->param.psi0 = randIn->psi0Min + epsilon1 * (randIn->psi0Max - randIn->psi0Min);
	    randIn->param.psi3 = randIn->psi3Min + epsilon2 * (randIn->psi3Max - randIn->psi3Min);
		  
		  
            while (((fend > randIn->param.tSampling/2.) || ( fend < randIn->param.fLower)) && (trials < 10))
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
	
	      trials++;
	    }
	    if (trials == 10 ) valid = 0;
	    else valid = 1;
	  }
	}	    
      }
    }
    else /* EOB , T1 and so on*/
    {
     if (randIn->param.approximant==SpinTaylor)
      {  
        /*Now, we randomize the spin parameters only*/
        INT4 temp;
        randIn->param.distance = 1;
        randIn->param.signalAmplitude = 1;
	    temp = randIn->param.massChoice;
/*        randIn->param.massChoice = spinOnly;*/
        randIn->spin1min = 0;
        randIn->spin1max = 1;
        randIn->spin2min = 0;
        randIn->spin2max = 1;
        randIn->inclinationMin = 0; /*inclination must be >0*/
        randIn->inclinationMax = 1;
        randIn->sourcePhiMin = 0.;
        randIn->sourcePhiMax = 1.;
        randIn->sourceThetaMin = 0.;
        randIn->sourceThetaMax = 1.;
        randIn->polarisationAngleMin = 0.;
        randIn->polarisationAngleMax = 1.;
        randIn->minDistance = 0.1;
        randIn->maxDistance = 1.;
        randIn->param.mass1 = 11; /*dummy but valid values overwritten later*/
        randIn->param.mass1 = 11;

        LALRandomInspiralSignal(status->statusPtr, signal, randIn);
        CHECKSTATUSPTR(status); 
        randIn->param.massChoice = temp;
      } 
    
    
      if (userParam.m1!=-1 && userParam.m2!=-1)
      { 
        randIn->param.mass1 = userParam.m1;
        randIn->param.mass2 = userParam.m2;
        randIn->param.massChoice = fixedMasses;
      }
      else if (userParam.tau0!=-1 && userParam.tau3!=-1) 
      {
        randIn->param.t0 = userParam.tau0;
        randIn->param.t3 = userParam.tau3;
        randIn->param.massChoice = fixedTau;
      }
      else if (userParam.binaryInjection == BHNS)
      {
        randIn->param.massChoice = bhns;
      }
      else if (userParam.minTotalMass > 0 && userParam.maxTotalMass > 0)
      {
        randIn->param.massChoice = minmaxTotalMass;
        randIn->MMax = userParam.maxTotalMass;
        randIn->MMin = userParam.minTotalMass;
      }
      else if (randIn->t0Min != 0)
      {
        randIn->param.massChoice = t03;
      }
      else
      {
        randIn->param.massChoice = m1Andm2;
        randIn->param.massChoice = totalMassUAndEta;	  
        randIn->param.distance = 1;
        randIn->param.signalAmplitude = 1;
      }
      
      /*Here, we  randomize the masses*/  
      LALRandomInspiralSignal(status->statusPtr, signal, randIn);
      CHECKSTATUSPTR(status);      
      

     
       
    
    
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
			     moments), status->statusPtr);
  
  
  DETATCHSTATUSPTR( status );
  RETURN( status );  

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
		  tmpltCurrent->psi0, tmpltCurrent->psi3, tmpltCurrent->f_final);
	} 
      else 
	{	  
	  fprintf(output, "%e %e %e %e\n", 	
		  tmpltCurrent->tau0, tmpltCurrent->tau3, tmpltCurrent->mass1, tmpltCurrent->mass2);
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
  LALSnprintf( fname, sizeof(fname), "TMPLTBANK.xml" ,
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
				    PROGRAM_NAME, CVS_REVISION_C, CVS_SOURCE_C, CVS_DATE_C ), &status );
  this_proc_param = processParamsTable.processParamsTable = 
    (ProcessParamsTable *) calloc( 1, sizeof(ProcessParamsTable) );

  BEFillProc(this_proc_param, coarseBankIn, randIn, userParam);

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
  
  /* close the output xml file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlStream ), &status );






}

void
 CreateListfromTmplt(InspiralTemplate  *params, 
		     SnglInspiralTable *tmplts)
{
  /* Right now only t03 and psi0psi3 bank exists so we 
     only need those information in principle and nothing else*/
  params->mass1 = tmplts->mass1; 
  params->mass2 = tmplts->mass2; 
  params->fFinal = tmplts->f_final; 
  params->t0 = tmplts->tau0;
  params->t3 = tmplts->tau3;
  params->psi0 = tmplts->psi0;
  params->psi3 = tmplts->psi3;
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
  coarseBankIn->MMax      	= -1;
  coarseBankIn->massRange       = MinMaxComponentMass; 
  coarseBankIn->etamin 	        = -1;
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
  coarseBankIn->computeMoments  = 0;
  coarseBankIn->NumFreqCut = 1;
  coarseBankIn->MaxFreqCut = SchwarzISCO;
  coarseBankIn->MinFreqCut = SchwarzISCO;
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
  randIn->t0Min                 = 0;                                           /* min tau0 to inject                   */
  randIn->t0Max                 = 0.;                                           /* max tau0 to inject                   */
  randIn->tnMin                 = 0.1;                                           /* min tau3 to inject                   */
  randIn->tnMax                 = 1.;                                           /* max tau3 to inject                   */
  randIn->etaMin                = BANKEFFICIENCY_MMIN
    * (BANKEFFICIENCY_MMAX - BANKEFFICIENCY_MMIN) 
    / (BANKEFFICIENCY_MMAX * 2) / (BANKEFFICIENCY_MMAX * 2);
  randIn->psi0Min  		= BANKEFFICIENCY_PSI0MIN;			/* psi0 range 				*/
  randIn->psi0Max  		= BANKEFFICIENCY_PSI0MAX;	
  randIn->psi3Min  		= BANKEFFICIENCY_PSI3MIN;	
  randIn->psi3Max  		= BANKEFFICIENCY_PSI3MAX;	
  randIn->param.approximant 	= BANKEFFICIENCY_SIGNAL;			/* approximant of h, the signal		*/
  randIn->param.tSampling 	= BANKEFFICIENCY_TSAMPLING;			/* sampling 				*/
  randIn->param.fCutoff 	= 0.; 			/*will be populate later once sampling is well defined*/
  randIn->param.startTime       = BANKEFFICIENCY_STARTTIME; 
  randIn->param.startPhase      = BANKEFFICIENCY_STARTPHASE; 
  randIn->param.nStartPad       = BANKEFFICIENCY_NSTARTPAD;
  randIn->param.signalAmplitude = BANKEFFICIENCY_SIGNALAMPLITUDE;
  randIn->param.nEndPad         = BANKEFFICIENCY_NENDPAD;
  randIn->NoiseAmp              = BANKEFFICIENCY_NOISEAMPLITUDE;
  randIn->param.eccentricity    = 0.;
}


/* ****************************************************************************
 * Other Param initialization							
 * ***************************************************************************/
void InitUserParametersIn(UserParametersIn *userParam)
{
  userParam->eccentricityMin   = 0.;
  userParam->eccentricityMax   = 0.;
  userParam->alphaFConstraint   = ALPHAFConstraint;
  userParam->extraFinalPrinting = 0; 
  userParam->eMatch = 0.5; 
  userParam->template		= BANKEFFICIENCY_TEMPLATE;
  /*By default, this value is empty and will be populate with the sampling frequency later*/
  userParam->signalfFinal       =  0.;
  userParam->signal       	= BANKEFFICIENCY_SIGNAL;
  userParam->m1           	= -1;
  userParam->m2           	= -1;
  userParam->numSeconds         = -1;
  userParam->psi0         	= -1;
  userParam->psi3         	= -1;
  userParam->tau0         	= -1;
  userParam->tau3         	= -1;
  userParam->t0FineMin         	= 0;
  userParam->t0FineMax         	= 0;
  userParam->t3FineMin         	= 0;
  userParam->t3FineMax         	= 0;
  userParam->t0FineBin         	= 0;
  userParam->t3FineBin         	= 0;
  userParam->printBestOverlap 	= BANKEFFICIENCY_PRINTBESTOVERLAP;
  userParam->printBestTemplate 	= BANKEFFICIENCY_PRINTBESTTEMPLATE;
  userParam->printSNRHisto      = BANKEFFICIENCY_PRINTSNRHISTO;
  userParam->printPsd           = BANKEFFICIENCY_PRINTPSD;
  userParam->printBank    	= BANKEFFICIENCY_PRINTBANK;
  userParam->printResultXml    	= BANKEFFICIENCY_PRINTRESULTXML;
  userParam->printPrototype    	= BANKEFFICIENCY_PRINTPROTOTYPE;
  userParam->faithfulness       = BANKEFFICIENCY_FAITHFULNESS;
  userParam->ntrials 	     	= BANKEFFICIENCY_NTRIALS;
  userParam->fastSimulation     = BANKEFFICIENCY_FASTSIMULATION;
  userParam->noiseModel         = LIGOI;
  userParam->binaryInjection    = NoUserChoice; 
  userParam->maxTotalMass       = -1;
  userParam->startPhase         = 1;
  userParam->numSeconds   	= -1;
  userParam->inputPSD     	= NULL;
  userParam->ambiguity     	= 0;
  userParam->eccentricBank     	= 0;
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

  if (*argc==1) 
  {
    Help();
    fprintf(stderr, "You must provide at least the template and signal argument");
    exit(0);
    }
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
	else if (!strcmp(argv[i], "HybridHexagonal")) { 
	  coarseBankIn->gridSpacing = HybridHexagonal;
	}
	else if (!strcmp(argv[i], "HexagonalNotOriented")) {
	  coarseBankIn->gridSpacing = HexagonalNotOriented;
	}
	else if (!strcmp(argv[i], "SquareNotOriented")) {
	  coarseBankIn->gridSpacing = SquareNotOriented;
	}
	else {
	  fprintf(stderr, "Unknown grid spacing after --bank-grid-spacing. Try --help.\n");
	  exit(0);
	}
      }        
      else if (!strcmp(argv[i],	"--bank-number-fcut")) {
        BEParseGetInt(argv,  &i, (INT4*)&(coarseBankIn->numFcutTemplates)); 
      }
      else if (!strcmp(argv[i],"--bank-mass-range")) {
	BEParseGetDouble2(argv,  &i, &(coarseBankIn->mMin), &(coarseBankIn->mMax));
      }
      else if (!strcmp(argv[i],"--bank-max-total-mass")) {
	BEParseGetDouble(argv,  &i, &(coarseBankIn->MMax));
      }
      else if (!strcmp(argv[i],"--bank-min-total-mass")) {
	BEParseGetDouble(argv,  &i, &(coarseBankIn->MMin));
      }
      else if (!strcmp(argv[i],	"--bank-psi0-range")) {
        BEParseGetDouble2(argv,  &i, &(coarseBankIn->psi0Min), &(coarseBankIn->psi0Max));
      }
      else if (!strcmp(argv[i],	"--bank-psi3-range")) {
        BEParseGetDouble2(argv,  &i, &(coarseBankIn->psi3Min), &(coarseBankIn->psi3Max));
      }
      else if (!strcmp(argv[i],"--debug")) {
        BEParseGetInt(argv,  &i, &(lalDebugLevel)); 
      }
      else if (!strcmp(argv[i],"--eccentricity-range")){
        BEParseGetDouble2(argv,  &i, &(userParam->eccentricityMin), &(userParam->eccentricityMax));
      }
      else if (!strcmp(argv[i], "--e-match")){
        BEParseGetDouble(argv,  &i, &(userParam->eMatch));     
      }
      else if (!strcmp(argv[i],"--signal-fl")) {
        BEParseGetDouble(argv,  &i, &(randIn->param.fLower));
      }
      else if  (!strcmp(argv[i],"--template-fl")) {
        BEParseGetDouble(argv,  &i, &(coarseBankIn->fLower));
      }
      else if  (!strcmp(argv[i],"--fl")) {
        BEParseGetDouble(argv,  &i, &(coarseBankIn->fLower));
	randIn->param.fLower = coarseBankIn->fLower;
      }  
      else if (!strcmp(argv[i], "--help") || !strcmp(argv[i],"--h")) {
	Help();
      }
      else if (!strcmp(argv[i],"--signal-max-total-mass")) {
	BEParseGetDouble(argv,  &i, &tmp1);
	userParam->maxTotalMass = (REAL4)tmp1;	
      }
      else if (!strcmp(argv[i],"--signal-min-total-mass")) {
	BEParseGetDouble(argv,  &i, &tmp1);
	userParam->minTotalMass = (REAL4)tmp1;	
      }
      else if (!strcmp(argv[i],"--ascii2xml")) {       
	ascii2xml = 1;
        /* we do not need to parse the other arguments */
        i = *argc+1;
      }
      else if (!strcmp(argv[i], "--m1")){
        BEParseGetDouble(argv,  &i, &(userParam->m1));
      }
      else if (!strcmp(argv[i], "--m2")) {
        BEParseGetDouble(argv,  &i, &(userParam->m2));
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
	else if (!strcmp(argv[i], "EGO")) {
	  userParam->noiseModel = EGO;
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
	  fprintf(stderr, "Unknown moise model after --noise-model. Try --help\n");
   	  exit(0);
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
      else if (!strcmp(argv[i],"--sampling")) {
        BEParseGetDouble(argv,  &i, &(coarseBankIn->tSampling));
        randIn->param.tSampling = coarseBankIn->tSampling;
      }   
      else if (!strcmp(argv[i],"--seed")){
        BEParseGetInt(argv,  &i, &(randIn->useed));
      }
      else if (!strcmp(argv[i], "--signal")) {
	BEParseGetString(argv, &i);
	
	if (strcmp(argv[i],	"TaylorT1")	==0) userParam->signal = TaylorT1;
	else if (strcmp(argv[i],"AmpCorPPN")	==0) userParam->signal = AmpCorPPN;
	else if (strcmp(argv[i],"Eccentricity")	==0) userParam->signal = Eccentricity;
    else if (strcmp(argv[i],"TaylorT2")	==0) userParam->signal = TaylorT2;
	else if (strcmp(argv[i],"TaylorT3")	==0) userParam->signal = TaylorT3;
	else if (strcmp(argv[i],"TaylorF1")	==0) userParam->signal = TaylorF1;
	else if (strcmp(argv[i],"TaylorF2")	==0) userParam->signal = TaylorF2;
	else if (strcmp(argv[i],"PadeT1")	==0) userParam->signal = PadeT1;
	else if (strcmp(argv[i],"PadeF1")	==0) userParam->signal = PadeF1;
	else if (strcmp(argv[i],"EOB")		==0) userParam->signal = EOB;
	else if (strcmp(argv[i],"BCV")		==0) userParam->signal = BCV;
	else if (strcmp(argv[i],"SpinTaylor")   ==0) userParam->signal = SpinTaylor;
	else  
	{
	  fprintf(stderr, "Wrong approximant. Use of of TaylorT1, TaylorT2, TaylorT3, EOB, PadeT1, BCV, SpinTaylor, AmpCorPPN, Eccentricity)\n");
	  exit( 1 );
	}
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
	else 
	{
	fprintf(stderr, "unknow simulation-type after --simulation-type. try --help \n");
	exit(0);
	}
      }
      else if (!strcmp(argv[i], "--t0-fine-range")){
        BEParseGetDouble2(argv,  &i, &(userParam->t0FineMin), &(userParam->t0FineMax));
      }
      else if (!strcmp(argv[i], "--t3-fine-range")){
        BEParseGetDouble2(argv,  &i, &(userParam->t3FineMin), &(userParam->t3FineMax));
      }
      else if (!strcmp(argv[i], "--t0-fine-bin")){
        BEParseGetDouble(argv,  &i, &(userParam->t0FineBin));
      }
      else if (!strcmp(argv[i], "--t3-fine-bin")){
        BEParseGetDouble(argv,  &i, &(userParam->t3FineBin));
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
	
	if (!strcmp(argv[i],	"TaylorT1")	)		userParam->template = TaylorT1;
	else if (!strcmp(argv[i],	"TaylorT2")	)	userParam->template = TaylorT2;
        else if (!strcmp(argv[i],	"Eccentricity")	)	userParam->template = Eccentricity;
        else if (!strcmp(argv[i],	"AmpCorPPN")	)	userParam->template = AmpCorPPN;
	else if (!strcmp(argv[i],	"TaylorT3")	)	userParam->template = TaylorT3;
	else if (!strcmp(argv[i],	"TaylorF1")	)	userParam->template = TaylorF1;
	else if (!strcmp(argv[i],	"TaylorF2")	)	userParam->template = TaylorF2;
	else if (!strcmp(argv[i],	"PadeT1")	)	userParam->template = PadeT1;
	else if (!strcmp(argv[i],	"PadeF1")	)	userParam->template = PadeF1;
	else if (!strcmp(argv[i],	"EOB")		)	userParam->template = EOB;
	else if (!strcmp(argv[i],	"BCV")		)       userParam->template = BCV;
	else 
	  {
      	    fprintf(stderr, "Wrong approximant. Use of of TaylorT1, TaylorT2, TaylorT3, EOB, PadeT1, BCV, AmpCorPPN, Eccentricity)\n");
	    exit( 1 );
	  }

	
	coarseBankIn->approximant = userParam->template;
	if ( coarseBankIn->approximant == BCV ) 	
	  coarseBankIn->space = Psi0Psi3; 
	else 
	  coarseBankIn->space = Tau0Tau3;
      }
      else if (!strcmp(argv[i],	"--template-order")) {
        BEParseGetInt(argv,  &i, (INT4*)(&(coarseBankIn->order)));
      }      
      /* no user parameters requeste, only flags below*/               
      else if (!strcmp(argv[i],"--ambiguity")) {
        userParam->ambiguity = 1;
      }
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
	fprintf(stderr, "BankEfficiency code"
		"Thomas Cokelaer, Thomas.Cokelaer@astro.cf.ac.uk\n"
		"CVS Version source:" CVS_ID_STRING_C "\n"
		"CVS Version include:" CVS_ID_STRING "\n"
		"CVS Tag: " CVS_NAME_STRING "\n");
	exit(0);
      }
      else if (!strcmp(argv[i],"--print-bank")) {
        userParam->printBank		= 1;
      }
      else if (!strcmp(argv[i],"--eccentric-bank")) {
        userParam->eccentricBank		= 1;
      }
      else if (!strcmp(argv[i],"--compute-moments")) {
        coarseBankIn->computeMoments = 1;
      }
      else if (!strcmp(argv[i],"--print-xml")) {
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
	  fprintf(stderr, msg);
  	  exit( 1 ); 
	}
    }
  else  
    {
      sprintf(msg, "Expect a float after option %s (got %s)\n ",
	      argv[*index],
	      argv[*index+1]);
      fprintf(stderr, msg);
      exit( 1 ); 
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
      fprintf(stderr, msg);
      exit( 1 );
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
	  fprintf(stderr, msg);
	}
    }
  else  
    {
      sprintf(msg, "Expect a float after option %s (got %s)\n ",
	      argv[*index],
	      argv[*index+1]);
      fprintf(stderr, msg);
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
          fprintf(stderr,msg);
	}
    }
  else  
    {
      sprintf(msg, "Expect 2 floats after option %s (got %s and %s)\n ",
	      argv[*index],
	      argv[*index+1],argv[*index+2]);
      fprintf(stderr, msg);
    }
  *index = *index +2 ;  
}  
  
    
    /* --- Function to check validity of some parsed parameters (TODO)--- */
void UpdateParams(InspiralCoarseBankIn *coarseBankIn,
		  RandomInspiralSignalIn *randIn,
		  UserParametersIn *userParam)
{
  CHAR   msg[2048];

  if (ascii2xml == 1) 
    return;
  
  userParam->useed = randIn->useed;

  /* -- If fupper is not set, then it is Nyquist frequency*/
  if (coarseBankIn->fUpper == -1 ){
    coarseBankIn->fUpper = coarseBankIn->tSampling/2. -1.;
  }
  /* -- If it is set, it must be below Nyquist*/
  if (coarseBankIn->fUpper >=coarseBankIn->tSampling/2. -1.){
    coarseBankIn->fUpper = coarseBankIn->tSampling/2 -1.;
  }
 
  /* -- the same for the fCutoff */
  if (randIn->param.fCutoff == 0){
      randIn->param.fCutoff = coarseBankIn->tSampling/2 -1.;
  }
  if (randIn->param.fCutoff >=coarseBankIn->tSampling/2. -1.)
  {
    randIn->param.fCutoff = coarseBankIn->tSampling/2 -1.;
  }

  /* -- Alpha_B for the bank must be positive*/
  if (coarseBankIn->alpha < 0 )
  {
    sprintf(msg, "--bank-alpha (%f) parameter must be positive in the range [0,1] \n",
	    coarseBankIn->alpha);
    fprintf(stderr, msg);
  }

  /* -- fLower must be below the cutoff frequencies */
  if (coarseBankIn->fUpper <= coarseBankIn->fLower 
      || coarseBankIn->fUpper >= coarseBankIn->tSampling/2)
  {
    sprintf(msg, "--bank-ffinal (%f) paramter must be greater than  bank-fl (%f) and less than sampling/2 %f\n",
	    coarseBankIn->fUpper,
	    coarseBankIn->fLower ,
	    coarseBankIn->tSampling/2); 
    fprintf(stderr,  msg);
  }
  
  
  /* -- mMin must be greater than mMax and positive */
  if  ((coarseBankIn->mMin >= coarseBankIn->mMax ) || (coarseBankIn->mMin <=0)){
    sprintf(msg, "--bank-mass-range (%f %f) parameter must be sorted and > 0 \n",
	    coarseBankIn->mMin, coarseBankIn->mMax);
    fprintf(stderr, msg);
  }
  
  /* -- the total mass is set for the bank, let us adjust the parameters*/
  if (coarseBankIn->MMax != -1)
  {
    coarseBankIn->mMax = coarseBankIn->MMax - coarseBankIn->mMin;
    if ( coarseBankIn->MMin != -1)
    {
      coarseBankIn->massRange = MinMaxComponentTotalMass;
    }
    else
    {
      coarseBankIn->massRange = MinMaxComponentMass;
    }
  }
  else
  {
    coarseBankIn->MMax = 2 * coarseBankIn->mMax;
    coarseBankIn->MMin = 2 * coarseBankIn->mMin;

    coarseBankIn->massRange = MinMaxComponentMass;
  }


       
  
  if (coarseBankIn->psi0Min <=0 || coarseBankIn->psi0Min > coarseBankIn->psi0Max) {
    sprintf(msg, "--bank-psi0-range (%f %f) paramter must be sorted and > 0 \n",
	    coarseBankIn->psi0Min, coarseBankIn->psi0Max);
    fprintf(stderr, msg);
  } 

  if (coarseBankIn->psi3Min >=0 || coarseBankIn->psi3Min > coarseBankIn->psi3Max) {
    sprintf(msg, "--bank-psi0-range (%f %f) paramter must be sorted and >= 0 \n",
	    coarseBankIn->psi3Min, coarseBankIn->psi3Max);
    fprintf(stderr, msg);
  }


  if (coarseBankIn->LowGM > coarseBankIn->HighGM)
    {
      sprintf(msg, "--bank-fcut-range (%f %f) expect sorted , typically 3 and 6",
	      coarseBankIn->LowGM , coarseBankIn->HighGM);
      fprintf(stderr, msg);
      exit(1);
    }
 
  if (coarseBankIn->fLower <10 || randIn->param.fLower < 10)
    {
      sprintf(msg, "--fl or --fl-signal or --fl-template must be >=10 Hz (%f %f)",
	      randIn->param.fLower , coarseBankIn->fLower);
      fprintf(stderr, msg);
      exit(1);
    }


  if (userParam->maxTotalMass != -1)
    {
      if (userParam->maxTotalMass < 2*randIn->mMin) 
	{	
	  sprintf(msg, "--signal-max-total-mass (%f) must be > twice the minimal mass (%f) ",
	      userParam->maxTotalMass , randIn->mMin );
          fprintf(stderr, msg);
          exit(1);
	}
    }
  if (userParam->minTotalMass != -1)
    {
      if (userParam->minTotalMass > 2*randIn->mMax) 
	{	
	  sprintf(msg, "--min-total-mass (%f) must be < twice the maximal mass (%f) ",
	      userParam->maxTotalMass , randIn->mMax );
          fprintf(stderr, msg);
          exit(1);
        }
    }



  if ((userParam->m1 != -1) && (userParam->m2 != -1)){
    randIn->param.mass1 = userParam->m1;
    randIn->param.mass2 = userParam->m2;
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
      fprintf(stderr, msg);
      exit(1);
    } 
  }
     
  if (userParam->psi0 != -1 && userParam->psi3 != -1){
  
   if (userParam->m1 != -1 ||userParam->m2 != -1 || userParam->tau0 != -1 || userParam->tau3 != -1){
     sprintf(msg, "--m1 --m2 --psi0 --psi3 --tau0 --tau3 error. If particular injection is requested,  you must choose either (--m1,--m2) options or (--psi0,--psi3) or (--tau0,--tau3)\n");
     fprintf(stderr, msg);
     exit(1);
   } 
  } 
  
 if (userParam->tau0 != -1 && userParam->tau3 != -1){   
   if (userParam->psi0 != -1 ||userParam->psi3 != -1 || userParam->m1 != -1 || userParam->m2 != -1){
     sprintf(msg, "--m1 --m2 --psi0 --psi3 --tau0 --tau3 error. If particular injection is requested,  you must choose either (--m1,--m2) options or (--psi0,--psi3) or (--tau0,--tau3)\n");
     fprintf(stderr, msg);  
     exit(1);
   } 
  } 


 
 if (coarseBankIn->mmCoarse <=0 || coarseBankIn->mmCoarse>=1){
   sprintf(msg, "--mm (%f) must be in the range ]0 1[\n", 
	   coarseBankIn->mmCoarse);
     fprintf(stderr, msg);  
     exit(1);
 } 

 /* -- BCV related */
 if (coarseBankIn->numFcutTemplates <= 0){
   sprintf(msg, "--bank-number-fcut (%d) must be > 0>\n", coarseBankIn->numFcutTemplates);
   fprintf(stderr, msg);  
   exit(1);
 } 
 

 if  ((randIn->mMin >= randIn->mMax ) || (randIn->mMin <=0)){
    sprintf(msg, "--signal-mass-range (%f %f) paramter must be sorted and > 0 \n",
	    randIn->mMin, randIn->mMax);
   fprintf(stderr, msg);  
   exit(1);
 }
 
 /* -- MMAX and MMIn must be redefined given the input parameters*/
 randIn->MMax = randIn->mMax * 2; 
 randIn->MMin = randIn->mMin * 2; 

 randIn->etaMin = randIn->mMin * (randIn->MMax - randIn->mMin) / 
 	randIn->MMax / randIn->MMax;
 
 /* -- However, if MMAx and MMIN are also defined by the user*/
 
 if (userParam->maxTotalMass != -1 && userParam->maxTotalMass <2*randIn->mMax)
  {
    randIn->MMax = userParam->maxTotalMass;
      
    if (userParam->minTotalMass != -1 && userParam->minTotalMass >2*randIn->mMin)
    {
       randIn->MMin = userParam->minTotalMass;
    }
  }
	
 
  if ((randIn->t0Min != 0 || randIn->t0Max != 0) && 
        ((randIn->t0Min >= randIn->t0Max ) || (randIn->t0Min <=0))){
    sprintf(msg, "--signal-tau0-range (%f %f) parameter must be sorted and > 0 \n",
            randIn->t0Min, randIn->t0Max);
    fprintf(stderr, msg);
    exit(1);
  }

  if ((randIn->tnMin != 0 || randIn->tnMax != 0) && 
        ((randIn->tnMin >= randIn->tnMax ) || (randIn->tnMin <=0))){
    sprintf(msg, "--signal-tau3-range (%f %f) paramter must be sorted and > 0 \n",
            randIn->tnMin, randIn->tnMax);
    fprintf(stderr, msg);
    exit(1);
  }

  if (userParam->faithfulness ==1 && randIn->type == 1)
    {
      fprintf(stderr, "No injection performed so the check option can not be used (change simulation-type option)\n");
      exit ( 1 );
    }
  if (userParam->binaryInjection == BHNS){
    if (randIn->mMin >3 || randIn->mMax < 3 ){
      fprintf(stderr, "BH-NS injection request to have the minimum mass less than 3 and maximum greater than 3.\n");
      exit( 1 );
    }
  }  

  if (userParam->signalfFinal>0 && userParam->signalfFinal< coarseBankIn->tSampling/2)
  {
    randIn->param.fCutoff = userParam->signalfFinal;
  }
  if (userParam->signalfFinal>0 && userParam->signalfFinal>=coarseBankIn->tSampling/2)
  {
    fprintf(stderr, "--signal-ffinal must be less than the sampling/2. Replaced it with %f\n", coarseBankIn->tSampling/2);
     exit(0)  ;
  }
  if (userParam->signalfFinal==0)
  userParam->signalfFinal = randIn->param.fCutoff ;
  

}




/* ****************************************************************************
 *  Documenation  on line
 *  **************************************************************************/
void Help(void)
{
  CHAR msg[2048];


  fprintf(stderr, "[NAME %s]\n ", CVS_NAME_STRING);
  fprintf(stderr, "[VERSION %s]\n ", CVS_ID_STRING);
  fprintf(stderr, "[VERSION %s]\n ", CVS_ID_STRING_C);
  fprintf(stderr, "[DESCRIPTION]\n");
  sprintf(msg, "\t lalapps_BankEfficiency is a standalone code that can be used to test the efficiency of\n"
	  "\t an inpiral template bank, in the framework of compact binary coalescence searches.\n"
	  "\t This code aimed at providing a test code to the template bank that are used by the search\n"
	  "\t pipeline\n\n"
	  "\t For the LAL fans, you should know that BankEfficiency uses the same functions as\n"
	  "\t lalapps_tmpltbank uses. However, the filtering is independant of the findchirp package.\n"
	  "\t It uses the noisemodel and inspiral packages instead.\n\n");
  fprintf(stderr, "%s",msg);
  sprintf(msg, 
	  "\t BankEfficiency allows to test the template banks in 3 different ways:\n"  
	  "\t  * injection only, in other words, to compute matches \n"
	  "\t  * noise only, to get the background distribution\n"
	  "\t  * noise and signal \n\n"
	  "\t You can use any design sensitivity curve provided in lal/noisemodesl such as LIGOI, VIRGO ...\n\n"
	  "\t Signal to be injected, or used for filtering are either physical template families (e.g., TaylorT3) \n"
	  "\t or phenomenological (BCV).\n");
  fprintf(stderr, "%s",msg);
  sprintf(msg, "\t Concerning computation, there is a naive fast option for physical template families, that sets \n"
  	  "\t a square of 0.3 seconds in tau_0 and 0.15 seconds in tau_3 around the injection\n\n"
	  "\t The code can be easily parallelised. See bep.py for an example of dag and sub condor file\n\n"
	  "\t A brief overview is given here. However, look at the option below, for more details. \n"
	  "\t Testing the bank requires to inject many signals. The number of trial is set with --n.\n" );
  fprintf(stderr, "%s",msg);
  sprintf(msg,"\t The --template and --signal options must be provided to specify the type of approximant to be used. The\n"
	  "\t randomization is done over the total mass, and within the --signal-mass-range and --bank-mass-range ranges.\n"
	  "\t By default xml results are not stored (although is could be a default option, usage in parallel would create\n"
	  "\t many xml files to be joined afterwards). The --check option replaces the template bank by a single template\n"
	  "\t which parameters are identical to the injection.\n\n\n");	  
  fprintf(stderr, "%s",msg);
  sprintf(msg,	"\t The parameter are randomized on the initial phase (can be removed using --no-start-pahse), and the masses\n"
		"\t By default the massesare randomized such  that total mass is uniform. One can set --signal-tau0-range\n"
		"\t which makes the injections uniformly distributed in the tau0/tau3 plane.\n");

  fprintf(stderr, "%s",msg);
  
  sprintf(msg, "[EXAMPLE]\n \t ./lalapps_BankEfficiency --template EOB --signal EOB --check --print-xml\n\n\n");
  fprintf(stderr, "%s",msg);


  
  fprintf(stderr, "[SYNOPSIS]\n");
  fprintf(stderr, "\t[--help]\n");
  fprintf(stderr, "\t[--verbose]  \t\t\t extra information on the screen \n" );
  fprintf(stderr, "\t[--ascii2xml]\t\t\t read the file Trigger.dat, concatenation of one or several outputs\n");
  fprintf(stderr, "\t             \t\t\t BankEfficiency and creates a unique XML outputs.\n");
  fprintf(stderr, "\t[--bank-alpha<float>]\t\t set the BCV alpha value in the moments computation\n");
  fprintf(stderr, "\t[--bank-fcut-range<float float>] set the range of BCV fcut (in units of GM) \n");
  fprintf(stderr, "\t[--bank-ffinal<float>]\t\t set the final frequency to be used in the BCV moments computation\n");
  fprintf(stderr, "\t[--bank-grid-spacing <gridSpacing>]\t set the grid type of the BCV bank (Square, SquareNotOriented,\n\t\t HybridHexagonal,Hexagonal, HexagonalNotOriented\t\n");
  fprintf(stderr, "\t[--bank-number-fcut<integer>]\t number of BCV fcut layers\n");
  fprintf(stderr, "\t[--bank-mass-range<float float>] mass range to be covered by the SPA bank\n");
  fprintf(stderr, "\t[--bank-max-total-mass<float>]  max total mass of the bank\n");
  fprintf(stderr, "\t[--bank-min-total-mass<float>]  min total mass of the bank\n");
  fprintf(stderr, "\t[--bank-psi0-range<float float>] psi_0 range to be covered by the BCV bank\n");
  fprintf(stderr, "\t[--bank-psi3-range<float float>] psi_3 to be covered by the BCV bank\n");
  fprintf(stderr, "\t[--debug<integer>]\t\t set the debug level (same as in lal)\n");
  fprintf(stderr, "\t[--e-match<float>]\t\t set the e-match for the fast simulation\n");
  fprintf(stderr, "\t[--eccentricity-range<float, float>]\t\t set the e-match for the fast simulation\n");
  fprintf(stderr, "\t[--fl-signal<float>]\t\t set the lower cut off frequency of signal to inject\n");
  fprintf(stderr, "\t[--fl-template<float>]\t\t set the lower cut off frequnecy of template \n");
  fprintf(stderr, "\t[--fl<float>]\t\t\t set both template and signal lower cutoff frequency \n");
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
  fprintf(stderr, "\t[--sampling<float>]\t\t set sampling frequency (4096 by default).\n");
  fprintf(stderr, "\t[--seed<integer>]\t\t set seed for random generator.\n");
  fprintf(stderr, "\t[--signal<string>]\t\t set signal approximant (TaylorT1, TaylorT2, TaylorT3, TaylorF2, PadeT1, EOB, SpinTaylor)\n");
  fprintf(stderr, "\t[--signal-alpha<float>]\t\t set alpha parameter of BCV injection\n");
  fprintf(stderr, "\t[--signal-amplitude<float>]\t set SNR of injection in the case NoiseandSignal simulation\n");
  fprintf(stderr, "\t[--signal-max-total-mass<float>]\t set maximum total mass to be injected !! ");
  fprintf(stderr, "\t				\t does not work with spin injection\n");
  fprintf(stderr, "\t[--signal-min-total-mass<float>]\t set minimum total mass to be injected\n");
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
  fprintf(stderr, "\t[--template<string>]\t\t set signal approximant (TaylorT1, TaylorT2, TaylorT3, TaylorF2, PadeT1, EOB, SpinTaylor)\n");
  fprintf(stderr, "\t[--template-order<integer>]\t set PN order of template\n");
  fprintf(stderr, "\t[--alpha-constraint]\t\t set BCV code to be constrained \n");
  fprintf(stderr, "\t[--bhns-injection]\t\t set injection to be only bhbs systems\n");
  fprintf(stderr, "\t[--no-alpha-constraint]\t\t set BCV code to be unconstrained\n");
  fprintf(stderr, "\t[--faithfulness]\t\t Check the code. Template parameters are equal to injection  parameters\n");
  fprintf(stderr, "\t                \t\t size of the bank is therefore unity. It computed the faithfulness instead of effectualness\n");
  fprintf(stderr, "\t[--real-noise]\t\t\t use real data and real PSD.force simulaion type to be Noise Only\n");
  fprintf(stderr, "\t[--no-start-phase]\t\t\t unset random phase which is always set to zero.\n");
  fprintf(stderr, "\t[--print-psd]\t\t\t print the psd in  a file BE_PSD_type_gpstime.dat\n");
  fprintf(stderr, "\t[--print-best-overlap]\t\t print best overlap and other information\n");
  fprintf(stderr, "\t[--print-snr-histo]\t\t print histogram of the correlation output\n");
  fprintf(stderr, "\t[--print-bank]\t\t\t print the bank in ascii and xml format\n");
  fprintf(stderr, "\t[--print-xml]\t\t\t print the results in an xml file\n");
  fprintf(stderr, "\t[--print-prototype]\t\t print a prototype to be used by condor script\n");
  fprintf(stderr, "\t[--fast-simulation]\t\t perform fast simulation in the case of SPA abnk\n");
  exit(0);
}


void 
BEAscii2Xml(void)
{
  UINT4 countline = 0, nfast_max=0;
  UINT8  id = 0;
  ResultIn trigger;
  REAL4 tau0, tau3, tau0I, tau3I, psi0, psi3,phaseI, ecc, eccI
; 
  FILE *input1, *input2, *bank;
  FILE *output;     

  SnglInspiralTable     *inputData = NULL;
  INT4 numFileTriggers = 0, nStartPad;

  char sbuf[2048];
  
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
      while(fgets(sbuf, 2048, input2) !=NULL)
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
      LALSnglInspiralTableFromLIGOLw( &inputData, BEASCII2XML_BANK , 0, -1 );
      

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
                  inputData->kappa,
		  inputData->chi,
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
                  inputData->bank_chisq,
		  inputData->bank_chisq_dof,
		  inputData->cont_chisq_dof,
		  inputData->cont_chisq,
		  inputData->sigmasq,
		  inputData->rsqveto_duration,
		  inputData->Gamma[0],
		  inputData->Gamma[1],  
		  inputData->Gamma[2],
		  inputData->Gamma[3],
		  inputData->Gamma[4],
		  inputData->Gamma[5],
		  inputData->Gamma[6],
		  inputData->Gamma[7],
		  inputData->Gamma[8],
		  inputData->Gamma[9],
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
  while   ((fgets(sbuf, 2048, input1))!= NULL)
  {
    sscanf(sbuf,BANKEFFICIENCY_PARAMS_ROW_SPACE,
        &trigger.psi0_trigger, &trigger.psi3_trigger, 
	&psi0, &psi3,  &tau0, &tau3, &tau0I, &tau3I, &ecc,&eccI,
	&trigger.fend_trigger, &trigger.fend_inject,
	&trigger.mass1_inject, &trigger.mass2_inject,
	&phaseI, &trigger.rho_final, &trigger.snrAtCoaTime, &trigger.phase,
	&trigger.alphaF, &trigger.bin, &nStartPad, &trigger.nfast, &nfast_max); 

      
    fprintf(output, BANKEFFICIENCY_PARAMS_ROW,
        trigger.psi0_trigger, trigger.psi3_trigger,
	psi0, psi3, tau0, tau3, tau0I, tau3I,ecc,eccI,
        trigger.fend_trigger, trigger.fend_inject,
	trigger.mass1_inject, trigger.mass2_inject,
	phaseI, trigger.rho_final, trigger.snrAtCoaTime, trigger.phase, 
        trigger.alphaF, trigger.bin, nStartPad, trigger.nfast, nfast_max); 
    fprintf(output,",\n");
      
    countline++;
  }


  fprintf(stderr,"read %d lines...done\n", countline);
  PRINT_LIGOLW_XML_TABLE_FOOTER(output);
  PRINT_LIGOLW_XML_FOOTER(output);

  fclose(output);
  fprintf(stderr,"closing xml file\n");
}



/* 
 * These two next functions are hack version of InspiralBankGeneration and
 * FineBank so that we can call a hybrid version of fineBank as if it was
 * available in BankGenerarion. Therefore, hardly any change is needed both in
 * lal and BankEfficiency. Once the code with this fine option is fully
 * tested, one can incorporate the change within lal. 
 * */

void
LALInspiralBankGeneration2(
     LALStatus *status,
     InspiralCoarseBankIn *input,
     SnglInspiralTable **first,
     INT4 *ntiles, 
     UserParametersIn *userParam)
{
  InspiralTemplateList *coarseList = NULL;
  InspiralTemplateList *fineList = NULL;
  InspiralFineBankIn *fineIn = NULL;
  SnglInspiralTable *bank;
  INT4 cnt = 0, flist =0;
  
  INITSTATUS(status, "LALInspiralBankGeneration2", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);
    
  ASSERT( input != NULL, status, LALINSPIRALBANKH_ENULL,
          LALINSPIRALBANKH_MSGENULL );
  ASSERT( *first == NULL, status, LALINSPIRALBANKH_ENULL,
          LALINSPIRALBANKH_MSGENULL );

  /* For nonspinning approximants, call LALInspiralCreateCoarseBank(). */
  switch( input->approximant )
  {
  case BCV:
  case EOB:
  case PadeT1:
  case PadeF1:
  case TaylorF1:
  case TaylorF2:
  case TaylorT1:
  case TaylorT2:
  case TaylorT3:
  case AmpCorPPN:
  case Eccentricity:

    /* Use LALInspiralCreateCoarseBank(). */
    TRY( LALInspiralCreateCoarseBank( status->statusPtr, &coarseList, ntiles,
         *input ), status );
    
    fineIn = (InspiralFineBankIn *)LALMalloc(sizeof(InspiralFineBankIn));
    fineIn->coarseIn = *input;
    fineIn->templateList = coarseList[0];
    fineIn->templateList.params.t0 = 0.9384;
    fineIn->templateList.params.t3 = 0.235;
    
    TRY( LALInspiralCreateFineBank2(status->statusPtr, &fineList, &flist, *fineIn, userParam), status);
    if (fineIn!=NULL) 
      LALFree(fineIn);           

    /* Convert output data structure. */
    bank = (SnglInspiralTable *) LALCalloc(1, sizeof(SnglInspiralTable));
    if (bank == NULL){
      ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
    }
    *first = bank;
    for( cnt = 0; cnt < flist; cnt++ )
    {

      /* do we want a template bank in the eccentyricity dimension (uniform
       * distribution ? */
      bank = bank->next = (SnglInspiralTable *) LALCalloc( 1, sizeof(
           SnglInspiralTable ) );
      if (bank == NULL)
      {
        ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
      }
      bank->mass1 = fineList[cnt].params.mass1;
      bank->mass2 = fineList[cnt].params.mass2;
      bank->mchirp = fineList[cnt].params.chirpMass;
      bank->mtotal = fineList[cnt].params.totalMass;
      bank->eta = fineList[cnt].params.eta;
      bank->tau0 = fineList[cnt].params.t0;
      bank->tau2 = fineList[cnt].params.t2;
      bank->tau3 = fineList[cnt].params.t3;
      bank->tau4 = fineList[cnt].params.t4;
      bank->tau5 = fineList[cnt].params.t5;
      bank->ttotal = fineList[cnt].params.tC;
      bank->psi0 = fineList[cnt].params.psi0;
      bank->psi3 = fineList[cnt].params.psi3;
      bank->f_final = fineList[cnt].params.fFinal;
      bank->eta = fineList[cnt].params.eta;
      bank->beta = fineList[cnt].params.beta;
      /* Copy the 10 metric co-efficients ... */
      memcpy (bank->Gamma, fineList[cnt].metric.Gamma, 10*sizeof(REAL4));
          
    }
    /* Free first template, which is blank. */
    bank = (*first)->next;
    LALFree( *first );
    *first = bank;
    /* free the coarse list returned by create coarse bank */
    LALFree( fineList );
    LALFree( coarseList );
    break;

  default:
    ABORT( status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE );

  }

  *ntiles = flist; 
  DETATCHSTATUSPTR(status);
  RETURN(status); 
}



void LALInspiralCreateFineBank2(LALStatus  *status,
    InspiralTemplateList **outlist,
    INT4                 *nlist,
    InspiralFineBankIn   fineIn,
    UserParametersIn *userParam)
{ /*
     </lalVerbatim> */

  
  REAL8  x0FineMin, x1FineMin, x0FineMax, x1FineMax;
  INT4  i, j, validPars, bins0, bins1;
  InspiralTemplate   *tempPars=NULL;
  InspiralBankParams *bankPars=NULL;
  
  
  INITSTATUS (status, "LALInspiralCreateFineBank", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);
  ASSERT ((INT4)fineIn.coarseIn.space>=0,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
  ASSERT ((INT4)fineIn.coarseIn.space<=1,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
  ASSERT ((REAL8)fineIn.templateList.params.t0 > 0, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
  
  *nlist = 0;
  
  tempPars = (InspiralTemplate *) LALCalloc(1, sizeof(InspiralTemplate));
  bankPars = (InspiralBankParams *) LALCalloc(1, sizeof(InspiralBankParams));
  *tempPars = fineIn.templateList.params;
  switch (fineIn.coarseIn.space)
  {
    case Tau0Tau2:
      ASSERT (fineIn.templateList.params.t2 > 0, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
      bankPars->x0 = fineIn.templateList.params.t0;
      bankPars->x1 = fineIn.templateList.params.t2;
      break;
    case Tau0Tau3:
      ASSERT (fineIn.templateList.params.t3 > 0, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
      bankPars->x0 = fineIn.templateList.params.t0;
      bankPars->x1 = fineIn.templateList.params.t3;
      break;
    default: /* JC: DEFAULT CASE ADDED HERE */
      ABORT( status, 9999, "Default case in switch." );
  }
  
  x0FineMin = userParam->t0FineMin;
  x0FineMax = userParam->t0FineMax;
  x1FineMin = userParam->t3FineMin;
  x1FineMax = userParam->t3FineMax;
  bins0 = userParam->t0FineBin;
  bins1 = userParam->t3FineBin;
  bankPars->dx0 = (x0FineMax - x0FineMin)/((float)bins0 -1);
  bankPars->dx1 = (x1FineMax - x1FineMin)/((float)bins1 -1);
  
  
  
  bankPars->x1 = x1FineMin;
  for(i=0; i<=bins1; i++) {
    bankPars->x0 = x0FineMin;
    for(j=0; j<=bins0; j++) {
      LALInspiralValidTemplate(status->statusPtr, &validPars, *bankPars, fineIn.coarseIn);
      CHECKSTATUSPTR(status);
      if (validPars) {
        LALInspiralComputeParams(status->statusPtr, tempPars, *bankPars, fineIn.coarseIn);
        CHECKSTATUSPTR(status);
        
        if (!(*outlist = (InspiralTemplateList*)
              LALRealloc(*outlist, sizeof(InspiralTemplateList)*(*nlist+1)))) {
          ABORT(status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
          outlist = NULL;
        }
        memset( *outlist + *nlist, 0, sizeof(InspiralTemplateList) );
        (*outlist)[*nlist].params = *tempPars;
        ++(*nlist);
      }
      bankPars->x0+=bankPars->dx0;
    }
    bankPars->x1+=bankPars->dx1;
  }
  
  if (tempPars!=NULL) LALFree(tempPars);
  if (bankPars!=NULL) LALFree(bankPars);
  
  DETATCHSTATUSPTR(status);
}
    


