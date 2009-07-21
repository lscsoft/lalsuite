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
/* --- version information ------------------------------------------------ */
NRCSID( BANKEFFICIENCYC, "$Id$");
RCSID(  "$Id$");

#define CVS_ID_STRING_C      "$Id$"
#define CVS_REVISION_C       "$Revision$"
#define CVS_NAME_STRING_C    "$Name$"
#define CVS_SOURCE_C         "$Source$"
#define CVS_DATE_C           "$Date$"
#define PROGRAM_NAME         "BankEfficiency"

/* --- global variables --------------------------------------------------- */
RandomParams  *randParams=NULL;
INT4          randnStartPad = 0;  /* injections always at the same time if zero*/
INT4          ascii2xml = 0;

/* --- Main program ------------------------------------------------------- */
int
main (INT4 argc, CHAR **argv )
{
  /* --- Local variables -------------------------------------------------- */
  INT4                   i;
  INT4                   j;
  INT4                   thisTemplateIndex;
  LALPNOrder             tempOrder;  /* temporary phase order */

  /* --- the general input structure --- */
  UserParametersIn       userParam;

  /* --- signal related --- */
  REAL4Vector            signal;
  static RandomInspiralSignalIn randIn;

  /* --- template bank related --- */
  Mybank                 mybank;
  static InspiralTemplate insptmplt;
  InspiralCoarseBankIn   coarseBankIn;
  INT4                   sizeBank = 0;
  MetadataTable          templateBank;
  SnglInspiralTable      *tmpltHead = NULL;

  /* --- filtering related --- */
  REAL4Vector               correlation;
  BankEfficiencyBCV         bankefficiencyBCV;
  InspiralWaveOverlapIn     overlapin;

  /* --- results and data mining --- */
  OverlapOutputIn        overlapOutputThisTemplate;
  OverlapOutputIn        overlapOutputBestTemplate;

  /* --- fft related --- */
  RealFFTPlan           *fwdp = NULL;
  RealFFTPlan           *revp = NULL;

  /* --- others ---*/
  LALStatus              status = blank_status;
  FILE                  *Foutput;

  /* --- fast option related --- */
  UINT4                   fast_index;

  /* --- eccentricity related to the bank */
  REAL4                  tempEccentricity = 0;

  /* --- ambiguity function and statistics --- */
  gsl_histogram         *histogramNoise = gsl_histogram_alloc (200);
  gsl_matrix            *amb1;

  /* --- parameter for the simulation --- */
  BankEfficiencySimulation   simulation;

  /* --- Some initialization --- */
  gsl_histogram_set_ranges_uniform (histogramNoise, 0.,20.);
  lal_errhandler = LAL_ERR_EXIT;
  lalDebugLevel = 0;
  templateBank.snglInspiralTable = NULL;
  memset( &templateBank, 0, sizeof( MetadataTable ) );
  memset( &randIn, 0, sizeof( RandomInspiralSignalIn ) );
  memset( &signal, 0, sizeof( REAL4Vector ) );
  memset( &mybank, 0, sizeof( Mybank ) );
  memset( &insptmplt, 0, sizeof( InspiralTemplate ) );
  memset( &coarseBankIn, 0, sizeof( InspiralCoarseBankIn ) );
  memset( &correlation, 0, sizeof( REAL4Vector ) );
  memset( &bankefficiencyBCV, 0, sizeof( BankEfficiencyBCV ) );
  memset( &overlapin, 0, sizeof( InspiralWaveOverlapIn ) );
  memset( &overlapOutputThisTemplate, 0, sizeof( OverlapOutputIn ) );
  memset( &overlapOutputBestTemplate, 0, sizeof( OverlapOutputIn ) );
  memset( &simulation, 0, sizeof( BankEfficiencySimulation ) );


  /* --- Initialization of structures related to all the user parameters --- */
  BankEfficiencyParametersInitialization(&coarseBankIn, &randIn, &userParam);

  /* --- Read user parameters --- */
  BankEfficiencyParseParameters(&argc, argv,&coarseBankIn, &randIn, &userParam);

  /* --- Check input parameters --- */
  BankEfficiencyUpdateParams(&coarseBankIn, &randIn, &userParam);

  /* eccentric bank initialisation */
  BankEfficiencyEccentricBankInit(&userParam);

  /* --- Init a random structure using the user seed --- */
  randParams = XLALCreateRandomParams(randIn.useed );

  /* --- Ascii2XML option --- */
  /* This is a call to the function that converts the ASCII output of
   * BankEfficiency into a standard XML output. Useful when several ASCII
   * output are available and needs to be converted all together into a single
   * XML file.
   * */
  if (ascii2xml == 1){
    BankEfficiencyAscii2Xml();
  }

  /* --- Condor related --- */
  /* If one want to use condor script to keep track of the input parameters,
   * we first need to create a prototype containing all parameters, which can
   * be read later using --ascii2xml option. Use with bep.dag DAGMan script.
  */
  if (userParam.printPrototype){
    BankEfficiencyPrintProtoXml(coarseBankIn, randIn, userParam);
  }

  /* --- Allocate memory -------------------------------------------------- */
  /* --- First, estimate the size of the signal --- */
  LAL_CALL(BankEfficiencyGetMaximumSize(&status,
      randIn, coarseBankIn, userParam, &(signal.length)),
      &status);

  /* --- Set size of the other vectors --- */
  randIn.psd.length     = signal.length/2 + 1;
  correlation.length    = signal.length;
  /* --- Allocate memory for some vectors --- */
  signal.data      = (REAL4*)LALCalloc(1, sizeof(REAL4) * signal.length);
  correlation.data = (REAL4*)LALCalloc(1, sizeof(REAL4) * correlation.length);
  randIn.psd.data  = (REAL8*)LALCalloc(1, sizeof(REAL8) * randIn.psd.length);
  /* --- Allocate memory for BCV filtering only --- */
  if (userParam.template == BCV)
  {
  	INT4 n;
  	n = signal.length;
    bankefficiencyBCV.FilterBCV1.length = n;
    bankefficiencyBCV.FilterBCV2.length = n;
    bankefficiencyBCV.FilterBCV1.data = (REAL4*)LALCalloc(1, sizeof(REAL4) * n);
    bankefficiencyBCV.FilterBCV2.data = (REAL4*)LALCalloc(1, sizeof(REAL4) * n);
  }

  /* --- Create the PSD noise ----------------------------------------------- */
  LAL_CALL(BankEfficiencyCreatePsd(&status, &coarseBankIn, &randIn, userParam),
       &status);

  /* --- Create the template bank ------------------------------------------- */
  LAL_CALL(BankEfficiencyCreateTemplateBank(&status, &coarseBankIn,
      &templateBank, &tmpltHead, userParam, &randIn, &sizeBank), &status);

  /* --- populate MyBank strucuture --- */
  BankEfficiencyInitMyBank(&mybank, &sizeBank, tmpltHead, userParam);

  /* --- Allocate memory for the ambiguity function ------------------------- */
  /* For each bank, we have t0,t3,and the maximum */
  amb1 = gsl_matrix_alloc(4,sizeBank); /*allocated t0,t3,max*/
  for (i=0;i<4; i++)
  {
    for (j=0;j<sizeBank; j++)
    {
      gsl_matrix_set(amb1, i, j, 0.); /* set it to zero */
    }
  }

  /* --- Estimate the fft's plans ----------------------------------------- */
  LALCreateForwardRealFFTPlan(&status, &fwdp, signal.length, 0);
  LALCreateReverseRealFFTPlan(&status, &revp, signal.length, 0);

  /* --- The overlap structure -------------------------------------------- */
  overlapin.nBegin      = 0;
  overlapin.nEnd        = 0;
  overlapin.psd         = randIn.psd;
  overlapin.fwdp        = randIn.fwdp = fwdp;
  overlapin.revp        = revp;
  overlapin.ifExtOutput = 0;  /* new option from Anand's WaveOverlap */
  overlapin.ifCorrelationOutput = 1; /* output of WaveOverlap is sqrt(x^2+y^2)*/

  /* --- Create Matrix for BCV Maximization once for all ------------------ */
  if ((userParam.template == BCV) )
  {
    LAL_CALL( BankEfficiencyCreatePowerVector(&status,
       &(bankefficiencyBCV.powerVector), randIn, signal.length),  &status);

    BankEfficiencyCreateBCVMomentVector(
        &(bankefficiencyBCV.moments), &coarseBankIn.shf,randIn.param.tSampling,
        randIn.param.fLower, signal.length);
  }

  /* ------------------------------------------------------------------------ */
  /* --- Main loop ---------------------------------------------------------- */
  /* ------------------------------------------------------------------------ */
  /* --- ntrials is used to keep track of the simulation trials           --- */
  simulation.ntrials = 0;
  simulation.N = userParam.ntrials;
  /* --- loop over the number of simulations (ntrials) --- */
  while (++simulation.ntrials <= userParam.ntrials)
  {
    if (vrbflg){
      fprintf(stderr,"Simulation number %d/%d\n", simulation.ntrials, userParam.ntrials);
    }

    /* --- initialisation for each simulation --- */
    BankEfficiencyInitOverlapOutputIn(&overlapOutputBestTemplate);

    /* --- set the fcutoff (signal,overlap) --- */
    randIn.param.fCutoff = userParam.signalfFinal;
    randIn.param.inclination = 0.;
    randIn.param.distance = 1.;

    /* --- generate the signal waveform --- */
    LAL_CALL(BankEfficiencyGenerateInputData(&status,
        &signal, &randIn, userParam), &status);

    /* --- populate the main structure of the overlap ---*/
    overlapin.signal = signal;

    /*  --- populate the insptmplt with the signal parameter --- */
    insptmplt = randIn.param; /* set the sampling and other common parameters */
    insptmplt.approximant  = userParam.template;        /* set the template   */
    insptmplt.order        = coarseBankIn.order;        /* set the phaseorder */
    insptmplt.ampOrder     = coarseBankIn.ampOrder;     /* set the amp order  */

    /* --- reset the bank (none of the template have been used) --- */
    for (i=0; i<sizeBank; i++){
      mybank.used[i] = 0;
    }

    /* --- set the general parameter for each simulation --- */
    simulation.filter_processed = 0; /* count number of templates being used  */
    simulation.lastEMatch       = 1; /* ematch of the last template being used*/
    simulation.filteringIndex   = 0; /* current filter number                 */
    simulation.stop             = 0; /* flag to stop this simulation          */
    simulation.bestSNR          = 0; /* value of the best SNR                 */
    simulation.bestSNRIndex     = 0; /* index of the template that has the best SNR*/
    simulation.SNRMax           = 0; /* For how long have we the best SNR     */
    simulation.fastParam1       = userParam.fastParam1;
    simulation.fastParam1       *= 1.+log10(userParam.eccentricBank.bins);
    simulation.eMatch           = userParam.eMatch; /* reset the ematch       */
    simulation.bestEMatch       = -1e16; /* reset ematch to small value*/
    simulation.randIn           = randIn;

    /* --- loop over the bank until counter reaches size bank --- */
    while ( simulation.filteringIndex<sizeBank && simulation.stop == 0 )
    {
      /* --- iterate the counter --- */
      simulation.filteringIndex++;

      /* -- reset the output structure --- */
      BankEfficiencyInitOverlapOutputIn(&overlapOutputThisTemplate);

      GetClosestValidTemplate(mybank, simulation.randIn, &fast_index);

      LAL_CALL(BankEfficiencyCreateListFromTmplt(&status,
             &insptmplt, mybank, fast_index), &status);
      thisTemplateIndex = fast_index;
      insptmplt.inclination = 0.;

      switch(userParam.template)
      {
        case BCV:

          insptmplt.massChoice = psi0Andpsi3;
          LAL_CALL(LALInspiralParameterCalc(&status, &(insptmplt)),&status);
          overlapin.param = insptmplt;

          LAL_CALL(LALInspiralParameterCalc( &status,  &(overlapin.param) ),
              &status);

          /* if faithfulness is required, the template parameters are
           * identical to the signal except for the name of the approximant.
           *  */
          if (userParam.faithfulness)
          {
            insptmplt                    = randIn.param;
            overlapin.param              = randIn.param;
            overlapin.param.approximant  = userParam.template;
            insptmplt.fFinal             = randIn.param.fFinal;
            /* --- the faithfulness requires only one filtering, so we stop - */
            simulation.stop              = 1;
           }

          LAL_CALL(BankEfficiencyInspiralOverlapBCV(&status,
                    &insptmplt, userParam, &randIn,
                    &overlapin, &overlapOutputThisTemplate,
                    &correlation, &bankefficiencyBCV),
                    &status);

          /* keep track of what has been done in the overlap function */

          overlapOutputThisTemplate.freq =  overlapin.param.fFinal;
          overlapOutputThisTemplate.templateNumber = thisTemplateIndex;
          simulation.filter_processed++;
          break;

        case AmpCorPPN:
        case TaylorT1:
        case Eccentricity:
        case TaylorT2:
        case TaylorT3:
        case TaylorF1:
        case TaylorF2:
        case EOB:
        case EOBNR:
        case PadeT1:
        case PadeF1:
        case SpinTaylor:
        case TaylorEt:
        case TaylorT4:
        case TaylorN:
	  /*
          if (vrbflg){
            fprintf(stderr,"closest template is tau0= %f tau3= %f  index=%d filterintIndex= %d\n",
              mybank.tau0[fast_index], mybank.tau3[fast_index],fast_index,simulation.filteringIndex);
          }
	  */
          /* --- first we create the template --- */
          insptmplt.massChoice = t03;
          LAL_CALL(LALInspiralParameterCalc( &status,  &(insptmplt) ),
              &status);
          /* --- set up the overlapin strucutre --- */
          overlapin.param = insptmplt;
          LAL_CALL(LALInspiralParameterCalc( &status,  &(overlapin.param) ),
              &status);
          overlapin.param.fCutoff = coarseBankIn.fUpper;
          overlapin.param.fLower  = coarseBankIn.fLower;
          overlapin.param.fFinal  = randIn.param.tSampling/2. - 1;

          /* --- the fast simulation option. We compute the ematch --- */
          simulation.lastEMatch = BankEfficiencyComputeEMatch(&randIn, mybank,
              thisTemplateIndex);
          if (simulation.lastEMatch > simulation.bestEMatch)
            simulation.bestEMatch = simulation.lastEMatch;

          /* --- if check option is on --- */
          if (userParam.faithfulness)
          {
            simulation.lastEMatch = 1; /* if faithfulness, masses are equal so match is 1*/
            tempOrder = insptmplt.order;
            tempEccentricity = insptmplt.eccentricity;
            /* --- now, we overwrite insptmplt with the signal params---*/
            insptmplt = randIn.param;
            overlapin.param = randIn.param;
            /* --- but get back some parameters --- */
            LAL_CALL(LALInspiralParameterCalc( &status,  &(overlapin.param) ), &status);
            overlapin.param.fCutoff = coarseBankIn.fUpper;
            overlapin.param.fFinal = randIn.param.tSampling/2. - 1;
            overlapin.param.approximant        = userParam.template;
            overlapin.param.fLower  = coarseBankIn.fLower;
            /* we want to keep the original order except for eccentricity
             * where order must be zero anyway*/
            overlapin.param.order              = tempOrder;
            overlapin.param.eccentricity       = tempEccentricity;
            insptmplt.eccentricity             = tempEccentricity;
            /* --- the faithfulness requires only one filtering, so we stop - */
            simulation.stop = 1;
          }
          /* if we want to cut integration before the fFinal */
          /*              overlapin.param.fCutoff = 1023;    */


          if (userParam.fastSimulation == Heuristic1)
          {
            /* --- decrease the ematch parameter if needed --- */
            if (simulation.lastEMatch < simulation.eMatch){
              simulation.eMatch = simulation.lastEMatch-1;
            }
            /* --- if SNRMax has reached the max iteration, we stop --- */
            if (simulation.SNRMax >= simulation.fastParam1
               && userParam.fastSimulation){
              simulation.stop = 1;
            }
          }

          /* --- filteing!=1 ensure that at least 1 filter will be used    ---*/
          /* --- while lastEMatch is greater than the ematch --- */
          if ((simulation.filteringIndex !=1) &&
              (userParam.fastSimulation != None)
              && (simulation.lastEMatch < simulation.eMatch ))
          {
            gsl_matrix_set(amb1,2,thisTemplateIndex, 0);
            gsl_matrix_set(amb1,3,thisTemplateIndex, 0);
            gsl_matrix_set(amb1,0,thisTemplateIndex, insptmplt.t0);
            gsl_matrix_set(amb1,1,thisTemplateIndex, insptmplt.t3);
            /* --- if we enter here, then we should stop the simulation --- */
            simulation.stop = 1;
          }
          else
          {
            LAL_CALL(BankEfficiencyWaveOverlap(&status,
                  &correlation,
                  &overlapin,
                  &overlapOutputThisTemplate,
                  randIn.param.nStartPad), &status);


            overlapOutputThisTemplate.templateNumber = thisTemplateIndex;
/*
	    fprintf(stderr, "%d\n", thisTemplateIndex);
*/
            /* we compute the averaged ambiguity function a t=ta and
             * the averaged maximizaed ambiguity function over time*/
            BankEfficiencyPopulateAmbiguityFunction(
                amb1,correlation,thisTemplateIndex,
                overlapOutputThisTemplate, randIn.param.nStartPad,
                insptmplt);

            /* is it correct to do so ? */
            insptmplt.fFinal = overlapin.param.fFinal;
            simulation.filter_processed++;

            if (vrbflg){
            fprintf(stderr, "snr=%f m1T=%f m2T=%f m1S=%f m2S=%f SigPad=%d Maxbin=%d\n",
                overlapOutputThisTemplate.rhoMax, insptmplt.mass1,
		insptmplt.mass2, randIn.param.mass1, randIn.param.mass2,
                randIn.param.nStartPad, overlapOutputThisTemplate.rhoBin);

            fflush(stderr);
            }

            if (overlapOutputThisTemplate.rhoMax > simulation.bestSNR)
            {
              simulation.bestSNR  = overlapOutputThisTemplate.rhoMax;
              simulation.bestSNRIndex = simulation.filteringIndex;
              simulation.SNRMax = 0;
              if (userParam.fastSimulation == Heuristic1)
              {
                simulation.randIn.param.t0 = insptmplt.t0;
                simulation.randIn.param.t3 = insptmplt.t3;
              }
              simulation.bestEMatch      = simulation.lastEMatch;
            }
            else
            {
              simulation.SNRMax++;
            }

          }
        break;
      } /*end of the switch over template*/



      /* accumulates histogram of the correlations over all templates and all
       * simulations. */
      if (userParam.printSNRHisto){
        BankEfficiencyUpdateSNRHistogram(&correlation, histogramNoise);
      }

      /* for the final results, we keep only the loudest SNR over the entire
       * template bank, for each simulations.*/
      BankEfficiencyKeepHighestValues( overlapOutputThisTemplate,
          &overlapOutputBestTemplate, insptmplt);

    }  /* --- end of  bank process --- */

    /* --- print the results on stdout and xml file if requested --- */
    LAL_CALL(BankEfficiencyFinalise(&status,mybank,
        overlapOutputBestTemplate,randIn,userParam,simulation,
        coarseBankIn), &status);

    if (userParam.template == BCV)
    {
      if (userParam.printBestOverlap || userParam.printBestTemplate)
      {
        userParam.extraFinalPrinting = 1;
        LAL_CALL(BankEfficiencyInspiralOverlapBCV(&status,
                    &insptmplt, userParam, &randIn,
                     &overlapin, &overlapOutputThisTemplate, &correlation,
                     &bankefficiencyBCV),
                     &status);
      }
    }

    /* --- print correlation of the last computed correlation --- */
    if (userParam.printBestOverlap ){
    	BankEfficiencySaveVector("BankEfficiency-correlation.dat", correlation,
    	    randIn.param.tSampling);
    }

  }  /* --- end while(trial) --- */


  /* --- save histrogram of all the SNR output correlations --- */
  if (userParam.printSNRHisto)
  {
    Foutput = fopen(BANKEFFICIENCY_SNRHISTOGRAM,"w");
    gsl_histogram_fprintf(Foutput, histogramNoise, "%f", "%g");
    fclose(Foutput);
  }


  /* --- we save the ambiguity in a file --- */
  BankEfficiencyPrintAmbiguity(userParam,sizeBank,amb1);

  /* free memory */
  while ( templateBank.snglInspiralTable )
  {
    tmpltHead = templateBank.snglInspiralTable;
    templateBank.snglInspiralTable = templateBank.snglInspiralTable->next;
    LALFree( tmpltHead );
  }

  /* --- destroy the plans, correlation and signal --- */

  XLALDestroyRandomParams(randParams);

  if (userParam.template == BCV)
  {
    LALFree(bankefficiencyBCV.powerVector.fm5_3.data);
    LALFree(bankefficiencyBCV.powerVector.fm2_3.data);
    LALFree(bankefficiencyBCV.powerVector.fm1_2.data);
    LALFree(bankefficiencyBCV.powerVector.fm7_6.data);
    LALFree(bankefficiencyBCV.moments.a11.data);
    LALFree(bankefficiencyBCV.moments.a21.data);
    LALFree(bankefficiencyBCV.moments.a22.data);

    if (userParam.template == BCV)
    {
      LALFree(bankefficiencyBCV.FilterBCV2.data);
      LALFree(bankefficiencyBCV.FilterBCV1.data);
    }
  }

  LALFree(randIn.psd.data);
  LALDDestroyVector( &status, &(coarseBankIn.shf.data) );

  LALFree(signal.data);
  LALFree(correlation.data);

  LALDestroyRealFFTPlan(&status,&fwdp);
  LALDestroyRealFFTPlan(&status,&revp);

  gsl_matrix_free(amb1);
  gsl_histogram_free(histogramNoise);

  free(mybank.mass1);
  free(mybank.mass2);
  /*todo  free other varaible in mybank*/

  LALCheckMemoryLeaks();

  return 0;
}




/* *************************************************************************
 * Some output Results
 * ********************************************************************** */
void BankEfficiencyKeepHighestValues(
  OverlapOutputIn  resultThisTemplate,
  OverlapOutputIn *resultBestTemplate,
  InspiralTemplate insptmplt)
{

  if (resultThisTemplate.rhoMax > resultBestTemplate->rhoMax)
  {
    resultBestTemplate->rhoMax         = resultThisTemplate.rhoMax;
    resultBestTemplate->phase          = resultThisTemplate.phase;
    resultBestTemplate->alpha          = resultThisTemplate.alpha;
    resultBestTemplate->rhoBin         = resultThisTemplate.rhoBin;
    resultBestTemplate->freq           = resultThisTemplate.freq;
    resultBestTemplate->templateNumber = resultThisTemplate.templateNumber;
    resultBestTemplate->snrAtCoaTime   = resultThisTemplate.snrAtCoaTime;
    resultBestTemplate->eccentricity   = insptmplt.eccentricity;
  }
}


/* ************************************************************************
 * Some output Results
 * ********************************************************************** */
void BankEfficiencyGetResult(
  LALStatus        *status,
  InspiralTemplate *list,
  InspiralTemplate  injected,
  OverlapOutputIn   bestOverlap,
  ResultIn         *result,
  UserParametersIn  userParam
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
    result->psi0_trigger = trigger.psi0;
    result->psi3_trigger = trigger.psi3;
    result->tau0_trigger = 0.;
    result->tau3_trigger = 0.;
    result->tau0_inject  = 0.;
    result->tau3_inject  = 0.;
  }
  else
  {
    LALInspiralParameterCalc( status->statusPtr,  &trigger );
    CHECKSTATUSPTR(status);
    result->psi0_inject  = 0.;
    result->psi3_inject  = 0.;
    result->psi0_trigger = 0.;
    result->psi3_trigger = 0.;
    result->tau0_trigger = trigger.t0;
    result->tau3_trigger = trigger.t3;
    result->tau0_inject  = injected.t0;
    result->tau3_inject  = injected.t3;
    result->polarisationAngle  = injected.polarisationAngle;
    result->inclination  = injected.inclination;
  }

  result->mass1_trigger = trigger.mass1;
  result->mass2_trigger = trigger.mass2;
  result->mass1_inject = injected.mass1;
  result->mass2_inject = injected.mass2;
  result->fend_inject  = injected.fFinal;
  result->fend_trigger = bestOverlap.freq;
  result->rho_final    = bestOverlap.rhoMax;
  result->alphaF       = bestOverlap.alpha*pow(bestOverlap.freq,2./3.);
  result->bin          = bestOverlap.rhoBin;
  result->phase        = bestOverlap.phase;
  result->layer        = bestOverlap.layer;
  result->phase        = bestOverlap.phase;
  result->snrAtCoaTime = bestOverlap.snrAtCoaTime;
  result->eccentricity = bestOverlap.eccentricity;

  DETATCHSTATUSPTR(status);
  RETURN (status);
}


void BankEfficiencyPrintResults(
  ResultIn                 result,
  RandomInspiralSignalIn   randIn,
  BankEfficiencySimulation simulation)
{
  FILE *fs;
  fprintf(stdout,
  "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %d\n",
      result.mass1_trigger, result.mass2_trigger,
/*
      randIn.param.psi0, randIn.param.psi3,
      result.tau0_trigger, result.tau3_trigger,
      randIn.param.t0, randIn.param.t3,
      result.eccentricity, randIn.param.eccentricity,randIn.param.alpha1,
*/
      randIn.param.mass1,randIn.param.mass2,
      result.fend_trigger, randIn.param.fFinal,
/*
      randIn.param.inclination,randIn.param.polarisationAngle,
      randIn.param.startPhase,
      result.snrAtCoaTime, result.phase,
      result.alphaF,result.bin, randIn.param.nStartPad, result.nfast,
*/
      result.rho_final,
      result.bin);
  fflush(stdout);
}


/* ****************************************************************************
 * Function to generate and stored the moments in  REAL4vectors.
 * ************************************************************************* */
void BankEfficiencyCreateBCVMomentVector(
  BankEfficiencyMoments *moments,
  REAL8FrequencySeries  *psd,
  REAL8                 sampling,
  REAL8                 fLower,
  INT4                  length)
{
  REAL8 m7 = 0.;                            /* the moments */
  REAL8 m5 = 0.;
  REAL8 m3 = 0.;
  REAL8 f;

  UINT4 kMin;
  UINT4 kMax;
  UINT4 k;

  InspiralMomentsIn     in;

  moments->a11.length =  length / 2 ;
  moments->a22.length =  length / 2 ;
  moments->a21.length =  length / 2 ;
  moments->a11.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) * moments->a11.length);
  moments->a22.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) * moments->a22.length);
  moments->a21.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) * moments->a21.length);

  /* inspiral structure */
  in.shf    = psd;          /* The spectrum             */
  in.xmin   = fLower;       /* Lower frequency of integral? is it correct or
                               should we put zero? */
  in.xmax   = sampling/2.;  /* We dont' need to cut here */
  in.norm   = 1./4.;        /* Some normalization */

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
  BankEfficiencySaveVector("moment1.dat",moments->a11,sampling);
}


/* ****************************************************************************
 * Create an orthogonal filter. by taking its complex conjuguate
 * ************************************************************************** */
void BankEfficiencyGetOrthogonalFilter(
  REAL4Vector *filter)
{
  UINT4   i;
  UINT4   n      = filter->length;
  UINT4   nby2   = filter->length / 2;
  REAL4   temp;

  for (i = 1; i < nby2; i++)
  {
    temp              =  filter->data[i];
    filter->data[i]   = -filter->data[n-i];
    filter->data[n-i] = temp;
  }
}

/* *****************************************************************************
 * Create a vector f^{a/b}
 * ************************************************************************** */
void BankEfficiencyCreateVectorFreqPower(
  REAL4Vector      *vector,
  InspiralTemplate  params,
  INT4              a,
  INT4              b)
{
  INT4  i;
  INT4  n = vector->length; /* Length of the vector to create */

  REAL8 power = (REAL8)a / (REAL8)b;       /* the frequency power */
  REAL8 f;                                 /* frequency */
  REAL8 df = params.tSampling/(REAL8)n/2.; /* sampling frequency */

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
void BankEfficiencyCreateBCVFilters(
  BankEfficiencyBCV *bankefficiencyBCV,
  UINT4              kMin,           /* position of the integration   */
  UINT4              kMax,
  REAL4              psi0,
  REAL4              psi3)
{
  UINT4 i;
  UINT4 n = bankefficiencyBCV->FilterBCV1.length;
  UINT4 nby2 = bankefficiencyBCV->FilterBCV1.length / 2;

  REAL8 amplitude, phase;
  REAL8 cos_phase;
  REAL8 sin_phase;

  REAL4 a11,a21,a22;

  a11 = bankefficiencyBCV->moments.a11.data[kMax];
  a22 = bankefficiencyBCV->moments.a22.data[kMax];
  a21 = bankefficiencyBCV->moments.a21.data[kMax];

  /* Create the templates */
  for (i = 0; i< kMin-1; i++)
  {
    bankefficiencyBCV->FilterBCV1.data[i] = 0;
    bankefficiencyBCV->FilterBCV2.data[i] = 0;
  }
  /* Should add the upper fendBCV here as well */

  for (i = kMin; i< nby2; i++)
  {
    phase  =   psi0 * bankefficiencyBCV->powerVector.fm5_3.data[i]    /* The same for the two filters     */
        + psi3 * bankefficiencyBCV->powerVector.fm2_3.data[i];

    cos_phase  = cos(phase);    /* The same for the two filters     */
    sin_phase  = sin(phase);

    /* Fill the first filter here */
    amplitude  = a11 * bankefficiencyBCV->powerVector.fm7_6.data[i];
    bankefficiencyBCV->FilterBCV1.data[i]   =  amplitude * cos_phase;
    bankefficiencyBCV->FilterBCV1.data[n-i] = -amplitude * sin_phase;

    /* Fill the second one here */
    amplitude =  a21 * bankefficiencyBCV->powerVector.fm7_6.data[i] +
        a22 * bankefficiencyBCV->powerVector.fm1_2.data[i];
    bankefficiencyBCV->FilterBCV2.data[i]   =  amplitude * cos_phase;
    bankefficiencyBCV->FilterBCV2.data[n-i] = -amplitude * sin_phase;
  }
}


/*  **************************************************************************
 *  The Overlap function for BCV templates
 *  **************************************************************************/
void BankEfficiencyWaveOverlapBCV(
  LALStatus               *status,
  REAL4Vector             *correlation,
  InspiralWaveOverlapIn   *overlapin,
  REAL4Vector             *FilterBCV1,
  REAL4Vector             *FilterBCV2,
  UserParametersIn         userParam ,
  OverlapOutputIn         *overlapOutput,
  BankEfficiencyMoments   *moments)
     /*  </lalVerbatim>  */
{
  INT4 rhoBinUnconstraint = 0.;
  INT4 rhoBinConstraint = 0.;

  REAL4 sampling;
  REAL4 rhoMaxConstraint = 0.;
  REAL4 alphaConstraint = 0.;
  REAL4 phaseConstraint = 0.;
  REAL4 rhoMaxUnconstraint = 0.;
  REAL4 alphaUnconstraint = 0.;
  REAL4 phaseUnconstraint = 0.;
  REAL4 rhoConstraint = 0.;
  REAL4 rhoUnconstraint = 0.;
  REAL4 df = 0.;
  REAL4 BestPhase = 0.;
  REAL4 thetab = 0.;
  REAL4 alphaMax = 0.;
  REAL4 a11, a22, a21;
  REAL4 rho = 0;
  REAL4 alphaFU, fm23, fp23;

  REAL4Vector templatee;
  REAL4Vector x1;
  REAL4Vector x2;
  REAL4Vector x3;
  REAL4Vector x4;
  REAL4Vector phaseV;
  REAL4Vector thetavVector;
  REAL4Vector rho1;
  REAL4Vector rho2;
  REAL4Vector rho3;
  REAL4Vector v0;
  REAL4Vector v1;
  REAL4Vector v2;                   /* output of the correlations       */
  REAL4Vector alphaf;               /* output of the correlations       */

  InspiralWaveCorrelateIn     corrin;/* the correlation input structure     */

  UINT4 k;
  UINT4 i;
  UINT4 nBegin;     /* beginning of valid correlation   */
  UINT4 nEnd;       /* end of valid correlation         */
  UINT4 n = correlation->length; /* length of the vectors       */

  REAL4 phi, phase, cphase, sphase, cphi, sphi;

  INITSTATUS (status, "BankEfficiencyWaveOverlapBCV", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);


  /* --- size affectation --- */
  x1.length = x2.length = x3.length = x4.length = n;
  v1.length = v2.length = v0.length = alphaf.length = n;
  phaseV.length = thetavVector.length = n;
  rho1.length = rho2.length = rho3.length = templatee.length = n;

  /* --- memory allocation --- */
  x1.data = (REAL4*) LALMalloc(sizeof(REAL4) * x1.length);
  x2.data = (REAL4*) LALMalloc(sizeof(REAL4) * x2.length);
  x3.data = (REAL4*) LALMalloc(sizeof(REAL4) * x3.length);
  x4.data = (REAL4*) LALMalloc(sizeof(REAL4) * x4.length);

  /* --- more memory allocation --- */
  phaseV.data    = (REAL4*) LALMalloc(sizeof(REAL4) * phaseV.length);
  templatee.data = (REAL4*) LALMalloc(sizeof(REAL4) * templatee.length);
  rho1.data      = (REAL4*) LALMalloc(sizeof(REAL4) * rho1.length);
  rho2.data      = (REAL4*) LALMalloc(sizeof(REAL4) * rho2.length);
  rho3.data      = (REAL4*) LALMalloc(sizeof(REAL4) * rho3.length);
  v1.data        = (REAL4*) LALMalloc(sizeof(REAL4) * v1.length);
  v2.data        = (REAL4*) LALMalloc(sizeof(REAL4) * v2.length);
  v0.data        = (REAL4*) LALMalloc(sizeof(REAL4) * v0.length);
  alphaf.data    = (REAL4*) LALMalloc(sizeof(REAL4) * v0.length);
  thetavVector.data = (REAL4*) LALMalloc(sizeof(REAL4) * v0.length);

  /* --- --- */
  sampling = overlapin->param.tSampling;

  overlapin->param.nStartPad  = 0;
  overlapin->param.nEndPad    = 1;
  overlapin->param.startPhase = 0;

  /*  We want to compute <x|h1> and <x|h2>; let's prepare the correlation
   *  process; let's get the input data namely x and the spectrum namely psd
   *  */

  BankEfficiencySaveVector("signal.dat", overlapin->signal, sampling);
  BankEfficiencySaveVector("filter1.dat", *FilterBCV1, sampling);

  corrin.fCutoff      = overlapin->param.fFinal;    /* the integral's cutoff frequency */
  corrin.samplingRate = overlapin->param.tSampling; /* the sampling          */
  corrin.df           = overlapin->param.tSampling /
      (REAL8) overlapin->signal.length;               /* the resolution in frequency */
  corrin.psd          = overlapin->psd;             /* the spectrum           */
  corrin.signal1      = overlapin->signal;          /* the signal x           */
  corrin.revp         = overlapin->revp;            /* fourier transform flag */

  /* --- Filter h1 and its orthonormal h1 --- */
  corrin.signal2        = *FilterBCV1;                      /* The first template    */
  LALInspiralWaveCorrelate(status->statusPtr, &x1, corrin); /* the correlation       */
  CHECKSTATUSPTR(status);
  BankEfficiencyGetOrthogonalFilter(FilterBCV1); /* get the orthonormalized template */
  corrin.signal2        = *FilterBCV1;                      /* the second template   */
  LALInspiralWaveCorrelate(status->statusPtr, &x3, corrin); /* the correlation       */
  CHECKSTATUSPTR(status);


  /* --- Filter h2 and its orthornormal filter --- */
  corrin.signal2        = *FilterBCV2;                      /* The first template */
  LALInspiralWaveCorrelate(status->statusPtr, &x2, corrin); /* The correlation    */
  CHECKSTATUSPTR(status);
  BankEfficiencyGetOrthogonalFilter( FilterBCV2); /* get the orthonormailzed templates*/
  corrin.signal2        = *FilterBCV2;                      /* the second template          */
  LALInspiralWaveCorrelate(status->statusPtr, &x4, corrin); /* the correlation */

  /* nbegin and end for the correlation processus
   * as well as overlapout affectation */
  nBegin = overlapin->nBegin;                  /* starting bin     */
  nEnd = FilterBCV1->length - overlapin->nEnd; /* ending valid bin */


  /* --- some output for debugging, checking --- */
  if (userParam.printBestOverlap && userParam.extraFinalPrinting)
  {
    BankEfficiencySaveVector("BankEfficiency_Filter.dat", x1, sampling);
    BankEfficiencySaveVectorAppend("BankEfficiency_Filter.dat", x2, sampling);
    BankEfficiencySaveVectorAppend("BankEfficiency_Filter.dat", x3, sampling);
    BankEfficiencySaveVectorAppend("BankEfficiency_Filter.dat", x4, sampling);
  }

  /* --- the BCV filtering --- */
  df = overlapin->param.tSampling / (REAL8)(n/2.) / 2.;
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
    REAL4 thetav;
    REAL4 x1_2, x2_2, x3_2, x4_2;
    REAL4 Vector0, Vector1, Vector2;

    x1_2 = x1.data[i] * x1.data[i];
    x2_2 = x2.data[i] * x2.data[i];
    x3_2 = x3.data[i] * x3.data[i];
    x4_2 = x4.data[i] * x4.data[i];

    Vector0 = x1_2 + x2_2 + x3_2 + x4_2;
    Vector1 = x1_2 + x3_2 - x2_2 - x4_2;
    Vector2 = 2 * (x1.data[i]*x2.data[i] + x3.data[i]*x4.data[i]);

    rhoUnconstraint = sqrt((Vector0 + sqrt(Vector1*Vector1+Vector2*Vector2))/2.);

    /* get thetav first in order to constraint alpha_F*/
    thetav = atan2(Vector2, Vector1);

    if(userParam.printBestOverlap && userParam.extraFinalPrinting)
    {
      rho1.data[i] = rhoUnconstraint;
      rho2.data[i] = sqrt((Vector0 + Vector1)/2.);
      rho3.data[i] = sqrt((Vector0+Vector1*cos(2*thetab)+Vector2*sin(2*thetab))/2.);
      v0.data[i]   = Vector0;
      v1.data[i]   = Vector1;
      v2.data[i]   = Vector2;
      phaseV.data[i] = 0.5 * thetav; /*used to compute alpha*/
      thetavVector.data[i] = thetav;
    }

    if (thetab >= 0)
    {
      if ( 0  <= thetav && thetav <= 2 * thetab)
      {
        rhoConstraint = rhoUnconstraint;
      }
      else if (thetab-LAL_PI <= thetav && thetav < 0)
      {
        rhoConstraint = sqrt((Vector0 + Vector1)/2.);
      }
      else if( (2*thetab  < thetav && thetav<=LAL_PI+1e-4 )
         || (-LAL_PI-1e-4<=thetav && thetav < -LAL_PI+thetab))
      {
        rhoConstraint =sqrt((Vector0+Vector1*cos(2*thetab)+Vector2*sin(2*thetab))/2.);;
      }
      else
      {
      fprintf(stderr,"must not enter here  thetav = %e thetab=%e\n ",
              thetav , thetab);
      exit( 1 );
      }
    }
    else
    {
      if ( 2*thetab  <= thetav && thetav <= 0)
      {
        rhoConstraint = rhoUnconstraint;
      }
      else if (0 < thetav &&  thetav  <= LAL_PI + thetab )
      {
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
        exit( 1 );
      }
    }

    if (userParam.alphaFConstraint == ALPHAFConstraint)
    {
      if(userParam.printBestOverlap || userParam.printSNRHisto)
      {
        correlation->data[i] = rhoConstraint;
      }
    }
    else
    {
      if(userParam.printBestOverlap || userParam.printSNRHisto)
      {
        alphaFU =   -(a22 * tan(.5*thetav))
                / (a11 + a21* tan(.5*thetav)) *fp23;

        alphaf.data[i] = alphaFU;

        if (alphaFU <=1 )
        {
          correlation->data[i] = rhoUnconstraint;
        }
        else
        {
          correlation->data[i] = -1;
        }
      }


      if ( rhoConstraint > rhoMaxConstraint)
      /* Keep the position of the max only  */
      {
        rhoMaxConstraint = rhoConstraint;
        rhoBinConstraint = i;
        phaseConstraint  = 0.5 * thetav;
        alphaConstraint  = -(a22 * tan(phaseConstraint))
          / (a11 + a21* tan(phaseConstraint));
      }

      alphaFU =   -(a22 * tan(0.5*thetav))
        / (a11 + a21 * tan(0.5*thetav)) * fp23;

      alphaf.data[i] = alphaFU;

      if ( rhoUnconstraint > rhoMaxUnconstraint && alphaFU<=1)
      /* Keep the position of the max only  */
      {
        rhoMaxUnconstraint  = rhoUnconstraint;
        rhoBinUnconstraint  = i;
        phaseUnconstraint   = 0.5*thetav;
        alphaUnconstraint = -(a22 * tan(phaseUnconstraint))
        / (a11 + a21 * tan(phaseUnconstraint));
      }
    }
  }

  /* --- Finally get the alpha value corresponding to the best rho --- */
  overlapOutput->rhoMax = rhoMaxConstraint;
  overlapOutput->rhoBin = rhoBinConstraint;
  overlapOutput->alpha  = alphaConstraint;
  overlapOutput->phase  = phaseConstraint;

  /* --- The final template to be print  ---------------------------------- */
  if ( userParam.printBestOverlap && userParam.extraFinalPrinting )
  {
    FILE *Foutput;
    /* print overlap, combinaison of the x_i filters */
    BankEfficiencySaveVector("BankEfficiency_Phase.dat", thetavVector, sampling);

    BankEfficiencySaveVector("BankEfficiency_rho1.dat", rho1, sampling);
    BankEfficiencySaveVector("BankEfficiency_rho2.dat", rho2, sampling);
    BankEfficiencySaveVector("BankEfficiency_rho3.dat", rho3, sampling);

    BankEfficiencySaveVector("BankEfficiency_v0.dat", v1, sampling);
    BankEfficiencySaveVector("BankEfficiency_v1.dat", v1, sampling);
    BankEfficiencySaveVector("BankEfficiency_v2.dat", v1, sampling);

    BankEfficiencySaveVector("BankEfficiency_alpha.dat",alphaf,sampling);

    Foutput=fopen("BankEfficiency_Overlap.dat","w");
    for (i=0; i< correlation->length; i++)
    {
      fprintf(Foutput, "%e\n", fabs(correlation->data[i]));
    }
    fclose(Foutput);
  }


  if (userParam.extraFinalPrinting && userParam.printBestOverlap)
  {
    /* print the final template given phase, alpha ans so on*/

    phase       = BestPhase;
    phi         = overlapin->param.startPhase; /*how to get it ohterwise ? */

    for (i=0; i<(UINT4)correlation->length; i++)
    {
      cphase = cos(phase);
      cphi   = cos(phi);
      sphase = sqrt(1-cphase*cphase);
      sphi   = sqrt(1-cphi*cphi);

      correlation->data[i] = fabs(
        x1.data[i] * cphase * cphi
        + x2.data[i] * sphase * cphi
        + x3.data[i] * cphase * sphi
        + x4.data[i] * sphase * sphi);
    }

    BankEfficiencySaveVector("BankEfficiency_Correlation.dat", *correlation, sampling);

    /* get the orthonormalized template     */
    BankEfficiencyGetOrthogonalFilter(FilterBCV1);
    corrin.signal2        = *FilterBCV1;

    for (i=0; i<(UINT4)correlation->length; i++)
    {
      cphase = cos(phaseV.data[i]);
      cphi   = cos(phi);
      sphase = sqrt(1-cphase*cphase);
      sphi   = sqrt(1-cphi*cphi);
      templatee.data[i]=FilterBCV1->data[i]*cphase*cphi;
      templatee.data[i]+=corrin.signal2.data[i] * sphase*cphi;
    }
    /* get the orthonormalized template     */
    BankEfficiencyGetOrthogonalFilter(FilterBCV2);
    corrin.signal2        = *FilterBCV2;
    for (i=0; i<(UINT4)correlation->length; i++)
    {
      cphase = cos(phaseV.data[i]);
      cphi   = cos(phi);
      sphase = sqrt(1-cphase*cphase);
      sphi   = sqrt(1-cphi*cphi);
      templatee.data[i]=FilterBCV2->data[i]*sphase*cphi;
      templatee.data[i]+=corrin.signal2.data[i]*sphase*sphi;
    }
    BankEfficiencySaveVector("BankEfficiency_BestTemplate.dat", templatee, sampling);
  }/*end if printTemplate*/

  /* Free memory */
  LALFree(x1.data);
  LALFree(x2.data);
  LALFree(x3.data);
  LALFree(x4.data);

  LALFree(phaseV.data);
  LALFree(templatee.data);
  LALFree(rho1.data);
  LALFree(rho2.data);
  LALFree(rho3.data);
  LALFree(v0.data);
  LALFree(v1.data);
  LALFree(v2.data);
  LALFree(alphaf.data);
  LALFree(thetavVector.data);

  DETATCHSTATUSPTR(status);
  RETURN(status);


}



/*  ****************************************************************************
 *  The Overlap function for BCV templates
 *  ***************************************************************************/



void BankEfficiencyFillProc(
  ProcessParamsTable    *this_proc_param,
  InspiralCoarseBankIn   coarseBankIn,
  RandomInspiralSignalIn randIn,
  UserParametersIn       userParam)
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
  if (coarseBankIn.numFcutTemplates>0)
  {
    ADD_PROCESS_PARAM("int","%d","--bank-number-fcut",
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
      BankEfficiencyGetStringFromGridType(coarseBankIn.gridSpacing));
  ADD_2PROCESS_PARAM("float","%f %f","--bank-eccentricity-range",
      userParam.eccentricBank.min, userParam.eccentricBank.max);
  ADD_PROCESS_PARAM("int","%d","--bank-eccentricity-bins",
      userParam.eccentricBank.bins);
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
      BankEfficiencyGetStringFromNoiseModel(userParam.noiseModel));
  ADD_PROCESS_PARAM("int","%d","--num-seconds",
      userParam.numSeconds);
  ADD_PROCESS_PARAM("int","%d","--multiply-length-by",
      userParam.increaseVector);
  ADD_PROCESS_PARAM("float","%f","--psi0",
      userParam.psi0);
  ADD_PROCESS_PARAM("float","%f","--psi3",
      userParam.psi3);
  ADD_PROCESS_PARAM("float","%f","--sampling",
      coarseBankIn.tSampling);
  ADD_PROCESS_PARAM("string","%s","--simulation-type",
      BankEfficiencyGetStringFromSimulationType(randIn.type));
  ADD_PROCESS_PARAM("float","%f","--signal-amplitude",
      randIn.SignalAmp);
  ADD_PROCESS_PARAM("float","%f","--signal-alpha",
      randIn.param.alpha);
  ADD_PROCESS_PARAM("float","%f","--signal-alpha1",
      randIn.param.alpha1);
  ADD_PROCESS_PARAM("float","%f","--signal-alpha2",
      randIn.param.alpha2);
  ADD_PROCESS_PARAM("int","%d","--signal-amp-order",
      randIn.param.ampOrder);
  ADD_2PROCESS_PARAM("float","%f %f","--signal-eccentricity-range",
      userParam.eccentricSignal.min, userParam.eccentricSignal.max);
  ADD_2PROCESS_PARAM("float","%f %f","--signal-inclination-range",
      randIn.inclinationMin, randIn.inclinationMax);
  ADD_2PROCESS_PARAM("float","%f %f","--signal-polarisation-range",
      randIn.polarisationAngleMin, randIn.polarisationAngleMax);
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
  ADD_2PROCESS_PARAM("float","%f %f","--signal-tau0-range",
      randIn.t0Min, randIn.t0Max);
  ADD_2PROCESS_PARAM("float","%f %f","--signal-tau3-range",
      randIn.tnMin, randIn.tnMax);
  ADD_PROCESS_PARAM("int","%d","--seed",
      userParam.useed);
  ADD_PROCESS_PARAM("string","%s","--signal",
      BankEfficiencyGetStringFromTemplate(userParam.signal));
  ADD_PROCESS_PARAM("int","%d","--signal-order",
      randIn.param.order);
  ADD_PROCESS_PARAM("string","%s","--template",
      BankEfficiencyGetStringFromTemplate(userParam.template));
  ADD_PROCESS_PARAM("int","%d","--template-order",
      coarseBankIn.order);
  ADD_PROCESS_PARAM("int","%d","--template-amp-order",
      coarseBankIn.ampOrder);
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
  ADD_PROCESS_PARAM("string", "%s", "--fast-simulation",
        BankEfficiencyGetStringFromFastSimulation(userParam.fastSimulation));
  ADD_PROCESS_PARAM("int", "%d", "--fast-param1",
        userParam.fastParam1);


  if (userParam.printResultXml){
    ADD_PROCESS_PARAM("float", "%s", "--xml-output"," ");
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
  if (userParam.binaryInjection == BHNS){
    ADD_PROCESS_PARAM("float", "%s", "--bhns-injection", " ");
  }



#undef ADD_PROCESS_PARAM
#undef ADD_2PROCESS_PARAM
}

/* --- convert a number to a string (random case) --- */
CHAR* BankEfficiencyGetStringFromSimulationType(
  INT4 input)
{
  static CHAR this[256];

  switch (input){
    case 0:
      snprintf(this, sizeof(this), "SignalOnly");
      break;
    case 1:
      snprintf(this, sizeof(this), "NoiseOnly");
      break;
    case 2:
      snprintf(this, sizeof(this), "NoiseAndSignal");
      break;
  }

  return this;
}

CHAR* BankEfficiencyGetStringFromFastSimulation(INT4 input)
{
 static CHAR this[256];

  switch (input){
    case None:
      snprintf(this, sizeof(this), "None");
      break;
    case EMatch:
      snprintf(this, sizeof(this), "EMatch");
      break;
    case Heuristic1:
      snprintf(this, sizeof(this), "Heuristic1");
      break;
  }
  return this;
}

/* --- convert a number to a string (template) --- */
CHAR* BankEfficiencyGetStringFromTemplate(INT4 input)
{
  static CHAR this[256];

  switch (input){
    case EOBNR:
      snprintf(this, sizeof(this), "EOBNR");
      break;
    case EOB:
      snprintf(this, sizeof(this), "EOB");
      break;
    case AmpCorPPN:
      snprintf(this, sizeof(this),"AmpCorPPN");
      break;
    case TaylorT1:
      snprintf(this, sizeof(this),"TaylorT1");
      break;
    case Eccentricity:
      snprintf(this, sizeof(this),"Eccentricity");
      break;
    case TaylorT2:
      snprintf(this,  sizeof(this),"TaylorT2");
      break;
    case TaylorT3:
      snprintf(this, sizeof(this), "TaylorT3");
      break;
    case PadeT1:
      snprintf(this, sizeof(this), "PadeT1");
      break;
    case TaylorF2:
      snprintf(this, sizeof(this), "TaylorF2");
      break;
    case BCV:
      snprintf(this, sizeof(this), "BCV");
      break;
    case SpinTaylor:
      snprintf(this, sizeof(this), "SpinTaylor");
      break;
    case TaylorEt:
      snprintf(this, sizeof(this), "TaylorEt");
      break;
    case TaylorN:
      snprintf(this, sizeof(this), "TaylorN");
      break;
    case TaylorT4:
      snprintf(this, sizeof(this), "TaylorT4");
      break;
  }
  return this;
}

/* --- convert a number to a string (bank) --- */
CHAR* BankEfficiencyGetStringFromGridType(INT4 input)
{
  static CHAR this[256];

  switch (input){
  case Square:
    snprintf(this, sizeof(this), "Square");
    break;
  case SquareNotOriented:
    snprintf(this, sizeof(this), "SquareNotOriented");
    break;
  case HexagonalNotOriented:
    snprintf(this,  sizeof(this), "HexagonalNotOriented");
    break;
  case Hexagonal:
    snprintf(this, sizeof(this),"Hexagonal");
    break;
  case HybridHexagonal:
    snprintf(this, sizeof(this),"HybridHexagonal");
    break;
  }

  return this;
}

/* --- convert a number to a string (noise model) --- */
CHAR* BankEfficiencyGetStringFromNoiseModel(
  INT4 input)
{
  static CHAR this[256];

  switch (input){
  case LIGOI:
    snprintf(this, sizeof(this),"LIGOI");
    break;
  case LIGOA:
    snprintf(this, sizeof(this),"LIGOA");
    break;
  case VIRGO:
    snprintf(this, sizeof(this),"VIRGO");
    break;
  case EGO:
    snprintf(this, sizeof(this),"EGO");
    break;
  case GEO:
    snprintf(this, sizeof(this),"GEO");
    break;
  case TAMA:
    snprintf(this, sizeof(this),"TAMA");
    break;
  case UNITY:
    snprintf(this, sizeof(this),"UNITY");
    break;
  case REALPSD:
    snprintf(this, sizeof(this),"REALPSD");
    break;
  }
  return this;
}

/* xml file for the standalone code */
void BankEfficiencyPrintResultsXml(
  InspiralCoarseBankIn         coarseBankIn,
  RandomInspiralSignalIn       randIn,
  UserParametersIn             userParam,
  ResultIn                     trigger,
  BankEfficiencySimulation     simulation)
{
  LALStatus             status = blank_status;
  MetadataTable         templateBank;
  CHAR                  ifo[3];               /* two character ifo code       */
  LIGOLwXMLStream       xmlStream;
  CHAR                  fname[256];
  LIGOTimeGPS           gpsStartTime    = { 0, 0 }; /* input data GPS start time*/
  LIGOTimeGPS           gpsEndTime  = { 0, 0 };     /* input data GPS end time*/
  CHAR                  comment[LIGOMETA_COMMENT_MAX];
  CHAR                  ifoName[MAXIFO][LIGOMETA_IFO_MAX];

  MetadataTable         processParamsTable;
  ProcessParamsTable   *this_proc_param = NULL;

  snprintf( fname, sizeof(fname),
           BANKEFFICIENCY_XMLRESULTS ,
           "simulation",
           gpsStartTime.gpsSeconds,
           gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );

  if (trigger.ntrial == 1)
  {
    strncpy( ifoName[0], "no", LIGOMETA_IFO_MAX * sizeof(CHAR) );
    strncpy( ifoName[1], "ne", LIGOMETA_IFO_MAX * sizeof(CHAR) );
    memset( ifo, 0, sizeof(ifo) );
    memcpy( ifo, "MC", sizeof(ifo) - 1 );

    /* --- we start to fill the xml file here --- */
    memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
    LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlStream, fname), &status );

    /* --- create the process and process params tables --- */
    templateBank.processTable =
        (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
    XLALGPSTimeNow ( &(templateBank.processTable->start_time)) ;

    if (strcmp(CVS_REVISION_C, "$Revi" "sion$"))
    {
      XLALPopulateProcessTable(templateBank.processTable, \
          PROGRAM_NAME, CVS_REVISION_C, CVS_SOURCE_C, CVS_DATE_C, 0);
    }
    else
    {
      XLALPopulateProcessTable(templateBank.processTable, \
          PROGRAM_NAME, lalappsGitCommitID, lalappsGitGitStatus, \
          lalappsGitCommitDate, 0);
    }


    this_proc_param = processParamsTable.processParamsTable =
      (ProcessParamsTable *) calloc( 1, sizeof(ProcessParamsTable) );

    BankEfficiencyFillProc(this_proc_param, coarseBankIn, randIn, userParam);

    memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );

    /* --- write process table --- */
    snprintf( templateBank.processTable->ifos, LIGOMETA_IFOS_MAX, "%s%s",
         ifoName[0], ifoName[1] );
    XLALGPSTimeNow ( &(templateBank.processTable->end_time) );

    /* --- write the template bank --- */
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

    /* finally write sngl inspiral table column's names */

    PRINT_LIGOLW_XML_BANKEFFICIENCY(xmlStream.fp->fp);
  }
  else
  {
    /* --- from the second simulation we just need to append the file --- */
    xmlStream.fp = XLALFileOpenAppend( fname, 0);
  }

  /* --- print the results of one simulation into the xml file --- */
  fprintf(xmlStream.fp->fp,BANKEFFICIENCY_PARAMS_ROW,
    trigger.mass1_trigger,
    trigger.mass2_trigger,
    randIn.param.psi0,
    randIn.param.psi3,
    trigger.tau0_trigger,
    trigger.tau3_trigger,
    randIn.param.t0,
    randIn.param.t3,
    trigger.eccentricity,
    randIn.param.eccentricity,
    randIn.param.alpha1,
    trigger.fend_trigger,
    trigger.fend_inject,
    trigger.mass1_inject,
    trigger.mass2_inject,
    trigger.inclination,
    trigger.polarisationAngle,
    randIn.param.startPhase,
    trigger.rho_final,
    trigger.snrAtCoaTime,
    trigger.phase,
    trigger.alphaF,
    trigger.bin,
    randIn.param.nStartPad,
    trigger.nfast,
    simulation.bestSNRIndex,
    simulation.bestEMatch
    );

  /* --- if we reached the last simulations, we close the file --- */
  if (trigger.ntrial == (UINT4)userParam.ntrials)
  {
    PRINT_LIGOLW_XML_TABLE_FOOTER(xmlStream.fp->fp);
    PRINT_LIGOLW_XML_FOOTER(xmlStream.fp->fp);
    XLALFileClose(xmlStream.fp);
    xmlStream.fp = NULL;
  }
  else
  {
    fprintf(xmlStream.fp->fp, ",\n");
    XLALFileClose(xmlStream.fp);
    xmlStream.fp = NULL;
  }



}


/* print prototype for condor code (option --print-result-prototype)*/
void
BankEfficiencyPrintProtoXml(
  InspiralCoarseBankIn   coarseBankIn,
  RandomInspiralSignalIn randIn,
  UserParametersIn       userParam
        )
{
  LALStatus             status = blank_status;

  MetadataTable         templateBank;
  CHAR                  ifo[3];            /* two character ifo code       */
  LIGOLwXMLStream       xmlStream;
  CHAR                  fname[256];
  /*LIGOTimeGPS           gpsStartTime    = { 0, 0 };
  LIGOTimeGPS           gpsEndTime  = { 0, 0 };    */
  CHAR                  comment[LIGOMETA_COMMENT_MAX];
  CHAR                  ifoName[MAXIFO][LIGOMETA_IFO_MAX];

  MetadataTable         processParamsTable;
  ProcessParamsTable   *this_proc_param = NULL;

  fprintf(stderr,"---%s\n",fname);
  sprintf(fname, BANKEFFICIENCY_PROTOTYPE);
  /*strcat(fname, "_");
  strcat(fname, userParam.tag);
  */
  strcat(fname, ".xml");

  strncpy( ifoName[0], "no", LIGOMETA_IFO_MAX * sizeof(CHAR) );
  strncpy( ifoName[1], "ne", LIGOMETA_IFO_MAX * sizeof(CHAR) );
  memset( ifo, 0, sizeof(ifo) );
  memcpy( ifo, "MC", sizeof(ifo) - 1 );

  /* -- we start to fill the xml file here --- */
  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlStream, fname), &status );


  /* create the process and process params tables */
  templateBank.processTable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
  XLALGPSTimeNow (&(templateBank.processTable->start_time));

  if (strcmp(CVS_REVISION_C, "$Revi" "sion$"))
  {
  XLALPopulateProcessTable(templateBank.processTable, \
	PROGRAM_NAME, CVS_REVISION_C, CVS_SOURCE_C, CVS_DATE_C, 0);
  }
  else
  {
    XLALPopulateProcessTable(templateBank.processTable, \
      PROGRAM_NAME, lalappsGitCommitID, lalappsGitGitStatus, \
	  lalappsGitCommitDate, 0);
  }

  this_proc_param = processParamsTable.processParamsTable =
      (ProcessParamsTable *) calloc( 1, sizeof(ProcessParamsTable) );

  BankEfficiencyFillProc(this_proc_param, coarseBankIn, randIn, userParam);

  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );

  /* write process table */
  snprintf( templateBank.processTable->ifos, LIGOMETA_IFOS_MAX, "%s%s",
         ifoName[0], ifoName[1] );
  XLALGPSTimeNow(&(templateBank.processTable->end_time));

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
  exit( 1 );

}

void BankEfficiencyInitOverlapOutputIn(
  OverlapOutputIn *this)
{
  BankEfficiencyPrintMessage(".");
  this->rhoMax         = -1.;
  this->phase          = -1.;
  this->rhoBin         =  0;
  this->templateNumber =  0;
  this->alpha          = -1.;
  this->layer          =  0;
  this->freq           = -1.;
}

/* Estimate the size of the longest template */
void BankEfficiencyGetMaximumSize(
  LALStatus              *status,
  RandomInspiralSignalIn  randIn,
  InspiralCoarseBankIn    coarseBankIn,
  UserParametersIn        userParam,
  UINT4                  *length)
{
  /* Use these for the longest template from the bank */
  InspiralTemplate params;
  UINT4 maxTmpltLength = 0;

  INITSTATUS( status, "BankEfficiencyGetMaximumSize", BANKEFFICIENCYC );
  ATTATCHSTATUSPTR( status );

  *length = 0;

  /* first the longest signal */
  params            = randIn.param;
  params.massChoice = m1Andm2;
  params.mass1 = params.mass2 = coarseBankIn.mMin;
  LAL_CALL(LALInspiralWaveLength(status->statusPtr, length, params),
       status->statusPtr);

  /* then the longest signal*/
  params.mass1 = params.mass2 = randIn.mMin;
  params.fLower = coarseBankIn.fLower;
  LAL_CALL(LALInspiralWaveLength(status->statusPtr, &maxTmpltLength, params),
           status->statusPtr);

  /* keep longest one */
  if (maxTmpltLength > *length)
    *length = maxTmpltLength;

  /* --m1 --m2 may have been requested. */
  if (userParam.m1 != -1)
  {
  	params.mass1 = userParam.m1;
    params.mass2 = userParam.m2;
    params.fLower = randIn.param.fLower;
    LAL_CALL(LALInspiralWaveLength(status->statusPtr, &maxTmpltLength, params),
           status->statusPtr);
    /* keep longest one */
    if (maxTmpltLength > *length)
      *length = maxTmpltLength;
  }

  if (userParam.tau0 != -1)
  {
  	params.t0 = userParam.tau0;
    params.t3 = userParam.tau3;
    params.massChoice = t03;
    params.fLower = randIn.param.fLower;
    LAL_CALL(LALInspiralParameterCalc( status->statusPtr,  &params ),
        status->statusPtr);

    LAL_CALL(LALInspiralWaveLength(status->statusPtr, &maxTmpltLength, params),
           status->statusPtr);
    /* keep longest one */
    if (maxTmpltLength > *length)
      {
      	*length = maxTmpltLength;
      }
  }



  /* ideally this should be implemented within lal.*/
  if( randIn.param.approximant  == AmpCorPPN || userParam.template == AmpCorPPN)
  {
    *length *= pow(1./(1.+randIn.param.ampOrder/2.), -8./3.);
  }

  /* if --num-seconds is provided, we want a specified length of data.
   * Otherwise, we stick to a minimal size given by twice the length of the
   * maximal template length.*/

  if (userParam.numSeconds * randIn.param.tSampling >= *length) {
    *length = userParam.numSeconds * randIn.param.tSampling;
  }
  else if (userParam.numSeconds != -1) {
    fprintf(stderr,
        "you asked for a signal duration of  %d seconds (--num-seconds) but a (%f %f) system  might be longer than that this duration (%f)...quitting \n",
        userParam.numSeconds, randIn.mMin,
        randIn.mMin, (REAL4)(*length)/(REAL4)randIn.param.tSampling);
    exit(0);
  }
  else {
    userParam.numSeconds = (INT4)(*length/randIn.param.tSampling);
  }

  if (userParam.increaseVector>1)
  {
  	userParam.numSeconds *= userParam.increaseVector;
  	*length *= userParam.increaseVector;
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );

}



void BankEfficiencyCreatePsd(
  LALStatus                *status,
  InspiralCoarseBankIn     *coarseBankIn,
  RandomInspiralSignalIn   *randIn,
  UserParametersIn          userParam)
{
  UINT4 i=0;
  REAL4 df;
  FILE  *Foutput;

  INITSTATUS( status, "BankEfficiencyCreatePsd", BANKEFFICIENCYC );
  ATTATCHSTATUSPTR( status );

  BankEfficiencyPrintMessage("Generating PSD ...");

  memset( &(coarseBankIn->shf), 0, sizeof(REAL8FrequencySeries) );
  coarseBankIn->shf.f0  = 0;
  LAL_CALL(LALDCreateVector( status->statusPtr, &(coarseBankIn->shf.data),
      randIn->psd.length ),
       status->statusPtr);

  coarseBankIn->shf.deltaF  = randIn->param.tSampling /
      ((randIn->psd.length-1)*2);

  /* --- Compute Noise Spectral Density --- */
  df = randIn->param.tSampling/(REAL4) ((randIn->psd.length-1)*2);

  switch (userParam.noiseModel)
  {
    case UNITY:
      LAL_CALL(LALNoiseSpectralDensity (status->statusPtr,
          coarseBankIn->shf.data, &LALLIGOIPsd, df),
          status->statusPtr);
      for (i = 0; i < coarseBankIn->shf.data->length; i++)
    coarseBankIn->shf.data->data[i] = 1;
      break;
    case LIGOI:
      LAL_CALL(LALNoiseSpectralDensity (status->statusPtr,
          coarseBankIn->shf.data, &LALLIGOIPsd, df),
          status->statusPtr);
      break;
    case LIGOA:
      LAL_CALL(LALNoiseSpectralDensity (status->statusPtr,
          coarseBankIn->shf.data, &LALAdvLIGOPsd, df),
          status->statusPtr);
      break;
    case VIRGO:
      LAL_CALL(LALNoiseSpectralDensity (status->statusPtr,
          coarseBankIn->shf.data, &LALVIRGOPsd, df),
          status->statusPtr);;
      break;
    case EGO:
      LAL_CALL(LALNoiseSpectralDensity (status->statusPtr,
          coarseBankIn->shf.data, &LALEGOPsd, df),
          status->statusPtr);;
      break;
    case GEO:
      LAL_CALL(LALNoiseSpectralDensity (status->statusPtr,
          coarseBankIn->shf.data, &LALGEOPsd, df),
          status->statusPtr);
      break;
    case TAMA:
      LAL_CALL(LALNoiseSpectralDensity (status->statusPtr,
          coarseBankIn->shf.data, &LALTAMAPsd, df),
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
      fprintf(Foutput, "%f %e\n",(REAL4)i*df, coarseBankIn->shf.data->data[i]);
    fclose(Foutput);
  }


  /* copy the psd in RandIn.psd */
  for (i = 0; i< coarseBankIn->shf.data->length; i++){
    randIn->psd.data[i] = coarseBankIn->shf.data->data[i];
  }

  BankEfficiencyPrintMessage(" ... done\n");

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

void BankEfficiencyGenerateInputData(
  LALStatus               *status,
  REAL4Vector             *signal,
  RandomInspiralSignalIn  *randIn,
  UserParametersIn         userParam)
{
  UINT4 trial ;
  UINT4 success ;

  INITSTATUS( status, "BankEfficiencyGenerateInputData", BANKEFFICIENCYC );
  ATTATCHSTATUSPTR( status );

  if(vrbflg) fprintf(stderr, "Signal mass1=%e mass2=%e\n", randIn->param.mass1, randIn->param.mass2);
  /* randomize the input start time */
  if (randnStartPad==1)
  {
    randIn->param.nStartPad = (int)( (float)rand() /
        (float)RAND_MAX*signal->length/2);
  }

  trial = 0 ;
  success = 0 ;

  /* we might force to compute a non random waveform giving the two
     input parameters m1 and m2 */
  randIn->param.approximant = userParam.signal;

  if (userParam.startPhase == 1){
    randIn->param.startPhase = BankEfficiencyRandom(0., LAL_PI);
  }
  else
  {
    randIn->param.startPhase = 0.;
  }

  if (userParam.signal==Eccentricity)
  {
    randIn->param.eccentricity = BankEfficiencyRandom(
      userParam.eccentricSignal.min, userParam.eccentricSignal.max);
  }
  else
  randIn->param.eccentricity=0;


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
      }
      else /* random parameters*/
      {
        randIn->param.massChoice = fixedPsi;
        {
          INT4 valid = 0 ;

          while (!valid)
          {
            REAL8 fend = 10000000.;
            INT4 trials = 0;

            randIn->param.psi0 = BankEfficiencyRandom(randIn->psi0Min, randIn->psi0Max);
            randIn->param.psi3 = BankEfficiencyRandom(randIn->psi3Min, randIn->psi3Max);


            while (((fend > randIn->param.tSampling/2.) ||
                ( fend < randIn->param.fLower)) && (trials < 10))
            {
              REAL8 fLR, fLSO;

              randIn->param.totalMass = -randIn->param.psi3/(16.L*LAL_PI
                  * LAL_PI * randIn->param.psi0);
              randIn->param.totalMass *= 2 / LAL_MTSUN_SI;

              fLR = 1.L/(LAL_PI * pow (3.L,1.5) * randIn->param.totalMass
                  * LAL_MTSUN_SI);
              fLSO = 1.L/(LAL_PI * pow (6.L,1.5) * randIn->param.totalMass
                  * LAL_MTSUN_SI);

              fend = BankEfficiencyRandom(fLSO, fLR);

              randIn->param.fFinal = fend;
              randIn->param.fCutoff = fend;

              trials++;
            }
            if (trials == 10 ) valid = 0;
            else valid = 1;
          }
        }
      }
      LALRandomInspiralSignal(status->statusPtr, signal, randIn);
      CHECKSTATUSPTR(status);
    }
    else /* EOB , T1 and so on*/
    {

      randIn->param.inclination = BankEfficiencyRandom(
          randIn->inclinationMin, randIn->inclinationMax);
      randIn->param.polarisationAngle = BankEfficiencyRandom(
          randIn->polarisationAngleMin,randIn->polarisationAngleMax);

      if (randIn->param.approximant == EOBNR)
        randIn->param.order = LAL_PNORDER_PSEUDO_FOUR;

      if (randIn->param.approximant == SpinTaylor)
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
        randIn->sourcePhiMin = 0.;
        randIn->sourcePhiMax = 1.;
        randIn->sourceThetaMin = 0.;
        randIn->sourceThetaMax = 1.;
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


  /* I guess that we can not whatever we want here, because
   * randIn->type==1 means no signal or more precisely
   * noise level is much greater than the signal
   * */
  if (randIn->type == 1)
  {
    randIn->param.massChoice = m1Andm2;
    LALRandomInspiralSignal(status->statusPtr, signal, randIn);
    CHECKSTATUSPTR(status);
  }

  /* --- print the signal waveform --- */
  if (vrbflg)
  {
  BankEfficiencyPrintMessage(" ... done\n");
  BankEfficiencySaveVector(BANKEFFICIENCY_ASCIISIGNAL, *signal,
     randIn->param.tSampling);
  }

  DETATCHSTATUSPTR(status);
  RETURN (status);

}



void BankEfficiencyCreatePowerVector(
  LALStatus                 *status,
  BankEfficiencyPowerVector *powerVector,
  RandomInspiralSignalIn     randIn,
  INT4                       length)
{
  INITSTATUS (status, "BankEfficiencyCreatePowerVector", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);

  powerVector->fm2_3.length = length / 2 ;
  powerVector->fm5_3.length = length / 2;
  powerVector->fm7_6.length = length / 2 ;
  powerVector->fm1_2.length = length / 2;

  powerVector->fm2_3.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) *
      powerVector->fm2_3.length);
  powerVector->fm5_3.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) *
      powerVector->fm5_3.length);
  powerVector->fm7_6.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) *
      powerVector->fm7_6.length);
  powerVector->fm1_2.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) *
      powerVector->fm1_2.length);

  /* --- Create useful vector once for all --- */
  LAL_CALL(BankEfficiencyCreateVectorFreqPower(&powerVector->fm5_3,
      randIn.param, -5, 3), status->statusPtr);
  LAL_CALL(BankEfficiencyCreateVectorFreqPower(&powerVector->fm2_3,
      randIn.param, -2, 3), status->statusPtr);
  LAL_CALL(BankEfficiencyCreateVectorFreqPower(&powerVector->fm7_6,
      randIn.param, -7, 6), status->statusPtr);
  LAL_CALL(BankEfficiencyCreateVectorFreqPower(&powerVector->fm1_2,
      randIn.param, -1, 2), status->statusPtr);

  DETATCHSTATUSPTR( status );
  RETURN( status );
}



void BankEfficiencyInspiralOverlapBCV(
  LALStatus  *status,
  InspiralTemplate          *list,
  UserParametersIn           userParam,
  RandomInspiralSignalIn    *randIn,
  InspiralWaveOverlapIn     *overlapin,
  OverlapOutputIn           *overlapout,
  REAL4Vector               *correlation,
  BankEfficiencyBCV         *bankefficiencyBCV)
{
  REAL4 fendBCV ;
  REAL4 df;
  INT4 n, kMin, kMax;

  INITSTATUS (status, "BankEfficiencyInspiralOverlapBCV", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);

  n       = bankefficiencyBCV->FilterBCV1.length;
  df      = randIn->param.tSampling / (REAL8)(n/2) / 2.;
  kMin    = floor( randIn->param.fLower /df );
  fendBCV = list->fFinal; /* ending frequency of template */

  /* if we want to test the BCV metric by injectd BCV signal, we
   * have to get the same final frequnecy for both template and signal */
  if (userParam.signal == BCV){
    fendBCV = list->fFinal = randIn->param.fFinal;
  }


  if ( fendBCV < randIn->param.fLower*1.5 )
  {
  	/*fprintf(stderr, "Warning. fendBCV less than 1.5*fLower. Forced to be 1.5*Lower");
  	fprintf(stderr, "psi0=%f psi3=%f ffinal=%f\n",list->psi0, list->psi3,fendBCV);
  	*/fendBCV = 1.5 * randIn->param.fLower;
  	list->fFinal = 1.5 * randIn->param.fLower;
  }
  overlapin->ifExtOutput    = 0;
  overlapin->param.fFinal   =  fendBCV;
  overlapin->param.fCutoff  =  fendBCV;
  overlapin->param          = *list;
  overlapout->rhoMax  = -1;

  kMax = floor(fendBCV / df);

  BankEfficiencyCompare(kMin, kMax,
      "of index used in BCV final frequency below fl. Change psi's ranges");

  /* Template creation */
  LAL_CALL(BankEfficiencyCreateBCVFilters(bankefficiencyBCV,
                kMin, kMax, list->psi0, list->psi3
                ), status->statusPtr);

  /* The overlap given two filters and the input signal*/
  /*BankEfficiencyInitOverlapOutputIn(&OverlapOutputThisTemplate); */
  /*be sure rhomax = 0 and so on*/
  LAL_CALL(BankEfficiencyWaveOverlapBCV(status->statusPtr,
                 correlation, overlapin,
                 &(bankefficiencyBCV->FilterBCV1),
                 &(bankefficiencyBCV->FilterBCV2),
                 userParam,
                 overlapout,
                 &(bankefficiencyBCV->moments)), status->statusPtr);

  /*overlapoutput->snrAtCoaTime = correlation->data[startPad];*/


  DETATCHSTATUSPTR( status );
  RETURN( status );

}

void BankEfficiencyBankPrintAscii(
  MetadataTable        templateBank ,
  UINT4                numCoarse,
  InspiralCoarseBankIn coarseBankIn )
{
  FILE *output;
  SnglInspiralTable *tmpltCurrent  = NULL;

  output = fopen(BANKEFFICIENCY_PRINTBANK_FILEASCII, "w");

  fprintf(output,"#Number of Coarse Bank Templates=%d\n",numCoarse);
  if (coarseBankIn.approximant == BCV)
  {
    fprintf(output, "#psi0Min=%e, psi0Max=%e, psi3Min=%e, psi3Max=%e\n",
        coarseBankIn.psi0Min, coarseBankIn.psi0Max,
        coarseBankIn.psi3Min, coarseBankIn.psi3Max);
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
          tmpltCurrent->tau0, tmpltCurrent->tau3,
          tmpltCurrent->mass1, tmpltCurrent->mass2);

    }
  }

  fprintf(output,"&\n");
  fclose(output);

}


void BankEfficiencyBankPrintXML(
  MetadataTable          templateBank ,
  InspiralCoarseBankIn   coarseBankIn,
  RandomInspiralSignalIn randIn,
  UserParametersIn       userParam)

{
  LALStatus             status = blank_status;
  MetadataTable         proctable;
  CHAR                  ifo[3];         /* two character ifo code       */
  LIGOLwXMLStream       xmlStream;
  CHAR                  fname[256];
  LIGOTimeGPS gpsStartTime  = { 0, 0 }; /* input data GPS start time    */
  LIGOTimeGPS gpsEndTime    = { 0, 0 }; /* input data GPS end time      */
  CHAR  comment[LIGOMETA_COMMENT_MAX];
  CHAR  ifoName[3][LIGOMETA_IFO_MAX];
  MetadataTable         processParamsTable;
  ProcessParamsTable   *this_proc_param = NULL;


  strncpy( ifoName[0], "no", LIGOMETA_IFO_MAX * sizeof(CHAR) );
  strncpy( ifoName[1], "ne", LIGOMETA_IFO_MAX * sizeof(CHAR) );
  memset( ifo, 0, sizeof(ifo) );
  memset( comment, 0, sizeof( comment ) );
  memset( &proctable, 0, sizeof( proctable ) );
  memset( &processParamsTable, 0, sizeof( processParamsTable ) );
  memcpy( ifo, "MC", sizeof(ifo) - 1 );
  

  /* --- first we create the filename --- */
  snprintf( fname, sizeof(fname), BANKEFFICIENCY_XMLBANK ,
           ifo, gpsStartTime.gpsSeconds,
           gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );

  /* -- we start to fill the xml file here --- */
  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );

  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlStream, fname), &status );

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );

  XLALGPSTimeNow (&(proctable.processTable->start_time));

  if (strcmp(CVS_REVISION_C, "$Revi" "sion$"))
  {
    XLALPopulateProcessTable(proctable.processTable, \
    PROGRAM_NAME, CVS_REVISION_C, CVS_SOURCE_C, CVS_DATE_C, 0);
  }
  else
  {
    XLALPopulateProcessTable(proctable.processTable, \
      PROGRAM_NAME, lalappsGitCommitID, lalappsGitGitStatus, \
      lalappsGitCommitDate, 0);
  }

  this_proc_param = processParamsTable.processParamsTable =
    (ProcessParamsTable *) calloc( 1, sizeof(ProcessParamsTable) );

  BankEfficiencyFillProc(this_proc_param, coarseBankIn, randIn, userParam);

  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );



  /* write process table */
  snprintf( proctable.processTable->ifos, LIGOMETA_IFOS_MAX, "%s%s",
           ifoName[0], ifoName[1] );
  XLALGPSTimeNow(&(proctable.processTable->end_time));

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
      LAL_CALL(LALBeginLIGOLwXMLTable(&status, &xmlStream, sngl_inspiral_table),
        &status );
      LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, templateBank,
                    sngl_inspiral_table ), &status );
      LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );
    }

  /* close the output xml file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlStream ), &status );
}

void BankEfficiencyCreateListFromTmplt(
  LALStatus         *status,
  InspiralTemplate  *params,
  Mybank             mybank,
  INT4               index)
{
  INITSTATUS(status, "BankEfficiencyCreateListFromTmplt", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);
  /* Right now only t03 and psi0psi3 bank exists so we
     only need those information in principle and nothing else*/
  params->mass1  = mybank.mass1[index];
  params->mass2  = mybank.mass2[index];
  params->fFinal = mybank.fFinal[index];
  params->t0     = mybank.tau0[index];
  params->t3     = mybank.tau3[index];
  params->psi0   = mybank.psi0[index];
  params->psi3   = mybank.psi3[index];
  params->eccentricity  = mybank.eccentricity[index];

  DETATCHSTATUSPTR(status);
}


/* ****************************************************************************
 * Initialization of the parameters fromn random signal, bank and others
 * structures
 * ***************************************************************************/
void BankEfficiencyParametersInitialization(
  InspiralCoarseBankIn   *coarseBankIn,
  RandomInspiralSignalIn *randIn,
  UserParametersIn       *userParam)
{
  BankEfficiencyInitInspiralCoarseBankIn(coarseBankIn);
  BankEfficiencyInitRandomInspiralSignalIn(randIn);
  BankEfficiencyInitUserParametersIn(userParam);
}


/* ****************************************************************************
 * CoarseIn initialization (bank and template)
 ***************************************************************************/
void BankEfficiencyInitInspiralCoarseBankIn(
  InspiralCoarseBankIn *coarseBankIn)
{
  coarseBankIn->fLower           = 40.;
  coarseBankIn->fUpper           = -1.;
  coarseBankIn->tSampling        = 2048;
  coarseBankIn->space            = Tau0Tau3;
  coarseBankIn->mmCoarse         = 0.80;
  coarseBankIn->mmFine           = 0.90;
  coarseBankIn->iflso            = 0.;
  coarseBankIn->mMin             = 5.;
  coarseBankIn->mMax             = 20.;
  coarseBankIn->MMax             = -1;
  coarseBankIn->massRange        = MinMaxComponentMass;
  coarseBankIn->etamin           = -1;
  coarseBankIn->psi0Min          = 10.;
  coarseBankIn->psi0Max          = 250000.;
  coarseBankIn->psi3Min          = -2200.;
  coarseBankIn->psi3Max          = -10.;
  coarseBankIn->alpha            = 0.01;
  coarseBankIn->numFcutTemplates = 5;
  coarseBankIn->approximant      = EOB;
  coarseBankIn->order            = LAL_PNORDER_TWO;
  coarseBankIn->LowGM            = 3;
  coarseBankIn->insidePolygon    = 1;
  coarseBankIn->HighGM           = 6;
  coarseBankIn->gridSpacing      = SquareNotOriented;
  coarseBankIn->computeMoments   = 0;
  coarseBankIn->numFreqCut       = 1;
  coarseBankIn->maxFreqCut       = FreqCut_SchwarzISCO;
  coarseBankIn->minFreqCut       = FreqCut_SchwarzISCO;
}


/* ****************************************************************************
 * Random initialization
 * ****************************************************************************/
void BankEfficiencyInitRandomInspiralSignalIn(
  RandomInspiralSignalIn *randIn)
{
  randIn->useed                 = 122888;  /* seed for MonteCarlo             */
  randIn->type                  = 0;       /* type of simulation: signal only */
  randIn->SignalAmp             = 10;      /* SNR of the signal if noise exits*/
  randIn->param.order           = LAL_PNORDER_TWO; /* and its order           */
  randIn->param.alpha           = 0;       /* alpha paramBCV                  */
  randIn->param.ieta            = 1;       /*                                 */
  randIn->param.mass1           =-1;       /* To allocate memory              */
  randIn->param.mass2           =-1;       /* idem                            */
  randIn->param.fLower          = 40.;     /* Lower freq.(templ)              */
  randIn->param.OmegaS          = 0.;      /* EOB parameter                   */
  randIn->param.Theta           = 0.;      /* EOB parameter                   */
  randIn->mMin                  = 5;       /* min mass to inject              */
  randIn->mMax                  = 20;      /* max mass to inject              */
  randIn->MMax                  = randIn->mMax * 2;  /* total mass max        */
  randIn->t0Min                 = 0;       /* min tau0 to inject              */
  randIn->t0Max                 = 0.;      /* max tau0 to inject              */
  randIn->tnMin                 = 0.1;     /* min tau3 to inject              */
  randIn->tnMax                 = 1.;      /* max tau3 to inject              */
  randIn->etaMin                = randIn->mMin* (randIn->mMax - randIn->mMax)
    / (randIn->MMax ) / (randIn->MMax);
  randIn->psi0Min               = 10.;     /* psi0 range                      */
  randIn->psi0Max               = 250000.;
  randIn->psi3Min               = -2200;
  randIn->psi3Max               = -10;
  randIn->param.approximant     = EOB;     /* signal approximant              */
  randIn->param.tSampling       = 2048;    /* sampling                        */
  randIn->param.fCutoff         = 0.;
  randIn->param.startTime       = 0;
  randIn->param.startPhase      = 0.;
  randIn->param.nStartPad       = 1000;
  randIn->param.signalAmplitude = 10;
  randIn->param.nEndPad         = 0;
  randIn->NoiseAmp              = 1;
  randIn->param.eccentricity    = 0.;
  randIn->param.ampOrder        = 5;
  randIn->inclinationMin        = 0.;
  randIn->inclinationMax        = LAL_PI;
  randIn->polarisationAngleMin  = 0.;
  randIn->polarisationAngleMax  = LAL_PI;
  randIn->param.alpha1                = 0; /*used to populate eccentriciyt at freq=fLower*/
}


/* ****************************************************************************
 * Other Param initialization
 * ***************************************************************************/
void BankEfficiencyInitUserParametersIn(
  UserParametersIn *userParam)
{
  userParam->BankFromFile        = 0;
  userParam->eccentricBank.min   = 0.;
  userParam->eccentricBank.max   = 0.1;
  userParam->eccentricBank.bins  = 1;
  userParam->eccentricSignal.min = 0.;
  userParam->eccentricSignal.max = 0.4;
  userParam->alphaFConstraint    = 1;
  userParam->extraFinalPrinting  = 0;
  userParam->eMatch              = 0.5;
  userParam->template            = EOB;
  /*By default, this value is empty and will be populate with the sampling
   * frequency later*/
  userParam->signalfFinal        =  0.;
  userParam->signal              = EOB;
  userParam->m1                  = -1;
  userParam->m2                  = -1;
  userParam->numSeconds          = -1;
  userParam->increaseVector      =  1;
  userParam->psi0                = -1;
  userParam->psi3                = -1;
  userParam->tau0                = -1;
  userParam->tau3                = -1;
  userParam->t0FineMin           = 0;
  userParam->t0FineMax           = 0;
  userParam->t3FineMin           = 0;
  userParam->t3FineMax           = 0;
  userParam->t0FineBin           = 0;
  userParam->t3FineBin           = 0;
  userParam->printBestOverlap    = BANKEFFICIENCY_PRINTBESTOVERLAP;
  userParam->printBestTemplate   = BANKEFFICIENCY_PRINTBESTTEMPLATE;
  userParam->printSNRHisto       = BANKEFFICIENCY_PRINTSNRHISTO;
  userParam->printPsd            = BANKEFFICIENCY_PRINTPSD;
  userParam->printBank           = BANKEFFICIENCY_PRINTBANK;
  userParam->printResultXml      = BANKEFFICIENCY_PRINTRESULTXML;
  userParam->printPrototype      = BANKEFFICIENCY_PRINTPROTOTYPE;
  userParam->faithfulness        = BANKEFFICIENCY_FAITHFULNESS;
  userParam->ntrials             = 1;
  userParam->fastSimulation      = None;
  userParam->noiseModel          = LIGOI;
  userParam->binaryInjection     = NoUserChoice;
  userParam->maxTotalMass        = -1;
  userParam->minTotalMass        = -1;
  userParam->startPhase          = 1;
  userParam->inputPSD            = NULL;
  userParam->ambiguity           = 0;
  userParam->fastParam1          = 50;
  sprintf(userParam->tag, "EOBNR");
}



/* ***************************************************************************
 *  Function to parse user parameters
 *  ************************************************************************ */
void BankEfficiencyParseParameters(
  INT4                    *argc,
  CHAR                   **argv,
  InspiralCoarseBankIn    *coarseBankIn,
  RandomInspiralSignalIn  *randIn,
  UserParametersIn        *userParam)
{
  INT4 i = 1;
  REAL8 tmp1, tmp2;

  if (*argc==1)
  {
    BankEfficiencyHelp();
    fprintf(stderr,
      "You must provide at least the template and signal argument");
    exit(0);
  }

  while(i < *argc)
  {
    if (!strcmp(argv[i],    "--bank-alpha")) {
      BankEfficiencyParseGetDouble(argv, &i, &(coarseBankIn->alpha));
    }
    else if (!strcmp(argv[i],"--bank-fcut-range")) {
      BankEfficiencyParseGetDouble2(argv, &i, &tmp1, &tmp2);
      coarseBankIn->LowGM  = (REAL4)tmp1;
      coarseBankIn->HighGM = (REAL4)tmp2;
    }
    else if (!strcmp(argv[i],   "--bank-ffinal")) {
      BankEfficiencyParseGetDouble(argv, &i, &(coarseBankIn->fUpper));
    }
    else if (!strcmp(argv[i],   "--bank-inside-polygon")) {
      INT4 temp;
      BankEfficiencyParseGetInt(argv, &i, &temp);
      coarseBankIn->insidePolygon = (InsidePolygon)temp;
    }
    else if (!strcmp(argv[i], "--bank-grid-spacing")) {
      BankEfficiencyParseGetString(argv, &i);
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
        fprintf(stderr,
            "Unknown grid spacing after --bank-grid-spacing. Try --help.\n");
        exit(0);
      }
    }
    else if (!strcmp(argv[i],   "--bank-number-fcut")) {
      BankEfficiencyParseGetInt(argv, &i,
          (INT4*)&(coarseBankIn->numFcutTemplates));
    }
    else if (!strcmp(argv[i],"--bank-mass-range")) {
      BankEfficiencyParseGetDouble2(argv, &i,
           &(coarseBankIn->mMin), &(coarseBankIn->mMax));
    }
    else if (!strcmp(argv[i],"--bank-max-total-mass")) {
      BankEfficiencyParseGetDouble(argv, &i, &(coarseBankIn->MMax));
    }
    else if (!strcmp(argv[i],"--bank-min-total-mass")) {
      BankEfficiencyParseGetDouble(argv, &i, &(coarseBankIn->MMin));
    }
    else if (!strcmp(argv[i],"--bank-from-file")) {
      userParam->BankFromFile = 1;
      BankEfficiencyParseGetString(argv, &i);
      userParam->BankFile = argv[i];
    }
    else if (!strcmp(argv[i],   "--bank-psi0-range")) {
      BankEfficiencyParseGetDouble2(argv, &i,
          &(coarseBankIn->psi0Min), &(coarseBankIn->psi0Max));
    }
    else if (!strcmp(argv[i],   "--bank-psi3-range")) {
      BankEfficiencyParseGetDouble2(argv, &i,
          &(coarseBankIn->psi3Min), &(coarseBankIn->psi3Max));
    }
    else if (!strcmp(argv[i],"--debug")) {
      BankEfficiencyParseGetInt(argv, &i, &(lalDebugLevel));
    }
    else if (!strcmp(argv[i],"--bank-eccentricity-range")){
      BankEfficiencyParseGetDouble2(argv, &i,
          &(userParam->eccentricBank.min), &(userParam->eccentricBank.max));
    }
    else if (!strcmp(argv[i],"--bank-eccentricity-bins")){
      BankEfficiencyParseGetInt(argv, &i ,&(userParam->eccentricBank.bins));
    }
    else if (!strcmp(argv[i], "--e-match")){
      BankEfficiencyParseGetDouble(argv, &i, &(userParam->eMatch));
    }
    else if (!strcmp(argv[i],"--fast-simulation")) {
      BankEfficiencyParseGetString(argv, &i);
      if (!strcmp(argv[i], "EMatch"))
        userParam->fastSimulation = EMatch;
      else if (!strcmp(argv[i], "Heuristic1"))
        userParam->fastSimulation = Heuristic1;
      else if (!strcmp(argv[i], "None"))
        userParam->fastSimulation = None;

      else
      {
        fprintf(stderr, "Wrong argument. Use either EMatch or Heuristic1 or None)\n");
        exit( 1 );
      }
    }
    else if (!strcmp(argv[i],"--signal-fl")) {
      BankEfficiencyParseGetDouble(argv, &i, &(randIn->param.fLower));
    }
    else if  (!strcmp(argv[i],"--template-fl")) {
      BankEfficiencyParseGetDouble(argv, &i, &(coarseBankIn->fLower));
    }
    else if  (!strcmp(argv[i],"--fast-param1")) {
      BankEfficiencyParseGetInt(argv, &i, &(userParam->fastParam1));
    }
    else if  (!strcmp(argv[i],"--fl")) {
      BankEfficiencyParseGetDouble(argv, &i, &(coarseBankIn->fLower));
      randIn->param.fLower = coarseBankIn->fLower;
    }
    else if (!strcmp(argv[i], "--help") || !strcmp(argv[i],"--h")) {
      BankEfficiencyHelp();
    }
    else if (!strcmp(argv[i],"--signal-max-total-mass")) {
      BankEfficiencyParseGetDouble(argv, &i, &tmp1);
      userParam->maxTotalMass = (REAL4)tmp1;
    }
    else if (!strcmp(argv[i],"--signal-min-total-mass")) {
      BankEfficiencyParseGetDouble(argv, &i, &tmp1);
      userParam->minTotalMass = (REAL4)tmp1;
    }
    else if (!strcmp(argv[i],"--ascii2xml")) {
      ascii2xml = 1;
      /* we do not need to parse the l10
       * other arguments */
      i = *argc+1;
    }
    else if (!strcmp(argv[i], "--m1")){
      BankEfficiencyParseGetDouble(argv, &i, &(userParam->m1));
    }
    else if (!strcmp(argv[i], "--m2")) {
      BankEfficiencyParseGetDouble(argv, &i, &(userParam->m2));
    }
    else if (!strcmp(argv[i],   "--mm")){
      BankEfficiencyParseGetDouble(argv, &i, &(coarseBankIn->mmCoarse));
    }
    else if (!strcmp(argv[i], "--multiply-length-by")) {
      BankEfficiencyParseGetInt(argv, &i, &(userParam->increaseVector));
    }
    else if (!strcmp(argv[i],   "--n") || !strcmp(argv[i], "--ntrial")){
      BankEfficiencyParseGetInt(argv, &i, &(userParam->ntrials));
    }
    else if (!strcmp(argv[i],   "--noise-amplitude")) {
      BankEfficiencyParseGetDouble(argv, &i, &(randIn->NoiseAmp));
    }
    else if (!strcmp(argv[i],   "--noise-model")){
      BankEfficiencyParseGetString(argv, &i);
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
        fprintf(stderr,
            "Unknown moise model after --noise-model. Try --help\n");
        exit(1);
      }
    }
    else if (!strcmp(argv[i], "--num-seconds")) {
      BankEfficiencyParseGetInt(argv, &i, &(userParam->numSeconds));
    }
    else if (!strcmp(argv[i], "--psi0")) {
      BankEfficiencyParseGetDouble(argv, &i, &(userParam->psi0));
    }
    else if (!strcmp(argv[i], "--psi3")) {
      BankEfficiencyParseGetDouble(argv, &i, &(userParam->psi3));
    }
    else if (!strcmp(argv[i],"--sampling")) {
      BankEfficiencyParseGetDouble(argv, &i, &(coarseBankIn->tSampling));
      randIn->param.tSampling = coarseBankIn->tSampling;
    }
    else if (!strcmp(argv[i],"--seed")){
      BankEfficiencyParseGetInt(argv, &i, &(randIn->useed));
    }
    else if (!strcmp(argv[i], "--signal")) {
      BankEfficiencyParseGetString(argv, &i);
      if (strcmp(argv[i], "TaylorT1") == 0)
        userParam->signal = TaylorT1;
      else if (strcmp(argv[i],"AmpCorPPN") == 0)
        userParam->signal = AmpCorPPN;
      else if (strcmp(argv[i],"Eccentricity") == 0)
        userParam->signal = Eccentricity;
      else if (strcmp(argv[i],"TaylorT2") == 0)
        userParam->signal = TaylorT2;
      else if (strcmp(argv[i],"TaylorT3") == 0)
        userParam->signal = TaylorT3;
      else if (strcmp(argv[i],"TaylorF1") == 0)
        userParam->signal = TaylorF1;
      else if (strcmp(argv[i],"TaylorF2") == 0)
        userParam->signal = TaylorF2;
      else if (strcmp(argv[i],"PadeT1") == 0)
        userParam->signal = PadeT1;
      else if (strcmp(argv[i],"PadeF1") == 0)
        userParam->signal = PadeF1;
      else if (strcmp(argv[i],"EOB") == 0)
        userParam->signal = EOB;
      else if (strcmp(argv[i],"EOBNR") == 0)
        userParam->signal = EOBNR;
      else if (strcmp(argv[i],"BCV") == 0)
        userParam->signal = BCV;
      else if (strcmp(argv[i],"SpinTaylor") == 0)
        userParam->signal = SpinTaylor;
      else if (strcmp(argv[i],"TaylorEt") == 0)
        userParam->signal = TaylorEt;
      else if (strcmp(argv[i],"TaylorN") == 0)
        userParam->signal = TaylorN;
      else if (strcmp(argv[i],"TaylorT4") == 0)
        userParam->signal = TaylorT4;
      else
      {
        fprintf(stderr,
         "Wrong approximant. use TaylorT1, TaylorT2, TaylorT3, EOB, PadeT1, BCV,\n\t SpinTaylor, AmpCorPPN, TaylorEt TaylorT4 TaylorN, Eccentricity)\n");
        exit( 1 );
      }
      randIn->param.approximant = userParam->signal;
    }
    else if (!strcmp(argv[i],   "--signal-alpha")){
      BankEfficiencyParseGetDouble(argv, &i, &(randIn->param.alpha));
    }
    else if (!strcmp(argv[i],   "--signal-amp-order")){
      INT4 temp;
      BankEfficiencyParseGetInt(argv, &i, &temp);
      randIn->param.ampOrder = (LALPNOrder)temp;
    }
    else if (!strcmp(argv[i],   "--signal-alpha1")){
      BankEfficiencyParseGetDouble(argv, &i, &(randIn->param.alpha1));
    }
    else if (!strcmp(argv[i],   "--signal-alpha2")){
      BankEfficiencyParseGetDouble(argv, &i, &(randIn->param.alpha2));
    }
    else if (!strcmp(argv[i],   "--signal-amplitude")) {
      BankEfficiencyParseGetDouble(argv, &i, &(randIn->SignalAmp));
    }
    else if (!strcmp(argv[i],"--signal-eccentricity-range")){
      BankEfficiencyParseGetDouble2(argv, &i,
          &(userParam->eccentricSignal.min), &(userParam->eccentricSignal.max));
    }
    else if (!strcmp(argv[i],   "--signal-ffinal")) {
      BankEfficiencyParseGetDouble(argv, &i, &(randIn->param.fCutoff));
      userParam->signalfFinal = randIn->param.fCutoff ;
    }
    else if (!strcmp(argv[i],"--signal-inclination-range")){
      BankEfficiencyParseGetDouble2(argv, &i,
          &(randIn->inclinationMin), &(randIn->inclinationMax));
    }
    else if (!strcmp(argv[i],"--signal-polarisation-range")){
      BankEfficiencyParseGetDouble2(argv, &i,
          &(randIn->polarisationAngleMin), &(randIn->polarisationAngleMax));
    }
    else if (!strcmp(argv[i],   "--signal-mass-range")){
      BankEfficiencyParseGetDouble2(argv, &i,
          &(randIn->mMin), &(randIn->mMax));
      randIn->param.mass1 = randIn->param.mass2 = randIn->mMin;
      /*default values for injection*/
    }
    else if (!strcmp(argv[i], "--signal-tau0-range")){
      BankEfficiencyParseGetDouble2(argv, &i,
          &(randIn->t0Min), &(randIn->t0Max));
    }
    else if (!strcmp(argv[i], "--signal-tau3-range")){
      BankEfficiencyParseGetDouble2(argv, &i,
          &(randIn->tnMin), &(randIn->tnMax));
    }
    else if (!strcmp(argv[i],   "--signal-order")) {
      INT4 temp;
      BankEfficiencyParseGetInt(argv, &i, &temp);
      randIn->param.order = (LALPNOrder)temp;
    }
    else if (!strcmp(argv[i],   "--signal-nstartpad")) {
      BankEfficiencyParseGetInt(argv, &i, (INT4*)&(randIn->param.nStartPad));
    }
    else if (!strcmp(argv[i],   "--signal-nendpad")) {
      BankEfficiencyParseGetInt(argv, &i, (INT4*)&(randIn->param.nEndPad));
    }
    else if (!strcmp(argv[i],   "--signal-psi0-range")) {
      BankEfficiencyParseGetDouble2(argv, &i,
          &(randIn->psi0Min), &(randIn->psi0Max));
    }
    else if (!strcmp(argv[i],   "--signal-psi3-range")){
      BankEfficiencyParseGetDouble2(argv, &i,
          &(randIn->psi3Min), &(randIn->psi3Max));
    }
    else if (!strcmp(argv[i],   "--signal-random-nstartpad")) {
     randnStartPad = 1;
    }
    else if (!strcmp(argv[i],   "--simulation-type")) {
      BankEfficiencyParseGetString(argv, &i);
      if (!strcmp(argv[i],"NoiseOnly"))
        randIn->type = 1;
      else if (!strcmp(argv[i],"SignalOnly"))
        randIn->type = 0;
      else if (!strcmp(argv[i],"NoiseAndSignal"))
        randIn->type = 2;
      else
      {
        fprintf(stderr,
            "unknow simulation-type after --simulation-type. try --help \n");
        exit(1);
      }
    }
    else if (!strcmp(argv[i], "--t0-fine-range")){
      BankEfficiencyParseGetDouble2(argv, &i,
          &(userParam->t0FineMin), &(userParam->t0FineMax));
    }
    else if (!strcmp(argv[i], "--t3-fine-range")){
      BankEfficiencyParseGetDouble2(argv, &i,
          &(userParam->t3FineMin), &(userParam->t3FineMax));
    }
    else if (!strcmp(argv[i], "--t0-fine-bin")){
      BankEfficiencyParseGetDouble(argv, &i, &(userParam->t0FineBin));
    }
    else if (!strcmp(argv[i], "--t3-fine-bin")){
      BankEfficiencyParseGetDouble(argv, &i, &(userParam->t3FineBin));
    }
    else if (!strcmp(argv[i], "--no-start-phase")){
      userParam->startPhase = 0;
    }
    else if (!strcmp(argv[i], "--tau0")) {
      BankEfficiencyParseGetDouble(argv, &i, &(userParam->tau0));
    }
    else if (!strcmp(argv[i], "--tau3")) {
      BankEfficiencyParseGetDouble(argv, &i, &(userParam->tau3));
    }
    else if (!strcmp(argv[i],"--template")) {
      BankEfficiencyParseGetString(argv, &i);

      if (!strcmp(argv[i], "TaylorT1"))
        userParam->template = TaylorT1;
      else if (!strcmp(argv[i], "TaylorT2"))
        userParam->template = TaylorT2;
      else if (!strcmp(argv[i] ,"Eccentricity"))
        userParam->template = Eccentricity;
      else if (!strcmp(argv[i], "AmpCorPPN"))
        userParam->template = AmpCorPPN;
      else if (!strcmp(argv[i], "TaylorT3"))
        userParam->template = TaylorT3;
      else if (!strcmp(argv[i], "TaylorF1"))
        userParam->template = TaylorF1;
      else if (!strcmp(argv[i], "TaylorF2"))
        userParam->template = TaylorF2;
      else if (!strcmp(argv[i], "PadeT1"))
        userParam->template = PadeT1;
      else if (!strcmp(argv[i], "PadeF1"))
        userParam->template = PadeF1;
      else if (!strcmp(argv[i], "EOB"))
        userParam->template = EOB;
      else if (!strcmp(argv[i], "TaylorEt"))
        userParam->template = TaylorEt;
      else if (!strcmp(argv[i], "TaylorN"))
        userParam->template = TaylorN;
      else if (!strcmp(argv[i], "TaylorT4"))
        userParam->template = TaylorT4;
      else if (!strcmp(argv[i], "EOBNR"))
        {
          userParam->template = EOBNR;
          coarseBankIn->order = LAL_PNORDER_PSEUDO_FOUR;
        }
      else if (!strcmp(argv[i], "BCV"))
        userParam->template = BCV;
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
    else if (!strcmp(argv[i],   "--template-order")) {
      INT4 temp;
      BankEfficiencyParseGetInt(argv, &i, &temp);
      coarseBankIn->order = (LALPNOrder)temp;
    }
    else if (!strcmp(argv[i],   "--template-amp-order")) {
      INT4 temp;
      BankEfficiencyParseGetInt(argv, &i, &temp);
      coarseBankIn->ampOrder = (LALPNOrder)temp;
    }
    else if (!strcmp(argv[i], "--user-tag")) {
      BankEfficiencyParseGetString(argv, &i);
      sprintf(userParam->tag,argv[i]);
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
      userParam->alphaFConstraint   = ALPHAFUnconstraint;
    }
    else if (!strcmp(argv[i],"--print-best-overlap")) {
      userParam->printBestOverlap   = 1;
      }
    else if (!strcmp(argv[i],"--faithfulness")) {
      userParam->faithfulness   = 1;
    }
    else if (!strcmp(argv[i],"--check")) {
      userParam->faithfulness   = 1;
    }
    else if (!strcmp(argv[i],"--print-psd")) {
      userParam->printPsd = 1;
    }
    else if (!strcmp(argv[i],"--print-snr-histo")) {
      userParam->printSNRHisto = 1;
    }
    else if (!strcmp(argv[i],"--verbose")) {
      vrbflg = 1;
    }
    else if (!strcmp(argv[i],"--version")) {
      fprintf(stderr, "BankEfficiency code"
          "Thomas Cokelaer, Thomas.Cokelaer@astro.cf.ac.uk\n"
          "CVS Version source:" CVS_ID_STRING_C "\n"
          "CVS Version include:" CVS_ID_STRING "\n"
          "CVS Tag: " CVS_NAME_STRING "\n");
      exit(1);
    }
    else if (!strcmp(argv[i],"--print-bank")) {
      userParam->printBank = 1;
    }
    else if (!strcmp(argv[i],"--compute-moments")) {
      coarseBankIn->computeMoments = 1;
    }
    else if (!strcmp(argv[i],"--xml-output")) {
      userParam->printResultXml = 1;
    }
    else if (!strcmp(argv[i],"--print-prototype")) {
      userParam->printPrototype = 1;
    }
    else {
      fprintf(stderr,
          "%s option does not exist. Type --help for options\n", argv[i]);
      exit(1);
    }
    i++;
  }
}


void BankEfficiencyParseGetInt(
  CHAR    **argv,
  INT4    *index,
  INT4    *data)
{
  CHAR *tmp;
  CHAR msg[2048];
  CHAR *tmp1;

  tmp1 = argv[*index+1];

  if ( argv[*index+1] != NULL )
  {
    *data = strtol(argv[*index+1], &tmp  , 10);
    if (*data==0 && tmp1[0] != '0')
    {
      sprintf(msg, "Expect a int after option %s (got %s)\n ",
          argv[*index],
          argv[*index+1]);
      fprintf(stderr, msg);
      exit( 1 );
    }
  }
  else
  {
    sprintf(msg, "Expect a int after option %s (got %s)\n ",
        argv[*index],
        argv[*index+1]);
    fprintf(stderr, msg);
    exit( 1 );
  }
  *index =*index + 1;
}

void BankEfficiencyParseGetString(
  CHAR    **argv,
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


void BankEfficiencyParseGetDouble(
  CHAR    **argv,
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


void BankEfficiencyParseGetDouble2(
  CHAR    **argv,
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
void BankEfficiencyUpdateParams(
  InspiralCoarseBankIn   *coarseBankIn,
  RandomInspiralSignalIn *randIn,
  UserParametersIn       *userParam)
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

  if (randIn->param.fCutoff >=coarseBankIn->tSampling/2. -1.){
    randIn->param.fCutoff = coarseBankIn->tSampling/2 -1.;
  }


  /* BCV related */
  BankEfficiencyValidity(coarseBankIn->alpha, 0, 1, "--bank-alpha");
  BankEfficiencyCompare(coarseBankIn->psi0Min, coarseBankIn->psi0Max, "--bank-psi0-range");
  BankEfficiencyCompare(coarseBankIn->psi3Min, coarseBankIn->psi3Max, "--bank-psi3-range");
  BankEfficiencyValidity(coarseBankIn->psi0Min, 0, 1e13, "--bank-psi0-range min");
  BankEfficiencyValidity(coarseBankIn->psi0Max, 0, 1e13, "--bank-psi0-range max");
  BankEfficiencyValidity(coarseBankIn->psi3Min, -1e13, 100, "--bank-psi3-range min");
  BankEfficiencyValidity(coarseBankIn->psi3Max, -1e13, 100, "--bank-psi3-range max");
  BankEfficiencyCompare(coarseBankIn->LowGM, coarseBankIn->HighGM, "--bank-fcut-range");

  /* --- other --- */
  BankEfficiencyCompare( coarseBankIn->fUpper, coarseBankIn->tSampling/2,
      "--bank-ffinal and --sampling don't agree");

  /* --- eccentricity related --- */
  BankEfficiencyCompare( userParam->eccentricSignal.min,
      userParam->eccentricSignal.max, "--signal-eccentricity-range");
  BankEfficiencyValidity( userParam->eccentricSignal.min, 0,1,
      "--signal-eccentricity-range min");
  BankEfficiencyValidity( userParam->eccentricSignal.max, 0,1,
      "--signal-eccentricity-range max");
  BankEfficiencyCompare( userParam->eccentricBank.min,
      userParam->eccentricSignal.max, "--bank-eccentricity-range");
  BankEfficiencyValidity( userParam->eccentricBank.min, 0,1,
      "--bank-eccentricity-range min");
  BankEfficiencyValidity( userParam->eccentricBank.max, 0,1,
      "--bank-eccentricity-range max");


  BankEfficiencyCompare(coarseBankIn->mMin, coarseBankIn->mMax, "--bank-mass-range");
  BankEfficiencyValidity(coarseBankIn->mMin, 0, 1e6, "--bank-mass-range min");
  BankEfficiencyValidity(coarseBankIn->mMax, 0, 1e6, "--bank-mass-range max");

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

  /* --- fLower related --- */
  BankEfficiencyCompare( coarseBankIn->fLower, coarseBankIn->fUpper,
       "--bank-ffinal and --fl don't agree");
  BankEfficiencyValidity(coarseBankIn->fLower, 10, 1e6, "--fl --template-fl");
  BankEfficiencyValidity(randIn->param.fLower, 10, 1e6, "--fl --signal-fl");

  if (userParam->maxTotalMass != -1)
  {
      if (userParam->maxTotalMass < 2*randIn->mMin)
    {
      sprintf(msg,
          "--signal-max-total-mass (%f) must be > twice the minimal mass (%f) ",
          userParam->maxTotalMass, randIn->mMin );
          fprintf(stderr, msg);
          exit(1);
    }
  }

  if (userParam->minTotalMass != -1)
  {
    if (userParam->minTotalMass > 2*randIn->mMax)
    {
      sprintf(msg,
          "--min-total-mass (%f) must be < twice the maximal mass (%f) ",
          userParam->maxTotalMass , randIn->mMax );
          fprintf(stderr, msg);
          exit(1);
      }
    }



  if ((userParam->m1 != -1) && (userParam->m2 != -1))
  {
    randIn->param.mass1 = userParam->m1;
    randIn->param.mass2 = userParam->m2;
    /* we need to update the mMin which is used to estimate the length of
     * the vectors*/
    if (userParam->m1  > userParam->m2)
    {
      randIn->mMin = userParam->m2;
      randIn->mMax = userParam->m1+1e-2;
    }
    else
    {
      randIn->mMin = userParam->m1;
      randIn->mMax = userParam->m2+1e-2;
    }
    if (userParam->psi0 != -1 ||userParam->psi3 != -1
        || userParam->tau0 != -1 || userParam->tau3 != -1)
    {
      sprintf(msg,
      "--m1 --m2 --psi0 --psi3 --tau0 --tau3 error. "
      "If particular injection is requested,  you must choose either"
      " (--m1,--m2) options or (--psi0,--psi3) or (--tau0,--tau3)\n");
      fprintf(stderr, msg);
      exit(1);
    }
  }

  if ( (userParam->psi0 != -1 && userParam->psi3 != -1))
  {
   if (userParam->m1 != -1 || userParam->m2 != -1
       || userParam->tau0 != -1 || userParam->tau3 != -1)
   {
     sprintf(msg, "--m1 --m2 --psi0 --psi3 --tau0 --tau3 error. "
     "If particular injection is requested,  you must choose either "
     "(--m1,--m2) options or (--psi0,--psi3) or (--tau0,--tau3)\n");
     fprintf(stderr, msg);
     exit(1);
   }
  }

 if (userParam->tau0 != -1 && userParam->tau3 != -1)
 {
 	randIn->t0Min = userParam->tau0/1.1;
 	randIn->t0Max = userParam->tau0 *1.1;
   if (userParam->psi0 != -1 ||userParam->psi3 != -1
       || userParam->m1 != -1 || userParam->m2 != -1)
   {
     sprintf(msg,
         "--m1 --m2 --psi0 --psi3 --tau0 --tau3 error."
         " If particular injection is requested, you must choose"
         " either (--m1,--m2) or (--psi0,--psi3) or (--tau0,--tau3)\n");
     fprintf(stderr, msg);
     exit(1);
   }
  }


  BankEfficiencyValidity(coarseBankIn->mmCoarse, 0, 1, "--mm");
  BankEfficiencyValidity(coarseBankIn->numFcutTemplates, 0, 100, "--bank-number-fcut");
  BankEfficiencyCompare(randIn->mMin, randIn->mMax, "--signal-mass-range");


 /* -- MMAX and MMIn must be redefined given the input parameters*/
 randIn->MMax = randIn->mMax * 2;
 randIn->MMin = randIn->mMin * 2;

 randIn->etaMin = randIn->mMin * (randIn->MMax - randIn->mMin) /
    randIn->MMax / randIn->MMax;

 /* -- However, if MMAx and MMIN are also defined by the user*/

 if (userParam->maxTotalMass != -1 && userParam->maxTotalMass <2*randIn->mMax)
  {
    randIn->MMax = userParam->maxTotalMass;

    if (userParam->minTotalMass != -1
        && userParam->minTotalMass >2*randIn->mMin)
    {
       randIn->MMin = userParam->minTotalMass;
    }
  }

  BankEfficiencyCompare(randIn->t0Min, randIn->t0Max, "--signal-tau0-range");
  BankEfficiencyCompare(randIn->tnMin, randIn->tnMax, "--signal-tau3-range");
  BankEfficiencyValidity(randIn->t0Min, 0 ,1e6, "--signal-tau0-range min");
  BankEfficiencyValidity(randIn->t0Max, 0 ,1e6, "--signal-tau0-range max");
  BankEfficiencyValidity(randIn->tnMin, 0 ,1e6, "--signal-tau3-range min");
  BankEfficiencyValidity(randIn->tnMax, 0 ,1e6, "--signal-tau3-range max");


  if (userParam->faithfulness ==1 && randIn->type == 1)
    {
      fprintf(stderr, "No injection performed so the check option can not be used (change simulation-type option)\n");
      exit ( 1 );
    }
  if (userParam->binaryInjection == BHNS){
    if (randIn->mMin >3 || randIn->mMax < 3 ){
      fprintf(stderr,
      "BH-NS injection must have a mass less than 3 and the other greater than 3.\n");
      exit( 1 );
    }
  }

  if (userParam->signalfFinal>0)
  {
    randIn->param.fCutoff = userParam->signalfFinal;
    BankEfficiencyCompare(userParam->signalfFinal, coarseBankIn->tSampling/2.,
       "--signal-ffinal must be less than sampling by 2");
    BankEfficiencyValidity(userParam->signalfFinal, randIn->param.fLower ,1e6, "--signal-ffinal");
  }

  if (userParam->signalfFinal==0){
    userParam->signalfFinal = randIn->param.fCutoff ;
  }
}




/* ****************************************************************************
 *  Documenation  on line
 *  **************************************************************************/
void BankEfficiencyHelp(void)
{
  CHAR msg[2048];

  fprintf(stderr, "[NAME %s]\n ", CVS_NAME_STRING);
  fprintf(stderr, "[VERSION %s]\n ", CVS_ID_STRING);
  fprintf(stderr, "[VERSION %s]\n ", CVS_ID_STRING_C);
  fprintf(stderr, "[DESCRIPTION]\n");
  sprintf(msg, ""
      "\t For a detailled documentation, please see BankEfficiency.tex.\n");
  fprintf(stderr, "%s",msg);

  sprintf(msg, "[EXAMPLE]\n \t ./lalapps_BankEfficiency --template EOB --signal EOB --check --print-xml\n\n\n");
  fprintf(stderr, "%s",msg);

  fprintf(stderr, "[SYNOPSIS]\n");
  fprintf(stderr, "\t[--help]\n");
  fprintf(stderr, "\t[--verbose]  \t\t\t extra information on the screen \n" );
  fprintf(stderr, "\t[--ascii2xml]\t\t\t read the file Trigger.dat, concatenation of one or several outputs\n");
  fprintf(stderr, "\t             \t\t\t BankEfficiency and creates a unique XML outputs.\n");
  fprintf(stderr, "\t[--bank-alpha<float>]\t\t set the BCV alpha value in the moments computation\n");
  fprintf(stderr, "\t[--bank-eccentricity-range]\t set the range of eccentricity of the eccentric bank\n");
  fprintf(stderr, "\t[--bank-eccentricity-bins]\t set the number of layers for eccentric bank\n");
  fprintf(stderr, "\t[--bank-fcut-range<float float>] set the range of BCV fcut (in units of GM) \n");
  fprintf(stderr, "\t[--bank-ffinal<float>]\t\t set the final frequency to be used in the BCV moments computation\n");
  fprintf(stderr, "\t[--bank-grid-spacing <gridSpacing>]\t set the grid type of the BCV bank (Square, SquareNotOriented,\n\t\t HybridHexagonal,Hexagonal, HexagonalNotOriented\t\n");
  fprintf(stderr, "\t[--bank-number-fcut<integer>]\t number of BCV fcut layers\n");
  fprintf(stderr, "\t[--bank-mass-range<float float>] mass range to be covered by the SPA bank\n");
  fprintf(stderr, "\t[--bank-max-total-mass<float>]  max total mass of the bank\n");
  fprintf(stderr, "\t[--bank-min-total-mass<float>]  min total mass of the bank\n");
  fprintf(stderr, "\t[--bank-psi0-range<float float>] psi_0 range to be covered by the BCV bank\n");
  fprintf(stderr, "\t[--bank-psi3-range<float float>] psi_3 to be covered by the BCV bank\n");
  fprintf(stderr, "\t[--bank-from-file\t\t <ascii filename containing a list of masses m1 m2>\n");
  fprintf(stderr, "\t[--debug<integer>]\t\t set the debug level (same as in lal)\n");
  fprintf(stderr, "\t[--e-match<float>]\t\t set the e-match for the fast simulation\n");
  fprintf(stderr, "\t[--signal-fl<float>]\t\t set the lower cut off frequency of signal to inject\n");
  fprintf(stderr, "\t[--template-fl<float>]\t\t set the lower cut off frequnecy of template \n");
  fprintf(stderr, "\t[--fl<float>]\t\t\t set both template and signal lower cutoff frequency \n");
  fprintf(stderr, "\t[--m1<float>]\t\t\t force injection first individual mass to be equal to m1. needs to set m2 as well then\n");
  fprintf(stderr, "\t[--m2<float>]\t\t\t force injection second individual mass to be equal to m2. needs to set m1 as well then\n");
  fprintf(stderr, "\t[--mm<float>]\t\t\t set minimal match of the bank\n");
  fprintf(stderr, "\t[--n<float>]\t\t\t set number of trial in the simulation\n");
  fprintf(stderr, "\t[--nStartpad<int>]\t\t number of leading zero in injection\n");
  fprintf(stderr, "\t[--nEndpad<int>]\t\t minimum number of trailing zeroes in injection (currently used for turning tapering on/off)\n");
  fprintf(stderr, "\t[--ntrial<float>]\t\t same as --n\n");
  fprintf(stderr, "\t[--noise-amplitude<float>]\t set noise amplitude when using NoiseAndSignal flag simulation\n");
  fprintf(stderr, "\t[--noise-model<string>]\t\t set noise model curve to be <LIGOI, LIGOA, VIRGO, GEO, TAMA, REALPSD>\n");
  fprintf(stderr, "\t[--num-seconds<integer>]\t set number of seconds of data to look at.\n");
  fprintf(stderr, "\t[--psi0<float>]\t\t\t force injection psi0  value; request to psi3 as well. \n");
  fprintf(stderr, "\t[--psi3<float>]\t\t\t force injection psi3 value; request psi0 as well\n");
  fprintf(stderr, "\t[--sampling<float>]\t\t set sampling frequency (4096 by default).\n");
  fprintf(stderr, "\t[--seed<integer>]\t\t set seed for random generator.\n");
  fprintf(stderr, "\t[--signal<string>]\t\t set signal approximant (TaylorT1, TaylorT2, TaylorT3, TaylorF2,\n\t PadeT1, EOB, TaylorEt,  TaylorT4, TaylorN, SpinTaylor)\n");
  fprintf(stderr, "\t[--signal-alpha<float>]\t\t set alpha parameter of BCV injection\n");
  fprintf(stderr, "\t[--signal-amplitude<float>]\t set SNR of injection in the case NoiseandSignal simulation\n");
  fprintf(stderr, "\t[--signal-eccentricity-range]\t\t set the range of eccentricity of the injection\n");
  fprintf(stderr, "\t[--signal-max-total-mass<float>]\t set maximum total mass to be injected !! ");
  fprintf(stderr, "\t               \t does not work with spin injection\n");
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
  fprintf(stderr, "\t[--t0-fine-range<float>]\t\t\t force the tau0 of the templates/bank to be in this range \n");
  fprintf(stderr, "\t[--t3-fine-range<float>]\t\t\t force the tau3 of the templates/bank to be in this range \n");
  fprintf(stderr, "\t[--t0-fine-bin<float>]\t\t\t set number of templates in t0 direction \n");
  fprintf(stderr, "\t[--t3-fine-bin<float>]\t\t\t set number of templates in t3 direction\n");

  fprintf(stderr, "\t[--template<string>]\t\t set signal approximant (TaylorT1, TaylorT2, TaylorT3, TaylorF2\n\t, PadeT1, EOB, TaylorEt, TaylorT4, TaylorN, SpinTaylor)\n");
  fprintf(stderr, "\t[--template-order<integer>]\t set PN order of template\n");
  fprintf(stderr, "\t[--alpha-constraint]\t\t set BCV code to be constrained \n");
  fprintf(stderr, "\t[--bhns-injection]\t\t set injection to be only bhbs systems\n");
  fprintf(stderr, "\t[--no-alpha-constraint]\t\t set BCV code to be unconstrained\n");
  fprintf(stderr, "\t[--faithfulness]\t\t Check the code. Template parameters are equal to injection  parameters\n");
  fprintf(stderr, "\t                \t\t size of the bank is therefore unity. It computed the faithfulness instead of effectualness\n");
  fprintf(stderr, "\t[--real-noise]\t\t\t use real data and real PSD.force simulaion type to be Noise Only\n");
  fprintf(stderr, "\t[--no-start-phase]\t\t\t unset random phase which is always set to zero.\n");
  fprintf(stderr, "\t[--print-psd]\t\t\t print the psd in  a file BankEfficiency_PSD_type_gpstime.dat\n");
  fprintf(stderr, "\t[--print-best-overlap]\t\t print best overlap and other information\n");
  fprintf(stderr, "\t[--print-snr-histo]\t\t print histogram of the correlation output\n");
  fprintf(stderr, "\t[--print-bank]\t\t\t print the bank in ascii and xml format\n");
  fprintf(stderr, "\t[--print-xml]\t\t\t print the results in an xml file\n");
  fprintf(stderr, "\t[--print-prototype]\t\t print a prototype to be used by condor script\n");
  fprintf(stderr, "\t[--fast-simulation]\t\t perform fast simulation [None, EMatch,Heuristic1]\n");
  fprintf(stderr, "\t[--fast-param1]\t\t set maximum number of matches to compute without finding a greater match (Heuristic1 method)\n");
  exit(1);
}


void BankEfficiencyAscii2Xml(void)
{
  UINT4 countline = 0, nfast_max=0;

  UINT8  id = 0;

  ResultIn trigger;

  REAL4 tau0, tau3, tau0I, tau3I, psi0, psi3,phaseI, ecc, eccI,eccI_fl;
  REAL4 polarisation,inclination;
  REAL4 bestEMatch;

  FILE *input1;
  FILE *input2;
  FILE *bank;
  FILE *output;


  SnglInspiralTable     *inputData = NULL;

  INT4 numFileTriggers = 0;
  INT4 nStartPad;

  char sbuf[2048];
  CHAR fname[256]="";

  /* Main program starts here */
  /* First we open the file containing the ascii results */
  fprintf(stderr,"opening the xml output data file -- %s",
      BANKEFFICIENCY_XMLRESULTS);
  output = fopen(BANKEFFICIENCY_XMLRESULTS,"w");
  fprintf(stderr,"done\n");


  fprintf(stderr,"opening the input ascii file -- %s",
      BANKEFFICIENCY_ASCIIRESULTS);
  if  ( (input1  = fopen(BANKEFFICIENCY_ASCIIRESULTS,"r"))==NULL) {
    fprintf(stderr,"error while opening input file %s\n",
        BANKEFFICIENCY_ASCIIRESULTS);
    exit( 1 );
  }
  fprintf(stderr,"done\n");

  fprintf(stderr,"opening the xml prototype (argument of BankEfficiency code)");
  sprintf(fname,BANKEFFICIENCY_PROTOTYPE);
  strcat(fname,".xml");


  if  ( (input2  = fopen(fname,"r"))==NULL)
    {
      fprintf(stderr,"error while opening input file %s\n", fname);
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
  if  ( (bank  = fopen(BANKEFFICIENCY_XMLBANK,"r"))==NULL)
    {
      fprintf(stderr,"Error while opening input file %s\n",
          BANKEFFICIENCY_XMLBANK);
      fprintf(stderr,"The XML file will not contain the bank table\n");
    }
  else
    {
      /* read prototype and save in outputfile */
      fprintf(stderr,"parsing the bank  -- ");
      fprintf( stderr, "reading triggers from file: %s\n",
          BANKEFFICIENCY_XMLBANK);
      numFileTriggers = LALSnglInspiralTableFromLIGOLw( &inputData,
          BANKEFFICIENCY_XMLBANK , 0, -1 );


      fprintf(stderr," done %d\n", numFileTriggers);
      myfprintf(output, LIGOLW_XML_SNGL_INSPIRAL );
      while(inputData)
      {
      /*      id = inputData->event_id->id;*/

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
          (double)inputData->cont_chisq_dof,
          (int)inputData->cont_chisq,
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
        &trigger.mass1_trigger, &trigger.mass2_trigger,
    &psi0, &psi3,  &tau0, &tau3, &tau0I, &tau3I, &ecc,&eccI,&eccI_fl,
    &trigger.fend_trigger, &trigger.fend_inject,
    &trigger.mass1_inject, &trigger.mass2_inject,&inclination,&polarisation,
    &phaseI, &trigger.rho_final, &trigger.snrAtCoaTime, &trigger.phase,
    &trigger.alphaF, &trigger.bin, &nStartPad, &trigger.nfast, &nfast_max,
    &bestEMatch);


    fprintf(output, BANKEFFICIENCY_PARAMS_ROW,
        trigger.mass1_trigger, trigger.mass2_trigger,
        psi0, psi3, tau0, tau3, tau0I, tau3I,ecc,eccI,eccI_fl,
        trigger.fend_trigger, trigger.fend_inject,
        trigger.mass1_inject, trigger.mass2_inject,inclination, polarisation,
        phaseI, trigger.rho_final, trigger.snrAtCoaTime, trigger.phase,
        trigger.alphaF, trigger.bin, nStartPad, trigger.nfast, nfast_max,
        bestEMatch);
    fprintf(output,",\n");

    countline++;
  }


  fprintf(stderr,"read %d lines...done\n", countline);
  PRINT_LIGOLW_XML_TABLE_FOOTER(output);
  PRINT_LIGOLW_XML_FOOTER(output);

  fclose(output);
  fprintf(stderr,"closing xml file\n");

  exit(1);
}

/*
 * These two next functions are hack version of InspiralBankGeneration and
 * FineBank so that we can call a hybrid version of fineBank as if it was
 * available in BankGenerarion. Therefore, hardly any change is needed both in
 * lal and BankEfficiency. Once the code with this fine option is fully
 * tested, one can incorporate the change within lal.
 * */

void BankEfficiencyInspiralBankGeneration(
  LALStatus *status,
  InspiralCoarseBankIn *input,
  SnglInspiralTable **first,
  INT4 *ntiles,
  UserParametersIn userParam)
{
  InspiralTemplateList *coarseList = NULL;
  InspiralTemplateList *fineList = NULL;
  InspiralFineBankIn *fineIn = NULL;
  SnglInspiralTable *bank;
  INT4 cnt = 0, flist =0;
  InspiralMomentsEtc moments;

  INITSTATUS(status, "BankEfficiencyInspiralBankGeneration", BANKEFFICIENCYC);
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
  case EOBNR:
  case PadeT1:
  case PadeF1:
  case TaylorF1:
  case TaylorF2:
  case TaylorT1:
  case TaylorT2:
  case TaylorT3:
  case AmpCorPPN:
  case Eccentricity:

    /* Use LALInspiralCreateCoarseBank(). to init some variables */
    TRY( LALInspiralCreateCoarseBank( status->statusPtr, &coarseList, ntiles,
         *input ), status );

    /* --- then re-compute the moments, used later for the metric/gammas---*/
    LALGetInspiralMoments(status->statusPtr, &moments, &(input->shf), &(coarseList[0].params));

    fineIn = (InspiralFineBankIn *)LALMalloc(sizeof(InspiralFineBankIn));
    fineIn->coarseIn = *input;
    fineIn->templateList = coarseList[0];
    fineIn->templateList.params.t0 = 0.9384;
    fineIn->templateList.params.t3 = 0.235;



    TRY( BankEfficiencyInspiralCreateFineBank(status->statusPtr,
        &fineList, &flist, *fineIn, userParam), status);
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
      bank->mass1   = fineList[cnt].params.mass1;
      bank->mass2   = fineList[cnt].params.mass2;
      bank->mchirp  = fineList[cnt].params.chirpMass;
      bank->mtotal  = fineList[cnt].params.totalMass;
      bank->eta     = fineList[cnt].params.eta;
      bank->tau0    = fineList[cnt].params.t0;
      bank->tau2    = fineList[cnt].params.t2;
      bank->tau3    = fineList[cnt].params.t3;
      bank->tau4    = fineList[cnt].params.t4;
      bank->tau5    = fineList[cnt].params.t5;
      bank->ttotal  = fineList[cnt].params.tC;
      bank->psi0    = fineList[cnt].params.psi0;
      bank->psi3    = fineList[cnt].params.psi3;
      bank->f_final = fineList[cnt].params.fFinal;
      bank->eta     = fineList[cnt].params.eta;
      bank->beta    = fineList[cnt].params.beta;
      /* Copy the 10 metric co-efficients ... */
      /* the gamma are empty for the fine bank, so we need to recompute them*/
      LALInspiralComputeMetric(status->statusPtr,&(fineList[cnt].metric),
      &(fineList[cnt].params), &moments);
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



void BankEfficiencyInspiralCreateFineBank(
  LALStatus             *status,
  InspiralTemplateList **outlist,
  INT4                  *nlist,
  InspiralFineBankIn     fineIn,
  UserParametersIn       userParam)
{ /*
     </lalVerbatim> */

  REAL8 x0FineMin;
  REAL8 x1FineMin;
  REAL8 x0FineMax;
  REAL8 x1FineMax;
  INT4  i, j, validPars, bins0, bins1;
  InspiralTemplate   *tempPars=NULL;
  InspiralBankParams *bankPars=NULL;

  INITSTATUS (status, "BankEfficiencyInspiralCreateFineBank", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);
  ASSERT ((INT4)fineIn.coarseIn.space>=0,  status,
      LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
  ASSERT ((INT4)fineIn.coarseIn.space<=1,  status,
      LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
  ASSERT ((REAL8)fineIn.templateList.params.t0 > 0, status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);

  *nlist = 0;

  tempPars = (InspiralTemplate *) LALCalloc(1, sizeof(InspiralTemplate));
  bankPars = (InspiralBankParams *) LALCalloc(1, sizeof(InspiralBankParams));
  *tempPars = fineIn.templateList.params;
  switch (fineIn.coarseIn.space)
  {
    case Tau0Tau2:
      ASSERT (fineIn.templateList.params.t2 > 0, status,
          LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
      bankPars->x0 = fineIn.templateList.params.t0;
      bankPars->x1 = fineIn.templateList.params.t2;
      break;
    case Tau0Tau3:
      ASSERT (fineIn.templateList.params.t3 > 0, status,
          LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
      bankPars->x0 = fineIn.templateList.params.t0;
      bankPars->x1 = fineIn.templateList.params.t3;
      break;
    default: /* JC: DEFAULT CASE ADDED HERE */
      ABORT( status, 9999, "Default case in switch." );
  }

  x0FineMin = userParam.t0FineMin;
  x0FineMax = userParam.t0FineMax;
  x1FineMin = userParam.t3FineMin;
  x1FineMax = userParam.t3FineMax;
  bins0 = userParam.t0FineBin;
  bins1 = userParam.t3FineBin;
  bankPars->dx0 = (x0FineMax - x0FineMin)/((REAL4)bins0 -1);
  bankPars->dx1 = (x1FineMax - x1FineMin)/((REAL4)bins1 -1);


  bankPars->x1 = x1FineMin;
  for(i=0; i<=bins1; i++) {
    bankPars->x0 = x0FineMin;
    for(j=0; j<=bins0; j++) {
      LALInspiralValidTemplate(status->statusPtr,
          &validPars, *bankPars, fineIn.coarseIn);
      CHECKSTATUSPTR(status);

      if (validPars) {
        LALInspiralComputeParams(status->statusPtr,
            tempPars, *bankPars, fineIn.coarseIn);
        CHECKSTATUSPTR(status);

        if (!(*outlist = (InspiralTemplateList*)
              LALRealloc(*outlist, sizeof(InspiralTemplateList)*(*nlist+1))))
        {
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
  if  (*nlist==0){
    fprintf(stderr, "No valid templates found in the range mMin=%f and mMax=%f you may want to extedn --bank-mass-range\n",
      fineIn.coarseIn.mMin, fineIn.coarseIn.mMax);

  }
  if (tempPars!=NULL) LALFree(tempPars);
  if (bankPars!=NULL) LALFree(bankPars);

  DETATCHSTATUSPTR(status);
}

/* --- compute ematch between two points in a t0/t3 plane --- */
REAL4 BankEfficiencyComputeEMatch(
  RandomInspiralSignalIn *randIn,
  Mybank                   mybank,
  INT4                     index)
{
  REAL8 dt0;
  REAL8 dt3;
  REAL8 g00, g01, g11;
  REAL8 match = 0;

  BankEfficiencyValidity(index, 0,mybank.size,
      "Requested template index larger than banksize!");



  dt0 = -(randIn->param.t0 - mybank.tau0[index]);
  dt3 = -(randIn->param.t3 - mybank.tau3[index]);

  g00 = mybank.gamma3[index] - mybank.gamma1[index] *
      mybank.gamma1[index] / mybank.gamma0[index];
  g01 = mybank.gamma4[index] - mybank.gamma1[index] *
      mybank.gamma2[index] / mybank.gamma0[index];
  g11 = mybank.gamma5[index] - mybank.gamma2[index] *
      mybank.gamma2[index] / mybank.gamma0[index];
  /*fprintf(stderr, "dt0=%f dt3=%f g00=%f g01=%f g11=%f\n",dt0,dt3,g00,g11,g01);*/
  match = 1 - (g00*dt0*dt0 + 2*g01*dt0*dt3 + g11*dt3*dt3);

  return (REAL4)match;
}

/* ---   --- */
void BankEfficiencyUpdateSNRHistogram(
  REAL4Vector   *correlation,
  gsl_histogram *histogramNoise)
{
  INT4 i;


  /* fill histogram of the correlation output. Useful to get a flavour
   * of the output distribution in presence of noise only for
   * instance. */

  for (i=0; i<(INT4)correlation->length; i++)
  {
    /* in the unconstraint case, if alphafCut is applied,
     * correlation =-1 if alphaF > 1. Therefore we do not count it
     * in the correlation histogram*/
    if (correlation->data[i]!=-1)
    {
    gsl_histogram_increment(histogramNoise, correlation->data[i]);
    }
  }
}


/* --- function to create the template bank and related parameters --- */
void BankEfficiencyCreateTemplateBank(
  LALStatus              *status,
  InspiralCoarseBankIn   *coarseBankIn,
  MetadataTable          *templateBank,
  SnglInspiralTable      **tmpltHead,
  UserParametersIn        userParam,
  RandomInspiralSignalIn *randIn,
  INT4                   *sizeBank
  )
{

  INITSTATUS (status, "BankEfficiencyCreateTemplateBank", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);

    /* --- make sure the pointer to the first template is null --- */
  templateBank->snglInspiralTable = NULL;

  /* --- We create a fine template bank if requested --- */

  /* Note that the call will fail if the specified order is less than 2PN.
   * Yet, we want to create a bank, so let us force the order to be a 2PN
   * template bank whatever the template/signal order are.
   * */
  INT4 temp_order = coarseBankIn->order;
  /* ---  keep track of the order --- */
  if (coarseBankIn->order < 4){
    coarseBankIn->order = 4;
  }

  if (userParam.BankFromFile)
  {
    BankEfficiencyReadBankFromFile (status->statusPtr, tmpltHead,
      sizeBank, coarseBankIn, userParam);
  }
  else if (userParam.t0FineBin > 0 && userParam.t3FineBin > 0)
  {
    BankEfficiencyInspiralBankGeneration(status->statusPtr,
        coarseBankIn, tmpltHead, sizeBank, userParam);
    fprintf(stderr, "Number of templates=%d\n", *sizeBank);
  }
  else
  {
    /* call to the standard LAL template bank */
    LALInspiralBankGeneration(status->statusPtr,
        coarseBankIn, tmpltHead, sizeBank);
    fprintf(stderr, "Number of templates=%d\n", *sizeBank);
  }
  /* --- get back the order ---*/
  coarseBankIn->order = temp_order;  /* get back the order*/


  /* --- sanity check --- */
  if ( sizeBank ){
    templateBank->snglInspiralTable = *tmpltHead;
  }
  else
  {
    fprintf(stderr, "The template bank is empty. Quitting.");
    exit( 1 );
  }

  /* --- print the bank in an ascii file and/or XML format --- */
  if (userParam.printBank)
  {
    BankEfficiencyBankPrintXML(*templateBank, *coarseBankIn, *randIn, userParam);
    BankEfficiencyBankPrintAscii(*templateBank, *sizeBank, *coarseBankIn);
  }

  /*BankEfficiencyPrintMessage(sprintf("done. Got %d templates\n", *sizeBank)) ;*/

  DETATCHSTATUSPTR(status);
}


REAL4 GetMaxVector(REAL4 *vect, INT4 n)
{
  INT4 i;
  REAL4 max = -1e16;
  for (i=0; i<n; i++){
    if (vect[i] > max){
      max = vect[i];
    }
  }
  return max;
}

REAL4 GetMinVectorNonZero(REAL4 *vect, INT4 n)
{
  INT4 i;
  REAL4 min = 1e16;
  for (i=0; i<n; i++){
    if ((vect[i] < min) && (vect[i]!=0)){
      min = vect[i];
    }
  }
  return min;
}



void BankEfficiencyWaveOverlap(
  LALStatus              *status,
  REAL4Vector            *correlation,
  InspiralWaveOverlapIn  *overlapin,
  OverlapOutputIn        *overlapOutputThisTemplate,
  INT4                    startPad)
{
  InspiralWaveOverlapOut    overlapout,overlapout2;

  INITSTATUS (status, "BankEfficiencyWaveOverlap", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);

  /* --- just to be sure a value has been computed --- */
  overlapout.max = -1;
  overlapout2.max = -1;

  /* --- compute the output correlation between the signal (in overlapin)
   * and the template in overlapin. --- */
  overlapin->param.startPhase = 0.;
  LAL_CALL(LALInspiralWaveOverlap(status->statusPtr,
      correlation, &overlapout, overlapin), status->statusPtr);
/*BankEfficiencySaveVector("corr.dat",*correlation,2048);  */
/*      fprintf(stderr,"phase= %f and overlap(0)= %f and thisoverlap= %f e= %f\n",
        0., overlapout.max, overlapout2.max, overlapin->param.eccentricity);
*/
#if 0
  if (overlapin->param.approximant == Eccentricity)
  {
    INT4 N;
    REAL4 step;
    REAL4 phase;
    N=5;
    step = LAL_PI/(N-1.);
   	for (phase=step; phase<=LAL_PI+1e-6; phase+=step)
    {
      overlapin->param.startPhase = phase;
      LAL_CALL(LALInspiralWaveOverlap(status->statusPtr,
        correlation, &overlapout2, overlapin), status->statusPtr);

      fprintf(stderr,"phase= %f and overlap(0)= %f and thisoverlap= %f e= %f step=%f\n",
        phase, overlapout.max, overlapout2.max, overlapin->param.eccentricity,step);

      if (overlapout2.max>overlapout.max) {overlapout=overlapout2;}
    }
  }
#endif


  overlapOutputThisTemplate->rhoMax       = overlapout.max;
  overlapOutputThisTemplate->phase       = overlapout.phase;
  overlapOutputThisTemplate->rhoBin       = overlapout.bin;
  overlapOutputThisTemplate->freq         = overlapin->param.fFinal;
  overlapOutputThisTemplate->snrAtCoaTime = correlation->data[startPad];


  DETATCHSTATUSPTR(status);
}


void BankEfficiencySaveVectorAppend(
  const char  *filename,
  REAL4Vector  correlation,
  REAL4        tSampling)
{
  FILE *Foutput;
  INT4 i;

  Foutput=  fopen(filename, "a+");
  fprintf(Foutput, "&\n");
  for (i=0; i<(INT4)correlation.length; i++){
    fprintf(Foutput, "%e %e\n",
        (REAL4)i/tSampling,correlation.data[i]);
  }
  fclose(Foutput);
}

void BankEfficiencySaveVector(
  const char  *filename,
  REAL4Vector  correlation,
  REAL4        tSampling)
{
  FILE *Foutput;
  INT4 i;

  Foutput=  fopen(filename, "w");
  for (i=0; i<(INT4)correlation.length; i++){
    fprintf(Foutput, "%e %e\n",
        (REAL4)i/tSampling,correlation.data[i]);
  }
  fclose(Foutput);
}


void  BankEfficiencyFinalise(
  LALStatus              *status,
  Mybank                  bank,
  OverlapOutputIn         overlapOutputBestTemplate,
  RandomInspiralSignalIn  randIn,
  UserParametersIn        userParam,
  BankEfficiencySimulation simulation,
  InspiralCoarseBankIn    coarseBankIn)
{
  INT4               thisTemplateIndex = 0;
  InspiralTemplate   insptmplt;
  ResultIn           result;
  INT4 ntrials = simulation.ntrials;
  INT4 filter_processed = simulation.filter_processed;
  INITSTATUS (status, "BankEfficiencyWaveOverlap", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);

  thisTemplateIndex = overlapOutputBestTemplate.templateNumber;
  /* set the insptmplt structure */
  insptmplt = randIn.param;
  if (!userParam.faithfulness){
    LAL_CALL(BankEfficiencyCreateListFromTmplt(status->statusPtr, &insptmplt,
      bank, thisTemplateIndex), status->statusPtr);
  }

  if (userParam.template==BCV)
    insptmplt.massChoice=psi0Andpsi3;
  else
    insptmplt.massChoice=m1Andm2;

  LAL_CALL(LALInspiralParameterCalc( status->statusPtr,
      &(insptmplt) ), status->statusPtr);


  LAL_CALL(BankEfficiencyGetResult( status->statusPtr, &insptmplt,
     randIn.param, overlapOutputBestTemplate, &result, userParam),
      status->statusPtr);
  /* keep track of the number of templates actually used.*/
  result.nfast  = filter_processed;
  result.ntrial = ntrials;

  BankEfficiencyPrintResults(result, randIn, simulation);

  if (userParam.printResultXml){
    BankEfficiencyPrintResultsXml(coarseBankIn,randIn,userParam,result,simulation);
  }

  if (filter_processed == 0)
  {
    fprintf(stderr, "Warning : no filter processed. This could be related to\n");
    fprintf(stderr, "          the --fast-simulation option. Check that the \n");
    fprintf(stderr, "          template bank ranges are wide enough and cover\n");
    fprintf(stderr, "          the values of the signal parameters. \n");
    fprintf(stderr, "          Try to change --bank-mass-range or \n");
  }
    DETATCHSTATUSPTR(status);
}

/* --- alias function to populate the ambiguity function --- */
void BankEfficiencyPopulateAmbiguityFunction(
  gsl_matrix      *amb1,
  REAL4Vector      correlation,
  INT4             tmpltIndex,
  OverlapOutputIn  outputTemplate,
  INT4             startPad,
  InspiralTemplate insptmplt
)
{
  /* save t0/t3 coordinates of the current template insptmplt */
  gsl_matrix_set(amb1,0,tmpltIndex, insptmplt.t0);
  gsl_matrix_set(amb1,1,tmpltIndex, insptmplt.t3);
  /* --- save best SNR of this template (accumulate SNR over simulations)--- */
  gsl_matrix_set(amb1, 2, tmpltIndex, gsl_matrix_get(amb1, 2, tmpltIndex) +
      outputTemplate.rhoMax);
  /* --- save the SNR at t=0 (accumulate SNR over simulations) --- */
  gsl_matrix_set(amb1, 3, tmpltIndex, gsl_matrix_get(amb1,3,tmpltIndex) +
      correlation.data[startPad]);

}


/* --- create a mybank structure based on the original template
 * bank structure (SnglinspiralTable) --- */
void BankEfficiencyInitMyBank(
  Mybank            *mybank,
  INT4              *sizeBank,
  SnglInspiralTable *tmpltHead,
  UserParametersIn   userParam)
{
  SnglInspiralTable      *tmpltCurrent = NULL;
  INT4 i,j;
  INT4 eccentricBins;
  REAL4 eccentricMin,eccentricMax, eccentricStep;
  INT4 thisIndex;

  /* --- first get the eccentric parameters --- */
  eccentricBins = userParam.eccentricBank.bins;
  eccentricStep = userParam.eccentricBank.step;
  eccentricMax  = userParam.eccentricBank.max;
  eccentricMin  = userParam.eccentricBank.min;

  /* ---  the new bank size equals nbins times the original bank size --- */
  *sizeBank *= eccentricBins;

  /* Allocate memory for mybank vectors */
  mybank->mass1  = (REAL4 *) malloc(*sizeBank * sizeof(REAL4*));
  mybank->mass2  = (REAL4 *) malloc(*sizeBank * sizeof(REAL4*));
  mybank->tau0   = (REAL4 *) malloc(*sizeBank * sizeof(REAL4*));
  mybank->tau3   = (REAL4 *) malloc(*sizeBank * sizeof(REAL4*));
  mybank->psi0   = (REAL4 *) malloc(*sizeBank * sizeof(REAL4*));
  mybank->psi3   = (REAL4 *) malloc(*sizeBank * sizeof(REAL4*));
  mybank->alpha  = (REAL4 *) malloc(*sizeBank * sizeof(REAL4*));
  mybank->fFinal = (REAL4 *) malloc(*sizeBank * sizeof(REAL4*));
  mybank->index  = (UINT4 *) malloc(*sizeBank * sizeof(UINT4*));
  mybank->used   = (UINT4 *) malloc(*sizeBank * sizeof(UINT4*));
  mybank->eccentricity = (REAL4 *) malloc(*sizeBank * sizeof(REAL4*));
  mybank->snr    = (REAL4 *) malloc(*sizeBank * sizeof(REAL4*));
  mybank->gamma0 = (REAL4 *) malloc(*sizeBank * sizeof(REAL4*));
  mybank->gamma1 = (REAL4 *) malloc(*sizeBank * sizeof(REAL4*));
  mybank->gamma2 = (REAL4 *) malloc(*sizeBank * sizeof(REAL4*));
  mybank->gamma3 = (REAL4 *) malloc(*sizeBank * sizeof(REAL4*));
  mybank->gamma4 = (REAL4 *) malloc(*sizeBank * sizeof(REAL4*));
  mybank->gamma5 = (REAL4 *) malloc(*sizeBank * sizeof(REAL4*));


  /* --- populate the new bank structure ---*/
  for (tmpltCurrent = tmpltHead, i=0;
       tmpltCurrent ;
       tmpltCurrent = tmpltCurrent->next, i++)
  {

    for (j=0; j<(INT4)eccentricBins; j++)
    {
      thisIndex = i * eccentricBins + j;
      mybank->mass1[thisIndex]  = tmpltCurrent->mass1;
      mybank->mass2[thisIndex]  = tmpltCurrent->mass2;
      mybank->tau0[thisIndex]   = tmpltCurrent->tau0;
      mybank->tau3[thisIndex]   = tmpltCurrent->tau3;
      mybank->psi0[thisIndex]   = tmpltCurrent->psi0;
      mybank->psi3[thisIndex]   = tmpltCurrent->psi3;
      mybank->alpha[thisIndex]  = tmpltCurrent->alpha;
      mybank->fFinal[thisIndex] = tmpltCurrent->f_final;
      mybank->gamma0[thisIndex] = tmpltCurrent->Gamma[0];
      mybank->gamma1[thisIndex] = tmpltCurrent->Gamma[1];
      mybank->gamma2[thisIndex] = tmpltCurrent->Gamma[2];
      mybank->gamma3[thisIndex] = tmpltCurrent->Gamma[3];
      mybank->gamma4[thisIndex] = tmpltCurrent->Gamma[4];
      mybank->gamma5[thisIndex] = tmpltCurrent->Gamma[5];
      mybank->index[thisIndex]  = thisIndex;
      mybank->used[thisIndex]   = (INT4) 0;
      mybank->eccentricity[thisIndex] = eccentricMin + (REAL4)j*eccentricStep;
      mybank->snr[thisIndex]    = 0.;
    }
  }

  mybank->size = (INT4) (*sizeBank);
  mybank->approximant = userParam.template;
  if ( vrbflg )
  {
    fprintf(stderr, "done. Using %d layers of eccenticity.\n", eccentricBins);
    fprintf(stderr, "The new bank size is Got %d templates\n", *sizeBank );
    fflush(stderr);
  }

}


/* --- Set the eccentric values of the eccentric template bank--- */
void BankEfficiencyEccentricBankInit(
UserParametersIn *userParam)
{
  /* --- init eccentric bank parameters --- */
  /* if bins > 1, we compute the step */
  if (userParam->eccentricBank.bins > 1)
  {
    userParam->eccentricBank.step = (userParam->eccentricBank.max -
        userParam->eccentricBank.min) / (userParam->eccentricBank.bins-1);
  }
  else
  {
  	/* otherwise, only 1 bin is required, so let us return something large ---*/
  	userParam->eccentricBank.step = 2*(userParam->eccentricBank.max -
        userParam->eccentricBank.min);
  }

  if ( vrbflg && userParam->eccentricBank.bins>=1)
  {
    fprintf(stderr, "Eccentric bank requested.\n");
    fprintf(stderr, "--- minimum value = %f\n",userParam->eccentricBank.min);
    fprintf(stderr, "--- maximum value = %f\n",userParam->eccentricBank.max);
    fprintf(stderr, "--- number of bins = %d\n",userParam->eccentricBank.bins);
    fprintf(stderr, "--- step between layers = %f\n",userParam->eccentricBank.step);
  }
}

void BankEfficiencyPrintAmbiguity(
  UserParametersIn userParam,
  INT4             sizeBank,
  gsl_matrix       *amb1
)
{
  FILE *Foutput;
  CHAR str[512];
  INT4 i;

  if (userParam.ambiguity)
  {
    if (vrbflg)
    {
      fprintf(stderr,"------->%s\n",userParam.tag);
      sprintf(str, "BankEfficiency-ambiguity_%d_%s.dat",
          userParam.useed, userParam.tag);
    }
    Foutput=  fopen(str,"w");

    for (i=0; i<sizeBank; i++)
    {
      fprintf(Foutput, "%e %e %e %e\n",
        gsl_matrix_get(amb1, 0, i),
        gsl_matrix_get(amb1, 1, i),
        gsl_matrix_get(amb1, 2, i) / userParam.ntrials,
        gsl_matrix_get(amb1, 3, i) / userParam.ntrials
        );
    }
    fclose(Foutput);
  }
}

void BankEfficiencyError(const char *str)
{
  fprintf(stderr,"%s", str);
  exit(1);
}

void BankEfficiencyCompare(REAL4 a, REAL4 b, const char *str)
{
  CHAR str2[1024];
  sprintf(str2,"Error: the range provided (%s) is not sorted.\n",str);
  if (a > b ){
	BankEfficiencyError(str2);
  }
}

void BankEfficiencyValidity(
  REAL4 a,
  REAL4 min,
  REAL4 max,
  const char * str)
{
   CHAR str2[1024];
   sprintf(str2,"Error: the range provided %s is not correct.\n",str);
   if ( (a < min) || (a > max) ){
	  BankEfficiencyError(str2);
   }
 }

void BankEfficiencyPrintMessage(const char *str)
{
  if (vrbflg){
  fprintf(stderr, "%s",str);
  }
}

REAL4 BankEfficiencyRandom(REAL4 min, REAL4 max)
{
  REAL4 u;
  u = XLALUniformDeviate(randParams);
  return min + u * (max-min);
}


void GetClosestValidTemplate(
  Mybank bank,
  RandomInspiralSignalIn randIn,
  UINT4 *fast_index
  )
{
  REAL4 distance = -1e16;
  UINT4 index =  0;
  UINT4 *index_backup;
  REAL4 ematch=-1e16;

  if (bank.approximant == BCV)
  {
    UINT4 i = 0;
    index = 0;
    while (i < bank.size)
    {
      if (bank.used[i]==0)
      {
        bank.used[i] = 1;
        index = i;
        i = bank.size+1;
      }
      else
      {
        i++;
      }
    }
    *fast_index = index;
  }
  else
  {
    UINT4 i;
    REAL4 tau0,tau3;
    REAL4 eccentricity;

    GetTau03FromMasses(randIn, &tau0, &tau3);

    index_backup = (UINT4 *) malloc(bank.size * sizeof(UINT4 *));
    for (i=0; i<bank.size; i++){
      index_backup[i] = bank.used[i];
    }
    eccentricity = randIn.param.eccentricity;

    for (i=0; i<bank.size; i++)
    {
      ematch = BankEfficiencyComputeEMatch(&randIn,bank, i);
      if ((ematch > distance) && bank.used[i] == 0)
      {
        distance = ematch;
        index = i;
      }
    }

    for (i=0; i<bank.size; i++)
      bank.used[i] = index_backup[i];

    bank.used[index] = 1;
 *fast_index = index;
     free(index_backup);
  }

}


void GetTau03FromMasses(
  RandomInspiralSignalIn randIn,
  REAL4 *tau0,
  REAL4 *tau3)
{
  REAL4 totalMass, eta, piFl,mass1,mass2;
  REAL4 ratio=1;
  REAL4 eccentricity;

  eccentricity = randIn.param.eccentricity;
  mass1 = randIn.param.mass1;
  mass2 = randIn.param.mass2;
  eccentricity = randIn.param.eccentricity;
  /*(Empirical equations based on a full simulation with masses > 3) and fl=40*/

  if (eccentricity >=0 && eccentricity <=0.35)
  {
    ratio = (eccentricity -8.8)/-8.8;
  }
  else if (eccentricity >=0.35 && eccentricity <=0.55)
  {
    ratio = (eccentricity - 2.2)/-1.9;
  }
  else if (eccentricity >=0.55 && eccentricity <=1)
  {
    ratio = (eccentricity - 0.93)/-0.43;
  }
  else
  {
    ratio = 1;
  }

  totalMass = (mass1+mass2);
  eta = mass1*mass2/totalMass/totalMass;
  piFl = LAL_PI * randIn.param.fLower;

  totalMass *= LAL_MTSUN_SI;
  *tau0     = 5.0/(256.0 * eta * pow(totalMass,5./3.)*pow(piFl,8./3.));
  *tau3     = LAL_PI/(8.0 * eta * pow(totalMass,2./3.)*pow(piFl,5./3.));

  *tau0 *= ratio;
}

void BankEfficiencyReadBankFromFile (
  LALStatus            *status,
  SnglInspiralTable    **first,
  INT4                 *nlist,
  InspiralCoarseBankIn *coarseIn,
  UserParametersIn     userParam)
{

  UINT4 i=0, cnt;
  FILE *bankfile;
  static InspiralTemplateList *list;
  CHAR str[1024];
  REAL8 m1, m2;
  SnglInspiralTable *bank;
  InspiralMomentsEtc moments;

  INITSTATUS (status, "ReadBankFromFile", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);

  if ( (bankfile = fopen(userParam.BankFile, "r")) == NULL)
  {
    sprintf(str, "No file by name %s; aborting\n", userParam.BankFile);
    BankEfficiencyError(str);
  }

  while( fscanf(bankfile, "%le %le", &m1, &m2) != EOF)
  {
    if (m1 < coarseIn->mMax && m2 < coarseIn->mMax &&
        m1 > coarseIn->mMin && m2 > coarseIn->mMin)
    {
      if (i==0)
        list = (InspiralTemplateList *) LALMalloc(sizeof(InspiralTemplateList));
      else
        list = (InspiralTemplateList *) LALRealloc( list, sizeof(InspiralTemplateList) * (i+1) );
      if (!list)
      {
        ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
      }
      memset (list + i, 0, sizeof(InspiralTemplateList) );
      list[i].ID = i;
      list[i].params.mass1 = m1;
      list[i].params.mass2 = m2;
      list[i].params.fLower = coarseIn->fLower;
      list[i].params.fCutoff = coarseIn->fUpper;
      list[i].params.fFinal = coarseIn->fUpper;
      list[i].params.tSampling = coarseIn->tSampling;
      list[i].params.approximant = coarseIn->approximant;
      list[i].params.order = coarseIn->order;
      list[i].params.ampOrder = coarseIn->ampOrder;
      list[i].params.massChoice = m1Andm2;
      LALInspiralParameterCalc(status->statusPtr, &(list[i].params));
      /*
      fprintf(stderr, "%e %e %e %e %e %e\n", m1, m2, list[i].params.totalMass, list[i].params.eta, list[i].params.t0, list[i].params.t3);
      */
      i++;
    }
  }
  fclose(bankfile);
  *nlist = i;
  fprintf(stderr, "Number of templates=%d\n", *nlist);

  /* --- then re-compute the moments, used later for the metric/gammas---*/
  /*
  */
  LALGetInspiralMoments(status->statusPtr, &moments, &(coarseIn->shf), &(list[0].params));

  /* Convert output data structure. */


  bank = (SnglInspiralTable *) LALCalloc(1, sizeof(SnglInspiralTable));
  if (bank == NULL){
    ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
  }
  *first = bank;
  for( cnt = 0; cnt < i; cnt++ )
  {
     bank = bank->next = (SnglInspiralTable *) LALCalloc( 1, sizeof(
     SnglInspiralTable ) );
     if (bank == NULL)
     {
       ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
     }
     bank->mass1   = list[cnt].params.mass1;
     bank->mass2   = list[cnt].params.mass2;
     bank->mchirp  = list[cnt].params.chirpMass;
     bank->mtotal  = list[cnt].params.totalMass;
     bank->eta     = list[cnt].params.eta;
     bank->tau0    = list[cnt].params.t0;
     bank->tau2    = list[cnt].params.t2;
     bank->tau3    = list[cnt].params.t3;
     bank->tau4    = list[cnt].params.t4;
     bank->tau5    = list[cnt].params.t5;
     bank->ttotal  = list[cnt].params.tC;
     bank->psi0    = list[cnt].params.psi0;
     bank->psi3    = list[cnt].params.psi3;
     bank->f_final = list[cnt].params.fFinal;
     bank->eta     = list[cnt].params.eta;
     bank->beta    = list[cnt].params.beta;
     /* Copy the 10 metric co-efficients ... */
     /* the gamma are empty for the fine bank, so we need to recompute them*/
     LALInspiralComputeMetric(status->statusPtr,&(list[cnt].metric),
     &(list[cnt].params), &moments);
     memcpy (bank->Gamma, list[cnt].metric.Gamma, 10*sizeof(REAL4));
   }
    /* Free first template, which is blank. */
    bank = (*first)->next;
    LALFree( *first );
    *first = bank;
    LALFree(list);
    return;
  }
