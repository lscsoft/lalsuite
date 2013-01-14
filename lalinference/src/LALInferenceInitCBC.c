/*
 *  LALInferenceCBCInit.c:  Bayesian Followup initialisation routines.
 *
 *  Copyright (C) 2012 Vivien Raymond and John Veitch
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


#include <stdio.h>
#include <lal/Date.h>
#include <lal/GenerateInspiral.h>
#include <lal/LALInference.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/StringInput.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/TimeSeries.h>
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALInferenceProposal.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceReadData.h>
#include <lal/LALInferenceInit.h>

void LALInferenceInitCBCTemplate(LALInferenceRunState *runState)
{
  char help[]="(--template [LAL,PhenSpin,LALGenerateInspiral,LALSim]\tSpecify template (default LAL)\n";
  ProcessParamsTable *ppt=NULL;
  ProcessParamsTable *commandLine=runState->commandLine;
  /* Print command line arguments if help requested */
  //Help is taken care of in LALInferenceInitCBCVariables
  //ppt=LALInferenceGetProcParamVal(commandLine,"--help");
  //if(ppt)
  //{
  //	fprintf(stdout,"%s",help);
  //	return;
  //}
  /* This is the LAL template generator for inspiral signals */
  runState->template=&LALInferenceTemplateXLALSimInspiralChooseWaveform;
  ppt=LALInferenceGetProcParamVal(commandLine,"--template");
  if(ppt) {
    if(!strcmp("LALSTPN",ppt->value)){
      fprintf(stderr,"ERROR: --template LALSTPN is deprecated. Try LALGenerateInspiral instead...\n");
      exit(1);
    }
    else if(!strcmp("PhenSpin",ppt->value))
      runState->template=&LALInferenceTemplatePSTRD;
    else if(!strcmp("LALGenerateInspiral",ppt->value))
      runState->template=&LALInferenceTemplateLALGenerateInspiral;
    else if(!strcmp("SpinTaylor",ppt->value))
      runState->template=&LALInferenceTemplateLALGenerateInspiral;
    else if(!strcmp("LAL",ppt->value)){
      fprintf(stderr,"Warning: --template LAL is deprecated. Please use --template LALSim in future runs.\n");
      runState->template=&LALInferenceTemplateLAL;
    }
    else if(!strcmp("LALSim",ppt->value))
      runState->template=&LALInferenceTemplateXLALSimInspiralChooseWaveform;
    else {
      XLALPrintError("Error: unknown template %s\n",ppt->value);
      XLALPrintError(help);
      XLAL_ERROR_VOID(XLAL_EINVAL);
    }
  }
  else if(LALInferenceGetProcParamVal(commandLine,"--LALSimulation")){
    fprintf(stderr,"Warning: --LALSimulation is deprecated, the LALSimulation package is now the default. To use LALInspiral specify:\n\
                    --template LALGenerateInspiral (for time-domain templates)\n\
                    --template LAL (for frequency-domain templates)\n");
  }
  else {
    fprintf(stdout,"Template function called is \"LALInferenceTemplateXLALSimInspiralChooseWaveform\"\n");
  }
  return;
}


/* Setup the variables to control template generation */
/* Includes specification of prior ranges */

void LALInferenceInitCBCVariables(LALInferenceRunState *state)
{

  char help[]="\
               \n\
               ------------------------------------------------------------------------------------------------------------------\n\
               --- Injection Arguments ------------------------------------------------------------------------------------------\n\
               ------------------------------------------------------------------------------------------------------------------\n\
               (--inj injections.xml)          Injection XML file to use.\n\
               (--event N)                     Event number from Injection XML file to use.\n\
               \n\
               ------------------------------------------------------------------------------------------------------------------\n\
               --- Template Arguments -------------------------------------------------------------------------------------------\n\
               ------------------------------------------------------------------------------------------------------------------\n\
               (--symMassRatio)                Jump in symmetric mass ratio eta, instead of q=m2/m1.\n\
               (--template)                    Specify template [LAL,PhenSpin,LALGenerateInspiral,LALSim] (default LALSim).\n\
               (--approx)                      Specify a template approximant and phase order to use.\n\
                                               (default TaylorF2threePointFivePN). Available approximants:\n\
                                               default modeldomain=\"time\": GeneratePPN, TaylorT1, TaylorT2, TaylorT3, TaylorT4, \n\
                                                                           EOB, EOBNR, EOBNRv2, EOBNRv2HM, SpinTaylor, \n\
                                                                           SpinQuadTaylor, SpinTaylorFrameless, SpinTaylorT4, \n\
                                                                           PhenSpinTaylorRD, NumRel.\n\
                                               default modeldomain=\"frequency\": TaylorF1, TaylorF2, TaylorF2RedSpin, \n\
                                                                                TaylorF2RedSpinTidal, IMRPhenomA, IMRPhenomB.\n\
               (--amporder PNorder)            Specify a PN order in amplitude to use (defaults: LALSimulation: max available; LALInspiral: newtownian).\n\
               (--fref fRef)                   Specify a reference frequency at which parameters are defined (default 0).\n\
               (--tidal)                       Enables tidal corrections, only with LALSimulation.\n\
               (--interactionFlags)            intercation flags, only with LALSimuation (LAL_SIM_INSPIRAL_INTERACTION_ALL).\n\
               (--modeldomain)                 domain the waveform template will be computed in (\"time\" or \"frequency\").\n\
               (--spinAligned or --aligned-spin)  template will assume spins aligned with the orbital angular momentum.\n\
                                               *Enables* spins for TaylorF2, TaylorF2RedSpin, TaylorF2RedSpinTidal, IMRPhenomB.\n\
               (--singleSpin)                  template will assume only the spin of the most massive binary component exists.\n\
               (--noSpin, --disable-spin)      template will assume no spins.\n\
               \n\
               ------------------------------------------------------------------------------------------------------------------\n\
               --- Starting Parameters ------------------------------------------------------------------------------------------\n\
               ------------------------------------------------------------------------------------------------------------------\n\
               (--trigtime time)               Trigger time to use.\n\
               (--time time)                   Waveform time (overrides random about trigtime).\n\
               (--mc mchirp)                   Trigger chirpmass to use.\n\
               (--eta eta)                     Trigger eta (symmetric mass ratio) to use.\n\
               (--q q)                         Trigger q (asymmetric mass ratio) to use.\n\
               (--phi phase)                   Trigger phase to use.\n\
               (--iota inclination)            Trigger inclination to use.\n\
               (--dist dist)                   Trigger distance.\n\
               (--ra ra)                       Trigger RA.\n\
               (--dec dec)                     Trigger declination.\n\
               (--psi psi)                     Trigger psi.\n\
               (--a1 a1)                       Trigger a1.\n\
               (--theta1 theta1)               Trigger theta1.\n\
               (--phi1 phi1)                   Trigger phi1.\n\
               (--a2 a2)                       Trigger a2.\n\
               (--theta2 theta2)               Trigger theta2.\n\
               (--phi2 phi2)                   Trigger phi2.\n\
               (--lambda1)                     Trigger lambda1.\n\
               (--lambda2)                     Trigger lambda2.\n\
               \n\
               ------------------------------------------------------------------------------------------------------------------\n\
               --- Prior Arguments ----------------------------------------------------------------------------------------------\n\
               ------------------------------------------------------------------------------------------------------------------\n\
               (--mc-min mchirp)               Minimum chirp mass.\n\
               (--mc-max mchirp)               Maximum chirp mass.\n\
               (--eta-min etaMin)              Minimum eta.\n\
               (--eta-max etaMax)              Maximum eta.\n\
               (--q-min qMin)                  Minimum q.\n\
               (--q-max qMax)                  Maximum q.\n\
               (--comp-min min)                Minimum component mass (1.0).\n\
               (--comp-max max)                Maximum component mass (30.0).\n\
               (--mtotalmin min)               Minimum total mass (2.0).\n\
               (--mtotalmax max)               Maximum total mass (35.0).\n\
               (--iota-max max)                Maximum inclination angle (pi).\n\
               (--Dmin dist)                   Minimum distance in Mpc (1).\n\
               (--Dmax dist)                   Maximum distance in Mpc (100).\n\
               (--lambda1-min)                 Minimum lambda1 (0).\n\
               (--lambda1-max)                 Maximum lambda1 (3000).\n\
               (--lambda2-min)                 Minimum lambda2 (0).\n\
               (--lambda2-max)                 Maximum lambda2 (3000).\n\
               (--dt time)                     Width of time prior, centred around trigger (0.1s).\n\
               \n\
               ------------------------------------------------------------------------------------------------------------------\n\
               --- Fix Parameters -----------------------------------------------------------------------------------------------\n\
               ------------------------------------------------------------------------------------------------------------------\n\
               (--fixMc)                       Do not allow chirpmass to vary.\n\
               (--fixEta)                      Do not allow mass ratio to vary.\n\
               (--fixQ)                        Do not allow mass ratio to vary.\n\
               (--fixPhi)                      Do not allow phase to vary.\n\
               (--fixIota)                     Do not allow inclination to vary.\n\
               (--fixDist)                     Do not allow distance to vary.\n\
               (--fixRa)                       Do not allow RA to vary.\n\
               (--fixDec)                      Do not allow declination to vary.\n\
               (--fixPsi)                      Do not allow polarization to vary.\n\
               (--fixA1)                       Do not allow spin to vary.\n\
               (--fixTheta1)                   Do not allow spin 1 colatitude to vary.\n\
               (--fixPhi1)                     Do not allow spin 1 longitude to vary.\n\
               (--fixA2)                       Do not allow spin 2 to vary.\n\
               (--fixTheta2)                   Do not allow spin 2 colatitude to vary.\n\
               (--fixPhi2)                     Do not allow spin 2 longitude to vary.\n\
               (--fixTime)                     Do not allow coalescence time to vary.\n\
               (--fixLambda1)                  Do not allow lambda1 to vary.\n\
               (--fixLambda2)                  Do not allow lambda2 to vary.\n\
               (--pinparams)                   List of parameters to set to injected values [mchirp,asym_massratio,etc].\n";


  /* Print command line arguments if state was not allocated */
  if(state==NULL)
    {
      fprintf(stdout,"%s",help);
      return;
    }

  /* Print command line arguments if help requested */
  if(LALInferenceGetProcParamVal(state->commandLine,"--help"))
    {
      fprintf(stdout,"%s",help);
      return;
    }


  LALStatus status;
  memset(&status,0,sizeof(status));
  SimInspiralTable *injTable=NULL;
  LALInferenceVariables *priorArgs=state->priorArgs;
  state->currentParams=XLALCalloc(1,sizeof(LALInferenceVariables));
  LALInferenceVariables *currentParams=state->currentParams;
  ProcessParamsTable *commandLine=state->commandLine;
  ProcessParamsTable *ppt=NULL;
  ProcessParamsTable *ppt_order=NULL;
  //INT4 AmpOrder=0;
  LALPNOrder PhaseOrder=LAL_PNORDER_THREE_POINT_FIVE;
  LALPNOrder AmpOrder=-1;//LAL_PNORDER_THREE_POINT_FIVE;//LAL_PNORDER_NEWTONIAN;
  Approximant approx=TaylorF2;
  REAL8 fRef = 0.0;
  LALInferenceApplyTaper bookends = LALINFERENCE_TAPER_NONE;
  UINT4 analytic=0;
  LALInferenceIFOData *dataPtr;
  LALInferenceDomain modelDomain;
  UINT4 event=0;
  UINT4 i=0;
  REAL8 m1=0;
  REAL8 m2=0;
  REAL8 Dmin=1.0;
  REAL8 Dmax=100.0;
  REAL8 mcMin=1.0;
  REAL8 mcMax=15.3;
  REAL8 mMin=1.0,mMax=30.0;
  REAL8 MTotMax=35.0;
  REAL8 MTotMin=2.0;
  REAL8 etaMin=0.0312;
  REAL8 etaMax=0.25;
  REAL8 qMin=mMin/mMax;
  REAL8 qMax=1.0;
  REAL8 m1min=mMin,m1max=mMax;
  REAL8 m2min=mMin,m2max=mMax;
  REAL8 iotaMin=0.0,iotaMax=LAL_PI;
  REAL8 psiMin=0.0,psiMax=LAL_PI;
  REAL8 decMin=-LAL_PI/2.0,decMax=LAL_PI/2.0;
  REAL8 raMin=0.0,raMax=LAL_TWOPI;
  REAL8 phiMin=0.0,phiMax=LAL_TWOPI;
  REAL8 a1min=0.0,a1max=1.0;
  REAL8 a2min=0.0,a2max=1.0;
  REAL8 theta1min=0.0,theta1max=LAL_PI;
  REAL8 theta2min=0.0,theta2max=LAL_PI;
  REAL8 phi1min=0.0,phi1max=LAL_TWOPI;
  REAL8 phi2min=0.0,phi2max=LAL_TWOPI;
  REAL8 dt=0.1;            /* Width of time prior */
  REAL8 lambda1Min=0.0;
  REAL8 lambda1Max=3000.0;
  REAL8 lambda2Min=0.0;
  REAL8 lambda2Max=3000.0;  
  REAL8 tmpMin,tmpMax;//,tmpVal;
  gsl_rng * GSLrandom=state->GSLrandom;
  REAL8 endtime=0.0, timeParam=0.0;
  REAL8 timeMin=endtime-dt,timeMax=endtime+dt;
  //REAL8 start_mc			=4.82+gsl_ran_gaussian(GSLrandom,0.025);
  REAL8 start_mc			=mcMin+gsl_rng_uniform(GSLrandom)*(mcMax-mcMin);
  REAL8 start_eta			=etaMin+gsl_rng_uniform(GSLrandom)*(etaMax-etaMin);
  REAL8 start_q           =qMin+gsl_rng_uniform(GSLrandom)*(qMax-qMin);
  REAL8 start_phase		=0.0+gsl_rng_uniform(GSLrandom)*(LAL_TWOPI-0.0);
  //REAL8 start_dist		=8.07955+gsl_ran_gaussian(GSLrandom,1.1);
  REAL8 start_dist		=Dmin+gsl_rng_uniform(GSLrandom)*(Dmax-Dmin);
  REAL8 start_ra			=0.0+gsl_rng_uniform(GSLrandom)*(LAL_TWOPI-0.0);
  REAL8 start_dec			=-LAL_PI/2.0+gsl_rng_uniform(GSLrandom)*(LAL_PI_2-(-LAL_PI_2));
  REAL8 start_psi			=0.0+gsl_rng_uniform(GSLrandom)*(LAL_PI-0.0);
  REAL8 start_iota		=0.0+gsl_rng_uniform(GSLrandom)*(LAL_PI-0.0);
  REAL8 start_a_spin1		=0.0+gsl_rng_uniform(GSLrandom)*(1.0-0.0);
  REAL8 start_theta_spin1 =0.0+gsl_rng_uniform(GSLrandom)*(LAL_PI-0.0);
  REAL8 start_phi_spin1	=0.0+gsl_rng_uniform(GSLrandom)*(LAL_TWOPI-0.0);
  REAL8 start_a_spin2		=0.0+gsl_rng_uniform(GSLrandom)*(1.0-0.0);
  REAL8 start_theta_spin2 =0.0+gsl_rng_uniform(GSLrandom)*(LAL_PI-0.0);
  REAL8 start_phi_spin2	=0.0+gsl_rng_uniform(GSLrandom)*(LAL_TWOPI-0.0);
  REAL8 start_lambda1 =lambda1Min+gsl_rng_uniform(GSLrandom)*(lambda1Max-lambda1Min);
  REAL8 start_lambda2 =lambda2Min+gsl_rng_uniform(GSLrandom)*(lambda2Max-lambda2Min);
  UINT4 spinAligned=0;
  UINT4 singleSpin=0;
  UINT4 noSpin=0;
  memset(currentParams,0,sizeof(LALInferenceVariables));

  /* Over-ride prior bounds if analytic test */
  if (LALInferenceGetProcParamVal(commandLine, "--correlatedGaussianLikelihood")) {
    analytic  = 1;
    m1min     = 14.927715;
    m1max     = 17.072285;
    m2min     = 5.829675;
    m2max     = 8.170325;
    iotaMin   = 1.4054428267948966;
    iotaMax   = 1.7361498267948965;
    phiMin    = 2.8701521535897934;
    phiMax    = 3.413033153589793;
    psiMin    = 1.3885563267948966;
    psiMax    = 1.7530363267948965;
    raMin     = 2.813050153589793;
    raMax     = 3.4701351535897933;
    decMin    = -0.300699;
    decMax    = 0.300699;
    Dmin      = 37.986000000000004;
    Dmax      = 62.013999999999996;
    timeMin   = -0.1073625;
    timeMax   = 0.1073625;
    a1min     = 0.3784565;
    a1max     = 0.6215435;
    a2min     = 0.421869;
    a2max     = 0.578131;
    theta1min = 1.3993998267948966;
    theta1max = 1.7421928267948965;
    theta2min = 1.4086158267948965;
    theta2max = 1.7329768267948966;
    phi1min   = 2.781852653589793;
    phi1max   = 3.501332653589793;
    phi2min   = 2.777215653589793;
    phi2max   = 3.5059696535897933;
  } else if (LALInferenceGetProcParamVal(commandLine, "--bimodalGaussianLikelihood")) {
    analytic  = 1;
    m1min     = 14.927715;
    m1max     = 18.787941;
    m2min     = 5.829675;
    m2max     = 10.042845;
    iotaMin   = 0.6200446634;
    iotaMax   = 1.2153172634;
    phiMin    = 1.2993558268;
    phiMax    = 2.2765416268;
    psiMin    = 0.6031581634;
    psiMax    = 1.2592221634;
    raMin     = 1.2422538268;
    raMax     = 2.4250068268;
    decMin    = -1.0860971634;
    decMax    = -0.0035807634;
    Dmin      = 12.986;
    Dmax      = 56.2364;
    timeMin   = -0.1373625;
    timeMax   = 0.2491425;
    a1min     = 0.0784565;
    a1max     = 0.5160131;
    a2min     = 0.121869;
    a2max     = 0.4031406;
    theta1min = 0.6140016634;
    theta1max = 1.2310290634;
    theta2min = 0.6232176634;
    theta2max = 1.2070674634;
    phi1min   = 1.2110563268;
    phi1max   = 2.5061203268;
    phi2min   = 1.2064193268;
    phi2max   = 2.5181765268;
  } else if (LALInferenceGetProcParamVal(commandLine, "--rosenbrockLikelihood")) {
    analytic  = 1;
    m1min     = 14.0;
    m1max     = 18.0;
    m2min     = 5.0;
    m2max     = 9.0;
    iotaMin   = -0.429203673;
    iotaMax   = 3.570796327;
    phiMin    = 1.141592654;
    phiMax    = 5.141592654;
    psiMin    = -0.429203673;
    psiMax    = 3.570796327;
    raMin     = 1.141592654;
    raMax     = 5.141592654;
    decMin    = -2.0;
    decMax    = 2.0;
    Dmin      = 48.0;
    Dmax      = 52.0;
    timeMin   = -2.0;
    timeMax   = 2.0;
    a1min     = -1.5;
    a1max     = 2.5;
    a2min     = -1.5;
    a2max     = 2.5;
    theta1min = -0.429203673;
    theta1max = 3.570796327;
    theta2min = -0.429203673;
    theta2max = 3.570796327;
    phi1min   = 1.141592654;
    phi1max   = 5.141592654;
    phi2min   = 1.141592654;
    phi2max   = 5.141592654;
  }


  if(LALInferenceGetProcParamVal(commandLine,"--skyLocPrior")){
    MTotMax=20.0;
    mMin=1.0;
    mMax=15.0;
    qMin=mMin/mMax;
    Dmin=10.0;
    Dmax=40.0;
    REAL8 densityVNR=1000.0;
    LALInferenceAddVariable(state->priorArgs,"densityVNR", &densityVNR , LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  }

  /* Read injection XML file for parameters if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--inj");
  if(ppt){
    SimInspiralTableFromLIGOLw(&injTable,ppt->value,0,0);
    if(!injTable){
      fprintf(stderr,"Unable to open injection file %s\n",ppt->value);
      exit(1);
    }
    ppt=LALInferenceGetProcParamVal(commandLine,"--event");
    if(ppt){
      event= atoi(ppt->value);
      fprintf(stderr,"Reading event %d from file\n",event);
      i=0;
      while(i<event) {i++; injTable=injTable->next;} /* select event */
    }
  }

  /* See if there are any parameters pinned to injection values */
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--pinparams"))){
    char *pinned_params=ppt->value;
    LALInferenceVariables tempParams;
    memset(&tempParams,0,sizeof(tempParams));
    char **strings=NULL;
    UINT4 N;
    LALInferenceParseCharacterOptionString(pinned_params,&strings,&N);
    LALInferenceInjectionToVariables(injTable,&tempParams);
    LALInferenceVariableItem *node=NULL;
    while(N>0){
      N--;
      char *name=strings[N];
      node=LALInferenceGetItem(&tempParams,name);
      if(node) LALInferenceAddVariable(currentParams,node->name,node->value,node->type,node->vary);
      else {fprintf(stderr,"Error: Cannot pin parameter %s. No such parameter found in injection!\n",node->name);}
    }
  }
  
  
  /* Over-ride approximant if user specifies */
  ppt=LALInferenceGetProcParamVal(commandLine,"--approx");
  if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--approximant"); //FIXME
  if(ppt){
    approx = XLALGetApproximantFromString(ppt->value);
    if( (int) approx == XLAL_FAILURE)
      ABORTXLAL(&status);
    ppt_order=LALInferenceGetProcParamVal(commandLine,"--order");
    if(ppt_order) PhaseOrder = XLALGetOrderFromString(ppt_order->value);
    else PhaseOrder = XLALGetOrderFromString(ppt->value);
    if( (int) PhaseOrder == XLAL_FAILURE)
      ABORTXLAL(&status);
  }
  ppt=LALInferenceGetProcParamVal(commandLine,"--amporder");
  if(ppt) AmpOrder=atoi(ppt->value);
  ppt=LALInferenceGetProcParamVal(commandLine,"--ampOrder");
  if(ppt) AmpOrder = XLALGetOrderFromString(ppt->value);
  fprintf(stdout,"Templates will run using Approximant %i, phase order %i, amp order %i\n",approx,PhaseOrder,AmpOrder);
  
  /* Set the modeldomain appropriately */
  switch(approx)
  {
    case GeneratePPN:
    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
    case TaylorT4:
    case EOB:
    case EOBNR:
    case EOBNRv2:
    case EOBNRv2HM:
    case SpinTaylor:
    case SpinTaylorT4:
    case SpinQuadTaylor:
    case SpinTaylorFrameless:
    case PhenSpinTaylorRD:
    case NumRel:
      modelDomain=LALINFERENCE_DOMAIN_TIME;
      break;
    case TaylorF1:
    case TaylorF2:
    case TaylorF2RedSpin:
    case TaylorF2RedSpinTidal:
    case IMRPhenomA:
    case IMRPhenomB:
      modelDomain=LALINFERENCE_DOMAIN_FREQUENCY;
      break;
    default:
      fprintf(stderr,"ERROR. Unknown approximant number %i. Unable to choose time or frequency domain model.",approx);
      exit(1);
      break;
  }
  
  
  ppt=LALInferenceGetProcParamVal(commandLine, "--fref");
  if (ppt) fRef = atof(ppt->value);

  ppt=LALInferenceGetProcParamVal(commandLine,"--modeldomain");
  if(ppt){
    if ( ! strcmp( "time", ppt->value ) )
    {
      modelDomain = LALINFERENCE_DOMAIN_TIME;
    }
    else if ( ! strcmp( "frequency", ppt->value ) )
    {
      modelDomain = LALINFERENCE_DOMAIN_FREQUENCY;
    }
    else
    {
      fprintf( stderr, "invalid argument to --modeldomain:\n"
              "unknown domain specified: "
              "domain must be one of: time, frequency\n");
      exit( 1 );
    }
  }
  
  dataPtr = state->data;
  while (dataPtr != NULL) {
    dataPtr->modelDomain = modelDomain;
    dataPtr = dataPtr->next;
  }
  /* This flag was added to account for the broken Big Dog
     injection, which had the opposite sign in H and L compared
     to Virgo. */
  if (LALInferenceGetProcParamVal(commandLine, "--crazyInjectionHLSign")) {
    INT4 flag = 1;
    LALInferenceAddVariable(currentParams, "crazyInjectionHLSign", &flag, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
  } else {
    INT4 flag = 0;
    LALInferenceAddVariable(currentParams, "crazyInjectionHLSign", &flag, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
  }

  /* Over-ride taper if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--taper");
  if(ppt){
    if(strstr(ppt->value,"STARTEND")) bookends=LALINFERENCE_TAPER_STARTEND;
    if(strstr(ppt->value,"STARTONLY")) bookends=LALINFERENCE_TAPER_START;
    if(strstr(ppt->value,"ENDONLY")) bookends=LALINFERENCE_TAPER_END;
    if(strstr(ppt->value,"RING")) bookends=LALINFERENCE_RING;
    if(strstr(ppt->value,"SMOOTH")) bookends=LALINFERENCE_SMOOTH;
  }

  /* Over-ride end time if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--trigtime");
  if(ppt && !analytic){
    endtime=atof(ppt->value);
    timeMin=endtime-dt; timeMax=endtime+dt;
    printf("Read end time %f\n",endtime);

  }

  /* Over-ride chirp mass if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--mc");
  if(ppt){
    start_mc=atof(ppt->value);
  }

  /* Over-ride eta if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--eta");
  if(ppt){
    start_eta=atof(ppt->value);
    LALInferenceMcEta2Masses(start_mc, start_eta, &m1, &m2);
    start_q=m2/m1;
  }

  /* Over-ride q if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--q");
  if(ppt){
    start_q=atof(ppt->value);
    LALInferenceQ2Eta(start_q, &start_eta);
  }

  /* Over-ride phase if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--phi");
  if(ppt){
    start_phase=atof(ppt->value);
  }

  /* Over-ride inclination if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--iota");
  if(ppt){
    start_iota=atof(ppt->value);
  }

  /* Over-ride distance if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--dist");
  if (ppt) {
    start_dist = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--ra");
  if (ppt) {
    start_ra = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--dec");
  if (ppt) {
    start_dec = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--psi");
  if (ppt) {
    start_psi = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--a1");
  if (ppt) {
    start_a_spin1 = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--theta1");
  if (ppt) {
    start_theta_spin1 = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--phi1");
  if (ppt) {
    start_phi_spin1 = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--a2");
  if (ppt) {
    start_a_spin2 = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--theta2");
  if (ppt) {
    start_theta_spin2 = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--phi2");
  if (ppt) {
    start_phi_spin2 = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--lambda1");
  if (ppt) {
    start_lambda1 = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--lambda2");
  if (ppt) {
    start_lambda2 = atof(ppt->value);
  }
  
  /* Over-ride time prior if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--dt");
  if(ppt){
    dt=atof(ppt->value);
  }

  /* Over-ride Distance min if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--Dmin");
  if(ppt){
    Dmin=atof(ppt->value);
  }

  /* Over-ride Distance max if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--Dmax");
  if(ppt){
    Dmax=atof(ppt->value);
  }

  /* Over-ride lambda1 min if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--lambda1-min");
  if(ppt){
    lambda1Min=atof(ppt->value);
  }
  
  /* Over-ride lambda1 max if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--lambda1-max");
  if(ppt){
    lambda1Max=atof(ppt->value);
  }
  
  /* Over-ride lambda2 min if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--lambda2-min");
  if(ppt){
    lambda2Min=atof(ppt->value);
  }
  
  /* Over-ride lambda2 max if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--lambda2-max");
  if(ppt){
    lambda2Max=atof(ppt->value);
  }
  
  /* Over-ride component masses */
  ppt=LALInferenceGetProcParamVal(commandLine,"--comp-min");
  if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--compmin");
  if(ppt)	mMin=atof(ppt->value);
  LALInferenceAddVariable(priorArgs,"component_min",&mMin,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--comp-max");
  if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--compmax");
  if(ppt)	mMax=atof(ppt->value);
  LALInferenceAddVariable(priorArgs,"component_max",&mMax,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
  
  
  /* Over-ride Mass priors if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--mc-min");
  if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--Mmin");
  if(ppt)
    {
      mcMin=atof(ppt->value);
      if (mcMin < 0)
        {
          fprintf(stderr,"ERROR: Minimum value of mchirp must be > 0");
          exit(1);
        }
    }
  else mcMin=pow(mMin*mMin,0.6)/pow(2.0*mMin,0.2);

  ppt=LALInferenceGetProcParamVal(commandLine,"--mc-max");
  if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--Mmax");
  if(ppt)
    {
      mcMax=atof(ppt->value);
      if (mcMax <= 0)
        {
          fprintf(stderr,"ERROR: Maximum value of mchirp must be > 0");
          exit(1);
        }
    }
  else mcMax=pow(mMax*mMax,0.6)/pow(2.0*mMax,0.2);

  ppt=LALInferenceGetProcParamVal(commandLine,"--eta-min");
  if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--etamin");
  if(ppt)
    {
      etaMin=atof(ppt->value);
      if (etaMin < 0.0)
        {
          fprintf(stderr,"ERROR: Minimum value of eta must be > 0");
          exit(1);
        }
    }

  ppt=LALInferenceGetProcParamVal(commandLine,"--eta-max");
  if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--etamax");
  if(ppt)
    {
      etaMax=atof(ppt->value);
      if (etaMax > 0.25 || etaMax <= 0.0)
        {
          fprintf(stderr,"ERROR: Maximum value of eta must be between 0 and 0.25\n");
          exit(1);
        }
    }


  ppt=LALInferenceGetProcParamVal(commandLine,"--MTotMax");
  if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--mtotalmax");
  if(ppt)	MTotMax=atof(ppt->value);
  LALInferenceAddVariable(priorArgs,"MTotMax",&MTotMax,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--MTotMin");
  if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--mtotalmin");
  if(ppt)	MTotMin=atof(ppt->value);
  LALInferenceAddVariable(priorArgs,"MTotMin",&MTotMin,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
  
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--q-min");
  if(ppt)
    {
      qMin=atof(ppt->value);
      if (qMin <= 0.0 || qMin < mMin/mMax || qMin < mMin/(MTotMax-mMin) || qMin > 1.0)
        {
          fprintf(stderr,"ERROR: invalid qMin ( max{0,mMin/mMax,mMin/(MTotMax-mMin) < q < 1.0} )");
          exit(1);
        }
    }

  ppt=LALInferenceGetProcParamVal(commandLine,"--q-max");
  if(ppt)
    {
      qMax=atof(ppt->value);
      if (qMax > 1.0 || qMax <= 0.0 || qMax < mMin/mMax || qMax < mMin/(MTotMax-mMin))
        {
          fprintf(stderr,"ERROR: invalid qMax ( max{0,mMin/mMax,mMin/(MTotMax-mMin) < q < 1.0} )");
          exit(1);
        }
    }

  ppt=LALInferenceGetProcParamVal(commandLine,"--iota-max");
  if (ppt) {
    iotaMax = atof(ppt->value);
  }


  LALInferenceAddVariable(currentParams, "LAL_APPROXIMANT", &approx,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(currentParams, "LAL_PNORDER",     &PhaseOrder,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
  if(LALInferenceGetProcParamVal(commandLine,"--ampOrder")) 
    LALInferenceAddVariable(currentParams, "LAL_AMPORDER",     &AmpOrder,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);

  if (LALInferenceGetProcParamVal(commandLine, "--fref"))
    LALInferenceAddVariable(currentParams, "fRef", &fRef, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);

  ppt=LALInferenceGetProcParamVal(commandLine,"--taper");
  if(ppt){
    LALInferenceAddVariable(currentParams, "LALINFERENCE_TAPER",     &bookends,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
  }
  ppt=LALInferenceGetProcParamVal(commandLine,"--newswitch");
  int newswitch=0;
  if(ppt){
    newswitch=1;
    LALInferenceAddVariable(currentParams, "newswitch", &newswitch, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
  }
  /* Set up the variable parameters */
  /* Jump in component masses if anaytic test */
  if (LALInferenceGetProcParamVal(commandLine, "--correlatedGaussianLikelihood") || 
      LALInferenceGetProcParamVal(commandLine, "--bimodalGaussianLikelihood") ||
      LALInferenceGetProcParamVal(commandLine, "--rosenbrockLikelihood") ||
      LALInferenceGetProcParamVal(commandLine, "--analyticnullprior")) {
    LALInferenceMcEta2Masses(start_mc, start_eta, &m1, &m2);
    LALInferenceAddVariable(currentParams, "m1",    &m1,    LALINFERENCE_REAL8_t,	LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddMinMaxPrior(priorArgs,	"m1",	&m1min,	&m1max,		LALINFERENCE_REAL8_t);
    LALInferenceAddVariable(currentParams, "m2",    &m2,    LALINFERENCE_REAL8_t,	LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddMinMaxPrior(priorArgs,	"m2",	&m2min,	&m2max,		LALINFERENCE_REAL8_t);
  } else {

    ppt=LALInferenceGetProcParamVal(commandLine,"--fixMc");
    if(ppt){
      LALInferenceAddVariable(currentParams, "chirpmass",    &start_mc,    LALINFERENCE_REAL8_t,	LALINFERENCE_PARAM_FIXED);
      if(lalDebugLevel>0) fprintf(stdout,"chirpmass fixed and set to %f\n",start_mc);
    }else{
      LALInferenceAddVariable(currentParams, "chirpmass",    &start_mc,    LALINFERENCE_REAL8_t,	LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs,	"chirpmass",	&mcMin,	&mcMax,		LALINFERENCE_REAL8_t);

    /* Check if running with symmetric (eta) or asymmetric (q) mass ratio.*/
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixQ");
    if(ppt){
      LALInferenceAddVariable(currentParams, "asym_massratio", &start_q, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      LALInferenceAddMinMaxPrior(priorArgs,	"asym_massratio",	&qMin,	&qMax,	LALINFERENCE_REAL8_t);
      if(lalDebugLevel>0) fprintf(stdout,"q fixed and set to %f\n",start_q);
    }else{
      ppt=LALInferenceGetProcParamVal(commandLine,"--fixEta");
      if(ppt){
        LALInferenceAddVariable(currentParams, "massratio", &start_eta, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        LALInferenceAddMinMaxPrior(priorArgs,	"massratio", &etaMin, &etaMax, LALINFERENCE_REAL8_t);
        if(lalDebugLevel>0) fprintf(stdout,"eta fixed and set to %f\n",start_eta);
      }else{
        ppt=LALInferenceGetProcParamVal(commandLine,"--symMassRatio");
        if(ppt){
          LALInferenceAddVariable(currentParams, "massratio", &start_eta, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
          LALInferenceAddMinMaxPrior(priorArgs,	"massratio", &etaMin, &etaMax, LALINFERENCE_REAL8_t);
        }else{
          LALInferenceAddVariable(currentParams, "asym_massratio", &start_q, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
          LALInferenceAddMinMaxPrior(priorArgs,	"asym_massratio",	&qMin,	&qMax,	LALINFERENCE_REAL8_t);
        }
      }
    }
  }

  /* Set up start time. */
  ppt=LALInferenceGetProcParamVal(commandLine, "--time");
  if (ppt) {
    /* User has specified start time. */
    timeParam = atof(ppt->value);
  } else {
    timeParam = endtime+gsl_ran_gaussian(GSLrandom,0.01);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--fixTime");
  if(ppt){
    LALInferenceAddVariable(currentParams, "time",            &timeParam   ,           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    if(lalDebugLevel>0) fprintf(stdout,"time fixed and set to %f\n",timeParam);
  }else{
    LALInferenceAddVariable(currentParams, "time",            &timeParam   ,           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
  }
  LALInferenceAddMinMaxPrior(priorArgs, "time",     &timeMin, &timeMax,   LALINFERENCE_REAL8_t);

  ppt=LALInferenceGetProcParamVal(commandLine,"--fixPhi");
  if(ppt){
    LALInferenceAddVariable(currentParams, "phase",           &start_phase,        LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    if(lalDebugLevel>0) fprintf(stdout,"phase fixed and set to %f\n",start_phase);
  }else{
    LALInferenceAddVariable(currentParams, "phase",           &start_phase,        LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
  }
  LALInferenceAddMinMaxPrior(priorArgs, "phase",     &phiMin, &phiMax,   LALINFERENCE_REAL8_t);

  ppt=LALInferenceGetProcParamVal(commandLine,"--fixDist");
  if(ppt){
    LALInferenceAddVariable(currentParams,"distance", &start_dist, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    if(lalDebugLevel>0) fprintf(stdout,"distance fixed and set to %f\n",start_dist);
  }else{
    LALInferenceAddVariable(currentParams,"distance", &start_dist, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
  }
  LALInferenceAddMinMaxPrior(priorArgs, "distance",     &Dmin, &Dmax,   LALINFERENCE_REAL8_t);


  ppt=LALInferenceGetProcParamVal(commandLine,"--fixRa");
  if(ppt){
    LALInferenceAddVariable(currentParams, "rightascension",  &start_ra,      LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    if(lalDebugLevel>0) fprintf(stdout,"R.A. fixed and set to %f\n",start_ra);
  }else{
    LALInferenceAddVariable(currentParams, "rightascension",  &start_ra,      LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
  }
  LALInferenceAddMinMaxPrior(priorArgs, "rightascension",     &raMin, &raMax,   LALINFERENCE_REAL8_t);

  ppt=LALInferenceGetProcParamVal(commandLine,"--fixDec");
  if(ppt){
    LALInferenceAddVariable(currentParams, "declination",     &start_dec,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    if(lalDebugLevel>0) fprintf(stdout,"declination fixed and set to %f\n",start_dec);
  }else{
    LALInferenceAddVariable(currentParams, "declination",     &start_dec,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
  }
  LALInferenceAddMinMaxPrior(priorArgs, "declination",     &decMin, &decMax,   LALINFERENCE_REAL8_t);

  ppt=LALInferenceGetProcParamVal(commandLine,"--fixPsi");
  if(ppt){
    LALInferenceAddVariable(currentParams, "polarisation",    &start_psi,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    if(lalDebugLevel>0) fprintf(stdout,"polarisation fixed and set to %f\n",start_psi);
  }else{
    LALInferenceAddVariable(currentParams, "polarisation",    &start_psi,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
  }
  LALInferenceAddMinMaxPrior(priorArgs, "polarisation",     &psiMin, &psiMax,   LALINFERENCE_REAL8_t);

  ppt=LALInferenceGetProcParamVal(commandLine,"--fixIota");
  if(ppt){
    LALInferenceAddVariable(currentParams, "inclination",     &start_iota,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    if(lalDebugLevel>0) fprintf(stdout,"iota fixed and set to %f\n",start_iota);
  }else{
    LALInferenceAddVariable(currentParams, "inclination",     &start_iota,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
  }
  LALInferenceAddMinMaxPrior(priorArgs, "inclination",     &iotaMin, &iotaMax,   LALINFERENCE_REAL8_t);

  /* Check for spin disabled */
  ppt=LALInferenceGetProcParamVal(commandLine, "--noSpin");
  if (!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--diable-spin");
  if (ppt) noSpin=1;
  
  /* Check for aligned spin */
  ppt=LALInferenceGetProcParamVal(commandLine, "--spinAligned");
  if (!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--aligned-spin");
  if(ppt) spinAligned=1;
  
  /* Check for single spin */
  ppt=LALInferenceGetProcParamVal(commandLine,"--singleSpin");
  if(ppt) singleSpin=1;
  
  if((approx==SpinTaylor || approx==SpinTaylorFrameless || approx==PhenSpinTaylorRD || approx==SpinTaylorT4) && noSpin){
    if(spinAligned) a1min=-1.0;
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixA1");
    if(ppt){
      LALInferenceAddVariable(currentParams, "a_spin1",     &start_a_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      if(lalDebugLevel>0) fprintf(stdout,"spin 1 fixed and set to %f\n",start_a_spin1);
    }else{
      LALInferenceAddVariable(currentParams, "a_spin1",     &start_a_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "a_spin1",     &a1min, &a1max,   LALINFERENCE_REAL8_t);

    if(spinAligned) fprintf(stdout,"Running with spin1 aligned to the orbital angular momentum.\n");
    else {
      ppt=LALInferenceGetProcParamVal(commandLine,"--fixTheta1");
      if(ppt){
        LALInferenceAddVariable(currentParams, "theta_spin1",     &start_theta_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        if(lalDebugLevel>0) fprintf(stdout,"theta 1 fixed and set to %f\n",start_theta_spin1);
      }else{
        LALInferenceAddVariable(currentParams, "theta_spin1",     &start_theta_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
      }
      LALInferenceAddMinMaxPrior(priorArgs, "theta_spin1",     &theta1min, &theta1max,   LALINFERENCE_REAL8_t);

      ppt=LALInferenceGetProcParamVal(commandLine,"--fixPhi1");
      if(ppt){
        LALInferenceAddVariable(currentParams, "phi_spin1",     &start_phi_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        if(lalDebugLevel>0) fprintf(stdout,"phi 1 fixed and set to %f\n",start_phi_spin1);
      }else{
        LALInferenceAddVariable(currentParams, "phi_spin1",     &start_phi_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
      }
      LALInferenceAddMinMaxPrior(priorArgs, "phi_spin1",     &phi1min, &phi1max,   LALINFERENCE_REAL8_t);
    }
    if(singleSpin) fprintf(stdout,"Running with first spin set to 0\n");
    else {
      if(spinAligned) a2min=-1.0;
      ppt=LALInferenceGetProcParamVal(commandLine,"--fixA2");
      if(ppt){
        LALInferenceAddVariable(currentParams, "a_spin2",     &start_a_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        if(lalDebugLevel>0) fprintf(stdout,"spin 2 fixed and set to %f\n",start_a_spin2);
      }else{
        LALInferenceAddVariable(currentParams, "a_spin2",     &start_a_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
      }
      LALInferenceAddMinMaxPrior(priorArgs, "a_spin2",     &a2min, &a2max,   LALINFERENCE_REAL8_t);

      if(spinAligned) fprintf(stdout,"Running with spin2 aligned to the orbital angular momentum.\n");
      else {
        ppt=LALInferenceGetProcParamVal(commandLine,"--fixTheta2");
        if(ppt){
          LALInferenceAddVariable(currentParams, "theta_spin2",     &start_theta_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
          if(lalDebugLevel>0) fprintf(stdout,"theta spin 2 fixed and set to %f\n",start_theta_spin2);
        }else{
          LALInferenceAddVariable(currentParams, "theta_spin2",     &start_theta_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
        }
        LALInferenceAddMinMaxPrior(priorArgs, "theta_spin2",     &theta2min, &theta2max,   LALINFERENCE_REAL8_t);

        ppt=LALInferenceGetProcParamVal(commandLine,"--fixPhi2");
        if(ppt){
          LALInferenceAddVariable(currentParams, "phi_spin2",     &start_phi_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
          if(lalDebugLevel>0) fprintf(stdout,"phi 2 fixed and set to %f\n",start_phi_spin2);
        }else{
          LALInferenceAddVariable(currentParams, "phi_spin2",     &start_phi_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
        }
        LALInferenceAddMinMaxPrior(priorArgs, "phi_spin2",     &phi2min, &phi2max,   LALINFERENCE_REAL8_t);
      }
    }
  }

  /* Freq domain spinning templates */
  if((approx==TaylorF2 || approx==TaylorF2RedSpin || approx==TaylorF2RedSpinTidal || approx==IMRPhenomB) && spinAligned){

    a1min=-1.0;
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixA1");
    if(ppt){
      LALInferenceAddVariable(currentParams, "spin1",     &start_a_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      if(lalDebugLevel>0) fprintf(stdout,"spin 1 fixed and set to %f\n",start_a_spin1);
    }else{
      LALInferenceAddVariable(currentParams, "spin1",     &start_a_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "spin1",     &a1min, &a1max,   LALINFERENCE_REAL8_t);

    a2min=-1.0;
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixA2");
    if(ppt){
      LALInferenceAddVariable(currentParams, "spin2",     &start_a_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      if(lalDebugLevel>0) fprintf(stdout,"spin 2 fixed and set to %f\n",start_a_spin2);
    }else{
      LALInferenceAddVariable(currentParams, "spin2",     &start_a_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "spin2",     &a2min, &a2max,   LALINFERENCE_REAL8_t);

  }

  ppt=LALInferenceGetProcParamVal(commandLine, "--TaylorF2ppE");
  if(approx==TaylorF2 && ppt){

    REAL8 start_alpha, start_A, start_a, start_beta, start_B, start_b;

    tmpMin = -1000;
    tmpMax = 1000;
    start_alpha = tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
    ppt=LALInferenceGetProcParamVal(commandLine,"--ppealpha");
    if (ppt) {
      start_alpha = atof(ppt->value);
    }
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixppealpha");
    if(ppt){
      LALInferenceAddVariable(currentParams, "ppealpha",     &start_alpha,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      if(lalDebugLevel>0) fprintf(stdout,"ppE alpha fixed and set to %f\n",start_alpha);
    }else{
      LALInferenceAddVariable(currentParams, "ppealpha",     &start_alpha,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "ppealpha",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

    start_beta = tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
    ppt=LALInferenceGetProcParamVal(commandLine,"--ppebeta");
    if (ppt) {
      start_beta = atof(ppt->value);
    }
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixppebeta");
    if(ppt){
      LALInferenceAddVariable(currentParams, "ppebeta",     &start_beta,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      if(lalDebugLevel>0) fprintf(stdout,"ppE beta fixed and set to %f\n",start_beta);
    }else{
      LALInferenceAddVariable(currentParams, "ppebeta",     &start_beta,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "ppebeta",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

    tmpMin = -3;
    tmpMax = 3;
    start_A = tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
    ppt=LALInferenceGetProcParamVal(commandLine,"--ppeuppera");
    if (ppt) {
      start_A = atof(ppt->value);
    }
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixppeuppera");
    if(ppt){
      LALInferenceAddVariable(currentParams, "ppeuppera",     &start_A,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      if(lalDebugLevel>0) fprintf(stdout,"ppE A fixed and set to %f\n",start_A);
    }else{
      LALInferenceAddVariable(currentParams, "ppeuppera",     &start_A,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "ppeuppera",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

    start_B = tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
    ppt=LALInferenceGetProcParamVal(commandLine,"--ppeupperb");
    if (ppt) {
      start_B = atof(ppt->value);
    }
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixppeupperb");
    if(ppt){
      LALInferenceAddVariable(currentParams, "ppeupperb",     &start_B,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      if(lalDebugLevel>0) fprintf(stdout,"ppE B fixed and set to %f\n",start_B);
    }else{
      LALInferenceAddVariable(currentParams, "ppeupperb",     &start_B,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "ppeupperb",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

    tmpMin = -3.0;
    tmpMax = 2.0/3.0;
    start_a = tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
    ppt=LALInferenceGetProcParamVal(commandLine,"--ppelowera");
    if (ppt) {
      start_a = atof(ppt->value);
    }
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixppelowera");
    if(ppt){
      LALInferenceAddVariable(currentParams, "ppelowera",     &start_a,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      if(lalDebugLevel>0) fprintf(stdout,"ppE a fixed and set to %f\n",start_a);
    }else{
      LALInferenceAddVariable(currentParams, "ppelowera",     &start_a,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "ppelowera",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

    tmpMin = -4.5;
    tmpMax = 1.0;
    start_b = tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
    ppt=LALInferenceGetProcParamVal(commandLine,"--ppelowerb");
    if (ppt) {
      start_b = atof(ppt->value);
    }
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixppelowerb");
    if(ppt){
      LALInferenceAddVariable(currentParams, "ppelowerb",     &start_b,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      if(lalDebugLevel>0) fprintf(stdout,"ppE b fixed and set to %f\n",start_b);
    }else{
      LALInferenceAddVariable(currentParams, "ppelowerb",     &start_b,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "ppelowerb",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

  }
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--tidal");
  if(ppt){
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixLambda1");
    if(ppt){
      LALInferenceAddVariable(currentParams, "lambda1",           &start_lambda1,        LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      if(lalDebugLevel>0) fprintf(stdout,"phase fixed and set to %f\n",start_lambda1);
    }else{
      LALInferenceAddVariable(currentParams, "lambda1",           &start_lambda1,        LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "lambda1",     &lambda1Min, &lambda1Max,   LALINFERENCE_REAL8_t);
  
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixLambda2");
    if(ppt){
    LALInferenceAddVariable(currentParams, "lambda2",           &start_lambda2,        LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      if(lalDebugLevel>0) fprintf(stdout,"phase fixed and set to %f\n",start_lambda2);
    }else{
      LALInferenceAddVariable(currentParams, "lambda2",           &start_lambda2,        LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "lambda2",     &lambda2Min, &lambda2Max,   LALINFERENCE_REAL8_t);
  }
  
  LALSimInspiralInteraction interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_ALL;
  ppt=LALInferenceGetProcParamVal(commandLine,"--interactionFlags");
  if(ppt){
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_NONE")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_NONE;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_15PN")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_15PN;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_2PN")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_2PN;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_SELF_2PN")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_SELF_2PN;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_QUAD_MONO_2PN")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_QUAD_MONO_2PN;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_25PN")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_25PN;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_TIDAL_5PN")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_TIDAL_5PN;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_TIDAL_6PN")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_TIDAL_6PN;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_ALL_SPIN")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_ALL_SPIN;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_ALL")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_ALL;
    LALInferenceAddVariable(currentParams, "interactionFlags", &interactionFlags,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED); 
 }
 
}



/* Setup the variable for the evidence calculation test for review */
/* 5-sigma ranges for analytic likeliood function */
/* https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/LALInferenceReviewAnalyticGaussianLikelihood */
void LALInferenceInitVariablesReviewEvidence(LALInferenceRunState *state)
{
    ProcessParamsTable *commandLine=state->commandLine;
    ProcessParamsTable *ppt=NULL;
    char **strings=NULL;
    char *pinned_params=NULL;
    UINT4 N=0,i,j;
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--pinparams"))){
            pinned_params=ppt->value;
            LALInferenceVariables tempParams;
            memset(&tempParams,0,sizeof(tempParams));
            LALInferenceParseCharacterOptionString(pinned_params,&strings,&N);
    }
	LALInferenceVariables *priorArgs=state->priorArgs;
        state->currentParams=XLALCalloc(1,sizeof(LALInferenceVariables));
        LALInferenceVariables *currentParams=state->currentParams;
	i=0;

	struct varSettings {const char *name; REAL8 val, min, max;};
	
	struct varSettings setup[]=
	{
		{.name="time", .val=0.0, .min=-0.1073625, .max=0.1073625},
		{.name="m1", .val=16., .min=14.927715, .max=17.072285},
		{.name="m2", .val=7., .min=5.829675, .max=8.170325},
		{.name="distance", .val=50., .min=37.986000000000004, .max=62.013999999999996},
		{.name="inclination", .val=LAL_PI/2., .min=1.4054428267948966, .max=1.7361498267948965},
		{.name="phase", .val=LAL_PI, .min=2.8701521535897934, .max=3.413033153589793},
		{.name="polarisation", .val=LAL_PI/2., .min=1.3885563267948966, .max=1.7530363267948965},
		{.name="rightascension", .val=LAL_PI, .min=2.813050153589793, .max=3.4701351535897933},
		{.name="declination", .val=0., .min=-0.300699, .max=0.300699},
		{.name="a_spin1", .val=0.5, .min=0.3784565, .max=0.6215435},
		{.name="a_spin2", .val=0.5, .min=0.421869, .max=0.578131},
		{.name="theta_spin1", .val=LAL_PI/2., .min=1.3993998267948966, .max=1.7421928267948965},
		{.name="theta_spin2", .val=LAL_PI/2., .min=1.4086158267948965, .max=1.7329768267948966},
		{.name="phi_spin1", .val=LAL_PI, .min=2.781852653589793, .max=3.501332653589793},
		{.name="phi_spin2", .val=LAL_PI, .min=2.777215653589793, .max=3.5059696535897933},
		{.name="END", .val=0., .min=0., .max=0.}
	};

	while(strcmp("END",setup[i].name))
	{
        LALInferenceParamVaryType type=LALINFERENCE_PARAM_CIRCULAR;
        /* Check if it is to be fixed */
        for(j=0;j<N;j++) if(!strcmp(setup[i].name,strings[j])) {type=LALINFERENCE_PARAM_FIXED; printf("Fixing parameter %s\n",setup[i].name); break;}
		LALInferenceAddVariable(currentParams,setup[i].name, &(setup[i].val) ,LALINFERENCE_REAL8_t, type);
	    LALInferenceAddMinMaxPrior(priorArgs, setup[i].name,    &(setup[i].min),    &(setup[i].max),    LALINFERENCE_REAL8_t);
		i++;
	}
	return;
}


void LALInferenceInitVariablesReviewEvidence_bimod(LALInferenceRunState *state)
{
  ProcessParamsTable *commandLine=state->commandLine;
  ProcessParamsTable *ppt=NULL;
  char **strings=NULL;
  char *pinned_params=NULL;
  UINT4 N=0,i,j;
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--pinparams"))){
    pinned_params=ppt->value;
    LALInferenceVariables tempParams;
    memset(&tempParams,0,sizeof(tempParams));
    LALInferenceParseCharacterOptionString(pinned_params,&strings,&N);
  }
  LALInferenceVariables *priorArgs=state->priorArgs;
  state->currentParams=XLALCalloc(1,sizeof(LALInferenceVariables));
  LALInferenceVariables *currentParams=state->currentParams;
  i=0;
  
  struct varSettings {const char *name; REAL8 val, min, max;};
  
  struct varSettings setup[]=
  {
    {.name="time", .val=0.05589, .min=-0.1373625, .max=0.2491425},
    {.name="m1", .val=16.857828, .min=14.927715, .max=18.787941},
    {.name="m2", .val=7.93626, .min=5.829675, .max=10.042845},
    {.name="distance", .val=34.6112, .min=12.986, .max=56.2364},
    {.name="inclination", .val=0.9176809634, .min=0.6200446634, .max=1.2153172634},
    {.name="phase", .val=1.7879487268, .min=1.2993558268, .max=2.2765416268},
    {.name="polarisation", .val=0.9311901634, .min=0.6031581634, .max=1.2592221634},
    {.name="rightascension", .val=1.8336303268, .min=1.2422538268, .max=2.4250068268},
    {.name="declination", .val=-0.5448389634, .min=-1.0860971634, .max=-0.0035807634},
    {.name="a_spin1", .val=0.2972348, .min=0.0784565, .max=0.5160131},
    {.name="a_spin2", .val=0.2625048, .min=0.121869, .max=0.4031406},
    {.name="theta_spin1", .val=0.9225153634, .min=0.6140016634, .max=1.2310290634},
    {.name="theta_spin2", .val=0.9151425634, .min=0.6232176634, .max=1.2070674634},
    {.name="phi_spin1", .val=1.8585883268, .min=1.2110563268, .max=2.5061203268},
    {.name="phi_spin2", .val=1.8622979268, .min=1.2064193268, .max=2.5181765268},
    {.name="END", .val=0., .min=0., .max=0.}
  };
  
  while(strcmp("END",setup[i].name))
  {
    LALInferenceParamVaryType type=LALINFERENCE_PARAM_CIRCULAR;
    /* Check if it is to be fixed */
    for(j=0;j<N;j++) if(!strcmp(setup[i].name,strings[j])) {type=LALINFERENCE_PARAM_FIXED; printf("Fixing parameter %s\n",setup[i].name); break;}
    LALInferenceAddVariable(currentParams,setup[i].name, &(setup[i].val) ,LALINFERENCE_REAL8_t, type);
    LALInferenceAddMinMaxPrior(priorArgs, setup[i].name,    &(setup[i].min),    &(setup[i].max),    LALINFERENCE_REAL8_t);
    i++;
  }
  return;
}

void LALInferenceInitVariablesReviewEvidence_banana(LALInferenceRunState *state)
{
  ProcessParamsTable *commandLine=state->commandLine;
  ProcessParamsTable *ppt=NULL;
  char **strings=NULL;
  char *pinned_params=NULL;
  UINT4 N=0,i,j;
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--pinparams"))){
    pinned_params=ppt->value;
    LALInferenceVariables tempParams;
    memset(&tempParams,0,sizeof(tempParams));
    LALInferenceParseCharacterOptionString(pinned_params,&strings,&N);
  }
  LALInferenceVariables *priorArgs=state->priorArgs;
  state->currentParams=XLALCalloc(1,sizeof(LALInferenceVariables));
  LALInferenceVariables *currentParams=state->currentParams;
  i=0;
  
  struct varSettings {const char *name; REAL8 val, min, max;};
  
  struct varSettings setup[]=
  {
    {.name="time", .val=0.0, .min=-2., .max=2.},
    {.name="m1", .val=16., .min=14., .max=18.},
    {.name="m2", .val=7., .min=5., .max=9.},
    {.name="distance", .val=50., .min=48., .max=52.},
    {.name="inclination", .val=LAL_PI/2., .min=-0.429203673, .max=3.570796327},
    {.name="phase", .val=LAL_PI, .min=1.141592654, .max=5.141592654},
    {.name="polarisation", .val=LAL_PI/2., .min=-0.429203673, .max=3.570796327},
    {.name="rightascension", .val=LAL_PI, .min=1.141592654, .max=5.141592654},
    {.name="declination", .val=0., .min=-2., .max=2.},
    {.name="a_spin1", .val=0.5, .min=-1.5, .max=2.5},
    {.name="a_spin2", .val=0.5, .min=-1.5, .max=2.5},
    {.name="theta_spin1", .val=LAL_PI/2., .min=-0.429203673, .max=3.570796327},
    {.name="theta_spin2", .val=LAL_PI/2., .min=-0.429203673, .max=3.570796327},
    {.name="phi_spin1", .val=LAL_PI, .min=1.141592654, .max=5.141592654},
    {.name="phi_spin2", .val=LAL_PI, .min=1.141592654, .max=5.141592654},
    {.name="END", .val=0., .min=0., .max=0.}
  };
  
  while(strcmp("END",setup[i].name))
  {
    LALInferenceParamVaryType type=LALINFERENCE_PARAM_CIRCULAR;
    /* Check if it is to be fixed */
    for(j=0;j<N;j++) if(!strcmp(setup[i].name,strings[j])) {type=LALINFERENCE_PARAM_FIXED; printf("Fixing parameter %s\n",setup[i].name); break;}
    LALInferenceAddVariable(currentParams,setup[i].name, &(setup[i].val) ,LALINFERENCE_REAL8_t, type);
    LALInferenceAddMinMaxPrior(priorArgs, setup[i].name,    &(setup[i].min),    &(setup[i].max),    LALINFERENCE_REAL8_t);
    i++;
  }
  return;
}


