/*------------------------------------------------------------------------------
 *    
 *    File Name: TappRpnTdomFreq.c
 *
 *    Created 25/08/9 
 *
 *    Authors: B. S. Sathyaprakash and D. K. Churches
 *
 *    Purpose: To calculate the waveform expected from an inspiralling binary
 *             system up to second post-Newtonian order, using an expansion
 *             in terms of frequency	
 *
 *    Diagnostics: Null pointer, invalid input size
 *
 *    Inputs:	*status	: pointer to our status structure
 *
 *              *params : pointer to an input parameter structure
 *
 *    Output:   *output : pointer to the output structure which
 *                        will contain the chirp waveform
 *
 *    Calls: InspiralParameterCalc,
 *           TappRpnTdomFreqPhase
 *
 *    Comments: See Sathyaprakash, PRD, 50, R7111, 1994, Eq.(5) or document for
 *		Tapp_RPN_tdom_freq.c for further details.
 *
 *    Notes: 
 *
 *    Acknowledgements: 
 *
 *    Revision History:  
 */



#include <lal/LALStdlib.h>
#include <lal/Inspiral.h>
#include <lal/FindRoot.h>



NRCSID (TAPPRPNTDOMFREQC, "$Id$");

void (*TappRpnTdomFreqPhase) (LALStatus *status,
                              REAL8 *phase, 
                              InspiralPhasesInput *params);

void (*TappRpnTdomFreqTofF) (LALStatus *status,
                             REAL8 *toff,
			     REAL8 f,
                             void *params);


void LALTappRpnTdomFreq(LALStatus *status, 
                     REAL8Vector *output, 
                     InspiralTemplate *params)
{

  REAL8 dt, fs, fu, fsPi, fHigh, phase0;
  REAL8 phase, v, x, q, totalMass, fLso, freq;
  INT4 i, n, startShift, endShift, m, count;
  InspiralParamsOutput out1;
  InspiralPhasesInput phaseIn;
  InspiralParamsInput paramCalc;
  DFindRootIn rootIn;
  InspiralToffInput toffIn;
  void *funcParams;


  INITSTATUS (status, "LALTappRpnTdomFreq", TAPPRPNTDOMFREQC);
  ATTATCHSTATUSPTR(status);


  ASSERT(output,status,TAPPRPNTDOMFREQ_ENULL,TAPPRPNTDOMFREQ_MSGENULL);
  ASSERT(output->data,status,TAPPRPNTDOMFREQ_ENULL,TAPPRPNTDOMFREQ_MSGENULL);
  ASSERT(params,status,TAPPRPNTDOMFREQ_ENULL,TAPPRPNTDOMFREQ_MSGENULL);


  ASSERT(params->nStartPad >= 0, status, TAPPRPNTDOMFREQ_ESIZE, TAPPRPNTDOMFREQ_MSGESIZE);
  ASSERT(params->nEndPad >= 0, status, TAPPRPNTDOMFREQ_ESIZE, TAPPRPNTDOMFREQ_MSGESIZE);
  ASSERT(params->fLower > 0, status, TAPPRPNTDOMFREQ_ESIZE, TAPPRPNTDOMFREQ_MSGESIZE);
  ASSERT(params->fCutoff >= 0, status, TAPPRPNTDOMFREQ_ESIZE, TAPPRPNTDOMFREQ_MSGESIZE);
  ASSERT(params->fCutoff > params->fLower, status, TAPPRPNTDOMFREQ_ESIZE, TAPPRPNTDOMFREQ_MSGESIZE);
  ASSERT(params->tSampling > 0, status, TAPPRPNTDOMFREQ_ESIZE, TAPPRPNTDOMFREQ_MSGESIZE);


  dt = 1.0/(params->tSampling);
  fs = params->fLower;
  fu = params->fCutoff;
  fsPi = fs*LAL_PI;
  startShift = params->nStartPad;
  endShift = params->nEndPad;
  phase0 = params->startPhase;



/* Initialize the members of the structure which will be used as the input
   to the function which will calculate the chirptimes etc   */

  paramCalc.m1 = params->mass1;
  paramCalc.m2 = params->mass2;
  paramCalc.totalMass = params->totalMass;
  paramCalc.eta = params->eta;
  paramCalc.mu = params->mu;
  paramCalc.fLower = params->fLower;
  paramCalc.massChoice = params->massChoice;
  paramCalc.order = params->order;
  
  switch (params->order) {
     case newtonian:
     case oneHalfPN:
          TappRpnTdomFreqPhase = &LALTappRpnTdomFreqPhase0PN;
          TappRpnTdomFreqTofF = &LALTappRpnTdomFreqTofF0PN;
          rootIn.function = LALTappRpnTdomFreqTofF0PN;
          break;
     case onePN:
          TappRpnTdomFreqPhase = &LALTappRpnTdomFreqPhase2PN;
          TappRpnTdomFreqTofF = &LALTappRpnTdomFreqTofF2PN;
          rootIn.function = LALTappRpnTdomFreqTofF2PN;
          break;
     case onePointFivePN:
          TappRpnTdomFreqPhase = &LALTappRpnTdomFreqPhase3PN;
          TappRpnTdomFreqTofF = &LALTappRpnTdomFreqTofF3PN;
          rootIn.function = LALTappRpnTdomFreqTofF3PN;
          break;
     case twoPN:
          TappRpnTdomFreqPhase = &LALTappRpnTdomFreqPhase4PN;
          TappRpnTdomFreqTofF = &LALTappRpnTdomFreqTofF4PN;
          rootIn.function = LALTappRpnTdomFreqTofF4PN;
          break;
     default:
          fprintf(stderr, "No order selected in LALTappRpnTdomTime ... exiting\n");
          exit(0);
     }
/* Call the function which calculates all the chirptimes etc. The output
   from this function is in the structure out1 */


  LALInspiralParameterCalc (status->statusPtr, &out1, &paramCalc);
  CHECKSTATUSPTR(status);


/* Now check that the outputted values from inspiralParameterCalc are within
   the allowed limits  */

  ASSERT(out1.totalMass > 0.4, status, TAPPRPNTDOMFREQ_ESIZE, TAPPRPNTDOMFREQ_MSGESIZE);
  ASSERT(out1.totalMass < 100, status, TAPPRPNTDOMFREQ_ESIZE, TAPPRPNTDOMFREQ_MSGESIZE);
  ASSERT(out1.eta >= 0, status, TAPPRPNTDOMFREQ_ESIZE, TAPPRPNTDOMFREQ_MSGESIZE);
  ASSERT(out1.eta <=0.25, status, TAPPRPNTDOMFREQ_ESIZE, TAPPRPNTDOMFREQ_MSGESIZE);
  ASSERT(out1.mu >= 0, status, TAPPRPNTDOMFREQ_ESIZE, TAPPRPNTDOMFREQ_MSGESIZE);



/* Initialize the members of the structure which will be used as the input to
   the function which will calculate t as a function of f  */


  toffIn.t0 = out1.tau0;
  toffIn.t2 = out1.tau2;
  toffIn.t3 = out1.tau3;
  toffIn.t4 = out1.tau4;
  toffIn.tc = out1.tauC;



/* Initialize the members of the structure which will be used as the input to
   the function which will calcute the phase of the wave  */

  phaseIn.p0 = 3.2 * fsPi * out1.tau0;
  phaseIn.p2 = 4.0 * fsPi * out1.tau2;
  phaseIn.p3 = 5.0 * fsPi * out1.tau3;
  phaseIn.p4 = 8.0 * fsPi * out1.tau4;
  phaseIn.pc = phaseIn.p0 + phaseIn.p2 - phaseIn.p3 + phaseIn.p4;


/* Calculate how many data points will be needed to specify the
   waveform.  */

  x = (out1.tauC/dt)+startShift+endShift;
  q = log10(x)/log10(2.0);
  m = ceil(q);
  n = pow(2.0,m);

  output->data = (REAL8 *)LALMalloc(sizeof(REAL8)*n);


  totalMass = out1.totalMass*LAL_MTSUN_SI;

  fLso = 1.0/(LAL_PI*totalMass*pow(6.0, (3.0/2.0)));
  if (fu < fLso) fHigh = fu/fs;
  else fHigh = fLso/fs;

  ASSERT(fHigh*fs < 0.5/dt, status, TAPPRPNTDOMFREQ_ESIZE, TAPPRPNTDOMFREQ_MSGESIZE);
  ASSERT(fHigh*fs > params->fLower, status, TAPPRPNTDOMFREQ_ESIZE, TAPPRPNTDOMFREQ_MSGESIZE);


  toffIn.t = 0.0;

/* Initialize the members of the structure which will be used as the input to
   the function which will calcute the frequency of the wave  */


  rootIn.xmax = 1.1*fu/fs;
  rootIn.xacc = 1.0e-8;
  rootIn.xmin = 0.999999;


  i=0;
  while (i<startShift) output->data[i++] = 0.0;



/* Now cast the input structure to argument 4 of BisectionFindRoot so that it 
  of type void * rather than InspiralToffInput  */

  funcParams = (void *) &toffIn;

  LALDBisectionFindRoot(status->statusPtr, &freq, &rootIn, funcParams);
  CHECKSTATUSPTR(status);

  
  phaseIn.f=freq;
  count=1;

    while (freq < fHigh) {
    v = pow(freq*LAL_PI*totalMass, oneby3);
    TappRpnTdomFreqPhase(status->statusPtr, &phase, &phaseIn);
    CHECKSTATUSPTR(status);
    output->data[i++] = params->signalAmplitude * v*v * cos(phase+phase0);
    toffIn.t=count*dt;
    ++count;
    funcParams = (void *) &toffIn;
    LALDBisectionFindRoot(status->statusPtr, &freq, &rootIn, funcParams);
    CHECKSTATUSPTR(status);
    phaseIn.f=freq;
  }

  while (i<n) output->data[i++]=0.0;
  output->length=n;

  DETATCHSTATUSPTR(status);
  RETURN(status);

}
