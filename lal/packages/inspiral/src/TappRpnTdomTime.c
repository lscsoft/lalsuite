/*------------------------------------------------------------------------------
 *    
 *    File Name: TappRpnTdomTime.c
 *
 *    Created 25/08/9 
 *
 *    Authors: B. S. Sathyaprakash and D. K. Churches
 *
 *    Purpose: To calculate the waveform expected from an inspiralling binary
 *             system up to second post-Newtonian order, using an expansion
 *             in terms of time	
 *
 *    Diagnostics: 
 *
 *    Inputs:	*status	: pointer to our status structure
 *
 *              *params : pointer to an input parameter structure
 *
 *    Output:   *output : pointer to the output structure which
 *                        will contain the chirp waveform
 *
 *    Calls: Tapp_PRN_tdom_time_phase, Tapp_PRN_tdom_time_frequency
 *
 *    Comments: 
 *
 *    Notes: 
 *
 *    Acknowledgements: 
 *
 *    Revision History:  
 */

#include <lal/LALStdlib.h>
#include <lal/Inspiral.h>



NRCSID (TAPPRPNTDOMTIMEC, "$Id$");


void (*TappRpnTdomTimeFrequency) (LALStatus *status,
                                  InspiralwaveFrequencyOutput *output,
			          InspiralwaveFrequencyInput *params);
void (*TappRpnTdomTimePhase) (LALStatus *status,
                              InspiralwavePhaseOutput *output,
			      InspiralwavePhaseInput *params);

void LALTappRpnTdomTime (LALStatus *status,
                      REAL8Vector *output, 
		      InspiralTemplate *params)
{

  INT4 n, m, i, startShift, endShift, count;
  REAL8 dt, fu, eta, tc, totalMass, t, c1, phi0, phi;
  REAL8 v, x, q, fLso, fHigh, amp, tmax;

/* Define the structures which will be used */

  InspiralwavePhaseInput input1;       /*Input structure for TappRpnTdomTimePhase()*/
  InspiralwavePhaseOutput output1;     /*Output structure for TappRpnTdomTimePhase()*/
  InspiralwaveFrequencyInput input2;   /*Input structure for TappRpnTdomTimeFrequency()*/
  InspiralwaveFrequencyOutput output2; /*Output structure for TappRpnTdomTimeFrequency()*/ 
  InspiralParamsInput paramCalc;       /*Input structure for inspiralParameterCalc()*/
  InspiralParamsOutput out1;           /*Output structure for inspiralParameterCalc()*/




  INITSTATUS (status, "LALTappRpnTdomTime", TAPPRPNTDOMTIMEC);
  ATTATCHSTATUSPTR(status);

  ASSERT(output, status, TAPPRPNTDOMTIME_ENULL, TAPPRPNTDOMTIME_MSGENULL);
  ASSERT(output->data, status, TAPPRPNTDOMTIME_ENULL, TAPPRPNTDOMTIME_MSGENULL);
  ASSERT(params, status, TAPPRPNTDOMTIME_ENULL, TAPPRPNTDOMTIME_MSGENULL); 
  ASSERT(params->nStartPad >= 0, status, TAPPRPNTDOMTIME_ESIZE, TAPPRPNTDOMTIME_MSGESIZE);
  ASSERT(params->nEndPad >= 0, status, TAPPRPNTDOMTIME_ESIZE, TAPPRPNTDOMTIME_MSGESIZE);
  ASSERT(params->fLower > 0, status, TAPPRPNTDOMTIME_ESIZE, TAPPRPNTDOMTIME_MSGESIZE);
  ASSERT(params->fCutoff > 0, status, TAPPRPNTDOMTIME_ESIZE, TAPPRPNTDOMTIME_MSGESIZE);
  ASSERT(params->fCutoff > params->fLower, status, TAPPRPNTDOMTIME_ESIZE, TAPPRPNTDOMTIME_MSGESIZE);
  ASSERT(params->tSampling > 0, status, TAPPRPNTDOMTIME_ESIZE, TAPPRPNTDOMTIME_MSGESIZE);


  

  dt = 1.0/(params->tSampling);          /* The sampling rate  */
  fu = params->fCutoff;            /* The upper frequency cutoff  */
  phi = params->startPhase;        /* The initial phase  */
  startShift = params->nStartPad;  /* The number of zeros at the start of the wave  */
  endShift = params->nEndPad;      /* The number of zeros at the end of the wave  */

/* Initialize the members of the structure which will be used as the input
   to the function inspiralParameterCalc which will calculate the chirptimes etc  */

  paramCalc.m1 = params->mass1;
  paramCalc.m2 = params->mass2;
  paramCalc.totalMass = params->totalMass;
  paramCalc.eta = params->eta;
  paramCalc.mu = params->mu;
  paramCalc.fLower = params->fLower;
  paramCalc.massChoice = params->massChoice;
  paramCalc.order = params->order;

  LALInspiralParameterCalc (status->statusPtr, &out1, &paramCalc);
  CHECKSTATUSPTR(status);


/* Now check that the outputted values from inspiralParameterCalc are within
   the allowed limits  */

  ASSERT(out1.totalMass > 0.4, status, TAPPRPNTDOMTIME_ESIZE, TAPPRPNTDOMTIME_MSGESIZE);
  ASSERT(out1.totalMass < 100, status, TAPPRPNTDOMTIME_ESIZE, TAPPRPNTDOMTIME_MSGESIZE);
  ASSERT(out1.eta >= 0, status, TAPPRPNTDOMTIME_ESIZE, TAPPRPNTDOMTIME_MSGESIZE);
  ASSERT(out1.eta <=0.25, status, TAPPRPNTDOMTIME_ESIZE, TAPPRPNTDOMTIME_MSGESIZE);
  ASSERT(out1.mu >= 0, status, TAPPRPNTDOMTIME_ESIZE, TAPPRPNTDOMTIME_MSGESIZE);


  switch (params->order) {
     case newtonian:
     case oneHalfPN:
          TappRpnTdomTimePhase = &LALTappRpnTdomTimePhase0PN;
          TappRpnTdomTimeFrequency = &LALTappRpnTdomTimeFrequency0PN;
          break;
     case onePN:
          TappRpnTdomTimePhase = &LALTappRpnTdomTimePhase2PN;
          TappRpnTdomTimeFrequency = &LALTappRpnTdomTimeFrequency2PN;
          break;
     case onePointFivePN:
          TappRpnTdomTimePhase = &LALTappRpnTdomTimePhase3PN;
          TappRpnTdomTimeFrequency = &LALTappRpnTdomTimeFrequency3PN;
          break;
     case twoPN:
          TappRpnTdomTimePhase = &LALTappRpnTdomTimePhase4PN;
          TappRpnTdomTimeFrequency = &LALTappRpnTdomTimeFrequency4PN;
          break;
     default:
          fprintf(stderr, "No order selected in LALTappRpnTdomTime ... exiting\n");
          exit(0);
     }
  
/* Call the function inspiralParameterCalc  */
 
  tc=out1.tauC;              /* Instant of coalescence of the compact objects */
  eta = out1.eta;                              /* Symmetric mass ratio  */
  totalMass = (out1.totalMass)*LAL_MTSUN_SI;   /* The mass of the system in seconds */


/* Calculate the number of data points needed to define the waveform. 
   The actual number used is the next highest power of 2  */

  x = (tc/dt)+startShift+endShift;      /* Actual number of points needed  */
  q = log10(x)/log10(2.0);
  m = ceil(q);
  n = pow(2.0,m);                       /* The next highest power of 2  */


/* Allocate memory for the waveform */

  output->data = (REAL8 *)LALMalloc(sizeof(REAL8)*n);


/* Calculate the frequency of the last stable orbit flso. If flso is less 
   than the user inputted upper frequency cutoff fu, then the waveforn is 
   truncated at f=flso.  If fu is less than flso, then we truncate the 
   waveform at fu. */

  fLso = 1.0/(LAL_PI*totalMass*pow(6.0, 1.5));
  if (fu < fLso) 
     fHigh = fu;
  else 
     fHigh = fLso;

/* Check that the highest frequency is less than half the sampling frequency - 
   the Nyquist theorum */

  ASSERT(fHigh < 0.5/dt, status, TAPPRPNTDOMTIME_ESIZE, TAPPRPNTDOMTIME_MSGESIZE);
  ASSERT(fHigh > params->fLower, status, TAPPRPNTDOMTIME_ESIZE, TAPPRPNTDOMTIME_MSGESIZE);

/* Initialize members of the structures which get fed into TappRpnTdomTimePhase()
   and TappRpnTdomTimeFrequency(). */

  input1.etaby2 = 0.5 * eta;
  input1.a2 = 3715./8064. + 55.*eta/96.;
  input1.a3 = 0.75*LAL_PI;
  input1.a4 = 9275495./14450688. + 284875.*eta/258048. + 1855.*pow(eta,2.0)/2048.;

  input2.eightPiM = 8.*LAL_PI*totalMass;
  input2.a2 = 743./2688.+(11.*eta)/32.;
  input2.a3 = 0.3*LAL_PI;
  input2.a4 =  1855099./14450688. +  56975.*eta/258048. + 371.*pow(eta,2.0)/2048;

/* Here's the part which calculates the waveform */

  c1 = eta/(5.*totalMass);
  i=0; while (i<startShift) output->data[i++] = 0.0;

  t=0.0;
  input1.td = input2.td = c1*(tc-t);
  TappRpnTdomTimePhase(status->statusPtr, &output1, &input1);
  CHECKSTATUSPTR(status);
  TappRpnTdomTimeFrequency(status->statusPtr, &output2, &input2);
  CHECKSTATUSPTR(status);
  phi0=-output1.phase+phi;
  i = 0;
  while (i < params->nStartPad) *(output->data + i++) = 0.;

  count = 0;
/*
  c2 = 2. * out1.eta * out1.totalMass;
*/
  tmax = tc - dt;
  while (output2.frequency<fHigh && t<tmax) 
  {
    v = pow(output2.frequency*LAL_PI*totalMass, oneby3);
    amp = v*v; 
    output->data[i++]=params->signalAmplitude * amp * cos(output1.phase+phi0);
    ++count;
    t=count*dt;
    input1.td = input2.td = c1*(tc-t);
    TappRpnTdomTimeFrequency(status->statusPtr, &output2, &input2);
    CHECKSTATUSPTR(status); 
    TappRpnTdomTimePhase(status->statusPtr, &output1, &input1);
    CHECKSTATUSPTR(status);
  }

  while (i<n) output->data[i++]=0.0;
  output->length=i;


  DETATCHSTATUSPTR(status);
  RETURN(status);
}


