#include "LALStdlib.h"
#include "Inspiral.h"



NRCSID (TAPPRPNTDOMTIMEFREQUENCYC, "$Id$");


void TappRpnTdomTimeFrequency0PN (Status *status,
                                  InspiralwaveFrequencyOutput *output,
			          InspiralwaveFrequencyInput *params) 
{

  INITSTATUS (status, "TappRpnTdomTimeFrequency", TAPPRPNTDOMTIMEFREQUENCYC);

  ASSERT(output, status, TAPPRPNTDOMTIMEFREQUENCY_ENULL, TAPPRPNTDOMTIMEFREQUENCY_MSGENULL);
  ASSERT(params, status, TAPPRPNTDOMTIMEFREQUENCY_ENULL, TAPPRPNTDOMTIMEFREQUENCY_MSGENULL);
  ASSERT(params->td > 0, status, TAPPRPNTDOMTIMEFREQUENCY_ESIZE, TAPPRPNTDOMTIMEFREQUENCY_MSGESIZE);


  output->frequency = pow(params->td,-threeby8)/params->eightPiM;

  RETURN(status);
}

void TappRpnTdomTimeFrequency1PN (Status *status,
                                  InspiralwaveFrequencyOutput *output,
			          InspiralwaveFrequencyInput *params) 
{

  INITSTATUS (status, "TappRpnTdomTimeFrequency", TAPPRPNTDOMTIMEFREQUENCYC);

  ASSERT(output, status, TAPPRPNTDOMTIMEFREQUENCY_ENULL, TAPPRPNTDOMTIMEFREQUENCY_MSGENULL);
  ASSERT(params, status, TAPPRPNTDOMTIMEFREQUENCY_ENULL, TAPPRPNTDOMTIMEFREQUENCY_MSGENULL);
  ASSERT(params->td > 0, status, TAPPRPNTDOMTIMEFREQUENCY_ESIZE, TAPPRPNTDOMTIMEFREQUENCY_MSGESIZE);


  output->frequency = pow(params->td,-threeby8)/params->eightPiM;

  RETURN(status);
}

void TappRpnTdomTimeFrequency2PN (Status *status,
                                  InspiralwaveFrequencyOutput *output,
			          InspiralwaveFrequencyInput *params) 
{

  INITSTATUS (status, "TappRpnTdomTimeFrequency", TAPPRPNTDOMTIMEFREQUENCYC);

  ASSERT(output, status, TAPPRPNTDOMTIMEFREQUENCY_ENULL, TAPPRPNTDOMTIMEFREQUENCY_MSGENULL);
  ASSERT(params, status, TAPPRPNTDOMTIMEFREQUENCY_ENULL, TAPPRPNTDOMTIMEFREQUENCY_MSGENULL);
  ASSERT(params->td > 0, status, TAPPRPNTDOMTIMEFREQUENCY_ESIZE, TAPPRPNTDOMTIMEFREQUENCY_MSGESIZE);

  output->frequency = (pow(params->td,-threeby8) 
                    + params->a2*pow(params->td,-fiveby8))
                    /params->eightPiM;

  RETURN(status);
}

void TappRpnTdomTimeFrequency3PN (Status *status,
                                  InspiralwaveFrequencyOutput *output,
			          InspiralwaveFrequencyInput *params) 
{

  INITSTATUS (status, "TappRpnTdomTimeFrequency", TAPPRPNTDOMTIMEFREQUENCYC);

  ASSERT(output, status, TAPPRPNTDOMTIMEFREQUENCY_ENULL, TAPPRPNTDOMTIMEFREQUENCY_MSGENULL);
  ASSERT(params, status, TAPPRPNTDOMTIMEFREQUENCY_ENULL, TAPPRPNTDOMTIMEFREQUENCY_MSGENULL);
  ASSERT(params->td > 0, status, TAPPRPNTDOMTIMEFREQUENCY_ESIZE, TAPPRPNTDOMTIMEFREQUENCY_MSGESIZE);

  output->frequency = (pow(params->td,-threeby8) 
                    + params->a2*pow(params->td,-fiveby8)
                    - params->a3*pow(params->td,-threeby4))
                    /params->eightPiM;

  RETURN(status);
}

void TappRpnTdomTimeFrequency4PN (Status *status,
                                  InspiralwaveFrequencyOutput *output,
			          InspiralwaveFrequencyInput *params) 
{


  INITSTATUS (status, "TappRpnTdomTimeFrequency", TAPPRPNTDOMTIMEFREQUENCYC);

  ASSERT(output, status, TAPPRPNTDOMTIMEFREQUENCY_ENULL, TAPPRPNTDOMTIMEFREQUENCY_MSGENULL);
  ASSERT(params, status, TAPPRPNTDOMTIMEFREQUENCY_ENULL, TAPPRPNTDOMTIMEFREQUENCY_MSGENULL);
  ASSERT(params->td > 0, status, TAPPRPNTDOMTIMEFREQUENCY_ESIZE, TAPPRPNTDOMTIMEFREQUENCY_MSGESIZE);

  output->frequency = (pow(params->td,-threeby8) 
                    + params->a2*pow(params->td,-fiveby8)
                    - params->a3*pow(params->td,-threeby4)
                    + params->a4*pow(params->td,-sevenby8))
                    /params->eightPiM;

  RETURN(status);
}
