#include <lal/LALStdlib.h>
#include <lal/Inspiral.h>



NRCSID (TAPPRPNTDOMTIMEFREQUENCYC, "$Id$");


void LALTappRpnTdomTimeFrequency0PN (LALStatus *status,
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

void LALTappRpnTdomTimeFrequency1PN (LALStatus *status,
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

void LALTappRpnTdomTimeFrequency2PN (LALStatus *status,
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

void LALTappRpnTdomTimeFrequency3PN (LALStatus *status,
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

void LALTappRpnTdomTimeFrequency4PN (LALStatus *status,
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
