/*-----------------------------------------------------------------------*
 * File Name: PowerETG.c
 *
 * Author: Julien Sylvestre
 *
 * Revision: $Id$ 
 *
 *-----------------------------------------------------------------------*/



/******** <lalVerbatim file="PowerETGCV"> ********
Author: Sylvestre, J
$Id$
********* </lalVerbatim> ********/

#include <config.h>
#include <lal/StdBurstSearch.h>
#include <lal/Sort.h>
#include <lal/LIGOMetadataTables.h>

#include <lal/AVFactories.h>
#include <lal/EPSearch.h>
#include <lal/LALRCSID.h>
#include <string.h>

#define SetStringParameter(par) if(!(params)) {ABORT(status, STDBURSTSEARCHH_ENULLPI, STDBURSTSEARCHH_MSGENULLPI);} \
  ASSERT ( par, status, STDBURSTSEARCHH_ENULLP, STDBURSTSEARCHH_MSGENULLP); \
  if(params->char_) { \
    strcpy(par,params->char_); \
  } else if(params->int4_) { \
    sprintf(par,"%i",*(params->int4_)); \
  } else if(params->real4_) { \
    sprintf(par,"%g",*(params->real4_)); \
  } else {ABORT(status, STDBURSTSEARCHH_ENULLPI, STDBURSTSEARCHH_MSGENULLPI);} \
  params = params->next
 
NRCSID (POWERETGC, "$Id$");

/******** <lalLaTeX file="PowerETGC"> ********
\noindent
Implement the POWER event trigger generator.
\subsubsection*{Prototype}
********* </lalLaTeX> ********/
/* <lalVerbatim> */
void
LALPowerETG(
		 LALStatus *status, 
		 EventIDColumn *output, 
		 REAL4TimeVectorSeries *input, 
		 BurstParameter *params
	       ) {
/* </lalVerbatim> */
/******** <lalLaTeX file="PowerETGC"> ********
\subsubsection*{Description}
Description of the parameters: \\
same vector of string as the power dso.

\subsubsection*{Uses}
\begin{verbatim}
...a bunch of stuff.
\end{verbatim}
********* </lalLaTeX> ********/
  CHAR buffer[1024];
  INT4 i, argc = 0;
  CHAR **argv = NULL;
  EPSearchParams *EPparams;
  REAL4TimeSeries tmp;
  REAL4Vector tv;
  SnglBurstTable *burstEvent = NULL;

  INITSTATUS (status, "LALPowerETG", POWERETGC);
  ATTATCHSTATUSPTR (status);

  /* trivial input check */
  ASSERT ( output, status, STDBURSTSEARCHH_ENULLP, STDBURSTSEARCHH_MSGENULLP);
  ASSERT ( input, status, STDBURSTSEARCHH_ENULLP, STDBURSTSEARCHH_MSGENULLP);
  ASSERT ( params, status, STDBURSTSEARCHH_ENULLP, STDBURSTSEARCHH_MSGENULLP);
  /* parse parameters */
  params = params->next;

  while(params) {
    SetStringParameter(buffer);
    
    argc++;
    
    argv = (CHAR **)LALRealloc(argv,argc*sizeof(CHAR *));
    if(!argv) {ABORT(status, STDBURSTSEARCHH_EMEM, STDBURSTSEARCHH_MSGEMEM);}

    argv[argc-1] = (CHAR *)LALCalloc(1+strlen(buffer),sizeof(CHAR));
    if(!(argv[argc-1])) {ABORT(status, STDBURSTSEARCHH_EMEM, STDBURSTSEARCHH_MSGEMEM);}

    strcpy(argv[argc-1],buffer);
  }

  /* pack input in one time series */
  tmp.epoch = input->epoch;
  tmp.deltaT = input->deltaT;
  tmp.f0 = input->f0;
  tmp.data = &tv;
  tv.length = input->data->vectorLength;
  tv.data = input->data->data;
  
  /* run EP ETG */
  EPInitSearch(status->statusPtr, (void *)(&EPparams), argv, argc);
  CHECKSTATUSPTR (status);

  EPConditionData(status->statusPtr, &tmp, EPparams);
  CHECKSTATUSPTR (status);
  
  EPSearch(status->statusPtr, EPparams, &burstEvent, EPparams->initParams->segDutyCycle);
  CHECKSTATUSPTR (status);
  
  EPFinalizeSearch(status->statusPtr, (void *)(&EPparams));
  CHECKSTATUSPTR (status);

  /* generate output */
  bzero(output, sizeof(EventIDColumn));

  while(burstEvent) {

    output->next = (EventIDColumn *)LALCalloc(1,sizeof(EventIDColumn));
    if(!(output->next)) {ABORT(status, STDBURSTSEARCHH_EMEM, STDBURSTSEARCHH_MSGEMEM);}
    
    output = output->next;

    output->snglBurstTable = burstEvent;

    burstEvent = burstEvent->next;

  }

  /* clean up */
  for(i=0;i<argc;i++) {
    LALFree(argv[i]);
  }
  LALFree(argv);

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}
