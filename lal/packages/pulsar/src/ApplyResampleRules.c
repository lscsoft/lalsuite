/*************************** <lalVerbatim file="ApplyResampleRulesCV">
Author: Creighton, T. D.
Revision: $Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{ApplyResampleRules.c}}
\label{ss:ApplyResampleRules.c}

Resamples a time series according to a set of resampling rules.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{ApplyResampleRulesCP}
\index{\texttt{LALApplyResampleRules()}}

\subsubsection*{Description}

This function sets \verb@output->deltaT@ and fills \verb@output->data@
with data from \verb@*input@, using the resampling rules specified in
\verb@*rules@.  If the timespan required to fill \verb@output->data@
is not a subset of the timespan covered by \verb@*input@ or
\verb@*rules@, the data at the nonintersecting times are set to zero.

\subsubsection*{Algorithm}

At present this routine is just a stub.  It does not apply or even
check \verb@*rules@, and instead simply makes \verb@*output@
equivalent to (a subset of) \verb@*input@.

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{ApplyResampleRulesCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/Resample.h>

NRCSID(APPLYRESAMPLERULESC,"$Id$");

/* <lalVerbatim file="ApplyResampleRulesCP"> */
void
LALApplyResampleRules( LALStatus       *stat,
		       REAL4TimeSeries *output,
		       REAL4TimeSeries *input,
		       ResampleRules   *rules )
{ /* </lalVerbatim> */
  INT4 nStart, nStop; /* output domain for which we can get data */

  INITSTATUS(stat,"LALApplyResampleRules",APPLYRESAMPLERULESC);

  /* Check that the inputs all exist. */
  ASSERT(rules,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(output,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(output->data,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(output->data->data,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(input,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(input->data,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(input->data->data,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);

  /* Set the output sampling time. */
  output->deltaT=input->deltaT;

  /* Find the difference between the input and output start times, in
     samples. */
  nStart=(INT4)((input->epoch.gpsSeconds-output->epoch.gpsSeconds)
		/output->deltaT);
  nStart+=(INT4)((input->epoch.gpsNanoSeconds
		  -output->epoch.gpsNanoSeconds)
		 /(1e9*output->deltaT));
  if(nStart>(INT4)(output->data->length))
    nStart=output->data->length;

  /* Ditto for stop times. */
  nStop=nStart+input->data->length; /* since deltaT's are equal */
  if(nStop>(INT4)(output->data->length))
    nStop=output->data->length;

  if(nStart>0)
    memset(output->data->data,0,nStart*sizeof(REAL4));
  if(nStop>nStart)
    memcpy(output->data->data+nStart,input->data->data,
	   (nStop-nStart)*sizeof(REAL4));
  if((INT4)(output->data->length)>nStop)
    memset(output->data->data+nStop,0,
	   (output->data->length-nStop)*sizeof(REAL4));

  /* That's all for the current stub. */
  RETURN(stat);
}
