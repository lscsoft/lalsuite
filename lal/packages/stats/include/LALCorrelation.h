/********************************* <lalVerbatim file="LALCorrelationHV">
Author: Yakushin, Igor
$Id$
********************************** </lalVerbatim> */

/********************************* <lalLaTeX>

\section{Header \texttt{LALCorrelation.h}}

[One sentence briefly defining scope of the header]

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALCorrelation.h>
\end{verbatim}

[Generic documentation on the header; this is the main place to
document any stuff not specific to the module]

\subsection*{Error conditions}
\input{LALCorrelationHE}

\subsection*{Structures}

[Document here any structures defined in the header.  
Also include any of them in the index; e.g.:]
% \index{\texttt{LALCorrelationOutput}}
% \index{\texttt{LALCorrelationInput}}
% \index{\texttt{LALCorrelationParams}}

\vfill{\footnotesize\input{LALCorrelationHV}}
\newpage\input{LALCorrelationC}
\newpage\input{LALCorrelationTestC}

********************************** </lalLaTeX> */

#ifndef _LALCORRELATION_H
#define _LALCORRELATION_H

#include <lal/LALStdlib.h>
/******* INCLUDE ANY OTHER LAL HEADERS needed for header (NOT module) ****/

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (LALCORRELATIONH, "$Id$");

/******************************** <lalErrTable file="LALCorrelationHE"> */

#define LALCORRELATIONH_ENULLP        1
#define LALCORRELATIONH_ESTART        2
#define LALCORRELATIONH_ESAMPLING     3

#define LALCORRELATIONH_MSGENULLP     "Null pointer"
#define LALCORRELATIONH_MSGESTART    "Time series do not start simultaneously" 
#define LALCORRELATIONH_MSGESAMPLING "Time series are not sampled with the same rate"

/* *********************************** </lalErrTable> */

/****** DEFINE OTHER GLOBAL CONSTANTS OR MACROS ************/

/****** DEFINE NEW STRUCTURES AND TYPES ************/

typedef struct 
tagCorrelationParams
{
  REAL4 maxTimeShiftNan;
}
CorrelationParams;

typedef struct
tagInputCorrelation
{
  REAL4TimeSeries *one;
  REAL4TimeSeries *two;
}
InputCorrelation;

typedef struct
tagOutputCorrelation
{
  REAL4 *timeShiftedCorrelation; /* cor(x(t-T),y(t)) or cor(x(t),y(t-T)) as a function of T*/

  INT4 maxCorrelationTimeShift;  /* time step (from shift=0) at which correlation is maximum */
  REAL4 maxCorrelationValue;     /* maximum value of the correlation for the range of shifts */

  INT4 minCorrelationTimeShift;
  REAL4 minCorrelationValue;

  LIGOTimeGPS start; /* start time for the time series under consideration */
  INT4 length; /* number of time steps in the time series */
  REAL8 deltaT; /* time separation in seconds between different time steps */

  INT4 shift; /* max time step shift between time series */
}
OutputCorrelation;

/****** INCLUDE EXTERNAL GLOBAL VARIABLES ************/

void
LALCorrelation( LALStatus                      *status,
		OutputCorrelation              **output,
		const InputCorrelation         *input,
		const CorrelationParams        *params);

#ifdef  __cplusplus
}
#endif

#endif /* _LDASCAMPMOMENT_H */
