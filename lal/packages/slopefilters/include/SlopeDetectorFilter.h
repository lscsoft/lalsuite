/*************************** <lalVerbatim file="SlopeDetectorFilterHV">
Author: Daw, E. J.
$Id$
**************************** </lalVerbatim> */

/*************************** <lalLaTeX>

\section{Header \texttt{SlopeDetectorFilter.h}}

[One sentence briefly defining scope of the header]

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/SlopeDetectorFilter.h>
\end{verbatim}

[Generic documentation on the header; this is the main place to
document any stuff not specific to the module]

\subsection*{Error conditions}
\input{SlopeDetectorFilterHE}

\subsection*{Structures}

[Document here any structures defined in the header.  
Also include any of them in the index; e.g.:]
% \index{\texttt{SlopeDetectorFilterOutput}}
% \index{\texttt{SlopeDetectorFilterInput}}
% \index{\texttt{SlopeDetectorFilterParams}}

\vfill{\footnotesize\input{SlopeDetectorFilterHV}}
\newpage\input{SlopeDetectorFilterC}
\newpage\input{SlopeDetectorFilterTestC}

********************************** </lalLaTeX> */

#ifndef _SLOPEDETECTORFILTER_H
#define _SLOPEDETECTORFILTER_H

#include <lal/LALStdlib.h>
/******* INCLUDE ANY OTHER LAL HEADERS needed for header (NOT module) ****/

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (SLOPEDETECTORFILTERH, "$Id$");

/******************************** <lalErrTable file="SlopeDetectorFilterHE"> */

#define SLOPEDETECTORFILTERH_EINPUTNULLP        1
#define SLOPEDETECTORFILTERH_EOUTPUTNULLP       2
#define SLOPEDETECTORFILTERH_EDATATOOSHORT      3

#define SLOPEDETECTORFILTERH_MSGEINPUTNULLP   "Null input pointer"
#define SLOPEDETECTORFILTERH_MSGEOUTPUTNULLP  "Null output pointer"
#define SLOPEDETECTORFILTERH_MSGEDATATOOSHORT "Data segment too short"

/************************************ </lalErrTable> */

/****** DEFINE OTHER GLOBAL CONSTANTS OR MACROS ************/

/****** DEFINE NEW STRUCTURES AND TYPES ************/

/****** INCLUDE EXTERNAL GLOBAL VARIABLES ************/

void
LALSlopeDetectorFilter( LALStatus            *status,
			REAL4Vector*         output_data,
			const REAL4Vector*   input_data,
			const UINT4          ntaps );

#ifdef  __cplusplus
}
#endif

#endif /* _SLOPEDETECTORFILTER_H */



