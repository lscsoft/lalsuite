/* <lalVerbatim file="LALMomentHV">

Author: Tibbits, M. M.
$Id$

</lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LALMoment.h}}
\label{s:LALMoment.h}

\begin{verbatim}
The LALDMoment() and LALSMoment() associated header file.
(S - single precision )
(D - double precision )
\end{verbatim}

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALMoment.h>
\end{verbatim}

\noindent This header provides the prototype for the LALDMoment() and LALSMoment() function.

 </lalLaTeX> */

/* Double Include Protection */
#ifndef _LALMOMENT_H
#define _LALMOMENT_H

#include <lal/LALStdlib.h>
#include <lal/LALDatatypes.h>


/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif


NRCSID( LALMOMENTH, "$Id$");


/*  <lalLaTeX>

\subsection*{Error codes}
\input{LALMomentHE}

\vfill{\footnotesize\input{LALMomentHV}}
</lalLaTeX>  */

/* <lalErrTable file="LALMomentHE"> */

#define LALMOMENTH_ENULL 1
#define	LALMOMENTH_ENNUL 2
#define LALMOMENTH_ELNTH 3
#define LALMOMENTH_ESEGZ 4
#define LALMOMENTH_ENUMZ 5
#define LALMOMENTH_EALOC 6

#define LALMOMENTH_MSGENULL "NULL pointer."
#define	LALMOMENTH_MSGENNUL "Non-NULL pointer."
#define LALMOMENTH_MSGELNTH "Must have more than one data point."
#define LALMOMENTH_MSGESEGZ "Invalid number of segments"
#define LALMOMENTH_MSGENUMZ "Invalid number of points in segment"
#define LALMOMENTH_MSGEALOC "Memory Allocation Error"

/* </lalErrTable> */


/* Function prototypes */


/*  <lalLaTeX>
\newpage\input{LALMomentC}
</lalLaTeX>  */

void LALSMoment 
(
	LALStatus		*status,
	REAL4			*result,
	REAL4Sequence		*data,
	INT4			whichMoment
);


void LALDMoment 
(
	LALStatus		*status,
	REAL8			*result,
	REAL8Sequence		*data,
	INT4			whichMoment
);

/*  <lalLaTeX>
\newpage\input{LALMomentTestC}
</lalLaTeX>  */

#ifdef  __cplusplus
}
#endif  /* C++ protection. */
#endif  /* Double-include protection. _LALMOMENT_H_ */
