/* <lalVerbatim file="LALRingDownHV">

Author: Tibbits, M. M.
$Id$

</lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LALRingDown.h}}
\label{s:LALRingDown.h}

\begin{verbatim}
The LALDRingDown() and LALSRingDown() associated header file.
(S - single precision )
(D - double precision )
\end{verbatim}

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALRingDown.h>
\end{verbatim}

\noindent This header provides the prototype for the LALDRingDown() and LALSRingDown() functions.

 </lalLaTeX> */

/* Double Include Protection */
#ifndef _LALRINGDOWN_H
#define _LALRINGDOWN_H

#include <lal/LALStdlib.h>
#include <lal/LALDatatypes.h>
#include <lal/LALConstants.h>

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif


NRCSID( LALRINGDOWNH, "$Id$");


/*  <lalLaTeX>

\subsection*{Error codes}
\input{LALRingDownHE}

\vfill{\footnotesize\input{LALRingDownHV}}
</lalLaTeX>  */

/* <lalErrTable file="LALRingDownHE"> */

#define LALRINGDOWNH_ENULL 1
#define	LALRINGDOWNH_ENNUL 2
#define LALRINGDOWNH_ELNTH 3
#define LALRINGDOWNH_ESEGZ 4
#define LALRINGDOWNH_ENUMZ 5
#define LALRINGDOWNH_EALOC 6
#define LALRINGDOWNH_EDATA 7
#define LALRINGDOWNH_ENEG  8
#define LALRINGDOWNH_ENPS  9
#define LALRINGDOWNH_EARG 10

#define LALRINGDOWNH_MSGENULL "NULL pointer."
#define	LALRINGDOWNH_MSGENNUL "Non-NULL pointer."
#define LALRINGDOWNH_MSGELNTH "Must have more than one data point."
#define LALRINGDOWNH_MSGESEGZ "Invalid number of segments"
#define LALRINGDOWNH_MSGENUMZ "Invalid number of points in segment"
#define LALRINGDOWNH_MSGEALOC "Memory Allocation Error"
#define LALRINGDOWNH_MSGEDATA "Output->data must be NULL"
#define LALRINGDOWNH_MSGENEG  "f0 must be >= 0 "
#define LALRINGDOWNH_MSGENPS  "n,df,and t must be positive"
#define LALRINGDOWNH_MSGEARG  "Improper argument"

/* </lalErrTable> */


/* <lalLaTeX>

\section*{Structures}
\input{LALRingDownHS}

</lalLaTeX>  */

/*  Structures  */
/* <lalVerbatim file="LALRingDownHS">  */
typedef struct
tagCOMPLEX8RingDownParams
{
	REAL4	f0;
	REAL4	df;
	REAL4	t;
	INT4	n;
}
COMPLEX8RingDownParams;
/* </lalVerbatim>  */

/* <lalVerbatim file="LALRingDownHS">  */
typedef struct
tagCOMPLEX16RingDownParams
{
	REAL8	f0;
	REAL8	df;
	REAL8	t;
	INT4	n;
}
COMPLEX16RingDownParams;
/* </lalVerbatim>  */


/* Function prototypes */

/*  <lalLaTeX>
\newpage\input{LALRingDownC}
</lalLaTeX>  */

void LALCOMPLEX8RingDown
(
	LALStatus		*status,
	COMPLEX8FrequencySeries	*result,
	COMPLEX8RingDownParams		*params
);

void LALCOMPLEX16RingDown
(
	LALStatus		*status,
	COMPLEX16FrequencySeries	*result,
	COMPLEX16RingDownParams		*params
);

/*  <lalLaTeX>
\newpage\input{LALRingDownTestC}
</lalLaTeX>  */

#ifdef  __cplusplus
}
#endif  /* C++ protection. */
#endif  /* Double-include protection. _LALRINGDOWN_H_ */
