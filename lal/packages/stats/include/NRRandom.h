/* <lalVerbatim file="NRRandomHV">

Author: Tibbits, M. M.
$Id$

</lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{NRRandom.h}}
\label{s:NRRandom.h}

\begin{verbatim}
The Numerical Recipes Random Functions associated header file.
\end{verbatim}

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/NRRandom.h>
\end{verbatim}

\noindent This header provides the prototype for the Numerical Recipes
 Random Functions as listed below.

 </lalLaTeX> */

/* Double Include Protection */
#ifndef _NRRANDOM_H
#define _NRRANDOM_H

#include <lal/LALStdlib.h>
#include <lal/LALDatatypes.h>
#include <lal/LALConstants.h>

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif


NRCSID( NRRANDOMH, "$Id$");


/*  <lalLaTeX>

\subsection*{Error codes}
\input{NRRandomHE}

</lalLaTeX>  */

/* <lalErrTable file="NRRandomHE"> */

#define NRRANDOMH_ENULL 1
#define	NRRANDOMH_ENNUL 2
#define NRRANDOMH_ELNTH 3
#define NRRANDOMH_ESEGZ 4
#define NRRANDOMH_ENUMZ 5
#define NRRANDOMH_EALOC 6
#define NRRANDOMH_EDATA 7
#define NRRANDOMH_EARG  8
#define NRRANDOMH_ETEST 9
#define NRRANDOMH_ESEED 10

#define NRRANDOMH_MSGENULL "NULL pointer"
#define	NRRANDOMH_MSGENNUL "Non-NULL pointer"
#define NRRANDOMH_MSGELNTH "Must have more than one data point"
#define NRRANDOMH_MSGESEGZ "Invalid number of segments"
#define NRRANDOMH_MSGENUMZ "Invalid number of points in segment"
#define NRRANDOMH_MSGEALOC "Memory Allocation Error"
#define NRRANDOMH_MSGEDATA "Output->data must be NULL"
#define NRRANDOMH_MSGEARG  "Improper argument"
#define NRRANDOMH_MSGETEST "Improper value for random test number"
#define NRRANDOMH_MSGESEED "Improper seed value"

/* </lalErrTable> */

/*  Constants  */
/* <lalVerbatim file="NRRandomHConstants">  */
#define IA     16807
#define IM     2147483647
#define AM     1.0/IM
#define IQ     127773
#define IR     2836
#define MASK   123459876
#define NTAB   32
#define NDIV   (1.0 + (IM - 1.0) / NTAB)
#define EPS    1.2e-7
#define RNMX   (1.0 - EPS)
#define IM1    2147483563
#define IM2    2147483399
#define AM1    (1.0 / IM1)
#define IMM1   (IM1 - 1)
#define IA1    40014
#define IA2    40692
#define IQ1    53668
#define IQ2    52774
#define IR1    12211
#define IR2    3791
#define NDIV1  (1.0+((REAL8)((IMM1)/NTAB))
#define MBIG   1000000000
#define MSEED  161803398
#define MZ     0
#define FAC    (1.0/MBIG)
#define NITER  4

/*  </lalVerbatim>
 <lalLaTeX>
\section*{Structures}
\input{NRRandomHS}

\vfill{\footnotesize\input{NRRandomHV}}
</lalLaTeX>  */

/*  Structures  */
/* <lalVerbatim file="NRRandomHS">  */
typedef float (*RanFunction)
(
	long	*idum
);
/* </lalVerbatim>  */

enum RandNumGen { nrran3, nrran4 };

/* Function prototypes */

/*  <lalLaTeX>
\newpage\input{NRRandomC}
</lalLaTeX>  */

float ran3(long *idum);

void psdes(unsigned long* lword, unsigned long *irword);

float ran4(long *idum);

float gasdev(long *idum, RanFunction gen);

void NRRan3
(
	LALStatus	*status,
	REAL4		*output,
	INT4		*seed
);

void NRRan4
(
	LALStatus	*status,
	REAL4		*output,
	INT4		*seed
);

void NRGasDev
(
	LALStatus	*status,
	REAL4		*output,
	INT4		*seed,
	INT2		*param
);

/*  <lalLaTeX>
\newpage\input{NRRandomTestC}
</lalLaTeX>  */

#ifdef  __cplusplus
}
#endif  /* C++ protection. */
#endif  /* Double-include protection. _NRRANDOM_H_ */
