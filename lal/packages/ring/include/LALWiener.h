/************************************ <lalVerbatim file="LALWienerHV">
Author: Dwyer, Steven J.
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\section{Header \texttt{LALWiener.h}}
\label{s:LALWiener.h}

This Package passes input through a wiener filter

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALWiener.h>
\end{verbatim}

\noindent 

*********************************************************** </lalLaTeX> */

#ifndef _LALWIENER_H  /* Double-include protection. */
#define _LALWIENER_H

#include <lal/LALDatatypes.h>
#include <lal/RealFFT.h>
#include <math.h>
#include <lal/LALMalloc.h>
#include <lal/AVFactories.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

NRCSID( LALWienerH, "$Id$" );

/**********************************************************  <lalLaTeX>

\subsection*{Error codes}
\input{LALWienerHE}

\vfill{\footnotesize\input{LALWienerHV}}

*********************************************************** </lalLaTeX> */

/********************************************************** <lalErrTable file="LALWienerHE" > */

#define LALWIENERH_ENULL 1
#define LALWIENERH_ESIZE 2
#define LALWIENERH_EDATA 3
#define LALWIENERH_EMEM  4

#define LALWIENERH_MSGENULL "Arguments contained an unexpected null pointer"
#define LALWIENERH_MSGESIZE "Template is smaller than source vector"
#define LALWIENERH_MSGEDATA "Member Data Not Allocated"
#define LALWIENERH_MSGEMEM  "Memory Allocation Error"

/********************************************************** </lalErrTable> */

/* <lalLaTeX>

\section*{Structures}
\input{LALWienerHS}
</lalLaTeX>  */





/* <lalVerbatim file="LALWienerHS">  */
typedef struct
tagWienerOutput
{
  REAL4Vector     *q;
  REAL4Vector     *z;
}
WienerOutput;

/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{WienerOutput}}
</lalLaTeX>  */




/* <lalVerbatim file="LALWienerHS">  */
typedef struct
tagWienerUnFormattedInput
{
  REAL4Vector     *h;
  REAL4Vector     *s;
}
WienerUnFormattedInput;

/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{WienerUnFormattedOutput}}
</lalLaTeX>  */





/* <lalVerbatim file="LALWienerHS">  */
typedef struct
tagWienerFormattedInput
{
  REAL4Vector        *h;
  COMPLEX8Vector     *st;
  INT2       *TotalLength; 
}
WienerFormattedInput;

/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{WienerFormattedOutput}}
</lalLaTeX>  */

/*  <lalLaTeX>
\vfill{\footnotesize\input{LALWienerHV}}
</lalLaTeX>  */




/********************************************************** <lalLaTeX>
\newpage\input{LALWienerC}
******************************************************* </lalLaTeX> */

/* Function prototypes */

void
LALUnFormattedWiener( 
     LALStatus *status, 
     WienerOutput *output, 
     WienerUnFormattedInput *input 
     );

void
LALFormattedWiener( 
     LALStatus *status, 
     WienerOutput *output, 
     WienerFormattedInput *input 
     );

/********************************************************** <lalLaTeX>
\newpage\input{LALWienerTestC}
******************************************************* </lalLaTeX> */

#ifdef  __cplusplus
}
#endif  /* C++ protection. */
#endif  /* Double-include protection. */











