/* <lalVerbatim file="InspiralBCVSpinBankTestCV">
 * Author: Van Den Broeck, C and Jones, G.
 * $Id$
 * </lalVerbatim> */ 

/* <lalLaTeX>
 * \subsection{Program \texttt{InspiralBCVSpinBankTest}}
 * \label{s:InspiralBCVSpinBankTest.c}
 * Tests InpsiralBCVSpinBankTest(). 
 *
 * \subsubsection*{Usage}
 *
 * \begin{verbatim}
 * InspiralBCVSpinBankTest
 * \end{verbatim} 
 *
 * This program uses InspiralBCVSpinBank() to generate a template bank from
 * command line parameter input.  
 *
 * \subsubsection{Command line options}
 * \begin{description}
 * \end{description}
 *
 * \subsubsection*{Exit codes}
 * \input{InspiralBCVSpinBankTestCE}
 *
 * \subsubsection*{Notes}
 *
 * \vfill{\footnotesize\input{InspiralBCVSpinBankTestCV}}
 * </lalLaTeX> */

/* <lalErrTable file="InspiralBCVSpinBankTestCE"> */
#define INSPIRALBCVSPINBANKTESTC_ENORM     0
#define INSPIRALBCVSPINBANKTESTC_EMEM      1
#define INSPIRALBCVSPINBANKTESTC_ESUB      2
#define INSPIRALBCVSPINBANKTESTC_EFILE     4

#define INSPIRALBCVSPINBANKTESTC_MSGENORM  "Normal exit"
#define INSPIRALBCVSPINBANKTESTC_MSGEMEM   "Memory allocation error"
#define INSPIRALBCVSPINBANKTESTC_MSGESUB   "Subroutine error"
#define INSPIRALBCVSPINBANKTESTC_MSGEFILE  "File I/O error"
/* </lalErrTable> */

int lalDebugLevel=33;
#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALGSL.h>
#include <gsl/gsl_matrix.h>

void testfunction(LALStatus *status);

int main (void)
{
  static LALStatus status;
  testfunction(&status);
  return 0;
}

void testfunction(LALStatus *status)
{
  int i, j; 
  float x;
  gsl_matrix *m;

  ATTATCHSTATUSPTR( status);
  CALLGSL((m = gsl_matrix_alloc (10, 3)), status);
  for (i = 0; i < 10; i++)
    for (j = 0; j < 3; j++)
      CALLGSL ( gsl_matrix_set (m, i, j, 0.23 + 100*i + j), status);
  
  for (i = 0; i < 10; i++)
    for (j = 0; j < 3; j++)
{
      CALLGSL((x = gsl_matrix_get (m, i, j)), status);
      printf ("m(%d,%d) = %g\n", i, j, x);
}

  CHECKSTATUSPTR (status);
  REPORTSTATUS (status);
  DETATCHSTATUSPTR( status);
  RETURN (status);
}


