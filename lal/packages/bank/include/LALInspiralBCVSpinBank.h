/* <lalVerbatim file="LALInspiralBCVSpinBankHV">

Author: Tagoshi, H., Van Den Broeck, C., Jones, G.
$Id$

</lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LALInspiralBCVSpinBank.h}}
\label{s:LALInspiralBCVSpinBank.h}

Header file for the template placement codes.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALInspiralBCVSpinBank.h>
\end{verbatim}

\noindent This header file covers routines that are used in template placement
for spinning black hole binaries using phenomenological BCV templates.

</lalLaTeX> */

#include <lal/LIGOMetadataTables.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>


void
LALInspiralBCVSpinBank(
    LALStatus         	 *status,
    SnglInspiralTable   **tiles,
    INT4      		 *ntiles,
    InspiralCoarseBankIn *coarseIn
    );

int BCVspin_metric(
   /*input*/
   int N,double *Sn,double fmin,double fmax,double beta,
   /*output*/
   double **bcv2metric,int dbg);

int BCVspin_spacing(
    double **metric3,
    double **a,
    double *deltax);

int BCVspin_effmetric(
    /* input */
    double **metric3,
    double **a,
    /* output */
    double **effmetric);
