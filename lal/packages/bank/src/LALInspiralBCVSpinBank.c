/*  <lalVerbatim file="LALInspiralBCVSpinBankCV">
Author: Tagoshi, H, Van Den Broeck, C, Jones, G., Sathyaprakash, B.S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>
\subsection{Module \texttt{LALInspiralBCVSpinBank.c}}
\subsubsection*{Prototypes}
\subsubsection*{Description}
\subsubsection*{Algorithm}
\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}
\subsubsection*{Notes}
\clearpage
</lalLaTeX>  */


#include <math.h>
#include <lal/AVFactories.h>
#include <lal/FlatMesh.h>
#include <lal/LALConfig.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALMalloc.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALStdlib.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/MatrixUtils.h>
#include <lal/SeqFactories.h>
#include <lal/LALInspiralBCVSpinBank.h>

int BCVspin_beta_placement (
		double beta_min,
		double beta_max,
		int n,
		double *Sn,
		double fmin,
		double fmax,
		double *beta_list,
		int *nbeta);

/* <lalVerbatim file="LALInspiralBCVSpinBankCP"> */
void
LALInspiralBCVSpinBank(
    LALStatus         	 *status,
    SnglInspiralTable   **tiles,
    INT4      		 *ntiles,
    InspiralCoarseBankIn *coarseIn
    )
/* </lalVerbatim> */

{

	/* coarseIn structure knows about the range of (1)beta, (2)psi0,
	 * and (3) psi3. It also provides the power spectrum of noise,
	 * minimal match (in mmCoarse), 
	 */
	double betaMin, betaMax, *dpsi0Min, *psi0Max, *beta_list;
	int n, nbeta; 
	
	betaMin = coarseIn->betaMin;
	betaMax = coarseIn->betaMax;
	BCVspin_beta_placement (beta_min, beta_max, n, Sn, fmin, fmax, &beta_list, &nbeta);
	beta=betamin; 
	for (i=0; i<beta; i++)
	{
		beta = beta_list[i];
		BCVspin_effmetric(&bcvspinmetric, *&a, &effmetric);
	}
}
