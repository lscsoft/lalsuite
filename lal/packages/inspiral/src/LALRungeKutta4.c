/*  <lalVerbatim file="LALRungeKutta4CV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALRungeKutta4.c}}

The code \texttt{LALRungeKutta4.c} solves a system of $n$ coupled first--order differential equations.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALRungeKutta4CP}
\index{\verb&LALRungeKutta4()&}

\subsubsection*{Description}


\subsubsection*{Algorithm}


\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALRungeKutta4CV}}

</lalLaTeX>  */







#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>

NRCSID (LALRUNGEKUTTA4C, "$Id$");

/*  <lalVerbatim file="LALRungeKutta4CP"> */
void LALRungeKutta4(LALStatus *status,
	            REAL8Vector *yout,
	            rk4In *input,
	            void *params)
{ /* </lalVerbatim>  */

	INT4 i;
	REAL8 xh,hh,h6;

   INITSTATUS(status, "LALRungeKutta4", LALRUNGEKUTTA4C);


   ASSERT (yout, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (yout->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (input, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);


	hh = input->h*0.5;
	h6 = input->h/6.0;
	xh = input->x+hh;


	for (i=0;i<=input->n-1;i++) input->yt->data[i]=input->y->data[i]+hh*input->dydx->data[i];

	input->function(input->yt, input->dyt, params);


 
	for (i=0;i<=input->n-1;i++) input->yt->data[i]=input->y->data[i]+hh*input->dyt->data[i];

	input->function(input->yt, input->dym, params);	

	for (i=0;i<=input->n-1;i++) {
		input->yt->data[i]=input->y->data[i]+input->h*input->dym->data[i];
		input->dym->data[i] += input->dyt->data[i];
	}
	input->function(input->yt, input->dyt, params);


	for (i=0;i<=input->n-1;i++)
		yout->data[i]=input->y->data[i]+h6*(input->dydx->data[i]+input->dyt->data[i]+2.0*input->dym->data[i]);



   RETURN (status);




}
