#include "Inspiral.h"
#include "LALStdlib.h"

NRCSID (RUNGEKUTTA4C, "$Id$");

void RungeKutta4(Status *status,
	         REAL8Vector *yout,
	         rk4In *input,
	         void *params)
{
	INT4 i;
	REAL8 xh,hh,h6;

   INITSTATUS(status, "RungeKutta4", RUNGEKUTTA4C);
   ATTATCHSTATUSPTR(status);

   ASSERT (yout, status, RUNGEKUTTA4_ENULL, RUNGEKUTTA4_MSGENULL);
   ASSERT (yout->data, status, RUNGEKUTTA4_ENULL, RUNGEKUTTA4_MSGENULL);
   ASSERT (input, status, RUNGEKUTTA4_ENULL, RUNGEKUTTA4_MSGENULL);
   ASSERT (params, status, RUNGEKUTTA4_ENULL, RUNGEKUTTA4_MSGENULL);


	hh = input->h*0.5;
	h6 = input->h/6.0;
	xh = input->x+hh;


	for (i=0;i<=input->n-1;i++) input->yt->data[i]=input->y->data[i]+hh*input->dydx->data[i];

	input->function(status->statusPtr, input->yt, input->dyt, params);
        CHECKSTATUSPTR(status);

 
	for (i=0;i<=input->n-1;i++) input->yt->data[i]=input->y->data[i]+hh*input->dyt->data[i];

	input->function(status->statusPtr, input->yt, input->dym, params);	
	CHECKSTATUSPTR(status);


	for (i=0;i<=input->n-1;i++) {
		input->yt->data[i]=input->y->data[i]+input->h*input->dym->data[i];
		input->dym->data[i] += input->dyt->data[i];
	}
	input->function(status->statusPtr, input->yt, input->dyt, params);
	CHECKSTATUSPTR(status);

	for (i=0;i<=input->n-1;i++)
		yout->data[i]=input->y->data[i]+h6*(input->dydx->data[i]+input->dyt->data[i]+2.0*input->dym->data[i]);





   DETATCHSTATUSPTR(status);
   RETURN (status);




}
