/************************************ <lalVerbatim file="LALMomentCV">
Author: Tibbits, M M
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{LALMoment.c}}
\label{s:LALMoment.c}

Routine to compute various moments of data.

\subsubsection*{Prototypes}
\input{LALMomentCP}

\subsubsection*{Description}
The data structure passed in is a REAL8 Sequence.  The only parameter is which moment to calculate.
The function the sums the data, calculates the average and then uses the algorithm below to find the moment
that was asked for.

\subsubsection*{Algorithm}
\begin{itemize}
\item \textit{Find the mean (here referred to as $ \overline{x} $).}
\item \textit{Sum, over all the elements, the quantity: \((x[k] - \overline{x})^{n}. \)}
\item \textit{Divide the sum just made by N-1. Call it moment-n}
\item \textit{If n is greater than 2:}
\begin{itemize}
\item \textit{Sum, over all the elements, the quantity \((x[k]-xbar)^{2}\)}
\item \textit{Divide the sum just made by N-1.}
\item \textit{Take the square root of this quantity just made and call it sigma.}
\item \textit{Divide moment-n by sigma$^{n}$.}
\end{itemize}

\item \textit{Return moment-n}
\end{itemize}

\subsubsection*{Uses}

Determination of a specific moment of a set of data.

\subsubsection*{Notes}

\begin{itemize}
\item \textit{Moments less than two are not allowed.}
\item \textit{The result structure must be Non-NULL when passed in.}
\item \textit{The function assumes that the length member of the data passed in is correct.}
\end{itemize}

\vfill{\footnotesize\input{LALMomentCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <lal/LALMoment.h>


NRCSID( LALMOMENTC, "$Id$");


/* <lalVerbatim file="LALMomentCP"> */
void LALMoment(
	LALStatus		*status,
	REAL8			*result,
	REAL8Sequence		*data,
	INT4			*whichMoment
)
{	/* </lalVerbatim> */

	/*  Turns error messages on/off */
	/* OMITTED -- JC
         * extern	INT4	lalDebugLevel;
         */


	/*  Constant for Moments Higher than two.  */
	const	INT2	two = 2;


	/*  Variable Declarations  */
	INT4	momentToCalculate;
	INT4	iterator;
	INT4	length;
	/* REAL8	temp	= 0.0; */
	REAL8	ave 	= 0.0;
	REAL8	momentn	= 0.0;
	REAL8	sigma	= 0.0;
	REAL8	sum	= 0.0;
	REAL8	base	= 0.0;

	INITSTATUS( status, "LALMoment", LALMOMENTC);

	/*  Check input for existence.  */
	/*  Result should come in Allocated  */
	ASSERT ( result, status, LALMOMENTH_ENULL, LALMOMENTH_MSGENULL);

	/*  whichMoment should be defined  */
	ASSERT ( whichMoment, status, LALMOMENTH_ENULL, LALMOMENTH_MSGENULL);

	/*  whichMoment must be greater than 1  */
	ASSERT ( *whichMoment > 1, status, LALMOMENTH_ENULL, LALMOMENTH_MSGENULL);

	/*  data must be defined  */
	ASSERT ( data, status, LALMOMENTH_ENULL, LALMOMENTH_MSGENULL);

	length = ((REAL8Sequence*)(data))->length;
	momentToCalculate = (*whichMoment);

	/*  length must be greater than one  */
	ASSERT ( length > 1, status, LALMOMENTH_ELNTH, LALMOMENTH_MSGELNTH);


	for (iterator = 0; iterator < length; iterator++)
	{
		sum += data->data[iterator];
	}

	ave = ( sum / length );


	/*  Only print debug info if told to do so  */
        /* OMITTED -- JC
	 * if (lalDebugLevel == 3 )
	 * {
         *   for ( iterator = 0; iterator < length; iterator++ )
         *   {
         *     printf("data[%d] = %f \n", iterator, data->data[iterator]);
         *   }
         *
         *   printf("avg     := %f\n",ave);
         *   printf("length	:= %d\n",length);
	 * }
         */

	/*  Return the Mean if whichMoment = 1  */
	if ( momentToCalculate == 1 )
	{
		*result = ave;
		RETURN(status);
	}


	for (iterator=0; iterator < length; iterator++)
	{
		base = (data->data[iterator] - ave);
		sum = pow( base, momentToCalculate );
		momentn += sum;
	}

	momentn /= ((REAL8)(length - 1));


	if ( momentToCalculate > 2 )
	{
		for (iterator = 0; iterator < length; iterator++)
		{
			base = (data->data[iterator] - ave);
			sum = pow( base, two );
			sigma += sum;
		}

		sum = sigma/(length - 1);
		sum = sqrt(sum);
		sum = pow( sum, momentToCalculate );

		if( sigma != 0.0)
		{
			momentn /= sigma;
		}
	}


	*result = momentn;

        /* OMITTED -- JC
	 * if ( lalDebugLevel == 3 )
	 * {
         *   printf("MESG:  momentn	:= %f\n",*result);
         *   printf("MESG:  LALMoment() complete. \n");
         * }
         */

	RETURN (status);
}
