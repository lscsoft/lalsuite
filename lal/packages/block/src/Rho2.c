/***************************** <lalVerbatim file="BlockRho2CV">
Author: Tibbits, Matthew M.
$Id$
 **************************************************** </lalVerbatim> */

/******************************* <lalLaTeX file="BlockRho2C">
\subsection{Module \texttt{BlockRho2.c}}
\label{s:BlockRho2.c}

\subsubsection*{Prototypes}
\input{Rho2CP}
\index{\texttt{BlockRho2()}}

\subsection*{Description}

To be written

\subsection*{Algorithms}

Matlab Code:

does the following:
%
\begin{enumerate}
\item 

\item 

\item 

\item
\end{enumerate}

\subsection*{Uses}
 
\begin{itemize}
\item \texttt{LALDDestroyVector}
\end{itemize}

\subsection*{Notes}

\vfill{\footnotesize\input{BlockRho2CV}}

 *************************************************** </lalLaTeX> */

#include <lal/BlockRho.h>

NRCSID( BLOCKRHO2, "$Id$");

/************************ <lalVerbatim file="Rho2CP"> */
static REAL8	nextOrder( INT4 arg )
/*********************************************** </lalVerbatim> */
{
	INT4	rem;
	REAL8	temp;

	temp = 0;

	if ( arg == 0 || arg == 1 )
		return ( 0 );

	rem = arg % 2;

	if ( rem == 0 )
	{
		temp = (arg / 2) * (( log(arg)) - 1 ) + ( log( 2 )/2 ) + log( LAL_PI / 2 ) / 2 + log(( 12 * arg ) + 1 ) - log (( 12 * arg ) + 2 );
		return (temp);
	}
	else
	{
		temp = (( arg - 1) / 2)*( log( arg - 1 ) -1) +  log( LAL_PI * ( arg - 1 )) / 2 + log(( 6 * arg ) - 5 ) - log(( 6 * arg ) - 6 );
		return (temp);
	}
}

/************************ <lalVerbatim file="Rho2CP"> */
static REAL8  iFunct( INT4 arg )
/*********************************************** </lalVerbatim> */
{
	REAL8	sum;
        INT4    rem;
	INT4	iterator;

        rem = arg % 2;

	sum = 0;

        if ( rem == 0 )
	{
		for ( iterator = (arg - 1); iterator >= 1; )
		{
			sum += log ( iterator );
			iterator -= 2;
		}
		sum += log ( LAL_PI/2 ) / 2;
	}
        else
	{
		for ( iterator = (arg - 1); iterator >= 2; )
		{
			sum += log ( iterator );
			iterator -= 2;
		}
	}

	return ( sum );
}

/************************ <lalVerbatim file="Rho2CP"> */
static void	lnIRatio(
	LALStatus	*status,
	REAL8Vector	**result,
	INT4		N,
	INT4		T
)
/*********************************************** </lalVerbatim> */
{
	REAL8Vector	*INk2;
	REAL8Vector	*Ik2;
	REAL8		IN2;
	INT4		iterator;

        INITSTATUS( status, "lnIRatio" , BLOCKRHO2);
        ATTATCHSTATUSPTR( status );

	INk2 = Ik2 = NULL;

	LALDCreateVector( status->statusPtr, &INk2, N - 3);
	LALDCreateVector( status->statusPtr, &Ik2,  N - 3);
	LALDCreateVector( status->statusPtr, result, N - 3);
	CHECKSTATUSPTR( status );

	for( iterator = 0; iterator < N - 3; iterator++)
		INk2->data[iterator] = Ik2->data[iterator] = 0;

	IN2 = nextOrder ( N - 2 );

	for( iterator = 2; iterator <= (N - 2); iterator++ )
	{
		if(( iterator - 2 ) < T )
			Ik2->data[ iterator - 2 ] = iFunct( iterator - 2 );
		else
			Ik2->data[ iterator - 2 ] = nextOrder( iterator - 2 );

		if(( N - iterator - 2 ) < T )
			INk2->data[ iterator - 2 ] = iFunct( N - iterator - 2 );
		else
			INk2->data[ iterator - 2 ] = nextOrder( N - iterator - 2 );
	}

	for( iterator = 0; iterator < (INT4)(*result)->length; iterator++)
		(*result)->data[iterator] = Ik2->data[iterator] + INk2->data[iterator] - IN2;

	LALDDestroyVector( status->statusPtr, &INk2);
	LALDDestroyVector( status->statusPtr, &Ik2);

        DETATCHSTATUSPTR( status );
        RETURN (status);
}

/************************ <lalVerbatim file="Rho2CP"> */
void LALBlockRho2 (
	LALStatus		*status,
	REAL8			*result,
	REAL8			*rpeak,
	INT4			*myindex,
	REAL8Sequence		*data,
	UINT4			*marginOfExclusion			
)
/*********************************************** </lalVerbatim> */
{
	/*  Variable Declarations  */
	REAL8Vector	*indexVector;
	REAL8Vector	*xbar0k;
	REAL8Vector	*x2bar0k;
        REAL8Vector     *xbarkNk;
        REAL8Vector     *x2barkNk;
	REAL8Vector	*temp;
        REAL8Vector     *temp2;
	REAL8Vector	*Y0k;
        REAL8Vector     *YkNk;
	REAL8Vector	*fom;
	REAL8Vector	*fomWithoutMargins;
	REAL8		Y0N;
	REAL8		tempReal;
	REAL8		xbar0N;
	REAL8		x2bar0N;
	REAL8		rootRatio;
const	INT4		two = 2;
	INT4		k;
        INT4		iterator;
        INT4		N;

	INITSTATUS( status, "LALBlockRho2" , BLOCKRHO2);
	ATTATCHSTATUSPTR( status );

	/*  Check input for existence.  */
	/*  data must be defined  */
	ASSERT ( data, status, LALMOMENTH_ENULL, LALMOMENTH_MSGENULL);

	/*  Variable Initializations  */
	N		= ((REAL8Sequence*)(data))->length;
	fom		= NULL;
	Y0k		= NULL;
	temp		= NULL;
	YkNk		= NULL;
	temp2		= NULL;
	xbar0k		= NULL;
	x2bar0k		= NULL;
	indexVector	= NULL;

	fomWithoutMargins = NULL;

	/*  length must be greater than one  */
	ASSERT ( N > 3, status, LALMOMENTH_ELNTH, LALMOMENTH_MSGELNTH);

	/*  Create a constant vector values from 1 to N  */
	LALDCreateVector( status->statusPtr, &indexVector, N);

	for (iterator = 1; iterator <= N; iterator++)
		indexVector->data[iterator - 1] = ((REAL8)(iterator));

	/* xbar0k = cumsum(data)./(1:N); */

	LALDCumSum( status->statusPtr, &xbar0k, data);

	LALDVectorDotSlashDVector( status->statusPtr, &temp, xbar0k, indexVector);
        CHECKSTATUSPTR( status );

	for ( iterator = 0; iterator < N; iterator++)
	  {
		xbar0k->data[iterator] = temp->data[iterator];
	  }

	LALDDestroyVector( status->statusPtr, &temp);
        CHECKSTATUSPTR( status );
        temp            = NULL;



        /* x2bar0k = cumsum(data.^2)./(1:N); */

	LALDVectorDotPowerI4(status->statusPtr, &temp, data, two);
        CHECKSTATUSPTR( status );

	LALDCumSum( status->statusPtr, &x2bar0k, temp);

        LALDDestroyVector( status->statusPtr, &temp);
        CHECKSTATUSPTR( status );
        temp            = NULL;

        LALDVectorDotSlashDVector( status->statusPtr, &temp, x2bar0k, indexVector);
        CHECKSTATUSPTR( status );

        for ( iterator = 0; iterator < N; iterator++)
                x2bar0k->data[iterator] = temp->data[iterator];

        LALDDestroyVector( status->statusPtr, &temp);
        CHECKSTATUSPTR( status );
        temp            = NULL;



	/* Y0k = x2bar0k - (xbar0k).^2 */

	LALDCreateVector( status->statusPtr, &Y0k, N);
	CHECKSTATUSPTR( status );

	for ( iterator = 0; iterator < N; iterator++)
		Y0k->data[iterator] = x2bar0k->data[iterator] - pow(xbar0k->data[iterator], two);

	/* xbar0N = xbar0k(end) */
	xbar0N	= xbar0k->data[xbar0k->length - 1];

	/* x2bar0N = x2bar0k(end) */
        x2bar0N  = x2bar0k->data[x2bar0k->length - 1];

	/* Y0N = x2bar0N - (xbar0N)^2 */
	Y0N = x2bar0N - pow(xbar0N,two);

	temp		= NULL;
	temp2		= NULL;
	xbarkNk		= NULL;
	x2barkNk	= NULL;




	/* xbarkNk=cumsum(data(end:-1:1))./(1:N) */

	LALDFlipVector( status->statusPtr, &temp, data);
        CHECKSTATUSPTR( status );

	LALDCumSum( status->statusPtr, &temp2, temp);
        CHECKSTATUSPTR( status );

        LALDVectorDotSlashDVector( status->statusPtr, &xbarkNk, temp2, indexVector);
        CHECKSTATUSPTR( status );

        LALDDestroyVector( status->statusPtr, &temp2);
        CHECKSTATUSPTR( status );

        temp2           = NULL;




	/* x2barkNk=cumsum(data(end:-1:1).^2)./(1:N) */

        LALDVectorDotPowerI4( status->statusPtr, &temp2, temp, two);
        CHECKSTATUSPTR( status );

        LALDDestroyVector( status->statusPtr, &temp);
        CHECKSTATUSPTR( status );

        temp            = NULL;

        LALDCumSum( status->statusPtr, &temp, temp2);
        CHECKSTATUSPTR( status );

	LALDVectorDotSlashDVector( status->statusPtr, &x2barkNk, temp, indexVector);
        CHECKSTATUSPTR( status );

        LALDDestroyVector( status->statusPtr, &temp2);
        LALDDestroyVector( status->statusPtr, &temp);
        CHECKSTATUSPTR( status );

        temp            = NULL;
        temp2           = NULL;




	/* YkNk=x2barkNk-(xbarkNk).^2 */

        LALDVectorDotPowerI4(status->statusPtr, &temp, xbarkNk, two);
        CHECKSTATUSPTR( status );

	LALDCreateVector( status->statusPtr, &YkNk, xbarkNk->length);
        CHECKSTATUSPTR( status );

        for ( iterator = 0; iterator < N; iterator++)
                YkNk->data[iterator] = x2barkNk->data[iterator] - temp->data[iterator];

        LALDDestroyVector( status->statusPtr, &temp);
        CHECKSTATUSPTR( status );

        temp            = NULL;



	/* fom = log(N*Y0N)*(N-1)/2 */

	tempReal = log(N*Y0N)*(N-1)/2;

	LALDCreateVector( status->statusPtr, &fom, N - 3);
        CHECKSTATUSPTR( status );

	/*  lnIRatio is computed here as the results do not depend on data other than the length N  */
	lnIRatio( status->statusPtr, &temp, N, 50 );
        CHECKSTATUSPTR( status );

	for (k = 2; k <= N - 2; k++)
	{
		/*  fom = log(N*Y0N)*(N-1)/2  - log(k.*Y0k(k)).*(k-1)/2  */
		fom->data[k - 2] = tempReal - (log(k * Y0k->data[k - 1])) * (k-1)/2;

		/*  fom = fom - log((N-k).*YkNk(N-k+1)).*(N-k-1)/2  */
		fom->data[k - 2] = fom->data[k - 2] - log((N - k) * YkNk->data[N - k - 1]) * ( N - k - 1) / 2;

		/*  fom = fom + lnIratio(N,50) */
		fom->data[k - 2] = fom->data[k - 2] + temp->data[k - 2];

		/*  fom = exp(fom)*sqrt(2*pi)  */
		fom->data[k - 2] = fom->data[k - 2] + log(2 * LAL_PI)/2.0;

		/*  fom = fom.*rootRatio  */
		rootRatio = ((REAL8)(N))/(k * (N - k));
		fom->data[k - 2] = fom->data[k - 2] + log( rootRatio ) / 2.0;
       	}

	LALDCreateVector( status->statusPtr, &fomWithoutMargins, fom->length - (2 * (*marginOfExclusion)));
	CHECKSTATUSPTR( status );

	for( k = (*marginOfExclusion); k < (INT4)fom->length - 2 * (INT4)((*marginOfExclusion)); k++ )
	  {
	    fomWithoutMargins->data[k - (*marginOfExclusion)] = fom->data[k];
	  }

	/*  [rpeak,ndx] = max(fom)  */
	LALDMax( status->statusPtr, rpeak, fomWithoutMargins, myindex);
        CHECKSTATUSPTR( status );

	for( k = 0; k < (INT4)(fomWithoutMargins->length); k++ )
	  {
	    fomWithoutMargins->data[k] = exp(fomWithoutMargins->data[k]);
	  }	

	/*  r = sum(fom)  */
	LALDSum( status->statusPtr, result, fomWithoutMargins);
	CHECKSTATUSPTR( status );

	(*result) = log((*result));

	/*  ndx = ndx + 2  */
	(*myindex) = (*myindex) + 2;

	LALDDestroyVector( status->statusPtr, &indexVector);
	LALDDestroyVector( status->statusPtr, &temp);
	LALDDestroyVector( status->statusPtr, &xbar0k);
	LALDDestroyVector( status->statusPtr, &x2bar0k);
	LALDDestroyVector( status->statusPtr, &YkNk);
	LALDDestroyVector( status->statusPtr, &Y0k);
	LALDDestroyVector( status->statusPtr, &fom);
        LALDDestroyVector( status->statusPtr, &xbarkNk);
        LALDDestroyVector( status->statusPtr, &x2barkNk);
	LALDDestroyVector( status->statusPtr, &fomWithoutMargins);


	DETATCHSTATUSPTR( status );
	RETURN (status);
}
