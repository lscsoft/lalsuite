/*********************************** <lalVerbatim file="ComputeSkyCV">
Author:  Berukoff, S.J., Papa, M.A.
$Id$
************************************ </lalVerbatim> */

/* <lalLaTeX> 
\subsection{Module \texttt{ComputeSky.c}}\label{ss:ComputeSky.c}
Computes the phase model coefficients necessary for a successful demodulation.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{ComputeSkyCP}
\idx{ComputeSky()}

\subsubsection*{Description}
Given an input index which refers to the sky patch under consideration, this
routine returns the phase model coefficients $A_{s\alpha}$ and $B_{s\alpha}$
which are needed to correctly account for the phase variance of a signal over
time.  The \verb@CSParams@ parameter structure contains relevant information
for this routine to properly run.  In particular, it contains an array of
timestamps in \verb@LIGOTimeGPS@ format, which are the GPS times of the first
data from each SFT.  The \verb@input@ is an \verb@INT4@ variable
\verb@iSkyCoh@, which is the index of the sky patch under consideration.  For
each sky patch, this code needs to be run once; the necessary phase model
coefficients are calculated, and can then be applied to the relevant spindown
parameter sets one is using in their search.  

\subsubsection*{Algorithm}
The routine uses a simplistic nested for-loop structure.  The outer loop is
over the number of SFTs in the observation timescale; this accounts for the
temporal variability of the phase model coefficients.  The inner loop is over
the number of spindown parameters in one set.  Inside the inner loop, the
values are calculated using the analytical formulae given in the
\verb@ComputeSky.h@ documentation.

\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{ComputeSkyCV}}

</lalLaTeX> */

#include <math.h>
#include <lal/LALConstants.h>
#include <lal/ComputeSky.h>

NRCSID (COMPUTESKYC, "$Id ComputeSky.c $");


/* <lalVerbatim file="ComputeSkyCP"> */
void ComputeSky	(LALStatus	*status, 
			REAL8 	*skyConst, 
			INT4 		iSkyCoh, 
			CSParams 	*params)
{  /* </lalVerbatim> */

	INT4	m, n;
	REAL8	t;
	REAL8	basedTbary;
	REAL8	dTbary;
	REAL8	tBary;
	REAL8	tDot;
	REAL8	t0;
	REAL8	tB0;
	

INITSTATUS (status, "ComputeSky", COMPUTESKYC);
ATTATCHSTATUSPTR(status);
	
/* Check for non-negativity of sky positions in SkyCoh[] */
ASSERT(iSkyCoh>=0, status, COMPUTESKYH_ENEGA, COMPUTESKYH_MSGENEGA);

/* Check to make sure sky positions are loaded */
ASSERT(params->skyPos!=NULL, status, COMPUTESKYH_ENULL, COMPUTESKYH_MSGENULL);
ASSERT(params->skyPos!=NULL, status, COMPUTESKYH_ENULL, COMPUTESKYH_MSGENULL);

/* Check to make sure parameters are loaded and reasonable */
ASSERT(params->spinDwnOrder>=0, status, COMPUTESKYH_ENEGA, COMPUTESKYH_MSGENEGA);
ASSERT(params->mObsSFT>=0, status, COMPUTESKYH_ENEGA, COMPUTESKYH_MSGENEGA);
ASSERT(params->tSFT>=0, status, COMPUTESKYH_ENEGA, COMPUTESKYH_MSGENEGA);

for(n=0;n<params->mObsSFT;n++)
{
	ASSERT(params->tGPS[n].gpsSeconds>=0, status, COMPUTESKYH_ENEGA, 	COMPUTESKYH_MSGENEGA);
}

/* Check to make sure pointer to output is not NULL */
ASSERT(skyConst!=NULL, status, COMPUTESKYH_ENNUL, COMPUTESKYH_MSGENNUL);

/* Check to make sure a tdb-like function is being used */
ASSERT((*(params->funcName))!=NULL, status, COMPUTESKYH_ENNUL, COMPUTESKYH_MSGENNUL);

/* calculation of Barycenter.c parameters */
		t0=(REAL8)params->tGPS[0].gpsSeconds*1.0+(REAL8)params->tGPS[0].gpsNanoSeconds*1.0E-9;
	
	(*(params->funcName))(params->skyPos[iSkyCoh], params->skyPos[iSkyCoh+1], t0, &tB0, &tDot);
	
	for (n=0; n<params->mObsSFT; n++) 
	{
		t=(REAL8)(params->tGPS[n].gpsSeconds)+(REAL8)(params->tGPS[n].gpsNanoSeconds)*1.0E-9+0.5*params->tSFT; 
	
		(*(params->funcName))(params->skyPos[iSkyCoh], params->skyPos[iSkyCoh+1], t, &tBary, &tDot);
		
		dTbary = tBary-tB0;		

		for (m=0; m<params->spinDwnOrder+1; m++) 
		{
			basedTbary = pow(dTbary, (REAL8)m);
			skyConst[2*n*(params->spinDwnOrder+1)+2*(INT4)m]=1.0/((REAL8)m+1.0)*basedTbary*dTbary-0.5*params->tSFT*tDot*basedTbary;
			skyConst[2*n*(params->spinDwnOrder+1)+2*(INT4)m+1]= params->tSFT*tDot*basedTbary;
			
 		}
	} 
	  /* writes the length of the vector of sky constants */
	/* Normal Exit */
	DETATCHSTATUSPTR(status);
	RETURN(status);
}

