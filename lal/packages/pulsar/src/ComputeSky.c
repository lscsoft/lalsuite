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
\index{\verb&ComputeSky()&}

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

static void tdb(REAL8 alpha, REAL8 delta, REAL8 t_GPS, REAL8 *T, REAL8 *Tdot, const CHAR *sw);


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


/* calculation of Barycenter.c parameters */
		t0=(REAL8)params->tGPS[0].gpsSeconds*1.0+(REAL8)params->tGPS[0].gpsNanoSeconds*1.0E-9;
	
	tdb(params->skyPos[iSkyCoh], params->skyPos[iSkyCoh+1], t0, &tB0, &tDot, params->sw);
	
	for (n=0; n<params->mObsSFT; n++) 
	{
		t=(REAL8)(params->tGPS[n].gpsSeconds)+(REAL8)(params->tGPS[n].gpsNanoSeconds)*1.0E-9+0.5*params->tSFT; 
	
		tdb(params->skyPos[iSkyCoh], params->skyPos[iSkyCoh+1], t, &tBary, &tDot,params->sw);
		
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


/* internal routine to convert detector time to barycentric time */
/* Note: this routine will be replaced by Barycenter(), written  */
/* by C. Cutler, Albert-Einstein-Institut.  			 */

static void tdb(REAL8 alpha, REAL8 delta, REAL8 t_GPS, REAL8 *T, REAL8 *Tdot, const CHAR *sw) 
{

	/*RA alpha and DEC delta have their usual meanings, but units 
  are radians, not hours or degrees */

	REAL8 tdiff;
	REAL8 tdiff_day;
	REAL8 tdiff_yr; 
	REAL8 beta;
	REAL8 x;
	REAL8 y; 
	REAL8 lam;
	REAL8 phi0;
	REAL8 psi0; 
	REAL8 phi = 0;
	REAL8 psi = 0;
	REAL8 cose;
	REAL8 sine;
	REAL8 sind; 
	REAL8 cosd;
	REAL8 sina;
	REAL8 cosa;
	REAL8 cosb;

	/*phi0 and psi0 given location and orientation, resp., of Earth at t_GPS =0 */
	REAL8 T_bary, Tdot_bary = 0;
	/* T_bary (what we want to know) is pulse emission time in TDB coords.*/
	/* The lab measures GPS time, but we set t_GPS =0 at first instant of 1998*/
	REAL8 AU = 1.49597870e13/2.99792458e10;
	REAL8 RE = 6.37814e8/2.99792458e10;
	REAL8 eps = (23.0 + 26.0/60 + 21.488/3600)*LAL_PI/180;
	/*    printf("eps is %lf \n",eps); */
	phi0=2.86e0;
	psi0=4.12e0;
	/* Those values of phi0 and psi0 are completely made up */

	/* Now converting from equatorial to ecliptic coords. */
	sine=sin(eps);
	cose=cos(eps);
	sind=sin(delta);
	cosd=cos(delta);
	sina=sin(alpha);
	cosa=cos(alpha);
	beta=cose*sind-sine*cosd*sina;
	cosb=sqrt(1.0-beta*beta);
	x=cosd*cosa/cosb;
	y=(sine*sind+cose*cosd*sina)/cosb;
	lam=atan2(y,x);
	/*  printf("beta, lam are (rad) %lf  %lf \n", beta, lam);  */
	phi=phi0+2*LAL_PI*t_GPS/3.15576e7;
	psi=psi0+2*LAL_PI*t_GPS/8.64e4;
	tdiff_day=RE*cosd*cos(psi-alpha);
	tdiff_yr=AU*cosb*cos(phi-lam);
	tdiff=tdiff_day+tdiff_yr;
	T_bary=t_GPS+tdiff;
	Tdot_bary=1.e0-RE*cosd*sin(psi-alpha)*2*LAL_PI/8.64e4 		
		-AU*cosb*sin(phi-lam)*2*LAL_PI/3.15576e7;
	
	/* switch to turn mod off/on */
	/* note that 'off' negates the function of the code! */
	if(sw!=NULL){
	*T=T_bary; 
 	*Tdot=Tdot_bary;}
	
	else {
	*T=t_GPS;
	*Tdot=1.0;}
	
}
