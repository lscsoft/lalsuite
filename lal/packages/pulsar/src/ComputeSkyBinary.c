/*********************************** <lalVerbatim file="ComputeSkyBinaryCV">
Author:  Messenger, C.J., Berukoff, S.J., Papa, M.A.
$Id$
************************************ </lalVerbatim> */

/* <lalLaTeX> 
\subsection{Module \texttt{ComputeSkyBinary.c}}\label{ss:ComputeSkyBinary.c}
Computes the phase model coefficients necessary for a successful demodulation
for the case of a continuous wave source in a binary system.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{ComputeSkyBinaryCP}
\index{\texttt{ComputeSkyBinary()}}

\subsubsection*{Description}
Given a set of input parameters defining a source location in the sky and 
the binary system in which the source resides, this routine returns the phase 
model coefficients $A_{s\alpha}$ and $B_{s\alpha}$ which are needed to 
correctly account for the phase variance of a signal over time. The 
\verb@CSBParams@ parameter structure contains relevant information
for this routine to properly run.  In particular, it contains an array of
timestamps in \verb@LIGOTimeGPS@ format, which are the GPS times of the first
data from each SFT.  The \verb@input@ is an \verb@INT4@ variable
\verb@iSkyCoh@, which is the index of the sky location under consideration.  For
each sky location and set of orbital parameters, this code needs to be run 
once; the necessary phase model
coefficients are calculated, and can then be applied to the relevant spindown
parameter sets one is using in their search.  

\subsubsection*{Algorithm}
The routine uses a simplistic nested for-loop structure.  The outer loop is
over the number of SFTs in the observation timescale; this accounts for the
temporal variability of the phase model coefficients.  The inner loop is over
the number of spindown parameters in one set.  Inside the inner loop, the
values are calculated using the analytical formulae given in the
\verb@ComputeSkyBinary.h@ documentation.

\subsubsection*{Uses}
\begin{verbatim}
LALBarycenter()
LALBarycenterEarth()
LALDBracketRoot()
LALDBisectionFindRoot()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{ComputeSkyBinaryCV}}

</lalLaTeX> */

#include <math.h>
#include <lal/LALConstants.h>
#include "ComputeSkyBinary.h"
#include <lal/FindRoot.h>
#include "LALBarycenter.h"

NRCSID (COMPUTESKYBINARYC, "$Id$");

static void TimeToFloat(REAL8 *f, LIGOTimeGPS *tgps);
static void FloatToTime(LIGOTimeGPS *tgps, REAL8 *f);
static void Ft(LALStatus *status, REAL8 *tr, REAL8 t, void *tr0);
REAL8 a;      /* semi major axis */
REAL8 Period;    /* Period */
REAL8 ecc;    /* eccentricity */
REAL8 parg;   /* argument of periapse */
LIGOTimeGPS Tperi;  /* Time of periapse passage as measured in SSB frame */
REAL8 p,q;    /* coeficients of phase model */
REAL8 E;      /* eccentric anomoly */

/* <lalVerbatim file="ComputeSkyBinaryCP"> */
void ComputeSkyBinary	(LALStatus	*status, 
		 REAL8 	*skyConst, 
		 INT8 		iSkyCoh, 
		 CSBParams 	*params)
{  /* </lalVerbatim> */
  

	INT4	m, n, nP;
	REAL8	t;
	REAL8   dTsource;
	REAL8	dTbary;
	REAL8   dTbarySP;
	REAL8	tBary;
	REAL8	tB0;
	REAL8   Tperifloat;
	REAL8   Tperi0new;
	REAL8   Tperi0;
	REAL8   dTperi;
	REAL8   dTbin;
	REAL8   dTcoord;
	REAL8   Tdotbin;
	REAL8   basedTperi;
	BarycenterInput baryinput;
	EmissionTime  emit;  
	LIGOTimeGPS tGPS;
	EarthState earth;
	EphemerisData *edat = NULL;
	DFindRootIn input;
	REAL8 tr0, t0;
	REAL8 T0Source;

 INITSTATUS (status, "ComputeSkyBinary", COMPUTESKYBINARYC);
 ATTATCHSTATUSPTR(status);
 
/* Check for non-negativity of sky positions in SkyCoh[] */
 ASSERT(iSkyCoh>=0, status, COMPUTESKYBINARYH_ENEGA, COMPUTESKYBINARYH_MSGENEGA);
 
/* Check to make sure sky positions are loaded */
 ASSERT(params->skyPos!=NULL, status, COMPUTESKYBINARYH_ENULL, COMPUTESKYBINARYH_MSGENULL);
 ASSERT(params->skyPos!=NULL, status, COMPUTESKYBINARYH_ENULL, COMPUTESKYBINARYH_MSGENULL);
 
/* Check to make sure parameters are loaded and reasonable */
 ASSERT(params->spinDwnOrder>=0, status, COMPUTESKYBINARYH_ENEGA, COMPUTESKYBINARYH_MSGENEGA);
 ASSERT(params->mObsSFT>=0, status, COMPUTESKYBINARYH_ENEGA, COMPUTESKYBINARYH_MSGENEGA);
 ASSERT(params->tSFT>=0, status, COMPUTESKYBINARYH_ENEGA, COMPUTESKYBINARYH_MSGENEGA);

/* Check to make sure orbital parameters are loaded and reasonable */
 ASSERT(params->SemiMajorAxis>=0, status, COMPUTESKYBINARYH_ENEGA, COMPUTESKYBINARYH_MSGENEGA);
 ASSERT(params->OrbitalPeriod>0, status, COMPUTESKYBINARYH_ENEGA, COMPUTESKYBINARYH_MSGENEGA);
 ASSERT(params->OrbitalEccentricity>=0, status, COMPUTESKYBINARYH_ENEGA, COMPUTESKYBINARYH_MSGENEGA);
 ASSERT((params->ArgPeriapse>=0)&&(params->ArgPeriapse<=LAL_TWOPI), status, COMPUTESKYBINARYH_ENEGA, COMPUTESKYBINARYH_MSGENEGA);
 /* ASSERT(params->TperiapseSSB!=NULL, status, COMPUTESKYBINARYH_ENEGA, COMPUTESKYBINARYH_MSGENEGA); */

 /* set orbital variables */
 a=params->SemiMajorAxis;
 Period=params->OrbitalPeriod;
 ecc=params->OrbitalEccentricity;
 parg=params->ArgPeriapse;
 Tperi.gpsSeconds=params->TperiapseSSB.gpsSeconds;
 Tperi.gpsNanoSeconds=params->TperiapseSSB.gpsNanoSeconds;
 
 for(n=0;n<params->mObsSFT;n++)
   {
     ASSERT(params->tGPS[n].gpsSeconds>=0, status, COMPUTESKYBINARYH_ENEGA, 	COMPUTESKYBINARYH_MSGENEGA);
   }
 
 /* Check to make sure pointer to output is not NULL */
 ASSERT(skyConst!=NULL, status, COMPUTESKYBINARYH_ENNUL, COMPUTESKYBINARYH_MSGENNUL);
 
 /* prepare params input time structure with first timestamp */
 params->baryinput->tgps.gpsSeconds=params->tGPS[0].gpsSeconds;
 params->baryinput->tgps.gpsNanoSeconds=params->tGPS[0].gpsNanoSeconds;
 
 /* prepare params input sky position structure */
 params->baryinput->alpha=params->skyPos[iSkyCoh];
 params->baryinput->delta=params->skyPos[iSkyCoh+1];
 
 /* Setup barycentric time at start of observation (first timestamp) */
 LALBarycenterEarth(status->statusPtr, params->earth, &(params->baryinput->tgps),params->edat);   
 LALBarycenter(status->statusPtr, params->emit, params->baryinput, params->earth);
 TimeToFloat(&tB0, &(params->emit->te));
 /* printf("first time stamp time in SSB frame tB0 = %lf\n",tB0); */ 
 
 /* calculate phase model coefficients explained in documentation at some point */
 p=((LAL_TWOPI/(Period*1.0))*a*sqrt(1-(ecc*ecc))*cos(parg))-ecc;
 q=(LAL_TWOPI/(Period*1.0))*a*sin(parg);
 /* printf("p = %le, q = %le\n",p,q); */

 /* convert input periapse passage time to a float */
 TimeToFloat(&Tperi0, &(Tperi));
 /* printf("input time of periapse passage in SSB frame Tperi0 = %lf\n",Tperi0); */ 

 /* begin loop over SFT's */
 for (n=0; n<params->mObsSFT; n++) 
   {
     /* Calculate the detector time at the mid point of current SFT */
     t=(REAL8)(params->tGPS[n].gpsSeconds)+(REAL8)(params->tGPS[n].gpsNanoSeconds)*1.0E-9+0.5*params->tSFT;    
     FloatToTime(&(params->baryinput->tgps), &t);

     /* Convert this mid point detector time into barycentric time */
     LALBarycenterEarth(status->statusPtr, params->earth, &(params->baryinput->tgps), params->edat);    
     LALBarycenter(status->statusPtr, params->emit, params->baryinput, params->earth);       

     /* Calculate the time difference since the orbit epoch (observed periapse passage) in barycentric time */    
     TimeToFloat(&tBary, &(params->emit->te));	
     dTbary=tBary-Tperi0;  
     /* printf("SSB time of timestamp tbary = %lf\n",tBary); */
     /* printf("SSB time since periapse passage dTbary = %lf\n",dTbary); */ 

     /* make sure that it is less than a single period (SP) for use in the root finding procedure */
     dTbarySP=Period*((dTbary/(1.0*Period))-floor(dTbary/(1.0*Period)));
     
     /* Calculate number of full orbits completed since orbit epoch */
     nP=(INT4)floor(dTbary/(1.0*Period));
     /* printf("New SSB time since periapse passage dTbary = %lf\n",dTbary); */

     /* begin root finding procedure */
     tr0 = dTbarySP;
     input.function = Ft;
     input.xmin = 0.0;
     input.xmax = LAL_TWOPI;
     input.xacc = 1E-12;

     /* expand domain until a root is bracketed */
     LALDBracketRoot(status->statusPtr,&input,&tr0); 

     /* bisect domain to find eccentric anomoly E corresponding to the current timestamp */
     LALDBisectionFindRoot(status->statusPtr,&E,&input,&tr0); 
    
     /* retarded time interval since input orbit epoch */ 
     dTperi=(Period/LAL_TWOPI)*(E-(ecc*sin(E)))+((REAL8)nP*Period); 
     /* printf("time since orbit epoch in retarded frame dTperi %lf\n",dTperi); */

     /* calculate dtbin/dtssb */
     dTcoord=(1.0-(ecc*cos(E)))/(1.0+(p*cos(E))-(q*sin(E)));  

     /* calculate dtbin/dtdet */
     Tdotbin = params->emit->tDot*dTcoord;   
     /* printf("dTcoord = %le, Tdotbin = %le\n",dTcoord,Tdotbin); */

     /* Loop over all spin down orders plus 0th order (f0) */
     for (m=0; m<params->spinDwnOrder+1; m++)     
       {
	 basedTperi = pow(dTperi, (REAL8)m);
	 /* Calculate A coefficients */
	 skyConst[2*n*(params->spinDwnOrder+1)+2*(INT4)m]=1.0/((REAL8)m+1.0)*basedTperi*dTperi-0.5*params->tSFT*basedTperi*Tdotbin;
	 /* Calculate B coefficients */
	 skyConst[2*n*(params->spinDwnOrder+1)+2*(INT4)m+1]= params->tSFT*basedTperi*Tdotbin;
	 /* printf("A = %le, B = %le\n",skyConst[2*n*(params->spinDwnOrder+1)+2*(INT4)m],skyConst[2*n*(params->spinDwnOrder+1)+2*(INT4)m+1]); */
       }    
   } 

 /* Normal Exit */
 DETATCHSTATUSPTR(status);
 RETURN(status);
}


/* Internal routines */
static void TimeToFloat(REAL8 *f, LIGOTimeGPS *tgps)
{
  INT4 x, y;

  x=tgps->gpsSeconds;
  y=tgps->gpsNanoSeconds;
  *f=(REAL8)x+(REAL8)y*1.e-9;
}

static void FloatToTime(LIGOTimeGPS *tgps, REAL8 *f)
{
  REAL8 temp0, temp2, temp3;
  REAL8 temp1, temp4;
  
  temp0 = floor(*f);     /* this is tgps.S */
  temp1 = (*f) * 1.e10;
  temp2 = fmod(temp1, 1.e10);
  temp3 = fmod(temp1, 1.e2); 
  temp4 = (temp2-temp3) * 0.1;

  tgps->gpsSeconds = (INT4)temp0;
  tgps->gpsNanoSeconds = (INT4)temp4;
}

static void Ft(LALStatus *status, REAL8 *tr, REAL8 E, void *tr0)
{
  INITSTATUS(status, "Ft", "Function Ft()");
  ASSERT(tr0,status, 1, "Null pointer");
  /* this is the function relating the observed time since periapse in the SSB to the eccentric anomoly E */
  *tr = *(REAL8 *)tr0*(-1.0) + (Period/LAL_TWOPI)*(E+(p*sin(E))+q*(cos(E)-1.0));
  RETURN(status);
}
                       
