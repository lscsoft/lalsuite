/*  <lalVerbatim file="LALEOBWaveformCV">
Author: Sathyaprakash, B. S., Cokelaer T.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALEOBWaveform.c} and 
\texttt{LALEOBWaveformTemplates.c}}

Module to generate effective-one-body waveforms.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALEOBWaveformCP}
\index{\verb&LALEOBWaveform()&}
\begin{itemize}
\item {\tt signal:} Output containing the inspiral waveform.
\item {\tt params:} Input containing binary chirp parameters.
\end{itemize}

\input{LALEOBWaveformTemplatesCP}
\index{\verb&LALEOBWaveformTemplates()&}
\begin{itemize}
\item {\tt signal1:} Output containing the 0-phase inspiral waveform.
\item {\tt signal2:} Output containing the $\pi/2$-phase inspiral waveform.
\item {\tt params:} Input containing binary chirp parameters.
\end{itemize}

\input{LALEOBWaveformForInjectionCP}
\index{\verb&LALEOBWaveformForInjection()&}
\begin{itemize}
\item {\tt inject\_hc:} Output containing the 0-phase inspiral waveform.
\item {\tt inject\_hp:} Output containing the $\pi/2$-phase inspiral waveform.
\item {\tt inject\_phase:} Output containing the phase of inspiral waveform.
\item {\tt inject\_freq:} Output containing the frequency of inspiral waveform.
\item {\tt params:} Input containing binary chirp parameters.
\end{itemize}


\subsubsection*{Description}
By solving four coupled ordinary differential equations in
Eq.~(\ref{eq:3.28})-(\ref{3.31}) this module computes the
waveform in Eq.~(\ref{4.1}) (see discussion in Sec.~\ref{sec:EOB}
for details on how the initial conditions are chosen, when the
waveform is terminated and so on). 
No quasi-normal mode oscillations are added to the plunge signal
so the waveform is terminated around $2.8\,M$.
\subsection*{3PN vs 2PN}
At 3PN, two additional parameters exist namely OmegaS and Zeta2. 
The first parameters should be set to zero. If the  second parameter
is also set to zero then the waveform correponds to the standard 
waveforms. 
\subsubsection*{Algorithm}
A fourth order Runge-Kutta is used to solve the differential equations.

\subsubsection*{Uses}
\begin{verbatim}
   LALInspiralSetup 
   LALInspiralChooseModel
   LALInspiralVelocity
   LALInspiralPhasing1
   LALDBisectionFindRoot
   LALRungeKutta4
   LALHCapDerivatives
   LALHCapDerivatives3PN
   LALlightRingRadius
   LALlightRingRadius3PN
   LALpphiInit
   LALpphiInit3PN
   LALprInit
   LALprInit3PN
   LALrOfOmega
   LALrOfOmega3PN
\end{verbatim}

\subsubsection*{Notes}
The length of the waveform returned by {\tt LALInspiralWaveLength} is 
occassionally smaller than what is required to hold an EOB waveform.
This is because EOB goes beyond the last stable orbit up to the light
ring while {\tt LALInspiralWaveLength} assumes that the waveform terminates
at the last stable orbit. It is recommended that a rather generous 
{\tt params->nEndPad} be used to prevent the code from crashing.

\vfill{\footnotesize\input{LALEOBWaveformCV}}

</lalLaTeX>  */
#include <lal/Units.h>
#include <lal/LALInspiral.h>
#include <lal/FindRoot.h>
#include <lal/SeqFactories.h>

typedef struct tagrOfOmegaIn {
   REAL8 eta, omega;
} rOfOmegaIn;

typedef struct tagPr3In {
  REAL8 eta, zeta2, omegaS, omega, vr,r,q;
  InspiralDerivativesIn in3copy; 
} pr3In;

void 
omegaofr3PN (
	     REAL8 *x,
	     REAL8 r, 
	     void *params) ;

static 
void LALHCapDerivatives(	REAL8Vector *values, 
				REAL8Vector *dvalues, 
				void 		*funcParams);

static 
void LALprInit(	REAL8 		*pr, 
				REAL8 		r, 
				InspiralDerivativesIn 	*ak);

static 
void LALpphiInit(	REAL8 *phase, 
			REAL8 r, 
			REAL8 eta);

static 
void LALlightRingRadius(	LALStatus 	*status,
				REAL8 		*x, 
				REAL8 		r, 
				void 	*params);
							
static void LALrOfOmega (	LALStatus 	*status, 
				REAL8 		*x,
				REAL8 		r, 
				void 		*params);

static
void LALHCapDerivatives3PN(	REAL8Vector 	*values,
				REAL8Vector 	*dvalues, 
				void 			*funcParams);
							
static
void LALprInit3PN(LALStatus *status, REAL8 *pr , REAL8 , void  *params);

static
void LALpphiInit3PN(REAL8 *phase, REAL8 r, REAL8 eta, REAL8 omegaS);

static
void LALlightRingRadius3PN(LALStatus *status, REAL8 *x, REAL8 r, void *params);

static
void LALrOfOmega3PN (LALStatus *status, REAL8 *x, REAL8 r, void *params);

static
void LALvr3PN(REAL8 *vr, void *params);


NRCSID (LALEOBWAVEFORMC, 
"$Id$");

/*--------------------------------------------------------------------*/

static void 
LALprInit(
   REAL8 *pr, 
   REAL8 r, 
   InspiralDerivativesIn *ak
   ) 
{ 

/* 
	This combines Eq. (4.13), (4.14) (4.16) of BD2 to get 
	the initial value of pr 
*/

   REAL8 	eta,
			z,
			omega,
			jadiab,
			djadiabdr, 
			H0cap, 
			v,
			FDIS,
			A,
			r2,  /* temp variable */
			r3,  /* temp variable */
			cr;

   eta = ak->coeffs->eta;
   r2 = r*r;
   r3 = r2*r;

   A = 1. - 2./r + 2.*eta/r3;
   z = r3 * A * A/(r3 - 3.*r2 + 5.*eta);
   omega = sqrt((1.-3*eta/r2)/(r3 * (1. + 2*eta*(sqrt(z)-1.))));
   jadiab = sqrt (r2 * (r2 - 3.*eta)/(r3 - 3.*r2 + 5.*eta));
   djadiabdr = (r3*r3 - 6.*r3*r2 + 3.*eta*r2*r2 +20.*eta*r3-30.*eta*eta*r)/
               (pow(r3 - 3.*r2 + 5.*eta, 2.)*2.*jadiab);
   H0cap = sqrt(1. + 2.*eta*(-1. + sqrt(z)))/eta;
   cr = A*A/((1.-6.*eta/r2) * eta * H0cap * sqrt(z));
   v = pow(omega,oneby3);
   FDIS = -ak->flux(v, ak->coeffs)/(eta * v*v*v); 
   *pr = FDIS/(djadiabdr*cr);
}

/*--------------------------------------------------------------------*/
static void 
LALpphiInit(
   REAL8 *phase, 
   REAL8 r, 
   REAL8 eta
   )
{ 
   REAL8 r2, r3;
   r2 = r*r;
   r3 = r2*r;
   *phase = pow(r2 * (r2 - 3.*eta) / (r3 - 3.* r2 + 5.*eta), 0.5);
}

/*--------------------------------------------------------------------*/
NRCSID (LALROFOMEGAC,
"$Id$");

static void 
LALrOfOmega (
   LALStatus *status, 
   REAL8 *x, 
   REAL8 r, 
   void *params
   ) 
{ 

   REAL8 r2, r3, A, z, eta, omega;
   rOfOmegaIn *rofomegain;

   INITSTATUS(status, "LALrOfOmega", LALROFOMEGAC);
   ATTATCHSTATUSPTR(status);
   rofomegain = (rOfOmegaIn *) params;
   eta = rofomegain->eta;
   omega = rofomegain->omega;
   r2 = r*r;
   r3 = r2*r;

   A = 1. - 2./r + 2.*eta/r3;
   z = r3 * A * A/(r3 - 3.*r2 + 5.*eta);
   *x = omega - sqrt((1.-3*eta/r2)/(r3 * (1. + 2*eta*(sqrt(z)-1.))));
   DETATCHSTATUSPTR(status);
   RETURN(status);
}

/*--------------------------------------------------------------------*/
NRCSID (LALLIGHTRINGRADIUSC, 
"$Id$");
static void 
LALlightRingRadius(
   LALStatus *status, 
   REAL8 *x, 
   REAL8 r, 
   void *params
   ) 
{ 
   REAL8 eta;
   rOfOmegaIn *rofomegain;

   INITSTATUS(status, "LALlightRingRadius", LALLIGHTRINGRADIUSC);
   ATTATCHSTATUSPTR(status);
   rofomegain = (rOfOmegaIn *) params;
   eta = rofomegain->eta;
   *x = pow(r,3.) - 3.*r*r + 5.* eta;
   DETATCHSTATUSPTR(status);
   RETURN(status);
}


static void 
LALHCapDerivatives(
   REAL8Vector *values, 
   REAL8Vector *dvalues, 
   void *funcParams
   ) 
{ 
   REAL8 r, s, p, q, r2, r3, p2, q2, A, B, dA, dB, hcap, Hcap, etahH;
   REAL8 omega, v, eta;
   InspiralDerivativesIn *ak;

   ak = (InspiralDerivativesIn *) funcParams;

   eta = ak->coeffs->eta;

   r = values->data[0];
   s = values->data[1];
   p = values->data[2];
   q = values->data[3];

   r2 = r*r;
   r3 = r2*r;
   p2 = p*p;
   q2 = q*q;
   A = 1. - 2./r + 2.*eta/r3;
   B = (1. - 6.*eta/r2)/A;
   dA = 2./r2 - 6.*eta/pow(r,4.);
   dB = (-dA * B + 12.*eta/r3)/A;
   hcap = pow (A*(1. + p2/B + q2/r2), 0.5);
   Hcap = pow (1. + 2.*eta*(hcap - 1.), 0.5) / eta;
   etahH = eta*hcap*Hcap;

   dvalues->data[0] = A*A*p/((1. - 6.*eta/r2) * etahH);
   dvalues->data[1] = omega = A * q / (r2 * etahH);

   v = pow(omega,oneby3);

   dvalues->data[2] = -0.5 * (dA * hcap * hcap/A - p2 * A * dB/(B*B) 
   									- 2. * A * q2/r3) / etahH;
   dvalues->data[3] = -ak->flux(v, ak->coeffs)/(eta * v*v*v); 
   /*
   printf("%e %e %e %e %e %e %e %e\n", r, s, p, q, A, B, hcap, Hcap);
   */
}



/*--------------------------------------------------------------------*/

  void
LALpphiInit3PN(
	    REAL8 *phase,
	    REAL8 r,
	    REAL8 eta,
	    REAL8 omegaS
	    )
{
  REAL8 u, u2, u3,  a4, a4p4eta, a4peta2, NA, DA, A, dA;
  

  u = 1./r;
  u2 = u*u;
  u3 = u2*u;
  a4 = (ninty4by3etc - 2. * omegaS) * eta;
  a4p4eta = a4 + 4. * eta;
  a4peta2 = a4 + eta * eta;
  NA = 2.*(4.-eta) + (a4 - 16. + 8. * eta) * u;
  DA = 2.*(4.-eta) + a4p4eta * u + 2. * a4p4eta * u2 + 4.*a4peta2 * u3;
  A = NA/DA;
  dA = ( (a4 - 16. + 8. * eta) * DA - NA * (a4p4eta + 4. * a4p4eta * u 
			+ 12. * a4peta2  * u2))/(DA*DA);

  *phase = sqrt(-dA/(2.*u*A + u2 * dA));

}
/*---------------------------------------------------------------*/

NRCSID (LALPRINIT3PN, 
"$Id$");
 void 
LALprInit3PN(
	     LALStatus *status, 
	     REAL8 *pr,
	     REAL8 p,
	     void *params
	     ) 
{
  REAL8   u, u2, u3, u4, p2, p3, p4, q2, A, DA, NA;
  REAL8  onebyD, AbyD, Heff, HReal, etahH;
  REAL8 omegaS, eta, a4, a4p4eta, a4peta2, z3, r, vr, q;
  pr3In *ak;
  
  INITSTATUS(status, "LALprInit3PN", LALPRINIT3PN);
  ATTATCHSTATUSPTR(status);
  ak = (pr3In *) params;
  
  eta = ak->eta;
  vr = ak->vr;
  r = ak->r;
  q = ak->q;
  omegaS = ak->omegaS;
  
     
   p2 = p*p;
   p3 = p2*p;
   p4 = p2*p2;
   q2 = q*q;
   u = 1./ r;
   u2 = u*u;
   u3 = u2 * u;
   u4 = u2 * u2;
   z3 = 2. * (4. - 3. * eta) * eta;
   a4 = (ninty4by3etc - 2. * omegaS) * eta;
   a4p4eta = a4 + 4. * eta;
   a4peta2 = a4 + eta * eta;  

/* From DJS 14, 15 */
   NA = 2.*(4.-eta) + (a4 - 16. + 8. * eta) * u;
   DA = 2.*(4.-eta) + a4p4eta * u + 2. * a4p4eta * u2 + 4.*a4peta2 * u3;
   A = NA/DA;
   onebyD = 1. + 6.*eta*u2 + 2. * ( 26. - 3. * eta) * eta * u3;
   AbyD = A * onebyD;   

   Heff = pow (A*(1. + AbyD * p2 + q*q * u2 + z3 * p4 * u2), 0.5);
   HReal = pow (1. + 2.*eta*(Heff - 1.), 0.5) / eta;
   etahH = eta*Heff*HReal;
   
   *pr = -vr +  A*(AbyD*p + 2. * z3 * u2 * p3)/etahH; 
   
   DETATCHSTATUSPTR(status);
   RETURN(status);
}

void 
omegaofr3PN (
	     REAL8 *x,
	     REAL8 r, 
	     void *params) 
{
   REAL8 u, u2, u3, a4, a4p4eta, a4peta2, eta, NA, DA, A, dA;
   REAL8   omegaS;

   /*include a status here ?*/
   pr3In *ak;
   ak = (pr3In *) params;
   omegaS = ak->omegaS;
   eta = ak->eta;

   u = 1./r;
   u2 = u*u;
   u3 = u2*u;
   a4 = (ninty4by3etc - 2. * omegaS) * eta;
   
   a4p4eta = a4 + 4. * eta;
   a4peta2 = a4 + eta * eta;
   NA = 2.*(4.-eta) + (a4 - 16. + 8. * eta) * u;
   DA = 2.*(4.-eta) + a4p4eta * u + 2. * a4p4eta * u2 + 4.*a4peta2 * u3;
   A = NA/DA;
   dA = ( (a4 - 16. + 8. * eta) * DA - NA 
			* (a4p4eta + 4. * a4p4eta * u + 12. * a4peta2  * u2))/(DA*DA);
   *x = pow(u,1.5) 
   			* sqrt ( -0.5 * dA /(1. + 2.*eta * (A/sqrt(A+0.5 * u*dA)-1.)));

}

void 
LALrOfOmega3PN(
	    LALStatus *status, 
	    REAL8 *x, 
	    REAL8 r, 
	    void *params)
{
  REAL8  omega1,omega2,eta ;
  pr3In *pr3in;
       
  status = NULL;
  pr3in = (pr3In *) params;
  eta = pr3in->eta;

  omega1 = pr3in->omega;
  omegaofr3PN(&omega2,r, params);
  *x = -omega1 + omega2;

}
/*--------------------------------------------------------------------*/
NRCSID (LALLIGHTRINGRADIUS3PNC,
"$Id$");
 void 
LALlightRingRadius3PN(
		      LALStatus *status, 
		      REAL8 *x, 
		      REAL8 r, 
		      void *params
		      ) 
{ 
  REAL8 eta, u, u2, u3, a4, a4p4eta,a4peta2, NA, DA, A, dA;
  rOfOmegaIn *rofomegain;
  REAL8 omegaS=0;
  status = NULL;
  rofomegain = (rOfOmegaIn *) params;
  eta = rofomegain->eta;


  u = 1./r;
  u2 = u*u;
  u3 = u2*u;
  a4 = (ninty4by3etc - 2. * omegaS) * eta;
  a4p4eta = a4 + 4. * eta;
  a4peta2 = a4 + eta * eta;
  NA = 2.*(4.-eta) + (a4 - 16. + 8. * eta) * u;
  DA = 2.*(4.-eta) + a4p4eta * u + 2. * a4p4eta * u2 + 4.*a4peta2 * u3;
  A = NA/DA;
  dA = ( (a4 - 16. + 8. * eta) * DA - NA * (a4p4eta + 4. 
	* a4p4eta * u + 12. * a4peta2  * u2))/(DA*DA);
  *x = 2 * A + dA * u;
}
/*--------------------------------------------------------------------*/

 void 
LALHCapDerivatives3PN(
		  REAL8Vector *values,
		  REAL8Vector *dvalues,
		  void *funcParams
		  )
{
   REAL8 r, s, p, q, u, u2, u3, u4, p2, p3, p4, q2, Apot, DA, NA;
   REAL8  dA, onebyD, DonebyD, AbyD, Heff, HReal, etahH;
   REAL8 omega, v, eta, a4, a4p4eta, a4peta2, z2, z30, z3, zeta2;
   REAL8 a5, n1, c1, d1, d2, d3, oneby4meta, vu, omegaS;
   REAL8    flexNonAdiab = 0;
   REAL8    flexNonCirc = 0;

   InspiralDerivativesIn *ak;

   ak = (InspiralDerivativesIn *) funcParams;
   eta = ak->coeffs->eta;
   zeta2 = ak->coeffs->zeta2;
   omegaS = ak->coeffs->omegaS;
   



   r = values->data[0];
   s = values->data[1];
   p = values->data[2];
   q = values->data[3];

   p2 = p*p;
   p3 = p2*p;
   p4 = p2*p2;
   q2 = q*q;
   u = 1./r;
   u2 = u*u;
   u3 = u2 * u;
   u4 = u2 * u2;
   z30 = 2.L * (4.L - 3.L * eta) * eta;
   z2 = 0.75L * z30 * zeta2,
   z3 = z30 * (1.L - zeta2);

   a4 = (ninty4by3etc - 2. * omegaS) * eta;
   a5 = 0;
   a4p4eta = a4 + 4. * eta;
   a4peta2 = a4 + eta * eta;
   
   /* From DJS 3PN Hamiltonian */
   oneby4meta = 1./(4.-eta);
   n1 = 0.5 * (a4 - 16. + 8.*eta) * oneby4meta;
   d1 = 0.5 * a4p4eta * oneby4meta;
   d2 = a4p4eta * oneby4meta;
   d3 = 2. * a4peta2 * oneby4meta;
   NA = 1. + n1 * u;
   DA = 1 + d1*u + d2*u2 + d3*u3;
   Apot = NA/DA;
   
   onebyD = 1. + 6.*eta*u2 + (2. * ( 26. - 3. * eta) * eta - z2)* u3;
   AbyD = Apot * onebyD;
   Heff = pow (Apot*(1. + AbyD * p2 + q*q * u2 + z30 * (p4 + zeta2*(-0.25*p4
        + 0.75  * p2 * q2 * u2 )) * u2), 0.5);
   HReal = pow (1. + 2.*eta*(Heff - 1.), 0.5) / eta;

   dA = -u2/(DA*DA) * (n1*DA - NA * (d1 + 2.*d2*u + 3.*d3*u2));

   DonebyD = -12.*eta*u3 - (6.*(26. - 3.*eta)*eta - z2)*u4;
   etahH = eta*Heff*HReal;

   dvalues->data[0] = Apot*(AbyD*p +  z30 * u2 *(2.* p3 
              + zeta2*(-0.5*p3 + 0.75*p*q2*u2)))/etahH;
   dvalues->data[1] = omega = Apot * q * u2 * (1. + 0.75*z30*zeta2*p2*u2)/ etahH;
   v = pow(omega,oneby3);

   dvalues->data[2] = -0.5 * Apot * (dA*Heff*Heff/(Apot*Apot) - 2.*q2*u3 
              + (dA * onebyD + Apot * DonebyD) * p2
              + z30 * u3 *(-2.* p4+zeta2*(0.5*p4 - 3.0*p2*q2*u2))) / etahH;

   c1 = 1.+(u2 - 2.*u3*Apot/dA) * q2;
   dvalues->data[3] = -(1. - flexNonAdiab*c1) * (1. + flexNonCirc*p2/(q2*u2)) 
   					* ak->flux(v,ak->coeffs)/(eta * v*v*v);
   vu = pow(u,0.5);
}



/*----------------------------------------------------------------------*/
 void LALvr3PN(REAL8 *vr, void *params ) 
{
  REAL8 A, dA, d2A, NA, DA, dDA, dNA, d2DA;
  REAL8 u, u2, u3, v, x1; 
  REAL8 eta,a4, a4p4eta, a4peta2, FDIS;
  
  pr3In *pr3in;
  pr3in = (pr3In *)params;

  eta = pr3in->eta;     
  u = 1./ pr3in->r;
  
  u2 = u*u;
  u3 = u2*u;

  
  a4 = (ninty4by3etc - 2. * pr3in->omegaS) * eta;
  a4p4eta = a4 + 4. * eta;
  a4peta2 = a4 + eta * eta;
  NA = 2.*(4.-eta) + (a4 - 16. + 8. * eta) * u;
  DA = 2.*(4.-eta) + a4p4eta * u + 2. * a4p4eta * u2 + 4.*a4peta2 * u3;
  A = NA/DA;
  dNA = (a4 - 16. + 8. * eta);
  dDA = (a4p4eta + 4. * a4p4eta * u + 12. * a4peta2  * u2);
  d2DA = 4. * a4p4eta + 24. * a4peta2 * u;
  
  dA = (dNA * DA - NA * dDA)/ (DA*DA);
  d2A = (-NA * DA * d2DA - 2. * dNA * DA * dDA + 2. * NA * dDA * dDA)/pow(DA,3.);
  v = pow(pr3in->omega,oneby3);
  FDIS = -pr3in->in3copy.flux(v, pr3in->in3copy.coeffs)/(eta* pr3in->omega);
  x1 = -1./u2 * sqrt (-dA * pow(2.* u * A + u2 * dA, 3.) ) 
  		/ (2.* u * dA * dA + A*dA - u * A * d2A);
  *vr = FDIS * x1;
}


/*-------------------------------------------------------------------*/








/*  <lalVerbatim file="LALEOBWaveformCP"> */
void 
LALEOBWaveform (
   LALStatus        *status,
   REAL4Vector      *signal,
   InspiralTemplate *params
   ) 
{ /* </lalVerbatim> */

   INT4 count, nn=4;
   REAL8 eta, m, rn, r, rOld, s, p, q, dt, t, h, v, omega, f, dR;
   REAL8Vector dummy, values, dvalues, newvalues, yt, dym, dyt;
   TofVIn in1;
   InspiralPhaseIn in2;
   InspiralDerivativesIn in3;
   rk4In in4;
   pr3In pr3in; 
   void *funcParams;
   expnCoeffs ak;
   expnFunc func;
   rOfOmegaIn rofomegain;
   DFindRootIn rootIn;
   
   INITSTATUS(status, "LALEOBWaveform", LALEOBWAVEFORMC);
   ATTATCHSTATUSPTR(status);

   ASSERT(signal,  status, 
	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signal->data,  status, 
   	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(params,  status, 
   	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(params->nStartPad >= 0, status, 
   	LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->nEndPad >= 0, status, 
   	LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->fLower > 0, status, 
   	LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->tSampling > 0, status, 
   	LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->totalMass > 0., status, 
   	LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   LALInspiralSetup (status->statusPtr, &ak, params);
   CHECKSTATUSPTR(status);
   LALInspiralChooseModel(status->statusPtr, &func, &ak, params);
   CHECKSTATUSPTR(status);

/*   params->nStartPad=0.0;
   params->nEndPad=0.0;
  */ 
   ASSERT(ak.totalmass/LAL_MTSUN_SI > 0.4, status, 
   	LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
/*   ASSERT(ak.totalmass/LAL_MTSUN_SI < 100, status, 
		LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);*/

/* Allocate all the memory required to dummy and then point the various
   arrays to dummy - this makes it easier to handle memory failures */

   dummy.length = nn * 6;

   values.length = dvalues.length = newvalues.length =
   yt.length = dym.length = dyt.length = nn;

   if (!(dummy.data = (REAL8 * ) LALMalloc(sizeof(REAL8) * nn * 6))) {
      ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
   }

   values.data = &dummy.data[0];
   dvalues.data = &dummy.data[nn];
   newvalues.data = &dummy.data[2*nn];
   yt.data = &dummy.data[3*nn];
   dym.data = &dummy.data[4*nn];
   dyt.data = &dummy.data[5*nn];

   dt = 1./params->tSampling;
   eta = ak.eta;
   m = ak.totalmass;

   t = 0.0;
   in1.t = t;
   in1.t0 = ak.t0;
   in1.v0 = ak.v0;
   in1.vlso = ak.vlso;
   in1.totalmass = ak.totalmass;
   in1.dEnergy = func.dEnergy;
   in1.flux = func.flux;
   in1.coeffs = &ak;

   LALInspiralVelocity(status->statusPtr, &v, &in1);
   CHECKSTATUSPTR(status);

   omega = pow(v,3.);
   f = omega/(LAL_PI*m);

   in2.v0 = ak.v0;
   in2.phi0 = params->startPhase;
   in2.dEnergy = func.dEnergy;
   in2.flux = func.flux;
   in2.coeffs = &ak;
   LALInspiralPhasing1(status->statusPtr, &s, v, &in2);
   CHECKSTATUSPTR(status);
   s = s/2.;

   rofomegain.eta = eta;
   rofomegain.omega = omega;

   rootIn.xacc = 1.0e-16;
   switch (params->order)
     {
     case twoPN:
     case twoPointFivePN:
       rootIn.function = LALlightRingRadius;
       break;
     case threePN:
       rootIn.function = LALlightRingRadius3PN;
       break;
     default:
       fprintf(stderr, 
	 			"There are no EOB waveforms implemented at order %d\n", 
	 			params->order);
       exit(0);
     }
   rootIn.xmax = 2.;
   rootIn.xmin = 4.;
   funcParams = (void *) &rofomegain;
/*-------------------------------------------------------------------
Userful for debugging: Make sure a solution for r exists.
--------------------------------------------------------
   for (r=0.01; r<10.; r+=.01) {
      LALlightRingRadius(status->statusPtr, &x, r, funcParams);
      CHECKSTATUSPTR(status);
      printf("%e %e\n", r, x);
   }
   printf("&\n");
   for (r=1; r<100.; r+=.1) {
      LALrOfOmega(status->statusPtr, &x, r, funcParams);
      CHECKSTATUSPTR(status);
      printf("%e %e\n", r, x);
   }
-------------------------------------------------------------------*/
   LALDBisectionFindRoot(status->statusPtr, &rn, &rootIn, funcParams);
   CHECKSTATUSPTR(status);
	
   
   
    /*rofOmega */
    switch (params->order)
     {
     case twoPN:
     case twoPointFivePN:
       rootIn.function = LALrOfOmega;
       rootIn.xmax = 300.;
       rootIn.xmin = 3.;
       LALDBisectionFindRoot(status->statusPtr, &r, &rootIn, funcParams);
       CHECKSTATUSPTR(status);
       break;
     case threePN:
       rootIn.function = LALrOfOmega3PN;
       pr3in.eta = eta;
       pr3in.omegaS = params->OmegaS;
       pr3in.zeta2 = params->Zeta2;
       pr3in.omega = omega;
       rootIn.xmax = 300.;
       rootIn.xmin = 3.;
       LALDBisectionFindRoot(status->statusPtr, &r, &rootIn, (void *)&pr3in);
       CHECKSTATUSPTR(status);
       break;
     default:
       fprintf(stderr, 
	 		"There are no EOB waveforms implemented at order %d\n", 
	 		params->order);
       exit(0);
     }

   /*params->rInitial = r;
   params->vInitial = v;
   params->rLightRing = rn;
   */
/* 
   LALInspiralPhasing1(v) gives the GW phase (= twice the orbital phase).
   The ODEs we solve give the orbital phase. Therefore, set the
   initial phase to be half the GW pahse.
*/
   in3.totalmass = ak.totalmass;
   in3.dEnergy = func.dEnergy;
   in3.flux = func.flux;
   in3.coeffs = &ak;
   funcParams = (void *) &in3;

switch (params->order)
     {
     case twoPN:
     case twoPointFivePN:
       LALpphiInit(&q, r, eta);
       LALprInit(&p, r, &in3);
       in4.function = LALHCapDerivatives;
       break;
     case threePN:
     case threePointFivePN:
       LALpphiInit3PN(&q,r,eta, params->OmegaS);
       rootIn.function = LALprInit3PN;
       rootIn.xmax = 5;
       rootIn.xmin = -10;
       /* first we compute vr (we need coeef->Fp6) */
       pr3in.in3copy = in3; 
       pr3in.eta = eta;
       pr3in.omegaS = params->OmegaS;
       pr3in.zeta2 = params->Zeta2;
       pr3in.r = r;
       pr3in.q = q; 
       pr3in.omega = omega;
       LALvr3PN(&pr3in.vr,(void *) &pr3in);
       /* then we compute the initial value of p */
       LALDBisectionFindRoot(status->statusPtr, &p, &rootIn,(void *) &pr3in );
       CHECKSTATUSPTR(status);
       in4.function = LALHCapDerivatives3PN;
       break;
     default:
       fprintf(stderr, 
	 			"There are no EOB waveforms at order %d\n",
				params->order);
       exit(0);
     }

   values.data[0] = r;
   values.data[1] = s;
   values.data[2] = p;
   values.data[3] = q;


   in4.y = &values;
   in4.h = dt/m;
   in4.n = nn;
   in4.yt = &yt;
   in4.dym = &dym;
   in4.dyt = &dyt;

   t = 0.0;
   count = 0;
   while (count < params->nStartPad)
   {
      *(signal->data + count) = 0.;
      count++;
   }

   t = 0.0;
   /* Choose rOld larger than r now; it doesn't matter how large; just 
    * to kickoff the while statement below
    */
   rOld = r+0.1;
   dR = 0.0;
   while ((r-dR)>=rn && r<rOld) 
   {
      if (count >= (INT4)signal->length)
      {
         LALFree(dummy.data);
         ABORT(status, LALINSPIRALH_EVECTOR, LALINSPIRALH_MSGEVECTOR);
      }

      rOld = r;
      if (params->order<threePN){
		  LALHCapDerivatives(&values, &dvalues, funcParams);
		  }
      else  LALHCapDerivatives3PN(&values, &dvalues, funcParams);


      v = pow(dvalues.data[1], oneby3);
      h = params->signalAmplitude * v*v * cos(2.*s);
      *(signal->data + count) = (REAL4) h;

      omega = dvalues.data[1];
      in4.dydx = &dvalues;
      in4.x = t/m;
      LALRungeKutta4(status->statusPtr, &newvalues, &in4, funcParams);
      CHECKSTATUSPTR(status);
     
      r = values.data[0] = newvalues.data[0];
      s = values.data[1] = newvalues.data[1];
      p = values.data[2] = newvalues.data[2];
      q = values.data[3] = newvalues.data[3];

      dR = rOld-r;
      t = (++count-params->nStartPad) * dt;
/*----------------------------------------------------------
      printf("%e %e %e %e %e %e %e\n", t, r, v, s, p, q, dR);
      if (v>ak->vlso) printf("TLSO=%e\n", t);
      printf("&\n");
----------------------------------------------------------*/
   }  

/*----------------------------------------------------------------- 
Record the final cutoff frequency of BD Waveforms for record keeping 
-----------------------------------------------------------------*/
   /*params->rFinal = rOld;
   params->vFinal = v;*/
   params->fFinal = pow(v,3.)/(LAL_PI*m);
   params->tC = t;
   while (count < (INT4)signal->length) 
   {
       *(signal->data + count) = 0.;
       count++;
   }

   LALFree(dummy.data);
   DETATCHSTATUSPTR(status);
   RETURN(status);
}


NRCSID (LALEOBWAVEFORMTEMPLATESC,
"$Id$");

/*  <lalVerbatim file="LALEOBWaveformTemplatesCP"> */

void 
LALEOBWaveformTemplates (
   LALStatus        *status,
   REAL4Vector      *signal1,
   REAL4Vector      *signal2,
   InspiralTemplate *params
   ) 
{ /* </lalVerbatim> */

   INT4 count, nn=4;
   REAL8 amp, eta, m, rn, r, rOld, s, p, q, dt, t, h1, h2, v, omega, f;
   REAL8Vector dummy, values, dvalues, newvalues, yt, dym, dyt;
   TofVIn in1;
   InspiralPhaseIn in2;
   InspiralDerivativesIn in3;
   rk4In in4;
   pr3In pr3in; 
   void *funcParams;
   expnCoeffs ak;
   expnFunc func;
   rOfOmegaIn rofomegain;
   DFindRootIn rootIn;
   
   INITSTATUS(status, "LALEOBWaveformTemplates", LALEOBWAVEFORMTEMPLATESC);
   ATTATCHSTATUSPTR(status);

   ASSERT(signal1,  status,
	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signal2,  status,
   	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signal1->data,  status, 
   	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signal2->data,  status, 
   	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(params,  status, LALINSPIRALH_ENULL, 
   	LALINSPIRALH_MSGENULL);
   ASSERT(params->nStartPad >= 0, status, 
   	LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->nEndPad >= 0, status, 
   	LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->fLower > 0, status, 
   	LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->tSampling > 0, status,
   	LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->totalMass > 0., status, 
   	LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   LALInspiralSetup (status->statusPtr, &ak, params);
   CHECKSTATUSPTR(status);
   LALInspiralChooseModel(status->statusPtr, &func, &ak, params);
   CHECKSTATUSPTR(status);

   ASSERT(ak.totalmass/LAL_MTSUN_SI > 0.4, status, 
   	LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(ak.totalmass/LAL_MTSUN_SI < 100, status, 
   	LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

/* Allocate all the memory required to dummy and then point the various
   arrays to dummy - this makes it easier to handle memory failures */

   dummy.length = nn * 6;

   values.length = dvalues.length = newvalues.length =
   yt.length = dym.length = dyt.length = nn;

   if (!(dummy.data = (REAL8 * ) LALMalloc(sizeof(REAL8) * nn * 6))) {
      ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
   }

   values.data = &dummy.data[0];
   dvalues.data = &dummy.data[nn];
   newvalues.data = &dummy.data[2*nn];
   yt.data = &dummy.data[3*nn];
   dym.data = &dummy.data[4*nn];
   dyt.data = &dummy.data[5*nn];

   dt = 1./params->tSampling;
   eta = ak.eta;
   m = ak.totalmass;

   t = 0.0;
   in1.t = t;
   in1.t0 = ak.t0;
   in1.v0 = ak.v0;
   in1.vlso = ak.vlso;
   in1.totalmass = ak.totalmass;
   in1.dEnergy = func.dEnergy;
   in1.flux = func.flux;
   in1.coeffs = &ak;

   LALInspiralVelocity(status->statusPtr, &v, &in1);
   CHECKSTATUSPTR(status);

   omega = pow(v,3.);
   f = omega/(LAL_PI*m);

   in2.v0 = ak.v0;
   in2.phi0 = params->startPhase;
   in2.dEnergy = func.dEnergy;
   in2.flux = func.flux;
   in2.coeffs = &ak;
   LALInspiralPhasing1(status->statusPtr, &s, v, &in2);
   CHECKSTATUSPTR(status);
   s = s/2.;

   rofomegain.eta = eta;
   rofomegain.omega = omega;

   
   rootIn.xacc = 1.0e-16;
   switch (params->order)
     {
     case twoPN:
     case twoPointFivePN:
       rootIn.function = LALlightRingRadius;
       break;
     case threePN:
       rootIn.function = LALlightRingRadius3PN;
       break;
     default:
       fprintf(stderr, 
	 	"There are no EOB waveforms implemented at order %d\n", 
	 	params->order);
       exit(0);
     }
   rootIn.xmax = 2.;
   rootIn.xmin = 4.;
   funcParams = (void *) &rofomegain;

/*-------------------------------------------------------------------
Userful for debugging: Make sure a solution for r exists.
--------------------------------------------------------
   for (r=0.01; r<10.; r+=.01) {
      LALlightRingRadius(status->statusPtr, &x, r, funcParams);
      CHECKSTATUSPTR(status);
      printf("%e %e\n", r, x);
   }
   printf("&\n");
   for (r=1; r<100.; r+=.1) {
      LALrOfOmega(status->statusPtr, &x, r, funcParams);
      CHECKSTATUSPTR(status);
      printf("%e %e\n", r, x);
   }
-------------------------------------------------------------------*/

   LALDBisectionFindRoot(status->statusPtr, &rn, &rootIn, funcParams);
   CHECKSTATUSPTR(status);


    /*rofOmega */
    switch (params->order)
     {
     case twoPN:
     case twoPointFivePN:
       rootIn.function = LALrOfOmega;
       rootIn.xmax = 100.;
       rootIn.xmin = 6.;
       LALDBisectionFindRoot(status->statusPtr, &r, &rootIn, funcParams);
       CHECKSTATUSPTR(status);
       break;
     case threePN:
       rootIn.function = LALrOfOmega3PN;
       pr3in.eta = eta;
       pr3in.omegaS = params->OmegaS;
       pr3in.zeta2 = params->Zeta2;
       pr3in.omega = omega;
       rootIn.xmax = 100.;
       rootIn.xmin = 6.;
       LALDBisectionFindRoot(status->statusPtr, &r, &rootIn, (void *)&pr3in);
       CHECKSTATUSPTR(status);
       break;
     default:
       fprintf(stderr, 
	 	"There are no EOB waveforms implemented at order %d\n", 
	 	params->order);
       exit(0);
     }

   /*params->rInitial = r;
   params->vInitial = v;
   params->rLightRing = rn;
*/
/* 
   LALInspiralPhasing1(v) gives the GW phase (= twice the orbital phase).
   The ODEs we solve give the orbital phase. Therefore, set the
   initial phase to be half the GW pahse.
*/
   in3.totalmass = ak.totalmass;
   in3.dEnergy = func.dEnergy;
   in3.flux = func.flux;
   in3.coeffs = &ak;
   funcParams = (void *) &in3;

   switch (params->order)
     {
     case twoPN:
     case twoPointFivePN:
       LALpphiInit(&q, r, eta);
       LALprInit(&p, r, &in3);
       in4.function = LALHCapDerivatives;
       break;
     case threePN:
     case threePointFivePN:
       LALpphiInit3PN(&q,r,eta, params->OmegaS);
       rootIn.function = LALprInit3PN;
       rootIn.xmax = 5;
       rootIn.xmin = -10;
       /* first we compute vr (we need coeef->Fp6) */
       pr3in.in3copy = in3; 
       pr3in.eta = eta;
       pr3in.omegaS = params->OmegaS;
       pr3in.r = r;
       pr3in.q = q;
       pr3in.zeta2 = params->Zeta2;
       pr3in.omega = omega;
       LALvr3PN(&pr3in.vr,(void *) &pr3in);
       /* then we compute the initial value of p */
       LALDBisectionFindRoot(status->statusPtr, &p, &rootIn,(void *) &pr3in );
       CHECKSTATUSPTR(status);
       in4.function = LALHCapDerivatives3PN;
       break;
     default:
       fprintf(stderr, 
	 	"There are no EOB waveforms at order %d\n", 
	 	params->order);
       exit(0);
     }

   values.data[0] = r;
   values.data[1] = s;
   values.data[2] = p;
   values.data[3] = q;

   
   in4.function = LALHCapDerivatives;
   in4.y = &values;
   in4.h = dt/m;
   in4.n = nn;
   in4.yt = &yt;
   in4.dym = &dym;
   in4.dyt = &dyt;

   t = 0.0;
   count = 0;
   while (count < params->nStartPad) 
   {
      *(signal1->data + count) = *(signal2->data + count) = 0.;
      count++;
   }

   t = 0.0;
   rOld = r+0.1;
   while (r>=rn && r<rOld) {
      if (count >= (INT4)signal1->length) 
      {
	LALFree(dummy.data);
	ABORT(status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
      }

      rOld = r;
      
     if (params->order<threePN) {
		 LALHCapDerivatives(&values, &dvalues, funcParams);
	 }
     else {  
		 LALHCapDerivatives3PN(&values, &dvalues, funcParams);
	 }
      CHECKSTATUSPTR(status);


      v = pow(omega, oneby3);
      amp = params->signalAmplitude *v*v;
      h1 = amp * cos(2.*s);
      h2 = amp * cos(2.*s + LAL_PI_2);
      *(signal1->data + count) = (REAL4) h1;
      *(signal2->data + count) = (REAL4) h2;
      
      omega = dvalues.data[1];
      in4.dydx = &dvalues;
      in4.x = t/m;
      LALRungeKutta4(status->statusPtr, &newvalues, &in4, funcParams);
      CHECKSTATUSPTR(status);


      r = values.data[0] = newvalues.data[0];
      s = values.data[1] = newvalues.data[1];
      p = values.data[2] = newvalues.data[2];
      q = values.data[3] = newvalues.data[3];

      t = (++count-params->nStartPad) * dt;
/*----------------------------------------------------------
      printf("%e %e %e %e %e %e %e\n", t, r, v, s, p, q, h);
      if (v>ak->vlso) printf("TLSO=%e\n", t);
      printf("&\n");
----------------------------------------------------------*/
   }  

/*----------------------------------------------------------------- 
Record the final cutoff frequency of BD Waveforms for record keeping 
-----------------------------------------------------------------*/
   /*params->rFinal = rOld;
   params->vFinal = v;*/
   params->fFinal = pow(v,3.)/(LAL_PI*m);
   while (count < (INT4)signal1->length) 
   {
      *(signal1->data + count) = *(signal2->data + count) = 0.;
      count++;
   }

   LALFree(dummy.data);
   DETATCHSTATUSPTR(status);
   RETURN(status);
}


/*=========================================================*/
/*======INJECTION =========================================*/
/*=========================================================*/

/*  <lalVerbatim file="LALEOBWaveformForInjectionCP"> */
void 
LALEOBWaveformForInjection (
			    LALStatus        *status,
			    CoherentGW       *waveform,
			    InspiralTemplate *params,
			    PPNParamStruc    *ppnParams
			    ) 
{ 
  /* </lalVerbatim> */
  UINT4 count, nn=4, i;
  REAL8  eta, m, rn, r, rOld, s, p, q, dt, t,  v, omega ;
  REAL8Vector dummy, values, dvalues, newvalues, yt, dym, dyt;

  REAL4Vector *a=NULL;/* pointers to generated amplitude  data */
  REAL4Vector *ff=NULL ;/* pointers to generated  frequency data */
  REAL8Vector *phi=NULL;/* pointer to generated phase data */

  CreateVectorSequenceIn in;
  TofVIn in1;
  InspiralPhaseIn in2;
  InspiralDerivativesIn in3;
  rk4In in4;
  pr3In pr3in; 
  void *funcParams;
  expnCoeffs ak;
  expnFunc func;
  rOfOmegaIn rofomegain;
  DFindRootIn rootIn;
  /* 2 Additional lines of code to convert the EOB frequency 
     into physical frequency for the inject package */
  REAL8 unitHz;
  REAL8 f2a;
  REAL8 mu; 
  REAL8 mTot;
  REAL8 cosI;/* cosine of system inclination */
  REAL8 etab;
  REAL8 fFac; /* SI normalization for f and t */
  REAL8 f2aFac;/* factor multiplying f in amplitude function */
  REAL8 apFac, acFac;/* extra factor in plus and cross amplitudes */
  REAL8 phiC;/* phase at coalescence */
  
 CHAR message[256];
  InspiralInit paramsInit;   
    
  INITSTATUS(status, "LALEOBWaveformForInjection", LALEOBWAVEFORMTEMPLATESC);
  ATTATCHSTATUSPTR(status);


  /* Make sure parameter and waveform structures exist. */
  ASSERT( params, status, LALINSPIRALH_ENULL,  LALINSPIRALH_MSGENULL );
  ASSERT(waveform, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);  
  /* Make sure waveform fields don't exist. */
  ASSERT( !( waveform->a ), status, 
  	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->f ), status, 
  	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->phi ), status, 
  	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->shift ), status, 
  	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );

  /* Compute some parameters*/
  LALInspiralInit(status->statusPtr, params, &paramsInit);
  CHECKSTATUSPTR(status);   
  if (paramsInit.nbins==0)
    {
      DETATCHSTATUSPTR(status);
      RETURN (status);      
    }
  func = paramsInit.func;
  ak   = paramsInit.ak;
  
  dt = 1./params->tSampling;
  

  mTot   =  params->mass1 + params->mass2;
  etab   =  params->mass1 * params->mass2;
  etab  /= mTot;
  etab  /= mTot;
  unitHz = (mTot) *LAL_MTSUN_SI*(REAL8)LAL_PI;
  cosI   = cos( params->inclination );
  mu     = etab * mTot;  
  fFac   = 1.0 / ( 4.0*LAL_TWOPI*LAL_MTSUN_SI*mTot );
  f2aFac = LAL_PI*LAL_MTSUN_SI*mTot*fFac;   
  apFac  = acFac = -2.0 * mu * LAL_MRSUN_SI/params->distance;
  apFac *= 1.0 + cosI*cosI;
  acFac *= 2.0*cosI;
    params->nStartPad = 0;

  /* Now we can allocate memory and vector for coherentGW structure*/     
  LALSCreateVector(status->statusPtr, &ff, paramsInit.nbins);
  CHECKSTATUSPTR(status);
  LALSCreateVector(status->statusPtr, &a, 2*paramsInit.nbins);
  CHECKSTATUSPTR(status);
  LALDCreateVector(status->statusPtr, &phi, paramsInit.nbins);
  CHECKSTATUSPTR(status);
  
   /* By default the waveform is empty */
  for (count = 0; count < paramsInit.nbins; count++) 
    {
      ff->data[count]           = 0.;
      a->data[2*count+1]        = 0.;
      phi->data[count]          = 0.;
      a->data[2*count]          = 0.;
    }
  count = 0;


 /* Allocate all the memory required to dummy and then point the various
     arrays to dummy - this makes it easier to handle memory failures */
  dummy.length = nn * 6;
  values.length 	= dvalues.length	= newvalues.length = nn;
  yt.length 		= dym.length		= dyt.length = nn;
  
  if (!(dummy.data = (REAL8 * ) LALMalloc(sizeof(REAL8) * nn * 6))) {
    ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
  }
  
  values.data    = &dummy.data[0];
  dvalues.data   = &dummy.data[nn];
  newvalues.data = &dummy.data[2*nn];
  yt.data    = &dummy.data[3*nn];
  dym.data   = &dummy.data[4*nn];
  dyt.data   = &dummy.data[5*nn];
  
  /* -- Some aliases -- */
  dt 	= 1./params->tSampling;
  eta	= ak.eta;
  m 	= ak.totalmass;
  
 
  /* let's find the initial velocity first fiven the lower frequency */
  t 		= 0.0;
  in1.t 	= t;
  in1.t0 	= ak.t0;
  in1.v0 	= ak.v0;
  in1.vlso 	= ak.vlso;
  in1.totalmass = ak.totalmass;
  in1.dEnergy 	= func.dEnergy;
  in1.flux 	= func.flux;
  in1.coeffs 	= &ak;
  
  LALInspiralVelocity(status->statusPtr, &v, &in1);
  CHECKSTATUSPTR(status);
  
  omega = pow(v,3.);
  
  /* -- then the initial phase -- */
  in2.v0 	= ak.v0;
  in2.phi0 	= params->startPhase;
  in2.dEnergy 	= func.dEnergy;
  in2.flux	= func.flux;
  in2.coeffs 	= &ak;
  LALInspiralPhasing1(status->statusPtr, &s, v, &in2);
  CHECKSTATUSPTR(status);
  s = s/2.;

  /* -- as well as the light ring value where to stop the waveform -- */
  rofomegain.eta = eta;
  rofomegain.omega = omega;
  
  
  rootIn.xacc = 1.0e-16;
  switch (params->order)
    {
    case twoPN:
    case twoPointFivePN:
      rootIn.function = LALlightRingRadius;
      break;
    case threePN:
      rootIn.function = LALlightRingRadius3PN;
      break;
    default:
      fprintf(stderr, 
		"There are no EOB waveforms implemented at order %d\n",
		params->order);
      exit(0);
    }
  rootIn.xmax = 2.;
  rootIn.xmin = 4.;
  funcParams = (void *) &rofomegain;

  LALDBisectionFindRoot(status->statusPtr, &rn, &rootIn, funcParams);
  CHECKSTATUSPTR(status);
  
 /* -- Finally, let's find the initial radius -- 
   * Here the code might crash in the BisectionFindRoot. First we have 
	to check that the values of the brackets are negative for the minimal value
	of r and positive for its maximum value. */

  /* For information, rofOmega of a (0.2, 0.2)Mo system with a Flower =
	 10Hz is equal to 459 therefore rootIn.xmax =1000 is large enough for most 
	 of the  systems we deal with . cokelaer, May 2004. */
  
  /* rootIn.min = 3  which corresponds to the local minima of the 
	 function. Below,  the function might return nan values. That's why we 
	 use the finite function of the standard math library.*/
	 
  rootIn.function = LALrOfOmega;
  rootIn.xmax 	= 1000.;
  rootIn.xmin 	= 3.;
  
  /* No check status here. We check manually since we want to get the
   *  hand back from those two calls.
   *  First, let's check what's happen with root.xmax */
  LALrOfOmega(status->statusPtr, &r, rootIn.xmax, funcParams);
  sprintf(message, 
	 	"EOB:initialCondition:For omega=rootIn.xmax=%f, LALrOfOmega return r=%f", 
	 	rootIn.xmax, r);
  LALInfo(status->statusPtr, message);
  if (r < 0 || finite(r)==0 ) {   
    sprintf(message, LALINSPIRALH_MSGEROOTINIT);
    LALWarning(status->statusPtr, message);
    sprintf(message, LALINSPIRALH_MSGESTOPPED);
    LALWarning(status->statusPtr, message);		
    
    LALSDestroyVector(status->statusPtr, &ff);
    CHECKSTATUSPTR(status);
    LALSDestroyVector(status->statusPtr, &a);
    CHECKSTATUSPTR(status);
    LALDDestroyVector(status->statusPtr, &phi);
    CHECKSTATUSPTR(status);
    
    LALFree(dummy.data);
    DETATCHSTATUSPTR(status);
    RETURN(status);
  }
  
  /* Then, let's check what's happen  with rootIn.xmin*/
  LALrOfOmega(status->statusPtr, &r, rootIn.xmin, funcParams);
  sprintf(message, 
  	"EOB:initialCondition:For omega=rootIn.xmin=%f, LALrOfOmega return r=%f", 
  	rootIn.xmin, r);
  
  LALInfo(status->statusPtr, message);
  if (r > 0 || finite(r)==0){
    sprintf(message, LALINSPIRALH_MSGEROOTINIT);
    LALWarning(status->statusPtr, message);		
    sprintf(message, LALINSPIRALH_MSGESTOPPED);
    LALWarning(status->statusPtr, message);	
    
    LALSDestroyVector(status->statusPtr, &ff);
    CHECKSTATUSPTR(status);
    LALSDestroyVector(status->statusPtr, &a);
    CHECKSTATUSPTR(status);
    LALDDestroyVector(status->statusPtr, &phi);
    CHECKSTATUSPTR(status);
    
    LALFree(dummy.data);
    DETATCHSTATUSPTR(status);
    RETURN(status);
  }


  /* -- the inital radius might be found since the brackets seem 
   * to give a positive  and a negative value. let's try to find it */   
  switch (params->order)
    {
    case twoPN:
    case twoPointFivePN:       
      funcParams 	= (void *) &rofomegain;
      rootIn.function 	= LALrOfOmega;
      rootIn.xmax 	= 1000.; /*should be add in the header*/
      rootIn.xmin 	= 3.;    /* should be add in the header */
      
      LALDBisectionFindRoot(status->statusPtr, &r, &rootIn, funcParams);
      CHECKSTATUSPTR(status);
      break;
    case threePN:
      rootIn.function 	= LALrOfOmega3PN;
      pr3in.eta 	= eta;
      pr3in.omegaS 	= params->OmegaS;
      pr3in.zeta2 	= params->Zeta2;
      pr3in.omega 	= omega;
      rootIn.xmax 	= 1000.;
      rootIn.xmin 	= 3.;
      
      LALDBisectionFindRoot(status->statusPtr, &r, &rootIn, (void *)&pr3in);
      CHECKSTATUSPTR(status);
	
    
      break;
    default: /*should not reach that point. it should be 
				checked before I guess*/
      fprintf(stderr, 
		"There are no EOB waveforms implemented at order %d\n", 
		params->order);
      exit(0);
    }
  
  sprintf(message, "EOB:initialCondition:Initial r found = %f\n", r);
  LALInfo(status->statusPtr, message); 

        
  if (r < 6){
    sprintf(message, 
	  	"EOB:initialCondition:Initial r found = %f ",r);
	sprintf(message, 	
	  	"too small (below 6 no waveform is generated)\n");
    LALWarning(status->statusPtr, message); 
    
    LALSDestroyVector(status->statusPtr, &ff);
    CHECKSTATUSPTR(status);
    LALSDestroyVector(status->statusPtr, &a);
    CHECKSTATUSPTR(status);
    LALDDestroyVector(status->statusPtr, &phi);
    CHECKSTATUSPTR(status);
    
    LALFree(dummy.data);
    DETATCHSTATUSPTR(status);
    RETURN(status);    
  }


  /* -- let's find the initial momentum -- */
  /* LALInspiralPhasing1(v) gives the GW phase (twice the orbital phase)
     The ODEs we solve give the orbital phase. Therefore, set the
     initial phase to be half the GW pahse.   */
  in3.totalmass 	= ak.totalmass;
  in3.dEnergy		= func.dEnergy;
  in3.flux 		= func.flux;
  in3.coeffs		= &ak;
  funcParams 		= (void *) &in3;
  
  switch (params->order)
    {
    case twoPN:
    case twoPointFivePN:
      LALpphiInit(&q, r, eta);
      LALprInit(&p, r, &in3);       
      in4.function = LALHCapDerivatives;
      break;
    case threePN:
    case threePointFivePN:
      LALpphiInit3PN(&q, r, eta, params->OmegaS);
      rootIn.function 	= LALprInit3PN;
      rootIn.xmax 	= 5;  /* have to check validity of those brackets.*/
      rootIn.xmin 	= -10; /* idem */
      /* first we compute vr (we need coeff->Fp6) */
      pr3in.in3copy 	= in3; 
      pr3in.eta 	= eta;
      pr3in.omegaS 	= params->OmegaS;
      pr3in.zeta2 	= params->Zeta2;
      pr3in.r		= r;
      pr3in.q		= q; 
      pr3in.omega 	= omega;
      
      LALvr3PN(&pr3in.vr, (void *) &pr3in);
      /* then we compute the initial value of p */
      LALDBisectionFindRoot(status->statusPtr, &p, &rootIn,(void *) &pr3in );
      CHECKSTATUSPTR(status);
      
      in4.function = LALHCapDerivatives3PN;
      break;
    default:
      break;
    }
  
  sprintf( message, 
		"rInitial= %f, vinitial= %f, lightring= %f, initial phase= %f, pphase%f, pr=%f ", 
	   r, v , rn, s, q, p);
  LALInfo(status->statusPtr, message);
  
  /** -------------------------------------------------------------- */
  /** -- actual computation started here -- */   
  values.data[0] = r;
  values.data[1] = s;
  values.data[2] = p;
  values.data[3] = q;
  
  in4.y 	= &values;
  in4.h 	= dt/m;
  in4.n 	= nn;
  in4.yt 	= &yt;
  in4.dym 	= &dym;
  in4.dyt 	= &dyt;
  
  t     = 0.0;
  count = 0;
  rOld  = r+0.1;
  
  
  do
    {
      if (count >= (INT4)ff->length)
      {
	LALSDestroyVector(status->statusPtr, &ff);
	CHECKSTATUSPTR(status);
	LALSDestroyVector(status->statusPtr, &a);
	CHECKSTATUSPTR(status);
	LALDDestroyVector(status->statusPtr, &phi);
	CHECKSTATUSPTR(status);
 	LALFree(dummy.data);
	ABORT(status, LALINSPIRALH_EVECTOR, LALINSPIRALH_MSGEVECTOR);
      }

      rOld = r;
      
      if (params->order<threePN) {
	LALHCapDerivatives(&values, &dvalues, funcParams);
      }
      else {
	LALHCapDerivatives3PN(&values, &dvalues, funcParams);
      }
      
      v = pow(omega, oneby3);
      
      /*    ff->data[count]           = (REAL4)  (omega * fFac ); */
      ff->data[count]= (REAL4)(omega/unitHz);
      f2a = pow (f2aFac * omega, 2./3.);
      a->data[2*count]          = (REAL4)(4.*apFac * f2a);
      a->data[2*count+1]        = (REAL4)(4.*acFac * f2a);
      phi->data[count]          = (REAL8)(2* s);
      
      omega = dvalues.data[1];
      in4.dydx = &dvalues;
      in4.x = t/m;
      LALRungeKutta4(status->statusPtr, &newvalues, &in4, funcParams);
      CHECKSTATUSPTR(status);
      
      r = values.data[0] = newvalues.data[0];
      s = values.data[1] = newvalues.data[1];
      p = values.data[2] = newvalues.data[2];
      q = values.data[3] = newvalues.data[3];
      
      count++;
    }
  while (r>=rn && r<rOld) ;
    
  
  /*----------------------------------------------------------------- 
    Record the final cutoff frequency of BD Waveforms for record keeping 
    -----------------------------------------------------------------*/
  /*params->rFinal = rOld;
  params->vFinal = v;*/
  params->fFinal = pow(v,3.)/(LAL_PI*m);
sprintf(message, "fFinal = %f", params->fFinal);
  LALInfo(status, message);

 sprintf(message, "cycles = %f", s/3.14159);
  LALInfo(status, message);

  sprintf( message, "final coalescence phase with respet to actual data =%f ",
  	(ff->data[count]-ff->data[count-1])/2/3.14159); 
  LALInfo(status->statusPtr, message);
  


  if ( (s/LAL_PI) < 2 ){
    sprintf(message, "The waveform has only %f cycles; we don't keep waveform with less than 2 cycles.", 
	      (double) s/ (double)LAL_PI );	
    LALWarning(status, message);
  }
  else
    {
      phiC =  phi->data[count-1] ;
      
      for (i=0; i<count;i++)
	{
	  phi->data[i] =  -phiC + phi->data[i] +ppnParams->phi;
	}
      
      /* Allocate the waveform structures. */
      if ( ( waveform->a = (REAL4TimeVectorSeries *)
	     LALMalloc( sizeof(REAL4TimeVectorSeries) ) ) == NULL ) {
	ABORT( status, LALINSPIRALH_EMEM,
	       LALINSPIRALH_MSGEMEM );
      }
      memset( waveform->a, 0, sizeof(REAL4TimeVectorSeries) );
      if ( ( waveform->f = (REAL4TimeSeries *)
	     LALMalloc( sizeof(REAL4TimeSeries) ) ) == NULL ) {
	LALFree( waveform->a ); waveform->a = NULL;
	ABORT( status, LALINSPIRALH_EMEM,
	       LALINSPIRALH_MSGEMEM );
      }
      memset( waveform->f, 0, sizeof(REAL4TimeSeries) );
      if ( ( waveform->phi = (REAL8TimeSeries *)
	     LALMalloc( sizeof(REAL8TimeSeries) ) ) == NULL ) {
	LALFree( waveform->a ); waveform->a = NULL;
	LALFree( waveform->f ); waveform->f = NULL;
	ABORT( status, LALINSPIRALH_EMEM,
	       LALINSPIRALH_MSGEMEM );
      }
      memset( waveform->phi, 0, sizeof(REAL8TimeSeries) );
      
      
      
      in.length = (UINT4)count;
      in.vectorLength = 2;
      LALSCreateVectorSequence( status->statusPtr,
				&( waveform->a->data ), &in );
      CHECKSTATUSPTR(status);      
      LALSCreateVector( status->statusPtr,
			&( waveform->f->data ), count);
      CHECKSTATUSPTR(status);      
      LALDCreateVector( status->statusPtr,
			&( waveform->phi->data ), count );
      CHECKSTATUSPTR(status);        
      
      
      
      memcpy(waveform->f->data->data , ff->data, count*(sizeof(REAL4)));
      memcpy(waveform->a->data->data , a->data, 2*count*(sizeof(REAL4)));
      memcpy(waveform->phi->data->data ,phi->data, count*(sizeof(REAL8)));
      
      
      waveform->a->deltaT = waveform->f->deltaT = waveform->phi->deltaT
	= 1./params->tSampling;
      
      waveform->a->sampleUnits = lalStrainUnit;
      waveform->f->sampleUnits = lalHertzUnit;
      waveform->phi->sampleUnits = lalDimensionlessUnit;
      waveform->position = ppnParams->position;
      waveform->psi = ppnParams->psi;
      
      
      LALSnprintf( waveform->a->name, 
	  	LALNameLength, "EOB inspiral amplitudes");
      LALSnprintf( waveform->f->name,
		  LALNameLength, "EOB inspiral frequency");
      LALSnprintf( waveform->phi->name, 
	  	LALNameLength, "EOB inspiral phase");
      
      /* --- fill some output ---*/
      ppnParams->tc     = (double)(count-1) / params->tSampling ;
      ppnParams->length = count;
      ppnParams->dfdt   = ((REAL4)(waveform->f->data->data[count-1] 
				   - waveform->f->data->data[count-2]))
	* ppnParams->deltaT;
      ppnParams->fStop  = params->fFinal;
      ppnParams->termCode        = GENERATEPPNINSPIRALH_EFSTOP;
      ppnParams->termDescription = GENERATEPPNINSPIRALH_MSGEFSTOP;
      
      ppnParams->fStart   = ppnParams->fStartIn;
    } /* end phase condition*/

  /* --- free memory --- */


  LALSDestroyVector(status->statusPtr, &ff);
  CHECKSTATUSPTR(status);
  LALSDestroyVector(status->statusPtr, &a);
  CHECKSTATUSPTR(status);
  LALDDestroyVector(status->statusPtr, &phi);
  CHECKSTATUSPTR(status);

 
  /*on peut utiliser tSampling pour dfdt*/
  
   LALFree(dummy.data);
   DETATCHSTATUSPTR(status);
   RETURN(status);
}
