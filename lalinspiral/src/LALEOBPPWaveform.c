/*
*  Copyright (C) 2007 Stas Babak, David Churches, Duncan Brown, David Chin, Jolien Creighton,
*                     B.S. Sathyaprakash, Craig Robinson , Thomas Cokelaer, Evan Ochsner
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/*  <lalVerbatim file="LALEOBPPWaveformCV">
Author: Sathyaprakash, B. S., Cokelaer T.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALEOBPPWaveform.c} and
\texttt{LALEOBPPWaveformTemplates.c}}

Module to generate effective-one-body waveforms.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALEOBPPWaveformCP}
\index{\verb&LALEOBPPWaveform()&}
\begin{itemize}
\item {\tt signalvec:} Output containing the inspiral waveform.
\item {\tt params:} Input containing binary chirp parameters.
\end{itemize}

\input{LALEOBPPWaveformTemplatesCP}
\index{\verb&LALEOBPPWaveformTemplates()&}
\begin{itemize}
\item {\tt signalvec1:} Output containing the 0-phase inspiral waveform.
\item {\tt signalvec2:} Output containing the $\pi/2$-phase inspiral waveform.
\item {\tt params:} Input containing binary chirp parameters.
\end{itemize}

\input{LALEOBPPWaveformForInjectionCP}
\index{\verb&LALEOBPPWaveformForInjection()&}
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
   LALHCapDerivativesP4PN
   LALlightRingRadius
   LALlightRingRadius3PN
   LALlightRingRadiusP4PN
   LALpphiInit
   LALpphiInit3PN
   LALpphiInitP4PN
   LALprInit
   LALprInit3PN
   LALprInitP4PN
   LALrOfOmega
   LALrOfOmega3PN
   LALrOfOmegaP4PN
\end{verbatim}

\subsubsection*{Notes}
The length of the waveform returned by {\tt LALInspiralWaveLength} is
occassionally smaller than what is required to hold an EOB waveform.
This is because EOB goes beyond the last stable orbit up to the light
ring while {\tt LALInspiralWaveLength} assumes that the waveform terminates
at the last stable orbit. It is recommended that a rather generous
{\tt params->nEndPad} be used to prevent the code from crashing.

\vfill{\footnotesize\input{LALEOBPPWaveformCV}}

</lalLaTeX>  */
#include <lal/Units.h>
#include <lal/LALInspiral.h>
#include <lal/FindRoot.h>
#include <lal/SeqFactories.h>
#include <lal/NRWaveInject.h>
#include <lal/LALComplex.h>

#include <gsl/gsl_sf_gamma.h>

typedef struct tagrOfOmegaIn {
   REAL8 eta, omega;
} rOfOmegaIn;

typedef struct tagPr3In {
  REAL8 eta, zeta2, omegaS, omega, vr,r,q;
  InspiralDerivativesIn in3copy;
} pr3In;


static inline REAL8
XLALCalculateA5( REAL8 eta );

static void
omegaofrP4PN (
             REAL8 *x,
             REAL8 r,
             void *params) ;


static
void LALHCapDerivativesP4PN(   REAL8Vector *values,
                               REAL8Vector *dvalues,
                               void        *funcParams);

static
void LALHCapDerivativesP4PNFF(   REAL8Vector *values,
                                 REAL8Vector *dvalues,
                                 void        *funcParams);

static
REAL8 XLALCalculateEOBA( REAL8    r,
                         REAL8  eta);
static
REAL8 XLALCalculateEOBD( REAL8    r,
                         REAL8	eta);
static
REAL8 XLALEffectiveHamiltonian( REAL8Vector           *values,
                                InspiralDerivativesIn *ak);
static
void LALprInitP4PN(LALStatus *status, REAL8 *pr , REAL8 , void  *params);

static
void LALpphiInitP4PN(REAL8 *phase, REAL8 r, REAL8 eta, REAL8 omegaS);

static
void LALlightRingRadiusP4PN(LALStatus *status, REAL8 *x, REAL8 r, void *params);

static
REAL8 XLALrOfOmegaP4PN (REAL8 r, void *params);

static
void LALvrP4PN(REAL8 *vr, void *params);

static void
LALEOBPPWaveformEngine (
                LALStatus        *status,
                REAL4Vector      *signalvec1,
                REAL4Vector      *signalvec2,
                REAL4Vector      *h,
                REAL4Vector      *a,
                REAL4Vector      *ff,
                REAL8Vector      *phi,
                UINT4            *countback,
                InspiralTemplate *params,
                InspiralInit     *paramsInit
                );

NRCSID (LALEOBPPWAVEFORMC,
"$Id$");


INT4 XLALGetFactorizedWaveform( COMPLEX16		*hlm,
				REAL8Vector           *values,
				REAL8Vector		  *dvalues,
                                InspiralDerivativesIn *ak,
                                const INT4            l,
                                const INT4            m
                                )
{
	static const char func[] = "XLALGetFactorizedWaveform";

        /* Status of function calls */
        INT4 status;
	
	REAL8 eta, eta2, eta3, dM, chiS, chiA, a, a2, a3;
	REAL8 r, pr, pp, Omega, v, v2, vh, vh3, k, hathatk, eulerlogxabs;
	REAL8 NQC, Hreal, Heff, Slm, deltalm, rholm;
        REAL8 aNQC1, aNQC2, aNQC3;
	COMPLEX16 Tlm;
        COMPLEX16 hNewton;
	gsl_sf_result lnr1, arg1, lnr2, arg2;

        /* Non-Keplerian velocity */
        REAL8 vPhi;

        pr3In  pr3in;
        

	if ( abs(m) > (INT4) l )
	{
	  XLAL_ERROR( func, XLAL_EINVAL );
	}
	
	eta	= ak->coeffs->eta;
	eta2	= eta * eta;
	eta3	= eta * eta2;
	dM	= sqrt(1 - 4*eta);
	
	chiS	= 0.0;
	chiA	= 0.0;
	a	= 0.0;
	a2	= a * a;
	a3	= a * a2;

	r	= values -> data[0];
        pr	= values -> data[2];
	pp	= values -> data[3];

	Heff	= XLALEffectiveHamiltonian( values, ak ); 
	Hreal	= sqrt( 1.0 + 2.0 * eta * ( Heff - 1.0) );
	Omega	= dvalues -> data[1];
	v	= pow( Omega, 1./3.);
	v2	= v * v;
	vh	= pow( Hreal * Omega, 1./3.);
	vh3	= vh * vh * vh;
	eulerlogxabs = LAL_GAMMA + log( 2*m ) + log( v );
	
        /* Calculate the non-Keplerian velocity */
        memset( &pr3in, 0, sizeof( pr3In ) );
        pr3in.omega  = Omega;
        pr3in.omegaS = ak->coeffs->omegaS;
        pr3in.zeta2  = ak->coeffs->zeta2;
        pr3in.eta    = eta;
  
        omegaofrP4PN( &vPhi, values->data[0], &pr3in );

        vPhi  = pow( vPhi, - twoby3 );
        vPhi *= Omega;

        /* Calculate the NQC term. aNQC values will be given after calibration. */
	aNQC1 = 0.0;
	aNQC2 = 0.0;
	aNQC3 = 0.0;
        NQC = 1.0 + pr*pr/r/r/Omega/Omega*( aNQC1 + aNQC2 / sqrt(r) + aNQC3 / r );
        
        /* Calculate the newtonian multipole */
        status = XLALCalculateNewtonianMultipole( &hNewton, vPhi * vPhi,
                         values->data[1], (UINT4)l, m, ak );
        if ( status == XLAL_FAILURE )
        {
          XLAL_ERROR( func, XLAL_EFUNC );
        }

        /* Calculate the source term */
	if ( ( (l+m)%2 ) == 0)
	{ 
	  Slm = Heff;
	}
	else
	{
	  Slm = v * pp;
	}

        /* Calculate the Tail term */	
	k	= m * Omega;
	hathatk = Hreal * k;
	XLAL_CALLGSL( status = gsl_sf_lngamma_complex_e( l+1.0, -2.0*hathatk, &lnr1, &arg1 ) );
	if (status != GSL_SUCCESS)
	{
	  XLALPrintError("Error in GSL function\n" );
	  XLAL_ERROR( func, XLAL_EFUNC );
	}
	XLAL_CALLGSL( status = gsl_sf_lngamma_complex_e( l+1.0, 0.0, &lnr2, &arg2 ) );
	if ( status != GSL_SUCCESS)
	{
	  XLALPrintError("Error in GSL function\n" );
	  XLAL_ERROR( func, XLAL_EFUNC );
	}
	Tlm = XLALCOMPLEX16Exp( XLALCOMPLEX16Rect( lnr1.val - lnr2.val + LAL_PI * hathatk, 
				arg1.val - arg2.val + 2.0 * hathatk * log(4.0*k/sqrt(LAL_E)) ) );

        /* Calculate the residue phase and amplitude terms */
	switch( l )
	{
	  case 2:
	    switch( abs(m) )
	    {
	      case 2:
	        deltalm = vh3*(7./3. + vh3*((-4*a)/3. + (428*LAL_PI)/105. 
			+ vh*vh*((20*a)/63. + (-2203./81. + (1712*LAL_PI*LAL_PI)/315.)*vh))) 
			- 24*eta*v*v2*v2;
		rholm	= 1 + v2*(-43./42. + (55*eta)/84. + v*((-2*(chiS + chiA*dM 
			- chiS*eta))/3. + v*(-20555./10584. + (chiS*chiS + 2*chiA*chiS*dM 
			+ chiA*chiA*(1 - 4*eta))/2. - (33025*eta)/21168. + (19583*eta2)/42336. 
			+ v*((-34*a)/21. + v*(1556919113./122245200. + (89*a2)/252. 
			- (48993925*eta)/9779616. - (6292061*eta2)/3259872. 
			+ (10620745*eta3)/39118464. - (428*eulerlogxabs)/105. 
			+ (41*eta*LAL_PI*LAL_PI)/192. + v*((18733*a)/15876. + a*a2/3. 
			+ v*(-387216563023./160190110080. + (18353*a2)/21168. - a2*a2/8. 
			+ (9202*eulerlogxabs)/2205. + (-16094530514677./533967033600. 
			+ (439877*eulerlogxabs)/55566.)*v2)))))));
	        break;
	      case 1:
	        deltalm = vh3*(2./3. + vh3*((-17*a)/35. + (107*LAL_PI)/105. 
			+ vh*((3*a2)/140. + (-272./81. + (214*LAL_PI*LAL_PI)/315.)*vh*vh))) 
			- (493*eta*v*v2*v2)/42.;
		rholm	= 1 + v*((-3*(chiA + chiS*dM))/(4.*dM) 
			+ v*(-59./56 - (9*pow(chiA + chiS*dM,2))/(32.*(1 - 4*eta)) 
			+ (23*eta)/84. + v*((-567*chiA*chiA*chiA - 1701*chiA*chiA*chiS*dM 
			+ chiA*(-4708 + 1701*chiS*chiS - 2648*eta)*(-1 + 4*eta) 
			+ chiS*pow(1 - 4*eta,1.5)*(4708 - 567*chiS*chiS 
			+ 1816*eta))/(2688.*pow(1 - 4*eta,1.5)) + v*(-47009./56448. 
			- (865*a2)/1792. - (405*a2*a2)/2048. - (10993*eta)/14112. 
			+ (617*eta2)/4704. + v*((-98635*a)/75264. + (2031*a*a2)/7168. 
			- (1701*a2*a3)/8192. + v*(7613184941./2607897600. 
			+ (9032393*a2)/1806336. + (3897*a2*a2)/16384. 
			- (15309*a3*a3)/65536. - (107*eulerlogxabs)/105. 
			+ v*((-3859374457*a)/1159065600. - (55169*a3)/16384. 
			+ (18603*a2*a3)/65536. - (72171*a2*a2*a3)/262144. 
			+ (107*a*eulerlogxabs)/140. + v*(-1168617463883./911303737344. 
			+ (6313*eulerlogxabs)/5880. + (-63735873771463./16569158860800. 
			+ (5029963*eulerlogxabs)/5927040.)*v2))))))));
	        break;
	      default:
                XLAL_ERROR( func, XLAL_EINVAL );
                break;
	    }
	    break;
	  case 3:
	    switch (m)
	    {
	      case 3:
	        deltalm = vh3*(13./10. + vh3*((-81*a)/20. + (39*LAL_PI)/7. + (-227827./3000. 
			+ (78*LAL_PI*LAL_PI)/7.)*vh3)) - (80897*eta*v*v2*v2)/2430.;
		rholm	= 1 + v2*(-7./6. + (2*eta)/3. + v*((chiS*dM*(-4 + 5*eta) 
			+ chiA*(-4 + 19*eta))/(6.*dM) + v*(-6719./3960. + a2/2. 
			- (1861*eta)/990. + (149*eta2)/330. + v*((-4*a)/3. 
			+ v*(3203101567./227026800. + (5*a2)/36. - (26*eulerlogxabs)/7. 
			+ v*((5297*a)/2970. + a*a2/3. + (-57566572157./8562153600. 
			+ (13*eulerlogxabs)/3.)*v))))));
	        break;
	      case 2:
		deltalm = vh3*((10 + 33*eta)/(15.*(1 - 3*eta)) + vh*(4*a + vh*vh*((-136.*a)/45. 
			+ (52.*LAL_PI)/21. + (-9112./405. + (208.*LAL_PI*LAL_PI)/63.)*vh3)));
		rholm	= 1 + v*((4*chiS*eta)/(3.*(1 - 3*eta)) 
			+ v*((-4*a2*eta2)/(9.*(1.-3.*eta)*(1.-3.*eta)) + (328. - 1115.*eta 
			+ 320.*eta2)/(270.*(-1 + 3*eta)) + v*((2*(45*a*pow(-1 + 3*eta,3) 
			- a*eta*(328 - 2099*eta + 5*(733 + 20*a2)*eta2 
			- 960*eta3)))/(405.*pow(-1 + 3*eta,3)) + v*(a2/3. + (-1444528. 
			+ 8050045.*eta - 4725605.*eta2 - 20338960.*eta3 
			+ 3085640.*eta2*eta2)/(1603800.*(1.-3.*eta)*(1.-3.*eta)) + v*((-2788*a)/1215. 
			+ v*(5849948554./940355325. + (488*a2)/405. - (104*eulerlogxabs)/63. 
			+ (-10607269449358./3072140846775. + (17056.*eulerlogxabs)/8505.)*v2))))));
		break;
	      case 1:
		deltalm = vh3*(13./30. + vh3*((61*a)/20. + (13*LAL_PI)/21. 
			+ vh*((-24*a2)/5. + (-227827./81000. + (26*LAL_PI*LAL_PI)/63.)*vh*vh))) 
			- (17*eta*v*v2*v2)/10.;
		rholm	= 1 + v2*(-13./18. - (2*eta)/9. + v*((chiA*(-4 + 11*eta) 
			+ chiS*dM*(-4 + 13*eta))/(6.*dM) + v*(101./7128. 
			- (5*a2)/6. - (1685*eta)/1782. - (829*eta2)/1782. + v*((4*a)/9. 
			+ v*(11706720301./6129723600. - (49*a2)/108. - (26*eulerlogxabs)/63. 
			+ v*((-2579*a)/5346. + a*a2/9. + (2606097992581./4854741091200. 
			+ (169*eulerlogxabs)/567.)*v))))));
		break;
              default:
                XLAL_ERROR( func, XLAL_EINVAL );
                break;
	    }
	    break;
	  case 4:
	    switch (m)
	    {
	      case 4:
	        deltalm = vh3*((112 + 219*eta)/(120.*(1 - 3*eta)) + ((-464*a)/75. 
			+ (25136*LAL_PI)/3465.)*vh3);
		rholm	= 1 + v2*((1614 - 5870*eta + 2625*eta2)/(1320.*(-1 + 3*eta)) 
			+ v*((chiA*(10 - 39*eta)*dM + chiS*(10 - 41*eta 
			+ 42*eta2))/(-15 + 45*eta) + v*(a2/2. + (-511573572 
			+ 2338945704*eta - 313857376*eta2 - 6733146000*eta3 
			+ 1252563795*eta2*eta2)/(317116800.*(1.-3.*eta)*(1.-3.*eta)) 
			+ v*((-69*a)/55. + (16600939332793./1098809712000. + (217*a2)/3960. 
			- (12568*eulerlogxabs)/3465.)*v))));
	        break;
	      case 3:
	        deltalm = vh3*((486 + 4961*eta)/(810.*(1 - 2*eta)) + vh*((11*a)/4. 
			+ (1571*LAL_PI*vh*vh)/385.));
		rholm	= 1 + v*((5*(chiA - chiS*dM)*eta)/(8.*dM*(-1 + 2*eta)) 
			+ v*((222 - 547*eta + 160*eta2)/(176.*(-1 + 2*eta)) 
			+ v2*(-6894273./7047040. + (3*a2)/8. + v*((-12113*a)/6160. 
			+ (1664224207351./195343948800. - (1571*eulerlogxabs)/770.)*v))));
	        break;
	      case 2:
		deltalm = vh3*((7*(1 + 6*eta))/(15.*(1 - 3*eta)) + ((212*a)/75. 
			+ (6284*LAL_PI)/3465.)*vh3);
		rholm	= 1 + v2*((1146 - 3530*eta + 285*eta2)/(1320.*(-1 + 3*eta)) 
			+ v*((chiA*(10 - 21*eta)*dM + chiS*(10 - 59*eta 
			+ 78*eta2))/(-15 + 45*eta) + v*(a2/2. + (-114859044 
			+ 295834536*eta + 1204388696*eta2 - 3047981160*eta3 
			- 379526805*eta2*eta2)/(317116800.*(1.-3.*eta)*(1.-3.*eta)) + v*((-7*a)/110. 
			+ (848238724511./219761942400. + (2323*a2)/3960. 
			- (3142*eulerlogxabs)/3465.)*v))));
		break;
	      case 1:
		deltalm = vh3*((2 + 507*eta)/(10.*(1 - 2*eta)) + vh*((11*a)/12. 
			+ (1571*LAL_PI*vh*vh)/3465.));
		rholm	= 1 + v*((5*(chiA - chiS*dM)*eta)/(8.*dM*(-1 + 2*eta)) 
			+ v*((602 - 1385*eta + 288*eta2)/(528.*(-1 + 2*eta)) 
			+ v2*(-7775491./21141120. + (3*a2)/8. + v*((-20033*a)/55440. 
			- (5*a*a2)/6. + (1227423222031./1758095539200. 
			- (1571*eulerlogxabs)/6930.)*v))));
		break;
	      default:
                XLAL_ERROR( func, XLAL_EINVAL );
                break;
	    }
	    break;
	  case 5:
	    switch (m)
	    {
	      case 5:
	        deltalm = ((96875. + 857528.*eta)*vh3)/(131250.*(1 - 2*eta));
		rholm	= 1 + v2*((487 - 1298*eta + 512*eta2)/(390.*(-1 + 2*eta)) 
			+ v*((-2*a)/3. + v*(-3353747./2129400. + a2/2. - (241*a*v)/195.)));
	        break;
	      case 4:
		deltalm = vh3*(8./15. + (12*a*vh)/5.);
		rholm	= 1 + v2*((-17448 + 96019*eta - 127610*eta2 
			+ 33320*eta3)/(13650.*(1 - 5*eta + 5*eta2)) + v*((-2*a)/15. 
			+ (-16213384./15526875. + (2*a2)/5.)*v));
		break;
	      case 3:
	        deltalm = (31*vh3)/70.;
		rholm	= 1 + v2*((375 - 850*eta + 176*eta2)/(390.*(-1 + 2*eta)) 
			+ v*((-2*a)/3. + v*(-410833./709800. + a2/2. - (103*a*v)/325.)));
	        break;
	      case 2:
		deltalm = vh3*(4./15. + (6*a*vh)/5.);
		rholm	= 1 + v2*((-15828 + 84679*eta - 104930*eta2 
			+ 21980*eta3)/(13650.*(1 - 5*eta + 5*eta2)) + v*((-2*a)/15. 
			+ (-7187914./15526875. + (2*a2)/5.)*v));
		break;
	      case 1:
		deltalm = (31*vh3)/210.;
		rholm	= 1 + v2*((319 - 626*eta + 8*eta2)/(390.*(-1 + 2*eta)) 
			+ v*((-2*a)/3. + v*(-31877./304200. + a2/2. + (139*a*v)/975.)));
		break;
	      default:
                XLAL_ERROR( func, XLAL_EINVAL );
                break;
	    }
	    break;
	  case 6:
	    switch (m)
	    {
	      case 6:
	        deltalm = (43*vh3)/70.;
		rholm	= 1 + v2*((-106 + 602*eta - 861*eta2 
			+ 273*eta3)/(84.*(1 - 5*eta + 5*eta2)) + v*((-2*a)/3. 
			+ (-1025435./659736. + a2/2.)*v));
	        break;
	      case 5:
		deltalm = (10*vh3)/21.;
		rholm	= 1 + v2*((-185 + 838*eta - 910*eta2 
			+ 220*eta3)/(144.*(1 - 4*eta + 3*eta2)) - (2*a*v)/9.);
		break;
	      case 4:
		deltalm = (43*vh3)/105.;
		rholm	= 1 + v2*((-86 + 462*eta - 581*eta2 
			+ 133*eta3)/(84.*(1 - 5*eta + 5*eta2)) + v*((-2*a)/3. 
			+ (-476887./659736. + a2/2.)*v));
		break;
	      case 3:
	        deltalm = (2*vh3)/7.;
		rholm	= 1 + v2*((-169 + 742*eta - 750*eta2 
			+ 156*eta3)/(144.*(1 - 4*eta + 3*eta2)) - (2*a*v)/9.);
	        break;
	      case 2:
		deltalm = (43*vh3)/210.;
		rholm	= 1 + v2*((-74 + 378*eta - 413*eta2 
			+ 49*eta3)/(84.*(1 - 5*eta + 5*eta2)) + v*((-2*a)/3. 
			+ (-817991./3298680. + a2/2.)*v));
		break;
	      case 1:
		deltalm = (2*vh3)/21.;
		rholm	= 1 + v2*((-161 + 694*eta - 670*eta2 
			+ 124*eta3)/(144.*(1 - 4*eta + 3*eta2)) - (2*a*v)/9.);
		break;
	      default:
                XLAL_ERROR( func, XLAL_EINVAL );
                break;
	    }
	    break;
	  case 7:
	    switch (m)
	    {
	      case 7:
	        deltalm = (19*vh3)/36.;
		rholm	= 1 + v2*((-906 + 4246*eta - 4963*eta2 
			+ 1380*eta3)/(714.*(1 - 4*eta + 3*eta2)) - (2*a*v)/3.);
	        break;
	      case 6:
	        deltalm = 0.0;
		rholm	= 1 + ((2144 - 16185*eta + 37828*eta2 - 29351*eta3 
			+ 6104*eta2*eta2)*v2)/(1666.*(-1 + 7*eta - 14*eta2 
			+ 7*eta3));
	        break;
	      case 5:
		deltalm = (95*vh3)/252.;
		rholm	= 1 + v2*((-762 + 3382*eta - 3523*eta2 
			+ 804*eta3)/(714.*(1 - 4*eta + 3*eta2)) - (2*a*v)/3.);
		break;
	      case 4:
		deltalm = 0.0;
		rholm	= 1 + ((17756 - 131805*eta + 298872*eta2 - 217959*eta3 
			+ 41076*eta2*eta2)*v2)/(14994.*(-1 + 7*eta - 14*eta2 
			+ 7*eta3));
		break;
	      case 3:
	        deltalm = (19*vh3)/84.;
		rholm	= 1 + v2*((-666 + 2806*eta - 2563*eta2 
			+ 420*eta3)/(714.*(1 - 4*eta + 3*eta2)) - (2*a*v)/3.);
	        break;
	      case 2:
		deltalm = 0.0;
		rholm	= 1 + ((16832 - 123489*eta + 273924*eta2 - 190239*eta3 
			+ 32760*eta2*eta2)*v2)/(14994.*(-1 + 7*eta - 14*eta2 
			+ 7*eta3));
		break;
	      case 1:
		deltalm = (19*vh3)/252.;
		rholm	= 1 + v2*((-618 + 2518*eta - 2083*eta2 
			+ 228*eta3)/(714.*(1 - 4*eta + 3*eta2)) - (2*a*v)/3.);
		break;
	      default:
                XLAL_ERROR( func, XLAL_EINVAL );
                break;
	    }
	    break;
	  case 8:
	    switch (m)
	    {
	      case 8:
	        deltalm = 0.0;
		rholm	= 1 + ((3482 - 26778*eta + 64659*eta2 - 53445*eta3 
			+ 12243*eta2*eta2)*v2)/(2736.*(-1 + 7*eta - 14*eta2 
			+ 7*eta3));
	        break;
	      case 7:
		deltalm = 0.0;
		rholm	= 1 + ((23478 - 154099*eta + 309498*eta2 - 207550*eta3 
			+ 38920*eta2*eta2)*v2)/(18240.*(-1 + 6*eta - 10*eta2 
			+ 4*eta3));
		break;
	      case 6:
	        deltalm = 0.0;
		rholm	= 1 + ((1002 - 7498*eta + 17269*eta2 - 13055*eta3 
			+ 2653*eta2*eta2)*v2)/(912.*(-1 + 7*eta - 14*eta2 
			+ 7*eta3));
	        break;
	      case 5:
		deltalm = 0.0;
		rholm	= 1 + ((4350 - 28055*eta + 54642*eta2 - 34598*eta3 
			+ 6056*eta2*eta2)*v2)/(3648.*(-1 + 6*eta - 10*eta2 
			+ 4*eta3));
		break;
	      case 4:
		deltalm = 0.0;
		rholm	= 1 + ((2666 - 19434*eta + 42627*eta2 - 28965*eta3 
			+ 4899*eta2*eta2)*v2)/(2736.*(-1 + 7*eta - 14*eta2 
			+ 7*eta3));
		break;
	      case 3:
	        deltalm = 0.0;
		rholm	= 1 + ((20598 - 131059*eta + 249018*eta2 - 149950*eta3 
			+ 24520*eta2*eta2)*v2)/(18240.*(-1 + 6*eta - 10*eta2 
			+ 4*eta3));
	        break;
	      case 2:
		deltalm = 0.0;
		rholm	= 1 + ((2462 - 17598*eta + 37119*eta2 - 22845*eta3 
			+ 3063*eta2*eta2)*v2)/(2736.*(-1 + 7*eta - 14*eta2 
			+ 7*eta3));
		break;
	      case 1:
		deltalm = 0.0;
		rholm	= 1 + ((20022 - 126451*eta + 236922*eta2 - 138430*eta3 
			+ 21640*eta2*eta2)*v2)/(18240.*(-1 + 6*eta - 10*eta2 
			+ 4*eta3));
		break;
	      default:
                XLAL_ERROR( func, XLAL_EINVAL );
                break;
	    }
	    break;
	  default:
            XLAL_ERROR( func, XLAL_EINVAL );
            break; 
	}

	*hlm = XLALCOMPLEX16MulReal( XLALCOMPLEX16Mul( Tlm, XLALCOMPLEX16Polar( 1.0, deltalm) ), 
				     Slm*pow(rholm,l) );
        *hlm = XLALCOMPLEX16Mul( *hlm, hNewton );
	return XLAL_SUCCESS;
}

static inline
REAL8 XLALCalculateA5( REAL8 eta )
{
  return -12.9499 + 204.779 * eta - 206.319 *eta*eta;
}

static
REAL8 XLALCalculateEOBA( REAL8   r,
                         REAL8 eta)
{
	REAL8 eta2, u, u2, u3, u4, a4, a5, NA, DA;

	eta2 = eta*eta;
	u = 1./r;
	u2 = u*u;
	u3 = u2*u;
	u4 = u2*u2;
	a4 = ninty4by3etc * eta;
	a5 = XLALCalculateA5( eta );
	NA = (32. - 24.*eta - 4.*a4 - a5*eta)*u + (a4 - 16. + 8.*eta);
	DA = a4 - 16. + 8.*eta - (2.*a4 + a5*eta + 8.*eta)*u - (4.*a4 + 2.*a5*eta + 16.*eta)*u2
	   - (8.*a4 + 4.*a5*eta + 2.*a4*eta + 16.*eta2)*u3
	   + (-a4*a4 - 8.*a5*eta - 8.*a4*eta + 2.*a5*eta2 - 16.*eta2)*u4;
	
	return NA/DA;	
}

static
REAL8 XLALCalculateEOBD( REAL8   r,
                         REAL8 eta)
{
	REAL8  u, u2, u3;

	u = 1./r;
	u2 = u*u;
	u3 = u2*u;
	
	return 1./(1.+6.*eta*u2+2.*eta*(26.-3.*eta)*u3);	
}

static
REAL8 XLALEffectiveHamiltonian( REAL8Vector           *values,
                                InspiralDerivativesIn *ak)
{
	REAL8 eta, r, phi, pr, pp, r2, pr2, pp2, z3, eoba, eobd;

   	eta = ak->coeffs->eta;

   	r   = values->data[0];
   	phi = values->data[1];
   	pr  = values->data[2];
   	pp  = values->data[3];

	r2   = r * r;
	pr2  = pr * pr;
	pp2  = pp * pp;
	
	eoba = XLALCalculateEOBA( r, eta );
	eobd = XLALCalculateEOBD( r, eta );
	z3   = 2. * ( 4. - 3. * eta ) * eta;
	return sqrt( eoba * ( 1. + eoba/eobd*pr2 + pp2/r2 + z3*pr2*pr2/r2 ) );	
}

/*-------------------------------------------------------------------*/
/*                      pseudo-4PN functions                         */
/*-------------------------------------------------------------------*/

void
LALpphiInitP4PN(
            REAL8 *phase,
            REAL8 r,
            REAL8 eta,
            REAL8 omegaS
            )
{
  REAL8 u, u2, u3, u4, a4, a5, eta2, NA, DA, A, dA;

  eta2 = eta*eta;
  u = 1./r;
  u2 = u*u;
  u3 = u2*u;
  u4 = u2*u2;
  a4 = (ninty4by3etc - 2. * omegaS) * eta;
  a5 = XLALCalculateA5( eta );
  NA = (32. - 24.*eta - 4.*a4 - a5*eta)*u + (a4 - 16. + 8.*eta);
  DA = a4 - 16. + 8.*eta - (2.*a4 + a5*eta + 8.*eta)*u - (4.*a4 + 2.*a5*eta + 16.*eta)*u2
       - (8.*a4 + 4.*a5*eta + 2.*a4*eta + 16.*eta2)*u3
       + (-a4*a4 - 8.*a5*eta - 8.*a4*eta + 2.*a5*eta2 - 16.*eta2)*u4;
  A = NA/DA;
  dA = ( (32. - 24.*eta - 4.*a4 - a5*eta) * DA - NA *
         ( -(2.*a4 + a5*eta + 8.*eta) - 2.*(4.*a4 + 2.*a5*eta + 16.*eta)*u
          - 3.*(8.*a4 + 4.*a5*eta + 2.*a4*eta + 16.*eta2)*u2
          + 4.*(-a4*a4 - 8.*a5*eta - 8.*a4*eta + 2.*a5*eta2 - 16.*eta2)*u3))/(DA*DA);
  *phase = sqrt(-dA/(2.*u*A + u2 * dA));
/* why is it called phase? This is initial j!? */
}

/*-------------------------------------------------------------------*/
  void
LALprInitP4PN(
             LALStatus *status,
             REAL8 *pr,
             REAL8 p,
             void *params
             )
{
  REAL8   u, u2, u3, u4, p2, p3, p4, q2, A, DA, NA;
  REAL8  onebyD, AbyD, Heff, HReal, etahH;
  REAL8 eta, eta2, a4, a5, z3, r, vr, q;
  pr3In *ak;

  ak = (pr3In *) params;

  eta = ak->eta;
  vr = ak->vr;
  r = ak->r;
  q = ak->q;
  eta2 = eta*eta;


   p2 = p*p;
   p3 = p2*p;
   p4 = p2*p2;
   q2 = q*q;
   u = 1./ r;
   u2 = u*u;
   u3 = u2 * u;
   u4 = u2 * u2;
   z3 = 2. * (4. - 3. * eta) * eta;
   a4 = ninty4by3etc * eta;
   a5 = XLALCalculateA5( eta );

   NA = (32. - 24.*eta - 4.*a4 - a5*eta)*u + (a4 - 16. + 8.*eta);
   DA = a4 - 16. + 8.*eta - (2.*a4 + a5*eta + 8.*eta)*u - (4.*a4 + 2.*a5*eta + 16.*eta)*u2
        - (8.*a4 + 4.*a5*eta + 2.*a4*eta + 16.*eta2)*u3
        + (-a4*a4 - 8.*a5*eta - 8.*a4*eta + 2.*a5*eta2 - 16.*eta2)*u4;
   A = NA/DA;
   onebyD = 1. + 6.*eta*u2 + 2. * ( 26. - 3. * eta) * eta * u3 + 36.*eta2*u4;
   AbyD = A * onebyD;

   Heff = pow (A*(1. + AbyD * p2 + q*q * u2 + z3 * p4 * u2), 0.5);
   HReal = pow (1. + 2.*eta*(Heff - 1.), 0.5) / eta;
   etahH = eta*Heff*HReal;

   *pr = -vr +  A*(AbyD*p + 2. * z3 * u2 * p3)/etahH;
/* This sets pr = dH/dpr - vr, calls rootfinder,
   gets value of pr s.t. dH/pr = vr */
   RETURN(status);
}


/*-------------------------------------------------------------------*/
static void
omegaofrP4PN (
             REAL8 *x,
             REAL8 r,
             void *params)
{
   REAL8 u, u2, u3, u4, a4, a5, eta, eta2, NA, DA, A, dA;

   /*include a status here ?*/
   pr3In *ak;
   ak = (pr3In *) params;
   eta = ak->eta;
   eta2 = eta*eta;

   u = 1./r;
   u2 = u*u;
   u3 = u2*u;
   u4 = u2*u2;
   a4 = ninty4by3etc * eta;
   a5 = XLALCalculateA5( eta );
   NA = (32. - 24.*eta - 4.*a4 - a5*eta)*u + (a4 - 16. + 8.*eta);
   DA = a4 - 16. + 8.*eta - (2.*a4 + a5*eta + 8.*eta)*u - (4.*a4 + 2.*a5*eta + 16.*eta)*u2
        - (8.*a4 + 4.*a5*eta + 2.*a4*eta + 16.*eta2)*u3
        + (-a4*a4 - 8.*a5*eta - 8.*a4*eta + 2.*a5*eta2 - 16.*eta2)*u4;
   A = NA/DA;
   dA = ( (32. - 24.*eta - 4.*a4 - a5*eta) * DA - NA *
          ( -(2.*a4 + a5*eta + 8.*eta) - 2.*(4.*a4 + 2.*a5*eta + 16.*eta)*u
           - 3.*(8.*a4 + 4.*a5*eta + 2.*a4*eta + 16.*eta2)*u2
           + 4.*(-a4*a4 - 8.*a5*eta - 8.*a4*eta + 2.*a5*eta2 - 16.*eta2)*u3))/(DA*DA);

   *x = pow(u,1.5) * sqrt ( -0.5 * dA /(1. + 2.*eta * (A/sqrt(A+0.5 * u*dA)-1.)));

}


/*-------------------------------------------------------------------*/

REAL8 
XLALrOfOmegaP4PN(
            REAL8 r,
            void *params )
{
  REAL8  omega1,omega2;
  pr3In *pr3in;

#ifndef LAL_NDEBUG
  if ( !params )
    XLAL_ERROR_REAL8( "XLALrOfOmegaP4PN", XLAL_EFAULT );
#endif

  pr3in = (pr3In *) params;

  omega1 = pr3in->omega;
  omegaofrP4PN(&omega2,r, params);
  return ( -omega1 + omega2 );
}

/*-------------------------------------------------------------------*/
static void
LALlightRingRadiusP4PN(
                      LALStatus *status,
                      REAL8 *x,
                      REAL8 r,
                      void *params
                      )
{
  REAL8 eta, eta2, u, u2, u3, u4, a4, a5, NA, DA, A, dA;
  rOfOmegaIn *rofomegain;
  status = NULL;
  rofomegain = (rOfOmegaIn *) params;
  eta = rofomegain->eta;
  eta2 = eta*eta;


  u = 1./r;
  u2 = u*u;
  u3 = u2*u;
  u4 = u2*u2;
  a4 = ninty4by3etc * eta;
  a5 = XLALCalculateA5( eta );
  NA = (32. - 24.*eta - 4.*a4 - a5*eta)*u + (a4 - 16. + 8.*eta);
  DA = a4 - 16. + 8.*eta - (2.*a4 + a5*eta + 8.*eta)*u - (4.*a4 + 2.*a5*eta + 16.*eta)*u2
       - (8.*a4 + 4.*a5*eta + 2.*a4*eta + 16.*eta2)*u3
       + (-a4*a4 - 8.*a5*eta - 8.*a4*eta + 2.*a5*eta2 - 16.*eta2)*u4;
  A = NA/DA;
  dA = ( (32. - 24.*eta - 4.*a4 - a5*eta) * DA - NA *
         ( -(2.*a4 + a5*eta + 8.*eta) - 2.*(4.*a4 + 2.*a5*eta + 16.*eta)*u
          - 3.*(8.*a4 + 4.*a5*eta + 2.*a4*eta + 16.*eta2)*u2
          + 4.*(-a4*a4 - 8.*a5*eta - 8.*a4*eta + 2.*a5*eta2 - 16.*eta2)*u3))/(DA*DA);
  *x = 2 * A + dA * u;
}

/*-------------------------------------------------------------------*/
 void
LALHCapDerivativesP4PN(
                  REAL8Vector *values,
                  REAL8Vector *dvalues,
                  void        *funcParams
                  )
{

   InspiralDerivativesIn *ak;

   REAL8 r, s, p, q, u, u2, u3, u4, u5, p2, p3, p4, q2, Apot, DA, NA;
   REAL8  dA, onebyD, DonebyD, AbyD, Heff, HReal, etahH;
   REAL8 omega, v, eta, eta2, a4, z2, z30, z3, zeta2;
   REAL8 a5, c1;
   double dr, ds, dp, dq;

   ak = (InspiralDerivativesIn *) funcParams;

   eta = ak->coeffs->eta;
   zeta2 = ak->coeffs->zeta2;

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
   u5 = u*u4;
   z30 = 2.L * (4.L - 3.L * eta) * eta;
   z2 = 0.75L * z30 * zeta2,
   z3 = z30 * (1.L - zeta2);
   eta2 = eta*eta;

   a4 = ninty4by3etc * eta;
   a5 = XLALCalculateA5( eta );

   NA = (32. - 24.*eta - 4.*a4 - a5*eta)*u + (a4 - 16. + 8.*eta);
   DA = a4 - 16. + 8.*eta - (2.*a4 + a5*eta + 8.*eta)*u - (4.*a4 + 2.*a5*eta + 16.*eta)*u2
        - (8.*a4 + 4.*a5*eta + 2.*a4*eta + 16.*eta2)*u3
        + (-a4*a4 - 8.*a5*eta - 8.*a4*eta + 2.*a5*eta2 - 16.*eta2)*u4;
   Apot = NA/DA; /* This A(u) assume zeta2=0 (the default value) */

   onebyD = 1. + 6.*eta*u2 + (2.*eta * ( 26. - 3.*eta ) - z2)* u3 + 36.*eta2*u4;
   AbyD = Apot * onebyD;
   Heff = pow (Apot*(1. + AbyD * p2 + q*q * u2 + z3 * p4 * u2), 0.5);
   HReal = pow (1. + 2.*eta*(Heff - 1.), 0.5) / eta;
   dA = -u2 * ( (32. - 24.*eta - 4.*a4 - a5*eta) * DA - NA *
          ( -(2.*a4 + a5*eta + 8.*eta) - 2.*(4.*a4 + 2.*a5*eta + 16.*eta)*u
          - 3.*(8.*a4 + 4.*a5*eta + 2.*a4*eta + 16.*eta2)*u2
          + 4.*(-a4*a4 - 8.*a5*eta - 8.*a4*eta + 2.*a5*eta2 - 16.*eta2)*u3))/(DA*DA);

   DonebyD = -12.*eta*u3 - (6.*eta*(26. - 3.*eta) - z2)*u4 - 144.*eta2*u5;
   etahH = eta*Heff*HReal;

   dr = dvalues->data[0] = Apot*(AbyD*p +  z30 * u2 *(2.* p3
              + zeta2*(-0.5*p3 + 0.75*p*q2*u2)))/etahH;
   ds = dvalues->data[1] = omega = Apot * q * u2 * (1. + 0.75*z30*zeta2*p2*u2)/ etahH;
   v = pow(omega,oneby3);

   dp = dvalues->data[2] = -0.5 * Apot * (dA*Heff*Heff/(Apot*Apot) - 2.*q2*u3
              + (dA * onebyD + Apot * DonebyD) * p2
              + z30 * u3 *(-2.* p4+zeta2*(0.5*p4 - 3.0*p2*q2*u2))) / etahH;
   c1 = 1.+(u2 - 2.*u3*Apot/dA) * q2;/*below:dpphi/dt = F_RR*/
   dq = dvalues->data[3] = - ak->flux(v,ak->coeffs)/(eta * omega);
   /*
   fprintf(stdout, "%e %e %e %e %e %e %e %e %e %e %e\n", r, s, p, q, Heff, v, Apot, dr, ds, dp, dq);
   */
}

/* Version which uses the factorized flux */
void
LALHCapDerivativesP4PNFF(
					   REAL8Vector *values,
					   REAL8Vector *dvalues,
					   void        *funcParams
					   )
{

  InspiralDerivativesIn *ak;

  /* Max l to sum up to in the factorized flux */
  const INT4 lMax = 8;

  REAL8 r, s, p, q, u, u2, u3, u4, u5, p2, p3, p4, q2, Apot, DA, NA;
  REAL8  dA, onebyD, DonebyD, AbyD, Heff, HReal, etahH;
  REAL8 omega, v, eta, eta2, a4, z2, z30, z3, zeta2;
  REAL8 a5, c1;
  double dr, ds, dp, dq;

  ak = (InspiralDerivativesIn *) funcParams;

  eta = ak->coeffs->eta;
  zeta2 = ak->coeffs->zeta2;

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
  u5 = u*u4;
  z30 = 2.L * (4.L - 3.L * eta) * eta;
  z2 = 0.75L * z30 * zeta2,
  z3 = z30 * (1.L - zeta2);
  eta2 = eta*eta;

  a4 = ninty4by3etc * eta;
  a5 = XLALCalculateA5( eta );

  NA = (32. - 24.*eta - 4.*a4 - a5*eta)*u + (a4 - 16. + 8.*eta);
  DA = a4 - 16. + 8.*eta - (2.*a4 + a5*eta + 8.*eta)*u - (4.*a4 + 2.*a5*eta + 16.*eta)*u2
       - (8.*a4 + 4.*a5*eta + 2.*a4*eta + 16.*eta2)*u3
       + (-a4*a4 - 8.*a5*eta - 8.*a4*eta + 2.*a5*eta2 - 16.*eta2)*u4;
  Apot = NA/DA; /* This A(u) assume zeta2=0 (the default value) */

  onebyD = 1. + 6.*eta*u2 + (2.*eta * ( 26. - 3.*eta ) - z2)* u3 + 36.*eta2*u4;
  AbyD = Apot * onebyD;
  Heff = pow (Apot*(1. + AbyD * p2 + q*q * u2 + z3 * p4 * u2), 0.5);
  HReal = pow (1. + 2.*eta*(Heff - 1.), 0.5) / eta;
  dA = -u2 * ( (32. - 24.*eta - 4.*a4 - a5*eta) * DA - NA *
       ( -(2.*a4 + a5*eta + 8.*eta) - 2.*(4.*a4 + 2.*a5*eta + 16.*eta)*u
        - 3.*(8.*a4 + 4.*a5*eta + 2.*a4*eta + 16.*eta2)*u2
        + 4.*(-a4*a4 - 8.*a5*eta - 8.*a4*eta + 2.*a5*eta2 - 16.*eta2)*u3))/(DA*DA);
	
  DonebyD = -12.*eta*u3 - (6.*eta*(26. - 3.*eta) - z2)*u4 - 144.*eta2*u5;
  etahH = eta*Heff*HReal;

  dr = dvalues->data[0] = Apot*(AbyD*p +  z30 * u2 *(2.* p3
       + zeta2*(-0.5*p3 + 0.75*p*q2*u2)))/etahH;
  ds = dvalues->data[1] = omega = Apot * q * u2 * (1. + 0.75*z30*zeta2*p2*u2)/ etahH;
  v = pow(omega,oneby3);

  dp = dvalues->data[2] = -0.5 * Apot * (dA*Heff*Heff/(Apot*Apot) - 2.*q2*u3
       + (dA * onebyD + Apot * DonebyD) * p2
       + z30 * u3 *(-2.* p4+zeta2*(0.5*p4 - 3.0*p2*q2*u2))) / etahH;
  c1 = 1.+(u2 - 2.*u3*Apot/dA) * q2;/*below:dpphi/dt = F_RR*/

  dq = XLALInspiralFactorizedFlux( values, dvalues, ak, lMax );
  /* This function can fail */
  /* TODO: Implement proper error checking */
  dq = dvalues->data[3] = - dq /(eta * omega);
  
   fprintf(stdout, "Hcapd: %e %e %e %e %e %e %e %e %e %e %e\n", r, s, p, q, Heff, v, Apot, dr, ds, dp, dq);
}


/*-------------------------------------------------------------------*/
 void LALvrP4PN(REAL8 *vr, void *params )
{
  REAL8 A, dA, d2A, NA, DA, dDA, dNA, d2DA;
  REAL8 u, u2, u3, u4, v, x1;
  REAL8 eta, eta2, a4, a5, FDIS;

  pr3In *pr3in;
  pr3in = (pr3In *)params;

  eta = pr3in->eta;
  u = 1./ pr3in->r;

  u2 = u*u;
  u3 = u2*u;
  u4 = u2*u2;
  eta2 = eta*eta;


  a4 = (ninty4by3etc - 2. * pr3in->omegaS) * eta;
  a5 = XLALCalculateA5( eta );
  NA = (32. - 24.*eta - 4.*a4 - a5*eta)*u + (a4 - 16. + 8.*eta);
  DA = a4 - 16. + 8.*eta - (2.*a4 + a5*eta + 8.*eta)*u - (4.*a4 + 2.*a5*eta + 16.*eta)*u2
       - (8.*a4 + 4.*a5*eta + 2.*a4*eta + 16.*eta2)*u3
       + (-a4*a4 - 8.*a5*eta - 8.*a4*eta + 2.*a5*eta2 - 16.*eta2)*u4;
  A = NA/DA;
  dNA = (32. - 24.*eta - 4.*a4 - a5*eta);
  dDA = - (2.*a4 + a5*eta + 8.*eta) - 2.*(4.*a4 + 2.*a5*eta + 16.*eta)*u
       - 3.*(8.*a4 + 4.*a5*eta + 2.*a4*eta + 16.*eta2)*u2
       + 4.*(-a4*a4 - 8.*a5*eta - 8.*a4*eta + 2.*a5*eta2 - 16.*eta2)*u3;
  d2DA = - 2.*(4.*a4 + 2.*a5*eta + 16.*eta)
       - 6.*(8.*a4 + 4.*a5*eta + 2.*a4*eta + 16.*eta2)*u
       + 12.*(-a4*a4 - 8.*a5*eta - 8.*a4*eta + 2.*a5*eta2 - 16.*eta2)*u2;

  dA = (dNA * DA - NA * dDA)/ (DA*DA);
  d2A = (-NA * DA * d2DA - 2. * dNA * DA * dDA + 2. * NA * dDA * dDA)/pow(DA,3.);
  v = pow(pr3in->omega,oneby3);
  FDIS = -pr3in->in3copy.flux(v, pr3in->in3copy.coeffs)/(eta* pr3in->omega);
  x1 = -1./u2 * sqrt (-dA * pow(2.* u * A + u2 * dA, 3.) )
                / (2.* u * dA * dA + A*dA - u * A * d2A);
  *vr = FDIS * x1;
}


/*-------------------------------------------------------------------*/

/*  <lalVerbatim file="LALEOBPPWaveformCP"> */
void
LALEOBPPWaveform (
   LALStatus        *status,
   REAL4Vector      *signalvec,
   InspiralTemplate *params
   )
{ /* </lalVerbatim> */

   UINT4 count;
   InspiralInit paramsInit;
   INITSTATUS(status, "LALEOBPPWaveform", LALEOBPPWAVEFORMC);
   ATTATCHSTATUSPTR(status);

   ASSERT(signalvec,  status,
	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signalvec->data,  status,
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

   LALInspiralSetup (status->statusPtr, &(paramsInit.ak), params);
   CHECKSTATUSPTR(status);
   LALInspiralChooseModel(status->statusPtr, &(paramsInit.func),
					 &(paramsInit.ak), params);
   CHECKSTATUSPTR(status);

   memset(signalvec->data, 0, signalvec->length * sizeof( REAL4 ));

   /* Call the engine function */
   LALEOBPPWaveformEngine(status->statusPtr, signalvec, NULL, NULL, NULL,
			NULL, NULL, &count, params, &paramsInit);
   CHECKSTATUSPTR( status );

   DETATCHSTATUSPTR(status);
   RETURN(status);
}


NRCSID (LALEOBPPWAVEFORMTEMPLATESC,
"$Id$");

/*  <lalVerbatim file="LALEOBPPWaveformTemplatesCP"> */

void
LALEOBPPWaveformTemplates (
   LALStatus        *status,
   REAL4Vector      *signalvec1,
   REAL4Vector      *signalvec2,
   InspiralTemplate *params
   )
{ /* </lalVerbatim> */

   UINT4 count;

   InspiralInit paramsInit;

   INITSTATUS(status, "LALEOBPPWaveformTemplates", LALEOBPPWAVEFORMTEMPLATESC);
   ATTATCHSTATUSPTR(status);

   ASSERT(signalvec1,  status,
	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signalvec2,  status,
   	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signalvec1->data,  status,
   	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signalvec2->data,  status,
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

   LALInspiralSetup (status->statusPtr, &(paramsInit.ak), params);
   CHECKSTATUSPTR(status);
   LALInspiralChooseModel(status->statusPtr, &(paramsInit.func),
					&(paramsInit.ak), params);
   CHECKSTATUSPTR(status);

   memset(signalvec1->data, 0, signalvec1->length * sizeof( REAL4 ));
   memset(signalvec2->data, 0, signalvec2->length * sizeof( REAL4 ));

   /* Call the engine function */
   LALEOBPPWaveformEngine(status->statusPtr, signalvec1, signalvec2, NULL, NULL,
			   NULL, NULL, &count, params, &paramsInit);
   CHECKSTATUSPTR( status );

   DETATCHSTATUSPTR(status);
   RETURN(status);
}


/*=========================================================*/
/*======INJECTION =========================================*/
/*=========================================================*/

/*  <lalVerbatim file="LALEOBPPWaveformForInjectionCP"> */
void
LALEOBPPWaveformForInjection (
			    LALStatus        *status,
			    CoherentGW       *waveform,
			    InspiralTemplate *params,
			    PPNParamStruc    *ppnParams
			    )
{
  /* </lalVerbatim> */
  UINT4 count, i;

  REAL4Vector *a=NULL;/* pointers to generated amplitude  data */
  REAL4Vector *h=NULL;/* pointers to generated polarization data */
  REAL4Vector *ff=NULL ;/* pointers to generated  frequency data */
  REAL8Vector *phi=NULL;/* pointer to generated phase data */

  REAL8 s;

  REAL8 phiC;/* phase at coalescence */
  CHAR message[256];
  InspiralInit paramsInit;

  CreateVectorSequenceIn in;

  INITSTATUS(status, "LALEOBPPWaveformForInjection", LALEOBPPWAVEFORMTEMPLATESC);
  ATTATCHSTATUSPTR(status);


  /* Make sure parameter and waveform structures exist. */
  ASSERT( params, status, LALINSPIRALH_ENULL,  LALINSPIRALH_MSGENULL );
  ASSERT(waveform, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  /* Make sure waveform fields don't exist. */
  ASSERT( !( waveform->a ), status,
  	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->h ), status,
  	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->f ), status,
  	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->phi ), status,
  	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->shift ), status,
  	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );

  params->ampOrder = 0;
  sprintf(message, "WARNING: Amp Order has been reset to %d", params->ampOrder);
  LALInfo(status, message);

  /* Compute some parameters*/
  LALInspiralInit(status->statusPtr, params, &paramsInit);
  CHECKSTATUSPTR(status);

  if (paramsInit.nbins==0)
    {
      DETATCHSTATUSPTR(status);
      RETURN (status);
    }
  /* Now we can allocate memory and vector for coherentGW structure*/
  LALSCreateVector(status->statusPtr, &ff, paramsInit.nbins);
  CHECKSTATUSPTR(status);
  LALSCreateVector(status->statusPtr, &a, 2*paramsInit.nbins);
  CHECKSTATUSPTR(status);
  LALDCreateVector(status->statusPtr, &phi, paramsInit.nbins);
  CHECKSTATUSPTR(status);
  LALSCreateVector(status->statusPtr, &h, 2*paramsInit.nbins);
  CHECKSTATUSPTR(status);

  /* By default the waveform is empty */
  memset(ff->data, 0, paramsInit.nbins * sizeof(REAL4));
  memset(a->data, 0, 2 * paramsInit.nbins * sizeof(REAL4));
  memset(phi->data, 0, paramsInit.nbins * sizeof(REAL8));
  memset(h->data, 0, 2 * paramsInit.nbins * sizeof(REAL4));

  /* Call the engine function */
  params->startPhase = ppnParams->phi;
  LALEOBPPWaveformEngine(status->statusPtr, NULL, NULL, h, a, ff,
			   phi, &count, params, &paramsInit);
  BEGINFAIL( status )
  {
     LALSDestroyVector(status->statusPtr, &ff);
     CHECKSTATUSPTR(status);
     LALSDestroyVector(status->statusPtr, &a);
     CHECKSTATUSPTR(status);
     LALDDestroyVector(status->statusPtr, &phi);
     CHECKSTATUSPTR(status);
     if( params->approximant == EOBNR_PP )
     {
       LALSDestroyVector(status->statusPtr, &h);
       CHECKSTATUSPTR(status);
     }
  }
  ENDFAIL( status );

  /* Check an empty waveform hasn't been returned */
  for (i = 0; i < phi->length; i++)
  {
    if (phi->data[i] != 0.0) break;
    if (i == phi->length - 1)
    {
      LALSDestroyVector(status->statusPtr, &ff);
      CHECKSTATUSPTR(status);
      LALSDestroyVector(status->statusPtr, &a);
      CHECKSTATUSPTR(status);
      LALDDestroyVector(status->statusPtr, &phi);
      CHECKSTATUSPTR(status);
      LALSDestroyVector(status->statusPtr, &h);
      CHECKSTATUSPTR(status);

      DETATCHSTATUSPTR( status );
      RETURN( status );
    }
  }

  s = 0.5 * phi->data[count - 1];

  sprintf(message, "fFinal = %f", params->fFinal);
  LALInfo(status, message);

  sprintf(message, "cycles = %f", s/3.14159);
  LALInfo(status, message);

  sprintf( message, "final coalescence phase with respet to actual data =%f ",
  	(ff->data[count]-ff->data[count-1])/2/3.14159);
  LALInfo(status, message);



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
	  phi->data[i] =  -phiC + phi->data[i] + ppnParams->phi;
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
      if ( ( waveform->h = (REAL4TimeVectorSeries *)
             LALMalloc( sizeof(REAL4TimeVectorSeries) ) ) == NULL )
      {
        ABORT( status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM );
      }
      memset( waveform->h, 0, sizeof(REAL4TimeVectorSeries) );


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
      LALSCreateVectorSequence( status->statusPtr,
                                &( waveform->h->data ), &in );
      CHECKSTATUSPTR(status);


      memcpy(waveform->f->data->data , ff->data, count*(sizeof(REAL4)));
      memcpy(waveform->a->data->data , a->data, 2*count*(sizeof(REAL4)));
      memcpy(waveform->phi->data->data ,phi->data, count*(sizeof(REAL8)));
      memcpy(waveform->h->data->data , h->data, 2*count*(sizeof(REAL4)));


      waveform->a->deltaT = waveform->f->deltaT = waveform->phi->deltaT
        = waveform->h->deltaT = 1./params->tSampling;

      waveform->a->sampleUnits = lalStrainUnit;
      waveform->f->sampleUnits = lalHertzUnit;
      waveform->phi->sampleUnits = lalDimensionlessUnit;
      waveform->h->sampleUnits = lalStrainUnit;
      waveform->position = ppnParams->position;
      waveform->psi = ppnParams->psi;


      snprintf( waveform->a->name,
	  	LALNameLength, "EOB inspiral amplitudes");
      snprintf( waveform->f->name,
		  LALNameLength, "EOB inspiral frequency");
      snprintf( waveform->phi->name,
	  	LALNameLength, "EOB inspiral phase");
      snprintf( waveform->h->name,
                LALNameLength, "EOB inspiral polarizations");

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
  LALSDestroyVector(status->statusPtr, &h);
  CHECKSTATUSPTR(status);

  /*on peut utiliser tSampling pour dfdt*/

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/* Engine function for generating waveform
   Craig Robinson 15/07/05 */
static void
LALEOBPPWaveformEngine (
                LALStatus        *status,
                REAL4Vector      *signalvec1,
                REAL4Vector      *signalvec2,
                REAL4Vector      *h,
                REAL4Vector      *a,
                REAL4Vector      *ff,
                REAL8Vector      *phi,
                UINT4            *countback,
                InspiralTemplate *params,
                InspiralInit     *paramsInit
                )
{


   UINT4                   count, nn=4, length = 0, hiSRndx=0, ndx=0, higherSR=0;

   REAL4Vector             *sig1, *sig2, *ampl, *freq;
   REAL8Vector             *phse;


   REAL8                   v2, eta, m, rn, r, rOld, s, p, q, dt, t, v, omega, f, ampl0;
   REAL8                   omegaOld;

   void                    *funcParams1, *funcParams2, *funcParams3;

   REAL8Vector             *values, *dvalues, *newvalues, *yt, *dym, *dyt;
   TofVIn                  in1;
   InspiralPhaseIn         in2;
   InspiralDerivativesIn   in3;
   rk4In                   in4;
   rk4GSLIntegrator        *integrator = NULL;
   pr3In                   pr3in;
   expnCoeffs              ak;
   expnFunc                func;
   rOfOmegaIn              rofomegain;
   DFindRootIn             rootIn1, rootIn2, rootIn3;

   REAL8 (*rOfOmegaFunc)(REAL8, void *); /* Function to be used in root finding later */

   /* Variables to allow the waveform to be generated */
   /* from a specific fLower */
   REAL8                   fCurrent;                /* The current frequency of the waveform */
   BOOLEAN                 writeToWaveform = 0;     /* Set to true when the current frequency
						     * crosses fLower */
   REAL8                   sInit, s0 = 0.0;  /* Initial phase, and phase to subtract */
   REAL8                   rmin = 20;        /* Smallest value of r at which to generate the waveform */
   COMPLEX16  MultSphHarmP;    /* Spin-weighted spherical harmonics */
   COMPLEX16  MultSphHarmM;    /* Spin-weighted spherical harmonics */
   COMPLEX16  hLM;             /* Factorized waveform */
   REAL4      x1, x2;
   UINT4      i, j, k, modeL;
   INT4       modeM;          /* number of modes required */
   REAL4      inclination;    /* binary inclination       */
   REAL4      coa_phase;      /* binary coalescence phase */
   REAL8      y_1, y_2, z1, z2; /* (2,2) and (2,-2) spherical harmonics needed in (h+,hx) */


   /* Used for EOBNR */
   COMPLEX8Vector *modefreqs;
   UINT4 resampFac;
   REAL8 resampEstimate;

   CHAR message[256];

   /* For checking XLAL return codes */
   INT4 xlalStatus;

   /* Variables used in injection */
   REAL8 unitHz;
   REAL8 cosI;/* cosine of system inclination */
   REAL8 apFac, acFac;/* extra factor in plus and cross amplitudes */

   /* Accuracy of root finding algorithms */
   const REAL8 xacc = 1.0e-12;

   REAL8 tStepBack; /* We need to step back 6M to attach ringdown */
   UINT4 nStepBack; /* Num points to step back */

   /* Keeping track of previous points so we can step back */
   REAL8Vector *pPrev = NULL;
   REAL8Vector *qPrev = NULL;
   REAL8Vector *rPrev = NULL;
   REAL8Vector *sPrev = NULL;

   /* Stuff at higher sample rate */
   REAL4Vector             *sig1Hi, *sig2Hi, *amplHi, *freqHi;
   REAL8Vector             *phseHi;
   UINT4                   lengthHiSR;

   /* Inidices of the ringdown matching points */
   /* peakIdx is the index where omega is a maximum */
   /* finalIdx is the index of the last point before the */
   /* integration breaks */
   UINT4Vector             *rdMatchPoint;
   UINT4                   peakIdx  = 0;
   UINT4                   finalIdx = 0;

   INITSTATUS(status, "LALEOBPPWaveformEngine", LALEOBPPWAVEFORMC);
   ATTATCHSTATUSPTR(status);


   ak   = paramsInit->ak;
   func = paramsInit->func;

   ASSERT(ak.totalmass/LAL_MTSUN_SI > 0.4, status,
   	LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   /* Check order is consistent */
   if ( params->order != LAL_PNORDER_PSEUDO_FOUR )
   {
     snprintf( message, 256, "Order must be LAL_PNORDER_PSEUDO_FOUR for approximant EOBNR_PP." );
     LALError( status, message );
     ABORT( status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE );
   }

   if (signalvec1) length = signalvec1->length; else if (ff) length = ff->length;

   /* Allocate some memory */
   values    = XLALCreateREAL8Vector( nn );
   dvalues   = XLALCreateREAL8Vector( nn );
   newvalues = XLALCreateREAL8Vector( nn );
   yt        = XLALCreateREAL8Vector( nn );
   dym       = XLALCreateREAL8Vector( nn );
   dyt       = XLALCreateREAL8Vector( nn );

   if ( !values || !dvalues || !newvalues || !yt || !dym || !dyt )
   {
     XLALDestroyREAL8Vector( values );
     XLALDestroyREAL8Vector( dvalues );
     XLALDestroyREAL8Vector( newvalues );
     XLALDestroyREAL8Vector( yt );
     XLALDestroyREAL8Vector( dym );
     XLALDestroyREAL8Vector( dyt );
     ABORTXLAL( status );
   }

   /* Set dt to sampling interval specified by user */
   dt =1./params->tSampling;
   eta = ak.eta;
   m = ak.totalmass;

   /* only used in injection case */
   unitHz = m*(REAL8)LAL_PI;
   cosI   = cos( params->inclination );
   apFac  = -2.0 * (1.0 + cosI*cosI) * eta * params->totalMass * LAL_MRSUN_SI/params->distance;
   acFac  = -4.0 * cosI * eta * params->totalMass * LAL_MRSUN_SI/params->distance;

   /* Set the amplitude depending on whether the distance is given */
   if ( params->distance > 0.0 )
     ampl0  = eta * params->totalMass * LAL_MRSUN_SI/params->distance;
   else
     ampl0  = params->signalAmplitude;

   /* Check we get a sensible answer */
   if ( ampl0 == 0.0 )
   {
     snprintf( message, 256, "Generating waveform of zero amplitude!!" );
     LALWarning( status, message );
   }

   /* Check that the 220 QNM freq. is less than the Nyquist freq. */
   /* Get QNM frequencies */
   modefreqs = XLALCreateCOMPLEX8Vector( 3 );
   xlalStatus = XLALGenerateQNMFreq( modefreqs, params, 2, 2, 3 );
   if ( xlalStatus != XLAL_SUCCESS )
   {
     XLALDestroyCOMPLEX8Vector( modefreqs );
     ABORTXLAL( status );
   }

   /* If Nyquist freq. <  220 QNM freq., exit */
   /* Note that we cancelled a factor of 2 occuring on both sides */
   if ( params->tSampling < modefreqs->data[0].re / LAL_PI )
   {
     XLALDestroyCOMPLEX8Vector( modefreqs );
     snprintf( message, 256, "Ringdown freq less than Nyquist freq. "
           "Increase sample rate or consider using EOB approximant.\n" );
     LALError(status->statusPtr, message);
     ABORT( status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
   }
   XLALDestroyCOMPLEX8Vector( modefreqs );

   /* Calculate the time we will need to step back for ringdown */
   tStepBack = 6.0 * params->totalMass * LAL_MTSUN_SI;
   nStepBack = ceil( tStepBack * params->tSampling );
   XLALPrintInfo( " We step back %d points\n", nStepBack );


   /* Allocate vectors to keep track of previous values */
   pPrev = XLALCreateREAL8Vector( nStepBack );
   qPrev = XLALCreateREAL8Vector( nStepBack );
   rPrev = XLALCreateREAL8Vector( nStepBack );
   sPrev = XLALCreateREAL8Vector( nStepBack );

   /* Calculate the resample factor for attaching the ringdown */
   /* We want it to be a power of 2 */
   resampEstimate = 20. / ( params->totalMass * LAL_MTSUN_SI * params->tSampling );
   resampFac = (UINT4) pow( 2.0, ceil(log2(resampEstimate) ) );

   /* The length of the vectors at the higher sample rate will be */
   /* the step back time plus the ringdown */
   lengthHiSR = ( nStepBack + (UINT4)(20.0 / modefreqs->data[0].im / dt) ) * resampFac;

   /* Find the initial velocity given the lower frequency */
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

   /* Then the initial phase */
   in2.v0 = ak.v0;
   in2.phi0 = params->startPhase;
   in2.dEnergy = func.dEnergy;
   in2.flux = func.flux;
   in2.coeffs = &ak;
   LALInspiralPhasing1(status->statusPtr, &s, v, &in2);
   CHECKSTATUSPTR(status);
/*
   LALInspiralPhasing1(v) gives the GW phase (= twice the orbital phase).
   The ODEs we solve give the orbital phase. Therefore, set the
   initial phase to be half the GW pahse.
*/
   s = s/2.;
   sInit = s;

   /* light ring value - where to stop evolution */
   rofomegain.eta = eta;
   rofomegain.omega = omega;
   rootIn1.xacc = xacc;
   rootIn3.xacc = xacc;
   rootIn1.xmax = 1.;
   rootIn1.xmin = 4.;
   rootIn2.xmax = 1000.;
   rootIn2.xmin = 3.;
   pr3in.eta = eta;
   pr3in.omegaS = params->OmegaS;
   pr3in.zeta2 = params->Zeta2;

   /* We will be changing the starting r if it is less than rmin */
   /* Therefore, we should reset pr3in.omega later if necessary. */
   /* For now we need it so that we can see what the initial r is. */

   pr3in.omega = omega;
   in3.totalmass = ak.totalmass;
   in3.dEnergy = func.dEnergy;
   in3.flux = func.flux;
   in3.coeffs = &ak;
   funcParams3 = (void *) &in3;

   funcParams1 = (void *) &rofomegain;

   switch (params->order)
   {
     case LAL_PNORDER_PSEUDO_FOUR:
       rootIn1.function = LALlightRingRadiusP4PN;
       rOfOmegaFunc = XLALrOfOmegaP4PN;
       funcParams2 = (void *) &pr3in;
       break;
     default:
       snprintf(message, 256, "There are no EOBNR_PP waveforms implemented at order %d\n", params->order);
       LALError( status, message );
       XLALDestroyREAL8Vector( values );
       XLALDestroyREAL8Vector( dvalues );
       XLALDestroyREAL8Vector( newvalues );
       XLALDestroyREAL8Vector( yt );
       XLALDestroyREAL8Vector( dym );
       XLALDestroyREAL8Vector( dyt );
       ABORT( status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
   }
   LALDBisectionFindRoot(status->statusPtr, &rn, &rootIn1, funcParams1);
   CHECKSTATUSPTR(status);
   r = XLALDBisectionFindRoot( rOfOmegaFunc, rootIn2.xmin,
              rootIn2.xmax, xacc, funcParams2);
   if ( XLAL_IS_REAL8_FAIL_NAN( r ) )
   {
     ABORTXLAL( status );
   }

   /* We want the waveform to generate from a point which won't cause */
   /* problems with the initial conditions. Therefore we force the code */
   /* to start at least at r = rmin (in units of M). */

   r = (r<rmin) ? rmin : r;

   rootIn3.xmax = 5;
   rootIn3.xmin = -10;
   pr3in.in3copy = in3;
   pr3in.r = r;

   /* Now that r is changed recompute omega corresponding */
   /* to that r and only then compute initial pr and pphi */

   switch (params->order)
   {
     case LAL_PNORDER_PSEUDO_FOUR:
       omegaofrP4PN (&omega, r, &pr3in);
       pr3in.omega = omega;
       LALpphiInitP4PN(&q, r, eta, params->OmegaS);
       rootIn3.function = LALprInitP4PN;
       /* first we compute vr (we need coeef->Fp6) */
       pr3in.q = q;
       funcParams2 = (void *) &pr3in;
       LALvrP4PN(&pr3in.vr, (void *) &pr3in);
       /* then we compute the initial value of p */
       LALDBisectionFindRoot(status->statusPtr, &p, &rootIn3, funcParams2);
       CHECKSTATUSPTR(status);
       if ( params->approximant == EOBNR_PF )
       {
         in4.function = LALHCapDerivativesP4PNFF;
       }
       else
       {
         in4.function = LALHCapDerivativesP4PN;
       }
       break;
     default:
       snprintf(message, 256, "There are no EOB/EOBNR waveforms implemented at order %d\n", params->order);
       LALError( status, message );
       XLALDestroyREAL8Vector( values );
       XLALDestroyREAL8Vector( dvalues );
       XLALDestroyREAL8Vector( newvalues );
       XLALDestroyREAL8Vector( yt );
       XLALDestroyREAL8Vector( dym );
       XLALDestroyREAL8Vector( dyt );
       ABORT( status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
   }

   values->data[0] = r;
   values->data[1] = s;
   values->data[2] = p;
   values->data[3] = q;
#if 0
   sprintf(message, "In EOB Initial values of r=%10.5e p=%10.5e q=%10.5e\n", r, p, q);
   LALInfo(status, message);
#endif

   in4.y = values;
   in4.h = dt/m;
   in4.n = nn;
   in4.yt = yt;
   in4.dym = dym;
   in4.dyt = dyt;

   /* Allocate memory for temporary arrays */
   sig1 = XLALCreateREAL4Vector ( length );
   sig2 = XLALCreateREAL4Vector ( length );
   ampl = XLALCreateREAL4Vector ( length*2 );
   freq = XLALCreateREAL4Vector ( length );
   phse = XLALCreateREAL8Vector ( length );

   if ( !sig1 || !sig2 || !ampl || !freq || !phse )
   {
     if ( sig1 ) XLALDestroyREAL4Vector( sig1 );
     if ( sig2 ) XLALDestroyREAL4Vector( sig2 );
     if ( ampl ) XLALDestroyREAL4Vector( ampl );
     if ( freq ) XLALDestroyREAL4Vector( freq );
     if ( phse ) XLALDestroyREAL8Vector( phse );
     XLALDestroyREAL8Vector( values );
     XLALDestroyREAL8Vector( dvalues );
     XLALDestroyREAL8Vector( newvalues );
     XLALDestroyREAL8Vector( yt );
     XLALDestroyREAL8Vector( dym );
     XLALDestroyREAL8Vector( dyt );
     ABORT( status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM );
   }

   memset(sig1->data, 0, sig1->length * sizeof( REAL4 ));
   memset(sig2->data, 0, sig2->length * sizeof( REAL4 ));
   memset(ampl->data, 0, ampl->length * sizeof( REAL4 ));
   memset(freq->data, 0, freq->length * sizeof( REAL4 ));
   memset(phse->data, 0, phse->length * sizeof( REAL8 ));

   /* And their higher sample rate counterparts */
   /* Allocate memory for temporary arrays */
   sig1Hi = XLALCreateREAL4Vector ( length );
   sig2Hi = XLALCreateREAL4Vector ( length );
   amplHi = XLALCreateREAL4Vector ( length*2 );
   freqHi = XLALCreateREAL4Vector ( length );
   phseHi = XLALCreateREAL8Vector ( length );

   if ( !sig1Hi || !sig2Hi || !amplHi || !freqHi || !phseHi )
   {
     if ( sig1 ) XLALDestroyREAL4Vector( sig1 );
     if ( sig2 ) XLALDestroyREAL4Vector( sig2 );
     if ( ampl ) XLALDestroyREAL4Vector( ampl );
     if ( freq ) XLALDestroyREAL4Vector( freq );
     if ( phse ) XLALDestroyREAL8Vector( phse );
     XLALDestroyREAL8Vector( values );
     XLALDestroyREAL8Vector( dvalues );
     XLALDestroyREAL8Vector( newvalues );
     XLALDestroyREAL8Vector( yt );
     XLALDestroyREAL8Vector( dym );
     XLALDestroyREAL8Vector( dyt );
     ABORT( status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM );
   }

   memset(sig1Hi->data, 0, sig1Hi->length * sizeof( REAL4 ));
   memset(sig2Hi->data, 0, sig2Hi->length * sizeof( REAL4 ));
   memset(amplHi->data, 0, amplHi->length * sizeof( REAL4 ));
   memset(freqHi->data, 0, freqHi->length * sizeof( REAL4 ));
   memset(phseHi->data, 0, phseHi->length * sizeof( REAL8 ));

   /* Initialize the GSL integrator */
   if (!(integrator = XLALRungeKutta4Init(nn, &in4)))
   {
     XLALDestroyREAL4Vector( sig1 );
     XLALDestroyREAL4Vector( sig2 );
     XLALDestroyREAL4Vector( ampl );
     XLALDestroyREAL4Vector( freq );
     XLALDestroyREAL8Vector( phse );
     XLALDestroyREAL8Vector( values );
     XLALDestroyREAL8Vector( dvalues );
     XLALDestroyREAL8Vector( newvalues );
     XLALDestroyREAL8Vector( yt );
     XLALDestroyREAL8Vector( dym );
     XLALDestroyREAL8Vector( dyt );
     ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
   }

   count = 0;
   if (a || signalvec2)
      params->nStartPad = 0; /* must be zero for templates and injection */

   count = params->nStartPad;

   /* Calculate the initial value of omega */
   in4.function(values, dvalues, funcParams3);
   omega = dvalues->data[1];

   /* Make sure old omega is < omega to get into the loop */
   omegaOld = omega - 0.1;

   /* Begin integration loop here */
   t = 0.0;
   rOld = r+0.1;

/*
   omegamatch = -0.05 -0.01 + 0.133 + 0.183 * params->eta + 1.161 * params->eta * params->eta;
*/

   while ( ( omega > omegaOld || !isnan(hLM.re) ) && r < rOld)
   {
      if (count > length)
      {
        XLALRungeKutta4Free( integrator );
        XLALDestroyREAL4Vector( sig1 );
        XLALDestroyREAL4Vector( sig2 );
        XLALDestroyREAL4Vector( ampl );
        XLALDestroyREAL4Vector( freq );
        XLALDestroyREAL8Vector( phse );
        XLALDestroyREAL8Vector( values );
        XLALDestroyREAL8Vector( dvalues );
        XLALDestroyREAL8Vector( newvalues );
        XLALDestroyREAL8Vector( yt );
        XLALDestroyREAL8Vector( dym );
        XLALDestroyREAL8Vector( dyt );
	ABORT(status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
      }

      rOld = r;

      if ( omega <= omegaOld && !peakIdx )
      {
        peakIdx = count - 1;
      }
        

      fCurrent = omega / (LAL_PI*m);
      if (!writeToWaveform)
      {
        s0 = s - sInit;
        if (r > rmin || fCurrent > f || fabs(fCurrent - f) < 1.0e-5)
        {
          writeToWaveform = 1;
        }
      }

      v = pow(omega, oneby3);
      v2 = v*v;

      if (writeToWaveform)
      {
	double st, amp;
	i = count;
	j = 2*count;
	k = j+1;
        st = 2.*(s - s0);
	amp = ampl0 * v2;
	/*--------------------------------------------------------
	   First we generate the real and imagninary parts of h22
	  --------------------------------------------------------*/
        xlalStatus = XLALGetFactorizedWaveform( &hLM, values, dvalues, &in3, 2, 2 );
        if ( xlalStatus == XLAL_FAILURE )
        {
          ABORTXLAL( status );
        }
        if ( !higherSR )
        {
         
	  sig1->data[i] =  (REAL4) ampl0 * hLM.re;
          sig2->data[i] =  (REAL4) ampl0 * hLM.im;
	  /*----------------------------------------------------------
	     ... then the frequency, amplitude of h+ and hx and phase
	    ----------------------------------------------------------*/
          freq->data[i] =  (REAL4)( omega );
          ampl->data[j] =  (REAL4)( apFac * v2 );
          ampl->data[k] =  (REAL4)( acFac * v2 );
          phse->data[i] =  (REAL8)( st );
        }
        else if ( !isnan( hLM.re) )
        {

          sig1Hi->data[i] =  (REAL4) ampl0 * hLM.re;
          sig2Hi->data[i] =  (REAL4) ampl0 * hLM.im;
          /*----------------------------------------------------------
             ... then the frequency, amplitude of h+ and hx and phase
            ----------------------------------------------------------*/
          freqHi->data[i] =  (REAL4)( omega );
          amplHi->data[j] =  (REAL4)( apFac * v2 );
          amplHi->data[k] =  (REAL4)( acFac * v2 );
          phseHi->data[i] =  (REAL8)( st );
        }
        else
        {
          finalIdx = --count;
        }
      }

      /* Integrate one step forward */
      in4.dydx = dvalues;
      in4.x = t/m;
      LALRungeKutta4(status->statusPtr, newvalues, integrator, funcParams3);
      BEGINFAIL( status )
      {
        XLALRungeKutta4Free( integrator );
        XLALDestroyREAL4Vector( sig1 );
        XLALDestroyREAL4Vector( sig2 );
        XLALDestroyREAL4Vector( ampl );
        XLALDestroyREAL4Vector( freq );
        XLALDestroyREAL8Vector( phse );
        XLALDestroyREAL8Vector( values );
        XLALDestroyREAL8Vector( dvalues );
        XLALDestroyREAL8Vector( newvalues );
        XLALDestroyREAL8Vector( yt );
        XLALDestroyREAL8Vector( dym );
        XLALDestroyREAL8Vector( dyt );
      }
      ENDFAIL( status );

      /* We need to track the dynamical variables prior to the current step */
      UINT4 startIndex = ndx < nStepBack ? nStepBack - ndx : 1;
      
      for ( i = startIndex; i < nStepBack; i++ )
      {
        rPrev->data[i-1] = rPrev->data[i];
        sPrev->data[i-1] = sPrev->data[i];
        pPrev->data[i-1] = pPrev->data[i];
        qPrev->data[i-1] = qPrev->data[i];
      }

      /* These are the current values of the dynamical variables */
      rPrev->data[nStepBack-1]=r;
      sPrev->data[nStepBack-1]=s;
      pPrev->data[nStepBack-1]=p;
      qPrev->data[nStepBack-1]=q;

      /* Update the values of the dynamical variables */
      r = values->data[0] = newvalues->data[0];
      s = values->data[1] = newvalues->data[1];
      p = values->data[2] = newvalues->data[2];
      q = values->data[3] = newvalues->data[3];

      /* Compute the derivaties at the new location */
      in4.function(values, dvalues, funcParams3);
      omegaOld = omega;
      omega = dvalues->data[1];


      /*----------------------------------------------------------------------*/
      /* We are going to terminate waveform generation if omega is greater    */
      /* than omegamatch - the frequency at which the ringdown is matched to  */
      /* merger waveform                                                      */
      /*----------------------------------------------------------------------*/
      if ( omega <= omegaOld && !higherSR )
      {
	/* We are now going to work with a higher sampling rate */
	/* Sometime in the future we might change code so that  */
	/* a higher sampling rate is used only if required */
	higherSR = 1;
        /*------------------------------------------------------------- */
	/* We are going to decrease the number of points by nStepBack+1 */
	/* In reality, note that we are really using the nStepBack      */
	/* points from the current step; the extra 1 is needed          */
	/* only because count is incremented before returning to the    */
	/* continuing the integration; the same is true with dt         */
        /*------------------------------------------------------------- */
	hiSRndx = count - nStepBack;
        count   = -1;
        t -= nStepBack * dt;
        dt /= (double) resampFac;
        t -= dt;
        in4.h = dt/m;

        r = values->data[0] = rPrev->data[0];
        s = values->data[1] = sPrev->data[0];
        p = values->data[2] = pPrev->data[0];
        q = values->data[3] = qPrev->data[0];

        /*----------------------------------------------------------------------*/
	/* Integration will stop if rOld is not reset to a value greater than r */
        /*----------------------------------------------------------------------*/

        rOld = r+0.1;

        in4.function(values, dvalues, funcParams3);
        omega = dvalues->data[1];

        /* Make sure omegaOld < omega so integration continues */
        omegaOld = omega - 0.1;
        fCurrent = omega/(LAL_PI*m);
      }

      if (writeToWaveform)
      {
	if (!higherSR)
          t = (++count-params->nStartPad) * dt;
	else
	{
	  t += dt;
	  count++;
	}

      }
      ndx++;
   }

   /*----------------------------------------------------------------------*/
   /* Record the final cutoff frequency of BD Waveforms for record keeping */
   /* ---------------------------------------------------------------------*/
   params->vFinal = v;
   if (signalvec1 && !signalvec2) params->tC = t;

   /* For now we just set fFinal to Nyquist */
   /* This is a hack to get the filtering code to work properly */
   params->fFinal = params->tSampling/2.;

   XLALRungeKutta4Free( integrator );
   XLALDestroyREAL8Vector( values );
   XLALDestroyREAL8Vector( dvalues );
   XLALDestroyREAL8Vector( newvalues );
   XLALDestroyREAL8Vector( yt );
   XLALDestroyREAL8Vector( dym );
   XLALDestroyREAL8Vector( dyt );
   XLALDestroyREAL8Vector( rPrev );
   XLALDestroyREAL8Vector( sPrev );
   XLALDestroyREAL8Vector( pPrev );
   XLALDestroyREAL8Vector( qPrev );


   /*--------------------------------------------------------------
    * Attach the ringdown waveform to the end of inspiral
     -------------------------------------------------------------*/
   REAL8 tmpSamplingRate = params->tSampling;
   params->tSampling *= resampFac;

   rdMatchPoint = XLALCreateUINT4Vector( 3 );

   /* Check the first matching point is sensible */
   if ( ceil( tStepBack / ( 2.0 * params->tSampling )) > peakIdx )
   {
     XLALPrintError( "Invalid index for first ringdown matching point.\n" );
     ABORT( status, LALINSPIRALH_ESIZE , LALINSPIRALH_MSGESIZE );
   }

   rdMatchPoint->data[0] = peakIdx - ceil( tStepBack * params->tSampling / 2.0 ) - 50;
   rdMatchPoint->data[1] = peakIdx - 50;
   rdMatchPoint->data[2] = finalIdx - 50;

   XLALPrintInfo( "Ringdown matching points: %e, %e\n", 
           (REAL8)rdMatchPoint->data[0]/resampFac + (REAL8)hiSRndx, 
           (REAL8)rdMatchPoint->data[1]/resampFac + (REAL8)hiSRndx );

   xlalStatus = XLALInspiralHybridAttachRingdownWave(sig1Hi, sig2Hi,
                   rdMatchPoint, params);
   if (xlalStatus != XLAL_SUCCESS )
   {
     XLALDestroyREAL4Vector( sig1 );
     XLALDestroyREAL4Vector( sig2 );
     XLALDestroyREAL4Vector( ampl );
     XLALDestroyREAL4Vector( freq );
     XLALDestroyREAL8Vector( phse );
     ABORTXLAL( status );
   }
   params->tSampling = tmpSamplingRate;
   count = hiSRndx;
   for(j=0; j<sig1Hi->length; j+=resampFac)
   {
     sig1->data[count] = sig1Hi->data[j];
     sig2->data[count] = sig2Hi->data[j];
     freq->data[count] = freqHi->data[j];
     if (sig1->data[count] == 0)
     {
       break;
     }
     count++;
   }
   *countback = count;

   /*-------------------------------------------------------------------
    * Compute the spherical harmonics required for constructing (h+,hx).
    * We are going to choose coa_phase to be zero. This perhaps should be
    * made compatible with the wave CoherentGW handles the phase at
    * coalecence. I have no idea how I (i.e., Sathya) might be able to
    * do this for EOBNR as there is no such thing as "phase at merger".
    -------------------------------------------------------------------*/
   inclination = (REAL4)params->inclination;
   coa_phase = 0.;
   /* -----------------------------------------------------------------
    * Attaching the (2,2) Spherical Harmonic
    * need some error checking
    *----------------------------------------*/
   modeL = 2;
   modeM = 2;
   xlalStatus = XLALSphHarm( &MultSphHarmP, modeL, modeM, inclination, coa_phase );
   if (xlalStatus != XLAL_SUCCESS )
   {
     XLALDestroyREAL4Vector( sig1 );
     XLALDestroyREAL4Vector( sig2 );
     XLALDestroyREAL4Vector( ampl );
     XLALDestroyREAL4Vector( freq );
     XLALDestroyREAL8Vector( phse );
     ABORTXLAL( status );
   }

   modeM = -2;
   xlalStatus = XLALSphHarm( &MultSphHarmM, modeL, modeM, inclination, coa_phase );
   if (xlalStatus != XLAL_SUCCESS )
   {
     XLALDestroyREAL4Vector( sig1 );
     XLALDestroyREAL4Vector( sig2 );
     XLALDestroyREAL4Vector( ampl );
     XLALDestroyREAL4Vector( freq );
     XLALDestroyREAL8Vector( phse );
     ABORTXLAL( status );
   }

   y_1 =   MultSphHarmP.re + MultSphHarmM.re;
   y_2 =   MultSphHarmM.im - MultSphHarmP.im;
   z1 = - MultSphHarmM.im - MultSphHarmP.im;
   z2 =   MultSphHarmM.re - MultSphHarmP.re;

#if 0
   sprintf(message, "MultSphHarm2,+2 re=%10.5e im=%10.5e\n", MultSphHarmP.re, MultSphHarmP.im);
   LALInfo(status, message);
   sprintf(message, "MultSphHarm2,-2 re=%10.5e im=%10.5e\n", MultSphHarmP.re, MultSphHarmM.im);
   LALInfo(status, message);
#endif

   /* Next, compute h+ and hx from h22, h22*, Y22, Y2-2 */
   for ( i = 0; i < sig1->length; i++)
   {
     freq->data[i] /= unitHz;
     x1 = sig1->data[i];
     x2 = sig2->data[i];
     sig1->data[i] = (x1 * y_1) + (x2 * y_2);
     sig2->data[i] = (x1 * z1) + (x2 * z2);
     if (x1 || x2)
     {
   /*
    * If the ringdown modes were added then artificially increase the phasing so
    * that it is nonzero until the end of the ringdown. When ringdown modes are
    * added the phase information is not used in injetions, only hplus and hcross
    * and therefore it shouldn't matter what phasing is as long as it is nonzero.
    */
       if ( i >= hiSRndx )
       {
	 phse->data[i] = phse->data[i-1]+LAL_PI/20.;
       }
     }
   }

   /*------------------------------------------------------
    * If required by the user copy other data sets to the
    * relevant arrays
    ------------------------------------------------------*/
   if (h)
   {
     for(i = 0; i < length; i++)
     {
       j = 2*i;
       k = j+1;
       h->data[j] = sig1->data[i];
       h->data[k] = sig2->data[i];
     }
   }
   if (signalvec1) memcpy(signalvec1->data , sig1->data, length * (sizeof(REAL4)));
   if (signalvec2) memcpy(signalvec2->data , sig2->data, length * (sizeof(REAL4)));
   if (ff)         memcpy(ff->data      , freq->data, length * (sizeof(REAL4)));
   if (a)          memcpy(a->data       , ampl->data, 2*length*(sizeof(REAL4)));
   if (phi)        memcpy(phi->data     , phse->data, length * (sizeof(REAL8)));

#if 0
   sprintf(message, "fFinal=%10.5e count=%d\n", params->fFinal, *countback);
   LALInfo(status, message);
#endif

   /* Clean up */
   XLALDestroyREAL4Vector ( sig1 );
   XLALDestroyREAL4Vector ( sig2 );
   XLALDestroyREAL4Vector ( ampl );
   XLALDestroyREAL4Vector ( freq );
   XLALDestroyREAL8Vector ( phse );


   DETATCHSTATUSPTR(status);
   RETURN(status);
}
