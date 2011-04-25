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
#include <lal/LALEOBNRv2Waveform.h>
#include <lal/LALAdaptiveRungeKutta4.h>
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
  EOBACoefficients *aCoeffs;
  InspiralDerivativesIn in3copy;
} pr3In;


static inline REAL8
XLALCalculateA5( REAL8 eta );

static inline REAL8
XLALCalculateA6( REAL8 eta );

static REAL8
omegaofrP4PN (
             const REAL8 r,
             const REAL8 eta,
             void *params);

static REAL8
nonKeplerianCoefficient(
                   REAL8Vector * restrict values,
                   const REAL8       eta,
                   EOBACoefficients *coeffs );

static
int LALHCapDerivativesP4PN(   double t,
                               const double values[],
                               double dvalues[],
                               void   *funcParams);

static
REAL8 XLALCalculateOmega (   REAL8 eta,
                             REAL8 r,
                             REAL8 pr,
                             REAL8 pPhi,
                             EOBACoefficients *aCoeffs );

static
int XLALFirstStoppingCondition(double t,
                               const double values[],
                               double dvalues[],
                               void   *funcParams);

static
int XLALHighSRStoppingCondition(double t,
                                const double values[],
                                double dvalues[],
                                void   *funcParams);

static
REAL8 XLALCalculateEOBA( const REAL8 r,
                         EOBACoefficients * restrict coeffs );

static
REAL8 XLALCalculateEOBD( REAL8    r,
                         REAL8	eta);

static
REAL8 XLALCalculateEOBdAdr( const REAL8 r,
                            EOBACoefficients * restrict coeffs );

static
REAL8 XLALEffectiveHamiltonian( const REAL8 eta,
                                const REAL8 r,
                                const REAL8 pr,
                                const REAL8 pp,
                                EOBACoefficients *aCoeffs );

static
void LALprInitP4PN(LALStatus *status, REAL8 *pr , REAL8 , void  *params);

static
REAL8 LALpphiInitP4PN( const REAL8 r,
                       EOBACoefficients * restrict coeffs );

static
REAL8 XLALrOfOmegaP4PN (REAL8 r, void *params);

static
REAL8 LALvrP4PN(const REAL8 r, const REAL8 omega, pr3In *params);

static void
LALEOBPPWaveformEngine (
                LALStatus        *status,
                REAL4Vector      *signalvec1,
                REAL4Vector      *signalvec2,
                REAL4Vector      *h,
                UINT4            *countback,
                InspiralTemplate *params,
                InspiralInit     *paramsInit
                );

NRCSID (LALEOBPPWAVEFORMC,
"$Id$");


INT4 XLALGetFactorizedWaveform( COMPLEX16             * restrict hlm,
				REAL8Vector           * restrict values,
                                const REAL8           v,
                                const INT4            l,
                                const INT4            m,
                                EOBParams             * restrict params
                                )
{
	static const char func[] = "XLALGetFactorizedWaveform";

        /* Status of function calls */
        INT4 status;
        INT4 i;

        REAL8 eta;	
	REAL8 r, pr, pp, Omega, v2, vh, vh3, k, hathatk, eulerlogxabs;
	REAL8 Hreal, Heff, Slm, deltalm, rholm, rholmPwrl;
	COMPLEX16 Tlm;
        COMPLEX16 hNewton;
	gsl_sf_result lnr1, arg1, z2;

        /* Non-Keplerian velocity */
        REAL8 vPhi;

        /* Non-quasicircular correction in dynamics */
        REAL8 a1, a2, a3;
        REAL8 hNQC = 1.0;

        /* Pre-computed coefficients */
        FacWaveformCoeffs *hCoeffs = params->hCoeffs;

	if ( abs(m) > (INT4) l )
	{
	  XLAL_ERROR( func, XLAL_EINVAL );
	}
	

        eta = params->eta;

        /* Check our eta was sensible */
        if ( eta > 0.25 )
        {
          XLALPrintError("Eta seems to be > 0.25 - this isn't allowed!\n" );
          XLAL_ERROR( func, XLAL_EINVAL );
        }
        else if ( eta == 0.25 && m % 2 )
        {
          /* If m is odd and dM = 0, hLM will be zero */
          memset( hlm, 0, sizeof( COMPLEX16 ) );
          return XLAL_SUCCESS;
        }
        
	r	= values->data[0];
        pr	= values->data[2];
	pp	= values->data[3];

	Heff	= XLALEffectiveHamiltonian( eta, r, pr, pp, params->aCoeffs ); 
	Hreal	= sqrt( 1.0 + 2.0 * eta * ( Heff - 1.0) );
	v2	= v * v;
        Omega   = v2 * v;
        vh3     = Hreal * Omega;
	vh	= cbrt(vh3);
	eulerlogxabs = LAL_GAMMA + log( 2.0 * (REAL8)m * v );


        /* Calculate the non-Keplerian velocity */
        vPhi = nonKeplerianCoefficient( values, eta, params->aCoeffs );

        vPhi  = r * cbrt(vPhi);
        vPhi *= Omega;

        /* Calculate the newtonian multipole */
        status = XLALCalculateNewtonianMultipole( &hNewton, vPhi * vPhi, r,
                         values->data[1], (UINT4)l, m, params );
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
	XLAL_CALLGSL( status = gsl_sf_fact_e( l, &z2 ) );
	if ( status != GSL_SUCCESS)
	{
	  XLALPrintError("Error in GSL function\n" );
	  XLAL_ERROR( func, XLAL_EFUNC );
	}
	Tlm = XLALCOMPLEX16Exp( XLALCOMPLEX16Rect( lnr1.val + LAL_PI * hathatk, 
				arg1.val + 2.0 * hathatk * log(4.0*k/sqrt(LAL_E)) ) );
        Tlm = XLALCOMPLEX16DivReal( Tlm, z2.val );

        /* Calculate the residue phase and amplitude terms */
	switch( l )
	{
	  case 2:
	    switch( abs(m) )
	    {
	      case 2:

                /* For 2,2 mode, calculate the NQC correction */
                a1 = -4.55919 + 18.761 * eta - 24.226 * eta*eta;
                a2 = 37.683 - 201.468 * eta + 324.591 * eta*eta;
                a3 = - 39.6024 + 228.899 * eta - 387.222 * eta * eta;

                hNQC = 1. + (pr*pr / (r*r*Omega*Omega)) * ( a1
                       + a2 / r + a3 / (r*sqrt(r)) );

	        deltalm = vh3*(hCoeffs->delta22vh3 + vh3*(hCoeffs->delta22vh6 
			+ vh*vh*(hCoeffs->delta22vh8 + hCoeffs->delta22vh9*vh))) 
			+ hCoeffs->delta22v5 *v*v2*v2;
		rholm	= 1. + v2*(hCoeffs->rho22v2 + v*(hCoeffs->rho22v3
			+ v*(hCoeffs->rho22v4
			+ v*(hCoeffs->rho22v5 + v*(hCoeffs->rho22v6 
			+ hCoeffs->rho22v6l*eulerlogxabs + v*(hCoeffs->rho22v7 
			+ v*(hCoeffs->rho22v8 + hCoeffs->rho22v8l*eulerlogxabs 
			+ (hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs)*v2)))))));
	        break;
	      case 1:
                {
	        deltalm = vh3*(hCoeffs->delta21vh3 + vh3*(hCoeffs->delta21vh6
			+ vh*(hCoeffs->delta21vh7 + (hCoeffs->delta21vh9)*vh*vh))) 
			+ hCoeffs->delta21v5*v*v2*v2;
		rholm	= 1. + v*(hCoeffs->rho21v1
			+ v*( hCoeffs->rho21v2 + v*(hCoeffs->rho21v3 + v*(hCoeffs->rho21v4 
			+ v*(hCoeffs->rho21v5 + v*(hCoeffs->rho21v6 + hCoeffs->rho21v6l*eulerlogxabs 
			+ v*(hCoeffs->rho21v7 + hCoeffs->rho21v7l * eulerlogxabs 
			+ v*(hCoeffs->rho21v8 + hCoeffs->rho21v8l * eulerlogxabs 
			+ (hCoeffs->rho21v10 + hCoeffs->rho21v10l * eulerlogxabs)*v2))))))));
                }
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
	        deltalm = vh3*(hCoeffs->delta33vh3 + vh3*(hCoeffs->delta33vh6 + hCoeffs->delta33vh9*vh3)) 
                        + hCoeffs->delta33v5*v*v2*v2;
		rholm	= 1. + v2*(hCoeffs->rho33v2 + v*(hCoeffs->rho33v3 + v*(hCoeffs->rho33v4 
			+ v*(hCoeffs->rho33v5 + v*(hCoeffs->rho33v6 + hCoeffs->rho33v6l*eulerlogxabs
			+ v*(hCoeffs->rho33v7 + (hCoeffs->rho33v8 + hCoeffs->rho33v8l*eulerlogxabs)*v))))));
	        break;
	      case 2:
		deltalm = vh3*(hCoeffs->delta32vh3 + vh*(hCoeffs->delta32vh4 + vh*vh*(hCoeffs->delta32vh6
			+ hCoeffs->delta32vh9*vh3)));
		rholm	= 1. + v*(hCoeffs->rho32v 
			+ v*(hCoeffs->rho32v2 + v*(hCoeffs->rho32v3 + v*(hCoeffs->rho32v4 + v*(hCoeffs->rho32v5
			+ v*(hCoeffs->rho32v6 + hCoeffs->rho32v6l*eulerlogxabs
			+ (hCoeffs->rho32v8 + hCoeffs->rho32v8l*eulerlogxabs)*v2))))));
		break;
	      case 1:
		deltalm = vh3*(hCoeffs->delta31vh3 + vh3*(hCoeffs->delta31vh6
			+ vh*(hCoeffs->delta31vh7 + hCoeffs->delta31vh9*vh*vh))) 
			+ hCoeffs->delta31v5*v*v2*v2;
		rholm	= 1. + v2*(hCoeffs->rho31v2 + v*(hCoeffs->rho31v3 + v*(hCoeffs->rho31v4 
			+ v*(hCoeffs->rho31v5 + v*(hCoeffs->rho31v6 + hCoeffs->rho31v6l*eulerlogxabs 
			+ v*(hCoeffs->rho31v7 + (hCoeffs->rho31v8 + hCoeffs->rho31v8l*eulerlogxabs)*v))))));
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
                
	        deltalm = vh3*(hCoeffs->delta44vh3 + hCoeffs->delta44vh6 *vh3);
		rholm	= 1. + v2*(hCoeffs->rho44v2
			+ v*( hCoeffs->rho44v3 + v*(hCoeffs->rho44v4
			+ v*(hCoeffs->rho44v5 + (hCoeffs->rho44v6
			+ hCoeffs->rho44v6l*eulerlogxabs)*v))));
	        break;
	      case 3:
	        deltalm = vh3*(hCoeffs->delta43vh3 + vh*(hCoeffs->delta43vh4 
			+ hCoeffs->delta43vh6*vh*vh));
		rholm	= 1. + v*(hCoeffs->rho43v
			+ v*(hCoeffs->rho43v2
			+ v2*(hCoeffs->rho43v4 + v*(hCoeffs->rho43v5
			+ (hCoeffs->rho43v6 + hCoeffs->rho43v6l*eulerlogxabs)*v))));
	        break;
	      case 2:
		deltalm = vh3*(hCoeffs->delta42vh3 + hCoeffs->delta42vh6*vh3);
		rholm	= 1. + v2*(hCoeffs->rho42v2
			+ v*(hCoeffs->rho42v3 + v*(hCoeffs->rho42v4 + v*(hCoeffs->rho42v5
			+ (hCoeffs->rho42v6 + hCoeffs->rho42v6l*eulerlogxabs)*v))));
		break;
	      case 1:
		deltalm = vh3*(hCoeffs->delta41vh3 + vh*(hCoeffs->delta41vh4
			+ hCoeffs->delta41vh6*vh*vh));
		rholm	= 1. + v*(hCoeffs->rho41v 
			+ v*(hCoeffs->rho41v2
			+ v2*(hCoeffs->rho41v4 + v*(hCoeffs->rho41v5 
			+ (hCoeffs->rho41v6 +  hCoeffs->rho41v6l*eulerlogxabs)*v))));
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
	        deltalm = hCoeffs->delta55vh3*vh3;
		rholm	= 1. + v2*( hCoeffs->rho55v2 
			+ v*(hCoeffs->rho55v3 + v*(hCoeffs->rho55v4 + hCoeffs->rho55v5*v)));
	        break;
	      case 4:
		deltalm = vh3*(hCoeffs->delta54vh3 + hCoeffs->delta54vh4*vh);
		rholm	= 1. + v2*(hCoeffs->rho54v2 + v*(hCoeffs->rho54v3
			+ hCoeffs->rho54v4*v));
		break;
	      case 3:
	        deltalm = hCoeffs->delta53vh3 * vh3;
		rholm	= 1. + v2*(hCoeffs->rho53v2 
			+ v*(hCoeffs->rho53v3 + v*(hCoeffs->rho53v4 + hCoeffs->rho53v5*v)));
	        break;
	      case 2:
		deltalm = vh3*(hCoeffs->delta52vh3 + hCoeffs->delta52vh4*vh);
		rholm	= 1. + v2*(hCoeffs->rho52v2 + v*(hCoeffs->rho52v3
			+ hCoeffs->rho52v4*v));
		break;
	      case 1:
		deltalm = hCoeffs->delta51vh3 * vh3;
		rholm	= 1. + v2*(hCoeffs->rho51v2 
			+ v*(hCoeffs->rho51v3 + v*(hCoeffs->rho51v4 + hCoeffs->rho51v5*v)));
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
	        deltalm = hCoeffs->delta66vh3*vh3;
		rholm	= 1. + v2*(hCoeffs->rho66v2 + v*(hCoeffs->rho66v3
			+ hCoeffs->rho66v4*v));
	        break;
	      case 5:
		deltalm = hCoeffs->delta65vh3*vh3;
		rholm	= 1. + v2*(hCoeffs->rho65v2 + hCoeffs->rho65v3*v);
		break;
	      case 4:
		deltalm = hCoeffs->delta64vh3 * vh3;
		rholm	= 1. + v2*(hCoeffs->rho64v2 + v*(hCoeffs->rho64v3
			+ hCoeffs->rho64v4*v));
		break;
	      case 3:
	        deltalm = hCoeffs->delta63vh3 * vh3;
		rholm	= 1. + v2*(hCoeffs->rho63v2 + hCoeffs->rho63v3*v);
	        break;
	      case 2:
		deltalm = hCoeffs->delta62vh3 * vh3;
		rholm	= 1. + v2*(hCoeffs->rho62v2 + v*(hCoeffs->rho62v3
			+ hCoeffs->rho62v4 * v));
		break;
	      case 1:
		deltalm = hCoeffs->delta61vh3 * vh3;
		rholm	= 1. + v2*(hCoeffs->rho61v2 + hCoeffs->rho61v3*v);
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
	        deltalm = hCoeffs->delta77vh3 * vh3;
		rholm	= 1. + v2*(hCoeffs->rho77v2 + hCoeffs->rho77v3 * v);
	        break;
	      case 6:
	        deltalm = 0.0;
		rholm	= 1. + hCoeffs->rho76v2 * v2;
	        break;
	      case 5:
		deltalm = hCoeffs->delta75vh3 * vh3;
		rholm	= 1. + v2*(hCoeffs->rho75v2 + hCoeffs->rho75v3*v);
		break;
	      case 4:
		deltalm = 0.0;
		rholm	= 1. + hCoeffs->rho74v2 * v2;
		break;
	      case 3:
	        deltalm = hCoeffs->delta73vh3 *vh3;
		rholm	= 1. + v2*(hCoeffs->rho73v2 + hCoeffs->rho73v3 * v);
	        break;
	      case 2:
		deltalm = 0.0;
		rholm	= 1. + hCoeffs->rho72v2 * v2;
		break;
	      case 1:
		deltalm = hCoeffs->delta71vh3 * vh3;
		rholm	= 1. + v2*(hCoeffs->rho71v2 +hCoeffs->rho71v3 * v);
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
		rholm	= 1. + hCoeffs->rho88v2 * v2;
	        break;
	      case 7:
		deltalm = 0.0;
		rholm	= 1. + hCoeffs->rho87v2 * v2;
		break;
	      case 6:
	        deltalm = 0.0;
		rholm	= 1. + hCoeffs->rho86v2 * v2;
	        break;
	      case 5:
		deltalm = 0.0;
		rholm	= 1. + hCoeffs->rho85v2 * v2;
		break;
	      case 4:
		deltalm = 0.0;
		rholm	= 1. + hCoeffs->rho84v2 * v2;
		break;
	      case 3:
	        deltalm = 0.0;
		rholm	= 1. + hCoeffs->rho83v2 * v2;
	        break;
	      case 2:
		deltalm = 0.0;
		rholm	= 1. + hCoeffs->rho82v2 * v2;
		break;
	      case 1:
		deltalm = 0.0;
		rholm	= 1. + hCoeffs->rho81v2 * v2;
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

        /* Raise rholm to the lth power */
        rholmPwrl = 1.0;
        i = l;
        while ( i-- )
        {
          rholmPwrl *= rholm;
        }

	*hlm = XLALCOMPLEX16MulReal( XLALCOMPLEX16Mul( Tlm, XLALCOMPLEX16Polar( 1.0, deltalm) ), 
				     Slm*rholmPwrl*hNQC );
        *hlm = XLALCOMPLEX16Mul( *hlm, hNewton );

	return XLAL_SUCCESS;
}

static inline
REAL8 XLALCalculateA5( const REAL8 eta )
{
  return - 5.82827 - 143.486 * eta + 447.045 * eta * eta;
}

static inline 
REAL8 XLALCalculateA6( const REAL8 UNUSED eta )
{
  return 184.0;
}


/**
 * Function to pre-compute the coefficients in the EOB A potential function
 */

static
int XLALCalculateEOBACoefficients(
          EOBACoefficients * const coeffs,
          const REAL8              eta
          )
{
  REAL8 eta2, eta3;
  REAL8 a4, a5, a6;

  eta2 = eta*eta;
  eta3 = eta2 * eta;

  /* Note that the definitions of a5 and a6 DO NOT correspond to those in the paper */
  /* Therefore we have to multiply the results of our a5 and a6 finctions by eta. */

  a4 = ninty4by3etc * eta;
  a5 = XLALCalculateA5( eta ) * eta;
  a6 = XLALCalculateA6( eta ) * eta;

  coeffs->n4 =  -64. + 12.*a4 + 4.*a5 + a6 + 64.*eta - 4.*eta2;
  coeffs->n5 = 32. -4.*a4 - a5 - 24.*eta;
  coeffs->d0 = 4.*a4*a4 + 4.*a4*a5 + a5*a5 - a4*a6 + 16.*a6 
             + (32.*a4 + 16.*a5 - 8.*a6) * eta + 4.*a4*eta2 + 32.*eta3;
  coeffs->d1 = 4.*a4*a4 + a4*a5 + 16.*a5 + 8.*a6 + (32.*a4 - 2.*a6)*eta + 32.*eta2 + 8.*eta3;
  coeffs->d2 = 16.*a4 + 8.*a5 + 4.*a6 + (8.*a4 + 2.*a5)*eta + 32.*eta2;
  coeffs->d3 = 8.*a4 + 4.*a5 + 2.*a6 + 32.*eta - 8.*eta2;
  coeffs->d4 = 4.*a4 + 2.*a5 + a6 + 16.*eta - 4.*eta2;
  coeffs->d5 = 32. - 4.*a4 - a5 - 24. * eta;

  return XLAL_SUCCESS;
}

static
REAL8 XLALCalculateEOBA( const REAL8 r,
                         EOBACoefficients * restrict coeffs )
{

  REAL8 r2, r3, r4, r5;
  REAL8 NA, DA;

  /* Note that this function uses pre-computed coefficients,
   * and assumes they have been calculated. Since this is a static function,
   * so only used here, I assume it is okay to neglect error checking 
   */

  r2 = r*r;
  r3 = r2 * r;
  r4 = r2*r2;
  r5 = r4*r;


  NA = r4 * coeffs->n4
     + r5 * coeffs->n5;

  DA = coeffs->d0
     + r  * coeffs->d1
     + r2 * coeffs->d2
     + r3 * coeffs->d3
     + r4 * coeffs->d4
     + r5 * coeffs->d5;

  return NA/DA;
}
  

static
REAL8 XLALCalculateEOBdAdr( const REAL8 r,
                            EOBACoefficients * restrict coeffs )
{
  REAL8 r2, r3, r4, r5;

  REAL8 NA, DA, dNA, dDA, dA;

  r2 = r*r;
  r3 = r2 * r;
  r4 = r2*r2;
  r5 = r4*r;

  NA = r4 * coeffs->n4
     + r5 * coeffs->n5;

  DA = coeffs->d0
     + r  * coeffs->d1
     + r2 * coeffs->d2
     + r3 * coeffs->d3
     + r4 * coeffs->d4
     + r5 * coeffs->d5;

  dNA = 4. * coeffs->n4 * r3
      + 5. * coeffs->n5 * r4;

  dDA = coeffs->d1
      + 2. * coeffs->d2 * r
      + 3. * coeffs->d3 * r2
      + 4. * coeffs->d4 * r3
      + 5. * coeffs->d5 * r4;

  dA = dNA * DA - dDA * NA;

  return dA / (DA*DA);
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


/**
 * Function to calculate the EOB effective Hamiltonian for the
 * given values of the dynamical variables. The coefficients in the
 * A potential function should already have been computed.
 * Note that the pr used here is the tortoise co-ordinate.
 */
static
REAL8 XLALEffectiveHamiltonian( const REAL8 eta,
                                const REAL8 r,
                                const REAL8 pr,
                                const REAL8 pp,
                                EOBACoefficients *aCoeffs )
{

        /* The pr used in here is the tortoise co-ordinate */
        REAL8 r2, pr2, pp2, z3, eoba;

        r2   = r * r;
        pr2  = pr * pr;
        pp2  = pp * pp;

        eoba = XLALCalculateEOBA( r, aCoeffs );
        z3   = 2. * ( 4. - 3. * eta ) * eta;
        return sqrt( pr2 + eoba * ( 1.  + pp2/r2 + z3*pr2*pr2/r2 ) );
}

/*-------------------------------------------------------------------*/
/*                      pseudo-4PN functions                         */
/*-------------------------------------------------------------------*/

REAL8
LALpphiInitP4PN(
            const REAL8 r,
            EOBACoefficients * restrict coeffs
            )
{
  REAL8 u, u2, A, dA;

  u  = 1. / r;
  u2 = u*u; 

  /* Comparing this expression with the old one, I *think* I will need to multiply */
  /* dA by - r^2. TODO: Check that this should be the case */

  A  = XLALCalculateEOBA( r, coeffs );
  dA = - r*r * XLALCalculateEOBdAdr( r, coeffs );
  return sqrt(-dA/(2.*u*A + u2 * dA));
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
  REAL8   u, u2, u3, u4, p2, p3, p4, q2, A;
  REAL8  onebyD, AbyD, Heff, HReal, etahH;
  REAL8 eta, eta2, z3, r, vr, q;
  pr3In *ak;

  /* TODO: Change this to use tortoise coord */

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

   A = XLALCalculateEOBA( r, ak->aCoeffs );
   onebyD = 1. / XLALCalculateEOBD( r, eta );
   AbyD = A * onebyD;

   Heff = sqrt(A*(1. + AbyD * p2 + q*q * u2 + z3 * p4 * u2));
   HReal = sqrt(1. + 2.*eta*(Heff - 1.)) / eta;
   etahH = eta*Heff*HReal;

   *pr = -vr +  A*(AbyD*p + 2. * z3 * u2 * p3)/etahH;
/* This sets pr = dH/dpr - vr, calls rootfinder,
   gets value of pr s.t. dH/pr = vr */
   RETURN(status);
}


/*-------------------------------------------------------------------*/

static REAL8
omegaofrP4PN (
             const REAL8 r,
             const REAL8 eta,
             void *params)
{
   REAL8 u, u3, A, dA, omega;

   EOBACoefficients *coeffs;
   coeffs = (EOBACoefficients *) params;

   u = 1./r;
   u3 = u*u*u;

   A  = XLALCalculateEOBA( r, coeffs );
   dA = XLALCalculateEOBdAdr( r, coeffs );

   /* Comparing this expression with the old one, I *think* I will need to multiply */
   /* dA by - r^2. TODO: Check that this should be the case */
   dA = - r * r * dA;

   omega = sqrt(u3) * sqrt ( -0.5 * dA /(1. + 2.*eta * (A/sqrt(A+0.5 * u*dA)-1.)));
   return omega;
}

static REAL8
nonKeplerianCoefficient(
                   REAL8Vector * restrict values,
                   const REAL8       eta,
                   EOBACoefficients *coeffs )
{

  REAL8 r    = values->data[0];
  REAL8 pphi = values->data[3];

  REAL8 A  = XLALCalculateEOBA( r, coeffs );
  REAL8 dA = XLALCalculateEOBdAdr( r, coeffs );

  return 2. * (1. + 2. * eta * ( -1. + sqrt( (1. + pphi*pphi/(r*r)) * A ) ) )
          / ( r*r * dA );
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
  omega2 = omegaofrP4PN(r, pr3in->eta, pr3in->aCoeffs);
  return ( -omega1 + omega2 );
}

/*-------------------------------------------------------------------*/

/* Version which uses the factorized flux */
int
LALHCapDerivativesP4PN( double UNUSED t,
                        const REAL8 values[],
                        REAL8       dvalues[],
                        void        *funcParams
                      )
{

  EOBParams *params = NULL;

  /* Max l to sum up to in the factorized flux */
  const INT4 lMax = 8;

  REAL8 eta, omega;

  double r, s, p, q;
  double dr, ds, dp, dq;
  REAL8 r2, p2, p4, q2;
  REAL8 u, u2, u3;
  REAL8 A, AoverSqrtD, dAdr, Heff, Hreal;
  REAL8 HeffHreal;
  REAL8 flux;
  REAL8 z3;

  /* Factorized flux function takes a REAL8Vector */
  /* We will wrap our current array for now */
  REAL8Vector valuesVec;

  params = (EOBParams *) funcParams;

  valuesVec.length = 4;
  memcpy( &(valuesVec.data), &(values), sizeof(REAL8 *) );

  eta = params->eta;

  z3   = 2. * ( 4. - 3. * eta ) * eta;

  r = values[0];
  s = values[1];
  p = values[2];
  q = values[3];

  u  = 1.0 / r;
  u2 = u * u;
  u3 = u2 * u;
  r2 = r * r;
  p2 = p*p;
  p4 = p2 * p2;
  q2 = q * q;

  A          = XLALCalculateEOBA(r, params->aCoeffs );
  dAdr       = XLALCalculateEOBdAdr(r, params->aCoeffs );
  AoverSqrtD = A / sqrt( XLALCalculateEOBD(r, eta) );

  /* Note that Hreal as given here is missing a factor of 1/eta */
  /* This is because it only enters into the derivatives in     */
  /* the combination eta*Hreal*Heff, so the eta would get       */
  /* cancelled out anyway. */

  Heff  = XLALEffectiveHamiltonian( eta, r, p, q, params->aCoeffs );
  Hreal = sqrt( 1. + 2.*eta*(Heff - 1.) );

  HeffHreal = Heff * Hreal;

  dr = dvalues[0] = AoverSqrtD * u2 * p * (r2 + 2. * p2 * z3 * A ) / HeffHreal;
  ds = dvalues[1] = omega = q * A * u2 / HeffHreal;

  /* Note that the only field of dvalues used in the flux is dvalues->data[1] */
  /* which we have already calculated. */
  flux = XLALInspiralFactorizedFlux( &valuesVec, omega, params, lMax );

  dp = dvalues[2] = 0.5 * AoverSqrtD * u3 * (  2.0 * ( q2 + p4 * z3) * A
                     - r * ( q2 + r2 + p4 * z3 ) * dAdr ) / HeffHreal
                     - (p / q) * (flux / (eta * omega));

  dq = flux;
  /* This function can fail */
  /* TODO: Implement proper error checking */

  dq = dvalues[3] = - dq / (eta * omega);

  return GSL_SUCCESS;
}

/* Function for calculating only omega */
static
REAL8 XLALCalculateOmega(   REAL8 eta,
                            REAL8 r,
                            REAL8 pr,
                            REAL8 pPhi,
                            EOBACoefficients *aCoeffs )
{

  REAL8 A, Heff, Hreal, HeffHreal;

  A    = XLALCalculateEOBA( r, aCoeffs );
  Heff = XLALEffectiveHamiltonian( eta, r, pr, pPhi, aCoeffs );
  Hreal = sqrt( 1. + 2.*eta*(Heff - 1.) );

  HeffHreal = Heff * Hreal;
  return pPhi * A / (HeffHreal*r*r);
}

/**
 * Function which will calculate the stopping condition for the 
 * initial sampling rate
 */
static int
XLALFirstStoppingCondition(double UNUSED t,
                           const double UNUSED values[],
                           double dvalues[],
                           void *funcParams
                          )
{

  EOBParams *params = (EOBParams *)funcParams;
  double omega = dvalues[1];

  if ( omega < params->omega )
  {
    return 1;
  }

  params->omega = omega;
  return GSL_SUCCESS;
}

/**
 * Function which will calculate the stopping condition for the 
 * initial sampling rate
 */
static int
XLALHighSRStoppingCondition(double UNUSED t,
                           const double values[],
                           double UNUSED dvalues[],
                           void UNUSED *funcParams
                          )
{

  if ( values[0] <= 1.0 )
  {
    return 1;
  }
  return 0;
}

/*-------------------------------------------------------------------*/
REAL8 LALvrP4PN( const REAL8 r,
                 const REAL8 omega,
                 pr3In *params )
{
  REAL8 A, dAdr, d2Adr2, dA, d2A, NA, DA, dDA, dNA, d2DA, d2NA;
  REAL8 r2, r3, r4, r5, u, u2, u3, u4, v, x1;
  REAL8 eta, FDIS;
  REAL8 twoUAPlusu2dA;

  EOBACoefficients *aCoeffs = params->aCoeffs;

  eta = params->eta;
  r2 = r*r;
  r3 = r2*r;
  r4 = r2*r2;
  r5 = r2*r3;

  u = 1./ r;

  u2 = u*u;
  u3 = u2*u;
  u4 = u2*u2;

  NA = r4 * aCoeffs->n4
     + r5 * aCoeffs->n5;

  DA = aCoeffs->d0
     + r  * aCoeffs->d1
     + r2 * aCoeffs->d2
     + r3 * aCoeffs->d3
     + r4 * aCoeffs->d4
     + r5 * aCoeffs->d5;

  dNA = 4. * aCoeffs->n4 * r3
      + 5. * aCoeffs->n5 * r4;

  dDA = aCoeffs->d1
      + 2. * aCoeffs->d2 * r
      + 3. * aCoeffs->d3 * r2
      + 4. * aCoeffs->d4 * r3
      + 5. * aCoeffs->d5 * r4;

  d2NA = 12. * aCoeffs->n4 * r2
       + 20. * aCoeffs->n5 * r3;

  d2DA = 2. * aCoeffs->d2
       + 6. * aCoeffs->d3 * r
       + 12. * aCoeffs->d4 * r2
       + 20. * aCoeffs->d5 * r3;

  A      = NA/DA;
  dAdr   = ( dNA * DA - dDA * NA ) / (DA*DA); 
  d2Adr2 = (d2NA*DA - d2DA*NA - 2.*dDA*DA*dAdr)/(DA*DA);

  /* The next derivatives are with respect to u, not r */

  dA  = - r2 * dAdr;
  d2A =  r4 * d2Adr2 + 2.*r3 * dAdr;

  v = cbrt(omega);
  FDIS = -params->in3copy.flux(v, params->in3copy.coeffs)/(eta*omega);

  twoUAPlusu2dA = 2.* u * A + u2 * dA;
  x1 = -r2 * sqrt (-dA * twoUAPlusu2dA * twoUAPlusu2dA * twoUAPlusu2dA )
                / (2.* u * dA * dA + A*dA - u * A * d2A);
  return (FDIS * x1);
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
   LALEOBPPWaveformEngine(status->statusPtr, signalvec, NULL, 
			NULL, &count, params, &paramsInit);
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
   LALEOBPPWaveformEngine(status->statusPtr, signalvec1, signalvec2, 
			   NULL, &count, params, &paramsInit);
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

  REAL4Vector *h=NULL;/* pointers to generated polarization data */

  REAL8 UNUSED phiC;/* phase at coalescence */
  InspiralInit paramsInit;

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
  XLALPrintWarning( "WARNING: Amp Order has been reset to %d\n", params->ampOrder);

  /* Compute some parameters*/
  LALInspiralInit(status->statusPtr, params, &paramsInit);
  CHECKSTATUSPTR(status);

  if (paramsInit.nbins==0)
    {
      DETATCHSTATUSPTR(status);
      RETURN (status);
    }
  /* Now we can allocate memory and vector for coherentGW structure*/
  if ( (h = XLALCreateREAL4Vector(2*paramsInit.nbins)) == NULL )
  {
    ABORTXLAL( status );
  }

  /* By default the waveform is empty */
  memset(h->data, 0, 2 * paramsInit.nbins * sizeof(REAL4));

  /* Call the engine function */
  params->startPhase = ppnParams->phi;
  LALEOBPPWaveformEngine(status->statusPtr, NULL, NULL, h, 
			   &count, params, &paramsInit);
  BEGINFAIL( status )
  {
     XLALDestroyREAL4Vector( h );
  }
  ENDFAIL( status );

  /* Check an empty waveform hasn't been returned */
  for (i = 0; i < h->length; i++)
  {
    if (h->data[i] != 0.0) break;
    if (i == h->length - 1)
    {
      XLALDestroyREAL4Vector(h);
      XLALPrintError( "An empty waveform has been generated!\n" );
      ABORT( status, LALINSPIRALH_ENOWAVEFORM, LALINSPIRALH_MSGENOWAVEFORM );
    }
  }

  /* Allocate the waveform structures. */
  if ( ( waveform->h = (REAL4TimeVectorSeries *)
         LALMalloc( sizeof(REAL4TimeVectorSeries) ) ) == NULL )
  {
    XLALDestroyREAL4Vector( h );
    ABORT( status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM );
  }
  memset( waveform->h, 0, sizeof(REAL4TimeVectorSeries) );

  if ( (waveform->h->data = XLALCreateREAL4VectorSequence( count, 2 )) == NULL )
  {
    LALFree( waveform->h );
    XLALDestroyREAL4Vector( h );
    ABORTXLAL( status );
  }


  memcpy(waveform->h->data->data , h->data, 2*count*(sizeof(REAL4)));

  waveform->h->deltaT = 1./params->tSampling;

  waveform->h->sampleUnits = lalStrainUnit;
  waveform->position = ppnParams->position;
  waveform->psi = ppnParams->psi;

  snprintf( waveform->h->name,
            LALNameLength, "EOB inspiral polarizations");

  /* --- fill some output ---*/
  /* Note that dfdt used to be set here. It has no relevance    */
  /* to this code, and was in fact probably wrong.              */
  /* However, certain codes use it to check for sampling issues */
  /* so if things fail, we may need to do something.            */
  ppnParams->tc     = (double)(count-1) / params->tSampling ;
  ppnParams->length = count;
  ppnParams->fStop  = params->fFinal;
  ppnParams->termCode        = GENERATEPPNINSPIRALH_EFSTOP;
  ppnParams->termDescription = GENERATEPPNINSPIRALH_MSGEFSTOP;

  ppnParams->fStart   = ppnParams->fStartIn;

  /* --- free memory --- */

  XLALDestroyREAL4Vector(h);

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
                UINT4            *countback,
                InspiralTemplate *params,
                InspiralInit     *paramsInit
                )
{


   UINT4                   count, nn=4, length = 0, hiSRndx=0;

   REAL4Vector             *sig1, *sig2, *freq;
   REAL8Vector             rVec, phiVec, prVec, pPhiVec;
   REAL8Vector             rVecHi, phiVecHi, prVecHi, pPhiVecHi;

   REAL8                   eta, m, r, s, p, q, dt, t, v, omega, f, ampl0;
   REAL8                   omegaOld;

   void                    *funcParams1, *funcParams2, *funcParams3;

   REAL8Vector             *values, *dvalues;
   InspiralDerivativesIn   in3;

   /* Variables for the integrator */
   ark4GSLIntegrator       *integrator = NULL;
   REAL8Array              *dynamics   = NULL;
   REAL8Array              *dynamicsHi = NULL;  
   INT4                    retLen;
   REAL8                   tMax;

   pr3In                   pr3in;
   expnCoeffs              ak;
   expnFunc                func;
   rOfOmegaIn              rofomegain;
   DFindRootIn             rootIn2, rootIn3;

   /* Stuff for pre-computed EOB values */
   EOBParams eobParams;
   EOBACoefficients        aCoeffs;
   FacWaveformCoeffs       hCoeffs;
   NewtonMultipolePrefixes prefixes;

   REAL8 (*rOfOmegaFunc)(REAL8, void *); /* Function to be used in root finding later */

   /* Variables to allow the waveform to be generated */
   /* from a specific fLower */
   REAL8                   sInit;  /* Initial phase, and phase to subtract */
   
   REAL8                   rmin = 20;        /* Smallest value of r at which to generate the waveform */
   COMPLEX16  MultSphHarmP;    /* Spin-weighted spherical harmonics */
   COMPLEX16  MultSphHarmM;    /* Spin-weighted spherical harmonics */
   COMPLEX16  hLM;             /* Factorized waveform */
   COMPLEX16  hNQC;            /* Non-quasicircular correction */
   REAL4      x1, x2;
   UINT4      i, j, k, modeL;
   INT4       modeM;          /* number of modes required */
   REAL4      inclination;    /* binary inclination       */
   REAL4      coa_phase;      /* binary coalescence phase */
   REAL8      y_1, y_2, z1, z2; /* (2,2) and (2,-2) spherical harmonics needed in (h+,hx) */


   /* Used for EOBNR */
   COMPLEX8Vector *modefreqs;
   UINT4 resampFac;
   UINT4 resampPwr; /* Power of 2 for resampling */
   REAL8 resampEstimate;

   /* For checking XLAL return codes */
   INT4 xlalStatus;

   /* Variables used in injection */
   REAL8 unitHz;
   REAL8 cosI;/* cosine of system inclination */
   REAL8 apFac, acFac;/* extra factor in plus and cross amplitudes */

   /* Accuracy of root finding algorithms */
   const REAL8 xacc = 1.0e-12;

   /* Accuracies of adaptive Runge-Kutta integrator */
   const REAL8 EPS_ABS = 1.0e-12;
   const REAL8 EPS_REL = 1.0e-10;

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

   /* Used in the calculation of the non-quasicircular correctlon */
   REAL8Vector             *ampNQC, *q1, *q2, *q3, *p1, *p2;
   EOBNonQCCoeffs           nqcCoeffs;

   /* Inidices of the ringdown matching points */
   /* peakIdx is the index where omega is a maximum */
   /* finalIdx is the index of the last point before the */
   /* integration breaks */
   UINT4Vector             *rdMatchPoint;
   UINT4                   peakIdx  = 0;
   UINT4                   finalIdx = 0;

   /* Counter used in unwrapping the waveform phase for NQC correction */
   INT4 phaseCounter;

   INITSTATUS(status, "LALEOBPPWaveformEngine", LALEOBPPWAVEFORMC);
   ATTATCHSTATUSPTR(status);


   ak   = paramsInit->ak;
   func = paramsInit->func;

   ASSERT(ak.totalmass/LAL_MTSUN_SI > 0.4, status,
   	LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   /* Check order is consistent */
   if ( params->order != LAL_PNORDER_PSEUDO_FOUR )
   {
     XLALPrintError( "Order must be LAL_PNORDER_PSEUDO_FOUR for approximant EOBNRv2.\n" );
     ABORT( status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE );
   }

   if (signalvec1) length = signalvec1->length; else if (h) length = h->length;

   /* Allocate some memory */
   values    = XLALCreateREAL8Vector( nn );
   dvalues   = XLALCreateREAL8Vector( nn );

   if ( !values || !dvalues )
   {
     XLALDestroyREAL8Vector( values );
     XLALDestroyREAL8Vector( dvalues );
     ABORTXLAL( status );
   }

   /* Set dt to sampling interval specified by user */
   dt =1./params->tSampling;
   eta = ak.eta;
   m = ak.totalmass;

   /* The maximum allowed time is determined by the length of the vectors */
   tMax = length * dt;

   /* only used in injection case */
   unitHz = m*(REAL8)LAL_PI;
   cosI   = cos( params->inclination );
   apFac  = -2.0 * (1.0 + cosI*cosI) * eta * params->totalMass * LAL_MRSUN_SI/params->distance;
   acFac  = -4.0 * cosI * eta * params->totalMass * LAL_MRSUN_SI/params->distance;

   /* Set the amplitude depending on whether the distance is given */
   if ( params->distance > 0.0 )
     ampl0  = params->totalMass * LAL_MRSUN_SI/params->distance;
   else
     ampl0  = params->signalAmplitude;

   /* Check we get a sensible answer */
   if ( ampl0 == 0.0 )
   {
     XLALPrintWarning( "Generating waveform of zero amplitude!!\n" );
   }

   /* Check that the 220 QNM freq. is less than the Nyquist freq. */
   /* Get QNM frequencies */
   modefreqs = XLALCreateCOMPLEX8Vector( 3 );
   xlalStatus = XLALGenerateQNMFreqV2( modefreqs, params, 2, 2, 3 );
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
     XLALPrintError( "Ringdown freq greater than Nyquist freq. "
           "Increase sample rate or consider using EOB approximant.\n" );
     ABORT( status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
   }

   /* Calculate the time we will need to step back for ringdown */
   tStepBack = 20.0 * params->totalMass * LAL_MTSUN_SI;
   nStepBack = ceil( tStepBack * params->tSampling );

   /* Set up structures for pre-computed EOB coefficients */
   memset( &eobParams, 0, sizeof(eobParams) );
   eobParams.eta = eta;
   eobParams.m1  = params->mass1;
   eobParams.m2  = params->mass2;
   eobParams.aCoeffs  = &aCoeffs;
   eobParams.hCoeffs  = &hCoeffs;
   eobParams.prefixes = &prefixes;

   if ( XLALCalculateEOBACoefficients( &aCoeffs, eta ) == XLAL_FAILURE )
   {
     ABORTXLAL( status );
   }

  if ( XLALCalcFacWaveformCoefficients( &hCoeffs, eta) == XLAL_FAILURE )
  {
    ABORTXLAL( status );
  }

  if ( XLALComputeNewtonMultipolePrefixes( &prefixes, eobParams.m1, eobParams.m2 )
         == XLAL_FAILURE )
  {
    ABORTXLAL( status );
  }

  funcParams3 = (void *) &eobParams;

   /* Allocate vectors to keep track of previous values */
   pPrev = XLALCreateREAL8Vector( nStepBack );
   qPrev = XLALCreateREAL8Vector( nStepBack );
   rPrev = XLALCreateREAL8Vector( nStepBack );
   sPrev = XLALCreateREAL8Vector( nStepBack );

   /* Calculate the resample factor for attaching the ringdown */
   /* We want it to be a power of 2 */
   /* Of course, we only want to do this if the required SR > current SR... */
   /* The form chosen for the resampleEstimate will essentially set */
   /* deltaT = M / 20. ( or less taking into account the power of 2 stuff */
   resampEstimate = 20. / ( params->totalMass * LAL_MTSUN_SI * params->tSampling );
   resampFac      = 1;

   if ( resampEstimate > 1 )
   {
     resampPwr = (UINT4)ceil(log2(resampEstimate));
     while ( resampPwr-- )
     {
       resampFac *= 2u;
     }
   }

   /* The length of the vectors at the higher sample rate will be */
   /* the step back time plus the ringdown */
   lengthHiSR = ( nStepBack + (UINT4)(20.0 / modefreqs->data[0].im / dt) ) * resampFac;

   /* Double it for good measure */
   lengthHiSR *= 2;

   /* We are now done with the ringdown modes - destroy them here */
   XLALDestroyCOMPLEX8Vector( modefreqs );
   modefreqs = NULL;

   /* Find the initial velocity given the lower frequency */
   f     = params->fLower;
   omega = f * LAL_PI * m;
   v     = cbrt( omega );

   /* Then the initial phase */
   s = params->startPhase;
   sInit = s;

   /* initial r as a function of omega - where to start evolution */
   rofomegain.eta = eta;
   rofomegain.omega = omega;
   rootIn3.xacc = xacc;
   rootIn2.xmax = 1000.;
   rootIn2.xmin = 3.;
   pr3in.eta = eta;
   pr3in.omegaS = params->OmegaS;
   pr3in.zeta2 = params->Zeta2;
   pr3in.aCoeffs = &aCoeffs;

   /* We will be changing the starting r if it is less than rmin */
   /* Therefore, we should reset pr3in.omega later if necessary. */
   /* For now we need it so that we can see what the initial r is. */

   pr3in.omega = omega;
   in3.totalmass = ak.totalmass;
   in3.dEnergy = func.dEnergy;
   in3.flux = func.flux;
   in3.coeffs = &ak;
   in3.nqcCoeffs = &nqcCoeffs;

   funcParams1 = (void *) &rofomegain;

   switch (params->order)
   {
     case LAL_PNORDER_PSEUDO_FOUR:
       rOfOmegaFunc = XLALrOfOmegaP4PN;
       funcParams2 = (void *) &pr3in;
       break;
     default:
       XLALPrintError( "There are no EOBNRv2 waveforms implemented at order %d\n", params->order);
       XLALDestroyREAL8Vector( values );
       XLALDestroyREAL8Vector( dvalues );
       ABORT( status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
   }
   r = XLALDBisectionFindRoot( rOfOmegaFunc, rootIn2.xmin,
              rootIn2.xmax, xacc, funcParams2);
   if ( XLAL_IS_REAL8_FAIL_NAN( r ) )
   {
     ABORTXLAL( status );
   }

   /* We want the waveform to generate from a point which won't cause */
   /* problems with the initial conditions. Therefore we force the code */
   /* to start at least at r = rmin (in units of M). */

   if (r < rmin) printf("Starting r will be reset.\n");
 
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
       omega = omegaofrP4PN( r, eta, &aCoeffs );
       pr3in.omega = omega;
       q = LALpphiInitP4PN(r, &aCoeffs );
       rootIn3.function = LALprInitP4PN;
       /* first we compute vr (we need coeef->Fp6) */
       pr3in.q = q;
       funcParams2 = (void *) &pr3in;
       pr3in.vr = LALvrP4PN(r, omega, (void *) &pr3in);
       /* then we compute the initial value of p */
       LALDBisectionFindRoot(status->statusPtr, &p, &rootIn3, funcParams2);
       CHECKSTATUSPTR(status);
       /* We need to change P to be the tortoise co-ordinate */
       /* TODO: Change prInit to calculate this directly */
       p = p * XLALCalculateEOBA(r, &aCoeffs);
       p = p / sqrt( XLALCalculateEOBD( r, eta ) );
       break;
     default:
       XLALPrintError( "There are no EOB/EOBNR waveforms implemented at order %d\n", params->order );
       XLALDestroyREAL8Vector( values );
       XLALDestroyREAL8Vector( dvalues );
       ABORT( status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
   }

   values->data[0] = r;
   values->data[1] = s;
   values->data[2] = p;
   values->data[3] = q;

   /* Allocate memory for temporary arrays */
   sig1 = XLALCreateREAL4Vector ( length );
   sig2 = XLALCreateREAL4Vector ( length );
   freq = XLALCreateREAL4Vector ( length );

   if ( !sig1 || !sig2 || !freq )
   {
     if ( sig1 ) XLALDestroyREAL4Vector( sig1 );
     if ( sig2 ) XLALDestroyREAL4Vector( sig2 );
     if ( freq ) XLALDestroyREAL4Vector( freq );
     XLALDestroyREAL8Vector( values );
     XLALDestroyREAL8Vector( dvalues );
     ABORT( status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM );
   }

   memset(sig1->data, 0, sig1->length * sizeof( REAL4 ));
   memset(sig2->data, 0, sig2->length * sizeof( REAL4 ));
   memset(freq->data, 0, freq->length * sizeof( REAL4 ));

   /* And their higher sample rate counterparts */
   /* Allocate memory for temporary arrays */
   sig1Hi = XLALCreateREAL4Vector ( lengthHiSR );
   sig2Hi = XLALCreateREAL4Vector ( lengthHiSR );
   amplHi = XLALCreateREAL4Vector ( lengthHiSR*2 );
   freqHi = XLALCreateREAL4Vector ( lengthHiSR );
   phseHi = XLALCreateREAL8Vector ( lengthHiSR );

   /* Allocate NQC vectors */
   ampNQC = XLALCreateREAL8Vector ( lengthHiSR );
   q1     = XLALCreateREAL8Vector ( lengthHiSR );
   q2     = XLALCreateREAL8Vector ( lengthHiSR );
   q3     = XLALCreateREAL8Vector ( lengthHiSR );
   p1     = XLALCreateREAL8Vector ( lengthHiSR );
   p2     = XLALCreateREAL8Vector ( lengthHiSR );

   if ( !sig1Hi || !sig2Hi || !amplHi || !freqHi || !phseHi )
   {
     if ( sig1 ) XLALDestroyREAL4Vector( sig1 );
     if ( sig2 ) XLALDestroyREAL4Vector( sig2 );
     if ( freq ) XLALDestroyREAL4Vector( freq );
     XLALDestroyREAL8Vector( values );
     XLALDestroyREAL8Vector( dvalues );
     ABORT( status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM );
   }

   memset(sig1Hi->data, 0, sig1Hi->length * sizeof( REAL4 ));
   memset(sig2Hi->data, 0, sig2Hi->length * sizeof( REAL4 ));
   memset(amplHi->data, 0, amplHi->length * sizeof( REAL4 ));
   memset(freqHi->data, 0, freqHi->length * sizeof( REAL4 ));
   memset(phseHi->data, 0, phseHi->length * sizeof( REAL8 ));
   memset(ampNQC->data, 0, ampNQC->length * sizeof( REAL8 ));
   memset(q1->data, 0, q1->length * sizeof( REAL8 ));
   memset(q2->data, 0, q2->length * sizeof( REAL8 ));
   memset(q3->data, 0, q3->length * sizeof( REAL8 ));
   memset(p1->data, 0, p1->length * sizeof( REAL8 ));
   memset(p2->data, 0, p2->length * sizeof( REAL8 ));

   /* Initialize the GSL integrator */
   if (!(integrator = XLALAdaptiveRungeKutta4Init(nn, LALHCapDerivativesP4PN, XLALFirstStoppingCondition, EPS_ABS, EPS_REL)))
   {
     XLALDestroyREAL4Vector( sig1 );
     XLALDestroyREAL4Vector( sig2 );
     XLALDestroyREAL4Vector( freq );
     XLALDestroyREAL8Vector( values );
     XLALDestroyREAL8Vector( dvalues );
     ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
   }

   count = 0;
   if (h || signalvec2)
      params->nStartPad = 0; /* must be zero for templates and injection */

   count = params->nStartPad;

  printf( "ICs = %e, %e, %e, %e\n", values->data[0], values->data[1], values->data[2], values->data[3] );

   /* Use the new adaptive integrator */
   /* TODO: Implement error checking */
   retLen = XLALAdaptiveRungeKutta4( integrator, &eobParams, values->data, 0., tMax/m, dt/m, &dynamics );

   /* We should have integrated to the peak of the frequency by now */
   hiSRndx = retLen - nStepBack;

   /* Set up the vectors, and re-initialize everything for the high sample rate */
   rVec.length  = phiVec.length = prVec.length = pPhiVec.length = retLen;
   rVec.data    = dynamics->data+retLen;
   phiVec.data  = dynamics->data+2*retLen;
   prVec.data   = dynamics->data+3*retLen;
   pPhiVec.data = dynamics->data+4*retLen;

   dt = dt/(REAL8)resampFac;
   values->data[0] = rVec.data[hiSRndx];
   values->data[1] = phiVec.data[hiSRndx];
   values->data[2] = prVec.data[hiSRndx];
   values->data[3] = pPhiVec.data[hiSRndx];

   /* We want to use a different stopping criterion for the higher sample rate */
   integrator->stop = XLALHighSRStoppingCondition;

   retLen = XLALAdaptiveRungeKutta4( integrator, &eobParams, values->data, 
     0, (lengthHiSR-1)*dt/m, dt/m, &dynamicsHi );

   rVecHi.length  = phiVecHi.length = prVecHi.length = pPhiVecHi.length = retLen;
   rVecHi.data    = dynamicsHi->data+retLen;
   phiVecHi.data  = dynamicsHi->data+2*retLen;
   prVecHi.data   = dynamicsHi->data+3*retLen;
   pPhiVecHi.data = dynamicsHi->data+4*retLen;

   /* We can now start calculating things for NQCs, and hiSR waveform */
   omegaOld = 0.0;
   phaseCounter = 0;

   for ( i=0; i < (UINT4)retLen; i++ )
   {
     omega = XLALCalculateOmega( eta, rVecHi.data[i], prVecHi.data[i], pPhiVecHi.data[i], &aCoeffs );
     /* For now we re-populate values - there may be a better way to do this */
     values->data[0] = r = rVecHi.data[i];
     values->data[1] = s = phiVecHi.data[i];
     values->data[2] = p = prVecHi.data[i];
     values->data[3] = q = pPhiVecHi.data[i];

     v = cbrt( omega ); 

     /* Only 2,2 mode for now */
     xlalStatus = XLALGetFactorizedWaveform( &hLM, values, v, 2, 2, &eobParams );

     ampNQC->data[i] = XLALCOMPLEX16Abs( hLM );
     sig1Hi->data[i] = (REAL4) ampl0 * hLM.re;
     sig2Hi->data[i] = (REAL4) ampl0 * hLM.im;
     phseHi->data[i] = XLALCOMPLEX16Arg( hLM ) + phaseCounter * LAL_TWOPI;
     if ( i && phseHi->data[i] > phseHi->data[i-1] )
     {
       phaseCounter--;
       phseHi->data[i] -= LAL_TWOPI;
     }
     q1->data[i] = p*p / (r*r*omega*omega);
     q2->data[i] = q1->data[i] / r;
     q3->data[i] = q2->data[i] / sqrt(r);
     p1->data[i] = p / ( r*omega );
     p2->data[i] = p1->data[i] * p*p;

     if ( omega <= omegaOld && !peakIdx )
     {
       peakIdx = i-1;
     }
     omegaOld = omega;
   }
   finalIdx = retLen - 1;

  /* Set the coalescence time */
  t = m * (dynamics->data[hiSRndx] + dynamicsHi->data[peakIdx]);

  /* Calculate the NQC correction */
  XLALCalculateNQCCoefficients( ampNQC, phseHi, q1,q2,q3,p1,p2, peakIdx, dt/m, eta, &nqcCoeffs );

  /* We can now calculate the waveform */
  i = 0;

  /* Find the point where we reach the low frequency cutoff */
  REAL8 lfCut = f * LAL_PI*m;

  while ( i < hiSRndx )
  {
    omega = XLALCalculateOmega( eta, rVec.data[i], prVec.data[i], pPhiVec.data[i], &aCoeffs );
    if ( omega > lfCut || fabs( omega - lfCut ) < 1.0e-5 )
    {
      break;
    }
    i++;
  }

  if ( i == hiSRndx )
  {
    XLALPrintError( "We don't seem to have crossed the low frequency cut-off\n" );
    ABORT( status, LALINSPIRALH_ENOWAVEFORM, LALINSPIRALH_MSGENOWAVEFORM );
  }

  /* Now create the (low-sampled) part of the waveform */
  FILE *out = fopen( "dynamics.dat", "w" );
  while ( i < hiSRndx )
  {
     omega = XLALCalculateOmega( eta, rVec.data[i], prVec.data[i], pPhiVec.data[i], &aCoeffs );
     /* For now we re-populate values - there may be a better way to do this */
     values->data[0] = r = rVec.data[i];
     values->data[1] = s = phiVec.data[i];
     values->data[2] = p = prVec.data[i];
     values->data[3] = q = pPhiVec.data[i];

    fprintf( out, "%.12e %.12e %.12e %.12e %.12e\n", dynamics->data[i], r, s, p, q );

     v = cbrt( omega );

     /* Only 2,2 mode for now */
     xlalStatus = XLALGetFactorizedWaveform( &hLM, values, v, 2, 2, &eobParams );

     xlalStatus = XLALEOBNonQCCorrection( &hNQC, values, omega, &nqcCoeffs );
     hLM = XLALCOMPLEX16Mul( hNQC, hLM );

     sig1->data[i] = (REAL4) ampl0 * hLM.re;
     sig2->data[i] = (REAL4) ampl0 * hLM.im;

     i++;
  }
  fclose( out ); out = NULL;

  /* Now apply the NQC correction to the high sample part */
  out = fopen( "dynamics_hi.dat", "w" );
  FILE *out2 = fopen( "waveform_hi.dat", "w" );
  FILE *out3 = fopen( "waveform_hi_nonqc.dat", "w" );
  for ( i = 0; i <= finalIdx; i++ )
  {
    omega = XLALCalculateOmega( eta, rVecHi.data[i], prVecHi.data[i], pPhiVecHi.data[i], &aCoeffs );

    /* For now we re-populate values - there may be a better way to do this */
    values->data[0] = r = rVecHi.data[i];
    values->data[1] = s = phiVecHi.data[i];
    values->data[2] = p = prVecHi.data[i];
    values->data[3] = q = pPhiVecHi.data[i];

    fprintf( out, "%.12e %.12e %.12e %.12e %.12e\n", dynamics->data[hiSRndx] + dynamicsHi->data[i], r, s, p, q );

    xlalStatus = XLALEOBNonQCCorrection( &hNQC, values, omega, &nqcCoeffs );

    hLM.re = sig1Hi->data[i];
    hLM.im = sig2Hi->data[i];

    fprintf( out3, "%.12e %.12e %.12e\n", dynamics->data[hiSRndx] + dynamicsHi->data[i], hLM.re, hLM.im );

    hLM = XLALCOMPLEX16Mul( hNQC, hLM );
    sig1Hi->data[i] = (REAL4) hLM.re;
    sig2Hi->data[i] = (REAL4) hLM.im;

    fprintf( out2, "%.12e %.12e %.12e\n", dynamics->data[hiSRndx] + dynamicsHi->data[i], hLM.re, hLM.im );

  }
  fclose( out ); out = NULL;
  fclose( out2 );
   /*----------------------------------------------------------------------*/
   /* Record the final cutoff frequency of BD Waveforms for record keeping */
   /* ---------------------------------------------------------------------*/
   params->vFinal = v;
   if (signalvec1 && !signalvec2) params->tC = t;

   /* For now we just set fFinal to Nyquist */
   /* This is a hack to get the filtering code to work properly */
   params->fFinal = params->tSampling/2.;

   XLALAdaptiveRungeKutta4Free( integrator );
   XLALDestroyREAL8Vector( values );
   XLALDestroyREAL8Vector( dvalues );

   /*--------------------------------------------------------------
    * Attach the ringdown waveform to the end of inspiral
     -------------------------------------------------------------*/
   REAL8 tmpSamplingRate = params->tSampling;
   params->tSampling *= resampFac;

   rdMatchPoint = XLALCreateUINT4Vector( 3 );

   /* Check the first matching point is sensible */
   if ( ceil( tStepBack * params->tSampling / 2.0 ) > peakIdx )
   {
     XLALPrintError( "Invalid index for first ringdown matching point.\n" );
     ABORT( status, LALINSPIRALH_ESIZE , LALINSPIRALH_MSGESIZE );
   }

   /*rdMatchPoint->data[0] = peakIdx - ceil( tStepBack * params->tSampling / 2.0 );*/
   /* Experiment */
   rdMatchPoint->data[0] = ceil( 5. * params->totalMass * LAL_MTSUN_SI * params->tSampling ) < peakIdx ?
       peakIdx - ceil( 5. * params->totalMass * LAL_MTSUN_SI * params->tSampling ) : 0;
   rdMatchPoint->data[1] = peakIdx;
   rdMatchPoint->data[2] = finalIdx;

   printf(" ringdown matching time = %e, next time  = %e\n", dynamics->data[hiSRndx] + dynamicsHi->data[peakIdx], dynamics->data[hiSRndx] + dynamicsHi->data[peakIdx+1] );

   XLALPrintInfo("Matching points: %u %u %u; last idx in array = %d\n", rdMatchPoint->data[0], rdMatchPoint->data[1], rdMatchPoint->data[2], sig1Hi->length);
   XLALPrintInfo("Time after peak that integration blows up: %e M\n", (rdMatchPoint->data[2] - rdMatchPoint->data[1]) * dt / (params->totalMass * LAL_MTSUN_SI ));
   
   XLALPrintInfo( "Ringdown matching points: %e, %e\n", 
           (REAL8)rdMatchPoint->data[0]/resampFac + (REAL8)hiSRndx, 
           (REAL8)rdMatchPoint->data[1]/resampFac + (REAL8)hiSRndx );

   xlalStatus = XLALInspiralHybridAttachRingdownWave(sig1Hi, sig2Hi,
                   rdMatchPoint, params);
   if (xlalStatus != XLAL_SUCCESS )
   {
     XLALDestroyREAL4Vector( sig1 );
     XLALDestroyREAL4Vector( sig2 );
     XLALDestroyREAL4Vector( freq );
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
     XLALDestroyREAL4Vector( freq );
     ABORTXLAL( status );
   }

   modeM = -2;
   xlalStatus = XLALSphHarm( &MultSphHarmM, modeL, modeM, inclination, coa_phase );
   if (xlalStatus != XLAL_SUCCESS )
   {
     XLALDestroyREAL4Vector( sig1 );
     XLALDestroyREAL4Vector( sig2 );
     XLALDestroyREAL4Vector( freq );
     ABORTXLAL( status );
   }

   y_1 =   MultSphHarmP.re + MultSphHarmM.re;
   y_2 =   MultSphHarmM.im - MultSphHarmP.im;
   z1 = - MultSphHarmM.im - MultSphHarmP.im;
   z2 =   MultSphHarmM.re - MultSphHarmP.re;

   /* Next, compute h+ and hx from h22, h22*, Y22, Y2-2 */
   out = fopen( "realAndImag.dat", "w" );
   for ( i = 0; i < sig1->length; i++)
   {
     freq->data[i] /= unitHz;
     x1 = sig1->data[i];
     x2 = sig2->data[i];
     fprintf( out, "%.12e %.12e %.12e\n", (REAL8)i/(m*params->tSampling), x1, x2);
     sig1->data[i] = (x1 * y_1) + (x2 * y_2);
     sig2->data[i] = (x1 * z1) + (x2 * z2);     
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

   /* Clean up */
   XLALDestroyREAL8Array( dynamics );
   XLALDestroyREAL8Array( dynamicsHi );
   XLALDestroyREAL4Vector ( sig1 );
   XLALDestroyREAL4Vector ( sig2 );
   XLALDestroyREAL4Vector ( freq );

   DETATCHSTATUSPTR(status);
   RETURN(status);
}
