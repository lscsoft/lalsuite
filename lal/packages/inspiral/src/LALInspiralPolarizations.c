/*
*  Copyright (C) 2008 Evan Ochsner
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

/*  <lalVerbatim file="LALInspiralPolarizationsCV">

Author: Ochsner, E.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralPolarizations.c} and \texttt{LALInspiralPolarizationsTemplates.c}}
This code generates the two GW polarizations of the inspiral waveform for a given phase model.


\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralPolarizationsCP}
\index{\verb&LALInspiralHplus()&}
\begin{itemize}
\item {\tt output:} Outputs either $h_+$ or $h_\times.$
\item {\tt params:} Input containing binary chirp parameters.
\item {\tt phase:} Input containing binary chirp phase.
\item {\tt v:} PN expansion parameter
\end{itemize}

\subsubsection*{Description}

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralPolarizationsCV}}

</lalLaTeX>  */

#include <lal/LALStdlib.h>
#include<lal/LALConstants.h>
#include<lal/LALInspiral.h>
/*-------------------------------------*/
/* Need to add a field ampOrder to the */
/*       InspiralTemplate struct       */
/*-------------------------------------*/
REAL4 LALInspiralHPlusPolarization( REAL8 phase, REAL8 v, InspiralTemplate *params )
{
	REAL8 hPlus = 0.0;
	REAL8 M = params->totalMass;
	REAL8 eta = params->eta;
	REAL8 eta2 = eta*eta;
	REAL8 dM = abs(params->mass1 - params->mass2) / params->totalMass;
	REAL8 cI = cos(params->inclination);
	REAL8 sI = sin(params->inclination);
	REAL8 sI2 = sI*sI;
	REAL8 cI2 = cI*cI;
	REAL8 cI4 = cI2*cI2;
	REAL8 cI6 = cI4*cI2;
	switch( params->ampOrder )
	{
		case LAL_PNORDER_TWO_POINT_FIVE:
		hPlus = hPlus + v * v * v * v * v * ( sI*dM*
                        ( (1771./5120. - cI2*1667./5120. + cI4*217./9216. - cI6/9216.)
                        + eta*(681./256. + cI2*13./768. - cI4*35./768. + cI6/2304.)
                        + eta2*(-3451./9216. + cI2*673./3072. - cI4*5./9216. - cI6/3072.) )
                        *cos(phase) + LAL_PI*( (19./3. + 3.*cI2 - cI4*2./3.)
                        + eta*(-16./3. + cI2*14./3. + 2.*cI4) )*cos(2.*phase)
                        + sI*dM* ( (3537./1024. - cI2*22977./5120. - cI4*15309./5120. + cI6*729./5120.)
                        + eta*(-23829./1280. + cI2*5529./1280. + cI4*7749./1280. - cI6*729./1280.)
                        + eta2*(29127./5120. - cI2*27267./5120. - cI4*1647./5120. + cI6*2187./5120.) )
                        *cos(3.*phase) - (16./3.)*LAL_PI*(1.+cI2)*sI2*(1.-3.*eta)*cos(4.*phase)
                        + sI*dM*( (-108125./9216. + cI2*40625./9216. + cI4*83125./9216. - cI6*15625./9216.)
                        + eta*(8125./256. - cI2*40625./2304. - cI4*48125./2304. + cI6*15625./2304.)
                        + eta2*(-119375./9216. + cI2*40625./3072. + cI4*44375./9216. - cI6*15625./3072.) )
                        *cos(5.*phase) + dM*sI2*sI2*sI*(1.+cI2)*(117649./46080.)*(1. - 4.*eta + 3.*eta2)
                        *cos(7.*phase) + ( (-9./5. + cI2*14./5. + cI4*7./5.)
                        + eta*(32. + cI2*56./5. - cI4*28./5.) )*sin(2.*phase) + sI2*(1.+cI2)
                        *( (56./5. - log(2.)*32./3.) - eta*(1193./30. - 32.*log(2.)) ) *sin(4.*phase)  );
		case LAL_PNORDER_TWO:
		hPlus = hPlus + v * v * v * v * (  (1./120.)*( (22.+396.*cI2 + 145.*cI4 - 5.*cI6)
                        + (5./3.)*eta*(706. - 216.*cI2 - 251.*cI4 + 15.*cI6)
                        - 5.*eta2* (98. - 108.*cI2 + 7.*cI4 + 5.*cI6) )* cos(2.*phase)
                        + (2./15.)*sI2* ( (59. + 35.*cI2 - 8.*cI4) - (5./3.)*eta*(131.+59.*cI2 - 24.*cI4)
                        + 5.*eta2* (21. - 3.*cI2 - 8.*cI4) )*cos(4.*phase)
                        - (81./40.)*(1. - 5.*eta + 5.*eta2)*sI2*sI2*(1.+cI2)* cos(6.*phase)
                        + (sI/40.)*dM* ( (11.+7.*cI2 + 10.*(5. + cI2)*log(2))*sin(phase)
                        - 5.*LAL_PI*(5. + cI2)*cos(phase) - 27.*(7. - 10.*log(3./2.))*(1 + cI2)*sin(3.*phase)
                        + 135.*LAL_PI*(1.+cI2)*cos(3.*phase) )  );
		case LAL_PNORDER_ONE_POINT_FIVE:
		hPlus = hPlus + v * v * v * (  (sI/192.)* dM*
                        (  ( (57.+60.*cI2 - cI4) - 2.*eta*(49. - 12.*cI2 - cI4) )* cos(phase)
                        - (27./2.)*( (73. + 40.*cI2 - 9.*cI4) - 2.*eta*(25. - 8.*cI2 - 9.*cI4) )
                        *cos(3.*phase) + (625./2.)*(1. - 2.*eta)*sI2*(1+cI2)* cos(5.*phase)  )
                        - 2.*LAL_PI* (1.+cI2)*cos(2.*phase)  );
		case LAL_PNORDER_ONE:
		hPlus = hPlus + v * v *(   (1./6.)*( (19.+9.*cI2 - 2.*cI4) - eta*(19. - 11.*cI2 - 6.*cI4) )
                        * cos( 2.*phase) - (4./3.)*sI2*(1.+cI2)*(1. - 3*eta)* cos(4.*phase)  );
		case LAL_PNORDER_HALF:
		hPlus = hPlus - v * (sI * dM / 8.) * ( (5. + cI2) * cos(phase) - 9.*(1. + cI2)*cos(3.*phase));
		case LAL_PNORDER_NEWTONIAN:
		hPlus = hPlus - (1.+cI2) * cos(2.*phase);
		break;
		default: fprintf(stderr, "There are no EOB waveforms at order %d in amplitude\n", params->ampOrder);
	}
        return 2.* M * LAL_MTSUN_SI* eta * v * v * hPlus * LAL_C_SI / ( LAL_PC_SI * 1.e6 );
}

REAL4 LALInspiralHCrossPolarization( REAL8 phase, REAL8 v, InspiralTemplate *params )
{
	REAL8 hCross = 0.0;
	REAL8 M = params->totalMass;
	REAL8 eta = params->eta;
	REAL8 eta2 = eta*eta;
	REAL8 dM = abs(params->mass1 - params->mass2) / params->totalMass;
	REAL8 cI = cos(params->inclination);
	REAL8 sI = sin(params->inclination);
	REAL8 sI2 = sI*sI;
	REAL8 cI2 = cI*cI;
	REAL8 cI4 = cI2*cI2;

	switch( params->ampOrder )
	{
		case LAL_PNORDER_TWO_POINT_FIVE:
		hCross = hCross + v * v * v * v * v * ( cI*( (2. - cI2*22./5.) + eta*(-282./5. + cI2*94./5.) )
                         *cos(2.*phase) + cI*sI2*( (-112./5. + log(2.)*64./3.) + eta*(1193./15. - 64.*log(2.)) )
                         *cos(4.*phase) + sI*cI*dM*( (-913./7680. + cI2*1891./11520. - cI4*7./4608.)
                         + eta*(1165./384. - cI2*235./576. + cI4*7./1152.)
                         + eta2*(-1301./4608. + cI2*301./2304. - cI4*7./1536.) )*sin(phase)
                         + cI*LAL_PI*( (34./3. - cI2*8./3.) - eta*(20./3. - 8.*cI2) )*sin(2.*phase)
                         + sI*cI*dM*( (12501./2560. - cI2*12069./1280. + cI4*1701./2560.)
                         + eta*(-19581./640. + cI2*7821./320. - cI4*1701./640.)
                         + eta2*(18903./2560. - cI2*11403./1280. + cI4*5103./2560.) )*sin(3.*phase)
                         - sI2*cI*LAL_PI*(32./3.)*(1. - 3.*eta)*sin(4.*phase)
                         + sI*cI*dM*( (-101875./4608. + cI2*6875./256. - cI4*21875./4608.)
                         + eta*(66875./1152. - cI2*44375./576. + cI4*21875./1152.)
                         + eta2*(-100625./4608. + cI2*83125./2304. - cI4*21875./1536.) )*sin(5.*phase)
                         + sI2*sI2*sI*cI*dM*(117649./23040.)*(1. - 4.*eta + 3.*eta2)*sin(7.*phase)  );
		case LAL_PNORDER_TWO:
		  hCross = hCross + v * v * v * v * (  (cI/60.)*( (68. + 226.*cI2 - 15.*cI4)
                           + (5./3.)*eta*(572. - 490.*cI2 + 45.*cI4)
                           - 5.*eta2*(56. - 70.*cI2 + 15.*cI4) )*sin(2.*phase)
                           + (4./15.)*cI*sI2*( (55. - 12.*cI2) - (5./3.)*eta*(119. - 36.*cI2)
                           + 5.*eta2*(17. - 12.*cI2) )*sin(4.*phase)
                           - (81./20.)*(1. - 5.*eta + 5.*eta2)*cI*sI2*sI2* sin(6.*phase)
                           - (3./20.)*sI*cI*dM* ( (3.+10.*log(2))*cos(phase) + 5.*LAL_PI*sin(phase)
                           - 9.*(7. - 10.*log(3./2.))*cos(3.*phase) - 45.*LAL_PI*sin(3.*phase) )  );
		case LAL_PNORDER_ONE_POINT_FIVE:
		  hCross = hCross + v * v * v * (  (sI*cI/96.)*dM* ( ( (63. - 5.*cI2) - 2.*eta*(23. - 5.*cI2)  )
                           *sin(phase) - (27./2.)*(  (67. - 15.*cI2) - 2.*eta*(19. - 15.*cI2)  )*sin(3.*phase)
                           + (625./2.)*(1. - 2.*eta)*sI2*sin(5.*phase)  ) - 4.*LAL_PI*cI* sin(2.*phase) ) ;
		case LAL_PNORDER_ONE:
		  hCross = hCross + v * v * (  (cI/3.)*( (17. - 4.*cI2) - eta* (13. - 12.*cI2) )
                           *sin(2.*phase) - (8./3.)*(1. - 3.*eta)*cI*sI2*sin(4.*phase)  );
		case LAL_PNORDER_HALF:
		hCross = hCross - v * (3./4.)*sI*cI*dM* ( sin(phase) - 3.* sin(3.*phase) );
		case LAL_PNORDER_NEWTONIAN:
		hCross = hCross - 2.* cI * sin(2.*phase);
		break;
		default:  fprintf(stderr, "There are no EOB waveforms at order %d in amplitude\n", params->ampOrder);
	}
        return 2.* M * LAL_MTSUN_SI * eta * v * v * hCross * LAL_C_SI/ ( LAL_PC_SI * 1.e6 );

}
