/*
*  Copyright (C) 2007 Jolien Creighton, Teviet Creighton
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

/**
\author Creighton, T. D.
\file
\ingroup pulsarTODO

\heading{Header \ref GenerateTaylorCW.h}
\latexonly\label{s_GenerateTaylorCW_h}\endlatexonly

Provides routines to generate Taylor-parameterized continuous
waveforms.

\heading{Synopsis}
\code
#include <lal/GenerateTaylorCW.h>
\endcode

This header covers routines to generate continuous quasiperiodic
waveforms whose frequency varies slowly and smoothly with time.  For
such sources the frequency function is normally described by its
Taylor "spindown" (or spin-up) coefficients.  This type of waveform
may be typical of objects such as neutron stars that are gradually
shedding angular momentum, or are accelerating in the gravitational
potential of a star cluster.  The Taylor expansion is likely
\e not suitable for long-term modelling of the frequency of waves
from glitching neutron stars, neutron stars in close binary orbits, or
neutron stars that are accreting or shedding angular momentum in a
stochastic manner.

The frequency and phase of such slowly-varying quasiperiodic sources
are given by their Taylor series:
\anchor eq_taylorcw-freq \anchor eq_taylorcw-phi \f{eqnarray}{
\label{eq_taylorcw-freq}
f(t)    & = & f_0 \left[ 1 + \sum_{k=1}^n f_k(t-t_0)^k \right] \;, \\
\label{eq_taylorcw-phi}
\phi(t) & = & \phi_0 + 2\pi f_0 \left[ (t-t_0) +
		\sum_{k=1}^n \frac{f_k}{k+1}(t-t_0)^{k+1} \right] \;,
\f}
where \f$f_k\f$ are the spin-normalized Taylor coefficients.  If the
source's spin is varying over some timescale \f$\tau\f$, one typically
expects that \f$f_k\sim\tau^{-k}\f$.  Note that in this and later
discussions, \f$f\f$ and \f$\phi\f$ refer to the frequency and phase of the
gravitational wave, which are typically some constant multiple of
(often twice) the frequency and phase of the rotating source.

The \c CoherentGW structure allows for a very general
description of waveforms with modulations in the amplitudes or
relative phases of the wave polarizations, as described in
\ref SimulateCoherentGW.h.  However, in this simplest model of
quasiperiodic waveforms, we neglect such phenomena as precession that
would produce these effects.  Thus for any given source one can choose
a polarization basis (described by some polarization angle \f$\psi\f$) in
which the wave has a constant elliptical polarization of the form:
\anchor eq_taylorcw-hplus \anchor eq_taylorcw-hcross \f{eqnarray}{
\label{eq_taylorcw-hplus}
h_+(t)      & = & A_+      \cos\phi(t) \;, \\
\label{eq_taylorcw-hcross}
h_\times(t) & = & A_\times \sin\phi(t) \;.
\f}



\heading{Types}

\heading{Structure \c TaylorCWParamStruc}


This structure stores the parameters for constructing a gravitational
waveform with a Taylor-polynomial frequency and phase.  As with the
\c PPNParamStruc type in \ref GeneratePPNInspiral.h, we divide
the fields into passed fields (which are supplied to the final
\c CoherentGW structure but not used in any calculations), input
fields (that are used by the waveform generator), and output fields
(that are set by the waveform generator).  They are:

\bigskip<em>Passed fields:</em>
<dl>
<dt><tt>SkyPosition position</tt></dt><dd> The location of the source on the
sky, normally in equatorial coordinates.</dd>

<dt><tt>REAL4 psi</tt></dt><dd> The polarization angle of the source, in
radians.</dd>

<dt><tt>LIGOTimeGPS epoch</tt></dt><dd> The start time \f$t_0\f$ of the output
series.</dd>
</dl>

<em>Input fields:</em>
<dl>
<dt><tt>REAL8 deltaT</tt></dt><dd> The requested sampling interval of the
waveform, in s.</dd>

<dt><tt>UINT4 length</tt></dt><dd> The number of samples in the generated
waveform.</dd>

<dt><tt>REAL4 aPlus, aCross</tt></dt><dd> The polarization amplitudes \f$A_+\f$,
\f$A_\times\f$, in dimensionless strain units.</dd>

<dt><tt>REAL8 phi0</tt></dt><dd> The wave phase at time \f$t_0\f$, in radians.</dd>

<dt><tt>REAL8 f0</tt></dt><dd> The wave frequency at time \f$t_0\f$, in Hz.</dd>

<dt><tt>REAL8Vector *f</tt></dt><dd> The spin-normalized Taylor parameters
\f$f_k\f$, as defined in Eq.\eqref{eq_taylorcw-freq}, above.  If
\c f=\c NULL, a monochromatic wave is generated.</dd>
</dl>

<em>Output fields:</em>
<dl>
<dt><tt>REAL4 dfdt</tt></dt><dd> The maximum value of \f$\Delta f\Delta t\f$
encountered over any timestep \f$\Delta t\f$ used in generating the
waveform.</dd>
</dl>

*/

#ifndef _GENERATETAYLORCW_H
#define _GENERATETAYLORCW_H

#include <lal/LALStdlib.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/SkyCoordinates.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

NRCSID( GENERATETAYLORCWH, "$Id$" );

/**
 \name Error Codes */ /*@{*/
#define GENERATETAYLORCWH_ENUL 1
#define GENERATETAYLORCWH_EOUT 2
#define GENERATETAYLORCWH_EMEM 3

#define GENERATETAYLORCWH_MSGENUL "Unexpected null pointer in arguments"
#define GENERATETAYLORCWH_MSGEOUT "Output field a, f, phi, or shift already exists"
#define GENERATETAYLORCWH_MSGEMEM "Out of memory"
/*@}*/


typedef struct tagTaylorCWParamStruc {
  /* Passed parameters. */
  SkyPosition position; /* location of source on sky */
  REAL4 psi;            /* polarization angle (radians) */
  LIGOTimeGPS epoch;    /* start time of output time series */

  /* Input parameters. */
  REAL8 deltaT;         /* requested sampling interval (s) */
  UINT4 length;         /* length of time series */
  REAL4 aPlus, aCross;  /* polarization amplitudes */
  REAL8 phi0;           /* initial phase */
  REAL8 f0;             /* initial frequency */
  REAL8Vector *f;       /* f0-normalized Taylor parameters */

  /* Output parameters. */
  REAL4 dfdt;           /* maximum value of df*dt over any timestep */
} TaylorCWParamStruc;


/* Function prototypes. */

void
LALGenerateTaylorCW( LALStatus          *,
		     CoherentGW         *output,
		     TaylorCWParamStruc *params );


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _GENERATETAYLORCW_H */
