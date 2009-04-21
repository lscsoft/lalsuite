/*
 * Copyright (C) 2004, 2005 Reinhard Prix
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
 * \author Reinhard Prix
 * \date 2005
 * \file 
 * \brief Routines to simulate pulsar-signals "exactly".

The motivation for this module is to provide functions to
simulate pulsar signals <em>with the best possible accuracy</em>, 
i.e. using no approximations, contrary to LALGeneratePulsarSignal(). 

Obviously this is not meant as a fast code to be used in a Monte-Carlo
simulation, but rather as a <em>reference</em> to compare other (faster)
functions agains, in order to be able to gauge the quality of a given
signal-generation routine.

We want to calculate \f$h(t)\f$, given by
\f[
\label{eq:1}
	h(t) = F_+(t)\, h_+(t) + F_\times(t) \,h_\times(t)\,,
\f]
where \f$F_+\f$ and \f$F_x\f$ are called the <em>beam-pattern</em> functions, 
which depend of the wave polarization \f$\psi\f$,
the source position \f$\alpha\f$, \f$\delta\f$ and the detector position and
orientation (\f$\gamma\f$, \f$\lambda\f$, \f$L\f$ and \f$\xi\f$). The expressions for 
the beam-pattern functions are given in \ref JKS98 "[JKS98]", which we write as
\f{eqnarray}
F_+(t) = \sin \zeta \cos 2\psi \, a(t)  + \sin \zeta \sin 2\psi \, b(t)\,,\\
F_\times(t) = \sin\zeta  \cos 2\psi \,b(t) - \sin\zeta \sin 2\psi \, a(t) \,,
\f}
where \f$\zeta\f$ is the angle between the interferometer arms, and 
\f{eqnarray}
a(t) &=& a_1 \cos[ 2 (\alpha - T)) ] + a_2 \sin[ 2(\alpha - T)]
+ a_3 \cos[ \alpha - T ] + a_4 \sin [ \alpha - T ] + a_5\,,\\
b(t) &=& b_1 \cos[ 2(\alpha - T)] + b_2 \sin[ 2(\alpha - T) ]
+ b_3 \cos[ \alpha - T ] + b_4 \sin[ \alpha - T]\,,
\f}
where \f$T\f$ is the local (mean) sidereal time of the detector, and the 
time-independent coefficients \f$a_i\f$ and \f$b_i\f$ are given by
\f{eqnarray}
a_1 &=& {1\over 16} \sin 2\gamma \,(3- \cos 2\lambda)\,(3 - \cos 2\delta)\,,\\
a_2 &=& -{1\over 4}\cos 2\gamma \,\sin \lambda \,(3 - \cos 2\delta) \,,\\
a_3 &=& {1\over 4} \sin 2\gamma \,\sin 2\lambda \,\sin 2\delta  \,\\
a_4 &=& -{1\over2} \cos 2\gamma \,\cos \lambda \,\sin 2 \delta\,,\\
a_5 &=& {3\over4} \sin 2\gamma \, \cos^2 \lambda \,\cos^2 \delta\,,
\f}
and 
\f{eqnarray}
b_1 &=& \cos 2\gamma \,\sin \lambda \,\sin \delta\,,\\
b_2 &=& {1\over 4} \sin 2\gamma \,(3-\cos 2\lambda)\, \sin \delta\,,\\
b_3 &=& \cos 2\gamma \,\cos \lambda \,\cos\delta \,, \\
b_4 &=& {1\over 2} \sin2\gamma \,\sin 2\lambda \,\cos\delta\,,
\f}

The source model considered is a plane-wave
\f{eqnarray}
h_+(t) &=& A_+\, \cos \Psi(t)\,,\\
h_\times(t) &=& A_\times \, \sin \Psi(t)\,,
\f}
where the wave-phase is \f$\Psi(t) = \Phi_0 + \Phi(t)\f$, and for an
isolated pulsar we have
\f{equation}
\Phi(t) = 2\pi \left[\sum_{s=0} {f^{(s)}(\tau_\mathrm{ref}) \over
(s+1)!} \left( \tau(t) - \tau_\mathrm{ref} \right)^{s+1} \right]\,,
\f}
where \f$\tau_\mathrm{ref}\f$ is the "reference time" for the definition
of the pulsar-parameters \f$f^{(s)}\f$ in the solar-system barycenter
(SSB), and \f$\tau(t)\f$ is the SSB-time of the phase arriving at the
detector at UTC-time \f$t\f$, which depends on the source-position
(\f$\alpha\f$, \f$\delta\f$) and the detector-position, namely
\f{equation}
  \tau (t) = t + { \vec{r}(t)\cdot\vec{n} \over c}\,,
\f}
where \f$\vec{r}(t)\f$ is the vector from SSB to the detector, and \f$\vec{n}\f$ 
is the unit-vector pointing \emph{to} the source.

This is a standalone "clean-room" implementation using no other 
outside-functions <em>except</em> for LALGPStoLMST1() to calculate 
the local (mean) sidereal time at the detector for given GPS-time, 
(which I double-checked with an independent Mathematica script),
and and LALBarycenter() to calculate \f$\tau(t)\f$.

 * $Id$
 *
 */

#include <lal/AVFactories.h>
#include <lal/TimeSeries.h>
#include <lal/GeneratePulsarSignal.h>
#include <lal/ComputeFstat.h>
#include <lal/ExtrapolatePulsarSpins.h>

NRCSID( SIMULATEPULSARSIGNALC, "$Id$");

extern INT4 lalDebugLevel;


/*----- Macros ----- */

#define SQ(x) ((x) * (x))

/** convert GPS-time to REAL8 */
#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

/** copy 3 components of Euklidean vector */
#define COPY_VECT(dst,src) do { (dst)[0] = (src)[0]; (dst)[1] = (src)[1]; (dst)[2] = (src)[2]; } while(0)

/** Simple Euklidean scalar product for two 3-dim vectors in cartesian coords */
#define SCALAR(u,v) ((u)[0]*(v)[0] + (u)[1]*(v)[1] + (u)[2]*(v)[2])

#define INIT_MEM(x) memset(&(x), 0, sizeof((x)))

#define TRUE  (1==1)
#define FALSE (1==0)

#define oneBillion 1000000000L

#define NUM_SPINDOWNS 	3

/* error-codes */
#define SIMULATEPULSARSIGNAL_ENULL 		1
#define SIMULATEPULSARSIGNAL_ENONULL		2
#define SIMULATEPULSARSIGNAL_EMEM		3
#define SIMULATEPULSARSIGNAL_ESYS		4
#define SIMULATEPULSARSIGNAL_EINPUT		5
#define SIMULATEPULSARSIGNAL_EFUN		6


#define SIMULATEPULSARSIGNAL_MSGENULL 		"Arguments contained an unexpected null pointer"
#define SIMULATEPULSARSIGNAL_MSGENONULL		"Output pointer is not NULL"
#define SIMULATEPULSARSIGNAL_MSGEMEM		"Out of memory"
#define SIMULATEPULSARSIGNAL_MSGESYS		"System error, probably while File I/O"
#define SIMULATEPULSARSIGNAL_MSGEINPUT		"Invalid input-arguments to function"
#define SIMULATEPULSARSIGNAL_MSGEFUN		"Subroutine failed"

static LALUnit emptyUnit;



/** Simulate a pulsar signal to best accuracy possible.
 *
 * NOTE: currently only isolated pulsars are supported
 *
 * NOTE2: we don't really use the highest possible accuracy right now,
 *   as we blatently neglect all relativistic timing effects (i.e. using dT=v.n/c)
 *
 * NOTE3: no heterodyning is performed here, the time-series is generated and sampled
 * at the given rate, that's all! ==> the caller needs to make sure about the 
 * right sampling rate to use (->aliasing) and do the proper post-treatment...
 *
 */
void
LALSimulateExactPulsarSignal (LALStatus *status, 
			      REAL4TimeSeries **timeSeries, 
			      const PulsarSignalParams *params)
{
  LALFrDetector *site = &(params->site->frDetector);
  REAL8 Delta, Alpha;
  UINT4 numSteps, i;

  REAL8 refTime, startTimeSSB;
  DetectorStateSeries *detStates = NULL;
  LIGOTimeGPSVector *timestamps = NULL;
  REAL8 dt;
  REAL8 vn[3];
  REAL8 A1, A2, A3, A4;
  REAL8 phi0, f0, f1dot, f2dot, f3dot;
  AMCoeffs *amcoe;
  REAL8 xAzi, yAzi;
  REAL8 Zeta, sinZeta;
  UINT4 numSpins = PULSAR_MAX_SPINS;

  CHAR *channel;

  INITSTATUS( status, "LALSimulatePulsarSignal", SIMULATEPULSARSIGNALC );
  ATTATCHSTATUSPTR(status);

  ASSERT ( timeSeries, status, SIMULATEPULSARSIGNAL_ENULL, SIMULATEPULSARSIGNAL_MSGENULL);
  ASSERT ( (*timeSeries)==NULL, status, SIMULATEPULSARSIGNAL_ENONULL, SIMULATEPULSARSIGNAL_MSGENONULL);
  /* don't accept heterodyning frequency */
  ASSERT ( params->fHeterodyne==0, status, SIMULATEPULSARSIGNAL_EINPUT, SIMULATEPULSARSIGNAL_MSGEINPUT);

  /* get timestamps of timeseries plus detector-states */
  dt = 1.0 / params->samplingRate;
  TRY ( LALMakeTimestamps(status->statusPtr, &timestamps, params->startTimeGPS, params->duration, dt),  status);

  numSteps = timestamps->length;

  TRY(LALGetDetectorStates(status->statusPtr, &detStates,timestamps,params->site,params->ephemerides,0), status );
  
  TRY ( LALDestroyTimestampVector (status->statusPtr, &timestamps), status );
  timestamps = NULL;

  amcoe = LALCalloc ( 1,  sizeof( *(amcoe)));
  amcoe->a = XLALCreateREAL4Vector ( numSteps );
  amcoe->b = XLALCreateREAL4Vector ( numSteps );
  TRY ( LALGetAMCoeffs (status->statusPtr, amcoe, detStates, params->pulsar.position ), status );

  /* create output timeseries (FIXME: should really know *detector* here, not just site!!) */
  if ( (channel = XLALGetChannelPrefix ( site->name )) == NULL )
    {
      LALPrintError ("\nXLALGetChannelPrefix() Failed to extract channel-prefix from site '%s'\n\n", site->name );
      ABORT (status, GENERATEPULSARSIGNALH_EDETECTOR, GENERATEPULSARSIGNALH_MSGEDETECTOR );
    }

  if ( NULL == ((*timeSeries) = XLALCreateREAL4TimeSeries( channel, &(detStates->data[0].tGPS), 0, dt, &emptyUnit, numSteps) ) )  {
    ABORT ( status, SIMULATEPULSARSIGNAL_EMEM, SIMULATEPULSARSIGNAL_MSGEMEM );
  }
  LALFree ( channel );

  /* orientation of detector arms */
  xAzi = site->xArmAzimuthRadians;
  yAzi = site->yArmAzimuthRadians;
  Zeta =  xAzi - yAzi;
  if (Zeta < 0) Zeta = -Zeta;
  if(params->site->type == LALDETECTORTYPE_CYLBAR) Zeta = LAL_PI_2;
  sinZeta = sin(Zeta);

  /* get source skyposition */
  Alpha = params->pulsar.position.longitude;
  Delta = params->pulsar.position.latitude;
  
  vn[0] = cos(Delta) * cos(Alpha);
  vn[1] = cos(Delta) * sin(Alpha);
  vn[2] = sin(Delta);

  /* get spin-parameters (restricted to maximally 3 spindowns right now) */
  phi0 = params->pulsar.phi0;
  f0   = params->pulsar.f0;

  if ( params->pulsar.spindown && (params->pulsar.spindown->length > numSpins) )
    {
      LALPrintError ("Sorry, SimulatePulsarSignal() only supports up to %d spindowns!\n", numSpins );
      ABORT (status,  SIMULATEPULSARSIGNAL_EINPUT,  SIMULATEPULSARSIGNAL_MSGEINPUT);
    }
  if ( params->pulsar.spindown && (params->pulsar.spindown->length >= 3 ) )
    f3dot = params->pulsar.spindown->data[2];
  else
    f3dot = 0.0;
  if ( params->pulsar.spindown && (params->pulsar.spindown->length >= 2 ) )
    f2dot = params->pulsar.spindown->data[1];
  else
    f2dot = 0.0;
  if ( params->pulsar.spindown && (params->pulsar.spindown->length >= 1 ) )
    f1dot = params->pulsar.spindown->data[0];
  else
    f1dot = 0.0;

  /* internally we always work with refTime = startTime->SSB, therefore
   * we need to translate the pulsar spin-params and initial phase to the
   * startTime
   */
  startTimeSSB = GPS2REAL8(detStates->data[0].tGPS) + SCALAR(vn, detStates->data[0].rDetector );
  if ( params->pulsar.refTime.gpsSeconds != 0 )
    {
      REAL8 refTime0 = GPS2REAL8(params->pulsar.refTime);
      REAL8 deltaRef = startTimeSSB - refTime0; 
      LIGOTimeGPS newEpoch;
      PulsarSpins fkdotOld, fkdotNew;
      
      XLALGPSSetREAL8( &newEpoch, startTimeSSB );

      INIT_MEM ( fkdotOld );
      fkdotOld[0] = f0;
      fkdotOld[1] = f1dot;
      fkdotOld[2] = f2dot;
      fkdotOld[3] = f3dot;

      TRY ( LALExtrapolatePulsarSpins ( status->statusPtr, fkdotNew, newEpoch, fkdotOld, params->pulsar.refTime ), status );

      /* Finally, need to propagate phase */
      phi0 += LAL_TWOPI * (               f0    * deltaRef 
			    + (1.0/2.0) * f1dot * deltaRef * deltaRef 
			    + (1.0/6.0) * f2dot * deltaRef * deltaRef * deltaRef
			    + (1.0/24.0)* f3dot * deltaRef * deltaRef * deltaRef * deltaRef 
			    );

      f0    = fkdotNew[0];
      f1dot = fkdotNew[1];
      f2dot = fkdotNew[2];
      f3dot = fkdotNew[3];

      refTime = startTimeSSB;

    } /* if refTime given */
  else /* if not given: use startTime -> SSB */
    refTime = startTimeSSB;

  /* get 4 amplitudes A_\mu */
  {
    REAL8 aPlus  = sinZeta * params->pulsar.aPlus;
    REAL8 aCross = sinZeta * params->pulsar.aCross;
    REAL8 twopsi = 2.0 * params->pulsar.psi;
  
    A1 =  aPlus * cos(phi0) * cos(twopsi) - aCross * sin(phi0) * sin(twopsi);
    A2 =  aPlus * cos(phi0) * sin(twopsi) + aCross * sin(phi0) * cos(twopsi);
    A3 = -aPlus * sin(phi0) * cos(twopsi) - aCross * cos(phi0) * sin(twopsi);
    A4 = -aPlus * sin(phi0) * sin(twopsi) + aCross * cos(phi0) * cos(twopsi);
  }

  /* main loop: generate time-series */
  for (i=0; i < detStates->length; i++)
    {
      REAL8 ai, bi;
      LIGOTimeGPS *tiGPS = &(detStates->data[i].tGPS);
      REAL8 ti, deltati, dT, taui;
      REAL8 phi_i, cosphi_i, sinphi_i;
      REAL8 hi;

      ti = GPS2REAL8 ( (*tiGPS) );
      deltati = ti - refTime;
      dT = SCALAR(vn, detStates->data[i].rDetector );
      taui = deltati + dT;

      phi_i = LAL_TWOPI * ( f0 * taui 
			    + (1.0/2.0) * f1dot * taui*taui
			    + (1.0/6.0) * f2dot * taui*taui*taui
			    + (1.0/24.0)* f3dot * taui*taui*taui*taui
			    );

      cosphi_i = cos(phi_i);
      sinphi_i = sin(phi_i);

      ai = amcoe->a->data[i];
      bi = amcoe->b->data[i];
      
      hi = A1 * ai * cosphi_i 
	+  A2 * bi * cosphi_i 
	+  A3 * ai * sinphi_i 
	+  A4 * bi * sinphi_i;

      (*timeSeries)->data->data[i] = (REAL4)hi;

    } /* for i < Nsteps */

  TRY ( LALDestroyDetectorStateSeries(status->statusPtr, &detStates ), status );
  XLALDestroyREAL4Vector(  amcoe->a );
  XLALDestroyREAL4Vector(  amcoe->b );
  LALFree ( amcoe );

  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* LALSimulateExactPulsarSignal() */


