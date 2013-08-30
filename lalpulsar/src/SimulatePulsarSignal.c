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

#include <lal/AVFactories.h>
#include <lal/TimeSeries.h>
#include <lal/GeneratePulsarSignal.h>
#include <lal/ComputeFstat.h>
#include <lal/ExtrapolatePulsarSpins.h>
#include <lal/SimulatePulsarSignal.h>

/*----- Macros ----- */

/** Simple Euklidean scalar product for two 3-dim vectors in cartesian coords */
#define SCALAR(u,v) ((u)[0]*(v)[0] + (u)[1]*(v)[1] + (u)[2]*(v)[2])
#define INIT_MEM(x) memset(&(x), 0, sizeof((x)))

// ----- global variables, initializers
static LALUnit emptyUnit;


// ----- function definitions


/**
 * Simulate a pulsar signal to best accuracy possible.
 * \author Reinhard Prix
 * \date 2005
 *
 * The motivation for this function is to provide functions to
 * simulate pulsar signals <em>with the best possible accuracy</em>,
 * i.e. using no approximations, contrary to LALGeneratePulsarSignal().
 *
 * Obviously this is not meant as a fast code to be used in a Monte-Carlo
 * simulation, but rather as a <em>reference</em> to compare other (faster)
 * functions agains, in order to be able to gauge the quality of a given
 * signal-generation routine.
 *
 * We want to calculate \f$h(t)\f$, given by
 * \f[
 * h(t) = F_+(t)\, h_+(t) + F_\times(t) \,h_\times(t)\,,
 * \f]
 * where \f$F_+\f$ and \f$F_x\f$ are called the <em>beam-pattern</em> functions,
 * which depend of the wave polarization \f$\psi\f$,
 * the source position \f$\alpha\f$, \f$\delta\f$ and the detector position and
 * orientation (\f$\gamma\f$, \f$\lambda\f$, \f$L\f$ and \f$\xi\f$). The expressions for
 * the beam-pattern functions are given in [\ref JKS98], which we write as
 * \f{eqnarray}{
 * F_+(t) = \sin \zeta \cos 2\psi \, a(t)  + \sin \zeta \sin 2\psi \, b(t)\,,\\
 * F_\times(t) = \sin\zeta  \cos 2\psi \,b(t) - \sin\zeta \sin 2\psi \, a(t) \,,
 * \f}
 * where \f$\zeta\f$ is the angle between the interferometer arms, and
 * \f{eqnarray}{
 * a(t) &=& a_1 \cos[ 2 (\alpha - T)) ] + a_2 \sin[ 2(\alpha - T)]
 * + a_3 \cos[ \alpha - T ] + a_4 \sin [ \alpha - T ] + a_5\,,\\
 * b(t) &=& b_1 \cos[ 2(\alpha - T)] + b_2 \sin[ 2(\alpha - T) ]
 * + b_3 \cos[ \alpha - T ] + b_4 \sin[ \alpha - T]\,,
 * \f}
 * where \f$T\f$ is the local (mean) sidereal time of the detector, and the
 * time-independent coefficients \f$a_i\f$ and \f$b_i\f$ are given by
 * \f{eqnarray}{
 * a_1 &=& \frac{1}{16} \sin 2\gamma \,(3- \cos 2\lambda)\,(3 - \cos 2\delta)\,,\\
 * a_2 &=& -\frac{1}{4}\cos 2\gamma \,\sin \lambda \,(3 - \cos 2\delta) \,,\\
 * a_3 &=& \frac{1}{4} \sin 2\gamma \,\sin 2\lambda \,\sin 2\delta  \,\\
 * a_4 &=& -\frac{1}{2} \cos 2\gamma \,\cos \lambda \,\sin 2 \delta\,,\\
 * a_5 &=& \frac{3}{4} \sin 2\gamma \, \cos^2 \lambda \,\cos^2 \delta\,,
 * \f}
 * and
 * \f{eqnarray}{
 * b_1 &=& \cos 2\gamma \,\sin \lambda \,\sin \delta\,,\\
 * b_2 &=& \frac{1}{4} \sin 2\gamma \,(3-\cos 2\lambda)\, \sin \delta\,,\\
 * b_3 &=& \cos 2\gamma \,\cos \lambda \,\cos\delta \,, \\
 * b_4 &=& \frac{1}{2} \sin2\gamma \,\sin 2\lambda \,\cos\delta\,,
 * \f}
 *
 * The source model considered is a plane-wave
 * \f{eqnarray}{
 * h_+(t) &=& A_+\, \cos \Psi(t)\,,\\
 * h_\times(t) &=& A_\times \, \sin \Psi(t)\,,
 * \f}
 * where the wave-phase is \f$\Psi(t) = \Phi_0 + \Phi(t)\f$, and for an
 * isolated pulsar we have
 * \f{equation}{
 * \Phi(t) = 2\pi \left[\sum_{s=0} \frac{f^{(s)}(\tau_\mathrm{ref})}{
 * (s+1)!} \left( \tau(t) - \tau_\mathrm{ref} \right)^{s+1} \right]\,,
 * \f}
 * where \f$\tau_\mathrm{ref}\f$ is the "reference time" for the definition
 * of the pulsar-parameters \f$f^{(s)}\f$ in the solar-system barycenter
 * (SSB), and \f$\tau(t)\f$ is the SSB-time of the phase arriving at the
 * detector at UTC-time \f$t\f$, which depends on the source-position
 * (\f$\alpha\f$, \f$\delta\f$) and the detector-position, namely
 * \f{equation}{
 * \tau (t) = t + \frac{ \vec{r}(t)\cdot\vec{n}}{c}\,,
 * \f}
 * where \f$\vec{r}(t)\f$ is the vector from SSB to the detector, and \f$\vec{n}\f$
 * is the unit-vector pointing <em>to</em> the source.
 *
 * This is a standalone "clean-room" implementation using no other
 * outside-functions <em>except</em> for LALGPStoLMST1() to calculate
 * the local (mean) sidereal time at the detector for given GPS-time,
 * (which I double-checked with an independent Mathematica script),
 * and and LALBarycenter() to calculate \f$\tau(t)\f$.
 *
 * NOTE: currently only isolated pulsars are supported
 *
 * NOTE2: we don't really use the highest possible accuracy right now,
 * as we blatently neglect all relativistic timing effects (i.e. using dT=v.n/c)
 *
 * NOTE3: no heterodyning is performed here, the time-series is generated and sampled
 * at the given rate, that's all! ==\> the caller needs to make sure about the
 * right sampling rate to use (-\>aliasing) and do the proper post-treatment...
 *
 */
REAL4TimeSeries *
XLALSimulateExactPulsarSignal ( const PulsarSignalParams *params )
{
  XLAL_CHECK_NULL ( params != NULL, XLAL_EINVAL, "Invalid NULL input 'params'\n");
  XLAL_CHECK_NULL ( params->samplingRate > 0, XLAL_EDOM, "Sampling rate must be positive, got samplingRate = %g\n", params->samplingRate );

  /* don't accept heterodyning frequency */
  XLAL_CHECK_NULL ( params->fHeterodyne == 0, XLAL_EINVAL, "Heterodyning frequency must be set to 0, got params->fHeterodyne = %g\n", params->fHeterodyne );

  UINT4 numSpins = 3;

  /* get timestamps of timeseries plus detector-states */
  REAL8 dt = 1.0 / params->samplingRate;
  LIGOTimeGPSVector *timestamps;
  XLAL_CHECK_NULL ( (timestamps = XLALMakeTimestamps ( params->startTimeGPS, params->duration, dt, 0 )) != NULL, XLAL_EFUNC );

  UINT4 numSteps = timestamps->length;

  DetectorStateSeries *detStates = XLALGetDetectorStates ( timestamps, params->site, params->ephemerides, 0 );
  XLAL_CHECK_NULL ( detStates != NULL, XLAL_EFUNC, "XLALGetDetectorStates() failed.\n");

  XLALDestroyTimestampVector ( timestamps );
  timestamps = NULL;

  AMCoeffs *amcoe = XLALComputeAMCoeffs ( detStates, params->pulsar.position );
  XLAL_CHECK_NULL ( amcoe != NULL, XLAL_EFUNC, "XLALComputeAMCoeffs() failed.\n");

  /* create output timeseries (FIXME: should really know *detector* here, not just site!!) */
  const LALFrDetector *site = &(params->site->frDetector);
  CHAR *channel = XLALGetChannelPrefix ( site->name );
  XLAL_CHECK_NULL ( channel != NULL, XLAL_EFUNC, "XLALGetChannelPrefix( %s ) failed.\n", site->name );

  REAL4TimeSeries *ts = XLALCreateREAL4TimeSeries ( channel, &(detStates->data[0].tGPS), 0, dt, &emptyUnit, numSteps );
  XLAL_CHECK_NULL ( ts != NULL, XLAL_EFUNC, "XLALCreateREAL4TimeSeries() failed.\n");
  XLALFree ( channel );
  channel = NULL;

  /* orientation of detector arms */
  REAL8 xAzi = site->xArmAzimuthRadians;
  REAL8 yAzi = site->yArmAzimuthRadians;
  REAL8 Zeta =  xAzi - yAzi;
  if (Zeta < 0) {
    Zeta = -Zeta;
  }
  if ( params->site->type == LALDETECTORTYPE_CYLBAR ) {
    Zeta = LAL_PI_2;
  }
  REAL8 sinZeta = sin(Zeta);

  /* get source skyposition */
  REAL8 Alpha = params->pulsar.position.longitude;
  REAL8 Delta = params->pulsar.position.latitude;
  REAL8 vn[3];
  vn[0] = cos(Delta) * cos(Alpha);
  vn[1] = cos(Delta) * sin(Alpha);
  vn[2] = sin(Delta);

  /* get spin-parameters (restricted to maximally 3 spindowns right now) */
  REAL8 phi0 = params->pulsar.phi0;
  REAL8 f0   = params->pulsar.f0;

  REAL8 f1dot = 0, f2dot = 0, f3dot = 0;
  if ( params->pulsar.spindown && (params->pulsar.spindown->length > numSpins) ) {
    XLAL_ERROR_NULL ( XLAL_EDOM, "Currently only supports up to %d spindowns!\n", numSpins );
  }
  if ( params->pulsar.spindown && (params->pulsar.spindown->length >= 3 ) ) {
    f3dot = params->pulsar.spindown->data[2];
  }
  if ( params->pulsar.spindown && (params->pulsar.spindown->length >= 2 ) ) {
    f2dot = params->pulsar.spindown->data[1];
  }
  if ( params->pulsar.spindown && (params->pulsar.spindown->length >= 1 ) ) {
    f1dot = params->pulsar.spindown->data[0];
  }

  /* internally we always work with refTime = startTime->SSB, therefore
   * we need to translate the pulsar spin-params and initial phase to the
   * startTime
   */
  REAL8 startTimeSSB = XLALGPSGetREAL8 ( &(detStates->data[0].tGPS) ) + SCALAR ( vn, detStates->data[0].rDetector );
  REAL8 refTime;
  if ( params->pulsar.refTime.gpsSeconds != 0 )
    {
      REAL8 refTime0 = XLALGPSGetREAL8 ( &(params->pulsar.refTime) );
      REAL8 deltaRef = startTimeSSB - refTime0;
      LIGOTimeGPS newEpoch;
      PulsarSpins fkdotOld, fkdotNew;

      XLALGPSSetREAL8( &newEpoch, startTimeSSB );

      INIT_MEM ( fkdotOld );
      fkdotOld[0] = f0;
      fkdotOld[1] = f1dot;
      fkdotOld[2] = f2dot;
      fkdotOld[3] = f3dot;
      REAL8 DeltaTau = XLALGPSDiff ( &newEpoch, &(params->pulsar.refTime) );

      int ret = XLALExtrapolatePulsarSpins ( fkdotNew, fkdotOld, DeltaTau );
      XLAL_CHECK_NULL ( ret == XLAL_SUCCESS, XLAL_EFUNC, "XLALExtrapolatePulsarSpins() failed.\n");

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
  else  { /* if not given: use startTime -> SSB */
    refTime = startTimeSSB;
  }

  /* get 4 amplitudes A_\mu */
  REAL8 aPlus  = sinZeta * params->pulsar.aPlus;
  REAL8 aCross = sinZeta * params->pulsar.aCross;
  REAL8 twopsi = 2.0 * params->pulsar.psi;

  REAL8 A1 =  aPlus * cos(phi0) * cos(twopsi) - aCross * sin(phi0) * sin(twopsi);
  REAL8 A2 =  aPlus * cos(phi0) * sin(twopsi) + aCross * sin(phi0) * cos(twopsi);
  REAL8 A3 = -aPlus * sin(phi0) * cos(twopsi) - aCross * cos(phi0) * sin(twopsi);
  REAL8 A4 = -aPlus * sin(phi0) * sin(twopsi) + aCross * cos(phi0) * cos(twopsi);

  /* main loop: generate time-series */
  for ( UINT4 i = 0; i < numSteps; i++)
    {
      LIGOTimeGPS *tiGPS = &(detStates->data[i].tGPS);

      REAL8 ti = XLALGPSGetREAL8 ( tiGPS );
      REAL8 deltati = ti - refTime;
      REAL8 dT = SCALAR(vn, detStates->data[i].rDetector );
      REAL8 taui = deltati + dT;

      REAL8 phi_i = LAL_TWOPI * ( f0 * taui
			    + (1.0/2.0) * f1dot * taui*taui
			    + (1.0/6.0) * f2dot * taui*taui*taui
			    + (1.0/24.0)* f3dot * taui*taui*taui*taui
			    );

      REAL8 cosphi_i = cos(phi_i);
      REAL8 sinphi_i = sin(phi_i);

      REAL8 ai = amcoe->a->data[i];
      REAL8 bi = amcoe->b->data[i];

      REAL8 hi = A1 * ai * cosphi_i
	+  A2 * bi * cosphi_i
	+  A3 * ai * sinphi_i
	+  A4 * bi * sinphi_i;

      ts->data->data[i] = (REAL4)hi;

    } /* for i < numSteps */

  XLALDestroyDetectorStateSeries( detStates );
  XLALDestroyAMCoeffs ( amcoe );

  return ts;

} /* XLALSimulateExactPulsarSignal() */
