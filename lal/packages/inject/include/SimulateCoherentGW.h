/*
 *  Copyright (C) 2007 Jolien Creighton, Reinhard Prix, Teviet Creighton, John Whelan
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

#ifndef _SIMULATECOHERENTGW_H
#define _SIMULATECOHERENTGW_H

#include <lal/LALStdlib.h>
#include <lal/DetectorSite.h>
#include <lal/SkyCoordinates.h>
#include <lal/LALBarycenter.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
   \addtogroup SimulateCoherentGW_h
   \author Creighton, T. D.

   \brief Provides routines to simulate generic gravitational waveforms
   originating from a particular source.

   \heading{Synopsis}
   \code
   #include <lal/SimulateCoherentGW.h>
   \endcode

   This header covers generic routines and structures to represent and
   simulate the effects of a plane gravitational wave propagating from a
   distinct point on the sky.

   Any plane gravitational wave is specified by a direction
   \f$\mathbf{\hat{n}}\f$ to its apparent source (i.e.\ opposite to its direction
   of propagation), and by the inistantaneous values \f$h_+(t)\f$,
   \f$h_\times(t)\f$ of its plus and cross polarizations as functions of
   (retarded) time \f$t=t_0+\mathbf{\hat{n}}\cdot(\mathbf{x}-\mathbf{x}_0)\f$, where
   \f$t_0\f$ is the time measured at some local reference point \f$\mathbf{x}_0\f$,
   and \f$t\f$ is the time measured by a synchronized clock at \f$\mathbf{x}\f$.  We
   adopt the standard meaning of the instantaneous strain amplitudes
   \f$h_{+,\times}\f$: in some reference transverse \f$x\f$-\f$y\f$ coordinate system
   oriented such that \f$\mathbf{\hat{x}}\times\mathbf{\hat{y}}=-\mathbf{\hat{n}}\f$
   points in the direction of propagation, two free observers originally
   separated by a displacement \f$(x,y)\f$ will experience an additional
   tidal displacement \f$\delta x=(xh_+ + yh_\times)/2\f$, \f$\delta
   y=(xh_\times - yh_+)/2\f$.

   \heading{Quasiperiodic waves:} Most astrophysical sources of
   gravitational radiation are described as \e quasiperiodic (or,
   less accurately, as "adiabatic"), in that they can be said to have
   an instantaneous frequency, amplitude, and polarization, all of which
   vary on timescales much longer than a wave period. Mathematically we
   write this as:
   \f[
   h_{+,\times}(t) = A_{+,\times}(t)\cos\phi(t)
   + B_{+,\times}(t)\sin\phi(t) \; ,
   \f]

   \wrapfig{r,0.47\textwidth,fig_phase_diagram}
   \image html inject_phase_diagram.png "Fig. [fig_phase_diagram]: Polarization phase diagram for a quasiperiodic gravitational wave. The phase point p(t) traces out the indicated ellipse in the h_+,h_x plane; the parameters A1, A2 and Phi remain roughly constant over many cycles in phi."
   \image latex inject_phase_diagram.eps "Polarization phase diagram for a quasiperiodic gravitational wave. The phase point p(t) traces out the indicated ellipse in the h_+,h_x plane; the parameters A1, A2 and Phi remain roughly constant over many cycles in phi." width=0.42\textwidth

   where \f$\phi(t)=2\pi\int f(t)\,dt\f$, and the <em>evolution timescale</em>
   \f$\tau=\min\{A/\dot{A},B/\dot{B},f/\dot{f}\}\f$ is much greater than
   \f$h/\dot{h}\sim1/f\f$.  Obviously it is mathematically impossible for the
   physical functions \f$h_{+,\times}(t)\f$ to specify uniquely more than two
   other functions of time; we rely on the notion of quasiperiodicity to
   define "natural" choices of instantaneous frequency and amplitude.
   The ambiguity in this choice is on the order of the amount that these
   quantities change over a cycle.

   While the above formula appears to have five degrees of freedom (two
   quadrature amplitudes \f$A\f$ and \f$B\f$ for each polarization, plus a common
   phase function \f$\phi\f$), there is a degeneracy between the two
   quadrature amplitudes and a shift in phase.  One could simply treat
   each polarization independently and represent the system with two
   amplitude functions \f$A_{+,\times}\f$ and two phase functions
   \f$\phi_{+,\times}\f$, but we would like to preserve the notion that the
   phases of the two waveforms derive from a single underlying
   instantaneous frequency.  We therefore write the waveforms in terms of
   two polarization amplitudes \f$A_1(t)\f$ and \f$A_2(t)\f$, a single phase
   function \f$\phi(t)\f$, and a polarization shift \f$\Phi(t)\f$:
   \anchor eq_quasiperiodic_hpluscross \f{eqnarray}{
   \label{eq_quasiperiodic_hpluscross}
   h_+(t) & = & A_1(t)\cos\Phi(t)\cos\phi(t)
   - A_2(t)\sin\Phi(t)\sin\phi(t) \; , \\
   h_\times(t) & = & A_1(t)\sin\Phi(t)\cos\phi(t)
   + A_2(t)\cos\Phi(t)\sin\phi(t) \; .
   \f}
   The physical meaning of these functions is shown in Fig.\figref{fig_phase_diagram}.

   There is a close relationship between the polarization shift \f$\Phi\f$
   and the orientation of the \f$x\f$-\f$y\f$ coordinates used to define our
   polarization basis: if we rotate the \f$x\f$ and \f$y\f$ axes by an angle
   \f$\Delta\psi\f$, we change \f$\Phi\f$ by an amount \f$-2\Delta\psi\f$.  (The
   factor of 2 comes from the fact that the + and \f$\times\f$ modes are
   quadrupolar: a + mode rotated \f$45^\circ\f$ is a \f$\times\f$ mode, and a
   mode rotated \f$90^\circ\f$ is the opposite of itself.)  We use the
   <em>polarization angle</em> \f$\psi\f$ to define the orientation of the
   \f$x\f$-axis of the polarization basis relative to an Earth-fixed
   reference frame (see the coordinate conventions below).  If \f$\Phi\f$ is
   constant, one can redefine \f$\psi\f$ such that \f$\Phi=0\f$; however, when
   \f$\Phi\f$ changes with time, we would nonetheless like our polarization
   basis to remain fixed.  We therefore retain the constant \f$\psi\f$ and
   the function \f$\Phi(t)\f$ as distinct quantities.

   The advantage of this quasiperiodic representation of a gravitational
   wave is that a physical sampling of the parameters \f$A_1\f$, \f$A_2\f$,
   \f$\phi\f$, and \f$\Phi\f$ need only be done on timescales \f$\Delta
   t\lesssim\tau\f$, whereas the actual wave functions \f$h_{+,\times}\f$ need
   to be sampled on timescales \f$\Delta t\lesssim1/f\f$.

   The following coordinate conventions are assumed:
   <ol>
   <li> Fig. 7 of [\ref Will_C_1996] defines standard coordinate
   conventions for nonprecessing binaries, and by extension, for any
   fixed-axis rotating source: If \f$\mathbf{\hat{z}}\f$ points in the direction
   of wave propagation (away from the source), and \f$\mathbf{\hat{l}}\f$ points
   in the (constant) direction of the source's angular momentum vector,
   then the \f$x\f$-\f$y\f$ coordinates used to define the + and \f$\times\f$
   polarizations are given by \f$\mathbf{\hat{x}}=|\csc
   i|\mathbf{\hat{z}}\times\mathbf{\hat{l}}\f$ and
   \f$\mathbf{\hat{y}}=\mathbf{\hat{z}}\times\mathbf{\hat{x}}\f$, where
   \f$i=\arccos(\mathbf{\hat{z}}\cdot\mathbf{\hat{l}})\f$ is the inclination angle
   between \f$\mathbf{\hat{l}}\f$ and \f$\mathbf{\hat{z}}\f$.  Such a system will
   generically have \f$A_1(t)=A(t)(1+\cos^2i)\f$, \f$A_2(t)=2A(t)\cos i\f$,
   \f$\Phi(t)=0\f$, and \f$f(t)>0\f$ (i.e.\ \f$\phi(t)\f$ increasing with time).  For
   precessing systems, prescriptions for \f$\mathbf{\hat{x}}\f$ and
   \f$\mathbf{\hat{y}}\f$ become ambiguous, but they \e must be fixed; the
   relations for \f$A_1\f$, \f$A_2\f$, and \f$\Phi\f$ will no longer be maintained.</li>

   <li> Appendix B of [\ref Anderson_W2000] defines a convention for
   the overal polarization angle \f$\psi\f$: Let \f$\mathbf{\hat{N}}\f$ be the
   direction of the Earth's north celestial pole, and define the
   direction of the <em>ascending node</em>
   \f$\mathbf{\hat{\Omega}}=|\csc\alpha|\mathbf{\hat{N}}\times\mathbf{\hat{z}}\f$, where
   \f$\alpha\f$ is the right ascension of the source.  Then \f$\psi\f$ is the
   angle, right-handed about \f$\mathbf{\hat{z}}\f$, from \f$\mathbf{\hat{\Omega}}\f$ to
   \f$\mathbf{\hat{x}}\f$.</li>

   <li> The direction of propagation of the wave is defined by the right
   ascension \f$\alpha\f$ and declination \f$\delta\f$ of the \e source, as
   seen from the point of measurement.  See \ref SkyCoordinates_h for a
   definition of these quantities.  We expect that these will be
   effectively constant for almost any gravitational wave source of
   interest.</li>
   </ol>

   \heading{The polarization response:} The relative strain induced in
   the test masses of a detector by a passing gravitational wave depends
   not only on the amplitudes \f$h_{+,\times}\f$ of the gravitational wave,
   but also on the design of the detector and its orientation with
   relative to the \f$x\f$-\f$y\f$ coordinate system used to define the + and
   \f$\times\f$ polarizations.  For a given detector, the response to each
   polarization thus depends on the right ascension \f$\alpha\f$, declination
   \f$\delta\f$, and polarization angle \f$\psi\f$ of the source (which define
   the orientation of the + and \f$\times\f$ polarization axes relative to
   the Earth), and on the time \f$t\f$ (which determines the orientation of
   the detector as the Earth rotates).  The strain \f$h(t)\f$ induced in the
   detector is thus given by two polarization response functions
   \f$F_{+,\times}(\alpha,\delta,\psi;t)\f$ by:
   \f[
   h(t) = h_+(t)F_+(\alpha,\delta,\psi;t) +
   h_\times(t)F_\times(\alpha,\delta,\psi;t) \; .
   \f]
   We will not discuss the computation of these functions \f$F_{+,\times}\f$,
   as these are covered under the header \ref DetResponse_h.

   \heading{The transfer function:} All gravitational wave detectors
   incorporate a set of analog and digital filters that convert a
   gravitational excitation on the test masses into a measurable output
   time series.  The effects of these functions are aggregated into a
   complex-valued <em>transfer function</em> \f${\cal T}(f)\f$, which gives the
   instrumental response (in units of "counts" from an
   analog\f$\rightarrow\f$digital converter) to gravitational waves of unit
   amplitued at the frequency \f$f\f$.  Specifically, if the strain exerted
   on the antenna is given by \f$h(t)=\mathrm{Re}[{\cal H}e^{2\pi ift}]\f$
   (where the complex amplitude \f$\cal H\f$ includes the phase of the wave),
   then the ADC output of the instrument is given by:
   \f[
   o(t) = \mathrm{Re}\left[ {\cal T}(f) {\cal H}e^{2\pi ift} \right] \; .
   \f]
   The transfer function has a strong frequency dependence in order to
   "whiten" the highly-coloured instrumental noise, and thus preserve
   instrumental sensitivity across a broad band of frequencies.

   We note that although the transfer function measures the response of
   the instrument to a gravitational wave, the term <em>response
   function</em> refers to inverse transformation of taking an instrumental
   response and computing a gravitational waveform; that is, \f${\cal
   R}(f)=1/{\cal T}(f)\f$.  This confusing bit of nomenclature arises from
   the fact that most data analysis deals with extracting gravitational
   waveforms from the instrumental output, rather than injecting
   waveforms into the output.

   For quasiperiodic waveforms with a well-defined instantaneous
   frequency \f$f(t)\f$ and phase \f$\phi(t)\f$, we can compute the response of
   the instrument entirely in the time domain in the adiabatic limit: if
   our instrumental excitation is a linear superposition of waveforms
   \f$h(t)=\mathrm{Re}[{\cal H}(t)e^{i\phi(t)}]\f$, then the output is a
   superposition of waves of the form
   \f[
   o(t) \approx \mathrm{Re}\left[ {\cal T}\{f(t)\}
   {\cal H}(t)e^{i\phi(t)} \right] \; .
   \f]
   This expression is approximate to the extent that \f${\cal T}(f)\f$ varies
   over the range \f$f\pm1/\tau\f$, where \f$\tau\f$ is the evolution timescale
   of \f${\cal H}(t)\f$ and \f$f(t)\f$.  Since the transfer function and
   polarization response (above) are linear operators, we can apply them
   in either order.

   \heading{A note on terminology:} We use the word "coherent" in the
   name of this header in the loosest possible sense, refering to any
   wave with a well-defined direction of propagation, whose wave
   amplitudes \f$h_{+,\times}\f$ are deterministic functions of retarded
   time.  Given a knowledge of these parameters, such a waveform is
   amenable to "coherent" detection in a network of detectors, through
   time-shifted matched filtering.

   However, coherence is often used to refer to a more restricted class
   of waveforms that are "effectively monochromatic" over some
   coherence timescale \f$t_\mathrm{coh}\f$; i.e.\ in any timespan
   \f$t_\mathrm{coh}\f$ there is a fixed-frequency sinusoid that is never
   more than \f$90^\circ\f$ out of phase with the waveform.  This is more
   retrictive even than our concept of quasiperiodic waves; for
   smoothly-varying waveforms one has \f$t_\mathrm{coh}\sim\dot{f}^{-1/2}\f$,
   which is much shorter than the evolution timescale \f$\tau\sim
   f/\dot{f}\f$ (provided \f$\tau\gg1/f\f$, as we have assumed).

*/
/*@{*/

/** This structure stores a representation of a plane
 * gravitational wave propagating from a particular point on the sky.
 * Several alternate representations are permitted to allow a more
 * natural characterization of quasiperiodic waveforms.
 *
 * \note It is permissible to set only some of the
 * \c REAL4TimeSeries or \c REAL4TimeVectorSeries fields above,
 * but the waveform is treated as being zero except during those times
 * when either \c h, or both \c a and \c phi, are defined.
 * Where \c shift is not specified, it is assumed that \f$\Phi\f$ is
 * zero; where \c f is not specified but \c phi is, \f$f(t)\f$ can be
 * computed as \f$\dot{\phi}(t)/2\pi\f$.  Where \c f and \c phi
 * overlap, or where \c h and any other time series overlap, they
 * must be defined consistently.
 *
 */
typedef struct tagCoherentGW {
  SkyPosition position;     /**< The location of the source in the sky; this should be in equatorial celestial coordinates, but routines may be able to do the conversion */
  REAL4 psi;                /**< The polarization angle \f$\psi\f$, in radians, as defined in Appendix B of [\ref Anderson_W2000] */
  REAL4TimeVectorSeries *h; /**< A time-sampled two-dimensional vector storing the waveforms \f$h_+(t)\f$ and \f$h_\times(t)\f$, in dimensionless strain */
  REAL4TimeVectorSeries *a; /**< A time-sampled two-dimensional vector storing the amplitudes \f$A_1(t)\f$ and \f$A_2(t)\f$, in dimensionless strain */
  REAL4TimeSeries *f;       /**< A time-sampled sequence storing the instantaneous frequency \f$f(t)\f$, in Hz. */
  REAL8TimeSeries *phi;     /**< A time-sampled sequence storing the phase function \f$\phi(t)\f$, in radians */
  REAL4TimeSeries *shift;   /**< A time-sampled sequence storing the polarization shift \f$\Phi(t)\f$, in radians */
  UINT4 dtDelayBy2;         /**< A user specified half-interval time step for the Doppler delay look-up table (will default to 400s if set to 0) */
  UINT4 dtPolBy2;           /**< A user defined half-interval time step for the polarisation response look-up table (will default to 300s if set to 0) */
} CoherentGW;

/** This structure contains information required to determine the response
 * of a detector to a gravitational waveform.
 */
typedef struct tagDetectorResponse {
  COMPLEX8FrequencySeries *transfer;    /**< The frequency-dependent transfer function of the interferometer, in ADC counts per unit strain amplitude at any given frequency;
                                         * if absent, the response will be given in raw strain rather than ADC output */
  LALDetector *site;                    /**< A structure storing site and polarization information, used to compute the polarization response and the propagation delay;
                                         * if absent, the response will be computed to the plus mode waveform with no time delay */
  EphemerisData *ephemerides;           /**< A structure storing the positions, velocities, and accelerations of the Earth and Sun centres of mass, used to compute
                                         * the propagation delay to the solar system barycentre;
                                         * if absent, the propagation delay will be computed to the Earth centre (rather than a true barycentre) */
  LIGOTimeGPS heterodyneEpoch;          /**< A reference time for heterodyned detector output time series, where the phase of the mixing signal is zero.
                                         * This parameter is only used when generating detector output time series with nonzero heterodyne frequency \c f0.
                                         * (Note: This should really be a parameter stored in the \c TimeSeries structure along with \c f0, but it isnt, so we
                                         * have to add it here.)
                                         */
} DetectorResponse;

/* Function prototypes. */

int XLALSimulateCoherentGW ( REAL4TimeSeries  *output, CoherentGW *CWsignal, DetectorResponse *detector );

void
LALSimulateCoherentGW( LALStatus        *status,
                       REAL4TimeSeries  *output,
                       CoherentGW       *input,
                       DetectorResponse *detector );

/*@}*/

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _SIMULATECOHERENTGW_H */
