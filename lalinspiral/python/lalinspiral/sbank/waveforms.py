# Copyright (C) 2011  Nickolas Fotopoulos
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

from __future__ import division

from math import isnan
import numpy as np
np.seterr(all="ignore")

import lal
import lalsimulation as lalsim
from lal import MSUN_SI, MTSUN_SI, PC_SI, PI, CreateREAL8Vector, CreateCOMPLEX8FrequencySeries
from lalsimulation import SimInspiralTaylorF2RedSpinComputeNoiseMoments, SimInspiralTaylorF2RedSpinMetricChirpTimes, SimIMRPhenomPCalculateModelParameters
from lalinspiral import InspiralSBankComputeMatch
from lalinspiral.sbank.psds import get_neighborhood_PSD
from lalinspiral.sbank.tau0tau3 import m1m2_to_tau0tau3

from pylal import spawaveform  # XXX: Remove when everything is ported to lalsim
from pylal.xlal.datatypes.snglinspiraltable import SnglInspiralTable

def compute_mchirp(m1, m2):
    return (m1 * m1 * m1 * m2 * m2 * m2 / (m1 + m2))**0.2

def ceil_pow_2( number ):
    return int(2**(np.ceil(np.log2( number ))))

#
# Represent template waveforms
#

def project_hplus_hcross(hplus, hcross, theta, phi, psi):
    # compute antenna factors Fplus and Fcross
    Fp = 0.5*(1 + np.cos(theta)**2)*np.cos(2*phi)*np.cos(2*psi) - np.cos(theta)*np.sin(2*phi)*np.sin(2*psi)
    Fc = 0.5*(1 + np.cos(theta)**2)*np.cos(2*phi)*np.sin(2*psi) + np.cos(theta)*np.sin(2*phi)*np.cos(2*psi)

    # form strain signal in detector
    # hoft = lal.CreateCOMPLEX16FrequncySeries("h(t)", hplus.epoch, hplus.f0, hplus.deltaF, lal.lalSecondUnit, hplus.data.length)
    hplus.data.data = Fp*hplus.data.data + Fc*hcross.data.data

    return hplus


def compute_sigmasq(htilde, deltaF):
    """
    Find norm of whitened h(f) array.
    """
    # vdot is dot with complex conjugation
    return float(np.vdot(htilde, htilde).real * 4 * deltaF)


def FrequencySeries_to_COMPLEX8FrequencySeries(fs):
    """
    Create a new COMPLEX8FrequencySeries and copy the results of the input FrequencySeries to it.
    """
    new = CreateCOMPLEX8FrequencySeries(fs.name, fs.epoch, fs.f0, fs.deltaF, fs.sampleUnits, fs.data.length)
    new.data.data[:] = fs.data.data[:]  # numpy automagic conversion
    return new

def create_moments(df, flow, len_PSD):
    # moments go from flow to fmax unlike the PSD, which starts at 0
    n = len_PSD - int(flow / df)
    momI_0 = lal.CreateREAL8Vector(n)
    momI_2 = lal.CreateREAL8Vector(n)
    momI_3 = lal.CreateREAL8Vector(n)
    momI_4 = lal.CreateREAL8Vector(n)
    momI_5 = lal.CreateREAL8Vector(n)
    momI_6 = lal.CreateREAL8Vector(n)
    momI_7 = lal.CreateREAL8Vector(n)
    momI_8 = lal.CreateREAL8Vector(n)
    momI_9 = lal.CreateREAL8Vector(n)
    momI_10 = lal.CreateREAL8Vector(n)
    momI_11 = lal.CreateREAL8Vector(n)
    momI_12 = lal.CreateREAL8Vector(n)
    momI_13 = lal.CreateREAL8Vector(n)
    momI_14 = lal.CreateREAL8Vector(n)
    momI_15 = lal.CreateREAL8Vector(n)
    momI_16 = lal.CreateREAL8Vector(n)
    momJ_5 = lal.CreateREAL8Vector(n)
    momJ_6 = lal.CreateREAL8Vector(n)
    momJ_7 = lal.CreateREAL8Vector(n)
    momJ_8 = lal.CreateREAL8Vector(n)
    momJ_9 = lal.CreateREAL8Vector(n)
    momJ_10 = lal.CreateREAL8Vector(n)
    momJ_11 = lal.CreateREAL8Vector(n)
    momJ_12 = lal.CreateREAL8Vector(n)
    momJ_13 = lal.CreateREAL8Vector(n)
    momJ_14 = lal.CreateREAL8Vector(n)
    momK_10 = lal.CreateREAL8Vector(n)
    momK_11 = lal.CreateREAL8Vector(n)
    momK_12 = lal.CreateREAL8Vector(n)
    return (momI_0, momI_2, momI_3, momI_4, momI_5, momI_6, momI_7, momI_8, momI_9, momI_10, momI_11, momI_12, momI_13, momI_14, momI_15, momI_16, momJ_5, momJ_6, momJ_7, momJ_8, momJ_9, momJ_10, momJ_11, momJ_12, momJ_13, momJ_14, momK_10, momK_11, momK_12)

PI_p5 = PI**5
def compute_chirptimes(mchirp, eta, chi, flow):
    theta0 = (125. / 2. / (16. * PI * flow * mchirp * MTSUN_SI)**5)**(1./3)
    theta3 = (16. * PI_p5 / 25. * theta0 * theta0 / (eta * eta * eta))**0.2
    theta3s = 113. / (48 * PI) * chi * theta3
    return theta0, theta3, theta3s

def compute_tau0(mc, flow):
    return 5. * mc * MTSUN_SI / (256 * (PI * flow * mc * MTSUN_SI)**(8./3))


class Template(object):
    """
    Base class that handles whitening, normalization, and
    zero-padding/unpadding. Derived classes should:

    * Add __slots__ for whatever new data they need to store.
    * Add param_names and param_formats tuples to the class describing
      the parameters that uniquely specify a template.
    * Store template information during __init__.
    * Call Template.__init__(self) during __init__.
    * Provide a _compute_waveform method that takes df and f_final to generate
      the frequency-domain waveform as a COMPLEX16FrequencySeries.
    * Provide a classmethod from_sngl, which creates an instance based on
      a sngl_inspiral object.
    * Provide a classmethod to_sngl, which creates a sngl_inspiral object
      from an instance of the subclass.
    * Provide a classmethod from_sim, which creates an instance based on
      a sim_inspiral object.
    """
    __slots__ = ("m1", "m2", "bank", "_mchirp", "_tau0", "_wf", "_metric", "sigmasq", "is_seed_point", "_f_final", "_fhigh_max")
    param_names = ("m1", "m2")
    param_formats = ("%.2f", "%.2f")

    def __init__(self, m1, m2, bank):

        self.m1 = float(m1)
        self.m2 = float(m2)
        self.bank = bank

        self._wf = {}
        self._metric = None
        self.sigmasq = 0.
        self._mchirp = compute_mchirp(m1, m2)
        self._tau0 = compute_tau0( self._mchirp, bank.flow)
        self._f_final = None
        self._fhigh_max = bank.fhigh_max


    @property
    def params(self):
        return tuple(getattr(self, k) for k in self.param_names)

    @property
    def f_final(self):

        if self._f_final is None:
            f_final = self._get_f_final()
            if self._fhigh_max:
                f_final = min(f_final, self._fhigh_max)
            self._f_final = ceil_pow_2(f_final)

        return self._f_final

    def _get_f_final(self):
        err_msg = "Template classes must provide a function _get_f_final "
        err_msg += "to compute the appropriate final frequency."
        raise NotImplementedError(err_msg)

    def __repr__(self):
        return "(%s)" % ", ".join(self.param_formats) % self.params

    def finalize_as_template(self):
        """
        This template is being added to a bank. Compute metric coefficients
        or do whatever other expensive things you didn't want to do when
        it was just a proposal.
        """
        pass

    def get_whitened_normalized(self, df, ASD=None, PSD=None):
        """
        Return a COMPLEX8FrequencySeries of the waveform, whitened by the
        given ASD and normalized. The waveform is not zero-padded to
        match the length of the ASD, so its normalization depends on
        its own length.
        """
        if not self._wf.has_key(df):
            wf = self._compute_waveform(df, self.f_final)
            if ASD is None:
                ASD = PSD**0.5
            if wf.data.length > len(ASD):
                raise ValueError("waveform has length greater than ASD; cannot whiten")
            arr_view = wf.data.data

            # whiten
            arr_view[:] /= ASD[:wf.data.length]
            arr_view[:int(self.bank.flow / df)] = 0.
            arr_view[int(self.f_final/df) : wf.data.length] = 0.

            # normalize
            self.sigmasq = compute_sigmasq(arr_view, df)
            arr_view[:] /= self.sigmasq**0.5

            # down-convert to single precision
            self._wf[df] = FrequencySeries_to_COMPLEX8FrequencySeries(wf)
        return self._wf[df]

    def metric_match(self, other, df, **kwargs):
        raise NotImplementedError

    def brute_match(self, other, df, workspace_cache, **kwargs):
        return InspiralSBankComputeMatch(self.get_whitened_normalized(df, **kwargs), other.get_whitened_normalized(df, **kwargs), workspace_cache)

    def clear(self):
        self._wf = {}


class AlignedSpinTemplate(Template):
    """
    A convenience class for aligned spin templates. Specific
    implementations of aligned spin templates should sub-class this
    class.
    """
    param_names = ("m1", "m2", "spin1z", "spin2z")
    param_formats = ("%.2f", "%.2f", "%.2f", "%.2f")
    __slots__ = ("spin1z", "spin2z", "chieff")

    def __init__(self, m1, m2, spin1z, spin2z, bank):

        Template.__init__(self, m1, m2, bank)
        self.spin1z = float(spin1z)
        self.spin2z = float(spin2z)
        self.chieff = lalsim.SimIMRPhenomBComputeChi(m1, m2, spin1z, spin2z)

    @classmethod
    def from_sim(cls, sim, bank):
        return cls(sim.mass1, sim.mass2, sim.spin1z, sim.spin2z, bank)

    @classmethod
    def from_sngl(cls, sngl, bank):
        return cls(sngl.mass1, sngl.mass2, sngl.spin1z, sngl.spin2z, bank)

    def to_sngl(self):
        # note that we use the C version; this causes all numerical
        # values to be initiated as 0 and all strings to be '', which
        # is nice
        row = SnglInspiralTable()

        row.mass1 = self.m1
        row.mass2 = self.m2
        row.mtotal = self.m1 + self.m2
        row.mchirp = self._mchirp
        row.eta = row.mass1 * row.mass2 / (row.mtotal * row.mtotal)
        row.tau0, row.tau3 = m1m2_to_tau0tau3(self.m1, self.m2, self.bank.flow)
        row.f_final = self.f_final
        row.template_duration = self._dur
        row.spin1z = self.spin1z
        row.spin2z = self.spin2z
        row.sigmasq = self.sigmasq

        return row


class IMRPhenomBTemplate(AlignedSpinTemplate):

    param_names = ("m1", "m2", "chieff")
    param_formats = ("%.2f", "%.2f", "%+.2f")
    slots = ("_dur")

    def __init__(self, m1, m2, spin1z, spin2z, bank):

        AlignedSpinTemplate.__init__(self, m1, m2, spin1z, spin2z, bank)
        self._dur = self._imrdur()

    def _get_f_final(self):
        return spawaveform.imrffinal(self.m1, self.m2, self.chieff)  # ISCO

    def _compute_waveform(self, df, f_final):
        return lalsim.SimIMRPhenomBGenerateFD(0, df,
            self.m1 * MSUN_SI, self.m2 * MSUN_SI,
            self.chieff, self.bank.flow, f_final, 1000000 * PC_SI)

    def _imrdur(self):
        """
        Ajith gave us the heuristic that chirp + 1000 M is a conservative
        estimate for the length of a full IMR waveform.
        """
        return lalsim.SimInspiralTaylorF2ReducedSpinChirpTime(self.bank.flow,
            self.m1 * MSUN_SI, self.m2 * MSUN_SI, self.chieff,
            7) + 1000 * (self.m1 + self.m2) * MTSUN_SI


class IMRPhenomCTemplate(IMRPhenomBTemplate):

    def _compute_waveform(self, df, f_final):
        return lalsim.SimIMRPhenomCGenerateFD(0, df,
            self.m1 * MSUN_SI, self.m2 * MSUN_SI,
            self.chieff, self.bank.flow, f_final, 1000000 * PC_SI)


class SEOBNRv2Template(AlignedSpinTemplate):

    param_names = ("m1", "m2", "spin1z", "spin2z")
    param_formats = ("%.2f", "%.2f", "%.2f", "%.2f")
    __slots__ = ("_dur")

    def __init__(self, m1, m2, spin1z, spin2z, bank):

        AlignedSpinTemplate.__init__(self, m1, m2, spin1z, spin2z, bank)
        self._dur = self._imrdur()

    def _get_f_final(self):
        return lalsim.SimInspiralGetFrequency(self.m1*lal.MSUN_SI,
                                self.m2*lal.MSUN_SI, 0., 0., self.spin1z, 0, 0,
                                self.spin2z, lalsim.fSEOBNRv2RD)

    def _compute_waveform(self, df, f_final):
        """
        Since SEOBNRv1 is a time domain waveform, we have to generate it,
        then FFT to frequency domain.
        """
        # need to compute dt from df, duration
        sample_rate = 2**np.ceil(np.log2(2*f_final))
        dt = 1. / sample_rate
        # get hplus
        hplus, hcross = lalsim.SimIMRSpinAlignedEOBWaveform(
            0., dt, self.m1 * MSUN_SI, self.m2 * MSUN_SI,
            self.bank.flow, 1e6*PC_SI, 0., self.spin1z, self.spin2z, 2) # last argument is SEOBNR version
        # zero-pad up to 1/df
        N = int(sample_rate / df)
        hplus = lal.ResizeREAL8TimeSeries(hplus, 0, N)
        # taper
        lalsim.SimInspiralREAL8WaveTaper(hplus.data, lalsim.SIM_INSPIRAL_TAPER_START)

        # create vector to hold output and plan
        htilde = lal.CreateCOMPLEX16FrequencySeries("h(f)", hplus.epoch, hplus.f0, df, lal.HertzUnit, int(N/2 + 1))
        fftplan = lal.CreateForwardREAL8FFTPlan(N, 0)

        # do the fft
        lal.REAL8TimeFreqFFT(htilde, hplus, fftplan)

        return htilde

    def _imrdur(self):
        """
        Following Jolien's suggestion for the duration of an IMR
        waveform, we compute the chirp time, scale by 10% and pad by
        one second.
        """
        # FIXME: This should be done better when possible!
        dur = lalsim.SimInspiralTaylorF2ReducedSpinChirpTime(self.bank.flow,
            self.m1 * MSUN_SI, self.m2 * MSUN_SI, self.chieff,
            7) + 1000 * (self.m1 + self.m2) * MTSUN_SI
        # FIXME: Is a minimal time of 1.0s too long?
        dur = 1.1 * dur + 1.0
        return dur


class SEOBNRv2ROMDoubleSpinTemplate(SEOBNRv2Template):
    def _compute_waveform(self, df, f_final):
        """
        The ROM is a frequency domain waveform, so easier.
        """
        # get hptilde
        approx_enum = lalsim.GetApproximantFromString('SEOBNRv2_ROM_DoubleSpin')
        flags = lalsim.SimInspiralCreateWaveformFlags()
        htilde, _ = lalsim.SimInspiralChooseFDWaveform(0, df, self.m1 * MSUN_SI,
            self.m2 * MSUN_SI, 0., 0., self.spin1z, 0., 0., self.spin2z,
            self.bank.flow, f_final, self.bank.flow, 1e6*PC_SI, 0., 0., 0.,
            flags, None, 1, 8, approx_enum)
        return htilde


class EOBNRv2Template(SEOBNRv2Template):

    param_names = ("m1", "m2")
    param_formats = ("%.2f", "%.2f")

    def __init__(self, m1, m2, bank):
        # Use everything from SEOBNRv2Template class except for the
        # _compute_waveform method; call parent __init__ with spins
        # set to zero
        SEOBNRv2Template.__init__(self, m1, m2, 0, 0, bank)

    def _compute_waveform(self, df, f_final):
        """
        Since EOBNRv2 is a time domain waveform, we have to generate it,
        then FFT to frequency domain.
        """
        # need to compute dt from df, duration
        sample_rate = 2**np.ceil(np.log2(2*f_final))
        dt = 1. / sample_rate
        # get hplus
        hplus, hcross = lalsim.SimIMREOBNRv2DominantMode(
            0., dt, self.m1 * MSUN_SI, self.m2 * MSUN_SI,
            self.bank.flow, 1e6*PC_SI, 0.)
        # zero-pad up to 1/df
        N = int(sample_rate / df)
        hplus = lal.ResizeREAL8TimeSeries(hplus, 0, N)
        # taper
        lalsim.SimInspiralREAL8WaveTaper(hplus.data, lalsim.SIM_INSPIRAL_TAPER_START)

        # create vector to hold output and plan
        htilde = lal.CreateCOMPLEX16FrequencySeries("h(f)", hplus.epoch, hplus.f0, df, lal.HertzUnit, int(N/2 + 1))
        fftplan = lal.CreateForwardREAL8FFTPlan(N, 0)

        # do the fft
        lal.REAL8TimeFreqFFT(htilde, hplus, fftplan)

        return htilde


class TaylorF2RedSpinTemplate(AlignedSpinTemplate):

    param_names = ("m1", "m2", "chired")
    param_formats = ("%.5f", "%.5f", "%+.4f")

    __slots__ = ("chired", "_dur", "_mchirp", "_tau0", "_eta", "_theta0", "_theta3", "_theta3s")

    def __init__(self, m1, m2, spin1z, spin2z, bank):

        AlignedSpinTemplate.__init__(self, m1, m2, spin1z, spin2z, bank)
        self.chired = lalsim.SimInspiralTaylorF2ReducedSpinComputeChi(m1, m2, spin1z, spin2z)
        self._dur = self._get_dur()
        self._eta = m1*m2/(m1+m2)**2
        self._theta0, self._theta3, self._theta3s = compute_chirptimes(self._mchirp, self._eta, self.chired, self.bank.flow)

    def _get_f_final(self):
        return 6**-1.5 / (PI * (self.m1 + self.m2) * MTSUN_SI)  # ISCO

    def _get_dur(self):
        return lalsim.SimInspiralTaylorF2ReducedSpinChirpTime(self.bank.flow,
            self.m1 * MSUN_SI, self.m2 * MSUN_SI, self.chired,
            7)

    def finalize_as_template(self):
        if not self.bank.use_metric: return

        df, PSD = get_neighborhood_PSD([self], self.bank.flow, self.bank.noise_model)

        if df not in self.bank._moments or len(PSD) - self.bank.flow // df > self.bank._moments[df][0].length:
            real8vector_psd = CreateREAL8Vector(len(PSD))
            real8vector_psd.data[:] = PSD
            self.bank._moments[df] = create_moments(df, self.bank.flow, len(PSD))
            SimInspiralTaylorF2RedSpinComputeNoiseMoments(*(self.bank._moments[df] + (real8vector_psd, self.bank.flow, df)))

        self._metric = SimInspiralTaylorF2RedSpinMetricChirpTimes(self._theta0, self._theta3, self._theta3s, self.bank.flow, df, *self.bank._moments[df])
        if isnan(self._metric[0]):
            raise ValueError("g00 is nan")

    def _compute_waveform(self, df, f_final):

        wf = lalsim.SimInspiralTaylorF2ReducedSpin(
            0, df, self.m1 * MSUN_SI, self.m2 * MSUN_SI, self.chired,
            self.bank.flow, 0, 1000000 * PC_SI, 7, 7)
        # have to resize wf to next pow 2 for FFT plan caching
        wf = lal.ResizeCOMPLEX16FrequencySeries( wf, 0, ceil_pow_2(wf.data.length) )
        return wf

    def metric_match(self, other, df, **kwargs):
        g00, g01, g02, g11, g12, g22 = self._metric
        dx0 = other._theta0 - self._theta0
        dx1 = other._theta3 - self._theta3
        dx2 = other._theta3s - self._theta3s
        match = 1 - (
            g00 * dx0 * dx0
          + 2 * g01 * dx0 * dx1
          + 2 * g02 * dx0 * dx2
          + g11 * dx1 * dx1
          + 2 * g12 * dx1 * dx2
          + g22 * dx2 * dx2
        )

        return match


class PrecessingTemplate(Template):
    """
    A generic class for precessing templates. These models require the
    full fifteen-dimensional parameter space to specify the observed
    signal in the detector.
    """
    param_names = ("m1", "m2", "spin1x", "spin1y", "spin1z", "spin2x", "spin2y", "spin2z", "theta", "phi", "iota", "psi")
    param_formats = ("%.2f","%.2f","%.2f","%.2f","%.2f","%.2f","%.2f","%.2f","%.2f","%.2f","%.2f","%.2f")
    __slots__ = param_names + ("bank", "_dur","_mchirp","_tau0")

    def __init__(self, m1, m2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, theta, phi, iota, psi, bank):

        Template.__init__(self, m1, m2, bank)
        self.m1 = m1
        self.m2 = m2
        self.spin1x = spin1x
        self.spin1y = spin1y
        self.spin1z = spin1z
        self.spin2x = spin2x
        self.spin2y = spin2y
        self.spin2z = spin2z
        self.theta = theta
        self.phi = phi
        self.iota = iota
        self.psi = psi
        self.bank = bank

        # derived quantities
        self._mchirp = compute_mchirp(m1, m2)

        # FIXME: What is an appropriate f_final for generic waveforms?
        self.chieff,_ = SimIMRPhenomPCalculateModelParameters(self.m1, self.m2,
                           self.bank.flow, np.sin(self.iota), float(0),
                           np.cos(self.iota), float(self.spin1x),
                           float(self.spin1y), float(self.spin1z),
                           float(self.spin2x), float(self.spin2y),
                           float(self.spin2z))[:2]

        # FIXME: Is this appropriate?
        self._dur = lalsim.SimInspiralTaylorF2ReducedSpinChirpTime(\
                       self.bank.flow, self.m1 * MSUN_SI, self.m2 * MSUN_SI,
                       self.chieff, 7) + 1000 * (self.m1 + self.m2) * MTSUN_SI

    def _get_f_final(self):
        return spawaveform.imrffinal(self.m1, self.m2, self.chieff)

    @classmethod
    def from_sim(cls, sim, bank):
        # theta = polar angle wrt overhead
        #       = pi/2 - latitude (which is 0 on the horizon)
        return cls(sim.mass1, sim.mass2, sim.spin1x, sim.spin1y, sim.spin1z, sim.spin2x, sim.spin2y, sim.spin2z, np.pi/2 - sim.latitude, sim.longitude, sim.inclination, sim.polarization, bank)


class IMRPhenomPTemplate(PrecessingTemplate):
    """
    IMRPhenomP precessing IMR model.
    """
    __slots__ = PrecessingTemplate.param_names + ("bank", "chieff", "chipre", "_dur", "_mchirp", "_tau0")
    param_names = ("m1", "m2", "chieff", "chipre")
    param_formats = ("%.2f", "%.2f", "%.2f", "%.2f")

    def __init__(self, m1, m2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, theta, phi, iota, psi, bank):

        PrecessingTemplate.__init__(self, m1, m2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, theta, phi, iota, psi, bank)
        # derived quantities
        self.chieff, self.chipre = SimIMRPhenomPCalculateModelParameters(self.m1, self.m2, self.bank.flow, np.sin(self.iota), float(0), np.cos(self.iota), float(self.spin1x), float(self.spin1y), float(self.spin1z), float(self.spin2x), float(self.spin2y), float(self.spin2z))[:2] # FIXME are the other four parameters informative?

        self._dur = self._imrdur()
        # FIXME: is this ffinal and dur appropriate for PhenomP?

    def _get_f_final(self):
        return spawaveform.imrffinal(self.m1, self.m2, self.chieff)

    def _imrdur(self):
        """
        Ajith gave us the heuristic that chirp + 1000 M is a conservative
        estimate for the length of a full IMR waveform.
        """
        return lalsim.SimInspiralTaylorF2ReducedSpinChirpTime(self.bank.flow,
            self.m1 * MSUN_SI, self.m2 * MSUN_SI, self.chieff,
            7) + 1000 * (self.m1 + self.m2) * MTSUN_SI

    def _compute_waveform(self, df, f_final):

        approx = lalsim.GetApproximantFromString( "IMRPhenomP" )
        phi0 = 0  # what is phi0?
        lmbda1 = lmbda2 = 0
        ampO = 3
        phaseO = 7 # are these PN orders correct for PhenomP?
        hplus_fd, hcross_fd = lalsim.SimInspiralChooseFDWaveform(
            phi0, df,
            self.m1*MSUN_SI, self.m2*MSUN_SI,
            self.spin1x, self.spin1y, self.spin1z, self.spin2x, self.spin2y, self.spin2z,
            self.bank.flow, f_final,
            40.0, # reference frequency, want it to be fixed to a constant value always and forever
            1e6*PC_SI, # irrelevant parameter for banks/banksims
            self.iota,
            lmbda1, lmbda2, # irrelevant parameters for BBH
            None, None, # non-GR parameters
            ampO, phaseO, approx)

        # project onto detector
        return project_hplus_hcross(hplus_fd, hcross_fd, self.theta, self.phi, self.psi)


class SpinTaylorT4Template(Template):

    param_names = ("m1","m2","s1x","s1y","s1z","s2x","s2y","s2z","inclination","theta","phi","psi")
    param_formats = ("%.2f","%.2f","%.2f","%.2f","%.2f","%.2f","%.2f","%.2f","%.2f","%.2f","%.2f","%.2f")

    __slots__ = ("m1","m2","s1x","s1y","s1z","s2x","s2y","s2z","inclination","theta","phi","psi","bank", "_dur","_mchirp", "_tau0")

    def __init__(self,m1,m2,s1x,s1y,s1z,s2x,s2y,s2z,inclination,theta,phi,psi,bank):
        Template.__init__(self, m1, m2, bank)
        self.m1 = float(m1)
        self.m2 = float(m2)
        self.s1x = float(s1x)
        self.s1y = float(s1y)
        self.s1z = float(s1z)
        self.s2x = float(s2x)
        self.s2y = float(s2y)
        self.s2z = float(s2z)
        self.inclination=float(inclination)
        self.theta = float(theta)
        self.phi = float(phi)
        self.psi = float(psi)
        self.bank = bank

        # derived quantities
        self._dur = lalsim.SimInspiralTaylorF2ReducedSpinChirpTime(bank.flow, m1 * MSUN_SI, m2 * MSUN_SI, lalsim.SimInspiralTaylorF2ReducedSpinComputeChi(self.m1*MSUN_SI,self.m2*MSUN_SI,self.s1z,self.s2z), 7)
        self._mchirp = compute_mchirp(m1, m2)

    def _get_f_final(self):
        return 6**-1.5 / (PI * (self.m1 + self.m2) * MTSUN_SI)

    def _compute_waveform(self, df, f_final):
        # Time domain, so compute then FFT
        # need to compute dt from df, duration
        sample_rate = 2**np.ceil(np.log2(2*f_final))
        dt = 1. / sample_rate

        # Generate waveform in time domain
        hplus, hcross = lalsim.SimInspiralSpinTaylorT4(
            0,				# GW phase at reference freq (rad)
            1,				# tail gauge term (default = 1)
            dt,				# sampling interval (s)
            self.m1 * lal.MSUN_SI,	# mass of companion 1 (kg)
            self.m2 * lal.MSUN_SI,	# mass of companion 2 (kg)
            self.bank.flow,			# start GW frequency (Hz)
            self.bank.flow,			# reference GW frequency at which phase is set (Hz)
            1e6*PC_SI,			# distance of source (m)
            self.s1x,			# initial value of S1x
            self.s1y,			# initial value of S1y
            self.s1z,			# initial value of S1z
            self.s2x,			# initial value of S2x
            self.s2y,			# initial value of S2y
            self.s2z,			# initial value of S2z
            np.sin(self.inclination),	# initial value of LNhatx
            0,				# initial value of LNhaty
            np.cos(self.inclination),	# initial value of LNhatz
            np.cos(self.inclination),	# initial value of E1x
            0,				# initial value of E1y
            -np.sin(self.inclination),	# initial value of E1z
            0,				# tidal deformability of mass 1
            0,				# tidal deformability of mass 2
            1, # phenom. parameter describing induced quad. moment of body 1 (=1 for BHs, ~2-12 for NSs)
            1, # phenom. parameter describing induced quad. moment of body 2 (=1 for BHs, ~2-12 for NSs)
            7, # twice PN spin order
            0, # twice PN tidal order
            7, # twice PN phase order
            0 # twice PN amplitude order
        )

        # project onto detector
        hoft = project_hplus_hcross(hplus, hcross, self.theta, self.phi, self.psi)

        # zero-pad up to 1/df
        N = ceil_pow_2(int(sample_rate / df))
        hoft = lal.ResizeREAL8TimeSeries(hoft, 0, N)

        # taper
        lalsim.SimInspiralREAL8WaveTaper(hoft.data, lalsim.SIM_INSPIRAL_TAPER_STARTEND)

        # create vector to hold output and plan
        htilde = lal.CreateCOMPLEX16FrequencySeries("h(f)", hoft.epoch, hoft.f0, df, lal.HertzUnit, int(N/2 + 1))
        fftplan = lal.CreateForwardREAL8FFTPlan(N, 0)

        # do the fft
        lal.REAL8TimeFreqFFT(htilde, hoft, fftplan)

        return htilde

    @classmethod
    def from_sim(cls, sim, bank):
        # theta = polar angle wrt overhead
        #       = pi/2 - latitude (which is 0 on the horizon)
        return cls(sim.mass1, sim.mass2, sim.spin1x, sim.spin1y, sim.spin1z, sim.spin2x, sim.spin2y, sim.spin2z, sim.inclination, np.pi/2 - sim.latitude, sim.longitude, sim.polarization, bank)


class SpinTaylorT5Template(Template):
    param_names = ("m1","m2","s1x","s1y","s1z","s2x","s2y","s2z","inclination","theta","phi","psi")
    param_formats = ("%.2f","%.2f","%.2f","%.2f","%.2f","%.2f","%.2f","%.2f","%.2f","%.2f","%.2f","%.2f")

    __slots__ = ("m1","m2","s1x","s1y","s1z","s2x","s2y","s2z","inclination","theta","phi","psi","bank", "_dur","_mchirp", "_tau0")

    def __init__(self,m1,m2,s1x,s1y,s1z,s2x,s2y,s2z,inclination,theta,phi,psi,bank):
        Template.__init__(self, m1, m2, bank)
        self.m1 = float(m1)
        self.m2 = float(m2)
        self.s1x = float(s1x)
        self.s1y = float(s1y)
        self.s1z = float(s1z)
        self.s2x = float(s2x)
        self.s2y = float(s2y)
        self.s2z = float(s2z)
        self.inclination=float(inclination)
        self.theta = float(theta)
        self.phi = float(phi)
        self.psi = float(psi)
        self.bank = bank

        # derived quantities
        self._dur = lalsim.SimInspiralTaylorF2ReducedSpinChirpTime(bank.flow, m1 * MSUN_SI, m2 * MSUN_SI, lalsim.SimInspiralTaylorF2ReducedSpinComputeChi(self.m1*MSUN_SI,self.m2*MSUN_SI,self.s1z,self.s2z), 7)
        self._mchirp = compute_mchirp(m1, m2)

    def _get_f_final(self):
        return 6**-1.5 / (PI * (self.m1 + self.m2) * MTSUN_SI)

    def _compute_waveform(self, df, f_final):
        # Time domain, so compute then FFT
        # need to compute dt from df, duration
        sample_rate = 2**np.ceil(np.log2(2*f_final))
        dt = 1. / sample_rate
        # Generate waveform in time domain
        hplus, hcross = lalsim.SimInspiralSpinTaylorT5(
            0,				# GW phase at reference freq (rad)
            dt,				# sampling interval (s)
            self.m1 * lal.MSUN_SI,	# mass of companion 1 (kg)
            self.m2 * lal.MSUN_SI,	# mass of companion 2 (kg)
            self.bank.flow,			# start GW frequency (Hz)
            1e6*PC_SI,			# distance of source (m)
            self.s1x,			# initial value of S1x
            self.s1y,			# initial value of S1y
            self.s1z,			# initial value of S1z
            self.s2x,			# initial value of S2x
            self.s2y,			# initial value of S2y
            self.s2z,			# initial value of S2z
            self.inclination,		# inclination angle - careful with definition (line of sight to total vs orbital angular momentum)
            7,				# twice PN phase order
            0
        )

        # project onto detector
        hoft = project_hplus_hcross(hplus, hcross, self.theta, self.phi, self.psi)

        # zero-pad up to 1/df
        N = int(sample_rate / df)
        hoft = lal.ResizeREAL8TimeSeries(hoft, 0, N)

        # taper
        lalsim.SimInspiralREAL8WaveTaper(hoft.data, lalsim.SIM_INSPIRAL_TAPER_STARTEND)

        # create vector to hold output and plan
        htilde = lal.CreateCOMPLEX16FrequencySeries("h(f)", hoft.epoch, hoft.f0, df, lal.HertzUnit, int(N/2 + 1))
        fftplan = lal.CreateForwardREAL8FFTPlan(N, 0)

        # do the fft
        lal.REAL8TimeFreqFFT(htilde, hoft, fftplan)

        return htilde

    @classmethod
    def from_sim(cls, sim, bank):
        # theta = polar angle wrt overhead
        #       = pi/2 - latitude (which is 0 on the horizon)
        return cls(sim.mass1, sim.mass2, sim.spin1x, sim.spin1y, sim.spin1z, sim.spin2x, sim.spin2y, sim.spin2z, sim.inclination, np.pi/2 - sim.latitude, sim.longitude, sim.polarization, bank)


waveforms = {
    "TaylorF2RedSpin": TaylorF2RedSpinTemplate,
    "IMRPhenomB": IMRPhenomBTemplate,
    "IMRPhenomC": IMRPhenomCTemplate,
    "IMRPhenomP": IMRPhenomPTemplate,
    "SEOBNRv2": SEOBNRv2Template,
    "SEOBNRv2_ROM_DoubleSpin": SEOBNRv2ROMDoubleSpinTemplate,
    "EOBNRv2": EOBNRv2Template,
    "SpinTaylorT4": SpinTaylorT4Template,
    "SpinTaylorT5": SpinTaylorT5Template,
}
