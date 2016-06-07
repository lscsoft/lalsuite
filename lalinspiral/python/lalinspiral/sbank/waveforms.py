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
from lalinspiral import InspiralSBankComputeMatch, InspiralSBankComputeMatchMaxSkyLoc
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

def compute_correlation(htilde1, htilde2, deltaF):
    """
    Find the real component of correlation between htilde1 and htilde2.
    """
    # vdot is dot with complex conjugation
    return float(np.vdot(htilde1,htilde2).real * 4 * deltaF)

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


class AlignedSpinTemplate(object):
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
    approximant = None
    __slots__ = ("m1", "m2", "spin1z", "spin2z", "chieff", "bank", "_dur", "_mchirp", "_tau0", "_wf", "_metric", "sigmasq", "is_seed_point", "_f_final", "_fhigh_max")
    param_names = ("m1", "m2", "spin1z", "spin2z")
    param_formats = ("%.2f", "%.2f", "%.2f", "%.2f")

    def __init__(self, m1, m2, spin1z, spin2z, bank):

        self.m1 = float(m1)
        self.m2 = float(m2)
        self.spin1z = float(spin1z)
        self.spin2z = float(spin2z)
        self.chieff = lalsim.SimIMRPhenomBComputeChi(m1, m2, spin1z, spin2z)
        self.bank = bank

        self._wf = {}
        self._metric = None
        self.sigmasq = 0.
        self._mchirp = compute_mchirp(m1, m2)
        self._tau0 = compute_tau0( self._mchirp, bank.flow)
        self._dur = self._get_dur()
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

    def _get_isco_f_final(self):
        return 6**-1.5 / (PI * (self.m1 + self.m2) * MTSUN_SI)  # ISCO

    def _get_imr_f_final(self):
        return lalsim.SimInspiralGetFrequency(self.m1*lal.MSUN_SI,
                                self.m2*lal.MSUN_SI, 0., 0., self.spin1z, 0, 0,
                                self.spin2z, lalsim.fSEOBNRv2RD)

    def _get_chirp_dur(self):
        return lalsim.SimInspiralTaylorF2ReducedSpinChirpTime(self.bank.flow,
            self.m1 * MSUN_SI, self.m2 * MSUN_SI, self.chieff, 7)

    def _get_imr_dur(self):
        """
        Following Jolien's suggestion for the duration of an IMR
        waveform, we compute the chirp time, scale by 10% and pad by
        one second.
        """
        # FIXME: Is a minimal time of 1.0s too long?
        return 1.1 * self._get_chirp_dur() + 1.0

    def _get_f_final(self):
        err_msg = "Template classes must provide a function _get_f_final "
        err_msg += "to compute the appropriate final frequency."
        raise NotImplementedError(err_msg)

    def _get_dur(self):
        err_msg = "Template classes must provide a function _get_dur "
        err_msg += "to compute template neighborhoods."
        raise NotImplementedError(err_msg)

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

    def __repr__(self):
        return "(%s)" % ", ".join(self.param_formats) % self.params

    def finalize_as_template(self):
        """
        This template is being added to a bank. Compute metric coefficients
        or do whatever other expensive things you didn't want to do when
        it was just a proposal.
        """
        pass


    def _compute_waveform(self, df, f_final):

        flags = lalsim.SimInspiralCreateWaveformFlags()
        phi0 = 0  # This is a reference phase, and not an intrinsic parameter
        lmbda1 = lmbda2 = 0 # No tidal terms here
        ampO = -1 # Are these the correct values??
        phaseO = -1 # Are these the correct values??
        approx = lalsim.GetApproximantFromString( self.approximant )

        if lalsim.SimInspiralImplementedFDApproximants(approx):
            hplus_fd, hcross_fd = lalsim.SimInspiralChooseFDWaveform(
                0, df,
                self.m1 * MSUN_SI, self.m2 * MSUN_SI,
                0., 0., self.spin1z, 0., 0., self.spin2z,
                self.bank.flow, f_final, self.bank.flow,
                1e6*PC_SI, 0., 0., 0., flags, None, ampO, phaseO, approx)

        else:
            hplus_fd, hcross_fd = lalsim.SimInspiralFD(
                phi0, df, self.m1*MSUN_SI, self.m2*MSUN_SI, 0,
                0, self.spin1z, 0, 0, self.spin2z,
                self.bank.flow, f_final, 40.0, 1e6*PC_SI, 0, 0,
                lmbda1, lmbda2, # irrelevant parameters for BBH
                None, None, # non-GR parameters
                ampO, phaseO, approx)

        return hplus_fd


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
                ASD2 = np.ones(wf.data.length) * np.inf
                ASD2[:len(ASD)] = ASD
                ASD = ASD2
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
        return InspiralSBankComputeMatch(self.get_whitened_normalized(df, **kwargs), other.get_whitened_normalized(df, **kwargs), workspace_cache[0])

    def clear(self):
        self._wf = {}


class IMRAlignedSpinTemplate(AlignedSpinTemplate):
    """
    A convenience class for IMR aligned spin templates. Specific
    implementations of IMR aligned spin templates should sub-class
    this class.
    """
    def _get_dur(self):
        return self._get_imr_dur()

    def _get_f_final(self):
        return self._get_imr_f_final()


class InspiralAlignedSpinTemplate(AlignedSpinTemplate):
    """
    A convenience class for inspiral-only aligned spin
    templates. Specific implementations of inspiral-only aligned spin
    templates should sub-class this class.
    """
    __slots__ = ("chired")
    def __init__(self, m1, m2, spin1z, spin2z, bank):

        self.chired = lalsim.SimInspiralTaylorF2ReducedSpinComputeChi(m1, m2, spin1z, spin2z)
        AlignedSpinTemplate.__init__(self, m1, m2, spin1z, spin2z, bank)

    def _get_dur(self):
        return self._get_chirp_dur()

    def _get_f_final(self):
        return self._get_isco_f_final()

#
# IMRPhenom*Template(IMRAlignedSpinTemplate)
#
class IMRPhenomBTemplate(IMRAlignedSpinTemplate):
    approximant = "IMRPhenomB"
    param_names = ("m1", "m2", "chieff")
    param_formats = ("%.2f", "%.2f", "%+.2f")
    def _compute_waveform(self, df, f_final):
        return lalsim.SimIMRPhenomBGenerateFD(0, df,
            self.m1 * MSUN_SI, self.m2 * MSUN_SI,
            self.chieff, self.bank.flow, f_final, 1000000 * PC_SI)


class IMRPhenomCTemplate(IMRPhenomBTemplate):
    approximant = "IMRPhenomC"
    def _compute_waveform(self, df, f_final):
        return lalsim.SimIMRPhenomCGenerateFD(
            0, df,
            self.m1 * MSUN_SI, self.m2 * MSUN_SI,
            self.chieff, self.bank.flow, f_final, 1000000 * PC_SI)

class IMRPhenomDTemplate(IMRPhenomBTemplate):
    approximant = "IMRPhenomD"
    def _compute_waveform(self, df, f_final):
        return lalsim.SimIMRPhenomDGenerateFD(
            0, 0, df, # ref phase, ref frequency, df
            self.m1 * MSUN_SI, self.m2 * MSUN_SI,
            self.spin1z, self.spin2z,
            self.bank.flow, f_final, 1000000 * PC_SI)


class SEOBNRv2Template(IMRAlignedSpinTemplate):
    approximant = "SEOBNRv2"

class SEOBNRv2ROMDoubleSpinTemplate(SEOBNRv2Template):
    approximant = "SEOBNRv2_ROM_DoubleSpin"

class SEOBNRv2ROMDoubleSpinHITemplate(SEOBNRv2Template):
    approximant = "SEOBNRv2_ROM_DoubleSpin_HI"

class EOBNRv2Template(SEOBNRv2Template):
    approximant = "EOBNRv2"
    param_names = ("m1", "m2")
    param_formats = ("%.2f", "%.2f")

    def __init__(self, m1, m2, bank):
        # Use everything from SEOBNRv2Template class except call
        # parent __init__ with spins set to zero
        SEOBNRv2Template.__init__(self, m1, m2, 0, 0, bank)

class TaylorF2RedSpinTemplate(InspiralAlignedSpinTemplate):
    approximant = "TaylorF2RedSpin"
    param_names = ("m1", "m2", "chired")
    param_formats = ("%.5f", "%.5f", "%+.4f")

    __slots__ = ("chired", "_dur", "_mchirp", "_tau0", "_eta", "_theta0", "_theta3", "_theta3s")

    def __init__(self, m1, m2, spin1z, spin2z, bank):

        AlignedSpinTemplate.__init__(self, m1, m2, spin1z, spin2z, bank)
        self.chired = lalsim.SimInspiralTaylorF2ReducedSpinComputeChi(m1, m2, spin1z, spin2z)
        self._dur = self._get_dur()
        self._eta = m1*m2/(m1+m2)**2
        self._theta0, self._theta3, self._theta3s = compute_chirptimes(self._mchirp, self._eta, self.chired, self.bank.flow)

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


class TaylorF2Template(InspiralAlignedSpinTemplate):
    approx_name = "TaylorF2"

    def _compute_waveform(self, df, f_final):
        phi0 = 0  # This is a reference phase, and not an intrinsic parameter
        lmbda1 = lmbda2 = 0 # No tidal terms here
        ampO = 0
        phaseO = 7
        approx = lalsim.GetApproximantFromString( self.approx_name )
        wave_flags = lalsim.SimInspiralCreateWaveformFlags()
        lalsim.SimInspiralSetSpinOrder(wave_flags, 5)

        hplus_fd, hcross_fd = lalsim.SimInspiralChooseFDWaveform(
                0, df, self.m1*MSUN_SI, self.m2*MSUN_SI, 0, 0, self.spin1z,
                0, 0, self.spin2z, self.bank.flow, f_final, 40.0, 1e6*PC_SI, 0,
                lmbda1, lmbda2, # irrelevant parameters for BBH
                wave_flags, None, # non-GR parameters
                ampO, phaseO, approx)

        # Must set values greater than _get_f_final to 0
        act_f_max = self._get_f_final()
        f_max_idx = int(act_f_max / df + 0.999)
        hplus_fd.data.data[f_max_idx:] = 0

        return hplus_fd


class PrecessingSpinTemplate(AlignedSpinTemplate):
    """
    A generic class for precessing templates. These models require the
    full fifteen-dimensional parameter space to specify the observed
    signal in the detector.
    """
    param_names = ("m1", "m2", "spin1x", "spin1y", "spin1z", "spin2x", "spin2y", "spin2z", "theta", "phi", "iota", "psi")
    param_formats = ("%.2f","%.2f","%.2f","%.2f","%.2f","%.2f","%.2f","%.2f","%.2f","%.2f","%.2f","%.2f")
    __slots__ = param_names + ("bank", "chieff", "chipre", "_dur","_mchirp","_tau0", "_wf_hp", "_wf_hc", "_hpsigmasq", "_hcsigmasq", "_hphccorr")

    def __init__(self, m1, m2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, theta, phi, iota, psi, bank):

        AlignedSpinTemplate.__init__(self, m1, m2, spin1z, spin2z, bank)
        self.spin1x = spin1x
        self.spin1y = spin1y
        self.spin2x = spin2x
        self.spin2y = spin2y

        self.theta = theta
        self.phi = phi
        self.iota = iota
        self.psi = psi

        self.chieff, self.chipre = SimIMRPhenomPCalculateModelParameters(self.m1, self.m2, self.bank.flow, np.sin(iota), float(0), np.cos(iota), spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, 1)[:2]

        self._wf = {}
        self._metric = None
        self.sigmasq = 0.
        self._wf_hp = {}
        self._wf_hc = {}
        self._hpsigmasq = {}
        self._hcsigmasq = {}
        self._hphccorr = {}

    @classmethod
    def from_sim(cls, sim, bank):
        # theta = polar angle wrt overhead
        #       = pi/2 - latitude (which is 0 on the horizon)
        return cls(sim.mass1, sim.mass2, sim.spin1x, sim.spin1y, sim.spin1z, sim.spin2x, sim.spin2y, sim.spin2z, np.pi/2 - sim.latitude, sim.longitude, sim.inclination, sim.polarization, bank)

    def _compute_waveform_comps(self, df, f_final):
        approx = lalsim.GetApproximantFromString( self.approximant )
        phi0 = 0  # This is a reference phase, and not an intrinsic parameter
        lmbda1 = lmbda2 = 0 # No tidal terms here
        ampO = -1 # Are these the correct values??
        phaseO = -1 # Are these the correct values??

        if lalsim.SimInspiralImplementedFDApproximants(approx):
            hplus_fd, hcross_fd = lalsim.SimInspiralChooseFDWaveform(
                phi0, df, self.m1*MSUN_SI, self.m2*MSUN_SI, self.spin1x,
                self.spin1y, self.spin1z, self.spin2x, self.spin2y, self.spin2z,
                self.bank.flow, f_final, self.bank.flow, 1e6*PC_SI, self.iota,
                lmbda1, lmbda2, # irrelevant parameters for BBH
                None, None, # non-GR parameters
                ampO, phaseO, approx)
        else:
            hplus_fd, hcross_fd = lalsim.SimInspiralFD(
                phi0, df, self.m1*MSUN_SI, self.m2*MSUN_SI,
                self.spin1x, self.spin1y, self.spin1z,
                self.spin2x, self.spin2y, self.spin2z,
                self.bank.flow, f_final, self.bank.flow, 1e6*PC_SI, 0, self.iota,
                lmbda1, lmbda2, # irrelevant parameters for BBH
                None, None, # non-GR parameters
                ampO, phaseO, approx)

        return hplus_fd, hcross_fd

    def _compute_waveform(self, df, f_final):
        hplus_fd, hcross_fd = self._compute_waveform_comps(df, f_final)

        # project onto detector
        return project_hplus_hcross(hplus_fd, hcross_fd, self.theta, self.phi,
                                    self.psi)

    def get_whitened_normalized_comps(self, df, ASD=None, PSD=None):
        """
        Return a COMPLEX8FrequencySeries of h+ and hx, whitened by the
        given ASD and normalized. The waveform is not zero-padded to
        match the length of the ASD, so its normalization depends on
        its own length.
        """
        if not self._wf_hp.has_key(df):
            # Clear self._wf as it won't be needed any more if calling here
            self._wf = {}
            # Generate a new wf
            hp, hc = self._compute_waveform_comps(df, self.f_final)
            if ASD is None:
                ASD = PSD**0.5
            if hp.data.length > len(ASD):
                err_msg = "waveform has length greater than ASD; cannot whiten"
                raise ValueError(err_msg)
            arr_view_hp = hp.data.data
            arr_view_hc = hc.data.data

            # Whiten
            arr_view_hp[:] /= ASD[:hp.data.length]
            arr_view_hp[:int(self.bank.flow / df)] = 0.
            arr_view_hp[int(self.f_final/df) : hp.data.length] = 0.

            arr_view_hc[:] /= ASD[:hc.data.length]
            arr_view_hc[:int(self.bank.flow / df)] = 0.
            arr_view_hc[int(self.f_final/df) : hc.data.length] = 0.

            # Get normalization factors
            self._hpsigmasq[df] = compute_sigmasq(arr_view_hp, df)
            self._hcsigmasq[df] = compute_sigmasq(arr_view_hc, df)
            self._hphccorr[df] = compute_correlation(arr_view_hp, arr_view_hc, df)

            self._wf_hp[df] = FrequencySeries_to_COMPLEX8FrequencySeries(hp)
            self._wf_hc[df] = FrequencySeries_to_COMPLEX8FrequencySeries(hc)


        return self._wf_hp[df], self._wf_hc[df], self._hpsigmasq[df],\
               self._hcsigmasq[df], self._hphccorr[df]

    def brute_match(self, other, df, workspace_cache, **kwargs):

        # Template generates hp and hc
        hp, hc, hpsigmasq, hcsigmasq, hphccorr =  self.get_whitened_normalized_comps(df, **kwargs)

        # Proposal generates h(t), sky loc is later discarded.
        proposal = other.get_whitened_normalized(df, **kwargs)

        # maximize over sky position of template
        return InspiralSBankComputeMatchMaxSkyLoc(hp, hc, hpsigmasq,
                                        hcsigmasq, hphccorr, proposal,
                                        workspace_cache[0], workspace_cache[1])

    @classmethod
    def from_sngl(cls, sngl, bank):
        # FIXME: Using alpha columns to hold theta, phi, iota, psi
        return cls(sngl.mass1, sngl.mass2, sngl.spin1x, sngl.spin1y,
                   sngl.spin1z, sngl.spin2x, sngl.spin2y, sngl.spin2z,
                   sngl.alpha1, sngl.alpha2, sngl.alpha3, sngl.alpha4, bank)

    def to_sngl(self):
        # note that we use the C version; this causes all numerical values to
        # be initiated as 0 and all strings to be '', which is nice
        row = SnglInspiralTable()
        row.mass1 = self.m1
        row.mass2 = self.m2
        row.mtotal = self.m1 + self.m2
        row.mchirp = self._mchirp
        row.eta = row.mass1 * row.mass2 / (row.mtotal * row.mtotal)
        row.tau0, row.tau3 = m1m2_to_tau0tau3(self.m1, self.m2, self.bank.flow)
        row.f_final = self.f_final
        row.template_duration = self._dur
        row.spin1x = self.spin1x
        row.spin1y = self.spin1y
        row.spin1z = self.spin1z
        row.spin2x = self.spin2x
        row.spin2y = self.spin2y
        row.spin2z = self.spin2z
        row.alpha1 = self.theta
        row.alpha2 = self.phi
        row.alpha3 = self.iota
        row.alpha4 = self.psi
        row.sigmasq = self.sigmasq
        return row


class IMRPrecessingSpinTemplate(PrecessingSpinTemplate):
    """
    A convenience class for IMR precessing spin templates. Specific
    implementations of IMR precessing spin templates should sub-class
    this class.
    """
    def _get_dur(self):
        return self._get_imr_dur()

    def _get_f_final(self):
        return self._get_imr_f_final()


class InspiralPrecessingSpinTemplate(PrecessingSpinTemplate):
    """
    A convenience class for inpsiral-only precessing spin templates. Specific
    implementations of inspiral-only precessing spin templates should sub-class
    this class.
    """
    def _get_dur(self):
        return self._get_chirp_dur()

    def _get_f_final(self):
        return self._get_isco_f_final()


class SpinTaylorF2Template(InspiralPrecessingSpinTemplate):
    approximant = "SpinTaylorF2"
    def __init__(self, m1, m2, spin1x, spin1y, spin1z,
                 theta, phi, iota, psi, bank):
        super(SpinTaylorF2Template,self).__init__(m1, m2,
                                    spin1x, spin1y, spin1z, 0, 0, 0,
                                    theta, phi, iota, psi, bank)

    def _compute_waveform_comps(self, df, f_final):
        hplus_fd, hcross_fd = \
          super(SpinTaylorF2Template,self)._compute_waveform_comps(df, f_final)
        # Must set values greater than _get_f_final to 0
        act_f_max = self._get_f_final()
        f_max_idx = int(act_f_max / df + 0.999)
        hplus_fd.data.data[f_max_idx:] = 0
        hcross_fd.data.data[f_max_idx:] = 0

        return hplus_fd, hcross_fd

    @classmethod
    def from_sngl(cls, sngl, bank):
        # FIXME: Using alpha columns to hold theta, phi, iota, psi
        assert sngl.spin2x == sngl.spin2y == sngl.spin2z == 0
        return cls(sngl.mass1, sngl.mass2, sngl.spin1x, sngl.spin1y,
                   sngl.spin1z, sngl.alpha1, sngl.alpha2, sngl.alpha3, sngl.alpha4, bank)

    def to_sngl(self):
        # note that we use the C version; this causes all numerical values to
        # be initiated as 0 and all strings to be '', which is nice
        row = SnglInspiralTable()
        row.mass1 = self.m1
        row.mass2 = self.m2
        row.mtotal = self.m1 + self.m2
        row.mchirp = self._mchirp
        row.eta = row.mass1 * row.mass2 / (row.mtotal * row.mtotal)
        row.tau0, row.tau3 = m1m2_to_tau0tau3(self.m1, self.m2, self.bank.flow)
        row.f_final = self.f_final
        row.template_duration = self._dur
        row.spin1x = self.spin1x
        row.spin1y = self.spin1y
        row.spin1z = self.spin1z
        row.spin2x = 0
        row.spin2y = 0
        row.spin2z = 0
        row.alpha1 = self.theta
        row.alpha2 = self.phi
        row.alpha3 = self.iota
        row.alpha4 = self.psi
        row.sigmasq = self.sigmasq
        return row


class SpinTaylorT2FourierTemplate(InspiralPrecessingSpinTemplate):
    approximant = "SpinTaylorT2Fourier"

class SpinTaylorT4Template(InspiralPrecessingSpinTemplate):
    approximant = "SpinTaylorT4"

class SpinTaylorT5Template(InspiralPrecessingSpinTemplate):
    approximant = "SpinTaylorT5"

class SEOBNRv3Template(IMRPrecessingSpinTemplate):
    approximant = "SEOBNRv3"

class IMRPhenomPTemplate(IMRPrecessingSpinTemplate):
    approximant = "IMRPhenomP"

class IMRPhenomPv2Template(IMRPrecessingSpinTemplate):
    approximant = "IMRPhenomPv2"


waveforms = {
    "TaylorF2RedSpin": TaylorF2RedSpinTemplate,
    "TaylorF2" : TaylorF2Template,
    "IMRPhenomB": IMRPhenomBTemplate,
    "IMRPhenomC": IMRPhenomCTemplate,
    "IMRPhenomD": IMRPhenomDTemplate,
    "IMRPhenomP": IMRPhenomPTemplate,
    "IMRPhenomPv2": IMRPhenomPv2Template,
    "SEOBNRv2": SEOBNRv2Template,
    "SEOBNRv2_ROM_DoubleSpin": SEOBNRv2ROMDoubleSpinTemplate,
    "SEOBNRv2_ROM_DoubleSpin_HI": SEOBNRv2ROMDoubleSpinHITemplate,
    "EOBNRv2": EOBNRv2Template,
    "SpinTaylorT4": SpinTaylorT4Template,
    "SpinTaylorT5": SpinTaylorT5Template,
    "SpinTaylorF2": SpinTaylorF2Template,
    "SpinTaylorT2Fourier": SpinTaylorT2FourierTemplate,
    "SEOBNRv3":SEOBNRv3Template
}
