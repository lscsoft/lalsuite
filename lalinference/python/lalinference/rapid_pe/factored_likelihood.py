# Copyright (C) 2013  Evan Ochsner, R. O'Shaughnessy
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

"""
Code to compute the log likelihood of parameters of a gravitational waveform. Precomputes terms that depend only on intrinsic parameters and computes the log likelihood for given values of extrinsic parameters

Requires python SWIG bindings of the LIGO Algorithms Library (LAL)
"""

import lal
import lalsimulation as lalsim
from lalinference.rapid_pe import lalsimutils as lsu
import numpy as np
from scipy import interpolate, integrate
from scipy import special
from itertools import product
from common_cl import distRef

__author__ = "Evan Ochsner <evano@gravity.phys.uwm.edu>, R. O'Shaughnessy"

#
# Main driver functions
#
def precompute_likelihood_terms(event_time_geo, t_window, P, data_dict,
        psd_dict, Lmax, fMax, analyticPSD_Q=False,
        inv_spec_trunc_Q=False, T_spec=0., remove_zero=True, verbose=True):
    """
    Compute < h_lm(t) | d > and < h_lm | h_l'm' >

    Returns:
        - Dictionary of interpolating functions, keyed on detector, then (l,m)
          e.g. rholms_intp['H1'][(2,2)]
        - Dictionary of "cross terms" <h_lm | h_l'm' > keyed on (l,m),(l',m')
          e.g. crossTerms[((2,2),(2,1))]
        - Dictionary of discrete time series of < h_lm(t) | d >, keyed the same
          as the interpolating functions.
          Their main use is to validate the interpolating functions
    """
    assert data_dict.keys() == psd_dict.keys()
    detectors = data_dict.keys()
    rholms = {}
    rholms_intp = {}
    crossTerms = {}

    # Compute hlms at a reference distance, distance scaling is applied later
    P.dist = distRef

    P.print_params()
    # Compute all hlm modes with l <= Lmax
    detectors = data_dict.keys()
    # Zero-pad to same length as data - NB: Assuming all FD data same resolution
    P.deltaF = data_dict[detectors[0]].deltaF
    hlms_list = lsu.hlmoff(P, Lmax) # a linked list of hlms
    hlms = lsu.SphHarmFrequencySeries_to_dict(hlms_list, Lmax) # a dictionary

    # If the hlm time series is identically zero, remove it
    zero_modes = []
    for mode in hlms.iterkeys():
        if remove_zero and np.sum(np.abs(hlms[mode].data.data)) == 0:
            zero_modes.append(mode)

    for mode in zero_modes:
        del hlms[mode]

    for det in detectors:
        # This is the event time at the detector
        t_det = compute_arrival_time_at_detector(det, P.phi, P.theta,event_time_geo)
        # The is the difference between the time of the leading edge of the
        # time window we wish to compute the likelihood in, and
        # the time corresponding to the first sample in the rholms
        rho_epoch = data_dict[det].epoch - hlms[hlms.keys()[0]].epoch
        t_shift =  float(t_det - t_window - rho_epoch)
        assert t_shift > 0
        # tThe leading edge of our time window of interest occurs
        # this many samples into the rholms
        N_shift = int( t_shift / P.deltaT )
        # Number of samples in the window [t_ref - t_window, t_ref + t_window]
        N_window = int( 2 * t_window / P.deltaT )
        # Compute cross terms < h_lm | h_l'm' >
        crossTerms[det] = compute_mode_cross_term_ip(hlms, psd_dict[det], P.fmin,
                fMax, 1./2./P.deltaT, P.deltaF, analyticPSD_Q,
                inv_spec_trunc_Q, T_spec)
        # Compute rholm(t) = < h_lm(t) | d >
        rholms[det] = compute_mode_ip_time_series(hlms, data_dict[det],
                psd_dict[det], P.fmin, fMax, 1./2./P.deltaT, N_shift, N_window,
                analyticPSD_Q, inv_spec_trunc_Q, T_spec)
        rhoXX = rholms[det][rholms[det].keys()[0]]
        # The vector of time steps within our window of interest
        # for which we have discrete values of the rholms
        # N.B. I don't do simply rho_epoch + t_shift, b/c t_shift is the
        # precise desired time, while we round and shift an integer number of
        # steps of size deltaT
        t = np.arange(N_window) * P.deltaT\
                + float(rho_epoch + N_shift * P.deltaT )
        if verbose:
            print("For detector", det, "...")
            print("\tData starts at %.20g" % float(data_dict[det].epoch))
            print("\trholm starts at %.20g" % float(rho_epoch))
            print("\tEvent time at detector is: %.18g" % float(t_det))
            print("\tInterpolation window has half width %g" % t_window)
            print("\tComputed t_shift = %.20g" % t_shift)
            print("\t(t_shift should be t_det - t_window - t_rholm = %.20g)" %\
                    (t_det - t_window - float(rho_epoch)))
            print("\tInterpolation starts at time %.20g" % t[0])
            print("\t(Should start at t_event - t_window = %.20g)" %\
                    (float(rho_epoch + N_shift * P.deltaT)))
        # The minus N_shift indicates we need to roll left
        # to bring the desired samples to the front of the array
        rholms_intp[det] =  interpolate_rho_lms(rholms[det], t)

    return rholms_intp, crossTerms, rholms

def factored_log_likelihood(extr_params, rholms_intp, crossTerms, Lmax):
    """
    Compute the log-likelihood = -1/2 < d - h | d - h > from:
        - extr_params is an object containing values of all extrinsic parameters
        - rholms_intp is a dictionary of interpolating functions < h_lm(t) | d >
        - crossTerms is a dictionary of < h_lm | h_l'm' >
        - Lmax is the largest l-index of any h_lm mode considered

    N.B. rholms_intp and crossTerms are the first two outputs of the function
    'precompute_likelihood_terms'
    """
    # Sanity checks
    assert rholms_intp.keys() == crossTerms.keys()
    detectors = rholms_intp.keys()

    RA = extr_params.phi
    DEC =  extr_params.theta
    tref = extr_params.tref # geocenter time
    phiref = extr_params.phiref
    incl = extr_params.incl
    psi = extr_params.psi
    dist = extr_params.dist

    # use EXTREMELY many bits
    i = 0
    _lnL = np.zeros(np.shape(extr_params.phi), dtype=np.float128)

    for RA, DEC, phiref, incl, psi, tref, dist in zip(extr_params.phi, \
                                        extr_params.theta, extr_params.phiref, \
                                        extr_params.incl, extr_params.psi, \
                                        extr_params.tref, extr_params.dist):
        # N.B.: The Ylms are a function of - phiref b/c we are passively rotating
        # the source frame, rather than actively rotating the binary.
        # Said another way, the m^th harmonic of the waveform should transform as
        # e^{- i m phiref}, but the Ylms go as e^{+ i m phiref}, so we must give
        # - phiref as an argument so Y_lm h_lm has the proper phiref dependence
        # FIXME: Strictly speaking, this should be inside the detector loop because
        # there *could* be different l,m pairs for different detectors. This never
        # happens in practice, so it's pulled out here, and we use the first
        # detector as a reference.
        Ylms = compute_spherical_harmonics(Lmax, incl, -phiref, rholms_intp[rholms_intp.keys()[0]])

        lnL = 0.
        for det in detectors:
            CT = crossTerms[det]

            # This is the GPS time at the detector
            t_det = compute_arrival_time_at_detector(det, RA, DEC, tref)
            F = complex_antenna_factor(det, RA, DEC, psi, t_det)

            det_rholms = {}  # rholms evaluated at time at detector
            for key in rholms_intp[det]:
                func = rholms_intp[det][key]
                det_rholms[key] = func(float(t_det))

            lnL += single_detector_log_likelihood(det_rholms, CT, Ylms, F, dist)

        _lnL[i] = lnL
        i += 1

    return _lnL

def factored_log_likelihood_time_marginalized(tvals, extr_params, rholms_intp, rholms, crossTerms, det_epochs, Lmax, interpolate=False):
    """
    Compute the log-likelihood = -1/2 < d - h | d - h > from:
        - extr_params is an object containing values of all extrinsic parameters
        - rholms_intp is a dictionary of interpolating functions < h_lm(t) | d >
        - crossTerms is a dictionary of < h_lm | h_l'm' >
        - Lmax is the largest l-index of any h_lm mode considered

    tvals is an array of timeshifts relative to the detector,
    used to compute the marginalized integral.
    It provides both the time prior and the sample points used for the integral.

    N.B. rholms_intp and crossTerms are the first two outputs of the function
    'precompute_likelihood_terms'
    """
    # Sanity checks
    assert rholms_intp.keys() == crossTerms.keys()
    detectors = rholms_intp.keys()

    RA = extr_params.phi
    DEC =  extr_params.theta
    tref = extr_params.tref # geocenter time
    phiref = extr_params.phiref
    incl = extr_params.incl
    psi = extr_params.psi
    dist = extr_params.dist

    # use EXTREMELY many bits
    i = 0
    _lnL = np.zeros(np.shape(extr_params.phi), dtype=np.float128)

    for RA, DEC, phiref, incl, psi, dist in zip(extr_params.phi, \
                                        extr_params.theta, extr_params.phiref, \
                                        extr_params.incl, extr_params.psi, \
                                        extr_params.dist):

        # N.B.: The Ylms are a function of - phiref b/c we are passively rotating
        # the source frame, rather than actively rotating the binary.
        # Said another way, the m^th harmonic of the waveform should transform as
        # e^{- i m phiref}, but the Ylms go as e^{+ i m phiref}, so we must give
        # - phiref as an argument so Y_lm h_lm has the proper phiref dependence
        # FIXME: Strictly speaking, this should be inside the detector loop because
        # there *could* be different l,m pairs for different detectors. This never
        # happens in practice, so it's pulled out here, and we use the first
        # detector as a reference.
        Ylms = compute_spherical_harmonics(Lmax, incl, -phiref, rholms[rholms.keys()[0]])

        lnL = 0.
        delta_t = tvals[1] - tvals[0]
        for det in detectors:
            CT = crossTerms[det]
            F = complex_antenna_factor(det, RA, DEC, psi, tref)
            rho_epoch = float(det_epochs[det])

            # This is the GPS time at the detector
            t_det = compute_arrival_time_at_detector(det, RA, DEC, tref)
            det_rholms = {}  # rholms evaluated at time at detector
            if ( interpolate ):
                # use the interpolating functions. 
                for key, func in rholms_intp[det].iteritems():
                    det_rholms[key] = func(float(t_det)+tvals)
            else:
                # do not interpolate, just use nearest neighbors.
                for key, rhoTS in rholms[det].iteritems():
                # PRB: these can be moved outside this loop to after t_det
                    tfirst = float(t_det)+tvals[0]
                    ifirst = int((tfirst - rho_epoch) / delta_t + 0.5)
                    ilast = ifirst + len(tvals)
                    det_rholms[key] = rhoTS[ifirst:ilast]

            lnL += single_detector_log_likelihood(det_rholms, CT, Ylms, F, dist)

        maxlnL = np.max(lnL)
        _lnL[i] = maxlnL + np.log(np.sum(np.exp(lnL - maxlnL)) * (tvals[1]-tvals[0]))
        i += 1

    return _lnL

def single_detector_log_likelihood(rholm_vals, crossTerms, Ylms, F, dist):
    """
    Compute the value of the log-likelihood at a single detector from
    several intermediate pieces of data.

    Inputs:
      - rholm_vals: A dictionary of values of inner product between data
            and h_lm modes, < h_lm(t*) | d >, at a single time of interest t*
      - crossTerms: A dictionary of inner products between h_lm modes:
            < h_lm | h_l'm' >
      - Ylms: Dictionary of values of -2-spin-weighted spherical harmonic modes
            for a certain inclination and ref. phase, Y_lm(incl, - phiref)
      - F: Complex-valued antenna pattern depending on sky location and
            polarization angle, F = F_+ + i F_x
      - dist: The distance from the source to detector in meters

    Outputs: The value of ln L for a single detector given the inputs.
    """

    invDistMpc = distRef/dist
    Fstar = np.conj(F)

    term1, term20, term21 = 0., 0., 0.
    # PRB: I think this loop can be vectorized with some work
    for pair1, Ylm1 in Ylms.iteritems():
        l1, m1 = pair1
        n_one_l1 = (-1)**l1
        Ylm1_conj = np.conj(Ylm1)
        term1 += Ylm1_conj * rholm_vals[pair1]
        tmp_term20, tmp_term21 = 0., 0.
	    # PRB: should also re-pack the crossterms into arrays
        for pair2, Ylm2 in Ylms.iteritems():
            tmp_term20 += crossTerms[(pair1, pair2)] * Ylm2
            tmp_term21 += Ylm2 * crossTerms[((l1, -m1), pair2)]
        term20 += tmp_term20 * Ylm1_conj
        term21 += tmp_term21 * n_one_l1 * Ylm1
    term1 = np.real( Fstar * term1 ) * invDistMpc 
    term1 += -0.25 * np.real( F * ( Fstar * term20 + F * term21 ) ) * invDistMpc * invDistMpc 

    return term1

def compute_mode_ip_time_series(hlms, data, psd, fmin, fMax, fNyq,
        N_shift, N_window, analyticPSD_Q=False,
        inv_spec_trunc_Q=False, T_spec=0.):
    """
    Compute the complex-valued overlap between
    each member of a SphHarmFrequencySeries 'hlms'
    and the interferometer data COMPLEX16FrequencySeries 'data',
    weighted the power spectral density REAL8FrequencySeries 'psd'.

    The integrand is non-zero in the range: [-fNyq, -fmin] union [fmin, fNyq].
    This integrand is then inverse-FFT'd to get the inner product
    at a discrete series of time shifts.

    Returns a SphHarmTimeSeries object containing the complex inner product
    for discrete values of the reference time tref.  The epoch of the
    SphHarmTimeSeries object is set to account for the transformation
    """
    rholms = {}
    assert data.deltaF == hlms[hlms.keys()[0]].deltaF
    assert data.data.length == hlms[hlms.keys()[0]].data.length
    deltaT = data.data.length/(2*fNyq)

    # Create an instance of class to compute inner product time series
    IP = lsu.ComplexOverlap(fmin, fMax, fNyq, data.deltaF, psd,
            analyticPSD_Q, inv_spec_trunc_Q, T_spec, full_output=True)

    # Loop over modes and compute the overlap time series
    for pair in hlms.keys():
        rho, rhoTS, rhoIdx, rhoPhase = IP.ip(hlms[pair], data)
        rhoTS.epoch = data.epoch - hlms[pair].epoch
        rholms[pair] = lal.CutCOMPLEX16TimeSeries(rhoTS, N_shift, N_window)

    return rholms

def interpolate_rho_lm(rholm, t):
    h_re = np.real(rholm.data.data)
    h_im = np.imag(rholm.data.data)
    # spline interpolate the real and imaginary parts of the time series
    h_real = interpolate.InterpolatedUnivariateSpline(t, h_re, k=3)
    h_imag = interpolate.InterpolatedUnivariateSpline(t, h_im, k=3)
    return lambda ti: h_real(ti) + 1j*h_imag(ti)

    # Little faster
    #def anon_intp(ti):
        #idx = np.searchsorted(t, ti)
        #return rholm.data.data[idx]
    #return anon_intp

    #from pygsl import spline
    #spl_re = spline.cspline(len(t))
    #spl_im = spline.cspline(len(t))
    #spl_re.init(t, np.real(rholm.data.data))
    #spl_im.init(t, np.imag(rholm.data.data))
    #@profile
    #def anon_intp(ti):
        #re = spl_re.eval_e_vector(ti)
        #return re + 1j*im
    #return anon_intp

    # Doesn't work, hits recursion depth
    #from scipy.signal import cspline1d, cspline1d_eval
    #re_coef = cspline1d(np.real(rholm.data.data))
    #im_coef = cspline1d(np.imag(rholm.data.data))
    #dx, x0 = rholm.deltaT, float(rholm.epoch)
    #return lambda ti: cspline1d_eval(re_coef, ti) + 1j*cspline1d_eval(im_coef, ti)


def interpolate_rho_lms(rholms, t):
    """
    Return a dictionary keyed on mode index tuples, (l,m)
    where each value is an interpolating function of the overlap against data
    as a function of time shift:
    rholm_intp(t) = < h_lm(t) | d >

    'rholms' is a dictionary keyed on (l,m) containing discrete time series of
    < h_lm(t_i) | d >
    't' is an array of the discrete times:
    [t_0, t_1, ..., t_N]
    """
    rholm_intp = {}
    for mode in rholms.keys():
        rholm = rholms[mode]
        # The mode is identically zero, don't bother with it
        if sum(abs(rholm.data.data)) == 0.0:
            continue
        rholm_intp[ mode ] = interpolate_rho_lm(rholm, t)

    return rholm_intp

def compute_mode_cross_term_ip(hlms, psd, fmin, fMax, fNyq, deltaF,
        analyticPSD_Q=False, inv_spec_trunc_Q=False, T_spec=0., verbose=True):
    """
    Compute the 'cross terms' between waveform modes, i.e.
    < h_lm | h_l'm' >.
    The inner product is weighted by power spectral density 'psd' and
    integrated over the interval [-fNyq, -fmin] union [fmin, fNyq]

    Returns a dictionary of inner product values keyed by tuples of mode indices
    i.e. ((l,m),(l',m'))
    """
    # Create an instance of class to compute inner product
    IP = lsu.ComplexIP(fmin, fMax, fNyq, deltaF, psd, analyticPSD_Q,
            inv_spec_trunc_Q, T_spec)

    crossTerms = {}

    for mode1 in hlms.keys():
        for mode2 in hlms.keys():
            crossTerms[ (mode1,mode2) ] = IP.ip(hlms[mode1], hlms[mode2])
            if verbose:
                print("       : U populated ", (mode1, mode2), "  = ",\
                        crossTerms[(mode1,mode2) ])

    return crossTerms


def complex_antenna_factor(det, RA, DEC, psi, tref):
    """
    Function to compute the complex-valued antenna pattern function:
    F+ + i Fx

    'det' is a detector prefix string (e.g. 'H1')
    'RA' and 'DEC' are right ascension and declination (in radians)
    'psi' is the polarization angle
    'tref' is the reference GPS time
    """
    detector = lalsim.DetectorPrefixToLALDetector(det)
    Fp, Fc = lal.ComputeDetAMResponse(detector.response, RA, DEC, psi, lal.GreenwichMeanSiderealTime(tref))

    return Fp + 1j * Fc

def compute_spherical_harmonics(Lmax, theta, phi, selected_modes=None):
    """
    Return a dictionary keyed by tuples
    (l,m)
    that contains the values of all
    -2Y_lm(theta,phi)
    with
    l <= Lmax
    -l <= m <= l
    """
    # PRB: would this be faster if we made it a 2d numpy array?  
    Ylms = {}
    for l in range(2,Lmax+1):
        for m in range(-l,l+1):
            if selected_modes is not None and (l,m) not in selected_modes:
                continue
            Ylms[ (l,m) ] = lal.SpinWeightedSphericalHarmonic(theta, phi,-2, l, m)

    return Ylms

def compute_arrival_time_at_detector(det, RA, DEC, tref):
    """
    Function to compute the time of arrival at a detector
    from the time of arrival at the geocenter.

    'det' is a detector prefix string (e.g. 'H1')
    'RA' and 'DEC' are right ascension and declination (in radians)
    'tref' is the reference time at the geocenter.  It can be either a float (in which case the return is a float) or a GPSTime object (in which case it returns a GPSTime)
    """
    detector = lalsim.DetectorPrefixToLALDetector(det)
    # if tref is a float or a GPSTime object,
    # it shoud be automagically converted in the appropriate way
    return tref + lal.TimeDelayFromEarthCenter(detector.location, RA, DEC, tref)

def compute_mode_iterator(Lmax):  # returns a list of (l,m) pairs covering all modes, as a list.  Useful for building iterators without nested lists
    mylist = []
    for L in np.arange(2, Lmax+1):
        for m in np.arange(-L, L+1):
            mylist.append((L,m))
    return mylist
