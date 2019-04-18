# Copyright (C) 2012 Chris Pankow, Evan Ochsner, Richard O'Shaughnessy
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
Integrate the extrinsic parameters of the prefactored likelihood function.
"""

from __future__ import print_function

# Basic stuff
import sys
import functools
from argparse import ArgumentParser

# analysis stuff
import numpy

# LAL analysis stuff
import lal
from glue.ligolw import utils, lsctables, table, ligolw
lsctables.use_in(ligolw.LIGOLWContentHandler)
from glue.ligolw.utils import process
import glue.lal

# our analysis stuff
from lalinference.rapid_pe import lalsimutils, factored_likelihood, mcsampler, xmlutils, common_cl
from lalinference.rapid_pe.common_cl import param_limits

from lalinference.io import fits as bfits

__author__ = "Evan Ochsner <evano@gravity.phys.uwm.edu>, Chris Pankow <pankow@gravity.phys.uwm.edu>, R. O'Shaughnessy"

#
# Option parsing
#
argp = ArgumentParser()
common_cl.add_datasource_params(argp)
common_cl.add_output_params(argp)

#
# Add the integration options
#
common_cl.add_integration_params(argp)

#
# Add the intrinsic parameters
#
common_cl.add_intrinsic_params(argp)

#
# Add the pinnable parameters
#
common_cl.add_pinnable_params(argp)

opts = argp.parse_args()

# Check both or neither of --data-start/end-time given
if opts.data_start_time is None and opts.data_end_time is not None:
    raise ValueError("You must provide both or neither of --data-start-time and --data-end-time.")
if opts.data_end_time is None and opts.data_start_time is not None:
    raise ValueError("You must provide both or neither of --data-start-time and --data-end-time.")

#
# Hardcoded variables
#
t_ref_wind = 50e-3 # Interpolate in a window +/- this width about event time. 
T_safety = 2. # Safety buffer (in sec) for wraparound corruption

#
# Set distance maximum
#
param_limits["distance"] = (param_limits["distance"][0], opts.distance_maximum)

#
# Inverse spectrum truncation control
#
T_spec = opts.inv_spec_trunc_time
if T_spec == 0.: # Do not do inverse spectrum truncation
    inv_spec_trunc_Q = False
    T_safety += 8. # Add a bit more safety buffer in this case
else:
    inv_spec_trunc_Q = True

#
# Initialize the RNG, if needed
#
# TODO: Do we seed a given instance of the integrator, or set it for all
# or both?
if opts.seed is not None:
    numpy.random.seed(opts.seed)

#
# Gather information about a injection put in the data
#
if opts.pin_to_sim is not None:
    xmldoc = utils.load_filename(opts.pin_to_sim, contenthandler=ligolw.LIGOLWContentHandler)
    sim_table = lsctables.SimInspiralTable.get_table(xmldoc)
    assert len(sim_table) == 1
    sim_row = sim_table[0]

#
# Gather information from the detection pipeline
#
if opts.coinc_xml is not None:
    xmldoc = utils.load_filename(opts.coinc_xml, contenthandler=ligolw.LIGOLWContentHandler)
    coinc_table = table.get_table(xmldoc, lsctables.CoincInspiralTable.tableName)
    assert len(coinc_table) == 1
    coinc_row = coinc_table[0]
    event_time = coinc_row.get_end()
    print("Coinc XML loaded, event time: %s" % str(coinc_row.get_end()))
elif opts.event_time is not None:
    event_time = glue.lal.LIGOTimeGPS(opts.event_time)
    print("Event time from command line: %s" % str(event_time))
else:
    raise ValueError("Either --coinc-xml or --event-time must be provided to parse event time.")

#
# Set masses
#
sngl_inspiral_table = None
if opts.mass1 is not None and opts.mass2 is not None:
    m1, m2 = opts.mass1, opts.mass2
elif opts.coinc_xml is not None:
    sngl_inspiral_table = table.get_table(xmldoc, lsctables.SnglInspiralTable.tableName)
    assert len(sngl_inspiral_table) == len(coinc_row.ifos.split(","))
    m1, m2 = None, None
    for sngl_row in sngl_inspiral_table:
        # NOTE: gstlal is exact match, but other pipelines may not be
        assert m1 is None or (sngl_row.mass1 == m1 and sngl_row.mass2 == m2)
        m1, m2 = sngl_row.mass1, sngl_row.mass2
else:
    raise ValueError("Need either --mass1 --mass2 or --coinc-xml to retrieve masses.")

lambda1, lambda2 = 0, 0
if opts.eff_lambda is not None:
    lambda1, lambda2 = lalsimutils.tidal_lambda_from_tilde(opts.mass1, opts.mass2, opts.eff_lambda, opts.deff_lambda or 0)

print("Performing integration for intrinsic parameters mass 1: %f, mass 2 %f, lambda1: %f, lambda2: %f, spin1z: %1.3f, spin2z: %1.3f" % (m1, m2, lambda1, lambda2, opts.spin1z or 0, opts.spin2z or 0))

#
# Template descriptors
#

fiducial_epoch = lal.LIGOTimeGPS(event_time.gpsSeconds, event_time.gpsNanoSeconds)

# Struct to hold template parameters
P = lalsimutils.ChooseWaveformParams(
	approx = lalsimutils.lalsim.GetApproximantFromString(opts.approximant),
    radec = False,
    fmin = opts.fmin_template, # minimum frequency of template
    fref = opts.reference_freq,
    ampO = opts.amp_order,
    tref = fiducial_epoch,

    dist = common_cl.distRef,
    incl = 0.0,
    phiref = 0.0,
    theta = 0.0,
    phi = 0.0,
    psi = 0.0,

    m1 = m1 * lal.MSUN_SI,
    m2 = m2 * lal.MSUN_SI,
    lambda1 = lambda1,
    lambda2 = lambda2
)

# Set spin values if appropriate
if opts.spin1z is not None or opts.spin2z is not None:
    for a in ("spin1z", "spin2z"):
        if getattr(opts, a) is not None:
            setattr(P, a, getattr(opts, a))
        else:
            setattr(P, a, 0.0)
# FIXME: We need to somehow be consistent here. Either use SpinTaylorT4, or set spins explicitly to zero.
"""
elif opts.coinc_xml is not None:
    xmldoc = utils.load_filename(opts.coinc_xml, contenthandler=ligolw.LIGOLWContentHandler)
    sngl_inspiral_table = table.get_table(xmldoc, lsctables.SnglInspiralTable.tableName)
    for a in ("spin1z", "spin2z"):
        if hasattr(sngl_inspiral_table[0], a):
            setattr(P, a, getattr(sngl_inspiral_table[0], a))
        else:
            setattr(P, a, 0.0)
"""

# User requested bounds for data segment
if opts.data_start_time is not None and opts.data_end_time is not None:
    start_time =  opts.data_start_time
    end_time =  opts.data_end_time
    print("Fetching data segment with start=", start_time)
    print("                             end=", end_time)

# Automatically choose data segment bounds so region of interest isn't corrupted
else:
    htmplt = lalsimutils.hoft(P)
    T_tmplt = - float(htmplt.epoch)
    T_seg = T_tmplt + T_spec + T_safety # Amount before and after event time
    start_time = float(event_time) - T_seg
    end_time = float(event_time) + T_seg
    print("Fetching data segment with start=", start_time)
    print("                             end=", end_time)
    print("\t\tEvent time is: ", float(event_time))
    print("\t\tT_seg is: ", T_seg)

#
# Load in data and PSDs
#
data_dict, psd_dict = {}, {}

if opts.zero_noise:
    if not opts.pin_to_sim:
        exit("--pin-to-sim must be provided with zero noise to do an injection.")
    for inst, _ in map(lambda c: c.split("="), opts.psd_file):
        print("Adding injections from sim inspiral row for instrument %s" % inst)
        P_inj = lalsimutils.ChooseWaveformParams()
        P_inj.copy_lsctables_sim_inspiral(sim_row)
        P_inj.detector = inst
        P_inj.radec = True # Get the right antenna patterns
        P_inj.deltaF, P_inj.fmax = 0.125/16, 2048
        data_dict[inst] = lalsimutils.non_herm_hoff(P_inj)
        print("Frequency binning: %f, length %d" % (data_dict[inst].deltaF,
                data_dict[inst].data.length))
else:
    for inst, chan in map(lambda c: c.split("="), opts.channel_name):
        print("Reading channel %s from cache %s" % (inst+":"+chan, opts.cache_file))
        data_dict[inst] = lalsimutils.frame_data_to_non_herm_hoff(opts.cache_file,
                inst+":"+chan, start=start_time, stop=end_time,
                window_shape=opts.window_shape)
        print("Frequency binning: %f, length %d" % (data_dict[inst].deltaF,
                data_dict[inst].data.length))

for inst, psdf in map(lambda c: c.split("="), opts.psd_file):
    print("Reading PSD for instrument %s from %s" % (inst, psdf))
    psd_dict[inst] = lalsimutils.get_psd_series_from_xmldoc(psdf, inst)

    deltaF = data_dict[inst].deltaF
    # Highest freq. at which PSD is defined
    fmax = psd_dict[inst].f0 + psd_dict[inst].deltaF * (psd_dict[inst].data.length - 1)

    # Assert upper limit of IP integral does not go past where PSD defined
    assert opts.fmax is None or opts.fmax<= fmax
    # Allow us to target a smaller upper limit than provided by the PSD. Important for numerical PSDs that turn over at high frequency
    if opts.fmax and opts.fmax < fmax:
        fmax = opts.fmax # fmax is now the upper freq. of IP integral

    psd_dict[inst] = lalsimutils.resample_psd_series(psd_dict[inst], deltaF)
    print("PSD deltaF after interpolation %f" % psd_dict[inst].deltaF)

    assert psd_dict[inst].deltaF == deltaF

# Ensure data and PSDs keyed to same detectors
if sorted(psd_dict.keys()) != sorted(data_dict.keys()):
    print("Got a different set of instruments based on data and PSDs provided.", file=sys.stderr)

# Ensure waveform has same sample rate, padded length as data
#
# N.B. This assumes all detector data has same sample rate, length
#
# data_dict holds 2-sided FrequencySeries, so their length is the same as
# that of the original TimeSeries that was FFT'd = Nsamples
# Also, deltaF = 1/T, with T = the duration (in sec) of the original TimeSeries
# Therefore 1/(data.length*deltaF) = T/Nsamples = deltaT
P.deltaT = 1./ (data_dict[data_dict.keys()[0]].data.length * deltaF)
P.deltaF = deltaF

#
# Perform the Precompute stage
#

# N.B. There is an implicit assumption all detectors use the same
# upper frequency limit for their inner product integrals
# N.B. P.fmin is being used to set the lower freq. limit of the IP integrals
# while in principal we may want to set it separately

t_window = 0.15
rholms_intp, cross_terms, rholms = factored_likelihood.precompute_likelihood_terms(fiducial_epoch, t_window, P, data_dict, psd_dict, opts.l_max, fmax, False, inv_spec_trunc_Q, T_spec)

if opts.pin_to_sim and not opts.zero_noise:
    P.copy_lsctables_sim_inspiral(sim_row)
    print("Pinned parameters from sim_inspiral")
    print("\tRA", P.phi, sim_row.longitude) 
    print("\tdec", P.theta, sim_row.latitude) 
    print("\tt ref %d.%d" % (P.tref.gpsSeconds, P.tref.gpsNanoSeconds), sim_row.get_time_geocent())
    print("\torb phase", P.phiref, sim_row.coa_phase) # ref. orbital phase
    print("\tinclination", P.incl, sim_row.inclination) # inclination
    print("\tpsi", P.psi, sim_row.polarization) # polarization angle
    print("\tdistance", P.dist/(1e6 * lal.PC_SI), sim_row.distance)  # luminosity distance

    logL = factored_likelihood.factored_log_likelihood(P, rholms_intp, cross_terms, opts.l_max)
    print("Pinned log likelihood: %g, (%g in \"SNR\")" % (logL, numpy.sqrt(2*logL)))
    tref = float(P.tref)
    tvals = numpy.arange(tref-0.01, tref+0.01, 0.00001)
    logLs = []
    for t in tvals:
        P.tref = lal.LIGOTimeGPS(t)
        logLs.append(factored_likelihood.factored_log_likelihood(P, rholms_intp, cross_terms, opts.l_max))
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot
    print("Maximum logL is %g, (%g in \"SNR\")" % (max(logLs), numpy.sqrt(2*max(logLs))))
    print("Which occurs at sample", numpy.argmax(logLs))
    print("This corresponds to time %.20g" % tvals[numpy.argmax(logLs)])
    print("The data event time is:  %.20g" % sim_row.get_time_geocent())
    print("Difference from geocenter t_ref is %.20g" %\
            (tvals[numpy.argmax(logLs)] - sim_row.get_time_geocent()))
    print("This difference in discrete samples: %.20g" %\
            ((tvals[numpy.argmax(logLs)]-sim_row.get_time_geocent())/P.deltaT))
    pyplot.plot(tvals-tref, logLs)
    pyplot.ylabel("log Likelihood")
    pyplot.xlabel("time (relative to %10.5f)" % tref)
    pyplot.axvline(0, color="k")
    pyplot.title("lnL(t),\n value at event time: %f" % logL)
    pyplot.grid()
    pyplot.savefig("logL.png")
    integral = numpy.sum( numpy.exp(logLs) * P.deltaT )
    print("Integral over t of likelihood is:", integral)
    print("The log of the integral is:", numpy.log(integral))
    exit()

#
# Parameter integral sampling strategy
#
params = {}
sampler = mcsampler.MCSampler()

#
# Psi -- polarization angle
# sampler: uniform in [0, pi)
#
psi_sampler = functools.partial(mcsampler.uniform_samp_vector, param_limits["psi"][0], param_limits["psi"][1])
sampler.add_parameter("psi", pdf = psi_sampler, cdf_inv = None, left_limit = param_limits["psi"][0], right_limit = param_limits["psi"][1], prior_pdf = mcsampler.uniform_samp_psi)

#
# Phi - orbital phase
# sampler: uniform in [0, 2*pi)
#
phi_sampler = functools.partial(mcsampler.uniform_samp_vector, param_limits["phi_orb"][0], param_limits["phi_orb"][1])
sampler.add_parameter("phi_orb", pdf = phi_sampler, cdf_inv = None, left_limit = param_limits["phi_orb"][0], right_limit = param_limits["phi_orb"][1], prior_pdf = mcsampler.uniform_samp_phase)

#
# inclination - angle of system angular momentum with line of sight
# sampler: cos(incl) uniform in [-1, 1)
#
incl_sampler = mcsampler.cos_samp_vector # this is NOT dec_samp_vector, because the angular zero point is different!
sampler.add_parameter("inclination", pdf = incl_sampler, cdf_inv = None, left_limit = param_limits["inclination"][0], right_limit = param_limits["inclination"][1], prior_pdf = mcsampler.uniform_samp_theta, adaptive_sampling = not opts.no_adapt)

#
# Distance - luminosity distance to source in parsecs
# sampler: uniform distance over [dmin, dmax), adaptive sampling
#
dist_sampler = functools.partial(mcsampler.uniform_samp_vector, param_limits["distance"][0], param_limits["distance"][1])
sampler.add_parameter("distance", pdf = dist_sampler, cdf_inv = None, left_limit = param_limits["distance"][0], right_limit = param_limits["distance"][1], prior_pdf = numpy.vectorize(lambda x: x**2/(param_limits["distance"][1]**3/3. - param_limits["distance"][0]**3/3.)), adaptive_sampling = not opts.no_adapt)

if opts.skymap_file is not None:
    #
    # Right ascension and declination -- use a provided skymap
    #
    smap, _ = bfits.read_sky_map(opts.skymap_file)
    ss_sampler = mcsampler.HealPixSampler(smap)
    isotropic_bstar_sampler = numpy.vectorize(lambda dec, ra: 1.0/len(smap))

    # FIXME: Should the left and right limits be modified?
    sampler.add_parameter(("declination", "right_ascension"), \
        pdf = ss_sampler.pseudo_pdf, \
        cdf_inv = ss_sampler.pseudo_cdf_inverse, \
        left_limit = (param_limits["declination"][0], param_limits["right_ascension"][0]), \
        right_limit = (param_limits["declination"][1], param_limits["right_ascension"][1]), \
        prior_pdf = isotropic_bstar_sampler)

else:
    #
    # Right ascension - angle in radians from prime meridian plus hour angle
    # sampler: uniform in [0, 2pi), adaptive sampling
    #
    ra_sampler = functools.partial(mcsampler.uniform_samp_vector, param_limits["right_ascension"][0], param_limits["right_ascension"][1])
    sampler.add_parameter("right_ascension", \
        pdf = ra_sampler, \
        cdf_inv = None, \
        left_limit = param_limits["right_ascension"][0], \
        right_limit =  param_limits["right_ascension"][1], \
        prior_pdf = mcsampler.uniform_samp_phase, \
        adaptive_sampling = not opts.no_adapt)

    #
    # declination - angle in radians from the north pole piercing the celestial
    # sky sampler: cos(dec) uniform in [-1, 1), adaptive sampling
    #
    dec_sampler = mcsampler.dec_samp_vector
    sampler.add_parameter("declination", \
        pdf = dec_sampler, \
        cdf_inv = None, \
        left_limit = param_limits["declination"][0], \
        right_limit = param_limits["declination"][1], \
        prior_pdf = mcsampler.uniform_samp_dec, \
        adaptive_sampling = not opts.no_adapt)

#
# Determine pinned and non-pinned parameters
#

pinned_params = common_cl.get_pinned_params(opts)
unpinned_params = common_cl.get_unpinned_params(opts, sampler.params)
print("{0:<25s} {1:>5s} {2:>5s} {3:>20s} {4:<10s}".format("parameter", "lower limit", "upper limit", "pinned?", "pin value"))
plen = len(sorted(sampler.params, key=lambda p: len(p))[-1])
for p in sampler.params:
    if pinned_params.has_key(p):
        pinned, value = True, "%1.3g" % pinned_params[p]
    else:
        pinned, value = False, ""

    if isinstance(p, tuple):
        for subp, subl, subr in zip(p, sampler.llim[p], sampler.rlim[p]):
            subp = subp + " "*min(0, plen-len(subp))
            print("|{0:<25s} {1:>1.3g}   {2:>1.3g} {3:>20s} {4:<10s}".format(subp, subl, subr, str(False), ""))
    else:
        p = p + " "*min(0, plen-len(p))
        print("{0:<25s} {1:>1.3g}   {2:>1.3g} {3:>20s} {4:<10s}".format(p, sampler.llim[p], sampler.rlim[p], str(pinned), value))

# Special case: t_ref is assumed to be relative to the epoch
if pinned_params.has_key("t_ref"):
    pinned_params["t_ref"] -= float(fiducial_epoch)

#
# Merge options into one big ol' kwargs dict
#

pinned_params.update({ 
    # Iteration settings and termination conditions
    "n": min(opts.n_chunk, opts.n_max), # Number of samples in a chunk
    "nmax": opts.n_max, # Total number of samples to draw before termination
    "neff": opts.n_eff, # Total number of effective samples to collect before termination

    # Adaptive sampling settings
    "tempering_exp": opts.adapt_weight_exponent if not opts.no_adapt else 1.0, # Weights will be raised to this power to prevent overconvergence
    "floor_level": opts.adapt_floor_level if not opts.no_adapt else 0.0, # The new sampling distribution at the end of each chunk will be floor_level-weighted average of a uniform distribution and the (L^tempering_exp p/p_s)-weighted histogram of sampled points.
    "history_mult": 10, # Multiplier on 'n' - number of samples to estimate marginalized 1-D histograms
    "n_adapt": 100 if not opts.no_adapt else 0, # Number of chunks to allow adaption over

    # Verbosity settings
    "verbose": True, 
    "extremely_verbose": False, 

    # Sample caching
    "save_intg": opts.save_samples, # Cache the samples (and integrand values)?
    "igrand_threshold_deltalnL": opts.save_deltalnL, # Threshold on distance from max L to save sample
    "igrand_threshold_p": opts.save_P # Threshold on cumulative probability contribution to cache sample
})

#
# Call the likelihood function for various extrinsic parameter values
#
if not opts.time_marginalization:

    #
    # tref - GPS time of geocentric end time
    # sampler: uniform in +/-2 ms window around estimated end time 
    #
    tref_sampler = functools.partial(mcsampler.uniform_samp_vector, param_limits["t_ref"][0], param_limits["t_ref"][1])

    sampler.add_parameter("t_ref", \
                          pdf = tref_sampler, \
                          cdf_inv = None, \
                          left_limit = param_limits["t_ref"][0], \
                          right_limit = param_limits["t_ref"][1], \
                          prior_pdf = functools.partial(mcsampler.uniform_samp_vector, param_limits["t_ref"][0], param_limits["t_ref"][1]))

    #
    # A note of caution:
    # In order to make the pinning interface work consistently, the names of 
    # parameters given to the sampler must match the argument names in the
    # called function. This is because the sampler has to reconstruct the
    # argument order to pass the right values, and it can only do that by
    # comparing the parameter names it knows to the arguments that are passed
    # to it.
    #
    def likelihood_function(right_ascension, declination, t_ref, phi_orb, inclination, psi, distance):

        P.phi = right_ascension # right ascension
        P.theta = declination # declination
        P.tref = t_ref + fiducial_epoch # ref. time (rel to epoch for data taking)
        P.phiref = phi_orb # ref. orbital phase
        P.incl = inclination # inclination
        P.psi = psi # polarization angle
        P.dist = distance * 1.e6 * lal.PC_SI # luminosity distance
    
        lnL = factored_likelihood.factored_log_likelihood(P, rholms_intp, cross_terms, opts.l_max)
    
        return numpy.exp(lnL)

    res, var, neff, dict_return = sampler.integrate(likelihood_function, *unpinned_params, **pinned_params)

else: # Sum over time for every point in other extrinsic params

    # Un-LALify the rholms, we don't need the extra information and we end up
    # invoking lots of data copies in slicing later on.
    det_epochs = {}
    for det, rho_ts in rholms.iteritems():
        ts_epoch = None
        for mode, mode_ts in rho_ts.iteritems():
            if ts_epoch:
                # Be super careful that we're looking at the *same* start times
                assert ts_epoch == mode_ts.epoch
            else:
                ts_epoch = mode_ts.epoch
            rholms[det][mode] = mode_ts.data.data

        det_epochs[det] = ts_epoch

    # Construct the time grid and advance it to the epoch of the mode SNR
    # time series
    tvals = numpy.linspace(-t_ref_wind, t_ref_wind, int(t_ref_wind * 2 / P.deltaT))

    def likelihood_function(right_ascension, declination, phi_orb, inclination,
            psi, distance):

        P.phi = right_ascension # right ascension
        P.theta = declination # declination
        P.tref = fiducial_epoch  # see 'tvals', above
        P.phiref = phi_orb # ref. orbital phase
        P.incl = inclination # inclination
        P.psi = psi # polarization angle
        P.dist = distance * 1.e6 * lal.PC_SI # luminosity distance

        lnL = factored_likelihood.factored_log_likelihood_time_marginalized(tvals, P, rholms_intp, rholms, cross_terms, det_epochs, opts.l_max,interpolate=opts.interpolate_time)
    
        return numpy.exp(lnL)

    res, var, neff, dict_return = sampler.integrate(likelihood_function, *unpinned_params, **pinned_params)

print(" lnLmarg is ", numpy.log(res), " with expected relative error ", numpy.sqrt(var)/res)
print(" note neff is ", neff)

# FIXME: This is turning into a mess
if opts.output_file:
    xmldoc = ligolw.Document()
    xmldoc.appendChild(ligolw.LIGO_LW())
    process.register_to_xmldoc(xmldoc, sys.argv[0], opts.__dict__)
    result_dict = {"mass1": opts.mass1, "mass2": opts.mass2, "event_duration": numpy.sqrt(var)/res, "ttotal": sampler.ntotal}
    if opts.spin1z is not None or sngl_inspiral_table:
        result_dict["spin1z"] = opts.spin1z or 0.0
    if opts.spin2z is not None or sngl_inspiral_table:
        result_dict["spin2z"] = opts.spin2z or 0.0
    if opts.eff_lambda is not None:
        result_dict["psi0"] = opts.eff_lambda
    if opts.deff_lambda is not None:
        result_dict["psi3"] = opts.deff_lambda

    xmlutils.append_likelihood_result_to_xmldoc(xmldoc, numpy.log(res), neff=neff, **result_dict)
    utils.write_filename(xmldoc, opts.output_file, gz=opts.output_file.endswith(".gz"))
    if opts.save_samples:
        samples = sampler._rvs
        samples["distance"] = samples["distance"]
        if not opts.time_marginalization:
            samples["t_ref"] += float(fiducial_epoch)
        else:
            samples["t_ref"] = float(fiducial_epoch)*numpy.ones(len(sampler._rvs["psi"]))
        samples["polarization"] = samples["psi"]
        samples["coa_phase"] = samples["phi_orb"]
        if ("declination", "right_ascension") in sampler.params:
            samples["latitude"], samples["longitude"] = samples[("declination", "right_ascension")]
        else:
            samples["latitude"] = samples["declination"]
            samples["longitude"] = samples["right_ascension"]
        samples["loglikelihood"] = numpy.log(samples["integrand"])
        samples["mass1"] = numpy.empty(samples["psi"].shape)
        samples["mass2"] = numpy.empty(samples["psi"].shape)
        samples["mass1"].fill(opts.mass1)
        samples["mass2"].fill(opts.mass2)
        samples["spin1z"] = numpy.empty(samples["psi"].shape)
        samples["spin2z"] = numpy.empty(samples["psi"].shape)
        samples["spin1z"].fill(opts.spin1z or 0.0)
        samples["spin2z"].fill(opts.spin2z or 0.0)
        if opts.eff_lambda is not None or opts.deff_lambda is not None:
            # FIXME: the column mapping isn't working right, we need to fix that
            # rather than give these weird names
            samples["psi0"] = numpy.ones(samples["psi"].shape)
            samples["psi3"] = numpy.ones(samples["psi"].shape)
            samples["psi0"].fill(opts.eff_lambda)
            samples["psi3"].fill(opts.deff_lambda)
        xmlutils.append_samples_to_xmldoc(xmldoc, samples)
        sample_output = opts.output_file.split(".")
        sample_output = sample_output[0] + ".samples." + ".".join(sample_output[1:])
        utils.write_filename(xmldoc, sample_output, gz=opts.output_file.endswith(".gz"))
