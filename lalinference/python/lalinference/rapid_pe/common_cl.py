# Copyright (C) 2013, 2015 Chris Pankow
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
Modularized command line options for rapidpe based programs
"""

from argparse import ArgumentParser

import numpy

import lal

#
# Command line parsing utilities
#
def parse_cl_key_value(params):
    """
    Convenience in parsing out parameter arrays in the form of something=something_else:
        --channel-name "H1=FAKE-STRAIN"
    """
    return dict([val.split("=") for val in params])

#
# Pinnable parameters -- for command line processing
#
LIKELIHOOD_PINNABLE_PARAMS = ["right_ascension", "declination", "psi", "distance", "phi_orb", "t_ref", "inclination"]

def get_pinned_params(opts):
    """
    Retrieve a dictionary of user pinned parameters and their pin values.
    """
    return dict([(p,v) for p, v in opts.__dict__.iteritems() if p in LIKELIHOOD_PINNABLE_PARAMS and v is not None]) 

def get_unpinned_params(opts, params):
    """
    Retrieve a set of unpinned parameters.
    """
    return params - set([p for p, v in opts.__dict__.iteritems() if p in LIKELIHOOD_PINNABLE_PARAMS and v is not None])

#
# Add the pinnable parameters
#
def add_pinnable_params(argp, include=None, exclude=None):
    """
    Add a set of command line options with the ability to pin a parameter to a value. Pinnable parameters are held in the LIKELIHOOD_PINNABLE_PARAMS list in the common_cl module.
    Specifying include will include only those parameters both in 'include' and LIKELIHOOD_PINNABLE_PARAMS. Specifying exclude will exclude parameters present in both 'exclude' in LIKELIHOOD_PINNABLE_PARAMS.
    Exclude and include cannot both be non-None. In theory, this isn't a problem, just a question of the order in which the operations are applied.
    """
    # FIXME: In theory, this isn't a problem, just a question of the order in
    # which the operations are applied
    if exclude is not None and include is not None:
        raise ValueError("Can't specify both include and exclude")

    pinned_params = LIKELIHOOD_PINNABLE_PARAMS
    if include is not None:
        pinned_params = list(set(include) & set(LIKELIHOOD_PINNABLE_PARAMS))
    if exclude is not None:
        pinned_params = list(set(pinnable_params)- set(exclude))

    pinnable = argp.add_argument_group("Pinnable Parameters", "Specifying these command line options will pin the value of that parameter to the specified value with a probability of unity.")
    for pin_param in pinned_params:
        option = "--" + pin_param.replace("_", "-")
        pinnable.add_argument(option, type=float, help="Pin the value of %s." % pin_param)
    return argp

#
# Data source options for command line processing
#
def add_datasource_params(argp):
    argp.add_argument("-c", "--cache-file", default=None, help="LIGO cache file containing all data needed.")
    argp.add_argument("-C", "--channel-name", action="append", help="instrument=channel-name, e.g. H1=FAKE-STRAIN. Can be given multiple times for different instruments.")
    argp.add_argument("-p", "--psd-file", action="append", help="instrument=psd-file, e.g. H1=H1_PSD.xml.gz. Can be given multiple times for different instruments.")
    argp.add_argument("-k", "--skymap-file", help="Use skymap stored in given FITS file.")
    argp.add_argument("-x", "--coinc-xml", help="gstlal_inspiral XML file containing coincidence information.")
    argp.add_argument("-f", "--reference-freq", type=float, default=100.0, help="Waveform reference frequency. Required, default is 100 Hz.")
    argp.add_argument("-a", "--approximant", default="TaylorT4", help="Waveform family to use for templates. Any approximant implemented in LALSimulation is valid.")
    argp.add_argument("-A", "--amp-order", type=int, default=0, help="Include amplitude corrections in template waveforms up to this e.g. (e.g. 5 <==> 2.5PN), default is Newtonian order.")
    argp.add_argument("--l-max", type=int, default=2, help="Include all (l,m) modes with l less than or equal to this value.")
    argp.add_argument("-s", "--data-start-time", type=float, default=None, help="GPS start time of data segment. If given, must also give --data-end-time. If not given, sane start and end time will automatically be chosen.")
    argp.add_argument("-e", "--data-end-time", type=float, default=None, help="GPS end time of data segment. If given, must also give --data-start-time. If not given, sane start and end time will automatically be chosen.")
    argp.add_argument("-F", "--fmax", type=float, help="Upper frequency of signal integration. Default is use PSD's maximum frequency.")
    argp.add_argument("-t", "--event-time", type=float, help="GPS time of the event --- probably the end time. Required if --coinc-xml not given.")
    argp.add_argument("-i", "--inv-spec-trunc-time", type=float, default=8., help="Timescale of inverse spectrum truncation in seconds (Default is 8 - give 0 for no truncation)")
    argp.add_argument("-w", "--window-shape", type=float, default=0, help="Shape of Tukey window to apply to data (default is no windowing)")
    return argp

def add_output_params(argp):
    argp.add_argument("-o", "--output-file", help="Save result to this file.")
    argp.add_argument("-S", "--save-samples", type=int, default=0, help="Save this number of sample points to output-file. Requires --output-file to be defined.")
    argp.add_argument("-L", "--save-deltalnL", type=float, default=None, help="Threshold on deltalnL for points preserved in output file.  Requires --output-file to be defined")
    argp.add_argument("-P", "--save-P", type=float,default=0, help="Threshold on cumulative probability for points preserved in output file.  Requires --output-file to be defined")
    return argp

#
# Add the integration options
#
def add_integration_params(argp):
    integration_params = argp.add_argument_group("Integration Parameters", "Control the integration with these options.")
    integration_params.add_argument("--distance-maximum", default=300.0, type=float, help="Override the maximum distance in the prior. Default is 300 Mpc.")
    integration_params.add_argument("-m", "--time-marginalization", action="store_true", help="Perform marginalization over time via direct numerical integration. Default is false.")
    # Default is actually None, but that tells the integrator to go forever or until n_eff is hit.
    integration_params.add_argument("--zero-noise", action="store_true", help="Do not use input data as noise. Use with --pin-to-sim to make an injection")
    integration_params.add_argument("--n-max", type=int, help="Total number of samples points to draw. If this number is hit before n_eff, then the integration will terminate. Default is 'infinite'.",default=None)
    integration_params.add_argument("--n-eff", type=int, default=100, help="Total number of effective samples points to calculate before the integration will terminate. Default is 100")
    integration_params.add_argument("--n-chunk", type=int, help="Chunk'.",default=100)
    integration_params.add_argument("--convergence-tests-on",default=False,action='store_true')
    integration_params.add_argument("--seed", type=int, help="Random seed to use. Default is to not seed the RNG.")
    integration_params.add_argument("--no-adapt", action="store_true", help="Turn off adaptive sampling. Adaptive sampling is on by default.")
    integration_params.add_argument("--adapt-weight-exponent", type=float, default=1.0, help="Exponent to use with weights (likelihood integrand) when doing adaptive sampling. Used in tandem with --adapt-floor-level to prevent overconvergence. Default is 1.0.")
    integration_params.add_argument("--adapt-floor-level", type=float, default=0.1, help="Floor to use with weights (likelihood integrand) when doing adaptive sampling. This is necessary to ensure the *sampling* prior is non zero during adaptive sampling and to prevent overconvergence. Default is 0.1 (no floor)")
    integration_params.add_argument("--interpolate-time", default=False,help="If using time marginalization, compute using a continuously-interpolated array. (Default=false)")
    integration_params.add_argument("--fmin-template", dest='fmin_template', type=float, default=40, help="Waveform starting frequency.  Default is 40 Hz.")
    return argp

#
# Add the intrinsic parameters
#
def add_intrinsic_params(argp):
    intrinsic_params = argp.add_argument_group("Intrinsic Parameters", "Intrinsic parameters (e.g component mass) to use.")
    intrinsic_params.add_argument("--pin-to-sim", help="Pin values to sim_inspiral table entry.")
    intrinsic_params.add_argument("--mass1", type=float, help="Value of first component mass, in solar masses. Required if not providing coinc tables.")
    intrinsic_params.add_argument("--mass2", type=float, help="Value of second component mass, in solar masses. Required if not providing coinc tables.")
    intrinsic_params.add_argument("--spin1z", type=float, help="Value of first component spin (aligned with angular momentum), dimensionless.")
    intrinsic_params.add_argument("--spin2z", type=float, help="Value of second component spin (aligned with angular momentum), dimensionless.")
    intrinsic_params.add_argument("--eff-lambda", type=float, help="Value of effective tidal parameter. Optional, ignored if not given.")
    intrinsic_params.add_argument("--deff-lambda", type=float, help="Value of second effective tidal parameter. Optional, ignored if not given.")
    return argp

def parse_param(popts):
    """
    Parse out the specification of the intrinsic space. Examples:

    >>> parse_param(["mass1=1.4", "mass2", "spin1z=-1.0,10"])
    {'mass1': 1.4, 'mass2': None, 'spin1z': (-1.0, 10.0)}
    """
    if popts is None:
        return {}, {}
    intr_prms, expand_prms = {}, {}
    for popt in popts:
        popt = popt.split("=")
        if len(popt) == 1:
            # Implicit expand in full parameter space -- not yet completely
            # implemented
            intr_prms[popt[0]] = None
        elif len(popt) == 2:
            popt[1] = popt[1].split(",")
            if len(popt[1]) == 1:
                # Fix intrinsic point
                intr_prms[popt[0]] = float(popt[1][0])
            else:
                expand_prms[popt[0]] = tuple(map(float, popt[1]))
    return intr_prms, expand_prms


#
# DAG workflow related dictionaries
#

#
# Default priorities for categories within subdags
#
JOB_PRIORITIES = { "ILE": 10,
    "SQL": 1,
    "PLOT": 1
}

MAXJOBS = { 
    # Subdag maxjobs
    # ILE runs: breadth first (high priority) but ensure that no more than
    # this many *events* are in the queue
    "ANALYSIS": 20,
    # Postprocessing DAGs
    "POST": 4,

    # Intra DAG maxjobs settings
    # sqlite conversion
    "SQL": 10,
    # plotting jobs
    "PLOT": 10,

    # analysis jobs, no throttle
    # NOTE: This is disabled in favor of having the ILE subdag be parented by
    # the postprocessing
    "ILE": None
}

#
# Set up parameters and bounds
#

# FIXME: These need to be imported in ILE
t_ref_wind = 50e-3 # Interpolate in a window +/- this width about event time. 
dmin = 1.    # min distance
dmax = 300.  # max distance FOR ANY SOURCE EVER. EUCLIDEAN
distRef = 100*1e6*lal.PC_SI # a fiducial distance for the template source.

param_limits = { "psi": (0, 2*numpy.pi),
    "phi_orb": (0, 2*numpy.pi),
    "distance": (dmin, dmax),
    "right_ascension": (0, 2*numpy.pi),
    "declination": (-numpy.pi/2, numpy.pi/2),
    "t_ref": (-t_ref_wind, t_ref_wind),
    "inclination": (0, numpy.pi),
    "lam_tilde": (0, 5000), # FIXME: Needs reference
    "dlam_tilde": (-500, 500)
}

