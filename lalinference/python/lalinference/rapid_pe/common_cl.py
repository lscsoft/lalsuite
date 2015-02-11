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

from optparse import OptionParser, OptionGroup

import numpy

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
def add_pinnable_params(optp, include=None, exclude=None):
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

    pinnable = OptionGroup(optp, "Pinnable Parameters", "Specifying these command line options will pin the value of that parameter to the specified value with a probability of unity.")
    for pin_param in pinned_params:
        option = "--" + pin_param.replace("_", "-")
        pinnable.add_option(option, type=float, help="Pin the value of %s." % pin_param)
    optp.add_option_group(pinnable)
    return optp

#
# Data source options for command line processing
#
def add_datasource_params(optp):
    optp.add_option("-c", "--cache-file", default=None, help="LIGO cache file containing all data needed.")
    optp.add_option("-C", "--channel-name", action="append", help="instrument=channel-name, e.g. H1=FAKE-STRAIN. Can be given multiple times for different instruments.")
    optp.add_option("-p", "--psd-file", action="append", help="instrument=psd-file, e.g. H1=H1_PSD.xml.gz. Can be given multiple times for different instruments.")
    optp.add_option("-k", "--skymap-file", help="Use skymap stored in given FITS file.")
    optp.add_option("-x", "--coinc-xml", help="gstlal_inspiral XML file containing coincidence information.")
    optp.add_option("-f", "--reference-freq", type=float, default=100.0, help="Waveform reference frequency. Required, default is 100 Hz.")
    optp.add_option("-a", "--approximant", default="TaylorT4", help="Waveform family to use for templates. Any approximant implemented in LALSimulation is valid.")
    optp.add_option("-A", "--amp-order", type=int, default=0, help="Include amplitude corrections in template waveforms up to this e.g. (e.g. 5 <==> 2.5PN), default is Newtonian order.")
    optp.add_option("--l-max", type=int, default=2, help="Include all (l,m) modes with l less than or equal to this value.")
    optp.add_option("-s", "--data-start-time", type=float, default=None, help="GPS start time of data segment. If given, must also give --data-end-time. If not given, sane start and end time will automatically be chosen.")
    optp.add_option("-e", "--data-end-time", type=float, default=None, help="GPS end time of data segment. If given, must also give --data-start-time. If not given, sane start and end time will automatically be chosen.")
    optp.add_option("-F", "--fmax", type=float, help="Upper frequency of signal integration. Default is use PSD's maximum frequency.")
    optp.add_option("-t", "--event-time", type=float, help="GPS time of the event --- probably the end time. Required if --coinc-xml not given.")
    optp.add_option("-i", "--inv-spec-trunc-time", type=float, default=8., help="Timescale of inverse spectrum truncation in seconds (Default is 8 - give 0 for no truncation)")
    optp.add_option("-w", "--window-shape", type=float, default=0, help="Shape of Tukey window to apply to data (default is no windowing)")
    return optp

def add_output_params(optp):
    optp.add_option("-o", "--output-file", help="Save result to this file.")
    optp.add_option("-S", "--save-samples", action="store_true", help="Save sample points to output-file. Requires --output-file to be defined.")
    optp.add_option("-L", "--save-deltalnL", type=float, default=float("Inf"), help="Threshold on deltalnL for points preserved in output file.  Requires --output-file to be defined")
    optp.add_option("-P", "--save-P", type=float,default=0, help="Threshold on cumulative probability for points preserved in output file.  Requires --output-file to be defined")
    return optp

#
# Add the integration options
#
def add_integration_params(optp):
    integration_params = OptionGroup(optp, "Integration Parameters", "Control the integration with these options.")
    integration_params.add_option("-m", "--time-marginalization", action="store_true", help="Perform marginalization over time via direct numerical integration. Default is false.")
    # Default is actually None, but that tells the integrator to go forever or until n_eff is hit.
    integration_params.add_option("--n-max", type=int, help="Total number of samples points to draw. If this number is hit before n_eff, then the integration will terminate. Default is 'infinite'.",default=1e7)
    integration_params.add_option("--n-eff", type=int, default=100, help="Total number of effective samples points to calculate before the integration will terminate. Default is 100")
    integration_params.add_option("--n-chunk", type=int, help="Chunk'.",default=100)
    integration_params.add_option("--convergence-tests-on",default=False,action='store_true')
    integration_params.add_option("--seed", type=int, help="Random seed to use. Default is to not seed the RNG.")
    integration_params.add_option("--no-adapt", action="store_true", help="Turn off adaptive sampling. Adaptive sampling is on by default.")
    integration_params.add_option("--adapt-weight-exponent", type=float, default=1.0, help="Exponent to use with weights (likelihood integrand) when doing adaptive sampling. Used in tandem with --adapt-floor-level to prevent overconvergence. Default is 1.0.")
    integration_params.add_option("--adapt-floor-level", type=float, default=0.1, help="Floor to use with weights (likelihood integrand) when doing adaptive sampling. This is necessary to ensure the *sampling* prior is non zero during adaptive sampling and to prevent overconvergence. Default is 0.1 (no floor)")
    integration_params.add_option("--interpolate-time", default=False,help="If using time marginalization, compute using a continuously-interpolated array. (Default=false)")
    optp.add_option_group(integration_params)
    return optp

#
# Add the intrinsic parameters
#
def add_intrinsic_params(optp):
    intrinsic_params = OptionGroup(optp, "Intrinsic Parameters", "Intrinsic parameters (e.g component mass) to use.")
    intrinsic_params.add_option("--pin-to-sim", help="Pin values to sim_inspiral table entry.")
    intrinsic_params.add_option("--mass1", type=float, help="Value of first component mass, in solar masses. Required if not providing coinc tables.")
    intrinsic_params.add_option("--mass2", type=float, help="Value of second component mass, in solar masses. Required if not providing coinc tables.")
    optp.add_option_group(intrinsic_params)
    return optp

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

param_limits = { "psi": (0, 2*numpy.pi),
    "phi_orb": (0, 2*numpy.pi),
    "distance": (dmin, dmax),
    "right_ascension": (0, 2*numpy.pi),
    "declination": (-numpy.pi/2, numpy.pi/2),
    "t_ref": (-t_ref_wind, t_ref_wind),
    "inclination": (0, numpy.pi)
}

