#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       cbcBayesPosToSimInspiral.py
#
#       Copyright 2013
#       Benjamin Farr <bfarr@u.northwestern.edu>,
#       Will M. Farr <will.farr@ligo.org>,
#       John Veitch <john.veitch@ligo.org>
#
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.
"""
Populate a sim_inspiral table with random draws from an ASCII table.
"""
from optparse import Option, OptionParser
import numpy as np
from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import ilwd
import matplotlib
matplotlib.use("Agg") # Needed to run on the CIT cluster
from pylal import bayespputils as bppu

# Create a datatype for all relavent fields to be filled in the sim_inspiral table
sim_inspiral_dt = [
        ('waveform','|S64'),
        ('taper','|S64'),
        ('f_lower', 'f8'),
        ('mchirp', 'f8'),
        ('eta', 'f8'),
        ('mass1', 'f8'),
        ('mass2', 'f8'),
        ('geocent_end_time', 'f8'),
        ('geocent_end_time_ns', 'f8'),
        ('distance', 'f8'),
        ('longitude', 'f8'),
        ('latitude', 'f8'),
        ('inclination', 'f8'),
        ('coa_phase', 'f8'),
        ('polarization', 'f8'),
        ('spin1x', 'f8'),
        ('spin1y', 'f8'),
        ('spin1z', 'f8'),
        ('spin2x', 'f8'),
        ('spin2y', 'f8'),
        ('spin2z', 'f8'),
        ('amp_order', 'i4'),
        ('numrel_data','|S64')
]

def get_input_filename(parser, args):
    """Determine name of input: either the sole positional command line argument,
    or /dev/stdin."""
    if len(args) == 0:
        infilename = '/dev/stdin'
    elif len(args) == 1:
        infilename = args[0]
    else:
        parser.error("Too many command line arguments.")
    return infilename

def standardize_param_name(params, possible_names, desired_name):
    for name in possible_names:
        if name in params: params[params.index(name)] = desired_name

def standardize_param_names(params):
    standardize_param_name(params, ['m1'], 'm1')
    standardize_param_name(params, ['m2'], 'm2')
    standardize_param_name(params, ['mc', 'chirpmass'], 'mchirp')
    standardize_param_name(params, ['massratio'], 'eta')
    standardize_param_name(params, ['d', 'dist'], 'distance')
    standardize_param_name(params, ['ra'], 'longitude')
    standardize_param_name(params, ['dec'], 'latitude')
    standardize_param_name(params, ['iota'], 'inclination')
    standardize_param_name(params, ['phi', 'phase', 'phi0'], 'phi_orb')
    standardize_param_name(params, ['psi', 'polarisation'], 'polarization')


def compute_mass_parameterizations(samples):
    params = samples.dtype.names
    has_mc = 'mchirp' in params
    has_eta = 'eta' in params
    has_q = 'q' in params
    has_ms = 'mass1' in params and 'mass2' in params
    has_mtotal = 'mtotal' in params

    if has_mc:
        mc = samples['mchirp']
        if not has_eta:
            if has_q:
                eta = bppu.q2eta(mc, samples['q'])
            else:
                raise ValueError("Chirp mass given with no mass ratio.")
        else:
            eta = samples['eta']

        if not has_ms:
            m1, m2 = bppu.mc2ms(mc, eta)

        mtotal = m1 + m2

    elif has_ms:
        m1 = samples['mass1']
        m2 = samples['mass2']
        mtotal = m1 + m2
        eta = m1 * m2 / (mtotal * mtotal)
        mc = mtotal * np.power(eta, 3./5.)

    elif has_mtotal:
        mtotal = samples['mtotal']
        if has_eta:
            eta = samples['eta']
            mc = mtotal * np.power(eta, 3./5.)
            m1, m2 = bppu.mc2ms(mc, eta)
        elif has_q:
            m1 = mtotal / (1 + samples['q'])
            m2 = mtotal - m1
        else:
            raise ValueError("Chirp mass given with no mass ratio.")
    return mc, eta, m1, m2, mtotal


if __name__ == "__main__":
    parser = OptionParser(
            description = __doc__,
            usage = "%prog [options] [INPUT]",
            option_list = [
                Option("-o", "--output", metavar="FILE.xml",
                    help="name of output XML file"),
                Option("--num-of-injs", metavar="NUM", type=int, default=200,
                    help="number of injections"),
                Option("--approx", metavar="APPROX", default="SpinTaylorT4threePointFivePN",
                    help="approximant to be injected"),
                Option("--taper", metavar="TAPER", default="TAPER_NONE",
                    help="Taper methods for injections"),
                Option("--flow", metavar="FLOW", type=float, default=None,
                    help="Taper methods for injections"),
                Option("--amporder", metavar="AMPORDER", type=int, default=0,
                    help="pN order in amplitude for injection")
            ]
    )

    opts, args = parser.parse_args()
    infilename = get_input_filename(parser, args)

    # Read in ASCII table, assuming column names are the first line
    with open(infilename, 'r') as inp:
        params = inp.readline().split()
        standardize_param_names(params)
        samples = np.loadtxt(inp, dtype=[(p, np.float) for p in params])

    N = opts.num_of_injs
    if len(samples) < N:
        raise ValueError("{} injections requested, but {} samples were provided.".format(N, len(samples)))

    # Choose subset for sim_inspiral_table
    selection = np.arange(len(samples))
    np.random.shuffle(selection)
    samples = samples[selection[:N]]

    # Create an empty structured array with names indentical to sim_inspiral fields
    injections = np.zeros((N,), dtype=sim_inspiral_dt)

    # Determine all mass parameterizations
    mc, eta, m1, m2, mtotal = compute_mass_parameterizations(samples)

    # Get cycle numbers as simulation_ids
    ids = range(N)

    # Compute cartesian spins
    if 'a1' in params and 'theta1' in params and 'phi1' in params:
        s1x, s1y, s1z = bppu.sph2cart(samples['a1'], samples['theta1'], samples['phi1'])
    elif 'a1z' in params:
        s1z = samples['a1z']
        s1x = np.zeros_like(s1z)
        s1y = np.zeros_like(s1z)
    else:
        s1x = np.zeros_like(m1)
        s1y = np.zeros_like(m1)
        s1z = np.zeros_like(m1)


    if 'a2' in params and 'theta2' in params and 'phi2' in params:
        s2x, s2y, s2z = bppu.sph2cart(samples['a2'], samples['theta2'], samples['phi2'])
    elif 'a2z' in params:
        s2z = samples['a2z']
        s2x = np.zeros_like(s2z)
        s2y = np.zeros_like(s2z)
    else:
        s2x = np.zeros_like(m2)
        s2y = np.zeros_like(m2)
        s2z = np.zeros_like(m2)

    system_frame_params = set([ \
            'costheta_jn', \
            'phi_jl', \
            'tilt1', 'tilt2', \
            'phi12', \
            'a1','a2', \
            'f_ref' \
    ])
    theta_jn=np.array([np.arccos(i) for i in samples['costheta_jn']])
    if set(params).intersection(system_frame_params) == system_frame_params:
        inclination, theta1, phi1, theta2, phi2, _ = bppu.physical2radiationFrame(
                theta_jn,
                samples['phi_jl'],
                samples['tilt1'],
                samples['tilt2'],
                samples['phi12'],
                samples['a1'],
                samples['a2'],
                m1, m2,
                samples['f_ref'])
        inclination = inclination.flatten()
        theta1 = theta1.flatten()
        phi1 = phi1.flatten()
        theta2 = theta2.flatten()
        phi2 = phi2.flatten()
        s1x, s1y, s1z = bppu.sph2cart(samples['a1'], theta1, phi1)
        s2x, s2y, s2z = bppu.sph2cart(samples['a2'], theta2, phi2)
    else:
        inclination = theta_jn

    print s1x.shape
    # Check if f_low is given on the command line. If not, try to take if from 'samples'.
    if opts.flow is None:
        try:
            flow = samples['flow']
        except:
            raise ValueError("No f_low found in input file or command line arguments.")
    else:
        try:
            samples['flow']
        except:
            pass
        else:  # executed if no exception is raised
            print('f_low given in both input file and command line.'
                  ' Using command line argument: %r' % opts.flow)
        flow = [opts.flow for i in xrange(N)]

    # Populate structured array
    injections['waveform'] = [opts.approx for i in xrange(N)]
    injections['taper'] = [opts.taper for i in xrange(N)]
    injections['f_lower'] = flow
    injections['mchirp'] = mc
    injections['eta'] = eta
    injections['mass1'] = m1
    injections['mass2'] = m2
    injections['geocent_end_time'] = np.modf(samples['time'])[1]
    injections['geocent_end_time_ns'] = np.modf(samples['time'])[0] * 10**9
    injections['distance'] = samples['distance']
    injections['longitude'] = samples['longitude']
    injections['latitude'] = samples['latitude']
    injections['inclination'] = inclination
    injections['coa_phase'] = samples['phi_orb']
    injections['polarization'] = samples['polarization']
    injections['spin1x'] = s1x
    injections['spin1y'] = s1y
    injections['spin1z'] = s1z
    injections['spin2x'] = s2x
    injections['spin2y'] = s2y
    injections['spin2z'] = s2z
    injections['amp_order'] = [opts.amporder for i in xrange(N)]
    injections['numrel_data'] = [ "" for _ in xrange(N)]

    # Create a new XML document
    xmldoc = ligolw.Document()
    xmldoc.appendChild(ligolw.LIGO_LW())
    sim_table = lsctables.New(lsctables.SimInspiralTable)
    xmldoc.childNodes[0].appendChild(sim_table)

    # Add empty rows to the sim_inspiral table
    for inj in xrange(N):
        row = sim_table.RowType()
        for slot in row.__slots__: setattr(row, slot, 0)
        sim_table.append(row)

    # Fill in IDs
    for i,row in enumerate(sim_table):
        row.process_id = ilwd.ilwdchar("process:process_id:{0:d}".format(i))
        row.simulation_id = ilwd.ilwdchar("sim_inspiral:simulation_id:{0:d}".format(ids[i]))

    # Fill rows
    for field in injections.dtype.names:
        vals = injections[field]
        for row, val in zip(sim_table, vals): setattr(row, field, val)

    # Write file
    output_file = open(opts.output, 'w')
    xmldoc.write(output_file)
    output_file.close()

