#!/usr/bin/env python
# -*- coding: utf-8 -*-
#       cbcBayesPosToSimBurst.py
#       C2014 Salvatore Vitale <salvatore.vitale@ligo.org> 
#       
#       based on
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
from pylal import bayespputils as bppu

# Create a datatype for all relavent fields to be filled in the sim_inspiral table
sim_inspiral_dt = [
        ('waveform','|S64'),
        ('frequency', 'f8'),
        ('q', 'f8'),
        ('bandwidth', 'f8'),
        ('duration', 'f8'),
        ('hrss', 'f8'),
        ('time_geocent_gps', 'i4'),
        ('time_geocent_gps_ns', 'i4'),
        ('ra','f8'),
        ('dec', 'f8'),
        ('pol_ellipse_angle', 'f8'),
        ('pol_ellipse_e', 'f8'),
        ('psi', 'f8'),
        ('amplitude', 'f8'),
        ('egw_over_rsquared', 'f8'),
        ('waveform_number', 'i8'),
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
    standardize_param_name(params, ['frequency','freq',], 'frequency')
    standardize_param_name(params, ['q','Q','quality'], 'q')
    standardize_param_name(params, ['duration', 'tau'], 'duration')
    standardize_param_name(params, ['ra','longitude','rightascension'],'ra')
    standardize_param_name(params, ['dec','latitude','declination'],'dec')
    standardize_param_name(params, ['alpha','polar_angle','pol_ellipse_angle'], 'pol_ellipse_angle')
    standardize_param_name(params, ['polar_eccentricity','eccentricity','pol_ellipse_e'], 'pol_ellipse_e')
    standardize_param_name(params, ['psi', 'polarisation','polarization'],'psi')
    

def compute_duration_parameterizations(samples):
    from math import sqrt
    params = samples.dtype.names
    has_q = 'q' in params
    has_dur='duration' in params
    has_frequency = 'frequency' in params
    if has_q:
        q = samples['q']
        if has_frequency:
          duration=[i/sqrt(2)/np.pi/j for (i,j) in zip(q,samples['frequency'])]
        else:
          duration=np.nan
    elif has_dur:
        duration = samples['duration']
        if has_frequency:
          q=[sqrt(2)*np.pi*i*j for (i,j) in zip(duration,samples['frequency'])]
        else:
          q=np.nan
    return q,duration 


if __name__ == "__main__":
    parser = OptionParser(
            description = __doc__,
            usage = "%prog [options] [INPUT]",
            option_list = [
                Option("-o", "--output", metavar="FILE.xml",
                    help="name of output XML file"),
                Option("--num-of-injs", metavar="NUM", type=int, default=200,
                    help="number of injections"),
                Option("--approx", metavar="APPROX", default="SineGaussianF",
                    help="approximant to be injected"),
                Option("--taper", metavar="TAPER", default="TAPER_NONE",
                    help="Taper methods for injections"),
                Option("--flow", metavar="FLOW", type=float, default=None,
                    help="Taper methods for injections"),
            ]
    )

    opts, args = parser.parse_args()
    infilename = get_input_filename(parser, args)

    # Read in ASCII table, assuming column names are the first line
    with open(infilename, 'r') as inp:
        params = inp.readline().split()
        standardize_param_names(params)
        samples = np.loadtxt(inp, dtype=[(p, np.float) for p in params])

    # Choose subset for sim_inspiral_table
    N = opts.num_of_injs
    selection = np.arange(N)
    np.random.shuffle(selection)
    samples = samples[selection]

    # Create an empty structured array with names indentical to sim_inspiral fields
    injections = np.zeros((N,), dtype=sim_inspiral_dt)

    # Determine all mass parameterizations
    q, dur = compute_duration_parameterizations(samples)

    # Get cycle numbers as simulation_ids
    ids = range(N)

    # Compute polar parameters from alpha if given
    if 'alpha' in params:
      # LIB may use use the convention that pol_ellipse_angle=alpha).
        injections['pol_ellipse_angle'] = samples['alpha']
    elif 'pol_ellipse_angle' in params:
        injections['pol_ellipse_angle'] = samples['pol_ellipse_angle']
    else:
        injections['pol_ellipse_angle'] = None

    if 'polar_eccentricity' in params:
        injections['pol_ellipse_e']= samples['polar_eccentricity']
    elif 'pol_ellipse_e' in params:
        injections['pol_ellipse_e']= samples['pol_ellipse_e']
    else:
        injections['pol_ellipse_e']=None
    print params
    print samples['pol_ellipse_e']
    print injections['pol_ellipse_e']
    if 'bandwidth' in params:
      injections['bandwidth'] = samples['bandwidth']
    else:
      injections['bandwidth'] = np.nan
    if 'time' in params:
      injections['time_geocent_gps'] = np.modf(samples['time'])[1]
      injections['time_geocent_gps_ns'] = np.modf(samples['time'])[0] * 10**9
    elif 'time_min' in params and 'time_max' in params:
      est_time=samples['time_min']+0.5*(samples['time_max']-samples['time_min'])
      injections['time_geocent_gps'] = np.modf(est_time)[1]
      injections['time_geocent_gps_ns'] = np.modf(est_time)[0] * 10**9
      
    # Populate structured array
    injections['waveform'] = [opts.approx for i in xrange(N)]
    try:
      injections['frequency'] = samples['frequency']
    except:
      injections['frequency']=[np.nan for i in samples['ra']]
    injections['duration'] = dur
    injections['q'] = q
    try:
      injections['hrss'] = samples['hrss']
    except:
      injections['hrss'] = np.exp(samples['loghrss'])
    injections['ra'] = samples['ra']
    injections['dec'] = samples['dec']
    injections['psi'] = samples['psi']

    # Create a new XML document
    xmldoc = ligolw.Document()
    xmldoc.appendChild(ligolw.LIGO_LW())
    #create timeslide table and set offsets to 0
    timeslide_table = lsctables.New(lsctables.TimeSlideTable)
    p=lsctables.Process
    p.process_id=ilwd.ilwdchar("process:process_id:{0:d}".format(0))
    timeslide_table.append_offsetvector({'H1':0,'V1':0,'L1':0,'H2':0},p)

    sim_table = lsctables.New(lsctables.SimBurstTable)
    xmldoc.childNodes[0].appendChild(timeslide_table)
    xmldoc.childNodes[0].appendChild(sim_table)

    # Add empty rows to the sim_inspiral table
    for inj in xrange(N):
        row = sim_table.RowType()
        for slot in row.__slots__: setattr(row, slot, 0)
        sim_table.append(row)

    # Fill in IDs
    for i,row in enumerate(sim_table):
        row.process_id = ilwd.ilwdchar("process:process_id:{0:d}".format(i))
        row.simulation_id = ilwd.ilwdchar("sim_burst:simulation_id:{0:d}".format(ids[i]))
        row.time_slide_id = ilwd.ilwdchar("time_slide:time_slide_id:{0:d}".format(0))
    # Fill rows
    for field in injections.dtype.names:
        vals = injections[field]
        for row, val in zip(sim_table, vals): setattr(row, field, val)

    # Write file
    output_file = open(opts.output, 'w')
    xmldoc.write(output_file)
    output_file.close()

