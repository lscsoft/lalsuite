#
# Copyright (C) 2013  Leo Singer
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
#
"""
Create a template bank that samples a regular lattice in (tau0, tau3), starting
from an initial (mass1, mass2), with lattice points spaced according to the
metric at the initial point.
"""
from __future__ import division
__author__ = "Leo Singer <leo.singer@ligo.org>"


# Command line interface.
from optparse import Option, OptionParser
from lalinference.bayestar import command

parser = OptionParser(
    formatter = command.NewlinePreservingHelpFormatter(),
    description = __doc__,
    usage="%prog [options] [INPUT.xml[.gz]] [-o OUTPUT.xml[.gz]]",
    option_list = [
        Option("--initial-mass1", metavar="Msun", type=float, default=1.4,
            help="Mass of first component of an initial lattice point [default: %default]"),
        Option("--initial-mass2", metavar="Msun", type=float, default=1.4,
            help="Mass of second component of an initial lattice point [default: %default]"),
        Option("--min-mass", metavar="Msun", type=float, default=1.,
            help="Minimum component mass [default: %default]"),
        Option("--max-mass", metavar="Msun", type=float, default=3.,
            help="Maximum component mass [default: %default]"),
        Option("-o", "--output", metavar="OUTPUT.xml[.gz]", default="/dev/stdout",
            help="Name of output file (default=stdout)")
    ]
)
opts, args = parser.parse_args()
infilename = command.get_input_filename(parser, args)


# Python standard library imports.
import os

# LIGO-LW XML imports.
from glue.ligolw import ligolw
from glue.ligolw.utils import process as ligolw_process
from glue.ligolw import table as ligolw_table
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import lsctables

# glue and LAL imports.
import lal
import lalinspiral.sbank.tau0tau3
import lalsimulation

# BAYESTAR imports.
from lalinference.bayestar import timing

# Other imports.
import numpy as np
import scipy.linalg


# Open output file.
xmldoc = ligolw.Document()
xmldoc.appendChild(ligolw.LIGO_LW())

# Write process metadata to output file.
process = ligolw_process.register_to_xmldoc(xmldoc, parser.get_prog_name(),
    opts.__dict__)

# Create a SnglInspiral table and initialize its row ID counter.
sngl_inspiral_table = lsctables.New(lsctables.SnglInspiralTable)
xmldoc.childNodes[0].appendChild(sngl_inspiral_table)
sngl_inspiral_table.set_next_id(lsctables.SnglInspiralID(0))

f_low = 10.
f_high = 2048.
df = 0.1

initial_mchirp = lalinspiral.sbank.tau0tau3.m1m2_to_mchirp(opts.initial_mass1, opts.initial_mass2)
initial_eta = opts.initial_mass1 * opts.initial_mass2 / (opts.initial_mass1 + opts.initial_mass2)**2
initial_chi = 0.
initial_chirp_times = lalsimulation.SimInspiralTaylorF2RedSpinChirpTimesFromMchirpEtaChi(initial_mchirp, initial_eta, initial_chi, f_low)
initial_theta0_theta3 = initial_chirp_times[:2]

# Sampled PSD.
S = lal.CreateREAL8Vector(int(f_high // df))
S.data = [lalsimulation.SimNoisePSDaLIGOZeroDetHighPower(i * df)
    for i in range(len(S.data))]

# Allocate noise moments.
moments = [lal.CreateREAL8Vector(int((f_high - f_low) // df)) for _ in range(29)]

# Compute noise moments.
lalsimulation.SimInspiralTaylorF2RedSpinComputeNoiseMoments(
    *(moments + [S, f_low, df]))

# Compute Fisher information matrix. Note factor of 2:
# The function SimInspiralTaylorF2RedSpinFisherMatrixChirpTimes
# returns the Fisher matrix for an SNR of 1/sqrt(2). The factor
# of 2 makes this the Fisher matrix at unit SNR.
I = lalsimulation.SimInspiralTaylorF2RedSpinFisherMatrixChirpTimes(
    *(initial_chirp_times + [f_low, df] + moments)).data * 2

# Blockwise separation of Fisher matrix. Parameters are in the following order:
# theta0, theta3, theta3S, t0, phi0
IA = I[0:2, 0:2] # (theta0, theta3) block
IB = I[0:2, 3:5] # cross block
ID = I[3:5, 3:5] # (time, phase) block

# Metric. We are dropping the theta3S terms completely and projecting out
# time and phase.
metric = IA - np.dot(IB, scipy.linalg.solve(ID, IB.T, sym_pos=True))

# Eigendecomposition of metric
metric_eigenvalues, metric_eigenvectors = np.linalg.eigh(metric)

# Shift between adjacent lattice points.
# FIXME: square root or no?
delta_theta0_theta3 = np.dot(metric_eigenvectors, np.diag(1 / np.sqrt(metric_eigenvalues)))

# FIXME: Determine appropriate boundaries to looping over lots of points that
# we are going to skip.
#
# T. Cokelaer (2007, http://dx.doi.org/10.1103/PhysRevD.76.102004) describes
# relationships between the component mass limits and the (tau0, tau3)
# boundaries.
n = 800
i0, i1 = np.mgrid[-n:n+1, -n:n+1]
i = np.column_stack((i0.flatten(), i1.flatten()))

# FIXME: Come up with a more natural way to specify the template spacing.
skip = 10
theta0_theta3 = np.dot(i, skip * delta_theta0_theta3.T) + initial_theta0_theta3

for th0, th3 in theta0_theta3:

    th3S = 0
    tau0 = th0 / (2 * np.pi * f_low)
    tau3 = -th3 / (2 * np.pi * f_low)

    mchirp, eta, chi = lalsimulation.SimInspiralTaylorF2RedSpinMchirpEtaChiFromChirpTimes(
        th0, th3, th3S, f_low
    )

    # Skip if either mchirp, eta, or chi are unphysical, unless this is the
    # initial point, which may be slightly unphysical just due to roundoff
    if np.all([th0, th3] == initial_theta0_theta3):
        mchirp = initial_mchirp
        eta = initial_eta
        mass1 = opts.initial_mass1
        mass2 = opts.initial_mass2
        mtotal = mass1 + mass2
    elif not (mchirp >= 0 and 0 <= eta <= 0.25 and -1 <= chi <= 1):
        continue
    else:
        mtotal = mchirp * eta ** -0.6
        mass1, mass2 = sorted(np.roots([1, -mtotal, eta * mtotal**2]))

    # Skip this one unless both component masses are in the appropriate range.
    if not(opts.min_mass <= mass1 <= opts.max_mass and opts.min_mass <= mass2 <= opts.max_mass):
        continue

    # Create new sngl_inspiral row and initialize its columns to None,
    # which produces an empty field in the XML output.
    sngl_inspiral = lsctables.SnglInspiral()
    for validcolumn in sngl_inspiral_table.validcolumns.iterkeys():
        setattr(sngl_inspiral, validcolumn, None)

    # Populate the row's fields.
    sngl_inspiral.event_id = sngl_inspiral_table.get_next_id()
    sngl_inspiral.mass1 = mass1
    sngl_inspiral.mass2 = mass2
    sngl_inspiral.tau0 = tau0
    sngl_inspiral.tau3 = tau3
    sngl_inspiral.mtotal = mtotal
    sngl_inspiral.mchirp = mchirp
    sngl_inspiral.eta = eta
    sngl_inspiral.chi = chi
    sngl_inspiral.f_final = timing.get_f_lso(mass1, mass2)

    # Add the row to the table in the document.
    sngl_inspiral_table.append(sngl_inspiral)

# Record process end time.
ligolw_process.set_process_end_time(process)

# Write output file.
ligolw_utils.write_filename(xmldoc, opts.output,
    gz=(os.path.splitext(opts.output)[-1]==".gz"))
