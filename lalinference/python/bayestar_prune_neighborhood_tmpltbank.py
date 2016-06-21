#
# Copyright (C) 2013-2016  Leo Singer
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
Remove all entries from a template bank except those that lie within a 1-sigma
error ellipse of a (mass1, mass2, chi=0) at a given SNR. Uses
TaylorF2ReducedSpin metric.
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"


# Command line interface.
import argparse
from lalinference.bayestar import command

parser = command.ArgumentParser()
parser.add_argument("--mass1", type=float, metavar="MSUN", required=True)
parser.add_argument("--mass2", type=float, metavar="MSUN", required=True)
parser.add_argument("--snr", type=float, metavar="SNR", required=True)
parser.add_argument(
    'input', metavar='IN.xml[.gz]', type=argparse.FileType('rb'),
    default='-', help='Name of input file [default: stdin]')
parser.add_argument(
    '-o', '--output', metavar='OUT.xml[.gz]', type=argparse.FileType('wb'),
    default='-', help='Name of output file [default: stdout]')
opts = parser.parse_args()


# Python standard library imports.
import os

# LIGO-LW XML imports.
from glue.ligolw.utils import process as ligolw_process
from glue.ligolw import table as ligolw_table
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import lsctables

# lal imports.
import lal
import lalinspiral.sbank.tau0tau3
import lalsimulation

# Other imports.
import numpy as np
from scipy import linalg

# BAYESTAR imports.
from lalinference.bayestar import ligolw as ligolw_bayestar


# Read input file.
xmldoc, _ = ligolw_utils.load_fileobj(
    opts.input, contenthandler=ligolw_bayestar.LSCTablesContentHandler)

# Write process metadata to output file.
process = command.register_to_xmldoc(xmldoc, parser, opts)

# Determine the low frequency cutoff from the template bank file.
f_low = ligolw_bayestar.get_template_bank_f_low(xmldoc)

# Get the SnglInspiral table.
sngl_inspiral_table = ligolw_table.get_table(xmldoc,
    lsctables.SnglInspiralTable.tableName)

# Determine central values of intrinsic parameters.
mchirp0 = lalinspiral.sbank.tau0tau3.m1m2_to_mchirp(opts.mass1, opts.mass2)
eta0 = opts.mass1 * opts.mass2 / np.square(opts.mass1 + opts.mass2)
chi0 = 0.

# Transform to chirp times.
thetas_0 = lalsimulation.SimInspiralTaylorF2RedSpinChirpTimesFromMchirpEtaChi(
    mchirp0, eta0, chi0, f_low)

# Sampled PSD.
f_high = 2048
df = 0.1
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
    *(thetas_0 + [f_low, df] + moments)).data * 2

# Now amplify to the requested SNR.
I *= np.square(opts.snr)

# Blockwise separation of Fisher matrix. Parameters are in the following order:
# theta0, theta3, theta3S, t0, phi0
IA = I[0:3, 0:3] # intrinsic block
IB = I[0:3, 3:5] # cross block
ID = I[3:5, 3:5] # extrinsic block
metric = IA - np.dot(IB, linalg.solve(ID, IB.T, sym_pos=True))


def predicate(sngl):
    """Return True if a template is within a 1-sigma radius of the central
    template, False otherwise"""
    thetas = lalsimulation.SimInspiralTaylorF2RedSpinChirpTimesFromMchirpEtaChi(
        sngl.mchirp, sngl.eta, sngl.chi, f_low)
    dtheta = np.asarray(thetas) - thetas_0
    distance = np.dot(dtheta, np.dot(metric, dtheta))
    return distance <= 1

# Grab the templates that are at most 1 sigma from the central (mass1, mass2).
rows_to_keep = filter(predicate, sngl_inspiral_table)
del sngl_inspiral_table[:]
sngl_inspiral_table.extend(rows_to_keep)


# Record process end time.
ligolw_process.set_process_end_time(process)

# Write output.
ligolw_utils.write_fileobj(xmldoc, opts.output,
    gz=(os.path.splitext(opts.output.name)[-1] == '.gz'))
