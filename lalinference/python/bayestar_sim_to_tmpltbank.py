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
from __future__ import division
__doc__ = """
Create a template bank from an inspinj-style injection file that has a
sngl_inspiral record for each unique set of intrinsic parameters
(mass1, mass2, spin1z, spin2z) described by the sim_inspiral rows in the
injection file.
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"


# Command line interface.
from optparse import Option, OptionParser
from lalinference.bayestar import command

parser = OptionParser(
    formatter = command.NewlinePreservingHelpFormatter(),
    description = __doc__,
    usage="%prog [options] [INPUT.xml[.gz]] [-o OUTPUT.xml[.gz]]",
    option_list = [
        Option("-o", "--output", metavar="OUTPUT.xml[.gz]", default="/dev/stdout",
            help="Name of output file [default: %default]"),
        Option("--f-low", type=float,
            help="Override low frequency cutoff found in sim_inspiral table.")
    ]
)
opts, args = parser.parse_args()
infilename = command.get_input_filename(parser, args)


# Python standard library imports.
import os
from collections import namedtuple

# LIGO-LW XML imports.
from glue.ligolw import ligolw
from glue.ligolw.utils import process as ligolw_process
from glue.ligolw import table as ligolw_table
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import lsctables

# glue and LAL imports.
import lalsimulation

# BAYESTAR imports.
from lalinference.bayestar import timing


_fields = 'mass1 mass2 mchirp eta spin1z spin2z'

class TaylorF2RedSpinIntrinsicParams(namedtuple('intrinsic_params', _fields)):
    """Immutable container for just the intrinsic parameters belonging to a
    sim_inspiral record."""

    def __new__(cls, sim_inspiral):
        for attr in 'spin1x spin1y spin2x spin2y'.split():
            if getattr(sim_inspiral, attr):
                raise NotImplementedError('sim_inspiral:{} column is nonzero,'
                'but only aligned-spin templates are supported'.format(attr))
        return super(TaylorF2RedSpinIntrinsicParams, cls).__new__(cls,
            *(getattr(sim_inspiral, field) for field in cls._fields)
        )

    @property
    def chi(self):
        return lalsimulation.SimInspiralTaylorF2ReducedSpinComputeChi(
            self.mass1, self.mass2, self.spin1z, self.spin2z)


# Read injection file.
xmldoc = ligolw_utils.load_filename(infilename)

# Extract simulation table from injection file.
sim_inspiral_table = ligolw_table.get_table(xmldoc,
    lsctables.SimInspiralTable.tableName)

# Get just the intrinsic parameters from the sim_inspiral table.
sim_inspiral_intrinsic_params = set(TaylorF2RedSpinIntrinsicParams(sim_inspiral)
    for sim_inspiral in sim_inspiral_table)

if opts.f_low is None:
    # Get the low-frequency cutoffs from the sim_inspiral table.
    f_lows = set(sim_inspiral.f_lower for sim_inspiral in sim_inspiral_table)

    # There can be only one!
    try:
        f_low, = f_lows
    except ValueError:
        raise ValueError(
            "sim_inspiral:f_lower columns are not unique, got values: "
            + ' '.join(f_lows))
else:
    f_low = opts.f_low

# Open output file.
out_xmldoc = ligolw.Document()
out_xmldoc.appendChild(ligolw.LIGO_LW())

# Write process metadata to output file. Masquerade as lalapps_tmpltbank and
# encode low frequency cutoff in command line arguments.
process = ligolw_process.register_to_xmldoc(out_xmldoc, "tmpltbank",
    dict(opts.__dict__, low_frequency_cutoff=f_low), ifos="H1",
    comment="Exact-match template bank")

# Record low-frequency cutoff in the SearchSummVars table.
search_summvars_table = lsctables.New(lsctables.SearchSummVarsTable)
out_xmldoc.childNodes[0].appendChild(search_summvars_table)
search_summvars = lsctables.SearchSummVars()
search_summvars.search_summvar_id = search_summvars_table.get_next_id()
search_summvars.process_id = process.process_id
search_summvars.name = "low-frequency cutoff"
search_summvars.string = None
search_summvars.value = f_low
search_summvars_table.append(search_summvars)

# Create a SnglInspiral table and initialize its row ID counter.
sngl_inspiral_table = lsctables.New(lsctables.SnglInspiralTable)
out_xmldoc.childNodes[0].appendChild(sngl_inspiral_table)
sngl_inspiral_table.set_next_id(lsctables.SnglInspiralID(0))

for sim_inspiral_intrinsic_param in sim_inspiral_intrinsic_params:

    # Create new sngl_inspiral row and initialize its columns to None,
    # which produces an empty field in the XML output.
    sngl_inspiral = lsctables.SnglInspiral()
    for validcolumn in sngl_inspiral_table.validcolumns.iterkeys():
        setattr(sngl_inspiral, validcolumn, None)

    # Populate the row's fields.
    sngl_inspiral.event_id = sngl_inspiral_table.get_next_id()
    sngl_inspiral.mass1 = sim_inspiral_intrinsic_param.mass1
    sngl_inspiral.mass2 = sim_inspiral_intrinsic_param.mass2
    sngl_inspiral.mtotal = sngl_inspiral.mass1 + sngl_inspiral.mass2
    sngl_inspiral.mchirp = sim_inspiral_intrinsic_param.mchirp
    sngl_inspiral.eta = sim_inspiral_intrinsic_param.eta
    sngl_inspiral.f_final = timing.get_f_lso(sngl_inspiral.mass1, sngl_inspiral.mass2)
    sngl_inspiral.chi = sim_inspiral_intrinsic_param.chi

    # Add the row to the table in the document.
    sngl_inspiral_table.append(sngl_inspiral)

# Record process end time.
ligolw_process.set_process_end_time(process)

# Write output file.
ligolw_utils.write_filename(out_xmldoc, opts.output,
    gz=(os.path.splitext(opts.output)[-1]==".gz"))
