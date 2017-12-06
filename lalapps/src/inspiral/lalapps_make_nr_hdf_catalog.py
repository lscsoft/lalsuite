# Copyright (C) 2015/2017 Patricia Schmidt, Ian W. Harry
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
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
This program constructs an NR waveform catalog in xml format by reading the relevant NR metadata from the HDF5 file. 
It also holds the pointers to the NR data directories. Any catalog generated that way can directly be used by 
lalapps_inspinj.
"""

import sys, argparse

import h5py
import lal
from lalsimulation import SimInspiralNRWaveformGetSpinsFromHDF5File

from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import ilwd
from glue.ligolw import utils
from glue.ligolw.utils import process as ligolw_process

import lalapps.git_version

__author__  = "Patricia Schmidt <patricia.schmidt@ligo.org>, "
__author__  += "Ian Harry <ian.harry@ligo.org>"
__version__ = lalapps.git_version.verbose_msg
__date__    = lalapps.git_version.date
__program__ = "lalapps_make_nr_hdf_catalog"


cols = lsctables.SimInspiralTable.validcolumns

def fill_missing_columns(sim):
    for entry in cols.keys():
        if not(hasattr(sim,entry)):
            if cols[entry] in ['real_4','real_8']:
                setattr(sim,entry,0.)
            elif cols[entry] == 'int_4s':
                setattr(sim,entry,0)
            elif cols[entry] == 'lstring':
                setattr(sim,entry,'')
            elif entry == 'simulation_id' or entry == 'process_id':
                continue
            else:
                print >> sys.stderr, "Column %s not recognized" %(entry)
                raise ValueError

_desc = __doc__[1:]
parser = argparse.ArgumentParser(description=_desc)

parser.add_argument('--version', action='version',
                    version=lalapps.git_version.verbose_msg)
parser.add_argument("-V", "--verbose", action="store_true",
                    help="print extra debugging information", default=False )
parser.add_argument("-o", "--output-file", action="store", type=str,
                    required=True, help="Output file name")
parser.add_argument("-i", "--input-files", nargs= '*', dest='inputs',
                    action="store", type=str, required=True,
                    help="Path(s) to HDF5 input files")

args = parser.parse_args()

# Prepare xml document
xmldoc = ligolw.Document()
xmldoc.appendChild(ligolw.LIGO_LW())

proc_id = ligolw_process.register_to_xmldoc\
    (xmldoc, "nr_catalog", args.__dict__, comment="",
     version=lalapps.git_version.version,
     cvs_repository='lalsuite/' + lalapps.git_version.branch,
     cvs_entry_time=lalapps.git_version.date).process_id

sim_table = lsctables.New(lsctables.SimInspiralTable)

inj_list = args.inputs

for count, inj in enumerate(inj_list):
    curr_sim = lsctables.SimInspiral()
    # Add the empty columns
    fill_missing_columns(curr_sim)
    # Set id columns
    curr_sim.process_id = proc_id
    curr_sim.simulation_id = ilwd.ilwdchar("sim_inspiral:simulation_id:%d"\
                                           %(count))
    curr_sim.numrel_data = inj
    f = h5py.File(inj, 'r')  
    curr_sim.eta = f.attrs['eta']
    if curr_sim.eta > 0.25 and curr_sim.eta < 0.2501:
        curr_sim.eta = 0.25
    # Populate spins columns with spins in LAL frame! Need to be
    # transformed from NR frame
    spins = SimInspiralNRWaveformGetSpinsFromHDF5File(inj)
    curr_sim.spin1x = spins[0]
    curr_sim.spin1y = spins[1]
    curr_sim.spin1z = spins[2]
    curr_sim.spin2x = spins[3]
    curr_sim.spin2y = spins[4]
    curr_sim.spin2z = spins[5]

    curr_sim.f_lower = f.attrs['f_lower_at_1MSUN']
    f.close()
    # These were the old columns used to specify min and max *l* modes in NINJA
    # not using them here as they ignore *m* modes.
    #curr_sim.numrel_mode_max = 0
    #curr_sim.numrel_mode_min = 0
    
    sim_table.append(curr_sim)

xmldoc.childNodes[-1].appendChild(sim_table)
utils.write_filename(xmldoc, args.output_file,
                     gz=args.output_file.endswith('gz'))
