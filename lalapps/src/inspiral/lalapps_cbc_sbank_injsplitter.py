# Copyright (C) 2015  Surabhi Sachdev, Tjonnie Li
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

###############################################################################
#
# IMPORT MODULES
#
###############################################################################

from time import strftime
from operator import attrgetter
from optparse import OptionParser
from pylal import spawaveform
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw import ligolw
from glue.ligolw.utils import process as ligolw_process
from pylal.datatypes import LIGOTimeGPS
import numpy

###############################################################################
#
# IMPORT MODULES
#
###############################################################################

class ContentHandler(ligolw.LIGOLWContentHandler):
    pass
lsctables.use_in(ContentHandler)

###############################################################################
#
# COMMAND LINE PARSING
#
###############################################################################

def parse_command_line():
    parser = OptionParser()
    parser.add_option("-o", "--output-path", metavar = "path", default = ".", help = "Set the path to the directory where output files will be written.  Default is \".\".")
    parser.add_option("-u", "--usertag", metavar = "usertag", help = "Number of splits",default="INJSPLITTER")
    parser.add_option("-n", "--nsplit", metavar = "count", type = "int", help = "Number of splits",default=1)
    parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
    options, filenames = parser.parse_args()

    if not filenames:
        raise ValueError, "must provide list of filenames"

    return options, filenames

options, filenames = parse_command_line()

###############################################################################
#
# MAIN CODE
#
###############################################################################

# GETTING COMMAND LINE OPTIONS FOR PRINTING INTO THE TABLE
opts_dict = dict((k, v) for k, v in options.__dict__.iteritems() if v is not False and v is not None)

# LOAD INJECTION TABLE
xmldoc=utils.load_filename(filenames[0], gz=filenames[0].endswith(".gz"), verbose = options.verbose, contenthandler=ContentHandler)
sim_inspiral_table=lsctables.table.get_table(xmldoc, lsctables.SimInspiralTable.tableName)
process_params_table = lsctables.table.get_table(xmldoc, lsctables.ProcessParamsTable.tableName)

# PREPARE PROCESS TABLE WITH INFORMATION ABOUT THE CURRENT PROGRAM
process = ligolw_process.register_to_xmldoc(xmldoc,
"lalapps_cbc_injsplitter", opts_dict,
version="no version", cvs_repository="sbank",
cvs_entry_time=strftime('%Y/%m/%d %H:%M:%S'))

sim_inspiral_table_split = lsctables.table.new_from_template(sim_inspiral_table)
sim_inspiral_table.parentNode.replaceChild(sim_inspiral_table_split, sim_inspiral_table)

# CALCULATE THE NUMBER OF TEMPLATES IN EACH SPLIT. SPLIT AS EVENLY AS POSSIBLE
#ninj = len(sim_inspiral_table);
#nquo,nrem = divmod(ninj,options.nsplit);
#nperbank = numpy.add(numpy.ones(nsplit)*nquo,1.0) ;
evensplit = numpy.array_split(sim_inspiral_table,options.nsplit);
for i in xrange(options.nsplit):
	sim_inspiral_table_split[:] = evensplit[i]
	ligolw_process.set_process_end_time(process)
	utils.write_filename(xmldoc, "%s_INJ_SPLIT_%04d.xml"%(options.usertag,i), gz = False, verbose = options.verbose)

#for i in xrange(options.nsplit):
	#first_row = i*ninjSub;
	#if i == options.nsplit-1:
		#last_row = length;
	#else:
		#last_row = first_row + ninjSub;
	#sim_inspiral_table_split[:] = sim_inspiral_table[first_row:last_row]
	#ligolw_process.set_process_end_time(process)
	#utils.write_filename(xmldoc, "%s_INJ_SPLIT_%04d.xml"%(options.usertag,i), gz = False, verbose = options.verbose)
