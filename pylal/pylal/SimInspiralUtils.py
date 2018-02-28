# Copyright (C) 2006  Duncan A. Brown
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
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw import ligolw
#
# =============================================================================
#
#                                   Input
#
# =============================================================================
#

class ExtractSimInspiralTableLIGOLWContentHandler(ligolw.PartialLIGOLWContentHandler):
  """
  LIGOLWContentHandler that will extract only the SimInspiralTable from a document.
  See glue.ligolw.LIGOLWContentHandler help for more info.
  """
  def __init__(self,document):
    def filterfunc(name,attrs):
      if name==ligolw.Table.tagName and attrs.has_key('Name'):
        return 0==table.CompareTableNames(attrs.get('Name'), lsctables.SimInspiralTable.tableName)
      else:
        return False
    ligolw.PartialLIGOLWContentHandler.__init__(self,document,filterfunc)


def ReadSimInspiralFromFiles(fileList, verbose=False):
  """
  Read the simInspiral tables from a list of files

  @param fileList: list of input files
  @param verbose: print ligolw_add progress
  """
  simInspiralTriggers = None

  lsctables.use_in(ExtractSimInspiralTableLIGOLWContentHandler)
  for thisFile in fileList:
    doc = utils.load_filename(thisFile, gz=(thisFile or "stdin").endswith(".gz"), verbose=verbose, contenthandler=ExtractSimInspiralTableLIGOLWContentHandler)
    # extract the sim inspiral table
    try: simInspiralTable = \
      table.get_table(doc, lsctables.SimInspiralTable.tableName)
    except: simInspiralTable = None
    if simInspiralTriggers and simInspiralTable: 
      simInspiralTriggers.extend(simInspiralTable)
    elif not simInspiralTriggers:
      simInspiralTriggers = simInspiralTable

  return simInspiralTriggers

