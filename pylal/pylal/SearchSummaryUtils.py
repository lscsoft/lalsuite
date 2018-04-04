# Copyright (C) 2006  Craig Robinson and Anand Sengupta
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

import gzip

from glue import segments

from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw.utils import search_summary as ligolw_search_summary

try:
  lsctables.use_in(ligolw.PartialLIGOLWContentHandler)
except AttributeError:
  # old glue did not allow .use_in().
  # FIXME:  remove when we can require the latest version of glue
  pass

def _table_filter_generator(tables):
  """
  Return a function suitable for a PartialLIGOLWContentHandler that filters
  for any table in the given list of tables.
  """
  def filter_func(name, attrs):
    for table in tables:
      if lsctables.IsTableProperties(table, name, attrs):
        return True
    return False
  return filter_func

def ReadTablesFromFiles(fileList, table_list, verbose=False):
  """
  Return a parsed document containing only the tables enumerated in
  table_list that appear in the files in fileList.
  """
  _is_table = _table_filter_generator(table_list)
  doc = ligolw.Document()
  searchsummary_handler = ligolw.PartialLIGOLWContentHandler(doc, _is_table)

  for thisFile in fileList:
    if thisFile.endswith(".gz"):
      fileobj = gzip.open(thisFile)
    else:
      fileobj = open(thisFile)
    ligolw.make_parser(searchsummary_handler).parse(fileobj)
  return doc

def GetSegListFromSearchSummaries(fileList, verbose=False):
  """
  Read segment lists from search summary tables
  @param fileList: list of input files.
  """
  required_tables = [lsctables.SearchSummaryTable, lsctables.ProcessTable]

  segList = segments.segmentlistdict()
  for thisFile in fileList:
    doc = ReadTablesFromFiles([thisFile], required_tables, verbose)
    try:
      segs = ligolw_search_summary.segmentlistdict_fromsearchsummary(doc)
    except:
      raise ValueError, "Cannot extract segments from the SearchSummaryTable of %s" % thisFile

    #Now add these segments to the existing list
    segList.extend(segs)

  for value in segList.values():
    value.sort()

  return segList

