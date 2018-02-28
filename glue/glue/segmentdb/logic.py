#
# Copyright (C) 2009  Larne Pekowsky
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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


"""
Utility methods for doing logical operations on sets of segments
"""


import sys
import os

import glue.segments

from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables

from glue.ligolw.utils import ligolw_add

from glue.segmentdb.segmentdb_utils import add_to_segment_definer
from glue.segmentdb.segmentdb_utils import add_to_segment
from glue.segmentdb.segmentdb_utils import add_to_segment_summary
from glue.segmentdb.segmentdb_utils import find_segments

from glue import git_version
__date__ = git_version.date
__version__ = git_version.id
__author__  = "Larne Pekowsky <lppekows@physics.syr.edu>"


# "enum" indicating operation
INTERSECT = 'intersect'
UNION     = 'union'
DIFF      = 'diff'



def run_file_operation(outdoc, filenames, use_segment_table, operation, preserve = True):
    """
    Performs an operation (intersect or union) across a set of files.
    That is, given a set of files each with segment definers DMT-FLAG1,
    DMT-FLAG2 etc the result is a file where 

    DMT-FLAG1 = (file 1's DMT-FLAG1 operation file 2's DMT-FLAG1 operation ...)
    DMT-FLAG2 = (file 1's DMT-FLAG2 operation file 2's DMT-FLAG2 operation ...)
    
    etc
    """

    proc_id = table.get_table(outdoc, lsctables.ProcessTable.tableName)[0].process_id

    # load up the files into individual documents
    xmldocs = [ligolw_add.ligolw_add(ligolw.Document(), [fname]) for fname in filenames]


    # Get the list of dinstinct segment_definers across all docs
    segment_definers = {}

    def register_definer(seg_def):
        key = (seg_def.ifos, seg_def.name, seg_def.version)
        segment_definers[key] = True
        return key

    for xmldoc in xmldocs:
        seg_def_table = table.get_table(xmldoc, lsctables.SegmentDefTable.tableName)
        map (register_definer, seg_def_table)

    # For each unique segment definer, find the intersection
    for ifo, name, version in segment_definers:
        if operation == INTERSECT:
            # If I were feeling especially functional-ist I'd write this
            # with reduce()
            result = glue.segments.segmentlist([glue.segments.segment(-glue.segments.infinity(), glue.segments.infinity())])
            for xmldoc in xmldocs:
                result &= find_segments(xmldoc, '%s:%s:%d' % (ifo, name, version), use_segment_table)
        elif operation == UNION:
            result = glue.segments.segmentlist([])

            for xmldoc in xmldocs:
                result |= find_segments(xmldoc, '%s:%s:%d' % (ifo, name, version), use_segment_table)
        elif operation == DIFF:
            result = find_segments(xmldocs[0], '%s:%s:%d' % (ifo, name, version), use_segment_table)
            
            for xmldoc in xmldocs[1:]:
                result -= find_segments(xmldoc, '%s:%s:%d' % (ifo, name, version), use_segment_table)
        else:
            raise NameError ("%s is not a known operation (intersect, union or diff)" % operation)


        # Add a segment definer for the result
        seg_def_id = add_to_segment_definer(outdoc, proc_id, ifo, name, version)

        # Add the segments
        if use_segment_table:
            add_to_segment(outdoc, proc_id, seg_def_id, result)
        else:
            add_to_segment_summary(outdoc, proc_id, seg_def_id, result)

    # If we're preserving, also load up everything into the output document.
    if preserve:
        # Add them to the output document
        map(lambda x: outdoc.appendChild(x.childNodes[0]), xmldocs)

        # Merge the ligolw elements and tables
        ligolw_add.merge_ligolws(outdoc)
        ligolw_add.merge_compatible_tables(outdoc)

    return outdoc, abs(result)



def run_segment_operation(outdoc, filenames, segments, use_segment_table, operation, result_name = 'RESULT', preserve = True):
    """
    Performs an operation (intersect or union) across a set of segments.
    That is, given a set of files each with segment definers DMT-FLAG1,
    DMT-FLAG2 etc and a list of segments DMT-FLAG1,DMT-FLAG1 this returns

    RESULT = (table 1's DMT-FLAG1 union table 2's DMT-FLAG1 union ...)
             operation
             (table 1's DMT-FLAG2 union table 2's DMT-FLAG2 union ...)
             operation
    etc
    """

    proc_id = table.get_table(outdoc, lsctables.ProcessTable.tableName)[0].process_id

    if preserve:
        indoc = ligolw_add.ligolw_add(outdoc, filenames)
    else:
        indoc = ligolw_add.ligolw_add(ligolw.Document(), filenames)

    # Start with a segment covering all of time, then
    # intersect with each of the fields of interest
    keys = segments.split(',')

    if operation == INTERSECT:
        sgmntlist = glue.segments.segmentlist([glue.segments.segment(-glue.segments.infinity(), glue.segments.infinity())])

        for key in keys:
            sgmntlist &= find_segments(indoc, key, use_segment_table)

    elif operation == UNION:
        sgmntlist = glue.segments.segmentlist([])

        for key in keys:
            sgmntlist |= find_segments(indoc, key, use_segment_table)
    elif operation == DIFF:
        sgmntlist = find_segments(indoc, keys[0], use_segment_table)

        for key in keys[1:]:
            sgmntlist -= find_segments(indoc, key, use_segment_table)
    else:
        raise NameError("%s is not a known operation (intersect, union or diff)" % operation)


    # Add a segment definer and segments
    seg_def_id = add_to_segment_definer(outdoc, proc_id, '', result_name, 1)

    if use_segment_table:
        add_to_segment(outdoc, proc_id, seg_def_id, sgmntlist)
    else:
        add_to_segment_summary(outdoc, proc_id, seg_def_id, sgmntlist)

    return outdoc, abs(sgmntlist)

