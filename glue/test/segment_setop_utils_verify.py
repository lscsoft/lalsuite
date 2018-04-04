#!/usr/bin/env python
#
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
#								   Preamble
#
# =============================================================================
#


"""
Tests ligolw_segment_intersect, ligolw_segment_union and ligolw_segment_diff
"""

import os
import sys
import pickle
import tempfile

from random import uniform

from glue.ligolw import utils
from glue.ligolw import ligolw
from glue.ligolw import lsctables

from glue.segments import segment
from glue.segments import segmentlist


from glue.segmentdb.segmentdb_utils import add_to_segment
from glue.segmentdb.segmentdb_utils import add_to_segment_summary
from glue.segmentdb.segmentdb_utils import add_to_segment_definer

from glue.ligolw.utils.process import append_process



#
# =============================================================================
#
#								 Library Code
#
# =============================================================================
#

def random_segments(start_time, end_time, num_segments, max_len):
	"""Creates a random segmentlist composed of num_segments segments,
	each of which is entirely contained in the given interval and 
	has length no more than max_len"""

	def random_segment():
		seg_len   = int(uniform(1, max_len))
		seg_start = int(uniform(start_time, end_time - seg_len))
		return segment(seg_start, seg_start + seg_len)

	return segmentlist([random_segment() for x in range(num_segments)])



def make_random_document(num_procs, num_seg_defs, num_segs, num_seg_sums, start_time, end_time, max_len):
	"""Create a ligolw document with random segments and segment_summary"""
	doc = ligolw.Document()
	doc.appendChild(ligolw.LIGO_LW())

	# Add some processes
	proc_ids = []
	seg_def_ids = {}
	segment_map = {}
	segment_sum_map = {}

	for count in range(num_procs):
		proc_id = append_process(doc, "Test program %d" % count).process_id
		proc_ids.append(proc_id)
		seg_def_ids[proc_id] = []

	# Add some segment definers
	for count in range(num_seg_defs):
		proc_id	= proc_ids[int(uniform(0,num_procs))]
		seg_def_id = add_to_segment_definer(doc, proc_id, "H1", "TEST_SEG_%d" % count, 1) 
		seg_def_ids[proc_id].append(seg_def_id)

		# Add some segments
		sgmntlst = random_segments(start_time, end_time, num_segs, max_len)
		add_to_segment(doc, proc_id, seg_def_id, sgmntlst)
		sgmntlst.coalesce()
		segment_map["TEST_SEG_%d" % count] = sgmntlst

		# Add some segment summaries
		sgmntlst = random_segments(start_time, end_time, num_segs, max_len)
		add_to_segment_summary(doc, proc_id, seg_def_id, sgmntlst)
		sgmntlst.coalesce()
		segment_sum_map["TEST_SEG_%d" % count] = sgmntlst

	return doc, segment_map, segment_sum_map
	

def segments_from_cmd(cmd):
	"""Run the specified command, the output of which should be lines of the form
	start_time, end_time
	The output is parsed and used to construct a (coalesced) segment list"""
	resultsegs = []

	pipe = os.popen(cmd)

	for l in pipe:
		start_t, end_t = l.strip().split(',')
		resultsegs.append(segment(int(start_t), int(end_t)))

	resultlst = segmentlist(resultsegs)
	resultlst.coalesce()

	return resultlst


def do_test(cmd, should_be, error_msg):
	"""Do a single test involving comparing the output of cmd to what it should be"""

	print 'Testing ' + error_msg

	resultlst = segments_from_cmd(cmd)
	
	if resultlst != should_be:
		print >>sys.stderr, 'Failure in ' + error_msg
		print resultlst
		print should_be

		sys.exit(-1)


#
# =============================================================================
#
#									 Main
#
# =============================================================================
#

if __name__ == '__main__':
	del lsctables.SegmentTable.validcolumns['start_time_ns']
	del lsctables.SegmentTable.validcolumns['end_time_ns']
	del lsctables.ProcessTable.validcolumns['domain']
	del lsctables.ProcessTable.validcolumns['jobid']
	del lsctables.ProcessTable.validcolumns['is_online']


	segment_map_map = {}
	segment_sum_map_map = {}

	# Create two random files, each with 2 procs, 4 segment definers each with 8 segs and 8 seg summaries
	tmp_dir = tempfile.mkdtemp()

	for filename in ['file1.xml', 'file2.xml']:
		doc, segment_map, segment_sum_map = make_random_document(2, 4, 8, 8, 800000000, 800000200, 10)
		segment_map_map[filename] = segment_map
		segment_map_map[filename] = segment_map
		segment_sum_map_map[filename] = segment_sum_map
		utils.write_filename(doc, tmp_dir + "/" + filename)

	# find the intersection between two segment definers
	do_test('cat %s/file1.xml | ligolw_segment_intersect -i H1:TEST_SEG_1,H1:TEST_SEG_2 --segment | ligolw_print -t segment -c start_time -c end_time' % tmp_dir,
			segment_map_map['file1.xml']['TEST_SEG_1'] & segment_map_map['file1.xml']['TEST_SEG_2'],
			'intersect between segments')


	# union
	do_test('cat %s/file1.xml | ligolw_segment_union -i H1:TEST_SEG_1,H1:TEST_SEG_2 --segment | ligolw_print -t segment -c start_time -c end_time' % tmp_dir,
			segment_map_map['file1.xml']['TEST_SEG_1'] | segment_map_map['file1.xml']['TEST_SEG_2'],
			'union between segments')


	# difference
	do_test('cat %s/file1.xml | ligolw_segment_diff -i H1:TEST_SEG_1,H1:TEST_SEG_2 --segment | ligolw_print -t segment -c start_time -c end_time' % tmp_dir,
			segment_map_map['file1.xml']['TEST_SEG_1'] - segment_map_map['file1.xml']['TEST_SEG_2'],
			'difference between segments')




	# find the intersection between two segment summaries
	do_test('cat %s/file1.xml | ligolw_segment_intersect -i H1:TEST_SEG_1,H1:TEST_SEG_2 --segment-sum | ligolw_print -t segment_summary -c start_time -c end_time' % tmp_dir,
			segment_sum_map_map['file1.xml']['TEST_SEG_1'] & segment_sum_map_map['file1.xml']['TEST_SEG_2'],
			'intersect between segment summaries')


	# union
	do_test('cat %s/file1.xml | ligolw_segment_union -i H1:TEST_SEG_1,H1:TEST_SEG_2 --segment-sum | ligolw_print -t segment_summary -c start_time -c end_time' % tmp_dir,
			segment_sum_map_map['file1.xml']['TEST_SEG_1'] | segment_sum_map_map['file1.xml']['TEST_SEG_2'],
			'union between segment summaries')


	# difference
	do_test('cat %s/file1.xml | ligolw_segment_diff -i H1:TEST_SEG_1,H1:TEST_SEG_2 --segment-sum | ligolw_print -t segment_summary -c start_time -c end_time' % tmp_dir,
			segment_sum_map_map['file1.xml']['TEST_SEG_1'] - segment_sum_map_map['file1.xml']['TEST_SEG_2'],
			'difference between segment summaries')



	# find the intersection between two segment definers across files
	do_test('ligolw_segment_intersect --segment %s/file1.xml %s/file2.xml | ' % (tmp_dir, tmp_dir) + 
			'ligolw_segment_intersect --segment -i H1:TEST_SEG_1 | ' +
			'ligolw_print -t segment -c start_time -c end_time', 
			segment_map_map['file1.xml']['TEST_SEG_1'] & segment_map_map['file2.xml']['TEST_SEG_1'],
			'intersect between files')

	# union
	do_test('ligolw_segment_union --segment %s/file1.xml %s/file2.xml | ' % (tmp_dir, tmp_dir) + 
			'ligolw_segment_intersect --segment -i H1:TEST_SEG_1 | ' +
			'ligolw_print -t segment -c start_time -c end_time', 
			segment_map_map['file1.xml']['TEST_SEG_1'] | segment_map_map['file2.xml']['TEST_SEG_1'],
			'union between files')

	# diff
	do_test('ligolw_segment_diff --segment %s/file1.xml %s/file2.xml | ' % (tmp_dir, tmp_dir) + 
			'ligolw_segment_intersect --segment -i H1:TEST_SEG_1 | ' +
			'ligolw_print -t segment -c start_time -c end_time', 
			segment_map_map['file1.xml']['TEST_SEG_1'] - segment_map_map['file2.xml']['TEST_SEG_1'],
			'difference between files')

	for filename in ['file1.xml', 'file2.xml']:
		os.remove(tmp_dir + "/" + filename)
	os.rmdir(tmp_dir)


