#!/usr/bin/env python
#
# Copyright (C) 2015 Chris Pankow
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

"""
Given a set of extrinsic evidence calculations on a given set of intrinsic parameters, refines the grid to do additional calculations.
"""

__author__ = "Chris Pankow <chris.pankow@ligo.org>"

import os
import sys
import glob
import json
import bisect
from collections import defaultdict
from optparse import OptionParser, OptionGroup

import numpy
from scipy.special import binom

from glue.ligolw import utils, ligolw, lsctables, ilwd
lsctables.use_in(ligolw.LIGOLWContentHandler)

from glue.ligolw.utils import process

def midpoint(pt1, pt2):
    diff = pt2 - pt1
    mid = diff / 2
    return pt1 + mid

def prune_duplicate_pts(pts):
    return numpy.array(list(set([tuple(pt) for pt in pts])))

#
# Option parsing
#

optp = OptionParser()
# Options needed by this program only.

optp.add_option("-o", "--output-glob", help="Where to find the output of the inttrinisic integrator.")
optp.add_option("-i", "--intrinsic-param", action="append", help="Adapt in this intrinsic parameter.")
optp.add_option("-O", "--use-overlap", help="Use overlap information to define 'closeness'.")
#optp.add_option("-m", "--refine-method", help="Midpoint or cell division")

opts, args = optp.parse_args()

class OverlapLookup(object):
    def __init__(self, npts):
        # This could get ugly.
        self._odata = numpy.zeros((npts, npts))
        self._intr_info = [None] * npts

    def add_overlap_row(self, idx, intr_prms, overlap_data, strict=False):
        # Overlap data could be incomplete, so we account for that here
        for row in overlap_data:
            h2_idx = row[-1]
            h2_overlap = row[-2]
            # FIXME:
            #intr_prms = row[:-2]
            intr_prms = row[:-2][:2]

            if self._intr_info[int(h2_idx)] is None:
                self._intr_info[int(h2_idx)] = intr_prms
            elif strict:
                assert all(numpy.isclose(self._intr_info[int(h2_idx)], intr_prms))
            self._intr_info[int(h2_idx)] = intr_prms
            self._odata[idx,int(h2_idx)] = h2_overlap

    def query(self, idx, npts, thresh=None):
        """
        Query a row (indexed with idx) of the overlap matrix and get the indexes of the highest npts.
        """
        # Sort by overlap value for this row
        sort_idx = self._odata[idx].argsort()

        # If we instead want a threshold, we have to find it
        if thresh is not None:
            # Find the value in the sorted overlaps which divides the indexes
            # to overlap values less or greater than the value
            tidx = bisect.bisect(self._odata[idx][sort_idx], thresh)
            return self._odata[idx][sort_idx][-tidx:]
        #return self._odata[idx][sort_idx][-npts:]
        return sort_idx[-npts:]

    def intr_prms_from_idx(self, idx):
        return self._intr_info[idx]

tree = OverlapLookup(5762)

# If asked, retrieve bank overlap
if opts.use_overlap is not None:
    with open(opts.use_overlap) as fin:
        overlap_toc = json.load(fin)["types"]
    # Just get the first type...
    overlap_toc = overlap_toc[overlap_toc.keys()[0]]

    # FIXME:
    bpath = os.path.dirname(opts.use_overlap)
    overlap = {}
    for idx, info in enumerate(overlap_toc):
        print os.path.join(bpath, info["filename"])
        with open(os.path.join(bpath, info["filename"])) as ovrlpfile:
            contents = json.load(ovrlpfile)
            overlap[idx] = numpy.array(contents["overlap"])
            tree.add_overlap_row(idx, (info["mass1"], info["mass2"]), numpy.array(contents["overlap"]))

tree.query(0, 1)
# TODO: Map an input point onto the template bank intrinsic space
# Hopefully the point is already present and we can just get it, otherwise it
# could incur an overlap calculation, or suffer from the effects of being close
# only in Euclidean terms

#
# Step 1: retrieve previous marginalized likelihood mapping
#

xmlfiles = glob.glob(opts.output_glob)
if not xmlfiles:
    exit("No files found with output glob.")

result_table = None
for xfile in xmlfiles:
    xmldoc = utils.load_filename(xfile, contenthandler=ligolw.LIGOLWContentHandler)
    if result_table is None:
        result_table = lsctables.SnglInspiralTable.get_table(xmldoc)
    else:
        result_table.extend(lsctables.SnglInspiralTable.get_table(xmldoc))

if not result_table:
    exit("Empty result table, nothing to process.")

#
# Step 2: Determine points which need refinement
#
from sklearn.neighbors import BallTree
# intr_prms = tuple(opts.intrinsic_params)
intr_prms = ("mass1", "mass2")
pts = numpy.array([tuple(getattr(t, a) for a in intr_prms) for t in result_table])
tree = BallTree(pts)

#
# Step 3: Do refinement
#

# The number of points to retrieve for an arbitrary scheme is not well defined.
# The number was chosen since that's the number of nearest neighbors including
# oneself and bisector points in all dimensions for a regular nD grid
ndim = len(intr_prms)
npts = int(numpy.sum(2**(ndim-m) * binom(ndim, m) for m in range(ndim)) + 1)
npts = min(npts, len(pts))

verbose = True

log_evid_thresh = 0
new_mid_pts = []
for pt in result_table:
    if pt.snr <= log_evid_thresh:
        continue

    pt_info = numpy.array([getattr(pt, a) for a in intr_prms])
    if verbose:
        print "Examining point with ID %d" % int(pt.event_id)
        print "and values (%s)" % ", ".join(map(str, pt_info))

    if verbose:
        sys.stderr.write("Querying for %d neighbors... " % (npts-1))
    dist, idx = tree.query(pt_info, k=npts, return_distance=True)
        
    # Remove the selected point
    # FIXME: Also remove things beyond a maximum distance if necessary
    dist = dist[0,1:]
    idx = idx[0,1:]
    if verbose:
        print "retrieved %d" % len(idx)
        print "\t" + "\n\t".join( ("%d %f" % a for a in zip(idx, dist)) )

    # get bisectors
    new_mid_pts.append(numpy.array(midpoint(pt_info, pts[idx])).T)

# Remove redundant points
# NOTE: This is only likely to happen in the case where two dimensions have a
# regular (and I think identically spaced) grid -- this could be modified to be
# more like "remove points which are too close".
new_pts = prune_duplicate_pts(numpy.hstack(new_mid_pts))

#
# Step 4: Set up next set of runs
#

##################################
xmldoc = ligolw.Document()
xmldoc.appendChild(ligolw.LIGO_LW())
procrow = process.append_process(xmldoc, program=sys.argv[0])
procid = procrow.process_id
process.append_process_params(xmldoc, procrow, process.process_params_from_dict(opts.__dict__))

rows = ["simulation_id", "process_id", "numrel_data"] + list(intr_prms)
sim_insp_tbl = lsctables.New(lsctables.SimInspiralTable, rows)
for itr, intr_prm in enumerate(new_pts):
    sim_insp = sim_insp_tbl.RowType()
    # FIXME: Need better IDs
    sim_insp.numrel_data = "MASS_SET_%d" % itr
    sim_insp.simulation_id = ilwd.ilwdchar("sim_inspiral:sim_inspiral_id:%d" % itr)
    sim_insp.process_id = procid
    for p, v in zip(intr_prms, intr_prm):
        setattr(sim_insp, p, v)
    sim_insp_tbl.append(sim_insp)

xmldoc.childNodes[0].appendChild(sim_insp_tbl)
channel_name = ["H=H", "L=L"]
ifos = "".join([o.split("=")[0][0] for o in channel_name])
#start = int(event_time)
start = 0
fname = "%s-MASS_POINTS-%d-1.xml.gz" % (ifos, start)
utils.write_filename(xmldoc, fname, gz=True, verbose=verbose)
