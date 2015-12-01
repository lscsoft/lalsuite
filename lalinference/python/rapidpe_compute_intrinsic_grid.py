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
import re
from collections import defaultdict
from optparse import OptionParser, OptionGroup
from copy import copy

import h5py
import numpy
from scipy.special import binom
from sklearn.neighbors import BallTree

from glue.ligolw import utils, ligolw, lsctables, ilwd
lsctables.use_in(ligolw.LIGOLWContentHandler)
from glue.ligolw.utils import process

import lalsimulation
from lalinference.rapid_pe import amrlib, lalsimutils, common_cl

def get_cr_from_grid(cells, weight, cr_thr=0.9):
    """
    Given a set of cells and the weight of that cell, calculate a N% CR including cells which contribute to that probability mass.
    """
    if cr_thr == 0.0:
        return numpy.empty((0,))

    # Arrange them all with their respective weight
    cell_sort = numpy.hstack( (weight[:,numpy.newaxis], cells) )

    # Sort and form the CDF
    cell_sort = cell_sort[cell_sort[:,0].argsort()]
    cell_sort[:,0] = cell_sort[:,0].cumsum()
    cell_sort[:,0] /= cell_sort[-1,0]

    # find the CR probability
    idx = cell_sort[:,0].searchsorted(1-cr_thr)

    return cell_sort[idx:,1:]

# FIXME: Should this be moved to amrlib --- seems a little too specialized?
ovrlp_obj = lalsimutils.Overlap(fLow=40, fMax=2000, deltaF=0.125, psd=lalsimulation.SimNoisePSDeLIGOModel)
def eval_grid(t1, cells, intr_prms):

    t1_norm = None
    t1_copy = copy(t1)

    overlaps = []
    for cell in cells:

        try:
            prms = intr_prms
            tmp = numpy.array(cell._center)[numpy.newaxis,:]
            tmp = amrlib.apply_inv_transform(tmp, prms, opts.distance_coordinates).T
            for p, pval in zip(prms, tmp):
                #if p == "lam_tilde":
                    #t1_copy.psi0 = pval[0]
                    #t1_copy.psi3 = 0
                #else:
                setattr(t1_copy, p, pval[0])
        except AssertionError:
            rejected.append((cell, 0.0))
            continue

        # FIXME: Hardcoded waveforms
        olap, t1_norm, _ = lalsimutils.overlap(t1, t1_copy, ovrlp_obj, 0.125, 40, "TaylorF2", "TaylorF2", t1_norm=t1_norm)
        overlaps.append(olap)

    return overlaps

def determine_region(pt, pts, ovrlp, ovrlp_thresh, expand_prms={}):
    """
    Given a point (pt) in a set of points (pts), with a function value at those points (ovrlp), return a rectangular hull such that the function exceeds the value ovrlp_thresh.
    """
    sidx = bisect.bisect(ovrlp, ovrlp_thresh)
    #print "Found %d neighbors with overlap >= %f" % (len(ovrlp[sidx:]), ovrlp_thresh)

    cell = amrlib.Cell.make_cell_from_boundaries(pt, pts[sidx:])
    for k, lim in expand_prms.iteritems():
        cell._bounds = numpy.vstack((cell._bounds, lim))
        # FIXME: Need to do center?
    return cell, sidx

def find_olap_index(tree, exact=True, **kwargs):
    """
    Given an object that can retrieve distance via a 'query' function (e.g. KDTree or BallTree), find the index of a point closest to the input point. Note that kwargs is used to get the current known values of the event. E.g.

    intr_prms = {'mass1': 1.4, 'mass2': 1.35}
    find_olap_index(tree, **intr_prms)
    """
    pt = numpy.array([kwargs[k] for k in intr_prms])

    # FIXME: Replace with standard function
    dist, m_idx = tree.query(pt, k=1)
    dist, m_idx = dist[0][0], int(m_idx[0][0])

    # FIXME: There's still some tolerance from floating point conversions
    if exact and dist > 0.000001:
        exit("Could not find template in bank, closest pt was %f away" % dist)
    return m_idx, pt

def write_to_xml(cells, intr_prms, fvals=None, fname=None, verbose=False):
    """
    Write a set of cells, with dimensions corresponding to intr_prms to an XML file as sim_inspiral rows.
    """
    xmldoc = ligolw.Document()
    xmldoc.appendChild(ligolw.LIGO_LW())
    procrow = process.append_process(xmldoc, program=sys.argv[0])
    procid = procrow.process_id
    process.append_process_params(xmldoc, procrow, process.process_params_from_dict(opts.__dict__))

    rows = ["simulation_id", "process_id", "numrel_data"] + list(intr_prms)
    if fvals is not None:
        rows.append("alpha1")
    sim_insp_tbl = lsctables.New(lsctables.SimInspiralTable, rows)
    for itr, intr_prm in enumerate(cells):
        sim_insp = sim_insp_tbl.RowType()
        # FIXME: Need better IDs
        sim_insp.numrel_data = "INTR_SET_%d" % itr
        sim_insp.simulation_id = ilwd.ilwdchar("sim_inspiral:sim_inspiral_id:%d" % itr)
        sim_insp.process_id = procid
        if fvals:
            sim_insp.alpha1 = fvals[itr]
        for p, v in zip(intr_prms, intr_prm._center):
            setattr(sim_insp, p, v)
        sim_insp_tbl.append(sim_insp)

    xmldoc.childNodes[0].appendChild(sim_insp_tbl)
    if fname is None:
        channel_name = ["H=H", "L=L"]
        ifos = "".join([o.split("=")[0][0] for o in channel_name])
        #start = int(event_time)
        start = 0
        fname = "%s-MASS_POINTS-%d-1.xml.gz" % (ifos, start)
    utils.write_filename(xmldoc, fname, gz=True, verbose=verbose)

#
# Plotting utilities
#
def plot_grid_cells(cells, color, axis1=0, axis2=1):
    from matplotlib.patches import Rectangle
    from matplotlib import pyplot
    ax = pyplot.gca()
    for cell in cells:
        ext1 = cell._bounds[axis1][1] - cell._bounds[axis1][0]
        ext2 = cell._bounds[axis2][1] - cell._bounds[axis2][0]

        ax.add_patch(Rectangle((cell._bounds[axis1][0], cell._bounds[axis2][0]), ext1, ext2, edgecolor = color, facecolor='none'))

#
# Option parsing
#

def parse_param(popts):
    """
    Parse out the specification of the intrinsic space. Examples:

    >>> parse_param(["mass1=1.4", "mass2", "spin1z=-1.0,10"])
    {'mass1': 1.4, 'mass2': None, 'spin1z': (-1.0, 10.0)}
    """
    intr_prms, expand_prms = {}, {}
    for popt in popts:
        popt = popt.split("=")
        if len(popt) == 1:
            # Implicit expand in full parameter space -- not yet completely
            # implemented
            intr_prms[popt[0]] = None
        elif len(popt) == 2:
            popt[1] = popt[1].split(",")
            if len(popt[1]) == 1:
                # Fix intrinsic point
                intr_prms[popt[0]] = float(popt[1][0])
            else:
                expand_prms[popt[0]] = tuple(map(float, popt[1]))
    return intr_prms, expand_prms

optp = OptionParser()

optp.add_option("-d", "--distance-coordinates", default="tau0_tau3", help="Coordinate system in which to calculate 'closeness'. Default is tau0_tau3.")
optp.add_option("-n", "--no-exact-match", action="store_true", help="Loosen criteria that the input intrinsic point must be a member of the input template bank.")
optp.add_option("-v", "--verbose", action='store_true', help="Be verbose.")

# FIXME: These two probably should only be for the initial set up. While it
# could work, in theory, for refinement, the procedure would be a bit more
# tricky.
# FIXME: This could be a single value (lock a point in) or a range (adapt across
# this is range). No argument given implies use entire known range (if
# available).
optp.add_option("-i", "--intrinsic-param", action="append", help="Adapt in this intrinsic parameter. If a pre-existing value is known (e.g. a search template was identified), specify this parameter as -i mass1=1.4 . This will indicate to the program to choose grid points which are commensurate with this value.")
# FIXME: This option to be replaced with semantics from above
#optp.add_option("-e", "--expand-param", action="append", help="Expand in this intrinsic parameter.")

grid_section = OptionGroup(optp, "initial gridding options", "Options for setting up the initial grid.")
grid_section.add_option("--setup", help="Set up the initial grid based on template bank overlaps. The new grid will be saved to this argument, e.g. --setup grid will produce a grid.npy file.")
grid_section.add_option("-t", "--tmplt-bank", help="XML file with template bank.")
grid_section.add_option("-O", "--use-overlap", help="Use overlap information to define 'closeness'.")
grid_section.add_option("-T", "--overlap-threshold", type=float, help="Threshold on overlap value.")
grid_section.add_option("-N", "--no-deactivate", action="store_true", help="Do not deactivate cells initially which have no template within them.")
optp.add_option_group(grid_section)

refine_section = OptionGroup(optp, "refine options", "Options for refining a pre-existing grid.")
refine_section.add_option("--refine", help="Refine a prexisting grid. Pass this option the grid points from previous levels (or the --setup) option.")
refine_section.add_option("-r", "--result-file", help="XML file containing newest result to refine.")
optp.add_option_group(refine_section)

opts, args = optp.parse_args()

if not (opts.setup or opts.refine):
    exit("Either --setup or --refine must be chosen")

# If asked, retrieve bank overlap
if opts.use_overlap is not None:
    h5file = h5py.File(opts.use_overlap, "r")

    # FIXME:
    #wfrm_fam = args.waveform_type
    # Just get the first one
    wfrm_fam = h5file.keys()[0]

    odata = h5file[wfrm_fam]
    m1, m2, ovrlp = odata["mass1"], odata["mass2"], odata["overlaps"]
    if opts.verbose:
        print "Using overlap data from %s" % wfrm_fam

# Hopefully the point is already present and we can just get it, otherwise it
# could incur an overlap calculation, or suffer from the effects of being close
# only in Euclidean terms

intr_prms, expand_prms = parse_param(opts.intrinsic_param)
intr_pt = numpy.array([intr_prms[k] for k in intr_prms])
intr_prms = intr_prms.keys()

#
# Step 1: retrieve templates / result
#
xmldoc = utils.load_filename(opts.tmplt_bank, contenthandler=ligolw.LIGOLWContentHandler)
tmplt_bank = lsctables.SnglInspiralTable.get_table(xmldoc)

results = []
if opts.result_file:
    for arg in glob.glob(opts.result_file):
        # FIXME: Bad hardcode
        if "samples" in arg:
            continue
        xmldoc = utils.load_filename(arg, contenthandler=ligolw.LIGOLWContentHandler)
        results.extend(lsctables.SnglInspiralTable.get_table(xmldoc))

    # Normalize
    # FIXME: If we have more than 1 copies -- This is tricky because we need to
    # pare down the duplicate sngl rows too
    maxlnevid = numpy.max([s.snr for s in results])
    total_evid = numpy.exp([s.snr - maxlnevid for s in results]).sum()
    for res in results:
        res.snr = numpy.exp(res.snr - maxlnevid)/total_evid

    res_pts = numpy.array([tuple(getattr(t, a) for a in intr_prms) for t in results])
    res_pts = amrlib.apply_transform(res_pts, intr_prms, opts.distance_coordinates)
    # FIXME: this needs to be done in a more consistent way
    results = numpy.array([res.snr for res in results])

#
# Step 2: Set up metric space
#

# NOTE: We use the template bank here because the overlap results might not have
# all the intrinsic information stored (e.g.: no spins, even though the bank is
# aligned-spin).
# FIXME: This is an oversight in the overlap calculator that should be
# rectified.
pts = numpy.array([tuple(getattr(t, a) for a in intr_prms) for t in tmplt_bank])
pts = amrlib.apply_transform(pts, intr_prms, opts.distance_coordinates)

# Transform and repack initial point
intr_pt = amrlib.apply_transform(intr_pt[numpy.newaxis,:], intr_prms, opts.distance_coordinates)[0]
intr_pt = dict(zip(intr_prms, intr_pt))

# FIXME: Can probably be moved to point index identification function -- it's
# not used again
tree = BallTree(pts)

#
# Step 3: Get the row of the overlap matrix to work with
#
m_idx, pt = find_olap_index(tree, not opts.no_exact_match, **intr_pt)

# Save the template for later use as well
t1 = tmplt_bank[m_idx]

#
# Rearrange data to correspond to input point
#

sort_order = ovrlp[m_idx].argsort()
ovrlp = numpy.array(ovrlp[m_idx])[sort_order]

# DANGEROUS: This assumes the (template bank) points are the same order as the
# overlaps. While we've taken every precaution to ensure this is true, it may
# not always be.
pts = pts[sort_order]
m_idx = sort_order[m_idx]

"""
if results is None:
    # Initial gridding
    init_region, idx = determine_region(pt, pts, ovrlp, opts.overlap_threshold)
else:
    # FIXME: This needs to be recalculated for the spacing, we should just save
    # and retrieve it in the grid file as necessary
    init_region, idx = determine_region(pt, res_pts, results, 0)
"""

intr_prms = list(intr_prms) + expand_prms.keys()
init_region, idx = determine_region(pt, pts, ovrlp, opts.overlap_threshold, expand_prms)

####### BEGIN INITIAL GRID CODE #########

# FIXME: To be reimplemented in a different way
#if opts.expand_param is not None:
    #expand_param(init_region, opts.expand_param)

extent_str = " ".join("(%f, %f)" % bnd for bnd in map(tuple, init_region._bounds))
center_str = " ".join(map(str, init_region._center))
print "Created cell with center " + center_str + " and extent " + extent_str

# TODO: Alternatively, check density of points in the region to determine the
# points to a side
grid, spacing = amrlib.create_regular_grid_from_cell(init_region, side_pts=5, return_cells=True)

# "Deactivate" cells not close to template points
# FIXME: This gets more and more dangerous in higher dimensions
# FIXME: Move to function
tree = BallTree(grid)
if not opts.no_deactivate:
    get_idx = set()
    for pt in pts[idx:]:
        get_idx.add(tree.query(pt, k=1, return_distance=False)[0][0])
    selected = grid[numpy.array(list(get_idx))]
else:
    selected = grid

#### BEGIN REFINEMENT OF RESULTS #########

if opts.result_file is not None:
    # FIXME: temporary -- we shouldn't need to know the spacing of the initial
    # region
    prev_cells, spacing = amrlib.deserialize_grid_cells(opts.refine)
    selected = numpy.array([c._center for c in prev_cells])
    # refinement already occurred once in the setup stage
    #spacing /= 2
    selected = amrlib.apply_transform(selected, intr_prms, opts.distance_coordinates)

    grid_tree = BallTree(selected)
    grid_idx = []
    # Reorder the grid points to match their weight indices
    for res in res_pts:
        dist, idx = grid_tree.query(res, k=1)
        # Stupid floating point inexactitude...
        #print res, selected[idx[0][0]]
        #assert numpy.allclose(res, selected[idx[0][0]])
        grid_idx.append(idx[0][0])

    if opts.refine:
        selected = get_cr_from_grid(selected[grid_idx], results, cr_thr=opts.overlap_threshold)
        print "Selected %d cells from %3.2f%% confidence region" % (len(selected), opts.overlap_threshold*100)

grid, spacing = amrlib.refine_regular_grid(selected, spacing, return_cntr=opts.setup)
print "%d cells after refinement" % len(grid)
grid = amrlib.prune_duplicate_pts(grid, init_region._bounds, spacing)

#
# Clean up
#

grid = numpy.array(grid)
bounds_mask = amrlib.check_grid(grid, intr_prms, opts.distance_coordinates)
grid = grid[bounds_mask]
print "%d cells after bounds checking" % len(grid)

if len(grid) == 0:
    exit("All cells would be removed by physical boundaries.")

# Convert back to physical mass
grid = amrlib.apply_inv_transform(grid, intr_prms, opts.distance_coordinates)

cells = amrlib.grid_to_cells(grid, spacing)
if opts.setup:
    npy_grid = amrlib.pack_grid_cells(cells, spacing)
    # FIXME: specify level better
    amrlib.serialize_grid_cells({1: npy_grid}, opts.setup)
else:
    npy_grid = amrlib.pack_grid_cells(cells, spacing)
    # FIXME: -1 will mean "whatever the next level is"
    #amrlib.serialize_grid_cells({-1: npy_grid}, opts.refine)
    amrlib.serialize_grid_cells({1: npy_grid}, opts.refine)

overlaps = eval_grid(t1, cells, intr_prms)

print "Selected %d cells for further analysis." % len(cells)
if opts.setup:
    fname = "HL-MASS_POINTS_LEVEL_0-0-1.xml.gz"
    write_to_xml(cells, intr_prms, overlaps, fname, verbose=opts.verbose)
else:
    m = re.search("LEVEL_(\d+)", opts.result_file)
    if m is not None:
        level = int(m.group(1)) + 1
        fname = "HL-MASS_POINTS_LEVEL_%d-0-1.xml.gz" % level
    else:
        fname = "HL-MASS_POINTS_LEVEL_X-0-1.xml.gz"
    write_to_xml(cells, intr_prms, overlaps, fname, verbose=opts.verbose)
