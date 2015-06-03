#!/usr/bin/env python
import os
import json
import argparse

import numpy
from sklearn.neighbors import BallTree

import lal
import lalsimulation
import lalinspiral
from lalinference.rapid_pe import lalsimutils, amrlib

from glue.ligolw import ligolw, utils, lsctables
from glue.ligolw.utils import process

lsctables.use_in(ligolw.LIGOLWContentHandler)
from pylal.series import read_psd_xmldoc, LIGOLWContentHandler

VALID_TMPLT_GENS = {"lalapps_cbc_sbank": "--flow", "tmpltbank": "--low-frequency-cutoff", "pycbc_geom_aligned_bank": "--f-low"}
def infer_flow(xmldoc):
    """
    Attempt to infer the low frequency by combing through the process table and trying to pick out the low frequency option given to that program. If you trust this, you will, for sure, be disappointed at some point in using this program.
    """

    proctable = lsctables.ProcessTable.get_table(xmldoc)
    # FIXME: ...but really, I don't you think can fix this...
    procs = set([p.program for p in proctable if VALID_TMPLT_GENS.has_key(p.program)])

    if len(procs) == 0:
        return None

    # FIXME: You askin' for trouble, son.
    try:
        return min([min(process.get_process_params(xmldoc, prog, VALID_TMPLT_GENS[prog], False)) for prog in procs])
    except ValueError:
        pass # No flow found. Bad luck for you.

    return None 

def parse_psd_file(filestr, fvals):
    """
    Map the user-provided PSD file string into a function to be called as PSD(f).
    """ 
  
    if not os.path.isfile(filestr):
        try:
            psd_func = getattr(lalsimulation, filestr)
            return numpy.array(map(psd_func, fvals))
        except AttributeError:
            pass

    try:
        xmldoc = utils.load_filename(filestr, contenthandler=LIGOLWContentHandler)
        psd = read_psd_xmldoc(xmldoc).values()[0]
        f = numpy.arange(0, len(psd.data)*psd.deltaF, psd.deltaF)
        psd = psd.data
    except:
        # FIXME: ugh!
        try:
            f, psd = numpy.loadtxt(filestr, unpack=True)
        except:
           exit("Can't parse PSD specifier %s as function or file." % filestr)

    def anon_interp(newf):
        return numpy.interp(newf, f, psd)
    return numpy.array(map(anon_interp, fvals))

argp = argparse.ArgumentParser()
argp.add_argument("-s", "--tmplt-start-index", type=int, help="Start at this index of the template bank.")
argp.add_argument("-e", "--tmplt-end-index", type=int, help="End at this index of the template bank.")
argp.add_argument("-t", "--tmplt-bank-file", help="File name of the template bank. Required.")
argp.add_argument("-d", "--distance-coordinates", default="mchirp_eta", help="Coordinate system in which to calculate 'closness'. Default is mchirp_eta.")
argp.add_argument("-p", "--psd-file", help="Name of PSD XML file. Required.")
argp.add_argument("-f", "--f-low", type=float, help="Lowest frequency component of template. Will attempt to infer from template bank, else must be provided.")
argp.add_argument("-F", "--delta-f", type=float, default=0.125, help="Frequency binning of the FD waveform. Default is 0.125.")
argp.add_argument("-a", "--approximant1", default="TaylorF2", help="Approximant to use for target waveform. Default is TaylorF2.")
argp.add_argument("-b", "--approximant2", default="TaylorF2", help="Approximant to use for overlapped waveform. Default is TaylorF2.")
argp.add_argument("-v", "--verbose", action="store_true", help="Be verbose.")
argp.add_argument("-V", "--too-verbose", action="store_true", help="Be absolutely, obnoxiously loquacious.")
args = argp.parse_args()

## DEFAULTS ##
f_high = 2048.0

if not args.tmplt_bank_file or not os.path.exists(args.tmplt_bank_file):
    exit("Template bank file either not specified or has an invalid path")

#
# Generate discrete PSD
#
delta_f = args.delta_f
fvals = numpy.arange(0, f_high, delta_f)
psd = parse_psd_file(args.psd_file, fvals)
    
#
# Extract and prepare template bank
#
xmldoc = utils.load_filename(args.tmplt_bank_file, contenthandler=ligolw.LIGOLWContentHandler)
tmplt_bank = lsctables.SnglInspiralTable.get_table(xmldoc)

if args.f_low is None:
    f_low = infer_flow(xmldoc)
    if args.verbose:
        print "Low frequency inferred from template bank is %f" % f_low
else:
    f_low = args.f_low
    if args.verbose:
        print "Low frequency from command line is %f" % f_low

if f_low is None:
    exit("Low frequency cutoff could not be inferred from template bank, and none was given.")

# lalapps_tmpltbank assigns 0 ID to all events, so we remap
# FIXME: Check for tmplt_bank: All others do assign IDs
for tmplt in tmplt_bank:
    tmplt.event_id = tmplt_bank.get_next_id()

# FIXME: Unhardcode
wtype = "%s_%s" % (args.approximant1, args.approximant2)
toc = {"types": {wtype: []}}
bdir = "%s/" % wtype
if not os.path.exists(bdir):
    os.mkdir(bdir)

# FIXME: This code should probably be libized
# FIXME: Unhardcode
intr_prms = ("mass1", "mass2", "spin1z", "spin2z")
pts = numpy.array([tuple(getattr(t, a) for a in intr_prms) for t in tmplt_bank])

# Dump full m1/m2 bank to JSON
with open("bank.json", "w") as fout:
    json.dump([list(a) + ["%s/%s_%d.json" % (bdir, wtype, i)] for i, a in enumerate(pts)], fout)

pts = amrlib.apply_transform(pts, intr_prms, args.distance_coordinates)

#
# Construct objects needed to identify neighbors and do overlaps
#
tree = BallTree(pts)

ovrlp = lalsimutils.Overlap(fLow=f_low, fMax=2000, deltaF=delta_f, psd=psd, analyticPSD_Q=False)

idx_range = range(args.tmplt_start_index or 0, args.tmplt_end_index or len(tmplt_bank))

# FIXME:
npts = len(tmplt_bank)
npts = 1000
for i1, pt in enumerate(pts):
    opt = amrlib.apply_inv_transform(pts[i1,numpy.newaxis].copy(), intr_prms, "mchirp_eta")[0]
    fname = "%s/%s_%d.json" % (bdir, wtype, i1)
    # FIXME: This makes my eyes bleed...
    toc["types"][wtype].append({})
    toc["types"][wtype][-1]["mass1"] = opt[0]
    toc["types"][wtype][-1]["mass2"] = opt[1]
    toc["types"][wtype][-1]["filename"] = fname

    if i1 not in idx_range:
        continue

    if os.path.exists(fname):
        continue

    dist, idx = tree.query(pt, k=npts, return_distance=True)

    t1 = tmplt_bank[i1]
    h1 = lalsimutils.generate_waveform_from_tmplt(t1, args.approximant1, delta_f, f_low)
    h1_norm = ovrlp.norm(h1)
    if args.verbose:
        print "--- (%f, %f) / (%f, %f)" % (t1.mass1, t1.mass2, t1.mchirp, t1.eta)

    ovrlps = []
    for d, i2 in numpy.vstack((dist, idx)).T:
        i2 = int(i2)
        t2 = tmplt_bank[i2]

        o12, _, _ = lalsimutils.overlap(h1, t2, ovrlp, delta_f, f_low, args.approximant1, args.approximant2, t1_norm=h1_norm)
        ovrlps.append(o12)

        if args.too_verbose:
            print d, t2.mass1, t2.mass2, t2.mchirp, t2.eta, o12

    opts = amrlib.apply_inv_transform(pts[idx][0], intr_prms, "mchirp_eta")
    opts = numpy.vstack((opts.T, ovrlps, idx[0]))

    ovrlps = {"mass1": opt[0], "mass2": opt[1], "overlap": [list(a) for a in opts.T]}
    with open(fname, "w") as fout:
        json.dump(ovrlps, fout)

with open("tmplt_bank.json", "w") as fout:
    json.dump(toc, fout)
