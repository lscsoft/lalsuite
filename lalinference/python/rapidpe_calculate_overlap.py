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
lsctables.use_in(ligolw.LIGOLWContentHandler)
from pylal.series import read_psd_xmldoc, LIGOLWContentHandler

# Adapted from similar code in gstlal.cbc_template_fir
def generate_waveform_from_tmplt_old(tmplt, approximant, delta_f=0.125, sample_rate=16384, f_low=40, amporder=-1, phaseorder=7):
    hfp, hfx = lalsimulation.SimInspiralChooseFDWaveform(
        0., # phase
        delta_f,
        lal.MSUN_SI * tmplt.mass1, lal.MSUN_SI * tmplt.mass2,
        tmplt.spin1x, tmplt.spin1y, tmplt.spin1z,
        tmplt.spin2x, tmplt.spin2y, tmplt.spin2z,
        f_low,
        2048.0, #FIXME
        0, #FIXME chosen until suitable default value for f_ref is defined
        1.e6 * lal.PC_SI, # distance
        0, # inclination
        0., # tidal deformability lambda 1
        0., # tidal deformability lambda 2
        None, # waveform flags
        None, # Non GR params
        amporder, phaseorder,
        lalsimulation.GetApproximantFromString(str(approximant))
    )
    hfp.data.data += 1j*hfx.data.data
    return hfp

# FIXME: Check for discrepancies between this and the old waveform generator
def generate_waveform_from_tmplt(tmplt, approximant, delta_f=0.125, f_low=40, amporder=-1, phaseorder=7):

    params = lalsimutils.ChooseWaveformParams(
        deltaF = delta_f,
        m1 = lal.MSUN_SI * tmplt.mass1, m2 = lal.MSUN_SI * tmplt.mass2,
        s1x = tmplt.spin1x, s1y = tmplt.spin1y, s1z = tmplt.spin1z,
        s2x = tmplt.spin2x, s2y = tmplt.spin2y, s2z = tmplt.spin2z,
        fmin = f_low, fref = 0,
        dist = 1.e6 * lal.PC_SI, # distance
        ampO = amporder, phaseO = phaseorder,
        approx = lalsimulation.GetApproximantFromString(str(approximant)),
        taper = lalsimulation.SIM_INSPIRAL_TAPER_START # FIXME: don't hardcode
    )
    return lalsimutils.hoff(params, Fp=1, Fc=1)

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

    def anon_interp(fvals):
        return numpy.interp(fvals, f, psd)
    return numpy.array(map(anon_interp, f))

argp = argparse.ArgumentParser()
argp.add_argument("-s", "--tmplt-start-index", type=int, help="Start at this index of the template bank.")
argp.add_argument("-e", "--tmplt-end-index", type=int, help="End at this index of the template bank.")
argp.add_argument("-t", "--tmplt-bank-file", help="File name of the template bank. Required.")
argp.add_argument("-d", "--distance-coordinates", default="mchirp_eta", help="Coordinate system in which to calculate 'closness'. Default is mchirp_eta.")
argp.add_argument("-p", "--psd-file", help="Name of PSD XML file. Required.")
argp.add_argument("-f", "--f-low", type=float, default=40., help="Lowest frequency component of template. Default is 40 Hz.")
argp.add_argument("-F", "--delta-f", type=float, default=0.125, help="Frequency binning of the FD waveform. Default is 0.125.")
argp.add_argument("-a", "--approximant1", default="TaylorF2", help="Approximant to use for target waveform. Default is TaylorF2.")
argp.add_argument("-b", "--approximant2", default="TaylorF2", help="Approximant to use for overlapped waveform. Default is TaylorF2.")
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

# lalapps_tmpltbank assigns 0 ID to all events, so we remap
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

# FIXME: add psd
ovrlp = lalsimutils.Overlap(fLow=args.f_low, fMax=2000, deltaF=delta_f, psd=psd, analyticPSD_Q=False)

idx_range = range(args.tmplt_start_index or 0, args.tmplt_end_index or len(tmplt_bank))

# FIXME:
npts = len(tmplt_bank)
#npts = 100
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

    t1 = tmplt_bank[i1]
    h1 = generate_waveform_from_tmplt(t1, args.approximant1, delta_f, args.f_low)
    h1_norm = ovrlp.norm(h1)
    dist, idx = tree.query(pt, k=npts, return_distance=True)
    print "--- (%f, %f) / (%f, %f)" % (t1.mass1, t1.mass2, t1.mchirp, t1.eta)

    ovrlps = []
    for d, i2 in numpy.vstack((dist, idx)).T:
        i2 = int(i2)
        t2 = tmplt_bank[i2]
        h2 = generate_waveform_from_tmplt(t2, args.approximant2, delta_f, args.f_low)
        h2_norm = ovrlp.norm(h2)
        o12 = ovrlp.ip(h1, h2) / h1_norm / h2_norm
        #print d, t2.mass1, t2.mass2, t2.mchirp, t2.eta, o12
        ovrlps.append(o12)

    opts = amrlib.apply_inv_transform(pts[idx][0], intr_prms, "mchirp_eta")
    opts = numpy.vstack((opts.T, ovrlps, idx[0]))

    ovrlps = {"mass1": opt[0], "mass2": opt[1], "overlap": [list(a) for a in opts.T]}
    with open(fname, "w") as fout:
        json.dump(ovrlps, fout)

with open("tmplt_bank.json", "w") as fout:
    json.dump(toc, fout)
