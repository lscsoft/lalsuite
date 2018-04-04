#!/usr/bin/env python
# read in the summary files and make coincidences etc.
from glue.ligolw import ligolw, lsctables, utils
from pylal.tools import XLALCalculateEThincaParameter
from glue import segments, lal, git_version
from glue.ligolw.utils import process
import copy
import lal
import numpy
import math
import glob
import optparse
import sys

def max_dt (row, eMatch):
    # Calculate the timing uncertainty + earth travel time
    # We use this as a conservative (but inexpensive to calculate)
    # ethinca test.
    a11 = row.Gamma0 / eMatch;
    a12 = row.Gamma1 / eMatch;
    a13 = row.Gamma2 / eMatch;
    a22 = row.Gamma3 / eMatch;
    a23 = row.Gamma4 / eMatch;
    a33 = row.Gamma5 / eMatch;

    x = (a23 * a23 - a22 * a33) * a22;
    denom = (a12*a23 - a22*a13) * (a12*a23 - a22*a13) - (a23*a23 - a22*a33) * (a12*a12 - a22*a11);

    return math.sqrt( x / denom ) + 2. * lal.REARTH_SI / lal.C_SI;

def find_slide_coincs(htrigs, ltrigs, min_mchirp, max_mchirp,
                      ethinca, slide, threshold):
    """ Calculate coincs from single detector triggers. We window based on
    mchirp, apply a new snr threshold, and only keep triggers that lie on a
    slide boundary.
    """

    coinc_sngls = lsctables.New(lsctables.SnglInspiralTable)
    num_trig = 0
    t_snrsq = float(threshold) ** 2

    # Collect the mchirp, new snr, and end time of the second detector's
    # triggers into numpy arrays for faster math operations later.
    l_ends = []
    lt_mc = []
    lt_sn = []
    lt_en = []
    for l in ltrigs:
        lt_mc.append(l.mchirp)
        lt_sn.append(l.get_new_snr()**2)
        lt_en.append(float(l.get_end()))
    lt_mc = numpy.array(lt_mc)
    lt_sn = numpy.array(lt_sn)
    lt_en = numpy.array(lt_en)

    # The list of index location of "good" triggers in the second detector
    # good triggers are ones that might possibly still form coincs with the
    # first detector's triggers. Initially, this is all of them.
    loc = numpy.arange(0, len(ltrigs))

    # Sort the first detector's triggers by newsnr ascending
    stats = []
    htrigs.sort(lambda a, b: cmp(a.get_new_snr(), b.get_new_snr()), reverse=True)

    # Iterate over each of the first detector's triggers and calculate what
    # coincs it will generate with the second detector's triggers
    for h in htrigs:
        h_end = float(h.get_end())
        h_mchirp = h.mchirp
        h_max_td = max_dt(h, ethinca)
        h_snrsq = h.get_new_snr() ** 2
        r_snrsq = t_snrsq - h_snrsq
        lind = 0

        # Grab only the information about triggers that could still form
        # coincs
        lt_snl = lt_sn[loc]
        lt_mcl = lt_mc[loc]
        lt_enl = lt_en[loc]

        # Determine the time difference between these triggers and the
        # considered "h" trigger
        num_slides = numpy.round((h_end - lt_enl)/slide)
        lt_en_tmp = lt_enl + num_slides * slide
        td = abs(h_end - lt_en_tmp)

        # Get location of coincs that are above the coinc new snr threshold
        # Triggers that fall below this threshold are never considered again
        goods = (lt_snl >= r_snrsq)

        # Get the location of coincs that pass the rough time coincidence
        # and are within the mchirp bin. We calculate ethinca only for these.
        good = goods & (td <=h_max_td) & (lt_mcl >= 2*min_mchirp - h_mchirp) & (lt_mcl <= 2*max_mchirp - h_mchirp)

        # For each remaining trigger we calculate the ethinca test and add
        # to the list of coincs if they pass
        for lind, ns in zip(loc[good], num_slides[good]):

            # don't allow zerolag or offsets too close to it!
            if abs(slide * ns) < 0.1:
                continue

            # Get the single inspiral trigger and apply the time offset
            # so that we can perform coincidence
            l = ltrigs[lind]
            l_end_old = l.get_end()
            l_end_tmp = l_end_old + float( slide * ns )

            l.set_end(l_end_tmp)
            epar = XLALCalculateEThincaParameter(h, l)

            # If we pass coincidence we add the first and second detector
            # trigger pair to the output list
            if epar < ethinca:
                # This should never happen
                if abs(slide*ns) < 0.1:
                    print 'I USED ZEROLAG!!'
                    exit(1)
                hcopy = copy.deepcopy(h)
                l.set_end(h.get_end())
                lcopy = copy.deepcopy(l)
                hcopy.event_id = lsctables.SnglInspiralID(num_trig)
                lcopy.event_id = lsctables.SnglInspiralID(num_trig)
                num_trig += 1
                coinc_sngls.append(hcopy)
                coinc_sngls.append(lcopy)
            l.set_end(l_end_old)

        # Finally set the list of second triggers to consider to the reduced
        # set that are still above the coinc new snr threshold
        loc = loc[goods]

        # No more triggers to consider
        if len(loc) == 0:
            break
    return( coinc_sngls )

def parse_mchirp_bins(param_range):
    slots = param_range.split(";")
    bins = []
    for slot in slots:
        minb, maxb = slot[1:-1].split(',')
        bins.append((float(minb), float(maxb)))
    return bins

def get_mchirp_bin(mchirp, mchirp_bins):
    for mchirp_low, mchirp_high in mchirp_bins:
        if mchirp >= mchirp_low and mchirp <= mchirp_high:
            return mchirp_low, mchirp_high
    return None, None

def get_veto_coincs(coincs, veto_start, veto_end):
    # Create a segment to veto the time given
    veto_seg = segments.segmentlist()
    veto_seg.append(segments.segment(veto_start, veto_end))

    # Remove the events around that time
    veto_coincs = coincs.veto(veto_seg)
    return veto_coincs

def get_background_livetime(lt1, lt2, slide_step, veto_window=0):
    """
    Returns the expected total time-slid coincident livetime for two detectors
    having single-ifo livetimes lt1, lt2, when doing all possible slides with
    time offsets being multiples of slide_step.
    Specifying veto_window > 0 removes a given duration from each of the
    single-ifo livetimes.

    Derivation: Consider an analysis of length T with two ifos having duty
    factors d1, d2. Then the expected coinc livetime per slide is T*d1*d2.
    The maximum possible number of slides is T/s where s is the step. The
    total time is the product which can be expressed as
    (T*d1)(T*d2)/s = lt1*lt2/s.
    """
    lt1 = float(lt1)
    lt2 = float(lt2)
    slide_step = float(slide_step)
    veto_window = float(veto_window)
    return (lt1 - veto_window*2) * (lt2 - veto_window*2) / slide_step

