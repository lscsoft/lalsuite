#!/usr/bin/env python

##
## Copyright (C) 2009 Reinhard Prix
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with with program; see the file COPYING. If not, write to the
##  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
##  MA  02111-1307  USA
##

## check metric codes 'getMetric' 'FstatMetric' and 'FstatMetric_v2' by
## comparing them against each other.
## Given that they represent 3 very different implementations of
## metric calculations, this provides a very powerful consistency test

import os
from pylab import *
from optparse import OptionParser
import subprocess

debug = False

## ----- test if LAL_DATA_PATH has been set ... needed to locate ephemeris-files
if not "LAL_DATA_PATH" in os.environ.keys():
    if "LAL_PREFIX" in os.environ.keys():
        os.environ["LAL_DATA_PATH"] = ( ".:%s/share/lal" % os.environ["LAL_PREFIX"] )
    else:
        print "Need environment-variable LAL_PREFIX, or LAL_DATA_PATH to be set"
        print "to your ephemeris-directory (e.g. /usr/local/share/lal)"
        print "This might indicate an incomplete LAL installation"
        sys.exit(1)

## ----- allow 'make test' to work from builddir != srcdir
if not "builddir" in os.environ.keys():
    builddir = "./"
else:
    builddir = os.environ["builddir"]

## ----- names of codes and input/output files
code0 = builddir + "lalapps_getMetric"
code1 = builddir + "lalapps_FstatMetric"
code2 = builddir + "lalapps_FstatMetric_v2"

## ---------- commandline input parsing ----------
parser = OptionParser()

parser = OptionParser(usage="usage: %prog [options]")

parser.add_option("--Alpha", help="skyposition Alpha in radians, equatorial coords", type="float", default=1.0 )
parser.add_option("--Delta", help="skyposition Delta in radians, equatorial coords", type="float", default=0.5 )
parser.add_option("--Freq",  help="Frequency", type="float", default=100 )
parser.add_option("--startTime",  help="GPS start time of observation", type="int", default=792576013 )
parser.add_option("--duration",   help="Duration of observation in seconds", type="float", default=60000 )
parser.add_option("--ephemYear",   help="Year-string to use as ephemeris-year identifier", type="string", default="05-09" )

parser.add_option("--IFOs",  help="Comma-separated list of detectors, eg. 'H1,L1, ...'", default="H1,L1,V1" )
parser.add_option("--IFOweights",  help="Comma-separated list of detectors, eg. 'H1,L1, ...'", default="1.0,0.9,1.2" )

parser.add_option("--h0",  help="GW amplitude h0", type="float", default=0.03 )
parser.add_option("--cosi",  help="Pulsar orientation-angle cos(iota) [-1,1]", type="float", default=-0.3 )
parser.add_option("--psi",  help="Wave polarization-angle psi [-pi/4, pi/4]", type="float", default=0.0 )
parser.add_option("--phi0",  help="GW initial phase phi_0 [0, 2pi]", type="float", default=0.0 )
parser.add_option("--coords",  help="Doppler-coordinates to compute metric in", default="Freq_Nat,Alpha,Delta" )
parser.add_option("--projection",  help="Project onto subspace orthogonal to this axis: 0=none, 1=1st-coord, 2=2nd-coord etc", type="int", default=0 )

(options, args) = parser.parse_args()


## ----- function to run lalapps codes
def run_code ( code, args, ps = [] ):
    ## run given code with commandline constructed from args-dictionary, appending optional 'ps' list argument
    ## returns tuple (retval, stdin, stderr)
    cmdline = ""
    arglist = [ code ]
    for opt in args.keys():
        val = args[opt]
        if type(val) == float:
            newarg = "--%s=%.17g" % (opt, val)
        elif type(val) == int:
            newarg = "--%s=%d" % (opt, val)
        elif type(val) == str:
            newarg = "--%s=%s" % (opt, val)
        else:
            print "Error: encountered argument '%s=%s' of unkown type '%s'." % ( opt, val, type(val) )
            sys.exit(1)	## FIXME: raise an exception

        arglist.append (newarg)

    if ps: arglist += ps

    if debug: print "Calling '%s'" % str(arglist)

    try:
        p1 = subprocess.Popen(arglist, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except OSError:
        print "Failed to launch command+args '%s'." % repr(arglist)
        print "Error-code = %d" % e
        sys.exit(1)

    (stdoutdata, stderrdata) = p1.communicate()

    retval = p1.returncode
    if retval != 0:
        print "\n\nError: program '%s' failed with exit code '%d'\n" % ( code, retval )
        print "cmdline: %s" % arglist
        print stderrdata
        raise subprocess.CalledProcessError(retval, code)

    return (stdoutdata, stderrdata)

## ----- function to convert a python list into a square matrix ('array')
def ConvertList2SquareMatrix ( slist ):
    """Converts the given list into a square matrix, return None on error"""
    mm = array ( slist )
    dim = sqrt( len(mm) )
    if ( int(dim) != dim ):
        print ("Error: input list of length %d cannot be turned into a square matrix." % len(mm) )
        return

    mm.shape = ( dim, dim )

    return mm


## ----- function to parse the octave-format metric files written by FstatMetric[_v2]
def parse_octave_metrics ( octstring ):
    """Parse octave-format metrics from the given string, as written by FstatMetric[_v2]
    This returns a dictionary with the metric names as keys
    """
    content = octstring.replace("];","]").replace(";", ",").replace("%", "#").lstrip()

    exec ( content )	## this provides python *lists* for all metrics found in the file

    ret = {}	## return-dictionary

    try:
        gPh_ij = ConvertList2SquareMatrix( gPh_ij )
        ret["gPh_ij"] = gPh_ij
    except NameError:
        pass

    try:
        g_ij = ConvertList2SquareMatrix( g_ij )
        ret["g_ij"] = g_ij
    except NameError:
        pass

    try:
        gOrb_ij = ConvertList2SquareMatrix( gOrb_ij )
        ret["gOrb_ij"] = gOrb_ij
    except NameError:
        pass

    try:
        gPtole_ij = ConvertList2SquareMatrix( gPtole_ij )
        ret["gPtole_ij"] = gPtole_ij
    except NameError:
        pass

    try:
        gF_ij   = ConvertList2SquareMatrix( gF_ij )
        gFav_ij = ConvertList2SquareMatrix( gFav_ij )
        m1_ij   = ConvertList2SquareMatrix( m1_ij )
        m2_ij   = ConvertList2SquareMatrix( m2_ij )
        m3_ij   = ConvertList2SquareMatrix( m3_ij )

        ret["gF_ij"] =  gF_ij
        ret["gFav_ij"] = gFav_ij
        ret["m1_ij"] = m1_ij
        ret["m2_ij"] = m2_ij
        ret["m3_ij"] = m3_ij
    except NameError:
        pass

    return ret


## ----- mini-function to load a text-file into a string
def load_file ( fname ):
    """mini-function to load a text-file into a string"""
    try:
        mfile = open( fname, "rb" )
        filestring = mfile.read()
        mfile.close()
    except IOError:
        print ("Error: failed to open '%s' for reading." % fname )
        sys.exit(1)

    return filestring


## ----- function to compare two metrics and return measures of aggreement
def compare_metrics ( g1_ij, g2_ij ):
    """Compare two metrics and return a measure of aggreement, namely the 3-tuple
    ( maxrel, rel2norm, reldet ), where
    - maxrel: is the largest relative difference in any component
    - rel2norm: is the 2-norm of the difference-matrix divided by the 2-norm of the second matrix
    - reldet: is the relative difference in determinants
    """

    dg_ij = g1_ij - g2_ij

    maxrel =  abs ( dg_ij / g2_ij ).max()

    rel2norm = sqrt ( norm( dg_ij, 2 ) / norm ( g2_ij, 2 ) );

    detg1 = det(g1_ij)
    detg2 = det(g2_ij)

    if debug:
        print ("det(g1) = %g, det(g2) = %g" % ( detg1, detg2 ) )
        print ("cond(g1)= %g, cond(g2) = %g" % ( cond(g1_ij), cond(g2_ij) ) )

    reldet = abs(detg1 - detg2) / abs ( detg2 )

    return ( maxrel, rel2norm, reldet )


## ----- set run parameters --------------------
outfile1 = "Fmetric.dat"
outfile2 = "Fmetric_v2.dat"

## we want to run a number of different comparisons:
## 1) compare 'phase-metrics' across all 3 codes, and
## 2) compare multi-IFO F-stat metric between FstatMetric and FstatMetric_v2

firstIFO = options.IFOs.split(",")[0]

## ========== 1) compare phase- and ptole- metrics
coords = "Freq,Alpha,Delta,f1dot"

common_args = { "Alpha" : options.Alpha,
                "Delta" : options.Delta,
                "Freq"  : options.Freq,
                "startTime" : options.startTime,
                "duration"  : options.duration,
                "ephemYear" : options.ephemYear
                }

for mettype in ["PHASE", "PTOLE"]:

    print """
## --------------------------------------------------------------------------------
## Comparing %s-METRICS between [0]getMesh, [1]FstatMetric and [2]FstatMetric_v2
## using coordinates '%s', startTime=%d, duration=%d, IFO=%s
## and Doppler-position: Freq=%f, Alpha=%f, Delta=%f
## --------------------------------------------------------------------------------""" \
        % (mettype, coords, options.startTime, options.duration, firstIFO, options.Freq, options.Alpha, options.Delta )

    ## ----- run getMetric
    args0 = common_args.copy()
    args0["IFO"] = firstIFO

    if mettype == "PHASE":
        args0["metricType"] = 3	## full-motion numerical phase metric
    elif mettype == "PTOLE":
        args0["metricType"] = 1	## analytic Ptole-metric
    else:
        print ("Invalid metric type '%s' encountered ..." % mettype )
        sys.exit(1)

    (stdout, stderr) = run_code ( code0, args0 )

    met0 = parse_octave_metrics ( stdout )
    if debug: print "getMetric output:\ng_ij = %s" % str(met0["g_ij"])

    ## ----- run FstatMetric
    args1 = common_args.copy()
    args1["IFOs"] = firstIFO
    args1["outputMetric"] = outfile1
    args1["unitsType"] = 0	## SI units

    if mettype == "PHASE":
        args1["metricType"] = 2	## full-motion numerical phase metric
    elif mettype == "PTOLE":
        args1["metricType"] = 4	## numerical Ptole-metric
    else:
        print ("Invalid metric type '%s' encountered ..." % mettype )
        sys.exit(1)


    (stdout, stderr) = run_code ( code1, args1 )

    octstr = load_file ( outfile1 )

    met1 =  parse_octave_metrics ( octstr )
    if debug: print "FstatMetric output:\ngPh_ij = %s" % str(met1["gPh_ij"])


    ## ----- run FstatMetric_v2
    args2 = common_args.copy()
    args2["IFOs"] = firstIFO
    args2["metricType"] = 0	## full-motion numerical phase metric
    args2["outputMetric"] = outfile2
    args2["coords"] = coords

    if mettype == "PHASE":
        args2["detMotionType"] = 0	## full ephemeris-based spin+orbit motion
        v1Name = "gPh_ij"
    elif mettype == "PTOLE":
        args2["detMotionType"] = 3	## spin + Ptole-orbit detector motion
        v1Name = "gPtole_ij"
    else:
        print ("Invalid metric type '%s' encountered ..." % mettype )
        sys.exit(1)

    (stdout, stderr) = run_code ( code2, args2 )

    octstr = load_file ( outfile2 )

    met2 =  parse_octave_metrics ( octstr )
    if debug: print "FstatMetric_v2 output:\ng_ij = %s" % str(met2["g_ij"])

    ## ---------- compare metrics against each other:
    relerr01 = compare_metrics ( met0["g_ij"],  met1[v1Name] )
    relerr02 = compare_metrics ( met0["g_ij"],  met2["g_ij"] )
    relerr12 = compare_metrics ( met1[v1Name],  met2["g_ij"] )

    print "relerr     = (       maxrel,              rel2norm,               reldet )"
    print "relerr 0-1 = " + str(relerr01)
    print "relerr 0-2 = " + str(relerr02)
    print "relerr 1-2 = " + str(relerr12)

    tolPh = 0.01
    if ( relerr01[0] > tolPh or relerr02[0] > tolPh or relerr12[0] > tolPh ):
        print ("\nRelative difference 'maxrel' in matrix-components exceeded tolerance %g!\n" % tolPh );
        sys.exit(1)

## end of loop over metric types

## ========== 2) compare F-metrics
coords = "Freq,Alpha,Delta,f1dot"
print """
## --------------------------------------------------------------------------------
## Comparing FSTAT-METRICS between [1]FstatMetric and [2]FstatMetric_v2
## using coordinates '%s', startTime=%d, duration=%d
## for IFOs = %s, with weights = %s
## and Doppler-position: Freq=%f, Alpha=%f, Delta=%f
## --------------------------------------------------------------------------------""" \
    % (coords, options.startTime, options.duration, options.IFOs, options.IFOweights, options.Freq, options.Alpha, options.Delta )

## ----- run FstatMetric
args1 = common_args.copy()
args1["IFOs"] = options.IFOs
args1["IFOweights"] = options.IFOweights
args1["metricType"] = 1	## full-motion numerical phase metric
args1["outputMetric"] = outfile1
args1["unitsType"] = 0	## SI units

(stdout, stderr) = run_code ( code1, args1 )

octstr = load_file ( outfile1 )

met1 =  parse_octave_metrics ( octstr )
if debug: print "FstatMetric output:\ngF_ij = %s" % str(met1["gF_ij"])

## ----- run FstatMetric_v2
args2 = common_args.copy()
args2["IFOs"] = options.IFOs
args2["IFOweights"] = options.IFOweights
args2["metricType"] = 1	## full-motion numerical F-stat metric
args2["outputMetric"] = outfile2
args2["coords"] = "Freq,Alpha,Delta,f1dot"

(stdout, stderr) = run_code ( code2, args2 )

octstr = load_file ( outfile2 )

met2 =  parse_octave_metrics ( octstr )
if debug: print "FstatMetric_v2 output:\ng_ij = %s" % str(met2["gF_ij"])


## ---------- compare F-metrics against each other:
relerr12 = compare_metrics ( met1["gF_ij"], met2["gF_ij"] )
print "relerr     = (       maxrel,              rel2norm,               reldet )"
print "relerr 1-2 = " + str(relerr12)

tolF = 0.05
if ( relerr12[0] > tolF ):
    print ("\nRelative difference 'maxrel' in matrix-components exceeded tolerance %g!\n" % tolF );
    sys.exit(1)


print ""


## ========== 3) compare FstatMetric_v2 vs analytic solutions

