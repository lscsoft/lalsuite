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

import os, sys
##from pylab import *
from numpy import *
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
def parse_octave_metrics ( octstring, metric_name="g_ij" ):
    """Parse octave-format metric of given name from the given octave-string, as written by FstatMetric[_v2]
    Returns a square array containing the metric, raises NameError if that metric_name could not be parsed.
    """

    content = octstring.replace("];","]").replace(";", ",").replace("%", "#").lstrip()

    exec ( content )	## this provides python *lists* for all metrics found in the file

    stmt = "g_ij = ConvertList2SquareMatrix( %s )" % metric_name

    try:
        exec(stmt)
    except NameError:
        print ("Failed to read metric named '%s' from input octave string." % metric_name )
        raise

    return g_ij


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
    - maxrel: is the largest relative difference in any component.
    - rel2norm: is the 2-norm of the difference-matrix divided by the 2-norm of the second matrix
    - reldet: is the relative difference in determinants
    """

    dg_ij = g1_ij - g2_ij

    ## Note: for the relative component-error, we need to be careful if
    ## any component is exactly zero: in that case we use the absolute error
    ## instead:
    inds0 = ( g2_ij == 0 )
    g2_ijRel = g2_ij.copy()
    g2_ijRel[inds0] = 1

    maxrel =  abs ( dg_ij / g2_ijRel ).max()

    rel2norm = sqrt ( linalg.norm( dg_ij, 2 ) / linalg.norm ( g2_ij, 2 ) );

    detg1 = linalg.det(g1_ij)
    detg2 = linalg.det(g2_ij)

    if debug:
        print ("det(g1) = %g, det(g2) = %g" % ( detg1, detg2 ) )
        print ("cond(g1)= %g, cond(g2) = %g" % ( linalg.cond(g1_ij), linalg.cond(g2_ij) ) )

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

    g0_ij = parse_octave_metrics ( stdout )
    if debug: print "getMetric output:g_ij =\n%s" % str(g0_ij)

    ## ----- run FstatMetric
    args1 = common_args.copy()
    args1["IFOs"] = firstIFO
    args1["outputMetric"] = outfile1
    args1["unitsType"] = 0	## SI units

    if mettype == "PHASE":
        args1["metricType"] = 2	## full-motion numerical phase metric
        v1Name = "gPh_ij"
    elif mettype == "PTOLE":
        args1["metricType"] = 4	## numerical Ptole-metric
        v1Name = "gPtole_ij"
    else:
        print ("Invalid metric type '%s' encountered ..." % mettype )
        sys.exit(1)


    (stdout, stderr) = run_code ( code1, args1 )

    octstr = load_file ( outfile1 )

    g1_ij =  parse_octave_metrics ( octstr, v1Name )
    if debug: print "FstatMetric output: %s =\n%s" % ( v1Name, str(g1_ij) )


    ## ----- run FstatMetric_v2
    args2 = common_args.copy()
    args2["IFOs"] = firstIFO
    args2["metricType"] = 0	## full-motion numerical phase metric
    args2["outputMetric"] = outfile2
    args2["coords"] = coords

    if mettype == "PHASE":
        args2["detMotionType"] = 0	## full ephemeris-based spin+orbit motion
    elif mettype == "PTOLE":
        args2["detMotionType"] = 3	## spin + Ptole-orbit detector motion
    else:
        print ("Invalid metric type '%s' encountered ..." % mettype )
        sys.exit(1)

    (stdout, stderr) = run_code ( code2, args2 )

    octstr = load_file ( outfile2 )

    g2_ij =  parse_octave_metrics ( octstr, "g_ij" )
    if debug: print "FstatMetric_v2 output:g_ij = \n%s" % str(g2_ij)

    ## ---------- compare metrics against each other:
    relerr01 = compare_metrics ( g0_ij,  g1_ij )
    relerr02 = compare_metrics ( g0_ij,  g2_ij )
    relerr12 = compare_metrics ( g1_ij,  g2_ij )

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

gF1_ij =  parse_octave_metrics ( octstr, "gF_ij" )
if debug: print "FstatMetric output: gF_ij = \n%s" % str(gF1_ij)

## ----- run FstatMetric_v2
args2 = common_args.copy()
args2["IFOs"] = options.IFOs
args2["IFOweights"] = options.IFOweights
args2["metricType"] = 1	## full-motion numerical F-stat metric
args2["outputMetric"] = outfile2
args2["coords"] = "Freq,Alpha,Delta,f1dot"

(stdout, stderr) = run_code ( code2, args2 )

octstr = load_file ( outfile2 )

gF2_ij =  parse_octave_metrics ( octstr, "gF_ij" )
if debug: print "FstatMetric_v2 output: g_ij = \n%s" % str(gF2_ij)


## ---------- compare F-metrics against each other:
relerr12 = compare_metrics ( gF1_ij, gF2_ij )
print "relerr     = (       maxrel,              rel2norm,               reldet )"
print "relerr 1-2 = " + str(relerr12)

tolF = 0.10
if ( relerr12[0] > tolF ):
    print ("\nRelative difference 'maxrel' in matrix-components exceeded tolerance %g!\n" % tolF );
    sys.exit(1)


## ========== 3) compare FstatMetric_v2 vs analytic solutions
coords = "Freq_Nat,f1dot_Nat,f2dot_Nat,f3dot_Nat"
print """
## --------------------------------------------------------------------------------
## Comparing fkdot PHASE-METRICS between [2]FstatMetric_v2 and [3]analytic solution
## using coordinates '%s', startTime=%d, duration=%d
## --------------------------------------------------------------------------------""" \
    % (coords, options.startTime, options.duration )

## ----- run FstatMetric_v2 with refTime == startTime
args2 = common_args.copy()
args2["IFO"] = firstIFO
args2["metricType"] = 0	## only compute phase-metric
args2["outputMetric"] = outfile2
args2["coords"] = coords

args2["refTime"] = 0	## ==> refTime == startTime
(stdout, stderr) = run_code ( code2, args2 )
octstr = load_file ( outfile2 )
gStart2_ij =  parse_octave_metrics ( octstr, "g_ij" )
if debug: print "refTime=startTime: FstatMetric_v2 output: g_ij = \n%s" % str(gStart2_ij)

## ----- run FstatMetric_v2 with refTime == mid-time of observation
args2["refTime"] = -1
(stdout, stderr) = run_code ( code2, args2 )
octstr = load_file ( outfile2 )
gMid2_ij =  parse_octave_metrics ( octstr, "g_ij" )
if debug: print "refTime=midTime: FstatMetric_v2 output: g_ij = \n%s" % str(gMid2_ij)

## analytic spin-metric for comparison
gStart3_ij = matrix ( [ [ 1.0/12, 1.0/12, 3.0/40, 1.0/15 ], \
                      [ 1.0/12, 4.0/45, 1.0/12, 8.0/105], \
                      [ 3.0/40, 1.0/12, 9.0/112,3.0/40 ], \
                      [ 1.0/15, 8.0/105,3.0/40,16.0/225]  \
                          ] );
if debug: print "refTime=startTime: analytic spin-metric: g_ij = \n%s" % str(gStart3_ij)

gMid3_ij = matrix ( [ [ 1.0/12,       0,   1.0/80,       0 ],  \
                    [      0, 1.0/180,        0,  1.0/840],	\
                    [ 1.0/80,       0,  1.0/448,       0 ],	\
                    [      0, 1.0/840,        0, 1.0/3600]   \
                        ] );
if debug: print "refTime=midTime: analytic spin-metric: g_ij = \n%s" % str(gMid3_ij)

## compare metrics
relerrStart_23 = compare_metrics ( gStart2_ij, gStart3_ij )
relerrMid_23 = compare_metrics ( gMid2_ij, gMid3_ij )
print "relerr          = (       maxrel,              rel2norm,               reldet )"
print "relerrStart 2-3 = " + str(relerrStart_23)
print "relerrMid 2-3   = " + str(relerrMid_23)

tolSpin = 0.01
if ( relerrStart_23[0] > tolSpin or relerrMid_23[0] > tolSpin ):
    print ("\nRelative difference 'maxrel' in matrix-components exceeded tolerance %g!\n" % tolSpin );
    sys.exit(1)


print ""
sys.exit(0)
