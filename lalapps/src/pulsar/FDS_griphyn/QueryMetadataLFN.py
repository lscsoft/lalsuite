#!/usr/bin/env python2

"""
Query a metadata catalog using metadata for the list of 
relevant logical filenames (LFNs).

Presently this script assumes we are querying version
two of MCS from the Globus team at ISI using the MCS tools,
and that the MCS client tools are already configured to
point to the correct MCS instance.
"""

versionN = 1.0

import os
import sys
import getopt
import tempfile
import popen2
import time

def usage():
        """
        Print a usage message to stderr.
        """

        msg = """\
NAME
        QueryMetadataLFN

SYNOPSIS
        QueryMetadataLFN --calibration=NAME --instrument=NAME
                --gps-start-time=GPS --gps-end-time=GPS
                --calibration-version=NAME [ --output=FILE ]

        QueryMetadataLFN --version

        QueryMetadataLFN --help

DESCRIPTION
        Query a Globus Metadata Catalog Service (MCS) to obtain
        the list of logical filenames (LFNs) described by
        the calibration, instrument, GPS start and end time, and
        calibration version given as inputs. Write the list of
        LFNs to stdout or to a file.

        -c, --calibration
                calibration tag for the data in the files

        -i, --instrument
                instrument or observatory that generated the data
                such as H1, H2, L1, G

        -s, --gps-start-time
                start of the GPS time range for the data

        -e, --gps-end-time
                end of the GPS time range for the data

        -v, --calibration-version
                version tag of the calibrated data 

        -o, --output
                file into which to write the list of LFNs
                stout is used by default

        -V, --version
                print version number of this script to stderr and exit

        -h, --help
                print this usage message

EXAMPLE

$ QueryMetadataLFN --calibration=Funky --instrument=H1 
        --gps-start-time=729277151 --gps-end-time=729278951 
        --calibration-version=3

"""

        print >>sys.stderr, msg


longopt = [
        "calibration=",
        "instrument=",
        "gps-start-time=",
        "gps-end-time=",
        "calibration-version=",
        "output=",
        "version",
        "help"
        ]

shortopt = "c:i:s:e:v:Vho:"


try:
        opts, args = getopt.getopt(sys.argv[1:], shortopt, longopt)
except getopt.GetoptError:
        print >>sys.stderr, "Error parsing command line"
        print >>sys.stderr, "Use --help for usage"
        sys.exit(1)

# set default values
calibration = None
instrument = None
gpsStart = None
gpsEnd = None
calversion = None
outputPath = None


for o, a in opts:
        if o in ("-h", "--help"):
                usage()
                sys.exit(0)
        elif o in ("-c", "--calibration"):
                calibration = a
        elif o in ("-i", "--instrument"):
                instrument = a
        elif o in ("-s", "--gps-start-time"):
                gpsStart = a
        elif o in ("-e", "--gps-end-time"):
                gpsEnd = a
        elif o in ("-v", "--calibration-version"):
                calversion = a
        elif o in ("-o", "--output"):
                outputPath = a
        elif o in ("-V", "--version"):
                print >>sys.stderr, versionN
                sys.exit(0)

         
# sanity checking on inputs
if not calibration:
        print >>sys.stderr, "Calibration must be specified"
        print >>sys.stderr, "Use --help for usage"
        sys.exit(1)

if not instrument:
        print >>sys.stderr, "Instrument must be specified"
        print >>sys.stderr, "Use --help for usage"
        sys.exit(1)

if instrument not in ('H1', 'H2', 'L1', 'G'):
        print >>sys.stderr, "Insrument %s not recognized by this script" % instrument
        print >>sys.stderr, "Use --help for usage"
        sys.exit(1)

if not gpsStart:
        print >>sys.stderr, "GPS start time must be specified"
        print >>sys.stderr, "Use --help for usage"
        sys.exit(1)

try:
        a = int(gpsStart)
        assert len(gpsStart) == 9
except:
        print >>sys.stderr, "GSP start time must be a 9 digit number"
        print >>sys.stderr, "Use --help for usage"


if not gpsEnd:
        print >>sys.stderr, "GPS end time must be specified"
        print >>sys.stderr, "Use --help for usage"
        sys.exit(1)
try:
        a = int(gpsEnd)
        assert len(gpsEnd) == 9
except:
        print >>sys.stderr, "GSP end time must be a 9 digit number"
        print >>sys.stderr, "Use --help for usage"


if not calversion:
        print >>sys.stderr, "Calibration version must be specified"
        print >>sys.stderr, "Use --help for usage"
        sys.exit(1)


# make sure that JAVA_HOME is in our environment
# currently necessary for MCS client tools

try:
        JAVA_HOME = os.environ["JAVA_HOME"]
except:
        print >>sys.stderr, "JAVA_HOME not defined in environment"
        sys.exit(1)
        

# make sure that we can find the MCS client tools
# and that MCS_HOME is set, and also set CLASSPATH appropriately
# for MCS tools

try:
        MCS_HOME = os.environ["MCS_HOME"]
except:
        print >>sys.stderr, "MCS_HOME not definied in environment"
        sys.exit(1)


queryPath = MCS_HOME + "/bin/query"

if not os.access(queryPath, os.X_OK):
        print >>sys.stderr, "Unable to execute %s" % queryPath
        sys.exit(1)


jarFileList = [ "%s/lib/%s" % (MCS_HOME, s) for s in os.listdir("%s/lib" % MCS_HOME)]


try:
        CLASSPATH=os.environ["CLASSPATH"]
except:
        CLASSPATH="."

for path in jarFileList:
        CLASSPATH += ":" + path

os.environ["CLASSPATH"] = CLASSPATH

# prepare a tempory file with the necessary attribute search
# syntax

try:
        fpath = tempfile.mktemp()
        f = open(fpath, "w")

        template = """\
calibration&string&=&%s
instrument&string&=&%s
gpsStart&integer&<=&%s
gpsEnd&integer&>=&%s
version&integer&=&%s
"""

        text = template % (calibration, instrument, gpsEnd, gpsStart, calversion)

        f.write(text)

        f.close()
except Exception, e:
        print >>sys.stderr, "Error: %s" % e
        sys.exit(1)

# query MCS for the list of LFNs and send the list to 
# a file or stdout

template = "%s/bin/query -l -f %s" 
cmd = template % (MCS_HOME, fpath)

capturestderr = 1
largebuf = 1024*1024
job = popen2.Popen3(cmd, capturestderr, largebuf)

ret = job.poll()

while ret == -1:
        time.sleep(1)
        ret = job.poll()

if ret != 0:
        errorList = job.childerr.readlines()
        print >>sys.stderr, "Error querying MCS"
        for s in errorList:
                print >>sys.stderr, s
        sys.exit(1)

outputList = job.fromchild.readlines()

# the last line is a summary we can delete
del outputList[-1]

if not outputPath:
        outf = sys.stdout
else:
        try:
                outf = open(outputPath, "w")
        except:
                print >>sys.stderr, "Error opening %s for writing" % outputPath

for s in outputList:
        print >>outf, s.strip()

# try to cleanup nicely
try:
        os.unlink(fpath)
except:
        pass

# exit nicely
sys.exit(0)


