import sys
from getopt import getopt
from ClusterComputeFUtils import *

# Constants
basename = "nbSFT"

# Parameter defaults
start = 0
end = 86400
instrument = 1
frequency = 1200.0
bandwidth = 1.0
spindown  = 0.0
spinband  = 0.0
metric = 1
mismatch = 0.02
threshold = 10.0
liststart = 0
num = -1

# Command-line options
short_args = "hs:e:f:b:t:i:l:n:p:"
long_args = [ "help", "start=", "end=", "frequency=", "bandwidth=", \
	      "threshold=", "instrument=", "liststart=", "num=", \
	      "spindown=" ]
usage = "\nusage: python %s [options] listfile\n\n" % sys.argv[0] \
	+ "arguments:\n" \
	+ "  listfile                list of RA and DEC pairs to search\n\n" \
	+ "options:\n" \
	+ "  -h, --help              print this help message\n" \
	+ "  -s, --start=START       start at GPS time START (s) [0]\n" \
	+ "  -e, --end=END           end at GPS time END (s) [86400]\n" \
	+ "  -i, --instrument=CODE   0 = GEO, 1 = LLO 2 = LHO, 3 = Rome\n" \
	+ "  -f, --frequency=FREQ    start at frequency FREQ (Hz) [1200.0]\n" \
	+ "  -b, --bandwidth=BAND    cover frequency bandwidth BAND (Hz) [1.0]\n" \
	+ "  -d, --spindown=SPIN     search up from SPIN spindown (Hz^2) [0.0]\n" \
	+ "  -p, --spinband=SBAND    search up to SPIN+SBAND (Hz^2) [0.0]\n" \
	+ "  -m, --metric=MCODE      0 = none, 1 = PtoleMetric, 2 = CoherentMetric [1]\n" \
	+ "  -x, --mismatch=MAX      maximum mismatch is MAX [0.02]\n" \
	+ "  -t, --threshold=THRESH  set F-statistic threshold to THRESH [10.0]\n" \
	+ "  -l, --liststart=LSTART  skip first LSTART points in list [0]\n" \
	+ "  -n, --num=NUM           search NUM points from list [-1]\n"
options, args = getopt( sys.argv[1:], short_args, long_args )
for opt in options:
	if opt[0] == "-h" or opt[0] == "--help":
		print usage
		sys.exit( 0 )
	elif opt[0] == "-s" or opt[0] == "--start":
		start = long( opt[1] )
	elif opt[0] == "-e" or opt[0] == "--end":
		end = long( opt[1] )
	elif opt[0] == "-i" or opt[0] == "--instrument":
		instrument = int( opt[1] )
	elif opt[0] == "-f" or opt[0] == "--frequency":
		frequency = double( opt[1] )
	elif opt[0] == "-b" or opt[0] == "--bandwidth":
		bandwidth = double( opt[1] )
	elif opt[0] == "-d" or opt[0] == "--spindown":
		spindown = float( opt[1] )
	elif opt[0] == "-p" or opt[0] == "--spinband":
		spinband = float( opt[1] )
	elif opt[0] == "-m" or opt[0] == "--metric":
		metric = int( opt[1] )
	elif opt[0] == "-x" or opt[0] == "--mismatch":
		mismatch = float( opt[1] )
	elif opt[0] == "-t" or opt[0] == "--threshold":
		threshold = float( opt[1] )
	elif opt[0] == "-l" or opt[0] == "--liststart":
		liststart = int( opt[1] )
	elif opt[0] == "-n" or opt[0] == "--num":
		num = int( opt[1] )
	else:
		print "%s: Unrecognized option %s" % ( sys.argv[0], opt[0] )
		print "%s --help          for usage information" % sys.argv[0]
		sys.exit( 2 )
if len( args ) < 1:
	print "%s: Missing argument" % sys.argv[0]
	print "%s --help          for usage information" % sys.argv[0]
	sys.exit( 2 )



### Get names of cluster shared directory and node local temp directory
sharedDir = "/home/teviet/scratch/"
localDir = "/home/teviet/tmp/"


### Create script to query metadata catalogue for names of SFT files
### for timespan and instrument in question.

### Create script to query local data replicator for locations of
### local SFTs and/or URLs of remote SFTs.

### Create script to move remote SFTs to sharedDir, and make list of
### all SFTs.
sftList = "sftList.txt"


### Get name of directory where ephemerides are stored, and year of ephemeris
ephemDir = "."
year = getyearstr( ephemDir, start, end )
if year == "":
	print "%s - no ephemerides for GPS times %i-%i" % ( sys.argv[0], start, end )
	sys.exit(1)

# Create script to distribute data to nodes
fp = open( "copyData.sh", "w" )
fp.write( "#! /bin/sh\n" )
fp.write( "cp " + sharedDir + basename + "* " + localDir + "\n" )
fp.close()

# Create script to return results from nodes
fp = open( "returnData.sh", "w" )
fp.write( "#! /bin/sh\n" )
fp.write( "rm " + localDir + basename + "*\n" )
fp.write( "mv " + localDir + "*.out" + " " + sharedDir + "\n" )
fp.close()

# Create SFT extraction Condor script
fp = open( "extractSFTband.sub", "w" )
fp.write( "executable = extractSFTband\n" )
fp.write( "arguments =" )
fp.write( " -i " + sftList )
fp.write( " -o " + basename )
fp.write( " -f %f" % frequency )
fp.write( " -b %f" % bandwidth )
fp.write( "\n" )
fp.write( "universe = standard\n" )
fp.write( "notification = Never\n" )
fp.write( "log = ClusterComputeF.log\n" )
fp.write( "queue\n" )
fp.close()

# Create F-statistic computation Condor script
fp = open( "ComputeFStatistic.sub", "w" )
fp.write( "executable = ComputeFStatistic\n" )
fp.write( "arguments =" )
fp.write( " -I %i" % instrument )
fp.write( " -f %f" % frequency )
fp.write( " -b %f" % bandwidth )
fp.write( " -s %f" % spindown )
fp.write( " -m %f" % spinband )
fp.write( " -M %i" % metric )
fp.write( " -X %f" % mismatch )
fp.write( " $(area)" )
fp.write( " -D " + localDir )
fp.write( " -i " + basename )
fp.write( " -E " + ephemDir )
fp.write( " -y " + year )
fp.write( " -F %f" % threshold )
fp.write( "\n" )
fp.write( "universe = standard\n" )
fp.write( "notification = Never\n" )
fp.write( "output = $(number).out\n" )
fp.write( "log = ClusterComputeF.log\n" )
fp.write( "queue\n" )
fp.close()

# Get list of areas to search
areas = []
fpList = open( args[0], "r" )
if liststart < 0:
	liststart = 0
for i in range( 0, liststart ):
	fpList.readline()
if num > 0:
	for i in range( 0, num ):
		line = fpList.readline()
		if len( line ) > 0:
			areas = areas + [line]
		else:
			i = num
	if num > len( areas ):
		print "%s: Only %i areas read from %s" % \
		      ( sys.argv[0], len( areas ), args[0] )
else:
	i = 1
	while i:
		line = fpList.readline()
		if len( line ) > 0:
			areas = areas + [line]
		else:
			i = 0
fpList.close()
num = len( areas )

# Create DAG
fp = open( "ClusterComputeF.dag", "w" )
### DAG jobs to migrate data here?
fp.write( "JOB extract extractSFTband.sub\n" )
for i in range( liststart, liststart + num ):
	fp.write( "JOB compute%i ComputeFStatistic.sub\n" % i )
for i in range( liststart, liststart + num ):
	fp.write( "SCRIPT PRE  compute%i copyData.sh\n" % i )
	fp.write( "SCRIPT POST compute%i returnData.sh\n" % i )
fp.write( "\nPARENT extract CHILD " )
for i in range( liststart, liststart + num ):
	fp.write( "compute%i " % i )
fp.write( "\n\n" )
for i in range( 0, num ):
	fp.write( "VARS compute%i" % ( i + liststart ) )
	fp.write( " number=\"%i\"" % ( i + liststart ) )
	fp.write( " area=\"" )
	for c in areas[i]:
		if c == '"':
			fp.write( "\\\"" )
		elif c == '\n':
			pass
		else:
			fp.write( c )
	fp.write( "\"\n" )
#	fp.write( " area=\"%s\"" % areas[i] )
fp.close()

# Execute Condor job
del areas
print "%s: Submitting jobs %i through %i" % \
      ( sys.argv[0], liststart, liststart + num - 1 )
#execlp( "condor_submit_dag", "ClusterComputeF.dag" )
