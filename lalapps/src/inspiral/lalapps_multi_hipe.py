"""
multi_hipe.in - multi inspiral pipeline driver script

This script uses master segment lists to determine a set of segment
lists appropriate to running inspiral_hipe on multiple epochs of fixed
length between some start and some end time.  At present, the script
only sets up directories and writes appropriate segment files to those
directories.  

It uses the same configuration file as the inspiral_hipe script to
determine various parameters and then set up analysis and injection
runs.

Coming soon:

  1.  read in the standard hipe config file to determine parameters,
  etc.  Modify as needed for particular epoch and write back out in
  the relevant directory

  2.  set up dags for zero-lag and injection runs in each of the
  epochs. 

  3.  write a log file

  4.  write a super-dag which would submit the other sub-dags.

"""

__author__ = 'Patrick Brady <patrick@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'

##############################################################################
# import standard modules
import os
import sys 
import pylab
import random
import string
import ConfigParser
from optparse import *
from glue import segments
from glue import segmentsUtils
from pylal import readMeta

##############################################################################
# define a few utility functions that make this job easier

# remove all segments of size less than min_length
def cleanlist(seglist, min_length):
  removals = segments.segmentlist()
  for seg in seglist:
    if seg.__abs__() < min_length:
      removals.append(seg)
  seglist = seglist - removals

  return seglist

# return those segments which intersect the interval
def getSegments ( seglistin, interval):

  seglistout = segments.segmentlist([s for s in seglistin \
      if (s[1] > interval[0] and s[0] < interval[1]) ])

  return seglistout

##############################################################################
# define usage and command line options and arguments - parse
usage = """usage: %prog ...

Construct a set of intervals of fixed duration spanning some larger
epoch; determine the segments appropriate to the inspiral code for
each of these intervals and write them out to approriate files in each
of the directories.

As of now, the code appears to do the right thing.  The segments are
constructed to allow a 1 second overlap between searched times from
one interval to the next.  BEWARE: This can result is some strange
behaviour if people are not careful when generating summary files for
the whole search.

This code also uses the inspiral_hipe config file to determine
information about segments and to insure that appropriate overlaps,
etc are being done.  

In the future, it will also construct dags to cover the particular
interval. Moreover, it will also generate dags for injection runs
associated with each interval. This should speed up the whole process
of bulk processing of the data. 


DIRECTORY HIERARCHY:

The directory hierarchy that a search would have then follows:

searchdir
  tstart-tend
    analysis
    injection001
    injection002
    ....
    injection00n
  .
  .
  .


METADATA FILES:

The script currently writes out a segwizard format segment file which
contains a listing of the intervals analyzed which is named

searchdir/multi_hipe_selectedsegs.txt

This can be used to add a new feature a little later on which would
read in this segment list, use it to determine the fixed analysis
interval, and construct the next interval in the sequence, generate
the dag, and launch it.

It would also be nice to have a master database that would include
information about each of the jobs required to complete the analysis
of the given interval of data noting it's success or failure.


RELATED TOOLS AND REQUIRED TOOLS:

With a structure like this, a host of other tools can be developed to
make the whole analysis engine work well. Here is a list of things
that we need with a note about its current status:

* inspiral_hipe:  exists and meta-stable
* multi_hipe: exists, but developmental
* follow_hipe: under development
* hipecoire: exists, stable, needs refinement
* upperlimit: exists, metastable
* hipe_likelihooh: under development
* hipe_summary_page: under development
* multi_hipe_summary: does not exist

"""
parser = OptionParser( usage )

parser.add_option("-q", "--verbose", action="store_true",default=False,\
  help="make things verbose" )
parser.add_option("-H","--h1-segments",action="store",type="string",\
  default=None, metavar=" H1_SEGMENTS", help="H1 input segment to read" )
parser.add_option("-K","--h2-segments",action="store",type="string",\
  default=None, metavar=" H2_SEGMENTS", help="H2 input segment to read" )
parser.add_option("-L","--l1-segments",action="store",type="string",\
  default=None, metavar=" L1_SEGMENTS", help="L1 input segment to read" )
parser.add_option("-Z","--v1-segments",action="store",type="string",\
  default=None, metavar=" V1_SEGMENTS", help="V1 input segment to read" )
parser.add_option("-m","--start-time",action="store",type="int",\
    default=None, metavar=" START TIME",\
    help="start time of search")
parser.add_option("-e","--end-time",action="store",type="int",\
    default=None, metavar=" END TIME",\
    help="end time of search")
parser.add_option("-n","--interval",action="store",type="int",\
    default=None, metavar=" INTERVAL",\
    help="length of each interval to be analyzed separately")
parser.add_option("-y","--ninjections",action="store",type="int",\
    default=0, metavar=" NINJ",\
    help="Number of injection runs to set up")
# read in the config file
parser.add_option("-f","--multi-hipe-config-file",action="store",type="string",\
  default=None, metavar=" FILE", help="use configuration file FILE" )
#parser.add_option("-p", "--log-path",action="store",type="string",\
    #metavar=" PATH",help="directory to write condor log file")
# Add some plotting capabilities to check things
parser.add_option("-Y", "--plot-segments", action="store_true",default=False,\
  help="plot segments for each interval with original segments" )
parser.add_option("-v", "--version",action="store_true",default=False,\
    help="print version information and exit")

parser.add_option("-u", "--user-tag",action="store",type="string",\
    default=None,metavar=" USERTAG",\
    help="tag the jobs with USERTAG (overrides value in ini file)")

parser.add_option("-g", "--g1-data",action="store_true",default=False,\
    help="analyze g1 data")
parser.add_option("-a", "--h1-data",action="store_true",default=False,\
    help="analyze h1 data")
parser.add_option("-b", "--h2-data",action="store_true",default=False,\
    help="analyze h2 data")
parser.add_option("-l", "--l1-data",action="store_true",default=False,\
    help="analyze l1 data")
parser.add_option("-X", "--v1-data",action="store_true",default=False,\
    help="analyze v1 data")

parser.add_option("-S", "--one-ifo",action="store_true",default=False,\
    help="analyze single ifo data (not usable for GEO)")
parser.add_option("-D", "--two-ifo",action="store_true",default=False,\
    help="analyze two interferometer data")
parser.add_option("-T", "--three-ifo",action="store_true",default=False,\
    help="analyze three interferometer data")
parser.add_option("-Q", "--four-ifo",action="store_true",default=False,\
    help="analyze four intereferometer data")

parser.add_option("-A", "--analyze-all",action="store_true",default=False,\
    help="analyze all ifos and all data (over-rides above)")

parser.add_option("-d", "--datafind",action="store_true",default=False,\
    help="run LSCdataFind to create frame cache files")
parser.add_option("-t", "--template-bank",action="store_true",default=False,\
    help="run lalapps_tmpltbank to generate template banks")
parser.add_option("-i", "--inspiral" ,action="store_true",default=False,\
    help="run lalapps_inspiral (or lalapps_ring) to generate triggers")
parser.add_option("-r", "--ringdown" ,action="store_true",default=False,\
    help="run lalapps_ring insted of lalapps_inspiral")
parser.add_option("-c", "--coincidence",action="store_true",default=False,\
    help="run lalapps_thinca to test for coincidence")
parser.add_option("-U", "--td-follow-bank", action="store_true",default=False,\
    help="Run tmpltbank for TD follow-up")
parser.add_option("-B", "--trigbank",action="store_true",default=False,\
    help="run lalapps_trigbank for banks of coinc triggers")
parser.add_option("-V", "--inspiral-veto",action="store_true",default=False,\
    help="run lalapps_inspiral with vetos")
parser.add_option("-W", "--td-follow-inspiral",action="store_true",\
    default=False,help="run lalapps_inspiral for TD follow-up")
parser.add_option("-C", "--second-coinc" ,action="store_true",default=False,\
    help="run lalapps_thinca on the inspiral veto triggers")
parser.add_option("-j", "--coherent-bank",action="store_true",default=False,\
    help="run lalapps_coherentbank to make coherent bank")
parser.add_option("-k", "--coherent-inspiral",action="store_true",
    default=False,help="run lalapps_coherent_inspiral for coherent analysis")
parser.add_option("-s", "--sire",action="store_true",default=False,\
    help="do sires to sweep up triggers")

parser.add_option("-R", "--read-cache",action="store_true",default=False,\
    help="read cache file from ini-file (if LSCDataFind is broken)")

parser.add_option("-P", "--priority",action="store",type="int",\
    metavar=" PRIO",help="run jobs with condor priority PRIO")


#parser.add_option("-f", "--config-file",action="store",type="string",\
   # metavar=" FILE",help="use configuration file FILE")

parser.add_option("-p", "--log-path",action="store",type="string",\
    metavar=" PATH",help="directory to write condor log file")

parser.add_option("-o", "--output-segs",action="store_true",default=False,\
    help="output the segment lists of analyzed data")

parser.add_option("-x","--dax", action="store_true", default=False,\
    help="create a dax instead of a dag")


( opts , args ) = parser.parse_args()

##############################################################################
# write command line to a file
print sys.argv[0:]


##############################################################################
# create the config parser object and read in the ini file
cp = ConfigParser.ConfigParser()
cp.read(opts.multi_hipe_config_file)
numslides = cp.get('input','num-slides')

##############################################################################
# get the pad and chunk lengths from the values in the ini file
paddata = int(cp.get('data', 'pad-data'))
n = int(cp.get('data', 'segment-length'))
s = int(cp.get('data', 'number-of-segments'))
r = int(cp.get('data', 'sample-rate'))
o = int(cp.get('inspiral', 'segment-overlap'))
length = ( n * s - ( s - 1 ) * o ) / r
overlap = o / r
minsciseg = length + 2 * paddata

##############################################################################
# Based on the start and end time, generate a list of epochs to
# analyze. An entire hipe dag will be run for each of these epochs.
search_epochs = segments.segmentlist()
istart = opts.start_time
while ( istart < opts.end_time ):
  iend = istart + opts.interval
  if iend > opts.end_time:
    iend = opts.end_time
  search_epochs.append(segments.segment(istart,iend))
  istart += opts.interval
# FIXME:  the writing out of the segments should be done at the end so
# that successfully generated dags, etc can be maintained from run to
# run
segmentsUtils.tosegwizard(file("multi_hipe_selectedsegs.txt",'w'),search_epochs)

##############################################################################
# Read in all the segment lists
ifolist = []
segdict = {}
if opts.h1_segments:
  tmplist = segmentsUtils.fromsegwizard(file(opts.h1_segments)).coalesce()
  segdict["H1"] = cleanlist(tmplist, minsciseg)
  ifolist.append("H1")

if opts.h2_segments:
  tmplist = segmentsUtils.fromsegwizard(file(opts.h2_segments)).coalesce()
  segdict["H2"] = cleanlist(tmplist, minsciseg)
  ifolist.append("H2")

if opts.l1_segments:
  tmplist = segmentsUtils.fromsegwizard(file(opts.l1_segments))
  segdict["L1"] = cleanlist(tmplist, minsciseg)
  ifolist.append("L1")

if opts.v1_segments:
  tmplist = segmentsUtils.fromsegwizard(file(opts.v1_segments))
  segdict["V1"] = cleanlist(tmplist, minsciseg)
  ifolist.append("V1")
################################################################################
# Setting the options for lalapps_inspiral_hipe 
################################################################################
multi_hipe_options = ["verbose", "h1_segments", "h2_segments", "l1_segments", "v1_segments", "start_time", "end_time", "interval", "ninjections", "multi_hipe_config_file", "plot_segments"]


hipe_arguments_tmp=" ".join("--%s %s" % (name.replace("_","-"), value) for name, value in opts.__dict__.items() if value is not None and value is not False and name not in multi_hipe_options)

hipe_arguments=hipe_arguments_tmp + " --config-file config.ini"
print hipe_arguments

##############################################################################
# loop over the intervals, constructing overlapping segment lists,
# making directories, and writing output to them
for interval in search_epochs:

  # name and the directory
  idirectory = str(interval[0])+"-"+str(interval[1])
  os.mkdir(idirectory)

  # extract the segmentlist for each of the interferometers
  for ifo in ifolist:
    tmplist = getSegments(segdict[ifo], interval)

    # now we need to fix segments that overlap the start and end of
    # the interval.  This is where the logic can get ugly, so read
    # with care.  First, handle the segment overlapping the start
    # .......
    try:
      segindex = tmplist.find(interval[0])
      tmpseg = tmplist[segindex]
      if ( tmpseg[1] - interval[0] >= minsciseg ):
        modifiedstart = max( interval[0] - overlap/2 - paddata - 1 , tmpseg[0] )
      else:
        modifiedstart = max( \
            min( tmpseg[1] - minsciseg, interval[0] - overlap/2 - paddata - 1 ),\
            tmpseg[0] )
    except ValueError, e:
      modifiedstart = interval[0]
      if opts.verbose:
        print ifo + ": No segment containing interval start " + str(e)

    # ....... and now the one overlapping the end time .......
    try:
      segindex = tmplist.find(interval[1])
      tmpseg = tmplist[segindex]
      if ( interval[1] - tmpseg[0] >= minsciseg ):
        modifiedend = min( interval[1] + overlap/2 + paddata + 1, tmpseg[1] )
      else:
        modifiedend = min( \
            max( tmpseg[0] + minsciseg, interval[1] + overlap/2 + paddata + 1 ),\
            tmpseg[1] )
    except ValueError, e:
      modifiedend = interval[1]
      if opts.verbose:
        print ifo + ": No segment containing interval end " + str(e)

    modifiedinterval = segments.segmentlist(\
        [segments.segment(modifiedstart,modifiedend)])

    tmplist = tmplist & modifiedinterval

    # write out the segment list to a segwizard file
    tmpoutfile = idirectory+"/"+ifo+"_selectedsegs.txt"
    segmentsUtils.tosegwizard(file(tmpoutfile,'w'),tmplist)

    # plot the segment lists
    if opts.plot_segments:
      pylab.figure()
      pylab.hold(True)
      y = pylab.asarray([0,0])
      y = y + 0.1
      for seg in tmplist:
        pylab.plot(seg,y,'b',linewidth=4)
      y = y + 0.1
      for seg in segdict[ifo]:
        pylab.plot(seg,y,'k',linewidth=4)
      pylab.axvline(interval[0], color='g')
      pylab.axvline(interval[1], color='g')
      pylab.axvline(interval[0]-overlap/2-paddata, color='r')
      pylab.axvline(interval[1]+overlap/2+paddata, color='r')
      pylab.ylim([0.0,0.5])
      pylab.xlim([interval[0]-2*minsciseg,interval[1]+2*minsciseg])
      pylab.savefig(ifo+"-"+idirectory+".png")

  # Next thing is to generate the dag for this interval of time.
  # The steps here are:
  #   1. make dir for zero-lag and playground
  #   2. copy in ini file (and modify if needed)
  #   3. generate dag
  #   4. make dir for injections and run inspinj
  #   5. repeat 2 & 3
  #   6. repeat 4 & 5 as needed

  ############################################################################
  # set up the analysis dag for this interval
  #
  # In doing this, we simply copy the configuration file into the
  # sub-directory and then run the dag generation script.  The exact
  # arguments to that script could be a problem, but we'll deal with
  # that later.  For now, we want to try it.
  #hipe_arguments = " --h1-data --h2-data --l1-data --two-ifo --three-ifo"
  #hipe_arguments += " --output-segs --log-path " + opts.log_path
  #hipe_arguments += " --config-file config.ini --datafind --template-bank"
  #hipe_arguments += " --inspiral --coincidence --trigbank --inspiral-veto"
  #hipe_arguments += " --second-coinc"

  cp.set('input','h1-segments',"../H1_selectedsegs.txt")
  cp.set('input','h2-segments',"../H2_selectedsegs.txt")
  cp.set('input','l1-segments',"../L1_selectedsegs.txt")
  cp.set('input','v1-segments',"../V1_selectedsegs.txt")
  cp.set('input','num-slides',numslides)
  cp.set('input','injection-file','')
  analysisdir = idirectory+"/analysis"
  os.mkdir(analysisdir)
  os.chdir(analysisdir)
  cp.write(file("config.ini",'w'))
  os.system("lalapps_inspiral_hipe " + hipe_arguments)
  os.chdir("../../")

  ############################################################################
  # set up the injection dag for this interval
  for inj in range(1,opts.ninjections + 1):

    # read in the configuration file for injections
    cpinj = ConfigParser.ConfigParser()
    cpinj.read(("injections%02d.ini" % opts.ninjections))

    # make the injection directory
    injectiondir = idirectory+"/injections%02d" % inj
    os.mkdir(injectiondir)
    os.chdir(injectiondir)

    # generate the list of injections
    injectionfile = "HL-INJECTION-" + str(modifiedstart) +\
      "-" + str(abs(modifiedinterval)) + ".xml"
    inspinjcmd = "lalapps_inspinj "
    inspinjcmd += " --gps-start-time " + str(modifiedstart)
    inspinjcmd += " --gps-end-time " + str(modifiedend)
    inspinjcmd += " --output " + injectionfile 
    inspinjcmd += " --seed " + str(random.randint(1,1000))
    for opt in cpinj.options('parameters'):
      arg = string.strip(cpinj.get('parameters',opt))
      inspinjcmd += " --" + opt + " " + arg
    os.system(inspinjcmd)
      
    # update parameters in the config file and write out
    cp.set('input','num-slides','')
    cp.set('input','injection-file',injectionfile)
    cp.write(file("config.ini",'w'))

    # generate the dag
    os.system("lalapps_inspiral_hipe " + hipe_arguments)
    os.chdir("../../")

  # Ultimately, we would like to construct a dag to allow running
  # these sub dags. For now,  we'll do it by hand .....

sys.exit(0)


