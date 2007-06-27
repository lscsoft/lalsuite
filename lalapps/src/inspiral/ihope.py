#!/usr/bin/env @PYTHONPROG@
import os, sys, copy
import ConfigParser
import optparse
import tempfile
from glue import segments
from glue import segmentsUtils
from glue import pipeline


##############################################################################
def make_external_call(command, show_stdout=False, show_command=False):
    """
    Run a program on the shell and print informative messages on failure.
    """
    if show_command: print command
    
    stdin, out, err = os.popen3(command)
    pid, status = os.wait()

    if status != 0:
        print >>sys.stderr, "External call failed."
        print >>sys.stderr, "  status: %d" % status
        print >>sys.stderr, "  stdout: %s" % out.read()
        print >>sys.stderr, "  stderr: %s" % err.read()
        print >>sys.stderr, "  command: %s" % command
        sys.exit(status)
    if show_stdout:
        print out.read()
    stdin.close()
    out.close()
    err.close()

##############################################################################
#
#  MAIN PROGRAM
#
##############################################################################
usage = """usage: %prog [options] 
"""

parser = optparse.OptionParser( usage )

# arguments
parser.add_option("-f", "--config-file",action="store",type="string",\
    metavar=" FILE",help="use configuration file FILE")

parser.add_option("-p", "--log-path",action="store",type="string",\
    metavar=" PATH",help="directory to write condor log file")

parser.add_option("-s", "--gps-start-time",action="store",type="int",\
    metavar=" GPS_START",help="begin analysis at GPS_START")

parser.add_option("-e", "--gps-end-time",action="store",type="int",\
    metavar=" GPS_END",help="end analysis at GPS_END")

(opts,args) = parser.parse_args()

##############################################################################
# Sanity check of input arguments
if not opts.config_file:
  print >> sys.stderr, "No configuration file specified."
  print >> sys.stderr, "Use --config-file FILE to specify location."
  sys.exit(1)

if not opts.log_path:
  print >> sys.stderr, "No log file path specified."
  print >> sys.stderr, "Use --log-path PATH to specify a location."
  sys.exit(1)

##############################################################################
# parse the ini file:
cp = ConfigParser.ConfigParser()
cp.read(opts.config_file)

# set up the analysis directory
analysisDirectory = str(opts.gps_start_time) + "-" + str(opts.gps_end_time)
os.mkdir(analysisDirectory)
os.chdir(analysisDirectory)

# work out the hipe call:
hipeCommand = cp.get("condor","hipe")
hipeCommand += " --log-path " + opts.log_path
for item in cp.items("hipe-arguments"):
  hipeCommand += " --" + item[0] + " " + item[1]

# create the hipe ini file
hipecp = copy.deepcopy(cp)
hipecp.remove_section("injections")
hipecp.remove_section("hipe-arguments")
hipecp.remove_option("condor","hipe")
hipecp.remove_section("segments")

##############################################################################
# Set up the IFOs and get the appropriate segments

ligoIfos = ["H1","H2","L1"]
ifos = []
for option in ["g1-data","h1-data","h2-data","l1-data"]:
  if cp.has_option("hipe-arguments",option): ifos.append(option[0:2].upper() )
if cp.has_option("hipe-arguments","analyze-all"): 
  ifos = ["G1", "H1", "H2", "L1"]

print "Setting up an analysis for " + str(ifos) + " from " + \
    str(opts.gps_start_time) + " to "  + str(opts.gps_end_time)
print

# Run LSCsegFind and LSCdataFind to determine the segments to analyze
for ifo in ifos:
  if ifo in ligoIfos: type = cp.get("input","ligo-type")
  elif ifo == "G1": type =   cp.get("input","geo-type")
 
  ifo_type = ifo + "_" + type
  dataFindFile = ifo_type + "-" + str(opts.gps_start_time) + "-" + \
      str(opts.gps_end_time - opts.gps_start_time) + ".txt"
  segFindFile = ifo + "_SCI_SEGS-" + str(opts.gps_start_time) + "-" + \
      str(opts.gps_end_time - opts.gps_start_time) + ".txt"
  segFile = ifo + "_selectedsegs.txt"
  missedFile = ifo + "_missedsegs.txt"

  print "Running LSCsegFind to determine science segments for " + ifo
  segFindCall = "LSCsegFind --interferometer=" + ifo + \
      " --type=" + cp.get("segments","analyze") + \
      " --gps-start-time=" + str(opts.gps_start_time) + \
      " --gps-end-time=" + str(opts.gps_end_time) + " > " + \
      segFindFile
  make_external_call(segFindCall)
  sfSegs = segmentsUtils.fromsegwizard(file(segFindFile)).coalesce()

  print "Running LSCdataFind to determine available data from " + type + \
      " frames for " + ifo
  dataFindCall = "LSCdataFind --observatory=" + ifo[0] + \
      " --type=" + ifo_type + \
      " --gps-start-time=" + str(opts.gps_start_time) + \
      " --gps-end-time=" + str(opts.gps_end_time) + " --show-times > " + \
      dataFindFile
  make_external_call(dataFindCall)
  dfSegs = segmentsUtils.fromsegwizard(file(dataFindFile)).coalesce()

  analyzedSegs = sfSegs.__and__(dfSegs)
  analyzedSegs = dfSegs
  segmentsUtils.tosegwizard(file(segFile,"w"), analyzedSegs)
  hipecp.set("input",ifo.lower() + "-segments","../" + segFile)
  print "Writing " + ifo + " segments of total time " + \
      str(analyzedSegs.__abs__()) + "s to file: " + \
      segFile

  missedSegs = sfSegs.__and__(dfSegs.__invert__()) 
  segmentsUtils.tosegwizard(file(missedFile,"w"), missedSegs)
  print "Writing " + ifo + " segments which cannot be analyzed to file " + \
      missedFile
  print "Not analyzing %d s, representing %.2f percent of time" %  \
     (missedSegs.__abs__(), 
     100. * missedSegs.__abs__() / analyzedSegs.__abs__() )
  print

##############################################################################
# create a log file that the Condor jobs will write to
basename = opts.config_file.rstrip(".ini")
logname = basename + '.dag.log.'
tempfile.tempdir = opts.log_path
logfile = tempfile.mktemp(prefix=logname)
fh = open( logfile, "w" )
fh.close()

##############################################################################
# create the DAG writing the log to the specified directory
dag = pipeline.CondorDAG(logfile)
dag.set_dag_file(basename)



##############################################################################
# Set up the directories for each run and run lalapps_inspiral_hipe
hipeRuns = cp.items("injections")
hipeRuns.append( ("slide-zero", "") )
     
print "Running inspiral hipe for analysis and injection runs"
for (hipeDir, injFile) in hipeRuns:

  os.mkdir(hipeDir)
  os.chdir("..")

  print "Running hipe in directory " + hipeDir 
  if injFile: 
    print "Injection file: " + injFile 
  else: print "No injections, " + str(cp.get("input","num-slides")) + \
      " time slides"
  print

  # link the executables 
  for (job, executable) in cp.items("condor"):
    if job != "universe":
      if executable[0] != "/": executable = "../../" + executable
      os.symlink(executable, 
          analysisDirectory + "/" + hipeDir + "/" + executable.split("/")[-1])
      cp.set("condor", job, executable.split("/")[-1] )

  # link the injection file
  if injFile:
    if injFile[0] != "/": injFile = "../../" + injFile
    os.symlink(injFile, 
        analysisDirectory + "/" + hipeDir + "/" + injFile.split("/")[-1])
    cp.set("condor", job, injFile.split("/")[-1] )
    hipecp.set("input", "num-slides", "")

  else: 
    hipecp.set("input","num-slides", cp.get("input","num-slides") )

  os.chdir(analysisDirectory + "/" + hipeDir)
  iniFile = "inspiral_hipe." + hipeDir + ".ini"
  hipecp.write(file(iniFile,"w"))

  make_external_call(hipeCommand + " --config-file " + iniFile) 

  hipeJob = pipeline.CondorDAGManJob(
      hipeDir + "/" + iniFile.rstrip("ini") + "dag")
  hipeNode = pipeline.CondorDAGNode(hipeJob)
  hipeNode.set_post_script( cp.get("condor", "rescue-script") )

  hipeNode.add_post_script_arg( 
      hipeNode.job().get_sub_file().rstrip(".condor.sub") )
  dag.add_node(hipeNode)

  os.chdir("..")

dag.write_sub_files()
dag.write_dag()

print 
print "  Created a DAG file which can be submitted by executing"
print "\n    condor_submit_dag " + analysisDirectory + "/" + dag.get_dag_file()
print "\n  from a condor submit machine"
print "\n  Before submitting the dag, you must execute"
print "\n    export _CONDOR_DAGMAN_LOG_ON_NFS_IS_ERROR"
print """
  If you are running LSCdataFind jobs, do not forget to initialize your grid
  proxy certificate on the condor submit machine by running the commands
  
    unset X509_USER_PROXY
    grid-proxy-init -hours 72
  
  Enter your pass phrase when promted. The proxy will be valid for 72 hours.
  If you expect the LSCdataFind jobs to take longer to complete, increase the
  time specified in the -hours option to grid-proxy-init. You can check that
  the grid proxy has been sucessfully created by executing the command:
  
    grid-cert-info -all -file /tmp/x509up_u`id -u`
  
  This will also give the expiry time of the proxy."""
