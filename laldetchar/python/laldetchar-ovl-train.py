usage = """ovl-train.py [--options] gps_start gps_stop"""
description = """written to generate an ovl vetolist over the specified range"""
__author__= 'Reed Essick (reed.essick@ligo.org)'

#=================================================

import os

from laldetchar.idq import idq
from laldetchar.idq import ovl

from glue.ligolw import ligolw
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import lsctables

import ConfigParser
import subprocess

from optparse import OptionParser

#=================================================

parser=OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-c", "--config", default="./config.ini", type="string")

parser.add_option("", "--ignore-science-segments", default=False, action="store_true")

parser.add_option("", "--skip-collect-sngl-chan", default=False, action="store_true", help="if you already have the correct single channel summary files, skip that step here")

parser.add_option("", "--redundancies", default=False, action="store_true")
parser.add_option("", "--safety", default=False, action="store_true")

opts, args = parser.parse_args()

if len(args) != 2:
	raise ValueError("please supply exactly 2 arguments")

gpsstart = int(args[0])
gpsstop = int(args[1])

### check config
if not os.path.exists(opts.config):
	raise ValueError("--config=%s does not exist"%opts.config)
config = ConfigParser.SafeConfigParser()
config.read(opts.config)

#=================================================

### make output directory
output_dir = "%s/%d_%d/"%(config.get("general","traindir"), gpsstart, gpsstop)
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

#========================
### get science segments
if not opts.ignore_science_segments:
	### query the database
	args = '%s -t %s -q -a %s -s %d -e %d' % (config.get("get_science_segments", "program"), 
					config.get("get_science_segments", "segdb"), 
					config.get("get_science_segments", "include"), 
					gpsstart, 
					gpsstop)

	if opts.verbose:
		print "getting science segments"
		print "\t", args

	segfile = "%s/science_segments-%d-%d.xml"%(output_dir, gpsstart, int(gpsstop-gpsstart))
	segfile_obj = open(segfile, "w")
	p = subprocess.Popen(args.split(), stdout=segfile_obj)
	p.wait()
	segfile_obj.close()

	### read in segments from xml file
	(scisegs, coveredseg) = idq.extract_dq_segments(segfile, config.get('get_science_segments', 'include'))
else:
	scisegs = [ [gpsstart, gpsstop] ]

### write segments to ascii list
if opts.verbose:
	print "writing segments to ascii list"

sciseg_path = "%s/science_segments-%d-%d.seg"%(output_dir, int(gpsstart), int(gpsstop-gpsstart))
f = open(sciseg_path, 'w')
for line in scisegs:
	print >> f, line[0], line[1]
f.close()

#========================
### create single channel summary files
snglchndir = config.get("general","snglchndir")
if opts.skip_collect_sngl_chan:
	if opts.verbose:
		print "skipping single channel summary file generation. Expecting all data to already exist in : %s"%snglchndir
else:
	if opts.verbose:
		print "writing single channel summary files into : %s"%snglchndir

	if not os.path.exists(snglchndir):
		os.makedirs(snglchndir)

	new_dirs = []
	for kwconfig, kwtrgdir in eval(config.get("general","kw")).items(): ### iterate over all kw directories
		if opts.verbose:
			print "\tkwconfig : %s\n\tkwtrgdir : %s"%(kwconfig, kwtrgdir)
		new_dirs += idq.collect_sngl_chan_kw(gpsstart, gpsstop, kwconfig, source_dir=kwtrgdir, output_dir=snglchndir)

	new_dirs = set(new_dirs)
	if opts.verbose:
		print "created %d new directories :"%len(new_dirs)
		for new_dir in new_dirs:
			print "\t", new_dir

#=================================================
### launch training job

#========================
### build params object
if opts.verbose:
	print "generating params object"

analysis_range = [gpsstart, gpsstop]

### specify ovl pointer
ovlsegs = [sciseg_path]

### load from config object
auxdir = snglchndir
gwdir = snglchndir

gwchans = eval(config.get('general', 'gwchannels'))
gwthr = config.getfloat('general', 'gw_kwsignif_thr')

ifos = [config.get('general', 'ifo')]

gwsets = eval(config.get('ovl_train', 'gwsets'))

safety = config.get('ovl_train', 'safety')

windows = eval(config.get('ovl_train', 'windows'))
thresholds = eval(config.get('ovl_train', 'thresholds'))

Psigthr = config.getfloat('ovl_train', 'Psigthr')
effbydtthr = config.getfloat('ovl_train', 'effbydtthr')

if config.has_option('general', 'selected-channels'):
	channels = config.get('general', 'selected-channels')
else:
	channels = ''

if config.has_option('general', 'unsafe-channels'):
	notused = config.get('general', 'unsafe-channels')
else:
	notused = ''

metric = config.get('ovl_train', 'metric')

### parse channels and not used
if channels == '':
	channels = False
if notused != '':
	notused = [l.strip('\n') for l in open(notused, 'r').readlines() if l.strip("\n")]
else:
	notused = []

params = ovl.params(analysis_range, 
		auxdir, 
		gwdir, 
		gwchans, 
		gwthr, 
		ifos, 
		gwsets, 
		metric=metric,
		scisegs=[sciseg_path], 
		vetosegs=None,
		channels=channels, 
		notused=notused, 
		windows=windows, 
		thresholds=thresholds, 
		Psigthr=Psigthr, 
		effbydtthr=effbydtthr, 
		safety=safety)

if opts.safety: ### a safety study
	if opts.verbose:
		print "launching safety study"

	vetolists = ovl.safety(params, output_dir=output_dir, verbose=opts.verbose, write_channels=True)

elif not config.has_option('ovl_train', 'convergent'): # "normal" training
	num_runs = config.getint('ovl_train', 'num_runs')
	incremental = config.getint('ovl_train', 'incremental')

	if opts.verbose:
		print "launching \"normal\" training with:\n\t%d runs\n\t%d incremental"%(num_runs, incremental)

	### launch training job
	vetolists = ovl.train(params, num_runs=num_runs, incremental=incremental, output_dir=output_dir, verbose=opts.verbose, write_channels=True)

else: # "convergent" training
	if opts.verbose:
		print "launching \"convergent\" training"

	vetolists = ovl.convergent_train(params, output_dir=output_dir, verbose=opts.verbose, write_channels=True)

### compute redundancies between channels 
if opts.redundancies:
	if opts.verbose: 
		print "launching redundancies job"
	ovl.redundancies(params, output_dir=output_dir, verbose=opts.verbose, write_channels=True)

if opts.verbose:
	print "Done"
