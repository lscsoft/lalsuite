usage = "laldetchar-ovl-evaluate.py [--options] vetolist patfile"
description = "generate predictions from OVL based on the vetolist and patfile. This is provided to allow for apples-to-apples comparison with MLA implementations"
author = "Reed Essick (reed.essick@ligo.org)"

#-------------------------------------------------

import os

import numpy as np

from laldetchar.idq import ovl
from laldetchar.idq import idq
from laldetchar.idq import event

from ConfigParser import SafeConfigParser

from optparse import OptionParser

#-------------------------------------------------

parser = OptionParser(usage=usage, description=description)

parser.add_option('-v', '--verbose', default=False, action='store_true')

parser.add_option('-c', '--config', default=None, type='string',
    help='the config file used with laldetchar-ovl-train')

parser.add_option('', '--output-filename', default=None, type='string',
    help='the name of the output file. If not supplied, a name will be automatically generated \
based on the gps times in patfiles.')

opts, args = parser.parse_args()

assert len(args)>=2, "please supply at least 2 input arguments\n%s"%usage
vetolist = args[0]
patfiles = args[1:]

if not opts.config:
    opts.config = raw_input('--config=')
config = SafeConfigParser()
config.read(opts.config)

#-------------------------------------------------

### find maximum vwin used in vetolist
if opts.verbose:
    print( "finding maximum veto window used in vetolist" )
win=0.001 ### use this as an absolute mimimum. Will almost certainly be replaced, unless vetolist is empty
file_obj = open(vetolist, 'r')
for line in file_obj:
    if line[0]=="#":
        pass
    else:
        win = max(win, float(line.strip().split()[ovl.vD['vwin']]))
file_obj.close()

#---

### find time ranges of interest from patfile
if opts.verbose:
    print( "defining segments in which we need KW triggers" )
gps = idq.slim_load_datfiles(patfiles, columns=['GPS_s', 'GPS_ms'])
gps = np.array(gps['GPS_s'], dtype=float) + 1e-6*np.array(gps['GPS_ms'], dtype=float)

Ngps = len(gps)
if Ngps==0:
    raise ValueError, 'please supply at least one GPS time within : '+patfile

elif opts.verbose:
    print( "found %d times"%Ngps )

segs = event.fixsegments([[t-win, t+win] for t in gps]) ### the segments in which we need KW triggers

#---

### look up KW trg files that intersect segs
if opts.verbose:
    print( "finding relevant kw_trgfiles" )
kw_trgfiles = []
### iterate over different configurations used in training
for kwconf, dirname in eval(config.get('general', 'kw')).items(): ### this is kinda ugly...
    if opts.verbose:
        print( "  searching for KW trgfiles corresponding to %s in %s within [%.3f, %.3f]"%(kwconf, dirname, segs[0][0], segs[-1][1]) )

    ### iterate over all trg files found in that directory
    for trgfile in idq.get_all_files_in_range(dirname, segs[0][0], segs[-1][1], pad=0, suffix='.trg'):
        ### check whether there is some overlap 
        ### not gauranteed if there are gaps between min and max gps times
        if event.livetime(event.andsegments([[idq.extract_start_stop(trgfile, suffix='.trg')], segs])): 
            if opts.verbose:
                print( "    kept : "+trgfile )
            kw_trgfiles.append( trgfile )

        elif opts.verbose:
            print( "    discarded : "+trgfile )

#---

if opts.verbose:
    print( "evaluating %d times using %d KW trgfiles"%(Ngps, len(kw_trgfiles) ) )
### set up output pointers
if opts.output_filename:
    datfile = os.path.basename(opts.output_filename)
    output_dir = os.path.dirname(opts.output_filename)

else:
    datfile = os.path.basename(patfile).replace(".pat", ".ovldat")
    output_dir = os.path.dirname(patfile)

if not output_dir: ### an empty string, meaning it is in the current directory
    output_dir = '.'

### actually run the evaluation
idq.ovl_evaluate(
    vetolist, 
    patfiles=patfiles, 
    kw_trgfiles=kw_trgfiles, 
    filename=datfile, 
    output_dir=output_dir,
)

if opts.verbose:
    print( "predictions written to : %s/%s"%(output_dir, datfile) )
