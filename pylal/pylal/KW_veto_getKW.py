#!/usr/bin/env python
#
# Copyright (C) 2009  Tomoki Isogai
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

"""
%prog --channel_name=channel_name [--segment_file=File | --gps_start_time=GPSTime --gps_end_time=GPSTime] [options]

Tomoki Isogai (isogait@carleton.edu)

This program gets all the KW triggers that fall into a specified segment list above a certain threshold for the channels specified.
Supported ifos are H1, H2, L1, and V1.
If you are running at the right cluster specified below, the code will use the following path as KW daily dump directory as default.

For S5 LIGO triggers:
/archive/home/lindy/public_html/triggers/s5/ at CIT

For post S5, trigger locations are:
H1:
/archive/home/lindy/public_html/triggers/s6/ at LHO
H2:
/archive/home/lindy/public_html/triggers/s6/ at LLO
V1:
/archive/home/mabizoua/public_html/KW/ at CIT

If you want to use triggers at other location, you can specify --KW_location.

For post S5 triggers, channel name must follow the notation:
(ifo)_(channel name)_(min freq)_(max freq) in capital letters.
For example,
H1_LSC-DARM_ERR_64_1024
V1_Pr_B1_ACp_40_1250

For S5 LIGO triggers, channel name follows the notation:
s5_(ifo)_(channel name) in small letters.
For example,
s5_l1_pobi

You can omit the frequency part if there is no ambiguity, but the code will give you an error if there are several possibilities.
"ls" the above directory to see channel names and max/min frequency available.

If --out_format is given, a file will be created for each channel with the name
(channel)-(start_time)-(duration)_KWtrigs.(specified extention)
If --name_tag is given, (name_tag)_ will be added to the name as prefix.
If --out_format is not given, the code prints out the result in stdout.
You can specify --order_by to sort the output. Supported options are 'GPSTime asc' for ascending time, 'GPSTime desc' for descending time, 'KWSignificance asc' for ascending KW Significance, and 'KWSignificance desc' for descending KW Significance. Default is 'GPSTime asc'.
"""

# =============================================================================
#
#                               PREAMBLE
#
# =============================================================================


from __future__ import division

import os
import sys
import shutil
import re
import optparse

try:
    import sqlite3
except ImportError:
   # pre 2.5.x
   from pysqlite2 import dbapi2 as sqlite3

from glue.segments import segment, segmentlist
from glue import segmentsUtils

from pylal import git_version
from pylal import KW_veto_utils as utils

__author__ = "Tomoki Isogai <isogait@carleton.edu>"
__date__ = "7/10/2009"
__version__ = "2.0"

def parse_commandline():
    """
    Parse the options given on the command-line.
    """
    parser = optparse.OptionParser(usage = __doc__,version=git_version.verbose_msg)

    parser.add_option("-K", "--KW_location", default=None,
                      help="Location of KW trigger folder if you are not using the folder specified in --help.")
    parser.add_option("-c", "--channel_name", action="append", default=[],
                      help="Channel names you want KW triggers from. See --help for the naming format. Can be provided multiple times to specify more than one channel. At least one channel is required.")
    parser.add_option("-m", "--min_thresh", type="int", default=0,
                      help="Minimum KW significance threshold for KW triggers. (Default: 0)")
    parser.add_option("-S", "--segment_file", default=None,
                      help="Segment file on which KW triggers are retrieved. This option or --gps_start_time and --gps_end_time are required.")
    parser.add_option("-s", "--gps_start_time", type="int",
                      help="GPS start time on which the KW triggers are retrieved. Required unless --segment_file is specified.")
    parser.add_option("-e", "--gps_end_time", type="int",
                      help="GPS end time on which KW triggers are retrieved. Required unless --segment_file is specified.")
    parser.add_option("-n", "--name_tag", default=None,
                      help="If given, this will be added as prefix in the output file name.")
    parser.add_option("-f", "--out_format", default='stdout',
                      help="Save in this format if specified, otherwise output in stdout. See --help for the supported format and output file name.")
    parser.add_option("-o","--order_by",default="GPSTime asc",
                      help="Order of the output. See -help for supported order. (Default: GPSTime asc)")
    parser.add_option("-l", "--scratch_dir", default=".",
                      help="Scratch directory to be used for database engine. Specify local scratch directory for better performance and less fileserver load.")
    parser.add_option("-v", "--verbose", action="store_true", default=False,
                      help="Run verbosely")
    
    opts, args = parser.parse_args()
    
    ############################ Sanity Checks ################################

    # at least 1 channel is necessary
    if len(opts.channel_name) < 1:
      parser.error("Error: at least 1 channel for --channel_name is required.")

    # check if necessary input exists
    if opts.segment_file is not None:    
      if not os.path.isfile(opts.segment_file):
        parser.error("Error: segment file %s not found"%opts.segment_file)
    else:
      for t in ('gps_start_time','gps_end_time'):
        if getattr(opts,t) is None:
          parser.error("Error: either --segment_file or --gps_start/end_time must be specified.")
    
    if opts.channel_name == []: #
        parser.error("Error: --channel_name is a required parameter")
    
    # check if scratch directory exists
    if not os.path.exists(opts.scratch_dir):
        parser.error("Error: %s does not exist"%opts.scratch_dir)
    
    # standardize the format (starting from period .) and see it's supported    
    if opts.out_format[0] != ".":
      opts.out_format = "." + opts.out_format
      if opts.out_format not in (\
         ".stdout",".pickle",".pickle.gz",".mat",".txt",".txt.gz",".db"):
        parser.error("Error: %s is not a supported format"%opts.out_format)
       
    # check order_by is valid and avoid SQL injection attack
    if opts.order_by not in ("GPSTime asc","GPSTime desc","KWSignificance asc","KWSignificance desc"):
      parser.error("Error: %s is not valid. See -help for the supported order."%opts.order_by)
 
    # show parameters
    if opts.verbose:
      print >> sys.stderr, "running KW_veto_getKW.py..."
      print >> sys.stderr, git_version.verbose_msg
      print >> sys.stderr, ""
      print >> sys.stderr, "******************* PARAMETERS ********************"
      print >> sys.stderr, "KW trigger directory:"
      if opts.KW_location is None:
        print >> sys.stderr, "/archive/home/lindy/putlic_html/triggers/s5 for S5 LIGO triggers"
        print >> sys.stderr, "/archive/home/lindy/public_html/triggers/s6 for post S5 LIGO triggers"
        print >> sys.stderr, "/archive/home/mabizoua/public_html/KW for Virgo triggers"
      else:
        print >> sys.stderr, opts.KW_location;
      print >> sys.stderr,'channels:'
      print >> sys.stderr, opts.channel_name; 
      print >> sys.stderr,'minimum threshold:'
      print >> sys.stderr, opts.min_thresh;
      if opts.segment_file is not None:
        print >> sys.stderr, 'segment file:'
        print >> sys.stderr, opts.segment_file; 
      else:
        print >> sys.stderr, 'start/end GPS time:'
        print >> sys.stderr, "%d - %d"%(opts.gps_start_time,opts.gps_end_time)
      print >> sys.stderr, 'name tag:'
      print >> sys.stderr, opts.name_tag;
      print >> sys.stderr,'output format:'
      print >> sys.stderr, opts.out_format[1:];
      print >> sys.stderr, 'order by:'
      print >> sys.stderr, opts.order_by;
      print >> sys.stderr,'scratch directory:'
      print >> sys.stderr, opts.scratch_dir;
      print >> sys.stderr, ""
        
    return opts

def find_file_from_channel(channel,daily_dir):
    """
    From the files in daily_dir, find the file that contains KW triggers
    for the channel, and return the full channels name and the file path.
    """
    # this function standardize the channel names
    def norm(channel_name):
      return channel_name.upper().replace("-","_")

    all_files = \
         [f for f in os.listdir(daily_dir) if os.path.splitext(f)[1] == '.trg']
    candidates = [f for f in all_files if re.match(norm(channel), norm(f)) != None]
    if len(candidates) == 1:
      full_name = norm(os.path.splitext(candidates[0])[0])
      return os.path.join(daily_dir,candidates[0]), full_name
    elif len(candidates) < 1:
      print >> sys.stderr, "Warning: no match found for %s in %s. See --help for channel name format. Attempting to ignore..."%(channel, daily_dir)
      return None, None
    # When there are more than two possibilities see if one of their name
    # mathes perfectly. If so, warn user and use it, otherwise give an error
    # and show the possibilities.
    else:
      refined_candidates = [f for f in candidates if re.match(norm(os.path.splitext(f)[0]),norm(channel)) != None]
      if len(refined_candidates) == 1:
        print >> sys.stderr, """
        Warning: Multiple possible files with the channel name %s:
                 %s
                 Using %s and ignoring...
        """%(channel,", ".join(candidates),refined_candidates[0])
        full_name = norm(os.path.splitext(refined_candidates[0])[0])
        return os.path.join(daily_dir,refined_candidates[0]), full_name
      else:
        print >> sys.stderr, """
        Error: Multiple possible files with the channel name %s:
               %s
               Please change the channels name so that it is unique.
               See --help for channel name format.
        """%(channel, ", ".join(candidates))
        sys.exit(1)

def get_trigs(channel, segs, min_thresh, trigs_loc=None,name_tag=None,\
              scratch_dir=".",verbose=True):
    """
    Get time and KW significance of KW triggers for a particular channel that
    occured in the specified segments and above specified KW significance
    threshold. 
    ifo has to be one of H1, H2, L1, V1.
    """
    if verbose: print >> sys.stderr, "getting data for %s..."%channel
  
    ## initialize SQLite database
    start_time = segs[0][0]
    end_time = segs[-1][1]
    duration = end_time - start_time
    prefix = "%s-%d-%d-KWtrigs"%(channel, start_time, duration)
    if name_tag != None:
      prefix = name_tag + "_" + prefix
    dbname = prefix+".db"
    # if the name already exists, rename the old one to avoid collision
    utils.rename(dbname)

    global KW_working_filename # so that it can be erased when error occurs
    KW_working_filename = utils.get_connection_filename(\
                         dbname,tmp_path=scratch_dir,verbose=verbose)

    KWconnection = sqlite3.connect(KW_working_filename)
    KWcursor = KWconnection.cursor()

    ## create a table for retrieved triggers
    KWcursor.execute('create table KWtrigs (GPSTime double, KWSignificance double, frequency int)')

    ## determine the KW trigger file we need
    ifo = channel.split("_")[0].upper()

    # Virgo case
    if trigs_loc is None and ifo == "V1":
      trigs_loc = "/archive/home/mabizoua/public_html/KW/"
    # LIGO case
    elif trigs_loc is None and ifo in ("H0","H1","H2","L0","L1"):
        trigs_loc = "/archive/home/lindy/public_html/triggers/s6/"
    # for S5, first two letters were not ifo but 's5'
    elif trigs_loc is None and ifo == "S5":   
      trigs_loc = "/archive/home/lindy/public_html/triggers/s5/"
    elif trigs_loc is None:
      print >> sys.stderr, "Error: channel name %s is unsupported. See --help for name format."%channel
      sys.exit(1)

    # sanity check
    if not os.path.exists(trigs_loc):
      print >> sys.stderr, "Error: KW daily dump %s not found."%trigs_loc

    ## select daily dump directories in the folder
    daily_dump_dirs = \
            [d for d in os.listdir(trigs_loc) if \
             re.match(r"(?P<start>\d{9,10}?)_(?P<end>\d{9,10}?)",d) != None \
             and len(d) < 22]
    # sort by time
    daily_dump_dirs.sort()

    # get the necessary daily dump folder
    analyze_dumps = [d for d in daily_dump_dirs \
                    if int(d.split("_")[0]) < end_time and \
                       int(d.split("_")[1]) > start_time]

    # for S5, some channels are not in certain dumps: do a little trick to keep
    # the name
    full_channel_name = ""
    for dump in analyze_dumps:
      trigs_file, tmp_name = \
        find_file_from_channel(channel,os.path.join(trigs_loc,dump))
      if full_channel_name == "" and tmp_name != None:
        full_channel_name = tmp_name
      if trigs_file != None:
        if verbose: print "retreiving data from %s..."%trigs_file

        for line in  open(trigs_file):
          # get central time and KW significance
          trig = line.split()
          t = float(trig[2]) 
          s = float(trig[7])
          f = float(trig[3])

          # check if KW trig is in the given segment and if its significance 
          # is above the minimum specified
          # FIXME: check inf/nan just in case
          if t in segs and s > min_thresh:
            # insert into the database
            # micro second for GPS time is accurate enough
            KWcursor.execute("insert into KWtrigs values (?, ?, ?)", ("%.3f"%t, "%.2f"%s,"%d"%f))

    if full_channel_name == "": # means there is no KW trigger
      full_channel_name = channel # better than nothing...

    # commit the insertions to the database
    KWconnection.commit()

    return dbname, KW_working_filename, KWconnection, KWcursor, full_channel_name

    
# =============================================================================
#
#                                    Main
#
# =============================================================================

def main():
    ## parse the command line
    opts = parse_commandline()
    
    ## get the segments on which triggers are retrieved
    # case 1: segment file is given   
    if opts.segment_file is not None:
      # read segment file
      if opts.segment_file.endswith(".txt"):
        seg_list = utils.read_segfile(opts.segment_file)
      elif opts.segment_file.endswith(".xml") or opts.segment_file.endswith(".xml.gz"):
        seg_list = utils.read_segfile_xml(opts.segment_file,opts.verbose)
    # case 2: start and end GPS time are given
    else:
      seg_list = segmentlist([segment(opts.gps_start_time,opts.gps_end_time)])
   
    ## loop over each channels and get KW triggers
    for chan in opts.channel_name:      
      # wrap in try/finally clause so that the code can erase temporary database
      # when it encounters an error and had to stop in the middle
      try:
        ## get triggers
        dbname, KW_working_filename, KWconnection, KWcursor, full_chan_name = \
            get_trigs(chan,seg_list,opts.min_thresh,trigs_loc=opts.KW_location,name_tag=opts.name_tag,scratch_dir=opts.scratch_dir,verbose=opts.verbose)
    
        # save/display the result
        outname = dbname.replace(chan, full_chan_name).replace(".db", opts.out_format)
        utils.save_db(KWcursor, "KWtrigs", outname, KW_working_filename,
                      order_by=opts.order_by, verbose=opts.verbose)
     
        # close the connection to database
        KWconnection.close()
           
      finally:
        # erase temporal database
        if globals().has_key('KW_working_filename'):
          db = globals()['KW_working_filename']
          if opts.verbose:
            print >> sys.stderr, "removing temporary workspace '%s'..." % db
          os.remove(db)

if __name__ == "__main__":
    main()    
