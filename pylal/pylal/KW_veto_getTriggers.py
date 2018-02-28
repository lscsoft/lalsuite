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
%prog --trigger_files=Files [--segment_file=File | --gps_start_time=GPSTime --gps_end_time=GPSTime] [options]

Tomoki Isogai (isogait@carleton.edu)

This program gets GPS time and SNR of GW triggers that fall into a specified segments above a given SNR threshold.
Input file must be ligolw xml or ascii format. For ligolw xml files, the code deals with sngl_inspiral, sngl_ringdown or singl_burst table. For ascii files, it is assumed that the first column represents GPS time of triggers and the second column represents SNR. The code ignores any empty lines and lines starting with "#" or "%", and gives an error if more that two columns are found.

If --out_format is given, a file will be created for each channel with the name
(start_time)_(end_time)_(SNR threshold)_GWtrigs.(specified extention)
If --name_tag is given, (name_tag)_ will be added to the name as prefix.
If --out_format is not given, the code prints out the result in stdout.
You can specify --order_by to sort the output. Supported options are 'GPSTime asc' for ascending time, 'GPSTime desc' for descending time, 'SNR asc' for ascending SNR, and 'SNR desc' for descending SNR. Default is 'GPSTime asc'.

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
import optparse
import cPickle

try:
    import sqlite3
except ImportError:
   # pre 2.5.x
   from pysqlite2 import dbapi2 as sqlite3

from glue import segmentsUtils
from glue.segments import segment, segmentlist

from pylal import git_version
from pylal import KW_veto_utils

__author__ = "Tomoki Isogai <isogait@carleton.edu>"
__date__ = "7/10/2009"
__version__ = "2.0"

def parse_commandline():
    """
    Parse the options given on the command-line.
    """
    parser = optparse.OptionParser(usage=__doc__,version=git_version.verbose_msg)

    parser.add_option("-t", "--trigger_files", action="append", default=[],
                      help="File containing triggers. Required.")
    parser.add_option("-S", "--segment_file", default=None,
                      help="Segments file by which triggers are filtered. This option or --gps_start_time and --gps_end_time are required.")
    parser.add_option("-s", "--gps_start_time", type="int",
                      help="GPS start time on which triggers are retrieved. Required unless --segment_file is specified.")
    parser.add_option("-e", "--gps_end_time", type="int",
                      help="GPS end time on which triggers are retrieved. Required unless --segment_file is specified.")
    parser.add_option("-m", "--min_thresh", default=8, type="float",
                      help="Code filters triggers below this minimum SNR value. (Default: 8)")
    parser.add_option("-n", "--name_tag", default=None,
                      help="If given, this will be added as prefix in the output file name.")
    parser.add_option("-f", "--out_format", default='.stdout',
                      help="Save in this format if specified, otherwise output in stdout.")
    parser.add_option("-o", "--order_by", default="GPSTime asc",
                      help="Order of the output. See --help for the supported format and output file name. (Default: GPSTime asc)")
    parser.add_option('-l', "--scratch_dir", default='.',
                      help="Scratch directory to be used for database engine. Specify local scratch directory for better performance and less fileserver load. (Default: current directory)")
    parser.add_option("-v", "--verbose", action="store_true",default=False, 
                      help="Run verbosely. (Default: False)")
    
    opts, args = parser.parse_args()
    
    
    ############################## Sanity Checks ##############################

    ## check if necessary input exists
    # trigger file
    if opts.trigger_files == []:
        parser.error("Error: --trigger_files is a required parameter")
    for t in opts.trigger_files:
      if not os.path.isfile(t):
        parser.error("Error: --trigger_files %s not found"%t)

    # segments
    if opts.segment_file is not None:
      if not os.path.isfile(opts.segment_file):
        parser.error("Error: segment file %s not found"%opts.segment_file)
    else:
      for t in ('gps_start_time','gps_end_time'):
        if getattr(opts,t) is None:
          parser.error("Error: either --segment_file or --gps_start/end_time must be specified.")
        
    # check if scratch directory exists
    if not os.path.exists(opts.scratch_dir):
        parser.error("Error: %s does not exist"%opts.scratch_dir)

    # standardize the format (starting from period .) and see if it's supported 
    if opts.out_format[0] != ".":
      opts.out_format = "." + opts.out_format
      if opts.out_format not in (\
         ".stdout",".pickle",".pickle.gz",".mat",".txt",".txt.gz",".db"):
        parser.error("Error: %s is not a supported format"%opts.out_format)

    # check order_by is valid and avoid SQL injection attack
    if opts.order_by not in ("GPSTime asc","GPSTime desc","SNR asc","SNR desc"):
      parser.error("Error: %s is not valid. See -help for the supported order."%opes.order_by)
        
    ## show parameters
    if opts.verbose:
        print >> sys.stderr, "running KW_veto_getTriggers..."
        print >> sys.stderr, git_version.verbose_msg
        print >> sys.stderr, ""
        print >> sys.stderr, "******************** PARAMETERS *****************"
        print >> sys.stderr, 'trigger file:'
        print >> sys.stderr, opts.trigger_files; 
        if opts.segment_file is not None:
          print >> sys.stderr, 'segment file:'
          print >> sys.stderr, opts.segment_file; 
        else:
          print >> sys.stderr, 'start/end GPS time:'
          print >> sys.stderr, "%d - %d"%(opts.gps_start_time,opts.gps_end_time)
        print >> sys.stderr, 'minimum SNR:'
        print >> sys.stderr, opts.min_thresh;
        print >> sys.stderr, 'name tag:'
        print >> sys.stderr, opts.name_tag; 
        print >> sys.stderr, 'output format:'
        print >> sys.stderr, opts.out_format[1:];
        print >> sys.stderr, 'order by:'
        print >> sys.stderr, opts.order_by;
        print >> sys.stderr, 'scratch directory:'
        print >> sys.stderr, opts.scratch_dir;
        print >> sys.stderr, ''
    
    return opts

def get_trigs_txt(GWcursor,trigger_file,segs,min_thresh,tracker,verbose):
    """
    Read trigger data from text file.
    It is assumed that the first column is GPS time and the second is SNR.
    """
    # read lines avoiding white space and comment lines
    trigs = [line.split() for line in open(trigger_file) if (line.split() != [] and not line.startswith("%") and not line.startswith("#"))]

    # check if there is triggers
    if len(trigs) == 0:
      print >> sys.stderr, "Error: no triggers found. Please check your trigger file %s."%trigger_file
      sys.exit(1)
    
    # check the number of columns
    if len(trigs[0]) != 2:
        print >> sys.stderr, "Error: two columns are assumed (1. GPS time 2. SNR), found %d columns. Please check the trigger file."%len(trigs[0])
        sys.exit(1)       
    
    ## get times and snrs
    for trig in trigs:
        t = float(trig[0])
        s = float(trig[1])
        if t in segs:
          if s >= min_thresh and s != float('inf'):
            # insert into the database
            # micro second for GPS time is accurate enough
            GWcursor.execute("insert into GWtrigs values (?, ?)",("%.3f"%t,"%.2f"%s))
            tracker['counter'] += 1
        else:
          # keep track if some triggers are outside of segments and later
          # warn the user (for veto code use)
          tracker['outside'] = True
    return GWcursor

def get_trigs_xml(GWcursor,trigger_file,segs,min_thresh,tracker,verbose):
    """
    Read trigger data from xml file that has ligolw xml format.
    Many lines in this function are adapted from ligolw_print by Kipp Cannon
    and modified to suit the need.
    """
    from glue.ligolw import ligolw
    from glue.ligolw import table
    from glue.ligolw import lsctables
    from glue.ligolw import utils
    # speed hacks
    # replace Glue's pure Python LIGOTimeGPS class with pyLAL's C version
    from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
    lsctables.LIGOTimeGPS = LIGOTimeGPS 

    # Enable column interning to save memory
    table.RowBuilder = table.InterningRowBuilder

    # for now, hardcode the table/column names
    # FIXME: assuming there is only one of those tables in the file
    tables = ['sngl_burst','sngl_inspiral','sngl_ringdown']
    columns = ['end_time','end_time_ns','snr']
    

    # don't parse other tables so as to improve parsing speed and reduce memory 
    # requirements.
    
    def ContentHandler(xmldoc):
        return ligolw.PartialLIGOLWContentHandler(xmldoc, lambda name, attrs:\
                           (name == ligolw.Table.tagName) and\
                           (table.StripTableName(attrs["Name"]) in tables))
    try:
      lsctables.use_in(ligolw.PartialLIGOLWContentHandler)
    except AttributeError:
      # old glue did not allow .use_in().
      # FIXME:  remove when we can require the latest version of glue
      pass

    # FIXME: find a way to load only columns necessary
    # something like lsctables.SnglInspiral.loadcolumns = columns?

    xmldoc = utils.load_url(trigger_file, verbose = verbose,\
                                 gz = trigger_file.endswith(".gz"),
                                 contenthandler = ContentHandler)
    table.InterningRowBuilder.strings.clear()
    for table_elem in xmldoc.getElements(lambda e:\
                                        (e.tagName == ligolw.Table.tagName)):
      # trigger specific time retrieval functions
      if table_elem.tableName[:-6] in ('sngl_inspiral'):
        get_time = lambda row: row.get_end()
      elif table_elem.tableName[:-6] in ('sngl_burst'):
        get_time = lambda row: row.get_peak()
      elif table_elem.tableName[:-6] in ('sngl_ringdown'):
        get_time = lambda row: row.get_start()
      else:
        print >> sys.stderr, "Error: This should not be happening. Please contact to the author with the error trace."
        sys.exit(1)

      for row in table_elem:
        t = get_time(row)
        if t in segs:
          if row.snr > min_thresh:
            # insert into the database
            # micro second for GPS time is accurate enough
            GWcursor.execute("insert into GWtrigs values (?, ?)",("%.3f"%t, "%.2f"%row.snr))
            tracker['counter'] += 1
        else:
          # some triggers are outside of segments
          tracker['outside'] = True
    xmldoc.unlink()
    return GWcursor
      
def get_trigs(trigger_files,segs,min_thresh,name_tag=None,scratch_dir=".",verbose=True):
    """
    prepare the SQL database and retrieve the triggers    
    """
    # initialize SQLite database
    start_time = segs[0][0]
    end_time = segs[-1][1]
    prefix = "%d_%d_%d_GWtrigs"%(start_time,end_time,min_thresh)
    if name_tag != None:
      prefix = name_tag + "_" + prefix
    dbname = prefix+".db"
    # if the name already exists, rename the old one to avoid collisions
    KW_veto_utils.rename(dbname)
 
    global GW_working_filename # so that it can be erased when error occurs
    GW_working_filename = KW_veto_utils.get_connection_filename(\
                         dbname,tmp_path=scratch_dir,verbose=verbose)
   
    GWconnection = sqlite3.connect(GW_working_filename)
    GWcursor = GWconnection.cursor()

    GWcursor.execute('create table GWtrigs (GPSTime double, SNR double)')

    ## find out file format and read in the data
    # tracker tracks 1) number of triggers retrieved and 2) if any triggers outside of specified segments
    tracker = {'counter': 0,'outside': False}
    for t in trigger_files:
      ext=os.path.splitext(t)[-1]
      if ext == '.txt' or ext == '.dat':
        if verbose: print >> sys.stderr, "getting triggers from txt/dat file..."
        GWcursor = get_trigs_txt(GWcursor,t,segs,min_thresh,tracker,verbose)
      elif ext == '.xml':
        if verbose: print >> sys.stderr, "getting triggers from xml file..."
        GWcursor = get_trigs_xml(GWcursor,t,segs,min_thresh,tracker,verbose)
      else:
        print >> sys.stderr, """
        Error: unrecognized file format: please see --help and modify the 
               file to a supported format
        """
        sys.exit(1)
        
    GWconnection.commit()

    if verbose: print >> sys.stderr, "times and SNR for triggers retrieved!"

    # if some triggers were outside of given segments, warn the user
    if tracker['outside']:
      print >> sys.stderr, """
      Warning: Some of the triggers are outside of the segment list.
               Unless intentional (using DQ flags etc.), make sure you
               are using the right segment list.
               Ignoring...
      """

    # check how many triggers remained after cuts
    if tracker['counter'] == 0:
      print >> sys.stderr, """
      Error : No triggers remained after cut. 
              Please check trigger files and segments.
      """
      sys.exit(1)   
        
    if verbose:
      print >> sys.stderr, "%d triggers are retrieved."%tracker['counter']
 
    return dbname, GW_working_filename, GWconnection, GWcursor, tracker['counter']

# =============================================================================
#
#                                   Main
#
# =============================================================================

def main():
    # parse commandline
    opts = parse_commandline()

    ## get the segments on which triggers are retrieved
    # case 1: segment file is given   
    if opts.segment_file is not None:
      # read segment file
      seg_list = KW_veto_utils.read_segfile(opts.segment_file)
    # case 2: start and end GPS time are given
    else:
      seg_list = segmentlist([segment(opts.gps_start_time,opts.gps_end_time)])
    
    try:
      dbname, GW_working_filename, GWconnection, GWcursor, GWnum =\
         get_trigs(opts.trigger_files, seg_list, opts.min_thresh, name_tag = opts.name_tag, scratch_dir=opts.scratch_dir, verbose=opts.verbose)

      # save/display the result
      outname = os.path.splitext(dbname)[0]+opts.out_format
      KW_veto_utils.save_db(GWcursor,"GWtrigs",outname,GW_working_filename,
                      order_by=opts.order_by,verbose=opts.verbose)
   
      # close the connection to the database
      GWconnection.close()
    finally:
      # erase temporal database
      if globals().has_key('GW_working_filename'):
        db = globals()["GW_working_filename"]
        if opts.verbose:
          print >> sys.stderr, "removing temporary workspace '%s'..." % db
        os.remove(db)

if __name__=="__main__":
    main()

