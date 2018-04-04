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
Tomoki Isogai (isogait@carleton.edu)

Utility functions for KW_veto codes

"""

from __future__ import division

import os
import sys
from math import floor, ceil
import cPickle
import glob
import shutil
import tempfile

try:
    import sqlite3
except ImportError:
    # pre 2.5.x
    from pysqlite2 import dbapi2 as sqlite3

from glue.segments import segment, segmentlist
from glue import segmentsUtils


# =============================================================================
# 
#                           Database Utility
#
# =============================================================================

def get_connection_filename(filename, tmp_path = None, replace_file = False, verbose = False):
  """
  copied from glue.ligolw.dbtables by Kipp Cannon
  I can't directly import from glue because dbtables collides with other tables
 
  Utility code for moving database files to a (presumably local)
  working location for improved performance and reduced fileserver
  load.
  """
  def mktmp(path, verbose = False):
    if not os.path.isdir(path):
      os.makedirs(path)
    fd, filename = tempfile.mkstemp(suffix = ".sqlite", dir = path)
    os.close(fd)
    if verbose:
      print >>sys.stderr, "using '%s' as workspace" % filename
    return filename

  def truncate(filename, verbose = False):
    if verbose:
      print >>sys.stderr, "'%s' exists, truncating ..." % filename
    try:
      fd = os.open(filename, os.O_WRONLY | os.O_TRUNC)
    except:
      if verbose:
        print >>sys.stderr, "cannot truncate '%s': %s" % (filename, str(e))
      return
    os.close(fd)
    if verbose:
      print >>sys.stderr, "done."

  def cpy(srcname, dstname, verbose = False):
    if verbose:
      print >>sys.stderr, "copying '%s' to '%s' ..." % (srcname, dstname)
    shutil.copy(srcname, dstname)

  database_exists = os.access(filename, os.F_OK)

  if tmp_path is not None:
    target = mktmp(tmp_path, verbose)
    if database_exists:
      if replace_file:
        # truncate database so that if this job
        # fails the user won't think the database
        # file is valid
        truncate(filename, verbose = verbose)
      else:
        # need to copy existing database to work
        # space for modifications
        i = 1
        while True:
          try:
            cpy(filename, target, verbose)
          except IOError, e:
            import errno
            import time
            if e.errno == errno.ENOSPC:
              if i < 5:
                if verbose:
                  print >>sys.stderr, "Warning: attempt %d: no space left on device, sleeping and trying again ..." % i
                time.sleep(10)
                i += 1
                continue
              else:
                if verbose:
                  print >>sys.stderr, "Warning: attempt %d: no space left on device: working with original file" % i
              os.remove(target)
              target = filename
            else:
              raise e
          break
  else:
    target = filename
    if database_exists and replace_file:
      truncate(target, verbose = verbose)

  del mktmp
  del truncate
  del cpy

  return target

# =============================================================================
#
#                         Reading and Writing Utility
#
# =============================================================================

def save_db(cursor, table_name, filename, working_filename, order_by = None, verbose=True):
    """
    Originally written by Nickolas Fotopoulos, modified.
    Convert a dictionary to SQLite database for later use and also save to
    specified format file.
    Supported file extensions are:
    * .pickle - Python pickle file (dictionary serialized unchanged)
    * .pickle.gz - gzipped pickle (dictionary serialized unchanged)
    * .mat - Matlab v4 file (dictionary keys become variable names; requires
             Scipy)
    * .txt - ASCII text 
    * .txt.gz - gzipped ASCII text
    * .bd - sqlite database

    FIXME: add xml
    """
    if verbose:
      print >> sys.stderr, "saving data in %s..." % filename

    if filename == '':
        raise ValueError, "Empty filename"

    ext = os.path.splitext(filename)[-1]

    ## .db (SQLite) file (save all table no matter what specified)
    if ext == '.db':
      shutil.copy(working_filename,filename) 
      return

    # FIXME: in theory it's suceptible to SQL injection attack
    cmd = "select * from %s"%table_name
    if order_by is not None:
      cmd += " order by %s"%order_by
    table = cursor.execute(cmd)

    ## no save -> print out in stdout
    if ext == '.stdout':
      print "\n".join(["# column %d: "%(i+1)+c[0] for i, c in enumerate(cursor.description)])
      for row in table:
        print "\t".join(map(str,row))
      return
    
    ## Set up file_handle
    file_handle = file(filename, 'wb')
 
    # For gz files, bind file_handle to a gzip file and find the new extension
    if ext == '.gz':
        import gzip
        file_handle = gzip.GzipFile(fileobj=file_handle, mode='wb')
        ext = os.path.splitext(filename[:-len(ext)])[-1]
 
    if ext == '.txt':
      lines = []
      lines.append("\n".join(["# column %d: "%(i+1)+c[0] for i, c in enumerate(cursor.description)])+"\n")
      for row in table:
        lines.append("\t".join(map(str,row)))
      file_handle.write("\n".join(lines))
      return 

    if ext == '.pickle' or '.mat':
      # create dictionary whose keys correspond to columns
      columns = [c[0] for c in cursor.description]
      # initialize
      data_dict = {}
      for c in columns:
        data_dict[c] = []
      # insert data
      for row in table:
        for c, r in zip(columns,row):
          data_dict[c].append(r)
        
      if ext == '.mat':
        import scipy.io
        scipy.io.savemat(filename, data_dict)
        return
 
      if ext == '.pickle':
        import cPickle
        output = cPickle.dumps(data_dict, protocol = 2)
        file_handle.write(output)
        return    

      else:
        raise ValueError, "Unrecognized file extension"

def load_db(dbFile, tmp_path, verbose):
  """
  Connect to database.
  Also, return a dictionary having variables from KW_veto_calc in a form:
  {variable name: variable value}
  """
  if verbose:
    print >> sys.stderr, "connecting to database %s..."%dbFile
  working_filename = get_connection_filename(\
                       dbFile, tmp_path=tmp_path, verbose=verbose)
  connection = sqlite3.connect(working_filename)
  cursor = connection.cursor()
  params = {}
  for v in cursor.execute("select var_name, value from params").fetchall():
    params[v[0]] = v[1]

  return cursor, connection, working_filename, params

def get_candidate(cursor, critical_usedPer):
  """
  This function returns  a tuple:
  (threshold, used percentage, # of coincident KW triggers above the threshold, # of total KW triggers above the threshold, # of vetoed GW triggers, veto efficiency, dead time, dead time percentage, (used percentage) / (random used percentage), (veto efficiency) / (dead time percentage)) that corresponds to threshold to be used for veto.
  This returns just None for non candidate channel.
  """
  return cursor.execute("select * from result where usedPercentage > ? order by threshold asc",(critical_usedPer,)).fetchone()
   

def rename(src):
  """
  If src alaready exists, this function rename it so that new file/directory
   won't overwite
  """
  index = 0; dst = os.path.normpath(src)
  while os.path.exists(dst):
    index += 1
    dst = "%s.%d"%(os.path.normpath(src),index)
  if index > 0:
    print >> sys.stderr, "Warning: %s already exists, renaming the existing one to %s and continuing..."%(src,dst)
    os.renames(src,dst)

def load_segs_db(cur,tableName):
  """
  Load segments in glue.segments.segmentlist format
  FIXME: possibly SQL Injection attack
  """
  seg_list = segmentlist([segment(s[0], s[1]) for s in\
                cur.execute("select start_time, end_time from %s"%tableName)])

  return seg_list.coalesce()

def write_segs_db(cur,seg_list,tableName):
  """
  Write glue.segments.segmentlist format segments into database.
  FIXME: possibly SQL Injection attack
  """
  seg_list.coalesce()
  cur.execute("create table %s (start_time int, end_time int)"%tableName)
  cur.executemany("insert into %s values (?,?)"%tableName,seg_list)

  return cur

def write_segs(segs,file_name):
  """
  Take glue.segments.segmentlist and write in a ascii file in a form compatible
  to the segment database.
  """
  open(file_name,'w').write("\n".join(["%d %d"%(s[0],s[1]) for s in segs]))

# read segment file
def read_segfile(segfile):
  """
  Read segments in segwizard form with sanity check and return
  glue.segments.segmentlist
  """
  try:
    seg_list =\
          segmentsUtils.fromsegwizard(open(segfile),coltype=float,strict=False)
  except(IOError):
    print >> sys.stderr, """
  Error: file contains no segments or glue.segmentsUtils.fromsegwizard is
         not reading segments correctly. Please check the seg file. 
         (possibly comments in the file is causing this)
         %s
    """%segfile
    raise

  seg_list.coalesce()
  if seg_list == [] or abs(seg_list) == 0:
    print >> sys.stderr, """
    Warning: File contains no segments or glue.segmentsUtils.fromsegwizard is
             not reading segments correctly. Please check the seg file.   
             (Possibly comments in the file is causing this.)
             Attempting to ignore.
    """
  return seg_list

def read_segfile_xml(segfile,verbose):
  """
  Read segment file in ligolw xml type and return in glue.segments.segmentlist
  format.
  """
  from glue.ligolw import ligolw
  from glue.ligolw import table
  from glue.ligolw import utils

  def ContentHandler(xmldoc):
    return ligolw.PartialLIGOLWContentHandler(xmldoc, lambda name, attrs:\
               (name == ligolw.Table.tagName) and\
               (table.StripTableName(attrs["Name"]) in ["segment"]))
  try:
    table.use_in(ligolw.PartialLIGOLWContentHandler)
  except AttributeError:
    # old glue did not allow .use_in().
    # FIXME:  remove when we can require the latest version of glue
    pass

  xmldoc = utils.load_url(segfile, verbose = verbose,gz = segfile.endswith(".gz"), contenthandler = ContentHandler)
  seg_list = segmentlist()
  for table_elem in xmldoc.getElements(lambda e:\
                                       (e.tagName == ligolw.Table.tagName)):
    for row in table_elem:
      seg_list.append(segment(row.start_time, row.end_time))
  xmldoc.unlink()
  return seg_list

def find_version_xml(segfile,seg,verbose):
  """
  Find out the version of the flag for the given seg.
  """
  from glue.ligolw import ligolw
  from glue.ligolw import table
  from glue.ligolw import utils

  def ContentHandler(xmldoc):
    return ligolw.PartialLIGOLWContentHandler(xmldoc, lambda name, attrs:\
               (name == ligolw.Table.tagName) and\
               (table.StripTableName(attrs["Name"]) in ["segment_definer","segment_summary"]))
  try:
    table.use_in(ligolw.PartialLIGOLWContentHandler)
  except AttributeError:
    # old glue did not allow .use_in().
    # FIXME:  remove when we can require the latest version of glue
    pass

  xmldoc = utils.load_url(segfile, verbose = verbose,gz = segfile.endswith(".gz"), contenthandler = ContentHandler)
  for n, table_elem in enumerate(xmldoc.getElements(lambda e:\
                                       (e.tagName == ligolw.Table.tagName))):
    if n == 0:
      definer = {}
      for row in table_elem:
        if row.name != "RESULT":
          definer[str(row.segment_def_id).split(":")[-1]] = row.version

    if n == 1:
      for row in table_elem:
        if seg[0] >= row.start_time and seg[1] <= row.end_time:
          if str(row.segment_def_id).split(":")[-1] in definer.keys():
            xmldoc.unlink()
            return definer[str(row.segment_def_id).split(":")[-1]]


def get_segment(eventTime, pos_window, neg_window):
  """ 
  Given the time of trigger, this aplies +/- window and round to integers
  in a way it broadens the window
  """ 
  return segment(floor(eventTime - neg_window), ceil(eventTime + pos_window))

def get_result_files(result_dir,verbose):
  """
  Make a list of the result files from KW_veto_calc and check sanity.
  """
  files_list = [os.path.join(result_dir,f) for f in os.listdir(result_dir) if f.endswith("-data.db")]
  if verbose:
    print >> sys.stderr, "result files:", files_list
  if files_list==[]: # check if there is at least one file
    print >> sys.stderr, "Error: no result files found in %s"%result_dir
    sys.exit(1)

  # check if they have the same name tag
  if os.path.commonprefix(files_list)=="":
    print >> sys.stderr, """
    Error: Possibly files from different runs are coexisting. 
           Please check the result_dir.
    """
    sys.exit(1)

  return files_list

def get_files_from_globs(globs):
  """
  Get individual file names from names with globs.
  Input could be a mixture of glob and ',' seperated file list, something like
  '../segs/*,H1_segs.txt' 
  also there might be white space between file names 
  """
  file_globs = [f.strip() for f in globs.split(",") if f.strip()!=""]
  all_files = [] 
  for f in file_globs: 
    files = glob.glob(f) 
    if files == []: 
      print >> sys.stderr, "Error: %s not found"%f 
      sys.exit(1) 
    else: 
      all_files.extend(files) 
               
  return all_files
