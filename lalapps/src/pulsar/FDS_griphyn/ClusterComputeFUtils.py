"""
This module contains some utility routines used by the
ClusterComputeF.py script.
"""

__author__ = 'Teviet Creighton <tdcreigh@ligo.caltech.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

import os
import re

def listyears( dirlist, base ):
  """
  Determine which years or year ranges are covered by ephemeris files
  within a list of files.

  dirlist = list of files in a directory (or directories)
  base = basename for ephemeris files (either "earth" or "sun")

  Returns a pair of lists ( years, ranges ), where:

  years = list of integers [ y, y, ... ] representing the years of
          ephemeris filenames of the form baseYY.dat

  ranges = list of pairs of integers [ [y1,y2], [y1,y2], ... ]
           representing the years of ephemeris filenames of the form
           baseYY-YY.dat

  In each case, the two-digit year specifier YY, ranging from 00 to
  99, is converted into an integer year in the range 1950 to 2049
  (where YY matches the last two digits of the actual year number).
  """
  years = []
  ranges = []
  for name in dirlist:
    # Match single-year names and append year to list
    thematch = re.match( r"%s(\d\d)\.dat" % base, name )
    if thematch:
      y1 = int( thematch.group(1) ) + 1900
      if y1 < 1950:
	y1 = y1 + 100
      years = years + [y1]
    # Match year-range names and append year to list
    thematch = re.match( r"%s(\d\d)-(\d\d)\.dat" % base, name )
    if thematch:
      y1 = int( thematch.group(1) ) + 1900
      y2 = int( thematch.group(2) ) + 1900
      if y1 < 1950:
	y1 = y1 + 100
      if y2 < 1950:
	y2 = y2 + 100
      ranges = ranges + [ [y1,y2] ]
  return years, ranges


def gps2year( gps ):
  """
  Converts a GPS time into an integer year containing that time.  Leap
  seconds are ignored, so this routine may give incorrect results for
  times within a minute of a year change.

  gps = GPS time in seconds

  Returns an integer year of the Common Era.
  """
  t = gps/86400 + 36      # complete days since start of 1980
  t4 = t / 1461           # complete 4-year intervals since start of 1980
  t1 = ( t % 1461 ) / 365 # complete years since last leap year
  return 1980 + 4*t4 + t1 # year number of GPS start time


def getyearstr( dirname, start, end ):
  """
  Determines the correct year or year range string to use to get the
  right ephemeris file from a directory, corresponding to the
  specified start and end times.  If start and end times lie within a
  given year, the routine gives preference to a single-year ephemeris
  (if one is available for that year) over an ephemeris covering a
  range of years.

  dirname = name of directory containing ephemeris files
  start = initial GPS time in seconds
  end = final GPS time in seconds

  Returns a string of the form "YY" or "YY-YY", where YY is the last
  two digits of a year or range of years covering the requested
  interval, for which both earth and sun ephemeris files are
  available.  If no such string could be found, "" is returned.
  """
  # Get start and stop years
  ti = gps2year( start )
  tf = gps2year( end )
    
  # Get lists of years and ranges available for Earth and Sun
  dirlist = os.listdir( dirname )
  earthyears, earthranges = listyears( dirlist, "earth" )
  sunyears, sunranges = listyears( dirlist, "sun" )
  del dirlist
    
  # Restrict to those years or ranges common to Earth and Sun
  years = []
  ranges = []
  for y1 in earthyears:
    for y2 in sunyears:
      if y1 == y2:
	years = years + [y1]
  for y1 in earthranges:
    for y2 in sunranges:
      if y1 == y2:
	ranges = ranges + [y1]
  del earthyears
  del sunyears
  del earthranges
  del sunranges
    
  # Look for suitable year or range, giving precedence to
  # single-year ephemeris files.
  y = ""
  if ti == tf:
    for year in years:
      if ti == year:
	y = "%02i" % ( year % 100 )
  if y == "":
    for year in ranges:
      if ti >= year[0] and tf <= year[1]:
	y = "%02i-%02i" % ( year[0] % 100, year[1] % 100 )
  del years
  del ranges
  return y


class skyPatches:
  """
  This class stores a list of sky patch descriptors to be passed as
  command-line arguments to ComputeFStatistic.  These descriptors
  should be strings of the form:

    -R (alpha1,delta1),...,(alphaN,deltaN)

  or

    -a alpha -d delta -z Dalpha -c Ddelta

  In the former case, the sequence of coordinate pairs cannot contain
  whitespace, even if surrounded by quotation marks: in fact,
  double-quotes " will cause errors when the patch is passed to the
  command line of a ComputeFStatistic Condor job.
  """
  #"
  def __init__( self ):
    self.list = []

  def read( self, listfile, liststart, num ):
    """
    Reads (non-blank) lines from the named file, storing them in the
    sky patch list (self.list).

    listfile = name of file containing list of patches
    liststart = number of patches to skip from start of file
    num = number of patches to read from file

    If num is negative, all the remaining lines in the file are read.
    If num is greater than the number of lines remaining, an EOFError
    is raised, but self.list contains all the lines read.
    """
    fp = open( listfile, "r" )
    if liststart < 0:
      liststart = 0
    for i in range( 0, liststart ):
      fp.readline()
    if num >= 0:
      for i in range( 0, num ):
	line = fp.readline()
	if len( line ) > 0:
	  if line[-1] == '\n':
	    line = line[:-1]
	  self.list = self.list + [line]
	else:
	  i = num
      if num > len( self.list ):
	fp.close()
	raise EOFError, "EOF before %i lines read" % num
    else:
      i = 1
      while i:
	line = fp.readline()
	if len( line ) > 0:
	  if line[-1] == '\n':
	    line = line[:-1]
	  self.list = self.list + [line]
	else:
	  i = 0
    fp.close()

  def check( self ):
    """
    Checks the syntax of all patch descriptors in self.list.  At
    present, it simply checks for quotation marks, and raises a
    UserWarning exception if it finds any.
    """
    quotes = False
    for line in self.list:
      if not quotes:
	for c in line:
	  if c == '"' or c == "'":
	    quotes = True
    if quotes:
      raise UserWarning, "Quotation marks found in sky patch descriptors"
