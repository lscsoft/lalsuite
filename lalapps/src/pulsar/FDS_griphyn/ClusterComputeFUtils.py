import os
import re

# Routine to list years and year ranges of ephemeris files: dirlist is
# a listing of all filenames in a directory, and base is the part of
# the filename preceding year information (usually "earth" or "sun").
# Ephemeris files are assumed to have names of the form "%s%02i.dat"
# or "%s%02i-%02i.dat", where the first string is the base, and the
# one or two integers are years or year ranges.  Years in the range 00
# to 49 are assumed to be 2000 to 2049, while years in the range 50 to
# 99 are assumed to be 1950 to 1999.  Two objects are returned: a list
# of years [ y, y, ... ] corresponding to single-year ephemeris files,
# and a list of year pairs [ [y1,y2], [y1,y2], ... ] corresponding to
# year range ephemeris files.

def listyears( dirlist, base ):
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


# Routine to convert a GPS time into a year number.  Works for years
# from 1901 to 2099.  Ignores leap seconds.

def gps2year( gps ):
	t = gps/86400 + 36      # complete days since start of 1980
	t4 = t / 1461           # complete 4-year intervals since start of 1980
	t1 = ( t % 1461 ) / 365 # complete years since last leap year
	return 1980 + 4*t4 + t1 # year number of GPS start time


# Routine to compute the correct year or year range string to use to
# get the right ephemeris files: dirname is the directory to look for
# ephemeris files, and start and end are GPS times specifying the
# duration of observation.  The routine returns a string, which can be
# applied to the format "earth%s.dat" or "sun%s.dat" to get the
# ephemeris file names.  Uses listyears to identify the year or year
# ranges available; returns an empty string "" if no suitable
# ephemeris files were found.

def getyearstr( dirname, start, end ):
    # Get start and stop years
    ti = gps2year( start )
    tf = gps2year( end )
    
    # Get lists of years and year ranges available for Earth and Sun
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
