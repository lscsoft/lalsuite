"""
  Matt Pitkin - Copyright 2013

  A script to convert a time of arrival of a signal from a given source location
  at the solar system barycentre (given in MJD) to the GPS time of arrival at a
  given site on Earth (these are assumed to be either standard GW detectors or a
  subset of radio telescopes.)
"""

import lal
import lalpulsar

from optparse import OptionParser

import sys
import math
import os
import numpy as np

__author__="Matthew Pitkin <matthew.pitkin@ligo.org>"
__version__="1.0"

# set some radio telescopes
radscopenames = ['PKS', # Parkes
                 'GBT', # Green Bank Telescope
                 'JBO', # Lovell Telescope at Jodrell Bank Observatory
                 'AO', # Arecibo Telescope
                 'EFF', # Effelsberg
                 'NCY', # Nancay
                 'HOB', # Hobart
                 'HART', # Hartebeesthoek
                 'VLA', # Very Large Array
                 'WSRT'] # Westerbork

# x, y, z coordinates of telescopes in metres taken from TEMPO2 observatories.dat file
radscopepos = {'PKS': np.array([ -4554231.5, 2816759.1, -3454036.3 ]), \
               'GBT': np.array([ 882589.65, -4924872.32, 3943729.348 ]), \
               'JBO': np.array([ 3822626.04, -154105.65, 5086486.04 ]), \
               'AO': np.array([ 2390490.0, -5564764.0, 1994727.0 ]), \
               'EFF': np.array([ 4033949.5, 486989.4, 4900430.8 ]), \
               'NCY': np.array([ 4324165.81, 165927.11, 4670132.83]), \
               'HOB': np.array([ -3950077.96, 2522377.31, -4311667.52 ]), \
               'HART': np.array([ 5085442.780, 2668263.483, -2768697.034 ]), \
               'VLA': np.array([ -1601192.0, -5041981.4, 3554871.4 ]), \
               'WSRT': np.array([ 3828445.659, 445223.600, 5064921.5677 ])}

if __name__=='__main__':
  description = """This script will return the GPS time of arrival at a given Earth-based detector given a \
MJD time of arrival of a signal from a particular sky position at the solar system barycentre"""

  usage = "Usage: %prog [options]"

  parser = OptionParser( usage = usage, description = description,
                         version = __version__ )

  parser.add_option("-t", "--mjdtime", dest="mjdtime", help="The MJD time of arrival at the SSB", type="float")

  parser.add_option("-r", "--right-ascension", dest="ra", help="The right ascension of the source in the format \
hh:mm:ss.s")

  parser.add_option("-d", "--declination", dest="dec", help="The declination of the source in the format \
dd:mm:ss.s")

  parser.add_option("-z", "--detector", dest="det", help="The gravitational wave detector or radio telescope acronym \
e.g. H1")

  # parse input options
  (opts, args) = parser.parse_args()

  # check options
  if not opts.__dict__['mjdtime']:
    print >> sys.stderr, "Must specify an MJD time of arrival"
    parser.print_help()
    sys.exit(1)
  else:
    mjd = opts.mjdtime

  if not opts.__dict__['ra']:
    print >> sys.stderr, "Must specify a source right ascension"
    parser.print_help()
    sys.exit(1)
  else:
    ra = opts.ra

  if not opts.__dict__['dec']:
    print >> sys.stderr, "Must specify a source declination"
    parser.print_help()
    sys.exit(1)
  else:
    dec = opts.dec

  if not opts.__dict__['det']:
    print >> sys.stderr, "Must specify a detector"
    parser.print_help()
    sys.exit(1)
  else:
    det = opts.det

  # get the path for the ephemeris files
  try:
    ephpath = os.environ['LALPULSAR_PREFIX']
  except:
    print >> sys.stderr, "Need environment variables 'LALPULSAR_PREFIX' to be set"
    sys.exit(1)

  earth = ephpath + '/share/lalpulsar/earth00-19-DE405.dat.gz'
  sun = ephpath + '/share/lalpulsar/sun00-19-DE405.dat.gz'

  # check files exist
  if not os.path.isfile(earth) or not os.path.isfile(sun):
    print >> sys.stderr, "Earth or Sun ephemeris file does not exist"
    sys.exit(1)

  # get ephemerides
  edat = lal.InitBarycenter( earth, sun )

  # convert MJD time at SSB to GPS time at SSB
  gpstime = lalpulsar.TTMJDtoGPS( mjd )

  # check time is in range of ephemeris files
  ephemstart = 630720013 # GPS time of Jan 1, 2000, 00:00:00 UTC
  ephemend = 1261872015  # GPS time of Jan 1, 2020, 00:00:00 UTC

  if gpstime < ephemstart or gpstime > ephemend:
    print >> sys.stderr, "Times (GPS %.9f) is outside the ephemeris file ranges!\n" % gpstime
    sys.exit(1)

  # convert time to LIGOTimeGPS structure
  gpsin = lal.LIGOTimeGPS()
  lal.GPSSetREAL8( gpsin, gpstime )

  # convert RA and dec to radians
  rarad = lalpulsar.hmsToRads( ra )
  decrad = lalpulsar.dmsToRads( dec )

  # setup a pulsar signal params structure
  params = lalpulsar.PulsarSignalParams()
  params.pulsar.position.latitude = decrad
  params.pulsar.position.longitude = rarad
  params.pulsar.position.system = lal.COORDINATESYSTEM_EQUATORIAL
  params.ephemerides = edat

  # set the site information
  if det.upper() in radscopenames:
    detector = lal.Detector()
    detector.location = radscopepos[det.upper()]
    params.site = detector
  else:
    params.site = lalpulsar.GetSiteInfo( det.upper() );

  # get time of arrival at detector
  toadet = lal.LIGOTimeGPS()
  lalpulsar.ConvertSSB2GPS( toadet, gpsin, params )

  # output the time
  #print >> sys.stdout, "%.9f" % (toadet.gpsSeconds + toadet.gpsNanoSeconds*1e-9)
  print >> sys.stdout, "%d.%09d" % (toadet.gpsSeconds, toadet.gpsNanoSeconds)

  sys.exit(0)