#!/usr/bin/env python

# Copyright (C) 2011 Duncan M. Macleod
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

# =============================================================================
# Preamble
# =============================================================================

from __future__ import division
import numpy as np,re,sys

import lal

from glue import git_version
from glue.ligolw import table

from pylal import inject,date
from pylal import inject
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS

__author__ = "Duncan M. Macleod <duncan.macleod@astro.cf.ac.uk>"
__version__ = git_version.id
__date__ = git_version.date

# useful regex
_comment_regex = re.compile(r"\s*([#;].*)?\Z", re.DOTALL)
_delim_regex   = re.compile('[,\s]')

# =============================================================================
# Sky Positions structure
# =============================================================================

class SkyPositionTable(table.Table):
    """
    glue.ligolw.table.Table holding SkyPosition objects.
    """

    tableName     = "sky_positions:table"
    valid_columns = {"process_id":  "ilwd:char",\
                     "latitude":    "real_4",\
                     "longitude":   "real_4",\
                     "probability": "real_4",\
                     "solid_angle": "real_4",\
                     "system":      "lstring"}

    def rotate(self, R):
        rotated = table.new_from_template(self)
        for row in self:
            rotated.append(row.rotate(R))
        return rotated

    def normalize(self):
        for row in self:
            row.normalize()

    def parseTimeDelayDegeneracy(self, ifos, gpstime=lal.GPSTimeNow(),\
                                 dt=0.0005):
     
        # get detectors
        detectors = [inject.cached_detector.get(inject.prefix_to_name[ifo])\
                         for ifo in ifos]
    
        timeDelays = []
        degenerate = []
    
        new = table.new_from_template(self)
        for n,row in enumerate(self):
            # get all time delays for this point
            timeDelays.append([])
            for i in xrange(len(ifos)):
                for j in xrange(i+1, len(ifos)):
                  timeDelays[n].append(date.XLALArrivalTimeDiff(\
                                           detectors[i].location,\
                                           detectors[j].location,\
                                           row.longitude,\
                                           row.latitude,\
                                           LIGOTimeGPS(gpstime)))
            # skip the first point
            if n==0:
                degenerate.append(False)
                new.append(row)
                continue
            else:
                degenerate.append(True)
            # test this point against all others
            for i in xrange(0,n):
                # if check point is degenerate, skip
                if degenerate[i]:
                    continue
                # check each time delay individually
                for j in xrange(0, len(timeDelays[n])):
                    # if this time delay is non-degenerate the point is valid
                    if np.fabs(timeDelays[i][j]-timeDelays[n][j]) >= dt:
                        degenerate[n] = False
                        break
                    else:
                        degenerate[n] = True
                if degenerate[n]:
                    break
    
            if not degenerate[n]:
                new.append(row)
    
        return new


class SkyPosition(object):

    __slots__ = SkyPositionTable.valid_columns.keys()

    def __iter__(self):
        return iter(tuple((self.longitude, self.latitude)))

    def __repr__(self):
        return "SkyPosition(%s, %s)" % (self.longitude, self.latitude)

    def __str__(self):
        return "(%s, %s)" % (self.longitude, self.latitude)

    def rotate(self, R):
        """
        Apply the 3x3 rotation matrix R to this SkyPosition.
        """

        cart = SphericalToCartesian(self)
        rot  = np.dot(np.asarray(R), cart)
        new  = CartesianToSpherical(rot, system=self.system)
        new.normalize()
        try:
            new.probability = self.probability
        except AttributeError:
            new.probability = None
        try:
            new.solid_angle = self.solid_angle
        except AttributeError:
            new.solid_angle = None
            
        return new

    def normalize(self):
        """
        Normalize this SkyPosition to be within standard limits
        [0 <= longitude < 2pi] and [-pi < latitude <= pi]
        """

        # first unwind positions, i.e. make sure that:
        # [0 <= alpha < 2pi] and [-pi < delta <= pi]

        # normalise longitude
        while self.longitude < 0:
            self.longitude += lal.TWOPI
        while self.longitude >= lal.TWOPI:
            self.longitude -= lal.TWOPI

        # normalise latitude
        while self.latitude <= -lal.PI:
            self.latitude += lal.TWOPI
        while self.latitude > lal.TWOPI:
            self.latitude -= lal.TWOPI

        #
        # now get latitude into canonical interval [-pi/2, pi/2)
        #

        if self.latitude > lal.PI_2:
            self.latitude = lal.PI - self.latitude
            if self.latitude < lal.PI:
                self.longitude += lal.PI
            else:
                self.longitude -= lal.PI

        if self.latitude < -lal.PI_2:
            self.latitude = -lal.PI - self.latitude
            if self.longitude < lal.PI:
                self.longitude += lal.PI
            else:
                self.longitude -= lal.PI

    def opening_angle(self, other):
        """
        Calcalate the opening angle between this SkyPosition and the other
        SkyPosition
        """

        if self == other:
            return 0.0

        s = np.sin
        c = np.cos

        return np.arccos(s(self.latitude) * s(other.latitude) +\
                         c(self.latitude) * c(other.latitude) *\
                         c(self.longitude - other.longitude))

    def get_time_delay(self, ifo1, ifo2, gpstime):
        """
        Return the time delay for this SkyPosition between arrival at ifo1
        relative to ifo2, for the given gpstime.
        """
        det1 = inject.cached_detector.get(inject.prefix_to_name[ifo1])
        det2 = inject.cached_detector.get(inject.prefix_to_name[ifo2])

        return date.XLALArrivalTimeDiff(det1.location, det2.location,\
                                        self.longitude, self.latitude,\
                                        LIGOTimeGPS(gpstime))



# =============================================================================
# Convert between coordinate systems
# =============================================================================

# a few last definitions (taken from SkyCoordinates.c)
lal.LAL_ALPHAGAL = 3.366032942
lal.LAL_DELTAGAL = 0.473477302
lal.LAL_LGAL     = 0.576

def ConvertSkyPosition(skyPoint, system, zenith=None, gpstime=None):
    """
    Convert the SkyPosition object skyPoint from it's current system to the
    new system.

    Valid systems are : 'horizon', 'geographic', 'equatorial', 'ecliptic',\
                        'galactic'

    SkyPosition zenith and gpstime should be given as appropriate for 'horizon'
    and 'geographic' conversions respectively. 
    """

    valid_systems = ['horizon', 'geographic', 'equatorial', 'ecliptic',\
                           'galactic']

    system = system.lower()
    skyPoint.system = skyPoint.system.lower()

    if system not in valid_systems:
        raise AttributeError("Unrecognised system='%s'" % system)

    while skyPoint.system != system:

        if skyPoint.system == 'horizon':
            skyPoint = HorizonToSystem(skyPoint, zenith)

        elif skyPoint.system == 'geographic':
            if system == 'horizon':
                if zenith.system == 'geographic':
                    skyPoint = SystemToHorizon(skyPoint, zenith)
                elif zenith.system == 'equatorial':
                    skyPoint = GeographicToEquatorial(skyPoint, gpstime)
                else:
                    raise AttibuteError("Cannot convert from geographic to "+\
                                        "horizon with zenith point not in "+\
                                        "geogrphic of equatorial systems.")
            else:
                skyPoint = GeographicToEquatorial(skyPoint, gpstime)

            if system == 'horizon' and zenith.system == 'equatorial':
                skyPoint = SystemToHorizon(skyPoint, zenith)
            elif system == 'ecliptic':
                skyPoint = EquatorialToEcliptic(skyPoint)
            elif system == 'galactic':
                skyPoint = EquatorialToGalactic(skyPoint)
            else:
                skyPoint = EquatorialToGeographic(skyPoint, gpstime)

        elif skyPoint.system == 'ecliptic':
            skyPoint = EclipticToEquatorial(skyPoint)

        elif skyPoint.system == 'galactic':
            skyPoint = GalacticToEquatorial(skyPoint)

    return skyPoint

def HorizonToSystem(input, zenith):
    """
    Convert the SkyPosition object input from 'horizon' to the inherited system
    using the SkyPosition zenith
    """

    if input.system.lower() != 'horizon':
        raise AttributeError("Cannot convert from horizon for point in "+\
                             "system='%s'" % input.system.lower())

    if zenith.system.lower() != 'equatorial'\
    and zenith.system.lower() != 'geographic':
        raise AttributeError("zenith must have coordinate system = "+\
                             "'equatorial' or 'geographic'")

    # intermediates
    sinP = np.sin(zenith.latitude)
    cosP = np.cos(zenith.latitude)
    sina = np.sin(input.latitude)
    cosa = np.cos(input.latitude)
    sinA = np.sin(input.longitude)
    cosA = np.cos(input.longitude)

    # final components
    sinD = sina*sinP + cosa*cosA*cosP
    sinH = cosa*sinA
    cosH = sina*cosP - cosa*cosA*sinP

    # output
    output             = SkyPosition()
    output.system      = zenith.system
    output.latitude    = np.arcsin(sinD)
    output.longitude   = zenith.longitude - np.arctan2(sinH, cosH)
    output.probability = input.probability
    output.solid_angle = input.solid_angle
    output.normalize()

    return output

def SystemToHorizon(input, zenith):
    """
    Convert the SkyPosition object input from the inherited system to 'horizon'
    using the SkyPosition zenith
    """

    if input.system.lower() != zenith.system.lower():
        raise AttributeError("input coordinate system must equal zenith system")

    if zenith.system.lower() != 'equatorial'\
    and zenith.system.lower() != 'geographic':
        raise AttributeError("zenith must have coordinate system = "+\
                             "'equatorial' or 'geographic'")

    # intermediates
    h = zenith.longitude - input.longitude
    sinH = np.sin(h)
    cosH = np.cos(h)
    sinP = np.sin(zenith.latitude)
    cosP = np.cos(zenith.latitude)
    sinD = np.sin(input.latitude)
    cosD = np.cos(input.latitude)

    # final components
    sina = sinD*sinP + cosD*cosP*cosH
    sinA = cosD*sinH
    cosA = sinD*cosP - cosD*sinP*cosH

    # output
    output = SkyPosition()
    output.system    = 'horizon'
    output.latitude  = np.arcsin(sina)
    output.longitude = np.arctan2(sinA, cosA)
    output.probability = input.probability
    output.solid_angle = input.solid_angle
    output.normalize()

    return output

def GeographicToEquatorial(input, gpstime):
    """
    Convert the SkyPosition object input from the inherited 'geographical'
    system to  to 'equatorial'.
    """

    if input.system.lower() != 'geographic':
        raise AttributeError("Input system!='geographic'")

    # get GMST
    gmst = np.fmod(date.XLALGreenwichSiderealTime(gpstime, 0),\
                   lal.TWOPI)

    # output
    output = SkyPosition()
    output.system = 'equatorial'
    output.longitude = np.fmod(input.longitude + gmst, lal.TWOPI)
    output.latitude = input.latitude
    output.probability = input.probability
    output.solid_angle = input.solid_angle
    output.normalize()

    return output

def EquatorialToEcliptic(input):
    """
    Convert the SkyPosition object input from the inherited 'equatorial'
    system to  to 'ecliptic'.
    """

    # intermediates
    sinD = np.sin(input.latitude)
    cosD = np.cos(input.latitude)
    sinA = np.sin(input.longitude)
    cosA = np.cos(input.longitude)

    # components
    sinB = sinD*np.cos(lal.IEARTH)\
                - cosD*sinA*np.sin(lal.IEARTH)
    sinL = cosD*sinA*np.cos(lal.IEARTH)\
                + sinD*np.sin(lal.IEARTH)
    cosL = cosD*cosA

    # output
    output.system    = 'ecliptic'
    output.latitude  = np.arcsin(sinB)
    output.longitude = np.arctan2(sinL, cosL)
    output.normalize()

    return output

def EquatorialToGalactic(input):
    """
    Convert the SkyPosition object input from the inherited 'equatorial'
    system to  to 'galactic'.
    """

    # intermediates. */
    a    = -lal.LAL_ALPHAGAL + input.longitude
    sinD = np.sin(input.latitude)
    cosD = np.cos(input.latitude)
    sinA = np.sin(a)
    cosA = np.cos(a)

    # components. */
    sinB = cosD*np.cos(lal.LAL_DELTAGAL)*cosA\
                + sinD*np.sin(lal.LAL_DELTAGAL)
    sinL = sinD*np.cos(lal.LAL_DELTAGAL)\
                - cosD*cosA*np.sin(lal.LAL_DELTAGAL)
    cosL = cosD*sinA

    # output
    output = SkyPosition()
    output.system    = 'galactic'
    output.latitude  = np.arcsin(sinB)
    output.longitude = np.arctan2(sinL, cosL) + lal.LAL_LGAL
    output.probability = input.probability
    output.solid_angle = input.solid_angle
    output.normalize()

    return output

def EquatorialToGeographic(input, gpstime):
    """
    Convert the SkyPosition object input from the inherited 'equatorial'
    system to  to 'geographic'.
    """

    if input.system.lower() != 'equatorial':
        raise AttributeError("Input system is not 'equatorial'")

    # get GMST
    gmst = np.fmod(date.XLALGreenwichSiderealTime(gpstime, 0),\
                   lal.TWOPI)

    # output
    output = SkyPosition()
    output.system = 'geographic'
    output.longitude = np.fmod(input.longitude - gmst, lal.TWOPI)
    output.latitude = input.latitude
    output.probability = input.probability
    output.solid_angle = input.solid_angle
    output.normalize()

    return output

def EclipticToEquatorial(input):
    """
    Convert the SkyPosition object input from the inherited 'eliptic'
    system to  to 'equatorial'.
    """

    # intermediates
    sinB = np.sin(input.latitude)
    cosB = np.cos(input.latitude)
    sinL = np.sin(input.longitude)
    cosL = np.cos(input.longitude)

    # components
    sinD = cosB*sinL*np.sin(lal.IEARTH)\
                + sinB*np.cos(lal.IEARTH)
    sinA = cosB*sinL*np.cos(lal.IEARTH)\
                - sinB*np.sin(lal.IEARTH)
    cosA = cosB*cosL

    # output
    output.system    = 'equatorial'
    output.latitude  = np.arcsin(sinD)
    output.longitude = np.arctan2(sinA, cosA)
    output.normalize()

    return output

# =============================================================================
# Read grid from file
# =============================================================================

def fromfile(fobj, loncol=0, latcol=1, probcol=None, angcol=None, system=None,\
             degrees=False):
    """
    Returns a SkyPositionTable as read from the file object fobj.
    """

    l = SkyPositionTable()

    if degrees:
        func = np.radians
    else:
        func = float

    for line in fobj:
        line = _comment_regex.split(line)[0]
        if not line: continue
        line = re.split(_delim_regex, line)
        p = SkyPosition()
        p.longitude   = func(float(line[loncol]))
        p.latitude    = func(float(line[latcol]))
        if probcol:
            p.probability = line[probcol]
        else:
            p.probability = None
        if angcol:
            p.solid_angle = line[angcol]
        else:
            p.solid_angle = None
        p.system      = system
        l.append(p)
        
    return l

def tofile(fobj, grid, delimiter=" ", degrees=False):
    """
    Writes the SkyPositionTable object grid to the file object fobj
    """

    for p in grid:
        if degrees:
            lon = np.degrees(p.longitude)
            lat = np.degrees(p.latitude)
        else:
            lon = p.longitude
            lat = p.latitude
        fobj.write("%s%s%s\n" % (lon, delimiter, lat))

# =============================================================================
# Generate grids
# =============================================================================

def SkyPatch(ifos, ra, dec, radius, gpstime, dt=0.0005, sigma=1.65,\
             grid='circular'):
    """
    Returns a SkyPositionTable of circular rings emanating from a given
    central ra and dec. out to the maximal radius.
    """

    # form centre point
    p = SkyPosition()
    p.longitude = ra
    p.latitude  = dec

    # get detectors
    ifos.sort()
    detectors  = []
    for ifo in ifos:
        if ifo not in inject.prefix_to_name.keys():
            raise ValueError("Interferometer '%s' not recognised." % ifo)
        detectors.append(inject.cached_detector.get(\
                             inject.prefix_to_name[ifo]))

    alpha = 0
    for i in xrange(len(ifos)):
        for j in xrange(i+1,len(ifos)):
            # get opening angle
            baseline = date.XLALArrivalTimeDiff(detectors[i].location,\
                                                detectors[j].location,\
                                                ra, dec, LIGOTimeGPS(gpstime))
            ltt      = inject.light_travel_time(ifos[i], ifos[j])
            angle    = np.arccos(baseline/ltt)
            
            # get window
            lmin = angle-radius
            lmax = angle+radius
            if lmin < lal.PI_2 and lmax > lal.PI_2:
                l = lal.PI_2
            elif np.fabs(lal.PI_2-lmin) <\
                     np.fabs(lal.PI_2-lmax):
                l = lmin
            else:
                l = lmax
 
            # get alpha
            dalpha = ltt * np.sin(l)
            if dalpha > alpha:
                alpha = dalpha
            
    # get angular resolution
    angRes = 2 * dt/alpha

    # generate grid
    if grid.lower() == 'circular':
        grid = CircularGrid(angRes, radius)
    else:
        raise RuntimeError("Must use grid='circular', others not coded yet")

    #
    # Need to rotate grid onto (ra, dec)
    #

    # calculate opening angle from north pole
    north = [0, 0, 1]
    angle = np.arccos(np.dot(north, SphericalToCartesian(p)))
    #angle = north.opening_angle(p)

    # get rotation axis
    axis = np.cross(north, SphericalToCartesian(p))
    axis = axis / np.linalg.norm(axis)

    # rotate grid
    R = _rotation(axis, angle)
    grid = grid.rotate(R)

    # 
    # calculate probability density function
    #

    # assume Fisher distribution in angle from centre
    kappa = (0.66*radius/sigma)**(-2)

    # compute probability
    for p in grid:
        overlap = np.cos(p.opening_angle(grid[0]))
        p.probability = np.exp(kappa*(overlap-1))

    probs = [p.probability for p in grid]
    for p in grid:
        p.probability = p.probability/max(probs)

    return grid

def CircularGrid(resolution, radius):
    """
    Generates a SkyPositionTable of circular grids around the North Pole with
    given angular resolution and a maximal radius.
    """

    l = SkyPositionTable()

    theta = 0
    while theta < radius+resolution:
        # get number of points in ring
        n  = max(1, int(np.ceil(lal.TWOPI\
                                * np.sin(theta)/resolution)))
        # get phi step
        dp = lal.TWOPI / n
        # get solid angle
        if theta == 0 or theta == lal.PI:
            dO = 4*lal.PI * np.sin(0.25*resolution)**2
        else:
            dO = 4*lal.PI/n * np.sin(theta) * np.sin(0.5*resolution)
 
        # assign points in ring
        for j in xrange(n):
            p = SkyPosition()
            p.system = 'equatorial'
            p.longitude = (-lal.PI + dp/2) + dp*j
            p.latitude  = lal.PI_2 - theta
            p.solid_angle = dO
            p.probability = 0
            p.normalize()
            l.append(p)
        theta += resolution

    # assign uniform probability distribution
    totSolidAngle = sum([p.solid_angle for p in l])
    for i,p in enumerate(l):
        l[i].probability = p.solid_angle/totSolidAngle

    return l

def TwoSiteSkyGrid(ifos, gpstime, dt=0.0005, sky='half', point=None):
    """
    Generates a SkyPositionTable for a two-site all-sky grid, covering either
    a hemisphere (sky='half') or full sphere (sky='full').

    The grid can be forced to pass through the given SkyPosition point if
    required.
    """

    # remove duplicate site references
    ifos = parse_sites(ifos)
    assert len(ifos)==2, "More than two sites in given list of detectors"

    # get detectors
    detectors = [inject.cached_detector.get(inject.prefix_to_name[ifo])\
                 for ifo in ifos]

    # get light travel time
    ltt = inject.light_travel_time(ifos[0], ifos[1])

    # get number of points
    max = int(np.floor(ltt/dt))

    # get baseline
    baseline = detectors[1].location - detectors[0].location
    baseline /= np.linalg.norm(baseline)

    #
    # construct rotation
    #

    # xaxis points overhead of detector 1 or given point
    # zaxis is baseline, or something else
    if point:
        xaxis = SphericalToCartesian(point)
        xaxis /= np.linalg.norm(xaxis)
        zaxis = np.cross(xaxis, baseline)
        zaxis /= np.linalg.norm(zaxis)
    else:
        xaxis = detectors[0].location
        xaxis /= np.linalg.norm(xaxis)
        zaxis = baseline
        yaxis = np.cross(zaxis, xaxis)
        yaxis /= np.linalg.norm(yaxis)

    yaxis = np.cross(xaxis, zaxis)
    yaxis /= np.linalg.norm(yaxis)

    # construct Euler rotation
    lineofnodes = np.cross([0,0,1], zaxis)
    lineofnodes /= np.linalg.norm(lineofnodes)

    # work out Euler angles
    alpha = np.arccos(np.dot([1,0,0], lineofnodes))
    beta  = np.arccos(np.dot([0,0,1], zaxis))
    gamma = np.arccos(np.dot(lineofnodes, xaxis))

    # get rotation
    R = _rotation_euler(alpha, beta, gamma)

    # construct sky table
    l = SkyPositionTable()
    if sky=='half': longs = [0]
    elif sky=='full': longs = [0,lal.PI]
    else: raise AttributeError("Sky type \"%s\" not recognised, please use "
                               "'half' or 'full'" % sky)

    for long in longs:
        for i in xrange(-max, max+1):
            # get time-delay
            t = i * dt
            # if time-delay greater that light travel time, skip
            if abs(t) >= ltt: continue
            # construct object
            e = SkyPosition()
            e.system = 'geographic'
            # project time-delay onto prime meridian
            e.longitude   = long
            e.latitude    = lal.PI_2-np.arccos(-t/ltt)
            e.probability = None
            e.solid_angle = None
            e.normalize()
            e = e.rotate(R)
            e = GeographicToEquatorial(e, LIGOTimeGPS(gpstime))
            if e not in l:
                l.append(e)

    return l

def ThreeSiteSkyGrid(ifos, gpstime, dt=0.0005, tiling='square', sky='half'):
    """
    Generates a SkyPositionTable for a three-site all-sky grid, covering either
    a hemisphere (sky='half') or full sphere (sky='full'), using either
    'square' or 'hexagonal tiling.
    """

    # remove duplicate site references
    pop = []
    for i,ifo in enumerate(ifos):
        if i==0:  continue
        site = ifo[0]
        if site in [x[0] for x in ifos[0:i]]:
            sys.stderr.write("Duplicate site reference, ignoring %s.\n" % ifo)
            pop.append(i)
    for p in pop[::-1]:
        ifos.pop(p)

    assert len(ifos)==3

    detectors  = []
    T          = []
    baseline   = []

    # load detectors
    for i,ifo in enumerate(ifos):
        detectors.append(inject.cached_detector.get(inject.prefix_to_name[ifo]))
        # get light travel time and baseline for other detectors to first
        if i>=1:
            T.append(inject.light_travel_time(ifos[0], ifos[i]))
            baseline.append(detectors[i].location - detectors[0].location)
            baseline[i-1] /= np.linalg.norm(baseline[i-1])

    # get angle between baselines
    angle = np.arccos(np.dot(baseline[0], baseline[1]))
    
    #
    # construct tiling vectors
    #

    tiling = tiling.lower()
    t  = np.zeros(len(baseline))
    tdgrid = []

    # square grid
    dx = np.array([1,0])
    dy = np.array([0,1])
    nx = int(np.floor(T[0]/dt))
    ny = int(np.floor(T[1]/dt))

    # hexagonal grid
    if tiling == 'square':
        dm = dx
        dn = dy
        nm = nx
        nn = ny
        x = np.linspace(-nx, nx, 2*nx+1)
        y = np.linspace(-ny, ny, 2*ny+1)
        grid = np.array([[i,j] for i in x for j in y])
    # hexagonal grid
    elif re.match('hex', tiling):
        dm = (2*dx+dy)
        dn = (-dx+dy)
        dx = dm - dn
        # tile grid so that origin gets tiled
        grid = []
        orig = np.array([0,0])

        # find bottom left point on grid
        while orig[1] >= -ny: orig -= dn

        # loop over y
        for y in xrange(-ny,ny+1):
            orig[0] = orig[0] % dx[0]
            # work out minimal x value
            minx = -nx
            while minx % 3 != orig[0]: minx+=1
            # tile x
            x = np.arange(minx, nx+1, dx[0])
            # append to grid
            grid.extend([[i,y] for i in x])
            orig += dn
            orig[0] = orig[0] % dx[0]
        grid = np.asarray(grid)

    #
    # Following Rabaste, Chassande-Motin, Piran (2009), construct an ellipse
    # of physical values in 2D time-delay space for 3 detectors
    #

    #
    # (theta, phi) coordinate system is as follows:
    # detector 1 is origin
    # z-axis joins positively to detector 2
    # x-axis is orthogonal to z-axis in plane joining D1, D2, and D3
    # y-axis forms a right-handed coordinate system
    # 

    # so need to rotate from Rabaste network coordinates into
    # earth-fixed coordinates at the end
    
    # set z-axis in Rabaste frame
    zaxis = baseline[0]

    # work out y-axis in Rabaste frame
    yaxis = np.cross(zaxis, baseline[1])
    yaxis /= np.linalg.norm(yaxis)

    # work out x-axis in Rabaste frame
    xaxis = np.cross(yaxis, zaxis)
    xaxis /= np.linalg.norm(xaxis)

    # work out lines of nodes
    north  = np.array([0, 0, 1])
    lineofnodes = np.cross(baseline[0], north)
    lineofnodes /= np.linalg.norm(lineofnodes)

    # work out Euler angles
    alpha = np.arccos(np.dot(xaxis, lineofnodes))
    beta  = np.arccos(np.dot(baseline[0], north))
    gamma = np.arccos(np.dot(lineofnodes, [1,0,0]))

    # construct rotation
    R     = _rotation_euler(alpha, beta, gamma)

    #
    # construct sky position table
    #

    l = SkyPositionTable()
    # loop over time-delay between D1 and D2
    k = 0
    for i in grid:
        t = dt*i
        # construct coefficients of elliptical equation in quadratic form
        A      = -T[1]/T[0] * t[0] * np.cos(angle)
        B      = T[1]**2 * ((t[0]/T[0])**2 - np.sin(angle)**2)
        # test condition for inside ellipse and place point
        condition = t[1]**2 + 2*A*t[1] + B
        if condition <= 0:
            ntheta = np.arccos(-t[0]/T[0])
            nphi   = np.arccos(-(T[0]*t[1]-T[1]*t[0]*np.cos(angle))/\
                                (T[1]*np.sqrt(T[0]**2-t[0]**2)*np.sin(angle)))
            p = SkyPosition()
            p.system = 'geographic'
            p.longitude   = nphi
            p.latitude    = lal.PI_2 - ntheta
            p.probability = None
            p.solid_angle = None
            p.normalize()
            p = p.rotate(R)
            p = GeographicToEquatorial(p, LIGOTimeGPS(gpstime))
            l.append(p)
            if sky=='full':
                p2 = SkyPosition()
                p2.longitude = p.longitude
                p2.latitude = p.latitude + lal.PI
                p2.system = 'equatorial'
                p2.normalize()
                l.append(p2)
 
    return l

def ISOTimeDelayLine(ifos, ra, dec, gpstime=None, n=3):
    """
    Returns the n-point SkyPositionTable describing a line of constant time
    delay through the given ra and dec. for the given 2-tuple ifos.
    """

    if gpstime:
        gpstime = LIGOTimeGPS(gpstime)

    p = SkyPosition()
    p.longitude = ra
    p.latitude  = dec
    p.system    = 'equatorial'
    p.probability = None
    p.solid_angle = None
    if gpstime:
        p = EquatorialToGeographic(p, gpstime)
    cart = SphericalToCartesian(p)

    # get baseline
    detectors = [inject.cached_detector.get(inject.prefix_to_name[ifo])\
                 for ifo in ifos]
    baseline = detectors[0].location - detectors[1].location
    baseline = baseline / np.linalg.norm(baseline)

    # get angle to baseline
    angle = np.arccos(np.dot(cart, baseline))

    # get evenly spaced ring over north pole
    longitudes = np.linspace(0, lal.TWOPI, n, endpoint=False)
    latitudes  = [lal.PI_2-angle]*len(longitudes)
    # get axis and angle of rotation
    north = np.array([0,0,1])
    axis  = np.cross(north, baseline)
    axis  = axis / np.linalg.norm(axis)
    angle = np.arccos(np.dot(north, baseline))
    R     = _rotation(axis, angle)

    # construct sky table
    iso = SkyPositionTable()
    for lon,lat in zip(longitudes, latitudes):
        e             = SkyPosition()
        e.longitude   = lon
        e.latitude    = lat
        e.probability = None
        e.solid_angle = None
        e.system      = 'geographic'
        e.normalize()
        e             = e.rotate(R)
        if gpstime:
            e         = GeographicToEquatorial(e, gpstime)
        iso.append(e)

    return iso

def MaxTimeDelayLine(ifos, ra, dec, gpstime=None, n=3):
    """
    Return the n-point SkyPositionTable describing the line perpendicular to
    the line of constant time delay through the given ra and dec for the
    2-tuple ifos.
    """

    # construct sky point
    p = SkyPosition()
    p.longitude = ra
    p.latitude  = dec
    p.system    = 'equatorial'

    # remove duplicate site references
    ifos = parse_sites(ifos)
    assert len(ifos)==2, "More than two sites in given list of detectors"

    # get light travel time
    ltt = inject.light_travel_time(ifos[0], ifos[1])
    # get number of points
    dt = ltt/n

    return TwoSiteSkyGrid(ifos, gpstime, dt=dt, point=p, sky='full')

def MaxTimeDelayLine3(ifo1, ifo2, ra, dec, gpstime=None, n=3):
    """
        Alternative implementation of MaxTimeDelayLine.
    """

    if gpstime:
        gpstime = LIGOTimeGPS(gpstime)

    p = SkyPosition()
    p.longitude = ra
    p.latitude  = dec
    p.system    = 'equatorial'
    p.probability = None
    p.solid_angle = None
    if gpstime:
        p = EquatorialToGeographic(p, gpstime)
    cart = SphericalToCartesian(p)

    # get baseline
    detectors = [inject.cached_detector.get(inject.prefix_to_name[ifo])\
                 for ifo in [ifo1,ifo2]]
    baseline = detectors[0].location - detectors[1].location
    baseline = baseline / np.linalg.norm(baseline)

    # get angular spacing
    dtheta = lal.TWOPI/n

    # get axis and angle of rotation
    north = np.array([0,0,1])
    axis  = np.cross(cart, baseline)
    axis  = axis / np.linalg.norm(axis)
    R     = _rotation(axis, dtheta)

    # set up list
    l = SkyPositionTable()
    # append central point
    l.append(p)

    for i in xrange(1, n):
        l.append(l[i-1].rotate(R))

    if gpstime:
        for i,p in enumerate(l):
            l[i] = GeographicToEquatorial(p, gpstime)
    return l
            
# =============================================================================
# Helpers
# =============================================================================

def SphericalToCartesian(skyPoint):
    """
    Convert SkyPosition object into Cartesian 3-vector
    """ 
 
    p     = np.zeros(3)
    phi   = skyPoint.longitude
    theta = lal.PI_2 - skyPoint.latitude
    a     = np.sin(phi)
    b     = np.cos(phi)
    c     = np.sin(theta)
    d     = np.cos(theta)

    p[0] = c*b
    p[1] = c*a
    p[2] = d

    return p 

def CartesianToSpherical(x, system='equatorial'):
    """
    Convert 3-tuple Cartesian sky position x to SkyPosition object in the
    given system
    """

    assert len(x)==3

    p = SkyPosition()
    p.system = system
    p.longitude   = np.arctan2(x[1], x[0])
    p.latitude    = lal.PI_2 - np.arccos(x[2])
    p.probability = None
    p.solid_angle = None
    p.normalize()

    return p

def _rotation(axis, angle):
    """
    Form 3x3 rotation matrix to rotate about a given 3-tuple axis by a given
    angle
    """

    R = np.zeros((3,3))

    ux,uy,uz = axis
    c        = np.cos(angle)
    s        = np.sin(angle)

    R[0][0] = c + ux**2*(1-c)
    R[0][1] = ux*uy*(1-c) - uz*s
    R[0][2] = ux*uz*(1-c) + uy*s
    
    R[1][0] = uy*ux*(1-c) + uz*s
    R[1][1] = c + uy**2*(1-c)
    R[1][2] = uy*uz*(1-c) - ux*s

    R[2][0] = uz*ux*(1-c) - uy*s
    R[2][1] = uz*uy*(1-c) + ux*s
    R[2][2] = c + uz**2*(1-c)

    return R

def _rotation_euler(alpha, beta, gamma):
    """
    Form rotation matrix from Euler angles.
    """
    c = np.cos
    s = np.sin

    # Define rotation matrix
    R = np.array([[c(alpha) * c(gamma) + s(alpha) * s(beta) * s(gamma),
                   c(beta) * s(alpha),
                   c(gamma) * s(alpha) * s(beta) - c(alpha) * s(gamma)],
                  [c(alpha) * s(beta) * s(gamma) - c(gamma) * s(alpha),
                   c(alpha) * c(beta),
                   s(alpha) * s(gamma) + c(alpha) * c(gamma) * s(beta)],
                  [c(beta) * s(gamma),
                   -s(beta),
                   c(beta) * c(gamma)]],
                  dtype=float)

    return R

def parse_sites(ifos):
    """
    Returns a new list of interferometers containing one only per site.
    I.e. this removes H2 if included.
    """

    # rebind ifos
    ifos = list(ifos)
    # remove duplicate site references
    ifos2 = []
    for ifo in ifos:
        if len([i for i in ifos2 if i.startswith(ifo[0])]) == 0:
            ifos2.append(ifo)

    return ifos2
