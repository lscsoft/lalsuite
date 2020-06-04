# Copyright (C) 2019 Matthew Pitkin
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
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

## \defgroup lalpulsar_py_simulateHeterodynedCW SimulateHeterodynedCW
## \ingroup lalpulsar_python
"""
The module provides the HeterodynedCWSimulator() class for simulating a
signal from a continuous wave source after application of a heterodyned as
described in Equations 7 and 8 of @cite Pitkin2017 . An example usage to
generate the complex heterodyned signal time series is:

~~~
from lalpulsar.simulateHeterodynedCW import HeterodynedCWSimulator
from lalpulsar.PulsarParametersWrapper import PulsarParametersPy
import lal
import numpy as np

# set the pulsar parameters
par = PulsarParametersPy()
par['RAJ'] = lal.TranslateHMStoRAD('01:23:34.5')
par['DECJ'] = lal.TranslateDMStoRAD('-45:01:23.4')
par['F'] = [123.456789, -9.87654321e-12]  # frequency and first derivative
pepoch = lal.TranslateStringMJDTTtoGPS('58000')   # frequency epoch
par['PEPOCH'] = pepoch.gpsSeconds + 1e-9*pepoch.gpsNanoSeconds
par['H0'] = 5.6e-26     # GW amplitude
par['COSIOTA'] = -0.2   # cosine of inclination angle
par['PSI'] = 0.4        # polarization angle (rads)
par['PHI0'] = 2.3       # initial phase (rads)

# set the GPS times of the data
times = np.arange(1000000000.0, 1000086400., 3600)

# set the detector
det = 'H1'  # the LIGO Hanford Observatory

# create the HeterodynedCWSimulator object
het = HeterodynedCWSimulator(par, det, times=times)

# get the model complex strain time series
model = het.model(usephase=True)
~~~

An example of getting the time series for a signal that has phase parameters
that are not identical to the heterodyned parameters would be:

~~~
from lalpulsar.simulateHeterodynedCW import HeterodynedCWSimulator
from lalpulsar.PulsarParametersWrapper import PulsarParametersPy
import lal
import numpy as np

# set the "heterodyne" pulsar parameters
par = PulsarParametersPy()
par['RAJ'] = lal.TranslateHMStoRAD('01:23:34.6')
par['DECJ'] = lal.TranslateDMStoRAD('-45:01:23.5')
par['F'] = [123.4567, -9.876e-12]  # frequency and first derivative
pepoch = lal.TranslateStringMJDTTtoGPS('58000')   # frequency epoch
par['PEPOCH'] = pepoch.gpsSeconds + 1e-9*pepoch.gpsNanoSeconds

# set the times
times = np.arange(1000000000., 1000000600., 60)

# set the detector
det = 'H1'  # the LIGO Hanford Observatory

# create the HeterodynedCWSimulator object
het = HeterodynedCWSimulator(par, det, times=times)

# set the updated parameters
parupdate = PulsarParametersPy()
parupdate['RAJ'] = lal.TranslateHMStoRAD('01:23:34.5')
parupdate['DECJ'] = lal.TranslateDMStoRAD('-45:01:23.4')
parupdate['F'] = [123.456789, -9.87654321e-12]  # frequency and first derivative
pepoch = lal.TranslateStringMJDTTtoGPS('58000')   # frequency epoch
parupdate['PEPOCH'] = pepoch.gpsSeconds + 1e-9*pepoch.gpsNanoSeconds
parupdate['H0'] = 5.6e-26     # GW amplitude
parupdate['COSIOTA'] = -0.2   # cosine of inclination angle
parupdate['PSI'] = 0.4        # polarization angle (rads)
parupdate['PHI0'] = 2.3       # initial phase (rads)

# get the model complex strain time series
model = het.model(parupdate, usephase=True, updateSSB=True)
~~~

"""
## @{

from __future__ import (division, print_function)

import numpy as np
from six import string_types

try:
    import lal
except ImportError:
    raise ImportError("SWIG wrappings of LAL cannot be imported")

try:
    import lalpulsar
except ImportError:
    raise ImportError("SWIG wrappings of LALPulsar cannot be imported")

try:
    from .PulsarParametersWrapper import PulsarParametersPy
except ImportError:
    raise ImportError("Cannot import PulsarParametersPy class")

from . import git_version


__author__ = "Matthew Pitkin <matthew.pitkin@ligo.org>"
__version__ = git_version.id
__date__ = git_version.date


DOWNLOAD_URL = 'https://git.ligo.org/lscsoft/lalsuite/raw/master/lalpulsar/lib/{}'

class HeterodynedCWSimulator(object):

    def __init__(self, par, det, times=None, earth_ephem=None,
                 sun_ephem=None, time_corr=None, ephem='DE405', units='TCB',
                 t0=None, dt=None):
        """
        A class to simulate strain data for a continuous gravitational-wave
        signal after the data has been heterodyned, i.e., after multiplying
        the data by a complex phase vector. This uses the Equations 7 and 8
        from @cite Pitkin2017 accessed via the XLALHeterodynedPulsarGetModel()
        function.

        @param par: a TEMPO-style text file, or a PulsarParametersPy()
            structure, containing the parameters for the source, in particular
            the phase parameters at which the data is "heterodyned".
        @param det: the name of a detector.
        @param times: an array of GPS times at which to calculate the
            heterodyned strain.
        @param t0: a time epoch in GPS seconds at which to calculate the
            detector response function. If not given and @b times is set,
            then the first value of @b times will be used.
        @param dt: the time steps (in seconds) in the data over which to
            average the detector response. If not given and @b times is set,
            the the time difference between the first two values in @b times
            will be used.
        @param earth_ephem: a file containing the Earth ephemeris information.
            If not set then a default file will be used.
        @param sun_ephem: a file containing the Earth ephemeris information.
            If not set then a default file will be used.
        @param time_corr: a file containing information on the time system
            corrections for, e.g., the TCB or TDB system. If not set then
            a default file will be used.
        @param ephem: The solar system ephemeris system to use for the Earth
            and Sun ephemeris, i.e., @c 'DE200', @c 'DE405', @c 'DE421', or
            @c 'DE430'. By default the @c 'EPHEM' value from the supplied
            @b par will be used, but if not found, and if this value is not
            set, it will default to @c 'DE405'.
        @param units: The time system used, i.e., @c 'TDB' or @c 'TCB'. By default
            the @c 'UNITS' value from the @b par will be used, but if not
            found, and if this value is not set, it will (like TEMPO2) default
            to @c 'TCB'.
        """

        self.__hetpar = self._read_par(par)
        self._set_detector(det)
        self.times = times

        # set default ephemeris strings
        self.__earthstr = 'earth00-40-{}.dat.gz'
        self.__sunstr = 'sun00-40-{}.dat.gz'
        self.__timecorrstr = '{}_2000-2040.dat.gz'

        # mapping between time units and time correction file prefix
        self.__units_map = {'TCB': 'te405',
                            'TDB': 'tdb'}

        self.ephem = ephem
        self.units = units

        # initialise the solar system ephemeris files
        self.__edat, self.__tdat = self._initialise_ephemeris(earth_ephem,
                                                              sun_ephem,
                                                              time_corr)

        # set the "heterodyne" SSB time delay
        if self.times is not None:
            self.__hetSSBdelay = lalpulsar.HeterodynedPulsarGetSSBDelay(self.hetpar.PulsarParameters(),
                                                                        self.gpstimes,
                                                                        self.detector,
                                                                        self.__edat,
                                                                        self.__tdat,
                                                                        self.__units_type)
        else:
            self.__hetSSBdelay = None

        # set the "heterodyne" BSB time delay
        if self.times is not None and self.hetpar['BINARY'] is not None:
            self.__hetBSBdelay = lalpulsar.HeterodynedPulsarGetBSBDelay(self.hetpar.PulsarParameters(),
                                                                        self.gpstimes,
                                                                        self.__hetSSBdelay,
                                                                        self.__edat)
        else:
            self.__hetBSBdelay = None

        # set the response function
        if self.times is None and t0 is None:
            raise ValueError("Must supply either 'times' or 't0' to calculate "
                             "the response function")
        else:
            self.__t0 = t0 if t0 is not None else self.times[0]

        if dt is None and self.times is None:
            raise ValueError("Must supply either 'times' or 'dt' to calculate "
                             "the response function")
        else:
            if self.times is not None and dt is None:
                if len(self.times) == 1:
                    raise ValueError("Must supply a 'dt' value")
                else:
                    self.__dt = self.times[1] - self.times[0]
            else:
                self.__dt = dt

        ra = self.hetpar['RA'] if self.hetpar['RAJ'] is None else self.hetpar['RAJ']
        dec = self.hetpar['DEC'] if self.hetpar['DECJ'] is None else self.hetpar['DECJ']
        if ra is None or dec is None:
            raise ValueError("Right ascension and/or declination have not "
                             "been set!")

        self.__resp = lalpulsar.DetResponseLookupTable(self.__t0,
                                                       self.detector,
                                                       ra,
                                                       dec,
                                                       2880,
                                                       self.__dt)

    @property
    def hetpar(self):
        return self.__hetpar

    def _set_detector(self, det):
        if isinstance(det, lal.Detector):
            # value is already a lal.Detector
            self.__detector = det
        else:
            if not isinstance(det, string_types):
                raise TypeError('Detector name must be a string')
            else:
                try:
                    self.__detector = lalpulsar.GetSiteInfo(det)
                except RuntimeError:
                    raise ValueError("Detector '{}' was not a valid detector "
                                     "name.".format(det))

        self.__detector_name = self.__detector.frDetector.name

    @property
    def detector(self):
        return self.__detector

    @property
    def resp(self):
        """
        Return the response function look-up table.
        """

        return self.__resp

    @property
    def times(self):
        return self.__times

    @property
    def gpstimes(self):
        return self.__gpstimes

    @times.setter
    def times(self, times):
        """
        Set an array of times, and also a @b LIGOTimeGPSVector() containing the
        times.
        """

        if times is None:
            self.__times = None
            self.__gpstimes = None
            return
        elif isinstance(times, lal.LIGOTimeGPS):
            self.__times = np.array([times.gpsSeconds + 1e-9*times.gpsNanoSeconds], dtype='float64')
            self.__gpstimes = lalpulsar.CreateTimestampVector(1)
            self.__gpstimes.data[0] = times
            return
        elif isinstance(times, lalpulsar.LIGOTimeGPSVector):
            self.__gpstimes = times
            self.__times = np.zeros(len(times.data), dtype='float64')
            for i, gpstime in enumerate(times.data):
                self.__times[i] = times.data[i].gpsSeconds + 1e-9*times.data[i].gpsNanoSeconds
            return
        elif isinstance(times, (int, float)):
            self.__times = np.array([times], dtype='float64')
        elif isinstance(times, (list, np.ndarray)):
            self.__times = np.array(times, dtype='float64')
        else:
            raise TypeError("Unknown data type for times")

        self.__gpstimes = lalpulsar.CreateTimestampVector(len(self.__times))
        for i, time in enumerate(self.__times):
            self.__gpstimes.data[i] = lal.LIGOTimeGPS(time)

    @property
    def ephem(self):
        return self.__ephem

    @ephem.setter
    def ephem(self, ephem):
        """
        Set the heterodyne solar system ephemeris version. This will attempt to
        use the value set in the heterodyne source parameters, but otherwise
        defaults to DE405.
        """

        if self.hetpar['EPHEM'] is not None:
            self.__ephem = self.hetpar['EPHEM']
        else:
            self.__ephem = 'DE405' if ephem is None else ephem

    @property
    def units(self):
        return self.__units

    @units.setter
    def units(self, units):
        """
        Set the time system units, i.e., either 'TDB' or 'TCB'. This will
        attempt to use the value set in the heterodyne source parameters, but
        otherwise defaults to 'TCB'.
        """

        if self.hetpar['UNITS'] is not None:
            self.__units = self.hetpar['UNITS']
        else:
            self.__units = 'TCB' if units is None else units

        if self.__units not in ['TCB', 'TDB']:
            raise ValueError("Unknown time system '{}' has been "
                             "given.".format(self.__units))

        if self.__units == 'TCB':
            self.__units_type = lalpulsar.lalpulsar.TIMECORRECTION_TCB
        else:
            self.__units_type = lalpulsar.lalpulsar.TIMECORRECTION_TDB

    def model(self, newpar=None, updateSSB=False, updateBSB=False,
              freqfactor=2., usephase=False, roq=False):
        """
        Compute the heterodyned strain model using
        XLALHeterodynedPulsarGetModel().

        @param newpar: A text parameter file, or PulsarParameterPy() object,
            containing a set of parameter at which to calculate the strain
            model. If this is @c None then the "heterodyne" parameters are used.
        @param updateSSB: set to @c True to update the solar system barycentring
            time delays compared to those used in heterodyning, i.e., if the
            @b newpar contains updated positional parameters.
        @param updateBSB: set to @c True to update the binary system barycentring
            time delays compared to those used in heterodying, i.e., if the
            @b newpar contains updated binary system parameters
        @param freqfactor: the factor by which the frequency evolution is
            multiplied for the source model. This defaults to 2 for emission
            from the \f$l=m=2\f$ quadrupole mode.
        @param usephase: set to @c True if the model is to include the phase
            evolution, i.e., if phase parameters are being updated, otherwise
            only two (six for non-GR sources) values giving the amplitides
            will be output.
        @param roq: a boolean value to set to @c True requiring the output for
            a ROQ model.

        @return a complex array called @b compstrain
        """

        if newpar is not None:
            parupdate = self._read_par(newpar)
        else:
            parupdate = self.hetpar

        # get frequency differences
        if usephase:
            if parupdate['DELTAF'] is None:
                try:
                    parupdate['DELTAF'] = parupdate['F'] - self.hetpar['F']
                except Exception as e:
                    raise ValueError("Frequencies are not set in parameter "
                                     "objects: {}".format(e))

        self.__nonGR = self._check_nonGR(parupdate)
        compstrain = lalpulsar.HeterodynedPulsarGetModel(parupdate.PulsarParameters(),
                                                         freqfactor,
                                                         int(usephase),  # phase is varying between par files
                                                         int(roq),       # using ROQ?
                                                         self.__nonGR,   # using non-tensorial modes?
                                                         self.gpstimes,
                                                         self.__hetSSBdelay,
                                                         int(updateSSB),  # the SSB delay should be updated compared to hetSSBdelay
                                                         self.__hetBSBdelay,
                                                         int(updateBSB),  # the BSB delay should be updated compared to hetBSBdelay
                                                         self.resp,
                                                         self.__edat,
                                                         self.__tdat,
                                                         self.__units_type)

        return compstrain.data.data

    def _read_par(self, par):
        """
        Read a TEMPO-style parameter file into a PulsarParameterPy object.
        """

        if isinstance(par, PulsarParametersPy):
            return par

        if isinstance(par, string_types):
            try:
                return PulsarParametersPy(par)
            except IOError:
                raise IOError("Could not read in parameter file: '{}'".format(par))
        else:
            raise TypeError("The parameter file must be a string")

    def _check_nonGR(self, par):
        """
        Check if the source parameters are for a non-GR model, i.e., are any of
        the amplitude/phase parameters for a non-GR model set
        """

        # non-GR amplitude parameters
        nonGRparams = ['HPLUS',
                       'HCROSS',
                       'HVECTORX',
                       'HVECTORY',
                       'HSCALARB',
                       'HSCALARL',
                       'HPLUS_F',
                       'HCROSS_F',
                       'HVECTORX_F',
                       'HVECTORY_F',
                       'HSCALARB_F',
                       'HSCALARL_F']

        for param in nonGRparams:
            if param in par.keys():
                return 1

        return 0

    def _initialise_ephemeris(self, earth_ephem, sun_ephem, time_corr):
        """
        Initialise the solar system ephemeris.
        """

        if earth_ephem is not None:
            earthfile = earth_ephem
        else:
            earthfile = self.__earthstr.format(self.ephem)

        if sun_ephem is not None:
            sunfile = sun_ephem
        else:
            sunfile = self.__sunstr.format(self.ephem)

        if time_corr is not None:
            timefile = time_corr
        else:
            timefile = self.__timecorrstr.format(self.__units_map[self.units])

        try:
            edat = lalpulsar.InitBarycenter(earthfile, sunfile)
        except RuntimeError:
            try:
                # try downloading the ephemeris files
                from astropy.utils.data import download_file

                efile = download_file(DOWNLOAD_URL.format(earthfile), cache=True)
                sfile = download_file(DOWNLOAD_URL.format(sunfile), cache=True)
                edat = lalpulsar.InitBarycenter(efile, sfile)
            except Exception as e:
                raise IOError("Could not read in ephemeris files: {}".format(e))

        try:
            tdat = lalpulsar.InitTimeCorrections(timefile)
        except RuntimeError:
            try:
                # try downloading the ephemeris files
                from astropy.utils.data import download_file

                tfile = download_file(DOWNLOAD_URL.format(timefile), cache=True)
                tdat = lalpulsar.InitTimeCorrections(tfile)
            except Exception as e:
                raise IOError("Could not read in time correction file: {}".format(e))

        return edat, tdat

## @}
