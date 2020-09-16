# Copyright (C) 2017 Karl Wette
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

## \defgroup lalpulsar_py_simulateCW SimulateCW
## \ingroup lalpulsar_python
"""
Generate strain time series of a continuous-wave signal in the detector frame,
given a function which computes signal phase and amplitudes as functions of time.

The given function should compute quantities \f$\Delta\phi(t)\f$, \f$a_+(t)\f$,
and \f$a_\times(t)\f$ such that the plus and cross gravitational-wave
polarisations \f$h_+(t)\f$ and \f$h_\times(t)\f$ are given by:
\f{eqnarray}{
h_+(t)      & = & a_+(t)     \cos[\phi_0 + \Delta\phi(t)] \; , \\
h_\times(t) & = & a_\times(t)\sin[\phi_0 + \Delta\phi(t)] \; .
\f}

This module provides a class CWSimulator() to generate the strain time series.
CWSimulator() can also directly write frame files and SFT files. Example usage:

~~~
import lal
from lalpulsar import simulateCW

def waveform(h0, cosi, freq, f1dot):
    def wf(dt):
        dphi = lal.TWOPI * (freq * dt + f1dot * 0.5 * dt**2)
        ap = h0 * (1.0 + cosi**2) / 2.0
        ax = h0 * cosi
        return dphi, ap, ax
    return wf

tref     = 900043200
tstart   = 900000000
Tdata    = 86400
h0       = 1e-24
cosi     = 0.123
psi      = 2.345
phi0     = 3.210
freq     = 10.0
f1dot    = -1.35e-8
dt_wf    = 5
alpha    = 6.12
delta    = 1.02
detector = 'H1'

wf = waveform(h0, cosi, freq, f1dot)
S = simulateCW.CWSimulator(tref, tstart, Tdata, wf, dt_wf, phi0, psi, alpha, delta, detector)

# To write SFT files
for file, i, N in S.write_sft_files(fmax=32, Tsft=1800, comment="simCW"):
    print('Generated SFT file %s (%i of %i)' % (file, i+1, N))

# To write frame files
for file, i, N in S.write_frame_files(fs=1, Tframe=1800, comment="simCW"):
    print('Generated frame file %s (%i of %i)' % (file, i+1, N))

~~~
"""
## @{

from __future__ import (division, print_function)

import os
import re
import math

import lal
import lalpulsar

from . import git_version

__author__ = "Karl Wette <karl.wette@ligo.org>"
__version__ = git_version.id
__date__ = git_version.date


class CWSimulator(object):

    def __init__(self, tref, tstart, Tdata, waveform, dt_wf, phi0, psi, alpha, delta, det_name,
                 earth_ephem_file="earth00-40-DE405.dat.gz", sun_ephem_file="sun00-40-DE405.dat.gz", tref_at_det=False, extra_comment=None):
        """
        Initialise a continuous-wave signal simulator.

        @param tref: reference time of signal phase at Solar System barycentre, in GPS seconds
            (but see @b tref_at_det)
        @param tstart: start time of signal, in GPS seconds
        @param Tdata: total duration of signal, in seconds
        @param waveform: function which computes signal phase and amplitudes as functions of time:
            @b dphi, @b aplus, @b across = @b waveform(dt), where:
                @b dt = time since reference time @b tref;
                @b dphi = phase of signal at time @b dt relative to reference time @b tref, in radians;
                @b aplus = strain amplitude of plus polarisation at time @b dt;
                @b across = strain amplitude of cross polarisation at time @b dt
        @param dt_wf: sampling time of the function @c waveform; this need only be small enough to ensure
            that @c dphi, @c aplus, and @c across are smoothly interpolated, and does not need to make
            the desired sampling frequency of the output strain time series
        @param phi0: initial phase of the gravitational-wave signal at @b tstart, in radians
        @param psi: polarisation angle of the gravitational-wave source, in radians
        @param alpha: right ascension of the gravitational-wave source, in radians
        @param delta: declination  of the gravitational-wave source, in radians
        @param det_name: name of the gravitational-wave detector to simulate a response for;
            e.g. @c "H1" for LIGO Hanford, @c "L1" for LIGO Livingston, @c "V1" for Virgo
        @param earth_ephem_file: name of file to load Earth ephemeris from
        @param sun_ephem_file: name of file to load Sun ephemeris from
        @param tref_at_det: default False; if True, shift reference time @b tref so that @b dt = 0 is
            @b tref in @e detector frame instead of Solar System barycentre frame, useful if e.g. one
            wants to turn on signal only for @b dt > 0, one can use @b tref as the turn-on time
        @param extra_comment: additional text to add to comment string in frame/SFT headers
               (not filenames), e.g. for wrapper script commandlines
        """

        self.__origin_str = "Generated by %s, %s-%s (%s)" % (__file__, git_version.id, git_version.status, git_version.date)
        if extra_comment:
            self.__origin_str += ", "+extra_comment

        # store arguments
        self.__tstart = tstart
        self.__Tdata = Tdata

        # parse detector name
        try:
            _, self.__site_index = lalpulsar.FindCWDetector(det_name, True)
            assert(self.__site_index >= 0)
        except:
            raise ValueError("Invalid detector name det_name='%s'" % det_name)
        self.__site = lal.CachedDetectors[self.__site_index]

        # load Earth and Sun ephemerides
        self.__ephemerides = lalpulsar.InitBarycenter(earth_ephem_file, sun_ephem_file)

        if tref_at_det:

            # calculate barycentric delay at reference time
            bary_state = lalpulsar.EarthState()
            bary_input = lalpulsar.BarycenterInput()
            bary_emit = lalpulsar.EmissionTime()
            bary_input.tgps = tref
            bary_input.site = self.__site
            for i in range(0, 3):
                bary_input.site.location[i] /= lal.C_SI
            bary_input.alpha = alpha
            bary_input.delta = delta
            bary_input.dInv = 0.0
            lalpulsar.BarycenterEarth(bary_state, bary_input.tgps, self.__ephemerides)
            lalpulsar.Barycenter(bary_emit, bary_input, bary_state)

            # adjust reference time so that dt = 0 is tref in detector frame
            tref += bary_emit.deltaT

        # start signal time series 'Tpad_wf' before/after output strain time series
        # add sufficient padding to signal for maximum Doppler modulation time shifts, otherwise
        # lalpulsar.PulsarSimulateCoherentGW() will output zeros without complaint (unless
        # you run with LAL_DEBUG_LEVEL=warning)
        Tpad_wf = 2.0 * lal.AU_SI / lal.C_SI
        tstart_wf = tstart - Tpad_wf
        Nwf = int(math.ceil(float(Tdata + 2*Tpad_wf) / float(dt_wf)))

        # create REAL8TimeSeries to store signal phase
        self.__phi = lal.CreateREAL8TimeSeries("phi", tstart_wf, 0, dt_wf, lal.DimensionlessUnit, Nwf)

        # create REAL4TimeVectorSeries to store signal amplitudes
        # - LAL provides no creator function for this type, so must be done manually
        self.__a = lal.REAL4TimeVectorSeries()
        self.__a.name = "a+,ax"
        self.__a.epoch = tstart_wf
        self.__a.deltaT = dt_wf
        self.__a.f0 = 0
        self.__a.sampleUnits = lal.StrainUnit
        self.__a.data = lal.CreateREAL4VectorSequence(Nwf, 2)

        # call waveform() to fill time series of signal phase and amplitudes
        dt = float(tstart_wf - tref)
        for i in range(0, Nwf):
            dphi, aplus, across = waveform(dt)
            self.__phi.data.data[i] = phi0 + dphi
            self.__a.data.data[i][0] = aplus
            self.__a.data.data[i][1] = across
            dt += dt_wf

        # create and initialise PulsarCoherentGW struct
        self.__wf = lalpulsar.PulsarCoherentGW()
        self.__wf.position.system = lal.COORDINATESYSTEM_EQUATORIAL
        self.__wf.position.longitude = alpha
        self.__wf.position.latitude = delta
        self.__wf.psi = psi
        self.__wf.phi = self.__phi
        self.__wf.a = self.__a

        # create and initialise PulsarDetectorResponse struct
        self.__detector = lalpulsar.PulsarDetectorResponse()
        self.__detector.site = self.__site
        self.__detector.ephemerides = self.__ephemerides

    def _simulate_coherent_gw(self, h, noise_sqrt_Sh, noise_seed):

        # generate strain time series
        lalpulsar.PulsarSimulateCoherentGW(h, self.__wf, self.__detector)

        # add Gaussian noise, if requested
        if noise_sqrt_Sh > 0:
            noise_sigma = math.sqrt(0.5 / h.deltaT) * noise_sqrt_Sh
            lalpulsar.AddGaussianNoise(h, noise_sigma, noise_seed)

    def get_strain(self, fs, tmin=0, tmax=None, noise_sqrt_Sh=0, noise_seed=0):
        """
        Generate strain time series of a continuous-wave signal in the detector frame.

        @param fs: sampling frequency of strain time series, in Hz
        @param tmin: start time for strain time series, as offsets from self.__tstart
        @param tmax: start time for strain time series, as offsets from self.__tstart
        @param noise_sqrt_Sh: if >0, add Gaussian noise with square-root single-sided power
            spectral density given by this value, in Hz^(-1/2)
        @param noise_seed: use this need for the random number generator used to create noise

        @return (@b t, @b h), where:
            @b t = start of time strain time series, in GPS seconds;
            @b h = strain time series
        """

        # process tmin/tmax range (interpreted relative to self.__tstart)
        if (tmin < 0) or (tmin>=self.__Tdata):
            raise ValueError('tmin must be within [0,{}).'.format(self.__Tdata))
        if tmax is None:
            tmax = self.__Tdata
        elif (tmax <= 0) or (tmax>self.__Tdata):
            raise ValueError('tmax must be within (0,{}].'.format(self.__Tdata))
        tspan = tmax-tmin

        # create REAL4TimeSeries to store output time series
        Nh = int(fs*tspan)
        h = lal.CreateREAL4TimeSeries("h", self.__tstart+tmin, 0, 1.0 / fs,
                                      lal.DimensionlessUnit, Nh)

        # generate strain time series
        self._simulate_coherent_gw(h, noise_sqrt_Sh, noise_seed)
        return (h.epoch, h.data.data)

    def get_strain_blocks(self, fs, Tblock, noise_sqrt_Sh=0, noise_seed=0):
        """
        Generate strain time series of a continuous-wave signal in the detector frame, in contiguous blocks.

        @param fs: sampling frequency of strain time series, in Hz
        @param Tblock: length of each block, in seconds; should divide evenly into @b Tdata
        @param noise_sqrt_Sh: if >0, add Gaussian noise with square-root single-sided power
            spectral density given by this value, in Hz^(-1/2)
        @param noise_seed: use this need for the random number generator used to create noise

        @return (@b t, @b h, @b i, @b N), where:
            @b t = start of time strain time series, in GPS seconds;
            @b h = strain time series;
            @b i = block index, starting from zero;
            @b N = number of blocks

        This is a Python generator function and so should be called as follows:
        ~~~
        S = CWSimulator(...)
        for t, h, i, N in S.get_strain_blocks(...):
            ...
        ~~~
        """

        # work out number of blocks
        Nblock = int(round(self.__Tdata / Tblock))
        if Tblock * Nblock > self.__Tdata:
            raise ValueError("Length of block Tblock=%g does not divide evenly into Tdata=%g" % (Tblock, self.__Tdata))

        # generate strain time series in blocks of length 'Tblock'
        tmin = 0
        for iblock in range(0, Nblock):
            epoch, hoft = self.get_strain(fs, tmin=tmin, tmax=tmin+Tblock,
                                          noise_sqrt_Sh=noise_sqrt_Sh,
                                          noise_seed=noise_seed + iblock)
            yield epoch, hoft, iblock, Nblock
            tmin += Tblock

    def write_frame_files(self, fs, Tframe, comment, out_dir=".", noise_sqrt_Sh=0, noise_seed=0):
        """
        Write frame files [1] containing strain time series of a continuous-wave signal.

        The strain time series is written as double-precision post-processed data (ProcData) channel named
        <tt>&lt;detector&gt;:SIMCW-STRAIN</tt>,
        where <tt>&lt;detector&gt;</tt> is the 2-character detector prefix (e.g. <tt>H1</tt> for LIGO Hanford,
        <tt>L1</tt> for LIGO Livingston, <tt>V1</tt> for Virgo).

        @param fs: sampling frequency of strain time series, in Hz
        @param Tframe: length of each frame, in seconds; should divide evenly into @b Tdata
        @param comment: frame file name comment, may only contain A-Z, a-z, 0-9, _, +, # characters
        @param out_dir: output directory to write frame files into
        @param noise_sqrt_Sh: if >0, add Gaussian noise with square-root single-sided power
            spectral density given by this value, in Hz^(-1/2)
        @param noise_seed: use this need for the random number generator used to create noise

        @return (@b file, @b i, @b N), where:
            @b file = name of frame file just written;
            @b i = frame file index, starting from zero;
            @b N = number of frame files

        This is a Python generator function and so should be called as follows:
        ~~~
        S = CWSimulator(...)
        for file, i, N in S.write_frame_files(...):
            ...
        ~~~

        [1] https://dcc.ligo.org/LIGO-T970130/public
        """

        try:
            import lalframe
        except ImportError:
            raise ImportError("SWIG wrappings of LALFrame cannot be imported")

        # check for valid frame filename comment (see LIGO-T010150)
        valid_comment = re.compile(r'^[A-Za-z0-9_+#]+$')
        if not valid_comment.match(comment):
            raise ValueError("Frame file comment='%s' may only contain A-Z, a-z, 0-9, _, +, # characters" % comment)

        # generate strain time series in blocks of length 'Tframe'
        frame_h = None
        for t, h, i, N in self.get_strain_blocks(fs, Tframe, noise_sqrt_Sh=noise_sqrt_Sh, noise_seed=noise_seed):

            # create and initialise REAL8TimeSeries to write to frame files
            if frame_h is None:
                frame_h_channel = '%s:SIMCW-STRAIN' % self.__site.frDetector.prefix
                frame_h = lal.CreateREAL8TimeSeries(frame_h_channel, t, 0, 1.0 / fs, lal.DimensionlessUnit, len(h))
            frame_h.epoch = t
            frame_h.data.data = h

            # create standard frame file name (see LIGO-T010150)
            frame_src = self.__site.frDetector.prefix[0]
            frame_desc = 'simCW_%s' % comment
            frame_t0 = int(t.gpsSeconds)
            frame_Tdata = int(math.ceil(float(t + Tframe)) - math.floor(float(t)))
            frame_name = '%s-%s-%u-%u.gwf' % (frame_src, frame_desc, frame_t0, frame_Tdata)
            frame_path = os.path.join(out_dir, frame_name)

            # create frame
            frame_det_bits = 2 * self.__site_index
            frame = lalframe.FrameNew(t, Tframe, "simCW", -1, i, frame_det_bits)

            # add strain time series to frame
            lalframe.FrameAddREAL8TimeSeriesProcData(frame, frame_h)

            # add history
            lalframe.FrameAddFrHistory(frame, "origin", self.__origin_str)

            # write frame
            lalframe.FrameWrite(frame, frame_path)

            # yield current file name for e.g. printing progress
            yield frame_path, i, N

    def get_sfts(self, fmax, Tsft, noise_sqrt_Sh=0, noise_seed=0, window=None, window_param=0):
        """
        Generate SFTs [2] containing strain time series of a continuous-wave signal.

        @param fmax: maximum SFT frequency, in Hz
        @param Tsft: length of each SFT, in seconds; should divide evenly into @b Tdata
        @param noise_sqrt_Sh: if >0, add Gaussian noise with square-root single-sided power
            spectral density given by this value, in Hz^(-1/2)
        @param noise_seed: use this need for the random number generator used to create noise
        @param window: if not None, window the time series before performing the FFT, using
            the named window function; see XLALCreateNamedREAL8Window()
        @param window_param: parameter for the window function given by @b window, if needed

        @return (@b sft, @b i, @b N), where:
            @b sft = SFT;
            @b i = SFT file index, starting from zero;
            @b N = number of SFTs

        This is a Python generator function and so should be called as follows:
        ~~~
        S = CWSimulator(...)
        for sft, i, N in S.get_sfts(...):
            ...
        ~~~

        [2] https://dcc.ligo.org/LIGO-T040164/public
        """

        # create timestamps for generating one SFT per time series
        sft_ts = lalpulsar.CreateTimestampVector(1)
        sft_ts.deltaT = Tsft

        # generate strain time series in blocks of length 'Tsft'
        sft_h = None
        sft_fs = 2 * fmax
        for t, h, i, N in self.get_strain_blocks(sft_fs, Tsft, noise_sqrt_Sh=noise_sqrt_Sh, noise_seed=noise_seed):

            # create and initialise REAL8TimeSeries to write to SFT files
            if sft_h is None:
                sft_name = self.__site.frDetector.prefix
                sft_h = lal.CreateREAL8TimeSeries(sft_name, t, 0, 1.0 / sft_fs, lal.DimensionlessUnit, len(h))
            sft_h.epoch = t
            sft_h.data.data = h

            # create SFT, possibly with windowing
            sft_ts.data[0] = t
            sft_vect = lalpulsar.MakeSFTsFromREAL8TimeSeries(sft_h, sft_ts, window, window_param)

            # yield current SFT
            yield sft_vect.data[0], i, N

    def write_sft_files(self, fmax, Tsft, comment, out_dir=".", noise_sqrt_Sh=0, noise_seed=0, window=None, window_param=0):
        """
        Write SFT files [2] containing strain time series of a continuous-wave signal.

        @param fmax: maximum SFT frequency, in Hz
        @param Tsft: length of each SFT, in seconds; should divide evenly into @b Tdata
        @param comment: SFT file name comment, may only contain A-Z, a-z, 0-9, _, +, # characters
        @param out_dir: output directory to write SFT files into
        @param noise_sqrt_Sh: if >0, add Gaussian noise with square-root single-sided power
            spectral density given by this value, in Hz^(-1/2)
        @param noise_seed: use this need for the random number generator used to create noise
        @param window: if not None, window the time series before performing the FFT, using
            the named window function; see XLALCreateNamedREAL8Window()
        @param window_param: parameter for the window function given by @b window, if needed

        @return (@b file, @b i, @b N), where:
            @b file = name of SFT file just written;
            @b i = SFT file index, starting from zero;
            @b N = number of SFT files

        This is a Python generator function and so should be called as follows:
        ~~~
        S = CWSimulator(...)
        for file, i, N in S.write_sft_files(...):
            ...
        ~~~

        [2] https://dcc.ligo.org/LIGO-T040164/public
        """

        # check for valid SFT filename comment (see LIGO-T040164)
        valid_comment = re.compile(r'^[A-Za-z0-9_+#]+$')
        if not valid_comment.match(comment):
            raise ValueError("SFT file comment='%s' may only contain A-Z, a-z, 0-9, _, +, # characters" % comment)

        # generate SFTs
        for sft, i, N in self.get_sfts(fmax, Tsft, noise_sqrt_Sh=noise_sqrt_Sh, noise_seed=noise_seed, window=window, window_param=window_param):

            # create standard SFT file name (see LIGO-T040164)
            sft_desc = 'simCW_%s' % comment
            sft_name = lalpulsar.GetOfficialName4SFT(sft, sft_desc)
            sft_path = os.path.join(out_dir, sft_name)

            # write SFT
            lalpulsar.WriteSFT2file(sft, sft_path, self.__origin_str)

            # yield current file name for e.g. printing progress
            yield sft_path, i, N

## @}
