from astropy.coordinates import Angle, SkyCoord
from astropy.time import Time
from astropy import units as u
from typing import Dict, NamedTuple, Union
from gwpy.timeseries import TimeSeries
from gwpy.frequencyseries import FrequencySeries
import lal
import lalsimulation as lalsim
import numpy as np
import warnings

class GravitationalWavePolarizations(NamedTuple):
    hp: Union[TimeSeries, FrequencySeries]
    hc: Union[TimeSeries, FrequencySeries]
    """
    GravitationalWavePolarizations class that takes hp,hc GwPy Time/FrequencySeries objects as inputs. Provides method to compute
    zero noise strain given a detector, sky-position, polarization and time of arrival values.
    Parameters
    ----------
    hp, hc : GwPy `TimeSeries` or `FrequencySeries` objects
    """
    def domain(self):
        if isinstance(self.hp, TimeSeries) and isinstance(self.hc, TimeSeries):
            return 'time'
        elif isinstance(self.hp, FrequencySeries) and isinstance(self.hc, FrequencySeries):
            return 'frequency'
        else:
            return 'mixed'

    def strain(self, det, ra, dec, psi, tgps):
        """
        Compute the detector strain in zero-noise

        Parameters
        ----------
        det : `str`
            Name of the detector, for eg: 'H1' or 'L1'
        ra : `~astropy.units.Quantity`
            Right-ascension of the binary in units of rad
        dec : `~astropy.units.Quantity`
            Declination of the binary in units of rad
        psi : `~astropy.units.Quantity`
            Polarization in units of rad
        tgps : `float`, `~astropy.units.Quantity`
            GPS Time of time-of-arrival of signal. Either as float or in units of s

        Returns
        -------
        h : `TimeSeries` or `FrequencySeries` object
            Waveform recomposed with detector response (Fp*hp + Fc*hc)
        """

        # Might change it later; for now keeping it as used by ligo.
        warnings.warn("This code is currently UNREVIEWED, use with caution!!")
        pos = SkyCoord(ra = ra, dec = dec)
        psi = Angle(psi)
        t = Time(tgps, format='gps', scale='utc')

        if isinstance(det, str):
            det = lalsim.DetectorPrefixToLALDetector(det)

        if self.domain() == 'time':
            hp = self.hp.to_lal()
            hc = self.hc.to_lal()
            hp.epoch += t.gps
            hc.epoch += t.gps
            h = lalsim.SimDetectorStrainREAL8TimeSeries(hp, hc, pos.ra.rad, pos.dec.rad, psi.rad, det)
            h = TimeSeries.from_lal(h)
            return h

        elif self.domain() == 'frequency':
            # WARNING: does not account for Earth rotation or non-LWL effects

            epoch = lal.LIGOTimeGPS(t.gps)
            dt = lal.TimeDelayFromEarthCenter(det.location, pos.ra.rad, pos.dec.rad, epoch)
            fp, fc = lal.ComputeDetAMResponse(det.response, pos.ra.rad, pos.dec.rad, psi.rad, lal.GreenwichMeanSiderealTime(epoch))

            # Shouldn't this have a overall time shift due to time-delay from earth center? or does epoch take care of it?
            # TODO: adress during phase 2 of review

            h = (fp * self.hp + fc * self.hc)
            h.epoch = Time(t.gps + dt, format='gps', scale='utc')
            return h
        else:
            return ValueError('hp and hc must both be either TimeSeries or FrequencySeries')

    def add_strain(self, data, det, ra, dec, psi, tgps, response=None):
         """
         Add strain to some LAL strain as required

         Parameters:
         data : LAL TimeSeries data
            Data in which to add detector response strain
        det : `str`
            Name of the detector, for eg: 'H1' or 'L1'
        ra : `~astropy.units.Quantity`
            Right-ascension of the binary in units of rad
        dec : `~astropy.units.Quantity`
            Declination of the binary in units of rad
        psi : `~astropy.units.Quantity`
            Polarization in units of rad
        tgps : `float`, `~astropy.units.Quantity`
            GPS Time of time-of-arrival of signal. Either as float or in units of s
        response : `COMPLEX16FrequencySeries`
            Response function passed to lalsimulation.SimAddInjectionREAL8TimeSeries,
            transforming strain to detector output units. Use None for unit response.

        Returns
        -------
        data : `TimeSeries` object
            Detector response injected in strain data
         """
         warnings.warn("This code is currently UNREVIEWED, use with caution!")
         h = self.strain(det, ra, dec, psi, tgps)
         # Adding condition to convert FrequencySeries to TimeSeries if strain is in freq domain
         if isinstance(h, FrequencySeries):
             h = h.ifft()
         h = h.to_lal()
         data = data.to_lal()
         lalsim.SimAddInjectionREAL8TimeSeries(data, h, response)
         data = TimeSeries.from_lal(data)
         return data


class ModePair(NamedTuple):
    """
    Returns a named tuple given the l,m values
    """
    l: int
    m: int


class SpinWeightedSphericalHarmonicMode(ModePair):
    """
    Class to return spin `s` weighted spherical harmonic given ModePair
    """

    def __new__(cls, s, l, m):
        new = super().__new__(cls, l, m)
        if new.l < abs(s):
            raise ValueError('Require l >= |s|')
        if abs(new.m) > new.l:
            raise ValueError('Require |m| <= l')
        new.s = s
        return new

    def __call__(self, theta, phi):
        theta = Angle(theta, u.rad)
        phi = Angle(phi, u.rad)
        return lal.SpinWeightedSphericalHarmonic(theta.rad, phi.rad, self.s, self.l, self.m)


class TensorSphericalHarmonicMode(SpinWeightedSphericalHarmonicMode):
    """
    Class to return spin `-2` weighted spherical harmonic given ModePair
    """
    def __new__(cls, l, m):
        return super().__new__(cls, -2, l, m)


class GravitationalWaveModes(Dict[SpinWeightedSphericalHarmonicMode, Union[TimeSeries, FrequencySeries]]):
    """
    Class for gravitational wave modes which returns the  waveform recomposed with -2 spin weighted spherical harmonics
    given a (theta, phi) value
    """
    def __call__(self, theta, phi):
        """
        Return plus and cross polarizations from the gravitational wave modes.

        Parameters
        ----------
        theta : `~astropy.units.Quantity`
            Theta angle (inclination, theta_jn) in units of rad.
        phi : `~astropy.units.Quantity`
            Phi angle (inclination, theta_jn) in units of rad.
        """
        if isinstance(theta, float) and isinstance(phi, float):
            raise TypeError("Theta and phi should be in units of astropy.units.rad, not float")
        elif theta.unit.is_equivalent(u.rad) and phi.unit.is_equivalent(u.rad):
            pass
        else:
            raise TypeError("Theta and phi should be in units of astropy.units.rad")


        hpc = 0.
        hp = 0.
        hc = 0.
        for key, hlm in self.items():
            if isinstance(key, str):
                pass
            elif isinstance(hlm, TimeSeries):
                mode = TensorSphericalHarmonicMode(int(key[0]), int(key[1]))
                hpc += mode(theta, phi) * hlm
                hp = hpc.real
                hc = -hpc.imag

            elif isinstance(hlm, FrequencySeries):
                # Recombination of freq domain hp/hc from here :
                # https://git.ligo.org/lscsoft/lalsuite/-/blob/master/lalsimulation/lib/LALSimInspiral.c#L3477
                l, m = int(key[0]), int(key[1])
                mode = TensorSphericalHarmonicMode(l, m)
                ylm = mode(theta, phi)
                hp  += 0.5*( ylm * hlm + np.conj(ylm) * np.conj(hlm)[::-1])
                hc  += 1j*0.5*( ylm * hlm - np.conj(ylm)* np.conj(hlm)[::-1])

        return GravitationalWavePolarizations(hp, hc)
