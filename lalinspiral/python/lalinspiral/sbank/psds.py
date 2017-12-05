# Copyright (C) 2011  Nickolas Fotopoulos
# Copyright (C) 2014-2017  Stephen Privitera
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

from __future__ import division

from math import log, ceil

from numpy import vectorize, arange, seterr, inf, ones_like
seterr(over="ignore")  # the PSD overflows frequently, but that's OK

from glue.ligolw import param
from glue.ligolw import utils
from lal import series as lalseries
import lal
import lalsimulation as lalsim


#
# Analytical PSDs
#
noise_models = {
    "LIGOIPsd": vectorize(lambda f: lal.LIGOIPsd(f)),  # what most lalapps_* uses
    "AdvLIGOPsd": vectorize(lambda f: 1e-49 * lal.AdvLIGOPsd(f)),  # what most lalapps_* uses; hasn't been XLALified, so doesn't have scaling factor. Gah!
    "iLIGOSRD": vectorize(lambda f: lalsim.SimNoisePSDiLIGOSRD(f)),
    "aLIGONoSRMLowPower": vectorize(lambda f: lalsim.SimNoisePSDaLIGONoSRMLowPower(f)),
    "aLIGONoSRMHighPower": vectorize(lambda f: lalsim.SimNoisePSDaLIGONoSRMHighPower(f)),
    "aLIGONSNSOpt": vectorize(lambda f: lalsim.SimNoisePSDaLIGONSNSOpt(f)),
    "aLIGOBHBH20Deg": vectorize(lambda f: lalsim.SimNoisePSDaLIGOBHBH20Deg(f)),
    "aLIGOHighFrequency": vectorize(lambda f: lalsim.SimNoisePSDaLIGOHighFrequency(f)),
    "aLIGOZeroDetLowPower": vectorize(lambda f: lalsim.SimNoisePSDaLIGOZeroDetLowPower(f)),
    "aLIGOZeroDetHighPower": vectorize(lambda f: lalsim.SimNoisePSDaLIGOZeroDetHighPower(f)),
}


psd_cache = {}  # keyed by df, flow, f_max
def get_PSD(df, flow, f_max, noise_model):
    """
    Return the frequency vector and sampled PSD using the noise_model,
    which are vectorized wrappings of the noise models in lalsimulation.
    flow sets the first non-zero frequency.
    """
    if (df, flow, f_max) in psd_cache:
        return psd_cache[df, flow, f_max]
    f = arange(f_max / df + 1) * df
    LIGO_PSD = inf * ones_like(f)
    ind_low = int(flow / df)
    LIGO_PSD[ind_low:] = noise_model(f[ind_low:])
    return psd_cache.setdefault((df, flow, f_max), LIGO_PSD)

asd_cache = {}
def get_ASD(df, flow, f_max, noise_model):
    """
    Some routines prefer ASDs over PSDs. Keep cache of ASDs, but we may
    be able to speed things up by drawing upon cache of PSDs.
    """
    if (df, flow, f_max) in asd_cache:
        return asd_cache[df, flow, f_max]
    return asd_cache.setdefault((df, flow, f_max), get_PSD(df, flow, f_max, noise_model)**0.5)

#
# Determine PSD for neighborhood
#

def next_pow2(n):
    return 1 << int(ceil(log(n, 2)))

def prev_pow2(n):
    return 1 << int(log(n, 2))

def get_neighborhood_df_fmax(waveforms, flow):
    """
    Return PSD that is optimized for this neighborhood, with small enough
    df and big enough f_max to cover all waveforms.
    """
    max_dur = max(w.dur for w in waveforms)
    assert 16384 * max_dur > 1   # chirp lasts long enough for one LIGO sample
    if max_dur >= 1:
        df = 1 / next_pow2(max_dur)
    else:
        df = prev_pow2(1 / max_dur)
    max_ffinal = max(w.f_final for w in waveforms)
    f_max = next_pow2(max_ffinal)  # will always be greater than 1
    assert f_max - flow >= 2 * df  # need a few frequencies at least!
    return df, f_max

def get_neighborhood_PSD(waveforms, flow, noise_model):
    """
    Return PSD that is optimized for this neighborhood, with small enough
    df and big enough f_max to cover all waveforms.
    """
    max_dur = max(w.dur for w in waveforms)
    assert 16384 * max_dur > 1   # chirp lasts long enough for one LIGO sample
    if max_dur >= 1:
        df = 1 / next_pow2(max_dur)
    else:
        df = prev_pow2(1 / max_dur)
    max_ffinal = max(w.f_final for w in waveforms)
    f_max = next_pow2(max_ffinal)  # will always be greater than 1
    assert f_max - flow >= 2 * df  # need a few frequencies at least!
    return df, get_PSD(df, flow, f_max, noise_model)

def get_neighborhood_ASD(waveforms, flow, noise_model):
    """
    Return ASD that is optimized for this neighborhood, with small enough
    df and big enough f_max to cover all waveforms.
    """
    max_dur = max(w.dur for w in waveforms)
    assert 16384 * max_dur > 1   # chirp lasts long enough for one LIGO sample
    if max_dur >= 1:
        df = 1 / next_pow2(max_dur)
    else:
        df = prev_pow2(1 / max_dur)
    max_ffinal = max(w.f_final for w in waveforms)
    f_max = next_pow2(max_ffinal)  # will always be greater than 1
    assert f_max - flow >= 2 * df  # need a few frequencies at least!
    return df, get_ASD(df, flow, f_max, noise_model)

#
# Read PSD from arrays in XML documents
#

def psd_instrument_dict(elem):
    out = {}
    for lw in elem.getElementsByTagName(u"LIGO_LW"):
        if not lw.hasAttribute(u"Name"):
            continue
        if lw.getAttribute(u"Name") != u"REAL8FrequencySeries":
            continue
        ifo = param.get_pyvalue(lw, u"instrument")
        out[ifo] = lalseries.parse_REAL8FrequencySeries(lw)
    return out

def read_psd(filename, verbose = False):
    return psd_instrument_dict(utils.load_filename(filename, verbose = verbose, contenthandler=lalseries.PSDContentHandler))
