# Copyright (C) 2012 Chris Pankow
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
import math

import numpy
from lal.rate import BinnedDensity, NDBins, IrregularBins, LinearBins

__author__ = "Chris Pankow <chris.pankow@ligo.org>"

#
# Stat utilities
#
def cumvar(arr, mean=None, var=None, n=0):
	"""
	Numerically stable running sample variance measure. If mean and var are supplied, they will be used as the history values. See 

    http://www.johndcook.com/standard_deviation.html

    for algorithm details.
	"""
	if mean and var:
		m, s = numpy.zeros(len(arr)+1), numpy.zeros(len(arr)+1, dtype=numpy.float128)
		m[0] = mean
		s[0] = var*(n-1)
		buf = numpy.array([0])
	else:
		m, s = numpy.zeros(arr.shape), numpy.zeros(arr.shape, dtype=numpy.float128)
		m[0] = arr[0]
		buf = numpy.array([])

	for i, x in enumerate(numpy.concatenate((buf, arr))):
		if mean is None:
			k = i+1+n
		else:
			k = i+n
		if i == 0: continue
		m[i] = m[i-1] + (x-m[i-1])/k
		s[i] = s[i-1] + (x-m[i-1])*(x-m[i])

	if mean and var:
		return s[1:]/numpy.arange(n, n + len(s)-1)
	else:
		norm = numpy.arange(n, n + len(s))
		norm[0] = 1 # avoid a warning about zero division
		return s/norm

def int_var(samples):
    mean = numpy.mean(samples)
    sq_mean = numpy.mean(samples**2)
    return (sq_mean-mean**2)/(len(samples)-1)

def neff(weights):
    neff = reduce(numpy.logaddexp, w) - max(w)
    return math.exp(neff)

def neff_frac(weights, ntotal):
    return neff / ntotal

#
# Histogramming utilities
#
def get_adaptive_binning(samples, edges, nbins=100, bintype='linear'):

    if bintype == "linear":
        return BinnedDensity(NDBins((LinearBins(edges[0], edges[-1], nbins),)))

    ordering = numpy.argsort(samples)
    stride = len(samples) / nbins
    bins = [samples[ordering[i]] for i in xrange(stride, nbins*stride, stride)]
    bins.insert(0, edges[0])
    bins.append(edges[1])
    return BinnedArray(NDBins((IrregularBins(bins),)))

def construct_alias_table(probs):
    """
    FIXME: instead of unity, use the normalization of the probabilities.
    Based on an expansive exposition and tutorial at http://www.keithschwarz.com/darts-dice-coins/

    Using the array of discrete, and sorted probabilities, Construct an alias table. This is accomplished by normalizing by the averge probability (typically 1/N), and removing probability mass from the last (largest) probability and "giving" it to the smallest so as the total of the removed mass and the current smallest probability is now unity. The alias table then records the index of the larger mass as the "alias" for the biased coin flip.

    The next smallest probability is chosen and the procedure repeated until the (formerly) largest probability also has unity mass, then the next largest probability is selected and the process continues untill all of the renomralized probabilities would have unity mass.

    What is returned is the alias table (of indexes) and the biased coin flip probabilities.
    """

    avg_prob = sum(probs)/len(probs)
    probs = [p / avg_prob for p in probs]
    #print "Average prob: %f" % avg_prob
    alias = numpy.linspace(0, len(probs)-1, len(probs), dtype=numpy.int64)

    #print probs
    # index running from the beginning (i) and end (j)
    i, j = 0, len(probs) - 1
    # FIXME: I *think* this can be i <= j or something like that.
    while i < len(probs) and j >= 0:
        #print "i, j: %f %f" % (probs[i], probs[j])
        # We "take away" mass from probs[j] to equalize to one probs[i], that
        # gives us the 'alias':
        # for i, e.g. alias[i] = j (iff probs[i] < 1)
        while probs[j] > 1 and i < len(probs):
            # This usually happens when i >= j and we need to redistribute the
            # last bit of prob mass between what used to be the larger values.
            if probs[i] >= 1:
                i += 1
                continue
            alias[i] = j
            #print "aliasing %d (%f) with %f from %d" % (i, probs[i], (1 - probs[i]), j)
            assert 1 - probs[i] > 0
            probs[j] -= 1 - probs[i]
            #print "Probs %d now %f" % (j, probs[j])
            #print probs
            i += 1
        j -= 1
    return numpy.array(probs), numpy.array(alias)

def draw_from_alias_table(cprobs, alias, ndraws=1):
    """
    For N draws, make N 'coin flips' and N index draws --- if the coin flip is under the biased number, use that index, otherwise, use the aliased index.
    """
    idxs, coins = numpy.random.randint(0, len(cprobs), ndraws), numpy.random.uniform(0, 1, ndraws)
    return numpy.where(coins >= cprobs[idxs], alias[idxs], idxs)

class AliasSampler(object):
    """
    A sampler designed to emulate an n-dimensional histogram, with amortized O(1) per sample draw times (compare to rejection sampling). The cost of initialization is something like N log N (not confirmed) where N is the number of histogram bins. This is effectively a memory for speed tradeoff, so high dimensional binning make take a long time to initialize and a lot of memory to carry out.
    """
    def __init__(self, binning):
        """
        The 'binning' argument is a pylal BinnedArray object, which is assumed to be a n-dimensional histogram.
        """

        # pylal.rate "Bins" object
        self._bins = binning
        # Number of actual bins (product of 1-d binning)
        self._nbins = numpy.prod(self._bins.bins.shape)
        # Unravelled bin indexing
        self._flat_binning = numpy.zeros(self._nbins)
        # N-dimensional lower and upper bin edges
        # These are used to produce uniformly distributed numbers *within* a bin
        self._nd_low = numpy.array(self._bins.bins.lower())
        self._nd_high = numpy.array(self._bins.bins.upper())

        vols = self._bins.bins.volumes()
        # The flat binning is filled with the probability *mass* within a bin
        for i, itr in enumerate(numpy.ndindex(*self._bins.array.shape)):
            self._flat_binning[i] = self._bins.array[itr] * vols[itr]

        #  ...and then normalized
        # FIXME: If the distribution is not frozen, _finalize would need to be
        # called manually before sampling again
        self._finalize()

    def _finalize(self):
        """
        This calculates the current normalization constant and reconstructs the alias table. It is not recommended to call this often --- only when the bin values change. Note that this is called at the end of the __init__ method, so the user should not have to call it for a frozen distribution.
        """
        norm = float(self._flat_binning.sum())
        if norm != 0:
            self._flat_binning[i] /= norm

        # We keep the bins sorted by probability to speed up the alias table
        # generation -- you have to remap the indices returned by
        # draw_from_alias_table
        self._sort = numpy.argsort(self._flat_binning)
        self._cprobs, self._alias = construct_alias_table(self._flat_binning[self._sort])

    def add_samples(self, samps):
        """
        Add samples to the distribution. Not implemented at this time so as to keep the distribution nominally frozen.
        """
        raise NotImplementedError("Distributions must be frozen.")
        #for s in samps:
            #self._bins[s] += 1

    def draw(self, ndraws):
        """
        Obtain 'ndraws' output samples from the histogram. Output is a ndim x ndraws numpy array.
        """
        # Make ndraws from the alias table --- the flattened indices are in
        # sorted order though, so we have to remap them back to the proper bins
        idx = self._sort[draw_from_alias_table(self._cprobs, self._alias, ndraws)]

        # Preconstruct the output draws
        draws = numpy.random.uniform(size=(self._bins.array.ndim, ndraws))

        # Convert 'idx' to a multi-dimensional index to get the bin we're
        # drawing from along the axis 'dim'
        for dim, nd_idx in enumerate(numpy.unravel_index(idx, self._bins.array.shape)):
            # Iterate over each dimension, shifting the random number to the bin
            # FIXME: We *could* precalculate this in _finalize
            ab = self._nd_high[dim] - self._nd_low[dim]
            draws[dim,:] = ab[nd_idx] * draws[dim,:] + self._nd_low[dim][nd_idx]

        return draws

#
# Test suite
#
if __name__ == "__main__":
    from collections import Counter
    from numpy.random import multivariate_normal

    #print construct_alias_table(sorted([1/4., 1/4., 1/4., 1/4.]))
    probs = sorted([1/12., 1/3., 1/2., 1/12.])
    cprobs, alias = construct_alias_table(probs)

    idx = draw_from_alias_table(cprobs, alias, 100000)

    bin_chosen = Counter()
    for i in idx:
        bin_chosen[i] += 1
    bin_chosen = numpy.array(bin_chosen.values(), dtype=numpy.float64)
    bin_chosen /= bin_chosen.sum()
    print sorted(probs)
    print sorted(bin_chosen)

    rvs = multivariate_normal(mean=(0,0), cov=((1, 2), (0.5, 2)), size=(100000,)).T
    nbin_side = 20
    binning = BinnedArray(NDBins((LinearBins(-3, 3, nbin_side), LinearBins(-3, 3, nbin_side))))
    for rv in rvs.T:
        if any(numpy.abs(rv) > 3):
            continue
        binning[tuple(rv)] += 1
    samp = AliasSampler(binning)
    new_rvs = samp.draw(100000)

    import matplotlib
    matplotlib.use("agg")
    from matplotlib import pyplot
    pyplot.subplot(1, 2, 1)
    h, x, y = numpy.histogram2d(rvs[0], rvs[1], bins=(20, 20), range=((-3, 3), (-3, 3)))
    pyplot.pcolor(h.T)
    pyplot.subplot(1, 2, 2)
    h, x, y = numpy.histogram2d(new_rvs[0], new_rvs[1], bins=(20, 20), range=((-3, 3), (-3, 3)))
    pyplot.pcolor(h.T)
    pyplot.savefig("alias_sampler_test.png")
