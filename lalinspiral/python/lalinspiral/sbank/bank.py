#!/usr/bin/env python
# Copyright (C) 2012  Nickolas Fotopoulos
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

import bisect
from operator import attrgetter

try:
    from glue.iterutils import inorder, uniq
except ImportError:
    raise ImportError("The sbank subpackage of lalinspiral depends on the glue and pylal packages.")

from lal import PI, MTSUN_SI
from lalinspiral import CreateSBankWorkspaceCache
from lalinspiral.sbank.psds import get_neighborhood_ASD, get_neighborhood_PSD

class lazy_nhoods(object):
    __slots__ = ("seq", "nhood_param")
    def __init__(self, seq, nhood_param="tau0"):
        self.seq = seq
        self.nhood_param = nhood_param

    def __getitem__(self, idx):
        return getattr(self.seq[idx], self.nhood_param)

    def __len__(self):
        return len(self.seq)

class Bank(object):

    def __init__(self, tmplt_class, noise_model, flow, use_metric=False, cache_waveforms=False, nhood_size=1.0, nhood_param="tau0"):
        self.tmplt_class = tmplt_class
        self.noise_model = noise_model
        self.flow = flow
        self.use_metric = use_metric
        self.cache_waveforms = cache_waveforms
        self.nhood_size = nhood_size
        self.nhood_param = "_" + nhood_param

        self._templates = []
        self._nmatch = 0
        self._nhoods = lazy_nhoods(self._templates, self.nhood_param)
        if use_metric:
            self._moments = {}
            self.compute_match = self._metric_match
        else:
            self._workspace_cache = CreateSBankWorkspaceCache()
            self.compute_match = self._brute_match

    def __len__(self):
        return len(self._templates)

    def __iter__(self):
        return iter(self._templates)

    def __repr__(self):
        return repr(self._templates)

    def insort(self, new):
        ind = bisect.bisect_left(self._nhoods, getattr(new, self.nhood_param))
        self._templates.insert(ind, new)
        new.finalize_as_template()

    @classmethod
    def merge(cls, *banks):
        if not banks:
            raise ValueError("cannot merge zero banks")
        cls_args = list(uniq((b.tmplt_class, b.noise_model, b.flow, b.use_metric) for b in banks))
        if len(cls_args) > 1:
            return ValueError("bank parameters do not match")
        merged = cls(*cls_args)
        merged._templates[:] = list(inorder(banks, key=attrgetter(bank.nhood_param)))
        merged._nmatch = sum(b._nmatch for b in banks)
        return merged

    @classmethod
    def from_sngls(cls, sngls, tmplt_class, *args, **kwargs):
        bank = cls(*((tmplt_class,) + args), **kwargs)
        bank._templates.extend([tmplt_class.from_sngl(s, bank=bank) for s in sngls])
        bank._templates.sort(key=attrgetter(bank.nhood_param))
        # Mark all templates as seed points
        for template in bank._templates:
            template.is_seed_point = True
        return bank

    @classmethod
    def from_sims(cls, sims, tmplt_class, *args):
        bank = cls(*((tmplt_class,) + args))
        bank._templates.extend([tmplt_class.from_sim(s, bank=bank) for s in sims])
        return bank

    def _metric_match(self, tmplt, proposal, f, **kwargs):
        return tmplt.metric_match(proposal, f, **kwargs)

    def _brute_match(self, tmplt, proposal, f, **kwargs):
        match = tmplt.brute_match(proposal, f,self._workspace_cache, **kwargs)
        if not self.cache_waveforms:
            tmplt.clear()
        return match

    def covers(self, proposal, min_match):
        """
        Return True if any template in the bank has match with proposal
        greater than min_match.
        """
        # find templates in the bank "near" this tmplt
        prop_nhd = getattr(proposal, self.nhood_param)
        low, high = _find_neighborhood(self._nhoods, prop_nhd, self.nhood_size)
        tmpbank = self._templates[low:high]
        if not tmpbank: return False

        # sort the bank by its nearness to tmplt in mchirp
        # NB: This sort comes up as a dominating cost if you profile,
        # but it cuts the number of match evaluations by 80%, so turns out
        # to be worth it even for metric match, where matches are cheap.
        tmpbank.sort(key=lambda b: abs( getattr(b, self.nhood_param) - prop_nhd))

        # set parameters of match calculation that are optimized for this block
        df, PSD = get_neighborhood_PSD(tmpbank + [proposal], self.flow, self.noise_model)

        # find and test matches
        for tmplt in tmpbank:
            self._nmatch += 1
            match = self.compute_match(tmplt, proposal, df, PSD=PSD)
            if match > min_match:
                return True
        return False

    def max_match(self, proposal):
        match, best_tmplt_ind = self.argmax_match(proposal)
        if not match: return (0., 0)
        return match, self._templates[best_tmplt_ind]

    def argmax_match(self, proposal):
        # find templates in the bank "near" this tmplt
        prop_nhd = getattr(proposal, self.nhood_param)
        low, high = _find_neighborhood(self._nhoods, prop_nhd, self.nhood_size)
        tmpbank = self._templates[low:high]
        if not tmpbank: return (0., 0)

        # set parameters of match calculation that are optimized for this block
        df, ASD = get_neighborhood_ASD(tmpbank + [proposal], self.flow, self.noise_model)

        # compute matches
        match, best_tmplt_ind = max((self.compute_match(tmplt, proposal, df, ASD=ASD), ind) for ind, tmplt in enumerate(tmpbank))
        self._nmatch += len(tmpbank)

        # best_tmplt_ind indexes tmpbank; add low to get index of full bank
        return match, best_tmplt_ind + low

    def clear(self):
        if hasattr(self, "_workspace_cache"):
            self._workspace_cache = CreateSBankWorkspaceCache()
        for tmplt in self._templates:
            tmplt.clear()


def _find_neighborhood(tmplt_locs, prop_loc, nhood_size=0.25):
    """
    Return the min and max indices of templates that cover the given
    template at prop_loc within a parameter difference of nhood_size (seconds).
    tmplt_locs should be a sequence of neighborhood values in sorted order.
    """
    low_ind = bisect.bisect_left(tmplt_locs, prop_loc - nhood_size)
    high_ind = bisect.bisect_right(tmplt_locs, prop_loc + nhood_size)
    return low_ind, high_ind
