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
    from pylal.xlal.constants import LAL_PI, LAL_MTSUN_SI
except ImportError:
    raise ImportError("The sbank subpackage of lalinspiral depends on the glue and pylal packages.")

from lalinspiral import CreateSBankWorkspaceCache
from lalinspiral.sbank.psds import get_neighborhood_ASD, get_neighborhood_PSD

class lazy_mchirps(object):
    __slots__ = ("seq",)

    def __init__(self, seq):
        self.seq = seq

    def __getitem__(self, ind):
        return self.seq[ind]._mchirp

    def __len__(self):
        return len(self.seq)

class Bank(object):
    __slots__ = ("tmplt_class", "noise_model", "flow", "use_metric", "_templates", "_nmatch", "_mchirps", "compute_match", "_workspace_cache", "_moments")

    def __init__(self, tmplt_class, noise_model, flow, use_metric=False):
        self.tmplt_class = tmplt_class
        self.noise_model = noise_model
        self.flow = flow
        self.use_metric = use_metric

        self._templates = []
        self._nmatch = 0
        self._mchirps = lazy_mchirps(self._templates)
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
        ind = bisect.bisect_left(self._mchirps, new._mchirp)
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
        merged._templates[:] = list(inorder(banks, key=attrgetter("_mchirp")))
        merged._nmatch = sum(b._nmatch for b in banks)
        return merged

    @classmethod
    def from_sngls(cls, sngls, tmplt_class, *args):
        bank = cls(*((tmplt_class,) + args))
        bank._templates.extend([tmplt_class.from_sngl(s, bank=bank) for s in sngls])
        bank._templates.sort(key=attrgetter("_mchirp"))
        return bank

    @classmethod
    def from_sims(cls, sims, tmplt_class, *args):
        bank = cls(*((tmplt_class,) + args))
        bank._templates.extend([tmplt_class.from_sim(s, bank=bank) for s in sims])
        bank._templates.sort(key=attrgetter("_mchirp"))
        return bank

    def _metric_match(self, tmplt, proposal, f, **kwargs):
        return tmplt.metric_match(proposal, f, **kwargs)

    def _brute_match(self, tmplt, proposal, f, **kwargs):
        return tmplt.brute_match(proposal, f,self._workspace_cache, **kwargs)

    def covers(self, proposal, min_match):
        """
        Return True if any template in the bank has match with proposal
        greater than min_match.
        """
        # find templates in the bank "near" this tmplt
        low, high = _find_neighborhood(self._mchirps, proposal._mchirp, self.flow)
        tmpbank = self._templates[low:high]
        if not tmpbank: return False

        # sort the bank by its nearness to tmplt in mchirp
        # NB: This sort comes up as a dominating cost if you profile,
        # but it cuts the number of match evaluations by 80%, so turns out
        # to be worth it even for metric match, where matches are cheap.
        pm = proposal._mchirp
        tmpbank.sort(key=lambda b: abs(b._mchirp - pm))

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
        # find templates near the proposal
        low, high = _find_neighborhood(self._mchirps, proposal._mchirp, self.flow)
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


def _find_neighborhood(tmplt_mcs, mc0, flow, dt0=0.25):
    """
    Return the min and max indices of templates that cover the given
    chirp mass (mc0) within a tau0 difference of dt0 (seconds).
    tmplt_mcs should be a sequence of mchirp values in sorted order.
    """
    tau0 = 5. * mc0 * LAL_MTSUN_SI / (256 * (LAL_PI * flow * mc0 * LAL_MTSUN_SI)**(8./3))
    low = mc0 * (1 - 0.6 * dt0 / tau0)  # Taylor expand mchirp = (tau0/A0)**-0.6
    high = mc0 * (1 + 0.6 * dt0 / tau0)
    low_ind = bisect.bisect_left(tmplt_mcs, low)
    high_ind = bisect.bisect_right(tmplt_mcs, high)
    return low_ind, high_ind
