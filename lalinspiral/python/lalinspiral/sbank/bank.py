#!/usr/bin/env python
# Copyright (C) 2012  Nickolas Fotopoulos
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

import bisect
from operator import attrgetter

import numpy as np

from glue.iterutils import inorder, uniq
from lal import PI, MTSUN_SI
from lalinspiral import CreateSBankWorkspaceCache
from lalinspiral.sbank.psds import get_neighborhood_ASD, get_neighborhood_PSD, get_PSD, get_neighborhood_df_fmax
from . import waveforms

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

    def __init__(self, noise_model, flow, use_metric=False, cache_waveforms=False, nhood_size=1.0, nhood_param="tau0", coarse_match_df=None, iterative_match_df_max=None, fhigh_max=None, optimize_flow=None, flow_column=None):

        self.noise_model = noise_model
        self.flow = flow
        self.use_metric = use_metric
        self.cache_waveforms = cache_waveforms
        self.coarse_match_df = coarse_match_df
        self.iterative_match_df_max = iterative_match_df_max
        self.optimize_flow = optimize_flow
        self.flow_column = flow_column

        if self.coarse_match_df and self.iterative_match_df_max and self.coarse_match_df < self.iterative_match_df_max:
            # If this case occurs coarse_match_df offers no improvement, turn off
            self.coarse_match_df = None

        if fhigh_max is not None:
            self.fhigh_max = (2**(np.ceil(np.log2( fhigh_max ))))
        else:
            self.fhigh_max = fhigh_max

        self.nhood_size = nhood_size
        self.nhood_param = nhood_param

        self._templates = []
        self._nmatch = 0
        self._nhoods = lazy_nhoods(self._templates, self.nhood_param)
        if use_metric:
            self._moments = {}
            self.compute_match = self._metric_match
        else:
            # The max over skyloc stuff needs a second cache
            self._workspace_cache = [CreateSBankWorkspaceCache(), CreateSBankWorkspaceCache()]
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

    def add_from_sngls(self, sngls, tmplt_class):
        newtmplts = [tmplt_class.from_sngl(s, bank=self) for s in sngls]
        for template in newtmplts:
            # Mark all templates as seed points
            template.is_seed_point = True
        self._templates.extend(newtmplts)
        self._templates.sort(key=attrgetter(self.nhood_param))

    def add_from_hdf(self, hdf_fp):
        num_points = len(hdf_fp['mass1'])
        newtmplts=[]
        for idx in xrange(num_points):
            if not idx % 100000:
                tmp = {}
                end_idx = min(idx+100000, num_points)
                for name in hdf_fp:
                    tmp[name] = hdf_fp[name][idx:end_idx]
            c_idx = idx % 100000
            approx = tmp['approximant'][c_idx]
            tmplt_class = waveforms.waveforms[approx]
            newtmplts.append(tmplt_class.from_dict(tmp, c_idx, self))
            newtmplts[-1].is_seed_point=True
        self._templates.extend(newtmplts)
        self._templates.sort(key=attrgetter(self.nhood_param))

    @classmethod
    def from_sngls(cls, sngls, tmplt_class, *args, **kwargs):
        bank = cls(*args, **kwargs)
        bank._templates.extend([tmplt_class.from_sngl(s, bank=bank) for s in sngls])
        bank._templates.sort(key=attrgetter(bank.nhood_param))
        # Mark all templates as seed points
        for template in bank._templates:
            template.is_seed_point = True
        return bank

    @classmethod
    def from_sims(cls, sims, tmplt_class, *args):
        bank = cls(*args)
        bank._templates.extend([tmplt_class.from_sim(s, bank=bank) for s in sims])
        return bank

    def _metric_match(self, tmplt, proposal, f, **kwargs):
        return tmplt.metric_match(proposal, f, **kwargs)

    def _brute_match(self, tmplt, proposal, f, **kwargs):
        match = tmplt.brute_match(proposal, f,self._workspace_cache, **kwargs)
        if not self.cache_waveforms:
            tmplt.clear()
        return match

    def covers(self, proposal, min_match, nhood=None):
        """
        Return (max_match, template) where max_match is either (i) the
        best found match if max_match < min_match or (ii) the match of
        the first template found with match >= min_match.  template is
        the Template() object which yields max_match.
        """
        max_match = 0
        template = None

        # find templates in the bank "near" this tmplt
        prop_nhd = getattr(proposal, self.nhood_param)
        if not nhood:
            low, high = _find_neighborhood(self._nhoods, prop_nhd, self.nhood_size)
            tmpbank = self._templates[low:high]
        else:
            tmpbank = nhood
        if not tmpbank: return (max_match, template)

        # sort the bank by its nearness to tmplt in mchirp
        # NB: This sort comes up as a dominating cost if you profile,
        # but it cuts the number of match evaluations by 80%, so turns out
        # to be worth it even for metric match, where matches are cheap.
        tmpbank.sort(key=lambda b: abs( getattr(b, self.nhood_param) - prop_nhd))

        # set parameters of match calculation that are optimized for this block
        df_end, f_max = get_neighborhood_df_fmax(tmpbank + [proposal], self.flow)
        if self.fhigh_max:
            f_max = min(f_max, self.fhigh_max)
        df_start = max(df_end, self.iterative_match_df_max)

        # find and test matches
        for tmplt in tmpbank:

            self._nmatch += 1
            df = df_start
            match_last = 0

            if self.coarse_match_df:
                # Perform a match at high df to see if point can be quickly
                # ruled out as already covering the proposal
                PSD = get_PSD(self.coarse_match_df, self.flow, f_max, self.noise_model)
                match = self.compute_match(tmplt, proposal, self.coarse_match_df,
                                           PSD=PSD)
                if match == 0:
                    err_msg = "Match is 0. This might indicate that you have "
                    err_msg += "the df value too high. Please try setting the "
                    err_msg += "coarse-value-df value lower."
                    # FIXME: This could be dealt with dynamically??
                    raise ValueError(err_msg)

                if (1 - match) > 0.05 + (1 - min_match):
                    continue

            while df >= df_end:

                PSD = get_PSD(df, self.flow, f_max, self.noise_model)
                match = self.compute_match(tmplt, proposal, df, PSD=PSD)
                if match == 0:
                    err_msg = "Match is 0. This might indicate that you have "
                    err_msg += "the df value too high. Please try setting the "
                    err_msg += "iterative-match-df-max value lower."
                    # FIXME: This could be dealt with dynamically??
                    raise ValueError(err_msg)

                # if the result is a really bad match, trust it isn't
                # misrepresenting a good match
                if (1 - match) > 0.05 + (1 - min_match):
                    break

                # calculation converged
                if match_last > 0 and abs(match_last - match) < 0.001:
                    break

                # otherwise, refine calculation
                match_last = match
                df /= 2.0

            if match > min_match:
                return (match, tmplt)

            # record match and template params for highest match
            if match > max_match:
                max_match = match
                template = tmplt

        return (max_match, template)

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
        matches = [self.compute_match(tmplt, proposal, df, ASD=ASD) for tmplt in tmpbank]
        best_tmplt_ind = np.argmax(matches)
        self._nmatch += len(tmpbank)

        # best_tmplt_ind indexes tmpbank; add low to get index of full bank
        return matches[best_tmplt_ind], best_tmplt_ind + low

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
