#  Copyright (C) 2013 Chris Pankow
#
#  This program is free software; ynu can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Fnundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with with program; see the file COPYING. If not, write to the
#  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#  MA  02111-1307  USA


import sys

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot
# FIXME: gridspec...
#import breakout

## \addtogroup pkg_py_laldetchar_hveto
"""hierarchical veto utilities"""
# \author Chris Pankow (<chris.pankow@ligo.org>)
# \heading{Synopsis}
# ~~~
# from laldetchar.hveto import hveto_utils
# ~~~

import sys
import re
import json
import itertools
from collections import defaultdict

class HVetoResult(object):
	@staticmethod
	def make_key(n, snr_t, wind):
		return str((n, snr_t, wind))

	@staticmethod
	def decode_key(kstr):
		m = re.match(HVetoResult.ROUND_KEY, kstr)
		try:
			return (int(m.group(1)), float(m.group(2)), float(m.group(3)))
		except IndexError, AttributeError:
			raise "Could not parse round values from key string."

	def __init__(self):
		# Yep. That is in fact a dictionary of dictionaries that return 0 by default. Confused yet?
		self.counts = defaultdict(dict)
		self.coincs = defaultdict(dict)

	def to_json(self, fobj):
		json.dump({"count": self.counts, "coinc": self.coincs}, fobj)

	def from_json(self, fobj):
		tmpdata = json.load(fobj)
		self.coincs.update(tmpdata["coinc"])
		self.counts.update(tmpdata["count"])

	def add_round(self, n, snr_t, wind, chan, chancount, chancoinc):
		"""
		Add a coinc and count matrix (a dictionary, keyed by channel to integer values, with a key corresponding to the round number (n), SNR threshold (snr_t) and time window (wind).
		"""
		self.counts[HVetoResult.make_key(n, snr_t, wind)][chan] = chancount
		self.coincs[HVetoResult.make_key(n, snr_t, wind)][chan] = chancoinc
	def get_keys(self):
		return set(self.counts.keys() + self.coincs.keys())

	def get_chan_list(self, n, snr_t, wind):
		"""
		Returns the union of the channels present in the internal coincidence and count dictionaries. The result is sorted so to preserve order.
		FIXME: This will only get the channels that are defined in the rows, not columns.
		"""
		return sorted(set(self.counts[HVetoResult.make_key(n, snr_t, wind)].keys() + self.coincs[HVetoResult.make_key(n, snr_t, wind)].keys()))

	def get_round_info(self, n, snr_t, wind, chan=None, crosschan=None, transpose=False):
		"""
		Retrieve information from a specific correlation matrix set, keyed by round number (n), SNR threshold (snr_t), and time window (wind). If neither chan nor crosschan is specified, then the entire matrix will be returned. If only chan is specified, a row will be returned. If both chan and crosschan are given, then the entry for that matrix will be returned. Optionally, setting transpose to True will couple the count and coinc values, e.g. they will be returned as a set, rather than independently.
		NOTE: Uses self.get_chan_list internally to retrieve values, so deterministic order is guaranteed.
		"""
		subkey = HVetoResult.make_key(n, snr_t, wind)
		count, coinc = self.counts[subkey], self.coincs[subkey]
		if count is None or coinc is None:
			raise ValueError("Could not retrieve count or coinc for round specifier %d %f %f" % (n, snr_t, wind))
		if chan is not None:
			if count.has_key(chan):
				if crosschan is not None and count[chan].has_key(crosschan):
					count = count[chan][crosschan]
				elif crosschan is not None:
					# We don't have a value for this
					count = 0
				else:
					count = [count[chan][cchan] for cchan in sorted(count[chan].keys())]
			if coinc.has_key(chan):
				if crosschan is not None and coinc[chan].has_key(crosschan):
					coinc = coinc[chan][crosschan]
				elif crosschan is not None:
					# We don't have a value for this
					coinc = 0
				else:
					coinc = [coinc[chan][cchan] for cchan in sorted(coinc[chan].keys())]
		else:
			chanlist = self.get_chan_list(n, snr_t, wind)
			count = [count[c1][c2] for c2 in sorted(count[c1].keys()) for c1 in chanlist]
			coinc = [coinc[c1][c2] for c2 in sorted(coinc[c1].keys()) for c1 in chanlist]

		if transpose:
			# Matrix?
			result = []
			if type(coinc) is list and type(coinc[0]) is list and type(count) is list and type(count[0]) is list:
				for c1row, c2row in zip(count, coinc):
					result.append(zip(c1row, c2ro2))
			# Row?
			elif type(coinc) is list and type(coinc[0]) is not list and type(count) is list and type(count[0]) is not list:
				result = zip(count, coinc)
			# Values
			else:
				result = (count, coinc)
			return result

		return count, coinc

# close doxygen
##
#	\defgroup	pkg_py_laldetchar_triggers_hveto_utils   HVeto Utils
#@}
