# Copyright (C) 2010 Nickolas Fotopoulos
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

__author__ = "Nickolas Fotopoulos <nickolas.fotopoulos@ligo.org>"

import collections

import numpy as np

from pylal import rate

class MovingHistogramFixedN(object):
	"""
	This class is a histogram with a bounded size such that as new
	elements come in, old elements are removed in chronological order.
	The timestamp is retained only so that some rough idea of live time
	is retained and so that we can enforce monotonicity.
	"""
	def __init__(self, bins, max_len):
		self.bins = bins
		self.max_len = max_len

		self.timestamps = collections.deque()
		self.counts = np.zeros(len(bins), dtype=int)
		self._hist_ind = collections.deque()
		if hasattr(bins, "volumes"):  # 1-D Bins don't have volumes() methods
			self._volumes = bins.volumes()
		else:
			self._volumes = bins.upper() - bins.lower()

	def __len__(self):
		"""
		Return the current number of elements in the histogram.
		"""
		return len(self._hist_ind)

	def _flush_oldest(self):
		"""
		Remove the oldest element from the histogram.
		"""
		self.counts[self._hist_ind.popleft()] -= 1
		self.timestamps.popleft() # FIXME: deques can work as ring buffers in Python >= 2.6

	def update(self, timestamp, stat):
		"""
		Push the stat's bin onto the queue and update the histogram.
		"""
		if len(self) and timestamp < self.timestamps[-1]:
			raise ValueError, "timestamp non-monotonic: %s" % str(timestamp)
		ind = self.bins[stat]
		self.counts[ind] += 1
		while len(self) >= self.max_len:  # get down to max - 1
			self._flush_oldest()
		self._hist_ind.append(ind)
		self.timestamps.append(timestamp)

	def get_sf(self, stat):
		"""
		Return the fraction of events with stats than stat. This is formally
		the survival function (sf), which is 1 - CDF(stat).
		"""
		# FIXME: This may by slow with a deque. Must profile.
		ind = self.bins[stat]
		return self.counts[ind:].sum() / len(self)

	def get_pdf(self, stat):
		"""
		Return the PDF at stat.
		"""
		ind = self.bins[stat]
		return float(self.counts[ind]) / self._volumes[ind] / len(self)

	def get_cdf(self, stat):
		"""
		Return the CDF of the stat. (To get the "CDF from the right", see get_sf().)
		"""
		# FIXME: This may by slow with a deque. Must profile.
		ind = self.bins[stat]
		return self.counts[:ind + 1].sum() / len(self)


	def get_oldest_timestamp(self):
		return self.timestamps[0]
