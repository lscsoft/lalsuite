#!/usr/bin/env python
# Copyright (C) 2010  Leo Singer
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
#
"""
Unit test suite for pylal.low_latency_clustering.
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"


import unittest
import numpy as np
from collections import deque
from pylal.low_latency_clustering import SlidingMaxHeap, clustered

try:
	all
except NameError:
	# pre Python 2.5
	from glue.iterutils import all


def all_distinct(iterable):
	l = list(iterable)
	s = set(l)
	return len(l) == len(s)


class TestHeapMonotonic(unittest.TestCase):

	def setUp(self):
		count_trials = 1000
		keys = np.cumsum(np.random.random_integers(1, 10000, count_trials))
		vals = range(count_trials)
		self.indata = zip(keys, vals)
		self.keyfunc = lambda x: x[0]
		self.valfunc = lambda x: x[1]
		self.heap = SlidingMaxHeap(value = self.valfunc)

	def test1(self):
		"""Insert monotonic time series, and then delete all elements."""
		for item in self.indata:
			self.heap.append(item)
			self.assertEqual(item, self.heap.root)
			self.assertTrue(self.heap.is_heap())
			self.assertTrue(all_distinct(self.heap.history))
		for item in self.indata:
			self.assertTrue(self.heap.is_heap())
			self.assertTrue(all_distinct(self.heap.history))
			self.heap._drop_oldest()

	def test2(self):
		"""Maintain fixed window history of a monotonic time series."""
		for i, item in enumerate(self.indata):
			self.heap.append(item)
			if i > 100:
				self.heap._drop_oldest()
			self.assertEqual(item, self.heap.root)
			self.assertTrue(self.heap.is_heap())


class TestHeapRandom(unittest.TestCase):

	@staticmethod
	def sliding_max_dumb(iterable, window_len, valfunc):
		"""Simpler, slower sliding max algorithm, but easier to verify correctness."""
		history = deque()
		for item in iterable:
			history.append(item)
			while len(history) > window_len:
				history.popleft()
			yield max((valfunc(x), x) for x in history)[1]

	@staticmethod
	def sliding_max(iterable, window_len, valfunc):
		heap = SlidingMaxHeap(value = valfunc)
		for i, item in enumerate(iterable):
			heap.append(item)
			heap.drop_while(lambda x: len(heap._heap) > window_len)
			yield heap.root

	def setUp(self):
		count_trials = 1000
		keys = np.cumsum(np.random.random_integers(1, 10000, count_trials))
		vals = np.random.rand(1, 100000, count_trials)
		self.indata = zip(keys, vals)
		self.keyfunc = lambda x: x[0]
		self.valfunc = lambda x: x[1]
		self.window_len = 50
		self.heap = SlidingMaxHeap(value = self.valfunc)

	def test(self):
		"""Test fixed length history sliding maximum on random input."""
		for output, expected_output in zip(self.sliding_max(self.indata, self.window_len, self.valfunc), self.sliding_max_dumb(self.indata, self.window_len, self.valfunc)):
			self.assertEqual(output, expected_output)


class TestCluster(unittest.TestCase):

	@staticmethod
	def cluster_dumb(list, window, key, value):
		"""The full O(N^2) clustering algorithm, but easier to verify correctness."""
		return [item for item in list if all(abs(key(item) - key(other)) >= window or value(item) >= value(other) for other in list)]

	def setUp(self):
		count_trials = 1000
		keys = np.cumsum(np.random.random_integers(1, 10, count_trials))
		vals = np.random.rand(count_trials)
		self.indata = zip(keys, vals)
		self.keyfunc = lambda x: x[0]
		self.valfunc = lambda x: x[1]
		self.window = 40
		latest_key = self.keyfunc(self.indata[-1]) - self.window
		self.required_outdata = filter(
			lambda x: self.window <= self.keyfunc(x) <= latest_key,
			self.cluster_dumb(self.indata, self.window, self.keyfunc, self.valfunc)
		)

	def test(self):
		"""Test one-dimensional clustering on random input."""
		clustered_items = list(clustered(self.indata, self.keyfunc, self.valfunc, self.window))
		self.assertEqual(clustered_items, self.required_outdata)


if __name__ == '__main__':
	suite = unittest.main()
