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
Algorithms for low latency clustering of event-like datasets.
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"


import array

try:
	all
except NameError:
	# pre Python 2.5
	from glue.iterutils import all


class SlidingMaxHeap(object):
	"""A max heap that retains the order in which samples were added to it.
	The history capacity grows dynamically by expanding to the smallest power of
	two required to hold the contents.

	This is essentially an ordinary max heap, storing values and their sample
	offsets from the originating stream, coupled with a dynamically sized ring
	buffer storing the heap ranks of each sample.

	This procedure should take about O(1) time per input sample on average."""

	# ring buffer stores heap ranks
	# heap stores lists of the form [value, offset, item]

	def __init__(self, value = None):
		if value is None:
			value = lambda x: x
		self._value = value
		self._history = array.array('I', (0, 0))
		self._offset = 0
		self._heap = []

	def _drop_oldest(self):

		# Get heap index of oldest sample.
		old_index = self._history[(self._offset - len(self._heap)) % len(self._history)]

		# Get the heap index of the lower-right leaf of the heap.
		new_index = len(self._heap) - 1

		# _swap the lower-right leaf with the oldest sample, and restore the heap invariant.
		if new_index != old_index:
			self._swap(old_index, new_index)
			self._sift(old_index)

		# Remove the oldest sample, which is now at the lower-right leaf.
		self._heap.pop()

	def drop_while(self, predicate):
		while len(self._heap) > 0 and predicate(self._heap[self._history[(self._offset - len(self._heap)) % len(self._history)]][2]):
			self._drop_oldest()

	def append(self, x):
		# If adding this element would overwrite the tail of the ring buffer, then
		# expand the ring buffer.
		if len(self._heap) == len(self._history):

			end_index = self._offset % len(self._history)
			begin_index = end_index - len(self._heap)

			if begin_index >= 0:
				old_history = self._history[begin_index:end_index]
			else:
				old_history = self._history[begin_index:] + self._history[:end_index]

			# Allocate new ring buffer with twice the capacity as before.
			new_capacity = len(self._history) << 1
			self._history = array.array('I', (0,) * new_capacity)

			end_index = self._offset % len(self._history)
			begin_index = end_index - len(self._heap)
			if begin_index >= 0:
				self._history[begin_index:end_index] = old_history
			else:
				self._history[begin_index:] = old_history[:-begin_index]
				self._history[:end_index] = old_history[-end_index:]

		# Add the new datum to the ring buffer and the heap.
		k = len(self._heap)
		self._history[self._offset % len(self._history)] = k
		self._heap.append( (self._value(x), self._offset, x) )
		self._offset += 1
		self._sift_up(k)

	@property
	def root(self):
		"""Return the maximum element in the history (the 'root' of the binary
		tree that is the heap)."""
		return self._heap[0][2]

	@property
	def history(self):
		end_index = self._offset % len(self._history)
		begin_index = end_index - len(self._heap)

		if begin_index >= 0:
			old_history = self._history[begin_index:end_index]
		else:
			old_history = self._history[begin_index % len(self._history):] + self._history[:end_index]

		return old_history

	@property
	def count(self):
		return len(self._heap)

	def is_heap(self):
		return all(self._cmp(i, (i - 1) / 2) == -1 for i in range(1, len(self._heap)))

	def _swap(self, i, j):
		offset_i = self._heap[i][1] % len(self._history)
		offset_j = self._heap[j][1] % len(self._history)
		self._history[offset_i], self._history[offset_j] = j, i
		self._heap[i], self._heap[j] = self._heap[j], self._heap[i]

	def _cmp(self, i, j):
		return cmp(self._heap[i][0], self._heap[j][0])

	def _sift_up(self, i):
		while i > 0:
			j = (i - 1) / 2
			if self._cmp(i, j) == -1:
				break
			self._swap(i, j)
			i = j

	def _sift_down(self, i):
		n = len(self._heap)
		while True:
			j = 2 * i + 1
			k = j + 1
			if k < n - 1:
				if self._cmp(i, j) == -1:
					if self._cmp(j, k) == -1:
						self._swap(i, k)
						i = k
					else:
						self._swap(i, j)
						i = j
				elif self._cmp(i, k) == -1:
					self._swap(i, k)
					i = k
				else:
					break
			elif j < n - 1 and self._cmp(i, j) == -1:
				self._swap(i, j)
				i = j
			else:
				break

	def _sift(self, i):
		if i == 0:
			self._sift_down(i)
		else:
			parent = (i - 1) / 2
			if self._heap[i] > self._heap[parent]:
				self._sift_up(i)
			else:
				self._sift_down(i)


class ClusterBuilder(object):
	"""Provides one-dimensional time symmetric clustering.  Every input that is
	greater than all samples that come before or after it within a certain window
	is reported.

	key is a function that, when applied to the items being clustered, yields a
	monotonically increasing quantity that is treated as the 'time' of the item.

	value is a function that, when applied to the items being clustered, is the
	quantity that is maximized -- for example, SNR.

	Note: This implementation is tuned for event-like datasets in which sample
	times are widely and irregularly spaced.  For dense, regularly sampled data,
	there is a conceptually identical but practically much simpler algorithm.

	The cluster() and flush() routines both return generators.  Both may be
	called multiple times.  History is retained between calls."""

	def __init__(self, key, value, window):
		self._key = key
		self._delta = window
		self._2delta = 2 * window
		self._heap = SlidingMaxHeap(lambda x: value(x[1]))
		self._min_key = None

	def flush(self, max_key):
		"""Proceed with clustering items in the history under the assumption that
		there are no further items with a key less than max_key.  This method
		is used internally, but the user may call this as well if information
		about the absence of triggers is available independently from the
		triggers.  This may be useful in low latency applications."""
		if self._min_key is None:
			self._min_key = max_key
			return

		_2delta = self._2delta
		_delta = self._delta
		min_key = self._min_key

		if max_key < min_key:
			raise ValueError("Input keys are not increasing monotonically.")

		while max_key - min_key >= _2delta:
			if self._heap.count > 0:
				root_key, root = self._heap.root
				if root_key > max_key - _delta:
					min_key = root_key - _delta
				elif root_key < min_key + _delta:
					min_key = root_key
				else:
					yield root
					min_key = root_key
				self._heap.drop_while(lambda x: x[0] <= min_key)
			else:
				min_key = max_key - _delta

		self._min_key = min_key

	def cluster(self, iterable):
		key = self._key

		for item in iterable:
			max_key = key(item)
			for clustered_item in self.flush(max_key):
				yield clustered_item
			self._heap.append((max_key, item))


def clustered(iterable, key, value, window):
	"""Utility function for clustering an iterable without having to instantiate
	a ClusterBuilder object."""
	cluster_builder = ClusterBuilder(key, value, window)
	for clustered_item in cluster_builder.cluster(iterable):
		yield clustered_item
