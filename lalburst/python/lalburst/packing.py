# Copyright (C) 2006  Kipp Cannon
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
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


"""
This module provides bin packing utilities.
"""


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
from .git_version import date as __date__
from .git_version import version as __version__


#
# =============================================================================
#
#                                     Bins
#
# =============================================================================
#


# FIXME:  this should probably be a subclass of list, but no doubt that will
# break things so make this change some time it's convenient.


class Bin(object):
	"""
	Bin object for use in packing algorithm implementations.  A Bin
	instance has two attributes:  size, which is the total "size" of
	the contents of the Bin, and objects, which is a list of the Bins
	contents.

	Example:

	>>> strings = Bin()
	>>> s = "hello"
	>>> strings.add(s, len(s))
	>>> s.objects
	['hello']
	>>> s.size
	5
	"""
	def __init__(self):
		"""
		Initialize a new Bin instance.
		"""
		self.objects = []
		self.size = 0

	def add(self, obj, size):
		"""
		Add the object, whose size is as given, to the bin.
		"""
		self.objects.append(obj)
		self.size += size
		return self

	def __iadd__(self, other):
		"""
		Add the contents of another Bin object to this one.
		"""
		self.objects.extend(other.objects)
		self.size += other.size
		return self

	def __cmp__(self, other):
		"""
		Compare two Bin objects by their sizes.
		"""
		return cmp(self.size, other.size)

	def __repr__(self):
		"""
		A representation of the Bin object.
		"""
		return "Bin(size=%s, %s)" % (str(self.size), str(self.objects))

	__str__ = __repr__


#
# =============================================================================
#
#                              Packing Algorithms
#
# =============================================================================
#


# FIXME:  this should probably also be a subclass of list.

class Packer(object):
	"""
	Parent class for packing algorithms.  Specific packing algorithms
	should sub-class this, providing implementations of the pack() and
	packlist() methods.
	"""
	def __init__(self, bins):
		"""
		Set the list of bins on which we shall operate
		"""
		self.bins = bins

	def pack(self, size, obj):
		"""
		Pack an object of given size into the bins.
		"""
		raise NotImplementedError

	def packlist(self, size_object_pairs):
		"""
		Pack a list of (size, object) tuples into the bins.

		Default implementation invokes self.pack() on each pair in
		the the order given.
		"""
		for size, obj in size_object_pairs:
			self.pack(size, obj)


class BiggestIntoEmptiest(Packer):
	"""
	Packs the biggest object into the emptiest bin.
	"""
	def pack(self, size, obj):
		min(self.bins).add(obj, size)

	def packlist(self, size_object_pairs):
		for size, obj in sorted(size_object_pairs, reverse = True):
			self.pack(size, obj)
