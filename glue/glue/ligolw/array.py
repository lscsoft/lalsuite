# Copyright (C) 2006--2015  Kipp Cannon
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
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
While the ligolw module provides classes and parser support for reading and
writing LIGO Light Weight XML documents, this module supplements that code
with classes and parsers that add intelligence to the in-RAM document
representation.

In particular, the document tree associated with an Array element is
enhanced.  During parsing, the Stream element in this module converts the
character data contained within it into the elements of a numpy array
object.  The array has the appropriate dimensions and type.  When the
document is written out again, the Stream element serializes the array back
into character data.

The array is stored as an attribute of the Array element.
"""


import itertools
import numpy
import re
import sys
import warnings
from xml.sax.saxutils import escape as xmlescape
from xml.sax.xmlreader import AttributesImpl as Attributes


from glue import git_version
from . import ligolw
from . import tokenizer
from . import types as ligolwtypes


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                                  Utilities
#
# =============================================================================
#


def getArraysByName(elem, name):
	"""
	Return a list of arrays with name name under elem.
	"""
	name = Array.ArrayName(name)
	return elem.getElements(lambda e: (e.tagName == ligolw.Array.tagName) and (e.Name == name))


def get_array(xmldoc, name):
	"""
	Scan xmldoc for an array named name.  Raises ValueError if not
	exactly 1 such array is found.
	"""
	arrays = getArraysByName(xmldoc, name)
	if len(arrays) != 1:
		raise ValueError("document must contain exactly one %s array" % Array.ArrayName(name))
	return arrays[0]


#
# =============================================================================
#
#                               Element Classes
#
# =============================================================================
#


class ArrayStream(ligolw.Stream):
	"""
	High-level Stream element for use inside Arrays.  This element
	knows how to parse the delimited character stream into the parent's
	array attribute, and knows how to turn the parent's array attribute
	back into a character stream.
	"""

	Delimiter = ligolw.attributeproxy(u"Delimiter", default = u" ")

	def __init__(self, *args):
		super(ArrayStream, self).__init__(*args)
		try:
			self.Encoding
		except AttributeError:
			pass
		else:
			raise ligolw.ElementError("non-default encoding '%s' not supported.  if this is critical, please report." % self.Encoding)
		self._tokenizer = tokenizer.Tokenizer(self.Delimiter)

	def config(self, parentNode):
		# some initialization that can only be done once parentNode
		# has been set.
		self._tokenizer.set_types([ligolwtypes.ToPyType[parentNode.Type]])
		parentNode.array = numpy.zeros(parentNode.get_shape(), ligolwtypes.ToNumPyType[parentNode.Type])
		self._array_view = parentNode.array.T.flat
		self._index = 0
		return self

	def appendData(self, content):
		# tokenize buffer, and assign to array
		tokens = tuple(self._tokenizer.append(content))
		next_index = self._index + len(tokens)
		self._array_view[self._index : next_index] = tokens
		self._index = next_index

	def endElement(self):
		# stream tokenizer uses delimiter to identify end of each
		# token, so add a final delimiter to induce the last token
		# to get parsed.
		self.appendData(self.Delimiter)
		if self._index != len(self._array_view):
			raise ValueError("length of Stream (%d elements) does not match array size (%d elements)" % (self._index, len(self._array_view)))
		del self._array_view
		del self._index

	def write(self, fileobj = sys.stdout, indent = u""):
		# avoid symbol and attribute look-ups in inner loop
		linelen = self.parentNode.array.shape[0]
		lines = self.parentNode.array.size / linelen
		tokens = itertools.imap(ligolwtypes.FormatFunc[self.parentNode.Type], self.parentNode.array.T.flat)
		islice = itertools.islice
		join = self.Delimiter.join
		w = fileobj.write

		w(self.start_tag(indent))
		if lines:
			newline = u"\n" + indent + ligolw.Indent
			w(newline)
			w(xmlescape(join(islice(tokens, linelen))))
			newline = self.Delimiter + newline
			for i in xrange(lines - 1):
				w(newline)
				w(xmlescape(join(islice(tokens, linelen))))
		w(u"\n" + self.end_tag(indent) + u"\n")


class Array(ligolw.Array):
	"""
	High-level Array element.
	"""
	class ArrayName(ligolw.LLWNameAttr):
		dec_pattern = re.compile(r"(?P<Name>[a-zA-Z0-9_:]+):array\Z")
		enc_pattern = u"%s:array"

	Name = ligolw.attributeproxy(u"Name", enc = ArrayName.enc, dec = ArrayName)

	def __init__(self, *args):
		"""
		Initialize a new Array element.
		"""
		super(Array, self).__init__(*args)
		self.array = None

	def get_shape(self):
		"""
		Return a tuple of this array's dimensions.  This is done by
		querying the Dim children.  Note that once it has been
		created, it is also possible to examine an Array object's
		.array attribute directly, and doing that is much faster.
		"""
		return tuple(c.n for c in self.getElementsByTagName(ligolw.Dim.tagName))[::-1]

	@classmethod
	def build(cls, name, array, dim_names = None):
		"""
		Construct a LIGO Light Weight XML Array document subtree
		from a numpy array object.

		Example:

		>>> import numpy, sys
		>>> a = numpy.arange(12, dtype = "double")
		>>> a.shape = (4, 3)
		>>> from_array(u"test", a).write(sys.stdout)	# doctest: +NORMALIZE_WHITESPACE
		<Array Type="real_8" Name="test:array">
			<Dim>3</Dim>
			<Dim>4</Dim>
			<Stream Delimiter=" " Type="Local">
				0 3 6 9
				1 4 7 10
				2 5 8 11
			</Stream>
		</Array>
		"""
		# Type must be set for .__init__();  easier to set Name
		# afterwards to take advantage of encoding handled by
		# attribute proxy
		elem = cls(Attributes({u"Type": ligolwtypes.FromNumPyType[str(array.dtype)]}))
		elem.Name = name
		if dim_names is None:
			dim_names = [None] * len(array.shape)
		elif len(dim_names) != len(array.shape):
			raise ValueError("dim_names must be same length as number of dimensions")
		for name, n in reversed(zip(dim_names, array.shape)):
			child = elem.appendChild(ligolw.Dim())
			if name is not None:
				child.Name = name
			child.n = n
		elem.appendChild(ArrayStream(Attributes({u"Type": ArrayStream.Type.default, u"Delimiter": ArrayStream.Delimiter.default})))
		elem.array = array
		return elem

	#
	# Element methods
	#

	def unlink(self):
		"""
		Break internal references within the document tree rooted
		on this element to promote garbage collection.
		"""
		super(Array, self).unlink()
		self.array = None


# FIXME:  remove when no longer used
from_array = Array.build


#
# =============================================================================
#
#                               Content Handler
#
# =============================================================================
#


#
# Override portions of a ligolw.LIGOLWContentHandler class
#


def use_in(ContentHandler):
	"""
	Modify ContentHandler, a sub-class of
	glue.ligolw.LIGOLWContentHandler, to cause it to use the Array and
	ArrayStream classes defined in this module when parsing XML
	documents.

	Example:

	>>> from glue.ligolw import ligolw
	>>> class MyContentHandler(ligolw.LIGOLWContentHandler):
	...	pass
	...
	>>> use_in(MyContentHandler)
	<class 'glue.ligolw.array.MyContentHandler'>
	"""
	def startStream(self, parent, attrs, __orig_startStream = ContentHandler.startStream):
		if parent.tagName == ligolw.Array.tagName:
			return ArrayStream(attrs).config(parent)
		return __orig_startStream(self, parent, attrs)

	def startArray(self, parent, attrs):
		return Array(attrs)

	ContentHandler.startStream = startStream
	ContentHandler.startArray = startArray

	return ContentHandler
