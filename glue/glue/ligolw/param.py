# Copyright (C) 2006--2009,2012--2016  Kipp Cannon
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
High-level support for Param elements.
"""


import re
import sys
from xml.sax.saxutils import escape as xmlescape
from xml.sax.xmlreader import AttributesImpl as Attributes
try:
	import yaml
except ImportError:
	# yaml serialization is optional
	pass


from glue import git_version
from . import ligolw
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


def getParamsByName(elem, name):
	"""
	Return a list of params with name name under elem.
	"""
	name = Param.ParamName(name)
	return elem.getElements(lambda e: (e.tagName == ligolw.Param.tagName) and (e.Name == name))


def get_param(xmldoc, name):
	"""
	Scan xmldoc for a param named name.  Raises ValueError if not
	exactly 1 such param is found.
	"""
	params = getParamsByName(xmldoc, name)
	if len(params) != 1:
		raise ValueError("document must contain exactly one %s param" % Param.ParamName(name))
	return params[0]


def get_pyvalue(xml, name):
	"""
	Convenience wrapper for get_param() that recovers an instance of a
	Python builtin type from a Param element.
	"""
	# Note:  the Param is automatically parsed into the correct Python
	# type, so this function is mostly a no-op.
	return get_param(xml, name).pcdata


#
# =============================================================================
#
#                               Element Classes
#
# =============================================================================
#


#
# FIXME: params of type string should be quoted in order to correctly
# delimit their extent.  If that were done, then the pcdata in a Param
# element could be parsed using the Stream tokenizer (i.e., as though it
# were a single-token stream), which would guarantee that Stream data and
# Param data is parsed using the exact same rules.  Unfortunately, common
# practice is to not quote Param string values, so we parse things
# differently here.  In particular, we strip whitespace from the start and
# stop of all Param pcdata.  If this causes your string Param values to be
# corrupted (because you need leading and trailing white space preserved),
# then you need to make everyone switch to quoting their string Param
# values, and once that is done then this code will be changed.  Perhaps a
# warning should be emitted for non-quoted strings to encourage a
# transition?
#


class Param(ligolw.Param):
	"""
	High-level Param element.  The value is stored in the pcdata
	attribute as the native Python type rather than as a string.
	"""
	class ParamName(ligolw.LLWNameAttr):
		dec_pattern = re.compile(r"(?P<Name>[a-z0-9_:]+):param\Z")
		enc_pattern = u"%s:param"

	Name = ligolw.attributeproxy(u"Name", enc = ParamName.enc, dec = ParamName)
	Scale = ligolw.attributeproxy(u"Scale", enc = ligolwtypes.FormatFunc[u"real_8"], dec = ligolwtypes.ToPyType[u"real_8"])
	Type = ligolw.attributeproxy(u"Type", default = u"lstring")

	def endElement(self):
		if self.pcdata is not None:
			# convert pcdata from string to native Python type
			if self.Type == u"yaml":
				try:
					yaml
				except NameError:
					raise NotImplementedError("yaml support not installed")
				self.pcdata = yaml.load(self.pcdata)
			else:
				self.pcdata = ligolwtypes.ToPyType[self.Type](self.pcdata.strip())

	def write(self, fileobj = sys.stdout, indent = u""):
		fileobj.write(self.start_tag(indent))
		for c in self.childNodes:
			if c.tagName not in self.validchildren:
				raise ligolw.ElementError("invalid child %s for %s" % (c.tagName, self.tagName))
			c.write(fileobj, indent + ligolw.Indent)
		if self.pcdata is not None:
			if self.Type == u"yaml":
				try:
					yaml
				except NameError:
					raise NotImplementedError("yaml support not installed")
				fileobj.write(xmlescape(yaml.dump(self.pcdata).strip()))
			else:
				# we have to strip quote characters from
				# string formats (see comment above)
				fileobj.write(xmlescape(ligolwtypes.FormatFunc[self.Type](self.pcdata).strip(u"\"")))
		fileobj.write(self.end_tag(u"") + u"\n")

	@classmethod
	def build(cls, name, Type, value, start = None, scale = None, unit = None, dataunit = None, comment = None):
		"""
		Construct a LIGO Light Weight XML Param document subtree.
		FIXME: document keyword arguments.
		"""
		elem = cls()
		elem.Name = name
		elem.Type = Type
		elem.pcdata = value
		# FIXME:  I have no idea how most of the attributes should be
		# encoded, I don't even know what they're supposed to be.
		if dataunit is not None:
			elem.DataUnit = dataunit
		if scale is not None:
			elem.Scale = scale
		if start is not None:
			elem.Start = start
		if unit is not None:
			elem.Unit = unit
		if comment is not None:
			elem.appendChild(ligolw.Comment()).pcdata = comment
		return elem

	@classmethod
	def from_pyvalue(cls, name, value, **kwargs):
		"""
		Convenience wrapper for .build() that constructs a Param
		element from an instance of a Python builtin type.  See
		.build() for a description of the valid keyword arguments.
		"""
		return cls.build(name, ligolwtypes.FromPyType[type(value)], value, **kwargs)


# FIXME: delete when nothing uses
new_param = Param.build
from_pyvalue = Param.from_pyvalue


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
	glue.ligolw.LIGOLWContentHandler, to cause it to use the Param
	class defined in this module when parsing XML documents.

	Example:

	>>> from glue.ligolw import ligolw
	>>> def MyContentHandler(ligolw.LIGOLWContentHandler):
	...	pass
	...
	>>> use_in(MyContentHandler)
	"""
	def startParam(self, parent, attrs):
		return Param(attrs)

	ContentHandler.startParam = startParam

	return ContentHandler
