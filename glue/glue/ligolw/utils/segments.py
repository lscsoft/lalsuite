# Copyright (C) 2008-2010,2012-2015  Kipp Cannon
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
Ask Kipp to document this!
"""

from lal import LIGOTimeGPS

from glue import git_version
from glue import iterutils
from glue import segments
from glue import segmentsUtils
from .. import lsctables


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                                 Segment List
#
# =============================================================================
#


class LigolwSegmentList(object):
	"""
	A description of a LIGO Light-Weight XML segment list.  Instances
	of this class carry all the metadata associated with a LIGO Light-
	Weight XML segment list including its name, version number, a
	comment, and so on.

	LIGO Light-Weight XML segment lists are three-state objects.  A
	segment list can be on, off, or undefined.  Two separate sequences
	of segments are used for this:  the "valid" list defines the
	intervals when the state of the segment list is known, and the
	"active" list defines the intervals when the segment list is on.
	It is not an error for the active list to be on during times when
	the segment lists state is unknown, this code does not impose any
	policy in that regard, but it should be expected that consumers of
	the segment list will treat all times when the segment list's state
	is unknown the same way.

	Example:

	>>> from glue.segments import *
	>>> segs = segmentlist([segment(0, 10), segment(20, 30)])
	>>> validity = segmentlist([segment(0, 10), segment(25, 100)])
	>>> x = LigolwSegmentList(active = segs, valid = validity, instruments = set(("H1",)), name = "test")
	>>> # x made copies of arguments
	>>> del segs[:]
	>>> segs
	[]
	>>> x.active
	[segment(0, 10), segment(20, 30)]
	>>> # some typical operations
	>>> x.active & x.valid	# known true
	[segment(0, 10), segment(25, 30)]
	>>> ~x.active & x.valid	# known false
	[segment(30, 100)]
	>>> x.active & ~x.valid	# not an error for this to be non-null
	[segment(20, 25)]
	>>> # make a copy
	>>> y = LigolwSegmentList(x)
	>>> del y.active[:]
	>>> y.active
	[]
	>>> x.active
	[segment(0, 10), segment(20, 30)]

	The arithmetic operators on this class implement Kleene's strong
	ternary logic, taking "true" to be (active & valid), "false" to be
	(~active & valid), and "unknown" to be ~valid.

	Example:

	>>> from glue.segments import *
	>>> segs = segmentlist([segment(0, 10), segment(20, 30)])
	>>> validity = segmentlist([segment(0, 35)])
	>>> x = LigolwSegmentList(active = segs, valid = validity, instruments = set(("H1",)), name = "test")
	>>> segs = segmentlist([segment(40, 50), segment(60, 70)])
	>>> validity = segmentlist([segment(35, 100)])
	>>> y = LigolwSegmentList(active = segs, valid = validity, instruments = set(("H1",)), name = "test")
	>>> (x | y).active
	[segment(0, 10), segment(20, 30), segment(40, 50), segment(60, 70)]
	>>> (x | y).valid
	[segment(0, 10), segment(20, 30), segment(40, 50), segment(60, 70)]
	>>> (x & y).active
	[]
	>>> (x & y).valid
	[segment(10, 20), segment(30, 40), segment(50, 60), segment(70, 100)]
	>>> (~x).active
	[segment(10, 20), segment(30, 35)]
	>>> (~x).valid
	[segment(0, 35)]

	With ternary logic the three basic Boolean operations AND, OR, and
	NOT, do not form a complete set of operations.  That is, there
	exist algebraic functions that cannot be implemented using
	combinations of these three operators alone.  One additional
	operator is required to construct a complete basis of logic
	operations, and we provide one:  .isfalse().  This operation
	inverts intervals of known state, and maps intervals of unknown
	state to false.

	>>> x.isfalse().active
	[segment(10, 20), segment(30, 35)]
	>>> x.isfalse().valid
	[segment(-infinity, infinity)]

	Unfortunately, one example of a function that cannot be constructed
	from the three basic Boolean operators is perhaps the most common
	operation we wish to perform with our tri-state segment lists.
	Often we wish to construct a tri-state list from two tri-state
	lists such that the final list's interval of validity is the union
	of the intervals of validity of the two source lists, and the state
	of the final list in that interval is the union the states of the
	source lists in that interval.  For example if from one source we
	know the state of some process spanning some time, and from another
	source we know the state of the same process spanning some other
	time, taken together we know the state of that process over the
	union of those times.  This function is given by

	>>> z = ~(x.isfalse() | y.isfalse() | (x & ~x & y & ~y))
	>>> z.active
	[segment(0, 10), segment(20, 30), segment(40, 50), segment(60, 70)]
	>>> z.valid
	[segment(0, 100)]

	Because this is inconvenient to type, slow, and not readable, a
	special in-place arithmetic operation named .update() is provided
	to implement this operation.

	>>> z = LigolwSegmentList(x).update(y)
	>>> z.active
	[segment(0, 10), segment(20, 30), segment(40, 50), segment(60, 70)]
	>>> z.valid
	[segment(0, 100)]

	The .update() method is not exactly equivalent to the operation
	above.  The .update() method demands that the two input lists'
	states be identical where their intervals of validity intersect.
	"""
	#
	# the columns the segment_definer, segment_summary and segment
	# tables need to have
	#

	segment_def_columns = (u"process_id", u"segment_def_id", u"ifos", u"name", u"version", u"comment")
	segment_sum_columns = (u"process_id", u"segment_sum_id", u"start_time", u"start_time_ns", u"end_time", u"end_time_ns", u"segment_def_id", u"comment")
	segment_columns = (u"process_id", u"segment_id", u"start_time", u"start_time_ns", u"end_time", u"end_time_ns", u"segment_def_id")

	def __init__(self, active = (), valid = (), instruments = (), name = None, version = None, comment = None):
		"""
		Initialize a new LigolwSegmentList instance.  active and
		valid are sequences that will be cast to
		segments.segmentlist objects.  They can be generator
		expressions.  The "active" sequence is what is usually
		thought of as the segment list, the "valid" sequence
		identifies the intervals of time for which the segment
		list's state is defined.
		"""
		# if we've only been passed an argument for active, see if
		# it's an object with the same attributes as ourselves and
		# if so initialize ourself as a copy of it.
		if not valid and not instruments and name is None and version is None and comment is None:
			try:
				self.valid = segments.segmentlist(active.valid)
				self.active = segments.segmentlist(active.active)
				self.instruments = set(active.instruments)
				self.name = active.name
				self.version = active.version
				self.comment = active.comment
				return
			except AttributeError:
				pass
		# we had more than one argument or it didn't have the
		# correct attributes, so do a normal initialization.  make
		# copies of mutable objects to avoid confusion
		self.valid = segments.segmentlist(valid)
		self.active = segments.segmentlist(active)
		self.instruments = set(instruments)
		self.name = name
		self.version = version
		self.comment = comment

	def sort(self, *args):
		"""
		Sort the internal segment lists.  The optional args are
		passed to the .sort() method of the segment lists.  This
		can be used to control the sort order by providing an
		alternate comparison function.  The default is to sort by
		start time with ties broken by end time.
		"""
		self.valid.sort(*args)
		self.active.sort(*args)

	def coalesce(self):
		"""
		Coalesce the internal segment lists.  Returns self.
		"""
		self.valid.coalesce()
		self.active.coalesce()
		return self

	def __ior__(self, other):
		"""
		If either is true the result is true, if both are false the
		result is false, otherwise the result is unknown.
		"""
		if self.instruments != other.instruments:
			raise ValueError("incompatible metadata")
		# times when one or the other is known to be true
		self.active &= self.valid
		self.active |= other.active & other.valid
		# times when result is defined
		self.valid = (self.valid & other.valid) | self.active
		return self

	def __iand__(self, other):
		"""
		If either is false the result is false, if both are true
		the result is true, otherwise the result is unknown.
		"""
		if self.instruments != other.instruments:
			raise ValueError("incompatible metadata")
		# times when one or the other is known to be false
		false = (self.valid & ~self.active) | (other.valid & ~other.active)
		# times when both are known to be true
		self.active &= self.valid
		self.active &= other.active & other.valid
		# times when result is defined
		self.valid = false | self.active
		return self

	def __or__(self, other):
		"""
		If either is true the result is true, if both are false the
		result is false, otherwise the result is unknown.
		"""
		result = type(self)(self)
		result |= other
		return result

	def __and__(self, other):
		"""
		If either is false the result is false, if both are true
		the result is true, otherwise the result is unknown.
		"""
		result = type(self)(self)
		result &= other
		return result

	def __invert__(self):
		"""
		If unknown the result is unknown, otherwise the state is
		inverted.
		"""
		result = type(self)(self)
		result.active = ~result.active & result.valid
		return result

	def isfalse(self):
		"""
		If unknown the result is false, otherwise the state is
		inverted.
		"""
		result = type(self)(self)
		result.active = ~result.active & result.valid
		result.valid = segments.segmentlist([segments.segment(-segments.infinity(), +segments.infinity())])
		return result

	def update(self, other):
		if self.instruments != other.instruments:
			raise ValueError("incompatible metadata")
		if (self.valid & other.valid).intersects(self.active ^ other.active):
			raise ValueError("result over-determined")
		self.active &= self.valid
		self.active |= other.active & other.valid
		self.valid |= other.valid
		return self


#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#


class LigolwSegments(set):
	"""
	An interface shim between code that makes use of segments in
	glue.segments form, and LIGO Light-Weight XML I/O code.

	This class is "attached" to an XML document object, at which time
	it parses and extracts the segment lists from the document, and
	clears the document's segment tables (preventing a second
	LigolwSegments object from being meaningfully attached to the same
	document).  When the application is finished manipulating the
	segment lists, they can be inserted into the XML document at which
	time the contents of the LigolwSegments object are cleared
	(preventing any further manipulations).

	This class is a subclass of the Python set builtin.  Each element
	of the set is a LigolwSegmentList instance describing one of the
	segment lists in the original XML document.

	This class may be used as a context manager to automate the
	replacement of segments back into the XML document, including in
	the event of an untrapped exception.  When used as a context
	manager, the process parameter of the .__init__() method is not
	optional.

	Example:

	>>> import sys
	>>> from glue.segments import *
	>>> from lal import LIGOTimeGPS
	>>> from glue.ligolw import ligolw, lsctables
	>>> xmldoc = ligolw.Document()
	>>> xmldoc.appendChild(ligolw.LIGO_LW())	# doctest: +ELLIPSIS
	<glue.ligolw.ligolw.LIGO_LW object at ...>
	>>> process = lsctables.Process()
	>>> process.process_id = lsctables.ProcessTable.get_next_id()
	>>> with LigolwSegments(xmldoc, process) as xmlsegments:
	...	h1segs = segmentlist([segment(LIGOTimeGPS(0), LIGOTimeGPS(10))])
	...	xmlsegments.insert_from_segmentlistdict({"H1": h1segs}, "test")
	>>> xmldoc.write(sys.stdout)		# doctest: +NORMALIZE_WHITESPACE
	<?xml version='1.0' encoding='utf-8'?>
	<!DOCTYPE LIGO_LW SYSTEM "http://ldas-sw.ligo.caltech.edu/doc/ligolwAPI/html/ligolw_dtd.txt">
	<LIGO_LW>
		<Table Name="segment_definer:table">
			<Column Type="ilwd:char" Name="segment_definer:process_id"/>
			<Column Type="ilwd:char" Name="segment_definer:segment_def_id"/>
			<Column Type="lstring" Name="segment_definer:ifos"/>
			<Column Type="lstring" Name="segment_definer:name"/>
			<Column Type="int_4s" Name="segment_definer:version"/>
			<Column Type="lstring" Name="segment_definer:comment"/>
			<Stream Delimiter="," Type="Local" Name="segment_definer:table">
				"process:process_id:0","segment_definer:segment_def_id:0","H1","test",,,
			</Stream>
		</Table>
		<Table Name="segment_summary:table">
			<Column Type="ilwd:char" Name="segment_summary:process_id"/>
			<Column Type="ilwd:char" Name="segment_summary:segment_sum_id"/>
			<Column Type="int_4s" Name="segment_summary:start_time"/>
			<Column Type="int_4s" Name="segment_summary:start_time_ns"/>
			<Column Type="int_4s" Name="segment_summary:end_time"/>
			<Column Type="int_4s" Name="segment_summary:end_time_ns"/>
			<Column Type="ilwd:char" Name="segment_summary:segment_def_id"/>
			<Column Type="lstring" Name="segment_summary:comment"/>
			<Stream Delimiter="," Type="Local" Name="segment_summary:table">
			</Stream>
		</Table>
		<Table Name="segment:table">
			<Column Type="ilwd:char" Name="segment:process_id"/>
			<Column Type="ilwd:char" Name="segment:segment_id"/>
			<Column Type="int_4s" Name="segment:start_time"/>
			<Column Type="int_4s" Name="segment:start_time_ns"/>
			<Column Type="int_4s" Name="segment:end_time"/>
			<Column Type="int_4s" Name="segment:end_time_ns"/>
			<Column Type="ilwd:char" Name="segment:segment_def_id"/>
			<Stream Delimiter="," Type="Local" Name="segment:table">
				"process:process_id:0","segment:segment_id:0",0,0,10,0,"segment_definer:segment_def_id:0"
			</Stream>
		</Table>
	</LIGO_LW>
	>>> xmlsegments = LigolwSegments(xmldoc)
	>>> xmlsegments.get_by_name("test")
	{u'H1': [segment(0.000000000, 10.000000000)]}
	>>> xmlsegments.get_by_name("wrong name")
	Traceback (most recent call last):
		...
	KeyError: "no segmentlists named 'wrong name'"

	NOTE:  the process of extracting and re-inserting the contents of
	the segment tables will, in general, randomize the IDs assigned to
	the rows of these tables.  If there are references to segment,
	segment_summary, or segment_definer row IDs in other tables in the
	document, those references will be broken by this process.
	"""
	def __init__(self, xmldoc, process = None):
		#
		# Find tables
		#

		try:
			self.segment_def_table = lsctables.SegmentDefTable.get_table(xmldoc)
		except ValueError:
			self.segment_def_table = xmldoc.childNodes[0].appendChild(lsctables.New(lsctables.SegmentDefTable, LigolwSegmentList.segment_def_columns))

		try:
			self.segment_sum_table = lsctables.SegmentSumTable.get_table(xmldoc)
		except ValueError:
			self.segment_sum_table = xmldoc.childNodes[0].appendChild(lsctables.New(lsctables.SegmentSumTable, LigolwSegmentList.segment_sum_columns))

		try:
			self.segment_table = lsctables.SegmentTable.get_table(xmldoc)
		except ValueError:
			self.segment_table = xmldoc.childNodes[0].appendChild(lsctables.New(lsctables.SegmentTable, LigolwSegmentList.segment_columns))

		#
		# Transform segment tables into a collection of
		# LigolwSegmentList objects for more convenient
		# manipulation
		#

		# construct empty LigolwSegmentList objects, one for each
		# entry in the segment_definer table, indexed by
		# segment_definer id
		segment_lists = dict((row.segment_def_id, LigolwSegmentList(instruments = row.instruments, name = row.name, version = row.version, comment = row.comment)) for row in self.segment_def_table)
		if len(segment_lists) != len(self.segment_def_table):
			raise ValueError("duplicate segment_def_id in segment_definer table")
		del self.segment_def_table[:]

		# populate LigolwSegmentList objects from segment table and
		# segment_summary table
		for row in self.segment_sum_table:
			try:
				segment_lists[row.segment_def_id].valid.append(row.segment)
			except KeyError as e:
				raise ValueError("invalid segment_def_id in segment_summary table: %s" % e)
		del self.segment_sum_table[:]
		for row in self.segment_table:
			try:
				segment_lists[row.segment_def_id].active.append(row.segment)
			except KeyError as e:
				raise ValueError("invalid segment_def_id in segment table: %s" % e)
		del self.segment_table[:]

		#
		# transcribe LigolwSegmentList objects into self.  segment
		# definer IDs are now meaningless
		#

		self.update(segment_lists.values())

		#
		# reset ID generators
		#

		self.segment_def_table.set_next_id(type(self.segment_def_table.next_id)(0))
		self.segment_table.set_next_id(type(self.segment_table.next_id)(0))
		self.segment_sum_table.set_next_id(type(self.segment_sum_table.next_id)(0))

		#
		# Save process row for later
		#

		self.process = process

		#
		# Done
		#


	def insert_from_segwizard(self, fileobj, instruments, name, version = None, comment = None):
		"""
		Parse the contents of the file object fileobj as a
		segwizard-format segment list, and insert the result as a
		new list of "active" segments into this LigolwSegments
		object.  A new entry will be created in the segment_definer
		table for the segment list, and instruments, name and
		comment are used to populate the entry's metadata.  Note
		that the "valid" segments are left empty, nominally
		indicating that there are no periods of validity.  Returns
		the newly created LigolwSegmentList object.
		"""
		ligolw_segment_list = LigolwSegmentList(active = segmentsUtils.fromsegwizard(fileobj, coltype = LIGOTimeGPS), instruments = instruments, name = name, version = version, comment = comment)
		self.add(ligolw_segment_list)
		return ligolw_segment_list


	def insert_from_segmentlistdict(self, seglists, name, version = None, comment = None):
		"""
		Insert the segments from the segmentlistdict object
		seglists as a new list of "active" segments into this
		LigolwSegments object.  The dictionary's keys are assumed
		to provide the instrument name for each segment list.  A
		new entry will be created in the segment_definer table for
		the segment lists, and the dictionary's keys, the name, and
		comment will be used to populate the entry's metadata.
		"""
		for instrument, segments in seglists.items():
			self.add(LigolwSegmentList(active = segments, instruments = set([instrument]), name = name, version = version, comment = comment))


	def coalesce(self):
		"""
		Coalesce the segment lists.  Returns self.
		"""
		for ligolw_segment_list in self:
			ligolw_segment_list.coalesce()
		return self


	def sort(self, *args):
		"""
		Sort the segment lists.  The optional args are passed to
		the .sort() methods of the segment lists.  This can be used
		to control the sort order by providing an alternate
		comparison function (the default is to sort all lists by
		segment start time with ties broken by end time).
		"""
		for ligolw_segment_list in self:
			ligolw_segment_list.sort(*args)


	def optimize(self):
		"""
		Identifies segment lists that differ only in their
		instruments --- they have the same valid and active
		segments, the same name, version and the same comment ---
		and then deletes all but one of them, leaving just a single
		list having the union of the instruments.
		"""
		self.sort()
		segment_lists = dict(enumerate(self))
		for target, source in [(idx_a, idx_b) for (idx_a, seglist_a), (idx_b, seglist_b) in iterutils.choices(segment_lists.items(), 2) if seglist_a.valid == seglist_b.valid and seglist_a.active == seglist_b.active and seglist_a.name == seglist_b.name and seglist_a.version == seglist_b.version and seglist_a.comment == seglist_b.comment]:
			try:
				source = segment_lists.pop(source)
			except KeyError:
				continue
			segment_lists[target].instruments |= source.instruments
		self.clear()
		self.update(segment_lists.values())


	def get_by_name(self, name, clip_to_valid = False):
		"""
		Retrieve the active segmentlists whose name equals name.
		The result is a segmentlistdict indexed by instrument.  All
		segmentlist objects within it will be copies of the
		contents of this object, modifications will not affect the
		contents of this object.  If clip_to_valid is True then the
		segmentlists will be intersected with their respective
		intervals of validity, otherwise they will be the verbatim
		active segments.

		NOTE:  the intersection operation required by clip_to_valid
		will yield undefined results unless the active and valid
		segmentlist objects are coalesced.
		"""
		result = segments.segmentlistdict()
		for seglist in self:
			if seglist.name != name:
				continue
			segs = seglist.active
			if clip_to_valid:
				# do not use in-place intersection
				segs = segs & seglist.valid
			for instrument in seglist.instruments:
				if instrument in result:
					raise ValueError("multiple '%s' segmentlists for instrument '%s'" % (name, instrument))
				# make copy so that instruments do not
				# share a single list, so that subsequent
				# arithmetic operations do not corrupt the
				# wrong instrument's segments
				result[instrument] = segments.segmentlist(segs)
		if not result:
			raise KeyError("no segmentlists named '%s'" % name)
		return result


	def finalize(self, process_row = None):
		"""
		Restore the LigolwSegmentList objects to the XML tables in
		preparation for output.  All segments from all segment
		lists are inserted into the tables in time order, but this
		is NOT behaviour external applications should rely on.
		This is done simply in the belief that it might assist in
		constructing well balanced indexed databases from the
		resulting files.  If that proves not to be the case, or for
		some reason this behaviour proves inconvenient to preserve,
		then it might be discontinued without notice.  You've been
		warned.
		"""
		if process_row is not None:
			process_id = process_row.process_id
		elif self.process is not None:
			process_id = self.process.process_id
		else:
			raise ValueError("must supply a process row to .__init__()")

		#
		# ensure ID generators are synchronized with table contents
		#

		self.segment_def_table.sync_next_id()
		self.segment_table.sync_next_id()
		self.segment_sum_table.sync_next_id()

		#
		# put all segment lists in time order
		#

		self.sort()

		#
		# generator function to convert segments into row objects,
		# each paired with the table to which the row is to be
		# appended
		#

		def row_generator(segs, target_table, process_id, segment_def_id):
			id_column = target_table.next_id.column_name
			for seg in segs:
				row = target_table.RowType()
				row.segment = seg
				row.process_id = process_id
				row.segment_def_id = segment_def_id
				setattr(row, id_column, target_table.get_next_id())
				if hasattr(row, "comment"):
					row.comment = None
				yield row, target_table

		#
		# populate the segment_definer table from the list of
		# LigolwSegmentList objects and construct a matching list
		# of table row generators.  empty ourselves to prevent this
		# process from being repeated
		#

		row_generators = []
		while self:
			ligolw_segment_list = self.pop()
			segment_def_row = self.segment_def_table.RowType()
			segment_def_row.process_id = process_id
			segment_def_row.segment_def_id = self.segment_def_table.get_next_id()
			segment_def_row.instruments = ligolw_segment_list.instruments
			segment_def_row.name = ligolw_segment_list.name
			segment_def_row.version = ligolw_segment_list.version
			segment_def_row.comment = ligolw_segment_list.comment
			self.segment_def_table.append(segment_def_row)

			row_generators.append(row_generator(ligolw_segment_list.valid, self.segment_sum_table, process_id, segment_def_row.segment_def_id))
			row_generators.append(row_generator(ligolw_segment_list.active, self.segment_table, process_id, segment_def_row.segment_def_id))

		#
		# populate segment and segment_summary tables by pulling
		# rows from the generators in time order
		#

		for row, target_table in iterutils.inorder(*row_generators):
			target_table.append(row)


	def __enter__(self):
		if self.process is None:
			raise ValueError("must supply a process row to .__init__()")
		return self


	def __exit__(self, *args):
		self.finalize()
		return False


#
# =============================================================================
#
#                                  Utilities
#
# =============================================================================
#


def has_segment_tables(xmldoc, name = None):
	"""
	Return True if the document contains a complete set of segment
	tables.  Returns False otherwise.  If name is given and not None
	then the return value is True only if the document's segment
	tables, if present, contain a segment list by that name.
	"""
	try:
		names = lsctables.SegmentDefTable.get_table(xmldoc).getColumnByName("name")
		lsctables.SegmentTable.get_table(xmldoc)
		lsctables.SegmentSumTable.get_table(xmldoc)
	except (ValueError, KeyError):
		return False
	return name is None or name in names


def segmenttable_get_by_name(xmldoc, name):
	"""
	Retrieve the segmentlists whose name equals name.  The result is a
	segmentlistdict indexed by instrument.

	The output of this function is not coalesced, each segmentlist
	contains the segments as found in the segment table.

	NOTE:  this is a light-weight version of the .get_by_name() method
	of the LigolwSegments class intended for use when the full
	machinery of that class is not required.  Considerably less
	document validation and error checking is performed by this
	version.  Consider using that method instead if your application
	will be interfacing with the document via that class anyway.
	"""
	#
	# find required tables
	#

	def_table = lsctables.SegmentDefTable.get_table(xmldoc)
	seg_table = lsctables.SegmentTable.get_table(xmldoc)

	#
	# segment_def_id --> instrument names mapping but only for
	# segment_definer entries bearing the requested name
	#

	instrument_index = dict((row.segment_def_id, row.instruments) for row in def_table if row.name == name)

	#
	# populate result segmentlistdict object from segment_def_map table
	# and index
	#

	instruments = set(instrument for instruments in instrument_index.values() for instrument in instruments)
	result = segments.segmentlistdict((instrument, segments.segmentlist()) for instrument in instruments)

	for row in seg_table:
		if row.segment_def_id in instrument_index:
			seg = row.segment
			for instrument in instrument_index[row.segment_def_id]:
				result[instrument].append(seg)

	#
	# done
	#

	return result
