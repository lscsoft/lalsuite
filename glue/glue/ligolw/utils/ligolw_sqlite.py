#
# Copyright (C) 2006-2014  Kipp Cannon
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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


"""
Convert tabular data in LIGO LW XML files to and from SQL databases.

This module provides a library interface to the machinery used by the
ligolw_sqlite command-line tool, facilitating it's re-use in other
applications.  The real XML<-->database translation machinery is
implemented in the glue.ligolw.dbtables module.  The code here wraps the
machinery in that mdoule in functions that are closer to the command-line
level operations provided by the ligolw_sqlite program.
"""


try:
	import sqlite3
except ImportError:
	# pre 2.5.x
	from pysqlite2 import dbapi2 as sqlite3
import sys


from glue import git_version
from .. import ligolw
from .. import dbtables
from .. import utils as ligolw_utils


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                                 Library Code
#
# =============================================================================
#


#
# Open database
#


def setup(target, check_same_thread = True):
	connection = sqlite3.connect(target, check_same_thread = check_same_thread)
	dbtables.idmap_sync(connection)
	return connection


#
# How to insert
#


def update_ids(connection, xmldoc, verbose = False):
	"""
	For internal use only.
	"""
	# NOTE:  it's critical that the xmldoc object be retrieved *before*
	# the rows whose IDs need to be updated are inserted.  The xml
	# retrieval resets the "last max row ID" values inside the table
	# objects, so if retrieval of the xmldoc is deferred until after
	# the rows are inserted, nothing will get updated.  therefore, the
	# connection and xmldoc need to be passed separately to this
	# function, even though it seems this function could reconstruct
	# the xmldoc itself from the connection.
	table_elems = xmldoc.getElementsByTagName(ligolw.Table.tagName)
	for i, tbl in enumerate(table_elems):
		if verbose:
			sys.stderr.write("updating IDs: %d%%\r" % (100.0 * i / len(table_elems)))
		tbl.applyKeyMapping()
	if verbose:
		sys.stderr.write("updating IDs: 100%\n")

	# reset ID mapping for next document
	dbtables.idmap_reset(connection)


def insert_from_url(url, preserve_ids = False, verbose = False, contenthandler = None):
	"""
	Parse and insert the LIGO Light Weight document at the URL into the
	database with which the content handler is associated.  If
	preserve_ids is False (default), then row IDs are modified during
	the insert process to prevent collisions with IDs already in the
	database.  If preserve_ids is True then IDs are not modified;  this
	will result in database consistency violations if any of the IDs of
	newly-inserted rows collide with row IDs already in the database,
	and is generally only sensible when inserting a document into an
	empty database.  If verbose is True then progress reports will be
	printed to stderr.  See glue.ligolw.dbtables.use_in() for more
	information about constructing a suitable content handler class.
	"""
	#
	# enable/disable ID remapping
	#

	orig_DBTable_append = dbtables.DBTable.append

	if not preserve_ids:
		try:
			dbtables.idmap_create(contenthandler.connection)
		except sqlite3.OperationalError:
			# assume table already exists
			pass
		dbtables.idmap_sync(contenthandler.connection)
		dbtables.DBTable.append = dbtables.DBTable._remapping_append
	else:
		dbtables.DBTable.append = dbtables.DBTable._append

	try:
		#
		# load document.  this process inserts the document's contents into
		# the database.  the XML tree constructed by this process contains
		# a table object for each table found in the newly-inserted
		# document and those table objects' last_max_rowid values have been
		# initialized prior to rows being inserted.  therefore, this is the
		# XML tree that must be passed to update_ids in order to ensure (a)
		# that all newly-inserted tables are processed and (b) all
		# newly-inserted rows are processed.  NOTE:  it is assumed the
		# content handler is creating DBTable instances in the XML tree,
		# not regular Table instances, but this is not checked.
		#

		xmldoc = ligolw_utils.load_url(url, verbose = verbose, contenthandler = contenthandler)

		#
		# update references to row IDs and cleanup ID remapping
		#

		if not preserve_ids:
			update_ids(contenthandler.connection, xmldoc, verbose = verbose)

	finally:
		dbtables.DBTable.append = orig_DBTable_append

	#
	# done.  unlink the document to delete database cursor objects it
	# retains
	#

	contenthandler.connection.commit()
	xmldoc.unlink()


def insert_from_xmldoc(connection, source_xmldoc, preserve_ids = False, verbose = False):
	"""
	Insert the tables from an in-ram XML document into the database at
	the given connection.  If preserve_ids is False (default), then row
	IDs are modified during the insert process to prevent collisions
	with IDs already in the database.  If preserve_ids is True then IDs
	are not modified;  this will result in database consistency
	violations if any of the IDs of newly-inserted rows collide with
	row IDs already in the database, and is generally only sensible
	when inserting a document into an empty database.  If verbose is
	True then progress reports will be printed to stderr.
	"""
	#
	# enable/disable ID remapping
	#

	orig_DBTable_append = dbtables.DBTable.append

	if not preserve_ids:
		try:
			dbtables.idmap_create(connection)
		except sqlite3.OperationalError:
			# assume table already exists
			pass
		dbtables.idmap_sync(connection)
		dbtables.DBTable.append = dbtables.DBTable._remapping_append
	else:
		dbtables.DBTable.append = dbtables.DBTable._append

	try:
		#
		# create a place-holder XML representation of the target
		# document so we can pass the correct tree to update_ids().
		# note that only tables present in the source document need
		# ID ramapping, so xmldoc only contains representations of
		# the tables in the target document that are also in the
		# source document
		#

		xmldoc = ligolw.Document()
		xmldoc.appendChild(ligolw.LIGO_LW())

		#
		# iterate over tables in the source XML tree, inserting
		# each into the target database
		#

		for tbl in source_xmldoc.getElementsByTagName(ligolw.Table.tagName):
			#
			# instantiate the correct table class, connected to the
			# target database, and save in XML tree
			#

			name = tbl.Name
			try:
				cls = dbtables.TableByName[name]
			except KeyError:
				cls = dbtables.DBTable
			dbtbl = xmldoc.childNodes[-1].appendChild(cls(tbl.attributes, connection = connection))

			#
			# copy table element child nodes from source XML tree
			#

			for elem in tbl.childNodes:
				if elem.tagName == ligolw.Stream.tagName:
					dbtbl._end_of_columns()
				dbtbl.appendChild(type(elem)(elem.attributes))

			#
			# copy table rows from source XML tree
			#

			for row in tbl:
				dbtbl.append(row)
			dbtbl._end_of_rows()

		#
		# update references to row IDs and clean up ID remapping
		#

		if not preserve_ids:
			update_ids(connection, xmldoc, verbose = verbose)

	finally:
		dbtables.DBTable.append = orig_DBTable_append

	#
	# done.  unlink the document to delete database cursor objects it
	# retains
	#

	connection.commit()
	xmldoc.unlink()


def insert_from_urls(urls, contenthandler, **kwargs):
	"""
	Iterate over a sequence of URLs, calling insert_from_url() on each,
	then build the indexes indicated by the metadata in lsctables.py.
	See insert_from_url() for a description of the additional
	arguments.
	"""
	verbose = kwargs.get("verbose", False)

	#
	# load documents
	#

	for n, url in enumerate(urls, 1):
		if verbose:
			sys.stderr.write("%d/%d:" % (n, len(urls)))
		insert_from_url(url, contenthandler = contenthandler, **kwargs)

	#
	# done.  build indexes
	#

	dbtables.build_indexes(contenthandler.connection, verbose)


#
# How to extract
#


def extract(connection, filename, table_names = None, verbose = False, xsl_file = None):
	"""
	Convert the database at the given connection to a tabular LIGO
	Light-Weight XML document.  The XML document is written to the file
	named filename.  If table_names is not None, it should be a
	sequence of strings and only the tables in that sequence will be
	converted.  If verbose is True then progress messages will be
	printed to stderr.
	"""
	xmldoc = ligolw.Document()
	xmldoc.appendChild(dbtables.get_xml(connection, table_names))
	ligolw_utils.write_filename(xmldoc, filename, gz = (filename or "stdout").endswith(".gz"), verbose = verbose, xsl_file = xsl_file)

	# delete cursors
	xmldoc.unlink()
