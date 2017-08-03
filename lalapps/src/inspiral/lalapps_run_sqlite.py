#
# Copyright (C) 2010 Chad Hanna
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


from optparse import OptionParser
import sqlite3
import sys

from glue.ligolw import ligolw
from glue.ligolw import dbtables
from glue.ligolw.utils import ligolw_sqlite
from lalapps import git_version

__author__ = "Chad Hanna <channa@ligo.caltech.edu>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


def parse_command_line():
	parser = OptionParser(
		version = "Name: %%prog\n%s" % git_version.verbose_msg, usage = "%prog (--sql CODE|--sql-file FILENAME) [options] database1.sqlite database2.xml ..."
	)
	parser.add_option("-t", "--tmp-space", metavar = "path", help = "Path to a directory suitable for use as a work area while manipulating the database file.  The database file will be worked on in this directory, and then moved to the final location when complete.  This option is intended to improve performance when running in a networked environment, where there might be a local disk with higher bandwidth than is available to the filesystem on which the final output will reside.")
	parser.add_option("-c", "--sql", metavar = "code", help = "Execute this SQL code.")
	parser.add_option("-s", "--sql-file", metavar = "filename", help = "Execute this SQL file.")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, filenames = parser.parse_args()
	if bool(options.sql) + bool(options.sql_file) != 1:
		raise ValueError("must set exactly one of --sql or --sql-file")
	return options, (filenames or [])


options, databases = parse_command_line()


if options.sql_file:
	# Parse the sql file into something that can be executed in sequence
	sql = [line.strip() for line in open(options.sql_file)]
	# Remove comments and pragmas
	sql = [s for s in sql if not s.startswith("--") and not s.startswith("PRAGMA")]
	sql = "\n".join(sql)
	sql = [statement.strip() for statement in sql.split(";\n")]
elif options.sql:
	sql = [statement.strip() for statement in options.sql.split(";")]
else:
	raise NotImplemented
# remove no-ops
sql = [statement for statement in sql if statement]

# Apply the SQL to all the databases sequentially
for n, filename in enumerate(databases, 1):
	#
	# access the document
	#

	if options.verbose:
		print >>sys.stderr, "%d/%d:" % (n, len(databases)),
	if filename.endswith(".xml") or filename.endswith(".xml.gz"):
		# load XML file into in-ram database for processing
		fileformat = "xml"
		working_filename = None
		connection = sqlite3.connect(":memory:")
		@dbtables.use_in
		class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
			connection = connection
		ligolw_sqlite.insert_from_url(filename, contenthandler = LIGOLWContentHandler, preserve_ids = True, verbose = options.verbose)
		del LIGOLWContentHandler
	else:
		# assume it's an SQLite file
		fileformat = "sqlite"
		working_filename = dbtables.get_connection_filename(filename, tmp_path = options.tmp_space, verbose = options.verbose)
		connection = sqlite3.connect(working_filename)
	# also use the scratch space for sqlite's temp store, but don't try
	# to do so if the working filename is the same as the original
	# filename as that most likely indicates some problem was found
	# with the scratch space like that it's full
	if options.tmp_space is not None and working_filename != filename:
		dbtables.set_temp_store_directory(connection, options.tmp_space, verbose = options.verbose)

	#
	# apply SQL
	#

	if options.verbose:
		print >>sys.stderr, "Executing SQL ..."
	cursor = connection.cursor()
	for statement in sql:
		if options.verbose:
			print >>sys.stderr, statement
		cursor.execute(statement)
		connection.commit()
	cursor.close()
	if options.verbose:
		print >>sys.stderr, "... Done."

	#
	# commit changes
	#

	if fileformat == "xml":
		# overwrite XML file with in-ram database' contents
		ligolw_sqlite.extract(connection, filename, verbose = options.verbose)
		connection.close()
	elif fileformat == "sqlite":
		# return SQLite file to its home
		connection.close()
		dbtables.put_connection_filename(filename, working_filename, verbose = options.verbose)
	else:
		raise Exception("internal error")
