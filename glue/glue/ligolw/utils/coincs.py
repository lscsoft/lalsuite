# Copyright (C) 2006-2010,2012-2013  Kipp Cannon
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
A collection of utilities to assist with the use of coinc tables.
"""


from glue import git_version
from .. import lsctables


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                                    Tables
#
# =============================================================================
#


def get_coinc_def_id(xmldoc, search, coinc_type, create_new = True, description = u""):
	"""
	Wrapper for the get_coinc_def_id() method of the CoincDefiner table
	class in glue.ligolw.lsctables.  This wrapper will optionally
	create a new coinc_definer table in the document if one does not
	already exist.
	"""
	try:
		coincdeftable = lsctables.CoincDefTable.get_table(xmldoc)
	except ValueError:
		# table not found
		if not create_new:
			raise
		# FIXME:  doesn't work if the document is stored in a
		# database.
		coincdeftable = lsctables.New(lsctables.CoincDefTable)
		xmldoc.childNodes[0].appendChild(coincdeftable)
	# make sure the next_id attribute is correct
	coincdeftable.sync_next_id()
	# get the id
	return coincdeftable.get_coinc_def_id(search, coinc_type, create_new = create_new, description = description)
