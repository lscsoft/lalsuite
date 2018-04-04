#
# Copyright (C) 2009  Larne Pekowsky
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
#				   Preamble
#
# =============================================================================
#


"""
Provides a unified interface for querying sqlite databases and a DB2 database
behind a ldbd interface
"""

import os
try:
    import pyRXP
except ImportError:
    import pyRXPU as pyRXP
from glue import ldbd

from glue import git_version
__date__ = git_version.date
__version__ = git_version.id
__author__  = "Larne Pekowsky <lppekows@physics.syr.edu>"


class QueryEngine:
	"""Abstract class.  Provides query() method that returns something that
	behaves like an array of tuples of results and a close() method that
	does any needed cleanup."""

	def query(sql):
		return None

	def close():
		pass


class SqliteQueryEngine(QueryEngine):
	"""QueryEngine for sqlite databases.  Really just a minimal wrapper
	around cursor"""

	connection = None
	cursor     = None

	def __init__(self, connection):
		self.connection = connection
	
	def query(self, sql):
		self.cursor = self.connection.cursor().execute(sql)
		ret         = []

		for row in self.cursor:
			ret.append(row)

		return ret

	def close(self):
		self.cursor.close()
		del self.cursor


class LdbdQueryEngine(QueryEngine):
	"""QueryEngine for databses behind ldbd.  Parses ligolw that a query
	returns into rows"""

	xmlparser = None
	lwtparser = None
	ligomd    = None

	rows      = None
	


	def __init__(self, client):
		def dtd_uri_callback(uri):
			if uri == 'http://ldas-sw.ligo.caltech.edu/doc/ligolwAPI/html/ligolw_dtd.txt':
				return 'file://localhost' + os.path.join( os.environ["GLUE_PREFIX"], 'etc/ligolw_dtd.txt' )
			else:
				return uri

		self.client	 = client
		self.xmlparser      = pyRXP.Parser()
		self.xmlparser.eoCB = dtd_uri_callback		
		self.lwtparser      = ldbd.LIGOLwParser()
		self.ligomd	 = ldbd.LIGOMetadata(self.xmlparser, self.lwtparser, None)


	def query(self, sql):
		xml = self.client.query(sql)
		
		# This is a kludge around bug 2317
                try:
                   self.client.__disconnect__()
		   self.client.__connect__(self.client.host, self.client.port, self.client.identity)
                except:
                   pass
 
		self.ligomd.parse(xml)
		res = self.ligomd.table
		self.rows = self.ligomd.table[res.keys()[0]]['stream']

		return self.rows
	
	def close(self):
		del self.rows

