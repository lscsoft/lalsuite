
# Copyright (C) 2012  LIGO Scientific Collaboration
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

"""
Some convenience logging classes courtesy of Leo Singer, provided as is.

Usage:
import logging
import ligo.gracedb.rest
import ligo.gracedb.logging
 
logging.basicConfig()
log = logging.getLogger('testing')
 
gracedb = ligo.gracedb.rest.GraceDb()
graceid = 'T62829'
 
log.addHandler(ligo.gracedb.logging.GraceDbLogHandler(gracedb, graceid))

# The following will create a log entry on the gracedb server
# (if the log level permits)
#
log.warn("this is a warning")
"""
 
# Perform explicit absolute import of Python standard library logging module.
# 'import logging' would not work here, because it would be interpreted as this
# module itself.
# module itself.g

logging = __import__('logging', level=0)

class GraceDbLogStream(object):
    def __init__(self, gracedb, graceid):
        self.gracedb = gracedb
        self.graceid = graceid
    def flush(self):
        pass
    def write(self, text):
        self.gracedb.writeLog(self.graceid, text)
 
class GraceDbLogFormatter(logging.Formatter):
    def __init__(self):
        logging.Formatter.__init__(self, logging.BASIC_FORMAT)
    def format(self, record):
        s = logging.Formatter.format(self, record)
        return '<div style="white-space:pre-wrap">' + s.strip("\n") + '</div>'
 
class GraceDbLogHandler(logging.StreamHandler):
    def __init__(self, gracedb, graceid):
        stream = GraceDbLogStream(gracedb, graceid)
        logging.StreamHandler.__init__(self, stream)
        self.formatter = GraceDbLogFormatter()
