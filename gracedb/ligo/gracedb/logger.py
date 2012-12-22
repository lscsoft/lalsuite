
"""
Some convenience logging classes courtesy of Leo Singer, provided as is.

Usage:
import logging
import ligo.gracedb.rest
import ligo.gracedb.logger
 
logging.basicConfig()
log = logging.getLogger('testing')
 
gracedb = ligo.gracedb.rest.GraceDb()
graceid = 'T62829'
 
log.addHandler(ligo.gracedb.logger.GraceDbLogHandler(gracedb, graceid))

# The following will create a log entry on the gracedb server
# (if the log level permits)
#
log.warn("this is a warning")
"""
 
import logging

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
