try:
        import sqlite3
except ImportError:
        # pre 2.5.x
        from pysqlite2 import dbapi2 as sqlite3
from glue.ligolw import dbtables

filename = "EOBNR_INJ_ONE.sqlite"
local_disk = "/tmp"        # use None to disable
working_filename = dbtables.get_connection_filename(filename, tmp_path = local_disk, verbose = True)

connection = sqlite3.connect(working_filename)

dbtables.DBTable_set_connection(connection)
xmldoc = dbtables.get_xml(connection)


xmldoc.write()

