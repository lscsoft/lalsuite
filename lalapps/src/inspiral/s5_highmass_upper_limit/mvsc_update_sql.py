try:
  import sqlite3
except ImportError:
  # pre 2.5.x
  from pysqlite2 import dbapi2 as sqlite3
from glue.ligolw import dbtables 
from pylal import SnglInspiralUtils
from pylal import db_thinca_rings
from time import clock,time
from optparse import *
import glob

usage="""
put the MVSC rank back into the SQL tables, after it's been turned into a likelihood, we put it in the likelihood column of the coinc_inspiral table
"""

time1=time()

parser=OptionParser(usage=usage, version="kari on 08/24/09")
parser.add_option("", "--files", default="*.dat", help="glob of .dat files that are output by SprOutputWriterApp, these should contain the MVSC rank in the last column")
parser.add_option("", "--infofiles", default="*_info.pat", help="glob of .pat files that are output by mvsc_get_doubles.py, these should contain the coinc_id (first column) for and sqlite database (second column) containing the trigger")
parser.add_option("", "--databases", default="*.sqlite", help="glob of SQLite databases that contain the tables that you're putting the MVSC likelhihood into")
(opts,args)=parser.parse_args()

files = glob.glob(opts.files)
infofiles = glob.glob(opts.infofiles)
databases = glob.glob(opts.databases)

def files_to_dict(f1, f1_info, databases, ldict):
  for f in databases: 
    try: ldict[f]
    except: ldict[f] = {}
  f1lines = open(f1).readlines()
  #the first line of the first file is a header
  f1lines.pop(0)
  f1_infolines = open(f1_info).readlines()
  if len(f1lines) != len(f1_infolines): 
    print>> sys.stderr, "file " + f1 + " and " +f1_info + " have different lengths"
    sys.exit(1)
  for i, f1row in enumerate(f1lines):
    likelihood = float(f1row.split()[-1])
    if likelihood == 1: likelihood = float("inf")
    else: likelihood = likelihood / (1.0 - likelihood)
    cid, dfile = f1_infolines[i].split()[0:2]
    ldict[dfile][cid] = likelihood
  return ldict
  
# first, initialize the databases, so that all of the likelihood values are set to 1
for database in databases:
  filename = database
  local_disk = None #"/tmp"
  working_filename = dbtables.get_connection_filename(filename, tmp_path = local_disk, verbose = True)
  connection = sqlite3.connect(working_filename)
  xmldoc = dbtables.get_xml(connection)
  cursor = connection.cursor()
  connection.cursor().execute("""
  UPDATE coinc_event
    SET likelihood = 1
  """)
  connection.commit()
  dbtables.put_connection_filename(filename, working_filename, verbose = True)

# pair up each file with the appropriate info file
files.sort()
#print files
infofiles.sort()
#print infofiles
if len(files) != len(infofiles): print 'error'
file_map=[]
for i in range(len(files)):
  file_map.append([files[i],infofiles[i]])
print file_map

ldict = {}
for pair in file_map:
  files_to_dict(str(pair[0]),str(pair[1]),databases, ldict)

for f in databases:
  filename = f
  local_disk = None #"/tmp"
  working_filename = dbtables.get_connection_filename(filename, tmp_path = local_disk, verbose = True)
  connection = sqlite3.connect(working_filename)
  xmldoc = dbtables.get_xml(connection)
  cursor = connection.cursor()
  def likelihood_from_coinc_id(id, ldict_f = ldict[f]):
    try: 
      return ldict_f[id]
    except: return 1
  if ldict[f]:  
    connection.create_function("likelihood_from_coinc_id", 1, likelihood_from_coinc_id)    
    cursor.execute("""
    UPDATE coinc_event
      SET likelihood = likelihood * likelihood_from_coinc_id(coinc_event.coinc_event_id)
    """)
    connection.commit() 
  dbtables.put_connection_filename(filename, working_filename, verbose = True)

time2=time()
elapsed_time=time2-time1
print elapsed_time


