#!/usr/bin/env @PYTHONPROG@
"""
This program makes the S5 high mass post processing dag

$Id$

This program creates cache files for the output of inspiral hipe
"""

__author__ = 'Chad Hanna <channa@phys.lsu.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

##############################################################################
# import standard modules and append the lalapps prefix to the python path
import sys, os, copy, math
import math
import socket, time
import re, string
from optparse import *
import tempfile
import ConfigParser
import urlparse
from UserDict import UserDict
sys.path.append('@PYTHONLIBDIR@')
import subprocess

##############################################################################
# import the modules we need to build the pipeline
from glue import pipeline
from glue import lal
from glue import segments
from glue import segmentsUtils
from pylal.webUtils import *
from pylal.webCondor import *
from lalapps import inspiral
from pylal import fu_utils
from glue.ligolw import lsctables

class hm_post_DAG(pipeline.CondorDAG):

  def __init__(self, config_file, log_path):
    self.basename = re.sub(r'\.ini',r'', config_file)
    tempfile.tempdir = log_path
    tempfile.template = self.basename + '.dag.log.'
    logfile = tempfile.mktemp()
    fh = open( logfile, "w" )
    fh.close()
    pipeline.CondorDAG.__init__(self,logfile)
    self.set_dag_file(self.basename)
    self.jobsDict = {}

class sqlite_job(pipeline.CondorDAGJob):
  """
  A sqlite3 job
  """
  def __init__(self, cp, tag_base='SQLITE3'):
    """
    """
    self.__prog__ = 'sqlite3'
    self.__executable = string.strip(cp.get('condor','sqlite3'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd("input","$(macroinput)")
    self.add_condor_cmd('getenv','True')
    self.tag_base = tag_base
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
    self.set_sub_file(tag_base+'.sub')
    self.set_stdout_file('logs/'+tag_base+'-$(macroid)-$(process).out')
    self.set_stderr_file('logs/'+tag_base+'-$(macroid-$(process)).err')

class ligolw_sqlite_job(pipeline.CondorDAGJob):
  """
  A ligolw_sqlite job
  """
  def __init__(self, cp, tag_base='LIGOLW_SQLITE'):
    """
    """
    self.__prog__ = 'ligolw_sqlite'
    self.__executable = string.strip(cp.get('condor','bash'))
    self.ligolw_sqlite = string.strip(cp.get('condor','ligolw_sqlite'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.tag_base = tag_base
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
    self.set_sub_file(tag_base+'.sub')
    self.set_stdout_file('logs/'+tag_base+'-$(macroid)-$(process).out')
    self.set_stderr_file('logs/'+tag_base+'-$(macroid)-$(process).err')


class ligolw_inspinjfind_job(pipeline.CondorDAGJob):
  """
  A ligolw_inspinjfind_job
  """
  def __init__(self, cp, tag_base='LIGOLW_INSPINJFIND'):
    """
    """
    self.__prog__ = 'ligolw_inspinjfind'
    self.__executable = string.strip(cp.get('condor','ligolw_inspinjfind'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.tag_base = tag_base
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
    self.set_sub_file(tag_base+'.sub')
    self.set_stdout_file('logs/'+tag_base+'-$(macroid)-$(process).out')
    self.set_stderr_file('logs/'+tag_base+'-$(macroid)-$(process).err')


class lalapps_newcorse_job(pipeline.CondorDAGJob):
  """
  A lalapps_newcorse_job
  """
  def __init__(self, cp, tag_base='LALAPPS_NEWCORSE'):
    """
    """
    self.__prog__ = 'lalapps_newcorse'
    self.__executable = string.strip(cp.get('condor','lalapps_newcorse'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.tag_base = tag_base
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
    self.set_sub_file(tag_base+'.sub')
    self.set_stdout_file('logs/'+tag_base+'-$(macroid)-$(process).out')
    self.set_stderr_file('logs/'+tag_base+'-$(macroid)-$(process).err')

class ligolw_segments_job(pipeline.CondorDAGJob):
  """
  A ligolw_segments_job
  """
  def __init__(self, cp, tag_base='LIGOLW_SEGMENTS'):
    """
    """
    self.__prog__ = 'ligolw_segments'
    self.__executable = string.strip(cp.get('condor','ligolw_segments'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.tag_base = tag_base
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
    self.set_sub_file(tag_base+'.sub')
    self.set_stdout_file('logs/'+tag_base+'-$(macroid)-$(process).out')
    self.set_stderr_file('logs/'+tag_base+'-$(macroid)-$(process).err')

class ligolw_thinca_to_coinc_job(pipeline.CondorDAGJob):
  """
  A ligolw_thinca_to_coinc_job
  """
  def __init__(self, cp, tag_base='LIGOLW_THINCA_TO_COINC'):
    """
    """
    self.__prog__ = 'ligolw_thinca_to_coinc'
    self.__executable = string.strip(cp.get('condor','ligolw_thinca_to_coinc'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.tag_base = tag_base
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
    self.set_sub_file(tag_base+'.sub')
    self.set_stdout_file('logs/'+tag_base+'-$(macroid)-$(process).out')
    self.set_stderr_file('logs/'+tag_base+'-$(macroid)-$(process).err')

class hm_upperlimit_job(pipeline.CondorDAGJob):
  """
  A hm_upperlimit_job
  """
  def __init__(self, cp, tag_base='HM_UPPERLIMIT'):
    """
    """
    self.__prog__ = 'hm_upperlimit'
    self.__executable = string.strip(cp.get('condor','hm_upperlimit'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.tag_base = tag_base
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
    self.set_sub_file(tag_base+'.sub')
    self.set_stdout_file('logs/'+tag_base+'-$(macroid)-$(process).out')
    self.set_stderr_file('logs/'+tag_base+'-$(macroid)-$(process).err')

class ligolw_sqlite_node(pipeline.CondorDAGNode):
  """
  """
  def __init__(self, job, dag, database, xml, id, p_node=[], replace=True, extract=False):

    pipeline.CondorDAGNode.__init__(self,job)
    cline = job.ligolw_sqlite + ' --database ' + database + ' --tmp-space /tmp --verbose '
    if replace: cline += " --replace "
    if extract: cline += " --extract " 
    cline += xml
    fn = "bash_scripts/ligolw_sqlite"+str(id)+".sh"
    f = open(fn,"w")
    f.write(cline)
    f.close
    self.add_macro("macroid", id)
    self.add_file_arg(fn)

    for p in p_node:
      self.add_parent(p)
    dag.add_node(self)

class ligolw_thinca_to_coinc_node(pipeline.CondorDAGNode):
  """
  """
  def __init__(self, job, dag, cache, vetoes, veto_name, prefix, id, effsnrfac=250.0, p_node=[], instruments='H1,H2,L1'):
    pipeline.CondorDAGNode.__init__(self,job)
    self.add_var_opt("ihope-cache", cache)
    self.add_var_opt("veto-segments", vetoes)
    self.add_var_opt("veto-segments-name",veto_name)
    self.add_var_opt("output-prefix",prefix)
    self.add_var_opt("effective-snr-factor",effsnrfac)
    self.add_var_opt("instruments",instruments)

    self.add_macro("macroid", id)
    for p in p_node:
      self.add_parent(p)
    dag.add_node(self)

class sqlite_node(pipeline.CondorDAGNode):
  """
  """
  def __init__(self, job, dag, database, sqlfile, id, p_node=[]):

    pipeline.CondorDAGNode.__init__(self,job)
    self.add_macro("macroid", id)
    self.add_file_arg(database)
    self.add_macro("macroinput", sqlfile)
    for p in p_node:
      self.add_parent(p)
    dag.add_node(self)

class ligolw_inspinjfind_node(pipeline.CondorDAGNode):
  """
  """
  def __init__(self, job, dag, xml, id, p_node=[]):

    pipeline.CondorDAGNode.__init__(self,job)
    self.add_var_arg(xml)
    self.add_macro("macroid", id)
    for p in p_node:
      self.add_parent(p)
    dag.add_node(self)

class ligolw_segments_node(pipeline.CondorDAGNode):
  """
  """
  def __init__(self, job, dag, ifodict, name, output, id, p_node=[], coalesce=True):
    pipeline.CondorDAGNode.__init__(self,job)
    # HA HA, we win!
    self.add_var_opt("insert-from-segwizard", " --insert-from-segwizard ".join(["%s=%s" % (instrument.upper(), filename) for instrument, filename in ifodict.items()]))
    self.add_var_opt("name",name)
    self.add_var_opt("output",output)
    self.add_macro("macroid", id)
    if coalesce: self.add_var_opt("coalesce","")
    for p in p_node:
      self.add_parent(p)
    dag.add_node(self)

class lalapps_newcorse_node(pipeline.CondorDAGNode):
  """
  """
  def __init__(self, job, dag, veto_segments, veto_segments_name, database, id, p_node=[],instruments = "H1,H2,L1", mass_bins="0,50,85,inf", live_time_program="thinca"):
    pipeline.CondorDAGNode.__init__(self,job)
    #FIXME make temp space?
    #self.add_var_opt("tmp-space","/tmp")
    self.add_var_opt("mass-bins",mass_bins)
    self.add_var_opt("live-time-program",live_time_program)
    self.add_var_opt("veto-segments",veto_segments)
    self.add_var_opt("veto-segments-name",veto_segments_name)
    self.add_var_arg(database)
    self.add_macro("macroid", id)
    for p in p_node:
      self.add_parent(p)
    dag.add_node(self)

class hm_upperlimit_node(pipeline.CondorDAGNode):
  """
  hm_upperlimit.py --instruments --output-name-tag --full-data-file --inj-data-glob --bootstrap-iterations --veto-segments --veto-segments-name
  """
  def __init__(self, job, dag, ifos, output_name_tag, full_data_file, inj_data_glob, bootstrap_iterations, veto_segments, veto_segments_name, id, p_node=[]):
    pipeline.CondorDAGNode.__init__(self,job)
    #FIXME make temp space?
    #self.add_var_opt("tmp-space","/tmp")
    self.add_var_opt("instruments",ifos)
    self.add_var_opt("output-name-tag",output_name_tag)
    self.add_var_opt("full-data-file",full_data_file)
    self.add_var_opt("inj-data-glob",inj_data_glob)
    self.add_var_opt("bootstrap-iterations",bootstrap_iterations)
    self.add_var_opt("veto-segments",veto_segments)
    self.add_var_opt("veto-segments-name",veto_segments_name)
    self.add_macro("macroid", id)
    for p in p_node:
      self.add_parent(p)
    dag.add_node(self)

def ifo_seg_dict(cp):
  out = {}
  out["H1"] = string.strip(cp.get('input','h1vetosegments'))
  out["H2"] = string.strip(cp.get('input','h2vetosegments'))
  out["L1"] = string.strip(cp.get('input','l1vetosegments'))
  return out

def grep(string, inname, outname):
    o = open(outname, "w")
    print "grepping " + inname + " for " + string + " and sending it to " + outname
    expr = re.compile(string)
    o.write(''.join(filter(expr.search,open(inname).readlines())))

###############################################################################
# MAIN PROGRAM
###############################################################################

cp = ConfigParser.ConfigParser()
cp.read("hm_post.ini")

try: os.mkdir("logs")
except: pass

try: os.mkdir("bash_scripts")
except: pass

cats = ["CAT_2","CAT_3"]
types = ["FULL_DATA"]
FULLDATACACHE = string.strip(cp.get('input','fulldatacache'))
INJCACHE = string.strip(cp.get('input','injcache'))
dag = hm_post_DAG("hm_post.ini", string.strip(cp.get('output','logpath')))
# to get injection file entries from the cache

#break down the cache to save on parsing
grep('HL', INJCACHE, "inj.cache")

#Setup jobs
sqliteJob = sqlite_job(cp)
ligolwSqliteJob = ligolw_sqlite_job(cp)
ligolwInspinjfindJob = ligolw_inspinjfind_job(cp)
lalappsNewcorseJob = lalapps_newcorse_job(cp)
ligolwSegmentsJob = ligolw_segments_job(cp)
ligolwThincaToCoincJob =  ligolw_thinca_to_coinc_job(cp)
hmUpperlimitJob = hm_upperlimit_job(cp)

n = 0
#Do the segments node
segNode = {}
for cat in cats:
  segNode[cat] = ligolw_segments_node(ligolwSegmentsJob, dag, ifo_seg_dict(cp), "vetoes", "vetoes_"+cat+".xml.gz", n); n+=1

#Some initialization
ligolwThincaToCoincNode = {}
sqliteNodeSimplify = {}
sqliteNodeRemoveH1H2 = {}
sqliteNodeCluster = {}
ligolwSqliteNodeVetoes = {}
ligolwSqliteNode = {}
ligolwSqliteNode2 = {}
ligolwSqliteNode3 = {}
ligolwSqliteNode4 = {}
ligolwInspinjfindNode = {}
lallappsNewcorseNode = {}
hmUpperlimitNode = {}
db = {}

# to get injection file entries from the cache
injcache = map(lal.CacheEntry, file("inj.cache"))
inj = injcache[0]
start_time = inj.segment[0]
end_time = inj.segment[1]

for type in types:
  for cat in cats:
    #break down the cache to save on parsing
    grep(type + ".*" + cat, FULLDATACACHE, type + cat + ".cache")
    try: os.mkdir(type+cat)
    except: pass
    ligolwThincaToCoincNode[type+cat] = ligolw_thinca_to_coinc_node(ligolwThincaToCoincJob, dag, type+cat+".cache", "vetoes_"+cat+".xml.gz", "vetoes", type+cat+"/S5_HM", n, effsnrfac=50, p_node=[segNode[cat]]); n+=1
    database = type+cat+".sqlite"
    try: db[cat].append(database) 
    except: db[cat] = [database]
    ligolwSqliteNode[type+cat] = ligolw_sqlite_node(ligolwSqliteJob, dag, database, type+cat+"/S5_HM_*"+type+"*"+cat+"*.xml.gz", n, p_node=[ligolwThincaToCoincNode[type+cat]], replace=True); n+=1
    ligolwSqliteNodeVetoes[type+cat] = ligolw_sqlite_node(ligolwSqliteJob, dag, database, "vetoes_"+cat+".xml.gz", n, p_node=[ligolwSqliteNode[type+cat]], replace=False); n+=1
    sqliteNodeSimplify[type+cat] = sqlite_node(sqliteJob, dag, database, string.strip(cp.get('input',"simplify")), n, p_node=[ligolwSqliteNodeVetoes[type+cat]]); n+=1
    sqliteNodeRemoveH1H2[type+cat] = sqlite_node(sqliteJob, dag, database, string.strip(cp.get('input',"remove_h1h2")),n, p_node=[sqliteNodeSimplify[type+cat]]); n+=1
    sqliteNodeCluster[type+cat] = sqlite_node(sqliteJob, dag, database, string.strip(cp.get('input',"cluster")),n, p_node=[sqliteNodeRemoveH1H2[type+cat]]); n+=1


for inj in injcache:
  for cat in cats:
    type = "_".join(inj.description.split("_")[2:])
    url = inj.url
    cachefile = type + cat + ".cache"
    try: os.mkdir(type+cat)
    except: pass
    #break down the cache
    grep(type + '.*' + cat, INJCACHE, cachefile)
    ligolwThincaToCoincNode[type+cat] = ligolw_thinca_to_coinc_node(ligolwThincaToCoincJob, dag, cachefile, "vetoes_"+cat+".xml.gz", "vetoes", type+cat+"/S5_HM_INJ", n, effsnrfac=50, p_node=[segNode[cat]]);n+=1
    database = type+cat+".sqlite"
    try: db[cat].append(database)
    except: db[cat] = [database]
    ligolwSqliteNode[type+cat] = ligolw_sqlite_node(ligolwSqliteJob, dag, database, type+cat+"/S5_HM_INJ*"+type+"*"+cat+"*.xml.gz",n, p_node=[ligolwThincaToCoincNode[type+cat]], replace=True);n+=1
    ligolwSqliteNode2[type+cat] = ligolw_sqlite_node(ligolwSqliteJob, dag, database, url, n, p_node=[ligolwSqliteNode[type+cat]], replace=False);n+=1
    ligolwSqliteNodeVetoes[cat] = ligolw_sqlite_node(ligolwSqliteJob, dag, database, "vetoes_"+cat+".xml.gz", n, p_node=[ligolwSqliteNode2[type+cat]], replace=False); n+=1
    sqliteNodeSimplify[type+cat] = sqlite_node(sqliteJob, dag, database, string.strip(cp.get('input',"simplify")), n, p_node=[ligolwSqliteNode2[type+cat]]);n+=1
    sqliteNodeRemoveH1H2[type+cat] = sqlite_node(sqliteJob, dag, database, string.strip(cp.get('input',"remove_h1h2")),n, p_node=[sqliteNodeSimplify[type+cat]]);n+=1
    sqliteNodeCluster[type+cat] = sqlite_node(sqliteJob, dag, database, string.strip(cp.get('input',"cluster")),n, p_node=[sqliteNodeRemoveH1H2[type+cat]]);n+=1

    ligolwSqliteNode3[type+cat] = ligolw_sqlite_node(ligolwSqliteJob, dag, database, database+".xml.gz", n, p_node=[sqliteNodeCluster[type+cat]], replace=False, extract=True); n+=1
    ligolwInspinjfindNode[type+cat] = ligolw_inspinjfind_node(ligolwInspinjfindJob, dag, database+".xml.gz", n, p_node=[ligolwSqliteNode3[type+cat]]);n+=1
    ligolwSqliteNode4[type+cat] = ligolw_sqlite_node(ligolwSqliteJob, dag, database, database+".xml.gz", n, p_node=[ligolwInspinjfindNode[type+cat]], replace=True);n+=1

# New corse
# First work out the parent/child relationships
p_nodes = []
for k in ligolwInspinjfindNode.keys():
  p_nodes.append(ligolwInspinjfindNode[k])
for k in sqliteNodeCluster.keys():
  p_nodes.append(sqliteNodeCluster[k])

for cat in cats:
  lallappsNewcorseNode[cat] = lalapps_newcorse_node(lalappsNewcorseJob, dag, "vetoes_%s.xml.gz" % cat, "vetoes", " ".join(db[cat]), n, p_nodes);n+=1

for cat in ['CAT_3']: 
  hmUpperlimitNode[cat] = hm_upperlimit_node(hmUpperlimitJob, dag, "H1,H2,L1", "", "FULL_DATACAT_3.sqlite", "*INJCAT_3.sqlite", 10000, "vetoes_%s.xml.gz" % cat, "vetoes", n, p_node=[lallappsNewcorseNode[cat]]);n+=1
  hmUpperlimitNode[cat] = hm_upperlimit_node(hmUpperlimitJob, dag, "H1,L1", "", "FULL_DATACAT_3.sqlite", "*INJCAT_3.sqlite", 10000, "vetoes_%s.xml.gz" % cat, "vetoes", n, p_node=[lallappsNewcorseNode[cat]]);n+=1
  hmUpperlimitNode[cat] = hm_upperlimit_node(hmUpperlimitJob, dag, "H2,L1", "", "FULL_DATACAT_3.sqlite", "*INJCAT_3.sqlite", 10000, "vetoes_%s.xml.gz" % cat, "vetoes", n, p_node=[lallappsNewcorseNode[cat]]);n+=1

dag.write_sub_files()
dag.write_dag()
dag.write_script()

print "\nYour database output should be...\n"
for cat in cats:
  print cat + ":\n", " ".join(db[cat]) + "\n"

