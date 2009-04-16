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
    self.add_condor_cmd('getenv','True')
    self.tag_base = tag_base
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
    self.set_sub_file(tag_base+'.sub')
    self.set_stdout_file('logs/'+tag_base+'-$(macroid).out')
    self.set_stderr_file('logs/'+tag_base+'-$(macroid).err')

class ligolw_sqlite_job(pipeline.CondorDAGJob):
  """
  A ligolw_sqlite job
  """
  def __init__(self, cp, tag_base='LIGOLW_SQLITE'):
    """
    """
    self.__prog__ = 'ligolw_sqlite'
    self.__executable = string.strip(cp.get('condor','ligolw_sqlite'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.tag_base = tag_base
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
    self.set_sub_file(tag_base+'.sub')
    self.set_stdout_file('logs/'+tag_base+'-$(macroid).out')
    self.set_stderr_file('logs/'+tag_base+'-$(macroid).err')


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
    self.set_stdout_file('logs/'+tag_base+'-$(macroid).out')
    self.set_stderr_file('logs/'+tag_base+'-$(macroid).err')


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
    self.set_stdout_file('logs/'+tag_base+'-$(macroid).out')
    self.set_stderr_file('logs/'+tag_base+'-$(macroid).err')

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
    self.set_stdout_file('logs/'+tag_base+'-$(macroid).out')
    self.set_stderr_file('logs/'+tag_base+'-$(macroid).err')

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
    self.set_stdout_file('logs/'+tag_base+'-$(macroid).out')
    self.set_stderr_file('logs/'+tag_base+'-$(macroid).err')

class ligolw_sqlite_node(pipeline.CondorDAGNode):
  """
  """
  def __init__(self, job, dag, database, xml, p_node=[], replace=False, extract=False):

    pipeline.CondorDAGNode.__init__(self,job)
    self.add_var_opt("database", database)
    self.add_var_opt("tmp-space", '\\tmp')
    self.add_var_opt("verbose","")
    self.add_file_arg(xml)
    if replace: self.add_var_opt("replace","")
    if extract: self.add_var_opt("extract","")

    for p in p_node:
      self.add_parent(p)
    dag.add_node(self)

class ligolw_thinca_to_coinc_node(pipeline.CondorDAGNode):
  """
  """
  def __init__(self, job, dag, cache, vetoes, veto_name, prefix, effsnrfac=250.0, p_node=[], replace=False):

    pipeline.CondorDAGNode.__init__(self,job)
    self.add_var_opt("ihope-cache", cache)
    self.add_var_opt("veto-segments", vetoes)
    self.add_var_opt("veto-segments-name",veto_name)
    self.add_var_opt("output-prefix",prefix)
    self.add_var_opt("effective-snr-factor",effsnrfac)
    for p in p_node:
      self.add_parent(p)
    dag.add_node(self)

class sqlite_node(pipeline.CondorDAGNode):
  """
  """
  def __init__(self, job, dag, database, sqlfile, p_node=[]):

    pipeline.CondorDAGNode.__init__(self,job)
    self.add_var_arg(database+" < " + sqlfile)
    for p in p_node:
      self.add_parent(p)
    dag.add_node(self)

class ligolw_inspinjfind_node(pipeline.CondorDAGNode):
  """
  """
  def __init__(self, job, dag, xml, p_node=[]):

    pipeline.CondorDAGNode.__init__(self,job)
    self.add_var_arg(xml)
    for p in p_node:
      self.add_parent(p)
    dag.add_node(self)

class ligolw_segments_node(pipeline.CondorDAGNode):
  """
  """
  def __init__(self, job, dag, ifodict, name, output, p_node=[], coalesce=True):
    pipeline.CondorDAGNode.__init__(self,job)
    for k in ifodict.keys():
      print ifodict[k]
      if ifodict[k]: self.add_var_opt("insert-from-segwizard="+k.upper()+"=",ifodict[k])
    self.add_var_opt("name",name)
    self.add_var_opt("output",output)
    if coalesce: self.add_var_opt("coalesce","")
    for p in p_node:
      self.add_parent(p)
    dag.add_node(self)


def ifo_seg_dict(cp):
  out = {}
  out["H1"] = string.strip(cp.get('input','h1vetosegments'))
  out["H2"] = string.strip(cp.get('input','h2vetosegments'))
  out["L1"] = string.strip(cp.get('input','l1vetosegments'))
  return out

###############################################################################
# MAIN PROGRAM
###############################################################################

cp = ConfigParser.ConfigParser()
cp.read("hm_post.ini")

try: os.mkdir("logs")
except: pass

cats = ["CAT_2","CAT_3"]
types = ["FULL_DATA"]
FULLDATACACHE = string.strip(cp.get('input','fulldatacache'))
INJCACHE = string.strip(cp.get('input','injcache'))
dag = hm_post_DAG("hm_post.ini", string.strip(cp.get('output','logpath')))
# to get injection file entries from the cache
command = 'grep HL ' + INJCACHE + " > " + "inj.cache"
print command
popen = os.popen(command)

#Setup jobs
sqliteJob = sqlite_job(cp)
ligolwSqliteJob = ligolw_sqlite_job(cp)
ligolwInspinjfindJob = ligolw_inspinjfind_job(cp)
lalappsNewcorseJob = lalapps_newcorse_job(cp)
ligolwSegmentsJob = ligolw_segments_job(cp)
ligolwThincaToCoincJob =  ligolw_thinca_to_coinc_job(cp)

#Do the segments node
segNode = {}
for cat in cats:
  segNode[cat] = ligolw_segments_node(ligolwSegmentsJob, dag, ifo_seg_dict(cp), "vetoes", "vetoes_"+cat+".xml.gz") 

#Run thinca_to_coinc on zero lag and time slides
for type in types:
  for cat in cats:
    command = 'grep "'  + type + ".*" + cat + '" ' + FULLDATACACHE + " > " + type + cat + ".cache"
    print command
    popen = os.popen(command)
    ligolwThincaToCoincNode = ligolw_thinca_to_coinc_node(ligolwThincaToCoincJob, dag, type+cat+".cache", "vetoes_"+cat+".xml.gz", "vetoes", "S5_HM", effsnrfac=50,p_node=[segNode[cat]])
    database = type+cat+".sqlite"
    ligolwSqliteNode = ligolw_sqlite_node(ligolwSqliteJob, dag, database, "S5_HM_*"+type+"*"+cat+"*.xml.gz", p_node=[ligolwThincaToCoincNode], replace=True)
    sqliteNodeSimplify = sqlite_node(sqliteJob, dag, database, string.strip(cp.get('input',"simplify")), p_node=[ligolwSqliteNode])
    sqliteNodeRemoveH1H2 = sqlite_node(sqliteJob, dag, database, string.strip(cp.get('input',"remove_h1h2")),p_node=[sqliteNodeSimplify])
    sqliteNodeCluster = sqlite_node(sqliteJob, dag, database, string.strip(cp.get('input',"cluster")),p_node=[sqliteNodeRemoveH1H2])


# to get injection file entries from the cache
injcache = map(lal.CacheEntry, file("inj.cache"))

for inj in injcache:
  for cat in cats:
    type = "_".join(inj.description.split("_")[2:])
    url = inj.url
    cachefile = type + cat + ".cache"
    command = 'grep "' + type + '.*' + cat + '" ' + INJCACHE +" > " + cachefile
    print command
    popen = os.popen(command)
    ligolwThincaToCoincNode = ligolw_thinca_to_coinc_node(ligolwThincaToCoincJob, dag, cachefile, "vetoes_"+cat+".xml.gz", "vetoes", "S5_HM_INJ", effsnrfac=50)
    database = type+cat+".sqlite"
    ligolwSqliteNode = ligolw_sqlite_node(ligolwSqliteJob, dag, database, "S5_HM_INJ*"+type+"*"+cat+"*.xml.gz", p_node=[ligolwThincaToCoincNode], replace=True)
    ligolwSqliteNode2 = ligolw_sqlite_node(ligolwSqliteJob, dag, database, url, p_node=[ligolwSqliteNode], replace=True)
    sqliteNodeSimplify = sqlite_node(sqliteJob, dag, database, string.strip(cp.get('input',"simplify")), p_node=[ligolwSqliteNode2])
    sqliteNodeRemoveH1H2 = sqlite_node(sqliteJob, dag, database, string.strip(cp.get('input',"remove_h1h2")),p_node=[sqliteNodeSimplify])
    sqliteNodeCluster = sqlite_node(sqliteJob, dag, database, string.strip(cp.get('input',"cluster")),p_node=[sqliteNodeRemoveH1H2])
    ligolwSqliteNode3 = ligolw_sqlite_node(ligolwSqliteJob, dag, database, database+".xml.gz", p_node=[sqliteNodeCluster], replace=False, extract=True)
    ligolwInspinjfindNode = ligolw_inspinjfind_node(ligolwInspinjfindJob, dag, database+".xml.gz", p_node=[ligolwSqliteNode3])
    ligolwSqliteNode4 = ligolw_sqlite_node(ligolwSqliteJob, dag, database, database+".xml.gz", p_node=[ligolwInspinjfindNode], replace=True)

dag.write_sub_files()
dag.write_dag()

 
