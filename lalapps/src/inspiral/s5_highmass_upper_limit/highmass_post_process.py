#!/usr/bin/env @PYTHONPROG@
"""
This program makes the S5 high mass post processing dag
"""

__author__ = 'Chad Hanna <channa@caltech.edu>'
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
    self.__universe = string.strip(cp.get('condor','sqlite3_universe'))
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd("input","$(macroinput)")
    self.add_condor_cmd('getenv','True')
    self.tag_base = tag_base
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
    self.set_sub_file(tag_base+'.sub')
    self.set_stdout_file('logs/'+tag_base+'-$(macroid)-$(process).out')
    self.set_stderr_file('logs/'+tag_base+'-$(macroid)-$(process).err')

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


class lalapps_newcorse_combined_job(pipeline.CondorDAGJob):
  """
  A lalapps_newcorse_job
  """
  def __init__(self, cp, tag_base='LALAPPS_NEWCORSE_COMBINED'):
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

class far_plot_job(pipeline.CondorDAGJob):
  """
  A far_plot Job
  """
  def __init__(self, cp, tag_base='FAR_PLOT'):
    """
    """
    self.__prog__ = 'far_plot'
    self.__executable = string.strip(cp.get('condor','far_plot'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.tag_base = tag_base
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
    self.set_sub_file(tag_base+'.sub')
    self.set_stdout_file('logs/'+tag_base+'-$(macroid)-$(process).out')
    self.set_stderr_file('logs/'+tag_base+'-$(macroid)-$(process).err')

class ul_plot_job(pipeline.CondorDAGJob):
  """
  A ul_plot Job
  """
  def __init__(self, cp, tag_base='UL_PLOT'):
    """
    """
    self.__prog__ = 'ul_plot'
    self.__executable = string.strip(cp.get('condor','ul_plot'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.tag_base = tag_base
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
    self.set_sub_file(tag_base+'.sub')
    self.set_stdout_file('logs/'+tag_base+'-$(macroid)-$(process).out')
    self.set_stderr_file('logs/'+tag_base+'-$(macroid)-$(process).err')

class summary_page_job(pipeline.CondorDAGJob):
  """
  A summary page job
  """
  def __init__(self, cp, tag_base='SUMMARY_PAGE'):
    """
    """
    self.__prog__ = 'summary_page'
    self.__executable = string.strip(cp.get('condor','summary_page'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.tag_base = tag_base
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
    self.set_sub_file(tag_base+'.sub')
    self.add_arg(string.strip(cp.get('output','web_page')))
    self.set_stdout_file('logs/'+tag_base+'-$(macroid)-$(process).out')
    self.set_stderr_file('logs/'+tag_base+'-$(macroid)-$(process).err')

class ligolw_sqlite_node(pipeline.CondorDAGNode):
  """
  """
  def __init__(self, job, dag, database, xml_list, id, p_node=[], replace=True, extract=False,cache_pat=None):

    pipeline.CondorDAGNode.__init__(self,job)
    #FIXME add tmp file space
    cline = job.ligolw_sqlite + ' --database ' + database + ' --verbose '
    if replace: cline += " --replace "
    if extract: cline += " --extract " 
    if cache_pat: cline += " --input-cache ligolw_sqlite_" + cache_pat + ".cache "
    for xml in xml_list: cline += xml + " "
    fn = "bash_scripts/ligolw_sqlite"+str(id)+".sh"
    f = open(fn,"w")
    if cache_pat: 
      f.write("find $PWD -name '*" + cache_pat + "*.xml.gz' -print | sed -e 's?^?- - - - file://localhost?' > ligolw_sqlite_" + cache_pat + ".cache\n")
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
  def __init__(self, job, dag, cache, vetoes, veto_name, prefix, id, start_time, end_time, effsnrfac=50.0, p_node=[], instruments='H1,H2,L1'):
    pipeline.CondorDAGNode.__init__(self,job)
    self.add_var_opt("ihope-cache", cache)
    self.add_var_opt("veto-segments", vetoes)
    self.add_var_opt("veto-segments-name",veto_name)
    self.add_var_opt("output-prefix",prefix)
    self.add_var_opt("effective-snr-factor",effsnrfac)
    self.add_var_opt("instruments",instruments)
    self.add_var_opt("experiment-start-time", start_time)
    self.add_var_opt("experiment-end-time", end_time)

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
  def __init__(self, job, dag, veto_segments_name, database, id, p_node=[],instruments = "H1,H2,L1", mass_bins="0,50,85,inf", live_time_program="thinca", categories="mtotal-ifos-oninstruments", rank="snr"):
    pipeline.CondorDAGNode.__init__(self,job)
    #FIXME make temp space?
    #self.add_var_opt("tmp-space","/tmp")
    self.add_var_opt("categories", categories)
    if mass_bins: self.add_var_opt("mass-bins", mass_bins)
    self.add_var_opt("live-time-program",live_time_program)
    self.add_var_opt("veto-segments-name",veto_segments_name)
    self.add_var_opt("rank-by", rank)
    self.add_var_arg(database)
    self.add_macro("macroid", id)
    for p in p_node:
      self.add_parent(p)
    dag.add_node(self)

class hm_upperlimit_node(pipeline.CondorDAGNode):
  """
  hm_upperlimit.py --instruments --output-name-tag --full-data-file --inj-data-glob --bootstrap-iterations --veto-segments-name
  """
  def __init__(self, job, dag, ifos, output_name_tag, full_data_file, inj_data_glob, bootstrap_iterations, veto_segments_name, id, p_node=[]):
    pipeline.CondorDAGNode.__init__(self,job)
    #FIXME make temp space?
    #self.add_var_opt("tmp-space","/tmp")
    self.add_var_opt("instruments",ifos)
    self.add_var_opt("output-name-tag",output_name_tag)
    self.add_var_opt("full-data-file",full_data_file)
    self.add_var_opt("inj-data-glob",inj_data_glob)
    self.add_var_opt("bootstrap-iterations",bootstrap_iterations)
    self.add_var_opt("veto-segments-name",veto_segments_name)
    self.add_macro("macroid", id)
    for p in p_node:
      self.add_parent(p)
    dag.add_node(self)

class far_plot_node(pipeline.CondorDAGNode):
  """
  """
  def __init__(self, job, dag, database, id, p_node):
    pipeline.CondorDAGNode.__init__(self,job)
    self.add_macro("macroid", id)
    self.add_var_arg(database)
    for p in p_node:
      self.add_parent(p)
    dag.add_node(self)
    
class ul_plot_node(pipeline.CondorDAGNode):
  """
  """
  def __init__(self, job, dag, xml_list, id, p_node):
    pipeline.CondorDAGNode.__init__(self,job)
    self.add_macro("macroid", id)
    self.add_var_arg(xml_list)
    for p in p_node:
      self.add_parent(p)
    dag.add_node(self)

class summary_page_node(pipeline.CondorDAGNode):
  """
  """
  def __init__(self, job, dag, open_box, id, p_node):
    pipeline.CondorDAGNode.__init__(self,job)
    self.add_macro("macroid", id)
    if open_box: self.add_var_arg("open")
    for p in p_node:
      self.add_parent(p)
    dag.add_node(self)


def ifo_seg_dict(cp):
  out = {}
  
  out["H1"] = string.strip(cp.get('input','h1vetosegments'))
  out["H2"] = string.strip(cp.get('input','h2vetosegments'))
  out["L1"] = string.strip(cp.get('input','l1vetosegments'))
  cat = "_".join(os.path.basename(out["H1"]).split("_")[1:3])
  return out, [cat]


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

# get the segments for a given category veto
seg_dict, cats = ifo_seg_dict(cp)

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
lalappsNewcorseJobCombined = lalapps_newcorse_combined_job(cp)
ligolwSegmentsJob = ligolw_segments_job(cp)
ligolwThincaToCoincJob =  ligolw_thinca_to_coinc_job(cp)
hmUpperlimitJob = hm_upperlimit_job(cp)
hmUpperlimitPlotJob = ul_plot_job(cp)
farPlotJob = far_plot_job(cp)
summaryPageJob = summary_page_job(cp)

n = 0
#Do the segments node
segNode = {}
for cat in cats:
  segNode[cat] = ligolw_segments_node(ligolwSegmentsJob, dag, seg_dict, "vetoes", "vetoes_"+cat+".xml.gz", n); n+=1

#Some initialization
ligolwThincaToCoincNode = {}
sqliteNodeSimplify = {}
sqliteNodeRemoveH1H2 = {}
sqliteNodeCluster = {}
ligolwSqliteNode = {}
ligolwSqliteNodeInjDBtoXML = {}
ligolwSqliteNodeInjXMLtoDB = {}
ligolwInspinjfindNode = {}
lallappsNewcorseNode = {}
lallappsNewcorseNodeCombined = {}
hmUpperlimitNode = {}
hmUpperlimitPlotNode = {}
farPlotNode = {}
summaryPageNode = {}
db = {}

# to get injection file entries from the cache
injcache = map(lal.CacheEntry, file("inj.cache"))
inj = injcache[0]
start_time = inj.segment[0]
end_time = inj.segment[1]

timestr = str(start_time) + "-" + str(end_time)

for type in types:
  for cat in cats:
    #break down the cache to save on parsing
    tag = type + "_" + cat
    grep('THINCA.*_SECOND_.*'+type + ".*" + cat, FULLDATACACHE, tag + ".cache")
    try: os.mkdir(tag)
    except: pass
    ligolwThincaToCoincNode[type+cat] = ligolw_thinca_to_coinc_node(ligolwThincaToCoincJob, dag, tag+".cache", "vetoes_"+cat+".xml.gz", "vetoes", tag+"/S5_HM_"+timestr, n, start_time, end_time, effsnrfac=50, p_node=[segNode[cat]]); n+=1
    database = tag+"_"+timestr+".sqlite"
    try: db[cat].append(database) 
    except: db[cat] = [database]
    #xml_list = [tag+"/S5_HM_"+timestr+"*"+type+"*"+cat+"*.xml.gz", "vetoes_"+cat+".xml.gz"]
    xml_list = ["vetoes_"+cat+".xml.gz"]
    ligolwSqliteNode[type+cat] = ligolw_sqlite_node(ligolwSqliteJob, dag, database, xml_list, n, p_node=[ligolwThincaToCoincNode[type+cat]], replace=True, cache_pat=tag); n+=1
    sqliteNodeSimplify[type+cat] = sqlite_node(sqliteJob, dag, database, string.strip(cp.get('input',"simplify")), n, p_node=[ligolwSqliteNode[type+cat]]); n+=1
    sqliteNodeRemoveH1H2[type+cat] = sqlite_node(sqliteJob, dag, database, string.strip(cp.get('input',"remove_h1h2")),n, p_node=[sqliteNodeSimplify[type+cat]]); n+=1
    sqliteNodeCluster[type+cat] = sqlite_node(sqliteJob, dag, database, string.strip(cp.get('input',"cluster")),n, p_node=[sqliteNodeRemoveH1H2[type+cat]]); n+=1

for inj in injcache:
  for cat in cats:
    type = "_".join(inj.description.split("_")[2:])
    tag = type + "_" + cat
    url = inj.url
    cachefile = tag + ".cache"
    try: os.mkdir(tag)
    except: pass
    #break down the cache
    grep('THINCA_SECOND_.*'+type + '.*' + cat, INJCACHE, cachefile)
    ligolwThincaToCoincNode[type+cat] = ligolw_thinca_to_coinc_node(ligolwThincaToCoincJob, dag, cachefile, "vetoes_"+cat+".xml.gz", "vetoes", tag+"/S5_HM_INJ_"+timestr, n, start_time, end_time, effsnrfac=50, p_node=[segNode[cat]]);n+=1
    database = tag+"_"+timestr+".sqlite"
    db_to_xml_name = tag +"_"+timestr+".xml.gz"
    try: db[cat].append(database)
    except: db[cat] = [database]
    #xml_list = [type+cat+"/S5_HM_INJ_"+timestr+"*"+type+"*"+cat+"*.xml.gz", url, "vetoes_"+cat+".xml.gz"]
    xml_list = [url, "vetoes_"+cat+".xml.gz"]
    ligolwSqliteNode[type+cat] = ligolw_sqlite_node(ligolwSqliteJob, dag, database, xml_list, n, p_node=[ligolwThincaToCoincNode[type+cat]], replace=True,cache_pat=tag);n+=1
    sqliteNodeSimplify[type+cat] = sqlite_node(sqliteJob, dag, database, string.strip(cp.get('input',"simplify")), n, p_node=[ligolwSqliteNode[type+cat]]);n+=1
    sqliteNodeRemoveH1H2[type+cat] = sqlite_node(sqliteJob, dag, database, string.strip(cp.get('input',"remove_h1h2")),n, p_node=[sqliteNodeSimplify[type+cat]]);n+=1
    sqliteNodeCluster[type+cat] = sqlite_node(sqliteJob, dag, database, string.strip(cp.get('input',"cluster")),n, p_node=[sqliteNodeRemoveH1H2[type+cat]]);n+=1
    ligolwSqliteNodeInjDBtoXML[type+cat] = ligolw_sqlite_node(ligolwSqliteJob, dag, database, [db_to_xml_name], n, p_node=[sqliteNodeCluster[type+cat]], replace=False, extract=True); n+=1
    ligolwInspinjfindNode[type+cat] = ligolw_inspinjfind_node(ligolwInspinjfindJob, dag, db_to_xml_name, n, p_node=[ligolwSqliteNodeInjDBtoXML[type+cat]]);n+=1
    ligolwSqliteNodeInjXMLtoDB[type+cat] = ligolw_sqlite_node(ligolwSqliteJob, dag, database, [db_to_xml_name], n, p_node=[ligolwInspinjfindNode[type+cat]], replace=True);n+=1

# New corse
# First work out the parent/child relationships
p_nodes = []
for k in ligolwInspinjfindNode.keys():
  p_nodes.append(ligolwSqliteNodeInjXMLtoDB[k])
for k in sqliteNodeCluster.keys():
  p_nodes.append(sqliteNodeCluster[k])

#New corse jobs 
for cat in cats:
  #to compute uncombined far
  lallappsNewcorseNode[cat] = lalapps_newcorse_node(lalappsNewcorseJob, dag, "vetoes", " ".join(db[cat]), n, p_nodes, mass_bins="0,50,85,inf", categories="mtotal-ifos-oninstruments", rank="snr");n+=1
  #to compute combined far 
  lallappsNewcorseNodeCombined[cat] = lalapps_newcorse_node(lalappsNewcorseJobCombined, dag, "vetoes", " ".join(db[cat]), n, [lallappsNewcorseNode[cat]], mass_bins=None, categories="oninstruments", rank="uncombined-ifar");n+=1

#Upper limit jobs
for cat in cats: 
  hmUpperlimitNode[cat+"H1,H2,L1"] = hm_upperlimit_node(hmUpperlimitJob, dag, "H1,H2,L1",timestr, "FULL_DATA_CAT_3_"+timestr+".sqlite", "*INJ_CAT_3_"+timestr+".sqlite", 10000, "vetoes", n, p_node=[lallappsNewcorseNodeCombined[cat]]);n+=1
  hmUpperlimitPlotNode[cat+"H1,H2,L1"] = ul_plot_node(hmUpperlimitPlotJob, dag, '2Dsearchvolume-' + timestr + '-H1H2L1.xml', n, [hmUpperlimitNode[cat+"H1,H2,L1"]]);n+=1

  hmUpperlimitNode[cat+"H1,L1"] = hm_upperlimit_node(hmUpperlimitJob, dag, "H1,L1", timestr, "FULL_DATA_CAT_3_"+timestr+".sqlite", "*INJ_CAT_3_"+timestr+".sqlite", 10000, "vetoes", n, p_node=[lallappsNewcorseNodeCombined[cat]]);n+=1
  hmUpperlimitPlotNode[cat+"H1,L1"] = ul_plot_node(hmUpperlimitPlotJob, dag, '2Dsearchvolume-' + timestr + '-H1L1.xml', n, [hmUpperlimitNode[cat+"H1,L1"]]);n+=1

  hmUpperlimitNode[cat+"H2,L1"] = hm_upperlimit_node(hmUpperlimitJob, dag, "H2,L1", timestr, "FULL_DATA_CAT_3_"+timestr+".sqlite", "*INJ_CAT_3_"+timestr+".sqlite", 10000, "vetoes", n, p_node=[lallappsNewcorseNodeCombined[cat]]);n+=1
  hmUpperlimitPlotNode[cat+"H2,L1"] = ul_plot_node(hmUpperlimitPlotJob, dag, '2Dsearchvolume-' + timestr + '-H2L1.xml', n, [hmUpperlimitNode[cat+"H2,L1"]]);n+=1

#IFAR plots and combined upper limit plots and summary page
for cat in cats:
  farPlotNode[cat] = far_plot_node(farPlotJob, dag, " ".join(db[cat]), n, [lallappsNewcorseNodeCombined[cat]]);n+=1
  fstr = '2Dsearchvolume-' + timestr + '-H1H2L1.xml 2Dsearchvolume-' + timestr + '-H1L1.xml 2Dsearchvolume-' + timestr + '-H2L1.xml'
  hmUpperlimitPlotNode[cat] = ul_plot_node(hmUpperlimitPlotJob, dag, fstr, n, [hmUpperlimitNode[cat+"H1,H2,L1"],hmUpperlimitNode[cat+"H1,L1"],hmUpperlimitNode[cat+"H2,L1"]]);n+=1
  # Summary pages (open and closed box)
  summaryPageNode[cat] = summary_page_node(summaryPageJob, dag, False, n, [hmUpperlimitPlotNode[cat]]);n+=1
  summaryPageNode[cat+"open"] = summary_page_node(summaryPageJob, dag, True, n, [hmUpperlimitPlotNode[cat]]);n+=1


dag.write_sub_files()
dag.write_dag()
dag.write_script()

print "\nYour database output should be...\n"
for cat in cats:
  print cat + ":\n", " ".join(db[cat]) + "\n"

