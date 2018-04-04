"""
This module contains condor jobs / node classes for the followup dag

This program creates cache files for the output of inspiral hipe
"""

__author__ = 'Chad Hanna <channa@phys.lsu.edu>'

##############################################################################
# import standard modules and append the lalapps prefix to the python path
import sys, os, copy, math
from subprocess import *
import socket, time
import re, string
from optparse import *
import tempfile
import ConfigParser 
import urlparse
from UserDict import UserDict
sys.path.append('@PYTHONLIBDIR@')

##############################################################################
# import the modules we need to build the pipeline
from glue import pipeline
from glue import lal
from glue import segments
from glue import segmentsUtils
from pylal.webUtils import *
from lalapps import inspiral
  
###### GENERIC CLASSES TO HELP WEBIFY A CONDOR DAG OUTPUT ####################
##############################################################################

class webTheJob:
  """
  webTheJob is a class intended to be inherited by a class that subclasses
  the condor DAG Job.  It is useful for setting up a standard structure to
  webify the output of a dag.  You'll want to use webTheNode and webTheDAG
  for your condor subclasses too. 
  """

  def __init__(self):
    pass

  def setupJobWeb(self, name, tag_base=None, cp=None):
    # Give this job a name.  Then make directories for the log files and such
    # This name is important since these directories will be included in
    # the web tree.
    self.name = name
    if not os.path.exists(name):
       os.mkdir(name)
    if not os.path.exists(name+'/logs'):
       os.mkdir(name+'/logs')
    if not os.path.exists(name+'/Images'):
       os.mkdir(name+'/Images')
    # Set up the usual stuff and name the log files appropriately
    self.tag_base = tag_base
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
    self.set_sub_file(name+'.sub')
    self.relPath = name + '/'
    self.outputPath = os.getcwd() + '/' + name + '/'
    self.set_stdout_file(self.outputPath+'/logs/'+name+'-$(macroid).out')
    self.set_stderr_file(self.outputPath+'/logs/'+name+'-$(macroid).err')
    if cp:
      if cp.has_section("condor-memory-requirement") and \
      		cp.has_option("condor-memory-requirement",name):
        requirement = cp.getint("condor-memory-requirement",name)
        self.add_condor_cmd("Requirements", \
		"(Memory > " + str(requirement) + ")")

class webTheNode:
  """
  webTheNode is a class intended to be inherited by a class that subclasses
  the condor DAG Node .  It is useful for setting up a standard structure to
  webify the output of a dag.  You'll want to use webTheJob and webTheDAG
  for your condor subclasses too. 
  """

  def __init__(self):
    pass

  def setupNodeWeb(self, job, passItAlong=False, content=None, page=None,webOverride=None,cache=None):
    # setup the node id
    self.add_macro("macroid", self.id)
    # determine the output web file name for the job
    #self.webFileName = job.outputPath + self.id + '.html'
    self.jobName = job.name
    #if page:
    #  self.webLink = page+'/'+job.relPath+self.id+'.html'
    #if webOverride:
    #  self.webLink = webOverride
    # standardize the I/O for executables that themselves write web pages.
    # this is great for developers since they can in some sense treat
    # the DAG as a black box for their executable.   They just have to 
    # comply with these few command line arguments to tie everything together
    if passItAlong:
      #self.add_var_opt("output-web-file",self.webFileName)
      self.add_var_opt("output-path",job.outputPath)
      self.add_var_opt("enable-output","")
      #self.add_var_opt("page-rel-path",job.relPath)
      #self.add_var_opt("page", page)
    #if content: self.writeContent(content)
    if cache:
      cache.appendCache(job.name,job.outputPath)
      try:
        if self.outputCache:
          cache.appendSubCache(job.name,self.outputCache)
      except: pass
      #if passItAlong:
        #output_cache = self.id.replace('-',' ') + ' ' + self.webFileName + '\n' + self.id.replace('-',' ') + ' ' + self.webFileName.replace('html','ini') + '\n'

        #cache.appendSubCache(job.name,output_cache)

  def writeContent(self,content):
    # The talkBack class is a way for the users job to provide information
    # back to the DAG web.  It is done through the reading and writing of
    # a config file (.ini) that has the same naming convention as the
    # web file
    #self.talkBack = talkBack(self.webFileName)
    #self.talkBack.read()
    content.appendTable(1,2,0,700)
    self.webTable = content.lastTable;
    content.lastTable.row[0].cell[0].link(self.webLink,self.friendlyName)
    #self.talkBack.readAppend(content.lastTable.row[0].cell[1])
    # Each time the dag is generated it checks for the existance of the
    # appropriate config file to include the contents in the web page
    #if 0 and self.talkBack.summaryText:
    #  content.lastTable.row[0].cell[0].linebreak()
    #  content.lastTable.row[0].cell[0].text(self.talkBack.summaryText)
    #if 0 and self.talkBack.summaryPlot:
    #  content.lastTable.row[0].cell[1].image(self.talkBack.summaryPlot)
    #if 0 and self.talkBack.summaryPlotCaption:
    #  content.lastTable.row[0].cell[1].linebreak()
    #  content.lastTable.row[0].cell[1].text(self.talkBack.summaryPlotCaption)
 
  def set_id(self, id):
    self.id = id
  
  def validate(self):
    self.validNode = True

  def invalidate(self):
    self.validNode = False

class webTheDAG:
  """
  webTheDAG is a class intended to be inherited by a class that subclasses
  the condor DAG Node .  It is useful for setting up a standard structure to
  webify the output of a dag.  You'll want to use webTheJob and webTheDAG
  for your condor subclasses too. 
  """

  def __init__(self):
    pass

  def setupDAGWeb(self,title,filename,cp,opts):
    self.publish_path = string.strip(cp.get('followup-output','page'))
    self.page = string.strip(cp.get('followup-output','url'))
    self.webPage = WebPage(title,filename,self.page)
    self.webDirs = {}
    self.cache = cacheStructure()
    try:
       os.mkdir('DAGWeb')
    except: pass
    if not opts.disable_dag_categories:
        for cp_opt in cp.options('condor-max-jobs'):
          self.add_maxjobs_category(cp_opt,cp.getint('condor-max-jobs',cp_opt))


  def writeDAGWeb(self,type):
    #self.webPage.cleanWrite(type)
    self.cache.writeCache()

  def appendSection(self,name):
    self.webPage.appendSection(name)
    inifile = name.replace(" ","_").replace("@","-").replace("=",'-') + '.ini'
    file = open('DAGWeb/'+inifile,'a')
    file.close()
#    talkback = talkBack('DAGWeb/'+inifile)
#    talkback.readAppend(self.webPage.lastSection)
#    talkback.read()
#    self.webPage.lastSection.appendTable(1,2,0,700)

#    if talkback.summaryText:
#      self.webPage.lastSection.lastTable.row[0].cell[0].linebreak()
#      self.webPage.lastSection.lastTable.row[0].cell[0].text(talkback.summaryText)
#    if talkback.summaryPlot:
#      self.webPage.lastSection.lastTable.row[0].cell[1].image(talkback.summaryPlot)
#    if talkback.summaryPlotCaption:
#      self.webPage.lastSection.lastTable.row[0].cell[1].linebreak()
#      self.webPage.lastSection.lastTable.row[0].cell[1].text(talkback.summaryPlotCaption)

  def appendSubSection(self,name):
    self.webPage.lastSection.appendSubSection(name)
    inifile = name.replace(" ","_").replace("@","-").replace("=",'-') + '.ini'
    file = open('DAGWeb/'+inifile,'a')
    file.close()
#    talkback = talkBack('DAGWeb/'+inifile)
#    talkback.readAppend(self.webPage.lastSection.lastSub)
#    talkback.read()
#    self.webPage.lastSection.lastSub.appendTable(1,2,0,700)

#    if talkback.summaryText:
#      self.webPage.lastSection.lastSub.lastTable.row[0].cell[0].linebreak()
#      self.webPage.lastSection.lastSub.lastTable.row[0].cell[0].text(talkback.summaryText)
#    if talkback.summaryPlot:
#      self.webPage.lastSection.lastSub.lastTable.row[0].cell[1].image(talkback.summaryPlot)
#    if talkback.summaryPlotCaption:
#      self.webPage.lastSection.lastSub.lastTable.row[0].cell[1].linebreak()
#      self.webPage.lastSection.lastSub.lastTable.row[0].cell[1].text(talkback.summaryPlotCaption)

  def addNode(self, node,jobType):
    try:
      self.jobsDict[jobType] = self.jobsDict[jobType] + 1
      self.webDirs[node.jobName] = node.jobName
    except:
      self.jobsDict[jobType] = 1
    self.add_node(node)


  def publishToHydra(self):
    dirStr = ''
    for dir in self.webDirs:
      dirStr += dir + ' '
    dirStr = 'rsync -vrz '+dirStr+' DAGWeb index.html '
    print dirStr
    copying_results = call(dirStr+self.publish_path, shell=True)
    if copying_results != 0:
      print >> sys.stderr, "the followup results could not be copied to "+self.publish_path
      sys.exit(1)

  def printNodeCounts(self):
    for jobs in self.jobsDict:
      print "\nFound " + str(self.jobsDict[jobs]) + " " + str(jobs) + " Jobs"

  def writeAll(self, type='IUL'):
    self.printNodeCounts()
    print "\n\n.......Writing DAG"
    self.write_sub_files()
    self.write_dag()
    self.write_script()
    self.writeDAGWeb(type)
    print "\n\n  Created a DAG file which can be submitted by executing"
    print "    condor_submit_dag " + self.get_dag_file()
    print """\n  from a condor submit machine
  Before submitting the dag, you must execute

    export _CONDOR_DAGMAN_LOG_ON_NFS_IS_ERROR=FALSE

  If you are running LSCdataFind jobs, do not forget to initialize your grid
  proxy certificate on the condor submit machine by running the commands

    unset X509_USER_PROXY
    grid-proxy-init -hours 72

  Enter your pass phrase when prompted. The proxy will be valid for 72 hours.
  If you expect the LSCdataFind jobs to take longer to complete, increase the
  time specified in the -hours option to grid-proxy-init. You can check that
  the grid proxy has been sucessfully created by executing the command:

    grid-cert-info -all -file /tmp/x509up_u`id -u`

  This will also give the expiry time of the proxy."""

