"""
This modules contains objects that make it simple for the user to 
create python scripts that build Condor DAGs to run code on the LSC
Data Grid.
"""

__author__ = 'Duncan Brown <duncan@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

import os
import sys
import string, re
import exceptions
import time
import random
import md5
import math

def s2play(t):
  """
  Return 1 if t is in the S2 playground, 0 otherwise
  t = GPS time to test if playground
  """
  if ((t - 729273613) % 6370) < 600:
    return 1
  else:
    return 0


class CondorError(exceptions.Exception):
  """Error thrown by Condor Jobs"""
  def __init__(self, args=None):
    self.args = args
class CondorJobError(CondorError):
  pass
class CondorSubmitError(CondorError):
  pass
class CondorDAGError(CondorError):
  pass
class CondorDAGNodeError(CondorError):
  pass
class SegmentError(exceptions.Exception):
  def __init__(self, args=None):
    self.args = args


class CondorJob:
  """
  Generic condor job class. Provides methods to set the options in the
  condor submit file for a particular executable
  """
  def __init__(self, universe, executable, queue):
    """
    @param universe: the condor universe to run the job in.
    @param executable: the executable to run.
    @param queue: number of jobs to queue.
    """
    self.__universe = universe
    self.__executable = executable
    self.__queue = queue

    # These are set by methods in the class
    self.__options = {}
    self.__short_options = {}
    self.__arguments = []
    self.__condor_cmds = {}
    self.__notification = None
    self.__log_file = None
    self.__err_file = None
    self.__out_file = None
    self.__sub_file_path = None

  def add_condor_cmd(self, cmd, value):
    """
    Add a Condor command to the submit file (e.g. a class add or evironment).
    @param cmd: Condor command directive.
    @param value: value for command.
    """
    self.__condor_cmds[cmd] = value

  def add_arg(self, arg):
    """
    Add an argument to the executable. Arguments are appended after any
    options and their order is guaranteed.
    @param arg: argument to add.
    """
    self.__arguments.append(arg)

  def get_args(self):
    """
    Return the list of arguments that are to be passed to the executable.
    """

    return self.__arguments

  def add_opt(self, opt, value):
    """
    Add a command line option to the executable. The order that the arguments
    will be appended to the command line is not guaranteed, but they will
    always be added before any command line arguments. The name of the option
    is prefixed with double hyphen and the program is expected to parse it
    with getopt_long().
    @param opt: command line option to add.
    @param value: value to pass to the option (None for no argument).
    """
    self.__options[opt] = value

  def get_opts(self):
    """
    Return the dictionary of opts for the job.
    """

    return self.__options

  def add_short_opt(self, opt, value):
    """
    Add a command line option to the executable. The order that the arguments
    will be appended to the command line is not guaranteed, but they will
    always be added before any command line arguments. The name of the option
    is prefixed with single hyphen and the program is expected to parse it
    with getopt() or getopt_long() (if a single character option), or
    getopt_long_only() (if multiple characters).  Long and (single-character)
    short options may be mixed if the executable permits this.
    @param opt: command line option to add.
    @param value: value to pass to the option (None for no argument).
    """
    self.__short_options[opt] = value

  def get_short_opts(self):
    """
    Return the dictionary of short options for the job.
    """

    return self.__short_options

  def add_ini_opts(self, cp, section):
    """
    Parse command line options from a given section in an ini file and
    pass to the executable.
    @param cp: ConfigParser object pointing to the ini file.
    @param section: section of the ini file to add to the options.
    """
    for opt in cp.options(section):
      arg = string.strip(cp.get(section,opt))
      self.__options[opt] = arg

  def set_notification(self, value):
    """
    Set the email address to send notification to.
    @param value: email address or never for no notification.
    """
    self.__notification = value

  def set_log_file(self, path):
    """
    Set the Condor log file.
    @param path: path to log file.
    """
    self.__log_file = path
    
  def set_stderr_file(self, path):
    """
    Set the file to which Condor directs the stderr of the job.
    @param path: path to stderr file.
    """
    self.__err_file = path

  def get_stderr_file(self):
    """
    Get the file to which Condor directs the stderr of the job.
    """
    return self.__err_file

  def set_stdout_file(self, path):
    """
    Set the file to which Condor directs the stdout of the job.
    @param path: path to stdout file.
    """
    self.__out_file = path
  
  def get_stdout_file(self):
    """
    Get the file to which Condor directs the stdout of the job.
    """
    return self.__out_file

  def set_sub_file(self, path):
    """
    Set the name of the file to write the Condor submit file to when
    write_sub_file() is called.
    @param path: path to submit file.
    """
    self.__sub_file_path = path

  def get_sub_file(self):
    """
    Get the name of the file which the Condor submit file will be
    written to when write_sub_file() is called.
    """
    return self.__sub_file_path

  def write_sub_file(self):
    """
    Write a submit file for this Condor job.
    """
    if not self.__log_file:
      raise CondorSubmitError, "Log file not specified."
    if not self.__err_file:
      raise CondorSubmitError, "Error file not specified."
    if not self.__out_file:
      raise CondorSubmitError, "Output file not specified."
  
    if not self.__sub_file_path:
      raise CondorSubmitError, 'No path for submit file.'
    try:
      subfile = open(self.__sub_file_path, 'w')
    except:
      raise CondorSubmitError, "Cannot open file " + self.__sub_file_path

    subfile.write( 'universe = ' + self.__universe + '\n' )
    subfile.write( 'executable = ' + self.__executable + '\n' )

    if self.__options.keys() or self.__short_options.keys() or self.arguments:
      subfile.write( 'arguments =' )
      for c in self.__options.keys():
        if self.__options[c]:
          subfile.write( ' --' + c + ' ' + self.__options[c] )
        else:
          subfile.write( ' --' + c )
      for c in self.__short_options.keys():
        if self.__short_options[c]:
          subfile.write( ' -' + c + ' ' + self.__short_options[c] )
        else:
          subfile.write( ' -' + c )
      for c in self.__arguments:
        subfile.write( ' ' + c )
      subfile.write( '\n' )

    for cmd in self.__condor_cmds.keys():
      subfile.write( cmd + " = " + self.__condor_cmds[cmd] + '\n' )

    subfile.write( 'log = ' + self.__log_file + '\n' )
    subfile.write( 'error = ' + self.__err_file + '\n' )
    subfile.write( 'output = ' + self.__out_file + '\n' )
    if self.__notification:
      subfile.write( 'notification = ' + self.__notification + '\n' )
    subfile.write( 'queue ' + str(self.__queue) + '\n' )

    subfile.close()



class CondorDAGJob(CondorJob):
  """
  A Condor DAG job never notifies the user on completion and can have variable
  options that are set for a particular node in the DAG. Inherits methods
  from a CondorJob.
  """
  def __init__(self, universe, executable):
    """
    universe = the condor universe to run the job in.
    executable = the executable to run in the DAG.
    """
    CondorJob.__init__(self, universe, executable, 1)
    CondorJob.set_notification(self, 'never')
    self.__var_opts = []
    self.__have_var_args = 0
    self.__bad_macro_chars = re.compile(r'[_-]')

  def add_var_opt(self, opt):
    """
    Add a variable (or macro) option to the condor job. The option is added 
    to the submit file and a different argument to the option can be set for
    each node in the DAG.
    @param opt: name of option to add.
    """
    if opt not in self.__var_opts:
      self.__var_opts.append(opt)
      macro = self.__bad_macro_chars.sub( r'', opt )
      self.add_opt(opt,'$(macro' + macro + ')')

  def add_var_arg(self):
    """
    Add a command to the submit file to allow variable (macro) arguments
    to be passed to the executable.
    """
    if not self.__have_var_args:
      self.add_arg('$(macroarguments)')
      self.__have_var_args = 1



class CondorDAGNode:
  """
  A CondorDAGNode represents a node in the DAG. It corresponds to a particular
  condor job (and so a particular submit file). If the job has variable
  (macro) options, they can be set here so each nodes executes with the
  correct options.
  """
  def __init__(self, job):
    """
    @param job: the CondorJob that this node corresponds to.
    """
    if not isinstance(job, CondorDAGJob):
      raise CondorDAGNodeError, "A DAG node must correspond to a Condor DAG job"
    self.__name = None
    self.__job = job
    self.__pre_script = None
    self.__pre_script_args = []
    self.__post_script = None
    self.__post_script_args = []
    self.__macros = {}
    self.__opts = {}
    self.__args = []
    self.__retry = 0
    self.__parents = []
    self.__bad_macro_chars = re.compile(r'[_-]')
    self.__output_files = []
    self.__input_files = []
    self.set_name()

  def __repr__(self):
    return self.__name

  def job(self):
    """
    Return the CondorJob that this node is associated with.
    """
    return self.__job
  
  def set_pre_script(self,script):
    """
    Sets the name of the pre script that is executed before the DAG node is
    run.
    @param script: path to script
    """
    self.__pre_script = script

  def add_pre_script_arg(self,arg):
    """
    Adds an argument to the pre script that is executed before the DAG node is
    run.
    """
    self.__pre_script_args.append(arg)

  def set_post_script(self,script):
    """
    Sets the name of the post script that is executed before the DAG node is
    run.
    @param script: path to script
    """
    self.__post_script = script

  def add_post_script_arg(self,arg):
    """
    Adds an argument to the post script that is executed before the DAG node is
    run.
    """
    self.__post_script_args.append(arg)

  def set_name(self):
    """
    Generate a unique name for this node in the DAG.
    """
    t = str( long( time.time() * 1000 ) )
    r = str( long( random.random() * 100000000000000000L ) )
    a = str( self.__class__ )
    self.__name = md5.md5(t + r + a).hexdigest()

  def add_input_file(self, filename):
    """
    Add filename as a necessary input file for this DAG node.

    @param filename: input filename to add
    """
    if filename not in self.__input_files:
        self.__input_files.append(filename)

  def add_output_file(self, filename):
    """
    Add filename as a output file for this DAG node.

    @param filename: output filename to add
    """
    if filename not in self.__output_files:
        self.__output_files.append(filename)

  def get_input_files(self):
    """
    Return list of input files for this DAG node.
    """

    return self.__input_files

  def get_output_files(self):
    """
    Return list of output files for this DAG node.
    """

    return self.__output_files


  def add_macro(self,name,value):
    """
    Add a variable (macro) for this node.  This can be different for
    each node in the DAG, even if they use the same CondorJob.  Within
    the CondorJob, the value of the macro can be referenced as
    '$(name)' -- for instance, to define a unique output or error file
    for each node.
    @param name: macro name.
    @param value: value of the macro for this node in the DAG
    """
    macro = self.__bad_macro_chars.sub( r'', name )
    self.__opts[macro] = value

  def get_opts(self):
    """
    Return the opts for this node. Note that this returns only
    the options for this instance of the node and not those
    associated with the underlying job template.
    """

    return self.__opts

  def add_var_opt(self,opt,value):
    """
    Add a variable (macro) option for this node. If the option
    specified does not exist in the CondorJob, it is added so the submit
    file will be correct when written.
    @param opt: option name.
    @param value: value of the option for this node in the DAG.
    """
    macro = self.__bad_macro_chars.sub( r'', opt )
    self.__opts['macro' + macro] = value
    self.__job.add_var_opt(opt)

  def add_var_arg(self, arg):
    """
    Add a variable (or macro) argument to the condor job. The argument is
    added to the submit file and a different value of the argument can be set
    for each node in the DAG.
    @param arg: name of option to add.
    """
    self.__args.append(arg)
    self.__job.add_var_arg()

  def get_args(self):
    """
    Return the arguments for this node. Note that this returns
    only the arguments for this instance of the node and not those
    associated with the underlying job template.
    """

    return self.__args

  def set_retry(self, retry):
    """
    Set the number of times that this node in the DAG should retry.
    @param retry: number of times to retry node.
    """
    self.__retry = retry

  def write_job(self,fh):
    """
    Write the DAG entry for this node's job to the DAG file descriptor.
    @param fh: descriptor of open DAG file.
    """
    fh.write( 'JOB ' + self.__name + ' ' + self.__job.get_sub_file() +  '\n' )
    fh.write( 'RETRY ' + self.__name + ' ' + str(self.__retry) + '\n' )

  def write_vars(self,fh):
    """
    Write the variable (macro) options and arguments to the DAG file
    descriptor.
    @param fh: descriptor of open DAG file.
    """
    if self.__macros.keys() or self.__opts.keys() or self.__args:
      fh.write( 'VARS ' + self.__name )
    for k in self.__macros.keys():
      fh.write( ' ' + str(k) + '="' + str(self.__macros[k]) + '"' )
    for k in self.__opts.keys():
      fh.write( ' ' + str(k) + '="' + str(self.__opts[k]) + '"' )
    if self.__args:
      fh.write( ' macroarguments="' + ' '.join(self.__args) + '"' )
    fh.write( '\n' )

  def write_parents(self,fh):
    """
    Write the parent/child relations for this job to the DAG file descriptor.
    @param fh: descriptor of open DAG file.
    """
    for parent in self.__parents:
      fh.write( 'PARENT ' + parent + ' CHILD ' + str(self) + '\n' )

  def write_pre_script(self,fh):
    """
    Write the pre script for the job, if there is one
    @param fh: descriptor of open DAG file.
    """
    if self.__pre_script:
      fh.write( 'SCRIPT PRE ' + str(self) + ' ' + self.__pre_script + ' ' +
        ' '.join(self.__pre_script_args) + '\n' )

  def write_post_script(self,fh):
    """
    Write the post script for the job, if there is one
    @param fh: descriptor of open DAG file.
    """
    if self.__post_script:
      fh.write( 'SCRIPT POST ' + str(self) + ' ' + self.__post_script + ' ' +
        ' '.join(self.__post_script_args) + '\n' )

  def write_input_files(self, fh):
    """
    Write as a comment into the DAG file the list of input files
    for this DAG node.

    @param fh: descriptor of open DAG file.
    """

    for f in self.__input_files:
       print >>fh, "## Job %s requires input file %s" % (self.__name, f)
 
  def write_output_files(self, fh):
    """
    Write as a comment into the DAG file the list of output files
    for this DAG node.

    @param fh: descriptor of open DAG file.
    """

    for f in self.__output_files:
       print >>fh, "## Job %s generates output file %s" % (self.__name, f)

  def set_log_file(self,log):
    """
    Set the Condor log file to be used by this CondorJob.
    @param log: path of Condor log file.
    """
    self.__job.set_log_file(log)

  def add_parent(self,node):
    """
    Add a parent to this node. This node will not be executed until the
    parent node has run sucessfully.
    @param node: CondorDAGNode to add as a parent.
    """
    if not isinstance(node, CondorDAGNode):
      raise CondorDAGNodeError, "Parent must be a Condor DAG node"
    self.__parents.append( str(node) )

  def get_cmd_line(self):
    """
    Return the full command line that will be used when this node
    is run by DAGman.
    """

    # pattern to find DAGman macros
    pat = re.compile(r'\$\((.+)\)')

    # first parse the options and replace macros with values
    options = self.job().get_opts()
    macros = self.get_opts()

    cmd = ""

    for k in options:
        val = options[k]
        m = pat.match(val)
        if m:
            key = m.group(1)
            value = macros[key]

            cmd += "--%s %s " % (k, value)
        else:
            cmd += "--%s %s " % (k, val)

    # second parse the short options and replace macros with values
    options = self.job().get_short_opts()

    for k in options:
        val = options[k]
        m = pat.match(val)
        if m:
            key = m.group(1)
            value = macros[key]

            cmd += "-%s %s " % (k, value)
        else:
            cmd += "-%s %s " % (k, val)

    # lastly parse the arguments and replace macros with values
    args = self.job().get_args()
    macros = self.get_args()

    for k in args:
        val = args[k]
        m = pat.match(val)
        if m:
            key = m.group(1)
            value = macros[key]

            cmd += "%s " % (k, value)
        else:
            cmd += "%s " % (k, val)

    return cmd
    


class CondorDAG:
  """
  A CondorDAG is a Condor Directed Acyclic Graph that describes a collection
  of Condor jobs and the order in which to run them. All Condor jobs in the
  DAG must write their Codor logs to the same file.
  NOTE: The log file must not be on an NFS mounted system as the Condor jobs
  must be able to get an exclusive file lock on the log file.
  """
  def __init__(self,log):
    """
    @param log: path to log file which must not be on an NFS mounted file system.
    """
    self.__log_file_path = log
    self.__dag_file_path = None
    self.__jobs = []
    self.__nodes = []

  def set_dag_file(self, path):
    """
    Set the name of the file into which the DAG is written.
    @param path: path to DAG file.
    """
    self.__dag_file_path = path

  def get_dag_file(self):
    """
    Return the path to the DAG file.
    """
    if not self.__log_file_path:
      raise CondorDAGError, "No path for DAG file"
    else:
      return self.__dag_file_path

  def add_node(self,node):
    """
    Add a CondorDAGNode to this DAG. The CondorJob that the node uses is 
    also added to the list of Condor jobs in the DAG so that a list of the
    submit files needed by the DAG can be maintained. Each unique CondorJob
    will be added once to prevent duplicate submit files being written.
    @param node: CondorDAGNode to add to the CondorDAG.
    """
    if not isinstance(node, CondorDAGNode):
      raise CondorDAGError, "Nodes must be class CondorDAGNode or subclass"
    node.set_log_file(self.__log_file_path)
    self.__nodes.append(node)
    if node.job() not in self.__jobs:
      self.__jobs.append(node.job())

  def write_sub_files(self):
    """
    Write all the submit files used by the dag to disk. Each submit file is
    written to the file name set in the CondorJob.
    """
    for job in self.__jobs:
      job.write_sub_file()

  def write_dag(self):
    """
    Write all the nodes in the DAG to the DAG file.
    """
    if not self.__dag_file_path:
      raise CondorDAGError, "No path for DAG file"
    try:
      dagfile = open( self.__dag_file_path, 'w' )
    except:
      raise CondorDAGError, "Cannot open file " + self.__dag_file_path
    for node in self.__nodes:
      node.write_job(dagfile)
      node.write_vars(dagfile)
      node.write_pre_script(dagfile)
      node.write_post_script(dagfile)
    for node in self.__nodes:
      node.write_parents(dagfile)
    dagfile.close()

  def write_dax(self):
    """
    Write all the nodes in the workflow to the DAX file.
    """
    if not self.__dag_file_path:
      raise CondorDAGError, "No path for DAX file"
    try:
      dagfile = open( self.__dag_file_path, 'w' )
    except:
      raise CondorDAGError, "Cannot open file " + self.__dag_file_path

    # write the preamble
    preamble = """\
<?xml version="1.0" encoding="UTF-8"?>
<adag xmlns="http://www.griphyn.org/chimera/DAX"
        xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" 
        xsi:schemaLocation="http://www.griphyn.org/chimera/DAX
        http://www.griphyn.org/chimera/dax-1.8.xsd"
        name="inspiral" index="0" count="1" version="1.8">
"""

    print >>dagfile, preamble

    # find unique input and output files from nodes
    input_file_dict = {}
    output_file_dict = {}
 
    for node in self.__nodes:
      input_files = node.get_input_files()
      output_files = node.get_output_files()
      for f in input_files:
         input_file_dict[f] = 1
      for f in output_files:
         output_file_dict[f] = 1

    # move union of input and output into inout
    inout_file_dict = {}

    for f in input_file_dict:
      if output_file_dict.has_key(f):
         inout_file_dict[f] = 1

    for f in inout_file_dict:
         del input_file_dict[f]
         del output_file_dict[f]

    # print input, inout, and output to dax
    filelist = input_file_dict.keys()
    filelist.sort()
    for f in filelist:
           msg = """\
    <filename file="%s" link="input"/>\
"""
           print >>dagfile, msg % f

    filelist = inout_file_dict.keys()
    filelist.sort()
    for f in filelist:
           msg = """\
    <filename file="%s" link="inout"/>\
"""
           print >>dagfile, msg % f

    filelist = output_file_dict.keys()
    filelist.sort()
    for f in filelist:
           msg = """\
    <filename file="%s" link="output"/>\
"""
           print >>dagfile, msg % f

    # write the jobs themselves to the DAX, making sure
    # to replace logical file references by the appropriate
    # xml, and adding the files used by each job both for
    # input and output

    # we save the ID number to DAG node name mapping so that
    # we can easily write out the child/parent relationship
    # later
    node_name_id_dict = {}

    id = 0
    for node in self.__nodes:
        executable = node.job()._CondorJob__executable
        node_name = node._CondorDAGNode__name

        id += 1
        id_tag = "ID%06d" % id
        node_name_id_dict[node_name] = id_tag

        cmd_line = node.get_cmd_line()
        
        # loop through all filenames looking for them in the command
        # line so that they can be replaced appropriately by xml tags
        for f in node.get_input_files():
            xml = '<filename file="%s" />' % f
            cmd_line = cmd_line.replace(f, xml)
        for f in node.get_output_files():      
            xml = '<filename file="%s" />' % f
            cmd_line = cmd_line.replace(f, xml)

        template = """\
<job id="%s" namespace="ligo" name="%s" version="1.0" level="1" dv-name="%s">
     <argument>%s
     </argument>\
"""
        xml = template % (id_tag, executable, node_name, cmd_line)

        print >>dagfile, xml

        for f in node.get_input_files():
                print >>dagfile, """\
     <uses file="%s" link="input" dontRegister="false" dontTransfer="false"/>\
""" % f

        for f in node.get_output_files():
                print >>dagfile, """\
     <uses file="%s" link="output" dontRegister="false" dontTransfer="false"/>\
""" % f

        print >>dagfile, "</job>"

    # print parent-child relationships to DAX
    
    for node in self.__nodes:
        child_id = node_name_id_dict[str(node)]
        if node._CondorDAGNode__parents:
                print >>dagfile, '<child ref="%s">' % child_id
                for parent in node._CondorDAGNode__parents:
                    parent_id = node_name_id_dict[parent]
                    print >>dagfile, '     <parent ref="%s"/>' % parent_id
                print >>dagfile, '</child>'


    print >>dagfile, "</adag>"

    dagfile.close()


class AnalysisJob:
  """
  Describes a generic analysis job that filters LIGO data as configured by
  an ini file.
  """
  def __init__(self,cp):
    """
    @param cp: ConfigParser object that contains the configuration for this job.
    """
    self.__cp = cp
    self.__channel = string.strip(self.__cp.get('input','channel'))

  def get_config(self,sec,opt):
    """
    Get the configration variable in a particular section of this jobs ini
    file.
    @param sec: ini file section.
    @param opt: option from section sec.
    """
    return string.strip(self.__cp.get(sec,opt))

  def channel(self):
    """
    Returns the name of the channel that this job is filtering. Note that 
    channel is defined to be IFO independent, so this may be LSC-AS_Q or
    IOO-MC_F. The IFO is set on a per node basis, not a per job basis.
    """
    return self.__channel



class AnalysisNode(CondorDAGNode):
  """
  Contains the methods that allow an object to be built to analyse LIGO
  data in a Condor DAG.
  """
  def __init__(self):
    self.__start = 0
    self.__end = 0
    self.__ifo = None
    self.__ifo_tag = None
    self.__input = None
    self.__output = None
    self.__calibration = None
    self.__calibration_cache_path = None
    self.__LHO2k = re.compile(r'H2')

  def set_start(self,time):
    """
    Set the GPS start time of the analysis node by setting a --gps-start-time
    option to the node when it is executed.
    @param time: GPS start time of job.
    """
    self.add_var_opt('gps-start-time',time)
    self.__start = time
    #if not self.__calibration and self.__ifo and self.__start > 0:
    #  self.calibration()

  def get_start(self):
    """
    Get the GPS start time of the node.
    """
    return self.__start
    
  def set_end(self,time):
    """
    Set the GPS end time of the analysis node by setting a --gps-end-time
    option to the node when it is executed.
    @param time: GPS end time of job.
    """
    self.add_var_opt('gps-end-time',time)
    self.__end = time

  def get_end(self):
    """
    Get the GPS end time of the node.
    """
    return self.__end

  def set_input(self,file):
    """
    Add an input to the node by adding a --input option.
    @param file: option argument to pass as input.
    """
    self.__input = file
    self.add_var_opt('input', file)
    self.add_input_file(file)

  def get_input(self):
    """
    Get the file that will be passed as input.
    """
    return self.__input

  def set_output(self, file):
    """
    Add an output to the node by adding a --output option.
    @param file: option argument to pass as output.
    """
    self.__output = file
    self.add_var_opt('output', file)
    self.add_output_file(file)

  def get_output(self):
    """
    Get the file that will be passed as output.
    """
    return self.__output

  def set_ifo(self,ifo):
    """
    Set the channel name to analyze and add a calibration file for that
    channel. The name of the ifo is prepended to the channel name obtained
    from the job configuration file and passed with a --channel-name option.
    A calibration file is obtained from the ini file and passed with a 
    --calibration-cache option.
    @param ifo: two letter ifo code (e.g. L1, H1 or H2).
    """
    self.__ifo = ifo
    self.add_var_opt('channel-name', ifo + ':' + self.job().channel())
    #if not self.__calibration and self.__ifo and self.__start > 0:
    #   self.calibration()

  def get_ifo(self):
    """
    Returns the two letter IFO code for this node.
    """
    return self.__ifo

  def set_ifo_tag(self,ifo_tag):
    """
    Set the ifo tag that is passed to the analysis code.
    @param ifo_tag: a string to identify one or more IFOs
    """
    self.__ifo_tag = ifo_tag
    self.add_var_opt('ifo-tag', ifo_tag)

  def get_ifo_tag(self):
    """
    Returns the IFO tag string
    """
    return self.__ifo_tag

  def set_cache(self,file):
    """
    Set the LAL frame cache to to use. The frame cache is passed to the job
    with the --frame-cache argument.
    @param file: calibration file to use.
    """
    self.add_var_opt('frame-cache', file)
    self.add_input_file(file)

  def calibration_cache_path(self):
    """
    Determine the path to the correct calibration cache file to use.
    """
    if self.__ifo and self.__start > 0:
        cal_path = self.job().get_config('calibration','path')

        if ( self.__LHO2k.match(self.__ifo) and 
          (self.__start >= 729273613) and (self.__start <= 734367613) ):
          if self.__start < int(
            self.job().get_config('calibration','H2-cal-epoch-boundary')):
            cal_file = self.job().get_config('calibration','H2-1')
          else:
            cal_file = self.job().get_config('calibration','H2-2')
        else:
          cal_file = self.job().get_config('calibration',self.__ifo)

        cal = os.path.join(cal_path,cal_file)
        self.__calibration_cache_path = cal
    else:
       msg = "__ifo and __start attributes must be set first"
       raise CondorDAGNodeError, msg 

  def calibration(self):
    """
    Set the path to the calibration cache file for the given IFO.
    During S2 the Hanford 2km IFO had two calibration epochs, so 
    if the start time is during S2, we use the correct cache file.
    """

    self.calibration_cache_path()
    self.add_var_opt('calibration-cache', self.__calibration_cache_path)
    self.__calibration = self.__calibration_cache_path
    self.add_input_file(self.__calibration)

  def get_calibration(self):
    """
    Return the calibration cache file to be used by the
    DAG.
    """

    return self.__calibration_cache_path



class AnalysisChunk:
  """
  An AnalysisChunk is the unit of data that a node works with, usually some
  subset of a ScienceSegment.
  """
  def __init__(self, start, end, trig_start = 0, trig_end = 0):
    """
    @param start: GPS start time of the chunk.
    @param end: GPS end time of the chunk.
    @param trig_start: GPS time at which to start generating triggers
    @param trig_end: GPS time at which to stop generating triggers
    """
    self.__start = start
    self.__end = end
    self.__length = end - start
    self.__trig_start = trig_start
    self.__trig_end = trig_end

  def __repr__(self):
    if self.__trig_start and self.__trig_end:
      return '<AnalysisChunk: start %d, end %d, trig_start %d, trig_end %d>' % (
        self.__start, self.__end, self.__trig_start, self.__trig_end)
    elif self.__trig_start and not self.__trig_end:
      return '<AnalysisChunk: start %d, end %d, trig_start %d>' % (
        self.__start, self.__end, self.__trig_start)
    elif not self.__trig_start and self.__trig_end:
      return '<AnalysisChunk: start %d, end %d, trig_end %d>' % (
        self.__start, self.__end, self.__trig_end)
    else:
      return '<AnalysisChunk: start %d, end %d>' % (self.__start, self.__end)

  def __len__(self):
    """
    Returns the length of data for which this AnalysisChunk will produce
    triggers (in seconds).
    """
    if self.__trig_start and self.__trig_end:
      x = self.__trig_end - self.__trig_start
    elif self.__trig_start and not self.__trig_end:
      x = self.__end - self.__trig_start
    elif not self.__trig_start and self.__trig_end:
      x = self.__trig_end - self.__start
    else:
      x = self.__end - self.__start

    if x < 0:
      raise SegmentError, self + 'has negative length'
    else:
      return x
    
  def start(self):
    """
    Returns the GPS start time of the chunk.
    """
    return self.__start

  def end(self):
    """
    Returns the GPS end time of the chunk.
    """
    return self.__end
    
  def dur(self):
    """
    Returns the length (duration) of the chunk in seconds.
    """
    return self.__length

  def trig_start(self):
    """
    Return the first GPS time at which triggers for this chunk should be
    generated.
    """
    return self.__trig_start

  def trig_end(self):
    """
    Return the last GPS time at which triggers for this chunk should be
    generated.
    """
    return self.__trig_end

  def set_trig_start(self,start):
    """
    Set the first GPS time at which triggers for this chunk should be
    generated.
    """
    self.__trig_start = start

  def set_trig_end(self,end):
    """
    Set the last GPS time at which triggers for this chunk should be
    generated.
    """
    self.__trig_end = end



class ScienceSegment:
  """
  A ScienceSegment is a period of time where the experimenters determine
  that the inteferometer is in a state where the data is suitable for 
  scientific analysis. A science segment can have a list of AnalysisChunks
  asscociated with it that break the segment up into (possibly overlapping)
  smaller time intervals for analysis.
  """
  def __init__(self,segment):
    """
    @param segment: a tuple containing the (segment id, gps start time, gps end
    time, duration) of the segment.
    """
    self.__id = segment[0]
    self.__start = segment[1]
    self.__end = segment[2]
    self.__dur = segment[3]
    self.__chunks = []
    self.__unused = self.dur()
    self.__ifo = None

  def __getitem__(self,i):
    """
    Allows iteration over and direct access to the AnalysisChunks contained
    in this ScienceSegment.
    """
    if i < 0: raise IndexError, "list index out of range"
    return self.__chunks[i]
    
  def __len__(self):
    """
    Returns the number of AnalysisChunks contained in this ScienceSegment.
    """
    return len(self.__chunks)

  def __repr__(self):
    return '<ScienceSegment: id %d, start %d, end %d, dur %d, unused %d>' % (
    self.id(),self.start(),self.end(),self.dur(),self.__unused)

  def __cmp__(self,other):
    """
    ScienceSegments are compared by the GPS start time of the segment.
    """
    return cmp(self.start(),other.start())

  def make_chunks(self,length=0,overlap=0,play=0,sl=0,excl_play=0):
    """
    Divides the science segment into chunks of length seconds overlapped by
    overlap seconds. If the play option is set, only chunks that contain S2
    playground data are generated. If the user has a more complicated way
    of generating chunks, this method should be overriden in a sub-class.
    Any data at the end of the ScienceSegment that is too short to contain a 
    chunk is ignored. The length of this unused data is stored and can be
    retrieved with the unused() method.
    @param length: length of chunk in seconds.
    @param overlap: overlap between chunks in seconds.
    @param play: 1 : only generate chunks that overlap with S2 playground data.
                 2 : as play = 1 plus compute trig start and end times to coincide
                 with the start/end of the playground
    @param sl: slide by sl seconds before determining playground data.
    @param excl_play: exclude the first excl_play second from the start and end
    of the chunk when computing if the chunk overlaps with playground.
    """
    time_left = self.dur()
    start = self.start()
    increment = length - overlap
    while time_left >= length:
      end = start + length
      if (not play) or (play and (((end-sl-excl_play-729273613) % 6370) < 
        (600+length-2*excl_play))):
	if (play == 2):
	  # calculate the start of the playground preceeding the chunk end
	  play_start = 729273613 + 6370 * \
		math.floor((end-sl-excl_play-729273613) / 6370)
 	  play_end = play_start + 600
	  trig_start = 0
	  trig_end = 0
	  if ( (play_end - 6370) > start ):
	    print "Two playground segments in this chunk:",
	    print "  Code to handle this case has not been implemented"
	    sys.exit(1)
          else:
	    if play_start > start:
	      trig_start = int(play_start)
	    if play_end < end:
	      trig_end = int(play_end)
	  self.__chunks.append(AnalysisChunk(start,end,trig_start,trig_end))
        else:
          self.__chunks.append(AnalysisChunk(start,end))
      start += increment
      time_left -= increment
    self.__unused = time_left - overlap

  def add_chunk(self,start,end,trig_start=0,trig_end=0):
    """
    Add an AnalysisChunk to the list associated with this ScienceSegment.
    @param start: GPS start time of chunk.
    @param end: GPS end time of chunk.
    @param trig_start: GPS start time for triggers from chunk
    """
    self.__chunks.append(AnalysisChunk(start,end,trig_start,trig_end))

  def unused(self):
    """
    Returns the length of data in the science segment not used to make chunks.
    """
    return self.__unused

  def set_unused(self,unused):
    """
    Set the length of data in the science segment not used to make chunks.
    """
    self.__unused = unused

  def id(self):
    """
    Returns the ID of this ScienceSegment.
    """
    return self.__id
    
  def start(self):
    """
    Returns the GPS start time of this ScienceSegment.
    """
    return self.__start

  def end(self):
    """
    Returns the GPS end time of this ScienceSegment.
    """
    return self.__end

  def set_start(self,t):
    """
    Override the GPS start time (and set the duration) of this ScienceSegment.
    @param t: new GPS start time.
    """
    self.__dur += self.__start - t
    self.__start = t

  def set_end(self,t):
    """
    Override the GPS end time (and set the duration) of this ScienceSegment.
    @param t: new GPS end time.
    """
    self.__dur -= self.__end - t
    self.__end = t

  def dur(self):
    """
    Returns the length (duration) in seconds of this ScienceSegment.
    """
    return self.__dur




class ScienceData:
  """
  An object that can contain all the science data used in an analysis. Can
  contain multiple ScienceSegments and has a method to generate these from
  a text file produces by the LIGOtools segwizard program.
  """
  def __init__(self):
    self.__sci_segs = []
    self.__file = None

  def __getitem__(self,i):
    """
    Allows direct access to or iteration over the ScienceSegments associated
    with the ScienceData.
    """
    return self.__sci_segs[i]

  def __repr__(self):
    return '<ScienceData: file %s>' % self.__file

  def __len__(self):
    """
    Returns the number of ScienceSegments associated with the ScienceData.
    """
    return len(self.__sci_segs)

  def read(self,file,min_length,slide_sec=0,buffer=0):
    """
    Parse the science segments from the segwizard output contained in file.
    @param file: input text file containing a list of science segments generated by
    segwizard.
    @param min_length: only append science segments that are longer than min_length.
    @param slide_sec: Slide each ScienceSegment by::

      delta > 0:
        [s,e] -> [s+delta,e].
      delta < 0:
        [s,e] -> [s,e-delta].

    @param buffer: shrink the ScienceSegment::

      [s,e] -> [s+buffer,e-buffer]
    """
    self.__file = file
    octothorpe = re.compile(r'\A#')
    for line in open(file):
      if not octothorpe.match(line) and int(line.split()[3]) >= min_length:
        (id,st,en,du) = map(int,line.split())

        # slide the data if doing a background estimation
        if slide_sec > 0:
          st += slide_sec
        elif slide_sec < 0:
          en += slide_sec
        du -= abs(slide_sec)

        # add a buffer
        if buffer > 0:
          st += buffer
          en -= buffer
          du -= 2*abs(buffer)

        x = ScienceSegment(tuple([id,st,en,du]))
        self.__sci_segs.append(x)

  def tama_read(self,file):
    """
    Parse the science segments from a tama list of locked segments contained in
	 	file.
    @param file: input text file containing a list of tama segments.
    """
    self.__file = file
    for line in open(file):
      columns = line.split()
      id = int(columns[0])
      start = int(math.ceil(float(columns[3])))
      end = int(math.floor(float(columns[4])))
      dur = end - start	
    
      x = ScienceSegment(tuple([id, start, end, dur]))
      self.__sci_segs.append(x)


  def make_chunks(self,length,overlap=0,play=0,sl=0,excl_play=0):
    """
    Divide each ScienceSegment contained in this object into AnalysisChunks.
    @param length: length of chunk in seconds.
    @param overlap: overlap between segments.
    @param play: if true, only generate chunks that overlap with S2 playground data.
    @param sl: slide by sl seconds before determining playground data.
    @param excl_play: exclude the first excl_play second from the start and end
    of the chunk when computing if the chunk overlaps with playground.
    """
    for seg in self.__sci_segs:
      seg.make_chunks(length,overlap,play,sl,excl_play)

  def make_chunks_from_unused(self,length,trig_overlap,play=0,min_length=0,sl=0,excl_play=0):
    """
    Create an extra chunk that uses up the unused data in the science segment.
    @param length: length of chunk in seconds.
    @param trig_overlap: length of time start generating triggers before the
    start of the unused data.
    @param play: 
                - 1 : only generate chunks that overlap with S2 playground data.
                - 2 : as 1 plus compute trig start and end times to coincide
                        with the start/end of the playground
    @param min_length: the unused data must be greater than min_length to make a
    chunk.
    @param sl: slide by sl seconds before determining playground data.
    @param excl_play: exclude the first excl_play second from the start and end
    of the chunk when computing if the chunk overlaps with playground.
    """
    for seg in self.__sci_segs:
      # if there is unused data longer than the minimum chunk length
      if seg.unused() > min_length:
        start = seg.end() - length
        end = seg.end()
        middle = start + length / 2
        if (not play) or (play and (((end-sl-excl_play-729273613)%6370) < 
          (600+length-2*excl_play))):
          trig_start = end - seg.unused() - trig_overlap
	  if (play == 2):
            # calculate the start of the playground preceeding the chunk end
	    play_start = 729273613 + 6370 * \
              math.floor((end-sl-excl_play-729273613) / 6370)
 	    play_end = play_start + 600
	    trig_end = 0
	    if ( (play_end - 6370) > start ):
	      print "Two playground segments in this chunk"
	      print "  Code to handle this case has not been implemented"
	      sys.exit(1)
            else:
	      if play_start > trig_start:
	        trig_start = int(play_start)
	      if (play_end < end):
	        trig_end = int(play_end)
	      if (trig_end == 0) or (trig_end > trig_start):
                seg.add_chunk(start, end, trig_start, trig_end)
          else:
            seg.add_chunk(start, end, trig_start)
        seg.set_unused(0)

  def make_short_chunks_from_unused(
    self,min_length,overlap=0,play=0,sl=0,excl_play=0):
    """
    Create a chunk that uses up the unused data in the science segment
    @param min_length: the unused data must be greater than min_length to make a
    chunk.
    @param overlap: overlap between chunks in seconds.
    @param play: if true, only generate chunks that overlap with S2 playground data.
    @param sl: slide by sl seconds before determining playground data.
    @param excl_play: exclude the first excl_play second from the start and end
    of the chunk when computing if the chunk overlaps with playground.
    """
    for seg in self.__sci_segs:
      if seg.unused() > min_length:
        start = seg.end() - seg.unused() - overlap
        end = seg.end()
        length = start - end
        if (not play) or (play and (((end-sl-excl_play-729273613)%6370) < 
        (600+length-2*excl_play))):
          seg.add_chunk(start, end, start)
        seg.set_unused(0)

  def intersection(self, other):
    """
    Replaces the ScienceSegments contained in this instance of ScienceData
    with the intersection of those in the instance other. Returns the number
    of segments in the intersection.
    @param other: ScienceData to use to generate the intersection
    """

    # we only deal with the case of two lists here
    length1 = len(self)
    length2 = len(other)

    # initialize list of output segments
    ostart = -1
    outlist = []
    iseg2 = -1
    start2 = -1
    stop2 = -1

    for seg1 in self:
      start1 = seg1.start()
      stop1 = seg1.end()
      id = seg1.id()

      # loop over segments from the second list which overlap this segment
      while start2 < stop1:
        if stop2 > start1:
          # these overlap

          # find the overlapping range
          if start1 < start2:
            ostart = start2
          else:
            ostart = start1
          if stop1 > stop2:
            ostop = stop2
          else:
            ostop = stop1

          x = ScienceSegment(tuple([id, ostart, ostop, ostop-ostart]))
          outlist.append(x)

          if stop2 > stop1:
            break

        # step forward
        iseg2 += 1
        if iseg2 < len(other):
          seg2 = other[iseg2]
          start2 = seg2.start()
          stop2 = seg2.end()
        else:
          # pseudo-segment in the far future
          start2 = 2000000000
          stop2 = 2000000000

    # save the intersection and return the length
    self.__sci_segs = outlist
    return len(self)

  

  def union(self, other):
    """
    Replaces the ScienceSegments contained in this instance of ScienceData
    with the union of those in the instance other. Returns the number of
    ScienceSegments in the union.
    @param other: ScienceData to use to generate the intersection
    """

    # we only deal with the case of two lists here
    length1 = len(self)
    length2 = len(other)

    # initialize list of output segments
    ostart = -1
    seglist = []

    i1 = -1
    i2 = -1
    start1 = -1
    start2 = -1
    id = -1
    
    while 1:
      # if necessary, get a segment from list 1
      if start1 == -1:
        i1 += 1
        if i1 < length1:
          start1 = self[i1].start()
          stop1 = self[i1].end()
          id = self[i1].id()
        elif i2 == length2:
          break

      # if necessary, get a segment from list 2
      if start2 == -1:
        i2 += 1
        if i2 < length2:
          start2 = other[i2].start()
          stop2 = other[i2].end()
        elif i1 == length1:
          break

      # pick the earlier segment from the two lists
      if start1 > -1 and ( start2 == -1 or start1 <= start2):
        ustart = start1
        ustop = stop1
        # mark this segment has having been consumed
        start1 = -1
      elif start2 > -1:
        ustart = start2
        ustop = stop2
        # mark this segment has having been consumed
        start2 = -1
      else:
        break

      # if the output segment is blank, initialize it; otherwise, see
      # whether the new segment extends it or is disjoint
      if ostart == -1:
        ostart = ustart
        ostop = ustop
      elif ustart <= ostop:
        if ustop > ostop:
          # this extends the output segment
          ostop = ustop
        else:
          # This lies entirely within the current output segment
          pass
      else:
         # flush the current output segment, and replace it with the
         # new segment
         x = ScienceSegment(tuple([id,ostart,ostop,ostop-ostart]))
         seglist.append(x)
         ostart = ustart
         ostop = ustop

    # flush out the final output segment (if any)
    if ostart != -1:
      x = ScienceSegment(tuple([id,ostart,ostop,ostop-ostart]))
      seglist.append(x)

    self.__sci_segs = seglist
    return len(self)


  def coalesce(self):
    """
    Coalesces any adjacent ScienceSegments. Returns the number of 
    ScienceSegments in the coalesced list.
    """

    # check for an empty list
    if len(self) == 0:
      return 0

    # sort the list of science segments
    self.__sci_segs.sort()

    # coalesce the list, checking each segment for validity as we go
    outlist = []
    ostop = -1

    for seg in self:
      start = seg.start()
      stop = seg.end()
      id = seg.id()
      if start > ostop:
        # disconnected, so flush out the existing segment (if any)
        if ostop >= 0:
          x = ScienceSegment(tuple([id,ostart,ostop,ostop-ostart]))
          outlist.append(x)
        ostart = start
        ostop = stop
      elif stop > ostop:
        # extend the current segment
        ostop = stop

    # flush out the final segment (if any)
    if ostop >= 0:
      x = ScienceSegment(tuple([id,ostart,ostop,ostop-ostart]))
      outlist.append(x)

    self.__sci_segs = outlist
    return len(self)


  def invert(self):
    """
    Inverts the ScienceSegments in the class (i.e. set NOT).  Returns the
    number of ScienceSegments after inversion.
    """

    # check for an empty list
    if len(self) == 0:
      # return a segment representing all time
      self.__sci_segs = ScienceSegment(tuple(0,0,1999999999,1999999999))

    # go through the list checking for validity as we go
    outlist = []
    ostart = 0
    for seg in self:
      start = seg.start()
      stop = seg.end()
      if start < 0 or stop < start or start < ostart:
        raise SegmentError, "Invalid list"
      if start > 0:
        x = ScienceSegment(tuple([0,ostart,start,start-ostart]))
        outlist.append(x)
      ostart = stop

    if ostart < 1999999999:
      x = ScienceSegment(tuple([0,ostart,1999999999,1999999999-ostart]))
      outlist.append(x)

    self.__sci_segs = outlist
    return len(self)

  
  def play(self):
    """
    Keep only times in ScienceSegments which are in the playground
    """

    length = len(self)

    # initialize list of output segments
    ostart = -1
    outlist = []
    begin_s2 = 729273613
    play_space = 6370
    play_len = 600

    for seg in self:
      start = seg.start()
      stop = seg.end()
      id = seg.id()
     
      # select first playground segment which ends after start of seg
      play_start = begin_s2+play_space*( 1 + 
	int((start - begin_s2 - play_len)/play_space) )

      while play_start < stop:
	if play_start > start:
	  ostart = play_start
	else:
	  ostart = start
        
        
	play_stop = play_start + play_len

	if play_stop < stop:
	  ostop = play_stop
	else:
	  ostop = stop

        x = ScienceSegment(tuple([id, ostart, ostop, ostop-ostart]))
        outlist.append(x)

        # step forward
	play_start = play_start + play_space 

    # save the playground segs and return the length
    self.__sci_segs = outlist
    return len(self)


  def intersect_3(self, second, third):
    """
    Intersection routine for three inputs.  Built out of the intersect, 
    coalesce and play routines
    """
    self.intersection(second)
    self.intersection(third)
    self.coalesce()
    return len(self)


  
class LSCDataFindJob(CondorDAGJob, AnalysisJob):
  """
  An LSCdataFind job used to locate data. The static options are
  read from the section [datafind] in the ini file. The stdout from
  LSCdataFind contains the paths to the frame files and is directed to a file
  in the cache directory named by site and GPS start and end times. The stderr
  is directed to the logs directory. The job always runs in the scheduler
  universe. The path to the executable is determined from the ini file.
  """
  def __init__(self,cache_dir,log_dir,config_file):
    """
    @param cache_dir: the directory to write the output lal cache files to.
    @param log_dir: the directory to write the stderr file to.
    @param config_file: ConfigParser object containing the path to the LSCdataFind
    executable in the [condor] section and a [datafind] section from which
    the LSCdataFind options are read.
    """
    self.__executable = config_file.get('condor','datafind')
    self.__universe = 'scheduler'
    CondorDAGJob.__init__(self,self.__universe,self.__executable)
    AnalysisJob.__init__(self,config_file)
    self.__cache_dir = cache_dir

    for sec in ['datafind']:
      self.add_ini_opts(config_file,sec)
    
    # we need a lal cache for files on the localhost
    self.add_opt('match','localhost')
    self.add_opt('lal-cache','')
    self.add_opt('url-type','file')

    self.add_condor_cmd('environment',
      """LD_LIBRARY_PATH=$ENV(LD_LIBRARY_PATH);PYTHONPATH=$ENV(PYTHONPATH);LSC_DATAFIND_SERVER=$ENV(LSC_DATAFIND_SERVER)""" )

    self.set_stderr_file(log_dir + '/datafind-$(macroobservatory)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_stdout_file(self.__cache_dir + '/$(macroobservatory)-$(macrogpsstarttime)-$(macrogpsendtime).cache')
    self.set_sub_file('datafind.sub')

  def get_cache_dir(self):
    """
    returns the directroy that the cache files are written to.
    """
    return self.__cache_dir


class LSCDataFindNode(CondorDAGNode, AnalysisNode):
  """
  A DataFindNode runs an instance of LSCdataFind in a Condor DAG.
  """
  def __init__(self,job):
    """
    @param job: A CondorDAGJob that can run an instance of LALdataFind.
    """
    CondorDAGNode.__init__(self,job)
    AnalysisNode.__init__(self)
    self.__start = 0
    self.__end = 0
    self.__observatory = None
    self.__output = None
    self.__job = job
   
  def __set_output(self):
    """
    Private method to set the file to write the cache to. Automaticaly set
    once the ifo, start and end times have been set.
    """
    if self.__start and self.__end and self.__observatory:
      self.__output = self.__job.get_cache_dir() + '/' + self.__observatory + '-'  
      self.__output += str(self.__start) + '-' + str(self.__end) + '.cache'
      self.add_output_file(self.__output)

  def set_start(self,time):
    """
    Set the start time of the datafind query.
    @param time: GPS start time of query.
    """
    self.add_var_opt('gps-start-time', time)
    self.__start = time
    self.__set_output()

  def set_end(self,time):
    """
    Set the end time of the datafind query.
    @param time: GPS end time of query.
    """
    self.add_var_opt('gps-end-time', time)
    self.__end = time
    self.__set_output()

  def set_observatory(self,obs):
    """
    Set the IFO to retrieve data for. Since the data from both Hanford 
    interferometers is stored in the same frame file, this takes the first 
    letter of the IFO (e.g. L or H) and passes it to the --observatory option
    of LSCdataFind.
    @param obs: IFO to obtain data for.
    """
    self.add_var_opt('observatory',obs)
    self.__observatory = obs
    self.__set_output()

  def get_output(self):
    """
    Return the output file, i.e. the file containing the frame cache data.
    """
    return self.__output


