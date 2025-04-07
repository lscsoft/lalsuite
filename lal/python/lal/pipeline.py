"""
This modules contains objects that make it simple for the user to
create python scripts that build Condor DAGs to run code on the LSC
Data Grid.

This file is part of the Grid LSC User Environment (GLUE)

GLUE is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from . import git_version
__author__ = 'Duncan Brown <duncan@gravity.phys.uwm.edu>'
__version__ = "git id %s" % git_version.id
__date__ = git_version.date

from collections import OrderedDict
import configparser
import igwn_segments as segments
import itertools
import math
import os
import random
import re
import stat
import sys
import time
from hashlib import md5


class CondorError(Exception):
  """Error thrown by Condor Jobs"""
  def __init__(self, args=None):
    self.args = args


class CondorJobError(CondorError):
  pass


class CondorSubmitError(CondorError):
  pass


class CondorDAGError(CondorError):
  pass


class CondorDAGJobError(CondorError):
  pass


class CondorDAGNodeError(CondorError):
  pass


class SegmentError(Exception):
  def __init__(self, args=None):
    self.args = args


class CondorJob(object):
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
    self.__condor_cmds = OrderedDict()
    self.__notification = None
    self.__log_file = None
    self.__in_file = None
    self.__err_file = None
    self.__out_file = None
    self.__sub_file_path = None
    self.__output_files = []
    self.__input_files = []
    self.__checkpoint_files = []
    self.__grid_type = None
    self.__grid_server = None
    self.__grid_scheduler = None
    self.__executable_installed = True

  def get_executable(self):
    """
    Return the name of the executable for this job.
    """
    return self.__executable

  def set_executable(self, executable):
    """
    Set the name of the executable for this job.
    """
    self.__executable = executable

  def get_universe(self):
    """
    Return the condor universe that the job will run in.
    """
    return self.__universe

  def set_universe(self, universe):
    """
    Set the condor universe for the job to run in.
    @param universe: the condor universe to run the job in.
    """
    self.__universe = universe

  def get_grid_type(self):
    """
    Return the grid type of the job.
    """
    return self.__grid_type

  def set_grid_type(self, grid_type):
    """
    Set the type of grid resource for the job.
    @param grid_type: type of grid resource.
    """
    self.__grid_type = grid_type

  def get_grid_server(self):
    """
    Return the grid server on which the job will run.
    """
    return self.__grid_server

  def set_grid_server(self, grid_server):
    """
    Set the grid server on which to run the job.
    @param grid_server: grid server on which to run.
    """
    self.__grid_server = grid_server

  def get_grid_scheduler(self):
    """
    Return the grid scheduler.
    """
    return self.__grid_scheduler

  def set_grid_scheduler(self, grid_scheduler):
    """
    Set the grid scheduler.
    @param grid_scheduler: grid scheduler on which to run.
    """
    self.__grid_scheduler = grid_scheduler

  def set_executable_installed(self,installed):
    """
    If executable installed is true, then no copying of the executable is
    done. If it is false, pegasus stages the executable to the remote site.
    Default is executable is installed (i.e. True).
    @param installed: true or fale
    """
    self.__executable_installed = installed

  def get_executable_installed(self):
    """
    return whether or not the executable is installed
    """
    return self.__executable_installed

  def add_condor_cmd(self, cmd, value):
    """
    Add a Condor command to the submit file (e.g. a class add or evironment).
    @param cmd: Condor command directive.
    @param value: value for command.
    """
    self.__condor_cmds[cmd] = value

  def get_condor_cmds(self):
    """
    Return the dictionary of condor keywords to add to the job
    """
    return self.__condor_cmds

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

  def add_checkpoint_file(self, filename):
    """
    Add filename as a checkpoint file for this DAG job.
    """
    if filename not in self.__checkpoint_files:
        self.__checkpoint_files.append(filename)

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

  def get_checkpoint_files(self):
    """
    Return a list of checkpoint files for this DAG node
    """
    return self.__checkpoint_files

  def add_arg(self, arg):
    """
    Add an argument to the executable. Arguments are appended after any
    options and their order is guaranteed.
    @param arg: argument to add.
    """
    self.__arguments.append(arg)

  def add_file_arg(self, filename):
    """
    Add a file argument to the executable. Arguments are appended after any
    options and their order is guaranteed. Also adds the file name to the
    list of required input data for this job.
    @param filename: file to add as argument.
    """
    self.__arguments.append(filename)
    if filename not in self.__input_files:
      self.__input_files.append(filename)

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

  def get_opt( self, opt):
    """
    Returns the value associated with the given command line option.
    Returns None if the option does not exist in the options list.
    @param opt: command line option
    """
    if opt in self.__options:
      return self.__options[opt]
    return None

  def add_file_opt(self, opt, filename):
    """
    Add a command line option to the executable. The order that the arguments
    will be appended to the command line is not guaranteed, but they will
    always be added before any command line arguments. The name of the option
    is prefixed with double hyphen and the program is expected to parse it
    with getopt_long().
    @param opt: command line option to add.
    @param filename: value to pass to the option (None for no argument).
    """
    self.__options[opt] = filename
    if filename not in self.__input_files:
      self.__input_files.append(filename)

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
      arg = str(cp.get(section,opt)).strip()
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

  def set_stdin_file(self, path):
    """
    Set the file from which Condor directs the stdin of the job.
    @param path: path to stdin file.
    """
    self.__in_file = path

  def get_stdin_file(self):
    """
    Get the file from which Condor directs the stdin of the job.
    """
    return self.__in_file

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
      raise CondorSubmitError("Log file not specified.")
    if not self.__err_file:
      raise CondorSubmitError("Error file not specified.")
    if not self.__out_file:
      raise CondorSubmitError("Output file not specified.")

    if not self.__sub_file_path:
      raise CondorSubmitError('No path for submit file.')
    try:
      subfile = open(self.__sub_file_path, 'w')
    except:
      raise CondorSubmitError("Cannot open file " + self.__sub_file_path)

    if self.__universe == 'grid':
      if self.__grid_type is None:
        raise CondorSubmitError('No grid type specified.')
      elif self.__grid_type == 'gt2':
        if self.__grid_server is None:
          raise CondorSubmitError('No server specified for grid resource.')
      elif self.__grid_type == 'gt4':
        if self.__grid_server is None:
          raise CondorSubmitError('No server specified for grid resource.')
        if self.__grid_scheduler is None:
          raise CondorSubmitError('No scheduler specified for grid resource.')
      else:
        raise CondorSubmitError('Unsupported grid resource.')

    subfile.write( 'universe = ' + self.__universe + '\n' )
    subfile.write( 'executable = ' + self.__executable + '\n' )

    if self.__universe == 'grid':
      if self.__grid_type == 'gt2':
        subfile.write('grid_resource = %s %s\n' % (self.__grid_type,
          self.__grid_server))
      if self.__grid_type == 'gt4':
        subfile.write('grid_resource = %s %s %s\n' % (self.__grid_type,
          self.__grid_server, self.__grid_scheduler))

    if self.__universe == 'grid':
      subfile.write('when_to_transfer_output = ON_EXIT\n')
      subfile.write('transfer_output_files = $(macrooutput)\n')
      subfile.write('transfer_input_files = $(macroinput)\n')

    if list(self.__options.keys()) or list(self.__short_options.keys()) or self.__arguments:
      subfile.write( 'arguments = "' )
      for c in self.__arguments:
        subfile.write( ' ' + c )
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
      subfile.write( ' "\n' )

    for cmd in self.__condor_cmds.keys():
      subfile.write( str(cmd) + " = " + str(self.__condor_cmds[cmd]) + '\n' )

    subfile.write( 'log = ' + self.__log_file + '\n' )
    if self.__in_file is not None:
      subfile.write( 'input = ' + self.__in_file + '\n' )
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
    super(CondorDAGJob,self).__init__(universe, executable, 1)
    CondorJob.set_notification(self, 'never')
    self.__var_opts = []
    self.__arg_index = 0
    self.__var_args = []
    self.__var_cmds = []
    self.__grid_site = None
    self.__bad_macro_chars = re.compile(r'[_-]')

  def create_node(self):
    """
    Create a condor node from this job. This provides a basic interface to
    the CondorDAGNode class. Most jobs in a workflow will subclass the
    CondorDAGNode class and overwrite this to give more details when
    initializing the node. However, this will work fine for jobs with very simp
    input/output.
    """
    return CondorDAGNode(self)

  def set_grid_site(self,site):
    """
    Set the grid site to run on. If not specified,
    will not give hint to Pegasus
    """
    self.__grid_site=str(site)
    if site != 'local':
      self.set_executable_installed(False)

  def get_grid_site(self):
    """
    Return the grid site for this node
    """
    return self.__grid_site

  def add_var_opt(self, opt, short=False):
    """
    Add a variable (or macro) option to the condor job. The option is added
    to the submit file and a different argument to the option can be set for
    each node in the DAG.
    @param opt: name of option to add.
    @param short: opt is short
    """
    if opt not in self.__var_opts:
      self.__var_opts.append(opt)
      macro = self.__bad_macro_chars.sub( r'', opt )
      if short:
        self.add_short_opt(opt,'$(macro' + macro + ')')
      else:
        self.add_opt(opt,'$(macro' + macro + ')')

  def add_var_condor_cmd(self, command):
    """
    Add a condor command to the submit file that allows variable (macro)
    arguments to be passes to the executable.
    """
    if command not in self.__var_cmds:
        self.__var_cmds.append(command)
        macro = self.__bad_macro_chars.sub( r'', command )
        self.add_condor_cmd(command, '$(macro' + macro + ')')

  def add_var_arg(self,arg_index,quote=False):
    """
    Add a command to the submit file to allow variable (macro) arguments
    to be passed to the executable.
    """
    try:
      self.__var_args[arg_index]
    except IndexError:
      if arg_index != self.__arg_index:
        raise CondorDAGJobError("mismatch between job and node var_arg index")
      if quote:
          self.__var_args.append("'$(macroargument%s)'" % str(arg_index))
      else:
          self.__var_args.append('$(macroargument%s)' % str(arg_index))
      self.add_arg(self.__var_args[self.__arg_index])
      self.__arg_index += 1


class CondorDAGManJob(object):
  """
  Condor DAGMan job class. Appropriate for setting up DAGs to run within a
  DAG.
  """
  def __init__(self, dag, dir=None):
    """
    dag = the name of the condor dag file to run
    dir = the diretory in which the dag file is located
    """
    self.__dag = dag
    self.__notification = None
    self.__dag_directory = dir

  def create_node(self):
    """
    Create a condor node from this job. This provides a basic interface to
    the CondorDAGManNode class. Most jobs in a workflow will subclass the
    CondorDAGManNode class and overwrite this to give more details when
    initializing the node. However, this will work fine for jobs with very simp
    input/output.
    """
    return CondorDAGManNode(self)

  def set_dag_directory(self, dir):
    """
    Set the directory where the dag will be run
    @param dir: the name of the directory where the dag will be run
    """
    self.__dag_directory = dir

  def get_dag_directory(self):
    """
    Get the directory where the dag will be run
    """
    return self.__dag_directory

  def set_notification(self, value):
    """
    Set the email address to send notification to.
    @param value: email address or never for no notification.
    """
    self.__notification = value

  def get_sub_file(self):
    """
    Return the name of the dag as the submit file name for the
    SUBDAG EXTERNAL command in the uber-dag
    """
    return self.__dag

  def write_sub_file(self):
    """
    Do nothing as there is not need for a sub file with the
    SUBDAG EXTERNAL command in the uber-dag
    """
    pass

  def get_dag(self):
    """
    Return the name of any associated dag file
    """
    return self.__dag


class CondorDAGNode(object):
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
    if not isinstance(job, CondorDAGJob) and \
        not isinstance(job,CondorDAGManJob):
      raise CondorDAGNodeError(
          "A DAG node must correspond to a Condor DAG job or Condor DAGMan job")
    self.__name = None
    self.__job = job
    self.__category = None
    self.__priority = None
    self.__pre_script = None
    self.__pre_script_args = []
    self.__post_script = None
    self.__post_script_args = []
    self.__macros = {}
    self.__opts = {}
    self.__args = []
    self.__arg_index = 0
    self.__retry = 0
    self.__parents = []
    self.__bad_macro_chars = re.compile(r'[_-]')
    self.__output_files = []
    self.__input_files = []
    self.__checkpoint_files = []
    self.__vds_group = None
    if isinstance(job, CondorDAGJob) and job.get_universe() == 'standard':
      self.__grid_start = 'none'
    else:
      self.__grid_start = None

    # generate the md5 node name
    t = str( int( time.time() * 1000 ) )
    r = str( int( random.random() * 100000000000000000 ) )
    a = str( self.__class__ )
    self.__name = md5((t + r + a).encode()).hexdigest()
    self.__md5name = self.__name

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

  def get_post_script(self):
    """
    returns the name of the post script that is executed before the DAG node is
    run.
    """
    return self.__post_script

  def add_post_script_arg(self,arg):
    """
    Adds an argument to the post script that is executed before the DAG node is
    run.
    """
    self.__post_script_args.append(arg)

  def get_post_script_arg(self):
    """
    Returns and array of arguments to the post script that is executed before
    the DAG node is run.
    """
    return self.__post_script_args

  def set_name(self,name):
    """
    Set the name for this node in the DAG.
    """
    self.__name = str(name)

  def get_name(self):
    """
    Get the name for this node in the DAG.
    """
    return self.__name

  def set_category(self,category):
    """
    Set the category for this node in the DAG.
    """
    self.__category = str(category)

  def get_category(self):
    """
    Get the category for this node in the DAG.
    """
    return self.__category

  def set_priority(self,priority):
    """
    Set the priority for this node in the DAG.
    """
    self.__priority = str(priority)

  def get_priority(self):
    """
    Get the priority for this node in the DAG.
    """
    return self.__priority

  def add_input_file(self, filename):
    """
    Add filename as a necessary input file for this DAG node.

    @param filename: input filename to add
    """
    if filename not in self.__input_files:
      self.__input_files.append(filename)
      if not isinstance(self.job(), CondorDAGManJob):
        if self.job().get_universe() == 'grid':
          self.add_input_macro(filename)

  def add_output_file(self, filename):
    """
    Add filename as a output file for this DAG node.

    @param filename: output filename to add
    """
    if filename not in self.__output_files:
      self.__output_files.append(filename)
      if not isinstance(self.job(), CondorDAGManJob):
        if self.job().get_universe() == 'grid':
          self.add_output_macro(filename)

  def add_checkpoint_file(self,filename):
    """
    Add filename as a checkpoint file for this DAG node
    @param filename: checkpoint filename to add
    """
    if filename not in self.__checkpoint_files:
        self.__checkpoint_files.append(filename)
        if not isinstance(self.job(), CondorDAGManJob):
            if self.job().get_universe() == 'grid':
                self.add_checkpoint_macro(filename)

  def get_input_files(self):
    """
    Return list of input files for this DAG node and its job.
    """
    input_files = list(self.__input_files)
    if isinstance(self.job(), CondorDAGJob):
      input_files = input_files + self.job().get_input_files()
    return input_files

  def get_output_files(self):
    """
    Return list of output files for this DAG node and its job.
    """
    output_files = list(self.__output_files)
    if isinstance(self.job(), CondorDAGJob):
      output_files = output_files + self.job().get_output_files()
    return output_files

  def get_checkpoint_files(self):
    """
    Return a list of checkpoint files for this DAG node and its job.
    """
    checkpoint_files = list(self.__checkpoint_files)
    if isinstance(self.job(), CondorDAGJob):
        checkpoint_files = checkpoint_files + self.job().get_checkpoint_files()
    return checkpoint_files

  def set_vds_group(self,group):
    """
    Set the name of the VDS group key when generating a DAX
    @param group: name of group for thus nore
    """
    self.__vds_group = str(group)

  def get_vds_group(self):
    """
    Returns the VDS group key for this node
    """
    return self.__vds_group

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

  def add_io_macro(self,io,filename):
    """
    Add a variable (macro) for storing the input/output files associated
    with this node.
    @param io: macroinput or macrooutput
    @param filename: filename of input/output file
    """
    io = self.__bad_macro_chars.sub( r'', io )
    if io not in self.__opts:
      self.__opts[io] = filename
    else:
      if filename not in self.__opts[io]:
        self.__opts[io] += ',%s' % filename

  def add_input_macro(self,filename):
    """
    Add a variable (macro) for storing the input files associated with
    this node.
    @param filename: filename of input file
    """
    self.add_io_macro('macroinput', filename)

  def add_output_macro(self,filename):
    """
    Add a variable (macro) for storing the output files associated with
    this node.
    @param filename: filename of output file
    """
    self.add_io_macro('macrooutput', filename)

  def add_checkpoint_macro(self,filename):
    self.add_io_macro('macrocheckpoint',filename)

  def get_opts(self):
    """
    Return the opts for this node. Note that this returns only
    the options for this instance of the node and not those
    associated with the underlying job template.
    """
    return self.__opts

  def add_var_condor_cmd(self, command, value):
    """
    Add a variable (macro) condor command for this node. If the command
    specified does not exist in the CondorJob, it is added so the submit file
    will be correct.
    PLEASE NOTE: AS with other add_var commands, the variable must be set for
    all nodes that use the CondorJob instance.
    @param command: command name
    @param value: Value of the command for this node in the DAG/DAX.
    """
    macro = self.__bad_macro_chars.sub( r'', command )
    self.__macros['macro' + macro] = value
    self.__job.add_var_condor_cmd(command)

  def add_var_opt(self,opt,value,short=False):
    """
    Add a variable (macro) option for this node. If the option
    specified does not exist in the CondorJob, it is added so the submit
    file will be correct when written.
    @param opt: option name.
    @param value: value of the option for this node in the DAG.
    @param short: opt is short
    """
    macro = self.__bad_macro_chars.sub( r'', opt )
    self.__opts['macro' + macro] = value
    self.__job.add_var_opt(opt,short)

  def add_file_opt(self,opt,filename,file_is_output_file=False):
    """
    Add a variable (macro) option for this node. If the option
    specified does not exist in the CondorJob, it is added so the submit
    file will be correct when written. The value of the option is also
    added to the list of input files for the DAX.
    @param opt: option name.
    @param filename: value of the option for this node in the DAG.
    @param file_is_output_file: A boolean if the file will be an output file
    instead of an input file.  The default is to have it be an input.
    """
    self.add_var_opt(opt,filename)
    if file_is_output_file: self.add_output_file(filename)
    else: self.add_input_file(filename)

  def add_var_arg(self, arg,quote=False):
    """
    Add a variable (or macro) argument to the condor job. The argument is
    added to the submit file and a different value of the argument can be set
    for each node in the DAG.
    @param arg: name of option to add.
    @param quote: quote
    """
    self.__args.append(arg)
    self.__job.add_var_arg(self.__arg_index,quote=quote)
    self.__arg_index += 1

  def add_file_arg(self, filename):
    """
    Add a variable (or macro) file name argument to the condor job. The
    argument is added to the submit file and a different value of the
    argument can be set for each node in the DAG. The file name is also
    added to the list of input files for the DAX.
    @param filename: name of option to add.
    """
    self.add_input_file(filename)
    self.add_var_arg(filename)

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

  def get_retry(self):
    """
    Return the number of times that this node in the DAG should retry.
    """
    return self.__retry

  def write_job(self,fh):
    """
    Write the DAG entry for this node's job to the DAG file descriptor.
    @param fh: descriptor of open DAG file.
    """
    if isinstance(self.job(),CondorDAGManJob):
      # create an external subdag from this dag
      fh.write( ' '.join(
        ['SUBDAG EXTERNAL', self.__name, self.__job.get_sub_file()]) )
      if self.job().get_dag_directory():
        fh.write( ' DIR ' + self.job().get_dag_directory() )
    else:
      # write a regular condor job
      fh.write( 'JOB ' + self.__name + ' ' + self.__job.get_sub_file() )
    fh.write( '\n')

    fh.write( 'RETRY ' + self.__name + ' ' + str(self.__retry) + '\n' )

  def write_category(self,fh):
    """
    Write the DAG entry for this node's category to the DAG file descriptor.
    @param fh: descriptor of open DAG file.
    """
    fh.write('CATEGORY ' + self.__name + ' ' + self.__category + '\n')

  def write_priority(self,fh):
    """
    Write the DAG entry for this node's priority to the DAG file descriptor.
    @param fh: descriptor of open DAG file.
    """
    fh.write('PRIORITY ' + self.__name + ' ' + self.__priority + '\n')

  def write_vars(self,fh):
    """
    Write the variable (macro) options and arguments to the DAG file
    descriptor.
    @param fh: descriptor of open DAG file.
    """
    if list(self.__macros.keys()) or list(self.__opts.keys()) or self.__args:
      fh.write( 'VARS ' + self.__name )
    for k in self.__macros.keys():
      fh.write( ' ' + str(k) + '="' + str(self.__macros[k]) + '"' )
    for k in self.__opts.keys():
      fh.write( ' ' + str(k) + '="' + str(self.__opts[k]) + '"' )
    if self.__args:
      for i in range(self.__arg_index):
        fh.write( ' macroargument' + str(i) + '="' + self.__args[i] + '"' )
    fh.write( '\n' )

  def write_parents(self,fh):
    """
    Write the parent/child relations for this job to the DAG file descriptor.
    @param fh: descriptor of open DAG file.
    """
    if len(self.__parents) > 0:
      fh.write( 'PARENT ' + " ".join((str(p) for p in self.__parents)) + ' CHILD ' + str(self) + '\n' )

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
        fh.write("## Job %s requires input file %s\n" % (self.__name, f))

  def write_output_files(self, fh):
    """
    Write as a comment into the DAG file the list of output files
    for this DAG node.

    @param fh: descriptor of open DAG file.
    """
    for f in self.__output_files:
        fh.write("## Job %s generates output file %s\n" % (self.__name, f))

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
    if not isinstance(node, (CondorDAGNode,CondorDAGManNode) ):
      raise CondorDAGNodeError("Parent must be a CondorDAGNode or a CondorDAGManNode")
    self.__parents.append( node )

  def get_cmd_tuple_list(self):
    """
    Return a list of tuples containg the command line arguments
    """

    # pattern to find DAGman macros
    pat = re.compile(r'\$\((.+)\)')
    argpat = re.compile(r'\d+')

    # first parse the arguments and replace macros with values
    args = self.job().get_args()
    macros = self.get_args()

    cmd_list = []

    for a in args:
      m = pat.search(a)
      if m:
        arg_index = int(argpat.findall(a)[0])
        try:
          cmd_list.append(("%s" % macros[arg_index], ""))
        except IndexError:
          cmd_list.append("")
      else:
        cmd_list.append(("%s" % a, ""))

    # second parse the options and replace macros with values
    options = self.job().get_opts()
    macros = self.get_opts()

    for k in options:
      val = options[k]
      m = pat.match(val)
      if m:
        key = m.group(1)
        value = macros[key]

        cmd_list.append(("--%s" % k, str(value)))
      else:
        cmd_list.append(("--%s" % k, str(val)))

    # lastly parse the short options and replace macros with values
    options = self.job().get_short_opts()

    for k in options:
      val = options[k]
      m = pat.match(val)
      if m:
        key = m.group(1)
        value = macros[key]

        cmd_list.append(("-%s" % k, str(value)))
      else:
        cmd_list.append(("-%s" % k, str(val)))

    return cmd_list

  def get_cmd_line(self):
    """
    Return the full command line that will be used when this node
    is run by DAGman.
    """

    cmd = ""
    cmd_list = self.get_cmd_tuple_list()
    for argument in cmd_list:
      cmd += ' '.join(argument) + " "

    return cmd

  def finalize(self):
    """
    The finalize method of a node is called before the node is
    finally added to the DAG and can be overridden to do any last
    minute clean up (such as setting extra command line arguments)
    """
    pass


class CondorDAGManNode(CondorDAGNode):
  """
  Condor DAGMan node class. Appropriate for setting up DAGs to run within a
  DAG. Adds the user-tag functionality to condor_dagman processes running in
  the DAG. May also be used to extend dagman-node specific functionality.
  """
  def __init__(self, job):
    """
    @param job: a CondorDAGNodeJob
    """
    super(CondorDAGManNode,self).__init__(job)
    self.__user_tag = None
    self.__maxjobs_categories = []
    self.__cluster_jobs = None

  def set_user_tag(self,usertag):
    """
    Set the user tag that is passed to the analysis code.
    @param usertag: the user tag to identify the job
    """
    self.__user_tag = str(usertag)

  def get_user_tag(self):
    """
    Returns the usertag string
    """
    return self.__user_tag

  def add_maxjobs_category(self,categoryName,maxJobsNum):
    """
    Add a category to this DAG called categoryName with a maxjobs of maxJobsNum.
    @param categoryName: category name
    @param maxJobsNum: max jobs num
    """
    self.__maxjobs_categories.append((str(categoryName),str(maxJobsNum)))

  def get_maxjobs_categories(self):
    """
    Return an array of tuples containing (categoryName,maxJobsNum)
    """
    return self.__maxjobs_categories

  def set_cluster_jobs(self,cluster):
    """
    Set the type of job clustering pegasus can use to collapse jobs
    @param cluster: clustering type
    """
    self.__cluster_jobs = str(cluster)

  def get_cluster_jobs(self):
    """
    Returns the usertag string
    """
    return self.__cluster_jobs


class CondorDAG(object):
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
    self.__maxjobs_categories = []
    self.__integer_node_names = 0
    self.__node_count = 0
    self.__nodes_finalized = 0

  def get_nodes(self):
    """
    Return a list containing all the nodes in the DAG
    """
    return self.__nodes

  def get_jobs(self):
    """
    Return a list containing all the jobs in the DAG
    """
    return self.__jobs

  def set_integer_node_names(self):
    """
    Use integer node names for the DAG
    """
    self.__integer_node_names = 1

  def set_dag_file(self, path):
    """
    Set the name of the file into which the DAG is written.
    @param path: path to DAG file.
    """
    self.__dag_file_path = path + '.dag'

  def get_dag_file(self):
    """
    Return the path to the DAG file.
    """
    if not self.__log_file_path:
      raise CondorDAGError("No path for DAG file")
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
      raise CondorDAGError("Nodes must be class CondorDAGNode or subclass")
    if not isinstance(node.job(), CondorDAGManJob):
      node.set_log_file(self.__log_file_path)
    self.__nodes.append(node)
    if self.__integer_node_names:
      node.set_name(str(self.__node_count))
    self.__node_count += 1
    if node.job() not in self.__jobs:
      self.__jobs.append(node.job())

  def add_maxjobs_category(self,categoryName,maxJobsNum):
    """
    Add a category to this DAG called categoryName with a maxjobs of maxJobsNum.
    @param categoryName: category name
    @param maxJobsNum: max jobs num
    """
    self.__maxjobs_categories.append((str(categoryName),str(maxJobsNum)))

  def get_maxjobs_categories(self):
    """
    Return an array of tuples containing (categoryName,maxJobsNum)
    """
    return self.__maxjobs_categories

  def write_maxjobs(self,fh,category):
    """
    Write the DAG entry for this category's maxjobs to the DAG file descriptor.
    @param fh: descriptor of open DAG file.
    @param category: tuple containing type of jobs to set a maxjobs limit for
        and the maximum number of jobs of that type to run at once.
    """
    fh.write('MAXJOBS ' + str(category[0]) + ' ' + str(category[1]) + '\n')

  def write_sub_files(self):
    """
    Write all the submit files used by the dag to disk. Each submit file is
    written to the file name set in the CondorJob.
    """
    if not self.__nodes_finalized:
      for node in self.__nodes:
        node.finalize()
    for job in self.__jobs:
      job.write_sub_file()

  def write_concrete_dag(self):
    """
    Write all the nodes in the DAG to the DAG file.
    """
    if not self.__dag_file_path:
      raise CondorDAGError("No path for DAG file")
    try:
      dagfile = open( self.__dag_file_path, 'w' )
    except:
      raise CondorDAGError("Cannot open file " + self.__dag_file_path)
    for node in self.__nodes:
      node.write_job(dagfile)
      node.write_vars(dagfile)
      if node.get_category():
        node.write_category(dagfile)
      if node.get_priority():
        node.write_priority(dagfile)
      node.write_pre_script(dagfile)
      node.write_post_script(dagfile)
      node.write_input_files(dagfile)
      node.write_output_files(dagfile)
    for node in self.__nodes:
      node.write_parents(dagfile)
    for category in self.__maxjobs_categories:
      self.write_maxjobs(dagfile, category)
    dagfile.close()

  def write_dag(self):
    """
    Write a dag.
    """
    if not self.__nodes_finalized:
      for node in self.__nodes:
        node.finalize()
    self.write_concrete_dag()

  def write_script(self):
    """
    Write the workflow to a script (.sh instead of .dag).

    Assuming that parents were added to the DAG before their children,
    dependencies should be handled correctly.
    """
    if not self.__dag_file_path:
      raise CondorDAGError("No path for DAG file")
    try:
      dfp = self.__dag_file_path
      outfilename = ".".join(dfp.split(".")[:-1]) + ".sh"
      outfile = open(outfilename, "w")
    except:
      raise CondorDAGError("Cannot open file " + self.__dag_file_path)

    for node in self.__nodes:
        outfile.write("# Job %s\n" % node.get_name())
        # Check if this is a DAGMAN Node
        if isinstance(node,CondorDAGManNode):
          outfile.write("condor_submit_dag %s\n\n" % (node.job().get_dag()))
        else:
          outfile.write("%s %s\n\n" % (node.job().get_executable(),
              node.get_cmd_line()))
    outfile.close()

    os.chmod(outfilename, os.stat(outfilename)[0] | stat.S_IEXEC)


class AnalysisJob(object):
  """
  Describes a generic analysis job that filters LIGO data as configured by
  an ini file.
  """
  def __init__(self,cp):
    """
    @param cp: ConfigParser object that contains the configuration for this job.
    """
    self.__cp = cp
    try:
      self.__channel = str(self.__cp.get('input','channel')).strip()
    except:
      self.__channel = None

  def get_config(self,sec,opt):
    """
    Get the configration variable in a particular section of this jobs ini
    file.
    @param sec: ini file section.
    @param opt: option from section sec.
    """
    return str(self.__cp.get(sec,opt)).strip()

  def set_channel(self,channel):
    """
    Set the name of the channel that this job is filtering.  This will
    overwrite the value obtained at initialization.
    """
    self.__channel = channel

  def channel(self):
    """
    Returns the name of the channel that this job is filtering. Note that
    channel is defined to be IFO independent, so this may be LSC-AS_Q or
    IOO-MC_F. The IFO is set on a per node basis, not a per job basis.
    """
    return self.__channel


class AnalysisNode(object):
  """
  Contains the methods that allow an object to be built to analyse LIGO
  data in a Condor DAG.
  """
  def __init__(self):
    self.__start = 0
    self.__end = 0
    self.__data_start = 0
    self.__pad_data = 0
    self.__data_end = 0
    self.__trig_start = 0
    self.__trig_end = 0
    self.__ifo = None
    self.__ifo_tag = None
    self.__input = None
    self.__output = None
    self.__calibration = None
    self.__calibration_cache = None
    self.__LHO2k = re.compile(r'H2')
    self.__user_tag = self.job().get_opts().get("user-tag", None)

  def set_start(self,time,pass_to_command_line=True):
    """
    Set the GPS start time of the analysis node by setting a --gps-start-time
    option to the node when it is executed.
    @param time: GPS start time of job.
    @param pass_to_command_line: add gps-start-time as variable option.
    """
    if pass_to_command_line:
      self.add_var_opt('gps-start-time',time)
    self.__start = time
    self.__data_start = time
    #if not self.__calibration and self.__ifo and self.__start > 0:
    #  self.calibration()

  def get_start(self):
    """
    Get the GPS start time of the node.
    """
    return self.__start

  def set_end(self,time,pass_to_command_line=True):
    """
    Set the GPS end time of the analysis node by setting a --gps-end-time
    option to the node when it is executed.
    @param time: GPS end time of job.
    @param pass_to_command_line: add gps-end-time as variable option.
    """
    if pass_to_command_line:
      self.add_var_opt('gps-end-time',time)
    self.__end = time
    self.__data_end = time

  def get_end(self):
    """
    Get the GPS end time of the node.
    """
    return self.__end

  def set_data_start(self,time):
    """
    Set the GPS start time of the data needed by this analysis node.
    @param time: GPS start time of job.
    """
    self.__data_start = time

  def get_data_start(self):
    """
    Get the GPS start time of the data needed by this node.
    """
    return self.__data_start

  def set_pad_data(self,pad):
    """
    Set the GPS start time of the data needed by this analysis node.
    @param pad: pad
    """
    self.__pad_data = pad

  def get_pad_data(self):
    """
    Get the GPS start time of the data needed by this node.
    """
    return self.__pad_data

  def set_data_end(self,time):
    """
    Set the GPS end time of the data needed by this analysis node.
    @param time: GPS end time of job.
    """
    self.__data_end = time

  def get_data_end(self):
    """
    Get the GPS end time of the data needed by this node.
    """
    return self.__data_end

  def set_trig_start(self,time,pass_to_command_line=True):
    """
    Set the trig start time of the analysis node by setting a
    --trig-start-time option to the node when it is executed.
    @param time: trig start time of job.
    @param pass_to_command_line: add trig-start-time as a variable option.
    """
    if pass_to_command_line:
      self.add_var_opt('trig-start-time',time)
    self.__trig_start = time

  def get_trig_start(self):
    """
    Get the trig start time of the node.
    """
    return self.__trig_start

  def set_trig_end(self,time,pass_to_command_line=True):
    """
    Set the trig end time of the analysis node by setting a --trig-end-time
    option to the node when it is executed.
    @param time: trig end time of job.
    @param pass_to_command_line: add trig-end-time as a variable option.
    """
    if pass_to_command_line:
      self.add_var_opt('trig-end-time',time)
    self.__trig_end = time

  def get_trig_end(self):
    """
    Get the trig end time of the node.
    """
    return self.__trig_end

  def set_input(self,filename,pass_to_command_line=True):
    """
    Add an input to the node by adding a --input option.
    @param filename: option argument to pass as input.
    @param pass_to_command_line: add input as a variable option.
    """
    self.__input = filename
    if pass_to_command_line:
      self.add_var_opt('input', filename)
    self.add_input_file(filename)

  def get_input(self):
    """
    Get the file that will be passed as input.
    """
    return self.__input

  def set_output(self,filename,pass_to_command_line=True):
    """
    Add an output to the node by adding a --output option.
    @param filename: option argument to pass as output.
    @param pass_to_command_line: add output as a variable option.
    """
    self.__output = filename
    if pass_to_command_line:
      self.add_var_opt('output', filename)
    self.add_output_file(filename)

  def get_output(self):
    """
    Get the file that will be passed as output.
    """
    return self.__output

  def set_ifo(self,ifo):
    """
    Set the ifo name to analyze. If the channel name for the job is defined,
    then the name of the ifo is prepended to the channel name obtained
    from the job configuration file and passed with a --channel-name option.
    @param ifo: two letter ifo code (e.g. L1, H1 or H2).
    """
    self.__ifo = ifo
    if self.job().channel():
      self.add_var_opt('channel-name', ifo + ':' + self.job().channel())

  def get_ifo(self):
    """
    Returns the two letter IFO code for this node.
    """
    return self.__ifo

  def set_ifo_tag(self,ifo_tag,pass_to_command_line=True):
    """
    Set the ifo tag that is passed to the analysis code.
    @param ifo_tag: a string to identify one or more IFOs
    @param pass_to_command_line: add ifo-tag as a variable option.
    """
    self.__ifo_tag = ifo_tag
    if pass_to_command_line:
      self.add_var_opt('ifo-tag', ifo_tag)

  def get_ifo_tag(self):
    """
    Returns the IFO tag string
    """
    return self.__ifo_tag

  def set_user_tag(self,usertag,pass_to_command_line=True):
    """
    Set the user tag that is passed to the analysis code.
    @param usertag: the user tag to identify the job
    @param pass_to_command_line: add user-tag as a variable option.
    """
    self.__user_tag = usertag
    if pass_to_command_line:
      self.add_var_opt('user-tag', usertag)

  def get_user_tag(self):
    """
    Returns the usertag string
    """
    return self.__user_tag

  def set_cache(self,filename):
    """
    Set the LAL frame cache to to use. The frame cache is passed to the job
    with the --frame-cache argument.
    @param filename: calibration file to use.
    """
    if isinstance( filename, str ):
      # the name of a lal cache file created by a datafind node
      self.add_var_opt('frame-cache', filename)
      self.add_input_file(filename)
    elif isinstance( filename, list ):
      # we have an LFN list
      self.add_var_opt('glob-frame-data',' ')
      # only add the LFNs that actually overlap with this job
      # XXX FIXME this is a very slow algorithm
      if len(filename) == 0:
        raise CondorDAGNodeError(
          "LDR did not return any LFNs for query: check ifo and frame type")
      for lfn in filename:
        a, b, c, d = lfn.split('.')[0].split('-')
        t_start = int(c)
        t_end = int(c) + int(d)
        if (t_start <= (self.get_data_end() + self.get_pad_data() + int(d) + 1)
          and t_end >= (self.get_data_start() - self.get_pad_data() - int(d) - 1)):
          self.add_input_file(lfn)
      # set the frame type based on the LFNs returned by datafind
      self.add_var_opt('frame-type', b)
    else:
      raise CondorDAGNodeError("Unknown LFN cache format")

  def calibration_cache_path(self):
    """
    Determine the path to the correct calibration cache file to use.
    """
    if self.__ifo and self.__start > 0:
        cal_path = self.job().get_config('calibration', 'path')

        # check if this is S2: split calibration epochs
        if ( self.__LHO2k.match(self.__ifo) and
          (self.__start >= 729273613) and (self.__start <= 734367613) ):
          if self.__start < int(
            self.job().get_config('calibration','H2-cal-epoch-boundary')):
            cal_file = self.job().get_config('calibration','H2-1')
          else:
            cal_file = self.job().get_config('calibration','H2-2')
        else:
            # if not: just add calibration cache
            cal_file = self.job().get_config('calibration',self.__ifo)

        cal = os.path.join(cal_path,cal_file)
        self.__calibration_cache = cal
    else:
       msg = "IFO and start-time must be set first"
       raise CondorDAGNodeError(msg)

  def calibration(self):
    """
    Set the path to the calibration cache file for the given IFO.
    During S2 the Hanford 2km IFO had two calibration epochs, so
    if the start time is during S2, we use the correct cache file.
    """
    # figure out the name of the calibration cache files
    # as specified in the ini-file
    self.calibration_cache_path()

    # old .calibration for DAG's
    self.add_var_opt('calibration-cache', self.__calibration_cache)
    self.__calibration = self.__calibration_cache
    self.add_input_file(self.__calibration)

  def get_calibration(self):
    """
    Return the calibration cache file to be used by the
    DAG.
    """
    return self.__calibration_cache


class AnalysisChunk(object):
  """
  An AnalysisChunk is the unit of data that a node works with, usually some
  subset of a ScienceSegment.
  """
  def __init__(self, start, end, trig_start=0, trig_end=0):
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
      raise SegmentError(self + 'has negative length')
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


class ScienceSegment(object):
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
    self.__df_node = None

  def __getitem__(self,i):
    """
    Allows iteration over and direct access to the AnalysisChunks contained
    in this ScienceSegment.
    """
    if i < 0: raise IndexError("list index out of range")
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

  def make_chunks(self,length=0,overlap=0,play=0,sl=0,excl_play=0,pad_data=0):
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
                 2 : as play = 1 plus compute trig start and end times to
                     coincide with the start/end of the playground
    @param sl: slide by sl seconds before determining playground data.
    @param excl_play: exclude the first excl_play second from the start and end
    of the chunk when computing if the chunk overlaps with playground.
    @param pad_data: exclude the first and last pad_data seconds of the segment
    when generating chunks
    """
    time_left = self.dur() - (2 * pad_data)
    start = self.start() + pad_data
    increment = length - overlap
    while time_left >= length:
      end = start + length
      if (not play) or (play and (((end - sl - excl_play - 729273613) % 6370)
        < (600 + length - 2 * excl_play))):
        if (play == 2):
        # calculate the start of the playground preceeding the chunk end
          play_start = 729273613 + 6370 * \
           math.floor((end - sl - excl_play - 729273613) / 6370)
          play_end = play_start + 600
          trig_start = 0
          trig_end = 0
          if ( (play_end - 6370) > start ):
            print("Two playground segments in this chunk:", end=' ')
            print("  Code to handle this case has not been implemented")
            sys.exit(1)
          else:
            if play_start > start:
              trig_start = int(play_start)
            if play_end < end:
              trig_end = int(play_end)
          self.__chunks.append(AnalysisChunk(start, end, trig_start, trig_end))
        else:
          self.__chunks.append(AnalysisChunk(start, end))
      start += increment
      time_left -= increment
    self.__unused = time_left - overlap

  def add_chunk(self, start, end, trig_start=0, trig_end=0):
    """
    Add an AnalysisChunk to the list associated with this ScienceSegment.
    @param start: GPS start time of chunk.
    @param end: GPS end time of chunk.
    @param trig_start: GPS start time for triggers from chunk
    @param trig_end: trig_end
    """
    self.__chunks.append(AnalysisChunk(start, end, trig_start, trig_end))

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

  def set_df_node(self,df_node):
    """
    Set the DataFind node associated with this ScienceSegment to df_node.
    @param df_node: the DataFind node for this ScienceSegment.
    """
    self.__df_node = df_node

  def get_df_node(self):
    """
    Returns the DataFind node for this ScienceSegment.
    """
    return self.__df_node


class ScienceData(object):
  """
  An object that can contain all the science data used in an analysis. Can
  contain multiple ScienceSegments and has a method to generate these from
  a text file produces by the LIGOtools segwizard program.
  """
  def __init__(self):
    self.__sci_segs = []
    self.__filename = None

  def __getitem__(self,i):
    """
    Allows direct access to or iteration over the ScienceSegments associated
    with the ScienceData.
    """
    return self.__sci_segs[i]

  def __repr__(self):
    return '<ScienceData: file %s>' % self.__filename

  def __len__(self):
    """
    Returns the number of ScienceSegments associated with the ScienceData.
    """
    return len(self.__sci_segs)

  def read(self,filename,min_length,slide_sec=0,buffer=0):
    """
    Parse the science segments from the segwizard output contained in file.
    @param filename: input text file containing a list of science segments generated by
    segwizard.
    @param min_length: only append science segments that are longer than min_length.
    @param slide_sec: Slide each ScienceSegment by::

      delta > 0:
        [s, e] -> [s+delta, e].
      delta < 0:
        [s, e] -> [s, e-delta].

    @param buffer: shrink the ScienceSegment::

      [s, e] -> [s+buffer, e-buffer]
    """
    self.__filename = filename
    octothorpe = re.compile(r'\A#')
    for line in open(filename):
      if not octothorpe.match(line) and int(line.split()[3]) >= min_length:
        (id, st, en, du) = list(map(int, line.split()))

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
          du -= 2 * abs(buffer)

        x = ScienceSegment(tuple([id, st, en, du]))
        self.__sci_segs.append(x)

  def append_from_tuple(self, seg_tuple):
    x = ScienceSegment(seg_tuple)
    self.__sci_segs.append(x)

  def tama_read(self, filename):
    """
    Parse the science segments from a tama list of locked segments contained in
                file.
    @param filename: input text file containing a list of tama segments.
    """
    self.__filename = filename
    for line in open(filename):
      columns = line.split()
      id = int(columns[0])
      start = int(math.ceil(float(columns[3])))
      end = int(math.floor(float(columns[4])))
      dur = end - start

      x = ScienceSegment(tuple([id, start, end, dur]))
      self.__sci_segs.append(x)

  def make_chunks(self, length, overlap=0, play=0, sl=0, excl_play=0, pad_data=0):
    """
    Divide each ScienceSegment contained in this object into AnalysisChunks.
    @param length: length of chunk in seconds.
    @param overlap: overlap between segments.
    @param play: if true, only generate chunks that overlap with S2 playground
    data.
    @param sl: slide by sl seconds before determining playground data.
    @param excl_play: exclude the first excl_play second from the start and end
    of the chunk when computing if the chunk overlaps with playground.
    @param pad_data: exclude the first and last pad_data seconds of the segment
    when generating chunks
    """
    for seg in self.__sci_segs:
      seg.make_chunks(length,overlap,play,sl,excl_play,pad_data)

  def make_chunks_from_unused(self,length,trig_overlap,play=0,min_length=0,
    sl=0,excl_play=0,pad_data=0):
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
    @param pad_data: exclude the first and last pad_data seconds of the segment
    when generating chunks

    """
    for seg in self.__sci_segs:
      # if there is unused data longer than the minimum chunk length
      if seg.unused() > min_length:
        end = seg.end() - pad_data
        start = end - length
        if (not play) or (play and (((end - sl - excl_play - 729273613) % 6370)
          < (600 + length - 2 * excl_play))):
          trig_start = end - seg.unused() - trig_overlap
          if (play == 2):
            # calculate the start of the playground preceeding the chunk end
            play_start = 729273613 + 6370 * \
              math.floor((end - sl - excl_play - 729273613) / 6370)
            play_end = play_start + 600
            trig_end = 0
            if ( (play_end - 6370) > start ):
              print("Two playground segments in this chunk")
              print("  Code to handle this case has not been implemented")
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
        if (not play) or (play and (((end - sl - excl_play - 729273613) % 6370)
          < (600 + length - 2 * excl_play))):
          seg.add_chunk(start, end, start)
        seg.set_unused(0)

  def make_optimised_chunks(self, min_length, max_length, pad_data=0):
    """
    Splits ScienceSegments up into chunks, of a given maximum length.
    The length of the last two chunks are chosen so that the data
    utilisation is optimised.
    @param min_length: minimum chunk length.
    @param max_length: maximum chunk length.
    @param pad_data: exclude the first and last pad_data seconds of the
    segment when generating chunks
    """
    for seg in self.__sci_segs:
      # pad data if requested
      seg_start = seg.start() + pad_data
      seg_end = seg.end() - pad_data

      if seg.unused() > max_length:
        # get number of max_length chunks
        N = (seg_end - seg_start) / max_length

        # split into chunks of max_length
        for i in range(N - 1):
          start = seg_start + (i * max_length)
          stop = start + max_length
          seg.add_chunk(start, stop)

        # optimise data usage for last 2 chunks
        start = seg_start + ((N - 1) * max_length)
        middle = (start + seg_end) / 2
        seg.add_chunk(start, middle)
        seg.add_chunk(middle, seg_end)
        seg.set_unused(0)
      elif seg.unused() > min_length:
        # utilise as single chunk
        seg.add_chunk(seg_start, seg_end)
      else:
        # no chunk of usable length
        seg.set_unused(0)

  def intersection(self, other):
    """
    Replaces the ScienceSegments contained in this instance of ScienceData
    with the intersection of those in the instance other. Returns the number
    of segments in the intersection.
    @param other: ScienceData to use to generate the intersection
    """

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

          x = ScienceSegment(tuple([id, ostart, ostop, ostop - ostart]))
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
         x = ScienceSegment(tuple([id, ostart, ostop, ostop - ostart]))
         seglist.append(x)
         ostart = ustart
         ostop = ustop

    # flush out the final output segment (if any)
    if ostart != -1:
      x = ScienceSegment(tuple([id, ostart, ostop, ostop - ostart]))
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
          # the following line produces a flake8 issue with ostart; see https://git.ligo.org/lscsoft/glue/-/issues/37
          x = ScienceSegment(tuple([id, ostart, ostop, ostop - ostart]))   # noqa: F821
          outlist.append(x)
        ostart = start  # noqa: F841
        ostop = stop  # noqa: F841
      elif stop > ostop:
        # extend the current segment
        ostop = stop

    # flush out the final segment (if any)
    if ostop >= 0:
      x = ScienceSegment(tuple([id, ostart, ostop, ostop - ostart]))
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
      self.__sci_segs = ScienceSegment(tuple([0,0,1999999999,1999999999]))

    # go through the list checking for validity as we go
    outlist = []
    ostart = 0
    for seg in self:
      start = seg.start()
      stop = seg.end()
      if start < 0 or stop < start or start < ostart:
        raise SegmentError("Invalid list")
      if start > 0:
        x = ScienceSegment(tuple([0, ostart, start, start - ostart]))
        outlist.append(x)
      ostart = stop

    if ostart < 1999999999:
      x = ScienceSegment(tuple([0, ostart, 1999999999, 1999999999 - ostart]))
      outlist.append(x)

    self.__sci_segs = outlist
    return len(self)

  def play(self):
    """
    Keep only times in ScienceSegments which are in the playground
    """

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
      play_start = begin_s2 + play_space * ( 1
        + int((start - begin_s2 - play_len) / play_space) )

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

        x = ScienceSegment(tuple([id, ostart, ostop, ostop - ostart]))
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

  def intersect_4(self, second, third, fourth):
    """
     Intersection routine for four inputs.
    """
    self.intersection(second)
    self.intersection(third)
    self.intersection(fourth)
    self.coalesce()
    return len(self)

  def split(self, dt):
    """
      Split the segments in the list is subsegments at least as long as dt
    """
    outlist = []
    for seg in self:
      start = seg.start()
      stop = seg.end()
      id = seg.id()

      while start < stop:
        tmpstop = start + dt
        if tmpstop > stop:
          tmpstop = stop
        elif tmpstop + dt > stop:
          tmpstop = int( (start + stop) / 2 )
        x = ScienceSegment(tuple([id, start, tmpstop, tmpstop - start]))
        outlist.append(x)
        start = tmpstop

    # save the split list and return length
    self.__sci_segs = outlist
    return len(self)


class LsyncCache(object):
  def __init__(self,path):
    # location of the cache file
    self.__path = path

    # dictionary where the keys are data types like 'gwf', 'sft', 'xml'
    # and the values are dictionaries
    self.cache = {'gwf': None, 'sft': None, 'xml': None}

    # for each type create a dictionary where keys are sites and values
    # are dictionaries
    for type in self.cache.keys():
      self.cache[type] = {}

  def group(self, lst, n):
    """
    Group an iterable into an n-tuples iterable. Incomplete
    tuples are discarded
    """
    return itertools.izip(*[itertools.islice(lst, i, None, n) for i in range(n)])

  def parse(self,type_regex=None):
    """
    Each line of the frame cache file is like the following:

    /frames/E13/LHO/frames/hoftMon_H1/H-H1_DMT_C00_L2-9246,H,H1_DMT_C00_L2,1,16 1240664820 6231 {924600000 924646720 924646784 924647472 924647712 924700000}

    The description is as follows:

    1.1) Directory path of files
    1.2) Site
    1.3) Type
    1.4) Number of frames in the files (assumed to be 1)
    1.5) Duration of the frame files.

    2) UNIX timestamp for directory modification time.

    3) Number of files that that match the above pattern in the directory.

    4) List of time range or segments [start, stop)

    We store the cache for each site and frameType combination
    as a dictionary where the keys are (directory, duration)
    tuples and the values are segment lists.

    Since the cache file is already coalesced we do not
    have to call the coalesce method on the segment lists.
    """
    path = self.__path
    cache = self.cache
    if type_regex:
      type_filter = re.compile(type_regex)
    else:
      type_filter = None

    f = open(path, 'r')

    # holds this iteration of the cache
    gwfDict = {}

    # parse each line in the cache file
    for line in f:
      # ignore lines that don't match the regex
      if type_filter and type_filter.search(line) is None:
        continue

      # split on spaces and then comma to get the parts
      header, modTime, fileCount, times = line.strip().split(' ', 3)
      dir, site, frameType, frameCount, duration = header.split(',')
      duration = int(duration)

      # times string has form { t1 t2 t3 t4 t5 t6 ... tN t(N+1) }
      # where the (ti, t(i+1)) represent segments
      #
      # first turn the times string into a list of integers
      times = [ int(s) for s in times[1:-1].split(' ') ]

      # group the integers by two and turn those tuples into segments
      segs = [ segments.segment(a) for a in self.group(times, 2) ]

      # initialize if necessary for this site
      if site not in gwfDict:
        gwfDict[site] = {}

      # initialize if necessary for this frame type
      if frameType not in gwfDict[site]:
        gwfDict[site][frameType] = {}

      # record segment list as value indexed by the (directory, duration) tuple
      key = (dir, duration)
      if key in gwfDict[site][frameType]:
        msg = "The combination %s is not unique in the frame cache file" \
          % str(key)
        raise RuntimeError(msg)

      gwfDict[site][frameType][key] = segments.segmentlist(segs)
    f.close()

    cache['gwf'] = gwfDict

  def get_lfns(self, site, frameType, gpsStart, gpsEnd):
    """
    """
    # get the cache from the manager
    cache = self.cache

    # if the cache does not contain any mappings for this site type return empty list
    if site not in cache['gwf']:
      return []

    # if the cache does nto contain any mappings for this frame type return empty list
    if frameType not in cache['gwf'][site]:
      return []

    # segment representing the search interval
    search = segments.segment(gpsStart, gpsEnd)

    # segment list representing the search interval
    searchlist = segments.segmentlist([search])

    # dict of LFNs returned that match the metadata query
    lfnDict = {}

    for key,seglist in cache['gwf'][site][frameType].items():
      dir, dur = key

      # see if the seglist overlaps with our search
      overlap = seglist.intersects(searchlist)

      if not overlap: continue

      # the seglist does overlap with search so build file paths
      # but reject those outside of the search segment

      for s in seglist:
        if s.intersects(search):
          t1, t2 = s
          times = range(t1, t2, dur)

          # loop through the times and create paths
          for t in times:
            if search.intersects(segments.segment(t, t + dur)):
              lfn = "%s-%s-%d-%d.gwf" % (site, frameType, t, dur)
              lfnDict[lfn] = None

    # sort the LFNs to deliver URLs in GPS order
    lfns = list(lfnDict.keys())
    lfns.sort()

    return lfns


class LSCDataFindJob(CondorDAGJob, AnalysisJob):
  """
  An LSCdataFind job used to locate data. The static options are
  read from the section [datafind] in the ini file. The stdout from
  LSCdataFind contains the paths to the frame files and is directed to a file
  in the cache directory named by site and GPS start and end times. The stderr
  is directed to the logs directory. The job always runs in the scheduler
  universe. The path to the executable is determined from the ini file.
  """
  def __init__(self,cache_dir,log_dir,config_file,lsync_cache_file=None,lsync_type_regex=None):
    """
    @param cache_dir: the directory to write the output lal cache files to.
    @param log_dir: the directory to write the stderr file to.
    @param config_file: ConfigParser object containing the path to the LSCdataFind executable in the [condor] section and a [datafind] section from which the LSCdataFind options are read.
    @param lsync_cache_file: lsync_cache_file
    @param lsync_type_regex: lsync_type_regex
    """
    self.__executable = config_file.get('condor','datafind')
    self.__universe = 'local'
    CondorDAGJob.__init__(self,self.__universe,self.__executable)
    AnalysisJob.__init__(self,config_file)
    self.__cache_dir = cache_dir
    self.__config_file = config_file
    self.__lsync_cache = None
    if config_file.has_option('condor','accounting_group'):
        self.add_condor_cmd('accounting_group',config_file.get('condor','accounting_group'))
    if lsync_cache_file:
      self.__lsync_cache = LsyncCache(lsync_cache_file)
      self.__lsync_cache.parse(lsync_type_regex)

    # we have to do this manually for backwards compatibility with type
    for o in self.__config_file.options('datafind'):
      opt = str(o).strip()
      if opt[:4] != "type":
        arg = str(self.__config_file.get('datafind',opt)).strip()
        self.add_opt(opt,arg)

    # we need a lal cache for file PFNs
    self.add_opt('lal-cache','')
    self.add_opt('url-type','file')

    self.add_condor_cmd('getenv','True')

    self.set_stderr_file(os.path.join(log_dir, 'datafind-$(macroobservatory)-$(macrotype)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err'))
    self.set_stdout_file(os.path.join(log_dir, 'datafind-$(macroobservatory)-$(macrotype)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out'))
    self.set_sub_file('datafind.sub')

  def get_cache_dir(self):
    """
    returns the directroy that the cache files are written to.
    """
    return self.__cache_dir

  def get_config_file(self):
    """
    return the configuration file object
    """
    return self.__config_file

  def lsync_cache(self):
    return self.__lsync_cache


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
    self.__lfn_list = None

    # try and get a type from the ini file and default to type None
    try:
      self.set_type(self.job().get_config_file().get('datafind','type'))
    except:
      self.__type = None

  def __set_output(self):
    """
    Private method to set the file to write the cache to. Automaticaly set
    once the ifo, start and end times have been set.
    """
    if self.__start and self.__end and self.__observatory and self.__type:
      self.__output = os.path.join(self.__job.get_cache_dir(), self.__observatory + '-' + self.__type + '_CACHE' + '-' + str(self.__start) + '-' + str(self.__end - self.__start) + '.lcf')
      self.set_output(self.__output)

  def set_start(self, time, pad=None):
    """
    Set the start time of the datafind query.
    @param time: GPS start time of query.
    @param pad: pad
    """
    if pad:
      self.add_var_opt('gps-start-time', int(time) - int(pad))
    else:
      self.add_var_opt('gps-start-time', int(time))
    self.__start = time
    self.__set_output()

  def get_start(self):
    """
    Return the start time of the datafind query
    """
    return self.__start

  def set_end(self,time):
    """
    Set the end time of the datafind query.
    @param time: GPS end time of query.
    """
    self.add_var_opt('gps-end-time', time)
    self.__end = time
    self.__set_output()

  def get_end(self):
    """
    Return the start time of the datafind query
    """
    return self.__end

  def set_observatory(self,obs):
    """
    Set the IFO to retrieve data for. Since the data from both Hanford
    interferometers is stored in the same frame file, this takes the first
    letter of the IFO (e.g. L or H) and passes it to the --observatory option
    of LSCdataFind.
    @param obs: IFO to obtain data for.
    """
    self.add_var_opt('observatory',obs)
    self.__observatory = str(obs)
    self.__set_output()

  def get_observatory(self):
    """
    Return the start time of the datafind query
    """
    return self.__observatory

  def set_type(self,type):
    """
    sets the frame type that we are querying
    """
    self.add_var_opt('type',str(type))
    self.__type = str(type)
    self.__set_output()

  def get_type(self):
    """
    gets the frame type that we are querying
    """
    return self.__type

  def get_output_cache(self):
    return self.__output

  def get_output(self):
    """
    Return the output file, i.e. the file containing the frame cache data.
    or the files itself as tuple (for DAX)
    """
    return self.__output


class LigolwAddJob(CondorDAGJob, AnalysisJob):
  """
  A ligolw_add job can be used to concatenate several ligo lw files
  """
  def __init__(self,log_dir,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','ligolw_add')
    self.__universe = 'vanilla'
    CondorDAGJob.__init__(self,self.__universe,self.__executable)
    AnalysisJob.__init__(self,cp)
    self.add_ini_opts(cp, "ligolw_add")

    self.add_condor_cmd('getenv','True')
    if cp.has_option('condor','accounting_group'):
        self.add_condor_cmd('accounting_group',cp.get('condor','accounting_group'))

    self.set_stdout_file(os.path.join( log_dir, 'ligolw_add-$(cluster)-$(process).out') )
    self.set_stderr_file(os.path.join( log_dir, 'ligolw_add-$(cluster)-$(process).err') )
    self.set_sub_file('ligolw_add.sub')


class LigolwAddNode(CondorDAGNode, AnalysisNode):
  """
  Runs an instance of ligolw_add in a Condor DAG.
  """
  def __init__(self,job):
    """
    @param job: A CondorDAGJob that can run an instance of ligolw_add
    """
    CondorDAGNode.__init__(self,job)
    AnalysisNode.__init__(self)


class LigolwCutJob(CondorDAGJob, AnalysisJob):
  """
  A ligolw_cut job can be used to remove parts of a ligo lw file
  """
  def __init__(self,log_dir,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','ligolw_cut')
    self.__universe = 'vanilla'
    CondorDAGJob.__init__(self,self.__universe,self.__executable)
    AnalysisJob.__init__(self,cp)

    self.add_condor_cmd('getenv','True')

    self.set_stdout_file(os.path.join( log_dir, 'ligolw_cut-$(cluster)-$(process).out') )
    self.set_stderr_file(os.path.join( log_dir, 'ligolw_cut-$(cluster)-$(process).err') )
    self.set_sub_file('ligolw_cut.sub')


class LigolwCutNode(CondorDAGNode, AnalysisNode):
  """
  Runs an instance of ligolw_cut in a Condor DAG.
  """
  def __init__(self,job):
    """
    @param job: A CondorDAGJob that can run an instance of ligolw_cut
    """
    CondorDAGNode.__init__(self,job)
    AnalysisNode.__init__(self)


class NoopJob(CondorDAGJob, AnalysisJob):
  """
  A Noop Job does nothing.
  """
  def __init__(self,log_dir,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = 'true'
    self.__universe = 'local'
    CondorDAGJob.__init__(self,self.__universe,self.__executable)
    AnalysisJob.__init__(self,cp)

    self.add_condor_cmd('getenv','True')
    self.add_condor_cmd('noop_job','True')
    if cp.has_option('condor','accounting_group'):
        self.add_condor_cmd('accounting_group',cp.get('condor','accounting_group'))

    self.set_stdout_file(os.path.join( log_dir, 'noop-$(cluster)-$(process).out') )
    self.set_stderr_file(os.path.join( log_dir, 'noop-$(cluster)-$(process).err') )
    self.set_sub_file('noop.sub')


class NoopNode(CondorDAGNode, AnalysisNode):
  """
  Run an noop job in a Condor DAG.
  """
  def __init__(self,job):
    """
    @param job: A CondorDAGJob that does nothing.
    """
    CondorDAGNode.__init__(self,job)
    AnalysisNode.__init__(self)
    self.__server = None
    self.__identity = None
    self.__insert = None
    self.__pfn = None
    self.__query = None


class SqliteJob(CondorDAGJob, AnalysisJob):
  """
  A cbc sqlite job adds to CondorDAGJob and AnalysisJob features common to jobs
  which read or write to a sqlite database. Of note, the universe is always set to
  local regardless of what's in the cp file, the extension is set
  to None so that it may be set by individual SqliteNodes, log files do not
  have macrogpsstarttime and endtime in them, and get_env is set to True.
  """
  def __init__(self, cp, sections, exec_name):
    """
    @param cp: a ConfigParser object from which options are read
    @param sections: list of sections in cp to get added options
    @param exec_name: the name of the sql executable
    """
    self.__exec_name = exec_name
    executable = cp.get('condor', exec_name)
    universe = 'vanilla'
    CondorDAGJob.__init__(self, universe, executable)
    AnalysisJob.__init__(self, cp)

    for sec in sections:
      if cp.has_section(sec):
        self.add_ini_opts(cp, sec)
      else:
        sys.stderr.write("warning: config file is missing section [" + sec + "]\n")

    self.add_condor_cmd('getenv', 'True')
    if cp.has_option('condor','accounting_group'):
        self.add_condor_cmd('accounting_group',cp.get('condor','accounting_group'))
    self.set_stdout_file('logs/' + exec_name + '-$(cluster)-$(process).out')
    self.set_stderr_file('logs/' + exec_name + '-$(cluster)-$(process).err')

  def set_exec_name(self, exec_name):
    """
    Set the exec_name name
    """
    self.__exec_name = exec_name

  def get_exec_name(self):
    """
    Get the exec_name name
    """
    return self.__exec_name


class SqliteNode(CondorDAGNode, AnalysisNode):
  """
  A cbc sqlite node adds to the standard AnalysisNode features common to nodes
  which read or write to a sqlite database. Specifically, it adds the set_tmp_space_path
  and set_database methods.
  """
  def __init__(self, job):
    """
    @param job: an Sqlite job
    """
    CondorDAGNode.__init__(self, job)
    AnalysisNode.__init__(self)
    self.__tmp_space = None
    self.__database = None

  def set_tmp_space(self, tmp_space):
    """
    Sets temp-space path. This should be on a local disk.
    @param tmp_space: tmp_space
    """
    self.add_var_opt('tmp-space', tmp_space)
    self.__tmp_space = tmp_space

  def get_tmp_space(self):
    """
    Gets tmp-space path.
    """
    return self.__tmp_space

  def set_database(self, database):
    """
    Sets database option.
    @param database: database
    """
    self.add_file_opt('database', database)
    self.__database = database

  def get_database(self):
    """
    Gets database option.
    """
    return self.__database


class LigolwSqliteJob(SqliteJob):
  """
  A LigolwSqlite job. The static options are read from the
  section [ligolw_sqlite] in the ini file.
  """
  def __init__(self, cp):
    """
    @param cp: ConfigParser object from which options are read.
    """
    exec_name = 'ligolw_sqlite'
    sections = ['ligolw_sqlite']
    super(LigolwSqliteJob,self).__init__(cp, sections, exec_name)

  def set_replace(self):
    """
    Sets the --replace option. This will cause the job
    to overwrite existing databases rather than add to them.
    """
    self.add_opt('replace','')


class LigolwSqliteNode(SqliteNode):
  """
  A LigolwSqlite node.
  """
  def __init__(self, job):
    """
    @param job: a LigolwSqliteJob
    """
    super(LigolwSqliteNode,self).__init__(job)
    self.__input_cache = None
    self.__xml_output = None
    self.__xml_input = None

  def set_input_cache(self, input_cache):
    """
    Sets input cache.
    @param input_cache: input_cache
    """
    self.add_file_opt('input-cache', input_cache)
    self.__input_cache = input_cache

  def get_input_cache(self):
    """
    Gets input cache.
    """
    return self.__input_cache

  def set_xml_input(self, xml_file):
    """
    Sets xml input file instead of cache
    @param xml_file: xml_file
    """
    self.add_var_arg(xml_file)

  def set_xml_output(self, xml_file):
    """
    Tell ligolw_sqlite to dump the contents of the database to a file.
    @param xml_file: xml_file
    """
    if self.get_database() is None:
      raise ValueError("no database specified")
    self.add_file_opt('extract', xml_file)
    self.__xml_output = xml_file

  def get_output(self):
    """
    Override standard get_output to return xml-file if xml-file is specified.
    Otherwise, will return database.
    """
    if self.__xml_output:
      return self.__xml_output
    elif self.get_database():
      return self.get_database()
    else:
      raise ValueError("no output xml file or database specified")


class DeepCopyableConfigParser(configparser.ConfigParser):
    """
    The standard SafeConfigParser no longer supports deepcopy() as of python
    2.7 (see http://bugs.python.org/issue16058). This subclass restores that
    functionality.
    """
    def __deepcopy__(self, memo):
        # http://stackoverflow.com/questions/23416370
        # /manually-building-a-deep-copy-of-a-configparser-in-python-2-7
        from io import StringIO
        config_string = StringIO()
        self.write(config_string)
        config_string.seek(0)
        new_config = self.__class__()
        new_config.readfp(config_string)
        return new_config
