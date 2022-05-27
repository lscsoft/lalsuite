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
import os
import re
import time
import random
import stat
from hashlib import md5
import warnings

warnings.warn(
	"all functionality within this module exists to support deprecated pipeline "
	"generation programs, and will be removed from gstlal in the future",
	DeprecationWarning,
)


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
      if self.__grid_type == None:
        raise CondorSubmitError('No grid type specified.')
      elif self.__grid_type == 'gt2':
        if self.__grid_server == None:
          raise CondorSubmitError('No server specified for grid resource.')
      elif self.__grid_type == 'gt4':
        if self.__grid_server == None:
          raise CondorSubmitError('No server specified for grid resource.')
        if self.__grid_scheduler == None:
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
    self.__dag_directory= dir

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
    if isinstance(job,CondorDAGJob) and job.get_universe()=='standard':
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
    fh.write( 'CATEGORY ' + self.__name + ' ' + self.__category +  '\n' )

  def write_priority(self,fh):
    """
    Write the DAG entry for this node's priority to the DAG file descriptor.
    @param fh: descriptor of open DAG file.
    """
    fh.write( 'PRIORITY ' + self.__name + ' ' + self.__priority +  '\n' )

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
    fh.write( 'MAXJOBS ' + str(category[0]) + ' ' + str(category[1]) +  '\n' )

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
