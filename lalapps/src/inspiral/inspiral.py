"""
Classes needed for the inspiral analysis pipeline.
This script produced the necessary condor submit and dag files to run
the standalone inspiral code on LIGO data
"""

__author__ = 'Duncan Brown <duncan@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

import string
import exceptions
from glue import pipeline


class InspiralError(exceptions.Exception):
  def __init__(self, args=None):
    self.args = args


class TmpltBankJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_tmpltbank job used by the inspiral pipeline. The static options
  are read from the sections [data] and [tmpltbank] in the ini file. The
  stdout and stderr from the job are directed to the logs directory. The job
  runs in the universe specfied in the ini file. The path to the executable
  is determined from the ini file.
  """
  def __init__(self,cp,dax=False,tag_base='TMPLTBANK'):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','tmpltbank')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp,dax)
    self.tag_base = tag_base

    for sec in ['data','tmpltbank']:
      try: self.add_ini_opts(cp,sec)
      except: pass
  
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")

    self.set_stdout_file('logs/tmpltbank-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/tmpltbank-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('tmpltbank.sub')


class InspInjJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_inspinj job used by the grb inspiral pipeline. The static options
  are read from the section [inspinj] in the ini file. The
  stdout and stderr from the job are directed to the logs directory. The
  job runs in the universe specified in the ini file. The path to the
  executable is determined from the ini file.
  """
  def __init__(self,cp,dax=False):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','inspinj')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp,dax)
    self.__listDone=[]
    self.__listNodes=[]

    for sec in ['inspinj']:
      try: self.add_ini_opts(cp,sec)
      except: pass

    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")

    self.set_stdout_file('logs/inspinj-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/inspinj-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')

  def set_done(self, number, node):
    self.__listDone.append(number)
    self.__listNodes.append(node)

  def check_node(self, number):
    if self.__listDone.count(number):
      index=self.__listDone.index(number)
      return self.__listNodes[index]
    return None    


class BbhInjJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_bbhinj job used by the online inspiral pipeline. The static options
  are read from the section [bbhinj] in the ini file. The
  stdout and stderr from the job are directed to the logs directory. The
  job runs in the universe specified in the ini file. The path to the 
  executable is determined from the ini file.
  """
  def __init__(self,cp,dax=False):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','bbhinj')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp,dax)

    for sec in ['bbhinj']:
      try: self.add_ini_opts(cp,sec)
      except: pass

    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")

    self.set_stdout_file('logs/bbhinj-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/bbhinj-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')


class RandomBankJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_randombank job used by the inspiral pipeline. The static options
  are read from the section [randombank] in the ini file. The stdout and
  stderr from the job are directed to the logs directory. The job runs in the
  universe specfied in the ini file. The path to the executable is determined
  from the ini file.
  """
  def __init__(self,cp,dax=False):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','randombank')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp,dax)

    try: self.add_ini_opts(cp,'randombank')
    except: pass

    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")

    self.set_stdout_file('logs/randombank-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/randombank-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('randombank.sub')


class SplitBankJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_splitbank job used by the inspiral pipeline. The static options
  are read from the section [splitbank] in the ini file. The stdout and stderr
  from the job are directed to the logs directory. The job runs in the
  universe specfied in the ini file. The path to the executable is determined
  from the ini file.
  """
  def __init__(self,cp,dax=False):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','splitbank')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp,dax)

    for sec in ['splitbank']:
      try: self.add_ini_opts(cp,sec)
      except: pass
  
    self.set_stdout_file('logs/splitbank-$(macrobankfile)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/splitbank-$(macrobankfile)-$(cluster)-$(process).err')
    self.set_sub_file('splitbank.sub')
    

class InspiralJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_inspiral job used by the inspiral pipeline. The static options
  are read from the sections [data] and [inspiral] in the ini file. The
  stdout and stderr from the job are directed to the logs directory. The job
  runs in the universe specfied in the ini file. The path to the executable
  is determined from the ini file.
  """
  def __init__(self,cp,dax=False,tag_base='INSPIRAL'):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','inspiral')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp,dax)
    self.tag_base = tag_base


    for sec in ['data','inspiral']:
      try: self.add_ini_opts(cp,sec)
      except: pass

    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")

    self.set_stdout_file('logs/inspiral-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/inspiral-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('inspiral.sub')


class TrigToTmpltJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_trigtotmplt job used by the inspiral pipeline. The static
  options are read from the section [trigtotmplt] in the ini file.  The
  stdout and stderr from the job are directed to the logs directory. The job
  always runs in the scheduler universe. The path to the executable is
  determined from the ini file.
  """
  def __init__(self,cp,dax=False):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','trigtotmplt')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp,dax)
    
    for sec in ['trigtotmplt']:
      try: self.add_ini_opts(cp,sec)
      except: pass

    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
    
    self.set_stdout_file('logs/trigtotmplt-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/trigtotmplt-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('trigtotmplt.sub')


class IncaJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_inca job used by the inspiral pipeline. The static options are
  read from the section [inca] in the ini file.  The stdout and stderr from
  the job are directed to the logs directory.  The path to the executable is 
  determined from the ini file.
  """
  def __init__(self,cp,dax=False):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','inca')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp,dax)
    
    for sec in ['inca']:
      try: self.add_ini_opts(cp,sec)
      except: pass

    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")

    self.set_stdout_file('logs/inca-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/inca-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('inca.sub')


class ThincaJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_thinca job used by the inspiral pipeline. The static options are
  read from the section [thinca] in the ini file.  The stdout and stderr from
  the job are directed to the logs directory.  The path to the executable is 
  determined from the ini file.
  """
  def __init__(self,cp,dax=False,tag_base='THINCA'):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','thinca')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp,False)
    self.tag_base = tag_base
    
    for sec in ['thinca']:
      try: self.add_ini_opts(cp,sec)
      except: pass

    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")

    self.set_stdout_file('logs/thinca-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/thinca-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('thinca.sub')


class SireJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_sire job used by the inspiral pipeline. The stdout and stderr from
  the job are directed to the logs directory. The path to the executable is 
  determined from the ini file.
  """
  def __init__(self,cp,dax=False):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','sire')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp,dax)

    for sec in ['sire']:
      try: self.add_ini_opts(cp,sec)
      except: pass

    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
    self.set_stdout_file('logs/sire-$(macroifo)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/sire-$(macroifo)-$(cluster)-$(process).err')
    self.set_sub_file('sire.sub')

class CoireJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_coire job used by the inspiral pipeline. The stdout and stderr from
  the job are directed to the logs directory. The path to the executable is
  determined from the ini file.
  """
  def __init__(self,cp,dax=False):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','coire')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp,dax)

    for sec in ['coire']:
      try: self.add_ini_opts(cp,sec)
      except: pass

    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
    self.set_stdout_file('logs/coire-$(macrocoinccut)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/coire-$(macrocoinccut)-$(cluster)-$(process).err')
    self.set_sub_file('coire.sub')
    
class FrJoinJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_frjoin job used by the inspiral pipeline. The path to the
  executable is determined from the ini file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','frjoin')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)
    
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")

    self.set_stdout_file('logs/frjoin-$(cluster)-$(process).out')
    self.set_stderr_file('logs/frjoin-$(cluster)-$(process).err')
    self.set_sub_file('frjoin.sub')

class CohBankJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_coherent_inspiral job used by the inspiral pipeline. The static
  options are read from the section [cohbank] in the ini file.  The stdout and
  stderr from the job are directed to the logs directory.  The path to the
  executable is determined from the ini file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','cohbank')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)
    
    for sec in ['cohbank']:
      try: self.add_ini_opts(cp,sec)
      except: pass

    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")

    self.set_stdout_file('logs/cohbank-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/cohbank-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('cohbank.sub')

class ChiaJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_coherent_inspiral job used by the inspiral pipeline. The static
  options are read from the section [chia] in the ini file.  The stdout and
  stderr from the job are directed to the logs directory.  The path to the
  executable is determined from the ini file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','chia')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)
    
    for sec in ['chia']:
      try: self.add_ini_opts(cp,sec)
      except: pass

    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")

    self.set_stdout_file('logs/chia-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/chia-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('chia.sub')   



class InspInjNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A InspInjNode runs an instance of the inspinj generation job in a
  Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_inspinj.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__usertag = job.get_config('pipeline','user-tag')

  def set_seed(self,seed):
    """
    Set the seed of the injection file by setting a --seed option to the
    node when it is executed. The seed is automatically the number of
    the injection 'round'.
    @param seed: seed of the job
    """
    self.add_var_opt('seed',seed)
    self.__seed = seed

  def set_output(self, outputName):
    """
    Set the output name of the injection file
    @param outputName: name of the injection file created
    """
    self.add_var_opt('output',outputName)
    self.__outputName = outputName

class BbhInjNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A BbhInjNode runs an instance of the bbhinj generation job in a 
  Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_bbhinj.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__usertag = job.get_config('pipeline','user-tag')

  def set_seed(self,seed):
    """
    Set the seed of the injection file by setting a --seed option to the
    node when it is executed.
    @param seed: seed of the job
    """
    self.add_var_opt('seed',seed)
    self.__seed = seed

  def get_output(self):
    """
    Returns the file name of output from the injection generation code. This 
    must be kept synchronized with the name of the output file in bbhinj.c.
    """
    if not self.get_start() or not self.get_end():
      raise InspiralError, "Start time or end time has not been set"
    if self.__usertag:
      bbhinject = 'HL-INJECTIONS_' + self.__usertag + '-'
      bbhinject = bbhinject + str(self.get_start()) + '-'
      bbhinject = bbhinject + str(self.get_end()-self.get_start()) + '.xml'
    elif self.__seed:
      bbhinject = 'HL-INJECTIONS_' + str(self.__seed) + '-'
      bbhinject = bbhinject + str(self.get_start()) + '-'
      bbhinject = bbhinject + str(self.get_end()-self.get_start()) + '.xml'
    else:
      bbhinject = 'HL-INJECTIONS-' + str(self.get_start()) + '-'
      bbhinject = bbhinject + str(self.get_end()-self.get_start()) + '.xml'

    self.add_output_file(bbhinject)

    return bbhinject


class TmpltBankNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A TmpltBankNode runs an instance of the template bank generation job in a
  Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_tmpltbank.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    try:
      self.__usertag = job.get_config('tmpltbank','user-tag')
    except:
      self.__usertag = job.get_config('pipeline','user-tag')

    try:
      self.__pad_data = int(self.job().get_opts()['pad-data'])
    except: 
      self.__pad_data = None

    try:
      self.__zip_output = job.get_config('tmpltbank','write-compress')
      self.__zip_output = True
    except:
      self.__zip_output = False

  def get_output(self):
    """
    Returns the file name of output from the template bank code. This must
    be kept synchronized with the name of the output file in tmpltbank.c.
    """
    tag_base = self.job().tag_base
    if not self.get_start() or not self.get_end() or not self.get_ifo():
      raise InspiralError, "Start time, end time or ifo has not been set"
    if self.__usertag and self.get_ifo_tag():
      bank = self.get_ifo() + '-' + tag_base + '_' + self.get_ifo_tag() + \
          "_" + self.__usertag + '-' 
      bank = bank + str(self.get_start())
    elif self.__usertag:
      bank = self.get_ifo() + '-' + tag_base + '_' + self.__usertag + '-'  
      bank = bank + str(self.get_start())
    elif self.get_ifo_tag():
      bank = self.get_ifo() + '-' + tag_base + '_' + self.get_ifo_tag() + '-'  
      bank = bank + str(self.get_start())
    else:
      bank = self.get_ifo() + '-' + tag_base + '-' + str(self.get_start())
    bank = bank + '-' + str(self.get_end() - self.get_start()) + '.xml'

    if self.__zip_output:
      bank += '.gz'

    self.add_output_file(bank)

    return bank

  def set_pad_data(self, pad):
    """
    Set the pad data value for this node 
    """
    self.__pad_data = pad
    self.add_var_opt('pad-data', pad)

  def get_pad_data(self):
    """
    Returns the injection file
    """
    return self.__pad_data

  def finalize(self):
    """
    set the data_start_time and data_end_time
    """
    if self.get_pad_data():
      pipeline.AnalysisNode.set_data_start(self,self.get_start() - \
          self.get_pad_data())
      pipeline.AnalysisNode.set_data_end(self,self.get_end() + \
          self.get_pad_data())

class RandomBankNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A RandomBankNode runs an instance of the random bank generation job in a
  Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_randombank.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    try:
      self.__usertag = job.get_config('randombank','user-tag')
    except:
      self.__usertag = job.get_config('pipeline','user-tag')

  def get_output(self):
    """
    Returns the file name of output from the template bank code. This must
    be kept synchronized with the name of the output file in randombank.c.
    """
    if not self.get_start() or not self.get_end():
      raise InspiralError, "Start time or end time has not been set"
    if self.__usertag:
      bank = 'P-TMPLTBANK_' + self.__usertag + '-' 
      bank = bank + str(self.get_start())
    else:
      bank = 'P-TMPLTBANK-' + str(self.get_start())
    bank = bank + '-' + str(self.get_end() - self.get_start()) + '.xml'

    self.add_output_file(bank)

    return bank


class SplitBankNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A SplitBankNode runs an instance of the split template bank job in a
  Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_tmpltbank.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    try:
      self.__usertag = job.get_config('splitbank','user-tag')
    except:
      self.__usertag = job.get_config('pipeline','user-tag')
    self.__bankfile = None
    self.__numbanks = None

  def set_bank(self,bank):
    self.add_var_opt('bank-file', bank)
    self.add_input_file(bank)
    self.__bankfile = bank

  def get_bank(self):
    return self.__bankfile

  def set_num_banks(self,numbanks):
    self.add_var_opt('number-of-banks',numbanks)
    self.__numbanks = int(numbanks)

  def get_num_banks(self):
    return self.__numbanks

  def get_output(self):
    """
    Returns a list of the file names of split banks. This must be kept
    synchronized with the name of the output files in splitbank.c.
    """
    if not self.get_bank() or not self.get_num_banks():
      raise InspiralError, "Bank file or number of banks has not been set"

    banks = []
    x = self.__bankfile.split('-')
    for i in range( 0, int(self.get_num_banks()) ):
      banks.append("%s-%s_%2.2d-%s-%s" % (x[0], x[1], i, x[2], x[3]))

    return banks


class InspiralNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  An InspiralNode runs an instance of the inspiral code in a Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_inspiral.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    try:
      self.__usertag = job.get_config('inspiral','user-tag')
    except:
      self.__usertag = job.get_config('pipeline','user-tag')

    self.__injections = None
    try:
      self.__pad_data = int(self.job().get_opts()['pad-data'])
    except: 
      self.__pad_data = None
    try:
      self.__zip_output = job.get_config('inspiral','write-compress')
      self.__zip_output = True
    except:
      self.__zip_output = False 

  def set_bank(self,bank):
    self.add_var_opt('bank-file', bank)
    self.add_input_file(bank)

  def set_injections(self, injections):
    """
    Set the injection file for this node
    """
    self.__injections = injections
    self.add_var_opt('injection-file', injections)
    self.add_input_file(injections)

  def get_injections(self):
    """
    Returns the injection file
    """
    return self.__injections

  def set_pad_data(self, pad):
    """
    Set the pad data value for this node 
    """
    self.__pad_data = pad
    self.add_var_opt('pad-data', pad)

  def get_pad_data(self):
    """
    Returns the injection file
    """
    return self.__pad_data

  def set_user_tag(self,usertag):
    self.__usertag = usertag
    self.add_var_opt('user-tag',usertag)

  def get_user_tag(self):
    return self.__usertag

  def get_output(self):
    """
    Returns the file name of output from the inspiral code. This must be kept
    synchronized with the name of the output file in inspiral.c.
    """
    if not self.get_start() or not self.get_end() or not self.get_ifo():
      raise InspiralError, "Start time, end time or ifo has not been set"

    tag_base = self.job().tag_base
    basename = self.get_ifo() + '-' + tag_base

    if self.get_ifo_tag():
      basename += '_' + self.get_ifo_tag()
    if self.__usertag:
      basename += '_' + self.__usertag

    filename = basename + '-' + str(self.get_start()) + '-' + \
      str(self.get_end() - self.get_start()) + '.xml'

    if self.__zip_output:
      filename += '.gz'

    self.add_output_file(filename)

    return filename

  def get_froutput(self):
    """
    Returns the file name of output frame from the inspiral code. This
    must be kept synchronized with the name of the output file in inspiral.c.
    """
    if not self.get_start() or not self.get_end() or not self.get_ifo():
      raise InspiralError, "Start time, end time or ifo has not been set"

    tag_base = self.job().tag_base
    basename = self.get_ifo() + '-' + tag_base

    if self.get_ifo_tag():
      basename += '_' + self.get_ifo_tag()
    if self.__usertag:
      basename += '_' + self.__usertag

    filename = basename + '-' + str(self.get_start()) + '-' + \
      str(self.get_end() - self.get_start()) + '.gwf'

    return filename  

  def finalize(self):
    """
    set the data_start_time and data_end_time
    """
    if self.get_pad_data():
      pipeline.AnalysisNode.set_data_start(self,self.get_start() - \
          self.get_pad_data())
      pipeline.AnalysisNode.set_data_end(self,self.get_end() + \
          self.get_pad_data())

class TrigToTmpltNode(pipeline.CondorDAGNode,pipeline.AnalysisNode):
  """
  A TrigToTmpltNode runs an instance of the triggered bank generator in a
  Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of inca in trigtotmplt mode.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__output = None
    self.__input_ifo = None
    self.__output_ifo = None
    try:
      self.__usertag = job.get_config('trigtotmplt','user-tag')
    except:
      self.__usertag = job.get_config('pipeline','user-tag')
    try:
      self.__zip_output = job.get_config('trigtotmplt','write-compress')
      self.__zip_output = True
    except:
      self.__zip_output = False

  def set_user_tag(self,usertag):
    self.__usertag = usertag
    self.add_var_opt('user-tag',usertag)

  def get_user_tag(self):
    return self.__usertag

  def set_input_ifo(self,ifo):
    self.add_var_opt('input-ifo', ifo)
    self.__input_ifo = ifo

  def get_input_ifo(self):
    return self.__input_ifo

  def set_output_ifo(self,ifo):
    self.add_var_opt('output-ifo', ifo)
    self.__output_ifo = ifo

  def get_output_ifo(self):
    return self.__output_ifo

  def get_output(self):
    """
    Returns the name of the output file from lalapps_trigbank
    """
    if not self.get_start() or not self.get_end() or not self.get_output_ifo():
      raise InspiralError, "Start time, end time or output ifo is not set"
      
    basename = self.get_output_ifo() + '-TRIGBANK'

    if self.get_ifo_tag():
      basename += '_' + self.get_ifo_tag()
    if self.__usertag:
      basename += '_' + self.__usertag 

    trigbank_name = basename + '-' + str(self.get_start()) + '-' + \
      str(self.get_end() - self.get_start()) + '.xml'

    if self.__zip_output:
      trigbank_name += '.gz'

    self.add_output_file(trigbank_name)
    return trigbank_name


class IncaNode(pipeline.CondorDAGNode,pipeline.AnalysisNode):
  """
  An IncaNode runs an instance of the inspiral coincidence code in a Condor
  DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_inca.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__ifo_a = None
    self.__ifo_b = None
    try:
      self.__usertag = job.get_config('inca','user-tag')
    except:
      self.__usertag = job.get_config('pipeline','user-tag')
    try:
      self.__zip_output = job.get_config('inca','write-compress')
      self.__zip_output = True
    except:
      self.__zip_output = False

  def set_ifo_a(self, ifo):
    """
    Set the interferometer code to use as IFO A.
    ifo = IFO code (e.g. L1, H1 or H2).
    """
    self.add_var_opt('ifo-a', ifo)
    self.__ifo_a = ifo

  def get_ifo_a(self):
    """
    Returns the IFO code of the primary interferometer.
    """
    return self.__ifo_a

  def set_ifo_b(self, ifo):
    """
    Set the interferometer code to use as IFO B.
    ifo = IFO code (e.g. L1, H1 or H2).
    """
    self.add_var_opt('ifo-b', ifo)
    self.__ifo_b = ifo

  def get_ifo_b(self):
    """
    Returns the IFO code of the primary interferometer.
    """
    return self.__ifo_b

  def set_user_tag(self,usertag):
    """
    Set the usertag for a given job
    """
    self.__usertag = usertag
    self.add_var_opt('user-tag',usertag)

  def get_user_tag(self):
    """
    Returns the usertag of the job
    """
    return self.__usertag

  def get_output_a(self):
    """
    Returns the file name of output from inca for ifo a. This must be kept
    synchronized with the name of the output file in inca.c.
    """
    if not self.get_start() or not self.get_end() or not self.get_ifo_a():
      raise InspiralError, "Start time, end time or ifo a has not been set"

    basename = self.get_ifo_a() + '-INCA'

    if self.get_ifo_tag():
      basename += '_' + self.get_ifo_tag()
    if self.__usertag:
      basename += '_' + self.__usertag 

    filename = basename + '-' + str(self.get_start()) + '-' + \
      str(self.get_end() - self.get_start()) + '.xml'

    if self.__zip_output:
      filename += '.gz'

    self.add_output_file(filename)
    return filename

  def get_output_b(self):
    """
    Returns the file name of output from inca for ifo b. This must be kept
    synchronized with the name of the output file in inca.c.
    """
    if not self.get_start() or not self.get_end() or not self.get_ifo_b():
      raise InspiralError, "Start time, end time or ifo a has not been set"

    basename = self.get_ifo_b() + '-INCA'

    if self.get_ifo_tag():
      basename += '_' + self.get_ifo_tag()
    if self.__usertag:
      basename += '_' + self.__usertag 

    filename = basename + '-' + str(self.get_start()) + '-' + \
      str(self.get_end() - self.get_start()) + '.xml'

    if self.__zip_output:
      filename += '.gz'

    self.add_output_file(filename)
    return filename


class ThincaNode(pipeline.CondorDAGNode,pipeline.AnalysisNode):
  """
  A ThincaNode runs an instance of the inspiral coincidence code in a Condor
  DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_inca.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__ifo_g1 = None
    self.__ifo_h1 = None
    self.__ifo_h2 = None
    self.__ifo_l1 = None
    self.__ifo_t1 = None
    self.__ifo_v1 = None
    self.__num_slides = None
    try:
      self.__usertag = job.get_config('thinca','user-tag')
    except:
      self.__usertag = job.get_config('pipeline','user-tag')
    self.__ifotag = None
    try:
      self.__zip_output = job.get_config('thinca','write-compress')
      self.__zip_output = True
    except:
      self.__zip_output = False

  def set_ifo(self, ifo):
    """
    Add the interferometer to the list of ifos
    ifo = IFO code (e.g. G1,L1, H1 or H2).
    """
    if ifo == 'G1':
      self.add_var_opt('g1-triggers','')
      self.__ifo_g1 = 'G1'
    elif ifo == 'H1':
      self.add_var_opt('h1-triggers','')
      self.__ifo_h1 = 'H1'
    elif ifo == 'H2':
      self.add_var_opt('h2-triggers','')
      self.__ifo_h2 = 'H2'
    elif ifo == 'L1':
      self.add_var_opt('l1-triggers','')
      self.__ifo_l1 = 'L1'
    elif ifo == 'T1':
      self.add_var_opt('t1-triggers','')
      self.__ifo_t1 = 'T1'
    elif ifo == 'V1':
      self.add_var_opt('v1-triggers','')
      self.__ifo_v1 = 'V1'
  
  def get_ifo_g1(self):
    """
    Returns the IFO code of g1.
    """
    return self.__ifo_g1
    
  def get_ifo_h1(self):
    """
    Returns the IFO code of h1.
    """
    return self.__ifo_h1

  def get_ifo_h2(self):
    """
    Returns the IFO code of h2.
    """
    return self.__ifo_h2

  def get_ifo_l1(self):
    """
    Returns the IFO code of l1.
    """
    return self.__ifo_l1

  def get_ifo_t1(self):
    """
    Returns the IFO code of t1.
    """
    return self.__ifo_t1

  def get_ifo_v1(self):
    """
    Returns the IFO code of v1.
    """
    return self.__ifo_v1

  def get_ifos(self):
    """
    Returns the ordered list of ifos.
    """
    ifos = ''
    if self.get_ifo_g1():
      ifos += self.get_ifo_g1()
    if self.get_ifo_h1():
      ifos += self.get_ifo_h1()
    if self.get_ifo_h2():
      ifos += self.get_ifo_h2()
    if self.get_ifo_l1():
      ifos += self.get_ifo_l1()
    if self.get_ifo_t1():
      ifos += self.get_ifo_t1()
    if self.get_ifo_v1():
      ifos += self.get_ifo_v1()

    return ifos

  def set_num_slides(self, num_slides):
    """
    Set number of time slides to undertake
    """
    self.add_var_opt('num-slides',num_slides)
    self.__num_slides = num_slides

  def get_num_slides(self):
    """
    Returns the num_slides from .ini (>0 => time slides desired)
    """
    return self.__num_slides

  def set_user_tag(self,usertag):
    """
    Set the usertag for a given job
    """
    self.__usertag = usertag
    self.add_var_opt('user-tag',usertag)

  def get_user_tag(self):
    """
    Returns the usertag of the job
    """
    return self.__usertag

  def set_ifo_tag(self,ifotag):
    """
    Set the ifotag for a given job (for second thinca)
    """
    self.__ifotag = ifotag
    self.add_var_opt('ifo-tag',ifotag)

  def get_ifo_tag(self):
    """
    Returns the ifo tag of the job
    """
    return self.__ifotag

  def get_output(self):
    """
    Returns the file name of output from thinca.  This must be kept
    synchronized with the name of the output file in thinca.c.
    """
    if not self.get_start() or not self.get_end() or not self.get_ifos():
      raise InspiralError, "Start time, end time or ifos have not been set"
    
    tag_base = self.job().tag_base
    if self.__num_slides:
      basename = self.get_ifos() + '-' + tag_base + '_SLIDE'
    else:
      basename = self.get_ifos() + '-' + tag_base

    if self.__ifotag:
      basename += '_' + self.__ifotag  

    if self.__usertag:
      basename += '_' + self.__usertag

    filename = basename + '-' + str(self.get_start()) + '-' + \
      str(self.get_end() - self.get_start()) + '.xml'

    if self.__zip_output:
      filename += '.gz'

    self.add_output_file(filename)
    return filename


class SireNode(pipeline.CondorDAGNode,pipeline.AnalysisNode):
  """
  A SireNode runs an instance of the single inspiral reader code in a Condor
  DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_sire.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__ifo  = None
    self.__ifotag = None
    self.__start = None
    self.__end   = None
    self.__injection_file = None
    try:
      self.__usertag = job.get_config('sire','user-tag')
    except:
      self.__usertag = job.get_config('pipeline','user-tag')

    try:
      self.__zip_output = job.get_config('coire','write-compress')
      self.__zip_output = True
    except:
      self.__zip_output = False

  def set_ifo(self, ifo):
    """
    Add the list of interferometers 
    """
    self.__ifo = ifo
    self.add_var_opt('ifo-cut',ifo)

  def get_ifo(self):
    """
    Returns the ifos
    """
    return self.__ifo

  def set_inj_file(self, file):
    """
    Sets the injection file
    """
    self.__injection_file = file
    self.add_var_opt('injection-file', file)

  def get_inj_file(self, file):
    """
    Sets the injection file
    """
    return self.__injection_file

  def set_ifo_tag( self, ifotag):
    """
    Set the ifotag
    """
    self.__ifotag = ifotag

  def get_ifo_tag( self ):
    """
    get the ifotag
    """
    return self.__ifotag

  def set_start(self, start):
    """
    Sets GPS start time
    """
    self.__start = start

  def get_start(self):
    """
    Returns GPS start time
    """
    return self.__start

  def set_end(self, end):
    """
    Sets GPS end time
    """
    self.__end = end

  def get_end(self):
    """
    Returns GPS end time
    """
    return self.__end

  def set_glob(self, file_glob):
    """
    Sets the glob name
    """
    self.add_var_opt('glob',file_glob)

  def set_input(self, input_file):
    """
    Sets the input file name
    """
    self.add_var_opt('input',input_file)

  def set_user_tag(self,usertag):
    """
    Set the usertag for a given job
    """
    self.__usertag = usertag
    self.add_var_opt('user-tag',usertag)

  def get_user_tag(self):
    """
    Returns the usertag of the job
    """
    return self.__usertag

  def get_output(self):
    """
    get the name of the output file
    """
    if not self.get_ifo():
      raise InspiralError, "ifos have not been set"

    fname = self.__ifo + "-SIRE"
    if self.__injection_file:
      fname += "_" + self.__injection_file.split("-")[1]
      fname += "_FOUND"

    if self.__ifotag: fname += "_" + self.__ifotag
    if self.__usertag: fname += "_" + self.__usertag

    if (self.__start and not self.__end) or (self.__end and not self.__start):
      raise InspiralError, "If one of start and end is set, both must be"

    if (self.__start):
      duration=self.__end - self.__start
      fname += "-" + str(self.__start) + "-" + str(duration)

    fname += ".xml"

    return fname

  def get_missed(self):
    """
    get the name of the missed file
    """
    if self.__injection_file:
      return self.get_output().replace("FOUND", "MISSED")
    else:
      return None

  def finalize(self):
    """
    set the output options
    """
    output = self.get_output()
 
    self.add_var_opt("output", output)
    self.add_var_opt("summary", output.replace("xml", "txt"))

    if self.__injection_file:
      self.add_var_opt('injection-file', self.__injection_file)
      self.add_var_opt('missed-injections', self.get_missed() )


class CoireNode(pipeline.CondorDAGNode,pipeline.AnalysisNode):
  """
  A CoireNode runs an instance of the inspiral coire code in a Condor
  DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_coire.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__ifos  = None
    self.__ifotag = None
    self.__start = None
    self.__end   = None
    self.__num_slides = None
    self.__injection_file = None
    self.__output_tag = None
    try:
      self.__usertag = job.get_config('coire','user-tag')
    except:
      self.__usertag = job.get_config('pipeline','user-tag')

  def set_ifos(self, ifos):
    """
    Add the list of interferometers 
    """
    self.__ifos = ifos
    self.add_var_opt('coinc-cut',ifos)

  def get_ifos(self):
    """
    Returns the ifos
    """
    return self.__ifos

  def set_slides(self, slides):
    """
    Add the number of time slides
    """
    self.__num_slides = slides 
    self.add_var_opt('num-slides',slides)

  def get_slides(self):
    """
    Returns the ifos
    """
    return self.__num_slides

  def set_ifo_tag( self, ifotag):
    """
    Set the ifotag
    """
    self.__ifotag = ifotag

  def get_ifo_tag( self ):
    """
    get the ifotag
    """
    return self.__ifotag

  def set_inj_file(self, file):
    """
    Sets the injection file
    """
    if file:
      self.__injection_file = file
      self.add_var_opt('injection-file', file)

  def get_inj_file(self, file):
    """
    Sets the injection file
    """
    return self.__injection_file

  def set_start(self, start):
    """
    Sets GPS start time
    """
    self.__start = start

  def get_start(self):
    """
    Returns GPS start time
    """
    return self.__start

  def set_end(self, end):
    """
    Sets GPS end time
    """
    self.__end = end

  def get_end(self):
    """
    Returns GPS end time
    """
    return self.__end

  def set_glob(self, file_glob):
    """
    Sets the glob name
    """
    self.add_var_opt('glob',file_glob)

  def set_input(self, input_file):
    """
    Sets the input file name
    """
    self.add_var_opt('input',input_file)

  def set_user_tag(self,usertag):
    """
    Set the usertag for a given job
    """
    self.__usertag = usertag
    self.add_var_opt('user-tag',usertag)

  def get_user_tag(self):
    """
    Returns the usertag of the job
    """
    return self.__usertag

  def set_output_tag(self):
    fname = "COIRE"
    if self.__num_slides: fname += "_SLIDE"
    if self.__injection_file:
      fname += "_" + self.__injection_file.split("-")[1]
      fname += "_FOUND"
    if self.__ifotag: fname += "_" + self.__ifotag
    if self.__usertag: fname += "_" + self.__usertag
    self.__output_tag = fname

  def get_output_tag(self):
    return self.__output_tag

  def get_output(self):
    """
    get the name of the output file
    """
    if not self.get_ifos():
      raise InspiralError, "ifos have not been set"

    self.set_output_tag()
    fname = self.__ifos + '-' + self.__output_tag

    if (self.__start and not self.__end) or \
           (self.__end and not self.__start):
      raise InspiralError, "If one of start and end is set, "\
            "both must be"

    if (self.__start):
      duration=self.__end - self.__start
      fname += "-" + str(self.__start) + "-" + str(duration)

    fname += ".xml"

    return fname

  def get_missed(self):
    """
    get the name of the missed file
    """
    if self.__injection_file:
      return self.get_output().replace("FOUND", "MISSED")
    else:
      return None

  def finalize(self):
    """
    set the output options
    """
    output = self.get_output()
 
    self.add_var_opt("output", output)
    self.add_var_opt("summary", output.replace("xml", "txt"))

    if self.__injection_file:
      self.add_var_opt('injection-file', self.__injection_file)
      self.add_var_opt('missed-injections', self.get_missed() )

class FrJoinNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A FrJoinNode runs an instance of lalapps_frjoin in a Condor DAG
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_frjoin.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)

  def set_output(self, outputName):
    """
    Set the output name of the frame file
    @param outputName: name of the injection file created
    """
    self.add_var_opt('output',outputName)
    self.__outputName = outputName
    
  def get_output(self):
    """
    Get the output name of the frame file
    """
    return self.__outputName



class CohBankNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A CohBankNode runs an instance of the coherent code in a Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_coherent_inspiral.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    try:
      self.__usertag = job.get_config('cohbank','user-tag')
    except:
      self.__usertag = job.get_config('pipeline','user-tag')
    self.__bank = None
    self.__ifos = None
    try:
      self.__zip_output = job.get_config('cohbank','write-compress')
      self.__zip_output = True
    except:
      self.__zip_output = False

  def set_user_tag(self,usertag):
    """
    Set the usertag for a given job
    """
    self.__usertag = usertag
    self.add_var_opt('user-tag',usertag)

  def get_user_tag(self):
    """
    Returns the usertag of the job
    """
    return self.__usertag     
    
  def set_bank(self,bank):
    self.add_var_opt('bank-file', bank)
    self.add_input_file(bank)
    self.__bank = bank

  def get_bank(self):
    return self.__bank

  def set_ifos(self,ifos):
    self.add_var_opt('ifos', ifos)
    self.__ifos = ifos
   
  def get_ifos(self):
    return self.__ifos
    
  def get_output(self):
    """
    Returns the file name of output from the coherent bank. 
    """
    
    if not self.get_ifos():
      raise InspiralError, "Ifos have not been set"
    
    basename = self.get_ifos() + '-COHBANK'

    if self.__usertag:
      basename += '_' + self.__usertag

    filename = basename + '-' + str(self.get_start()) + '-' + \
      str(self.get_end() - self.get_start()) + '.xml'

    if self.__zip_output:
      filename += '.gz'

    self.add_output_file(filename)

    return filename    


class ChiaNode(pipeline.CondorDAGNode,pipeline.AnalysisNode):
  """
  A ChiaNode runs an instance of the coherent_inspiral code in a Condor
  DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_coherent_inspiral.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    try:
      self.__zip_output = job.get_config('chia','write-compress')
      self.__zip_output = True
    except:
      self.__zip_output = False

  def set_bank(self,bank):
    self.add_var_opt('bank-file', bank)
    self.add_input_file(bank)

    
##############################################################################
# some functions to make life easier later

def overlap_test(interval1, interval2, slide_sec=0):
  """
  Test whether the two intervals could possibly overlap with one of them being
  slid by a maximum time of slide_sec.  Perform three tests:
  1)  Does the start of interval 1 lie within interval 2's range (with the 
    start decremented by slide_sec and the end incremented by slide_sec)
  2)  Does the end of interval 1 lie within interval 2's range (with the start 
    decremented by slide_sec and the end incremented by slide_sec)
  3)  Does interval 1 completely cover (the extended) interval 2, 
    ie is interval_1 start before (interval 2 start - slide_sec) AND 
    interval 1 end after (interval 2 end + slide_sec)
  If any of the above conditions are satisfied then return 1, else 0.
  """
  if ( 
    interval1.start() >= interval2.start() - slide_sec
    and interval1.start() <= interval2.end() + slide_sec
    ) or (
    interval1.end() >= interval2.start() - slide_sec 
    and interval1.end() <= interval2.end() + slide_sec
    ) or (
    interval1.start() <= interval2.start() - slide_sec    
    and interval1.end() >= interval2.end() + slide_sec ):
    return 1
  else:
    return 0



