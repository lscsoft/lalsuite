"""
inspiral.py - classes needed for the inspiral analysis pipeline

$Id$

This script produced the necessary condor submit and dag files to run
the standalone inspiral code on LIGO data
"""

__author__ = 'Duncan Brown <duncan@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

import string
import exceptions
import pipeline



class InspiralError(exceptions.Exception):
  def __init__(self, args=None):
    self.args = args



class DataFindJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  def __init__(self,cp):
    pipeline.CondorDAGJob.__init__(self,'LALDataFind','scheduler')
    pipeline.AnalysisJob.__init__(self,cp)

    for sec in ['datafind']:
      self.add_ini_args(cp,sec)

    self.set_stdout_file('cache/$(frame-cache)')
    self.set_stderr_file('logs/datafind-$(instrument)-$(start)-$(end).err')
    self.set_sub_file('datafind.sub')



class TmpltBankJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  def __init__(self,cp):
    pipeline.CondorDAGJob.__init__(self,'lalapps_tmpltbank','vanilla')
    pipeline.AnalysisJob.__init__(self,cp)

    for sec in ['data','tmpltbank']:
      self.add_ini_args(cp,sec)

    self.set_stdout_file(
      'logs/tmpltbank-$(channel-name)-$(gps-start-time)-$(gps-end-time).out')
    self.set_stderr_file(
      'logs/tmpltbank-$(channel-name)-$(gps-start-time)-$(gps-end-time).err')
    self.set_sub_file('tmpltbank.sub')
    


class InspiralJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  def __init__(self,cp):
    pipeline.CondorDAGJob.__init__(self,'lalapps_inspiral','vanilla')
    pipeline.AnalysisJob.__init__(self,cp)

    for sec in ['data','inspiral']:
      self.add_ini_args(cp,sec)

    self.set_stdout_file(
      'logs/inspiral-$(channel-name)-$(gps-start-time)-$(gps-end-time).out')
    self.set_stderr_file(
      'logs/inspiral-$(channel-name)-$(gps-start-time)-$(gps-end-time).err')
    self.set_sub_file('inspiral.sub')
    


class TrigToTmpltJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  def __init__(self,cp):
    pipeline.CondorDAGJob.__init__(self,'lalapps_trigtotmplt','scheduler')
    pipeline.AnalysisJob.__init__(self,cp)
    
    for sec in ['trigtotmplt']:
      self.add_ini_args(cp,sec)
    
    self.set_stdout_file(
      'logs/trigtotmplt-$(channel-name)-$(gps-start-time)-$(gps-end-time).out')
    self.set_stderr_file(
      'logs/trigtotmplt-$(channel-name)-$(gps-start-time)-$(gps-end-time).err')
    self.set_sub_file('trigtotmplt.sub')



class IncaJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  def __init__(self,cp):
    pipeline.CondorDAGJob.__init__(self,'lalapps_inca','scheduler')
    pipeline.AnalysisJob.__init__(self,cp)
    
    for sec in ['inca']:
      self.add_ini_args(cp,sec)

    self.set_stdout_file('logs/inca-$(gps-start-time)-$(gps-end-time).out')
    self.set_stderr_file('logs/inca-$(gps-start-time)-$(gps-end-time).err')
    self.set_sub_file('inca.sub')



class DataFindNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  def __init__(self,job):
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__start = 0
    self.__end = 0
    self.__instrument = None
    self.__cache = None
   
  def set_start(self,time):
    self.add_var('start', time)
    self.__start = time
    self.set_cache()

  def set_end(self,time):
    self.add_var('end', time)
    self.__end = time
    self.set_cache()

  def set_ifo(self,ifo):
    self.add_var('instrument',ifo[0])
    self.__instrument = ifo[0]
    self.set_cache()

  def set_cache(self):
    if self.__start and self.__end and self.__instrument:
      self.__cache = self.__instrument + '-' + str(self.__start) 
      self.__cache = self.__cache + '-' + str(self.__end) + '.cache'
      self._CondorDAGNode__vars['frame-cache'] = self.__cache

  def get_output(self):
    return self.__cache



class TmpltBankNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  def __init__(self,job):
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)

  def get_output(self):
    if not self.get_start() or not self.get_end() or not self.get_ifo():
      raise InspiralError, "Start time, end time or ifo has not been set"
    bank = self.get_ifo() + '-TMPLTBANK-' + str(self.get_start())
    bank = bank + '-' + str(self.get_end() - self.get_start()) + '.xml'
    return bank



class InspiralNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  def __init__(self,job):
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__output = None

  def set_ifo(self,ifo):
    self._AnalysisNode__ifo = ifo
    self.add_var('channel-name', ifo + ':' + self.job().channel())
    self.add_var('calibration-cache', self.job().calibration(ifo))
    
  def set_bank(self,bank):
    self.add_var('bank-file', bank)
    out = self.get_ifo() + '-INSPIRAL-' + str(self.get_start()) + '-'
    self.__output = out + str(self.get_end() - self.get_start()) + '.xml'

  def get_output(self):
     return self.__output



class TrigToTmpltNode(pipeline.CondorDAGNode,pipeline.AnalysisNode):
  def __init__(self,job):
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)



class IncaNode(pipeline.CondorDAGNode,pipeline.AnalysisNode):
  def __init__(self,job):
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__input_a = []
    self.__input_b = []
    
  def add_input_a(self, file):
    self.__input_a.append(file)
    self.add_var('input-a', ' '.join(self.__input_a))

  def add_input_b(self, file):
    self.__input_b.append(file)
    self.add_var('input-b', ' '.join(self.__input_b))

  def set_output(self,string):
    self.__output = string + '-INCA-' + str(self.get_start()) + '-'
    self.__output = self.__output + str(self.get_end() - self.get_start())
    self.__output = self.__output + '.xml'

  def get_output(self):
    return self.__output
