"""
Classes needed for the cosmic string analysis pipeline.
"""

__author__ = 'Xavier Siemens<siemens@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'

import string,sys,os
from glue import pipeline
from pylab import *
import operator
from math import *
import numpy
import time

class StringError(Exception):
  def __init__(self, args=None):
    self.args = args


class DataFindJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A LSCdataFind job used by the string pipeline. The static options are
  read from the section [datafind] in the ini file. The stdout from
  LSCdataFind contains the paths to the frame files and is directed to a file
  in the cache directory named by site and GPS start and end times. The stderr
  is directed to the logs directory. The job always runs in the scheduler
  universe. The path to the executable is determined from the ini file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','datafind')
    self.__universe = 'scheduler'
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)

    for sec in ['datafind']:
      self.add_ini_opts(cp,sec)

    self.add_condor_cmd('environment',
      """LD_LIBRARY_PATH=$ENV(LD_LIBRARY_PATH);PYTHONPATH=$ENV(PYTHONPATH)""" )

    self.set_stderr_file('logs/datafind-$(macroinstrument)-$(macrostart)-$(macroend)-$(cluster)-$(process).err')
    self.set_stdout_file('cache/$(macroinstrument)-$(macrostart)-$(macroend).cache')
    self.set_sub_file('datafind.sub')


class StrainJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_StringSearch job used by the string pipeline. The static options
  are read from the section in the ini file. The
  stdout and stderr from the job are directed to the logs directory. The job
  runs in the universe specified in the ini file. The path to the executable
  is determined from the ini file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','strain')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)

    for sec in ['strain']:
      self.add_ini_opts(cp,sec)

    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")

    self.set_stdout_file('logs/strain-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/strain-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('strain.sub')

class NoiseJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_StringSearch job used by the string pipeline. The static options
  are read from the section in the ini file. The
  stdout and stderr from the job are directed to the logs directory. The job
  runs in the universe specified in the ini file. The path to the executable
  is determined from the ini file.
  """
  def __init__(self,cp,dax=False):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','noise')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp,dax)

    for sec in ['noisecomp']:
      self.add_ini_opts(cp,sec)

    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")

    self.set_stdout_file('logs/noisecomp-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/noisecomp-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('noisecomp.sub')

class EpochData(object):
  """
  Holds calibration data epochs
  """
  def __init__(self,cp,opts):
    self.file = open(cp.get('epochs','epochs_data'),'r')
    self.epoch_data = []
    for line in self.file:
      if line.strip()[0] != '#':
        self.epoch_data.append(line.split())
    self.epoch_data.sort()

  def epoch_segs(self):
    tmpList = []
    for line in self.epoch_data:
      tmpList.append(tuple(map(int,(line[0:4]))))
    return tmpList

class DataFindNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A DataFindNode runs an instance of datafind in a Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of LSCdataFind.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__start = 0
    self.__end = 0
    self.__instrument = None
    self.__output = None

  def __set_output(self):
    """
    Private method to set the file to write the cache to. Automaticaly set
    once the ifo, start and end times have been set.
    """
    if self.__start and self.__end and self.__instrument:
      self.__output = 'cache/' + self.__instrument + '-' + str(self.__start)
      self.__output = self.__output + '-' + str(self.__end) + '.cache'

  def set_start(self,time):
    """
    Set the start time of the datafind query.
    time = GPS start time of query.
    """
    self.add_var_opt('start', time)
    self.__start = time
    self.__set_output()

  def set_end(self,time):
    """
    Set the end time of the datafind query.
    time = GPS end time of query.
    """
    self.add_var_opt('end', time)
    self.__end = time
    self.__set_output()

  def set_ifo(self,ifo):
    """
    Set the IFO to retrieve data for. Since the data from both Hanford
    interferometers is stored in the same frame file, this takes the first
    letter of the IFO (e.g. L or H) and passes it to the --instrument option
    of LSCdataFind.
    ifo = IFO to obtain data for.
    """
    self.add_var_opt('instrument',ifo[0])
    self.__instrument = ifo[0]
    self.__set_output()

  def get_output(self):
    """
    Return the output file, i.e. the file containing the frame cache data.
    """
    return self.__output



class StrainNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A RingNode runs an instance of the ring code in a Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_StringSearch.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)


class NoiseNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A RingNode runs an instance of the ring code in a Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_StringSearch.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)

  def add_cache(self,file,frame):
    if isinstance( file, str ):
      # the name of a lal cache file created by a datafind node
      self.add_var_opt(frame, file)
      self.add_input_file(file)
    else:
      cacheFileName = 'cache/'+self.get_name()+'-'+frame+'-'+str(self.get_start())+'-'+str(self.get_end())+'.cache'
      cacheFile = open(cacheFileName,'w')
      self.add_var_opt(frame, cacheFileName)
      self.add_input_file(cacheFileName)
      # check we have an LFN list
      from glue import LDRdataFindClient
      if isinstance( file, LDRdataFindClient.lfnlist ):
        #self.add_var_opt('glob-frame-data',' ')
        # only add the LFNs that actually overlap with this job
        for lfn in file:
          a, b, c, d = lfn.split('.')[0].split('-')
          t_start = int(c)
          t_end = int(c) + int(d)
          if ( t_start <= self.get_end() and t_end >= self.get_start() ):
            self.add_input_file(lfn)
            cacheFile.write(a+' '+b+' '+c+' '+d+' '+lfn+'\n')
        # set the frame type based on the LFNs returned by datafind
        #self.add_var_opt('frame-type',b)
      else:
        raise CondorDAGNodeError("Unknown LFN cache format")

# Convenience functions to cat together noise output files.
def open_noise_cat_file(dir):
  outfilename = dir + '/' + dir + '.cat'
  try:  outfile = open(outfilename,'w')
  except:
    sys.stderr.write('Could not open '+outfilename+' for writing')
    sys.exit(1)
  return outfile

def cat_noise_jobs(file,node):
  out = node.get_output_files()
  tmpfile = open(out[0],'r')
  tmpstr = ''.join( tmpfile.readlines() )
  try:  file.write(tmpstr)
  except: pass
  return

def bin(binVec,histVec,number):
  if number < binVec[0]:
    histVec[0] += 1
    return
  if number > binVec[-1]:
    histVec[-1] += 1
    return
  histVec[int(ceil(number/(binVec[-2]-binVec[-1])))+int(len(binVec)/2-1)] += 1
  #for i in range(1,len(binVec)):
    #if number > binVec[i-1] and number < binVec[i]:
      #histVec[i] += 1

def plot_systematics(filelist,cp,dir,epoch,dag,opts):
  flist = []
  for file in filelist:
    if os.path.split(file)[-1][:8] == "out-bin-":
      flist.append(file)
  freq = []
  freqfile = cp.get('noisecomp','freq-file')
  for line in open(freqfile,'r').readlines():
    freq.append(float(line.strip()))

  hfr1 = {}
  hfi1 = {}
  htr1 = {}
  hti1 = {}
  htrfr = {}
  htifi = {}
  htrfi = {}
  htifr = {}
  Ai = {}
  Ar = {}
  N = {}
  amp = {}
  amphf  = {}
  Nt = 0

  for f in freq:
    hfr1[f] = 0.0
    hfi1[f] = 0.0
    htr1[f] = 0.0
    hti1[f] = 0.0
    htrfr[f] = 0
    htifi[f] = 0
    htrfi[f] = 0
    htifr[f] = 0
    Ai[f] = 0.0
    Ar[f] = 0.0
    N[f] = 0.0
    amp[f] = 0.0
    amphf[f] = 0.0
  freqcnt = 0;

  print("\tfirst pass through systematics files...")
  for file in flist:
    try: input = open(file,'r')
    except:
      print("WARNING: file " + file + " doesn't exist")
      continue
    for line in input.readlines():
      tmp = line.split()
      if len(tmp) == 1:
        Nt += 1
        continue

      ampt = float(tmp[0])**2 + float(tmp[1])**2
      ampf = float(tmp[2])**2 + float(tmp[3])**2
      N[freq[freqcnt]] += 1.0
      amp[freq[freqcnt]] += ampt/2.0 + ampf/2.0
      amphf[freq[freqcnt]] += ampf
      htr1[freq[freqcnt]] += float(tmp[0])
      hti1[freq[freqcnt]] += float(tmp[1])
      hfr1[freq[freqcnt]] += float(tmp[2])
      hfi1[freq[freqcnt]] += float(tmp[3])
      htrfr[freq[freqcnt]] += float(tmp[0])*float(tmp[2])#/ampf
      htifi[freq[freqcnt]] += float(tmp[1])*float(tmp[3])#/ampf
      htrfi[freq[freqcnt]] += float(tmp[0])*float(tmp[3])#/ampf
      htifr[freq[freqcnt]] += float(tmp[1])*float(tmp[2])#/ampf
#      htrfr[freq[freqcnt]] += float(tmp[0])*float(tmp[2])/ampf
#      htifi[freq[freqcnt]] += float(tmp[1])*float(tmp[3])/ampf
#      htrfi[freq[freqcnt]] += float(tmp[0])*float(tmp[3])/ampf
#      htifr[freq[freqcnt]] += float(tmp[1])*float(tmp[2])/ampf



      freqcnt += 1
      if freqcnt >= len(freq): freqcnt = 0
    #if N[freq[freqcnt]] > 100: break
  #Actually make it the mean
  for f in freq:
    #print "...found " + str(N[f]) + " of " + str(Nt) + " for "  + str(f)
    htr1[f] /= N[f]
    hti1[f] /= N[f]
    hfr1[f] /= N[f]
    hfi1[f] /= N[f]
    #Ai[f] = (htifr[f]-htrfi[f])/N[f]
    #Ar[f] = (htrfr[f]+htifi[f])/N[f]
    Ai[f] = (htifr[f]-htrfi[f])/(amp[f])
    Ar[f] = (htrfr[f]+htifi[f])/(amp[f])
    amp[f] = sqrt(amp[f]/N[f])


  fname = "Ar_Ai_"+epoch[1]+"-"+epoch[2]+".txt"
  fl = open(fname,'w')
  fl.write("#freq h(t) re sys\th(t) im sys\th(t) mag sys\th(t) phase sys\n")
  mag = {}
  phase = {}
  for f in freq:
    mag[f] = sqrt(Ar[f]*Ar[f]+Ai[f]*Ai[f])
    phase[f] = atan2(Ai[f],Ar[f])*180.0/3.14159
    fl.write(str(f) + "\t"+str(Ar[f])+"\t"+str(Ai[f])+"\t"+str(mag[f])+"\t"+str(phase[f])+"\n")
  fl.close()


  realHistVecs = {}
  imagHistVecs = {}
  binVec = []
  nbins = int(ceil(float(N.values()[0])/50.0));
  for f in freq:
    realHistVecs[f] = zeros(2*nbins)
    imagHistVecs[f] = zeros(2*nbins)
  for b in range(-nbins,nbins): binVec.append(float(b)/nbins)


  xr1 = {}
  xi1 = {}
  xr2 = {}
  xi2 = {}
  xr3 = {}
  xi3 = {}
  xr4 = {}
  xi4 = {}
  N = {}

  for f in freq:
    xr1[f] = 0.0
    xi1[f] = 0.0
    xr2[f] = 0.0
    xi2[f] = 0.0
    xr3[f] = 0.0
    xi3[f] = 0.0
    xr4[f] = 0.0
    xi4[f] = 0.0
    N[f] = 0.0

  freqcnt = 0;

  print("\tsecond pass through systematics files...")
  #Compute the moments of the distribution
  for file in flist:
    try: input = open(file,'r')
    except:
      print("WARNING: file " + file + " doesn't exist")
      continue
    for line in input.readlines():
      tmp = line.split()
      if len(tmp) == 1: continue

      N[freq[freqcnt]] += 1.0
      htr = float(tmp[0])
      hti = float(tmp[1])
      hfr = float(tmp[2])
      hfi = float(tmp[3])
      xr = htr-Ar[freq[freqcnt]]*hfr + Ai[freq[freqcnt]]*hfi
      xi = hti-Ar[freq[freqcnt]]*hfi - Ai[freq[freqcnt]]*hfr
      bin(binVec,realHistVecs[freq[freqcnt]],xr/amp[freq[freqcnt]])
      bin(binVec,imagHistVecs[freq[freqcnt]],xi/amp[freq[freqcnt]])
      xr1[freq[freqcnt]] += xr
      xi1[freq[freqcnt]] += xi
      xr2[freq[freqcnt]] += xr*xr
      xi2[freq[freqcnt]] += xi*xi
      xr3[freq[freqcnt]] += xr*xr*xr
      xi3[freq[freqcnt]] += xi*xi*xi
      xr4[freq[freqcnt]] += xr*xr*xr*xr
      xi4[freq[freqcnt]] += xi*xi*xi*xi
      freqcnt += 1
      if freqcnt >= len(freq): freqcnt = 0
    #if N[freq[freqcnt]] > 100: break


  #Put them in units of the noise amplitude
  for f in freq:
    xr1[f] /= N[f]*amp[f]
    xi1[f] /= N[f]*amp[f]
    xr2[f] = sqrt(xr2[f]/N[f])/amp[f]
    xi2[f] = sqrt(xi2[f]/N[f])/amp[f]
    if xr3[f]:
      xr3[f] = pow(abs(xr3[f]/N[f]),1.0/3.0)/amp[f]*xr3[f]/abs(xr3[f])
    else: xr3[f] = 0.0
    if xi3[f]:
      xi3[f] = pow(abs(xi3[f]/N[f]),1.0/3.0)/amp[f]*xi3[f]/abs(xi3[f])
    else: xi3[f] = 0.0
    xr4[f] = pow(abs(xr4[f]/N[f]-3.0*pow(xr2[f]*amp[f],4)),1.0/4.0)/amp[f]*(xr4[f]/N[f]-3.0*pow(xr2[f]*amp[f],4))/abs(xr4[f]/N[f]-3.0*pow(xr2[f]*amp[f],4))
    xi4[f] = pow(abs(xi4[f]/N[f]-3.0*pow(xi2[f]*amp[f],4)),1.0/4.0)/amp[f]*(xi4[f]/N[f]-3.0*pow(xi2[f]*amp[f],4))/abs(xi4[f]/N[f]-3.0*pow(xi2[f]*amp[f],4))

  fname = "x1_x2_x3_x4_"+epoch[1]+"-"+epoch[2]+".txt"
  fl = open(fname,'w')
  fl.write("#freq \t xr \t xi \t xr^2 \t xi^2 \t xr^3 \t xi^3 \t xr^4 \t xi^4 \n")
  for f in freq:
    fl.write(str(f) + '\t' + str(xr1[f]) + '\t' + str(xi1[f]) + '\t' + str(xr2[f]) + '\t' + str(xi2[f]) + '\t' + str(xr3[f]) + '\t' + str(xi3[f]) + '\t' + str(xr4[f]) + '\t' + str(xi4[f]) + '\n')

  rootN = []
  rootNdeg = []
  vals = N.values()
  vals.sort()
  for f in vals:
    rootN.append(1/sqrt(f))
    rootNdeg.append(180/sqrt(f)/3.14159)

  fl.close()
  # Plot the results
  print("\tplotting...")
  # Plot the systematic in magnitude
  magfigname = "sys_mag"+epoch[1]+"-"+epoch[2]+".png"
  figure(1)
  keys = mag.keys()
  keys.sort()
  errorbar(keys,[mag[k] for k in keys],rootN)
  title('h(t) and h(f) magnitude systematics '+epoch[1]+"-"+epoch[2]+'\n')
  xlabel('Freq')
  ylabel('Mag')
  grid()
  savefig(dir + '/'+ magfigname)
  thumb = 'thumb-'+magfigname
  savefig(dir + '/'+ thumb,dpi=20)
  clf()
  #close()

 # Plot the systematic in phase
  phasefigname = "sys_phase"+epoch[1]+"-"+epoch[2]+".png"
  #figure(1)
  keys = phase.keys()
  keys.sort()
  errorbar(keys,[phase[k] for k in keys],rootNdeg)
  #errorbar(phase.keys(),phase.values(),rootNdeg)
  title('h(t) and h(f) phase systematics '+epoch[1]+"-"+epoch[2]+'\n')
  xlabel('Freq')
  ylabel('Phase (degrees)')
  grid()
  savefig(dir + '/'+ phasefigname)
  thumb = 'thumb-'+phasefigname
  savefig(dir + '/'+ thumb,dpi=20)
  clf()
  #close()

 # Plot the residual moments
  x1figname = "sys_x1_"+epoch[1]+"-"+epoch[2]+".png"
  #figure(1)
  keys = xr1.keys()
  keys.sort()
  plot(keys,[xr1[k] for k in keys])
  keys = xi1.keys()
  keys.sort()
  plot(keys,[xi1[k] for k in keys])
  #plot(xi1.keys(),xi1.values(),'r')
  legend(['real','imaginary'])
  title('mean '+epoch[1]+"-"+epoch[2]+'\n')
  xlabel('Freq')
  ylabel('mean')
  grid()
  savefig(dir + '/'+ x1figname)
  thumb = 'thumb-'+x1figname
  savefig(dir + '/'+ thumb,dpi=20)
  clf()
  #close()

 # Plot the residual moments
  x2figname = "sys_x2_"+epoch[1]+"-"+epoch[2]+".png"
  #figure(1)
  keys = xr2.keys()
  keys.sort()
  plot(keys,[xr2[k] for k in keys])
  keys = xi2.keys()
  keys.sort()
  plot(keys,[xi2[k] for k in keys])
#  plot(xr2.keys(),xr2.values())
#  plot(xi2.keys(),xi2.values(),'r')
  legend(['real','imaginary'])
  title('residual noise sqrt of second moment '+epoch[1]+"-"+epoch[2]+'\n')
  xlabel('Freq')
  ylabel('sigma')
  grid()
  savefig(dir + '/'+ x2figname)
  thumb = 'thumb-'+x2figname
  savefig(dir + '/'+ thumb,dpi=20)
  clf()
  #close()

 # Plot the residual moments
  x3figname = "sys_x3_"+epoch[1]+"-"+epoch[2]+".png"
  #figure(1)
  keys = xr3.keys()
  keys.sort()
  plot(keys,[xr3[k] for k in keys])
  keys = xi3.keys()
  keys.sort()
  plot(keys,[xi3[k] for k in keys])
#  plot(xr3.keys(),xr3.values())
#  plot(xi3.keys(),xi3.values(),'r')
  legend(['real','imaginary'])
  title('residual noise cube root of third moment '+epoch[1]+"-"+epoch[2]+'\n')
  xlabel('Freq')
  ylabel('cube root of skew')
  grid()
  savefig(dir + '/'+ x3figname)
  thumb = 'thumb-'+x3figname
  savefig(dir + '/'+ thumb,dpi=20)
  clf()
  #close()

 # Plot the residual moments
  x4figname = "sys_x4_"+epoch[1]+"-"+epoch[2]+".png"
  #figure(1)
  keys = xr4.keys()
  keys.sort()
  plot(keys,[xr4[k] for k in keys])
  keys = xi4.keys()
  keys.sort()
  plot(keys,[xi4[k] for k in keys])
 # plot(xr4.keys(),xr4.values())
 # plot(xi4.keys(),xi4.values(),'r')
  legend(['real','imaginary'])
  title('residual noise fourth root of excess kurtosis '+epoch[1]+"-"+epoch[2]+'\n')
  xlabel('Freq')
  ylabel('fourth root of excess kurtosis')
  grid()
  savefig(dir + '/'+ x4figname)
  thumb = 'thumb-'+x4figname
  savefig(dir + '/'+ thumb,dpi=20)
  clf()
  #close()

  pgname = dir + '/' + "sys_plots"+epoch[1]+"-"+epoch[2]+".html"
  page = open(pgname,'w')
  page.write('<h2>Plots of systematic errors between h(t) and h(f) '+epoch[1]+"-"+epoch[2]+'</h2><hr><br><br>\n')
  page.write('<h3>Systematics in magnitude and phase and moments of the residual distributions</h3><hr>\n')
  page.write('<img src='+magfigname+' width=600>\n')
  page.write('<img src='+phasefigname+' width=600>\n')
  page.write('<img src='+x1figname+' width=600>\n')
  page.write('<img src='+x2figname+' width=600>\n')
  page.write('<img src='+x3figname+' width=600>\n')
  page.write('<img src='+x4figname+' width=600><br><br>\n')
  page.write('<h3>Raw distribution of residual noise</h3><hr><br>\n')
  for f in freq:
    #time.sleep(10)
    print("plotting "+str(f))
    figname = "n_hist_"+str(f)+'_'+epoch[1]+"-"+epoch[2]+".png"
    #figure(1)
    plot(binVec,realHistVecs[f])
    plot(binVec,imagHistVecs[f],'r')
    legend(['real','imaginary'])
    title('residual noise distribution '+epoch[1]+"-"+epoch[2]+'\n'+'freq = '+str(f))
    ylabel('Number')
    xlabel('n / amp of h(t)')
    thumb = 'thumb-'+figname
    savefig(dir + '/'+ thumb,dpi=20)
    grid()
    savefig(dir + '/'+ figname)
    clf()
    #close()
    page.write('<a href='+figname+'><img src='+thumb+'></a>\n')

  page.close

def plot_noise_jobs(filelist,cp,dir,epoch,dag,qjob,opts):
  flist = []
  for file in filelist:
    if os.path.split(file)[-1][:8] != "out-bin-":
      flist.append(file)
  filelist = flist
  qfile = write_qscan_conf(dir,cp)
  filelist.sort()
  ifo = cp.get('pipeline','ifo')
  duration = cp.get('plot','duration')
  specList = []

  # identify the first file that exists and go from there...
  for ix in range(0,len(filelist)):
    try:
      f = open(filelist[ix],'r')
      break
    except:
      pass
  del(filelist[0:ix])

  input = open(filelist[0],'r')
  start = input.readlines()[0].split()[0]
  input.close()
  time = []
  freq = []
  fignames = []
  vetolist = []
  timeFreqTuple = []
  freqfile = cp.get('noisecomp','freq-file')
  for line in open(freqfile,'r').readlines():
    freq.append(float(line.strip()))
  #STOP = 0;
  for file in filelist:
    try: input = open(file,'r')
    except:
      print("WARNING: file " + file + " doesn't exist")
      continue
    #if STOP > 100: break
    #STOP+=1
    for line in input.readlines():
      tmp = line.split()
      if (float(tmp[0]) - float(start)) > float(duration):
        if opts.veto_list: vetolist.extend(make_veto_list(timeFreqTuple,cp))
        if opts.plot:
          fn,ft = plot_noise_spec(specList,cp,dir,dag,qjob,qfile,timeFreqTuple)
          fignames.append([fn,ft])
        start = float(tmp[0])
        specList = []
        timeFreqTuple = []
      specCol = []
      for ix in range(0,len(freq)):
        #PREVIOUS VERSION HAD 12 COLUMNS!
        #timeFreqTuple.append((float(tmp[0]),float(freq[ix]),float(tmp[1+2*ix])/float(tmp[2+2*ix]) ))
        timeFreqTuple.append((float(tmp[0]),float(freq[ix]),float(tmp[1+2*ix])/float(tmp[2+2*ix]) ))
        specCol.append(float(tmp[1+2*ix])/float(tmp[2+2*ix]))
      specList.append(specCol)
    input.close()
  fignames.sort(reverse=True)
  if opts.plot: write_spec_web_page(dir,ifo,cp,epoch,freq,fignames)
  if opts.veto_list: write_veto_list(ifo,epoch,vetolist,cp)

def make_veto_list(tftuple,cp):
  thresh = float(cp.get("veto","threshold"))
  veto = []
  tftuple.sort(key=operator.itemgetter(2))
  for tf in tftuple:
    if tf[2] < 1.0-thresh:
      veto.append(tf[0])
    else:
      break
  tftuple.sort(key=operator.itemgetter(2),reverse=True)
  for tf in tftuple:
    if tf[2] > 1.0+thresh:
      veto.append(tf[0])
    else:
      break
  return list(set(veto))

def write_veto_list(ifo,epoch,veto,cp):
  duration = cp.get("noisecomp","time")
  thresh = cp.get("veto","threshold")
  vfname = "veto_times_"+thresh+ifo+"_"+epoch[1]+"-"+epoch[2]+".txt"
  vfile = open(vfname,'w')
  vfile.write("# VETO list generated by lalapps_noise_comp_pipe\n")
  vfile.write("# flagged " + duration + " second times when h(t)/h(f) exceeded " + thresh + "\n")
  for ix in range(0,len(veto)):
    v = veto[ix]
    vfile.write(str(ix+1)+'\t'+str(v)+'\t'+str(v+float(duration))+'\t'+duration+'\n')
  vfile.close()

def write_spec_web_page(dir,ifo,cp,epoch,freq,fignames):
  web = open(dir+'/'+'index.html','w')
  web.write('<table>')
  web.write('<tr><td><b>IFO=</b>'+ifo+'</td></tr>\n')
  web.write('<tr><td><b>BAND=</b>'+cp.get('noisecomp','band')+'</td></tr>\n')
  web.write('<tr><td><b>TIME=</b>'+cp.get('noisecomp','time')+'</td></tr>\n')
  web.write('<tr><td><b>FCAL=</b>'+cp.get('noisecomp','fcal')+'</td></tr>\n')
  web.write('<tr><td><b>GAMMA-FUDGE-FACTOR=</b>'+cp.get('noisecomp','gamma-fudge-factor')+'</td></tr>\n')
  web.write('<tr><td><b>HOFT=</b>'+cp.get('noisecomp','hoft-channel')+'</td></tr>\n')
  web.write('<tr><td><b>ASQ=</b>'+cp.get('noisecomp','asq-channel')+'</td></tr>\n')
  web.write('<tr><td><b>EXC=</b>'+cp.get('noisecomp','exc-channel')+'</td></tr>\n')
  web.write('<tr><td><b>DARM=</b>'+cp.get('noisecomp','darm-channel')+'</td></tr>\n')
  web.write('<tr><td><b>DERR=</b>'+cp.get('noisecomp','derr-channel')+'</td></tr>\n')
  web.write('<tr><td><b>EPOCH=</b>'+epoch[0]+'</td></tr>\n')
  web.write('<tr><td><b>START=</b>'+epoch[1]+'</td></tr>\n')
  web.write('<tr><td><b>STOP=</b>'+epoch[2]+'</td></tr>\n')
  web.write('<tr><td><b>DURATION=</b>'+epoch[3]+'</td></tr>\n')
  web.write('<tr><td><b>OPEN LOOP GAIN (RE)=</b>'+epoch[4]+'</td></tr>\n')
  web.write('<tr><td><b>OPEN LOOP GAIN (IM)=</b>'+epoch[5]+'</td></tr>\n')
  web.write('<tr><td><b>SERVO (RE)=</b>'+epoch[6]+'</td></tr>\n')
  web.write('<tr><td><b>SERVO (IM)=</b>'+epoch[7]+'</td></tr>\n')
  web.write('<tr><td><b>WHITENER (RE)=</b>'+epoch[8]+'</td></tr>\n')
  web.write('<tr><td><b>WHITENER (IM)=</b>'+epoch[9]+'</td></tr>\n')
  web.write('<tr><td><b>OLG FILE=</b>'+epoch[11]+'</td></tr>\n')
  web.write('<tr><td><b>SENSING FILE=</b>'+epoch[12]+'</td></tr>\n')
  web.write('</table><br><br>')
  web.write('<h1>SPEC GRAMS, starting with the bad stuff between '+str(freq[0])+' and '+str(freq[-1])+' Hz</h1>')
  web.write('<table><tr>\n')
  cnter = 1
  for fig in fignames:
    if cnter == 6:
      web.write('</tr><tr>')
      cnter = 0
    web.write('<td><a href='+fig[0]+'><img src=thumb-'+fig[0]+'></a>\n')
    web.write('<br><a href='+str(fig[1])+'>'+str(fig[1])+'</a></td>\n')
    cnter += 1
  web.write('</tr></table>\n')
  web.close()

def plot_noise_spec(specList,cp,dir,dag,qjob,qfile,tftuple):
  fignames = []
  time = list(set(map(operator.itemgetter(0),tftuple)))
  time.sort()
  freq = list(set(map(operator.itemgetter(1),tftuple)))
  freq.sort()
  freq = array(freq,typecode='f')
  Time = array(time,typecode='d') - time[0]
  X,Y = meshgrid(Time,freq)
  start = str(time[0])
  end = str(time[-1])
  flat_specList = []
  tftuple.sort(key=operator.itemgetter(2))
  MIN = tftuple[0][2]
  MINTime = tftuple[0][0]
  tftuple.sort(key=operator.itemgetter(2),reverse=True)
  MAX = tftuple[0][2]
  MAXTime = tftuple[0][0]
  OUTLIER = [1-MIN, MAX-1]
  if (1-MIN) > (MAX-1):
    qscanTime = MINTime
  else:
    qscanTime = MAXTime
  dag.add_node(qscanNode(qjob,qscanTime,qfile,cp.get('pipeline','ifo'),dir,OUTLIER))
  figname = str(max(OUTLIER))+'-'+dir + '-' + start + '-' + end + '.png'
  A = array(specList,typecode='f')
  figure(1)
  pcolor(X,Y,A.transpose(),shading='flat',vmin=0.95,vmax=1.05)
  print("...plotting qscan for " + start)
  title('h(t) and h(f) power ratios per freq bin GPS '+start + '\n min = '+str(MIN) + ' max = '+str(MAX) )
  xlabel('Time')
  ylabel('Frequency')
  colorbar()
  savefig(dir + '/'+ figname)
  thumb = 'thumb-'+figname
  savefig(dir + '/'+ thumb,dpi=20)
  clf()
  #close()
  return figname,qscanTime

def write_qscan_conf(epoch,cp):
  duration = cp.get('noisecomp','time')
  qdur = pow(2,ceil(log(float(duration),2)))
  qdur = str(qdur)
  stdArgs = '''  sampleFrequency:\t\t2048
  searchTimeRange:\t\t'''
  stdArgs += qdur
  stdArgs += '''
  searchFrequencyRange:\t\t[32 Inf]
  searchQRange:\t\t\t[4 64]
  searchMaximumEnergyLoss:\t0.2
  whiteNoiseFalseRate:\t\t1e2
  searchWindowDuration:\t\t0.5
  plotTimeRanges:\t\t['''
  stdArgs += qdur
  stdArgs +=''']
  plotFrequencyRange:\t\t[]
  plotMaximumEnergyLoss:\t0.2
  plotNormalizedEnergyRange:\t[0 25.5]
'''

  L1Channels = ["'L1:LSC-DARM_ERR'", "'L1:LSC-DARM_CTRL_EXC_DAQ'", "'L1:LSC-ETMX_CAL_EXC_DAQ'", "'L1:LSC-ETMX_EXC_DAQ'", "'L1:LSC-ETMY_CAL_EXC_DAQ'", "'L1:LSC-ETMX_CAL'", "'L1:LSC-ETMY_CAL'", "'L1:LSC-AS1I_CORR_OUT_DAQ'", "'L1:LSC-AS1_I_DAQ'", "'L1:LSC-AS1_Q_DAQ'", "'L1:LSC-AS2I_CORR_OUT_DAQ'", "'L1:LSC-AS2_I_DAQ'", "'L1:LSC-AS2_Q_DAQ'", "'L1:LSC-AS3I_CORR_OUT_DAQ'", "'L1:LSC-AS3_I_DAQ'", "'L1:LSC-AS3_Q_DAQ'", "'L1:LSC-AS4I_CORR_OUT_DAQ'", "'L1:LSC-AS4_I_DAQ'", "'L1:LSC-AS4_Q_DAQ'", "'L1:LSC-AS_AC'", "'L1:LSC-AS_DC'", "'L1:LSC-AS_I'", "'L1:LSC-POB_I'", "'L1:LSC-POB_Q'", "'L1:LSC-REFL_AC'", "'L1:LSC-REFL_DC'", "'L1:LSC-REFL_I'", "'L1:LSC-REFL_Q'", "'L1:LSC-SPOB_I'", "'L1:ASC-BS_P'", "'L1:ASC-BS_Y'", "'L1:ASC-ETMX_P'", "'L1:ASC-ETMX_Y'", "'L1:ASC-ETMY_P'", "'L1:ASC-ETMY_Y'", "'L1:ASC-ITMX_P'", "'L1:ASC-ITMX_Y'", "'L1:ASC-ITMY_P'", "'L1:ASC-ITMY_Y'", "'L1:ASC-QPDX_DC'", "'L1:ASC-QPDX_P'", "'L1:ASC-QPDX_Y'", "'L1:ASC-QPDY_DC'", "'L1:ASC-QPDY_P'", "'L1:ASC-QPDY_Y'", "'L1:ASC-RM_P'", "'L1:ASC-RM_Y'", "'L1:ASC-WFS1_QP'", "'L1:ASC-WFS1_QY'", "'L1:ASC-WFS2_IP'", "'L1:ASC-WFS2_IY'", "'L1:ASC-WFS2_QP'", "'L1:ASC-WFS2_QY'", "'L1:ASC-WFS3_IP'", "'L1:ASC-WFS3_IY'", "'L1:ASC-WFS4_IP'", "'L1:ASC-WFS4_IY'", "'L1:IOO-MC_F'", "'L1:PSL-FSS_MIXERM_F'", "'L1:SUS-BS_OPLEV_PERROR'", "'L1:SUS-BS_OPLEV_YERROR'", "'L1:SUS-ETMX_OPLEV_PERROR'", "'L1:SUS-ETMX_OPLEV_YERROR'", "'L1:SUS-ETMX_SENSOR_SIDE'", "'L1:SUS-ETMY_OPLEV_PERROR'", "'L1:SUS-ETMY_OPLEV_YERROR'", "'L1:SUS-ETMY_SENSOR_SIDE'", "'L1:SUS-ITMX_OPLEV_PERROR'", "'L1:SUS-ITMX_OPLEV_YERROR'", "'L1:SUS-ITMY_OPLEV_PERROR'", "'L1:SUS-ITMY_OPLEV_YERROR'", "'L1:SUS-MMT3_OPLEV_PERROR'", "'L1:SUS-MMT3_OPLEV_YERROR'", "'L1:SUS-RM_OPLEV_PERROR'", "'L1:SUS-RM_OPLEV_YERROR'", "'L1:SEI-BS_RX'", "'L1:SEI-BS_RY'", "'L1:SEI-BS_RZ'", "'L1:SEI-BS_X'", "'L1:SEI-BS_Y'", "'L1:SEI-BS_Z'", "'L1:SEI-ETMX_RX'", "'L1:SEI-ETMX_RY'", "'L1:SEI-ETMX_RZ'", "'L1:SEI-ETMX_STS2_X'", "'L1:SEI-ETMX_X'", "'L1:SEI-ETMX_Y'", "'L1:SEI-ETMX_Z'", "'L1:SEI-ETMY_RX'", "'L1:SEI-ETMY_RY'", "'L1:SEI-ETMY_RZ'", "'L1:SEI-ETMY_STS2_Y'", "'L1:SEI-ETMY_X'", "'L1:SEI-ETMY_Y'", "'L1:SEI-ETMY_Z'", "'L1:SEI-ITMX_RX'", "'L1:SEI-ITMX_RY'", "'L1:SEI-ITMX_RZ'", "'L1:SEI-ITMX_X'", "'L1:SEI-ITMX_Y'", "'L1:SEI-ITMX_Z'", "'L1:SEI-ITMY_RX'", "'L1:SEI-ITMY_RY'", "'L1:SEI-ITMY_RZ'", "'L1:SEI-ITMY_X'", "'L1:SEI-ITMY_Y'", "'L1:SEI-ITMY_Z'", "'L1:SEI-LVEA_STS2_X'", "'L1:SEI-LVEA_STS2_Y'", "'L1:SEI-LVEA_STS2_Z'", "'L1:SEI-MC1_RX'", "'L1:SEI-MC1_RY'", "'L1:SEI-MC1_RZ'", "'L1:SEI-MC1_X'", "'L1:SEI-MC1_Y'", "'L1:SEI-MC1_Z'", "'L1:SEI-MC2_RX'", "'L1:SEI-MC2_RY'", "'L1:SEI-MC2_RZ'", "'L1:SEI-MC2_X'", "'L1:SEI-MC2_Y'", "'L1:SEI-MC2_Z'", "'L1:SEI-OUT_RX'", "'L1:SEI-OUT_RY'", "'L1:SEI-OUT_RZ'", "'L1:SEI-OUT_X'", "'L1:SEI-OUT_Y'", "'L1:SEI-RM_RX'", "'L1:SEI-RM_RY'", "'L1:SEI-RM_RZ'", "'L1:SEI-RM_X'", "'L1:SEI-RM_Y'", "'L1:SEI-RM_Z'", "'L0:PEM-BSC1_ACCX'", "'L0:PEM-BSC1_ACCY'", "'L0:PEM-BSC1_ACCZ'", "'L0:PEM-BSC2_ACCX'", "'L0:PEM-BSC2_ACCY'", "'L0:PEM-BSC2_ACCZ'", "'L0:PEM-BSC3_ACCX'", "'L0:PEM-BSC3_ACCY'", "'L0:PEM-BSC3_ACCZ'", "'L0:PEM-BSC4_ACCX'", "'L0:PEM-BSC4_ACCY'", "'L0:PEM-BSC4_ACCZ'", "'L0:PEM-BSC4_MIC'", "'L0:PEM-BSC5_ACCX'", "'L0:PEM-BSC5_ACCY'", "'L0:PEM-BSC5_ACCZ'", "'L0:PEM-BSC5_MIC'", "'L0:PEM-COIL_MAGX'", "'L0:PEM-COIL_MAGZ'", "'L0:PEM-EX_BAYMIC'", "'L0:PEM-EX_MAGX'", "'L0:PEM-EX_MAGY'", "'L0:PEM-EX_MAGZ'", "'L0:PEM-EX_SEISX'", "'L0:PEM-EX_SEISY'", "'L0:PEM-EX_SEISZ'", "'L0:PEM-EX_V1'", "'L0:PEM-EY_BAYMIC'", "'L0:PEM-EY_MAGX'", "'L0:PEM-EY_MAGY'", "'L0:PEM-EY_MAGZ'", "'L0:PEM-EY_SEISX'", "'L0:PEM-EY_SEISY'", "'L0:PEM-EY_SEISZ'", "'L0:PEM-EY_V1'", "'L0:PEM-HAM1_ACCX'", "'L0:PEM-HAM1_ACCZ'", "'L0:PEM-HAM2_ACCX'", "'L0:PEM-HAM2_ACCZ'", "'L0:PEM-ISCT1_ACCX'", "'L0:PEM-ISCT1_ACCY'", "'L0:PEM-ISCT1_ACCZ'", "'L0:PEM-ISCT1_MIC'", "'L0:PEM-ISCT4_ACCX'", "'L0:PEM-ISCT4_ACCY'", "'L0:PEM-ISCT4_ACCZ'", "'L0:PEM-ISCT4_MIC'", "'L0:PEM-LVEA_BAYMIC'", "'L0:PEM-LVEA_MAGX'", "'L0:PEM-LVEA_MAGY'", "'L0:PEM-LVEA_MAGZ'", "'L0:PEM-LVEA_MIC'", "'L0:PEM-LVEA_SEISX'", "'L0:PEM-LVEA_SEISY'", "'L0:PEM-LVEA_SEISZ'", "'L0:PEM-LVEA_V1'", "'L0:PEM-RADIO_LVEA'"]


  H1Channels = ["'H1:LSC-DARM_ERR'", "'H1:LSC-DARM_CTRL'", "'H1:LSC-AS_Q'", "'H1:LSC-DARM_CTRL_EXC_DAQ'", "'H1:LSC-ETMX_CAL_EXC_DAQ'", "'H1:LSC-ETMX_EXC_DAQ'", "'H1:LSC-ETMY_CAL_EXC_DAQ'", "'H1:LSC-ETMX_CAL'", "'H1:LSC-ETMY_CAL'", "'H1:LSC-AS_AC'", "'H1:LSC-AS_DC'", "'H1:LSC-AS_I'", "'H1:LSC-AS_Q_0FSR'", "'H1:LSC-AS_Q_1FSR'", "'H1:LSC-MC_L'", "'H1:LSC-MICH_CTRL'", "'H1:LSC-POB_I'", "'H1:LSC-POB_Q'", "'H1:LSC-POY_DC'", "'H1:LSC-PRC_CTRL'", "'H1:LSC-REFL_AC'", "'H1:LSC-REFL_DC'", "'H1:LSC-REFL_I'", "'H1:LSC-REFL_Q'", "'H1:LSC-SPOB_I'", "'H1:ASC-BS_P'", "'H1:ASC-BS_Y'", "'H1:ASC-ETMX_P'", "'H1:ASC-ETMX_Y'", "'H1:ASC-ETMY_P'", "'H1:ASC-ETMY_Y'", "'H1:ASC-ITMX_P'", "'H1:ASC-ITMX_Y'", "'H1:ASC-ITMY_P'", "'H1:ASC-ITMY_Y'", "'H1:ASC-QPDX_DC'", "'H1:ASC-QPDX_P'", "'H1:ASC-QPDX_Y'", "'H1:ASC-QPDY_DC'", "'H1:ASC-QPDY_P'", "'H1:ASC-QPDY_Y'", "'H1:ASC-RM_P'", "'H1:ASC-RM_Y'", "'H1:ASC-WFS1_QP'", "'H1:ASC-WFS1_QY'", "'H1:ASC-WFS2_IP'", "'H1:ASC-WFS2_IY'", "'H1:ASC-WFS2_QP'", "'H1:ASC-WFS2_QY'", "'H1:ASC-WFS3_IP'", "'H1:ASC-WFS3_IY'", "'H1:ASC-WFS4_IP'", "'H1:ASC-WFS4_IY'", "'H1:IOO-MC_F'", "'H1:PSL-FSS_MIXERM_F'", "'H1:SUS-BS_OPLEV_PERROR'", "'H1:SUS-BS_OPLEV_YERROR'", "'H1:SUS-BS_SENSOR_SIDE'", "'H1:SUS-ETMX_OPLEV_PERROR'", "'H1:SUS-ETMX_OPLEV_YERROR'", "'H1:SUS-ETMX_SENSOR_SIDE'", "'H1:SUS-ETMY_OPLEV_PERROR'", "'H1:SUS-ETMY_OPLEV_YERROR'", "'H1:SUS-ETMY_SENSOR_SIDE'", "'H1:SUS-FMX_OPLEV_PERROR'", "'H1:SUS-FMX_OPLEV_YERROR'", "'H1:SUS-FMX_SENSOR_SIDE'", "'H1:SUS-FMY_OPLEV_PERROR'", "'H1:SUS-FMY_OPLEV_YERROR'", "'H1:SUS-FMY_SENSOR_SIDE'", "'H1:SUS-ITMX_OPLEV_PERROR'", "'H1:SUS-ITMX_OPLEV_YERROR'", "'H1:SUS-ITMX_SENSOR_SIDE'", "'H1:SUS-ITMY_OPLEV_PERROR'", "'H1:SUS-ITMY_OPLEV_YERROR'", "'H1:SUS-ITMY_SENSOR_SIDE'", "'H1:SUS-MMT3_OPLEV_PERROR'", "'H1:SUS-MMT3_OPLEV_YERROR'", "'H1:SUS-RM_OPLEV_PERROR'", "'H1:SUS-RM_OPLEV_YERROR'", "'H1:SUS-RM_SENSOR_SIDE'", "'H1:TCS-ITMX_PD1AC'", "'H1:TCS-ITMX_PD2AC'", "'H1:TCS-ITMY_PD1AC'", "'H1:TCS-ITMY_PD2AC'"]

  H2Channels = ["'H2:LSC-DARM_ERR'", "'H2:LSC-DARM_CTRL'", "'H2:LSC-AS_Q'", "'H2:LSC-DARM_CTRL_EXC_DAQ'", "'H2:LSC-ETMX_CAL_EXC_DAQ'", "'H2:LSC-ETMX_EXC_DAQ'", "'H2:LSC-ETMY_CAL_EXC_DAQ'", "'H2:LSC-ETMX_CAL'", "'H2:LSC-ETMY_CAL'", "'H2:LSC-AS_AC'", "'H2:LSC-AS_DC'", "'H2:LSC-AS_I'", "'H2:LSC-AS_Q_0FSR'", "'H2:LSC-AS_Q_1FSR'", "'H2:LSC-MC_L'", "'H2:LSC-MICH_CTRL'", "'H2:LSC-POB_I'", "'H2:LSC-POB_Q'", "'H2:LSC-POY_DC'", "'H2:LSC-PRC_CTRL'", "'H2:LSC-REFL_AC'", "'H2:LSC-REFL_DC'", "'H2:LSC-REFL_I'", "'H2:LSC-REFL_Q'", "'H2:LSC-SPOB_I'", "'H2:ASC-BS_P'", "'H2:ASC-BS_Y'", "'H2:ASC-ETMX_P'", "'H2:ASC-ETMX_Y'", "'H2:ASC-ETMY_P'", "'H2:ASC-ETMY_Y'", "'H2:ASC-ITMX_P'", "'H2:ASC-ITMX_Y'", "'H2:ASC-ITMY_P'", "'H2:ASC-ITMY_Y'", "'H2:ASC-QPDX_DC'", "'H2:ASC-QPDX_P'", "'H2:ASC-QPDX_Y'", "'H2:ASC-QPDY_DC'", "'H2:ASC-QPDY_P'", "'H2:ASC-QPDY_Y'", "'H2:ASC-RM_P'", "'H2:ASC-RM_Y'", "'H2:ASC-WFS1_QP'", "'H2:ASC-WFS1_QY'", "'H2:ASC-WFS2_IP'", "'H2:ASC-WFS2_IY'", "'H2:ASC-WFS2_QP'", "'H2:ASC-WFS2_QY'", "'H2:ASC-WFS3_IP'", "'H2:ASC-WFS3_IY'", "'H2:ASC-WFS4_IP'", "'H2:ASC-WFS4_IY'", "'H2:IOO-MC_F'", "'H2:PSL-FSS_MIXERM_F'", "'H2:SUS-BS_OPLEV_PERROR'", "'H2:SUS-BS_OPLEV_YERROR'", "'H2:SUS-BS_SENSOR_SIDE'", "'H2:SUS-ETMX_OPLEV_PERROR'", "'H2:SUS-ETMX_OPLEV_YERROR'", "'H2:SUS-ETMX_SENSOR_SIDE'", "'H2:SUS-ETMY_OPLEV_PERROR'", "'H2:SUS-ETMY_OPLEV_YERROR'", "'H2:SUS-ETMY_SENSOR_SIDE'", "'H2:SUS-FMX_OPLEV_PERROR'", "'H2:SUS-FMX_OPLEV_YERROR'", "'H2:SUS-FMX_SENSOR_SIDE'", "'H2:SUS-FMY_OPLEV_PERROR'", "'H2:SUS-FMY_OPLEV_YERROR'", "'H2:SUS-FMY_SENSOR_SIDE'", "'H2:SUS-ITMX_OPLEV_PERROR'", "'H2:SUS-ITMX_OPLEV_YERROR'", "'H2:SUS-ITMX_SENSOR_SIDE'", "'H2:SUS-ITMY_OPLEV_PERROR'", "'H2:SUS-ITMY_OPLEV_YERROR'", "'H2:SUS-ITMY_SENSOR_SIDE'", "'H2:SUS-MMT3_OPLEV_PERROR'", "'H2:SUS-MMT3_OPLEV_YERROR'", "'H2:SUS-RM_OPLEV_PERROR'", "'H2:SUS-RM_OPLEV_YERROR'", "'H2:SUS-RM_SENSOR_SIDE'", "'H2:TCS-ITMX_PD1AC'", "'H2:TCS-ITMX_PD2AC'", "'H2:TCS-ITMY_PD1AC'", "'H2:TCS-ITMY_PD2AC'"]


  H0Channels = ["'H0:PEM-BSC10_ACC1Y'", "'H0:PEM-BSC10_MAGX'", "'H0:PEM-BSC10_MAGY'", "'H0:PEM-BSC10_MAGZ'", "'H0:PEM-BSC10_MIC'", "'H0:PEM-BSC1_ACCY'", "'H0:PEM-BSC1_MAG1X'", "'H0:PEM-BSC1_MAG1Y'", "'H0:PEM-BSC1_MAG1Z'", "'H0:PEM-BSC2_ACCX'", "'H0:PEM-BSC2_ACCY'", "'H0:PEM-BSC3_ACCX'", "'H0:PEM-BSC4_ACCX'", "'H0:PEM-BSC4_ACCY'", "'H0:PEM-BSC5_ACCX'", "'H0:PEM-BSC5_MAGX'", "'H0:PEM-BSC5_MAGY'", "'H0:PEM-BSC5_MAGZ'", "'H0:PEM-BSC5_MIC'", "'H0:PEM-BSC6_ACCY'", "'H0:PEM-BSC6_MAGX'", "'H0:PEM-BSC6_MAGY'", "'H0:PEM-BSC6_MAGZ'", "'H0:PEM-BSC6_MIC'", "'H0:PEM-BSC7_ACCX'", "'H0:PEM-BSC7_MIC'", "'H0:PEM-BSC8_ACCY'", "'H0:PEM-BSC8_MIC'", "'H0:PEM-BSC9_ACC1X'", "'H0:PEM-BSC9_MAGX'", "'H0:PEM-BSC9_MAGY'", "'H0:PEM-BSC9_MAGZ'", "'H0:PEM-BSC9_MIC'", "'H0:PEM-COIL_MAGX'", "'H0:PEM-COIL_MAGZ'", "'H0:PEM-EX_SEISX'", "'H0:PEM-EX_SEISY'", "'H0:PEM-EX_SEISZ'", "'H0:PEM-EX_V1'", "'H0:PEM-EX_V2'", "'H0:PEM-EY_SEISX'", "'H0:PEM-EY_SEISY'", "'H0:PEM-EY_SEISZ'", "'H0:PEM-EY_V1'", "'H0:PEM-EY_V2'", "'H0:PEM-HAM1_ACCX'", "'H0:PEM-HAM1_ACCZ'", "'H0:PEM-HAM3_ACCX'", "'H0:PEM-HAM7_ACCX'", "'H0:PEM-HAM7_ACCZ'", "'H0:PEM-HAM9_ACCX'", "'H0:PEM-IOT1_MIC'", "'H0:PEM-IOT7_MIC'", "'H0:PEM-ISCT10_ACCX'", "'H0:PEM-ISCT10_ACCY'", "'H0:PEM-ISCT10_ACCZ'", "'H0:PEM-ISCT10_MIC'", "'H0:PEM-ISCT1_ACCX'", "'H0:PEM-ISCT1_ACCY'", "'H0:PEM-ISCT1_ACCZ'", "'H0:PEM-ISCT1_MIC'", "'H0:PEM-ISCT4_ACCX'", "'H0:PEM-ISCT4_ACCY'", "'H0:PEM-ISCT4_ACCZ'", "'H0:PEM-ISCT4_MIC'", "'H0:PEM-ISCT7_ACCX'", "'H0:PEM-ISCT7_ACCY'", "'H0:PEM-ISCT7_ACCZ'", "'H0:PEM-ISCT7_MIC'", "'H0:PEM-LVEA2_V1'", "'H0:PEM-LVEA2_V2'", "'H0:PEM-LVEA2_V3'", "'H0:PEM-LVEA_MAGX'", "'H0:PEM-LVEA_MAGY'", "'H0:PEM-LVEA_MAGZ'", "'H0:PEM-LVEA_MIC'", "'H0:PEM-LVEA_SEISX'", "'H0:PEM-LVEA_SEISY'", "'H0:PEM-LVEA_SEISZ'", "'H0:PEM-MX_SEISX'", "'H0:PEM-MX_SEISY'", "'H0:PEM-MX_SEISZ'", "'H0:PEM-MX_V1'", "'H0:PEM-MX_V2'", "'H0:PEM-MY_SEISX'", "'H0:PEM-MY_SEISY'", "'H0:PEM-MY_SEISZ'", "'H0:PEM-MY_V1'", "'H0:PEM-MY_V2'", "'H0:PEM-PSL1_ACCX'", "'H0:PEM-PSL1_ACCZ'", "'H0:PEM-PSL1_MIC'", "'H0:PEM-PSL2_ACCX'", "'H0:PEM-PSL2_ACCZ'", "'H0:PEM-PSL2_MIC'", "'H0:PEM-RADIO_CS_1'", "'H0:PEM-RADIO_CS_2'", "'H0:PEM-RADIO_LVEA'"]

  ifo = cp.get('pipeline','ifo')
  qfilename = str(ifo)+str(epoch)+qdur+'.conf'
  qfile = open(qfilename,'w')
  qfile.write('[Context,Context]')
  qfile.write('[Gravitational,Gravitational wave data]\n')
  qfile.write('{\n')
  qfile.write("  channelName:\t'"+cp.get('noisecomp','hoft-channel')+"'\n")
  qfile.write("  frameType:\t'"+cp.get('input','type-hoft')+"'\n")
  qfile.write(stdArgs)
  qfile.write('}\n')
  qfile.write('{\n')
  qfile.write("  channelName:\t'"+cp.get('noisecomp','derr-channel')+"'\n")
  qfile.write("  frameType:\t'"+cp.get('input','type-derr')+"'\n")
  qfile.write(stdArgs)
  qfile.write('}\n')
  qfile.close()
  qfilename = str(ifo)+str(epoch)+qdur+'.conf'
  qfile = open(qfilename+'FULL','w')
  qfile.write('[Context,Context]')
  qfile.write('[Gravitational,Gravitational wave data]\n')
  qfile.write('{\n')
  qfile.write("  channelName:\t'"+cp.get('noisecomp','hoft-channel')+"'\n")
  qfile.write("  frameType:\t'"+cp.get('input','type-hoft')+"'\n")
  qfile.write(stdArgs)
  qfile.write('}\n\n')
  if ifo == 'H1':
    for chan in H1Channels+H0Channels:
      qfile.write('{\n')
      qfile.write("  channelName:\t"+chan+"\n")
      qfile.write("  frameType:\t'RDS_R_L1'\n")
      qfile.write(stdArgs)
      qfile.write('}\n\n')
  if ifo == 'H2':
    for chan in H2Channels+H0Channels:
      qfile.write('{\n')
      qfile.write("  channelName:\t"+chan+"\n")
      qfile.write("  frameType:\t'RDS_R_L1'\n")
      qfile.write(stdArgs)
      qfile.write('}\n\n')
  if ifo == 'L1':
    for chan in L1Channels:
      qfile.write('{\n')
      qfile.write("  channelName:\t"+chan+"\n")
      qfile.write("  frameType:\t'RDS_R_L1'\n")
      qfile.write(stdArgs)
      qfile.write('}\n\n')
  qfile.close()

  return qfilename


class qscanJob(pipeline.CondorDAGJob):
  """
  A qscan job
  """
  def __init__(self, cp, tag_base='QSCAN'):
    """
    """
    self.__executable = string.strip(cp.get('condor','qscan'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
    self.add_condor_cmd('getenv','True')
    self.set_stdout_file('logs/qscan-$(macrochannelname)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/qscan-$(macrochannelname)-$(cluster)-$(process).err')
    self.set_sub_file('qscan.sub')



class qscanNode(pipeline.CondorDAGNode):
  """
  Runs an instance of a qscan job
  """
  def __init__(self,job,time,qfile,ifo,dir,OUTLIER):
    """
    job = A CondorDAGJob that can run an instance of qscan.
    """
    self.id = ifo + '-qscan-' + repr(time)

    pipeline.CondorDAGNode.__init__(self,job)
    self.add_var_arg(repr(time))
    if max(OUTLIER) < 0.10:
      # just look at darm and h(t) for puny outliers.
      self.add_file_arg(qfile)
    else:
      print(".....found 10% outlier running full qscan\n")
      # run the standard qscan on outliers greater than 10%
      self.add_file_arg(qfile+'FULL')
    self.add_var_arg('@default')
    self.add_var_arg(dir)

class MkdirJob(pipeline.CondorDAGJob):
  """
  Runs an instance of mkdir in a DAG/DAX. Useful for grid submission.
  """
  def __init__(self,log_dir, cp, dax=False):
    self.__executable = cp.get('condor','mkdir')
    self.__universe = 'local'
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
#    AnalysisJob.__init__(self,cp,dax)
    self.add_condor_cmd('getenv','True')
    self.set_stdout_file(os.path.join( log_dir, 'mkdir-$(cluster)-$(process).out') )
    self.set_stderr_file(os.path.join( log_dir, 'mkdir-$(cluster)-$(process).err') )
    self.set_sub_file('mkdir.sub')

class MkdirNode(pipeline.CondorDAGNode):
  """
  Runs an instance of mkdir in a DAG/DAX. Useful for grid submission.
  """
  def __init__(self,job,dir):
    """
    @param job: A CondorDAGJob that can run an instance of ligolw_add
    @param dir UNDOCUMENTED
    """
    pipeline.CondorDAGNode.__init__(self,job)
    self.set_name('mkdir_'+dir.replace(' ',''))
    self.add_var_arg(dir)

    for file in dir.split():
      self.add_output_file(file+'/'+'.log')
      try: os.mkdir(file)
      except: pass


