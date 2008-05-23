#!/usr/bin/python
    
__Id__ = "$Id$"
__author__ = "Thomas Cokelaer <Thomas.Cokelaer@astro.cf.ac.uk>"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]
__name__ = "plotinspmissed"
__title__ = "Found and Missed plots for triggers"
    


import sys
import os
import glob
from pylab import *

mtsun = 4.92549095e-6

def GetChirpMassEtaFromMasses(mass1, mass2):
  """
  Converts the individual masses into chirp mass and eta parameters
  
  @status: mature
  @param mass1: the first mass of the system
  @param mass2: the second mass of the system
  @return: the chirpMass and eta parameters
  """
  totalMass = mass1 + mass2
  eta = mass1 * mass2 / totalMass / totalMass
  chirpMass = totalMass * eta**(3./5.)
  return chirpMass,eta


def GetChirpMassEtaFromTaus(tau0, tau3, fl):
  """
  Converts the tau parameters intpo chirp mass and eta parameters. Requires the
  lower cut-off frequency.
   
  @status: mature
  @param tau0: the tau0 parameter in seconds
  @param tau3: the tau3 parameter in seconds
  @param fl: the lower cut off frequency in Hz
  @return the chirpMass and eta parameters.
  @type tau3: double
  @type tau0: double
  @type fl: double
  """
  m1, m2 = GetMassesFromTaus(tau0, tau3, fl)
  chirpMass, eta = GetChirpMassEtaFromMasses(m1, m2)
  return chirpMass,eta


def GetMassesFromTaus(tau0, tau3, fl):
  """
  Converts the tau parameters into components masses. Requires the lower 
  cut-off frequency.
  
  @status: mature
  @param tau0: the tau0 parameter in seconds
  @param tau3: the tau3 parameter in seconds
  @param fl: the lower cut off frequency in Hz
  @return the mass1 and mass2
  @type tau3: double
  @type tau0: double
  @type fl: double
  """
  M = 5./32. / pi / pi / fl * tau3 / tau0 / mtsun
  eta = 1./8.* (32./5. * pi * tau0 / tau3)**(2./3) / tau3 / fl
  m1 = M/2.* (1.-sqrt(1.-4.*eta) )
  m2 = M/2.* (1.-sqrt(1.-4.*eta) )
  return m1, m2


 #--------------------------------------------------------- ReadXMLFiles class
class ReadXMLFiles:
  """
  A class that allows to read a XML file and returns tables such as sngl_inspiral,
  bankefficiency,process_params as dictionaries
  """
  def __init__(self, opts):
    """
    """
    self.files = glob.glob(opts.glob)
    self.verbose = opts.verbose
    self.results = None
    self.bank = None
    self.params = None
    self.values = None
    self.fl = None

  def getvalue(self,param):      
    """
    param is a string corrsponding to the field "param" to be read within a
    table.
    """
    value = None
    
    for p,v in zip(self.params,self.values):      
      if param in p:
        value = v
    return str(value.strip("\""))
  
  # a method to read a table within an XML files
  def readXML(self, filename, table, type="real"):
    """
    tool to read an XML file and return the corresponding data 
    @param filename: the name of the input file.
    @param table: the table to read
    """
    file = open(filename, "r")
    line = file.readline()  
    
    # first, we search for the requested table
    if self.verbose : print "searching for the table "+table 
    while line and line.find(table) < 0:
      line = file.readline() 
    

    if self.verbose: print "rebuilding the dictionary"
    if not line:
      return None
    # create a dictionary using the parameters of the xml table as keys
    table_summ = {}
    line = file.readline()  
    sortDict = []
    while line.find("Column") >= 0:
      this = line.split("\"")
      key = this[1].split(":")[2]
      sortDict.append(key)
      table_summ[key] = []
      line = file.readline()  
    count = 0

    print table_summ
    # now parse the whole data and coutn the rows until "\Stream" is found 
    line = file.readline()  
    while line.find("Stream") < 0 :
      count = count + 1
      fields = line.split(",")
      for field,key in zip(fields,sortDict):
        try:
          table_summ[key].append(float(field))
        except:
          table_summ[key].append((field))
      print len(line)
      if len(line)==0: break
      line = file.readline()  
      
    # check that the file was properly closed with a </stream> tag
    if self.verbose: print "Found "+str(count )+" elements"
      
    #convert  
    if self.verbose : print "converting to arrays"
    for key in sortDict:
      try:
        table_summ[key] = array(table_summ[key])
      except:
        print 'skip array conversion for '+ key 
        pass
  
    file.close()
    
    if table_summ is None: 
        print 'Warning: '+ table +'table cannot be found in'+filename
    return table_summ


  def getTables(self):
    """
    Method to read bankefficiency,sngl_inspiral and process_params tables in a 
    file that matches a user glob.
    """    
    # right now only one file can be read.
    if len(self.files)>1:
        raise ValueError,"""More than 1 file match your --glob argument. Only 
    is required. Try again."""
    
    for file in self.files:
      # first we need to read the process_params, and extract Fl, the lower cut-
      # off frequency, which will be used later on  
      try:
        process = self.readXML(file, 'process_params')
        self.params = process['param']
        self.values = process['value']
        self.fl = float(self.getvalue('--fl'))
      except:
        self.params = None
        self.values = None
        raise ValueError, """ The XMl file must contain a process_params table 
    with the --fl option at least. It does not seem to be present in the file 
    provide."""

      
        
        
      try:
        results = self.readXML(file, 'bankefficiency')
        # some data products that are not in the original table.
        results['totalmass_sim'] = results['mass1_sim'] + results['mass2_sim']
        
        results['eta_sim'] = results['mass1_sim'] * results['mass2_sim'] / \
            results['totalmass_sim'] / \
            results['totalmass_sim']
        # note the dummy variable. Indeed, these functions returns 2 values and 
        # therefore we need 2 outputs. Otherwise, the function returns a tuple and
        # not an array.
        results['chirpmass'],dummy = GetChirpMassEtaFromTaus(\
                                   array(results['tau0']),\
                                   array(results['tau3']), self.fl)
        results['chirpmass_sim'],dummy = GetChirpMassEtaFromMasses(\
                                     array(results['mass1_sim']),\
                                     array(results['mass2_sim']))
        m1,m2 = GetMassesFromTaus(array(results['tau0']),\
                                  array(results['tau3']),self.fl)
        results['totalmass'] = m1+m2
        results['eta'] = m1*m2/(m1+m2)/(m1+m2)
        self.results = results
      except:
        self.results = None
        pass
        
      # reading the sngl_inspiral table
      try:
        self.bank = self.readXML(file, 'sngl_inspiral')
      except:
        self.bank = None
        pass
    
      if self.bank is None:
        print \
"""Warning: no sngl_inspiral table found. so,no information
         related to the template bank can be extracted"""
    
    return self.results, self.bank, self.params, self.values

class Plotting():
  def __init__(self, verbose=False):
    self.figure_num = 1
    self.options={}
    self.marker='bo'
    self.markersize=10
    self.alpha=0.6
    self.linewidth = 1
    self.hold = False
    self.title = None
    self.verbose = verbose

  def settitle(self, text=None):
    self.title = text
      
      
  def plot_histogram_and_fit(self,data, nbins=20,fit=True):
    """
    """
    fig = None
    handle = None
    if self.hold is False:
      if self.verbose is True:
        print 'Plotting '+ str(self.figure_num) +' in progress...'
      figure(self.figure_num)
      self.figure_num += 1
      
    try:
      n, bins, patches = hist(data, nbins, normed=1)
      if fit is True:
        y = normpdf( bins, mean(data), 1./max(n))
        plotting.hold = True
        plotting.plot(bins, y, 'r--', linewidth=2)
        plotting.hold = False
        if self.title is not None:
          title(self.title) 
      gca().grid(True)
    except: print """Error while creating this plot"""
    return fig,handle

  def scatter(self, xdata, ydata, zdata, markersize=None, alpha=None):
    """
    @param xdata: an x vector
    @param ydata: a y vector
    @param zdata; a z vector
    @param markersize: size of the marker
    @type markersize: integer
    @param alpha: the transparency
    @type alpha: float
    """
    fig = None
    handle = None
    if self.hold is False:
      if self.verbose is True:
        print 'Plotting '+ str(self.figure_num) +' in progress...'
    
      figure(self.figure_num)
      self.figure_num += 1
    if markersize==None:
      markersize = self.markersize
    if alpha==None:
      alpha = self.alpha
    try:
      scatter(xdata,ydata, c=zdata,s=markersize, alpha=alpha)
      xlim(min(xdata), max(xdata))
      ylim(min(ydata), max(ydata))
      gca().grid(True)
      colorbar()
      if self.title is not None:
          title(self.title) 
    except: print """Error while creating this plot"""
    return fig,handle
    
  def plot(self, xdata,ydata, 
           marker=None, markersize=None, alpha=None, linewidth=None):
    """
    simple call to plot function. 
    Iterate the figure number
    
    @param xdata: an x vector
    @param ydata: a y vector
    @param marker: symbol and color of the data 
    @type marker: combinaison such as 'ro' for red circles.
    @param markersize: size of the marker
    @type markersize: integer
    @param alpha: the transparency
    @type alpha: float
    @param linewidth: the transparency
    @type linewidth: float
    """
    fig = None
    handle = None
    if self.hold is False:
      if self.verbose is True:
        print 'Plotting '+ str(self.figure_num) +' in progress...'
 
      figure(self.figure_num)
      self.figure_num += 1
      
    if marker==None:
      marker = self.marker
    if alpha==None:
      alpha = self.alpha
    if markersize==None:
      markersize = self.markersize
    if linewidth==None:
      linewidth = self.linewidth
      
    try:  
      plot(xdata,ydata, marker, \
           markersize=markersize,alpha=alpha,linewidth=linewidth)
      if self.title is not None:
          title(self.title) 
    except: print """Error while creating this plot"""
    
    gca().grid(True)
    return fig,handle


  def vectors2image(self, xdata,ydata, N=50):
    """
    simple call to plot function. 
    Iterate the figure number
    
    @param xdata: an x vector
    @param ydata: a y vector
    @param marker: symbol and color of the data 
    @type marker: combinaison such as 'ro' for red circles.
    @param markersize: size of the marker
    @type markersize: integer
    @param alpha: the transparency
    @type alpha: float
    @param linewidth: the transparency
    @type linewidth: float
    """
    fig = None
    handle = None
    if self.hold is False:
      if self.verbose is True:
        print 'Plotting '+ str(self.figure_num) +' in progress...'
 
      fig = figure(self.figure_num)
      self.figure_num += 1
      
      
    xmin = floor(min(xdata)*100)/100
    xmax = ceil(max(xdata)*100)/100    
    
    ymin = floor(min(ydata)*100)/100
    ymax = ceil(max(ydata)*100)/100
    
    dx = (xmax-xmin) / (N+1.)
    dy = (ymax-ymin) / (N+1.)
    
    x = arange(xmin, xmax, dx)
    y = arange(ymin, ymax, dy)
    X,Y = meshgrid(x,y)
    Z = X.copy()
    Z = Z * 0

    for s,M in zip(xdata,ydata):
      i = N-(int)(floor((M-ymin)/dy))
      j = (int)(floor((s-xmin)/dx))
      Z[i,j] += 1
      
    imshow(Z)
    colorbar()
    locs = arange(0,N+0.01,N/5.)
    
    
    newlabels = (arange(xmin,xmax,(xmax-xmin)/(len(locs))))
    Dx = (xmax-xmin)/(len(locs)-1)
    test = []
    for i in range(0,len(locs)):
      this = xmin+i*Dx
      test.append(str(int(this)))
    xticks(locs,tuple(test))

    newlabels = (arange(ymin,ymax,(ymax-ymin)/(len(locs))))
    Dy = (ymax-ymin)/(len(locs)-1)
    test = []
    for i in range(0,len(locs)):
      this = ymin+i*Dy
      this = floor(this*1000)/1000
      test.append(str(this))
    
    yticks(locs,tuple(test))

    
    
    if self.title is not None:
      title(self.title) 
#   except: print """Error while creating this plot"""
    
    gca().grid(True)
    return fig,handle
