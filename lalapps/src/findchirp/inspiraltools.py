#!/usr/bin/python
    
__Id__ = "$Id$"
__author__ = "Thomas Cokelaer <Thomas.Cokelaer@astro.cf.ac.uk>"
__version__ = "$Revision$"
__date__ = "$Date$"
__name__ = "inspiraltools"
__title__ = "Found and Missed plots for triggers"
    

import sys
import os
import glob
from pylab import *
# for the delaunay and griddata
try:
    from scipy.sandbox import delaunay
except:
    print """"Warning: the delaunay package could not be imported. 
    Functionalities such as the contour or surf plots will not be available. 
    Fix your configuration and packages."""    
#import numpy.core.ma as ma
from numpy import *
from math import log10

mtsun = 4.92549095e-6

def GetChirpMassEtaFromMasses(mass1, mass2):
  """
  Converts the individual masses into chirp mass and eta parameters
  
  status: mature
  @param mass1: the first mass of the system
  @param mass2: the second mass of the system
  @return: the chirpMass and eta parameters
  @author: Thomas Cokelaer
  """
  totalMass = mass1 + mass2
  eta = mass1 * mass2 / totalMass / totalMass
  chirpMass = totalMass * eta**(3./5.)
  return chirpMass,eta


def GetChirpMassEtaFromTaus(tau0, tau3, fl):
  """
  Converts the tau parameters intpo chirp mass and eta parameters. Requires the
  lower cut-off frequency.
   
  status: mature
  @param tau0: the tau0 parameter in seconds
  @param tau3: the tau3 parameter in seconds
  @param fl: the lower cut off frequency in Hz
  @return the chirpMass and eta parameters.
  type tau3: double
  type tau0: double
  type fl: double
  @author: Thomas Cokelaer
  """
  m1, m2, M, eta, dummy = GetMassesFromTaus(tau0, tau3, fl)
  chirpMass, eta = GetChirpMassEtaFromMasses(m1, m2)
  return chirpMass,eta


def GetMassesFromTaus(tau0, tau3, fl):
  """
  Converts the tau parameters into components masses. Requires the lower 
  cut-off frequency.
  
  status: mature
  @param tau0: the tau0 parameter in seconds
  @param tau3: the tau3 parameter in seconds
  @param fl: the lower cut off frequency in Hz
  @return the mass1 and mass2 , total mass and eta
  type tau3: double
  type tau0: double
  type fl: double
  @author: Thomas Cokelaer
  """
  epsilon = 0.01
  
  M = 5./32. / pi / pi / fl * tau3 / tau0 / mtsun
  eta = 1./8.* (32./5. * pi * tau0 / tau3)**(2./3) / tau3 / fl
  chirpMass = M * eta**(3./5.)
  # we may want to round up the numerical errors because eta must be <=0.25
  # of course if we have 0.30, this is another problem. One may want to adjust 
  # the epsilonm which is set empirically for now.
  try: 
    size = len(eta)
    for this,i in zip(eta,range(0,size)):
      if (this<0.25+epsilon) & (this>0.25):
        eta[i] = 0.25
        
  except:
    size = len([eta])
    if (eta>0.25) & (eta<0.25+epsilon):
        eta = 0.25
  
  
  m1 = M/2.* (1.-sqrt(1.-4.*eta) )
  m2 = M/2.* (1.+sqrt(1.-4.*eta) )
  return m1, m2, M, eta, chirpMass

def GetTausFromMasses(m1,m2,fl):
  """
  Converts the mass parameters into components masses. Requires the lower 
  cut-off frequency.
  
  status: mature
  @param m1: the individual mass1
  @param m2: the individual mass2
  @param fl: the lower cut off frequency in Hz
  @return the tau parameters
  type m1: double
  type m2: double
  type fl: double
  @author: Thomas Cokelaer
  
  """
  M = (float)(m1 + m2)
  M2 = M * M
  eta = m1 * m2 / M2
  piFl = pi * fl
  
  
  M *= mtsun;
  tau0 = 5.0 / (256.0 * eta * M**(5./3.) * piFl**(8./3.))
  tau3 = pi / (8.0 * eta * M**(2./3.) * piFl**(5./3.))

  return tau0,tau3

 #--------------------------------------------------------- ReadXMLFiles class
class ReadXMLFiles:
  """
  read a XML file and returns tables (e.g., sngl_inspiral) as dictionaries
  
  @author: Thomas Cokelaer
  """
  def __init__(self, opts):
    """
    @param opts: field requireda are verbose, glob
    """
    self.files = glob.glob(opts.glob)    
    self.verbose = opts.verbose
    self.results = None
    self.bank = None
    self.params = None
    self.values = None
    self.fl = None    
    try:
      self.nlines_to_read = opts.n
    except:
      self.nlines_to_read = 1e16

  def getvalue(self, param):      
    """
    param is a string corrsponding to the field "param" to be read within a
    table.
    @author: Thomas Cokelaer
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
    @param type UNDOCUMENTED
    @author: Thomas Cokelaer
    """
        
    file = open(filename, "r")
    line = file.readline()  
    # first, we search for the requested table
    if self.verbose : print "searching for the table "+table 
    while line and line.find(table) < 0:
      line = file.readline() 
    

    if self.verbose: print "reading the data..."
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

    # now parse the whole data and coutn the rows until "\Stream" is found 
    line = file.readline()
    
    
    
    while line.find("Stream") < 0:
      count = count + 1
      fields = line.split(",")
      for field,key in zip(fields,sortDict):
        try:
          table_summ[key].append(float(field))
        except:
          table_summ[key].append((field))
      if len(line)==0: break
      # if the table is bank efficiency and the user provided a limited 
      # number of lines to be read, then we break
      if table=='bankefficiency':
        if count >= self.nlines_to_read: break
      
        
      line = file.readline()  
    
      
    # check that the file was properly closed with a </stream> tag
    if self.verbose: print "Found "+str(count-1)+" elements"
      
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


  def getTables(self, name=None):
    """
    Method to read bankefficiency,sngl_inspiral and process_params tables in a 
    file that matches a user glob.
    @author: Thomas Cokelaer
    """    
    # right now only one file can be read.
    if len(self.files)>1:
        print "Warning: More than 1 file match your --glob argument. use the last one. Fix me"
    
    if len(self.files)==0:
      print 'Error, no file to read. check the spelling'
      sys.exit(1)
    
    for file in self.files:
      if self.verbose is True:
          print 'Reading ' + file
      # first we need to read the process_params, and extract Fl, the lower cut-
      # off frequency, which will be used later on  
      if name is None:
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

        #try:
        results = self.readXML(file, 'bankefficiency')
          
          # some data products that are not in the original table.
        results['totalmass_sim'] = results['mass1_sim'] + results['mass2_sim']        
        results['eta_sim'] = results['mass1_sim'] * results['mass2_sim'] / \
          results['totalmass_sim'] / results['totalmass_sim']
          # note the dummy variable. Indeed, these functions returns 2 values and 
          # therefore we need 2 outputs. Otherwise, the function returns a tuple and
          # not an array.
        results['chirpmass'],dummy = GetChirpMassEtaFromTaus(\
                                   array(results['tau0']),\
                                   array(results['tau3']), self.fl)
        results['chirpmass_sim'],dummy = GetChirpMassEtaFromMasses(\
                                     array(results['mass1_sim']),\
                                     array(results['mass2_sim']))
        m1,m2,M,eta, chirpMass = GetMassesFromTaus(array(results['tau0']),\
                                  array(results['tau3']),self.fl)        
        results['totalmass'] = M
        results['eta'] = eta
        results['mass1'] = m1
        results['mass2'] = m2
        results['chirpmass'] = chirpMass
        # it is better to get chirpmass from M and eta instead of tau0/tau3
        
        self.results = results
        
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
    
    if name is not None:
      try:
        self.results = self.readXML(file, name)
      except:
        self.results = None
        pass
    
    return self.results, self.bank, self.params, self.values


class Plotting:
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
    
  def setfigure(self):
    fig = None
    if self.hold is False:
      if self.verbose is True:
        print 'Plotting '+ str(self.figure_num) +' in progress...',
      fig = figure(self.figure_num)
      self.figure_num += 1
    return fig
   
      
  def plot_histogram_and_fit(self,data, nbins=20,fit=True):
    """
    @param data: the data to look at
    @param nbins: number of bins for the histogram
    @param fit: create a fit on top of the histogram
    """
    fig = self.setfigure()
    handle = None
      
    try:
      n, bins, patches = hist(data, nbins, normed=1)
    except:
      print """Error inside histogram (hist function)""" 
        
    if fit is True:
      try:
        y = normpdf( bins, mean(data), 1./max(n))
        plotting.hold = True        
        try:
          plotting.plot(bins, y, 'r--', linewidth=2)          
          plotting.hold = False
        except:pass
        
      except:
        print """Error inside histogram (normpdf function)""" 
    if self.title is not None:
      title(self.title) 
    
    gca().grid(True)
    
    return fig,handle

  def scatter(self, xdata, ydata, zdata, markersize=10, alpha=1,\
              vmin=None,vmax=None,type='normal'):
    """
    @param xdata: an x vector
    @param ydata: a y vector
    @param zdata; a z vector
    @param markersize: size of the marker
    type markersize: integer
    @param alpha: the transparency
    type alpha: float
    @param vmin UNDOCUMENTED
    @param vmax UNDOCUMENTED
    @param type UNDOCUMENTED
    @author: Thomas Cokelaer
    """
    fig = self.setfigure()
    handle = None
    
    if markersize==None:
      markersize = self.markersize
    if alpha==None:
      alpha = self.alpha
    
    try:
      if type=='log':
        scatter(xdata,ydata, c=log10(zdata),s=markersize, alpha=alpha,vmin=vmin,vmax=vmax)
        colorbar()
      else:  
        scatter(xdata, ydata, c=zdata,s=markersize, alpha=alpha,vmin=vmin,vmax=vmax)
        colorbar()
        
      xlim(min(xdata), max(xdata))
      ylim(min(ydata), max(ydata))
      
      gca().grid(True)
      
      if self.title is not None:
          title(self.title) 
    except: 
      print """Error inside scatter plots"""      
      
    return fig,handle
    
  def plot(self, xdata,ydata, 
           marker=None, markersize=None, alpha=None, linewidth=None):
    """
    simple call to plot function. 
    Iterate the figure number
    
    @param xdata: an x vector
    @param ydata: a y vector
    @param marker: symbol and color of the data 
    type marker: combinaison such as 'ro' for red circles.
    @param markersize: size of the marker
    type markersize: integer
    @param alpha: the transparency
    type alpha: float
    @param linewidth: the transparency
    type linewidth: float
    @author: Thomas Cokelaer
    """
    
    fig = self.setfigure()
    handle = None
         
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
    except: print """Error inside plot"""
    
    gca().grid(True)
    return fig,handle

  def griddata(self,xdata,ydata,zdata,xbin=20,ybin=20):
    
    if self.verbose: print 'Entering griddata...',
    
    xmin = min(xdata)
    xmax = max(xdata)
    xstep = (xmax - xmin)/xbin

    ymin = min(ydata)
    ymax = max(ydata)
    ystep = (ymax - ymin)/ybin
    print xmin,xmax,xstep,ymin,ymax,ystep
    xi, yi = mgrid[xmin:xmax:xstep, ymin:ymax:ystep]
    # triangulate data
    
    try:
        
      print 'here1'  
      tri = delaunay.Triangulation(xdata,ydata)
      print 'here2'
      interp = tri.nn_interpolator(zdata)    
      print 'here3'
      zi = interp(xi,yi)
      print 'here4'
    except:
      print 'Problem in griddata'

    if self.verbose: print 'Done'
    return xi,yi,zi

  def surf(self,xdata,ydata,zdata,xbin=20,ybin=20,vmin=0,vmax=1):
    """
    """
    fig = self.setfigure() 
    handle = None
   
    xi,yi,zi = self.griddata(xdata,ydata,zdata,xbin,ybin) 
    zim = ma.masked_where(isnan(zi),zi)
    figure(figsize=(8,8))
    p = pcolor(xi,yi,zim,shading='interp',cmap=cm.jet,
               vmin=vmin,vmax=vmax)
    gca().grid(True)
    #c = contour(xi,yi,zim,cmap=cm.jet)
    colorbar()
    return p

  def contourf(self,xdata,ydata,zdata,\
               xbin=20,ybin=20,xmin=0, xmax=0.4,vmin=0,vmax=1,\
               colormap=None,colorbar_title='',fontsize=16):
    """
     function to call contour plot
     the vmin option is harcoded for the time being using [0.5 0.8 0.9 0.95 and 1]
     colormap: a list to set the colormap
    """
    fig = self.setfigure() 
    handle = None
    
    xi,yi,zi = self.griddata(xdata,ydata,zdata,xbin,ybin)       
    
    zim = ma.masked_where(isnan(zi),zi)
    figure(figsize=(8,8))
    

    
    if colormap is None:
        map = [0.5, 0.8, 0.9, 0.95 ,0.965,1]
        c = contourf(xi,yi,zim,map)
    elif colormap=='cube':    
        map = [0.5**3., 0.8**3., 0.9**3.,0.95**3., 0.965**3.,1]
        c = contourf(xi,yi,zim,map)        
    else:
        c = contourf(xi,yi,zim,colormap)
            
    ylims = ylim()
    axis([xmin, xmax, ylims[0], ylims[1]])

    
    gca().grid(True)
    h = colorbar(format='%.3f')
    h.set_label(colorbar_title, fontsize=fontsize)
        
    return c  


  def vectors2image(self, xdata,ydata, N=50):
    """
    simple call to plot function. 
    Iterate the figure number
    
    @param xdata: an x vector
    @param ydata: a y vector
    @param N UNDOCUMENTED
    @author: Thomas Cokelaer
    """
    fig = self.setfigure()
         
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
