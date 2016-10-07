#!/usr/bin/python

# combine_posteriors.py
# Copyright 2016 Jeroen Meidam <jmeidam@nikhef.nl>

import sys
py_version = sys.version_info[:2]
from os import path, makedirs
from pylab import figure,ioff,linspace
import matplotlib.patches as mpatches
from numpy import load,save,searchsorted,random,cumsum,savetxt,genfromtxt,array

from copy import deepcopy
from scipy.stats import gaussian_kde
ioff()

###########################################
#  Instructions
###########################################
"""
This script is used to plot combined posteriors for any number of sources
(depending on how crowded you are willing to make the plot).
The individual posteriors are also plotted.

The script can be used to create combined posteriors from single runs or batch runs
For single runs, the lalinference posterior files must have the following form:
  posterior_sourcename_paramname.dat
  e.g.: posterior_gw150914_dchi2.dat
For batch runs it must be as follows:
  posterior_sourcename_batchname.dat
  e.g.: posterior_gw150914_dchis.dat

To do a batch run, settings["singles"] must be set to False

When doing a batch run, violin plot compatible files are
automatically created.
These files are written to the settings["datadir"] path and have the form:
  posterior_samples_dchis.dat
"""

###########################################
#  Plotting parameters
###########################################

### As in the TGR paper:
from matplotlib import rc
import matplotlib
matplotlib.rc('text.latex', preamble = '\usepackage{txfonts}')
rc('text', usetex=True)
rc('font', family='serif')
rc('font', serif='times')
rc('mathtext', default='sf')
rc("lines", markeredgewidth=1)
rc("lines", linewidth=2)
rc('axes', labelsize=22) #24
rc("axes", linewidth=1.0) #2)
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('legend', fontsize=16) #10 in pltsettings.py
rc('xtick.major', pad=6) #8)
rc('ytick.major', pad=6) #8)
rc('xtick.minor', size=5) #8)
rc('ytick.minor', size=5) #8)

green = (0.17254901960784313, 0.6274509803921569, 0.17254901960784313)
gray = (0.4980392156862745, 0.4980392156862745, 0.4980392156862745)
lightgray = (0.8117647058823529, 0.8117647058823529, 0.8117647058823529)
blue = (0.12156862745098039, 0.4666666666666667, 0.7058823529411765)
red = (0.8392156862745098, 0.15294117647058825, 0.1568627450980392)
lightblue = (0.6823529411764706, 0.7803921568627451, 0.9098039215686274)
darkgrey = (0.34901960784313724, 0.34901960784313724, 0.34901960784313724)
purple = (0.5803921568627451, 0.403921568627451, 0.7411764705882353)
cyan = (0.09019607843137255, 0.7450980392156863, 0.8117647058823529)
orange = (1.0, 0.4980392156862745, 0.054901960784313725)
brown = (0.5490196078431373, 0.33725490196078434, 0.29411764705882354)


###########################################
#  User input
###########################################


settings = {"rescan_files":True,#false does not work for some reason
            "fill_all":False, #also fill the single posterior curves
            "combinedposterior_samples":20000, #number of samples to draw from the combined pdf
            "outpath":"/home/jeroen/Dropbox/combined_posteriors/", #plots will end up here
            "datadir":"/home/jeroen/Dropbox/combined_posteriors/", #data should be here
            "sources":["gw150914","gw151226"],
            "singles":True, #False if using parameters from batch runs
            "parameters":["dchi0","dchi1","dchi2","dchi3","dchi4","dchi5l","dchi6","dchi6l","dchi7",
                          "dsigma2","dsigma3","dsigma4",
                          "dbeta2","dbeta3",
                          "dalpha2","dalpha3","dalpha4"],
            "kde_cov_factor":0.25,
            "kde_points":1000}

xranges = {"dchi0":[-1.0,1.0],
          "dchi1":[-1.5,1.5],
          "dchi2":[-3,3],
          "dchi3":[-3,3],
          "dchi4":[-8,8],
          "dchi5l":[-3,3],
          "dchi6":[-5,5],
          "dchi6l":[-25,25],
          "dchi7":[-10,10],
          "dsigma2":[-8,8],
          "dsigma3":[-8,8],
          "dsigma4":[-25,25],
          "dbeta2":[-3,3],
          "dbeta3":[-3,3],
          "dalpha2":[-7,7],
          "dalpha3":[-7,7],
          "dalpha4":[-7,7]}

xranges_batch = {"dchi0":[-10.0,10.0],
          "dchi1":[-5,5],
          "dchi2":[-25,25],
          "dchi3":[-25,25],
          "dchi4":[-30,30],
          "dchi5l":[-30,30],
          "dchi6":[-30,30],
          "dchi6l":[-30,30],
          "dchi7":[-30,30],
          "dsigma2":[-40,40],
          "dsigma3":[-40,40],
          "dsigma4":[-40,40],
          "dbeta2":[-3,3],
          "dbeta3":[-3,3],
          "dalpha2":[-7,7],
          "dalpha3":[-7,7],
          "dalpha4":[-7,7]}

colors = {"gw150914":orange,
          "gw151226":blue,
          "combined":gray}

dchilatex = {"dchi0":r"$\delta\chi_0$","dchi1":r"$\delta\chi_1$","dchi2":r"$\delta\chi_2$",
             "dchi3":r"$\delta\chi_3$","dchi4":r"$\delta\chi_4$","dchi5l":r"$\delta\chi_5^l$",
             "dchi6":r"$\delta\chi_6$","dchi6l":r"$\delta\chi_6^l$","dchi7":r"$\delta\chi_7$",
             "dsigma2":r"$\delta\sigma_2$","dsigma3":r"$\delta\sigma_3$","dsigma4":r"$\delta\sigma_4$",
             "dbeta2":r"$\delta\beta_2$","dbeta3":r"$\delta\beta_3$",
             "dalpha2":r"$\delta\alpha_2$","dalpha3":r"$\delta\alpha_3$","dalpha4":r"$\delta\alpha_4$"}


###########################################
#  Classes and functions
###########################################

class posterior():
  def __init__(self,sourcename,paramname,posteriorfile,kde_min,kde_max,kde_points=1000,kde_cov_factor=0.25):
    """
    Container for the posterior samples and the kde
    """

    self.file = posteriorfile
    self.paramname = paramname
    self.sourcename = sourcename
    self.posterior=[]

    self.xaxis = linspace(kde_min,kde_max,kde_points)

    with open(posteriorfile,'r') as f:
      line = f.readline() #get header
      col = 0
      for p in line.strip().split():
        if p == paramname:
          break
        col+=1
      data = f.readlines()
      for line in data:
        self.posterior.append(float(line.split()[col]))

    #calculate kde
    density = gaussian_kde(self.posterior)
    density.covariance_factor = lambda : settings["kde_cov_factor"]
    density._compute_covariance()
    #self.density = density
    self.kde = density(self.xaxis)



class parameter():
  def __init__(self,paramname,sourcenames,settings,batchname=""):
    """
    Contains the combined posterior (kde) and all individual samples and kdes per source
    Provide batchname if parameters are part of a batch
    """

    self.name = paramname
    self.sourcenames = sourcenames
    self.sources = {}
    if settings["singles"]:
      self.kde_min = xranges[paramname][0]
      self.kde_max = xranges[paramname][1]
    else:
      self.kde_min = xranges_batch[paramname][0]
      self.kde_max = xranges_batch[paramname][1]

    print "working on",paramname

    for sourcename in self.sourcenames:
      if batchname == "":
        posteriorfile = path.join(settings["datadir"],"posterior_%s_%s.dat"%(sourcename,paramname))
      else:
        posteriorfile = path.join(settings["datadir"],"posterior_%s_%s.dat"%(sourcename,batchname))
      self.sources[sourcename] = posterior(sourcename,paramname,posteriorfile,
                                           self.kde_min,self.kde_max,
                                           kde_points=settings["kde_points"],kde_cov_factor=settings["kde_cov_factor"])

    self.combined_kde = deepcopy(self.sources[self.sourcenames[0]].kde)
    self.kde_xaxis = self.sources[self.sourcenames[0]].xaxis
    #Calculate combined
    for i in range(1,len(sourcenames)):
      newkde = self.sources[self.sourcenames[i]].kde
      for n in range(settings["kde_points"]):
        self.combined_kde[n] = self.combined_kde[n]*newkde[n]

    self.combined_kde = normalize_kde(self.combined_kde,self.kde_xaxis)

  def get_combined_kde(self):
    return self.combined_kde

  def get_source_kde(self,source):
    return self.sources[source].kde

  def get_kde_xaxis(self):
    return self.kde_xaxis

  def write_combined_posterior_samples(self,Nsamples):
    thiskde = self.combined_kde
    thisxaxis = self.kde_xaxis
    samples = generate_rand_from_pdf(thiskde, thisxaxis, Nsamples)
    if settings["singles"]:
      filename = path.join(settings["datadir"],"combined_posterior_%s.dat"%self.name)
    else:
      filename = path.join(settings["datadir"],"combined_posterior_batch_%s.dat"%self.name)
    savetxt(filename,samples)



def generate_rand_from_pdf(pdf, x_axis, N):
  cdf = cumsum(pdf)
  cdf = cdf / cdf[-1]
  values = random.rand(N)
  value_bins = searchsorted(cdf, values)
  random_from_cdf = x_axis[value_bins]
  return random_from_cdf

def batchname_from_parameter(param):
  if "dchi" in param:
    batchname = "dchis"
  elif "dbeta" in param:
    batchname = "dbetas"
  elif "dalpha" in param:
    batchname = "dalphas"
  elif "dsigma" in param:
    batchname = "dsigmas"
  else:
    print "could not set batchname"
    return ""
  return batchname

def ensure_dir(f):
  """
  Create folder if it does not exist
  """
  if not path.exists(f):
      makedirs(f)

def init(settings):
  """
  Check paths and files used in postprocessing
  """

  ##SANITY CHECKS


  ##FOLDER STRUCTURE:
  ensure_dir(settings["outpath"])
  ensure_dir(settings["datadir"])

  if settings["singles"]:
    datafile = path.join(settings["datadir"],"combined_posteriors_data.npy")
  else:
    datafile = path.join(settings["datadir"],"combined_posteriors_data_batches.npy")

  #get sources
  if (path.isfile(datafile) and not settings["rescan_files"]):
    print "reading data from", datafile
    parameters = load(datafile)
  else:
    parameters = {}
    if settings["singles"]:
      for p in settings["parameters"]:
        parameters[p] = parameter(p,settings["sources"],settings)
    else:
      for p in settings["parameters"]:
        batchname = batchname_from_parameter(p)
        parameters[p] = parameter(p,settings["sources"],settings,batchname=batchname)
    print "storing data in", datafile
    save(datafile,parameters)

  return parameters

def normalize_kde(kde,x):

  area = 0.0
  dx = (max(x) - min(x))/len(x)
  for i in kde:
    area += i*dx

  norm = 1.0/area
  for n in range(len(kde)):
    kde[n] *= norm

  return kde

def plot_combined_kde(parameters,param,settings):

  print "plotting",param
  sourcenames = settings["sources"]

  ## Start plotting
  fig = figure()
  ax = fig.add_subplot(111)
  ax.set_xlabel(dchilatex[param])
  #ax.set_ylabel("$ $")

  param_obj = parameters[param]

  xaxis = param_obj.get_kde_xaxis()
  ax.set_xlim(param_obj.kde_min,param_obj.kde_max)

  if settings["fill_all"]:
    for source in sourcenames:
      ax.fill_between(xaxis, 0, param_obj.get_source_kde(source), alpha=0.7, color=colors[source])
  ax.fill_between(xaxis, 0, param_obj.get_combined_kde(), alpha=1.0, color=lightgray)

  for source in sourcenames:
    ax.plot(xaxis,param_obj.get_source_kde(source),label=r"${\rm %s}$"%source,c=colors[source])
  ax.plot(xaxis,param_obj.get_combined_kde(),label=r"${\rm combined}$",c=colors["combined"],lw=1)

#  ax.legend(loc='upper left', fancybox=True, framealpha=0.8)

  handles, labels = ax.get_legend_handles_labels()
  if settings["fill_all"]:
    for i in range(len(sourcenames)):
      nh = mpatches.Patch(facecolor=colors[sourcenames[i]], alpha=0.7, edgecolor=colors[sourcenames[i]], label='')
      handles[i] = nh
  nh = mpatches.Patch(facecolor=lightgray, alpha=1.0, edgecolor=colors["combined"], label='')
  handles[-1] = nh
  ax.legend(handles, labels, loc='upper left', fancybox=True, framealpha=0.8)

  if settings["singles"]:
    figname = path.join(settings["outpath"],"posteriors_%s"%param)
  else:
    figname = path.join(settings["outpath"],"posteriors_%s_batch"%param)

  fig.savefig(figname+".pdf", bbox_inches = 'tight')
  fig.savefig(figname+".png", bbox_inches = 'tight')




###########################################
#  Do stuff
###########################################

parameters = init(settings)
for param in settings["parameters"]:
  plot_combined_kde(parameters,param,settings)

for param in settings["parameters"]:
  p = parameters[param]
  p.write_combined_posterior_samples(settings["combinedposterior_samples"])


###############################################################################
## Create violin plotting script compatible files
###############################################################################

#The violin plotting script expects a batch file to contain collumns for all the
#parameters, not one file per parameter.
#This ugly routine creates such batch files.
if not settings["singles"]:
  import fileinput
  batches = ["dchi","dalpha","dsigma","dbeta"]
  for b in batches:
    batchfilename_out = path.join(settings["datadir"],"posterior_samples_%ss.dat"%b)
    thisbatch = []
    paramnamesinbatch = []
    for param in settings["parameters"]:
      if b in param:
        paramnamesinbatch.append(param)
        samples = genfromtxt(path.join(settings["datadir"],"combined_posterior_batch_%s.dat"%param))
        thisbatch.append(samples)
    header = " ".join(paramnamesinbatch)
    savetxt(batchfilename_out,array(thisbatch).T)
    #add the header in a tremendously ugly, but effective fashion:
    for line in fileinput.input([batchfilename_out], inplace=True):
        if fileinput.isfirstline():
            print header
        print line,



