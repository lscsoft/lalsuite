#/bin/python

'''

Modified Dispersion Relation post-processing module for PE 

'''

import sys
import os
import argparse

import matplotlib
matplotlib.use('Agg')
from pylab import *
from scipy import stats
from scipy import integrate
from scipy.optimize import newton 
from numpy import random 
import h5py
import lal
from pylal.bayespputils import calculate_redshift, DistanceMeasure, lambda_a, amplitudeMeasure

###########################################
#
#  Define hardcoded
#
###########################################

fig_width_pt = 2*246.0/1.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]

plot_params = {'axes.grid' : True,
               'savefig.dpi' : 150,
               'axes.labelsize': 16, 
               'axes.titlesize': 20, 
               'font.size': 16, 
               #'font.family': 'serif', 
               'legend.fontsize': 10, 
               'xtick.labelsize': 12, 
               'ytick.labelsize': 0, 
               'xtick.major.size': 8,
               'ytick.major.size': 8,
               'xtick.major.pad': 6,
               'xtick.minor.size': 5,
               'ytick.minor.size': 5,
               'axes.grid' : True,
               'text.usetex': False,
               'lines.markersize' : 4, 
#               'lines.linewidth' : 2, 
               'lines.markeredgewidth' : 1, 
               'figure.figsize': fig_size
               }
          
rcParams.update(plot_params)

parlabels = {'log10lambda_a':r'$\log \lambda_\mathbb{A}$',
                    'log10lambda_eff':r'$\log \lambda_{\mathrm{eff}}$',
                    'lambda_a':r'$\lambda_\mathbb{A}$',
                    'lambda_eff':r'$\lambda_{\mathrm{eff}}$'}

lambda_a_min = 1.0e13
lambda_a_max = 1.0e26
KDE_bandwidth = 0.02
x_min = log10(lambda_a_min)
x_max = log10(lambda_a_max)

def mirrorEdges(data, left, right, margin=0.2, weights=None):
  '''Mirror data left and right of the edges to reduce kde edge effects'''
  lmirr = left - (data - left)
  rmirr = right - (data - right)
  mwidth = (right - left)*margin
  ledge = left - mwidth
  redge = right + mwidth
  lidx = where(lmirr > ledge)[0]
  ldata = lmirr[lidx]
  ridx = where(rmirr < redge)[0]
  rdata = rmirr[ridx]
  extdata = hstack((ldata, data, rdata))
  if weights is not None:
    lw = weights[lidx]
    rw = weights[ridx]
    extweights = hstack((lw, weights, rw))
    return extdata, extweights
  return extdata

def EnergyScale(lambda_A):
    """
    Calculate mass/energy scale in eV as:    
    m c^2 = h c / \lambda_A
    Valid for alpha != 2
    """
    return 4.135667662E-15*lal.C_SI/lambda_A

def A_LIV(lambda_A, alpha):
    """
    Calculate \mathbb{A}*c^{2-alpha} in (eV)^{2-a} as: 
    A = (m_A*c)^{2-alpha} = (h/lambda_A)^{2-alpha}
    Valid for alpha != 2
    """
    return (4.135667662E-15*lal.C_SI/lambda_A)**(2-alpha)

def lambda_A_of_eff(leffedata, zdata, alphaLIV, cosmology):
    """
    Interface function for lambda_a
    """
    dist = vectorize(lal.LuminosityDistance, excluded=[0])(cosmology, zdata)
    return lambda_a(zdata, alphaLIV, leffedata, dist)

def weighted_1dQuantile(quantile, data, weights=None):
  '''Calculate quantiles for confidence intervals on sorted weighted posteriors'''

  if weights is None:
    norm = len(data)
    weights = ones(len(data))
  else:
    norm = sum(weights)

  ### Sort data and weights array
  wdata = array(zip(data, weights) ,dtype=[('data','<f8'),('weights','<f8')])

  swdata = sort(wdata, order='data')
  sdata = swdata['data']
  sweights = swdata['weights']

  cumw = cumsum(sweights)/norm
  iplus = where(cumw > quantile)[0][0]
  x = sdata[iplus]
  # Interpolate if possible
  if iplus > 0:
    x -= (cumw[iplus] - quantile) * (sdata[iplus] - sdata[iplus-1])/(cumw[iplus] - cumw[iplus-1])
  else:
    print "WARNING: Quantile interpolation not possible, edge value used."
  return x

def normalize(pdf, x):
    '''Normalize a 1-dim posterior function on a given grid'''
    return pdf/trapz(pdf,x)

def downsample(data, k):
  '''downsample a dataset'''
  n = len(data)
  s = random.choice(n, k, replace=False)
  ds_data = data[s]
  return ds_data
  
  
def merge_posterior_samples(post_list, weights_list, label_list, weighted=False):
  '''Merge posterior samples from different runs with e.g. different template approximants'''
  merged_label = "_".join(label_list)
  n_list = [len(post) for post in post_list]
  ntot = sum(array(n_list))
  N = len(post_list)
  print n_list, len(post_list), ntot

  # Check if columns for all parameters exist in all datasets
  # for data in post_list:
  #  if not set(param_list) < set(data.dtype.names):
  #      print "ERROR: No data for " + param_list
  #      sys.exit(-1)
    # Now only keep the wanted columns
  #  data = data[param_list]

  # Check weighted flag; if true keep all samples and apply appropriate weights
  if weighted:
    merged_post = concatenate(post_list)
    # Populate a list of arrays with the weights of each posterior sample.
    # Each run i must carry the same overall weight, \propto 1/n_i
    # which we normalize according to the total number of samples ntot.
    new_weights_list = [ntot/(len(post)*N)*ones(len(post))*weights for post,weights in zip(post_list,weights_list)]
    return merged_post, concatenate(new_weights_list), merged_label
  # If not weighted, then downsample each run to the smallest n_i in the list
  nmin = min(n_list)
  ds_post=[]
  ds_weight=[]
  for post, weight in zip(post_list, weights_list):
    s = random.choice(len(post), nmin, replace=False)
    ds_post.append(post[s])
    ds_weight.append(weight[s])
  merged_post = concatenate(ds_post)
  merged_weight = concatenate(ds_weight)
  return merged_post, merged_weight, merged_label 
  

def output_stats(outdir, post, weights, label):
    if alphaLIV < 2.0:
        """Calculating Posterior Quantiles (lower)"""
        PQ_68 = weighted_1dQuantile(0.32, post, weights)
        PQ_90 = weighted_1dQuantile(0.1, post, weights) 
        PQ_95 = weighted_1dQuantile(0.05, post, weights)
        PQ_99 = weighted_1dQuantile(0.01, post, weights)
    elif alphaLIV > 2.0:
        """Calculating Posterior Quantiles (upper)"""
        PQ_68 = -weighted_1dQuantile(0.32, -post,weights)
        PQ_90 = -weighted_1dQuantile(0.1, -post, weights)
        PQ_95 = -weighted_1dQuantile(0.05, -post, weights)
        PQ_99 = -weighted_1dQuantile(0.01, -post, weights)
    else:
        print "Cannot handle alpha=2 yet. Exiting..."
        sys.exit(-1)

    # print min(post)
    print " Summary ( alpha = ", str(alphaLIV), ")"
    print "=-=-=-=-=-=-=-=-=-=-=-=-="
    print "shape:", shape(post), " min:", min(post), " max:", max(post)
    print "log(lambda_A)\t68% PQ: ", PQ_68, "\t90% PQ: ", PQ_90, "\t95% PQ: ", PQ_95, "\t99% PQ: ", PQ_99
    print "lambda_A [m]\t68% PQ: ", 10**PQ_68, "\t90% PQ: ", 10**PQ_90, "\t95% PQ: ", 10**PQ_95, "\t99% PQ: ", 10**PQ_99
    print "E_A [eV]\t68% PQ: ", EnergyScale(10**PQ_68), "\t90% PQ: ", EnergyScale(10**PQ_90), "\t95% PQ: ", EnergyScale(10**PQ_95), "\t99% PQ: ", EnergyScale(10**PQ_99)
    print "A [(eV/c)^" + str(2-alphaLIV) + "]\t68% PQ: ", A_LIV(10**PQ_68, alphaLIV), "\t90% PQ: ", A_LIV(10**PQ_90, alphaLIV), "\t95% PQ: ", A_LIV(10**PQ_95, alphaLIV), "\t99% PQ: ", A_LIV(10**PQ_99, alphaLIV)
    
    savetxt(os.path.join(outdir, 'credible_regions_%s_%s.txt'%(combine_param, label)),
             array([PQ_68, PQ_90, PQ_95, PQ_99]).transpose(),
             header = "68%\t90%\t95%\t99%",
             fmt = '%10.7f')


if __name__ == "__main__":

###########################################
#
#  Parse arguments (argparser)
#
###########################################

  parser = argparse.ArgumentParser(description="This is a post-processing script for calculating bounds on parameterized modified dispersion relations.")
  parser.add_argument("-i", "--input", type=str, dest="datafiles", nargs='+', help="list of .hdf5 or .dat file(s) containing data points (one file per chain)", metavar="POSTERIOR FILES")
  parser.add_argument("-l", "--label", type=str, dest="labels", nargs='+', help="source-identifying string(s)", default=None)
  parser.add_argument("-o", "--output", type=str, dest="outputfolder", help="outputfolder", metavar="OUTPUT FOLDER",default=".")
  parser.add_argument("--prior", type=str, dest="prior", choices=['mass','A'], help="use prior uniform in {mass, A} (default is uniform in lalinference parameter used)")
  parser.add_argument("-a", "--alpha", type=float, dest="alphaLIV", help="Exponent of Lorentz invariance violating term", default=0.0)
  parser.add_argument("-s", "--single-source", action="store_true", dest="single_source", help="All posterior input files will be marginalized over as coming from a single source", default=False)
  parser.add_argument("-n", type=int, dest="nbins", help="number of interpolation points for the  gaussian KDE (default 256)", default=256)
  parser.add_argument("-r", type=int, dest="rseed", help="Random seed for re-sampling.", default=12345)


  args = parser.parse_args()
  
  datafiles = args.datafiles
  labels = args.labels
  outfolder = args.outputfolder
  prior = args.prior
  alphaLIV = args.alphaLIV
  # combine_param = args.combine_param
  combine_param = "log10lambda_a"
  merge_weights = True

  random.seed(args.rseed)
  
  cosmology = lal.CreateCosmologicalParameters(0.7,0.3,0.7,-1.0,0.0,0.0) ## these are dummy parameters that are being used for the initialization. they are going to be set to their defaults Planck 2015 values in the next line
  lal.SetCosmologicalParametersDefaultValue(cosmology) ## setting h, omega_matter and omega_lambda to their default Planck 2015 values available in LAL

  # Prior PDFs for log10lambda_A
  priorPDF = {}
  priorPDF['mass'] = lambda loglA,loglAmin,loglAmax: log(10)*(pow(10.0,loglAmax)-pow(10.0,loglAmin))/(pow(10.0,loglA+loglAmax+loglAmin))
  priorPDF['A'] = lambda loglA,loglAmin,loglAmax: log(10)*(alphaLIV - 2.0)*pow(10.0,loglA*(alphaLIV-2.0))/(pow(10.0,loglAmax*(alphaLIV-2.0))-pow(10.0,loglAmin*(alphaLIV-2.0)))

  if not os.path.exists(outfolder):
    os.makedirs(outfolder)
  
  if labels:
    if len(set(labels))!=len(datafiles):
      print "ERROR: need to give same number of datafiles and discrete labels"
      sys.exit(-1)

  print datafiles
  wpostlist = []
  postlist = []
  weightlist = []
  logweightlist = []
  for (dfile, lab) in zip(datafiles, labels):
    print "Post-processing " + lab
    if os.path.splitext(dfile)[1] == '.hdf5':
      #  with h5py.File(dfile, 'r') as h5file:
      h5file = h5py.File(dfile, 'r') 
      ns = h5file['lalinference']['lalinference_nest']['posterior_samples']
      data = array(ns[()])
    else:
      if os.path.splitext(dfile)[1] !=  '.dat':
        print 'WARNING: data format seems to be incompatible...'
      data = genfromtxt(dfile, names=True)


    """Converting (log)distance posterior to Mpc"""
    if "logdistance" in data.dtype.names:
      distdata = exp(data["logdistance"]) # calculate_redshift needs distances in Mpc. Use * 1e6 * lal.PC_SI to convert to meters
      print "Logarithmic distance parameter detected."
    elif "distance" in data.dtype.names:
      distdata = data["distance"] # calculate_redshift needs distances in Mpc. * 1e6 * lal.PC_SI to convert to meters
      print "Linear distance parameter detected."
    else:
      print "ERROR: No distance posterior! Exiting..."
      sys.exit(-1)

    """Calculating redshifts"""
    zdata = vectorize(calculate_redshift)(distdata, cosmology.h, cosmology.om, cosmology.ol, cosmology.w0)

    #logldata = data["logl"]

    """Calculating posteriors for lambda_{eff} parameters"""
    if "log10lambda_a" in data.dtype.names:
      loglamAdata = data["log10lambda_a"]
      lamAdata = pow(10, loglamAdata)
    elif "lambda_a" in data.dtype.names:
      lamAdata = data["lambda_a"]
      loglamAdata = log10(lamAdata)
    elif "lambda_eff" in data.dtype.names:
      leffdata = data["lambda_eff"]
      logleffdata = log10(leffdata)
      lamAdata = lambda_A_of_eff(leffdata, zdata, alphaLIV, cosmology)
      lameff = True
    elif "log10lambda_eff" in data.dtype.names:
      logleffdata = data["log10lambda_eff"]
      leffdata = pow(10, logleffdata)
      lamAdata = lambda_A_of_eff(leffdata, zdata, alphaLIV, cosmology)
      loglamAdata = log10(lamAdata)
      lameff = True
    else:
      print "No valid lambda parameter found in " + lab + " data! Exiting..."
      sys.exit(-1)

    mgdata = EnergyScale(lamAdata)
    if prior == 'mass':
        # apply uniform mass prior
        print 'Applying prior uniform in mass scale'
        print "Weighing loglambda_A posterior points by 1/\lambda_\mathbb{A}"
        logweights = 1.0/lamAdata
        print "Weighing lambda_A posterior points by 1/\lambda_\mathbb{A}^2"
        weights = lamAdata**(-2)
    elif prior == 'A':
        #apply uniform A prior
        print 'Applying prior uniform in A'
        print "Weighing loglambda_A posterior points by \lambda_\mathbb{A}^{alpha-2}"
        logweights = lamAdata**(alphaLIV-2)
        print "Weighing lambda_A posterior points by \lambda_\mathbb{A}^{alpha-3}"
        weights = lamAdata**(alphaLIV-3)
    else:
        print "Using default priors"
        weights = lamAdata*0+1.0
        logweights = lamAdata*0+1.0
        
    output_stats(outfolder, loglamAdata, logweights, lab)

    #wpostlist.append(mirrorEdges(loglamAdata, x_min, x_max, weights=logweights)) 
    #postlist.append(mirrorEdges(loglamAdata, x_min, x_max))
    postlist.append(loglamAdata)
    logweightlist.append(logweights)
    weightlist.append(weights)

  if args.single_source:
    print 'Merging posteriors'
    merged_loglamA, merged_logweights, merged_label = merge_posterior_samples(postlist, logweightlist, labels, weighted=True)
    print shape(merged_loglamA)
    print shape(merged_logweights)
    output_stats(outfolder, merged_loglamA, merged_logweights, merged_label)
    outdata = vstack((merged_loglamA, merged_logweights)).T
    savetxt(os.path.join(outfolder,'merged_loglambdaA_posteriors_%s.dat'%merged_label), outdata, header = "log10lambda_a\tweights")
    
  if len(postlist) == 1 or args.single_source:
    print 'DONE! (single source only)'
    sys.exit(0)

  # Get combined posterior for loglambda_A
  print "\n==============================\nCombining all of the above sources\n==============================\n"
  x_min = min(map(min,postlist))-1 # one decade wider for the KDE to fit the plot
  x_max = max(map(max,postlist))+1 # one decade wider for the KDE to fit the plot
  print "log10lambda_a in [" + str(x_min) + ',' + str(x_max) + ']'
  x = linspace(x_min, x_max, args.nbins)

  # KDE for each source
  kdes = map(lambda x : stats.gaussian_kde(x, bw_method=0.02), postlist)

  # Prior PDF on x
  yp = priorPDF[prior](x, x_min, x_max)

  # evaluate the KDEs on the bins
  pdfs = squeeze(array([map(kde, x) for kde in kdes]), axis=2)

  # compute the combined pdf accounting for the prior
  combined_pdf = sum(log(pdfs), axis=0) + log(yp)
  combined_pdf = normalize(exp(combined_pdf), x)

  # array of unique colors for the single pdfs
  colors = cm.rainbow(linspace(0,1,len(postlist)))

  # calculate combined CDF
  from scipy.integrate import cumtrapz
  cdf = cumtrapz(combined_pdf, x, initial=0.0)

  combined_90bound = x[abs(cdf-0.1).argmin()] if alphaLIV < 2.0 else x[abs(cdf - 0.9).argmin()]
  print "COMBINED 90pc BOUND : ", pow(10.0,combined_90bound)

  # Plot posteriors and combined PDF/CDF (all normalized)
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.set_title("$\\alpha = " + str(alphaLIV) + "$")
  ax.set_xlim([x_min, x_max])

  # Plot bounds as vertical lines
  if alphaLIV==0:
    lambdag_bounds = {"Solar system":(2.8e15,'r','-.'), "Binary Pulsars":(1.6e13,'g','--')}
    for lb in lambdag_bounds:
      (bound, lc, ls) = lambdag_bounds[lb]
      ax.axvline(log10(bound), linestyle=ls, linewidth=2, color=lc, label=lb)
  ax.axvline(combined_90bound, linewidth=2, color='k', label="GW combined")

  for c,pdf,lab in zip(colors, pdfs, labels):
    ax.plot(x, normalize(pdf*yp, x), color=c, label=lab)
  ax.plot(x, combined_pdf,color='k')
  ax.plot(x, cdf if alphaLIV < 2.0 else 1.0 - cdf,color='k')

  # Labels and legend
  try:
    ax.set_xlabel(parlabels[combine_param])
  except:
    ax.set_xlabel(combine_param)
  ax.set_ylabel(r"$\mathrm{probability}$ $\mathrm{density}$")
  ax.legend(loc='upper right')

  savefig(os.path.join(outfolder,'comb_post_%s.pdf'%combine_param),bbox_inches='tight')

  # Compute a few credible regions from the combined CDF
  regions = matrix([pow(10,x[abs(cdf-cl).argmin()]) if alphaLIV < 2.0 else pow(10,x[abs(cdf-(1.-cl)).argmin()]) for cl in [0.05,0.1,0.32,0.50,0.68,0.9,0.95]])

  savetxt(os.path.join(outfolder,'comb_credible_regions_%s.txt'%combine_param), regions, header = "5%\t10%\t32%\t50%\t68%\t90%\t95%")

  print "DONE!"
