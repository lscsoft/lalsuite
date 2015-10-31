"""
Perform the consistency check between the inspiral and ringdown estimates of the mass and spin of the final 
black hole in a binary black hole merger. 

P. Ajith, Abhirup Ghosh, Archisman Ghosh, 2015-09-18

$Id:$
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.signal as ss
from scipy import interpolate
from optparse import OptionParser
import time, os 
import imrtestgr as tgr 
import pickle, gzip
import sys 

from matplotlib import rc
import matplotlib
matplotlib.rc('text.latex', preamble = '\usepackage{txfonts}')

rc('text', usetex=True)
rc('font', family='serif')
rc('font', serif='times')
rc('mathtext', default='sf')
rc("lines", markeredgewidth=1)
rc("lines", linewidth=2)
rc('axes', labelsize=10) #24
rc("axes", linewidth=0.5) #2)
rc('xtick', labelsize=8)
rc('ytick', labelsize=8)
rc('legend', fontsize=10) #16
rc('xtick.major', pad=6) #8)
rc('ytick.major', pad=6) #8)
rc('xtick.minor', size=5) #8)
rc('ytick.minor', size=5) #8)

def set_tick_sizes(ax, major, minor):
  for l in ax.get_xticklines() + ax.get_yticklines():
    l.set_markersize(major)
  for tick in ax.xaxis.get_minor_ticks() + ax.yaxis.get_minor_ticks():
    tick.tick1line.set_markersize(minor)
    tick.tick2line.set_markersize(minor)
  ax.xaxis.LABELPAD=10.
  ax.xaxis.OFFSETTEXTPAD=10.

# Module for confidence calculations
class confidence(object):
  def __init__(self, counts):
    # Sort in descending order in frequency
    self.counts_sorted = np.sort(counts.flatten())[::-1]
    # Get a normalized cumulative distribution from the mode
    self.norm_cumsum_counts_sorted = np.cumsum(self.counts_sorted) / np.sum(counts)
    # Set interpolations between heights, bins and levels
    self._set_interp()
  def _set_interp(self):
    self._length = len(self.counts_sorted)
    # height from index
    self._height_from_idx = interpolate.interp1d(np.arange(self._length), self.counts_sorted, bounds_error=False, fill_value=0.)
    # index from height
    self._idx_from_height = interpolate.interp1d(self.counts_sorted[::-1], np.arange(self._length)[::-1], bounds_error=False, fill_value=self._length)
    # level from index
    self._level_from_idx = interpolate.interp1d(np.arange(self._length), self.norm_cumsum_counts_sorted, bounds_error=False, fill_value=1.)
    # index from level
    self._idx_from_level = interpolate.interp1d(self.norm_cumsum_counts_sorted, np.arange(self._length), bounds_error=False, fill_value=self._length)
  def level_from_height(self, height):
    return self._level_from_idx(self._idx_from_height(height))
  def height_from_level(self, level):
    return self._height_from_idx(self._idx_from_level(level))

######################################################################################
################################# MAIN PROGRAM #######################################
######################################################################################

if __name__ == '__main__':

  start_time = time.time()

  # read inputs from command line 
  parser = OptionParser()
  parser.add_option("-i", "--insp-post", dest="insp_post", help="file containing the posterior samples from the lalinference inspiral run")
  parser.add_option("-r", "--ring-post", dest="ring_post", help="file containing the posterior samples from the lalinference ringdown run")
  parser.add_option("-m", "--imr-post", dest="imr_post", help="file containing the posterior samples from the full lalinference IMR run")
  parser.add_option("-f", "--fit-formula", dest="fit_formula", help="fitting formula to be used for the calculation of final mass/spin [options: 'nospin_Pan2011', 'nonprecspin_Healy2014'")
  parser.add_option("-p", "--mf-af-prior", dest="prior_Mfaf_file", help="pickle file containing the interpolation object of the prior in (Mf,af) used in the lalinference runs")
  parser.add_option("-o", "--out-dir", dest="out_dir", help="output directory")
  parser.add_option("--insp-fhigh", dest="insp_fhigh", help="Upper cutoff freq for the inspiral analysis")
  parser.add_option("--ring-flow", dest="ring_flow", help="Lower cutoff freq for the ringdown analysis")
  parser.add_option("-M", "--mtot-inj", dest="M_inj", help="injected value of total mass (if this is an injection)")
  parser.add_option("-Q", "--mratio-inj", dest="q_inj", help="injected value of mass ratio (if this is an injection)")
  parser.add_option("--chi1-inj", dest="chi1_inj", help="injected value of spin of mass m1 (if this is an injection)")
  parser.add_option("--chi2-inj", dest="chi2_inj", help="injected value of spin of mass m2 (if this is an injection)")
  parser.add_option("-w", "--waveform", dest="waveform", help="waveform used for recovery")
  parser.add_option("-d", "--debug-plots", dest="debug_plots", help="debug plots")
  (options, args) = parser.parse_args()

  insp_post = options.insp_post
  ring_post = options.ring_post
  imr_post = options.imr_post
  prior_Mfaf_file = options.prior_Mfaf_file
  out_dir = options.out_dir
  fit_formula = options.fit_formula
  debug_plots = options.debug_plots
  insp_fhigh = float(options.insp_fhigh)
  ring_flow = float(options.ring_flow)
  waveform = options.waveform
  N_bins = 201

  # create output directory and copy the script 
  os.system('mkdir -p %s' %out_dir)
  os.system('mkdir -p %s/data' %out_dir)
  os.system('mkdir -p %s/img' %out_dir)
  os.system('cp %s %s' %(__file__, out_dir))
  os.system('cp ../src/imrtestgr.py %s' %(out_dir))
  os.system('cp ../templates/*.html %s' %(out_dir)) # this assumes that we are running the script from the scripts directory FIXME 
  os.system('cp ../templates/*.css %s' %(out_dir))
  os.system('cp ../templates/*.js %s' %(out_dir))
  os.system('cp %s %s/data'%(prior_Mfaf_file, out_dir))

  # creating file to save the run command
  run_command = open('%s/command.txt'%(out_dir),'w')
  for arg in sys.argv:
    run_command.write('%s ' %arg)
  run_command.close()

  # creating the parameter table
  param_table = [['Upper cutoff freq for the inspiral analysis: %s Hz'%insp_fhigh],
                ['Lower cutoff freq for the ringdown analysis: %s Hz'%ring_flow],
                ['Waveform approximant: %s'%(waveform)],
                ['Final mass/spin fitting formula: %s'%(fit_formula)]]
  np.savetxt('%s/summary_table.txt'%(out_dir), np.array(param_table), delimiter='\t', fmt='%s')

  # creating soft links for lalinference results
  insp_posplots = insp_post.replace("/posterior_samples.dat"," ")
  ring_posplots = ring_post.replace("/posterior_samples.dat"," ")
  imr_posplots = imr_post.replace("/posterior_samples.dat"," ")

  os.system('ln -s %s %s/lalinf_insp' %(insp_posplots, out_dir))
  os.system('ln -s %s %s/lalinf_ring' %(ring_posplots, out_dir))
  os.system('ln -s %s %s/lalinf_imr' %(imr_posplots, out_dir))

  # read the injection mass parameters if this is an injection
  M_inj = options.M_inj
  q_inj = options.q_inj
  chi1_inj = options.chi1_inj
  chi2_inj = options.chi2_inj

  if M_inj == None or q_inj == None or chi1_inj == None or chi2_inj == None:
    plot_injection_lines = False
  else:
    plot_injection_lines = True
    M_inj = float(options.M_inj)
    q_inj = float(options.q_inj)
    chi1_inj = float(options.chi1_inj)
    chi2_inj = float(options.chi2_inj)
    m2_inj = M_inj/(1.+q_inj)
    m1_inj = M_inj*q_inj/(1.+q_inj)
    eta_inj = q_inj/(1.+q_inj)**2.
    Mf_inj, af_inj = tgr.calc_final_mass_spin(m1_inj, m2_inj, chi1_inj, chi2_inj, fit_formula)

  ###############################################################################################
  # Read the posteriors from the inspiral, ringdown and imr lalinference runs (after post-processing) 
  ###############################################################################################
  insp_data = np.genfromtxt(insp_post, dtype=None, names=True)
  m1_i, m2_i = insp_data['m1'], insp_data['m2']
  if ('a1' in insp_data.dtype.names) and ('a2' in insp_data.dtype.names):
    chi1_i, chi2_i = insp_data['a1'], insp_data['a2']
  else:
    chi1_i, chi2_i = np.zeros(len(m1_i)), np.zeros(len(m2_i))
  Mf_i, af_i = tgr.calc_final_mass_spin(m1_i, m2_i, chi1_i, chi2_i, fit_formula)

  ring_data = np.genfromtxt(ring_post, dtype=None, names=True)
  m1_r, m2_r = ring_data['m1'], ring_data['m2']
  if ('a1' in ring_data.dtype.names) and ('a2' in ring_data.dtype.names):
    chi1_r, chi2_r = ring_data['a1'], ring_data['a2']
  else:
    chi1_r, chi2_r = np.zeros(len(m1_r)), np.zeros(len(m2_r))
  Mf_r, af_r = tgr.calc_final_mass_spin(m1_r, m2_r, chi1_r, chi2_r, fit_formula)

  imr_data = np.genfromtxt(imr_post, dtype=None, names=True)
  m1_imr, m2_imr = imr_data['m1'], imr_data['m2']
  if ('a1' in imr_data.dtype.names) and ('a2' in imr_data.dtype.names):
    chi1_imr, chi2_imr = imr_data['a1'], imr_data['a2']
  else:
    chi1_imr, chi2_imr = np.zeros(len(m1_imr)), np.zeros(len(m2_imr))
  Mf_imr, af_imr = tgr.calc_final_mass_spin(m1_imr, m2_imr, chi1_imr, chi2_imr, fit_formula)

  print '... read posteriors'
  ###############################################################################################

  ###############################################################################################
  # compute the limits of integration for computing delta_Mf and delta_af 
  ###############################################################################################
  Mf_lim = max(abs(np.append(np.append(Mf_i, Mf_r), Mf_imr)))
  af_lim = max(abs(np.append(np.append(af_i, af_r), af_imr)))

  Mf_bins = np.linspace(-Mf_lim, Mf_lim, N_bins)
  af_bins = np.linspace(-af_lim, af_lim, N_bins)

  dMf = np.mean(np.diff(Mf_bins))
  daf = np.mean(np.diff(af_bins))

  Mf_intp = (Mf_bins[:-1] + Mf_bins[1:])/2.
  af_intp = (af_bins[:-1] + af_bins[1:])/2.
  ###############################################################################################


  ###############################################################################################
  # Undo the effect of the prior from the lalinference posterior. Lalinference assumes a        #
  # uniform prior in component masses. We need to assume a uniform prior in Mf, af              #
  ###############################################################################################

  # read the interpolation object, reconstruct the data from the interpolation object 
  f = gzip.open(prior_Mfaf_file,'rb')
  P_Mfaf_pr_interp_obj = pickle.load(f)
  P_Mfaf_pr = P_Mfaf_pr_interp_obj(Mf_intp, af_intp)

  # compute the 2D posterior distributions for the inspiral, ringodwn and IMR analyses 
  P_Mfaf_i, Mf_bins, af_bins = np.histogram2d(Mf_i, af_i, bins=(Mf_bins, af_bins), normed=True)
  P_Mfaf_r, Mf_bins, af_bins = np.histogram2d(Mf_r, af_r, bins=(Mf_bins, af_bins), normed=True)
  P_Mfaf_imr, Mf_bins, af_bins = np.histogram2d(Mf_imr, af_imr, bins=(Mf_bins, af_bins), normed=True)

  P_Mfaf_i = P_Mfaf_i.T
  P_Mfaf_r = P_Mfaf_r.T
  P_Mfaf_imr = P_Mfaf_imr.T

  # compute the corrected 2D posteriors in Mf and af by dividing by the prior distribution 
  P_Mfaf_i = P_Mfaf_i/P_Mfaf_pr
  P_Mfaf_r = P_Mfaf_r/P_Mfaf_pr
  P_Mfaf_imr = P_Mfaf_imr/P_Mfaf_pr

  # removing nan's
  P_Mfaf_i[np.isnan(P_Mfaf_i)] = 0.
  P_Mfaf_r[np.isnan(P_Mfaf_r)] = 0.
  P_Mfaf_imr[np.isnan(P_Mfaf_imr)] = 0.

  # removing infinities
  P_Mfaf_i[np.isinf(P_Mfaf_i)] = 0.
  P_Mfaf_r[np.isinf(P_Mfaf_r)] = 0.
  P_Mfaf_imr[np.isinf(P_Mfaf_imr)] = 0.

  print '... computed (prior) corrected posteriors'
  ###############################################################################################

  ################################################################################################
  # compute the posterior of (delta_Mf, delta_af)
  ################################################################################################
  P_dMfdaf = dMf*daf*ss.correlate2d(P_Mfaf_i, P_Mfaf_r, boundary='fill', mode='same')

  print '... computed P(delta_Mf, delta_af)'
  ###############################################################################################

  ################################################################################################
  # compute the posterior of (delta_Mf/Mf, delta_af/af)
  ################################################################################################
  # compute interpolation objects for the Mf,af posterior and delta_Mf and delta_af posterior 
  Mf_intp = (Mf_bins[:-1] + Mf_bins[1:])/2.
  af_intp = (af_bins[:-1] + af_bins[1:])/2.
  P_dMfdaf_interp_object = scipy.interpolate.interp2d(Mf_intp, af_intp, P_dMfdaf, fill_value=0., bounds_error=False)
  P_Mfaf_imr_interp_object = scipy.interpolate.interp2d(Mf_intp, af_intp, P_Mfaf_imr, fill_value=0., bounds_error=False)

  # defining limits of delta_Mf/Mf and delta_af/af. limits are currently set arbitrarily FIXME 
  dMfbyMf_vec = np.linspace(-1.0, 1.0, N_bins)
  dafbyaf_vec = np.linspace(-1.0, 1.0, N_bins)

  # compute the P(dMf/Mf, daf/af) by evaluating the integral 
  dx = np.mean(np.diff(dMfbyMf_vec))
  dy = np.mean(np.diff(dafbyaf_vec))
  P_dMfbyMf_dafbyaf = np.zeros(shape=(N_bins,N_bins))

  for i, v2 in enumerate(dafbyaf_vec):
    for j, v1 in enumerate(dMfbyMf_vec):
      P_dMfbyMf_dafbyaf[i,j] = tgr.calc_sum(Mf_intp, af_intp, v1, v2, P_dMfdaf_interp_object, P_Mfaf_imr_interp_object)*dx*dy

  # normalization
  P_dMfbyMf_dafbyaf /= np.sum(P_dMfbyMf_dafbyaf) * dx * dy

  # Marginalization to one-dimensional joint_posteriors
  P_dMfbyMf = np.sum(P_dMfbyMf_dafbyaf, axis=0) * dy
  P_dafbyaf = np.sum(P_dMfbyMf_dafbyaf, axis=1) * dx
  
  # injection confidence
  conf_v1v2 = confidence(P_dMfbyMf_dafbyaf)
  inj_height = P_dMfbyMf_dafbyaf[np.argmin(abs(dMfbyMf_vec)), np.argmin(abs(dafbyaf_vec))]
  inj_level = conf_v1v2.level_from_height(inj_height)
  print 'GR is consistent with %.1f%% confidence level'%(100.*inj_level)

  # save results 
  np.savetxt(out_dir+'/data/Mfaf.dat', (Mf_bins,af_bins))
  np.savetxt(out_dir+'/data/P_Mfaf_i.dat', P_Mfaf_i)
  np.savetxt(out_dir+'/data/P_Mfaf_r.dat', P_Mfaf_r)
  np.savetxt(out_dir+'/data/P_Mfaf_imr.dat', P_Mfaf_imr)
  np.savetxt(out_dir+'/data/P_dMfdaf.dat', P_dMfdaf)
  np.savetxt(out_dir+'/data/dMfbyMf_vec.dat', dMfbyMf_vec)
  np.savetxt(out_dir+'/data/dafbyaf_vec.dat', dafbyaf_vec)
  np.savetxt(out_dir+'/data/P_dMfbyMf_dafbyaf.dat', P_dMfbyMf_dafbyaf)
  np.savetxt(out_dir+'/data/P_dMfbyMf.dat', P_dMfbyMf)
  np.savetxt(out_dir+'/data/P_dafbyaf.dat', P_dafbyaf)
  np.savetxt(out_dir+'/data/GR_confidence.txt', [inj_level])

  #########################################################################################

  #########################################################################################
  # plotting                  
  #########################################################################################
  #inspiral
  P_m1m2_i, m1_bins_i, m2_bins_i = np.histogram2d(m1_i, m2_i, bins=50, normed=True)
  P_chi1chi2_i, chi1_bins_i, chi2_bins_i = np.histogram2d(chi1_i, chi2_i, bins=50, normed=True)

  P_m1m2_i = P_m1m2_i.T
  P_chi1chi2_i = P_chi1chi2_i.T
  
  conf_m1m2_i = confidence(P_m1m2_i)
  s1_m1m2_i = conf_m1m2_i.height_from_level(0.68) 
  s2_m1m2_i = conf_m1m2_i.height_from_level(0.95)
  
  conf_chi1chi2_i = confidence(P_chi1chi2_i)
  s1_chi1chi2_i = conf_chi1chi2_i.height_from_level(0.68)
  s2_chi1chi2_i = conf_chi1chi2_i.height_from_level(0.95)
  
  conf_Mfaf_i = confidence(P_Mfaf_i)
  s1_Mfaf_i = conf_Mfaf_i.height_from_level(0.68)
  s2_Mfaf_i = conf_Mfaf_i.height_from_level(0.95)

  plt.figure(figsize=(5,5))
  plt.pcolormesh(m1_bins_i, m2_bins_i, tgr.gf(P_m1m2_i), cmap='YlOrBr')
  plt.contour(m1_bins_i[:-1], m2_bins_i[:-1], tgr.gf(P_m1m2_i), levels=(s1_m1m2_i,s2_m1m2_i), linewidths=(1,1.5))
  if plot_injection_lines == True:
    plt.axvline(x=m1_inj, ls='--', color='k')
    plt.axhline(y=m2_inj, ls='--', color='k')
  plt.xlabel('$m_1 [M_{\odot}]$')
  plt.ylabel('$m_2 [M_{\odot}]$')
  plt.xlim([min(m1_i), max(m1_i)])
  plt.ylim([min(m2_i), max(m2_i)])
  plt.grid()
  plt.savefig('%s/img/inspiral_m1m2_thumb.png'%(out_dir), dpi=72)
  plt.savefig('%s/img/inspiral_m1m2.png'%(out_dir), dpi=300)

  plt.figure(figsize=(5,5))
  plt.plot(m1_i, m2_i, 'k.', ms=0.1)
  if plot_injection_lines == True:
    plt.axvline(x=m1_inj, ls='--', color='k')
    plt.axhline(y=m2_inj, ls='--', color='k')
  plt.xlabel('$m_1 [M_{\odot}]$')
  plt.ylabel('$m_2 [M_{\odot}]$')
  plt.xlim([min(m1_i), max(m1_i)])
  plt.ylim([min(m2_i), max(m2_i)])
  plt.grid()
  plt.savefig('%s/img/inspiral_m1m2_scatter_thumb.png'%(out_dir), dpi=72)
  plt.savefig('%s/img/inspiral_m1m2_scatter.png'%(out_dir), dpi=300)

  plt.figure(figsize=(5,5))
  plt.pcolormesh(chi1_bins_i, chi2_bins_i, tgr.gf(P_chi1chi2_i), cmap='YlOrBr')
  plt.contour(chi1_bins_i[:-1], chi2_bins_i[:-1], tgr.gf(P_chi1chi2_i), levels=(s1_chi1chi2_i,s2_chi1chi2_i), linewidths=(1,1.5))
  if plot_injection_lines == True:
    plt.axvline(x=chi1_inj, ls='--', color='k')
    plt.axhline(y=chi2_inj, ls='--', color='k')
  plt.xlabel('$\chi _1$')
  plt.ylabel('$\chi _2$')
  plt.xlim([min(chi1_i), max(chi1_i)])
  plt.ylim([min(chi2_i), max(chi2_i)])
  plt.grid()
  plt.savefig('%s/img/inspiral_chi1chi2_thumb.png'%(out_dir), dpi=72)
  plt.savefig('%s/img/inspiral_chi1chi2.png'%(out_dir), dpi=300)

  plt.figure(figsize=(5,5))
  plt.plot(chi1_i, chi2_i, 'k.', ms=0.1)
  if plot_injection_lines == True:
    plt.axvline(x=chi1_inj, ls='--', color='k')
    plt.axhline(y=chi2_inj, ls='--', color='k')
  plt.xlabel('$\chi _1$')
  plt.ylabel('$\chi _2$')
  plt.xlim([min(chi1_i), max(chi1_i)])
  plt.ylim([min(chi2_i), max(chi2_i)])
  plt.grid()
  plt.savefig('%s/img/inspiral_chi1chi2_scatter_thumb.png'%(out_dir), dpi=72)
  plt.savefig('%s/img/inspiral_chi1chi2_scatter.png'%(out_dir), dpi=300)

  plt.figure(figsize=(5,5))
  plt.pcolormesh(Mf_bins, af_bins, tgr.gf(P_Mfaf_i), cmap='YlOrBr')
  plt.contour(Mf_bins[:-1], af_bins[:-1], tgr.gf(P_Mfaf_i), levels=(s1_Mfaf_i,s2_Mfaf_i), linewidths=(1,1.5))
  if plot_injection_lines == True:
    plt.axvline(x=Mf_inj, ls='--', color='k')
    plt.axhline(y=af_inj, ls='--', color='k')
  plt.xlabel('$M_f [M_{\odot}]$')
  plt.ylabel('$a_f/M_f$')
  plt.xlim([min(Mf_i), max(Mf_i)])
  plt.ylim([min(af_i), max(af_i)])
  plt.grid()
  plt.savefig('%s/img/inspiral_Mfaf_thumb.png'%(out_dir), dpi=72)
  plt.savefig('%s/img/inspiral_Mfaf.png'%(out_dir), dpi=300)


  #ringdown
  P_m1m2_r, m1_bins_r, m2_bins_r = np.histogram2d(m1_r, m2_r, bins=50, normed=True)
  P_chi1chi2_r, chi1_bins_r, chi2_bins_r = np.histogram2d(chi1_r, chi2_r, bins=50, normed=True)

  P_m1m2_r = P_m1m2_r.T
  P_chi1chi2_r = P_chi1chi2_r.T
  
  conf_m1m2_r = confidence(P_m1m2_r)
  s1_m1m2_r = conf_m1m2_r.height_from_level(0.68) 
  s2_m1m2_r = conf_m1m2_r.height_from_level(0.95)
  
  conf_chi1chi2_r = confidence(P_chi1chi2_r)
  s1_chi1chi2_r = conf_chi1chi2_r.height_from_level(0.68)
  s2_chi1chi2_r = conf_chi1chi2_r.height_from_level(0.95)
  
  conf_Mfaf_r = confidence(P_Mfaf_r)
  s1_Mfaf_r = conf_Mfaf_r.height_from_level(0.68)
  s2_Mfaf_r = conf_Mfaf_r.height_from_level(0.95)

  plt.figure(figsize=(5,5))
  plt.pcolormesh(m1_bins_r, m2_bins_r, tgr.gf(P_m1m2_r), cmap='YlOrBr')
  plt.contour(m1_bins_r[:-1], m2_bins_r[:-1], tgr.gf(P_m1m2_r), levels=(s1_m1m2_r,s2_m1m2_r), linewidths=(1,1.5))
  if plot_injection_lines == True:
    plt.axvline(x=m1_inj, ls='--', color='k')
    plt.axhline(y=m2_inj, ls='--', color='k')
  plt.xlabel('$m_1 [M_{\odot}]$')
  plt.ylabel('$m_2 [M_{\odot}]$')
  plt.xlim([min(m1_r), max(m1_r)])
  plt.ylim([min(m2_r), max(m2_r)])
  plt.grid()
  plt.savefig('%s/img/ringdown_m1m2.png'%(out_dir), dpi=300)
  plt.savefig('%s/img/ringdown_m1m2_thumb.png'%(out_dir), dpi=72)

  plt.figure(figsize=(5,5))
  plt.plot(m1_r, m2_r, 'k.', ms=0.1)
  if plot_injection_lines == True:
    plt.axvline(x=m1_inj, ls='--', color='k')
    plt.axhline(y=m2_inj, ls='--', color='k')
  plt.xlabel('$m_1 [M_{\odot}]$')
  plt.ylabel('$m_2 [M_{\odot}]$')
  plt.xlim([min(m1_r), max(m1_r)])
  plt.ylim([min(m2_r), max(m2_r)])
  plt.grid()
  plt.savefig('%s/img/ringdown_m1m2_scatter_thumb.png'%(out_dir), dpi=72)
  plt.savefig('%s/img/ringdown_m1m2_scatter.png'%(out_dir), dpi=300)

  plt.figure(figsize=(5,5))
  plt.pcolormesh(chi1_bins_r, chi2_bins_r, tgr.gf(P_chi1chi2_r), cmap='YlOrBr')
  plt.contour(chi1_bins_r[:-1], chi2_bins_r[:-1], tgr.gf(P_chi1chi2_r), levels=(s1_chi1chi2_r,s2_chi1chi2_r), linewidths=(1,1.5))
  if plot_injection_lines == True:
    plt.axvline(x=chi1_inj, ls='--', color='k')
    plt.axhline(y=chi2_inj, ls='--', color='k')
  plt.xlabel('$\chi _1$')
  plt.ylabel('$\chi _2$')
  plt.xlim([min(chi1_r), max(chi1_r)])
  plt.ylim([min(chi2_r), max(chi2_r)])
  plt.grid()
  plt.savefig('%s/img/ringdown_chi1chi2_thumb.png'%(out_dir), dpi=72)
  plt.savefig('%s/img/ringdown_chi1chi2.png'%(out_dir), dpi=300)

  plt.figure(figsize=(5,5))
  plt.plot(chi1_r, chi2_r, 'k.', ms=0.1)
  if plot_injection_lines == True:
    plt.axvline(x=chi1_inj, ls='--', color='k')
    plt.axhline(y=chi2_inj, ls='--', color='k')
  plt.xlabel('$\chi _1$')
  plt.ylabel('$\chi _2$')
  plt.xlim([min(chi1_r), max(chi1_r)])
  plt.ylim([min(chi2_r), max(chi2_r)])
  plt.grid()
  plt.savefig('%s/img/ringdown_chi1chi2_scatter_thumb.png'%(out_dir), dpi=72)
  plt.savefig('%s/img/ringdown_chi1chi2_scatter.png'%(out_dir), dpi=300)

  plt.figure(figsize=(5,5))
  plt.pcolormesh(Mf_bins, af_bins, tgr.gf(P_Mfaf_r), cmap='YlOrBr')
  plt.contour(Mf_bins[:-1], af_bins[:-1], tgr.gf(P_Mfaf_r), levels=(s1_Mfaf_r,s2_Mfaf_r), linewidths=(1,1.5))
  if plot_injection_lines == True:
    plt.axvline(x=Mf_inj, ls='--', color='k')
    plt.axhline(y=af_inj, ls='--', color='k')
  plt.xlabel('$M_f [M_{\odot}]$')
  plt.ylabel('$a_f/M_f$')
  plt.xlim([min(Mf_r), max(Mf_r)])
  plt.ylim([min(af_r), max(af_r)])
  plt.grid()
  plt.savefig('%s/img/ringdown_Mfaf.png'%(out_dir), dpi=300)
  plt.savefig('%s/img/ringdown_Mfaf_thumb.png'%(out_dir), dpi=72)

  #IMR
  P_m1m2_imr, m1_bins_imr, m2_bins_imr = np.histogram2d(m1_imr, m2_imr, bins=50, normed=True)
  P_chi1chi2_imr, chi1_bins_imr, chi2_bins_imr = np.histogram2d(chi1_imr, chi2_imr, bins=50, normed=True)

  P_m1m2_imr = P_m1m2_imr.T
  P_chi1chi2_imr = P_chi1chi2_imr.T
  
  conf_m1m2_imr = confidence(P_m1m2_imr)
  s1_m1m2_imr = conf_m1m2_imr.height_from_level(0.68) 
  s2_m1m2_imr = conf_m1m2_imr.height_from_level(0.95)
  
  conf_chi1chi2_imr = confidence(P_chi1chi2_imr)
  s1_chi1chi2_imr = conf_chi1chi2_imr.height_from_level(0.68)
  s2_chi1chi2_imr = conf_chi1chi2_imr.height_from_level(0.95)
  
  conf_Mfaf_imr = confidence(P_Mfaf_imr)
  s1_Mfaf_imr = conf_Mfaf_imr.height_from_level(0.68)
  s2_Mfaf_imr = conf_Mfaf_imr.height_from_level(0.95)

  plt.figure(figsize=(5,5))
  plt.pcolormesh(m1_bins_imr, m2_bins_imr, tgr.gf(P_m1m2_imr), cmap='YlOrBr')
  plt.contour(m1_bins_imr[:-1], m2_bins_imr[:-1], tgr.gf(P_m1m2_imr), levels=(s1_m1m2_imr,s2_m1m2_imr), linewidths=(1,1.5))
  if plot_injection_lines == True:
    plt.axvline(x=m1_inj, ls='--', color='k')
    plt.axhline(y=m2_inj, ls='--', color='k')
  plt.xlabel('$m_1 [M_{\odot}]$')
  plt.ylabel('$m_2 [M_{\odot}]$')
  plt.xlim([min(m1_imr), max(m1_imr)])
  plt.ylim([min(m2_imr), max(m2_imr)])
  plt.grid()
  plt.savefig('%s/img/imr_m1m2.png'%(out_dir), dpi=300)
  plt.savefig('%s/img/imr_m1m2_thumb.png'%(out_dir), dpi=72)

  plt.figure(figsize=(5,5))
  plt.plot(m1_imr, m2_imr, 'k.', ms=0.1)
  if plot_injection_lines == True:
    plt.axvline(x=m1_inj, ls='--', color='k')
    plt.axhline(y=m2_inj, ls='--', color='k')
  plt.xlabel('$m_1 [M_{\odot}]$')
  plt.ylabel('$m_2 [M_{\odot}]$')
  plt.xlim([min(m1_imr), max(m1_imr)])
  plt.ylim([min(m2_imr), max(m2_imr)])
  plt.grid()
  plt.savefig('%s/img/imr_m1m2_scatter_thumb.png'%(out_dir), dpi=72)
  plt.savefig('%s/img/imr_m1m2_scatter.png'%(out_dir), dpi=300)

  plt.figure(figsize=(5,5))
  plt.pcolormesh(chi1_bins_imr, chi2_bins_imr, tgr.gf(P_chi1chi2_imr), cmap='YlOrBr')
  plt.contour(chi1_bins_imr[:-1], chi2_bins_imr[:-1], tgr.gf(P_chi1chi2_imr), levels=(s1_chi1chi2_imr,s2_chi1chi2_imr), linewidths=(1,1.5))
  if plot_injection_lines == True:
    plt.axvline(x=chi1_inj, ls='--', color='k')
    plt.axhline(y=chi2_inj, ls='--', color='k')
  plt.xlabel('$\chi _1$')
  plt.ylabel('$\chi _2$')
  plt.xlim([min(chi1_imr), max(chi1_imr)])
  plt.ylim([min(chi2_imr), max(chi2_imr)])
  plt.grid()
  plt.savefig('%s/img/imr_chi1chi2_thumb.png'%(out_dir), dpi=72)
  plt.savefig('%s/img/imr_chi1chi2.png'%(out_dir), dpi=300)

  plt.figure(figsize=(5,5))
  plt.plot(chi1_imr, chi2_imr, 'k.', ms=0.1)
  if plot_injection_lines == True:
    plt.axvline(x=chi1_inj, ls='--', color='k')
    plt.axhline(y=chi2_inj, ls='--', color='k')
  plt.xlabel('$\chi _1$')
  plt.ylabel('$\chi _2$')
  plt.xlim([min(chi1_imr), max(chi1_imr)])
  plt.ylim([min(chi2_imr), max(chi2_imr)])
  plt.grid()
  plt.savefig('%s/img/imr_chi1chi2_scatter_thumb.png'%(out_dir), dpi=72)
  plt.savefig('%s/img/imr_chi1chi2_scatter.png'%(out_dir), dpi=300)

  plt.figure(figsize=(5,5))
  plt.pcolormesh(Mf_bins, af_bins, tgr.gf(P_Mfaf_imr), cmap='YlOrBr')
  plt.contour(Mf_bins[:-1], af_bins[:-1], tgr.gf(P_Mfaf_imr), levels=(s1_Mfaf_imr,s2_Mfaf_imr), linewidths=(1,1.5))
  if plot_injection_lines == True:
    plt.axvline(x=Mf_inj, ls='--', color='k')
    plt.axhline(y=af_inj, ls='--', color='k')
  plt.xlabel('$M_f [M_{\odot}]$')
  plt.ylabel('$a_f/M_f$')
  plt.xlim([min(Mf_imr), max(Mf_imr)])
  plt.ylim([min(af_imr), max(af_imr)])
  plt.grid()
  plt.savefig('%s/img/imr_Mfaf.png'%(out_dir), dpi=300)
  plt.savefig('%s/img/imr_Mfaf_thumb.png'%(out_dir), dpi=72)

  # IR overlap
  plt.figure(figsize=(5,5))
  CSi = plt.contour(Mf_bins[:-1], af_bins[:-1], tgr.gf(P_Mfaf_i), levels=(s1_Mfaf_i,s2_Mfaf_i), linewidths=(1,1.5), colors='orange')
  CSr = plt.contour(Mf_bins[:-1], af_bins[:-1], tgr.gf(P_Mfaf_r), levels=(s1_Mfaf_r,s2_Mfaf_r), linewidths=(1,1.5), colors='red')
  CSimr = plt.contour(Mf_bins[:-1], af_bins[:-1], tgr.gf(P_Mfaf_imr), levels=(s1_Mfaf_imr,s2_Mfaf_imr), linewidths=(1,1.5), colors='k')
  if plot_injection_lines == True:
    plt.axvline(x=Mf_inj, ls='--', color='k')
    plt.axhline(y=af_inj, ls='--', color='k')
  plt.xlim([min(np.append(Mf_i, Mf_r)), max(np.append(Mf_i, Mf_r))])
  plt.ylim([min(np.append(af_i, af_r)), max(np.append(af_i, af_r))])
  plt.xlabel('$M_f~[M_\odot]$')
  plt.ylabel('$a_f/M_f$')
  plt.grid()

  strs_i = [ 'inspiral', 'inspiral' ]
  strs_r = [ 'ringdown', 'ringdown' ]
  strs_imr = [ 'IMR', 'IMR' ]
  fmt_i = {}
  fmt_r = {}
  fmt_imr = {}
  for l,s in zip(CSi.levels, strs_i):
    fmt_i[l] = s
  for l,s in zip(CSr.levels, strs_r):
    fmt_r[l] = s
  for l,s in zip(CSimr.levels, strs_imr):
    fmt_imr[l] = s

  # Label every other level using strings
  plt.clabel(CSi,CSi.levels[::2],inline=True,fmt=fmt_i,fontsize=14, use_clabeltext=True)
  plt.clabel(CSr,CSr.levels[::2],inline=True,fmt=fmt_r,fontsize=14, use_clabeltext=True)
  plt.clabel(CSimr,CSimr.levels[::2],inline=True,fmt=fmt_imr,fontsize=10)

  plt.savefig('%s/img/IMR_overlap.png'%(out_dir), dpi=300)
  plt.savefig('%s/img/IMR_overlap_thumb.png'%(out_dir), dpi=72)

  #(dMf, daf)
  conf_dMfdaf = confidence(P_dMfdaf)
  s1_dMfdaf = conf_dMfdaf.height_from_level(0.68)
  s2_dMfdaf = conf_dMfdaf.height_from_level(0.95)

  plt.figure(figsize=(5,5))
  plt.pcolormesh(Mf_bins, af_bins, tgr.gf(P_dMfdaf), cmap='YlOrBr')
  plt.contour(Mf_bins[:-1], af_bins[:-1], tgr.gf(P_dMfdaf), levels=(s1_dMfdaf,s2_dMfdaf), linewidths=(1,1.5))
  plt.plot(0, 0, 'k+', ms=12, mew=2)
  plt.xlabel('$\Delta M_f~[M_\odot]$')
  plt.ylabel('$\Delta a_f$')
  plt.grid()
  plt.savefig('%s/img/dMfdaf.png'%(out_dir), dpi=300)
  plt.savefig('%s/img/dMfdaf_thumb.png'%(out_dir), dpi=72)

  #(dMf/Mf, daf/af)
  conf_v1v2 = confidence(P_dMfbyMf_dafbyaf)
  s1_v1v2 = conf_v1v2.height_from_level(0.68)
  s2_v1v2 = conf_v1v2.height_from_level(0.95)
  
  conf_v1 = confidence(P_dMfbyMf)
  s1_v1 = conf_v1.height_from_level(0.68)
  s2_v1 = conf_v1.height_from_level(0.95)

  conf_v2 = confidence(P_dafbyaf)
  s1_v2 = conf_v2.height_from_level(0.68)
  s2_v2 = conf_v2.height_from_level(0.95)

  # Calculation of condifence edges

  left1_v1 = min(dMfbyMf_vec[np.where(P_dMfbyMf>=s1_v1)[0]])
  right1_v1 = max(dMfbyMf_vec[np.where(P_dMfbyMf>=s1_v1)[0]])

  left2_v1 = min(dMfbyMf_vec[np.where(P_dMfbyMf>=s2_v1)[0]])
  right2_v1 = max(dMfbyMf_vec[np.where(P_dMfbyMf>=s2_v1)[0]])

  left1_v2 = min(dafbyaf_vec[np.where(P_dafbyaf>s1_v2)[0]])
  right1_v2 = max(dafbyaf_vec[np.where(P_dafbyaf>s1_v2)[0]])

  left2_v2 = min(dafbyaf_vec[np.where(P_dafbyaf>s2_v2)[0]])
  right2_v2 = max(dafbyaf_vec[np.where(P_dafbyaf>s2_v2)[0]])

  plt.figure(figsize=(5,5))
  plt.subplot2grid((3,3), (0,0), colspan=2)
  plt.plot(dMfbyMf_vec, P_dMfbyMf, color='k', lw=1)
  plt.axvline(x=left1_v1, color='c', lw=0.5, ls='-.')
  plt.axvline(x=right1_v1, color='c', lw=0.5, ls='-.')
  plt.axvline(x=left2_v1, color='b', lw=0.5, ls='-.')
  plt.axvline(x=right2_v1, color='b', lw=0.5, ls='-.')
  #plt.xlabel('$\Delta M_f/M_f$')
  plt.ylabel('$P(\Delta M_f/M_f)$')
  #plt.grid()

  plt.subplot2grid((3,3), (1,0), colspan=2, rowspan=2)
  plt.pcolormesh(dMfbyMf_vec,dafbyaf_vec,P_dMfbyMf_dafbyaf, cmap='YlOrBr')
  plt.contour(dMfbyMf_vec,dafbyaf_vec,tgr.gf(P_dMfbyMf_dafbyaf), levels=(s1_v1v2,s2_v1v2), linewidths=(1,1.5))
  plt.plot(0, 0, 'k+', ms=12, mew=2)
  plt.xlabel('$\Delta M_f/M_f$')
  plt.ylabel('$\Delta a_f/a_f$')
  plt.xlim([-1.,1.])
  plt.ylim([-1.,1.])
  plt.grid()

  plt.subplot2grid((3,3), (1,2), rowspan=2)
  plt.plot(P_dafbyaf, dafbyaf_vec,'k', lw=1)
  plt.axhline(y=left1_v2, color='c', lw=0.5, ls='-.')
  plt.axhline(y=right1_v2, color='c', lw=0.5, ls='-.')
  plt.axhline(y=left2_v2, color='b', lw=0.5, ls='-.')
  plt.axhline(y=right2_v2, color='b', lw=0.5, ls='-.')
  #plt.ylabel('$\Delta a_f/a_f$')
  plt.xlabel('$P(\Delta a_f/a_f)$')
  #plt.grid()

  plt.savefig('%s/img/dMfbyMfdafbyaf.png' %(out_dir), dpi=300)
  plt.savefig('%s/img/dMfbyMfdafbyaf_thumb.png' %(out_dir), dpi=72)

  print '... made summary plots' 

  print '... completed in %f seconds' %(time.time()-start_time)
  #########################################################################################

