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
from scipy import interpolate
from optparse import OptionParser
import time, os
import lalinference.imrtgr.imrtgrutils as tgr
import pickle, gzip
import sys
from lalinference import git_version

from scipy.stats import gaussian_kde   #rahul: for KDE implementation

from matplotlib import rc
import matplotlib
matplotlib.rc('text.latex', preamble=r'\usepackage{txfonts}')

rc('text', usetex=True)
rc('font', family='serif')
rc('font', serif='times')
rc('mathtext', default='sf')
rc("lines", markeredgewidth=1)
rc("lines", linewidth=2)
rc('axes', labelsize=10)
rc("axes", linewidth=0.5)
rc('xtick', labelsize=8)
rc('ytick', labelsize=8)
rc('legend', fontsize=10)
rc('xtick.major', pad=6)
rc('ytick.major', pad=6)
rc('xtick.minor', size=5)
rc('ytick.minor', size=5)

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
  parser.add_option("-i", "--insp-post", dest="insp_post", help="file containing the posterior samples from the PE inspiral run")
  parser.add_option("-r", "--ring-post", dest="ring_post", help="file containing the posterior samples from the PE ringdown run")
  parser.add_option("-m", "--imr-post", dest="imr_post", help="file containing the posterior samples from the full PE IMR run")
  parser.add_option("-f", "--fit-formula", dest="fit_formula", help="fitting formula to be used for the calculation of final mass/spin [options: 'nospin_Pan2011', 'nonprecspin_Healy2014', 'bbh_average_fits_precessing'", default="nonprecspin_Healy2014")
  parser.add_option("-p", "--mf-chif-prior", dest="prior_Mfchif_file", help="pickle file containing the interpolation object of the prior in (Mf, chif) used in the PE runs", default=None)
  parser.add_option("-o", "--out-dir", dest="out_dir", help="output directory")
  parser.add_option("--insp-fhigh", dest="insp_fhigh", help="Upper cutoff freq for the inspiral analysis")
  parser.add_option("--ring-flow", dest="ring_flow", help="Lower cutoff freq for the ringdown analysis")
  parser.add_option("--m1-inj", dest="m1_inj", help="injected value of component mass m1 (if this is an injection)")
  parser.add_option("--m2-inj", dest="m2_inj", help="injected value of component mass m2 (if this is an injection)")
  parser.add_option("--chi1-inj", dest="chi1_inj", help="injected value of spin magnitude of mass m1 (if this is an injection)")
  parser.add_option("--chi2-inj", dest="chi2_inj", help="injected value of spin magnitude of mass m2 (if this is an injection)")
  parser.add_option("--chi1z-inj", dest="chi1z_inj", help="injected value of z-component of spin of mass m1 (if this is an injection)")
  parser.add_option("--chi2z-inj", dest="chi2z_inj", help="injected value of z-component of spin of mass m2 (if this is an injection)")
  parser.add_option("--phi12-inj", dest="phi12_inj", help="injected value of the azimuth angle of (hats2 - hats1) from vecL (if this is an injection)", default=0.)
  parser.add_option("-w", "--waveform", dest="waveform", help="waveform used for recovery")
  parser.add_option("-d", "--debug-plots", dest="debug_plots", help="debug plots")
  parser.add_option("--N_bins", type="int", dest="N_bins", default=201, help="number of bins (default=201)")
  parser.add_option("--dMfbyMf_lim", type="float", dest="dMfbyMf_lim", default=1., help="absolute value of limit for range of dMfbyMf_vec, defined as [-dMfbyMf_lim, +dMfbyMf_lim]")
  parser.add_option("--dchifbychif_lim", type="float", dest="dchifbychif_lim", default=1., help="absolute value of limit for range of dchifbychif_vec, defined as [-dchifbychif_lim, +dchifbychif_lim]")
  parser.add_option("--use_KDE", type="int", dest="MfafKDE", help="use KDE or not after getting samples of Mf, af")
  parser.add_option("-s", "--pepostproc", dest="pepostproc", default= "cbcbayespostproc", help="package used to generate summary webpage to visualize the output [options: 'cbcbayespostproc', 'pesummary']")

  (options, args) = parser.parse_args()
  MfafKDE = options.MfafKDE
  insp_post = options.insp_post
  ring_post = options.ring_post
  imr_post = options.imr_post
  prior_Mfchif_file = options.prior_Mfchif_file
  out_dir = options.out_dir
  fit_formula = options.fit_formula
  debug_plots = options.debug_plots
  pepostproc = options.pepostproc
  if options.insp_fhigh is not None:
    insp_fhigh = float(options.insp_fhigh)
  else:
    print('Inspiral cutoff freq not provided. To have it displayed on the results page, please pass command line option --insp-fhigh.')
    insp_fhigh = np.nan
  if options.ring_flow is not None:
    ring_flow = float(options.ring_flow)
  else:
    print('Ringdown cutoff freq not provided. To have it displayed on the results page, please pass command line option --ring-flow.')
    ring_flow = np.nan
  waveform = options.waveform
  if waveform is None:
    print('Recovery approximant not provided. To have it displayed on the results page, please pass command line option --waveform.')

  N_bins = int(options.N_bins) # Number of grid points along either axis (dMfbyMf, dchifbychif) for computation of the posteriors
  dMfbyMf_lim = float(options.dMfbyMf_lim)
  dchifbychif_lim = float(options.dchifbychif_lim)
  lalinference_datadir = os.getenv('LALINFERENCE_DATADIR')

  # create output directory and copy the script
  os.system('mkdir -p %s' %out_dir)
  os.system('mkdir -p %s/data' %out_dir)
  os.system('mkdir -p %s/img' %out_dir)
  os.system('cp %s %s' %(__file__, out_dir))
  os.system('cp %s %s/'%(os.path.join(lalinference_datadir, 'imrtgr_webpage_templates/*.*'), out_dir))

  # creating file to save the run command
  run_command = open('%s/command.txt'%(out_dir),'w')
  for arg in sys.argv:
    run_command.write('%s\n' %arg)
  run_command.write("\n")
  run_command.write("\n")
  run_command.write("%s"%git_version.verbose_msg)
  run_command.close()

  # creating soft links for PE results
  if pepostproc == 'cbcbayespostproc':
    insp_post_link = insp_post
    ring_post_link = ring_post
    imr_post_link = imr_post
  # posterior samples data file generated by pesummary is stored in
  # <path to summary directory>/samples/<posterior_samples>.dat
  # need to go to the parent directory before creating soft link
  elif pepostproc == 'pesummary':
    insp_post_link = os.path.realpath(os.path.dirname(insp_post))
    ring_post_link = os.path.realpath(os.path.dirname(ring_post))
    imr_post_link = os.path.realpath(os.path.dirname(imr_post))
    os.system("sed -i 's/posplots/home/g' %s/result.html"%out_dir)
  else:
    raise ValueError('Package to generate summary webpage unknown !')

  insp_posplots = os.path.realpath(os.path.dirname(insp_post_link))
  ring_posplots = os.path.realpath(os.path.dirname(ring_post_link))
  imr_posplots = os.path.realpath(os.path.dirname(imr_post_link))

  insp_target = os.path.join(out_dir, 'pe_insp')
  ring_target = os.path.join(out_dir, 'pe_ring')
  imr_target = os.path.join(out_dir, 'pe_imr')

  if insp_posplots != insp_target:
    if os.path.islink(insp_target):
      print('... removing existing link %s'%(insp_target))
      os.system('rm %s'%(insp_target))
    print('... linking %s to %s' %(insp_posplots, insp_target))
    os.system('ln -s %s %s' %(insp_posplots, insp_target))
  if ring_posplots != ring_target:
    if os.path.islink(ring_target):
      print('... removing existing link %s'%(ring_target))
      os.system('rm %s'%(ring_target))
    print('... linking %s to %s' %(ring_posplots, ring_target))
    os.system('ln -s %s %s' %(ring_posplots, ring_target))
  if imr_posplots != imr_target:
    if os.path.islink(imr_target):
      print('... removing existing link %s'%(imr_target))
      os.system('rm %s'%(imr_target))
    print('... linking %s to %s' %(imr_posplots, imr_target))
    os.system('ln -s %s %s' %(imr_posplots, imr_target))

  # read the injection mass parameters if this is an injection
  m1_inj = options.m1_inj
  m2_inj = options.m2_inj
  chi1_inj = options.chi1_inj
  chi2_inj = options.chi2_inj
  chi1z_inj = options.chi1z_inj
  chi2z_inj = options.chi2z_inj
  phi12_inj = options.phi12_inj

  if m1_inj == None or m2_inj == None or chi1_inj == None or chi2_inj == None or chi1z_inj == None or chi2z_inj == None or phi12_inj == None:
    plot_injection_lines = False
  else:
    m1_inj = float(m1_inj)
    m2_inj = float(m2_inj)
    chi1_inj = float(chi1_inj)
    chi2_inj = float(chi2_inj)
    chi1z_inj = float(chi1z_inj)
    chi2z_inj = float(chi2z_inj)
    phi12_inj = float(phi12_inj)
    plot_injection_lines = True
    q_inj = m1_inj/m2_inj
    eta_inj = q_inj/(1.+q_inj)**2.
    Mf_inj, chif_inj = tgr.calc_final_mass_spin(m1_inj, m2_inj, chi1_inj, chi2_inj, chi1z_inj, chi2z_inj, phi12_inj, fit_formula)

  ###############################################################################################
  # Read the posteriors from the inspiral, ringdown and imr PE runs (after post-processing)
  ###############################################################################################

  # parameter names assigned according to the package used to generate summary webpage
  if pepostproc == 'cbcbayespostproc':
    m1_id, m2_id, a1_id, a2_id, a1z_id, a2z_id, phi12_id = 'm1', 'm2', 'a1', 'a2', 'a1z', 'a2z', 'phi12'
  elif pepostproc == 'pesummary':
    m1_id, m2_id, a1_id, a2_id, a1z_id, a2z_id, phi12_id = 'mass_1', 'mass_2', 'a_1', 'a_2', 'spin_1z', 'spin_2z', 'phi_12'
  else:
    raise ValueError('Package to generate summary webpage unknown !')

  # read data from the inspiral posterior file
  insp_data = np.genfromtxt(insp_post, dtype=None, names=True)

  m1_i, m2_i, chi1_i, chi2_i, chi1z_i, chi2z_i = insp_data[m1_id], insp_data[m2_id], insp_data[a1_id], insp_data[a2_id], insp_data[a1z_id], insp_data[a2z_id]
  # if there is phi12 in the posterior, read the values.
  if phi12_id in insp_data.dtype.names:
    phi12_i = insp_data[phi12_id]
  else:
    phi12_i = np.zeros(len(m1_i))
  # compute the final mass and spin
  Mf_i, chif_i = tgr.calc_final_mass_spin(m1_i, m2_i, chi1_i, chi2_i, chi1z_i, chi2z_i, phi12_i, fit_formula)

  # read data from the ringdown posterior file
  ring_data = np.genfromtxt(ring_post, dtype=None, names=True)

  m1_r, m2_r, chi1_r, chi2_r, chi1z_r, chi2z_r = ring_data[m1_id], ring_data[m2_id], ring_data[a1_id], ring_data[a2_id], ring_data[a1z_id], ring_data[a2z_id]
  # if there is phi12 in the posterior, read the values.
  if phi12_id in ring_data.dtype.names:
    phi12_r = ring_data[phi12_id]
  else:
    phi12_r = np.zeros(len(m1_r))
  # compute the final mass and spin
  Mf_r, chif_r = tgr.calc_final_mass_spin(m1_r, m2_r, chi1_r, chi2_r, chi1z_r, chi2z_r, phi12_r, fit_formula)

  # read data from the IMR posterior file
  imr_data = np.genfromtxt(imr_post, dtype=None, names=True)

  m1_imr, m2_imr, chi1_imr, chi2_imr, chi1z_imr, chi2z_imr = imr_data[m1_id], imr_data[m2_id], imr_data[a1_id], imr_data[a2_id], imr_data[a1z_id], imr_data[a2z_id]
  # if there is phi12 in the posterior, read the values.
  if phi12_id in imr_data.dtype.names:
    phi12_imr = imr_data[phi12_id]
  else:
    phi12_imr = np.zeros(len(m1_imr))
  # compute the final mass and spin
  Mf_imr, chif_imr = tgr.calc_final_mass_spin(m1_imr, m2_imr, chi1_imr, chi2_imr, chi1z_imr, chi2z_imr, phi12_imr, fit_formula)


  print('... read posteriors')
  ###############################################################################################

  ###############################################################################################
  # compute the limits of integration for computing delta_Mf and delta_chif
  ###############################################################################################
  Mf_lim = max(abs(np.append(np.append(Mf_i, Mf_r), Mf_imr)))
  chif_lim = max(abs(np.append(np.append(chif_i, chif_r), chif_imr)))

  # the integral used to compute (Delta Mf, Delta af) has limits from -infinity to +infinity. We
  # are approximating this by setting the limits to (-Mf_lim to Mf_lim) and (-chif_lim to chif_lim)
  # where Mf_lim and chif_lim are the max values of Mf and chif where the posteriors have nonzero
  # support. The scipy.signal.correlate2d function requires arguments x_bins and y_bins that need
  # to be symmetric around zero
  Mf_bins = np.linspace(-Mf_lim, Mf_lim, N_bins)
  chif_bins = np.linspace(-chif_lim, chif_lim, N_bins)

  dMf = np.mean(np.diff(Mf_bins))
  dchif = np.mean(np.diff(chif_bins))

  Mf_intp = (Mf_bins[:-1] + Mf_bins[1:])/2.
  chif_intp = (chif_bins[:-1] + chif_bins[1:])/2.


  print('useKDE=',MfafKDE)
  if MfafKDE==1:
    print('replacing lal P(Mfaf) with its KDE pdf')
    M_i,C_i=np.meshgrid(Mf_intp,chif_intp)

    joint_data=np.vstack([Mf_i,chif_i]);kernel=gaussian_kde(joint_data)
    f_i = lambda x,y:kernel.evaluate([x,y])
    print("for inspiral kernel",kernel.integrate_box([-Mf_lim,-chif_lim],[Mf_lim,chif_lim]))
    P_Mfchif_i = np.vectorize(f_i)(M_i,C_i)/kernel.integrate_box([-Mf_lim,-chif_lim],[Mf_lim,chif_lim])

    joint_data=np.vstack([Mf_r,chif_r]);kernel=gaussian_kde(joint_data)#;M_i,C_i=np.meshgrid(Mf_bins,chif_bins)
    f_r = lambda x,y:kernel.evaluate([x,y])
    print("for post-inspiral kernel",kernel.integrate_box([-Mf_lim,-chif_lim],[Mf_lim,chif_lim]))
    P_Mfchif_r = np.vectorize(f_r)(M_i,C_i)/kernel.integrate_box([-Mf_lim,-chif_lim],[Mf_lim,chif_lim])

    joint_data=np.vstack([Mf_imr,chif_imr]);kernel=gaussian_kde(joint_data)#;M_i,C_i=np.meshgrid(Mf_bins,chif_bins)
    f_imr = lambda x,y:kernel.evaluate([x,y])
    print("for imr kernel",kernel.integrate_box([-Mf_lim,-chif_lim],[Mf_lim,chif_lim]))
    P_Mfchif_imr = np.vectorize(f_imr)(M_i,C_i)/kernel.integrate_box([-Mf_lim,-chif_lim],[Mf_lim,chif_lim])


    #rahul: end KDE of Mf,af

  elif MfafKDE==0:
    print('using default samples, NOKDE')
    # compute the 2D posterior distributions for the inspiral, ringodwn and IMR analyses
    P_Mfchif_i, Mf_bins, chif_bins = np.histogram2d(Mf_i, chif_i, bins=(Mf_bins, chif_bins), density=True)
    P_Mfchif_r, Mf_bins, chif_bins = np.histogram2d(Mf_r, chif_r, bins=(Mf_bins, chif_bins), density=True)
    P_Mfchif_imr, Mf_bins, chif_bins = np.histogram2d(Mf_imr, chif_imr, bins=(Mf_bins, chif_bins), density=True)
    # transpose to go from (X,Y) indexing returned by np.histogram2d() to array (i,j) indexing for further
    # computations. From now onwards, different rows (i) correspond to different values of Mf and different
    # columns (j) correspond to different values of chif
    #rahul:Transpose is forbidden if using KDE
    P_Mfchif_i = P_Mfchif_i.T
    P_Mfchif_r = P_Mfchif_r.T
    P_Mfchif_imr = P_Mfchif_imr.T
    print('computed P_Mfchif without using KDE')


  ###############################################################################################


  ###############################################################################################
  # Undo the effect of the prior from the lalinference posterior. Lalinference assumes a        #
  # uniform prior in component masses. We need to assume a uniform prior in Mf, chif              #
  ###############################################################################################

  if prior_Mfchif_file is not None:

    os.system('cp %s %s/data'%(prior_Mfchif_file, out_dir))

    # read the interpolation object, reconstruct the data from the interpolation object
    f = gzip.open(prior_Mfchif_file,'rb')
    P_Mfchif_pr_interp_obj = pickle.load(f)
    P_Mfchif_pr = P_Mfchif_pr_interp_obj(Mf_intp, chif_intp)

    # compute the corrected 2D posteriors in Mf and chif by dividing by the prior distribution
    P_Mfchif_i = P_Mfchif_i/P_Mfchif_pr
    P_Mfchif_r = P_Mfchif_r/P_Mfchif_pr
    P_Mfchif_imr = P_Mfchif_imr/P_Mfchif_pr

    # removing nan's
    P_Mfchif_i[np.isnan(P_Mfchif_i)] = 0.
    P_Mfchif_r[np.isnan(P_Mfchif_r)] = 0.
    P_Mfchif_imr[np.isnan(P_Mfchif_imr)] = 0.

    # removing infinities
    P_Mfchif_i[np.isinf(P_Mfchif_i)] = 0.
    P_Mfchif_r[np.isinf(P_Mfchif_r)] = 0.
    P_Mfchif_imr[np.isinf(P_Mfchif_imr)] = 0.

    print('... computed (prior) corrected posteriors')

  ################################################################################################
  # compute the posterior of (delta_Mf/Mf, delta_chif/chif)
  ################################################################################################

  # compute interpolation objects for the Mf,chif posterior and delta_Mf and delta_chif posterior
  P_Mfchif_i_interp_object = scipy.interpolate.interp2d(Mf_intp, chif_intp, P_Mfchif_i, fill_value=0., bounds_error=False)
  P_Mfchif_r_interp_object = scipy.interpolate.interp2d(Mf_intp, chif_intp, P_Mfchif_r, fill_value=0., bounds_error=False)

  # defining limits of delta_Mf/Mf and delta_chif/chif.
  dMfbyMf_vec = np.linspace(-dMfbyMf_lim, dMfbyMf_lim, N_bins)
  dchifbychif_vec = np.linspace(-dchifbychif_lim, dchifbychif_lim, N_bins)

  # compute the P(dMf/Mf, dchif/chif) by evaluating the integral
  diff_dMfbyMf = np.mean(np.diff(dMfbyMf_vec))
  diff_dchifbychif = np.mean(np.diff(dchifbychif_vec))
  P_dMfbyMf_dchifbychif = np.zeros(shape=(N_bins,N_bins))

  # compute the posterior on the fractional deviation parameters (delta_Mf/Mf, delta_chif/chif).
  # Approximate the integral in Eq.(6) of the document LIGO-P1500185-v5 by a discrete sum
  for i, v2 in enumerate(dchifbychif_vec):
    for j, v1 in enumerate(dMfbyMf_vec):
      P_dMfbyMf_dchifbychif[i,j] = tgr.calc_sum(Mf_intp, chif_intp, v1, v2, P_Mfchif_i_interp_object, P_Mfchif_r_interp_object)

  # normalization
  P_dMfbyMf_dchifbychif /= np.sum(P_dMfbyMf_dchifbychif) * diff_dMfbyMf * diff_dchifbychif

  # Marginalization to one-dimensional joint_posteriors
  P_dMfbyMf = np.sum(P_dMfbyMf_dchifbychif, axis=0) * diff_dchifbychif
  P_dchifbychif = np.sum(P_dMfbyMf_dchifbychif, axis=1) * diff_dMfbyMf

  # compute the confidence region corresponding to the GR value (delta_Mf/Mf = 0, delta_chif/chif = 0).
  # the 'confidence' class is defined on top of this script
  conf_v1v2 = confidence(P_dMfbyMf_dchifbychif)
  gr_height = P_dMfbyMf_dchifbychif[np.argmin(abs(dMfbyMf_vec)), np.argmin(abs(dchifbychif_vec))] # taking value closest to (0,0)
  gr_conf_level = conf_v1v2.level_from_height(gr_height)
  print('... no deviation from GR above %.1f%% confidence level'%(100.*gr_conf_level))

  # creating the parameter table
  param_table = [['Upper cutoff freq for the inspiral analysis: %s Hz'%insp_fhigh],
                ['Lower cutoff freq for the ringdown analysis: %s Hz'%ring_flow],
                ['Waveform approximant: %s'%(waveform)],
                ['Final mass/spin fitting formula: %s'%(fit_formula)],
                ['No deviation from GR above %.1f%% confidence level'%(100.*gr_conf_level)]]
  np.savetxt('%s/summary_table.txt'%(out_dir), np.array(param_table), delimiter='\t', fmt='%s')

  # save results
  np.savetxt(out_dir+'/data/Mfchif.dat.gz', (Mf_bins,chif_bins))
  np.savetxt(out_dir+'/data/P_Mfchif_i.dat.gz', P_Mfchif_i)
  np.savetxt(out_dir+'/data/P_Mfchif_r.dat.gz', P_Mfchif_r)
  np.savetxt(out_dir+'/data/P_Mfchif_imr.dat.gz', P_Mfchif_imr)
  np.savetxt(out_dir+'/data/dMfbyMf_vec.dat.gz', dMfbyMf_vec)
  np.savetxt(out_dir+'/data/dchifbychif_vec.dat.gz', dchifbychif_vec)
  np.savetxt(out_dir+'/data/P_dMfbyMf_dchifbychif.dat.gz', P_dMfbyMf_dchifbychif)
  np.savetxt(out_dir+'/data/P_dMfbyMf.dat.gz', P_dMfbyMf)
  np.savetxt(out_dir+'/data/P_dchifbychif.dat.gz', P_dchifbychif)
  np.savetxt(out_dir+'/data/GR_confidence.txt', [gr_conf_level])

  #########################################################################################

  #########################################################################################
  # plotting
  #########################################################################################
  #inspiral
  P_m1m2_i, m1_bins_i, m2_bins_i = np.histogram2d(m1_i, m2_i, bins=50, density=True)
  P_chi1chi2_i, chi1_bins_i, chi2_bins_i = np.histogram2d(chi1_i, chi2_i, bins=50, density=True)

  P_m1m2_i = P_m1m2_i.T
  P_chi1chi2_i = P_chi1chi2_i.T

  conf_m1m2_i = confidence(P_m1m2_i)
  s1_m1m2_i = conf_m1m2_i.height_from_level(0.68)
  s2_m1m2_i = conf_m1m2_i.height_from_level(0.95)

  conf_chi1chi2_i = confidence(P_chi1chi2_i)
  s1_chi1chi2_i = conf_chi1chi2_i.height_from_level(0.68)
  s2_chi1chi2_i = conf_chi1chi2_i.height_from_level(0.95)

  conf_Mfchif_i = confidence(P_Mfchif_i)
  s1_Mfchif_i = conf_Mfchif_i.height_from_level(0.68)
  s2_Mfchif_i = conf_Mfchif_i.height_from_level(0.95)

  plt.figure(figsize=(5,5))
  plt.pcolormesh(m1_bins_i, m2_bins_i, tgr.gf(P_m1m2_i), cmap='YlOrBr')
  plt.contour(m1_bins_i[:-1], m2_bins_i[:-1], tgr.gf(P_m1m2_i), levels=(s2_m1m2_i,s1_m1m2_i), linewidths=(1,1.5))
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
  plt.contour(chi1_bins_i[:-1], chi2_bins_i[:-1], tgr.gf(P_chi1chi2_i), levels=(s2_chi1chi2_i,s1_chi1chi2_i), linewidths=(1,1.5))
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
  plt.pcolormesh(Mf_bins, chif_bins, tgr.gf(P_Mfchif_i), cmap='YlOrBr')
  plt.contour(Mf_bins[:-1], chif_bins[:-1], tgr.gf(P_Mfchif_i), levels=(s2_Mfchif_i,s1_Mfchif_i), linewidths=(1,1.5))
  if plot_injection_lines == True:
    plt.axvline(x=Mf_inj, ls='--', color='k')
    plt.axhline(y=chif_inj, ls='--', color='k')
  plt.xlabel('$M_f [M_{\odot}]$')
  plt.ylabel('$\chi_f$')
  plt.xlim([min(Mf_i), max(Mf_i)])
  plt.ylim([min(chif_i), max(chif_i)])
  plt.grid()
  plt.savefig('%s/img/inspiral_Mfchif_thumb.png'%(out_dir), dpi=72)
  plt.savefig('%s/img/inspiral_Mfchif.png'%(out_dir), dpi=300)


  #ringdown
  P_m1m2_r, m1_bins_r, m2_bins_r = np.histogram2d(m1_r, m2_r, bins=50, density=True)
  P_chi1chi2_r, chi1_bins_r, chi2_bins_r = np.histogram2d(chi1_r, chi2_r, bins=50, density=True)

  P_m1m2_r = P_m1m2_r.T
  P_chi1chi2_r = P_chi1chi2_r.T

  conf_m1m2_r = confidence(P_m1m2_r)
  s1_m1m2_r = conf_m1m2_r.height_from_level(0.68)
  s2_m1m2_r = conf_m1m2_r.height_from_level(0.95)

  conf_chi1chi2_r = confidence(P_chi1chi2_r)
  s1_chi1chi2_r = conf_chi1chi2_r.height_from_level(0.68)
  s2_chi1chi2_r = conf_chi1chi2_r.height_from_level(0.95)

  conf_Mfchif_r = confidence(P_Mfchif_r)
  s1_Mfchif_r = conf_Mfchif_r.height_from_level(0.68)
  s2_Mfchif_r = conf_Mfchif_r.height_from_level(0.95)

  plt.figure(figsize=(5,5))
  plt.pcolormesh(m1_bins_r, m2_bins_r, tgr.gf(P_m1m2_r), cmap='YlOrBr')
  plt.contour(m1_bins_r[:-1], m2_bins_r[:-1], tgr.gf(P_m1m2_r), levels=(s2_m1m2_r,s1_m1m2_r), linewidths=(1,1.5))
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
  plt.contour(chi1_bins_r[:-1], chi2_bins_r[:-1], tgr.gf(P_chi1chi2_r), levels=(s2_chi1chi2_r,s1_chi1chi2_r), linewidths=(1,1.5))
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
  plt.pcolormesh(Mf_bins, chif_bins, tgr.gf(P_Mfchif_r), cmap='YlOrBr')
  plt.contour(Mf_bins[:-1], chif_bins[:-1], tgr.gf(P_Mfchif_r), levels=(s2_Mfchif_r,s1_Mfchif_r), linewidths=(1,1.5))
  if plot_injection_lines == True:
    plt.axvline(x=Mf_inj, ls='--', color='k')
    plt.axhline(y=chif_inj, ls='--', color='k')
  plt.xlabel('$M_f [M_{\odot}]$')
  plt.ylabel('$\chi_f$')
  plt.xlim([min(Mf_r), max(Mf_r)])
  plt.ylim([min(chif_r), max(chif_r)])
  plt.grid()
  plt.savefig('%s/img/ringdown_Mfchif.png'%(out_dir), dpi=300)
  plt.savefig('%s/img/ringdown_Mfchif_thumb.png'%(out_dir), dpi=72)

  #IMR
  P_m1m2_imr, m1_bins_imr, m2_bins_imr = np.histogram2d(m1_imr, m2_imr, bins=50, density=True)
  P_chi1chi2_imr, chi1_bins_imr, chi2_bins_imr = np.histogram2d(chi1_imr, chi2_imr, bins=50, density=True)

  P_m1m2_imr = P_m1m2_imr.T
  P_chi1chi2_imr = P_chi1chi2_imr.T

  conf_m1m2_imr = confidence(P_m1m2_imr)
  s1_m1m2_imr = conf_m1m2_imr.height_from_level(0.68)
  s2_m1m2_imr = conf_m1m2_imr.height_from_level(0.95)

  conf_chi1chi2_imr = confidence(P_chi1chi2_imr)
  s1_chi1chi2_imr = conf_chi1chi2_imr.height_from_level(0.68)
  s2_chi1chi2_imr = conf_chi1chi2_imr.height_from_level(0.95)

  conf_Mfchif_imr = confidence(P_Mfchif_imr)
  s1_Mfchif_imr = conf_Mfchif_imr.height_from_level(0.68)
  s2_Mfchif_imr = conf_Mfchif_imr.height_from_level(0.95)

  plt.figure(figsize=(5,5))
  plt.pcolormesh(m1_bins_imr, m2_bins_imr, tgr.gf(P_m1m2_imr), cmap='YlOrBr')
  plt.contour(m1_bins_imr[:-1], m2_bins_imr[:-1], tgr.gf(P_m1m2_imr), levels=(s2_m1m2_imr,s1_m1m2_imr), linewidths=(1,1.5))
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
  plt.contour(chi1_bins_imr[:-1], chi2_bins_imr[:-1], tgr.gf(P_chi1chi2_imr), levels=(s2_chi1chi2_imr,s1_chi1chi2_imr), linewidths=(1,1.5))
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
  plt.pcolormesh(Mf_bins, chif_bins, tgr.gf(P_Mfchif_imr), cmap='YlOrBr')
  plt.contour(Mf_bins[:-1], chif_bins[:-1], tgr.gf(P_Mfchif_imr), levels=(s2_Mfchif_imr,s1_Mfchif_imr), linewidths=(1,1.5))
  if plot_injection_lines == True:
    plt.axvline(x=Mf_inj, ls='--', color='k')
    plt.axhline(y=chif_inj, ls='--', color='k')
  plt.xlabel('$M_f [M_{\odot}]$')
  plt.ylabel('$\chi_f$')
  plt.xlim([min(Mf_imr), max(Mf_imr)])
  plt.ylim([min(chif_imr), max(chif_imr)])
  plt.grid()
  plt.savefig('%s/img/imr_Mfchif.png'%(out_dir), dpi=300)
  plt.savefig('%s/img/imr_Mfchif_thumb.png'%(out_dir), dpi=72)

  # IR overlap
  plt.figure(figsize=(5,5))
  CSi = plt.contour(Mf_bins[:-1], chif_bins[:-1], tgr.gf(P_Mfchif_i), levels=(s2_Mfchif_i,s1_Mfchif_i), linewidths=(1,1.5), colors='orange')
  CSr = plt.contour(Mf_bins[:-1], chif_bins[:-1], tgr.gf(P_Mfchif_r), levels=(s2_Mfchif_r,s1_Mfchif_r), linewidths=(1,1.5), colors='red')
  CSimr = plt.contour(Mf_bins[:-1], chif_bins[:-1], tgr.gf(P_Mfchif_imr), levels=(s2_Mfchif_imr,s1_Mfchif_imr), linewidths=(1,1.5), colors='k')
  if plot_injection_lines == True:
    plt.axvline(x=Mf_inj, ls='--', color='k')
    plt.axhline(y=chif_inj, ls='--', color='k')
  plt.xlim([min(np.append(Mf_i, Mf_r)), max(np.append(Mf_i, Mf_r))])
  plt.ylim([min(np.append(chif_i, chif_r)), max(np.append(chif_i, chif_r))])
  plt.xlabel('$M_f~[M_\odot]$')
  plt.ylabel('$\chi_f$')
  plt.grid()

  CSi.levels = np.asarray(CSi.levels)
  CSr.levels = np.asarray(CSr.levels)
  CSimr.levels = np.asarray(CSimr.levels)

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

  ## Label every other level using strings
  plt.clabel(CSi,CSi.levels[::2],inline=True,fmt=fmt_i,fontsize=14, use_clabeltext=True)
  plt.clabel(CSr,CSr.levels[::2],inline=True,fmt=fmt_r,fontsize=14, use_clabeltext=True)
  plt.clabel(CSimr,CSimr.levels[::2],inline=True,fmt=fmt_imr,fontsize=10)

  plt.savefig('%s/img/IMR_overlap.png'%(out_dir), dpi=300)
  plt.savefig('%s/img/IMR_overlap_thumb.png'%(out_dir), dpi=72)

  #(dMf/Mf, dchif/chif)
  conf_v1v2 = confidence(P_dMfbyMf_dchifbychif)
  s1_v1v2 = conf_v1v2.height_from_level(0.68)
  s2_v1v2 = conf_v1v2.height_from_level(0.95)

  conf_v1 = confidence(P_dMfbyMf)
  s1_v1 = conf_v1.height_from_level(0.68)
  s2_v1 = conf_v1.height_from_level(0.95)

  conf_v2 = confidence(P_dchifbychif)
  s1_v2 = conf_v2.height_from_level(0.68)
  s2_v2 = conf_v2.height_from_level(0.95)

  # Calculation of condifence edges
  left1_v1 = min(dMfbyMf_vec[np.where(P_dMfbyMf>=s1_v1)[0]])
  right1_v1 = max(dMfbyMf_vec[np.where(P_dMfbyMf>=s1_v1)[0]])

  left2_v1 = min(dMfbyMf_vec[np.where(P_dMfbyMf>=s2_v1)[0]])
  right2_v1 = max(dMfbyMf_vec[np.where(P_dMfbyMf>=s2_v1)[0]])

  left1_v2 = min(dchifbychif_vec[np.where(P_dchifbychif>s1_v2)[0]])
  right1_v2 = max(dchifbychif_vec[np.where(P_dchifbychif>s1_v2)[0]])

  left2_v2 = min(dchifbychif_vec[np.where(P_dchifbychif>s2_v2)[0]])
  right2_v2 = max(dchifbychif_vec[np.where(P_dchifbychif>s2_v2)[0]])

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
  plt.pcolormesh(dMfbyMf_vec,dchifbychif_vec,P_dMfbyMf_dchifbychif, cmap='YlOrBr')
  plt.contour(dMfbyMf_vec,dchifbychif_vec,tgr.gf(P_dMfbyMf_dchifbychif), levels=(s2_v1v2,s1_v1v2), linewidths=(1,1.5))
  plt.plot(0, 0, 'k+', ms=12, mew=2)
  plt.xlabel('$\Delta M_f/M_f$')
  plt.ylabel('$\Delta \chi_f/\chi_f$')
  plt.xlim([-dMfbyMf_lim,dMfbyMf_lim])
  plt.ylim([-dchifbychif_lim,dchifbychif_lim])
  plt.grid()

  plt.subplot2grid((3,3), (1,2), rowspan=2)
  plt.plot(P_dchifbychif, dchifbychif_vec,'k', lw=1)
  plt.axhline(y=left1_v2, color='c', lw=0.5, ls='-.')
  plt.axhline(y=right1_v2, color='c', lw=0.5, ls='-.')
  plt.axhline(y=left2_v2, color='b', lw=0.5, ls='-.')
  plt.axhline(y=right2_v2, color='b', lw=0.5, ls='-.')
  #plt.ylabel('$\Delta \chi_f/\chi_f$')
  plt.xlabel('$P(\Delta \chi_f/\chi_f)$')
  #plt.grid()

  plt.savefig('%s/img/dMfbyMfdchifbychif.png' %(out_dir), dpi=300)
  plt.savefig('%s/img/dMfbyMfdchifbychif_thumb.png' %(out_dir), dpi=72)

  print('... made summary plots')

  print('... completed in %f seconds' %(time.time()-start_time))
  #########################################################################################

