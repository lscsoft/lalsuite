"""
Core python module for the IMR consistency test 

P. Ajith, Abhirup Ghosh, 2015-09-15

$Id:$
"""

import numpy as np 
import scipy.ndimage.filters as filter
from multiprocessing import Pool
from functools import partial
import time 
import nrutils as nr

""" calculate the mass and spin of the final black hole using an NR-inspired fitting formula """
def calc_final_mass_spin(m1, m2, chi1, chi2, fit_formula):
  
  if fit_formula == 'nospin_Pan2011':
    chif = nr.bbh_final_spin_non_spinning_Panetal(m1, m2)
    mf = nr.bbh_final_mass_non_spinning_Panetal(m1, m2)
  elif fit_formula == 'nonprecspin_Healy2014':
    chif = nr.bbh_final_spin_non_precessing_Healyetal(m1, m2, chi1, chi2)
    mf = nr.bbh_final_mass_non_precessing_Healyetal(m1, m2, chi1, chi2, chif)
  elif fit_formula == 'nonprecspin_Husa2015':
    chif = nr.bbh_final_spin_non_precessing_Husaetal(m1, m2, chi1, chi2)
    mf = nr.bbh_final_mass_non_precessing_Husaetal(m1, m2, chi1, chi2)
  else:
    raise ValueError("unknown spin fit formula")
  return mf, chif

""" compute the integrand of P(dMf/Mf, dchif/chif). """
def P_integrand(chif, Mf, v1, v2, P_dMfdchif_interp_object, P_Mfchif_imr_interp_object):

  Mf_mat, chif_mat = np.meshgrid(Mf, chif)

  # Create dMf and dchif vectors corresponding to the given v1 and v2. These vectors have to be 
  # monotonically increasing in order to evaluate the interpolated prob densities. Hence, for 
  # v1, v2 < 0, flip them, evaluate the prob density (in column or row) and flip it back 
  dMf = v1*Mf
  dchif = v2*chif

  if v1 < 0.:
    dMf = np.flipud(dMf)
  if v2 < 0.:
    dchif = np.flipud(dchif)
  P_delta = P_dMfdchif_interp_object(dMf, dchif)

  if v1 < 0.:
    P_delta = np.fliplr(P_delta)
  if v2 < 0.:
    P_delta = np.flipud(P_delta)

  P_imr = P_Mfchif_imr_interp_object(Mf, chif)
  return P_imr*P_delta*abs(Mf_mat*chif_mat), P_imr, P_delta

""" compute P(dMf/Mf, dchif/chif). """
def calc_sum(Mf, chif, v1, v2, P_dMfdchif_interp_object, P_Mfchif_imr_interp_object):

  Pintg, P_imr, P_delta = P_integrand(chif, Mf, v1, v2, P_dMfdchif_interp_object, P_Mfchif_imr_interp_object)
  return np.sum(Pintg)

""" gaussian filter of histogram """
def gf(P):
  return filter.gaussian_filter(P, sigma=2.0)

""" generate prior samples in (Mf, chif) assuming uniform prior in component masses """
def calc_Mfchif_prior_samples(comp_mass_prior_min, comp_mass_prior_max, comp_spin_min, comp_spin_max, Mf_bins, chif_bins, fit_formula, spin_angle_dist, N_sampl, thread):

  # generate random samples of the prior uniform in component masses 
  m1_pr = np.random.uniform(comp_mass_prior_min, comp_mass_prior_max, N_sampl)
  m2_pr = np.random.uniform(comp_mass_prior_min, comp_mass_prior_max, N_sampl)

  # generate samples of component spins, for non-precessing spins (aligned/anti-aligned with the orb ang momentum)  
  if spin_angle_dist == 'aligned':
    chi1_pr = np.random.uniform(comp_spin_min, comp_spin_max, N_sampl)
    chi2_pr = np.random.uniform(comp_spin_min, comp_spin_max, N_sampl)
  # generate samples of component spins, for uninform spin magnitudes and isotropic spin directions (uniform in cos_tilt_angle) 
  elif spin_angle_dist == 'isotropic':
    chi1_pr = np.random.uniform(comp_spin_min, comp_spin_max, N_sampl)*np.random.uniform(-1., 1., N_sampl)
    chi2_pr = np.random.uniform(comp_spin_min, comp_spin_max, N_sampl)*np.random.uniform(-1., 1., N_sampl)

  # return the corrresponding samples of prior in Mf, chif 
  return calc_final_mass_spin(m1_pr, m2_pr, chi1_pr, chi2_pr, fit_formula)

""" calculate the prior distribution in (Mf, chif) assuming uniform prior in component masses """
def calc_Mfchif_prior(comp_mass_prior_min, comp_mass_prior_max, comp_spin_min, comp_spin_max, Mf_bins, chif_bins, fit_formula, spin_angle_dist, N_sampl, num_threads):

  # check inputs 
  if comp_mass_prior_min < 1. or comp_mass_prior_min > 1000.: raise ValueError("comp_mass_prior_min should be in the interval [1, 1000]")
  if comp_mass_prior_max < 1. or comp_mass_prior_max > 1000.: raise ValueError("comp_mass_prior_max should be in the interval [1, 1000]")
  if comp_mass_prior_max <= comp_mass_prior_min : raise ValueError("comp_mass_prior_max should be greater than comp_mass_prior_min")
  if N_sampl < 1: raise ValueError("N_sampl should be greater than 1")
  if comp_spin_max < comp_spin_min: raise ValueError("comp_spin_max should be greater than comp_spin_min")
  if spin_angle_dist == 'aligned':
    if comp_spin_min < -1. or comp_spin_min > 1.: raise ValueError("comp_spin_min should be in the interval [-1, 1] for the case of aligned spin distributions")
    if comp_spin_max < -1. or comp_spin_max > 1.: raise ValueError("comp_spin_max should be in the interval [-1, 1] for the case of aligned spin distributions")
  elif spin_angle_dist == 'isotropic':
    if comp_spin_min < 0. or comp_spin_min > 1.: raise ValueError("comp_spin_min should be in the interval [0, 1] for the case of isotrpic spin distributions")
    if comp_spin_max < 0. or comp_spin_max > 1.: raise ValueError("comp_spin_max should be in the interval [0, 1] for the case of isotrpic spin distributions")
  else:
    raise ValueError("spin_angle_dist should be 'aligned' or 'isotropic'")


  # generate samples in Mf and chif corresponding to uniform samples in m1, m2. Parallelise the calculation over multiple threads 
  thread_vec = range(num_threads)
  N_sampl = int(N_sampl/num_threads)
  p = Pool(num_threads)
  func = partial(calc_Mfchif_prior_samples, comp_mass_prior_min, comp_mass_prior_max, comp_spin_min, comp_spin_max, Mf_bins, chif_bins, fit_formula, spin_angle_dist, N_sampl)
  final_params = p.map(func, thread_vec)
  p.close()
  p.join()
  
  final_params = np.array(final_params)
  final_params.reshape(N_sampl*num_threads, 2)
  Mf_pr = np.ravel(final_params[:,0])
  chif_pr = np.ravel(final_params[:,1])

  # compute the 2D prior distribution in Mf and chif 
  P_Mfchif_pr, Mf_bins, chif_bins = np.histogram2d(Mf_pr, chif_pr, bins=(Mf_bins, chif_bins), normed=True)

  return P_Mfchif_pr.T
