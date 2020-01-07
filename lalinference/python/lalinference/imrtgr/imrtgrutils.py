"""
Core python module for the IMR consistency test

P. Ajith, Abhirup Ghosh, 2015-09-15

$Id:$
"""

import numpy as np
import scipy.ndimage.filters as filter
from multiprocessing import Pool
from functools import partial
from . import nrutils as nr

""" calculate the mass and spin of the final black hole using an NR-inspired fitting formula """
def calc_final_mass_spin(m1, m2, chi1, chi2, chi1z, chi2z, phi12, fit_formula):
  """ Calculate the mass and spin of the final black hole using an NR-inspired fitting formula.

  inputs:
  m1, m2: initial masses
  chi1, chi2: initial spin magnitudes
  chi1z, chi2z: z-components of the initial spins
  phi12: in-plane angle between initial spins
  fit_formula: fitting formula to be used for the calculation of final mass/spin

  output:
  mf, chif: final mass, final spin
  """

  tilt1 = np.arccos(chi1z/chi1)
  tilt2 = np.arccos(chi2z/chi2)

  if fit_formula == 'nospin_Pan2011':
    chif = nr.bbh_final_spin_non_spinning_Panetal(m1, m2)
    mf = nr.bbh_final_mass_non_spinning_Panetal(m1, m2)
  elif fit_formula == 'nonprecspin_Healy2014':
    chif = nr.bbh_final_spin_non_precessing_Healyetal(m1, m2, chi1z, chi2z, version="2014")
    mf = nr.bbh_final_mass_non_precessing_Healyetal(m1, m2, chi1z, chi2z, version="2014", chif=chif)
  elif fit_formula == 'nonprecspin_Husa2015':
    chif = nr.bbh_final_spin_non_precessing_Husaetal(m1, m2, chi1z, chi2z)
    mf = nr.bbh_final_mass_non_precessing_Husaetal(m1, m2, chi1z, chi2z)
  elif fit_formula == 'bbh_average_fits_precessing':
    mf_fits = ["UIB2016", "HL2016"]
    chif_fits = ["UIB2016", "HBR2016", "HL2016"]
    mf = nr.bbh_average_fits_precessing(m1, m2, chi1, chi2, tilt1, tilt2, 0, "Mf", mf_fits)
    chif = nr.bbh_average_fits_precessing(m1, m2, chi1, chi2, tilt1, tilt2, phi12, "af", chif_fits)
  else:
    raise ValueError("unknown spin fit formula")
  return mf, chif

""" compute the integrand of P(dMf/Mfbar, dchif/chifbar). """
def P_integrand(chif, Mf, v1, v2, P_Mfchif_i_interp_object, P_Mfchif_r_interp_object):

  """ Compute the integrand of P(dMf/Mfbar, dchif/chifbar).

  inputs:
  chif: vector of values of final spin
  Mf: vector of values of final mass
  v1: dMf/Mfbar value
  v2: dchif/chifbar value
  P_Mfchif_i_interp_object: interpolation function of P_i(Mf, chif)
  P_Mfchif_r_interp_object: interpolation function of P_r(Mf, chif)

  output: integrand of P(dMf/Mfbar, dchif/chifbar)
  """


  Mf_mat, chif_mat = np.meshgrid(Mf, chif)

  # Create dMf and dchif vectors corresponding to the given v1 and v2. These vectors have to be
  # monotonically increasing in order to evaluate the interpolated prob densities. Hence, for
  # v1, v2 < 0, flip them, evaluate the prob density (in column or row) and flip it back
  dMf_i = (1.+v1/2.)*Mf
  dchif_i = (1.+v2/2.)*chif

  dMf_r = (1.-v1/2.)*Mf
  dchif_r = (1.-v2/2.)*chif

  if (1.+v1/2.) < 0.:
    dMf_i = np.flipud(dMf_i)
  if (1.+v2/2.) < 0.:
    dchif_i = np.flipud(dchif_i)
  P_i = P_Mfchif_i_interp_object(dMf_i, dchif_i)

  if (1.+v1/2.) < 0.:
    P_i = np.fliplr(P_i)
  if (1.+v2/2.) < 0.:
    P_i = np.flipud(P_i)

  if (1.-v1/2.) < 0.:
    dMf_r = np.flipud(dMf_r)
  if (1.-v2/2.) < 0.:
    dchif_r = np.flipud(dchif_r)
  P_r = P_Mfchif_r_interp_object(dMf_r, dchif_r)

  if (1.-v1/2.) < 0.:
    P_r = np.fliplr(P_r)
  if (1.-v2/2.) < 0.:
    P_r = np.flipud(P_r)

  return P_i*P_r*abs(Mf_mat*chif_mat), P_i, P_r

""" compute P(dMf/Mfbar, dchif/chifbar). """
def calc_sum(Mf, chif, v1, v2, P_Mfchif_i_interp_object, P_Mfchif_r_interp_object):

  Pintg, P_i, P_r = P_integrand(chif, Mf, v1, v2, P_Mfchif_i_interp_object, P_Mfchif_r_interp_object)
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
  P_Mfchif_pr, Mf_bins, chif_bins = np.histogram2d(Mf_pr, chif_pr, bins=(Mf_bins, chif_bins), density=True)

  return P_Mfchif_pr.T
