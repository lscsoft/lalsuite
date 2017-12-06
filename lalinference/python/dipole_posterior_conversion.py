"""
Transforms the posterior on a generic -1PN coefficient (dipolecoeff) from a LALInference run into other, less generic dipole radiation constraints
Possible prior effects in the posterior transformation are taken into account.

Maria Haney, 2016

$Id:$
"""
import numpy as np
import scipy.optimize as so
from optparse import OptionParser
from scipy import random
from scipy import interpolate as interp

# Functions for converting the posterior of generic -1PN coefficient 'dipolecoeff' into different dipole radiation constraints

def convert_dipolecoeff_into_dipolepar_Horbatsch_etal(m1, m2, dipolecoeff):
    """
    Converts the posterior of a generic PN coefficient at -1PN phase order (dipolecoeff) into that of the free parameter mu (mu_dipole), parametrizing certain 'scalarization'-type dipole radiation and quantifying the aquired scalar charge of the BBH.
    See: Section 2.4 of M.W. Horbatsch, C.P. Burgess, Cosmic black-hole hair growth and quasar OJ287, JPAC 5, 010 (2012), arXiv:1111.4009
    
    Input parameters
    ----------
    m1, m2 : component masses
    dipolecoeff : generic -1PN coefficient
    
    Returns
    ------
    scalar dipole radiation parameter, 'mu_dipole'
    """

    m1 = np.vectorize(float)(np.array(m1))
    m2 = np.vectorize(float)(np.array(m2))
    dchi = np.vectorize(float)(np.array(dipolecoeff))

    m = m1 + m2
    delta = (m1-m2)/m

    return np.sqrt(21./10.) * np.sqrt(-1.*dchi) / (m*delta)

def convert_dipolecoeff_into_dipolepar_Barausse_etal(dipolecoeff):
    """
    Converts the posterior of a generic PN coefficient at -1PN phase order (dipolecoeff) into that of the free parameter B (B_dipole), parametrizing theory-agnostic dipole radiation corrections to the flux.
    See: E. Barausse, N. Yunes, K. Chamberlain, Theory-Agnostic Constraints on Black-Hole Dipole Radiation with Multi-Band Gravitational-Wave Astrophysics, Phys. Rev. Lett. 116, 241104 (2016), arXiv:1603.04075

    Input parameters
    ----------
    dipolecoeff : generic -1PN coefficient

    Returns
    ------
    theory-agnostic dipole flux parameter, 'B_dipole'
    """
    dchi = np.vectorize(float)(np.array(dipolecoeff))

    return -(7./4.) * dchi

def convert_dipolecoeff_into_scalchargediff(dipolecoeff):
    """
    For comparison with Ben Farr's results
    Documentation goes here...
    
    Input parameters
    ----------
    dipolecoeff : generic -1PN coefficient
    
    Returns
    ------
    BBH scalar charge difference, 'deltaq'
    """
    dchi = np.vectorize(float)(np.array(dipolecoeff))

    return np.sqrt(42./5.) * np.sqrt(-1.*dchi)

# Utility function for selecting only those samples for which dipolecoeff <= 0

def physical_samples_only(m1, m2, dipolecoeff):

    m1_new = []
    m2_new = []
    dchi_new = []
    m1_val = 0.
    m2_val = 0.
    dchi_val = 0.

    for i in range(len(dipolecoeff)):
       if dipolecoeff[i] <= 0:
          dchi_val = dipolecoeff[i]
          m1_val = m1[i]
          m2_val = m2[i]
          dchi_new.append(dchi_val)
          m1_new.append(m1_val)
          m2_new.append(m2_val)

    return m1_new, m2_new, dchi_new




if __name__ == '__main__':

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rc
    matplotlib.rc('text.latex', preamble = '\usepackage{txfonts}')

    rc('text', usetex=True)
    rc('mathtext', default='rm')

    parser = OptionParser()
    parser.add_option("--input-file", dest="input_file", help="file containing the posterior samples from a LALInference run")
    parser.add_option("--out-dir", dest="out_dir", help="output directory")
    parser.add_option("--conv-formula", dest="conv_formula", default="dipolescalar", help="formula for chosen dipole parameter conversions (options: 'dipolescalar', 'dipolefluxpar', 'chargediff')")
    parser.add_option("--lower-range", type="float", dest="lower_range", default=0., help="lower limit of range for P(mu) plots (default=0)")
    parser.add_option("--upper-range", type="float", dest="upper_range", default=0.1, help="upper limit of range for P(mu) plots (default=0.1)")
    parser.add_option("--N-bins", type="int", dest="N_bins", default=201, help="number of bins (default=201)")
    parser.add_option("--N-samp", type="int", dest="N_samp", default=10000000, help="number of sampling points for prior computation")
    parser.add_option("--comp-mass-min", type="float", dest="comp_mass_min", default=10., help="lower limit of component mass prior (default=10)")
    parser.add_option("--comp-mass-max", type="float", dest="comp_mass_max", default=80., help="upper limit of component mass prior (default=80)")
    (options, args) = parser.parse_args()

    input_file = options.input_file
    out_dir = options.out_dir
    conv_formula = options.conv_formula
    lower_range = options.lower_range
    upper_range = options.upper_range
    N_bins = int(options.N_bins)
    N_samp = int(options.N_samp)
    comp_mass_min = options.comp_mass_min
    comp_mass_max = options.comp_mass_max

    # set up strings for file names and plotting
    if conv_formula == 'dipolescalar':
      name_strs = 'mu_dipole'
      plot_strs = '\mu_{dipole}'
    elif conv_formula == 'chargediff':
      name_strs = 'delta_q'
      plot_strs = '\delta q'
    elif conv_formula == 'dipolefluxpar':
      name_strs = 'B_dipole'
      plot_strs = 'B_{dipole}'
    else:
      raise ValueError("unknown dipole conversion formula")

    # uniform priors in component masses and -1PN dchi
    mass_range = (comp_mass_min, comp_mass_max)
    dchi_range = (-1., 1.)

    m1_pr =  random.uniform(*mass_range, size=N_samp)
    m2_pr =  random.uniform(*mass_range, size=N_samp)
    dchi_pr =  random.uniform(*dchi_range, size=N_samp)

    # for the conversion of m1, m2, dchi priors into mu and delta_q priors, only those samples for which dchi<=0 are physically meaningful (since both mu and delta_q ~ sqrt(-dchi) )
    m1_newpr, m2_newpr, dchi_newpr = physical_samples_only(m1_pr, m2_pr, dchi_pr)

    m1_newpr = np.vectorize(float)(np.array(m1_newpr))
    m2_newpr = np.vectorize(float)(np.array(m2_newpr))
    dchi_newpr = np.vectorize(float)(np.array(dchi_newpr))

    # convert uniform priors in m1, m2, dchi into a prior for a less generic -1PN parameter
    if conv_formula == 'dipolescalar':
      dipolepar_pr = convert_dipolecoeff_into_dipolepar_Horbatsch_etal(m1_newpr, m2_newpr, dchi_newpr)
      bin_range_pr = (lower_range,upper_range)
    elif conv_formula == 'chargediff':
      dipolepar_pr = convert_dipolecoeff_into_scalchargediff(dchi_newpr)
      bin_range_pr = (min(dipolepar_pr), max(dipolepar_pr))
    elif conv_formula == 'dipolefluxpar':
      dipolepar_pr = convert_dipolecoeff_into_dipolepar_Barausse_etal(dchi_pr)
      bin_range_pr = (min(dipolepar_pr), max(dipolepar_pr))
    else:
      raise ValueError("unknown dipole conversion formula")

    # plot prior distribution
    plt.figure()
    plt.hist(dipolepar_pr, range=bin_range_pr, bins=N_bins, normed=True, facecolor='grey')
    plt.xlabel(r'$%s$'%plot_strs)
    plt.ylabel(r'$Prior\,Probability\,Density$')
    plt.grid()
    plt.savefig(out_dir+'prior_%s.png'%name_strs)

    print '... computed transformed prior'

    # read posterior samples from LALInference run
    samp_data = np.genfromtxt(input_file, dtype=None, names=True)

    dchi_arr = samp_data['dipolecoeff']
    m1_arr = samp_data['m1']
    m2_arr = samp_data['m2']

    # for the conversion of m1, m2, dchi posterior samples into mu and delta_q posteriors, only those samples for which dchi<=0 are physically meaningful (since both mu and delta_q ~ sqrt(-dchi) )
    m1, m2, dchi = physical_samples_only(m1_arr, m2_arr, dchi_arr)

    m1 = np.vectorize(float)(np.array(m1))
    m2 = np.vectorize(float)(np.array(m2))
    dchi = np.vectorize(float)(np.array(dchi))

    # convert posterior samples in m1, m2, dchi into a posterior for a less generic -1PN parameter
    if conv_formula == 'dipolescalar':
      dipolepar = convert_dipolecoeff_into_dipolepar_Horbatsch_etal(m1, m2, dchi)
      bin_range = (lower_range,upper_range)
    elif conv_formula == 'chargediff':
      dipolepar = convert_dipolecoeff_into_scalchargediff(dchi)
      bin_range = (min(dipolepar), max(dipolepar))
    elif conv_formula == 'dipolefluxpar':
      dipolepar = convert_dipolecoeff_into_dipolepar_Barausse_etal(dchi_arr)
      bin_range = (min(dipolepar), max(dipolepar))
    else:
      raise ValueError("unknown dipole conversion formula")

    print '... computed transformed posterior'

    # plot posterior distribution
    plt.figure()
    plt.hist(dipolepar, range=bin_range, bins=N_bins, normed=True, facecolor='grey')
    plt.xlabel(r'$%s$'%plot_strs)
    plt.ylabel(r'$Probability\,Density$')
    plt.grid()
    plt.savefig(out_dir+'%s.png'%name_strs)

    # create interpolation object of prior distribution to undo prior effects in the posterior transformation
    dipolepar_range = bin_range  

    dipolepar_bins = np.linspace(min(dipolepar_range), max(dipolepar_range), N_bins)
    P_dipolepar_pr, dipolepar_bins = np.histogram(dipolepar_pr, bins=dipolepar_bins, normed=True)
    dipolepar_bins_intp = (dipolepar_bins[:-1] + dipolepar_bins[1:])/2.
    P_dipolepar_pr_interp_obj = interp.interp1d(dipolepar_bins_intp, P_dipolepar_pr, fill_value=0., bounds_error=False)

    # divide posterior distribution by interpolation object to compute prior-corrected posterior
    P_dipolepar, dipolepar_bins = np.histogram(dipolepar, bins=dipolepar_bins, normed=True)
    P_dipolepar_corr = P_dipolepar/P_dipolepar_pr_interp_obj(dipolepar_bins_intp)

    # check for nans and infinities
    num_nan = len(np.where(np.isnan(P_dipolepar_corr))[0])
    num_inf = len(np.where(np.isinf(P_dipolepar_corr))[0])

    # remove any nans and infinities (N-samp warning)
    if num_nan != 0.:
      P_dipolepar_corr[np.isnan(P_dipolepar_corr)] = 0.
      print "%d nans removed from prior-corrected posterior; insufficient N-samp"%num_nan
    if num_inf != 0.:
      P_dipolepar_corr[np.isinf(P_dipolepar_corr)] = 0.
      print "%d infinities removed from prior-corrected posterior; insufficient N-samp"%num_inf

    print '... computed prior-corrected posterior'

    #plot prior-corrected posterior distribution
    plt.figure()
    width = 1.0 * (dipolepar_bins[1] - dipolepar_bins[0])
    center = (dipolepar_bins[:-1] + dipolepar_bins[1:])/2.
    plt.bar(center, P_dipolepar_corr, align='center', width=width, facecolor='grey')
    plt.xlabel(r'$%s$'%plot_strs)
    plt.ylabel(r'$Prior-corrected\,Probability\,Density$')
    plt.grid()
    plt.savefig(out_dir+'%s_priorcorr.png'%name_strs)
