"""
Various fitting formulas provided by numerical relativity

Archisman Ghosh, Nathan K. Johnson-McDaniel, P. Ajith, 2015-04-09
"""

import numpy as np
import scipy.optimize as so
try:
    import lal
except ImportError:
    print('Cannot import lal SWIG bindings')

# Functions to check that input values are physical

def _check_m(m1,m2):
    if np.any(m1<=0):
        raise ValueError("m1 has to be positive")
    if np.any(m2<=0):
        raise ValueError("m2 has to be positive")

def _check_chi(chi1, chi2):
    if np.any(chi1>1):
      raise ValueError("chi1 has to be <= 1")
    if np.any(chi2>1):
      raise ValueError("chi2 has to be <= 1")
    if np.any(chi1<0):
      raise ValueError("chi1 has to be nonnegative")
    if np.any(chi2<0):
      raise ValueError("chi2 has to be nonnegative")

def _check_chi_aligned(chi1, chi2):
    if np.any(abs(chi1)>1):
      raise ValueError("chi1 has to be <= 1")
    if np.any(abs(chi2)>1):
      raise ValueError("chi2 has to be <= 1")

def _check_mchi(m1,m2,chi1,chi2):
    _check_m(m1,m2)
    _check_chi(chi1,chi2)

def _check_mchi_aligned(m1,m2,chi1,chi2):
    _check_m(m1,m2)
    _check_chi_aligned(chi1,chi2)

###########################
# Final mass and spin fits
###########################

# list of final mass and spin fit names
class bbh_final_state_fits:
    pan2011   = "Pan2011"   # Pan et al.             [Phys Rev D 84, 124052 (2011)] (mass+spin, non-spinning)
    hlz2014   = "HLZ2014"   # Healy et al.           [Phys Rev D 90, 104004 (2014)] (mass+spin, aligned)
    phenD     = "PhenomD"   # Husa et al.            [Phys Rev D 93, 044006 (2016)] (mass+spin, aligned)
    uib2016v1 = "UIB2016v1" # Jimenez-Forteza et al. [arXiv:1611.00332v1]           (mass+spin, aligned)
    uib2016   = "UIB2016"   # Jimenez-Forteza et al. [LIGO-P1600270-v4 / arXiv:1611.00332v2]           (mass+spin, aligned)
    hbr2016   = "HBR2016"   # Hofmann et al.         [ApJL 825:L19 (2016)]          (spin only, precessing)
    hl2016    = "HL2016"    # Healy and Lousto       [arXiv:1610.09713]             (mass+spin, aligned)

# list of Kerr truncation behaviours
class bbh_Kerr_trunc_opts:
    trunc        = "TRUNCATE"       # truncate fit value to 1.0, with an info message
    trunc_silent = "TRUNCATESILENT" # truncate fit value to 1.0, without info message
    keep         = "KEEP"           # keep fit value, even if superextremal, but still give the info message
    ign          = "IGNORE"         # keep fit value, even if superextremal, without info message
    err          = "ERROR"          # abort with an error if fit value is superextremal

# Function to truncate final spins at Kerr limit (with various options)

def _truncate_at_Kerr_limit(chif, behavior, fitname="this"):
    if behavior not in list(bbh_Kerr_trunc_opts.__dict__.values()):
      raise ValueError("Unknown user option '%s' for Kerr truncation behavior." % behavior)
    if behavior!=bbh_Kerr_trunc_opts.ign and np.any(abs(chif)>1.0):
      if isinstance(chif,(list,np.ndarray)):
        idx_over = np.where(abs(chif)>1.0)
        if behavior==bbh_Kerr_trunc_opts.trunc or behavior==bbh_Kerr_trunc_opts.trunc_silent:
          chif_trunc = np.sign(chif[idx_over])*1.0
          if behavior==bbh_Kerr_trunc_opts.trunc:
            print("Truncating %d excessive chif values from %s fit to Kerr limit of +-1.0" % (np.size(idx_over), fitname))
          chif[idx_over] = chif_trunc
        elif behavior==bbh_Kerr_trunc_opts.keep:
          print("Note: %s fit predicts %d chif values in excess of the Kerr limit." % (fitname, np.size(idx_over)))
        elif behavior==bbh_Kerr_trunc_opts.err:
          raise ValueError("%s fit predicts %d chif values in excess of the Kerr limit." % (fitname, np.size(idx_over)))
      else:
        if behavior==bbh_Kerr_trunc_opts.trunc or behavior==bbh_Kerr_trunc_opts.trunc_silent:
          chif_trunc = np.sign(chif)*1.0
          if behavior==bbh_Kerr_trunc_opts.trunc:
            print("Truncating excessive chif of %f from %s fit to Kerr limit of %f" % (chif, fitname, chif_trunc))
          chif = chif_trunc
        elif behavior==bbh_Kerr_trunc_opts.keep:
          print("Note: %s fit predicts chif of %f in excess of the Kerr limit." % (fitname, chif))
        elif behavior==bbh_Kerr_trunc_opts.err:
          raise ValueError("%s fit predicts chif of %f in excess of the Kerr limit." % (fitname, chif))
    return chif

# Final mass and spin functions

def bbh_final_mass_non_spinning_Panetal(m1, m2):
    """
    Calculate the mass of the final BH resulting from the merger of two non-spinning black holes using eqn. 29a of from Pan et al, Phys Rev D 84, 124052 (2011).

    Parameters
    ----------
    m1, m2 : component masses

    Returns
    ------
    final mass, mf
    """
    m1 = np.vectorize(float)(np.array(m1))
    m2 = np.vectorize(float)(np.array(m2))

    m = m1 + m2
    eta = m1*m2/(m1+m2)**2.
    return m*(1. + (np.sqrt(8./9.)-1.)*eta - 0.4333*(eta**2.) - (0.4392*(eta**3)))

def bbh_final_spin_non_spinning_Panetal(m1, m2):
    """
    Calculate the spin of the final BH resulting from the merger of two non-spinning black holes using eqn. 29b of Pan et al, Phys Rev D 84, 124052 (2011).

    Parameters
    ----------
    m1, m2 : component masses

    Returns
    ------
    final dimensionless spin, chif
    """
    m1 = np.vectorize(float)(np.array(m1))
    m2 = np.vectorize(float)(np.array(m2))

    eta = m1*m2/(m1+m2)**2.
    return np.sqrt(12.)*eta - 3.871*(eta**2.) + 4.028*(eta**3)

def calc_isco_radius(a):
    """
    Calculate the ISCO radius of a Kerr BH as a function of the Kerr parameter using eqns. 2.5 and 2.8 from Ori and Thorne, Phys Rev D 62, 24022 (2000)

    Parameters
    ----------
    a : Kerr parameter

    Returns
    -------
    ISCO radius
    """

    a = np.minimum(np.array(a),1.) # Only consider a <=1, to avoid numerical problems

    # Ref. Eq. (2.5) of Ori, Thorne Phys Rev D 62 124022 (2000)
    z1 = 1.+(1.-a**2.)**(1./3)*((1.+a)**(1./3) + (1.-a)**(1./3))
    z2 = np.sqrt(3.*a**2 + z1**2)
    a_sign = np.sign(a)
    return 3+z2 - np.sqrt((3.-z1)*(3.+z1+2.*z2))*a_sign


def _RIT_setup(m1, m2, chi1, chi2):
    """ Internal function to compute the combinations of masses and spins that are used in the RIT fits """

    _check_mchi_aligned(m1, m2, chi1, chi2)

    m = m1+m2
    eta = m1*m2/m/m
    delta_m = (m1 - m2)/m

    S = (m1*m1*chi1 + m2*m2*chi2)/m**2 # symmetric spin (dimensionless -- called \tilde{S} in the paper)
    Delta = (m2*chi2 - m1*chi1)/m # antisymmetric spin (dimensionless -- called tilde{Delta} in the paper

    return eta, delta_m, S, Delta

def _RIT_symm_express(eta, delta_m, S, Delta, coeffs):
    """
    Internal function to return the basic even-in-interchange-of-particles expression used in the RIT fits for the final mass, final spin, and peak luminosity

    The coefficients in the expression are passed as a vector as [0 1 2a 2b 2c 2d 3a 3b 3c 3d 4a 4b 4c 4d 4e 4f 4g 4h 4i] (giving the subscripts).
    """

    S2 = S*S

    dm2 = delta_m*delta_m

    Deltadm = Delta*delta_m

    Delta2 = Delta*Delta

    return 16.*eta*eta*(coeffs[0] + coeffs[1]*S + coeffs[2]*Deltadm + coeffs[3]*S2 + coeffs[4]*Delta2 \
        + coeffs[5]*dm2 + coeffs[6]*Deltadm*S + coeffs[7]*S*Delta2 + coeffs[8]*S2*S \
        + coeffs[9]*S*dm2 + coeffs[10]*Deltadm*S2 + coeffs[11]*Delta2*Deltadm \
        + coeffs[12]*Delta2*Delta2 + coeffs[13]*S2*S2 + coeffs[14]*Delta2*S2 + coeffs[15]*dm2*dm2 + coeffs[16]*Deltadm*dm2 \
        + coeffs[17]*Delta2*dm2 + coeffs[18]*S2*dm2)

def _final_spin_diff_Healyetal(a_f, eta, delta_m, S, Delta, version):
    """ Internal function: the final spin with the Healy et al. fits is determined by minimizing this function """

    # calculate ISCO radius
    r_isco = calc_isco_radius(a_f)

    # angular momentum at ISCO -- Eq.(2.8) of Ori, Thorne Phys Rev D 62 124022 (2000)
    J_isco = (3*np.sqrt(r_isco)-2*a_f)*2./np.sqrt(3*r_isco)

    # fitting coefficients
    if version == "2014": # From Table XI of Healy et al Phys Rev D 90, 104004 (2014) [fourth order fits]
        L0  = 0.686710
        L1  = 0.613247
        L2a = -0.145427
        L2b = -0.115689
        L2c = -0.005254
        L2d = 0.801838
        L3a = -0.073839
        L3b = 0.004759
        L3c = -0.078377
        L3d = 1.585809
        L4a = -0.003050
        L4b = -0.002968
        L4c = 0.004364
        L4d = -0.047204
        L4e = -0.053099
        L4f = 0.953458
        L4g = -0.067998
        L4h = 0.001629
        L4i = -0.066693
    elif version == "2016": # From Table III of Healy and Lousto arXiv:1610.09713 (values taken from Matlab implementation and thus slightly more precise than the ones in the table)
        L0  = 0.686732132
        L1  = 0.613284976
        L2a = -0.148530075
        L2b = -0.113826318
        L2c = -0.00323995784
        L2d = 0.798011319
        L3a = -0.0687823713
        L3b = 0.00129103641
        L3c = -0.0780143929
        L3d = 1.55728564
        L4a = -0.00571010557
        L4b = 0.005919799
        L4c = -0.00170575554
        L4d = -0.0588819084
        L4e = -0.0101866693
        L4f = 0.964444768
        L4g = -0.11088507
        L4h = -0.00682082169
        L4i = -0.0816482139
    else:
        raise ValueError('Unknown version--should be either "2014" or "2016".')

    dm2 = delta_m*delta_m
    dm4 = dm2*dm2
    dm6 = dm4*dm2

    a_f_new = _RIT_symm_express(eta, delta_m, S, Delta, [L0, L1, L2a, L2b, L2c, L2d, L3a, L3b, L3c, L3d, L4a, L4b, L4c, L4d, L4e, L4f, L4g, L4h, L4i]) \
        + S*(1. + 8.*eta)*dm4 + eta*J_isco*dm6

    return abs(a_f-a_f_new)

def bbh_final_spin_non_precessing_Healyetal(m1, m2, chi1, chi2, version="2014"):
    """
    Calculate the spin of the final BH resulting from the merger of two black holes with non-precessing spins using fit from Healy et al Phys Rev D 90, 104004 (2014) (version == "2014") or the small update from Healy and Lousto arXiv:1610.09713 (version == "2016")

    Parameters
    ----------
    m1, m2 : component masses
    chi1, chi2 : dimensionless spins of two BHs

    Returns
    -------
    dimensionless final spin, chif
    """
    m1 = np.vectorize(float)(np.array(m1))
    m2 = np.vectorize(float)(np.array(m2))
    chi1 = np.vectorize(float)(np.array(chi1))
    chi2 = np.vectorize(float)(np.array(chi2))

    # Vectorize the function if arrays are provided as input
    if np.size(m1) * np.size(m2) * np.size(chi1) * np.size(chi2) > 1:
        return np.vectorize(bbh_final_spin_non_precessing_Healyetal)(m1, m2, chi1, chi2, version)

    eta, delta_m, S, Delta = _RIT_setup(m1, m2, chi1, chi2)

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # compute the final spin
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    x, cov_x = so.leastsq(_final_spin_diff_Healyetal, 0., args=(eta, delta_m, S, Delta, version))

    # The first element returned by so.leastsq() is a scalar in early versions of scipy (like 0.7.2) while it is a tuple of length 1 in later versions of scipy (like 0.10.1). The following bit ensures that a scalar is returned for a set of scalar inputs in a version-independent way.
    if hasattr(x, '__len__'):
      chif = x[0]
    else:
      chif = x

    return chif

def bbh_final_mass_non_precessing_Healyetal(m1, m2, chi1, chi2, version="2014", chif=None):
    """
    Calculate the mass of the final BH resulting from the merger of two black holes with non-precessing spins using fit from Healy et al Phys Rev D 90, 104004 (2014) (version == "2014") or the small update from Healy and Lousto arXiv:1610.09713 (version == "2016")

    Parameters
    ----------
    m1, m2 : component masses
    chi1, chi2 : dimensionless spins of two BHs
    chif: final spin (optional), if already calculated

    Returns
    -------
    final mass, mf
    """
    m1 = np.vectorize(float)(np.array(m1))
    m2 = np.vectorize(float)(np.array(m2))
    chi1 = np.vectorize(float)(np.array(chi1))
    chi2 = np.vectorize(float)(np.array(chi2))

    eta, delta_m, S, Delta = _RIT_setup(m1, m2, chi1, chi2)

    if chif is None:
        chif = bbh_final_spin_non_precessing_Healyetal(m1, m2, chi1, chi2, version)
    else:
        chif = np.array(chif)

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # now compute the final mass
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    r_isco = calc_isco_radius(chif)

    # fitting coefficients
    if version == "2014": # From Table XI of Healy et al Phys Rev D 90, 104004 (2014) [fourth order fits]
        M0  = 0.951507
        K1  = -0.051379
        K2a = -0.004804
        K2b = -0.054522
        K2c = -0.000022
        K2d = 1.995246
        K3a = 0.007064
        K3b = -0.017599
        K3c = -0.119175
        K3d = 0.025000
        K4a = -0.068981
        K4b = -0.011383
        K4c = -0.002284
        K4d = -0.165658
        K4e = 0.019403
        K4f = 2.980990
        K4g = 0.020250
        K4h = -0.004091
        K4i = 0.078441
    elif version == "2016": # From Table III of Healy and Lousto arXiv:1610.09713 (values taken from Matlab implementation and thus slightly more precise than the ones in the table)
        M0  = 0.951659087
        K1  = -0.0511301363
        K2a = -0.00569897591
        K2b = -0.0580644933
        K2c = -0.00186732281
        K2d = 1.99570464
        K3a = 0.00499137602
        K3b = -0.00923776244
        K3c = -0.120577082
        K3d = 0.0164168385
        K4a = -0.0607207285
        K4b = -0.00179828653
        K4c = 0.000654388173
        K4d = -0.156625642
        K4e = 0.0103033606
        K4f = 2.97872857
        K4g = 0.00790433045
        K4h = 0.000631241195
        K4i = 0.0844776942
    else:
        raise ValueError('Unknown version--should be either "2014" or "2016".')

    # binding energy at ISCO -- Eq.(2.7) of Ori, Thorne Phys Rev D 62 124022 (2000)
    E_isco = (1. - 2./r_isco + chif/r_isco**1.5)/np.sqrt(1. - 3./r_isco + 2.*chif/r_isco**1.5)

    dm2 = delta_m*delta_m
    dm6 = dm2*dm2*dm2

    m = m1 + m2

    # final mass -- Eq. (14) of Healy et al Phys Rev D 90, 104004 (2014)
    mf = _RIT_symm_express(eta, delta_m, S, Delta, [M0, K1, K2a, K2b, K2c, K2d, K3a, K3b, K3c, K3d, K4a, K4b, K4c, K4d, K4e, K4f, K4g, K4h, K4i]) \
        + (1+eta*(E_isco+11.))*dm6

    return mf*m

def bbh_final_spin_projected_spin_Healyetal(m1, m2, chi1, chi2, tilt1, tilt2):
    """
    wrapper to bbh_final_spin_projected_spins() for backwards compatibility
    """
    return bbh_final_spin_projected_spins(m1, m2, chi1, chi2, tilt1, tilt2, "HLZ2014")

def bbh_final_mass_projected_spin_Healyetal(m1, m2, chi1, chi2, tilt1, tilt2, chif=None):
    """
    wrapper to bbh_final_mass_projected_spins() for backwards compatibility
    """
    return bbh_final_mass_projected_spins(m1, m2, chi1, chi2, tilt1, tilt2, "HLZ2014", chif)

def bbh_final_spin_precessing_Healyetal_extension_Minit(m1, m2, chi1, chi2, tilt1, tilt2, phi12):
   """
   wrapper to bbh_final_spin_precessing() for backwards compatibility
   """
   return bbh_final_spin_precessing(m1, m2, chi1, chi2, tilt1, tilt2, phi12, "HLZ2014")

def bbh_final_mass_projected_spins(m1, m2, chi1, chi2, tilt1, tilt2, fitname, chif=None):
    """
    Calculate the mass of the final BH resulting from the merger of two black holes,
    only using the projected spins along the total angular momentum
    and some aligned-spin fit from the literature

    Parameters
    ----------
    m1, m2 : component masses
    chi1, chi2 : dimensionless spins of two BHs
    tilt1, tilt2 : tilts (in radians) in the new spin convention
    fitname: fit selection currently supports Pan2011 (non-spinning), HLZ2014, PhenomD, UIB2016, HL2016
    chif: final spin (optional, only used for HLZ2014 and HL2016), if already calculated

    Returns
    -------
    final mass, mf
    """

    m1    = np.vectorize(float)(np.array(m1))
    m2    = np.vectorize(float)(np.array(m2))
    chi1  = np.vectorize(float)(np.array(chi1))
    chi2  = np.vectorize(float)(np.array(chi2))
    tilt1 = np.vectorize(float)(np.array(tilt1))
    tilt2 = np.vectorize(float)(np.array(tilt2))
    if chif is not None:
       chif = np.vectorize(float)(np.array(chif))

    _check_mchi(m1,m2,chi1,chi2) # Check that inputs are physical

    chi1proj = chi1*np.cos(tilt1)
    chi2proj = chi2*np.cos(tilt2)

    if fitname==bbh_final_state_fits.pan2011:
       if np.any(chi1!=0) or np.any(chi2!=0) or np.any(tilt1!=0) or np.any(tilt2!=0):
          print("Note: Pan2011 fit does not use spins.")
       if chif is not None:
          print("Note: Precomputed chif not used by this fit.")
       mf = bbh_final_mass_non_spinning_Panetal(m1, m2)
    elif fitname==bbh_final_state_fits.hlz2014:
       mf = bbh_final_mass_non_precessing_Healyetal(m1, m2, chi1proj, chi2proj, version="2014", chif=chif)
    elif fitname==bbh_final_state_fits.hl2016:
       mf = bbh_final_mass_non_precessing_Healyetal(m1, m2, chi1proj, chi2proj, version="2016", chif=chif)
    elif fitname==bbh_final_state_fits.phenD:
       if chif is not None:
          print("Note: Precomputed chif not used by this fit.")
       mf = bbh_final_mass_non_precessing_Husaetal(m1, m2, chi1proj, chi2proj)
    elif fitname==bbh_final_state_fits.uib2016:
       if chif is not None:
          print("Note: Precomputed chif not used by this fit.")
       mf = bbh_final_mass_non_precessing_UIB2016(m1, m2, chi1proj, chi2proj, version="v2")
    elif fitname==bbh_final_state_fits.uib2016v1:
       if chif is not None:
          print("Note: Precomputed chif not used by this fit.")
       mf = bbh_final_mass_non_precessing_UIB2016(m1, m2, chi1proj, chi2proj, version="v1")
    else:
       raise ValueError("Unrecognized fit name.")

    return mf

def bbh_final_spin_projected_spins(m1, m2, chi1, chi2, tilt1, tilt2, fitname, truncate=bbh_Kerr_trunc_opts.trunc):
    """
    Calculate the (signed) dimensionless spin parameter of the final BH resulting from the merger of two black holes,
    only using the projected spins along the total angular momentum
    and some aligned-spin fit from the literature

    Parameters
    ----------
    m1, m2 : component masses
    chi1, chi2 : dimensionless spins of two BHs
    tilt1, tilt2 : tilts (in radians) in the new spin convention
    fitname: fit selection currently supports Pan2011 (non-spinning), HLZ2014, PhenomD, UIB2016, HBR2016, HL2016

    Returns
    -------
    (signed) dimensionless final spin parameter, chif
    """

    m1    = np.vectorize(float)(np.array(m1))
    m2    = np.vectorize(float)(np.array(m2))
    chi1  = np.vectorize(float)(np.array(chi1))
    chi2  = np.vectorize(float)(np.array(chi2))
    tilt1 = np.vectorize(float)(np.array(tilt1))
    tilt2 = np.vectorize(float)(np.array(tilt2))

    _check_mchi(m1,m2,chi1,chi2) # Check that inputs are physical

    chi1proj = chi1*np.cos(tilt1)
    chi2proj = chi2*np.cos(tilt2)

    if fitname==bbh_final_state_fits.pan2011:
       if np.any(chi1!=0) or np.any(chi2!=0) or np.any(tilt1!=0) or np.any(tilt2!=0):
          print("Note: Pan2011 fit does not use spins.")
       chif = bbh_final_spin_non_spinning_Panetal(m1, m2)
    elif fitname==bbh_final_state_fits.hlz2014:
       chif = bbh_final_spin_non_precessing_Healyetal(m1, m2, chi1proj, chi2proj, version="2014")
    elif fitname==bbh_final_state_fits.hl2016:
       chif = bbh_final_spin_non_precessing_Healyetal(m1, m2, chi1proj, chi2proj, version="2016")
    elif fitname==bbh_final_state_fits.phenD:
       chif = bbh_final_spin_non_precessing_Husaetal(m1, m2, chi1proj, chi2proj)
    elif fitname==bbh_final_state_fits.uib2016:
       chif = bbh_final_spin_non_precessing_UIB2016(m1, m2, chi1proj, chi2proj, version="v2")
    elif fitname==bbh_final_state_fits.uib2016v1:
       chif = bbh_final_spin_non_precessing_UIB2016(m1, m2, chi1proj, chi2proj, version="v1")
    elif fitname==bbh_final_state_fits.hbr2016:
       chif = bbh_final_spin_non_precessing_HBR2016(m1, m2, chi1proj, chi2proj)
    else:
       raise ValueError("Unrecognized fit name.")

    chif = _truncate_at_Kerr_limit(chif, truncate, fitname) # optionally truncate and/or warn at Kerr limit, depending on user option

    return chif

def bbh_final_spin_precessing(m1, m2, chi1, chi2, tilt1, tilt2, phi12, fitname, truncate=bbh_Kerr_trunc_opts.trunc):
    """
    Calculate the magnitude of the dimensionless spin parameter of the final BH resulting from the merger of two black holes,
    including the in-plane spin components;
    by either using a precessing fit from the literature;
    or by first projecting the spins along the angular momentum and using an aligned-spin fit from the literature,
    then adding in quadrature the in-plane dimensionful spin scaled by the initial mass squared.

    Parameters
    ----------
    m1, m2 : component masses
    chi1, chi2 : dimensionless spins of two BHs
    tilt1, tilt2 : tilts (in radians) in the new spin convention
    phi12: angle (in radians) between in-plane spin components
    fitname: fit selection currently supports Pan2011 (non-spinning), HLZ2014 (aligned+augmentation),
                                              PhenomD (aligned+augmentation), UIB2016 (aligned+augmentation),
                                              HL2016 (aligned+augmentation), HBR2016 (precessing)

    Returns
    -------
    magnitude of the dimensionless final spin parameter, chif
    """

    m1    = np.vectorize(float)(np.array(m1))
    m2    = np.vectorize(float)(np.array(m2))
    chi1  = np.vectorize(float)(np.array(chi1))
    chi2  = np.vectorize(float)(np.array(chi2))
    tilt1 = np.vectorize(float)(np.array(tilt1))
    tilt2 = np.vectorize(float)(np.array(tilt2))
    phi12 = np.vectorize(float)(np.array(phi12))

    _check_mchi(m1,m2,chi1,chi2) # Check that inputs are physical

    if fitname==bbh_final_state_fits.pan2011 or fitname==bbh_final_state_fits.hlz2014 or fitname==bbh_final_state_fits.phenD or fitname==bbh_final_state_fits.uib2016 or fitname==bbh_final_state_fits.uib2016v1 or fitname==bbh_final_state_fits.hl2016:
       precfit = False
    elif fitname==bbh_final_state_fits.hbr2016:
       precfit = True
       chif = bbh_final_spin_precessing_HBR2016(m1, m2, chi1, chi2, tilt1, tilt2, phi12)
    else:
       raise ValueError("Unrecognized fit name.")

    if not precfit:
       # First compute the final mass and parallel component of the final spin using the aligned components of the initial spins
       chifaligned = bbh_final_spin_projected_spins(m1, m2, chi1, chi2, tilt1, tilt2, fitname, truncate)

       # Now compute the squared magnitude of the in-plane dimensionful spin, first computing the magnitudes of the initial in-plane spins
       S1perpmag = m1*m1*chi1*np.sin(tilt1)
       S2perpmag = m2*m2*chi2*np.sin(tilt2)
       Sperpmag2 = S1perpmag*S1perpmag + S2perpmag*S2perpmag + 2.*S1perpmag*S2perpmag*np.cos(phi12)

       # Combine together
       chif = (chifaligned*chifaligned + Sperpmag2/(m1+m2)**4.)**0.5

    chif = _truncate_at_Kerr_limit(chif, truncate, "augmented "+fitname) # optionally truncate and/or warn at Kerr limit, depending on user option

    return chif

def bbh_final_mass_non_precessing_Husaetal(m1, m2, chi1, chi2):
    """
    Calculate the mass of the final BH resulting from the
    merger of two black holes with non-precessing spins using the fits
    used by IMRPhenomD, given in Eqs. (3.6) and (3.8) of Husa et al.
    arXiv:1508.07250. Note that Eq. (3.8) gives the radiated energy, not
    the final mass directly

    m1, m2: component masses
    chi1, chi2: dimensionless spins of two BHs
    """
    m1 = np.vectorize(float)(np.array(m1))
    m2 = np.vectorize(float)(np.array(m2))
    chi1 = np.vectorize(float)(np.array(chi1))
    chi2 = np.vectorize(float)(np.array(chi2))

    if np.any(abs(chi1)>1):
      raise ValueError("chi1 has to be in [-1, 1]")
    if np.any(abs(chi2)>1):
      raise ValueError("chi2 has to be in [-1, 1]")

    # binary parameters
    m = m1+m2
    msq = m*m

    eta = m1*m2/msq
    eta2 = eta*eta
    eta3 = eta2*eta
    eta4 = eta3*eta

    S1 = chi1*m1**2/msq                   # spin angular momentum 1 (in m = 1 units)
    S2 = chi2*m2**2/msq                   # spin angular momentum 2 (in m = 1 units)
    S = S1+S2                        # total spin
    Sh = S/(1. - 2.*eta)                    # rescaled total spin

    # Expressions copied from LALSimIMRPhenomD_internals.c (except with two notation differences: S is capitalized in chif and s -> Sh in Mf, in addition to the "m*(1. - ...)" to obtain the final mass from the radiated mass in m = 1 units which is calculated in the LAL code)

    Mf = m*(1. - ((0.055974469826360077*eta + 0.5809510763115132*eta2 - 0.9606726679372312*eta3 + 3.352411249771192*eta4)*
    (1. + (-0.0030302335878845507 - 2.0066110851351073*eta + 7.7050567802399215*eta2)*Sh))/(1. + (-0.6714403054720589 - 1.4756929437702908*eta + 7.304676214885011*eta2)*Sh))

    return Mf

def bbh_final_spin_non_precessing_Husaetal(m1, m2, chi1, chi2):
    """
    Calculate the spin of the final BH resulting from the
    merger of two black holes with non-precessing spins using the fits
    used by IMRPhenomD, given in Eqs. (3.6) and (3.8) of Husa et al.
    arXiv:1508.07250. Note that Eq. (3.8) gives the radiated energy, not
    the final mass directly

    m1, m2: component masses
    chi1, chi2: dimensionless spins of two BHs
    """
    # Vectorize the function if arrays are provided as input
    m1 = np.vectorize(float)(np.array(m1))
    m2 = np.vectorize(float)(np.array(m2))
    chi1 = np.vectorize(float)(np.array(chi1))
    chi2 = np.vectorize(float)(np.array(chi2))

    if np.any(abs(chi1)>1):
      raise ValueError("chi1 has to be in [-1, 1]")
    if np.any(abs(chi2)>1):
      raise ValueError("chi2 has to be in [-1, 1]")

    # binary parameters
    m = m1+m2
    msq = m*m

    eta = m1*m2/msq
    eta2 = eta*eta
    eta3 = eta2*eta
    eta4 = eta3*eta

    S1 = chi1*m1**2/msq                   # spin angular momentum 1 (in m = 1 units)
    S2 = chi2*m2**2/msq                   # spin angular momentum 2 (in m = 1 units)
    S = S1+S2                        # total spin
    Sh = S/(1. - 2.*eta)                    # rescaled total spin

    Ssq = S*S
    Scu = Ssq*S
    Squ = Scu*S

    # Expressions copied from LALSimIMRPhenomD_internals.c (except with two notation differences: S is capitalized in chif and s -> Sh in Mf, in addition to the "m*(1. - ...)" to obtain the final mass from the radiated mass in m = 1 units which is calculated in the LAL code)

    chif = 3.4641016151377544*eta - 4.399247300629289*eta2 + 9.397292189321194*eta3 - 13.180949901606242*eta4 + (1 - 0.0850917821418767*eta - 5.837029316602263*eta2)*S + (0.1014665242971878*eta - 2.0967746996832157*eta2)*Ssq + (-1.3546806617824356*eta + 4.108962025369336*eta2)*Scu + (-0.8676969352555539*eta + 2.064046835273906*eta2)*Squ

    return chif

def bbh_UIBfits_setup(m1, m2, chi1, chi2):
    """
    common setup function for UIB final-state and luminosity fit functions
    """

    # Vectorize the function if arrays are provided as input
    m1   = np.vectorize(float)(np.array(m1))
    m2   = np.vectorize(float)(np.array(m2))
    chi1 = np.vectorize(float)(np.array(chi1))
    chi2 = np.vectorize(float)(np.array(chi2))

    if np.any(m1<0):
      raise ValueError("m1 must not be negative")
    if np.any(m2<0):
      raise ValueError("m2 must not be negative")

    if np.any(abs(chi1)>1):
      raise ValueError("chi1 has to be in [-1, 1]")
    if np.any(abs(chi2)>1):
      raise ValueError("chi2 has to be in [-1, 1]")

    # binary masses
    m    = m1+m2
    if np.any(m<=0):
      raise ValueError("m1+m2 must be positive")
    msq  = m*m
    m1sq = m1*m1
    m2sq = m2*m2

    # symmetric mass ratio
    eta  = m1*m2/msq
    if np.any(eta>0.25):
      print("Truncating eta from above to 0.25. This should only be necessary in some rounding corner cases, but better check your m1 and m2 inputs...")
      eta = np.minimum(eta,0.25)
    if np.any(eta<0.0):
      print("Truncating negative eta to 0.0. This should only be necessary in some rounding corner cases, but better check your m1 and m2 inputs...")
      eta = np.maximum(eta,0.0)
    eta2 = eta*eta
    eta3 = eta2*eta
    eta4 = eta2*eta2

    # spin variables (in m = 1 units)
    S1    = chi1*m1sq/msq # spin angular momentum 1
    S2    = chi2*m2sq/msq # spin angular momentum 2
    Stot  = S1+S2         # total spin
    Shat  = (chi1*m1sq+chi2*m2sq)/(m1sq+m2sq) # effective spin, = msq*Stot/(m1sq+m2sq)
    Shat2 = Shat*Shat
    Shat3 = Shat2*Shat
    Shat4 = Shat2*Shat2

    # asymmetric spin combination (spin difference), where the paper assumes m1>m2
    # to make our implementation symmetric under simultaneous exchange of m1<->m2 and chi1<->chi2,
    # we flip the sign here when m2>m1
    chidiff  = chi1 - chi2
    if np.any(m2>m1):
      chidiff = np.sign(m1-m2)*chidiff
    chidiff2 = chidiff*chidiff

    # typical squareroots and functions of eta
    sqrt2 = 2.**0.5
    sqrt3 = 3.**0.5
    sqrt1m4eta = (1. - 4.*eta)**0.5

    return m, eta, eta2, eta3, eta4, Stot, Shat, Shat2, Shat3, Shat4, chidiff, chidiff2, sqrt2, sqrt3, sqrt1m4eta

def bbh_final_mass_non_precessing_UIB2016(m1, m2, chi1, chi2, version="v2"):
    """
    Calculate the final mass with the aligned-spin NR fit
    by Xisco Jimenez Forteza, David Keitel, Sascha Husa et al.
    [LIGO-P1600270] [https://arxiv.org/abs/1611.00332]
    versions v1 and v2 use the same ansatz,
    with v2 calibrated to additional SXS and RIT data

    m1, m2: component masses
    chi1, chi2: dimensionless spins of two BHs
    Results are symmetric under simultaneous exchange of m1<->m2 and chi1<->chi2.
    """

    m, eta, eta2, eta3, eta4, Stot, Shat, Shat2, Shat3, Shat4, chidiff, chidiff2, sqrt2, sqrt3, sqrt1m4eta = bbh_UIBfits_setup(m1, m2, chi1, chi2)

    if version == "v1":
        # rational-function Pade coefficients (exact) from Eq. (22) of 1611.00332v1
        b10 = 0.487
        b20 = 0.295
        b30 = 0.17
        b50 = -0.0717
        # fit coefficients from Tables VII-X of 1611.00332v1
        # values at increased numerical precision copied from
        # https://gravity.astro.cf.ac.uk/cgit/cardiff_uib_share/tree/Luminosity_and_Radiated_Energy/UIBfits/LALInference/EradUIB2016_pyform_coeffs.txt
        # git commit 636e5a71462ecc448060926890aa7811948d5a53
        a2 = 0.5635376058169299
        a3 = -0.8661680065959881
        a4 = 3.181941595301782
        b1 = -0.15800074104558132
        b2 = -0.15815904609933157
        b3 = -0.14299315232521553
        b5 = 8.908772171776285
        f20 = 3.8071100104582234
        f30 = 25.99956516423936
        f50 = 1.552929335555098
        f10 = 1.7004558922558886
        f21 = 0.
        d10 = -0.12282040108157262
        d11 = -3.499874245551208
        d20 = 0.014200035799803777
        d30 = -0.01873720734635449
        d31 = -5.1830734185518725
        f11 = 14.39323998088354
        f31 = -232.25752840151296
        f51 = -0.8427987782523847

    elif version == "v2":
        # rational-function Pade coefficients (exact) from Eq. (22) of LIGO-P1600270-v4
        b10 = 0.346
        b20 = 0.211
        b30 = 0.128
        b50 = -0.212
        # fit coefficients from Tables VII-X of LIGO-P1600270-v4
        # values at increased numerical precision copied from
        # https://dcc.ligo.org/DocDB/0128/P1600270/004/FinalStateUIB2016_suppl_Erad_coeffs.txt
        a2 = 0.5609904135313374
        a3 = -0.84667563764404
        a4 = 3.145145224278187
        b1 = -0.2091189048177395
        b2 = -0.19709136361080587
        b3 = -0.1588185739358418
        b5 = 2.9852925538232014
        f20 = 4.271313308472851
        f30 = 31.08987570280556
        f50 = 1.5673498395263061
        f10 = 1.8083565298668276
        f21 = 0.
        d10 = -0.09803730445895877
        d11 = -3.2283713377939134
        d20 = 0.01118530335431078
        d30 = -0.01978238971523653
        d31 = -4.91667749015812
        f11 = 15.738082204419655
        f31 = -243.6299258830685
        f51 = -0.5808669012986468

    else:
        raise ValueError('Unknown version -- should be either "v1" or "v2".')

    # Calculate the radiated-energy fit from Eq. (28) of LIGO-P1600270-v4
    Erad = (((1. + -2.0/3.0*sqrt2)*eta + a2*eta2 + a3*eta3 + a4*eta4)*(1. + b10*b1*Shat*(f10 + f11*eta + (16. - 16.*f10 - 4.*f11)*eta2) + b20*b2*Shat2*(f20 + f21*eta + (16. - 16.*f20 - 4.*f21)*eta2) + b30*b3*Shat3*(f30 + f31*eta + (16. - 16.*f30 - 4.*f31)*eta2)))/(1. + b50*b5*Shat*(f50 + f51*eta + (16. - 16.*f50 - 4.*f51)*eta2)) + d10*sqrt1m4eta*eta2*(1. + d11*eta)*chidiff + d30*Shat*sqrt1m4eta*eta*(1. + d31*eta)*chidiff + d20*eta3*chidiff2

    # Convert to actual final mass
    Mf = m*(1.-Erad)

    return Mf

def bbh_final_spin_non_precessing_UIB2016(m1, m2, chi1, chi2, version="v2"):
    """
    Calculate the final spin with the aligned-spin NR fit
    by Xisco Jimenez Forteza, David Keitel, Sascha Husa et al.
    [LIGO-P1600270] [https://arxiv.org/abs/1611.00332]
    versions v1 and v2 use the same ansatz,
    with v2 calibrated to additional SXS and RIT data

    m1, m2: component masses
    chi1, chi2: dimensionless spins of two BHs
    Results are symmetric under simultaneous exchange of m1<->m2 and chi1<->chi2.
    """

    m, eta, eta2, eta3, eta4, Stot, Shat, Shat2, Shat3, Shat4, chidiff, chidiff2, sqrt2, sqrt3, sqrt1m4eta = bbh_UIBfits_setup(m1, m2, chi1, chi2)

    if version == "v1":
        # rational-function Pade coefficients (exact) from Eqs. (7) and (8) of 1611.00332v1
        a20 = 5.28
        a30 = 1.27
        a50 = 2.89
        b10 = -0.194
        b20 = 0.075
        b30 = 0.00782
        b50 = -0.527
        # fit coefficients from Tables I-IV of 1611.00332v1
        # values at increased numerical precision copied from
        # https://gravity.astro.cf.ac.uk/cgit/cardiff_uib_share/tree/Luminosity_and_Radiated_Energy/UIBfits/LALInference/FinalSpinUIB2016_pyform_coeffs.txt
        # git commit 636e5a71462ecc448060926890aa7811948d5a53
        a2 = 3.772362507208651
        a3 = -9.627812453422376
        a5 = 2.487406038123681
        b1 = 1.0005294518146604
        b2 = 0.8823439288807416
        b3 = 0.7612809461506448
        b5 = 0.9139185906568779
        f21 = 8.887933111404559
        f31 = 23.927104476660883
        f50 = 1.8981657997557002
        f11 = 4.411041530972546
        f52 = 0.
        d10 = 0.2762804043166152
        d11 = 11.56198469592321
        d20 = -0.05975750218477118
        d30 = 2.7296903488918436
        d31 = -3.388285154747212
        f12 = 0.3642180211450878
        f22 = -40.35359764942015
        f32 = -178.7813942566548
        f51 = -5.556957394513334

    elif version == "v2":
        # rational-function Pade coefficients (exact) from Eqs. (7) and (8) of LIGO-P1600270-v4
        a20 = 5.24
        a30 = 1.3
        a50 = 2.88
        b10 = -0.194
        b20 = 0.0851
        b30 = 0.00954
        b50 = -0.579
        # fit coefficients from Tables I-IV of LIGO-P1600270-v4
        # values at increased numerical precision copied from
        # https://dcc.ligo.org/DocDB/0128/P1600270/004/FinalStateUIB2016_suppl_spin_coeffs.txt
        a2 = 3.8326341618708577
        a3 = -9.487364155598392
        a5 = 2.5134875145648374
        b1 = 1.0009563702914628
        b2 = 0.7877509372255369
        b3 = 0.6540138407185817
        b5 = 0.8396665722805308
        f21 = 8.77367320110712
        f31 = 22.830033250479833
        f50 = 1.8804718791591157
        f11 = 4.409160174224525
        f52 = 0.
        d10 = 0.3223660562764661
        d11 = 9.332575956437443
        d20 = -0.059808322561702126
        d30 = 2.3170397514509933
        d31 = -3.2624649875884852
        f12 = 0.5118334706832706
        f22 = -32.060648277652994
        f32 = -153.83722669033995
        f51 = -4.770246856212403

    else:
        raise ValueError('Unknown version -- should be either "v1" or "v2".')

    # Calculate the fit for the Lorb' quantity from Eq. (16) of LIGO-P1600270-v4
    Lorb = (2.*sqrt3*eta + a20*a2*eta2 + a30*a3*eta3)/(1. + a50*a5*eta) + (b10*b1*Shat*(f11*eta + f12*eta2 + (64. - 16.*f11 - 4.*f12)*eta3) + b20*b2*Shat2*(f21*eta + f22*eta2 + (64. - 16.*f21 - 4.*f22)*eta3) + b30*b3*Shat3*(f31*eta + f32*eta2 + (64. - 16.*f31 - 4.*f32)*eta3))/(1. + b50*b5*Shat*(f50 + f51*eta + f52*eta2 + (64. - 64.*f50 - 16.*f51 - 4.*f52)*eta3)) + d10*sqrt1m4eta*eta2*(1. + d11*eta)*chidiff + d30*Shat*sqrt1m4eta*eta3*(1. + d31*eta)*chidiff + d20*eta3*chidiff2

    # Convert to actual final spin
    chif = Lorb + Stot

    return chif

def _bbh_HBR2016_setup(m1, m2, chi1, chi2):
    """
    Setup function for the Hofmann, Barausse, and Rezzolla final spin fits to vectorize the masses and spins and calculate the mass ratio.
    """

    # Vectorize if arrays are provided as input
    m1 = np.vectorize(float)(np.array(m1))
    m2 = np.vectorize(float)(np.array(m2))
    chi1 = np.vectorize(float)(np.array(chi1))
    chi2 = np.vectorize(float)(np.array(chi2))

    return m1, m2, chi1, chi2, m2/m1

def _bbh_HBR2016_ell(m1, m2, chi1z, chi2z, version):
    """
    Compute the orbital angular momentum ell used by the Hofmann, Barausse, and Rezzolla final spin fit [from ApJL 825, L19 (2016)], henceforth HBR.
    The three versions available correspond to the three choices of fit coefficients given in Table 1 of that paper, with 6, 16, and 20 coefficients, respectively.

    version can thus be "M1J2", "M3J3", or "M3J4"
    """

    # Set upper bounds for sums and coefficients from Table 1 in HBR; k00 is calculated from Eq. (11) in the paper

    if version == "M1J2":
        nM = 1
        nJ = 2

        k00 = -3.82116
        k01 = -1.2019
        k02 = -1.20764
        k10 = 3.79245
        k11 = 1.18385
        k12 = 4.90494
        xi = 0.41616

    elif version == "M3J3":
        nM = 3
        nJ = 3

        k00 = -5.83164
        k01 = 2.87025
        k02 = -1.53315
        k03 = -3.78893
        k10 = 32.9127
        k11 = -62.9901
        k12 = 10.0068
        k13 = 56.1926
        k20 = -136.832
        k21 = 329.32
        k22 = -13.2034
        k23 = -252.27
        k30 = 210.075
        k31 = -545.35
        k32 = -3.97509
        k33 = 368.405
        xi = 0.463926

    elif version == "M3J4":
        nM = 3
        nJ = 4

        k00 = -5.97723
        k01 = 3.39221
        k02 = 4.48865
        k03 = -5.77101
        k04 = -13.0459
        k10 = 35.1278
        k11 = -72.9336
        k12 = -86.0036
        k13 = 93.7371
        k14 = 200.975
        k20 = -146.822
        k21 = 387.184
        k22 = 447.009
        k23 = -467.383
        k24 = -884.339
        k30 = 223.911
        k31 = -648.502
        k32 = -697.177
        k33 = 753.738
        k34 = 1166.89
        xi = 0.474046

    else:
        raise ValueError('Unknown version--should be either "M1J2", "M3J3", or "M3J4".')

    # Calculate eta, atot, and aeff; note that HBR call the symmetric mass ratio nu instead of eta

    m = m1 + m2
    q = m2/m1
    eta = m1*m2/(m*m)

    atot = (chi1z + q*q*chi2z)/((1.+q)*(1.+q)) # Eq. (12) in HBR
    aeff = atot + xi*eta*(chi1z + chi2z) # Inline equation below Eq. (7); see also Eq. (15) for the precessing analogue

    # Calculate ISCO energy and angular momentum
    r_isco = calc_isco_radius(aeff)
    e_isco = (1. - 2./3./r_isco)**0.5 # Eq. (2) in HBR
    l_isco = 2./(3.*3.**0.5)*(1. + 2.*(3.*r_isco - 2.)**0.5) # Eq. (3) in HBR

    # The following expressions give the double sum in Eq. (13) in HBR (without the overall factor of eta that is put in below) specialized to the nJ values for the three versions, i.e., nJ = 2 => M1J2, nJ = 3 => M3J3, nJ = 4 => M3J4
    if nJ >= 2:
        aeff2 = aeff*aeff
        ksum = k00 + k01*aeff + k02*aeff2 + eta*(k10 + k11*aeff + k12*aeff2)
    if nJ >= 3:
        eta2 = eta*eta
        eta3 = eta2*eta
        aeff3 = aeff2*aeff
        ksum = ksum + (k03 + eta*k13)*aeff3 + eta2*(k20 + k21*aeff + k22*aeff2 + k23*aeff3) + eta3*(k30 + k31*aeff + k32*aeff2 + k33*aeff3)
    if nJ >= 4:
        ksum = ksum + (k04 + eta*k14 + eta2*k24 + eta3*k34)*aeff3*aeff

    # Calculate the absolute value of ell
    ell = abs((l_isco - 2.*atot*(e_isco - 1.)) + eta*ksum) # Eq. (13) in HBR

    return ell

def bbh_final_spin_non_precessing_HBR2016(m1, m2, chi1z, chi2z, version="M3J3"):
        """
        Calculate the (signed) dimensionless spin of the final BH resulting from the
        merger of two black holes with aligned spins using the fit from Hofmann, Barausse, and Rezzolla ApJL 825, L19 (2016), henceforth HBR.

        The three versions available correspond to the three choices of fit coefficients given in Table 1 of that paper, with 6, 16, and 20 coefficients, respectively.

        version can thus be "M1J2", "M3J3", or "M3J4"

        m1, m2: component masses
        chi1z, chi2z: components of the dimensionless spins of the two BHs along the orbital angular momentum
        """

        # Calculate q and vectorize the masses and spins if arrays are provided as input

        m1, m2, chi1z, chi2z, q = _bbh_HBR2016_setup(m1, m2, chi1z, chi2z)

        # Calculate the final spin

        atot = (chi1z + chi2z*q*q)/((1.+q)*(1.+q)) # Eq. (12) in HBR

        ell = _bbh_HBR2016_ell(m1, m2, chi1z, chi2z, version)

        return atot + ell/(1./q + 2. + q) # Eq. (12) in HBR, writing the symmetric mass ratio in terms of q

def bbh_final_spin_precessing_HBR2016(m1, m2, chi1, chi2, tilt1, tilt2, phi12, version="M3J3"):
        """
        Calculate the dimensionless spin of the final BH resulting from the
        merger of two black holes with precessing spins using the fit from Hofmann, Barausse, and Rezzolla ApJL 825, L19 (2016), henceforth HBR.

        The three versions available correspond to the three choices of fit coefficients given in Table 1 of that paper, with 6, 16, and 20 coefficients, respectively.

        version can thus be "M1J2", "M3J3", or "M3J4"

        m1, m2: component masses
        chi1, chi2: dimensionless spins of two BHs
        tilt1, tilt2: tilt angles of the spins from the orbital angular momentum
        phi12: angle between in-plane spin components
        """

        # Vectorize the function if arrays are provided as input
        if np.size(m1) * np.size(m2) * np.size(chi1) * np.size(chi2) * np.size(tilt1) * np.size(tilt2) * np.size(phi12) > 1:
            return np.vectorize(bbh_final_spin_precessing_HBR2016)(m1, m2, chi1, chi2, tilt1, tilt2, phi12, version)

        # Calculate q and vectorize the masses and spins if arrays are provided as input

        m1, m2, chi1, chi2, q = _bbh_HBR2016_setup(m1, m2, chi1, chi2)

        # Vectorize the spin angles if arrays are provided as input
        tilt1 = np.vectorize(float)(np.array(tilt1))
        tilt2 = np.vectorize(float)(np.array(tilt2))
        phi12 = np.vectorize(float)(np.array(phi12))

        # Set eps (\epsilon_\beta or \epsilon_\gamma) to the value given below Eq. (18) in HBR

        eps = 0.024

        # Computing angles defined in Eq. (17) of HBR. The betas and gammas expressions are for the starred quantities computed using the second (approximate) equality in Eq. (18) in HBR
        cos_beta = np.cos(tilt1)
        cos_betas = np.cos(tilt1 + eps*np.sin(tilt1))
        cos_gamma = np.cos(tilt2)
        cos_gammas = np.cos(tilt2 + eps*np.sin(tilt2))
        cos_alpha = ((1 - cos_beta*cos_beta)*(1 - cos_gamma*cos_gamma))**0.5*np.cos(phi12) + cos_beta*cos_gamma # This rewrites the inner product definition of cos_alpha in terms of cos_beta, cos_gamma, and phi12

        # Define a shorthand and compute the final spin

        q2 = q*q

        ell = _bbh_HBR2016_ell(m1, m2, chi1*cos_betas, chi2*cos_gammas, version)

        # Compute the final spin value [Eq. (16) in HBR], truncating the argument of the square root at zero if it becomes negative
        sqrt_arg = chi1*chi1 + chi2*chi2*q2*q2 + 2.*chi1*chi2*q2*cos_alpha + 2.*(chi1*cos_betas + chi2*q2*cos_gammas)*ell*q + ell*ell*q2
        if sqrt_arg < 0.:
            print("bbh_final_spin_precessing_HBR2016(): The argument of the square root is %f; truncating it to zero."%sqrt_arg)
            sqrt_arg = 0.

        # Return the final spin value [Eq. (16) in HBR]
        return sqrt_arg**0.5/((1.+q)*(1.+q))

#######################
# Peak luminosity fits
#######################

# list of peak luminosity fit names
class bbh_Lpeak_fits:
    t1600018  = "T1600018" # Jimenez Forteza et al. [LIGO-T1600018 (2016)]         (aligned)
    uib2016   = "UIB2016"  # Keitel et al.          [arXiv:1612.09566]             (aligned)
    hl2016    = "HL2016"   # Healy and Lousto       [arXiv:1610.09713]             (aligned)

def _rescale_Lpeak(Lpeak):
    """
    convert from geometric units ("Planck luminosity" of c^5/G) to multiples of 10^56 erg/s = 10^49 J/s
    """
    LumPl_ergs_per_sec = lal.LUMPL_SI*1e-49 # Approximate value = 3628.505
    return LumPl_ergs_per_sec*Lpeak

def bbh_aligned_Lpeak_6mode_SHXJDK(q, chi1, chi2):
    """
    wrapper to bbh_peak_luminosity_non_precessing_T1600018() for backwards compatibility

    q: mass ratio (here m2/m1, where m1>m2)
    chi1: the component of the dimensionless spin of m1 along the angular momentum (z)
    chi2: the component of the dimensionless spin of m2 along the angular momentum (z)
    """
    # Vectorize the function if arrays are provided as input
    q = np.vectorize(float)(np.array(q))

    # from bayespputils.py convention is q = m2/m1, where m1>m2.
    if np.any(q<=0.):
      raise ValueError("q has to be > 0.")
    if np.any(q>1.):
      raise ValueError("q has to be <= 1.")

    # T1600018 fit expects m1>m2;
    # implementation here should be symmetric,
    # but we follow that convention to be sure
    m1 = 1./(1.+q)
    m2 = q/(1.+q)

    return bbh_peak_luminosity_non_precessing_T1600018(m1, m2, chi1, chi2)

def bbh_peak_luminosity_non_precessing_T1600018(m1, m2, chi1, chi2):
    """
    Calculate the peak luminosity (using modes 22, 21, 33, 32, 44, and 43) of a binary black hole with aligned spins using the fit made by Sascha Husa, Xisco Jimenez Forteza, David Keitel [LIGO-T1500598] using 5th order in chieff and return results in units of 10^56 ergs/s

    m1, m2: component masses
    chi1: the component of the dimensionless spin of m1 along the angular momentum (z)
    chi2: the component of the dimensionless spin of m2 along the angular momentum (z)
    Results are symmetric under simultaneous exchange of m1<->m2 and chi1<->chi2.
    """

    # Vectorize the function if arrays are provided as input
    m1 = np.vectorize(float)(np.array(m1))
    m2 = np.vectorize(float)(np.array(m2))
    chi1 = np.vectorize(float)(np.array(chi1))
    chi2 = np.vectorize(float)(np.array(chi2))

    # Calculate powers of eta and the effective spin S (not used in this fit)
    m, eta, eta2, eta3, eta4, Stot, Shat, Shat2, Shat3, Shat4, chidiff, chidiff2, sqrt2, sqrt3, sqrt1m4eta = bbh_UIBfits_setup(m1, m2, chi1, chi2)

    dm2 = 1. - 4.*eta # equivalent to sqrt1m4eta**2

    chi_eff = (m1*chi1+m2*chi2)/m # equivalent to (q_inv*chi1 + chi2)/(1. + q_inv)
    chi_eff2 = chi_eff*chi_eff
    chi_eff3 = chi_eff2*chi_eff
    chi_eff4 = chi_eff3*chi_eff
    chi_eff5 = chi_eff4*chi_eff

    # Calculate best fit (from [https://dcc.ligo.org/T1500598-v4])
    Lpeak = (0.012851338846828302 + 0.007822265919928252*chi_eff + 0.010221856361035788*chi_eff2 + 0.015805535732661396*chi_eff3 + 0.0011356206806770043*chi_eff4 - 0.009868152529667197*chi_eff5)*eta2 + (0.05681786589129071 - 0.0017473702709303457*chi_eff - 0.10150706091341818*chi_eff2 - 0.2349153289253309*chi_eff3 + 0.015657737820040145*chi_eff4 + 0.19556893194885075*chi_eff5)*eta4 + 0.026161288241420833*dm2**0.541825641769908*eta**3.1629576945611757*chidiff + 0.0007771032100485481*dm2**0.4499151697918658*eta**1.7800346166040835*chidiff2

    # Convert to 10^56 ergs/s units
    return _rescale_Lpeak(Lpeak)

def bbh_peak_luminosity_non_precessing_UIB2016(m1, m2, chi1, chi2):
    """
    Calculate the peak luminosity with the aligned-spin NR fit
    by David Keitel, Xisco Jimenez Forteza, Sascha Husa, Lionel London et al.
    [LIGO-P1600279-v5] [https://arxiv.org/abs/1612.09566v1]
    using modes up to lmax=6,
    and return results in units of 10^56 ergs/s

    m1, m2: component masses
    chi1, chi2: dimensionless spins of two BHs
    Results are symmetric under simultaneous exchange of m1<->m2 and chi1<->chi2.
    """

    m, eta, eta2, eta3, eta4, Stot, Shat, Shat2, Shat3, Shat4, chidiff, chidiff2, sqrt2, sqrt3, sqrt1m4eta = bbh_UIBfits_setup(m1, m2, chi1, chi2)
    eta5 = eta3*eta2

    # fit coefficients corresponding to Table I, II, IV,
    # exact values corresponding to https://arxiv.org/src/1612.09566v1/anc/LpeakUIB2016_suppl_coeffs.txt
    # fi2 coefficients are replaced by functions of fi0, fi1 as in Eq. (10)
    a0 = 0.8742169580717333
    a1 = -2.111792574893241
    a2 = 35.214103272783646
    a3 = -244.94930678226913
    a4 = 877.1061892200927
    a5 = -1172.549896493467
    b1 = 0.9800204548606681
    b2 = -0.1779843936224084
    b4 = 1.7859209418791981
    d10 = 3.789116271213293
    d20 = 0.40214125006660567
    d30 = 4.273116678713487
    f10 = 1.6281049269810424
    f11 = -3.6322940180721037
    f20 = 31.710537408279116
    f21 = -273.84758785648336
    f30 = -0.23470852321351202
    f31 = 6.961626779884965
    f40 = 0.21139341988062182
    f41 = 1.5255885529750841
    f60 = 3.0901740789623453
    f61 = -16.66465705511997
    f70 = 0.8362061463375388
    f71 = 0.

    # calculate the Lpeak/(eta2*L0) fit from Eq. (14), using the constraints for fi2 from Eq. (10)
    Lpeak = a0 + a1*eta + a2*eta2 + a3*eta3 + a4*eta4 + a5*eta5 + (0.465*b1*Shat*(f10 + f11*eta + (16. - 16.*f10 - 4.*f11)*eta2) + 0.107*b2*Shat2*(f20 + f21*eta + (16. - 16.*f20 - 4.*f21)*eta2) + Shat3*(f30 + f31*eta + (-16.*f30 - 4.*f31)*eta2) + Shat4*(f40 + f41*eta + (-16.*f40 - 4.*f41)*eta2))/(1. - 0.328*b4*Shat*(f60 + f61*eta + (16. - 16.*f60 - 4.*f61)*eta2) + Shat2*(f70 + f71*eta + (-16.*f70 - 4.*f71)*eta2)) + eta3*((d10+d30*Shat)*sqrt1m4eta*chidiff + d20*chidiff2)

    # Lpeak(eta=0.25,chi1=chi2=0)/0.25^2
    L0 = 0.016379197203103536

    # Convert to actual peak luminosity
    Lpeak = Lpeak*eta2*L0

    # Convert to 10^56 ergs/s units
    return _rescale_Lpeak(Lpeak)

def bbh_peak_luminosity_non_precessing_Healyetal(m1, m2, chi1z, chi2z):
    """
    Calculate the peak luminosity of an aligned-spin binary black hole coalescence using the fit from Healy and Lousto arXiv:1610.09713 (henceforth HL)

    Parameters
    ----------
    m1, m2 : component masses
    chi1z, chi2z : components of the dimensionless spins of the two BHs along their orbital angular momentum

    Returns
    -------
    peak luminosity in units of 10^56 ergs/s

    """

    m1 = np.vectorize(float)(np.array(m1))
    m2 = np.vectorize(float)(np.array(m2))
    chi1z = np.vectorize(float)(np.array(chi1z))
    chi2z = np.vectorize(float)(np.array(chi2z))

    eta, delta_m, S, Delta = _RIT_setup(m1, m2, chi1z, chi2z)

    # fitting coefficients (from the RIT group's private Matlab implementation; versions to fewer digits are given in Table IV of HL)
    coeffs = [0.00102101737,0.000897428902,-9.77467189e-05,0.000920883841,1.86982704e-05,-0.000391316975,-0.000120214144,0.000148102239,0.00137901473,-0.000493730555,0.000884792724,3.29254224e-07,1.70170729e-05,0.00151484008,-0.000148456828,0.000136659593,0.000160343115,-6.18530577e-05,-0.00103602521]

    # Return the peak luminosity in units of 10^56 ergs/s [cf. Eq. (3) in HL]
    Lpeak = _RIT_symm_express(eta, delta_m, S, Delta, coeffs)

    return _rescale_Lpeak(Lpeak)

def bbh_peak_luminosity_projected_spins(m1, m2, chi1, chi2, tilt1, tilt2, fitname):
    """
    Calculate the peak luminosity of the merger of two black holes,
    only using the projected spins along the total angular momentum
    and some aligned-spin fit from the literature

    Parameters
    ----------
    m1, m2 : component masses
    chi1, chi2 : dimensionless spins of two BHs
    tilt1, tilt2 : tilts (in radians) in the new spin convention
    fitname: fit selection currently supports T1600018, UIB2016, HL2016

    Returns
    -------
    peak luminosity, Lpeak, in units of 10^56 ergs/s
    """

    m1    = np.vectorize(float)(np.array(m1))
    m2    = np.vectorize(float)(np.array(m2))
    chi1  = np.vectorize(float)(np.array(chi1))
    chi2  = np.vectorize(float)(np.array(chi2))
    tilt1 = np.vectorize(float)(np.array(tilt1))
    tilt2 = np.vectorize(float)(np.array(tilt2))

    _check_mchi(m1,m2,chi1,chi2) # Check that inputs are physical

    chi1proj = chi1*np.cos(tilt1)
    chi2proj = chi2*np.cos(tilt2)

    if fitname==bbh_Lpeak_fits.t1600018:
       Lpeak = bbh_peak_luminosity_non_precessing_T1600018(m1, m2, chi1proj, chi2proj)
    elif fitname==bbh_Lpeak_fits.uib2016:
       Lpeak = bbh_peak_luminosity_non_precessing_UIB2016(m1, m2, chi1proj, chi2proj)
    elif fitname==bbh_Lpeak_fits.hl2016:
       Lpeak = bbh_peak_luminosity_non_precessing_Healyetal(m1, m2, chi1proj, chi2proj)
    else:
       raise ValueError("Unrecognized fit name.")

    return Lpeak

#########################
# Convenience functions
#########################

def bbh_average_fits_precessing(m1, m2, chi1, chi2, tilt1, tilt2, phi12, quantity, fits):
    """
    Calculate the average of the predictions of various fits, currently just for the final mass, final spin, or peak luminosity of a quasicircular
    binary black hole coalescence. The final spin calculation returns the magnitude of the final spin, including the contribution
    from in-plane spins; the other two cases just apply aligned-spin fits to the components of the spins along the orbital angular momentum.

    Parameters
    ----------
    m1, m2 : component masses
    chi1, chi2 : dimensionless spins of two BHs
    tilt1, tilt2 : tilts (in radians) in the new spin convention
    phi12: angle (in radians) between in-plane spin components (only used for the final spin)
    quantity: "Mf", "af", "afz", or "Lpeak"
    fits: An array of fit names to be used. The possible fit names are those known by bbh_final_mass_projected_spins, bbh_final_spin_precessing,
          and bbh_peak_luminosity_projected_spins

    The shape of m1 is used to determine the shape of the output.

    Returns
    -------
    Average of the results for the given fits for the chosen quantity
    """

    if quantity != "af" and max(abs(np.atleast_1d(phi12))) != 0:
        print("Note: phi12 is only used for the full final spin calculation.")

    if quantity not in ["Mf", "af", "afz", "Lpeak"]:
        raise ValueError("Unknown quantity: %s"%quantity)

    # Define function to return the appropriate quantity

    def _return_quantity(m1, m2, chi1, chi2, tilt1, tilt2, phi12, fitname):
        if quantity == "Mf":
            return bbh_final_mass_projected_spins(m1, m2, chi1, chi2, tilt1, tilt2, fitname)
        elif quantity == "af":
            return bbh_final_spin_precessing(m1, m2, chi1, chi2, tilt1, tilt2, phi12, fitname)
        elif quantity == "afz":
            return bbh_final_spin_projected_spins(m1, m2, chi1, chi2, tilt1, tilt2, fitname)
        elif quantity == "Lpeak":
            return bbh_peak_luminosity_projected_spins(m1, m2, chi1, chi2, tilt1, tilt2, fitname)

    # Define a function that returns the length of an array when passed an array and 1 when not passed an array
    def _len_smart(x):
        if hasattr(x, '__len__'):
            return len(x)
        else:
            return 1

    # Initialize

    fits = np.atleast_1d(fits)

    if not fits.size: # Check to make sure that fits is not an empty array
        raise ValueError("The list of fits passed cannot be an empty array.")

    num_fits = len(fits)
    num_data = max(_len_smart(m1), _len_smart(m2), _len_smart(chi1), _len_smart(chi2), _len_smart(tilt1), _len_smart(tilt2), _len_smart(phi12))

    data = np.zeros([num_fits, num_data])

    # Convert the input from bayespputils into a format that is compatible with the following loop

    if hasattr(m1, 'shape'):
        data_shape = m1.shape
    else:
        data_shape = -1

    # Loop over the fits

    for k, fit in enumerate(fits):
        data_portion = _return_quantity(m1, m2, chi1, chi2, tilt1, tilt2, phi12, fit)
        data[k] = data_portion.reshape(-1)

    # Calculate average:
    data_avg = np.mean(data,axis=0)

    return data_avg.reshape(data_shape)
