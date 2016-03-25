""" 
Various fitting formulas provided by numerical relativity 

P. Ajith, 2015-04-09; updated by NKJ-M, March 2016
"""

import numpy as np 
import scipy.optimize as so

def bbh_final_mass_non_spinning_Panetal(m1, m2):
    """
    Calculate the mass of the final BH resulting from the merger of two non-spinning black holes using fit from Pan et al, Phys Rev D 84, 124052 (2011).
    
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
    Calculate the spin of the final BH resulting from the merger of two non-spinning black holes using fit from Pan et al, Phys Rev D 84, 124052 (2011).
    
    Parameters
    ----------
    m1, m2 : component masses
    
    Returns
    ------
    final spin, sf
    """
    m1 = np.vectorize(float)(np.array(m1))
    m2 = np.vectorize(float)(np.array(m2))
    
    eta = m1*m2/(m1+m2)**2.
    return np.sqrt(12.)*eta - 3.871*(eta**2.) + 4.028*(eta**3)

def calc_isco_radius(a):
    """
    Calculate the ISCO radius of a Kerr BH as a function of the Kerr parameter
    
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

def _final_spin_diff_Healyetal(a_f, eta, delta_m, S, Delta):
    """ Internal function: the final spin is determined by minimizing this function """
    
    # calculate ISCO radius
    r_isco = calc_isco_radius(a_f)
    
    # angular momentum at ISCO -- Eq.(2.8) of Ori, Thorne Phys Rev D 62 124022 (2000)
    J_isco = (3*np.sqrt(r_isco)-2*a_f)*2./np.sqrt(3*r_isco)
    
    # fitting coefficients - Table XI of Healy et al Phys Rev D 90, 104004 (2014)
    # [fourth order fits]
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
    
    a_f_new = (4.*eta)**2.*(L0  +  L1*S +  L2a*Delta*delta_m + L2b*S**2. + L2c*Delta**2 \
        + L2d*delta_m**2. + L3a*Delta*S*delta_m + L3b*S*Delta**2. + L3c*S**3. \
        + L3d*S*delta_m**2. + L4a*Delta*S**2*delta_m + L4b*Delta**3.*delta_m \
        + L4c*Delta**4. + L4d*S**4. + L4e*Delta**2.*S**2. + L4f*delta_m**4 + L4g*Delta*delta_m**3. \
        + L4h*Delta**2.*delta_m**2. + L4i*S**2.*delta_m**2.) \
        + S*(1. + 8.*eta)*delta_m**4. + eta*J_isco*delta_m**6.
    
    return abs(a_f-a_f_new)

def bbh_final_spin_non_precessing_Healyetal(m1, m2, chi1, chi2):
    """
    Calculate the spin of the final BH resulting from the merger of two black holes with non-precessing spins using fit from Healy et al Phys Rev D 90, 104004 (2014)
    
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
        return np.vectorize(bbh_final_spin_non_precessing_Healyetal)(m1, m2, chi1, chi2)
    
    # binary parameters
    m = m1+m2
    q = m1/m2
    eta = q/(1.+q)**2.
    delta_m = (m1-m2)/m
    
    S1 = chi1*m1**2 # spin angular momentum 1
    S2 = chi2*m2**2 # spin angular momentum 2
    S = (S1+S2)/m**2 # symmetric spin (dimensionless -- called \tilde{S} in the paper)
    Delta = (S2/m2-S1/m1)/m # antisymmetric spin (dimensionless -- called tilde{Delta} in the paper
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # compute the final spin
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    
    x, cov_x = so.leastsq(_final_spin_diff_Healyetal, 0., args=(eta, delta_m, S, Delta))
    chif = x[0]
    
    return chif

def bbh_final_mass_non_precessing_Healyetal(m1, m2, chi1, chi2, chif=None):
    """
    Calculate the mass of the final BH resulting from the merger of two black holes with non-precessing spins using fit from Healy et al Phys Rev D 90, 104004 (2014)
    
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
    
    # binary parameters
    m = m1+m2
    q = m1/m2
    eta = q/(1.+q)**2.
    delta_m = (m1-m2)/m
    
    S1 = chi1*m1**2 # spin angular momentum 1
    S2 = chi2*m2**2 # spin angular momentum 2
    S = (S1+S2)/m**2 # symmetric spin (dimensionless -- called \tilde{S} in the paper)
    Delta = (S2/m2-S1/m1)/m # antisymmetric spin (dimensionless -- called tilde{Delta} in the paper
    
    if chif is None:
        chif = bbh_final_spin_non_precessing_Healyetal(m1, m2, chi1, chi2)
    else:
        chif = np.array(chif)
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # now compute the final mass
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    r_isco = calc_isco_radius(chif)
    
    # fitting coefficients - Table XI of Healy et al Phys Rev D 90, 104004 (2014)
    # [fourth order fits]
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
    
    # binding energy at ISCO -- Eq.(2.7) of Ori, Thorne Phys Rev D 62 124022 (2000)
    E_isco = (1. - 2./r_isco + chif/r_isco**1.5)/np.sqrt(1. - 3./r_isco + 2.*chif/r_isco**1.5)
    
    # final mass -- Eq. (14) of Healy et al Phys Rev D 90, 104004 (2014)
    mf = (4.*eta)**2*(M0 + K1*S + K2a*Delta*delta_m + K2b*S**2 + K2c*Delta**2 + K2d*delta_m**2 \
        + K3a*Delta*S*delta_m + K3b*S*Delta**2 + K3c*S**3 + K3d*S*delta_m**2 \
        + K4a*Delta*S**2*delta_m + K4b*Delta**3*delta_m + K4c*Delta**4 + K4d*S**4 \
        + K4e*Delta**2*S**2 + K4f*delta_m**4 + K4g*Delta*delta_m**3 + K4h*Delta**2*delta_m**2 \
        + K4i*S**2*delta_m**2) + (1+eta*(E_isco+11.))*delta_m**6.
    
    return mf*m

def bbh_final_spin_projected_spin_Healyetal(m1, m2, chi1, chi2, tilt1, tilt2):
    """
    Calculate the spin of the final BH resulting from the merger of two black holes projecting spin along angular momentum and using fit from Healy et al Phys Rev D 90, 104004 (2014)
    
    Parameters
    ----------
    m1, m2 : component masses
    chi1, chi2 : dimensionless spins of two BHs
    tilt1, tilt2 : tilts (in radians) in the new spin convention
    
    Returns
    -------
    final spin, chif
    """
    return bbh_final_spin_non_precessing_Healyetal(m1, m2, chi1*np.cos(tilt1), chi2*np.cos(tilt2))

def bbh_final_mass_projected_spin_Healyetal(m1, m2, chi1, chi2, tilt1, tilt2, chif=None):
    """
    Calculate the mass of the final BH resulting from the merger of two black holes projecting spin along angular momentum and using fit from Healy et al Phys Rev D 90, 104004 (2014)
    
    Parameters
    ----------
    m1, m2 : component masses
    chi1, chi2 : dimensionless spins of two BHs
    tilt1, tilt2 : tilts (in radians) in the new spin convention
    chif: final spin (optional), if already calculated
    
    Returns
    -------
    final mass, mf
    """
    return bbh_final_mass_non_precessing_Healyetal(m1, m2, chi1*np.cos(tilt1), chi2*np.cos(tilt2), chif=chif)

def bbh_final_spin_precessing_Healyetal_extension_Mf(m1, m2, chi1, chi2, tilt1, tilt2, phi12):
    """
    Calculate the spin of the final BH resulting from the merger of two black holes including the in-plane spine by projecting the spins along the angular momentum and using fit from Healy et al Phys Rev D 90, 104004 (2014) and then adding in quadrature the in-plane dimensionful spin scaled by the final mass squared.
    
    Parameters
    ----------
    m1, m2 : component masses
    chi1, chi2 : dimensionless spins of two BHs
    tilt1, tilt2 : tilts (in radians) in the new spin convention
    phi12: angle (in radians) between in-plane spin components

    Returns
    -------
    final spin, chif
    """

    # First compute the final mass and parallel component of the final spin using the aligned components of the initial spins
    chifpara = bbh_final_spin_non_precessing_Healyetal(m1, m2, chi1*np.cos(tilt1), chi2*np.cos(tilt2))
    Mf = bbh_final_mass_non_precessing_Healyetal(m1, m2, chi1*np.cos(tilt1), chi2*np.cos(tilt2), chifpara)

    # Now compute the squared magnitude of the in-plane dimensionful spin, first computing the magnitudes of the initial in-plane spins

    S1perpmag = m1*m1*chi1*np.sin(tilt1)
    S2perpmag = m2*m2*chi2*np.sin(tilt2)

    Sperpmag2 = S1perpmag*S1perpmag + S2perpmag*S2perpmag + 2.*S1perpmag*S2perpmag*np.cos(phi12)

    # Combine together and return

    return (chifpara*chifpara + Sperpmag2/Mf**4.)**0.5

def bbh_final_spin_precessing_Healyetal_extension_Minit(m1, m2, chi1, chi2, tilt1, tilt2, phi12):
    """
    Calculate the spin of the final BH resulting from the merger of two black holes including the in-plane spine by projecting the spins along the angular momentum and using fit from Healy et al Phys Rev D 90, 104004 (2014) and then adding in quadrature the in-plane dimensionful spin scaled by the initial mass squared.
    
    Parameters
    ----------
    m1, m2 : component masses
    chi1, chi2 : dimensionless spins of two BHs
    tilt1, tilt2 : tilts (in radians) in the new spin convention
    phi12: angle (in radians) between in-plane spin components

    Returns
    -------
    final spin, chif
    """

    # First compute the final mass and parallel component of the final spin using the aligned components of the initial spins
    chifpara = bbh_final_spin_non_precessing_Healyetal(m1, m2, chi1*np.cos(tilt1), chi2*np.cos(tilt2))

    # Now compute the squared magnitude of the in-plane dimensionful spin, first computing the magnitudes of the initial in-plane spins

    S1perpmag = m1*m1*chi1*np.sin(tilt1)
    S2perpmag = m2*m2*chi2*np.sin(tilt2)

    Sperpmag2 = S1perpmag*S1perpmag + S2perpmag*S2perpmag + 2.*S1perpmag*S2perpmag*np.cos(phi12)

    # Combine together and return

    return (chifpara*chifpara + Sperpmag2/(m1+m2)**4.)**0.5

def bbh_final_mass_non_precessing_Husaetal(m1, m2, chi1, chi2): 
    """ 
    Calculate the mass and spin of the final BH resulting from the 
    merger of two black holes with non-precessing spins using the fits
    used by IMRPhenomD, given in Eqs. (3.6) and (3.8) of Husa et al.
    Phys Rev D 93, 044006 (2016). Note that Eq. (3.8) gives the radiated energy, not
    the final mass directly

    m1, m2: component masses
    chi1, chi2: dimensionless spins of two BHs
    """
    m1 = np.vectorize(float)(np.array(m1))
    m2 = np.vectorize(float)(np.array(m2))
    chi1 = np.vectorize(float)(np.array(chi1))
    chi2 = np.vectorize(float)(np.array(chi2))
    
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
    Calculate the mass and spin of the final BH resulting from the 
    merger of two black holes with non-precessing spins using the fits
    used by IMRPhenomD, given in Eqs. (3.6) and (3.8) of Husa et al.
    Phys Rev D 93, 044006 (2016). Note that Eq. (3.8) gives the radiated energy, not
    the final mass directly

    m1, m2: component masses
    chi1, chi2: dimensionless spins of two BHs
    """
    # Vectorize the function if arrays are provided as input
    m1 = np.vectorize(float)(np.array(m1))
    m2 = np.vectorize(float)(np.array(m2))
    chi1 = np.vectorize(float)(np.array(chi1))
    chi2 = np.vectorize(float)(np.array(chi2))
    
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

def bbh_final_spin_precessing_Husaetal_extension_Minit_chi_p(m1, m2, chi1z, chi2z, chi_p):
    """
    Calculate the spin of the final BH resulting from the merger of two black holes including the in-plane spine by projecting the spins along the angular momentum and using fit from Husa et al Phys Rev D 93, 044006 (2016) and then adding in quadrature the in-plane dimensionful spin scaled by the final mass squared. Here we compute the in-plane spin from chi_p as is done in FinalSpinIMRPhenomD_all_in_plane_spin_on_larger_BH in LALSimIMRPhenomP.c
    
    Parameters
    ----------
    m1, m2 : component masses
    chi1z, chi2z : components of the dimensionless spins of the two BHs along the orbital angular momentum
    chi_p: IMRPhenomP in-plane spin parameter
 
    Returns
    -------
    final spin, chif
    """

    # Compute the component of the spin parallel to the orbital angular momentum using the IMRPhenomD fit

    chifpara = bbh_final_spin_non_precessing_Husaetal(m1, m2, chi1z, chi2z)

    # Compute the ratio of the more massive black hole to the total mass

    qf = np.maximum(m1, m2)/(m1 + m2)

    # Compute the (already scaled) in-plane spin from qf and chi_p

    Splane = qf*qf*chi_p

    # Combine together and return

    return (chifpara*chifpara + Splane*Splane)**0.5

def bbh_final_spin_precessing_Husaetal_extension_Minit(m1, m2, chi1, chi2, tilt1, tilt2, phi12):
    """
    Calculate the spin of the final BH resulting from the merger of two black holes including the in-plane spine by projecting the spins along the angular momentum and using the fit from Husa et al Phys Rev D 93, 044006 (2016) and then adding in quadrature the in-plane dimensionful spin scaled by the initial mass squared.
    
    Parameters
    ----------
    m1, m2 : component masses
    chi1, chi2 : dimensionless spins of two BHs
    tilt1, tilt2 : tilts (in radians) in the new spin convention
    phi12: angle (in radians) between in-plane spin components

    Returns
    -------
    final spin, chif
    """

    # First compute the final mass and parallel component of the final spin using the aligned components of the initial spins
    chifpara = bbh_final_spin_non_precessing_Husaetal(m1, m2, chi1*np.cos(tilt1), chi2*np.cos(tilt2))

    # Now compute the squared magnitude of the in-plane dimensionful spin, first computing the magnitudes of the initial in-plane spins

    S1perpmag = m1*m1*chi1*np.sin(tilt1)
    S2perpmag = m2*m2*chi2*np.sin(tilt2)

    Sperpmag2 = S1perpmag*S1perpmag + S2perpmag*S2perpmag + 2.*S1perpmag*S2perpmag*np.cos(phi12)

    # Combine together and return

    return (chifpara*chifpara + Sperpmag2/(m1+m2)**4.)**0.5

def bbh_final_spin_precessing_Barausse_and_Rezzolla(m1, m2, a1, a2, tilt1, tilt2, phi12): 
	""" 
	Calculate the dimensionless spin of the final BH resulting from the 
	merger of two black holes with precessing spins using the fit from Barausse and Rezzolla ApJL 704, L40 (2009). We base our implementation on the IMRPhenomPv2 one in FinalSpinBarausse2009 in LALSimIMRPhenomP.c. (We duplicate this here as a stop-gap until that function can be made public for direct use.)

	m1, m2: component masses (with m1 > m2)
	a1, a2: dimensionless spins of two BHs
	tilt1, tilt2: tilt angles of the spins from the orbital angular momentum
	phi12: angle between in-plane spin components 
	"""

        # Vectorize the function if arrays are provided as input
        m1 = np.vectorize(float)(np.array(m1))
        m2 = np.vectorize(float)(np.array(m2))
        a1 = np.vectorize(float)(np.array(a1))
        a2 = np.vectorize(float)(np.array(a2))
        tilt1 = np.vectorize(float)(np.array(tilt1))
        tilt2 = np.vectorize(float)(np.array(tilt2))
        phi12 = np.vectorize(float)(np.array(phi12))

	# Coefficients

	s4 = -0.1229
	s5 = 0.4537
	t0 = -2.8904
	t2 = -3.5171
	t3 = 2.5763

	# Computing angles
	cos_beta_tilde = np.cos(tilt1)
	cos_gamma_tilde = np.cos(tilt2)
	cos_alpha = ((1 - cos_beta_tilde*cos_beta_tilde)*(1 - cos_gamma_tilde*cos_gamma_tilde))**0.5*np.cos(phi12) + cos_beta_tilde*cos_gamma_tilde

	# Definitions
	q = m2/m1
	nu = m1*m2/(m1+m2)**2

	# Shorthands
	nu2 = nu*nu
	q2 = q*q
	q4 = q2*q2
	q2p = 1. + q2
	q2p2 = q2p*q2p
	qp = 1. + q
	qp2 = qp*qp
	a1_2 = a1*a1
	a2_2 = a2*a2

	# Compute the final spin and return it

	l = 2.*3.**0.5 + t2*nu + t3*nu2 + (s4 / q2p2) * (a1_2 + a2_2*q4 + 2.*a1*a2*q2*cos_alpha) + ((s5*nu + t0 + 2.)/q2p) * (a1*cos_beta_tilde + a2*cos_gamma_tilde*q2)
	
	return (1. / qp2) * (a1_2 + a2_2*q4 + 2.*a1*a2*q2*cos_alpha + 2.*(a1*cos_beta_tilde + a2*q2*cos_gamma_tilde)*l*q + l*l*q2)**0.5

def qnmfreqs_berti(a, l, m, n):
    """
    compute the (complex) QNM frequencies for different 
    overtones using the fits provided by Berti, Cardoso, Will 
    
    a         : Kerr parameter of the BH 
    l, m, n    : indices of spheroidal harnonics modes 
    """

    # load the data file containing the fits (Berti et al,  gr-qc/0512160)
    lVec, mVec, nVec, f1Vec, f2Vec, f3Vec, q1Vec, q2Vec, q3Vec = np.loadtxt('../src/Berti_QNMfitcoeffsWEB.dat', unpack=True)

    idx = np.logical_and(np.logical_and(lVec == l, mVec == m), nVec == n)

    # evaluate the Berti et al fits to the complex frequencies 
    if len(lVec[idx]) == 1:

         f1 = f1Vec[idx]
         f2 = f2Vec[idx]
         f3 = f3Vec[idx]
         q1 = q1Vec[idx]
         q2 = q2Vec[idx]
         q3 = q3Vec[idx]

         omega = (f1 + f2*(1.-a)**f3) # fit to omega
         Q = q1 + q2*(1.-a)**q3       # fit to quality factor 
         tau = 2.*Q/omega

         return omega - 1j/tau  # complex frequency 

    elif len(lVec[idx]) < 1:
            print '# No data matching this l, m, n combination (l = #d m = #d n = #d)' %(l, m, n)
            exit()
    else:
            print '# More than on fit point corresponding to this l, m, n combination'
            exit()

def calc_fqnm_dominant_mode(chif): 
    """
    calculate the (dimensionless) freq of the dominant QNM of a Kerr BH with spin chif 
    """

    Omega = qnmfreqs_berti(chif, 2, 2, 0)
    return np.real(Omega)/(2*np.pi) 


