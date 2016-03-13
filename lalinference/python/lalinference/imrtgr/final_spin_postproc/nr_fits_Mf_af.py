""" 
Various fitting formulas provided by numerical relativity 

P. Ajith, 2015-04-09
N. K. Johnson-McDaniel, 10.2015 - 03.2016 
"""

import numpy as np 
import scipy.optimize as so

def calc_isco_radius(a): 
	""" 
	Calculate the ISCO radius of a Kerr BH as a function of the Kerr parameter 
	
	a : Kerr parameter (numpy array) 
	"""

	# Ref. Eq. (2.5) of Ori, Thorne Phys Rev D 62 124022 (2000)
	z1 = 1.+(1.-a**2.)**(1./3)*((1.+a)**(1./3) + (1.-a)**(1./3))
	z2 = np.sqrt(3.*a**2 + z1**2)
	a_sign = np.sign(a)
	return 3+z2 - np.sqrt((3.-z1)*(3.+z1+2.*z2))*a_sign

def calc_isco_freq(a): 
	""" 
	Calculate the ISCO frequency of a Kerr BH as a function of the Kerr parameter 
	
	a : Kerr parameter (numpy array) 
	"""
	
	r_isco = calc_isco_radius(a)
	u_isco = r_isco**-0.5	
	v_isco = u_isco*(1.-a*u_isco**3.+a**2.*u_isco**6)**(1./3.)
	return v_isco**3./np.pi
	
def _final_spin_diff(a_f, eta, delta_m, S, Delta): 
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

	return np.ravel(abs(a_f-a_f_new))


def bbh_final_mass_and_spin_non_precessing(m1, m2, chi1, chi2): 
	""" 
	Calculate the mass and spin of the final BH resulting from the 
	merger of two black holes with non-precessing spins using the fit from
	Healy, Lousto, and Zlochower PRD 90, 104004 (2014)

	m1, m2: component masses
	chi1, chi2: dimensionless spins of two BHs
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(m1, '__len__'):
		return np.vectorize(bbh_final_mass_and_spin_non_precessing)(m1, m2, chi1, chi2)
	
	# binary parameters 
	m = m1+m2  
	q = m1/m2 
	eta = q/(1.+q)**2. 
	delta_m = (m1-m2)/m

	S1 = chi1*m1**2				# spin angular momentum 1 
	S2 = chi2*m2**2				# spin angular momentum 2 
	S = (S1+S2)/m**2 			# symmetric spin (dimensionless -- called \tilde{S} in the paper) 
	Delta = (S2/m2-S1/m1)/m		# antisymmetric spin (dimensionless -- called tilde{Delta} in the paper

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
	# compute the final spin 
	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
	# res = so.minimize_scalar(_final_spin_diff, bounds=(-0.99999, 0.99999), args=(eta, delta_m, S, Delta), method='Bounded', tol=1e-6, options={'maxiter':100, 'disp':False})
	# a_f = res.x

	x, cov_x = so.leastsq(_final_spin_diff, 0., args=(eta, delta_m, S, Delta))
	a_f = x[0]

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
	# now compute the final mass  
	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
	r_isco = calc_isco_radius(a_f)

	# fitting coefficients - Table X1 of Healy et al Phys Rev D 90, 104004 (2014) 
	# [forth order fits]
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
	E_isco = (1. - 2./r_isco + a_f/r_isco**1.5)/np.sqrt(1. - 3./r_isco + 2.*a_f/r_isco**1.5)

	# final mass -- Eq. (14) of Healy et al Phys Rev D 90, 104004 (2014)
	mf = (4.*eta)**2*(M0 + K1*S + K2a*Delta*delta_m + K2b*S**2 + K2c*Delta**2 + K2d*delta_m**2 \
		+ K3a*Delta*S*delta_m + K3b*S*Delta**2 + K3c*S**3 + K3d*S*delta_m**2 \
		+ K4a*Delta*S**2*delta_m + K4b*Delta**3*delta_m + K4c*Delta**4 + K4d*S**4 \
		+ K4e*Delta**2*S**2 + K4f*delta_m**4 + K4g*Delta*delta_m**3 + K4h*Delta**2*delta_m**2 \
		+ K4i*S**2*delta_m**2) + (1+eta*(E_isco+11.))*delta_m**6.
	
	return mf*m, a_f  

def bbh_final_mass_and_spin_non_precessing_Husaetal(m1, m2, chi1, chi2): 
	""" 
	Calculate the mass and spin of the final BH resulting from the 
	merger of two black holes with non-precessing spins using the fits
	used by IMRPhenomD, given in Eqs. (3.6) and (3.8) of Husa et al.
	arXiv:1508.07250. Note that Eq. (3.8) gives the radiated energy, not
	the final mass directly

	m1, m2: component masses
	chi1, chi2: dimensionless spins of two BHs
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(m1, '__len__'):
		return np.vectorize(bbh_final_mass_and_spin_non_precessing_Husaetal)(m1, m2, chi1, chi2)
	
	# binary parameters 
	m = m1+m2  
        msq = m*m

	eta = m1*m2/msq
	eta2 = eta*eta
	eta3 = eta2*eta
	eta4 = eta3*eta

	S1 = chi1*m1**2/msq	       		# spin angular momentum 1 (in m = 1 units) 
	S2 = chi2*m2**2/msq		       	# spin angular momentum 2 (in m = 1 units)
	S = S1+S2		    	        # total spin
	Sh = S/(1. - 2.*eta)                    # rescaled total spin

	S2 = S*S
	S3 = S2*S
	S4 = S3*S

	# Expressions copied from LALSimIMRPhenomD_internals.c (except with two notation differences: S is capitalized in chif and s -> Sh in Mf, in addition to the "m*(1. - ...)" to obtain the final mass from the radiated mass in m = 1 units which is calculated in the LAL code)

	Mf = m*(1. - ((0.055974469826360077*eta + 0.5809510763115132*eta2 - 0.9606726679372312*eta3 + 3.352411249771192*eta4)*
    (1. + (-0.0030302335878845507 - 2.0066110851351073*eta + 7.7050567802399215*eta2)*Sh))/(1. + (-0.6714403054720589 - 1.4756929437702908*eta + 7.304676214885011*eta2)*Sh))

        chif = 3.4641016151377544*eta - 4.399247300629289*eta2 + 9.397292189321194*eta3 - 13.180949901606242*eta4 + (1 - 0.0850917821418767*eta - 5.837029316602263*eta2)*S + (0.1014665242971878*eta - 2.0967746996832157*eta2)*S2 + (-1.3546806617824356*eta + 4.108962025369336*eta2)*S3 + (-0.8676969352555539*eta + 2.064046835273906*eta2)*S4 

	return Mf, chif

def bbh_final_mass_and_spin_precessing_IMRPhenomPv2(m1, m2, chi1z, chi2z, chi_p): 
	""" 
	Calculate the mass and dimensionless spin of the final BH resulting from the 
	merger of two black holes with precessing spins using the modification
	of the IMRPhenomD fit used by IMRPhenomPv2, given in
	FinalSpinIMRPhenomD_all_in_plane_spin_on_larger_BH in LALSimIMRPhenomP.c
	plus the final mass fit with just the aligned components of the spins,
	again as in IMRPhenomPv2

	m1, m2: component masses
	chi1z, chi2z: components of the dimensionless spins of the two BHs along the orbital angular momentum
	chi_p: IMRPhenomP in-plane spin parameter
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(m1, '__len__'):
		return np.vectorize(bbh_final_mass_and_spin_precessing_IMRPhenomPv2)(m1, m2, chi1z, chi2z, chi_p)

	# Compute the ratio of the more massive black hole to the total mass

	qf = max(m1, m2)/(m1 + m2)

	# Compute the component of the spin parallel to the orbital angular momentum using the IMRPhenomD fit

	mf, chif_par = bbh_final_mass_and_spin_non_precessing_Husaetal(m1, m2, chi1z, chi2z)

	# Compute the in-plane spin with the scaling from the mass ratio

	Splane = qf*qf*chi_p

	return mf, (Splane*Splane + chif_par*chif_par)**0.5

def bbh_final_mass_and_spin_precessing_IMRPhenomP_tilts(m1, m2, chi1, chi2, tilt1, tilt2, phi12):
	""" 
	Calculate the spin of the final BH resulting from the 
	merger of two black holes with precessing spins using
	the Healy, Lousto, and Zlochower fit for the aligned part
	and including the contribution from the in-plane spin components a la IMRPhenomP (same as in bbh_final_mass_and_spin_HLZ_extension_precessing except with Mf -> m1 + m2)

	m1, m2: component masses
	chi1, chi2: dimensionless spins of two BHs
	tilt1, tilt2: tilt angles of the spins from the orbital angular momentum
	phi12: angle between in-plane spin components 
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(m1, '__len__'):
		return np.vectorize(bbh_final_mass_and_spin_precessing_IMRPhenomP_tilts)(m1, m2, chi1, chi2, tilt1, tilt2, phi12)

	# First compute the final mass and parallel component of the final spin using the aligned components of the initial spins
	Mf, afpara = bbh_final_mass_and_spin_non_precessing_Husaetal(m1, m2, chi1*np.cos(tilt1), chi2*np.cos(tilt2))

	# Now compute the magnitude of the in-plane final dimensionful spin, first computing the magnitudes of the initial in-plane spins

	S1perpmag = m1*m1*chi1*np.sin(tilt1)
	S2perpmag = m2*m2*chi2*np.sin(tilt2)

	Sperpmag = (S1perpmag*S1perpmag + S2perpmag*S2perpmag + 2.*S1perpmag*S2perpmag*np.cos(phi12))**0.5

	# Put everything together

	return Mf, (afpara*afpara + Sperpmag*Sperpmag/(m1+m2)**4.)**0.5

def bbh_final_mass_and_spin_HLZ_extension_precessing(m1, m2, chi1, chi2, tilt1, tilt2, phi12):
	""" 
	Calculate the spin of the final BH resulting from the 
	merger of two black holes with precessing spins using
	the Healy, Lousto, and Zlochower fit for the aligned part
	and including the contribution from the in-plane spin components

	m1, m2: component masses
	chi1, chi2: dimensionless spins of two BHs
	tilt1, tilt2: tilt angles of the spins from the orbital angular momentum
	phi12: angle between in-plane spin components 
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(m1, '__len__'):
		return np.vectorize(bbh_final_mass_and_spin_HLZ_extension_precessing)(m1, m2, chi1, chi2, tilt1, tilt2, phi12)

	# First compute the final mass and parallel component of the final spin using the aligned components of the initial spins
	Mf, afpara = bbh_final_mass_and_spin_non_precessing(m1, m2, chi1*np.cos(tilt1), chi2*np.cos(tilt2))

	# Now compute the magnitude of the in-plane final dimensionful spin, first computing the magnitudes of the initial in-plane spins

	S1perpmag = m1*m1*chi1*np.sin(tilt1)
	S2perpmag = m2*m2*chi2*np.sin(tilt2)

	Sperpmag = (S1perpmag*S1perpmag + S2perpmag*S2perpmag + 2.*S1perpmag*S2perpmag*np.cos(phi12))**0.5

	# Put everything together

	return Mf, (afpara*afpara + Sperpmag*Sperpmag/Mf**4.)**0.5

def bbh_final_mass_and_spin_HLZ_extension_precessing_IMRPhenomP(m1, m2, chi1, chi2, tilt1, tilt2, phi12):
	""" 
	Calculate the spin of the final BH resulting from the 
	merger of two black holes with precessing spins using
	the Healy, Lousto, and Zlochower fit for the aligned part
	and including the contribution from the in-plane spin components a la IMRPhenomP (same as in bbh_final_mass_and_spin_HLZ_extension_precessing except with Mf -> m1 + m2)

	m1, m2: component masses
	chi1, chi2: dimensionless spins of two BHs
	tilt1, tilt2: tilt angles of the spins from the orbital angular momentum
	phi12: angle between in-plane spin components 
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(m1, '__len__'):
		return np.vectorize(bbh_final_mass_and_spin_HLZ_extension_precessing_IMRPhenomP)(m1, m2, chi1, chi2, tilt1, tilt2, phi12)

	# First compute the final mass and parallel component of the final spin using the aligned components of the initial spins
	Mf, afpara = bbh_final_mass_and_spin_non_precessing(m1, m2, chi1*np.cos(tilt1), chi2*np.cos(tilt2))

	# Now compute the magnitude of the in-plane final dimensionful spin, first computing the magnitudes of the initial in-plane spins

	S1perpmag = m1*m1*chi1*np.sin(tilt1)
	S2perpmag = m2*m2*chi2*np.sin(tilt2)

	Sperpmag = (S1perpmag*S1perpmag + S2perpmag*S2perpmag + 2.*S1perpmag*S2perpmag*np.cos(phi12))**0.5

	# Put everything together

	return Mf, (afpara*afpara + Sperpmag*Sperpmag/(m1+m2)**4.)**0.5

def bbh_final_spin_precessing_Barausse_and_Rezzolla(m1, m2, a1, a2, tilt1, tilt2, phi12): 
	""" 
	Calculate the dimensionless spin of the final BH resulting from the 
	merger of two black holes with precessing spins using the fit from Barausse and Rezzolla ApJL 704, L40 (2009). We base our implementation on the IMRPhenomPv2 one in FinalSpinBarausse2009 in LALSimIMRPhenomP.c.

	m1, m2: component masses (with m1 > m2)
	chi1, chi2: dimensionless spins of two BHs
	tilt1, tilt2: tilt angles of the spins from the orbital angular momentum
	phi12: angle between in-plane spin components 
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(m1, '__len__'):
		return np.vectorize(bbh_final_spin_precessing_Barausse_and_Rezzolla)(m1, m2, a1, a2, tilt1, tilt2, phi12)

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

def bbh_final_mass_and_spin_precessing_HLZ_mass_Barausse_and_Rezzolla_spin(m1, m2, a1, a2, tilt1, tilt2, phi12): 
	""" 
	Wrapper function to combine together the HLZ final mass with the Barausse and Rezzolla final spin 

	m1, m2: component masses (with m1 > m2)
	chi1, chi2: dimensionless spins of two BHs
	tilt1, tilt2: tilt angles of the spins from the orbital angular momentum
	phi12: angle between in-plane spin components 
	"""

	mf_HLZ, af_HLZ = bbh_final_mass_and_spin_non_precessing(m1, m2, a1*np.cos(tilt1), a2*np.cos(tilt2))

	af_BR = bbh_final_spin_precessing_Barausse_and_Rezzolla(m1, m2, a1, a2, tilt1, tilt2, phi12)

	return mf_HLZ, af_BR
