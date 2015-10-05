""" 
Various fitting formulas provided by numerical relativity 

P. Ajith, 2015-04-09 
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

	# fitting coefficients - Table X1 of Healy et al Phys Rev D 90, 104004 (2014) 
	# [forth order fits]
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
	merger of two black holes with non-precessing spins

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

def qnmfreqs_berti(a, l, m, n):
	"""
	compute the (complex) QNM frequencies for different 
	overtones using the fits provided by Berti, Cardoso, Will 
	
	a     	: Kerr parameter of the BH 
	l, m, n	: indices of spheroidal harnonics modes 
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

def calc_fqnm_dominant_mode(af): 
	"""
	calculate the (dimensionless) freq of the dominant QNM of a Kerr BH with spin af 
	"""

	Omega = qnmfreqs_berti(af, 2, 2, 0)
	return np.real(Omega)/(2*np.pi) 


