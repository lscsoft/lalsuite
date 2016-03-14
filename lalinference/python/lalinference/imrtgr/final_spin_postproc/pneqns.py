import dE_by_Flux
from scipy.integrate import ode
import numpy as np 

MTSUN_SI = 4.926e-6


def precession_eqns(t, y_vec, m1, m2): 
	""" The coupled set of ODEs containing the PN precession eqns as well as the evolution of dv/dt 
	    All these equations are listed in Ajith (2011) (http://arxiv.org/pdf/1107.1267v2.pdf).
	"""

	# unpack the vector of variables 
	v, S1x, S1y, S1z, S2x, S2y, S2z, Lx, Ly, Lz = y_vec 

	m1_sqr = (m1*MTSUN_SI)**2. #(m1*lal.MTSUN_SI)**2.
	m2_sqr = (m2*MTSUN_SI)**2. #(m2*lal.MTSUN_SI)**2.
	m = (m1+m2)*MTSUN_SI #*lal.MTSUN_SI
	q = m2/m1
	q1 = 1./q
	eta = m1*m2/(m1+m2)**2 
	delta = (m1-m2)/(m1+m2)

	# spin and angular momentum vectors 
	Ln = np.array([Lx, Ly, Lz])
	S1 = np.array([S1x, S1y, S1z])
	S2 = np.array([S2x, S2y, S2z])
	chi1 = S1/m1_sqr 
	chi2 = S2/m2_sqr 
	chi_s = (chi1+chi2)/2.
	chi_a = (chi1-chi2)/2.

	#print 'v = ', v , ' chi1 = ', chi1, ' chi2 = ', chi2 

	# dot products 
	chiadL = np.dot(chi_a, Ln) 
	chisdL = np.dot(chi_s, Ln) 
	chissqr = np.dot(chi_s, chi_s) 
	chiasqr = np.dot(chi_a, chi_a) 
	chisdchia = np.dot(chi_s, chi_a) 
	S1dL = np.dot(S1, Ln) 
	S2dL = np.dot(S2, Ln) 

	# magnitude of the orb ang momentum. The prefactor of dLN/dt in Eq.(3.6) of Ajith (2011)
	Lmag = (eta*m**2./v)*(1.+(1.5+eta/6.)*v**2.)

	# precession freqs. Eqs.(3.8) of Ajith (2011) 
	Omega1=(v**5./m)*((3./4+eta/2.-3.*delta/4.)*Ln+0.5*v/m**2.*(-3.0*(S2dL +q*S1dL)*Ln+S2)+v**2.*(9./16.+5.*eta/4.-eta**2./24.-9.*delta/16.+5.*delta*eta/8.)*Ln)
	Omega2=(v**5./m)*((3./4+eta/2.+3.*delta/4.)*Ln+0.5*v/m**2.*(-3.0*(q1*S2dL+S1dL)*Ln+S1)+v**2.*(9./16.+5.*eta/4.-eta**2./24.+9.*delta/16.-5.*delta*eta/8.)*Ln)

	# spin precession eqns. Eqs.(3.7) of Ajith (2011)
	dS1_dt = np.cross(Omega1, S1)
	dS2_dt = np.cross(Omega2, S2)

	# ang momentum precession eqn. Eqs.(3.6) of Ajith (2011)
	dL_dt =-(dS1_dt+dS2_dt)/Lmag

	# orb frequency evolution. Eqs.(3.5) of Ajith (2011)
	dv_dt = -1./(m*dE_by_Flux.denergy_by_flux(v, eta, delta, chiadL, chisdL, chiasqr, chissqr, chisdchia, 7))

	return  [dv_dt, dS1_dt[0], dS1_dt[1], dS1_dt[2], dS2_dt[0], dS2_dt[1], dS2_dt[2], dL_dt[0], dL_dt[1], dL_dt[2]]
	

def evolve_spins_dt(v0, m1, m2, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, Lnx, Lny, Lnz, v_final, dt): 
	""" evolve the spins and orb ang momentum according to the PN precession eqns."""

	m1_sqr = (m1*MTSUN_SI)**2. #(m1*lal.MTSUN_SI)**2.
	m2_sqr = (m2*MTSUN_SI)**2. #(m2*lal.MTSUN_SI)**2.

	# pack the inputs 
	y_vec0 = [v0, chi1x*m1_sqr, chi1y*m1_sqr, chi1z*m1_sqr, chi2x*m2_sqr, chi2y*m2_sqr, chi2z*m2_sqr, Lnx, Lny, Lnz]
	t0 = 0.

	T_MAX = 1e7
	R_TOL = 1e-6 
	A_TOL = 1e-6 

	# initialize the ODE solver 
	backend = "dopri5"
	solver = ode(precession_eqns)
	solver.set_integrator(backend, atol=R_TOL, rtol=R_TOL)  # nsteps=1
	solver.set_initial_value(y_vec0, t0)
	solver.set_f_params(m1, m2)

	y_result = []
	t_output = []
	y_result.append(y_vec0)
	t_output.append(t0)

	# evolve the eqns 
	while solver.successful() and solver.t < T_MAX and solver.y[0] <= v_final: 
			
		solver.integrate(solver.t + dt, step=1)
		y_result.append(solver.y)
		t_output.append(solver.t)

		#print '... t = ', solver.t, ' y = ', solver.y 
			
	Y = np.array(y_result)
	t = np.array(t_output)
	v_v, S1x_v, S1y_v, S1z_v, S2x_v, S2y_v, S2z_v, Lx_v, Ly_v, Lz_v  = Y.T

	return v_v, S1x_v/m1_sqr, S1y_v/m1_sqr, S1z_v/m1_sqr, S2x_v/m2_sqr, S2y_v/m2_sqr, S2z_v/m2_sqr, Lx_v, Ly_v, Lz_v

def find_tilts_and_phi12_at_freq(v0, m1, m2, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, Lnx, Lny, Lnz, v_final, dt):
        """ given the spins and ang momentum at a given frequency, find the tilt and in-plane spin angles at a later frequency """

	print("v0 = %f, m1 = %f, m2 = %f, chi1x = %f, chi1y = %f, chi1z = %f, chi2x = %f, chi2y = %f, chi2z = %f"%(v0, m1, m2, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z))

        # evolve the spins                                                                                                                                                           
        v_v, chi1x_v, chi1y_v, chi1z_v, chi2x_v, chi2y_v, chi2z_v, Lnx_v, Lny_v, Lnz_v = evolve_spins_dt(v0, m1, m2, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, Lnx, Lny, Lnz, v_final, dt)

	chi1x_v = chi1x_v[-1]
	chi1y_v = chi1y_v[-1]
	chi1z_v = chi1z_v[-1]

	chi2x_v = chi2x_v[-1]
	chi2y_v = chi2y_v[-1]
	chi2z_v = chi2z_v[-1]

	Lnx_v = Lnx_v[-1]
	Lny_v = Lny_v[-1]
	Lnz_v = Lnz_v[-1]

        # norms and normalizing                                                                                                                    
        chi1_norm = np.sqrt(chi1x_v**2+chi1y_v**2+chi1z_v**2)
        chi2_norm = np.sqrt(chi2x_v**2+chi2y_v**2+chi2z_v**2)
        L_norm = np.sqrt(Lnx_v**2+Lny_v**2+Lnz_v**2)

	Lnx_v /= L_norm
	Lny_v /= L_norm
	Lnz_v /= L_norm

   # dot products                                                                                                                                  
        chi1dL_v = chi1x_v*Lnx_v + chi1y_v*Lny_v + chi1z_v*Lnz_v
        chi2dL_v = chi2x_v*Lnx_v + chi2y_v*Lny_v + chi2z_v*Lnz_v

	# in-plane spins

	chi1inplanex = chi1x_v - chi1dL_v*Lnx_v
	chi1inplaney = chi1y_v - chi1dL_v*Lny_v
	chi1inplanez = chi1z_v - chi1dL_v*Lnz_v
	chi2inplanex = chi2x_v - chi2dL_v*Lnx_v
	chi2inplaney = chi2y_v - chi2dL_v*Lny_v
	chi2inplanez = chi2z_v - chi2dL_v*Lnz_v

	print("cos tilt1 = %f, cos tilt2 = %f, cos phi12 = %f"%(chi1dL_v/chi1_norm, chi2dL_v/chi2_norm, (chi1inplanex*chi2inplanex + chi1inplaney*chi2inplaney + chi1inplanez*chi2inplanez)/((chi1inplanex*chi1inplanex + chi1inplaney*chi1inplaney + chi1inplanez*chi1inplanez)*(chi2inplanex*chi2inplanex + chi2inplaney*chi2inplaney + chi2inplanez*chi2inplanez))**0.5))

	return np.arccos(chi1dL_v/chi1_norm), np.arccos(chi2dL_v/chi2_norm), np.arccos((chi1inplanex*chi2inplanex + chi1inplaney*chi2inplaney + chi1inplanez*chi2inplanez)/((chi1inplanex*chi1inplanex + chi1inplaney*chi1inplaney + chi1inplanez*chi1inplanez)*(chi2inplanex*chi2inplanex + chi2inplaney*chi2inplaney + chi2inplanez*chi2inplanez))**0.5)

def find_S_and_L_at_freq_dt(v0, m1, m2, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, Lnx, Lny, Lnz, v_final, dt):
        """ given the spins and ang momentum at a given frequency, find the spins and L at a later frequency """

        # evolve the spins                                                                                                                                                           
        v_v, chi1x_v, chi1y_v, chi1z_v, chi2x_v, chi2y_v, chi2z_v, Lnx_v, Lny_v, Lnz_v = evolve_spins_dt(v0, m1, m2, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, Lnx, Lny, Lnz, v_final, dt)

        # dot products                                                                                                                                                               
        chi1dL_v = chi1x_v*Lnx_v + chi1y_v*Lny_v + chi1z_v*Lnz_v
        chi2dL_v = chi2x_v*Lnx_v + chi2y_v*Lny_v + chi2z_v*Lnz_v

        # norms                                                                                                                                                                      
        chi1_norm = np.sqrt(chi1x_v**2+chi1y_v**2+chi1z_v**2)
        chi2_norm = np.sqrt(chi2x_v**2+chi2y_v**2+chi2z_v**2)
        L_norm = np.sqrt(Lnx_v**2+Lny_v**2+Lnz_v**2)

        return chi1x_v[-1], chi1y_v[-1], chi1z_v[-1], chi2x_v[-1], chi2y_v[-1], chi2z_v[-1], Lnx_v[-1], Lny_v[-1], Lnz_v[-1], v_v[-1]
