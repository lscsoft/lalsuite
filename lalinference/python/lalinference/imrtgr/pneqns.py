"""
Code to use post-Newtonian equations [from Ajith, PRD 84, 084037 (2011), arXiv:1107.1267] to evolve spins in a binary black hole

P. Ajith, A. Gupta, and N. K. Johnson-McDaniel, 04.2016, based on earlier code
"""

from scipy.integrate import ode
import numpy as np
import lalsimulation as lalsim
from lal import PI, MTSUN_SI, MSUN_SI, GAMMA
from numpy import log
from numpy.linalg import norm


"""Preccesion frequency spins """
def Omega(v, m1, m2, S1, S2, Ln): #S1, S2 are full spin vector while Ln is a unit vector

    ''' Eqs.(3.8) of Ajith (2011) (http://arxiv.org/pdf/1107.1267v2.pdf)'''

    m = m1+m2
    deltam = m1-m2
    eta = (m1*m2)/(m1+m2)**2.
    return (v**5./m)*((3./4+eta/2.-3.*deltam/m/4.)*Ln+0.5*v/m**2.*(-3.0*np.dot((S2 +(m2/m1)*S1), Ln)*Ln+S2)+v**2.*(9./16.+5.*eta/4.-eta**2./24.-9.*deltam/m/16.+5.*deltam*eta/8./m)*Ln)


""" compute the re-expanded dEnergy/flux """
def denergy_by_flux(v, eta, delta, chiadL, chisdL, chiasqr, chissqr, chisdchia, order):

    ''' Eqs.(3.2) of Ajith (2011) http://arxiv.org/pdf/1107.1267v2.pdf '''

    # different powers of v
    v2 = v*v; v3 = v2*v; v4 = v3*v; v5 = v4*v; v6 = v5*v; v7 = v6*v; v9 = v7*v2

    # initialize the cofficients
    dEbF0 = dEbF2 = dEbF3 = dEbF4 = dEbF5 = dEbF6 = dEbF6L = dEbF7 = 0.

    if order >= 0:
            dEbF0 = -5./(32.*eta)
    if order >= 2:
            dEbF2 = 2.2113095238095237 + (11*eta)/4.
    if order >= 3:
            dEbF3 = (113*chiadL*delta)/12. + chisdL*(9.416666666666666 - (19.*eta)/3.) - 4.*PI
    if order >= 4:
            dEbF4 = 3.010315295099521 + (233*chisdchia*delta)/48. - (719.*chiadL*chisdL*delta)/48. + \
                    chiasqr*(2.4270833333333335 - 10.*eta) + chisdL*chisdL*(-7.489583333333333 - eta/24.) + \
                    chissqr*(2.4270833333333335 + (7.*eta)/24.) + (5429.*eta)/1008. + (617*eta*eta)/144. + \
                    chiadL*chiadL*(-7.489583333333333 + 30.*eta)
    if order >= 5:
                    dEbF5 = chiadL*delta*(72.71676587301587 + (7*eta)/2.) + chisdL*(72.71676587301587 - \
                                            (1213*eta)/18. - 17*eta*eta/2.) - (7729*PI)/672. + (13*eta*PI)/8.
    if order >= 6:
            dEbF6 = -115.2253249962622 - 15211*eta*eta/6912. + 25565*eta*eta*eta/5184. + \
                                    32*PI*PI/3. + eta*(258.1491854023631 - 451*PI*PI/48.) + (1712*GAMMA)/105.
            dEbF6L = 1712./105.
    if order >= 7:
            dEbF7 = (-15419335.*PI)/1.016064e6 - (75703.*eta*PI)/6048. + (14809.*eta*eta*PI)/3024.

    return (dEbF0/v9)*(1. + dEbF2*v2 + dEbF3*v3 + dEbF4*v4 + dEbF5*v5 + (dEbF6+dEbF6L*log(4.*v))*v6 + dEbF7*v7)


def precession_eqns(t, y_vec, m1, m2):
	""" The coupled set of ODEs containing the PN precession eqns as well as the evolution of dv/dt
	    All these equations are listed in Ajith (2011) (http://arxiv.org/pdf/1107.1267v2.pdf).
	"""

	# unpack the vector of variables
	v, S1x, S1y, S1z, S2x, S2y, S2z, Lx, Ly, Lz = y_vec

	m1_sqr = m1*m1
	m2_sqr = m2*m2
	m = m1+m2
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

	# magnitude of the orb ang momentum. The prefactor of dLN/dt in Eq.(3.6) of Ajith (2011)
	Lmag = (eta*m**2./v)*(1.+(1.5+eta/6.)*v**2.)

	# precession freqs. Eqs.(3.8) of Ajith (2011)
        Omega1= Omega(v, m1, m2, S1, S2, Ln)
        Omega2= Omega(v, m2, m1, S2, S1, Ln)

	# spin precession eqns. Eqs.(3.7) of Ajith (2011)
	dS1_dt = np.cross(Omega1, S1)
	dS2_dt = np.cross(Omega2, S2)

	# ang momentum precession eqn. Eqs.(3.6) of Ajith (2011)
	dL_dt =-(dS1_dt+dS2_dt)/Lmag

	# orb frequency evolution. Eqs.(3.5) of Ajith (2011)
	dv_dt = -1./(m*denergy_by_flux(v, eta, delta, chiadL, chisdL, chiasqr, chissqr, chisdchia, 7))

	return  [dv_dt, dS1_dt[0], dS1_dt[1], dS1_dt[2], dS2_dt[0], dS2_dt[1], dS2_dt[2], dL_dt[0], dL_dt[1], dL_dt[2]]


def evolve_spins_dt(v0, m1, m2, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, Lnx, Lny, Lnz, v_final, dt):
	""" evolve the spins and orb ang momentum according to the PN precession eqns."""

        # some sanity checks
        if m1<0:
            raise ValueError("ERROR: mass1 is negative")
        if m2<0:
            raise ValueError("ERROR: mass2 is negative")
        if (chi1x**2+chi1y**2+chi1z**2)>1.:
            raise ValueError("ERROR: magnitude of spin1 is greater than 1")
        if (chi2x**2+chi2y**2+chi2z**2)>1.:
            raise ValueError("ERROR: magnitude of spin2 is greater than 1")
        if dt<0:
            raise ValueError("ERROR: time step is negative")

	m1_sqr = m1*m1
	m2_sqr = m2*m2

	# pack the inputs
	y_vec0 = [v0, chi1x*m1_sqr, chi1y*m1_sqr, chi1z*m1_sqr, chi2x*m2_sqr, chi2y*m2_sqr, chi2z*m2_sqr, Lnx, Lny, Lnz]
	t0 = 0.

        # Defining chi_eff and Newtonian chirp time. T_MAX is set to twice of chirp time.
        chi_eff = lalsim.SimInspiralTaylorF2ReducedSpinComputeChi(m1, m2, chi1z, chi2z)
	T_MAX = 2.* lalsim.SimInspiralTaylorF2ReducedSpinChirpTime(v0**3./(PI*(m1+m2)*MTSUN_SI), m1*MSUN_SI, m2*MSUN_SI, chi_eff, 0)/MTSUN_SI
	R_TOL = 1e-6
	A_TOL = 1e-6

        #print("T_MAX = %f"%(T_MAX))

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
	while solver.successful() and solver.t < 2.*T_MAX and solver.y[0] <= v_final and abs(np.sqrt((solver.y[7])**2+(solver.y[8])**2+(solver.y[9])**2)-1.)<1e-8:

		solver.integrate(solver.t + dt, step=1)
		y_result.append(solver.y)
		t_output.append(solver.t)

		#print '... t = ', solver.t, ' y = ', solver.y

	Y = np.array(y_result)
	t = np.array(t_output)
	v_v, S1x_v, S1y_v, S1z_v, S2x_v, S2y_v, S2z_v, Lx_v, Ly_v, Lz_v  = Y.T

        if abs(v_v[-1]-v_final)>1e-2:
            raise ValueError("v_max = {0} is more than 0.01 different from v_final = {1}".format(v_v[-1], v_final))

        if abs(np.sqrt((solver.y[7])**2+(solver.y[8])**2+(solver.y[9])**2)-1.)>=1e-8:
            raise ValueError("norm of Ln is more than 10^{-8} different from 1, with distance {0}".format(abs(np.sqrt((solver.y[7])**2+(solver.y[8])**2+(solver.y[9])**2)-1.)))

	return v_v, S1x_v/m1_sqr, S1y_v/m1_sqr, S1z_v/m1_sqr, S2x_v/m2_sqr, S2y_v/m2_sqr, S2z_v/m2_sqr, Lx_v, Ly_v, Lz_v


def find_tilts_and_phi12_at_freq(v0, m1, m2, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, Lnx, Lny, Lnz, v_final, dt):
        """ given the spins and ang momentum at a given frequency, find the tilt and in-plane spin angles at a later frequency """

	print("v0 = %f, m1 = %f, m2 = %f, chi1x = %f, chi1y = %f, chi1z = %f, chi2x = %f, chi2y = %f, chi2z = %f"%(v0, m1, m2, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z))

        # evolve the spins
        v_v, chi1x_v, chi1y_v, chi1z_v, chi2x_v, chi2y_v, chi2z_v, Lnx_v, Lny_v, Lnz_v = evolve_spins_dt(v0, m1, m2, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, Lnx, Lny, Lnz, v_final, dt)

	chi1_v = np.array([chi1x_v[-1], chi1y_v[-1], chi1z_v[-1]])
	chi2_v = np.array([chi2x_v[-1], chi2y_v[-1], chi2z_v[-1]])

	Ln_v = np.array([Lnx_v[-1], Lny_v[-1], Lnz_v[-1]])

        # norms and normalizing
        chi1_norm = norm(chi1_v)
        chi2_norm = norm(chi2_v)

        Ln_v /= norm(Ln_v)

        # dot products
        chi1dL_v = np.dot(chi1_v, Ln_v)
        chi2dL_v = np.dot(chi2_v, Ln_v)

	# in-plane spins

	chi1inplane = chi1_v - chi1dL_v*Ln_v
	chi2inplane = chi2_v - chi2dL_v*Ln_v

        # Defining cosine of tilts and phi12
        cos_tilt1 = chi1dL_v/chi1_norm
        cos_tilt2 = chi2dL_v/chi2_norm
        cos_phi12 = np.dot(chi1inplane,chi2inplane)/(norm(chi1inplane)*norm(chi2inplane))

	print("cos tilt1 = %f, cos tilt2 = %f, cos phi12 = %f"%(cos_tilt1, cos_tilt2, cos_phi12))

	return np.arccos(cos_tilt1), np.arccos(cos_tilt2), np.arccos(cos_phi12)
