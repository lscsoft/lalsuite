"""
This script computes the tilt angles at infinity from a given reference frequency, by combining orbit-averaged evolution at higher frequencies
until a transition frequency, with precession-averaged evolution until infinite separation.
There is also an option to compute the bounds on the tilts and an average value at a finite separation, though this has not been tested extensively.
This implementation is described in the paper, <https://dcc.ligo.org/P2100029>, arXiv:2107.11902

Sumeet Kulkarni, 2021
"""

import numpy as np
import lal
from .calc_tilts_prec_avg_regularized import prec_avg_tilt_comp
from .tilts_at_infinity_utils import *
import lalsimulation as lalsim
from warnings import warn

# Define a function to transform spin basis:
def get_tilts(chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, Lnx, Lny, Lnz):

        """
        Given the spins and ang momentum at a given frequency, find the tilt and in-plane spin angles at that frequency.

        Inputs:
        chi1x, chi1y, chi1z: Cartesian spin-magnitude components for the primary object (m1) of the binary
        chi2x, chi2y, chi2z: Cartesian spin-magnitude components for the secondary object (m2) of the binary
        Lnx, Lny, Lnz: Cartesian components of the direction of the Newtonian orbital angular momentum (will be normalized) of the binary

        Output:
        tilt1, tilt2: tilt angles of the binary spin-vectors w.r.t. the z-axis: the direction of the Newtonian orbital angular momentum
        phi12: angle between projection of the two spin-vectors onto the xy plane
        """

        chi1 = np.array([chi1x, chi1y, chi1z])
        chi2 = np.array([chi2x, chi2y, chi2z])

        Ln = np.array([Lnx, Lny, Lnz])

        # norms and normalizing
        chi1_norm = np.linalg.norm(chi1)
        chi2_norm = np.linalg.norm(chi2)

        Ln /= np.linalg.norm(Ln)

        # dot products
        chi1dL = np.dot(chi1, Ln)
        chi2dL = np.dot(chi2, Ln)

        # in-plane spins
        chi1inplane = chi1 - chi1dL*Ln
        chi2inplane = chi2 - chi2dL*Ln

        # Defining cosine of tilts and phi12
        cos_tilt1 = np.clip(chi1dL/chi1_norm, -1., 1.)
        cos_tilt2 = np.clip(chi2dL/chi2_norm, -1., 1.)
        cos_phi12 = np.clip(np.dot(chi1inplane,chi2inplane)/(np.linalg.norm(chi1inplane) * np.linalg.norm(chi2inplane)),-1.,1.)

        # set quadrant of phi12
        phi12_evol_i = np.arccos(cos_phi12)
        if np.sign(np.dot(Ln,np.cross(chi1, chi2))) < 0:
            phi12_evol_i = 2.*np.pi - phi12_evol_i

        return np.arccos(cos_tilt1), np.arccos(cos_tilt2), phi12_evol_i


# Define a function to calculate v_trans based on the fitting curve, for a given mass ratio
def calc_v_trans(q):
    """
    Calculates the transition orbital speed (v_trans) to shift from orbit-averaged to precession-averaged evolution in this
    hybrid spin evolution code. v_trans depends on the mass ratio, (q), and is determined using the fitting curve from the
    paper [Eq. (31) in the paper, <https://dcc.ligo.org/P2100029>, arXiv:2107.11902]

    Input:

    q, the binary mass ratio defined as m2/m1, with m1 being the primary (heavier) object; range (0,1].

    Output:
    v_transition (float)
    """
    if q <= 0. or q > 1.:
        raise ValueError("The mass ratio must be a float between 0 and 1, defined as m2/m1 where m1 is the heavier component")
    # Coefficients of the fit:
    a = -0.05
    b = 0.06
    return(a*q**2 + b)


# The number of steps in integration, based on the v_trans:
def get_nsteps(v_trans):
    """
    Determine the number of frequency steps to use in each segment of the integration while performing orbit-averaged evolution
    from fref to v_trans

    Input:
    v_trans (float): the transition orbital speed determined by the fit in calc_v_trans()
    The value for v_trans must be >= 0.01, which is the least transition orbital speed given in our fit in the paper [Eq. (31) in
    <https://dcc.ligo.org/P2100029>, arXiv:2107.11902]

    Output:
    n_steps: (int), the number of steps for evolution to be used in calc_tilts_at_infty_hybrid_evolve()
    """
    if v_trans >= 0.05:
        n = 200
    elif 0.03 < v_trans <= 0.05:
        n = 800
    elif 0.02 < v_trans <= 0.03:
        n = 3200
    elif 0.01 <= v_trans <= 0.02:
        n = 12800
    else:
        raise ValueError("The number of steps has not been calibrated for v_trans < 0.01 and might lead to memory errors")
    return n


# The timesteps based on how low to go in v:
def get_dt_constant(v_trans):
    """
    Determine the step-size (dt) to use in each frequency interval used for orbit-averaged evolution from fref to v_trans.
    This function returns the denominator X in the constant c = 1/X, such that f_high*dt = c in determining the timestep
    of orbital evolution. Here, f_high is the higher frequency in the frequency interval.

    Input:
    v_trans (float): the transition orbital speed determined by the fit in calc_v_trans()

    Output:
    dt_constant: (int) to be used in calc_tilts_at_infty_hybrid_evolve()
    """
    if v_trans > 0.05:
        dt_c = 256
    elif 0.03 < v_trans <= 0.05:
        dt_c = 128
    elif 0.01 <= v_trans <= 0.03:
        dt_c = 64
    else:
        raise ValueError("Evolution to v_trans lower than 0.01 not possible with current configuration")
    return dt_c


def calc_tilts_at_infty_hybrid_evolve(m1, m2, chi1, chi2, tilt1, tilt2, phi12, fref, approx="SpinTaylorT5", spinO=6, lscorr=1, verbose=False, prec_only=False, version='v1', failure_mode='None', **kwargs):
    """
    Calculate tilts at infinity with hybrid orbit-averaged and precession-averaged evolution
    Evolves tilt1 and tilt2 for a given binary from a reference frequency, fref, to infinite separation by first evolving using orbit-averaged evolution
    until a transition orbital speed (v_trans) given by calc_v_trans(), and from that point until infinity using precession-averaged evolution.

    There is also an option to compute the bounds on the tilts and average values at a finite separation (given by the value of the orbital angular momentum, Lf)
    but the application of the hybrid evolution to this case is not so well tested. In particular, for small enough values of Lf, the hybrid evolution will end
    at a larger value of the orbital angular momentum than Lf and the precession-averaged evolution will fail because of this. This function does not check for
    this case.

    Inputs:

    --Required--
    m1, m2: Detector frame masses of the binary, in kg
    chi1, chi2: The dimensionless spin-magnitudes of the binary
    tilt1, tilt2: tilt angles of the binary spin-vectors w.r.t. the z-axis: the direction of the Newtonian orbital angular momentum
    phi12: angle between projection of the two spin-vectors onto the xy plane
    fref: Reference frequency, in Hz; must be <= 20 Hz, since the stepwise evolution is currently only calibrated for such values

    --Optional--
    approx: The approximant to use (options: "SpinTaylorT1", "SpinTaylorT4", "SpinTaylorT5"), default: "SpinTaylorT5"
    spinO: The order of the spin contributions included in the post-Newtonian expansion (possible values:[4,5,6]), default: 6
    lscorr: activates spin contributions to the orbital angular momentum in the precession equations, default: 1 (active)
    verbose (bool): display details of integration steps and print output of orbit-averaged evolution, default: False
    prec_only: Use only the precession-averaged evolution for a quick but less accurate result, default: False
    version: Version of the calculation to use--currently, two versions are available: v1 divides the orbit-averaged part of the hybrid evolution into a series of multiple integration steps, while in v2 it proceeds in a single step using a modified integrator that only outputs the final spin values, default: "v1", though this will be changed to "v2" in a near future release, since the "v2" evolution is much faster and somewhat more accurate
    failure_mode:  How the code behaves when the evolution fails. 'Error' means that the code raises a RuntimeError, while 'NAN' and 'None' mean that the code raises a RuntimeWarning and returns np.nan or None for the output, respectively, default: 'None'

    **kwargs: dict, optional: precession-averaged evolution settings, passed through **kwargs to prec_avg_tilt_comp()
    Please refer to the prec_avg_tilt_comp() documentation in calc_tilts_prec_avg_regularized.py for the list of settings

    NOTE: The keyword argument Lf: Final magnitude of orbital angular momentum, if set to None, gives the output at infinity; this is the default.
          To obtain tilts at a finite separation, provide the corresponding Lf in total mass = 1 units

    NOTE: prec_avg_tilt_comp() options LPNorder and LPNspins are not available through here, but are set automatically by the choice of spinO and lscorr.

    Output:

    Dictionary with entries 'tilt1_inf', 'tilt2_inf' for evolution to infinity, and entries 'tilt1_transition', 'tilt2_transition',
    'phi12_transition' for the tilts at the transition orbital speed (v_trans).

    NOTE: If kwarg Lf is not set to None in prec_avg_tilt_comp(), the output gives bounds and average values for the tilts at a final separation determined by Lf.
          In this case, the entries 'tilt1_inf' and 'tilt2_inf' will be replaced by 'tilt1_sep_min', 'tilt1_sep_max', 'tilt1_sep_avg', 'tilt2_sep_min', 'tilt2_sep_max',
          and 'tilt2_sep_avg'.
    """

    # Check version
    if version in ['v1','v2']:
        version_prec_avg = 'v1'

        if version == 'v1':
            warn("v1 of the hybrid tilts at infinity calculation is deprecated, since v2 is much faster and somewhat more accurate. In a near future release, the default version will be changed to v2.", FutureWarning)
    else:
        raise ValueError("Only versions ['v1', 'v2'] are available at the moment, while version = %s"%version)

    # Set the failure output and the string corresponding to it
    # These default to None, since they are not used if failure_mode == 'Error'

    if failure_mode == 'NAN':
        failure_output = np.nan
        failure_output_string = 'np.nan'
    else:
        failure_output = None
        failure_output_string = 'None'

    # Check if Lf is defined in kwargs:
    if 'Lf' in kwargs:
        Lf = kwargs['Lf']
    else:
        Lf = None

    ## Run checks for input parameter values of the  masses, spin magnitudes, tilts and reference frequency:

    # Masses
    check_masses(m1, m2)
    # Exactly equal masses:
    eq_mass_check(m1, m2, Lf)
    # Spin Magnitudes:
    check_spin_mags(chi1, chi2)
    # Tilt inclinations:
    check_tilts(tilt1, tilt2)
    # Reference Frequency:
    check_fref(fref, m1, m2, evol_type="hybrid")

    # Return tilts as they are for non-spinning and exactly aligned/anti-aligned cases:
    if (tilt1 == np.pi and tilt2 in [0., np.pi]) or (tilt1 == 0. and tilt2 == 0.):
        return package_tilts(tilt1, tilt2, Lf, swap=False)

    # Note: returning the input tilts for single non-spinning cases not appropriate if the instantaneous terms in the orbit-averaged evolution are activated.
    if chi1 == 0. or chi2 == 0.:
        return package_tilts(tilt1, tilt2, Lf, swap=False)

    # If prec_only is set to True, omit hybrid evolution and return tilts at infinity using only precession-averaged evolution for fast results.
    if prec_only:
        spin_angles_output  = prec_avg_tilt_comp(m1, m2, chi1, chi2, tilt1, tilt2, phi12, fref, LPNorder = 2, LPNspins = False, **kwargs)
        spin_angles_output["tilt1_transition"] = None
        spin_angles_output["tilt2_transition"] = None
        spin_angles_output["phi12_transition"] = None
        spin_angles_output["f_transition"] = None
        return spin_angles_output


    M = m1 + m2 #total mass, in kg.
    q = min(m1/m2,m2/m1) #mass ratio (make sure this is in the range [0,1) )

    # Save input masses in kg.
    m1_kg = m1
    m2_kg = m2

    # Rescale masses and frequency to the stellar-mass scale of total mass = 200 M_sun units:
    fac = 200*MSUN_SI/M
    m1 *= fac
    m2 *= fac
    fref /= fac

    # Define the rescaled total mass in seconds, for usage in converting from frequencies to orbital speeds, and vice versa:
    MT_s = (m1 + m2) * kg_to_s

    # Transition orbital speed:
    v_trans = calc_v_trans(q)
    if verbose:
        print("v_trans = {}".format(v_trans))
    f_trans = v_trans**3/np.pi/MT_s

    # Set parameters:
    v_ref = np.power((fref*(np.pi*MT_s)),1./3)

    phaseO = 7
    tideO = 0


    if lscorr == 0:
        LPNspins = False
    elif lscorr == 1:
        LPNspins = True
    else:
        raise ValueError("The value of lscorr should either be 0 or 1")
    inst = 0
    approx = lalsim.GetApproximantFromString(approx)

    if spinO == 4:
        LPNorder = 0
    elif spinO == 5 and lscorr == 0:
        LPNorder = 1
    elif spinO == 5 and lscorr == 1:
        LPNorder = 1.5
    elif spinO == 6 and lscorr == 0:
        LPNorder = 1
    elif spinO == 6 and lscorr == 1:
        LPNorder = 1.5
    # spinO == 7 will only be available when the instantaneous evolution is implemented
    #elif spinO == 7 and lscorr == 0:
    #    LPNorder = 2
    #elif spinO == 7 and lscorr == 1:
    #    LPNorder = 2.5
    else:
        raise ValueError("Check your SpinO input, which must be one from [4,5,6]; given input = {}".format(spinO))

   # Transform from spherical to cartesian for initial values of the tilts:
    chi1x_high, chi1y_high, chi1z_high = chi1*np.sin(tilt1), 0.0, chi1*np.cos(tilt1)
    chi2x_high, chi2y_high, chi2z_high = chi2*np.sin(tilt2)*np.cos(phi12), chi2*np.sin(tilt2)*np.sin(phi12), chi2*np.cos(tilt2)

    Lnx_high = 0.0
    Lny_high = 0.0
    Lnz_high = 1.0

    E1x_high = 1.0
    E1y_high = 0.0
    E1z_high = 0.0

    lalpars=lal.CreateDict()
    lalsim.SimInspiralWaveformParamsInsertFinalFreq(lalpars, f_trans)
    lalsim.SimInspiralWaveformParamsInsertPNSpinOrder(lalpars, spinO)
    lalsim.SimInspiralWaveformParamsInsertPNPhaseOrder(lalpars, phaseO)
    lalsim.SimInspiralWaveformParamsInsertLscorr(lalpars, lscorr)

    if f_trans < fref:

        ###################### v2: single-step evolution ############################################

        if version == "v2":

            # Use a uniform timestep (with dt_constant = 64) for all binary parameters:
            dt = 1./(64 * fref)
            # Activate option to save only final values in Spin Taylor evolution:
            lalsim.SimInspiralWaveformParamsInsertOnlyFinal(lalpars, 1)

            dictparams = {'phiRef':0., 'deltaT': dt, 'm1_SI': m1, 'm2_SI': m2, 'fStart': fref, 'fRef': fref, 's1x': chi1x_high, 's1y': chi1y_high, 's1z': chi1z_high,
                          's2x': chi2x_high, 's2y': chi2y_high, 's2z': chi2z_high, 'lnhatx': Lnx_high, 'lnhaty': Lny_high, 'lnhatz': Lnz_high,
                          'e1x': E1x_high, 'e1y': E1y_high, 'e1z': E1z_high, 'LALparams': lalpars, 'approx': approx}

            try:
                _,_,c1x, c1y, c1z, c2x, c2y, c2z, Lnx, Lny, Lnz, E1x, E1y, E1z = lalsim.SimInspiralSpinTaylorOrbitalDriver(**dictparams)
            except:
                failure_message = "The orbit-averaged evolution failed."
                return evolution_error_handling(failure_mode, failure_message, failure_output, failure_output_string, Lf, swap=False, hybrid_evol=True)

            c1x = c1x.data.data[-1]
            c2x = c2x.data.data[-1]
            c1y = c1y.data.data[-1]
            c2y = c2y.data.data[-1]
            c2z = c2z.data.data[-1]
            c1z = c1z.data.data[-1]
            Lnx = Lnx.data.data[-1]
            Lny = Lny.data.data[-1]
            Lnz = Lnz.data.data[-1]

            tilt1_transition, tilt2_transition, phi12_transition = get_tilts(c1x, c1y, c1z, c2x, c2y, c2z, Lnx, Lny, Lnz)

            if verbose:
                print("The tilts at transition are: tilt1 = {0}, tilt2 = {1}, phi12 = {2}".format(tilt1_transition, tilt2_transition, phi12_transition))

            spin_angles_output = prec_avg_tilt_comp(m1, m2, chi1, chi2, tilt1_transition, tilt2_transition, phi12_transition, f_trans, LPNorder = LPNorder, LPNspins = LPNspins,
                                                        version=version_prec_avg, **kwargs)
            spin_angles_output["tilt1_transition"] = tilt1_transition
            spin_angles_output["tilt2_transition"] = tilt2_transition
            spin_angles_output["phi12_transition"] = phi12_transition
            spin_angles_output["f_transition"] = v_trans**3/np.pi/(M*kg_to_s) #return the transition frequency based on the original input masses in kg.


    ################### v1: multi-step evolution #############################

        elif version == 'v1':

            # Use a timestep optimized by the get_dt_constant() function:
            dt_constant = get_dt_constant(v_trans)
            dt = 1./(dt_constant * fref)

            dictparams = {'phiRef':0., 'deltaT': dt, 'm1_SI': m1, 'm2_SI': m2, 'fStart': fref, 'fRef': fref, 's1x': chi1x_high, 's1y': chi1y_high, 's1z': chi1z_high,
                          's2x': chi2x_high, 's2y': chi2y_high, 's2z': chi2z_high, 'lnhatx': Lnx_high, 'lnhaty': Lny_high, 'lnhatz': Lnz_high,
                          'e1x': E1x_high, 'e1y': E1y_high, 'e1z': E1z_high, 'LALparams': lalpars, 'approx': approx}


    ####################### Step 1: Use LALSIM.SimInspiralSpinTaylorPNEvolveOrbit to get orbit-averaged spins at the transition frequency #########################

            n_steps = get_nsteps(v_trans)

            if verbose:
                print("Starting orbital evolution in nsteps = ",n_steps)

            v_start = 0.5*(v_ref + v_trans) # This is the v-value until which all binaries can be evolved from a given fref without segfaulting # NEEDS CALIBRATION FOR fref > 20 Hz
            v_arr = np.geomspace(v_start,v_trans,n_steps+1)
            v_arr = np.insert(v_arr, 0, v_ref) # start the orbital speed integration steps array with v_ref

            for i in range(len(v_arr)-1):
                    if verbose:
                        print("starting step {}".format(i+1))
                    v_high = v_arr[i]
                    f_high = v_high**3/np.pi/MT_s
                    dictparams['fStart'] = f_high
                    dictparams['fRef'] = f_high

                    v_low = v_arr[i+1]
                    f_low = v_low**3/np.pi/MT_s
                    lalsim.SimInspiralWaveformParamsInsertFinalFreq(lalpars, f_low)
                    dictparams['LALparams'] = lalpars
                    dictparams['deltaT'] = 1./(dt_constant * f_high)
                    if verbose:
                        print("Evolving to v_low = {0} corresponding to f_low = {1} with dt = {2}".format(v_low, f_low, dictparams['deltaT']))

                    try:
                        v_out,_,chi1x_low, chi1y_low, chi1z_low, chi2x_low, chi2y_low, chi2z_low, Lnx_low, Lny_low, Lnz_low, E1x, E1y, E1z = lalsim.SimInspiralSpinTaylorOrbitalDriver(**dictparams)
                    except:
                        failure_message = "The orbit-averaged evolution failed."
                        return evolution_error_handling(failure_mode, failure_message, failure_output, failure_output_string, Lf, swap=False, hybrid_evol=True)

                    if abs(v_out.data.data[0] - v_low) > 1.e-4:
                        failure_message = f"The orbit-averaged evolution failed to reached the specified final velocity, v={v_low} at step {i+1}. The input binary parameters cannot be evolved from the given frequency {fref*fac} Hz."
                        return evolution_error_handling(failure_mode, failure_message, failure_output, failure_output_string, Lf, swap=False, hybrid_evol=True)


                    dictparams['s1x'] = chi1x_low.data.data[0]
                    dictparams['s1y'] = chi1y_low.data.data[0]
                    dictparams['s1z'] = chi1z_low.data.data[0]

                    dictparams['s2x'] = chi2x_low.data.data[0]
                    dictparams['s2y'] = chi2y_low.data.data[0]
                    dictparams['s2z'] = chi2z_low.data.data[0]

                    dictparams['lnhatx'] = Lnx_low.data.data[0]
                    dictparams['lnhaty'] = Lny_low.data.data[0]
                    dictparams['lnhatz'] = Lnz_low.data.data[0]

                    dictparams['e1x'] = E1x.data.data[0]
                    dictparams['e1y'] = E1y.data.data[0]
                    dictparams['e1z'] = E1z.data.data[0]

            c1x = dictparams['s1x']
            c1y = dictparams['s1y']
            c1z = dictparams['s1z']

            c2x = dictparams['s2x']
            c2y = dictparams['s2y']
            c2z = dictparams['s2z']

            Lnx = dictparams['lnhatx']
            Lny = dictparams['lnhaty']
            Lnz = dictparams['lnhatz']

            tilt1_transition, tilt2_transition, phi12_transition = get_tilts(c1x, c1y, c1z, c2x, c2y, c2z, Lnx, Lny, Lnz)

            if verbose:
                print("The tilts at transition are: tilt1 = {0}, tilt2 = {1}, phi12 = {2}".format(tilt1_transition, tilt2_transition, phi12_transition))

         ################################# Step 2: Use precession-averaged evolution to compute tilts at infinity ###################################################
            spin_angles_output = prec_avg_tilt_comp(m1, m2, chi1, chi2, tilt1_transition, tilt2_transition, phi12_transition, f_trans, LPNorder = LPNorder, LPNspins = LPNspins,
                                                    version=version_prec_avg, **kwargs)
            spin_angles_output["tilt1_transition"] = tilt1_transition
            spin_angles_output["tilt2_transition"] = tilt2_transition
            spin_angles_output["phi12_transition"] = phi12_transition
            spin_angles_output["f_transition"] = v_trans**3/np.pi/(M*kg_to_s) #return the transition frequency based on the original input masses in kg.

    else:
            if verbose:
                print("The fitted transition frequency is higher than fref; Computing tilts at infinity from {} Hz. instead using precession-averaged evolution only".format(fref))
            spin_angles_output  = prec_avg_tilt_comp(m1, m2, chi1, chi2, tilt1, tilt2, phi12, fref, LPNorder = LPNorder, LPNspins = LPNspins, version=version_prec_avg,
                                                     **kwargs)
            spin_angles_output["tilt1_transition"] = None
            spin_angles_output["tilt2_transition"] = None
            spin_angles_output["phi12_transition"] = None
            spin_angles_output["f_transition"] = None

    return spin_angles_output
