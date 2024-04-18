"""
Class for NRSur7dq4Remnant model for the mass, spin and kick of the final BH of
generically precessing binary black holes.

Vijay Varma, 2019.
"""

import numpy as np
import lalsimulation as lalsim
import lal

from . import pn_spin_evolution_wrapper
from . import quaternion_utils
from .nrfits import NRFits

class NRSur7dq4Remnant(NRFits):
    r"""
    Class for NRSur7dq4Remnant model for the mass, spin and kick of the final BH
    of generically precessing binary black holes.

    Paper: arxiv:1905.09300

    Parameter ranges for usage: \n
        q = [1, 6.01] \n
        \f$\chi_{A}, \chi_{B}\f$ = [-1, 1]

    Training parameter ranges: \n
        q = [1, 4.01] \n
        \f$\chi_{A}, \chi_{B}\f$ = [-0.81, 0.81]

    But extrapolates reasonably to the above mass ratios and spins. However,
    if a guarantee of accuracy is required, this model should be used within
    the training parameter range.
    """


    #-------------------------------------------------------------------------
    def _get_coorbital_frame_spins_at_idx(self, chiA, chiB, omega, lNhat, phi, \
            idx):
        """ Computes PN spins and dynamics at a given idx.

        Inputs:
            chiA:       Dimless spin evolution of BhA in inertial frame.
            chiB:       Dimless spin evolution of BhB in inertial frame.
            omega:      Orbital frequency evolution (this is total-mass times
                        the angular orbital frequency).
            lNhat:      Orbital angular momentum direction evolution.
            phi:        Orbital phase evolution.
            idx:        Index for output.

        Outputs (all are time series):
            chiA_at_idx_coorb:  Spin of BhA at idx, in coorbital frame.
            chiB_at_idx_coorb:  Spin of BhB at idx, in coorbital frame.
            quat_copr_at_idx:   Coprecessing frame quaternion at idx.
            phi_at_idx:         Orbital phase in the coprecessing frame at idx.
            omega_at_idx        Orbital frequency at idx (this is total-mass
                                times the angular orbital frequency).

        The inertial frame is assumed to be aligned to the coorbital frame at
        the first index.
        """

        # Get omega, inertial spins, angular momentum direction and orbital
        # phase at idx
        omega_at_idx = omega[idx]
        chiA_at_idx = chiA[idx]
        chiB_at_idx = chiB[idx]
        lNhat_at_idx = lNhat[idx]
        phi_at_idx = phi[idx]

        # Align the z-direction along orbital angular momentum direction
        # at idx. This moves us into the coprecessing frame.
        quat_copr_at_idx = quaternion_utils.align_vec_quat(lNhat_at_idx)
        # The zeroth element is taken because transform_time_dependent_vector
        # is designed for arrays, but here we only have one element
        chiA_at_idx_copr = quaternion_utils.transform_time_dependent_vector(
                np.array([quat_copr_at_idx]).T,
                np.array([chiA_at_idx]).T, inverse=1).T[0]
        chiB_at_idx_copr = quaternion_utils.transform_time_dependent_vector(
                np.array([quat_copr_at_idx]).T,
                np.array([chiB_at_idx]).T, inverse=1).T[0]

        # get coorbital frame spins at idx
        chiA_at_idx_coorb = quaternion_utils.rotate_in_plane(
                chiA_at_idx_copr, phi_at_idx)
        chiB_at_idx_coorb = quaternion_utils.rotate_in_plane(
                chiB_at_idx_copr, phi_at_idx)

        return chiA_at_idx_coorb, chiB_at_idx_coorb, quat_copr_at_idx, \
            phi_at_idx, omega_at_idx

    #-------------------------------------------------------------------------
    def _get_surrogate_dynamics(self, q, chiA0, chiB0, init_quat, \
            init_orbphase, omega_ref, unlimited_extrapolation):
        """ A wrapper for NRSur7dq4 dynamics.

        Inputs:
            q:          Mass ratio, mA/mB >= 1.
            chiA0:      Dimless spin of BhA in the coorbital frame at omega_ref.
            chiB0:      Dimless spin of BhB in the coorbital frame at omega_ref.
            init_quat:  Coprecessing frame quaternion at omega_ref.
            init_orbphase:
                        Orbital phase in the coprecessing frame at omega_ref.
            omega_ref:  Orbital frequency in the coprecessing frame at the
                        reference epoch where the input spins are given. Note:
                        This is total-mass times the angular orbital frequency.
            unlimited_extrapolation:
                        If True, allows unlimited extrapolation to regions well
                        outside the surrogate's training region. Else, raises
                        an error.

        Outputs:
            t_dyn:      Time values at which the dynamics are returned. These
                        are nonuniform and sparse.
            copr_quat:  Time series of coprecessing frame quaternions.
            orbphase:   Orbital phase time series in the coprecessing frame.
            chiA_copr:  Time series of spin of BhA in the coprecessing frame.
            chiB_copr:  Time series of spin of BhB in the coprecessing frame.
        """

        approxTag = lalsim.SimInspiralGetApproximantFromString("NRSur7dq4")
        LALParams = lal.CreateDict()
        if unlimited_extrapolation:
            lal.DictInsertUINT4Value(LALParams, "unlimited_extrapolation", 1)

        t_dyn, quat0, quat1, quat2, quat3, orbphase, chiAx, chiAy, chiAz, \
            chiBx, chiBy, chiBz = lalsim.PrecessingNRSurDynamics(q, \
            chiA0[0], chiA0[1], chiA0[2], chiB0[0], chiB0[1], chiB0[2], \
            omega_ref, init_quat[0], init_quat[1], init_quat[2], init_quat[3], \
            init_orbphase, LALParams, approxTag)

        t_dyn = t_dyn.data
        orbphase = orbphase.data
        copr_quat = np.array([quat0.data, quat1.data, quat2.data, quat3.data])
        chiA_copr = np.array([chiAx.data, chiAy.data, chiAz.data]).T
        chiB_copr = np.array([chiBx.data, chiBy.data, chiBz.data]).T

        return t_dyn, copr_quat, orbphase, chiA_copr, chiB_copr

    #-------------------------------------------------------------------------
    def _get_pn_spins_at_surrogate_start(self, q, \
            chiA0, chiB0, omega0, omega_switch_IG, t_sur_switch,
            unlimited_extrapolation):
        """ Gets PN spins and frame dynamics at a time close to the start
            of the surrogate waveform model.

        Inputs:
            q:          Mass ratio, mA/mB >= 1.
            chiA0:      Initial dimless spin of BhA in the LAL inertial frame.
            chiB0:      Initial dimless spin of BhB in the LAL inertial frame.
            omega0:     Initial orbital frequency. The spins are specified at
                        this frequency. Note: This is total-mass times the
                        angular orbital frequency.
            omega_switch_IG:
                        Initial guess for the orbital frequency at which to
                        switch from PN to surrogate dynamics. We use PN to
                        evolve the spins until the surrogate becomes valid. The
                        surrogate is in time domain, so it is not
                        straightforward to determine the frequency at which it
                        is valid. So, we evolve with PN until omega_switch_IG.
                        However, this is chosen to be overly conservative such
                        that this initial guess works for all q<=6 cases.
                        Ideally, we would like to use as much as possible of
                        the surrogate spin evolution. So, using the PN spins at
                        omega_switch_IG, the surrogate dynamics are evolved
                        both backwards and forwards in time. This lets us get
                        surrogate orbital frequency at t_sur_switch, as well as
                        the corresponding PN spins at this time.
            t_sur_switch:
                        The time at which the PN spins are returned so that
                        surrogate dynamics can be regenerated starting as early
                        as possible.
            unlimited_extrapolation:
                        If True, Allow unlimited extrapolation to regions well
                        outside the surrogate's training region. Else, raises
                        an error.

        Outputs:
            PN spin and frame dynamics at t_sur_switch, as well as the spins
            and orbital frequency evolution.

            chiA_PN_at_idx_coorb:
                        PN spin of BhA at t_sur_switch, in the coorbital frame.
            chiB_PN_at_idx_coorb:
                        PN spin of BhB at t_sur_switch, in the coorbital frame.
            quat_PN_copr_at_idx:
                        PN coprecessing frame quaternions at t_sur_switch.
            phi_PN_at_idx:
                        PN orbital phase at t_sur_switch.
            omega_PN_at_idx:
                        PN orbital frequency at t_sur_switch.
            chiA_PN:    PN Spin evolution of BhA in the inertial frame.
            chiB_PN:    PN Spin evolution of BhB in the inertial frame.
            omega_PN:   PN orbital frequency evolution. Note: This is total-mass
                        times the angular orbital frequency.
        """

        # Get PN spin evolution starting at omega0
        omega_PN, phi_PN, chiA_PN, chiB_PN, lNhat_PN \
            = pn_spin_evolution_wrapper.spin_evolution(q, chiA0, chiB0, omega0)

        # Get PN coorbital frame spins and frame dynamics at
        # omega_PN=omega_switch_IG
        idx = np.argmin(np.abs(omega_PN - omega_switch_IG))
        chiA_PN_at_idx_coorb, chiB_PN_at_idx_coorb, quat_PN_copr_at_idx, \
            phi_PN_at_idx, omega_PN_at_idx \
            = self._get_coorbital_frame_spins_at_idx(chiA_PN, chiB_PN, \
            omega_PN, lNhat_PN, phi_PN, idx)

        # Now evaluate the surrogate dynamics (both forwards and backwards)
        # using PN spins at omega_switch_IG
        dyn_times, quat_sur, orbphase_sur, chiA_copr_sur, chiB_copr_sur \
            = self._get_surrogate_dynamics(q, chiA_PN_at_idx_coorb, \
            chiB_PN_at_idx_coorb, quat_PN_copr_at_idx, phi_PN_at_idx, \
            omega_switch_IG, unlimited_extrapolation)
        omega_sur = np.gradient(orbphase_sur, dyn_times)

        # Get surrogate orbital frequency at t_sur_switch, which is
        # close to the start of the surrogate data
        omega_init_sur = omega_sur[np.argmin(np.abs(dyn_times - t_sur_switch))]

        # Get PN coorbital frame spins and frame dynamics at omega_init_sur.
        # These must be in the coorbital frame as the surrogate dynamics
        # expects the inputs in the LAL convention which agrees with the
        # coorbital frame.
        idx = np.argmin(np.abs(omega_PN - omega_init_sur))
        chiA_PN_at_idx_coorb, chiB_PN_at_idx_coorb, quat_PN_copr_at_idx, \
            phi_PN_at_idx, omega_PN_at_idx \
            = self._get_coorbital_frame_spins_at_idx(chiA_PN, chiB_PN, \
            omega_PN, lNhat_PN, phi_PN, idx)

        return chiA_PN_at_idx_coorb, chiB_PN_at_idx_coorb, \
            quat_PN_copr_at_idx, phi_PN_at_idx, omega_PN_at_idx, \
            chiA_PN, chiB_PN, omega_PN

    #-------------------------------------------------------------------------
    def _evolve_spins(self, q, chiA0, chiB0, omega0, omega_switch_IG=0.03, \
            t_sur_switch=-4000, return_spin_evolution=False, \
            unlimited_extrapolation=False):
        """ Uses PN and surrogate spin evolution to evolve spins of the
        component BHs from an initial orbital frequency = omega0 until t=-100 M
        from the peak of the waveform.

        Inputs:
            q:          Mass ratio, mA/mB >= 1.
            chiA0:      Initial dimless spin of BhA in the inertial frame.
            chiB0:      Initial dimless spin of BhB in the inertial frame. Note
                        that the LAL inertial frame coincides with the
                        coorbital frame used in the surrogate models.
            omega0:     Initial orbital frequency in the coprecessing frame.
                        The spins are specified at this frequency. Note: This
                        is total-mass times the angular orbital frequency.
            omega_switch_IG:
                        Initial guess for orbital frequency used to determine
                        if PN spin evolution is required, and where to switch
                        from PN spin evolution to the surrogate.
                If omega0 >= omega_switch_IG:
                        PN evolution is not required and the surrogate dynamics
                        are used to evolve the spin until t=-100M. Default
                        value for omega_switch_IG is 0.03 which is expected to
                        be large enough for q<=6.  NOTE: If you get errors
                        about omega0 being too small for the surrogate, try
                        increasing omega_switch_IG.
                If omega0 < omega_switch_IG:
                        Use PN to evolve the spins until the surrogate becomes
                        valid. The surrogate is in time domain, so it is not
                        straightforward to determine the frequency at which it
                        is valid. So, we evolve with PN until omega_switch_IG.
                        However, this is chosen to be overly conservative such
                        that this initial guess works for all q<=6 cases.
                        Ideally, we would like to use as much as possible of
                        the surrogate spin evolution. So, using the PN spins at
                        omega_switch_IG, the surrogate dynamics are evolved
                        both backwards and forwards in time.  Then the PN spins
                        at omega = omega_sur[t_sur_switch] are used to
                        reevaluate the surrogate dynamics. This way,
                        omega_switch_IG is only used as an initial guess for
                        where to switch to the surrogate, and we switch to the
                        surrogate as early as possible.
            t_sur_switch:
                        The time at which we switch to the surrogate evolution.
                        This should be slightly larger than the surrogate start
                        time of -4300M. Since we use the PN spins at
                        omega_sur[t_sur_switch] to reevaluate the surrogate
                        dynamics, it is possible that for the given PN spins at
                        omega_sur[t_sur_switch], the surrogate is not long
                        enough. Default: -4000.
            return_spin_evolution:
                        Also return the spin evolution and dynamics as a
                        dictionary.
            unlimited_extrapolation:
                        If True, Allow unlimited extrapolation to regions well
                        outside the surrogate's training region. Else, raises
                        an error.
        Outputs:
            chiA_coorb_fitnode:
                        Spin of BhA in the coorbital frame at t=-100M.
            chiB_coorb_fitnode:
                        Spin of BhB in the coorbital frame at t=-100M.
            quat_fitnode:
                        Coprecessing frame quaternions at t=-100M.
            orbphase_fitnode:
                        Orbital phase in the coprecessing frame at t=-100M.
            spin_evolution:
                        Dictionary containing the PN and surrogate spin
                        evolution (in the inertial frame) and dynamics.
                        Default: None, unless return_spin_evolution=True.
        """

        # If omega0 is below the NRSur7dq4 initial guess frequency, we use PN
        # to evolve the spins. We get the initial PN spins and omega_init_sur
        # that should go into the surrogate such that the initial time is
        # t_sur_switch.
        if omega0 < omega_switch_IG:
            chiA0_nrsur_coorb, chiB0_nrsur_coorb, quat0_nrsur_copr, \
                phi0_nrsur, omega_init_sur, chiA_PN, chiB_PN, omega_PN \
                = self._get_pn_spins_at_surrogate_start(q, \
                chiA0, chiB0, omega0, omega_switch_IG, t_sur_switch,
                unlimited_extrapolation)

        # If omega0 >= omega_switch_IG, we evolve spins directly with NRSur7dq4
        # waveform model. We set the coprecessing frame quaternion to identity
        # and orbital phase to 0 at omega=omega0, hence the coorbital frame is
        # the same as the inertial frame here.
        else:
            # Note that here we set omega_init_sur to omega0
            chiA0_nrsur_coorb, chiB0_nrsur_coorb, quat0_nrsur_copr, \
                phi0_nrsur, omega_init_sur, chiA_PN, chiB_PN, omega_PN \
                = chiA0, chiB0, [1,0,0,0], 0, omega0, None, None, None

        # Now evaluate the surrogate dynamics using PN spins at omega_init_sur
        dyn_times, quat_sur, orbphase_sur, chiA_copr_sur, chiB_copr_sur \
            = self._get_surrogate_dynamics(q, chiA0_nrsur_coorb, \
            chiB0_nrsur_coorb, quat0_nrsur_copr, phi0_nrsur, \
            omega_init_sur, unlimited_extrapolation)

        # get data at time node where remnant fits are done
        fitnode_time = -100
        nodeIdx = np.argmin(np.abs(dyn_times - fitnode_time))
        quat_fitnode = quat_sur.T[nodeIdx]
        orbphase_fitnode = orbphase_sur[nodeIdx]

        # get coorbital frame spins at the time node
        chiA_coorb_fitnode = quaternion_utils.rotate_in_plane(
                chiA_copr_sur[nodeIdx], orbphase_fitnode)
        chiB_coorb_fitnode = quaternion_utils.rotate_in_plane(
                chiB_copr_sur[nodeIdx], orbphase_fitnode)

        if return_spin_evolution:
            # Transform spins to the reference inertial frame
            chiA_sur = quaternion_utils.transform_time_dependent_vector(
                    quat_sur, chiA_copr_sur.T).T
            chiB_sur = quaternion_utils.transform_time_dependent_vector(
                    quat_sur, chiB_copr_sur.T).T
            spin_evolution = {
                    't_sur': dyn_times,
                    'chiA_sur': chiA_sur,
                    'chiB_sur': chiB_sur,
                    'orbphase_sur': orbphase_sur,
                    'quat_sur': quat_sur,
                    'omega_PN': omega_PN,
                    'chiA_PN': chiA_PN,
                    'chiB_PN': chiB_PN,
                    'omega_init_sur': omega_init_sur,
                }
        else:
            spin_evolution = None

        return chiA_coorb_fitnode, chiB_coorb_fitnode, quat_fitnode, \
                orbphase_fitnode, spin_evolution

    # ------------------------------------------------------------------------
    def _get_fit_params(self, m1, m2, chiA_vec, chiB_vec, f_ref,
            extra_params_dict):
        """
        If f_ref = -1, assumes the reference epoch is at t=-100M and
            just returns the input spins, which are expected to be in
            the coorbital frame at t=-100M.
        Else, takes the spins at f_ref and evolves them using a combination of
        PN and NRSur7dq4 until t=-100M. The returned spins are in the coorbital
        frame at t=-100M.

        See eval_fits.eval_nrfit() for the definitions of the arguments of
        this function.
        """

        unlimited_extrapolation = extra_params_dict["unlimited_extrapolation"]

        q = m1/m2
        M = m1 + m2     # m1, m2 are assumed to be in kgs

        if f_ref == -1:
            # If f_ref = -1, assume the reference epoch is at t=-100M from the
            # total waveform amplitude peak. All input and output quantities
            # are in the coorbital frame at this time.
            chiA_coorb_fitnode = chiA_vec
            chiB_coorb_fitnode = chiB_vec
            # We set this to None, and this will be used in _eval_fit() to
            # determine if the remnant vectors need any further transformations.
            quat_fitnode = None
            orbphase_fitnode = None
        else:
            # Else, evolve the spins from f_ref to t = -100 M from the peak.

            # reference orbital angular frequency times the total mass (M*omega)
            omega_ref = f_ref*np.pi*M*lal.G_SI/lal.C_SI**3

            chiA_coorb_fitnode, chiB_coorb_fitnode, quat_fitnode, \
                orbphase_fitnode, spin_evolution = self._evolve_spins(q, \
                chiA_vec, chiB_vec, omega_ref, \
                unlimited_extrapolation=unlimited_extrapolation)

        fit_params = [q, chiA_coorb_fitnode, chiB_coorb_fitnode, quat_fitnode, \
            orbphase_fitnode]

        return fit_params

    # ------------------------------------------------------------------------
    def _eval_fit(self, fit_params, fit_type, extra_params_dict):
        """ Evaluates a particular fit for NRSur7dq4Remnant using the fit_params
        returned by _get_fit_params().

        fit_type can be one of "FinalMass", "FinalSpin" or "RecoilKick".

        Passing extra_params_dict = {"unlimited_extrapolation": True}, will
        ignore any extrapolation errors. USE AT YOUR OWN RISK!! Default: False.
        """

        q, chiA_vec, chiB_vec, quat_fitnode, orbphase_fitnode = fit_params
        LALParams = lal.CreateDict()
        if extra_params_dict["unlimited_extrapolation"]:
            lal.DictInsertUINT4Value(LALParams, "unlimited_extrapolation", 1)

        if fit_type == "FinalMass":
            # FinalMass is given as a fraction of total mass
            val = lalsim.NRSur7dq4Remnant(q,
                chiA_vec[0], chiA_vec[1], chiA_vec[2],
                chiB_vec[0], chiB_vec[1], chiB_vec[2], "mf",
                LALParams).data[0]
        else:
            if fit_type == "FinalSpin":
                val = lalsim.NRSur7dq4Remnant(q,
                    chiA_vec[0], chiA_vec[1], chiA_vec[2],
                    chiB_vec[0], chiB_vec[1], chiB_vec[2], "chif",
                    LALParams).data
            elif fit_type == "RecoilKick":
                val = lalsim.NRSur7dq4Remnant(q,
                    chiA_vec[0], chiA_vec[1], chiA_vec[2],
                    chiB_vec[0], chiB_vec[1], chiB_vec[2], "vf",
                    LALParams).data
            else:
                raise ValueError("Invalid fit_type=%s. This model only allows "
                    "'FinalMass', 'FinalSpin' and 'RecoilKick'."%fit_type)

            # The fits are constructed in the coorbital frame at t=-100M.
            # If quat_fitnode is None, this means that no spin evolution
            # was done and the return values should be in the coorbital frame
            # at t=-100. So, we don't need any further transformations.
            if quat_fitnode is not None:
                # If quat_fitnode is not None, we transform the remnant vectors
                # into the LAL inertial frame, which is the same as the
                # coorbital frame at omega0.  This is done using the
                # coprecessing frame quaternion and orbital phase at t=-100M.
                # These are defined w.r.t. the reference LAL inertial frame
                # defined at omega0.
                val = quaternion_utils.transform_vector_coorb_to_inertial(val, \
                    orbphase_fitnode, quat_fitnode)

        return np.atleast_1d(val)
