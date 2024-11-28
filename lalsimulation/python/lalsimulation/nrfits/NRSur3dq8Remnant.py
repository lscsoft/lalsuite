"""
Class for NRSur3dq8Remnant model for the remnant mass, spin and kick velocity
for nonprecessing BBH systems. This model was called surfinBH3dq8 in the paper.

Vijay Varma, 2019.
"""

import numpy as np
import lalsimulation as lalsim
import lal

from .nrfits import NRFits

class NRSur3dq8Remnant(NRFits):
    r"""
    Class for NRSur3dq8Remnant model for the remnant mass, spin and kick
    velocity for nonprecessing BBH systems. This model was called surfinBH3dq8
    in the paper.

    Paper: arxiv:1809.09125. The model is referred to as surfinBH3dq8 in the
    paper.

    Parameter ranges for usage: \n
        q = [1, 9.1] \n
        \f$\chi_{1z}, \chi_{2z}\f$ = [-0.91, 0.91] \n
        OR \n
        q = [1, 10.1] \n
        \f$\chi_{1z}, \chi_{2z}\f$ = [-0.81, 0.81]

    Training parameter ranges:  \n
        q = [1, 8] \n
        \f$\chi_{1z}, \chi_{2z}\f$ = [-0.81, 0.81]

    But extrapolates reasonably to the above mass ratios and spins. However,
    if a guarantee of accuracy is required, this model should be used within
    the training parameter range.
    """

    # ------------------------------------------------------------------------
    def _get_fit_params(self, m1, m2, chiA_vec, chiB_vec, f_ref,
            extra_params_dict):
        """ No processing or spin evolution is required here, just returns the
        mass ratio, total mass and spins along the z-direction.

        See eval_fits.eval_nrfit() for the definitions of the arguments of
        this function.
        """

        # f_ref is used to set the reference frame. For this model, only the
        # final kick components depend on the reference frame. However, to set
        # this at an arbitrary f_ref, we will need to know the orbital phase
        # difference between f_ref and the time at which the fits are
        # constructed, i.e. t=-100M from the peak of the total waveform
        # amplitude. While this could be done using a waveform model, we don't
        # expect the kick direction of this model to be very useful.  So, we
        # only allow f_ref = -1, which assumes that the reference epoch is at
        # t=-100M. The kick direction can still be used, but it is always
        # retured in this frame.
        if f_ref != -1:
            raise ValueError("This model only works for f_ref=-1.")

        q = m1/m2
        fit_params = [q, chiA_vec[2], chiB_vec[2]]
        return fit_params

    # ------------------------------------------------------------------------
    def _eval_fit(self, fit_params, fit_type, extra_params_dict):
        """ Evaluates a particular fit for NRSur3dq8Remnant using the fit_params
        returned by _get_fit_params().
        """
        q, chiAz, chiBz = fit_params
        LALParams = lal.CreateDict()
        if extra_params_dict["unlimited_extrapolation"]:
            lal.DictInsertUINT4Value(LALParams, "unlimited_extrapolation", 1)

        if fit_type == "FinalMass":
            # FinalMass is given as a fraction of total mass
            val = lalsim.NRSur3dq8Remnant(q, chiAz, chiBz, "mf", LALParams)
        elif fit_type == "FinalSpin":
            # chifx and chify are zero for aligned-spin systems
            chifz = lalsim.NRSur3dq8Remnant(q, chiAz, chiBz, "chifz", LALParams)
            val = [0,0,chifz]
        elif fit_type == "RecoilKick":
            # vfz is zero for aligned-spin systems
            vfx = lalsim.NRSur3dq8Remnant(q, chiAz, chiBz, "vfx", LALParams)
            vfy = lalsim.NRSur3dq8Remnant(q, chiAz, chiBz, "vfy", LALParams)
            val = [vfx, vfy, 0]
        else:
            raise ValueError("Invalid fit_type=%s. This model only allows "
                "'FinalMass', 'FinalSpin' and 'RecoilKick'."%fit_type)

        return np.atleast_1d(val)
