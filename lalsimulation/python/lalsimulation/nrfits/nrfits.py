"""
Base class for Numerical Relativity fits such as remnant BH mass, spin, etc.

Vijay Varma, 2019.
"""

#=============================================================================
class NRFits(object):
    """
    Base class to evaluate NR fits.

    For each new fit, you need to do the following:
        1. Add a new derived class for this class in a separate file and add the
            new class to fits_collection in eval_fits.py.
        2. override _get_fit_params() and _eval_fit() in the new derived class.
        3. Add the new filename to Makefile.am

    See NRSur7dq4Remnant.py for an example.
    """

    # ------------------------------------------------------------------------
    def _get_fit_params(self, m1, m2, chiA_vec, chiB_vec, f_ref,
            extra_params_dict):
        """ Not all fits will require all of these parameters, this function
        should be used to get the required params, and if necessary, process
        them to get the parameters that are used to evaluate the actual fits.

        For example: chiA_vec/chiB_vec are defined at f_ref, this function
        could take these initial spins and evolve them to get the spins at
        ISCO, which are then used to evaluate the fits.

        See eval_fits.eval_nrfit() for the definitions of the arguments of
        this function.
        """
        raise NotImplementedError("Please override me.")
        return fit_params

    # ------------------------------------------------------------------------
    def _eval_fit(self, fit_params, fit_type, extra_params_dict):
        """ Evaluates a particular fit for a given model using the fit_params
        returned by _get_fit_params().
        """
        raise NotImplementedError("Please override me.")

    # ------------------------------------------------------------------------
    def __call__(self, m1, m2, chiA_vec, chiB_vec, f_ref, fit_types_list,
            extra_params_dict):
        """ Evaluates all fits given in fit_types_list and returns them as
            a dictionary.

            See eval_fits.eval_nrfit() for the definitions of the arguments of
            this function.
        """

        # Get the parameters used to evaluate the fits. This includes any
        # spin evolution that may be required
        fit_params = self._get_fit_params(m1, m2, chiA_vec, chiB_vec,
            f_ref, extra_params_dict)

        # Add more if needed and update the fit_types_list documentation in
        # eval_fits.py
        allowed_fit_types = [
                "FinalMass",
                "FinalSpin",
                "RecoilKick",
                "PeakLuminosity",
                ]

        # Loop over fit_types_list and return a dictionary of values
        return_dict = {}
        for fit_type in fit_types_list:
            if fit_type not in allowed_fit_types:
                raise ValueError("Invalid fit_type=%s. "%fit_type \
                    + "Should be one of ["+ ", ".join(allowed_fit_types) + "].")
            return_dict[fit_type] = self._eval_fit(fit_params, fit_type,
                    extra_params_dict)

        return return_dict


