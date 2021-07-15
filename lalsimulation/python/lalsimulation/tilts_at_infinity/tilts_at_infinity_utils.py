"""
This contains utility functions and definitions used in calc_tilts_prec_avg_regularized.py and
hybrid_spin_evolution.py

N. K. Johnson-McDaniel, 2021
"""

import numpy as np

from lal import G_SI, C_SI, MSUN_SI

import warnings
from warnings import warn

# Define the kg to m and s conversions

kg_to_m = G_SI/C_SI**2

kg_to_s = G_SI/C_SI**3

# Setup for error handling

class Error(Exception):
    """Base class for exceptions in this module"""
    pass

class NonprecessingError(Error):
    """Exception raised when the evolution finds that the system is nonprecessing to numerical accuracy."""

    def __init__(self, message):
        self.message = message

# Setup for warnings

class Warn(UserWarning):
    """Base class for warnings in this module"""
    pass

class ValueWarning(Warn):
    """Warning raised when an input is dubious"""
    pass

# Make warnings cleaner

def warning_formatter(message, category, filename, lineno, file=None, line=None):
    return '%s: %s\n'%(category.__name__, message)

warnings.formatwarning = warning_formatter

# Define functions


def format_error(err):
        return "%s: %s"%(type(err).__name__, err)


def check_masses(m1, m2):
    if m1 <= 0. or m2 <= 0.:
        raise ValueError(
            "The input masses must both be positive, while they are m1, m2 = %e, %e kg" % (m1, m2))

    if m1 < 0.09*MSUN_SI or m2 < 0.09*MSUN_SI:
        warn("One or both of the masses is rather small: m1 = %e kg = %e Msun, m2 = %e kg = %e Msun. "
             "The masses must be input in kg." % (m1, m1/MSUN_SI, m2, m2/MSUN_SI), ValueWarning)


def eq_mass_check(m1, m2, Lf):
    if m1 == m2:
        if Lf is None:
            raise ValueError(
                "Cannot compute tilts at infinity for exactly equal-mass binaries, as those quantities are not well "
                "defined in this case.")
        else:
            raise ValueError(
                "The computation of the bounds and average value of the tilts is not yet implemented for exactly "
                "equal-mass binaries.")


def check_spin_mags(chi1, chi2):
    if chi1 < 0. or chi1 > 1. or chi2 < 0. or chi2 > 1.:
        raise ValueError(
            "The magnitudes of the spins must both be between 0 and 1, while they are chi1, chi2 = %f, %f" % (chi1, chi2))


def check_tilts(tilt1, tilt2):
    if tilt1 < 0. or tilt1 > np.pi or tilt2 < 0. or tilt2 > np.pi:
        raise ValueError(
            "The tilt angles must both be between 0 and pi, while they are tilt1, tilt2 = %f, %f" % (tilt1, tilt2))


def check_fref(fref, m1, m2, evol_type):
    if fref <= 0.:
        raise ValueError(
            "The reference frequency must be positive, while it is fref = %f Hz" % fref)

    f_ISCO = 6.**(-1.5)/(np.pi*(m1 + m2)*kg_to_s)

    if fref > f_ISCO:
        warn("The reference frequency should not be close to merger, where the %s evolution is highly inaccurate, "
             "while it is fref = %f Hz, which is greater than the Schwarzschild ISCO frequency associated "
             "with the binary's total mass of %f Hz"%(evol_type, fref, f_ISCO), ValueWarning)


def package_tilts(tilt1, tilt2, Lf, swap):
    """
    Package tilts to be returned by prec_avg_tilt_comp_vec_inputs() or prec_avg_tilt_comp() depending on whether
    Lf is None or not. Also swaps the tilts if necessary. Only works in the case when Lf is not None when the
    min, max, and average values are all the same.

    Inputs:

    tilt1, tilt2: Tilts
    Lf: Final orbital angular momentum (here just acts as a switch depending on whether it is None or not)
    swap: Whether to swap tilt1 and tilt2 before returning (True) or not (False)

    Output: dictionary with entries 'tilt1_inf', 'tilt2_inf' for evolution to infinity--Lf is None--and entries
            'tilt1_sep_min', 'tilt1_sep_max', 'tilt1_sep_avg', 'tilt2_sep_min', 'tilt2_sep_max', 'tilt2_sep_avg'
            for evolution to a finite separation (i.e., a finite orbital angular momentum), when Lf is not None
    """

    if swap:
        tilt1, tilt2 = tilt2, tilt1

    if Lf is None:
        return {'tilt1_inf': tilt1, 'tilt2_inf': tilt2}
    else:
        return {'tilt1_sep_min': tilt1, 'tilt1_sep_max': tilt1, 'tilt1_sep_avg': tilt1, 'tilt2_sep_min': tilt2,
                'tilt2_sep_max': tilt2, 'tilt2_sep_avg': tilt2}

def evolution_error_handling(failure_mode, failure_message, failure_output, failure_output_string, Lf, swap):
    """
    Take care of the error message or warning and returning something for the tilts when the
    precession-averaged evolution fails

    Inputs:
    failure_mode: The failure mode (either "Error", "NAN", or "None")
    failure_message: The message to print in the error or warning
    failure output: What to output for the tilts in the event of a failure
    failure_output_string: The string associated with the failure output
    Lf: Final orbital angular momentum (here just acts as a switch depending on whether it is None or not)
    swap: Whether to swap tilt1 and tilt2 before returning (True) or not (False)

    Output: dictionary with entries 'tilt1_inf', 'tilt2_inf' for evolution to infinity--Lf is None--and entries
            'tilt1_sep_min', 'tilt1_sep_max', 'tilt1_sep_avg', 'tilt2_sep_min', 'tilt2_sep_max', 'tilt2_sep_avg'
            for evolution to a finite separation (i.e., a finite orbital angular momentum), when Lf is not None.
            The entries of the dictionary are failure_output
    """

    if failure_mode == 'Error':
        raise RuntimeError(failure_message)
    else:
        failure_message += " Returning %s for the tilts."%failure_output_string

        warn(failure_message, RuntimeWarning)

        return package_tilts(failure_output, failure_output, Lf, swap)
