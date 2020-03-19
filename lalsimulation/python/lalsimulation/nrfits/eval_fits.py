""" Evaluates Numerical Relativity fits such as remnant BH mass, spin, etc, for
    various models.

Vijay Varma, 2019.
"""
import numpy as np
import warnings

from lal import MSUN_SI

from .NRSur7dq4Remnant import NRSur7dq4Remnant
from .NRSur3dq8Remnant import NRSur3dq8Remnant
from . import quaternion_utils

#=============================================================================
class FitAttributes(object):
    """ Saves attributes of a particular fit.
    """
    #-------------------------------------------------------------------------
    def __init__(self, fit_class, desc, refs, model_keywords):
        """ Initializes a FitAttributes object. See usage examples below.

            - fit_class: Class to load for a given model. type: should be a
                daughter class of nrfits.NRFits.
            - desc: A brief description of the model, type: str.
            - refs: Paper reference, type: str.
            - model_keywords: List of strings, with one or more of
                allowed_model_keywords.
        """

        # Types of model keywords allowed (add more if required)
        allowed_model_keywords = ['nonspinning', 'aligned_spin', 'precessing', \
            'eccentric', 'tidal']

        for key in model_keywords:
            if key not in allowed_model_keywords:
                raise Exception("Invalid model_keyword: %s"%key)

        self.fit_class = fit_class
        self.desc = desc
        self.refs = refs
        self.model_keywords = model_keywords


#=============================================================================

# Collection of all implemented fits
fits_collection = {}

fits_collection['NRSur3dq8Remnant'] = FitAttributes( \
    fit_class = NRSur3dq8Remnant(),
    desc = 'Fits for remnant mass, spin and kick velocity for nonprecessing'
        ' BBH systems, trained on mass ratios q <= 8 and dimensionless spin '
        'magnitudes <= 0.81. This model was called surfinBH3dq8 in the paper.',
    refs = 'arxiv:1809.09125',
    model_keywords = ['aligned_spin'],
    )

fits_collection['NRSur7dq4Remnant'] = FitAttributes( \
    fit_class = NRSur7dq4Remnant(),
    desc = 'Fits for remnant mass, spin and kick velocity for generically'
        ' precessing BBH systems, trained on mass ratios q <= 4.01 and '
        'dimensionless spin magnitudes <= 0.81.',
    refs = 'arxiv:1905.09300',
    model_keywords = ['precessing'],
    )


#=============================================================================
def truncate_output_to_physical_limits(val_dict, behavior):
    """
    Function to truncate fit values when they go beyond the following physical
    limits:
        - FinalMass: 1, FinalMass is expected as a fraction of total mass here.
        - FinalSpin: magnitude = 1, the Kerr limit for dimensionless spins.
        - RecoilKick: magnitude = 1, as this is in units of c.

    Inputs: \n
    - val_dict: A dictionary of fits with one or more of the above keys. See
        eval_nrfit.

    - behavior (str) can be one of the following:
        * "TRUNCATE":         Rescale output magnitude to maximum limit, with an
                            info message.
        * "TRUNCATESILENT":   Rescale output magnitude to maximum limit, without
                            info message
        * "KEEP":             Keep output magnitude, even if unphysical, but
                            still give the info message
        * "IGNORE":           Keep output magnitude, even if unphysical, without
                            info message
        * "ERROR":            Abort with an error if output magnitude is
                            unphysical.
    """

    # Don't need to do anything
    if behavior == "IGNORE":
        return val_dict

    allowed_behaviors = ["ERROR","IGNORE","KEEP","TRUNCATE","TRUNCATESILENT"]
    if behavior not in allowed_behaviors:
        raise ValueError("Invalid physical_limit_violation_behavior=%s"%behavior
            + ". Should be one of ["+ ", ".join(allowed_behaviors) + "].")

    physical_limits = {
        "FinalMass": 1,     # Remnant can't be heavier than the total mass
        "FinalSpin": 1,     # Kerr limit
        "RecoilKick": 1,    # Speed of light limit
        }

    for fit_type in val_dict.keys():
        limit = physical_limits[fit_type]
        magnitude = np.linalg.norm(val_dict[fit_type])
        if magnitude > limit:
            # Error/Warning message
            message = ""
            if fit_type in ["FinalSpin", "RecoilKick"]:
                message += "%s magnitude=%.7e"%(fit_type, magnitude)
            else:
                message += "Dimensionless %s=%.7e"%(fit_type, magnitude)
            message += " is over max limit=%.2e."%(limit)

            if behavior == "ERROR":             # raise error
                raise RuntimeError(message + " Adapt " \
                "extra_params_dict['physical_limit_violation_behavior'] to "\
                "continue without an error.")

            elif behavior in ["TRUNCATE", "TRUNCATESILENT"]: # rescale magnitude
                message += " Setting to max limit."
                val_dict[fit_type] *= limit/magnitude

            # Only these cases get a warning
            if behavior in ["KEEP", "TRUNCATE"]:
                warnings.warn(message)

    return val_dict

#=============================================================================
def check_extra_params_and_set_defaults(extra_params_dict):
    """ Does some sanity checks on extra_params_dict.
    If any of the default_keywords are not specified, sets them to the
    default values.
    """

    # NOTE: To add a new key to extra_params_dict, set a default value here
    # and update the documentation of eval_nrfit
    default_keywords = {
        'Lambda1': None,
        'Lambda2': None,
        'eccentricity': None,
        'mean_anomaly': None,
        'unlimited_extrapolation': False,
        'physical_limit_violation_behavior': "ERROR",
        }

    if extra_params_dict is None:
        extra_params_dict = {}

    # Sanity checks
    for key in extra_params_dict.keys():
        if key not in default_keywords.keys():
            raise ValueError('Invalid key %s in extra_params_dict. '%(key)
            + "Should be one of ["+ ", ".join(default_keywords.keys()) + "].")

    # set to default if keyword is not specified
    for key in default_keywords:
        if key not in extra_params_dict.keys():
            extra_params_dict[key] = default_keywords[key]

    return extra_params_dict



#=============================================================================
def eval_nrfit(m1, m2, chiA_vec, chiB_vec, model_name, fit_types_list, f_ref=-1,
        extra_params_dict=None):
    """
    Evaluates Numerical Relativity fits for a given model.

    - m1:
        mass of object 1 in kg. \n
        This needs to be in kg for models that allow f_ref != -1. When
        f_ref = -1, if other units are used for m1/m2, the remnant mass will be
        returned in the same units. \n
    - m2:
        mass of the object 2 in kg. \n
    - chiA_vec:
        dimensionless spin (3-vector) of object 1 at the reference epoch. \n
    - chiB_vec:
        dimensionless spin (3-vector) of object 2 at the reference epoch. \n
        This follows the same convention as the waveform interface (see
        https://dcc.ligo.org/T1800226/public), where the spin components are
        defined as:
        * \f$\chi_z = \chi \cdot \hat{L}\f$, where \f$L\f$ is the orbital
            angular momentum vector at the reference epoch.
        * \f$\chi_x = \chi \cdot \hat{n}\f$, where \f$n =\f$ body2 -> body1
            is the separation vector pointing from body2 to body1 at the
            reference epoch.
        * \f$\chi_y = \chi \cdot (\hat{L} \times \hat{n})\f$. \n
        These spin components are frame-independent as they are defined using
        vector inner products. One can also think of the above spins as being
        defined in the following reference frame: The positive z-axis is along
        the orbital angular momentum at the reference epoch. The separation
        vector from the object 2 to object 1 at the reference epoch is along
        the positive x-axis. The y-axis completes the right-handed triad. The
        returned vectors such as final spin and kick velocity are also defined
        in the same frame.
    - model_name:
        model used for fit. \n
        See lalsimulation.nrfits.eval_fits.fits_collection.keys() for all
        available fits. Details about a particular model, including its
        validity, can be found by doing the following: \n
            >> from lalsimulation.nrfits.eval_fits import fits_collection \n
            >> help(fits_collection[model_name].fit_class)
    - fit_types_list:
        List of fits to evaluate. \n
        Allowed values for elements are
        "FinalMass", "FinalSpin", "RecoilKick", and "PeakLuminosity". \n
        Example: fit_types_list = ["FinalMass", "FinalSpin"]
    - f_ref:
        reference frequency (in Hz) used to set the reference epoch. \n
        The reference epoch is set such that the orbital frequency is equal to
        f_ref/2. If f_ref = -1, the reference epoch is taken to be the
        time/frequency at which the fits were constructed. This can be
        different for different models, see the documentation for the specific
        model. Default: f_ref = -1.
    - extra_params_dict:
        Any additional parameters required for a specific model.
        Default: None. \n
        Allowed keys for extra_params_dict are:
        * Lambda1:
            Tidal deformability of object 1. Default: None.
        * Lambda2:
            Tidal deformability of object 2. Default: None.
        * eccentricity:
            Eccentricity at the reference epoch. Default: None.
        * mean_anomaly:
            Mean anomaly at the reference epoch. Default: None.
        * unlimited_extrapolation:
            Some models will raise an error if evaluated well outside
            their calibration range. If unlimited_extrapolation = True,
            then these errors are not raised and the model will still
            produce an output. USE AT YOUR OWN RISK!! Default: False.
        * physical_limit_violation_behavior:
            What to do if the fit values violate physical limits
            FinalMass > total mass, |FinalSpin| > 1, or |RecoilKick| > 1. \n
            Allowed options are (Default: "ERROR"):
            - "ERROR":
                Abort with an error message.
            - "TRUNCATE":
                Rescale output magnitude to maximum limit, with a warning.
            - "TRUNCATESILENT":
                Rescale output magnitude to maximum limit, without
                warning.
            - "KEEP":
                Keep output magnitude, but with a warning.
            - "IGNORE":
                Keep output magnitude, without warning.


    - Returns: return_dict. \n
        Dictionary of values corresponding to the keys in fit_types_list. \n
        Example: mf = return_dict["FinalMass"]
    """

    ### Sanity checks
    if model_name not in fits_collection.keys():
        raise ValueError("Invalid model_name=%s. "%model_name \
            + "Should be one of ["+ ", ".join(fits_collection.keys()) + "].")

    if m1 < 0.09 * MSUN_SI and f_ref != -1:
        warnings.warn("Small value of m1 = %e (kg) = %e (Msun) requested. "
            "When f_ref != -1, component masses must be in kgs, perhaps you "
            "are using different units?"%(m1, m1/MSUN_SI))

    if m2 < 0.09 * MSUN_SI and f_ref != -1:
        warnings.warn("Small value of m2 = %e (kg) = %e (Msun) requested. "
            "When f_ref != -1, component masses must be in kgs, perhaps you "
            "are using different units?"%(m2, m2/MSUN_SI))

    if m1 <= 0 or m2 <= 0:
        raise ValueError("Got nonpositive mass: m1=%.3e, m2=%.3e"%(m1,m2))

    chiA_vec = np.atleast_1d(chiA_vec)
    chiB_vec = np.atleast_1d(chiB_vec)
    if len(chiA_vec) != 3 or len(chiB_vec) != 3:
        raise TypeError("Expected input spins to be 3-vectors.")

    if np.linalg.norm(chiA_vec) > 1:
        raise ValueError("Invalid spin magnitude |chiA_vec|=%.3f."%( \
            np.linalg.norm(chiA_vec)))

    if np.linalg.norm(chiB_vec) > 1:
        raise ValueError("Invalid spin magnitude |chiB_vec|=%.3f."%( \
            np.linalg.norm(chiB_vec)))

    if not type(fit_types_list) == list:
        raise TypeError("fit_types_list should be a list.")


    # do sanity checks on extra_params_dict and set default values if
    # required.
    extra_params_dict = check_extra_params_and_set_defaults(extra_params_dict)

    # Some further sanity checks to makes sure extra_params_dict is
    # compatible with the given model
    if (extra_params_dict['Lambda1'] is not None) \
            or (extra_params_dict['Lambda2'] is not None):
        if 'tidal' not in fits_collection[model_name].model_keywords:
            raise ValueError("This model does not allow Lambda1/Lambda2.")

    if (extra_params_dict['eccentricity'] is not None) \
            or (extra_params_dict['mean_anomaly'] is not None):
        if 'eccentric' not in fits_collection[model_name].model_keywords:
            raise ValueError("This model does not allow eccentricity or "
                    "mean_anomaly.")

    if 'aligned_spin' in fits_collection[model_name].model_keywords:
        if np.linalg.norm(chiA_vec[:2]) > 0 or np.linalg.norm(chiB_vec[:2]) > 0:
            raise ValueError("This model only allows nonprecessing spins")

    if 'nonspinning' in fits_collection[model_name].model_keywords:
        if np.linalg.norm(chiA_vec) > 0 or np.linalg.norm(chiB_vec) > 0:
            raise ValueError("This model only allows zero spins")


    swapped_labels = False
    # If m1 < m2, we switch the labels of the two objects, and then rotate the
    # spins by pi about the z-direction. This amounts to a rigid rotation of
    # the full system about the z-axis by pi. After computing the remnant
    # properties, the final spin and kick vectors are rotated by pi to undo
    # this switch.
    if m1 < m2:
        temp = m1
        m1 = m2
        m2 = temp
        temp = chiA_vec
        chiA_vec = chiB_vec
        chiB_vec = temp
        chiA_vec = quaternion_utils.rotate_in_plane(chiA_vec, np.pi)
        chiB_vec = quaternion_utils.rotate_in_plane(chiB_vec, np.pi)
        swapped_labels = True

    # Compute NR fit quantities
    val_dict = fits_collection[model_name].fit_class(m1, m2, chiA_vec, chiB_vec,
            f_ref, fit_types_list, extra_params_dict)

    # If the output breaks physical limits, truncate it if requested.
    val_dict = truncate_output_to_physical_limits(val_dict, \
        extra_params_dict['physical_limit_violation_behavior'])

    # So far FinalMass was dimless, now rescale to same units as total mass
    if 'FinalMass' in val_dict.keys():
        val_dict['FinalMass'] *= (m1 + m2)

    # Rotate remnant vectors by pi if the component lables were swapped
    if swapped_labels:
        if 'FinalSpin' in val_dict.keys():
            val_dict['FinalSpin'] = quaternion_utils.rotate_in_plane( \
                val_dict['FinalSpin'], np.pi)
        if 'RecoilKick' in val_dict.keys():
            val_dict['RecoilKick'] = quaternion_utils.rotate_in_plane( \
                val_dict['RecoilKick'], np.pi)

    return val_dict
