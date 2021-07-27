"""
Implement the computation of the tilt angles at infinity or bounds and average values of tilts at a finite separation
from tilt angles at some frequency, taken to be sufficiently small so that precession-averaged evolution is accurate.
This implementation is described in the paper, <https://dcc.ligo.org/P2100029>, arXiv:2107.11902.

N. K. Johnson-McDaniel, 2021
"""

import numpy as np

from numpy.linalg import norm

from numpy.polynomial import Polynomial as P

from scipy.special import ellipe, ellipk

from scipy.integrate import ode

from .tilts_at_infinity_utils import *

from warnings import warn

try:
    from scipy.integrate import solve_ivp
except ImportError:
    solve_ivp_avail = False
else:
    solve_ivp_avail = True

try:
    import mpmath as mp
except ImportError:
    mpmath_avail = False
else:
    mpmath_avail = True

# Set available PN orders for the orbital angular momentum expression in terms of frequency

LPN_orders = [0, 1, 1.5, 2, 2.5]

LPN_order_string = ', '.join(map(str, LPN_orders))

# Define functions


def L_vec_from_masses_spins_and_f(m1, m2, f, PNorder=2, S1_vec=np.zeros(3), S2_vec=np.zeros(3),
                                  no_mass_spin_input_check=False):
    """
    Compute the PN orbital angular momentum from the binary's masses and spins and m = 2 GW frequency
    using the expression in Eq. (4) of <https://dcc.ligo.org/public/0122/T1500554/023/dLdS.pdf>.

    We orbit-average the contributions from the spins as in Eq. (17) of that reference, but the spins are
    set to zero by default.

    Inputs:

    - Required:

    m1, m2: Detector frame masses in kg
    f: m = 2 GW frequency in Hz

    - Optional:

    PNorder: Order of PN corrections to include--choices are 0, 1, 1.5, 2, or 2.5, where the half-integer orders are
             the spin-orbit contributions, default: 2
    S1_vec, S2_vec: 3d numpy arrays giving the components of the dimensionful spins (with total mass = 1) in a coordinate
                    system where the Newtonian orbital angular momentum is in the z-direction, default: [0., 0., 0.]
    no_mass_spin_input_check: Turn off mass and spin input checking for when this is called by other functions that have
                              already performed this check, default: False

    Output: 3d array giving the PN orbital angular momentum in total mass = 1 units in a coordinate system where the
            Newtonian orbital angular momentum is in the z-direction
    """

    # Check inputs

    if not no_mass_spin_input_check:
        check_masses(m1, m2)

    if f <= 0.:
        raise ValueError(
            "The frequency must be positive, while it is f = %f" % f)

    if PNorder not in LPN_orders:
        raise ValueError("The PN order to use to compute the orbital angular momentum must be one of %s, while it is "
                         "PNorder = %f" % (LPN_order_string, PNorder))

    M = m1 + m2
    Msq = M*M

    S1 = norm(S1_vec)
    S2 = norm(S2_vec)

    if not no_mass_spin_input_check and (S1 > m1*m1/Msq or S2 > m2*m2/Msq):
        raise ValueError(
            "The magnitudes of the spins must be at most (m1/M)^2 = %f and (m2/M)^2 = %f, respectively, "
            "where M is the total mass, while they are |S1|, |S2| = %f, %f" %
            (m1*m1/Msq, m2*m2/Msq, S1, S2))

    # Compute useful quantities

    X1 = m1/M
    X2 = m2/M

    X1sq = X1*X1
    X2sq = X2*X2

    eta = X1*X2

    zhat = np.array([0., 0., 1.])

    S1z_vec = S1_vec[2]*zhat
    S2z_vec = S2_vec[2]*zhat

    v = (np.pi*M*kg_to_s*f)**(1./3.)

    v2 = v*v
    v3 = v2*v

    # Put together L to the desired PN order, starting by giving the magnitude and direction of the Newtonian expression

    LN = eta/v

    L = 1.*zhat

    if PNorder >= 1:
        L += (1.5+eta/6.)*v2*zhat
    if PNorder >= 1.5:
        L += (-0.25*((1. + 3./X1)*S1_vec + (1. + 3./X2)*S2_vec)
              - ((1. + 27./X1)*S1z_vec + (1. + 27./X2)*S2z_vec)/12.
              )*v3
    if PNorder >= 2:
        L += (3.375 - 2.375*eta + eta*eta/24.)*v2*v2*zhat
    if PNorder >= 2.5:
        L += ((-71. + 21./X1 - 13.*X1 - 9.*X1sq)*S1_vec + (-71. + 21./X2 - 13.*X2 - 9.*X2sq)*S2_vec
              + (1. - 147./X1 - 367.*X1/3. + 13.*X1sq/3.)*S1z_vec
              + (1. - 147./X2 - 367.*X2/3. + 13.*X2sq/3.)*S2z_vec
              )*v3*v2/48.

    return LN*L


def L_from_masses_spins_and_f(m1, m2, f, PNorder=2, S1_vec=np.zeros(3), S2_vec=np.zeros(3), no_mass_spin_input_check=False):
    """
    Compute the magnitude of the PN orbital angular momentum from the binary's masses and spins and m = 2 GW frequency
    using the expression in Eq. (4) of <https://dcc.ligo.org/public/0122/T1500554/023/dLdS.pdf>.

    We orbit-average the contributions from the spins as in Eq. (17) of that reference, but the spins are
    set to zero by default.

    This wraps L_vec_from_masses_spins_and_f()

    Inputs:

    - Required:

    m1, m2: Detector frame masses in kg
    f: m = 2 GW frequency in Hz

    - Optional:

    PNorder: Order of PN corrections to include--choices are 0, 1, 1.5, 2, or 2.5, where the half-integer orders are
             the spin-orbit contributions, default: 2
    S1_vec, S2_vec: 3d numpy arrays giving the components of the dimensionful spins (with total mass = 1) in a coordinate
                    system where the Newtonian orbital angular momentum is in the z-direction, default: [0., 0., 0.]
    no_mass_spin_input_check: Turn off mass and spin input checking for when this is called by other functions that have
                              already performed this check, default: False

    Output: Magnitude of the PN orbital angular momentum in total mass = 1 units
    """

    return norm(L_vec_from_masses_spins_and_f(m1, m2, f, PNorder=PNorder, S1_vec=S1_vec, S2_vec=S2_vec,
                                              no_mass_spin_input_check=no_mass_spin_input_check))


def L_from_masses_a_and_ecc(m1, m2, a, e):
    """
    Compute the magnitude of the Newtonian orbital angular momentum from the binary's masses, semimajor axis, and
    eccentricity

    Inputs:

    m1, m2: Source frame masses in kg
    a: Semimajor axis in m
    e: eccentricity

    Output: Magnitude of orbital angular momentum in total mass = 1 units
    """

    # Check inputs

    check_masses(m1, m2)

    if a <= 0.:
        raise ValueError(
            "The semimajor axis must be positive, while it is a = %e" % a)

    if e < 0. or e > 1.:
        raise ValueError(
            "The eccentricity must be between 0 and 1, while it is e = %f" % e)

    # Compute

    M = m1 + m2

    eta = m1*m2/M/M

    return eta*(a*(1. - e*e)/(M*kg_to_m))**0.5


def E_ratio_term(m):
    """
    Compute m*E_ratio(m) [cf. Eq. (B13) in the paper <https://dcc.ligo.org/P2100029>, arXiv:2107.11902, which gives m^2*E_ratio(m)]
    """

    omm = 1. - m
    fpm = 4. + m

    return (3.*(1. + 3.*(4. - m)/(omm*fpm))/(16.*omm**1.5) - 0.5)*m/fpm


def eq_coeffs(u, kappaxiq, q, S1, S2, xi, S0sq):
    """
    Compute coefficients for the cubic equation in bar{S}^2 to be solved [Eq. (20) in the paper <https://dcc.ligo.org/P2100029>, arXiv:2107.11902]

    Inputs:
    u: 1/(2L), where L is the orbital angular momentum (i.e., the integration variable)
    kappaxiq: The rescaled variable defined in Eq. (14) of the paper
    q: Mass ratio of the binary, with the < 1 convention
    S1, S2: Magnitudes of the dimensionful spins of the binary (with total mass = 1)
    xi: Conserved effective spin, defined in Eq. (2) of the paper
    S0sq: Square of magnitude of the initial total spin (at u = u0, and with total mass = 1)

    Output: The coefficients of the equation (vareps, Btransf, Ctransf, Dtransf)
    """

    # Define useful quantities

    q2 = q*q

    eps = 1. - q

    opq = 1. + q

    u2 = u*u

    S1sq = S1*S1
    S2sq = S2*S2

    kappaxiq2 = kappaxiq*kappaxiq

    # Compute the coefficients of the transformed equation q (1 - q^2) u^2 bar{S}^6 + bar{B} bar{S}^4 +
    # bar{C} bar{S}^2 + bar{D} = 0 [Eqs. (20-22) of the paper <https://dcc.ligo.org/P2100029>, arXiv:2107.11902].

    # Here we refer to the barred coefficients with the suffix "transf" and the coefficient of bar{S}^6 as
    # vareps, as in Appendix C of the paper <https://dcc.ligo.org/P2100029>, arXiv:2107.11902

    Sigma = 0.5*(S0sq - S1sq - S2sq)
    Upsilon = opq*(2.*q*Sigma + q2*S1sq + S2sq)
    zeta = q*(Sigma + q*S1sq)*xi

    vareps = q*opq*eps*u2

    Btransf = 0.25*opq*eps*eps + q*eps * \
        (xi - 2.*opq*kappaxiq)*u + Upsilon*u2

    Ctransf = eps*(opq*(Sigma + q*kappaxiq2) - q*kappaxiq *
                   xi) + 2.*(zeta - Upsilon*kappaxiq)*u

    Dtransf = opq*(Sigma*Sigma - S1sq*S2sq) + q2*S1sq*xi * \
        xi/opq - 2.*zeta*kappaxiq + Upsilon*kappaxiq2

    return vareps, Btransf, Ctransf, Dtransf


def solve_cubic(coeffs, u, kappaxiq, q, S1, S2, imag_tol, use_mpmath=False, root_tol=1.e-6, polyroots_maxsteps=50, polyroots_extraprec=50):
    """
    Solve the cubic equation, check for complex roots, and order the roots.
    The inputs besides the coefficients are just needed in order to be printed
    if one of the checks is triggered

    Inputs:

    - Required:

    coeffs: Array of coefficients of cubic [vareps, Btransf, Ctransf, Dtransf]
    u: 1/(2L), where L is the orbital angular momentum (i.e., the integration variable)
    kappaxiq: The rescaled variable defined in Eq. (14) of the paper <https://dcc.ligo.org/P2100029>, arXiv:2107.11902
    q: Mass ratio of the binary, with the < 1 convention
    S1, S2: Magnitudes of the dimensionful spins of the binary (with total mass = 1)
    imag_tol: Tolerance for the imaginary part of roots of the cubic

    - Optional:

    use_mpmath: If True, use mpmath.polyroots to solve the cubic instead of numpy's Polynomial roots method
                (slower but more accurate), default: False
    root_tol: Tolerance for the error estimate from mpmath.polyroots, if use_mpmath == True, default: 1.e-6
    polyroots_maxsteps: The maximum number of steps allowed in mpmath.polyroots if use_mpmath == True, default: 50
    polyroots_extraprec: The number of extra digits to use in mpmath.polyroots if use_mpmath == True, default: 50

    Output: The roots of the equation ordered from smallest to largest (Ssq3, Ssqm, Ssqp)
    """

    # Extract coefficients and solve

    vareps, Btransf, Ctransf, Dtransf = coeffs

    if use_mpmath:
        Btransf_out, Ctransf_out, Dtransf_out = Btransf, Ctransf, Dtransf

        roots, err = mp.polyroots([vareps, Btransf, Ctransf, Dtransf], error=True, maxsteps=polyroots_maxsteps, extraprec=polyroots_extraprec)

        SsqI, SsqII, SsqIII = roots

        # Check that the error on the roots isn't too large
        if err > root_tol:
            warn("polyroots error of %e greater than root_tol = %e"%(err, root_tol), RuntimeWarning)
    else:
        # The [0] indices here and elsewhere are necessary to avoid a warning
        Btransf_out, Ctransf_out, Dtransf_out = Btransf[0], Ctransf[0], Dtransf[0]

        cubic = P([Dtransf_out, Ctransf_out, Btransf_out, vareps])

        SsqI, SsqII, SsqIII = cubic.roots()

    # Check that all roots are real
    if max(abs(SsqI.imag), abs(SsqII.imag), abs(SsqIII.imag)) > imag_tol:
        raise RuntimeError("At least one of the roots of the cubic in S^2 has an imaginary part larger than the "
                           "tolerance of %.2g set in the code.\nThe roots are (%.10f, %.10f), (%.10f, %.10f), "
                           "(%.10f, %.10f), and the coefficients are %.10e, %.10e, %.10e, %.10e at u = %.10f, "
                           "kappaxiq = %.10f. q = %f, S1, S2 = %f, %f" %
                           (imag_tol, SsqI.real, SsqI.imag, SsqII.real, SsqII.imag, SsqIII.real, SsqIII.imag,
                            vareps, Btransf_out, Ctransf_out, Dtransf_out, u, kappaxiq, q, S1, S2))

    # Order the roots

    Ssq3, Ssqm, Ssqp = np.sort([SsqI.real, SsqII.real, SsqIII.real])

    # Check for complex conjugates that passed the previous check, but allow for roots that are real and equal,
    # which occurs to numerical accuracy for close-to-aligned-spin cases
    if Ssq3.real == Ssqm.real:
        if Ssq3.imag != 0. or Ssqm.imag != 0.:
            warn(
                r"\bar{S}_3^2 and \bar{S}_-^2 have the same real parts, with values of: %e, %e"%(Ssq3, Ssqm), RuntimeWarning)
        elif Ssq3 == 0. and Ssqm == 0.:
            raise NonprecessingError(r"The system is nonprecessing: \bar{S}^2_3 = \bar(S}^2_- = 0") # This system is nonprecessing to numerical accuracy

    if Ssqm.real == Ssqp.real:
        if Ssqm.imag != 0. or Ssqp.imag != 0.:
            warn(
                r"\bar{S}_-^2 and \bar{S}_+^2 have the same real parts, with values of: %e, %e"%(Ssqm, Ssqp), RuntimeWarning)
        elif Ssqm == 0. and Ssqp == 0.:
            raise NonprecessingError(r"The system is nonprecessing: \bar{S}^2_- = \bar(S}^2_+ = 0") # This system is nonprecessing to numerical accuracy

    return Ssq3, Ssqm, Ssqp


def dkappaxiq_du(u, kappaxiq, q, S1, S2, xi, S0sq, imag_tol, root_tol, lin_tol, root_tol_raise_err, use_mpmath=False,
                 polyroots_maxsteps=50, polyroots_extraprec=50):
    """
    Compute (<S^2>_prec - S0^2)/(1 - q) [= dkappa_{xi q}/du] using Eq. (42) in Chatziioannou, Klein, Yunes, and Cornish PRD 95, 104004 (2017)
    [rewritten as bar{S}_+^2 + [E(m)/K(m) - 1](bar{S}_+^2 - bar{S}_3^2), where bar{S}_*^2 := (S_*^2 - S_0^2)/(1 - q)].
    Here kappa_{xi q} := (kappa - S_0^2/(2 L) - q xi/(1 + q))/(1 - q) - xi/4, and S_0 is the magnitude of the
    initial total spin (i.e., at u = u0).

    Inputs:

    - Required:

    u: 1/(2L), where L is the orbital angular momentum (i.e., the integration variable)
    kappaxiq: The rescaled variable defined above and in Eq. (14) of the paper <https://dcc.ligo.org/P2100029>, arXiv:2107.11902
    q: Mass ratio of the binary, with the < 1 convention
    S1, S2: Magnitudes of the dimensionful spins of the binary (with total mass = 1)
    xi: Conserved effective spin, defined in Eq. (2) of the paper <https://dcc.ligo.org/P2100029>, arXiv:2107.11902
    S0sq: Square of magnitude of the initial total spin (at u = u0, and with total mass = 1)
    imag_tol: Tolerance for the imaginary part of roots of the cubic
    root_tol: Tolerance for the error in the roots of the cubic (is an absolute tolerance for small values of the roots, as
              generally occurs during the beginning of the evolution, and is a fractional tolerance when the roots are large,
              which generally occurs towards the end of the evolution; also is reduced automatically when the roots are close
              together)
    lin_tol: Tolerance for errors on the r.h.s. of dkappa/du due to linearization in m (is a fractional tolerance when m is
             large, as generally occurs during the beginning of the evolution, and is an absolute tolerance when m is small,
             as generally occurs towards the end of the evolution)
    root_tol_raise_err: Whether to raise an error (True) or just print a warning (False) when the error in the
                        roots of the cubic is larger than root_tol

    - Optional:

    use_mpmath: If True, use mpmath functions instead of numpy functions to solve the cubic and evaluate the elliptic integrals
                (slower but more accurate), default: False
    polyroots_maxsteps: The maximum number of steps allowed in mpmath.polyroots if use_mpmath == True, default: 50
    polyroots_extraprec: The number of extra digits to use in mpmath.polyroots if use_mpmath == True, default: 50

    Output: (<S^2>_prec - S0^2)/(1 - q) [= dkappa_{xi q}/du]

    """

    # Compute the desired quantity for u > 0, otherwise set to 0

    if u > 0.:
        # Compute coefficients of equation

        vareps, Btransf, Ctransf, Dtransf = eq_coeffs(
            u, kappaxiq, q, S1, S2, xi, S0sq)

        # Compute the roots of the equation

        # First, check whether vareps is small enough that we can linearize <S^2>_pr in mm and compute
        # S_+^2 + S_-^2 to good accuracy by solving the quadratic one obtains by setting eps = 0 in the cubic
        # (so S_+^2 + S_-^2 = -Ctransf/Btransf).

        # Check if vareps is small enough that we can just solve the quadratic for S_+^2 + S_-^2
        # This is currently just implemented for the case where Btransf, Ctransf, Dtransf > 0
        # This is the case where there are issues in practice and also is the one that is easy to analyze

        only_quadratic = False

        # The second check is that vareps is small enough that the square root in z_eps_abs makes sense
        if min(Btransf, Ctransf, Dtransf) > 0 and 4.*Ctransf*vareps < Btransf*Btransf:
            z_eps_abs = abs(
                (-Btransf + (Btransf*Btransf - 4.*Ctransf*vareps)**0.5)/(2.*Ctransf))

            # Error on -Ctransf/Btransf expression for S_+^2 + S_-^2, E_\varepsilon^\text{sum} [Eq. (C5) in the paper <https://dcc.ligo.org/P2100029>, arXiv:2107.11902]
            sum_err_bound = z_eps_abs * \
                (abs(Ctransf*Ctransf - Btransf*Dtransf) + Ctransf*Dtransf*z_eps_abs) / \
                (Btransf*(Btransf - Ctransf*z_eps_abs))

            # Check that the additional inequality is satisfied
            if z_eps_abs*(Ctransf/Btransf + sum_err_bound) < 1:
                m_bound = z_eps_abs*(Ctransf/Btransf + sum_err_bound) / \
                    (1. - z_eps_abs*(Ctransf/Btransf + sum_err_bound))

                # This bounds S_+^2 - S_-^2 by abs(S_+^2 + S_-^2) in the first part, which is not sharp
                # We check that m_bound < 1 since E_ratio_term is only real for m_bound < 1
                # It should be possible to restrict computing E_ratio_term only for smaller m_bound values,
                # but it's not clear if this is necessary
                if m_bound < 1.:
                    if E_ratio_term(m_bound) <= 0.5*lin_tol*max(1., 1./(Ctransf/Btransf + sum_err_bound)) \
                       and sum_err_bound < lin_tol:
                        only_quadratic = True

        if only_quadratic:
            Ssq_pr = -0.5*Ctransf/Btransf

            # Check that Ssq_pr [= (bar{S}_+^2 + bar{S}_-^2)/2] gives an unbarred value that is nonnegative

            Ssq_pr_unbarred = (1. - q)*Ssq_pr + S0sq

            if Ssq_pr_unbarred < 0.:
                raise RuntimeError("The unbarred S_+^2 + S_-^2 value (which was all that was necessary to compute to this accuracy) is negative, "
                                   "indicating a problem with the evolution, with a value of %e"%(2.*Ssq_pr_unbarred))
        else:
            Ssq3, Ssqm, Ssqp = solve_cubic(
                [vareps, Btransf, Ctransf, Dtransf], u, kappaxiq, q, S1, S2, imag_tol, use_mpmath=use_mpmath, root_tol=root_tol,
                polyroots_maxsteps=polyroots_maxsteps, polyroots_extraprec=polyroots_extraprec)

            # Check that Ssqp and Ssqm give unbarred values that are nonnegative

            Ssqp_unbarred = (1. - q)*Ssqp + S0sq
            Ssqm_unbarred = (1. - q)*Ssqm + S0sq

            if Ssqp_unbarred < 0. or Ssqm_unbarred < 0.:
                raise RuntimeError("At least one of the unbarred S_+^2 and S_-^2 values are negative, indicating a problem with the evolution, "
                                   "with values of %e, %e"%(Ssqp_unbarred, Ssqm_unbarred))

            # Compute <S^2>_pr

            if Ssqp == Ssq3:
                raise RuntimeError(
                    "Degenerate nonprecessing case with all roots equal. Cannot continue with computation.")

            # When abs(Ssq3) is large, mm is small, and one can write Ssq_pr to good accuracy in a form that does not
            # involve Ssq3 (which has a large residual when put back into the equation when it is large, likely due to
            # catastrophic cancellation), but the residual is still smaller than Ssq3, so we take the value of Ssq3
            # output by the solver to give us an accurate enough value to use in checking the size of mm

            mm = (Ssqp - Ssqm)/(Ssqp - Ssq3)  # Eq. (25) in Chatziioannou, Klein, Yunes, and Cornish PRD 95, 104004 (2017)

            # Check that the equation is satisfied to good accuracy (by checking that the residual changes sign over an
            # interval given by root_tol) and compute Ssq_pr

            def cubic(Ssq):
                if use_mpmath:
                    return mp.polyval([vareps, Btransf, Ctransf, Dtransf], Ssq)
                else:
                    return ((vareps*Ssq + Btransf)*Ssq + Ctransf)*Ssq + Dtransf

            def cubic_deriv(Ssq):
                if use_mpmath:
                    return mp.polyval([3*vareps, 2*Btransf, Ctransf], Ssq)
                else:
                    return (3.*vareps*Ssq + 2.*Btransf)*Ssq + Ctransf

            # If the difference between Ssqp and Ssqm or Ssqm and Ssq3 is smaller than the difference set by root_tol,
            # shrink root_tol appropriately
            if Ssqp - Ssqm < 1.1*max(0.5*abs(Ssqp + Ssqm), 1.)*root_tol:
                root_tol_pm = 0.5*(Ssqp - Ssqm)/max(0.5*abs(Ssqp + Ssqm), 1.)
            else:
                root_tol_pm = root_tol

            if Ssqm - Ssq3 < 1.1*max(0.5*abs(Ssqm + Ssq3), 1.)*root_tol:
                root_tol_m3 = 0.5*(Ssqm - Ssq3)/max(0.5*abs(Ssqm + Ssq3), 1.)
            else:
                root_tol_m3 = root_tol

            root_tol_message_extra_bit_flag = False

            root_tol_orig = root_tol

            if root_tol_pm < root_tol and root_tol_pm <= root_tol_m3:
                root_tol = root_tol_pm

                root_tol_message_extra_bit_flag = True

                root_tol_message_extra_bit_insert = "S_+^2 - S_-^2"
            elif root_tol_m3 < root_tol:
                root_tol = root_tol_m3

                root_tol_message_extra_bit_flag = True

                root_tol_message_extra_bit_insert = "S_-^2 - S_3^2"

            if root_tol_message_extra_bit_flag:
                root_tol_message_extra_bit = " the internally modified (since %s is small) version "\
                                             "of"%root_tol_message_extra_bit_insert
            else:
                root_tol_message_extra_bit = ""

            Ssqp_resp, Ssqp_resm = cubic(
                Ssqp + max(abs(Ssqp), 1.)*root_tol), cubic(Ssqp - max(abs(Ssqp), 1.)*root_tol)
            Ssqm_resp, Ssqm_resm = cubic(
                Ssqm + max(abs(Ssqm), 1.)*root_tol), cubic(Ssqm - max(abs(Ssqm), 1.)*root_tol)

            # Check for special cases and if we can linearize

            if use_mpmath:
                Btransf_out, Ctransf_out, Dtransf_out = Btransf, Ctransf, Dtransf
            else:
                Btransf_out, Ctransf_out, Dtransf_out = Btransf[0], Ctransf[0], Dtransf[0]

            if Ssqp == Ssqm:
                # Since here we have a double root, we will check that the derivative of the cubic changes sign
                # If the difference between Ssqm and Ssq3 is smaller than the difference set by root_tol,
                # shrink root_tol appropriately
                if Ssqm - Ssq3 < 1.1*max(0.5*abs(Ssqm + Ssq3), 1.)*root_tol_orig:
                    root_tol_m0 = 0.5*(Ssqm - Ssq3)/max(0.5*abs(Ssqm + Ssq3), 1.)

                    root_tol_message_extra_bit = " the internally modified (since S_-^2 - S_3^2 is small) "\
                                                 "version of"
                else:
                    root_tol_m0 = root_tol_orig

                    root_tol_message_extra_bit = ""
                Ssqp_resp, Ssqp_resm = cubic_deriv(
                    Ssqp + max(abs(Ssqp), 1.)*root_tol_m0), cubic_deriv(Ssqp - max(abs(Ssqp), 1.)*root_tol_m0)

                if np.sign(Ssqp_resp) == np.sign(Ssqp_resm):
                    root_tol_message = "Error in S_+^2 is greater than%s the tolerance root_tol = %e--the "\
                                       "residual does not change sign over the interval determined by root_tol "\
                                       "(S_+^2 = S_-^2, so m = 0 and we do not need to consider S_3^2)\n"\
                                       "Printing the residuals for S^2 = S_+^2 + max(abs(S_+^2),1)*root_tol, S_+^2 "\
                                       "- max(abs(S_+^2),1)*root_tol, as well as the values of u, S_+^2, "\
                                       "q*(1 - q^2)*u^2, bar{B}, bar{C}, bar{D}\n%e %e %e %f %e %e %e %e" % \
                                       (root_tol_message_extra_bit, root_tol_m0, Ssqp_resp, Ssqp_resm, u, Ssqp, vareps,
                                        Btransf_out, Ctransf_out, Dtransf_out)

                    if root_tol_raise_err:
                        raise RuntimeError(root_tol_message)
                    else:
                        warn(root_tol_message, RuntimeWarning)

                Ssq_pr = Ssqp # Here m = 0
            elif Ssqm == Ssq3:
                # Since here we have a double root, we will check that the derivative of the cubic changes sign
                # If the difference between Ssqp and Ssqm is smaller than the difference set by root_tol,
                # shrink root_tol appropriately
                if Ssqp - Ssqm < 1.1*max(0.5*abs(Ssqp + Ssqm), 1.)*root_tol_orig:
                    root_tol_m1 = 0.5*(Ssqp - Ssqm)/max(0.5*abs(Ssqp + Ssqm), 1.)

                    root_tol_message_extra_bit = " the internally modified (since S_+^2 - S_-^2 is small) "\
                                                 "version of"
                else:
                    root_tol_m1 = root_tol_orig

                    root_tol_message_extra_bit = ""
                Ssqm_resp, Ssqm_resm = cubic_deriv(
                    Ssqm + max(abs(Ssqm), 1.)*root_tol_m1), cubic_deriv(Ssqm - max(abs(Ssqm), 1.)*root_tol_m1)

                if np.sign(Ssqm_resp) == np.sign(Ssqm_resm):
                    root_tol_message = "Error in S_-^2 is greater than%s the tolerance root_tol = %e--the "\
                                       "residual does not change sign over the interval determined by root_tol "\
                                       "(S_-^2 = S_3^2, so m = 1 and we do not need to consider S_+^2)\n"\
                                       "Printing the residuals for S^2 = S_-^2 + max(abs(S_-^2),1)*root_tol, S_-^2 "\
                                       "- max(abs(S_-^2),1)*root_tol, as well as the values of u, S_-^2, "\
                                       "q*(1 - q^2)*u^2, bar{B}, bar{C}, bar{D}\n%e %e %e %f %e %e %e %e" % \
                                       (root_tol_message_extra_bit, root_tol_m1, Ssqm_resp, Ssqm_resm, u, Ssqm, vareps,
                                        Btransf_out, Ctransf_out, Dtransf_out)

                    if root_tol_raise_err:
                        raise RuntimeError(root_tol_message)
                    else:
                        warn(root_tol_message, RuntimeWarning)

                Ssq_pr = Ssq3 # Here m = 1
            elif E_ratio_term(mm) <= lin_tol*max(1., 1./(Ssqp - Ssqm)): # Linearization check
                if np.sign(Ssqp_resp) == np.sign(Ssqp_resm) or np.sign(Ssqm_resp) == np.sign(Ssqm_resm):
                    root_tol_message = "Error in S_+^2 and/or S_-^2 is greater than%s the tolerance root_tol = %e--the "\
                                       "residual does not change sign over the interval determined by root_tol "\
                                       "(m is small enough to be linearized in, so we do not consider S_3^2)\n"\
                                       "Printing the residuals for S^2 = S_+^2 + max(abs(S_+^2),1)*root_tol, S_+^2 "\
                                       "- max(abs(S_+^2),1)*root_tol, S_-^2 + max(abs(S_-^2),1)*root_tol, "\
                                       "S_-^2 - max(abs(S_-^2),1)*root_tol, as well as the values of u, S_+^2, S_-^2, "\
                                       "q*(1 - q^2)*u^2, bar{B}, bar{C}, bar{D}\n%e %e %e %e %e %f %f %e %e %e %e" % \
                                       (root_tol_message_extra_bit, root_tol, Ssqp_resp, Ssqp_resm, Ssqm_resp, Ssqm_resm,
                                        u, Ssqp, Ssqm, vareps, Btransf_out, Ctransf_out, Dtransf_out)

                    if root_tol_raise_err:
                        raise RuntimeError(root_tol_message)
                    else:
                        warn(root_tol_message, RuntimeWarning)

                # Eq. (42) in Chatziioannou, Klein, Yunes, and Cornish PRD 95, 104004 (2017), with the elliptic functions linearized in m
                Ssq_pr = 0.5*(Ssqp + Ssqm)
            else:
                Ssq3_resp, Ssq3_resm = cubic(
                    Ssq3 + max(abs(Ssq3), 1.)*root_tol), cubic(Ssq3 - max(abs(Ssq3), 1.)*root_tol)

                if np.sign(Ssqp_resp) == np.sign(Ssqp_resm) or np.sign(Ssqm_resp) == np.sign(Ssqm_resm) \
                   or np.sign(Ssq3_resp) == np.sign(Ssq3_resm):
                    root_tol_message = "Error in S_+^2, S_-^2, and/or S_3^2 is greater than%s the tolerance root_tol = "\
                                       "%e--the residual does not change sign over the interval determined by "\
                                       "root_tol\nPrinting the residuals for S^2 = S_+^2 + max(abs(S_+^2),1)*root_tol, "\
                                       "S_+^2 - max(abs(S_+^2),1)*root_tol, S_-^2 + max(abs(S_-^2),1)*root_tol, "\
                                       "S_-^2 - max(abs(S_-^2),1)*root_tol, S_3^2 + max(abs(S_3^2),1)*root_tol, "\
                                       "S_3^2 - max(abs(S_-^2),1)*root_tol, as well as the values of u, S_+^2, S_-^2, "\
                                       "S_3^2, q*(1 - q^2)*u^2, bar{B}, bar{C}, bar{D}\n"\
                                       "%e %e %e %e %e %e %e %f %f %f %e %e %e %e" % \
                                       (root_tol_message_extra_bit, root_tol, Ssqp_resp, Ssqp_resm, Ssqm_resp, Ssqm_resm,
                                        Ssq3_resp, Ssq3_resm, u, Ssqp, Ssqm, Ssq3, vareps, Btransf_out, Ctransf_out,
                                        Dtransf_out)

                    if root_tol_raise_err:
                        raise RuntimeError(root_tol_message)
                    else:
                        warn(root_tol_message, RuntimeWarning)

                if use_mpmath:
                    ellip_ratio = mp.ellipe(mm)/mp.ellipk(mm)
                else:
                    ellip_ratio = ellipe(mm)/ellipk(mm)

                # Eq. (17) in the paper <https://dcc.ligo.org/P2100029>, arXiv:2107.11902
                Ssq_pr = Ssqp + (ellip_ratio - 1.)*(Ssqp - Ssq3)
    else:
        Ssq_pr = 0.

    return Ssq_pr


def kappaxiq_evol(q, S1, S2, S1L, xi, S0sq, u0, uf, du, method, imag_tol, root_tol, lin_tol, ode_atol, ode_rtol,
                  root_tol_raise_err, use_solve_ivp=False, solve_ivp_lsoda_min_step=0., use_mpmath=False,
                  polyroots_maxsteps=50, polyroots_extraprec=50, odefun_degree=None):
    """
    Compute kappaxiq when the binary has u = uf starting from u = u0

    Inputs:

    - Required:

    q: Mass ratio of the binary, with the < 1 convention
    S1, S2: Magnitudes of the dimensionful spins of the binary (with total mass = 1)
    S1L: Component of the dimensionful spin of hole 1 along the orbital angular momentum (with total mass = 1)
    xi: Conserved effective spin, defined in Eq. (2) of the paper <https://dcc.ligo.org/P2100029>, arXiv:2107.11902
    S0sq: Square of initial dimensionful total spin (in total mass = 1 units)
    u0: Initial value of u
    uf: Final value of u
    du: Step in u for the evolution
    method: method to use for scipy integration ("vode", "lsoda", "dopri5", or "dop853" for the ode interface and
            "RK45", "RK23", "DOP853", or "LSODA" for the solve_ivp interface); not used if use_mpmath == True
    imag_tol: Tolerance for the imaginary part of roots of the cubic
    root_tol: Tolerance for the error in the roots of the cubic (is an absolute tolerance for small values of the roots, as
              generally occurs during the beginning of the evolution, and is a fractional tolerance when the roots are large,
              which generally occurs towards the end of the evolution; also is reduced automatically when the roots are close
              together)
    lin_tol: Tolerance for errors on the r.h.s. of dkappa/du due to linearization in m (is a fractional tolerance when m is
             large, as generally occurs during the beginning of the evolution, and is an absolute tolerance when m is small,
             as generally occurs towards the end of the evolution)
    ode_atol: Absolute tolerance for the default scipy ode solver, also the single tolerance setting for the mpmath ode solver
              if it is used
    ode_rtol: Relative tolerance for the default scipy ode solver
    root_tol_raise_err: Whether to raise an error (True) or just print a warning (False) when the error in the roots of
                        the cubic is larger than root_tol

    - Optional:

    use_solve_ivp: If True, use scipy's solve_ivp function for the ODE integration, which gives access to a different set of
                   integrators and allows tighter tolerances for LSODA than the default ode function, but can hang for certain
                   parameters, default: False
    solve_ivp_lsoda_min_step: Minimum step to use in the solve_ivp LSODA integrator, default: 0.

    use_mpmath: If True, use mpmath functions and ode integration instead of numpy/scipy (slower but more accurate), default: False
    polyroots_maxsteps: The maximum number of steps allowed in mpmath.polyroots if use_mpmath == True, default: 50
    polyroots_extraprec: The number of extra digits to use in mpmath.polyroots if use_mpmath == True, default: 50
    odefun_degree: Degree of polynomial used by mpmath.odefun if use_mpmath == True; calculated automatically by
                   mpmath.odefun if None, default: None

    Output: kappaxiq at Lf
    """

    # Compute initial value of kappaxiq

    kappaxiq0 = S1L

    if use_solve_ivp:
        # Use scipy.integrate.solve_ivp

        # Set up kwargs for solver (necessary since one gets a warning if one passes min_step with an integrator other than LSODA)
        solve_ivp_kwargs = {'atol': ode_atol, 'rtol': ode_rtol}

        if method == "LSODA":
            solve_ivp_kwargs['min_step'] = solve_ivp_lsoda_min_step

        # Solve
        soln = solve_ivp(dkappaxiq_du, [u0, uf], [kappaxiq0], method=method, t_eval=[uf], args=(q, S1, S2, xi, S0sq, imag_tol,
                         root_tol, lin_tol, root_tol_raise_err), **solve_ivp_kwargs)

        # Return kappaxiq at u = uf

        return soln.y[0]
    elif use_mpmath:
        # Use mpmath.odefun

        # Solve (we have to write it in terms of ub := -u so that the final value of ub is greater than the initial value, as
        # required by mpmath.odefun)
        soln = mp.odefun(lambda ub, kappaxiq: -dkappaxiq_du(-ub, kappaxiq, q, S1, S2, xi, S0sq, imag_tol, root_tol,
                                                             lin_tol, root_tol_raise_err, use_mpmath=True,
                                                             polyroots_maxsteps=polyroots_maxsteps,
                                                             polyroots_extraprec=polyroots_extraprec),
                         -u0, kappaxiq0, ode_atol, odefun_degree)

        # Return kappaxiq at u = uf

        return soln(-uf)
    else:
        # Use scipy.integrate.ode

        # initialize the ODE solver
        solver = ode(dkappaxiq_du)
        solver.set_integrator(method, atol=ode_atol, rtol=ode_rtol)
        solver.set_initial_value(kappaxiq0, u0)
        solver.set_f_params(q, S1, S2, xi, S0sq, imag_tol,
                            root_tol, lin_tol, root_tol_raise_err)

        # Solve

        end = False

        while solver.successful() and not end:
            new_t = solver.t - du

            if new_t <= uf:
                new_t = uf
                end = True

            solver.integrate(new_t, step=1)

        if abs(solver.t - uf) > 0.75*du:
            raise RuntimeError(
                "Final u of %e is more than 0.75*du different from uf of %e" % (solver.t, uf))

        # Return kappaxiq at u = uf

        return solver.y


def prec_avg_tilt_comp_vec_inputs(m1, m2, chi1_vec, chi2_vec, fref, L0_vec=None, Lf=None, du=1.e-3, ode_atol_base=1.e-8,
                                  ode_atol_inplane_spin_scaling=True, ode_atol_floor=1.e-13, ode_rtol=-1, method="lsoda",
                                  use_solve_ivp=False, solve_ivp_lsoda_min_step=0., use_fallback=True, du_fallback=1.e-5,
                                  ode_atol_fallback=1.e-13, ode_rtol_fallback=-1, method_fallback="lsoda",
                                  use_mpmath_fallback=True, mpmath_dps=30, ode_tol_mpmath=1.e-15, odefun_degree=None,
                                  polyroots_maxsteps=50, polyroots_extraprec=50, LPNorder=2, LPNspins=True, imag_tol=1.e-6,
                                  root_tol=1.e-6, lin_tol=-1, lin_tol_fallback=-1, root_tol_raise_err=False,
                                  failure_mode='None', version='v1'):
    """
    Compute tilt angles at infinite separation or bounds and average value of tilt angles at finite separation
    from masses, spins, tilt angles, and in-plane spin angle at some frequency, taken to be small enough that
    the precession averaged spin evolution from Gerosa et al. PRD 92, 064016 (2015) [implemented following
    Chatziioannou, Klein, Yunes, and Cornish PRD 95, 104004 (2017)] is accurate, using the regularized expressions
    from <https://dcc.ligo.org/P2100029>, arXiv:2107.11902. This version takes vectors of spins as input (as opposed
    to spin angles) and allows one optionally to specify the initial orbital angular momentum instead of computing it
    from fref (and the binary's parameters).

    Inputs:

    - Required:

    m1, m2: Detector frame masses of the binary, in kg
    chi1_vec, chi2_vec: 3d numpy arrays with the dimensionless spins of the binary; these must be in coordinates where
                        the Newtonian orbital angular momentum is in the z-direction if L0_vec is None (the default),
                        otherwise all that is necessary is that chi1_vec, chi2_vec, and L0_vec are all in the same
                        coordinate system
    fref: Reference frequency, in Hz; ignored when L0_vec is not None

    - Optional initial and final L settings:

    L0_vec: 3d numpy array with the initial orbital angular momentum in units where the binary's total mass = 1 and in
            the same coordinates as the spins; this is computed from fref and the binary parameters if None is passed
            and in that case given in coordinates where the Newtonian orbital angular momentum is in the z-direction,
            default: None
    Lf: Final magnitude of orbital angular momentum in total mass = 1 units for output at finite separation,
        None gives the output at infinity, default: None

    - Optional settings for the evolution:

    du: Step in u for the evolution (external to the ODE integrator), default: 1.e-3
    ode_atol_base: Base value of absolute tolerance for ode solver, default: 1.e-8
    ode_atol_inplane_spin_scaling: Scale ode_atol_base by min(chi1*np.sin(tilt1),chi2*np.sin(tilt2)), with a floor of
                                   ode_atol_floor, to give good accuracy for small in-plane spin cases, default: True
    ode_atol_floor: Floor value for absolute tolerance for ode solver used when ode_atol_inplane_spin_scaling == True;
                    should be several orders of magnitude less than ode_atol_base, default: 1.e-13
    ode_rtol: Relative tolerance for ode solver, -1 sets this to 0.1*ode_atol, default: -1
    method: method to use for integration (either from scipy's ode interface, where "vode", "lsoda", "dopri5", and
            "dop853" have been tested, though "lsoda" is the recommended one; "RK45", "RK23", "DOP853", or "LSODA" for
            scipy's solve_ivp interface, if use_solve_ivp == True ("lsoda" and "dop853" are equivalent to "LSODA" and
            "DOP853" here, for compatibility); or "odefun" to use mpmath.odefun with its default Taylor method, which
            is all that is available as of v1.2.0; this uses the mpmath fallback settings and disables the fallback
            evolution--it is also quite slow), default: "lsoda"
    use_solve_ivp: If True, use scipy's solve_ivp interface and its integrators instead of the default ode interface; note
                   that the solve_ivp integrators apparently hang for certain difficult cases, default: False
    solve_ivp_lsoda_min_step: Minimum step to use in the solve_ivp LSODA integrator, default: 0.
    use_fallback: Whether to use the fallback integration if the primary integration fails (True) or not (False),
                  not available if method == "odefun", default: True
    du_fallback: Step in u for the evolution if the first evolution fails, default: 1.e-5
    ode_atol_fallback: Absolute tolerance for ode solver if the first evolution fails, default: 1.e-13
    ode_rtol_fallback: Relative tolerance for ode solver if the first evolution fails, -1 sets this to
                       0.1*ode_atol_fallback, default: -1
    method_fallback: method to use for integration if the first evolution fails (any of the scipy integrators available for method),
                     default: "lsoda"
    use_mpmath_fallback: Use the slow but accurate mpmath integration if the first fallback integration fails (True) or not (False),
                         default: True
    mpmath_dps: The number of decimal places to use for the mpmath computations, if the mpmath integration is triggered or
                method == "odefun", default: 30
    ode_tol_mpmath: Tolerance for the mpmath ode solver, only used if the mpmath fallback evolution is triggered or
                    method == "odefun", default: 1.e-15
    odefun_degree: Degree of polynomial used by mpmath.odefun if the mpmath fallback evolution is triggered or method == "odefun";
                   calculated automatically by mpmath.odefun if None, default: None
    polyroots_maxsteps: The maximum number of steps allowed in mpmath.polyroots if the mpmath fallback evolution is triggered or
                        method == "odefun", default: 50
    polyroots_extraprec: The number of extra digits to use in mpmath.polyroots if the mpmath fallback evolution is triggered or
                         method == "odefun", default: 50
    LPNorder: PN order to use in computing the intial orbital angular momentum from fref--choices are
              0, 1, 1.5, 2, or 2.5, where the half-integer orders are the spin contributions; ignored when L0_vec is
              not None, default: 2
    LPNspins: Whether to include the orbit-averaged spin contributions in the computation of the initial orbital
              angular momentum; ignored when L0_vec is not None, default: True
    imag_tol: Tolerance for the imaginary part of roots of the cubic, default: 1.e-6
    root_tol: Tolerance for the error in the roots of the cubic (is an absolute tolerance for small values of the roots, as
              generally occurs during the beginning of the evolution, and is a fractional tolerance when the roots are large,
              which generally occurs towards the end of the evolution; also is reduced automatically when the roots are close
              together), default: 1.e-6
    lin_tol: Tolerance for errors on the r.h.s. of dkappa/du due to linearization in m (is a fractional tolerance when m is
             large, as generally occurs during the beginning of the evolution, and is an absolute tolerance when m is small,
             as generally occurs towards the end of the evolution); -1 sets this to be the same as ode_atol, except for the
             mpmath evolution, where it is always set to the ode_tol_mpmath value, default: -1
    lin_tol_fallback: Tolerance for errors on the r.h.s. of dkappa/du due to linearization in m if the first evolution
                      fails, -1 sets this to be the same as ode_atol_fallback, default: -1
    root_tol_raise_err: Whether to raise an error (True) or just print a warning (False) when the error in the roots
                        of the cubic is larger than root_tol, default: False
    failure_mode: How the code behaves when the evolution fails (after the fallback evolution, if use_fallback is True, the
                  default). 'Error' means that the code raises a RuntimeError, while 'NAN' and 'None' mean that the code
                  raises a RuntimeWarning and returns np.nan or None for the output, respectively, default: 'None'
    version: Version of the calculation to use--currently only the initial version, v1, is available, default: "v1"

    Output: dictionary with entries 'tilt1_inf', 'tilt2_inf' for evolution to infinity and entries 'tilt1_sep_min',
            'tilt1_sep_max', 'tilt1_sep_avg', 'tilt2_sep_min', 'tilt2_sep_max', 'tilt2_sep_avg' for evolution to a
            finite separation (i.e., a finite orbital angular momentum)
    """

    # Set scipy integrators for which this has been tested

    if use_solve_ivp:
        tested_scipy_integrators = ["RK45", "RK23", "DOP853", "LSODA"]
    else:
        tested_scipy_integrators = ["vode", "lsoda", "dopri5", "dop853"]

    tested_scipy_integrator_string = ', '.join(tested_scipy_integrators)

    # Set known failure modes

    known_failure_modes = ['Error', 'NAN', 'None']

    # Check version

    if version != 'v1':
        raise ValueError("Only version v1 is available at the moment, while version = %s"%version)

    # Make "lsoda" and "dop853" uppercase when passed with use_solve_ivp == True, for compatibility

    methods_to_uppercase = ["lsoda", "dop853"]

    if use_solve_ivp:
        if method in methods_to_uppercase:
            method = method.upper()

        if method_fallback in methods_to_uppercase:
            method_fallback = method_fallback.upper()

    # Check that inputs are physical and otherwise make sense

    check_masses(m1, m2)

    if type(chi1_vec) is not np.ndarray or type(chi2_vec) is not np.ndarray or len(np.atleast_1d(chi1_vec)) != 3 or \
       len(np.atleast_1d(chi2_vec)) != 3:
        raise ValueError(
            "chi1_vec and chi2_vec must both be numpy arrays of length 3, while they are", chi1_vec, chi2_vec)

    if L0_vec is not None:
        if len(L0_vec) != 3:
            raise ValueError(
                "L0_vec must be either None or a numpy array of length 3, while it is", L0_vec)

    chi1 = norm(chi1_vec)
    chi2 = norm(chi2_vec)

    if chi1 == 0. or chi1 > 1. or chi2 == 0. or chi2 > 1.:
        raise ValueError(
            "The magnitudes of the spins must both be positive and at most 1, while they are chi1, chi2 = %f, %f" %
            (chi1, chi2))

    if fref is not None:
        check_fref(fref, m1, m2, "precession-averaged")
    elif L0_vec is None:
        raise ValueError("One of fref and L0_vec must not be None")

    if du <= 0. or du_fallback <= 0.:
        raise ValueError(
            "The step in u and its fallback value must both be positive, while they are du, du_fallback = %e, %e" %
            (du, du_fallback))

    if use_solve_ivp and not solve_ivp_avail:
        import scipy

        raise ValueError("solve_ivp is not available in the version of scipy being used (v%s), so use_solve_ivp == True "
                         "is not allowed."%scipy.__version__)

    if use_solve_ivp:
        interface_str = "solve_ivp"
    else:
        interface_str = "ode"

    if method not in np.append(tested_scipy_integrators, "odefun") or method_fallback not in tested_scipy_integrators:
        raise ValueError("The integration method must be odefun or one of the tested scipy integrators available with "
                         "the %s interface, %s, and the fallback integration method must be one of the tested scipy "
                         "integrators, while they are method, method_fallback = %s, %s" %
                         (interface_str, tested_scipy_integrator_string, method, method_fallback))

    if method == "odefun" and not mpmath_avail:
        raise ValueError("mpmath is not available, so odefun is not a possible choice for the integration method, "
                         "which must in this case be one of the tested scipy integrators: %s"%tested_scipy_integrator_string)

    if mpmath_dps <= 0 or not isinstance(mpmath_dps, int):
        raise ValueError("The number of decimal places to use for the mpmath computations must be a natural number, "
                         "while it is mpmath_dps = %f"%mpmath_dps)

    if odefun_degree is not None and (odefun_degree <= 0 or not isinstance(odefun_degree, int)):
        raise ValueError("The degree for the mpmath.odefun must be a natural number, while it is odefun_degree = %f"%odefun_degree)

    if polyroots_maxsteps <= 0 or polyroots_extraprec <= 0 or not isinstance(polyroots_maxsteps, int) \
            or not isinstance(polyroots_extraprec, int):
        raise ValueError("The maximum number of steps and extra precision for mpmath.polyroots must both be natural "
                         "numbers, while they are polyroots_maxsteps, polyroots_extraprec = %f, %f"%(polyroots_maxsteps,
                                                                                                     polyroots_extraprec))

    if L0_vec is None and LPNorder not in LPN_orders:
        raise ValueError("The PN order to use to compute the initial orbital angular momentum must be one of %s, "
                         "while it is LPNorder = %f" % (LPN_order_string, LPNorder))

    if failure_mode not in known_failure_modes:
        known_failure_modes_string = ', '.join(known_failure_modes)
        raise ValueError("The failure mode must be one of %s, while it is failure_mode = %s" %
                         (known_failure_modes_string, failure_mode))

    # Print a warning if both fref and L0_vec are not None

    if fref is not None and L0_vec is not None:
        warn("fref and L0_vec are both not None. In this case, the L0_vec value is used to initialize the "
             "evolution and fref is not used.", ValueWarning)

    # Set the failure output and the string corresponding to it
    # These default to None, since they are not used if failure_mode == 'Error'

    if failure_mode == 'NAN':
        failure_output = np.nan
        failure_output_string = 'np.nan'
    else:
        failure_output = None
        failure_output_string = 'None'

    # Set failure message for mpmath evolution

    def failure_message_mpmath(err):
        return "Encountered the exception\n\n%s\n\nusing the mpmath evolution and mpmath_dps = %d, " \
               "ode_tol_mpmath = %e, odefun_degree = %s, polyroots_maxsteps = %d, polyroots_extraprec " \
               "= %d. However, it may be possible to evolve this case successfully using more digits " \
               "and/or a tighter tolerance "%(format_error(err), mpmath_dps, ode_tol_mpmath, str(odefun_degree),
                                              polyroots_maxsteps, polyroots_extraprec)

    # Check that m1 > m2, and swap the objects if m1 < m2; raise an error if m1 == m2
    # This does not rotate the spins, since we do not need to track the orientation of the binary

    swap = False

    if m1 < m2:
        swap = True
        m1, m2 = m2, m1
        chi1, chi2 = chi2, chi1
        chi1_vec, chi2_vec = chi2_vec, chi1_vec

    eq_mass_check(m1, m2, Lf)

    # Save the masses in kg for later use

    m1_kg = m1
    m2_kg = m2

    # Compute derived quantities (scaling the masses so that their sum is 1)

    q = m2/m1

    M = m1 + m2

    m1 /= M
    m2 /= M

    S1 = m1*m1*chi1
    S2 = m2*m2*chi2

    S1_vec = m1*m1*chi1_vec
    S2_vec = m2*m2*chi2_vec

    # Compute initial orbital angular momentum (in total mass = 1 units)
    # We need to pass the masses in kg here to convert the frequency when computing from fref

    if L0_vec is not None:
        L0 = norm(L0_vec)
    else:
        if LPNspins:
            S1_vec_for_L = S1_vec
            S2_vec_for_L = S2_vec
        else:
            S1_vec_for_L = np.zeros(3)
            S2_vec_for_L = np.zeros(3)

        L0_vec = L_vec_from_masses_spins_and_f(
            m1_kg, m2_kg, fref, PNorder=LPNorder, S1_vec=S1_vec_for_L, S2_vec=S2_vec_for_L, no_mass_spin_input_check=True)

        L0 = norm(L0_vec)

    if L0 == 0.:
        raise ValueError("L0_vec must not be the zero vector")

    if Lf is not None and Lf <= L0:
        raise ValueError(
            "The magnitude of the final orbital angular momentum must be greater than that of the "
            "initial orbital angular momentum, %f, while it is Lf = %f" % (L0, Lf))

    # Setting ode_atol for the case when ode_atol_inplane_spin_scaling == True

    L0_hat = L0_vec/L0

    chi1_inplane = norm(chi1_vec - np.dot(chi1_vec, L0_hat)*L0_hat)
    chi2_inplane = norm(chi2_vec - np.dot(chi2_vec, L0_hat)*L0_hat)

    if ode_atol_inplane_spin_scaling:
        ode_atol = max(min(chi1_inplane, chi2_inplane)
                       * ode_atol_base, ode_atol_floor)
    else:
        ode_atol = ode_atol_base

    # Take care of the -1 cases for ode_rtol and lin_tol and their fallback versions

    if ode_rtol == -1:
        ode_rtol = 0.1*ode_atol

    if ode_rtol_fallback == -1:
        ode_rtol_fallback = 0.1*ode_atol_fallback

    if lin_tol == -1:
        lin_tol = ode_atol

    if lin_tol_fallback == -1:
        lin_tol_fallback = ode_atol_fallback

    # Check the tolerance settings

    if ode_atol_base <= 0. or ode_atol_floor <= 0. or ode_rtol <= 0. or ode_atol_fallback <= 0. \
       or ode_rtol_fallback <= 0. or ode_tol_mpmath <= 0.:
        raise ValueError(
            "The ode solver tolerance values and their fallback values must all be positive, while they are "
            "ode_atol_base, ode_atol_floor, ode_rtol, ode_atol_fallback, ode_rtol_fallback, "
            "ode_tol_mpmath = %e, %e, %e, %e, %e, %e" %
            (ode_atol_base, ode_atol_floor, ode_rtol, ode_atol_fallback, ode_rtol_fallback, ode_tol_mpmath))

    if ode_atol_inplane_spin_scaling and ode_atol_floor >= ode_atol_base:
        raise ValueError("ode_atol_floor should be less than ode_atol_base, likely by several orders of magnitude, "
                         "while they are ode_atol_floor, ode_atol_base = %e, %e" %
                         (ode_atol_floor, ode_atol_base))

    if imag_tol < 0. or root_tol < 0. or lin_tol < 0. or lin_tol_fallback < 0.:
        raise ValueError("The tolerance values must all be nonnegative, while they are imag_tol, root_tol, lin_tol, "
                         "lin_tol_fallback = %e, %e, %e, %e" % (imag_tol, root_tol, lin_tol, lin_tol_fallback))

    # Compute effective spin [Eq. (12) of Gerosa et al.]

    xi = np.dot((1. + q)*S1_vec + (1. + 1./q)*S2_vec, L0_hat)

    # Compute projection of S1_vec onto L0_vec

    S1L = np.dot(S1_vec, L0_hat)
    S2L = np.dot(S2_vec, L0_hat)

    # Return the initial tilts when spins are aligned and not in an unstable configuration

    cos_tilt1 = np.clip(S1L/S1, -1., 1.)
    cos_tilt2 = np.clip(S2L/S2, -1., 1.)

    tilt1 = np.arccos(cos_tilt1)
    tilt2 = np.arccos(cos_tilt2)

    if (tilt1 == np.pi and tilt2 in [0., np.pi]) or (tilt1 == 0. and tilt2 == 0.):
        return package_tilts(tilt1, tilt2, Lf, swap=swap)

    # Compute initial total dimensionful spin and angular momentum in total mass = 1 units

    S0sq = norm(S1_vec + S2_vec)**2

    # Compute kappa at infinity or Lf

    u0 = 0.5/L0

    if Lf is not None:
        uf = 0.5/Lf
    else:
        uf = 0.

    mpmath_tilts = False

    if method == "odefun":
            try:
                mpmath_tilts = True

                with mp.workdps(mpmath_dps):
                    k_evol = kappaxiq_evol(q, S1, S2, S1L, xi, S0sq, u0, uf, du, method, imag_tol, root_tol, ode_tol_mpmath,
                                           ode_tol_mpmath, ode_rtol, root_tol_raise_err, use_mpmath=True,
                                           polyroots_maxsteps=polyroots_maxsteps, polyroots_extraprec=polyroots_extraprec,
                                           odefun_degree=odefun_degree)

            except (RuntimeError, mp.ctx_base.StandardBaseContext.NoConvergence) as err:
                return evolution_error_handling(failure_mode, failure_message_mpmath(err), failure_output, failure_output_string, Lf, swap)

            except NonprecessingError as ne:
                return package_tilts(tilt1, tilt2, Lf, swap)
    else:
        try:
            k_evol = kappaxiq_evol(q, S1, S2, S1L, xi, S0sq, u0, uf, du, method, imag_tol, root_tol, lin_tol,
                                   ode_atol, ode_rtol, root_tol_raise_err, use_solve_ivp=use_solve_ivp,
                                   solve_ivp_lsoda_min_step=solve_ivp_lsoda_min_step)

        except RuntimeError as re:
            if use_fallback:
                warn("Encountered the exception\n\n%s\n\nusing method %s and du = %e, ode_atol = %e, ode_rtol = %e, "
                     "lin_tol = %e, trying with %s with du = %e, ode_atol = %e, ode_rtol = %e, lin_tol = %e" %
                     (format_error(re), method, du, ode_atol, ode_rtol, lin_tol, method_fallback, du_fallback,
                      ode_atol_fallback, ode_rtol_fallback, lin_tol_fallback), RuntimeWarning)

                try:
                    k_evol = kappaxiq_evol(q, S1, S2, S1L, xi, S0sq, u0, uf, du_fallback, method_fallback, imag_tol,
                                           root_tol, lin_tol_fallback, ode_atol_fallback, ode_rtol_fallback,
                                           root_tol_raise_err, solve_ivp_lsoda_min_step=solve_ivp_lsoda_min_step)

                except RuntimeError as re2:
                    if use_mpmath_fallback and mpmath_avail:
                        warn("Encountered the exception\n\n%s\n\nin the first fallback evolution. Now trying with the mpmath "
                             "fallback with mpmath_dps = %d, ode_tol_mpmath = %e, odefun_degree = %s, "
                             "polyroots_maxsteps = %d, polyroots_extraprec = %d"%(format_error(re2), mpmath_dps, ode_tol_mpmath, str(odefun_degree),
                                                                                  polyroots_maxsteps, polyroots_extraprec), RuntimeWarning)

                        try:
                            mpmath_tilts = True

                            with mp.workdps(mpmath_dps):
                                k_evol = kappaxiq_evol(q, S1, S2, S1L, xi, S0sq, u0, uf, du, method, imag_tol, root_tol, ode_tol_mpmath,
                                                       ode_tol_mpmath, ode_rtol, root_tol_raise_err, use_mpmath=True,
                                                       polyroots_maxsteps=polyroots_maxsteps, polyroots_extraprec=polyroots_extraprec,
                                                       odefun_degree=odefun_degree)

                        except (RuntimeError, mp.ctx_base.StandardBaseContext.NoConvergence) as err:
                            return evolution_error_handling(failure_mode, failure_message_mpmath(err), failure_output, failure_output_string, Lf, swap)

                        except NonprecessingError as ne:
                            return package_tilts(tilt1, tilt2, Lf, swap)
                    else:
                        if not use_mpmath_fallback:
                            mpmath_bit = "it is disabled"
                        else:
                            mpmath_bit = "mpmath is not available"

                        failure_message = "Fallback evolution failed with the exception\n\n%s\n\nThe mpmath fallback evolution might work in this " \
                                          "case, but %s. Alternatively, you can try adjusting the tolerance settings for the evolution, but this " \
                                          "case may not be possible to evolve accurately with the scipy integration." % (format_error(re2), mpmath_bit)

                        return evolution_error_handling(failure_mode, failure_message, failure_output, failure_output_string, Lf, swap)

                except NonprecessingError as ne:
                    return package_tilts(tilt1, tilt2, Lf, swap)

            else:
                failure_message = "Evolution failed with the exception\n\n%s\n\nYou can try adjusting the tolerance settings " \
                                  "for the evolution, but this case may not be possible to evolve accurately with the current " \
                                  "version of the precession-averaged evolution."%format_error(re)

                return evolution_error_handling(failure_mode, failure_message, failure_output, failure_output_string, Lf, swap)

        except NonprecessingError as ne:
            return package_tilts(tilt1, tilt2, Lf, swap)

    # Compute the desired output. We have already accounted for the case where at least one spin is zero,
    # so there are no concerns about division by zero here.

    if Lf is None:
        # Compute cosines of tilt angles at infinity [Eqs. (15) in the paper <https://dcc.ligo.org/P2100029>, arXiv:2107.11902]

        cos_tilt1inf = k_evol/S1
        cos_tilt2inf = q*(xi/(1. + q) - k_evol)/S2

        # Clip to avoid issues with taking the arccos if there are rounding issues

        cos_tilt1inf = np.clip(cos_tilt1inf, -1., 1.)
        cos_tilt2inf = np.clip(cos_tilt2inf, -1., 1.)

        # Return tilt angles at infinity

        if mpmath_tilts:
            tilt1inf = float(mp.acos(cos_tilt1inf))
            tilt2inf = float(mp.acos(cos_tilt2inf))
        else:
            tilt1inf = np.arccos(cos_tilt1inf[0])
            tilt2inf = np.arccos(cos_tilt2inf[0])

        return package_tilts(tilt1inf, tilt2inf, None, swap)
    else:
        vareps, Btransf, Ctransf, Dtransf = eq_coeffs(
            uf, k_evol, q, S1, S2, xi, S0sq)

        _, Ssqm, Ssqp = solve_cubic(
            [vareps, Btransf, Ctransf, Dtransf], uf, k_evol, q, S1, S2, imag_tol, use_mpmath=mpmath_tilts,
            root_tol=root_tol, polyroots_maxsteps=polyroots_maxsteps, polyroots_extraprec=polyroots_extraprec)

        if mpmath_tilts:
            lin_tol_use = ode_tol_mpmath
        else:
            lin_tol_use = lin_tol

        Ssqavg = dkappaxiq_du(uf, k_evol, q, S1, S2, xi, S0sq, imag_tol, root_tol, lin_tol_use, root_tol_raise_err,
                              use_mpmath=mpmath_tilts, polyroots_maxsteps=polyroots_maxsteps,
                              polyroots_extraprec=polyroots_extraprec)

        # Compute cosines of tilt angles at separation [Eqs. (23) in paper <https://dcc.ligo.org/P2100029>, arXiv:2107.11902]

        cos_tilt1sep_Sm = (k_evol - uf*Ssqm)/S1
        cos_tilt2sep_Sm = q*(xi/(1. + q) - S1*cos_tilt1sep_Sm)/S2

        cos_tilt1sep_Sp = (k_evol - uf*Ssqp)/S1
        cos_tilt2sep_Sp = q*(xi/(1. + q) - S1*cos_tilt1sep_Sp)/S2

        cos_tilt1sep_Savg = (k_evol - uf*Ssqavg)/S1
        cos_tilt2sep_Savg = q*(xi/(1. + q) - S1*cos_tilt1sep_Savg)/S2

        # Clip to avoid issues with taking the arccos if there are rounding issues

        cos_tilt1sep_Sm = np.clip(cos_tilt1sep_Sm, -1., 1.)
        cos_tilt2sep_Sp = np.clip(cos_tilt2sep_Sp, -1., 1.)

        cos_tilt1sep_Sm = np.clip(cos_tilt1sep_Sm, -1., 1.)
        cos_tilt2sep_Sp = np.clip(cos_tilt2sep_Sp, -1., 1.)

        cos_tilt1sep_Savg = np.clip(cos_tilt1sep_Savg, -1., 1.)
        cos_tilt2sep_Savg = np.clip(cos_tilt2sep_Savg, -1., 1.)

        # Return tilt angles at separation

        if swap:
            cos_tilt1sep_Sm, cos_tilt2sep_Sm, cos_tilt1sep_Sp, cos_tilt2sep_Sp = cos_tilt2sep_Sp, cos_tilt1sep_Sp, \
                                                                                 cos_tilt2sep_Sm, cos_tilt1sep_Sm
            cos_tilt1sep_Savg, cos_tilt2sep_Savg = cos_tilt2sep_Savg, cos_tilt1sep_Savg

        if mpmath_tilts:
            tilt1sepmin = float(mp.acos(cos_tilt1sep_Sm))
            tilt1sepmax = float(mp.acos(cos_tilt1sep_Sp))
            tilt1sepavg = float(mp.acos(cos_tilt1sep_Savg))

            tilt2sepmin = float(mp.acos(cos_tilt2sep_Sp))
            tilt2sepmax = float(mp.acos(cos_tilt2sep_Sm))
            tilt2sepavg = float(mp.acos(cos_tilt2sep_Savg))
        else:
            tilt1sepmin = np.arccos(cos_tilt1sep_Sm[0])
            tilt1sepmax = np.arccos(cos_tilt1sep_Sp[0])
            tilt1sepavg = np.arccos(cos_tilt1sep_Savg[0])

            tilt2sepmin = np.arccos(cos_tilt2sep_Sp[0])
            tilt2sepmax = np.arccos(cos_tilt2sep_Sm[0])
            tilt2sepavg = np.arccos(cos_tilt2sep_Savg[0])

        return {'tilt1_sep_min': tilt1sepmin, 'tilt1_sep_max': tilt1sepmax, 'tilt1_sep_avg': tilt1sepavg,
                'tilt2_sep_min': tilt2sepmin, 'tilt2_sep_max': tilt2sepmax, 'tilt2_sep_avg': tilt2sepavg}


def prec_avg_tilt_comp(m1, m2, chi1, chi2, tilt1, tilt2, phi12, fref, Lf=None, **kwargs):
    """
    Compute tilt angles at infinite separation or bounds and average value of tilt angles at finite separation from
    masses, spins, tilt angles, and in-plane spin angle at some frequency, taken to be small enough that the precession
    averaged spin evolution from Gerosa et al. PRD 92, 064016 (2015) [implemented following Chatziioannou, Klein, Yunes,
    and Cornish PRD 95, 104004 (2017)] is accurate, using the regularized expressions. This takes in spin angles and a
    reference frequency to initialize the evolution.

    Inputs:

    - Required:

    m1, m2: Detector frame masses of the binary, in kg
    chi1, chi2: Dimensionless spin magnitudes of the binary
    tilt1, tilt2: Tilt angles of the binary's spins (w.r.t. the Newtonian orbital angular momentum) at fref
    phi12: Angle between the in-plane components of the spins at fref
    fref: Reference frequency, in Hz

    - Optional settings [passed to prec_avg_tilt_comp_vec_inputs(); all but Lf passed through kwargs]:

    Lf: Final magnitude of orbital angular momentum in total mass = 1 units for output at finite separation
        None gives the output at infinity, default: None

    du: Step in u for the evolution (external to the ODE integrator), default: 1.e-3
    ode_atol_base: Base value of absolute tolerance for ode solver, default: 1.e-8
    ode_atol_inplane_spin_scaling: Scale ode_atol_base by min(chi1*np.sin(tilt1),chi2*np.sin(tilt2)), with a floor of
                                   ode_atol_floor, to give good accuracy for small in-plane spin cases, default: True
    ode_atol_floor: Floor value for absolute tolerance for ode solver used when ode_atol_inplane_spin_scaling == True;
                    should be several orders of magnitude less than ode_atol_base, default: 1.e-13
    ode_rtol: Relative tolerance for ode solver, -1 sets this to 0.1*ode_atol, default: -1
    method: method to use for integration (either from scipy's ode interface, where "vode", "lsoda", "dopri5", and
            "dop853" have been tested, though "lsoda" is the recommended one; "RK45", "RK23", "DOP853", or "LSODA" for
            scipy's solve_ivp interface, if use_solve_ivp == True ("lsoda" and "dop853" are equivalent to "LSODA" and
            "DOP853" here, for compatibility); or "odefun" to use mpmath.odefun with its default Taylor method, which
            is all that is available as of v1.2.0; this uses the mpmath fallback settings and disables the fallback
            evolution--it is also quite slow), default: "lsoda"
    use_solve_ivp: If True, use scipy's solve_ivp interface and its integrators instead of the default ode interface; note
                   that the solve_ivp integrators apparently hang for certain difficult cases, default: False
    solve_ivp_lsoda_min_step: Minimum step to use in the solve_ivp LSODA integrator, default: 0.
    use_fallback: Whether to use the fallback integration if the primary integration fails (True) or not (False),
                  not available if method == "odefun", default: True
    du_fallback: Step in u for the evolution if the first evolution fails, default: 1.e-5
    ode_atol_fallback: Absolute tolerance for ode solver if the first evolution fails, default: 1.e-13
    ode_rtol_fallback: Relative tolerance for ode solver if the first evolution fails, -1 sets this to
                       0.1*ode_atol_fallback, default: -1
    method_fallback: method to use for integration if the first evolution fails (any of the scipy integrators available for method),
                     default: "lsoda"
    use_mpmath_fallback: Use the slow but accurate mpmath integration if the first fallback integration fails (True) or not (False),
                         default: True
    mpmath_dps: The number of decimal places to use for the mpmath computations, if the mpmath integration is triggered or
                method == "odefun", default: 30
    ode_tol_mpmath: Tolerance for the mpmath ode solver, only used if the mpmath fallback evolution is triggered or
                    method == "odefun", default: 1.e-15
    odefun_degree: Degree of polynomial used by mpmath.odefun if the mpmath fallback evolution is triggered or method == "odefun";
                   calculated automatically by mpmath.odefun if None, default: None
    polyroots_maxsteps: The maximum number of steps allowed in mpmath.polyroots if the mpmath fallback evolution is triggered or
                        method == "odefun", default: 50
    polyroots_extraprec: The number of extra digits to use in mpmath.polyroots if the mpmath fallback evolution is triggered or
                         method == "odefun", default: 50
    LPNorder: PN order to use in computing the intial orbital angular momentum from fref--choices are
              0, 1, 1.5, 2, or 2.5, where the half-integer orders are the spin contributions, default: 2
    LPNspins: Whether to include the orbit-averaged spin contributions in the computation of the initial orbital
              angular momentum, default: True
    imag_tol: Tolerance for the imaginary part of roots of the cubic, default: 1.e-6
    root_tol: Tolerance for the error in the roots of the cubic (is an absolute tolerance for small values of the roots, as
              generally occurs during the beginning of the evolution, and is a fractional tolerance when the roots are large,
              which generally occurs towards the end of the evolution; also is reduced automatically when the roots are close
              together), default: 1.e-6
    lin_tol: Tolerance for errors on the r.h.s. of dkappa/du due to linearization in m (is a fractional tolerance when m is
             large, as generally occurs during the beginning of the evolution, and is an absolute tolerance when m is small,
             as generally occurs towards the end of the evolution); -1 sets this to be the same as ode_atol, except for the
             mpmath evolution, where it is always set to the ode_tol_mpmath value, default: -1
    lin_tol_fallback: Tolerance for errors on the r.h.s. of dkappa/du due to linearization in m if the first evolution
                      fails, -1 sets this to be the same as ode_atol_fallback, default: -1
    root_tol_raise_err: Whether to raise an error (True) or just print a warning (False) when the error in the roots
                        of the cubic is larger than root_tol, default: False
    failure_mode: How the code behaves when the evolution fails (after the fallback evolution, if use_fallback is True, the
                  default). 'Error' means that the code raises a RuntimeError, while 'NAN' and 'None' mean that the code
                  raises a RuntimeWarning and returns np.nan or None for the output, respectively, default: 'None'
    version: Version of the calculation to use--currently only the initial version, v1, is available, default: "v1"

    Output: dictionary with entries 'tilt1_inf', 'tilt2_inf' for evolution to infinity and entries 'tilt1_sep_min',
            'tilt1_sep_max', 'tilt1_sep_avg', 'tilt2_sep_min', 'tilt2_sep_max', 'tilt2_sep_avg' for evolution to a
            finite separation (i.e., a finite orbital angular momentum)
    """

    # Return initial tilts for single-spin or nonspinning cases

    if chi1 == 0. or chi2 == 0.:
        return package_tilts(tilt1, tilt2, Lf, False)

    # Check that spins and tilt angles are physical/make sense [other checks carried out in prec_avg_tilt_comp_vec_inputs()]

    check_spin_mags(chi1, chi2)

    check_tilts(tilt1, tilt2)

    # Convert spin magnitudes and angles to vectors

    chi1_vec = chi1*np.array([np.sin(tilt1), 0., np.cos(tilt1)])

    chi2_vec = chi2*np.array([np.sin(tilt2)*np.cos(phi12),
                              np.sin(tilt2)*np.sin(phi12), np.cos(tilt2)])

    # Pass to prec_avg_tilt_comp_vec_inputs() to do calculation

    return prec_avg_tilt_comp_vec_inputs(m1, m2, chi1_vec, chi2_vec, fref, Lf=Lf, **kwargs)
