import numpy as np
import lalsimulation as lalsim
from lal import MSUN_SI, MTSUN_SI

def spin_evolution(q, chiA0, chiB0, omega0, approximant='SpinTaylorT4',
        dt=0.1, spinO=6, phaseO=7):
    """
    Wrapper for PN spin and dynamics evolution in LAL.
    Inputs:
        - q:              Mass ratio (q>=1)
        - chiA0:          Dimensionless spin of BhA at initial freq.
        - chiB0:          Dimensionless spin of BhB at initial freq.
        - omega0:         Initial orbital frequency in dimensionless units.
        - approximant:    'SpinTaylorT1/T4/T5'. Default: SpinTaylorT4.
        - dt:             Dimensionless step time for evolution. Default: 0.1
        - spinO:          Twice PN order of spin effects. Default: 7.
        - phaseO:         Twice PN order in phase. Default: 7.

    Outputs (all are time series):
        - Omega:          Dimensionless orbital frequency.
        - Phi:            Orbital phase (radians)
        - ChiA:           Dimensionless spin of BhA
        - ChiB:           Dimensionless spin of BhB
        - LNhat:          Orbital angular momentum direction
        - E1:             Orbital plane basis vector

    The frame is defined at the initial frequency omega0, as follows: \n
        - z-axis is set by the orbital angular momentum direction.
        - x-axis is the separation vector from the lighter BH to the heavier BH.
        - y-axis completes the triad by right-hand rule. \n
        All quantities are defined in this fixed frame, including initial spins,
        returned spins, other vectors like LNhat, etc.
    """

    approxTag = lalsim.SimInspiralGetApproximantFromString(approximant)

    # Total mass in solar masses
    M = 100     # This does not affect the returned values as they are
                # dimensionless

    # time step and initial GW freq in SI units
    MT = M*MTSUN_SI
    deltaT = dt*MT
    fStart = omega0/np.pi/MT

    # component masses of the binary
    m1_SI =  M*MSUN_SI*q/(1.+q)
    m2_SI =  M*MSUN_SI/(1.+q)

    # spins at fStart
    s1x, s1y, s1z = chiA0
    s2x, s2y, s2z = chiB0

    # integrate as far forward as possible
    fEnd = 0

    # initial value of orbital angular momentum unit vector, i.e at fStart
    lnhatx, lnhaty, lnhatz = 0,0,1

    # initial value of orbital plane basis vector, i.e at fStart
    e1x, e1y, e1z = 1, 0, 0

    # tidal deformability parameters
    lambda1, lambda2 = 0, 0

    # spin-induced quadrupole moments
    quadparam1, quadparam2 = 1, 1

    # twice PN order of tidal effects
    tideO = 0

    # include some known L-S terms
    lscorr = 1

    # evolve spins and collect data into a nice format. The input values start
    # with lower-case letters, while output values start with upper-case
    # letters.
    V, Phi, S1x, S1y, S1z, S2x, S2y, S2z, LNhatx, LNhaty, LNhatz, \
        E1x, E1y, E1z = lalsim.SimInspiralSpinTaylorPNEvolveOrbit(deltaT, \
        m1_SI, m2_SI, fStart, fEnd, s1x, s1y, s1z, s2x, s2y, s2z, \
        lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2, \
        quadparam1, quadparam2, spinO, tideO, phaseO, lscorr, approxTag)
    V = np.array(V.data.data)
    Phi = np.array(Phi.data.data)
    ChiA = np.array([S1x.data.data, S1y.data.data, S1z.data.data]).T
    ChiB = np.array([S2x.data.data, S2y.data.data, S2z.data.data]).T
    LNhat = np.array([LNhatx.data.data, LNhaty.data.data, LNhatz.data.data]).T
    E1 = np.array([E1x.data.data, E1y.data.data, E1z.data.data]).T

    # orbital frequency
    Omega = V**3
    return Omega, Phi, ChiA, ChiB, LNhat
