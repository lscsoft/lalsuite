# Copyright (C) 2011  Nickolas Fotopoulos, Stephen Privitera
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

from __future__ import division

# standard
from math import sqrt

import numpy
from numpy.random.mtrand import uniform
from scipy.optimize import fsolve

# local
try:
    from glue.iterutils import choices
except ImportError:
    raise ImportError("The sbank subpackage of lalinspiral depends on the glue and pylal packages.")

from lal import PI, MTSUN_SI
#
# functions for converting between m1-m2 and tau0-tau3 coords
#

def A0(flow):
    '''
    A0 is a conversion factor between tau0 and mchirp, defined in eqn
    B3 of the appendix of arxiv.org/abs/0706.4437
    '''
    return 5./(256*(PI * flow)**(8./3)) # eqn B3

def A3(flow):
    '''
    A3 is a conversion factor between tau3 and tau0*M, defined in eqn
    B3 of the appendix of arxiv.org/abs/0706.4437
    '''
    return PI/(8*(PI*flow)**(5./3)) # eqn B3

def tau0tau3_to_m1m2(tau0, tau3, flow):
    '''
    Convert tau0-tau3 coordinates to m1-m2 coordinates.
    '''
    # compute mtotal in seconds
    mtotal = 5. / (32 * PI * PI * flow) * (tau3 / tau0)
    eta = ( A0(flow) / tau0 ) * mtotal**(-5./3)

    # convert to solar masses
    mtotal /= MTSUN_SI
    m1 = mtotal * (0.5 + (0.25 - eta)**0.5)

    return m1, mtotal - m1

def m1m2_to_tau0tau3(m1, m2, flow):
    '''
    Convert m1-m2 coordinates to tau0-tau3 coordinates.
    '''
    # convert solar masses to seconds
    m1_s = m1 * MTSUN_SI
    m2_s = m2 * MTSUN_SI

    # compute mass-dependent terms
    mtotal = m1_s + m2_s
    eta = m1_s * m2_s / (mtotal * mtotal)

    # eqn B4 for tau0,tau3
    tau0 = ( A0(flow) / eta ) * mtotal**(-5./3)
    tau3 = ( A3(flow) / eta ) * mtotal**(-2./3)

    return tau0, tau3


#
# functions for handling m1-m2 constraints
#

def set_default_constraints(constraints):
    '''
    Check that the constraints dict does not containt unknown
    constraints. Check that required constraints are given. To define
    a 2D parameter space, one needs at least three and possibly four
    points of intersection. This function preps the constraint dict so
    that it defines a bounded region (if possible) or crashes (if
    not). By convention, we require mass1 limits to be set.
    '''
    # complain about unknown constraints; don't silently ignore!
    known_constraints = ["mass1","mass2","mratio","mtotal","mchirp","spin1","spin2"]
    unknown_constraints = [k for k in constraints.keys() if k not in known_constraints]
    if len(unknown_constraints):
        raise ValueError("unknown constraints %s" % ', '.join(unknown_constraints))

    # component mass1 required
    mass1_min, mass1_max = constraints.setdefault('mass1', (None, None))
    if mass1_min is None or mass1_max is None:
        raise ValueError("Must specify both minimum and maximum mass1.")

    # component mass2 can be specified independently
    # if not specified, symmetry will be assumed
    mass2_min, mass2_max = constraints.setdefault('mass2', (None, None))
    if mass2_min is None:
        mass2_min = mass1_min
    if mass2_max is None:
        mass2_max = mass1_max
    constraints['mass2'] = (mass2_min, mass2_max)

    # mtotal can be given or inferred from component mass limits
    mtotal_min, mtotal_max = constraints.setdefault('mtotal', (None, None))
    if mtotal_min is None or mtotal_min < mass1_min + mass2_min:
        mtotal_min = mass1_min + mass2_min
    if mtotal_max is None or mass1_max + mass2_max < mtotal_max:
        mtotal_max = mass1_max + mass2_max

    # mratio can be given or inferred from component mass limits
    qmin, qmax = constraints.setdefault('mratio', (None, None))
    if qmin is None or qmin < mass1_min/mass2_max:
        qmin = max(mass1_min/mass2_max, 1) # q = m1/m2 > 1 by convention
    if qmin < 1:
        raise ValueError("We use the convention that q = m1/m2 > 1.")
    if qmax is None or min(mass1_max, mtotal_max - mass2_min) / mass2_min < qmax:
        qmax = min(mass1_max, mtotal_max - mass2_min) / mass2_min # q = m1/m2 by convention
    constraints['mratio'] = (qmin, qmax)

    # mchirp can be given or inferred from component mass limits
    mcmin, mcmax = constraints.setdefault('mchirp', (None, None))
    if mcmin is None:
        mcmin = min(m1m2_to_mchirp(m1, m2) for m1 in (mass1_min, mass1_max) for m2 in (mass2_min, mass2_max))
    if mcmax is None:
        mcmax = max(m1m2_to_mchirp(m1, m2) for m1 in (mass1_min, mass1_max) for m2 in (mass2_min, mass2_max))
    constraints['mchirp'] = (mcmin, mcmax)

    # update mtotal constraints based on mchirp cuts one can show that
    # for fixed mchirp, dM/dq > 0 provided q>1 so just set q=qmin. on
    # the other hand, mchirp constraints do not influence bounds on q
    if mtotal_min < mcmin * ((1+qmin)**2/qmin)**(3./5):
        mtotal_min = mcmin * ((1+qmin)**2/qmin)**(3./5)
    if mtotal_max > mcmax * ((1+qmax)**2/qmax)**(3./5):
        mtotal_max = mcmax * ((1+qmax)**2/qmax)**(3./5)
    constraints['mtotal'] = (mtotal_min, mtotal_max)

    return constraints

def m1m2_to_mratio(m1,m2):
    return m1/m2

def m1m2_to_m1(m1,m2):
    return m1

def m1m2_to_m2(m1,m2):
    return m2

def m1m2_to_mtotal(m1,m2):
    return m1+m2

def m1m2_to_mchirp(m1,m2):
    return (m1 * m1 * m1 * m2 * m2 * m2 / (m1 + m2))**0.2

m1m2_mappings = { \
    "mratio":m1m2_to_mratio,
    "mass1":m1m2_to_m1,
    "mass2":m1m2_to_m2,
    "mtotal":m1m2_to_mtotal,
    "mchirp":m1m2_to_mchirp}

def allowed_m1m2(m1, m2, m1m2_constraints, tol=1e-10):
    '''
    Return those values of m1,m2 from the input arrays
    that are consistent with the specified constraints.
    '''
    for param in m1m2_constraints.keys():

        # get the constraints
        mn,mx = m1m2_constraints[param]

        # don't bother if there aren't any constraints
        if mn is None and mx is None:
            continue

        # otherwise compute surviving pairs
        params = m1m2_mappings[param](m1,m2)
        if mn is None:
            include = (params <= mx + tol)
        elif mx is None:
            include = (mn - tol <= params)
        else:
            include = (mn - tol <= params)*(params <= mx + tol)

        m1 = m1[include]
        m2 = m2[include]

    return m1,m2


def mchirpm1_to_m2(mc, m1, tol=1e-6):
    # solve cubic for m2
    p = -mc**5/m1**3
    q = p*m1
    rts = numpy.roots([1, 0, p, q])

    # remove complex and negative roots
    rts = [r.real for r in rts if abs(r.imag) < tol and r.real > 0]

    if len(rts) == 0:
        m2 = float('nan')
    elif len(rts) == 1:
        m2 = rts[0]
    else:
        raise ValueError("No unique real solution for m2 found for mchirp=%f and m1=%f"%(mc, m1))

    return m2

def tau0tau3_bound(flow, **constraints):
    '''
    For a given set of constraints on the m1-m2 parameter space,
    this function returns the corners of a box which bounds the
    same region (tightly) above and below in tau0-tau3 space.

    Supported constraints: mass1, mass2, mtotal, mratio, mchirp
    '''
    # ensure we can at least bound the region in m1-m2
    constraints = set_default_constraints( constraints )

    # controls delta-m resolution
    npts = 1e4

    # draw constant component mass lines
    m1min, m1max = constraints['mass1']
    m2min, m2max = constraints['mass2']
    m1m2_points = [(m1, m2min) for m1 in numpy.linspace(m1min, m1max, npts)]
    m1m2_points += [(m1, m2max) for m1 in numpy.linspace(m1min, m1max, npts)]
    m1m2_points += [(m1min, m2) for m2 in numpy.linspace(m2min, m2max, npts)]
    m1m2_points += [(m1max, m2) for m2 in numpy.linspace(m2min, m2max, npts)]

    # find where total mass meets mass1, mass2
    # only need the corners here because these constraints
    # are also lines in tau0-tau3 space
    for val in constraints['mtotal']:
        m1m2_points += [(m1, val-m1) for m1 in numpy.linspace(m1min, m1max, npts)]
        m1m2_points += [(val-m2, m2) for m2 in numpy.linspace(m2min, m2max, npts)]

    # draw constant mratio lines
    for val in constraints['mratio']:
        m1m2_points += [(m1, m1/val) for m1 in numpy.linspace(m1min, m1max, npts)]
        m1m2_points += [(val*m2, m2) for m2 in numpy.linspace(m2min, m2max, npts)]

    # draw constant mchirp lines
    mcmin, mcmax = constraints.setdefault('mchirp', (None, None))
    if mcmax is not None:
        m1m2_points += [(m1, mchirpm1_to_m2(mcmax, m1)) for m1 in numpy.linspace(m1min, m1max, npts)]
    if mcmin is not None:
        m1m2_points += [(m1, mchirpm1_to_m2(mcmin, m1)) for m1 in numpy.linspace(m1min, m1max, npts)]

    # filter these down to only those that satisfy ALL constraints
    m1 = numpy.array([i[0] for i in m1m2_points])
    m2 = numpy.array([i[1] for i in m1m2_points])
    m1, m2 = allowed_m1m2(m1, m2, constraints)

    if len(m1) == 0:
        raise ValueError("The requested parameter space is empty.")

    # infer lower and upper tau0-tau3 constraints
    lims_tau0, lims_tau3 = m1m2_to_tau0tau3(m1, m2, flow)
    lims_tau0 = [min( lims_tau0 ), max( lims_tau0 )]
    lims_tau3 = [min( lims_tau3 ), max( lims_tau3 )]

    return lims_tau0, lims_tau3


def urand_tau0tau3_generator(flow, **constraints):
    """
    This is a generator for random (m1, m2) pairs that are uniformly
    distributed in (tau0, tau3) space, subject to the constraints given in
    the inputs.

    @param flow UNDOCUMENTED
    @param constraints: must specify a mass1 range; mtotal, q, and tau0
    ranges are optional. The arguments for each of these keywords should be a
    tuple of (min, max). E.g., mtotal = (50, 100) constrains
    50 <= mtotal < 100.

    Example:
    >>> urand_tau0tau3 = urand_tau0tau3_generator(40, mass1 = (1, 99), mtotal = (55, 100), mratio = (1, 10))
    >>> urand_tau0tau3.next()
    (49.836271184652254, 6.6152675269639829)
    >>> urand_tau0tau3.next()
    (46.701700599815645, 14.07462196069671)

    Equivalently, use a dictionary with **:
    >>> urand_tau0tau3 = urand_tau0tau3_generator(40, **{"mass1": (1, 99), "mtotal": (55, 100), "mratio": (1, 10)})
    """
    # check that minimal constraints are set and set defaults
    constraints = set_default_constraints(constraints)

    # draw a box around the parameter space in tau0-tau3 coords
    lims_tau0, lims_tau3 = tau0tau3_bound(flow, **constraints)
    tau0_min, tau0_max = lims_tau0
    tau3_min, tau3_max = lims_tau3

    # avoid repetitive lookups
    mass1_min, mass1_max = constraints['mass1']
    mass2_min, mass2_max = constraints['mass2']
    mtotal_min, mtotal_max = constraints['mtotal']
    qmin, qmax = constraints['mratio']

    # precompute useful coefficients
    _A0 = A0(flow)
    _A3 = A3(flow)
    A0_A3 = _A0 / _A3

    # The first part of the while loop can be the tight inner loop for
    # high-mass banks. Let's go crazy optimizing.
    # FIXME: This can be coded without discards if someone is willing to do
    # the math on the non-linear shape of the tau0, tau3 boundaries.
    from numpy.random.mtrand import uniform
    minus_five_thirds = -5. / 3.

    while 1:   # NB: "while 1" is inexplicably much faster than "while True"
        tau0 = uniform(tau0_min, tau0_max)

        mtot = A0_A3 * uniform(tau3_min, tau3_max) / tau0  # seconds
        eta = _A0 / tau0 * mtot**minus_five_thirds
        if eta > 0.25: continue

        mtot /= MTSUN_SI  # back to solar masses
        mass1 = mtot * (0.5 + sqrt(0.25 - eta)) # mass1 is the larger component
        mass2 = mtot - mass1

        if mtotal_min < mtot < mtotal_max and \
                mass1_min <= mass1 <= mass1_max and \
                mass2_min <= mass2 <= mass2_max and \
                qmin < mass1/mass2 < qmax:
            yield mass1, mass2


def IMRPhenomB_param_generator(flow, **kwargs):
    """
    Specify the min and max mass of the bigger component, then
    the min and max mass of the total mass. This function includes
    restrictions on q and chi based on IMRPhenomB's range of believability.
    Ref: http://arxiv.org/pdf/0908.2090

    @param flow UNDOCUMENTED
    @param kwargs: must specify a component_mass range; mtotal, q, chi, and tau0
    ranges are optional. If no chi is specified, the IMRPhenomB limits will be used.
    See urand_tau0tau3_generator for more usage help.
    """
    # get args

    # FIXME: PhenomB ignores spin2 and therefore requires symmetry in
    # the spins. In BBH use cases, this is OK, but for NSBH
    # applications this is undesired. The weird chi-q bounds make
    # implementing this trick
    smin, smax = kwargs.pop('spin1', (-1.0, 1.0))
    kwargs.pop('spin2')
    Warning("PhenomB: spin2 limits not implemented. Using spin1 limits for both components.")
    # the rest will be checked in the call to urand_tau0tau3_generator

    # IMRPhenomB has special bounds on chi, so we will silently truncate
    chi_low_bounds = (max(-0.85, smin), min(0.85, smax))
    chi_high_bounds = (max(-0.5, smin), min(0.75, smax))
    for mass1, mass2 in urand_tau0tau3_generator(flow, **kwargs):
        q = max(mass1/mass2, mass2/mass1)
        if 4 < q <= 10:
            spin1 = uniform(*chi_high_bounds) #guaranteed to give chi in correct range
            spin2 = uniform(*chi_high_bounds)
        elif 1 <= q <= 4:
            spin1 = uniform(*chi_low_bounds) #guaranteed to give chi in correct range
            spin2 = uniform(*chi_low_bounds)
        else:
            raise ValueError("mass ratio out of range")
        yield mass1, mass2, spin1, spin2


def IMRPhenomC_param_generator(flow, **kwargs):
    """
    Generate random parameters for the IMRPhenomC waveform model.
    Specify the min and max mass of the bigger component, then the min
    and max mass of the total mass. This function includes
    restrictions on q and chi based on IMRPhenomC's range of
    believability, namely q <=20 and |chi| <= 0.9.

    @param flow: low frequency cutoff
    @param kwargs: constraints on waveform parameters. See urand_tau0tau3_generator for more usage help. If no spin limits are specified, the IMRPhenomC limits will be used.
    """

    # get spin limits. IMRPhenomC has special bounds on chi, so we
    # will silently truncate
    s1min, s1max = kwargs.pop('spin1', (-0.9, 0.9))
    s2min, s2max = kwargs.pop('spin2', (s1min, s1max))
    s1min, s1max = (max(-0.9, s1min), min(0.9, s1max))
    s2min, s2max = (max(-0.9, s2min), min(0.9, s2max))

    for mass1, mass2 in urand_tau0tau3_generator(flow, **kwargs):

        q = max(mass1/mass2, mass2/mass1)
        if q <= 20:
            spin1 = uniform(s1min, s1max)
            spin2 = uniform(s2min, s2max)
        else:
            raise ValueError("mass ratio out of range")

        yield mass1, mass2, spin1, spin2


def aligned_spin_param_generator(flow, **kwargs):
    """
    Specify the min and max mass of the bigger component, the min and
    max mass of the total mass and the min and max values for the
    z-axis spin angular momentum.
    """
    # get args
    spin1_bounds = kwargs.pop('spin1', (-1., 1.))
    spin2_bounds = kwargs.pop('spin2', (-1., 1.))
    # the rest will be checked in the call to urand_tau0tau3_generator

    for mass1, mass2 in urand_tau0tau3_generator(flow, **kwargs):
        spin1 = uniform(*spin1_bounds)
        spin2 = uniform(*spin2_bounds)
        yield mass1, mass2, spin1, spin2


def SpinTaylorT4_param_generator(flow, **kwargs):
    # FIXME implement!
    raise NotImplementedError

proposals = {"IMRPhenomB":IMRPhenomB_param_generator,
             "IMRPhenomC":IMRPhenomC_param_generator,
             "TaylorF2RedSpin":aligned_spin_param_generator,
             "EOBNRv2":urand_tau0tau3_generator,
             "SEOBNRv1":aligned_spin_param_generator,
             "SpinTaylorT4":SpinTaylorT4_param_generator}
