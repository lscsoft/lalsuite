import numpy
import lal
from lalinference.rapid_pe import lalsimutils
m1m2 = numpy.vectorize(lalsimutils.m1m2)

def midpoint(pt1, pt2):
    diff = pt2 - pt1
    mid = diff / 2
    return pt1 + mid

def prune_duplicate_pts(pts):
    return numpy.array(list(set([tuple(pt) for pt in pts])))

def transform_m1m2_mceta(m1, m2):
    return lalsimutils.Mceta(m1, m2)

def transform_mceta_m1m2(mc, eta):
    return m1m2(mc, eta)

__prefac_0 = 5. / 256 / numpy.pi
__prefac_3 = 1. / 8
__dim_mass = lal.G_SI / lal.C_SI**3 * lal.MSUN_SI
def transform_m1m2_tau0tau3(m1, m2, flow=40.):
    mt = m1 + m2
    eta = m1 * m2 / mt**2
    mt *= numpy.pi * flow * __dim_mass
    return (__prefac_0 / flow / eta * mt**(-5./3), __prefac_3 / flow / eta * mt**(-2./3))

__prefac_tau = 5. / 32 / numpy.pi
def transform_tau0tau3_m1m2(tau0, tau3, flow=40.):
    mt = __prefac_tau / numpy.pi / flow * tau3 / tau0
    eta = 1.0 / 8 / flow / tau3 * (tau0 / __prefac_tau / tau3)**(2./3)
    m1, m2 = m1m2(mt*eta**(3./5), eta)
    return m1 / __dim_mass, m2 / __dim_mass

VALID_TRANSFORMS_MASS = { \
    "mchirp_eta": transform_m1m2_mceta,
    "tau0_tau3": transform_m1m2_tau0tau3,
    None: None
}

INVERSE_TRANSFORMS_MASS = { \
    transform_m1m2_mceta: transform_mceta_m1m2,
    transform_m1m2_tau0tau3: transform_tau0tau3_m1m2,
    None: None
}

def apply_transform(pts, intr_prms, mass_transform=None):
    # You know what... recarrays are dumb, and so's your face.
    # FIXME: Why does numpy want me to repack this into tuples!?
    #tpts = numpy.array([tuple(pt) for pt in pts.T], dtype = numpy.dtype([(a ,"float64") for a in intr_prms]))
    m1_idx, m2_idx = intr_prms.index("mass1"), intr_prms.index("mass2")
    if mass_transform:
       pts[:,m1_idx], pts[:,m2_idx] = VALID_TRANSFORMS_MASS[mass_transform](pts[:,m1_idx], pts[:,m2_idx])

    # Independent transforms go here

    return pts

def apply_inv_transform(pts, intr_prms, mass_transform=None):
    m1_idx, m2_idx = intr_prms.index("mass1"), intr_prms.index("mass2")
    if mass_transform:
       pts[:,m1_idx], pts[:,m2_idx] = INVERSE_TRANSFORMS_MASS[VALID_TRANSFORMS_MASS[mass_transform]](pts[:,m1_idx], pts[:,m2_idx])

    # Independent transforms go here

    return pts
