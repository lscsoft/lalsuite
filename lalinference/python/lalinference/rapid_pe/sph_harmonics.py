import math
from math import pi
import numpy as np
import lal

_a2m2 = math.sqrt(5.0 / (64.0 * math.pi))
_a2m1 = math.sqrt(5.0 / (16.0 * math.pi)) 
_a20  = math.sqrt(15.0 /(32.0 * math.pi))
_a21  = math.sqrt(5.0 / (16.0 * math.pi))
_a22  = math.sqrt(5.0 / (64.0 * math.pi))
def _compute_sph_l_eq_2(theta, phi, selected_modes=None):
    """
    Compute spherical harmonics for various m with l=2.
    """
    Ylms = {}
    one_m_theta = 1.0 - np.cos(theta)
    one_p_theta = 1.0 + np.cos(theta)
    snTheta = np.sin(theta)
    for m in range(-2, 3):
        if selected_modes is not None and (2,m) not in selected_modes:
            continue
        if m == -2:
            Ylms[(2,m)] = _a2m2 * one_m_theta * one_m_theta * np.exp(1j*m*phi)
        elif m == -1:
            Ylms[(2,m)] = _a2m1 * snTheta * one_m_theta * np.exp(1j*m*phi)
        elif m == 0:
            Ylms[(2,m)] = _a20 * snTheta * snTheta #* np.exp(1j*m*phi)
        elif m == 1:
            Ylms[(2,m)] = _a21 * snTheta * one_p_theta * np.exp(1j*m*phi)
        elif m == 2:
            Ylms[(2,m)] = _a22 * one_p_theta * one_p_theta * np.exp(1j*m*phi)

    return Ylms

_a3m3 = math.sqrt(21.0 / (2.0 * math.pi))
_a3m2 = math.sqrt(7.0 / (4.0 * math.pi))
_a3m1 = math.sqrt(35.0 / (2.0 * math.pi)) 
_a30  = math.sqrt(105.0 /(2.0 * math.pi))
_a31  = -math.sqrt(35.0 / (2.0 * math.pi))
_a32  = math.sqrt(7.0 / math.pi)
_a33  = -math.sqrt(21.0 / (2.0 * math.pi))

def _compute_sph_l_eq_3(theta, phi, selected_modes=None):
    """
    Compute spherical harmonics for various m with l=3.
    """
    Ylms = {}
    # Precompute the most used (or most likely to be used) ones
    cos2th = np.cos(theta / 2.0)
    costh = np.cos(theta)
    sin2th = np.sin(theta / 2.0)
    sinth = np.sin(theta)
    for m in range(-3, 4):
        if selected_modes is not None and (3,m) not in selected_modes:
            continue
        if m == -3:
            Ylms[(3,m)] = _a3m3 * cos2th * sin2th**5.0 * np.exp(1j*m*phi)
        if m == -2:
            Ylms[(3,m)] = _a3m2 * (2.0 + 3.0 * costh) * sin2th**4.0 * np.exp(1j*m*phi)
        if m == -1:
            Ylms[(3,m)] = _a3m1 * (sinth + 4.0 * np.sin(2.0 * theta) \
                    - 3.0 * np.sin(3.0 * theta)) / 32.0 * np.exp(1j*m*phi)
        if m == 0:
            Ylms[(3,m)] = _a30 * costh * sinth**2.0 / 4.0 #deleted an ) at the end of the 2.0
        if m == 1:
            Ylms[(3,m)] = _a31 * (sinth - 4.0 * np.sin(2.0 * theta) - 3.0 * np.sin(3.0 * theta)) / 32.0 * np.exp(1j*m*phi)
        if m == 2:
            Ylms[(3,m)] = _a32 * cos2th**4.0 * \
                  (-2.0 + 3.0 * costh) / 2.0 * np.exp(1j*m*phi)
        if m == 3:
            Ylms[(3,m)] = _a33 * cos2th**5.0 * sin2th * np.exp(1j*m*phi)

    return Ylms


_a4m4 = 3.0 * math.sqrt(7.0 / math.pi)
_a4m3 = 3.0 * math.sqrt(7.0 / 2.0 / math.pi)
_a4m2 =	3.0 / (4.0 * math.sqrt(math.pi)) 
_a4m1 =	3.0 / (32.0 * math.sqrt(2.0 * math.pi)) 
_a40  =	3.0 * math.sqrt(5.0 / (2.0 * math.pi)) / 16.0 
_a41  =	3.0 / (32.0 * math.sqrt(2.0 * math.pi))			 
_a42  = 3.0 / (4.0 * math.sqrt(math.pi)) 
_a43  =	-3.0 * math.sqrt(7.0 / (2.0 * math.pi)) 
_a44  = 3.0 * math.sqrt(7.0 / math.pi)

def _compute_sph_l_eq_4(theta, phi, selected_modes=None):
    """
    Compute spherical harmonics for various m with l=4.
    """
    Ylms = {}
    # precompute some terms
    cos2th = np.cos(theta / 2.0)
    costh = np.cos(theta)
    sin2th = np.sin(theta / 2.0)
    sinth = np.sin(theta)
    for m in range(-4, 5):
        if selected_modes is not None and (4,m) not in selected_modes:
            continue
        if m == -4:			
            Ylms [(4,m)] = _a4m4 * (cos2th**2.0 * sin2th**6.0) * np.exp(1j*m*phi)
        if m == -3:			
            Ylms [(4,m)] = _a4m3 * (cos2th * (1.0 + 2.0 * costh) * sin2th**5.0) * np.exp(1j*m*phi)
        if m == -2:			
            Ylms [(4,m)] = _a4m2 * ((9.0 + 14.0 * costh + 7.0*np.cos(2.0*theta)) * sin2th**4.0) * np.exp(1j*m*phi) 
        if m == -1:			
            Ylms [(4,m)] = _a4m1 * (3.0 * sinth + 2.0 * np.sin(2.0 * theta) + 7.0 * np.sin(3.0 * theta) - 7.0 * np.sin(4.0 * theta)) * np.exp(1j*m*phi) 
        if m ==  0:			
            Ylms [(4,m)] = _a40  * ((5.0 + 7.0 * np.cos(2.0 * theta)) * sinth**2.0) * np.exp(1j*m*phi)
        if m ==  1:			
            Ylms [(4,m)] = _a41  * (3.0 * sinth - 2.0 * np.sin(2.0 * theta) + 7.0 * np.sin(3.0 * theta) + 7.0 * np.sin(4.0 * theta)) * np.exp(1j*m*phi) 
        if m ==  2:			
            Ylms [(4,m)] = _a42  * (cos2th**4.0 * (9.0 - 14.0*costh + 7.0 * np.cos(2.0 * theta))) * np.exp(1j*m*phi)
        if m ==  3:			
            Ylms [(4,m)] = _a43  * (cos2th**5.0 * (-1.0 + 2.0 * costh) * sin2th) * np.exp(1j*m*phi)
        if m ==  4:			
            Ylms [(4,m)] = _a44  * (cos2th**6.0 * sin2th**2.0) * np.exp(1j*m*phi)

    return Ylms


_a5m5 = math.sqrt(330.0 / math.pi)
_a5m4 = math.sqrt(33.0 / math.pi) 
_a5m3 = math.sqrt(33.0 / (2.0 * math.pi)) / 4.0
_a5m2 = math.sqrt(11.0 / math.pi) / 8.0 
_a5m1 = math.sqrt(77.0 / math.pi) / 256.0
_a50  = math.sqrt(1155.0 / (2.0 * math.pi)) / 32.0 
_a51  = math.sqrt(77.0 / math.pi) / 256.0 
_a52  = math.sqrt(11.0 / math.pi) / 8.0 
_a53  = -math.sqrt(33.0 / (2.0 * math.pi)) / 4.0 
_a54  = math.sqrt(33.0 / math.pi) 
_a55  = -math.sqrt(330.0 / math.pi) 


def _compute_sph_l_eq_5(theta, phi, selected_modes=None):
    """
    Compute spherical harmonics for various m with l=5.
    """
    Ylms = {}
    # precompute some terms
    cos2th = np.cos(theta / 2.0)
    costh = np.cos(theta)
    sin2th = np.sin(theta / 2.0)
    sinth = np.sin(theta)

    for m in range(-5, 6):
        if selected_modes is not None and (5,m) not in selected_modes:
            continue
        if m == -5: 
            Ylms[(5,m)] = _a5m5 * (cos2th**3.0 * sin2th**7.0) * np.exp(1j*m*phi)
        if m == -4: 
            Ylms[(5,m)] = _a5m4 * (cos2th**2.0 * (2.0 + 5.0 * costh) * sin2th**6.0) * np.exp(1j*m*phi)
        if m == -3: 
            Ylms[(5,m)] = _a5m3 * (17.0 + 24.0 * np.cos(theta) + 15.0 * np.cos(2.0 * theta)) * cos2th * sin2th**5.0 * np.exp(1j*m*phi)
        if m == -2: 
            Ylms[(5,m)] = _a5m2 * (32.0 + 57.0 * costh + 36.0 * np.cos(2.0 * theta) + 15.0 * np.cos(3.0 * theta)) * sin2th**4.0 * np.exp(1j*m*phi)
        if m == -1: 
            Ylms[(5,m)] = _a5m1 * (2.0 * sinth + 8.0 * np.sin(2.0 * theta) + 3.0 * np.sin(3.0*theta) + 12.0 * np.sin(4.0 * theta) - 15.0 * np.sin(5.0 * theta)) * np.exp(1j*m*phi)
        if m == 0: 
            Ylms[(5,m)] = _a50 * (5.0 * costh + 3.0 * np.cos(3.0 * theta)) * sinth**2.0 * np.exp(1j*m*phi)
        if m ==  1: 
            Ylms[(5,m)] = _a51 * (-2.0 * sinth + 8.0 * np.sin(2.0 * theta) - 3.0 * np.sin(3.0 * theta) + 12.0 * np.sin(4.0*theta) + 15.0 * np.sin(5.0*theta)) * np.exp(1j*m*phi)
        if m ==  2: 
            Ylms[(5,m)] = _a52 * (-32.0 + 57.0 * costh - 36.0 * np.cos(2.0 * theta) + 15.0 * np.cos(3.0 * theta))* cos2th**4.0 * np.exp(1j*m*phi)
        if m ==  3: 
            Ylms[(5,m)] = _a53 * (17.0 - 24.0 * costh + 15.0 * np.cos(2.0 * theta)) * cos2th**5.0 * sin2th * np.exp(1j*m*phi)
        if m ==  4: 
            Ylms[(5,m)] = _a54 * (-2.0 + 5.0 * costh) * cos2th**6.0 * sin2th**2  * np.exp(1j*m*phi)
        if m ==  5: 
            Ylms[(5,m)] = _a55 * cos2th**7.0 * sin2th **3.0 * np.exp(1j*m*phi)


    return Ylms

_a6m6 = 3.0 * math.sqrt(715.0 / math.pi) / 2.0 
_a6m5 = math.sqrt(2145.0 / math.pi) / 2.0 
_a6m4 = math.sqrt(195.0 / (2.0 * math.pi)) / 8.0
_a6m3 = 3.0 * math.sqrt(13.0 / math.pi) / 32.0 
_a6m2 = math.sqrt(13.0 / math.pi) / 256.0 
_a6m1 = math.sqrt(65.0 / (2.0 * math.pi)) / 64.0 
_a60  = math.sqrt(1365.0 / math.pi) / 512.0 
_a61  = math.sqrt(65.0 / (2.0 * math.pi)) / 64.0
_a62  = math.sqrt(13.0 / math.pi) / 256.0 
_a63  = -3.0 * math.sqrt(13.0 / math.pi) / 32.0
_a64  = math.sqrt(195.0 / (2.0 * math.pi)) / 8.0
_a65  = -math.sqrt(2145.0 / math.pi) / 2.0
_a66  = 3.0 * math.sqrt(715.0 / math.pi) / 2.0

def _compute_sph_l_eq_6(theta, phi, selected_modes=None):
    """
    Compute spherical harmonics for various m with l=6.
    """
    Ylms = {}
    # precompute some terms
    cos2th = np.cos(theta / 2.0)
    costh = np.cos(theta)
    sin2th = np.sin(theta / 2.0)
    sinth = np.sin(theta)
    for m in range(-6, 7):
        if selected_modes is not None and (6,m) not in selected_modes:
            continue
        if m == -6:
            Ylms[(6,m)] = _a6m6 * cos2th**4.0 * sin2th**8.0 * np.exp(1j*m*pi)
        if m == -5:
            Ylms[(6,m)] = _a6m5 * (1.0 + 3.0 * costh) * cos2th**3 * sin2th**7.0 * np.exp(1j*m*pi)
        if m == -4:
            Ylms[(6,m)] = _a6m4 * (35.0 + 44.0 * costh + 33.0 * np.cos(2.0 * theta)) * cos2th**2 * sin2th**6 * np.exp(1j*m*pi)
        if m == -3:
            Ylms[(6,m)] = _a6m3 * (98.0 + 185.0 * costh + 110.0 * np.cos(2.0 * theta) + 55.0 * np.cos(3.0 * theta)) * cos2th * cos2th**5.0 * np.exp(1j*m*pi)
        if m == -2:
            Ylms[(6,m)] = _a6m2 * (1709.0 + 3096.0 * costh + 2340 * np.cos(2.0 * theta) + 1320.0 * np.cos(3.0 *theta) + 495.0 * np.cos(4.0 * theta)) * sin2th**4.0 * np.exp(1j*m*pi)
        if m == -1:
            Ylms[(6,m)] = _a6m1 * (161.0 + 252.0 * costh + 252.0 * np.cos(2.0 * theta) + 132.0 * np.cos(3.0 * theta) + 99.0 * np.cos(4.0 * theta)) * cos2th * sin2th**3.0 * np.exp(1j*m*pi)
        if m ==  0:
            Ylms[(6,m)] = _a60  * (35.0 + 60.0 * np.cos(2.0 * theta) + 33.0 * np.cos(4.0 * theta)) * sin2th**2.0 * np.exp(1j*m*pi)
        if m ==  1:
            Ylms[(6,m)] = _a61  * (161.0 - 252.0 * costh + 252.0 * np.cos(2.0 * theta) - 132.0 * np.cos(3.0 * theta) + 99.0 * np.cos(4.0 * theta)) * cos2th**3.0 * sin2th * np.exp(1j*m*pi)
        if m ==  2:
            Ylms[(6,m)] = _a62  * (1709.0 - 3096.0 * costh + 2340.0 * np.cos(2.0*theta) - 1320.0 * np.cos(3.0 * theta) + 495.0 * np.cos(4.0 * theta)) * cos2th**4.0 * np.exp(1j*m*pi)
        if m ==  3:
            Ylms[(6,m)] = _a63  * (-98.0 + 185.0 * costh - 110.0 * np.cos(2.0 * theta) + 55.0 * np.cos(3.0*theta)) * cos2th**5.0 * sin2th * np.exp(1j*m*pi)
        if m ==  4:
            Ylms[(6,m)] = _a64  * (35.0 - 44.0 * costh + 33.0 * np.cos(2.0 * theta)) * cos2th**6.0 * sin2th**2.0 * np.exp(1j*m*pi)
        if m ==  5:
            Ylms[(6,m)] = _a65  * (-1.0 + 3.0 * costh) * cos2th**7.0 * sin2th**3.0 * np.exp(1j*m*pi)
        if m ==  6:
            Ylms[(6,m)] = _a66  * (cos2th**8.0 * sin2th**4.0) * np.exp(1j*m*pi)

    return Ylms

def compute_spherical_harmonics(Lmax, theta, phi, selected_modes=None):
    """
    Return a dictionary keyed by tuples
    (l,m)
    that contains the values of all
    -2Y_lm(theta,phi)
    with
    l <= Lmax
    -l <= m <= l
    """
    Ylms = _compute_sph_l_eq_2(theta, phi, selected_modes)
    Ylms.update(_compute_sph_l_eq_3(theta, phi, selected_modes))
    Ylms.update(_compute_sph_l_eq_4(theta, phi, selected_modes))
    Ylms.update(_compute_sph_l_eq_5(theta, phi, selected_modes))
    Ylms.update(_compute_sph_l_eq_6(theta, phi, selected_modes))

    return Ylms

if __name__ == "__main__":
    theta, phi = np.random.uniform(0, pi, 2)

    # Do unit tests
    Ylm = compute_spherical_harmonics(2, theta, phi, None)
    for (l, m), val in Ylm.iteritems():
        if np.isclose(lal.SpinWeightedSphericalHarmonic(theta, phi, -2, l, m), val):	
            print("Test successful for (l,m) = (%d, %d)" % (l, m))
        else:
            print("Test unsucessful for (l,m) = (%d, %d)" % (l, m))
