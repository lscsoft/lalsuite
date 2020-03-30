#!/usr/bin/env python3

"""
This script checks that the list of approximants
has not been changed in a backward-incompatible way
by reordering the Approximant enum.
"""

import pytest
import lalsimulation as lalsim

# This list was last updated 2020-03-30 by john.veitch
known_approximants={
    0: 'TaylorT1',
    1: 'TaylorT2',
    2: 'TaylorT3',
    3: 'TaylorF1',
    4: 'EccentricFD',
    5: 'TaylorF2',
    6: 'TaylorF2Ecc',
    7: 'TaylorF2NLTides',
    8: 'TaylorR2F4',
    9: 'TaylorF2RedSpin',
    10: 'TaylorF2RedSpinTidal',
    11: 'PadeT1',
    12: 'PadeF1',
    13: 'EOB',
    14: 'BCV',
    15: 'BCVSpin',
    16: 'SpinTaylorT1',
    17: 'SpinTaylorT2',
    18: 'SpinTaylorT3',
    19: 'SpinTaylorT4',
    20: 'SpinTaylorT5',
    21: 'SpinTaylorF2',
    22: 'SpinTaylorFrameless',
    23: 'SpinTaylor',
    24: 'PhenSpinTaylor',
    25: 'PhenSpinTaylorRD',
    26: 'SpinQuadTaylor',
    27: 'FindChirpSP',
    28: 'FindChirpPTF',
    29: 'GeneratePPN',
    30: 'BCVC',
    31: 'FrameFile',
    32: 'AmpCorPPN',
    33: 'NumRel',
    34: 'NumRelNinja2',
    35: 'Eccentricity',
    36: 'EOBNR',
    37: 'EOBNRv2',
    38: 'EOBNRv2HM',
    39: 'EOBNRv2_ROM',
    40: 'EOBNRv2HM_ROM',
    41: 'TEOBResum_ROM',
    42: 'SEOBNRv1',
    43: 'SEOBNRv2',
    44: 'SEOBNRv2_opt',
    45: 'SEOBNRv3',
    46: 'SEOBNRv3_pert',
    47: 'SEOBNRv3_opt',
    48: 'SEOBNRv3_opt_rk4',
    49: 'SEOBNRv4',
    50: 'SEOBNRv4_opt',
    51: 'SEOBNRv2T',
    52: 'SEOBNRv4T',
    53: 'SEOBNRv1_ROM_EffectiveSpin',
    54: 'SEOBNRv1_ROM_DoubleSpin',
    55: 'SEOBNRv2_ROM_EffectiveSpin',
    56: 'SEOBNRv2_ROM_DoubleSpin',
    57: 'SEOBNRv2_ROM_DoubleSpin_HI',
    58: 'Lackey_Tidal_2013_SEOBNRv2_ROM',
    59: 'SEOBNRv4_ROM',
    60: 'SEOBNRv4_ROM_NRTidal',
    61: 'SEOBNRv4T_surrogate',
    62: 'HGimri',
    63: 'IMRPhenomA',
    64: 'IMRPhenomB',
    65: 'IMRPhenomFA',
    66: 'IMRPhenomFB',
    67: 'IMRPhenomC',
    68: 'IMRPhenomD',
    69: 'IMRPhenomD_NRTidal',
    70: 'IMRPhenomHM',
    71: 'IMRPhenomP',
    72: 'IMRPhenomPv2',
    73: 'IMRPhenomPv2_NRTidal',
    74: 'IMRPhenomFC',
    75: 'TaylorEt',
    76: 'TaylorT4',
    77: 'EccentricTD',
    78: 'TaylorN',
    79: 'SpinTaylorT4Fourier',
    80: 'SpinTaylorT2Fourier',
    81: 'SpinDominatedWf',
    82: 'NR_hdf5',
    83: 'NRSur4d2s',
    84: 'NRSur7dq2',
    85: 'SEOBNRv4HM',
    86: 'NRHybSur3dq8'
}

@pytest.mark.parametrize("i, name", known_approximants.items())
def test_approximant(i, name):
    a = lalsim.GetStringFromApproximant(i)
    assert a == name, (
        "The Approximant enum is modified in an non-backward-compatible way, "
        "Approximant {i} is {a} but it should be {name}; "
        "please append {a} to the Approximants enum, rather than reordering the list".format(
            i=i,
            a=a,
            name=name,
        )
    )


def regenerate_known_approximants():
    """
    Function here for convenience if you need to regenerate the list
    Not used in test
    """
    appx = {}
    for i in range(lalsim.NumApproximants):
        try:
            appx[i]=lalsim.GetStringFromApproximant(i)
        except RuntimeError:
            # Integer which is not used in enum?
            print('No approximant numbered {}'.format(i))
    print(appx)
    return appx

if __name__ == '__main__':
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-approx-enum.xml"]
    sys.exit(pytest.main(args=[__file__] + args))

