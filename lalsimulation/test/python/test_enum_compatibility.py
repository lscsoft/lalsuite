#!/usr/bin/env python3

"""
This script checks that the list of approximants
has not been changed in a backward-incompatible way
by reordering the Approximant enum.
"""
import sys
import pytest
import lalsimulation as lalsim

# This list was last updated 2020-03-30 by john.veitch
known_approximants = {
 'TaylorT1': 0,
 'TaylorT2': 1,
 'TaylorT3': 2,
 'TaylorF1': 3,
 'EccentricFD': 4,
 'TaylorF2': 5,
 'TaylorF2Ecc': 6,
 'TaylorF2NLTides': 7,
 'TaylorR2F4': 8,
 'TaylorF2RedSpin': 9,
 'TaylorF2RedSpinTidal': 10,
 'PadeT1': 11,
 'PadeF1': 12,
 'EOB': 13,
 'BCV': 14,
 'BCVSpin': 15,
 'SpinTaylorT1': 16,
 'SpinTaylorT3': 18,
 'SpinTaylorT4': 19,
 'SpinTaylorT5': 20,
 'SpinTaylorF2': 21,
 'SpinTaylorFrameless': 22,
 'SpinTaylor': 23,
 'PhenSpinTaylor': 24,
 'PhenSpinTaylorRD': 25,
 'SpinQuadTaylor': 26,
 'FindChirpSP': 27,
 'FindChirpPTF': 28,
 'GeneratePPN': 29,
 'BCVC': 30,
 'FrameFile': 31,
 'AmpCorPPN': 32,
 'NumRel': 33,
 'NumRelNinja2': 34,
 'Eccentricity': 35,
 'EOBNR': 36,
 'EOBNRv2': 37,
 'EOBNRv2HM': 38,
 'EOBNRv2_ROM': 39,
 'EOBNRv2HM_ROM': 40,
 'TEOBResum_ROM': 41,
 'SEOBNRv1': 42,
 'SEOBNRv2': 43,
 'SEOBNRv2_opt': 44,
 'SEOBNRv3': 45,
 'SEOBNRv3_pert': 46,
 'SEOBNRv3_opt': 47,
 'SEOBNRv3_opt_rk4': 48,
 'SEOBNRv4': 49,
 'SEOBNRv4_opt': 50,
 'SEOBNRv4P': 51,
 'SEOBNRv4PHM': 52,
 'SEOBNRv2T': 53,
 'SEOBNRv4T': 54,
 'SEOBNRv1_ROM_EffectiveSpin': 55,
 'SEOBNRv1_ROM_DoubleSpin': 56,
 'SEOBNRv2_ROM_EffectiveSpin': 57,
 'SEOBNRv2_ROM_DoubleSpin': 58,
 'SEOBNRv2_ROM_DoubleSpin_HI': 59,
 'Lackey_Tidal_2013_SEOBNRv2_ROM': 60,
 'SEOBNRv4_ROM': 61,
 'SEOBNRv4_ROM_NRTidal': 62,
 'SEOBNRv4T_surrogate': 63,
 'HGimri': 64,
 'IMRPhenomA': 65,
 'IMRPhenomB': 66,
 'IMRPhenomFA': 67,
 'IMRPhenomFB': 68,
 'IMRPhenomC': 69,
 'IMRPhenomD': 70,
 'IMRPhenomD_NRTidal': 71,
 'IMRPhenomHM': 72,
 'IMRPhenomP': 73,
 'IMRPhenomPv2': 74,
 'IMRPhenomPv2_NRTidal': 75,
 'IMRPhenomFC': 76,
 'TaylorEt': 77,
 'TaylorT4': 78,
 'EccentricTD': 79,
 'TaylorN': 80,
 'SpinTaylorT4Fourier': 81,
 'SpinTaylorT5Fourier': 82,
 'SpinDominatedWf': 83,
 'NR_hdf5': 84,
 'NRSur4d2s': 85,
 'NRSur7dq2': 86,
 'NRSur7dq4': 87,
 'SEOBNRv4HM': 88,
 'NRHybSur3dq8': 89,
 'SEOBNRv4HM_ROM': 90,
 'SEOBNRv4_ROM_NRTidalv2': 91,
 'SEOBNRv4_ROM_NRTidalv2_NSBH': 92,
 'IMRPhenomD_NRTidalv2': 93,
 'IMRPhenomNSBH': 94,
 'IMRPhenomPv2_NRTidalv2': 95,
 'IMRPhenomXAS': 96,
 'IMRPhenomXHM': 97,
 'IMRPhenomPv3': 98,
 'IMRPhenomPv3HM': 99
}

@pytest.mark.parametrize("name, i", known_approximants.items())
def test_approximant(name, i):
    a = lalsim.GetApproximantFromString(name)
    assert a == i, (
        "The Approximant enum is modified in an non-backward-compatible way, "
        "Approximant {name} is {a} but it should be {i}; "
        "please append new approximants to the Approximants enum in LALSimInspiral.h, "
        "rather than reordering the list".format(
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
            appx[lalsim.GetStringFromApproximant(i)]=i
        except RuntimeError:
            # Integer which is not used in enum?
            print('No approximant numbered {}'.format(i))
    print(appx)
    return appx

if __name__ == '__main__':
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-approx-enum.xml"]
    sys.exit(pytest.main(args=[__file__] + args))

