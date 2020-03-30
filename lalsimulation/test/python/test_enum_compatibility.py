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
known_approximants={
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
 'SpinTaylorT2': 17,
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
 'SEOBNRv2T': 51,
 'SEOBNRv4T': 52,
 'SEOBNRv1_ROM_EffectiveSpin': 53,
 'SEOBNRv1_ROM_DoubleSpin': 54,
 'SEOBNRv2_ROM_EffectiveSpin': 55,
 'SEOBNRv2_ROM_DoubleSpin': 56,
 'SEOBNRv2_ROM_DoubleSpin_HI': 57,
 'Lackey_Tidal_2013_SEOBNRv2_ROM': 58,
 'SEOBNRv4_ROM': 59,
 'SEOBNRv4_ROM_NRTidal': 60,
 'SEOBNRv4T_surrogate': 61,
 'HGimri': 62,
 'IMRPhenomA': 63,
 'IMRPhenomB': 64,
 'IMRPhenomFA': 65,
 'IMRPhenomFB': 66,
 'IMRPhenomC': 67,
 'IMRPhenomD': 68,
 'IMRPhenomD_NRTidal': 69,
 'IMRPhenomHM': 70,
 'IMRPhenomP': 71,
 'IMRPhenomPv2': 72,
 'IMRPhenomPv2_NRTidal': 73,
 'IMRPhenomFC': 74,
 'TaylorEt': 75,
 'TaylorT4': 76,
 'EccentricTD': 77,
 'TaylorN': 78,
 'SpinTaylorT4Fourier': 79,
 'SpinTaylorT2Fourier': 80,
 'SpinDominatedWf': 81,
 'NR_hdf5': 82,
 'NRSur4d2s': 83,
 'NRSur7dq2': 84,
 'SEOBNRv4HM': 85,
 'NRHybSur3dq8': 86
}

@pytest.mark.parametrize("name, i", known_approximants.items())
def test_approximant(name, i):
    a = lalsim.GetApproximantFromString(name)
    assert a == i, (
        "The Approximant enum is modified in an non-backward-compatible way, "
        "Approximant {name} is {a} but it should be {i}; "
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
            appx[lalsim.GetStringFromApproximant(i)]=i
        except RuntimeError:
            # Integer which is not used in enum?
            print('No approximant numbered {}'.format(i))
    print(appx)
    return appx

if __name__ == '__main__':
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-approx-enum.xml"]
    sys.exit(pytest.main(args=[__file__] + args))

