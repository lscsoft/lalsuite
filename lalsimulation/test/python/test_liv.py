# -*- coding: utf-8 -*-
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http: //www.gnu.org/licenses/>.

"""Simple test to see if the LIV parameters can be read and inserted correctly
"""

import sys
import pytest
import lal
import lalsimulation
import numpy as np

# -- utility functions ---------------------

def read_liv_params(LALparams):
    """
    Reads LIV parameters
    """

    logLeff = lalsimulation.SimInspiralWaveformParamsLookupNonGRLIVLogLambdaEff(LALparams)
    signOfA = lalsimulation.SimInspiralWaveformParamsLookupNonGRLIVASign(LALparams)
    alpha = lalsimulation.SimInspiralWaveformParamsLookupNonGRLIVAlpha(LALparams)
    return logLeff, signOfA, alpha

def set_liv_pars(LALparams, logLeff, Asign, alpha):
    """
    Sets LIV parameters according to input values.

    - logLeff: LIV wavelength-like parameter log10lambda_eff to be set.
    - Asign: LIV sign in dispersion relation, +1.0 or -1.0.
    - alpha: Exponent of momentum in dispersion relation, between 0 and 4.
    """

    lalsimulation.SimInspiralWaveformParamsInsertNonGRLIVLogLambdaEff(LALparams, logLeff)
    lalsimulation.SimInspiralWaveformParamsInsertNonGRLIVASign(LALparams, Asign)
    lalsimulation.SimInspiralWaveformParamsInsertNonGRLIVAlpha(LALparams, alpha)
    return None

def is_liv_enabled_by_default(LALparams):
    """
    This checks if LIV flag is disabled by default.
    """

    return lalsimulation.SimInspiralWaveformParamsLookupEnableLIV(LALparams)

def enable_liv(LALparams):
    """
    This enables the LIV flag, by enabling it,
    the LIV parameters will be sampled upon.
    """

    lalsimulation.SimInspiralWaveformParamsInsertEnableLIV(LALparams, 1)
    return lalsimulation.SimInspiralWaveformParamsLookupEnableLIV(LALparams)

# -- test functions ---------------------

def test_correct_liv_pars():
    """
    This tests if the default LIV parameters are correct.
    Additionally it tests LIV parameters are being inserted correctly.

    `expected_result = np.array([100.0,1.0,0.0])` are the default values of log10lambda_eff,
     A_sign and alpha respectively
    """
    LALpars = lal.CreateDict()
    expected_result = np.array([100.0, 1.0, 0.0])
    actual_result = np.array(read_liv_params(LALpars))
    np.testing.assert_almost_equal(actual_result, expected_result, 7, "Default LIV values are not set correctly")
    ## Checking if parameters can be inserted properly
    set_liv_pars(LALpars, 44.,-1.,1.5)
    expected_result = np.array([44., -1., 1.5])
    actual_result = np.array(read_liv_params(LALpars))
    np.testing.assert_almost_equal(actual_result, expected_result, 7, "LIV values are not inserted correctly")

def test_liv_flag_disabled_by_default():
    """
    This tests that the liv flag is disabled by default.
    Additionally it checks the flag may be enabled properly.

    `expected_result = 0` is the default value of the LIV flag
    """
    LALpars = lal.CreateDict()
    expected_result = 0
    actual_result = is_liv_enabled_by_default(LALpars)
    np.testing.assert_approx_equal(actual_result, expected_result, 7, "Incorrect setting of LIV flag by default")
    ## Now check if it can be inserted correctly
    expected_result = 1
    actual_result = enable_liv(LALpars)
    np.testing.assert_approx_equal(actual_result, expected_result, 7, "LIV flag not inserted correctly")

# -- run the tests ------------------------------

if __name__ == '__main__':
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-liv.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
