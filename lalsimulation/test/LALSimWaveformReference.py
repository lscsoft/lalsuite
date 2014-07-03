#!/usr/bin/env python

# Copyright (C) 2014 Frank Ohme, Evan Ochsner
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
# A copy of the GNU General Public License may be found at
# http://www.gnu.org/copyleft/gpl.html
# or write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

import lal
import lalsimulation as lalsim
# import math
import numpy as np
import unittest
from optparse import OptionParser
import ConfigParser
import io
import sys, os

NEW_DATA_STR = '######### NEW DATASET #############\n'
DEFAULT_FILE = 'reviewed_waveforms.asc'

usage = ('usage: %prog [options]\nChecks the waveform generation'
        + ' in LALSimulation against reference data.')

parser = OptionParser(usage = usage)
parser.add_option('-r', '--reference-file', action = 'store', type = 'string',
        dest = 'reffilename', default = DEFAULT_FILE,
        metavar = "FILE", help = 'location of the file containing '
        + 'reference waveform data [default: %default]')
parser.add_option('-t', '--tolerance', action = 'store', type = float,
        metavar = "TOL", dest = 'tolerance', default = 1.e-10,
        help = 'test whether data are consistent to within this ' +
        'fractional difference [default: %default]')
parser.add_option('-a', '--approximant', action = 'store', type = 'string',
        dest = 'approx', default = 'all',
        help = 'waveform approximant [default: %default]')

(options, args) = parser.parse_args()

def diff(v1, v2, meanval):
    """
    Return
        abs( (v1[i]-v2[i])/v1[i] )
    Unless v1[i]==0, then return
        abs( v1[i] - v2[i]/ meanval )
    where meanval is the expected average absolute value of v1, v2
    """
    out = np.array( [ np.abs((v1[i]-v2[i])/v1[i]) if v1[i]!=0\
            else np.abs((v1[i]-v2[i])/meanval) for i in range(len(v1)) ] )
    return out


def generatePolarizationTest(datasets):
    '''Generates test method to be added to the CheckReferenceWaveforms class.
    The input is a list of datasets extracted from the reference file.'''
    def test_approx(self):
        '''check for consistent waveform polarisations'''
        for conf in datasets:
            approx = lalsim.GetApproximantFromString(conf.get('approximant',
                    'approximant'))
            domain = conf.get('approximant', 'domain')
            approxstr = conf.get('approximant', 'approximant') + ' / ' + domain
            epochref = conf.getfloat('waveform-data', 'epoch')
            if domain=='TD':
                hpref, hcref = [np.array(map(float, (conf.get('waveform-data',
                        l)).split())) for l in ['hp', 'hc']]
            if domain=='FD':
                hpRref, hpIref, hcRref, hcIref = [np.array(map(float,
                        (conf.get('waveform-data', l)).split()))
                        for l in ['hp_real', 'hp_imag', 'hc_real', 'hc_imag']]
                hpref = hpRref + 1j * hpIref
                hcref = hcRref + 1j * hcIref

            names = self.paramnames[domain]
            parstring =' / '.join([name + ': '
                    + str(conf.get('parameters', name)) for name in names])
            parDict = {p: self.paramtype[p](conf.get('parameters', p)) for p in conf.options('parameters')}
            parDict['m1'] *= lal.LAL_MSUN_SI
            parDict['m2'] *= lal.LAL_MSUN_SI
            parDict['distance'] *= (1.e6 * lal.LAL_PC_SI)

            params = [parDict[name] for name in names]
            params.append(approx)

            hp, hc = self.waveformgenerator(domain, params)
            epoch = float(hp.epoch)

            # Actual test starts here
            self.assertTrue( np.abs(epochref-epoch) < options.tolerance,
                    self.errmsg('epoch', approxstr, parstring))
            self.assertEqual(hp.data.data.size, hpref.size,
                             self.errmsg('length of generated hplus array',
                             approxstr, parstring))
            self.assertEqual(hc.data.data.size, hcref.size,
                             self.errmsg('length of generated hcross array',
                             approxstr, parstring))
            ampmean = np.sqrt(np.abs(hpref * hpref + hcref * hcref)).mean()
            hpdiff = diff(hpref, hp.data.data, ampmean)
            hcdiff = diff(hcref, hc.data.data, ampmean)
            self.assertTrue( (hpdiff < options.tolerance).all(),
                    self.errmsg('hplus', approxstr, parstring))
            self.assertTrue( (hcdiff < options.tolerance).all(),
                    self.errmsg('hcross', approxstr, parstring))

    return test_approx

def addApproxTestToClass(approx, dataset):
    '''adds test method as "test_approx" (where approx is the approximant name)
    to the class CheckReferenceWaveforms'''
    test_method = generatePolarizationTest(dataset)
    test_method.__name__ = 'test_' + approx
    test_method.__doc__ = test_method.__doc__ + ' for ' + approx
    setattr(CheckReferenceWaveforms, test_method.__name__, test_method)



class ReferenceFile:
    '''When initialized with the filename for the reference file,
    this class contains in the information given in the file in a
    ConfigParser compatible way. Various datasets are constructed
    that represent blocks of data between the string
    ######### NEW DATASET #############'''
    def __init__(self, filename):
        infile = open(filename, 'r')
        self.name = filename
        self.content = infile.readlines()
        infile.close()
        self.size = len(self.content)
        self.newapproxindex = [i for i in range(self.size)
                if self.content[i] == NEW_DATA_STR]
        defaultparams = {'waveformFlags': None, 'nonGRparams': None}
        #possibly add more above
        ConfigParser.RawConfigParser.optionxform = str
        # prevent ConfigParser to use lower case version of option
        self.dataset = [ConfigParser.RawConfigParser(defaults = defaultparams)
                for i in range(len(self.newapproxindex))]
        self.newapproxindex.append(self.size)
        for i in range(len(self.newapproxindex) - 1):
            begin, end = self.newapproxindex[i:(i+2)]
            filepart = io.BytesIO('\n'.join(self.content[(begin+1):end]))
            self.dataset[i].readfp(filepart)


class CheckReferenceWaveforms(unittest.TestCase):
    def waveformgenerator(self, domain, arg):
        func = {'TD': lalsim.SimInspiralChooseTDWaveform,
                'FD': lalsim.SimInspiralChooseFDWaveform}[domain]
        return func(*arg)

    paramnames = {'TD': ['phiref', 'deltaT', 'm1', 'm2', 'spin1x', 'spin1y', 'spin1z',
                         'spin2x', 'spin2y', 'spin2z', 'fmin', 'fref', 'distance', 'inclination',
                         'lambda1', 'lambda2', 'waveformFlags', 'nonGRparams',
                         'ampOrder', 'phaseOrder'],
                  'FD': ['phiref', 'deltaF', 'm1', 'm2', 'spin1x', 'spin1y', 'spin1z',
                         'spin2x', 'spin2y', 'spin2z', 'fmin', 'fmax', 'fref', 'distance',
                         'inclination', 'lambda1', 'lambda2', 'waveformFlags',
                         'nonGRparams', 'ampOrder', 'phaseOrder']}

    paramtype = {'phiref':float, 'deltaT':float, 'deltaF':float,
                 'm1':float, 'm2':float,
                 'spin1x':float, 'spin1y':float, 'spin1z':float,
                 'spin2x':float, 'spin2y':float, 'spin2z':float,
                 'fmin':float, 'fref':float, 'distance':float, 'fmax':float,
                 'inclination':float, 'lambda1':float, 'lambda2':float,
                 'waveformFlags':lambda x: x, 'nonGRparams':lambda x: x,
                 'ampOrder':int, 'phaseOrder':int}
    #TODO: introduce function that properly handles waveformFlags and nonGRparams

    def errmsg(self, obj, approxstr, par):
        return ('{1} fails consistency test of {0} for the following '
        + 'parameters:\n{2}').format(obj, approxstr, par)



if __name__ == '__main__':
    if options.reffilename == DEFAULT_FILE:
        filepath = (os.path.dirname(os.path.abspath(__file__)))
        absfilename = filepath + '/' + options.reffilename
    else:
        absfilename = options.reffilename
    reffile = ReferenceFile(absfilename)
    allapprox = list(set([conf.get('approximant', 'approximant')
            for conf in reffile.dataset]))

    if options.approx == 'all':
        datasets = [[conf for conf in reffile.dataset
                if conf.get('approximant', 'approximant') == approx]
                for approx in allapprox]
        for approx, dataset in zip(allapprox, datasets):
            addApproxTestToClass(approx, dataset)
    else:
        dataset = [conf for conf in reffile.dataset
                if conf.get('approximant', 'approximant') == options.approx]
        addApproxTestToClass(options.approx, dataset)

    suite = unittest.TestLoader().loadTestsFromTestCase(CheckReferenceWaveforms)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if result.wasSuccessful():
        sys.exit(0)
    else:
        sys.exit(1)

