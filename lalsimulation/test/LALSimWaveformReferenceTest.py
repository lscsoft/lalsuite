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
parser.add_option('-a', '--approximant', action = 'store', type = 'string',
        dest = 'approx', default = 'all',
        help = 'waveform approximant [default: %default]')
parser.add_option('-p', '--plot', action = 'store_true', dest = 'plot',
                  default = False, help = 'save debugging plots if tests ' + 
                  'fail [default: %default] WORKS FOR TIME-DOMAIN APPROXIMANTS ONLY!')
parser.add_option('-s', '--separate', action = 'store_true', dest = 'sep',
                  default = False, help = 'seperate all data sets into ' +
                  'individual test, even if they belong to the same ' + 
                  'approximant [default: %default]')

(options, args) = parser.parse_args()


if options.plot:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

def waveformgenerator(domain, arg):
    func = {'TD': lalsim.SimInspiralChooseTDWaveform,
            'FD': lalsim.SimInspiralChooseFDWaveform}[domain]
    return func(*arg)
    
def recalculate_waveform(conf):
    '''Reads single data set and returns both the reference waveform and the newly calculated
    waveform with the same parameters.
    Returns: [hpref, hcref, hpnew, hcnew]
    '''
    approx = lalsim.GetApproximantFromString(conf.get('approximant',
            'approximant'))
    domain = conf.get('approximant', 'domain')
    if domain=='TD':
        hpref, hcref = [np.array(map(float, (conf.get('waveform-data',
                l)).split())) for l in ['hp', 'hc']]
    if domain=='FD':
        hpRref, hpIref, hcRref, hcIref = [np.array(map(float,
                (conf.get('waveform-data', l)).split()))
                for l in ['hp_real', 'hp_imag', 'hc_real', 'hc_imag']]
        hpref = hpRref + 1j * hpIref
        hcref = hcRref + 1j * hcIref
                                                                                                               
    names = CheckReferenceWaveforms.paramnames[domain]
    parDict = dict([ (p, CheckReferenceWaveforms.paramtype[p](conf.get('parameters', p)) ) for p in conf.options('parameters') ])
    parDict['m1'] *= lal.MSUN_SI
    parDict['m2'] *= lal.MSUN_SI
    parDict['distance'] *= (1.e6 * lal.PC_SI)
                                                                                                               
    params = [parDict[name] for name in names]
    params.append(approx)
                                                                                                               
    hp, hc = waveformgenerator(domain, params)
    
    return [hpref, hcref, hp, hc]

if options.plot:
    def waveformplots(labels, xnew, ynew, xref=None, yref=None, name='', counter=0, \
                      plot_func = plt.plot):
        '''create plots of input arrays'''
        plt.figure()
        plt.title(' / '.join(labels))
        if yref is not None:
            plot_func(xnew, ynew, 'b', label='current state')
            plot_func(xref, yref, 'r--', label='reference')
            plt.legend(loc='lower left')
        else:
            plot_func(xnew, ynew, 'b')

        plt.xlabel('time [s] or frequency [Hz]')
        plt.ylabel(name)
        plt.savefig(labels[0] + '_' + str(counter).zfill(2) + '_' + name + '.png')
        plt.clf()


def generateAttributes(datasets, counter = ''):
    '''Generates test and plot methods to be added to the CheckReferenceWaveforms class.
    The input is a list of datasets extracted from the reference file.'''
    waveforms = [recalculate_waveform(conf) for conf in datasets]
        
    def test_approx(self):
        '''check for consistent waveform polarisations'''
        for conf, wfs in zip(datasets, waveforms):
            hpref, hcref, hp, hc = wfs
            domain = conf.get('approximant', 'domain')
            approxstr = conf.get('approximant', 'approximant') + ' / ' + domain
            epochref = conf.getfloat('waveform-data', 'epoch')

            names = self.paramnames[domain]
            parstring =' / '.join([name + ': '
                    + str(conf.get('parameters', name)) for name in names])
            epoch = float(hp.epoch)

            # Actual test starts here
            self.assertTrue(np.allclose(epochref, epoch),
                    self.errmsg('epoch', approxstr, parstring))
            self.assertEqual(hp.data.data.size, hpref.size,
                             self.errmsg('length of generated hplus array',
                             approxstr, parstring))
            self.assertEqual(hc.data.data.size, hcref.size,
                             self.errmsg('length of generated hcross array',
                             approxstr, parstring))
            hpmean = np.abs(hpref).mean()
            hcmean = np.abs(hcref).mean()
            self.assertTrue(np.allclose(hp.data.data / hpmean, hpref / hpmean),
                            self.errmsg('hplus', approxstr, parstring))
            self.assertTrue(np.allclose(hc.data.data / hcmean, hcref / hcmean),
                            self.errmsg('hcross', approxstr, parstring))
    if options.plot:                        
        def plot_approx(self):
            '''plot the amplitude and phase'''
            i = 1
            for conf, wfs in zip(datasets, waveforms):
                hpref, hcref, hp, hc = wfs
                approx = conf.get('approximant', 'approximant')
                domain = conf.get('approximant', 'domain')
                m1, m2, s1z, s2z = ['%.2g' % conf.getfloat('parameters', name)
                                for name in ['m1', 'm2', 'spin1z', 'spin2z']]
                iota = conf.getfloat('parameters', 'inclination')
                cfac = np.cos(iota)
                pfac = 0.5 * (1. + cfac*cfac);

                href = hpref / pfac + 1.j * hcref / cfac
                h = hp.data.data / pfac + 1.j * hc.data.data / cfac
                steppar = {'TD': 'deltaT', 'FD': 'deltaF'}
                dX = conf.getfloat('parameters', steppar[domain])
                xvals, xvals_ref = [np.linspace(0., dX * (x.size - 1), x.size) \
                                   for x in [h, href]]
                
                if counter is not '':
                    num = counter
                else: 
                    num = i

                if domain=='TD':
                    xvals_ref += conf.getfloat('waveform-data','epoch')
                    xvals += hp.epoch
                else:
                    #remove zero frequency entry
                    xvals = xvals[1:]
                    xvals_ref = xvals_ref[1:]
                    h = h[1:]
                    href = href[1:]

                plotfunc = {'TD': plt.plot, 'FD': plt.loglog}
                waveformplots([approx, m1, m2, s1z, s2z], xvals, np.abs(h), \
                              xvals_ref, np.abs(href), 'amplitude', num, \
                              plotfunc[domain])
                plotfunc['FD'] = plt.semilogx
                waveformplots([approx, m1, m2, s1z, s2z], xvals, np.unwrap(np.angle(h)), \
                              xvals_ref, np.unwrap(np.angle(href)), 'phase', num, \
                              plotfunc[domain])
                if np.allclose(xvals, xvals_ref, atol = 1e-6):
                    waveformplots([approx, m1, m2, s1z, s2z], xvals, \
                      np.unwrap(np.angle(h)) - np.unwrap(np.angle(href)), \
                              name = 'phase_diff', counter = num, \
                              plot_func = plotfunc[domain])
                    sel = (np.abs(href) > 0)
                    waveformplots([approx, m1, m2, s1z, s2z], xvals[[sel]], \
                      np.abs(h[[sel]]) / np.abs(href[[sel]]), name = 'amp_quot', \
                      counter = num, plot_func = plotfunc[domain])
                
                i += 1
                
        return test_approx, plot_approx

    else:
        return test_approx

def addApproxTestToClass(approx, dataset, counter=''):
    '''adds test method as "test_approx" (where approx is the approximant name)'''
    if options.plot:
        test_method, plot_method = generateAttributes(dataset, counter)
        plot_method.__name__ = 'plot_' + approx + counter
        plot_method.__doc__ = plot_method.__doc__ + ' of ' + approx
        setattr(CheckReferenceWaveforms, plot_method.__name__, plot_method)
    else:
        test_method = generateAttributes(dataset, counter)
        
    test_method.__name__ = 'test_' + approx + counter
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
        # possibly add more above
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
    
    if options.approx == 'all' and not options.sep:
        datasets = [[conf for conf in reffile.dataset
                if conf.get('approximant', 'approximant') == approx]
                for approx in allapprox]
        for approx, dataset in zip(allapprox, datasets):
            addApproxTestToClass(approx, dataset)
    elif options.approx == 'all' and options.sep:
        i=1
        for conf in reffile.dataset:
            addApproxTestToClass(conf.get('approximant', 'approximant'), [conf],\
            str(i).zfill(2))
            i += 1
    else:
        dataset = [conf for conf in reffile.dataset
                if conf.get('approximant', 'approximant') == options.approx]
        if options.sep:
            i = 1
            for conf in dataset:
                addApproxTestToClass(options.approx, [conf], str(i).zfill(2))
                i += 1
        else:
            addApproxTestToClass(options.approx, dataset)
        


    suite = unittest.TestLoader().loadTestsFromTestCase(CheckReferenceWaveforms)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if result.wasSuccessful():
        sys.exit(0)
    else:
        if options.plot:
            for failure in result.failures:
                testinst = failure[0]
                teststr = testinst.__dict__['_testMethodName']
                testinst.__getattribute__('plot' + teststr[4:])()
        sys.exit(1)
