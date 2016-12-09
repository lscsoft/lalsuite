#!/usr/bin/env python
# Copyright (C) 2014 Evan Ochsner, Frank Ohme
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
import copy
import unittest
from optparse import OptionParser, OptionGroup
import ConfigParser
from lalsimulation import git_version
import io

NEW_DATA_STR = '######### NEW DATASET #############\n'

# write all parameter information consistently
# [long name, type, default, help, metavar]
defstr = ' [default: %default]'
parameterInfo = [['domain', 'string', 'TD', 'Generate either a time (TD) or ' +
        'frequency (FD) domain waveform'+defstr , 'DOM'],
        ['approximant', 'string', 'TaylorT1', 'Approximant name'+defstr, 'APPROX'],
        ['phaseOrder', 'int', -1, 'Twice PN order of phase, [default -1 ' +
        'for highest available order]', 'ORD'],
        ['ampOrder', 'int', -1, 'Twice PN order of amplitude, [default -1 ' +
        'for highest available order]', 'ORD'],
        ['sampleRate', 'float', 16384., 'Sample rate in Hz for TD waveforms. ' +
        '(If specified different from default, SRATE will be used to ' +
        'calculate DELTAT)'+defstr, 'SRATE'],
        ['deltaT', 'float', 1./16384., 'Time steps in seconds for TD ' +
        'waveforms. (Will be overwritten if SRATE is defined.)'+defstr,
        'DELTAT'],
        ['deltaF', 'float', 1./8., 'Frequency bin size in Hz for FD waveforms' +
        defstr, 'DF'],
        ['m1', 'float', 10., 'Mass of first object in solar masses'+defstr,
        'M1'],
        ['m2', 'float', 1.4, 'Mass of second object in solar masses'+defstr,
        'M2'],
        ['distance', 'float', 100., 'Distance in Mpc'+defstr, 'DIST'],
        ['inclination', 'float', lal.PI/3., 'Inclination of binary to ' +
        'line of sight [default: PI/3]', 'INCL'],
        ['spin1x', 'float', 0., 'Vector components for spin of mass1 ' +
        '(default all 0)', 'S1X'],
        ['spin1y', 'float', 0., 'z-axis=line of sight, L in x-z plane at ' +
        'reference', 'S1Y'],
        ['spin1z', 'float', 0., 'Kerr limit: S1X^2 + S1Y^2 + S1Z^2 <= 1','S1Z'],
        ['spin2x', 'float', 0., 'Vector components for spin of mass2 ' +
        '(default all 0)', 'S2X'],
        ['spin2y', 'float', 0., 'z-axis=line of sight, L in x-z plane at ' +
        'reference', 'S2Y'],
        ['spin2z', 'float', 0., 'Kerr limit: S2X^2 + S2Y^2 + S2Z^2 <= 1','S2Z'],
        ['lambda1', 'float', 0., '(tidal deformability of mass 1) / ' +
        '(mass of body 1)^5 Reasonable values (~128-2560 for NS, 0 for BH)' +
        defstr, 'L1'],
        ['lambda2', 'float', 0., '(tidal deformability of mass 2) / ' +
        '(mass of body 2)^5 Reasonable values (~128-2560 for NS, 0 for BH)' +
        defstr, 'L2'],
        ['f_min', 'float', 40., 'Lowest GW frequency in Hz'+defstr, 'FMIN'],
        ['f_max', 'float', 0., 'Highest GW frequency in Hz - only used by FD ' +
        'waveforms [default 0 (as high as possible)]', 'FMAX'],
        ['f_ref', 'float', 0., 'Reference GW frequency in Hz at which phase, ' +
        'binary orientation are defined [default 0 (coalescence)]', 'FREF'],
        ['phiref', 'float', 1.2, 'Orbital phase at FREF'+defstr, 'PHIREF'],
        ['longAscNodes', 'float', 0., 'Longitude of ascending nodes', 'LONGASCNODES'],
        ['eccentricity', 'float', 0., 'eccentricity', 'ECCENTRICITY'],
        ['meanPerAno', 'float', 0., 'Mean periastron anomaly', 'MEANPERANO'],
        ['eps', 'float', 1.e-2, 'Tolerance check', 'TOLERANCE']]

paramnames = {'TD': ['m1', 'm2', 'spin1x', 'spin1y', 'spin1z', 'spin2x', 'spin2y', 'spin2z',
                     'distance', 'inclination', 'phiref', 'longAscNodes', 'eccentricity', 'meanPerAno',
                     'deltaT', 'f_min', 'f_ref'],
              'FD': ['m1', 'm2', 'spin1x', 'spin1y', 'spin1z', 'spin2x', 'spin2y', 'spin2z',
                         'distance', 'inclination', 'phiref', 'longAscNodes', 'eccentricity', 'meanPerAno',
                         'deltaF', 'f_min', 'f_max', 'f_ref']}

parameterDicts = [{'name': p[0], 'type': p[1], 'default': p[2], 'help': p[3],
        'metavar': p[4]} for p in parameterInfo]

defaultDict = {p['name']: p['default'] for p in parameterDicts}
defSampleRate = defaultDict['sampleRate']  # needed later
partypes = {p['name']: p['type'] for p in parameterDicts}

usage = ('usage: %prog [options]\n' +
    'Generate a file containing a waveform (h+, hx), the parameters used' +
    ' to generate it, and the git version of the code used. This is' +
    ' useful for generating reference versions of reviewed waveforms.')

# Input/output file related options
parser = OptionParser(usage=usage)
parser.add_option('-o', '--output', metavar='FILE',
        default='waveformreference.dat',
        help='Write/append output to FILE [default: %default]')
parser.add_option('-i', '--input', metavar='FILE', default=None,
        dest = 'infile', help='Specify waveform parameters in FILE; if ' +
        'parameters are specified in both FILE and by a command-line option, ' +
        'the latter will be used.')

# Approximant options
group = OptionGroup(parser, 'Approximant details', 'Options that can be ' +
        'specified in input file under section [approximant]')

typeget = {'string': lambda x: x, 'float': float, 'int': int}

for p in parameterDicts[:2]:
    group.add_option('--' + p['name'],
            default = typeget[p['type']](p['default']),
            metavar = p['metavar'], type = p['type'], help = p['help'])
parser.add_option_group(group)

# Waveform parameters
group = OptionGroup(parser, 'Waveform parameters', 'Options that can be ' +
        'specified in input file under section [parameters]')
for p in parameterDicts[2:]:
    group.add_option('--' + p['name'],
            default = typeget[p['type']](p['default']),
            metavar = p['metavar'], type = p['type'], help = p['help'])
parser.add_option_group(group)

# Actually reading the options
(opts, args) = parser.parse_args()
# Sanity check args related to sample rate
if opts.sampleRate != defaultDict['sampleRate']\
        and opts.deltaT != defaultDict['deltaT']:
    raise ValueError('Specified both --sampleRate and --deltaT. ' +
            'Provide one or the other.')


class InputFile:
    '''When initialized with the filename for the input file,
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
        ConfigParser.RawConfigParser.optionxform = str
        # prevent ConfigParser to use lower case version of option
        self.dataset = [ConfigParser.RawConfigParser(defaults = defaultDict)
                for i in range(len(self.newapproxindex))]
        self.newapproxindex.append(self.size)
        for i in range(len(self.newapproxindex) - 1):
            begin, end = self.newapproxindex[i:(i+2)]
            filepart = io.BytesIO('\n'.join(self.content[(begin+1):end]))
            self.dataset[i].readfp(filepart)


def updateOptions(config):
    '''Re-define defaults according to dictionary in config [ConfigParser
    object] and read options again'''
    typeget = {'string': config.get, 'float': config.getfloat, 'int': config.getint}
    pardict1 = {p['name']: typeget[p['type']]('approximant', p['name'])
            for p in parameterDicts[:2]}
    pardict2 = {p['name']: typeget[p['type']]('parameters', p['name'])
            for p in parameterDicts[2:]}
    defaultDict = dict(pardict1.items() + pardict2.items())
    parser.set_defaults(**defaultDict)
    (opts, args) = parser.parse_args()
    return opts


def generateWfparameters(opts):
    '''Turn parameters specified by input file and/or command line options into
    parameter dictionary that can be used for waveform generation.'''
    optsDict = vars(opts)
    inputpar = copy.deepcopy(optsDict)
    inputpar['m1'] *= lal.MSUN_SI
    inputpar['m2'] *= lal.MSUN_SI
    inputpar['distance'] *= (1.e6 * lal.PC_SI)
    inputpar['approximant'] = lalsim.GetApproximantFromString(
        inputpar['approximant'])

    if (inputpar['sampleRate'] != defSampleRate):
        inputpar['deltaT'] = 1./inputpar['sampleRate']
        optsDict['deltaT'] = inputpar['deltaT']

    domain = optsDict['domain']
    wfinputpar=[inputpar[name] for name in paramnames[domain]]
    LALpars=lal.CreateDict()
    if (inputpar['lambda1']):
        lalsim.SimInspiralWaveformParamsInsertTidalLambda1(LALpars,inputpar['lambda1'])
    if (inputpar['lambda2']):
        lalsim.SimInspiralWaveformParamsInsertTidalLambda2(LALpars,inputpar['lambda2'])
    if (inputpar['ampOrder']!=-1):
        lalsim.SimInspiralWaveformParamsInsertPNAmplitudeOrder(LALpars,inputpar['ampOrder'])
    if (inputpar['phaseOrder']!=-1):
        lalsim.SimInspiralWaveformParamsInsertPNPhaseOrder(LALpars,inputpar['phaseOrder'])
    wfinputpar.append(LALpars)
    wfinputpar.append(inputpar['approximant'])

    return domain, wfinputpar, optsDict


def generateWaveformData(domain,inputpar):
    '''Generate waveform data from input parameter dictionary'''

    waveformgenerator={'TD': lalsim.SimInspiralChooseTDWaveform,
                       'FD': lalsim.SimInspiralChooseFDWaveform}[domain]
    hp, hc = waveformgenerator(*inputpar)
    assert hp.epoch==hc.epoch

    return hp, hc


def writeData(filename, polarisations, optsDict):
    '''Write waveform data to file'''
    hp, hc = polarisations

    # Open output file
    fp = open(filename, 'a')
    fp.write(NEW_DATA_STR)

    # Write metadata sections of output file
    fp.write('[approximant]\n')
    for p in ['approximant', 'domain']:
        fp.write(p + ' = ' + optsDict[p] + '\n')

    fp.write('\n[auxiliary]\n')
    fp.write('LALSuite-git = ' + git_version.version + '\n')
    if optsDict['approximant']=="SEOBNRv3":
        fp.write('EPS = ' + "{:.1e}".format(optsDict['eps']) + '\n')

    fp.write('\n[parameters]\n')

    formattext = {'float': lambda val: '%.16e' % val, 'string': lambda x: x,
                   'int': lambda val: str(val)}

    # collate parameters that are either (1) standard waveform
    # interface parameters, or (2) non-standard additional parameters
    # (e.g., lambda1, ampOrder, ...)
    outputpars = set(paramnames[optsDict['domain']]).intersection(
            set(optsDict.keys()))
    remainingpars = (set(optsDict.keys()) - outputpars).intersection(
            set(defaultDict.keys()) - {'approximant', 'domain'})
    non_std_pars = []
    for par in remainingpars:
        if optsDict[par] != defaultDict[par]:
            non_std_pars.append(par)
    outputpars = list(outputpars) + non_std_pars

    for p in outputpars:
        line = p + ' = ' + formattext[partypes[p]](optsDict[p]) + '\n'
        fp.write(line)

    # Write the waveform data section
    fp.write('\n[waveform-data]\n')
    line = 'epoch = %.16e\n' % hp.epoch; fp.write(line)
    if optsDict['domain']=='TD':
        fp.write('hp =')
        for x in hp.data.data:
            fp.write(' %.16e' % x)
        fp.write('\nhc =')
        for x in hc.data.data:
            fp.write(' %.16e' % x)
    elif optsDict['domain']=='FD':
        hpreal = np.real(hp.data.data)
        hpimag = np.imag(hp.data.data)
        hcreal = np.real(hc.data.data)
        hcimag = np.imag(hc.data.data)
        fp.write('hp_real =')
        for x in hpreal:
            fp.write(' %.16e' % x)
        fp.write('\nhp_imag =')
        for x in hpimag:
            fp.write(' %.16e' % x)
        fp.write('\nhc_real =')
        for x in hcreal:
            fp.write(' %.16e' % x)
        fp.write('\nhc_imag =')
        for x in hcimag:
            fp.write(' %.16e' % x)
    fp.write('\n\n\n')

    fp.close()

    return


if __name__ == '__main__':
    if opts.infile:
        infile = InputFile(opts.infile)
        # if parameters are specified by input file, loop through all
        # datasets of input file, set default accordingly and generate wfs.
        # NOTE: command line options overwrite parameters in input file
        # for every dataset
        for conf in infile.dataset:
            domain, par, optsdict = generateWfparameters(updateOptions(conf))
            polarisations = generateWaveformData(domain, par)
            writeData(opts.output, polarisations, optsdict)
    else:
        domain, par, optsdict = generateWfparameters(opts)
        polarisations = generateWaveformData(domain, par)
        writeData(opts.output, polarisations, optsdict)
