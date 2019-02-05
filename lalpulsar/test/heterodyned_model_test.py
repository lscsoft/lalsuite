#!/usr/bin/env python

# Copyright (C) 2019  Matthew Pitkin
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

"""
Test script for the heterodyned pulsar model functions.

These functions are tested against fiducial outputs from the reviewed
lalapps_pulsar_parameter_estimation_nested code.
"""

from __future__ import print_function, division

import os
import sys
import numpy as np
import lal
import lalpulsar
from lalpulsar.PulsarParametersWrapper import PulsarParametersPy

"""
The first test output is to create a signal model for a source assuming that
the heterodyne precisely matches the signal, i.e., it only varies due to the
antenna pattern. lalapps_pulsar_parameter_estimation_nested has been run with
the following pulsar parameter "inj.par" file:

PSRJ     TEST
RAJ      01:23:34.5
DECJ     -45:01:23.4
F0       123.456789
F1       -9.87654321e-12
PEPOCH   58000
H0       5.6e-26
COSIOTA  -0.2
PSI      0.4
PHI0     2.3
EPHEM    DE405
UNITS    TCB

lalapps_pulsar_parameter_estimation_nested --par-file inj.par --outfile test.hdf --harmonics 2 --inject-file inj.par --inject-only --inject-output inj.txt --fake-data H1,L1 --fake-starts 1000000000,1000000000 --fake-lengths 86400,86400 --fake-dt 3600,3600

The outputs of "inj.txt_H1_2.0_signal_only" and "inj.txt_L1_2.0_signal_only"
are given below.
"""

t1output = {'H1': np.array([[1000000000.0, -4.365015078242e-27, -4.216860020522e-27],
                            [1000003600.0, -3.123840987225e-27, -7.303345528632e-27],
                            [1000007200.0, -1.392894194743e-27, -8.343779137631e-27],
                            [1000010800.0, 3.907780411954e-28, -7.371234643743e-27],
                            [1000014400.0, 1.803192456526e-27, -4.935204629063e-27],
                            [1000018000.0, 2.543802234234e-27, -1.933962480155e-27],
                            [1000021600.0, 2.510964710484e-27, 6.441088928009e-28],
                            [1000025200.0, 1.822605255722e-27, 1.999009861526e-27],
                            [1000028800.0, 7.770155400880e-28, 1.741944189379e-27],
                            [1000032400.0, -2.352775438513e-28, 1.624267692628e-30],
                            [1000036000.0, -8.441413982604e-28, -2.614474888096e-27],
                            [1000039600.0, -8.062469575161e-28, -5.193209436791e-27],
                            [1000043200.0, -7.603746349482e-29, -6.776100558983e-27],
                            [1000046800.0, 1.178187576640e-27, -6.635532319641e-27],
                            [1000050400.0, 2.617633050484e-27, -4.491492229213e-27],
                            [1000054000.0, 3.824370674045e-27, -6.087150052676e-28],
                            [1000057600.0, 4.415959728830e-27, 4.253199979849e-27],
                            [1000061200.0, 4.152430242107e-27, 9.024733948035e-27],
                            [1000064800.0, 3.006434304359e-27, 1.259817232053e-26],
                            [1000068400.0, 1.177337218278e-27, 1.411375249637e-26],
                            [1000072000.0, -9.549597694961e-28, 1.318431214347e-26],
                            [1000075600.0, -2.924662644700e-27, 9.998085875592e-27],
                            [1000079200.0, -4.298094465681e-27, 5.272226410224e-27],
                            [1000082800.0, -4.783899643376e-27, 6.937006559007e-29]]),
            'L1': np.array([[1000000000.0, 2.354774544244e-27, 2.989220110321e-27],
                            [1000003600.0, 1.102264220569e-27, 4.794624473630e-27],
                            [1000007200.0, -3.784783044226e-28, 4.901350258365e-27],
                            [1000010800.0, -1.726096919602e-27, 3.558770900936e-27],
                            [1000014400.0, -2.638664308831e-27, 1.375241458539e-27],
                            [1000018000.0, -2.950862387455e-27, -8.626400678749e-28],
                            [1000021600.0, -2.672975553462e-27, -2.415588636184e-27],
                            [1000025200.0, -1.981652053751e-27, -2.800065246359e-27],
                            [1000028800.0, -1.165367276485e-27, -1.923079648198e-27],
                            [1000032400.0, -5.396781233382e-28, -1.064506005032e-28],
                            [1000036000.0, -3.556885019889e-28, 2.005730381569e-27],
                            [1000039600.0, -7.267224431594e-28, 3.631191457942e-27],
                            [1000043200.0, -1.593447257663e-27, 4.074861141635e-27],
                            [1000046800.0, -2.737073733735e-27, 2.933672841813e-27],
                            [1000050400.0, -3.837271308977e-27, 2.246600367826e-28],
                            [1000054000.0, -4.559321270321e-27, -3.599805614295e-27],
                            [1000057600.0, -4.646886209647e-27, -7.754852546804e-27],
                            [1000061200.0, -3.995330812549e-27, -1.131770424887e-26],
                            [1000064800.0, -2.685441479913e-27, -1.346283635637e-26],
                            [1000068400.0, -9.681500316052e-28, -1.367509409753e-26],
                            [1000072000.0, 7.960559086451e-28, -1.188425812350e-26],
                            [1000075600.0, 2.227270382987e-27, -8.484795106371e-27],
                            [1000079200.0, 3.021952545474e-27, -4.235613315750e-27],
                            [1000082800.0, 3.029100989189e-27, -6.642344227215e-29]])}


"""
The second test output is to create a signal model for a source assuming that
the heterodyne and injected signal do not precisely matches. 
lalapps_pulsar_parameter_estimation_nested has been run with the following
pulsar parameter "inj.par" file:

PSRJ     TEST
RAJ      01:23:34.5
DECJ     -45:01:23.4
F0       123.456789
F1       -9.87654321e-12
PEPOCH   58000
H0       5.6e-26
COSIOTA  -0.2
PSI      0.4
PHI0     2.3
EPHEM    DE405
UNITS    TCB

and "heterodyne" "het.par" file:

PSRJ     TEST
RAJ      01:23:34.6
DECJ     -45:01:23.5
F0       123.4567
F1       -9.876e-12
PEPOCH   58000
H0       5.6e-26
COSIOTA  -0.2
PSI      0.4
PHI0     2.3
EPHEM    DE405
UNITS    TCB

lalapps_pulsar_parameter_estimation_nested --par-file het.par --outfile test.hdf --harmonics 2 --inject-file inj.par --inject-only --inject-output inj.txt --fake-data H1 --fake-starts 1000000000 --fake-lengths 600 --fake-dt 60

The output of "inj.txt_H1_2.0_signal_only" is given below.
"""

t2output = np.array([[1000000000.0, -5.941104012924e-27, 1.561111119568e-27],
                     [1000000060.0, -6.078674025313e-27, 1.108663867367e-27],
                     [1000000120.0, -6.181603080801e-27, 6.451136978413e-28],
                     [1000000180.0, -6.249044892549e-27, 1.731192817987e-28],
                     [1000000240.0, -6.280357307660e-27, -3.045989420574e-28],
                     [1000000300.0, -6.275106653741e-27, -7.852744278380e-28],
                     [1000000360.0, -6.233070880130e-27, -1.266110481373e-27],
                     [1000000420.0, -6.154241482804e-27, -1.744296420398e-27],
                     [1000000480.0, -6.038824191671e-27, -2.217023814819e-27],
                     [1000000540.0, -5.887238436697e-27, -2.681502811182e-27]])


"""
The third test output is to create a signal model for a source assuming that
the heterodyne and injected signal do not precisely matches, but with emission
at both once and twice the rotation frequency.
lalapps_pulsar_parameter_estimation_nested has been run with
the following pulsar parameter "inj.par" file:

PSRJ     TEST
RAJ      01:23:34.5
DECJ     -45:01:23.4
F0       123.456789
F1       -9.87654321e-12
PEPOCH   58000
C22      5.6e-26
C21      1.4e-25
COSIOTA  -0.2
PSI      0.4
PHI21    2.3
PHI22    4.5
EPHEM    DE405
UNITS    TCB

and "heterodyne" "het.par" file:

PSRJ     TEST
RAJ      01:23:34.6
DECJ     -45:01:23.5
F0       123.4567
F1       -9.876e-12
PEPOCH   58000
C22      5.6e-26
C21      1.4e-25
COSIOTA  -0.2
PSI      0.4
PHI21    2.3
PHI22    4.5
EPHEM    DE405
UNITS    TCB

lalapps_pulsar_parameter_estimation_nested --par-file het.par --outfile test.hdf --harmonics 1,2 --inject-file inj.par --inject-only --inject-output inj.txt --fake-data H1 --fake-starts 1000000000,1000000000 --fake-lengths 600,600 --fake-dt 60,60

The outputs of "inj.txt_H1_1.0_signal_only" and "inj.txt_H1_2.0_signal_only"
are given below.
"""

t3output = {'1': np.array([[1000000000.0, -2.281620347081e-26, -7.267059337744e-27],
                           [1000000060.0, -2.242012648444e-26, -8.024920903732e-27],
                           [1000000120.0, -2.199802615415e-26, -8.763862546799e-27],
                           [1000000180.0, -2.155075521735e-26, -9.482964910868e-27],
                           [1000000240.0, -2.107919965732e-26, -1.018134596439e-26],
                           [1000000300.0, -2.058427708343e-26, -1.085816221335e-26],
                           [1000000360.0, -2.006693507621e-26, -1.151260986472e-26],
                           [1000000420.0, -1.952814948068e-26, -1.214392590347e-26],
                           [1000000480.0, -1.896892272393e-26, -1.275138910868e-26],
                           [1000000540.0, -1.839028196029e-26, -1.333432099023e-26]]),
            '2': np.array([[1000000000.0, 1.151114436476e-26, -4.292865557392e-27],
                           [1000000060.0, 1.187524854552e-26, -3.419959925105e-27],
                           [1000000120.0, 1.217263381782e-26, -2.518042744682e-27],
                           [1000000180.0, 1.240108521541e-26, -1.592235817764e-27],
                           [1000000240.0, 1.255878166730e-26, -6.478246234004e-28],
                           [1000000300.0, 1.264430777434e-26, 3.097719790379e-28],
                           [1000000360.0, 1.265666324682e-26, 1.275032881006e-27],
                           [1000000420.0, 1.259526996162e-26, 2.242366499355e-27],
                           [1000000480.0, 1.245997657263e-26, 3.206142957363e-27],
                           [1000000540.0, 1.225106070777e-26, 4.160726677163e-27]])}


"""
The fourth test output is to create a signal model for a source assuming that
the heterodyne and injected signal do not precisely matches, adding in some
binary system parameters.
lalapps_pulsar_parameter_estimation_nested has been run with
the following pulsar parameter "inj.par" file:

PSRJ     TEST
RAJ      01:23:34.5
DECJ     -45:01:23.4
F0       123.456789
F1       -9.87654321e-12
PEPOCH   58000
H0       5.6e-26
COSIOTA  -0.2
PSI      0.4
PHI0     2.3
BINARY   BT
T0       58121.3
ECC      0.0001
OM       1.2
A1       8.9
PB       0.54
EPHEM    DE405
UNITS    TCB

and "heterodyne" "het.par" file:

PSRJ     TEST
RAJ      01:23:34.6
DECJ     -45:01:23.5
F0       123.4567
F1       -9.876e-12
PEPOCH   58000
H0       5.6e-26
COSIOTA  -0.2
PSI      0.4
PHI0     2.3
BINARY   BT
T0       58121.3
ECC      0.0001
OM       2.2
A1       8.9
PB       0.54
EPHEM    DE405
UNITS    TCB

lalapps_pulsar_parameter_estimation_nested --par-file het.par --outfile test.hdf --harmonics 2 --inject-file inj.par --inject-only --inject-output inj.txt --fake-data H1 --fake-starts 1000000000 --fake-lengths 600 --fake-dt 60

The output of "inj.txt_H1_2.0_signal_only" is given below.
"""

t4output = np.array([[1000000000.0, -2.005382688010e-27, -5.806222962878e-27],
                     [1000000060.0, 6.157890501156e-27, 5.097038871699e-28],
                     [1000000120.0, -2.949373852705e-27, 5.470793562975e-27],
                     [1000000180.0, -3.871712827206e-27, -4.908194386786e-27],
                     [1000000240.0, 6.071606930851e-27, -1.634398269891e-27],
                     [1000000300.0, -8.720670796956e-28, 6.263634603914e-27],
                     [1000000360.0, -5.468817487024e-27, -3.247498058720e-27],
                     [1000000420.0, 5.125788383325e-27, -3.826689389048e-27],
                     [1000000480.0, 1.601886962341e-27, 6.230292960469e-27],
                     [1000000540.0, -6.411302897381e-27, -8.632664101744e-28]])


"""
The fifth test output is to create a signal model for a source assuming that
the heterodyne precisely matches the signal, i.e., it only varies due to the
antenna pattern, but including vector and scalar polarisation modes.
lalapps_pulsar_parameter_estimation_nested has been run with
the following pulsar parameter "inj.par" file:

PSRJ       TEST
RAJ        01:23:34.5
DECJ       -45:01:23.4
F0         123.456789
F1         -9.87654321e-12
PEPOCH     58000
HPLUS      5.6e-26
HCROSS     1.3e-26
HVECTORX   1.4e-26
HVECTORY   2.3e-26
HSCALARB   4.5e-26
HSCALARL   3.1e-26
PHI0TENSOR 0.4
PSITENSOR  1.2
PHI0SCALAR 3.1
PSISCALAR  0.2
PHI0VECTOR 4.5
PSIVECTOR  2.4
EPHEM      DE405
UNITS      TCB

lalapps_pulsar_parameter_estimation_nested --par-file inj.par --outfile test.hdf --harmonics 2 --inject-file inj.par --inject-only --inject-output inj.txt --fake-data H1 --fake-starts 1000000000 --fake-lengths 86400 --fake-dt 3600 --nonGR --inject-nonGR

The output of "inj.txt_H1_2.0_signal_only" is given below.
"""

t5output = np.array([[1000000000.0, 2.362531579953e-26, 1.176813780671e-26],
                     [1000003600.0, 2.152194509402e-26, 1.336412201064e-26],
                     [1000007200.0, 1.516382224064e-26, 1.177246654454e-26],
                     [1000010800.0, 6.398562101048e-27, 7.456172010976e-27],
                     [1000014400.0, -2.399506614349e-27, 1.582538413351e-27],
                     [1000018000.0, -8.969877836658e-27, -4.291514941390e-27],
                     [1000021600.0, -1.176734778693e-26, -8.637246923757e-27],
                     [1000025200.0, -1.036272117124e-26, -1.036209884816e-26],
                     [1000028800.0, -5.536465938605e-27, -9.098830531150e-27],
                     [1000032400.0, 9.593356581661e-28, -5.297636085703e-27],
                     [1000036000.0, 6.903808455425e-27, -9.736075872324e-29],
                     [1000039600.0, 1.023786500728e-26, 4.989077888210e-27],
                     [1000043200.0, 9.650277837673e-27, 8.490030277625e-27],
                     [1000046800.0, 4.959079544959e-27, 9.378384683966e-27],
                     [1000050400.0, -2.815874543073e-27, 7.354154735066e-27],
                     [1000054000.0, -1.171113927621e-26, 2.929341510260e-27],
                     [1000057600.0, -1.933951129658e-26, -2.707704965036e-27],
                     [1000061200.0, -2.353253161094e-26, -8.012087385831e-27],
                     [1000064800.0, -2.293243052000e-26, -1.150043241736e-26],
                     [1000068400.0, -1.737338292421e-26, -1.215410941899e-26],
                     [1000072000.0, -7.948860127300e-27, -9.699169393802e-27],
                     [1000075600.0, 3.255444398910e-27, -4.686765876987e-27],
                     [1000079200.0, 1.369234859571e-26, 1.647785967195e-27],
                     [1000082800.0, 2.100706686852e-26, 7.709825036865e-27]])


# set ephemeris files
earthephem = os.path.join(os.environ['LAL_TEST_PKGDATADIR'], 'earth00-19-DE405.dat.gz')
sunephem = os.path.join(os.environ['LAL_TEST_PKGDATADIR'], 'sun00-40-DE405.dat.gz')
timefile = os.path.join(os.environ['LAL_TEST_PKGDATADIR'], 'te405_2000-2040.dat.gz')

# get ephemeris files
edat = lalpulsar.InitBarycenter(earthephem, sunephem)
tdat = lalpulsar.InitTimeCorrections(timefile)


def test_one():
    par = PulsarParametersPy()
    par['F'] = [123.456789, -9.87654321e-12]  # set frequency
    par['DELTAF'] = [0.0, 0.0]  # frequency difference
    par['RAJ'] = lal.TranslateHMStoRAD('01:23:34.5')  # set right ascension
    par['DECJ'] = lal.TranslateDMStoRAD('-45:01:23.4')  # set declination
    pepoch = lal.TranslateStringMJDTTtoGPS('58000')
    par['PEPOCH'] = pepoch.gpsSeconds + 1e-9*pepoch.gpsNanoSeconds
    par['H0'] = 5.6e-26
    par['COSIOTA'] = -0.2
    par['PSI'] = 0.4
    par['PHI0'] = 2.3

    freqfactor = 2.  # set frequency factor

    for det in t1output.keys():
        # convert into GPS times
        gpstimes = lalpulsar.CreateTimestampVector(len(t1output[det]))
        for i, time in enumerate(t1output[det][:,0]): 
            gpstimes.data[i] = lal.LIGOTimeGPS(time)

        detector = lalpulsar.GetSiteInfo(det)

        # set the response function look-up table
        dt = t1output[det][1,0] - t1output[det][0,0]  # time step
        resp = lalpulsar.DetResponseLookupTable(t1output[det][0,0],
                                                detector,
                                                par['RAJ'],
                                                par['DECJ'],
                                                2880,
                                                dt)

        # get the heterodyned file SSB delay
        hetSSBdelay = lalpulsar.HeterodynedPulsarGetSSBDelay(par.PulsarParameters(),
                                                             gpstimes,
                                                             detector,
                                                             edat,
                                                             tdat,
                                                             lalpulsar.TIMECORRECTION_TCB)

        fullsignal = lalpulsar.HeterodynedPulsarGetModel(par.PulsarParameters(),
                                                         freqfactor,
                                                         1,
                                                         0,
                                                         0,
                                                         gpstimes,
                                                         hetSSBdelay,
                                                         0,
                                                         None,
                                                         0,
                                                         resp,
                                                         edat,
                                                         tdat,
                                                         lalpulsar.TIMECORRECTION_TCB)

        # check output matches that from lalapps_pulsar_parameter_estimation_nested
        if np.any(np.abs(fullsignal.data.data.real - t1output[det][:,1]) > 1e-34):
            return False
        elif np.any(np.abs(fullsignal.data.data.imag - t1output[det][:,2]) > 1e-34):
            return False

    return True


def test_two():
    parhet = PulsarParametersPy()
    parhet['F'] = [123.4567, -9.876e-12]  # set frequency
    parhet['RAJ'] = lal.TranslateHMStoRAD('01:23:34.6')  # set right ascension
    parhet['DECJ'] = lal.TranslateDMStoRAD('-45:01:23.5')  # set declination
    pepoch = lal.TranslateStringMJDTTtoGPS('58000')
    parhet['PEPOCH'] = pepoch.gpsSeconds + 1e-9*pepoch.gpsNanoSeconds
    parhet['H0'] = 5.6e-26
    parhet['COSIOTA'] = -0.2
    parhet['PSI'] = 0.4
    parhet['PHI0'] = 2.3
    
    parinj = PulsarParametersPy()
    parinj['F'] = [123.456789, -9.87654321e-12]  # set frequency
    parinj['DELTAF'] = parinj['F'] - parhet['F']  # frequency difference
    parinj['RAJ'] = lal.TranslateHMStoRAD('01:23:34.5')  # set right ascension
    parinj['DECJ'] = lal.TranslateDMStoRAD('-45:01:23.4')  # set declination
    pepoch = lal.TranslateStringMJDTTtoGPS('58000')
    parinj['PEPOCH'] = pepoch.gpsSeconds + 1e-9*pepoch.gpsNanoSeconds
    parinj['H0'] = 5.6e-26
    parinj['COSIOTA'] = -0.2
    parinj['PSI'] = 0.4
    parinj['PHI0'] = 2.3

    freqfactor = 2.  # set frequency factor
    det = 'H1'  # the detector

    # convert into GPS times
    gpstimes = lalpulsar.CreateTimestampVector(len(t2output))
    for i, time in enumerate(t2output[:,0]): 
        gpstimes.data[i] = lal.LIGOTimeGPS(time)

    detector = lalpulsar.GetSiteInfo(det)

    # set the response function look-up table
    dt = t2output[1,0] - t2output[0,0]  # time step
    resp = lalpulsar.DetResponseLookupTable(t2output[0,0],
                                            detector,
                                            parhet['RAJ'],
                                            parhet['DECJ'],
                                            2880,
                                            dt)

    # get the heterodyned file SSB delay
    hetSSBdelay = lalpulsar.HeterodynedPulsarGetSSBDelay(parhet.PulsarParameters(),
                                                         gpstimes,
                                                         detector,
                                                         edat,
                                                         tdat,
                                                         lalpulsar.TIMECORRECTION_TCB)

    fullsignal = lalpulsar.HeterodynedPulsarGetModel(parinj.PulsarParameters(),
                                                     freqfactor,
                                                     1,
                                                     0,
                                                     0,
                                                     gpstimes,
                                                     hetSSBdelay,
                                                     1,
                                                     None,
                                                     0,
                                                     resp,
                                                     edat,
                                                     tdat,
                                                     lalpulsar.TIMECORRECTION_TCB)

    # check output matches that from lalapps_pulsar_parameter_estimation_nested
    if np.any(np.abs(fullsignal.data.data.real - t2output[:,1]) > 1e-34):
        return False
    elif np.any(np.abs(fullsignal.data.data.imag - t2output[:,2]) > 1e-34):
        return False
    else:
        return True


def test_three():
    parhet = PulsarParametersPy()
    parhet['F'] = [123.4567, -9.876e-12]  # set frequency
    parhet['RAJ'] = lal.TranslateHMStoRAD('01:23:34.6')  # set right ascension
    parhet['DECJ'] = lal.TranslateDMStoRAD('-45:01:23.5')  # set declination
    pepoch = lal.TranslateStringMJDTTtoGPS('58000')
    parhet['PEPOCH'] = pepoch.gpsSeconds + 1e-9*pepoch.gpsNanoSeconds
    parhet['C22'] = 5.6e-26
    parhet['C21'] = 1.4e-25
    parhet['COSIOTA'] = -0.2
    parhet['PSI'] = 0.4
    parhet['PHI21'] = 2.3
    parhet['PHI22'] = 4.5

    parinj = PulsarParametersPy()
    parinj['F'] = [123.456789, -9.87654321e-12]  # set frequency
    parinj['DELTAF'] = parinj['F'] - parhet['F']  # frequency difference
    parinj['RAJ'] = lal.TranslateHMStoRAD('01:23:34.5')  # set right ascension
    parinj['DECJ'] = lal.TranslateDMStoRAD('-45:01:23.4')  # set declination
    pepoch = lal.TranslateStringMJDTTtoGPS('58000')
    parinj['PEPOCH'] = pepoch.gpsSeconds + 1e-9*pepoch.gpsNanoSeconds
    parinj['C22'] = 5.6e-26
    parinj['C21'] = 1.4e-25
    parinj['COSIOTA'] = -0.2
    parinj['PSI'] = 0.4
    parinj['PHI21'] = 2.3
    parinj['PHI22'] = 4.5

    det = 'H1'  # the detector
    detector = lalpulsar.GetSiteInfo(det)

    for harmonic in t3output.keys():
        freqfactor = float(harmonic)  # set frequency factor

        # convert into GPS times
        gpstimes = lalpulsar.CreateTimestampVector(len(t3output[harmonic]))
        for i, time in enumerate(t3output[harmonic][:,0]): 
            gpstimes.data[i] = lal.LIGOTimeGPS(time)

        # set the response function look-up table
        dt = t3output[harmonic][1,0] - t3output[harmonic][0,0]  # time step
        resp = lalpulsar.DetResponseLookupTable(t3output[harmonic][0,0],
                                                detector,
                                                parhet['RAJ'],
                                                parhet['DECJ'],
                                                2880,
                                                dt)

        # get the heterodyned file SSB delay
        hetSSBdelay = lalpulsar.HeterodynedPulsarGetSSBDelay(parhet.PulsarParameters(),
                                                             gpstimes,
                                                             detector,
                                                             edat,
                                                             tdat,
                                                             lalpulsar.TIMECORRECTION_TCB)

        fullsignal = lalpulsar.HeterodynedPulsarGetModel(parinj.PulsarParameters(),
                                                         freqfactor,
                                                         1,
                                                         0,
                                                         0,
                                                         gpstimes,
                                                         hetSSBdelay,
                                                         1,
                                                         None,
                                                         0,
                                                         resp,
                                                         edat,
                                                         tdat,
                                                         lalpulsar.TIMECORRECTION_TCB)

        # check output matches that from lalapps_pulsar_parameter_estimation_nested
        if np.any(np.abs(fullsignal.data.data.real - t3output[harmonic][:,1]) > 1e-34):
            return False
        elif np.any(np.abs(fullsignal.data.data.imag - t3output[harmonic][:,2]) > 1e-34):
            return False

    return True


def test_four():
    parhet = PulsarParametersPy()
    parhet['F'] = [123.4567, -9.876e-12]  # set frequency
    parhet['RAJ'] = lal.TranslateHMStoRAD('01:23:34.6')  # set right ascension
    parhet['DECJ'] = lal.TranslateDMStoRAD('-45:01:23.5')  # set declination
    pepoch = lal.TranslateStringMJDTTtoGPS('58000')
    parhet['PEPOCH'] = pepoch.gpsSeconds + 1e-9*pepoch.gpsNanoSeconds
    parhet['H0'] = 5.6e-26
    parhet['COSIOTA'] = -0.2
    parhet['PSI'] = 0.4
    parhet['PHI0'] = 2.3
    parhet['BINARY'] = 'BT'
    T0 = lal.TranslateStringMJDTTtoGPS('58121.3')
    parhet['T0'] = T0.gpsSeconds + 1e-9*T0.gpsNanoSeconds
    parhet['OM'] = np.deg2rad(2.2)
    parhet['A1'] = 8.9
    parhet['PB'] = 0.54*86400.
    parhet['ECC'] = 0.0001
    
    parinj = PulsarParametersPy()
    parinj['F'] = [123.456789, -9.87654321e-12]  # set frequency
    parinj['DELTAF'] = parinj['F'] - parhet['F']  # frequency difference
    parinj['RAJ'] = lal.TranslateHMStoRAD('01:23:34.5')  # set right ascension
    parinj['DECJ'] = lal.TranslateDMStoRAD('-45:01:23.4')  # set declination
    pepoch = lal.TranslateStringMJDTTtoGPS('58000')
    parinj['PEPOCH'] = pepoch.gpsSeconds + 1e-9*pepoch.gpsNanoSeconds
    parinj['H0'] = 5.6e-26
    parinj['COSIOTA'] = -0.2
    parinj['PSI'] = 0.4
    parinj['PHI0'] = 2.3
    parinj['BINARY'] = 'BT'
    T0 = lal.TranslateStringMJDTTtoGPS('58121.3')
    parinj['T0'] = T0.gpsSeconds + 1e-9*T0.gpsNanoSeconds
    parinj['OM'] = np.deg2rad(1.2)
    parinj['A1'] = 8.9
    parinj['PB'] = 0.54*86400.
    parinj['ECC'] = 0.0001

    freqfactor = 2.  # set frequency factor
    det = 'H1'  # the detector

    # convert into GPS times
    gpstimes = lalpulsar.CreateTimestampVector(len(t4output))
    for i, time in enumerate(t4output[:,0]): 
        gpstimes.data[i] = lal.LIGOTimeGPS(time)

    detector = lalpulsar.GetSiteInfo(det)

    # set the response function look-up table
    dt = t4output[1,0] - t4output[0,0]  # time step
    resp = lalpulsar.DetResponseLookupTable(t4output[0,0],
                                            detector,
                                            parhet['RAJ'],
                                            parhet['DECJ'],
                                            2880,
                                            dt)

    # get the heterodyned file SSB delay
    hetSSBdelay = lalpulsar.HeterodynedPulsarGetSSBDelay(parhet.PulsarParameters(),
                                                         gpstimes,
                                                         detector,
                                                         edat,
                                                         tdat,
                                                         lalpulsar.TIMECORRECTION_TCB)

    # get the heterodyned file BSB delay
    hetBSBdelay = lalpulsar.HeterodynedPulsarGetBSBDelay(parhet.PulsarParameters(),
                                                         gpstimes,
                                                         hetSSBdelay,
                                                         edat)

    fullsignal = lalpulsar.HeterodynedPulsarGetModel(parinj.PulsarParameters(),
                                                     freqfactor,
                                                     1,  # phase is varying between par files
                                                     0,  # not using ROQ
                                                     0,  # not using non-tensorial modes
                                                     gpstimes,
                                                     hetSSBdelay,
                                                     1,  # the SSB delay should be updated compared to hetSSBdelay
                                                     hetBSBdelay,
                                                     1,  # the BSB delay should be updated compared to hetBSBdelay
                                                     resp,
                                                     edat,
                                                     tdat,
                                                     lalpulsar.TIMECORRECTION_TCB)

    # check output matches that from lalapps_pulsar_parameter_estimation_nested
    if np.any(np.abs(fullsignal.data.data.real - t4output[:,1]) > 1e-33):
        return False
    elif np.any(np.abs(fullsignal.data.data.imag - t4output[:,2]) > 1e-33):
        return False
    else:
        return True


def test_five():
    par = PulsarParametersPy()
    par['F'] = [123.456789, -9.87654321e-12]  # set frequency
    par['DELTAF'] = [0.0, 0.0]  # frequency difference
    par['RAJ'] = lal.TranslateHMStoRAD('01:23:34.5')  # set right ascension
    par['DECJ'] = lal.TranslateDMStoRAD('-45:01:23.4')  # set declination
    pepoch = lal.TranslateStringMJDTTtoGPS('58000')
    par['PEPOCH'] = pepoch.gpsSeconds + 1e-9*pepoch.gpsNanoSeconds
    par['HPLUS'] = 5.6e-26
    par['HCROSS'] = 1.3e-26
    par['HVECTORX'] = 1.4e-26
    par['HVECTORY'] = 2.3e-26
    par['HSCALARB'] = 4.5e-26
    par['HSCALARL'] = 3.1e-26
    par['PHI0TENSOR'] = 0.4
    par['PSITENSOR'] = 1.2
    par['PHI0SCALAR'] = 3.1
    par['PSISCALAR'] = 0.2
    par['PHI0VECTOR'] = 4.5
    par['PSIVECTOR'] = 2.4
    par['PHI0'] = 2.3

    freqfactor = 2.  # set frequency factor

    det = 'H1'  # detector
    detector = lalpulsar.GetSiteInfo(det)

    gpstimes = lalpulsar.CreateTimestampVector(len(t5output))
    for i, time in enumerate(t5output[:,0]): 
        gpstimes.data[i] = lal.LIGOTimeGPS(time)


    # set the response function look-up table
    dt = t5output[1,0] - t5output[0,0]  # time step
    resp = lalpulsar.DetResponseLookupTable(t5output[0,0],
                                            detector,
                                            par['RAJ'],
                                            par['DECJ'],
                                            2880,
                                            dt)

    # get the heterodyned file SSB delay
    hetSSBdelay = lalpulsar.HeterodynedPulsarGetSSBDelay(par.PulsarParameters(),
                                                         gpstimes,
                                                         detector,
                                                         edat,
                                                         tdat,
                                                         lalpulsar.TIMECORRECTION_TCB)

    fullsignal = lalpulsar.HeterodynedPulsarGetModel(par.PulsarParameters(),
                                                     freqfactor,
                                                     1,
                                                     0,
                                                     1,  # use non-GR modes
                                                     gpstimes,
                                                     hetSSBdelay,
                                                     0,
                                                     None,
                                                     0,
                                                     resp,
                                                     edat,
                                                     tdat,
                                                     lalpulsar.TIMECORRECTION_TCB)


    # check output matches that from lalapps_pulsar_parameter_estimation_nested
    if np.any(np.abs(fullsignal.data.data.real - t5output[:,1]) > 1e-34):
        return False
    elif np.any(np.abs(fullsignal.data.data.imag - t5output[:,2]) > 1e-34):
        return False

    return True


# run tests
t1 = test_one()
t2 = test_two()
t3 = test_three()
t4 = test_four()
t5 = test_five()

if np.any(np.invert(np.array([t1, t2, t3, t4, t5]))):
    sys.exit(1)
else:
    sys.exit(0)
