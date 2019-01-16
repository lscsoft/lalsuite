#
#  Copyright (C) 2007 Xavier Siemens
#  Copyright (C) 2018 Daichi Tsuna
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with with program; see the file COPYING. If not, write to the
#  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#  MA  02111-1307  USA
#

#/*********************************************************************************/
#/*            Cosmic string burst rate computation code for large loops          */
#/*                                                                               */
#/*                  Xavier Siemens, Jolien Creighton, Irit Maor                  */
#/*                                                                               */
#/*                         UWM/Caltech - September 2006                          */
#/*********************************************************************************/

#Port to Python from the C program cs_gammaLargeLoops.c
#Original C source code by Xavier Siemens
#Port by Daichi Tsuna, November 2018

from __future__ import print_function

import math
import numpy
from optparse import OptionParser
import scipy.integrate
import sys

from lalburst import cs_gamma
from lalburst import git_version

#Constants from cs_lambda_cosmo.h
LAMBDA_Z_EQ = 5440.0
LAMBDA_H_0 = 2.27e-18
LAMBDA_OMEGA_M = 0.279
LAMBDA_OMEGA_R = 8.5e-5
LOOPS_RAD_POWER = 50.0	# Gamma
CUSPS_PER_LOOP = 1.0	# c


def parse_command_line():
    parser = OptionParser(
        version = "Name: %%prog\n%s" % git_version.verbose_msg
    )
    parser.add_option("-a", "--frequency", type="float", help="Lowest frequency.")
    parser.add_option("-b", "--Gmustart", type="float", help="Lowest Gmu.")
    parser.add_option("-c", "--Gmuend", type="float", help="Largest Gmu.")
    parser.add_option("-d", "--nGmu", type="int", help="Number of Gmu bins to do.")
    parser.add_option("-e", "--pstart", type="float", help="Lowest p.")
    parser.add_option("-f", "--pend", type="float", help="Largest p.")
    parser.add_option("-g", "--np", type="int", help="Number of p bins to do.")
    parser.add_option("-m", "--model", metavar="Siemens06|Blanco-Pillado14|Ringeval07", type="string", help="Model of loop distribution. There are three models that can be taken in this code, the 'Siemens06' model, the 'Blanco-Pillado14' model, and the 'Ringeval07' model. See arXiv:1712.01168 for details on each model.")
    parser.add_option("-n", "--efficiency-file", type="string", help="File with efficiency values and errors.")
    options, filenames = parser.parse_args()

    required_options = ["efficiency_file", "frequency", "Gmustart", "Gmuend", "nGmu", "pstart", "pend", "np", "model"]
    missing_options = [option for option in required_options if getattr(options, option) is None]
    if missing_options:
        raise ValueError("missing required option(s) %s" % ", ".join("--%s" % option.replace("_", "-") for option in missing_options))
    if options.model not in ("Siemens06", "Blanco-Pillado14", "Ringeval07"):
    	raise ValueError("--model \"%s\" not recognized" % options.model)
    assert options.nGmu >= 2
    assert options.np >= 2

    return options, (filenames or [None])

#------------------------------------------------------------------------------
#              Main
#------------------------------------------------------------------------------


ops, files = parse_command_line()

#Open the efficiency file and read the three columns into three arrays
amp, eff, Deff = numpy.loadtxt(ops.efficiency_file, dtype="double", unpack=True)
dlnA = numpy.log(amp[1:]) - numpy.log(amp[:-1])

#Open the output file and print the header
outfile = open("gamma.dat", 'w')
outfile.write('%     p           Gmu       gammaAverage    gammaMin      gammaMax\n')

#Decide the redshift range to integrate over
lnz_min = -50.0
lnz_max = 50.0
dlnz = 0.1
z = numpy.logspace(lnz_min, lnz_max, int((lnz_max-lnz_min)/dlnz)+1, base = math.e)
dRdA = numpy.zeros(len(amp), dtype = float)

for i in range(ops.np):
	P = math.exp(math.log(ops.pstart) + i * (math.log(ops.pend) - math.log(ops.pstart)) / (ops.np - 1))
	for j in range(ops.nGmu):
		Gmu = math.exp(math.log(ops.Gmustart) + j * (math.log(ops.Gmuend) - math.log(ops.Gmustart)) / (ops.nGmu - 1))
		print("%.1f%%: Gmu=%10.4g, p=%4.2g\r" % (100.0 * (i * ops.nGmu + j) / (ops.np * ops.nGmu), Gmu, P), file=sys.stderr)

		dRdzdA = cs_gamma.finddRdzdA(Gmu, ops.frequency, LOOPS_RAD_POWER, amp, z, ops.model)

		#Integrate over z
		for k in range(len(amp)):
			dRdA[k] = scipy.integrate.simps(dRdzdA[k,:-1] * z[:-1] *  dlnz)
		#Integrate over A
		gammaAverage = scipy.integrate.simps(eff[:-1] * dRdA[:-1] * amp[:-1] * dlnA) * CUSPS_PER_LOOP / P
		gammaMin = scipy.integrate.simps(numpy.clip(eff[:-1] - Deff[:-1], 0.0, 1.0) * dRdA[:-1] * amp[:-1] * dlnA) * CUSPS_PER_LOOP / P
		gammaMax = scipy.integrate.simps(numpy.clip(eff[:-1] + Deff[:-1], 0.0, 1.0) * dRdA[:-1] * amp[:-1] * dlnA) * CUSPS_PER_LOOP / P

		outfile.write("%.17g  %.17g  %.17g  %.17g  %.17g\n" % (P, Gmu, gammaAverage, gammaMin, gammaMax))
print("100.0%%: Gmu=%10.4g, p=%4.2g" % (Gmu, P), file=sys.stderr)
