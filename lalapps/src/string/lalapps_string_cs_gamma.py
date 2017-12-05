#
#  Copyright (C) 2007 Xavier Siemens
#  Copyright (C) 2010 Andrew Mergl
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
#/*            Cosmic string burst rate computation code for small loops          */
#/*                                                                               */
#/*                  Xavier Siemens, Jolien Creighton, Irit Maor                  */
#/*                                                                               */
#/*                         UWM/Caltech - September 2006                          */
#/*********************************************************************************/

#Port to Python from the C program cs_gamma.c
#Original C source code by Xavier Seimens
#Port by Andrew Mergl, June 2010

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
P = 1.0	# reconnection probability


def parse_command_line():
    parser = OptionParser(
        version = "Name: %%prog\n%s" % git_version.verbose_msg
    )
    parser.add_option("-a", "--frequency", type="float", help="Lowest frequency.")
    parser.add_option("-b", "--Gmustart", type="float", help="Lowest Gmu.")
    parser.add_option("-c", "--Gmuend", type="float", help="Largest Gmu.")
    parser.add_option("-d", "--nGmu", type="int", help="Nubmer of Gmu bins to do.")
    parser.add_option("-e", "--epsilonstart", type="float", help="Lowest epsilon")
    parser.add_option("-f", "--epsilonend", type="float", help="Largest epsilon.")
    parser.add_option("-g", "--nepsilon", type="int", help="Number of epsilon bins to do")
    parser.add_option("-i", "--index", type="float", help="Index for alpha as function of Gmu")
    parser.add_option("-m", "--efficiency-file", type="string", help="File with efficiency values and errors.")
    options, filenames = parser.parse_args()

    required_options = ["efficiency_file", "frequency", "Gmustart", "Gmuend", "nGmu", "epsilonstart", "epsilonend", "nepsilon", "index"]
    missing_options = [option for option in required_options if getattr(options, option) is None]
    if missing_options:
        raise ValueError("missing required option(s) %s" % ", ".join("--%s" % option.replace("_", "-") for option in missing_options))
    assert options.nGmu >= 2
    assert options.nepsilon >= 2

    return options, (filenames or [None])

#------------------------------------------------------------------------------
#              Main
#------------------------------------------------------------------------------


ops, files = parse_command_line()

#Open the efficiency file and read the three columns into three arrays
amp, eff, Deff = numpy.loadtxt(ops.efficiency_file, dtype="double", unpack=True)

#Open the output file and print the header
outfile = open("gamma.dat", 'w')
outfile.write('%     p           n           epsilon         Gmu       gammaAverage    gammaMin      gammaMax\n')
for i in range(ops.nepsilon):
    epsilon = math.exp(math.log(ops.epsilonstart) + i * (math.log(ops.epsilonend) - math.log(ops.epsilonstart)) / (ops.nepsilon - 1))
    for j in range(ops.nGmu):
        Gmu = math.exp(math.log(ops.Gmustart) + j * (math.log(ops.Gmuend) - math.log(ops.Gmustart)) / (ops.nGmu - 1))
        print >>sys.stderr, "%.1f%%: Gmu=%10.4g, epsilon=%10.4g, p=%4.2g\r" % (100.0 * (i * ops.nGmu + j) / (ops.nepsilon * ops.nGmu), Gmu, epsilon, P),

        alpha = epsilon * (LOOPS_RAD_POWER * Gmu)**ops.index

        zofA = cs_gamma.findzofA(Gmu, alpha, amp) 
        dRdz = cs_gamma.finddRdz(Gmu, alpha, ops.frequency, LOOPS_RAD_POWER, zofA)

        Dlnz = numpy.log(zofA[1:]) - numpy.log(zofA[:-1])

        gammaAverage = scipy.integrate.simps(eff[:-1] * zofA[:-1] * dRdz[:-1] * -Dlnz) * CUSPS_PER_LOOP / P
        gammaMin = scipy.integrate.simps(numpy.clip(eff[:-1] - Deff[:-1], 0.0, 1.0) * zofA[:-1] * dRdz[:-1] * -Dlnz) * CUSPS_PER_LOOP / P
        gammaMax = scipy.integrate.simps(numpy.clip(eff[:-1] + Deff[:-1], 0.0, 1.0) * zofA[:-1] * dRdz[:-1] * -Dlnz) * CUSPS_PER_LOOP / P

        outfile.write("%.17g  %.17g  %.17g  %.17g  %.17g  %.17g  %.17g\n" % (P, ops.index, epsilon, Gmu, gammaAverage, gammaMin, gammaMax))
print >>sys.stderr, "100.0%%: Gmu=%10.4g, epsilon=%10.4g, p=%4.2g" % (Gmu, epsilon, P)
