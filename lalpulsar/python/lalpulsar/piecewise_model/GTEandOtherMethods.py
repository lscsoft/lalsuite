# Copyright (C) 2019--2023 Benjamin Grace
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

import numpy as np
from . import BasisFunctions as bf
import matplotlib.pyplot as plt

# The general torque equation we have formulated
def gte(t, f0, ngte, kgte):
    return f0 * (1 + (ngte - 1) * f0 ** (ngte - 1) * kgte * t) ** (1 / (1 - ngte))

# Inverse of the GTE
def gteinv(freq, f0, ngte, kgte):
    return ((freq / f0) ** (1 - ngte) - 1) / ((ngte - 1) * f0 ** (ngte - 1) * kgte)

# Arbitrary derivative of the gte
def gtederivs(t, f0, ngte, kgte, deriv):
    if deriv == 0:
        return gte(t, f0, ngte, kgte)
    elif deriv == 1:
        return -kgte * (f0 ** ngte) * (1 + (ngte - 1) * kgte * (f0 ** (ngte - 1)) * t) ** (1 / (1 - ngte) - 1)
    elif deriv == 2:
        return ngte * (kgte ** 2) * (f0 ** (2 * ngte - 1)) * (1 + (ngte - 1) * kgte * (f0 ** (ngte - 1)) * t) ** (1 / (1 - ngte) - 2)

    # The below expression may not correct
    else:
        print("Derivatives calculated here have not been thoroughly checked!")
        prefactor = (kgte ** deriv) * (f0 ** (deriv * ngte - (deriv - 1)))

        factorialfactor = ngte - 1

        for i in range(deriv):
            factorialfactor *= ((1 / (1 - ngte)) - i)

        return prefactor * factorialfactor * (1 + (ngte - 1) * f0 ** (ngte - 1) * kgte * t) ** (1 / (ngte - 1) - deriv)

# Returns the minimum and maximum values of the GTE at a time t with initial frequency f0 where ngte and kgte can vary
# between the given ranges
def gteminmax(t, fmin, fmax, nmin, nmax, kmin, kmax):
    if t < 0:
        return [gte(t, fmin, nmin, kmin), gte(t, fmax, nmax, kmax)]
    else:
        return [gte(t, fmin, nmax, kmax), gte(t, fmax, nmin, kmin)]

# Returns the minimum and maximum values of the derivatives of the GTE given the parameter ranges
def gtederivminmax(t, fmin, fmax, nmin, nmax, kmin, kmax, deriv):
    if deriv == 0:
        return gteminmax(t, fmin, fmax, nmin, nmax, kmin, kmax)

    if deriv % 2 == 0:
        minderiv = gtederivs(t, fmin, nmin, kmin, deriv)
        maxderiv = gtederivs(t, fmax, nmax, kmax, deriv)
    else:
        minderiv = gtederivs(t, fmax, nmax, kmax, deriv)
        maxderiv = gtederivs(t, fmin, nmin, kmin, deriv)

    return [minderiv, maxderiv]

# Rounds number to certain number of significant figures
def roundsig(x, sig=4):
    return round(x, sig - int(np.floor(np.log10(abs(x)))) - 1)

# Returns the braking index range for the next knot given ngte and our tolerance for ngte between knots.
def nminmax(ngte, ntol, seglength):
    nmin = ngte * (1 - ntol * seglength)
    nmax = ngte  # * (1 + ntol * seglength)

    return [nmin, nmax]

# Written to give the braking index range on the previous knot given ngte and out tolerance on some knot. Should be written
# to be the inverse of the above method.
def nminmaxinverse(ngte, ntol, seglength):
    nmin = ngte  # / (1 + ntol * seglength)
    nmax = ngte / (1 - ntol * seglength)

    return [nmin, nmax]

# Returns the allowed braking index range if we start with an allowed braking index range at t = 0 of [n0min, n0max]
# and a range that we are allowed to evole into of [nmin, nmax]. Method based on the result of recursively using the
# nminmax method.
def nminmaxatknot(n0min, n0max, nmin, nmax, ntol, knotnum):
    thisnmin = n0min
    thisnmax = n0max

    for i in range(1, knotnum):
        seglength = bf.knotslist[i] - bf.knotslist[i - 1]
        thisnmin = nminmax(thisnmin, ntol, seglength)[0]
        thisnmax = nminmax(thisnmax, ntol, seglength)[1]

    if thisnmin < nmin:
        thisnmin = nmin
    if thisnmax > nmax:
        thisnmax = nmax

    return [thisnmin, thisnmax]

# Returns the allowed range for the braking index on the next knot if we have an initial range of n0min and n0max and
# the value of the braking index at knotnum - 1
def nextnminmax(ngte, n0min, n0max, nmin, nmax, ntol, knotnum, seglength):
    nextmin, nextmax = nminmax(ngte, ntol, seglength)

    allowedmin, allowedmax = nminmaxatknot(n0min, n0max, nmin, nmax, ntol, knotnum)

    if nextmin < allowedmin:
        nextmin = allowedmin
    if nextmax > allowedmax:
        nextmax = allowedmax

    return [nextmin, nextmax]

# Gives the next allowed ranges of the kgte on the next knot.
def kminmax(kgte, ktol, seglength):
    kmin = kgte * (1 - ktol * seglength)
    kmax = kgte  # * (1 + ktol * seglength)

    return [kmin, kmax]

# Gives the minimum and maximum allowed values of the value kgte at a given knot given its initial starting values. Written
# as applying the kminmax method recursively
def kminmaxatknot(k0min, k0max, kmin, kmax, ktol, knotnum):
    thiskmin = k0min
    thiskmax = k0max

    for i in range(1, knotnum):
        seglength = bf.knotslist[i] - bf.knotslist[i - 1]
        thiskmin = kminmax(thiskmin, ktol, seglength)[0]
        thiskmax = kminmax(thiskmax, ktol, seglength)[1]

    if thiskmin < kmin:
        thiskmin = kmin
    if thiskmax > kmax:
        thiskmax = kmax

    return [thiskmin, thiskmax]

# As the nextnminmax method but now for the constant kgte
def nextkminmax(kgte, k0min, k0max, kmin, kmax, ktol, knotnum, seglength):
    nextmin, nextmax = kminmax(kgte, ktol, seglength)

    allowedmin, allowedmax = kminmaxatknot(k0min, k0max, kmin, kmax, ktol, knotnum)

    if nextmin < allowedmin:
        nextmin = allowedmin
    if nextmax > allowedmax:
        nextmax = allowedmax

    return [nextmin, nextmax]

# Returns the time at which the GTE is at half its initial frequency value
def timeatwhichfreqishalved(f0, ngte, kgte):
    return roundsig((2 ** (ngte - 1) - 1) / ((ngte - 1) * f0 ** (ngte - 1) * kgte))

# TODO: Return the value of kgte for a purely EM (ngte=3) source
def kforEMsource():
    return 1.9e-17   # units: s

# TODO: Return the value of kgte for a purely GW (ngte=5) source
def kforGWsource():
    return 1.7e-20   # units: s^3

# TODO: Return the value of kgte for a purely r-mode (ngte=7) source
def kforRmodesource():
    pass

# Returns a 2D list containing the running braking indices and kgte values of a given template in the form [ngte's, kgte's].
def runningbrakingindexandk(template):
    s = len(template[0])
    knotnum = len(template)

    ns = []
    ks = []

    if s == 3:
        for knottemplate in template:
            f0 = knottemplate[0]
            fd0 = knottemplate[1]
            fdd0 = knottemplate[2]

            thisn = fdd0 * f0 / fd0 ** 2
            print(f0)
            print(fd0)
            print(fdd0)
            print(thisn)
            print(f0 ** thisn)
            thisk = -fd0 / f0 ** thisn

            ns.append(thisn)
            ks.append(thisk)

    if s == 2:
        for i, knottemplate in enumerate(template[:-1]):
            f0 = knottemplate[0]
            fd0 = knottemplate[1]
            f1 = template[i + 1][0]
            fd1 = template[i + 1][1]

            thisn = 1 + (f0 * fd1 - f1 * fd0) / (fd0 * fd1 * (bf.knotslist[i + 1] - bf.knotslist[i]))
            thisk = - fd0 / f0 ** thisn

            ns.append(thisn)
            ks.append(thisk)

    if s == 1 and knotnum >= 2:
        for i, knottemplate in enumerate(template[:-2]):
            f0 = knottemplate[0]
            f1 = template[i + 1][0]
            f2 = template[i + 2][0]

            t0 = bf.knotslist[i]
            t1 = bf.knotslist[i + 1]
            t2 = bf.knotslist[i + 2]

            fd0 = (f1 - f0) / (t1 - t0)
            fd1 = (f2 - f1) / (t2 - t1)

            thisn = 1 + (f0 * fd1 - f1 * fd0) / (fd0 * fd1 * (t1 - t0))

            try:
                thisk = - fd0 / f0 ** thisn
            except OverflowError as error:
                print(error)
                print("Setting thisk value to 1")
                thisk = 1

            ns.append(thisn)
            ks.append(thisk)

    return [ns, ks]
"""
#tempstr = "1.43879014E+02 -1.20883749E-05  5.07818350E-12  1.43872629E+02 -1.20856456E-05  5.07611515E-12  1.43864383E+02 -1.20821159E-05  5.07344041E-12  1.43849039E+02 -1.20755063E-05  5.06842882E-12  1.43830442E+02 -1.20672990E-05  5.06218941E-12  1.43805340E+02 -1.20562765E-05  5.05381968E-12  1.43761220E+02 -1.20374432E-05  5.03958354E-12  1.43717172E+02 -1.20182481E-05  5.02505130E-12  1.43659214E+02 -1.19936977E-05  5.00655636E-12  1.43619126E+02 -1.19762575E-05  4.99338843E-12"
#tempstr = "1.43863460E+02  1.54326187E-04 -4.03184188E-06  1.43857764E+02 -1.91187060E-05  2.41020111E-07  1.43851640E+02 -5.85225334E-06 -1.23023207E-07  1.43837184E+02 -1.41526903E-05  2.01676539E-08  1.43820610E+02 -1.13064990E-05  2.90098000E-09  1.43798825E+02 -1.06402984E-05 -3.16688115E-09  1.43759297E+02 -1.08720903E-05  2.67286913E-09  1.43721129E+02 -1.04594391E-05 -1.52489024E-09  1.43669461E+02 -1.03484484E-05  1.20869587E-09  1.43634732E+02 -1.06857562E-05 -3.53909665E-09"
tempstr = "1.24675130E+02 -5.20883350E-06  1.08810580E-12  1.24672379E+02 -5.20822650E-06  1.08787595E-12  1.24668825E+02 -5.20744305E-06  1.08757933E-12  1.24662211E+02 -5.20598580E-06  1.08702774E-12  1.24654192E+02 -5.20421958E-06  1.08635936E-12  1.24643364E+02 -5.20183537E-06  1.08545741E-12  1.24624322E+02 -5.19764456E-06  1.08387284E-12  1.24605294E+02 -5.19345974E-06  1.08229157E-12  1.24580235E+02 -5.18795221E-06  1.08021206E-12  1.24562889E+02 -5.18414247E-06  1.07877464E-12"

tempsplit = tempstr.split()

temp = []
knottemp = []

for string in tempsplit:
        if len(knottemp) < 3:
                knottemp.append(float(string))
        else:
                temp.append(knottemp)
                knottemp = [float(string)]

nsandks = runningbrakingindexandk(temp)
print(nsandks[0])
print()
print(nsandks[1])
"""

# Determines whether rangea fits inside rangeb
def rangeinsiderange(rangea, rangeb):
    if rangea[0] >= rangeb[0] and rangea[1] <= rangeb[1]:
        return True
    else:
        return False

# Quick plotting for the GTE and Taylor Expansions

def texp(tx, f0, ngte, kgte):
    t = tx * 365 * 24 * 3600
    fd0 = gtederivs(t, f0, ngte, kgte, 1)
    fdd0 = gtederivs(t, f0, ngte, kgte, 2)

    return f0 + fd0 * t + 1 / 2 * fdd0 * t ** 2

# A component of the matrix which converts a list of knot template params to taylor expansion params. BUT we force p0 to be zero.
# The reason being is that this matrix is very ill-conditioned/numerically unstable. By substituting p0 = 0 and making p1
# be the difference between them, the condition number is reduced, as well as this matrix being simplified (it contains
# many zeros if p0 = 0). A few things to consider, because we force p0 = 0, it may be important in other places where we
# are translating between PW parameters and doppler parameters that all times have tstart subtracted from them. So long as
# all times have this substraction, the calculated doppler parameters should be the same as if the times were left as they
# are. Be careful out there!
def knotmatrixdash(p0orig, p1orig):

    p0 = 0
    p1 = p1orig - p0orig

    # While in this example we have specifically set p0 = 0 and p1 = p1 - p0, we keep the full form of the below matrix.
    # The below matrix is the more general case for any arbitrary knot p0 and p1. We set p0 = 0 in the above as this
    # matrix is extremely ill-conditioned, and by having this p0 = 0 the matrix has a lowered condition number. Note
    # that if this matrix is used, certain translations on time may be required in order to produce the correct taylor
    # expansion model values.

    kmd = [[-10 * p0 ** 2 * p1 ** 3 + 5 * p0 * p1 ** 4 - p1 ** 5, 4 * p0 ** 2 * p1 ** 3 - p0 * p1 ** 4,
            -(p0 ** 2 * p1 ** 3) / 2., p0 ** 5 - 5 * p0 ** 4 * p1 + 10 * p0 ** 3 * p1 ** 2,
            -(p0 ** 4 * p1) + 4 * p0 ** 3 * p1 ** 2, (p0 ** 3 * p1 ** 2) / 2.],
           [30 * p0 ** 2 * p1 ** 2, -12 * p0 ** 2 * p1 ** 2 - 4 * p0 * p1 ** 3 + p1 ** 4,
            (3 * p0 ** 2 * p1 ** 2) / 2. + p0 * p1 ** 3, -30 * p0 ** 2 * p1 ** 2,
            p0 ** 4 - 4 * p0 ** 3 * p1 - 12 * p0 ** 2 * p1 ** 2, -(p0 ** 3 * p1) - (3 * p0 ** 2 * p1 ** 2) / 2.],
           [-60 * p0 ** 2 * p1 - 60 * p0 * p1 ** 2, 24 * p0 ** 2 * p1 + 36 * p0 * p1 ** 2,
            -3 * p0 ** 2 * p1 - 6 * p0 * p1 ** 2 - p1 ** 3, 60 * p0 ** 2 * p1 + 60 * p0 * p1 ** 2,
            36 * p0 ** 2 * p1 + 24 * p0 * p1 ** 2, p0 ** 3 + 6 * p0 ** 2 * p1 + 3 * p0 * p1 ** 2],
           [60 * p0 ** 2 + 240 * p0 * p1 + 60 * p1 ** 2, -24 * p0 ** 2 - 120 * p0 * p1 - 36 * p1 ** 2,
            3 * p0 ** 2 + 18 * p0 * p1 + 9 * p1 ** 2, -60 * p0 ** 2 - 240 * p0 * p1 - 60 * p1 ** 2,
            -36 * p0 ** 2 - 120 * p0 * p1 - 24 * p1 ** 2, -9 * p0 ** 2 - 18 * p0 * p1 - 3 * p1 ** 2],
           [-360 * p0 - 360 * p1, 168 * p0 + 192 * p1, -24 * p0 - 36 * p1, 360 * p0 + 360 * p1, 192 * p0 + 168 * p1,
            36 * p0 + 24 * p1],
           [720, -360, 60, -720, -360, -60]]

    return kmd

# A component of the matrix which converts a list of knot template params to taylor expansion params
def dpmat(p0, p1):
    dp = p0 - p1

    mat = np.diag([dp ** 5, dp ** 4, dp ** 3, dp ** 5, dp ** 4, dp ** 3])

    return mat

"""
p0 = 9 * 10 ** 8
p1 = p0 + 528
print(np.matmul(knotmatrixdash(p0, p1), np.linalg.inv(dpmat(p0, p1))))
print(np.linalg.cond(knotmatrixdash(p0, p1), np.linalg.inv(dpmat(p0, p1))))
print()
print(np.matmul(knotmatrixdash(p0, p1, settozero=False), np.linalg.inv(dpmat(p0, p1))))
print(np.linalg.cond(np.matmul(knotmatrixdash(9 * 10 ** 8, 9 * 10 ** 8 + 528, settozero=False), np.linalg.inv(dpmat(p0, p1)))))
"""

# Converts a list of template params [f00, f01, f02, f10, f11, f12] to the equivalent Taylor expansion parameters. The
# Taylor expansion parameters are those of the form f0 + f1 t + 1/2 f2 t^2 + ... This method is only valid when S = 3.
# Again, this problem is ill-conditioned, by the same reasoning given for the knotmatrixdash method, subtracting the
# start time of your data (bf.knotslist[0]) will give a better result.
def KnotTempToTExpParams(temp, p0, p1):
    kmdash = knotmatrixdash(p0, p1)
    dp = dpmat(p0, p1)

    """
    kmdashcond = np.linalg.cond(kmdash)
    dpcond = np.linalg.cond(np.linalg.inv(dp))
    print("kmdash and dp condition numbers")
    print([kmdashcond, dpcond, kmdashcond * dpcond])
    """

    pmat = np.matmul(kmdash, np.linalg.inv(dp))

    #texpparams = np.matmul(pmat, temp)

    # Conditioning. Greatest reduction in condition number of matrix A which is not symmetric by scaling of diagonal
    # matrix D by AD. Source: http://ftp.demec.ufpr.br/CFD/bibliografia/Higham_2002_Accuracy%20and%20Stability%20of%20Numerical%20Algorithms.pdf

    pinv = np.linalg.inv(pmat)
    pmattrans = np.transpose(pmat)

    diagelems = []

    for i in range(len(pinv)):
        b = np.linalg.norm(pinv[i])
        a = np.linalg.norm(pmattrans[i])

        elem = (b / a) ** (1/2)
        diagelems.append(elem)

    diag = np.diag(diagelems)
    tempdash = np.matmul(np.linalg.inv(diag), temp)

    pmatdash = np.matmul(pmat, diag)
    """
    condpmat = np.linalg.cond(pmat)
    condpmatdash = np.linalg.cond(pmatdash)
    print("Old and new condition numbers:")
    print([condpmat, condpmatdash, np.log10(condpmat / condpmatdash)])
    """

    newtexpparams = np.matmul(pmatdash, tempdash)

    return newtexpparams

# Calculated the value of a taylor expansion given parameters [f0, f1, f2, f3, f4, f5] at a time t with tref
def TExpModelValue(texpparams, t, tref=0):
    f0 = texpparams[0]
    f1 = texpparams[1]
    f2 = texpparams[2]
    f3 = texpparams[3]
    f4 = texpparams[4]
    f5 = texpparams[5]

    ti = t - tref

    return f0 + f1 * ti + 1 / 2 * f2 * ti ** 2 + 1 / 6 * f3 * ti ** 3 + 1 / 24 * f4 * ti ** 4 + 1 / 120 * f5 * ti ** 5

# Converts taylor expansion parameters f0 + f1 * t + ... + 1/120 f5 t^5 to taylor parameters f0 + f1 (t - tref) + ...
# + 1/120 f5 (t - tref)^5. Once again, this is an ill-conditioned problem and it is worth subtracting the start time
# of your data (bf.knotlist[0]) to reduce the condition number of the problem, by the same reasoning given for the
# knotmatrixdash method
def TExpParamsToTExpRefParams(texpparams, tref):
    refmatinv = [[1, tref, (tref ** 2) / 2., (tref ** 3) / 6., (tref ** 4) / 24., (tref ** 5) / 120.],
                 [0,    1,             tref, (tref ** 2) / 2.,  (tref ** 3) / 6., (tref ** 4) / 24.],
                 [0,    0,                1,             tref,  (tref ** 2) / 2., (tref ** 3) / 6.],
                 [0,    0,                0,                1,              tref, (tref ** 2) / 2.],
                 [0,    0,                0,                0,                 1, tref],
                 [0,    0,                0,                0,                 0, 1]]

    refmat = np.linalg.inv(refmatinv)
    #texprefparams = np.matmul(refmatinv, texpparams)

    # Conditioning. Greatest reduction in condition number of matrix A which is not symmetric by scaling of diagonal
    # matrix D by AD. Source: http://ftp.demec.ufpr.br/CFD/bibliografia/Higham_2002_Accuracy%20and%20Stability%20of%20Numerical%20Algorithms.pdf
    refmatinvtrans = np.transpose(refmatinv)

    diagelems = []

    for i in range(len(refmat)):
        b = np.linalg.norm(refmat[i])
        a = np.linalg.norm(refmatinvtrans[i])

        elem = (np.abs(b) / np.abs(a)) ** (1/2)
        diagelems.append(elem)

    diag = np.diag(diagelems)
    refmatinvdash = np.matmul(refmatinv, diag)

    #trefparams = np.matmul(refmatinvdash, np.matmul(np.linalg.inv(diag), texpparams))
    trefparams = np.matmul(refmatinvdash, np.matmul(np.linalg.inv(diag), texpparams))
    """
    print("Tref old and new cond numbers:")
    refmatinvcond = np.linalg.cond(refmatinv)
    refmatinvdashcond = np.linalg.cond(refmatinvdash)
    print([refmatinvcond, refmatinvdashcond, np.log10(refmatinvcond / refmatinvdashcond)])
    """

    return trefparams

# The matrix which transforms from the PW parameters to tref parameters
def ParamTransformationMatrix(tstart, tend, reftime, s):

    dt   = tend - tstart
    dref = reftime - tstart

    if s == 3:
        matrix = [[1 -(6*dref**5)/dt**5 + (15*dref**4)/dt**4 - (10*dref**3)/dt**3, dref - (3*dref**5)/dt**4 + (8*dref**4)/dt**3 - (6*dref**3)/dt**2, dref**2/2. - dref**5/(2.*dt**3) + (3*dref**4)/(2.*dt**2) - (3*dref**3)/(2.*dt), (6*dref**5)/dt**5 - (15*dref**4)/dt**4 + (10*dref**3)/dt**3, (-3*dref**5)/dt**4 + (7*dref**4)/dt**3 - (4*dref**3)/dt**2,dref**5/(2.*dt**3) - dref**4/dt**2 + dref**3/(2.*dt)],
                  [ (-30*dref**4)/dt**5 + (60*dref**3)/dt**4 - (30*dref**2)/dt**3, 1 - (15*dref**4)/dt**4 + (32*dref**3)/dt**3 - (18*dref**2)/dt**2, dref - (5*dref**4)/(2.*dt**3) + (6*dref**3)/dt**2 - (9*dref**2)/(2.*dt), (30*dref**4)/dt**5 - (60*dref**3)/dt**4 + (30*dref**2)/dt**3, (-15*dref**4)/dt**4 + (28*dref**3)/dt**3 - (12*dref**2)/dt**2, (5*dref**4)/(2.*dt**3) - (4*dref**3)/dt**2 + (3*dref**2)/(2.*dt)],
                  [(-120*dref**3)/dt**5 + (180*dref**2)/dt**4 - (60*dref)/dt**3, (-60*dref**3)/dt**4 + (96*dref**2)/dt**3 - (36*dref)/dt**2, 1 - (10*dref**3)/dt**3 + (18*dref**2)/dt**2 - (9*dref)/dt, (120*dref**3)/dt**5 - (180*dref**2)/dt**4 + (60*dref)/dt**3, (-60*dref**3)/dt**4 + (84*dref**2)/dt**3 - (24*dref)/dt**2, (10*dref**3)/dt**3 - (12*dref**2)/dt**2 + (3*dref)/dt],
                  [(-360*dref**2)/dt**5 + (360*dref)/dt**4 - 60/dt**3, (-180*dref**2)/dt**4 + (192*dref)/dt**3 - 36/dt**2, (-30*dref**2)/dt**3 + (36*dref)/dt**2 - 9/dt, (360*dref**2)/dt**5 - (360*dref)/dt**4 + 60/dt**3, (-180*dref**2)/dt**4 + (168*dref)/dt**3 - 24/dt**2, (30*dref**2)/dt**3 - (24*dref)/dt**2 + 3/dt],
                  [(-720*dref)   /dt**5 + 360/dt**4, (-360*dref)/dt**4 + 192/dt**3, (-60*dref)/dt**3 + 36/dt**2, (720*dref)/dt**5 - 360/dt**4, (-360*dref)/dt**4 + 168/dt**3, (60*dref)/dt**3 - 24/dt**2],
                  [ -720         /dt**5, -360/dt**4, -60/dt**3, 720/dt**5, -360/dt**4, 60/dt**3]]

    elif s == 2:

        matrix = [[1 - (3*dref**2)/dt**2 - (2*dref**3)/dt**3,dref + (2*dref**2)/dt + dref**3/dt**2,(3*dref**2)/dt**2 + (2*dref**3)/dt**3,dref**2/dt + dref**3/dt**2],
                  [(-6*dref)/dt**2 - (6*dref**2)/dt**3,1 + (4*dref)/dt + (3*dref**2)/dt**2,(6*dref)/dt**2 + (6*dref**2)/dt**3,(2*dref)/dt + (3*dref**2)/dt**2],
                  [-6/dt**2 - (12*dref)/dt**3,4/dt + (6*dref)/dt**2,6/dt**2 + (12*dref)/dt**3,2/dt + (6*dref)/dt**2],
                  [-12/dt**3,6/dt**2,12/dt**3,6/dt**2]]

    return matrix

# The matrix used for conditioning the PW params to Tref params matrix. This matrix is included to hopefully avoid re-computing it over and over to
# increase efficiency. The matrix was worked out by computing the commented out section symbolically in Mathematica for the 'ParamTransformationMatrix(tstart, tend, reftime)'
# method.
def ParamTransformationConditioningMatrix(tstart, tend, reftime):
    # Algorithm used for determining conditioning. Source: http://ftp.demec.ufpr.br/CFD/bibliografia/Higham_2002_Accuracy%20and%20Stability%20of%20Numerical%20Algorithms.pdf
    """
    diagelems = []

    matrixinv = np.linalg.inv(matrix)

    matrixinvtrans = np.transpose(matrixinv)

    for i in range(len(matrix)):
        b = np.linalg.norm(matrix[i])
        a = np.linalg.norm(matrixinvtrans[i])

        elem = (np.abs(b) / np.abs(a)) ** (1/2)
        diagelems.append(elem)

    diagmat = np.diag(diagelems)
    """

    dt   = tend - tstart
    dref = reftime - tstart

    elems = [abs(abs(1 - (6*dref**5)/dt**5 + (15*dref**4)/dt**4 - (10*dref**3)/dt**3)**2 + abs((6*dref**5)/dt**5 - (15*dref**4)/dt**4 + (10*dref**3)/dt**3)**2 +
           abs(dref - (3*dref**5)/dt**4 + (8*dref**4)/dt**3 - (6*dref**3)/dt**2)**2 + abs((-3*dref**5)/dt**4 + (7*dref**4)/dt**3 - (4*dref**3)/dt**2)**2 +
           abs(dref**2/2. - dref**5/(2.*dt**3) + (3*dref**4)/(2.*dt**2) - (3*dref**3)/(2.*dt))**2 + abs(dref**5/(2.*dt**3) - dref**4/dt**2 + dref**3/(2.*dt))**2)**
         0.25/2**0.25,
     abs(abs((-30*dref**4)/dt**5 + (60*dref**3)/dt**4 - (30*dref**2)/dt**3)**2 +
           abs((30*dref**4)/dt**5 - (60*dref**3)/dt**4 + (30*dref**2)/dt**3)**2 + abs(1 - (15*dref**4)/dt**4 + (32*dref**3)/dt**3 - (18*dref**2)/dt**2)**2 +
           abs((-15*dref**4)/dt**4 + (28*dref**3)/dt**3 - (12*dref**2)/dt**2)**2 + abs(dref - (5*dref**4)/(2.*dt**3) + (6*dref**3)/dt**2 - (9*dref**2)/(2.*dt))**2 +
           abs((5*dref**4)/(2.*dt**3) - (4*dref**3)/dt**2 + (3*dref**2)/(2.*dt))**2)**0.25/
        abs(2 + abs(dref)**2 + abs(((-8640*dref)/dt**9 + 8640/dt**8)*dt**9)**2/7.46496e7)**0.25,

       abs(abs((-120*dref**3)/dt**5 + (180*dref**2)/dt**4 - (60*dref)/dt**3)**2 + abs((120*dref**3)/dt**5 - (180*dref**2)/dt**4 + (60*dref)/dt**3)**2 +
           abs((-60*dref**3)/dt**4 + (96*dref**2)/dt**3 - (36*dref)/dt**2)**2 + abs((-60*dref**3)/dt**4 + (84*dref**2)/dt**3 - (24*dref)/dt**2)**2 +
           abs(1 - (10*dref**3)/dt**3 + (18*dref**2)/dt**2 - (9*dref)/dt)**2 + abs((10*dref**3)/dt**3 - (12*dref**2)/dt**2 + (3*dref)/dt)**2)**0.25/
        abs(2 + abs(dref)**2 + abs(dref)**4/4. + abs(((-8640*dref)/dt**9 + 8640/dt**8)*dt**9)**2/7.46496e7 +
           abs(((4320*dref**2)/dt**9 - (8640*dref)/dt**8 + 4320/dt**7)*dt**9)**2/7.46496e7)**0.25,

       abs(abs((-360*dref**2)/dt**5 + (360*dref)/dt**4 - 60/dt**3)**2 + abs((360*dref**2)/dt**5 - (360*dref)/dt**4 + 60/dt**3)**2 +
           abs((-180*dref**2)/dt**4 + (192*dref)/dt**3 - 36/dt**2)**2 + abs((-180*dref**2)/dt**4 + (168*dref)/dt**3 - 24/dt**2)**2 +
           abs((-30*dref**2)/dt**3 + (36*dref)/dt**2 - 9/dt)**2 + abs((30*dref**2)/dt**3 - (24*dref)/dt**2 + 3/dt)**2)**0.25/
        abs(abs(dref)**2 + abs(dref)**4/4. + abs(dref)**6/36. + abs(((-8640*dref)/dt**9 + 8640/dt**8)*dt**9)**2/7.46496e7 +
           abs(((4320*dref**2)/dt**9 - (8640*dref)/dt**8 + 4320/dt**7)*dt**9)**2/7.46496e7 +
           abs(((-1440*dref**3)/dt**9 + (4320*dref**2)/dt**8 - (4320*dref)/dt**7 + 1440/dt**6)*dt**9)**2/7.46496e7)**0.25,

       abs(abs((720*dref)/dt**5 - 360/dt**4)**2 + abs((-720*dref)/dt**5 + 360/dt**4)**2 + abs((-360*dref)/dt**4 + 168/dt**3)**2 +
           abs((-360*dref)/dt**4 + 192/dt**3)**2 + abs((60*dref)/dt**3 - 24/dt**2)**2 + abs((-60*dref)/dt**3 + 36/dt**2)**2)**0.25/
        abs(abs(dref)**4/4. + abs(dref)**6/36. + abs(dref)**8/576. + abs(((4320*dref**2)/dt**9 - (8640*dref)/dt**8 + 4320/dt**7)*dt**9)**2/7.46496e7 +
           abs(((-1440*dref**3)/dt**9 + (4320*dref**2)/dt**8 - (4320*dref)/dt**7 + 1440/dt**6)*dt**9)**2/7.46496e7 +
           abs(((360*dref**4)/dt**9 - (1440*dref**3)/dt**8 + (2160*dref**2)/dt**7 - (1440*dref)/dt**6 + 360/dt**5)*dt**9)**2/7.46496e7)**0.25,

       abs(1036800/abs(dt)**10 + 259200/abs(dt)**8 + 7200/abs(dt)**6)**0.25/
        abs(abs(dref)**6/36. + abs(dref)**8/576. + abs(dref)**10/14400. +
           abs(((-1440*dref**3)/dt**9 + (4320*dref**2)/dt**8 - (4320*dref)/dt**7 + 1440/dt**6)*dt**9)**2/7.46496e7 +
           abs(((360*dref**4)/dt**9 - (1440*dref**3)/dt**8 + (2160*dref**2)/dt**7 - (1440*dref)/dt**6 + 360/dt**5)*dt**9)**2/7.46496e7 +
           abs(((-72*dref**5)/dt**9 + (360*dref**4)/dt**8 - (720*dref**3)/dt**7 + (720*dref**2)/dt**6 - (360*dref)/dt**5 + 72/dt**4)*dt**9)**2/7.46496e7)**0.25]

    return np.diag(elems)

previousconds = []
newconds = []

# An alternative method to transoform our PW parameters to Doppler parameters in a taylor expansion with (t - tref) terms.
# Like the above methods, this is ill-conditioned, however less so than the above. Again similar to the above methods,
# this method works best if the start time (bf.knotslist[0]) of your data is subtracted from all time elements, p0, p1
# and tref. If you are using the conditioning in this method, it slows down dramatically. Good luck
def PWParamstoTrefParams(pwparams, p0, p1, tref, s):

    matrix = ParamTransformationMatrix(p0, p1, tref, s)

    print(len(matrix))
    print(len(pwparams))
    print(pwparams)
    # If no conditioning is required
    return np.matmul(matrix, pwparams)

    # Uncomment for conditioning
    """
    diagmat = ParamTransformationConditioningMatrix(p0, p1, tref)

    pwparamsdash = np.matmul(np.linalg.inv(diagmat), pwparams)

    matrixdash = np.matmul(matrix, diagmat)
    trefparams = np.matmul(matrixdash, pwparamsdash)

    return trefparams
    """

# As the above method, but now we accept the transformation and conditioning matrices as input. In this way, if this method is needed a large
# number of times, these matrices can be computed once initially and then reused. Conditioning is not always required, so consider just
# multiplying the transformation matrix with your PW parameters vector directly.
def PWParamstoTrefParamsPreComputed(pwparams, PWtoTrefMat, CondMat):

    pwparamsdash = np.matmul(np.linalg.inv(CondMat), pwparams)

    matrixdash = np.matmul(PWtoTrefMat, CondMat)
    trefparams = np.matmul(matrixdash, pwparamsdash)

    return trefparams

# Plot a PW model. pwparams only needs to be a list of parameters. Knots should be defined before using this method
def PlotPWModel(pwparams, show=True, label="", linewidth=2):
        times = np.linspace(bf.knotslist[0], bf.knotslist[-1], 200)
        values = []

        s = int(len(pwparams) / len(bf.knotslist))

        basiscoeffs = bf.allcoeffs(s)

        for t in times:
                segment = 0

                for i in range(len(bf.knotslist)):
                        if bf.knotslist[i] <= t <= bf.knotslist[i + 1]:
                                segment = i
                                break

                value = 0
                s = int(len(pwparams)/len(bf.knotslist))

                for spindown in range(s):
                        value += pwparams[segment * s + spindown] * bf.basisfunctionvalue(t, segment, 0, spindown, basiscoeffs) + pwparams[(segment + 1) * s + spindown] * bf.basisfunctionvalue(t, segment, 1, spindown, basiscoeffs)

                values.append(value)

        plt.plot(times, values, label=label, linewidth=linewidth)

        if show:
                plt.legend()
                plt.show()

# Not finished
def isvalidtemplate(template, tbank):

        # Check initial frequency is within global range.
        # purposes
        if not (tbank.fmin <= template[0] <= tbank.fmax):
                return False

        kmin = tbank.kmin
        kmax = tbank.kmax

        # Check first spin down on first knot
        if not (-kmax * tbank.fmax ** tbank.nmax <= template[1] <= -kmin * tbank.fmin ** tbank.nmin):
                return False

        s = (len(template) / len(bf.knotslist))

        knottemplates = []
        thisknottemplate = []

        for elem in template:
                thisknottemplate.append(elem)
                if len(thisknottemplate) == s:
                        knottemplates.append(thisknottemplate)
                        thisknottempalte = []

        ns = [knot[2] * knot[0] / knot[1] ** 2 for knot in knottemplates]
        ks = [-knot[1] / knot[0] ** (knot[2] * knot[0] / knot[1] ** 2) for knot in knottemplates]

        ntol = tbank.ntol
        ktol = tbank.ktol

        # Check all braking indices are within global range
        if max(ns) > tbank.nmax or min(ns) < tbank.nmin:
                return False

        # Check all kgte values are within global range
        if max(ks) > kmax or min(ks) < kmin:
                return False

"""
day = 24 * 3600
year = 365 * day
ts = np.linspace(0, 150, 1000)
fs3 = gte(ts, 1000, 3, 10**-15)
fs5 = gte(ts, 1000, 5, 10**-15)
fs7 = gte(ts, 1000, 7, 10**-15)
#fts = texp(ts, 200, 3, 10**-16)
#plt.plot(ts, fs, label="Expected Spin Down")
plt.plot(ts, fs3, label="EM radiation: ngte = 3")
plt.plot(ts, fs5, label="GW radiation: ngte = 5")
plt.plot(ts, fs7, label="R-modes: ngte = 7")
plt.xlabel("Seconds from NS Birth")
plt.ylabel("Neutron Star Spin Frequency")
plt.legend()
plt.show()
"""
