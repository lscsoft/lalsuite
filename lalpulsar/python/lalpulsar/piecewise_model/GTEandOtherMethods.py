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

from . import BasisFunctions as bf
import lalpulsar as lp

import numpy as np
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
    return (2 ** (ngte - 1) - 1) / ((ngte - 1) * f0 ** (ngte - 1) * kgte)

# Returns value of kgte for a purely EM (ngte=3) source. Typical values: Izz \in [1, 3] * 10^38, mag_field (in Gauss) \in [, 10^15].
# Izz default value is the same as that from the long transient search paper. Defaults give kgte=4.47762e-18
def kforEMsource(Izz=4.34e38, mag_field=1e12, radius=1e4):
    
    c = 299792458
    mu0 = 1.2566e-6
    
    # Magnetic moment of the NS, given as magnetic field (in Tesla) * radius ** 3
    mp = (mag_field / 10 ** 4) * radius ** 3
    
    numerator   = 2 * mp ** 2 * np.pi ** 2
    denominator = 3 * c ** 3 * Izz * mu0
    
    return numerator / denominator   # units: s

# Returns value of kgte for a purely GW (ngte=5) source. Typical values: Izz \in [1, 3] * 10^38, ellip \in [1e-4, 1e-8].
# Izz default value is the same as that from the long transient search paper. Defaults give kgte=1.7182e-20
def kforGWsource(Izz=4.34e38, ellip=1e-4, radius=1e4):
    c = 299792458
    G = 6.6743e-11
    
    numerator = 32 * G * Izz * np.pi ** 4 * ellip ** 2
    denominator = 5 * c ** 5
    
    return numerator / denominator   # units: s^3

# TODO: Return the value of kgte for a purely r-mode (ngte=7) source
def kforRmodesource():
    pass # units: s ^ 5

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
    # Conditioning not currently implemented, not necessary at current
    
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

# Checks to see whether a given template is within the parameter space defined by tbank
def template_inside_p_space(template, tbank):
	
	f0s = []
	f1s = []
	f2s = []
	
	for i, param in enumerate(template):
		
		if i % tbank.s == 0:
			f0s.append(param)
		elif i % tbank.s == 1:
			f1s.append(param)
		elif i % tbank.s == 2:
			f2s.append(param)
	
	# Frequency checks:
	
	if not (tbank.fmin <= f0s[0] <= tbank.fmax):
		return False
	
	for i, f0 in enumerate(f0s[:-1]):
		dt = bf.knotslist[i + 1] - bf.knotslist[i]
		
		this_fmax = gte(dt, f0, tbank.nmin, tbank.kmin)
		this_fmin = gte(dt, f0, tbank.nmax, tbank.kmax)
		
		if not (this_fmin <= f0s[i + 1] <= this_fmax):
			return False
	
	# F1 checks:
	
	for i, f1 in enumerate(f1s):
		
		f0 = f0s[i]
		
		this_f1_max = gtederivs(0, f0, tbank.nmin, tbank.kmin, 1)
		this_f1_min = gtederivs(0, f0, tbank.nmax, tbank.kmax, 1)
		
		if not (this_f1_min <= f1 <= this_f1_max):
			return False
	
	if tbank.s == 3:
		for i, f2 in enumerate(f2s):
			f0 = f0s[i]
			
			this_f2_max = gtederivs(0, f0, tbank.nmax, tbank.kmax, 2)
			this_f2_min = gtederivs(0, f0, tbank.nmin, tbank.kmin, 2)
			
			if not (this_f2_min <= f2 <= this_f1_max):
				return False
	
	return True
	
		

	
	

