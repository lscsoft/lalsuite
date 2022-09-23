import numpy as np
from scipy import integrate
import BasisFunctions as bf
import SamplingMethods as sm
import GTEandOtherMethods as gom
import MyErrors
import matplotlib.pyplot as plt
import copy
import time


# Our b vector for MOLS
def bvec(points, f0, n, kgte):
    return [gom.gte(t - bf.knotslist[0], f0, n, kgte) for t in points]


# Constructs an individual row of the a matrix for the sample point 'point' and the coefficients of our basis functions
# 'coeffs'
def rowofa(point, coeffs, s):
    j = sm.thisint(point)
    ints = len(bf.knotslist) - 1
    firstzeros = [0] * j * s
    lastzeros = [0] * (ints - 1 - j) * s

    coeffsb = []
    coeffsc = []

    for thiss in range(s):
        coeffb = bf.basisfunctionvalue(point, j, 0, thiss, coeffs)
        coeffc = bf.basisfunctionvalue(point, j, 1, thiss, coeffs)

        coeffsb.append(coeffb)
        coeffsc.append(coeffc)

    return firstzeros + coeffsb + coeffsc + lastzeros


# Builds our a matrix
def a(points, coeffs, s):
    amat = []

    for point in points:
        amat.append(rowofa(point, coeffs, s))

    return amat


# Returns a diagonal matrix containing the diagonal elements of mat to power -1/2
def pmat(mat):
    return np.diag(bf.listpseudoinv(np.diagonal(mat), 1 / 2))


# Conditions mat by multiplying mat by pmat on each of its sides
def matdash(mat):
    pm = pmat(mat)

    return np.matmul(pm, np.matmul(mat, pm))


# Returns mat multiplied by its transpose
def ata(mat):
    return np.matmul(np.transpose(mat), mat)


# Returns the solutions for our LSTSQ problem
def sols(coeffs, ppint, s, f0, n, kgte, conditioning=True):
    while True:
        try:
            points = sm.samplepoints(ppint, f0, n, kgte)
            break
        except MyErrors.SegmentContainsNoSamplePoints:
            print("Reattempting MOLS fitting with greater sampling")
            ppint *= 2

    b = bvec(points, f0, n, kgte)
    amat = a(points, coeffs, s)

    atamat = ata(amat)

    if conditioning:
        pm = pmat(atamat)
    else:
        pm = np.identity(len(atamat))

    lhs = np.matmul(pm, np.matmul(atamat, pm))
    rhs = np.matmul(np.matmul(pm, np.transpose(amat)), b)

    # print("Condition number for matrix A  is: " + "{:.2E}".format(np.linalg.cond(amat)))
    # print("Condition number for matrix A' is: " + "{:.2E}".format(np.linalg.cond(lhs)))

    try:
        params = np.matmul(pm, np.linalg.solve(lhs, rhs))
    except np.linalg.LinAlgError as error:
        print(error)
        print("Error in calculating MOLS parameters, using Python LSTSQ method instead")
        print()
        params = np.matmul(pm, np.linalg.lstsq(lhs, rhs)[0])

    return params


# Partitions the 1D list params into a 3D list for our parameters. Extract parameter by [ int ][ B or C ][ k ]
def solsbyint(params, s):
    ints = int(len(params) / s) - 1

    partedsols = np.zeros((ints, 2, s))

    for i in range(ints):
        solsb = np.array(params[i * s: (i + 1) * s])
        solsc = np.array(params[(i + 1) * s: (i + 2) * s])

        partedsols[i][0] = solsb
        partedsols[i][1] = solsc

    return partedsols


def solsbetweenknots(knotnuma, knotnumb, coeffs, ppint, s, f0, n, kgte, conditioning=True):

    while True:
        try:
            points = sm.samplepointswithinknots(knotnuma, knotnumb, ppint, f0, n, kgte)
            break
        except MyErrors.SegmentContainsNoSamplePoints:
            print("Reattempting MOLS fitting with greater sampling")
            ppint *= 2

    b = bvec(points, f0, n, kgte)
    amat = a(points, coeffs, s)

    for i in range(knotnuma * s):
        [row.pop(0) for row in amat]

    atamat = ata(amat)

    if conditioning:
        pm = pmat(atamat)
    else:
        pm = np.identity(len(atamat))

    lhs = np.matmul(pm, np.matmul(atamat, pm))
    rhs = np.matmul(np.matmul(pm, np.transpose(amat)), b)

    # print("Condition number for matrix A  is: " + "{:.2E}".format(np.linalg.cond(amat)))
    # print("Condition number for matrix A' is: " + "{:.2E}".format(np.linalg.cond(lhs)))

    try:
        params = np.matmul(pm, np.linalg.solve(lhs, rhs))
    except np.linalg.LinAlgError as error:
        print(error)
        print("Error in calculating MOLS parameters, using Python LSTSQ method instead")
        print()
        params = np.matmul(pm, np.linalg.lstsq(lhs, rhs)[0])

    zerobuff = [0] * (s * knotnuma)
    paramszerobuff = zerobuff + list(params)

    return paramszerobuff


# Comparison of the two different 'sols' methods
"""
s = 3
f0 = 1000
n = 5
kgte = 10 ** -14
ppint = 30
bf.knotslist = [0, 20, 40, 60, 80, 100]
coeffs = bf.allcoeffs(s)
knotnuma = 0
knotnumb = 2

print(bf.knotslist)
print(list(solsbyint(list(sols(coeffs, ppint, s, f0, n, kgte, conditioning=True)), s)))
print(list(solsbyint(list(solsbetweenknots(knotnuma, knotnumb, coeffs, ppint, s, f0, n, kgte, conditioning=True)), s)))
"""
points = []

# Calculates model value at given point for basis function coeffs 'coeffs' and parameter values 'params'
def modelvalueatpoint(point, coeffs, params, ignoreintcheck=False, singleseg=False):
    """
    f0 = params[0][0][0]
    f1 = params[0][0][1]
    f2 = params[0][0][2]
    
    t = point - bf.knotslist[0]
    #print("Model t: " + str(t))
    
    return f0 + f1 * t + 1/2 * f2 * t ** 2
    """

    if ignoreintcheck:
        try:
            j = sm.thisint(point)
        except MyErrors.PointNotWithinKnotBoundaries:

            t = point - bf.knotslist[0]

            if t < 0:
                f0 = params[0][0][0]
                f1 = params[0][0][1]
                f2 = params[0][0][2]

                return f0 + f1 * t + 1/2 * f2 * t ** 2
            elif t > bf.knotslist[-1] - bf.knotslist[0]:
                f0 = params[-1][1][0]
                f1 = params[-1][1][1]
                f2 = params[-1][1][2]

                return f0 + f1 * t + 1/2 * f2 * t ** 2
            
            return 0
    else:
        j = sm.thisint(point)

    #if singleseg:
    #    j = 0

    points.append(point)
    
    s = np.shape(params)[2]

    val = 0

    for thiss in range(s):
        basisfuncs = bf.basisfunctionvalue(point, j, 0, thiss, coeffs)
        basisfunce = bf.basisfunctionvalue(point, j, 1, thiss, coeffs)


        if singleseg:
            val += params[0][0][thiss] * basisfuncs
            val += params[0][1][thiss] * basisfunce
        else:
            val += params[j][0][thiss] * basisfuncs
            val += params[j][1][thiss] * basisfunce

    return val


def phase(point, coeffs, params, ignoreintcheck=False):

    model = lambda t: modelvalueatpoint(t, coeffs, params, ignoreintcheck=ignoreintcheck)
    #print("Integral bounds: " + str([bf.knotslist[0], point]))
    phasemodel = 2 * np.pi * integrate.quad(model, bf.knotslist[0], point, epsabs=0)[0]

    return phasemodel
    """
    dt = point - bf.knotslist[0]
    freq = params[0][0][0]
    f1dot = params[0][0][1]
    f2dot = params[0][0][2]
    
    print("Dt is: " + str(dt))
    
    return 2 * np.pi * (freq * dt + f1dot * 0.5 * dt**2 + f2dot * 1/6 * dt**3) #phasemodel
    """


# Builds a list of the 'correct' physical parameters as by the general torque equation
def correctparams(s, f0, n, kgte):
    ints = len(bf.knotslist) - 1
    params = []

    for i in range(ints):
        paramsstart = []
        paramsend = []
        for thisk in range(s):
            ti = bf.knotslist[i] - bf.knotslist[0]
            tip1 = bf.knotslist[i + 1] - bf.knotslist[0]

            paramsstart.append(gom.gtederivs(ti, f0, n, kgte, thisk))
            paramsend.append(gom.gtederivs(tip1, f0, n, kgte, thisk))

        params.append([paramsstart, paramsend])

    return params


# Plots our model
def modelplotter(ppint, s, f0, n, kgte, trueparams=False):
    res = 1

    coeffs = bf.allcoeffs(s)
    params = sols(coeffs, ppint, s, f0, n, kgte)
    partparams = solsbyint(params, s)

    correctps = correctparams(s, f0, n, kgte)

    xpoints = np.linspace(bf.knotslist[0], bf.knotslist[-1], res * (bf.knotslist[-1] - bf.knotslist[0]))

    gtepoints = []
    correctpoints = []
    ypoints = []

    for x in xpoints:
        ypoints.append(modelvalueatpoint(x, coeffs, partparams))
        gtepoints.append(gom.gte(x - bf.knotslist[0], f0, n, kgte))

        if trueparams:
            correctpoints.append(modelvalueatpoint(x, coeffs, correctps))

    plt.plot(xpoints, ypoints, label="LSTSQ Model")
    plt.plot(xpoints, gtepoints, label="GTE")

    if trueparams:
        plt.plot(xpoints, correctpoints, label="Correct Coefficients")

    plt.xlabel("Time (s)")
    plt.ylabel("Rotational Frequency (Hz)")
    plt.ylim(0, f0 * 1.1)
    plt.legend()
    plt.show()


def modelplotterknotspecific(knotnuma, knotnumb, ppint, s, f0, n, kgte, trueparams=False):
    res = 1

    coeffs = bf.allcoeffs(s)
    params = solsbetweenknots(knotnuma, knotnumb, coeffs, ppint, s, f0, n, kgte)
    partparams = solsbyint(params, s)

    correctps = correctparams(s, f0, n, kgte)

    xpoints = np.linspace(bf.knotslist[knotnuma], bf.knotslist[knotnumb], res * (bf.knotslist[knotnumb] - bf.knotslist[knotnuma]))

    gtepoints = []
    correctpoints = []
    ypoints = []

    for x in xpoints:
        ypoints.append(modelvalueatpoint(x, coeffs, partparams))
        gtepoints.append(gom.gte(x - bf.knotslist[0], f0, n, kgte))
        if trueparams:
            correctpoints.append(modelvalueatpoint(x, coeffs, correctps))

    plt.plot(xpoints, ypoints, label="LSTSQ Model")
    plt.plot(xpoints, gtepoints, label="GTE")
    if trueparams:
        plt.plot(xpoints, correctpoints, label="Correct Coefficients")
    plt.xlabel("Time (s)")
    plt.ylabel("Rotational Frequency (Hz)")
    plt.legend()
    plt.show()


# Gives the value of the error bound for a particular point 't' for the interval it lies in
def errorbounds(t):
    j = sm.thisint(t)

    return 1 / (bf.knotslist[j + 1] - bf.knotslist[j])


# Plots our error
def errorplotter(ppint, s, f0, n, kgte):
    res = 1

    coeffs = bf.allcoeffs(s)
    params = sols(coeffs, ppint, s, f0, n, kgte)
    partparams = solsbyint(params, s)

    correctps = correctparams(s, f0, n, kgte)

    xpoints = []

    # Sometimes when plotting the error it is useful to have a higher resolution for plotting at earlier times as we
    # plot them on a log log plot. This for loop chooses plotting points such that there is the same number of points
    # per piecewise interval. This is not as important when just plotting the piecewise model.
    ints = len(bf.knotslist) - 1
    for i in range(ints):
        sublist = np.linspace(bf.knotslist[i], bf.knotslist[i + 1], res * (bf.knotslist[i + 1] - bf.knotslist[i]) / ints)
        for elem in sublist:
            xpoints.append(elem)

    errorpoints = []
    correctpoints = []
    ypoints = []

    for i, x in enumerate(xpoints):
        ypoints.append(np.abs(modelvalueatpoint(x, coeffs, partparams) - gom.gte(x - bf.knotslist[0], f0, n, kgte)))
        correctpoints.append(modelvalueatpoint(x, coeffs, correctps) - gom.gte(x - bf.knotslist[0], f0, n, kgte))
        errorpoints.append(errorbounds(x))

    title = "f0, n, kgte: " + str([f0, n, kgte]) + ". S = " + str(s) + ", Ints = " + str(ints)

    plt.loglog(xpoints, ypoints, label="LSTSQ Model")
    plt.loglog(xpoints, correctpoints, label="Correct Coefficients")
    plt.loglog(xpoints, errorpoints, label="Error Bounds")
    plt.legend()
    plt.title(title)
    plt.show()


def errorplotterknotspecific(knotnuma, knotnumb, ppint, s, f0, n, kgte):
    res = 1

    ints = knotnumb - knotnuma

    coeffs = bf.allcoeffs(s)
    params = solsbetweenknots(knotnuma, knotnumb, coeffs, ppint, s, f0, n, kgte)
    partparams = solsbyint(params, s)

    correctps = correctparams(s, f0, n, kgte)

    xpoints = []

    # Sometimes when plotting the error it is useful to have a higher resolution for plotting at earlier times as we
    # plot them on a log log plot. This for loop chooses plotting points such that there is the same number of points
    # per piecewise interval. This is not as important when just plotting the piecewise model.
    for i in range(knotnuma, knotnumb):
        sublist = np.linspace(bf.knotslist[i], bf.knotslist[i + 1], res * (bf.knotslist[i + 1] - bf.knotslist[i]) / ints)
        for elem in sublist:
            xpoints.append(elem)

    errorpoints = []
    correctpoints = []
    ypoints = []

    for i, x in enumerate(xpoints):
        ypoints.append(np.abs(modelvalueatpoint(x, coeffs, partparams) - gom.gte(x - bf.knotslist[0], f0, n, kgte)))
        correctpoints.append(modelvalueatpoint(x, coeffs, correctps) - gom.gte(x - bf.knotslist[0], f0, n, kgte))
        errorpoints.append(errorbounds(x))

    title = "f0, n, kgte: " + str([f0, n, kgte]) + ". S = " + str(s) + ", Ints = " + str(ints)

    plt.loglog(xpoints, ypoints, label="LSTSQ Model")
    plt.loglog(xpoints, correctpoints, label="Correct Coefficients")
    plt.loglog(xpoints, errorpoints, label="Error Bounds")
    plt.legend()
    plt.title(title)
    plt.show()


# Returns the value of the error of our model at a given point
def errorvalueatpoint(point, coeffs, params, f0, n, kgte):
    error = np.abs(modelvalueatpoint(point, coeffs, params) - gom.gte(point - bf.knotslist[0], f0, n, kgte))

    return error


# Plots a template against the gte with the parameters f0, n and k. Template should be a list of knot templates
def plotatemplate(template, f0, n, kgte):
    s = len(template[0])

    basiscoeffs = bf.allcoeffs(s)

    tempparamsbyint = []

    for i in range(len(template) - 1):
        tempparamsbyint.append([template[i], template[i + 1]])

    res = 600

    xpoints = np.linspace(bf.knotslist[0], bf.knotslist[-1], res)

    ypoints = []
    gtepoints = []

    for x in xpoints:
        modelpoint = modelvalueatpoint(x, basiscoeffs, tempparamsbyint)

        if modelpoint > 2 * 10 ** 3:
            modelpoint = 2 * 10 ** 3
        elif modelpoint < 0:
            modelpoint = 0

        ypoints.append(modelpoint)
        gtepoints.append(gom.gte(x - bf.knotslist[0], f0, n, kgte))

    plt.plot(xpoints, ypoints, label="Template Model")
    plt.plot(xpoints, gtepoints, label="GTE, n = " + str(n))
    plt.legend()
    plt.show()


def plotGTE(f0, n, kgte, ts, te, show=True, label=''):

    xpoints = np.linspace(ts, te, 50)

    ypoints = []

    for x in xpoints:
        ypoints.append(gom.gte(x - bf.knotslist[0], f0, n, kgte))

    plt.plot(xpoints, ypoints, label=label)
    if show:
        plt.show()


"""
tstart = 5000000000000
tdata = 200
tref = tstart + 1/2 * tdata
s = 3
f0 = 200
n = 5
tau = 21600
kgte = gom.kwhichresultsingivenhalflife(tau, f0, n)

bf.knotslist = np.array([tstart, tstart + 1/2 * tdata, tstart + tdata], dtype='float')
someparams = correctparams(s, f0, n, kgte)

ps = 1
pe = ps + 1
texpparams = gom.KnotTempToTExpParams(someparams[ps][0] + someparams[ps][1], bf.knotslist[ps], bf.knotslist[pe])
texprefparams = gom.TExpParamsToTExpRefParams(texpparams, tref - bf.knotslist[ps])
print()
print(bf.knotslist)
print("Basis function template")
print(someparams[ps][0] + someparams[ps][1])
print("Standard Taylor Expansion Template")
print(texpparams)
print("Tref Taylor Expansion template")
print(texprefparams)
print()
print("GTE at tref: " + str(gom.gte(tref - bf.knotslist[ps], f0, n, kgte)) + ", " + str(gom.gte(tref - bf.knotslist[ps], f0, n, kgte) - texprefparams[0]))

times = np.linspace(bf.knotslist[ps], bf.knotslist[pe], 50)
ytemp = []
ytexp = []
ytref = []
diffs = []
diffsref = []
gtevals = []
gtediffstexp = []
gtediffstref = []

for time in times:
    coeffs = bf.allcoeffs(s)

    tempval = modelvalueatpoint(time, coeffs, someparams)
    texpval = gom.TExpModelValue(texpparams, time - bf.knotslist[ps], tref=0)
    texprefval = gom.TExpModelValue(texprefparams, time, tref=tref)
    gteval = gom.gte(time - bf.knotslist[ps], f0, n, kgte)

    ytemp.append(tempval)
    ytexp.append(texpval)
    ytref.append(texprefval)
    gtevals.append(gteval)
    diffs.append(tempval - texpval)
    diffsref.append(tempval - texprefval)

    gtediffstexp.append(gteval - texpval)
    gtediffstref.append(gteval - texprefval)

plt.plot(times, ytemp, label="Standard Template")
plt.plot(times, ytexp, label="Taylor Template")
plt.plot(times, ytref, label="Taylor Ref Template")
#plt.plot(times, gtevals, label="GTE")
plt.legend()
plt.show()

plt.plot(times, diffs, label="Diffs Between Taylor and tref params")
plt.plot(times, gtediffstexp, label="Taylor Diffs with GTE")
#plt.plot(times, diffsref, label="Taylor Ref Diffs")
#plt.plot(times, gtediffstref, label="Taylor Ref Diffs with GTE")
plt.legend()
plt.show()
"""

"""
s = 3
bf.knotslist = np.array([0, 2000, 4000])
#bf.knotslist = np.linspace(0, 4000, 5)
knotnuma = 1
knotnumb = 4
ppint = 30
f0 = 1000
n = 5
kgte = 10 ** -14

print(bf.knotslist)
coeffs = bf.allcoeffs(s)
modelplotter(ppint, s, f0, n, kgte)
#modelplotterknotspecific(1, 4, ppint, s, f0, n, kgte)
#errorplotter(ppint, s, f0, n, kgte)
#errorplotterknotspecific(1, 4, ppint, s, f0, n, kgte)
"""

"""
times = np.linspace(0, 20 * np.pi, 200)
coss = np.cos(times)

noise = [np.random.randint(-100, 100) * 1/100 for i in range(len(times))]
noisedata = coss + noise

plt.plot(times, coss, label="GW Signal")
plt.plot(times, noisedata, label="Data")
plt.legend()
plt.show()
"""
