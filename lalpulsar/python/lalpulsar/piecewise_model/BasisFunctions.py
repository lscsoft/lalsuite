import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from . import MyErrors

# In this notebook we build the methods required for creating the basis functions of our piecewise model. It is worth
# noting here that initially conditioning of the matrices used here was undertaken. Initially this was necessary as
# when plotting the basis functions all that was returned was numerical fuzz, however many months after having written
# these methods, I cannot recover this behaviour and the conditioning used here seems to increase the condition number
# of the matrices built here. For now the conditioning methods in this notebook have been commented out until a later
# time, however for now without conditioning the results produced appear to be returned accurately despite the
# ill-conditionedness of the problem.

# Lists that the condition numbers of the many matrices required to calculate the condition numbers. Condition numbers
# are added to these lists as the matrices are built. These lists can be used for comparison of the conditioned and
# unconditioned matrices and whether the problems are ill-conditioned or not.
wcondnumbers = []
wdashcondnumbers = []

knotslist = [0.]

# Returns the value of the ith knot. If we are using the methods in the EstimatingKnots notebook, this method should
# simply extract the ith element from the knotslist variable above.
def p(i):
    return knotslist[i]

    if i < 0 or i >= len(knotslist):
        print("Invalid knot number. Length of knotslist is: " + str(len(knotslist)) + ". Knot number requested: " + str(i))
        raise MyErrors.InvalidKnotNumber

    return knotslist[i]

# Our basis function
def u(t, i):
    return (t - p(i)) / (p(i + 1) - p(i))

# Builds the W matrix whose inverse is the coefficients of our basis function
def w(i, s):
    t = Symbol('t')

    ps = p(i)
    pe = p(i + 1)

    matrixupper = []
    matrixlower = []

    thisrowderivs = [u(t, i) ** m for m in range(2 * s)]

    for row in range(s):

        if row != 0:
            for j, elem in enumerate(thisrowderivs):
                thisrowderivs[j] = diff(elem, t, 1)

        thisrows = []
        thisrowe = []

        for elem in thisrowderivs:
            thiselem = lambdify(t, elem, 'numpy')

            thisrows.append(thiselem(ps))
            thisrowe.append(thiselem(pe))

        matrixupper.append(thisrows)
        matrixlower.append(thisrowe)

    return np.array(matrixupper + matrixlower)

# Returns the pseudo inverse of a list to the given power
def listpseudoinv(lst, pwr=1):
    invlist = []

    for elem in lst:
        if elem != 0:
            invlist.append(elem ** -pwr)
        else:
            invlist.append(0)

    return invlist

# Builds a diagonal matrix with elements equal to the inverse of the diagonal elements of the given matrix
def d(mat):
    return np.diag(listpseudoinv(np.diagonal(mat)))

# Returns a conditioned matrix where each row is divided by its diagonal element
def dcond(mat):

    return np.matmul(d(mat), mat)

# Returns the coefficients of our basis function in a 3D list. Reference elements by [ B or C ][ k ][ specific coeff ]
def basiscoeffs(i, s, conditioning=True):
    wmat = w(i, s)
    #wcondnumbers.append(np.linalg.cond(wmat))

    if conditioning:
        dmat = d(wmat)
        dwmat = np.matmul(dmat, wmat)

    else:
        dmat = np.identity(len(wmat))
        dwmat = wmat

    #wmatcondnum = np.linalg.cond(wmat)
    #dwmatcondnum = np.linalg.cond(np.matmul(d(wmat), wmat))

    #wdashcondnumbers.append(np.linalg.cond(dwmat))

    try:
        coeffs = np.transpose(np.linalg.solve(dwmat, dmat))
    except np.linalg.LinAlgError as error:
        print(error)
        print("Error in calculating basis functions. Using Python LSTSQ method instead")
        coeffs = np.linalg.lstsq(dwmat, dmat)[0]

    blist = coeffs[0:s]
    clist = coeffs[s: 2 * s]

    return np.array([blist, clist])

# Returns the coefficients of all basis functions in a 4D list. Reference elements by [ int ][ B or C][ k ][ specific
# coeff ]
def allcoeffs(s):
    coeffs = [basiscoeffs(i, s) for i in range(len(knotslist) - 1)]

    #print("The largest condition number for a W j is: " + "{:.2E}".format(max(wcondnumbers)))
    #print("The largest condition number for a W'j is: " + "{:.2E}".format(max(wdashcondnumbers)))
    return np.array(coeffs)

#allcoeffs(3, 10, 10000)

# Returns the value of a specified basis function given the 4D list of coefficients coeffs.
def basisfunctionvalue(t, i, borc, s, coeffs):
    val = 0

    if t < p(i) or p(i + 1) < t:
        return 0

    thesecoeffs = coeffs[i][borc][s]

    for m, coeff in enumerate(thesecoeffs):
        val += coeff * u(t, i) ** m

    return val

# Plots all basis functions
def allfunctionplotter(s):
    coeffs = allcoeffs(s)
    ints = len(knotslist) - 1

    xpointsintervals = []

    res = 30

    for i in range(ints):
        xpoints = np.linspace(p(i), p(i + 1), res)
        xpointsintervals.append(xpoints)

    ypointsintervals = []

    for i, xpoints in enumerate(xpointsintervals):
        ypointsb = []
        ypointsc = []

        for thisk in range(s):
            ypointsbk = []
            ypointsck = []
            for t in xpoints:
                ypointsbk.append(basisfunctionvalue(t, i, 0, thisk, coeffs))
                ypointsck.append(basisfunctionvalue(t, i, 1, thisk, coeffs))

            ypointsb.append(ypointsbk)
            ypointsc.append(ypointsck)

        ypointsintervals.append(np.array([ypointsb, ypointsc]))

    fig, axs = plt.subplots(ints, s)
    for i in range(ints):
        xpoints = xpointsintervals[i]
        ypointsbck = ypointsintervals[i]
        for thisk in range(s):
            axs[i, thisk].plot(xpoints, ypointsbck[0][thisk])
            axs[i, thisk].plot(xpoints, ypointsbck[1][thisk])
            axs[i, thisk].set_title('B: ' + str(i) + ', ' + str(thisk))

    plt.tight_layout()
    plt.show()

# Plots the basis functions with the given coefficients. Coeffs parameter in the same form as the output of the
# allcoeffs method
def plotcoeffs(coeffs):

    ints = len(coeffs)
    s = len(coeffs[0][0])

    xpointsintervals = []

    res = 30

    for i in range(ints):
        xpoints = np.linspace(p(i), p(i + 1), res)
        xpointsintervals.append(xpoints)

    ypointsintervals = []

    for i, xpoints in enumerate(xpointsintervals):
        ypointsb = []
        ypointsc = []

        for thisk in range(s):
            ypointsbk = []
            ypointsck = []
            for t in xpoints:
                ypointsbk.append(basisfunctionvalue(t, i, 0, thisk, coeffs))
                ypointsck.append(basisfunctionvalue(t, i, 1, thisk, coeffs))

            ypointsb.append(ypointsbk)
            ypointsc.append(ypointsck)

        ypointsintervals.append(np.array([ypointsb, ypointsc]))

    fig, axs = plt.subplots(ints, s)
    for i in range(ints):
        xpoints = xpointsintervals[i]
        ypointsbck = ypointsintervals[i]
        for thisk in range(s):
            axs[i, thisk].plot(xpoints, ypointsbck[0][thisk])
            axs[i, thisk].plot(xpoints, ypointsbck[1][thisk])
            axs[i, thisk].set_title('B: ' + str(i) + ', ' + str(thisk))

    plt.tight_layout()
    plt.show()

"""
s = 3
ints = 2
knotslist = np.array([0, 43200, 86400]) + np.array([20000, 20000, 20000])

coeffseg = allcoeffs(s)
print(knotslist)
print(coeffseg)
plotcoeffs(coeffseg)
"""
