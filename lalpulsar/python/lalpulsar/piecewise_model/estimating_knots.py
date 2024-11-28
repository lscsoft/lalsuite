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

## \file
## \ingroup lalpulsar_python_piecewise_model
"""
Construct a list of knots which gives the maximum allowed spacing between
them while maintaining the desired accuracy.
"""

import ast
import logging
import os.path

import numpy as np

from . import basis_functions as bf
from . import mols_for_gte as mols
from . import errors
from . import semicoherent_metric_methods as scmm

knotarchivefile = None


def setknotarchivepath(path):
    global knotarchivefile
    knotarchivefile = os.path.join(path, "KnotArchive")


# In this module, we use the function, max(|f_GTE(t) - F_PP(t)|) - Delta f_i0,
# to choose the spacing of our knots. The Delta f_i0 is the spacing between the frequency parameter given by the
# metric, it is equal to sqrt(mu/g_ii). Ideally, we wish for the maximum difference between our piecewise model and
# the GTE to be minimised and smaller than the spacing of our frequency parameters. When finding the maximum
# difference between the GTE and our model, max(|f_GTE(t) - F_PP(t)|), we maximise over the value t which falls into
# some interval [p_i, p_i+1]. If we make interval sizes small the difference will go to zero, however this would
# dramatically increase the number of intervals we need, and hence most likely increase the template bank size.
# However, if we make them too large the error becomes greater than our parameter spacing. As such, we wish to
# find the values [p_i, p_i+1] such that the given function above is zero. This is the longest we can make a segment
# before the error becomes greater than our parameter spacing. We start by assuming p_0 = 0 and then inductively
# calculate all of the following knots.

# Working in this module it is important that the knots function, p(i, ints, dur) in the basis_functions file is
# written correctly. See that file for how it should be defined when using this file. As well as this, this module
# requires methods written in the mols_for_gte, basis_functions and semicoherent_metric_methods files, all of which rely
# on the p(i, ints, dur) method. This can be a problem, as for the methods in this module to work, p(i, ints,
# dur) should simply extract an element from a 'knots list' that we build inductively in this module. Hence,
# using the knots generated in this module can at times be tricky. The recommendation for using knots built by this
# file is to either run this file for given model parameter and then use the output as the knots by copying and
# pasting the output into the basis_functions file to be used as all knots, or in the module you wish to run, run
# the GlobalVariableDeclarations module which initialises our knot list (by importing this module and using the
# methods below) and then import the GlobalVariableDeclarations module into the module you wish to run. As an
# example, if you want to run module A which requires knots, import the GlobalVariableDeclarations module into A,
# making sure that in the GlobalVariableDecalartions module the knots list is initialised with one of the methods
# here.
#
# Knots generated in this module are saved to the KnotArchive text file.

# Another word of caution, some methods in this module are outdated, being replaced by more efficient methods later
# in the module. These methods typically run much slower than those which they are replaced by. Further more, many of
# these methods do not make use of the KnotArchive file which stores all generated values of knots. I have tried to
# note these outdated methods were appropriate and add an explanation to why they have been rewritten further in the
# module. They are kept in this module however for completeness and a reference in case they are needed in the
# future

# When using this module be sure to check that the method p and the knotslist parameter in the basis_functions module
# are defined correctly or are being properly updated.

# List that we append the maximum parameter spacing for each segment we build knots for in this module. Kept only for
# inspection, not used in any methods except for printing.
parameterspacings = []

# Returns the parameter spacing associated with the zeroth derivative at the lth knot. Note, if the knotslist variable
# in the basis_functions is never updated, then this method will only work for l = 0
def metricelementspacingonl(l, coeffs, mu):

    metvalue = scmm.metricelementonl(l, l, l, 1, 1, 0, 0, coeffs)

    stepsize = np.sqrt(mu / metvalue)

    return stepsize


# This method is an outdated method, replaced later in this module. Chronologiclly this was replaced first by the
# knotatleff method and then finally by the knotatleffusinggivenknotlist method. No other methods in this module use
# this method.
#
# Calcualtes the time (and hence knot) at which the function, max(|f_GTE(t) - F_PP(t)| - Delta f_i0 is zero. This is
# achieved by sampling all times t between durmin and durmax and seeing if the appropriate value of that function is
# within a given tolerance of zero. If not, we reuse this method with an extended duration of durmax to 2 * durmax. This
# has the obvious flaw that if the correct knot value is below durmin, then this method will not be able to find it.
# Because this method does not make use of faster algorithms, such as the bisection method which we use later on, it
# runs very slowly.
def knotatl(l, s, durmin, durmax, steps, f0, ngte, kgte, mu):
    errorlists = []
    paramspacinglist = []

    durlist = np.linspace(durmin, durmax, steps)

    for dur in durlist:
        basiscoeffs = bf.allcoeffs(s)
        modelparams = mols.solsbyint(mols.sols(basiscoeffs, 20, s, f0, ngte, kgte), s)
        paramspacinglist.append(metricelementspacingonl(l, basiscoeffs, mu))

        suberrorlist = []

        for t in np.linspace(durmin, dur, steps):
            suberrorlist.append(
                mols.errorvalueatpoint(t, basiscoeffs, modelparams, f0, ngte, kgte)
            )
        errorlists.append(suberrorlist)

    for i, dur in enumerate(durlist):

        if np.max(errorlists[i]) - paramspacinglist[i] > 0:
            return dur

    logging.debug("No root found for specified duration, increasing duration span")
    return knotatl(l, s, durmax, 2 * durmax, steps, f0, ngte, kgte)


# This method is not used in our finalised knot generating method. It is however used in the buildcompleteknotlist
# method. It was kept because the buildcompelteknotlist method (and this method) do not read or write to any files.
# As such, if we want to do some experimentation with how we build our knots but not interfere and risk altering our
# library of knots we can use these methods instead.
#
#
# The same as the knotatl method but now done more efficiently by using a bisection method. Instead of sampling all
# possible times/knots, we instead sample two possible durations, negdur and posdur. These knots should be chosen
# such that negdur results in the function, max(|f_GTE(t) - F_PP(t)| - Delta f_i0, being negative and posdur such
# that it is positive. We then select a time inbetween these two possible choices and calculated that function value.
# If the value is within a given tolerance of zero we return that time as the new knot value. If not, we recursively use
# this method by replacing negdur or posdur with this 'inbetween time' in accordance with the bisection method until we
# find a time that is close enough to zero for us to be confident we can use that as the new knot value. This method
# has a maximum recursion depth of 10.
#
# This method is somewhat sensitive to the initial choices of negdur and posdur. The checkfirst parameter is used
# such that negdur and posdur have appropriate values. We must check that negdur and posdur result in negative and
# positive values of the diff function respectively. In certain cases the ill-conditionedness of the basis functions
# we calculate also influences our results. We check for this too in certain cases.
#
# The checkfirst value should always be used when knotatleff is first called and set to True. Once this check has
# initially been cleared we do not need to check the negdur and posdur values again.
#
# In this method we refer to a boolean value fullMOLS, it is defined below.
def knotatleff(
    l,
    s,
    negdur,
    posdur,
    steps,
    f0,
    ngte,
    kgte,
    mu,
    checkfirst=False,
    recursiondepth=0,
    ps=None,
    negcounter=0,
):
    # logging.debug("Current recursion depth: " + str(counter))

    knotnuma = len(bf.knotslist) - 1
    knotnumb = knotnuma + 1

    # As negdur is a candidate knot, it should be greater than the last knot we currently have in bf.knotslist.
    if bf.knotslist[l] >= negdur:
        logging.debug(
            "negdur parameter not within appropriate range, adjusting with larger value"
        )
        logging.debug(
            "negdur = " + str(negdur) + ", previous knot = " + str(bf.knotslist[l])
        )
        return knotatleff(
            l, s, bf.knotslist[l] + 1, posdur, steps, f0, ngte, kgte, mu, True, ps=ps
        )

    if checkfirst:
        bf.knotslist.append(negdur)
        basiscoeffsneg = bf.allcoeffs(s)

        # Can't remember why but this sometimes throws a TypeError. I suspect it probably happens when negdur is too
        # large and the sample points in the sampling_methods module are chosen poorly. Corrected for by using a
        # smaller negdur value. Don't ever remember seeing this error for the posdur parameter though
        try:
            # logging.debug("Prev knot: " + str(ps))
            if fullMOLS:
                modelparamsneg = mols.solsbyint(
                    mols.sols(basiscoeffsneg, 20, s, f0, ngte, kgte), s
                )
            else:
                modelparamsneg = mols.solsbyint(
                    mols.solsbetweenknots(
                        knotnuma, knotnumb, basiscoeffsneg, 20, s, f0, ngte, kgte
                    ),
                    s,
                )
        except TypeError as error:
            logging.debug(error)
            logging.debug("Reducing negdur duration")
            diff = negdur - ps
            bf.knotslist.pop()
            return knotatleff(
                l, s, ps + diff / 10, posdur, steps, f0, ngte, kgte, mu, True, 0, ps=ps
            )

        # Calculating the diff function value for negdur and posdur
        paramspacingneg = metricelementspacingonl(l, basiscoeffsneg, mu)

        errorlistneg = []

        for t in np.linspace(bf.knotslist[-1], negdur, steps):
            errorlistneg.append(
                mols.errorvalueatpoint(
                    t, basiscoeffsneg, modelparamsneg, f0, ngte, kgte
                )
            )

        bf.knotslist.pop()

        bf.knotslist.append(posdur)
        basiscoeffspos = bf.allcoeffs(s)
        if fullMOLS:
            modelparamspos = mols.solsbyint(
                mols.sols(basiscoeffspos, 20, s, f0, ngte, kgte), s
            )
        else:
            modelparamspos = mols.solsbyint(
                mols.solsbetweenknots(
                    knotnuma, knotnumb, basiscoeffspos, 20, s, f0, ngte, kgte
                ),
                s,
            )

        paramspacingpos = metricelementspacingonl(l, basiscoeffspos, mu)

        errorlistpos = []

        for t in np.linspace(bf.knotslist[-1], posdur, steps):
            errorlistpos.append(
                mols.errorvalueatpoint(
                    t, basiscoeffspos, modelparamspos, f0, ngte, kgte
                )
            )

        bf.knotslist.pop()

        maxerrorneg = np.max(errorlistneg)
        maxerrorpos = np.max(errorlistpos)

        logging.debug("Neg and pos durs: " + str([negdur, posdur]))
        logging.debug("Max errors: " + str([maxerrorneg, maxerrorpos]))
        logging.debug("Metric values: " + str([paramspacingneg, paramspacingpos]))
        logging.debug()

        diffvalneg = maxerrorneg - paramspacingneg
        diffvalpos = maxerrorpos - paramspacingpos

        # Checking if posdur and negdur have the appropriate signs for the bisection method and readjusting if not
        if not (diffvalneg < 0 and diffvalpos > 0):
            # logging.debug("Diff vals are: " + str([diffvalneg, diffvalpos]))
            if diffvalpos < 0:
                # logging.debug("Negative pos val")
                diff = posdur - ps
                return knotatleff(
                    l, s, negdur, ps + 2 * diff, steps, f0, ngte, kgte, mu, True, ps=ps
                )
            if diffvalneg > 0:
                # logging.debug("Positive neg val")
                diff = negdur - ps

                return knotatleff(
                    l,
                    s,
                    ps + diff / 10,
                    posdur,
                    steps,
                    f0,
                    ngte,
                    kgte,
                    mu,
                    True,
                    ps=ps,
                    negcounter=negcounter + 1,
                )

    # Calculating the diff function value for the time between negdur and posdur
    halfdur = (posdur + negdur) / 2
    bf.knotslist.append(halfdur)
    basiscoeffs = bf.allcoeffs(s)
    if fullMOLS:
        modelparams = mols.solsbyint(mols.sols(basiscoeffs, 20, s, f0, ngte, kgte), s)
    else:
        modelparams = mols.solsbyint(
            mols.solsbetweenknots(
                knotnuma, knotnumb, basiscoeffs, 20, s, f0, ngte, kgte
            ),
            s,
        )
    paramspacing = metricelementspacingonl(l, basiscoeffs, mu)

    errorlist = []

    for t in np.linspace(bf.knotslist[-1], halfdur, steps):
        errorlist.append(
            mols.errorvalueatpoint(t, basiscoeffs, modelparams, f0, ngte, kgte)
        )

    bf.knotslist.pop()

    maxerror = np.max(errorlist)

    diffval = maxerror - paramspacing

    # Either return halfdur as our new knot value or recursing
    if -(2**-10) < diffval <= 0:
        parameterspacings.append(paramspacing)
        return halfdur
    elif recursiondepth > 10:
        logging.debug(
            "Recursion depth of 10 reached, terminating and returning current knot value"
        )
        logging.debug("Function value for this knot is " + str(diffval))
        parameterspacings.append(paramspacing)
        return halfdur
    else:
        if diffval < 0:
            return knotatleff(
                l,
                s,
                halfdur,
                posdur,
                steps,
                f0,
                ngte,
                kgte,
                mu,
                False,
                recursiondepth + 1,
            )
        elif diffval > 0:
            return knotatleff(
                l,
                s,
                negdur,
                halfdur,
                steps,
                f0,
                ngte,
                kgte,
                mu,
                False,
                recursiondepth + 1,
            )


# Same as the knotatleff method but instead of refering to the bf.knotslist parameter we instead include this list of
# previous knots as the parameter knotslist. The motivation behind this is that previous methods used to generate
# knots had to start from an initial knot p_0 = 0. This was inconvenient if we instead wanted to calculate p_100 and
# had already calculated p_99. Instead we carry all generated knots as a parameter in this function. We of course do
# not need to carry the entire knot list for the model specifications, only the last knot in that list,
# but we do this anyway.
#
# The logic in this method is identical to that of the knotatleff method with only the addition of the knotslist
# parameter. Read that method's description for more details about the process of this method. This method has
# a maximum recursion depth of 10.
#
# The fullMOLS parameter allows us to decide when considering the worse case scenario whether we want to consider the
# accuracy of the entire MOLS solution or just the MOLS solution which belongs to the segment we are currently trying
# to determine the knot for. Typically, if we consider the MOLS solution for the full signal, the knots are spread
# further apart.
fullMOLS = True


def knotatleffusinggivenknotlist(
    l,
    s,
    negdur,
    posdur,
    steps,
    f0,
    ngte,
    kgte,
    knotslist,
    mu,
    checkfirst=False,
    counter=0,
    ps=None,
    negcounter=0,
    prevdiffvalneg=np.inf,
):
    # logging.debug("Current recursion depth: " + str(counter))

    knotnuma = len(knotslist) - 1
    knotnumb = knotnuma + 1

    # As negdur is a candidate knot, it should be greater than the last knot we currently have in bf.knotslist.
    bf.knotslist = knotslist
    if knotslist[l] >= negdur:
        logging.debug(
            "negdur parameter not within appropriate range, adjusting with larger value"
        )
        logging.debug(
            "negdur = " + str(negdur) + ", previous knot = " + str(knotslist[l])
        )
        return knotatleffusinggivenknotlist(
            l,
            s,
            knotslist[l] + 1,
            posdur,
            steps,
            f0,
            ngte,
            kgte,
            knotslist,
            mu,
            True,
            ps=ps,
        )

    if checkfirst:
        bf.knotslist.append(negdur)
        basiscoeffsneg = bf.allcoeffs(s)

        # Can't remember why but this sometimes throws a TypeError. I suspect it probably happens when negdur is too
        # large and the sample points in the sampling_methods module are chosen poorly. Corrected for by using a
        # smaller negdur value. Don't ever remember seeing this error for the posdur parameter though
        try:
            if fullMOLS:
                modelparamsneg = mols.solsbyint(
                    mols.sols(basiscoeffsneg, 20, s, f0, ngte, kgte), s
                )
            else:
                modelparamsneg = mols.solsbyint(
                    mols.solsbetweenknots(
                        knotnuma, knotnumb, basiscoeffsneg, 20, s, f0, ngte, kgte
                    ),
                    s,
                )
        except TypeError as error:
            logging.debug(error)
            logging.debug("Reducing negdur duration")
            bf.knotslist.pop()
            diff = negdur - ps
            return knotatleffusinggivenknotlist(
                l,
                s,
                ps + diff / 10,
                posdur,
                steps,
                f0,
                ngte,
                kgte,
                knotslist,
                mu,
                True,
                0,
                ps=ps,
            )

        # Calculating the diff function value for negdur and posdur
        paramspacingneg = metricelementspacingonl(l, basiscoeffsneg, mu)

        errorlistneg = []

        for t in np.linspace(knotslist[-1], negdur, steps):
            errorlistneg.append(
                mols.errorvalueatpoint(
                    t, basiscoeffsneg, modelparamsneg, f0, ngte, kgte
                )
            )

        bf.knotslist.pop()

        bf.knotslist.append(posdur)
        basiscoeffspos = bf.allcoeffs(s)

        if fullMOLS:
            modelparamspos = mols.solsbyint(
                mols.sols(basiscoeffspos, 20, s, f0, ngte, kgte), s
            )
        else:
            modelparamspos = mols.solsbyint(
                mols.solsbetweenknots(
                    knotnuma, knotnumb, basiscoeffspos, 20, s, f0, ngte, kgte
                ),
                s,
            )

        paramspacingpos = metricelementspacingonl(l, basiscoeffspos, mu)

        errorlistpos = []

        for t in np.linspace(knotslist[-1], posdur, steps):
            errorlistpos.append(
                mols.errorvalueatpoint(
                    t, basiscoeffspos, modelparamspos, f0, ngte, kgte
                )
            )

        bf.knotslist.pop()

        maxerrorneg = np.max(errorlistneg)
        maxerrorpos = np.max(errorlistpos)

        diffvalneg = maxerrorneg - paramspacingneg
        diffvalpos = maxerrorpos - paramspacingpos

        # Checking if posdur and negdur have the appropriate signs for the bisection method and readjusting if not
        if not (diffvalneg < 0 and diffvalpos > 0):
            if diffvalpos < 0:
                diff = posdur - ps
                return knotatleffusinggivenknotlist(
                    l,
                    s,
                    negdur,
                    ps + 2 * diff,
                    steps,
                    f0,
                    ngte,
                    kgte,
                    knotslist,
                    mu,
                    True,
                    ps=ps,
                    prevdiffvalneg=diffvalneg,
                )
            if diffvalneg > 0:
                diff = negdur - ps

                if diffvalneg > prevdiffvalneg:
                    return knotatleffusinggivenknotlist(
                        l,
                        s,
                        ps + 2 * diff,
                        posdur,
                        steps,
                        f0,
                        ngte,
                        kgte,
                        knotslist,
                        mu,
                        True,
                        ps=ps,
                        negcounter=negcounter + 1,
                        prevdiffvalneg=diffvalneg,
                    )

                # if negcounter > 5 or (diff/10) < 1:
                #    logging.debug("Too many negdur reductions")
                #    logging.debug()
                #    raise errors.TooManyNegdurReductions

                return knotatleffusinggivenknotlist(
                    l,
                    s,
                    ps + diff / 10,
                    posdur,
                    steps,
                    f0,
                    ngte,
                    kgte,
                    knotslist,
                    mu,
                    True,
                    ps=ps,
                    negcounter=negcounter + 1,
                    prevdiffvalneg=diffvalneg,
                )

    # Calculating the diff function value for the time between negdur and posdur
    halfdur = (posdur + negdur) / 2
    bf.knotslist.append(halfdur)
    basiscoeffs = bf.allcoeffs(s)

    if fullMOLS:
        modelparams = mols.solsbyint(mols.sols(basiscoeffs, 20, s, f0, ngte, kgte), s)
    else:
        modelparams = mols.solsbyint(
            mols.solsbetweenknots(
                knotnuma, knotnumb, basiscoeffs, 20, s, f0, ngte, kgte
            ),
            s,
        )

    paramspacing = metricelementspacingonl(l, basiscoeffs, mu)

    errorlist = []

    for t in np.linspace(knotslist[-1], halfdur, steps):
        errorlist.append(
            mols.errorvalueatpoint(t, basiscoeffs, modelparams, f0, ngte, kgte)
        )
    bf.knotslist.pop()
    maxerror = np.max(errorlist)

    diffval = maxerror - paramspacing

    # Either return halfdur as our new knot value or recurse
    if -(2**-10) < diffval <= 0:
        parameterspacings.append(paramspacing)
        return halfdur
    elif counter > 10:
        logging.debug(
            "Recursion depth of 10 reached, terminating and returning current knot value"
        )
        logging.debug("Function value for this knot is " + str(diffval))
        parameterspacings.append(paramspacing)
        return halfdur
    else:
        if diffval < 0:
            return knotatleffusinggivenknotlist(
                l,
                s,
                halfdur,
                posdur,
                steps,
                f0,
                ngte,
                kgte,
                knotslist,
                mu,
                counter=counter + 1,
            )
        elif diffval > 0:
            return knotatleffusinggivenknotlist(
                l,
                s,
                negdur,
                halfdur,
                steps,
                f0,
                ngte,
                kgte,
                knotslist,
                mu,
                counter=counter + 1,
            )


# knotlist = [0, 2.107421875, 4.18438720703125]
# ints = len(knotlist)
# logging.debug(knotatleffusinggivenknotlist(ints - 1, 1, ints, 1, 10, 30, knotlist, checkfirst=True, ps=knotlist[-1]))

# Calculates and adds the next knot to the KnotArchive file with the appropriate specifications by using the
# knotatleffusinggivenknotlist method
def addnextknottofile(f0, nmax, kgte, s, mu, steps=30):
    spindownspecification = [f0, nmax, kgte, s, mu, fullMOLS]

    knotarchive = open(knotarchivefile, "r")
    alllines = knotarchive.readlines()
    knotarchive.close()

    foundspecifications = False
    newknot = 0

    for i, line in enumerate(alllines):
        thisline = ast.literal_eval(line)

        if thisline[0] == spindownspecification:
            foundspecifications = True

            knotslist = thisline[1]
            ps = thisline[1][-1]
            ints = len(knotslist)

            if ps == 0:
                negdurincrement = 1.1
                while True:
                    try:
                        newknot = knotatleffusinggivenknotlist(
                            ints - 1,
                            s,
                            negdurincrement,
                            10 * negdurincrement,
                            steps,
                            f0,
                            nmax,
                            kgte,
                            knotslist,
                            mu,
                            True,
                            ps=ps,
                        )

                        break
                    except errors.TooManyNegdurReductions:
                        negdurincrement *= 10
            else:
                negdurincrement = 0
                while True:
                    ordmag = np.floor(np.log10(ps)) + negdurincrement
                    try:
                        newknot = knotatleffusinggivenknotlist(
                            ints - 1,
                            s,
                            ps + 10**ordmag,
                            ps + 10 ** (ordmag + 0.5),
                            steps,
                            f0,
                            nmax,
                            kgte,
                            knotslist,
                            mu,
                            True,
                            ps=ps,
                        )
                        break
                    except errors.TooManyNegdurReductions:
                        negdurincrement += 0.5
            thisline[1].append(newknot)
            updatedline = str(thisline) + "\n"
            alllines[i] = updatedline
            logging.info("Knotnum and knot val: " + str(ints) + ", " + str(newknot))
            break

    if foundspecifications:
        knotarchive = open(knotarchivefile, "w")
        knotarchive.writelines(alllines)
        knotarchive.close()
    else:
        newline = [spindownspecification, [0]]

        knotarchive = open(knotarchivefile, "a")
        knotarchive.write(str(newline) + "\n")
        knotarchive.close()

        addnextknottofile(f0, nmax, kgte, s, mu, steps=steps)

    return newknot


# Calculates and writes all knots up the knotnum th knot to the KnotArchive file. E.g. If we want to know the first 10
# knots for certain model specifications, if we run this method with knotnum = 10, if the KnotArchive file does not
# already have those 10 knot values, this method will calculate and add those knot values to that file until all knots
# up to the 10th knot are present
def addknottoknotnum(s, knotnum, f0, nmax, kgte, mu, steps=30):
    spindownspecifications = [f0, nmax, kgte, s, mu, fullMOLS]

    knotarchive = open(knotarchivefile, "r")
    alllines = knotarchive.readlines()
    knotarchive.close()

    currentknotnum = 0

    for i, line in enumerate(alllines):
        thisline = ast.literal_eval(line)

        if thisline[0] == spindownspecifications:
            currentknotnum = len(thisline[1]) - 1
            break

    while currentknotnum <= knotnum:
        addnextknottofile(f0, nmax, kgte, s, mu, steps=steps)
        currentknotnum += 1


# Ad the above method but now instead of calculating up to a given knotnumber instead calculates all knots up to a given
# signal duration.
def addknottodur(s, dur, f0, nmax, kgte, mu, steps=30):
    spindownspecifications = [f0, nmax, kgte, s, mu, fullMOLS]

    knotarchive = open(knotarchivefile, "r")
    alllines = knotarchive.readlines()
    knotarchive.close()

    currentdur = 0

    for i, line in enumerate(alllines):
        thisline = ast.literal_eval(line)

        if thisline[0] == spindownspecifications:
            currentdur = thisline[1][-1]
            break

    while currentdur < dur:
        currentdur = addnextknottofile(f0, nmax, kgte, s, mu, steps=steps)


# Calculates and returns all knots either up to the time dur or given knot number from the KnotArchive file. By using
# this method we no longer need to always recalculate our knots, hopefully speeding up some of our other modules.
def allidealisedknots(s, dur, steps, f0, nmax, kgte, mu, knotnum=0):
    spindownspecifications = [f0, nmax, kgte, s, mu, fullMOLS]

    if os.path.exists(knotarchivefile):
        pass
    else:
        with open(knotarchivefile, "w") as _:
            pass

    if knotnum != 0:
        addknottoknotnum(s, knotnum, f0, nmax, kgte, mu, steps=steps)
    else:
        addknottodur(s, dur, f0, nmax, kgte, mu, steps=steps)

    knotarchive = open(knotarchivefile, "r")
    alllines = knotarchive.readlines()
    knotarchive.close()

    for i, line in enumerate(alllines):
        thisline = ast.literal_eval(line)

        if thisline[0] == spindownspecifications:
            bf.knotslist = thisline[1]

            knotsbelowdur = []

            if knotnum != 0:
                knotsbelowdur = thisline[1][0 : knotnum + 1]
            else:
                # In case the knots stored in KnotArchive go well past the specified duration, we make sure we only set the knots to go to the duration we require
                for i, knot in enumerate(bf.knotslist):
                    if knot < dur:
                        knotsbelowdur.append(knot)
                    else:
                        knotsbelowdur.append(dur)
                        break

            bf.knotslist = knotsbelowdur

            return bf.knotslist
