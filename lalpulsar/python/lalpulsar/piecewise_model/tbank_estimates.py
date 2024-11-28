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
Function to estimate the size of piecewise model template banks.
"""

import logging
import time

import numpy as np

import lal
import lalpulsar as lp

from . import basis_functions as bf
from . import semicoherent_metric_methods as scmm


# Returns tiling statistics from a tiling object
def tilingstatistics(tiling, dim, iterator=-1):

    if iterator == 0:
        totalpoints = []
        minpoints = []
        maxpoints = []
        minvals = []
        maxvals = []

        for i in range(dim):
            logging.debug("Dim :" + str(i))
            tilingstats = lp.LatticeTilingStatistics(tiling, i)
            totalpoints.append(tilingstats.total_points)
            minpoints.append(tilingstats.min_points)
            maxpoints.append(tilingstats.max_points)
            minvals.append(tilingstats.min_value)
            maxvals.append(tilingstats.max_value)

        logging.info("Total points up to dimension: %s", str(totalpoints))
        logging.info("Min points in dimension: %s", str(minpoints))
        logging.info("Max points in dimension: %s", str(maxpoints))
        logging.info("Min value in dimension: %s", str(minvals))
        logging.info("Max value in dimension: %s", str(maxvals))

        alltilingstats = [totalpoints, minpoints, maxpoints, minvals, maxvals]

        return alltilingstats

    else:
        return []


# Sets the bounds on the tiling lattice by the parameters given in a TBank object
def setbounds(tiling, tbank):
    s = tbank.s
    fmin = tbank.fmin
    fmax = tbank.fmax
    nmin = tbank.nmin
    nmax = tbank.nmax
    kmin = tbank.kmin
    kmax = tbank.kmax
    dur = tbank.dur
    mismatch = tbank.maxmismatch

    bbox = tbank.flags_bbox
    intbox = tbank.flags_intbox

    knots = bf.knotslist

    lp.SetLatticeTilingPiecewiseBounds(
        tiling, s, fmin, fmax, nmin, nmax, kmin, kmax, knots, bbox, intbox
    )


# Returns number of templates required to cover a parameter space with given knots
def recursivetbanksize(
    s,
    fmin,
    fmax,
    nmin,
    nmax,
    kmin,
    kmax,
    knots,
    mismatch,
    printouts=1000000,
    recdepth=1,
    prevalidtempsandnmin=[-1, -1],
    metric=[],
):

    filename = "BInd_" + str(fmin) + "_" + str(fmax) + ".txt"

    if nmin < nmin:
        return prevalidtempsandnmin

    finalknot = len(knots)

    # Create LatticeTiling object
    tiling = lp.CreateLatticeTiling(s * finalknot)

    bf.knotslist = knots
    logging.info("Knots list: %s", str(bf.knotslist))

    logging.info("Computing metric")

    if metric == []:
        metric = scmm.metric(s)

    logging.info("Metric calculated")

    stepsizes = [np.sqrt(2 / metric[i][i]) for i in range(s * finalknot)]
    logging.info("Maximum parameter step sizes: %s", str(stepsizes))

    # Set Bounds
    lp.SetLatticeTilingPiecewiseBounds(
        tiling, s, fmin, fmax, nmin, nmax, kmin, kmax, knots
    )

    # Set metric, mismatch and lattice type
    lp.SetTilingLatticeAndMetric(tiling, lp.TILING_LATTICE_ANSTAR, metric, mismatch)

    # Create Iterator
    iterator = lp.CreateLatticeTilingIterator(tiling, s * finalknot)
    lp.ResetLatticeTilingIterator(iterator)

    start = time.time()
    tiles = 0

    p = lal.gsl_vector(s * finalknot)

    while lp.NextLatticeTilingPoint(iterator, p) != 0:

        tiles += 1
        if tiles % printouts == 0:
            logging.info("Current number of tiles: %s", str(tiles))
            logging.info("Elapsed time: %s", str(time.time() - start))

        if tiles > 10**8:
            break

    if lp.NextLatticeTilingPoint(iterator, p) != 0:

        if recdepth >= 12:

            if prevalidtempsandnmin[0] != -1:
                logging.info(
                    "Final recursive step. Valid range found. Frequency range: %s, braking index range: %s",
                    str([fmin, fmax]),
                    str([nmin, nmax]),
                )

                f = open(filename, "a")
                f.write(str([nmin, tiles]) + "\n")
                f.close()

                return prevalidtempsandnmin

            else:
                logging.info(
                    "Final recdepth %s reached with no valid range found. Frequency range: %s, braking index range: %s",
                    str(recdepth),
                    str([fmin, fmax]),
                    str([nmin, nmax]),
                )

                f = open(filename, "a")
                f.write(str([nmin, tiles]) + "\n")
                f.close()

                return [-1, 10**8]

        else:
            logging.info(
                "Recursive step %s. No valid range found for this nmin. Frequency range: %s, braking index range: %s",
                str(recdepth),
                str([fmin, fmax]),
                str([nmin, nmax]),
            )

            f = open(filename, "a")
            f.write(str([nmin, tiles]) + "\n")
            f.close()

            newnmin = nmin + (nmax - nmin) * 2 ** -(recdepth + 1)
            return recursivetbanksize(
                s,
                fmin,
                fmax,
                newnmin,
                nmax,
                kmin,
                kmax,
                knots,
                mismatch,
                printouts=printouts,
                recdepth=recdepth + 1,
                metric=metric,
            )

    else:
        if recdepth >= 12:
            logging.info(
                "Statistics calculated: %s",
                str(
                    tilingstatistics(
                        tiling,
                        s * finalknot,
                        iterator=lp.NextLatticeTilingPoint(iterator, p),
                    )
                ),
            )
            logging.info(
                "Frequency range: %s, Braking index range: %s",
                str([fmin, fmax]),
                str([nmin, nmax]),
            )
            logging.info("Final tile count: %s", str(tiles))
            logging.info("Total elapsed time: %s", str(time.time() - start))

            logging.info(
                "Final recursive step. Valid range found. Frequency range: %s, braking index range: %s",
                str([fmin, fmax]),
                str([nmin, nmax]),
            )

            f = open(filename, "a")
            f.write(str([nmin, tiles]) + "\n")
            f.close()

            return [tiles, nmin]

        else:
            logging.info(
                "Statistics calculated: %s",
                str(
                    tilingstatistics(
                        tiling,
                        s * finalknot,
                        iterator=lp.NextLatticeTilingPoint(iterator, p),
                    )
                ),
            )
            logging.info(
                "Frequency range: %s, Braking index range: %s",
                str([fmin, fmax]),
                str([nmin, nmax]),
            )
            logging.info("Final tile count: %s", str(tiles))
            logging.info("Total elapsed time: %s", str(time.time() - start))

            logging.info(
                "Recursive step %s. Valid range found for this nmin. Frequency range: %s, braking index range: %s",
                str(recdepth),
                str([fmin, fmax]),
                str([nmin, nmax]),
            )

            f = open(filename, "a")
            f.write(str([nmin, tiles]) + "\n")
            f.close()

            newnmin = nmin - (nmax - nmin) * 2 ** -(recdepth + 1)
            return recursivetbanksize(
                s,
                fmin,
                fmax,
                newnmin,
                nmax,
                kmin,
                kmax,
                knots,
                mismatch,
                printouts=printouts,
                prevalidtempsandnmin=[tiles, nmin],
                recdepth=recdepth + 1,
                metric=metric,
            )
