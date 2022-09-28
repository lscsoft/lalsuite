import lalpulsar as lp
import lal

from . import EstimatingKnots as ek
from . import SemicoherentMetricMethods as scmm
from . import BasisFunctions as bf

import numpy as np
import time
import os
import logging
import cProfile, pstats, io

# Returns tiling statistics from a tiling object
def tilingstatistics(tiling, dim, iterator=-1):

        if iterator == 0:
                totalpoints = []
                minpoints = []
                maxpoints = []
                minvals = []
                maxvals = []

                for i in range(dim):
                        print("Dim :" + str(i))
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

# Writes tiling statistics to a given file
def tilingstatstofile(tilingstats, path, tbank):
        print("Path is: " + str(path))
        dimensions = tbank.s * len(bf.knotslist)

        if not os.path.exists(path):
                first_line_columns = "{:173}".format("Parameter Space Values") + "| "
                length_of_min_max_columns = "{:" + str(20 * dimensions) + "}"
                first_line_columns += length_of_min_max_columns.format("Total number of points in dimension") + "| "
                first_line_columns += length_of_min_max_columns.format("Minimum number of points in dimension") + "| "
                first_line_columns += length_of_min_max_columns.format("Maximum number of points in dimension") + "| "
                first_line_columns += length_of_min_max_columns.format("Minimum Values") + "| "
                first_line_columns += length_of_min_max_columns.format("Maximum Values") + "| "
                first_line_columns += "\n"

                second_line_columns = ""

                tbankstr = tbank.toString()
                splitstr = tbankstr.split()

                # Column titles for the template bank
                second_line_columns += "{:^5}".format("S")
                second_line_columns += "{:^8}".format("fmin")
                second_line_columns += "{:^8}".format("fmax")
                second_line_columns += "{:^8}".format("fmaxtrue")
                second_line_columns += "{:^9}".format("nmin")
                second_line_columns += "{:^9}".format("nmax")
                second_line_columns += "{:^9}".format("nmin0")
                second_line_columns += "{:^9}".format("nmax0")
                second_line_columns += "{:^20}".format("ntol")
                second_line_columns += "{:^20}".format("taumin")
                second_line_columns += "{:^20}".format("taumax")
                second_line_columns += "{:^20}".format("ktol")
                second_line_columns += "{:^8}".format("mismatch")
                second_line_columns += "{:^20}".format("dur")
                second_line_columns += "| "

                # Column titles for the total points in each dimension
                for i in range(dimensions):
                        second_line_columns += "{:^20.10G}".format(i)

                second_line_columns += "| "

                # Column Titles for min and max points of statistics
                for i in range(dimensions):
                        second_line_columns += "{:^20.10G}".format(i)

                second_line_columns += "| "

                for i in range(dimensions):
                        second_line_columns += "{:^20.10G}".format(i)

                second_line_columns += "| "

                # Column Titles for the min values and max values of statistics
                for i in range(dimensions):
                        second_line_columns += "{:^20.10G}".format(i)

                second_line_columns += "| "

                for i in range(dimensions):
                        second_line_columns += "{:^20.10G}".format(i)

                second_line_columns += "| "
                second_line_columns += "\n"

                new_file = open(path, 'a')
                new_file.write(first_line_columns)
                new_file.write(second_line_columns)
                new_file.close()

        new_statistics = ""
        new_statistics += "{:^5}".format(tbank.s)
        new_statistics += "{:^8}".format(tbank.fmin)
        new_statistics += "{:^8}".format(tbank.fmax)
        new_statistics += "{:^8}".format(tbank.fmaxtrue)
        new_statistics += "{:^9.3f}".format(tbank.nmin)
        new_statistics += "{:^9.3f}".format(tbank.nmax)
        new_statistics += "{:^9.3f}".format(tbank.nmin0)
        new_statistics += "{:^9.3f}".format(tbank.nmax0)
        new_statistics += "{:^20.10G}".format(tbank.ntol)
        new_statistics += "{:^20}".format(tbank.taumin)
        new_statistics += "{:^20}".format(tbank.taumax)
        new_statistics += "{:^20.10G}".format(tbank.ktol)
        new_statistics += "{:^8.2f}".format(tbank.mismatch)
        new_statistics += "{:^20}".format(tbank.dur)
        new_statistics += "| "

        for total_points in tilingstats[0]:
                new_statistics += "{:^20.10G}".format(total_points)

        new_statistics += "| "

        for min_points in tilingstats[1]:
                new_statistics += "{:^20.10G}".format(min_points)

        new_statistics += "| "

        for max_points in tilingstats[2]:
                new_statistics += "{:^20.10G}".format(max_points)

        new_statistics += "| "

        for min_vals in tilingstats[3]:
                new_statistics += "{:^20.10G}".format(min_vals)

        new_statistics += "| "

        for max_vals in tilingstats[4]:
                new_statistics += "{:^20.10G}".format(max_vals)

        new_statistics += "| \n"

        stats_file = open(path,'a');

        stats_file.write(new_statistics)
        stats_file.close()

# Returns number of templates required to cover a parameter space with given knots
def pwboundtbanksizecustomknots(s, fmin, fmax, fmaxtrue, nmin, nmax, nmin0, nmax0, ntol, taumin, taumax, ktol, knots, mismatch, reset, maxtemps=10000000, printouts=1000000, stats=False, dirname=None, filename=None, tbank=None, padding_flags_bbox=[], padding_flags_int=[]):

        if padding_flags_bbox == []:
                padding_flags_bbox = [0] * s * len(bf.knotslist)
        if padding_flags_int == []:
                padding_flags_int  = [0] * s * len(bf.knotslist)

        finalknot = len(knots)

        # Create LatticeTiling object
        tiling = lp.CreateLatticeTiling(s * finalknot)

        bf.knotslist = knots
        logging.info("Knots list: %s", str(bf.knotslist))

        logging.info("Computing metric")
        metric = scmm.metric(s)
        logging.info("Metric calculated")

        stepsizes = [np.sqrt(mismatch / metric[i][i]) for i in range(s * finalknot)]
        logging.info("Maximum parameter step sizes: %s", str(stepsizes))

        # Set Bounds
        if s == 3:
                lp.SetLatticeTilingPiecewiseBounds(  tiling, fmin, fmax, fmaxtrue, nmin, nmax, nmin0, nmax0, ntol, taumin, taumax, ktol, knots, finalknot, padding_flags_bbox, padding_flags_int, reset)
        elif s == 2:
                lp.SetLatticeTilingPiecewiseBoundsS2(tiling, fmin, fmax, fmaxtrue, nmin, nmax, nmin0, nmax0, ntol, taumin, taumax, ktol, knots, finalknot, padding_flags_bbox, padding_flags_int, reset)

        # Set metric, mismatch and lattice type
        lp.SetTilingLatticeAndMetric(tiling, lp.TILING_LATTICE_ANSTAR, metric, mismatch)

        # Create Iterator
        iterator = lp.CreateLatticeTilingIterator(tiling, s * finalknot)

        # Calculating template banke size. Two options, we calculatetemplate bank size using the lp.LatticeTilingStatistics method, however using this we do not get a regular print out
        # of how many templates we are currently at. The other option is we do not calculate the statistics, but a message is logged every 'printouts' number of templates to keep track
        # of where we are currently at.

        start_time = time.time()

        if stats:
                # In this section we calculate the entire number of templates required to cover a parameter space using lp.LatticeTilingStatistics. No prinouts are given to indicate
                # how big the template bank is as we tile the parameter space

                start = time.time()
                t = time.localtime()
                tilingstats = tilingstatistics(tiling, s * finalknot, iterator=0)
                logging.info("Time elapsed: %s", str(time.time() - start))

                if not dirname:
                        logging.warning("No directory given, saving generated stats to Stats_Dir")
                        dirname = "Stats_Dir"

                try:
                        os.mkdir(dirname)
                except FileExistsError:
                        pass

                if filename:
                        tilingstatstofile(tilingstats, [fmin, fmax], dirname + "/" + filename)
                if not filename and tbank:
                        dimensions = str(tbank.s * len(bf.knotslist))
                        filename = "Stats_Archive_With_" + dimensions + "_Dimensions.txt"
                else:
                        filename = "Default_file_name_" + str(fmin) + "_" + str(fmax) + ".txt"
                        logging.warning("No file name given, saving to file with name: " + str(filename))

                tilingstatstofile(tilingstats, dirname + "/" + filename, tbank)

                total_temps = tilingstats[0][-1]

                print("Templates per second: " + str(total_temps / (time.time() - start_time)))
                logging.info("Templates per second: " + str(total_temps / (time.time() - start_time)))

                return tilingstats[0][-1]

        else:
                # In this section we tile the parameter space without calculating any statistics. A message is logged every 'printouts' number of templates to indicate
                # how many templates are currently used to tile the parameter space.
                start = time.time()
                tiles = 0
                p = lal.gsl_vector(s * finalknot)
                while lp.NextLatticeTilingPoint(iterator, p) != 0:
                        tiles += 1

                        if tiles % printouts == 0:
                                logging.info("Current number of tiles: %s", str(tiles))
                                logging.info("Elapsed time: %s", str(time.time() - start))

                                print("Current number of tiles: %s", str(tiles))
                                print("Elapsed time: %s", str(time.time() - start))

                        if tiles > maxtemps:
                                logging.info("Maximum number of specified templates has been reached: %s", str(tiles))
                                break

                #if lp.NextLatticeTilingPoint(iterator, p) == 0:
                #       logging.info("Statistics calculated: %s", str(tilingstatistics(tiling, s * finalknot, iterator=lp.NextLatticeTilingPoint(iterator, p))))

                logging.info("Final tile count: %s", str(tiles))
                logging.info("Total elapsed time: %s", str(time.time() - start))

                print("Templates per second: " + str(tiles / (time.time() - start_time)))
                logging.info("Templates per second: " + str(tiles / (time.time() - start_time)))

                return tiles

# Returns number of templates required to cover a parameter space using knot algorithm
def pwboundtbanksize(s, fmin, fmax, fmaxtrue, nmin, nmax, nmin0, nmax0, ntol, taumin, taumax, ktol, dur, mismatch, reset, printouts=0, stats=False, dirname=None, filename=None):

        ek.allidealisedknots(s, dur, 40, fmaxtrue, nmax, taumin, mu=mismatch)
        knotnum = ek.getknotnum(s, dur, fmaxtrue, nmax, taumin, mu=mismatch)
        knots = bf.knotslist[:knotnum + 1]
        # This way our signal length is exactly dur and we don't have our last segment extending potentiall well beyond what we need it to be, and hence increasing computational cost
        knots[-1] = dur

        return pwboundtbanksizecustomknots(s, fmin, fmax, fmaxtrue, nmin, nmax, nmin0, nmax0, ntol, taumin, taumax, ktol, knots, mismatch, reset, printouts=printouts, stats=stats, dirname=dirname, filename=filename)

# Returns the values of the centres of each template required to cover a given parameter space. Will only return up to maxtemps
def pwboundtempscustomknots(s, fmin, fmax, fmaxtrue, nmin, nmax, nmin0, nmax0, ntol, taumin, taumax, ktol, knots, mismatch, maxtemps=10000, padding_flags_bbox=[], padding_flags_int=[]):

        if padding_flags_bbox == []:
                padding_flags_bbox = [0] * s * len(bf.knotslist)
        if padding_flags_int == []:
                padding_flags_int  = [0] * s * len(bf.knotslist)

        finalknot = len(knots)

        # Create LatticeTiling object
        tiling = lp.CreateLatticeTiling(s * finalknot)

        bf.knotslist = knots
        logging.info("Knots list: %s", str(bf.knotslist))

        metric = scmm.metric(s)

        # Set Bounds
        if s == 3:
                lp.SetLatticeTilingPiecewiseBounds(  tiling, fmin, fmax, fmaxtrue, nmin, nmax, nmin0, nmax0, ntol, taumin, taumax, ktol, knots, finalknot, padding_flags_bbox, padding_flags_int)
        elif s == 2:
                lp.SetLatticeTilingPiecewiseBoundsS2(tiling, fmin, fmax, fmaxtrue, nmin, nmax, nmin0, nmax0, ntol, taumin, taumax, ktol, knots, finalknot, padding_flags_bbox, padding_flags_int)

        # Set metric, mismatch and lattice type
        lp.SetTilingLatticeAndMetric(tiling, lp.TILING_LATTICE_ANSTAR, metric, mismatch)

        # Create Iterator
        iterator = lp.CreateLatticeTilingIterator(tiling, s * finalknot)

        p = lal.gsl_vector(s * finalknot)

        temps = []

        for i in range(0, maxtemps):
                fin = lp.NextLatticeTilingPoint(iterator, p)
                temps.append((p.data.copy()))

                if fin == 0:
                        break

        return temps

# Returns the templates to cover a given parameter space with knots built from the knot algorithm
def pwboundtemplates(s, fmin, fmax, fmaxtrue, nmin, nmax, nmin0, nmax0, ntol, taumin, taumax, ktol, dur, mismatch, maxtemps=1000, padding_flags_bbox=[], padding_flags_int=[]):

        ek.allidealisedknots(s, dur, 30, fmaxtrue, nmax, taumin, mu=mismatch)
        knotnum = ek.getknotnum(s, dur, fmaxtrue, nmax, taumin, mu=mismatch)
        knots = bf.knotslist[:knotnum + 1]
        knots[-1] = dur

        return pwboundtempscustomknots(s, fmin, fmax, fmaxtrue, nmin, nmax, nmin0, nmax0, ntol, taumin, taumax, ktol, knots, mismatch, maxtemps=maxtemps, padding_flags_bbox=padding_flags_bbox, padding_flags_int=padding_flags_int)

# Sets the bounds on the tiling lattice by the parameters given in a TBank object
def setbounds(tiling, padding_flags_bbox, padding_flags_int, tbank, reset):
        s        = tbank.s
        fmin     = tbank.fmin
        fmax     = tbank.fmax
        fmaxtrue = tbank.fmaxtrue
        nmin     = tbank.nmin
        nmax     = tbank.nmax
        nmin0    = tbank.nmin0
        nmax0    = tbank.nmax0
        ntol     = tbank.ntol
        taumin   = tbank.taumin
        taumax   = tbank.taumax
        ktol     = tbank.ktol
        dur      = tbank.dur
        mismatch = tbank.mismatch
        maxtemps = tbank.maxtemps

        knots = bf.knotslist

        finalknot = len(knots)

        flagoption = tbank.flagoption

        if s == 3:
                lp.SetLatticeTilingPiecewiseBounds(  tiling, fmin, fmax, fmaxtrue, nmin, nmax, nmin0, nmax0, ntol, taumin, taumax, ktol, knots, finalknot, padding_flags_bbox, padding_flags_int, reset)
        elif s == 2:
                lp.SetLatticeTilingPiecewiseBoundsS2(tiling, fmin, fmax, fmaxtrue, nmin, nmax, nmin0, nmax0, ntol, taumin, taumax, ktol, knots, finalknot, padding_flags_bbox, padding_flags_int)

# Calculates the size of a template bank for a TBank object
def PWTBankSizeWithObject(tbank, reset, stats=False, padding_flags_bbox=[], padding_flags_int=[]):

        pr = cProfile.Profile()
        pr.enable()

        s = tbank.s
        fmin = tbank.fmin
        fmax = tbank.fmax
        fmaxtrue = tbank.fmaxtrue
        nmin = tbank.nmin
        nmax = tbank.nmax
        nmin0 = tbank.nmin0
        nmax0 = tbank.nmax0
        ntol = tbank.ntol
        taumin = tbank.taumin
        taumax = tbank.taumax
        ktol = tbank.ktol
        mismatch = tbank.mismatch
        dur = tbank.dur
        maxtemps = tbank.maxtemps

        #ek.allidealisedknots(s, dur, 40, fmaxtrue, nmax, taumin, mismatch)
        #tbank.knots = bf.knotslist
        knots = tbank.knots

        temp = pwboundtbanksizecustomknots(s, fmin, fmax, fmaxtrue, nmin, nmax, nmin0, nmax0, ntol, taumin, taumax, ktol, knots, mismatch, reset, maxtemps=maxtemps, tbank=tbank, stats=stats, padding_flags_bbox=padding_flags_bbox, padding_flags_int=padding_flags_int)

        pr.disable()

        s = io.StringIO()
        sortby = pstats.SortKey.CUMULATIVE
        ps = pstats.Stats(pr, stream=s).sort_stats('cumtime')
        ps.print_stats()

        with open("TBankEstimatesProfile.txt", 'w+') as f:
                f.write(s.getvalue())

        return temp

# Gives the template bank for a given TBank object
def PWTBankWithObject(tbank, padding_flags_bbox=[], padding_flags_int=[]):
        s = tbank.s
        fmin = tbank.fmin
        fmax = tbank.fmax
        fmaxtrue = tbank.fmaxtrue
        nmin = tbank.nmin
        nmax = tbank.nmax
        nmin0 = tbank.nmin0
        nmax0 = tbank.nmax0
        ntol = tbank.ntol
        taumin = tbank.taumin
        taumax = tbank.taumax
        ktol = tbank.ktol
        mismatch = tbank.mismatch
        dur = tbank.dur
        maxtemps = tbank.maxtemps

        #ek.allidealisedknots(s, dur, 40, fmaxtrue, nmax, taumin, mismatch)
        #tbank.knots = bf.knotslist
        knots = tbank.knots

        temps = pwboundtempscustomknots(s, fmin, fmax, fmaxtrue, nmin, nmax, nmin0, nmax0, ntol, taumin, taumax, ktol, knots, mismatch, maxtemps=tbank.maxtemps, padding_flags_bbox=padding_flags_bbox, padding_flags_int=padding_flags_int)

        return temps

# Returns number of templates required to cover a parameter space with given knots
def recursivetbanksize(s, fmin, fmax, fmaxtrue, nmin, nmax, nmin0, nmax0, ntol, taumin, taumax, ktol, knots, mismatch, printouts=1000000, recdepth=1, prevalidtempsandnmin0=[-1, -1], metric=[]):

        filename = "BInd_" + str(fmin) + "_" + str(fmax) + ".txt"

        if nmin0 < nmin:
                return prevalidtempsandnmin0

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
        if s == 3:
                lp.SetLatticeTilingPiecewiseBounds(  tiling, fmin, fmax, fmaxtrue, nmin, nmax, nmin0, nmax0, ntol, taumin, taumax, ktol, knots, finalknot, padding_flags_bbox, padding_flags_int)
        elif s == 2:
                lp.SetLatticeTilingPiecewiseBoundsS2(tiling, fmin, fmax, fmaxtrue, nmin, nmax, nmin0, nmax0, ntol, taumin, taumax, ktol, knots, finalknot, padding_flags_bbox, padding_flags_int)

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

                if tiles > 10 ** 8:
                        break

        if lp.NextLatticeTilingPoint(iterator, p) != 0:

                if recdepth >= 12:

                        if prevalidtempsandnmin0[0] != -1:
                                logging.info("Final recursive step. Valid range found. Frequency range: %s, braking index range: %s", str([fmin, fmax]), str([nmin0, nmax0]))

                                f=open(filename,'a')
                                f.write(str([nmin0, tiles]) + "\n")
                                f.close()

                                return prevalidtempsandnmin0

                        else:
                                logging.info("Final recdepth %s reached with no valid range found. Frequency range: %s, braking index range: %s", str(recdepth), str([fmin, fmax]), str([nmin0, nmax0]))

                                f=open(filename,'a')
                                f.write(str([nmin0, tiles]) + "\n")
                                f.close()

                                return [-1, 10 ** 8]

                else:
                        logging.info("Recursive step %s. No valid range found for this nmin0. Frequency range: %s, braking index range: %s", str(recdepth), str([fmin, fmax]), str([nmin0, nmax0]))

                        f=open(filename,'a')
                        f.write(str([nmin0, tiles]) + "\n")
                        f.close()

                        newnmin0 = nmin0 + (nmax - nmin) * 2 ** -(recdepth + 1)
                        return recursivetbanksize(s, fmin, fmax, fmaxtrue, nmin, nmax, newnmin0, nmax0, ntol, taumin, taumax, ktol, knots, mismatch, printouts=printouts, recdepth=recdepth + 1, metric=metric)

        else:
                if recdepth >= 12:
                        logging.info("Statistics calculated: %s", str(tilingstatistics(tiling, s * finalknot, iterator=lp.NextLatticeTilingPoint(iterator, p))))
                        logging.info("Frequency range: %s, Braking index range: %s", str([fmin, fmax]), str([nmin0, nmax0]))
                        logging.info("Final tile count: %s", str(tiles))
                        logging.info("Total elapsed time: %s", str(time.time() - start))

                        logging.info("Final recursive step. Valid range found. Frequency range: %s, braking index range: %s", str([fmin, fmax]), str([nmin0, nmax0]))

                        f=open(filename,'a')
                        f.write(str([nmin0, tiles]) + "\n")
                        f.close()

                        return [tiles, nmin0]

                else:
                        logging.info("Statistics calculated: %s", str(tilingstatistics(tiling, s * finalknot, iterator=lp.NextLatticeTilingPoint(iterator, p))))
                        logging.info("Frequency range: %s, Braking index range: %s", str([fmin, fmax]), str([nmin0, nmax0]))
                        logging.info("Final tile count: %s", str(tiles))
                        logging.info("Total elapsed time: %s", str(time.time() - start))

                        logging.info("Recursive step %s. Valid range found for this nmin0. Frequency range: %s, braking index range: %s", str(recdepth), str([fmin, fmax]), str([nmin0, nmax0]))

                        f=open(filename,'a')
                        f.write(str([nmin0, tiles]) + "\n")
                        f.close()

                        newnmin0 = nmin0 - (nmax - nmin) * 2 ** -(recdepth + 1)
                        return recursivetbanksize(s, fmin, fmax, fmaxtrue, nmin, nmax, newnmin0, nmax0, ntol, taumin, taumax, ktol, knots, mismatch, printouts=printouts, prevalidtempsandnmin0=[tiles, nmin0], recdepth=recdepth + 1, metric=metric)

# Gives the template bank for a given TBank object
def PWTBankSizeWithObjectRecursive(tbank):
        s = tbank.s
        fmin = tbank.fmin
        fmax = tbank.fmax
        fmaxtrue = tbank.fmaxtrue
        nmin = tbank.nmin
        nmax = tbank.nmax
        nmin0 = tbank.nmin0
        nmax0 = tbank.nmax0
        ntol = tbank.ntol
        taumin = tbank.taumin
        taumax = tbank.taumax
        ktol = tbank.ktol
        mismatch = tbank.mismatch
        dur = tbank.dur
        maxtemps = tbank.maxtemps

        ek.allidealisedknots(s, dur, 40, fmaxtrue, nmax, taumin, mismatch)
        tbank.knots = bf.knotslist
        knots = tbank.knots

        return recursivetbanksize(s, fmin, fmax, fmaxtrue, nmin, nmax, nmin0, nmax0, ntol, taumin, taumax, ktol, knots, mismatch)

# Calculates the size of a template bank for a TBank object
def PWTBankSizeWithObjectEnMass(tbank, jobid, stats=False):

        freqs = np.linspace(100, 2000, 39)
        binds = [2, 4, 4.5, 4.75, 4.8, 4.85, 4.9, 4.95, 4.99, 4.995, 4.999, 5]

        freqranges = []
        bindranges = []

        for i in range(len(binds) - 1):
                thisrange = [binds[i], binds[i + 1]]
                bindranges.append(thisrange)

        freqsandbinds = []

        for freq in freqs:
                for bind in bindranges:
                        freqsandbinds.append([freq, bind])

        pr = cProfile.Profile()
        pr.enable()

        s = tbank.s
        fmin = freqsandbinds[jobid - 1][0]
        fmax = fmin + 0.01
        fmaxtrue = tbank.fmaxtrue
        nmin = tbank.nmin
        nmax = tbank.nmax
        nmin0 = freqsandbinds[jobid - 1][1][0]
        nmax0 = freqsandbinds[jobid - 1][1][1]
        ntol = tbank.ntol
        taumin = tbank.taumin
        taumax = tbank.taumax
        ktol = tbank.ktol
        mismatch = tbank.mismatch
        dur = tbank.dur
        maxtemps = tbank.maxtemps

        ek.allidealisedknots(s, dur, 40, fmaxtrue, nmax, taumin, mismatch)
        tbank.knots = bf.knotslist
        knots = tbank.knots

        temp = pwboundtbanksizecustomknots(s, fmin, fmax, fmaxtrue, nmin, nmax, nmin0, nmax0, ntol, taumin, taumax, ktol, knots, mismatch, stats=stats)

        filename = "TempsEnMass.txt"

        f=open(filename,'a')
        f.write(str([fmin, fmax, nmin0, nmax0, temp]) + "\n")
        f.close()

        pr.disable()

        s = io.StringIO()
        sortby = pstats.SortKey.CUMULATIVE
        ps = pstats.Stats(pr, stream=s).sort_stats('cumtime')
        ps.print_stats()

        with open("TBankEstimatesProfile.txt", 'w+') as f:
                f.write(s.getvalue())

        return temp

