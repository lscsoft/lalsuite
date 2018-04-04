# Copyright (C) 2012  Matthew West
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


#
# =============================================================================
#
#                                 Preamble
#
# =============================================================================
#

"""
Collection of functions to compute the efficiency and effective 4-volume
"""

import sqlite3
import math
import pdb
from operator import itemgetter
import numpy
from scipy import special

from glue import segments
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import dbtables

from pylal import antenna
from pylal import ligolw_sqlutils as sqlutils
from pylal import ligolw_cbc_compute_durations as compute_dur


#
# =============================================================================
#
#                                 
#
# =============================================================================
#


def chirp_dist(distance, mchirp):
        mchirp_DNS = (1.4+1.4) * (1./4)**(3.0/5.0)

        return distance * (mchirp_DNS/mchirp)**(5.0/6.0)

def decisive_dist(
    h_dist, l_dist, v_dist, 
    mchirp, weight_dist, ifos):
    
    dist_list = []
    if 'H1' in ifos or 'H2' in ifos:
        dist_list.append(h_dist)
    if 'L1' in ifos:
        dist_list.append(l_dist)
    if 'V1' in ifos:
        dist_list.append(v_dist) 

    if weight_dist:
        return chirp_dist(sorted(dist_list)[1], mchirp) 
    else:
        return sorted(dist_list)[1]

def end_time_with_ns(end_time, end_time_ns):
    time = end_time + 1e-9*end_time_ns
    return time

def get_livetime(connection, veto_cat, on_ifos, datatype):
    sqlquery = """
    SELECT duration
    FROM experiment_summary
        JOIN experiment ON (
            experiment_summary.experiment_id == experiment.experiment_id)
    WHERE
        datatype = :0
        AND veto_def_name = :1
        AND instruments = :2 """

    # total livetime in seconds 
    total_dur = numpy.sum(connection.execute(sqlquery, (datatype, veto_cat, on_ifos)).fetchall() )

    return total_dur

#
# =============================================================================
#
#                         Injections Functions
#
# =============================================================================
#

def inj_dist_range(dist_bounds, dist_scale = "linear", step = 4.0):

    if dist_scale == "linear":
        dist_bin_edges = numpy.arange(dist_bounds[0]-step, dist_bounds[1]+step, step)
    elif dist_scale == "log":
        log_limits = numpy.log10([dist_bounds[0], dist_bounds[1]])/numpy.log10(step)
        dist_bin_edges = numpy.power(
            step,
            numpy.arange(log_limits[0]-1, log_limits[1]+1)
        )

    return dist_bin_edges


def successful_injections(
    connection,
    tag,
    on_ifos,
    veto_cat,
    dist_type = "distance",
    weight_dist = False,
    verbose = False):

    """
    My attempt to get a list of the simulations that actually made
    it into some level of coincident time
    """

    xmldoc = dbtables.get_xml(connection)
    connection.create_function('end_time_with_ns', 2, end_time_with_ns)

    # Get the veto segments as dictionaries, keyed by veto category
    veto_segments = compute_dur.get_veto_segments(xmldoc, verbose)

    # ------------------------ Get List of Injections ------------------------ #
    sql_params_dict = {}
    sqlquery = """
        SELECT DISTINCT
            simulation_id,
            end_time_with_ns(geocent_end_time, geocent_end_time_ns),"""
    # add the desired distance measure to the SQL query
    if dist_type == "distance":
        connection.create_function('distance_func', 2, chirp_dist)
        sqlquery += """
            distance_func(distance, sim_inspiral.mchirp)
        FROM sim_inspiral """
    elif dist_type == "decisive_distance":
        connection.create_function('decisive_dist_func', 6, decisive_dist)
        sql_params_dict['ifos'] = on_ifos
        sql_params_dict['weight_dist'] = weight_dist
        sqlquery += """
            decisive_dist_func(
                eff_dist_h, eff_dist_l, eff_dist_v,
                sim_inspiral.mchirp, :weight_dist, :ifos)
        FROM sim_inspiral """

    if tag != 'ALL_INJ':
        # if a specific injection set is wanted
        sqlquery += """
        JOIN process_params ON (
            process_params.process_id == sim_inspiral.process_id)
        WHERE process_params.value = :usertag) """
        sql_params_dict["usertag"] = tag
    else:
        # for all injections
        tag = 'FULL_DATA'

    # Get segments that define which time was filtered
    ifo_segments = compute_dur.get_single_ifo_segments(
        connection,
        program_name = "inspiral",
        usertag = tag)

    zero_lag_dict = dict([(ifo, 0.0) for ifo in ifo_segments])

    successful_inj = []
    # determine coincident segments for that veto category 
    coinc_segs = compute_dur.get_coinc_segments(
        ifo_segments - veto_segments[veto_cat],
        zero_lag_dict)

    # Apply vetoes to single-ifo filter segments
    for injection in connection.execute(sqlquery, sql_params_dict):
        inj_segment = segments.segment(injection[1], injection[1])
        if coinc_segs[on_ifos].intersects_segment( inj_segment ):
            successful_inj.append( injection )

    return successful_inj


def found_injections(
    connection,
    tag,
    on_ifos,
    dist_type = "distance",
    weight_dist = False,
    verbose = False):

    connection.create_function('end_time_with_ns', 2, end_time_with_ns)

    sql_params_dict = {'ifos': on_ifos}
    sqlquery = """
    SELECT DISTINCT
        sim_inspiral.simulation_id,
        end_time_with_ns(geocent_end_time, geocent_end_time_ns), """
    # add the desired distance measure to the SQL query
    if dist_type == "distance":
        connection.create_function('distance_func', 2, chirp_dist)
        sqlquery += """
            distance_func(distance, sim_inspiral.mchirp), """
    elif dist_type == "decisive_distance":
        connection.create_function('decisive_dist_func', 6, decisive_dist)
        sql_params_dict['weight_dist'] = weight_dist
        sqlquery += """
            decisive_dist_func(
                eff_dist_h, eff_dist_l, eff_dist_v,
                sim_inspiral.mchirp, :weight_dist, :ifos), """

    sqlquery += """
        false_alarm_rate,
        coinc_inspiral.snr
    FROM
        coinc_event_map AS coincs
        JOIN coinc_event_map AS sims, coinc_inspiral, coinc_event, sim_inspiral ON (
            coincs.coinc_event_id == sims.coinc_event_id
            AND coinc_event.coinc_event_id == coincs.event_id
            AND coinc_inspiral.coinc_event_id == coincs.event_id
            AND sim_inspiral.simulation_id == sims.event_id)
        JOIN process_params ON (
            process_params.process_id == sim_inspiral.process_id)
    WHERE
        coincs.table_name = "coinc_event"
        AND sims.table_name = "sim_inspiral"
        AND coinc_event.instruments = :ifos """

    if tag != 'ALL_INJ':
        sqlquery += """
        AND process_params.value = :usertag
        """
        sql_params_dict["tag"] = tag

    injections = set(connection.execute(sqlquery, sql_params_dict).fetchall())

    # Get foreground coinc events
    sqlquery = """
        SELECT 
            end_time_with_ns(end_time, end_time_ns) AS trig_time,
            snr AS trig_snr
        FROM coinc_inspiral AS ci
            JOIN experiment_map AS em, experiment_summary AS es ON (
                ci.coinc_event_id == em.coinc_event_id
                AND em.experiment_summ_id == es.experiment_summ_id )
        WHERE es.datatype == 'all_data';
    """
    foreground = connection.executescript(sqlquery).fetchall()

    # Remove injection coincs that correspond closely in time and snr to foreground coincs 
    inj_snr = numpy.array([inj[3] for inj in injections])
    inj_time = numpy.array([inj[1] for inj in injections])
    idx2remove = []
    for time, snr in foreground:
        indices =  numpy.where(numpy.abs(inj_time - time) < 1.0)
        if len(indices[0]):
            idx = numpy.where(inj_snr[indices]/snr < 1.25)
            if len(idx[0]):
                idx2remove += list(indices[idx[0]])
    for i in sorted(idx2remove, reverse=True):
        del injections[i]
    
    # Sort found injections by FAR (largest to smallest)
    injections = sorted(injections, key=itemgetter(3), reverse=True)
    found_inj = [inj[0:3] for inj in injections]
    inj_fars = [inj[3] for inj in injections]
    inj_snrs = [inj[4] for inj in injections]

    return found_inj, inj_fars, inj_snrs


def binomial_confidence(K, N, eff_bin_edges, confidence):
    """
    Calculate the optimal Bayesian credible interval for p(eff|k,n)
    Posterior generated with binomial p(k|eff,n) and a uniform p(eff)
    is the beta function: Beta(eff|k+1,n-k+1) where n is the number
    of injected signals and k is the number of found signals.
    """
    eff_low = numpy.zeros(len(N))
    eff_high = numpy.zeros(len(N))
    for i, n in enumerate(N):
        if n!= 0:
             # construct the point-mass-function
             eff_cdf = special.betainc(K[i]+1, n-K[i]+1, eff_bin_edges)
             pmf = ( eff_cdf[1:] - eff_cdf[:-1] )

             # determine the indices for the highest density interval
             a = numpy.argsort(pmf)[::-1]
             sorted_cdf = numpy.cumsum(numpy.sort(pmf)[::-1])
             j = numpy.argmin( numpy.abs(sorted_cdf - confidence) )

             eff_low[i] = eff_bin_edges[:-1][ numpy.min(a[:(j+1)]) ]
             eff_high[i] = eff_bin_edges[:-1][ numpy.max(a[:(j+1)]) ]

    return eff_low, eff_high

def detection_efficiency(
    successful_inj,
    found_inj,
    found_fars,
    far_list,
    r,
    confidence):
    """
    This function determines the peak efficiency for a given bin and associated
    'highest density' confidence interval. The calculation is done for results
    from each false-alarm-rate threshold
    """
    # catching any edge cases were the injection end_time is nearly on a second boundary
    successful_inj = set(successful_inj) | set(found_inj)
    # histogram of successful injections into coincidence time post vetoes
    successful_dist = [inj[2] for inj in successful_inj]
    N, _ = numpy.histogram(successful_dist, bins = r)

    significant_dist = [inj[2] for inj in found_inj]

    eff_bin_edges = numpy.linspace(0, 1, 1e3+1)
    eff = {
        'mode': {},
        'low': {},
        'high': {}
    }
    for far in far_list:
        for idx, coinc_far in enumerate(found_fars):
            if coinc_far <= far:
                new_start = idx
                break
        # Histogram found injections with FAR < threshold
        K, _ = numpy.histogram(significant_dist[new_start:], bins = r)
        eff['mode'][far] = numpy.nan_to_num(numpy.float_(K)/N)

        # computes the confidence interval for each distance bin
        eff['low'][far], eff['high'][far] = binomial_confidence(K, N, eff_bin_edges, confidence)

    return eff


def rescale_dist(
    on_ifos, dist_type, weight_dist,
    phys_dist=None, param_dist=None
    ):

    N_signals = int(1e6)
    trigTime = 0.0

    # if decisive distance is desired, get the antenna responses for each signal
    if dist_type == 'decisive_distance':
        # sky position (right ascension & declination)
        ra = 360 * numpy.random.rand(N_signals)
        dec = 180 * numpy.random.rand(N_signals) - 90
        # additional angles
        inclination = 180 * numpy.random.rand(N_signals)
        polarization = 360 * numpy.random.rand(N_signals)
   
        f_q = {} 
        for ifo in on_ifos:
            f_q[ifo] = numpy.zeros(N_signals)
            for index in range(N_signals):
                _, _, _, f_q[ifo][index] = antenna.response(
                   trigTime,
                   ra[index], dec[index],
                   inclination[index], polarization[index],
                   'degree', ifo )
    
    prob_d_d = {}
    for j in range(len(phys_dist)-1):
        # for this physical distance range, create signals that are uniform in volume
        volume = 4*numpy.pi/3 * numpy.random.uniform(
            low = phys_dist[j]**3.0,
            high = phys_dist[j+1]**3.0,
            size = N_signals)
        dist = numpy.power(volume*(3/(4*numpy.pi)), 1./3)

        # create decisive distance (if desired)
        if dist_type == 'decisive_distance':
            dist_eff = {}
            for ifo in on_ifos:
                dist_eff[ifo] = dist / f_q[ifo]
            dist_dec = numpy.sort(dist_eff.values(), 0)[1]

        # weight distance measure by chirp mass (if desired)
        if weight_dist:
            # Component masses are Gaussian distributed around the Chandrasekar mass
            mass1, mass2 = 0.13 * numpy.random.randn(2, N_signals) + 1.40
            mchirp = numpy.power(mass1+mass2, -1./5) * numpy.power(mass1*mass2, 3./5)
            if dist_type == 'decisive_distance':
                dist_chirp = chirp_dist(dist_dec, mchirp)
            if dist_type == 'distance':
                dist_chirp = chirp_dist(dist, mchirp)
            N_d, _ = numpy.histogram(dist_chirp, bins=param_dist)
        else:
            N_d, _ = numpy.histogram(dist_dec, bins=param_dist)
    
        prob_d_d[phys_dist[j+1]] = numpy.float_(N_d)/numpy.sum(N_d)

    return prob_d_d

def eff_vs_dist(measured_eff, prob_dc_d):
    """
    This function creates a weighted average efficiency as a function of distance
    by computing eff_wavg(D) = \sum_dc eff_mode(dc)p(dc|d). p(dc|d) is the probability
    a signal has a parameterized distance dc if its physical distance is d.

    The confidence interval for eff_wavg(d) is constructed from the quadrature sum
    of the difference between the modes and the bounds, with each term again
    weighted by p(dc|d).

    This operation is done for each false-alarm-rate threshold.
    """
    eff_dist = {
        'wavg': {},
        'low': {},
        'high': {}
    }
    for far, modes in measured_eff['mode'].items():
        eff_dist['wavg'][far] = numpy.sum(modes * prob_dc_d.values(), 1)
        eff_dist['low'][far] = numpy.sqrt(numpy.sum(
            (measured_eff['low'][far] - modes)**2 * prob_dc_d.values(), 1)
        )
        eff_dist['high'][far] = numpy.sqrt(numpy.sum(
            (measured_eff['high'][far] - modes)**2 * prob_dc_d.values(), 1)
        )
    return eff_dist


def volume_efficiency(measured_eff, V_shell,  prob_dc_d):
    """
    This function creates a weighted average efficiency within a given volume
    by computing eff_wavg(D) = \sum_dc eff_mode(dc)p(dc|D). p(dc|D) is the
    probability a signal has a parameterized distance dc if it falls within
    physical distance D.

    The confidence interval for eff_wavg(D) is constructed from the quadrature sum
    of the difference between the modes and the bounds, with each term again
    weighted by p(dc|D).

    This operation is done for each false-alarm-rate threshold.
    """
    
    cumVol = numpy.cumsum(V_shell)
    p_dc_d = numpy.array(prob_dc_d.values())
    V_dc_D = numpy.cumsum((V_shell * p_dc_d.transpose()).transpose(), 0)
    p_dc_D = (V_dc_D.transpose()/cumVol).transpose()

    vol_eff = {
        'wavg': {},
        'low': {},
        'high': {}
    }
    for far, modes in measured_eff['mode'].items():
        vol_eff['wavg'][far] = numpy.sum(modes * p_dc_D, 1)
        vol_eff['low'][far] = numpy.sqrt(numpy.sum(
            (measured_eff['low'][far] - modes)**2 * p_dc_D, 1)
        )
        vol_eff['high'][far] = numpy.sqrt(numpy.sum(
            (measured_eff['high'][far] - modes)**2 * p_dc_D, 1)
        )

    return vol_eff
