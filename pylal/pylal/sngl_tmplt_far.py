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
A collection of utilities that allow one to estimate a false alarm rate from 
triggers of a single template without the need for doing time-slides.
"""

import sqlite3
import numpy as np
from scipy import stats
import bisect

from glue import iterutils
from glue import segments
from glue.ligolw import dbtables

from pylal import ligolw_sqlutils as sqlutils
from pylal import ligolw_cbc_compute_durations as compute_dur



#
# =============================================================================
#
#                                Library API 
#
# =============================================================================
#

def get_newsnr(snr, chisq, chisq_dof):
    fac = 6.0
    rchisq = chisq/(2*chisq_dof-2)
    newsnr = snr/(0.5*(1+rchisq**(fac/2.0)))**(1.0/fac)
    if rchisq < 1.0:
        newsnr = snr
    return newsnr

def get_effsnr(snr, chisq, chisq_dof):
    fac = 250.0
    rchisq = chisq/(2*chisq_dof-2)
    effsnr = snr/((1 + snr**2/fac)*rchisq)**(1./4)
    return effsnr

def get_snr(snr, chisq, chisq_dof):
    return snr

def get_snr_over_chi(snr, chisq, chisq_dof):
    return snr/chisq**(1./2)

def set_getsnr_function(connection, statistic):
    """
    Set the get_snr SQL function. Available options are:
    -- "rawsnr", "newsnr", "effsnr", "snroverchi"
    """
    if statistic == "rawsnr":
        connection.create_function('get_snr', 3, get_snr)
    elif statistic == "newsnr":
        connection.create_function('get_snr', 3, get_newsnr)
    elif statistic == "effsnr":
        connection.create_function('get_snr', 3, get_effsnr)
    elif statistic == "snroverchi":
        connection.create_function('get_snr', 3, get_snr_over_chi)
    else:
        raise ValueError, "%s is not a valid snr statistic" % statistic

def quadrature_sum(tuple):
    snrsq_array = np.array(tuple)**2.
    comb_snr = np.sum( snrsq_array )**(1/2.)
    return comb_snr

def end_time_w_ns(end_time, end_time_ns):
    time = end_time + 1e-9*end_time_ns
    return time

def compute_cumrate(hist, T):
    # create the survival function
    cum_hist = np.cumsum( hist[::-1] )[::-1]

    rate = {}
    # moments of rate pdf assuming that the bkgd coincs are a Poisson process
    rate["mode"] = cum_hist / T
    rate["mean"] = (cum_hist + 1) / T
    rate["stdev"] = (cum_hist + 1)**(1/2.) / T

    # compute 1-sigma peak density confidence interval for each snr-bin in 
    #     cumulative rate histogram
    sigma = 0.674490
    rate['lower_lim'] = np.zeros(len(rate['mode']))
    rate['upper_lim'] = np.zeros(len(rate['mode']))
    for i, R in enumerate(rate['mode']):
        # define limits to likely rate given peak and stdev
        max_range = rate['mode'][i]+5*rate['stdev'][i]
        if 5*rate['stdev'][i] < rate['mode'][i]:
            min_range = rate['mode'][i]-5*rate['stdev'][i]
        else:
            min_range = 0.0
        r = np.linspace(min_range, max_range, 1e4+1)

        #FIXME: using self-made gamma distribution because one in scipy
        #    version 0.7.2 is broken 
        cdf = stats.gamma.cdf(r, cum_hist[i]+1, scale=1./T)
        pdf = (cdf[1:]-cdf[:-1])/(r[1]-r[0])

        sort_idx = np.argsort(pdf)[::-1]
        norm_cumsum = np.cumsum(np.sort(pdf)[::-1])/np.sum(pdf)
        idx = np.argmin(np.abs(norm_cumsum - sigma))
        rate['lower_lim'][i] += R - r[np.min(sort_idx[:idx+1])]
        rate['upper_lim'][i] += r[np.max(sort_idx[:idx+1])] - R

    return rate

def calc_far( snrs, cn_snr, snr_bins, far_column, t_fgd, t_bkgd ):

    significance = []
    for snr in snrs:
        # which bin the coinc event falls into
        idx = bisect.bisect_right(snr_bins, snr) - 1
        if idx == (len(snr_bins) - 1):
            N = 0.0
        else:
            N = cn_snr[idx]

        # compute desired statistic
        if far_column == 'peak_far':
            significance.append( N/t_bkgd )
        elif far_column == 'mean_far':
            significance.append( (N+1)/t_bkgd )
        elif far_column == 'marginalized_fap':
            significance.append( 1-(1+t_fgd/t_bkgd)**(-N-1) )

    return significance

#
# =============================================================================
#
#                                SNR Histograms 
#
# =============================================================================
#

def sngl_snr_hist(
    connection,
    ifo,
    mchirp,
    eta,
    min_snr,
    snr_stat = None,
    sngls_width = None,
    usertag = "FULL_DATA",
    datatype = None,
    sngls_bins = None):
    """
    Creates a histogram of sngl_inspiral triggers and returns a list of counts
    and the associated snr bins.

    @connection: connection to a SQLite database with lsctables
    @ifo: the instrument one desires triggers from
    @mchirp: the chirp mass from the desired template
    @eta: the symmetric mass ratio from the desired template
    @min_snr: a lower threshold on the value of the snr_stat
    @sngls_width: the bin width for the histogram
    @usertag: the usertag for the triggers. The default is "FULL_DATA".
    @datatype: the datatype (all_data, slide, ...) if single-ifo triggers from
        coincident events is desired. The default is to collect all triggers.
    @sngls_bins: a list of bin edges for the snr-histogram
    """

    # create function for the desired snr statistic
    set_getsnr_function(connection, snr_stat)

    connection.create_function('end_time_w_ns', 2, end_time_w_ns)

    # set SQL statement parameters
    sql_params_dict = {
        "mchirp": mchirp, "eta": eta,
        "min_snr": min_snr, "ifo": ifo,
        "usertag": usertag}

    # SQLite query to get a list of (snr, gps-time) tuples for inspiral triggers
    sqlquery = """
    SELECT DISTINCT
        get_snr(snr, chisq, chisq_dof) as snr_stat,
        end_time_w_ns(end_time, end_time_ns)
    FROM sngl_inspiral
        JOIN process_params ON (
            process_params.process_id == sngl_inspiral.process_id) """
    # if a datatype is given, get only the inspiral triggers from coincs of that type
    if datatype:
        sqlquery += """
        JOIN coinc_event_map, experiment_map, experiment_summary ON (
            coinc_event_map.event_id == sngl_inspiral.event_id
            AND experiment_map.coinc_event_id == coinc_event_map.coinc_event_id
            AND experiment_summary.experiment_summ_id == experiment_map.experiment_summ_id) """
    sqlquery += """
    WHERE
        sngl_inspiral.ifo == :ifo 
        AND snr_stat >= :min_snr
        AND sngl_inspiral.mchirp == :mchirp
        AND sngl_inspiral.eta == :eta
        AND process_params.value == :usertag """
    if datatype:
        sqlquery += """
        AND experiment_summary.datatype == :type
        """
        sql_params_dict["type"] = datatype

    # get dq-veto segments
    xmldoc = dbtables.get_xml(connection)
    veto_segments = compute_dur.get_veto_segments(xmldoc, False)
    veto_segments = veto_segments[ veto_segments.keys()[0] ]

    snr_array = np.array([])
    # apply vetoes to the list of trigger times
    for snr, trig_time in connection.execute( sqlquery, sql_params_dict ):
        trig_segment = segments.segment(trig_time, trig_time)
        if not veto_segments[ifo].intersects_segment( trig_segment ):
            snr_array = np.append( snr_array, snr )

    if sngls_bins is None:
        sngls_bins = np.arange(min_snr, np.max(snr_array) + sngls_width, sngls_width)
    
    # make the binned snr histogram
    sngls_hist, _ = np.histogram(snr_array, bins=sngls_bins)

    return sngls_hist, sngls_bins


def all_sngl_snr_hist(
    connection,
    mchirp,
    eta,
    all_ifos,
    min_snr = 5.5,
    sngls_width = 0.01,
    snr_stat = None):
    """
    Creates a pair of dictionaries containing single-ifo snr histograms and
    their associated snr bins. The keys are the instruments filtered in the
    analysis.

    @connection: connection to a SQLite database with lsctables
    @min_snr: a lower threshold on the value of the snr_stat
    @sngls_width: the bin width for the histogram
    @mchirp: the chirp mass from the desired template
    @eta: the symmetric mass ratio from the desired template
    @all_ifos: a list containing the instruments used in the analysis 
    """

    sngl_ifo_hist = {}
    sngl_ifo_midbins = {}
    # loop over the analyzed ifos
    for ifo in all_ifos:
        sngl_ifo_hist[ifo], bins = sngl_snr_hist(
            connection,
            ifo,
            mchirp, eta,
            min_snr,
            sngls_width = sngls_width,
            snr_stat = snr_stat)
        # define the midpoint of each snr bin
        sngl_ifo_midbins[ifo] = 0.5*( bins[1:] + bins[:-1] )

    return sngl_ifo_hist, sngl_ifo_midbins


def coinc_snr_hist(
    connection,
    ifos,
    mchirp,
    eta,
    min_snr = 5.5,
    datatype = None,
    slide_id = None,
    combined_bins = None,
    snr_stat = None):
    """
    Creates a histogram of coinc_inspiral triggers and returns a list of counts
    in each of the snr bins.

    @connection: connection to a SQLite database with lsctables
    @sngls_threshold: a lower threshold on the value of the snr_stat
    @ifos: a tuple containing the ifos a coinc must come from 
    @mchirp: the chirp mass from the desired template
    @eta: the symmetric mass ratio from the desired template
    @datatype: the datatype (all_data, slide, ...) if single-ifo triggers from
        coincident events is desired. The default is to collect all triggers.
    @combined_bins: a list of bin edges for the snr-histogram
    """

    # create function for the desired snr statistic
    set_getsnr_function(connection, snr_stat)

    # set SQL statement parameters
    query_params_dict = {
        "ifo1": ifos[0], "ifo2": ifos[1],
        "type": datatype,
        "min_snr": min_snr,
        "mchirp": mchirp, "eta": eta
    }

    # get a list of the combined snrs from the coinc_inspiral table
    sqlquery = """
    SELECT
        coinc_inspiral.coinc_event_id,
        get_snr(si_ifo1.snr, si_ifo1.chisq, si_ifo1.chisq_dof),
        get_snr(si_ifo2.snr, si_ifo2.chisq, si_ifo2.chisq_dof)"""
    if len(ifos) > 2:
        sqlquery += """,
        get_snr(si_ifo3.snr, si_ifo3.chisq, si_ifo3.chisq_dof)"""
        query_params_dict["ifo3"] = ifos[2]

    sqlquery += """
    FROM coinc_inspiral
        JOIN sngl_inspiral AS si_ifo1, coinc_event_map AS cem_ifo1 ON (
            coinc_inspiral.coinc_event_id == cem_ifo1.coinc_event_id
            AND cem_ifo1.event_id == si_ifo1.event_id
            AND si_ifo1.ifo == :ifo1)
        JOIN sngl_inspiral AS si_ifo2, coinc_event_map AS cem_ifo2 ON (
            coinc_inspiral.coinc_event_id == cem_ifo2.coinc_event_id
            AND cem_ifo2.event_id == si_ifo2.event_id
            AND si_ifo2.ifo == :ifo2)"""
    if len(ifos) > 2:
        sqlquery += """
        JOIN sngl_inspiral AS si_ifo3, coinc_event_map AS cem_ifo3 ON (
            coinc_inspiral.coinc_event_id == cem_ifo3.coinc_event_id
            AND cem_ifo3.event_id == si_ifo3.event_id
            AND si_ifo3.ifo == :ifo3)"""

    sqlquery += """
        JOIN experiment_map AS expr_map, experiment_summary AS expr_summ ON (
            expr_map.coinc_event_id == coinc_inspiral.coinc_event_id
            AND expr_summ.experiment_summ_id == expr_map.experiment_summ_id)
    WHERE
        expr_summ.datatype == :type
        AND si_ifo1.mchirp == :mchirp
        AND si_ifo1.eta == :eta
        AND get_snr(si_ifo1.snr, si_ifo1.chisq, si_ifo1.chisq_dof) >= :min_snr
        AND get_snr(si_ifo2.snr, si_ifo2.chisq, si_ifo2.chisq_dof) >= :min_snr"""
    if len(ifos) > 2:
        sqlquery += """
        AND get_snr(si_ifo3.snr, si_ifo3.chisq, si_ifo3.chisq_dof) >= :min_snr"""
    if slide_id:
        query_params_dict['tsid'] = slide_id
        sqlquery += """
        AND expr_summ.time_slide_id == :tsid """

    # execute query
    events = connection.execute(sqlquery, query_params_dict).fetchall()

    coincs = {}
    coincs['snrs'] = np.array([ quadrature_sum(event[1:]) for event in events ])
    coincs['ids'] = [event[0] for event in events]

    # make histogram of coinc snr
    binned_snr, _ = np.histogram(coincs['snrs'], bins=combined_bins)
    coincs['snr_hist'] = map(float, binned_snr)

    return coincs


def all_possible_coincs(
    sngl_ifo_hist,
    sngl_ifo_midbins,
    combined_bins,
    ifos
):
    """
    Creates a histogram of all possible coincident events and returns a list of counts
    in each of the snr bins. This is made using the single-ifo snr histograms.

    @sngl_ifo_hist: a dictionary containing an snr histogram for each IFO
    @sngl_ifo_midbins: a dictionary the snr mid-bin values for sngl_ifo_hist
    @combined_bins: a list of bin edges for the combined-snr histogram
    @zerolag_coinc_hist: the snr histogram for zerolag coincident events
    @ifos: a list of analyzed instruments to generate a coinc
    """

    if len(ifos) > 3:
        raise ValueError, "Can only estimate FARs for doubles & triples"

    binWidth = combined_bins[1] - combined_bins[0]
    combined_counts = np.zeros(len(combined_bins)-1)

    # for computing doubles rate: N01_ij = n0_i*n1_j
    N01 = np.outer(sngl_ifo_hist[ifos[0]], sngl_ifo_hist[ifos[1]])
    len0, len1 = N01.shape
    for idx0, snr0 in enumerate(sngl_ifo_midbins[ifos[0]]):
        for idx1, snr1 in enumerate(sngl_ifo_midbins[ifos[1]]):
            if len(ifos) == 2:
                combined_snr = quadrature_sum( (snr0, snr1) )
                index =  int( np.floor((combined_snr - min(combined_bins))/binWidth) )
                combined_counts[index] += N01[idx0,idx1]
            else:
                # for computing triples rate: N012_ijk = n0_i*n1_j*n2_k
                N012 = np.outer(N01, sngl_ifo_hist[ifos[2]])
                N012 = N_012.reshape(len0, len1, N012.size/(len0*len1))
                for idx2, snr2 in enumerate(sngl_ifo_midbins[ifos[2]]):
                    combined_snr = quadrature_sum( (snr0, snr1, snr2) )
                    index =  int( np.floor((combined_snr - min(combined_bins))/binWidth) )
                    combined_counts[index] += N012[idx0,idx1,idx2]

    return combined_counts

#
# =============================================================================
#
#                                Time Functions 
#
# =============================================================================
#

def coinc_time(connection, datatype, inc_ifo_list, unit, tsid=None):
    params_dict = {'type': datatype}
    cursor = connection.cursor()

    sqlquery = """
    SELECT duration, instruments
    FROM experiment_summary
        JOIN experiment ON (
            experiment.experiment_id == experiment_summary.experiment_id)
    WHERE datatype = :type """
    if tsid:
        params_dict['tsid'] = tsid
        sqlquery += " AND time_slide_id = :tsid"
    cursor.execute( sqlquery, params_dict )

    # all elements of inclusive ifo_list must be found in ifos for T to be included
    livetime = np.sum([ T for T, ifos in cursor if set(ifos.split(',')).issuperset(set(inc_ifo_list)) ])

    return sqlutils.convert_duration(livetime, unit)


def get_singles_times( connection, verbose = False ):
    # get single-ifo filtered segments
    ifo_segments = compute_dur.get_single_ifo_segments(
        connection, program_name = "inspiral", usertag = "FULL_DATA")

    # get all veto segments
    xmldoc = dbtables.get_xml(connection)
    veto_segments = compute_dur.get_veto_segments(xmldoc, verbose)

    sngls_durations = {}

    for veto_def_name, veto_seg_dict in veto_segments.items():
        post_vetoes_ifosegs = ifo_segments - veto_seg_dict
        for ifo, segs_list in post_vetoes_ifosegs.items():
            sngls_durations[ifo] = float( abs(segs_list) )

    return sngls_durations


def get_coinc_window(connection, ifos):
    # determine the minimum time shift
    sqlquery = """
    SELECT MIN(ABS(offset))
    FROM time_slide
    WHERE offset != 0.0
    """
    shift = connection.execute( sqlquery ).fetchone()[0]

    # SQL query to get gps end-times for a given type of double
    sqlquery = """
    SELECT DISTINCT
        si_ifo1.end_time, si_ifo2.end_time,
        si_ifo1.end_time_ns, si_ifo2.end_time_ns
    FROM coinc_inspiral
        JOIN sngl_inspiral AS si_ifo1, coinc_event_map AS cem_ifo1 ON (
            coinc_inspiral.coinc_event_id == cem_ifo1.coinc_event_id 
            AND cem_ifo1.event_id == si_ifo1.event_id
            AND si_ifo1.ifo == ?)
        JOIN sngl_inspiral AS si_ifo2, coinc_event_map AS cem_ifo2 ON (
            coinc_inspiral.coinc_event_id == cem_ifo2.coinc_event_id 
            AND cem_ifo2.event_id == si_ifo2.event_id
            AND si_ifo2.ifo == ?)
    """
    # loop over pairs of instruments from the ifos list
    tau = {}
    for ifo_pair in iterutils.choices(ifos, 2):

        toa_diff = np.array([])
        # determine the difference in trigger end-times after sliding the times
        for coinc in connection.execute( sqlquery, tuple(ifos) ):
            dT_sec = (coinc[0]-coinc[1]) + (coinc[2]-coinc[3])*1e-9
            num_shifts = np.round(dT_sec/shift)
            toa_diff = np.append( toa_diff, dT_sec - num_shifts*shift )

        # the max-min of time differences defines the size of the coinc window
        tau[','.join(ifos)] = np.max(toa_diff) - np.min(toa_diff)

    return tau


def eff_bkgd_time(T_i, tau_ij, ifos, unit):
    # numerator is the product of the relevant single-ifo analyzed times
    numerator = np.prod([T for (ifo, T) in T_i.items() if ifo in ifos])

    # sorted list of the coincidence windows from smallest to largest
    taus = sorted(tau_ij.values())
    # denominator is a non-trivial combination of the coincident windows
    if len(ifos) == 2:
        denominator = taus[0]
    elif len(ifos) == 3:
        denominator = 0.5*tau[0]*(tau[2] + tau[1]) - \
                      0.25*( (tau[2]-tau[1])**2. + tau[0]**2. )
    else:
        raise ValueError, "Can only estimate background times for double & triples"

    return sqlutils.convert_duration(numerator/denominator, unit)

