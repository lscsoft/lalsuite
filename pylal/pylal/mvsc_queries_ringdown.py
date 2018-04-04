try:
	import sqlite3
except ImportError:
	# pre 2.5.x
	from pysqlite2 import dbapi2 as sqlite3
from glue.ligolw import dbtables
from glue.ligolw import table
from glue.ligolw import ilwd
from glue import segments
from pylal import SnglInspiralUtils
from pylal import db_thinca_rings
from pylal import git_version
from time import clock,time
from optparse import *
import glob
import sys
import random
import math

usage="""
this is a module for use in mvsc_get_doubles
"""

__author__ = "Kari Hodge <khodge@ligo.caltech.edu>, Paul T Baker <paul.baker@ligo.org>"

class CandidateEventQuery:
	# this is the list of parameters that will describe each event in the training and testing sets:
	parameters = "ds_sq delta_t df dQ gQQ gff gtt gQf gtf gtQ a_snr b_snr coinc_snr choppedl snr_sq snr_ratio a_eff_D b_eff_D eff_D_ratio delta_eff_D" #null_stat eff_coh_snr"
	# these are the sqlite queries used to extract these parameters (the dimensions to be considered in the multivariate statitical classification algorithm)
	select_count="""
		SELECT
			COUNT(coinc_ringdown.coinc_event_id)"""
	select_dimensions="""
		SELECT
			coinc_ringdown.coinc_event_id,
			snglA.*,
			snglB.*,
			insp_coinc_event.time_slide_id,
			calc_delta_t(snglA.ifo, snglA.start_time, snglA.start_time_ns, snglB.ifo, snglB.start_time, snglB.start_time_ns, insp_coinc_event.time_slide_id),
			abs(snglA.frequency - snglB.frequency),
			abs(snglA.Quality - snglB.Quality),
			gQQ(snglA.Quality, snglB.Quality),
			gff(snglA.frequency, snglB.frequency, snglA.Quality, snglB.Quality),
			gtt(snglA.frequency, snglB.frequency, snglA.Quality, snglB.Quality),
			gQf(snglA.frequency, snglB.frequency, snglA.Quality, snglB.Quality),
			gtf(snglA.Quality, snglB.Quality),
			gtQ(snglA.frequency, snglB.frequency, snglA.Quality, snglB.Quality),
			snglA.snr,
			snglB.snr,
			coinc_ringdown.snr,
			coinc_ringdown.choppedl_snr,
			coinc_ringdown.snr_sq,
			max(snglA.snr/snglB.snr,snglB.snr/snglA.snr),
			snglA.eff_dist,
			snglB.eff_dist,
			max(snglA.eff_dist/snglB.eff_dist,snglB.eff_dist/snglA.eff_dist),
			abs(snglA.eff_dist - snglB.eff_dist)"""
#			coinc_ringdown.null_stat,
#			coinc_ringdown.eff_coh_snr"""
	add_select_injections="""
		, coinc_ringdown.start_time+coinc_ringdown.start_time_ns*.000000001,
		process_params.value,
		sim_ringdown.distance"""
	add_from_injections="""
		FROM
			coinc_ringdown
			JOIN coinc_event_map AS mapA ON (mapA.coinc_event_id == coinc_ringdown.coinc_event_id)
			JOIN coinc_event_map AS mapB ON (mapB.coinc_event_id == coinc_ringdown.coinc_event_id)
			JOIN sngl_ringdown AS snglA ON (snglA.event_id == mapA.event_id)
			JOIN sngl_ringdown AS snglB ON (snglB.event_id == mapB.event_id)
			JOIN coinc_event_map AS mapC ON (mapC.event_id == coinc_ringdown.coinc_event_id)
			JOIN coinc_event_map AS mapD ON (mapD.coinc_event_id == mapC.coinc_event_id)
			JOIN sim_ringdown ON (sim_ringdown.simulation_id == mapD.event_id)
			JOIN coinc_event AS sim_coinc_event ON (sim_coinc_event.coinc_event_id == mapD.coinc_event_id)
			JOIN coinc_event AS insp_coinc_event ON (insp_coinc_event.coinc_event_id == mapA.coinc_event_id)
			JOIN coinc_definer ON (coinc_definer.coinc_def_id == sim_coinc_event.coinc_def_id)
			JOIN process_params ON (process_params.process_id == sim_ringdown.process_id)
		WHERE
			mapA.table_name == 'sngl_ringdown'
			AND mapB.table_name == 'sngl_ringdown'
			AND mapC.table_name == 'coinc_event'
			AND mapD.table_name == 'sim_ringdown'
			AND snglA.ifo == ?
			AND snglB.ifo == ?
			AND snglA.start_time > ?
			AND snglA.start_time < ?
			AND ( 
				(process_params.program == 'rinj' AND process_params.param == '--waveform' )
				OR (process_params.program == 'inspinj'AND process_params.param == '--d-distr') 
			)"""
	add_where_all="""
			AND coinc_definer.description == ?
		ORDER BY coinc_ringdown.start_time+coinc_ringdown.start_time_ns*.000000001"""
	add_where_exact="""
			AND coinc_definer.description == ?
		ORDER BY coinc_ringdown.start_time+coinc_ringdown.start_time_ns*.000000001"""
	add_select_fulldata="""
		, experiment_summary.datatype"""
	add_from_fulldata="""
		FROM
			coinc_ringdown
			JOIN coinc_event_map AS mapA ON (mapA.coinc_event_id == coinc_ringdown.coinc_event_id)
			JOIN coinc_event_map AS mapB ON (mapB.coinc_event_id == coinc_ringdown.coinc_event_id)
			JOIN sngl_ringdown AS snglA ON (snglA.event_id == mapA.event_id)
			JOIN sngl_ringdown AS snglB ON (snglB.event_id == mapB.event_id)
			JOIN coinc_event AS insp_coinc_event ON (mapA.coinc_event_id == insp_coinc_event.coinc_event_id)
			JOIN coinc_definer ON (coinc_definer.coinc_def_id == insp_coinc_event.coinc_def_id)
			JOIN experiment_map ON (experiment_map.coinc_event_id == coinc_ringdown.coinc_event_id)
			JOIN experiment_summary ON (experiment_summary.experiment_summ_id == experiment_map.experiment_summ_id)
		WHERE
			coinc_definer.search == 'ring'
			AND coinc_definer.search_coinc_type == 0
			AND mapA.table_name == 'sngl_ringdown'
			AND mapB.table_name == 'sngl_ringdown'
			AND snglA.ifo == ?
			AND snglB.ifo == ?
			AND snglA.start_time > ?
			AND snglA.start_time < ?"""

