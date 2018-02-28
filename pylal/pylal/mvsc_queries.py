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

__author__ = "Kari Hodge <khodge@ligo.caltech.edu>"

class CandidateEventQuery:
	# this is the list of parameters that will describe each event in the training and testing sets:
	parameters = "ethinca delta_t ab_dmchirp_rel ab_deta_rel a_snr b_snr a_chisq_red b_chisq_red a_effective_snr b_effective_snr a_rsq_veto_duration b_rsq_veto_duration a_cont_chisq_red b_cont_chisq_red coinc_inspiral_snr sngl_gps_time_a sngl_gps_time_b"
	# these are the sqlite queries used to extract these parameters (the dimensions to be considered in the multivariate statitical classification algorithm) 
	select_count="""
		SELECT
			COUNT(coinc_inspiral.coinc_event_id)"""
	select_dimensions="""
		SELECT
			coinc_inspiral.coinc_event_id,
			snglA.*,
			snglB.*,
			insp_coinc_event.time_slide_id,
			calc_delta_t(snglA.ifo, snglA.end_time, snglA.end_time_ns, snglB.ifo, snglB.end_time, snglB.end_time_ns, insp_coinc_event.time_slide_id),
			abs(2*(snglA.mchirp - snglB.mchirp)/(snglA.mchirp+snglB.mchirp)),
			abs(2*(snglA.eta - snglB.eta)/(snglA.eta+snglB.eta)),
			snglA.snr,
			snglB.snr,
			snglA.chisq/(2*snglA.chisq_dof-2),
			snglB.chisq/(2*snglB.chisq_dof-2),
			calc_effective_snr(snglA.snr, snglA.chisq, snglA.chisq_dof),
			calc_effective_snr(snglB.snr, snglB.chisq, snglB.chisq_dof),
			snglA.rsqveto_duration,
			snglB.rsqveto_duration,
			CASE snglA.cont_chisq_dof
				WHEN 0.0 THEN 1.0
				ELSE snglA.cont_chisq/snglA.cont_chisq_dof END,
			CASE snglB.cont_chisq_dof
				WHEN 0.0 THEN 1.0
				ELSE snglB.cont_chisq/snglB.cont_chisq_dof END,
			coinc_inspiral.snr,
			snglA.end_time+snglA.end_time_ns*.000000001,
			snglB.end_time+snglB.end_time_ns*.000000001"""
	add_select_injections="""
		, coinc_inspiral.end_time+coinc_inspiral.end_time_ns*.000000001,
		process_params.value,
		sim_inspiral.distance"""
	add_select_fulldata="""
		, experiment_summary.datatype"""
	add_from_injections="""
		FROM
			coinc_inspiral
			JOIN coinc_event_map AS mapA ON (mapA.coinc_event_id == coinc_inspiral.coinc_event_id)
			JOIN coinc_event_map AS mapB ON (mapB.coinc_event_id == coinc_inspiral.coinc_event_id)
			JOIN sngl_inspiral AS snglA ON (snglA.event_id == mapA.event_id)
			JOIN sngl_inspiral AS snglB ON (snglB.event_id == mapB.event_id)
			JOIN coinc_event_map AS mapC ON (mapC.event_id == coinc_inspiral.coinc_event_id)
			JOIN coinc_event_map AS mapD ON (mapD.coinc_event_id == mapC.coinc_event_id)
			JOIN sim_inspiral ON (sim_inspiral.simulation_id == mapD.event_id)
			JOIN coinc_event AS sim_coinc_event ON (sim_coinc_event.coinc_event_id == mapD.coinc_event_id)
			JOIN coinc_event AS insp_coinc_event ON (insp_coinc_event.coinc_event_id == mapA.coinc_event_id)
			JOIN coinc_definer ON (coinc_definer.coinc_def_id == sim_coinc_event.coinc_def_id)
			JOIN process_params ON (process_params.process_id == sim_inspiral.process_id)
		WHERE
			mapA.table_name == 'sngl_inspiral'
			AND mapB.table_name == 'sngl_inspiral'
			AND mapC.table_name == 'coinc_event'
			AND mapD.table_name == 'sim_inspiral'
			AND snglA.ifo == ?
			AND snglB.ifo == ?
			AND process_params.program == 'inspinj' AND process_params.param == '--d-distr'"""
	add_where_all="""
			AND coinc_definer.coinc_def_id == "coinc_definer:coinc_def_id:12"
		ORDER BY coinc_inspiral.end_time+coinc_inspiral.end_time_ns*.000000001"""
	add_where_exact="""
			AND coinc_definer.description == ?
		ORDER BY coinc_inspiral.end_time+coinc_inspiral.end_time_ns*.000000001"""
	add_from_fulldata="""
		FROM
			coinc_inspiral
			JOIN coinc_event_map AS mapA ON (mapA.coinc_event_id == coinc_inspiral.coinc_event_id)
			JOIN coinc_event_map AS mapB ON (mapB.coinc_event_id == coinc_inspiral.coinc_event_id)
			JOIN sngl_inspiral AS snglA ON (snglA.event_id == mapA.event_id)
			JOIN sngl_inspiral AS snglB ON (snglB.event_id == mapB.event_id)
			JOIN coinc_event AS insp_coinc_event ON (mapA.coinc_event_id == insp_coinc_event.coinc_event_id)
			JOIN coinc_definer ON (coinc_definer.coinc_def_id == insp_coinc_event.coinc_def_id)
			JOIN experiment_map ON (experiment_map.coinc_event_id == coinc_inspiral.coinc_event_id)
			JOIN experiment_summary ON (experiment_summary.experiment_summ_id == experiment_map.experiment_summ_id)
		WHERE
			coinc_definer.search == 'inspiral'
			AND coinc_definer.search_coinc_type == 0
			AND mapA.table_name == 'sngl_inspiral'
			AND mapB.table_name == 'sngl_inspiral'
			AND snglA.ifo == ?
			AND snglB.ifo == ?
			ORDER BY coinc_inspiral.end_time+coinc_inspiral.end_time_ns*.000000001"""
