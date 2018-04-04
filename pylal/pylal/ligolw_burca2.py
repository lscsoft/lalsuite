# Copyright (C) 2007--2014  Kipp Cannon
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
#                                   Preamble
#
# =============================================================================
#


import sys
import traceback


from glue.ligolw import lsctables
from glue.text_progress_bar import ProgressBar
from pylal import git_version


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                              Library Interface
#
# =============================================================================
#


#
# Core routine
#


def assign_likelihood_ratios(connection, coinc_def_id, offset_vectors, vetoseglists, events_func, veto_func, ln_likelihood_ratio_func, likelihood_params_func, verbose = False, params_func_extra_args = ()):
	"""
	Assigns likelihood ratio values to coincidences.
	"""
	#
	# Convert offset vector keys to strings so that we can use the
	# dictionary inside an SQL query (they might be
	# glue.ligolw.ilwd_char objects)
	#

	offset_vectors = dict((unicode(time_slide_id), offset_vector) for time_slide_id, offset_vector in offset_vectors.items())

	#
	# Create a cursor object for events_func() to reuse
	#

	cursor = connection.cursor()

	#
	# Construct the in-SQL likelihood ratio function.  Rely on Python's
	# closure mechanism to retain all local variables at the time of
	# this function's creation for use inside the function.
	#

	def ln_likelihood_ratio(coinc_event_id, time_slide_id):
		try:
			params = likelihood_params_func([event for event in events_func(cursor, coinc_event_id) if veto_func(event, vetoseglists)], offset_vectors[time_slide_id], *params_func_extra_args)
			return ln_likelihood_ratio_func(params) if params is not None else None
		except:
			traceback.print_exc()
			raise

	connection.create_function("ln_likelihood_ratio", 2, ln_likelihood_ratio)

	#
	# Iterate over all coincs, assigning likelihood ratios.
	#

	if verbose:
		print >>sys.stderr, "computing likelihood ratios ..."

	connection.cursor().execute("""
UPDATE
	coinc_event
SET
	likelihood = ln_likelihood_ratio(coinc_event_id, time_slide_id)
WHERE
	coinc_def_id == ?
	""", (unicode(coinc_def_id),))

	#
	# Done
	#

	connection.commit()
	cursor.close()


def assign_likelihood_ratios_xml(xmldoc, coinc_def_id, offset_vectors, vetoseglists, events_func, veto_func, ln_likelihood_ratio_func, likelihood_params_func, verbose = False, params_func_extra_args = ()):
	"""
	Assigns likelihood ratio values to coincidences (XML version).
	"""
	#
	# Iterate over all coincs, assigning likelihood ratios.
	#

	coinc_event_table = lsctables.CoincTable.get_table(xmldoc)

	if verbose:
		progressbar = ProgressBar("computing likelihood ratios", max = len(coinc_event_table))
	else:
		progressbar = None

	for coinc_event in coinc_event_table:
		if progressbar is not None:
			progressbar.increment()
		if coinc_event.coinc_def_id != coinc_def_id:
			continue
		params = likelihood_params_func([event for event in events_func(None, coinc_event.coinc_event_id) if veto_func(event, vetoseglists)], offset_vectors[coinc_event.time_slide_id], *params_func_extra_args)
		coinc_event.likelihood = ln_likelihood_ratio_func(params) if params is not None else None

	del progressbar

	#
	# Done
	#

	return


#
# Burst-specific interface
#


def sngl_burst_events_func(cursor, coinc_event_id, row_from_cols):
	return map(row_from_cols, cursor.execute("""
SELECT
	sngl_burst.*
FROM
	sngl_burst
	JOIN coinc_event_map ON (
		coinc_event_map.table_name == 'sngl_burst'
		AND coinc_event_map.event_id == sngl_burst.event_id
	)
WHERE
	coinc_event_map.coinc_event_id == ?
	""", (coinc_event_id,)))


def sngl_burst_veto_func(event, vetoseglists):
	# return True if event should be *retained*
	return event.ifo not in vetoseglists or event.peak not in vetoseglists[event.ifo]


def ligolw_burca2(database, ln_likelihood_ratio, params_func, verbose = False, params_func_extra_args = ()):
	"""
	Assigns likelihood ratio values to excess power coincidences.
	database is pylal.SnglBurstUtils.CoincDatabase instance, and
	ln_likelihood_ratio is a LnLikelihoodRatio class instance.
	"""
	#
	# Run core function
	#

	assign_likelihood_ratios(
		connection = database.connection,
		coinc_def_id = database.bb_definer_id,
		offset_vectors = database.time_slide_table.as_dict(),
		vetoseglists = database.vetoseglists,
		events_func = lambda cursor, coinc_event_id: sngl_burst_events_func(cursor, coinc_event_id, database.sngl_burst_table.row_from_cols),
		veto_func = sngl_burst_veto_func,
		ln_likelihood_ratio_func = ln_likelihood_ratio,
		likelihood_params_func = params_func,
		verbose = verbose,
		params_func_extra_args = params_func_extra_args
	)

	#
	# Done
	#

	return database
