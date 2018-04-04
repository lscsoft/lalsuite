#!/usr/bin/python
try:
  import sqlite3
except ImportError:
  # pre 2.5.x
  from pysqlite2 import dbapi2 as sqlite3
from glue.ligolw import dbtables
from optparse import *

import glob
import matplotlib
matplotlib.use('Agg')
import pylab
import numpy
import bisect
from pylal import db_rinca_rings
from glue import segments

parser=OptionParser(usage="""
%prog [sqlite database file(s)]
in S6, fulldata databases include timeslides, injections, and zerolag
so in general, you will only be providing one file, for example:
~/lalsuite/pylal/bin/mvsc_plot.py H1L1-FULL_DATA_CAT_2_VETO_CLUSTERED_CBC_RESULTS-951868815-1209600.sqlite
"""
, version='%prog')
(opts,files)=parser.parse_args()

timeslide_likelihood = []
timeslide_snr = []
zerolag_likelihood = []
zerolag_snr = []
injection_likelihood = []
injection_snr = []
for filename in files:
  local_disk = None #"/tmp"
  working_filename = dbtables.get_connection_filename(filename, tmp_path = local_disk, verbose = True)
  connection = sqlite3.connect(working_filename)
  dbtables.DBTable_set_connection(connection)
  xmldoc = dbtables.get_xml(connection)
  cursor = connection.cursor()
  for likelihood, snr, is_background in connection.cursor().execute("""
  SELECT
    insp_coinc_event.likelihood,
    coinc_ringdown.snr,
    EXISTS (
      SELECT
        *
      FROM
        time_slide
      WHERE
       time_slide.time_slide_id == insp_coinc_event.time_slide_id
       AND time_slide.offset != 0
     )
  FROM
    coinc_ringdown
    JOIN coinc_event AS insp_coinc_event ON (insp_coinc_event.coinc_event_id == coinc_ringdown.coinc_event_id)
  """):
    if is_background:
      timeslide_likelihood.append(likelihood)
      timeslide_snr.append(snr)
    else:
      zerolag_likelihood.append(likelihood)
      zerolag_snr.append(snr)
  for likelihood,snr in connection.cursor().execute("""
  SELECT
    insp_coinc_event.likelihood,
    coinc_ringdown.snr
  FROM
    coinc_ringdown
    JOIN coinc_event_map AS mapC ON (mapC.event_id == coinc_ringdown.coinc_event_id)
    JOIN coinc_event_map AS mapD ON (mapD.coinc_event_id == mapC.coinc_event_id)
    JOIN sim_inspiral ON (sim_inspiral.simulation_id == mapD.event_id)
    JOIN coinc_event AS sim_coinc_event ON (sim_coinc_event.coinc_event_id == mapD.coinc_event_id)
    JOIN coinc_event AS insp_coinc_event ON (insp_coinc_event.coinc_event_id == coinc_ringdown.coinc_event_id)
    JOIN coinc_definer AS insp_coinc_definer ON (insp_coinc_definer.coinc_def_id == insp_coinc_event.coinc_def_id)
    JOIN coinc_definer AS sim_coinc_definer ON (sim_coinc_definer.coinc_def_id == sim_coinc_event.coinc_def_id)
  WHERE
    insp_coinc_definer.search == 'ring'
    AND sim_coinc_definer.search == 'ring'
    AND insp_coinc_definer.search_coinc_type == 0
    AND sim_coinc_definer.search_coinc_type == 2
    AND mapC.table_name == 'coinc_event'
    AND mapD.table_name == 'sim_inspiral'
  """):
    injection_likelihood.append(likelihood)
    injection_snr.append(snr)
  #livetimes = db_rinca_rings.get_rinca_livetimes(db_rinca_rings.get_rinca_rings_by_available_instruments(connection,'rinca'), db_rinca_rings.get_veto_segments(connection, 'vetoes'), db_rinca_rings.get_background_offset_vectors(connection), verbose=True)
  # above is the correct calculation of livetimes, but it takes a while. You can comment it out and uncomment the line below for an approximate calculation of livetime.
  livetimes = db_rinca_rings.get_rinca_livetimes(db_rinca_rings.get_rinca_rings_by_available_instruments(connection,'rinca'), segments.segmentlistdict(), db_rinca_rings.get_background_offset_vectors(connection), verbose=True)
  dbtables.put_connection_filename(filename, working_filename, verbose = True)
print "number of timeslides:", len(timeslide_likelihood)
print "number of zerolags:", len(zerolag_likelihood)
print "number of injections:", len(injection_likelihood)

all_likelihoods = numpy.array(timeslide_likelihood + zerolag_likelihood + injection_likelihood, dtype=float)
timeslide_likelihoods = numpy.array(timeslide_likelihood)
zerolag_likelihoods = numpy.array(zerolag_likelihood)
injection_likelihoods = numpy.array(injection_likelihood)
timeslide_snrs = numpy.array(timeslide_snr)
zerolag_snrs = numpy.array(zerolag_snr)
injection_snrs = numpy.array(injection_snr)
# set all zero likelihoods to the next lowest calculated likelihood, and all infinity likelihoods to the next highest, so plotting can work
plottable_likelihoods = all_likelihoods[~numpy.isinf(all_likelihoods) & (all_likelihoods > 0)]
min_likelihood = plottable_likelihoods.min()
max_likelihood = plottable_likelihoods.max()
timeslide_likelihoods.clip(min_likelihood, max_likelihood, out=timeslide_likelihoods)
zerolag_likelihoods.clip(min_likelihood, max_likelihood, out=zerolag_likelihoods)
injection_likelihoods.clip(min_likelihood, max_likelihood, out=injection_likelihoods)
# to find the 'threshold,' which here is the likelihood or combined effective snr of the 100th loudest timeslide
sorted_timeslide_likelihoods = numpy.sort(timeslide_likelihoods)
sorted_timeslide_snrs = numpy.sort(timeslide_snrs)
pylab.figure(0)
pylab.loglog(injection_likelihoods,injection_snrs,'rx',label='injections')
pylab.hold(1)
pylab.loglog(timeslide_likelihoods,timeslide_snrs,'k.',label='timeslides')
pylab.axvline(x=sorted_timeslide_likelihoods[-100], color='g', label='100th loudest timeslide, by MVSCL')
pylab.axhline(y=sorted_timeslide_snrs[-100], label='100th loudest timeslide, by SNR')
pylab.xlabel('MVSC likelihood')
pylab.ylabel('combined effective SNR')
pylab.legend(loc='lower right')
pylab.hold(0)
pylab.savefig('MVSC_likelihood_scatterplot.png')

t=0
for comb in livetimes:
  t += sum(livetimes[comb])
print "total livetime:", t

bkg_likelihoods = timeslide_likelihood[:]
bkg_likelihoods.sort()
bkg_snrs = timeslide_snr[:]
bkg_snrs.sort()

injection_likelihood_ifar = []
timeslide_likelihood_ifar = []
zerolag_likelihood_ifar = []
injection_snr_ifar = []
timeslide_snr_ifar = []
zerolag_snr_ifar = []

def get_ifar(statistic_list,bkg_statistics,ifar_list):
  for i in range(len(statistic_list)):
    n = len(bkg_likelihoods) - bisect.bisect_left(bkg_statistics,statistic_list[i])
    if n != 0:
      ifar_list.append(t/float(n))
    else:
      ifar_list.append(float('inf'))

get_ifar(injection_likelihood,bkg_likelihoods,injection_likelihood_ifar)
get_ifar(timeslide_likelihood,bkg_likelihoods,timeslide_likelihood_ifar)
get_ifar(zerolag_likelihood,bkg_likelihoods,zerolag_likelihood_ifar)
get_ifar(injection_snr,bkg_snrs,injection_snr_ifar)
get_ifar(timeslide_snr,bkg_snrs,timeslide_snr_ifar)
get_ifar(zerolag_snr,bkg_snrs,zerolag_snr_ifar)

# note, infs are not plotted!
pylab.figure(1)
pylab.loglog(injection_likelihood_ifar,injection_snr_ifar,'rx',label='injections')
pylab.hold(1)
pylab.loglog(timeslide_likelihood_ifar,timeslide_snr_ifar,'k.',label='timeslides')
pylab.hold(0)
pylab.xlabel('IFAR from MVSC Likelihood')
pylab.ylabel('IFAR from combined effective SNR')
pylab.savefig('IFAR_IFAR_scatterplot.png')
ax = pylab.figure(0).gca()

pylab.figure(2)
pylab.hexbin(injection_likelihood_ifar,injection_snr_ifar,gridsize=50,xscale='log',yscale='log',bins='log',cmap=pylab.cm.jet,alpha=1)
pylab.gca = ax
pylab.colorbar()
pylab.xlabel('IFAR from MVSC Likelihood')
pylab.ylabel('IFAR from combined effective SNR')
pylab.title('log10 density of injections')
pylab.savefig('IFAR_densityfig_injection.png')

pylab.figure(3)
pylab.hexbin(timeslide_likelihood_ifar,timeslide_snr_ifar,gridsize=50,xscale='log',yscale='log',bins='log',cmap=pylab.cm.jet,alpha=1)
pylab.gca = ax
pylab.colorbar()
pylab.xlabel('IFAR from MVSC Likelihood')
pylab.ylabel('IFAR from combined effective SNR')
pylab.title('log10 density of timeslides')
pylab.savefig('IFAR_densityfig_timeslide.png')
