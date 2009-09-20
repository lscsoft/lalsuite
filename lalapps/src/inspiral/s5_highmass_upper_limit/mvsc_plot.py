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

parser=OptionParser(usage='make plots', version='%prog')
parser.add_option("", "--injections", default="*INJ*.sqlite", help="glob of injection sqlite databases")
parser.add_option("", "--fulldata", default ="FULL_DATA*.sqlite", help="glob of full data sqlite databases (these should already include the timeslides)")
(opts,args)=parser.parse_args()

inj_files = glob.glob(opts.injections)
fulldata_files = glob.glob(opts.fulldata)

timeslide_likelihood = []
timeslide_snr = []
zerolag_likelihood = []
zerolag_snr = []
for filename in fulldata_files:
  local_disk = None #"/tmp"
  working_filename = dbtables.get_connection_filename(filename, tmp_path = local_disk, verbose = True)
  connection = sqlite3.connect(working_filename)
  dbtables.DBTable_set_connection(connection)
  xmldoc = dbtables.get_xml(connection)
  cursor = connection.cursor()
  for likelihood, snr, is_background in connection.cursor().execute("""
  SELECT 
    insp_coinc_event.likelihood, 
    coinc_inspiral.snr,
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
    coinc_inspiral 
    JOIN coinc_event AS insp_coinc_event ON (insp_coinc_event.coinc_event_id == coinc_inspiral.coinc_event_id)
  """):
    if is_background:
      timeslide_likelihood.append(likelihood)
      timeslide_snr.append(snr)
    else:
      zerolag_likelihood.append(likelihood)
      zerolag_snr.append(snr)
  dbtables.put_connection_filename(filename, working_filename, verbose = True)
#print len(timeslide_likelihood)
#print len(zerolag_likelihood)

injection_likelihood = []
injection_snr = []
for filename in inj_files:
  local_disk = None #"/tmp"
  working_filename = dbtables.get_connection_filename(filename, tmp_path = local_disk, verbose = True)
  connection = sqlite3.connect(working_filename)
  dbtables.DBTable_set_connection(connection)
  xmldoc = dbtables.get_xml(connection)
  cursor = connection.cursor()
  for likelihood,snr in connection.cursor().execute("""
  SELECT
    insp_coinc_event.likelihood, 
    coinc_inspiral.snr
  FROM
    coinc_inspiral 
    JOIN coinc_event_map AS mapC ON (mapC.event_id == coinc_inspiral.coinc_event_id)
    JOIN coinc_event_map AS mapD ON (mapD.coinc_event_id == mapC.coinc_event_id)
    JOIN sim_inspiral ON (sim_inspiral.simulation_id == mapD.event_id)
    JOIN coinc_event AS sim_coinc_event ON (sim_coinc_event.coinc_event_id == mapD.coinc_event_id)
    JOIN coinc_event AS insp_coinc_event ON (insp_coinc_event.coinc_event_id == coinc_inspiral.coinc_event_id)
  WHERE
    sim_coinc_event.coinc_def_id == 'coinc_definer:coinc_def_id:2'
    AND insp_coinc_event.coinc_def_id == 'coinc_definer:coinc_def_id:0'
    AND mapC.table_name == 'coinc_event'
    AND mapD.table_name == 'sim_inspiral'
  """):
    injection_likelihood.append(likelihood)
    injection_snr.append(snr)
  dbtables.put_connection_filename(filename, working_filename, verbose = True)
#print len(injection_likelihood)
#print timeslide_likelihood
#print injection_likelihood

# map all 0 likelihoods to lowest non-zero likelihood and all 'inf' likelihoods to the highest non-infinity likelihood
all_likelihood = timeslide_likelihood + zerolag_likelihood + injection_likelihood
likelihood_set = set(all_likelihood)
likelihood_list = list(likelihood_set)
likelihood_list.sort()
#print likelihood_list
min_likelihood = likelihood_list[1]
max_likelihood = likelihood_list[-2]
#print min_likelihood
#print max_likelihood
for list in timeslide_likelihood, zerolag_likelihood, injection_likelihood:
  #print 'list'
  #print list
  for i in range(len(list)):
    if list[i] == 0:
      list[i] = min_likelihood
    if list[i] == float('inf'):
      list[i] = max_likelihood
  #print list
# find the "threshold," which is the likelihood of the 100th loudest timeslide
timeslide_likelihood.sort()
pylab.figure(0)
pylab.loglog(injection_likelihood,injection_snr,'r*',label='injections')
pylab.hold(1)
pylab.loglog(timeslide_likelihood,timeslide_snr,'.',label='time slides')
pylab.axvline(x=timeslide_likelihood[-100], label='100th loudest time slide')
pylab.xlabel('likelihood')
pylab.ylabel('combined effective SNR')
pylab.legend(loc='lower right') 
pylab.hold(0)
pylab.savefig('myfigure.png')
