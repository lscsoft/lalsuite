#!/usr/bin/env python2.2
"""
inspiral_pipeline.py - standalone inspiral pipeline driver script

$Id$

This script produced the necessary condor submit and dag files to run
the standalone inspiral code on LIGO data
"""

__author__ = 'Duncan Brown <duncan@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

import sys, os
import getopt, re
import tempfile
import pipeline, inspiral
import ConfigParser

# create the config parser object and read in the ini file
cp = ConfigParser.ConfigParser()
cp.read('pipe2.ini')

# create the DAG
dag = pipeline.CondorDAG('dag.log')
dag.set_dag_file('this.dag')

# create the jobs to use in the DAG
df_job = inspiral.DataFindJob(cp)
tmplt_job = inspiral.TmpltBankJob(cp)
insp_job = inspiral.InspiralJob(cp)
trig_job = inspiral.TrigToTmpltJob(cp)
inca_job = inspiral.IncaJob(cp)

# get the pad and chunk lengths from the values in the ini file
pad = int(cp.get('data', 'pad-data'))
n = int(cp.get('data', 'segment-length'))
s = int(cp.get('data', 'number-of-segments'))
r = int(cp.get('data', 'sample-rate'))
o = int(cp.get('inspiral', 'segment-overlap'))
length = ( n * s - ( s - 1 ) * o ) / r
overlap = o / r

# read the science segments from the input file and create the chunks
data = pipeline.ScienceData()
data.read(cp.get('input','segments'))
data.make_chunks(length,overlap,0)

vrbflg = 1
if vrbflg:
  print data
  for seg in data:
    print seg
    for chunk in seg:
      print chunk
    print seg.unused(), 'seconds remaining'

# get the order of the ifos to filter
ifo1 = cp.get('pipeline','ifo1')
ifo2 = cp.get('pipeline','ifo2')
ifo1_snr = cp.get('pipeline','ifo1-snr-threshold')
ifo2_snr = cp.get('pipeline','ifo2-snr-threshold')

# create all the LALdataFind jobs to run in sequence
prev_df1 = None
prev_df2 = None
first_df2 = None

insp_nodes = []

for seg in data:
  # find all the data
  df1 = inspiral.DataFindNode(df_job)

  df1.set_start(seg.start() - pad)
  df1.set_end(seg.end() + pad)
  df1.set_ifo(ifo1)
  if prev_df1: df1.add_parent(prev_df1)

  df2 = inspiral.DataFindNode(df_job)
  if not first_df2: first_df2 = df2

  df2.set_start(seg.start() - pad)
  df2.set_end(seg.end() + pad)
  df2.set_ifo(ifo2)
  if prev_df2: df2.add_parent(prev_df2)

  dag.add_node(df1)
  dag.add_node(df2)

  prev_df1 = df1
  prev_df1 = df2

  seg_insp_nodes = []

  for chunk in seg:
    bank = inspiral.TmpltBankNode(tmplt_job)
    bank.set_start(chunk.start())
    bank.set_end(chunk.end())
    bank.set_ifo(ifo1)
    bank.set_cache(df1.get_output())
    bank.add_parent(df1)
    dag.add_node(bank)

    insp1 = inspiral.InspiralNode(insp_job)
    insp1.set_start(chunk.start())
    insp1.set_end(chunk.end())
    insp1.set_ifo(ifo1)
    insp1.add_var('snr-threshold',ifo1_snr)
    insp1.set_cache(df1.get_output())
    insp1.set_bank(bank.get_output())
    insp1.add_parent(bank)

    trigbank = inspiral.TrigToTmpltNode(trig_job)
    trigbank.set_input(insp1.get_output())
    trigbank.set_output(ifo2 + '-TRIGBANK_' + ifo1 + '-' + str(chunk.start())
      + '-' + str(chunk.dur()) + '.xml')
    trigbank.add_parent(insp1)
    dag.add_node(trigbank)

    insp2 = inspiral.InspiralNode(insp_job)
    insp2.set_start(chunk.start())
    insp2.set_end(chunk.end())
    insp2.set_ifo(ifo2)
    insp2.add_var('snr-threshold',ifo2_snr)
    insp2.set_cache(df2.get_output())
    insp2.set_bank(trigbank.get_output())
    insp2.add_parent(df2)
    insp2.add_parent(trigbank)
    dag.add_node(insp2)

    seg_insp_nodes.append(tuple([insp1,insp2]))

  insp_nodes.append(seg_insp_nodes)
    
# now find coincidences between the two inspiral jobs
for i in range(len(data)):
  for j in range(len(data[i])):
    chunk = data[i][j]
    inca = inspiral.IncaNode(inca_job)
    inca.set_start(chunk.start())
    inca.set_end(chunk.end())
    inca.set_output( ifo1 + ifo2 + '-INCA-' + str(chunk.start())
      + '-' + str(chunk.dur()) + '.xml')
    
    try: 
      data[i][j-1]
      inca.add_parent(insp_nodes[i][j-1][0])
      inca.add_parent(insp_nodes[i][j-1][1])
      inca.add_input_a(insp_nodes[i][j-1][0].get_output())
      inca.add_input_b(insp_nodes[i][j-1][1].get_output())
    except IndexError:
      pass

    inca.add_parent(insp_nodes[i][j][0])
    inca.add_parent(insp_nodes[i][j][1])
    inca.add_input_a(insp_nodes[i][j][0].get_output())
    inca.add_input_b(insp_nodes[i][j][1].get_output())

    try:
      data[i][j+1]
      inca.add_parent(insp_nodes[i][j+1][0])
      inca.add_parent(insp_nodes[i][j+1][1])
      inca.add_input_a(insp_nodes[i][j+1][0].get_output())
      inca.add_input_b(insp_nodes[i][j+1][1].get_output())
    except IndexError:
      pass
      
    dag.add_node(inca)

# write out the DAG
dag.write_sub_files()
dag.write_dag()

sys.exit(0)
