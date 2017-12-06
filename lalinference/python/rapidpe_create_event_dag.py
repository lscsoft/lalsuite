# Copyright (C) 2012 Chris Pankow
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

"""
Creates a dag workflow to perform extrinsic marginalization calculation.
"""

import sys
import os
from argparse import ArgumentParser

import numpy as np

import lal

import glue.lal
from glue.ligolw import utils, ligolw, lsctables, table
lsctables.use_in(ligolw.LIGOLWContentHandler)
from glue.ligolw.utils import process

from glue.ligolw.utils import process
from glue import pipeline

from lalinference.rapid_pe import common_cl, dagutils
from lalinference.rapid_pe.common_cl import JOB_PRIORITIES, MAXJOBS

__author__ = "Evan Ochsner <evano@gravity.phys.uwm.edu>, Chris Pankow <pankow@gravity.phys.uwm.edu>, R. O'Shaughnessy"

#
# Option parsing
#

argp = ArgumentParser()
# Options needed by this program only.

argp.add_argument("-T", "--template-bank-xml", help="Input template bank as a sim_inspiral or sngl_inspiral table. Required.")
argp.add_argument("-D", "--working-directory", default="./", help="Directory in which to stage DAG components.")
argp.add_argument("-l", "--log-directory", default="./", help="Directory in which to place condor logs.")
argp.add_argument("-W", "--web-output", default="./", help="Directory to place web accessible plots and webpages.")
argp.add_argument("-O", "--output-name", default="marginalize_extrinsic_parameters", help="Filename (without extension) to write DAG to.")
argp.add_argument("--n-copies", default=1, help="Number of copies of each integrator instance to run per mass point. Default is one.")
argp.add_argument("--write-script", action="store_true", help="In addition to the DAG, write a script to this filename to execute the workflow.")
argp.add_argument("--write-eff-lambda", action="store_true", help="Use psi0 column of template bank XML as effective lambda point to calculate in DAG.")
argp.add_argument("--write-deff-lambda", action="store_true", help="Use psi3 column of template bank XML as delta effective lambda point to calculate in DAG.")
argp.add_argument("--condor-command", action="append", help="Append these condor commands to the submit files. Useful for account group information.")

for cat, val in MAXJOBS.iteritems():
    optname = "--maxjobs-%s" % cat.lower().replace("_", "-")
    argp.add_argument(optname, type=int, default=MAXJOBS[cat], help="Set MAXJOBS in DAGs for category %s. Default is %s" % (cat, str(val)))

# Options transferred to ILE
common_cl.add_datasource_params(argp)
common_cl.add_integration_params(argp)
common_cl.add_output_params(argp)
common_cl.add_intrinsic_params(argp)
common_cl.add_pinnable_params(argp)

opts = argp.parse_args()

if not opts.template_bank_xml:
    exit("Option --template-bank-xml is required.")

condor_commands = None
if opts.condor_command is not None:
    condor_commands = dict([c.split("=") for c in opts.condor_command])

#
# Get trigger information from coinc xml file
#
# FIXME: CML should package this for us

# Get end time from coinc inspiral table or command line
xmldoc = None
if opts.coinc_xml is not None:
    xmldoc = utils.load_filename(opts.coinc_xml, contenthandler=ligolw.LIGOLWContentHandler)
    coinc_table = table.get_table(xmldoc, lsctables.CoincInspiralTable.tableName)
    assert len(coinc_table) == 1
    coinc_row = coinc_table[0]
    event_time = coinc_row.get_end()
    print "Coinc XML loaded, event time: %s" % str(coinc_row.get_end())
elif opts.event_time is not None:
    event_time = glue.lal.LIGOTimeGPS(opts.event_time)
    print "Event time from command line: %s" % str(event_time)
else:
    raise ValueError("Either --coinc-xml or --event-time must be provided to parse event time.")

xmldoc, tmplt_bnk = utils.load_filename(opts.template_bank_xml, contenthandler=ligolw.LIGOLWContentHandler), None
try:
    tmplt_bnk = lsctables.SimInspiralTable.get_table(xmldoc)
except ValueError:
    print >>sys.stderr, "Exactly one sim_inspiral table was not found in %s, trying sngl_inspiral" % opts.template_bank_xml

if tmplt_bnk is None:
    tmplt_bnk = lsctables.SnglInspiralTable.get_table(xmldoc)

#
# Post processing options
#
# FIXME: Remove these entirely
use_ile_postproc = False
use_bayespe_postproc = False

# initialize the analysis subdag
dag = pipeline.CondorDAG(log=os.getcwd())

if opts.maxjobs_ile is not None:
    dag.add_maxjobs_category("ILE", opts.maxjobs_ile)

# This is a subdag used for all our plotting and postproc so they don't block
# completion of an individual event's ILEs
ppdag = pipeline.CondorDAG(log=os.getcwd())
ppdag.add_maxjobs_category("SQL", MAXJOBS["SQL"])
ppdag.add_maxjobs_category("PLOT", MAXJOBS["PLOT"])

if not os.path.exists(opts.log_directory):
    os.makedirs(opts.log_directory) # Make a directory to hold log files of jobs

# All the intrinsic parameters we're gridding in
intr_prms = set(("mass1", "mass2"))
for p in ("spin1z", "spin2z"): # FIXME: Add all
    if hasattr(tmplt_bnk[0], p):
        intr_prms.add(p)

# These have explicit options because they map to non-standard columns and I
# want the user to explicity use these columns if they've written them
if opts.write_eff_lambda:
    intr_prms.add("eff_lambda")
if opts.write_deff_lambda:
    intr_prms.add("deff_lambda")

ile_job_type, ile_sub_name = dagutils.write_integrate_likelihood_extrinsic_sub(
        tag='integrate',
        condor_commands=condor_commands,
        intr_prms=intr_prms,
        log_dir=opts.log_directory,
        cache_file=opts.cache_file,
        channel_name=opts.channel_name,
        psd_file=opts.psd_file,
        coinc_xml=opts.coinc_xml,
        reference_freq=opts.reference_freq,
        fmax=opts.fmax,
        fmin_template=opts.fmin_template,
		approximant=opts.approximant,
		amp_order=opts.amp_order,
		l_max=opts.l_max,
        event_time=event_time,
        time_marginalization=opts.time_marginalization,
        save_samples=opts.save_samples,
        output_file=opts.output_file,
        n_eff=opts.n_eff,
        n_max=opts.n_max,
        ncopies=opts.n_copies,
        save_deltalnL=opts.save_deltalnL,
        save_P=opts.save_P,
        n_chunk=opts.n_chunk,
        convergence_tests_on=opts.convergence_tests_on,
        adapt_floor_level=opts.adapt_floor_level,
        adapt_weight_exponent=opts.adapt_weight_exponent,
        skymap_file=opts.skymap_file,
        distance_maximum=opts.distance_maximum
        )
ile_job_type.write_sub_file()

if use_bayespe_postproc:
    if not os.path.exists(opts.web_output):
        os.makedirs(opts.web_output)
    bpp_plot_job_type, bpp_plot_job_name = dagutils.write_bayes_pe_postproc_sub(tag="bayes_pp_plot", log_dir=opts.log_directory, web_dir=opts.web_output)
    bpp_plot_job_type.write_sub_file()
    bpp_plot_node = pipeline.CondorDAGNode(bpp_plot_job_type)
    bpp_plot_node.set_category("PLOT")
    bpp_plot_node.set_pre_script(dagutils.which("bayes_pe_preprocess"))
    ppdag.add_node(bpp_plot_node)

#
# Make the posterior plot here since we need to make it the child of every sql
# node in the DAG
#
if use_ile_postproc:
    pos_plot_job_type, pos_plot_job_name = dagutils.write_posterior_plot_sub(tag="pos_plot", log_dir=opts.log_directory)
    pos_plot_job_type.write_sub_file()
    pos_plot_node = pipeline.CondorDAGNode(pos_plot_job_type)
    pos_plot_node.set_pre_script(dagutils.which("coalesce.sh"))
    pos_plot_node.set_category("PLOT")
    pos_plot_node.set_priority(JOB_PRIORITIES["PLOT"])
    ppdag.add_node(pos_plot_node)

sql_job_type, sql_job_name = dagutils.write_result_coalescence_sub(tag="coalesce", log_dir=opts.log_directory)
sql_job_type.write_sub_file()

# TODO: Mass index table
for i, tmplt in enumerate(tmplt_bnk):
    mass_grouping = "MASS_SET_%d" % i

    ile_node = pipeline.CondorDAGNode(ile_job_type)
    ile_node.set_priority(JOB_PRIORITIES["ILE"])
    ile_node.add_macro("macromass1", tmplt.mass1)
    ile_node.add_macro("macromass2", tmplt.mass2)
    if opts.write_eff_lambda:
        ile_node.add_macro("macroefflambda", tmplt.psi0)
    if opts.write_deff_lambda:
        ile_node.add_macro("macrodefflambda", tmplt.psi3)
    if hasattr(tmplt, "spin1z"):
        ile_node.add_macro("macrospin1z", tmplt.spin1z)
    if hasattr(tmplt, "spin2z"):
        ile_node.add_macro("macrospin2z", tmplt.spin2z)
    if use_bayespe_postproc:
        # If we're using the Bayesian PE post processing script, dump the data
        ile_node.set_post_script(dagutils.which("process_ile_output"))
        ile_node.add_post_script_arg("--output ILE_%s.txt" % mass_grouping)
        ile_node.add_post_script_arg("--glob *-%s-*.xml.gz" % mass_grouping)

    # This is to identify output from groupings of the same mass point
    ile_node.add_macro("macromassid", mass_grouping)

    ile_node.set_category("ILE")
    dag.add_node(ile_node)

    sql_node = pipeline.CondorDAGNode(sql_job_type)
    # FIXME: Uncomment this line and the next to reconnect the subdags
    #if use_bayespe_postproc:
        #bpp_plot_node.add_parent(ile_node)
    #sql_node.add_parent(ile_node)
    sql_node.set_priority(JOB_PRIORITIES["SQL"])

    # The sql node needs to run a PRE script in order to coalesce the data into
    # a cache
    sql_node.set_pre_script(dagutils.which("coalesce.sh"))
    sql_node.add_pre_script_arg(mass_grouping)

    # This is to identify output from groupings of the sane mass point
    sql_node.add_macro("macromassid", mass_grouping)

    sql_node.set_category("SQL")
    ppdag.add_node(sql_node)

    if use_ile_postproc:
        tri_plot_job_type, tri_plot_job_name = dagutils.write_tri_plot_sub(tag="tri_plot", injection_file=opts.sim_xml, log_dir=opts.log_directory)
        tri_plot_job_type.write_sub_file()
        tri_plot_node = pipeline.CondorDAGNode(tri_plot_job_type)
        tri_plot_node.add_macro("macromassid", mass_grouping)
        tri_plot_node.set_category("PLOT")
        tri_plot_node.set_priority(JOB_PRIORITIES["PLOT"])
        ppdag.add_node(tri_plot_node)
        tri_plot_node.add_parent(sql_node)
    
        # In the interest of not blocking later DAGs for completion in an
        # uberdag, this is now dependent on the SQL step
        pos_plot_node.add_parent(sql_node)

dag_name=opts.output_name
dag.set_dag_file(dag_name)
dag.write_concrete_dag()
if opts.write_script:
    dag.write_script()

print "Created a DAG named %s\n" % dag_name
print "This will run %i instances of %s in parallel\n" % (len(tmplt_bnk), ile_sub_name)

# FIXME: Adjust name on command line
if use_bayespe_postproc:
    ppdag_name="posterior_pp"
    ppdag.set_dag_file(ppdag_name)
    ppdag.add_maxjobs_category("ANALYSIS", MAXJOBS["ANALYSIS"])
    ppdag.add_maxjobs_category("POST", MAXJOBS["POST"])
    ppdag.write_concrete_dag()
    if opts.write_script:
        ppdag.write_script()

    print "Created a postprocessing DAG named %s\n" % ppdag_name

xmldoc = ligolw.Document()
xmldoc.appendChild(ligolw.LIGO_LW())
process.register_to_xmldoc(xmldoc, sys.argv[0], opts.__dict__)
utils.write_filename(xmldoc, opts.output_name + ".xml.gz", gz=True)
