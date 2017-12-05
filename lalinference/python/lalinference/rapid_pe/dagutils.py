# Copyright (C) 2013, 2015  Evan Ochsner, Chris Pankow
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
A collection of routines to manage Condor workflows (DAGs).
"""

import os

from glue import pipeline

__author__ = "Evan Ochsner <evano@gravity.phys.uwm.edu>, Chris Pankow <pankow@gravity.phys.uwm.edu>"

# Taken from
# http://pythonadventures.wordpress.com/2011/03/13/equivalent-of-the-which-command-in-python/
def is_exe(fpath):
    return os.path.exists(fpath) and os.access(fpath, os.X_OK)

def which(program):
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file): return exe_file

    return None

# FIXME: Keep in sync with arguments of integrate_likelihood_extrinsic
def write_integrate_likelihood_extrinsic_sub(tag='integrate', exe=None, log_dir=None, intr_prms=("mass1", "mass2"), ncopies=1, condor_commands=None, **kwargs):
    """
    Write a submit file for launching jobs to marginalize the likelihood over
    extrinsic parameters.

    Inputs:
        - 'tag' is a string to specify the base name of output files. The output
          submit file will be named tag.sub, and the jobs will write their
          output to tag-ID.out, tag-ID.err, tag.log, where 'ID' is a unique
          identifier for each instance of a job run from the sub file.
        - 'cache' is the path to a cache file which gives the location of the
          data to be analyzed.
        - 'coinc' is the path to a coincident XML file, from which masses and
          times will be drawn FIXME: remove this once it's no longer needed.
        - 'channelH1/L1/V1' is the channel name to be read for each of the
          H1, L1 and V1 detectors.
        - 'psdH1/L1/V1' is the path to an XML file specifying the PSD of
          each of the H1, L1, V1 detectors.
        - 'ncopies' is the number of runs with identical input parameters to
          submit per condor 'cluster'

    Outputs:
        - An instance of the CondorDAGJob that was generated for ILE
    """

    assert len(kwargs["psd_file"]) == len(kwargs["channel_name"])

    exe = exe or which("rapidpe_integrate_extrinsic_likelihood")
    ile_job = pipeline.CondorDAGJob(universe="vanilla", executable=exe)
    # This is a hack since CondorDAGJob hides the queue property
    ile_job._CondorJob__queue = ncopies

    ile_sub_name = tag + '.sub'
    ile_job.set_sub_file(ile_sub_name)

    #
    # Logging options
    #
    uniq_str = "$(macromassid)-$(cluster)-$(process)"
    ile_job.set_log_file("%s%s-%s.log" % (log_dir, tag, uniq_str))
    ile_job.set_stderr_file("%s%s-%s.err" % (log_dir, tag, uniq_str))
    ile_job.set_stdout_file("%s%s-%s.out" % (log_dir, tag, uniq_str))

    if kwargs.has_key("output_file") and kwargs["output_file"] is not None:
        #
        # Need to modify the output file so it's unique
        #
        ofname = kwargs["output_file"].split(".")
        ofname, ext = ofname[0], ".".join(ofname[1:])
        ile_job.add_file_opt("output-file", "%s-%s.%s" % (ofname, uniq_str, ext))
        del kwargs["output_file"]
        if kwargs.has_key("save_samples") and kwargs["save_samples"] is True:
            ile_job.add_opt("save-samples", '')
            del kwargs["save_samples"]

    #
    # Add normal arguments
    # FIXME: Get valid options from a module
    #
    for opt, param in kwargs.iteritems():
        if isinstance(param, list) or isinstance(param, tuple):
            # NOTE: Hack to get around multiple instances of the same option
            for p in param:
                ile_job.add_arg("--%s %s" % (opt.replace("_", "-"), str(p)))
        elif param is True or param is None:
            ile_job.add_opt(opt.replace("_", "-"), '')
        # Explcitly check for False to turn it off
        elif param is False:
            continue
        else:
            ile_job.add_opt(opt.replace("_", "-"), str(param))

    #
    # Macro based options
    #
    ile_job.add_var_opt("mass1")
    ile_job.add_var_opt("mass2")
    for p in intr_prms:
        ile_job.add_var_opt(p.replace("_", "-"))

    ile_job.add_condor_cmd('getenv', 'True')
    if condor_commands is not None:
        for cmd, value in condor_commands.iteritems():
            ile_job.add_condor_cmd(cmd, value)
    ile_job.add_condor_cmd('request_memory', '2048')
    
    return ile_job, ile_sub_name

def write_result_coalescence_sub(tag='coalesce', exe=None, log_dir=None, output_dir="./", use_default_cache=True):
    """
    Write a submit file for launching jobs to coalesce ILE output
    """

    exe = exe or which("ligolw_sqlite")
    sql_job = pipeline.CondorDAGJob(universe="vanilla", executable=exe)

    sql_sub_name = tag + '.sub'
    sql_job.set_sub_file(sql_sub_name)

    #
    # Logging options
    #
    uniq_str = "$(cluster)-$(process)"
    sql_job.set_log_file("%s%s-%s.log" % (log_dir, tag, uniq_str))
    sql_job.set_stderr_file("%s%s-%s.err" % (log_dir, tag, uniq_str))
    sql_job.set_stdout_file("%s%s-%s.out" % (log_dir, tag, uniq_str))

    if use_default_cache:
        sql_job.add_opt("input-cache", "ILE_$(macromassid).cache")
    else:
        sql_job.add_arg("$(macrofiles)")
    #sql_job.add_arg("*$(macromassid)*.xml.gz")
    sql_job.add_opt("database", "ILE_$(macromassid).sqlite")
    #if os.environ.has_key("TMPDIR"):
        #tmpdir = os.environ["TMPDIR"]
    #else:
        #print >>sys.stderr, "WARNING, TMPDIR environment variable not set. Will default to /tmp/, but this could be dangerous."
        #tmpdir = "/tmp/"
    tmpdir = "/dev/shm/"
    sql_job.add_opt("tmp-space", tmpdir)
    sql_job.add_opt("verbose", '')

    sql_job.add_condor_cmd('getenv', 'True')
    sql_job.add_condor_cmd('request_memory', '1024')
    
    return sql_job, sql_sub_name

def write_posterior_plot_sub(tag='plot_post', exe=None, log_dir=None, output_dir="./"):
    """
    Write a submit file for launching jobs to coalesce ILE output
    """

    exe = exe or which("plot_like_contours")
    plot_job = pipeline.CondorDAGJob(universe="vanilla", executable=exe)

    plot_sub_name = tag + '.sub'
    plot_job.set_sub_file(plot_sub_name)

    #
    # Logging options
    #
    uniq_str = "$(cluster)-$(process)"
    plot_job.set_log_file("%s%s-%s.log" % (log_dir, tag, uniq_str))
    plot_job.set_stderr_file("%s%s-%s.err" % (log_dir, tag, uniq_str))
    plot_job.set_stdout_file("%s%s-%s.out" % (log_dir, tag, uniq_str))

    plot_job.add_opt("show-points", '')
    plot_job.add_opt("dimension1", "mchirp")
    plot_job.add_opt("dimension2", "eta")
    plot_job.add_opt("input-cache", "ILE_all.cache")
    plot_job.add_opt("log-evidence", '')

    plot_job.add_condor_cmd('getenv', 'True')
    plot_job.add_condor_cmd('request_memory', '1024')
    
    return plot_job, plot_sub_name

def write_tri_plot_sub(tag='plot_tri', injection_file=None, exe=None, log_dir=None, output_dir="./"):
    """
    Write a submit file for launching jobs to coalesce ILE output
    """

    exe = exe or which("make_triplot")
    plot_job = pipeline.CondorDAGJob(universe="vanilla", executable=exe)

    plot_sub_name = tag + '.sub'
    plot_job.set_sub_file(plot_sub_name)

    #
    # Logging options
    #
    uniq_str = "$(cluster)-$(process)"
    plot_job.set_log_file("%s%s-%s.log" % (log_dir, tag, uniq_str))
    plot_job.set_stderr_file("%s%s-%s.err" % (log_dir, tag, uniq_str))
    plot_job.set_stdout_file("%s%s-%s.out" % (log_dir, tag, uniq_str))

    plot_job.add_opt("output", "ILE_triplot_$(macromassid).png")
    if injection_file is not None:
        plot_job.add_opt("injection", injection_file)
    plot_job.add_arg("ILE_$(macromassid).sqlite")

    plot_job.add_condor_cmd('getenv', 'True')
    #plot_job.add_condor_cmd('request_memory', '2048')
    
    return plot_job, plot_sub_name

def write_1dpos_plot_sub(tag='1d_post_plot', exe=None, log_dir=None, output_dir="./"):
    """
    Write a submit file for plotting 1d posterior cumulants.
    """

    exe = exe or which("postprocess_1d_cumulative")
    plot_job = pipeline.CondorDAGJob(universe="vanilla", executable=exe)

    plot_sub_name = tag + '.sub'
    plot_job.set_sub_file(plot_sub_name)

    #
    # Logging options
    #
    uniq_str = "$(cluster)-$(process)"
    plot_job.set_log_file("%s%s-%s.log" % (log_dir, tag, uniq_str))
    plot_job.set_stderr_file("%s%s-%s.err" % (log_dir, tag, uniq_str))
    plot_job.set_stdout_file("%s%s-%s.out" % (log_dir, tag, uniq_str))

    plot_job.add_opt("save-sampler-file", "ILE_$(macromassid).sqlite")
    plot_job.add_opt("disable-triplot", '')
    plot_job.add_opt("disable-1d-density", '')

    plot_job.add_condor_cmd('getenv', 'True')
    plot_job.add_condor_cmd('request_memory', '2048')
    
    return plot_job, plot_sub_name

def write_bayes_pe_postproc_sub(tag='bayespe_post_plot', exe=None, log_dir=None, web_dir="./", inj_xml=None):
    """
    Write a submit file for postprocessing output and pushing it through cbcBayesPostProc.py
    """

    exe = exe or which("cbcBayesPostProc.py")
    plot_job = pipeline.CondorDAGJob(universe="vanilla", executable=exe)

    plot_sub_name = tag + '.sub'
    plot_job.set_sub_file(plot_sub_name)

    #
    # Logging options
    #
    uniq_str = "$(cluster)-$(process)"
    plot_job.set_log_file("%s%s-%s.log" % (log_dir, tag, uniq_str))
    plot_job.set_stderr_file("%s%s-%s.err" % (log_dir, tag, uniq_str))
    plot_job.set_stdout_file("%s%s-%s.out" % (log_dir, tag, uniq_str))

    #
    # Injection options
    #
    plot_job.add_opt("outpath", web_dir)
    if inj_xml:
        plot_job.add_opt("inj", inj_xml)
        # FIXME: Since we put individual sim entries into their own XML, this is
        # always zero. We might need to tweak this if we use a bigger one
        plot_job.add_opt("eventnum", 0)

    # Calculate evidence (just to compare)
    plot_job.add_opt("dievidence", '')

    plot_job.add_opt("header", "header.txt")
    plot_job.add_opt("data", "tmp")

    plot_job.add_condor_cmd('getenv', 'True')
    plot_job.add_condor_cmd('request_memory', '1024')
    
    return plot_job, plot_sub_name
