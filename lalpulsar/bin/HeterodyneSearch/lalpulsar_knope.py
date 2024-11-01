# -*- coding: utf-8 -*-
#
#       lalpulsar_knope.py
#
#       Copyright 2015
#       Matthew Pitkin <matthew.pitkin@ligo.org>,
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

## \file
## \ingroup lalpulsar_bin_HeterodyneSearch
"""
The KNOwn Pulsar pipelinE - lalpulsar_knope

DAG generation code for running the known pulsar search pipeline
"""

from __future__ import print_function, division

from lalpulsar import knope_utils as knope
import argparse
from configparser import ConfigParser
import sys
import pickle

description = """Setup a Condor DAG file to run the known pulsar search pipeline based on information given in config.ini.
The user must specify the configuration file for the script to run.
"""

parser = argparse.ArgumentParser(description=description)

parser.add_argument(
    "inifile", help="The configuation (.ini) file"
)  # the positional argument for the configuration file
parser.add_argument(
    "--condor-submit",
    action="store_true",
    default=False,
    help="Automatically submit the Condor DAG",
)
parser.add_argument(
    "-r",
    "--run-path",
    dest="runpath",
    default=None,
    help="Set the directory to run the pipeline in (overwrites any value in the config.ini file)",
)
parser.add_argument(
    "-p",
    "--pulsar",
    dest="pulsarlist",
    action="append",
    default=None,
    help="A pulsar name to search for rather than all pulsars given in a parameter file directory (this can be specified multiple times to search for more than one pulsar).",
)

opts = parser.parse_args()

# check that at least the ini file has been given
inifile = opts.inifile

# parser .ini file
try:
    cp = ConfigParser()
    cp.optionxform = str
    cp.readfp(open(inifile))
except:
    print("Error... problem parsing '%s' configuration file" % inifile, file=sys.stderr)
    sys.exit(1)

if opts.runpath is not None:
    cp.set("analysis", "run_dir", opts.runpath)

# Check if we're running in automated mode or not
try:
    automated = cp.getboolean("analysis", "autonomous")
except:
    automated = False

# Check if configuration file says to submit the DAG
submitdag = opts.condor_submit
if not submitdag:
    try:
        submitdag = cp.getboolean("analysis", "submit_dag")
    except:
        submitdag = False

# Create DAG from ConfigParser object
dag = knope.knopeDAG(cp, inifile, pulsarlist=opts.pulsarlist)
if dag.error_code != 0:  # check for any errors that occurred
    if dag.error_code in knope.KNOPE_ERROR.keys():
        print(knope.KNOPE_ERROR[dag.error_code], file=sys.stderr)
    else:
        print("Error... unrecognised error code!", file=sys.stderr)

    # only exit if not in autonomous mode and the error message
    if not automated or dag.error_code != knope.KNOPE_ERROR_NO_SEGMENTS:
        sys.exit(dag.error_code)

# write out DAG and submit files (unless in automated mode and no new segment files were found)
if not automated or dag.error_code != knope.KNOPE_ERROR_NO_SEGMENTS:
    dag.write_sub_files()
    dag.write_dag()

    print("Successfully created DAG file: '%s'" % dag.get_dag_file())

    if submitdag:
        from subprocess import Popen

        x = Popen(["condor_submit_dag", dag.get_dag_file()])
        x.wait()
        if x.returncode == 0:
            print("Submitted DAG file")
        else:
            print("Unable to submit DAG file")
    else:
        print("Run 'condor_submit_dag %s' to submit DAG file" % dag.get_dag_file())

    # output DAG class to pickle file if given
    if cp.has_option("analysis", "pickle_file"):
        try:
            pfile = cp.get("analysis", "pickle_file")
            fp = open(pfile, "wb")
            pickle.dump(dag, fp)
            fp.close()
        except:
            print("Warning... could not output analysis class to pickle file")
else:
    print("No new science segments found in current time frame")

sys.exit(0)
