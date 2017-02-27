# Copyright (C) 2016  Surabhi Sachdev, Tjonnie Li, Stephen Privitera
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


###############################################################################
#
# IMPORT MODULES
#
###############################################################################

import os
import socket, tempfile
from optparse import OptionParser
from glue.pipeline import DeepCopyableConfigParser as dcConfigParser
from glue import pipeline
from lalapps import inspiral

###############################################################################
#
# CLASSES AND FUNCTIONS
#
###############################################################################

def log_path():
    host = socket.getfqdn()
    username = os.environ['USER']
    #FIXME add more hosts as you need them
    if 'caltech.edu' in host: return '/usr1/' + username
    if 'phys.uwm.edu' in host: return '/localscratch/' + username
    if 'aei.uni-hannover.de' in host: return '/local/user/' + username
    if 'phy.syr.edu' in host: return '/usr1/' + username
    else: return os.environ["TMPDIR"]

class bank_DAG(pipeline.CondorDAG):

    def __init__(self, name, logpath = log_path()):
        self.basename = name
        tempfile.tempdir = logpath
        tempfile.template = self.basename + '.dag.log.'
        logfile = tempfile.mktemp()
        fh = open( logfile, "w" )
        fh.close()
        pipeline.CondorDAG.__init__(self,logfile, dax=False)
        self.set_dag_file(self.basename)
        self.set_dax_file(self.basename)
        self.jobsDict = {}
        self.node_id = 0
        self.output_cache = []

    def add_node(self, node):
        node.set_retry(3)
        self.node_id += 1
        node.add_macro("macroid", self.node_id)
        pipeline.CondorDAG.add_node(self, node)

    def write_cache(self):
        out = self.basename + ".cache"
        f = open(out,"w")
        for c in self.output_cache:
            f.write(str(c)+"\n")
        f.close()

class InspinjJob(inspiral.InspiralAnalysisJob):
    def __init__(self, cp, dax=False, tag_base="inspinj"):
        exec_name = 'lalapps_inspinj'
        extension = 'xml'
        sections = ['inspinj']
        inspiral.InspiralAnalysisJob.__init__(self,cp,sections,exec_name,extension,dax)
        self.set_universe("local")
        self.set_sub_file(tag_base+'.sub')
        self.tag_base = tag_base
        self.set_stdout_file('logs/'+tag_base+'-$(macroid)-$(process).out')
        self.set_stderr_file('logs/'+tag_base+'-$(macroid)-$(process).err')
        if cp.has_section("accounting"):
            self.add_condor_cmd('accounting_group', cp.get("accounting", "accounting-group"))


class InspinjNode(pipeline.CondorDAGNode):
    def __init__(self, job, dag, tag=None, seed=0, p_node=[]):
        pipeline.CondorDAGNode.__init__(self,job)
        self.add_var_opt("seed", seed)
        outputFile = "HL-INJECTIONS_" + str(seed)
        if tag is not None:
            outputFile += "_" + tag
        opts = self.job().get_opts()
        outputFile += "-" + opts["gps-start-time"] + "-" + str(int(opts["gps-end-time"]) - int(opts["gps-start-time"])) + ".xml"
        self.add_var_opt("output", outputFile)
        self.add_output_file(outputFile)
        for p in p_node:
            self.add_parent(p)
        dag.add_node(self)

class BankSimJob(inspiral.InspiralAnalysisJob):
    def __init__(self, cp, dax=False, tag_base ="banksim"):
        exec_name = 'lalapps_cbc_sbank_sim'
        extension = 'xml'
        sections = ['banksim']
        inspiral.InspiralAnalysisJob.__init__(self,cp,sections,exec_name,extension,dax)
        self.add_condor_cmd('request_memory', '3999')
        self.tag_base = tag_base
        self.add_condor_cmd('getenv','True')
        self.set_stdout_file('logs/'+tag_base+'-$(macroid)-$(process).out')
        self.set_stderr_file('logs/'+tag_base+'-$(macroid)-$(process).err')
        if cp.has_section("accounting"):
            self.add_condor_cmd('accounting_group', cp.get("accounting", "accounting-group"))
        if "OMP_NUM_THREADS" in os.environ:
            self.add_condor_cmd('request_cpus', os.environ["OMP_NUM_THREADS"])


class BankSimNode(pipeline.CondorDAGNode):
    def __init__(self, job, dag, injs, approx, injmin=0, injmax=-1, tag="SIM", p_node=[]):
        pipeline.CondorDAGNode.__init__(self,job)
        self.add_file_opt("injection-file", injs)
        self.add_var_opt("injection-approx", approx)
        self.add_var_opt("injection-min", injmin)
        self.add_var_opt("injection-max", injmax)
        self.add_var_opt("user-tag", tag)
        self.add_output_file("%s.h5" % tag)
        for p in p_node:
            if p is not None:
                self.add_parent(p)
        dag.add_node(self)


class PlotSimJob(inspiral.InspiralAnalysisJob):
    def __init__(self, cp, dax=False, tag_base="lalapps_cbc_sbank_plot_sim"):
        exec_name = 'lalapps_cbc_sbank_plot_sim'
        extension = 'h5'
        sections = []
        inspiral.InspiralAnalysisJob.__init__(self,cp,sections,exec_name,extension,dax)
        self.tag_base = tag_base
        self.set_universe("local")
        self.set_sub_file(tag_base+'.sub')
        self.set_stdout_file('logs/'+tag_base+'-$(macroid)-$(process).out')
        self.set_stderr_file('logs/'+tag_base+'-$(macroid)-$(process).err')
        self.add_condor_cmd('getenv','True')
        if cp.has_section("accounting"):
            self.add_condor_cmd('accounting_group', cp.get("accounting", "accounting-group"))
        self.add_condor_cmd('request_memory', '3999')


class PlotSimNode(pipeline.CondorDAGNode):
    def __init__(self, job, dag, h5files=[],p_node=[]):
        pipeline.CondorDAGNode.__init__(self,job)
        [self.add_file_arg(h5f) for h5f in h5files]
        for p in p_node:
            self.add_parent(p)
        dag.add_node(self)


usage = """

lalapps_cbc_sbank_sim_pipe generates a DAG for parallelizing the
 verification of template banks via lalapps_cbc_sbank_sim. Its primary
 input is an ini file, which takes form given below. Once you have an
 ini file, create the DAG by

lalapps_cbc_sbank_sim_pipe --config-file INIFILE.ini --user-tag SBANK_PIPE

;;; BEGIN EXAMPLE INI ;;;

[condor] ; This section points to the executables, and provides condor options
universe = vanilla
lalapps_cbc_sbank_sim = /home/sprivite/local/opt/master/bin/lalapps_cbc_sbank_sim
lalapps_cbc_sbank_plot_sim = /home/sprivite/local/opt/master/bin/lalapps_cbc_sbank_plot_sim
lalapps_inspinj = /home/sprivite/local/opt/master/bin/lalapps_inspinj

;[accounting]
;accounting-group = ???

[banksim]
; This section contains the parameters of the entire bank parameter
; space you wish to cover. sbank_pipe will divide the space for you.
template-bank = MYBANK1[:MYAPPROX1]
;template-bank = MYBANK2[:MYAPPROX2]
noise-model = aLIGOZeroDetHighPower
; reference-psd = PSD.xml.gz
flow = 30
; extra filters on injection parameters produced by inspinj, if needed
instrument = L1
mtotal-max = 400
mtotal-min = 4
mratio-max = 97
duration-min = 0.2
neighborhood-size = 2.0

;
; FOR BANK SIMS
;
[injections]
; As in the classic IHOPE style, this section is a list of injection
; runs in the format "name" = "seed".
seobnr = 1234
imrc = 1234

[inspinj]
; This section gives options that are common to all injection runs,
; regardless of name or seed. Be careful, since flag options (those
; that set a value to true) cannot be overridden below.
f-lower = 25.0
waveform = SEOBNRv2_ROM_DoubleSpin
l-distr = random
d-distr = volume
i-distr = uniform
max-inc = 179.99
polarization = uniform
min-distance = 10000
max-distance = 1000000
m-distr = componentMass
taper-injection = startend
min-mass1 = 2.0
max-mass1 = 296.0
min-mass2 = 2.0
max-mass2 = 296.0

[seobnr]
; This section overrides the generic inspinj options from the section
; above. If an option is not included here, it will be taken from the
; [inspinj] section. The options in this section apply only to the
; [seobnr] run, in this example.
waveform = SEOBNRv4_ROM
enable-spin =
min-spin1 = 0.0
max-spin1 = 0.98
min-spin2 = 0.0
max-spin2 = 0.98
aligned =

[imrd]
; This section overrides the generic inspinj options from the section
; above. If an option is not included here, it will be taken from the
; [inspinj] section. The options in this section apply only to the
; [imrc] run, in this example.
waveform = IMRPhenomD
enable-spin =
min-spin1 = 0.0
max-spin1 = 0.9
min-spin2 = 0.0
max-spin2 = 0.9
aligned =

;;; END EXAMPLE INI ;;;

"""

###############################################################################
#
# OPTION PARSER
#
###############################################################################

def parse_command_line():
    parser = OptionParser(usage=usage)
    parser.add_option("--config-file", default=None, help="Read options for generating template bank placement pipeline from configuration ini file..")
    parser.add_option("--user-tag", default="SBANK", help="Make your results feel special and give them a unique name.")
    options, filenames = parser.parse_args()

    return options, filenames

options, filenames = parse_command_line()


###############################################################################
#
# CREATE DAG
#
###############################################################################

# READ IN CONFIGPARSER
try: os.mkdir("logs")
except: pass

# read config file
cp = dcConfigParser()
cp.read(options.config_file)

# CREATE TOP LEVEL DAG AND JOBS
dag = bank_DAG(options.user_tag)
banksimJob = BankSimJob(cp)
plotsimJob = PlotSimJob(cp)

# for each injection run and every subbank, perform a banksim
for inj_run in cp.options("injections"):
    # set up inspinj jobs
    # each option in the [injections] section is a separate injection run
    # the inspinj section sets options that are common to all injection runs
    injperjob = int(cp.get("sbank_pipe", "injections-per-job"))
    njobs = int(cp.get("sbank_pipe", "jobs"))

    gpsstart = 1000000000 # irrelevant for banksim.
    # because we will prune injections, we generate 50 times more
    # than we think we need. this is fast.
    # FIXME: this is still no guarantee we have enough injections!!
    gpsend = gpsstart + 50 * (njobs * injperjob)
    cp.set("inspinj", "gps-start-time", str(gpsstart))
    cp.set("inspinj", "gps-end-time", str(gpsend))
    cp.set("inspinj", "time-step", str(1))

    # override default inspinj options with those given in the inj_run section
    for opt in cp.options(inj_run):
        cp.set('inspinj', opt, cp.get(inj_run, opt))

    # set random seed
    seed = int(cp.get("injections", inj_run))
    cp.set('inspinj', 'seed', str(seed))
    tag = options.user_tag + "_INJECTIONS_" + inj_run.upper()

    # inspinj nodes
    inspinjJob = InspinjJob(cp, tag_base = tag.lower())
    inj_node = InspinjNode(inspinjJob, dag, tag=tag,seed = seed, p_node = [])

    # reset the config parser option back to default values
    for opt in cp.options(inj_run):
        cp.remove_option('inspinj',opt)
    cp.read(options.config_file)

    # add sim nodes
    inj_name, = inj_node.get_output_files()
    waveform = cp.get(inj_run, "waveform")
    sim_nodes = []
    for j in range(njobs):

        tag = "%05d" % j
        sim_name = "SBANK_SIM_%s_%s_%05d" % (options.user_tag, inj_run.upper(), j)
        sim_nodes.append(BankSimNode(banksimJob, dag, inj_name, waveform, j*injperjob, (j+1)*injperjob, sim_name, [inj_node]))

    # merge and plot the partial sims
    sim_names = [node.get_output_files()[0] for node in sim_nodes]
    inputs = []
    [inputs.extend(sim_node.get_output_files()[:]) for sim_node in sim_nodes]
    PlotSimNode(plotsimJob, dag, inputs, sim_nodes)

# write the dag
dag.write_sub_files()
dag.write_script()
dag.write_cache()
dag.write_dag()
