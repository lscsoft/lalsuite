# Copyright (C) 2014  Stephen Privitera
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


##############################################################################
# import standard modules and append the lalapps prefix to the python path
import sys, os, shutil
from itertools import cycle, izip
import subprocess, socket, tempfile

##############################################################################
# import the modules we need to build the pipeline
from optparse import OptionParser
from glue.pipeline import DeepCopyableConfigParser as dcConfigParser

from glue import pipeline
from lalapps import inspiral

def which(prog):
    which = subprocess.Popen(['/usr/bin/which', prog], stdout=subprocess.PIPE)
    out = which.stdout.read().strip()
    if not out:
        print >>sys.stderr, "ERROR: could not find %s in your path, have you built the proper software and source the proper env. scripts?" % (prog,prog)
        raise ValueError
    return out

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


class SBankJob(inspiral.InspiralAnalysisJob):
    def __init__(self,cp,dax=False, tag_base="lalapps_cbc_sbank"):
        exec_name = 'lalapps_cbc_sbank'
        extension = 'xml'
        sections = ['sbank']
        inspiral.InspiralAnalysisJob.__init__(self,cp,sections,exec_name,extension,dax)
        self.tag_base = tag_base
        self.set_sub_file(tag_base+'.sub')
        self.set_stdout_file('logs/'+tag_base+'-$(macroid)-$(process).out')
        self.set_stderr_file('logs/'+tag_base+'-$(macroid)-$(process).err')
        if cp.has_section("accounting"):
            self.add_condor_cmd('accounting_group', cp.get("accounting", "accounting-group"))
        self.add_condor_cmd('getenv','True')
        self.add_condor_cmd('request_memory', '3999')
        if "OMP_NUM_THREADS" in os.environ:
            self.add_condor_cmd('request_cpus', os.environ["OMP_NUM_THREADS"])


class SBankNode(pipeline.CondorDAGNode):
    def __init__(self, job, dag, tag="SBANK", seed=0, bank_seed=[], mchirp_boundaries_file=None, mchirp_boundaries_index=None, p_node=[]):
        pipeline.CondorDAGNode.__init__(self,job)
        if (mchirp_boundaries_file is None) ^ (mchirp_boundaries_index is None):
            raise ValueError("must supply both mchirp_boundaries_file and mchirp_boundaries_index or neither")
        if mchirp_boundaries_file is not None:
            self.add_var_arg("--mchirp-boundaries-file %s --mchirp-boundaries-index %s" % (mchirp_boundaries_file, mchirp_boundaries_index))
        for bseed in bank_seed:
            self.add_var_arg("--bank-seed %s" % bseed)
        self.add_var_opt("seed", seed)

        fout = "%s.xml.gz" % tag
        self.add_var_opt("output-filename", fout)
        self.add_output_file(fout)

        for p in p_node:
            self.add_parent(p)

        dag.add_node(self)


class SBankChooseMchirpBoundariesJob(inspiral.InspiralAnalysisJob):
    def __init__(self,cp,dax=False, tag_base="lalapps_cbc_sbank_choose_mchirp_boundaries"):
        exec_name = 'lalapps_cbc_sbank_choose_mchirp_boundaries'
        extension = 'txt'
        sections = ['split']
        inspiral.InspiralAnalysisJob.__init__(self,cp,sections,exec_name,extension,dax)
        self.set_universe("local")
        self.set_sub_file(tag_base+'.sub')
        self.tag_base = tag_base
        self.set_stdout_file('logs/'+tag_base+'-$(macroid)-$(process).out')
        self.set_stderr_file('logs/'+tag_base+'-$(macroid)-$(process).err')
        if cp.has_section("accounting"):
            self.add_condor_cmd('accounting_group', cp.get("accounting", "accounting-group"))
        self.add_condor_cmd('getenv','True')


class SBankChooseMchirpBoundariesNode(pipeline.CondorDAGNode):
    def __init__(self, job, dag, bank, tag = None, p_node=[]):
        pipeline.CondorDAGNode.__init__(self,job)
        self.add_file_arg(bank)
        nbanks = int(self.job().get_opts()["nbanks"])
        self.add_var_opt("nbanks", nbanks)
        if tag:
            output_fname = "%s_mchirp_boundaries.txt" % tag
        else:
            output_fname = "mchirp_boundaries.txt"
        self.add_var_opt("output-file", output_fname)
        self.add_output_file(output_fname)
        for p in p_node:
            self.add_parent(p)
        dag.add_node(self)


class LWAddJob(pipeline.CondorDAGJob):
    """
    A ligolw_add node. This node is used to combine the split template banks
    into one aggregate bank.
    """
    def __init__(self, cp, executable=which('ligolw_add'), tag_base='ligolw_add'):
        """
        """
        self.__prog__ = 'ligolw_add'
        self.__executable = executable
        self.__universe = 'vanilla'
        pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
        self.add_condor_cmd('request_memory', '3999')
        self.add_condor_cmd('getenv','True')
        if cp.has_section("accounting"):
            self.add_condor_cmd('accounting_group', cp.get("accounting", "accounting-group"))
        self.tag_base = tag_base
        self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
        self.set_sub_file(tag_base+'.sub')
        self.set_stdout_file('logs/'+tag_base+'-$(macroid)-$(process).out')
        self.set_stderr_file('logs/'+tag_base+'-$(macroid)-$(process).err')


class LWAddNode(pipeline.CondorDAGNode):
    """
    """
    def __init__(self, job, dag, xmls,output,p_node=[]):

        pipeline.CondorDAGNode.__init__(self,job)
        for x in xmls: self.add_file_arg(x)
        self.add_var_opt("output", output)
        self.add_output_file(output)
        for p in p_node:
            if p is not None:
                self.add_parent(p)
        dag.add_node(self)


usage = """

lalapps_cbc_sbank_pipe generates a DAG for parallelizing the
construction of template banks via lalapps_cbc_sbank. Its primary
input is an ini file, which takes form given below. Once you have an
ini file, create the DAG by

lalapps_cbc_sbank_pipe --config-file INIFILE.ini --user-tag SBANK_PIPE

;;; BEGIN EXAMPLE INI ;;;

[condor] ; This section points to the executables, and provides condor options
universe = vanilla
lalapps_cbc_sbank = /home/sprivite/local/opt/master/bin/lalapps_cbc_sbank
lalapps_cbc_sbank_choose_mchirp_boundaries = /home/sprivite/local/opt/master/bin/lalapps_cbc_sbank_choose_mchirp_boundaries

;[accounting]
;accounting-group = ???

[sbank]
; This section contains the parameters of the entire bank parameter
; space you wish to cover. sbank_pipe will divide the space for you.
approximant = SEOBNRv2_ROM_DoubleSpin
match-min = 0.97
flow = 30.0
noise-model = aLIGOZeroDetHighPower
instrument = H1
mass1-min = 2.0
mass1-max = 296.0
mass2-min = 2.0
mass2-max = 296.0
mtotal-min = 4.0
mtotal-max = 350.0
mratio-min = 1.0
mratio-max = 14.0
spin1-max = 0.98
spin1-min = -0.98
spin2-max = 0.98
spin2-min = -0.98
aligned-spin =
gps-start-time = 1000000000
gps-end-time =  1000050000
neighborhood-param = dur
neighborhood-size = 0.5
checkpoint = 500
;cache-waveforms =
; please check your memory requirements before setting this option
; if jobs use too much memory, they will swap and run forever
coarse-match-df = 4.0
iterative-match-df-max= 1.0
; If you want to add precomputed banks to the workflow
;bank-seed = FILENAME[:APPROX]

[coarse-sbank]
; This section is for planning the splitting of the parameter
; space. To do so, we generate a "coarse" bank, i.e., a bank on the
; same parameter space but with a weaker convergence criteria. This
; process gives a rough measure of the density of templates the final
; bank will require. Use the max-new-templates option to prevent the job
; from running forever, but the more templates you have in the coarse
; bank, the less "over-coverage" you will incur from the bank
; splitting process. A good rule of thumb is that you want ~1000
; templates per split bank.
max-new-templates = 15000
match-min = 0.93
convergence-threshold = 1000

[split]
; This section configures the parallelization. nbanks are how many
; splits (in chirp mass) you want. You can crank it to infinity at the
; cost of overcoverage.
nbanks = 15

;;; END EXAMPLE INI ;;;

"""

def parse_command_line():
    parser = OptionParser(usage=usage)
    parser.add_option("--config-file", default=None, help="Read options for generating template bank placement pipeline from configuration ini file..")
    parser.add_option("--user-tag", default="SBANK", help="Make your results feel special and give them a unique name.")
    options, filenames = parser.parse_args()

    return options, filenames

options, filenames = parse_command_line()

try: os.mkdir("logs")
except: pass

# create top level dag object
dag = bank_DAG(options.user_tag)

# read config file
cp = dcConfigParser()
cp.read(options.config_file)

# initialize sbank job objects
sbankJob = SBankJob(cp)
mm = cp.get("sbank", "match-min")
cp.remove_option("sbank", "match-min")  # don't want it entering via add_ini_opts

# set up bank generation
# Two modes:
#   1. generate coarse, partition mchirp space, generate sub-banks, ligolw_add
#   2. take given bank, split it, and sim it
nbanks = int(cp.get("split","nbanks"))

# for layering, we prefer nbanks to be odd so that no second-stage job
# works on the bank boundary
if not (nbanks % 2):
    nbanks += 1
    cp.set("split", "nbanks", str(nbanks))

bank_names = []
bank_nodes = []

# set up sole coarse node to plan out the mini-sbank nodes
coarse_sbank_node = SBankNode(sbankJob, dag, "SBANK_COARSE")
coarse_mm = cp.get("coarse-sbank", "match-min")
coarse_sbank_node.add_var_opt("match-min", coarse_mm)
coarse_thresh = cp.get("coarse-sbank", "convergence-threshold")
coarse_sbank_node.add_var_arg("--convergence-threshold %s" % coarse_thresh)

if cp.has_option("coarse-sbank", "max-new-templates"):
    templates_max = cp.get("coarse-sbank", "max-new-templates")
    assert templates_max >= 3*nbanks # you need at least a few templates per bank
else:
    templates_max = 15*nbanks  # to prevent the coarse bank from running forever
coarse_sbank_node.add_var_arg("--max-new-templates %s" % templates_max)

xmlCoarse, = coarse_sbank_node.get_output_files()
pnode = [coarse_sbank_node]
bank_names.append(xmlCoarse)
bank_nodes.append(coarse_sbank_node)

# use coarse bank to choose mchirp regions of roughly equal template number
sbankChooseMchirpBoundariesJob = SBankChooseMchirpBoundariesJob(cp)
sbankChooseMchirpBoundariesNode = SBankChooseMchirpBoundariesNode(sbankChooseMchirpBoundariesJob, dag, xmlCoarse, options.user_tag, pnode)
mchirp_boundaries_fname, = sbankChooseMchirpBoundariesNode.get_output_files()

# generate a bank for each mchirp region
# first compute even numbered split banks
for j in xrange(0, nbanks, 2):

    bank_node = SBankNode(sbankJob, dag, "SBANK_SPLIT_BANK_%04d"%j, seed="%d" % (j*nbanks+1), mchirp_boundaries_file=mchirp_boundaries_fname, mchirp_boundaries_index=str(j), p_node=[sbankChooseMchirpBoundariesNode], bank_seed=[xmlCoarse])
    bank_node.add_var_opt("match-min", mm)
    bank_node.set_priority(1)  # want complete bank before sims
    bank_nodes.append(bank_node)
    bank_name, = bank_node.get_output_files()
    bank_names.append(bank_name)

# then compute odd numbered split banks using even banks as seeds
for j in xrange(1, nbanks, 2):

    if j < nbanks - 1:
        p_node = [bank_nodes[(j+1)/2], bank_nodes[(j+3)/2]]
        bank_seed = [xmlCoarse, bank_names[(j+1)/2], bank_names[(j+3)/2]]
    else:
        p_node = [bank_nodes[(j+1)/2]]
        bank_seed = [xmlCoarse, bank_names[(j+1)/2]]

    bank_node = SBankNode(sbankJob, dag, "SBANK_SPLIT_BANK_%04d"%j, seed="%d" % (j*nbanks+1), mchirp_boundaries_file=mchirp_boundaries_fname, mchirp_boundaries_index=str(j), p_node=p_node, bank_seed=bank_seed)
    bank_node.add_var_opt("match-min", mm)
    bank_node.set_priority(1)  # want complete bank before sims
    bank_nodes.append(bank_node)
    bank_name, = bank_node.get_output_files()
    bank_names.append(bank_name)

# recombine bank fragments under a common name
lwaddJob = LWAddJob(cp, tag_base=options.user_tag + "_lwadd")
lwaddNode = LWAddNode(lwaddJob, dag, bank_names, "SBANK_COMBINED-%s.xml.gz" % options.user_tag, bank_nodes)

# write the dag
dag.write_sub_files()
dag.write_script()
dag.write_cache()
dag.write_dag()
