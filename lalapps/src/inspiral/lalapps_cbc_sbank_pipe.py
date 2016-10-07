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
        self.add_condor_cmd('request_memory', '1999')
        if "OMP_NUM_THREADS" in os.environ:
            self.add_condor_cmd('request_cpus', os.environ["OMP_NUM_THREADS"])


class SBankNode(pipeline.CondorDAGNode):
    def __init__(self, job, dag, tag=None, seed=0, bank_seed=None, mchirp_boundaries_file=None, mchirp_boundaries_index=None, p_node=[]):
        pipeline.CondorDAGNode.__init__(self,job)
        if (mchirp_boundaries_file is None) ^ (mchirp_boundaries_index is None):
            raise ValueError("must supply both mchirp_boundaries_file and mchirp_boundaries_index or neither")
        if mchirp_boundaries_file is not None:
            self.add_var_arg("--mchirp-boundaries-file %s --mchirp-boundaries-index %s" % (mchirp_boundaries_file, mchirp_boundaries_index))
        if bank_seed is not None:
            self.add_var_arg("--bank-seed %s" % bank_seed)
        self.add_var_opt("seed", seed)

        if tag:
            self.add_var_opt("user-tag", tag)
            self.add_output_file("%s-SBANK_%s-%d-%d.xml.gz" % (job.get_opt("instrument"), tag, int(job.get_opt("gps-start-time")), int(job.get_opt("gps-end-time")) - int(job.get_opt("gps-start-time"))))
        else:
            self.add_output_file("%s-SBANK-%d-%d.xml.gz" % (job.get_opt("instrument"), int(job.get_opt("gps-start-time")), int(job.get_opt("gps-end-time")) - int(job.get_opt("gps-start-time"))))

        for p in p_node:
            self.add_parent(p)

        dag.add_node(self)


class SBankSplitJob(inspiral.InspiralAnalysisJob):
    def __init__(self,cp,dax=False, tag_base="ligolw_cbc_sbank_splitter"):
        exec_name = 'lalapps_cbc_sbank_splitter'
        extension = 'xml'
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


class SBankSplitNode(pipeline.CondorDAGNode):
    def __init__(self, job, dag, bank, tag = None, p_node=[]):
        pipeline.CondorDAGNode.__init__(self,job)
        self.add_file_arg(bank)
        if tag:
            self.add_var_opt("user-tag", tag)
        for p in p_node:
            self.add_parent(p)
        nbanks = int(self.job().get_opts()["nbanks"])
        for i in xrange(nbanks):
            self.add_output_file("H1-SBANK_SPLIT_%04d-%s.xml" % (i+1, tag))
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

class InspinjJob(inspiral.InspiralAnalysisJob):
    def __init__(self,cp,dax=False, tag_base="inspinj"):
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
    def __init__(self,job, dag, tag=None,seed=0, p_node=[]):
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
    def __init__(self,cp,dax=False, tag_base ="banksim"):
        exec_name = 'lalapps_cbc_sbank_sim'
        extension = 'xml'
        sections = ['banksim']
        inspiral.InspiralAnalysisJob.__init__(self,cp,sections,exec_name,extension,dax)
        self.add_condor_cmd('request_memory', '1999')
        self.tag_base = tag_base
        self.add_condor_cmd('getenv','True')
        self.set_stdout_file('logs/'+tag_base+'-$(macroid)-$(process).out')
        self.set_stderr_file('logs/'+tag_base+'-$(macroid)-$(process).err')
        if cp.has_section("accounting"):
            self.add_condor_cmd('accounting_group', cp.get("accounting", "accounting-group"))
        if "OMP_NUM_THREADS" in os.environ:
            self.add_condor_cmd('request_cpus', os.environ["OMP_NUM_THREADS"])


class BankSimNode(pipeline.CondorDAGNode):
    def __init__(self, job, dag, injs, approx, bank, tag, p_node=[]):
        pipeline.CondorDAGNode.__init__(self,job)
        self.add_file_opt("injection-file", injs)
        self.add_var_opt("injection-approx", approx)
        self.add_file_opt("template-bank", bank)
        self.add_var_opt("user-tag", tag)
        self.add_output_file("%s.h5" % tag)
        for p in p_node:
            if p is not None:
                self.add_parent(p)
        dag.add_node(self)


class LWAddJob(pipeline.CondorDAGJob):
    """
    A ligolw_add node. This node is used to combine the split template banks
    into one aggregate bank.
    """
    def __init__(self, executable=which('ligolw_add'), tag_base='ligolw_add'):
        """
        """
        self.__prog__ = 'ligolw_add'
        self.__executable = executable
        self.__universe = 'vanilla'
        pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
        self.add_condor_cmd('getenv','True')
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


class MergeSimsJob(inspiral.InspiralAnalysisJob):
    def __init__(self,cp,dax=False, tag_base="lalapps_cbc_sbank_merge_sims"):
        exec_name = 'lalapps_cbc_sbank_merge_sims'
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
        self.add_condor_cmd('request_memory', '1999')


class MergeSimsNode(pipeline.CondorDAGNode):
    def __init__(self,job, dag, tag='sbank_merge_sims', h5files=[], p_node=[]):
        pipeline.CondorDAGNode.__init__(self,job)
        output_fname = "SBANK_MERGED_SIMS_%s.h5" % tag
        self.add_var_opt("output", output_fname)
        for f in h5files:
            self.add_file_arg(f)
        for p in p_node:
            self.add_parent(p)
        self.add_output_file(output_fname)
        dag.add_node(self)


class PlotSimJob(inspiral.InspiralAnalysisJob):
    def __init__(self,cp,dax=False, tag_base="lalapps_cbc_sbank_plot_sim"):
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
        self.add_condor_cmd('request_memory', '1999')


class PlotSimNode(pipeline.CondorDAGNode):
    def __init__(self,job, dag, h5files=[],p_node=[]):
        pipeline.CondorDAGNode.__init__(self,job)
        if len(h5files) != 1:
            raise ValueError("Can only plot one h5 file at a time")
        self.add_file_arg(h5files[0])
        for p in p_node:
            self.add_parent(p)
        dag.add_node(self)


usage = """

lalapps_cbc_sbank_pipe generates a DAG for parallelizing the
construction and verification of template banks via lalapps_cbc_sbank
and lalapps_cbc_sbank_sim. Its primary input is an ini file, which
takes form given below. Once you have an ini file, create the DAG by

lalapps_cbc_sbank_pipe --config-file INIFILE.ini --user-tag SBANK_PIPE

If already have a bank, but want to patch holes, or extend the
parameter space it covers, you can provide the option --bank-seed
BANK_SEED.xml. The bank BANK_SEED.xml can be generated by any means,
and the resulting DAG will generate a bank with that bank as a subset.

[FIXME: In bank-seed mode, the *banksim* portion is not parallelized
along the bank seed, so if your bank seed is already most of the final
bank, your bank sim will spend most of its time on the input bank.]

You can also add --template-bank BANK.xml to run only the sims. To do
so, you have to remove the "template-weight" option from the [split]
section.

;;; BEGIN EXAMPLE INI ;;;

[condor] ; This section points to the executables, and provides condor options
universe = vanilla
lalapps_cbc_sbank = /home/sprivite/local/opt/master/bin/lalapps_cbc_sbank
lalapps_cbc_sbank_splitter = /home/sprivite/local/opt/master/bin/lalapps_cbc_sbank_splitter
lalapps_cbc_sbank_choose_mchirp_boundaries = /home/sprivite/local/opt/master/bin/lalapps_cbc_sbank_choose_mchirp_boundaries
lalapps_cbc_sbank_sim = /home/sprivite/local/opt/master/bin/lalapps_cbc_sbank_sim
lalapps_cbc_sbank_plot_sim = /home/sprivite/local/opt/master/bin/lalapps_cbc_sbank_plot_sim
lalapps_cbc_sbank_merge_sims = /home/sprivite/local/opt/master/bin/lalapps_cbc_sbank_merge_sims
lalapps_inspinj = /home/sprivite/local/opt/master/bin/lalapps_inspinj

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
template-weight = equal

;
; FOR BANK SIMS ONLY
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
gps-start-time = 1000000000
gps-end-time =  1000010000
time-step = 1

[seobnr]
; This section overrides the generic inspinj options from the section
; above. If an option is not included here, it will be taken from the
; [inspinj] section. The options in this section apply only to the
; [seobnr] run, in this example.
waveform = SEOBNRv2_ROM_DoubleSpin
enable-spin =
min-spin1 = 0.0
max-spin1 = 0.98
min-spin2 = 0.0
max-spin2 = 0.98
aligned =

[imrc]
; This section overrides the generic inspinj options from the section
; above. If an option is not included here, it will be taken from the
; [inspinj] section. The options in this section apply only to the
; [imrc] run, in this example.
waveform = SEOBNRv2_ROM_DoubleSpin
enable-spin =
min-spin1 = 0.0
max-spin1 = 0.9
min-spin2 = 0.0
max-spin2 = 0.9
aligned =

[banksim]
; In this section, you can globally override some inspinj options, and
; reset some of the template bank options. For example, you may want
; to generate the template bank using the TaylorF2RedSpin metric, but
; run the bank sim using IMRPhenomC templates. See
; lalapps_cbc_sbank_sim for available options. Any option set here
; applies to every injection run.
noise-model = aLIGOZeroDetHighPower
mtotal-max = 350
mtotal-min = 4
mratio-max = 14

;;; END EXAMPLE INI ;;;

"""

def parse_command_line():
    parser = OptionParser(usage=usage)
    parser.add_option("--config-file", default=None, help="Read options for generating template bank placement pipeline from configuration ini file..")
    parser.add_option("--user-tag", default="SBANK", help="Make your results feel special and give them a unique name.")
    parser.add_option("--template-bank", default=None, help="Skip the template bank generation stage. Use given template bank for bank sim.")
    parser.add_option("--bank-seed", default=None, help="Use the provided template bank as a start bank. sBank will fill in any holes in this bank.")
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

# using the bank-seed option in the ini file can only ever lead to
# problems, so don't let the user do it
if cp.has_option("sbank", "bank-seed"):
    raise ValueError("You have indicated a bank seed in the sbank section of the ini file. Please specify instead via the command line option --bank-seed.")

# set up bank generation
# Two modes:
#   1. generate coarse, partition mchirp space, generate sub-banks, ligolw_add
#   2. take given bank, split it, and sim it
nbanks = int(cp.get("split","nbanks"))
bank_names = []
bank_nodes = []
if options.template_bank:
    # skip planning stage if bank provided :)
    xmlCoarse = os.path.basename(options.template_bank)
    if not os.path.isfile(xmlCoarse):
        shutil.copy(options.template_bank, xmlCoarse)

    cp.set("split", "instrument", cp.get("sbank", "instrument")) # these parameters need to agree or the DAG will crash
    sbankSplitJob = SBankSplitJob(cp)
    sbankSplitNode = SBankSplitNode(sbankSplitJob, dag, xmlCoarse, options.user_tag)

    bank_names = sbankSplitNode.get_output_files()
    bank_nodes = [sbankSplitNode]
else:
    if options.bank_seed:
        xmlCoarse = options.bank_seed
        pnode = []
        bank_names.append(options.bank_seed)
        bank_nodes.append(None)
    else:
        # set up sole coarse node to plan out the mini-sbank nodes
        coarse_sbank_node = SBankNode(sbankJob, dag, "COARSE")
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
    for j in xrange(nbanks):
        bank_node = SBankNode(sbankJob, dag, "%04d"%j, seed="%d" % (j*nbanks+1), mchirp_boundaries_file=mchirp_boundaries_fname, mchirp_boundaries_index=str(j), p_node=[sbankChooseMchirpBoundariesNode], bank_seed=xmlCoarse)
        bank_node.add_var_opt("match-min", mm)
        bank_node.set_priority(1)  # want complete bank before sims
        bank_nodes.append(bank_node)
        bank_name, = bank_node.get_output_files()
        bank_names.append(bank_name)

# recombine bank fragments under a common name
if not options.template_bank:
    lwaddJob = LWAddJob(tag_base=options.user_tag + "_lwadd")
    lwaddNode = LWAddNode(lwaddJob, dag, bank_names, "H1-SBANK_COMBINED-%s.xml.gz"%options.user_tag,bank_nodes)

# set up banksim parameters according to sbank parameters (if not provided)
if not cp.has_option("banksim","flow"):
    cp.set("banksim","flow",cp.get("sbank","flow"))

# you could want change the template waveform I guess...
if not cp.has_option("banksim","template-approx"):
    cp.set("banksim","template-approx",cp.get("sbank","approximant"))

# copy over PSD if not specified
if not cp.has_option("banksim","reference-psd") and not cp.has_option("banksim","noise-model"):
    if cp.has_option("sbank","reference-psd"):
        cp.set("banksim","reference-psd",cp.get("sbank","reference-psd"))
        cp.set("banksim","instrument",cp.get("sbank","instrument"))
    else:
        cp.set("banksim","noise-model",cp.get("sbank","noise-model"))

# sim nodes
banksimJob = BankSimJob(cp)

# add sim fragments together nodes
merge_simsJob = MergeSimsJob(cp)
plotsimJob = PlotSimJob(cp)

# for each injection run and every subbank, perform a banksim
for inj_run in cp.options("injections"):
    # set up inspinj jobs
    # each option in the [injections] section is a separate injection run
    # the inspinj section sets options that are common to all injection runs
    for opt in cp.options(inj_run):
        # override default inspinj options with those given in the inj_run section
        cp.set('inspinj',opt,cp.get(inj_run,opt))

    # set random seed
    seed = int(cp.get("injections",inj_run))
    cp.set('inspinj','seed',str(seed))
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
    for bank_node in bank_nodes:
        for bank_name in bank_node.get_output_files():
            ind = bank_name.index("-", 3)  # start searching for hyphens after ind 3
            base, _ = os.path.splitext(bank_name)
            sim_name = "%s_%s-%s" % (base[:ind], "SIM_%s"%inj_run.upper(), base[ind + 1:])
            sim_nodes.append(BankSimNode(banksimJob, dag, inj_name, waveform, bank_name, sim_name.replace(".h5",""), [inj_node, bank_node]))

    # merge and plot the partial sims
    sim_names = [node.get_output_files()[0] for node in sim_nodes]
    merge_sims_node = MergeSimsNode(merge_simsJob, dag, tag, sim_names, sim_nodes )
    PlotSimNode(plotsimJob, dag, merge_sims_node.get_output_files(), [merge_sims_node])

# write the dag
dag.write_sub_files()
dag.write_script()
dag.write_cache()
dag.write_dag()
