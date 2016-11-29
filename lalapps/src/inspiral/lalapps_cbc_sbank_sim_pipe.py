# Copyright (C) 2015  Surabhi Sachdev, Tjonnie Li, Stephen Privitera
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

import sys, os, shutil
from itertools import cycle, izip
import subprocess, socket, tempfile
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

class InjSplitJob(inspiral.InspiralAnalysisJob):
    def __init__(self,cp,dax=False, tag_base ="injsplitter"):
        exec_name = 'lalapps_cbc_sbank_injsplitter'
        extension = 'xml'
        sections = ['injsplitter']
        inspiral.InspiralAnalysisJob.__init__(self,cp,sections,exec_name,extension,dax)
        self.add_condor_cmd('request_memory', '1999')
        self.tag_base = tag_base
        self.add_condor_cmd('getenv','True')
        self.set_stdout_file('logs/'+tag_base+'-$(macroid)-$(process).out')
        self.set_stderr_file('logs/'+tag_base+'-$(macroid)-$(process).err')


class InjSplitNode(pipeline.CondorDAGNode):
    def __init__(self, job, dag, injfile, nsplit, usertag, p_node=[]):
        pipeline.CondorDAGNode.__init__(self,job)
        outfiles = ["%s_INJ_SPLIT_%04d.xml"%(usertag,i) for i in xrange(nsplit)]
        self.add_var_opt("nsplit", nsplit)
        self.add_var_opt("usertag", usertag)
        self.add_file_arg(injfile)
        for f in outfiles:
            self.add_output_file(f)
        for p in p_node:
            if p is not None:
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


class BankSimNode(pipeline.CondorDAGNode):
    def __init__(self, job, dag, injs, approx, bank, tag, p_node=[]):
        pipeline.CondorDAGNode.__init__(self,job)
        self.add_file_opt("injection-file", injs)
        #self.add_var_opt("injection-approx", approx)
        self.add_file_opt("template-bank", bank)
        self.add_var_opt("user-tag", tag)
        self.add_output_file("%s.h5" % tag)
        for p in p_node:
            if p is not None:
                self.add_parent(p)
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
        self.add_condor_cmd('request_memory', '3999')


class PlotSimNode(pipeline.CondorDAGNode):
    def __init__(self,job, dag, h5files=[],p_node=[]):
        pipeline.CondorDAGNode.__init__(self,job)
        [self.add_file_arg(h5f) for h5f in h5files]
        for p in p_node:
            self.add_parent(p)
        dag.add_node(self)


###############################################################################
#
# OPTION PARSER
#
###############################################################################

def parse_command_line():
    parser = OptionParser()
    parser.add_option("--config-file", default=None, help="Read options for generating template bank placement pipeline from configuration ini file..")
    options, filenames = parser.parse_args()

    return options, filenames

options, filenames = parse_command_line()


###############################################################################
#
# CREATE DAG
#
###############################################################################

# READ IN CONFIGPARSER
cp = dcConfigParser()
cp.read(options.config_file)

# FIXME: HARDCODED FOR TESTING PURPOSES
inj_name = cp.get("run","injection");
inj_approx = cp.get("banksim",'injection-approx');
bank_name = cp.get("run",'bank');
usertag = cp.get("run",'tag');
psdfile = cp.get("banksim",'reference-psd');
nsplit = int(cp.get("injsplitter","nsplit"));

# CREATE TOP LEVEL DAG AND JOBS
dag = bank_DAG(usertag)
injsplitJob = InjSplitJob(cp);
banksimJob = BankSimJob(cp)
plotsimJob = PlotSimJob(cp)

# SETUP SPLIT INJECTION JOB
injsplitNode = InjSplitNode(injsplitJob, dag, inj_name, nsplit, usertag);

# SETUP BANKSIM NODE
sim_nodes = [BankSimNode(banksimJob, dag, f, inj_approx, bank_name, "%s_BANKSIM_SPLIT_%04d"%(usertag,i), [injsplitNode]) for i,f in zip(xrange(nsplit),injsplitNode.get_output_files())]

# SETUP MERGESIM NODE
sim_names = [node.get_output_files()[0] for node in sim_nodes]

# SETUP PLOTSIM NODE
PlotSimNode(plotsimJob, dag, sim_names, sim_nodes)

# CREATE DAG
dag.write_sub_files()
dag.write_script()
dag.write_cache()
dag.write_dag()

# CREATE LOG DIRECTORY
if not os.path.exists('logs'):
	os.makedirs('logs');
