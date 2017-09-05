# Copyright (C) 2013 Ruslan Vaulin
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

## \defgroup laldetchar_py_idq_auxmvc AuxMVC Module
## \ingroup laldetchar_py_idq
## Synopsis
# ~~~
# from laldetchar.idq import auxmvc
# ~~~
# \author Ruslan Vaulin (<ruslan.vaulin@ligo.org>)

"""Module to keep definitions of condor job and node classes for auxmvc algorithms.
"""

import os
import sys
import tempfile
from glue import pipeline
import tempfile
from laldetchar import git_version

__author__ = 'Ruslan Vaulin <ruslan.vaulin@ligo.org>'
__version__ = git_version.id
__date__ = git_version.date


## \addtogroup laldetchar_py_idq_auxmvc
# @{

# definitions of condor job and node classes for auxmvc algroithms.

def construct_command(node):
    command_string = node.job().get_executable() + ' ' \
        + node.get_cmd_line()
    return [opt for opt in command_string.split(' ') if opt.strip()]


class auxmvc_DAG(pipeline.CondorDAG):

    def __init__(self, basename, log_path):
        self.basename = basename
        tempfile.tempdir = log_path
        tempfile.template = self.basename + '.dag.log.'
        logfile = tempfile.mktemp()
        fh = open(logfile, 'w')
        fh.close()
        pipeline.CondorDAG.__init__(self, logfile)
        self.set_dag_file(self.basename)
        self.jobsDict = {}

    # self.id = 0

    def add_node(self, node):

    # self.id+=1

        pipeline.CondorDAG.add_node(self, node)


#####################  JOB and NODE classes for auxmvc pipeline  #################################

class auxmvc_analysis_job(pipeline.AnalysisJob, pipeline.CondorDAGJob):

    """
  A basic auxmvc job class. Sets common atributes needed for any auxmvc job. It uses config parser object to 
  set the options. 
  """

    def __init__(
        self,
        cp,
        sections,
        exec_name,
        tag_base='',
        id='',
        extension='',
        dax=False,
        short_opts=False,
        ):
        """
    cp = ConfigParser object from which options are read.
    sections = sections of the ConfigParser that get added to the opts
    exec_name = exec_name name in ConfigParser
    """

        self.__exec_name = exec_name
        self.__extension = extension
        self.tag_base = tag_base
        universe = cp.get('condor', 'universe')
        executable = cp.get('condor', exec_name)
        pipeline.CondorDAGJob.__init__(self, universe, executable)
        pipeline.AnalysisJob.__init__(self, cp, dax)
        self.add_condor_cmd('copy_to_spool', 'False')
        self.add_condor_cmd('getenv', 'True')
        self.add_condor_cmd('environment',
                            'KMP_LIBRARY=serial;MKL_SERIAL=yes')
        self.__use_gpus = cp.has_option('condor', 'use-gpus')

        ### accounting tags
        if cp.has_option('condor', 'accounting_group'):
            self.add_condor_cmd('accounting_group', cp.get('condor', 'accounting_group') )

        if cp.has_option('condor', 'accounting_group_user'):
            self.add_condor_cmd('accounting_group_user', cp.get('condor', 'accounting_group_user') )

        for sec in sections:
            if cp.has_section(sec):
                if short_opts:
                    self.add_short_ini_opts(cp, sec)
                else:
                    self.add_ini_opts(cp, sec)
            else:
                print >> sys.stderr, \
                    'warning: config file is missing section [' + sec \
                    + ']'

        self.set_stdout_file('logs/' + tag_base + id + '.out')
        self.set_stderr_file('logs/' + tag_base + id + '.err')
        self.set_sub_file(tag_base + '.sub')

    def set_exec_name(self, exec_name):
        """
    Set the exec_name name
    """

        self.__exec_name = exec_name

    def set_exec_name(self, exec_name):
        """
    Set the exec_name name
    """

        self.__exec_name = exec_name

    def get_exec_name(self):
        """
    Get the exec_name name
    """

        return self.__exec_name

    def set_extension(self, extension):
        """
    Set the file extension
    """

        self.__extension = extension

    def get_extension(self):
        """
    Get the extension for the file name
    """

        return self.__extension

    def get_use_gpus(self):
        """
    Get whether this job was requested to run on a GPU node
    """

        return self.__use_gpus

    def add_short_ini_opts(self, cp, section):
        """
    Parse command line options from a given section in an ini file and
    pass to the executable as short options.
    @param cp: ConfigParser object pointing to the ini file.
    @param section: section of the ini file to add to the options.
    """

        for opt in cp.options(section):
            arg = cp.get(section, opt)
            self.add_short_opt(opt, arg)


class build_auxmvc_vectors_job(auxmvc_analysis_job):

    """
  Job for building auxmvc feature vectors. 
  """

    def __init__(
        self,
        cp,
        main_channel,
        channels=None,
        unsafe_channels=None,
        ):
        """
    """

        sections = ['build_auxmvc_vectors']
        exec_name = 'idq_build_auxmvc_vectors'
        tag_base = 'build_vectors'
        auxmvc_analysis_job.__init__(self, cp, sections, exec_name,
                tag_base=tag_base)
        self.add_opt('main-channel', main_channel)
        if channels:
            self.add_opt('channels', channels)
        if unsafe_channels:
            self.add_opt('unsafe-channels', unsafe_channels)


class build_auxmvc_vectors_node(pipeline.CondorDAGNode):

    """
  Dag node for building auxmvc feature vector.
  """

    def __init__(
        self,
        job,
        trigdir,
        gps_start_time,
        gps_end_time,
        output_file,
        dq_segments=None,
        dq_segments_name='',
        p_node=[],
        ):
        job.set_stdout_file('logs/' + job.tag_base + '-'
                            + str(gps_start_time) + '-'
                            + str(gps_end_time - gps_start_time)
                            + '.out')
        job.set_stderr_file('logs/' + job.tag_base + '-'
                            + str(gps_start_time) + '-'
                            + str(gps_end_time - gps_start_time)
                            + '.err')
        pipeline.CondorDAGNode.__init__(self, job)
        self.add_output_file(output_file)
        self.add_var_opt('trigger-dir', trigdir)
        self.add_var_opt('gps-start-time', gps_start_time)
        self.add_var_opt('gps-end-time', gps_end_time)
        self.add_var_opt('output-file', output_file)
        if dq_segments:
            self.add_var_opt('dq-segments', dq_segments)
            self.add_var_opt('dq-segments-name', dq_segments_name)
        for p in p_node:
            self.add_parent(p)


class prepare_training_auxmvc_samples_job(auxmvc_analysis_job):

    """
  Job for preparing training auxmvc samples. 
  """

    def __init__(self, cp):
        """
    """

        sections = ['prepare_training_auxmvc_samples']
        exec_name = 'idq_prepare_training_auxmvc_samples'
        tag_base = 'training_auxmvc'
        auxmvc_analysis_job.__init__(self, cp, sections, exec_name,
                tag_base=tag_base)


class prepare_training_auxmvc_samples_node(pipeline.CondorDAGNode):

    """
  Node for preparing training auxmvc samples job. 
  """

    def __init__(
        self,
        job,
        source_dir,
        gps_start_time,
        gps_end_time,
        output_file,
        dq_segments='',
        dq_segments_name='',
        p_node=[],
        ):
        job.set_stdout_file('logs/'
                            + os.path.split(output_file)[1].replace('.pat'
                            , '.out'))
        job.set_stderr_file('logs/'
                            + os.path.split(output_file)[1].replace('.pat'
                            , '.err'))
        pipeline.CondorDAGNode.__init__(self, job)
        self.add_output_file(output_file)
        self.add_var_opt('source-directory', source_dir)
        self.add_var_opt('gps-start-time', gps_start_time)
        self.add_var_opt('gps-end-time', gps_end_time)
        if dq_segments and dq_segments_name:
            self.add_var_opt('dq-segments', dq_segments)
            self.add_var_opt('dq-segments-name', dq_segments_name)
        self.add_var_opt('output-file', output_file)
        for p in p_node:
            self.add_parent(p)


class add_file_to_cache_job(auxmvc_analysis_job):

    """
  Job for preparing training auxmvc samples. 
  """

    def __init__(self, cp):
        """
    """

#        sections = ['add_file_to_cache']
        sections = [] ### never actually used, only references "condor" section
        exec_name = 'add_file_to_cache'
        tag_base = 'add_file_to_cache'
        auxmvc_analysis_job.__init__(self, cp, sections, exec_name,
                tag_base=tag_base)


class add_file_to_cache_node(pipeline.CondorDAGNode):

    """
  Node for preparing training auxmvc samples job. 
  """

    def __init__(
        self,
        job,
        files,
        cache,
        p_node=[],
        ):
        job.set_stdout_file('logs/'
                            + os.path.split(files[0])[1].split('.')[0]
                            + '_adding_to_cache.out')
        job.set_stderr_file('logs/'
                            + os.path.split(files[0])[1].split('.')[0]
                            + '_adding_to_cache.err')
        pipeline.CondorDAGNode.__init__(self, job)
        self.add_var_arg(cache)
        for file in files:
            self.add_var_arg(file)
        for p in p_node:
            self.add_parent(p)


class train_forest_job(auxmvc_analysis_job):

    """
  Training job for random forest (MVSC). 
  """

    def __init__(self, cp):
        """
    """

        sections = ['train_forest']
        exec_name = 'SprBaggerDecisionTreeApp'
        tag_base = 'train_forest'
        auxmvc_analysis_job.__init__(
            self,
            cp,
            sections,
            exec_name,
            tag_base=tag_base,
            short_opts=True,
            )


class train_forest_node(pipeline.CondorDAGNode):

    """
  Dag node for training the random forest (MVSC).
  """

    def __init__(
        self,
        job,
        training_data_file,
        trainedforest_filename,
        p_node=[],
        ):
        job.set_stdout_file('logs/'
                            + os.path.split(training_data_file)[1].replace('.pat'
                            , '.out'))
        job.set_stderr_file('logs/'
                            + os.path.split(training_data_file)[1].replace('.pat'
                            , '.err'))
        pipeline.CondorDAGNode.__init__(self, job)
        self.add_input_file(training_data_file)
        self.training_data_file = self.get_input_files()[0]
        self.trainedforest = trainedforest_filename
        self.add_output_file(self.trainedforest)
        self.add_var_opt('f', self.trainedforest, short=True)
        self.add_file_arg(' %s' % self.training_data_file)
        for p in p_node:
            self.add_parent(p)


class use_forest_job(auxmvc_analysis_job):

    """
  Job using random forest to evaluate unclassified data.
  """

    def __init__(self, cp):
        """
    """

        sections = ['forest_evaluate']
        exec_name = 'SprOutputWriterApp'
        tag_base = 'forest_evaluate'
        auxmvc_analysis_job.__init__(
            self,
            cp,
            sections,
            exec_name,
            tag_base=tag_base,
            short_opts=True,
            )
        self.add_short_opt('A', '')


class use_forest_node(pipeline.CondorDAGNode):

    """
  Node for radnom forest evaluation job. 
  """

    def __init__(
        self,
        job,
        trainedforest,
        file_to_rank,
        ranked_file,
        p_node=[],
        ):
        job.set_stdout_file('logs/'
                            + os.path.split(ranked_file)[1].replace('.dat'
                            , '.out'))
        job.set_stderr_file('logs/'
                            + os.path.split(ranked_file)[1].replace('.dat'
                            , '.err'))
        pipeline.CondorDAGNode.__init__(self, job)
        self.add_input_file(trainedforest)
        self.add_input_file(file_to_rank)
        self.add_output_file(ranked_file)
        self.trainedforest = self.get_input_files()[0]
        self.file_to_rank = self.get_input_files()[1]
        self.ranked_file = ranked_file
        self.add_file_arg(' %s %s %s' % (self.trainedforest,
                          self.file_to_rank, self.ranked_file))
        for p in p_node:
            self.add_parent(p)


class forest_add_excluded_vars_job(auxmvc_analysis_job):

    """
  A simple fix job that adds the variables excluded by forest (MVSC) from classification into the output file.
  Need to be run right after use_forest job. 
  """

    def __init__(self, cp):
        """
    """

#        sections = ['forest_add_excluded_vars']
        sections = [] ### not used...
        exec_name = 'forest_add_excluded_vars'
        tag_base = 'forest_add_excluded_vars'
        auxmvc_analysis_job.__init__(self, cp, sections, exec_name, tag_base=tag_base)
        self.add_opt('excluded-variables', cp.get('forest_evaluate', 'z'))


class forest_add_excluded_vars_node(pipeline.CondorDAGNode):

    """
  Node for forest_add_excluded_vars_job.
  """

    def __init__(
        self,
        job,
        patfile,
        datfile,
        p_node=[],
        ):
        job.set_stdout_file('logs/'
                            + os.path.split(datfile)[1].replace('.dat',
                            'faev.out'))
        job.set_stderr_file('logs/'
                            + os.path.split(datfile)[1].replace('.dat',
                            'faev.err'))
        pipeline.CondorDAGNode.__init__(self, job)
        self.add_input_file(patfile)
        self.add_input_file(datfile)
        self.add_var_opt('pat-file', patfile)
        self.add_var_opt('dat-file', datfile)
        for p in p_node:
            self.add_parent(p)


class plot_channels_significance_job(pipeline.CondorDAGJob):

    """   
  Job that makes verious plots and histograms using significance of the axuiloary channels.      
  """

    def __init__(self, cp):
        sections = ['plot-forest-channels-significance']
        exec_name = 'auxmvc_plot_mvsc_channels_significance'
        tag_base = 'plot_channels_signif'
        auxmvc_analysis_job.__init__(self, cp, sections, exec_name,
                tag_base=tag_base)


class plot_channels_significance_node(pipeline.CondorDAGNode):

    """
  Node for plot_channels_significance_job.
  """

    def __init__(
        self,
        job,
        input,
        p_node=[],
        ):
        pipeline.CondorDAGNode.__init__(self, job)
        self.add_var_opt('input', input)
        for p in p_node:
            self.add_parent(p)


class result_plots_job(pipeline.CondorDAGJob):

    """
  Job that makes plots based on results of evaluation e.g. ROC curves.
  """

    def __init__(self, cp, tag_base='RESULT_PLOTS'):
        """
    """

        sections = ['result_plots']
        exec_name = 'auxmvc_result_plots'
        tag_base = 'auxmvc_result_plots'
        auxmvc_analysis_job.__init__(self, cp, sections, exec_name,
                tag_base=tag_base)


class result_plots_node(pipeline.CondorDAGNode):

    """
  Node for result_plots_job.
  """

    def __init__(
        self,
        job,
        datfiles,
        p_node=[],
        ):
        pipeline.CondorDAGNode.__init__(self, job)
        for file in datfiles:
            self.add_file_arg(file[0])
        for p in p_node:
            self.add_parent(p)


########################  svm for idq  ############################

class use_svm_job(auxmvc_analysis_job):

    """
    """

    def __init__(self, cp):
        """
        """

        sections = ['svm_evaluate']
        exec_name = 'svm_evaluate_cmd'
        tag_base = 'svm_evaluate'
        auxmvc_analysis_job.__init__(self, cp, sections, exec_name,
                tag_base=tag_base)


class use_svm_node(pipeline.CondorDAGNode):

    """
    Node for SVM evaluation job.
    """

    def __init__(
        self,
        job,
        cp,
        test_file,
        range_file,
        svm_model,
        predict_file,
        p_node=[],
        ):
        job.set_stdout_file('logs/'
                            + os.path.split(predict_file)[1].replace('.dat'
                            , '.out'))
        job.set_stderr_file('logs/'
                            + os.path.split(predict_file)[1].replace('.dat'
                            , '.err'))
        pipeline.CondorDAGNode.__init__(self, job)
        self.add_input_file(test_file)
        self.add_input_file(range_file)
        self.add_input_file(svm_model)
        self.add_output_file(predict_file)

        # self.scale_cmd = cp.get('svm_evaluate','svm_scale_cmd')
        # self.predict_cmd = cp.get('svm_evaluate', 'svm_predict_cmd')

        self.test_file = self.get_input_files()[0]
        self.range_file = self.get_input_files()[1]
        self.svm_model = self.get_input_files()[2]
        self.predict_file = self.get_output_files()[0]
        self.add_var_opt('i', self.test_file, short=True)
        self.add_var_opt('r', self.range_file, short=True)
        self.add_var_opt('m', self.svm_model, short=True)
        self.add_var_opt('o', self.predict_file, short=True)

        # self.add_file_arg(" --scale %s --predict %s -i %s -r %s -m %s -o %s" % (self.scale_cmd, self.predict_cmd, self.test_file, self.range_file, self.svm_model, self.predict_file))

        for p in p_node:
            self.add_parent(p)


class train_svm_job(auxmvc_analysis_job):

    """
  Training job for svm.
  """

    def __init__(self, cp):
        """
    """

        sections = ['train_svm']  # no section in configuration yet
        exec_name = 'svm_train_cmd'
        tag_base = 'train_svm'
        auxmvc_analysis_job.__init__(self, cp, sections, exec_name,
                tag_base=tag_base)


class train_svm_node(pipeline.CondorDAGNode):

    """
    Node for SVM train job.
    """

    def __init__(
        self,
        job,
        cp,
        train_file,
        range_file,
        model_file,
        p_node=[],
        ):
        job.set_stdout_file('logs/'
                            + os.path.split(train_file)[1].replace('.pat'
                            , '.out'))
        job.set_stderr_file('logs/'
                            + os.path.split(train_file)[1].replace('.pat'
                            , '.err'))
        pipeline.CondorDAGNode.__init__(self, job)
        self.add_input_file(train_file)
        self.add_output_file(range_file)
        self.add_output_file(model_file)

        # self.scale_cmd = cp.get('svm_evaluate','svm_scale_cmd')
        # self.train_cmd = cp.get('svm_evaluate','svm_train_cmd')
        # self.gamma = cp.get('svm_evaluate','svm_gamma')
        # self.cost  = cp.get('svm_evaluate','svm_cost')

        self.train_file = self.get_input_files()[0]
        self.range_file = self.get_output_files()[0]
        self.model_file = self.get_output_files()[1]

        self.train_file_svm = os.path.abspath(self.train_file) + '.mid'
        self.scale_file = os.path.abspath(self.train_file) + '.scale'
        self.add_var_opt('train-file', self.train_file)
        self.add_var_opt('train-file-svm', self.train_file_svm)
        self.add_var_opt('scale-file', self.scale_file)
        self.add_var_opt('range-file', self.range_file)
        self.add_var_opt('model-file', self.model_file)

        # self.add_file_arg(" --scale %s --train %s -g %s -c %s " % (self.scale_cmd, self.train_cmd, self.train_file, self.train_file_svm, self.scale_file, self.range_file, self.model_file, self.gamma, self.cost))

        self.set_post_script('/bin/rm ')
        self.add_post_script_arg(self.train_file_svm)
        self.add_post_script_arg(self.scale_file)

        for p in p_node:
            self.add_parent(p)


########################  ann for idq  ############################

class convert_annfile_job(auxmvc_analysis_job):

    """
  Job for converting pat files to FANN type files.
  """

    def __init__(self, cp):
        """
    """

        sections = ['ann_convert']
        exec_name = 'ConvertANNData'
        tag_base = 'ann_convert'
        auxmvc_analysis_job.__init__(
            self,
            cp,
            sections,
            exec_name,
            tag_base=tag_base,
            short_opts=False,
            )


class convert_annfile_node(pipeline.CondorDAGNode):

    """
  Dag node for converting pat files to FANN type files.
  """

    def __init__(
        self,
        job,
        training_data_file,
        p_node=[],
        ):
        job.set_stdout_file('logs/'
                            + os.path.split(training_data_file)[1].replace('.pat'
                            , '.out'))
        job.set_stderr_file('logs/'
                            + os.path.split(training_data_file)[1].replace('.pat'
                            , '.err'))
        pipeline.CondorDAGNode.__init__(self, job)
        self.add_input_file(training_data_file)
        self.training_data_file = self.get_input_files()[0]
        self.fann_file = training_data_file.replace('.pat', '.ann')
        self.add_output_file(self.fann_file)
        self.add_file_arg(' %s' % self.training_data_file)
        for p in p_node:
            self.add_parent(p)


class train_ann_job(auxmvc_analysis_job):

    """
  Training job for Artifical Neural Networks (ANN) with iRPROP- algorithm.
  """

    def __init__(self, cp):
        """
    """

        sections = ['train_ann']
        exec_name = 'TrainNeuralNet'
        tag_base = 'train_ann'
        auxmvc_analysis_job.__init__(
            self,
            cp,
            sections,
            exec_name,
            tag_base=tag_base,
            short_opts=False,
            )


class train_ann_node(pipeline.CondorDAGNode):

    """
  Dag node for training the Artificial Neural Networks (ANN) with iRPROP- algorithm.
  """

    def __init__(
        self,
        job,
        training_data_file,
        trained_ann_filename,
        p_node=[],
        ):
        job.set_stdout_file('logs/'
                            + os.path.split(training_data_file)[1].replace('.ann'
                            , '.out'))
        job.set_stderr_file('logs/'
                            + os.path.split(training_data_file)[1].replace('.ann'
                            , '.err'))
        pipeline.CondorDAGNode.__init__(self, job)
        self.add_input_file(training_data_file)
        self.training_data_file = self.get_input_files()[0]
        self.trained_ann = trained_ann_filename
        self.add_output_file(self.trained_ann)
        self.add_file_arg(' -t %s -s %s' % (self.training_data_file,
                          self.trained_ann))
        for p in p_node:
            self.add_parent(p)


class use_ann_job(auxmvc_analysis_job):

    """
  Job using ANN to evaluate unclassified data.
  """

    def __init__(self, cp):
        """
    """

        sections = ['ann_evaluate']
        exec_name = 'EvaluateNeuralNet'
        tag_base = 'ann_evaluate'
        auxmvc_analysis_job.__init__(
            self,
            cp,
            sections,
            exec_name,
            tag_base=tag_base,
            short_opts=False,
            )
        self.add_short_opt('', '')


class use_ann_node(pipeline.CondorDAGNode):

    """
  Node for ANN evaluation job. 
  """

    def __init__(
        self,
        job,
        trained_ann,
        file_to_rank,
        ranked_file,
        p_node=[],
        ):
        job.set_stdout_file('logs/'
                            + os.path.split(ranked_file)[1].replace('.dat'
                            , '.out'))
        job.set_stderr_file('logs/'
                            + os.path.split(ranked_file)[1].replace('.dat'
                            , '.err'))
        pipeline.CondorDAGNode.__init__(self, job)
        self.add_input_file(trained_ann)
        self.add_input_file(file_to_rank)
        self.add_output_file(ranked_file)
        self.trained_ann = self.get_input_files()[0]
        self.file_to_rank = self.get_input_files()[1]
        self.ranked_file = ranked_file
        self.add_file_arg(' -n %s -e %s -s %s' % (self.trained_ann,
                          self.file_to_rank, self.ranked_file))
        for p in p_node:
            self.add_parent(p)


##@}
