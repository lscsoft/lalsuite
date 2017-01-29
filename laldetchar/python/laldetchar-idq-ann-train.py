# Copyright (C) 2015 Young-Min Kim
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the# Free Software Foundation; either version 3 of the License, or (at your# option
# ) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


###########################################
### This code uses FANN library 
### for Aritifical Neural Network analysis
###########################################
from optparse import *
import glob
import sys
import os
import numpy
import random
import time
###########################################
from pyfann import libfann
###########################################
import GAnnUtils
import G1DConnections
import GAnnMutators
import GAnnEvaluators
import GAnnCrossovers
import GAnnConsts
import GAnnSelectors
import GAnnInitializators
from pyevolve import GSimpleGA, Consts
import GAnnGA
###########################################
from laldetchar import git_version

__author__ = 'Young-Min Kim <young-min.kim@ligo.org>'
__version__ = git_version.id
__date__ = git_version.date

description = """This program runs ANN training tasks."""

parser = OptionParser(version='Name: %%prog\n%s'%git_version.verbose_msg, 
                                usage='%prog [options]', 
                                description=description)
parser.add_option("-v","--verbose",action="store_true",help="verbose mode for printing run processing.")
parser.add_option("-t","--training-file",action="store",type="string",help="training data file such as *. ann  with FANN type data format.")
parser.add_option("-s","--saving-results",action="store",type="string",default=False,help="saving result network")
############### network configuration parameters######################################
parser.add_option("-n","--hidden-neurons",action="store",type="string",help="number of neurons per each layer. 5,4,3 --> layer 1,2,3 have 5,4,3 neurons, respectively.First number is for input layer and last number is for output layer. ")
parser.add_option("-c","--connection-rate",action="store",type="float",default=1.0,help="connection rate. Default is 1 which gives full connections")
parser.add_option("-d","--hidden-activation",action="store", type="string", default="SIGMOID", help="Activation fuction for hidden layers. Default is SIGMOID.")
parser.add_option("-o","--output-activation",action="store", type="string", default="SIGMOID", help="Activation fuction for output layer. Default is SIGMOID.")
parser.add_option("-f","--steep-hidden",action="store",type="float",default=0.5,help="steepness of hidden layer activation function")
parser.add_option("-g","--steep-out",action="store",type="float",default=0.5,help="steepness of output layer activation function")
############### FANN parameters######################################
parser.add_option("","--desired-error",action="store",type="float",default=0.0001,help="training terminate when MSE < desired_error. Default is 0.0001")
parser.add_option("","--learning-rate",action="store",type="float",default=0.7,help="connection rate. Default is 1 which gives full connections")
parser.add_option("-m","--max-epochs",action="store",type="int",default=1000,help="max iterations")
parser.add_option("-w","--weights-min",action="store",type="float",default=-0.1,help="minimum weight")
parser.add_option("-x","--weights-max",action="store",type="float",default=0.1,help="max weight")
parser.add_option("-y","--increase-factor",action="store",type="float",default=1.2,help="rprop increase factor")
parser.add_option("-z","--decrease-factor",action="store",type="float",default=0.5,help="rprop decrease factor")
parser.add_option("-a","--delta-min",action="store",type="float",default=0.0,help="rprop delta minimum")
parser.add_option("-b","--delta-max",action="store",type="float",default=50.0,help="rprop delta maximum")
############### GA parameters######################################
parser.add_option("","--mutation-rate",action="store",type="float",default=0.2,help="Generations for GA run")
parser.add_option("","--generations",action="store",type="int",default=20,help="Generations for GA run")
parser.add_option("","--population",action="store",type="int",default=20,help="Population Size for GA run")
parser.add_option("","--range-min",action="store",type="float",default=-1.0,help="minimum weight for GA run")
parser.add_option("","--range-max",action="store",type="float",default=1.0,help="maximum weight for GA run")
parser.add_option("","--gauss-mu",action="store",type="float",default=0.0,help="mean value of Gaussian distribution for GA run")
parser.add_option("","--gauss-sigma",action="store",type="float",default=1.0,help="standard deviation of Gaussian distribution for GA run")
parser.add_option("-I","--initializator",action="store",type="string",default="uniform",help="Initializator for initial population of GA")
parser.add_option("","--import-network",action="store",type="string",default="",help="Import network from existing file")
parser.add_option("","--set-new-parameters",action="store",type="string", default="True",help="set new paramters on imported network")
parser.add_option("","--training-machine",action="store",type="string",default="ga,fann",help="use ga,fann if you want to use ga aided backpropagation algorithm for training network. If not, just use a single machine, ga or fann.")
(opts,files)=parser.parse_args()

############# MAIN ###########################
start_time=time.time()

########## Network Configuration ###########
if opts.verbose:
    print "Creating network."

########## Import training data #####################
if opts.verbose:
    print "Getting training data : %s" % opts.training_file
train_data = libfann.training_data()
train_data.read_train_from_file(opts.training_file)
train_data.shuffle_train_data()
input_number = train_data.num_input_train_data()
#train_data.scale_train_data(0.0,1.0)


# import exiting network or generate neural network
ann = libfann.neural_net()
if opts.import_network:
    if opts.verbose:
        print "%% Import form existing network."
    ann.create_from_file(opts.import_network)
    layers = ann.get_layer_array()
    bias = ann.get_bias_array()
else:
    layers = [input_number] + map(int,opts.hidden_neurons.split(",")) + [1]
    bias = [1]+[1 for i in range(len(layers[1:-1]))]+[0]
    ann.create_sparse_array(opts.connection_rate, layers)
if opts.verbose:
    # define paramters for total layers : [num_neurons[0],num_neurons[1],.....,num_neurons[-1]]. num_neurons[0] is the number of input variables and num_neurons[-1] is the number of output variables.
    print "The layer structure is following:"
    print layers
    
if (opts.set_new_parameters == "True") or (not opts.import_network):
    if opts.verbose:
        print "Setting newtork parameters"
    # activation functions : SIGMOID,SIGMOID_STEPWISE,SIGMOID_SYMMETRIC,SIGMOID_SYMMETRIC_STEPWISE,LINEAR,THERESHOLD,THRESHOLD_SYMMETRIC,GAUSSIAN,GAUSSIAN_SYMMETRIC,ELLIOT,ELLIOT_SYMMETRIC,LINEAR_PIECE,LINEAR_PIECE_SYMMETRIC,SIN_SYMMETRIC,COS_SYMMETRIC,SIN,COS
    if opts.hidden_activation == "SIGMOID_SYMMETRIC_STEPWISE":
	    ann.set_activation_function_hidden(libfann.SIGMOID_SYMMETRIC_STEPWISE)
    elif opts.hidden_activation == "GAUSSIAN":
	    ann.set_activation_function_hidden(libfann.GAUSSIAN)
    elif opts.hidden_activation == "GAUSSSIAN_SYMMETRIC":
	    ann.set_activation_function_hidden(libfann.GAUSSIAN_SYMMETRIC)
    elif opts.hidden_activation == "SIGMOID":
	    ann.set_activation_function_hidden(libfann.SIGMOID)
    elif opts.hidden_activation == "SIGMOID_STEPWISE":
	    ann.set_activation_function_hidden(libfann.SIGMOID_STEPWISE)
    ann.set_activation_steepness_hidden(opts.steep_hidden)

    if opts.output_activation == "SIGMOID_SYMMETRIC_STEPWISE":
	    ann.set_activation_function_output(libfann.SIGMOID_SYMMETRIC_STEPWISE)
    elif opts.output_activation == "GAUSSIAN":
	    ann.set_activation_function_output(libfann.GAUSSIAN)
    elif opts.output_activation == "GAUSSIAN_SYMMETRIC":
	    ann.set_activation_function_output(libfann.GAUSSIAN_SYMMETRIC)
    elif opts.output_activation == "SIGMOID":
	    ann.set_activation_function_output(libfann.SIGMOID)
    elif opts.output_activation == "SIGMOID_STEPWISE":
	    ann.set_activation_function_output(libfann.SIGMOID_STEPWISE)
    ann.set_activation_steepness_output(opts.steep_out)


########## Training Machine #####################
machine = opts.training_machine.split(",")

########## GA Training #####################
if "ga" in machine:
    if opts.verbose:
        print "Setting GA training parameters"

    # setting Initializator
    if opts.initializator == "gauss":
        gaInit=GAnnInitializators.G1DConnInitializatorGaussian
    elif opts.initializator == "uniform":
        gaInit=GAnnInitializators.G1DConnInitializatorUniform
    else:
        raise exceptions.ValueError("Wrong argument {0} for {1} ".format(opts.initializator, '--initializator'))

    # generating genome
    if opts.verbose:
        print "number of total connections : %i" % ann.get_total_connections()
    totalConnections = ann.get_total_connections()
    genome = G1DConnections.G1DConnections(totalConnections)
#if opts.import_network:
    genome.genomeList = GAnnUtils.ToGAnnConn(ann.get_connection_array())

    genome.evaluator.set(GAnnEvaluators.evaluateMSE)

    genome.setParams(rangemin=opts.range_min, rangemax=opts.range_max, layers=layers, bias=bias, gauss_mu=opts.gauss_mu, gauss_sigma=opts.gauss_sigma,connectionrate=opts.connection_rate)
    genome.initializator.set(gaInit)
    genome.mutator.set(GAnnMutators.G1DConnUnbiasedMutateWeights)
    genome.crossover.set(GAnnCrossovers.G1DConnCrossoverWeights)
    ga = GAnnGA.GAnnGA(genome, ann, train_data)
    ga.setMutationRate(opts.mutation_rate)
    ga.setPopulationSize(opts.population)
    ga.setGenerations(opts.generations)
    if opts.verbose:
        print "Start running GA"
        print "initial GA MSE : %f" % ann.test_data(train_data)
    start_ga_time = time.time()
    ga.evolve(freq_stats=1)
    end_ga_time = time.time()
    if opts.verbose:
        print "Time elpased for Pyevolve training: %f seconds" % (end_ga_time - start_ga_time)

    # choose best initial connection weights
    best=ga.bestIndividual()
    if opts.verbose:
        print "number of total connections in best genome : %i" % len(best.toList())
    # set initial connection weights on ann
    ann.set_weight_array(best.toList())
    ann.save(opts.saving_results)
    if opts.verbose:
        print "Best initial connection wetighs are set on Network"
        print "current GA MSE : %f" % ann.test_data(train_data)
        print "GA Training with Pyevolve is done."
    if "fann" not in machine:
        if opts.verbose:
            print "Trained Network by Genetic Algorithm is saved in \n%s." % opts.saving_results

########## FANN Training #####################
if "fann" in machine:
    if opts.verbose:
        print "Setting FANN training parameters"
    learning_rate = opts.learning_rate
    desired_error = opts.desired_error
    max_iterations = opts.max_epochs
    iterations_between_reports = 100

    if (opts.set_new_parameters == "True") or (not opts.import_network):
        #ann.set_input_scaling_params(train_data,-100,100)

        # training algorithm : INCREMENTAL, BATCH, RPROP, QUICKPROP
        # learning rate and learning momentum are not used during RPROP training
        ann.set_training_algorithm(libfann.TRAIN_RPROP)
        ann.set_learning_rate(learning_rate)
        #ann.set_training_algorithm(libfann.TRAIN_INCREMENTAL)
        #ann.set_learning_momentum(learning_momentum)

        #rpop parameters
        ann.set_rprop_increase_factor(opts.increase_factor)# >1, default 1.2
        ann.set_rprop_decrease_factor(opts.decrease_factor)# <1, default 0.5
        ann.set_rprop_delta_min(opts.delta_min)# small positive number, default 0.0
        ann.set_rprop_delta_max(opts.delta_max)# positive number, default 50.0

    # initial connection weights for FANN  
    if "ga" not in machine:
        #ann.init_weights(train_data)
        ann.randomize_weights(opts.weights_min,opts.weights_max)

    if opts.verbose:
        ann.print_parameters()
        print "Start training network with %s-algorithm in FANN" % ann.get_training_algorithm()

    start_fann_time = time.time()
    ann.train_on_data(train_data, max_iterations, iterations_between_reports, desired_error)
    end_fann_time = time.time()

    if opts.verbose:
        print "Time elpased for FANN training: %f seconds" % (end_fann_time - start_fann_time)
        if "ga" not in machine:
            print "Trained Network by only BackPropagation Algorithm in FANN library is saved in \n%s." % opts.saving_results
        else:
            print "Trained Network by GA+FANN is saved in \n%s." % opts.saving_results
    ann.save(opts.saving_results)
    end_time = time.time()
    if opts.verbose:
        print "Total Running Time: %f seconds" % (end_time - start_time)
