# -*- coding: utf-8 -*-
#
#       cbcBayesMCMC2pos.py
#
#       Copyright 2016
#       Carl-Johan Haster <carl-johan.haster@ligo.org>
#
#       Following methodology from cbcBayesThermoInt.py and lalapps_nest2pos.py
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

from six.moves import range

import matplotlib
matplotlib.use('Agg') #sets backend to not need to open windows

import numpy as np
import astropy.table as apt
import scipy.integrate as si
from optparse import OptionParser

import h5py
from lalinference import git_version
from lalinference import bayespputils as bppu
from lalinference.io import read_samples, write_samples

from lalinference import LALINFERENCE_PARAM_LINEAR as LINEAR
from lalinference import LALINFERENCE_PARAM_CIRCULAR as CIRCULAR
from lalinference import LALINFERENCE_PARAM_FIXED as FIXED
from lalinference import LALINFERENCE_PARAM_OUTPUT as OUTPUT

__author__="Carl-Johan Haster <carl-johan.haster@ligo.org>>"
__version__= "git id %s"%git_version.id
__date__= git_version.date

def multipleFileCB(opt, opt_str, value, parser):
    args=[]

    def floatable(str):
      try:
        float(str)
        return True
      except ValueError:
        return False

    for arg in parser.rargs:
      # stop on --foo like options
      if arg[:2] == "--" and len(arg) > 2:
        break
      # stop on -a, but not on -3 or -3.0
      if arg[:1] == "-" and len(arg) > 1 and not floatable(arg):
        break
      args.append(arg)

    del parser.rargs[:len(args)]
    #Append new files to list if some already specified
    if getattr(parser.values, opt.dest):
        oldargs = getattr(parser.values, opt.dest)
        oldargs.extend(args)
        args = oldargs
    setattr(parser.values, opt.dest, args)

mcmc_group_id = '/lalinference/lalinference_mcmc'

def reassign_metadata(new_posterior, original_hdf5):

	# Make sure output file has same metadata as original
	# input hdf5 file

	base_file = read_samples(original_hdf5)

	mcmc_diagnostics_params = ['nLocalTemps','randomSeed']

	meta_dict = {}

	for colname, column in base_file.columns.items():
		meta_dict[colname] = column.meta

	for colname, column in new_posterior.columns.items():
		if colname in mcmc_diagnostics_params:
			column.meta = {'vary': OUTPUT}
			# these parameters are fixed within a run,
			# but doesn't have to be equal between runs.
		elif colname in meta_dict:
			column.meta = meta_dict[colname]
		elif 'cos'+colname in meta_dict:
			column.meta = meta_dict['cos'+colname]
		elif 'sin'+colname in meta_dict:
			column.meta = meta_dict['sin'+colname]
		elif 'log'+colname in meta_dict:
			column.meta = meta_dict['log'+colname]
		elif colname.startswith('chain_'):
			column.meta = {'vary': OUTPUT}
			# same argument as with mcmc_diagnostics_params
		else:
			column.meta = {'vary': FIXED}

	return new_posterior

def downsample_and_evidence(data_hdf5, deltaLogP=None, fixedBurnin=None, nDownsample=None, verbose=False):

    # Remove burnin from beginning of MCMC-chains, downsample by the chain's respective autocorrelation length.
    # Compute the evidence for the set of parallel tempered chains through a thermodynamic integral

	if not data_hdf5.lower().endswith(('.hdf5', '.hdf', '.h5')):
		print('cbcBayesMCMC2pos only suports hdf5 input, for older file formats plese revert to cbcBayesThermoInt and cbcBayesPosProc')
		sys.exit(1)

	peparser = bppu.PEOutputParser('hdf5')

	ps, samps = peparser.parse(data_hdf5, deltaLogP=deltaLogP, fixedBurnins=fixedBurnin,
		nDownsample=nDownsample, tablename=None)
	posterior_samples = apt.Table(samps, names=ps)

	highTchains = []
	for i in range(1,int(posterior_samples['nTemps'][0])):
		ps, samps = peparser.parse(data_hdf5, deltaLogP=deltaLogP, fixedBurnins=fixedBurnin,
			nDownsample=nDownsample, tablename='chain_'+str('%02.f' %i))
		highTchains.append(apt.Table(samps, names=ps))
		if verbose: print('chain_'+str('%02.f' %i)+' at a temperature '+str(highTchains[i-1]['temperature'].mean()))

	betas = np.zeros(len(highTchains)+1)
	logls = np.zeros_like(betas)

	betas[0] = 1./np.median(posterior_samples['temperature'])
	logls[0] = np.median(posterior_samples['logl'])

	for i in range(len(highTchains)):
		betas[i+1] = 1./np.median(highTchains[i]['temperature'])
		logls[i+1] = np.median(highTchains[i]['logl'])

	inds = np.argsort(betas)[::-1]

	betas = betas[inds]
	logls = logls[inds]

	# Now extend to infinite temperature by copying the last <log(L)>.
	# This works as long as the chains have extended to high enough
	# temperature to sample the prior.
	# If infinite temperature is already included, this 'duplicate'
	# will not change the final evidence.
	ebetas = np.concatenate((betas, [0.0]))
	elogls = np.concatenate((logls, [logls[-1]]))

	ebetas2 = np.concatenate((betas[::2], [0.0]))
	elogls2 = np.concatenate((logls[::2], [logls[::2][-1]]))

	evidence = -si.trapz(elogls, ebetas)
	evidence2 = -si.trapz(elogls2, ebetas2)

	posterior_samples['chain_log_evidence'] = evidence
	posterior_samples['chain_delta_log_evidence'] = np.absolute(evidence - evidence2)
	posterior_samples['chain_log_noise_evidence'] = posterior_samples['nullLogL']
	posterior_samples['chain_log_bayes_factor'] = posterior_samples['chain_log_evidence'] - posterior_samples['chain_log_noise_evidence']

	if verbose:
		print('logZ = '+str(posterior_samples['chain_log_evidence'][0])+'+-'+str(posterior_samples['chain_delta_log_evidence'][0]))
		print('logB_SN = '+str(posterior_samples['chain_log_bayes_factor'][0]))

	posterior_samples = reassign_metadata(posterior_samples, data_hdf5)

	return posterior_samples


def weight_and_combine(pos_chains, verbose=False):

	# Combine several posterior chains into one
	# weighted by their relative evidence

	log_evs = np.zeros(len(pos_chains))
	log_noise_evs = np.zeros_like(log_evs)

	for i in range(len(pos_chains)):
		log_evs[i] = pos_chains[i]['chain_log_evidence'][0]
		log_noise_evs[i] = pos_chains[i]['chain_log_noise_evidence'][0]
	if verbose: print('Computed log_evidences: %s'%(str(log_evs)))

	max_log_ev = log_evs.max()

	fracs=[np.exp(log_ev-max_log_ev) for log_ev in log_evs]
	if verbose: print('Relative weights of input files: %s'%(str(fracs)))

	Ns=[fracs[i]/len(pos_chains[i]) for i in range(len(fracs))]
	Ntot=max(Ns)
	fracs=[n/Ntot for n in Ns]
	if verbose: print('Relative weights of input files taking into account their length: %s'%(str(fracs)))

	final_posterior = pos_chains[0][np.random.uniform(size=len(pos_chains[0]))<fracs[0]]

	for i in range(1,len(pos_chains)):
		final_posterior = apt.vstack([final_posterior,
			pos_chains[i][np.random.uniform(size=len(pos_chains[i]))<fracs[i]]])

	final_log_evidence = reduce(np.logaddexp, log_evs) - np.log(len(log_evs))
	final_log_noise_evidence = reduce(np.logaddexp, log_noise_evs) - np.log(len(log_noise_evs))

	run_level= 'lalinference/lalinference_mcmc'
	metadata = {}
	metadata[run_level] = {}
	metadata[run_level]['log_bayes_factor'] = final_log_evidence - final_log_noise_evidence
	metadata[run_level]['log_evidence'] = final_log_evidence
	metadata[run_level]['log_noise_evidence'] = final_log_noise_evidence
	metadata[run_level]['log_max_likelihood'] = final_posterior['logl'].max()

	# This has already been burned-in and downsampled,
	# remove the cycle column to stop cbcBayesPosProc
	# from doing it again.
	final_posterior.remove_column('cycle')

	return final_posterior, metadata

USAGE='''%prog [options] PTMCMC_datafile.hdf5 [PTMCMC_datafile2.hdf5 ...]
Compute the evidence for a set of parallel tempered MCMC chains
thourgh thermodynamical integration. If using several input PTMCMC files,
combine them into one set of posterior samples, weighted by their relative
evidences.
'''

if __name__ == '__main__':
	parser = OptionParser(USAGE)
	parser.add_option(
		'-p', '--pos', action='store', type='string', default=None,
		help='Output file for posterior samples', metavar='posterior.hdf5')
	parser.add_option(
		'-v', '--verbose', action='store_true', default=False,
		help='Print some additional information')
	parser.add_option('-d','--data',dest='data',action='callback',
		callback=multipleFileCB,help='datafile')
	parser.add_option(
		'-s','--downsample',action='store',default=None,
		help='Approximate number of samples to record in the posterior',type='int')
	parser.add_option(
		'-l','--deltaLogP',action='store',default=None,
		help='Difference in logPosterior to use for convergence test.',type='float')
	parser.add_option(
		'-b','--fixedBurnin',dest='fixedBurnin',action="callback",
		callback=multipleFileCB,help='Fixed number of iteration for burnin.')
	opts, args = parser.parse_args()

	datafiles=[]
	if args:
		datafiles = datafiles + args
	if opts.data:
		datafiles = datafiles + opts.data

	if opts.fixedBurnin:
	# If only one value for multiple chains, assume it's to be applied to all chains
		if len(opts.fixedBurnin) == 1:
			fixedBurnins = [int(opts.fixedBurnin[0]) for df in datafiles]
		else:
			fixedBurnins = [int(fixedBurnin) for fixedBurnin in opts.fixedBurnin]
	else:
		fixedBurnins = [None]*len(datafiles)

	chain_posteriors = []

	for i in range(len(datafiles)):
		chain_posteriors.append(downsample_and_evidence(datafiles[i],
			deltaLogP=opts.deltaLogP, fixedBurnin=fixedBurnins[i], nDownsample=opts.downsample, verbose=opts.verbose))

	final_posterior, metadata = weight_and_combine(chain_posteriors, verbose=opts.verbose)

	write_samples(final_posterior, opts.pos,
		path='/'.join(['','lalinference','lalinference_mcmc','posterior_samples']), metadata=metadata)
