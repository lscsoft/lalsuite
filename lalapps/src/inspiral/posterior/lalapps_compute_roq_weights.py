import numpy as np
import sys
from math import ceil
from optparse import OptionParser
import os

data_dir = './'
# load in data from file

parser = OptionParser(usage="usage: %prog [options]",
                          version="%prog")
parser.add_option("-d", "--data", type='string',
                      action="append", 
                      dest="data_file",
                      help="data file",)
parser.add_option("-p", "--psd", type='string',
                      action="append", # optional because action defaults to "store"
                      dest="psd_file",
                      help="psd file",)
parser.add_option("-t", "--time_prior", type=float,
                      action="store", # optional because action defaults to "store"
                      dest="dt",
                      help="width of time prior",)
parser.add_option("-i", "--ifo", type='string',
                      action="append", 
                      dest="IFOs",
                      help="list of ifos",)
parser.add_option("-s", "--seglen", type=float,
                      action="store", 
                      dest="seglen",
                      help="segment length",)
parser.add_option("-f", "--fLow", type=float,
                      action="store",
                      dest="fLow",
                      help="low frequency cut off",)
parser.add_option("-T", "--delta_tc", type=float,
                      action="store",
                      dest="delta_tc",
                      help="width of tc subdomain",)
parser.add_option("-B", "--basis-set", type='string',
                      action="store",
                      dest="b_matrix_path",
                      help="B matrix",)
parser.add_option("-o", "--out", type='string',
                  action="store",
                  dest="outpath",
                  default="./",
                  help="output path",)

(options, args) = parser.parse_args()

B = np.load(options.b_matrix_path)

def BuildWeights(data, B, deltaF):

        ''' for a data array and reduced basis compute roq weights
        
        B: (reduced basis element)*invV (the inverse Vandermonde matrix)
        data: data set 
        PSD: detector noise power spectral density (must be same shape as data)
        deltaF: integration element df

        '''

        weights = np.dot(B, data.conjugate()) * deltaF * 4.

        return weights.T
##################################

relative_tc_shift = options.seglen - 2. 

# loop over ifos 

i=0

for ifo in options.IFOs:

	dat_file = np.column_stack( np.loadtxt(options.data_file[i]) )
	data = dat_file[1] + 1j*dat_file[2]
	fseries = dat_file[0]
        deltaF = fseries[1] - fseries[0]
	fseries = fseries[int(options.fLow/deltaF):len(fseries)]
	data = data[int(options.fLow/deltaF):len(data)]


	psdfile = np.column_stack( np.loadtxt(options.psd_file[i]) )
	psd = psdfile[1]
	psd = psd[int(options.fLow/deltaF):len(psd)]
	data /= psd


	# print len(data),len(psd),len(basis_set)

	assert len(data) == len(psd) == B.shape[1]

	for k in range(len(data)):
		if np.isnan(data[k].real):
			data[k] = 0+0j

	tc_shifted_data = []  # array to be filled with data, shifted by discrete time tc

	tcs = np.linspace(relative_tc_shift - options.dt - 0.022, relative_tc_shift + options.dt + 0.022, ceil(2.*(options.dt+0.022) / options.delta_tc) )# array of relative time shifts to be applied to the data
	print "time steps = "+str(len(tcs))
	for j in range(len(tcs)):

	        tc = tcs[j]

	        exp_2pi_i_tc = np.array(np.exp(1j*2.*np.pi*fseries*tc))
        	data_tc = np.array(data*exp_2pi_i_tc)
	        tc_shifted_data.append(data_tc)

	tc_shifted_data = np.array(tc_shifted_data).T

	#*************************************************************************** #

	weights_path = os.path.join(options.outpath,"weights_%s.dat"%ifo)

	weights_file = open(weights_path, "wb")

	print "Computing weights for "+ifo
	weights = BuildWeights(tc_shifted_data, B, deltaF)
	print "Weights have been computed for "+ifo

	(weights.T).tofile(weights_file)

	weights_file.close()
	i += 1

size_file_path = os.path.join(options.outpath,"roq_sizes.dat")
np.savetxt(size_file_path,np.array((len(tcs),B.shape[0],B.shape[1])),fmt='%u')
