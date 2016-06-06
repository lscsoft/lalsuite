import numpy as np
import sys
from math import ceil
from optparse import OptionParser
import os

#####################################################################

def _block_slices(dim_size, block_size):
    """Generator that yields slice objects for indexing into 
    sequential blocks of an array along a particular axis
    """
    count = 0
    while True:
        yield slice(count, count + block_size, 1)
        count += block_size
        if count > dim_size:
            raise StopIteration

def blockwise_dot(A, B, deltaF, max_elements=int(2**27), out=None):
    """
    Computes the dot product of two matrices in a block-wise fashion. 
    Only blocks of `A` with a maximum size of `max_elements` will be 
    processed simultaneously.
    """

    m,  n = A.shape
    n1, o = B.shape

    if n1 != n:
        raise ValueError('matrices are not aligned')

    if A.flags.f_contiguous:
        # prioritize processing as many columns of A as possible
        max_cols = max(1, max_elements / m)
        max_rows =  max_elements / max_cols

    else:
        # prioritize processing as many rows of A as possible
        max_rows = max(1, max_elements / n)
        max_cols =  max_elements / max_rows

    if out is None:
        out = np.empty((m, o), dtype=np.result_type(A, B))
    elif out.shape != (m, o):
        raise ValueError('output array has incorrect dimensions')

    for mm in _block_slices(m, max_rows):
        out[mm, :] = 0
        for nn in _block_slices(n, max_cols):
            A_block = A[mm, nn].copy()  # copy to force a read
            out[mm, :] += np.dot(A_block, B[nn, :]) * 4 * deltaF
            del A_block

    return out
######################################################################

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
                      dest="b_matrix_directory",
                      help="B matrix directory",)
parser.add_option("-o", "--out", type='string',
                  action="store",
                  dest="outpath",
                  default="./",
                  help="output path",)

(options, args) = parser.parse_args()

B_linear = np.load(options.b_matrix_directory + "/B_linear.npy")
B_quadratic = np.load(options.b_matrix_directory + "/B_quadratic.npy")
basis_params = np.loadtxt(options.b_matrix_directory + "/params.dat").T

def BuildWeights(data, B, deltaF):

        ''' for a data array and reduced basis compute roq weights
        
        B: (reduced basis element)*invV (the inverse Vandermonde matrix)
        data: data set 
        PSD: detector noise power spectral density (must be same shape as data)
        deltaF: integration element df

        '''
        weights = np.dot(data, B.conjugate()) * deltaF * 4.
	#weights = np.einsum('ij,jk->ik', B, data.conjugate() * deltaF * 4)
        return weights
##################################

relative_tc_shift = options.seglen - 2. 

# loop over ifos 

i=0
scale_factor = 0

for ifo in options.IFOs:

	dat_file = np.column_stack( np.loadtxt(options.data_file[i]) )
	data = dat_file[1] + 1j*dat_file[2]
	fseries = dat_file[0]
        deltaF = fseries[1] - fseries[0]
	fHigh = fseries[-1]
	
	if options.fLow: 
		fLow = options.fLow
		scale_factor = int(basis_params[0] / fLow)

	else:
		fLow = basis_params[0]

		assert fHigh == basis_params[1]

	fseries = fseries[int(fLow/deltaF):len(fseries)+1]
	data = data[int(fLow/deltaF):len(data)+1]


	psdfile = np.column_stack( np.loadtxt(options.psd_file[i]) )
	psd = psdfile[1]

	psd[-1] = psd[-1 -1 ]

	psd = psd[int(fLow/deltaF):len(psd)+1]
	data /= psd

	assert len(data) == len(psd) == B_linear.shape[1] == B_quadratic.shape[1]

	#for the dot product, it's convenient to work with transpose of B:
	B_linear = B_linear.T 
	B_quadratic = B_quadratic.T

	for k in range(len(data)):
		if np.isnan(data[k].real):
			data[k] = 0+0j


	tcs = np.linspace(relative_tc_shift - options.dt - 0.026, relative_tc_shift + options.dt + 0.026, ceil(2.*(options.dt+0.026) / options.delta_tc) )# array of relative time shifts to be applied to the data

	
	tc_shifted_data = np.zeros([len(tcs), len(fseries)], dtype=complex)  # array to be filled with data, shifted by discrete time tc

	print "time steps = "+str(len(tcs))
	for j in range(len(tcs)):

	        tc = tcs[j]

	        #exp_2pi_i_tc = np.array(np.exp(1j*2.*np.pi*fseries*tc))
        	#data_tc = np.array(data*exp_2pi_i_tc)
	        tc_shifted_data[j] = data * np.exp(1j*2.*np.pi*fseries*tc)

	#tc_shifted_data = tc_shifted_data.T#np.array(tc_shifted_data).T



	#*************************************************************************** #
	print "Computing weights for "+ifo
	weights_path_linear = os.path.join(options.outpath,"weights_linear_%s.dat"%ifo)
	weights_file_linear = open(weights_path_linear, "wb")
	
	max_block_gigabytes = 4
	max_elements = int((max_block_gigabytes * 2 ** 30) / 8) # max number of double complex elements

	weights_linear = blockwise_dot(tc_shifted_data, B_linear.conjugate(), deltaF, max_elements).T

	(weights_linear).tofile(weights_file_linear)
        weights_file_linear.close()

	del tc_shifted_data

	#*************************************************************************** #	
	weights_path_quadratic = os.path.join(options.outpath,"weights_quadratic_%s.dat"%ifo)
        weights_file_quadratic = open(weights_path_quadratic, "wb")
	weights_quadratic = (BuildWeights(1./psd, B_quadratic, deltaF).T).real

	(weights_quadratic).tofile(weights_file_quadratic)
        weights_file_quadratic.close()
	size_file_path = os.path.join(options.outpath,"roq_sizes.dat")

	#*************************************************************************** #
	
	B_linear = B_linear.T
	B_quadratic = B_quadratic.T
	
        np.savetxt(size_file_path,np.array((len(tcs),B_linear.shape[0],B_quadratic.shape[0],B_linear.shape[1])),fmt='%u')
	print "Weights have been computed for "+ifo

	i += 1

	#save the fnodes as a dat file if they're not already:


fnodes_linear = np.load(options.b_matrix_directory + "/fnodes_linear.npy")
fnodes_quadratic = np.load(options.b_matrix_directory + "/fnodes_quadratic.npy")

if scale_factor:

	fnodes_linear /= scale_factor
	fnodes_quadratic  /= scale_factor


fnodes_linear_path = os.path.join(options.outpath,"fnodes_linear.dat")
fnodes_linear_file = open(fnodes_linear_path, "wb")

fnodes_linear.tofile(fnodes_linear_file)
fnodes_linear_file.close()

fnodes_quadratic_path = os.path.join(options.outpath,"fnodes_quadratic.dat")
fnodes_quadratic_file = open(fnodes_quadratic_path, "wb")

fnodes_quadratic.tofile(fnodes_quadratic_file)
fnodes_quadratic_file.close()

tcs_file_path = os.path.join(options.outpath,"tcs.dat")
tcs_file = open(tcs_file_path, "wb")

tcs.tofile(tcs_file)
tcs_file.close()
