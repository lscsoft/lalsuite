import numpy as np
from math import ceil
from optparse import OptionParser
import os
from scipy.linalg import solve_triangular
import lal
import lalsimulation

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
            return

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
        max_cols = max(1, max_elements // m)
        max_rows =  max_elements // max_cols

    else:
        # prioritize processing as many rows of A as possible
        max_rows = max(1, max_elements // n)
        max_cols =  max_elements // max_rows

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

def ehat(j, length):
    """Return a unit vector whose j-th component is 1"""
    ehat = np.zeros(length)
    ehat[j] = 1.0
    return ehat

def DEIM(basis):
    """Calculate interpolation nodes following the Algorithm 5 in J Sci Comput
    (2013) 57:604-637.

    Parameters
    ----------
    basis : ndarray
        ndarray whose i-th row is the i-th basis vector

    Return
    ------
    nodes : array
        interpolation nodes
    B : ndarray
        ndarray whose i-th row is the weight for the i-th frequency node.
   """
    vec_num, vec_len = basis.shape
    # Step (1)
    j = abs(basis[0]).argmax()
    # Step (2)
    UT = np.zeros(shape=(vec_num, vec_len), dtype=complex)
    PT = np.zeros(shape=(vec_num, vec_len))
    UT[0] = basis[0]
    PT[0] = ehat(j, vec_len)
    PTU = np.zeros(shape=(vec_num, vec_num), dtype=complex)
    nodes = [j]
    # Step (3)
    for i in range(vec_num - 1):
        i += 1
        ei = basis[i]
        # Step (4)
        PTU[i-1][:i:] = np.array([UT[k][nodes[i-1]] for k in range(i)])
        c = solve_triangular(
            PTU[:i:, :i:], np.array([ei[node] for node in nodes]),
            lower=True, check_finite=False)
        # Step (5)
        r = ei - c.dot(UT[:i:])
        # Step (6)
        j = abs(r).argmax()
        # Step (7)
        UT[i] = r
        PT[i] = ehat(j, vec_len)
        nodes.append(j)

    # Calculate B
    VT = basis
    B = np.linalg.inv(VT.dot(PT.transpose())).dot(VT)

    # ordering
    B = np.array([B[i] for i in np.argsort(nodes)])
    nodes.sort()

    return nodes, B

def construct_nodes(
    selected_params, flow, fhigh, deltaF, approximant, quadratic
):
    """Construct frequency nodes and weights from parameter values selected by greedy algorithm.

    Parameters
    ----------
    selected_params : ndarray
	ndarray whose i-th row is the i-th parameter set selected by greedy
        algorithm
    flow : float
    fhigh : float
    deltaF : float
    approximant : string
    quadratic : bool
        construct nodes for products of waveforms when this is True

    Return
    ------
    f_nodes : ndarray
        interpolation frequency nodes
    B : ndarray
        ndarray whose i-th roq is the weight for the i-th frequency node.
    """
    # generate waveforms at selected parameters
    start_index = int(flow / deltaF)
    def _waveform(param):
        # FIXME: Need to be fixed when we take into account precession.
        m1, m2, a1z, a2z = param
        hp, _  = lalsimulation.SimInspiralChooseFDWaveform(
            m1 * lal.MSUN_SI, m2 * lal.MSUN_SI,
            0.0, 0.0, a1z, 0.0, 0.0, a2z,
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            deltaF, flow, fhigh, flow, None,
            lalsimulation.GetApproximantFromString(approximant)
        )
        hp = hp.data.data[start_index::]
        if quadratic:
            return np.conj(hp) * hp
        else:
            return hp
    basis = np.array([_waveform(param) for param in selected_params]).transpose()

    # orthogonalize waveforms
    basis, _ = np.linalg.qr(basis)
    basis = basis.transpose()

    # discrete empirical interpolation
    nodes, B = DEIM(basis)
    return np.array([flow + deltaF * node for node in nodes]), B

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
parser.add_option("-a", "--approx", type='string',
                      action="store",
                      dest="approximant",
                      help="approximant name",)
parser.add_option("-s", "--seglen", type=float,
                      action="store",
                      dest="seglen",
                      help="segment length",)
parser.add_option("-f", "--fLow", type=float,
                      action="store",
                      dest="fLow",
                      help="low frequency cut off",)
parser.add_option("-u", "--fHigh", type=float,
                      action="store",
                      dest="fHigh",
                      help="high frequency cut off",)
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

basis_params = np.loadtxt(os.path.join(options.b_matrix_directory, "params.dat")).T
flow_in_params, fhigh_in_params, deltaF_in_params = basis_params[0], basis_params[1], 1. / basis_params[2]

B_linear_path = os.path.join(options.b_matrix_directory, "B_linear.npy")
B_quadratic_path = os.path.join(options.b_matrix_directory, "B_quadratic.npy")
fnodes_linear_path = os.path.join(options.b_matrix_directory, "fnodes_linear.npy")
fnodes_quadratic_path = os.path.join(options.b_matrix_directory, "fnodes_quadratic.npy")
params_linear_path = os.path.join(options.b_matrix_directory, "selected_params_linear.npy")
params_quadratic_path = os.path.join(options.b_matrix_directory, "selected_params_quadratic.npy")
if os.path.exists(B_linear_path) and os.path.exists(B_quadratic_path) and \
    os.path.exists(fnodes_linear_path) and os.path.exists(fnodes_quadratic_path):
    B_linear = np.load(B_linear_path)
    B_quadratic = np.load(B_quadratic_path)
    fnodes_linear = np.load(fnodes_linear_path)
    fnodes_quadratic = np.load(fnodes_quadratic_path)
elif os.path.exists(params_linear_path) and os.path.exists(params_quadratic_path):
    selected_params_linear = np.load(params_linear_path)
    selected_params_quadratic = np.load(params_quadratic_path)
    fnodes_linear, B_linear = construct_nodes(selected_params_linear, flow_in_params, fhigh_in_params, deltaF_in_params, options.approximant, False)
    fnodes_quadratic, B_quadratic = construct_nodes(selected_params_linear, flow_in_params, fhigh_in_params, deltaF_in_params, options.approximant, True)
else:
    print("No ROQ data found. Please make sure that you have B_(linear, quadratic).npy "
          "and fnodes_(linear, quadratic).npy, or selected_params_(linear, quadratic).npy")
    raise

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
    if options.seglen:
        deltaF = 1./options.seglen
    else:
        deltaF = fseries[1] - fseries[0]
    if options.fHigh:
        fHigh = options.fHigh
    else:
        fHigh = fseries[-1]
    fHigh_index = int(fHigh / deltaF)

    print('Desired fHigh is', fHigh,', actual fHigh is', fseries[fHigh_index])

    if options.fLow:
        fLow = options.fLow
        scale_factor = flow_in_params / fLow

    else:
        fLow = flow_in_params
        assert fHigh == fhigh_in_params
    fLow_index = int(fLow / deltaF)

    print('Desired fLow is', fLow,', actual fLow is', fseries[fLow_index])

    fseries = fseries[fLow_index:fHigh_index]
    data = data[fLow_index:fHigh_index]

    psdfile = np.column_stack( np.loadtxt(options.psd_file[i]) )
    psd = psdfile[1]

    psd[-1] = psd[-1 -1 ]
    psd = psd[fLow_index:fHigh_index]
    data /= psd

    # only get frequency components up to fHigh
    B_linear = B_linear.T[0:(fHigh_index - fLow_index)][:].T
    B_quadratic = B_quadratic.T[0:(fHigh_index-fLow_index)][:].T
    print(B_linear.shape[1], B_quadratic.shape[1], len(data), len(psd))
    assert len(data) == len(psd) == B_linear.shape[1] == B_quadratic.shape[1]



    #for the dot product, it's convenient to work with transpose of B:
    B_linear = B_linear.T
    B_quadratic = B_quadratic.T

    for k in range(len(data)):
        if np.isnan(data[k].real):
            data[k] = 0+0j

    #0.045 comes from the diameter of the earth in light seconds: the maximum time-delay between earth-based observatories
    tcs = np.linspace(relative_tc_shift - options.dt - 0.045, relative_tc_shift + options.dt + 0.045, ceil(2.*(options.dt+0.045) / options.delta_tc) )# array of relative time shifts to be applied to the data


    tc_shifted_data = np.zeros([len(tcs), len(fseries)], dtype=complex)  # array to be filled with data, shifted by discrete time tc

    print("time steps = "+str(len(tcs)))
    for j in range(len(tcs)):

        tc = tcs[j]

        #exp_2pi_i_tc = np.array(np.exp(1j*2.*np.pi*fseries*tc))
        #data_tc = np.array(data*exp_2pi_i_tc)
        tc_shifted_data[j] = data * np.exp(1j*2.*np.pi*fseries*tc)

    #tc_shifted_data = tc_shifted_data.T#np.array(tc_shifted_data).T



    #*************************************************************************** #
    print("Computing weights for "+ifo)
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
    print("Weights have been computed for "+ifo)

    i += 1

    #save the fnodes as a dat file if they're not already:

if scale_factor:
    print("scale factor = %f"%scale_factor)
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
