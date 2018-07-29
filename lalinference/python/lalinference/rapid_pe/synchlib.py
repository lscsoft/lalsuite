import socket as pysocket
import errno
import numpy

def init_socket(address, port):
    """
    Lower level socket initialization routine.
    """
    sock = pysocket.socket(pysocket.AF_INET, pysocket.SOCK_STREAM)
    sock.connect((address, port))
    return sock

def prepare_data(data, dims, offset):
    """
    Pack the samples starting from 'offset' in 'data' enumerated from 'dims' into an array.
    """
    data = numpy.vstack([data[d][offset:] for d in dims]).T
    return data.astype(numpy.float64)

def send_samples(data, address=None, port=1890, socket=None, verbose=False):
    """
    Send a nsamp x dim array of samples ('data') to the 'address' and 'port'. If 'socket' is provided, it will be used instead of generating a new one.
    """

    #
    # Generate a socket if necessary
    #
    if socket is None:
        sender = pysocket.socket(pysocket.AF_INET, pysocket.SOCK_STREAM)
    else:
        sender = socket

    msg_len = numpy.prod(data.shape) * 8
    if verbose:
        print("Sending %d bytes to %s:%d" % (msg_len, address, port))

    try:
        if socket is None:
            sender.connect((address, port))
        sender.send(data.tobytes("C"))
    except pysocket.error as sockerr:
        if verbose:
            print("Unable to send data.")
        # FIXME: Reenable this
        #if sockerr.errno != errno.ECONNREFUSED or sockerr.errno:
            #raise sockerr

    return sender
