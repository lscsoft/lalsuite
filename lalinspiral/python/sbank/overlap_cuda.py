# Copyright (C) 2011  Nickolas Fotopoulos
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

from numpy import float32, complex64
import pycuda.autoinit  # initializes upon import
import pycuda.driver as cuda
from pycuda.reduction import ReductionKernel
from pycuda.elementwise import ElementwiseKernel


class create_workspace_cache(dict):
    def __init__(self, *args, **kwargs):
        self.stream = cuda.Stream()

        self.conj_knl = ElementwiseKernel(
            "pycuda::complex<float> *h",
            "h[i]._M_im = -h[i].imag()",
            preamble="#include <pycuda-complex.hpp>")

        self.imult_knl = ElementwiseKernel(
            "pycuda::complex<float> *a, pycuda::complex<float> *b",
            "a[i] = a[i]*b[i]",
            preamble="#include <pycuda-complex.hpp>")

        self.max_abs_knl = ReductionKernel(dtype_out=float32(0.).dtype, neutral="0",
            reduce_expr="max(a,b)",
            map_expr="x[i].real()*x[i].real()+x[i].imag()*x[i].imag()",
            arguments="pycuda::complex<float> *x",
            preamble="#include <pycuda-complex.hpp>")
        dict.__init__(self, *args, **kwargs)

    def get_workspace(self, n):
        from pyfft.cuda import Plan as pycufftplan
        import pycuda.gpuarray as gpuarray

        ws = self.get(n)
        if ws: return ws
        return self.setdefault(n,
            (pycufftplan(int(n), stream=self.stream, normalize=False),
             gpuarray.empty(n, dtype=complex64(0.).dtype),
             gpuarray.empty(n, dtype=complex64(0.).dtype)))

def destroy_workspace_cache(ws):
    pass

def compute_match(tmplt, inj, workspace):
    min_len = min(inj.data.length, tmplt.data.length)
    max_len = max(inj.data.length, tmplt.data.length)
    n = 2 * (min_len - 1)
    plan, a_gpu, b_gpu = workspace.get_workspace(n)

    # Note that findchirp paper eq 4.2 defines a positive-frequency integral,
    # so we should only fill the positive frequencies (first half of zf).
    a_gpu.fill(0.)  # reusing memory, so must reset it
    if tmplt.data.data.dtype == complex64:
        a_gpu[:min_len].set(tmplt.data.data[:min_len])  # a = tmplt
    else:
        a_gpu[:min_len].set(tmplt.data.data[:min_len].astype(complex64))  # a = tmplt
    b_gpu.fill(0.)  # reusing memory, so must reset it
    if inj.data.data.dtype == complex64:
        b_gpu[:min_len].set(inj.data.data[:min_len])  # b = inj
    else:
        b_gpu[:min_len].set(inj.data.data[:min_len].astype(complex64))  # b = inj

    workspace.conj_knl(b_gpu)  # b = conj(b)
    workspace.imult_knl(a_gpu, b_gpu)  # a *= b

    plan.execute(a_gpu, data_out=b_gpu, inverse=True)

    return 4. * inj.deltaF * (workspace.max_abs_knl(b_gpu).get() * (min_len - 1.) / (max_len - 1.))**0.5 / n
