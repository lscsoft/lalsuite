//
// Copyright (C) 2020 Karl Wette
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA  02110-1301  USA
//

///
/// \file
/// \ingroup lalpulsar_bin_Weave
///

#include "config.h"

#include <cuda.h>
#include <cuda_runtime_api.h>

#include <lal/LALStdlib.h>

#define CUDA_BLOCK_SIZE 512

#define XLAL_CHECK_CUDA_CALL(...) do { \
  cudaError_t retn; \
  XLAL_CHECK ( ( retn = (__VA_ARGS__) ) == cudaSuccess, XLAL_EERR, "%s failed with return code %i", #__VA_ARGS__, retn ); \
  } while(0)

///
/// CUDA kernel to find the maximum of \a nvec vectors
///
__global__ void VectorsMaxREAL4CUDA( REAL4 *max, const REAL4 **vec, const size_t nvec, const size_t nbin ) {
  int k = threadIdx.x + blockDim.x*blockIdx.x;
  if (k >= nbin) {
    return;
  }
  max[k] = vec[0][k];
  for ( size_t j = 1; j < nvec; ++j ) {
    if ( vec[j][k] > max[k] ) {
      max[k] = vec[j][k];
    }
  }
}

///
/// Find the maximum of \a nvec vectors in \a vec[], of length \a nbin, and return the result in \a max
///
extern "C" int XLALVectorsMaxREAL4CUDA( REAL4 *max, const REAL4 **vec, const size_t nvec, const size_t nbin ) {
  XLAL_CHECK( max != NULL, XLAL_EFAULT );
  XLAL_CHECK( vec != NULL, XLAL_EFAULT );
  XLAL_CHECK( vec[0] != NULL, XLAL_EFAULT );
  XLAL_CHECK( nvec > 0, XLAL_EINVAL );
  XLAL_CHECK( nbin > 0, XLAL_EINVAL );

  VectorsMaxREAL4CUDA<<<(nbin + CUDA_BLOCK_SIZE - 1)/CUDA_BLOCK_SIZE, CUDA_BLOCK_SIZE>>>( max, vec, nvec, nbin );
  XLAL_CHECK_CUDA_CALL ( cudaGetLastError() );
  XLAL_CHECK_CUDA_CALL ( cudaDeviceSynchronize() );

  return XLAL_SUCCESS;
}

///
/// CUDA kernel to add \a nvec vectors
///
__global__ void VectorsAddREAL4CUDA( REAL4 *sum, const REAL4 **vec, const size_t nvec, const size_t nbin ) {
  int k = threadIdx.x + blockDim.x*blockIdx.x;
  if (k >= nbin) {
    return;
  }
  sum[k] = vec[0][k];
  for ( size_t j = 1; j < nvec; ++j ) {
    sum[k] += vec[j][k];
  }
}

///
/// Add \a nvec vectors in \a vec[], of length \a nbin, and return the result in \a sum
///
extern "C" int XLALVectorsAddREAL4CUDA( REAL4 *sum, const REAL4 **vec, const size_t nvec, const size_t nbin ) {
  XLAL_CHECK( sum != NULL, XLAL_EFAULT );
  XLAL_CHECK( vec != NULL, XLAL_EFAULT );
  XLAL_CHECK( vec[0] != NULL, XLAL_EFAULT );
  XLAL_CHECK( nvec > 0, XLAL_EINVAL );
  XLAL_CHECK( nbin > 0, XLAL_EINVAL );

  VectorsAddREAL4CUDA<<<(nbin + CUDA_BLOCK_SIZE - 1)/CUDA_BLOCK_SIZE, CUDA_BLOCK_SIZE>>>( sum, vec, nvec, nbin );
  XLAL_CHECK_CUDA_CALL ( cudaGetLastError() );
  XLAL_CHECK_CUDA_CALL ( cudaDeviceSynchronize() );

  return XLAL_SUCCESS;
}

// Local Variables:
// mode: c
// c-file-style: "linux"
// c-basic-offset: 2
// End:
