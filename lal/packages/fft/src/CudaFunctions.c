#include <lal/LALDatatypes.h>
#include <lal/XLALError.h>
#include <cuda_runtime.h>
#include "CudaFunctions.h"


REAL4 *XLALCudaMallocReal(UINT4 size)
{
    static const char *func = "XLALCudaMallocReal";
    REAL4 *d_data;
    cudaMalloc( (void **)&d_data, sizeof(REAL4) * size );

    if( !d_data )
        XLAL_ERROR_NULL( func, XLAL_ENOMEM );
    return d_data;
}

COMPLEX8 *XLALCudaMallocComplex(UINT4 size)
{
    static const char *func = "XLALCudaMallocComplex";
    COMPLEX8 *d_data;
    cudaMalloc( (void **)&d_data, sizeof(COMPLEX8) * size );

    if( !d_data )
        XLAL_ERROR_NULL( func, XLAL_ENOMEM );
    return d_data;
}

void XLALCudaFree(void *d_data)
{
    static const char *func = "XLALCudaFree";
    if( !d_data )
        XLAL_ERROR_VOID( func, XLAL_EFAULT );
    cudaFree(d_data);
}
