/*
 * Copyright (C) 2017  Leo Singer
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 */


#include "bayestar_moc.h"
#include <math.h>
#include <gsl/gsl_math.h>
#include <stdlib.h>
#include <string.h>
#include <chealpix.h>
#include <lal/XLALError.h>


uint64_t nest2uniq64(uint8_t order, uint64_t nest)
{
    return nest + ((int64_t) 1 << 2 * (order + 1));
}


int8_t uniq2order64(uint64_t n)
{
    if (n < 4)
        return -1;

    int8_t order;
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || defined(_M_IX86)
    int64_t o;
    asm("bsrq %1, %0\n\t"
        : "=r" (o)
        : "rm" (n));
    order = o;
#else
    for (order = -1; n; n >>= 1, order ++)
        /* noop */;
#endif
    return (order >> 1) - 1;
}


double uniq2pixarea64(uint64_t n)
{
    return ldexp(M_PI / 3, -2 * uniq2order64(n));
}


int8_t uniq2nest64(uint64_t *n)
{
    int8_t order = uniq2order64(*n);
    *n -= (uint64_t) 1 << 2 * (order + 1);
    return order;
}


void pix2ang_uniq64(uint64_t ipix, double *theta, double *phi)
{
    int8_t order = uniq2nest64(&ipix);
    int64_t nside = (int64_t) 1 << order;
    pix2ang_nest64(nside, ipix, theta, phi);
}


void *moc_rasterize64(const void *pixels, size_t offset, size_t itemsize, size_t len, size_t *npix)
{
    const size_t pixelsize = offset + itemsize;
    const void *pixel = (const char *) pixels + (len - 1) * pixelsize;
    uint64_t ipix = *(const uint64_t *) pixel;
    const int8_t max_order = uniq2order64(ipix);
    *npix = 12 * ((size_t) 1 << 2 * max_order);
    void *ret = calloc(*npix, itemsize);
    if (!ret)
        XLAL_ERROR_NULL(XLAL_ENOMEM, "not enough memory to allocate image");

    for (size_t i = 0; i < len; i ++)
    {
        pixel = (const char *) pixels + i * pixelsize;
        ipix = *(const uint64_t *) pixel;
        const int8_t order = uniq2nest64(&ipix);
        const size_t reps = (size_t) 1 << 2 * (max_order - order);
        for (size_t j = 0; j < reps; j ++)
            memcpy((char *) ret + (ipix * reps + j) * itemsize, (const char *) pixel + offset, itemsize);
    }

    return ret;
}
