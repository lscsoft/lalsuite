/*
 * Copyright (C) 2012  Nickolas Fotopoulos
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <lal/XLALError.h>
#include <lal/FrequencySeries.h>
#include <lal/LALInspiralSBankOverlap.h>
#include <sys/types.h>

#define MAX_NUM_WS 32  /* maximum number of workspaces */
#define CHECK_OOM(ptr, msg) if (!(ptr)) { XLALPrintError((msg)); XLAL_ERROR_NULL(XLAL_ENOMEM); }

/*
 * set up workspaces
 */

WS *XLALCreateSBankWorkspaceCache(void) {
    WS *workspace_cache = calloc(MAX_NUM_WS, sizeof(WS));
    CHECK_OOM(workspace_cache, "unable to allocate workspace\n");
    return workspace_cache;
}

void XLALDestroySBankWorkspaceCache(WS *workspace_cache) {
    size_t k = MAX_NUM_WS;
    for (;k--;) {
        if (workspace_cache[k].n) {
            fftwf_destroy_plan(workspace_cache[k].plan);
            fftwf_free(workspace_cache[k].zf);
            fftwf_free(workspace_cache[k].zt);
        }
    }
    free(workspace_cache);
}

static WS *get_workspace(WS *workspace_cache, const size_t n) {
    if (!n) {
        fprintf(stderr, "Zero size workspace requested\n");
        abort();
    }

    /* if n already in cache, return it */
    WS *ptr = workspace_cache;
    while (ptr->n) {
        if (ptr->n == n) return ptr;
        if (++ptr - workspace_cache > MAX_NUM_WS) return NULL;  /* out of space! */
    }

    /* if n not in cache, ptr now points at first blank entry */
    ptr->n = n;
    ptr->zf = fftwf_malloc(sizeof(fftwf_complex) * n);
    CHECK_OOM(ptr->zf, "unable to allocate workspace array zf\n");
    memset(ptr->zf, 0, n * sizeof(*(ptr->zf)));
    ptr->zt = fftwf_malloc(sizeof(fftwf_complex) * n);
    CHECK_OOM(ptr->zf, "unable to allocate workspace array zt\n");
    memset(ptr->zt, 0, n * sizeof(*(ptr->zt)));
    ptr->plan = fftwf_plan_dft_1d(n, ptr->zf, ptr->zt, FFTW_BACKWARD, FFTW_MEASURE);
    CHECK_OOM(ptr->plan, "unable to allocate plan");

    return ptr;
}

/* by default, complex arithmetic will call built-in function __muldc3, which does a lot of error checking for inf and nan; just do it manually */
static void multiply_conjugate(complex float * restrict out, complex float *a, complex float *b, const size_t size) {
    size_t k = 0;
    for (;k < size; ++k) {
        const float ar = crealf(a[k]);
        const float br = crealf(b[k]);
        const float ai = cimagf(a[k]);
        const float bi = cimagf(b[k]);
        __real__ out[k] = ar * br + ai * bi;
        __imag__ out[k] = ar * -bi + ai * br;
    }
}

static double abs2(const complex float x) {
    const REAL8 re = crealf(x);
    const REAL8 im = cimagf(x);
    return re * re + im * im;
}

/* interpolate the peak with a parabolic interpolation */
static double vector_peak_interp(const double ym1, const double y, const double yp1) {
    const double dy = 0.5 * (yp1 - ym1);
    const double d2y = 2. * y - ym1 - yp1;
    return y + 0.5 * dy * dy / d2y;
}

/*
 * Returns the match for whitened, normalized, positive-frequency
 * COMPLEX8FrequencySeries inputs.
 */
REAL8 XLALInspiralSBankComputeMatch(const COMPLEX8FrequencySeries *inj, const COMPLEX8FrequencySeries *tmplt, WS *workspace_cache) {
    size_t min_len = (inj->data->length <= tmplt->data->length) ? inj->data->length : tmplt->data->length;
    size_t max_len = (inj->data->length >= tmplt->data->length) ? inj->data->length : tmplt->data->length;

    /* get workspace for + and - frequencies */
    size_t n = 2 * (min_len - 1);   /* no need to integrate implicit zeros */
    WS *ws = get_workspace(workspace_cache, n);
    if (!ws) {
        XLALPrintError("out of space in the workspace_cache\n");
        XLAL_ERROR_REAL8(XLAL_ENOMEM);
    }

    /* compute complex SNR time-series in freq-domain, then time-domain */
    /* Note that findchirp paper eq 4.2 defines a positive-frequency integral,
       so we should only fill the positive frequencies (first half of zf). */
    multiply_conjugate(ws->zf, inj->data->data, tmplt->data->data, min_len);
    fftwf_execute(ws->plan); /* plan is reverse */

    /* maximize over |z(t)|^2 */
    complex float *zdata = ws->zt;
    size_t k = n;
    ssize_t argmax = -1;
    REAL8 max = 0.;
    for (;k--;) {
        REAL8 temp = abs2(zdata[k]);
        if (temp > max) {
            argmax = k;
            max = temp;
        }
    }
    if (max == 0.) return 0.;

    /* refine estimate of maximum */
    REAL8 result;
    if (argmax == 0 || argmax == (ssize_t) n - 1)
        result = max;
    else
        result = vector_peak_interp(abs2(zdata[argmax - 1]), abs2(zdata[argmax]), abs2(zdata[argmax + 1]));

    /* normalization depends on vector length; correct for differing normalizations */
    result *= (min_len - 1.) / (max_len - 1.);

    /* compute match */
    /* return 4. * inj->deltaF * sqrt(result) / n; */  /* inverse FFT = reverse / n */
    return 4. * inj->deltaF * sqrt(result);  /* Ajith, 2012-11-09: The division by n is inconsistent with the revised convention of the normalization in compute_sigmasq in SBank. I check this by computing the match of a normalized, whitened template with itself */
}
