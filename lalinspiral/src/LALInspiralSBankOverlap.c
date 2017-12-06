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
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>
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
            XLALDestroyCOMPLEX8FFTPlan(workspace_cache[k].plan);
            XLALDestroyCOMPLEX8Vector(workspace_cache[k].zf);
            XLALDestroyCOMPLEX8Vector(workspace_cache[k].zt);
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
    ptr->zf = XLALCreateCOMPLEX8Vector(n);
    CHECK_OOM(ptr->zf->data, "unable to allocate workspace array zf\n");
    memset(ptr->zf->data, 0, n * sizeof(COMPLEX8));

    ptr->zt = XLALCreateCOMPLEX8Vector(n);
    CHECK_OOM(ptr->zf->data, "unable to allocate workspace array zt\n");
    memset(ptr->zt->data, 0, n * sizeof(COMPLEX8));

    ptr->n = n;
    ptr->plan = XLALCreateReverseCOMPLEX8FFTPlan(n, 1);
    CHECK_OOM(ptr->plan, "unable to allocate plan");

    return ptr;
}

/* by default, complex arithmetic will call built-in function __muldc3, which does a lot of error checking for inf and nan; just do it manually */
static void multiply_conjugate(COMPLEX8 * restrict out, COMPLEX8 *a, COMPLEX8 *b, const size_t size) {
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

static double abs2(const COMPLEX8 x) {
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
    multiply_conjugate(ws->zf->data, inj->data->data, tmplt->data->data, min_len);
    XLALCOMPLEX8VectorFFT(ws->zt, ws->zf, ws->plan); /* plan is reverse */

    /* maximize over |z(t)|^2 */
    COMPLEX8 *zdata = ws->zt->data;
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

    /* compute match */
    /* return 4. * inj->deltaF * sqrt(result) / n; */  /* inverse FFT = reverse / n */
    return 4. * inj->deltaF * sqrt(result);  /* Ajith, 2012-11-09: The division by n is inconsistent with the revised convention of the normalization in compute_sigmasq in SBank. I check this by computing the match of a normalized, whitened template with itself */
}

REAL8 XLALInspiralSBankComputeMatchMaxSkyLoc(const COMPLEX8FrequencySeries *hp, const COMPLEX8FrequencySeries *hc, const REAL8 hpsigmasq, const REAL8 hcsigmasq, const REAL8 hphccorr, const COMPLEX8FrequencySeries *proposal, WS *workspace_cache1, WS *workspace_cache2) {
    /* FIXME: Add sanity checking for consistency of lengths in input */

    /* What does this do? */
    size_t min_len = (hp->data->length <= proposal->data->length) ? hp->data->length : proposal->data->length;

    /* get workspace for + and - frequencies */
    size_t n = 2 * (min_len - 1);   /* no need to integrate implicit zeros */
    WS *ws1 = get_workspace(workspace_cache1, n);
    if (!ws1) {
        XLALPrintError("out of space in the workspace_cache\n");
        XLAL_ERROR_REAL8(XLAL_ENOMEM);
    }
    WS *ws2 = get_workspace(workspace_cache2, n);
    if (!ws2) {
        XLALPrintError("out of space in the workspace_cache\n");
        XLAL_ERROR_REAL8(XLAL_ENOMEM);
    }


    /* compute complex SNR time-series in freq-domain, then time-domain */
    /* Note that findchirp paper eq 4.2 defines a positive-frequency integral,
       so we should only fill the positive frequencies (first half of zf). */
    multiply_conjugate(ws1->zf->data, hp->data->data, proposal->data->data, min_len);
    XLALCOMPLEX8VectorFFT(ws1->zt, ws1->zf, ws1->plan); /* plan is reverse */
    multiply_conjugate(ws2->zf->data, hc->data->data, proposal->data->data, min_len);
    XLALCOMPLEX8VectorFFT(ws2->zt, ws2->zf, ws2->plan);


    /* COMPUTE DETECTION STATISTIC */

    /* First start with constant values */
    REAL8 delta = 2 * (hphccorr) / hpsigmasq;
    REAL8 gamma = hcsigmasq / hpsigmasq;
    REAL8 denom = 4 * gamma - delta * delta;
    if (denom < 0)
    {
        fprintf(stderr, "DANGER WILL ROBINSON: CODE IS BROKEN!!\n");
    }

    /* Now the tricksy bit as we loop over time*/
    COMPLEX8 *hpdata = ws1->zt->data;
    COMPLEX8 *hcdata = ws2->zt->data;
    size_t k = n;
    /* FIXME: This is needed if we turn back on peak refinement. */
    /*ssize_t argmax = -1;*/
    REAL8 max = 0.;
    for (;k--;) {
        COMPLEX8 ratio = hcdata[k] / hpdata[k];
        REAL8 ratio_real = creal(ratio);
        REAL8 ratio_imag = cimag(ratio);
        REAL8 beta = 2 * ratio_real;
        REAL8 alpha = ratio_real * ratio_real + ratio_imag * ratio_imag;
        REAL8 sqroot = alpha*alpha + alpha * (delta*delta - 2 * gamma) + gamma*gamma;
        sqroot += beta * (beta * gamma - delta * (gamma + alpha));
        sqroot = sqrt(sqroot);
        REAL8 brckt = 2*alpha + 2*gamma - beta*delta + 2*sqroot;
        brckt = brckt / denom;
        REAL8 det_stat_sq = abs2(hpdata[k]) * brckt / hpsigmasq;

        if (det_stat_sq > max) {
            /*argmax = k;*/
            max = det_stat_sq;
        }
    }
    if (max == 0.) return 0.;

    /* FIXME: For now do *not* refine estimate of peak. */
    /* REAL8 result;
    if (argmax == 0 || argmax == (ssize_t) n - 1)
        result = max;
    else
        result = vector_peak_interp(abs2(zdata[argmax - 1]), abs2(zdata[argmax]), abs2(zdata[argmax + 1])); */

    /* Return match */
    return 4. * proposal->deltaF * sqrt(max);
}
