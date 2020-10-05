
/*
 * This file is part of TEOBResumS
 *
 * Copyright (C) 2017-2018  Alessandro Nagar, Sebastiano Bernuzzi,
 * Sarp Ackay, Gregorio Carullo, Walter Del Pozzo, Ka Wa Tsang, Michalis Agathos
 * LALSimulation implementation by Michalis Agathos
 *
 * TEOBResumS is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * TEOBResumS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see http://www.gnu.org/licenses/.
 *
 */

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

//#ifndef _OPENMP
//#define omp ignore
//#endif

#include <complex.h>
#include <math.h>

#include <gsl/gsl_const.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>

#include "LALSimTEOBResumS.h"

#include <stdio.h>


/* Return symm mass ratio from q */
REAL8 q_to_nu(const REAL8 q)
{
   REAL8 nu = 0;
   XLAL_CHECK_REAL8(q>0.0, XLAL_EDOM, "Mass ratio q cannot be negative!");
   if (q>0.)
       nu = q/((q+1.)*(q+1.));
   return nu;
}

/* Return mass ratio M1/M from nu */
REAL8 nu_to_X1(const REAL8 nu)
{
    return 0.5*(1.+sqrt(1.-4.*nu));
}

/* Eulerlog function (constants are defined in header) */
static const REAL8 Logm[] = {0.,Log1,Log2,Log3,Log4,Log5,Log6,Log7};
REAL8 Eulerlog(const REAL8 x,const INT4 m)
{
    REAL8 logm = 0.;
    if ((m>0) & (m<8)) logm = Logm[m];
    else logm = log((REAL8)m);
    return LAL_GAMMA + LAL_LN2 + logm + 0.5*log(x);
}

/* Spline interpolation with GSL routines */
void interp_spline(double *t, double *y, int n, double *ti, int ni, REAL8 *yi)
{
  //#pragma omp parallel
  //{
        gsl_interp_accel *acc = gsl_interp_accel_alloc ();
        gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n);
        gsl_spline_init (spline, t, y, n);

  //#pragma omp for
        for (INT4 k = 0; k < ni; k++)
        {
            yi[k] = (REAL8) gsl_spline_eval (spline, ti[k], acc);
        }
        gsl_spline_free (spline);
        gsl_interp_accel_free (acc);
  //}
}

/* Find nearest point index in 1d array */
INT4 find_point_bisection(REAL8 x, INT4 n, REAL8 *xp, INT4 o)
{
    INT4 i0 = o-1, i1 = n-o;
    INT4 i;
    if (n < 2*o)
    {
        XLAL_ERROR(XLAL_EDOM, "not enough points to interpolate");
    }
    if (x <= xp[i0]) return 0;
    if (x >  xp[i1]) return n-2*o;
    while (i0 != i1-1)
    {
        i = (i0+i1)/2;
        if (x < xp[i]) i1 = i; else i0 = i;
    }
    return i0-o+1;
}

/* Barycentric Lagrange interpolation at xx with n points of f(x),
 equivalent to standard Lagrangian interpolation */
#define tiny 1e-12
REAL8 baryc_f(REAL8 xx, INT4 n, REAL8 *f, REAL8 *x)
{
    REAL8 omega[n];
    REAL8 o, num, den, div, ci;
    INT4 i, j;
    for (i = 0; i < n; i++)
    {
        if (fabs(xx - x[i]) <= tiny) return f[i];
        o = 1.;
        for (j = 0; j < n; j++)
        {
            if (j != i)
            {
                o /= (x[i] - x[j]);
            }
        }
        omega[i] = o;
    }
    num = den = 0.;
    for (i = 0; i < n; i++)
    {
        div  = xx - x[i];
        ci   = omega[i]/div;
        den += ci;
        num += ci * f[i];
    }
    return( num/den );
}

/* Barycentric Lagrange interpolation at xx with n points of f(x),
 compute weights */
void baryc_weights(INT4 n, REAL8 *x, REAL8 *omega)
{
    REAL8 o;
    INT4 i, j;
    for (i = 0; i < n; i++)
    {
        o = 1.;
        for (j = 0; j < n; j++)
        {
            if (j != i)
            {
                o /= (x[i] - x[j]);
            }
        }
        omega[i] = o;
    }
}

/* Barycentric Lagrange interpolation at xx with n points of f(x),
 use precomputed weights */
REAL8 baryc_f_weights(REAL8 xx, INT4 n, REAL8 *f, REAL8 *x, REAL8 *omega)
{
    INT4 i;
    REAL8 num, den, div, ci;
    num = den = 0.;
    for (i = 0; i < n; i++)
    {
        div  = xx - x[i];
        if (fabs(div) <= tiny) return f[i];
        ci   = omega[i]/div;
        den += ci;
        num += ci * f[i];
    }
    return( num/den );
}

/* 1d Lagrangian barycentric interpolation */
REAL8 interp1d (const INT4 order, REAL8 xx, INT4 nx, REAL8 *f, REAL8 *x)
{
    REAL8 ff;
    INT4 ix;
    INT4 ox = order > nx ? nx : order;
    ix = find_point_bisection(xx, nx, x, ox/2);
    ff = baryc_f(xx, ox, &f[ix], &x[ix]);
    return( ff );
}

/* Find max location by poynomial interpolation around x0 (uniform grid) */
REAL8 find_max (const INT4 n, REAL8 dx, REAL8 x0, REAL8 *f, REAL8 *fmax)
{
    const INT4 i = (n-1)/2; /* To centre the grid for n = N, e.g., n=7, i-3 = 0. */
    REAL8 xmax = x0;
    REAL8 d1f = 0., d2f = 0.;
    if (n==3)
    {
        d1f = 0.5*(f[i+1]-f[i-1]);
        d2f = (f[i-1]-2*f[i]+f[i+1]);
    }
    else if (n==5)
    {
        d1f = (8.*(f[i+1]-f[i-1]) - f[i+2] + f[i-2]);
        d2f = (-30*f[i]+16*(f[i+1]+f[i-1])-(f[i+2]+f[i-2]));
    }
    else if (n==7)
    {
        d1f = ( 45.0*(f[i+1]-f[i-1]) - 9.0*(f[i+2] - f[i-2]) + f[i+3] - f[i-3] )/(60.0);
        d2f = ( -490.0*f[i]+270.0*(f[i+1]+f[i-1])-27.0*(f[i+2]+f[i-2])+2.0*(f[i+3]+f[i-3]) )/(180.0);
    }
    else XLAL_ERROR_REAL8(XLAL_EINVAL, "Implemented only n = 3,5,7");

    if (d2f != 0.)
    {
        xmax -= dx*d1f/d2f;
    }
    /* Eval function */
    if (fmax!=NULL)
    {
        if (n==3)
        {
            *fmax = ( -((dx + x0 - xmax)*(-x0 + xmax)*f[-1 + i])
                     + (dx - x0 + xmax)*(2*(dx + x0 - xmax)*f[i] + (-x0 + xmax)*f[1 + i]) );
            *fmax /= (2.*SQ(dx));
        }
        else if (n==5)
        {
            *fmax = ((dx + x0 - xmax)*(2*dx + x0 - xmax)*(-x0 + xmax)*(dx - x0 + xmax)*f[-2 + i]
                     + (2*dx - x0 + xmax)*(-4*(dx + x0 - xmax)*(2*dx + x0 - xmax)*(-x0 + xmax)*f[-1 + i]+(dx - x0 + xmax)*(6*(dx + x0 - xmax)*(2*dx + x0 - xmax)*f[i] + (-x0 + xmax)*(4*(2*dx + x0 - xmax)*f[1 + i]- (dx + x0 - xmax)*f[2 + i]))));
            *fmax /= (24.*pow(dx,4));

        }
        else
        {
          XLAL_ERROR_REAL8(XLAL_EINVAL, "Implemented only n = 3,5");
        }
    }

    return xmax;
}

///* Factorial */
//double fact(int n)
//{
//    double f[] = {1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880.,
//        3628800., 39916800., 479001600., 6227020800., 87178291200.};
//    if (n < 0){
//        errorexit(" computing a negative factorial.\n");
//    } else if (n <= 14){
//        return f[n];
//    } else {
//        return n*fact(n-1);
//    }
//}

///* Wigner d-function */
//double wigner_d_function(int l, int m, int s, double i)
//{
//    double dWig = 0.;
//    double costheta = cos(i*0.5);
//    double sintheta = sin(i*0.5);
//    int ki = MAX( 0  , m-s );
//    int kf = MIN( l+m, l-s );
//    int k;
//    for( k = ki; k <= kf; k++ ){
//        dWig +=
//        ( pow(-1.,k) * pow(costheta,2*l+m-s-2*k) * pow(sintheta,2*k+s-m) )/
//        ( fact(k) * fact(l+m-k) * fact(l-s-k) * fact(s-m+k) );
//    }
//    return (sqrt(fact(l+m) * fact(l-m) * fact(l+s) * fact(l-s)) * dWig);
//}

///* Spin-weighted spherical harmonic
// Ref: https://arxiv.org/pdf/0709.0093.pdf */
//int spinsphericalharm(double *rY, double *iY, int s, int l, int m, double phi, double i)
//{
//    if ((l<0) || (m<-l) || (m>l)) {
//        errorexit(" wrong (l,m) inside spinspharmY\n");
//    }
//    double c = pow(-1.,-s) * sqrt( (2.*l+1.)/(4.*M_PI) );
//    double dWigner = c * wigner_d_function(l,m,-s,i);
//    *rY = cos((double)(m)*phi) * dWigner;
//    *iY = sin((double)(m)*phi) * dWigner;
//    return XLAL_SUCCESS;
//}

/*
 * ModeArray is a structure which allows to select the modes to include
 * in the waveform.
 * This function will create a structure with the given selection of modes
 * as enumerated in LALSimTEOBResumS.h.
 * Default values for different types of binaries are hardcoded as macros.
 */
INT4 XLALSetup_TEOB_mode_array(LALValue *ModeArray, INT4 modeType)
{
    switch (modeType)
    {
        /* Only include (2,2) and (2,-2) mode */
        case TEOB_MODES_22:
            XLALSimInspiralModeArrayDeactivateAllModes(ModeArray);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 2);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -2);
            break;
        /* For BNS all modes up to l=8 are available */
        case TEOB_MODES_ALL:
            XLALSimInspiralModeArrayActivateAllModes(ModeArray);
            break;

        default:
            XLAL_ERROR(XLAL_EINVAL, "Unknown TEOB mode-list code.");
            break;
    }

    return XLAL_SUCCESS;
}

/*
 * ModeArray is a structure which allows to select the modes to include.
 * This function checks if the selected modes are available for a given type
 * of binary with TEOBResumS. It also checks that the symmetric modes are
 * consistently active.
 */
INT4 XLALCheck_TEOB_mode_array(LALValue *ModeArray,
                                      UINT4 use_tidal)
{
    /* If one select BNS all the modes up to (8,8) are available */
    /* If not, only the (2,2) mode is available */

    /*Loop over all the possible modes
     *we only check +m modes, when one select (l,m) mode is actually
     *selecting (l,|m|) mode
     */
    for (UINT4 ELL = 2; ELL <= LAL_SIM_L_MAX_MODE_ARRAY; ELL++)
    {
        for (INT4 EMM = 0; EMM <= (INT4) ELL; EMM++)
        {
            if (XLALSimInspiralModeArrayIsModeActive(ModeArray, ELL, EMM))
            {
                /* Make sure that no modes exist with l>8 */
                XLAL_CHECK(ELL<=8, XLAL_ENOSYS, "Modes with l>8 are not available in TEOBResumS.\n");

                /* Make sure that for non-BNS systems, only \f$(2,\pm 2)\f$ modes are active */
                if (!use_tidal) XLAL_CHECK((ELL==2) && (abs(EMM)==2), XLAL_EDOM, "Modes beyond (2,\\pm 2) are only available for BNS in TEOBResumS.");

                /* Make sure (l,-m) mode is activated if (l,m) is active. Does not
                 * exit with error but silently adds the symmetric mode if not present.
                 * Bitwise operation in implementation ensures that no check is necessary.
                 */
                XLALSimInspiralModeArrayActivateMode(ModeArray, ELL, -EMM);
            }
            else if (XLALSimInspiralModeArrayIsModeActive(ModeArray, ELL, -EMM))
            {
                /* Make sure (l,m) mode is active if (l,-m) is active */
                XLAL_CHECK(XLALSimInspiralModeArrayIsModeActive(ModeArray, ELL, EMM), XLAL_EDOM, "Symmetric (l,-m) mode cannot be included without the (l,m) mode being active.");
            }
        }
    }

    return XLAL_SUCCESS;
}

/* (h+, hx) polarizations from the multipolar waveform in the source frame */
void XLALSimIMRComputePolarisations(REAL8Sequence *hplus_out,
                                    REAL8Sequence *hcross_out,
                                    SphHarmTimeSeries *hlm,
                                    LALValue *modeArray,
                                    REAL8 amplitude_prefactor,
                                    REAL8 theta,
                                    REAL8 phi)
{
    INT4 mneg = 1; /* m>0 modes only, add m<0 modes afterwards */
    COMPLEX16 Y, Ym=0.0, hpc;
    UINT4 l;
    INT4 m;
    SphHarmTimeSeries *this_mode = hlm;

    while (this_mode)
    {
        l = this_mode->l;
        m = this_mode->m;
        if (!XLALSimInspiralModeArrayIsModeActive(modeArray, l, m))
        {
            this_mode = this_mode->next;
            continue;
        }
        mneg = XLALSimInspiralModeArrayIsModeActive(modeArray, l, - m);

        Y = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, m);
        if ( (mneg) && (m != 0) )
        {
            Ym = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, - m);
            if ( l % 2 ) Ym = -Ym; /* l is odd */
        }

        for (UINT4 i=0; i < this_mode->mode->data->length; i++)
        {
            hpc=0.0;
            hpc += this_mode->mode->data->data[i] * Y;
            /* add m<0 modes */
            if ( (mneg) && (m != 0) )
                hpc += conj(this_mode->mode->data->data[i]) * Ym;

            hpc *= amplitude_prefactor;
            hplus_out->data[i]  += creal(hpc);
            hcross_out->data[i] -= cimag(hpc);
        }

        this_mode = this_mode->next;
    }

    return;
}

/* (h+, hx) polarizations from the multipolar waveform in the source frame;
    polar representation of time series */
void XLALSimIMRComputePolarisationsPolar(REAL8Sequence *hplus_out,
                                         REAL8Sequence *hcross_out,
                                         SphHarmPolarTimeSeries *hlm,
                                         LALValue *modeArray,
                                         REAL8 amplitude_prefactor,
                                         REAL8 theta,
                                         REAL8 phi)
{
    INT4 mneg = 1; /* m>0 modes only, add m<0 modes afterwards */
    COMPLEX16 Y, Ym = 0.0, hpc;
    UINT4 l;
    INT4 m;
    SphHarmPolarTimeSeries *this_mode = hlm;

    while (this_mode)
    {
        l = this_mode->l;
        m = this_mode->m;
        if (!XLALSimInspiralModeArrayIsModeActive(modeArray, l, m))
        {
            this_mode = this_mode->next;
            continue;
        }
        mneg = XLALSimInspiralModeArrayIsModeActive(modeArray, l, - m);
        //#pragma omp parallel
        //{
            /* Convention (master) theta = (-)incl and phi = pi/2 - phiRef */
            Y = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, m);
            if ( (mneg) && (m != 0) )
            {
                Ym = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, -m);
                if ( l % 2 ) Ym = -Ym; /* l is odd */
            }
            //#pragma omp for
            for (UINT4 i=0; i < this_mode->ampl->data->length; i++)
            {
                hpc = cpolar(this_mode->ampl->data->data[i], -this_mode->phase->data->data[i]) * Y;

                /* add m<0 modes */
                if ( (mneg) && (m != 0) )
                    hpc += cpolar(this_mode->ampl->data->data[i], this_mode->phase->data->data[i]) * Ym;

                hpc *= amplitude_prefactor;
                hplus_out->data[i] += creal(hpc);
                hcross_out->data[i] -= cimag(hpc);
            }
      //}
      this_mode = this_mode->next;
    }

    return;
}


/* 4th order centered stencil first derivative, uniform grids */
INT4 D0(REAL8 *f, REAL8 dx, INT4 n, REAL8 *df)
{
    const REAL8 oo12dx  = 1./(12*dx);
    INT4 i;
    for (i=2; i<n-2; i++)
    {
        df[i] = (8.*(f[i+1]-f[i-1]) - f[i+2] + f[i-2])*oo12dx;
    }
    i = 0;
    df[i] = (-25.*f[i] + 48.*f[i+1] - 36.*f[i+2] + 16.*f[i+3] - 3.*f[i+4])*oo12dx;
    i = 1;
    df[i] = (-3.*f[i-1] - 10.*f[i] + 18.*f[i+1] - 6.*f[i+2] + f[i+3])*oo12dx;
    i = n-2;
    df[i] = - (-3.*f[i+1] - 10.*f[i] + 18.*f[i-1] - 6.*f[i-2] + f[i-3])*oo12dx;
    i = n-1;
    df[i] = - (-25.*f[i] + 48.*f[i-1] - 36.*f[i-2] + 16.*f[i-3] - 3.*f[i-4])*oo12dx;
    return XLAL_SUCCESS;
}

/* 4th order centered stencil second derivative, uniform grids */
INT4 D2(REAL8 *f, REAL8 dx, INT4 n, REAL8 *d2f)
{
    const REAL8 oo12dx2  = 1./(dx*dx*12);
    INT4 i;
    for (i=2; i<n-2; i++)
    {
        d2f[i] = (-30*f[i]+16*(f[i+1]+f[i-1])-(f[i+2]+f[i-2]))*oo12dx2;
    }
    i= 0;
    d2f[i] = (45*f[i]-154*f[i+1]+214*f[i+2]-156*f[i+3]+61*f[i+4]-10*f[i+5])*oo12dx2;
    i= 1;
    d2f[i] = (10*f[i-1]-15*f[i]-4*f[i+1]+14*f[i+2]-6*f[i+3]+f[i+4])*oo12dx2;
    i = n-2;
    d2f[i] = (10*f[i+1]-15*f[i]-4*f[i-1]+14*f[i-2]-6*f[i-3]+f[i-4])*oo12dx2;
    i = n-1;
    d2f[i] = (45*f[i]-154*f[i-1]+214*f[i-2]-156*f[i-3]+61*f[i-4]-10*f[i-5])*oo12dx2;
    return XLAL_SUCCESS;
}

/* Third-order polynomial integration */
REAL8 cumint3(REAL8 *f, REAL8 *x, const INT4 n, REAL8 *sum)
{
    REAL8 xe[n+2], fe[n+2];
    REAL8 *x0,*x1,*x2,*x3;
    REAL8 *f0,*f1,*f2,*f3;
    REAL8 a,b,c,d,e,h,g,z;
    const REAL8 oo12 = 0.08333333333333333;

    for (INT4 i=1; i < n+1; i++)
    {
        xe[i] = x[i-1];
        fe[i] = f[i-1];
    }
    xe[0]   = x[3];
    xe[n+1] = x[n-4];
    fe[0]   = f[3];
    fe[n+1] = f[n-4];

    x0 = &xe[0];
    x1 = &xe[1];
    x2 = &xe[2];
    x3 = &xe[3];

    f0 = &fe[0];
    f1 = &fe[1];
    f2 = &fe[2];
    f3 = &fe[3];

    sum[0] = 0.;
    for (INT4 i=0; i < n-1; i++)
    {
        a = x1[i]-x0[i];
        b = x2[i]-x1[i];
        c = x3[i]-x2[i];
        d = f1[i]-f0[i];
        e = f2[i]-f1[i];
        h = f3[i]-f2[i];
        g = 0.5*(f1[i]+f2[i]);
        z = b*g + oo12*b*b*(c*b*(2*c+b)*(c+b)*d-a*c*(c-a)*(2*c+2*a+3*b)*e-a*b*(2*a+b)*(a+b)*h)/(a*c*(a+b)*(c+b)*(c+a+b));
        sum[i+1] = sum[i] + z;
    }

    return sum[n-1];
}

/* Simple unwrap for phase angles */ // Do NOT mess with it
void unwrap(REAL8 *p, const INT4 size)
{
    if (size < 1) return;
    INT4 j;
    INT4 fact = 0;  // For making the initial negative phase angles positive
    REAL8 curr, prev;
    REAL8 corr = 0.0;
    REAL8 dphi = 0.0;

    prev = p[0];
    if( p[0] < 0 ) fact = 1;
    if( p[1] < p[0] )
        dphi = LAL_TWOPI;

    for (j = 1; j < size; j++)
    {
        p[j] += fact*LAL_TWOPI;
        curr = p[j];
        if( curr < prev )
            dphi = LAL_TWOPI;
        corr += dphi;
        p[j] += corr - fact*LAL_TWOPI;
        prev = curr;
        dphi = 0.0;
    }
}

/* Unwrap unsign number of cycles from reference phase as proxy */
void unwrap_proxy(REAL8 *p, REAL8 *r, const INT4 size, const INT4 shift0)
{
    if (size < 1) return;

    INT4 *n = (INT4 *) calloc(size, sizeof(INT4));
    XLAL_CHECK_VOID(n, XLAL_ENOMEM, "Out of memory");

    const REAL8 ooTwoPi = 1.0/(LAL_TWOPI);
    const REAL8 r0 = r[0];

    /* no cycles from reference phase, r>0 */
    for (INT4 i = 0; i < size; i++)
        n[i] = floor ( (r[i] - r0) * ooTwoPi );

    if (shift0)
    {
        /* shift phase : p(0) = r(0) */
        const REAL8 dp = (r0 - p[0]);
        for (INT4 i = 0; i < size; i++)
            p[i] += dp;
    }

    /* unwrap based on no cycles */
    const REAL8 p0 = p[0];
    INT4 np = 0;
    for (INT4 i = 0; i < size; i++)
    {
        p[i] += n[i]*LAL_TWOPI;
        /* correct cases p = - Pi + epsilon */
        np = floor ( ( p[i] - p0 )*ooTwoPi );
        p[i] += (n[i]-np)*LAL_TWOPI;
    }

    XLALFree(n);
}

/* Compute size of a uniform grid t0:dt:tf */
INT4 get_uniform_size(const REAL8 tN, const REAL8 t0, const REAL8 dt)
{
    return ((INT4)((tN - t0)/dt + 1));
}

/* Compute optimized timestep after merger */
REAL8 get_mrg_timestep(REAL8 q, REAL8 chi1, REAL8 chi2)
{
    REAL8 dt = q+chi1+chi2; // MA: just to get rid of UNUSED errors
    dt = 0.1;
    // ...
    return dt;
}

/* Multipolar waveform at time point (complex) */
void Waveform_lm_t_alloc (LALTEOBResumSWaveformModeSingleTime **wav)
{
    *wav = (LALTEOBResumSWaveformModeSingleTime *) calloc(1, sizeof(LALTEOBResumSWaveformModeSingleTime));
    XLAL_CHECK_VOID(wav, XLAL_ENOMEM, "Could not allocate mode instance.\n");
    (*wav)->time = 0.;

    return;
}

void Waveform_lm_t_free (LALTEOBResumSWaveformModeSingleTime *wav)
{
    free(wav);
    return;
}

/* NQC data */
void NQCdata_alloc (NQCdata **nqc)
{
  *nqc = (NQCdata *) calloc(1, sizeof(NQCdata));
  XLAL_CHECK_VOID(nqc, XLAL_ENOMEM, "Out of memory");
  (*nqc)->flx = (NQCcoefs *) calloc(1, sizeof(NQCcoefs));
  XLAL_CHECK_VOID((*nqc)->flx, XLAL_ENOMEM, "Out of memory");
  (*nqc)->hlm = (NQCcoefs *) calloc(1, sizeof(NQCcoefs));
  XLAL_CHECK_VOID((*nqc)->hlm, XLAL_ENOMEM, "Out of memory");
}

void NQCdata_free (NQCdata *nqc)
{
  if(nqc->flx) free (nqc->flx);
  if(nqc->hlm) free (nqc->hlm);
  if(nqc) free (nqc);
}


/* Dynamics */
void XLALTEOBDynamicsInit (LALTEOBResumSDynamics **dyn, INT4 size, const char *name)
{
    (*dyn) = calloc(1, sizeof(LALTEOBResumSDynamics));
    XLAL_CHECK_VOID(dyn, XLAL_ENOMEM, "Could not allocate TEOB Dynamics.\n");

    strcpy((*dyn)->name,name);
    (*dyn)->size = size;
    (*dyn)->time = calloc ((size_t) size, sizeof(REAL8) );   // TODO: speedup: calloc < malloc+memset

    for (INT4 v = 0; v < TEOB_DYNAMICS_NVARS; v++)
    {
        (*dyn)->data[v] = calloc ((size_t) size, sizeof(REAL8) );
    }
    NQCdata_alloc(&(*dyn)->NQC);
    return;
}

void XLALTEOBDynamicsPush (LALTEOBResumSDynamics **dyn, INT4 size)
{
    (*dyn)->time = realloc( (*dyn)->time, size * sizeof(REAL8) );
    for (INT4 v = 0; v < TEOB_DYNAMICS_NVARS; v++)
    {
        (*dyn)->data[v] = realloc ( (*dyn)->data[v], size * sizeof(REAL8) );
        if ((*dyn)->data[v] == NULL) XLAL_ERROR_VOID(XLAL_ENOMEM, "Could not allocate TEOB Dynamics.\n");
        /* if (dn>0) memset( (*dyn)->data[v] + n, 0, dn * sizeof(REAL8) ); */ // TODO: MA: use this?
    }
    (*dyn)->size = size;
}

/* Interp and overwrite a multipolar waveform */
void XLALTEOBDynamicsInterp (LALTEOBResumSDynamics *dyn, const INT4 size, const REAL8 t0, const REAL8 dt, const char *name)
{
    /* Alloc and init aux memory */
    LALTEOBResumSDynamics *dyn_aux;
    const INT4 oldsize = dyn->size;
    XLALTEOBDynamicsInit(&dyn_aux, oldsize, "");
    memcpy(dyn_aux->time, dyn->time, oldsize * sizeof(REAL8));
    for (int v = 0; v < TEOB_DYNAMICS_NVARS; v++)
        memcpy(dyn_aux->data[v], dyn->data[v], oldsize * sizeof(REAL8));

    /* Overwrite and realloc arrays */
    dyn->dt   = dt;
    dyn->size = size;

    if (strcmp(name, "")) strcpy(dyn->name, name);
    if(dyn->time) free(dyn->time);

    dyn->time = malloc ( size * sizeof(REAL8) );

    for (int v = 0; v < TEOB_DYNAMICS_NVARS; v++)
    {
        if (dyn->data[v]) free(dyn->data[v]);
        dyn->data[v] = malloc ( size * sizeof(REAL8) );
        /* memset(dyn->data[v], 0., size*sizeof(double)); */
    }

    /* Fill new time array */
    for (int i = 0; i < size; i++)
        dyn->time[i] = i*dt + t0;

    /* Interp */
    for (int k = 0; k < TEOB_DYNAMICS_NVARS; k++)
        interp_spline(dyn_aux->time, dyn_aux->data[k], dyn_aux->size, dyn->time, size, dyn->data[k]);

    /* Free aux memory */
    XLALFreeTEOBDynamics (dyn_aux);
}

/* Extract Dynamics at times t >= to and t < tn, Alloc a new Dynamics var */
void XLALTEOBDynamicsExtract (LALTEOBResumSDynamics *dyna, const REAL8 to, const REAL8 tn, LALTEOBResumSDynamics **dynb, const char *name)
{
    /* Check limits */
    if (tn<to)
        XLAL_ERROR_VOID(XLAL_EINVAL, "Bad choice of times: tn < to");
    if (to > dyna->time[dyna->size-1])
        XLAL_ERROR_VOID(XLAL_EINVAL, "Nothing to extract, to > time[size-1]");
    if (tn < dyna->time[0])
        XLAL_ERROR_VOID(XLAL_EINVAL, "Nothing to extract, tn < time[0]");

    /* Find indexes of closest elements to  (to, tn) */
    INT4 io = 0;
    INT4 in = dyna->size-1;
    if (to > dyna->time[0])
        io = find_point_bisection(to, dyna->size, dyna->time, 1);
    if (tn < dyna->time[dyna->size-1])
        in = find_point_bisection(tn, dyna->size, dyna->time, 0);

    /* Calculate the new size */
    const INT4 N  = in-io;


    /* Alloc output waveform b */
    XLALTEOBDynamicsInit (dynb, N, name);

    /* TODO: Parameters are not copied in the new wf !*/
    /*
     (*dynb) = (Dynamics *) calloc(1, sizeof(Dynamics));
     if (dynb == NULL)
     errorexit("Out of memory");
     strcpy((*dynb)->name,name);
     (*dynb)->size = N;
     memcpy(*dynb, dyna, sizeof(Dynamics)); // copy parameters
     (*dynb)->time = malloc ( N * sizeof(double) );
     for (int v = 0; v < TEOB_DYNAMICS_NVARS; v++)
     (*dynb)->data[v] = malloc ( N * sizeof(double) );
     */

    /* Copy the relevant part of a into b */
    for (int i = 0; i < N; i++)
        (*dynb)->time[i] = dyna->time[io + i];
    for (int v = 0; v < TEOB_DYNAMICS_NVARS; v++)
    {
        for (int i = 0; i < N; i++)
        {
            (*dynb)->data[v][i] = dyna->data[v][io + i];
        }
    }

}

/* Join two dynamics time series at t = to */
void XLALTEOBDynamicsJoin (LALTEOBResumSDynamics *dyna, LALTEOBResumSDynamics *dynb, REAL8 to)
{
    /* Time arrays are suppose to be ordered as
     dyna->time:  x x x x x x x x x
     dynb->time:       o o o o o o o o o
     to        :                |
     But they do not need to overlap or be uniformly spaced.
     Note to can be
     to > dyna->time[hlma->size-1] => extend the dynamics data
     to < dynb->time[0]            => join the whole b dynamics
     Following checks enforce the above structure, if possible.
     */
    if (dyna->time[0] > dynb->time[0]) SWAPTRS( dyna, dynb );
    XLAL_CHECK_VOID (to <= dynb->time[dynb->size-1], XLAL_EINVAL, "Joining time outside range. Dynamics not joined.");
    XLAL_CHECK_VOID (to > dyna->time[0], XLAL_EINVAL, "Joining time outside range. Dynamics not joined.");

    /* Find indexes of closest elements to to */
    const INT4 ioa = find_point_bisection(to, dyna->size, dyna->time, 1);
    INT4 iob = find_point_bisection(to, dynb->size, dynb->time, 1);
    if ( DEQUAL(dyna->time[ioa], dynb->time[iob], 1e-10) ) iob++;

    /* Calculate the new size */
    const INT4 Nb = dynb->size - iob;
    const INT4 N  = ioa + Nb;

    /* Resize a */
    XLALTEOBDynamicsPush (&dyna, N);

    /* Copy the relevant part of b into a */
    for (int i = 0; i < Nb; i++)
        dyna->time[ioa + i] = dynb->time[iob + i];
    for (int v = 0; v < TEOB_DYNAMICS_NVARS; v++)
    {
        for (int i = 0; i < Nb; i++)
        {
            dyna->data[v][ioa + i] = dynb->data[v][iob + i];
        }
    }
}


void XLALFreeTEOBDynamics (LALTEOBResumSDynamics *dyn)
{
    if(dyn->time) free(dyn->time);
    for (INT4 v = 0; v < TEOB_DYNAMICS_NVARS; v++)
        if (dyn->data[v]) free(dyn->data[v]);
    NQCdata_free(dyn->NQC);
    free(dyn);
    return;
}


/* Sync some quick access parameters in dyn with parameter database
 to be used carefully */
void XLALTEOBDynamicsSetParams(LALTEOBResumSDynamics *dyn,                 /*< pointer to LALTEOBResumSDynamics struct to be filled in */
                               UNUSED LALValue **pModeArray,                /*< Mode array to be returned */
                               const REAL8 m1,                             /*< mass of companion 1 (kg) */
                               const REAL8 m2,                             /*< mass of companion 2 (kg) */
                               const REAL8 UNUSED S1x,                            /*< x-component of the dimensionless spin of object 1 */
                               const REAL8 UNUSED S1y,                            /*< y-component of the dimensionless spin of object 1 */
                               const REAL8 S1z,                            /*< z-component of the dimensionless spin of object 1 */
                               const REAL8 UNUSED S2x,                            /*< x-component of the dimensionless spin of object 2 */
                               const REAL8 UNUSED S2y,                            /*< y-component of the dimensionless spin of object 2 */
                               const REAL8 S2z,                            /*< z-component of the dimensionless spin of object 2 */
                               const REAL8 dt,                         /*< sampling interval (s) */
                               const REAL8 lambda1,                    /*< l=2 tidal polarizability for m1 */
                               const REAL8 lambda2,                    /*< l=2 tidal polarizability for m2 */
                               const REAL8 lambda1oct,                    /*< l=3 tidal polarizability for m1 */
                               const REAL8 lambda2oct,                    /*< l=3 tidal polarizability for m2 */
                               const REAL8 lambda1hex,                    /*< l=4 tidal polarizability for m1 */
                               const REAL8 lambda2hex,                    /*< l=4 tidal polarizability for m2 */
                               LALDict *LALparams                  /*< LALDict dictionary for extra waveform options */
                               )
{
    LALSimInspiralSpinOrder spinO = XLALSimInspiralWaveformParamsLookupPNSpinOrder(LALparams);
    LALSimInspiralTidalOrder tidalO = XLALSimInspiralWaveformParamsLookupPNTidalOrder(LALparams);
    LALSimInspiralGETides geTides = TEOB_GE_TIDES_DEFAULT;
    LALSimInspiralGMTides gmTides = TEOB_GM_TIDES_DEFAULT;
    if (XLALDictContains(LALparams, "GEtideO"))
        geTides = XLALSimInspiralWaveformParamsLookupGETides(LALparams);
    if (XLALDictContains(LALparams, "GMtideO"))
        gmTides = XLALSimInspiralWaveformParamsLookupGMTides(LALparams);


    /* Set auxiliary parameters */
    dyn->dt = dt;
    dyn->M  = (m1+m2)/LAL_MSUN_SI; /* Msun */
    dyn->q  =  m1/m2;
    dyn->nu = q_to_nu(dyn->q);
    dyn->X1 = nu_to_X1(dyn->nu);
    dyn->X2 = 1. - dyn->X1;

    dyn->chi1  = S1z;
    dyn->chi2  = S2z;
    dyn->S1    = SQ(dyn->X1) * dyn->chi1;
    dyn->S2    = SQ(dyn->X2) * dyn->chi2;
    dyn->a1    = dyn->X1*dyn->chi1;
    dyn->a2    = dyn->X2*dyn->chi2;
    REAL8 aK   = dyn->a1 + dyn->a2;
    dyn->aK2   = SQ(aK);
    dyn->S     = dyn->S1 + dyn->S2;                     /* in the EMRL this becomes the spin of the BH */
    dyn->Sstar = dyn->X2*dyn->a1 + dyn->X1*dyn->a2;  /* in the EMRL this becomes the spin of the particle */

    INT4 ssorder;
    INT4 rc_order;
    REAL8 c3 = 0.;
    INT4 teobModeChoice = TEOB_MODES_NOPT; /* initialize to unavailable value */


    /* Reading tidal order, spin order and GM tides from LALDict */

    if (geTides > LAL_SIM_INSPIRAL_GETIDES_NOPT) XLAL_ERROR_VOID(XLAL_EDATA, "Gravitoelectric tides option not supported.");

    switch (tidalO)
    {
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL:
            dyn->use_tidal = geTides;
            break;
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
            dyn->use_tidal = LAL_SIM_INSPIRAL_GETIDES_OFF;
            break;
        default:
            XLAL_ERROR_VOID(XLAL_EDATA, "Tidal order not supported.");
    }

    /* Determine the type of binary based on lambda being 0 */
    dyn->bhns = 0;
    if (lambda1 == 0.0 && lambda2 ==0.0) dyn->use_tidal = LAL_SIM_INSPIRAL_GETIDES_OFF;
    else if (lambda1 == 0.0) dyn->bhns = 1;
    else if (lambda2 == 0.0) dyn->bhns = 2;

    INT4 use_Yagi_fits = 1;
    if (XLALDictContains(LALparams, "use_Yagi_fits"))
        use_Yagi_fits = XLALDictLookupINT4Value(LALparams, "use_Yagi_fits");

    if ((dyn->use_tidal != LAL_SIM_INSPIRAL_GETIDES_OFF) && (!use_Yagi_fits))
    {
        use_Yagi_fits = 0;
        XLAL_CHECK_VOID(lambda1 != 0 && lambda2 != 0 && lambda1oct != 0 && lambda2oct != 0 && lambda1hex != 0 && lambda2hex != 0, XLAL_EDOM, "If Yagi fits are not used, then all tidal parameters (quad, oct, hex) need to be positive!");
    }

    if (gmTides > LAL_SIM_INSPIRAL_GMTIDES_NOPT) XLAL_ERROR_VOID(XLAL_EDATA, "Gravitomagnetic tides option not supported.");
    dyn->use_tidal_gravitomagnetic = gmTides;

    switch (spinO)
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
            dyn->use_spins = 1;
            break;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            dyn->use_spins = 0;
            break;
        default:
            XLAL_ERROR_VOID(XLAL_EDATA, "Spin order not supported.");
    }


    if (dyn->use_tidal!=LAL_SIM_INSPIRAL_GETIDES_OFF)
    {
        REAL8 LambdaAl2 = lambda1;
        REAL8 LambdaBl2 = lambda2;
        REAL8 LambdaAl3, LambdaAl4, LambdaBl3, LambdaBl4, SigmaAl2, SigmaBl2;
        LambdaAl3 = lambda1oct;
        LambdaAl4 = lambda1hex;
        LambdaBl3 = lambda2oct;
        LambdaBl4 = lambda2hex;
        SigmaAl2 = SigmaBl2 = 0.0;

        if (use_Yagi_fits)
        {
            /* Allow for manual override of oct and hex Lambdas */
            if (LambdaAl3==0.0) LambdaAl3 = Yagi13_fit_barlamdel(LambdaAl2, 3);
            if (LambdaBl3==0.0) LambdaBl3 = Yagi13_fit_barlamdel(LambdaBl2, 3);
            if (LambdaAl4==0.0) LambdaAl4 = Yagi13_fit_barlamdel(LambdaAl2, 4);
            if (LambdaBl4==0.0) LambdaBl4 = Yagi13_fit_barlamdel(LambdaBl2, 4);

            /* Write new values to LALDict so that they are globally accessible */
            XLALSimInspiralWaveformParamsInsertTidalOctupolarLambda1(LALparams, LambdaAl3);
            XLALSimInspiralWaveformParamsInsertTidalOctupolarLambda2(LALparams, LambdaBl3);
            XLALSimInspiralWaveformParamsInsertTidalHexadecapolarLambda1(LALparams, LambdaAl4);
            XLALSimInspiralWaveformParamsInsertTidalHexadecapolarLambda2(LALparams, LambdaBl4);

            /* Old fits [ref] */
            // dyn->SigmaAl2  = Yagi13_fit_barsigmalambda(dyn->LambdaAl2);
            // dyn->SigmaBl2  = Yagi13_fit_barsigmalambda(dyn->LambdaBl2);

            /* New fits [ref] */
            SigmaAl2 = JFAPG_fit_Sigma_Irrotational(LambdaAl2);
            SigmaBl2 = JFAPG_fit_Sigma_Irrotational(LambdaBl2);
        }

        if (dyn->use_tidal_gravitomagnetic != LAL_SIM_INSPIRAL_GMTIDES_OFF)
        {
            SigmaAl2 = JFAPG_fit_Sigma_Irrotational(LambdaAl2);
            SigmaBl2 = JFAPG_fit_Sigma_Irrotational(LambdaBl2);
        }

        /* Tidal coupling constants */
        REAL8 XA    = dyn->X1;
        REAL8 XB    = dyn->X2;
        dyn->kapA2  = 3.   * LambdaAl2 * XA*XA*XA*XA*XA / dyn->q;
        dyn->kapA3  = 15.  * LambdaAl3 * XA*XA*XA*XA*XA*XA*XA / dyn->q;
        dyn->kapA4  = 105. * LambdaAl4 * XA*XA*XA*XA*XA*XA*XA*XA*XA / dyn->q;

        dyn->kapB2  = 3.   * LambdaBl2 * XB*XB*XB*XB*XB * dyn->q;
        dyn->kapB3  = 15.  * LambdaBl3 * XB*XB*XB*XB*XB*XB*XB * dyn->q;
        dyn->kapB4  = 105. * LambdaBl4 * XB*XB*XB*XB*XB*XB*XB*XB*XB * dyn->q;

        /* gravitomagnetic tidal coupling constants el = 2 only */
        dyn->kapA2j = 24.   * SigmaAl2 * XA*XA*XA*XA*XA / dyn->q;
        dyn->kapB2j = 24.   * SigmaBl2 * XB*XB*XB*XB*XB * dyn->q;

        dyn->kapT2  = dyn->kapA2 + dyn->kapB2;
        dyn->kapT3  = dyn->kapA3 + dyn->kapB3;
        dyn->kapT4  = dyn->kapA4 + dyn->kapB4;
        dyn->kapT2j = dyn->kapA2j + dyn->kapB2j;

        XLAL_CHECK_VOID(dyn->kapT2 > 0., XLAL_EDOM, "kappaT2 must be >= 0");
        XLAL_CHECK_VOID(dyn->kapT3 > 0., XLAL_EDOM, "kappaT3 must be >= 0");
        XLAL_CHECK_VOID(dyn->kapT4 > 0., XLAL_EDOM, "kappaT4 must be >= 0");

        /* Tidal coefficients cons dynamics
         \f$\bar{\alpha}_n^{(\ell)}\f$, Eq.(37) of Damour&Nagar, PRD 81, 084016 (2010) */
        dyn->bar_alph2_1 = (5./2.*XA*dyn->kapA2 + 5./2.*XB*dyn->kapB2)/dyn->kapT2;
        dyn->bar_alph2_2 = ((3.+XA/8.+ 337./28.*XA*XA)*dyn->kapA2 + (3.+XB/8.+ 337./28.*XB*XB)*dyn->kapB2)/dyn->kapT2;
        dyn->bar_alph3_1 = ((-2.+15./2.*XA)*dyn->kapA3 + (-2.+15./2.*XB)*dyn->kapB3)/dyn->kapT3;
        dyn->bar_alph3_2 = ((8./3.-311./24.*XA+110./3.*XA*XA)*dyn->kapA3 + (8./3.-311./24.*XB+110./3.*XB*XB)*dyn->kapB3)/dyn->kapT3;
        /* Gravitomagnetic term, see Eq.(6.27) if Bini-Damour-Faye 2012 */
        dyn->bar_alph2j_1 = ( dyn->kapA2j*(1. + (11./6.)*XA + XA*XA) + dyn->kapB2j*(1. + (11./6.)*XB + XB*XB) )/dyn->kapT2j;

        /* Tidal coefficients for the amplitude */
        dyn->khatA2  = 3./2. * LambdaAl2 * XB/XA * (REAL8) gsl_pow_int((double) XA,5);
        dyn->khatB2  = 3./2. * LambdaBl2 * XA/XB * (REAL8) gsl_pow_int((double) XB,5);

        /* Self-spin coefficients */
        dyn->C_Q1   = 1.;
        dyn->C_Q2   = 1.;
        dyn->C_Oct1 = 1.;
        dyn->C_Oct2 = 1.;
        dyn->C_Hex1 = 1.;
        dyn->C_Hex2 = 1.;
        if (LambdaAl2>0.)
        {
            REAL8 logC_Q1 = logQ(log(LambdaAl2));
            dyn->C_Q1           = exp(logC_Q1);
            dyn->C_Oct1         = Yagi14_fit_Coct(dyn->C_Q1);
            dyn->C_Hex1         = Yagi14_fit_Chex(dyn->C_Q1);
        }
        if (LambdaBl2>0.)
        {
            REAL8 logC_Q2 = logQ(log(LambdaBl2));
            dyn->C_Q2           = exp(logC_Q2);
            dyn->C_Oct2         = Yagi14_fit_Coct(dyn->C_Q2);
            dyn->C_Hex2         = Yagi14_fit_Chex(dyn->C_Q2);
        }

        /* Set more as needed ... */
        ssorder = SS_NLO;
        rc_order = RC_NNLO;
        c3 = 0.;

        /* Use default modes for BNS */
        teobModeChoice = TEOB_MODES_BNS_DEFAULT;
    }
    else
    {
        ssorder = SS_LO;
        rc_order = RC_NLO;
        c3 = eob_c3_fit_global(dyn->nu, dyn->chi1, dyn->chi2, dyn->X1, dyn->X2, dyn->a1, dyn->a2);
        /* Use default modes for BBH */
        teobModeChoice = TEOB_MODES_BBH_DEFAULT;
    }
    dyn->cN3LO = c3;

    if (*pModeArray == NULL)
    {
        *pModeArray = XLALSimInspiralCreateModeArray();
        XLALSetup_TEOB_mode_array(*pModeArray, teobModeChoice);
    }
    XLALCheck_TEOB_mode_array(*pModeArray, dyn->use_tidal);

    /* Default NQC settings:
     *    - if tides are on no NQC for flux and hlm (BNS/BHNS)
     *    - if tides are off and spins are on, calculate NQC for hlm but not for flux (BBH spins)
     *    - if tides and spins are off, then calculate NQC for flux and hlm using nrfit_nospin201602 (BBH no spins)
     * TODO: parse optional flags in LALDict
     */

    if (dyn->use_tidal)
    {
        /* BNS setup */
        dyn->nqc_coeffs_hlm = NQC_OFF;
        dyn->nqc_coeffs_flx = NQC_OFF;
    }
    else
    {
        if (dyn->use_spins)
        {
            /* BBH spins setup */
            dyn->nqc_coeffs_hlm = NQC_COMPUTE;
            dyn->nqc_coeffs_flx = NQC_OFF;
        }
        else
        {
            /* BBH no spins setup */
            dyn->nqc_coeffs_hlm = NQC_NR_NOSPIN;
            dyn->nqc_coeffs_flx = NQC_NR_NOSPIN;
        }
    }

    /* NQC data */
    eob_nqc_setcoefs(dyn);

    if (dyn->use_tidal)
    {
        /* back up */
        INT4 utidal_tmp, uspin_tmp;
        utidal_tmp = dyn->use_tidal;
        uspin_tmp  = dyn->use_spins;

        /* Temporarily set tidal to NNLO and switch off spin for LR estimation */
        dyn->use_tidal = LAL_SIM_INSPIRAL_GETIDES_NNLO;
        dyn->use_spins = 0;
        int status;

        // TODO: tidesFlag usage?
        XLAL_CHECK_VOID((status = eob_dyn_adiabLR(dyn, &(dyn->rLR_tidal), utidal_tmp))==XLAL_SUCCESS, status, "Light ring not found! Exiting...");

        /* restore */
        dyn->use_tidal = utidal_tmp;
        dyn->use_spins = uspin_tmp;
    }

    /* p-order 2GSF pole in the tidal interaction.

        See Eq.(29) of TEOBResumS paper
        A. Nagar et al. Phys. Rev. D 98, 104052 (2018), arXiv:1806.01772 [gr-qc]
        originating from calculations in
        D. Bini and T. Damour, Phys.Rev. D90, 124037 (2014), arXiv:1409.6933 [gr-qc]
        and
        S. R. Dolan, P. Nolan, A. C. Ottewill, N. Warburton, and B. Wardell,
        Phys. Rev. D91, 023009 (2015), arXiv:1406.4890 [gr-qc]

        The interested user may set a custom value
        AT THEIR OWN RISK by passing the appropriate
        REAL8 parameter "pGSF_tidal" through LALDict.
    */
    dyn->pGSF_tidal = 4.0;
    if (XLALDictContains(LALparams, "pGSF_tidal"))
        dyn->pGSF_tidal = XLALDictLookupREAL8Value(LALparams, "pGSF_tidal");

    if (!(dyn->use_tidal))
    {
        HealyBBHFitRemnant(S1z, S2z, dyn->q, &(dyn->Mbhf), &(dyn->abhf));
        dyn->abhf = JimenezFortezaRemnantSpin(dyn->nu, dyn->X1, dyn->X2, S1z, S2z);
    }

    dyn->ode_timestep = ODE_TSTEP_ADAPTIVE;
    dyn->t_stop = 1e10;    // HARDCODED
    dyn->r_stop = 1.0;     // HARDCODED

    /* Function pointers */
    EOBWavFlmSFunc eob_wav_flm_s;
    EOBDynSGetRCFunc eob_dyn_s_get_rc;

    /* Set f_lm fun pointer */
    if (ssorder==SS_LO) {
        /* eob_wav_flm_s = &eob_wav_flm_s_old; */
        eob_wav_flm_s = &eob_wav_flm_s_SSLO;
    } else if (ssorder==SS_NLO) {
        eob_wav_flm_s = &eob_wav_flm_s_SSNLO;
    } else XLAL_ERROR_VOID(XLAL_EINVAL, "SS order not recognized.\n");
    dyn->eob_wav_flm_s = eob_wav_flm_s;

    /* Set rc fun pointer */
    if (rc_order==RC_LO)
    {
        eob_dyn_s_get_rc = &eob_dyn_s_get_rc_LO;
    }
    else if (rc_order==RC_NLO)
    {
        eob_dyn_s_get_rc = &eob_dyn_s_get_rc_NLO;
    }
    else if (rc_order==RC_NNLO)
    {
        eob_dyn_s_get_rc = &eob_dyn_s_get_rc_NNLO;
    }
    else if (rc_order==RC_NNLOS4)
    {
        eob_dyn_s_get_rc = &eob_dyn_s_get_rc_NNLO_S4;
    }
    else if (rc_order==RC_NOSPIN)
    {
        eob_dyn_s_get_rc = &eob_dyn_s_get_rc_NOSPIN;
    }
    else if (rc_order==RC_NOTIDES)
    {
        eob_dyn_s_get_rc = &eob_dyn_s_get_rc_NOTIDES;
    }
    else XLAL_ERROR_VOID(XLAL_EINVAL, "Centrifugal radius option not recognized.\n");
    dyn->eob_dyn_s_get_rc = eob_dyn_s_get_rc;

    /* Pre-compute coefficients for resummed amplitudes */
    eob_wav_flm_coeffs(dyn->nu, dyn->clm);

    return;
}

/* Convert time in sec to dimensionless and mass-rescaled units */
REAL8 time_units_factor(REAL8 M)
{
    return 1./(M*LAL_MTSUN_SI);
}

REAL8 time_units_conversion(REAL8 M, REAL8 t)
{
    return t/(M*LAL_MTSUN_SI);
}

/* Convert frequency in Hz to dimensionless radius */
REAL8 radius0(REAL8 M, REAL8 fHz)
{
    REAL8 x = (M*LAL_MTSUN_SI*fHz*2.*LAL_PI)/2.;
    return cbrt(1/(x*x));
}

/* Mass and angular momentum of the final black hole
 Healey, Lousto and Zochlower (HLZ),
 arXiv: 1406.7295, published as PRD 90, 104004 (2014)
 WARNING: the formula uses the convention that M2 > M1, so that
 chi2 should refer to the black hole with the largest
 mass. In the EOB code, this is given by chi1, since
 in EOB code we use the convention that M1 > M2

 Here it is q=M2/M1, with M2>M1

 Improved with (Eisco, Jisco) + iterative procedure 23/02/2016
 parameters (TABLE VI)
 */
void HealyBBHFitRemnant(REAL8 chi1, REAL8 chi2, REAL8 q, REAL8 *mass, REAL8 *spin)
{

    /* Final mass:                    Angular momentum: */

    REAL8 M0  =  0.951507;            REAL8 L0  =  0.686710;
    REAL8 K1  = -0.051379;            REAL8 L1  =  0.613247;
    REAL8 K2a = -0.004804;            REAL8 L2a = -0.145427;
    REAL8 K2b = -0.054522;            REAL8 L2b = -0.115689;
    REAL8 K2c = -0.000022;            REAL8 L2c = -0.005254;
    REAL8 K2d =  1.995246;            REAL8 L2d =  0.801838;
    REAL8 K3a =  0.007064;            REAL8 L3a = -0.073839;
    REAL8 K3b = -0.017599;            REAL8 L3b =  0.004759;
    REAL8 K3c = -0.119175;            REAL8 L3c = -0.078377;
    REAL8 K3d =  0.025000;            REAL8 L3d =  1.585809;
    REAL8 K4a = -0.068981;            REAL8 L4a = -0.003050;
    REAL8 K4b = -0.011383;            REAL8 L4b = -0.002968;
    REAL8 K4c = -0.002284;            REAL8 L4c =  0.004364;
    REAL8 K4d = -0.165658;            REAL8 L4d = -0.047204;
    REAL8 K4e =  0.019403;            REAL8 L4e = -0.053099;
    REAL8 K4f =  2.980990;            REAL8 L4f =  0.953458;
    REAL8 K4g =  0.020250;            REAL8 L4g = -0.067998;
    REAL8 K4h = -0.004091;            REAL8 L4h =  0.001629;
    REAL8 K4i =  0.078441;            REAL8 L4i = -0.066693;

    /* Parameters */
    REAL8 nu      = q/((1.+q)*(1.+q));

    /* Masses: convention here is that m2>m1 */
    REAL8 X2      = 0.5*(1.+sqrt(1.-4*nu));
    REAL8 X1      = 1.-X2;

    /* Spin variables */
    REAL8 s1      = X1*X1*chi1;
    REAL8 s2      = X2*X2*chi2;
    REAL8 S       = s1 + s2;
    REAL8 S2      = S*S;
    REAL8 S3      = S*S2;
    REAL8 S4      = S2*S2;
    REAL8 Delta   = X1/X2*s2 - X2/X1*s1 + s2 - s1;
    REAL8 Delta2  = Delta*Delta;
    REAL8 Delta3  = Delta*Delta2;
    REAL8 Delta4  = Delta2*Delta2;

    /* Mass ratio variables */
    REAL8 deltam  = -sqrt(1-4*nu); // X1 - X2
    REAL8 deltam2 =  deltam*deltam;
    REAL8 deltam3 =  deltam*deltam2;
    REAL8 deltam4 =  deltam*deltam3;
    REAL8 deltam6 =  deltam2*deltam4;

    /* Initialize the angular momentum */
    REAL8 a0 = s1 + s2;
    INT4 a0_sign = 0;

    if (a0==0)
    {
        a0_sign=0;
    } else if (a0>0)
    {
        a0_sign=1;
    } else
    {   // if (a0<0) {
        a0_sign=-1;
    }

    /* Set-up an interative procedure to compute properly the "isco" quantities */
    REAL8 a2;
    REAL8 Z1;
    REAL8 Z2;
    REAL8 risco;
    REAL8 uisco;
    REAL8 Eisco;
    REAL8 Jisco;
    REAL8 abh;
    REAL8 Mbh=0.;

    UINT4 i;
    for(i=0; i<20; i++)
    {
        a2     = a0*a0;
        Z1     = 1 + cbrt(1-a2)*(cbrt(1+a0) + cbrt(1-a0));
        Z2     = sqrt(3*a2 + Z1*Z1);
        risco  = 3 + Z2 - a0_sign*sqrt((3-Z1)*(3+Z1+2.*Z2));
        uisco  = 1./risco;
        Eisco  = (1 - 2.*uisco + a0*sqrt(uisco*uisco*uisco))/sqrt(1-3*uisco + 2*a0*sqrt(uisco*uisco*uisco));
        Jisco  = 2./(sqrt(3.*risco))*(3.*sqrt(risco)-2.*a0);

        /* Dimensionless spin: J/Mbh^2 */
        abh = (4*nu)*(4*nu)*(L0 + L1*S + L2a*Delta*deltam + L2b*S2 + L2c*Delta2 + L2d*deltam2 + L3a*Delta*S*deltam + L3b*S*Delta2 + L3c*S3 + L3d*S*deltam2 + L4a*Delta*S2*deltam + L4b*Delta3*deltam + L4c*Delta4 + L4d*S4 + L4e*Delta2*S2 + L4f*deltam4 + L4g*Delta*deltam3 + L4h*Delta2*deltam2 + L4i*S2*deltam2) + S*(1+8*nu)*deltam4 + nu*Jisco*deltam6;

        Mbh = (4*nu)*(4*nu)*(M0 + K1*S + K2a*Delta*deltam + K2b*S2 + K2c*Delta2 + K2d*deltam2 + K3a*Delta*S*deltam + K3b*S*Delta2 + K3c*S3 + K3d*S*deltam2 + K4a*Delta*S2*deltam + K4b*Delta3*deltam + K4c*Delta4 + K4d*S4 + K4e*Delta2*S2 + K4f*deltam4 + K4g*Delta*deltam3 + K4h*Delta2*deltam2 + K4i*S2*deltam2) + (1 + nu*(Eisco + 11))*deltam6;

        a0 = abh;
    }

    *mass = Mbh;
    *spin = abh;
}

/* Final spin fit of */
// TODO: Refs needed
REAL8 JimenezFortezaRemnantSpin(REAL8 nu, REAL8 X1, REAL8 X2, REAL8 chi1, REAL8 chi2)
{

    const REAL8 xnu     = sqrt(1.0-4.0*nu);
    const REAL8 Dchi    = chi1-chi2;
    const REAL8 S       = (X1*X1*chi1+X2*X2*chi2)/(X1*X1+X2*X2);
    const REAL8 a2      = 3.833;
    const REAL8 a3      = -9.49;
    const REAL8 a5      = 2.513;

    /* The functional form is taken from eq. (7), page 5. */
    REAL8 Lorb_spin_zero  = (1.3*a3*nu*nu*nu + 5.24*a2*nu*nu + 2.*sqrt(3)*nu)/(2.88*a5*nu + 1);

    /* Coeffcients taken from Table II, page 6: */
    REAL8 b1      = 1.00096;
    REAL8 b2      = 0.788;
    REAL8 b3      = 0.654;
    REAL8 b5      = 0.840;

    /* These values are taken from Table III, page 7: */
    REAL8 f21     = 8.774;
    REAL8 f31     = 22.83;
    REAL8 f50     = 1.8805;
    REAL8 f11     = 0.345225*f21 + 0.0321306*f31 - 3.66556*f50 + 7.5397;

    /* These values are taken from Table IV, page 10 */
    REAL8 f12     = 0.512;
    REAL8 f22     = -32.1;
    REAL8 f32     = -154;
    REAL8 f51     = -4.77;

    /* The following quantities were taken from the relation given in eq. (11), */
    /* page 7: fi3 = 64 - 64.*fi0 - 16.*fi1 - 4.*fi2; */
    REAL8 f13     = 64 - 16.*f11 - 4.*f12;
    REAL8 f23     = 64 - 16.*f21 - 4.*f22;
    REAL8 f33     = 64 - 16.*f31 - 4.*f32;
    REAL8 f53     = 64 - 64.*f50 - 16.*f51;

    /* this transformation is given in eq. (9), page (7) */
    REAL8 b1t     = b1*(f11*nu + f12*nu*nu + f13*nu*nu*nu);
    REAL8 b2t     = b2*(f21*nu + f22*nu*nu + f23*nu*nu*nu);
    REAL8 b3t     = b3*(f31*nu + f32*nu*nu + f33*nu*nu*nu);
    REAL8 b5t     = b5*(f50 + f51*nu + f53*nu*nu*nu);

    /* The functional form is taken from eq. (8), page 6. */
    REAL8 Lorb_eq_spin  = (0.00954*b3t*S*S*S + 0.0851*b2t*S*S - 0.194*b1t*S)/(1 - 0.579*b5t*S);

    /* These values are taken from Table IV, page 10: */
    REAL8 d10     = 0.322;
    REAL8 d11     = 9.33;
    REAL8 d20     = -0.0598;
    REAL8 d30     = 2.32;
    REAL8 d31     = -3.26;

    /* The functional form is taken from eq. (19a-c), page 10.*/
    REAL8 A1      = d10*xnu*nu*nu*(d11*nu+1);
    REAL8 A2      = d20*nu*nu*nu;
    REAL8 A3      = d30*xnu*nu*nu*nu*(d31*nu+1);

    /* The functional form is taken from eq. (15), page 9. */
    REAL8 Lorb_uneq_mass  = A1*Dchi + A2*Dchi*Dchi + A3*S*Dchi;

    return X1*X1*chi1+X2*X2*chi2 + Lorb_spin_zero + Lorb_eq_spin + Lorb_uneq_mass;
}

/* logQ-vs-log(lambda) fit of Table I of Yunes-Yagi
    here x = log(lambda) and the output is the log of the coefficient
    that describes the quadrupole deformation due to spin. */
REAL8 logQ(REAL8 x)
{
  const REAL8 ai = 0.194;
  const REAL8 bi = 0.0936;
  const REAL8 ci = 0.0474;
  const REAL8 di = -4.21e-3;
  const REAL8 ei = 1.23e-4;
  const REAL8 x2 = x*x;
  const REAL8 x3 = x*x2;
  const REAL8 x4 = x*x3;
  return ai + bi*x + ci*x2 + di*x3 + ei*x4;
}

/* Yagi 2013 fits for NS multipolar
 \f$\bar{\lambda}_\ell = 2 k_\ell/(C^{2\ell+1} (2\ell-1)!!)\f$
 Eq.(9,10),(61); Tab.I; Fig.8 http://arxiv.org/abs/1311.0872 */
REAL8 Yagi13_fit_barlamdel(REAL8 barlam2, int ell)
{
    if (barlam2<=0.) return 0.;
    REAL8 lnx = log(barlam2);
    REAL8 coeffs[5];
    if (ell == 3)
    {
        coeffs[0] =  2.52e-5;
        coeffs[1] = -1.31e-3;
        coeffs[2] =  2.51e-2;
        coeffs[3] =  1.18;
        coeffs[4] = -1.15;
    }
    else if (ell == 4)
    {
        coeffs[0] =  2.8e-5;
        coeffs[1] = -1.81e-3;
        coeffs[2] =  3.95e-2;
        coeffs[3] =  1.43;
        coeffs[4] = -2.45;
    }
    else
        XLAL_ERROR_REAL8(XLAL_EINVAL, "Yagi fits are for ell=3,4.");
    REAL8 lny = coeffs[0]*lnx*lnx*lnx*lnx+coeffs[1]*lnx*lnx*lnx+coeffs[2]*lnx*lnx+coeffs[3]*lnx+coeffs[4];
    return exp(lny);
}

/* Yagi 2013 fits for NS multipolar
 \f$\bar{\sigma_2}( \bar{\lambda}_2 )\f$
 Eq.(9,11),(61); Tab.I; Fig.9 http://arxiv.org/abs/1311.0872
 See also later erratum */
REAL8 Yagi13_fit_barsigmalambda(REAL8 barlam2)
{
    if (barlam2<=0.) return 0.;
    REAL8 lnx = log(barlam2);
    REAL8 coeffs[5];
    /*
     coeffs[4] = 0.126;
     coeffs[3] = 0.617;
     coeffs[2] = 2.81e-2;
     coeffs[1] = 3.59e-4;
     coeffs[0] = -3.61e-5;
     */
    coeffs[4] = -2.01;
    coeffs[3] =  0.462;
    coeffs[2] =  1.68e-2;
    coeffs[1] = -1.58e-4;
    coeffs[0] = -6.03e-6;
    REAL8 lny = coeffs[0]*lnx*lnx*lnx*lnx+coeffs[1]*lnx*lnx*lnx+coeffs[2]*lnx*lnx+coeffs[3]*lnx+coeffs[4];

    return -1.0*exp(lny);
}

/* Yagi et al. fits for C_Oct
 Eq. (90) and Table I of https://arxiv.org/abs/1403.6243 */
REAL8 Yagi14_fit_Coct(REAL8 C_Q)
{
    REAL8 A0  = -0.925;
    REAL8 B1  =  1.98;
    REAL8 nu1 =  0.273;

    REAL8 cubrootCoct = A0 + B1*pow(C_Q,nu1);

    return cubrootCoct*cubrootCoct*cubrootCoct;
}

/* Yagi et al. fits for C_Hex
 Eq. (90) and Table I of https://arxiv.org/abs/1403.6243 */
REAL8 Yagi14_fit_Chex(REAL8 C_Q)
{
    REAL8 A0  = -0.413;
    REAL8 B1  =  1.5;
    REAL8 nu1 =  0.466;

    REAL8 fourthrootChex = A0 + B1*pow(C_Q,nu1);

    return SQ(SQ(fourthrootChex));
}

// TODO: Refs needed
REAL8 JFAPG_fit_Sigma_Irrotational(REAL8 barlam2)
{
    if (barlam2<=0.) return 0.;
    REAL8 lnx = log(barlam2);
    REAL8 coeffs[6];

    coeffs[5] = -2.03;
    coeffs[4] =  0.487;
    coeffs[3] =  9.69e-3;
    coeffs[2] =  1.03e-3;
    coeffs[1] = -9.37e-5;
    coeffs[0] =  2.24e-6;
    REAL8 lny = coeffs[0]*lnx*lnx*lnx*lnx*lnx+coeffs[1]*lnx*lnx*lnx*lnx+coeffs[2]*lnx*lnx*lnx+coeffs[3]*lnx*lnx+coeffs[4]*lnx+coeffs[5];

    return -1.0*exp(lny);
}

// TODO: Refs needed
REAL8 JFAPG_fit_Sigma_Static(REAL8 barlam2)
{
    if (barlam2<=0.) return 0.;
    REAL8 lnx = log(barlam2);
    REAL8 coeffs[6];

    coeffs[5] = -2.66;
    coeffs[4] =  0.786;
    coeffs[3] =  -0.01;
    coeffs[2] =  1.28e-3;
    coeffs[1] = -6.37e-5;
    coeffs[0] =  1.18e-6;
    REAL8 lny = coeffs[0]*lnx*lnx*lnx*lnx*lnx+coeffs[1]*lnx*lnx*lnx*lnx+coeffs[2]*lnx*lnx*lnx+coeffs[3]*lnx*lnx+coeffs[4]*lnx+coeffs[5];

    return exp(lny);
}

/* Join two multipolar waveforms at t = to */
void XLALSphHarmPolarJoin (SphHarmPolarTimeSeries *hlma, SphHarmPolarTimeSeries *hlmb, REAL8 to)
{
    /* Time arrays are suppose to be ordered as
     hlma->tdata->data:  x x x x x x x x x
     hlmb->tdata->data:       o o o o o o o o o
     t0               :                |
     But they do not need to overlap or be uniformly spaced.
     Note to can be
     t0 > hlma->time[hlma->size-1] => extend the a waveform
     t0 < hlmb->time[0]            => join the whole b waveform
     Following checks enforce the above structure, if possible.
     */
    if (hlma->tdata->data[0] > hlmb->tdata->data[0])
    {
        SWAPTRS( hlma, hlmb );
    }
    if ((to > hlmb->tdata->data[hlmb->tdata->length-1]) || (to <= hlma->tdata->data[0])) return;

    /* Find indexes of closest elements to "to"
     If the two arrays exactly overlap at "to" this should give:
     time_a[ioa] = time_b[iob] */
    const int ioa = find_point_bisection(to, hlma->tdata->length, hlma->tdata->data, 1);
    int iob       = find_point_bisection(to, hlmb->tdata->length, hlmb->tdata->data, 1);

    if ( DEQUAL(hlma->tdata->data[ioa], hlmb->tdata->data[iob], 1e-10) ) iob++;

    /* Calculate the new size N */
    const int N  = ioa + hlmb->tdata->length - iob;

    /* Resize a */
    XLALResizeSphHarmPolarTimeSeries(hlma, 0, N);
    XLALResizeREAL8Sequence(hlma->tdata, 0, N);

    /* Copy the relevant part of b into a */
    for (int i = ioa; i < N; i++)
    {
        hlma->tdata->data[i] = hlmb->tdata->data[iob + i - ioa];
    }
    SphHarmPolarTimeSeries *this = hlma;
    SphHarmPolarTimeSeries *that = hlmb;
    while (this && that)
    {
        for (int i = ioa; i < N; i++)
        {
            XLAL_CHECK_VOID((this->l==that->l) && (this->m==that->m), XLAL_EFAULT, "Mismatch of l and m when joining modes.");
            this->ampl->data->data[i]  = that->ampl->data->data[iob  + i - ioa];
            this->phase->data->data[i] = that->phase->data->data[iob + i - ioa];
        }
        this = this->next;
        that = that->next;
    }
    XLAL_CHECK_VOID(!this && !that, XLAL_EFAULT, "SphHarmTimeSeries to be joined must have the same set of modes");
}
