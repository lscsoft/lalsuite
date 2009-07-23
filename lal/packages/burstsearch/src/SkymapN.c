/*
 * $Id$
 *
 * Copyright (C) 2008  Antony Searle
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

static void dummy(void);

#ifdef SKYMAP_N

#define CAT2(A, B) A##B
#define CAT(A, B) CAT2(A, B)


//static void CAT(diag, SKYMAP_N)(double a[SKYMAP_N][SKYMAP_N], double b[SKYMAP_N])
//{
//    int i;
//    int j;
//    for (i = 0; i != SKYMAP_N; ++i)
//        for (j = 0; j != SKYMAP_N; ++j)
//            a[i][j] = (i == j) ? b[i] : 0.;
//}

void CAT(XLALSkymap2DirectionPropertiesConstruct, SKYMAP_N)(
    CAT(XLALSkymap2PlanType, SKYMAP_N)* plan,
    XLALSkymap2SphericalPolarType* directions,
    CAT(XLALSkymap2DirectionPropertiesType, SKYMAP_N)* properties
    )
{
    double x[SKYMAP_N];
    int j;
    XLALSkymapCartesianFromSpherical(x, *directions);
    for (j = 0; j != SKYMAP_N; ++j)
    {
        properties->delay[j] = floor(site_time(plan->site + j, x) * plan->sampleFrequency + 0.5);
        site_response(properties->f[j], plan->site + j, x);
    }
}

void CAT(XLALSkymap2KernelConstruct, SKYMAP_N)(
    CAT(XLALSkymap2DirectionPropertiesType, SKYMAP_N)* properties,
    double wSw[SKYMAP_N],
    CAT(XLALSkymap2KernelType, SKYMAP_N)* kernel
    )
{
    
    int i, j, k, l;

    {

        //
        // W = diag(w.S_j^{-1}.w)
        //
        // F (F^T W F + I) F^T
        //
        
        double a[2][2]; // F^T W F
        double b[2][2]; // inv(A)

        // Compute the kernel

        // F^T W F + I
                
        for (i = 0; i != 2; ++i)
        {
            for (j = 0; j != 2; ++j)
            {
                a[i][j] = ((i == j) ? 1.0 : 0.0);
                for (k = 0; k != SKYMAP_N; ++k)
                {
                    a[i][j] += properties->f[k][i] * wSw[k] * properties->f[k][j];
                }
            }
        }
        
        // (F^T W F + I)^{-1}

        inv22(b, a);
        
        // F (F^T W F + I)^{-1} F^T
        
        for (i = 0; i != SKYMAP_N; ++i)
        {
            for (j = 0; j != SKYMAP_N; ++j)
            {
                kernel->k[i][j] = 0.0;
                for (k = 0; k != 2; ++k)
                {
                    for (l = 0; l != 2; ++l)
                    {
                         kernel->k[i][j] += properties->f[i][k] * b[k][l] * properties->f[j][l];
                    }
                }
            }
        }
        
        kernel->logNormalization = 0.5 * log(det22(b));
        
    }



    //{

        // Compute the normalization

    //    double a;
    //    double b[SKYMAP_N][SKYMAP_N];
    //    double c;

    //    a = wSw[0] * wSw[1] * wSw[2];

    //    for (i = 0; i != SKYMAP_N; ++i)
    //    {
    //        for (j = 0; j != SKYMAP_N; ++j)
    //        {
    //            b[i][j] = -kernel->k[i][j];
    //        }
    //        b[i][i] += 1. / wSw[i];
    //    }
    //    c = CAT(det, SKYMAP_N)(b);

    //    kernel->logNormalization = 0.5 * log(a * c);

    //}
}

static double CAT(CAT(ip, SKYMAP_N), SKYMAP_N)(double a[SKYMAP_N], double b[SKYMAP_N][SKYMAP_N], double c[SKYMAP_N])
{
    double d = 0.;
    int i;
    int j;
    for (i = 0; i != SKYMAP_N; ++i)
        for (j = 0; j != SKYMAP_N; ++j)
            d += a[i] * b[i][j] * c[j];
    return d;
}

void CAT(XLALSkymap2Apply, SKYMAP_N)(
    CAT(XLALSkymap2DirectionPropertiesType, SKYMAP_N)* properties,
    CAT(XLALSkymap2KernelType, SKYMAP_N)* kernel,
    double* xSw[SKYMAP_N],
    int tau,
    double* posterior
    )
{
    double x[SKYMAP_N];
    int j;
    for (j = 0; j != SKYMAP_N; ++j)
        x[j] = xSw[j][tau + properties->delay[j]];
    *posterior = 0.5 * CAT(CAT(ip, SKYMAP_N), SKYMAP_N)(x, kernel->k, x) + kernel->logNormalization;
}

void CAT(XLALSkymap2PlanConstruct, SKYMAP_N)(int sampleFrequency, int siteNumbers[SKYMAP_N], CAT(XLALSkymap2PlanType, SKYMAP_N)* plan)
{
    int i;
    static const char func[] = "XLALSkymap2ConstructPlanN";

    if (sampleFrequency <= 0)
    {
        XLALPrintError("%s(): invalid sample frequency %d\n", func, sampleFrequency);
        XLAL_ERROR_VOID(func, XLAL_EINVAL);
    }

    plan->sampleFrequency = sampleFrequency;

    for (i = 0; i != SKYMAP_N; ++i)
    {
        plan->site[i] = lalCachedDetectors[siteNumbers[i]];
    }
    
}

#undef CAT
#undef CAT2

#endif // SKYMAP_N

