/*
 * Copyright (C) 2013 M. Carney, L. Wade
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */
/**
 * @addtogroup LALSimNeutronStarEOS_c
 * @{
 */
/**
 * @name Creation routines for spectral decomposition equations of state
 * @{
 */

#include <lal/LALSimReadData.h>
#include <gsl/gsl_interp.h>
#include <lal/LALSimNeutronStar.h>

/** @cond */

/* Maps the abscissae from [-1,1] to [a,b] */
static void GLBoundConversion(double a, double b, double abscissae[], int nEval)
{
    int i;

    // Converting to new evalutation points
    for(i=0;i<nEval;i++)
    {
        abscissae[i]=((b-a)/2.0)*abscissae[i] + (a+b)/2.0;
    }
}

/* Resets 10 point array to standard abscissae */
static void resetAbscissae(double abscissae[])
{
       abscissae[0] = -0.9739065285171717;
       abscissae[1] = -0.8650633666889845;
       abscissae[2] = -0.6794095682990244;
       abscissae[3] = -0.4333953941292472;
       abscissae[4] = -0.1488743389816312;
       abscissae[5] =  0.1488743389816312;
       abscissae[6] =  0.4333953941292472;
       abscissae[7] =  0.6794095682990244;
       abscissae[8] =  0.8650633666889845;
       abscissae[9] =  0.9739065285171717;
}

/* Specral decomposition of eos's adiabatic index */
static double AdiabaticIndex(double gamma[], double x, int size)
{
        double Gamma, logGamma = 0;
        int i;
        for(i=0;i<size;i++)
        {
          logGamma += gamma[i]*pow(x,(double) i);
        }
        Gamma = exp(logGamma);
        return Gamma;
}

/* Returns energy density given a pressure and spectral decomposition parameters */
static double eos_e_of_p_spectral_decomposition(double x, double gamma[], int size, double p0, double e0)
{
    // Integration/Placeholder variables
    int i;
    int j;
    int nEval = 10;
    double Integral;
    double IPrime;
    double e, mu;

    // Declaring arrays needed for Gauss-Legendre abscissae and weights
    double abscissae[nEval], abscissaePrime[nEval];
    resetAbscissae(abscissae);
    double weights[] = {
       0.0666713443086881,
       0.1494513491505806,
       0.2190863625159820,
       0.2692667193099963,
       0.2955242247147529,
       0.2955242247147529,
       0.2692667193099963,
       0.2190863625159820,
       0.1494513491505806,
       0.0666713443086881 
       };

    /* Calculation of \mu(x) */
    /* Calculation based on integral in Eq. 8 of PRD 82 103011 (2010): 

       \mu(x) = \exp[ - Integral(x) ] 

       Intergral(x) = \int_{0}^{x} 1 / \Gamma(xp) dxp

                    = (x/2) \int_{-1}^{1} 1 / \Gamma[ y(t) ] dt,

       where y(t) = (x/2) t + x/2.

       10-point Gaussian quadrature integration
    */
    Integral = 0.0;

    /* Shift abscissae, anticipating a change of integration limits */
    GLBoundConversion(0.0,x,abscissae, nEval);
    for(i=0;i<nEval;i++)
    {
      /* Integral for mu(x) calculation  */
      Integral += weights[i]*pow(AdiabaticIndex(gamma,abscissae[i],size),-1.0);
    }

    Integral*=(x/2.0);
    mu = exp(-Integral);
    /* end \mu(x) calculation */

    /* Calculation of \epsilon(x) */
    /* Calculation based on integral in Eq. 7 of PRD 82 103011 (2010):

       \epsilon(x) = \epsilon_0 / \mu(x) + p_0 * Integral / \mu(x)

       Intergral(x) = \int_{0}^{x} \mu(xp) \exp(xp) / \Gamma(xp) dxp

                    = (x/2) \int_{-1}^{1} \mu[ y(t) ] \exp[ y(t) ] / \Gamma[ y(t) ] dt,

       where y(t) = (x/2) t + x/2.

       10-point Gaussian quadrature integration
    */
    Integral = 0.0;

    for(i=0;i<nEval;i++)
    {
      /* calculate \mu(xp) */
      IPrime = 0.0;
      resetAbscissae(abscissaePrime);
      GLBoundConversion(0.0,abscissae[i],abscissaePrime, nEval);

      for(j=0;j<nEval;j++)
      {
        IPrime += weights[j]*pow(AdiabaticIndex(gamma,abscissaePrime[j],size),-1.0);
      }

      IPrime *= (abscissae[i]/2.0);
      IPrime = exp(-IPrime);
      /* end mu(xp) calculation, and calculate integral in epsilon(x) */
      Integral += weights[i]*(exp(abscissae[i])*IPrime/AdiabaticIndex(gamma,abscissae[i],size));
    }

    Integral*=(x/2.0);

    e = e0/mu +(p0/mu)*Integral;
    /* end epsilon(x) calculation */

    return e;
}

/** @endcond */

/**
 * @brief Reads spectral decomposition eos parameters to make an eos.
 * @details Reads an array of spectral decomposition eos parameters and the
 * array length to construct an eos using a spectral decomposition of the 
 * adiabatic index, outlined in PRD 82 103011 (2010).
 * @param[in] gamma[] Array of spectral decomposition eos parameters.
 * @param[in] size The length of the gamma array.
 * @return A pointer to neutron star equation of state structure.
 */
LALSimNeutronStarEOS *XLALSimNeutronStarEOSSpectralDecomposition(double gamma[], int size)
{
    LALSimNeutronStarEOS * eos;
    size_t ndat_low = 69;
    size_t ndat = ndat_low + 500;
    size_t i;

    // Minimum pressure and energy density of core EOS geom
    double e0 = 9.54629006e-11; // For reference: 1.1555e35 [cgs]
    double p0 = 4.43784199e-13; // For reference: 5.3716e32 [cgs]

    double xmax = 12.3081;
    double pmax = p0*exp(xmax); // For reference: 9.82905e-8 [geom], 1.18973e38 [cgs]

    // Declaring array and allocating memory space
    double *edat;
    double *pdat;
    double xdat[ndat-ndat_low];

    pdat = XLALCalloc(ndat, sizeof(*pdat));
    edat = XLALCalloc(ndat, sizeof(*edat));

    // Low density EOS values to be stitched (SLy, in geom)
    double pdat_low[]={
         0.00000000e+00,   2.49730009e-31,   1.59235347e-30,
         1.01533235e-29,   6.47406376e-29,   4.12805731e-28,
         2.63217321e-27,   1.67835262e-26,   1.07016799e-25,
         6.82371225e-25,   4.35100369e-24,   1.91523482e-23,
         8.06001537e-23,   3.23144235e-22,   4.34521997e-22,
         1.18566090e-21,   3.16699528e-21,   8.31201995e-21,
         2.15154075e-20,   5.51600847e-20,   7.21972469e-20,
         1.34595234e-19,   2.50269468e-19,   3.41156366e-19,
         4.16096744e-19,   5.66803746e-19,   1.05098304e-18,
         1.94663211e-18,   3.60407863e-18,   4.67819652e-18,
         6.36373536e-18,   8.65904266e-18,   1.17739845e-17,
         1.60126190e-17,   2.06809005e-17,   2.81253637e-17,
         3.82385968e-17,   4.91532870e-17,   6.68349199e-17,
         9.08868981e-17,   1.19805457e-16,   1.23523557e-16,
         1.67975513e-16,   2.14575704e-16,   2.38949918e-16,
         2.71834450e-16,   3.69579177e-16,   4.80543818e-16,
         6.22823125e-16,   6.44883854e-16,   6.51906933e-16,
         6.90079430e-16,   7.51717272e-16,   7.84188682e-16,
         8.12280996e-16,   8.94822824e-16,   1.78030908e-15,
         2.83170525e-15,   4.35257355e-15,   6.44272433e-15,
         9.21014776e-15,   1.72635532e-14,   2.96134301e-14,
         4.76007735e-14,   7.28061891e-14,   1.06973879e-13,
         1.78634067e-13,   3.17897582e-13,   4.16939514e-13};

    double edat_low[]={
         0.00000000e+00,   9.76800363e-24,   3.08891410e-23,
         9.76800489e-23,   3.08891490e-22,   9.76801000e-22,
         3.08891816e-21,   9.76803078e-21,   3.08893141e-20,
         9.76811525e-20,   3.08898527e-19,   7.75906672e-19,
         1.94909879e-18,   4.89622085e-18,   6.16357861e-18,
         1.22969164e-17,   2.45420997e-17,   4.89823567e-17,
         9.77310869e-17,   1.95024387e-16,   2.45531498e-16,
         3.89264229e-16,   6.16926264e-16,   7.76837038e-16,
         9.00642353e-16,   1.19237569e-15,   1.89037060e-15,
         3.09452823e-15,   4.90635934e-15,   5.96458915e-15,
         7.51061114e-15,   9.79776989e-15,   1.23353809e-14,
         1.55335650e-14,   1.88188196e-14,   2.46271194e-14,
         3.10096915e-14,   3.74392368e-14,   4.91814886e-14,
         6.19454280e-14,   7.62224545e-14,   8.11155415e-14,
         1.05200677e-13,   1.26436151e-13,   1.37076948e-13,
         1.55799897e-13,   2.02900352e-13,   2.47139208e-13,
         3.11281590e-13,   3.19529901e-13,   3.22140999e-13,
         3.36213572e-13,   3.58537297e-13,   3.70113877e-13,
         3.80033593e-13,   4.24821517e-13,   1.57119195e-12,
         2.36083738e-12,   3.33152493e-12,   4.49715211e-12,
         5.87555551e-12,   9.35905600e-12,   1.39898961e-11,
         2.00067862e-11,   2.76245081e-11,   3.69196076e-11,
         5.35224936e-11,   7.81020115e-11,   9.19476188e-11};

    // Populating first 69 points with low density EOS
    for(i=0;i<ndat_low;i++)
    {
      pdat[i]=pdat_low[i];
      edat[i]=edat_low[i];
    }

    // Generating higher density portion of EOS with spectral decomposition
    double logpmax = log(pmax);
    double logp0 = log(p0);
    double dlogp = (logpmax-logp0)/(ndat - ndat_low);

    // Calculating pressure table
    for(i=0;i < ndat-ndat_low;i++)
    {
      pdat[i+ndat_low] = exp(logp0 + dlogp*(i));
    }

    // Calculating xdat table
    for(i = 0;i<ndat - ndat_low;i++) {
      xdat[i] = log(pdat[i+ndat_low]/p0);
    }
    for(i = ndat_low;i < ndat;i++)
    {
      // Calculate energy density for each dimensionless pressure value
      edat[i] = eos_e_of_p_spectral_decomposition(xdat[i-ndat_low], gamma, size, p0,e0);
    }

    /* Updating eos structure */
    eos = eos_alloc_tabular(pdat,edat,ndat);

    if(snprintf(eos->name, sizeof(eos->name), "4-Param Spec Decomp (g0=%.4g, g1=%.4g, g2=%.4g, g3=%.4g)",
    gamma[0], gamma[1], gamma[2], gamma[3]) >= (int) sizeof(eos->name))
        XLAL_PRINT_WARNING("EOS name too long");

    return eos;
}

/**
 * @brief Reads 4 spectral decomposition eos parameters to make an eos.
 * @details Reads 4 spectral decomposition eos parameters and to construct an
 * eos using a spectral decomposition of the adiabatic index, outlined in 
 * PRD 82 103011 (2010).
 * @param[in] SDgamma0,SDgamma1,SDgamma2,SDgamma3 The spectral decomposition
 * eos parameters.
 * @return A pointer to neutron star equation of state structure.
 */
LALSimNeutronStarEOS *XLALSimNeutronStarEOS4ParameterSpectralDecomposition(double SDgamma0, double SDgamma1, double SDgamma2, double SDgamma3){
    double gamma[] = {SDgamma0, SDgamma1, SDgamma2, SDgamma3};
    LALSimNeutronStarEOS * eos;
    eos = XLALSimNeutronStarEOSSpectralDecomposition(gamma, 4);
    return eos;
}

/** @} */
/** @} */
