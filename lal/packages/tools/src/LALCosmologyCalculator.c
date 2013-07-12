/* Copyright (C) 2012 Walter Del Pozzo, Tjonnie Li
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
#include <lal/LALMalloc.h>
#include <lal/LALCosmologyCalculator.h>

/** 
* The set of functions in this module implements the standard cosmological distance measures 
* defined a Friedmann-Robertson-Walker-Lemaitre metric.
* For a detailed reference to the various formulae, see Hogg 1999 ( http://arxiv.org/abs/astro-ph/9905116 ). 
**/ 

/** 
* Computes the luminosity distance as the angular distance divided by (1+z)^2
* Eq. 21 in Hogg 1999 ( http://arxiv.org/abs/astro-ph/9905116 ) 
**/

double XLALLuminosityDistance(
            LALCosmologicalParameters *omega, 
            double z)
{
    double x=1.0+z;
    return XLALComovingTransverseDistance(omega, z)*x;
}

/** 
* Computes the angular size distance by dividing the comoving transverse distance by 1+z
* Eq. 18 in Hogg 1999 ( http://arxiv.org/abs/astro-ph/9905116 )
**/

double XLALAngularDistance(
            LALCosmologicalParameters *omega, 
            double z)
{
    double x=1.0+z;
    return XLALComovingTransverseDistance(omega, z)/x;
}

 /** 
* Computes the comoving line-of-sight distance (usually called the proper distance) by integrating the Hubble parameter.  
* Eq. 15 in Hogg 1999 ( http://arxiv.org/abs/astro-ph/9905116 ).
**/

double XLALComovingLOSDistance(
            LALCosmologicalParameters *omega, 
            double z)
{
    double dH = XLALHubbleDistance(omega);
    return dH*XLALIntegrateHubbleParameter(omega,z);
}
/**
* Computes the comoving transverse distance from the comoving line-of-sight distance (proper distance) and
* Eq. 16 in Hogg 1999 ( http://arxiv.org/abs/astro-ph/9905116 )
**/
double XLALComovingTransverseDistance(
            LALCosmologicalParameters *omega, 
            double z)
{
    double ok = omega->ok;
    double dH = XLALHubbleDistance(omega);
    double dC = XLALComovingLOSDistance(omega,z);
    
    if (fabs(ok)<1e-7) 
    {
        
        return dC;
    }
    else if (ok>1e-7)
    {
        return dH*sinh(sqrt(ok)*dC/dH)/sqrt(ok);
    }
    else if (ok<1e-7) 
    {
        return dH*sin(sqrt(fabs(ok))*dC/dH)/sqrt(fabs(ok));
    }
    else 
    {
        fprintf(stderr,"Something funny happened. Aborting.\n");
        exit(-1);
    }

}
/**
* Computes the Hubble distance, independent of the redshift. This is the distance with sets the units for all others.
* It returns the distance in Mpc. 
* Eq. 4 in Hogg 1999 ( http://arxiv.org/abs/astro-ph/9905116 )
**/
double XLALHubbleDistance(
            LALCosmologicalParameters *omega
            )
{
    double x = GSL_CONST_MKSA_SPEED_OF_LIGHT/(100.0*omega->h);
    return 1e-3*x;
}

/**
* Computes the inverse of the Hubble parameter at redshift z 
* Eq. 14 in Hogg 1999 ( http://arxiv.org/abs/astro-ph/9905116 )
**/

double XLALHubbleParameter(double z,
            void *omega
            )
{
    LALCosmologicalParameters *p = (LALCosmologicalParameters *) omega;

    double E=0.0;
    double x = 1.0+z;
    double om=p->om;
    double ol=p->ol;
    double ok=p->ok;
    double w0=p->w0;
    double w1=p->w1;
    double w2=p->w2;
    E = sqrt(om*x*x*x+ok*x*x+ol*pow(x,3.*(1.0+w0+w1+w2))*exp(-3.0*((w1+w2)*z/x + w2*z*z/(2.0*x*x))));
    return  1.0/E; 
}
/**
* Computes the integral of inverse of the Hubble parameter at redshift z.
* The integral is computed using the built-in function gsl_integration_qag of the gsl library. 
* The integral is performed using a Gauss-Kronrod adaptive method (see http://www.gnu.org/software/gsl/manual/html_node/QAG-adaptive-integration.html ) 
**/
double XLALIntegrateHubbleParameter(LALCosmologicalParameters *omega, double z)
{
    double result = 0.0;
    double error;
    double epsabs = 5e-5;
    double epsrel = 1e-5;
    size_t limit = 1024;
    int key = 1;
    
    gsl_function F;
    F.function = &XLALHubbleParameter;
    F.params  = omega;
    
    gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1024);

    gsl_integration_qag (&F, 0.0, z, epsabs, epsrel, 
                    limit, key, w, &result, &error);

    gsl_integration_workspace_free (w);

    return result;
}
/** 
* This function computes the value of a uniform probability distribution over the comoving volume.
* Details of the derivation of these formulae can be found in Coward, Burman 2005 ( http://arxiv.org/abs/astro-ph/0505181 )
**/
double XLALUniformComovingVolumeDistribution(
            LALCosmologicalParameters *omega, 
            double z,
            double zmax)
{
    double norm = XLALIntegrateComovingVolumeDensity(omega, zmax);          
    double unnorm_density = XLALUniformComovingVolumeDensity(z,(void *)omega);
    return unnorm_density/norm;
}

/** 
* This function computes the value of a uniform probability density distribution over the comoving volume at redshift z.
**/

double XLALUniformComovingVolumeDensity( 
            double z,
            void *omega)
{
    
    LALCosmologicalParameters *p = (LALCosmologicalParameters *)omega;

    double x = 1.0+z;
    double dm = XLALComovingTransverseDistance(omega,z);
    double E = XLALHubbleParameter(z,omega);
    double unnorm_density = 4.0*M_PI*dm*dm*E*XLALHubbleDistance(p)/x;
    return unnorm_density;
}
/**
* Function that integrates the uniform in comoving volume density to compute the normalisation factor. 
* Consistently with Coward, Burman 2005 ( http://arxiv.org/abs/astro-ph/0505181 ) the integral should be always performed up to infinity.
* However, if the user specifies a value for zmax, the probability density will be normalised by integrating up to that value. 
* To let the upper limit of integration go to infinity, specify zmax < 0. 
* The integration is performed using gsl_integration_qagiu from the gsl library (see http://www.gnu.org/software/gsl/manual/html_node/QAGI-adaptive-integration-on-infinite-intervals.html )
**/
double XLALIntegrateComovingVolumeDensity(LALCosmologicalParameters *omega, double z)
{
    double result = 0.0;
    double error;
    double epsabs = 5e-5;
    double epsrel = 1e-5;
    size_t limit = 1024;
    int key = 1;
    
    gsl_function F;
    F.function = &XLALUniformComovingVolumeDensity;
    F.params  = omega;
    
    gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1024);

    if (z<0.0) gsl_integration_qagiu (&F, 0.0, epsabs, epsrel, 
                    limit, w, &result, &error);
    
    else gsl_integration_qag (&F, 0.0, z, epsabs, epsrel, 
                    limit, key, w, &result, &error);

    gsl_integration_workspace_free (w);

    return result;
}
/**
* Creates a LALCosmologicalParameters structure from the values of the cosmological parameters.
* Note that the boundary condition \Omega_m + \Omega_k + \Omega_\Lambda = 1 is imposed here.  
**/
LALCosmologicalParameters *XLALCreateCosmologicalParameters(double h, double om, double ol, double w0, double w1, double w2)
{
    LALCosmologicalParameters *p = (LALCosmologicalParameters *)malloc(sizeof(LALCosmologicalParameters));
    p->h = h;
    p->om=om;
    p->ol=ol;
    p->ok=1.0-om-ol;
    p->w0=w0;
    p->w1=w1;
    p->w2=w2;
    return p;
}
/**
* Destroys a LALCosmologicalParameters structure.
**/

void XLALDestroyCosmologicalParameters(LALCosmologicalParameters *omega)
{
    LALFree(omega);
}
/**
* The next set of functions return specific parameters from the LALCosmologicalParameters structure 
**/

double XLALGetHubbleConstant(LALCosmologicalParameters *omega)
{
    return omega->h;
}
double XLALGetOmegaMatter(LALCosmologicalParameters *omega)
{
    return omega->om;
}
double XLALGetOmegaLambda(LALCosmologicalParameters *omega)
{
    return omega->ol;
}
double XLALGetOmegaK(LALCosmologicalParameters *omega)
{
    return omega->ok;
}
double XLALGetW0(LALCosmologicalParameters *omega)
{
    return omega->w0;
}
double XLALGetW1(LALCosmologicalParameters *omega)
{
    return omega->w1;
}
double XLALGetW2(LALCosmologicalParameters *omega)
{
    return omega->w2;
}
/**
* Function to set a LALCosmologicalParameters structure to the default LambdaCDM value (see http://arxiv.org/abs/1303.5076 )
**/
void XLALSetCosmologicalParametersDefaultValue(LALCosmologicalParameters *omega)
{
    omega->h = 0.671;
    omega->om=0.3175;
    omega->ok=1.0-0.3175-0.68251;
    omega->ol=0.68251;
    omega->w0=-1.0;
    omega->w1=0.0;
    omega->w2=0.0;
}
/**
* Function to create a LALCosmologicalRateParameters structure
**/
LALCosmologicalRateParameters *XLALCreateCosmologicalRateParameters(double r0, double W, double Q, double R)
{
    LALCosmologicalRateParameters *p = (LALCosmologicalRateParameters *)malloc(sizeof(LALCosmologicalRateParameters));
    p->r0=r0;
    p->Q=Q;
    p->W=W;
    p->R=R;
    return p;
}
/**
* Destroys a LALCosmologicalParameters structure.
**/

void XLALDestroyCosmologicalRateParameters(LALCosmologicalRateParameters *rate)
{
    LALFree(rate);
}
/**
* Returns the local rate 
**/
double XLALGetLocalRate(LALCosmologicalRateParameters *rate)
{
    return rate->r0;
}
/**
* Implements the fit to the SFR in Eq.7 of Coward, Burman 2005 ( http://arxiv.org/abs/astro-ph/0505181 )
* See also Porciani & Madau ( http://arxiv.org/pdf/astro-ph/0008294v2.pdf ) and references therein
**/
double XLALStarFormationDensity(double z, void *rate)
{
    LALCosmologicalRateParameters *p = (LALCosmologicalRateParameters *)rate;
    return p->r0*(1.0+p->W)*exp(p->Q*z)/(exp(p->R*z)+p->W);
}
/** 
* Returns the Rate weighted uniform comoving volume density
**/
double XLALRateWeightedUniformComovingVolumeDensity(double z, void *params)
{
    LALCosmologicalParametersAndRate *p = (LALCosmologicalParametersAndRate *)params;
    double dvdz = XLALUniformComovingVolumeDensity(z,p->omega);
    double ez = XLALStarFormationDensity(z,p->rate);
    return ez*dvdz;
}

/**
* Function that integrates the uniform in comoving volume density to compute the normalisation factor. 
* Consistently with Coward, Burman 2005 ( http://arxiv.org/abs/astro-ph/0505181 ) the integral should be always performed up to infinity.
* However, if the user specifies a value for zmax, the probability density will be normalised by integrating up to that value. 
* To let the upper limit of integration go to infinity, specify zmax < 0. 
* The integration is performed using gsl_integration_qagiu from the gsl library (see http://www.gnu.org/software/gsl/manual/html_node/QAGI-adaptive-integration-on-infinite-intervals.html )
**/

double XLALIntegrateRateWeightedComovingVolumeDensity(LALCosmologicalParametersAndRate *p, double z)
{
    double result = 0.0;
    double error;
    double epsabs = 5e-5;
    double epsrel = 1e-5;
    size_t limit = 1024;
    int key = 1;
    
    gsl_function F;
    F.function = &XLALRateWeightedUniformComovingVolumeDensity;
    F.params  = p;
    
    gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1024);

    if (z<0.0) gsl_integration_qagiu (&F, 0.0, epsabs, epsrel, 
                    limit, w, &result, &error);
    
    else gsl_integration_qag (&F, 0.0, z, epsabs, epsrel, 
                    limit, key, w, &result, &error);

    gsl_integration_workspace_free (w);

    return result;
}
/**
*  Returns the source redshifts distribution function obtained by normalizing the differential rate dR/dz integrated throughout the cosmos, as seen in our frame. 
*  It gives the probability density distribution function of an event occurring in the redshift range 0 to z, see Eq.12 in Coward & Burman ( http://arxiv.org/pdf/astro-ph/0505181v1.pdf )
**/
double XLALRateWeightedComovingVolumeDistribution(LALCosmologicalParametersAndRate *p, double z, double zmax)
{
    double norm = XLALIntegrateRateWeightedComovingVolumeDensity(p, zmax);          
    double unnorm_density = XLALRateWeightedUniformComovingVolumeDensity(z,(void *)p);
    return unnorm_density/norm;
}
/**
* Creates a cosmological parameters and rate structure
**/ 
LALCosmologicalParametersAndRate *XLALCreateCosmologicalParametersAndRate(void)
{
    LALCosmologicalParametersAndRate *p = (LALCosmologicalParametersAndRate *)malloc(sizeof(LALCosmologicalParametersAndRate));
    p->omega=XLALCreateCosmologicalParameters(0.0,0.0,0.0,0.0,0.0,0.0);
    p->rate=XLALCreateCosmologicalRateParameters(0.0,0.0,0.0,0.0);
    return p;
}
/**
* Destroys a cosmological parameters and rate structure
**/ 
void XLALDestroyCosmologicalParametersAndRate(LALCosmologicalParametersAndRate *p)
{
    XLALDestroyCosmologicalRateParameters(p->rate);
    XLALDestroyCosmologicalParameters(p->omega);
    LALFree(p);
}
/**
* Function to set a LALCosmologicalRateParameters structure to the default independent of z value 
**/
void XLALSetCosmologicalRateParametersDefaultValue(LALCosmologicalRateParameters *params)
{
    params->r0=5E-12;
    params->W= 0.0;
    params->Q= 0.0;
    params->R= 0.0;
}
