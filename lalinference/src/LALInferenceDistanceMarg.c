//
//  LALInferenceDistanceMarg.c
//  
//
//  Created by John Veitch on 20/03/2018.
//

#include "LALInferenceDistanceMarg.h"
#include <math.h>
#include <gsl/gsl_integration.h>


struct integrand_args
{
    double A,B;
};

double dist_integral(double rho_opt, double rho_match, double dist_min, double dist_max)
{
    size_t limit=100;
    double result,abserr;
    gsl_function F;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc (limit);

    F.function = &dist_snr_pdf;
    struct integrand_args args = {rho_opt, rho_match};
    F.params = &args;
    
    gsl_integration_qag (&F, dist_min, dist_max, 0, 1e-6, limit, GSL_INTEG_GAUSS61, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    return(result / (dist_max-dist_min));
}


double dist_snr_pdf(double dL, void *args)
{
    struct integrand_args *a = (struct integrand_args *)args;
    return exp(-a->A/dL/dL + a->B/dL)*dL*dL;
}

