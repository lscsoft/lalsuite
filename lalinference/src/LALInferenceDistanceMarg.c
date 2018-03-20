//
//  LALInferenceDistanceMarg.c
//  
//
//  Created by John Veitch on 20/03/2018.
//

#include "LALInferenceDistanceMarg.h"
#include <math.h>
double dist_integral(double rho_opt, double rho_match);

double dist_snr_pdf(double rho_opt, double rho_match, double dL)
{
    double A = (rho_opt / dL)*(rho_opt / dL);
    double B = rho_match/dL;
    return exp(-A + B);
}

