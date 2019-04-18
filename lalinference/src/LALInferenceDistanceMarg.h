//
//  LALInferenceDistanceMarg.h
//  
//
//  Created by John Veitch on 20/03/2018.
//

#ifndef LALInferenceDistanceMarg_h
#define LALInferenceDistanceMarg_h

#include <stdio.h>

/**
 * Compute the integral
 Int p(d | D, rho_opt, rho_mf) dD
 * D is luminosity distance
 * rho_opt is optimal SNR at 1 Mpc
 * rho_mf is <d|h> at 1 Mpc
 */
double dist_integral(double rho_opt, double rho_match, double dist_min, double dist_max);

double dist_snr_pdf(double dL, void *args);

#endif /* LALInferenceDistanceMarg_h */
