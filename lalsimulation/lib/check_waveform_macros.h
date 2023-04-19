/* Internal utility macro to check all spin components are zero
   returns 1 if all spins zero, otherwise returns 0 */
#define checkSpinsZero(s1x, s1y, s1z, s2x, s2y, s2z) \
    (((s1x) != 0. || (s1y) != 0. || (s1z) != 0. || (s2x) != 0. || (s2y) != 0. || (s2z) != 0.) ? 0 : 1)

/* Internal utility macro to check that the second body's spin components are zero.
   Returns 1 if all components are zero, otherwise returns 0 */
#define checkCOSpinZero(s2x, s2y, s2z) \
    (((s2x) != 0. || (s2y) != 0. || (s2z) != 0.) ? 0 : 1)

/* Internal utility macro to check transverse spins are zero
   returns 1 if x and y components of spins are zero, otherwise returns 0 */
#define checkTransverseSpinsZero(s1x, s1y, s2x, s2y) \
    (((s1x) != 0. || (s1y) != 0. || (s2x) != 0. || (s2y) != 0. ) ? 0 : 1)

/* Internal utility macro to check aligned spins very close to equal
   returns 1 if z components of spins are very close to equal, otherwise returns 0 */
#define checkAlignedSpinsEqual(s1z, s2z) \
    ((fabs((s1z) - (s2z)) > 1e-6) ? 0 : 1)

/* Internal utility macro to check tidal parameters are zero
   returns 1 if both tidal parameters zero, otherwise returns 0 */
#define checkTidesZero(lambda1, lambda2) \
    (((lambda1) != 0. || (lambda2) != 0. ) ? 0 : 1)
    
#define SANITY_WARNINGS(m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, deltaT, fmin, params) \
    do { \
        if ((deltaT) > 1.) \
            XLALPrintWarning("XLAL Warning - %s: Large value of deltaT = %e requested.\nPerhaps sample rate and time step size were swapped?\n", __func__, (deltaT)); \
        if ((deltaT) < 1./16385.) \
            XLALPrintWarning("XLAL Warning - %s: Small value of deltaT = %e requested.\nCheck for errors, this could create very large time series.\n", __func__, (deltaT)); \
        if ((m1) < 0.09 * LAL_MSUN_SI) \
            XLALPrintWarning("XLAL Warning - %s: Small value of m1 = %e (kg) = %e (Msun) requested.\nPerhaps you have a unit conversion error?\n", __func__, (m1), (m1)/LAL_MSUN_SI); \
        if ((m2) < 0.09 * LAL_MSUN_SI) \
            XLALPrintWarning("XLAL Warning - %s: Small value of m2 = %e (kg) = %e (Msun) requested.\nPerhaps you have a unit conversion error?\n", __func__, (m2), (m2)/LAL_MSUN_SI); \
        if ((m1) + (m2) > 1000. * LAL_MSUN_SI) \
            XLALPrintWarning("XLAL Warning - %s: Large value of total mass m1+m2 = %e (kg) = %e (Msun) requested.\nSignal not likely to be in band of ground-based detectors.\n", __func__, (m1)+(m2), ((m1)+(m2))/LAL_MSUN_SI); \
        if ((S1x)*(S1x) + (S1y)*(S1y) + (S1z)*(S1z) > 1.000001) \
            XLALPrintWarning("XLAL Warning - %s: S1 = (%e,%e,%e) with norm > 1 requested.\nAre you sure you want to violate the Kerr bound?\n", __func__, (S1x), (S1y), (S1z)); \
        if ((S2x)*(S2x) + (S2y)*(S2y) + (S2z)*(S2z) > 1.000001) \
            XLALPrintWarning("XLAL Warning - %s: S2 = (%e,%e,%e) with norm > 1 requested.\nAre you sure you want to violate the Kerr bound?\n", __func__, (S2x), (S2y), (S2z)); \
        if ((f_min) < 1.) \
            XLALPrintWarning("XLAL Warning - %s: Small value of fmin = %e requested.\nCheck for errors, this could create a very long waveform.\n", __func__, (f_min)); \
        if ((f_min) > 40.000001) \
            XLALPrintWarning("XLAL Warning - %s: Large value of fmin = %e requested.\nCheck for errors, the signal will start in band.\n", __func__, (f_min)); \
    } while (0)
