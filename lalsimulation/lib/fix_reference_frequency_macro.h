/*This defines a macro to fix the reference frequency, depending on the waveform model */

/*
 * certain approximants adopt the convention that f_ref=0 refers to the start
 * of the waveform while other approximants adopt the convention that f_ref=0
 * refers to the end of the waveform. in the former case, this routine will
 * return the explicit value of f_ref, which will be f_min.
 */

#define FIX_REFERENCE_FREQUENCY(f_ref, f_min, approximant) \
    if (f_ref == 0) \
        switch (approximant) { \
        case SpinTaylorT1:  \
        case SpinTaylorT5:  \
        case SpinTaylorT3:  \
        case SpinTaylorT4:  \
        case SpinTaylorT5Fourier:  \
        case SpinTaylorT4Fourier:  \
        case SpinTaylorF2:  \
        case IMRPhenomP:  \
        case IMRPhenomPv2:  \
        case IMRPhenomPv3:  \
        case IMRPhenomPv3HM:  \
        case IMRPhenomPv2_NRTidal:  \
        case IMRPhenomPv2_NRTidalv2:  \
        case IMRPhenomXP:  \
        case IMRPhenomXP_NRTidalv2: \
        case IMRPhenomXPHM:  \
        case NRSur4d2s:  \
        case IMRPhenomT:  \
        case IMRPhenomTHM:  \
        case IMRPhenomTP:  \
        case IMRPhenomTPHM:  \
        case IMRPhenomXO4a: \
        case TEOBResumS:  \
            f_ref = f_min;   \
            break; \
        default:  \
            break;  \
        }
