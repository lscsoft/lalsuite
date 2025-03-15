from . import pyseobnr_model

import lalsimulation as lalsim
from ..core.waveform import LALCompactBinaryCoalescenceGenerator


def gwsignal_get_waveform_generator(waveform_approximant):
    if waveform_approximant == "SEOBNRv5HM":
        wf_gen = pyseobnr_model.SEOBNRv5HM()
    elif waveform_approximant == "SEOBNRv5EHM":
        wf_gen = pyseobnr_model.SEOBNRv5EHM()
    elif waveform_approximant == "SEOBNRv5PHM":
        wf_gen = pyseobnr_model.SEOBNRv5PHM()
    else:
        try:
            lal_approx = lalsim.SimInspiralGetApproximantFromString(
                waveform_approximant
            )
        except:
            raise ValueError("Approximant not implemented in GWSignal!")

        if lalsim.SimInspiralImplementedFDApproximants(
            lal_approx
        ) or lalsim.SimInspiralImplementedTDApproximants(lal_approx):
            wf_gen = LALCompactBinaryCoalescenceGenerator(waveform_approximant)
        else:
            # Should never get here
            raise ValueError("Approximant not implemented in GWSignal!")
    return wf_gen
