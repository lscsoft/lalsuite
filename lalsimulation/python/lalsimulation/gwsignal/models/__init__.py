import sys
from importlib import import_module
import lalsimulation as lalsim

from ..core.waveform import LALCompactBinaryCoalescenceGenerator
from . import pyseobnr_model

GWSIGNAL_APPROXIMANTS = {
    "SEOBNRv5HM": ('pyseobnr_model', 'SEOBNRv5HM'),
    "SEOBNRv5PHM": ('pyseobnr_model', 'SEOBNRv5PHM'),
    "SEOBNRv5EHM": ('pyseobnr_model', 'SEOBNRv5EHM'),
    "TEOBResumSDALI": ('teobresums', 'TEOBResumSDALI'),
}

def gwsignal_get_waveform_generator(waveform_approximant):
    if waveform_approximant in GWSIGNAL_APPROXIMANTS:
        module_name, class_name = GWSIGNAL_APPROXIMANTS[waveform_approximant]
        module = import_module(f'.{module_name}', package=sys.modules[__name__].__package__)
        wf_gen = module.__getattribute__(class_name)()
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
