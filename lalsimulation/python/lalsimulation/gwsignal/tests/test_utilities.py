from pycbc.filter import match, optimized_match
from pycbc.types import TimeSeries as PyCBCTimeSeries


def compute_match(h1, h2):
    h1_pycbc = PyCBCTimeSeries(h1.data, delta_t=h1.dt.value)
    h2_pycbc = PyCBCTimeSeries(h2.data, delta_t=h2.dt.value)

    m, _ = optimized_match(h1_pycbc, h2_pycbc)
    return m
