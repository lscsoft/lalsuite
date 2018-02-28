"""
This module is deprecated.  Do not use.
"""

import warnings
warnings.warn("pylal.llwapp.get_coinc_def_id() is deprecated, use glue.ligolw.utils.coincs.get_coinc_def_id() instead", DeprecationWarning)
from glue.ligolw.utils.coincs import get_coinc_def_id
