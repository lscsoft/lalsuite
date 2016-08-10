import warnings
warnings.warn(
    'The module `lalinference.fits` has been deprecated. Please import '
    '`lalinference.io` or  `lalinference.io.fits` instead.')
from .io.fits import *
