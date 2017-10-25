import doctest, sys
from lalburst import date

sys.exit(bool(doctest.testmod(date).failed))
