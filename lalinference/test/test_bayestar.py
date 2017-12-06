import sys
import doctest
import lalinference.bayestar.sky_map
import lalinference.bayestar.filter

failures, tests = doctest.testmod(lalinference.bayestar.filter, verbose=True)

print('Running C unit tests.')
failures += lalinference.bayestar.sky_map.test()

if failures > 0:
    sys.exit(1)
