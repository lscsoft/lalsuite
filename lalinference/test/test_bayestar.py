import sys
import doctest
import lalinference.bayestar.sky_map
import lalinference.bayestar.filter
import lalinference.fits

print('Running C unit tests.')
total_failures = lalinference.bayestar.sky_map.test()

print('Running Python unit tests.')
modules = [
    lalinference.bayestar.filter,
    lalinference.fits]
for module in modules:
    failures, tests = doctest.testmod(module, verbose=True)
    total_failures += failures

if total_failures > 0:
    sys.exit(1)
