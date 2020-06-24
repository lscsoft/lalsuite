import os, sys

# test that local TwoSpect executable can be called
os.system("lalapps_TwoSpect --version")

# return test status 'SKIP', since this test doesn't really do anything yet
sys.exit(77)
