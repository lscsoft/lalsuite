
import os

from distutils.core import setup

version = "1.0"

setup(
  name = "ligo-gracedb",
  version = version,
  maintainer = "Brian Moe",
  maintainer_email = "brian.moe@ligo.org",
  description = "Gravity Wave Candidate Event Database"
  long_description = "",

  url = "http://www.lsc-group.phys.uwm.edu/daswg/gracedb.html",
  license = 'GPL',
  provides = ['ligo.gracedb'],
  packages = [ 'ligo.gracedb'],

  requires = ['ligo', 'M2Crypto'],

  scripts = [
    os.path.join('bin','gracedb'),
  ],

)
