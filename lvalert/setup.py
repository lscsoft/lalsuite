# 
# setup script for lvalert

import os

from distutils.core import setup

version = "1.1"

setup(
  name = "ligo-lvalert",
  version = version,
  maintainer = "Branson Stephens",
  maintainer_email = "branson.stephens@ligo.org",
  description = "LIGO-Virgo Alert Network",
  long_description = "The LIGO-Virgo Alert Network (LVAlert) is a prototype notification service built XMPP to provide a basic notification tool which allows multiple producers and consumers of notifications.",

  url = "http://www.lsc-group.phys.uwm.edu/daswg/lvalert.html",
  license = 'GPL',
  provides = ['ligo.lvalert'],
  packages = [ 'ligo.lvalert'],

  requires = ['ligo', 'pyxmpp', 'libxml2', 'M2Crypto', 'dnspython'],

  scripts = [
    os.path.join('bin','lvalert_admin'),
    os.path.join('bin','lvalert_send'),
    os.path.join('bin','lvalert_listen'),
  ],

)

# Dependencies.

# Python name     RHEL Name       Debian name
# (easy_install)  (yum)           (apt-get)
# 
# pyxmpp          pyxmpp *        python-pyxmpp
# libxml2 **      libxml2-python  python-libxml2
# M2Crypto        m2crypto        python-m2crypto
# dnspython       python-dns      python-dnspython
#
#   * lscsoft package
# 
#  ** Not in normal easy_install/pip place.  Also requires swig if installing with pip.
#        sudo apt-get install swig
#        pip install ftp://xmlsoft.org/libxml2/python/libxml2-python-2.6.9.tar.gz

