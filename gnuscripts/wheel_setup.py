"""
LIGO Scientific Collaboration Algorithm Library - minimal Python package

LALSuite is the LIGO Scientific Collaboration Algorithm Library for
gravitational-wave analysis. Its primary purpose is searching for and
characterizing astrophysical signals in gravitational-wave time series data,
particularly data from ground-based detectors such as `LIGO
<https://www.ligo.org>`_ and `Virgo <http://www.virgo-gw.eu>`_.

LALSuite consists of a set of ``configure``-``make``-``install`` style software
packages organized by problem domain or source classification. This Python
package provides a standalone, dependency-free binary distribution of the
libraries and Python modules in LALSuite for Linux and macOS.

To install, simply run::

    $ pip install lalsuite
"""

from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext as _build_ext

# Recursive globbing like in Python 3.5
from glob2 import glob


class build_ext(_build_ext):
    """Override setuptools build_ext: just copy pre-built modules."""

    def build_extension(self, ext):
        src = ext.name.replace('.', '/') + '.so'
        dst = self.get_ext_fullpath(ext.name)
        self.copy_file(src, dst)


setup(
    name='lalsuite',
    license='GPL-2+',
    description=__doc__.strip().split('\n')[0],
    long_description='\n'.join(__doc__.strip().split('\n')[1:]).strip(),
    url='https://git.ligo.org/lscsoft/lalsuite',
    author='LIGO Scientific Collaboration',
    author_email='lal-discuss@ligo.org',
    maintainer='Adam Mercer',
    maintainer_email='adam.mercer@ligo.org',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics'],
    cmdclass={'build_ext': build_ext},
    packages=find_packages(),
    package_data={'': ['*.so']},
    ext_modules=[Extension(filename.replace('.so', '').replace('/', '.'), [])
                 for filename in glob('lal*/**/*.so')],
    install_requires=['numpy>=1.7']
)
