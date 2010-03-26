from distutils.core import setup
from distutils.extension import Extension
import os.path
import sys
import numpy

## {{{ http://code.activestate.com/recipes/502261/ (r1)
def pkgconfig(*packages, **kw):
    import commands
    flag_map = {'-I': 'include_dirs', '-L': 'library_dirs', '-l': 'libraries'}
    for token in commands.getoutput("pkg-config --libs --cflags %s" % ' '.join(packages)).split():
        kw.setdefault(flag_map.get(token[:2]), []).append(token[2:])
    return kw
## end of http://code.activestate.com/recipes/502261/ }}}

setup(
    name='pynds',
    author='Leo Singer',
    author_email='leo.singer@ligo.org',
    description='NDS1/NDS2 client',
    packages=['nds'],
    ext_modules=[
        Extension('nds_ext', ['nds_ext.cpp'],
            **pkgconfig('nds2-client',
                include_dirs=[numpy.get_include()],
                libraries=['boost_python']
            )
        )
    ]
)
