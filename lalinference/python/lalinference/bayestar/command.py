#
# Copyright (C) 2013-2016  Leo Singer
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
"""
Functions that support the command line interface.
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"


import argparse
import contextlib
from distutils.dir_util import mkpath
from distutils.errors import DistutilsFileError
import errno
import glob
import inspect
import itertools
import os
import shutil
import sys
import tempfile
from matplotlib import cm
from .. import cmap



@contextlib.contextmanager
def TemporaryDirectory(suffix='', prefix='tmp', dir=None, delete=True):
    try:
        dir = tempfile.mkdtemp(suffix=suffix, prefix=prefix, dir=dir)
        yield dir
    finally:
        if delete:
            shutil.rmtree(dir)


def chainglob(patterns):
    """Generate a list of all files matching a list of globs."""
    return itertools.chain.from_iterable(glob.iglob(s) for s in patterns)


waveform_parser = argparse.ArgumentParser(add_help=False)
group = waveform_parser.add_argument_group(
    'waveform options', 'Options that affect template waveform generation')
# FIXME: The O1 uberbank high-mass template, SEOBNRv2_ROM_DoubleSpin, does
# not support frequencies less than 30 Hz.
group.add_argument('--f-low', type=float, metavar='Hz', default=30,
    help='Low frequency cutoff [default: %(default)s]')
group.add_argument('--f-high-truncate', type=float, default=0.95,
    help='Truncate waveform at this fraction of the maximum frequency of the '
    'PSD [default: %(default)s]')
group.add_argument('--waveform', default='o1-uberbank',
    help='Template waveform approximant (e.g., TaylorF2threePointFivePN) '
    '[default: O1 uberbank mass-dependent waveform]')
del group


prior_parser = argparse.ArgumentParser(add_help=False)
group = prior_parser.add_argument_group(
    'prior options', 'Options that affect the BAYESTAR likelihood')
group.add_argument('--phase-convention', default='antifindchirp',
    choices=('findchirp', 'antifindchirp'),
    help='Phase convention [default: %(default)s]')
group.add_argument('--min-distance', type=float, metavar='Mpc',
    help='Minimum distance of prior in megaparsecs '
    '[default: infer from effective distance]')
group.add_argument('--max-distance', type=float, metavar='Mpc',
    help='Maximum distance of prior in megaparsecs '
    '[default: infer from effective distance]')
group.add_argument('--prior-distance-power', type=int, metavar='-1|2',
    default=2, help='Distance prior '
    '[-1 for uniform in log, 2 for uniform in volume, default: %(default)s]')
del group


skymap_parser = argparse.ArgumentParser(add_help=False)
group = skymap_parser.add_argument_group(
    'sky map output options', 'Options that affect sky map output')
group.add_argument('--nside', '-n', type=int, default=-1,
    help='HEALPix resolution [default: auto]')
group.add_argument('--chain-dump', default=False, action='store_true',
    help='For MCMC methods, dump the sample chain to disk [default: no]')
del group


class MatplotlibFigureType(argparse.FileType):
    def __init__(self):
        super(MatplotlibFigureType, self).__init__('wb')

    @staticmethod
    def __show():
        from matplotlib import pyplot as plt
        return plt.show()

    def __save(self):
        from matplotlib import pyplot as plt
        return plt.savefig(self.string)

    def __call__(self, string):
        if string == '-':
            return self.__show
        else:
            with super(MatplotlibFigureType, self).__call__(string):
                pass
            from matplotlib import pyplot as plt
            plt.switch_backend('agg')
            self.string = string
            return self.__save

class HelpChoicesAction(argparse.Action):
    def __init__(self,
                 option_strings,
                 choices=(),
                 dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS):
        name = option_strings[0].replace('--help-', '')
        super(HelpChoicesAction, self).__init__(
            option_strings=option_strings,
            dest=dest,
            default=default,
            nargs=0,
            help='show support values for --' + name + ' and exit')
        self._name = name
        self._choices = choices

    def __call__(self, parser, namespace, values, option_string=None):
        print('Supported values for --' + self._name + ':')
        for choice in self._choices:
            print(choice)
        parser.exit()

def type_with_sideeffect(type):
    def decorator(sideeffect):
        def func(value):
            ret = type(value)
            sideeffect(ret)
            return ret
        return func
    return decorator

@type_with_sideeffect(str)
def colormap(value):
    from matplotlib import rcParams
    rcParams['image.cmap'] = value

@type_with_sideeffect(float)
def figwith(value):
    from matplotlib import rcParams
    rcParams['figure.figsize'][0] = float(value)

@type_with_sideeffect(float)
def figheight(value):
    from matplotlib import rcParams
    rcParams['figure.figsize'][1] = float(value)

@type_with_sideeffect(int)
def dpi(value):
    from matplotlib import rcParams
    rcParams['figure.dpi'] = rcParams['savefig.dpi'] = float(value)

figure_parser = argparse.ArgumentParser(add_help=False)
colormap_choices = sorted(cm.cmap_d.keys())
group = figure_parser.add_argument_group(
    'figure options', 'Options that affect figure output format')
group.add_argument(
    '-o', '--output', metavar='FILE.{pdf,png}',
    default='-', type=MatplotlibFigureType(),
    help='name of output file [default: plot to screen]')
group.add_argument(
    '--colormap', default='cylon', choices=colormap_choices,
    type=colormap, metavar='CMAP',
    help='name of matplotlib colormap [default: %(default)s]')
group.add_argument(
    '--help-colormap', action=HelpChoicesAction, choices=colormap_choices)
group.add_argument(
    '--figure-width', metavar='INCHES', type=figwith, default='8',
    help='width of figure in inches [default: %(default)s]')
group.add_argument(
    '--figure-height', metavar='INCHES', type=figheight, default='6',
    help='height of figure in inches [default: %(default)s]')
group.add_argument(
    '--dpi', metavar='PIXELS', type=dpi, default=300,
    help='resolution of figure in dots per inch [default: %(default)s]')
del colormap_choices
del group


# Defer loading SWIG bindings until version string is needed.
class VersionAction(argparse._VersionAction):
    def __call__(self, parser, namespace, values, option_string=None):
        from .. import InferenceVCSVersion
        self.version = 'LALInference ' + InferenceVCSVersion
        super(VersionAction, self).__call__(
            parser, namespace, values, option_string)


class ArgumentParser(argparse.ArgumentParser):
    """
    An ArgumentParser subclass with some sensible defaults.

    - Any ``.py`` suffix is stripped from the program name, because the
      program is probably being invoked from the stub shell script.

    - The description is taken from the docstring of the file in which the
      ArgumentParser is created.

    - If the description is taken from the docstring, then whitespace in
      the description is preserved.

    - A ``--version`` option is added that prints the version of LALInference.
    """
    def __init__(self,
                 prog=None,
                 usage=None,
                 description=None,
                 epilog=None,
                 parents=[],
                 formatter_class=None,
                 prefix_chars='-',
                 fromfile_prefix_chars=None,
                 argument_default=None,
                 conflict_handler='error',
                 add_help=True):
        if prog is None:
            prog = os.path.basename(sys.argv[0]).replace('.py', '')
        if description is None:
            parent_frame = inspect.currentframe().f_back
            description = parent_frame.f_locals.get('__doc__', None)
            if formatter_class is None:
                formatter_class = argparse.RawDescriptionHelpFormatter
        if formatter_class is None:
            formatter_class = argparse.HelpFormatter
        super(ArgumentParser, self).__init__(
                 prog=prog,
                 usage=usage,
                 description=description,
                 epilog=epilog,
                 parents=parents,
                 formatter_class=argparse.RawDescriptionHelpFormatter,
                 prefix_chars=prefix_chars,
                 fromfile_prefix_chars=fromfile_prefix_chars,
                 argument_default=argument_default,
                 conflict_handler=conflict_handler,
                 add_help=add_help)
        self.add_argument('--version', action=VersionAction)


class DirType(object):
    """Factory for directory arguments."""

    def __init__(self, create=False):
        self._create = create

    def __call__(self, string):
        if self._create:
            try:
                mkpath(string)
            except DistutilsFileError as e:
                raise argparse.ArgumentTypeError(e.message)
        else:
            try:
                os.listdir(string)
            except OSError as e:
                raise argparse.ArgumentTypeError(e)
        return string


class SQLiteType(argparse.FileType):
    """Open an SQLite database, or fail if it does not exist.
    FIXME: use SQLite URI when we drop support for Python < 3.4.
    See: https://docs.python.org/3.4/whatsnew/3.4.html#sqlite3"""

    def __init__(self, mode='r'):
        super(SQLiteType, self).__init__(mode + 'b')

    def __call__(self, string):
        if string == '-':
            raise argparse.ArgumentTypeError(
                'Cannot open stdin/stdout as an SQLite database')
        with super(SQLiteType, self).__call__(string):
            import sqlite3
            return sqlite3.connect(string)


def sqlite_get_filename(connection):
    """Get the name of the file associated with an SQLite connection"""
    result = connection.execute('pragma database_list').fetchall()
    try:
        (_, _, filename), = result
    except ValueError:
        raise RuntimeError('Expected exactly one attached database')
    return filename


def rename(src, dst):
    """Like os.rename(src, dst), but works across different devices because it
    catches and handles EXDEV ('Invalid cross-device link') errors."""
    try:
        os.rename(src, dst)
    except OSError as e:
        if e.errno == errno.EXDEV:
            dir, suffix = os.path.split(dst)
            tmpfid, tmpdst = tempfile.mkstemp(dir=dir, suffix=suffix)
            try:
                os.close(tmpfid)
                shutil.copy2(src, tmpdst)
                os.rename(tmpdst, dst)
            except:
                os.remove(tmpdst)
                raise
        else:
            raise


def register_to_xmldoc(xmldoc, parser, opts, **kwargs):
    from glue.ligolw.utils import process
    return process.register_to_xmldoc(
        xmldoc, parser.prog,
        {key: (value.name if hasattr(value, 'read') else value)
        for key, value in opts.__dict__.items()})
