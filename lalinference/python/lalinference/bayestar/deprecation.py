import inspect
import os
import re
import sys
import warnings
import textwrap

__all__ = ('warn',)


# Python <3.7 hides all DeprecationWarnings by default, whereas Python 3.7
# shows DeprecationWarnings that occur in main script (__main__) or in
# an interactive session. This subclass produces the Python 3.7 behavior.
class MovedWarning(DeprecationWarning):
    pass


warnings.filterwarnings(
    'default', category=MovedWarning, module='__main__')


def rewrap(s):
    return textwrap.fill(re.sub(r'\s+', ' ', s))


def warn(new=None, old=None, kind='module'):
    if old is None:
        old = inspect.getmodule(inspect.currentframe().f_back.f_code).__name__
        kind = 'module'
        if old == '__main__':
            old = os.path.basename(sys.argv[0])
            old.replace('.py', '')
            kind = 'tool'

    if new is None:
        new = old.replace('lalinference', 'ligo.skymap')

    kwargs = dict(old=old, new=new, kind=kind)

    message = '"{}" has moved to ligo.skymap.\n\n'.format(old)

    message += rewrap(
        '''The {old} {kind} from LALInference has been replaced by the {new}
        {kind} from the ligo.skymap package. The old {kind} is no longer
        tested and supported and will be removed in a future version of
        LALInference. Please update to the ligo.skymap package.'''.format(
            **kwargs))

    python_version = sys.version_info[:2]
    min_python_version = (3, 5)
    if python_version < min_python_version:
        message += '\n\n'
        message += rewrap(
            '''The ligo.skymap package requires Python {min_major}.{min_minor}
            or newer. You have Python {major}.{minor}. Before installing
            ligo.skymap, please update to Python {min_major}.{min_minor} or
            newer.'''.format(
                major=python_version[0],
                minor=python_version[1],
                min_major=min_python_version[0],
                min_minor=min_python_version[1]))

    message += '\n\n'
    message += rewrap(
        '''You can install ligo.skymap with pip by running this command in your
        shell:''')

    message += '\n\n'
    message += '    pip install ligo.skymap'
    message += '\n'

    warnings.warn(message, MovedWarning, stacklevel=3)
