import os
import pwd


def existsInPath(name):
    for dir in os.getenv('PATH').split(os.pathsep):
        fullname = os.path.join(dir, name)
        if os.path.exists(fullname): return fullname


# this get_username() is more robust than os.getenv('USERNAME') or
# getpass.getuser(), because the latter rely on the user's environment
# and thus cannot be trusted, since the underlying variables can be
# (accidentally or intentionally) removed or changed.

def get_username():
    return pwd.getpwuid(os.getuid()).pw_name
