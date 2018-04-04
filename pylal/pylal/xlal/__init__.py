"""
This subpackage provides interfaces to LAL's functions and datatypes.
"""


from pylal import git_version


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date

__all__ = ["date"]
