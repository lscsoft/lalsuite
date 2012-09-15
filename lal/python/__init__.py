#
# Import the contents of the swig binding into our namespace.  The swig
# bindings are optional, so ignore the import error if they aren't
# installed.
#


try:
	from lal import *
except ImportError:
	pass
