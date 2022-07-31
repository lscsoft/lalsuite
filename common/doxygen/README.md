# LALSuite Doxygen documentation

Online HTML documentation for the LALSuite C libraries and Pure-python extras
is generated using [Doxygen](http://www.doxygen.nl/).

This directory (`/doxygen/`) contains the common files that are shared between
the doxygen builds for all of the LALSuite subpackages.

The build is controlled by the automake file `/gnuscripts/lalsuite_doxygen.am`,
customisations for each subpackage can be added in
`<package>/doxygen/Makefile.am`.
