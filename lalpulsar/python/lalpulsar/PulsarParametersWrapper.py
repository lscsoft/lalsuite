# Copyright (C) 2018 Matthew Pitkin
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
"""

Wrapper class for the PulsarParameters structure, so that it can be accessed
in a dictionary-like way.

"""

from __future__ import division, absolute_import, print_function

import os
import re
import numpy as np

try:
    import lal
except ImportError:
    raise ImportError("SWIG wrappings of LAL cannot be imported")

try:
    import lalpulsar
except ImportError:
    raise ImportError("SWIG wrappings of LALPulsar cannot be imported")

try:
    from six import string_types
except ImportError:
    raise ImportError("Could not import six")

class PulsarParametersPy(object):
    """
    A class to wrap the SWIG-wrapped lalpulsar.PulsarParameters structure.

    This class lets you access the structure in a more Pythonic way, as well as providing
    a nice format for holding pulsar (`.par`) parameter files.

    The class can be used to set numerical values (double precision, unsigned integers), strings,
    or vectors of floating point values, e.g.:

        >>> from lalpulsar.PulsarParametersWrapper import PulsarParametersPy
        >>> pppy = PulsarParametersPy() # an empty structure
        >>> pppy['DECJ'] = 0.23         # set a numerical value
        >>> pppy['BINARY'] = 'BT'       # set a string value
        >>> pppy['F'] = [10.2, 1.4e-11] # set a vector of float values

    Args:
        pp (PulsarParameters, str): a lalpulsar.PulsarParameters structure, or a string giving the
            path to a TEMPO-style (`.par`) pulsar parameter file. If nothing is given then an empty
            lalpulsar.PulsarParameters structure is created. The `read()` method can subsequently
            be used to read in a `.par` file, or parameters can be added.

    Examples:
        An example of initialising the class with a previously created `lalpulsar.PulsarParameters`
        structure is:

            >>> import lalpulsar
            >>> from lalpulsar.PulsarParametersWrapper import PulsarParametersPy
            >>> # read in a pulsar parameter file
            >>> pp = lalpulsar.ReadTEMPOParFile('apulsar.par')
            >>> # view as a PulsarParametersPy object
            >>> pppy = PulsarParametersPy(pp)

        The same thing could be achieved more directly using:

            >>> pppy = PulsarParametersPy('apulsar.par')

        or, equivalently with:

            >>> pppy = PulsarParametersPy()
            >>> pppy.read('apulsar.par')

         parameter can be set in the class using, e.g.
    """

    keynames = []   # parameter names in PulsarParameters structure
    length = 0      # number of parameters
    _pulsarparameters = None

    def __init__(self, pp=None):
        # if pp is None create empty PulsarParameters structure
        if pp is None:
            self._pulsarparameters = lalpulsar.PulsarParameters()
        else:
            # check if pp is a pulsar parameters type or a (par file)
            if not isinstance(pp, lalpulsar.PulsarParameters) and isinstance(pp, string_types):
                if os.path.isfile(pp):
                    # try reading in file
                    self.read(pp)
                else:
                    raise ValueError("Input string does not point to a file")
            elif isinstance(pp, lalpulsar.PulsarParameters):
                self._pulsarparameters = pp
            else:
                raise ValueError("Expected 'lalpulsar.PulsarParameters' type, string, or None")

    def __len__(self):
        return length

    def __getitem__(self, key):
        """
        Get value from pulsar parameters
        """

        if self._pulsarparameters is None:
            return None

        # check if key finishes with "_ERR", in which case check for error value
        geterr = False
        tkey = key
        if key[-4:].upper() == "_ERR":
            geterr = True
            tkey = key[:-4] # get the actual parameter key name

        # check if the key is asking for an individal parameter from a vector parameter
        # (e.g. 'F0' gets the first value from the 'F' vector)
        sname = re.sub(r'_\d', '', tkey) if '_' in tkey else re.sub(r'\d', '', tkey)
        sidx = None
        indkey = None
        if sname != tkey:
            # check additional index is an integer
            try:
                sidx = int(tkey.split('_')[-1]) if '_' in tkey else int(tkey[len(sname):])
            except ValueError:
                pass

            # change tkey for checking parameter exists
            if sidx is not None:
                indkey = tkey # key with index
                tkey = sname

        # check if parameter given by the key is present
        if not lalpulsar.PulsarCheckParam(self._pulsarparameters, tkey):
            return None

        # get type of parameter
        ptype = lalpulsar.PulsarGetParamType(self._pulsarparameters, tkey)

        if ptype == lalpulsar.PULSARTYPE_REAL8_t:
            if not geterr:
                value = lalpulsar.PulsarGetREAL8Param(self._pulsarparameters, tkey)
            else:
                value = lalpulsar.PulsarGetREAL8ParamErr(self._pulsarparameters, tkey)
        elif ptype == lalpulsar.PULSARTYPE_REAL8Vector_t:
            if not geterr:
                if sidx is None:
                    tmpvalue = lalpulsar.PulsarGetREAL8VectorParam(self._pulsarparameters, tkey)
                    value = tmpvalue.data # 'data' in a REAL8Vector gets returned as a numpy array
                else:
                    value = lalpulsar.PulsarGetREAL8VectorParamIndividual(self._pulsarparameters, indkey)
            else:
                if sidx is None:
                    tmpvalue = lalpulsar.PulsarGetREAL8VectorParamErr(self._pulsarparameters, tkey)
                    value = tmpvalue.data
                else:
                    value = lalpulsar.PulsarGetREAL8VectorParamErrIndividual(self._pulsarparameters, indkey)
        elif ptype == lalpulsar.PULSARTYPE_string_t:
            if not geterr:
                value = lalpulsar.PulsarGetStringParam(self._pulsarparameters, tkey)
            else:
                raise ValueError("String-type cannot have an error")
        elif ptype == lalpulsar.PULSARTYPE_UINT4_t:
            if not geterr:
                value = lalpulsar.PulsarGetUINT4Param(self._pulsarparameters, tkey)
            else:
                raise ValueError("UINT4-type cannot have an error")
        else:
            raise ValueError("Unrecognised type")

        return value

    def __setitem__(self, key, value):
        """
        Set the value of a key
        """

        # if parameter exists remove it
        if lalpulsar.PulsarCheckParam(self._pulsarparameters, key):
            lalpulsar.PulsarRemoveParam(self._pulsarparameters, key)

        if isinstance(value, float):
            lalpulsar.PulsarAddREAL8Param(self._pulsarparameters, key, value)
        elif isinstance(value, string_types):
            lalpulsar.PulsarAddStringParam(self._pulsarparameters, key, value)
        elif isinstance(value, int):
            if value < 0.: # store negative integers as floats
                lalpulsar.PulsarAddREAL8Param(self._pulsarparameters, key, float(value))
            else:
                lalpulsar.PulsarAddUINT4Param(self._pulsarparameters, key, value)
        elif isinstance(value, list) or isinstance(value, np.ndarray):
            tarray = lal.CreateREAL8Vector(len(value))
            for i, tv in enumerate(value):
                if isinstance(tv, float):
                    tarray.data[i] = tv
                else:
                    raise ValueError("Non-float value in list or array")
            lalpulsar.PulsarAddREAL8VectorParam(self._pulsarparameters, key, tarray)
        else:
            raise ValueError("Data-type not one of know types")

    def keys(self):
        """
        Return a list of the parameter names stored in the PulsarParameters structure
        """

        thisitem = self._pulsarparameters.head
        self.keynames = [] # clear any previous key names
        self.length = 0
        while thisitem:
            tname = thisitem.name
            self.keynames.append(tname)
            self.length += 1

            thisitem = thisitem.next # move on to next value

        return self.keynames

    def values(self):
        """
        Return the values of each parameter in the structure
        """

        tvalues = []

        keys = self.keys()
        for key in keys:
            if lalpulsar.PulsarCheckParam(self._pulsarparameters, key):
                # get type of parameter
                ptype = lalpulsar.PulsarGetParamType(self._pulsarparameters, key)

                if ptype == lalpulsar.PULSARTYPE_REAL8_t:
                    value = lalpulsar.PulsarGetREAL8Param(self._pulsarparameters, key)
                elif ptype == lalpulsar.PULSARTYPE_REAL8Vector_t:
                    tmpvalue = lalpulsar.PulsarGetREAL8VectorParam(self._pulsarparameters, key)
                    value = tmpvalue.data # 'data' in a REAL8Vector gets returned as a numpy array
                elif ptype == lalpulsar.PULSARTYPE_string_t:
                    value = lalpulsar.PulsarGetStringParam(self._pulsarparameters, key)
                elif ptype == lalpulsar.PULSARTYPE_UINT4_t:
                    value = lalpulsar.PulsarGetUINT4Param(self._pulsarparameters, key)
                else:
                    raise ValueError("UINT4-type cannot have an error")
            else:
                raise ValueError("Could not find {} in strcuture".format(key))

            tvalues.append(value)

        return tvalues    

    def as_dict(self):
        """
        Return the contents (not error at the moment) of the structure as a dictionary
        """

        tdict = {}

        for tpair in self.items():
            tdict[tpair[0]] = tpair[1]

        return tdict

    def items(self):
        """
        Return list of item tuples for each parameter in the structure
        """

        tkeys = self.keys()
        tvalues = self.values()

        titems = []

        for (tk, tv) in zip(tkeys, tvalues):
            titems.append((tk, tv))

        return titems

    def read(self, filename):
        """
        Read a TEMPO-style parameter file into a PulsarParameters structure
        """

        # remove existing pulsarparameters
        if self._pulsarparameters is not None:
            del self._pulsarparameters

        pp = lalpulsar.ReadTEMPOParFile(filename)

        if pp is None:
            raise IOError("Problem reading in pulsar parameter file '{}'".format(filename))

        self._pulsarparameters = pp

    def PulsarParameters(self):
        """
        Return the PulsarParameters structure
        """

        return self._pulsarparameters
