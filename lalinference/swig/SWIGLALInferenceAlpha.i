//
// Copyright (C) 2011--2014 Karl Wette
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA  02111-1307  USA
//

///
/// \defgroup SWIGLALInferenceAlpha_i Interface SWIGLALInferenceAlpha.i
/// \ingroup lalinference_swig
/// \brief SWIG code which must appear \e before the LALInference headers.
/// \author Karl Wette
///

// Local Variables:
// mode: c
// End:

%typemap(in) char *[] {
    // check if is a list
    if(PyList_Check($input))
    {
        int size = PyList_Size($input);
        int i = 0;
        $1 = (char **)malloc((size + 1)*sizeof(char *));
        for(i = 0; i < size; i++)
        {
            PyObject * o = PyList_GetItem($input, i);
            if(PyString_Check(o))
                $1[i] = PyString_AsString(o);
            else
            {
                PyErr_SetString(PyExc_TypeError, "list must contain strings");
                free($1);
                return NULL;
            }
        }
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "not a list");
        return NULL;
    }
}

// clean up the char ** array
%typemap(freearg) char *[] {
    free((char *) $1);
}

