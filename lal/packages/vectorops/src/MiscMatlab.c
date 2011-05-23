/**
\author Tibbits, M. M.
\file

\heading{Module \ref MiscMatlab.c}
\latexonly\label{s_MiscMatlab_c}\endlatexonly

This file reproduces the last few matlab functions that we needed for our purposes.
It creates useable forms of cumsum, sum, max, and finally an implemenation of the
array addressing in matlab.  Matlab has an easy of inverting a vector, (end: -1: 1)
and the final function, FlipVector returns a result vector that has been flipped in
that same manner.

\heading{Description}

This file reproduces the last few matlab functions that we needed for our purposes.
It creates useable forms of cumsum, sum, max, and finally an implemenation of the
array addressing in matlab.  Matlab has an easy of inverting a vector, (end: -1: 1)
and the final function, FlipVector returns a result vector that has been flipped in
that same manner.

\heading{Algorithms}

The algorithms are the same as in matlab.  Flip vector was discussed above.  Sum
takes the sum of all of the elements in a vector.  Cum sum takes an input vector:

vector input[25];
vector output[25];

output[0] = input[0];
output[1] = input[0] + input[1];
output[2] = input[0] + input[1] + input[2];

etc

\heading{Uses}

<ul>
<li> \c LALDCreateVector</li>
</ul>

\heading{Notes}

At the current time none of the operations have been specified for neither the
complex datatypes nor the unsigned datatypes.

Also, the prototypes are out of order as I have used m4 to create all of the
functions from one codebase.

\heading{Prototypes}

*/

#include <math.h>
#include <lal/Matrix.h>

NRCSID( MATLABMATRIXSUMC, "$Id$");

#define TYPECODE D
#define TYPE REAL8
#include "MiscMatlab_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE S
#define TYPE REAL4
#include "MiscMatlab_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I2
#define TYPE INT2
#include "MiscMatlab_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I4
#define TYPE INT4
#include "MiscMatlab_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I8
#define TYPE INT8
#include "MiscMatlab_source.c"
#undef TYPECODE
#undef TYPE
