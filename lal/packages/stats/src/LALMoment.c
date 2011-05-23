/**
\author Tibbits, M M
\file

\brief Routine to compute various moments of data.

\heading{Description}
The data structure passed in is either a REAL8 or a REAL4 Sequence.  The only parameter is which
moment to calculate. The function the sums the data, calculates the average and then it returns the
average for the first moment, it returns the variance for the second moment, and it returns the n-th
moment about the mean for higher order moments.

\heading{Algorithm}
<ul>
<li> Find the mean (here referred to as \f$ \overline{x} \f$).</li>
<li> Sum, over all the elements, the quantity: \f$(x[k] - \overline{x})^{n}\f$</li>
<li> Divide the sum just made by N-1. Call it moment-n</li>
<li> If n is greater than 2:</li>
   <ul>
   <li> Sum, over all the elements, the quantity: \f$(x[k] - \overline{x})^{n}\f$</li>
   <li> <em>Divide the sum just made by N. Call it moment-n</em></li>
   </ul>
<li>Return moment-n</li>
</ul>

\heading{Uses}

Determination of a specific moment of a set of data.

\heading{Notes}

<ul>
<li> <em>Moments less than two are not allowed.</em></li>
<li> <em>The result structure must be Non-NULL when passed in.</em></li>
<li> <em>The function assumes that the length member of the data passed in is correct.</em></li>
</ul>



*/

#include <math.h>
#include <lal/LALMoment.h>

NRCSID( LALMOMENTC, "$Id$");

#define TYPECODE D
#define TYPE REAL8
#include "LALMoment_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE S
#define TYPE REAL4
#include "LALMoment_source.c"
#undef TYPECODE
#undef TYPE
