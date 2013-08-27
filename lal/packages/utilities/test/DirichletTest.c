/*
*  Copyright (C) 2007 John Whelan
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/**
   \file
   \ingroup Dirichlet_h
   \author UTB Relativity Group; contact whelan@phys.utb.edu

   \brief Test suite for LALDirichlet().

   \heading{Usage}
   \code
   ./DirichletTest
   \endcode

\heading{Description}

This program tests the function <tt>LALDirichlet()</tt>.
It tests all error conditions listed in the Error codes table.
It also writes to files the values of the Dirichlet kernel for three
different valid test cases.
See Figs.\figref{utilities_f_dirichlet_fig1}-\figref{utilities_f_dirichlet_fig3}.

\image html  utilitiesDirichletFig1.png "Fig. [utilities_f_dirichlet_fig1]: Dirichlet kernel for N=10, Delta x =0.01, and 0 \<= x \<= 1."
\image latex utilitiesDirichletFig1.pdf "Dirichlet kernel for N=10, Delta x =0.01, and 0 \<= x \<= 1." width=4in


\image html  utilitiesDirichletFig2.png "Fig. [utilities_f_dirichlet_fig2]: Dirichlet kernel for N=11, Delta x =0.01, and 0\<= x \<= 1."
\image latex utilitiesDirichletFig2.pdf "Dirichlet kernel for N=11, Delta x =0.01, and 0\<= x \<= 1." width=4in

\image html  utilitiesDirichletFig3.png "Fig. [utilities_f_dirichlet_fig3]: Dirichlet kernel for N=10, Delta x =0.01, and 0\<= x \<= 2."
\image latex utilitiesDirichletFig3.pdf "Dirichlet kernel for N=10, Delta x =0.01, and 0\<= x \<= 2." width=4in

\heading{Exit codes}
<table><tr><th>Code</th><th>Explanation</th></tr>
<tr><td>0</td><td>Success, normal exit.</td></tr>
<tr><td>1</td><td>Subroutine failed.</td></tr>
</table>

*/

/** \cond DONT_DOXYGEN */

#include <lal/LALStdlib.h>
#include <math.h>
#include <string.h>
#include <lal/AVFactories.h>
#include <lal/PrintVector.h>
#include <lal/Dirichlet.h>


int check ( LALStatus*, INT4, const CHAR* );

int
main( void )
{

  static LALStatus            status;
  REAL4Vector*                poutput = NULL;
  DirichletParameters         parameters;



#ifndef LAL_NDEBUG
  REAL4Vector                 dummy;
  if ( ! lalNoDebug )
  {
    parameters.length = 11 ;
    dummy.length = 11;
    /* test behavior for null pointer to input parameters */
    LALDirichlet( &status, &dummy, NULL );

    if ( check( &status, DIRICHLETH_ENULLPIN, DIRICHLETH_MSGENULLPIN ) )
    {
      return 1;
    }
    printf("PASS: %s\n", DIRICHLETH_MSGENULLPIN);

    /* test behavior for LALDirichlet parameter N <= 0  */

    parameters.n = 0;
    LALDirichlet( &status, &dummy, &parameters);
    if ( check( &status, DIRICHLETH_ENVALUE, DIRICHLETH_MSGENVALUE ) )
    {
      return 1;
    }
    printf("PASS: %s\n", DIRICHLETH_MSGENVALUE);
}
#endif /* LAL_NDEBUG */

  /* define valid value of N */
  parameters.n = 10 ;

#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
    parameters.length = 0;
    LALDirichlet( &status, &dummy, &parameters);
    if ( check( &status, DIRICHLETH_ESIZE, DIRICHLETH_MSGESIZE ) )
    {
      return 1;
    }
    printf("PASS: %s\n", DIRICHLETH_MSGESIZE);
  }
#endif /* LAL_NDEBUG */

  /* define valid value for specified length of output vector */
  parameters.length = 11 ;

#ifndef LAL_NDEBUG
  dummy.length = 11;
  if ( ! lalNoDebug )
  {
    /* test behavior for x spacing <= 0 */
    parameters.deltaX = -4.0;
    LALDirichlet( &status, &dummy, &parameters);
    if ( check( &status, DIRICHLETH_EDELTAX, DIRICHLETH_MSGEDELTAX ) )
    {
      return 1;
    }
    printf("PASS: %s\n", DIRICHLETH_MSGEDELTAX );
    parameters.deltaX = 0.0;
    LALDirichlet( &status, &dummy, &parameters);
    if ( check( &status, DIRICHLETH_EDELTAX, DIRICHLETH_MSGEDELTAX ) )
    {
      return 1;
    }
    printf("PASS: %s\n", DIRICHLETH_MSGEDELTAX );
  }
#endif  /* LAL_NDEBUG */

  /* define valid delta x */
  parameters.deltaX = 0.1;

#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
    /* test behavior for null pointer to output vector */
    LALDirichlet( &status, NULL, &parameters );
    if ( check( &status, DIRICHLETH_ENULLPOUT, DIRICHLETH_MSGENULLPOUT))
    {
      return 1;
    }
    printf("PASS: %s\n", DIRICHLETH_MSGENULLPOUT);

    /* test behavior for length of output vector not equal to length  */
    /* specified in input parameters */
    dummy.length = 10;
    LALDirichlet( &status, &dummy, &parameters );
    if ( check( &status, DIRICHLETH_ESIZEMM, DIRICHLETH_MSGESIZEMM ) )
    {
      return 1;
    }
    printf( "PASS: %s\n", DIRICHLETH_MSGESIZEMM );
  }
#endif  /* LAL_NDEBUG */

  /* assign valid output vector length */

#ifndef LAL_NDEBUG
  dummy.length = parameters.length;
  if ( ! lalNoDebug )
  {
    /* test behavior for null pointer to data member of output vector */
    dummy.data = NULL;
    LALDirichlet( &status, &dummy, &parameters );
    if ( check( &status, DIRICHLETH_ENULLPDOUT, DIRICHLETH_MSGENULLPDOUT))
    {
      return 1;
    }
    printf("PASS: %s\n", DIRICHLETH_MSGENULLPDOUT);
  }
# endif  /* LAL_NDEBUG */

  /* VALID TEST DATA #1 */
  /* call Dirichet() with valid data (N=even) */
  parameters.n      = 10;
  parameters.length = 101;
  parameters.deltaX = 0.01;


  LALSCreateVector (&status, &poutput, parameters.length);
  if ( check( &status, 0 , "") )
  {
    return 1;
  }


  LALDirichlet( &status, poutput, &parameters );
  if ( check( &status, 0 , "") )
  {
    return 1;
  }


  LALPrintVector(poutput);
  if ( check( &status, 0 , "") )
  {
    return 1;
  }


  LALSDestroyVector( &status, &poutput );
  if ( check( &status, 0 , "") )
  {
    return 1;
  }

  /* VALID TEST DATA #2 */
  /* call Dirichet() with valid data (N=odd) */
  parameters.n      = 11;
  parameters.length = 101;
  parameters.deltaX = 0.01;


  LALSCreateVector(&status, &poutput, parameters.length);
  if ( check( &status, 0 , "") )
  {
    return 1;
  }


  LALDirichlet( &status, poutput, &parameters );
  if ( check( &status, 0 , "") )
  {
    return 1;
  }


  LALPrintVector(poutput);
  if ( check( &status, 0 , "") )
  {
    return 1;
  }


  LALSDestroyVector( &status, &poutput );
  if ( check( &status, 0 , "") )
  {
    return 1;
  }


  /* VALID TEST DATA #3 */
  /* call Dirichet() with valid data (x=0 to 2) */
  parameters.n      = 10;
  parameters.length = 201;
  parameters.deltaX = 0.01;


  LALSCreateVector(&status, &poutput, parameters.length);
  if ( check( &status, 0 , "") )
  {
    return 1;
  }


  LALDirichlet( &status, poutput, &parameters );
  if ( check( &status, 0 , "") )
  {
    return 1;
  }


  LALPrintVector(poutput);
  if ( check( &status, 0 , "") )
  {
    return 1;
  }


  LALSDestroyVector( &status, &poutput );
  if ( check( &status, 0 , "") )
  {
    return 1;
  }


  return 0;
}
/*------------------------------------------------------------------------*/

int
check( LALStatus* status, INT4 code, const CHAR* message )
{
  if ( status->statusCode!= code )
  {
    printf ( "FAIL: did not recognize \"%s\"\n", message );
    return 1;
  }
  else if (code && strcmp( message, status->statusDescription)) {
    printf("FAIL: incorrect warning message \"%s\" not \"%s\"\n",
	   status->statusDescription, message);

    return 1;
  }

  return 0;
}
/** \endcond */
