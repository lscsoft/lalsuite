/*
*  Copyright (C) 2007 Jolien Creighton
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
   \ingroup PrintVector_h

   \brief Tests the routines in \ref PrintVector.c.

   \heading{Usage}
   \code
PrintVectorTest
   \endcode

   \heading{Description}

This program generates and prints a sequence of
   <tt>\<datatype\></tt>Vectors; the program itself always
returns success, so the testing function is actually served by
examinaton of the output files.

\heading{Exit codes}
<table>
<tr><th>Code</th><th>Explanation</th></tr>
<tr><td>0</td><td>Always returned.</td></tr>
/table>

*/
/** \cond DONT_DOXYGEN */

#define LAL_USE_OLD_COMPLEX_STRUCTS

#ifndef _PRINTVECTOR_H
#include <lal/PrintVector.h>
#ifndef _PRINTVECTOR_H
#define _PRINTVECTOR_H
#endif
#endif

#ifndef _MATH_H
#include <math.h>
#ifndef _MATH_H
#define _MATH_H
#endif
#endif

#ifndef _LALSTDLIB_H
#include <lal/LALStdlib.h>
#ifndef _LALSTDLIB_H
#define _LALSTDLIB_H
#endif
#endif

#ifndef _AVFACTORIES_H
#include <lal/AVFactories.h>
#ifndef _AVFACTORIES_H
#define _AVFACTORIES_H
#endif
#endif


int main( void )
{
  static LALStatus status;

  COMPLEX16Vector   *zVector;
  COMPLEX8Vector    *cVector;
  REAL8Vector       *dVector;
  REAL4Vector       *sVector;
  INT2Vector        *i2Vector;
  INT4Vector        *i4Vector;
  INT8Vector        *i8Vector;
  UINT2Vector       *u2Vector;
  UINT4Vector       *u4Vector;
  UINT8Vector       *u8Vector;
  CHARVector        *charVector;

  COMPLEX16        *z;
  COMPLEX8         *c;
  REAL8            *d;
  REAL4            *s;
  INT2             *i2;
  INT4             *i4;
  INT8             *i8;
  UINT2             *u2;
  UINT4             *u4;
  UINT8             *u8;
  CHAR              *ch;
  INT2             n;

  zVector = NULL;

  LALZCreateVector( &status, &zVector, 8 );
  for ( n=zVector->length, z=zVector->data; n > 0 ; --n, ++z ) {
    z->real_FIXME = sinh(90.0*(4-n));
    z->imag_FIXME = - 1 / (1e-300 + creal(*z));
  }
  LALZPrintVector(zVector);

  cVector = NULL;

  LALCCreateVector( &status, &cVector, 8 );
  for ( n=cVector->length, c=cVector->data; n > 0 ; --n, ++c ) {
    c->realf_FIXME = sinh(9.0*(4-n));
    c->imagf_FIXME = - 1 / (1e-30 + crealf(*c));
  }
  LALCPrintVector(cVector);

  dVector = NULL;

  LALDCreateVector( &status, &dVector, 8 );
  for ( n=dVector->length, d=dVector->data; n > 0 ; --n, ++d ) {
    *d = sinh(90.0*(4-n));
  }
  LALDPrintVector(dVector);
  for ( n=dVector->length, d=dVector->data; n > 0 ; --n, ++d ) {
    *d = 1 / (1e-300 + *d);
  }
  LALDPrintVector(dVector);

  sVector = NULL;

  LALSCreateVector( &status, &sVector, 8 );
  for ( n=sVector->length, s=sVector->data; n > 0 ; --n, ++s ) {
    *s = sinh(9.0*(4-n));
  }
  LALSPrintVector(sVector);
  for ( n=sVector->length, s=sVector->data; n > 0 ; --n, ++s ) {
    *s = 1 / (1e-30 + *s);
  }
  LALSPrintVector(sVector);

  i2Vector = NULL;

  LALI2CreateVector( &status, &i2Vector, 10 );
  for ( n=i2Vector->length, i2=i2Vector->data; n > 0 ; --n, ++i2 ) {
    *i2 = pow((n-5),1+abs(n-5));
  }
  LALI2PrintVector(i2Vector);

  i4Vector = NULL;

  LALI4CreateVector( &status, &i4Vector, 10 );
  for ( n=i4Vector->length, i4=i4Vector->data; n > 0 ; --n, ++i4 ) {
    *i4 = pow((n-5),1+2*abs(n-5));
  }
  LALI4PrintVector(i4Vector);

  i8Vector = NULL;

  LALI8CreateVector( &status, &i8Vector, 10 );
  for ( n=i8Vector->length, i8=i8Vector->data; n > 0 ; --n, ++i8 ) {
    *i8 = pow((n-5),1+abs(n-5)*abs(n-5));
	/* *i8 = pow(10,2*n); */
  }
  LALI8PrintVector(i8Vector);

  u2Vector = NULL;

  LALU2CreateVector( &status, &u2Vector, 10 );
  for ( n=u2Vector->length, u2=u2Vector->data; n > 0 ; --n, ++u2 ) {
    *u2 = pow(abs(n-5),1+abs(n-5));
  }
  LALU2PrintVector(u2Vector);

  u4Vector = NULL;

  LALU4CreateVector( &status, &u4Vector, 10 );
  for ( n=u4Vector->length, u4=u4Vector->data; n > 0 ; --n, ++u4 ) {
    *u4 = pow(abs(n-5),1+2*abs(n-5));
  }
  LALU4PrintVector(u4Vector);

  u8Vector = NULL;

  LALU8CreateVector( &status, &u8Vector, 10 );
  for ( n=u8Vector->length, u8=u8Vector->data; n > 0 ; --n, ++u8 ) {
    *u8 = pow(abs(n-5),1+abs(n-5)*abs(n-5));
	/* *u8 = pow(10,2*n); */
  }
  LALU8PrintVector(u8Vector);


  charVector = NULL;

  LALCHARCreateVector( &status, &charVector, 10 );
  ch = charVector->data;
  *(ch++) = 'a';
  *(ch++) = 'A';
  *(ch++) = '%';
  *(ch++) = '@';
  *(ch++) = '2';
  *(ch++) = 'e';
  *(ch++) = '+';
  *(ch++) = '\a';
  *(ch++) = 127;
  *(ch++) = '\0';

  LALCHARPrintVector(charVector);

  return 0;
}

/** \endcond */
