/*
*  Copyright (C) 2007 Bernd Machenschalk, Chris Van Den Broeck
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
 * \author Van Den Broeck, C and Jones, G.
 * \file
 *
 * \brief Tests InpsiralBCVSpinBankTest().
 *
 * \heading{Usage}
 *
 * \code
 * InspiralBCVSpinBankTest
 * \endcode
 *
 * This program uses InspiralBCVSpinBank() to generate a template bank from
 * command line parameter input.
 *
 * \heading{Command line options}
 *
 * \heading{Notes}
 *
 */

/**\name Error Codes */ /*@{*/
#define INSPIRALBCVSPINBANKTESTC_ENORM     0
#define INSPIRALBCVSPINBANKTESTC_EMEM      1
#define INSPIRALBCVSPINBANKTESTC_ESUB      2
#define INSPIRALBCVSPINBANKTESTC_EFILE     4

#define INSPIRALBCVSPINBANKTESTC_MSGENORM  "Normal exit"
#define INSPIRALBCVSPINBANKTESTC_MSGEMEM   "Memory allocation error"
#define INSPIRALBCVSPINBANKTESTC_MSGESUB   "Subroutine error"
#define INSPIRALBCVSPINBANKTESTC_MSGEFILE  "File I/O error"
/*@}*/

#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALGSL.h>
#include <gsl/gsl_matrix.h>

void testfunction(LALStatus *status);

int main (void)
{
  static LALStatus status;
  testfunction(&status);
  return 0;
}

void testfunction(LALStatus *status)
{
  int i, j;
  float x = 0;
  gsl_matrix *m = NULL;

  ATTATCHSTATUSPTR( status);
  CALLGSL((m = gsl_matrix_alloc (10, 3)), status);
  for (i = 0; i < 10; i++)
    for (j = 0; j < 3; j++)
      CALLGSL ( gsl_matrix_set (m, i, j, 0.23 + 100*i + j), status);

  for (i = 0; i < 10; i++)
    for (j = 0; j < 3; j++)
{
      CALLGSL((x = gsl_matrix_get (m, i, j)), status);
      printf ("m(%d,%d) = %g\n", i, j, x);
}

  CHECKSTATUSPTR (status);
  REPORTSTATUS (status);
  DETATCHSTATUSPTR( status);
  RETURN (status);
}


