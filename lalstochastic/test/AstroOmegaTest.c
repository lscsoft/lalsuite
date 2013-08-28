/*
*  Copyright (C) 2007 Jolien Creighton, Tania Regimbau
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
 * \author Regimbau Tania
 * \file
 * \ingroup AstroOmega_h
 *
 * \brief This programs verifies that the routine LALAstroOmega() gives the expected results (computed separetly with mathematica) for a set of input parameters.
 *
 * \heading{Exit codes}
 * returns 0 on success, otherwise returns 1.
 *
 * \heading{Uses}
 * \code
 * lalDebugLevel
 * LALAstroOmega()
 * \endcode
 *
 */

#include <stdio.h>
#include <math.h>
#include <lal/LALConfig.h>
#include <lal/LALStdlib.h>
#include <lal/Integrate.h>
#include <lal/AstroOmega.h>

static void SDensity (REAL8 *dEgw, REAL8 nu)
 {
  *dEgw=9.27E35 * pow(nu,3.);
  return;
 }

int main ( void )
 {
  static LALStatus status;
  AstroOmegaParams p;
  AstroOmegaSourceParams sourcep;
  AstroOmegaCosmoParams cosmop;
  REAL8 nu, test;


  cosmop.ho = 0.65;
  cosmop.density_matter = 0.3;
  cosmop.density_vacuum = 0.7;
  cosmop.density_k = 0.;
  nu = 1100.;
  /*source*/
  sourcep.lambda = 4.84E-3;
  sourcep.numax = 4000.;
  sourcep.SDensitySource = SDensity;
  p.cosmoparams = cosmop;
  p.sourceparams = sourcep;


  p.extraparams = &nu;
  LALAstroOmega (&status, &test, nu,&p);
  if (fabs(test-1.16E-10) > 1.E-12)
   {
    printf("error! the right value is 1.16E-10 not %.2e\n",test);
   }
  else printf("omega(%f)= %.2e o.k\n", nu,test);

   /*
  for (nu = 0.;nu < sourcep.numax;nu = nu + 10.)
   {
    p.extraparams = &nu;
    LALAstroOmega (&status, &test, nu, &p);
    printf("%f\t %.2e\n", nu,test);}
  */
  return 0;
  }
