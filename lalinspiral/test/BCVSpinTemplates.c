/*
*  Copyright (C) 2007 B.S. Sathyaprakash, Thomas Cokelaer
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
 * \author B.S. Sathyaprakash
 * \file
 *
 * \brief Creates a template mesh for BCVSpin using the mismatch metric.
 *
 * ### Usage ###
 *
 *
 * ### Description ###
 *
 *
 * ### Algorithm ###
 *
 *
 * ### Uses ###
 *
 * \code
 * lalDebugLevel
 * \endcode
 *
 * ### Notes ###
 *
 */

#include <math.h>
#include <stdlib.h>
#include <lal/LALInspiralBank.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>

/* Default parameter settings. */

int
main(int argc, char **argv)
{
  /* top-level status structure */
 UINT4   numPSDpts=8193;
  static LALStatus status;
  static InspiralCoarseBankIn coarseIn;
  void *noisemodel = LALLIGOIPsd;
  double beta;
  SnglInspiralTable *tiles=NULL, *first=NULL;
  FILE *fpr;

/* Number of templates is nlist */
  INT4 nlist1, j;


  fpr = fopen("BCVSpinTemplates.out", "w");
  nlist1 = 0;
  coarseIn.HighGM = 6.;
  coarseIn.LowGM = -4.;
  coarseIn.fLower = 40.L;
  coarseIn.fUpper = 400.L;
  coarseIn.tSampling = 4096.L;
  coarseIn.order = twoPN;
  coarseIn.mmCoarse = 0.90;
  coarseIn.mmFine = 0.90;
  coarseIn.iflso = 0.0L;
  coarseIn.mMin = 3.0;
  coarseIn.mMax = 20.0;
  coarseIn.MMax = coarseIn.mMax * 2.;
  coarseIn.massRange = MinMaxComponentMass;
  /* coarseIn.massRange = MinComponentMassMaxTotalMass;*/
  /* minimum value of eta */
  coarseIn.etamin = coarseIn.mMin * ( coarseIn.MMax - coarseIn.mMin) / pow(coarseIn.MMax,2.);
  coarseIn.psi0Min = 1.0e4;
  coarseIn.psi0Max = 6.0e4;
  coarseIn.psi3Min = -5.0e2;
  coarseIn.psi3Max = 1.0e1;
  coarseIn.alpha = 0.L;
  coarseIn.numFcutTemplates = 1;
  coarseIn.betaMin = 10.0;
  coarseIn.betaMax = 700.;
  coarseIn.spinBank = 2;
  coarseIn.iseed = 9295883;
  coarseIn.nTIni = 10000;
  coarseIn.ShMaxSz = 1024;
  coarseIn.gridSpacing = Hexagonal;
  coarseIn.insidePolygon = True;

  memset( &(coarseIn.shf), 0, sizeof(REAL8FrequencySeries) );
  coarseIn.shf.f0 = 0;
  LALDCreateVector( &status, &(coarseIn.shf.data), numPSDpts );
  coarseIn.shf.deltaF = coarseIn.tSampling / (2.*(REAL8) coarseIn.shf.data->length + 1.L);
  LALNoiseSpectralDensity (&status, coarseIn.shf.data, noisemodel, coarseIn.shf.deltaF );

  coarseIn.approximant = BCVSpin;
  coarseIn.space       = Psi0Psi3;


  LALInspiralBCVSpinRandomBank (&status, &tiles, &nlist1, &coarseIn);
/*
  LALInspiralBankGeneration(&status, &coarseIn, &tiles, &nlist1);
*/
  LALDDestroyVector( &status, &(coarseIn.shf.data) );

  fprintf (fpr, "#numtemplaes=%d %e %e %e %e %e %e\n", nlist1, coarseIn.psi0Min, coarseIn.psi0Max, coarseIn.psi3Min, coarseIn.psi3Max, coarseIn.betaMin, coarseIn.betaMax);
  beta = tiles->beta;
  first = tiles;
  for (j=0; j<nlist1; j++)
  {
	  fprintf(fpr, "%7.3f %e %e\n", beta, tiles->psi0, tiles->psi3);
	  tiles = tiles->next;
	  if (tiles != NULL && beta != tiles->beta)
	  {
		  beta = tiles->beta;
	  }
  }
  fclose(fpr);

  tiles = first;
  for (j=0; j<nlist1; j++)
  {
          first = tiles ->next;
	  if (tiles != NULL) LALFree(tiles);
          tiles = first;
  }
  LALCheckMemoryLeaks();
  return 0;
}
