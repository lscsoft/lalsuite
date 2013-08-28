/*
*  Copyright (C) 2007 Bernd Machenschalk, David Churches, Duncan Brown, Jolien Creighton, B.S. Sathyaprakash
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
 * \author Sathyaprakash, B. S.
 * \ingroup LALNoiseModels_h
 * \file
 *
 * \brief This program can be used generate expected noise
 * NoiseSpectralDensity in various interferometers.
 * See the beginning of the NoiseModels module to see details on how
 * this test program works.
 *
 * ### Uses ###
 *
 * \code
 * LALDCreateVector
 * LALNoiseSpectralDensity
 * LALGEOPsd
 * LALLIGOIPsd
 * LALTAMAPsd
 * LALVIRGOPsd
 * LALAdvLIGOPsd
 * LALDDestroyVector
 * LALCheckMemoryLeaks
 * \endcode
 *
 */
#include <lal/AVFactories.h>
#include <lal/LALNoiseModels.h>

/** \cond DONT_DOXYGEN */

/** \endcond */

int main ( void )
{

   static LALStatus status;
   REAL8Vector *psd=NULL;
   REAL8       df;
   INT4 i, length;
   FILE *NoisePsdFile;

   fprintf(stderr, "This test code computes the amplitude spectrum \n");
   fprintf(stderr, "of GEO, LIGO, VIRGO, TAMA and AdvLIGO and writes them in \n");
   fprintf(stderr, "NoisePSDTest.out in a format suitable for \n");
   fprintf(stderr, "display with xmgr/xgrace\n");

   if ( (NoisePsdFile = fopen("NoisePSDTest.out", "w")) == NULL)
   {
      fprintf(stderr, "Can't open output file\n");
      exit(0);
   }
   df = 1.0;
   length = 8193;

   LALDCreateVector(&status, &psd, length);
   fprintf(stderr, "Length of vector=%d\n", psd->length);

   LALNoiseSpectralDensity(&status, psd, &LALGEOPsd, df);

   for (i=2; i<length; i++) {
        if (psd->data[i]) fprintf (NoisePsdFile, "%d %e\n", i, sqrt(psd->data[i]));
   }
   fprintf(NoisePsdFile, "&\n");
   LALNoiseSpectralDensity(&status, psd, &LALLIGOIPsd, df);
   for (i=2; i<length; i++) {
        if (psd->data[i]) fprintf (NoisePsdFile, "%d %e\n", i, sqrt(psd->data[i]));
   }
   fprintf(NoisePsdFile, "&\n");
   LALNoiseSpectralDensity(&status, psd, &LALVIRGOPsd, df);
   for (i=2; i<length; i++) {
        if (psd->data[i]) fprintf (NoisePsdFile, "%d %e\n", i, sqrt(psd->data[i]));
   }
   fprintf(NoisePsdFile, "&\n");
   LALNoiseSpectralDensity(&status, psd, &LALTAMAPsd, df);
   for (i=2; i<length; i++) {
        if (psd->data[i]) fprintf (NoisePsdFile, "%d %e\n", i, sqrt(psd->data[i]));
   }
   fprintf(NoisePsdFile, "&\n");
   LALNoiseSpectralDensity(&status, psd, &LALAdvLIGOPsd, df);
   for (i=2; i<length; i++) {
        if (psd->data[i]) fprintf (NoisePsdFile, "%d %e\n", i, sqrt(psd->data[i]));
   }
   LALDDestroyVector(&status, &psd);
   LALCheckMemoryLeaks();
   return 0;
}
