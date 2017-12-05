/*
*  Copyright (C) 2007 B.S. Sathyaprakash
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

/*
 * Given the noise amplitude/power spectral density as a two-column
 * ascii file this code returns the distance in Mpc as a function of the
 * total mass to binaries with equal masses and an rms orientation 
 * and produce an SNR of 5 when integrated up to the LSO at R=6M
 * or the light ring at R=2.8 M.
 *
 */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
int main (int argc, char **argv)
{

	double totalMass, eta, fOld, x, sevenBy3=7./3.;
	double pwrSpec, LO;
	double freq;
	double rhorms, rmsD, idealD, Crms, integral;  
        int    ampSpec;
	double lso_radius=6., light_ring_radius=2.8, MTSUN_SI=4.925e-6, pi=3.14159, MPC_SI=3.e22, C_SI=3.e8;
   

	if (argc!=3)
	{
                fprintf(stderr, "* Given the noise amplitude/power spectral density as a two-column\n");
                fprintf(stderr, "* ascii file this code returns the distance in Mpc as a function of the\n");
                fprintf(stderr, "* total mass to binaries with equal masses and an rms orientation \n");
                fprintf(stderr, "* and produce an SNR of 5 when integrated up to the LSO R=6M or \n");
                fprintf(stderr, "* the light ring R=2.8 M.\n");
		fprintf(stderr, "-----------------------------\n");
		fprintf (stderr, "Usage: %s ampSpec (0/1) LO (0/1)\n", argv[0]);
		fprintf (stderr, "ampSpec: should be 0 if input spectrum is power and 1 if amplitude\n");
		fprintf (stderr, "LO: should be 0 if the last orbit is at R=6M (Schwarzschild) and nonzero if at R=2.8M (as in EOB)\n");
		fprintf(stderr, "-----------------------------\n");
		exit(1);
	}

	/* Cast the inputs */
	ampSpec = atoi(argv[1]);

	if (atoi(argv[2]) == 0)
		LO = lso_radius;
	else
		LO = light_ring_radius;

	/* Equal masses means eta=1/4 */
	eta = 0.25;
	/* rms SNR = 5 */
	rhorms = 5.;
	/* initiate the SNR integral to zero */
	integral = 0.;

	/* If the file is empty write a warning message and exit */
	if ( scanf("%le %le\n", &freq, &x) == EOF)
	{
		fprintf(stderr, "Empty input file\n");
		exit(1);
	}
	fOld = freq;
	
	while (scanf("%le %le\n", &freq, &x) != EOF)
	{
		if (ampSpec)
			/* If the file contains amplitude spectrum convert it to power spectrum */
			pwrSpec = x*x;
		else
			/* If the file contains power spectrum use as is */
			pwrSpec = x;

		integral += (freq-fOld) * pow(freq, -sevenBy3)/pwrSpec;
			
		/* Compute the total mass of a binary whose last orbit is the current frequency */
	        totalMass = 1./(pow(LO, 1.5) * pi * freq);
		/* Normalisation to covert integral above to a distance */
		Crms = pow(totalMass, 5./6.) * pow(2.*eta/15., 0.5) / (rhorms * pow(pi, 2./3.));
		/* Distance corresponding to a source with RMS orientation producing an SNR equal to rhorms */
		rmsD = Crms * sqrt(integral) / MPC_SI * C_SI;
		/* Distance corresponding to an optimally oriented source with the same SNR */
		idealD = 2.5 * rmsD;
		/* Convert the mass in seconds to solar masses */
		totalMass /= MTSUN_SI;
	        if (totalMass < 5000.) fprintf(stdout, "%e %e %e %e\n", totalMass, rmsD, idealD, freq);
		fOld = freq;
	}
	fprintf(stdout, "&\n");
	return 0;
}
