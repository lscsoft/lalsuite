/* This programme expects six inputs at the command line:
 *
 * m1: mass of the first object in the binary source,
 * m2: mass of the second object in the binary source,
 * rhorms: rms SNR required from the source,
 * flow: lower frequency cutoff,
 * LO: Last orbit in units of M (either 6, for LSO, or 2.8, for light ring)
 * ispec: an integer indicating the nature of input data: power (0) or amplitude (1) spectrum is input
 *
 * reads a data from a file contianing two columns of frequency
 * and power/amplitude spectral density. It outputs the distance in Kpc
 * for a source with the given rms SNR. 
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

int main (int argc, char **argv)
{
	double m1, m2, totalMass, eta, flso, lso, fOld, sevenBy3 = 7.L/3.L;
	double x, pwrSpec;
	double freq, flow;
	double rhorms, rmsD, idealD, Crms, integral;  
	double MTSUN_SI=4.925e-6, pi=3.14159, PC_SI=3.e19, C_SI=3.e8;
	int ispec;
   
	if (argc !=7 )
	{
		fprintf(stderr, "-----------------------------\n");
		fprintf(stderr, "Usage:%s mass1/solarMass mass2/solarMass rhoRMS fLow LO ispec < PSD_File\n",argv[0]);
		fprintf(stderr, "-----------------------------\n");
		exit(1);
	}

	m1 = atof(argv[1]);
	m2 = atof(argv[2]);
	rhorms = atof(argv[3]);
	flow = atof(argv[4]);
	lso = atof(argv[5]);
	ispec = atoi(argv[6]);

	totalMass = (m1+m2);
	eta = m1*m2/(totalMass*totalMass);
	totalMass *= MTSUN_SI;
	Crms = pow(totalMass,5./6.) * pow(2.*eta/15.,0.5) / (rhorms * pow(pi,2./3.));

	flso = 1.L/ ( pow(lso, 1.5) * pi * totalMass);
	integral = 0.;
	if (scanf("%le %le\n", &freq, &x) == EOF) 
	{
		fprintf(stderr, "No data in file\n");
		return 1;
	}
	fOld = freq;
		
	while ( scanf("%le %le\n", &freq, &x) != EOF)
	{
	
		if (freq > flso) break;
		if (x && freq > flow)
		{
			if (ispec)
				pwrSpec = x*x;
			else
				pwrSpec = x;
			integral += (freq - fOld) * pow(freq, -sevenBy3)/pwrSpec;
		}
		fOld = freq;
	}

	rmsD = Crms * sqrt(integral) / PC_SI * C_SI;
	idealD = 2.5 * rmsD;
	fprintf(stdout, "%e %e %e %e %e %e\n", flow, rmsD, totalMass/MTSUN_SI, idealD, m1, m2);
	return 0;
}
