#include <stdlib.h>
#include <math.h>
#include <stdio.h>

int main (int argc, char **argv)
{
	double m1, m2, totalMass, eta, flso, lso=2.8, fOld, x, sevenBy3 = 7.L/3.L;
	double pwrSpec;
	double freq, flow;
	double rhorms, rmsD, idealD, Crms, integral;  
	double MTSUN_SI=4.925e-6, pi=3.14159, PC_SI=3.e19, C_SI=3.e8;
   
	if (argc !=5 )
	{
		fprintf(stderr, "\n*****************************\n");
		fprintf(stderr, "Usage:%s mass1/solarMass mass2/solarMass rhoRMS fLow < PSD_File\n",argv[0]);
		fprintf(stderr, "*****************************\n");
		exit(1);
	}

	m1 = atof(argv[1]);
	m2 = atof(argv[2]);
	rhorms = atof(argv[3]);
	flow = atof(argv[4]);

	totalMass = (m1+m2);
	eta = m1*m2/(totalMass*totalMass);
	totalMass *= MTSUN_SI;
	Crms = pow(totalMass,5./6.) * pow(2.*eta/15.,0.5) / (rhorms * pow(pi,2./3.));

	integral = 0.;
	scanf("%le %le\n", &freq, &x);
	fOld = freq;

	while ( scanf("%le %le\n", &freq, &x) != EOF)
	{
	
		flso = 1.L/ ( pow(lso, 1.5) * pi * totalMass);
		if (freq > flso) break;
		if (x && freq > flow)
		{
			pwrSpec = x;
			integral += (freq - fOld) * pow(freq, -sevenBy3)/pwrSpec;
		}
		fOld = freq;
	}

	rmsD = Crms * sqrt(integral) / PC_SI * C_SI;
	idealD = 2.5 * rmsD;
	fprintf(stdout, "%e %e %e %e %e %e\n", flow, rmsD, totalMass/MTSUN_SI, idealD, m1, m2, flow);
	return 0;
}
