#include <stdio.h>
#include <math.h>
#include <lal/LALDetectors.h>
#include <lal/DetResponse.h>


struct detectors {
	LALDetector *detector;
	int n;
};


static double point_source_average(struct detectors detectors, double ra, double dec)
{
	const int N = 1000;
	double integral = 0;
	int i, j;
	double gmst;

	for(gmst = i = 0; gmst < LAL_TWOPI; gmst = ++i * LAL_TWOPI / N) {
		double fsquared = 1;

		for(j = 0; j < detectors.n; j++) {
			double fp, fc;
			XLALComputeDetAMResponse(&fp, &fc, (const REAL4(*)[3])detectors.detector[j].response, ra, dec, 0, gmst);
			fsquared *= fp * fp + fc * fc;
		}

		integral += fsquared * LAL_TWOPI / N;
	}

	return integral / LAL_TWOPI;
}


static double all_sky_average(struct detectors detectors)
{
	const int N = 1000;
	double integral = 0;
	double ra, dec;
	int i, j, k;

	for(ra = i = 0; ra < LAL_TWOPI; ra = ++i * LAL_TWOPI / N)
		for(dec = (-N / 2 + (j = 0)) * LAL_PI / N; dec < LAL_PI_2; dec = (-N / 2 + ++j) * LAL_PI / N) {
			double fsquared = 1;

			for(k = 0; k < detectors.n; k++) {
				double fp, fc;
				XLALComputeDetAMResponse(&fp, &fc, (const REAL4(*)[3]) detectors.detector[k].response, ra, dec, 0, 0);
				fsquared *= fp * fp + fc * fc;
			}

			integral += fsquared * sin(LAL_PI_2 + dec) * (LAL_PI / N) * (LAL_TWOPI / N);
		}

	return integral / (4 * LAL_PI);
}


int main(void)
{
	struct detectors detectors;
	/* galactic core */
	double ra = 2.0318570464121519;
	double dec = -0.50628171572274738;
	/* sn1987a */
	/*double ra = 1.4637516476674761;
	double dec = -1.2089885547477417;*/
	int i;

	detectors.n = 3;
	detectors.detector = malloc(detectors.n * sizeof(*detectors.detector));

	detectors.detector[0] = lalCachedDetectors[LAL_LHO_4K_DETECTOR];
	detectors.detector[1] = lalCachedDetectors[LAL_LHO_2K_DETECTOR];
	detectors.detector[2] = lalCachedDetectors[LAL_LLO_4K_DETECTOR];

	printf("\nAverage of product of (F+^2 + Fx^2)\n-----------------------------------\n");
	printf("\nDetectors: ");
	for(i = 0; i < detectors.n; i++)
		printf(" %s", detectors.detector[i].frDetector.prefix);

	printf("\npoint source @ ra = %.16g, dec = %.16g:\n\t%.16g\n", ra, dec, point_source_average(detectors, ra, dec));
	printf("all sky:\n\t%.16g\n\n", all_sky_average(detectors));

	return 0;
}
