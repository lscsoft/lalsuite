/*
*  Copyright (C) 2011 Drew Keppel
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
\author Keppel, D.

\brief Create waveforms based on the TaylorT2 model.

Outputs a two files with three columns corresponding to time (in seconds),
\f$h_+\f$, and \f$h_\times\f$. The first output file is <tt>T2wave1.dat</tt> in
the current directory. This is generated using the
XLALSimInspiralTaylorT2PNRestricted routine. The seond output file is
<tt>T2wave2.dat</tt> in the current directory. This is generated using the
LALInspiralWave3Templates routine.

\heading{Usage}

<tt>LALInspiralTaylor22Test</tt>

*/

#include <time.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/LALSimInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/GenerateInspiral.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

int main(void) {
	static LALStatus    mystatus;

	/* variables for timing purposes */
	clock_t start, diff;
	int msec;

	REAL4TimeSeries *h_plus;
	REAL4TimeSeries *h_cross;
	REAL8TimeSeries *hplus;
	REAL8TimeSeries *hcross;
	InspiralTemplate params;

	FILE *outputfile;
	INT4 i,length,O;
	REAL8 dt, m, m1, m2, mu, nu, lambda1, lambda2;
	LIGOTimeGPS tc = LIGOTIMEGPSZERO;

	memset( &mystatus, 0, sizeof(LALStatus) );
	memset( &params, 0, sizeof(InspiralTemplate) );

	m1 = 5.;
	m2 = 5.;
	m = m1 + m2;
	mu = m1 * m2 / m;
	nu = mu / m;

	lambda1 = 0.;
	lambda2 = 0.;
	LALSimInspiralWaveformFlags *waveFlags;
	waveFlags = XLALSimInspiralCreateWaveformFlags();
	XLALSimInspiralSetInteraction( waveFlags, XLALGetInteractionFromString( "ALL" ) );

	O = 7;
	switch (O)
	{
		case 0:
		case 1:
			O = 0;
			params.order = LAL_PNORDER_NEWTONIAN;
			break;
		case 2:
			params.order = LAL_PNORDER_ONE;
			break;
		case 3:
			params.order = LAL_PNORDER_ONE_POINT_FIVE;
			break;
		case 4:
			params.order = LAL_PNORDER_TWO;
			break;
		case 5:
			params.order = LAL_PNORDER_TWO_POINT_FIVE;
			break;
		case 6:
			params.order = LAL_PNORDER_THREE;
			break;
		case 7:
			params.order = LAL_PNORDER_THREE_POINT_FIVE;
			break;
	}

	params.approximant = TaylorT2;
	params.ampOrder = LAL_PNORDER_NEWTONIAN;
	params.mass1 = m1;
	params.mass2 = m2;
	params.totalMass = m;
	params.mu = mu;
	params.eta = nu;
	params.tSampling = 4096;
	params.fCutoff = 2047.;
	params.fLower = 40.;
	params.distance = 1e6 * LAL_PC_SI;
	params.signalAmplitude = 4.0 * params.mu * LAL_MRSUN_SI/params.distance;
	params.ieta = 1;

	

	dt = 1. / params.tSampling;

	start = clock();
	length = XLALSimInspiralTaylorT2PNRestricted(&hplus, &hcross, 0., dt, params.mass1*LAL_MSUN_SI, params.mass2*LAL_MSUN_SI, params.fLower, 0., params.distance, lambda1, lambda2, XLALSimInspiralGetInteraction(waveFlags), 0, O);
	diff = clock() - start;
	msec = diff * 1000 / CLOCKS_PER_SEC;
	printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);

	fprintf(stderr, "length = %i\n", length);
	fprintf(stderr, "T = %f\n", (float) length * dt);

	outputfile = fopen("T2wave1.dat","w");

	length = hplus->data->length;

	for(i = 0; i < length; i++) {
		fprintf(outputfile,"%e\t%e\t%e\n",
			i*dt,
			hplus->data->data[i],
			hcross->data->data[i]);
	}

	fclose(outputfile);
	fprintf(stdout,"waveform saved in T2wave1.dat\n" );

	start = clock();
	h_plus = XLALCreateREAL4TimeSeries("", &tc, 0.0, dt, &lalDimensionlessUnit, params.tSampling*3);
	h_cross = XLALCreateREAL4TimeSeries("", &tc, 0.0, dt, &lalDimensionlessUnit, params.tSampling*3);

	fprintf(stderr, "Lower cut-off frequency used will be %fHz\n", params.fLower);

	/* --- now we can call the injection function --- */
	XLALInspiralWave2Templates(h_plus->data, h_cross->data, &params);
	diff = clock() - start;
	msec = diff * 1000 / CLOCKS_PER_SEC;
	printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);

	if ( mystatus.statusCode )
	{
		fprintf( stderr, "LALInspiralTaylorT2Test: error generating waveform\n" );
		exit( 1 );
	}

	/* --- and finally save in a file --- */

	outputfile = fopen("T2wave2.dat","w");

	length = h_plus->data->length;

	for(i = 0; i < length; i++) {
		fprintf(outputfile,"%e\t%e\t%e\n",
			i*dt,
			h_plus->data->data[i],
			h_cross->data->data[i]);
	}

	fclose(outputfile);
	fprintf(stdout,"waveform saved in T2wave2.dat\n" );
	return 0;
}
