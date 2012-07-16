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

\brief Create waveforms based on the chosen model from lalsimulation and lalinspiral routines.

Outputs a two files with three columns corresponding to time (in seconds),
\f$h_+\f$, and \f$h_\times\f$. The first output file is <tt>wave1.dat</tt> in
the current directory. This is generated using the lalsimulation routine. The
seond output file is <tt>wave2.dat</tt> in the current directory. This is
generated using the lalinspiral routine.

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

typedef struct{
  INT4 order;
  char waveformString[LIGOMETA_WAVEFORM_MAX];
} OtherParamIn;

static void
ParseParameters( UINT4            argc,
                 char             **argv,
                 OtherParamIn     *otherIn)
{
  UINT4 i = 1;

  while(i < argc)
    {
      if( strcmp(argv[i],"--approximant") == 0 )
      {
        strcpy( otherIn->waveformString, argv[++i] );
      }
      else if( strcmp(argv[i],"--order") == 0 )
      {
        otherIn->order = atoi(argv[++i]);
      }
      i++;
    }
}

int main(int argc , char **argv)
{
	/* set us up to fail hard */
	lalDebugLevel = 7;
	XLALSetErrorHandler(XLALAbortErrorHandler);
	static LALStatus    mystatus;

	/* variables for timing purposes */
	clock_t start, diff;
	int msec;

	REAL4TimeSeries *h_plus;
	REAL4TimeSeries *h_cross;
	REAL8TimeSeries *hplus;
	REAL8TimeSeries *hcross;
	InspiralTemplate params;
	InspiralInit paramsInit;

	const char *name = "wave1.dat";
	FILE *outputfile;
	int i,length;
	REAL8 dt;
	LIGOTimeGPS tc = LIGOTIMEGPSZERO;
	/* The two tidal lambda's and the WaveformFlags below can't be
	 * set from the command line. However, they aren't used in the old
	 * code, so to sanity check old vs new code, this is OK.
	 * We just set the lambda's to zero and create a WaveformFlags struct
	 * with default values (in particular, it turns on all interactions).
	 */  
	REAL8 lambda1 = 0., lambda2 = 0., fRef = 0.;
	LALSimInspiralWaveformFlags *waveFlags=XLALSimInspiralCreateWaveformFlags();
	LALSimInspiralTestGRParam *nonGRparams = NULL;

	memset( &mystatus, 0, sizeof(LALStatus) );
	memset( &params, 0, sizeof(InspiralTemplate) );

	static OtherParamIn otherIn; /* some extra parameters to parse */
        otherIn.order = 7; /* default order */
        strncpy(otherIn.waveformString, "TaylorT3", sizeof(otherIn.waveformString)); /* default approximant */
	ParseParameters(argc, argv, &otherIn);/* let's parse user parameters */
	LALInspiralITStructureSetDefault(&mystatus, &params);
	LALInspiralITStructureParseParameters(&mystatus, argc, argv, &params);

	switch (otherIn.order)
	{
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
		case 8:
			params.order = LAL_PNORDER_PSEUDO_FOUR;
			break;
		default:
			fprintf(stderr, "unimplemented post-Newtonian order: %d/2\n", otherIn.order);
			exit(1);
	}

	if (strcmp(otherIn.waveformString, "TaylorEt") == 0)
		params.approximant = TaylorEt;
	else if (strcmp(otherIn.waveformString, "TaylorT1") == 0)
		params.approximant = TaylorT1;
	else if (strcmp(otherIn.waveformString, "TaylorT2") == 0)
		params.approximant = TaylorT2;
	else if (strcmp(otherIn.waveformString, "TaylorT3") == 0)
		params.approximant = TaylorT3;
	else if (strcmp(otherIn.waveformString, "TaylorT4") == 0)
		params.approximant = TaylorT4;
	else if (strcmp(otherIn.waveformString, "IMRPhenomA") == 0)
		params.approximant = IMRPhenomA;
	else if (strcmp(otherIn.waveformString, "IMRPhenomB") == 0)
		params.approximant = IMRPhenomB;
	else if (strcmp(otherIn.waveformString, "SpinTaylorFrameless") == 0)
		params.approximant = SpinTaylorFrameless;
	else if (strcmp(otherIn.waveformString, "PhenSpinTaylorRD") == 0)
		params.approximant = PhenSpinTaylorRD;
	else if (strcmp(otherIn.waveformString, "EOBNRv2") == 0)
		params.approximant = EOBNRv2;
	else if (strcmp(otherIn.waveformString, "EOBNRv2HM") == 0)
		params.approximant = EOBNRv2HM;
	else
	{
		fprintf(stderr, "unrecognized or unimplemented approximant: %s\n", otherIn.waveformString);
		exit(1);
	}

	params.totalMass = params.mass1 + params.mass2;
        params.distance *= LAL_PC_SI*1e6;
	params.signalAmplitude = 4.0 * params.mu * LAL_MRSUN_SI/params.distance;
	dt = 1. / params.tSampling;

	XLALInspiralInit(&params, &paramsInit);

	start = clock();
	/* --- now we can call the lalsimulation function --- */
	length = XLALSimInspiralChooseTDWaveform(&hplus, &hcross, 0., dt, params.mass1*LAL_MSUN_SI, params.mass2*LAL_MSUN_SI, params.spin1[0], params.spin1[1], params.spin1[2], params.spin2[0], params.spin2[1], params.spin2[2], params.fLower, fRef, params.distance, params.inclination, lambda1, lambda2, waveFlags, nonGRparams, 0, otherIn.order, params.approximant);
	XLALSimInspiralDestroyWaveformFlags(waveFlags);
	XLALSimInspiralDestroyTestGRParam(nonGRparams);
	diff = clock() - start;
	msec = diff * 1000 / CLOCKS_PER_SEC;
	printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
	if ( length == XLAL_FAILURE )
	{
		fprintf( stderr, "error generating lalsimulation waveform\n" );
		exit( 1 );
	}
	length = hplus->data->length;

	fprintf(stderr, "length = %i\n", length);
	fprintf(stderr, "T = %f\n", (float) length * dt);

	name = "wave1.dat";
	outputfile = fopen(name,"w");

	for(i = 0; i < length; i++) {
		fprintf(outputfile,"%e\t%e\t%e\n",
			i*dt,
			hplus->data->data[i],
			hcross->data->data[i]);
	}

	fclose(outputfile);
	fprintf(stdout,"waveform saved in %s\n", name);

	start = clock();
	h_plus = XLALCreateREAL4TimeSeries("", &tc, 0.0, dt, &lalDimensionlessUnit, length * 1.1);
	h_cross = XLALCreateREAL4TimeSeries("", &tc, 0.0, dt, &lalDimensionlessUnit, length * 1.1);
	/* --- now we can call the lalinspiral function --- */
	LALInspiralWaveTemplates(&mystatus, h_plus->data, h_cross->data, &params);
	diff = clock() - start;
	msec = diff * 1000 / CLOCKS_PER_SEC;
	printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);

	if ( mystatus.statusCode )
	{
		fprintf( stderr, "error generating lalinspiral waveform\n" );
		exit( 1 );
	}

	/* --- and finally save in a file --- */

	name = "wave2.dat";
	outputfile = fopen(name,"w");

	length = h_plus->data->length;

	for(i = 0; i < length; i++) {
		fprintf(outputfile,"%e\t%e\t%e\n",
			i*dt,
			h_plus->data->data[i],
			h_cross->data->data[i]);
	}

	fclose(outputfile);
	fprintf(stdout,"waveform saved in %s\n", name);
	return 0;
}
