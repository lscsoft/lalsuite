/*
*  Copyright (C) 2007 Philip Charlton, Duncan Brown, Jolien Creighton
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

#include <config.h>
#if defined LAL_FFTW3_ENABLED
/* fftw3 not yet supported */
int main( void ) { return 77; }
#else /* fftw2 implementation */


#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
#include <fcntl.h>

#include <lal/LALFCTInterface.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>
#include <lal/Random.h>
#include <lal/LALNoiseModels.h>
#include <lal/LALInspiralBank.h>

NRCSID(LALFCTINTERFACETESTC, "$Id:");

#define USAGE "Usage: LALFCTInterfaceTest [ --full-test ] [ --verbose ]\n"

int lalDebugLevel = 1;

/*
#define USE_WRITE
*/

/* Global test parameters */
INT4 doingFullTest = 0;
INT4 doingVerboseTest = 0;

/* Global parameters for chirp functions */

const REAL8 f0 = 175.0L; /* Hz, for LIGO I */

/* Lower frequency cutoff */
const REAL8 fcutoff = 40.0L;  /* Hz, for LIGO I */

/* Sampling rate */
#define SRATE 2048.0L

const REAL8 srate = SRATE; /* Hz */

/* Nyquist frequency */
const REAL8 Nyq = SRATE/2.0L;  /* Hz */

#undef SRATE

/* Initial scaling factors - must be reset before creating the plan */
REAL8 newtPhiScale = 0.0;
REAL8 onePNPhiScale = 0.0;

/* Replacement for printf that only prints when verbose is on */
void
testprintf(const char* format, ...);

void
testprintf(const char* format, ...)
{
    if (doingVerboseTest)
    {
        va_list args;
        va_start(args, format);
        vprintf(format, args);
        va_end(args);
    }
}


/* Quadratic chirp phase function */
static
REAL4 quadPhi(const REAL4 x)
{
    return x*x;
}

/* Unit chirp amplitude */
static
REAL4 unitA(REAL4 x)
{
    x = 0.0;

    return 1.0;
}

/*
  Newtonian chirp phase, scaled so that it's maximum is 1

  The usual expression for the Newtonian part of the phase is

      6                      / f \ (-5/3)
      - * tau_0 * pi * f_0 * | - |
      5                      \f_0/

  however the FCT has a factor of 2*pi built in, so we
  write it as

                       3         / f \ (-5/3)
      2 * pi * tau_0 * - * f_0 * | - |
                       5         \f_0/

  and identify the phase function with

               3         / f \ (-5/3)
      phi(f) = - * f_0 * | - |
               5         \f_0/

  Then we pick a scaling of tau_0 and phi so that the maximum value
  we want to search for tau_0 corresponds to the integer N/2 where
  N is the length of the FCT in the frequency direction:

      tauhat_0 = ? * tau_0

      phihat_0 = ? * phi
*/
static
REAL4 newtPhi(const REAL4 x)
{
    const REAL8 pwr = -5.0/3.0;

    /* Convert x from fraction of Nyquist to Hz */
    const REAL8 f = x*Nyq;

    if (f <= fcutoff)
    {
        return 0.0;
    }
    else
    {
        return (3.0/5.0*f0*pow(f/f0, pwr)/newtPhiScale);
    }
}

/*
  1-PN chirp phase, scaled so that it's maximum is 1

  The usual expression for the 1-PN part of the phase is


      tau_2 * 2 * (pi * f_0)^2 * f^(-1)


  however the FCT has a factor of 2*pi built in, so we
  write it as

                                  / f \ (-1)
      2 * pi * tau_2 * pi * f_0 * | - |
                                  \f_0/

  and identify the phase function with

                          / f \ (-1)
      phi(f) = pi * f_0 * | - |
                          \f_0/

  or for simplicity

      phi(f) = pi * f_0^2 / f

  Then we pick a scaling of tau_0 and phi so that the maximum value
  we want to search for tau_0 corresponds to the integer N/2 where
  N is the length of the FCT in the frequency direction:

      tauhat_1 = ? * tau_2

      phihat_1 = ? * phi
*/
static
REAL4 onePNPhi(const REAL4 x)
{
    /* Convert x from fraction of Nyquist to Hz */
    const REAL8 f = x*Nyq;

    if (f <= fcutoff)
    {
        return 0.0;
    }
    else
    {
        return LAL_PI*f0*f0/f/onePNPhiScale;
    }
}

/*
  Inspiral chirp amplitude

  At this stage we neglect the overall scaling and just use the
  frequency-dependent part of the amplitude, f^(-7/6).

  Again, the amplitude is cutoff at 40 Hz. We also don't need to
  worry about f0 because that is only needed in converting the
  coefficients of the phase function into chirp times (tau's)
*/
static
REAL4 inspiralA(const REAL4 x)
{
    const REAL8 pwr = -7.0/6.0;

    /* Convert x from fraction of Nyquist to Hz */
    const REAL8 f = x*Nyq;

    if (f <= fcutoff)
    {
        return 0.0;
    }
    else
    {
        return pow(f, pwr);
    }
}

/*
  Set these to point to the phase function and amplitude that you want
*/
/* LALFCTPhaseFn phi0 = newtPhi; */
#define phi0 newtPhi
REAL8 phi0Scale = 0.0;

/* LALFCTPhaseFn phi1 = onePNPhi; */
#define phi1 onePNPhi
REAL8 phi1Scale = 0.0;

LALFCTPhaseFn A = unitA;

/*
  The number of dimensions of the FCT: 1 for the linear term plus 1 for
  each phase function
*/
#define NUMBER_OF_DIMENSIONS 3

/* The length (number of samples) of the initial time series */
#define TSERIES_LENGTH 65536

/*
  The number of data cubes to do in one FCT calculation. We only
  usually do one data cube at a time
*/
#define NUM_DATA_CUBES 1

/*
  This is the number of times a new time series is obtained or generated
*/
const INT4 numDataSegments = 1;

/* const INT4 tseries_length = TSERIES_LENGTH; */
enum { tseries_length = TSERIES_LENGTH };

/* Variables for generating the list of chirp parameters */

/* A pointer to the list of chirp parameters */
InspiralTemplateList* templateList = 0;

/* Actual size of the list gets returned in this variable */
INT4 NtemplateList = 0;

/*
  This is the length of the data cube in the dim0 direction.
*/
/* const INT4 data_length = TSERIES_LENGTH/2; */
enum { data_length = TSERIES_LENGTH/2 };

/*
  This is the length of the data cube in the dim1 direction
  The length can be 1: this corresponds to choosing a single
  value of tau0.
*/
/* const INT4 dim1_data_length = 1; */
/* const INT4 max_dim1_data_length = TSERIES_LENGTH/2; */
enum { dim1_data_length = 1 };
enum { max_dim1_data_length = TSERIES_LENGTH/2 };

/*
  This is the length of the data cube in the dim2 direction
  The length can be 1: this corresponds to choosing a single
  value of tau2.
*/
/* const INT4 dim2_data_length = 1; */
/* const INT4 max_dim2_data_length = TSERIES_LENGTH/2; */
enum { dim2_data_length = 1 };
enum { max_dim2_data_length = TSERIES_LENGTH/2 };

/* Parameters for setting up the fct plan */
/* const INT4 number_of_dimensions = NUMBER_OF_DIMENSIONS; */
/* const INT4 dimension_0_stride = 1; */
/* const REAL4 offset = 0.0; */
enum { number_of_dimensions = NUMBER_OF_DIMENSIONS };
enum { dimension_0_stride = 1 };
#define offset 0.0

/* delta = 1/N where N is the length of the FCT input */
const REAL4 delta = 2.0/TSERIES_LENGTH;

#undef TSERIES_LENGTH

LALFCTPlan* fctPlan = 0;
LALFCTDataCube cube[NUM_DATA_CUBES];
LALFCTSetDataCubesInput dataCubesIn;

ComplexFFTPlan* plan = 0;

COMPLEX8Vector* gfct_in = 0;
COMPLEX8Vector* gchirp = 0;
REAL4Vector* gtseries = 0;

COMPLEX8Vector* gtmp_in = 0;
COMPLEX8Vector* gtmp_out = 0;

COMPLEX8Vector* gfct_out = 0;

/*
  Get (or generate) a time series
*/
static
void getTimeSeriesData(LALStatus* const status, REAL4Vector* const ts)
{
#ifdef USE_WRITE
    const int fd_output = creat("tseries.out", 0644);
#endif
    RandomParams* randpar = 0;
    UINT4 i = 0;

    for (i = 0; i < ts->length; ++i)
    {
	ts->data[i] = 0.0;
    }

#if 1
    LALCreateRandomParams (status, &randpar, 0);
    LALNormalDeviates(status, ts, randpar);
    LALDestroyRandomParams (status, &randpar);
#endif

#ifdef USE_WRITE
    write(fd_output, ts->data, ts->length*sizeof(*(ts->data)));
#endif
}

/*
  The chirp function h_{tc,t0,t1}(f)

  Parameters are

  f - the value at which to evaluate the chirp function
  tc - coefficient of the linear part of the phase ie. f
  tau0 - coefficient of the first phase function ie. phi0(f)
  tau2 - coefficient of the second phase function ie. phi1(f)

*/
static
void chirpFn(const REAL8 f, const REAL8 tc,
	     const REAL8 tau0, const REAL8 tau2,
	     COMPLEX8* const val)
{
    const REAL8 t0 = tau0*phi0Scale;
    const REAL8 t1 = tau2*phi1Scale;
    const REAL8 phase_f = LAL_TWOPI*(tc*f + t0*phi0(f) + t1*phi1(f));

    val->re =  A(f)*cos(phase_f);
    val->im = -A(f)*sin(phase_f);
}

static
void
insertChirp(COMPLEX8Vector* in)
{
    /* Typical ranges for tau0 are 0.295 to 43.50 */
    const REAL8 tau0 = 1.0;

    /* Typical ranges for tau2 are 0.097 to 1.940 */
    const REAL8 tau2 = 0.2;

    const REAL8 tc = 20.0;
    /*
      t0 = tau0*phi0Scale;
      t2 = tau2*phi1Scale;
    */

    UINT8 l = 0;

    for (l = 0; l < gchirp->length; ++l)
    {
	const REAL8 f = offset + delta*l;
	COMPLEX8 val = { 0.0, 0.0 };

	chirpFn(f, tc, tau0, tau2, &val);

	gchirp->data[l].re = val.re;
	gchirp->data[l].im = val.im;
    }

    /* Add it - why is there no function for adding vectors in LAL??? */
    for (l = 0; l < in->length; ++l)
    {
	const REAL8 f = (offset + delta*l)*Nyq;

	in->data[l].re += gchirp->data[l].re;
 	in->data[l].im += gchirp->data[l].im;

	/* Zero out everything below the cutoff */
	if (f <= fcutoff)
	{
	    in->data[l].re = 0.0;
	    in->data[l].im = 0.0;
	}
    }

}

static
void
generateFCTInput(LALStatus* const status,
		 COMPLEX8Vector* out,
		 REAL4Vector* in)
{
#ifdef USE_WRITE
    const int fd_output = creat("fseries.out", 0644);
#endif

    UINT4 i = 0;

    for (i = 0; i < gtmp_in->length; ++i)
    {
	gtmp_in->data[i].re = in->data[i];
	gtmp_in->data[i].im = 0.0;
    }

    LALCOMPLEX8VectorFFT(status, gtmp_out, gtmp_in, plan);

    /* Only fill in the output vector up to the size requested in out */
    for (i = 0; i < out->length; ++i)
    {
	out->data[i].re = gtmp_out->data[i].re;
	out->data[i].im = gtmp_out->data[i].im;
    }

    /* Now insert the phase function into the background signal */
    insertChirp(out);

#ifdef USE_WRITE
    write(fd_output, out->data, out->length*sizeof(*(out->data)));
#endif
}

static
void
calculateFCTCubes(LALStatus* status, COMPLEX8Vector* fct_in)
{
    REAL8 overall_max = 0.0;

    INT4 i = 0;
    INT4 k_max = 0;
    INT4 k0_max = 0;
    INT4 k1_max = 0;
    INT4 Ntemplates = 0;

    REAL4 tau0_max = 0.0;
    REAL4 tau2_max = 0.0;

    /* Strides are always the same */
    cube[0].stride[0] = 1;
    cube[0].stride[1] = 1;
    cube[0].stride[2] = 1;

    /*
      Data cube dimensions for the frequency axis are always going
      to be the full length
    */
    cube[0].start_locations[0] = 0;
    cube[0].end_locations[0] = data_length;

    /* Initial dim1 data cube limits */
    cube[0].start_locations[1] = 0;
    cube[0].end_locations[1] = 0;

    dataCubesIn.data_cube = &cube[0];
    dataCubesIn.num_data_cubes = NUM_DATA_CUBES;

    /*
      To make the test shorter, set the number of templates to a
      small number unless --full-test is on
    */

    if (doingFullTest)
    {
	Ntemplates = NtemplateList;
    }
    else
    {
	Ntemplates = 4;
    }

    for (i = 0; i < Ntemplates; ++i)
    {
	const REAL4 tc_start = cube[0].start_locations[0];
	const REAL4 tc_end   = cube[0].end_locations[0] - 1;

	REAL8 pwr = 0.0;
	REAL8 max_pwr = 0.0;

	/*
	  The offset into the output where the max power for this
	  data cube occurs
	*/
	INT4 max_pwr_index = 0;

	/* Indices into the current data cube for the chirp parameters */
	INT4 K  = 0;
	INT4 K0 = 0;
	INT4 K1 = 0;

	/* Indices into the whole space */
	INT4 k = 0;
	INT4 k0 = 0;
	INT4 k1 = 0;

	UINT4 l = 0;
	INT4 m = 0;

	REAL4 tau0 = 0.0;
	REAL4 tau2 = 0.0;

	REAL4 tau0_start = 0.0;
	REAL4 tau0_end   = 0.0;

	REAL4 tau2_start = 0.0;
	REAL4 tau2_end   = 0.0;

	/*
	  Convert the tau's from the template list to the appropriate
	  integer index
	*/
	cube[0].start_locations[1]
	    = (templateList[i].params.t0*phi0Scale + 0.5);
	cube[0].end_locations[1] = cube[0].start_locations[1] + 1;

	tau0_start = cube[0].start_locations[1]/phi0Scale;
	tau0_end   = (cube[0].end_locations[1] - 1)/phi0Scale;

	/* Tau2 */
	cube[0].start_locations[2]
	    = (templateList[i].params.t2*phi1Scale + 0.5);
	cube[0].end_locations[2] = cube[0].start_locations[2] + 1;

	tau2_start = cube[0].start_locations[2]/phi1Scale;
	tau2_end   = (cube[0].end_locations[2] - 1)/phi1Scale;

	testprintf("Searching over data cube bounded by:\n");
	testprintf("    %5d <= k  <= %5d (%5f <= tc   <= %5f)\n",
	       cube[0].start_locations[0],
	       cube[0].end_locations[0] - 1,
	       tc_start,
	       tc_end);
	testprintf("    %5d <= k0 <= %5d (%5f <= tau0 <= %5f)\n",
	       cube[0].start_locations[1],
	       cube[0].end_locations[1] - 1,
	       tau0_start,
	       tau0_end);
	testprintf("    %5d <= k1 <= %5d (%5f <= tau2 <= %5f)\n",
	       cube[0].start_locations[2],
	       cube[0].end_locations[2] - 1,
	       tau2_start,
	       tau2_end);

	LALFCTSetDataCubes(status, &dataCubesIn, fctPlan);

	LALFCTCalculate(status, gfct_out, fct_in, fctPlan);

	for (l = 0; l < gfct_out->length; ++l)
	{
	    pwr = gfct_out->data[l].re*gfct_out->data[l].re
		+ gfct_out->data[l].im*gfct_out->data[l].im;
	    if (pwr > max_pwr)
	    {
		max_pwr = pwr;
		max_pwr_index = l;
	    }
	}

	/*
	  Using the formula

	  l = K0*data_length + K

	  we can solve for K and K0
	*/

	K = max_pwr_index % data_length;
	m = (max_pwr_index - K)/data_length;
	K0 = m % dim1_data_length;
	K1 = (m - K0)/dim1_data_length;

	k  = cube[0].start_locations[0] + K;
	k0 = cube[0].start_locations[1] + K0;
	k1 = cube[0].start_locations[2] + K1;

	tau0 = k0/phi0Scale;
	tau2 = k1/phi1Scale;

	testprintf("Max pwr: fct_out->data[%d](K = %d, K0 = %d, K1 = %d) = %g\n",
	       max_pwr_index, K, K0, K1, max_pwr);
	testprintf("    k  = %5d (tc   = %d)\n", k, k);
	testprintf("    k0 = %5d (tau0 = %f)\n", k0, tau0);
	testprintf("    k1 = %5d (tau2 = %f)\n", k1, tau2);
	testprintf("\n");

	if (max_pwr > overall_max)
	{
	    overall_max = max_pwr;
	    k_max  = k;
	    k0_max = k0;
	    k1_max = k1;
	}

    }

    tau0_max = k0_max/phi0Scale;
    tau2_max = k1_max/phi1Scale;

    testprintf("Overall max pwr: %g\n", overall_max);
    testprintf("    k  = %5d (tc   = %d)\n", k_max, k_max);
    testprintf("    k0 = %5d (tau0 = %f)\n", k0_max, tau0_max);
    testprintf("    k1 = %5d (tau2 = %f)\n", k1_max, tau2_max);

#ifdef USE_WRITE
    writefct((float*) fct_out->data,
	     data_length, fct_out->length/data_length, "fct.out");
#endif
}

/*
  This is the old version that I didn't want to change - it
  runs over some range of parameters
*/
#if 0
static
void
calculateFCTCubes2(LALStatus* status, COMPLEX8Vector* fct_in)
{
    REAL8 overall_max = 0.0;

    UINT4 l = 0;
    INT4 k_max = 0;
    INT4 k0_max = 0;

    REAL4 tau0_max = 0.0;

    /* Strides are always the same */
    cube[0].stride[0] = 1;
    cube[0].stride[1] = 1;

    /*
      Data cube dimensions for the frequency axis are always going
      to be the full length
    */
    cube[0].start_locations[0] = 0;
    cube[0].end_locations[0] = data_length;

    /* Initial dim1 data cube limits */
    cube[0].start_locations[1] = 0;
    cube[0].end_locations[1] = dim1_data_length;

    dataCubesIn.data_cube = &cube[0];
    dataCubesIn.num_data_cubes = NUM_DATA_CUBES;

    while (cube[0].start_locations[1] < max_dim1_data_length)
    {
	const REAL4 tc_start = cube[0].start_locations[0];
	const REAL4 tc_end   = cube[0].end_locations[0] - 1;

	const REAL4 tau0_start = cube[0].start_locations[1]/phi0Scale;
	const REAL4 tau0_end   = (cube[0].end_locations[1] - 1)/phi0Scale;

	REAL8 pwr = 0.0;
	REAL8 max_pwr = 0.0;

	/*
	  The offset into the output where the max power for this
	  data cube occurs
	*/
	INT4 max_pwr_index = 0;

	/* Indices into the current data cube for the chirp parameters */
	INT4 K  = 0;
	INT4 K0 = 0;

	/* Indices into the whole space */
	INT4 k = 0;
	INT4 k0 = 0;

	REAL4 tau0 = 0.0;

	testprintf("Searching over data cube bounded by:\n");
	testprintf("    %5d <= k  <= %5d (%5f <= tc   <= %5f)\n",
	       cube[0].start_locations[0],
	       cube[0].end_locations[0] - 1,
	       tc_start,
	       tc_end);
	testprintf("    %5d <= k0 <= %5d (%5f <= tau0 <= %5f)\n",
	       cube[0].start_locations[1],
	       cube[0].end_locations[1] - 1,
	       tau0_start,
	       tau0_end);

	LALFCTSetDataCubes(status, &dataCubesIn, fctPlan);

	LALFCTCalculate(status, gfct_out, fct_in, fctPlan);

	for (l = 0; l < gfct_out->length; ++l)
	{
	    pwr = gfct_out->data[l].re*gfct_out->data[l].re
		+ gfct_out->data[l].im*gfct_out->data[l].im;
	    if (pwr > max_pwr)
	    {
		max_pwr = pwr;
		max_pwr_index = l;
	    }
	}

	/*
	  Using the formula

	  l = K0*data_length + K

	  we can solve for K and K0
	*/
	K0 = max_pwr_index/data_length;
	K  = max_pwr_index - K0*data_length;

	k  = cube[0].start_locations[0] + K;
	k0 = cube[0].start_locations[1] + K0;

	tau0 = k0/phi0Scale;

	testprintf("Max pwr: fct_out->data[%d](K = %d, K0 = %d) = %g\n",
	       max_pwr_index, K, K0, max_pwr);
	testprintf("    k  = %5d (tc   = %d)\n", k, k);
	testprintf("    k0 = %5d (tau0 = %f)\n", k0, tau0);
	testprintf("\n");

	if (max_pwr > overall_max)
	{
	    overall_max = max_pwr;
	    k_max  = k;
	    k0_max = k0;
	}

	/* Go to the next data cube */
	cube[0].start_locations[1] = cube[0].end_locations[1];
	cube[0].end_locations[1]
	    = cube[0].start_locations[1] + dim1_data_length;
    }

    tau0_max = k0_max/phi0Scale;

    testprintf("Overall max pwr: %g\n", overall_max);
    testprintf("    k  = %5d (tc   = %d)\n", k_max, k_max);
    testprintf("    k0 = %5d (tau0 = %f)\n", k0_max, tau0_max);

#ifdef USE_WRITE
    writefct((float*) fct_out->data,
	     data_length, fct_out->length/data_length, "fct.out");
#endif
}
#endif /* exclusion of function */

static
void
InspiralSearch(LALStatus* const status)
{
    INT4 curr_data_segment = 0;

    testprintf("Commencing search...\n");

    while (curr_data_segment < numDataSegments)
    {
	getTimeSeriesData(status, gtseries);
	generateFCTInput(status, gfct_in, gtseries);
	calculateFCTCubes(status, gfct_in);
	++curr_data_segment;
    }
}

/*
  Use the LALInspiralBank code to generate a list of tau's
*/
static
void
setupParamList(LALStatus* const status)
{
    INT4 i = 0;

    const UINT4 numPSDpts = 1048576;
    void *noisemodel = LALLIGOIPsd;

    REAL8 t0_min = LAL_REAL8_MAX;
    REAL8 t0_max = 0.0;
    REAL8 t2_min = LAL_REAL8_MAX;
    REAL8 t2_max = 0.0;

    InspiralCoarseBankIn coarseIn;

    testprintf("Generating list of chirp parameters...\n");

    coarseIn.massRange = MinMaxComponentMass;
    coarseIn.mMin = 2.0;
    coarseIn.mMax = 3.0;
    coarseIn.MMax = 2.0*coarseIn.mMax; /* 10.6 */
    coarseIn.mmCoarse = 0.80;
    coarseIn.mmFine = 0.97;
    coarseIn.fLower = fcutoff;
    coarseIn.fUpper = Nyq;
    coarseIn.tSampling = srate;
    coarseIn.order = twoPN;
    coarseIn.approximant = TaylorT1;
    coarseIn.space = Tau0Tau3;

    /* minimum value of eta */
    coarseIn.etamin = coarseIn.mMin*(coarseIn.MMax - coarseIn.mMin)/
	pow(coarseIn.MMax,2.0);
    coarseIn.iflso = 0;

    /* the bank creation routines now take a frequency series containing */
    /* the psd */
    memset( &(coarseIn.shf), 0, sizeof(REAL8FrequencySeries) );
    coarseIn.shf.f0 = 0;
    LALDCreateVector( status, &(coarseIn.shf.data), numPSDpts );
    coarseIn.shf.deltaF = coarseIn.tSampling / (REAL8) coarseIn.shf.data->length;
    LALNoiseSpectralDensity (status, coarseIn.shf.data, noisemodel, coarseIn.shf.deltaF );


    /* LALInspiralCreateCoarseBank() now expects its second argument to */
    /* be a handle (pointer to a pointer) to a variable of type         */
    /* InspiralTemplateList. The pointer should be NULL on entry        */
    /*                                  - duncan                        */
    LALInspiralCreateCoarseBank(status, &templateList,
				&NtemplateList, coarseIn);

    LALDDestroyVector( status, &(coarseIn.shf.data) );

    testprintf("Generated %d sets of parameters\n", NtemplateList);

    testprintf("tau0 tau2 tau3 M eta m1 m2\n");
    for (i = 0; i < NtemplateList; ++i)
    {
	testprintf("%e %e %e %e %e %e %e\n",
	       templateList[i].params.t0,
	       templateList[i].params.t2,
	       templateList[i].params.t3,
	       templateList[i].params.totalMass,
	       templateList[i].params.eta,
	       templateList[i].params.mass1,
	       templateList[i].params.mass2);

	if (t0_max < templateList[i].params.t0)
	{
	    t0_max = templateList[i].params.t0;
	}

	if (t0_min > templateList[i].params.t0)
	{
	    t0_min = templateList[i].params.t0;
	}

	if (t2_max < templateList[i].params.t2)
	{
	    t2_max = templateList[i].params.t2;
	}

	if (t2_min > templateList[i].params.t2)
	{
	    t2_min = templateList[i].params.t2;
	}
    }

    testprintf("Min tau0 = %g\n", t0_min);
    testprintf("Max tau0 = %g\n", t0_max);

    testprintf("Min tau2 = %g\n", t2_min);
    testprintf("Max tau2 = %g\n", t2_max);
}

static
void
setupGlobals(LALStatus* const status, int argc, char **argv)
{
    const LALCreateFCTPlanInput planIn = {
	data_length,
	number_of_dimensions,
	dimension_0_stride,
	/* Array of phase fns */
	{
	    phi0,
	    phi1
	}
    };

    const LALFCTSetUnitsInput setUnitsIn = {
	offset,
	/* delta */
	2.0/tseries_length
    };

    LALFCTSetOversamplingFactorInput setOfac = { 0, 0 };

    UINT8 output_data_size = 0;

    /* I have to guess what the size of the list of chirp params will be */
    /* const INT4 Nlist = 5000; */

    INT4 i = 0;

    /* Do command line arguments */
    for (i = 1; i < argc; ++i)
    {
	if (strcmp("-h", argv[i]) == 0)
	{
	    printf(USAGE);
	    exit(0);
	}
	else if (strcmp("--full-test", argv[i]) == 0)
	{
	    doingFullTest = 1;
	}
	else if (strcmp("--verbose", argv[i]) == 0)
	{
	    doingVerboseTest = 1;
	}
	else
	{
	    printf("Invalid argument: %s\n", argv[i]);
	    printf(USAGE);
	    exit(1);
	}
    }

    /*
      Set the global scaling functions.

      The relationship is that

        t_i   = tau_i * scale factor
        phi_i = phase_fn_i / scale factor
    */
    newtPhiScale  = 3.0/5.0*f0*pow(fcutoff/f0, -5.0/3.0);
    onePNPhiScale = LAL_PI*f0*f0/fcutoff;

    /* Set phase function scaling factors */
    phi0Scale = newtPhiScale;
    phi1Scale = onePNPhiScale;

    /* Set up list of templates */
    /* this memory allocation is now done with a realloc() in         */
    /* LALCreateCoarseInspiralBank()                                  */
    /*                                          - duncan              */
    /* templateList = LALCalloc(Nlist, sizeof(InspiralTemplateList)); */
    setupParamList(status);

    /* Create the various arrays we use for storage */
    LALSCreateVector(status, &gtseries, tseries_length);

    LALCCreateVector(status, &gfct_in, data_length);
    LALCCreateVector(status, &gtmp_in, data_length);
    LALCCreateVector(status, &gtmp_out, data_length);
    LALCCreateVector(status, &gchirp, data_length);

    /*
      Plan for doing the FFT of the tseries. Note that this
      is a "backward" plan in LAL convention - LAL and FCT
      use the convention that a "forward" transform has a
      factor of +1 in the exponent and a "backward" transform
      has a -1.
    */
    /* LALEstimateInvComplexFFTPlan(status, &plan, gtmp_in->length); */
    /*
     * >> THIS HAS BEEN CHANGED <<
     *
     * Now "forward" has a -1 in LAL, so I am changing the above line to do a
     * transform with a -1 (as before), which is now the "forward" transform
     *
     * --Jolien
     */
    LALCreateForwardComplexFFTPlan(status, &plan, gtmp_in->length, 0);

    /* Set up the plan and output area */
    LALCreateFCTPlan(status, &fctPlan, &planIn);

    setOfac.ofac = 1;
    setOfac.dim = 1;
    LALFCTSetOversamplingFactor(status, &setOfac, fctPlan);

    setOfac.ofac = 1;
    setOfac.dim = 2;
    LALFCTSetOversamplingFactor(status, &setOfac, fctPlan);

    LALFCTSetUnits(status, &setUnitsIn, fctPlan);

    /*
      Initialise a data cube of the size we plan to use, so
      we can allocate the correct output size
    */
    cube[0].start_locations[0] = 0;
    cube[0].end_locations[0]   = data_length;

    cube[0].start_locations[1] = 0;
    cube[0].end_locations[1]   = dim1_data_length;

    cube[0].start_locations[2] = 0;
    cube[0].end_locations[2]   = dim2_data_length;

    cube[0].stride[0] = 1;
    cube[0].stride[1] = 1;
    cube[0].stride[2] = 1;

    dataCubesIn.data_cube = &cube[0];
    dataCubesIn.num_data_cubes = NUM_DATA_CUBES;

    /* Get the required output size */
    LALFCTSetDataCubes(status, &dataCubesIn, fctPlan);
    LALFCTOutputDataSize(status, &output_data_size, fctPlan);

    /* Allocate an output vector of this size */
    LALCCreateVector(status, &gfct_out, output_data_size);
}

static
void
cleanup(LALStatus* const status)
{
    /* Deallocate the arrays and plan */

    LALCDestroyVector(status, &gfct_out);

    LALDestroyFCTPlan(status, &fctPlan);

    LALDestroyComplexFFTPlan(status, &plan);

    LALCDestroyVector(status, &gchirp);
    LALCDestroyVector(status, &gtmp_out);
    LALCDestroyVector(status, &gtmp_in);
    LALCDestroyVector(status, &gfct_in);

    LALSDestroyVector(status, &gtseries);

    LALFree(templateList);
}

int
main(int argc, char **argv)
{
    LALStatus* status = LALCalloc(1, sizeof(*status));

    setupGlobals(status, argc, argv);

    InspiralSearch(status);

    cleanup(status);

    LALFree(status);

    LALCheckMemoryLeaks();

    return 0;
}
#endif /* fftw2 implementation */
