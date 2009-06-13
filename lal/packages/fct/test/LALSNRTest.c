/*
*  Copyright (C) 2007 Philip Charlton, Jolien Creighton
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
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <assert.h>

#include <lal/LALConstants.h>
#include <lal/LALNoiseModels.h>
#include <lal/LALFCTInterface.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/ComplexFFT.h>
#include <lal/GeneratePPNInspiral.h>

NRCSID(LALSNRTESTC, "$Id:");

#undef USE_WRITE

int lalDebugLevel = 1;

/* Global parameters for chirp functions */
/* Lower frequency cutoff */
/* NOTE - now using Sathya's definition of tau_0, not Owen's */
const REAL8 f_seismic = 40.0;

/* Sampling rate */
#define SRATE 2048
#define NSECS 64
#define TSERIES_LENGTH (SRATE*NSECS)
#define NUMBER_OF_DIMENSIONS 2
#define NUM_DATA_CUBES 1

/* Nyquist frequency */
const REAL8 NYQUIST = SRATE/2.0;  /* Hz */

/*
  This is the length of the data cube in the dim0 direction.
*/
enum { data_length = TSERIES_LENGTH/2 };

/*
  This is the length of the data cube in the dim1 direction
  The length can be 1: this corresponds to choosing a single
  value of tau0.
*/
const int dim1_data_length = 8;
const int max_dim1_data_length = TSERIES_LENGTH/2;

/* Parameters for setting up the fct plan */
enum { dimension_0_stride = 1 };

/* OFFSET from 0 and DELTA = Nyquist/data_length */
/*
const REAL8 OFFSET = 0.0;
const REAL8 DELTA  = SRATE/((REAL8) TSERIES_LENGTH);
*/
#define OFFSET 0.0
#define DELTA (SRATE/(REAL8)TSERIES_LENGTH)
#define SDELTA (SRATE/(REAL4)TSERIES_LENGTH)

LALFCTDataCube cube[NUM_DATA_CUBES];
LALFCTSetDataCubesInput dataCubesIn;

LALFCTPlan* gfctPlan = 0;
COMPLEX8Vector* gfct_out = 0;

REAL8Vector* gSh = 0;

#ifdef USE_WRITE
int writefct(const REAL4* const out, const int n, const int m,
             const char* const filename);

int writefct(const REAL4* const out, const int n, const int m,
             const char* const filename)
{
    /* The size of each element in the out array */
    const float out_size = sizeof(*out);

    const REAL4 M = m;
    const REAL4 N = n;

    /* Open the file */
    const int fd_output = creat(filename, 0644);

    if (fd_output != 0)
    {
        /* Write info about endian-ness and local sizeof(*out) */
        write(fd_output, &out_size, sizeof(out_size));

        /* Write the dimensions of the fct */
        write(fd_output, &M, sizeof(M));
        write(fd_output, &N, sizeof(N));

        /* Write the body of the fct */
        write(fd_output, out, sizeof(*out)*2*n*m);

        close(fd_output);
    }
    else
    {
        return -1;
    }

    return 0;
}
#endif

/*
  Return the inner product of two frequency series, *assuming* that they
  are the positive frequency components of the FFT of a real series
*/
static
REAL8 inner(const COMPLEX8Vector* const f, const COMPLEX8Vector* const g)
{
    UINT4 k = 0;
    REAL8 sum = 0.0;

    assert(f->length == g->length);
    assert(gSh != 0);
    assert(gSh->length == f->length);

    /* Want real(sum(f * conj(g))) */
    for (k = 0; k < f->length; ++k)
    {
        sum += (f->data[k].re*g->data[k].re + f->data[k].im*g->data[k].im)/
	    gSh->data[k];
    }

    /*
      This is being subtracted so that we get the equivalent of the
      2-sided inner product
    */
    sum = (2.0*sum - f->data[0].re*g->data[0].re/gSh->data[0])/f->length;

    return sum;
}

#if 0
static
void generateTimeDomainChirp(LALStatus* const status,
                             const REAL8 tau0,
                             REAL4Vector* const out,
                             REAL8* const tc,
                             REAL8* MTot,
                             REAL8* Eta)
{
#ifdef USE_WRITE
    const int fd_output = creat("tchirp.out", 0644);
#endif

    PPNParamStruc params;
    CoherentGW waveform;

    const REAL8 eta = 0.25;

    const REAL8 t1 = pow(LAL_PI*f_seismic, 8.0/3.0);
    const REAL8 mTot = pow(256.0*t1*eta*tau0/5.0, -3.0/5.0)/LAL_MTSUN_SI;

    const REAL8 dist = 8.5;

    UINT4 j = 0;

    assert(out != 0);
    assert(out->length > 0);

    params.position.latitude = 0.0;
    params.position.longitude = 0.0;
    params.position.system = COORDINATESYSTEM_EQUATORIAL;
    params.psi = 0.0;
    params.lengthIn = out->length;

    /* Variable parameters. */
    params.epoch.gpsSeconds = 6100000;
    params.epoch.gpsNanoSeconds = 0;

    params.deltaT = 1.0/SRATE;
    params.mTot = mTot;
    params.eta = eta;
    params.inc = 0;

    /*
      Tev's code uses a coalescence phase different by a factor of 2
      (so I should use -LAL_PI/4 here instead of -LAL_PI/8) BUT also
      has an overall sign of -1 on the amplitude of chirp, giving
      another +LAL_PI, so 3/LAL_PI/4 is the total
    */
    params.phi = 3.0*LAL_PI/4.0;
    params.d = dist*LAL_PC_SI*1.0e+03;
    params.fStartIn = f_seismic;
    params.fStopIn = NYQUIST;

    /* PPN parameter. */
    params.ppn = 0;
    LALSCreateVector(status, &(params.ppn), 1);
    params.ppn->data[0] = 1.0;

    /* Output parameters. */
    memset( &waveform, 0, sizeof(waveform) );

    /* Generate waveform. */
    LALGeneratePPNInspiral(status, &waveform, &params);

    for (j = 0; j < params.length; ++j)
    {
        REAL8 p = waveform.phi->data->data[j];
        REAL8 ap = waveform.a->data->data[2*j];

        out->data[j] = ap*cos(p);
    }

    for (j = params.length; j < out->length; ++j)
    {
        out->data[j] = 0.0;
    }

    LALSDestroyVector(status, &(params.ppn));
    LALSDestroyVectorSequence(status, &(waveform.a->data));
    LALSDestroyVector(status, &(waveform.f->data));
    LALDDestroyVector(status, &(waveform.phi->data));

    *tc = params.tc;
    *MTot = mTot;
    *Eta = eta;

#ifdef USE_WRITE
    write(fd_output, out->data, out->length*sizeof(*(out->data)));
#endif
}
#endif

static
REAL4 maxPhi(LALFCTPhaseFn phi, const INT4 N)
{
    REAL4 max = -LAL_REAL4_MAX;
    INT4 k = 0;

    for (k = 0; k < N; ++k)
    {
        const REAL4 f = OFFSET + DELTA*k;
        if (phi(f) > max)
        {
            max = phi(f);
        }
    }

    return max;
}

/*
  Newtonian chirp phase, scaled so that it's maximum is 1
*/
static
REAL4 newtPhi(const REAL4 f)
{
    const REAL8 pwr = -5.0/3.0;

    if (f < f_seismic)
    {
        return 0.0;
    }
    else
    {
        return pow(f/f_seismic, pwr);
    }
}

/*
  Inspiral chirp amplitude

  At this stage we neglect the overall scaling and just use the
  frequency-dependent part of the amplitude, f^(-7/6).
*/
static
REAL4 inspiralA(const REAL4 f)
{
    const REAL8 pwr = -7.0/6.0;

    if (f < f_seismic)
    {
        return 0.0;
    }
    else
    {
        return pow(f, pwr);
    }
}

/*
  Unit inspiral chirp amplitude
*/
static
REAL4 unitInspiralA(const REAL4 f)
{
    if (f < f_seismic)
    {
        return 0.0;
    }
    else
    {
        return 1.0;
    }
}

LALFCTPhaseFn A = inspiralA;

static
REAL4 phi0(const REAL4 f)
{
    static INT4 initialised = 0;
    static REAL4 max = 0.0;

    if (initialised == 0)
    {
        max = maxPhi(newtPhi, data_length);
        initialised = 1;
    }

    return newtPhi(f)/max;
}

/*
  The chirp function h_{tc,lam0}(f) where f is frequency in Hz

  Parameters are

  tc - coefficient of the linear part of the phase ie. f
  lam0 - coefficient of the first phase function ie. phi0(f)
*/
static
void
generateFrequencyDomainChirp(const REAL8 tc, const REAL8 lam0,
                             COMPLEX8Vector* const in)
{
#ifdef USE_WRITE
    const int fd_output = creat("fchirp.out", 0644);
#endif

    UINT8 l = 0;

    assert(in != 0);
    assert(in->length > 0);

    for (l = 0; l < in->length; ++l)
    {
        const REAL4 f = OFFSET + DELTA*l;
        const REAL4 Amp = A(f);

        REAL8 phase_f = LAL_TWOPI*(tc*f + lam0*phi0(f));

        phase_f = fmod(phase_f, LAL_TWOPI);
        if (phase_f < 0)
        {
            phase_f += LAL_TWOPI;
        }

        /* Overall sign of phase function is -1 */
        in->data[l].re =  Amp*cos(phase_f);
        in->data[l].im = -Amp*sin(phase_f);
    }

#ifdef USE_WRITE
    write(fd_output, in->data, in->length*sizeof(*(in->data)));
#endif
}

/*
  This function generates the template that is equivalent to that used
  by the FCT
*/
static
void
generateFCTTemplate(const INT4 k, const INT4 k0,
		    COMPLEX8Vector* const in)
{
#ifdef USE_WRITE
    const int fd_output = creat("template.out", 0644);
#endif

    const REAL4 M = in->length;
    const REAL4 N = in->length;

    UINT8 n = 0;

    assert(in != 0);
    assert(in->length > 0);

    for (n = 0; n < in->length; ++n)
    {
        const REAL4 f = OFFSET + DELTA*n;
        const REAL4 Amp = A(f);
	const INT8 mu_n = (M*phi0(f) + 0.5);

        REAL8 phase_f = LAL_TWOPI*((k*n)/N + (k0*mu_n)/M);

        phase_f = fmod(phase_f, LAL_TWOPI);
        if (phase_f < 0)
        {
            phase_f += LAL_TWOPI;
        }

        in->data[n].re =  Amp*cos(phase_f);
        in->data[n].im = -Amp*sin(phase_f);
    }

#ifdef USE_WRITE
    write(fd_output, in->data, in->length*sizeof(*(in->data)));
#endif
}

#if 0
static
void
transformTimeDomainChirp(LALStatus* const status,
                         COMPLEX8Vector* const out,
                         const REAL4Vector* const in)
{
    ComplexFFTPlan* plan = 0;

    COMPLEX8Vector* tmp_in = 0;
    COMPLEX8Vector* tmp_out = 0;

    UINT4 i = 0;

    assert(out->length == in->length/2);

    LALCCreateVector(status, &tmp_in, 2*out->length);
    LALCCreateVector(status, &tmp_out, 2*out->length);

    for (i = 0; i < tmp_in->length; ++i)
    {
        tmp_in->data[i].re = in->data[i];
        tmp_in->data[i].im = 0.0;
    }

    LALCreateForwardComplexFFTPlan(status, &plan, tmp_in->length, 0);

    LALCOMPLEX8VectorFFT(status, tmp_out, tmp_in, plan);

    for (i = 0; i < out->length; ++i)
    {
        out->data[i].re = tmp_out->data[i].re;
        out->data[i].im = tmp_out->data[i].im;
    }

    LALDestroyComplexFFTPlan(status, &plan);

    LALCDestroyVector(status, &tmp_out);
    LALCDestroyVector(status, &tmp_in);
}
#endif

static
void
FCTCalculate(LALStatus* status, COMPLEX8Vector* const fct_out,
             const COMPLEX8Vector* const fct_in,
             const LALFCTSetDataCubesInput* cubesIn)
{
    UINT4 l = 0;
    COMPLEX8Vector* tmp_fct_in = 0;

    assert(cubesIn != 0);
    assert(fct_in != 0);
    assert(fct_out != 0);

    assert(fct_in->length == data_length);

    LALCCreateVector(status, &tmp_fct_in, fct_in->length);

    /*
      The input needs to be scaled by the amplitude of the template,
      and by the inverse of the power spectrum, since this is NOT included
      in the FCT
    */

    for (l = 0; l < tmp_fct_in->length; ++l)
    {
        const REAL8 f = OFFSET + DELTA*l;

        tmp_fct_in->data[l].re = A(f)*fct_in->data[l].re/gSh->data[l];
        tmp_fct_in->data[l].im = A(f)*fct_in->data[l].im/gSh->data[l];
    }

    LALFCTSetDataCubes(status, cubesIn, gfctPlan);

    LALFCTCalculate(status, fct_out, tmp_fct_in, gfctPlan);

    LALCDestroyVector(status, &tmp_fct_in);
}

static
void
calculateFCTInner(LALStatus* const status,
                  const COMPLEX8Vector* const fct_in,
                  const INT4 dim1_start,
                  REAL4* const inner_est,
                  INT4* const tc_est,
                  INT4* const lam0_est)
{
    INT4 k = 0;
    INT4 k0 = 0;
    INT4 K = 0;
    INT4 K0 = 0;
    INT4 m = 0;

    UINT4 l = 0;

    REAL8 inner_prod = 0.0;
    REAL8 max_inner_prod = -1.0e+99;
    INT4  max_inner_prod_index = 0;

    assert(fct_in != 0);
    assert(fct_in->length > 0);

    cube[0].start_locations[0] = 0;
    cube[0].end_locations[0]   = data_length;

    cube[0].start_locations[1] = dim1_start;
    cube[0].end_locations[1]   = dim1_start + dim1_data_length;

    assert(gfct_out != 0);
    assert(gfct_out->length % fct_in->length == 0);

    FCTCalculate(status, gfct_out, fct_in, &dataCubesIn);

    for (l = 0; l < gfct_out->length; ++l)
    {
        /*
          pwr = gfct_out->data[l].re*gfct_out->data[l].re
          + gfct_out->data[l].im*gfct_out->data[l].im;
        */
        inner_prod = 2.0*gfct_out->data[l].re/fct_in->length;
        if (inner_prod > max_inner_prod)
        {
            max_inner_prod = inner_prod;
            max_inner_prod_index = l;
        }
    }

    /*
      Using the formula

      l = K0*data_length + K

      we can solve for K and K0
    */

    K = max_inner_prod_index % data_length;
    m = (max_inner_prod_index - K)/data_length;
    K0 = m % dim1_data_length;

    k  = cube[0].start_locations[0] + K;
    k0 = cube[0].start_locations[1] + K0;

    *inner_est = max_inner_prod;
    *tc_est = k;
    *lam0_est = k0;

#ifdef USE_WRITE
    writefct((float*) gfct_out->data,
             data_length, gfct_out->length/data_length, "fct.out");
#endif
}

static
void
calculateDCTInner(LALStatus* const status,
                  const COMPLEX8Vector* const dct_in,
                  const INT4 dim1_start,
                  REAL4* const inner_est,
                  INT4* const tc_est,
                  INT4* const lam0_est)
{
    INT4 k = 0;
    INT4 K = 0;
    INT4 K0 = 0;

    UINT4 l = 0;

    REAL8 inner_prod = 0.0;
    REAL8 max_inner_prod = -1.0e+99;

    COMPLEX8Vector* tmp_in = 0;
    COMPLEX8Vector* tmp_out = 0;
    ComplexFFTPlan* plan = 0;

    assert(dct_in != 0);
    assert(dct_in->length > 0);

    LALCCreateVector(status, &tmp_in, dct_in->length);
    LALCCreateVector(status, &tmp_out, dct_in->length);

    LALCreateReverseComplexFFTPlan(status, &plan, tmp_in->length, 0);

    for (k = 0; k < dim1_data_length; ++k)
    {
        for (l = 0; l < tmp_in->length; ++l)
        {
            const REAL4 f = OFFSET + DELTA*l;

            /* NOTE: must use opposite sign to the generated chirp */
            const REAL8 Amp = A(f);
            REAL8 phase_f = LAL_TWOPI*(dim1_start + k)*phi0(f);
            REAL8 cos_ph = 0.0;
            REAL8 sin_ph = 0.0;

	    phase_f = fmod(phase_f, LAL_TWOPI);
	    if (phase_f < 0)
	    {
		phase_f += LAL_TWOPI;
	    }

            cos_ph = cos(phase_f);
            sin_ph = sin(phase_f);

            tmp_in->data[l].re = Amp*(dct_in->data[l].re*cos_ph
                                   - dct_in->data[l].im*sin_ph)/gSh->data[l];
            tmp_in->data[l].im = Amp*(dct_in->data[l].re*sin_ph
                                   + dct_in->data[l].im*cos_ph)/gSh->data[l];
        }

        LALCOMPLEX8VectorFFT(status, tmp_out, tmp_in, plan);

        for (l = 0; l < tmp_out->length; ++l)
        {
            inner_prod = 2.0*tmp_out->data[l].re/dct_in->length;
            if (inner_prod > max_inner_prod)
            {
                max_inner_prod = inner_prod;

                K = l;
                K0 = dim1_start + k;
            }
        }
    }

    *inner_est = max_inner_prod;
    *tc_est = K;
    *lam0_est = K0;

    LALDestroyComplexFFTPlan(status, &plan);

    LALCDestroyVector(status, &tmp_in);
    LALCDestroyVector(status, &tmp_out);
}

static
void
CalculateSNRLoss(LALStatus* const status)
{
    const REAL4 deltaLam0 = 1000.5;

    const REAL8 tc = 3*NSECS/4.0;
    REAL4 tau0 = 0.0;

    REAL4 lambda = 0.0;
    REAL4 lam0 = 0.0;
    INT4 fct_lam_est = 0;
    INT4 fct_lam0_est = 0;

    INT4 dct_lam_est = 0;
    INT4 dct_lam0_est = 0;

    COMPLEX8Vector* chirp = 0;
    COMPLEX8Vector* template = 0;

    REAL4 inner_fct = 0.0;
    REAL4 inner_dct = 0.0;
    REAL4 tmp_inner = 0.0;

    REAL4 SNR_fct = 0.0;
    REAL4 SNR_dct = 0.0;

    REAL4 chirp_norm = 0.0;
    REAL4 template_norm = 0.0;

    INT4 dim1_start = 0;
    INT4 lam0Idx = 0;

    LALCCreateVector(status, &chirp, data_length);
    LALCCreateVector(status, &template, data_length);

    lam0 = dim1_data_length/2.0;

    while (lam0 < data_length)
    {
        /*
          Create the chirp with the given parameters, using the time-domain
          chirp generator
        */
        tau0 = (5.0/3.0)/f_seismic*lam0;
        generateFrequencyDomainChirp(tc, lam0, chirp);

        /* Calculate its norm */
        chirp_norm = sqrt(inner(chirp, chirp));

        /* Convert tc to sample number */
        lambda = SRATE*tc/2;

        printf("----------------------------------------------------------\n");
        printf("Chirp parameters:    tc     = %e, tau0 = %e\n",
               tc, tau0);
        printf("Unitless parameters: lambda = %e, lam0 = %e\n",
               lambda, lam0);

        /*
          Plug the chirp into the FCT, returning the value of
          max(<chirp, template>), where the inner product is approximated
          using the FCT and the maximum is over all available templates.
          Also returns the estimated template parameters.
        */
        dim1_start = lam0 - dim1_data_length/2;
	inner_fct = 0.0;
        fct_lam_est = 0;
        fct_lam0_est = 0;
        calculateFCTInner(status, chirp, dim1_start,
                          &inner_fct, &fct_lam_est, &fct_lam0_est);

        printf("FCT parameters:      lambda = %d, lam0 = %d\n",
               fct_lam_est, fct_lam0_est);

	inner_dct = 0.0;
        dct_lam_est = 0;
        dct_lam0_est = 0;
        calculateDCTInner(status, chirp, dim1_start,
                          &inner_dct, &dct_lam_est, &dct_lam0_est);

        printf("DCT parameters:      lambda = %d, lam0 = %d\n",
               dct_lam_est, dct_lam0_est);

	/* Now calculate the SNR's */

	/* FCT SNR */
        generateFCTTemplate(fct_lam_est, fct_lam0_est, template);
        template_norm = sqrt(inner(template, template));
        SNR_fct = inner_fct/(chirp_norm*template_norm);

	tmp_inner = inner(chirp, template);

	printf("FCT INNER = %e\n", inner_fct);
	printf("TMP INNER = %e\n", tmp_inner);

	/* DCT SNR */
        generateFrequencyDomainChirp(2.0*dct_lam_est/SRATE, dct_lam0_est,
				     template);
        template_norm = sqrt(inner(template, template));
        SNR_dct = inner_dct/(chirp_norm*template_norm);

#if 0
        printf("estimated inner product = %e\n", inner_fct);
        printf("actual    inner product = %e\n", inner_dct);

        printf("Fraction of SNR recovered from estimated "
               "inner product = %e\n", SNR_est);
        printf("Fraction of SNR recovered from actual "
               "inner product = %e\n", SNR_actual);

        printf("           lambda estimate = %5d\n", lam_est);
        printf("           lam0 estimate = %5d\n", lam0_est);
#endif

        printf("DCT overlap   FCT overlap\n");
        printf("%e  %e\n", SNR_dct, SNR_fct);

        ++lam0Idx;
        lam0 = dim1_data_length/2.0 + lam0Idx*deltaLam0;
    }

    LALCDestroyVector(status, &chirp);
    LALCDestroyVector(status, &template);
}

static
void
setupGlobals(LALStatus* const status)
{
    const LALCreateFCTPlanInput planIn = {
        data_length,
        NUMBER_OF_DIMENSIONS,
        dimension_0_stride,
        /* Array of phase fns */
        {
            phi0
        }
    };

    const LALFCTSetUnitsInput setUnitsIn = {
        OFFSET,
        SDELTA
    };

    UINT4 k = 0;
    UINT8 output_data_size = 0;

    /* Set up the plan and output area */
    LALCreateFCTPlan(status, &gfctPlan, &planIn);
    LALFCTSetUnits(status, &setUnitsIn, gfctPlan);

    /*
      Initialise a data cube of the size we plan to use, so
      we can allocate the correct output size
    */
    cube[0].start_locations[0] = 0;
    cube[0].end_locations[0]   = data_length;

    cube[0].start_locations[1] = 0;
    cube[0].end_locations[1]   = dim1_data_length;

    cube[0].stride[0] = 1;
    cube[0].stride[1] = 1;

    dataCubesIn.data_cube = &cube[0];
    dataCubesIn.num_data_cubes = NUM_DATA_CUBES;

    /* Get the required output size */
    LALFCTSetDataCubes(status, &dataCubesIn, gfctPlan);
    LALFCTOutputDataSize(status, &output_data_size, gfctPlan);

    /* Allocate an output vector of this size */
    LALCCreateVector(status, &gfct_out, output_data_size);

    /* Allocate a vector for the noise spectrum */
    LALDCreateVector(status, &gSh, TSERIES_LENGTH/2);

    /* Initialise the noise model */
    for (k = 0; k < gSh->length; ++k)
    {
        const REAL8 f = OFFSET + DELTA*k;
        REAL8 psd = 0.0;

	if (f < f_seismic)
	{
	    gSh->data[k] = 1.0;
	    gSh->data[k] = 1.0e+99;
	}
	else
	{
	    LALLIGOIPsd(status, &psd, f);
	    gSh->data[k] = 1.0;
	    gSh->data[k] = psd;
	}
    }
}

static
void
Cleanup(LALStatus* const status)
{
    /* Deallocate the arrays and plan */
    LALCDestroyVector(status, &gfct_out);
    LALDDestroyVector(status, &gSh);

    LALDestroyFCTPlan(status, &gfctPlan);
}

int main(void);

int
main()
{
    LALStatus* const status = LALCalloc(1, sizeof(*status));

    setupGlobals(status);

    CalculateSNRLoss(status);

    Cleanup(status);

    LALFree(status);
    LALCheckMemoryLeaks();

    return 0;
}
#endif /* fftw2 implementation */
