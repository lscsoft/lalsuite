#include <lal/LALInferenceGenerateROQ.h>
#include <lal/LALConstants.h>
#include <gsl/gsl_randist.h>

#include <time.h>
#include <math.h>

/* check whether to include omp.h for use of multiple cores */
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

/* test waveform length */
#define WL 4096

/* number of training set waveforms */
#define TSSIZE 1000

#define TOLERANCE 10e-12

/* tolerance allow for fractional percentage log likelihood difference */
#define LTOL 0.1

/* simple inspiral phase model */
double calc_phase(double frequency, double Mchirp);

/* model for a purely real frequency domain inspiral-like signal */
double real_model(double frequency, double Mchirp);

/* model for a complex frequency domain inspiral-like signal */
COMPLEX16 imag_model(double frequency, double Mchirp);

double calc_phase(double frequency, double Mchirp){
  return (-0.25*LAL_PI + ( 3./( 128. * pow(Mchirp*LAL_MTSUN_SI*LAL_PI*frequency, 5./3.) ) ) );
}

double real_model(double frequency, double Mchirp){
  return ( pow(frequency, -7./6.) * pow(Mchirp*LAL_MTSUN_SI,5./6.) * cos(calc_phase(frequency,Mchirp)) );
}

COMPLEX16 imag_model(double frequency, double Mchirp){
  return ( pow(frequency, -7./6.) * pow(Mchirp*LAL_MTSUN_SI,5./6.) * cexp(I*calc_phase(frequency,Mchirp)) );
}

int main(void) {
  gsl_matrix *TS; /* the training set of real waveforms */
  gsl_matrix_complex *cTS; /* the training set of complex waveforms */

  size_t TSsize;  /* the size of the training set (number of waveforms) */
  size_t wl;      /* the length of each waveform */
  size_t k = 0, j = 0, i = 0, nbases = 0, cnbases = 0;

  REAL8 *RB = NULL;        /* the real reduced basis set */
  COMPLEX16 *cRB = NULL;   /* the complex reduced basis set */
  LALInferenceREALROQInterpolant *interp = NULL;
  LALInferenceCOMPLEXROQInterpolant *cinterp = NULL;

  gsl_vector *freqs;

  double tolerance = TOLERANCE; /* tolerance for reduced basis generation loop */

  TSsize = TSSIZE;
  wl = WL;

  /* allocate memory for training set */
  TS = gsl_matrix_calloc(TSsize, wl);
  cTS = gsl_matrix_complex_calloc(TSsize, wl);

  /* the waveform model is just a simple chirp so set up chirp mass range for training set */
  double fmin0 = 48, fmax0 = 256, f0 = 0., m0 = 0.;
  double df = (fmax0-fmin0)/(wl-1.); /* model time steps */
  freqs = gsl_vector_alloc(wl); /* times at which to calculate the model */
  double Mcmax = 2., Mcmin = 1.5, Mc = 0.;

  gsl_vector_view fweights = gsl_vector_view_array(&df, 1);

  /* set up training sets (one real and one complex) */
  for ( k=0; k < TSsize; k++ ){
    Mc = pow(pow(Mcmin, 5./3.) + (double)k*(pow(Mcmax, 5./3.)-pow(Mcmin, 5./3.))/((double)TSsize-1), 3./5.);

    for ( j=0; j < wl; j++ ){
      f0 = fmin0 + (double)j*(fmax0-fmin0)/((double)wl-1.);

      gsl_complex gctmp;
      COMPLEX16 ctmp;
      m0 = real_model(f0, Mc);
      ctmp = imag_model(f0, Mc);
      GSL_SET_COMPLEX(&gctmp, creal(ctmp), cimag(ctmp));
      gsl_vector_set(freqs, j, f0);
      gsl_matrix_set(TS, k, j, m0);
      gsl_matrix_complex_set(cTS, k, j, gctmp);
    }
  }

  /* create reduced orthonormal basis from training set */
  if ( (RB = LALInferenceGenerateREAL8OrthonormalBasis(&fweights.vector, tolerance, TS, &nbases)) == NULL){
    fprintf(stderr, "Error... problem producing basis\n");
    return 1;
  }

  if ( (cRB = LALInferenceGenerateCOMPLEX16OrthonormalBasis(&fweights.vector, tolerance, cTS, &cnbases)) == NULL){
    fprintf(stderr, "Error... problem producing basis\n");
    return 1;
  }

  /* free the training set */
  gsl_matrix_free(TS);
  gsl_matrix_complex_free(cTS);

  gsl_matrix_view RBview = gsl_matrix_view_array(RB, nbases, wl);
  gsl_matrix_complex_view cRBview = gsl_matrix_complex_view_array((double*)cRB, cnbases, wl);

  fprintf(stderr, "No. nodes (real)  = %zu, %zu x %zu\n", nbases, RBview.matrix.size1, RBview.matrix.size2);
  fprintf(stderr, "No. nodes (complex)  = %zu, %zu x %zu\n", cnbases, cRBview.matrix.size1, cRBview.matrix.size2);

  /* get the interpolant */
  interp = LALInferenceGenerateREALROQInterpolant(&RBview.matrix);
  cinterp = LALInferenceGenerateCOMPLEXROQInterpolant(&cRBview.matrix);

  /* free the reduced basis */
  XLALFree(RB);
  XLALFree(cRB);

  /* now get the terms for the likelihood with and without the reduced order quadrature
   * and do some timing tests */

  /* create the model dot model weights */
  REAL8 varval = 1.;
  gsl_vector_view vars = gsl_vector_view_array(&varval, 1);

  gsl_matrix *mmw = LALInferenceGenerateREALModelModelWeights(interp->B, &vars.vector);
  gsl_matrix_complex *cmmw = LALInferenceGenerateCOMPLEXModelModelWeights(cinterp->B, &vars.vector);

  /* let's create some Gaussian random data */
  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  REAL8 *data = XLALCalloc(wl, sizeof(REAL8));
  COMPLEX16 *cdata = XLALCalloc(wl, sizeof(COMPLEX16));
  for ( i=0; i<wl; i++ ){
    data[i] = gsl_ran_gaussian(r, 1.0);                               /* real data */
    cdata[i] = gsl_ran_gaussian(r, 1.0) + I*gsl_ran_gaussian(r, 1.0); /* complex data */
  }

  /* create the data dot model weights */
  gsl_vector_view dataview = gsl_vector_view_array(data, wl);
  gsl_vector *dmw = LALInferenceGenerateREAL8DataModelWeights(interp->B, &dataview.vector, &vars.vector);

  gsl_vector_complex_view cdataview = gsl_vector_complex_view_array((double*)cdata, wl);
  gsl_vector_complex *cdmw = LALInferenceGenerateCOMPLEX16DataModelWeights(cinterp->B, &cdataview.vector, &vars.vector);

  /* pick a chirp mass and generate a model to compare likelihoods */
  double randMc = 1.873; /* a random frequency to create a model */

  gsl_vector *modelfull = gsl_vector_alloc(wl);
  gsl_vector *modelreduced = gsl_vector_alloc(nbases);
  gsl_vector_complex *cmodelfull = gsl_vector_complex_alloc(wl);
  gsl_vector_complex *cmodelreduced = gsl_vector_complex_alloc(cnbases);

  /* create models */
  for ( i=0; i<wl; i++ ){
    /* models at all frequencies */
    gsl_vector_set(modelfull, i, real_model(gsl_vector_get(freqs, i), randMc));

    COMPLEX16 cval = imag_model(gsl_vector_get(freqs, i), randMc);
    gsl_complex gcval;
    GSL_SET_COMPLEX(&gcval, creal(cval), cimag(cval));
    gsl_vector_complex_set(cmodelfull, i, gcval);
  }

  /* models at interpolant nodes */
  for ( i=0; i<nbases; i++ ){ /* real model */
    gsl_vector_set(modelreduced, i, real_model(gsl_vector_get(freqs, interp->nodes[i]), randMc));
  }
  for ( i=0; i<cnbases; i++ ){ /* complex model */
    COMPLEX16 cval = imag_model(gsl_vector_get(freqs, cinterp->nodes[i]), randMc);
    gsl_complex gcval;
    GSL_SET_COMPLEX(&gcval, creal(cval), cimag(cval));
    gsl_vector_complex_set(cmodelreduced, i, gcval);
  }

  /* timing variables */
  struct timeval t1, t2, t3, t4;
  double dt1, dt2;

  /* start with the real model */
  /* get the model model term with the full model */
  REAL8 mmfull, mmred;
  gettimeofday(&t1, NULL);
  XLAL_CALLGSL( gsl_blas_ddot(modelfull, modelfull, &mmfull) );        /* real model */
  gettimeofday(&t2, NULL);

  /* now get it with the reduced order quadrature */
  gettimeofday(&t3, NULL);
  mmred = LALInferenceROQREAL8ModelDotModel(mmw, modelreduced);
  gettimeofday(&t4, NULL);

  dt1 = (double)((t2.tv_sec + t2.tv_usec*1.e-6) - (t1.tv_sec + t1.tv_usec*1.e-6));
  dt2 = (double)((t4.tv_sec + t4.tv_usec*1.e-6) - (t3.tv_sec + t3.tv_usec*1.e-6));
  fprintf(stderr, "Real Signal:\n - M dot M (full) = %le [%.9lf s], M dot M (reduced) = %le [%.9lf s], time ratio = %lf\n", mmfull, dt1, mmred, dt2, dt1/dt2);

  /* get the data model term with the full model */
  REAL8 dmfull, dmred;
  gettimeofday(&t1, NULL);
  XLAL_CALLGSL( gsl_blas_ddot(&dataview.vector, modelfull, &dmfull) );
  gettimeofday(&t2, NULL);

  /* now get it with the reduced order quadrature */
  gettimeofday(&t3, NULL);
  dmred = LALInferenceROQREAL8DataDotModel(dmw, modelreduced);
  gettimeofday(&t4, NULL);

  dt1 = (double)((t2.tv_sec + t2.tv_usec*1.e-6) - (t1.tv_sec + t1.tv_usec*1.e-6));
  dt2 = (double)((t4.tv_sec + t4.tv_usec*1.e-6) - (t3.tv_sec + t3.tv_usec*1.e-6));
  fprintf(stderr, " - D dot M (full) = %le [%.9lf s], D dot M (reduced) = %le [%.9lf s], time ratio = %lf\n", dmfull, dt1, dmred, dt2, dt1/dt2);

  /* check difference in log likelihoods */
  double Lfull, Lred, Lfrac;

  Lfull = mmfull - 2.*dmfull;
  Lred = mmred - 2.*dmred;
  Lfrac = 100.*fabs(Lfull-Lred)/fabs(Lfull); /* fractional log likelihood difference (in %) */

  fprintf(stderr, " - Fractional difference in log likelihoods = %lf%%\n", Lfrac);

  XLALFree(data);
  gsl_vector_free(modelfull);
  gsl_vector_free(modelreduced);
  gsl_matrix_free(mmw);
  gsl_vector_free(dmw);

  /* check log likelihood difference is within tolerance */
  if ( Lfrac > LTOL ) { return 1; }

  /* now do the same with the complex model */
  /* get the model model term with the full model */
  COMPLEX16 cmmred, cmmfull;
  gsl_complex cmmfulltmp;
  gettimeofday(&t1, NULL);
  XLAL_CALLGSL( gsl_blas_zdotc(cmodelfull, cmodelfull, &cmmfulltmp) ); /* complex model */
  cmmfull = GSL_REAL(cmmfulltmp) + I*GSL_IMAG(cmmfulltmp);
  gettimeofday(&t2, NULL);

  gettimeofday(&t3, NULL);
  cmmred = LALInferenceROQCOMPLEX16ModelDotModel(cmmw, cmodelreduced);
  gettimeofday(&t4, NULL);

  dt1 = (double)((t2.tv_sec + t2.tv_usec*1.e-6) - (t1.tv_sec + t1.tv_usec*1.e-6));
  dt2 = (double)((t4.tv_sec + t4.tv_usec*1.e-6) - (t3.tv_sec + t3.tv_usec*1.e-6));
  fprintf(stderr, "Complex Signal:\n - M dot M (full) = %le [%.9lf s], M dot M (reduced) = %le [%.9lf s], time ratio = %lf\n", creal(cmmfull), dt1, creal(cmmred), dt2, dt1/dt2);

  COMPLEX16 cdmfull, cdmred;
  gsl_complex cdmfulltmp;
  gettimeofday(&t1, NULL);
  XLAL_CALLGSL( gsl_blas_zdotc(&cdataview.vector, cmodelfull, &cdmfulltmp) );
  cdmfull = GSL_REAL(cdmfulltmp) + I*GSL_IMAG(cdmfulltmp);
  gettimeofday(&t2, NULL);

  gettimeofday(&t3, NULL);
  cdmred = LALInferenceROQCOMPLEX16DataDotModel(cdmw, cmodelreduced);
  gettimeofday(&t4, NULL);

  dt1 = (double)((t2.tv_sec + t2.tv_usec*1.e-6) - (t1.tv_sec + t1.tv_usec*1.e-6));
  dt2 = (double)((t4.tv_sec + t4.tv_usec*1.e-6) - (t3.tv_sec + t3.tv_usec*1.e-6));
  fprintf(stderr, " - D dot M (full) = %le [%.9lf s], D dot M (reduced) = %le [%.9lf s], time ratio = %lf\n", creal(cdmfull), dt1, creal(cdmred), dt2, dt1/dt2);

  /* check difference in log likelihoods */
  Lfull = creal(cmmfull) - 2.*creal(cdmfull);
  Lred = creal(cmmred) - 2.*creal(cdmred);
  Lfrac = 100.*fabs(Lfull-Lred)/fabs(Lfull); /* fractional log likelihood difference (in %) */

  fprintf(stderr, " - Fractional difference in log likelihoods = %lf%%\n", Lfrac);

  XLALFree(cdata);
  gsl_vector_complex_free(cmodelfull);
  gsl_vector_complex_free(cmodelreduced);
  gsl_matrix_complex_free(cmmw);
  gsl_vector_complex_free(cdmw);
  LALInferenceRemoveREALROQInterpolant( interp );
  LALInferenceRemoveCOMPLEXROQInterpolant( cinterp );

  /* check log likelihood difference is within tolerance */
  if ( Lfrac > LTOL ) { return 1; }

  return 0;
}
