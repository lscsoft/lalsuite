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
#define WL 1024

/* number of training set waveforms */
#define TSSIZE 1000

#define TOLERANCE 10e-12

/* tolerance allow for fractional percentage log likelihood difference */
#define LTOL 0.1

/* simple inspiral phase model */
double calc_phase(double frequency, double Mchirp);

/* model for a purely real frequency domain inspiral-like signal */
double real_model(double frequency, double Mchirp, double modperiod);

/* model for a complex frequency domain inspiral-like signal */
COMPLEX16 imag_model(double frequency, double Mchirp, double modperiod);

double calc_phase(double frequency, double Mchirp){
  return (-0.25*LAL_PI + ( 3./( 128. * pow(Mchirp*LAL_MTSUN_SI*LAL_PI*frequency, 5./3.) ) ) );
}

double real_model(double frequency, double Mchirp, double modperiod){
  return ( pow(frequency, -7./6.) * pow(Mchirp*LAL_MTSUN_SI,5./6.) * cos(calc_phase(frequency,Mchirp)) )*sin(LAL_TWOPI*frequency/modperiod);
}

COMPLEX16 imag_model(double frequency, double Mchirp, double modperiod){
  return ( pow(frequency, -7./6.) * pow(Mchirp*LAL_MTSUN_SI,5./6.) * cexp(I*calc_phase(frequency,Mchirp)) )*sin(LAL_TWOPI*frequency/modperiod);
}

int main(void) {
  REAL8Array *TS = NULL, *TSquad = NULL, *cTSquad = NULL;  /* the training set of real waveforms (and quadratic model) */
  COMPLEX16Array *cTS = NULL;              /* the training set of complex waveforms */

  size_t TSsize;  /* the size of the training set (number of waveforms) */
  size_t wl;      /* the length of each waveform */
  size_t k = 0, j = 0, i = 0;

  REAL8Array *RBlinear = NULL, *RBquad = NULL, *cRBquad = NULL;       /* the real reduced basis set */
  COMPLEX16Array *cRBlinear = NULL;   /* the complex reduced basis set */

  LALInferenceREALROQInterpolant *interp = NULL, *interpQuad = NULL, *cinterpQuad = NULL;
  LALInferenceCOMPLEXROQInterpolant *cinterp = NULL;

  double tolerance = TOLERANCE; /* tolerance for reduced basis generation loop */

  TSsize = TSSIZE;
  wl = WL;

  /* allocate memory for training set */
  UINT4Vector *TSdims = XLALCreateUINT4Vector( 2 );
  TSdims->data[0] = TSsize;
  TSdims->data[1] = wl;
  TS = XLALCreateREAL8Array( TSdims );
  TSquad = XLALCreateREAL8Array( TSdims );
  cTS = XLALCreateCOMPLEX16Array( TSdims );
  cTSquad = XLALCreateREAL8Array( TSdims );

  gsl_matrix_view TSview, TSviewquad, cTSviewquad;
  TSview = gsl_matrix_view_array( TS->data, TSdims->data[0], TSdims->data[1] );
  TSviewquad = gsl_matrix_view_array( TSquad->data, TSdims->data[0], TSdims->data[1] );
  cTSviewquad = gsl_matrix_view_array( cTSquad->data, TSdims->data[0], TSdims->data[1] );
  gsl_matrix_complex_view cTSview;
  cTSview = gsl_matrix_complex_view_array( (double *)cTS->data, TSdims->data[0], TSdims->data[1] );

  XLALDestroyUINT4Vector( TSdims );

  /* the waveform model is just a simple chirp so set up chirp mass range for training set */
  double fmin0 = 48, fmax0 = 256, f0 = 0., m0 = 0.;
  double df = (fmax0-fmin0)/(wl-1.); /* model time steps */
  double Mcmax = 2., Mcmin = 1.5, Mc = 0.;
  double periodmax = 1./99.995, periodmin = 1./100., modperiod = 0.;

  REAL8Vector *fweights = XLALCreateREAL8Vector( 1 );
  fweights->data[0] = df;

  REAL8Vector *freqs = XLALCreateREAL8Vector( wl );

  /* random number generator setup */
  const gsl_rng_type *T;
  gsl_rng *r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  /* set up training sets (one real and one complex) */
  for ( k=0; k < TSsize; k++ ){
    Mc = pow(pow(Mcmin, 5./3.) + (double)k*(pow(Mcmax, 5./3.)-pow(Mcmin, 5./3.))/((double)TSsize-1), 3./5.);
    modperiod = gsl_ran_flat(r, periodmin, periodmax);

    for ( j=0; j < wl; j++ ){
      f0 = fmin0 + (double)j*(fmax0-fmin0)/((double)wl-1.);

      gsl_complex gctmp;
      COMPLEX16 ctmp;
      m0 = real_model(f0, Mc, modperiod);
      ctmp = imag_model(f0, Mc, modperiod);
      GSL_SET_COMPLEX(&gctmp, creal(ctmp), cimag(ctmp));
      freqs->data[j] = f0;
      gsl_matrix_set(&TSview.matrix, k, j, m0);
      gsl_matrix_set(&TSviewquad.matrix, k, j, m0*m0);
      gsl_matrix_complex_set(&cTSview.matrix, k, j, gctmp);
      gsl_matrix_set(&cTSviewquad.matrix, k, j, creal(ctmp*conj(ctmp)));
    }
  }

  /* create reduced orthonormal basis from training set for linear part */
  REAL8 maxprojerr = 0.;
  maxprojerr = LALInferenceGenerateREAL8OrthonormalBasis(&RBlinear, fweights, tolerance, TS);
  fprintf(stderr, "No. linear nodes (real) = %d, %d x %d; Maximum projection err. = %le\n", RBlinear->dimLength->data[0], RBlinear->dimLength->data[0], RBlinear->dimLength->data[1], maxprojerr);
  maxprojerr = LALInferenceGenerateCOMPLEX16OrthonormalBasis(&cRBlinear, fweights, tolerance, cTS);
  fprintf(stderr, "No. linear nodes (complex) = %d, %d x %d; Maximum projection err. = %le\n", cRBlinear->dimLength->data[0], cRBlinear->dimLength->data[0], cRBlinear->dimLength->data[1], maxprojerr);
  maxprojerr = LALInferenceGenerateREAL8OrthonormalBasis(&RBquad, fweights, tolerance, TSquad);
  fprintf(stderr, "No. quadratic nodes (real)  = %d, %d x %d; Maximum projection err. = %le\n", RBquad->dimLength->data[0], RBquad->dimLength->data[0], RBquad->dimLength->data[1], maxprojerr);
  maxprojerr = LALInferenceGenerateREAL8OrthonormalBasis(&cRBquad, fweights, tolerance, cTSquad);
  fprintf(stderr, "No. quadratic nodes (complex)  = %d, %d x %d; Maximum projection err. = %le\n", cRBquad->dimLength->data[0], cRBquad->dimLength->data[0], cRBquad->dimLength->data[1], maxprojerr);

  /* free the training set */
  XLALDestroyREAL8Array( TS );
  XLALDestroyCOMPLEX16Array( cTS );
  XLALDestroyREAL8Array( TSquad );
  XLALDestroyREAL8Array( cTSquad );
  XLALDestroyREAL8Vector( fweights );

  /* get the linear interpolant */
  interp = LALInferenceGenerateREALROQInterpolant(RBlinear);
  cinterp = LALInferenceGenerateCOMPLEXROQInterpolant(cRBlinear);

  /* get the quadratic interpolant */
  interpQuad = LALInferenceGenerateREALROQInterpolant(RBquad);
  cinterpQuad = LALInferenceGenerateREALROQInterpolant(cRBquad);

  /* free the reduced basis */
  XLALDestroyREAL8Array(RBlinear);
  XLALDestroyCOMPLEX16Array(cRBlinear);
  XLALDestroyREAL8Array(RBquad);
  XLALDestroyREAL8Array(cRBquad);

  /* now get the terms for the likelihood with and without the reduced order quadrature
   * and do some timing tests */

  /* create the model dot model weights */
  REAL8Vector *vars = XLALCreateREAL8Vector( 1 );
  vars->data[0] = 1.;

  REAL8Vector *mmw = LALInferenceGenerateQuadraticWeights(interpQuad->B, vars);
  REAL8Vector *cmmw = LALInferenceGenerateQuadraticWeights(cinterpQuad->B, vars);

  /* let's create some Gaussian random data */
  REAL8Vector *data = XLALCreateREAL8Vector( wl );
  COMPLEX16Vector *cdata = XLALCreateCOMPLEX16Vector( wl );
  for ( i=0; i<wl; i++ ){
    data->data[i] = gsl_ran_gaussian(r, 1.0);                               /* real data */
    cdata->data[i] = gsl_ran_gaussian(r, 1.0) + I*gsl_ran_gaussian(r, 1.0); /* complex data */
  }

  /* create the data dot model weights */
  REAL8Vector *dmw = LALInferenceGenerateREAL8LinearWeights(interp->B, data, vars);
  COMPLEX16Vector *cdmw = LALInferenceGenerateCOMPLEX16LinearWeights(cinterp->B, cdata, vars);

  XLALDestroyREAL8Vector( vars );

  /* pick a chirp mass and generate a model to compare likelihoods */
  double randMc = 1.873; /* a random frequency to create a model */
  double randmp = 1./99.9989;

  gsl_vector *modelfull = gsl_vector_alloc(wl);
  REAL8Vector *modelreduced = XLALCreateREAL8Vector( interp->B->dimLength->data[0] );
  REAL8Vector *modelreducedquad = XLALCreateREAL8Vector( interpQuad->B->dimLength->data[0] );
  gsl_vector_complex *cmodelfull = gsl_vector_complex_alloc(wl);
  COMPLEX16Vector *cmodelreduced = XLALCreateCOMPLEX16Vector( cinterp->B->dimLength->data[0] );
  REAL8Vector *cmodelreducedquad = XLALCreateREAL8Vector( cinterpQuad->B->dimLength->data[0] );

  /* create models */
  for ( i=0; i<wl; i++ ){
    /* models at all frequencies */
    gsl_vector_set(modelfull, i, real_model(freqs->data[i], randMc, randmp));

    COMPLEX16 cval = imag_model(freqs->data[i], randMc, randmp);
    gsl_complex gcval;
    GSL_SET_COMPLEX(&gcval, creal(cval), cimag(cval));
    gsl_vector_complex_set(cmodelfull, i, gcval);
  }

  /* models at interpolant nodes */
  for ( i=0; i<modelreduced->length; i++ ){ /* real model */
    REAL8 rm = real_model(freqs->data[interp->nodes[i]], randMc, randmp);
    modelreduced->data[i] = rm;
  }
  for ( i=0; i<modelreducedquad->length; i++ ){ /* real model */
    REAL8 rm = real_model(freqs->data[interpQuad->nodes[i]], randMc, randmp);
    modelreducedquad->data[i] = rm*rm;
  }
  for ( i=0; i<cmodelreduced->length; i++ ){ /* complex model */
    COMPLEX16 crm = imag_model(freqs->data[cinterp->nodes[i]], randMc, randmp);
    cmodelreduced->data[i] = crm;
  }
  for ( i=0; i<cmodelreducedquad->length; i++ ){ /* complex model */
    COMPLEX16 crm = imag_model(freqs->data[cinterpQuad->nodes[i]], randMc, randmp);
    cmodelreducedquad->data[i] = creal(crm*conj(crm));
  }

  XLALDestroyREAL8Vector( freqs );

  /* timing variables */
  struct timeval t1, t2, t3, t4;
  double dt1, dt2;

  /* start with the real model */
  /* get the model model term with the full model */
  REAL8 mmfull, mmred;
  gettimeofday(&t1, NULL);
  XLAL_CALLGSL( gsl_blas_ddot(modelfull, modelfull, &mmfull) ); /* real model */
  gettimeofday(&t2, NULL);

  /* now get it with the reduced order quadrature */
  gettimeofday(&t3, NULL);
  mmred = LALInferenceROQREAL8DotProduct(mmw, modelreducedquad);
  gettimeofday(&t4, NULL);

  dt1 = (double)((t2.tv_sec + t2.tv_usec*1.e-6) - (t1.tv_sec + t1.tv_usec*1.e-6));
  dt2 = (double)((t4.tv_sec + t4.tv_usec*1.e-6) - (t3.tv_sec + t3.tv_usec*1.e-6));
  fprintf(stderr, "Real Signal:\n - M dot M (full) = %le [%.9lf s], M dot M (reduced) = %le [%.9lf s], time ratio = %lf\n", mmfull, dt1, mmred, dt2, dt1/dt2);

  /* get the data model term with the full model */
  REAL8 dmfull, dmred;
  gsl_vector_view dataview = gsl_vector_view_array(data->data, wl);
  gettimeofday(&t1, NULL);
  XLAL_CALLGSL( gsl_blas_ddot(&dataview.vector, modelfull, &dmfull) );
  gettimeofday(&t2, NULL);

  /* now get it with the reduced order quadrature */
  gettimeofday(&t3, NULL);
  dmred = LALInferenceROQREAL8DotProduct(dmw, modelreduced);
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

  XLALDestroyREAL8Vector(data);
  gsl_vector_free(modelfull);
  XLALDestroyREAL8Vector(modelreduced);
  XLALDestroyREAL8Vector(modelreducedquad);
  XLALDestroyREAL8Vector(mmw);
  XLALDestroyREAL8Vector(dmw);

  /* check log likelihood difference is within tolerance */
  if ( Lfrac > LTOL ) { return 1; }

  /* now do the same with the complex model */
  /* get the model model term with the full model */
  REAL8 cmmred, cmmfull;
  gsl_complex cmmfulltmp;
  gettimeofday(&t1, NULL);
  XLAL_CALLGSL( gsl_blas_zdotc(cmodelfull, cmodelfull, &cmmfulltmp) ); /* complex model */
  cmmfull = GSL_REAL(cmmfulltmp);
  gettimeofday(&t2, NULL);

  gettimeofday(&t3, NULL);
  cmmred = LALInferenceROQREAL8DotProduct(cmmw, cmodelreducedquad);
  gettimeofday(&t4, NULL);
  dt1 = (double)((t2.tv_sec + t2.tv_usec*1.e-6) - (t1.tv_sec + t1.tv_usec*1.e-6));
  dt2 = (double)((t4.tv_sec + t4.tv_usec*1.e-6) - (t3.tv_sec + t3.tv_usec*1.e-6));
  fprintf(stderr, "Complex Signal:\n - M dot M (full) = %le [%.9lf s], M dot M (reduced) = %le [%.9lf s], time ratio = %lf\n", cmmfull, dt1, cmmred, dt2, dt1/dt2);

  COMPLEX16 cdmfull, cdmred;
  gsl_complex cdmfulltmp;
  gsl_vector_complex_view cdataview = gsl_vector_complex_view_array((double*)cdata->data, wl);
  gettimeofday(&t1, NULL);
  XLAL_CALLGSL( gsl_blas_zdotc(&cdataview.vector, cmodelfull, &cdmfulltmp) );
  cdmfull = GSL_REAL(cdmfulltmp) + I*GSL_IMAG(cdmfulltmp);
  gettimeofday(&t2, NULL);

  gettimeofday(&t3, NULL);
  cdmred = LALInferenceROQCOMPLEX16DotProduct(cdmw, cmodelreduced);
  gettimeofday(&t4, NULL);

  dt1 = (double)((t2.tv_sec + t2.tv_usec*1.e-6) - (t1.tv_sec + t1.tv_usec*1.e-6));
  dt2 = (double)((t4.tv_sec + t4.tv_usec*1.e-6) - (t3.tv_sec + t3.tv_usec*1.e-6));
  fprintf(stderr, " - D dot M (full) = %le [%.9lf s], D dot M (reduced) = %le [%.9lf s], time ratio = %lf\n", creal(cdmfull), dt1, creal(cdmred), dt2, dt1/dt2);

  /* check difference in log likelihoods */
  Lfull = cmmfull - 2.*creal(cdmfull);
  Lred = cmmred - 2.*creal(cdmred);
  Lfrac = 100.*fabs(Lfull-Lred)/fabs(Lfull); /* fractional log likelihood difference (in %) */

  fprintf(stderr, " - Fractional difference in log likelihoods = %lf%%\n", Lfrac);

  XLALDestroyCOMPLEX16Vector(cdata);
  gsl_vector_complex_free(cmodelfull);
  XLALDestroyCOMPLEX16Vector(cmodelreduced);
  XLALDestroyREAL8Vector(cmodelreducedquad);
  XLALDestroyREAL8Vector(cmmw);
  XLALDestroyCOMPLEX16Vector(cdmw);

  LALInferenceRemoveREALROQInterpolant( interp );
  LALInferenceRemoveCOMPLEXROQInterpolant( cinterp );
  LALInferenceRemoveREALROQInterpolant( interpQuad );
  LALInferenceRemoveREALROQInterpolant( cinterpQuad );

  /* check log likelihood difference is within tolerance */
  if ( Lfrac > LTOL ) { return 1; }

  LALCheckMemoryLeaks();

  return 0;
}
