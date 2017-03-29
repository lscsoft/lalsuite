/*
 *  Copyright (C) 2014 Michael Puerrer, John Veitch
 *  Reduced Order Model for SEOBNR
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

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <stdbool.h>
#include <alloca.h>
#include <string.h>
#include <libgen.h>


#include <gsl/gsl_errno.h>
#include <gsl/gsl_chebyshev.h>

// #include <gsl/gsl_blas.h>
// #include <gsl/gsl_min.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/FrequencySeries.h>
#include <lal/Date.h>
#include <lal/StringInput.h>
#include <lal/Sequence.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_spline.h>
//needed for read_vector and nudge
#include "LALSimIMRSEOBNRROMUtilities.c"

#include <lal/LALConfig.h>
#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
#endif


/************** Model parameters **************/

//Prepending G to make clear it is a global variable
// #define Gntimes 69072 //(will be reduced in the future)
// #define Gnamp 21
#define Gntimes 73624 //(will be reduced in the future)
#define Gnamp 12
#define Gnphase 7
#define Gnq 16
#define Gnlambda1 16
#define Gnlambda2 16

#ifdef LAL_PTHREAD_LOCK
static pthread_once_t TEOBResumROM_is_initialized = PTHREAD_ONCE_INIT;
#endif

static const double Gparams_min[] = {0.5,50.,50.}; //qmin,lambda1min,lambda2min
static const double Gparams_max[] = {1.0,5000.,5000.}; //qmax,lambda1max,lambda2max

/*************** type definitions ******************/

typedef struct tagTEOBResumROMdataDS_coeff
{
  gsl_vector* c_amp;
  gsl_vector* c_phi;
} TEOBResumROMdataDS_coeff;

struct tagTEOBResumROMdataDS_submodel
{
  gsl_vector* cvec_amp;      // amplitude projection coefficients
  gsl_vector* cvec_phi;      // phase projection coefficients
  gsl_matrix *Bamp;          // Reduced SVD basis for amplitude
  gsl_matrix *Bphi;          // Reduced SVD basis for phase
  gsl_vector* times;  // AMplitude prefactor coefficient
  const double *params_min;
  const double *params_max;
  int n_amp;                 // Number frequency points for amplitude
  int n_phi;                 // Number of frequency points for phase
  int nq, nl1, nl2, ntimes;         // Number of points in eta, chi1, chi2
};
typedef struct tagTEOBResumROMdataDS_submodel TEOBResumROMdataDS_submodel;

struct tagTEOBResumROMdataDS
{
  UINT4 setup;
  TEOBResumROMdataDS_submodel* sub1;
  TEOBResumROMdataDS_submodel* sub2;
  TEOBResumROMdataDS_submodel* sub3;
};
typedef struct tagTEOBResumROMdataDS TEOBResumROMdataDS;

static TEOBResumROMdataDS __lalsim_TEOBResumROMDS_data;

typedef int (*load_dataPtr)(const char*, gsl_vector *, gsl_vector *, gsl_matrix *, gsl_matrix *, gsl_vector *);

/**************** Internal functions **********************/

static void TEOBResumROM_Init_LALDATA(void);
static int TEOBResumROM_Init(const char dir[]);
static bool TEOBResumROM_IsSetup(void);

static int TEOBResumROMdataDS_Init(TEOBResumROMdataDS *romdata, const char dir[]);
static void TEOBResumROMdataDS_Cleanup(TEOBResumROMdataDS *romdata);

static int TEOBResumROMdataDS_Init_submodel(
  TEOBResumROMdataDS_submodel **submodel,
  const int n_amp,
  const int n_phi,
  const int nq,
  const int nl1,
  const int nl2,
  const int ntimes,
  const double *params_min,
  const double *params_max,
  const char dir[],
  load_dataPtr load_data
);


static double gsl_cheb_evaluate_polynomial(int n, double x);
static double gsl_cheb_eval_3d(gsl_vector *c_ijk, int nx, int ny, int nz, double x, double y, double z);
// static double chebev_3d(gsl_vector *c_ijk, int nx, int ny, int nz, double x, double y, double z, const double xyz_min[], const double xyz_max[]);
static int chebyshev_interpolation3d(
  double q,
  double lambda1,
  double lambda2,
  int nx, int ny, int nz,
  gsl_vector *cvec_amp,
  gsl_vector *cvec_phi,
  int nk_amp,
  int nk_phi,
  const double xyz_min[],
  const double xyz_max[],
  gsl_vector *interp_amp,
  gsl_vector *interp_phi);

static void TEOBResumROMdataDS_Cleanup_submodel(TEOBResumROMdataDS_submodel *submodel);

static int TEOBResumROMCore(
  REAL8TimeSeries **hPlus,
  REAL8TimeSeries **hCross,
  double phiRef,
  double deltaT,
  double fLow,
  double distance,
  double inclination,
  double Mtot_sec,
  double eta,
  double lambda1,
  double lambda2
);

static int load_data_romeos(const char dir[], gsl_vector *cvec_amp, gsl_vector *cvec_phi, gsl_matrix *Bamp, gsl_matrix *Bphi, gsl_vector *times);

/********************* Definitions begin here ********************/

/** Setup SEOBNRv2ROMDoubleSpin model using data files installed in dir
 */
static int TEOBResumROM_Init(const char dir[]) {
  if(__lalsim_TEOBResumROMDS_data.setup) {
    XLALPrintError("Error: DSTEOBResumROMdata was already set up!");
    XLAL_ERROR(XLAL_EFAILED);
  }

  TEOBResumROMdataDS_Init(&__lalsim_TEOBResumROMDS_data, dir);

  if(__lalsim_TEOBResumROMDS_data.setup) {
    return(XLAL_SUCCESS);
  }
  else {
    return(XLAL_EFAILED);
  }
}

/** Helper function to check if the SEOBNRv2ROMDoubleSpin model has been initialised */
static bool TEOBResumROM_IsSetup(void) {
  if(__lalsim_TEOBResumROMDS_data.setup)
    return true;
  else
    return false;
}

// Read binary ROM data for basis functions and coefficients for submodel 1
static int load_data_romeos(const char dir[], gsl_vector *cvec_amp, gsl_vector *cvec_phi, gsl_matrix *Bamp, gsl_matrix *Bphi, gsl_vector *times) {
  int ret = XLAL_SUCCESS;
  ret |= read_vector(dir, "TEOBResumROM_Amp_ciall.dat", cvec_amp);
  ret |= read_vector(dir, "TEOBResumROM_Phase_ciall.dat", cvec_phi);
  ret |= read_matrix(dir, "TEOBResumROM_Bamp_matrix.dat", Bamp);
  ret |= read_matrix(dir, "TEOBResumROM_Bphase_matrix.dat", Bphi);
  ret |= read_vector(dir, "TEOBResumROM_times.dat", times);
  return(ret);
}

/* Set up a new ROM submodel, using data contained in dir */
static int TEOBResumROMdataDS_Init_submodel(
  TEOBResumROMdataDS_submodel **submodel,
  const int n_amp,
  const int n_phi,
  const int nq,
  const int nl1,
  const int nl2,
  const int ntimes,
  const double *params_min,
  const double *params_max,
  const char dir[],
  load_dataPtr load_data
) {
  int ret = XLAL_FAILURE;

  if(!submodel) exit(1);
  /* Create storage for submodel structures */
  if (!*submodel)
    *submodel = XLALCalloc(1,sizeof(TEOBResumROMdataDS_submodel));
  else
    TEOBResumROMdataDS_Cleanup_submodel(*submodel);

  int N = nq*nl1*nl2;

  // Initalize actual ROM data
  (*submodel)->cvec_amp = gsl_vector_alloc(N*n_amp);
  (*submodel)->cvec_phi = gsl_vector_alloc(N*n_phi);
  (*submodel)->Bamp = gsl_matrix_alloc(n_amp, ntimes);
  (*submodel)->Bphi = gsl_matrix_alloc(n_phi, ntimes);
  (*submodel)->times = gsl_vector_alloc(ntimes);

  // Load ROM data for this submodel
  ret=load_data(dir, (*submodel)->cvec_amp, (*submodel)->cvec_phi, (*submodel)->Bamp, (*submodel)->Bphi, (*submodel)->times);

  // Initialize other members
  (*submodel)->n_amp = n_amp;
  (*submodel)->n_phi = n_phi;
  (*submodel)->nq = nq;
  (*submodel)->nl1 = nl1;
  (*submodel)->nl2 = nl2;
  (*submodel)->ntimes = ntimes;

  (*submodel)->params_min = params_min;
  (*submodel)->params_max = params_max;

  return ret;
}

/* Deallocate contents of the given TEOBResumROMdataDS_submodel structure */
static void TEOBResumROMdataDS_Cleanup_submodel(TEOBResumROMdataDS_submodel *submodel) {
  if(submodel->cvec_amp) gsl_vector_free(submodel->cvec_amp);
  if(submodel->cvec_phi) gsl_vector_free(submodel->cvec_phi);
  if(submodel->Bamp) gsl_matrix_free(submodel->Bamp);
  if(submodel->Bphi) gsl_matrix_free(submodel->Bphi);
  if(submodel->times) gsl_vector_free(submodel->times);
}

/* Set up a new ROM model, using data contained in dir */
int TEOBResumROMdataDS_Init(TEOBResumROMdataDS *romdata, const char dir[]) {
  int ret = XLAL_FAILURE;

  /* Create storage for structures */
  if(romdata->setup) {
    XLALPrintError("WARNING: You tried to setup the TEOBResumROM model that was already initialised. Ignoring\n");
    return (XLAL_FAILURE);
  }

  //gsl_set_error_handler(&err_handler);

  load_dataPtr load_data = &load_data_romeos;
  ret = TEOBResumROMdataDS_Init_submodel(&(romdata)->sub1, Gnamp, Gnphase, Gnq, Gnlambda1, Gnlambda2, Gntimes, Gparams_min, Gparams_max, dir, load_data);
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : submodel 1 loaded sucessfully.\n", __func__);

  if(XLAL_SUCCESS==ret)
    romdata->setup=1;
  else
    TEOBResumROMdataDS_Cleanup(romdata);

  return (ret);
}

/* Deallocate contents of the given TEOBResumROMdataDS structure */
static void TEOBResumROMdataDS_Cleanup(TEOBResumROMdataDS *romdata) {
  TEOBResumROMdataDS_Cleanup_submodel((romdata)->sub1);
  XLALFree((romdata)->sub1);
  (romdata)->sub1 = NULL;
  romdata->setup=0;
}


/* NOTE: this function does not give correct results yet */
/* Wrapper function to evaluate 3d chebyshev polynomial.
 * p(x,y,z) = \sum_{i,j,k} c_{ijk} T_i(x) T_j(y) T_k(z)
 * This is done by expressing the sum as follows:
 * f(x;y,z) = \sum_i h_i(y;z) T_i(x)
 *          = \sum_i [ \sum_j g_ij(z) T_j(y) ] T_i(x)
 *          = \sum_i \sum_j [ \sum_k c_ijk T_k(z) ] T_j(y) T_i(x)
 * The inner sum is evaluated with a custom implementation of Clenshaw's Recurrence to
 * avoid having to create the 1d coefficient vector for each ij. This implementation
 * picks out the coefficients from the 3d coefficient vector c_ijk.
 * The following sums are evaluated with gsl_cheb_eval at no additional cost.
 */
// static double chebev_3d(gsl_vector *c_ijk, int nx, int ny, int nz, double x, double y, double z, const double xyz_min[], const double xyz_max[]){
//
//   int i,j,k;
//   int index=0;
//   int Nyz = ny*nz ;
//
//   double gij_z,hi_y,f_xyz;
//   double d1=0.0,d2=0.0,sv,G,G2 ;
//
//   double a=xyz_min[2];
//   double b=xyz_max[2];
//   gsl_cheb_series *coeffs_1 = gsl_cheb_alloc(ny-1); //GSL assumes coeffs from c[0] to and including c[order]. So N = order+1
//   coeffs_1->a = xyz_min[1] ;
//   coeffs_1->b = xyz_max[1] ;
//   gsl_cheb_series *coeffs_2 = gsl_cheb_alloc(nx-1);
//   coeffs_2->a = xyz_min[0] ;
//   coeffs_2->b = xyz_max[0] ;
//
//   for (i=0;i<nx;i++){
//
//     for (j=0;j<ny;j++){
//
//       G = (2.0*z-a-b)/(b-a); //Change of variable
//       G2 = 2.0*G;
//       for (k=nz-1 ; k>=1 ; k--) { //Clenshaw's recurrence
//         sv = d1;
//         index = (i%Nyz)*Nyz + j*nz + k%nz ; //select coefficient corresponding to current ijk
//         d1 = G2*d1 - d2 + gsl_vector_get(c_ijk,index);
//         d2 = sv;
//       }
//       index = (i%Nyz)*Nyz + j*nz ; //c_ij0
//
//       //evaluate g_{ij}(z), which will be the coefficients for evaluation of h_i(y;z)
//       gij_z = G*d1 - d2 + gsl_vector_get(c_ijk,index); //NOTE: gsl returns G*d1 - d2 + 0.5 * cs->c[0]
//       coeffs_1->c[j] = gij_z ;
//     }
//
//     //evaluate h_i(y;z), which will be the coefficients for evaluation of f(x;y,z)
//     hi_y = gsl_cheb_eval(coeffs_1, y) + 0.5*coeffs_1->c[0] ;
//     coeffs_2->c[i] = hi_y ;
//   }
//
//   //evaluate f(x;y,z)
//   f_xyz = gsl_cheb_eval(coeffs_2, x) + 0.5*coeffs_2->c[0] ;
//
//   gsl_cheb_free(coeffs_1);
//   gsl_cheb_free(coeffs_2);
//
//   return f_xyz ;
//
// }

/* Evaluate the n-th Chebyshev polynomial at value x i.e. this is only T_n(x) without any coefficients or summing */
static double gsl_cheb_evaluate_polynomial(int n, double x){

  double Tnx = 0.0 ;
  //T_0(x)
  if (n==0){
    Tnx = 1.0 ;
  }
  //T_1(x)
  else if (n==1){
    Tnx = x ;
  }
  //T_2(x)
  else if (n==2){
    Tnx = 2.0*x*x - 1.0 ;
  }
  //T_n(x)
  else {
    double d1=0.0,d2=0.0,sv ;
    int k=0;

    sv = d1;
    d1 = 2.0*x*d1 - d2 + 1.0 ; //last coefficient is 1.0
    d2 = sv ;
    for (k=n-1 ; k>=1 ; k--) { //Clenshaw's recurrence
      sv = d1;
      d1 = 2.0*x*d1 - d2 ; //all coefficients except the last one are 0.0
      d2 = sv;
    }
    Tnx = x*d1 - d2 ; //c[0] is also 0.0
  }

  return Tnx ;

}

/* Temporary function that is about 5 times slower compared to the commented out chebev_3d. */
static double gsl_cheb_eval_3d(gsl_vector *c_ijk, int nx, int ny, int nz, double x, double y, double z){

  double sum=0.0;
  int i,j,k;
  int index=0;
  double Tix=0.,Tjy=0.,Tkz=0.;

  for (i=0;i<nx;i++){
    Tix=gsl_cheb_evaluate_polynomial(i,x);
    for (j=0;j<ny;j++){
      Tjy=gsl_cheb_evaluate_polynomial(j,y);
      for (k=0;k<nz;k++){
        Tkz=gsl_cheb_evaluate_polynomial(k,z);
        sum+=gsl_vector_get(c_ijk,index)*Tix*Tjy*Tkz;
        index+=1;
      }
    }
  }

  return sum ;

}


/* Evaluate the chebyshev polinomials for amplitude and phase at node T_j
 */
static int chebyshev_interpolation3d(
  double q,
  double lambda1,
  double lambda2,
  int nx, int ny, int nz,
  gsl_vector *cvec_amp,
  gsl_vector *cvec_phi,
  int nk_amp,
  int nk_phi,
  const double xyz_min[],
  const double xyz_max[],
  gsl_vector *interp_amp,   //return: A(T_j;q,lambda1,lambda2)    <=> p_j(q,lambda1,lambda2)
  gsl_vector *interp_phi)   //return: \Phi(T_j;q,lambda1,lambda2) <=> p_j(q,lambda1,lambda2)
{

  double sum = 0.0;
  int k = 0;
  int N=nx*ny*nz ;

  double xrescale = (q-0.5*(xyz_max[0]+xyz_min[0])) / (0.5*(xyz_max[0]-xyz_min[0]));
  double yrescale = (lambda1-0.5*(xyz_max[1]+xyz_min[1])) / (0.5*(xyz_max[1]-xyz_min[1]));
  double zrescale = (lambda2-0.5*(xyz_max[2]+xyz_min[2])) / (0.5*(xyz_max[2]-xyz_min[2]));

  /*-- Calculate interp_amp --*/
  for (k=0; k<nk_amp; k++) { // For each empirical node
    gsl_vector v = gsl_vector_subvector(cvec_amp, k*N, N).vector; // Pick out the coefficient matrix corresponding to the k-th node.
//     sum = chebev_3d(&v, nx, ny, nz, q, lambda1, lambda2, xyz_min, xyz_max) ;
    sum = gsl_cheb_eval_3d(&v, nx, ny, nz, xrescale, yrescale, zrescale) ;
    gsl_vector_set(interp_amp, k, sum); //write p_k(x,y,z)
  }

  /*-- Calculate interp_phi --*/
  for (k=0; k<nk_phi; k++) { // For each empirical node
    gsl_vector v = gsl_vector_subvector(cvec_phi, k*N, N).vector; // Pick out the coefficient matrix corresponding to the k-th node.
//     sum = chebev_3d(&v, nx, ny, nz, q, lambda1, lambda2, xyz_min, xyz_max) ;
    sum = gsl_cheb_eval_3d(&v, nx, ny, nz, xrescale, yrescale, zrescale) ;
    gsl_vector_set(interp_phi, k, sum); //write p_k(x,y,z)
  }

  return 0;
}


static int TEOBResumROMCore(
  REAL8TimeSeries **hPlus,
  REAL8TimeSeries **hCross,
  double phiRef, // orbital reference phase NOTE: unused
  double deltaT,
  double fLow,
  double distance,
  double inclination,
  double Mtot, // in Msol
  double eta,
  double lambda1, //in range 50 - 5000
  double lambda2  //in range 50 - 5000
)
{

  //NOTE: silly
  inclination = inclination + phiRef - phiRef ;

  /* Check output arrays */
  if(!hPlus || !hCross)
    XLAL_ERROR(XLAL_EFAULT);
  TEOBResumROMdataDS *romdata=&__lalsim_TEOBResumROMDS_data;
  if(*hPlus || *hCross)
  {
    XLALPrintError("(*hPlus) and (*hCross) are supposed to be NULL, but got %p and %p",(*hPlus),(*hCross));
    XLAL_ERROR(XLAL_EFAULT);
  }
  int retcode=0;

  REAL8TimeSeries *hp;
  REAL8TimeSeries *hc;

  /* Select ROM submodel */
  TEOBResumROMdataDS_submodel *submodel;
  submodel = romdata->sub1;

  double x = sqrt(1.0-4.0*eta) ;
  double q = (1-x)/(1+x);

  //Allocate space for the nodes
  gsl_vector *amp_at_nodes = gsl_vector_alloc(submodel->times->size);
  gsl_vector *phi_at_nodes = gsl_vector_alloc(submodel->times->size);

  double *amp_interp = calloc(Gntimes,sizeof(double));
  double *phi_interp = calloc(Gntimes,sizeof(double));
  double *freqs = calloc(Gntimes,sizeof(double));
  double *physical_times = calloc(Gntimes,sizeof(double));

  //calculate chebyshev interpolants (A(T_j),Phi(T_j))
  retcode = chebyshev_interpolation3d(q,lambda1,lambda2,
                                      Gnq, Gnlambda1, Gnlambda2,
                                      submodel->cvec_amp,submodel->cvec_phi,
                                      Gnamp,Gnphase,Gparams_min, Gparams_max,
                                      amp_at_nodes,phi_at_nodes);

  //calculate A(T_j) and Phi(T_j) at nodes
  double BjAmp_tn=0.0;
  double BjPhi_tn=0.0;
  int n,j;
  double c3 = LAL_C_SI*LAL_C_SI*LAL_C_SI ;
  double time_to_phys = LAL_G_SI*Mtot*LAL_MSUN_SI/c3 ;
  for (n=0;n<Gntimes;n++){
    BjAmp_tn=0.0 ;
    BjPhi_tn=0.0 ;
    for (j=0;j<Gnamp;j++){
      BjAmp_tn+=gsl_vector_get(amp_at_nodes,j)*gsl_matrix_get(submodel->Bamp,j,n);
    }
    for (j=0;j<Gnphase;j++){
      BjPhi_tn+=gsl_vector_get(phi_at_nodes,j)*gsl_matrix_get(submodel->Bphi,j,n);
    }
    //convert times in to seconds
    physical_times[n]=gsl_vector_get(submodel->times,n)*time_to_phys;
    amp_interp[n]=BjAmp_tn;
    phi_interp[n]=BjPhi_tn;
  }



  //Resampling A(t) and Phi(t) to arbitrary deltaT ---\n");
  gsl_interp_accel *acc_amp = gsl_interp_accel_alloc();
  gsl_interp_accel *acc_phi = gsl_interp_accel_alloc();
  gsl_interp_accel *acc_fot = gsl_interp_accel_alloc();
  gsl_spline *ampoft_spline = gsl_spline_alloc (gsl_interp_cspline, Gntimes);
  gsl_spline *phioft_spline = gsl_spline_alloc (gsl_interp_cspline, Gntimes);

  //initializing splines
  gsl_spline_init(ampoft_spline, physical_times, amp_interp, Gntimes);
  gsl_spline_init(phioft_spline, physical_times, phi_interp, Gntimes);

  double der ;
  //calculate frequencies at nodes (derivative of phi(t)/2pi)
  int i_end_mono = Gntimes;
  for (n=0;n<Gntimes;n++) {
    der = gsl_spline_eval_deriv (phioft_spline, physical_times[n], acc_phi);
    freqs[n] = 0.5*der/LAL_PI;//omegaoft(time_phys)/(2*np.pi)
    //determine up to where f is monotonically increasing
    if (n > 0 && i_end_mono == Gntimes) {
      if (freqs[n] < freqs[n-1]) i_end_mono = n ;
    }
  }

  //create t(f) spline
  double physical_times_end = physical_times[Gntimes-1];//save the last element before resizing
  //Resize freqs and physical_times to include only monotonically increasing values
  realloc(freqs,i_end_mono*sizeof(double));
  realloc(physical_times,i_end_mono*sizeof(double));
  //interpolate to get t(f)
  gsl_spline *toffreq_spline = gsl_spline_alloc (gsl_interp_cspline, i_end_mono);
  gsl_spline_init(toffreq_spline, freqs, physical_times, i_end_mono);

  //calculate parameters to resample with even spacing
  double tstart = gsl_spline_eval(toffreq_spline, fLow, acc_fot);
  int Ntimes_res = (int) ceil((physical_times_end-tstart)/deltaT);
  double *times_res = calloc(Ntimes_res,sizeof(double));
  double *amp_res = calloc(Ntimes_res,sizeof(double));
  double *phi_res = calloc(Ntimes_res,sizeof(double));
  double t=tstart;
  fprintf(stdout,"tstart=%.2f\n",tstart);
  //fprintf(stdout,"Ntimes=%d\n",Ntimes_res);

  //for scaling the amplitude
  double h22_to_h = 4.0*eta*sqrt(5.0/LAL_PI)/8.0;
  double amp_units = LAL_G_SI*Mtot*LAL_MSUN_SI/(LAL_C_SI*LAL_C_SI*distance) ;

  //Adjust for inclination angle [0,pi]
  double cosi = cos(inclination);
  double inc_plus = (1.0+cosi*cosi)/2.0;
  double inc_cross = cosi;


  //Generate h+(t) and hx(t)

  //XLALGPSAdd(&tC, -1 / deltaF);  /* coalesce at t=0 */
  LIGOTimeGPS tC = LIGOTIMEGPSZERO;
  XLALGPSAdd(&tC, tstart);
  //XLALGPSAdd(&tC, -1.0*j*deltaT);
  /* Allocate hplus and hcross */
  hp = XLALCreateREAL8TimeSeries("hplus: TD waveform", &tC, 0.0, deltaT, &lalStrainUnit, Ntimes_res);
  if (!hp) XLAL_ERROR(XLAL_EFUNC);
  memset(hp->data->data, 0, Ntimes_res * sizeof(REAL8));

  hc = XLALCreateREAL8TimeSeries("hcross: TD waveform", &tC, 0.0, deltaT, &lalStrainUnit, Ntimes_res);
  if (!hc) XLAL_ERROR(XLAL_EFUNC);
  memset(hc->data->data, 0, Ntimes_res * sizeof(REAL8));

  times_res[0] = t ;
  amp_res[0] = gsl_spline_eval(ampoft_spline, t, acc_amp)*amp_units*h22_to_h;
  double phi0 = gsl_spline_eval(phioft_spline, t, acc_phi);
  phi_res[0] = 0.0;
  hp->data->data[0] = inc_plus*amp_res[0]*cos(phi_res[0]);
  hc->data->data[0] = inc_cross*amp_res[0]*sin(phi_res[0]);
  t+=deltaT;
  for (n=1;n<Ntimes_res;n++) {
    times_res[n] = t;
    amp_res[n] = gsl_spline_eval(ampoft_spline, t, acc_amp)*amp_units*h22_to_h;
    //Zero the phase at the beginning (-phi0)
    phi_res[n] = gsl_spline_eval(phioft_spline, t, acc_phi)-phi0;

    hp->data->data[n] = inc_plus*amp_res[n]*cos(phi_res[n]);
    hc->data->data[n] = inc_cross*amp_res[n]*sin(phi_res[n]);

    t+=deltaT;
  }

  *hPlus = hp;
  *hCross = hc;

  gsl_spline_free (ampoft_spline);
  gsl_spline_free (phioft_spline);
  gsl_spline_free (toffreq_spline);
  gsl_interp_accel_free (acc_amp);
  gsl_interp_accel_free (acc_phi);
  gsl_interp_accel_free (acc_fot);

  gsl_vector_free(amp_at_nodes);
  gsl_vector_free(phi_at_nodes);

  free(amp_interp);
  free(phi_interp);
  free(freqs);
  free(physical_times);
  free(times_res);
  free(amp_res);
  free(phi_res);

  if (retcode==0){
    return(XLAL_SUCCESS);
  } else {
    return(retcode);
  }

}

/**
 * @addtogroup LALSimInspiralTEOBResumROM_c
 *
 * @{
 *
 * @name TEOBResum Reduced Order Model (Tidal effects)
 *
 * @author Jeroen Meidam, ... (Based on SEOBNRv2ROM code written by Michael Puerrer and John Veitch)
 *
 * @brief C code for TEOBResum reduced order model which includes tidal effects.
 * See ... for the basic approach.
 * Further details in ...
 *
 * This is a time domain model that approximates the time domain EOB model with tidal effects.
 *
 * The binary data files are available at https://github.com/benjaminlackey/cbcrom/tree/master/data.
 * Put the *.dat files into a location in your LAL_DATA_PATH.
 * They must have the names
 *  - TEOBResumROM_Amp_ciall.dat
 *  - TEOBResumROM_Phase_ciall.dat
 *  - TEOBResumROM_Bamp_matrix.dat
 *  - TEOBResumROM_Bphase_matrix.dat
 *  - TEOBResumROM_times.dat
 *
 * @note Approximant name is TEOBResum_ROM
 *
 * @note Parameter ranges:
 *   0.5 <= q <= 1.0
 *   50 <= lambda_i <= 5000
 *
 * @{
 */


//Function to test data reading
int XLALSimInspiralTEOBResumROM(
  REAL8TimeSeries **hPlus, /**< Output: Frequency-domain waveform h+ */
  REAL8TimeSeries **hCross, /**< Output: Frequency-domain waveform hx */
  REAL8 phiRef,                                 /**< Orbital phase at reference frequency*/
  REAL8 deltaT,                                 /**< Sampling frequency (Hz) */
  REAL8 fLow,                                   /**< Starting GW frequency (Hz) */
  REAL8 fRef,                                   /**< Reference frequency (Hz); 0 defaults to fLow */
  REAL8 distance,                               /**< Distance of source (m) */
  REAL8 inclination,                            /**< Inclination of source (rad) */
  REAL8 m1SI,                                   /**< Mass of companion 1 (kg) */
  REAL8 m2SI,                                   /**< Mass of companion 2 (kg) */
  REAL8 lambda1,                                /**< dimensionless tidal deformability of body 1 */
  REAL8 lambda2)                                /**< dimensionless tidal deformability of body 1 */
{
  /* Internally we need m1 > m2, so change around if this is not the case */
  if (m1SI < m2SI) {
    // Swap m1 and m2
    double m1temp = m1SI;
    double l1temp = lambda1;
    m1SI = m2SI;
    lambda1 = lambda2;
    m2SI = m1temp;
    lambda2 = l1temp;
  }

  /* Get masses in terms of solar mass */
  double mass1 = m1SI / LAL_MSUN_SI;
  double mass2 = m2SI / LAL_MSUN_SI;
  double Mtot = mass1+mass2;
  double eta = mass1 * mass2 / (Mtot*Mtot);    /* Symmetric mass-ratio */

  double x = sqrt(1.0-4.0*eta) ;
  double q = (1-x)/(1+x);

  // 'Nudge' q to allowed boundary values if close by (q could be e.g. 0.499999)
  if (q > Gparams_max[0])  nudge(&q, Gparams_max[0], 1e-6);
  if (q < Gparams_min[0])  nudge(&q, Gparams_min[0], 1e-6);

  /* Check that parameters are within the limits for this model */
  if ( q<Gparams_min[0] || q>Gparams_max[0] ) {
    fprintf(stderr,"Error: q not in range (%.2f not in [%.2f,%.2f])\n",q,Gparams_min[0],Gparams_max[0]);
    return(XLAL_EDOM);
  }
  if ( lambda1<Gparams_min[1] || lambda1>Gparams_max[1] ) {
    fprintf(stderr,"Error: lambda1 not in range (%.2f not in [%.2f,%.2f])\n",lambda1,Gparams_min[1],Gparams_max[1]);
    return(XLAL_EDOM);
  }
  if ( lambda2<Gparams_min[2] || lambda2>Gparams_max[2] ) {
    fprintf(stderr,"Error: lambda2 not in range (%.2f not in [%.2f,%.2f])\n",lambda2,Gparams_min[2],Gparams_max[2]);
    return(XLAL_EDOM);
  }

  if (fRef!=0.0) fprintf(stdout,"WARNING: fREf != 0.0 -> TEOBResum_ROM does not do anything with fRef. It will be evaluated from fLow.\n");

  // Load ROM data if not loaded already
  fprintf(stdout,"initializing with TEOBResumROM_Init_LALDATA()\n");
  #ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&TEOBResumROM_is_initialized, TEOBResumROM_Init_LALDATA);
  #else
  TEOBResumROM_Init_LALDATA();
  #endif

  if(!TEOBResumROM_IsSetup()) XLAL_ERROR(XLAL_EFAILED,"Error setting up TEOBResumROM data - check your $LAL_DATA_PATH\n");

  int retcode = TEOBResumROMCore(hPlus,hCross, phiRef, deltaT, fLow, distance, inclination, Mtot, eta, lambda1, lambda2);

  return(retcode);
}



/** Setup TEOBResum_ROM model using data files installed in $LAL_DATA_PATH
 */
static void TEOBResumROM_Init_LALDATA(void)
{
  if (TEOBResumROM_IsSetup()) return;

  // If we find one ROM datafile in a directory listed in LAL_DATA_PATH,
  // then we expect the remaining datafiles to also be there.
  char datafile[] = "TEOBResumROM_Phase_ciall.dat";

  char *path = XLALFileResolvePathLong(datafile, PKG_DATA_DIR);
  if (path==NULL)
    XLAL_ERROR_VOID(XLAL_EIO, "Unable to resolve data file %s in $LAL_DATA_PATH\n", datafile);
  char *dir = dirname(path);
  int ret = TEOBResumROM_Init(dir);
  XLALFree(path);

  if(ret!=XLAL_SUCCESS)
    XLAL_ERROR_VOID(XLAL_FAILURE, "Unable to find TEOBResumROM data files in $LAL_DATA_PATH\n");
}
