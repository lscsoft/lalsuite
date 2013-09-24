/*
 Compile with:

 gcc -O2 -o BayesLine BayesLine.c -lgsl -lm

 Usage:

 ./BayesLine --help

 */

/*********************************************************************************/
/*                                                                               */
/*                                                                               */
/*     BayesLine v1.0. Written by Neil J. Cornish between 5/20/13 - 6/5/13       */
/*                                                                               */
/*     NJC thanks Tyson Littenberg and Will Farr for comments and suggestions    */
/*                                                                               */
/*                                                                               */
/*                                                                               */
/*     BayesLine fits the LIGO/Virgo power spectra using a model made up         */
/*     of N Lorentzian lines (described by central frequency f, quality          */
/*     factor Q and amplitude A) and cubic spline with M control points.         */
/*     The number of terms in each model, N, M, are free to vary via RJMCMC      */
/*     updates. The code initializes the models in a non-Markovian fashion,      */
/*     then refines the models with a full Markovian RJMCMC subroutine. This     */
/*     subroutine (LorentzMCMC) can be called by other codes to update the       */
/*     spectral fit to the residuals (data - signal model). Doing this is        */
/*     highly recommended, as it ensures that the spectral model is not          */
/*     eating any of the gravitational wave signal. Since the code is            */
/*     transdimensional (and very free in its movement between dimensions)       */
/*     it will not "over-fit" the spectral model.                                */
/*                                                                               */
/*                                                                               */
/*********************************************************************************/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>

#include <gsl/gsl_sort.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_statistics.h>

#define SAmin  1.0e-48
#define SAmax  1.0e-30
#define Qmin  1.0e2  // allowing Q's below ~ 1e2 starts a fight with the spline model as lines get too fat
#define Qmax 1.0e8
#define Amin  1.0e-44
#define Amax  1.0e-30
#define lAwidth 2.3    // 2.3 is one decade
#define zeta 1.0
#define kappa 0.8
//#define fgrid 4.0   // the stencil separation in Hz for the spline model. Going below 2 Hz is dangerous - will fight with line model
// as we are getting down to the line width of the broadest lines

typedef struct
{
  int gnuplot;
  int verbose;
  int zerologL;

  char ifile[64];
  char ofile[64];

}bayesline_opt;

typedef struct
{
  int n;
  int size;

  int *larray;

  double *Q;
  double *A;
  double *f;

}lorentzianParams;

typedef struct
{
  int tmax;
  int ncut;
  int nmin;
  int tfull;
  int sgmts;

  double df;
  double fny;
  double Tobs;
  double fmin;
  double fmax;
  double flow;
  double fgrid;
  double fstep;
  double fhigh;
  double cadence;

}dataParams;

typedef struct
{
  int n;
  double *points;
  double *data;

}splineParams;

typedef struct
{
  int n;
  double *f;
  double *pow;

  double *Sn;
  double *Sbase;
  double *Sline;

}psdParams;

void KILL(char* Message);

int opt_parse(bayesline_opt *opts, int argc, char **argv );
int bayesline_usage( const char *program );

double loglike           (double *respow, double *Snf, int ncut);
double loglike_fit_spline(double *respow, double *Snf, int ncut);

double loglike_pm        (double *respow, double *Sn, double *Snx, int ilow, int ihigh);
double loglike_single    (double *respow, double *Sn, double *Snx, int ilowx, int ihighx, int ilowy, int ihighy);

double sample(double *fprop, double pmax, dataParams *data, gsl_rng *r);
double lprop(double f, double *fprop, dataParams *data);

void full_spectrum_single(double *Sn, double *Snx, double *Sbasex, double *sfreq, dataParams *data, lorentzianParams *line_x, lorentzianParams *line_y, int ii, int *ilowx, int *ihighx, int *ilowy, int *ihighy);
void full_spectrum_add_or_subtract(double *Snew, double *Sold, double *Sbase, double *sfreq, dataParams *data, lorentzianParams *lines, int ii, int *ilow, int *ihigh, int flag);
void full_spectrum_spline(double *Sline, double *Sbase, double *sfreq, dataParams *data, lorentzianParams *lines);

void spectrum_spline(double *Sn, double *Sbase, double *sfreq, dataParams *data, lorentzianParams *lines, splineParams *spline);

void SpecFitSpline    (bayesline_opt opts, int steps, double *freq, double *power, splineParams *spline, double *Snf, int ncut, gsl_rng *r);
void LorentzSplineFit (bayesline_opt opts, int steps, dataParams *data, lorentzianParams *lines, splineParams *spline, double *sfreq, double *spow, gsl_rng *r);
void LorentzSplineMCMC(bayesline_opt opts, int steps, dataParams *data, double *freq, double *power, lorentzianParams *lines, splineParams *spline, splineParams *spline_x, double *Snf, int focus, double *dan, gsl_rng *r);

void CubicSplineGSL(int N, double *x, double *y, int Nint, double *xint, double *yint);

void create_dataParams(dataParams *data, double *f, int n);

void create_lorentzianParams(lorentzianParams *lines, int size);
void destroy_lorentzianParams(lorentzianParams *lines);

void create_splineParams(splineParams *spline, int size);
void destroy_splineParams(splineParams *spline);

void create_psdParams(psdParams *psd, int size);
void destroy_psdParams(psdParams *psd);

void gnuplot_full_spectrum(int n, double *f, double *p, double *Sn, FILE *pipe);

int main(int argc, char *argv[])
{
  int n;
  int imin, imax;
  int i, j, k;
  int jj, kk;
  double *fa, *Sna;
  double dan;
  double *spow, *sfreq, *Sn, *Sbase;
  int sm;
  double mdn;
  double spread, deltafmax;
  double sdatRe,sdatIm;
  double *freq;
  double *power;
  double *wndw;
  double *Snf;
  double x, *y, z;
  double fsq, frs;
  int nspline, smodn;
  double *xint,*yint;

  dataParams *data = malloc( (size_t)sizeof(dataParams) );

  lorentzianParams *lines_x    = malloc( (size_t)sizeof(lorentzianParams) );
  lorentzianParams *lines_full = malloc( (size_t)sizeof(lorentzianParams) );

  splineParams *spline   = malloc( (size_t)sizeof(splineParams) );
  splineParams *spline_x = malloc( (size_t)sizeof(splineParams) );

  char filename[100];
  char burnrows[1000];
  FILE *pipe;
  FILE *infile;
	FILE *outfile;

  /* Parse command line */
  bayesline_opt opts;
  opt_parse( &opts, argc, argv );

  fprintf(stdout,"\n============ BayesLine ==============\n");
  fprintf(stdout," Bayesian PSD modeling of %s:\n",opts.ifile);
  fprintf(stdout," Final SSD model will be output to %s:\n",opts.ofile);

  fprintf(stdout," gnuplot status updates are ");
  if(opts.gnuplot) fprintf(stdout,"ENABLED\n");
  else             fprintf(stdout,"DISABLED\n");

  fprintf(stdout," Verbose settings are ");
  if(opts.verbose) fprintf(stdout,"MAXIMAL\n");
  else             fprintf(stdout,"MINIMAL \n");

  fprintf(stdout," likelihood functions are ");
  if(opts.zerologL) fprintf(stdout,"CONSTANT\n");
  else              fprintf(stdout,"COMPUTED\n");

  fprintf(stdout,"\n");

  /* set up GSL random number generator */
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *r = gsl_rng_alloc (T);
  gsl_rng_env_setup();

  // This section reads in data
  sprintf(filename,"%s",opts.ifile);
	infile = fopen(filename,"r");

	n = 0;
	while(!feof(infile))
	{
    fgets(burnrows,1000,infile);
    n++;
  }
  n -= 1;
  rewind(infile);

  freq  = malloc((size_t)(sizeof(double)*(n)));
  power = malloc((size_t)(sizeof(double)*(n)));
  Snf   = malloc((size_t)(sizeof(double)*(n)));

  for(i=0; i<n; i++)
	{
    fscanf(infile,"%lg%lg%lg\n", &freq[i], &sdatRe, &sdatIm);
    power[i] = 2.0*(sdatRe*sdatRe+sdatIm*sdatIm);
  }
  fclose(infile);


  // storage for data meta parameters
  create_dataParams(data,freq,n);

  // storage for line model for a segment. These get updated and passed back from the fitting routine
  create_lorentzianParams(lines_x,data->tmax);

  // storage for master line model
  create_lorentzianParams(lines_full,data->tfull);

  // start with a single line. This number gets updated and passed back
  lines_x->n = 1;

  //start with ~4 Hz resolution for the spline
  nspline = (int)(data->fstep/data->fgrid)+2;
  smodn   = nspline*data->sgmts;

  fa      = malloc((size_t)(sizeof(double)*(smodn+1)));
  Sna     = malloc((size_t)(sizeof(double)*(smodn+1)));
  xint    = malloc((size_t)(sizeof(double)*(nspline)));
  yint    = malloc((size_t)(sizeof(double)*(nspline)));

  create_splineParams(spline,nspline);

  kk = 0;
  jj = 0;

  if(opts.gnuplot)
  {
    pipe = popen("gnuplot -persist","w");
    fprintf(pipe, "set logscale y\n");
  }

  /******************************************************************************/
  /*                                                                            */
  /*  Rapid non-Markovian fit over small bandwidth blocks of data               */
  /*                                                                            */
  /******************************************************************************/

  // loop over the frequency segments
  data->flow = data->fmin;
  printf("\nLorentzSplineFit() loop over frequency segments...\n");
  do
  {
    data->fhigh = data->flow + data->fstep;
    imax = (int)(floor(data->fhigh * data->Tobs));
    imin = (int)(floor(data->flow  * data->Tobs));
    data->ncut = imax-imin;

    lines_x->n = 1;

    y     = malloc((size_t)(sizeof(double)*(data->ncut)));
    Sn    = malloc((size_t)(sizeof(double)*(data->ncut)));
    Sbase = malloc((size_t)(sizeof(double)*(data->ncut)));
    spow  = malloc((size_t)(sizeof(double)*(data->ncut)));
    sfreq = malloc((size_t)(sizeof(double)*(data->ncut)));
    wndw  = malloc((size_t)(sizeof(double)*(data->ncut)));

    sm = data->ncut/2;

    for(i=0; i<data->ncut; i++)
    {
      j = i+imin-data->nmin;

      spow[i]  = power[j];
      wndw[i]  = power[j];
      sfreq[i] = freq[j];
    }
    // The 1.47 factor converts from median to mean
    // We don't care about the even-odd issue since we have so many terms
    gsl_sort(wndw,1,data->ncut);
    mdn = log( 1.47*( gsl_stats_quantile_from_sorted_data(wndw,1,data->ncut,0.5) ) );

    for(j=0; j<nspline; j++)
    {
      spline->points[j] = data->flow+(data->fhigh-data->flow)*(double)(j)/(double)(nspline-1);
      spline->data[j]   = mdn;
    }

    LorentzSplineFit(opts, 40000, data, lines_x, spline, sfreq, spow, r);

    //Interpoloate {spoints,sdata} ==> {sfreq,y}
    CubicSplineGSL(spline->n,spline->points,spline->data,data->ncut,sfreq,y);

    for(i=0; i< data->ncut; i++)
    {
      Sbase[i] = exp(y[i]);
      Sn[i]    = Sbase[i];
      fsq      = sfreq[i]*sfreq[i];

      for(j=0; j< lines_x->n; j++)
      {
        spread = (1.0e-2*lines_x->Q[j]);
        if(spread < 50.0) spread = 50.0;  // maximum half-width is f_resonance/100
        deltafmax = lines_x->f[j]/spread;
        frs = lines_x->f[j]*lines_x->f[j];
        x = fabs(lines_x->f[j]-sfreq[i]);
        z = 1.0;
        if(x > deltafmax) z = exp(-(x-deltafmax)/deltafmax);
        Sn[i] += z*lines_x->A[j]*frs*frs/(frs*fsq+lines_x->Q[j]*lines_x->Q[j]*(fsq-frs)*(fsq-frs));
      }
    }

    if(opts.verbose)
    {
      outfile = fopen("model_final.dat","w");
      for(i=0; i<data->ncut; i++)
      {
        fprintf(outfile, "%f %e %e %e %e\n", sfreq[i], spow[i], Sn[i], Sbase[i], spow[i]/Sn[i]);
      }
      fclose(outfile);
    }
    if(opts.gnuplot)
    {
      fprintf(pipe, "plot 'model_final.dat' using 1:2 title 'data' with lines, 'model_final.dat' using 1:3 title 'model' with lines lt 3, 'model_final.dat' using 1:4 title 'base' with lines lt 2\n");
      fflush(pipe);
    }
    free(y);
    free(Sbase);
    free(Sn);
    free(spow);
    free(sfreq);
    free(wndw);

    for(j=0; j<nspline; j++) xint[j] = data->flow+(data->fhigh-data->flow)*(double)(j)/(double)(nspline-1);
    CubicSplineGSL(nspline,spline->points,spline->data,nspline,xint,yint);

    for(j=0; j<nspline; j++)
    {
      fa[jj]  = xint[j];
      Sna[jj] = yint[j];
      jj++;
    }
    if(opts.verbose)fprintf(stdout,"New lines = %d; Total lines = %d;  f = %f Hz, Sn(f)~%e\n", lines_x->n, kk+lines_x->n, fa[jj-1], Sna[jj-1]);

    for(j=0; j<lines_x->n; j++)
    {
      lines_full->f[kk] = lines_x->f[j];
      lines_full->Q[kk] = lines_x->Q[j];
      lines_full->A[kk] = lines_x->A[j];
      kk++;
    }

    data->flow += data->fstep;
  }while(data->fhigh+data->fstep < data->fmax);
  fprintf(stdout, "  Line model used %i Lorentzians\n",kk);

  // this is the initial fit to the smooth part of the spectrum
  outfile = fopen("spec.dat","w");
  for(i=0; i< jj; i++)
  {
    fprintf(outfile, "%f %e\n", fa[i], Sna[i]);
  }
  fclose(outfile);

  destroy_lorentzianParams(lines_x);
  destroy_splineParams(spline);
  free(xint);
  free(yint);


  /******************************************************************************/
  /*                                                                            */
  /*  Full-data spline fit (holding lines fixed from search phase)              */
  /*                                                                            */
  /******************************************************************************/

  // current number of terms in the Lorentzian model
  lines_full->n = kk-1;

  // maximum number of terms in the Lorentzian model
  data->tmax = 4*lines_full->n;

  if(data->tmax    < 20) data->tmax    = 20;
  if(lines_full->n < 1 ) lines_full->n = 1;

  //start with ~4 Hz resolution for the spline
  nspline = (int)((data->fhigh-data->fmin)/data->fgrid)+2;

  spline = malloc( (size_t)sizeof(splineParams) );
  create_splineParams(spline,nspline);
  create_splineParams(spline_x,nspline);

  //initialize spline grid & PSD values
  k=0;
  for(j=0; j<nspline; j++)
  {
    x = data->fmin+(data->fhigh-data->fmin)*(double)(j)/(double)(nspline-1);

    if     (x <= fa[0])    z = Sna[0];
    else if(x >= fa[jj-1]) z = Sna[jj-1];
    else
    {
      i=k-10;
      mdn = 0.0;
      do
      {
        if(i>=0) mdn = fa[i];
        i++;
      } while(x > mdn);

      k = i;
      z = Sna[i];
    }

    spline->points[j] = x;
    spline->data[j]   = z;
  }

  // produce an initial spline fit to the smooth part of the spectrum
  SpecFitSpline(opts, 100000, fa, Sna, spline, Snf, jj, r);

  outfile = fopen("spec_fit.dat","w");
  for(i=0; i< jj; i++) fprintf(outfile, "%f %e\n", fa[i], Snf[i]);
  fclose(outfile);

  outfile = fopen("spec_end.dat","w");
  for(j=0; j< nspline; j++) fprintf(outfile, "%f %e\n", spline->points[j], spline->data[j]);
  fclose(outfile);

  if(opts.gnuplot)
  {
    fprintf(pipe, "plot 'spec.dat' using 1:(exp($2)) title 'raw fit' with lines, 'spec_fit.dat' using 1:(exp($2)) title 'smooth fit' with lines lt 3\n");
    fflush(pipe);
  }


  /******************************************************************************/
  /*                                                                            */
  /*  Full spline/line MCMC stage                                               */
  /*                                                                            */
  /******************************************************************************/

  /*
   The spow and sfreq arrays hold the frequency and power of the full section of data to be whitened
   The model parameters are the number of Lorentzians, nx, and their frequency ff, amplitude AA and
   quality factor QQ, as well as the number of terms in the power law fit to the smooth part of the
   spectrum, segs, and their amplitudes at 100 Hz, SA, and their spectral slopes SP. The maximum number
   of terms in the Lorentzian model is tmax (around 200-300 should cover most spectra), and the
   maximum number of terms in the power law fit is segmax (aound 12-20 should be enough). When called
   from another MCMC code, the LorentzMCMC code needs to be passed the most recent values for the
   noise model. Arrays to hold them have to be declared somewhere in advance. The updated values are
   passed back. The first entry in the call to LoretnzMCMC are the number of iterations. It probably
   makes sense to do 100 or so each time it is called.
   */

  // Re-set dataParams structure to use full dataset
  data->flow  = data->fmin;
  imax  = (int)(floor(data->fhigh*data->Tobs));
  imin  = (int)(floor(data->flow*data->Tobs));
  data->ncut  = imax-imin;

  spow  = malloc((size_t)(sizeof(double)*(data->ncut)));
  sfreq = malloc((size_t)(sizeof(double)*(data->ncut)));

  for(i=0; i< data->ncut; i++)
  {
    j = i + imin - data->nmin;
    spow[i] = power[j];
    sfreq[i] = freq[j];
  }

  //copy spineParams structure spline to spline_x
  /*
   spline_x will be varied (and updated) by LorentzSplineMCMC.
   spline contains the original grid used in the RJ proposals
   */
  for(j=0; j<nspline; j++)
  {
    spline_x->points[j] = spline->points[j];
    spline_x->data[j]   = spline->data[j];
  }

  LorentzSplineMCMC(opts, 20000, data, sfreq, spow, lines_full, spline, spline_x, Snf, 0, &dan, r);
  if(opts.gnuplot) gnuplot_full_spectrum(data->ncut, sfreq, spow, Snf, pipe);

  //alternate between targeting outliers and MCMCing full solution until outliers are gone or max iterations are reached
  j = 0;
  do
  {
    LorentzSplineMCMC(opts, 5000, data, sfreq, spow, lines_full, spline, spline_x, Snf, 1, &dan, r);
    if(opts.gnuplot) gnuplot_full_spectrum(data->ncut, sfreq, spow, Snf, pipe);

    LorentzSplineMCMC(opts, 5000, data, sfreq, spow, lines_full, spline, spline_x, Snf, 0, &dan, r);
    if(opts.gnuplot) gnuplot_full_spectrum(data->ncut, sfreq, spow, Snf, pipe);

    j++;

  } while (dan > 25.0 && j < 8);

  //final full-model MCMC
  LorentzSplineMCMC(opts, 20000, data, sfreq, spow, lines_full, spline, spline_x, Snf, 0, &dan, r);

  if(opts.gnuplot) gnuplot_full_spectrum(data->ncut, sfreq, spow, Snf, pipe);

  /******************************************************************************/
  /*                                                                            */
  /*  Output PSD in format for LALInference                                     */
  /*                                                                            */
  /******************************************************************************/

  sprintf(filename,"%s",opts.ofile);
  outfile = fopen(filename,"w");

  //at frequencies below fmin output SAmax (killing lower frequencies in any future inner products)
  imax = (int)(floor(data->fhigh * data->Tobs));
  imin = (int)(floor(data->flow  * data->Tobs));

  for(i=0; i<2*imax; i++)
  {
  if(i>=imin && i<imax)
    fprintf(outfile,"%.12g %.12g\n",(double)i/data->Tobs, sqrt(Snf[i-imin]/data->Tobs) );
  else
    fprintf(outfile,"%.12g %.12g\n",(double)i/data->Tobs, sqrt(SAmax/data->Tobs)       );
  }
  fclose(outfile);

  free(fa);
  free(Sna);
  free(spow);
  free(sfreq);
  free(freq);
  free(power);
  free(Snf);

  return 0;
}

void KILL(char* Message)
{
  fprintf(stderr,"\a\n");
  fprintf(stderr,"%s",Message);
  fprintf(stderr,"Terminating the program.\n\n\n");
  exit(1);

  return;
}

double sample(double *fprop, double pmax, dataParams *data, gsl_rng *r)
{
  int i;
  double f, alpha;

  int ssize    = data->ncut;
  double flow  = data->flow;
  double fhigh = data->fhigh;
  double Tobs  = data->Tobs;

  do
  {
    f = flow+gsl_rng_uniform(r)*(fhigh-flow);
    i = (int)(floor((f-flow)*Tobs));
    if(i < 0) i=0;
    if(i > ssize-1) i=ssize-1;
    alpha = pmax*gsl_rng_uniform(r);
    //printf("%d %f %f %f\n", i, f, alpha, fprop[i]);
  } while(alpha > fprop[i]);

  return(f);
}


double lprop(double f, double *fprop, dataParams *data)
{
  int i;

  int ssize    = data->ncut;
  double flow  = data->flow;
  double Tobs  = data->Tobs;

  i = (int)((f-flow)*Tobs);
  if(i < 0) i=0;
  else if(i > ssize-1) i=ssize-1;

  return(log(fprop[i]));
}

double loglike_fit_spline(double *respow, double *Snf, int ncut)
{
  double lgl, x;
  int i;

  lgl = 0.0;
  for(i=0; i< ncut; i++)
  {
    x = (respow[i]-Snf[i])*(respow[i]-Snf[i])/0.1;
    lgl -= (x);
  }

  return(lgl);
}



double loglike(double *respow, double *Snf, int ncut)
{
  double lgl, x;
  int i;

  // leavimng out the log(2Pi) terms since they cancel in Hastings ratio
  lgl = 0.0;
  for(i=0; i< ncut; i++)
  {
    x = respow[i]/Snf[i];
    lgl -= (x+log(Snf[i]));
  }

  return(lgl);
}

void LorentzSplineFit(bayesline_opt opts, int steps, dataParams *data, lorentzianParams *lines_x, splineParams *spline, double *sfreq, double *spow, gsl_rng *r)
{
  int spass, passes, psmax;
  double logLx, logLy, logH;
  int i, j, k, ki, ii, jj, mc;
  int check;
  double alpha, heat;
  double SAmaxx, SAminn, lSAmax, lSAmin;
  double lQmin, lQmax, QQ;
  double Aminn, Amaxx, lAmin, lAmax;
  int ac0, ac1, cnt;
  int cc0, cc1, cc2;
  double *Sn, *Sbase;
  double e1, e2, e3, e4;
  double x2, x3, x4;
  double s1, s2, s3, s4;
  int typ;
  double xsm, pmax, y, z;
  double mdn, baseav;
  double logpx, logpy, x, beta;
  double *fprop;

  int ncut     = data->ncut;
  int tmax     = data->tmax;
  double flow  = data->flow;
  double fhigh = data->fhigh;

  lorentzianParams *lines_y = malloc(sizeof(lorentzianParams));
  create_lorentzianParams(lines_y,tmax);

  splineParams *spline_y = malloc(sizeof(splineParams));
  create_splineParams(spline_y,spline->n);

  int    nspline   = spline->n;
  double *spoints  = spline->points;
  double *sdatax   = spline->data;
  double *sdatay   = spline_y->data;
  double *spointsy = spline_y->points;

  Sn     = malloc((size_t)(sizeof(double)*(ncut)));
  Sbase  = malloc((size_t)(sizeof(double)*(ncut)));
  fprop  = malloc((size_t)(sizeof(double)*(ncut)));


  mdn    = exp(sdatax[0]);

  SAmaxx = 1.0e2*mdn;
  SAminn = 1.0e-2*mdn;
  Aminn = mdn;
  Amaxx = mdn*1.0e6;
  if(SAmaxx > SAmax) SAmaxx = SAmax;
  if(SAminn < SAmin) SAminn = SAmin;
  if(Amaxx  > Amax)  Amaxx  = Amax;
  if(Aminn  < Amin)  Aminn  = Amin;
  lQmin  = log(Qmin);
  lQmax  = log(Qmax);
  lAmin  = log(Aminn);
  lAmax  = log(Amaxx);
  lSAmin = log(SAminn);
  lSAmax = log(SAmaxx);

  // this is the fractional error estimate on the noise level
  s1 = 1.0/sqrt((double)(ncut));

  s2 = 0.01;
  s3 = 0.5;
  s4 = 0.5;

  // set up proposal for frequency jumps
  xsm =0.0;
  for(i=0; i< ncut; i++)
	{
    x = spow[i]/mdn;
    if(x <  100.0) x = 1.0;
    if(x >= 100.0) x = 100.0;
    fprop[i] = x;
    xsm += x;
  }

  pmax = -1.0;
  for(i=0; i< ncut; i++)
	{
    fprop[i] /= xsm;
    if(fprop[i] > pmax) pmax = fprop[i];
  }

  for(i=0; i< lines_x->n; i++)
  {
    lines_x->larray[i] = i;
    lines_y->larray[i] = i;
  }

  alpha = gsl_rng_uniform(r);

  for(i=0; i<lines_x->n; i++)
  {
    lines_x->Q[i] = exp(lQmin+(lQmax-lQmin)*gsl_rng_uniform(r));
    lines_x->f[i] = sample(fprop, pmax, data, r);
    lines_x->A[i] = mdn;
  }

  spectrum_spline(Sn, Sbase, sfreq, data, lines_x, spline);

  baseav = 0.0;
  for(i=0; i< ncut; i++)
	{
    baseav += log(Sbase[i]);
  }
  baseav /= (double)(ncut);


  if(!opts.zerologL) logLx = loglike(spow, Sn, ncut);
  else               logLx = 1.0;

  cnt = 0;
  passes = 0;
  spass = steps/4;

  psmax = 2;

  do
  {
    ac0 = 0;
    ac1 = 0;
    cc0 = 1;
    cc1 = 1;
    cc2 = 1;

    for(mc=0; mc < spass; mc++)
    {
      heat = 1.0;
      if(mc < spass/2) heat = pow(10.0,1.0*(double)(spass-2.0*mc)/(double)(spass));

      //copy over current state
      for(i=0; i< lines_x->n; i++)
      {
        lines_y->larray[i] = lines_x->larray[i];
        lines_y->Q[i] = lines_x->Q[i];
        lines_y->f[i] = lines_x->f[i];
        lines_y->A[i] = lines_x->A[i];
      }
      for(i=0; i< nspline; i++)
      {
        sdatay[i]   = sdatax[i];
        spointsy[i] = spoints[i];
      }
      alpha = gsl_rng_uniform(r);

      if(alpha > 0.5)  // try a transdimensional move
      {

        alpha = gsl_rng_uniform(r);
        if(alpha < 0.5)  // try and add a new term
        {
          lines_y->n = lines_x->n+1;
          typ = 2;
        }
        else // try and remove term
        {
          lines_y->n = lines_x->n-1;
          typ = 3;
        }

        check = 0;
        if(lines_y->n < 0 || lines_y->n > tmax) check = 1;


        if(check == 0)
        {

          lines_y->Q[i] = lines_x->Q[i];
          lines_y->f[i] = lines_x->f[i];
          lines_y->A[i] = lines_x->A[i];

          if(lines_y->n < lines_x->n)
          {
            i=(int)(gsl_rng_uniform(r)*(double)(lines_x->n)); // pick a term to try and kill
            k = 0;
            for(j=0; j< lines_x->n; j++)
            {
              if(j != i)
              {
                lines_y->larray[k] = lines_x->larray[j];
                jj = lines_x->larray[j];
                lines_y->A[jj] = lines_x->A[jj];
                lines_y->Q[jj] = lines_x->Q[jj];
                lines_y->f[jj] = lines_x->f[jj];
                k++;
              }
              if(j == i) ki = lines_x->larray[j];  // take note of who's demise is being proposed
            }

            logpx = lprop(lines_x->f[ki], fprop, data);
            logpy = -log((double)(ncut));    // corresponds to uniform density - just as likely to kill from alines_y->nwhere

            y = log(lines_x->A[ki])-baseav;
            if(y < 0.0)
            {
              z = 0.0;
            }
            else
            {
              z = 2.0*kappa*gsl_ran_gaussian_pdf(y, lAwidth);
            }
            logpx += log((1.0-kappa)/(lAmax-lAmin) + z);
            logpy += -log((lAmax-lAmin));
          }

          if(lines_y->n > lines_x->n)
          {
            // copy over current state
            for(j=0; j< lines_x->n; j++)
            {
              k = lines_x->larray[j];
              lines_y->A[k] = lines_x->A[k];
              lines_y->f[k] = lines_x->f[k];
              lines_y->Q[k] = lines_x->Q[k];
              lines_y->larray[j] = lines_x->larray[j];
            }

            // find a label that isn't in use

            i = -1;
            do
            {
              i++;
              k = 0;
              for(j=0; j< lines_x->n; j++)
              {
                if(i==lines_x->larray[j]) k = 1;
              }
            } while(k == 1);

            lines_y->larray[lines_x->n] = i;
            ii=i;

            // draw new term
            //lines_y->A[i] = exp(lAmin+(lAmax-lAmin)*gsl_rng_uniform(r));
            lines_y->Q[i] = exp(lQmin+(lQmax-lQmin)*gsl_rng_uniform(r));
            lines_y->f[i] = sample(fprop, pmax, data, r);
            logpy = lprop(lines_y->f[i], fprop, data);
            logpx = -log((double)(ncut));    // corresponds to uniform density - just as likely to kill from alines_y->nwhere

            alpha = gsl_rng_uniform(r);
            if(alpha < kappa)
            {
              y = fabs(gsl_ran_gaussian(r, lAwidth));
              lines_y->A[i] = exp(baseav+y);
            }
            else
            {
              lines_y->A[i] = exp(lAmin+(lAmax-lAmin)*gsl_rng_uniform(r));
              y = log(lines_y->A[i])-baseav;
            }
            if(y < 0.0)
            {
              z = 0.0;
            }
            else
            {
              z = 2.0*kappa*gsl_ran_gaussian_pdf(y, lAwidth);
            }
            logpy += log((1.0-kappa)/(lAmax-lAmin) + z);
            logpx += -log((lAmax-lAmin));

            if(lines_y->A[i] > Amaxx) check = 1;
            if(lines_y->A[i] < Aminn) check = 1;
          }
        }
      }
      else  // regular MCMC update
      {

        lines_y->n = lines_x->n;

        e1 = s1;

        for(i=0; i< nspline; i++) sdatay[i] = sdatax[i] + gsl_ran_gaussian(r, e1);

        if(lines_y->n > 0)
        {

          //pick a term to update
          jj=(int)(gsl_rng_uniform(r)*(double)(lines_x->n));
          //find label of who is geting updated
          ii = lines_x->larray[jj];

          // copy over current state
          for(i=0; i< lines_x->n; i++)
          {
            k = lines_x->larray[i];
            lines_y->f[k] = lines_x->f[k];
            lines_y->A[k] = lines_x->A[k];
            lines_y->Q[k] = lines_x->Q[k];
            lines_y->larray[i] = lines_x->larray[i];
          }


          alpha = gsl_rng_uniform(r);
          if(alpha > 0.8)  // 0.8
          {
            lines_y->Q[ii] = exp(lQmin+(lQmax-lQmin)*gsl_rng_uniform(r));
            lines_y->f[ii] = sample(fprop, pmax, data, r);
            logpy = lprop(lines_y->f[ii], fprop, data);
            logpx = lprop(lines_x->f[ii], fprop, data);
            typ = 0;


            alpha = gsl_rng_uniform(r);
            if(alpha < kappa)
            {
              y = fabs(gsl_ran_gaussian(r, lAwidth));
              lines_y->A[ii] = exp(baseav+y);
            }
            else
            {
              lines_y->A[ii] = exp(lAmin+(lAmax-lAmin)*gsl_rng_uniform(r));
              y = log(lines_y->A[ii])-baseav;
            }
            if(y < 0.0)
            {
              z = 0.0;
            }
            else
            {
              z = 2.0*kappa*gsl_ran_gaussian_pdf(y, lAwidth);
            }
            logpy += log((1.0-kappa)/(lAmax-lAmin) + z);


            y = log(lines_x->A[ii])-baseav;
            if(y < 0.0)
            {
              z = 0.0;
            }
            else
            {
              z = 2.0*kappa*gsl_ran_gaussian_pdf(y, lAwidth);
            }
            logpx += log((1.0-kappa)/(lAmax-lAmin) + z);
          }
          else
          {
            typ = 1;
            alpha = gsl_rng_uniform(r);

            if     (alpha > 0.9) beta = 1.0e+1;
            else if(alpha > 0.6) beta = 1.0e+0;
            else if(alpha > 0.3) beta = 1.0e-1;
            else                 beta = 1.0e-2;

            e2 = beta*s2;
            e3 = beta*s3;
            e4 = beta*s4;
            x2 = gsl_ran_gaussian(r, e2);
            x4 = gsl_ran_gaussian(r, e4);
            x3 = gsl_ran_gaussian(r, e3);

            lines_y->A[ii] = lines_x->A[ii]*exp(x3);
            lines_y->Q[ii] = lines_x->Q[ii]*exp(x4);
            lines_y->f[ii] = lines_x->f[ii]+x2;

            logpx = logpy = 0.0;
          }

          check =0;


          for(i=0; i<nspline; i++) if(sdatay[i] > lSAmax) check = 1;
          for(i=0; i<nspline; i++) if(sdatay[i] < lSAmin) check = 1;
          if(lines_y->A[ii] > Amaxx) check = 1;
          if(lines_y->A[ii] < Aminn) check = 1;
          if(lines_y->f[ii] < flow)  check = 1;
          if(lines_y->f[ii] > fhigh) check = 1;
          if(lines_y->Q[ii] < Qmin)  check = 1;
          if(lines_y->Q[ii] > Qmax)  check = 1;

        }
      }

      if(check == 0)
      {

        if(typ == 0) cc0++;
        if(typ == 1) cc1++;
        if(typ == 2) cc2++;

        if(!opts.zerologL)
        {
          spectrum_spline(Sn, Sbase, sfreq, data, lines_y, spline_y);
          logLy = loglike(spow, Sn, ncut);
        }
        else
        {
          logLy = 1.0;
        }

        // prior on line number e(-zeta * n).  (this is a prior, not a proposal, so opposite sign)
        // effectively an SNR cut on lines
        logpy += zeta*(double)(lines_y->n);
        logpx += zeta*(double)(lines_x->n);


        logH = (logLy - logLx)/heat - logpy + logpx;
        alpha = log(gsl_rng_uniform(r));

        if(logH > alpha)
        {
          if(mc%1000 < 10)
          {
            baseav = 0.0;
            for(i=0; i< ncut; i++)
            {
              baseav += log(Sbase[i]);
            }
            baseav /= (double)(ncut);
          }

          if(typ == 0) ac0++;
          if(typ == 1) ac1++;
          logLx = logLy;
          lines_x->n = lines_y->n;
          for(i=0; i< nspline; i++) sdatax[i] = sdatay[i];
          for(i=0; i< tmax; i++)
          {
            lines_x->larray[i] = lines_y->larray[i];
            lines_x->A[i] = lines_y->A[i];
            lines_x->f[i] = lines_y->f[i];
            lines_x->Q[i] = lines_y->Q[i];
          }
        }
      }
    }

    passes++;

    spass *= 2;

    if(lines_x->n > 20) psmax = 3;

    if(opts.verbose)fprintf(stdout,"pass = %d  n = %d\n", passes, lines_x->n);


  }while((lines_x->n > 5) && (passes < psmax));


  QQ = -1.0;
  x = 0.0;
  // re-map the array to 0..mx ordering
  for(i=0; i< lines_x->n; i++)
  {
    k = lines_x->larray[i];
    lines_y->f[i] = lines_x->f[k];
    lines_y->A[i] = lines_x->A[k];
    lines_y->Q[i] = lines_x->Q[k];
    if(lines_x->Q[k] > QQ)
    {
      QQ = lines_x->Q[k];
      x = lines_x->A[k];
    }
  }

  if(opts.verbose) fprintf(stdout,"Nlines=%i: Highest Q %e, A of this line  %e\n", lines_x->n, QQ, x);

  // return the last value of the chain
  for(i=0; i< lines_x->n; i++)
  {
    lines_x->f[i] = lines_y->f[i];
    lines_x->A[i] = lines_y->A[i];
    lines_x->Q[i] = lines_y->Q[i];
  }

  free(fprop);
  free(Sn);
  free(Sbase);
  destroy_splineParams(spline_y);
  destroy_lorentzianParams(lines_y);
}

void LorentzSplineMCMC(bayesline_opt opts, int steps, dataParams *data, double *freq, double *power, lorentzianParams *lines_x, splineParams *spline, splineParams *spline_x, double *Snf, int focus, double *dan, gsl_rng *r)
{
  fprintf(stdout,"\nLorentzSplineMCMC( focus=%i )\n",focus);
  int skip;
  int nsy, nsx;
  double logLx, logLy, logH;
  int ilowx, ihighx, ilowy, ihighy;
  int i, j, k, ki, ii, jj, mc;
  int check;
  double alpha;
  double lSAmax, lSAmin;
  double lQmin, lQmax;
  double lAmin, lAmax;
  int ac0, ac1, ac2;
  int cc0, cc1, cc2;
  double *Sn, *Sbase, *Sbasex, *Sline, *Snx;
  double *xint;
  double e1, e2, e3, e4;
  double x2, x3, x4;
  double s1, s2, s3, s4;
  int typ;
  double xsm, pmax, fcl, fch, dff;
  double baseav;
  double logpx, logpy, x, y, z, beta;
  double Ac;
  double *fprop;
  double *sdatay;
  double *spointsy;
  int *foc;
  FILE *outfile, *out;

  int ncut = data->ncut;
  int tmax = data->tmax;

  double flow  = data->flow;
  double fhigh = data->fhigh;

  int nspline   = spline->n;
  int *nsplinex = &spline_x->n;

  double *sdatax = spline_x->data;

  double *spoints  = spline->points;
  double *spointsx = spline_x->points;

  //control screen output based on verbose
  if(opts.verbose) skip=100;
  else            skip=1000;


  Snx    = malloc((size_t)(sizeof(double)*(ncut)));
  Sn     = malloc((size_t)(sizeof(double)*(ncut)));
  Sbasex = malloc((size_t)(sizeof(double)*(ncut)));
  Sbase  = malloc((size_t)(sizeof(double)*(ncut)));
  Sline  = malloc((size_t)(sizeof(double)*(ncut)));

  sdatay   = malloc((size_t)(sizeof(double)*(nspline)));
  spointsy = malloc((size_t)(sizeof(double)*(nspline)));

  foc     = malloc((size_t)(sizeof(int)*(tmax)));

  // This keeps track of whos who in the Lorentzian model
  // Necessary complication when using delta likelihood
  // calculations that only update a single line
  lorentzianParams *lines_y = malloc(sizeof(lorentzianParams));
  create_lorentzianParams(lines_y,tmax);

  fprop = malloc((size_t)(sizeof(double)*(ncut)));
  xint  = malloc((size_t)(sizeof(double)*(ncut)));

  // maxima and minima for the noise spectal slopes and amplitudes
  // uniform priors in slope and log amplitude
  lQmin = log(Qmin);
  lQmax = log(Qmax);
  lAmin = log(Amin);
  lAmax = log(Amax);
  lSAmin = log(SAmin);
  lSAmax = log(SAmax);

  nsx = *nsplinex;

  dff = 0.01;  // half-width of frequency focus region (used if focus == 1)

  // this is the fractional error estimate on the noise level
  s1 = 1.0/sqrt((double)(ncut));

  s2 = 0.01;
  s3 = 0.5;
  s4 = 0.5;

  for(i=0; i<lines_x->n; i++)
  {
    lines_x->larray[i] = i;
    lines_y->larray[i] = i;
  }

  baseav = 0.0;

  //Interpolate {spointsx,sdatax} ==> {freq,xint}
  CubicSplineGSL(nsx,spointsx,sdatax,ncut,freq,xint);

  for(i=0; i<ncut; i++)
  {
    Sbase[i] = exp(xint[i]);
    Sbasex[i] = Sbase[i];
    baseav += xint[i];
  }

  baseav /= (double)(ncut);

  full_spectrum_spline(Sline, Sbase, freq, data, lines_x);
  for(i=0; i< ncut; i++) Sn[i] = Sbase[i]+Sline[i];

  for(i=0; i<ncut; i++) Snx[i] = Sn[i];

  if(!opts.zerologL) logLx = loglike(power, Sn, ncut);
  else              logLx = 1.0;

  // set up proposal for frequency jumps
  xsm =0.0;
  pmax = -1.0;
  for(i=0; i< ncut; i++)
	{
    x = power[i]/Sbase[i];
    if(x > pmax)
    {
      pmax = x;
      k = i;
    }
    if(x < 10.0) x = 1.0;
    if(x >= 10.0) x = 100.0;
    fprop[i] = x;
    xsm += x;
  }

  // define the focus region (only used if focus flag = 1)
  fcl = freq[k]-dff;
  fch = freq[k]+dff;
  Ac = power[k];
  if(fcl < freq[0]) fcl = freq[0];
  if(fch > freq[ncut-1]) fch = freq[ncut-1];

  if(opts.verbose)
  {
    out = fopen("fprop.dat","w");
    pmax = -1.0;
    for(i=0; i< ncut; i++)
    {
      fprop[i] /= xsm;
      if(fprop[i] >  pmax) pmax = fprop[i];
      fprintf(out,"%f %e\n", freq[i], fprop[i]);
    }
    fclose(out);


    outfile = fopen("fullspectrum_start.dat","w");
    for(i=0; i< ncut; i++) fprintf(outfile,"%f %e %e %e %e %e\n", freq[i], power[i], Sn[i], power[i]/Sn[i], Sbase[i], Sline[i]);
    fclose(outfile);

    outfile = fopen("Fchain.dat","w");
  }

  ac0 = 0;
  ac1 = 0;
  ac2 = 0;
  cc0 = 1;
  cc1 = 1;
  cc2 = 1;


  for(mc=0; mc < steps; mc++)
  {

    //copy over current state
    lines_y->n = lines_x->n;
    nsy = nsx;
    for(i=0; i< tmax; i++)
    {
      lines_y->larray[i] = lines_x->larray[i];
      lines_y->Q[i] = lines_x->Q[i];
      lines_y->f[i] = lines_x->f[i];
      lines_y->A[i] = lines_x->A[i];
    }
    for(i=0; i<nspline; i++)
    {
      sdatay[i] = sdatax[i];
      spointsy[i] = spointsx[i];
    }

    beta = gsl_rng_uniform(r);

    if(beta > 0.9)  // update the smooth part of the spectrum
    {

      alpha = gsl_rng_uniform(r);

      logpx = logpy = 0.0;

      if(alpha > 0.8)  // try a transdimensional move
      {
        alpha = gsl_rng_uniform(r);
        if(alpha > 0.5)// || nsx<3)  // try and add a new term
        {
          nsy = nsx+1;
          typ = 5;
        }
        else // try and remove term
        {
          nsy = nsx-1;
          typ = 6;
        }

        if(nsy > 0 && nsy < nspline)
        {

          if(nsy < nsx)
          {

            ki=1+(int)(gsl_rng_uniform(r)*(double)(nsx-2)); // pick a term to try and kill - cant be first or last term
            k = 0;
            for(j=0; j<nsx; j++)
            {
              if(j != ki)
              {
                sdatay[k] = sdatax[j];
                spointsy[k] = spointsx[j];
                k++;
              }
            }

          }

          if(nsy > nsx)
          {

            // have to randomly pick a new point that isn't already in use
            do
            {
              ki=1+(int)(gsl_rng_uniform(r)*(double)(nspline-2));  // pick a point to add
              ii = 0;
              for(j=0; j<nsx; j++)
              {
                if(fabs(spointsx[j]-spoints[ki]) < 1.0e-3) ii = 1;  // this means the point is already in use
              }
            } while (ii == 1);
            ii = 0;
            for(j=0; j<nsx; j++)
            {

              if(spointsx[j] < spoints[ki])
              {
                sdatay[j] = sdatax[j];
                spointsy[j] = spointsx[j];
              }

              if((spointsx[j] > spoints[ki]) && ii == 0)  // found where to slot the new term in
              {
                sdatay[j] = lSAmin +(lSAmax - lSAmin)*gsl_rng_uniform(r);
                spointsy[j] = spoints[ki];
                ii = 1;
              }

              if((spointsx[j] > spoints[ki]) && ii == 1)
              {
                sdatay[j+1] = sdatax[j];
                spointsy[j+1] = spointsx[j];
              }

            }
          }
        }
      }
      else  // regular MCMC update
      {
        typ = 4;

        nsy = nsx;

        //pick a term to update
        ii=(int)(gsl_rng_uniform(r)*(double)(nsx));

        // use a variety of jump sizes by using a sum of gaussians of different width
        e1 = 0.0005;
        alpha = gsl_rng_uniform(r);
        if(alpha > 0.8)
        {
          e1 = 0.002;
        }
        else if(alpha > 0.6)
        {
          e1 = 0.005;
        }
        else if(alpha > 0.4)
        {
          e1 = 0.05;
        }

        // propose new value for the selected term
        sdatay[ii] = sdatax[ii]+gsl_ran_gaussian(r, e1);

      }

      check = 0;

      if(nsy < 1 || nsy > nspline-1) check = 1;

      for(i=0; i<nsy; i++)
      {
        if(sdatay[i] > lSAmax) check = 1;
        if(sdatay[i] < lSAmin) check = 1;
      }
    }
    else    // update the line model
    {

      alpha = gsl_rng_uniform(r);

      if(alpha > 0.8)  // try a transdimensional move
      {

        alpha = gsl_rng_uniform(r);
        if(alpha < 0.5)  // try and add a new term
        {
          lines_y->n = lines_x->n+1;
          typ = 2;
        }
        else // try and remove term
        {
          lines_y->n = lines_x->n-1;
          typ = 3;
        }

        check = 0;
        if(lines_y->n < 0 || lines_y->n > tmax) check = 1;


        if(check == 0)
        {
          if(lines_y->n < lines_x->n)
          {
            i=(int)(gsl_rng_uniform(r)*(double)(lines_x->n)); // pick a term to try and kill
            k = 0;
            for(j=0; j< lines_x->n; j++)
            {
              if(j != i)
              {
                lines_y->larray[k] = lines_x->larray[j];
                jj = lines_x->larray[j];
                lines_y->A[jj] = lines_x->A[jj];
                lines_y->Q[jj] = lines_x->Q[jj];
                lines_y->f[jj] = lines_x->f[jj];
                k++;
              }
              if(j == i) ki = lines_x->larray[j];  // take note of who's demise is being proposed
            }

            logpx = lprop(lines_x->f[ki], fprop, data);
            logpy = -log((double)(ncut));    // corresponds to uniform density - just as likely to kill from anywhere

            y = log(lines_x->A[ki])-baseav;
            if(y < 0.0)
            {
              z = 0.0;
            }
            else
            {
              z = 2.0*kappa*gsl_ran_gaussian_pdf(y, lAwidth);
            }
            logpx += log((1.0-kappa)/(lAmax-lAmin) + z);
            logpy += -log((lAmax-lAmin));
          }

          if(lines_y->n > lines_x->n)
          {
            // copy over current state
            for(j=0; j< lines_x->n; j++)
            {
              k = lines_x->larray[j];
              lines_y->A[k] = lines_x->A[k];
              lines_y->f[k] = lines_x->f[k];
              lines_y->Q[k] = lines_x->Q[k];
              lines_y->larray[j] = lines_x->larray[j];
            }

            // find a label that isn't in use

            i = -1;
            do
            {
              i++;
              k = 0;
              for(j=0; j< lines_x->n; j++)
              {
                if(i==lines_x->larray[j]) k = 1;
              }
            } while(k == 1);

            lines_y->larray[lines_x->n] = i;
            ii=i;

            // draw new term
            lines_y->A[i] = exp(lAmin+(lAmax-lAmin)*gsl_rng_uniform(r));
            lines_y->Q[i] = exp(lQmin+(lQmax-lQmin)*gsl_rng_uniform(r));
            lines_y->f[i] = sample(fprop, pmax, data, r);
            logpy = lprop(lines_y->f[i], fprop, data);
            logpx = -log((double)(ncut));    // corresponds to uniform density - just as likely to kill from anywhere

            alpha = gsl_rng_uniform(r);
            if(alpha < kappa)
            {
              y = fabs(gsl_ran_gaussian(r, lAwidth));
              lines_y->A[i] = exp(baseav+y);
            }
            else
            {
              lines_y->A[i] = exp(lAmin+(lAmax-lAmin)*gsl_rng_uniform(r));
              y = log(lines_y->A[i])-baseav;
            }
            if(y < 0.0)
            {
              z = 0.0;
            }
            else
            {
              z = 2.0*kappa*gsl_ran_gaussian_pdf(y, lAwidth);
            }
            logpy += log((1.0-kappa)/(lAmax-lAmin) + z);
            logpx += -log((lAmax-lAmin));


            if(focus==1) // using focused region (not Markovian - not paying penalty for proposing in such a small region)
            {
              lines_y->f[i] = fcl+(fch-fcl)*gsl_rng_uniform(r);
              lines_y->Q[i] = Qmax/10.0;
              lines_y->A[i] = Ac;

              logpy = 0.0;
              logpx = 0.0;
            }
            if(lines_y->A[i] > Amax) check = 1;
            if(lines_y->A[i] < Amin) check = 1;
          }
        }
      }
      else  // regular MCMC update
      {
        lines_y->n = lines_x->n;

        if(lines_y->n > 0)
        {

          //pick a term to update
          jj=(int)(gsl_rng_uniform(r)*(double)(lines_x->n));
          //find label of who is geting updated
          ii = lines_x->larray[jj];

          // copy over current state
          for(i=0; i< lines_x->n; i++)
          {
            k = lines_x->larray[i];
            lines_y->f[k] = lines_x->f[k];
            lines_y->A[k] = lines_x->A[k];
            lines_y->Q[k] = lines_x->Q[k];
            lines_y->larray[i] = lines_x->larray[i];
          }

          if(focus == 1)
          {
            // find if any lines are in the focus region
            j = 0;
            for(i=0; i< lines_x->n; i++)
            {
              k = lines_x->larray[i];
              if(lines_x->f[k] > fcl && lines_x->f[k] < fch)
              {
                foc[j] = k;
                j++;
              }
            }

            x = 0.0;
            if(j > 0)  // some lines are currently in the focus region
            {
              x = 0.8;
              jj=(int)(gsl_rng_uniform(r)*(double)(j));
              //find label of who is getting updated in the focus region
              ii = foc[jj];
            }

          }
          else
          {
            x = 0.8;
          }

          alpha = gsl_rng_uniform(r);
          if(alpha > x)
          {
            // here we try and move an exisiting line to a totally new location
            if(focus != 1)
            {
              lines_y->f[ii] = sample(fprop, pmax, data, r);
              logpy = lprop(lines_y->f[ii], fprop, data);
              logpx = lprop(lines_x->f[ii], fprop, data);

              lines_y->Q[ii] = exp(lQmin+(lQmax-lQmin)*gsl_rng_uniform(r));
              alpha = gsl_rng_uniform(r);
              if(alpha < kappa)
              {
                y = fabs(gsl_ran_gaussian(r, lAwidth));
                lines_y->A[ii] = exp(baseav+y);
              }
              else
              {
                lines_y->A[ii] = exp(lAmin+(lAmax-lAmin)*gsl_rng_uniform(r));
                y = log(lines_y->A[ii])-baseav;
              }
              if(y < 0.0)
              {
                z = 0.0;
              }
              else
              {
                z = 2.0*kappa*gsl_ran_gaussian_pdf(y, lAwidth);
              }
              logpy += log((1.0-kappa)/(lAmax-lAmin) + z);
              logpx += -log((lAmax-lAmin));


            }
            else  // using focused region (not Markovian - not paying penalty for proposing in such a small region)
            {
              lines_y->f[ii] = fcl+(fch-fcl)*gsl_rng_uniform(r);
              lines_y->Q[ii] = Qmax/10.0;
              lines_y->A[ii] = Ac;

              logpy = 0.0;
              logpx = 0.0;
            }
            typ = 0;
          }
          else
          {
            typ = 1;
            alpha = gsl_rng_uniform(r);

            if     (alpha > 0.9) beta = 1.0e+1;
            else if(alpha > 0.6) beta = 1.0e+0;
            else if(alpha > 0.3) beta = 1.0e-1;
            else                 beta = 1.0e-2;

            e2 = beta*s2;
            e3 = beta*s3;
            e4 = beta*s4;
            x2 = gsl_ran_gaussian(r, e2);
            x4 = gsl_ran_gaussian(r, e4);
            x3 = gsl_ran_gaussian(r, e3);

            lines_y->A[ii] = lines_x->A[ii]*exp(x3);
            lines_y->Q[ii] = lines_x->Q[ii]*exp(x4);
            lines_y->f[ii] = lines_x->f[ii]+x2;

            logpx = logpy = 0.0;
          }

          check =0;

          if(lines_y->A[ii] > Amax)  check = 1;
          if(lines_y->A[ii] < Amin)  check = 1;
          if(lines_y->f[ii] < flow)  check = 1;
          if(lines_y->f[ii] > fhigh) check = 1;
          if(lines_y->Q[ii] < Qmin)  check = 1;
          if(lines_y->Q[ii] > Qmax)  check = 1;

        }
      }
    }

    //If line parameters satisfy priors, continue with MCMC update
    if(check == 0)
    {

      if(typ == 0) cc0++;
      if(typ == 1) cc1++;
      if(typ == 4) cc2++;

      if(!opts.zerologL)
      {
        if(typ > 3)  // need to do a full update (slower)
        {
          if(nsy>2)
          {

            //Interpolate {spointsy,sdatay} ==> {freq,xint}
            CubicSplineGSL(nsy,spointsy,sdatay,ncut,freq,xint);

            for(i=0; i<ncut; i++)  Sbase[i] = exp(xint[i]);

            full_spectrum_spline(Sline, Sbase, freq, data, lines_y);
            for(i=0; i< ncut; i++) Sn[i] = Sbase[i]+Sline[i];
            logLy = loglike(power, Sn, ncut);
          }
          else logLy = -1e60;
        }

        if(typ == 1 || typ == 0)  // fixed dimension MCMC of line ii
        {
          full_spectrum_single(Sn, Snx, Sbasex, freq, data, lines_x, lines_y, ii, &ilowx, &ihighx, &ilowy, &ihighy);
          logLy = logLx + loglike_single(power, Sn, Snx, ilowx, ihighx, ilowy, ihighy);
        }

        if(typ == 2)  // add new line with label ii
        {
          full_spectrum_add_or_subtract(Sn, Snx, Sbasex, freq, data, lines_y, ii, &ilowy, &ihighy,  1);
          logLy = logLx + loglike_pm(power, Sn, Snx, ilowy, ihighy);
        }

        if(typ == 3)  // remove line with label ki
        {
          full_spectrum_add_or_subtract(Sn, Snx, Sbasex, freq, data, lines_x, ki, &ilowx, &ihighx, -1);
          logLy = logLx + loglike_pm(power, Sn, Snx, ilowx, ihighx);
        }

      }
      else
      {
        logLy = 1.0;
      }

      // prior on line number e(-zeta * n).  (this is a prior, not a proposal, so opposite sign)
      // effectively an SNR cut on lines
      logpy += zeta*(double)(lines_y->n);
      logpx += zeta*(double)(lines_x->n);

      logH  = (logLy - logLx) - logpy + logpx;
      alpha = log(gsl_rng_uniform(r));

      if(logH > alpha)
      {
        if(typ == 0) ac0++;
        if(typ == 1) ac1++;
        if(typ == 4) ac2++;
        logLx = logLy;
        lines_x->n = lines_y->n;
        nsx = nsy;
        for(i=0; i< ncut; i++)
        {
          Snx[i] = Sn[i];
          if(typ > 3) Sbasex[i] = Sbase[i];
        }
        for(i=0; i< tmax; i++)
        {
          lines_x->larray[i] = lines_y->larray[i];
          lines_x->A[i] = lines_y->A[i];
          lines_x->f[i] = lines_y->f[i];
          lines_x->Q[i] = lines_y->Q[i];
        }
        for(i=0; i<nspline; i++)
        {
          sdatax[i] = sdatay[i];
          spointsx[i] = spointsy[i];
        }
      }

    }//end prior check

    //Every 1000 steps update focus region
    if(mc%1000 == 0)
    {
      pmax = -1.0;
      for(i=0; i< ncut; i++)
      {
        x = power[i]/Snx[i];
        if(x > pmax)
        {
          pmax = x;
          k = i;
        }
      }

      // define the focus region (only used if focus flag = 1)
      fcl = freq[k]-dff;
      fch = freq[k]+dff;
      Ac = power[k];
      if(fcl < freq[0]) fcl = freq[0];
      if(fch > freq[ncut-1]) fch = freq[ncut-1];

      if(focus==1)fprintf(stdout,"Focusing on [%f %f] Hz...\n", fcl, fch);


      // set up proposal for frequency jumps
      xsm =0.0;
      baseav = 0.0;
      for(i=0; i< ncut; i++)
      {
        x = power[i]/Snx[i];
        if(x < 16.0) x = 1.0;
        if(x >= 16.0) x = 100.0;
        fprop[i] = x;
        xsm += x;
        baseav += log(Sbasex[i]);
      }
      baseav /= (double)(ncut);


      pmax = -1.0;
      for(i=0; i< ncut; i++)
      {
        fprop[i] /= xsm;
        if(fprop[i] > pmax) pmax = fprop[i];
      }

    }//end focus update


    //Output chain stats to screen (and file if verbose)
    if(mc%skip == 0)
    {

      if(mc==0)fprintf(stdout,"#cycle logLx logLy logH Nlines Nspline proptype check line_b/d_acc line_acc spline_acc\n");
      fprintf(stdout,"%d %.9e %.9e %.9e %d %d %d %d %f %f %f\n", mc, logLx, logLy, logH, lines_x->n, nsx, typ, check,
              (double)(ac0)/(double)(cc0), (double)(ac1)/(double)(cc1), (double)(ac2)/(double)(cc2));

      if(opts.verbose)
      {
        fprintf(outfile,"%d %.9e %d %d ", mc/10, logLx, lines_x->n, nsx);
        for(i=0; i< tmax/2; i++) fprintf(outfile,"%f %e %e ", lines_x->f[i], lines_x->A[i], lines_x->Q[i]);
        fprintf(outfile,"\n");
      }
    }

  }//End MCMC loop

  //Interpolate {spointsx,sdatax} ==> {freq,xint}
  CubicSplineGSL(nsx,spointsx,sdatax,ncut,freq,xint);
  for(i=0; i< ncut; i++) Sbase[i] = exp(xint[i]);

  full_spectrum_spline(Sline, Sbase, freq, data, lines_x);
  for(i=0; i< ncut; i++) Sn[i] = Sbase[i]+Sline[i];

  if(opts.verbose)
  {
    fclose(outfile);

    outfile = fopen("fullspectrum_end.dat","w");
    for(i=0; i< ncut; i++)
    {
      fprintf(outfile,"%f %e %e %e %e %e %e\n", freq[i], power[i], Sn[i], power[i]/Sn[i], Sbase[i], Sline[i], Snx[i]);
    }
    fclose(outfile);
  }

  // return updated spectral model
  for(i=0; i< ncut; i++) Snf[i] = Snx[i];

  // re-map the array to 0..mx ordering
  for(i=0; i< lines_x->n; i++)
  {
    k = lines_x->larray[i];
    lines_y->f[i] = lines_x->f[k];
    lines_y->A[i] = lines_x->A[k];
    lines_y->Q[i] = lines_x->Q[k];
  }

  // return the last value of the chain
  for(i=0; i< lines_x->n; i++)
  {
    lines_x->f[i] = lines_y->f[i];
    lines_x->A[i] = lines_y->A[i];
    lines_x->Q[i] = lines_y->Q[i];
  }

  // check for outliers
  pmax = -1.0;
  for(i=0; i< ncut; i++)
  {
    x = power[i]/Snx[i];
    if(x > pmax) pmax = x;
  }

  *dan = pmax;
  *nsplinex = nsx;

  free(foc);
  free(xint);
  free(fprop);
  free(Snx);
  free(Sn);
  free(Sbase);
  free(Sbasex);
  free(Sline);
  free(sdatay);
  free(spointsy);

  destroy_lorentzianParams(lines_y);

  fprintf(stdout, "  PSD model used %i points, %i lines, dan factor = %g\n",nsx,lines_x->n,pmax);
}

void spectrum_spline(double *Sn, double *Sbase, double *sfreq, dataParams *data, lorentzianParams *lines, splineParams *spline)
{
  int i, j, k, n;
  int istart, istop;

  int nspline     = spline->n;
  double *spoints = spline->points;
  double *sdata   = spline->data;

  double *x,*Stemp;

  n = data->ncut;

  x     = malloc((size_t)(sizeof(double)*(n)));
  Stemp = malloc((size_t)(sizeof(double)*(n)));

  //Interpolate {spoints,sdata} ==> {sfreq,x}
  CubicSplineGSL(nspline,spoints,sdata,n,sfreq,x);

  for(i=0; i< n; i++)
	{
    Sbase[i] = exp(x[i]);          // spline base model
    Sn[i] = 0.0;
  }

  for(k=0; k<lines->n; k++)
  {
    j = lines->larray[k];
    for(i=0; i<n; i++) Stemp[i]=Sn[i];
    full_spectrum_add_or_subtract(Sn, Stemp, Sbase, sfreq, data, lines, j, &istart, &istop, 1);
  }

  for(i=0; i< n; i++)
  {
    Sn[i] += Sbase[i];
  }

  free(x);
  free(Stemp);
}

double loglike_pm(double *respow, double *Sn, double *Snx, int ilow, int ihigh)
{
  double lgl, x;
  int i;

  // leavimng out the log(2Pi) terms since they cancel in Hastings ratio
  lgl = 0.0;
  for(i=ilow; i< ihigh; i++)
  {
    x = respow[i]/Sn[i]-respow[i]/Snx[i];
    lgl -= (x+log(Sn[i]/Snx[i]));
  }

  return(lgl);
}



double loglike_single(double *respow, double *Sn, double *Snx, int ilowx, int ihighx, int ilowy, int ihighy)
{
  double lgl, x;
  int i;
  int ilow, ihigh;
  int imid1, imid2;

  ilow = ilowx;
  if(ilowy < ilow) ilow = ilowy;

  if(ilow == ilowx)
  {
    if(ihighx <= ilowy)  // separate regions
    {
      imid1 = ihighx;
      imid2 = ilowy;
      ihigh = ihighy;
    }

    if(ihighx > ilowy) // overlapping regions
    {
      if(ihighx < ihighy)
      {
        imid1 = ihighx;
        imid2 = ihighx;
        ihigh = ihighy;
      }
      else
      {
        imid1 = ilowy;
        imid2 = ilowy;
        ihigh = ihighx;
      }
    }
  }

  if(ilow == ilowy)
  {
    if(ihighy <= ilowx)  // separate regions
    {
      imid1 = ihighy;
      imid2 = ilowx;
      ihigh = ihighx;
    }

    if(ihighy > ilowx) // overlapping regions
    {
      if(ihighy < ihighx)
      {
        imid1 = ihighy;
        imid2 = ihighy;
        ihigh = ihighx;
      }
      else
      {
        imid1 = ilowx;
        imid2 = ilowx;
        ihigh = ihighy;
      }
    }
  }

  // leavimng out the log(2Pi) terms since they cancel in Hastings ratio
  lgl = 0.0;
  for(i=ilow; i< imid1; i++)
  {
    x = respow[i]/Sn[i]-respow[i]/Snx[i];
    lgl -= (x+log(Sn[i]/Snx[i]));
  }
  for(i=imid2; i< ihigh; i++)
  {
    x = respow[i]/Sn[i]-respow[i]/Snx[i];
    lgl -= (x+log(Sn[i]/Snx[i]));
  }

  return(lgl);
}

void full_spectrum_add_or_subtract(double *Snew, double *Sold, double *Sbase, double *sfreq, dataParams *data, lorentzianParams *lines, int ii, int *ilow, int *ihigh, int flag)
{
  int i;
  double dS,f2,f4;
  double deltf;
  double fsq, x, z, deltafmax, spread;
  int istart, istop, imid, idelt;

  double A = lines->A[ii];
  double Q = lines->Q[ii];
  double f = lines->f[ii];

  int    ncut = data->ncut;
  double Tobs = data->Tobs;
  double flow = data->flow;

  // copy over current model
  for(i=0; i< ncut; i++) Snew[i] = Sold[i];

  // here we figure out how many frequency bins are needed for the line profile
  imid = (int)((f-flow)*Tobs);
  spread = (1.0e-2*Q);

  if(spread < 50.0) spread = 50.0;  // maximum half-width is f_resonance/50
  deltafmax = f/spread;
  deltf = 4.0*deltafmax;
  idelt = (int)(deltf*Tobs)+1;
  if(A < 10.0*Sbase[imid]) idelt = (int)(20.0*f*Tobs/Q)+1;

  istart = imid-idelt;
  istop = imid+idelt;
  if(istart < 0) istart = 0;
  if(istop > ncut) istop = ncut;

  *ilow  = istart;
  *ihigh = istop;


  // add or remove the old line
  f2=f*f;
  f4=f2*f2;
  for(i=istart; i<istop; i++)
  {
    fsq = sfreq[i]*sfreq[i];
    x = fabs(f-sfreq[i]);
    z = 1.0;
    if(x > deltafmax) z = exp(-(x-deltafmax)/deltafmax);
    dS = z*A*f4/(f2*fsq+Q*Q*(fsq-f2)*(fsq-f2));
    switch(flag)
    {
      case 1: //add new line
        Snew[i] += dS;
        break;
      case -1: //remove line
        Snew[i] -= dS;
        break;
      default:
        break;
    }
  }
}

void full_spectrum_single(double *Sn, double *Snx, double *Sbasex, double *sfreq, dataParams *data, lorentzianParams *line_x, lorentzianParams *line_y, int ii,
                          int *ilowx, int *ihighx, int *ilowy, int *ihighy)
{

  double *Stemp = malloc((size_t)(sizeof(double)*(data->ncut)));

  full_spectrum_add_or_subtract(Stemp, Snx, Sbasex, sfreq, data, line_x, ii, ilowx, ihighx,-1);
  full_spectrum_add_or_subtract(Sn, Stemp,  Sbasex, sfreq, data, line_y, ii, ilowy, ihighy, 1);

  free(Stemp);
}



void full_spectrum_spline(double *Sline, double *Sbase, double *sfreq, dataParams *data, lorentzianParams *lines)
{
  int i, j, k;
  int istart, istop;

  double *Stemp = malloc((size_t)(sizeof(double)*(data->ncut)));

  for(i=0; i<data->ncut; i++) Sline[i] = 0.0;
  for(k=0; k<lines->n; k++)
  {
    j = lines->larray[k];
    for(i=0; i<data->ncut; i++) Stemp[i]=Sline[i];
    full_spectrum_add_or_subtract(Sline, Stemp, Sbase, sfreq, data, lines, j, &istart, &istop, 1);
  }

  free(Stemp);
}

void SpecFitSpline(bayesline_opt opts, int steps, double *freq, double *power, splineParams *spline, double *Snf, int ncut, gsl_rng *r)
{
  printf("\nSpecFitSpline()\n");
  int i, j, k, ii, ki, mc;
  int nsy, nsx;
  int check;

  int nspline     = spline->n;
  double *sdata   = spline->data;
  double *spoints = spline->points;

  double *sdatax, *sdatay;
  double *spointsx, *spointsy;
  double *Snfx;

  double alpha;
  double logLx, logLy, logH;
  double lSAmax, lSAmin;
  double e1;

  nsx = nspline;

  // maxima and minima for the noise spectal amplitudes
  lSAmin = log(SAmin);
  lSAmax = log(SAmax);

  sdatax   = malloc((size_t)(sizeof(double)*(nspline)));
  sdatay   = malloc((size_t)(sizeof(double)*(nspline)));
  spointsx = malloc((size_t)(sizeof(double)*(nspline)));
  spointsy = malloc((size_t)(sizeof(double)*(nspline)));
  Snfx     = malloc((size_t)(sizeof(double)*(ncut+1)));

  for(i=0; i<nspline; i++)
  {
    sdatax[i]   = sdata[i];
    spointsx[i] = spoints[i];
  }


  // check that line amplitudes are within prior ranges
  for(i=0; i<nsx; i++)
  {
    if(sdatax[i] > lSAmax || sdatax[i] < lSAmin)
    {
      printf("%i: %lg < %lg < %lg ?\n",i, exp(lSAmin), exp(sdatax[i]), exp(lSAmax));
      KILL("spectrum priors not wide enough\n");
    }
  }

  //Interpolate {spointsx,sdatax} ==> {freq,Snf}
  CubicSplineGSL(nspline, spointsx, sdatax, ncut, freq, Snf);

  if(!opts.zerologL) logLx = loglike_fit_spline(power, Snf, ncut);
  else               logLx = 1.0;

  for(mc=0; mc < steps; mc++)
	{
    // copy over the current state
    nsy = nsx;
    for(i=0; i<nsx; i++)
    {
      sdatay[i]   = sdatax[i];
      spointsy[i] = spointsx[i];
    }

    alpha = gsl_rng_uniform(r);

    if(alpha > 0.8)  // try a transdimensional move
    {
      alpha = gsl_rng_uniform(r);

      if(alpha > 0.5) nsy = nsx+1; // try and add a new term
      else            nsy = nsx-1; // try and remove term

      if(nsy > 2 && nsy <= nspline)
      {

        if(nsy < nsx)
        {

          ki=1+(int)(gsl_rng_uniform(r)*(double)(nsx-2)); // pick a term to try and kill - cant be first or last term
          k = 0;
          for(j=0; j<nsx; j++)
          {
            if(j != ki)
            {
              sdatay[k]   = sdatax[j];
              spointsy[k] = spointsx[j];
              k++;
            }
          }

        }

        if(nsy > nsx)
        {
          // have to randomly pick a new point that isn't already in use
          do
          {
            ki=1+(int)(gsl_rng_uniform(r)*(double)(nspline-2));  // pick a point to add
            ii = 0;
            for(j=0; j<nsx; j++)
            {
              if(fabs(spointsx[j]-spoints[ki]) < 1.0e-3) ii = 1;  // this means the point is already in use
            }
          } while (ii == 1);

          ii = 0;
          for(j=0; j<nsx; j++)
          {
            if(spointsx[j] < spoints[ki])
            {
              sdatay[j]   = sdatax[j];
              spointsy[j] = spointsx[j];
            }

            if((spointsx[j] > spoints[ki]) && ii == 0)  // found where to slot the new term in
            {
              sdatay[j]   = lSAmin +(lSAmax - lSAmin)*gsl_rng_uniform(r);
              spointsy[j] = spoints[ki];
              ii = 1;
            }

            if((spointsx[j] > spoints[ki]) && ii == 1)
            {
              sdatay[j+1]   = sdatax[j];
              spointsy[j+1] = spointsx[j];
            }

          }
        }
      }
    }//end transdimensional move

    else  // regular MCMC update
    {
      nsy = nsx;

      //pick a term to update
      ii=(int)(gsl_rng_uniform(r)*(double)(nsx));

      // use a variety of jump sizes by using a sum of gaussians of different width
      e1 = 0.0005;
      alpha = gsl_rng_uniform(r);

      if(alpha > 0.8)      e1 = 0.002;
      else if(alpha > 0.7) e1 = 0.005;
      else if(alpha > 0.6) e1 = 0.05;

      // propose new value for the selected term
      sdatay[ii] = sdatax[ii]+gsl_ran_gaussian(r, e1);
    }


    //check that priors on dimension and line amplitude are satisfied
    check = 0;
    if(nsy < 3 || nsy > nspline-1) check = 1;
    for(i=0; i<nsy; i++) if(sdatay[i] > lSAmax || sdatay[i] < lSAmin) check = 1;

    //skip interpolation & Hastings ratio if parameters are out of bounds
    if(check == 0)
    {
      //interpolate {spointsy,sdatay} ==> {freq,Snf}
      CubicSplineGSL(nsy, spointsy, sdatay, ncut, freq, Snf);

      alpha = log(gsl_rng_uniform(r));

      if(!opts.zerologL) logLy = loglike_fit_spline(power, Snf, ncut);
      else               logLy = 1.0;

      logH  = logLy - logLx;

      if(logH > alpha)
      {
        logLx = logLy;
        nsx = nsy;
        for(i=0; i< ncut; i++) Snfx[i] = Snf[i];
        for(i=0; i<nsx; i++)
        {
          sdatax[i]   = sdatay[i];
          spointsx[i] = spointsy[i];
        }
      }

    }//end prior check

  }

  // return the most recent accepted estimate for the spectrum
  for(i=0; i< ncut; i++) Snf[i] = Snfx[i];

  fprintf(stdout,"  Spline used %d of %d points\n", nsx, nspline);

  // pass back fit

  //interpolate {spointsx,sdatax} ==> {spoints,sdata}
  CubicSplineGSL(nsx, spointsx, sdatax, nspline, spoints, sdata);

  free(Snfx);
  free(sdatax);
  free(sdatay);
  free(spointsx);
  free(spointsy);

}

void CubicSplineGSL(int N, double *x, double *y, int Nint, double *xint, double *yint)
{
  int n;
  double tol=1.0e-6;


  /* set up GSL cubic spline */
  gsl_spline       *cspline = gsl_spline_alloc(gsl_interp_cspline, N);
  gsl_interp_accel *acc    = gsl_interp_accel_alloc();

  /* get derivatives */
  gsl_spline_init(cspline,x,y,N);

  /* interpolate */
  for(n=0; n<Nint; n++)
  {
    /*
     GSL cubic spline throws errors if
     interpolated points are at end of
     spline control points
     */
    if     (fabs(xint[n]-x[0])<tol)
      yint[n] = y[0];

    else if(fabs(xint[n]-x[N-1])<tol)
      yint[n] = y[N-1];

    else
      yint[n]=gsl_spline_eval(cspline,xint[n],acc);
  }

  gsl_spline_free (cspline);
  gsl_interp_accel_free (acc);

}

void create_dataParams(dataParams *data, double *f, int n)
{

  // length of segment in seconds, this should be read in from the frame file
  data->Tobs = rint(1.0/(f[1]-f[0]));

  // frequency resolution
  data->df = 1.0/data->Tobs;

  // sample cadence in Hz, this should be read in from the frame file
  data->cadence = pow(2.0,rint(log((double)(n))/log(2.0))+1.0)/data->Tobs;

  // Nyquist frequency
  data->fny = 2.0/data->cadence;

  // size of segments in Hz
  // If frequency snippets are too large need longer initial runs to get convergence
  data->fstep = 9.0;

  // This sets the maximum number of Lorentzian lines per segment.
  // For segements shorter than 16 Hz this always seems to be enough
  data->tmax = 40;

  // approximate number of segments
  data->sgmts = (int)((f[n-1]-f[0])/data->fstep)+2;

  // Maximum number of lines for full data set
  data->tfull = data->tmax*data->sgmts+1;

  //minimum frequency
  data->fmin = f[0];

  //maximum frequency
  data->fmax = f[n-1];

  //minimum Fourier bin
  data->nmin = (int)(f[0]*data->Tobs);

  // the stencil separation in Hz for the spline model. Going below 2 Hz is dangerous - will fight with line model
  data->fgrid = 4.0;
}

void create_lorentzianParams(lorentzianParams *lines, int size)
{
  lines->n    = 0;
  lines->size = size;

  lines->larray = malloc((size_t)(sizeof(int)*size));

  lines->f = malloc((size_t)(sizeof(double)*size));
  lines->Q = malloc((size_t)(sizeof(double)*size));
  lines->A = malloc((size_t)(sizeof(double)*size));
}

void destroy_lorentzianParams(lorentzianParams *lines)
{
  free(lines->larray);
  free(lines->f);
  free(lines->Q);
  free(lines->A);
  free(lines);
}

void create_splineParams(splineParams *spline, int size)
{
  spline->n = size;

  spline->data   = malloc((size_t)(sizeof(double)*size));
  spline->points = malloc((size_t)(sizeof(double)*size));
}

void destroy_splineParams(splineParams *spline)
{
  free(spline->data);
  free(spline->points);
  free(spline);
}

void create_psdParams(psdParams *psd, int size)
{
  psd->n = size;

  psd->f     = malloc( (size_t)(sizeof(double)*size) );
  psd->pow   = malloc( (size_t)(sizeof(double)*size) );

  psd->Sn    = malloc( (size_t)(sizeof(double)*size) );
  psd->Sbase = malloc( (size_t)(sizeof(double)*size) );
  psd->Sline = malloc( (size_t)(sizeof(double)*size) );
}

void destroy_psdParams(psdParams *psd)
{
  free(psd->f);
  free(psd->pow);
  free(psd->Sn);
  free(psd->Sbase);
  free(psd->Sline);
  free(psd);
}

void gnuplot_full_spectrum(int n, double *f, double *p, double *Sn, FILE *pipe)
{
  int i;
  FILE *outfile = fopen("fullspectrum_final.dat","w");
  for(i=0; i<n; i++) fprintf(outfile,"%f %e %e %e\n", f[i], p[i], Sn[i], p[i]/Sn[i]);
  fclose(outfile);

  fprintf(pipe, "plot 'fullspectrum_final.dat' using 1:4 title 'whitened' with lines\n");
  fflush(pipe);
}


int opt_parse( bayesline_opt *opts, int argc, char **argv )
{
  int no_arg  = 0;   // no argument
  int req_arg = 1;  // required argument

  static bayesline_opt cmd_param;

  struct option cmd_opt[] =
  {
    { "gnuplot",  no_arg,  0, 'g' },
    { "help",     no_arg,  0, 'h' },
    { "infile",  req_arg,  0, 'i' },
    { "zeroLogL", no_arg,  0, 'l' },
    { "outfile", req_arg,  0, 'o' },
    { "verbose",  no_arg,  0, 'v' },
    { 0, 0, 0, 0 }
  };

  char args[] = "ghi:lo:v";
  char *program = argv[0];

  if(argc==1)
  {
    bayesline_usage(program);
    exit(0);
  }

  //set defaults
  cmd_param.gnuplot  = 0;
  cmd_param.ifile[0] = 0;
  cmd_param.zerologL = 0;
  cmd_param.ofile[0] = 0;
  cmd_param.verbose  = 0;

  while(1)
  {
    int opt_indx = 0;
    int c;

    c = getopt_long( argc, argv, args, cmd_opt, &opt_indx );
    if ( c == -1 ) // end of options
      break;

    switch ( c )
    {
      case 0: // if option set a flag, nothing else to do
        if ( cmd_opt[opt_indx].flag ) break;
        else
        {
          fprintf(stderr,"error parsing option %s with argument %s\n", cmd_opt[opt_indx].name, optarg );
          exit(1);
        }
      case 'g': // gnuplot output
        cmd_param.gnuplot = 1;
        cmd_param.verbose = 1;
        break;
      case 'h': // help
        bayesline_usage(program);
        exit(0);
      case 'i': // in-file
        sprintf(cmd_param.ifile,"%s",optarg);
        break;
      case 'l': // likelihood=constant test
        cmd_param.zerologL = 1;
      case 'o': // out-file
        sprintf(cmd_param.ofile,"%s",optarg);
        break;
      case 'v': // verbose output
        cmd_param.verbose = 1;
        break;
      case '?':
        fprintf(stderr,"unknown error while parsing options\n" );
        exit(1);
      default:
        fprintf(stderr,"unknown error while parsing options\n" );
        exit(1);
    }
  }
  
  if ( optind < argc )
  {
    fprintf( stderr, "extraneous command line arguments:\n" );
    while ( optind < argc ) fprintf( stderr, "%s\n", argv[optind++] );
    exit(1);
  }
  
  //check that all necessary command line arguments were supplied
  if(!cmd_param.ifile[0])
  {
    fprintf(stderr,"ERROR: required option '--infile=filename', must specify input data file\n");
    bayesline_usage(program);
    exit(0);
  }
  if(!cmd_param.ofile[0]) sprintf(cmd_param.ofile,"ssd_%s",cmd_param.ifile);
  
  
  *opts=cmd_param;
  return 0;
}

int bayesline_usage( const char *program)
{
  fprintf( stdout,"\n");
  fprintf( stdout, "REQUIRED:\n");
  fprintf( stdout, "      -i, --infile=FILE     filename for fourier domain data, expecting 3 columns: {f,Re,Im}\n");
  fprintf( stdout, "OPTIONAL:\n");
  fprintf( stdout, "     (-o, --outfile=FILE)   filename for output PSD model [default: ssd_infile] \n" );
  fprintf( stdout, "     (-v, --verbose)        thorough amount of screen output\n");
  fprintf( stdout, "     (-g, --gnuplot)        use gnuplot for progress reports\n");
  fprintf( stdout, "     (-l, --zerologL)       use logL=const to test detailed balance [use with --verbose] \n");
  fprintf( stdout, "     (-h, --help)           print this message and exit\n" );
  fprintf( stdout,"\n");
  return 0;
}
