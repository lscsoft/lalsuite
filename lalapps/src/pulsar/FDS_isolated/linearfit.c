
/* B. Krishnan, July 2005 */
/* small modification by M.A. Papa */
/* $Id$*/
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_cdf.h>


int main( int argc, char *argv[]){

  double *x=NULL; /* x values -- here h_0 values */
  double *y=NULL;  /* y values -- here Conf levels */
  double *w=NULL;  /* weights */
  FILE *fp=NULL;

  int n=0, r; 
  double temp1, temp2, temp3, temp4, temp5, temp6;
  double c0, c1, cov00, cov01, cov11, chisq;
  double h0, y0, y0_err;

  /* these variables are for producing plots */
  /*   int i; */
  /*   double  h0_fit, dh0; */
  /*   FILE *fpout=NULL; */


  /* open file for reading */
  fp = fopen( argv[1], "r");

  /* read the file to count lines */  
  do{
    r=fscanf(fp,"%lf%lf%lf%lf%lf\n", &temp1, &temp2, &temp3, &temp4, &temp5); 
    if (r==5) {
      n++;
    }
  } while (r != EOF);

  /* go back to start of file */
  rewind(fp);

  /* allocate sufficient memory for x, y and z */
  x = (double *)malloc(sizeof(double)*n);
  if (x==NULL) 
    {
      fprintf(stderr, "unable to allocate memory\n");
      return 1;
    }
  y = (double *)malloc(sizeof(double)*n);
    if (y==NULL) 
    {
      fprintf(stderr, "unable to allocate memory\n");
      return 1;
    }
  w = (double *)malloc(sizeof(double)* n);
  if (w==NULL) 
    {
      fprintf(stderr, "unable to allocate memory\n");
      return 1;
    }

  /* read file again and fill in values of vectors */  
  n = 0;
  do{
     r=fscanf(fp,"%lf%lf%lf%lf%lf\n", &temp1, &temp2, x+n, &temp4, y+n); 
    if (r==5) {
      temp6 = y[n];
      w[n] = temp1/(temp6*(1.0-temp6)); /* weight = 1/variance*/
      n++;
    }
  } while (r != EOF); 

  /* now we are done with the input file */  
  fclose(fp);
  
  /* call gsl fitting function */
  gsl_fit_wlinear(x, 1, w, 1, y, 1, n, &c0, &c1, &cov00, &cov01, &cov11, &chisq); 

  /* calculate h095 based on the fit */  
  h0 = (0.95 - c0)/c1;
  gsl_fit_linear_est(h0, c0, c1, cov00, cov01, cov11, &y0, &y0_err); 

  /* write output */
  fprintf(stdout, "%11.7E  %11.7E  %11.7E  %11.7E  ", h0*1.0e-19, c1*1.0e+19, y0_err, chisq );

  /*   fprintf(stdout, "best fit : Y = %g + %g X\n", c0, c1);   */
  /*   fprintf(stdout, "covariance matrix:\n");   */
  /*   fprintf(stdout, "[ %g,  %g\n  %g,  %g]\n", cov00, cov01, cov01, cov11);    */
  /*   fprintf(stdout, "chisq = %g\n", chisq);    */

  /*  dh0 = 0.1 * h0;  */
  /* fpout = fopen("fitout","w");   */
  /*   for (i = 0; i< 101; i++)   */
  /*     {   */
  /*       h0_fit = (h0-dh0) + i * 0.02 * dh0;   */
  /*       gsl_fit_linear_est(h0_fit, c0, c1, cov00, cov01, cov11, &y0, &y0_err);   */
  /*       fprintf(fpout, "%g  %g  %g\n", h0_fit, y0, y0_err);   */
  /*     }   */
  
  /*  fclose(fpout);  */

  /* calculate the error bounds on h0^95 */
  {  
    double h0_err_lo, h0_err_hi;  
    h0_err_lo = h0;  
    h0_err_hi = h0;  
    do{  
      gsl_fit_linear_est(h0_err_hi, c0, c1, cov00, cov01, cov11, &y0, &y0_err);   
      h0_err_hi += 0.0001*h0;  
    } while (y0 - y0_err < 0.95);  
    fprintf(stdout, "%11.7E  ", (h0_err_hi-h0)*1.0E-19);  
    do{  
      gsl_fit_linear_est(h0_err_lo, c0, c1, cov00, cov01, cov11, &y0, &y0_err);   
      h0_err_lo -= 0.0001*h0;  
    } while (y0 + y0_err > 0.95);  
    fprintf(stdout, "%11.7E  ", (h0 - h0_err_lo)*1.0E-19);  
    fprintf(stdout, "%g  ", 100.0*(h0_err_hi - h0_err_lo)/h0);  
  }  


  /* calculate chi-square quantities to test if the fit is acceptable */
  /* distribution is chi square with n-2 degrees of freedom */
  /* mean of chi square with k degrees of freedom is k and variance is 2k */
  fprintf(stdout , "%d  %d  ", n-2, 2*(n-2));

  /* also print area under chi-square distribution from chisq value */
  fprintf(stdout, "%g\n", gsl_cdf_chisq_Q( chisq, n-2));
  
  free(x); 
  free(y); 
  free(w);

  return 0;
}

