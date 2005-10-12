int BCVspin_metric(
   /*input*/
   int N,double *Sn,double fmin,double fmax,double beta,
   /*output*/
   double **bcv2metric,int dbg);

int BCVspin_spacing(
    double **metric3,
    double **a,
    double *deltax);

int BCVspin_effmetric(
    /* input */
    double **metric3,
    double **a,
    /* output */
    double **effmetric);
