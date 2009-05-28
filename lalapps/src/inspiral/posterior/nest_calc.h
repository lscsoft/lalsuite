/* Nested sampling algorithm definitions */
/* (C) John Veitch 2009 */
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>

extern gsl_matrix *cov_mat;

extern CHAR outfile[512];
extern double etawindow;

double logadd(double a,double b);

extern INT4 seed;

void NestInit2PN(LALMCMCParameter *parameter, void *iT);

REAL8 mean(REAL8 *array,int N);

REAL8 nestZ(INT4 Nruns, INT4 Nlive, LALMCMCParameter **Live,LALMCMCInput *MCMCinput);

REAL4 MCMCSampleLimitedPrior(LALMCMCParameter *sample, LALMCMCParameter *temp, LALMCMCInput *MCMCInput,REAL8 minL,gsl_matrix *covM,INT4 N);

void calcCVM(gsl_matrix *cvm, LALMCMCParameter **samples,int N);

REAL8 ang_dist(REAL8 a1, REAL8 a2);

REAL8 ang_var(LALMCMCParameter **list, const char *pname, int N);

void fprintSample(FILE *fp,LALMCMCParameter *sample);

REAL8 sample_logt(int Nlive);

void Inject2PN(LALMCMCParameter *parameter, LALMCMCInput *MCMCinput, double SNR);
