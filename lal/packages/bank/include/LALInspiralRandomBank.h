typedef struct
tagMCBankIn
{
	int   verbose; /* if !=0 print helpful information */
	int   dim;     /* dimension of the parameer space */
	int   nIni;    /* initial number of seed points */
	int   nFin;    /* number of points remaining after filtering */
	float MM;      /* Maximum mismatch */
	float *pMax;   /* max values in the range of parameters */
	float *pMin;   /* min values in the range of parameters */
	float *p;      /* coordinates of an arbitrary pt on the parameter space */
	int   seed;    /* random seed for random bank */
	float error;   /* error condition */
	double *Sn;    /* power spectral density of the detector noise background */
	int   Npsd;    /* length of the spectral density array, normally N/2+1, where */
                       /* N is the number of samples */
} MCBankIn;

typedef void (MetricFunc) (
		MCBankIn in,
		int *diag,
		float *gij);


void EstimateNumberOfTemplates(
		MetricFunc *metric,
		MCBankIn *in);

void CreateRandomBank(
		MCBankIn *in,
		float *bank);

void FilterRandomBank(
		MetricFunc *metric,
		MCBankIn *in,
		float *bank);

void MonteCarloBank(
		MetricFunc *metric,
		MCBankIn *in,
		float *bank);

void BankStats(
		MetricFunc *metric,
		MCBankIn *in,
		float *bank);

void MCComputeDistance(
		int dim,
		int diag,
		float *x,
		float *y,
		float *gij,
		float *dist);

void TwoSphereMetric(
		MCBankIn in,
		int *diag,
		float *gij);

void FlatSpaceMetric(
		MCBankIn in,
		int *diag,
		float *gij);

void TwoDFlatSpacePolarMetric(
		MCBankIn in,
		int *diag,
		float *gij);

void BCVSpinMetric(
		MCBankIn in,
		int *diag,
		float *gij);

void SchwarzschildMetric(
		MCBankIn in,
		int *diag,
		float *gij);
