

int LALInferenceRemoveLinesChiSquared(
    REAL8FrequencySeries        *spectrum,
    const REAL8TimeSeries       *tseries,
    UINT4                        seglen,
    UINT4                        stride,
    const REAL8Window           *window,
    const REAL8FFTPlan          *plan,
    REAL8                       *pvalues
    );

static void median_cleanup_REAL8( REAL8FrequencySeries *work, UINT4 n );
static int compare_REAL8( const void *p1, const void *p2 );
double chisqr(int Dof, double Cv);
static double igf(double S, double Z);
double gamma(double N);

int LALInferenceRemoveLinesKS(
    REAL8FrequencySeries        *spectrum,
    const REAL8TimeSeries       *tseries,
    UINT4                        seglen,
    UINT4                        stride,
    const REAL8Window           *window,
    const REAL8FFTPlan          *plan,
    REAL8                       *pvalues
    );

