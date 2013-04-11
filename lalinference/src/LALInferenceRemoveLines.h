

int LALInferenceRemoveLinesChiSquared(
    REAL8FrequencySeries        *spectrum,
    const REAL8TimeSeries       *tseries,
    UINT4                        seglen,
    UINT4                        stride,
    const REAL8Window           *window,
    const REAL8FFTPlan          *plan,
    REAL8                       *pvalues
    );

double chisqr(int Dof, double Cv);
double igf(double S, double Z);
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

