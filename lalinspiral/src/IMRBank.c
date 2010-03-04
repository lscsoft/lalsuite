#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/LALInspiralBank.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/IMRBank.h>
/* This function sets the pointers of the cumulative noise moment arrays to
 * zero */
int tmpltcnt = 0;

static REAL8 eta(REAL8 m1, REAL8 m2)
  {
  return m1*m2/(m1+m2)/(m1+m2);
  }

static REAL8 chirpmass(REAL8 m1, REAL8 m2)
  {
  return pow(m1*m2,0.6)/pow(m1+m2,0.2) / LAL_MTSUN_SI;
  }

static int XLALCreateIMRBankCumulativeNoiseMoments(IMRBankCumulativeNoiseMoments *moments, REAL8FrequencySeries *psd, REAL8 flow)
  {
  UINT4 length = psd->data->length;
  REAL8 deltaF = psd->deltaF;
  UINT4 i = 0;

  moments->psd = psd;
  moments->length = length;
  moments->deltaF = deltaF;
  moments->epoch = psd->epoch;
  moments->f0 = psd->f0;
  moments->sampleUnits = psd->sampleUnits;
  moments->flow = flow;
  /* Don't yet allocate memory for the moments */
  for (i=0; i < 25; i++)
    {
    moments->minus3[i] = NULL;
    moments->plus3[i] = NULL;
    moments->logminus3[i] = NULL;
    moments->logplus3[i] = NULL;
    moments->logsqplus3[i] = NULL;
    moments->logsqminus3[i] = NULL;
    }
  return 0;
  }

/* This function computes the negative noise moment integrals */
static int computeMinusMoment(IMRBankCumulativeNoiseMoments *moments, INT4 power, UINT4 logFlag)
  {
  REAL8FrequencySeries *psd = moments->psd;
  UINT4 length = psd->data->length;
  REAL8 flow = moments->flow;
  REAL8 deltaF = psd->deltaF / flow;
  UINT4 j = 0;
  REAL8 f;
  if (logFlag != 0 || power >= 0) return 1;
  if (moments->minus3[-power]) return 0;
  if ( !(moments->minus3[-power] = XLALCreateREAL8Vector(length)) ) return 1;
  moments->minus3[-power]->data[0] = 0;
  for (j = 1; j < length; j++)
    {
    if (psd->data->data[j] != 0)
      {
      f = deltaF*j;
      if (deltaF*j < 1)  moments->minus3[-power]->data[j] = 0;
      else moments->minus3[-power]->data[j] = moments->minus3[-power]->data[j-1]
        + 1./pow(f,-power / 3.0) / (psd->data->data[j]) * deltaF;
        /* pow(flow,-power/3.0);*/
      }
    else
      {
      moments->minus3[-power]->data[j]=moments->minus3[-power]->data[j-1];
      }
    }
  return 0;
  }

/* This function computes the positive noise moment integrals */

static int computePlusMoment(IMRBankCumulativeNoiseMoments *moments, INT4 power, UINT4 logFlag)
  {

   REAL8FrequencySeries *psd = moments->psd;
   UINT4 length = psd->data->length;
   REAL8 flow = moments->flow;
   REAL8 deltaF = psd->deltaF / flow;
   UINT4 j = 0;
   REAL8 f;
   if (logFlag != 0 || power < 0) return 1;
   if (moments->plus3[power]) return 0;
   if ( !(moments->plus3[power] = XLALCreateREAL8Vector(length)) ) return 1;
   moments->plus3[power]->data[0] = 0;

  for (j = 1; j < length; j++)
    {
    if (psd->data->data[j] != 0)
      {
      f = deltaF*j;
      if (deltaF*j < 1)  moments->plus3[power]->data[j] = 0;
      else moments->plus3[power]->data[j] = moments->plus3[power]->data[j-1]
        + pow(f,power / 3.0) / (psd->data->data[j]) * deltaF;
        /* pow(flow,power/3.0);*/
      }
    else
      {
      moments->plus3[power]->data[j]=moments->plus3[power]->data[j-1];
      }
    }
  return 0;
  }

/* This function computes the positive log noise moment integrals */

static int computePlusLogMoment(IMRBankCumulativeNoiseMoments *moments, INT4 power, UINT4 logFlag)
  {

   REAL8FrequencySeries *psd = moments->psd;
   UINT4 length = psd->data->length;
   REAL8 flow = moments->flow;
   REAL8 deltaF = psd->deltaF/flow;
   UINT4 j = 0;
   REAL8 f;
   if (logFlag != 1 || power < 0) return 1;
   if (moments->logplus3[power]) return 0;
   if ( !(moments->logplus3[power] = XLALCreateREAL8Vector(length)) ) return 1;
   moments->logplus3[power]->data[0] = 0;

  for (j = 1; j < length; j++)
    {
    if (psd->data->data[j] != 0)
      {
      f = deltaF*j;
      if (deltaF*j < 1)  moments->logplus3[power]->data[j] = 0;
      else moments->logplus3[power]->data[j] =
           moments->logplus3[power]->data[j-1]
	   + pow(f,power / 3.0)  * log(f)
        / (psd->data->data[j]) * deltaF ;
        /* pow(flow,power/3.0);*/
      }
    else
      {
      moments->logplus3[power]->data[j]=moments->logplus3[power]->data[j-1];
      }
    }
  return 0;
  }

/* This function computes the negative log noise moment integrals */

static int computeMinusLogMoment(IMRBankCumulativeNoiseMoments *moments, INT4 power, UINT4 logFlag)
  {

   REAL8FrequencySeries *psd = moments->psd;
   UINT4 length = psd->data->length;
   REAL8 flow = moments->flow;
   REAL8 deltaF = psd->deltaF/flow;
   UINT4 j = 0;
   REAL8 f;
   if (logFlag != 1 || power >= 0) return 1;
   if (moments->logminus3[-power]) return 0;
   if ( !(moments->logminus3[-power] = XLALCreateREAL8Vector(length)) ) return 1;
   moments->logminus3[-power]->data[0] = 0;

  for (j = 1; j < length; j++)
    {
    if (psd->data->data[j] != 0)
      {
      f = deltaF*j;
      if (deltaF*j < 1)  moments->logminus3[-power]->data[j] = 0;
      else moments->logminus3[-power]->data[j] =
            moments->logminus3[-power]->data[j-1]
            + 1./pow(f,-power / 3.0)  * log(f)
            / (psd->data->data[j]) * deltaF ;
            /* pow(flow,-power/3.0);*/
      }
    else
      {
      moments->logminus3[-power]->data[j]=moments->logminus3[-power]->data[j-1];
      }
    }
  return 0;
  }

/* This function computes the negative log squared noise moment integrals */

static int computeMinusLogsqMoment(IMRBankCumulativeNoiseMoments *moments, INT4 power, UINT4 logFlag)
  {

   REAL8FrequencySeries *psd = moments->psd;
   UINT4 length = psd->data->length;
   REAL8 flow = moments->flow;
   REAL8 deltaF = psd->deltaF/flow;
   UINT4 j = 0;
   REAL8 f;
   if (logFlag != 2 || power >= 0) return 1;
   if (moments->logsqminus3[-power]) return 0;
   if ( !(moments->logsqminus3[-power] = XLALCreateREAL8Vector(length)) ) return 1;
   moments->logsqminus3[-power]->data[0] = 0;

  for (j = 1; j < length; j++)
    {
    if (psd->data->data[j] != 0)
      {
      f = deltaF*j;
      if (deltaF*j < 1)  moments->logsqminus3[-power]->data[j] = 0;
      else
      moments->logsqminus3[-power]->data[j]
        = moments->logsqminus3[-power]->data[j-1]
        + 1./pow(f,-power / 3.0)  * log(f) * log(f)
        / (psd->data->data[j]) * deltaF ;
        /* pow(flow,-power/3.0);*/
      }
    else
      {
      moments->logsqminus3[-power]->data[j]=moments->logsqminus3[-power]->data[j-1];
      }
    }
  return 0;
  }

/* This function computes the positive log squared noise moment integrals */

static int computePlusLogsqMoment(IMRBankCumulativeNoiseMoments *moments, INT4 power, UINT4 logFlag)
  {

   REAL8FrequencySeries *psd = moments->psd;
   UINT4 length = psd->data->length;
   REAL8 flow = moments->flow;
   REAL8 deltaF = psd->deltaF/flow;
   UINT4 j = 0;
   REAL8 f;
   if (logFlag != 2 || power < 0) return 1;
   if (moments->logsqplus3[power]) return 0;
   if ( !(moments->logsqplus3[power] = XLALCreateREAL8Vector(length)) ) return 1;
   moments->logsqplus3[power]->data[0] = 0;

  for (j = 1; j < length; j++)
    {
    if (psd->data->data[j] != 0)
      {
      f = deltaF*j;
      if (deltaF*j < 1)  moments->logsqplus3[power]->data[j] = 0;
      else
      moments->logsqplus3[power]->data[j]
        = moments->logsqplus3[power]->data[j-1]
        + pow(f,power / 3.0)  * log(f) * log(f)
        / (psd->data->data[j]) * deltaF ;
        /* pow(flow,power/3.0);*/
      }
    else
      {
      moments->logsqplus3[power]->data[j]=moments->logsqplus3[power]->data[j-1];
      }
    }
  return 0;
  }

static int XLALComputeIMRBankCumulativeNoiseMoment(
                                     IMRBankCumulativeNoiseMoments *moments,
                                     INT4 power, UINT4 logFlag)
  {

  if (logFlag == 0 && power < 0) computeMinusMoment(moments,power,logFlag);
  if (logFlag == 0 && power >= 0) computePlusMoment(moments,power,logFlag);
  if (logFlag == 1 && power < 0) computeMinusLogMoment(moments,power,logFlag);
  if (logFlag == 1 && power >= 0) computePlusLogMoment(moments,power,logFlag);
  if (logFlag == 2 && power < 0) computeMinusLogsqMoment(moments,power,logFlag);
  if (logFlag == 2 && power >= 0) computePlusLogsqMoment(moments,power,logFlag);

  return 0;
  }




static int XLALDestroyIMRBankCumulativeNoiseMoments(IMRBankCumulativeNoiseMoments *moments)
  {
  UINT4 i = 0;
  for (i = 1; i < 25; i++)
    {
    if (moments->minus3[i]) XLALDestroyREAL8Vector(moments->minus3[i]);
    if (moments->plus3[i]) XLALDestroyREAL8Vector(moments->plus3[i]);
    if (moments->logminus3[i]) XLALDestroyREAL8Vector(moments->logminus3[i]);
    if (moments->logplus3[i]) XLALDestroyREAL8Vector(moments->logplus3[i]);
    if (moments->logsqminus3[i]) XLALDestroyREAL8Vector(moments->logsqminus3[i]);
    if (moments->logsqplus3[i]) XLALDestroyREAL8Vector(moments->logsqplus3[i]);
    }
  return 0;
  }

static REAL8 x(IMRBankCumulativeNoiseMoments *moments,
               REAL8 m, REAL8 flow, REAL8 fhigh, INT4 power, INT4 mpower)
  {
  /*printf("getting moment %d\n",power);*/
  REAL8FrequencySeries *psd = moments->psd;
  UINT4 fl = floor(flow / moments->deltaF);
  UINT4 fh = floor(fhigh / moments->deltaF);
  REAL8 newflow = flow;
  REAL8 output = 0;

  if (fh > psd->data->length)
    {
    fh = psd->data->length -1;
    fhigh = fh * moments->deltaF;
    }
  if (fl >= fh)
    {
    newflow = fhigh - 1.0;
    fl = floor(newflow / moments->deltaF);
    }

  XLALComputeIMRBankCumulativeNoiseMoment(moments, power,0);

  if (power < 0)
    {
    output
      = (moments->minus3[-power]->data[fh]-moments->minus3[-power]->data[fl]);
    }
  else
    {
    output
      = (moments->plus3[power]->data[fh]-moments->plus3[power]->data[fl]);
    }

  if (mpower < 0) output /= pow(m,-mpower/3.);
  if (mpower > 0) output *= pow(m,mpower/3.);
  /*printf("output of x function for m=%f fl=%f fh=%f pow=%d/3 is %e\n",
         m,flow,fhigh,power,output);*/
  /* FIXME: UH OH! Recursion! if the output is zero then bump up the frequency until it isn't */
  if (output <= 0) return 0.0; /*x(moments, m, (newflow-0.5), (fhigh+0.5), power, mpower);*/
  else return output;
  }

static REAL8 lx(IMRBankCumulativeNoiseMoments *moments,
                REAL8 m, REAL8 flow, REAL8 fhigh, INT4 power, INT4 mpower)
  {
  REAL8FrequencySeries *psd = moments->psd;
  REAL8 output = 0;
  UINT4 fl = floor(flow / moments->deltaF);
  UINT4 fh = floor(fhigh / moments->deltaF);
  REAL8 newflow = flow;

  if (fh > psd->data->length)
    {
    fh = psd->data->length -1;
    fhigh = fh * moments->deltaF;
    }
  if (fl >= fh)
    {
    newflow = fhigh - 1.0;
    fl = floor(newflow / moments->deltaF);
    }

  XLALComputeIMRBankCumulativeNoiseMoment(moments,power,1);

  if (power < 0)
    {
    /* FIXME check on this!!!!! Do I need log m too? */
    output = (moments->logminus3[-power]->data[fh] -  moments->logminus3[-power]->data[fl]);
    }
  else
    {
    /* FIXME check on this!!!!! Do I need log m too?*/
    output = (moments->logplus3[power]->data[fh] -  moments->logplus3[power]->data[fl] );
    }
  if (mpower < 0) output /= pow(m,-mpower/3.);
  if (mpower > 0) output *= pow(m,mpower/3.);
  /* FIXME: UH OH! Recursion! if the output is zero then bump up the frequency until it isn't */
  if (output <= 0) return 0.0;/* lx(moments, m, (newflow-0.5), (fhigh+0.5), power, mpower);*/
  else return output;
  }

static REAL8 lsqx(IMRBankCumulativeNoiseMoments *moments,
                REAL8 m, REAL8 flow, REAL8 fhigh, INT4 power, INT4 mpower)
  {
  REAL8FrequencySeries *psd = moments->psd;
  REAL8 output = 0;
  UINT4 fl = floor(flow / moments->deltaF);
  UINT4 fh = floor(fhigh / moments->deltaF);
  REAL8 newflow = flow;

  if (fh > psd->data->length)
    {
    fh = psd->data->length -1;
    fhigh = fh * moments->deltaF;
    }
  if (fl >= fh)
    {
    newflow = fhigh - 1.0;
    fl = floor(newflow / moments->deltaF);
    }

  XLALComputeIMRBankCumulativeNoiseMoment(moments,power,2);
  if (power < 0)
    {
    /* FIXME check on this!!!!! Do I need log m too? */
    output = (moments->logminus3[-power]->data[fh] -  moments->logminus3[-power]->data[fl]);
    }
  else
    {
    /* FIXME check on this!!!!! Do I need log m too?*/
    output = (moments->logplus3[power]->data[fh] -  moments->logplus3[power]->data[fl] );
    }
  if (mpower < 0) output /= pow(m,-mpower/3.);
  if (mpower > 0) output *= pow(m,mpower/3.);
  /* FIXME: UH OH! Recursion! if the output is zero then bump up the frequency until it isn't */
  if (output <= 0) return 0.0; /*lsqx(moments, m, (newflow-0.5), (fhigh+0.5), power, mpower);*/
  else return output;
  }

static REAL8 JPsiEta(REAL8 mass1, REAL8 mass2,
         REAL8 fl, REAL8 fh, INT4 xpow, IMRBankCumulativeNoiseMoments *I)
  {
  REAL8 output = 0;
  REAL8 n = eta(mass1,mass2);
  REAL8 m = (mass1+mass2);
  /*DOUBLE CHECKED */
  output = 1.0 / n / n * ( 0.0
    - 0.0366605 * x(I,m,fl,fh,-3+xpow,-3)
    - 3.76385 * x(I,m,fl,fh,0+xpow,0)

    - 0.00347799 * x(I,m,fl,fh,-5+xpow,-5)
    + 0.549222 * x(I,m,fl,fh,-2+xpow,-2)
    - 0.481733 * x(I,m,fl,fh,-1+xpow,-1)
    + 9.55141 * x(I,m,fl,fh,1+xpow,1)
    - 47.9368 * x(I,m,fl,fh,2+xpow,2)
    + 0.685673 * n * n * x(I,m,fl,fh,-1+xpow,-1)
    + 1.51082 * n * n * x(I,m,fl,fh,1+xpow,1)
    - 15.4692 * n * n * x(I,m,fl,fh,2+xpow,2)
    - 6.77125 * n * n * n * x(I,m,fl,fh,1+xpow,1)

    + 11.2916/3.0 * log(1./m/LAL_PI) * x(I,m,fl,fh,0+xpow,0)
    - 11.2916/3.0 * lx(I,m,fl,fh,0+xpow,0)

    + 288.855/3.0 * log(1./m/LAL_PI) * x(I,m,fl,fh,3+xpow,3)
    - 288.855/3.0 * lx(I,m,fl,fh,3+xpow,3)

    - 11.1937/3.0 * log(1./m/LAL_PI) * x(I,m,fl,fh,1+xpow,1)
    + 11.1937/3.0 * lx(I,m,fl,fh,1+xpow,1)
    );
  return output; /*/I->flow/2./LAL_PI;*/
  /*return output /2./LAL_PI/I->flow; */
  }


static REAL8 JPsiEtaPsiEta(REAL8 mass1, REAL8 mass2,
         REAL8 fl, REAL8 fh, INT4 xpow, IMRBankCumulativeNoiseMoments *I)
  {
  REAL8 output = 0;
  REAL8 n = eta(mass1,mass2);
  REAL8 m = (mass1+mass2);
  /*DOUBLE CHECKED*/
  output = 1.0 / n / n / n / n * ( 0.0

    + 83437.0/9.0 * ( log(1./LAL_PI/m) *log(1./LAL_PI/m)* x(I,m,fl,fh,6+xpow,6)
                    - 2.0 * log(1./LAL_PI/m) * lx(I,m,fl,fh,6+xpow,6)
                    + lsqx(I,m,fl,fh,6+xpow,6) )

    + 9.0/16384.0 / pow(LAL_PI,10.0/3.0) * x(I,m,fl,fh,-10+xpow,-10)
    + 3715./688128. / pow(LAL_PI,8./3.) * x(I,m,fl,fh,-8+xpow,-8)
    - 12096./688128. / pow(LAL_PI,4./3.) * x(I,m,fl,fh,-7+xpow,-7)
    + 0.046337 / LAL_PI/LAL_PI * x(I,m,fl,fh,-6+xpow,-6)
    - 0.139045 / LAL_PI/LAL_PI * x(I,m,fl,fh,-5+xpow,-5)
    + 2.66999 / LAL_PI/LAL_PI * x(I,m,fl,fh,-4+xpow,-4)
    - 0.0470734 * n*n / LAL_PI/LAL_PI * x(I,m,fl,fh,-6+xpow,-6)
    - 0.599909 * n*n / LAL_PI/LAL_PI * x(I,m,fl,fh, -4+xpow,-4)
    + 0.464865 * n*n*n / LAL_PI/LAL_PI * x(I,m,fl,fh,-4+xpow,-4)


    - 0.775197 / 3.0/LAL_PI/LAL_PI * log(1./m/LAL_PI) * x(I,m,fl,fh,-5+xpow,-5)
    + 0.775197 / 3.0/LAL_PI/LAL_PI * lx(I,m,fl,fh,-5+xpow,-5)

    + 0.768476 / 3.0/LAL_PI/LAL_PI * log(1./m/LAL_PI) * x(I,m,fl,fh,-4+xpow,-4)
    - 0.768476 / 3.0/LAL_PI/LAL_PI *  lx(I,m,fl,fh,-4+xpow,-4)


    + 0.0802598 * x(I,m,fl,fh,-3+xpow,-3)
    - 4.60263 * x(I,m,fl,fh,-2+xpow,-2)
    + 17.6328 * x(I,m,fl,fh,-1+xpow,-1)
    + 0.860777 * n*n * x(I,m,fl,fh,-3+xpow,-3)
    - 0.771398 * n*n * x(I,m,fl,fh,-2+xpow,-2)
    - 2.36778 * n*n * x(I,m,fl,fh,-1+xpow,-1)
    + 0.496475 * n*n*n * x(I,m,fl,fh,-3+xpow,-3)
    - 7.43784 * n*n*n * x(I,m,fl,fh,-2+xpow,-2)
    + 0.470148 * n*n*n*n * x(I,m,fl,fh,-1+xpow,-1)

    - 0.827909 * log(1./m/LAL_PI) * x(I,m,fl,fh,-3+xpow,-3)
    + 0.827909 * lx(I,m,fl,fh,-3+xpow,-3)

    + 11.2146 * log(1./m/LAL_PI) * x(I,m,fl,fh,-2+xpow,-2)
    - 11.2146 * lx(I,m,fl,fh,-2+xpow,-2)

    - 23.1737 * log(1./m/LAL_PI) * x(I,m,fl,fh,-1+xpow,-1)
    + 23.1737 * lx(I,m,fl,fh,-1+xpow,-1)

    + 15.4846 * n*n * log(1./m/LAL_PI) * x(I,m,fl,fh,-1+xpow,-1)
    - 15.4849 * n*n * lx(I,m,fl,fh,-1+xpow,-1)


    - 915.728 * x(I,m,fl,fh,3+xpow,3)
    + 2297.94 * x(I,m,fl,fh,4+xpow,4)
    - 440.353 * n*n * x(I,m,fl,fh,3+xpow,3)
    + 1483.009 * n*n * x(I,m,fl,fh,4+xpow,4)
    + 649.184 * n*n*n * x(I,m,fl,fh,3+xpow,3)
    - 46.7423 * n*n*n*n * x(I,m,fl,fh,3+xpow,3)
    + 239.297 * n*n*n*n * x(I,m,fl,fh,4+xpow,4)
    + 209.492 * n*n*n*n*n * x(I,m,fl,fh,3+xpow,3)


    - 1101.24/3.0 * log(1./m/LAL_PI) * x(I,m,fl,fh,3+xpow,3)
    + 1101.24/3.0 * lx(I,m,fl,fh,3+xpow,3)

    + 5517.94/3. * log(1./m/LAL_PI) * x(I,m,fl,fh,4+xpow,4)
    - 5517.94/3. * lx(I,m,fl,fh,4+xpow,4)

    - 27963.5/3. * log(1./m/LAL_PI) * x(I,m,fl,fh,5+xpow,5)
    + 27963.5/3. * lx(I,m,fl,fh,5+xpow,5)

    + 346.315/3. *n*n * log(1./m/LAL_PI) * x(I,m,fl,fh,3+xpow,3)
    - 346.315/3. *n*n * lx(I,m,fl,fh,3+xpow,3)

    + 872.813/3. *n*n * log(1./m/LAL_PI) * x(I,m,fl,fh,4+xpow,4)
    - 872.813/3. *n*n * lx(I,m,fl,fh,4+xpow,4)

    - 8936.71/3. *n*n * log(1./m/LAL_PI) * x(I,m,fl,fh,5+xpow,5)
    + 8936.71/3. *n*n  * lx(I,m,fl,fh,5+xpow,5)

    - 3911.81/3. *n*n*n * log(1./m/LAL_PI) * x(I,m,fl,fh,4+xpow,4)
    + 3911.81/3. * lx(I,m,fl,fh,4+xpow,4)

    + 6523.24/9. * log(1./m/LAL_PI) * log(1./m/LAL_PI) * x(I,m,fl,fh,3+xpow,3)
    - 6523.24/9. * 2.0 * log(1./m/LAL_PI) * lx(I,m,fl,fh,3+xpow,3)
    + 6523.24/9. * lsqx(I,m,fl,fh,3+xpow,3)

    - 6466.69/9. * log(1./m/LAL_PI) * log(1./m/LAL_PI) * x(I,m,fl,fh,4+xpow,4)
    + 6466.69/9. * 2.0 * log(1./m/LAL_PI) * lx(I,m,fl,fh,4+xpow,4)
    - 6466.69/9. * lsqx(I,m,fl,fh,4+xpow,4)


    - 47.6918 * x(I,m,fl,fh,0+xpow,0)
    - 25.7147 * x(I,m,fl,fh,1+xpow,1)
    + 452.083 * x(I,m,fl,fh,2+xpow,2)
    - 5.3494 * n*n * x(I,m,fl,fh,0+xpow,0)
    - 62.2069 * n*n * x(I,m,fl,fh,1+xpow,1)
    + 145.309 * n*n * x(I,m,fl,fh, 2+xpow,2)
    + 6.52387 * n*n*n * x(I,m,fl,fh,0+xpow,0)
    + 50.972 * n*n*n * x(I,m,fl,fh,1+xpow,1)
    - 129.35 * n*n*n * x(I,m,fl,fh,2+xpow,2)
    + 2.07185 * n*n*n*n * x(I,m,fl,fh,0+xpow,0)
    - 21.2137 * n*n*n*n * x(I,m,fl,fh,1+xpow,1)
    + 2.28257 * n*n*n*n * x(I,m,fl,fh,2+xpow,2)
    - 9.28573 * n*n*n*n*n * x(I,m,fl,fh,0+xpow,0)
    + 45.8498 * n*n*n*n*n*n * x(I,m,fl,fh,2+xpow,2)

    - 95.3939/3.0 * log(1./m/LAL_PI) * x(I,m,fl,fh,0+xpow,0)
    + 95.3939/3.0 * lx(I,m,fl,fh,0+xpow,0)

    + 617.254/3.0 * log(1./m/LAL_PI) * x(I,m,fl,fh,1+xpow,1)
    - 617.254/3.0 * lx(I,m,fl,fh,1+xpow,1)

    - 1574.69/3.0 * log(1./m/LAL_PI) * x(I,m,fl,fh,2+xpow,2)
    + 1574.69/3.0 * lx(I,m,fl,fh,2+xpow,2)

    - 15.3504/3.0 * n*n * log(1./m/LAL_PI) * x(I,m,fl,fh,0+xpow,0)
    + 15.3504/3.0 * n*n * lx(I,m,fl,fh,0+xpow,0)

    + 34.119/3.0 * n*n * log(1./m/LAL_PI) * x(I,m,fl,fh,1+xpow,1)
    - 34.119/3.0 * n*n * lx(I,m,fl,fh,1+xpow,1)

    + 12.9537/3.0 * n*n * log(1./m/LAL_PI) * x(I,m,fl,fh,2+xpow,2)
    - 12.9537/3.0 * n*n * lx(I,m,fl,fh,2+xpow,2)

    - 152.916/3.0 * n*n*n * log(1./m/LAL_PI) * x(I,m,fl,fh,1+xpow,1)
    + 152.916/3.0 * n*n*n * lx(I,m,fl,fh,1+xpow,1)

    + 151.59/3.0 * n*n*n * log(1./m/LAL_PI) * x(I,m,fl,fh,2+xpow,2)
    - 151.59/3.0 * n*n*n * lx(I,m,fl,fh,2+xpow,2)

    + 127.499/9.0 * log(1./m/LAL_PI)*log(1./m/LAL_PI) * x(I,m,fl,fh,0+xpow,0)
    - 127.499/9.0 * 2.0 * log(1./m/LAL_PI) * lx(I,m,fl,fh,0+xpow,0)
    + 127.499/9.0 * lsqx(I,m,fl,fh,0+xpow,0)

    - 252.788/9.0 * log(1./m/LAL_PI)*log(1./m/LAL_PI) * x(I,m,fl,fh,1+xpow,1)
    + 252.788/9.0 * 2.0 * log(1./m/LAL_PI) * lx(I,m,fl,fh,1+xpow,1)
    - 252.788/9.0 * lsqx(I,m,fl,fh,1+xpow,1)

    + 125.298/9.0 * log(1./m/LAL_PI)*log(1./m/LAL_PI) * x(I,m,fl,fh,2+xpow,2)
    - 125.298/9.0 * 2.0 * log(1./m/LAL_PI) * lx(I,m,fl,fh,2+xpow,2)
    + 125.298/9.0 * lsqx(I,m,fl,fh,2+xpow,2)

    );
  return output; /*/I->flow/I->flow/2./LAL_PI/2./LAL_PI;*/
  /*return output /2./LAL_PI/2./LAL_PI/I->flow/I->flow; */
  }



static REAL8 JPsiM(REAL8 mass1, REAL8 mass2,
         REAL8 fl, REAL8 fh, INT4 xpow, IMRBankCumulativeNoiseMoments *I)
  {
  REAL8 output = 0;
  REAL8 n = eta(mass1,mass2);
  REAL8 m = (mass1+mass2);
  /*printf("n = %e, m = %e \n",n,m);*/
  /* DOUBLE CHECKED */
  output = -1.0 / n / m * ( 0.0
    + 0.0366605 * x(I,m,fl,fh,-3+xpow,-3)
    - 3.76385 * x(I,m,fl,fh,0+xpow,0)
    - 96.2849 * x(I,m,fl,fh,+3+xpow,+3)
    + 0.00579665 * x(I,m,fl,fh,-5+xpow,-5)
    - 0.366148 * x(I,m,fl,fh,-2+xpow,-2)
    + 0.160578 * x(I,m,fl,fh,-1+xpow,-1)
    + 6.91502 * x(I,m,fl,fh,+1+xpow,+1)
    - 31.9579 * x(I,m,fl,fh,+2+xpow,+2)
    + 0.0455913 * n * x(I,m,fl,fh,-3+xpow,-3)
    + 0.53178 * n * x(I,m,fl,fh,0+xpow,0)
    + 0.287398 * n * x(I,m,fl,fh,-1+xpow,-1)
    + 37.8542 * n * x(I,m,fl,fh,1+xpow,+1)
    - 26.3593 * n * x(I,m,fl,fh,+2+xpow,+2)
    + 0.228558 * n * n * x(I,m,fl,fh,-1+xpow,-1)
    - 0.503606 * n * n * x(I,m,fl,fh,+1+xpow,+1)
    + 10.3128 * n * n * x(I,m,fl,fh,+2+xpow,+2)
    + 1.12854 * n * n * n * x(I,m,fl,fh,+1+xpow,+1)

    + 288.855/3. * log(1./m/LAL_PI) * x(I,m,fl,fh,+3+xpow,+3)
    - 288.855/3. * lx(I,m,fl,fh,+3+xpow,+3)
    - 3.73122/3. * log(1./m/LAL_PI) * x(I,m,fl,fh,+1+xpow,+1)
    + 3.73122 / 3.0 * lx(I,m,fl,fh,+1+xpow,+1)

  );

  /* Lets use solar masses here and stick with it */
  return output;
  /*return output  /2./LAL_PI/I->flow; */
  }

static REAL8 JPsiMPsiM(REAL8 mass1, REAL8 mass2,
         REAL8 fl, REAL8 fh, INT4 xpow, IMRBankCumulativeNoiseMoments *I)
  {
  REAL8 output = 0;
  REAL8 n = eta(mass1,mass2);
  REAL8 m = (mass1+mass2);
  /*DOUBLE CHECKED */
  output = 1.0 / m / m / n / n * (0.0
    + 42644.4 * x(I,m,fl,fh,6+xpow,6)
    + 0.0000336011 * x(I,m,fl,fh,-10+xpow,-10)
    + 0.000425016 * x(I,m,fl,fh,-8+xpow,-8)
    - 0.00424486 * x(I,m,fl,fh,-7+xpow,-7)
    + 0.000528553 * n * x(I,m,fl,fh, -8+xpow,-8)
    + 0.00320562 * x(I,m,fl,fh,-6+xpow,-6)
    - 0.0704818 * x(I,m,fl,fh,-5+xpow,-5)
    + 0.242512 * x(I,m,fl,fh,-4+xpow,-4)
    + 0.00667353 * n * x(I,m,fl,fh,-6+xpow,-6)
    - 0.0272212 * n * x(I,m,fl,fh,-5+xpow,-5)
    + 0.474562 * n * x(I,m,fl,fh,-4+xpow,-4)
    + 0.0047283 * n*n * x(I,m,fl,fh,-6+xpow,-6)
    + 0.0371162 * n*n * x(I,m,fl,fh,-4+xpow,-4)
    + 0.033924 * n*n*n * x(I,m,fl,fh,-4+xpow,-4)
    - 0.764057 * x(I,m,fl,fh,-3+xpow,-3)
    + 0.999372 * x(I,m,fl,fh,-2+xpow,-2)
    - 9.65841 * x(I,m,fl,fh,-1+xpow,-1)
    - 0.820186 * n * x(I,m,fl,fh,-3+xpow,-3)
    + 3.23871 * n * x(I,m,fl,fh,-2+xpow,-2)
    - 34.5591 * n * x(I,m,fl,fh,-1+xpow,-1)
    + 0.00067654 * n*n * x(I,m,fl,fh,-3+xpow,-3)
    + 3.57066 * n*n * x(I,m,fl,fh,-2+xpow,-2)
    - 2.69353 * n*n * x(I,m,fl,fh,-1+xpow,-1)
    + 0.168154 * n*n*n * x(I,m,fl,fh,-2+xpow,-2)
    + 0.357007 * n*n*n * x(I,m,fl,fh, -1+xpow,-1)
    + 0.155142 * n*n*n*n * x(I,m,fl,fh,-2+xpow,-2)
    + 1021.53 * x(I,m,fl,fh,3+xpow,3)
    - 2422.69 * x(I,m,fl,fh,4+xpow,4)
    + 13198.9 * x(I,m,fl,fh,5+xpow,5)
    - 3078.72 * n * x(I,m,fl,fh,3+xpow,3)
    - 13949.4 * n * x(I,m,fl,fh,4+xpow,4)
    + 10886.7 * n * x(I,m,fl,fh,5+xpow,5)
    - 1791.44 * n*n * x(I,m,fl,fh,3+xpow,3)
    + 243.658 * n*n * x(I,m,fl,fh,4+xpow,4)
    - 4259.3 * n*n * x(I,m,fl,fh,5+xpow,5)
    + 735.185 * n*n*n * x(I,m,fl,fh,3+xpow,3)
    - 1009.78 * n*n*n * x(I,m,fl,fh,4+xpow,4)
    - 69.8824 * n*n*n*n * x(I,m,fl,fh,3+xpow,3)
    + 106.354 * n*n*n*n * x(I,m,fl,fh, 4+xpow,4)
    + 23.2769 * n*n*n*n*n * x(I,m,fl,fh,3+xpow,3)
    + 25.1061 * x(I,m,fl,fh,0+xpow,0)
    + 78.1876 * x(I,m,fl,fh,1+xpow,1)
    + 243.784 * x(I,m,fl,fh,2+xpow,2)
    + 13.4186 * n * x(I,m,fl,fh,0+xpow,0)
    - 302.915 * n * x(I,m,fl,fh,1+xpow,1)
    + 677.095 * n * x(I,m,fl,fh,2+xpow,2)
    + 18.1317 * n*n * x(I,m,fl,fh,0+xpow,0)
    + 17.6088 * n*n * x(I,m,fl,fh,1+xpow,1)
    + 1224.48 * n*n * x(I,m,fl,fh,2+xpow,2)
    + 17.3768 * n*n*n * x(I,m,fl,fh,0+xpow,0)
    - 15.1545 * n*n*n * x(I,m,fl,fh,1+xpow,1)
    - 8.33761 * n*n*n * x(I,m,fl,fh,2+xpow,2)
    + 0.418249 * n*n*n*n * x(I,m,fl,fh,0+xpow,0)
    + 5.914442 * n*n*n*n * x(I,m,fl,fh,1+xpow,1)
    + 85.6937 * n*n*n*n * x(I,m,fl,fh,2+xpow,2)
    + 0.515874 * n*n*n*n*n * x(I,m,fl,fh,0+xpow,0)
    - 1.13668 * n*n*n*n*n * x(I,m,fl,fh,2+xpow,2)
    + 1.27361 * n*n*n*n*n*n * x(I,m,fl,fh,2+xpow,2)

/*This was a problem */
    - 119300./3. * log(1./m/LAL_PI) * x(I,m,fl,fh,6+xpow,6)
    + 119300./3. * lx(I,m,fl,fh,6+xpow,6)

    - 0.0432572/3. * log(1./m/LAL_PI) * x(I,m,fl,fh,-4+xpow,-4)
    + 0.0432572/3. * lx(I,m,fl,fh,-4+xpow,-4)

    + 3.0752/3. * log(1./m/LAL_PI) * x(I,m,fl,fh,-2+xpow,-2)
    - 3.0752/3. * lx(I,m,fl,fh,-2+xpow,-2)

    + 2.73236/3. * log(1./m/LAL_PI) * x(I,m,fl,fh,-1+xpow,-1)
    - 2.73236/3. * lx(I,m,fl,fh,-1+xpow,-1)

    - 0.34022/3. * log(1./m/LAL_PI) * n * x(I,m,fl,fh,-2+xpow,-2)
    + 0.34022/3. * lx(I,m,fl,fh,-2+xpow,-2)

    - 1935.93/3. * log(1./m/LAL_PI) * x(I,m,fl,fh,3+xpow,3)
    + 1935.93/3. * lx(I,m,fl,fh,3+xpow,3)

    + 6358.42/3. * log(1./m/LAL_PI) * x(I,m,fl,fh,4+xpow,4)
    - 6358.42/3. * lx(I,m,fl,fh,4+xpow,4)

    - 18462.4/3. * log(1./m/LAL_PI) * x(I,m,fl,fh,5+xpow,5)
    + 18462.4/3. * lx(I,m,fl,fh,5+xpow,5)

    + 503.919/3. * n * log(1./m/LAL_PI) * x(I,m,fl,fh,3+xpow,3)
    - 503.919/3. * n * lx(I,m,fl,fh,3+xpow,3)

    + 21868.7/3. * n * log(1./m/LAL_PI) * x(I,m,fl,fh,4+xpow,4)
    - 21868.7/3. * n * lx(I,m,fl,fh,4+xpow,4)

    - 15228./3.  * n * log(1./m/LAL_PI) * x(I,m,fl,fh,5+xpow,5)
    + 15228./3.  * n * lx(I,m,fl,fh,5+xpow,5)

    - 76.9588/3. * n*n * log(1./m/LAL_PI) * x(I,m,fl,fh,3+xpow,3)
    + 76.9588/3. * n*n * lx(I,m,fl,fh,3+xpow,3)

    - 290.938/3. * n*n * log(1./m/LAL_PI) * x(I,m,fl,fh,4+xpow,4)
    + 290.938/3. * n*n * lx(I,m,fl,fh,4+xpow,4)

    + 5957.81/3. * n*n * log(1./m/LAL_PI) * x(I,m,fl,fh,5+xpow,5)
    - 5957.81/3. * n*n * lx(I,m,fl,fh,5+xpow,5)

    + 651.969/3. *n*n*n * log(1./m/LAL_PI) * x(I,m,fl,fh,4+xpow,4)
    - 651.969/3. *n*n*n * lx(I,m,fl,fh,4+xpow,4)

    + 19.9808/3. * log(1./m/LAL_PI) * x(I,m,fl,fh,0+xpow,0)
    - 19.9808/3. * lx(I,m,fl,fh,0+xpow,0)

    - 183.44/3.  * log(1./m/LAL_PI) * x(I,m,fl,fh,1+xpow,1)
    + 183.44/3. * lx(I,m,fl,fh,1+xpow,1)

    + 30.5397/3. * log(1./m/LAL_PI) * x(I,m,fl,fh,2+xpow,2)
    - 30.5397/3. * lx(I,m,fl,fh,2+xpow,2)

    + 24.1946/3. * n * log(1./m/LAL_PI) * x(I,m,fl,fh,0+xpow,0)
    - 24.1946/3. * n * lx(I,m,fl,fh,0+xpow,0)

    - 3.96838/3. * n * log(1./m/LAL_PI) * x(I,m,fl,fh,1+xpow,1)
    + 3.96838/3. * n * lx(I,m,fl,fh,1+xpow,1)

    - 116.51/3. * n * log(1./m/LAL_PI) * x(I,m,fl,fh,2+xpow,2)
    + 116.51/3. * n * lx(I,m,fl,fh,2+xpow,2)

    - 1.7056/3. * n*n * log(1./m/LAL_PI) * x(I,m,fl,fh,0+xpow,0)
    + 1.7056/3. * n*n * lx(I,m,fl,fh,0+xpow,0)

    - 8.42168 * n*n*n * log(1./m/LAL_PI) * x(I,m,fl,fh,2+xpow,2)
    + 8.42168 * n*n*n * lx(I,m,fl,fh,2+xpow,2)

    + 83437.0/9.0 * log(1./m/LAL_PI)*log(1./m/LAL_PI) * x(I,m,fl,fh,6+xpow,6)
    - 83437.0/9.0 * 2.0 * log(1./m/LAL_PI) * lx(I,m,fl,fh,6+xpow,6)
    + 83437.0/9.0 * lsqx(I,m,fl,fh,6+xpow,6)

    + 13.922/9.0 * log(1./m/LAL_PI)*log(1./m/LAL_PI) * x(I,m,fl,fh,2+xpow,2)
    - 13.922/9.0 * 2.0 * log(1./m/LAL_PI) * lx(I,m,fl,fh,2+xpow,2)
    + 13.922/9.0 * lsqx(I,m,fl,fh,2+xpow,2)

    - 2155.56/9.0 * log(1./m/LAL_PI)*log(1./m/LAL_PI) * x(I,m,fl,fh,4+xpow,4)
    + 2155.56/9.0 * 2.0 * log(1./m/LAL_PI) * lx(I,m,fl,fh,4+xpow,4)
    - 2155.56/9.0 * lsqx(I,m,fl,fh,4+xpow,4)

  );
  return output;
  /*return output / 2./LAL_PI/2./LAL_PI/I->flow/I->flow; */

  }

/* FINISH CHECKING HERE DOWN! */
static REAL8 JPsiMPsiEta(REAL8 mass1, REAL8 mass2,
         REAL8 fl, REAL8 fh, INT4 xpow, IMRBankCumulativeNoiseMoments *I)
  {
  REAL8 output = 0;
  REAL8 n = eta(mass1,mass2);
  REAL8 m = (mass1+mass2);
  output = 1 / m * (0.0
    - 24.7587 * x(I,m,fl,fh,0+xpow,0)
    + 0.476432 * x(I,m,fl,fh,1+xpow,1)
    - 12.9204 * x(I,m,fl,fh,2+xpow,2)
    + 0.000340013 / n/n/n * x(I,m,fl,fh,-8+xpow,-8)
    - 0.0044571 / n/n/n * x(I,m,fl,fh,-7+xpow,-7)
    + 0.000158566 /n/n * x(I,m,fl,fh,-8+xpow,-8)
    + 15.899 /n/n/n * x(I,m,fl,fh,0+xpow,0)
    + 107.161 /n/n/n * x(I,m,fl,fh,1+xpow,1)
    + 0.0000201607 /n/n/n * x(I,m,fl,fh,-10+xpow,-10)
    - 413.144 /n/n/n * x(I,m,fl,fh,2+xpow,2)
    + 31.9702 /n/n * x(I,m,fl,fh,0+xpow,0)
    + 138.472 /n/n * x(I,m,fl,fh,1+xpow,1)
    - 435.282 /n/n * x(I,m,fl,fh,2+xpow,2)
    - 18.7378 /n * x(I,m,fl,fh,0+xpow,0)
    + 44.112 /n * x(I,m,fl,fh,1+xpow,1)
    + 40.9749 /n * x(I,m,fl,fh,2+xpow,2)
    + 1.94537 * n * x(I,m,fl,fh,0+xpow,0)
    + 0.0652036 * n * x(I,m,fl,fh,1+xpow,1)
    + 257.081 * n * x(I,m,fl,fh,2+xpow,2)
    + 0.773811 * n*n * x(I,m,fl,fh,0+xpow,0)
    - 5.11506 * n*n * x(I,m,fl,fh,2+xpow,2)
    + 7.64163 * n*n*n * x(I,m,fl,fh,2+xpow,2)

    + 10.5896/3. /n/n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,0+xpow,0)
    - 10.5896/3. /n/n/n * lx(I,m,fl,fh,0+xpow,0)

    - 158.645/3 /n/n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,1+xpow,1)
    + 158.645/3. /n/n/n * lx(I,m,fl,fh,1+xpow,1)

    + 139.151/3.  /n/n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,2+xpow,2)
    - 139.151/3.  /n/n/n * lx(I,m,fl,fh,2+xpow,2)

    - 198.06/3. /n * log(1./m/LAL_PI) * x(I,m,fl,fh,2+xpow,2)
    + 198.06/3. /n * lx(I,m,fl,fh,2+xpow,2)

    - 12.743/3. * log(1./m/LAL_PI) * x(I,m,fl,fh,1+xpow,1)
    + 12.743/3. * lx(I,m,fl,fh,1+xpow,1)

    - 12.6325/3. * log(1./m/LAL_PI) * x(I,m,fl,fh,2+xpow,2)
    + 12.6325/3. * lx(I,m,fl,fh,2+xpow,2)

    + 31.9102/3. /n/n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,0+xpow,0)
    - 31.9102/3. /n/n/n * lx(I,m,fl,fh,0+xpow,0)

    - 28.4929/3. /n/n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,1+xpow,1)
    + 28.4929/3. /n/n/n * lx(I,m,fl,fh,1+xpow,1)

    + 427.513/3. /n/n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,2+xpow,2)
    - 427.513/3. /n/n/n * lx(I,m,fl,fh,2+xpow,2)

    - 15.958/3. /n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,0+xpow,0)
    + 15.958/3. /n/n * lx(I,m,fl,fh,0+xpow,0)

    - 421.481/3. /n/n *  log(1./m/LAL_PI) * x(I,m,fl,fh,1+xpow,1)
    + 421.481/3. /n/n *  lx(I,m,fl,fh,1+xpow,1)

    + 638.378/3. /n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,2+xpow,2)
    - 638.378/3. /n/n * lx(I,m,fl,fh,2+xpow,2)

    + 5.1168/3. /n * log(1./m/LAL_PI) * x(I,m,fl,fh,0+xpow,0)
    - 5.1168/3. /n * lx(I,m,fl,fh,0+xpow,0)

    + 5.68649/3. /n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,1+xpow,1)
    - 5.68649/3. /n/n * lx(I,m,fl,fh,1+xpow,1)

    - 182.468/3. /n * log(1./m/LAL_PI) * x(I,m,fl,fh,2+xpow,2)
    + 182.468/3. /n * lx(I,m,fl,fh,2+xpow,2)

    + 42.1313/9. /n/n/n*log(1./m/LAL_PI)*log(1/m/LAL_PI) * x(I,m,fl,fh,1+xpow,1)
    - 42.1313/9. /n/n/n*log(1./m/LAL_PI)*lx(I,m,fl,fh,1+xpow,1)
    + 42.1313/9. /n/n/n*lsqx(I,m,fl,fh,1+xpow,1)

    - 41.7661/9. /n/n/n*log(1./m/LAL_PI)*log(1/m/LAL_PI) * x(I,m,fl,fh,2+xpow,2)
    + 41.7661/9. /n/n/n*log(1./m/LAL_PI)*lx(I,m,fl,fh,2+xpow,2)
    - 41.7661/9. /n/n/n*lsqx(I,m,fl,fh,2+xpow,2)

    + 0.0119149 * x(I,m,fl,fh,-4+xpow,-4)
    + 0.004694922 /n/n/n * x(I,m,fl,fh,-6+xpow,-6)
    - 0.0248308 /n/n/n * x(I,m,fl,fh,-5+xpow,-5)
    + 0.193328 /n/n/n * x(I,m,fl,fh,-4+xpow,-4)
    + 0.00267062 /n/n * x(I,m,fl,fh,-6+xpow,-6)
    - 0.0231902 /n/n * x(I,m,fl,fh,-5+xpow,-5)
    + 0.164152 /n/n * x(I,m,fl,fh,-4+xpow,-4)
    - 0.00317969 /n * x(I,m,fl,fh,-6+xpow,-6)
    - 0.0272673 /n * x(I,m,fl,fh,-4+xpow,-4)

    - 0.0654532/3. /n/n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,-5+xpow,-5)
    + 0.0654532/3. /n/n/n * lx(I,m,fl,fh,-5+xpow,-5)
    + 0.0519086/3. /n/n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,-4+xpow,-4)
    - 0.0519086/3. /n/n/n * lx(I,m,fl,fh,-4+xpow,-4)

    + 0.0237375 * x(I,m,fl,fh,-2+xpow,-2)
    - 2.75846 * x(I,m,fl,fh,-1+xpow,-1)
    - 0.0978549 /n/n/n * x(I,m,fl,fh,-3+xpow,-3)
    + 0.33489 /n/n/n * x(I,m,fl,fh,-2+xpow,-2)
    - 0.923641 /n/n/n * x(I,m,fl,fh,-1+xpow,-1)
    - 0.0583737 /n/n * x(I,m,fl,fh,-3+xpow,-3)
    + 0.79863 /n/n * x(I,m,fl,fh,-2+xpow,-2)
    - 18.2337 /n/n * x(I,m,fl,fh,-1+xpow,-1)
    + 0.251066 /n * x(I,m,fl,fh,-3+xpow,-3)
    - 0.0738498 /n * x(I,m,fl,fh,-2+xpow,-2)
    + 5.21599 /n * x(I,m,fl,fh,-1+xpow,-1)
    + 0.151994 * n * x(I,m,fl,fh,-2+xpow,-2)

    + 1.00463/3. /n/n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,-2+xpow,-2)
    - 1.00463/3. /n/n/n * lx(I,m,fl,fh,-2+xpow,-2)

    - 0.413954/3. /n/n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,-3+xpow,-3)
    + 0.413954/3. /n/n/n * lx(I,m,fl,fh,-3+xpow,-3)

    + 2.73357/3. /n/n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,-2+xpow,-2)
    - 2.73357/3. /n/n/n * lx(I,m,fl,fh,-2+xpow,-2)

    - 3.86244/3. /n/n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,-1+xpow,-1)
    + 3.86244/3. /n/n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,-1+xpow,-1)

    - 0.514796/3. /n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,-3+xpow,-3)
    + 0.514796/3. /n/n * lx(I,m,fl,fh,-3+xpow,-3)

    + 0.510333/3. /n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,-2+xpow,-2)
    - 0.510333/3. /n/n * lx(I,m,fl,fh,-2+xpow,-2)

    - 3.24404/3. /n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,-1+xpow,-1)
    + 3.24404/3. /n/n * lx(I,m,fl,fh,-1+xpow,-1)

    - 2.58077/3. /n * log(1./m/LAL_PI) * x(I,m,fl,fh,-1+xpow,-1)
    + 2.58077/3. /n * lx(I,m,fl,fh,-1+xpow,-1)

    + 27812.3/3. /n/n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,6+xpow,6)
    - 27812.3/3. /n/n/n * lx(I,m,fl,fh,6+xpow,6)

    - 83437./9. /n/n/n*log(1./m/LAL_PI)*log(1./m/LAL_PI)*x(I,m,fl,fh,6+xpow,6)
    + 83437./9. /n/n/n*log(1./m/LAL_PI)*lx(I,m,fl,fh,6+xpow,6)
    - 83437./9. /n/n/n*lsqx(I,m,fl,fh,6+xpow,6)

    + 46.103 * x(I,m,fl,fh,3+xpow,3)
    - 1059.73 * x(I,m,fl,fh,4+xpow,4)
    + 274.324 /n/n/n * x(I,m,fl,fh,3+xpow,3)
    - 612.301 /n/n/n * x(I,m,fl,fh,4+xpow,4)
    - 4615.59 /n/n/n * x(I,m,fl,fh,5+xpow,5)
    + 2066.38 /n/n * x(I,m,fl,fh,3+xpow,3)
    - 1263.58 /n/n * x(I,m,fl,fh,4+xpow,4)
    + 32.6094 /n * x(I,m,fl,fh,3+xpow,3)
    + 145.469 /n * x(I,m,fl,fh,4+xpow,4)
    - 1489.45 /n * x(I,m,fl,fh,5+xpow,5)
    - 201.857 * n * x(I,m,fl,fh,3+xpow,3)
    + 159.531 *n * x(I,m,fl,fh,4+xpow,4)
    + 87.2883 *n*n * x(I,m,fl,fh,3+xpow,3)

    + 1955.91/3. * log(1./m/LAL_PI) * x(I,m,fl,fh,4+xpow,4)
    - 1955.91/3. * lx(I,m,fl,fh,4+xpow,4)

    + 1087.21/3. /n/n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,3+xpow,3)
    - 1087.21/3. /n/n/n * lx(I,m,fl,fh,3+xpow,3)

    - 2758.97/3. /n/n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,4+xpow,4)
    + 2758.97/3. /n/n/n * lx(I,m,fl,fh,4+xpow,4)

    + 13846.8/3. /n/n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,5+xpow,5)
    - 13846.8/3. /n/n/n * lx(I,m,fl,fh,5+xpow,5)

    - 436.407/3. /n * log(1./m/LAL_PI) * x(I,m,fl,fh,4+xpow,4)
    + 436.407/3. /n * lx(I,m,fl,fh,4+xpow,4)

    - 4468.36/3. /n * log(1./m/LAL_PI) * x(I,m,fl,fh,5+xpow,5)
    + 4468.36/3. /n * lx(I,m,fl,fh,5+xpow,5)

    - 325.984/3. * log(1./m/LAL_PI) * x(I,m,fl,fh,4+xpow,4)
    + 325.984/3. * lx(I,m,fl,fh,4+xpow,4)

    + 1637.83/3. /n/n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,3+xpow,3)
    - 1637.83/3. /n/n/n * lx(I,m,fl,fh,3+xpow,3)

    - 3075.22/3. /n/n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,4+xpow,4)
    + 3075.22/3. /n/n/n * lx(I,m,fl,fh,4+xpow,4)

    + 9321.18/3. /n/n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,5+xpow,5)
    - 9321.18/3. /n/n/n * lx(I,m,fl,fh,5+xpow,5)

    - 448.655/3. /n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,3+xpow,3)
    + 448.655/3. /n/n * lx(I,m,fl,fh,3+xpow,3)

    - 10934.4/3. /n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,4+xpow,4)
    + 10934.4/3. /n/n * lx(I,m,fl,fh,4+xpow,4)

    + 7614.02/3. /n/n * log(1./m/LAL_PI) * x(I,m,fl,fh,5+xpow,5)
    - 7614.02/3. /n/n * lx(I,m,fl,fh,5+xpow,5)

    + 57.7191/3. /n * log(1./m/LAL_PI) * x(I,m,fl,fh,3+xpow,3)
    - 57.7191/3. /n * lx(I,m,fl,fh,3+xpow,3)

    + 145.469/3. /n * log(1./m/LAL_PI) * x(I,m,fl,fh,4+xpow,4)
    - 145.469/3. /n * lx(I,m,fl,fh,4+xpow,4)

    - 2978.9/3. /n * log(1./m/LAL_PI) * x(I,m,fl,fh,5+xpow,5)
    + 2978.9/3. /n * lx(I,m,fl,fh,5+xpow,5)

    - 3261.62/9. /n/n * log(1./m/LAL_PI)*log(1/m/LAL_PI)*x(I,m,fl,fh,3+xpow,3)
    + 3261.62/9. /n/n * log(1./m/LAL_PI)*lx(I,m,fl,fh,3+xpow,3)
    - 3261.62/9. /n/n/n * lsqx(I,m,fl,fh,3+xpow,3)

    + 3233.34/9. /n/n/n * log(1./m/LAL_PI)*log(1/m/LAL_PI)*x(I,m,fl,fh,4+xpow,4)
    - 3233.34/9. /n/n/n * log(1./m/LAL_PI)*lx(I,m,fl,fh,4+xpow,4)
    + 3233.34/9. /n/n/n * lsqx(I,m,fl,fh,4+xpow,4)

    + 1077.78/9. /n/n/n * log(1./m/LAL_PI)*log(1/m/LAL_PI)*x(I,m,fl,fh,4+xpow,4)
    - 1077.78/9. /n/n/n * log(1./m/LAL_PI)*lx(I,m,fl,fh,4+xpow,4)
    + 1077.78/9. /n/n/n * lsqx(I,m,fl,fh,4+xpow,4)

    );
  return output;
  /*return output /2./LAL_PI/2./LAL_PI/I->flow/I->flow; */
  }

static REAL8 JPsiTime(REAL8 mass1, REAL8 mass2,
         REAL8 fl, REAL8 fh, INT4 xpow, IMRBankCumulativeNoiseMoments *I)
  {
  /*printf("computing JPSITIME\n");*/
  REAL8 output = 0;
  REAL8 m = (mass1+mass2);
  /* FINISH THIS!!! */
  output = -2.0 * LAL_PI * x(I,m,fl,fh,3+xpow,0);
  return output;/*2./LAL_PI;*//*/I->flow/2./LAL_PI;*/
  /*return output /LAL_PI; */
  }


static REAL8 JPsiTimePsiTime(REAL8 mass1, REAL8 mass2,
         REAL8 fl, REAL8 fh, INT4 xpow, IMRBankCumulativeNoiseMoments *I)
  {
  /*printf("computing JPSITIMEPSITIME\n");*/
  REAL8 output = 0;
  REAL8 m = (mass1+mass2);
  /* FINISH THIS!!! */
  /*printf("setting output for %f %f %f %d\n",m,fl,fh,xpow);*/
  output = 4.0 * LAL_PI * LAL_PI * x(I,m,fl,fh,6+xpow,0);
  return output;/*4./LAL_PI/LAL_PI;*//*/I->flow/I->flow/2./LAL_PI/2./LAL_PI;*/
  /*return output /LAL_PI/LAL_PI;*/
  }


static REAL8 JPsiTimePsiM(REAL8 mass1, REAL8 mass2,
         REAL8 fl, REAL8 fh, INT4 xpow, IMRBankCumulativeNoiseMoments *I)
  {
  REAL8 output = 0;
  REAL8 n = eta(mass1,mass2);
  REAL8 m = (mass1+mass2);
  /* FINISH THIS!!! */
  output = 1.0 / n / m / m * (0.0
    + 0.230345 * x(I,m,fl,fh, 0+xpow,0)
    - 23.649 * x(I,m,fl,fh,3+xpow,3)
    - 604.976 * x(I,m,fl,fh,6+xpow,6)
    + 0.0364214 * x(I,m,fl,fh,-2+xpow,-2)
    - 2.30058 * x(I,m,fl,fh,1+xpow,1)
    + 1.00894 * x(I,m,fl,fh,2+xpow,2)
    + 43.4484 * x(I,m,fl,fh,4+xpow,4)
    - 200.797 * x(I,m,fl,fh,5+xpow,5)
    + 0.286458 * n * x(I,m,fl,fh,0+xpow,0)
    + 3.34127 * n * x(I,m,fl,fh,3+xpow,3)
    + 1.80515 * n * x(I,m,fl,fh,2+xpow,2)
    + 237.845 * n * x(I,m,fl,fh,4+xpow,4)
    - 165.621 * n * x(I,m,fl,fh,5+xpow,5)
    + 1.43607 * n*n * x(I,m,fl,fh,2+xpow,2)
    - 3.16425 * n*n * x(I,m,fl,fh,4+xpow,4)
    + 64.7973 * n*n * x(I,m,fl,fh,5+xpow,5)
    + 7.09083 * n*n*n * x(I,m,fl,fh,4+xpow,4)

    + 1814.93/3. * log(1./m/LAL_PI) * x(I,m,fl,fh,6+xpow,6)
    - 1814.93/3. * lx(I,m,fl,fh,6+xpow,6)
    - 23.444/3. * log(1./m/LAL_PI) * x(I,m,fl,fh,4+xpow,4)
    + 23.444/3. * lx(I,m,fl,fh,4+xpow,4)

  );
  return output;
  /*return output /LAL_PI/2./LAL_PI/I->flow;*/
  }

static REAL8 JPsiTimePsiEta(REAL8 mass1, REAL8 mass2,
         REAL8 fl, REAL8 fh, INT4 xpow, IMRBankCumulativeNoiseMoments *I)
  {
  REAL8 output = 0;
  REAL8 n = eta(mass1,mass2);
  REAL8 m = (mass1+mass2);
  /* FINISH THIS!!! */
  output = 1.0 / n / n / m * (0.0

    - 1814.93/3. * log(1./m/LAL_PI) * x(I,m,fl,fh,6+xpow,6)
    + 1814.93/3. * lx(I,m,fl,fh,6+xpow,6)

    + 0.230345/3. * x(I,m,fl,fh,0+xpow,0)
    + 23.649/3. *  x(I,m,fl,fh,3+xpow,3)
    + 0.0218528 *  x(I,m,fl,fh,-2+xpow,-2)
    - 3.45086 * x(I,m,fl,fh,-1+xpow,-1)
    + 3.20682 * x(I,m,fl,fh,2+xpow,2)
    - 60.0133 * x(I,m,fl,fh,4+xpow,4)
    + 301.196 * x(I,m,fl,fh,5+xpow,5)
    - 4.30821 * n*n * x(I,m,fl,fh,2+xpow,2)
    - 9.49274 * n*n * x(I,m,fl,fh,4+xpow,4)
    + 97.196 * n*n * x(I,m,fl,fh,5+xpow,5)
    + 42.545 * n*n*n * x(I,m,fl,fh,4+xpow,4)

    - 70.947/3. * log(1./m/LAL_PI) * x(I,m,fl,fh,3+xpow,3)
    + 70.947/3. * lx(I,m,fl,fh,3+xpow,3)
    + 70.3319/3. * log(1./m/LAL_PI) * x(I,m,fl,fh,4+xpow,4)
    - 70.3319/3. * lx(I,m,fl,fh,4+xpow,4)

  );
  return output;
  /*return output /LAL_PI/2./LAL_PI/I->flow;*/
  }


static REAL8 XLALComputeIMRBankMetricTimeTime(REAL8 mass1, REAL8 mass2,
         REAL8 fl, REAL8 fh, INT4 xpow, IMRBankCumulativeNoiseMoments *I)
  {
  return (JPsiTimePsiTime(mass1, mass2, fl, fh, xpow, I)
    - JPsiTime(mass1, mass2, fl, fh, xpow, I)
      * JPsiTime(mass1, mass2, fl, fh, xpow, I))/2.0;

  }

static REAL8 XLALComputeIMRBankMetricTimeMass(REAL8 mass1, REAL8 mass2,
         REAL8 fl, REAL8 fh, INT4 xpow, IMRBankCumulativeNoiseMoments *I)
  {

  return (JPsiTimePsiM(mass1, mass2, fl, fh, xpow, I)
    - JPsiTime(mass1, mass2, fl, fh, xpow, I)
      * JPsiM(mass1, mass2, fl, fh, xpow, I))/2.0;

  }

static  REAL8 XLALComputeIMRBankMetricTimeEta(REAL8 mass1, REAL8 mass2,
         REAL8 fl, REAL8 fh, INT4 xpow, IMRBankCumulativeNoiseMoments *I)
  {

  return (JPsiTimePsiEta(mass1, mass2, fl, fh, xpow, I)
    - JPsiTime(mass1, mass2, fl, fh, xpow, I)
      * JPsiEta(mass1, mass2, fl, fh, xpow, I))/2.0;
  }

static REAL8 XLALComputeIMRBankMetricMassMass(REAL8 mass1, REAL8 mass2,
         REAL8 fl, REAL8 fh, INT4 xpow, IMRBankCumulativeNoiseMoments *I)
  {

  return (JPsiMPsiM(mass1, mass2, fl, fh, xpow, I)
    - JPsiM(mass1, mass2, fl, fh, xpow, I)
      * JPsiM(mass1, mass2, fl, fh, xpow, I))/2.0;
  }


static REAL8 XLALComputeIMRBankMetricEtaEta(REAL8 mass1, REAL8 mass2,
         REAL8 fl, REAL8 fh, INT4 xpow, IMRBankCumulativeNoiseMoments *I)
  {
  return (JPsiEtaPsiEta(mass1, mass2, fl, fh, xpow, I)
    - JPsiEta(mass1, mass2, fl, fh, xpow, I)
    * JPsiEta(mass1, mass2, fl, fh, xpow, I))/2.0;

  }

static REAL8 XLALComputeIMRBankMetricMassEta(REAL8 mass1, REAL8 mass2,
         REAL8 fl, REAL8 fh, INT4 xpow, IMRBankCumulativeNoiseMoments *I)
  {

  return (JPsiMPsiEta(mass1, mass2, fl, fh, xpow, I)
    - JPsiM(mass1, mass2, fl, fh, xpow, I)
    * JPsiEta(mass1, mass2, fl, fh, xpow, I))/2.0;
  }

#if 0
static int printmetric(IMRBankMetric *metric,FILE *FP)
  {
  REAL8 MM = metric->data[1][1]-metric->data[0][1]*metric->data[0][1]/metric->data[0][0];
  REAL8 MN = metric->data[1][2]-metric->data[0][1]*metric->data[0][2]/metric->data[0][0];
  REAL8 NN = metric->data[2][2]-metric->data[0][2]*metric->data[0][2]/metric->data[0][0];
  REAL8 vol = sqrt(fabs( (MM*NN-MN*MN) ));
  fprintf(FP,"%e %e %e %e %e %e %e %e %e\n",metric->m1, metric->m2,metric->data[0][0],metric->data[0][1],metric->data[0][2],metric->data[1][1],metric->data[1][2],metric->data[2][2],vol);
  return 0;
  }
#endif


static double tau0fromm1m2(REAL8 mass1, REAL8 mass2, REAL8 flow)
  {
  REAL8 M = mass1+mass2;
  REAL8 n = eta(mass1,mass2);
  REAL8 c0 = 5. / 256. / pow(LAL_PI*flow,8./3.);
  return c0 / n / pow(M,5./3.);
  }

static double tau3fromm1m2(REAL8 mass1, REAL8 mass2, REAL8 flow)
  {
  REAL8 M = mass1+mass2;
  REAL8 n = eta(mass1,mass2);
  REAL8 c3 = LAL_PI / 8. / pow(LAL_PI*flow,5./3.);
  return c3 / n / pow(M,2./3.);
  }


static int IMRBankMetricToTau0Tau3(IMRBankMetric *metric,IMRBankCumulativeNoiseMoments *moments)
  {
  REAL8 dtdt, dtdt0, dtdt3, dmdt, dmdt0, dmdt3, dndt, dndt0, dndt3;
  REAL8 h00,h01,h02,h10,h20,h11,h12,h21,h22;
  /* fix this */

  REAL8 c8_3 = 2. * LAL_PI / pow(moments->flow,8./3.);
  REAL8 c5_3 = 2. * LAL_PI / pow(moments->flow,5./3.);
  REAL8 p2 = 2. * LAL_PI;

  REAL8 g00 = metric->data[0][0];/* /2./LAL_PI/flow/2./LAL_PI/flow;*/
  REAL8 g01 = metric->data[0][1];/* /2./LAL_PI/flow/2./LAL_PI/flow;*/
  REAL8 g02 = metric->data[0][2];/* /2./LAL_PI/flow/2./LAL_PI/flow;*/
  REAL8 g11 = metric->data[1][1];/* /2./LAL_PI/flow/2./LAL_PI/flow;*/
  REAL8 g12 = metric->data[1][2];/* /2./LAL_PI/flow/2./LAL_PI/flow;*/
  REAL8 g22 = metric->data[2][2];/* /2./LAL_PI/flow/2./LAL_PI/flow;*/

  metric->tau0 = tau0fromm1m2(metric->m1,metric->m2,moments->flow);
  metric->tau3 = tau3fromm1m2(metric->m1,metric->m2,moments->flow);

  dtdt = 1;
  dtdt0 = 0;
  dtdt3 = 0;
  dmdt = 0;
  dmdt0 = (0.0 - (metric->M)/metric->tau0);
  dmdt3 = ((metric->M)/metric->tau3);
  dndt = 0;
  dndt0 = (2./3.*metric->eta / metric->tau0);
  dndt3 = (0.0-5./3.*metric->eta / metric->tau3);


  h00 = dtdt*g00 + dmdt*g01 + dndt*g02;
  h01 = dtdt0*g00 + dmdt0*g01 + dndt0*g02;
  h02 = dtdt3*g00 + dmdt3*g01 + dndt3*g02;
  h10 = g01*dtdt + g11*dmdt + g12*dndt;
  h11 = dtdt0*g01 + dmdt0*g11 + dndt0*g12;
  h12 = dtdt3*g01 + dmdt3*g11 + dndt3*g12;
  h20 = dtdt*g02 + dmdt*g12 + dndt*g22;
  h21 = dtdt0*g02 + dmdt0*g12 + dndt0*g22;
  h22 = dtdt3*g02 + dmdt3*g12 + dndt3*g22;
  /* fix here down*/
  metric->data[0][0] = p2 * p2 * (h00*dtdt + h10*dmdt + h20*dndt);
  metric->data[0][1] = p2 * c8_3 * (h01*dtdt + h11*dmdt + h21*dndt);
  metric->data[0][2] =  p2 * c5_3 * (h02*dtdt + h12*dmdt + h22*dndt);
  metric->data[1][0] =  p2 * c5_3 * (h00*dtdt0 + h10*dmdt0 + h20*dndt0);
  metric->data[1][1] =  c8_3 * c8_3 * (h01*dtdt0 + h11*dmdt0 + h21*dndt0);
  metric->data[1][2] =  c5_3 * c8_3 * (h02*dtdt0 + h12*dmdt0 + h22*dndt0);
  metric->data[2][0] =  p2 * c5_3 * (h00*dtdt3 + h10*dmdt3 + h20*dndt3);
  metric->data[2][1] =  c8_3 * c5_3 * (h01*dtdt3 + h11*dmdt3 + h21*dndt3);
  metric->data[2][2] =  c5_3 * c5_3 * (h02*dtdt3 + h12*dmdt3 + h22*dndt3);

  /*printf("%e %e %e\n %e %e %e\n %e %e %e\n\n",metric->data[0][0],metric->data[0][1],metric->data[0][2],metric->data[1][0],metric->data[1][1],metric->data[1][2],metric->data[2][0],metric->data[2][1],metric->data[2][2]);*/
  return 0;
  }

static int XLALComputeIMRBankMetric(REAL8 mass1, REAL8 mass2, IMRBankCumulativeNoiseMoments *moments, IMRBankMetric *metric)
  {
  REAL8 q = 0;
  REAL8 m = (mass1+mass2);
  REAL8 fm = 0;
  REAL8 fr = 0;
  REAL8 fl = moments->flow;

  metric->m1 = mass1;
  metric->m2 = mass2;
  metric->M = mass1+mass2;
  metric->eta = eta(mass1,mass2);

  if (mass1 < mass2) q = mass1/mass2;
  else q = mass2/mass1;

  fm = (0.8 * q * q * q - 2.6 * q * q + 2.8 * q + 1.0) / LAL_PI / m / pow(6,3.0/2.);
  /*fm = 1.0/LAL_PI/m/pow(2.8,3./2.);*/
  if (fm < fl) fl = fm-10;
  /*fr = 1+1.0/LAL_PI/m/pow(3.0,3./2.);*/
  fr = (0.099*q*q*q -0.27*q*q + 0.25*q +0.19) / LAL_PI / m;
  if (fr<fm) fl = fr-10;

  /*printf("starting metric 00\n");*/
  metric->data[0][0] =
    ( XLALComputeIMRBankMetricTimeTime(mass1,mass2,fl,fm,-7,moments) +
    XLALComputeIMRBankMetricTimeTime(mass1,mass2,fm,fr,-4,moments) ) /
    ( x(moments,1,fl,fm,-7,0) + x(moments,1,fm,fr,-4,0) );
  /*printf("starting metric 01\n");*/

  metric->data[0][1] = metric->data[1][0] =
    ( XLALComputeIMRBankMetricTimeMass(mass1,mass2,fl,fm,-7,moments) +
    XLALComputeIMRBankMetricTimeMass(mass1,mass2,fm,fr,-4,moments) ) /
    ( x(moments,1,fl,fm,-7,0) + x(moments,1,fm,fr,-4,0) );

  /*printf("starting metric 02\n");*/
  metric->data[0][2] = metric->data[2][0] =
    ( XLALComputeIMRBankMetricTimeEta(mass1,mass2,fl,fm,-7,moments) +
    XLALComputeIMRBankMetricTimeEta(mass1,mass2,fm,fr,-4,moments) ) /
    ( x(moments,1,fl,fm,-7,0) + x(moments,1,fm,fr,-4,0) );

  /*printf("starting metric 11\n");*/
  metric->data[1][1] =
    ( XLALComputeIMRBankMetricMassMass(mass1,mass2,fl,fm,-7,moments) +
    XLALComputeIMRBankMetricMassMass(mass1,mass2,fm,fr,-4,moments) ) /
    ( x(moments,1,fl,fm,-7,0) + x(moments,1,fm,fr,-4,0) );

  /*printf("starting metric 22\n");*/
  metric->data[2][2] =
    ( XLALComputeIMRBankMetricEtaEta(mass1,mass2,fl,fm,-7,moments) +
    XLALComputeIMRBankMetricEtaEta(mass1,mass2,fm,fr,-4,moments) ) /
    ( x(moments,1,fl,fm,-7,0) + x(moments,1,fm,fr,-4,0) );

  /*printf("starting metric 12\n");*/
  metric->data[1][2] = metric->data[2][1] =
    ( XLALComputeIMRBankMetricMassEta(mass1,mass2,fl,fm,-7,moments) +
    XLALComputeIMRBankMetricMassEta(mass1,mass2,fm,fr,-4,moments) ) /
    ( x(moments,1,fl,fm,-7,0) + x(moments,1,fm,fr,-4,0) );
  /*printMetric(metric);*/
  return 0;

  }

static REAL8 jacobian(REAL8 m1, REAL8 m2)
  {
  return fabs(m2-m1)/(m2+m1)/(m2+m1);
  }

static REAL8 mDensity(REAL8 m1, REAL8 m2, IMRBankCumulativeNoiseMoments *I)
  {
  IMRBankMetric metric;
  REAL8 MM,MN,NN;
  XLALComputeIMRBankMetric(m1,m2,I,&metric);

  MM = metric.data[1][1]-metric.data[0][1]*metric.data[0][1]/metric.data[0][0];
  MN = metric.data[1][2]-metric.data[0][1]*metric.data[0][2]/metric.data[0][0];
  NN = metric.data[2][2]-metric.data[0][2]*metric.data[0][2]/metric.data[0][0];

  return sqrt(fabs( (MM*NN-MN*MN) ));
  }

#if 0
static REAL8 Mfromtau0tau3(REAL8 t0, REAL8 t3, REAL8 flow)
  {
  REAL8 c0 = 5. / 256. / pow(LAL_PI*flow,8./3.);
  REAL8 c3 = LAL_PI / 8. / pow(LAL_PI*flow,5./3.);
  return c0/c3*t3/t0;
  }

static REAL8 Etafromtau0tau3(REAL8 t0, REAL8 t3, REAL8 flow)
  {
  REAL8 c0 = 5. / 256. / pow(LAL_PI*flow,8./3.);
  REAL8 c3 = LAL_PI / 8. / pow(LAL_PI*flow,5./3.);
  return pow(c3/t3,5./3.)*pow(t0/c0,2./3.);
  }

static REAL8 m1fromMassEta(REAL8 M, REAL8 n)
  {
  return (M + sqrt(M*M-4.0*n*M*M))/2.0;
  }

static REAL8 m2fromMassEta(REAL8 M, REAL8 n)
  {
  return (M - sqrt(M*M-4.0*n*M*M))/2.0;
  }
#endif

#if 0
static REAL8 m1fromtau0tau3(REAL8 t0, REAL8 t3, REAL8 flow)
  {
  REAL8 M = Mfromtau0tau3(t0,t3,flow);
  REAL8 n = Etafromtau0tau3(t0,t3,flow);
  /*printf("M %f eta %f\n",M,n);*/
  if (n > 0.25) return 0.0;
  return m1fromMassEta(M,n);
  }

static REAL8 m2fromtau0tau3(REAL8 t0, REAL8 t3, REAL8 flow)
  {
  REAL8 M = Mfromtau0tau3(t0,t3,flow);
  REAL8 n = Etafromtau0tau3(t0,t3,flow);
  /*printf("M %f eta %f\n",M,n);*/
  if (n > 0.25) return 0.0;
  return m2fromMassEta(M,n);
  }
#endif

#if 0
static REAL8 tDensity(REAL8 t0, REAL8 t3, IMRBankCumulativeNoiseMoments *I)
  {
  IMRBankMetric metric;
  REAL8 MM,MN,NN;
  REAL8 m1 = m1fromtau0tau3(t0,t3,I->flow);
  REAL8 m2 = m2fromtau0tau3(t0,t3,I->flow);
  /*printf("m1 %f m2 %f \n", m1,m2);*/
  XLALComputeIMRBankMetric(m1,m2,I,&metric);
  IMRBankMetricToTau0Tau3(&metric,I);

  MM = metric.data[1][1]-metric.data[0][1]*metric.data[0][1]/metric.data[0][0];
  MN = metric.data[1][2]-metric.data[0][1]*metric.data[0][2]/metric.data[0][0];
  NN = metric.data[2][2]-metric.data[0][2]*metric.data[0][2]/metric.data[0][0];
  return sqrt(fabs( (MM*NN-MN*MN) ));
  }
#endif

static REAL8 eta_chirpmass_volume(REAL8 mbox[3])
  {
  REAL8 m1 = mbox[0];
  REAL8 m2 = mbox[1];
  REAL8 size = mbox[2];
  REAL8 e1 = eta(m1,m2);
  REAL8 e2 = eta(m1+size,m2);
  REAL8 e3 = eta(m1,m2+size);
  REAL8 e4 = eta(m1+size,m2+size);
  REAL8 c1 = chirpmass(m1,m2);
  REAL8 c2 = chirpmass(m1+size,m2);
  REAL8 c3 = chirpmass(m1,m2+size);
  REAL8 c4 = chirpmass(m1+size,m2+size);
  REAL8 maxe = e1;
  REAL8 maxc = c1;
  REAL8 mine = e1;
  REAL8 minc = c1;

  if (e2 > maxe) maxe = e2;
  if (e3 > maxe) maxe = e3;
  if (e4 > maxe) maxe = e4;

  if (c2 > maxc) maxc = c2;
  if (c3 > maxc) maxc = c3;
  if (c4 > maxc) maxc = c4;

  if (e2 < maxe) mine = e2;
  if (e3 < maxe) mine = e3;
  if (e4 < maxe) mine = e4;

  if (c2 < minc) minc = c2;
  if (c3 < minc) minc = c3;
  if (c4 < minc) minc = c4;

  return (maxe-mine) * (maxc-minc);

  }

static REAL8 integrateMassVolume(REAL8 mbox[3],
                       IMRBankCumulativeNoiseMoments *I)

  {
  REAL8 m1 = mbox[0];
  REAL8 m2 = mbox[1];
  REAL8 size = mbox[2];
  REAL8 volume = 0;
  REAL8 sf = 0.80; /* gives volume that is 2x overlapping */
  REAL8 g1 = mDensity(sf*m1,m2*sf,I);
  REAL8 g2 = mDensity(sf*m1+size/sf/sf,m2*sf,I);
  REAL8 g3 = mDensity(m1*sf,m2*sf+size/sf/sf,I);
  REAL8 g4 = mDensity(m1*sf+size/sf/sf,m2*sf+size/sf/sf,I);

  REAL8 j1 = jacobian(m1*sf,m2*sf);
  REAL8 j2 = jacobian(m1*sf+size/sf/sf,m2*sf);
  REAL8 j3 = jacobian(m1*sf,m2*sf+size/sf/sf);
  REAL8 j4 = jacobian(m1*sf+size/sf/sf,m2*sf+size/sf/sf);
  REAL8 maxg = g1;
  REAL8 maxj = j1;

  if (j2 > maxj) maxj = j2;
  if (j3 > maxj) maxj = j3;
  if (j4 > maxj) maxj = j4;

  if (g2 > maxg) maxg = g2;
  if (g3 > maxg) maxg = g3;
  if (g4 > maxg) maxg = g4;
  volume = maxg * maxj * size/sf * size/sf;
  /*fprintf(stderr,"volume %e\n",volume);*/
  return volume;
  }


static REAL8 XLALComputeNumberOfIMRTemplatesInSquareIMRBankMassRegion(
                      REAL8 mbox[3],
                      REAL8 mm, IMRBankCumulativeNoiseMoments *I)
  {
  REAL8 out;
  REAL8 vol;

  vol = integrateMassVolume(mbox,I);
  out = vol / mm / (4.0 * LAL_PI * LAL_PI * I->flow * I->flow)
      / (4.0 * LAL_PI * LAL_PI * I->flow * I->flow);
  return out;
  }


static IMRBankMassRegion * createIMRBankMassRegion(REAL8 mass1, REAL8 mass2, REAL8 size)
  {
  IMRBankMassRegion *region = calloc(1, sizeof(IMRBankMassRegion));
  region->mbox[0] = mass1;
  region->mbox[1] = mass2;
  region->mbox[2] = size;
  return region;
  }


static int appendtotailMass(IMRBankMassRegion *elem, IMRBankMassRegion **tail)
  {
  (*tail)->next = elem;
  *tail = elem;
  return 0;
  }

static int divideAndConquerMass(IMRBankMassRegion *list,
                            IMRBankMassRegion **tail)

  {
  REAL8 m1 = list->mbox[0];
  REAL8 m2 = list->mbox[1];
  REAL8 size = list->mbox[2];
  REAL8 newsize = size/2.;
  /* No checking if it is null  BAD */
  appendtotailMass( createIMRBankMassRegion(m1,m2,newsize), tail);
  appendtotailMass( createIMRBankMassRegion(m1+newsize,m2,newsize), tail);
  appendtotailMass( createIMRBankMassRegion(m1,m2+newsize,newsize), tail);
  appendtotailMass( createIMRBankMassRegion(m1+newsize,m2+newsize,newsize), tail);
  return 0;
  }


static int addtemplatesMass(REAL8 mbox[3],
                        InspiralCoarseBankIn *in, SnglInspiralTable **head,
			IMRBankCumulativeNoiseMoments *I
			)
  {
  /*Now this function always assumes one template to be added */
  REAL8 m1 = 0;
  REAL8 m2 = 0;
  REAL8 mtot = mbox[0]+mbox[1];
  REAL8 size = mbox[2];
  IMRBankMetric metric;
  REAL8 MM,MN,NN;
  /* check to see if it is the bottom half of the mass/mass plane */
  /* and that is falls within the boundaries */

  m1 = mbox[0]+size/2.;/* * gsl_rng_uniform_pos(r); */
  m2 = mbox[1]+size/2.; /* * gsl_rng_uniform_pos(r);*/
  /* Don't forget the input masses are in Kg */
  if ( (mtot <= 0.99 * in->MMin*LAL_MTSUN_SI)
        || (mtot >= 1.01 * in->MMax*LAL_MTSUN_SI)
        || (m2>m1)
	|| (m1 < 0.99 * in->mMin*LAL_MTSUN_SI)
	|| (m2 < 0.99 * in->mMin*LAL_MTSUN_SI)
	|| (m1 > 1.01 * in->mMax*LAL_MTSUN_SI)
	|| (m2 > 1.01 * in->mMax*LAL_MTSUN_SI)      )  return 0;

  fprintf(stderr, "template %d\n",tmpltcnt++);
  XLALComputeIMRBankMetric(m1-size/2.,m2-size/2.,I,&metric);
  IMRBankMetricToTau0Tau3(&metric,I);
  /* project out the time dimension */
  MM = metric.data[1][1]-metric.data[0][1]*metric.data[0][1]/metric.data[0][0];
  MN = metric.data[1][2]-metric.data[0][1]*metric.data[0][2]/metric.data[0][0];
  NN = metric.data[2][2]-metric.data[0][2]*metric.data[0][2]/metric.data[0][0];
  (*head)->mass1 = m1 / LAL_MTSUN_SI;
  (*head)->mass2 = m2 / LAL_MTSUN_SI;
  (*head)->tau0 = metric.tau0;
  (*head)->tau3 = metric.tau3;
  /*printmetric(&metric,FP);*/
  (*head)->eta = eta(m1,m2);
  (*head)->mchirp = pow(m1*m2,0.6)/pow(m1+m2,0.2) / LAL_MTSUN_SI;
  (*head)->Gamma[0] = metric.data[0][0];
  (*head)->Gamma[1] = metric.data[0][1];
  (*head)->Gamma[2] = metric.data[0][2];
  (*head)->Gamma[3] = metric.data[1][1];
  (*head)->Gamma[4] = metric.data[1][2];
  (*head)->Gamma[5] = metric.data[2][2];
  (*head)->next = (SnglInspiralTable *) LALCalloc(1,sizeof(SnglInspiralTable));
  (*head) = (*head)->next;
  return 0;
  }



static int checkNumberOfTemplatesMass(IMRBankMassRegion *list,
                                  IMRBankMassRegion **tail,
				  InspiralCoarseBankIn *in,
				  IMRBankCumulativeNoiseMoments *I,
				  SnglInspiralTable **head
				  )
  {
  REAL8 mm = 1.0-in->mmCoarse;
  REAL8 numTmps = XLALComputeNumberOfIMRTemplatesInSquareIMRBankMassRegion(
                  list->mbox, mm, I);

  if (numTmps >= 1.0)  divideAndConquerMass(list,tail);
  else addtemplatesMass(list->mbox,in,head,I);
  return 0;
  }


static int destroyregionlistMass(IMRBankMassRegion *head)
  {
  IMRBankMassRegion *tmp = NULL;

  while(head)
    {
    tmp = head;
    head = head->next;
    free(tmp);
    }
  return 0;
  }


static int normalize_psd(InspiralCoarseBankIn *in)
  {
  UINT4 i;
  double f = 0;
  double min = 1.0;
  UINT4 startIX;
  REAL8Vector *vec = in->shf.data;
  startIX = floor(in->fLower / in->shf.deltaF);
  vec->data[0] = vec->data[1];
  for (i=0; i < vec->length; i++)
    if (vec->data[i] < min && vec->data[i] != 0) min = vec->data[i];
  for (i=0; i<vec->length; i++)
    vec->data[i] = 2.0 / min; /*10e50;*/
  return 0;
  }


/* This is the main function that actually is called by the bank code */
int XLALTileIMRBankMassRegion(InspiralCoarseBankIn *in, SnglInspiralTable **first)
  {
  /* Convert all masses to geometrized units */
  REAL8 mass1 = LAL_MTSUN_SI*(0.80*in->mMin);
  REAL8 mass2 = LAL_MTSUN_SI*(0.80*in->mMin);
  REAL8 size = LAL_MTSUN_SI*(1.25*in->mMax - 0.80*in->mMin);

  REAL8 flow = in->fLower;
  SnglInspiralTable *tab = NULL;
  SnglInspiralTable *tmp = NULL;
  IMRBankCumulativeNoiseMoments moments;
  IMRBankMassRegion *head = createIMRBankMassRegion(mass1,mass2,size);
  IMRBankMassRegion *list = head;
  IMRBankMassRegion *tail = head;
  int cnt;

  *first = (SnglInspiralTable *) LALCalloc(1,sizeof(SnglInspiralTable));
  tab = *first;

  normalize_psd(in);

  XLALCreateIMRBankCumulativeNoiseMoments(&moments,&(in->shf),flow);
  while(list)
    {
    checkNumberOfTemplatesMass(list,&tail,in,&moments,first);
    list = list->next;
    }
  XLALDestroyIMRBankCumulativeNoiseMoments(&moments);
  destroyregionlistMass(head);
  /* send back the beginning of the list */
  *first = tab;
  cnt = 0;
  while (tab)
    {
    if ( (tab->next) && !(tab->next->next))
      {
      tmp = tab->next;
      tab->next = NULL;
      free(tmp);
      break;
      }
    cnt++;
    tab = tab->next;
    }
  return cnt;
}
