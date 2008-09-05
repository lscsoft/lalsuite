#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/LALInspiralBank.h>
#include <lal/LIGOMetadataTables.h>
#include <gsl/gsl_rng.h>
#include <lal/IMRBank.h>


int XLALCreateLALIMRBankCumulativeNoiseMoments(LALIMRBankCumulativeNoiseMoments *moments, REAL8FrequencySeries *psd, REAL8 flow)
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
  /* Allocate memory for the moments */
  for (i=0; i < 15; i++)
    {
    moments->minus3[i] = NULL;
    moments->plus3[i] = NULL;
    moments->logminus3[i] = NULL;
    moments->logplus3[i] = NULL;
    }
  return 0;
  }

int XLALComputeIMRBankCumulativeNoiseMoment(LALIMRBankCumulativeNoiseMoments *moments, INT4 power, UINT4 logFlag)
  {
  REAL8FrequencySeries *psd = moments->psd;
  UINT4 length = psd->data->length;
  REAL8 deltaF = psd->deltaF;
  UINT4 j = 0;
  if (power < 0) 
    {
    if (!logFlag)
      {
      moments->minus3[-power] = XLALCreateREAL8Vector(length);
      moments->minus3[-power]->data[0] = 0;
      }
    else
      {
      moments->logminus3[-power] = XLALCreateREAL8Vector(length);
      moments->logminus3[-power]->data[0] = 0;
      } 
    }
  else
    {
    if (!logFlag)
      {
      moments->plus3[power] = XLALCreateREAL8Vector(length);
      moments->plus3[power]->data[0] = 0;
      }
    else
      {
      moments->logplus3[power] = XLALCreateREAL8Vector(length);
      moments->logplus3[power]->data[0] = 0;
      }
    }
  for (j = 1; j < length; j++)
    {
    if (psd->data->data[j] > 0.0)
      {
      if (power < 0)
        if (!logFlag)
	  {
          moments->minus3[-power]->data[j] = moments->minus3[-power]->data[j-1]
             + 1./pow(deltaF*j,-power / 3.0) / (psd->data->data[j]) * deltaF;
	  }
	
        else
          moments->logminus3[-power]->data[j] = 
                           moments->logminus3[-power]->data[j-1]
                           + 1./pow(deltaF*j,-power/3.0) 
			   / (psd->data->data[j])
                           * log(deltaF*j) * deltaF;
      else 
        if (!logFlag)
          moments->plus3[power]->data[j] = moments->plus3[power]->data[j-1]
                           + pow(deltaF*j,power/3.0) / (psd->data->data[j])
                           * deltaF;
        else
          moments->logplus3[power]->data[j] = 
                           moments->logplus3[power]->data[j-1]
                           + pow(deltaF*j,power/ 3.0) / (psd->data->data[j])
                           * log(deltaF*j) * deltaF;
      }
    else 
      {
      if (power < 0) 
        if (!logFlag)
          moments->minus3[-power]->data[j]=moments->minus3[-power]->data[j-1];
        else 
          moments->logminus3[-power]->data[j]=moments->logminus3[-power]->data[j-1];
      else 
        if (!logFlag)
          moments->plus3[power]->data[j]=moments->plus3[power]->data[j-1];
        else
          moments->logplus3[power]->data[j]=moments->logplus3[power]->data[j-1];
      } 
    }
  return 0;
  }




int XLALDestroyLALIMRBankCumulativeNoiseMoments(LALIMRBankCumulativeNoiseMoments *moments)
  {
  UINT4 i = 0;
  for (i = 1; i < 15; i++)
    {
    if (moments->minus3[i]) XLALDestroyREAL8Vector(moments->minus3[i]);
    if (moments->plus3[i]) XLALDestroyREAL8Vector(moments->plus3[i]);
    if (moments->logminus3[i]) XLALDestroyREAL8Vector(moments->logminus3[i]);
    if (moments->logplus3[i]) XLALDestroyREAL8Vector(moments->logplus3[i]);
    }
  return 0;
  }

static REAL8 x(LALIMRBankCumulativeNoiseMoments *moments,
               REAL8 m, REAL8 flow, REAL8 fhigh, INT4 power, INT4 mpower)
  {
  /*printf("getting moment %d\n",power);*/
  REAL8FrequencySeries *psd = moments->psd;
  UINT4 fl = floor(flow / moments->deltaF);
  UINT4 fh = floor(fhigh / moments->deltaF);
  REAL8 output = 0;
  if (fl > psd->data->length) fl = psd->data->length -2;
  if (fh > psd->data->length) fh = psd->data->length -2;
  if (power < 0) 
    {
    /*printf("checking if the neg moment has already been computed...\n");*/
    if (!moments->minus3[-power])
      {
      /*printf("it hasn't, so I'll compute it...\n");*/
      XLALComputeIMRBankCumulativeNoiseMoment(moments, power,0);
      }
    output = 1./pow(m,-mpower/3.0) 
        * (moments->minus3[-power]->data[fh]-moments->minus3[-power]->data[fl]);
    }
  else
    {
    /*printf("checking if the pos moment has already been computed...\n");*/
    if (!moments->plus3[power])
      {
      /*printf("it hasn't, so I'll compute it...\n");*/
      XLALComputeIMRBankCumulativeNoiseMoment(moments, power,0);
      }
    output = pow(m,mpower/3.0)
           * (moments->plus3[power]->data[fh]-moments->plus3[power]->data[fl]);
    }
  /*printf("output of x function %e\n",output);*/
  return output;
  }

static REAL8 lx(LALIMRBankCumulativeNoiseMoments *moments,
                REAL8 m, REAL8 flow, REAL8 fhigh, INT4 power, INT4 mpower)
  {
  REAL8FrequencySeries *psd = moments->psd;
  UINT4 fl = floor(flow / moments->deltaF);
  UINT4 fh = floor(fhigh / moments->deltaF);
  if (fl > psd->data->length) fl = psd->data->length -1;
  if (fh > psd->data->length) fh = psd->data->length -1;
  if (power < 0)
    {
    if (!moments->logminus3[-power])
      XLALComputeIMRBankCumulativeNoiseMoment(moments,power,1);
    return 
    /* FIXME check on this!!!!! Do I need log m too? */
    pow(m,mpower/3.0) * (moments->logminus3[-power]->data[fh] -  moments->logminus3[-power]->data[fl]);
    }
  else
    {
    if (!moments->logplus3[power])
      XLALComputeIMRBankCumulativeNoiseMoment(moments,power,1);
    return 
    /* FIXME check on this!!!!! Do I need log m too?*/
    pow(m,mpower/3.0) * (moments->logplus3[power]->data[fh] -  moments->logplus3[power]->data[fl] );
    }
  }
       
REAL8 XLALComputeLALIMRBankMetricMassMass(REAL8 mass1, REAL8 mass2, 
         REAL8 fl, REAL8 fh, INT4 xpow, LALIMRBankCumulativeNoiseMoments *I)
  {
  REAL8 output = 0;
  REAL8 n = mass1*mass2/(mass1+mass2)/(mass1+mass2);
  REAL8 m = (mass1+mass2)*LAL_MTSUN_SI;
  /*printf("n = %e, m = %e \n",n,m);*/
  output = -1.0 / 2.0 / m / m * (0.0 
    - 96.2849 / n * x(I,m,fl,fh,3+xpow,3)  
    - 0.0154577 / n * x(I,m,fl,fh,-5+xpow,-5)  
    - 0.073321  / n * x(I,m,fl,fh,-3+xpow,-3) 
    + 0.610247 / n * x(I,m,fl,fh,-2+xpow,-2) 
    - 0.214104 / n * x(I,m,fl,fh,-1+xpow,-1) 
    - 0.0911825 * x(I,m,fl,fh,-3+xpow,-3) 
    - 0.383064 * x(I,m,fl,fh,-1+xpow,-1) 
    - 0.304744 * n * x(I,m,fl,fh,-1+xpow,-1) 
    + 3.76385 / n * x(I,m,fl,fh,xpow,0) 
    - 3.36628 / n * x(I,m,fl,fh,1+xpow,1) 
    + 10.6526 / n * x(I,m,fl,fh,2+xpow,2)
    - 0.53178 * x(I,m,fl,fh,0+xpow,0)
    - 25.2361 * x(I,m,fl,fh,1+xpow,1)
    + 8.78645 * x(I,m,fl,fh,2+xpow,2)
    + 0.335737 * n *x(I,m,fl,fh,1+xpow,1)
    - 3.4376 * n * x(I,m,fl,fh,2+xpow,2)
    - 0.752361 * n * n * x(I,m,fl,fh,1+xpow,1) 
    + 2.48748 / n * (log(1)-1./3.*log(m*LAL_PI)) * x(I,m,fl,fh,1+xpow,1)
    - 2.48748 / 3. / n * lx(I,m,fl,fh,1+xpow,1) 
    );
    
  /*printf("Unnormalized metric component %e\n",output);*/

  /* WARNING THIS RESCALES TO SOLAR MASS ! */
  return output*LAL_MTSUN_SI*LAL_MTSUN_SI;
 
  }

REAL8 XLALComputeLALIMRBankMetricEtaEta(REAL8 mass1, REAL8 mass2,
         REAL8 fl, REAL8 fh, INT4 xpow, LALIMRBankCumulativeNoiseMoments *I)
  {
  REAL8 output = 0;
  REAL8 n = mass1*mass2/(mass1+mass2)/(mass1+mass2);
  REAL8 m = (mass1+mass2)*LAL_MTSUN_SI;
  /*printf("n = %e, m = %e \n",n,m);*/
  /* FIXME I MIGHT HAVE THE SIGN WRONG */
  output = -1.0 / 2.0 * (0.0
    + 6.77125 * x(I,m,fl,fh,1+xpow,1)
    - 7.52771 / n / n / n * x(I,m,fl,fh,0+xpow,0)
    + 19.1028 / n / n / n * x(I,m,fl,fh,1+xpow,1)
    - 95.8736 / n / n / n * x(I,m,fl,fh,2+xpow,2) 
    - 3./64. / n / n / n / pow(LAL_PI,5./3.) * x(I,m,fl,fh,-5+xpow,-5)  
    - 2496480. / 10838016. / n / n / n / LAL_PI * x(I,m,fl,fh,-3+xpow,-3) 
    - 15293365. / 10838016./n/n/n/pow(LAL_PI,1./3.) * x(I,m,fl,fh,-1+xpow,-1)
    + 8128512. / 10838016./n/n/n*pow(LAL_PI,1./3.) * x(I,m,fl,fh,-2+xpow,-2)
    + 577.709 / n / n / n * (log(1)-log(LAL_PI*m)/3.) * x(I,m,fl,fh,3+xpow,3) 
    - 577.709/3. / n / n / n * lx(I,m,fl,fh,3+xpow,3) 
    + 22.5831 / n / n / n * (log(1)-log(LAL_PI*m)/3.) * x(I,m,fl,fh,xpow,0) 
    - 22.5831/3. / n / n / n * lx(I,m,fl,fh,xpow,0) 
    - 22.5831 / n / n / n * (log(1)-log(LAL_PI*m)/3.) * x(I,m,fl,fh,1+xpow,1)
    - 22.5831/3. / n / n / n * lx(I,m,fl,fh,1+xpow,1)

    );

  /*printf("Unnormalized metric component %e\n",output);*/
  return output ;

  }

REAL8 XLALComputeLALIMRBankMetricMassEta(REAL8 mass1, REAL8 mass2,
         REAL8 fl, REAL8 fh, INT4 xpow, LALIMRBankCumulativeNoiseMoments *I)
  {
  REAL8 output = 0;
  REAL8 n = mass1*mass2/(mass1+mass2)/(mass1+mass2);
  REAL8 m = (mass1+mass2)*LAL_MTSUN_SI;
  /*printf("n = %e, m = %e \n",n,m);*/
  output = -1.0 / 2.0 / m * (0.0
    - 0.503603 * x(I,m,fl,fh,1+xpow,1)
    + 10.3128 * x(I,m,fl,fh,2+xpow,2)
    + 3.76385 / n / n * x(I,m,fl,fh,xpow,0)
    - 6.91502 / n / n * x(I,m,fl,fh,1+xpow,1)
    + 31.9579 / n / n * x(I,m,fl,fh,2+xpow,2)
    + 2.25708 * n * x(I,m,fl,fh,1+xpow,1)
    - 5./128. / n / n / pow(LAL_PI,5./3.) * x(I,m,fl,fh,-5+xpow,-5)
    - 7489440./65028096/n/n/LAL_PI * x(I,m,fl,fh,-3+xpow,-3)
    - 5.*3058673/65028096/n/n/pow(LAL_PI,1./3.) * x(I,m,fl,fh,-1+xpow,-1)
    + 5.*4353552/65028096/pow(LAL_PI,1./3.) * x(I,m,fl,fh,-3+xpow,-3)
    + 16257024/65028096/n/n*pow(LAL_PI,1./3.) * x(I,m,fl,fh,-2+xpow,-2)
    + 96.2849/n/n * x(I,m,fl,fh,1+xpow,1)
    - 288.855/n/n*(log(0.682784)-log(m)/3.) * x(I,m,fl,fh,3+xpow,3)
    + 288.855/3./n/n * lx(I,m,fl,fh,3+xpow,3)
    + 3.73122/n/n*(log(1)-log(LAL_PI*m)/3.) * x(I,m,fl,fh,1+xpow,1)
    - 3.73122/3./n/n * lx(I,m,fl,fh,1+xpow,1)
    );
  /* WARNING THIS RESCALES TO SOLAR MASS! */
  return output*LAL_MTSUN_SI;
  }

int XLALComputeLALIMRBankMetric(REAL8 mass1, REAL8 mass2, LALIMRBankCumulativeNoiseMoments *moments, LALIMRBankMetric *metric)
  {
  REAL8 q = 0;
  REAL8 m = (mass1+mass2);
  REAL8 fm = 0;
  REAL8 fr = 0;
  REAL8 fl = moments->flow;

  if (mass1 < mass2) q = mass1/mass2;
  else q = mass2/mass1;

  fm = (0.8 * q * q * q - 2.6 * q * q + 2.8 * q + 1.0) / LAL_PI / m / LAL_MTSUN_SI / pow(6,3.0/2.0);
  fr = (0.099*q*q*q -0.27*q*q + 0.25*q +0.19) / LAL_PI / m / LAL_MTSUN_SI;
  /*fm = 1.0 /  LAL_PI / m / LAL_MTSUN_SI / pow(3,3.0/2.0);
  fr = 1.0 /  LAL_PI / m / LAL_MTSUN_SI / pow(3,3.0/2.0);*/
  /* fr = */
  /*printf("calling xomputeIMRmetricMAsMass with fl %f fm %f fr %f.... q %f m %f \n",fl,fm,fr,q,m);*/

  metric->data[0][0] = 
          ( XLALComputeLALIMRBankMetricMassMass(mass1,mass2,fl,fm,-7,moments) 
          +  XLALComputeLALIMRBankMetricMassMass(mass1,mass2,fm,fr,-4,moments) )
          / (x(moments,1,fl,fm,-7,1) + x(moments,1,fm,fr,-4,1));

  metric->data[1][1] =
          ( XLALComputeLALIMRBankMetricEtaEta(mass1,mass2,fl,fm,-7,moments)
          +  XLALComputeLALIMRBankMetricEtaEta(mass1,mass2,fm,fr,-4,moments) )
          / (x(moments,1,fl,fm,-7,1) + x(moments,1,fm,fr,-4,1));

  metric->data[0][1] = metric->data[1][0] = 
          ( XLALComputeLALIMRBankMetricMassEta(mass1,mass2,fl,fm,-7,moments)
          +  XLALComputeLALIMRBankMetricMassEta(mass1,mass2,fm,fr,-4,moments) )
          / (x(moments,1,fl,fm,-7,1) + x(moments,1,fm,fr,-4,1));

  return 0;

  }

static REAL8 eta(REAL8 m1, REAL8 m2)
  {
  return m1*m2/(m1+m2)/(m1+m2);
  }

static REAL8 XLALComputeNumberOfIMRTemplatesInSquareLALIMRBankMassRegion(REAL8 mbox[3], 
                      REAL8 mm, LALIMRBankCumulativeNoiseMoments *I)
  {
  LALIMRBankMetric metric;
  REAL8 volumeElem = 0;
  REAL8 volume = 0;
  REAL8 totalMassDiff = 2*mbox[2];
  REAL8 cm1 = mbox[0];/* + mbox[2]/2.0;*/
  REAL8 cm2 = mbox[1];/* + mbox[2]/2.0;*/
  REAL8 m1 = mbox[0];
  REAL8 m2 = mbox[1];
  REAL8 size = mbox[2];
  REAL8 eta1 = eta(m1,m2);
  REAL8 eta2 = eta(m1+size,m2);
  REAL8 eta3 = eta(m1,m2+size);
  REAL8 eta4 = eta(m1+size,m2+size);
  REAL8 maxEta = eta1;
  REAL8 minEta = eta1;
  REAL8 etaDiff = 0;

  if (eta2 > maxEta) maxEta = eta2;
  if (eta3 > maxEta) maxEta = eta3;
  if (eta4 > maxEta) maxEta = eta4;
  if (eta2 < minEta) minEta = eta2;
  if (eta3 < minEta) minEta = eta3;
  if (eta4 < minEta) minEta = eta4;

   etaDiff = maxEta-minEta;
  
  XLALComputeLALIMRBankMetric(cm1,cm2,I,&metric);
  volumeElem = sqrt(fabs(metric.data[0][0]*metric.data[1][1]-metric.data[0][1]*metric.data[0][1]));
  
  volume = etaDiff * totalMassDiff;
  /*printf("massdiff %f etadiff %f cm1 %f cm3 %f volumeElem %f volume %f\n",
         totalMassDiff,etaDiff,cm1,cm2,volumeElem,volume);*/
  /*FIXME what should this be? I have fudged it a bit */	 
  return volume * volumeElem / (20.0 * mm);
  }

static LALIMRBankMassRegion * createLALIMRBankMassRegion(REAL8 mass1, REAL8 mass2, REAL8 size)
  {
  LALIMRBankMassRegion *region = calloc(1, sizeof(LALIMRBankMassRegion));
  region->mbox[0] = mass1;
  region->mbox[1] = mass2;
  region->mbox[2] = size;
  return region;
  }

static int appendtotail(LALIMRBankMassRegion *elem, LALIMRBankMassRegion **tail)
  {
  (*tail)->next = elem;
  *tail = elem;
  return 0;
  }

static int divideAndConquer(LALIMRBankMassRegion *list, 
                            LALIMRBankMassRegion **tail)

  {
  REAL8 m1 = list->mbox[0];
  REAL8 m2 = list->mbox[1];
  REAL8 size = list->mbox[2];
  REAL8 newsize = size/2.;
  /* No checking if it is null  BAD */
  appendtotail( createLALIMRBankMassRegion(m1,m2,newsize), tail);
  appendtotail( createLALIMRBankMassRegion(m1+newsize,m2,newsize), tail);
  appendtotail( createLALIMRBankMassRegion(m1,m2+newsize,newsize), tail);
  appendtotail( createLALIMRBankMassRegion(m1+newsize,m2+newsize,newsize), tail);
  return 0;
  }

static int addtemplates(REAL8 mbox[3], REAL8 numTmps, 
                        InspiralCoarseBankIn *in, SnglInspiralTable **head,
			gsl_rng *r,
			UINT4 *ntiles)
  {
  UINT4 num = floor(numTmps+0.50);
  REAL8 m1 = 0;
  REAL8 m2 = 0;
  REAL8 mtot = mbox[0]+mbox[1];
  REAL8 size = mbox[2];
  UINT4 i = 0;
  /* check to see if it is the bottom half of the mass/mass plane */
  /* and that is falls within the boundaries */

  if ( (mbox[1] < mbox[0]) && (mtot >= in->MMin) && (mtot <= in->MMax) )
  {
    for (i = 0; i < num; i++)
      {
      /*head = (SnglInspiralTable *) LALCalloc(1,sizeof(SnglInspiralTable));*/
      m1 = mbox[0]+size * gsl_rng_uniform_pos(r); /*( (double)rand() / ((double)(RAND_MAX)+(double)(1)) );*/
      m2 = mbox[1]+size * gsl_rng_uniform_pos(r); /*( (double)rand()/((double)(RAND_MAX)+(double)(1)) );*/
      (*head)->mass1 = m1;
      (*head)->mass2 = m2;
      (*head)->eta = m1*m2/(m1+m2)/(m1+m2);
      (*head)->mchirp = pow(m1*m2,0.6)/pow(m1+m2,0.2);
      (*head)->next = (SnglInspiralTable *) LALCalloc(1,sizeof(SnglInspiralTable));
      (*head) = (*head)->next;
      (*ntiles)++;
      }
    }
  return 0;
  }


static int checkNumberOfTemplates(LALIMRBankMassRegion *list, 
                                  LALIMRBankMassRegion **tail, 
				  InspiralCoarseBankIn *in, 
				  LALIMRBankCumulativeNoiseMoments *I, 
				  SnglInspiralTable **head,
				  gsl_rng *r,
				  UINT4 *ntiles)
  {
  REAL8 mm = 1.0-in->mmCoarse;
  REAL8 numTmps = XLALComputeNumberOfIMRTemplatesInSquareLALIMRBankMassRegion(
                  list->mbox, mm, I);
  
  if (numTmps > 4)  divideAndConquer(list,tail);
  else addtemplates(list->mbox,numTmps,in,head,r,ntiles);
  return 0;
  }

static int destroyregionlist(LALIMRBankMassRegion *head)
  {
  LALIMRBankMassRegion *tmp = NULL;

  while(head)
    {
    tmp = head;
    head = head->next;
    free(tmp);
    }
  return 0;
  }

static int cleanuplast(SnglInspiralTable *head)
  {
  SnglInspiralTable *tmp = NULL;
  while(head)
    {
    tmp = head->next;
    if (tmp && !(tmp->next))
      {
      head->next = NULL;
      LALFree(tmp);
      }
    head = head->next;
    }
  }

/* This is the main function that */
INT4 XLALTileLALIMRBankMassRegion(InspiralCoarseBankIn *in, SnglInspiralTable **first)
  {

  /* Just a bit of padding for the big box.  but the templates will only be
   * in the specified region */

  REAL8 mass1 = 0.9*in->mMin;
  REAL8 mass2 = 0.9*in->mMin;
  REAL8 size = 1.1*in->mMax - 0.9*in->mMin;
  REAL8 mm = 1.0 - in->mmCoarse;
  REAL8 flow = in->fLower;
  UINT4 ntiles = 0;
  LALIMRBankCumulativeNoiseMoments moments;
  LALIMRBankMassRegion *head = createLALIMRBankMassRegion(mass1,mass2,size);
  LALIMRBankMassRegion *list = head;
  LALIMRBankMassRegion *tail = head;
  gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
  SnglInspiralTable *first_return = NULL;
  /*if (*first) LALfree(*first);*/
  /* *first = NULL;*/
  *first = (SnglInspiralTable *) LALCalloc(1,sizeof(SnglInspiralTable));
  first_return = *first;
  XLALCreateLALIMRBankCumulativeNoiseMoments(&moments,&(in->shf),flow);

  while(list)
    {
    checkNumberOfTemplates(list,&tail,in,&moments,first,r,&ntiles);
    list = list->next;
    }
  XLALDestroyLALIMRBankCumulativeNoiseMoments(&moments);
  destroyregionlist(head);
  gsl_rng_free(r);
  *first = first_return;
  /* remove pesky zero mass in last template */
  cleanuplast(*first);
  *first = first_return;
  return ntiles;
  }


