#include <lal/LALInspiralBank.h>

NRCSID (LALNORMALISEC, "$Id$");

void
LALNormalise (
   LALStatus *status, 
   REAL4Vector *in, 
   REAL8 *norm,
   REAL8Vector psd) 
{
  INT4 i, n;

  INITSTATUS (status, "LALNormalise", LALNORMALISEC);
  if (psd.length != in->length/2+1) {
  fprintf(stderr, "LALNormalise: Incompatible lengths %d, %d\n",
     psd.length, in->length);
     exit(1);
  }

  n = in->length;
  *norm = 0;
  for (i=2; i<n; i+=2)
     if (psd.data[i/2+1])
        *norm += (pow(in->data[i],2.)+pow(in->data[i+1],2.))/psd.data[i/2+1];
  *norm = sqrt(*norm);
  
  i=n;
  while (i--) *(in->data+i) = *(in->data+i) / *norm;
  RETURN(status);
}

