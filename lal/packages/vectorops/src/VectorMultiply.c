#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/VectorOps.h>

NRCSID (VECTORMULTIPLYC, "$Id$");

void
LALCCVectorDivide (
    LALStatus            *status,
    COMPLEX8Vector       *out,
    const COMPLEX8Vector *in1,
    const COMPLEX8Vector *in2
    )
{
  COMPLEX8 *a;
  COMPLEX8 *b;
  COMPLEX8 *c;
  INT4      n;

  INITSTATUS (status, "LALCCVectorDivide", VECTORMULTIPLYC);

  ASSERT (out, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in1, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (in2, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (out->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in1->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in2->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (out->length > 0, status, VECTOROPS_ESIZE, VECTOROPS_MSGESIZE);
  ASSERT (in1->length == out->length, status,
          VECTOROPS_ESZMM, VECTOROPS_MSGESZMM);
  ASSERT (in2->length == out->length, status,
          VECTOROPS_ESZMM, VECTOROPS_MSGESZMM);

  a = in1->data;
  b = in2->data;
  c = out->data;
  n = out->length;

  while (n-- > 0)
  {
    REAL4 ar = a->re;
    REAL4 ai = a->im;
    REAL4 br = b->re;
    REAL4 bi = b->im;

    if (fabs(br) > fabs(bi))
    {
      REAL4 ratio = bi/br;
      REAL4 denom = br + ratio*bi;

      c->re = (ar + ratio*ai)/denom;
      c->im = (ai - ratio*ar)/denom;
    }
    else
    {
      REAL4 ratio = br/bi;
      REAL4 denom = bi + ratio*br;

      c->re = (ar*ratio + ai)/denom;
      c->im = (ai*ratio - ar)/denom;
    }

    ++a;
    ++b;
    ++c;
  }

  RETURN (status);
}


void
LALZZVectorDivide (
    LALStatus             *status,
    COMPLEX16Vector       *out,
    const COMPLEX16Vector *in1,
    const COMPLEX16Vector *in2
    )
{
  COMPLEX16 *a;
  COMPLEX16 *b;
  COMPLEX16 *c;
  INT4       n;

  INITSTATUS (status, "LALZZVectorDivide", VECTORMULTIPLYC);

  ASSERT (out, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in1, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (in2, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (out->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in1->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in2->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (out->length > 0, status, VECTOROPS_ESIZE, VECTOROPS_MSGESIZE);
  ASSERT (in1->length == out->length, status,
          VECTOROPS_ESZMM, VECTOROPS_MSGESZMM);
  ASSERT (in2->length == out->length, status,
          VECTOROPS_ESZMM, VECTOROPS_MSGESZMM);

  a = in1->data;
  b = in2->data;
  c = out->data;
  n = out->length;

  while (n-- > 0)
  {
    REAL8 ar = a->re;
    REAL8 ai = a->im;
    REAL8 br = b->re;
    REAL8 bi = b->im;

    if (fabs(br) > fabs(bi))
    {
      REAL8 ratio = bi/br;
      REAL8 denom = br + ratio*bi;

      c->re = (ar + ratio*ai)/denom;
      c->im = (ai - ratio*ar)/denom;
    }
    else
    {
      REAL8 ratio = br/bi;
      REAL8 denom = bi + ratio*br;

      c->re = (ar*ratio + ai)/denom;
      c->im = (ai*ratio - ar)/denom;
    }

    ++a;
    ++b;
    ++c;
  }

  RETURN (status);
}



void
LALCCVectorMultiply (
    LALStatus            *status,
    COMPLEX8Vector       *out,
    const COMPLEX8Vector *in1,
    const COMPLEX8Vector *in2
    )
{
  COMPLEX8 *a;
  COMPLEX8 *b;
  COMPLEX8 *c;
  INT4      n;

  INITSTATUS (status, "LALCCVectorMultiply", VECTORMULTIPLYC);

  ASSERT (out, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in1, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (in2, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (out->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in1->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in2->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (out->length > 0, status, VECTOROPS_ESIZE, VECTOROPS_MSGESIZE);
  ASSERT (in1->length == out->length, status,
          VECTOROPS_ESZMM, VECTOROPS_MSGESZMM);
  ASSERT (in2->length == out->length, status,
          VECTOROPS_ESZMM, VECTOROPS_MSGESZMM);

  a = in1->data;
  b = in2->data;
  c = out->data;
  n = out->length;

  while (n-- > 0)
  {
    REAL4 ar = a->re;
    REAL4 ai = a->im;
    REAL4 br = b->re;
    REAL4 bi = b->im;

    c->re = ar*br - ai*bi;
    c->im = ar*bi + ai*br;

    ++a;
    ++b;
    ++c;
  }

  RETURN (status);
}


void
LALZZVectorMultiply (
    LALStatus             *status,
    COMPLEX16Vector       *out,
    const COMPLEX16Vector *in1,
    const COMPLEX16Vector *in2
    )
{
  COMPLEX16 *a;
  COMPLEX16 *b;
  COMPLEX16 *c;
  INT4       n;

  INITSTATUS (status, "LALZZVectorMultiply", VECTORMULTIPLYC);

  ASSERT (out, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in1, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (in2, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (out->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in1->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in2->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (out->length > 0, status, VECTOROPS_ESIZE, VECTOROPS_MSGESIZE);
  ASSERT (in1->length == out->length, status,
          VECTOROPS_ESZMM, VECTOROPS_MSGESZMM);
  ASSERT (in2->length == out->length, status,
          VECTOROPS_ESZMM, VECTOROPS_MSGESZMM);

  a = in1->data;
  b = in2->data;
  c = out->data;
  n = out->length;

  while (n-- > 0)
  {
    REAL8 ar = a->re;
    REAL8 ai = a->im;
    REAL8 br = b->re;
    REAL8 bi = b->im;

    c->re = ar*br - ai*bi;
    c->im = ar*bi + ai*br;

    ++a;
    ++b;
    ++c;
  }

  RETURN (status);
}


void
LALCCVectorMultiplyConjugate (
    LALStatus            *status,
    COMPLEX8Vector       *out,
    const COMPLEX8Vector *in1,
    const COMPLEX8Vector *in2
    )
{
  COMPLEX8 *a;
  COMPLEX8 *b;
  COMPLEX8 *c;
  INT4      n;

  INITSTATUS (status, "LALCCVectorMultiplyConjugate", VECTORMULTIPLYC);

  ASSERT (out, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in1, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (in2, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (out->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in1->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in2->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (out->length > 0, status, VECTOROPS_ESIZE, VECTOROPS_MSGESIZE);
  ASSERT (in1->length == out->length, status,
          VECTOROPS_ESZMM, VECTOROPS_MSGESZMM);
  ASSERT (in2->length == out->length, status,
          VECTOROPS_ESZMM, VECTOROPS_MSGESZMM);

  a = in1->data;
  b = in2->data;
  c = out->data;
  n = out->length;

  while (n-- > 0)
  {
    REAL4 ar = a->re;
    REAL4 ai = a->im;
    REAL4 br = b->re;
    REAL4 bi = b->im;

    c->re = ar*br + ai*bi;
    c->im = ai*br - ar*bi;

    ++a;
    ++b;
    ++c;
  }

  RETURN (status);
}


void
LALZZVectorMultiplyConjugate (
    LALStatus             *status,
    COMPLEX16Vector       *out,
    const COMPLEX16Vector *in1,
    const COMPLEX16Vector *in2
    )
{
  COMPLEX16 *a;
  COMPLEX16 *b;
  COMPLEX16 *c;
  INT4       n;

  INITSTATUS (status, "LALZZVectorMultiplyConjugate", VECTORMULTIPLYC);

  ASSERT (out, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in1, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (in2, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (out->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in1->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in2->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (out->length > 0, status, VECTOROPS_ESIZE, VECTOROPS_MSGESIZE);
  ASSERT (in1->length == out->length, status,
          VECTOROPS_ESZMM, VECTOROPS_MSGESZMM);
  ASSERT (in2->length == out->length, status,
          VECTOROPS_ESZMM, VECTOROPS_MSGESZMM);

  a = in1->data;
  b = in2->data;
  c = out->data;
  n = out->length;

  while (n-- > 0)
  {
    REAL8 ar = a->re;
    REAL8 ai = a->im;
    REAL8 br = b->re;
    REAL8 bi = b->im;

    c->re = ar*br + ai*bi;
    c->im = ai*br - ar*bi;

    ++a;
    ++b;
    ++c;
  }

  RETURN (status);
}


void
LALSCVectorMultiply (
    LALStatus            *status,
    COMPLEX8Vector       *out,
    const REAL4Vector    *in1,
    const COMPLEX8Vector *in2
    )
{
  REAL4    *a;
  COMPLEX8 *b;
  COMPLEX8 *c;
  INT4      n;

  INITSTATUS (status, "LALSCVectorMultiply", VECTORMULTIPLYC);

  ASSERT (out, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in1, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (in2, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (out->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in1->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in2->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (out->length > 0, status, VECTOROPS_ESIZE, VECTOROPS_MSGESIZE);
  ASSERT (in1->length == out->length, status,
          VECTOROPS_ESZMM, VECTOROPS_MSGESZMM);
  ASSERT (in2->length == out->length, status,
          VECTOROPS_ESZMM, VECTOROPS_MSGESZMM);

  a = in1->data;
  b = in2->data;
  c = out->data;
  n = out->length;

  while (n-- > 0)
  {
    REAL4 fac = *a;
    REAL4 br  = b->re;
    REAL4 bi  = b->im;

    c->re = fac*br;
    c->im = fac*bi;

    ++a;
    ++b;
    ++c;
  }

  RETURN (status);
}


void
LALDZVectorMultiply (
    LALStatus             *status,
    COMPLEX16Vector       *out,
    const REAL8Vector     *in1,
    const COMPLEX16Vector *in2
    )
{
  REAL8     *a;
  COMPLEX16 *b;
  COMPLEX16 *c;
  INT4       n;

  INITSTATUS (status, "LALDZVectorMultiply", VECTORMULTIPLYC);

  ASSERT (out, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in1, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (in2, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (out->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in1->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in2->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (out->length > 0, status, VECTOROPS_ESIZE, VECTOROPS_MSGESIZE);
  ASSERT (in1->length == out->length, status,
          VECTOROPS_ESZMM, VECTOROPS_MSGESZMM);
  ASSERT (in2->length == out->length, status,
          VECTOROPS_ESZMM, VECTOROPS_MSGESZMM);

  a = in1->data;
  b = in2->data;
  c = out->data;
  n = out->length;

  while (n-- > 0)
  {
    REAL8 fac = *a;
    REAL8 br  = b->re;
    REAL8 bi  = b->im;

    c->re = fac*br;
    c->im = fac*bi;

    ++a;
    ++b;
    ++c;
  }

  RETURN (status);
}


void
LALSSVectorMultiply (
    LALStatus            *status,
    REAL4Vector          *out,
    const REAL4Vector    *in1,
    const REAL4Vector    *in2
    )
{
  REAL4 *a;
  REAL4 *b;
  REAL4 *c;
  INT4   n;

  INITSTATUS (status, "LALSSVectorMultiply", VECTORMULTIPLYC);

  ASSERT (out, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in1, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (in2, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (out->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in1->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in2->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (out->length > 0, status, VECTOROPS_ESIZE, VECTOROPS_MSGESIZE);
  ASSERT (in1->length == out->length, status,
          VECTOROPS_ESZMM, VECTOROPS_MSGESZMM);
  ASSERT (in2->length == out->length, status,
          VECTOROPS_ESZMM, VECTOROPS_MSGESZMM);

  a = in1->data;
  b = in2->data;
  c = out->data;
  n = out->length;

  while (n-- > 0)
  {
    *c++ = (*a++)*(*b++);
  }

  RETURN (status);
}


void
LALDDVectorMultiply (
    LALStatus            *status,
    REAL8Vector          *out,
    const REAL8Vector    *in1,
    const REAL8Vector    *in2
    )
{
  REAL8 *a;
  REAL8 *b;
  REAL8 *c;
  INT4   n;

  INITSTATUS (status, "LALDDVectorMultiply", VECTORMULTIPLYC);

  ASSERT (out, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in1, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (in2, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (out->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in1->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in2->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (out->length > 0, status, VECTOROPS_ESIZE, VECTOROPS_MSGESIZE);
  ASSERT (in1->length == out->length, status,
          VECTOROPS_ESZMM, VECTOROPS_MSGESZMM);
  ASSERT (in2->length == out->length, status,
          VECTOROPS_ESZMM, VECTOROPS_MSGESZMM);

  a = in1->data;
  b = in2->data;
  c = out->data;
  n = out->length;

  while (n-- > 0)
  {
    *c++ = (*a++)*(*b++);
  }

  RETURN (status);
}

