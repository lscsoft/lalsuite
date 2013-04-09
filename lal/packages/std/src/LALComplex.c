/** \cond DONT_DOXYGEN */
/*
 *
 * NOTE: THIS CODE WAS AUTOMATICALLY GENERATED FROM GSL CODE
 * ORIGINAL COPYRIGHT IS BELOW.
 *
 *
 */

/* complex/math.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Jorma Olavi T<E4>htinen, Brian Gou
gh
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
.
 */


#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <math.h>
#define LAL_NO_COMPLEX_MACROS
#include <lal/LALComplex.h>


COMPLEX16
XLALCOMPLEX16Rect (REAL8 x, REAL8 y)
{
  COMPLEX16 z;
  do {(&z)->re=(x); (&z)->im=(y);} while(0);
  return z;
}

COMPLEX16
XLALCOMPLEX16Polar (REAL8 r, REAL8 theta)
{
  COMPLEX16 z;
  do {(&z)->re=(r * cos (theta)); (&z)->im=(r * sin (theta));} while(0);
  return z;
}





REAL8
XLALCOMPLEX16Arg (COMPLEX16 z)
{
  REAL8 x = ((z).re);
  REAL8 y = ((z).im);

  if (x == 0.0 && y == 0.0)
    {
      return 0;
    }

  return atan2 (y, x);
}

REAL8
XLALCOMPLEX16Abs (COMPLEX16 z)
{
  return hypot (((z).re), ((z).im));
}

REAL8
XLALCOMPLEX16Abs2 (COMPLEX16 z)
{
  REAL8 x = ((z).re);
  REAL8 y = ((z).im);

  return (x * x + y * y);
}

REAL8
XLALCOMPLEX16LogAbs (COMPLEX16 z)
{
  REAL8 xabs = fabs (((z).re));
  REAL8 yabs = fabs (((z).im));
  REAL8 max, u;

  if (xabs >= yabs)
    {
      max = xabs;
      u = yabs / xabs;
    }
  else
    {
      max = yabs;
      u = xabs / yabs;
    }



  return log (max) + 0.5 * log1p (u * u);
}






COMPLEX16
XLALCOMPLEX16Add (COMPLEX16 a, COMPLEX16 b)
{
  REAL8 ar = ((a).re), ai = ((a).im);
  REAL8 br = ((b).re), bi = ((b).im);

  COMPLEX16 z;
  do {(&z)->re=(ar + br); (&z)->im=(ai + bi);} while(0);
  return z;
}

COMPLEX16
XLALCOMPLEX16AddReal (COMPLEX16 a, REAL8 x)
{
  COMPLEX16 z;
  do {(&z)->re=(((a).re) + x); (&z)->im=(((a).im));} while(0);
  return z;
}

COMPLEX16
XLALCOMPLEX16AddImag (COMPLEX16 a, REAL8 y)
{
  COMPLEX16 z;
  do {(&z)->re=(((a).re)); (&z)->im=(((a).im) + y);} while(0);
  return z;
}


COMPLEX16
XLALCOMPLEX16Sub (COMPLEX16 a, COMPLEX16 b)
{
  REAL8 ar = ((a).re), ai = ((a).im);
  REAL8 br = ((b).re), bi = ((b).im);

  COMPLEX16 z;
  do {(&z)->re=(ar - br); (&z)->im=(ai - bi);} while(0);
  return z;
}

COMPLEX16
XLALCOMPLEX16SubReal (COMPLEX16 a, REAL8 x)
{
  COMPLEX16 z;
  do {(&z)->re=(((a).re) - x); (&z)->im=(((a).im));} while(0);
  return z;
}

COMPLEX16
XLALCOMPLEX16SubImag (COMPLEX16 a, REAL8 y)
{
  COMPLEX16 z;
  do {(&z)->re=(((a).re)); (&z)->im=(((a).im) - y);} while(0);
  return z;
}

COMPLEX16
XLALCOMPLEX16Mul (COMPLEX16 a, COMPLEX16 b)
{
  REAL8 ar = ((a).re), ai = ((a).im);
  REAL8 br = ((b).re), bi = ((b).im);

  COMPLEX16 z;
  do {(&z)->re=(ar * br - ai * bi); (&z)->im=(ar * bi + ai * br);} while(0);
  return z;
}

COMPLEX16
XLALCOMPLEX16MulReal (COMPLEX16 a, REAL8 x)
{
  COMPLEX16 z;
  do {(&z)->re=(x * ((a).re)); (&z)->im=(x * ((a).im));} while(0);
  return z;
}

COMPLEX16
XLALCOMPLEX16MulImag (COMPLEX16 a, REAL8 y)
{
  COMPLEX16 z;
  do {(&z)->re=(-y * ((a).im)); (&z)->im=(y * ((a).re));} while(0);
  return z;
}

COMPLEX16
XLALCOMPLEX16Div (COMPLEX16 a, COMPLEX16 b)
{
  REAL8 ar = ((a).re), ai = ((a).im);
  REAL8 br = ((b).re), bi = ((b).im);

  REAL8 s = 1.0 / XLALCOMPLEX16Abs (b);

  REAL8 sbr = s * br;
  REAL8 sbi = s * bi;

  REAL8 zr = (ar * sbr + ai * sbi) * s;
  REAL8 zi = (ai * sbr - ar * sbi) * s;

  COMPLEX16 z;
  do {(&z)->re=(zr); (&z)->im=(zi);} while(0);
  return z;
}

COMPLEX16
XLALCOMPLEX16DivReal (COMPLEX16 a, REAL8 x)
{
  COMPLEX16 z;
  do {(&z)->re=(((a).re) / x); (&z)->im=(((a).im) / x);} while(0);
  return z;
}

COMPLEX16
XLALCOMPLEX16DivImag (COMPLEX16 a, REAL8 y)
{
  COMPLEX16 z;
  do {(&z)->re=(((a).im) / y); (&z)->im=(- ((a).re) / y);} while(0);
  return z;
}

COMPLEX16
XLALCOMPLEX16Conjugate (COMPLEX16 a)
{
  COMPLEX16 z;
  do {(&z)->re=(((a).re)); (&z)->im=(-((a).im));} while(0);
  return z;
}

COMPLEX16
XLALCOMPLEX16Negative (COMPLEX16 a)
{
  COMPLEX16 z;
  do {(&z)->re=(-((a).re)); (&z)->im=(-((a).im));} while(0);
  return z;
}

COMPLEX16
XLALCOMPLEX16Inverse (COMPLEX16 a)
{
  REAL8 s = 1.0 / XLALCOMPLEX16Abs (a);

  COMPLEX16 z;
  do {(&z)->re=((((a).re) * s) * s); (&z)->im=(-(((a).im) * s) * s);} while(0);
  return z;
}





COMPLEX16
XLALCOMPLEX16Sqrt (COMPLEX16 a)
{
  COMPLEX16 z;

  if (((a).re) == 0.0 && ((a).im) == 0.0)
    {
      do {(&z)->re=(0); (&z)->im=(0);} while(0);
    }
  else
    {
      REAL8 x = fabs (((a).re));
      REAL8 y = fabs (((a).im));
      REAL8 w;

      if (x >= y)
        {
          REAL8 t = y / x;
          w = sqrt (x) * sqrt (0.5 * (1.0 + sqrt (1.0 + t * t)));
        }
      else
        {
          REAL8 t = x / y;
          w = sqrt (y) * sqrt (0.5 * (t + sqrt (1.0 + t * t)));
        }

      if (((a).re) >= 0.0)
        {
          REAL8 ai = ((a).im);
          do {(&z)->re=(w); (&z)->im=(ai / (2.0 * w));} while(0);
        }
      else
        {
          REAL8 ai = ((a).im);
          REAL8 vi = (ai >= 0) ? w : -w;
          do {(&z)->re=(ai / (2.0 * vi)); (&z)->im=(vi);} while(0);
        }
    }

  return z;
}

COMPLEX16
XLALCOMPLEX16SqrtReal (REAL8 x)
{
  COMPLEX16 z;

  if (x >= 0)
    {
      do {(&z)->re=(sqrt (x)); (&z)->im=(0.0);} while(0);
    }
  else
    {
      do {(&z)->re=(0.0); (&z)->im=(sqrt (-x));} while(0);
    }

  return z;
}

COMPLEX16
XLALCOMPLEX16Exp (COMPLEX16 a)
{
  REAL8 rho = exp (((a).re));
  REAL8 theta = ((a).im);

  COMPLEX16 z;
  do {(&z)->re=(rho * cos (theta)); (&z)->im=(rho * sin (theta));} while(0);
  return z;
}

COMPLEX16
XLALCOMPLEX16Pow (COMPLEX16 a, COMPLEX16 b)
{
  COMPLEX16 z;

  if (((a).re) == 0 && ((a).im) == 0.0)
    {
      do {(&z)->re=(0.0); (&z)->im=(0.0);} while(0);
    }
  else
    {
      REAL8 logr = XLALCOMPLEX16LogAbs (a);
      REAL8 theta = XLALCOMPLEX16Arg (a);

      REAL8 br = ((b).re), bi = ((b).im);

      REAL8 rho = exp (logr * br - bi * theta);
      REAL8 beta = theta * br + bi * logr;

      do {(&z)->re=(rho * cos (beta)); (&z)->im=(rho * sin (beta));} while(0);
    }

  return z;
}

COMPLEX16
XLALCOMPLEX16PowReal (COMPLEX16 a, REAL8 b)
{
  COMPLEX16 z;

  if (((a).re) == 0 && ((a).im) == 0)
    {
      do {(&z)->re=(0); (&z)->im=(0);} while(0);
    }
  else
    {
      REAL8 logr = XLALCOMPLEX16LogAbs (a);
      REAL8 theta = XLALCOMPLEX16Arg (a);
      REAL8 rho = exp (logr * b);
      REAL8 beta = theta * b;
      do {(&z)->re=(rho * cos (beta)); (&z)->im=(rho * sin (beta));} while(0);
    }

  return z;
}

COMPLEX16
XLALCOMPLEX16Log (COMPLEX16 a)
{
  REAL8 logr = XLALCOMPLEX16LogAbs (a);
  REAL8 theta = XLALCOMPLEX16Arg (a);

  COMPLEX16 z;
  do {(&z)->re=(logr); (&z)->im=(theta);} while(0);
  return z;
}

COMPLEX16
XLALCOMPLEX16Log10 (COMPLEX16 a)
{
  return XLALCOMPLEX16MulReal (XLALCOMPLEX16Log (a), 1 / log (10.));
}

COMPLEX16
XLALCOMPLEX16LogB (COMPLEX16 a, COMPLEX16 b)
{
  return XLALCOMPLEX16Div (XLALCOMPLEX16Log (a), XLALCOMPLEX16Log (b));
}





COMPLEX16
XLALCOMPLEX16Sin (COMPLEX16 a)
{
  REAL8 R = ((a).re), I = ((a).im);

  COMPLEX16 z;

  if (I == 0.0)
    {


      do {(&z)->re=(sin (R)); (&z)->im=(0.0);} while(0);
    }
  else
    {
      do {(&z)->re=(sin (R) * cosh (I)); (&z)->im=(cos (R) * sinh (I));} while(0);
    }

  return z;
}

COMPLEX16
XLALCOMPLEX16Cos (COMPLEX16 a)
{
  REAL8 R = ((a).re), I = ((a).im);

  COMPLEX16 z;

  if (I == 0.0)
    {


      do {(&z)->re=(cos (R)); (&z)->im=(0.0);} while(0);
    }
  else
    {
      do {(&z)->re=(cos (R) * cosh (I)); (&z)->im=(sin (R) * sinh (-I));} while(0);
    }

  return z;
}

COMPLEX16
XLALCOMPLEX16Tan (COMPLEX16 a)
{
  REAL8 R = ((a).re), I = ((a).im);

  COMPLEX16 z;

  if (fabs (I) < 1)
    {
      REAL8 D = pow (cos (R), 2.0) + pow (sinh (I), 2.0);

      do {(&z)->re=(0.5 * sin (2 * R) / D); (&z)->im=(0.5 * sinh (2 * I) / D);} while(0);
    }
  else
    {
      REAL8 u = exp (-I);
      REAL8 C = 2 * u / (1 - pow (u, 2.0));
      REAL8 D = 1 + pow (cos (R), 2.0) * pow (C, 2.0);

      REAL8 S = pow (C, 2.0);
      REAL8 T = 1.0 / tanh (I);

      do {(&z)->re=(0.5 * sin (2 * R) * S / D); (&z)->im=(T / D);} while(0);
    }

  return z;
}

COMPLEX16
XLALCOMPLEX16Sec (COMPLEX16 a)
{
  COMPLEX16 z = XLALCOMPLEX16Cos (a);
  return XLALCOMPLEX16Inverse (z);
}

COMPLEX16
XLALCOMPLEX16Csc (COMPLEX16 a)
{
  COMPLEX16 z = XLALCOMPLEX16Sin (a);
  return XLALCOMPLEX16Inverse(z);
}


COMPLEX16
XLALCOMPLEX16Cot (COMPLEX16 a)
{
  COMPLEX16 z = XLALCOMPLEX16Tan (a);
  return XLALCOMPLEX16Inverse (z);
}





COMPLEX16
XLALCOMPLEX16Arcsin (COMPLEX16 a)
{
  REAL8 R = ((a).re), I = ((a).im);
  COMPLEX16 z;

  if (I == 0)
    {
      z = XLALCOMPLEX16ArcsinReal (R);
    }
  else
    {
      REAL8 x = fabs (R), y = fabs (I);
      REAL8 r = hypot (x + 1, y), s = hypot (x - 1, y);
      REAL8 A = 0.5 * (r + s);
      REAL8 B = x / A;
      REAL8 y2 = y * y;

      REAL8 real, imag;

      const REAL8 A_crossover = 1.5, B_crossover = 0.6417;

      if (B <= B_crossover)
        {
          real = asin (B);
        }
      else
        {
          if (x <= 1)
            {
              REAL8 D = 0.5 * (A + x) * (y2 / (r + x + 1) + (s + (1 - x)));
              real = atan (x / sqrt (D));
            }
          else
            {
              REAL8 Apx = A + x;
              REAL8 D = 0.5 * (Apx / (r + x + 1) + Apx / (s + (x - 1)));
              real = atan (x / (y * sqrt (D)));
            }
        }

      if (A <= A_crossover)
        {
          REAL8 Am1;

          if (x < 1)
            {
              Am1 = 0.5 * (y2 / (r + (x + 1)) + y2 / (s + (1 - x)));
            }
          else
            {
              Am1 = 0.5 * (y2 / (r + (x + 1)) + (s + (x - 1)));
            }

          imag = log1p (Am1 + sqrt (Am1 * (A + 1)));
        }
      else
        {
          imag = log (A + sqrt (A * A - 1));
        }

      do {(&z)->re=((R >= 0) ? real : -real); (&z)->im=((I >= 0) ? imag : -imag);} while(0);
    }

  return z;
}

COMPLEX16
XLALCOMPLEX16ArcsinReal (REAL8 a)
{
  COMPLEX16 z;

  if (fabs (a) <= 1.0)
    {
      do {(&z)->re=(asin (a)); (&z)->im=(0.0);} while(0);
    }
  else
    {
      if (a < 0.0)
        {
          do {(&z)->re=(-1.57079632679489661923132169163975144); (&z)->im=(acosh (-a));} while(0);
        }
      else
        {
          do {(&z)->re=(1.57079632679489661923132169163975144); (&z)->im=(-acosh (a));} while(0);
        }
    }

  return z;
}

COMPLEX16
XLALCOMPLEX16Arccos (COMPLEX16 a)
{
  REAL8 R = ((a).re), I = ((a).im);
  COMPLEX16 z;

  if (I == 0)
    {
      z = XLALCOMPLEX16ArccosReal (R);
    }
  else
    {
      REAL8 x = fabs (R), y = fabs (I);
      REAL8 r = hypot (x + 1, y), s = hypot (x - 1, y);
      REAL8 A = 0.5 * (r + s);
      REAL8 B = x / A;
      REAL8 y2 = y * y;

      REAL8 real, imag;

      const REAL8 A_crossover = 1.5, B_crossover = 0.6417;

      if (B <= B_crossover)
        {
          real = acos (B);
        }
      else
        {
          if (x <= 1)
            {
              REAL8 D = 0.5 * (A + x) * (y2 / (r + x + 1) + (s + (1 - x)));
              real = atan (sqrt (D) / x);
            }
          else
            {
              REAL8 Apx = A + x;
              REAL8 D = 0.5 * (Apx / (r + x + 1) + Apx / (s + (x - 1)));
              real = atan ((y * sqrt (D)) / x);
            }
        }

      if (A <= A_crossover)
        {
          REAL8 Am1;

          if (x < 1)
            {
              Am1 = 0.5 * (y2 / (r + (x + 1)) + y2 / (s + (1 - x)));
            }
          else
            {
              Am1 = 0.5 * (y2 / (r + (x + 1)) + (s + (x - 1)));
            }

          imag = log1p (Am1 + sqrt (Am1 * (A + 1)));
        }
      else
        {
          imag = log (A + sqrt (A * A - 1));
        }

      do {(&z)->re=((R >= 0) ? real : 3.14159265358979323846264338327950288 - real); (&z)->im=((I >= 0) ? -imag : imag);} while(0);
    }

  return z;
}

COMPLEX16
XLALCOMPLEX16ArccosReal (REAL8 a)
{
  COMPLEX16 z;

  if (fabs (a) <= 1.0)
    {
      do {(&z)->re=(acos (a)); (&z)->im=(0);} while(0);
    }
  else
    {
      if (a < 0.0)
        {
          do {(&z)->re=(3.14159265358979323846264338327950288); (&z)->im=(-acosh (-a));} while(0);
        }
      else
        {
          do {(&z)->re=(0); (&z)->im=(acosh (a));} while(0);
        }
    }

  return z;
}

COMPLEX16
XLALCOMPLEX16Arctan (COMPLEX16 a)
{
  REAL8 R = ((a).re), I = ((a).im);
  COMPLEX16 z;

  if (I == 0)
    {
      do {(&z)->re=(atan (R)); (&z)->im=(0);} while(0);
    }
  else
    {




      REAL8 r = hypot (R, I);

      REAL8 imag;

      REAL8 u = 2 * I / (1 + r * r);




      if (fabs (u) < 0.1)
        {
          imag = 0.25 * (log1p (u) - log1p (-u));
        }
      else
        {
          REAL8 A = hypot (R, I + 1);
          REAL8 B = hypot (R, I - 1);
          imag = 0.5 * log (A / B);
        }

      if (R == 0)
        {
          if (I > 1)
            {
              do {(&z)->re=(1.57079632679489661923132169163975144); (&z)->im=(imag);} while(0);
            }
          else if (I < -1)
            {
              do {(&z)->re=(-1.57079632679489661923132169163975144); (&z)->im=(imag);} while(0);
            }
          else
            {
              do {(&z)->re=(0); (&z)->im=(imag);} while(0);
            };
        }
      else
        {
          do {(&z)->re=(0.5 * atan2 (2 * R, ((1 + r) * (1 - r)))); (&z)->im=(imag);} while(0);
        }
    }

  return z;
}

COMPLEX16
XLALCOMPLEX16Arcsec (COMPLEX16 a)
{
  COMPLEX16 z = XLALCOMPLEX16Inverse (a);
  return XLALCOMPLEX16Arccos (z);
}

COMPLEX16
XLALCOMPLEX16ArcsecReal (REAL8 a)
{
  COMPLEX16 z;

  if (a <= -1.0 || a >= 1.0)
    {
      do {(&z)->re=(acos (1 / a)); (&z)->im=(0.0);} while(0);
    }
  else
    {
      if (a >= 0.0)
        {
          do {(&z)->re=(0); (&z)->im=(acosh (1 / a));} while(0);
        }
      else
        {
          do {(&z)->re=(3.14159265358979323846264338327950288); (&z)->im=(-acosh (-1 / a));} while(0);
        }
    }

  return z;
}

COMPLEX16
XLALCOMPLEX16Arccsc (COMPLEX16 a)
{
  COMPLEX16 z = XLALCOMPLEX16Inverse (a);
  return XLALCOMPLEX16Arcsin (z);
}

COMPLEX16
XLALCOMPLEX16ArccscReal (REAL8 a)
{
  COMPLEX16 z;

  if (a <= -1.0 || a >= 1.0)
    {
      do {(&z)->re=(asin (1 / a)); (&z)->im=(0.0);} while(0);
    }
  else
    {
      if (a >= 0.0)
        {
          do {(&z)->re=(1.57079632679489661923132169163975144); (&z)->im=(-acosh (1 / a));} while(0);
        }
      else
        {
          do {(&z)->re=(-1.57079632679489661923132169163975144); (&z)->im=(acosh (-1 / a));} while(0);
        }
    }

  return z;
}

COMPLEX16
XLALCOMPLEX16Arccot (COMPLEX16 a)
{
  COMPLEX16 z;

  if (((a).re) == 0.0 && ((a).im) == 0.0)
    {
      do {(&z)->re=(1.57079632679489661923132169163975144); (&z)->im=(0);} while(0);
    }
  else
    {
      z = XLALCOMPLEX16Inverse (a);
      z = XLALCOMPLEX16Arctan (z);
    }

  return z;
}





COMPLEX16
XLALCOMPLEX16Sinh (COMPLEX16 a)
{
  REAL8 R = ((a).re), I = ((a).im);

  COMPLEX16 z;
  do {(&z)->re=(sinh (R) * cos (I)); (&z)->im=(cosh (R) * sin (I));} while(0);
  return z;
}

COMPLEX16
XLALCOMPLEX16Cosh (COMPLEX16 a)
{
  REAL8 R = ((a).re), I = ((a).im);

  COMPLEX16 z;
  do {(&z)->re=(cosh (R) * cos (I)); (&z)->im=(sinh (R) * sin (I));} while(0);
  return z;
}

COMPLEX16
XLALCOMPLEX16Tanh (COMPLEX16 a)
{
  REAL8 R = ((a).re), I = ((a).im);

  COMPLEX16 z;

  if (fabs(R) < 1.0)
    {
      REAL8 D = pow (cos (I), 2.0) + pow (sinh (R), 2.0);

      do {(&z)->re=(sinh (R) * cosh (R) / D); (&z)->im=(0.5 * sin (2 * I) / D);} while(0);
    }
  else
    {
      REAL8 D = pow (cos (I), 2.0) + pow (sinh (R), 2.0);
      REAL8 F = 1 + pow (cos (I) / sinh (R), 2.0);

      do {(&z)->re=(1.0 / (tanh (R) * F)); (&z)->im=(0.5 * sin (2 * I) / D);} while(0);
    }

  return z;
}

COMPLEX16
XLALCOMPLEX16Sech (COMPLEX16 a)
{
  COMPLEX16 z = XLALCOMPLEX16Cosh (a);
  return XLALCOMPLEX16Inverse (z);
}

COMPLEX16
XLALCOMPLEX16Csch (COMPLEX16 a)
{
  COMPLEX16 z = XLALCOMPLEX16Sinh (a);
  return XLALCOMPLEX16Inverse (z);
}

COMPLEX16
XLALCOMPLEX16Coth (COMPLEX16 a)
{
  COMPLEX16 z = XLALCOMPLEX16Tanh (a);
  return XLALCOMPLEX16Inverse (z);
}





COMPLEX16
XLALCOMPLEX16Arcsinh (COMPLEX16 a)
{
  COMPLEX16 z = XLALCOMPLEX16MulImag(a, 1.0);
  z = XLALCOMPLEX16Arcsin (z);
  z = XLALCOMPLEX16MulImag (z, -1.0);
  return z;
}

COMPLEX16
XLALCOMPLEX16Arccosh (COMPLEX16 a)
{
  COMPLEX16 z = XLALCOMPLEX16Arccos (a);
  z = XLALCOMPLEX16MulImag (z, ((z).im) > 0 ? -1.0 : 1.0);
  return z;
}

COMPLEX16
XLALCOMPLEX16ArccoshReal (REAL8 a)
{
  COMPLEX16 z;

  if (a >= 1)
    {
      do {(&z)->re=(acosh (a)); (&z)->im=(0);} while(0);
    }
  else
    {
      if (a >= -1.0)
        {
          do {(&z)->re=(0); (&z)->im=(acos (a));} while(0);
        }
      else
        {
          do {(&z)->re=(acosh (-a)); (&z)->im=(3.14159265358979323846264338327950288);} while(0);
        }
    }

  return z;
}

COMPLEX16
XLALCOMPLEX16Arctanh (COMPLEX16 a)
{
  if (((a).im) == 0.0)
    {
      return XLALCOMPLEX16ArctanhReal (((a).re));
    }
  else
    {
      COMPLEX16 z = XLALCOMPLEX16MulImag(a, 1.0);
      z = XLALCOMPLEX16Arctan (z);
      z = XLALCOMPLEX16MulImag (z, -1.0);
      return z;
    }
}

COMPLEX16
XLALCOMPLEX16ArctanhReal (REAL8 a)
{
  COMPLEX16 z;

  if (a > -1.0 && a < 1.0)
    {
      do {(&z)->re=(atanh (a)); (&z)->im=(0);} while(0);
    }
  else
    {
      do {(&z)->re=(atanh (1 / a)); (&z)->im=((a < 0) ? 1.57079632679489661923132169163975144 : -1.57079632679489661923132169163975144);} while(0);
    }

  return z;
}

COMPLEX16
XLALCOMPLEX16Arcsech (COMPLEX16 a)
{
  COMPLEX16 t = XLALCOMPLEX16Inverse (a);
  return XLALCOMPLEX16Arccosh (t);
}

COMPLEX16
XLALCOMPLEX16Arccsch (COMPLEX16 a)
{
  COMPLEX16 t = XLALCOMPLEX16Inverse (a);
  return XLALCOMPLEX16Arcsinh (t);
}

COMPLEX16
XLALCOMPLEX16Arccoth (COMPLEX16 a)
{
  COMPLEX16 t = XLALCOMPLEX16Inverse (a);
  return XLALCOMPLEX16Arctanh (t);
}


COMPLEX8
XLALCOMPLEX8Rect (REAL4 x, REAL4 y)
{
  COMPLEX8 z;
  do {(&z)->re=(x); (&z)->im=(y);} while(0);
  return z;
}

COMPLEX8
XLALCOMPLEX8Polar (REAL4 r, REAL4 theta)
{
  COMPLEX8 z;
  do {(&z)->re=(r * cos (theta)); (&z)->im=(r * sin (theta));} while(0);
  return z;
}





REAL4
XLALCOMPLEX8Arg (COMPLEX8 z)
{
  REAL4 x = ((z).re);
  REAL4 y = ((z).im);

  if (x == 0.0 && y == 0.0)
    {
      return 0;
    }

  return atan2 (y, x);
}

REAL4
XLALCOMPLEX8Abs (COMPLEX8 z)
{
  return hypot (((z).re), ((z).im));
}

REAL4
XLALCOMPLEX8Abs2 (COMPLEX8 z)
{
  REAL4 x = ((z).re);
  REAL4 y = ((z).im);

  return (x * x + y * y);
}

REAL4
XLALCOMPLEX8LogAbs (COMPLEX8 z)
{
  REAL4 xabs = fabs (((z).re));
  REAL4 yabs = fabs (((z).im));
  REAL4 max, u;

  if (xabs >= yabs)
    {
      max = xabs;
      u = yabs / xabs;
    }
  else
    {
      max = yabs;
      u = xabs / yabs;
    }



  return log (max) + 0.5 * log1p (u * u);
}






COMPLEX8
XLALCOMPLEX8Add (COMPLEX8 a, COMPLEX8 b)
{
  REAL4 ar = ((a).re), ai = ((a).im);
  REAL4 br = ((b).re), bi = ((b).im);

  COMPLEX8 z;
  do {(&z)->re=(ar + br); (&z)->im=(ai + bi);} while(0);
  return z;
}

COMPLEX8
XLALCOMPLEX8AddReal (COMPLEX8 a, REAL4 x)
{
  COMPLEX8 z;
  do {(&z)->re=(((a).re) + x); (&z)->im=(((a).im));} while(0);
  return z;
}

COMPLEX8
XLALCOMPLEX8AddImag (COMPLEX8 a, REAL4 y)
{
  COMPLEX8 z;
  do {(&z)->re=(((a).re)); (&z)->im=(((a).im) + y);} while(0);
  return z;
}


COMPLEX8
XLALCOMPLEX8Sub (COMPLEX8 a, COMPLEX8 b)
{
  REAL4 ar = ((a).re), ai = ((a).im);
  REAL4 br = ((b).re), bi = ((b).im);

  COMPLEX8 z;
  do {(&z)->re=(ar - br); (&z)->im=(ai - bi);} while(0);
  return z;
}

COMPLEX8
XLALCOMPLEX8SubReal (COMPLEX8 a, REAL4 x)
{
  COMPLEX8 z;
  do {(&z)->re=(((a).re) - x); (&z)->im=(((a).im));} while(0);
  return z;
}

COMPLEX8
XLALCOMPLEX8SubImag (COMPLEX8 a, REAL4 y)
{
  COMPLEX8 z;
  do {(&z)->re=(((a).re)); (&z)->im=(((a).im) - y);} while(0);
  return z;
}

COMPLEX8
XLALCOMPLEX8Mul (COMPLEX8 a, COMPLEX8 b)
{
  REAL4 ar = ((a).re), ai = ((a).im);
  REAL4 br = ((b).re), bi = ((b).im);

  COMPLEX8 z;
  do {(&z)->re=(ar * br - ai * bi); (&z)->im=(ar * bi + ai * br);} while(0);
  return z;
}

COMPLEX8
XLALCOMPLEX8MulReal (COMPLEX8 a, REAL4 x)
{
  COMPLEX8 z;
  do {(&z)->re=(x * ((a).re)); (&z)->im=(x * ((a).im));} while(0);
  return z;
}

COMPLEX8
XLALCOMPLEX8MulImag (COMPLEX8 a, REAL4 y)
{
  COMPLEX8 z;
  do {(&z)->re=(-y * ((a).im)); (&z)->im=(y * ((a).re));} while(0);
  return z;
}

COMPLEX8
XLALCOMPLEX8Div (COMPLEX8 a, COMPLEX8 b)
{
  REAL4 ar = ((a).re), ai = ((a).im);
  REAL4 br = ((b).re), bi = ((b).im);

  REAL4 s = 1.0 / XLALCOMPLEX8Abs (b);

  REAL4 sbr = s * br;
  REAL4 sbi = s * bi;

  REAL4 zr = (ar * sbr + ai * sbi) * s;
  REAL4 zi = (ai * sbr - ar * sbi) * s;

  COMPLEX8 z;
  do {(&z)->re=(zr); (&z)->im=(zi);} while(0);
  return z;
}

COMPLEX8
XLALCOMPLEX8DivReal (COMPLEX8 a, REAL4 x)
{
  COMPLEX8 z;
  do {(&z)->re=(((a).re) / x); (&z)->im=(((a).im) / x);} while(0);
  return z;
}

COMPLEX8
XLALCOMPLEX8DivImag (COMPLEX8 a, REAL4 y)
{
  COMPLEX8 z;
  do {(&z)->re=(((a).im) / y); (&z)->im=(- ((a).re) / y);} while(0);
  return z;
}

COMPLEX8
XLALCOMPLEX8Conjugate (COMPLEX8 a)
{
  COMPLEX8 z;
  do {(&z)->re=(((a).re)); (&z)->im=(-((a).im));} while(0);
  return z;
}

COMPLEX8
XLALCOMPLEX8Negative (COMPLEX8 a)
{
  COMPLEX8 z;
  do {(&z)->re=(-((a).re)); (&z)->im=(-((a).im));} while(0);
  return z;
}

COMPLEX8
XLALCOMPLEX8Inverse (COMPLEX8 a)
{
  REAL4 s = 1.0 / XLALCOMPLEX8Abs (a);

  COMPLEX8 z;
  do {(&z)->re=((((a).re) * s) * s); (&z)->im=(-(((a).im) * s) * s);} while(0);
  return z;
}





COMPLEX8
XLALCOMPLEX8Sqrt (COMPLEX8 a)
{
  COMPLEX8 z;

  if (((a).re) == 0.0 && ((a).im) == 0.0)
    {
      do {(&z)->re=(0); (&z)->im=(0);} while(0);
    }
  else
    {
      REAL4 x = fabs (((a).re));
      REAL4 y = fabs (((a).im));
      REAL4 w;

      if (x >= y)
        {
          REAL4 t = y / x;
          w = sqrt (x) * sqrt (0.5 * (1.0 + sqrt (1.0 + t * t)));
        }
      else
        {
          REAL4 t = x / y;
          w = sqrt (y) * sqrt (0.5 * (t + sqrt (1.0 + t * t)));
        }

      if (((a).re) >= 0.0)
        {
          REAL4 ai = ((a).im);
          do {(&z)->re=(w); (&z)->im=(ai / (2.0 * w));} while(0);
        }
      else
        {
          REAL4 ai = ((a).im);
          REAL4 vi = (ai >= 0) ? w : -w;
          do {(&z)->re=(ai / (2.0 * vi)); (&z)->im=(vi);} while(0);
        }
    }

  return z;
}

COMPLEX8
XLALCOMPLEX8SqrtReal (REAL4 x)
{
  COMPLEX8 z;

  if (x >= 0)
    {
      do {(&z)->re=(sqrt (x)); (&z)->im=(0.0);} while(0);
    }
  else
    {
      do {(&z)->re=(0.0); (&z)->im=(sqrt (-x));} while(0);
    }

  return z;
}

COMPLEX8
XLALCOMPLEX8Exp (COMPLEX8 a)
{
  REAL4 rho = exp (((a).re));
  REAL4 theta = ((a).im);

  COMPLEX8 z;
  do {(&z)->re=(rho * cos (theta)); (&z)->im=(rho * sin (theta));} while(0);
  return z;
}

COMPLEX8
XLALCOMPLEX8Pow (COMPLEX8 a, COMPLEX8 b)
{
  COMPLEX8 z;

  if (((a).re) == 0 && ((a).im) == 0.0)
    {
      do {(&z)->re=(0.0); (&z)->im=(0.0);} while(0);
    }
  else
    {
      REAL4 logr = XLALCOMPLEX8LogAbs (a);
      REAL4 theta = XLALCOMPLEX8Arg (a);

      REAL4 br = ((b).re), bi = ((b).im);

      REAL4 rho = exp (logr * br - bi * theta);
      REAL4 beta = theta * br + bi * logr;

      do {(&z)->re=(rho * cos (beta)); (&z)->im=(rho * sin (beta));} while(0);
    }

  return z;
}

COMPLEX8
XLALCOMPLEX8PowReal (COMPLEX8 a, REAL4 b)
{
  COMPLEX8 z;

  if (((a).re) == 0 && ((a).im) == 0)
    {
      do {(&z)->re=(0); (&z)->im=(0);} while(0);
    }
  else
    {
      REAL4 logr = XLALCOMPLEX8LogAbs (a);
      REAL4 theta = XLALCOMPLEX8Arg (a);
      REAL4 rho = exp (logr * b);
      REAL4 beta = theta * b;
      do {(&z)->re=(rho * cos (beta)); (&z)->im=(rho * sin (beta));} while(0);
    }

  return z;
}

COMPLEX8
XLALCOMPLEX8Log (COMPLEX8 a)
{
  REAL4 logr = XLALCOMPLEX8LogAbs (a);
  REAL4 theta = XLALCOMPLEX8Arg (a);

  COMPLEX8 z;
  do {(&z)->re=(logr); (&z)->im=(theta);} while(0);
  return z;
}

COMPLEX8
XLALCOMPLEX8Log10 (COMPLEX8 a)
{
  return XLALCOMPLEX8MulReal (XLALCOMPLEX8Log (a), 1 / log (10.));
}

COMPLEX8
XLALCOMPLEX8LogB (COMPLEX8 a, COMPLEX8 b)
{
  return XLALCOMPLEX8Div (XLALCOMPLEX8Log (a), XLALCOMPLEX8Log (b));
}





COMPLEX8
XLALCOMPLEX8Sin (COMPLEX8 a)
{
  REAL4 R = ((a).re), I = ((a).im);

  COMPLEX8 z;

  if (I == 0.0)
    {


      do {(&z)->re=(sin (R)); (&z)->im=(0.0);} while(0);
    }
  else
    {
      do {(&z)->re=(sin (R) * cosh (I)); (&z)->im=(cos (R) * sinh (I));} while(0);
    }

  return z;
}

COMPLEX8
XLALCOMPLEX8Cos (COMPLEX8 a)
{
  REAL4 R = ((a).re), I = ((a).im);

  COMPLEX8 z;

  if (I == 0.0)
    {


      do {(&z)->re=(cos (R)); (&z)->im=(0.0);} while(0);
    }
  else
    {
      do {(&z)->re=(cos (R) * cosh (I)); (&z)->im=(sin (R) * sinh (-I));} while(0);
    }

  return z;
}

COMPLEX8
XLALCOMPLEX8Tan (COMPLEX8 a)
{
  REAL4 R = ((a).re), I = ((a).im);

  COMPLEX8 z;

  if (fabs (I) < 1)
    {
      REAL4 D = pow (cos (R), 2.0) + pow (sinh (I), 2.0);

      do {(&z)->re=(0.5 * sin (2 * R) / D); (&z)->im=(0.5 * sinh (2 * I) / D);} while(0);
    }
  else
    {
      REAL4 u = exp (-I);
      REAL4 C = 2 * u / (1 - pow (u, 2.0));
      REAL4 D = 1 + pow (cos (R), 2.0) * pow (C, 2.0);

      REAL4 S = pow (C, 2.0);
      REAL4 T = 1.0 / tanh (I);

      do {(&z)->re=(0.5 * sin (2 * R) * S / D); (&z)->im=(T / D);} while(0);
    }

  return z;
}

COMPLEX8
XLALCOMPLEX8Sec (COMPLEX8 a)
{
  COMPLEX8 z = XLALCOMPLEX8Cos (a);
  return XLALCOMPLEX8Inverse (z);
}

COMPLEX8
XLALCOMPLEX8Csc (COMPLEX8 a)
{
  COMPLEX8 z = XLALCOMPLEX8Sin (a);
  return XLALCOMPLEX8Inverse(z);
}


COMPLEX8
XLALCOMPLEX8Cot (COMPLEX8 a)
{
  COMPLEX8 z = XLALCOMPLEX8Tan (a);
  return XLALCOMPLEX8Inverse (z);
}





COMPLEX8
XLALCOMPLEX8Arcsin (COMPLEX8 a)
{
  REAL4 R = ((a).re), I = ((a).im);
  COMPLEX8 z;

  if (I == 0)
    {
      z = XLALCOMPLEX8ArcsinReal (R);
    }
  else
    {
      REAL4 x = fabs (R), y = fabs (I);
      REAL4 r = hypot (x + 1, y), s = hypot (x - 1, y);
      REAL4 A = 0.5 * (r + s);
      REAL4 B = x / A;
      REAL4 y2 = y * y;

      REAL4 real, imag;

      const REAL4 A_crossover = 1.5, B_crossover = 0.6417;

      if (B <= B_crossover)
        {
          real = asin (B);
        }
      else
        {
          if (x <= 1)
            {
              REAL4 D = 0.5 * (A + x) * (y2 / (r + x + 1) + (s + (1 - x)));
              real = atan (x / sqrt (D));
            }
          else
            {
              REAL4 Apx = A + x;
              REAL4 D = 0.5 * (Apx / (r + x + 1) + Apx / (s + (x - 1)));
              real = atan (x / (y * sqrt (D)));
            }
        }

      if (A <= A_crossover)
        {
          REAL4 Am1;

          if (x < 1)
            {
              Am1 = 0.5 * (y2 / (r + (x + 1)) + y2 / (s + (1 - x)));
            }
          else
            {
              Am1 = 0.5 * (y2 / (r + (x + 1)) + (s + (x - 1)));
            }

          imag = log1p (Am1 + sqrt (Am1 * (A + 1)));
        }
      else
        {
          imag = log (A + sqrt (A * A - 1));
        }

      do {(&z)->re=((R >= 0) ? real : -real); (&z)->im=((I >= 0) ? imag : -imag);} while(0);
    }

  return z;
}

COMPLEX8
XLALCOMPLEX8ArcsinReal (REAL4 a)
{
  COMPLEX8 z;

  if (fabs (a) <= 1.0)
    {
      do {(&z)->re=(asin (a)); (&z)->im=(0.0);} while(0);
    }
  else
    {
      if (a < 0.0)
        {
          do {(&z)->re=(-1.57079632679489661923132169163975144); (&z)->im=(acosh (-a));} while(0);
        }
      else
        {
          do {(&z)->re=(1.57079632679489661923132169163975144); (&z)->im=(-acosh (a));} while(0);
        }
    }

  return z;
}

COMPLEX8
XLALCOMPLEX8Arccos (COMPLEX8 a)
{
  REAL4 R = ((a).re), I = ((a).im);
  COMPLEX8 z;

  if (I == 0)
    {
      z = XLALCOMPLEX8ArccosReal (R);
    }
  else
    {
      REAL4 x = fabs (R), y = fabs (I);
      REAL4 r = hypot (x + 1, y), s = hypot (x - 1, y);
      REAL4 A = 0.5 * (r + s);
      REAL4 B = x / A;
      REAL4 y2 = y * y;

      REAL4 real, imag;

      const REAL4 A_crossover = 1.5, B_crossover = 0.6417;

      if (B <= B_crossover)
        {
          real = acos (B);
        }
      else
        {
          if (x <= 1)
            {
              REAL4 D = 0.5 * (A + x) * (y2 / (r + x + 1) + (s + (1 - x)));
              real = atan (sqrt (D) / x);
            }
          else
            {
              REAL4 Apx = A + x;
              REAL4 D = 0.5 * (Apx / (r + x + 1) + Apx / (s + (x - 1)));
              real = atan ((y * sqrt (D)) / x);
            }
        }

      if (A <= A_crossover)
        {
          REAL4 Am1;

          if (x < 1)
            {
              Am1 = 0.5 * (y2 / (r + (x + 1)) + y2 / (s + (1 - x)));
            }
          else
            {
              Am1 = 0.5 * (y2 / (r + (x + 1)) + (s + (x - 1)));
            }

          imag = log1p (Am1 + sqrt (Am1 * (A + 1)));
        }
      else
        {
          imag = log (A + sqrt (A * A - 1));
        }

      do {(&z)->re=((R >= 0) ? real : 3.14159265358979323846264338327950288 - real); (&z)->im=((I >= 0) ? -imag : imag);} while(0);
    }

  return z;
}

COMPLEX8
XLALCOMPLEX8ArccosReal (REAL4 a)
{
  COMPLEX8 z;

  if (fabs (a) <= 1.0)
    {
      do {(&z)->re=(acos (a)); (&z)->im=(0);} while(0);
    }
  else
    {
      if (a < 0.0)
        {
          do {(&z)->re=(3.14159265358979323846264338327950288); (&z)->im=(-acosh (-a));} while(0);
        }
      else
        {
          do {(&z)->re=(0); (&z)->im=(acosh (a));} while(0);
        }
    }

  return z;
}

COMPLEX8
XLALCOMPLEX8Arctan (COMPLEX8 a)
{
  REAL4 R = ((a).re), I = ((a).im);
  COMPLEX8 z;

  if (I == 0)
    {
      do {(&z)->re=(atan (R)); (&z)->im=(0);} while(0);
    }
  else
    {




      REAL4 r = hypot (R, I);

      REAL4 imag;

      REAL4 u = 2 * I / (1 + r * r);




      if (fabs (u) < 0.1)
        {
          imag = 0.25 * (log1p (u) - log1p (-u));
        }
      else
        {
          REAL4 A = hypot (R, I + 1);
          REAL4 B = hypot (R, I - 1);
          imag = 0.5 * log (A / B);
        }

      if (R == 0)
        {
          if (I > 1)
            {
              do {(&z)->re=(1.57079632679489661923132169163975144); (&z)->im=(imag);} while(0);
            }
          else if (I < -1)
            {
              do {(&z)->re=(-1.57079632679489661923132169163975144); (&z)->im=(imag);} while(0);
            }
          else
            {
              do {(&z)->re=(0); (&z)->im=(imag);} while(0);
            };
        }
      else
        {
          do {(&z)->re=(0.5 * atan2 (2 * R, ((1 + r) * (1 - r)))); (&z)->im=(imag);} while(0);
        }
    }

  return z;
}

COMPLEX8
XLALCOMPLEX8Arcsec (COMPLEX8 a)
{
  COMPLEX8 z = XLALCOMPLEX8Inverse (a);
  return XLALCOMPLEX8Arccos (z);
}

COMPLEX8
XLALCOMPLEX8ArcsecReal (REAL4 a)
{
  COMPLEX8 z;

  if (a <= -1.0 || a >= 1.0)
    {
      do {(&z)->re=(acos (1 / a)); (&z)->im=(0.0);} while(0);
    }
  else
    {
      if (a >= 0.0)
        {
          do {(&z)->re=(0); (&z)->im=(acosh (1 / a));} while(0);
        }
      else
        {
          do {(&z)->re=(3.14159265358979323846264338327950288); (&z)->im=(-acosh (-1 / a));} while(0);
        }
    }

  return z;
}

COMPLEX8
XLALCOMPLEX8Arccsc (COMPLEX8 a)
{
  COMPLEX8 z = XLALCOMPLEX8Inverse (a);
  return XLALCOMPLEX8Arcsin (z);
}

COMPLEX8
XLALCOMPLEX8ArccscReal (REAL4 a)
{
  COMPLEX8 z;

  if (a <= -1.0 || a >= 1.0)
    {
      do {(&z)->re=(asin (1 / a)); (&z)->im=(0.0);} while(0);
    }
  else
    {
      if (a >= 0.0)
        {
          do {(&z)->re=(1.57079632679489661923132169163975144); (&z)->im=(-acosh (1 / a));} while(0);
        }
      else
        {
          do {(&z)->re=(-1.57079632679489661923132169163975144); (&z)->im=(acosh (-1 / a));} while(0);
        }
    }

  return z;
}

COMPLEX8
XLALCOMPLEX8Arccot (COMPLEX8 a)
{
  COMPLEX8 z;

  if (((a).re) == 0.0 && ((a).im) == 0.0)
    {
      do {(&z)->re=(1.57079632679489661923132169163975144); (&z)->im=(0);} while(0);
    }
  else
    {
      z = XLALCOMPLEX8Inverse (a);
      z = XLALCOMPLEX8Arctan (z);
    }

  return z;
}





COMPLEX8
XLALCOMPLEX8Sinh (COMPLEX8 a)
{
  REAL4 R = ((a).re), I = ((a).im);

  COMPLEX8 z;
  do {(&z)->re=(sinh (R) * cos (I)); (&z)->im=(cosh (R) * sin (I));} while(0);
  return z;
}

COMPLEX8
XLALCOMPLEX8Cosh (COMPLEX8 a)
{
  REAL4 R = ((a).re), I = ((a).im);

  COMPLEX8 z;
  do {(&z)->re=(cosh (R) * cos (I)); (&z)->im=(sinh (R) * sin (I));} while(0);
  return z;
}

COMPLEX8
XLALCOMPLEX8Tanh (COMPLEX8 a)
{
  REAL4 R = ((a).re), I = ((a).im);

  COMPLEX8 z;

  if (fabs(R) < 1.0)
    {
      REAL4 D = pow (cos (I), 2.0) + pow (sinh (R), 2.0);

      do {(&z)->re=(sinh (R) * cosh (R) / D); (&z)->im=(0.5 * sin (2 * I) / D);} while(0);
    }
  else
    {
      REAL4 D = pow (cos (I), 2.0) + pow (sinh (R), 2.0);
      REAL4 F = 1 + pow (cos (I) / sinh (R), 2.0);

      do {(&z)->re=(1.0 / (tanh (R) * F)); (&z)->im=(0.5 * sin (2 * I) / D);} while(0);
    }

  return z;
}

COMPLEX8
XLALCOMPLEX8Sech (COMPLEX8 a)
{
  COMPLEX8 z = XLALCOMPLEX8Cosh (a);
  return XLALCOMPLEX8Inverse (z);
}

COMPLEX8
XLALCOMPLEX8Csch (COMPLEX8 a)
{
  COMPLEX8 z = XLALCOMPLEX8Sinh (a);
  return XLALCOMPLEX8Inverse (z);
}

COMPLEX8
XLALCOMPLEX8Coth (COMPLEX8 a)
{
  COMPLEX8 z = XLALCOMPLEX8Tanh (a);
  return XLALCOMPLEX8Inverse (z);
}





COMPLEX8
XLALCOMPLEX8Arcsinh (COMPLEX8 a)
{
  COMPLEX8 z = XLALCOMPLEX8MulImag(a, 1.0);
  z = XLALCOMPLEX8Arcsin (z);
  z = XLALCOMPLEX8MulImag (z, -1.0);
  return z;
}

COMPLEX8
XLALCOMPLEX8Arccosh (COMPLEX8 a)
{
  COMPLEX8 z = XLALCOMPLEX8Arccos (a);
  z = XLALCOMPLEX8MulImag (z, ((z).im) > 0 ? -1.0 : 1.0);
  return z;
}

COMPLEX8
XLALCOMPLEX8ArccoshReal (REAL4 a)
{
  COMPLEX8 z;

  if (a >= 1)
    {
      do {(&z)->re=(acosh (a)); (&z)->im=(0);} while(0);
    }
  else
    {
      if (a >= -1.0)
        {
          do {(&z)->re=(0); (&z)->im=(acos (a));} while(0);
        }
      else
        {
          do {(&z)->re=(acosh (-a)); (&z)->im=(3.14159265358979323846264338327950288);} while(0);
        }
    }

  return z;
}

COMPLEX8
XLALCOMPLEX8Arctanh (COMPLEX8 a)
{
  if (((a).im) == 0.0)
    {
      return XLALCOMPLEX8ArctanhReal (((a).re));
    }
  else
    {
      COMPLEX8 z = XLALCOMPLEX8MulImag(a, 1.0);
      z = XLALCOMPLEX8Arctan (z);
      z = XLALCOMPLEX8MulImag (z, -1.0);
      return z;
    }
}

COMPLEX8
XLALCOMPLEX8ArctanhReal (REAL4 a)
{
  COMPLEX8 z;

  if (a > -1.0 && a < 1.0)
    {
      do {(&z)->re=(atanh (a)); (&z)->im=(0);} while(0);
    }
  else
    {
      do {(&z)->re=(atanh (1 / a)); (&z)->im=((a < 0) ? 1.57079632679489661923132169163975144 : -1.57079632679489661923132169163975144);} while(0);
    }

  return z;
}

COMPLEX8
XLALCOMPLEX8Arcsech (COMPLEX8 a)
{
  COMPLEX8 t = XLALCOMPLEX8Inverse (a);
  return XLALCOMPLEX8Arccosh (t);
}

COMPLEX8
XLALCOMPLEX8Arccsch (COMPLEX8 a)
{
  COMPLEX8 t = XLALCOMPLEX8Inverse (a);
  return XLALCOMPLEX8Arcsinh (t);
}

COMPLEX8
XLALCOMPLEX8Arccoth (COMPLEX8 a)
{
  COMPLEX8 t = XLALCOMPLEX8Inverse (a);
  return XLALCOMPLEX8Arctanh (t);
}

/** \endcond */
