/** the way of trimming x to the interval [0..1) for the sin_cos_LUT functions
    give significant differences in speed, so we provide various ways here.
    We also record the way we are using for logging */

#if EAH_SINCOS_ROUND == EAH_SINCOS_ROUND_PLUS2
/* this only makes sense for the linear sin/cos approximation */
#ifdef _MSC_VER /* no C99 */
#define SINCOS_TRIM_X(y,x) \
  { \
    __asm FLD     QWORD PTR x 	\
    __asm FRNDINT             	\
    __asm FSUBR   QWORD PTR x 	\
    __asm FLD1                	\
    __asm FADDP   ST(1),ST	\
    __asm FSTP    QWORD PTR y 	\
    }
#else
#define SINCOS_TRIM_X(y,x) \
  y = x - rint(x) + 1.0;
#endif
#elif EAH_SINCOS_ROUND == EAH_SINCOS_ROUND_FLOOR 
#define SINCOS_TRIM_X(y,x) \
  y = x - floor(x);
#elif EAH_SINCOS_ROUND == EAH_SINCOS_ROUND_INT4
#define SINCOS_TRIM_X(y,x) \
  y = x-(INT4)x; \
  if ( y < 0.0 ) { y += 1.0; }
#elif EAH_SINCOS_ROUND == EAH_SINCOS_ROUND_INT8
#define SINCOS_TRIM_X(y,x) \
  y = x-(INT8)x; \
  if ( y < 0.0 ) { y += 1.0; }
#elif EAH_SINCOS_ROUND == EAH_SINCOS_ROUND_MODF
#define SINCOS_TRIM_X(y,x) \
{ \
  REAL8 dummy; \
  y = modf(x, &dummy); \
  if ( y < 0.0 ) { y += 1.0; } \
}
#endif

/* sin/cos Lookup tables */
#if (EAH_SINCOS_VARIANT == EAH_SINCOS_VARIANT_LAL)
#define SINCOS_LUT_RES 64 /* resolution of lookup-table */
static REAL4 sinLUT[SINCOS_LUT_RES+1];
static REAL4 cosLUT[SINCOS_LUT_RES+1];
#elif (EAH_SINCOS_VARIANT == EAH_SINCOS_VARIANT_LINEAR)
#define SINCOS_LUT_RES 1024 /* should be multiple of 4 */
static REAL4 sincosLUTbase[SINCOS_LUT_RES+SINCOS_LUT_RES/4];
static REAL4 sincosLUTdiff[SINCOS_LUT_RES+SINCOS_LUT_RES/4];
#endif

#define SINCOS_LUT_RES_F	(1.0 * SINCOS_LUT_RES)
#define OO_SINCOS_LUT_RES	(1.0 / SINCOS_LUT_RES)

#define X_TO_IND	(1.0 * SINCOS_LUT_RES * OOTWOPI )
#define IND_TO_X	(LAL_TWOPI * OO_SINCOS_LUT_RES)

#if (EAH_SINCOS_VARIANT == EAH_SINCOS_VARIANT_LAL)
/* "traditional" version with 2 LUT */

static void local_sin_cos_2PI_LUT_init (void)
{
  UINT4 k;
  static REAL8 const oo_lut_res = OO_SINCOS_LUT_RES;
  for (k=0; k <= SINCOS_LUT_RES; k++) {
    sinLUT[k] = sin( LAL_TWOPI * k * oo_lut_res );
    cosLUT[k] = cos( LAL_TWOPI * k * oo_lut_res );
  }
}

static int local_sin_cos_2PI_LUT_trimmed (REAL4 *sin2pix, REAL4 *cos2pix, REAL8 x)
{
  INT4 i0;
  REAL8 d, d2;
  REAL8 ts, tc;
  static REAL8 const oo_lut_res = OO_SINCOS_LUT_RES;

#ifndef LAL_NDEBUG
  if ( x < 0.0 || x >= 1.0 )
    {
      XLALPrintError("\nFailed numerica in local_sin_cos_2PI_LUT(): x = %f not in [0,1)\n\n", x );
      return XLAL_FAILURE;
    }
#endif

  i0 = (INT4)( x * SINCOS_LUT_RES_F + 0.5 );	/* i0 in [0, SINCOS_LUT_RES ] */
  d = d2 = LAL_TWOPI * (x - oo_lut_res * i0);
  d2 *= 0.5 * d;

  ts = sinLUT[i0];
  tc = cosLUT[i0];
   
  /* use Taylor-expansions for sin/cos around LUT-points */
  (*sin2pix) = ts + d * tc - d2 * ts;
  (*cos2pix) = tc - d * ts - d2 * tc;

  return XLAL_SUCCESS;
} /* local_sin_cos_2PI_LUT() */



#elif (EAH_SINCOS_VARIANT == EAH_SINCOS_VARIANT_LINEAR)

/* linear approximation developed with Akos */

#define SINCOS_ADDS  402653184.0
#define SINCOS_MASK1 0xFFFFFF
#define SINCOS_MASK2 0x003FFF
#define SINCOS_SHIFT 14

static void local_sin_cos_2PI_LUT_init (void)
{
  static const REAL8 step = LAL_TWOPI / (REAL8)SINCOS_LUT_RES;
  static const REAL8 div  = 1.0 / ( 1 << SINCOS_SHIFT );
  REAL8 start, end, true_mid, linear_mid;
  int i;

  start = 0.0; /* sin(0 * step) */
  for( i = 0; i < SINCOS_LUT_RES + SINCOS_LUT_RES/4; i++ ) {
    true_mid = sin( ( i + 0.5 ) * step );
    end = sin( ( i + 1 ) * step );
    linear_mid = ( start + end ) * 0.5;
    sincosLUTbase[i] = start + ( ( true_mid - linear_mid ) * 0.5 );
    sincosLUTdiff[i] = ( end - start ) * div;
    start = end;
  }
}


static int local_sin_cos_2PI_LUT_trimmed (REAL4 *sin2pix, REAL4 *cos2pix, REAL8 x) {
  INT4  i, n, ix;
  union {
    REAL8 asreal;
    struct {
      INT4 high;
      INT4 low;
    } as2int;
  } ux;

  static const REAL4* cosbase = sincosLUTbase + (SINCOS_LUT_RES/4);
  static const REAL4* cosdiff = sincosLUTdiff + (SINCOS_LUT_RES/4);

#ifndef LAL_NDEBUG
  if(x > SINCOS_ADDS) {
    LogPrintf(LOG_DEBUG,"sin_cos_LUT: x too large: %22f > %f\n",x,SINCOS_ADDS);
    return XLAL_FAILURE;
  } else if(x < -SINCOS_ADDS) {
    LogPrintf(LOG_DEBUG,"sin_cos_LUT: x too small: %22f < %f\n",x,-SINCOS_ADDS);
    return XLAL_FAILURE;
  }
#endif

#if EAH_SINCOS_F2IBITS == EAH_SINCOS_F2IBITS_MEMCPY
  ux.asreal = x + SINCOS_ADDS;
  memcpy(&(ix), &ux.asreal, sizeof(ix));
#elif EAH_SINCOS_F2IBITS == EAH_SINCOS_F2IBITS_UNION
  ux.asreal = x + SINCOS_ADDS;
#ifdef __BIG_ENDIAN__
  ix = ux.as2int.low;
#else /* BIG_ENDIAN */
  ix = ux.as2int.high;
#endif /* BIG_ENDIAN */
#endif /* EAH_SINCOS_F2IBITS */

  i  = ix & SINCOS_MASK1;
  n  = ix & SINCOS_MASK2;
  i  = i >> SINCOS_SHIFT;
  
  (*sin2pix) = sincosLUTbase[i]  + n * sincosLUTdiff[i];
  (*cos2pix) = cosbase[i]        + n * cosdiff[i];

  return XLAL_SUCCESS;
}

#else  /* EAH_SINCOS_VARIANT */
#error no valid EAH_SINCOS_VARIANT specified
#endif /* EAH_SINCOS_VARIANT */
