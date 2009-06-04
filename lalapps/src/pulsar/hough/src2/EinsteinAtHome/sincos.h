/** the way of trimming x to the interval [0..1) for the sin_cos_LUT functions
    give significant differences in speed, so we provide various ways here.
    We also record the way we are using for logging */

#if EAH_SINCOS_ROUND == EAH_SINCOS_ROUND_PLUS2
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
#else /* C99, no MSC */
#define SINCOS_TRIM_X(y,x) \
  y = x - rint(x) + 1.0;
#endif /* C99 / MSC */

#elif EAH_SINCOS_ROUND == EAH_SINCOS_ROUND_FLOOR 
#define SINCOS_TRIM_X(y,x) \
  y = x - floor(x);
#endif

/* sin/cos Lookup tables */
#define SINCOS_LUT_RES 1024 /* should be multiple of 4 */

static REAL4 sincosLUTbase[SINCOS_LUT_RES+SINCOS_LUT_RES/4];
static REAL4 sincosLUTdiff[SINCOS_LUT_RES+SINCOS_LUT_RES/4];

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


/* x must already been trimmed to interval [0..2) */
static int local_sin_cos_2PI_LUT_trimmed (REAL4 *sin2pix, REAL4 *cos2pix, REAL8 x) {
  static const REAL4* cosbase = sincosLUTbase + (SINCOS_LUT_RES/4);
  static const REAL4* cosdiff = sincosLUTdiff + (SINCOS_LUT_RES/4);
  INT4  i, n;
  union {
    REAL8 asreal;
    struct {
      INT4 high;
      INT4 low;
    } as2int;
  } ux;

#ifndef LAL_NDEBUG
  if(x > SINCOS_ADDS) {
    LogPrintf(LOG_DEBUG,"sin_cos_LUT: x too large: %22f > %f\n",x,SINCOS_ADDS);
    return XLAL_FAILURE;
  } else if(x < -SINCOS_ADDS) {
    LogPrintf(LOG_DEBUG,"sin_cos_LUT: x too small: %22f < %f\n",x,-SINCOS_ADDS);
    return XLAL_FAILURE;
  }
#endif

  ux.asreal = x + SINCOS_ADDS;
#ifdef __BIG_ENDIAN__
  i  = ux.as2int.low & SINCOS_MASK1;
  n  = ux.as2int.low & SINCOS_MASK2;
#else
  i  = ux.as2int.high & SINCOS_MASK1;
  n  = ux.as2int.high & SINCOS_MASK2;
#endif
  i  = i >> SINCOS_SHIFT;
  
  (*sin2pix) = sincosLUTbase[i] + n * sincosLUTdiff[i];
  (*cos2pix) = cosbase[i]       + n * cosdiff[i];

  return XLAL_SUCCESS;
}
