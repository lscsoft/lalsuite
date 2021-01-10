
#include <stdio.h>
#include <math.h>
#include <lal/LatticeTiling.h>
#include <lal/PiecewiseModel.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif


typedef struct{
  double radius;
  int pm;
} CircleBoundInfo;


double CircleBound(
  const void *data,
  const size_t dim UNUSED,
  const gsl_matrix *cache UNUSED,
  const gsl_vector* point
  )
{
  /* I think this defines a YBoundInfo struct with the values in data */
  const CircleBoundInfo* info = (const CircleBoundInfo*)data;
  
  /* This extracts the value located at the 0th position in the gsl_vector point */
  const double x = gsl_vector_get(point, 0);
  
  /* Define important variables to determine the bounds */
  
  double radius = info->radius;
  int pm = info->pm;
	
  double y = pm * sqrt(radius * radius - x * x);

  return y;
}


int XLALSetLatticeTilingCircleBound(
  LatticeTiling* tiling,
  const double radius
  )
{
  /* Simple checks. Not sure what I assume the types of errors are, something to ask about later */
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(radius > 0.0, XLAL_EINVAL);
  
  /* Build parameter space bounds */
  CircleBoundInfo XLAL_INIT_DECL( info_lower );
  CircleBoundInfo XLAL_INIT_DECL( info_upper );
  info_lower.radius = info_upper.radius = radius;
  info_lower.pm = -1;
  info_upper.pm = 1;
  
  
  /* Set bounds and check for success? */
  XLAL_CHECK(XLALSetLatticeTilingBound(tiling, 1, CircleBound, sizeof( info_lower ), &info_lower, &info_upper) == XLAL_SUCCESS, XLAL_EFAILED);
  
  return XLAL_SUCCESS;
}

/* Data needed to define bounds for a piecewise segment */
typedef struct{
  double nmin;
  double nmax;
  double kmin;
  double kmax;
  double segmentlength; /* Length of piecewise segment of interest */
  double nf0; /* Previous frequency parameter */
  double nf1; /* Previous frequency derivative parameter */
  double nf2; /* Previous frequency double derivative parameter */
} PiecewiseSegmentBoundsInfo;



int XLALSetLatticeTilingPiecewiseBounds(
  LatticeTiling* tiling,
  const double fmin,
  const double fmax,
  const double nmin,
  const double nmax,
  const double kmin,
  const double kmax,
  const double knots[]
  )
{
  /* Simple checks. Not sure what I assume the types of errors are, something to ask about later */
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(radius > 0.0, XLAL_EINVAL);
  
  
  /* First knot bounds */
  SetLatticeTilingConstantBound(tiling, 0, fmin, fmax);
  FirstKnotBoundInfo XLAL_INIT_DECL( info_first_knot_lower );
  FirstKnotBoundInfo XLAL_INIT_DECL( info_first_knot_upper );
  
  info_first_knot_lower.n = nmax;
  info_first_knot_lower.k = kmax;
  
  info_first_knot_upper.n = nmin;
  info_first_knot_upper.k = kmin;
  
  XLAL_CHECK(XLALSetLatticeTilingBound(tiling, 1, FirstKnotDerivBound, sizeof( info_first_knot_lower ), &info_first_knot_lower, &info_first_knot_upper) == XLAL_SUCCESS, XLAL_EFAILED)
  XLAL_CHECK(XLALSetLatticeTilingBound(tiling, 2, FirstKnotDerivBound, sizeof( info_first_knot_lower ), &info_first_knot_lower, &info_first_knot_upper) == XLAL_SUCCESS, XLAL_EFAILED)
  
  int knot;
  int finalknot = sizeof(knots)/sizeof(knot);
  
  for (knot = 1; knot <= finalknot; ++knot){
    segmentlength = knots[knot] - knots[knot - 1];
    
    PiecewiseBoundInfo XLAL_INIT_DECL( info_knot_lower );
    PiecewiseBoundInfo XLAL_INIT_DECL( info_knot_upper );
    
    info_knot_lower.nmin = info_knot_upper.nmin = nmin;
    info_knot_lower.nmax = info_knot_upper.nmax = nmax;
    info_knot_lower.ntol = info_knot_upper.ntol = ntol;
    info_knot_lower.kmin = info_knot_upper.kmin = kmin;
    info_knot_lower.kmax = info_knot_upper.kmax = kmin;
    info_knot_lower.segmentlength = info_knot_upper.segmentlength = segmentlength;
    
    
    info_knot_lower.upperlower = -1;
    info_knot_upper.upperlower = 1;
    
    dimindex = 3 * knot
    
    XLAL_CHECK(XLALSetLatticeTilingBound(tiling, dimindex, F0Bound, sizeof( info_knot_lower ), &info_knot_lower, &info_knot_upper) == XLAL_SUCCESS, XLAL_EFAILED)
    XLAL_CHECK(XLALSetLatticeTilingBound(tiling, dimindex + 1, F1Bound, sizeof( info_knot_lower ), &info_knot_lower, &info_knot_upper) == XLAL_SUCCESS, XLAL_EFAILED)
    XLAL_CHECK(XLALSetLatticeTilingBound(tiling, dimindex + 2, F2Bound, sizeof( info_knot_lower ), &info_knot_lower, &info_knot_upper) == XLAL_SUCCESS, XLAL_EFAILED)
    
  }
  
  return XLAL_SUCCESS;
}

typedef struct{
  double n;
  double k;
} FirstKnotBoundInfo;


/* Return bounds on the f1 and f2 parameters on the first knot */
double FirstKnotDerivBound(
  const void *data,
  const size_t dim,
  const gsl_matrix *cache UNUSED,
  const gsl_vector* point
  )
{
  /* I think this defines a YBoundInfo struct with the values in data */
  const FirstKnotBoundInfo* info = (const FirstKnotBoundInfo*)data;
  
  /* This extracts the value located at the 0th position in the gsl_vector point */
  const double f0 = gsl_vector_get(point, 0);
  
  double n = info->n;
  double k = info->k;
  
  if (dim == 1){
    double f1 = -k * pow(f0, n);
    return f1;
  }
  else if (dim == 2){
    const double f1 = gsl_vector_get(point, 1);
    double f2 = n * f1 * f1 / f0;
    return f2;
  }
}


double GTEAndDerivs(
  double f0,
  double n,
  double k,
  double t,
  int d
  )
{
  double base = 1 + (n + 1) * k * pow(f0, n - 1) * t;
  double power = 1 / (1 - n);
  
  if (d == 0){
    return f0 * pow(base, power);
  }
  else if (d == 1){
    double factor = -1 * pow(f0, n) * k;
    return factor * pow(base, power - 1);
  }
  else if (d == 1){
    double factor = power(f0, 2 * n - 1) * k * k * n;
    return factor * pow(base, power - 2);
  }
}


double * GetLastElements(
  int initialindex,
  double elems[]
  )
{
  int lastindex = sizeof(elems)/sizeof(elems[0]) - 1;
  
  static double lastelems[initialindex];
  
  int i;
  int j = 0;
  for (i == lastindex - initialindex; i <= lastindex; ++i){
    lastelems[j] = elems[i];
    ++j;
  }
  
  return lastelems;
  
}


double F0BoundMinMax(
  double f0,
  double na, /* For calculating upper bound, na > nb. For calculating lower bound na < nb */
  double nb, /* For calculating upper bound, na > nb. For calculating lower bound na < nb */
  double ka, /* For calculating upper bound, ka > kb. For calculating lower bound ka < kb */
  double kb, /* For calculating upper bound, ka > kb. For calculating lower bound ka < kb */
  double seglength,
  int minmax
)
{
  double gted = GTEAndDerivs(f0, na, ka, seglength, 1)
  double gtedd = GTEAndDerivs(f0, nb, kb, seglength, 2)
  
  double brakingindex = na * gted / gtedd
  
  double paramrange = GTEAndDerivs(f0, nb, kb, seglength, 0)
  
  
  if (brakingindex < paramrange){
    if (minmax == 1){
      return brakingindex
    }
    else if (minmax == -1){
      return paramrange
    }
  }
  else if (brakingindex >= paramrange){
    if (minmax == 1){
      return paramrange
    }
    else if (minmax == -1){
      return brakingindex
    }
  }
}

double * NMinMax(
  double n,
  double ntol,
  double nmin,
  double nmax,
  double segmentlength UNUSED
  )
{
  static double nminmax[2];
  
  double nextnmin = n * (1 - ntol);
  double nextnmax = n * (1 + ntol);
  
  if (nextnmin < nmin){
    nextnmin = nmin;
  }
  if (nextnmax > nmax){
    nextnmax = nmax;
  }
  
  nminmax[0] = nextnmin;
  nminmax[1] = nextnmax;
  
  return nminmax;
}

typedef struct(
  double nmin;
  double nmax;
  double ntol;
  double kmin;
  double kmax;
  double segmentlength;
  int upperlower
) PiecewiseBoundInfo;


double F0Bound(
  const void *data,
  const size_t dim UNUSED
  const gsl_matrix *cache UNUSED,
  const gsl_vector* point
  )
{
  double prevknot = GetLastElements(point.data, 3);
  double nprev = prevknot[2] * prevknot[1] / prevknot[0];
  double f0 = prevknot[0]
  
  /* I think this defines a YBoundInfo struct with the values in data */
  const PiecewiseBoundInfo* info = (const PiecewiseBoundInfo*)data;
  
  double nmin = PiecewiseBoundInfo.nmin;
  double nmax = PiecewiseBoundInfo.nmax;
  double ntol = PiecewiseBoundInfo.ntol;
  double kmin = PiecewiseBoundInfo.kmin;
  double kmax = PiecewiseBoundInfo.kmax
  double segmentlength = PiecewiseBoundInfo.segmentlength;
  int upperlower = PiecewiseBoundInfo.upperlower;
  
  double nminmax[2] = NMinMax(nprev, ntol, nmin, nmax, segmentlength);
  
  if (upperlower == 1){
    double upperbound = F0BoundMinMax(f0, nminmax[2], nminmax[1], kmax, kmin, segmengthlength, 1);
    return upperbound;
  }
  else if (upperlower == 2){
    double lowerbound = F0BoundMinMax(f0, nminmax[1], nminmax[2], kmin, kmax, segmengthlength, -1);
    return lowerbound;
  }
}


double F1BoundMinMax(
  double f0,
  double na, /* For upper bound, na < nb */
  double nb
  double ka, /* For upper bound ka < kb. */
  double kb,
  double seglength,
  int minmax
)
{
  double prangecondition1 = -ka * pow(f0, na);
  double prangecondition2 = GTEAndDerivs(f0, na, ka, seglength, 1);
  
  double brakingindexcriteria = GTEAndDerivs(f0, nb, kb, seglength, 2) * f0 / nb;
  
  double bound = 0;
  
  if (minmax == 1){
    if (prangecondition1 > prangecondition2){
      bound = prangecondition1;
    }
    else {
      bound = prangecondition2;
    }
    if (bound > brakingindexcriteria){
      return bound
    }
    else {
      return brakingindexcriteria
    }
  }
  else if (minmax == -1){
    if (prangecondition1 < prangecondition2){
      bound = prangecondition1;
    }
    else {
      bound = prangecondition2;
    }
    if (bound < brakingindexcriteria){
      return bound
    }
    else {
      return brakingindexcriteria
    }
  }
}


double F1Bound(
  const void *data,
  const size_t dim UNUSED
  const gsl_matrix *cache UNUSED,
  const gsl_vector* point
  )
{
  double prevknot = GetLastElements(point.data, 4);
  double nprev = prevknot[2] * prevknot[1] / prevknot[0];
  double f0 = prevknot[3]
  
  /* I think this defines a YBoundInfo struct with the values in data */
  const PiecewiseBoundInfo* info = (const PiecewiseBoundInfo*)data;
  
  double nmin = PiecewiseBoundInfo.nmin;
  double nmax = PiecewiseBoundInfo.nmax;
  double ntol = PiecewiseBoundInfo.ntol;
  double kmin = PiecewiseBoundInfo.kmin;
  double kmax = PiecewiseBoundInfo.kmax
  double segmentlength = PiecewiseBoundInfo.segmentlength;
  int upperlower = PiecewiseBoundInfo.upperlower;
  
  double nminmax[2] = NMinMax(nprev, ntol, nmin, nmax, segmentlength);
  
  if (upperlower == 1){
    double upperbound = F1BoundMinMax(f0, nminmax[1], nminmax[2], kmin, kmax, segmengthlength, 1);
    return upperbound;
  }
  else if (upperlower == 2){
    double lowerbound = F1BoundMinMax(f0, nminmax[2], nminmax[1], kmax, kmin, segmengthlength, -1);
    return lowerbound;
  }
}


double F2BoundMinMax(
  double f0,
  double f1,
  double n, /* For upper bound n should be maximised */
  double k,
  double seglength,
  int minmax
)
{
  double brakingindexcriteria = n * pow(f1, 2) / f0;
  double prangecondition = GTEAndDerivs(f0, n, k, seglength, 2);
 
  if (minmax == 1){
    if (brakingindexcriteria > prangecriteria) {
      return brakingindexcriteria
    }
    else {
      return prangecriteria;
    }
  }
  else if (minmax == -1) {
    if (brakingindexcriteria < prangecriteria) {
      return brakingindexcriteria
    }
    else {
      return prangecriteria;
    }
  }
}

double F2Bound(
  const void *data,
  const size_t dim UNUSED
  const gsl_matrix *cache UNUSED,
  const gsl_vector* point
  )
{
  double prevknot = GetLastElements(point.data, 5);
  double nprev = prevknot[2] * prevknot[1] / prevknot[0];
  double f0 = prevknot[3];
  double f1 = prevknot[4]
  
  /* I think this defines a YBoundInfo struct with the values in data */
  const PiecewiseBoundInfo* info = (const PiecewiseBoundInfo*)data;
  
  double nmin = PiecewiseBoundInfo.nmin;
  double nmax = PiecewiseBoundInfo.nmax;
  double ntol = PiecewiseBoundInfo.ntol;
  double kmin = PiecewiseBoundInfo.kmin;
  double kmax = PiecewiseBoundInfo.kmax
  double segmentlength = PiecewiseBoundInfo.segmentlength;
  int upperlower = PiecewiseBoundInfo.upperlower;
  
  double nminmax[2] = NMinMax(nprev, ntol, nmin, nmax, segmentlength);
  
  if (upperlower == 1){
    double upperbound = F2BoundMinMax(f0, f1, nminmax[2], kmax, segmentlength, 1);
    return upperbound;
  }
  else if (upperlower == 2){
    double lowerbound = F2BoundMinMax(f0, f1, nminmax[1], kmin, segmentlength, -1);
    return lowerbound;
  }
}



















