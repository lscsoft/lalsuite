
#include <stdio.h>
#include <math.h>
#include <lal/LatticeTiling.h>
#include <lal/PiecewiseModel.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

///
/// Information required for creating circular bounds
///
typedef struct{
  double radius; /// Radius of circle
  int pm;        /// +1 for upper bound, -1 for lower bound
} CircleBoundInfo;

///
/// Function which describes the circular bound
///
static double CircleBound(
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
  
  if (x > 1 || x < 1){
    return 0.;
  }
  
  double y = pm * sqrt(radius * radius - x * x);
  
  return y;
}

///
/// Sets the bound for a circular parameter space
///
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


typedef struct tagFirstKnotBoundInfo{
  double n;  /// A braking index (should be of an optimised value)
  double k;  /// A k value (should be of an optimised value)
} FirstKnotBoundInfo;



///
/// Defines the bounds for the first and second derivative frequency parameters on the first knot (t = 0)
///
static double FirstKnotDerivBound(
  const void *data,
  const size_t dim,
  const gsl_matrix *cache UNUSED,
  const gsl_vector* point
  )
{
  const FirstKnotBoundInfo* info = (const FirstKnotBoundInfo*)data;
  
  const double f0 = gsl_vector_get(point, 0);
  
  double n = info->n;
  double k = info->k;
  
  if(f0 <= 0){
    return 0.;
  }
  
  if (dim == 1){
    double f1 = -k * pow(f0, n);
    
    return f1;
  }
  else if (dim == 2){
    const double f1 = gsl_vector_get(point, 1);
    
    if(f1 >= 0){
      return 0.;
    }
    
    double f2 = n * pow(f1, 2) / f0;
    
    double thisn = f2 * f0 / pow(f1, 2);
    printf("%f \n", thisn);
    
    return f2;
  }

  return NAN;
}





///
/// A struct containing the relevant information for calculating upper and lower bounds on a specific piecewise segment
///
typedef struct tagPiecewiseBoundInfo{
  double nmin;  /// Global minimum braking index
  double nmax;  /// Global maximum braking index
  double ntol;  /// Braking index tolerance (percentage) between adjacent knots
  double kmin;  /// Global minimum k value
  double kmax;  /// Global maximum k value
  double ktol;  /// k value tolerance (percentage) between adjacent knots
  double segmentlength;  /// The time difference between the two knots relevant to this piecewise segment
  int upperlower;  /// +1 for calculating upper bound, -1 for calculating lower bound
} PiecewiseBoundInfo;


///
/// The general torque equation and its first two derivatives 
///
static double GTEAndDerivs(
  double f0, /// Initial frequency
  double n,  /// Braking index
  double k,  /// k value
  double t,  /// Time at which to calculat evalue
  int d      /// Derivative order (d <= 2)
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
  else if (d == 2){
    double factor = pow(f0, 2 * n - 1) * k * k * n;
    double secondderiv = factor * pow(base, power - 2);
    
    return secondderiv;
    
  }
  
  return 0.0;
}

///
/// Calculates the valid braking index for the next knot based upon the previous braking index. Arbitrary definition
///
static double * NMinMax(
  double n,     /// Previous braking index
  double ntol,  /// Braking index tolerance (percentage)
  double nmin,  /// Minimum acceptable braking index
  double nmax,  /// Maximum acceptable braking index
  double segmentlength UNUSED /// Time difference between current knot and previous knot. Kept in case function definition changes at a later date
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


///
/// Calculates the valid k value for the next knot based upon the previous value. Arbitrary definition
///
static double * KMinMax(
  double k,     /// Previous braking index
  double ktol,  /// Braking index tolerance (percentage)
  double kmin,  /// Minimum acceptable braking index
  double kmax,  /// Maximum acceptable braking index
  double segmentlength UNUSED /// Time difference between current knot and previous knot. Kept in case function definition changes at a later date
  )
{
  static double kminmax[2];
  
  double nextkmin = k * (1 - ktol);
  double nextkmax = k * (1 + ktol);
  
  if (nextkmin < kmin){
    nextkmin = kmin;
  }
  if (nextkmax > kmax){
    nextkmax = kmax;
  }
  
  kminmax[0] = nextkmin;
  kminmax[1] = nextkmax;
  
  return kminmax;
}



///
/// Calculates the minimum and maximum bounds for a frequency frequency parameter
///
static double F0BoundMinMax(
  double f0, /// Frequency value of the previous knot
  double na, /// A braking index value. For calculating upper bound, na > nb. For calculating lower bound na < nb
  double nb, /// A braking index value. For calculating upper bound, na > nb. For calculating lower bound na < nb
  double ka, /// A k value. For calculating upper bound, ka > kb. For calculating lower bound ka < kb
  double kb, /// A k value. For calculating upper bound, ka > kb. For calculating lower bound ka < kb
  double seglength, /// Time difference between this knot and the previous knot
  int minmax /// +1 for calculating the maximum bound, -1 for calculating the minimum bound
)
{
  /// These combinations of gted and gtedd give one of the optimised values of f0 (either minimum or maximum) such that the braking index at the next knot still
  /// falls within an acceptable range. This acceptable range is between [na, nb] or [nb, na], depending upon the value of minmax. Do a double check of this maths
  
  printf("%f, %f, %f, %f, %f \n", f0, na, nb, ka, kb);
  
  double gted = GTEAndDerivs(f0, na, ka, seglength, 1);
  double gtedd = GTEAndDerivs(f0, nb, kb, seglength, 2);
  double brakingindex = na * pow(gted, 2) / gtedd;
  
  /// Parameter range optimisation for frequency parameter
  double paramrange = GTEAndDerivs(f0, nb, kb, seglength, 0);
  
  // Selecting the appropriate bound as the optimised value
  if (brakingindex < paramrange){
    if (minmax == 1){
      return brakingindex;
    }
    else if (minmax == -1){
      return paramrange;
    }
  }
  else if (brakingindex >= paramrange){
    if (minmax == 1){
      return paramrange;
    }
    else if (minmax == -1){
      return brakingindex;
    }
  }
  
  return 0.0;
}


///
/// Sets the bound on the frequency parameter
///
static double F0Bound(
  const void *data,
  const size_t dim,
  const gsl_matrix *cache UNUSED,
  const gsl_vector* point
  )
{
  
  double f0n1 = gsl_vector_get(point, dim - 3);
  double f1n1 = gsl_vector_get(point, dim - 2);
  double f2n1 = gsl_vector_get(point, dim - 1);
  
  if(f0n1 <= 0 || f1n1 >= 0 || f2n1 <= 0){
    return 0.;
  }
  
  const PiecewiseBoundInfo* info = (const PiecewiseBoundInfo*)data;
  
  double nmin = info->nmin;
  double nmax = info->nmax;
  double ntol = info->ntol;
  double kmin = info->kmin;
  double kmax = info->kmax;
  double ktol = info->ktol;
  double segmentlength = info->segmentlength;
  int upperlower = info->upperlower;
  
  double nprev = f2n1 * f0n1 / pow(f1n1, 2);
  
  if (nprev > nmax){
    nprev = nmax;
  }
  else if (nprev < nmin){
    nprev = nmin;
  }
  
  double kprev = - f1n1 / pow(f0n1, nprev);
  
  if (kprev > kmax){
    kprev = kmax;
  }
  else if (kprev < kmin){
    kprev = kmin;
  }
  
  double *nminmax = NMinMax(nprev, ntol, nmin, nmax, segmentlength);
  double *kminmax = KMinMax(kprev, ktol, kmin, kmax, segmentlength);
  
  ///printf("%f \n", kprev);
  
  if (upperlower == 1){
    double upperbound = F0BoundMinMax(f0n1, nminmax[1], nminmax[0], kminmax[1], kminmax[0], segmentlength, 1);
    return upperbound;
  }
  else if (upperlower == -1){
    double lowerbound = F0BoundMinMax(f0n1, nminmax[0], nminmax[1], kminmax[0], kminmax[1], segmentlength, -1);
    return lowerbound;
  }
  
  return 0.0;
}

///
/// Function which determines the minimum/maximum value of the first derivative parameters
///
static double F1BoundMinMax(
  double f0,        /// Frequency parameter at the corresponding knot
  double fprev,     /// Frequency parameter at the previous knot
  double na,        /// A braking index value. For calculating upper bound, na < nb. For calculating lower bound na > nb
  double nb,        /// A braking index value. For calculating upper bound, na < nb. For calculating lower bound na > nb
  double ka,        /// A k value. For calculating upper bound, ka < kb. For calculating lower bound ka > kb
  double kb,        /// A k value. For calculating upper bound, ka < kb. For calculating lower bound ka > kb
  double seglength, /// Time difference between this knot and the previous
  int minmax        /// +1 for maximum bound, -1 for minimum bound
)
{
  double prangecondition1 = -ka * pow(f0, na);
  double prangecondition2 = GTEAndDerivs(fprev, na, ka, seglength, 1);
  
  double brakingindexcriteria = sqrt(GTEAndDerivs(fprev, nb, kb, seglength, 2) * f0 / nb);
  double bound = 0;
  
  
  if (minmax == 1){
    if (prangecondition1 > prangecondition2){
      bound = prangecondition2;
    }
    else {
      bound = prangecondition1;
    }
    if (bound > brakingindexcriteria){
      return brakingindexcriteria;
    }
    else {
      return bound;
    }
  }
  else if (minmax == -1){
    if (prangecondition1 < prangecondition2){
      bound = prangecondition2;
    }
    else {
      bound = prangecondition1;
    }
    if (bound < brakingindexcriteria){
      return brakingindexcriteria;
    }
    else {
      return bound;
    }
  }
  
  return 0.0;
}

///
/// Sets the bound on the first derivative frequency parameter 
///
static double F1Bound(
  const void *data,
  const size_t dim,
  const gsl_matrix *cache UNUSED,
  const gsl_vector* point
  )
{
  double f0n1 = gsl_vector_get(point, dim - 4);
  double f1n1 = gsl_vector_get(point, dim - 3);
  double f2n1 = gsl_vector_get(point, dim - 2);
  
  if(f0n1 <= 0 || f1n1 >= 0 || f2n1 <= 0){
    return 0.;
  }
  
  double f0 = gsl_vector_get(point, dim - 1);
  
  if(f0 <= 0){
    return 0.;
  }
  
  ///printf("%f, %f, %f, %f \n", f0n1, f1n1, f2n1, f0);
  
  const PiecewiseBoundInfo* info = (const PiecewiseBoundInfo*)data;
  
  double nmin = info->nmin;
  double nmax = info->nmax;
  double ntol = info->ntol;
  double kmin = info->kmin;
  double kmax = info->kmax;
  double ktol = info->ktol;
  double segmentlength = info->segmentlength;
  int upperlower = info->upperlower;
  
  
  double nprev = f2n1 * f0n1 / pow(f1n1, 2);
  
  if (nprev > nmax){
    nprev = nmax;
  }
  else if (nprev < nmin){
    nprev = nmin;
  }
  
  double kprev = - f1n1 / pow(f0n1, nprev);
  
  if (kprev > kmax){
    kprev = kmax;
  }
  else if (kprev < kmin){
    kprev = kmin;
  }
  
  
  double *nminmax = NMinMax(nprev, ntol, nmin, nmax, segmentlength);
  double *kminmax = NMinMax(kprev, ktol, kmin, kmax, segmentlength);
  
  if (upperlower == 1){
    double upperbound = F1BoundMinMax(f0, f0n1, nminmax[0], nminmax[1], kminmax[0], kminmax[1], segmentlength, 1);
    
    return upperbound;
  }
  else if (upperlower == -1){
    double lowerbound = F1BoundMinMax(f0, f0n1, nminmax[1], nminmax[0], kminmax[1], kminmax[0], segmentlength, -1);
    
    return lowerbound;
  }
  
  return 0.0;
}

///
/// Calculates the bound for the second derivative frequency parameter 
///
static double F2BoundMinMax(
  double f0,        /// Frequency parameter on this knot
  double fprev,     /// Frequency parameter on the previous knot
  double f1,        /// First derivative frequency parameter on this knot
  double n,         /// Optimal braking index value. For upper bound n should be maximised, for lower bound n should be minimised
  double k,         /// Optimal k value. For upper bound k should be maximised, for lower bound k should be minimised
  double seglength, /// Time difference between this knot and the previous
  int minmax        /// +1 for maximum value, -1 for minimum value
)
{
  double brakingindexcriteria = n * pow(f1, 2) / f0;
  double prangecriteria = GTEAndDerivs(fprev, n, k, seglength, 2);
  
  if (minmax == 1){
    if (brakingindexcriteria > prangecriteria) {
      return prangecriteria;
    }
    else {
      return brakingindexcriteria;
    }
  }
  else if (minmax == -1) {
    if (brakingindexcriteria < prangecriteria) {
      return prangecriteria;
    }
    else {
      return brakingindexcriteria;
    }
  }
  
  return 0.0;
}

///
/// Sets the bounds on the second derivative frequency parameter 
///
static double F2Bound(
  const void *data,
  const size_t dim,
  const gsl_matrix *cache UNUSED,
  const gsl_vector* point
  )
{
  double f0n1 = gsl_vector_get(point, dim - 5);
  double f1n1 = gsl_vector_get(point, dim - 4);
  double f2n1 = gsl_vector_get(point, dim - 3);
  
  if(f0n1 <= 0 || f1n1 >= 0 || f2n1 <= 0){
    return 0.;
  }
  
  double f0 = gsl_vector_get(point, dim - 2);
  double f1 = gsl_vector_get(point, dim - 1);
  
  if(f0 <= 0 || f1 >= 0){
    return 0.;
  }
  
  const PiecewiseBoundInfo* info = (const PiecewiseBoundInfo*)data;
  
  double nmin = info->nmin;
  double nmax = info->nmax;
  double ntol = info->ntol;
  double kmin = info->kmin;
  double kmax = info->kmax;
  double ktol = info->ktol;
  double segmentlength = info->segmentlength;
  int upperlower = info->upperlower;
  
  double nprev = f2n1 * f0n1 / pow(f1n1, 2);
  
  if (nprev > nmax){
    nprev = nmax;
  }
  else if (nprev < nmin){
    nprev = nmin;
  }
  
  double kprev = - f1n1 / pow(f0n1, nprev);
  
  if (kprev > kmax){
    kprev = kmax;
  }
  else if (kprev < kmin){
    kprev = kmin;
  }
  
  double *nminmax = NMinMax(nprev, ntol, nmin, nmax, segmentlength);
  double *kminmax = NMinMax(kprev, ktol, kmin, kmax, segmentlength);
  
  if (upperlower == 1){
    double upperbound = F2BoundMinMax(f0, f0n1, f1, nminmax[1], kminmax[1], segmentlength, 1);
    return upperbound;
  }
  else if (upperlower == -1){
    double lowerbound = F2BoundMinMax(f0, f0n1, f1, nminmax[0], kminmax[0], segmentlength, -1);
    return lowerbound;
  }
  
  return 0.0;
}

///
/// Converts a given tau value to a k value (or a given k value to a tau value)
///
static double ktauconversion(
  double f0,  /// Frequency at time t = 0. For calculating kmax/kmin from initial tau values, f0 should be maximised. Likewise for calculating taumin/taumax from k values 
  double n,   /// A braking index. For calculating kmax/kmin from initial tau values, n should be maximised. Likewise for calculating taumin/taumax from k values 
  double kortau  /// The k or tau values we wish to convert. For calculating kmax, taumin should be used. For calculating kmin, taumax should be used.
  )
{
  double numerator = pow(2, n - 1) - 1;
  double denominator = (n - 1) * pow(f0, n - 1) * kortau;
  
  double ktau = numerator / denominator;
  
  return ktau;
}


///
/// Sets the bounds for the piecewise model
///
int XLALSetLatticeTilingPiecewiseBounds(
  LatticeTiling* tiling,
  const double fmin, /// Minimum spin frequency
  const double fmax, /// Maximum spin frequency
  const double nmin, /// Minimum braking index
  const double nmax, /// Maximum braking index
  const double ntol, /// Tolerance (percentage) between braking indices on adjacent knots
  const double taumin, /// Minimum tau value
  const double taumax, /// Maximum tau value
  const double ktol, /// Tolerance (percentage) between k values on adjacent knots
  const gsl_vector* knots, /// List of knots
  const int finalknot /// The number of the final knot
  )
{
  /* Simple checks. Not sure what I assume the types of errors are, something to ask about later */
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  
  /// Converting tau values to k values
  double kmin = ktauconversion(fmax, nmax, taumax);
  double kmax = ktauconversion(fmax, nmax, taumin);
  
  
  /// Setting the first knot bounds
  printf("Setting bounds on knot 0 \n");
  XLALSetLatticeTilingConstantBound(tiling, 0, fmin, fmax);
  FirstKnotBoundInfo XLAL_INIT_DECL( info_first_knot_lower );
  FirstKnotBoundInfo XLAL_INIT_DECL( info_first_knot_upper );
  
  info_first_knot_lower.n = nmax;
  info_first_knot_lower.k = kmax;
  
  info_first_knot_upper.n = nmin;
  info_first_knot_upper.k = kmin;
  
  XLAL_CHECK(XLALSetLatticeTilingBound(tiling, 1, FirstKnotDerivBound, sizeof( info_first_knot_lower ), &info_first_knot_lower, &info_first_knot_upper) == XLAL_SUCCESS, XLAL_EFAILED);
  XLAL_CHECK(XLALSetLatticeTilingBound(tiling, 2, FirstKnotDerivBound, sizeof( info_first_knot_lower ), &info_first_knot_upper, &info_first_knot_lower) == XLAL_SUCCESS, XLAL_EFAILED);
  
  /// Setting the bounds for all following knots
  
  for (int knot = 1; knot < finalknot; ++knot){
    printf("Setting bounds on knot %i \n", knot);
    double segmentlength = gsl_vector_get(knots, knot) - gsl_vector_get(knots, knot - 1);
    
    PiecewiseBoundInfo XLAL_INIT_DECL( info_knot_lower );
    PiecewiseBoundInfo XLAL_INIT_DECL( info_knot_upper );
    
    info_knot_lower.nmin = info_knot_upper.nmin = nmin;
    info_knot_lower.nmax = info_knot_upper.nmax = nmax;
    info_knot_lower.ntol = info_knot_upper.ntol = ntol;
    info_knot_lower.kmin = info_knot_upper.kmin = kmin;
    info_knot_lower.kmax = info_knot_upper.kmax = kmin;
    info_knot_lower.ktol = info_knot_upper.ktol = ktol;
    info_knot_lower.segmentlength = info_knot_upper.segmentlength = segmentlength;
    
    info_knot_lower.upperlower = -1;
    info_knot_upper.upperlower = 1;
    
    int dimindex = 3 * knot;
    
    XLAL_CHECK(XLALSetLatticeTilingBound(tiling, dimindex,     F0Bound, sizeof( info_knot_lower ), &info_knot_lower, &info_knot_upper) == XLAL_SUCCESS, XLAL_EFAILED);
    XLAL_CHECK(XLALSetLatticeTilingBound(tiling, dimindex + 1, F1Bound, sizeof( info_knot_lower ), &info_knot_lower, &info_knot_upper) == XLAL_SUCCESS, XLAL_EFAILED);
    XLAL_CHECK(XLALSetLatticeTilingBound(tiling, dimindex + 2, F2Bound, sizeof( info_knot_lower ), &info_knot_lower, &info_knot_upper) == XLAL_SUCCESS, XLAL_EFAILED);
    
    printf("Finished that knot \n");
    printf("\n");
    
  }
  return XLAL_SUCCESS;
}














