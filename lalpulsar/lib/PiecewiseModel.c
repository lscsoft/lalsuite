
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
  double fmin;  /// Global minimum starting frequency
  double fmax;  /// Global maximum starting frequency
  double n;     /// A braking index (should be of an optimised value)
  double k;     /// A k value (should be of an optimised value)
} FirstKnotBoundInfo;


///
/// A struct containing the relevant information for calculating upper and lower bounds on a specific piecewise segment
///
typedef struct tagPiecewiseBoundInfo{
  double fmin;  /// Minimum frequency at t = 0
  double fmax;  /// Maximum frequency at t = 0
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
  
  double gted = GTEAndDerivs(f0, na, ka, seglength, 1);
  double gtedd = GTEAndDerivs(f0, nb, kb, seglength, 2);
  double brakingindex = na * pow(gted, 2) / gtedd;
  
  /// Parameter range optimisation for frequency parameter
  double paramrange = GTEAndDerivs(f0, nb, kb, seglength, 0);
  
  printf("These are the f0 options %E, %E \n", brakingindex, paramrange);
  printf("Min max is %d \n", minmax);
  
  if (minmax == 1) {
    if (brakingindex <= paramrange){
      return brakingindex;
    }
    else if (brakingindex > paramrange){
      return paramrange;
    }
  }
  else if (minmax == -1){
    if (brakingindex <= paramrange){
      return paramrange;
    }
    else if (paramrange < brakingindex){
      return brakingindex;
    }
  }
  
  return NAN;
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
  
  double brakingindexcriteria = -sqrt(GTEAndDerivs(fprev, na, kb, seglength, 2) * f0 / nb);
  double bound;
  
  ///printf("These are the f1 options %E, %E, %E \n", prangecondition1, prangecondition2, brakingindexcriteria);
  
  
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
  
  ///printf("For f2, these are the two criteria %E, %E \n", brakingindexcriteria, prangecriteria);
  
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


/// Takes an input val, and returns val to its closest value within the range [valmin, valmax].
static double resetinsidebounds(
  double val,     /// Value
  double valmin,  /// Lower bound of range
  double valmax   /// Upper bound of range
  )
{
  if (valmin <= val && val <= valmax){
    return val;
  }
  else if (val < valmin){
    printf("I am activated! Too small %E, %E, %E \n", val, valmin, valmax);
    return valmin;
  }
  else if (val > valmax){
    printf("I am activated! Too big %E, %E, %E \n", val, valmin, valmax);
    return valmax;
  }
  
  return NAN;
}


static void resetdimonpoint(
  gsl_vector* point,
  int dim,
  double fmin,
  double fmax,
  double nmin,
  double nmax,
  double ntol,
  double kmin,
  double kmax,
  double ktol,
  double segmentlength
  )
{
  if (dim == 0){
    double f0 = gsl_vector_get(point, dim);
    double val = resetinsidebounds(f0, fmin, fmax);
    gsl_vector_set(point, dim, val);
    return;
  }
  
  else if (dim == 1){
    double f0 = gsl_vector_get(point, dim - 1);
    double f1 = gsl_vector_get(point, dim);
    double lower = -kmax * pow(f0, nmax);
    double upper = -kmin * pow(f0, nmin);
    
    double val = resetinsidebounds(f1, lower, upper);
    gsl_vector_set(point, dim, val);
    return;
  }
  else if (dim == 2){
    double f0 = gsl_vector_get(point, dim - 2);
    double f1 = gsl_vector_get(point, dim - 1);
    double f2 = gsl_vector_get(point, dim);
    
    double lower = nmin * pow(f1, 2) / f0;
    double upper = nmax * pow(f1, 2) / f0;
    
    double val = resetinsidebounds(f2, lower, upper);
    gsl_vector_set(point, dim, val);
    return;
  }
  else if (dim % 3 == 0){
    printf("Second knot f0 \n");
    double f0n1 = gsl_vector_get(point, dim - 3);
    double f1n1 = gsl_vector_get(point, dim - 2);
    double f2n1 = gsl_vector_get(point, dim - 1);
    double f0 =   gsl_vector_get(point, dim);
    
    double nprev = f2n1 * f0n1 / pow(f1n1, 2);
    double kprev = - f1n1 / pow(f0n1, nprev);
    
    double *nminmax = NMinMax(nprev, ntol, nmin, nmax, segmentlength);
    double *kminmax = KMinMax(kprev, ktol, kmin, kmax, segmentlength);
    
    ///printf("%f, %f, %f \n", nprev, nminmax[0], nminmax[1]);
    
    double lower = F0BoundMinMax(f0n1, nminmax[0], nminmax[1], kminmax[0], kminmax[1], segmentlength, -1);
    double upper = F0BoundMinMax(f0n1, nminmax[1], nminmax[0], kminmax[1], kminmax[0], segmentlength,  1);
    
    double val = resetinsidebounds(f0, lower, upper);
    gsl_vector_set(point, dim, val);
    
  }
  else if (dim % 3 == 1){
    printf("Second knot f1 \n");
    double f0n1 = gsl_vector_get(point, dim - 4);
    double f1n1 = gsl_vector_get(point, dim - 3);
    double f2n1 = gsl_vector_get(point, dim - 2);
    double f0 =   gsl_vector_get(point, dim - 1);
    double f1 =   gsl_vector_get(point, dim);
    
    double nprev = f2n1 * f0n1 / pow(f1n1, 2);
    double kprev = - f1n1 / pow(f0n1, nprev);
    
    double *nminmax = NMinMax(nprev, ntol, nmin, nmax, segmentlength);
    double *kminmax = KMinMax(kprev, ktol, kmin, kmax, segmentlength);
    
    double lower = F1BoundMinMax(f0, f0n1, nminmax[1], nminmax[0], kminmax[1], kminmax[0], segmentlength, -1);
    double upper = F1BoundMinMax(f0, f0n1, nminmax[0], nminmax[1], kminmax[0], kminmax[1], segmentlength,  1);
    printf("Lower upper %E, %E \n", lower, upper);
    double val = resetinsidebounds(f1, lower, upper);
    gsl_vector_set(point, dim, val);
  }
  else if (dim % 3 == 2){
    printf("Second knot f2 \n");
    double f0n1 = gsl_vector_get(point, dim - 5);
    double f1n1 = gsl_vector_get(point, dim - 4);
    double f2n1 = gsl_vector_get(point, dim - 3);
    double f0 =   gsl_vector_get(point, dim - 2);
    double f1 =   gsl_vector_get(point, dim - 1);
    double f2 =   gsl_vector_get(point, dim);
    
    double nprev = f2n1 * f0n1 / pow(f1n1, 2);
    double kprev = - f1n1 / pow(f0n1, nprev);
    
    double *nminmax = NMinMax(nprev, ntol, nmin, nmax, segmentlength);
    double *kminmax = KMinMax(kprev, ktol, kmin, kmax, segmentlength);
    
    double lower = F2BoundMinMax(f0, f0n1, f1, nminmax[0], kminmax[0], segmentlength, -1);
    double upper = F2BoundMinMax(f0, f0n1, f1, nminmax[1], kminmax[1], segmentlength,  1);
    
    double val = resetinsidebounds(f2, lower, upper);
    gsl_vector_set(point, dim, val);
  }
}


static void resetoutofboundspoint(
  gsl_vector* point,
  double fmin,
  double fmax,
  double nmin,
  double nmax,
  double ntol,
  double kmin,
  double kmax,
  double ktol,
  double segmentlength
  )
{
  int dim = point->size;
  /*
  printf("Before \n");
  for (int i = 0; i < dim; ++i){
    printf("%E, ", gsl_vector_get(point, i));
  }
  printf("\n");
  */
  for (int i = 0; i < dim; ++i){
    
    if (dim % 3 == 0){
      printf("Dim being reset is f0 \n");
    }
    else if (dim % 3 == 1){
      printf("Dim being reset is f1 \n");
    }
    else if (dim % 3 == 2){
      printf("Dim being reset is f2 \n");
    }
    
    resetdimonpoint(point, i, fmin, fmax, nmin, nmax, ntol, kmin, kmax, ktol, segmentlength);
  }
  /*
  printf("After \n");
  for (int i = 0; i < dim; ++i){
    printf("%E, ", gsl_vector_get(point, i));
  }
  printf("\n");
  */
}


///
/// Defines the bounds for the first and second derivative frequency parameters on the first knot (t = 0).
///
static double FirstKnotDerivBound(
  const void *data,
  const size_t dim,
  const gsl_matrix *cache UNUSED,
  const gsl_vector* point
  )
{
  const FirstKnotBoundInfo* info = (const FirstKnotBoundInfo*)data;
  
  double fmin = info->fmin;
  double fmax = info->fmax;
  double n = info->n;
  double k = info->k;
  
  double f0 = resetinsidebounds(gsl_vector_get(point, 0), fmin, fmax);
  
  if (dim == 1){
    double f1 = -k * pow(f0, n);
    
    return f1;
  }
  else if (dim == 2){
    
    double lower = -k * pow(fmax, n);
    double upper = -k * pow(fmin, n);
    
    double f1 = resetinsidebounds(gsl_vector_get(point, 1), lower, upper);
 
    double f2 = n * pow(f1, 2) / f0;
    
    
    return f2;
  }

  return NAN;
}


///
/// Sets the bound on the frequency parameter
///
static double F0Bound(
  const void *data,
  const size_t dim,
  const gsl_matrix *cache UNUSED,
  const gsl_vector* pointorig
  )
{
  size_t vectorlength = pointorig->size;
  gsl_vector* point = gsl_vector_alloc(vectorlength);
  memcpy(point, pointorig, sizeof(gsl_vector));
  
  const PiecewiseBoundInfo* info = (const PiecewiseBoundInfo*)data;
  
  double fmin = info->fmin;
  double fmax = info->fmax;
  double nmin = info->nmin;
  double nmax = info->nmax;
  double ntol = info->ntol;
  double kmin = info->kmin;
  double kmax = info->kmax;
  double ktol = info->ktol;
  double segmentlength = info->segmentlength;
  int upperlower = info->upperlower;
  
  resetoutofboundspoint(point, fmin, fmax, nmin, nmax, ntol, kmin, kmax, ktol, segmentlength);
  
  double f0n1 = gsl_vector_get(point, dim - 3);
  double f1n1 = gsl_vector_get(point, dim - 2);
  double f2n1 = gsl_vector_get(point, dim - 1);
  
  double nprev = f2n1 * f0n1 / pow(f1n1, 2);
  double kprev = - f1n1 / pow(f0n1, nprev);
  
  double *nminmax = NMinMax(nprev, ntol, nmin, nmax, segmentlength);
  double *kminmax = KMinMax(kprev, ktol, kmin, kmax, segmentlength);
  
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
/// Sets the bound on the first derivative frequency parameter 
///
static double F1Bound(
  const void *data,
  const size_t dim,
  const gsl_matrix *cache UNUSED,
  const gsl_vector* pointorig
  )
{
  size_t vectorlength = pointorig->size;
  gsl_vector* point = gsl_vector_alloc(vectorlength);
  memcpy(point, pointorig, sizeof(gsl_vector));
  
  const PiecewiseBoundInfo* info = (const PiecewiseBoundInfo*)data;
  
  double fmin = info->fmin;
  double fmax = info->fmax;
  double nmin = info->nmin;
  double nmax = info->nmax;
  double ntol = info->ntol;
  double kmin = info->kmin;
  double kmax = info->kmax;
  double ktol = info->ktol;
  double segmentlength = info->segmentlength;
  int upperlower = info->upperlower;
  
  resetoutofboundspoint(point, fmin, fmax, nmin, nmax, ntol, kmin, kmax, ktol, segmentlength);
  
  double f0n1 = gsl_vector_get(point, dim - 4);
  double f1n1 = gsl_vector_get(point, dim - 3);
  double f2n1 = gsl_vector_get(point, dim - 2);
  
  double f0 = gsl_vector_get(point, dim - 1);
  
  ///printf("%f, %f, %f, %f \n", f0n1, f1n1, f2n1, f0);
  
  double nprev = f2n1 * f0n1 / pow(f1n1, 2);
  double kprev = - f1n1 / pow(f0n1, nprev);
  
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
/// Sets the bounds on the second derivative frequency parameter 
///
static double F2Bound(
  const void *data,
  const size_t dim,
  const gsl_matrix *cache UNUSED,
  const gsl_vector* pointorig
  )
{
  size_t vectorlength = pointorig->size;
  gsl_vector* point = gsl_vector_alloc(vectorlength);
  memcpy(point, pointorig, sizeof(gsl_vector));
  
  const PiecewiseBoundInfo* info = (const PiecewiseBoundInfo*)data;
  
  double fmin = info->fmin;
  double fmax = info->fmax;
  double nmin = info->nmin;
  double nmax = info->nmax;
  double ntol = info->ntol;
  double kmin = info->kmin;
  double kmax = info->kmax;
  double ktol = info->ktol;
  double segmentlength = info->segmentlength;
  int upperlower = info->upperlower;
  
  resetoutofboundspoint(point, fmin, fmax, nmin, nmax, ntol, kmin, kmax, ktol, segmentlength);
  
  double f0n1 = gsl_vector_get(point, dim - 5);
  double f1n1 = gsl_vector_get(point, dim - 4);
  double f2n1 = gsl_vector_get(point, dim - 3);
  
  double f0 = gsl_vector_get(point, dim - 2);
  double f1 = gsl_vector_get(point, dim - 1);
  
  double nprev = f2n1 * f0n1 / pow(f1n1, 2);
  double kprev = - f1n1 / pow(f0n1, nprev);
  
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
  double f0,     /// Frequency at time t = 0. For calculating kmax/kmin from initial tau values, f0 should be maximised. Likewise for calculating taumin/taumax from k values 
  double n,      /// A braking index. For calculating kmax/kmin from initial tau values, n should be maximised. Likewise for calculating taumin/taumax from k values 
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
  const double fmin,   /// Minimum spin frequency
  const double fmax,   /// Maximum spin frequency
  const double nmin,   /// Minimum braking index
  const double nmax,   /// Maximum braking index
  const double ntol,   /// Tolerance (percentage) between braking indices on adjacent knots
  const double taumin, /// Minimum tau value
  const double taumax, /// Maximum tau value
  const double ktol,   /// Tolerance (percentage) between k values on adjacent knots
  const gsl_vector* knots, /// List of knots
  const int finalknot  /// The number of the final knot
  )
{
  /* Simple checks. Not sure what I assume the types of errors are, something to ask about later */
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  
  /// Converting tau values to k values
  double kmin = ktauconversion(fmax, nmax, taumax);
  double kmax = ktauconversion(fmax, nmax, taumin);
  
  
  /// Setting the first knot bounds
  XLALSetLatticeTilingConstantBound(tiling, 0, fmin, fmax);
  FirstKnotBoundInfo XLAL_INIT_DECL( info_first_knot_lower );
  FirstKnotBoundInfo XLAL_INIT_DECL( info_first_knot_upper );
  
  info_first_knot_lower.fmin = info_first_knot_upper.fmin = fmin;
  info_first_knot_lower.fmax = info_first_knot_upper.fmax = fmax;
  
  info_first_knot_lower.n = nmax;
  info_first_knot_lower.k = kmax;
  
  info_first_knot_upper.n = nmin;
  info_first_knot_upper.k = kmin;
  
  XLAL_CHECK(XLALSetLatticeTilingBound(tiling, 1, FirstKnotDerivBound, sizeof( info_first_knot_lower ), &info_first_knot_lower, &info_first_knot_upper) == XLAL_SUCCESS, XLAL_EFAILED);
  XLAL_CHECK(XLALSetLatticeTilingBound(tiling, 2, FirstKnotDerivBound, sizeof( info_first_knot_lower ), &info_first_knot_upper, &info_first_knot_lower) == XLAL_SUCCESS, XLAL_EFAILED);
  
  /// Setting the bounds for all following knots
  
  for (int knot = 1; knot < finalknot; ++knot){
    double segmentlength = gsl_vector_get(knots, knot) - gsl_vector_get(knots, knot - 1);
    
    PiecewiseBoundInfo XLAL_INIT_DECL( info_knot_lower );
    PiecewiseBoundInfo XLAL_INIT_DECL( info_knot_upper );
    
    info_knot_lower.fmin = info_knot_upper.fmin = fmin;
    info_knot_lower.fmax = info_knot_upper.fmax = fmax;
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
    
  }
  
  printf("Bounds set \n");
  printf("\n");
  return XLAL_SUCCESS;
}



/// resetoutofboundspoint(point, fmin, fmax, nmin, nmax, ntol, kmin, kmax, ktol, segmentlength);







