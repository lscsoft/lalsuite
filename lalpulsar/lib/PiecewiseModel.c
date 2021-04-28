
#include <stdio.h>
#include <math.h>
#include <lal/LatticeTiling.h>
#include <lal/PiecewiseModel.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/*
///
/// A method for printing the contents of a gsl vector
///
static void printvector(
  gsl_vector* vec
  )
{
  size_t len = vec->size;
  
  printf("{");
  for (size_t i = 0; i < len; ++i){
    double elem = gsl_vector_get(vec, i);
    printf("%E", elem);
    
    if(i < len - 1){
      printf(", ");
    } else{
      printf("} \n");
    }
  }
}
*/


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
  
  if (x > radius || x < - radius){
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
  //XLALSetLatticeTilingPaddingFlags(tiling, 1, LATTICE_TILING_PAD_NONE);
  
  return XLAL_SUCCESS;
}

///
/// Information required to determine the bounds on the parameters on the first knot
///
typedef struct tagFirstKnotBoundInfo{
  double fmin;  /// Global minimum starting frequency
  double fmax;  /// Global maximum starting frequency
  double nmin;  /// Initial minimum allowed braking
  double nmax;  /// Initial maximum allowed braking index
  double kmin;  /// Initial minimum k value
  double kmax;  /// Initial maximum k value
  int minmax;   /// +1 for upper bound, -1 for lower bound
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
  double t,  /// Time at which to evaluate the GTE
  int d      /// Derivative order (d <= 2)
  )
{
  double base = 1 + (n - 1) * k * pow(f0, n - 1) * t;
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
  
  return NAN;
}

///
/// Calculates the valid braking index for the next knot based upon the previous braking index. Can be defined in any manner. 
///
static void NMinMax(
  double * nminmax,           /// Nminmax array
  double n,                   /// Previous braking index
  double ntol,                /// Braking index tolerance (percentage)
  double nmin,                /// Minimum acceptable braking index
  double nmax,                /// Maximum acceptable braking index
  double segmentlength UNUSED /// Time difference between current knot and previous knot. Kept in case function definition changes at a later date
  )
{
  double nextnmin = n * (1 - ntol);
  double nextnmax = n; //* (1 + ntol);
  
  if (nmin > nmax){
    printf("nmin greater than nmax, %E, %E \n", nmin, nmax);
    nminmax[0] = NAN;
    nminmax[1] = NAN;
  }
  
  if (nextnmin < nmin){
    nextnmin = nmin;
  }
  if (nextnmax > nmax){
    nextnmax = nmax;
  }
  
  nminmax[0] = nextnmin;
  nminmax[1] = nextnmax;
  
  if (nextnmin > nextnmax){
    printf("Ranges from NMinMax are incorrect %E, %E \n", nextnmin, nextnmax);
    
    nminmax[0] = NAN;
    nminmax[1] = NAN;
  }
}


///
/// Calculates the valid k value for the next knot based upon the previous value. Can be defined in any manner. 
///
static void KMinMax(
  double * kminmax,           /// kminmax array
  double k,                   /// Previous braking index
  double ktol,                /// Braking index tolerance (percentage)
  double kmin,                /// Minimum acceptable braking index
  double kmax,                /// Maximum acceptable braking index
  double segmentlength UNUSED /// Time difference between current knot and previous knot. Kept in case function definition changes at a later date
  )
{
  double nextkmin = k * (1 - ktol);
  double nextkmax = k; //* (1 + ktol);
  
  if (nextkmax < kmin || nextkmin > kmax){
    printf("Calculated ranges outside of global range, %E, %E, %E, %E, %E \n", k, kmin, kmax, k - kmin, k - kmax);
  }
  
  if (kmin > kmax){
    printf("kmin greater than kmax, %E, %E \n", kmin, kmax);
    kminmax[0] = NAN;
    kminmax[1] = NAN;
  }
  
  
  if (nextkmin < kmin){
    nextkmin = kmin;
  }
  if (nextkmax > kmax){
    nextkmax = kmax;
  }
  
  kminmax[0] = nextkmin;
  kminmax[1] = nextkmax;
  
  if (nextkmin > nextkmax){
    printf("Ranges from KMinMax are incorrect %E, %E \n", nextkmin, nextkmax);
    kminmax[0] = NAN;
    kminmax[1] = NAN;
  }
}


///
/// Calculates the minimum and maximum bounds for a frequency parameter
///
static double F0BoundMinMax(
  double f0,        /// Frequency value of the previous knot
  double na,        /// A braking index value. For calculating upper bound, na > nb. For calculating lower bound na < nb
  double nb,        /// A braking index value. For calculating upper bound, na > nb. For calculating lower bound na < nb
  double ka UNUSED, /// A k value. For calculating upper bound, ka > kb. For calculating lower bound ka < kb
  double kb,        /// A k value. For calculating upper bound, ka > kb. For calculating lower bound ka < kb
  double seglength, /// Time difference between this knot and the previous knot
  int minmax        /// +1 for calculating upper bound, -1 for calculating lower bound
)
{
  /// These combinations of gted and gtedd give one of the optimised values of f0 (either minimum or maximum) such that the braking index at the next knot still
  /// falls within an acceptable range. This acceptable range is between [na, nb] or [nb, na], depending upon the value of minmax. Do a double check of this maths
  if (na > nb && minmax == -1){
    printf("F0BoundMinMax being called incorrectly, with %f, %f, %d \n", na, nb, minmax);
    return NAN;
  }
  else if (na < nb && minmax == 1){
    printf("F0BoundMinMax being called incorrectly, with %f, %f, %d \n", na, nb, minmax);
    return NAN;
  }
  
  /// These variables are used to calculate the f0 value using braking index criteria
  double gted = GTEAndDerivs(f0, na, ka, seglength, 1);
  double gtedd = GTEAndDerivs(f0, nb, kb, seglength, 2);
  double bindexcriteria = na * pow(gted, 2) / gtedd;
  
  /// Parameter range optimisation for frequency parameter
  double paramrange = GTEAndDerivs(f0, nb, kb, seglength, 0);
  
  gsl_vector* vals = gsl_vector_alloc(2);
  gsl_vector_set(vals, 0, bindexcriteria);
  gsl_vector_set(vals, 1, paramrange);
  
  if(minmax == 1){
    double max = gsl_vector_min(vals);
    return max;
  }
  else if (minmax == -1){
    double min = gsl_vector_max(vals);
    return min;
  }
  
  return NAN;
}


///
/// Function which determines the minimum/maximum value of the first derivative parameters
///
static double F1BoundMinMax(
  double f0 UNUSED,        /// Frequency parameter at the corresponding knot
  double fprev,     /// Frequency parameter at the previous knot
  double na,        /// A braking index value. For calculating upper bound, na < nb. For calculating lower bound na > nb
  double nb,        /// A braking index value. For calculating upper bound, na < nb. For calculating lower bound na > nb
  double ka,        /// A k value. For calculating upper bound, ka < kb. For calculating lower bound ka > kb
  double kb UNUSED, /// A k value. For calculating upper bound, ka < kb. For calculating lower bound ka > kb
  double seglength, /// Time difference between this knot and the previous
  int minmax        /// +1 for upper bound, -1 for lower bound
  )
{
  /// Checking that the method is being called with the appropriate values
  if (na < nb && minmax == -1){
    printf("F1BoundMinMax being called incorrectly, with %f, %f, %d \n", na, nb, minmax);
    return NAN;
  }
  if (na > nb && minmax == 1){
    printf("F1BoundMinMax being called incorrectly, with %f, %f, %d \n", na, nb, minmax);
    return NAN;
  }
  
  /// The two parameter range conditions
  double prangecondition1 = -ka * pow(f0, na);
  double prangecondition2 = GTEAndDerivs(fprev, na, ka, seglength, 1);
  
  gsl_vector* vals = gsl_vector_alloc(2);
  gsl_vector_set(vals, 0, prangecondition1);
  gsl_vector_set(vals, 1, prangecondition2);
  
  if(minmax == 1){
    double max = gsl_vector_min(vals);
    return max;
  }
  else if(minmax == -1){
    double min = gsl_vector_max(vals);
    return min;
  }
  
  return NAN; 
}

///
/// Calculates the bound for the second derivative frequency parameter 
///
static double F2BoundMinMax(
  double f0,        /// Frequency parameter on this knot
  double fprev,     /// Frequency parameter on the previous knot
  double f1,        /// First derivative frequency parameter on this knot
  double n,         /// Optimal braking index value. For upper bound n should be maximised, for lower bound n should be minimised
  double ka,        /// A k value. For calculating upper bound, ka > kb. For calculating lower bound ka < kb
  double kb,        /// A k value. For calculating upper bound, ka > kb. For calculating lower bound ka < kb 
  double seglength, /// Time difference between this knot and the previous
  int minmax        /// +1 for upper bound, -1 for lower bound
)
{
  /// Checking that the method is being called with the appropriate values
  if (ka > kb && minmax == -1){
    printf("F2BoundMinMax being called incorrectly, with %f, %f, %d \n", ka, kb, minmax);
    return NAN;
  }
  if (ka < kb && minmax == 1){
    printf("F2BoundMinMax being called incorrectly, with %f, %f, %d \n", ka, kb, minmax);
    return NAN;
  }
  
  /// Braking index and parameter range conditions.
  double bindexcriteria = n * pow(f1, 2) / f0;
  double prangecriteria = GTEAndDerivs(fprev, n, ka, seglength, 2);
  double kcriteria = - pow(f1, 2) * log(- kb / f1) / (f0 * log(f0));
  
  gsl_vector* vals = gsl_vector_alloc(3);
  gsl_vector_set(vals, 0, bindexcriteria);
  gsl_vector_set(vals, 1, prangecriteria);
  gsl_vector_set(vals, 2, kcriteria);
  
  if(minmax == 1){
    double max = gsl_vector_min(vals);
    return max;
  }
  else if(minmax == -1){
    double min = gsl_vector_max(vals);
    return min;
  }
  
  return NAN;
}

///
/// Takes an input val, and returns val to its closest value within the range [valmin, valmax].
///
static double resetinsidebounds(
  double val,     /// Value
  double valmin,  /// Lower bound of range
  double valmax   /// Upper bound of range
  )
{
  
  if (valmin <= val && val <= valmax){
    return val;
  }
  else if (valmin > valmax){
  
    if ((valmin - valmax) / valmax < pow(10, -10)){
      return valmax;
    }
    printf("Oh no, valmin is bigger than valmax! \n");
    printf("val, valmin and valmax: %E, %E, %E, %E \n", val, valmin, valmax, valmax - valmin);
    return NAN;
  }
  else if (val < valmin){
    return valmin;
  }
  else if (val > valmax){
    return valmax;
  }
  
  return NAN;
}


///
/// A method which will reset a given point to be within our parameter sapce
///
static void resetdimonpoint(
  gsl_vector* point,   /// The point which we are resetting to be within the parameter space
  int dim,             /// The dimension which we are resetting
  double fmin,         /// The global minimum frequency
  double fmax,         /// The global maximum frequency
  double nmin,         /// The minimum braking index allowed for the knot belonging to dim
  double nmax,         /// The maximum braking index allowed for the knot belonging to dim
  double ntol,         /// The percentage tolerance we allow for the braking index on adjacent knots
  double kmin,         /// The minimum k value allowed for the knot belonging to dim
  double kmax,         /// The maximum k value allowed for the knot belonging to dim
  double ktol,         /// The percentage tolerance we allow for the k value on adjacent knots
  double segmentlength /// The length of the segment we are working on
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
    
    /// We need a special case when f1 is equal to one of its bounds. The reason being a floating point error arises which results in calculated k values being
    /// off by only a small percentage. This is explained in further comments below.
    if (f1 == -kmax * pow(f0, nmax)){
      double val = nmax * pow(f1, 2) / f0;
      gsl_vector_set(point, dim, val);
      return;
    }
    else if (f1 == -kmin * pow(f0, nmin)){
      double val = nmin * pow(f1, 2) / f0;
      gsl_vector_set(point, dim, val);
      return;
    }
    
    double lower = F2BoundMinMax(f0, 0.0000000001, f1, nmin, kmin, kmax, 0, -1);
    double upper = F2BoundMinMax(f0, 100000000000, f1, nmax, kmax, kmin, 0, 1);
    
    double val = resetinsidebounds(f2, lower, upper);
    gsl_vector_set(point, dim, val);
    return;
    
  }
  else if (dim % 3 == 0){
  
    double f0n1 = gsl_vector_get(point, dim - 3);
    double f1n1 = gsl_vector_get(point, dim - 2);
    double f2n1 = gsl_vector_get(point, dim - 1);
    double f0 =   gsl_vector_get(point, dim);
    
    double nprev = f2n1 * f0n1 / pow(f1n1, 2);
    double kprev = - f1n1 / pow(f0n1, nprev);
    
    /// In the instant that f1n1 has a value at one of its extrema, the upper and lower bound on f2 are the same point (i.e. there is only one allowed value for f1)
    /// However, given the equations for f2 (either starting from braking index requirements or k value requirements), small floating point error arises between the 
    /// two end results. For this reason, we must choose one equation to use to calculate f2 (in this case we have chosen the braking index equation) and then manually
    /// edit the value of k and n we decide to use for our checks on subsequent knots from that point. This is the reason for the below conditions. The typically error
    /// between kprev and either kmax and kmin is on the order of ~10^-23, likewise for nprev and nmin/nmax the error is of the order ~10^-16.
    if (kprev > kmax && kprev * 0.99999 <= kmax){
      kprev = kmax;
    }
    if (kprev < kmin && kprev * 1.00001 >= kmin){
      kprev = kmin;
    }
    
    if (nprev > nmax && nprev * 0.99999 <= nmax){
      nprev = nmax;
    }
    if (nprev < nmin && nprev * 1.00001 >= nmin){
      nprev = nmin;
    }
    
    double nminmax[2];
    double kminmax[2];
    
    NMinMax(nminmax, nprev, ntol, nmin, nmax, segmentlength);
    KMinMax(kminmax, kprev, ktol, kmin, kmax, segmentlength);
    
    double lower = F0BoundMinMax(f0n1, nminmax[0], nminmax[1], kminmax[0], kminmax[1], segmentlength, -1);
    double upper = F0BoundMinMax(f0n1, nminmax[1], nminmax[0], kminmax[1], kminmax[0], segmentlength,  1);
    
    double val = resetinsidebounds(f0, lower, upper);
    gsl_vector_set(point, dim, val);
    
  }
  else if (dim % 3 == 1){
  
    double f0n1 = gsl_vector_get(point, dim - 4);
    double f1n1 = gsl_vector_get(point, dim - 3);
    double f2n1 = gsl_vector_get(point, dim - 2);
    double f0 =   gsl_vector_get(point, dim - 1);
    double f1 =   gsl_vector_get(point, dim);
    
    double nprev = f2n1 * f0n1 / pow(f1n1, 2);
    double kprev = - f1n1 / pow(f0n1, nprev);
    
    /// In the instant that f1n1 has a value at one of its extrema, the upper and lower bound on f2 are the same point (i.e. there is only one allowed value for f1)
    /// However, given the equations for f2 (either starting from braking index requirements or k value requirements), small floating point error arises between the 
    /// two end results. For this reason, we must choose one equation to use to calculate f2 (in this case we have chosen the braking index equation) and then manually
    /// edit the value of k and n we decide to use for our checks on subsequent knots from that point. This is the reason for the below conditions. The typically error
    /// between kprev and either kmax and kmin is on the order of ~10^-23, likewise for nprev and nmin/nmax the error is of the order ~10^-16.
    if (kprev > kmax && kprev * 0.99999 <= kmax){
      kprev = kmax;
    }
    if (kprev < kmin && kprev * 1.00001 >= kmin){
      kprev = kmin;
    }
    
    if (nprev > nmax && nprev * 0.99909 <= nmax){
      nprev = nmax;
    }
    if (nprev < nmin && nprev * 1.00001 >= nmin){
      nprev = nmin;
    }
    
    double nminmax[2];
    double kminmax[2];
    
    NMinMax(nminmax, nprev, ntol, nmin, nmax, segmentlength);
    KMinMax(kminmax, kprev, ktol, kmin, kmax, segmentlength);
    
    double lower = F1BoundMinMax(f0, f0n1, nminmax[1], nminmax[0], kminmax[1], kminmax[0], segmentlength, -1);
    double upper = F1BoundMinMax(f0, f0n1, nminmax[0], nminmax[1], kminmax[0], kminmax[1], segmentlength,  1);
    
    double val = resetinsidebounds(f1, lower, upper);
    gsl_vector_set(point, dim, val);
  }
  else if (dim % 3 == 2){
    
    double f0n1 = gsl_vector_get(point, dim - 5);
    double f1n1 = gsl_vector_get(point, dim - 4);
    double f2n1 = gsl_vector_get(point, dim - 3);
    double f0 =   gsl_vector_get(point, dim - 2);
    double f1 =   gsl_vector_get(point, dim - 1);
    double f2 =   gsl_vector_get(point, dim);
    
    double nprev = f2n1 * f0n1 / pow(f1n1, 2);
    double kprev = - f1n1 / pow(f0n1, nprev);
    
    /// In the instant that f1n1 has a value at one of its extrema, the upper and lower bound on f2 are the same point (i.e. there is only one allowed value for f1)
    /// However, given the equations for f2 (either starting from braking index requirements or k value requirements), small floating point error arises between the 
    /// two end results. For this reason, we must choose one equation to use to calculate f2 (in this case we have chosen the braking index equation) and then manually
    /// edit the value of k and n we decide to use for our checks on subsequent knots from that point. This is the reason for the below conditions. The typically error
    /// between kprev and either kmax and kmin is on the order of ~10^-23, likewise for nprev and nmin/nmax the error is of the order ~10^-16.
    if (kprev > kmax && kprev * 0.99999 <= kmax){
      kprev = kmax;
    }
    if (kprev < kmin && kprev * 1.00001 >= kmin){
      kprev = kmin;
    }
    
    if (nprev > nmax && nprev * 0.99999 <= nmax){
      nprev = nmax;
    }
    if (nprev < nmin && nprev * 1.00001 >= nmin){
      nprev = nmin;
    }
    
    double nminmax[2];
    double kminmax[2];
    
    NMinMax(nminmax, nprev, ntol, nmin, nmax, segmentlength);
    KMinMax(kminmax, kprev, ktol, kmin, kmax, segmentlength);
    
    double lower = F2BoundMinMax(f0, f0n1, f1, nminmax[0], kminmax[0], kminmax[1], segmentlength, -1);
    double upper = F2BoundMinMax(f0, f0n1, f1, nminmax[1], kminmax[1], kminmax[0], segmentlength,  1);
    
    double val = resetinsidebounds(f2, lower, upper);
    gsl_vector_set(point, dim, val);
  }
}

///
/// Resets a point such that it is within the bounds of our parameter space
///
static void resetoutofboundspoint(
  gsl_vector* point,   /// The point which we are resetting to be within the parameter space
  double fmin,         /// The global minimum frequency
  double fmax,         /// The global maximum frequency
  double nmin,         /// The minimum braking index allowed for the knot belonging to dim
  double nmax,         /// The maximum braking index allowed for the knot belonging to dim
  double ntol,         /// The percentage tolerance we allow for the braking index on adjacent knots
  double kmin,         /// The minimum k value allowed for the knot belonging to dim
  double kmax,         /// The maximum k value allowed for the knot belonging to dim
  double ktol,         /// The percentage tolerance we allow for the k value on adjacent knots
  double segmentlength /// The length of the segment we are working on
  )
{
  int dim = point->size;
  
  for (int i = 0; i < dim; ++i){
    resetdimonpoint(point, i, fmin, fmax, nmin, nmax, ntol, kmin, kmax, ktol, segmentlength);
  }
}


///
/// Defines the bounds for the first and second derivative frequency parameters on the first knot (t = 0).
///
static double FirstKnotDerivBound(
  const void *data,
  const size_t dim,
  const gsl_matrix *cache UNUSED,
  const gsl_vector* pointorig
  )
{
  size_t vectorlength = pointorig->size;
  gsl_vector* point = gsl_vector_alloc(vectorlength);
  memcpy(point, pointorig, sizeof(gsl_vector));
  
  const FirstKnotBoundInfo* info = (const FirstKnotBoundInfo*)data;
  
  double fmin = info->fmin;
  double fmax = info->fmax;
  double nmin = info->nmin;
  double nmax = info->nmax;
  double kmin = info->kmin;
  double kmax = info->kmax;
  int minmax = info->minmax;
  
  double f0 = resetinsidebounds(gsl_vector_get(point, 0), fmin, fmax);
  gsl_vector_set(point, 0, f0);
  
  if (dim == 1){
    
    if (minmax == 1){
      double f1 = -kmin * pow(f0, nmin);
      return f1;
    }
    else if (minmax == -1){
      double f1 = -kmax * pow(f0, nmax);
      return f1;
    }
  }
  else if (dim == 2){
    
    double lower = -kmax * pow(f0, nmax);
    double upper = -kmin * pow(f0, nmin);
    
    double f1 = resetinsidebounds(gsl_vector_get(point, 1), lower, upper);
    gsl_vector_set(point, 1, f1);
    
    if (minmax == 1){
      double f2 = F2BoundMinMax(f0, 100000000000, f1, nmax, kmax, kmin, 0, 1);
      return f2;
    }
    else if (minmax == -1){
      double f2 = F2BoundMinMax(f0, 0.0000000001, f1, nmin, kmin, kmax, 0, -1);
      return f2;
    }
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
  
  if (kprev > kmax && kprev * 0.99999 <= kmax){
      kprev = kmax;
    }
  if (kprev < kmin && kprev * 1.00001 >= kmin){
    kprev = kmin;
  }
    
  if (nprev > nmax && nprev * 0.99999 <= nmax){
    nprev = nmax;
  }
  if (nprev < nmin && nprev * 1.00001 >= nmin){
    nprev = nmin;
  }
  
  double nminmax[2];
  double kminmax[2];
  
  NMinMax(nminmax, nprev, ntol, nmin, nmax, segmentlength);
  KMinMax(kminmax, kprev, ktol, kmin, kmax, segmentlength);
  
  if (upperlower == 1){
    double upperbound = F0BoundMinMax(f0n1, nminmax[1], nminmax[0], kminmax[1], kminmax[0], segmentlength, 1);
    return upperbound;
  }
  else if (upperlower == -1){
    double lowerbound = F0BoundMinMax(f0n1, nminmax[0], nminmax[1], kminmax[0], kminmax[1], segmentlength, -1);
    return lowerbound;
  }
  
  return NAN;
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
  
  double nprev = f2n1 * f0n1 / pow(f1n1, 2);
  double kprev = - f1n1 / pow(f0n1, nprev);
  
  if (kprev > kmax && kprev * 0.99999 <= kmax){
      kprev = kmax;
    }
  if (kprev < kmin && kprev * 1.00001 >= kmin){
    kprev = kmin;
  }
    
  if (nprev > nmax && nprev * 0.99999 <= nmax){
    nprev = nmax;
  }
  if (nprev < nmin && nprev * 1.00001 >= nmin){
    nprev = nmin;
  }
  
  double nminmax[2];
  double kminmax[2];
  
  NMinMax(nminmax, nprev, ntol, nmin, nmax, segmentlength);
  KMinMax(kminmax, kprev, ktol, kmin, kmax, segmentlength);
  
  if (upperlower == 1){
    double upperbound = F1BoundMinMax(f0, f0n1, nminmax[0], nminmax[1], kminmax[0], kminmax[1], segmentlength, 1);
    return upperbound;
  }
  else if (upperlower == -1){
    double lowerbound = F1BoundMinMax(f0, f0n1, nminmax[1], nminmax[0], kminmax[1], kminmax[0], segmentlength, -1);
    return lowerbound;
  }
  
  return NAN;
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
  
  if (kprev > kmax && kprev * 0.99999 <= kmax){
      kprev = kmax;
    }
  if (kprev < kmin && kprev * 1.00001 >= kmin){
    kprev = kmin;
  }
    
  if (nprev > nmax && nprev * 0.99999 <= nmax){
    nprev = nmax;
  }
  if (nprev < nmin && nprev * 1.00001 >= nmin){
    nprev = nmin;
  }
  
  double nminmax[2];
  double kminmax[2];
  
  NMinMax(nminmax, nprev, ntol, nmin, nmax, segmentlength);
  KMinMax(kminmax, kprev, ktol, kmin, kmax, segmentlength);
  
  if (upperlower == 1){
    double upperbound = F2BoundMinMax(f0, f0n1, f1, nminmax[1], kminmax[1], kminmax[0], segmentlength, 1);
    return upperbound;
  }
  else if (upperlower == -1){
    double lowerbound = F2BoundMinMax(f0, f0n1, f1, nminmax[0], kminmax[0], kminmax[1], segmentlength, -1);
    return lowerbound;
  }
  
  return NAN;
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
  /// Converting tau values to k values
  double kmin = ktauconversion(fmax, nmax, taumax);
  double kmax = ktauconversion(fmax, nmax, taumin);
  
  printf("kmin and kmax %E, %E \n", kmin, kmax);
  
  /// Simple checks. Not sure what the types of errors are, something to ask about later
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(fmin < fmax, XLAL_EFAULT);
  XLAL_CHECK(nmin < nmax, XLAL_EFAULT);
  XLAL_CHECK(taumin < taumax, XLAL_EFAULT);
  XLAL_CHECK(kmin < kmax, XLAL_EFAULT);
  
  /// Setting the first knot bounds
  XLALSetLatticeTilingConstantBound(tiling, 0, fmin, fmax);
  
  FirstKnotBoundInfo XLAL_INIT_DECL( info_first_knot_lower );
  FirstKnotBoundInfo XLAL_INIT_DECL( info_first_knot_upper );
  
  info_first_knot_lower.fmin = info_first_knot_upper.fmin = fmin;
  info_first_knot_lower.fmax = info_first_knot_upper.fmax = fmax;
  info_first_knot_lower.nmin = info_first_knot_upper.nmin = nmin;
  info_first_knot_lower.nmax = info_first_knot_upper.nmax = nmax;
  info_first_knot_lower.kmin = info_first_knot_upper.kmin = kmin;
  info_first_knot_lower.kmax = info_first_knot_upper.kmax = kmax;
  
  info_first_knot_lower.minmax = -1;
  info_first_knot_upper.minmax = 1;
  
  XLAL_CHECK(XLALSetLatticeTilingBound(tiling, 1, FirstKnotDerivBound, sizeof( info_first_knot_lower ), &info_first_knot_lower, &info_first_knot_upper) == XLAL_SUCCESS, XLAL_EFAILED);
  XLAL_CHECK(XLALSetLatticeTilingBound(tiling, 2, FirstKnotDerivBound, sizeof( info_first_knot_lower ), &info_first_knot_upper, &info_first_knot_lower) == XLAL_SUCCESS, XLAL_EFAILED);
  
  /// Padding flags options. For reference.
  
  //LATTICE_TILING_PAD_NONE  = 0x00,      ///< Do not add padding, and generate points strictly within parameter space
  //LATTICE_TILING_PAD_LHBBX = 0x01,      ///< Add half-bounding-box padding to lower physical parameter-space bounds
  //LATTICE_TILING_PAD_UHBBX = 0x02,      ///< Add half-bounding-box padding to upper physical parameter-space bounds
  //LATTICE_TILING_PAD_LINTP = 0x04,      ///< Add integer point padding to lower integer parameter-space bounds
  //LATTICE_TILING_PAD_UINTP = 0x08,      ///< Add integer point padding to upper integer parameter-space bounds
  //LATTICE_TILING_PAD_MAX   = 0x20,
  
  LatticeTilingPaddingFlags flags = LATTICE_TILING_PAD_NONE;
  
  XLALSetLatticeTilingPaddingFlags(tiling, 0, flags);
  XLALSetLatticeTilingPaddingFlags(tiling, 1, flags);
  XLALSetLatticeTilingPaddingFlags(tiling, 2, flags);
  
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
    info_knot_lower.kmax = info_knot_upper.kmax = kmax;
    info_knot_lower.ktol = info_knot_upper.ktol = ktol;
    info_knot_lower.segmentlength = info_knot_upper.segmentlength = segmentlength;
    
    info_knot_lower.upperlower = -1;
    info_knot_upper.upperlower = 1;
    
    int dimindex = 3 * knot;
    
    XLAL_CHECK(XLALSetLatticeTilingBound(tiling, dimindex,     F0Bound, sizeof( info_knot_lower ), &info_knot_lower, &info_knot_upper) == XLAL_SUCCESS, XLAL_EFAILED);
    XLAL_CHECK(XLALSetLatticeTilingBound(tiling, dimindex + 1, F1Bound, sizeof( info_knot_lower ), &info_knot_lower, &info_knot_upper) == XLAL_SUCCESS, XLAL_EFAILED);
    XLAL_CHECK(XLALSetLatticeTilingBound(tiling, dimindex + 2, F2Bound, sizeof( info_knot_lower ), &info_knot_lower, &info_knot_upper) == XLAL_SUCCESS, XLAL_EFAILED);
    
    XLALSetLatticeTilingPaddingFlags(tiling, dimindex,     flags);
    XLALSetLatticeTilingPaddingFlags(tiling, dimindex + 1, flags);
    XLALSetLatticeTilingPaddingFlags(tiling, dimindex + 2, flags);
  }
  
  return XLAL_SUCCESS;
}


