
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
/// A method for printing the contents of a gsl vector
static void printvector(
  gsl_vector* vec
  )
{
  size_t len = vec->size;
  
  printf("{");
  for (size_t i = 0; i < len; ++i){
    double elem = gsl_vector_get(vec, i);
    printf("%f", elem);
    
    if(i < len - 1){
      printf(", ");
    } else{
      printf("} \n");
    }
  }
}



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

///
/// Information required to determine the bounds on for the parameters on the first knot
///
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
static double * NMinMax(
  double n,                   /// Previous braking index
  double ntol,                /// Braking index tolerance (percentage)
  double nmin,                /// Minimum acceptable braking index
  double nmax,                /// Maximum acceptable braking index
  double segmentlength UNUSED /// Time difference between current knot and previous knot. Kept in case function definition changes at a later date
  )
{
  static double nminmax[2];
  
  double nextnmin = n * (1 - ntol);
  double nextnmax = n * (1 + ntol);
  
  if (nmin > nmax){
    printf("nmin greater than nmax, %E, %E \n", nmin, nmax);
    nminmax[0] = NAN;
    nminmax[1] = NAN;
    return nminmax;
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
    return nminmax;
  }
  
  return nminmax;
}


///
/// Calculates the valid k value for the next knot based upon the previous value. Can be defined in any manner. 
///
static double * KMinMax(
  double k,                   /// Previous braking index
  double ktol,                /// Braking index tolerance (percentage)
  double kmin,                /// Minimum acceptable braking index
  double kmax,                /// Maximum acceptable braking index
  double segmentlength UNUSED /// Time difference between current knot and previous knot. Kept in case function definition changes at a later date
  )
{
  static double kminmax[2];
  double nextkmin = k * (1 - ktol);
  double nextkmax = k * (1 + ktol);
  
  //printf("Using KMinMax %E, %E, %E, %E, %E \n", k, kmin, kmax, nextkmin, nextkmax);
  
  if (nextkmax < kmin || nextkmin > kmax){
    printf("Calculated ranges outside of global range \n");
  }
  
  if (kmin > kmax){
    printf("kmin greater than kmax, %E, %E \n", kmin, kmax);
    kminmax[0] = NAN;
    kminmax[1] = NAN;
    return kminmax;
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
    return kminmax;
  }
  
  return kminmax;
}



///
/// Calculates the minimum and maximum bounds for a frequency frequency parameter
///
static double F0BoundMinMax(
  double f0,        /// Frequency value of the previous knot
  double na,        /// A braking index value. For calculating upper bound, na > nb. For calculating lower bound na < nb
  double nb,        /// A braking index value. For calculating upper bound, na > nb. For calculating lower bound na < nb
  double ka UNUSED, /// A k value. For calculating upper bound, ka > kb. For calculating lower bound ka < kb
  double kb,        /// A k value. For calculating upper bound, ka > kb. For calculating lower bound ka < kb
  double seglength, /// Time difference between this knot and the previous knot
  int minmax        /// +1 for calculating the maximum bound, -1 for calculating the minimum bound
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
  
  ///printf("These are the f0 options %E, %E \n", brakingindex, paramrange);
  ///printf("Min max is %d \n", minmax);
  
  if (minmax == 1) {
    if (bindexcriteria <= paramrange){
    
      return bindexcriteria;
    }
    else if (bindexcriteria > paramrange){
      return paramrange;
    }
  }
  else if (minmax == -1){
    if (bindexcriteria <= paramrange){
      return paramrange;
    }
    else if (paramrange < bindexcriteria){
      return bindexcriteria;
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
  double kb UNUSED, /// A k value. For calculating upper bound, ka < kb. For calculating lower bound ka > kb
  double seglength, /// Time difference between this knot and the previous
  int minmax        /// +1 for maximum bound, -1 for minimum bound
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
  
  /// The braking index criteria
  double bindexcriteria = -sqrt(GTEAndDerivs(fprev, na, ka, seglength, 2) * f0 / nb);
  double bound;
  
  if (minmax == 1){
    if (prangecondition1 > prangecondition2){
      bound = prangecondition2;
    }
    else {
      bound = prangecondition1;
    }
    if (bound > bindexcriteria){
      return bindexcriteria;
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
    if (bound < bindexcriteria){
      return bindexcriteria;
    }
    else {
      return bound;
    }
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
  double k,         /// Optimal k value. For upper bound k should be maximised, for lower bound k should be minimised
  double seglength, /// Time difference between this knot and the previous
  int minmax        /// +1 for maximum value, -1 for minimum value
)
{
  /// Braking index and parameter range conditions.
  double bindexcriteria = n * pow(f1, 2) / f0;
  double prangecriteria = GTEAndDerivs(fprev, n, k, seglength, 2);
  
  if (minmax == 1){
    if (bindexcriteria > prangecriteria) {
      return prangecriteria;
    }
    else {
      return bindexcriteria;
    }
  }
  else if (minmax == -1) {
    if (bindexcriteria < prangecriteria) {
      return prangecriteria;
    }
    else {
      return bindexcriteria;
    }
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
  
  if (valmin > valmax){
    printf("Oh no, valmin is bigger than valmax! \n");
    printf("val, valmin and valmax: %E, %E, %E \n", val, valmin, valmax);
    return NAN;
  }
  
  if (valmin <= val && val <= valmax){
    return val;
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
/// The Lambert W function (the same as the ProductLog function in Mathematica). Taken from https://stackoverflow.com/questions/60211021/lambert-w-function-in-c-sharp.
/// It may be that other better numerical estimates for the Lambert W function exist. For now will work with this approximation.
///
static double LambertW(
  double x
  )
{
  // LambertW is not defined in this case
  if (x < -exp(-1)){
    printf("Lambert W function not defined for x = %f \n", x);
    return NAN;
  }

  // computes the first branch for real values only

  // amount of iterations (empirically found)
  int iterations;
  if (10 > (int) ceil(log10(x) / 3)){
    iterations = 10;
  } else {
    iterations = (int) ceil(log10(x) / 3);
  }
  // initial guess is based on 0 < ln(a) < 3
  double w = 3 * log(x + 1) / 4;

  // Halley's method via eqn (5.9) in Corless et al (1996). (Obviously there are other methods besides Halley's method, this however seems the simplest to implement for my purposes).
  for (int i = 0; i < iterations; i++)
    w = w - (w * exp(w) - x) / (exp(w) * (w + 1) - (w + 2) * (w * exp(w) - x) / (2 * w + 2));
    
  //printf("Lambert W function results are for x = %f are %f \n", x, w);

  return w;
}

///
/// Due to some interdependence on the conditions between the f1 and f2 parameters on the first knot, we need to define a special method for how these parameters can be 'reset' if 
/// a point is selected outside the parameter space. This method calculates the optimised value of f1 allowed given the calculate k values on that knot. If they are unacceptable, the f1
/// value is reset appropriately using the solution to f1 = - thisk f0 ^ thisn, where thisk = - f1/f0^thisn and thisn = f2 * f0 / f1^2.
///
static void resetfirstdim(
  gsl_vector* point,    /// The vector of which we are resetting the variables f1 and f2 on
  double kmin,          /// The minimum allowed value of k for the first knot
  double kmax,          /// The maximum allowed value of k for the first knot
  double tol            /// How close we desire the current value of k for this knot to be within the range [kmin, kmax]. Defined as a percentage, e.g. if tol = 0.1, we expect k for this knot to be no more than 10% lower or greater than kmin or kmax respectively.
  )
{
  double f0 = gsl_vector_get(point, 0);
  double f1 = gsl_vector_get(point, 1);
  double f2 = gsl_vector_get(point, 2);
  
  double n = f2 * f0 / pow(f1, 2);
  double thisk = - f1 / pow(f0, n);
  
  //printf("Resetting first dim. %f, %f, %f \n", f0, f1, f2);
  if (kmin <= thisk && thisk <= kmax){
    //printf("It's all good brother \n");
    return;
  }
  /// Solutions for the appropriate bounds on f1 were calculated using mathematica and the inequality that kmin < - f1 / f0^n < kmax, where the value of n used is calculated using the values
  /// on that specific knot.
  else if (thisk < kmin){
    //printf("Not all good, thisk < kmin, using LambertW function\n");
    double numerator = 2 * f0 * f2 * log(f0);
    //printf("Lamb W elems %f, %f, %E \n", f0, f2, kmin);
    double denominator = LambertW(numerator / pow(kmin * (1 + tol), 2));
    
    double newf1 = - sqrt(numerator / denominator);
    gsl_vector_set(point, 1, newf1);
    return;
  }
  else if (kmax < thisk){
    //printf("Not all good, kmax < thisk using LambertW function \n");
    double numerator = 2 * f0 * f2 * log(f0);
    //printf("Lamb W elems %f, %f, %E \n", f0, f2, kmin);
    double denominator = LambertW(numerator / pow(kmax * (1 - tol), 2));
    
    double newf1 = - sqrt(numerator / denominator);
    gsl_vector_set(point, 1, newf1);
    return;
  }
  
  printf("Something strange happened, no cases caught, not using LambertW function \n");
  //printf("K values here are %E, %E, %E \n", thisk, kmin, kmax);
  printvector(point);
  if (!(f2 > 0)){
    printf("f2 has a weird value, %E \n", f2);
  }
  if (f2 == 0.){
    printf("f2 is 0! \n");
  }
  
  return; 
}

///
/// As the above method but now for the parameter f2 on the first knot. This value is instead reset however using the expression we have for the braking index, i.e. f2 = n f2^2 / f0
///
static void resetseconddim(
  gsl_vector* point, /// The vector of which we are resetting the variables f1 and f2 on
  double nmin,       /// The minimum allowed value of n for the first knot
  double nmax,       /// The maximum allowed value of n for the first knot
  double kmin,       /// The minimum allowed value of k for the first knot
  double kmax,       /// The maximum allowed value of k for the first knot
  double tol         /// The same tol as in the above method, however now for both n and k criteria
  )
{
  double f0 = gsl_vector_get(point, 0);
  double f1 = gsl_vector_get(point, 1);
  double f2 = gsl_vector_get(point, 2);
  
  if (f2 < 0){
    //printf("f2 below 0, resetting using GTE. Second dim only method. \n");
    f2 = GTEAndDerivs(f0, nmin, kmin, 0, 2);
    gsl_vector_set(point, 2, f2);
  }
  if (f2 > 1){
    //printf("f2 above 1, resetting using GTE. Second dim only method. \n");
    f2 = GTEAndDerivs(f0, nmax, kmax, 0, 2);
    gsl_vector_set(point, 2, f2);
  }
  
  
  if (nmin * pow(f1, 2) / f0 <= f2 && f2 <= nmax * pow(f1, 2) / f0){
    if (!(f2 > 0)){
    printf("f2 appropriate, but somehow, f2 !> 0. %E \n", f2);
  }
    return;
  }
  else if (f2 < nmin * pow(f1, 2) / f0){
    //printf("f2 too small, resetting. %f \n", f2);
    double newf2 = (1 + tol) * nmin * pow(f1, 2) / f0;
    if (!(f2 > 0)){
      printf("f2 was too small, but f2 !> 0. %E \n", f2);
    }
    gsl_vector_set(point, 2, newf2);
    return;
  }
  else if (f2 > nmax * pow(f1, 2) / f0){
    //printf("f2 too large, resetting. %f \n", f2);
    double newf2 = (1 - tol) * nmax * pow(f1, 2) / f0;
    if (!(f2 > 0)){
      printf("f2 was too big, but f2 !> 0. %E \n", f2);
    }
    gsl_vector_set(point, 2, newf2);
    return;
  }
  
  return; 
}

///
/// Due to the interdependence of the conditions on f1 and f2 on the first knot, if there values are chosen outside of the parameter space, resetting them can be challenging. Here I have written
/// a recursive numerical method which gives a way to reset such a parameter. First, if we find f1 to be in need of resetting, we use the resetfirstdim method defined above. Then, if the f2
/// parameter also needs resetting then we use the resetseconddim method above. If with both of these resets we find that the first knot is within the appropriate bounds, we terminate. Otherwise
/// we repeat this process for a maximum of 100 iterations. This method has so far worked successfully.
///
static void resetfirstandseconddim(
  gsl_vector* point, /// The vector of which we are resetting the variables f1 and f2 on
  double nmin,       /// The minimum allowed value of n for the first knot
  double nmax,       /// The maximum allowed value of n for the first knot
  double kmin,       /// The minimum allowed value of k for the first knot
  double kmax,       /// The maximum allowed value of k for the first knot
  double tol,        /// The same tol as in the above method, however now for both n and k criteria
  int itr
  )
{
  int maxitr = 100;
  if (itr >= maxitr){
    printf("Number of iterations exceeds %d, returning current values of point. (Perhaps I should throw an error here instead?)\n", maxitr);
    return;
  }
  
  double f0 = gsl_vector_get(point, 0);
  double f1 = gsl_vector_get(point, 1);
  double f2 = gsl_vector_get(point, 2);
  
  
  if (f2 <= 0){
    // We use nmin/kmin because if f2 is negative, it should be closest to its minimum values. Hence we use the extrema which result in the minimum values for f2.
    //printf("f2 below 0, resetting using GTE. Global only method. \n");
    f2 = GTEAndDerivs(f0, nmin, kmin, 0, 2);
    if (!(f2 > 0)){
      printf("lol, f2 reset is too small. nmin and kmin, f2 !> 0. %E \n", f2);
    }
    gsl_vector_set(point, 2, f2);
  }
  if (f2 > 1){
    //printf("f2 above 1, resetting using GTE. Global only method. \n");
    f2 = GTEAndDerivs(f0, nmax, kmax, 0, 2);
    if (!(f2 > 0)){
      printf("lol, f2 reset is too small. nmax and kmax, f2 !> 0. %E \n", f2);
    }
    gsl_vector_set(point, 2, f2);
  }
  
  if (!(f2 > 0)){
    printf("Global reset after initial checks, f2 !> 0. %E \n", f2);
  }
  
  
  double n = f2 * f0 / pow(f1, 2);
  double thisk = - f1 / pow(f0, n);
  
  if (kmin <= thisk && thisk <= kmax){
    if (nmin * pow(f1, 2) / f0 <= f2 && f2 <= nmax * pow(f1, 2) / f0){
      return;
    }
    resetseconddim(point, nmin, nmax, kmin, kmax, tol);
    itr++;
    resetfirstandseconddim(point, nmin, nmax, kmin, kmax, tol, itr);
    return;
  }
  else {
    resetfirstdim(point, kmin, kmax, tol);
    itr++;
    resetfirstandseconddim(point, nmin, nmax, kmin, kmax, tol, itr);
    return;
  }
  return;
  
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
  //printf("Variables for the resetting method \n");
  //printf("%f, %f, %f, %f, %f, %E, %E, %f, %f \n", fmin, fmax, nmin, nmax, ntol, kmin, kmax, ktol, segmentlength);
  if (dim == 0){
    double f0 = gsl_vector_get(point, dim);
    double val = resetinsidebounds(f0, fmin, fmax);
    gsl_vector_set(point, dim, val);
    return;
  }
  
  else if (dim == 1){
    
    resetfirstandseconddim(point, nmin, nmax, kmin, kmax, 0.1, 0);
    //double f1 = gsl_vector_get(point, 1);
    //XLAL_CHECK(!(f1 < 0), XLAL_EFAULT);
    return;
    
    /*
    
    double f0 = gsl_vector_get(point, dim - 1);
    double f1 = gsl_vector_get(point, dim);
    
    double lower = -kmax * pow(f0, nmax);
    double upper = -kmin * pow(f0, nmin);
    
    printf("For the first knot f1, %E, %E, %E \n", f1, lower, upper); 
    
    double val = resetinsidebounds(f1, lower, upper);
    gsl_vector_set(point, dim, val);
    return;
    */
    
    
  }
  else if (dim == 2){
    resetfirstandseconddim(point, nmin, nmax, kmin, kmax, 0.1, 0);
    return;
    
    /*
    double f0 = gsl_vector_get(point, dim - 2);
    double f1 = gsl_vector_get(point, dim - 1);
    double f2 = gsl_vector_get(point, dim);
    
    double lower = nmin * pow(f1, 2) / f0;
    double upper = nmax * pow(f1, 2) / f0;
    
    double val = resetinsidebounds(f2, lower, upper);
    gsl_vector_set(point, dim, val);
    return;
    */
  }
  else if (dim % 3 == 0){
  
    //printf("Second knot f0 \n");
    
    double f0n1 = gsl_vector_get(point, dim - 3);
    double f1n1 = gsl_vector_get(point, dim - 2);
    double f2n1 = gsl_vector_get(point, dim - 1);
    double f0 =   gsl_vector_get(point, dim);
    
    double nprev = f2n1 * f0n1 / pow(f1n1, 2);
    double kprev = - f1n1 / pow(f0n1, nprev);
    
    if (kprev > kmax || kprev < kmin){
      printf("Second knot f0. kprev outside global ranges, %E, %E, %E \n", kprev, kmin, kmax);
    }
    if (kprev <= kmax && kprev >= kmin){
      //printf("It made it inside!!! \n");
      //printf("\n");
    }
    
    double *nminmax = NMinMax(nprev, ntol, nmin, nmax, segmentlength);
    double *kminmax = KMinMax(kprev, ktol, kmin, kmax, segmentlength);
    
    ///printf("Variables used in F0BoundMinMax %f, %f, %f, %E, %E, %E \n", nprev, nminmax[0], nminmax[1], kprev, kminmax[0], kminmax[1]);
    
    double lower = F0BoundMinMax(f0n1, nminmax[0], nminmax[1], kminmax[0], kminmax[1], segmentlength, -1);
    double upper = F0BoundMinMax(f0n1, nminmax[1], nminmax[0], kminmax[1], kminmax[0], segmentlength,  1);
    
    double val = resetinsidebounds(f0, lower, upper);
    gsl_vector_set(point, dim, val);
    
  }
  else if (dim % 3 == 1){
    //printf("Second knot f1 \n");
    double f0n1 = gsl_vector_get(point, dim - 4);
    double f1n1 = gsl_vector_get(point, dim - 3);
    double f2n1 = gsl_vector_get(point, dim - 2);
    double f0 =   gsl_vector_get(point, dim - 1);
    double f1 =   gsl_vector_get(point, dim);
    
    double nprev = f2n1 * f0n1 / pow(f1n1, 2);
    double kprev = - f1n1 / pow(f0n1, nprev);
    
    if (kprev > kmax || kprev < kmin){
      printf("Second knot f1. kprev outside global ranges, %E, %E, %E \n", kprev, kmin, kmax);
    }
    
    double *nminmax = NMinMax(nprev, ntol, nmin, nmax, segmentlength);
    double *kminmax = KMinMax(kprev, ktol, kmin, kmax, segmentlength);
    
    double lower = F1BoundMinMax(f0, f0n1, nminmax[1], nminmax[0], kminmax[1], kminmax[0], segmentlength, -1);
    double upper = F1BoundMinMax(f0, f0n1, nminmax[0], nminmax[1], kminmax[0], kminmax[1], segmentlength,  1);
    
    
    ///printf("Lower upper %E, %E \n", lower, upper);
    
    
    double val = resetinsidebounds(f1, lower, upper);
    gsl_vector_set(point, dim, val);
  }
  else if (dim % 3 == 2){
  
    //printf("Second knot f2 \n");
    
    double f0n1 = gsl_vector_get(point, dim - 5);
    double f1n1 = gsl_vector_get(point, dim - 4);
    double f2n1 = gsl_vector_get(point, dim - 3);
    double f0 =   gsl_vector_get(point, dim - 2);
    double f1 =   gsl_vector_get(point, dim - 1);
    double f2 =   gsl_vector_get(point, dim);
    
    double nprev = f2n1 * f0n1 / pow(f1n1, 2);
    double kprev = - f1n1 / pow(f0n1, nprev);
    
    if (kprev > kmax || kprev < kmin){
      printf("kprev outside global ranges, %E, %E, %E \n", kprev, kmin, kmax);
    }
    
    double *nminmax = NMinMax(nprev, ntol, nmin, nmax, segmentlength);
    double *kminmax = KMinMax(kprev, ktol, kmin, kmax, segmentlength);
    
    double lower = F2BoundMinMax(f0, f0n1, f1, nminmax[0], kminmax[0], segmentlength, -1);
    double upper = F2BoundMinMax(f0, f0n1, f1, nminmax[1], kminmax[1], segmentlength,  1);
    
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
  
  //printf("Before and after F0 reset \n");
  //printvector(point);
  resetoutofboundspoint(point, fmin, fmax, nmin, nmax, ntol, kmin, kmax, ktol, segmentlength);
  //printvector(point);
  
  double f0n1 = gsl_vector_get(point, dim - 3);
  double f1n1 = gsl_vector_get(point, dim - 2);
  double f2n1 = gsl_vector_get(point, dim - 1);
  
  double nprev = f2n1 * f0n1 / pow(f1n1, 2);
  double kprev = - f1n1 / pow(f0n1, nprev);
  
  double *nminmax = NMinMax(nprev, ntol, nmin, nmax, segmentlength);
  double *kminmax = KMinMax(kprev, ktol, kmin, kmax, segmentlength);
  
  if (upperlower == 1){
    double upperbound = F0BoundMinMax(f0n1, nminmax[1], nminmax[0], kminmax[1], kminmax[0], segmentlength, 1);
    //printf("F0 upper bound %f, %f \n", f0n1, upperbound);
    return upperbound;
  }
  else if (upperlower == -1){
    double lowerbound = F0BoundMinMax(f0n1, nminmax[0], nminmax[1], kminmax[0], kminmax[1], segmentlength, -1);
    //printf("F0 lower bound %f, %f \n", f0n1, lowerbound);
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
  
  //printf("Before and after F1 reset \n");
  //printvector(point);
  resetoutofboundspoint(point, fmin, fmax, nmin, nmax, ntol, kmin, kmax, ktol, segmentlength);
  //printvector(point);
  
  double f0n1 = gsl_vector_get(point, dim - 4);
  double f1n1 = gsl_vector_get(point, dim - 3);
  double f2n1 = gsl_vector_get(point, dim - 2);
  
  double f0 = gsl_vector_get(point, dim - 1);
  
  ///printf("%f, %f, %f, %f \n", f0n1, f1n1, f2n1, f0);
  
  double nprev = f2n1 * f0n1 / pow(f1n1, 2);
  double kprev = - f1n1 / pow(f0n1, nprev);
  
  double *nminmax = NMinMax(nprev, ntol, nmin, nmax, segmentlength);
  double *kminmax = KMinMax(kprev, ktol, kmin, kmax, segmentlength);
  
  //printf("nminmax %f, %f \n", nminmax[0], nminmax[1]);
  //printf("kminmax %f, %f \n", kminmax[0], kminmax[1]);
  
  if (upperlower == 1){
    double upperbound = F1BoundMinMax(f0, f0n1, nminmax[0], nminmax[1], kminmax[0], kminmax[1], segmentlength, 1);
    //printf("F1 upper bound %f, %f \n", f1n1, upperbound);
    return upperbound;
  }
  else if (upperlower == -1){
    double lowerbound = F1BoundMinMax(f0, f0n1, nminmax[1], nminmax[0], kminmax[1], kminmax[0], segmentlength, -1);
    //printf("F1 lower bound %f, %f \n", f1n1, lowerbound);
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
  
  //printf("Before and after F2 reset \n");
  //printvector(point);
  resetoutofboundspoint(point, fmin, fmax, nmin, nmax, ntol, kmin, kmax, ktol, segmentlength);
  //printvector(point);
  
  double f0n1 = gsl_vector_get(point, dim - 5);
  double f1n1 = gsl_vector_get(point, dim - 4);
  double f2n1 = gsl_vector_get(point, dim - 3);
  
  double f0 = gsl_vector_get(point, dim - 2);
  double f1 = gsl_vector_get(point, dim - 1);
  
  double nprev = f2n1 * f0n1 / pow(f1n1, 2);
  double kprev = - f1n1 / pow(f0n1, nprev);
  
  double *nminmax = NMinMax(nprev, ntol, nmin, nmax, segmentlength);
  double *kminmax = KMinMax(kprev, ktol, kmin, kmax, segmentlength);
  
  if (upperlower == 1){
    double upperbound = F2BoundMinMax(f0, f0n1, f1, nminmax[1], kminmax[1], segmentlength, 1);
    printf("F2 upper bound %E, %E \n", f2n1, upperbound);
    return upperbound;
  }
  else if (upperlower == -1){
    double lowerbound = F2BoundMinMax(f0, f0n1, f1, nminmax[0], kminmax[0], segmentlength, -1);
    printf("F2 lower bound %E, %E \n", f2n1, lowerbound);
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
  XLALSetLatticeTilingPaddingFlags(tiling, 0, LATTICE_TILING_PAD_NONE);
  
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
  
  XLALSetLatticeTilingPaddingFlags(tiling, 1, LATTICE_TILING_PAD_NONE);
  XLALSetLatticeTilingPaddingFlags(tiling, 2, LATTICE_TILING_PAD_NONE);
  
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
    
    XLALSetLatticeTilingPaddingFlags(tiling, dimindex,     LATTICE_TILING_PAD_NONE);
    XLALSetLatticeTilingPaddingFlags(tiling, dimindex + 1, LATTICE_TILING_PAD_NONE);
    XLALSetLatticeTilingPaddingFlags(tiling, dimindex + 2, LATTICE_TILING_PAD_NONE);
    
  }
  
  printf("Bounds set \n");
  printf("\n");
  return XLAL_SUCCESS;
}


