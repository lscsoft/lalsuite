//
// Copyright (C) 2021, Ben Grace
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307 USA
//

#include <stdio.h>
#include <math.h>
#include <lal/LatticeTiling.h>
#include <lal/LogPrintf.h>
#include <lal/PiecewiseModel.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif


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


///
/// Information required to determine the bounds on the parameters on the first knot
///
typedef struct tagFirstKnotBoundInfo{
  double fmin;             ///< Minimum starting frequency
  double fmax;             ///< Maximum starting frequency
  double nmin;             ///< Minimum braking index
  double nmax;             ///< Maximum braking index
  double kmin;             ///< Minimum k value
  double kmax;             ///< Maximum k value
  int minmax;              ///< +1 for upper bound, -1 for lower bound
  int reset;               ///< +1 for resetting point methods to be used, -1 for not
} FirstKnotBoundInfo;

///
/// A struct containing the relevant information for calculating upper and lower bounds on a specific piecewise segment
///
typedef struct tagPiecewiseBoundInfo{
  double fmin;             ///< Minimum starting frequency
  double fmax;             ///< Maximum starting frequency
  double nmin;             ///< Minimum braking index
  double nmax;             ///< Maximum braking index
  double ntol;             ///< Braking index tolerance (percentage per second) between adjacent knots
  double kmin;             ///< Minimum k value
  double kmax;             ///< Maximum k value
  double ktol;             ///< k value tolerance (percentage per second) between adjacent knots
  double segmentlength;    ///< The duration of the current piecewise segment
  int upperlower;          ///< +1 for calculating upper bound, -1 for calculating lower bound
  int reset;               ///< +1 for resetting point methods to be used, -1 for not
} PiecewiseBoundInfo;

///
/// The general torque equation and its first two derivatives 
///
static double GTEAndDerivs(
  double f0,               ///< Initial frequency
  double n,                ///< Braking index
  double k,                ///< k value
  double t,                ///< Time at which to evaluate the GTE
  int d                    ///< Derivative order (d <= 2)
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
/// Alters the nminmax array to give the allowed range of the braking index for a knot given the braking index value on the previous knot
///
static int NMinMax(
  double * nminmax,        ///< Nminmax array
  double n,                ///< Previous braking index
  double ntol,             ///< Braking index tolerance (percentage per second)
  double nmin,             ///< Minimum allowed braking index
  double nmax,             ///< Maximum allowed braking index
  double segmentlength     ///< Time difference between current knot and previous knot
  )
{

  XLAL_CHECK(nmin <= nmax, XLAL_EINVAL, "nmin greater than nmax, [nmin, nmax] = [%f, %f]", nmin, nmax);
  //printf("Nminmax values are, %f, %f, %f \n", n, nmin, nmax);
  double nextnmin = n * (1 - ntol * segmentlength);
  double nextnmax = n;
  //printf("Next n values are, %f, %f \n", nextnmin, nextnmax);
  ///< If the code is running properly (at least with padding flags turned off) this shouldn't ever really come up, but if it does, it doesn't necessarily 
  ///< mean things will go wrong. It usually just means that nextnmax/nextnmin are very close to either nmin or nmax
  if (nextnmax < nmin || nextnmin > nmax){
    XLAL_PRINT_WARNING("Calculated n ranges outside of global range, %E, %E, %E, %E, %E, %E, %E \n", n, nmin, nmax, nextnmin, nextnmax, n - nmin, n - nmax);
  }
  
  if (nextnmin < nmin){
    nextnmin = nmin;
  }
  if (nextnmax > nmax){
    nextnmax = nmax;
  }
  
  nminmax[0] = nextnmin;
  nminmax[1] = nextnmax;
  
  XLAL_CHECK(nextnmin <= nextnmax, XLAL_EINVAL, "Calculated braking index ranges incorrect, [nextnmin, nextnmax] = [%f, %f]", nextnmin, nextnmax);
  
  return 1;
}

///
/// Alters the kminmax array to give the allowed range of the k value for a knot given the k value on the previous knot
///
static int KMinMax(
  double * kminmax,        ///< kminmax array
  double k,                ///< Previous k value index
  double ktol,             ///< k value tolerance (percentage per second)
  double kmin,             ///< Minimum allowed k value
  double kmax,             ///< Maximum allowed k value
  double segmentlength     ///< Time difference between current knot and previous knot
  )
{
  
  XLAL_CHECK(kmin <= kmax, XLAL_EINVAL, "kmin greater than kmax, [kmin, kmax] = [%E, %E]", kmin, kmax);
  
  //printf("Kminmax values are, %E, %E, %E \n", k, kmin, kmax);
  double nextkmin = k * (1 - ktol * segmentlength);
  double nextkmax = k;
  //printf("Next k values are, %E, %E \n", nextkmin, nextkmax);
  
  ///< If the code is running properly (at least with padding flags turned off) this shouldn't ever really come up, but if it does, it doesn't necessarily 
  ///< mean things will go wrong. It usually just means that nextkmax/nextkmin are very close to either kmin or kmax
  if (nextkmax < kmin || nextkmin > kmax){
    XLAL_PRINT_WARNING("Calculated k ranges outside of global range, %E, %E, %E, %E, %E, %E, %E \n", k, kmin, kmax, nextkmin, nextkmax, k - kmin, k - kmax);
  }
  
  if (nextkmin < kmin){
    nextkmin = kmin;
  }
  if (nextkmax > kmax){
    nextkmax = kmax;
  }
  
  kminmax[0] = nextkmin;
  kminmax[1] = nextkmax;
  
  XLAL_CHECK(nextkmin <= nextkmax, XLAL_EINVAL, "Calculated k value ranges incorrect, [nextkmin, nextkmax] = [%E, %E]", nextkmin, nextkmax);

  return 1;
}

///
/// Calculates the minimum and maximum bounds for a frequency parameter
///
static double F0BoundMinMax(
  double f0,               ///< Frequency value of the previous knot
  double UNUSED na,        ///< A braking index value. For calculating upper bound, na > nb. For calculating lower bound na < nb.
  double nb,               ///< A braking index value. For calculating upper bound, na > nb. For calculating lower bound na < nb. Revision: For upper bound, nb should be minimised
  double UNUSED ka,        ///< A k value. For calculating upper bound, ka > kb. For calculating lower bound ka < kb
  double kb,               ///< A k value. For calculating upper bound, ka > kb. For calculating lower bound ka < kb. Revision: For upper bound, kb should be minimised
  double seglength,        ///< Time difference between this knot and the previous knot
  int UNUSED minmax        ///< +1 for calculating upper bound, -1 for calculating lower bound
)
{
  
  ///< Parameter range optimisation for frequency parameter
  double paramrange = GTEAndDerivs(f0, nb, kb, seglength, 0);
  return paramrange;
}

///
/// Function which determines the minimum/maximum value of the first derivative parameters
///
static double F1BoundMinMax(
  double UNUSED f0,        ///< Frequency parameter at the corresponding knot
  double fprev,            ///< Frequency parameter at the previous knot
  double na,               ///< A braking index value. For calculating upper bound, na < nb. For calculating lower bound na > nb. Revision: For upper bound, na should be minimised
  double UNUSED nb,        ///< A braking index value. For calculating upper bound, na < nb. For calculating lower bound na > nb
  double ka,               ///< A k value. For calculating upper bound, ka < kb. For calculating lower bound ka > kb. Revision: For upper bound, ka should be minimised
  double UNUSED kb,        ///< A k value. For calculating upper bound, ka < kb. For calculating lower bound ka > kb
  double seglength,        ///< Time difference between this knot and the previous
  int minmax               ///< +1 for upper bound, -1 for lower bound
  )
{
  double prangecondition1 = -ka * pow(f0, na);
  double prangecondition2 = GTEAndDerivs(fprev, na, ka, seglength, 1);
  
  gsl_vector* vals = gsl_vector_alloc(2);
  gsl_vector_set(vals, 0, prangecondition1);
  gsl_vector_set(vals, 1, prangecondition2);
  
  double rtn = NAN;
  
  if(minmax == 1){
    double max = gsl_vector_min(vals);
    rtn = max;
  }
  else if(minmax == -1){
    double min = gsl_vector_max(vals);
    rtn = min;
  }
  
  gsl_vector_free(vals);
  return rtn; 
}

///
/// Calculates the bound for the second derivative frequency parameter 
///
static double F2BoundMinMax(
  double f0,               ///< Frequency parameter on this knot
  double fprev,            ///< Frequency parameter on the previous knot
  double f1,               ///< First derivative frequency parameter on this knot
  double n,                ///< Optimal braking index value. For upper bound n should be maximised, for lower bound n should be minimised
  double ka,               ///< A k value. For calculating upper bound, ka > kb. For calculating lower bound ka < kb
  double kb,               ///< A k value. For calculating upper bound, ka > kb. For calculating lower bound ka < kb 
  double seglength,        ///< Time difference between this knot and the previous
  int minmax               ///< +1 for upper bound, -1 for lower bound
)
{

  XLAL_CHECK((ka <= kb && minmax == -1) || (ka >= kb && minmax == 1), XLAL_EINVAL, "F2BoundMinMax being called incorrectly, with [ka, kb, minmax] = [%E, %E, %d]", ka, kb, minmax);
  
  ///< Braking index and parameter range conditions.
  double bindexcriteria = n * pow(f1, 2) / f0;
  double prangecriteria = GTEAndDerivs(fprev, n, ka, seglength, 2);
  double kcriteria = - pow(f1, 2) * log(- kb / f1) / (f0 * log(f0));
  
  gsl_vector* vals = gsl_vector_alloc(3);
  gsl_vector_set(vals, 0, bindexcriteria);
  gsl_vector_set(vals, 1, prangecriteria);
  gsl_vector_set(vals, 2, kcriteria);
  
  double rtn = NAN;
  
  if(minmax == 1){
    double max = gsl_vector_min(vals);
    rtn =  max;
  }
  else if(minmax == -1){
    double min = gsl_vector_max(vals);
    rtn = min;
  }
  
  gsl_vector_free(vals);
  return rtn;
}

///
/// Takes an input val, and returns val to its closest value within the range [valmin, valmax].
///
static double resetinsidebounds(
  double val,              ///< Value
  double valmin,           ///< Lower bound of range
  double valmax            ///< Upper bound of range
  )
{
  
  if (valmin <= val && val <= valmax){
    return val;
  }
  else if (valmin > valmax){
    
    XLAL_CHECK(((valmin - valmax) / valmax < pow(10, -4)), XLAL_EINVAL, "Valmin bigger than valmax by more than accepted tolerance: [val, valmin, valmax] = [%E, %E, %E]", val, valmin, valmax);
    return valmax;
    
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
/// Sets a value val to one of its extreme values if it lies just outside its acceptable range. 'Just outside' defined by valtol
///
static double resetvalwithintol(
  double val,              ///< Value
  double valmin,           ///< Lower bound of range
  double valmax,           ///< Upper bound of range
  double UNUSED valtol     ///< Percentage tolerance val may lie outside of range
  )
{
  if (val > valmax && val * (1 - 0.001) <= valmax){
    val = valmax;
  }
  else if (val < valmin && val * (1 + 0.001) >= valmin){
    val = valmin;
  }
  return val;
}

///
/// A method which will reset a given point to be within our parameter space. Method not used when LatticeTilingPaddingFlags = LATTICE_TILING_PADDING_FLAGS_NONE
///
static void resetdimonpoint(
  gsl_vector* point,       ///< The point which we are resetting to be within the parameter space
  int dim,                 ///< The dimension which we are resetting
  double fmin,             ///< Minimum starting frequency
  double fmax,             ///< Maximum starting frequency
  double nmin,             ///< Minimum allowed braking index at current knot
  double nmax,             ///< Maximum allowed braking index at current knot
  double ntol,             ///< Braking index tolerance (percentage per second)
  double kmin,             ///< Minimum allowed k value at current knot
  double kmax,             ///< Maximum allowed k value at current knot
  double ktol,             ///< k value tolerance (percentage per second)
  double segmentlength     ///< The length of the segment we are working on
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
    
    ///< We need a special case when f1 is equal to one of its bounds. The reason being a floating point error arises which results in calculated k values being
    ///< off by only a small percentage. This is explained in further comments below.
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
    
    ///< The parameter range criteria involving fprev is not valid for at the first knot (as there are no preceding knots). To then discount this criteria in the 
    ///< F2BoundMinMax function, we set fprev to extreme values to ensure this criteria is never selected without having to write a new method
    double lower = F2BoundMinMax(f0, 0.0000000001, f1, nmin, kmin, kmax, 0, -1);
    double upper = F2BoundMinMax(f0, 100000000000, f1, nmax, kmax, kmin, 0,  1);

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
    
    ///< In the instant that f1n1 has a value at one of its extrema, the upper and lower bound on f2 are the same point (i.e. there is only one allowed value for f1)
    ///< However, given the equations for f2 (either starting from braking index requirements or k value requirements), small floating point error arises between the 
    ///< two end results. For this reason, we must choose one equation to use to calculate f2 (in this case we have chosen the braking index equation) and then manually
    ///< edit the value of k and n we decide to use for our checks on subsequent knots from that point. This is the reason for the below conditions. The typically error
    ///< between kprev and either kmax and kmin is on the order of ~10^-23, likewise for nprev and nmin/nmax the error is of the order ~10^-16.
    
    nprev = resetvalwithintol(nprev, nmin, nmax, ntol);
    kprev = resetvalwithintol(kprev, kmin, kmax, ktol);
    
    double nminmax[2];
    double kminmax[2];
    
    NMinMax(nminmax, nprev, 0.001, nmin, nmax, segmentlength);
    KMinMax(kminmax, kprev, 0.001, kmin, kmax, segmentlength);
    
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
    
    ///< In the instant that f1n1 has a value at one of its extrema, the upper and lower bound on f2 are the same point (i.e. there is only one allowed value for f1)
    ///< However, given the equations for f2 (either starting from braking index requirements or k value requirements), small floating point error arises between the 
    ///< two end results. For this reason, we must choose one equation to use to calculate f2 (in this case we have chosen the braking index equation) and then manually
    ///< edit the value of k and n we decide to use for our checks on subsequent knots from that point. This is the reason for the below conditions. The typically error
    ///< between kprev and either kmax and kmin is on the order of ~10^-23, likewise for nprev and nmin/nmax the error is of the order ~10^-16.
    
    nprev = resetvalwithintol(nprev, nmin, nmax, ntol);
    kprev = resetvalwithintol(kprev, kmin, kmax, ktol);
    
    double nminmax[2];
    double kminmax[2];
    
    NMinMax(nminmax, nprev, 0.001, nmin, nmax, segmentlength);
    KMinMax(kminmax, kprev, 0.001, kmin, kmax, segmentlength);
    
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
    
    ///< In the instant that f1n1 has a value at one of its extrema, the upper and lower bound on f2 are the same point (i.e. there is only one allowed value for f1)
    ///< However, given the equations for f2 (either starting from braking index requirements or k value requirements), small floating point error arises between the 
    ///< two end results. For this reason, we must choose one equation to use to calculate f2 (in this case we have chosen the braking index equation) and then manually
    ///< edit the value of k and n we decide to use for our checks on subsequent knots from that point. This is the reason for the below conditions. The typically error
    ///< between kprev and either kmax and kmin is on the order of ~10^-23, likewise for nprev and nmin/nmax the error is of the order ~10^-16.
    
    nprev = resetvalwithintol(nprev, nmin, nmax, ntol);
    kprev = resetvalwithintol(kprev, kmin, kmax, ktol);
    
    double nminmax[2];
    double kminmax[2];
    
    NMinMax(nminmax, nprev, 0.001, nmin, nmax, segmentlength);
    KMinMax(kminmax, kprev, 0.001, kmin, kmax, segmentlength);
    
    double lower = F2BoundMinMax(f0, f0n1, f1, nminmax[0], kminmax[0], kminmax[1], segmentlength, -1);
    double upper = F2BoundMinMax(f0, f0n1, f1, nminmax[1], kminmax[1], kminmax[0], segmentlength,  1);

    double val = resetinsidebounds(f2, lower, upper);
    gsl_vector_set(point, dim, val);
  }
}

///
/// Resets a point such that it is within the bounds of our parameter space. Method not used when LatticeTilingPaddingFlags = LATTICE_TILING_PADDING_FLAGS_NONE
///
static void resetoutofboundspoint(
  gsl_vector* point,       ///< The point which we are resetting to be within the parameter space
  double fmin,             ///< Minimum starting frequency
  double fmax,             ///< Maximum starting frequency
  double nmin,             ///< Minimum allowed braking index at current knot
  double nmax,             ///< Maximum allowed braking index at current knot
  double ntol,             ///< Braking index tolerance (percentage per second)
  double kmin,             ///< Minimum allowed k value at current knot
  double kmax,             ///< Maximum allowed k value at current knot
  double ktol,             ///< k value tolerance (percentage per second)
  double segmentlength     ///< The length of the segment we are working on
  )
{
  int dim = point->size;
  printf("Performing reset \n");
  ///printf("Before: ");
  ///printvector(point);
  for (int i = 0; i < dim; ++i){
    resetdimonpoint(point, i, fmin, fmax, nmin, nmax, ntol, kmin, kmax, ktol, segmentlength);
  }
  ///printf("After:   ");
  ///printvector(point);
  ///printf("\n");
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
  gsl_vector_memcpy(point, pointorig);
  
  const FirstKnotBoundInfo* info = (const FirstKnotBoundInfo*)data;
  
  double fmin = info->fmin;
  double fmax = info->fmax;
  double nmin = info->nmin;
  double nmax = info->nmax;
  double kmin = info->kmin;
  double kmax = info->kmax;
  int minmax = info->minmax;
  int reset = info->reset;
  
  double f0 = gsl_vector_get(point, 0);
  
  if (reset == 1){
    f0 = resetinsidebounds(gsl_vector_get(point, 0), fmin, fmax);
    gsl_vector_set(point, 0, f0);
  }
  
  double rtn = NAN;
  
  if (dim == 1){
    
    if (minmax == 1){
      double f1 = -kmin * pow(f0, nmin);
      rtn = f1;
    }
    else if (minmax == -1){
      double f1 = -kmax * pow(f0, nmax);
      rtn = f1;
    }
  }
  else if (dim == 2){
    
    double f1 = gsl_vector_get(point, 1);
    
    if (reset == 1){
    
      double lower = -kmax * pow(f0, nmax);
      double upper = -kmin * pow(f0, nmin);
      
      f1 = resetinsidebounds(gsl_vector_get(point, 1), lower, upper);
      gsl_vector_set(point, 1, f1);
    }
    
    if (minmax == 1){
      double f2 = F2BoundMinMax(f0,  1E10, f1, nmax, kmax, kmin, 0, 1);
      rtn = f2;
    }
    else if (minmax == -1){
      double f2 = F2BoundMinMax(f0, 1E-10, f1, nmin, kmin, kmax, 0, -1);
      rtn = f2;
    }
  }
  
  gsl_vector_free(point);
  return rtn;
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
  gsl_vector_memcpy(point, pointorig);
  
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
  int reset = info->reset;
  
  if (reset == 1){
    resetoutofboundspoint(point, fmin, fmax, nmin, nmax, ntol, kmin, kmax, ktol, segmentlength);
  }
  
  double f0n1 = gsl_vector_get(point, dim - 3);
  double f1n1 = gsl_vector_get(point, dim - 2);
  double f2n1 = gsl_vector_get(point, dim - 1);
  
  double nprev = f2n1 * f0n1 / pow(f1n1, 2);
  double kprev = - f1n1 / pow(f0n1, nprev);
  
  nprev = resetvalwithintol(nprev, nmin, nmax, ntol);
  kprev = resetvalwithintol(kprev, kmin, kmax, ktol);
  
  double nminmax[2];
  double kminmax[2];
  
  NMinMax(nminmax, nprev, 0.001, nmin, nmax, segmentlength);
  KMinMax(kminmax, kprev, 0.001, kmin, kmax, segmentlength);
  
  double rtn = NAN;
  
  if (upperlower == 1){
    double upperbound = F0BoundMinMax(f0n1, nminmax[1], nminmax[0], kminmax[1], kminmax[0], segmentlength, 1);
    rtn = upperbound;
  }
  else if (upperlower == -1){
    double lowerbound = F0BoundMinMax(f0n1, nminmax[0], nminmax[1], kminmax[0], kminmax[1], segmentlength, -1);
    rtn = lowerbound;
  }
  
  gsl_vector_free(point);
  return rtn;
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
  gsl_vector_memcpy(point, pointorig);
  
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
  int reset = info->reset;
  
  if (reset == 1){
    resetoutofboundspoint(point, fmin, fmax, nmin, nmax, ntol, kmin, kmax, ktol, segmentlength);
  }
  
  double f0n1 = gsl_vector_get(point, dim - 4);
  double f1n1 = gsl_vector_get(point, dim - 3);
  double f2n1 = gsl_vector_get(point, dim - 2);
  
  double f0 = gsl_vector_get(point, dim - 1);
  
  double nprev = f2n1 * f0n1 / pow(f1n1, 2);
  double kprev = - f1n1 / pow(f0n1, nprev);
  
  nprev = resetvalwithintol(nprev, nmin, nmax, ntol);
  kprev = resetvalwithintol(kprev, kmin, kmax, ktol);
  
  double nminmax[2];
  double kminmax[2];
  
  NMinMax(nminmax, nprev, 0.001, nmin, nmax, segmentlength);
  KMinMax(kminmax, kprev, 0.001, kmin, kmax, segmentlength);
  
  double rtn = NAN;
  
  if (upperlower == 1){
    double upperbound = F1BoundMinMax(f0, f0n1, nminmax[0], nminmax[1], kminmax[0], kminmax[1], segmentlength, 1);
    rtn = upperbound;
  }
  else if (upperlower == -1){
    double lowerbound = F1BoundMinMax(f0, f0n1, nminmax[1], nminmax[0], kminmax[1], kminmax[0], segmentlength, -1);
    rtn = lowerbound;
  }
  
  gsl_vector_free(point);
  return rtn;
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
  gsl_vector_memcpy(point, pointorig);
  
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
  int reset = info->reset;
  
  if (reset == 1){
    resetoutofboundspoint(point, fmin, fmax, nmin, nmax, ntol, kmin, kmax, ktol, segmentlength);
  }
  
  double f0n1 = gsl_vector_get(point, dim - 5);
  double f1n1 = gsl_vector_get(point, dim - 4);
  double f2n1 = gsl_vector_get(point, dim - 3);
  
  double f0 = gsl_vector_get(point, dim - 2);
  double f1 = gsl_vector_get(point, dim - 1);
  
  double nprev = f2n1 * f0n1 / pow(f1n1, 2);
  double kprev = - f1n1 / pow(f0n1, nprev);
  
  nprev = resetvalwithintol(nprev, nmin, nmax, ntol);
  kprev = resetvalwithintol(kprev, kmin, kmax, ktol);
  
  double nminmax[2];
  double kminmax[2];
  
  NMinMax(nminmax, nprev, 0.001, nmin, nmax, segmentlength);
  KMinMax(kminmax, kprev, 0.001, kmin, kmax, segmentlength);
  
  double rtn = NAN;
  
  if (upperlower == 1){
    double upperbound = F2BoundMinMax(f0, f0n1, f1, nminmax[1], kminmax[1], kminmax[0], segmentlength, 1);
    rtn = upperbound;
  }
  else if (upperlower == -1){
    double lowerbound = F2BoundMinMax(f0, f0n1, f1, nminmax[0], kminmax[0], kminmax[1], segmentlength, -1);
    rtn = lowerbound;
  }
  
  gsl_vector_free(point);
  return rtn;
}

///
/// Converts a given tau value to a k value (or a given k value to a tau value)
///
static double ktauconversion(
  double f0,               ///< Frequency at time t = 0. For calculating kmax/kmin from tau values, f0 should be maximised. Likewise for calculating taumin/taumax from k values 
  double n,                ///< A braking index. For calculating kmax/kmin from tau values, n should be maximised. Likewise for calculating taumin/taumax from k values 
  double kortau            ///< The k or tau values we wish to convert. For calculating kmax, taumin should be used. For calculating kmin, taumax should be used.
  )
{
  double numerator = pow(2, n - 1) - 1;
  double denominator = (n - 1) * pow(f0, n - 1) * kortau;
  
  double ktau = numerator / denominator;
  
  return ktau;
}

///
/// Returns the upper or lower bound for the dimension 'dim' given values for all previous dimensions, point_up_to_dim, as well as the relevant parameter space information
///
double XLALPiecewiseParameterBounds(
  const size_t      dim,              /// The dimension of the parameter we wish to calculate the bounds for
  const gsl_vector* point_up_to_dim,  /// The point up the given dimension
  const int         upperlower,       /// +1 to return upper bound, -1 to return lower bound
  const double      fmin,             /// Global maximum frequency
  const double      fmax,             /// Global minimum frequency
  const double      nmin,             /// Minimum braking index
  const double      nmax,             /// Maximum braking index
  const double      ntol,             /// Tolerance (percentage per second) between braking indices on adjacent knots
  const double      kmin,             /// Minimum k value
  const double      kmax,             /// Maximum k value
  const double      ktol,             /// Tolerance (percentage per second) between k values on adjacent knots
  const double      seglength         /// The length of the segment. The time between the previous knot and the knot that the parameter 'dim' resides on
  )
{ 
  
  double bound;
  size_t zero = 0;
  gsl_matrix *cache = gsl_matrix_alloc(zero, zero);
  
  
  if (dim == 0){
    if (upperlower == 1){
      return fmax;
    }
    else {
      return fmin;
    }
  }
  
  if (dim < 3){
    FirstKnotBoundInfo XLAL_INIT_DECL( info_first_knot );
    
    info_first_knot.fmin = fmin;
    info_first_knot.fmax = fmax;
    info_first_knot.nmin = nmin;
    info_first_knot.nmax = nmax;
    info_first_knot.kmin = kmin;
    info_first_knot.kmax = kmax;
    
    info_first_knot.minmax = upperlower;
    
    info_first_knot.reset = 1;
    
    
    bound = FirstKnotDerivBound(&info_first_knot, dim, cache, point_up_to_dim);
  }
  else {
    PiecewiseBoundInfo XLAL_INIT_DECL( info_knot );
    
    info_knot.fmin = fmin;
    info_knot.fmax = fmax;
    info_knot.nmin = nmin;
    info_knot.nmax = nmax;
    info_knot.ntol = ntol;
    info_knot.kmin = kmin;
    info_knot.kmax = kmax;
    info_knot.ktol = ktol;
    info_knot.segmentlength = seglength;
    
    info_knot.upperlower = upperlower;
    
    info_knot.reset = 1;
    
    if (dim % 3 == 0){
      bound = F0Bound(&info_knot, dim, cache, point_up_to_dim);
    } else if (dim % 3 == 1){
      bound = F1Bound(&info_knot, dim, cache, point_up_to_dim);
    } else if (dim % 3 == 2){
      bound = F2Bound(&info_knot, dim, cache, point_up_to_dim);
    }
  }
  
  gsl_matrix_free(cache);
  return bound;
}



///
/// Sets the bounds for the piecewise model
///
int XLALSetLatticeTilingPiecewiseBounds(
  LatticeTiling* tiling,        ///< Lattice tiling
  const double fmin,            ///< Minimum initial frequency
  const double fmax,            ///< Maximum initial frequency
  const double fmaxtrue,        ///< Maximum spin frequency with which to calculate k value ranges with (useful when computing tiles in parrallel and fmax != fmaxtrue)
  const double nmin,            ///< Minimum braking index
  const double nmax,            ///< Maximum braking index
  const double nmin0,           ///< Minimum braking index for the first knot. Useful if you want to brake up a search into partitions separated by templates with braking indices within a certain range
  const double nmax0,           ///< Maximum braking index for the first knot. Useful if you want to brake up a search into partitions separated by templates with braking indices within a certain range
  const double ntol,            ///< Tolerance (percentage per second) between braking indices on adjacent knots
  const double taumin,          ///< Minimum tau value
  const double taumax,          ///< Maximum tau value
  const double ktol,            ///< Tolerance (percentage per second) between k values on adjacent knots
  const gsl_vector* knots,      ///< List of knots
  const int finalknot,          ///< The number of the final knot
  const gsl_vector* bboxpad,    ///< Vector containing fractional bounding box padding
  const gsl_vector_int* intpad ///< Vector containing number of integer points to use for padding
  )
{
  // Converting tau values to k values
  double kmin = ktauconversion(fmaxtrue, nmax, taumax);
  double kmax = ktauconversion(fmaxtrue, nmax, taumin);
  
  LogPrintf(LOG_DEBUG, "kmin and kmax %E, %E \n", kmin, kmax);
  
  XLAL_CHECK(tiling != NULL,     XLAL_EINVAL);
  XLAL_CHECK(fmin   < fmax,      XLAL_EINVAL, "Bad frequency range: [%f, %f]", fmin, fmax);
  XLAL_CHECK(nmin   < nmax,      XLAL_EINVAL, "Bad braking index range: [%f, %f]", nmin, nmax);
  XLAL_CHECK(fmax   <= fmaxtrue, XLAL_EINVAL, "fmax larger than fmaxtrue: [%f, %f]", fmax, fmaxtrue);
  XLAL_CHECK(nmin0  < nmax0,     XLAL_EINVAL, "Bad braking index knot0 range: [%f, %f]", nmin0, nmax0);
  XLAL_CHECK(nmin0  >= nmin,     XLAL_EINVAL, "nmin0 smaller than nmin: [%f, %f]", nmin0, nmin);
  XLAL_CHECK(nmax0  <= nmax,     XLAL_EINVAL, "nmax0 greater than nmax: [%f, %f]", nmax0, nmax);
  XLAL_CHECK(taumin < taumax,    XLAL_EINVAL, "Bad tau range: [%f, %f]", taumin, taumax);
  XLAL_CHECK(kmin   < kmax,      XLAL_EINVAL, "Bad k range: [%f, %f]", kmin, kmax);
  
  FirstKnotBoundInfo XLAL_INIT_DECL( info_first_knot_lower );
  FirstKnotBoundInfo XLAL_INIT_DECL( info_first_knot_upper );
  
  info_first_knot_lower.fmin = info_first_knot_upper.fmin = fmin;
  info_first_knot_lower.fmax = info_first_knot_upper.fmax = fmax;
  info_first_knot_lower.nmin = info_first_knot_upper.nmin = nmin0;
  info_first_knot_lower.nmax = info_first_knot_upper.nmax = nmax0;
  info_first_knot_lower.kmin = info_first_knot_upper.kmin = kmin;
  info_first_knot_lower.kmax = info_first_knot_upper.kmax = kmax;
  
  info_first_knot_lower.minmax = -1;
  info_first_knot_upper.minmax = 1;
  
  /// Whether we reset points to be within the parameter space or not. 1 is reset, -1 is not reset
  int reset_value = -1;
  
  printf("Padding flag vectors are: \n");
  
  size_t vectorlength = bboxpad->size;
  gsl_vector* bboxpad_copy = gsl_vector_alloc(vectorlength);
  gsl_vector_memcpy(bboxpad_copy, bboxpad);
  
  gsl_vector* intpad_copy = gsl_vector_alloc(vectorlength);
  gsl_vector_memcpy(intpad_copy, bboxpad);
  
  printvector(bboxpad_copy);
  printvector(intpad_copy);
  
  
  info_first_knot_lower.reset = info_first_knot_upper.reset = reset_value;
  
  ///< Setting the first knot bounds
  XLALSetLatticeTilingConstantBound(tiling, 0, fmin, fmax);
  XLALSetLatticeTilingPadding(tiling, 0, gsl_vector_get(bboxpad, 0), gsl_vector_get(bboxpad, 0), gsl_vector_int_get(intpad, 0), gsl_vector_int_get(intpad, 0));
  
  XLAL_CHECK(XLALSetLatticeTilingBound(tiling, 1, FirstKnotDerivBound, sizeof( info_first_knot_lower ), &info_first_knot_lower, &info_first_knot_upper) == XLAL_SUCCESS, XLAL_EFAILED);
  XLALSetLatticeTilingPadding(tiling, 1, gsl_vector_get(bboxpad, 1), gsl_vector_get(bboxpad, 1), gsl_vector_int_get(intpad, 1), gsl_vector_int_get(intpad, 1));
  
  XLAL_CHECK(XLALSetLatticeTilingBound(tiling, 2, FirstKnotDerivBound, sizeof( info_first_knot_lower ), &info_first_knot_upper, &info_first_knot_lower) == XLAL_SUCCESS, XLAL_EFAILED);
  XLALSetLatticeTilingPadding(tiling, 2, gsl_vector_get(bboxpad, 2), gsl_vector_get(bboxpad, 2), gsl_vector_int_get(intpad, 2), gsl_vector_int_get(intpad, 2));
  
  ///< Setting the bounds for all following knots
  
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
    
    info_knot_lower.reset = info_knot_upper.reset = reset_value;
    
    int dimindex = 3 * knot;
    
    XLAL_CHECK(XLALSetLatticeTilingBound(tiling, dimindex, F0Bound, sizeof( info_knot_lower ), &info_knot_lower, &info_knot_upper) == XLAL_SUCCESS, XLAL_EFAILED);
    XLALSetLatticeTilingPadding(tiling, dimindex, gsl_vector_get(bboxpad, dimindex), gsl_vector_get(bboxpad, dimindex), gsl_vector_int_get(intpad, dimindex), gsl_vector_int_get(intpad, dimindex));
    
    XLAL_CHECK(XLALSetLatticeTilingBound(tiling, dimindex + 1, F1Bound, sizeof( info_knot_lower ), &info_knot_lower, &info_knot_upper) == XLAL_SUCCESS, XLAL_EFAILED);
    XLALSetLatticeTilingPadding(tiling, dimindex + 1, gsl_vector_get(bboxpad, dimindex + 1), gsl_vector_get(bboxpad, dimindex), gsl_vector_int_get(intpad, dimindex + 1), gsl_vector_int_get(intpad, dimindex + 1));
  
    XLAL_CHECK(XLALSetLatticeTilingBound(tiling, dimindex + 2, F2Bound, sizeof( info_knot_lower ), &info_knot_lower, &info_knot_upper) == XLAL_SUCCESS, XLAL_EFAILED);
    XLALSetLatticeTilingPadding(tiling, dimindex + 2, gsl_vector_get(bboxpad, dimindex + 2), gsl_vector_get(bboxpad, dimindex + 2), gsl_vector_int_get(intpad, dimindex + 2), gsl_vector_int_get(intpad, dimindex + 2));
    
  }
  
  return XLAL_SUCCESS;
}






/// Methods below are for when S = 2








/*

///
/// The LambertW function. Taken from: https://stackoverflow.com/questions/60211021/lambert-w-function-in-c-sharp
///
static double LambertW(
  double x
  )
{
  
  // LambertW is not defined in this section
  XLAL_CHECK(x > 0, XLAL_EINVAL, "LambertW function not defined for x < 0. Given x: %f \n", x);
  
  // computes the first branch for real values only
  
  // amount of iterations (empirically found)
  int amountOfIterations = (int) ceill(log10(x) / 3);
  
  // initial guess is based on 0 < ln(a) < 3
  double w = 3 * log(x + 1) / 4;
  
  // Halley's method via eqn (5.9) in Corless et al (1996)
  for (int i = 0; i < amountOfIterations; i++)
    w = w - (w * exp(w) - x) / (exp(w) * (w + 1) - (w + 2) * (w * exp(w) - x) / (2 * w + 2));
  
  return w;
}
*/

///
/// Function which determines the minimum/maximum value of the first derivative parameters
///
static double F1BoundMinMaxS2(
  double f0,               ///< Frequency parameter at current knot
  double fprev,            ///< Frequency parameter on the previous knot
  double f1prev,           ///< Frequency derivative parameter at previous knot
  double na,               ///< A braking index value. For calculating upper bound, na < nb. For calculating lower bound na > nb
  double nb,               ///< A braking index value. For calculating upper bound, na < nb. For calculating lower bound na > nb
  double ka,               ///< A k value. For calculating upper bound, ka < kb. For calculating lower bound ka > kb
  double kb UNUSED,        ///< A k value. For calculating upper bound, ka < kb. For calculating lower bound ka > kb
  double seglength,        ///< Time difference between this knot and the previous
  int minmax               ///< +1 for upper bound, -1 for lower bound
  )
{

  XLAL_CHECK((na >= nb && minmax == -1) || (na <= nb && minmax == 1), XLAL_EINVAL, "F1BoundMinMaxS2 being called incorrectly, with [na, nb, minmax] = [%f, %f, %d]", na, nb, minmax);
  
  ///< The two parameter range conditions
  double prangecondition1 = -ka * pow(f0, na);
  double prangecondition2 = GTEAndDerivs(fprev, na, ka, seglength, 1);
  
  double bindexcondition = f0 * f1prev / (fprev - (nb - 1) * seglength * f1prev);
  
  ///double lambert_w_argument = - (f0 * fprev * f1prev * log(f0)) / (pow(seglength, 2) * kb);
  ///double kcondition  = f0 * log(f0) / (seglength * LambertW(lambert_w_argument));
  
  
  gsl_vector* vals = gsl_vector_alloc(3);
  gsl_vector_set(vals, 0, prangecondition1);
  gsl_vector_set(vals, 1, prangecondition2);
  gsl_vector_set(vals, 2, bindexcondition);
  ///gsl_vector_set(vals, 3, kcondition);
  
  
  double rtn = NAN;
  printf("Variables used in F1BoundMinMaxS2 are, %f, %f, %f, %f, %f, %E \n", f0, fprev, f1prev, na, nb, ka);
  if(minmax == 1){
    
    printf("F1 Max bounds are: ");
    printvector(vals);
    
    double max = gsl_vector_min(vals);
    rtn = max;
  }
  else if(minmax == -1){
    double min = gsl_vector_max(vals);
    
    printf("F1 Min bounds are: ");
    printvector(vals);
    
    rtn = min;
  }
  
  gsl_vector_free(vals);
  return rtn; 
}

///
/// Defines the bounds for the first and second derivative frequency parameters on the second knot (t = 0).
///
static double SecondKnotBoundS2(
  const void *data,
  const size_t dim,
  const gsl_matrix *cache UNUSED,
  const gsl_vector* pointorig
  )
{
  size_t vectorlength = pointorig->size;
  gsl_vector* point = gsl_vector_alloc(vectorlength);
  gsl_vector_memcpy(point, pointorig);
  
  const PiecewiseBoundInfo* info = (const PiecewiseBoundInfo*)data;
  
  double fmin = info->fmin;
  double fmax = info->fmax;
  double nmin = info->nmin;
  double nmax = info->nmax;
  double kmin = info->kmin;
  double kmax = info->kmax;
  double segmentlength = info->segmentlength;
  int upperlower = info->upperlower;
  int reset = info->reset;
  
  double f0prev = gsl_vector_get(point, 0);
  double f1prev = gsl_vector_get(point, 1);
  
  if (reset == 1){
    f0prev = resetinsidebounds(gsl_vector_get(point, 0), fmin, fmax);
    f1prev = resetinsidebounds(gsl_vector_get(point, 1), -kmax * pow(f0prev, nmax), -kmin * pow(f0prev, nmin));
    
    gsl_vector_set(point, 0, f0prev);
    gsl_vector_set(point, 1, f1prev); 
  }
  
  double rtn = NAN;
  
  if (dim == 2){
    
    if (upperlower == 1){
      double f0 = F0BoundMinMax(f0prev, nmax, nmin, kmax, kmin, segmentlength, upperlower);
      rtn = f0;
    }
    else if (upperlower == -1){
      double f0 = F0BoundMinMax(f0prev, nmin, nmax, kmin, kmax, segmentlength, upperlower);
      rtn = f0;
    }
  }
  else if (dim == 3){
    
    double f0 = gsl_vector_get(point, 2);
    
    if (reset == 1){
      
      double lower = F0BoundMinMax(f0prev, nmin, nmax, kmin, kmax, segmentlength, -1);
      double upper = F0BoundMinMax(f0prev, nmax, nmin, kmax, kmin, segmentlength,  1); 
      
      f0 = resetinsidebounds(gsl_vector_get(point, 2), lower, upper);
      gsl_vector_set(point, 2, f0);
    }
    
    printf("F1BoundMinMaxS2 from SecondKnotBoundS2 \n");
    
    if (upperlower == 1){
      double f1 = F1BoundMinMaxS2(f0, f0prev, f1prev, nmin, nmax, kmin, kmax, segmentlength, upperlower);
      rtn = f1;
    }
    else if (upperlower == -1){
      double f1 = F1BoundMinMaxS2(f0, f0prev, f1prev, nmax, nmin, kmax, kmin, segmentlength, upperlower);
      rtn = f1;
    }
  }
  
  gsl_vector_free(point);
  return rtn;
}

///
/// A method which will reset a given point to be within our parameter space. Method not used if LatticeTilingPaddingFlags = LATTICE_TILING_PADDING_FLAGS_NONE
///
static void resetdimonpointS2(
  gsl_vector* point,       ///< The point which we are resetting to be within the parameter space
  int dim,                 ///< The dimension which we are resetting
  double fmin,             ///< Minimum initial frequency
  double fmax,             ///< Maximum initial frequency
  double nmin,             ///< Minimum allowed braking index
  double nmax,             ///< Maximum allowed braking index
  double ntol,             ///< Braking index tolerance (percentage per second)
  double kmin,             ///< Minimum allowed k value
  double kmax,             ///< Maximum allowed k value
  double ktol,             ///< k value tolerance (percentage per second)
  double segmentlength     ///< The length of the segment we are working on
  )
{
  if (dim == 0){
  
    double f0  = gsl_vector_get(point, 0);
    double val = resetinsidebounds(f0, fmin, fmax);
    gsl_vector_set(point, dim, val);
    return;
  }
  
  else if (dim == 1){
  
    double f0 = gsl_vector_get(point, 0);
    double f1 = gsl_vector_get(point, 1);
    
    double lower = -kmax * pow(f0, nmax);
    double upper = -kmin * pow(f0, nmin);
    
    double val = resetinsidebounds(f1, lower, upper);
    gsl_vector_set(point, dim, val);
    return;
  }
  
  else if (dim == 2){
    
    double f0n1 = gsl_vector_get(point, 0);
    double f0   = gsl_vector_get(point, 2);
    
    double lower = F0BoundMinMax(f0n1, nmin, nmax, kmin, kmax, segmentlength, -1);
    double upper = F0BoundMinMax(f0n1, nmax, nmin, kmax, kmin, segmentlength,  1);
    
    double val = resetinsidebounds(f0, lower, upper);
    gsl_vector_set(point, dim, val);
    return;
  }
  
  else if (dim == 3){
  
    double f0n1 = gsl_vector_get(point, 0);
    double f1n1 = gsl_vector_get(point, 1);
    double f0   = gsl_vector_get(point, 2);
    double f1   = gsl_vector_get(point, 3);
    
    printf("F1BoundMinMaxS2 from resetdimonpointS2 second knot\n");
    
    double lower = F1BoundMinMaxS2(f0, f0n1, f1n1, nmax, nmin, kmax, kmin, segmentlength, -1);
    double upper = F1BoundMinMaxS2(f0, f0n1, f1n1, nmin, nmax, kmin, kmax, segmentlength,  1);
    
    double val = resetinsidebounds(f1, lower, upper);
    gsl_vector_set(point, dim, val);
    return;
  
  }
  else if (dim % 2 == 0){
    double f0nn1 = gsl_vector_get(point, dim - 4);
    double f1nn1 = gsl_vector_get(point, dim - 3);
  
    double f0n1 = gsl_vector_get(point, dim - 2);
    double f1n1 = gsl_vector_get(point, dim - 1);
    
    double f0   = gsl_vector_get(point, dim);
  
    double nprev = 1 + (f0nn1 * f1n1 - f0n1 * f1nn1) / (f1nn1 * f1n1 * segmentlength);
    double kprev = - f1n1 / pow(f0n1, nprev);
  
    nprev = resetvalwithintol(nprev, nmin, nmax, ntol);
    kprev = resetvalwithintol(kprev, kmin, kmax, ktol);
  
    double nminmax[2];
    double kminmax[2];
    
    printf("Running KMinMax in resetdimonpointS2 \n");
    NMinMax(nminmax, nprev, 0.001, nmin, nmax, segmentlength);
    KMinMax(kminmax, kprev, 0.001, kmin, kmax, segmentlength);
    printf("Completed KMinMax in resetdimonpointS2 \n");
    
    double lower = F0BoundMinMax(f0n1, nminmax[0], nminmax[1], kminmax[0], kminmax[1], segmentlength, -1);
    double upper = F0BoundMinMax(f0n1, nminmax[1], nminmax[0], kminmax[1], kminmax[0], segmentlength,  1);

    double val = resetinsidebounds(f0, lower, upper);
    gsl_vector_set(point, dim, val);
    return;
    
  }
  else if (dim % 2 == 1){
    double f0nn1 = gsl_vector_get(point, dim - 5);
    double f1nn1 = gsl_vector_get(point, dim - 4);
    
    double f0n1  = gsl_vector_get(point, dim - 3);
    double f1n1  = gsl_vector_get(point, dim - 2);
    
    double f0   = gsl_vector_get(point, dim - 1);
    double f1   = gsl_vector_get(point, dim);
  
    double nprev = 1 + (f0nn1 * f1n1 - f0n1 * f1nn1) / (f1nn1 * f1n1 * segmentlength);
    double kprev = - f1n1 / pow(f0n1, nprev);
    
    nprev = resetvalwithintol(nprev, nmin, nmax, ntol);
    kprev = resetvalwithintol(kprev, kmin, kmax, ktol);
    
    double nminmax[2];
    double kminmax[2];
    
    printf("Running KMinMax in rest points \n");
    NMinMax(nminmax, nprev, 0.001, nmin, nmax, segmentlength);
    KMinMax(kminmax, kprev, 0.001, kmin, kmax, segmentlength);
    printf("Completed KMinMax in rest points \n");
    
    printf("F1BoundMinMaxS2 from resetdimonpointS2 \n");
    
    double lower = F1BoundMinMaxS2(f0, f0n1, f1n1, nminmax[1], nminmax[0], kminmax[1], kminmax[0], segmentlength, -1);
    double upper = F1BoundMinMaxS2(f0, f0n1, f1n1, nminmax[0], nminmax[1], kminmax[0], kminmax[1], segmentlength,  1);
    
    double val = resetinsidebounds(f1, lower, upper);
    gsl_vector_set(point, dim, val);
    return;
  }
}

///
/// Resets a point such that it is within the bounds of our parameter space. Method not used if LatticeTilingPaddingFlags = LATTICE_TILING_PADDING_FLAGS_NONE
///
static void resetoutofboundspointS2(
  gsl_vector* point,       ///< The point which we are resetting to be within the parameter space
  double fmin,             ///< Minimum initial frequency
  double fmax,             ///< Maximum initial frequency
  double nmin,             ///< Minimum allowed braking index
  double nmax,             ///< Maximum allowed braking index
  double ntol,             ///< Braking index tolerance (percentage per second)
  double kmin,             ///< Minimum allowed k value
  double kmax,             ///< Maximum allowed k value
  double ktol,             ///< k value tolerance (percentage per second)
  double segmentlength     ///< The length of the segment we are working on
  )
{
  printvector(point);
  int dim = point->size;
  //printvector(point);
  for (int i = 0; i < dim; ++i){
    resetdimonpointS2(point, i, fmin, fmax, nmin, nmax, ntol, kmin, kmax, ktol, segmentlength);
  }
  printvector(point);
  //printvector(point);
  //printf("\n");
}

///
/// Sets the bound on the frequency parameter
///
static double F0BoundS2(
  const void *data,
  const size_t dim,
  const gsl_matrix *cache UNUSED,
  const gsl_vector* pointorig
  )
{
  XLAL_CHECK(dim >= 4, XLAL_EINVAL, "F0BoundS2 being called before there is an appropriate number of knots. Minimum dimension required: 4. Current dimension: %zu \n", dim);
  
  size_t vectorlength = pointorig->size;
  gsl_vector* point = gsl_vector_alloc(vectorlength);
  gsl_vector_memcpy(point, pointorig);
  
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
  int reset = info->reset;
  
  if (reset == 1){
    resetoutofboundspointS2(point, fmin, fmax, nmin, nmax, ntol, kmin, kmax, ktol, segmentlength);
  }
  
  double f0nn1 = gsl_vector_get(point, dim - 4);
  double f1nn1 = gsl_vector_get(point, dim - 3);
  
  double f0n1  = gsl_vector_get(point, dim - 2);
  double f1n1  = gsl_vector_get(point, dim - 1);
  
  double nprev = 1 + (f0nn1 * f1n1 - f0n1 * f1nn1) / (f1nn1 * f1n1 * segmentlength);
  double kprev = - f1n1 / pow(f0n1, nprev);
  printf("Nprev and kprev are: %f, %E \n", nprev, kprev);
  nprev = resetvalwithintol(nprev, nmin, nmax, ntol);
  kprev = resetvalwithintol(kprev, kmin, kmax, ktol);
  
  double nminmax[2];
  double kminmax[2];
  
  printf("Nprev and kprev are: %f, %E \n", nprev, kprev);
  
  printf("Running KMinMax in F0BoundS2 \n");
  NMinMax(nminmax, nprev, 0.001, nmin, nmax, segmentlength);
  KMinMax(kminmax, kprev, 0.001, kmin, kmax, segmentlength);
  printf("Completed KMinMax in F0BoundS2 \n");
  
  double rtn = NAN;
  
  if (upperlower == 1){
    double upperbound = F0BoundMinMax(f0n1, nminmax[1], nminmax[0], kminmax[1], kminmax[0], segmentlength, 1);
    rtn = upperbound;
  }
  else if (upperlower == -1){
    double lowerbound = F0BoundMinMax(f0n1, nminmax[0], nminmax[1], kminmax[0], kminmax[1], segmentlength, -1);
    rtn = lowerbound;
  }
  
  gsl_vector_free(point);
  return rtn;
}

///
/// Sets the bound on the first derivative frequency parameter 
///
static double F1BoundS2(
  const void *data,
  const size_t dim,
  const gsl_matrix *cache UNUSED,
  const gsl_vector* pointorig
  )
{
  XLAL_CHECK(dim >= 4, XLAL_EINVAL, "F1BoundS2 being called before there is an appropriate number of knots. Minimum dimension required: 4. Current dimension: %zu \n", dim);
  
  size_t vectorlength = pointorig->size;
  gsl_vector* point = gsl_vector_alloc(vectorlength);
  gsl_vector_memcpy(point, pointorig);
  
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
  int reset = info->reset;
  
  if (reset == 1){
    resetoutofboundspointS2(point, fmin, fmax, nmin, nmax, ntol, kmin, kmax, ktol, segmentlength);
  }
  
  double f0nn1 = gsl_vector_get(point, dim - 5);
  double f1nn1 = gsl_vector_get(point, dim - 4);
  
  double f0n1  = gsl_vector_get(point, dim - 3);
  double f1n1  = gsl_vector_get(point, dim - 2);
  
  double f0    = gsl_vector_get(point, dim - 1);
  
  double nprev = 1 + (f0nn1 * f1n1 - f0n1 * f1nn1) / (f1nn1 * f1n1 * segmentlength);
  double kprev = - f1n1 / pow(f0n1, nprev);
  
  nprev = resetvalwithintol(nprev, nmin, nmax, ntol);
  kprev = resetvalwithintol(kprev, kmin, kmax, ktol);
  
  double nminmax[2];
  double kminmax[2];
  
  printf("Running KMinMax in F1BoundS2 \n");
  NMinMax(nminmax, nprev, 0.001, nmin, nmax, segmentlength);
  KMinMax(kminmax, kprev, 0.001, kmin, kmax, segmentlength);
  printf("Completed KMinMax in F1BoundS2 \n");
  
  double rtn = NAN;
  
  printf("F1BoundMinMaxS2 from F1BoundS2 \n");
  
  if (upperlower == 1){
    double upperbound = F1BoundMinMaxS2(f0, f0n1, f1n1, nminmax[0], nminmax[1], kminmax[0], kminmax[1], segmentlength,  1);
    rtn = upperbound;
  }
  else if (upperlower == -1){
    double lowerbound = F1BoundMinMaxS2(f0, f0n1, f1n1, nminmax[1], nminmax[0], kminmax[1], kminmax[0], segmentlength, -1);
    rtn = lowerbound;
  }
  
  gsl_vector_free(point);
  return rtn;
}

///
/// Sets the bounds for the piecewise model
///
int XLALSetLatticeTilingPiecewiseBoundsS2(
  LatticeTiling* tiling,
  const double fmin,            /// Minimum spin frequency to search over
  const double fmax,            /// Maximum spin frequency to search over
  const double fmaxtrue,        /// Maximum spin frequency used to calculate k and knots (useful for computing tiles in parrallel and fmax != fmaxtrue)
  const double nmin,            /// Minimum braking index
  const double nmax,            /// Maximum braking index
  const double nmin0,           /// Minimum braking index for the first knot. Useful if you want to brake up a search into partitions separated by templates with braking indices within a certain range
  const double nmax0,           /// Maximum braking index for the first knot. Useful if you want to brake up a search into partitions separated by templates with braking indices within a certain range
  const double ntol,            /// Tolerance (percentage per second) between braking indices on adjacent knots
  const double taumin,          /// Minimum spin half life when n = nmax, f0 = fmaxtrue
  const double taumax,          /// Maximum spin half life when n = nmax, f0 = fmaxtrue
  const double ktol,            /// Tolerance (percentage per second) between k values on adjacent knots
  const gsl_vector* knots,      /// List of knots
  const int finalknot,          /// The number of the final knot
  const gsl_vector* bboxpad,    /// Vector containing fractional bounding box padding
  const gsl_vector_int* intpad  /// Vector containing number of integer points to use for padding
  )
{
  ///< Converting tau values to k values
  double kmin = ktauconversion(fmaxtrue, nmax, taumax);
  double kmax = ktauconversion(fmaxtrue, nmax, taumin);
  
  printf("K values are: %E, %E \n", kmin, kmax);
  
  LogPrintf(LOG_DEBUG, "kmin and kmax %E, %E \n", kmin, kmax);
  
  XLAL_CHECK(tiling != NULL,     XLAL_EINVAL);
  XLAL_CHECK(fmin   < fmax,      XLAL_EINVAL, "Bad frequency range: [%f, %f] \n", fmin, fmax);
  XLAL_CHECK(nmin   < nmax,      XLAL_EINVAL, "Bad braking index range: [%f, %f] \n", nmin, nmax);
  XLAL_CHECK(fmax   <= fmaxtrue, XLAL_EINVAL, "fmax larger than fmaxtrue: [%f, %f] \n", fmax, fmaxtrue);
  XLAL_CHECK(nmin0  < nmax0,     XLAL_EINVAL, "Bad braking index knot0 range: [%f, %f] \n", nmin0, nmax0);
  XLAL_CHECK(nmin0  >= nmin,     XLAL_EINVAL, "nmin0 smaller than nmin: [%f, %f] \n", nmin0, nmin);
  XLAL_CHECK(nmax0  <= nmax,     XLAL_EINVAL, "nmax0 greater than nmax: [%f, %f] \n", nmax0, nmax);
  XLAL_CHECK(taumin < taumax,    XLAL_EINVAL, "Bad tau range: [%f, %f] \n", taumin, taumax);
  XLAL_CHECK(kmin   < kmax,      XLAL_EINVAL, "Bad k range: [%f, %f] \n", kmin, kmax);
  
  FirstKnotBoundInfo XLAL_INIT_DECL( info_first_knot_lower );
  FirstKnotBoundInfo XLAL_INIT_DECL( info_first_knot_upper );
  
  info_first_knot_lower.fmin = info_first_knot_upper.fmin = fmin;
  info_first_knot_lower.fmax = info_first_knot_upper.fmax = fmax;
  info_first_knot_lower.nmin = info_first_knot_upper.nmin = nmin0;
  info_first_knot_lower.nmax = info_first_knot_upper.nmax = nmax0;
  info_first_knot_lower.kmin = info_first_knot_upper.kmin = kmin;
  info_first_knot_lower.kmax = info_first_knot_upper.kmax = kmax;
  
  info_first_knot_lower.minmax = -1;
  info_first_knot_upper.minmax =  1;
  
  /// Whether we reset points to be within the parameter space or not. 1 is reset, -1 is not reset
  int reset_value = -1;
  
  info_first_knot_lower.reset = info_first_knot_upper.reset = reset_value;
  
  ///< We only need to use the resetting methods if flags != LATTICE_TILING_PAD_NONE
  ///if (flags == LATTICE_TILING_PAD_NONE){
  ///  info_first_knot_lower.reset = info_first_knot_upper.reset = -1;
  ///}
  
  ///< Setting the bounds on the first knot
  XLALSetLatticeTilingConstantBound(tiling, 0, fmin, fmax);
  XLALSetLatticeTilingPadding(tiling, 0, gsl_vector_get(bboxpad, 0), gsl_vector_get(bboxpad, 0), gsl_vector_int_get(intpad, 0), gsl_vector_int_get(intpad, 0));
  
  XLAL_CHECK(XLALSetLatticeTilingBound(tiling, 1, FirstKnotDerivBound, sizeof( info_first_knot_lower ), &info_first_knot_lower, &info_first_knot_upper) == XLAL_SUCCESS, XLAL_EFAILED);
  XLALSetLatticeTilingPadding(tiling, 1, gsl_vector_get(bboxpad, 1), gsl_vector_get(bboxpad, 1), gsl_vector_int_get(intpad, 1), gsl_vector_int_get(intpad, 1));
  
  
  ///< Setting the bounds for all following knots
  
  for (int knot = 1; knot < finalknot; ++knot){
  
    double segmentlength = gsl_vector_get(knots, knot) - gsl_vector_get(knots, knot - 1);
    
    PiecewiseBoundInfo XLAL_INIT_DECL( info_knot_lower );
    PiecewiseBoundInfo XLAL_INIT_DECL( info_knot_upper );
    
    info_knot_lower.fmin = info_knot_upper.fmin = fmin;
    info_knot_lower.fmax = info_knot_upper.fmax = fmax;
    info_knot_lower.ntol = info_knot_upper.ntol = ntol;
    info_knot_lower.kmin = info_knot_upper.kmin = kmin;
    info_knot_lower.kmax = info_knot_upper.kmax = kmax;
    info_knot_lower.ktol = info_knot_upper.ktol = ktol;
    info_knot_lower.segmentlength = info_knot_upper.segmentlength = segmentlength;
    
    info_knot_lower.upperlower = -1;
    info_knot_upper.upperlower = 1;
    
    info_knot_lower.reset = info_knot_upper.reset = reset_value;
    
    int dimindex = 2 * knot;
    
    if (knot == 1){
      info_knot_lower.nmin = info_knot_upper.nmin = nmin0;
      info_knot_lower.nmax = info_knot_upper.nmax = nmax0;
    
      XLAL_CHECK(XLALSetLatticeTilingBound(tiling, dimindex,     SecondKnotBoundS2, sizeof( info_knot_lower ), &info_knot_lower, &info_knot_upper) == XLAL_SUCCESS, XLAL_EFAILED);
      XLALSetLatticeTilingPadding(tiling, dimindex, gsl_vector_get(bboxpad, dimindex), gsl_vector_get(bboxpad, dimindex), gsl_vector_int_get(intpad, dimindex), gsl_vector_int_get(intpad, dimindex));
      
      printf("F1BoundMinMaxS2 from SecondKnotBoundS2 second knot\n");
      
      XLAL_CHECK(XLALSetLatticeTilingBound(tiling, dimindex + 1, SecondKnotBoundS2, sizeof( info_knot_lower ), &info_knot_lower, &info_knot_upper) == XLAL_SUCCESS, XLAL_EFAILED);
      XLALSetLatticeTilingPadding(tiling, dimindex + 1, gsl_vector_get(bboxpad, dimindex + 1), gsl_vector_get(bboxpad, dimindex + 1), gsl_vector_int_get(intpad, dimindex + 1), gsl_vector_int_get(intpad, dimindex + 1));
    }
    else {
      info_knot_lower.nmin = info_knot_upper.nmin = nmin;
      info_knot_lower.nmax = info_knot_upper.nmax = nmax;
      
      XLAL_CHECK(XLALSetLatticeTilingBound(tiling, dimindex,     F0BoundS2, sizeof( info_knot_lower ), &info_knot_lower, &info_knot_upper) == XLAL_SUCCESS, XLAL_EFAILED);
      XLALSetLatticeTilingPadding(tiling, dimindex, gsl_vector_get(bboxpad, dimindex), gsl_vector_get(bboxpad, dimindex), gsl_vector_int_get(intpad, dimindex), gsl_vector_int_get(intpad, dimindex));
      
      XLAL_CHECK(XLALSetLatticeTilingBound(tiling, dimindex + 1, F1BoundS2, sizeof( info_knot_lower ), &info_knot_lower, &info_knot_upper) == XLAL_SUCCESS, XLAL_EFAILED);
      XLALSetLatticeTilingPadding(tiling, dimindex + 1, gsl_vector_get(bboxpad, dimindex + 1), gsl_vector_get(bboxpad, dimindex + 1), gsl_vector_int_get(intpad, dimindex + 1), gsl_vector_int_get(intpad, dimindex + 1));
    }
    
  }
  
  return XLAL_SUCCESS;
}


