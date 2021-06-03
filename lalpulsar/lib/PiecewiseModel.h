
#include <stdio.h>
#include <math.h>
#include <lal/LatticeTiling.h>


#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif


int XLALSetLatticeTilingCircleBound(
  LatticeTiling* tiling,
  const double radius
  );

///
/// Sets the bounds for the piecewise model
///
int XLALSetLatticeTilingPiecewiseBounds(
  LatticeTiling* tiling,
  const double fmin,       /// Minimum spin frequency to search over
  const double fmax,       /// Maximum spin frequency to search over
  const double fmaxtrue,   /// Maximum spin frequency used to calculate k and knots (useful for computing tiles in parrallel and fmax != fmaxtrue)
  const double nmin,       /// Minimum braking index
  const double nmax,       /// Maximum braking index
  const double ntol,       /// Tolerance (percentage per second) between braking indices on adjacent knots
  const double taumin,     /// Minimum spin half life when n = nmax, f0 = fmaxtrue
  const double taumax,     /// Maximum spin half life when n = nmax, f0 = fmaxtrue
  const double ktol,       /// Tolerance (percentage per second) between k values on adjacent knots
  const gsl_vector* knots, /// List of knots
  const int finalknot      /// The number of the final knot
  );
  
///
/// Sets the bounds for the piecewise model when we are using 2 spin down parameters for each knot
///
int XLALSetLatticeTilingPiecewiseBoundsS2(
  LatticeTiling* tiling,
  const double fmin,       /// Minimum spin frequency to search over
  const double fmax,       /// Maximum spin frequency to search over
  const double fmaxtrue,   /// Maximum spin frequency used to calculate k and knots (useful for computing tiles in parrallel and fmax != fmaxtrue)
  const double nmin,       /// Minimum braking index
  const double nmax,       /// Maximum braking index
  const double ntol,       /// Tolerance (percentage per second) between braking indices on adjacent knots
  const double taumin,     /// Minimum spin half life when n = nmax, f0 = fmaxtrue
  const double taumax,     /// Maximum spin half life when n = nmax, f0 = fmaxtrue
  const double ktol,       /// Tolerance (percentage per second) between k values on adjacent knots
  const gsl_vector* knots, /// List of knots
  const int finalknot      /// The number of the final knot
  );

















