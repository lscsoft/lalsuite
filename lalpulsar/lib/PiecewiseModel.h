
#include <stdio.h>
#include <math.h>
#include <lal/LatticeTiling.h>


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
  );


int XLALSetLatticeTilingCircleBound(
  LatticeTiling* tiling,
  const double radius
  );

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
  );

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
  );


double GTEAndDerivs(
  double f0,
  double n,
  double k,
  double t,
  int d
  );


double * GetLastElements(
  int initialindex,
  double elems[]
  );


double F0BoundMinMax(
  double f0,
  double na, /* For calculating upper bound, na > nb. For calculating lower bound na < nb */
  double nb, /* For calculating upper bound, na > nb. For calculating lower bound na < nb */
  double ka, /* For calculating upper bound, ka > kb. For calculating lower bound ka < kb */
  double kb, /* For calculating upper bound, ka > kb. For calculating lower bound ka < kb */
  double seglength,
  int minmax
);

double * NMinMax(
  double n,
  double ntol,
  double nmin,
  double nmax,
  double segmentlength UNUSED
  );

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
  );


double F1BoundMinMax(
  double f0,
  double na, /* For upper bound, na < nb */
  double nb
  double ka, /* For upper bound ka < kb. */
  double kb,
  double seglength,
  int minmax
);


double F1Bound(
  const void *data,
  const size_t dim UNUSED
  const gsl_matrix *cache UNUSED,
  const gsl_vector* point
  );


double F2BoundMinMax(
  double f0,
  double f1,
  double n, /* For upper bound n should be maximised */
  double k,
  double seglength,
  int minmax
);

double F2Bound(
  const void *data,
  const size_t dim UNUSED
  const gsl_matrix *cache UNUSED,
  const gsl_vector* point
  );



















