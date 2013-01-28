#include <lal/LALInference.h>
#include <gsl/gsl_rng.h>

#define NDIM 9

static void seedRng(gsl_rng *rng) {
  FILE *devRandom = fopen("/dev/urandom", "r");
  unsigned long int seed;

  if (devRandom == NULL) {
    seed = 314159;
  } else {
    if(!fread(&seed, sizeof(unsigned long), 1, devRandom))
    {
        exit(1);
    }
    fclose(devRandom);
  }
  
  gsl_rng_set(rng, seed);
}

static void generatePt(gsl_rng *rng, REAL8 *pt) {
  size_t i;

  for (i = 0; i < NDIM; i++) {
    pt[i] = gsl_rng_uniform(rng);
  }
}

static int inBounds(REAL8 *pt, REAL8 *low, REAL8 *high) {
  size_t i;
  
  for (i = 0; i < NDIM; i++) {
    if (pt[i] < low[i] || pt[i] > high[i]) return 0;
  }

  return 1;
}

static int checkNumberInCell(LALInferenceKDTree *cell) {
  size_t i, count = 0, subCount = 0;

  if (cell == NULL || cell->npts == 0) {
    return 1;
  }

  if (cell->left != NULL) subCount += cell->left->npts;
  if (cell->right != NULL) subCount += cell->right->npts;

  for (i = 0; i < cell->npts; i++) {
    if (inBounds(cell->pts[i], cell->lowerLeft, cell->upperRight)) {
      count++;
    }
  }

  if (cell->npts != 1 && subCount != cell->npts) {
    fprintf(stderr, "Failure: number of points in sub-cells (%zd) is not equal\
 to this cell's number of points (%zd).\n",
            subCount, cell->npts);
    return 0;
  }
  
  if (count != cell->npts) {
    fprintf(stderr, "Failure: number of points actually in cell (%zd) differs\
 from cell's count (%zd).\n",
            count, cell->npts);
    return 0;
  }

  return checkNumberInCell(cell->left) && checkNumberInCell(cell->right);
}

static int checkFindCell(LALInferenceKDTree *tree, gsl_rng *rng) {
  size_t i;

  /* Does findCell work for all the tree points? */
  for (i = 0; i < tree->npts; i++) {
    size_t nCell = (gsl_rng_uniform(rng) < 0.5 ? 0 : gsl_rng_uniform_int(rng, tree->npts));
    REAL8 *pt = tree->pts[i];
    LALInferenceKDTree *theCell = LALInferenceKDFindCell(tree, pt, nCell);

    if (nCell > 1 && theCell->npts >= nCell) {
      fprintf(stderr, "Found non-leaf node that has too many points for requested nCell.\n");
      return 0;
    }

    if (!inBounds(pt, theCell->lowerLeft, theCell->upperRight)) {
      fprintf(stderr, "Searched-for point is not in-bounds for cell returned by findCell.\n");
      return 0;
    }
  }

  /* How about for some points that we make up? */
  for (i = 0; i < 100; i++) {
    size_t nCell = (gsl_rng_uniform(rng) < 0.5 ? 0 : gsl_rng_uniform_int(rng, tree->npts));
    REAL8 pt[NDIM];
    LALInferenceKDTree *theCell;
    size_t j;
    
    for (j = 0; j < NDIM; j++) {
      pt[j] = gsl_rng_uniform(rng);
    }

    theCell = LALInferenceKDFindCell(tree, pt, nCell);
    
    if (theCell == NULL) {
      fprintf(stderr, "Random point not contained in some cell of tree.\n");
      return 0;
    }

    if (nCell > 1 && theCell->npts >= nCell) {
      fprintf(stderr, "Random point's cell is non-leaf and contains too many points for requested nCell.\n");
      return 0;
    }

    if (!inBounds(pt, theCell->lowerLeft, theCell->upperRight)) {
      fprintf(stderr, "Random point not contained in bounds of cell returned by findCell.\n");
      return 0;
    }
  }

  return 1;
}

int main(void);
int main() {
  REAL8 zero[NDIM];
  REAL8 one[NDIM];
  const size_t NPTS = 1000;
  size_t i;
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  LALInferenceKDTree *tree = NULL;

  seedRng(rng);
  
  for (i = 0; i < NDIM; i++) {
    zero[i] = 0.0;
    one[i] = 1.0;
  }

  tree=LALInferenceKDEmpty(zero, one, NDIM);

  for (i = 0; i < NPTS; i++) {
    REAL8 *pt = XLALCalloc(NDIM, sizeof(REAL8));
    generatePt(rng, pt);
    LALInferenceKDAddPoint(tree, pt);
  }

  if (tree->npts != NPTS) return 3;

  if (!checkNumberInCell(tree) && tree->dim == NDIM) {
    return 1;
  }

  if (!checkFindCell(tree, rng)) return 2;

  return 0;
}
