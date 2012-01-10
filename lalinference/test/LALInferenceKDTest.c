#include <../src/LALInference.h>
#include <gsl/gsl_rng.h>

#define NDIM 5

extern void LALInferencePrintCell(LALInferenceKDCell *cell, size_t dim, FILE *stream);
extern void LALInferencePrintKDTree(LALInferenceKDTree *tree, FILE *stream);

static void seedRng(gsl_rng *rng) {
  FILE *devRandom = fopen("/dev/urandom", "r");
  unsigned long int seed;

  if (devRandom == NULL) {
    seed = 314159;
  } else {
    fread(&seed, sizeof(unsigned long), 1, devRandom);
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

static int checkNumberInCell(LALInferenceKDCell *cell, REAL8 **pts, size_t npts) {
  size_t i, count = 0, countPts = 0, subCount = 0;

  if (cell == NULL || cell->npts == 0) {
    return 1;
  }

  if (cell->left != NULL) subCount += cell->left->npts;
  if (cell->right != NULL) subCount += cell->right->npts;

  for (i = 0; i < npts; i++) {
    if (inBounds(pts[i], cell->lowerLeft, cell->upperRight)) {
      count++;
    }

    if (inBounds(pts[i], cell->pointsLowerLeft, cell->pointsUpperRight)) {
      countPts++;
    }
  }

  if (cell->npts != 1 && subCount != cell->npts) {
    fprintf(stderr, "Failure: number of points in sub-cells (%ld) is not equal to this cell's number of points (%ld).\n",
            subCount, cell->npts);
    LALInferencePrintCell(cell, NDIM, stderr);
    return 0;
  }
  
  if (count != countPts) {
    fprintf(stderr, "Failure: number of points in cell bounds (%ld) and points bounds (%ld) doesn't match.\n",
            count, countPts);
    LALInferencePrintCell(cell, NDIM, stderr);
    return 0;
  }

  if (count != cell->npts) {
    fprintf(stderr, "Failure: number of points actually in cell (%ld) differs from cell's count (%ld).\n",
            count, cell->npts);
    LALInferencePrintCell(cell, NDIM, stderr);
    return 0;
  }

  return checkNumberInCell(cell->left, pts, npts) && checkNumberInCell(cell->right, pts, npts);
}

static int checkFindCell(LALInferenceKDTree *tree, gsl_rng *rng) {
  size_t i;

  /* Does findCell work for all the tree points? */
  for (i = 0; i < tree->npts; i++) {
    size_t nCell = (gsl_rng_uniform(rng) < 0.5 ? 0 : gsl_rng_uniform_int(rng, tree->npts));
    REAL8 *pt = tree->pts[i];
    LALInferenceKDCell *theCell = LALInferenceKDFindCell(tree, pt, nCell);

    if (nCell > 1 && theCell->npts >= nCell) return 0;
    if (!inBounds(pt, theCell->lowerLeft, theCell->upperRight)) return 0;
    if (!inBounds(pt, theCell->pointsLowerLeft, theCell->pointsUpperRight)) return 0;
  }

  /* How about for some points that we make up? */
  for (i = 0; i < 100; i++) {
    size_t nCell = (gsl_rng_uniform(rng) < 0.5 ? 0 : gsl_rng_uniform_int(rng, tree->npts));
    REAL8 pt[NDIM];
    LALInferenceKDCell *theCell;
    size_t j;
    
    for (j = 0; j < NDIM; j++) {
      pt[j] = gsl_rng_uniform(rng);
    }

    theCell = LALInferenceKDFindCell(tree, pt, nCell);
    
    if (theCell == NULL) return 0;
    if (nCell > 1 && theCell->npts >= nCell) return 0;
    if (!inBounds(pt, theCell->lowerLeft, theCell->upperRight)) return 0;
  }

  return 1;
}

int main(void);
int main() {
  REAL8 pt[NDIM];
  REAL8 zero[NDIM] = {0.0, 0.0, 0.0, 0.0, 0.0};
  REAL8 one[NDIM] = {1.0, 1.0, 1.0, 1.0, 1.0};
  const size_t NPTS = 100;
  size_t i;
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  LALInferenceKDTree *tree = LALInferenceKDEmpty(zero, one, NDIM);

  seedRng(rng);
  
  for (i = 0; i < NPTS; i++) {
    generatePt(rng, pt);
    LALInferenceKDAddPoint(tree, pt);
  }

  if (!checkNumberInCell(tree->topCell, tree->pts, tree->npts) && tree->ndim == NDIM) {
    return 1;
  }

  if (!checkFindCell(tree, rng)) return 2;

  return 0;
}
