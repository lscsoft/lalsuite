/*  denest.c: Density estimation from discrete samples. 
    Copyright (C) 2009 Will M. Farr <w-farr@northwestern.edu>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */


#include"denest.h"
#include<assert.h>
#include<stdlib.h>
#include<string.h>

void free_density_tree(tree *t) {
  if (t == NULL) {
    return;
  } else {
    free(t->lower_left);
    free(t->upper_right);
    free_density_tree(t->left);
    free_density_tree(t->right);
    free(t);
    return;
  }
}

static void swap_pts(double **pts, size_t i, size_t j) {
  double *tmp = pts[i];
  pts[i] = pts[j];
  pts[j] = tmp;
}

/* Make sure that the median of the three elements pts[0], pts[n-1],
   pts[n/2] is stored at pts[n/2]. */
static void ensure_median_of_three(double **pts, size_t n, size_t pdim) {
  size_t med_ind = n / 2; 

  if (pts[0][pdim] > pts[med_ind][pdim]) {
    swap_pts(pts, 0, med_ind);
  }

  if (pts[med_ind][pdim] > pts[n-1][pdim]) {
    swap_pts(pts, med_ind, n-1);
  }

  if (pts[0][pdim] > pts[med_ind][pdim]) {
    swap_pts(pts, 0, med_ind);
  }
}

/* Partition the array pts.  On exit, all the points with indices
   smaller than part_ind have smaller coordinates in the pdim
   dimension than the point at pts[part_ind]; all the points with
   larger indices have larger coordinates in this dimension. 

   The algorithm is basically quicksort, but the left- and
   right-sub-arrays need not be sorted by recursive calls.  The
   running time is O(npts). */
static void partition_pts(size_t npts, size_t pdim, double **pts, size_t part_ind) {
  if (npts <= 1) {
    return;
  } else if (npts == 2) {
    if (pts[0][pdim] > pts[1][pdim]) {
      /* Out of order, so swap. */
      swap_pts(pts, 0, 1);
    }
    return;
  } else {
    /* At least three points. */
    size_t low_ind = 0; /* Stores the index of the last point known to
                           be smaller than the pivot element, counting
                           from the start of the array. */
    size_t high_ind = npts - 1; /* Stores the index of the last point
                                   known to be larger than the pivot
                                   element, counting from the end of the array. */
    double *pvt; /* The pivot element. */

    ensure_median_of_three(pts, npts, pdim);

    /* We put the pivot in a known location (pts[1]) because later we
       will re-locate it to the appropriate position. */
    swap_pts(pts, 1, npts/2);
    pvt = pts[1]; 

    do {
      do {
        low_ind++;
      } while (pts[low_ind][pdim] < pvt[pdim]);

      do {
        high_ind--;
      } while (pts[high_ind][pdim] > pvt[pdim]);
      /* At this point, low_ind and high_ind each point to a elements
         that are "out of position" relative to the pivot. */

      /* As long as the low and high pointers have not crossed, swap
         the out-of-position elements. */
      if (low_ind < high_ind) swap_pts(pts, low_ind, high_ind); 

    } while(low_ind < high_ind);

    /* Now we know that pts[high_ind] < pvt.  Since the low and high
       pointers have crossed, we are done with the partitioning.  We
       can now put pvt (stored in pts[1]) in the proper location:
       pts[high_ind]. */
    swap_pts(pts, 1, high_ind);

    /* The pivot element ended up in high_ind; if that is equal to
       part_ind, then we're done.  Otherwise, we should partition one
       of the sub-arrays (so that, in the end, the partition element
       is in part_ind). */
    if (high_ind == part_ind) {
      return;
    } else if (high_ind > part_ind) {
      /* Partition only the beginning of the array again. */
      partition_pts(high_ind, pdim, pts, part_ind);
      return;
    } else {
      /* Partition the end of the array, adjusting part_ind
         appropriately to correspond to the same spot in the full
         array. */
      partition_pts(npts - high_ind, pdim, pts + high_ind, part_ind - high_ind);
      return;
    }
  }
}

static size_t largest_dim(size_t ndim, double *lower_left, double *upper_right) {
  size_t ld = 0, i;
  double max_delta = -1.0/0.0;

  for (i = 0; i < ndim; i++) {
    double delta = upper_right[i] - lower_left[i];

    if (delta > max_delta) {
      ld = i;
      max_delta = delta;
    }
  }

  return ld;
}

/* To split the bounds along pdim when the points are partitioned into
   larger and smaller sets than pts[part_ind], we find the maximum
   coordinate of the small points, and the minimum coordinate of the
   large points.  Then the boundary in pdim is half way between the
   two. */
static void split_bounds(size_t ndim, size_t npts, size_t pdim, size_t part_ind,
                         double **pts, double *lower_left, double *upper_right,
                         double *new_lower_left, double *new_upper_right) {
  size_t i;

  for (i = 0; i < ndim; i++) {
    if (i == pdim) {
      size_t j;
      double max_of_left = -1.0/0.0, min_of_right = 1.0/0.0;
      for (j = 0; j < part_ind; j++) {
        if (pts[j][pdim] > max_of_left) {
          max_of_left = pts[j][pdim];
        }
      }
      for (j = part_ind; j < npts; j++) {
        if (pts[j][pdim] < min_of_right) {
          min_of_right = pts[j][pdim];
        }
      }

      new_lower_left[i] = 0.5*(max_of_left+min_of_right);
      new_upper_right[i] = new_lower_left[i];
    } else {
      new_lower_left[i] = lower_left[i];
      new_upper_right[i] = upper_right[i];
    }
  }
}

tree *make_density_tree(size_t ndim, size_t npts, double **pts, 
                        double *lower_left, double *upper_right) {
  if (npts == 0) {
    /* No tree needed for 0 points. */
    return NULL;
  } else if (npts == 1) {
    /* A singleton tree has both left and right subtrees NULL. */
    tree *t = malloc(sizeof(tree));

    t->ndim = ndim;
    t->left = NULL;
    t->right = NULL;
    t->lower_left = malloc(ndim*sizeof(double));
    t->upper_right = malloc(ndim*sizeof(double));

    memcpy(t->lower_left, lower_left, ndim*sizeof(double));
    memcpy(t->upper_right, upper_right, ndim*sizeof(double));

    return t;
  } else {
    /* At least two points. */
    size_t part_ind = npts / 2;
    size_t pdim = largest_dim(ndim, lower_left, upper_right);
    double *new_lower_left, *new_upper_right;
    tree *t = malloc(sizeof(tree));

    new_lower_left = malloc(ndim*sizeof(double));
    new_upper_right = malloc(ndim*sizeof(double));

    /* Partition points along longest dimension. */
    partition_pts(npts, pdim, pts, part_ind);

    /* Compute the bounds needed for the sub-trees. */
    split_bounds(ndim, npts, pdim, part_ind, pts, lower_left, upper_right,
                 new_lower_left, new_upper_right);

    t->ndim = ndim;

    t->lower_left = malloc(ndim*sizeof(double));
    t->upper_right = malloc(ndim*sizeof(double));

    memcpy(t->lower_left, lower_left, ndim*sizeof(double));
    memcpy(t->upper_right, upper_right, ndim*sizeof(double));

    /* Construct the trees from the left and right elements. */
    t->left = make_density_tree(ndim, part_ind, pts, lower_left, new_upper_right);
    t->right = make_density_tree(ndim, npts-part_ind, pts+part_ind, new_lower_left, 
                                 upper_right);

    free(new_lower_left);
    free(new_upper_right);

    return t;
  }
}

static int in_cell(size_t ndim, double *pt, tree *t) {
  if (t == NULL) {
    return 0;
  } else {
    size_t i;

    for (i = 0; i < ndim; i++) {
      if (pt[i] < t->lower_left[i] || pt[i] > t->upper_right[i]) return 0;
    }

    return 1;
  }
}

static tree *cell_containing(size_t ndim, double *pt, tree *t) {
  assert(t != NULL);

  if (t->left == NULL || t->right == NULL) {
    /* We've got a singleton; make sure that pt is actually contained
       in it. */
    assert(in_cell(ndim, pt, t));
    return t;
  } else if (in_cell(ndim, pt, t->left)) {
    /* pt is in the left sub-tree. */
    return cell_containing(ndim, pt, t->left);
  } else {
    /* pt is in the right sub-tree. */
    return cell_containing(ndim, pt, t->right);
  }
}

void sample_density(size_t ndim, size_t npts, double **pts, tree *t,
                    uniform_random rng, void *rng_data, double *output_pt) {
  size_t i = (size_t) (rng(rng_data)*npts);
  double *pt = pts[i];
  tree *cell = cell_containing(ndim, pt, t);
  size_t j;

  /* Draw uniformly from within the cell that contains point. */
  for (j = 0; j < ndim; j++) {
    output_pt[j] = cell->lower_left[j] + 
      (rng(rng_data))*(cell->upper_right[j] - cell->lower_left[j]);
  }  
}
