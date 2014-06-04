//
// Copyright (C) 2014 Karl Wette
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

#ifndef _GSLHELPERS_H
#define _GSLHELPERS_H
/// \cond DONT_DOXYGEN

#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#define ALLOC_GSL_VAL(val, name, call) \
  name = (call); \
  XLAL_CHECK_VAL(val, name != NULL, XLAL_ENOMEM, #call " failed")

#define ALLOC_GSL_1D_VAL(val, type, name, n) \
  name = gsl_##type##_calloc(n); \
  XLAL_CHECK_VAL(val, name != NULL, XLAL_ENOMEM, "gsl_"#type"_calloc(%zu) failed", n)

#define ALLOC_GSL_2D_VAL(val, type, name, m, n) \
  name = gsl_##type##_calloc(m, n); \
  XLAL_CHECK_VAL(val, name != NULL, XLAL_ENOMEM, "gsl_"#type"_calloc(%zu,%zu) failed", m, n)

#define PRINT_GSL_1D(type, name, fmt) \
  do { \
    printf("%s:%i ", strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') : __FILE__, __LINE__); \
    printf("%s = [", #name); \
    for (size_t GH_i = 0; name != NULL && GH_i < name->size; ++GH_i) { \
      printf(" "fmt, gsl_##type##_get(name, GH_i)); \
    } \
    printf(" ]\n"); \
  } while (0)

#define PRINT_GSL_2D(type, name, fmt) \
  do { \
    printf("%s:%i ", strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') : __FILE__, __LINE__); \
    printf("%s = [\n", #name); \
    for (size_t GH_i = 0; name != NULL && GH_i < name->size1; ++GH_i) { \
      printf("  "); \
      for (size_t GH_j = 0; GH_j < name->size2; ++GH_j) { \
        printf(" "fmt, gsl_##type##_get(name, GH_i, GH_j)); \
      } \
      printf(";\n"); \
    } \
    printf("]\n"); \
  } while (0)

#define FREE_GSL(type, ...) \
  do { \
    gsl_##type *GH_ptrs[] = {__VA_ARGS__}; \
    for (size_t GH_i = 0; GH_i < sizeof(GH_ptrs)/sizeof(GH_ptrs[0]); ++GH_i) { \
      gsl_##type##_free(GH_ptrs[GH_i]); \
    } \
  } while (0)

#define CALL_GSL_VAL(val, ...)		CALL_GSL_VAL_(val, __VA_ARGS__, NULL, NULL)
#define CALL_GSL_VAL_(val, call, fmt, ...) \
  do { \
    int GH_retn = 0; \
    XLAL_CALLGSL(GH_retn = (call)); \
    if (GH_retn != 0) { \
      if (fmt != NULL) { \
        XLAL_ERROR_VAL(val, XLAL_EFAILED, fmt, __VA_ARGS__, NULL); \
      } else { \
        XLAL_ERROR_VAL(val, XLAL_EFAILED, #call " failed: %s", gsl_strerror(GH_retn)); \
      } \
    } \
  } while (0)

#define GALLOC(name, call)		ALLOC_GSL_VAL(XLAL_FAILURE, name, call)
#define GALLOC_NULL(type, name, n)	ALLOC_GSL_VAL(NULL, name, call)

#define GAVEC(name, n)			ALLOC_GSL_1D_VAL(XLAL_FAILURE, vector, name, n)
#define GAVEC_NULL(name, n)		ALLOC_GSL_1D_VAL(NULL, vector, name, n)
#define GPVEC(name, fmt)		PRINT_GSL_1D(vector, name, fmt)
#define GFVEC(...)			FREE_GSL(vector, __VA_ARGS__)

#define GAVECI(name, n)			ALLOC_GSL_1D_VAL(XLAL_FAILURE, vector_int, name, n)
#define GAVECI_NULL(name, n)		ALLOC_GSL_1D_VAL(NULL, vector_int, name, n)
#define GPVECI(name, fmt)		PRINT_GSL_1D(vector_int, name, fmt)
#define GFVECI(...)			FREE_GSL(vector_int, __VA_ARGS__)

#define GAVECU(name, n)			ALLOC_GSL_1D_VAL(XLAL_FAILURE, vector_uint, name, n)
#define GAVECU_NULL(name, n)		ALLOC_GSL_1D_VAL(NULL, vector_uint, name, n)
#define GPVECU(name, fmt)		PRINT_GSL_1D(vector_uint, name, fmt)
#define GFVECU(...)			FREE_GSL(vector_uint, __VA_ARGS__)

#define GAPERM(name, n)			ALLOC_GSL_1D_VAL(XLAL_FAILURE, permutation, name, n)
#define GAPERM_NULL(name, n)		ALLOC_GSL_1D_VAL(NULL, permutation, name, n)
#define GPPERM(name, fmt)		PRINT_GSL_1D(permutation, name, fmt)
#define GFPERM(...)			FREE_GSL(permutation, __VA_ARGS__)

#define GAMAT(name, m, n)		ALLOC_GSL_2D_VAL(XLAL_FAILURE, matrix, name, m, n)
#define GAMAT_NULL(name, m, n)		ALLOC_GSL_2D_VAL(NULL, matrix, name, m, n)
#define GPMAT(name, fmt)		PRINT_GSL_2D(matrix, name, fmt)
#define GFMAT(...)			FREE_GSL(matrix, __VA_ARGS__)

#define GCALL(...)			CALL_GSL_VAL(XLAL_FAILURE, __VA_ARGS__);
#define GCALL_NULL(...)			CALL_GSL_VAL(NULL, __VA_ARGS__);

/// \endcond
#endif // _GSLHELPERS_H
