static volatile const char *cpu_type_features_h_rcsid = "$Id$";

extern int global_cpu_type;

#if defined(__POWERPC__) && defined(__APPLE__)

#define FEAT_CPU_TYPE_1 FEAT_ALTIVEC
/* #define FEAT_CPU_TYPE_2 FEAT_G5 */

#elif defined(_MSC_VER)

#define FEAT_CPU_TYPE_1 FEAT_SSE

#else /* gcc x86 */

#define FEAT_CPU_TYPE_1 FEAT_SSE

#endif
