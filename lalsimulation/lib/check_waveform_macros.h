/* Internal utility macro to check all spin components are zero
   returns 1 if all spins zero, otherwise returns 0 */
#define checkSpinsZero(s1x, s1y, s1z, s2x, s2y, s2z) \
    (((s1x) != 0. || (s1y) != 0. || (s1z) != 0. || (s2x) != 0. || (s2y) != 0. || (s2z) != 0.) ? 0 : 1)

/* Internal utility macro to check that the second body's spin components are zero.
   Returns 1 if all components are zero, otherwise returns 0 */
#define checkCOSpinZero(s2x, s2y, s2z) \
    (((s2x) != 0. || (s2y) != 0. || (s2z) != 0.) ? 0 : 1)

/* Internal utility macro to check transverse spins are zero
   returns 1 if x and y components of spins are zero, otherwise returns 0 */
#define checkTransverseSpinsZero(s1x, s1y, s2x, s2y) \
    (((s1x) != 0. || (s1y) != 0. || (s2x) != 0. || (s2y) != 0. ) ? 0 : 1)

/* Internal utility macro to check aligned spins very close to equal
   returns 1 if z components of spins are very close to equal, otherwise returns 0 */
#define checkAlignedSpinsEqual(s1z, s2z) \
    ((fabs((s1z) - (s2z)) > 1e-6) ? 0 : 1)

/* Internal utility macro to check tidal parameters are zero
   returns 1 if both tidal parameters zero, otherwise returns 0 */
#define checkTidesZero(lambda1, lambda2) \
    (((lambda1) != 0. || (lambda2) != 0. ) ? 0 : 1)
