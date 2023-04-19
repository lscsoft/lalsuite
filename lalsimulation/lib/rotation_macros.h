#include <math.h>
#include <lal/LALConstants.h>


/* Macro functions to rotate the components of a vector about an axis */
#define ROTATEZ(angle, vx, vy, vz)\
	tmp1 = vx*cos(angle) - vy*sin(angle);\
	tmp2 = vx*sin(angle) + vy*cos(angle);\
	vx = tmp1;\
	vy = tmp2

#define ROTATEY(angle, vx, vy, vz)\
	tmp1 = vx*cos(angle) + vz*sin(angle);\
	tmp2 = - vx*sin(angle) + vz*cos(angle);\
	vx = tmp1;\
	vz = tmp2

/* here's the reference explaining why we perform this rotation https://dcc.ligo.org/LIGO-G1900275 */
#define ROT_HP_HC(hp, hc, omega) \
    if ((omega) != 0.0) { \
        REAL8 _c = cos(2.0 * omega); \
        REAL8 _s = sin(2.0 * omega); \
        for (UINT4 _ind = 0; _ind < (hp)->data->length; ++_ind) { \
            REAL8 _x = (hp)->data->data[_ind]; \
            REAL8 _y = (hc)->data->data[_ind]; \
            (hp)->data->data[_ind] = _c * _x + _s * _y; \
            (hc)->data->data[_ind] = _c * _y - _s * _x; \
        } \
    } else ((void)0)
