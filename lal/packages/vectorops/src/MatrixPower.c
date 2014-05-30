#include <math.h>
#include <lal/Matrix.h>

#define TYPECODE D
#define TYPE REAL8
#define TYPECODE2 D
#define TYPE2 REAL8
#include "MatrixPower_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE D
#define TYPE REAL8
#define TYPECODE2 S
#define TYPE2 REAL4
#include "MatrixPower_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE D
#define TYPE REAL8
#define TYPECODE2 I2
#define TYPE2 INT2
#include "MatrixPower_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE D
#define TYPE REAL8
#define TYPECODE2 I4
#define TYPE2 INT4
#include "MatrixPower_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE D
#define TYPE REAL8
#define TYPECODE2 I8
#define TYPE2 INT8
#include "MatrixPower_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE S
#define TYPE REAL4
#define TYPECODE2 S
#define TYPE2 REAL4
#include "MatrixPower_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE S
#define TYPE REAL4
#define TYPECODE2 I2
#define TYPE2 INT2
#include "MatrixPower_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE S
#define TYPE REAL4
#define TYPECODE2 I4
#define TYPE2 INT4
#include "MatrixPower_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE S
#define TYPE REAL4
#define TYPECODE2 I8
#define TYPE2 INT8
#include "MatrixPower_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE I2
#define TYPE INT2
#define TYPECODE2 I2
#define TYPE2 INT2
#include "MatrixPower_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE I4
#define TYPE INT4
#define TYPECODE2 I2
#define TYPE2 INT2
#include "MatrixPower_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE I4
#define TYPE INT4
#define TYPECODE2 I4
#define TYPE2 INT4
#include "MatrixPower_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE I8
#define TYPE INT8
#define TYPECODE2 I2
#define TYPE2 INT2
#include "MatrixPower_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE I8
#define TYPE INT8
#define TYPECODE2 I4
#define TYPE2 INT4
#include "MatrixPower_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE I8
#define TYPE INT8
#define TYPECODE2 I8
#define TYPE2 INT8
#include "MatrixPower_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2
