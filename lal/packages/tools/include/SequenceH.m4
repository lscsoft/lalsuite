dnl $Id$
define(`SEQUENCETYPE',DATATYPE`Sequence')

void `XLALDestroy'SEQUENCETYPE (
	SEQUENCETYPE *sequence
);

SEQUENCETYPE *`XLALCreate'SEQUENCETYPE (
	size_t length
);

SEQUENCETYPE *`XLALCut'SEQUENCETYPE (
	SEQUENCETYPE *sequence,
	size_t first,
	size_t length
);

SEQUENCETYPE *`XLALCopy'SEQUENCETYPE (
	SEQUENCETYPE *sequence
);

void `XLALShift'SEQUENCETYPE (
	SEQUENCETYPE *sequence,
	int count
);

SEQUENCETYPE *`XLALResize'SEQUENCETYPE (
	SEQUENCETYPE *sequence,
	int first,
	size_t length
);

SEQUENCETYPE *`XLALShrink'SEQUENCETYPE (
	SEQUENCETYPE *sequence,
	size_t first,
	size_t length
);

DATATYPE `XLAL'DATATYPE`Sum' (
	const DATATYPE *data,
	size_t first,
	size_t count
);

SQUAREDATATYPE `XLAL'DATATYPE`SumSquares' (
	const DATATYPE *data,
	size_t first,
	size_t count
);

DATATYPE `XLAL'SEQUENCETYPE`Sum' (
	const SEQUENCETYPE *sequence,
	size_t first,
	size_t count
);

SQUAREDATATYPE `XLAL'SEQUENCETYPE`SumSquares' (
	const SEQUENCETYPE *sequence,
	size_t first,
	size_t count
);

ifelse(DATATYPE, COMPLEX8,
SEQUENCETYPE *`XLALConjugate'SEQUENCETYPE (
        SEQUENCETYPE *series
);
, DATATYPE, COMPLEX16,

SEQUENCETYPE *`XLALConjugate'SEQUENCETYPE (
        SEQUENCETYPE *series
);
,)
