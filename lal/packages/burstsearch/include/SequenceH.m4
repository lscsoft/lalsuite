dnl $Id$
define(SEQUENCETYPE,DATATYPE`Sequence')
void `XLALDestroy'SEQUENCETYPE (
	SEQUENCETYPE *sequence
);

void `LALDestroy'SEQUENCETYPE (
	LALStatus *status,
	SEQUENCETYPE *sequence
);

SEQUENCETYPE *`XLALCreate'SEQUENCETYPE (
	size_t length
);

void `LALCreate'SEQUENCETYPE (
	LALStatus *status,
	SEQUENCETYPE **output,
	size_t length
);

SEQUENCETYPE *`XLALCut'SEQUENCETYPE (
	SEQUENCETYPE *sequence,
	size_t first,
	size_t length
);

void `LALCut'SEQUENCETYPE (
	LALStatus *status,
	SEQUENCETYPE **output,
	SEQUENCETYPE *input,
	size_t first,
	size_t length
);

size_t `XLALShrink'SEQUENCETYPE (
	SEQUENCETYPE *sequence,
	size_t first,
	size_t length
);

void `LALShrink'SEQUENCETYPE (
	LALStatus *status,
	SEQUENCETYPE *sequence,
	size_t first,
	size_t length
);
