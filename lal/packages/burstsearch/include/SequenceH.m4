dnl $Id$
define(`SEQUENCETYPE',DATATYPE`Sequence')
/* <lalVerbatim file="SequenceDestroyP"> */
void `XLALDestroy'SEQUENCETYPE (
	SEQUENCETYPE *sequence
);
/* </lalVerbatim> */

/* <lalVerbatim file="SequenceDestroyP"> */
void `LALDestroy'SEQUENCETYPE (
	LALStatus *status,
	SEQUENCETYPE *sequence
);
/* </lalVerbatim> */

/* <lalVerbatim file="SequenceCreateP"> */
SEQUENCETYPE *`XLALCreate'SEQUENCETYPE (
	size_t length
);
/* </lalVerbatim> */

/* <lalVerbatim file="SequenceCreateP"> */
void `LALCreate'SEQUENCETYPE (
	LALStatus *status,
	SEQUENCETYPE **output,
	size_t length
);
/* </lalVerbatim> */

/* <lalVerbatim file="SequenceCutP"> */
SEQUENCETYPE *`XLALCut'SEQUENCETYPE (
	SEQUENCETYPE *sequence,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalVerbatim file="SequenceCutP"> */
void `LALCut'SEQUENCETYPE (
	LALStatus *status,
	SEQUENCETYPE **output,
	SEQUENCETYPE *input,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalVerbatim file="SequenceShrinkP"> */
size_t `XLALShrink'SEQUENCETYPE (
	SEQUENCETYPE *sequence,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalVerbatim file="SequenceShrinkP"> */
void `LALShrink'SEQUENCETYPE (
	LALStatus *status,
	SEQUENCETYPE *sequence,
	size_t first,
	size_t length
);
/* </lalVerbatim> */
