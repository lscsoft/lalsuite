dnl $Id$
define(`SEQUENCETYPE',DATATYPE`Sequence')
/* <lalLaTeX file="SequenceDestroyP">
\idx{`XLALDestroy'SEQUENCETYPE ()}
</lalLaTeX> <lalVerbatim file="SequenceDestroyP"> */
void `XLALDestroy'SEQUENCETYPE (
	SEQUENCETYPE *sequence
);
/* </lalVerbatim> */

/* <lalLaTeX file="SequenceDestroyP">
\idx{`LALDestroy'SEQUENCETYPE ()}
</lalLaTeX> <lalVerbatim file="SequenceDestroyP"> */
void `LALDestroy'SEQUENCETYPE (
	LALStatus *status,
	SEQUENCETYPE *sequence
);
/* </lalVerbatim> */

/* <lalLaTeX file="SequenceCreateP">
\idx{`XLALCreate'SEQUENCETYPE ()}
</lalLaTeX> <lalVerbatim file="SequenceCreateP"> */
SEQUENCETYPE *`XLALCreate'SEQUENCETYPE (
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="SequenceCreateP">
\idx{`LALCreate'SEQUENCETYPE ()}
</lalLaTeX> <lalVerbatim file="SequenceCreateP"> */
void `LALCreate'SEQUENCETYPE (
	LALStatus *status,
	SEQUENCETYPE **output,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="SequenceCutP">
\idx{`XLALCut'SEQUENCETYPE ()}
</lalLaTeX> <lalVerbatim file="SequenceCutP"> */
SEQUENCETYPE *`XLALCut'SEQUENCETYPE (
	SEQUENCETYPE *sequence,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="SequenceCutP">
\idx{`LALCut'SEQUENCETYPE ()}
</lalLaTeX> <lalVerbatim file="SequenceCutP"> */
void `LALCut'SEQUENCETYPE (
	LALStatus *status,
	SEQUENCETYPE **output,
	SEQUENCETYPE *input,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="SequenceCutP">
\idx{`XLALCopy'SEQUENCETYPE ()}
</lalLaTeX> <lalVerbatim file="SequenceCutP"> */
SEQUENCETYPE *`XLALCopy'SEQUENCETYPE (
	SEQUENCETYPE *sequence
);
/* </lalVerbatim> */

/* <lalLaTeX file="SequenceCutP">
\idx{`LALCopy'SEQUENCETYPE ()}
</lalLaTeX> <lalVerbatim file="SequenceCutP"> */
void `LALCopy'SEQUENCETYPE (
	LALStatus *status,
	SEQUENCETYPE **output,
	SEQUENCETYPE *input
);
/* </lalVerbatim> */

/* <lalLaTeX file="SequenceShiftP">
\idx{`XLALShift'SEQUENCETYPE ()}
</lalLaTeX> <lalVerbatim file="SequenceShiftP"> */
void `XLALShift'SEQUENCETYPE (
	SEQUENCETYPE *sequence,
	int count
);
/* </lalVerbatim> */

/* <lalLaTeX file="SequenceShiftP">
\idx{`LALShift'SEQUENCETYPE ()}
</lalLaTeX> <lalVerbatim file="SequenceShiftP"> */
void `LALShift'SEQUENCETYPE (
	LALStatus *status,
	SEQUENCETYPE *sequence,
	int count
);
/* </lalVerbatim> */

/* <lalLaTeX file="SequenceResizeP">
\idx{`XLALResize'SEQUENCETYPE ()}
</lalLaTeX> <lalVerbatim file="SequenceResizeP"> */
size_t `XLALResize'SEQUENCETYPE (
	SEQUENCETYPE *sequence,
	int first,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="SequenceResizeP">
\idx{`LALResize'SEQUENCETYPE ()}
</lalLaTeX> <lalVerbatim file="SequenceResizeP"> */
void `LALResize'SEQUENCETYPE (
	LALStatus *status,
	SEQUENCETYPE *sequence,
	int first,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="SequenceResizeP">
\idx{`XLALShrink'SEQUENCETYPE ()}
</lalLaTeX> <lalVerbatim file="SequenceResizeP"> */
size_t `XLALShrink'SEQUENCETYPE (
	SEQUENCETYPE *sequence,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="SequenceResizeP">
\idx{`LALShrink'SEQUENCETYPE ()}
</lalLaTeX> <lalVerbatim file="SequenceResizeP"> */
void `LALShrink'SEQUENCETYPE (
	LALStatus *status,
	SEQUENCETYPE *sequence,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

