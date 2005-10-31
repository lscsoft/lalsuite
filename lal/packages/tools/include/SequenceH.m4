dnl $Id$
define(`SEQUENCETYPE',DATATYPE`Sequence')
/* <lalLaTeX file="SequenceDestroyP">
\idx{`XLALDestroy'SEQUENCETYPE ()}
</lalLaTeX> <lalVerbatim file="SequenceDestroyP"> */
void `XLALDestroy'SEQUENCETYPE (
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
\idx{`XLALCopy'SEQUENCETYPE ()}
</lalLaTeX> <lalVerbatim file="SequenceCutP"> */
SEQUENCETYPE *`XLALCopy'SEQUENCETYPE (
	SEQUENCETYPE *sequence
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
\idx{`XLALShrink'SEQUENCETYPE ()}
</lalLaTeX> <lalVerbatim file="SequenceResizeP"> */
size_t `XLALShrink'SEQUENCETYPE (
	SEQUENCETYPE *sequence,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="SequenceSumP">
\idx{`XLAL'DATATYPE`Sum' ()}
</lalLaTeX> <lalVerbatim file="SequenceSumP"> */
DATATYPE `XLAL'DATATYPE`Sum' (
	const DATATYPE *data,
	size_t first,
	size_t count
);
/* </lalVerbatim> */

/* <lalLaTeX file="SequenceSumP">
\idx{`XLAL'DATATYPE`SumSquares' ()}
</lalLaTeX> <lalVerbatim file="SequenceSumP"> */
SQUAREDATATYPE `XLAL'DATATYPE`SumSquares' (
	const DATATYPE *data,
	size_t first,
	size_t count
);
/* </lalVerbatim> */

/* <lalLaTeX file="SequenceSumP">
\idx{`XLAL'SEQUENCETYPE`Sum' ()}
</lalLaTeX> <lalVerbatim file="SequenceSumP"> */
DATATYPE `XLAL'SEQUENCETYPE`Sum' (
	const SEQUENCETYPE *sequence,
	size_t first,
	size_t count
);
/* </lalVerbatim> */

/* <lalLaTeX file="SequenceSumP">
\idx{`XLAL'SEQUENCETYPE`SumSquares' ()}
</lalLaTeX> <lalVerbatim file="SequenceSumP"> */
SQUAREDATATYPE `XLAL'SEQUENCETYPE`SumSquares' (
	const SEQUENCETYPE *sequence,
	size_t first,
	size_t count
);
/* </lalVerbatim> */

