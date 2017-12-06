#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)

#define SERIESTYPE CONCAT2(DATATYPE,FrequencySeries)
#define SEQUENCETYPE CONCAT2(DATATYPE,Sequence)

#define CSERIES CONCAT2(XLALConjugate,SERIESTYPE)
#define CSEQUENCE CONCAT2(XLALConjugate,SEQUENCETYPE)

SERIESTYPE *CSERIES (
	SERIESTYPE *series
)
{
	CSEQUENCE (series->data);
	return series;
}

#undef SERIESTYPE
#undef SEQUENCETYPE

#undef CSERIES
#undef CSEQUENCE
