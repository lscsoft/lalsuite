#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)

#define SEQUENCETYPE CONCAT2(DATATYPE,Sequence)
#define CSEQUENCE CONCAT2(XLALConjugate,SEQUENCETYPE)

SEQUENCETYPE *CSEQUENCE (
	SEQUENCETYPE *sequence
)
{
	unsigned int i;

	for(i = 0; i < sequence->length; i++)
		sequence->data[i] = CONJ(sequence->data[i]);

	return sequence;
}

#undef SEQUENCETYPE
#undef CSEQUENCE
