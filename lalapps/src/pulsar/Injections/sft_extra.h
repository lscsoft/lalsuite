#include <lal/LALDatatypes.h>
#include <lal/SFTfileIO.h>

void write_timeSeriesR4 (FILE *fp, const REAL4TimeSeries *series);
void write_timeSeriesR8 (FILE *fp, const REAL8TimeSeries *series);
void dump_SFT (FILE *fp, const SFTtype *sft, INT4 format);
