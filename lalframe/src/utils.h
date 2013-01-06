#include <lal/LALFrameU.h>
enum { ADC_CHAN_TYPE, PROC_CHAN_TYPE, SIM_CHAN_TYPE };
LALFrameUFrameH *framecpy(LALFrameUFrFile * frfile, size_t pos);
LALFrameUFrameH *framedup(LALFrameUFrameH * frame);
int copydetectors(LALFrameUFrameH * frame, LALFrameUFrFile * frfile);
int copydetector(LALFrameUFrameH * frame, LALFrameUFrDetector * detector);
int copychannels(LALFrameUFrameH * frame, LALFrameUFrFile * frfile, size_t pos, const char *match);
int copychannel(LALFrameUFrameH * frame, LALFrameUFrChan * channel, int chantype);
