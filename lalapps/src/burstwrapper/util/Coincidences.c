#include <Coincidences.h>
#include <math.h>

int Get2Coincidences(
		     SnglBurstTableC *output,
		     SnglBurstTableC *input1,
		     SnglBurstTableC *input2,
		     Coincidence2Params *cparams
		     ) {
  /* NOTE: events are ordered by GPS seconds */

  SnglBurstTableC *e1, *e2, *e2ls;

  if(!(input1->next) || !(input2->next)) {
    return 0;
  }

  e1 = input1->next;
  e2ls = input2->next;

  while(e1) {

    while(e2ls->start_time.gpsSeconds + (int)ceil(cparams->dt) < e1->start_time.gpsSeconds) {
      e2ls = e2ls->next;
      
      if(!e2ls) {
	break;
      }
    }

    if(!e2ls) {
      break;
    }

    e2 = e2ls;

    while(e2) {

      int gotc = 1;

      if(e2->start_time.gpsSeconds - e2ls->start_time.gpsSeconds > (int)ceil(cparams->dt)) {
	break;
      }

      if(cparams->dt > 0.0 &&
	 fabs((double)(e1->start_time.gpsSeconds - e2->start_time.gpsSeconds) + 1E-9*(double)(e1->start_time.gpsNanoSeconds - e2->start_time.gpsNanoSeconds)) > cparams->dt) {
	gotc = 0;
      }

      if(cparams->df > 0.0 &&
	 fabs(e1->central_freq - e2->central_freq) > cparams->df) {
	gotc = 0;
      }

      if(cparams->dbw > 0.0 &&
	 fabs(e1->bandwidth - e2->bandwidth) > cparams->dbw) {
	gotc = 0;
      }

      if(cparams->ddur > 0.0 &&
	 fabs(e1->duration - e2->duration) > cparams->ddur) {
	gotc = 0;
      }

      if(cparams->damp > 0.0 &&
	 fabs(e1->amplitude - e2->amplitude) > cparams->damp) {
	gotc = 0;
      }

      if(cparams->dconf > 0.0 &&
	 fabs(e1->confidence - e2->confidence) > cparams->dconf) {
	gotc = 0;
      }

      if(cparams->dsnr > 0.0 &&
	 fabs(e1->snr - e2->snr) > cparams->dsnr) {
	gotc = 0;
      }

      if(gotc) {
	output->next = (SnglBurstTableC *)calloc(1,sizeof(SnglBurstTableC));
	output = output->next;

	output->start_time.gpsSeconds = (e1->start_time.gpsSeconds+e2->start_time.gpsSeconds)/2;
	output->start_time.gpsNanoSeconds = (e1->start_time.gpsNanoSeconds+e2->start_time.gpsNanoSeconds)/2;
	output->central_freq = (e1->central_freq+e2->central_freq)/2.0;
	output->bandwidth = (e1->bandwidth+e2->bandwidth)/2.0;
	output->duration = (e1->duration+e2->duration)/2.0;
	output->confidence = (e1->confidence+e2->confidence)/2.0;
	output->amplitude = (e1->amplitude+e2->amplitude)/2.0;
	output->snr = (e1->snr+e2->snr)/2.0;
      }

      e2 = e2->next;
    }

    e1 = e1->next;
  }

  return 0;
}
