#include <lal/FileIO.h>
#include <lal/LALStdio.h>
#include <lal/SFTfileIO.h>
#include "sft_extra.h"
#include "sft_extra.h"



/* write a time-series into a text-file */
  void
write_timeSeriesR4 (FILE *fp, const REAL4TimeSeries *series)
{
  REAL8 timestamp; 
  UINT4 i;

  if (series == NULL)
  {
    printf ("\nempty input!\n");
    return; 
  }

  timestamp = 1.0*series->epoch.gpsSeconds + series->epoch.gpsNanoSeconds * 1.0e-9;

  for( i = 0; i < series->data->length; i++)
  {
    fprintf( fp, "%16.9f %e\n", timestamp, series->data->data[i] );
    timestamp += series->deltaT;
  }

  return;

} /* write_timeSeriesR4() */

  void
write_timeSeriesR8 (FILE *fp, const REAL8TimeSeries *series)
{
  REAL8 timestamp; 
  UINT4 i;

  timestamp = 1.0*series->epoch.gpsSeconds + series->epoch.gpsNanoSeconds * 1.0e-9;

  for( i = 0; i < series->data->length; i++)
  {
    fprintf( fp, "%f %e\n", timestamp, series->data->data[i] );
    timestamp += series->deltaT;
  }

  return;

} /* write_timeSeriesR4() */


/* dump an SFT into a text-file 
 *  * format: 0 = openDX (include header), 1 = xmgrace (no header)
 *   */
void dump_SFT (FILE *fp, const SFTtype *sft, INT4 format)
{

  REAL4 valre, valim;
  UINT4 i;
  REAL8 Tsft, freqBand;
  REAL8 f0, df, ff;
  UINT4 nsamples;
  REAL4 norm;
  REAL4 P_k;

  f0 = sft->f0;
  df = sft->deltaF;
  nsamples = sft->data->length;
  Tsft = 1.0 / sft->deltaF;
  freqBand = (nsamples-1.0) * df;

  norm = (REAL4)( Tsft / nsamples);

  /* if openDX format: add a header with number of points..*/
  if ( format == 0)
  {
    fprintf (fp, "points = %d\n", nsamples);
    fprintf (fp, "format = ascii\n");
    fprintf (fp, "field = field0\n");
    fprintf (fp, "structure = 2-vector\n");
    fprintf (fp, "type = float\n");
    fprintf (fp, "dependency = positions\n");
    fprintf (fp, "header = marker \"SFT-data\\n\"\n");
    fprintf (fp, "positions = regular, %f, %f \n", f0, df);
    fprintf(fp, "end\n\n");

    /* write some SFT header-info */
    fprintf (fp, "SFT-header\n");
    fprintf (fp, "Name = %s\n", sft->name);
    fprintf (fp, "Timestamp = %d s, %d ns\n", sft->epoch.gpsSeconds, sft->epoch.gpsNanoSeconds);
    fprintf (fp, "Tsft = %f\n", Tsft);
    fprintf (fp, "Start-frequency = %f Hz\n", f0);
    fprintf (fp, "Frequency-step = %f Hz\n", df);
    fprintf (fp, "Frequency-band = %f Hz\n", freqBand);
    fprintf (fp, "Number of frequency-bins nsamples = %d\n", nsamples);

    /* write SFT-data */
    fprintf (fp, "\nSFT-data\n");
  }

  for (i=0; i < nsamples; i++)
  {
    ff = f0 + i*df;
    valre = norm * sft->data->data[i].re;
    valim = norm * sft->data->data[i].im;
    if ( (i==0) && (i == nsamples-1) )
      P_k = sqrt(valre*valre + valim*valim);
    else
      P_k = 2.0 * sqrt(valre*valre + valim*valim);

    if (format == 1) /* xmgrace */
      fprintf (fp, "%f %e %e %e\n", ff, valre, valim, P_k );
    else if (format == 0) /* openDX */
      fprintf (fp, "%e %e\n", valre, valim);

  } /* for i < nsamples */

  return;

} /* dump_SFT() */
