/************************************************************************************/
/* The functions below deal with the reading and writing of the header information  */
/* stored at the start of a binary search template file.  The information in the    */
/* header defines all the information used to generate the file.                    */
/*                                                                                  */
/*			           C. Messenger                                     */
/*                                                                                  */
/*                         BIRMINGHAM UNIVERISTY -  2004                            */
/************************************************************************************/

#include "GenerateBinaryMesh_v1.h"

LALStatus status;

int WriteMeshFileHeader(FILE *fp,BinaryMeshFileHeader *BMFheader) 
{

  /* this function simply outputs the mesh header to file */

  /* output the header information */
  fprintf(fp,"Search_maximum_search_frequency_Hz                      %6.12f\n",BMFheader->fmax);
  fprintf(fp,"Search_T_span_sec                                       %lf\n",BMFheader->tspan);
  fprintf(fp,"Search_Tobs_start_GPS_sec                               %d\n",BMFheader->tstart.gpsSeconds);
  fprintf(fp,"Search_Tobs_start_GPS_nano                              %d\n",BMFheader->tstart.gpsNanoSeconds);
  fprintf(fp,"Search_number_of_filters                                %d\n",BMFheader->Nfilters);
  fprintf(fp,"Search_template_mismatch                                %lf\n",BMFheader->mismatch);
  fprintf(fp,"Search_detector                                         %s\n",BMFheader->det);
  fprintf(fp,"Search_right_ascension                                  %lf\n",BMFheader->RA);
  fprintf(fp,"Search_declination                                      %lf\n",BMFheader->dec);
  fprintf(fp,"Source_projected_orbital_semi_major_axis_CENTER_sec     %6.12f\n",BMFheader->sma_0);
  fprintf(fp,"Source_projected_orbital_semi_major_axis_MIN_sec        %6.12f\n",BMFheader->sma_MIN);
  fprintf(fp,"Source_projected_orbital_semi_major_axis_MAX_sec        %6.12f\n",BMFheader->sma_MAX);
  fprintf(fp,"Source_SSB_periapse_passage_CENTER_GPS_sec              %d\n",BMFheader->tperi_0.gpsSeconds);
  fprintf(fp,"Source_SSB_periapse_passage_CENTER_GPS_nanosec          %d\n",BMFheader->tperi_0.gpsNanoSeconds);
  fprintf(fp,"Source_SSB_periapse_passage_MIN_GPS_sec                 %d\n",BMFheader->tperi_MIN.gpsSeconds);
  fprintf(fp,"Source_SSB_periapse_passage_MIN_GPS_nanosec             %d\n",BMFheader->tperi_MIN.gpsNanoSeconds);
  fprintf(fp,"Source_SSB_periapse_passage_MAX_GPS_sec                 %d\n",BMFheader->tperi_MAX.gpsSeconds);
  fprintf(fp,"Source_SSB_periapse_passage_MAX_GPS_nanosec             %d\n",BMFheader->tperi_MAX.gpsNanoSeconds);
  fprintf(fp,"Source_orbital_eccentricity_MIN                         %6.12f\n",BMFheader->ecc_MIN);
  fprintf(fp,"Source_orbital_eccentricity_MAX                         %6.12f\n",BMFheader->ecc_MAX);
  fprintf(fp,"Source_argument_of_periapse_MIN_rad                     %6.12f\n",BMFheader->argp_MIN);
  fprintf(fp,"Source_argument_of_periapse_MAX_rad                     %6.12f\n",BMFheader->argp_MAX);
  fprintf(fp,"Source_orbital_period_MIN_sec                           %6.12f\n",BMFheader->period_MIN);
  fprintf(fp,"Source_orbital_period_MAX_sec                           %6.12f\n",BMFheader->period_MAX);
  fprintf(fp,"Search_XY_metric_XX_element                             %lf\n",BMFheader->metric_XX);
  fprintf(fp,"Search_XY_metric_XY_element                             %lf\n",BMFheader->metric_XY);
  fprintf(fp,"Search_XY_metric_YY_element                             %lf\n",BMFheader->metric_YY);
  fprintf(fp,"Search_template_generation_version                      v1\n");
  fprintf(fp,"\n");

  return 0;

}

/*********************************************************************************************/

int ReadMeshFileHeader(FILE *fp,BinaryMeshFileHeader *BMFheader) 
{

  /* this function simply reads in the mesh header to memory */

  char dmp[256];

  fscanf(fp,"%s%lf",dmp,&BMFheader->fmax);
  fscanf(fp,"%s%lf",dmp,&BMFheader->tspan);
  fscanf(fp,"%s%d",dmp,&BMFheader->tstart.gpsSeconds);
  fscanf(fp,"%s%d",dmp,&BMFheader->tstart.gpsNanoSeconds);
  fscanf(fp,"%s%d",dmp,&BMFheader->Nfilters);
  fscanf(fp,"%s%lf",dmp,&BMFheader->mismatch);
  fscanf(fp,"%s%s",dmp,BMFheader->det);
  fscanf(fp,"%s%lf",dmp,&BMFheader->RA);
  fscanf(fp,"%s%lf",dmp,&BMFheader->dec);
  fscanf(fp,"%s%lf",dmp,&BMFheader->sma_0);
  fscanf(fp,"%s%lf",dmp,&BMFheader->sma_MIN);
  fscanf(fp,"%s%lf",dmp,&BMFheader->sma_MAX);
  fscanf(fp,"%s%d",dmp,&BMFheader->tperi_0.gpsSeconds);
  fscanf(fp,"%s%d",dmp,&BMFheader->tperi_0.gpsNanoSeconds);
  fscanf(fp,"%s%d",dmp,&BMFheader->tperi_MIN.gpsSeconds);
  fscanf(fp,"%s%d",dmp,&BMFheader->tperi_MIN.gpsNanoSeconds);
  fscanf(fp,"%s%d",dmp,&BMFheader->tperi_MAX.gpsSeconds);
  fscanf(fp,"%s%d",dmp,&BMFheader->tperi_MAX.gpsNanoSeconds);
  fscanf(fp,"%s%lf",dmp,&BMFheader->ecc_MIN);
  fscanf(fp,"%s%lf",dmp,&BMFheader->ecc_MAX);
  fscanf(fp,"%s%lf",dmp,&BMFheader->argp_MIN);
  fscanf(fp,"%s%lf",dmp,&BMFheader->argp_MAX);
  fscanf(fp,"%s%lf",dmp,&BMFheader->period_MIN);
  fscanf(fp,"%s%lf",dmp,&BMFheader->period_MAX);
  fscanf(fp,"%s%lf",dmp,&BMFheader->metric_XX);
  fscanf(fp,"%s%lf",dmp,&BMFheader->metric_XY);
  fscanf(fp,"%s%lf",dmp,&BMFheader->metric_YY);
  fscanf(fp,"%s%s",dmp,BMFheader->version);
  fscanf(fp," ");

  return 0;

}

/*********************************************************************************************/

