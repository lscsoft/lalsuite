/*
*  Copyright (C) 2007 Chris Messenger
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

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

int WriteMeshFileHeader(FILE *fp,BinaryMeshFileHeader *BMFheader) 
{

  /* this function simply outputs the mesh header to file */

  /* output the header information */
  fprintf(fp,"Search_maximum_search_frequency_Hz                      %6.12f\n",BMFheader->f_max);
  fprintf(fp,"Search_T_span_sec                                       %f\n",BMFheader->tspan);
  fprintf(fp,"Search_Tobs_start_GPS_sec                               %d\n",BMFheader->tstart.gpsSeconds);
  fprintf(fp,"Search_Tobs_start_GPS_nano                              %d\n",BMFheader->tstart.gpsNanoSeconds);
  fprintf(fp,"Search_number_of_filters                                %d\n",BMFheader->Nfilters);
  fprintf(fp,"Search_template_mismatch                                %f\n",BMFheader->mismatch);
  fprintf(fp,"Search_detector                                         %s\n",BMFheader->det);
  fprintf(fp,"Search_right_ascension                                  %f\n",BMFheader->RA);
  fprintf(fp,"Search_declination                                      %f\n",BMFheader->dec);
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
  fprintf(fp,"Search_XY_metric_XX_element                             %f\n",BMFheader->metric_XX);
  fprintf(fp,"Search_XY_metric_XY_element                             %f\n",BMFheader->metric_XY);
  fprintf(fp,"Search_XY_metric_YY_element                             %f\n",BMFheader->metric_YY);
  fprintf(fp,"Search_template_generation_version                      v1\n");
  fprintf(fp,"\n");

  return 0;

}

/*********************************************************************************************/

int ReadMeshFileHeader(FILE *fp,BinaryMeshFileHeader *BMFheader) 
{

  /* this function simply reads in the mesh header to memory */

  char dmp[256];
  int count = 0;
  count += fscanf(fp,"%s%lf",dmp,&BMFheader->f_max);
  count += fscanf(fp,"%s%lf",dmp,&BMFheader->tspan);
  count += fscanf(fp,"%s%d",dmp,&BMFheader->tstart.gpsSeconds);
  count += fscanf(fp,"%s%d",dmp,&BMFheader->tstart.gpsNanoSeconds);
  count += fscanf(fp,"%s%u",dmp,&BMFheader->Nfilters);
  count += fscanf(fp,"%s%lf",dmp,&BMFheader->mismatch);
  count += fscanf(fp,"%s%s",dmp,BMFheader->det);
  count += fscanf(fp,"%s%lf",dmp,&BMFheader->RA);
  count += fscanf(fp,"%s%lf",dmp,&BMFheader->dec);
  count += fscanf(fp,"%s%lf",dmp,&BMFheader->sma_0);
  count += fscanf(fp,"%s%lf",dmp,&BMFheader->sma_MIN);
  count += fscanf(fp,"%s%lf",dmp,&BMFheader->sma_MAX);
  count += fscanf(fp,"%s%d",dmp,&BMFheader->tperi_0.gpsSeconds);
  count += fscanf(fp,"%s%d",dmp,&BMFheader->tperi_0.gpsNanoSeconds);
  count += fscanf(fp,"%s%d",dmp,&BMFheader->tperi_MIN.gpsSeconds);
  count += fscanf(fp,"%s%d",dmp,&BMFheader->tperi_MIN.gpsNanoSeconds);
  count += fscanf(fp,"%s%d",dmp,&BMFheader->tperi_MAX.gpsSeconds);
  count += fscanf(fp,"%s%d",dmp,&BMFheader->tperi_MAX.gpsNanoSeconds);
  count += fscanf(fp,"%s%lf",dmp,&BMFheader->ecc_MIN);
  count += fscanf(fp,"%s%lf",dmp,&BMFheader->ecc_MAX);
  count += fscanf(fp,"%s%lf",dmp,&BMFheader->argp_MIN);
  count += fscanf(fp,"%s%lf",dmp,&BMFheader->argp_MAX);
  count += fscanf(fp,"%s%lf",dmp,&BMFheader->period_MIN);
  count += fscanf(fp,"%s%lf",dmp,&BMFheader->period_MAX);
  count += fscanf(fp,"%s%lf",dmp,&BMFheader->metric_XX);
  count += fscanf(fp,"%s%lf",dmp,&BMFheader->metric_XY);
  count += fscanf(fp,"%s%lf",dmp,&BMFheader->metric_YY);
  count += fscanf(fp,"%s%s",dmp,BMFheader->version);
  count += fscanf(fp," ");

  if ( count != 28 ) fprintf (stderr, "\nSome fscanf() conversions failed in %s\n\n", __func__ );

  return 0;

}

/*********************************************************************************************/

