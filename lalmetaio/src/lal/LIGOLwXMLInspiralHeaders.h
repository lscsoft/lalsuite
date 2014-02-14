/*
*  Copyright (C) 2007 Sukanta Bose, Duncan Brown, Jolien Creighton, Kipp Cannon, Anand Sengupta
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

/*-----------------------------------------------------------------------
 *
 * File Name: LIGOLwXMLInspiralHeaders.h
 *
 * Author: Brown, D. A.
 *
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <lal/LALStdlib.h>

/**
 * \author Brown, D. A.
 * \file
 * \ingroup lalmetaio
 *
 * \brief This header provides provides <tt>\#define</tt>s for the inspiral related
 * tables of LIGO light weight XML files.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/LIGOLwXMLInspiralHeaders.h>
 * \endcode
 *
 * It will need to be kept up to date with changes in the LIGO database table definitions.
 * The quantities which are defined in this file are
 *
 * <ul>
 * <li> LIGOLW_XML_SNGL_INSPIRAL</li>
 * <li> SNGL_INSPIRAL_ROW</li>
 * <li> LIGOLW_XML_MULTI_INSPIRAL</li>
 * <li> MULTI_INSPIRAL_ROW</li>
 * <li> LIGOLW_XML_SIM_INSPIRAL</li>
 * <li> SIM_INSPIRAL_ROW</li>
 * </ul>
 *
 */

#ifndef _LIGOLWXMLINSPIRALHEADERS_H
#define _LIGOLWXMLINSPIRALHEADERS_H

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

#define PRINT_LIGOLW_XML_SNGL_INSPIRAL(fp) ( \
fputs( "   <Table Name=\"sngl_inspiralgroup:sngl_inspiral:table\">\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:ifo\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:search\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:channel\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:end_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:end_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:end_time_gmst\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:impulse_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:impulse_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:template_duration\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:event_duration\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:amplitude\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:eff_distance\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:coa_phase\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:mass1\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:mass2\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:mchirp\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:mtotal\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:eta\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:kappa\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:chi\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:tau0\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:tau2\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:tau3\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:tau4\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:tau5\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:ttotal\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:psi0\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:psi3\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:alpha\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:alpha1\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:alpha2\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:alpha3\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:alpha4\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:alpha5\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:alpha6\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:beta\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:f_final\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:snr\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:bank_chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:bank_chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:cont_chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:cont_chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:sigmasq\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:rsqveto_duration\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:Gamma0\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:Gamma1\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:Gamma2\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:Gamma3\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:Gamma4\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:Gamma5\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:Gamma6\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:Gamma7\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:Gamma8\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:Gamma9\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:spin1x\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:spin1y\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:spin1z\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:spin2x\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:spin2y\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:spin2z\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:event_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Stream Name=\"sngl_inspiralgroup:sngl_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SNGL_INSPIRAL_ROW \
"         \"process:process_id:0\",\"%s\",\"%s\",\"%s\",%d,%d,%22.16e,%d,%d,%22.16e,%22.16e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%e,%d,%e,%d,%22.16e,%e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,\"sngl_inspiral:event_id:%" LAL_INT8_FORMAT "\""

#define PRINT_LIGOLW_XML_SNGL_INSPIRAL_BNS(fp) ( \
fputs( "   <Table Name=\"sngl_inspiralgroup:sngl_inspiral:table\">\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:ifo\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:search\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:channel\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:end_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:end_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:end_time_gmst\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:template_duration\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:eff_distance\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:coa_phase\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:mass1\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:mass2\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:mchirp\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:mtotal\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:eta\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:tau0\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:tau3\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:ttotal\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:f_final\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:snr\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:bank_chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:bank_chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:cont_chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:cont_chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:sigmasq\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:rsqveto_duration\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Stream Name=\"sngl_inspiralgroup:sngl_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SNGL_INSPIRAL_ROW_BNS \
"         \"process:process_id:0\",\"%s\",\"%s\",\"%s\",%d,%d,%22.16e,%22.16e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%e,%d,%e,%d,%22.16e,%e"

#define PRINT_LIGOLW_XML_SNGL_INSPIRAL_BCV(fp) ( \
fputs( "   <Table Name=\"sngl_inspiralgroup:sngl_inspiral:table\">\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:ifo\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:search\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:channel\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:end_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:end_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:end_time_gmst\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:template_duration\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:eff_distance\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:coa_phase\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:mchirp\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:eta\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:psi0\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:psi3\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:alpha\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:f_final\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:snr\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:bank_chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:bank_chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:cont_chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:cont_chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:sigmasq\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:rsqveto_duration\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Stream Name=\"sngl_inspiralgroup:sngl_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SNGL_INSPIRAL_ROW_BCV \
"         \"process:process_id:0\",\"%s\",\"%s\",\"%s\",%d,%d,%22.16e,%22.16e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%e,%d,%e,%d,%22.16e,%e"

#define PRINT_LIGOLW_XML_MULTI_INSPIRAL(fp) ( \
fputs( "   <Table Name=\"multi_inspiralgroup:multi_inspiral:table\">\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:process_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:ifos\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:search\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:end_time_gmst\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:impulse_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:impulse_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:amplitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:distance\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:eff_dist_h1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:eff_dist_h2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:eff_dist_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:eff_dist_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:eff_dist_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:eff_dist_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:eff_dist_h1h2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:coa_phase\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:mass1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:mass2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:mchirp\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:eta\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:chi\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:kappa\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:tau0\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:tau2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:tau3\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:tau4\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:tau5\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:ttotal\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:snr\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:snr_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:chisq\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:bank_chisq\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:bank_chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:cont_chisq\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:cont_chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:trace_snr\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:snr_h1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:snr_h2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:snr_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:snr_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:snr_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:snr_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:amp_term_1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:amp_term_2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:amp_term_3\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:amp_term_4\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:amp_term_5\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:amp_term_6\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:amp_term_7\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:amp_term_8\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:amp_term_9\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:amp_term_10\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:sigmasq_h1\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:sigmasq_h2\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:sigmasq_l\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:sigmasq_g\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:sigmasq_t\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:sigmasq_v\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:chisq_h1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:chisq_h2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:chisq_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:chisq_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:chisq_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:chisq_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:sngl_chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:bank_chisq_h1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:bank_chisq_h2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:bank_chisq_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:bank_chisq_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:bank_chisq_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:bank_chisq_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:sngl_bank_chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:cont_chisq_h1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:cont_chisq_h2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:cont_chisq_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:cont_chisq_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:cont_chisq_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:cont_chisq_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:sngl_cont_chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:ra\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:dec\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:ligo_angle\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:ligo_angle_sig\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:inclination\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:polarization\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:null_statistic\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:null_stat_h1h2\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:null_stat_degen\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:event_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:h1quad_re\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:h1quad_im\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:h2quad_re\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:h2quad_im\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:l1quad_re\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:l1quad_im\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:g1quad_re\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:g1quad_im\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:t1quad_re\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:t1quad_im\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:v1quad_re\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:v1quad_im\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:coh_snr_h1h2\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:cohSnrSqLocal\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:autoCorrCohSq\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:crossCorrCohSq\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:autoCorrNullSq\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:crossCorrNullSq\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:ampMetricEigenVal1\"  Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:ampMetricEigenVal2\"  Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:time_slide_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Stream Name=\"multi_inspiralgroup:multi_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define MULTI_INSPIRAL_ROW \
  "         \"process:process_id:0\",\"%s\",\"%s\",%d,%d,%22.16e,%d,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%e,%d,%e,%d,%e,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%e,%e,%e,%e,%e,%e,%d,%e,%e,%e,%e,%e,%e,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,\"multi_inspiral:event_id:%" LAL_INT8_FORMAT "\",%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,\"time_slide:time_slide_id:%" LAL_INT8_FORMAT "\""

#define PRINT_LIGOLW_XML_SIM_INSPIRAL(fp) ( \
fputs( "   <Table Name=\"sim_inspiralgroup:sim_inspiral:table\">\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF ||  \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:waveform\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:geocent_end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:geocent_end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:h_end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:h_end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:l_end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:l_end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:g_end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:g_end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:t_end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:t_end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:v_end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:v_end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:end_time_gmst\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:source\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:mass1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:mass2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:mchirp\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:eta\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:distance\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:longitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:latitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:inclination\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:coa_phase\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:polarization\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:psi0\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:psi3\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:alpha\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:alpha1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:alpha2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:alpha3\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:alpha4\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:alpha5\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:alpha6\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:beta\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:spin1x\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:spin1y\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:spin1z\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:spin2x\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:spin2y\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:spin2z\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:theta0\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:phi0\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:f_lower\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:f_final\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:eff_dist_h\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:eff_dist_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:eff_dist_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:eff_dist_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:eff_dist_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:numrel_mode_min\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:numrel_mode_max\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:numrel_data\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:amp_order\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:taper\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:bandpass\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:simulation_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Stream Name=\"sim_inspiralgroup:sim_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SIM_INSPIRAL_ROW \
"         \"process:process_id:0\",\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%22.16e,\"%s\",%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%d,\"%s\",%d,\"%s\",%d,\"sim_inspiral:simulation_id:%" LAL_INT8_FORMAT "\""

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LIGOLWXMLINSPIRALHEADERS_H */
