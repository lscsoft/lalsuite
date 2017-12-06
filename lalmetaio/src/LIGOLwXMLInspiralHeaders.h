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
 * \ingroup lalmetaio_general
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
fputs( "   <Table Name=\"sngl_inspiral:table\">\n", fp ) == EOF || \
fputs( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"ifo\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"channel\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"end_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"end_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"end_time_gmst\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"impulse_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"impulse_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"template_duration\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"event_duration\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"amplitude\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"eff_distance\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"coa_phase\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"mass1\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"mass2\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"mchirp\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"mtotal\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"eta\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"kappa\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"chi\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"tau0\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"tau2\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"tau3\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"tau4\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"tau5\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"ttotal\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"psi0\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"psi3\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"alpha\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"alpha1\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"alpha2\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"alpha3\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"alpha4\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"alpha5\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"alpha6\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"beta\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"f_final\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"snr\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"bank_chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"bank_chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"cont_chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"cont_chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sigmasq\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"rsqveto_duration\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"Gamma0\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"Gamma1\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"Gamma2\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"Gamma3\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"Gamma4\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"Gamma5\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"Gamma6\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"Gamma7\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"Gamma8\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"Gamma9\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"spin1x\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"spin1y\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"spin1z\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"spin2x\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"spin2y\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"spin2z\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"event_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Stream Name=\"sngl_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SNGL_INSPIRAL_ROW \
"         \"process:process_id:%ld\",\"%s\",\"%s\",\"%s\",%d,%d,%22.16e,%d,%d,%22.16e,%22.16e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%e,%d,%e,%d,%22.16e,%e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,\"sngl_inspiral:event_id:%ld\""

#define PRINT_LIGOLW_XML_SNGL_INSPIRAL_BNS(fp) ( \
fputs( "   <Table Name=\"sngl_inspiral:table\">\n", fp ) == EOF || \
fputs( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"ifo\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"channel\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"end_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"end_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"end_time_gmst\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"template_duration\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"eff_distance\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"coa_phase\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"mass1\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"mass2\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"mchirp\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"mtotal\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"eta\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"tau0\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"tau3\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"ttotal\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"f_final\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"snr\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"bank_chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"bank_chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"cont_chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"cont_chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sigmasq\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"rsqveto_duration\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Stream Name=\"sngl_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SNGL_INSPIRAL_ROW_BNS \
"         \"process:process_id:0\",\"%s\",\"%s\",\"%s\",%d,%d,%22.16e,%22.16e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%e,%d,%e,%d,%22.16e,%e"

#define PRINT_LIGOLW_XML_SNGL_INSPIRAL_BCV(fp) ( \
fputs( "   <Table Name=\"sngl_inspiral:table\">\n", fp ) == EOF || \
fputs( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"ifo\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"channel\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"end_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"end_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"end_time_gmst\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"template_duration\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"eff_distance\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"coa_phase\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"mchirp\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"eta\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"psi0\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"psi3\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"alpha\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"f_final\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"snr\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"bank_chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"bank_chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"cont_chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"cont_chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sigmasq\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"rsqveto_duration\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Stream Name=\"sngl_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SNGL_INSPIRAL_ROW_BCV \
"         \"process:process_id:0\",\"%s\",\"%s\",\"%s\",%d,%d,%22.16e,%22.16e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%e,%d,%e,%d,%22.16e,%e"

#define PRINT_LIGOLW_XML_MULTI_INSPIRAL(fp) ( \
fputs( "   <Table Name=\"multi_inspiral:table\">\n" , fp ) == EOF || \
fputs( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"ifos\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"search\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"end_time_gmst\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"impulse_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"impulse_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"amplitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"distance\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"eff_dist_h1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"eff_dist_h2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"eff_dist_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"eff_dist_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"eff_dist_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"eff_dist_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"eff_dist_h1h2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"coa_phase\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"mass1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"mass2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"mchirp\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"eta\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"chi\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"kappa\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"tau0\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"tau2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"tau3\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"tau4\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"tau5\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"ttotal\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"snr\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"snr_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"chisq\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"bank_chisq\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"bank_chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"cont_chisq\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"cont_chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"trace_snr\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"snr_h1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"snr_h2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"snr_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"snr_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"snr_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"snr_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"amp_term_1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"amp_term_2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"amp_term_3\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"amp_term_4\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"amp_term_5\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"amp_term_6\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"amp_term_7\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"amp_term_8\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"amp_term_9\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"amp_term_10\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sigmasq_h1\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sigmasq_h2\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sigmasq_l\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sigmasq_g\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sigmasq_t\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sigmasq_v\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"chisq_h1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"chisq_h2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"chisq_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"chisq_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"chisq_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"chisq_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sngl_chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"bank_chisq_h1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"bank_chisq_h2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"bank_chisq_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"bank_chisq_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"bank_chisq_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"bank_chisq_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sngl_bank_chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"cont_chisq_h1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"cont_chisq_h2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"cont_chisq_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"cont_chisq_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"cont_chisq_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"cont_chisq_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sngl_cont_chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"ra\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"dec\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"ligo_angle\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"ligo_angle_sig\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"inclination\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"polarization\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"null_statistic\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"null_stat_h1h2\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"null_stat_degen\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"event_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"h1quad_re\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"h1quad_im\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"h2quad_re\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"h2quad_im\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"l1quad_re\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"l1quad_im\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"g1quad_re\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"g1quad_im\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"t1quad_re\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"t1quad_im\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"v1quad_re\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"v1quad_im\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"coh_snr_h1h2\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"cohSnrSqLocal\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"autoCorrCohSq\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"crossCorrCohSq\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"autoCorrNullSq\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"crossCorrNullSq\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"ampMetricEigenVal1\"  Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"ampMetricEigenVal2\"  Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"time_slide:time_slide_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Stream Name=\"multi_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define MULTI_INSPIRAL_ROW \
  "         \"process:process_id:0\",\"%s\",\"%s\",%d,%d,%22.16e,%d,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%e,%d,%e,%d,%e,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%e,%e,%e,%e,%e,%e,%d,%e,%e,%e,%e,%e,%e,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,\"multi_inspiral:event_id:%" LAL_INT8_FORMAT "\",%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,\"time_slide:time_slide_id:%" LAL_INT8_FORMAT "\""

#define PRINT_LIGOLW_XML_SIM_INSPIRAL(fp) ( \
fputs( "   <Table Name=\"sim_inspiral:table\">\n" , fp ) == EOF || \
fputs( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF ||  \
fputs( "      <Column Name=\"waveform\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"geocent_end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"geocent_end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"h_end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"h_end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"l_end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"l_end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"g_end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"g_end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"t_end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"t_end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"v_end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"v_end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"end_time_gmst\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"source\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"mass1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"mass2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"mchirp\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"eta\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"distance\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"longitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"latitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"inclination\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"coa_phase\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"polarization\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"psi0\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"psi3\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"alpha\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"alpha1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"alpha2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"alpha3\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"alpha4\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"alpha5\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"alpha6\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"beta\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"spin1x\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"spin1y\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"spin1z\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"spin2x\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"spin2y\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"spin2z\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"theta0\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"phi0\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"f_lower\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"f_final\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"eff_dist_h\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"eff_dist_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"eff_dist_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"eff_dist_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"eff_dist_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"numrel_mode_min\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"numrel_mode_max\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"numrel_data\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"amp_order\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"taper\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"bandpass\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"simulation_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Stream Name=\"sim_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SIM_INSPIRAL_ROW \
"         \"process:process_id:%ld\",\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%22.16e,\"%s\",%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%d,\"%s\",%d,\"%s\",%d,\"sim_inspiral:simulation_id:%" LAL_INT8_FORMAT "\""

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LIGOLWXMLINSPIRALHEADERS_H */
