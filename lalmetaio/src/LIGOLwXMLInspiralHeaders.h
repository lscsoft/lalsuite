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
XLALFilePuts( "   <Table Name=\"sngl_inspiral:table\">\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ifo\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"search\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"channel\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time_gmst\" Type=\"real_8\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"impulse_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"impulse_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"template_duration\" Type=\"real_8\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"event_duration\" Type=\"real_8\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amplitude\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_distance\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"coa_phase\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mass1\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mass2\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mchirp\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mtotal\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eta\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"kappa\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chi\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"tau0\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"tau2\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"tau3\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"tau4\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"tau5\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ttotal\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"psi0\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"psi3\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha1\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha2\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha3\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha4\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha5\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha6\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"beta\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"f_final\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"snr\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sigmasq\" Type=\"real_8\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"rsqveto_duration\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"Gamma0\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"Gamma1\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"Gamma2\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"Gamma3\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"Gamma4\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"Gamma5\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"Gamma6\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"Gamma7\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"Gamma8\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"Gamma9\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin1x\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin1y\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin1z\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin2x\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin2y\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin2z\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"event_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Stream Name=\"sngl_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SNGL_INSPIRAL_ROW \
"         \"process:process_id:%ld\",\"%s\",\"%s\",\"%s\",%d,%d,%22.16e,%d,%d,%22.16e,%22.16e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%e,%d,%e,%d,%22.16e,%e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,\"sngl_inspiral:event_id:%ld\""

#define PRINT_LIGOLW_XML_SNGL_INSPIRAL_BNS(fp) ( \
XLALFilePuts( "   <Table Name=\"sngl_inspiral:table\">\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ifo\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"search\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"channel\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time_gmst\" Type=\"real_8\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"template_duration\" Type=\"real_8\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_distance\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"coa_phase\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mass1\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mass2\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mchirp\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mtotal\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eta\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"tau0\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"tau3\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ttotal\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"f_final\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"snr\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sigmasq\" Type=\"real_8\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"rsqveto_duration\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Stream Name=\"sngl_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SNGL_INSPIRAL_ROW_BNS \
"         \"process:process_id:0\",\"%s\",\"%s\",\"%s\",%d,%d,%22.16e,%22.16e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%e,%d,%e,%d,%22.16e,%e"

#define PRINT_LIGOLW_XML_SNGL_INSPIRAL_BCV(fp) ( \
XLALFilePuts( "   <Table Name=\"sngl_inspiral:table\">\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ifo\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"search\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"channel\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time_gmst\" Type=\"real_8\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"template_duration\" Type=\"real_8\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_distance\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"coa_phase\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mchirp\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eta\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"psi0\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"psi3\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"f_final\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"snr\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sigmasq\" Type=\"real_8\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"rsqveto_duration\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Stream Name=\"sngl_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SNGL_INSPIRAL_ROW_BCV \
"         \"process:process_id:0\",\"%s\",\"%s\",\"%s\",%d,%d,%22.16e,%22.16e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%e,%d,%e,%d,%22.16e,%e"

#define PRINT_LIGOLW_XML_MULTI_INSPIRAL(fp) ( \
XLALFilePuts( "   <Table Name=\"multi_inspiral:table\">\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ifos\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"search\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time_gmst\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"impulse_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"impulse_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amplitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"distance\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_h1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_h2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_h1h2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"coa_phase\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mass1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mass2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mchirp\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eta\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chi\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"kappa\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"tau0\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"tau2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"tau3\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"tau4\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"tau5\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ttotal\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"snr\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"snr_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"trace_snr\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"snr_h1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"snr_h2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"snr_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"snr_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"snr_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"snr_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amp_term_1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amp_term_2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amp_term_3\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amp_term_4\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amp_term_5\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amp_term_6\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amp_term_7\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amp_term_8\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amp_term_9\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amp_term_10\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sigmasq_h1\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sigmasq_h2\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sigmasq_l\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sigmasq_g\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sigmasq_t\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sigmasq_v\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq_h1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq_h2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sngl_chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq_h1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq_h2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sngl_bank_chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq_h1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq_h2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sngl_cont_chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ra\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"dec\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ligo_angle\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ligo_angle_sig\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"inclination\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"polarization\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"null_statistic\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"null_stat_h1h2\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"null_stat_degen\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"event_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"h1quad_re\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"h1quad_im\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"h2quad_re\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"h2quad_im\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"l1quad_re\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"l1quad_im\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"g1quad_re\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"g1quad_im\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"t1quad_re\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"t1quad_im\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"v1quad_re\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"v1quad_im\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"coh_snr_h1h2\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cohSnrSqLocal\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"autoCorrCohSq\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"crossCorrCohSq\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"autoCorrNullSq\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"crossCorrNullSq\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ampMetricEigenVal1\"  Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ampMetricEigenVal2\"  Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"time_slide:time_slide_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Stream Name=\"multi_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define MULTI_INSPIRAL_ROW \
  "         \"process:process_id:0\",\"%s\",\"%s\",%d,%d,%22.16e,%d,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%e,%d,%e,%d,%e,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%e,%e,%e,%e,%e,%e,%d,%e,%e,%e,%e,%e,%e,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,\"multi_inspiral:event_id:%" LAL_INT8_FORMAT "\",%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,\"time_slide:time_slide_id:%" LAL_INT8_FORMAT "\""

#define PRINT_LIGOLW_XML_SIM_INSPIRAL(fp) ( \
XLALFilePuts( "   <Table Name=\"sim_inspiral:table\">\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF ||  \
XLALFilePuts( "      <Column Name=\"waveform\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"geocent_end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"geocent_end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"h_end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"h_end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"l_end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"l_end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"g_end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"g_end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"t_end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"t_end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"v_end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"v_end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time_gmst\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"source\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mass1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mass2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mchirp\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eta\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"distance\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"longitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"latitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"inclination\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"coa_phase\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"polarization\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"psi0\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"psi3\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha3\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha4\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha5\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha6\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"beta\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin1x\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin1y\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin1z\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin2x\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin2y\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin2z\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"theta0\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"phi0\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"f_lower\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"f_final\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_h\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"numrel_mode_min\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"numrel_mode_max\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"numrel_data\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amp_order\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"taper\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bandpass\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"simulation_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Stream Name=\"sim_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SIM_INSPIRAL_ROW \
"         \"process:process_id:%ld\",\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%22.16e,\"%s\",%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%d,\"%s\",%d,\"%s\",%d,\"sim_inspiral:simulation_id:%ld\""

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LIGOLWXMLINSPIRALHEADERS_H */
