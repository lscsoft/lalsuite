/*----------------------------------------------------------------------- 
 * 
 * File Name: LIGOLwXMLHeaders.h
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <lal/LALStdlib.h>

#if 0
<lalVerbatim file="LIGOLwXMLHeadersHV">
Author: Brown, D. A.
$Id$
</lalVerbatim> 
<lalLaTeX>
\section{Header \texttt{LIGOLwXMLHeaders.h}}
\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LIGOLwXMLHeaders.h>
\end{verbatim}

This header provides provides \verb|#define|s for the common elements of LIGO
light weight XML files.  It provides the XML header and footer, as well as table
definitions for the various metadata tables.  It will need to be kept up to date
with changes in the LIGO database table definitions.  The quantities which are
defined in this file are 

\begin{itemize}
\item LIGOLW\_XML\_HEADER
\item LIGOLW\_XML\_FOOTER 
\item LIGOLW\_XML\_TABLE\_FOOTER
\item LIGOLW\_XML\_PROCESS 
\item PROCESS\_ROW 
\item LIGOLW\_XML\_PROCESS\_PARAMS 
\item PROCESS\_PARAMS\_ROW 
\item LIGOLW\_XML\_SEARCH\_SUMMARY 
\item SEARCH\_SUMMARY\_ROW 
\item LIGOLW\_XML\_SEARCH\_SUMMVARS 
\item SEARCH\_SUMMVARS\_ROW 
\item LIGOLW\_XML\_SNGL\_BURST 
\item SNGL\_BURST\_ROW 
\item LIGOLW\_XML\_SNGL\_INSPIRAL 
\item SNGL\_INSPIRAL\_ROW 
\item LIGOLW\_XML\_MULTI\_INSPIRAL 
\item MULTI\_INSPIRAL\_ROW
\item LIGOLW\_XML\_SIM\_INSPIRAL
\item SIM\_INSPIRAL\_ROW 
\item LIGOLW\_XML\_SIM\_BURST 
\item SIM\_BURST\_ROW 
\item LIGOLW\_XML\_SUMM\_VALUE 
\item SUMM\_VALUE\_ROW 
\item LIGOLW\_XML\_SIM\_INST\_PARAMS
\item SIM\_INST\_PARAMS\_ROW 
\item LIGOLW\_XML\_STOCHASTIC
\item STOCHASTIC\_ROW
\item LIGOLW\_XML\_STOCH\_SUMM
\item STOCH\_SUMM\_ROW
\item LIGOLW\_XML\_EXT\_TRIGGERS
\item EXT\_TRIGGERS\_ROW
\end{itemize}

\vfill{\footnotesize\input{LIGOLwXMLHeadersHV}}
</lalLaTeX>
#endif

#ifndef _LIGOLWXMLHEADERS_H
#define _LIGOLWXMLHEADERS_H

NRCSID( LIFOLWXMLHEADERSH, "$Id$" );

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

#define PRINT_LIGOLW_XML_HEADER(fp) ( \
fputs( "<?xml version='1.0' encoding='utf-8' ?>\n", fp) == EOF || \
fputs( "<!DOCTYPE LIGO_LW [\n", fp ) == EOF || \
fputs( "<!ELEMENT LIGO_LW ((LIGO_LW|Comment|Param|Table|Array|Stream|IGWDFrame|AdcData|AdcInterval|Time|Detector)*)>\n", fp ) == EOF || \
fputs( "<!ATTLIST LIGO_LW\n", fp ) == EOF || \
fputs( "         Name CDATA #IMPLIED\n", fp ) == EOF || \
fputs( "         Type CDATA #IMPLIED>\n", fp ) == EOF || \
fputs( "\n", fp ) == EOF || \
fputs( "<!ELEMENT Comment (#PCDATA)>\n", fp ) == EOF || \
fputs( "\n", fp ) == EOF || \
fputs( "<!ELEMENT Param (#PCDATA|Comment)*>\n", fp ) == EOF || \
fputs( "<!ATTLIST Param \n", fp ) == EOF || \
fputs( "         Name CDATA #IMPLIED\n", fp ) == EOF || \
fputs( "         Type CDATA #IMPLIED\n", fp ) == EOF || \
fputs( "         Start CDATA #IMPLIED\n", fp ) == EOF || \
fputs( "         Scale CDATA #IMPLIED\n", fp ) == EOF || \
fputs( "         Unit CDATA #IMPLIED\n", fp ) == EOF || \
fputs( "         DataUnit CDATA #IMPLIED>\n", fp ) == EOF || \
fputs( "\n", fp ) == EOF || \
fputs( "<!ELEMENT Table (Comment?,Column*,Stream?)>\n", fp ) == EOF || \
fputs( "<!ATTLIST Table \n", fp ) == EOF || \
fputs( "         Name CDATA #IMPLIED\n", fp ) == EOF || \
fputs( "         Type CDATA #IMPLIED>\n", fp ) == EOF || \
fputs( "\n", fp ) == EOF || \
fputs( "<!ELEMENT Column EMPTY>\n", fp ) == EOF || \
fputs( "<!ATTLIST Column\n", fp ) == EOF || \
fputs( "         Name CDATA #IMPLIED\n", fp ) == EOF || \
fputs( "         Type CDATA #IMPLIED\n", fp ) == EOF || \
fputs( "         Unit CDATA #IMPLIED>\n", fp ) == EOF || \
fputs( "\n", fp ) == EOF || \
fputs( "<!ELEMENT Array (Dim*,Stream?)>\n", fp ) == EOF || \
fputs( "<!ATTLIST Array \n", fp ) == EOF || \
fputs( "         Name CDATA #IMPLIED\n", fp ) == EOF || \
fputs( "         Type CDATA #IMPLIED\n", fp ) == EOF || \
fputs( "         Unit CDATA #IMPLIED>\n", fp ) == EOF || \
fputs( "\n", fp ) == EOF || \
fputs( "<!ELEMENT Dim (#PCDATA)>\n", fp ) == EOF || \
fputs( "<!ATTLIST Dim \n", fp ) == EOF || \
fputs( "         Name  CDATA #IMPLIED\n", fp ) == EOF || \
fputs( "         Unit CDATA #IMPLIED\n", fp ) == EOF || \
fputs( "         Start CDATA #IMPLIED\n", fp ) == EOF || \
fputs( "         Scale CDATA #IMPLIED>\n", fp ) == EOF || \
fputs( "\n", fp ) == EOF || \
fputs( "<!ELEMENT Stream (#PCDATA)>\n", fp ) == EOF || \
fputs( "<!ATTLIST Stream \n", fp ) == EOF || \
fputs( "         Name      CDATA #IMPLIED\n", fp ) == EOF || \
fputs( "         Type      (Remote|Local) \"Local\"\n", fp ) == EOF || \
fputs( "         Delimiter CDATA \",\"\n", fp ) == EOF || \
fputs( "         Encoding  CDATA #IMPLIED\n", fp ) == EOF || \
fputs( "         Content   CDATA #IMPLIED>\n", fp ) == EOF || \
fputs( "\n", fp ) == EOF || \
fputs( "<!ELEMENT IGWDFrame ((Comment|Param|Time|Detector|AdcData|LIGO_LW|Stream?|Array|IGWDFrame)*)>\n", fp ) == EOF || \
fputs( "<!ATTLIST IGWDFrame \n", fp ) == EOF || \
fputs( "         Name CDATA #IMPLIED>\n", fp ) == EOF || \
fputs( "\n", fp ) == EOF || \
fputs( "<!ELEMENT Detector ((Comment|Param|LIGO_LW)*)>\n", fp ) == EOF || \
fputs( "<!ATTLIST Detector \n", fp ) == EOF || \
fputs( "         Name CDATA #IMPLIED>\n", fp ) == EOF || \
fputs( "\n", fp ) == EOF || \
fputs( "<!ELEMENT AdcData ((AdcData|Comment|Param|Time|LIGO_LW|Array)*)>\n", fp ) == EOF || \
fputs( "<!ATTLIST AdcData \n", fp ) == EOF || \
fputs( "         Name CDATA #IMPLIED>\n", fp ) == EOF || \
fputs( "\n", fp ) == EOF || \
fputs( "<!ELEMENT AdcInterval ((AdcData|Comment|Time)*)>\n", fp ) == EOF || \
fputs( "<!ATTLIST AdcInterval \n", fp ) == EOF || \
fputs( "         Name CDATA #IMPLIED\n", fp ) == EOF || \
fputs( "         StartTime CDATA #IMPLIED\n", fp ) == EOF || \
fputs( "         DeltaT CDATA #IMPLIED>\n", fp ) == EOF || \
fputs( "\n", fp ) == EOF || \
fputs( "<!ELEMENT Time (#PCDATA)>\n", fp ) == EOF || \
fputs( "<!ATTLIST Time \n", fp ) == EOF || \
fputs( "         Name CDATA #IMPLIED\n", fp ) == EOF || \
fputs( "         Type (GPS|Unix|ISO-8601) \"ISO-8601\">\n", fp ) == EOF || \
fputs( "]>\n", fp ) == EOF || \
fputs( "\n", fp ) == EOF || \
fputs( "<LIGO_LW>\n", fp ) == EOF || \
fputs( "   <Comment>metadata</Comment>\n", fp ) == EOF )


#define PRINT_LIGOLW_XML_FOOTER(fp) \
fputs( "</LIGO_LW>", fp )

#define PRINT_LIGOLW_XML_TABLE_FOOTER(fp) ( \
fputs( "\n", fp ) == EOF || \
fputs( "      </Stream>\n", fp ) == EOF || \
fputs( "   </Table>\n", fp ) )

#define PRINT_LIGOLW_XML_PROCESS(fp) ( \
fputs( "   <Table Name=\"processgroup:process:table\">\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:program\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:version\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:cvs_repository\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:cvs_entry_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:comment\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:is_online\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:node\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:username\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:unix_procid\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:start_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:end_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:jobid\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:domain\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:ifos\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Stream Name=\"processgroup:process:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define PROCESS_ROW \
"         \"%s\",\"%s\",\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",%d,%d,%d,%d,\"%s\",\"%s\",\"process:process_id:0\""

#define PRINT_LIGOLW_XML_PROCESS_PARAMS(fp) ( \
fputs( "   <Table Name=\"process_paramsgroup:process_params:table\">\n", fp ) == EOF || \
fputs( "      <Column Name=\"process_paramsgroup:process_params:program\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"process_paramsgroup:process_params:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"process_paramsgroup:process_params:param\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"process_paramsgroup:process_params:type\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"process_paramsgroup:process_params:value\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Stream Name=\"process_paramsgroup:process_params:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define PROCESS_PARAMS_ROW \
"         \"%s\",\"process:process_id:0\",\"%s\",\"%s\",\"%s\""

#define PRINT_LIGOLW_XML_SEARCH_SUMMARY(fp) ( \
fputs( "   <Table Name=\"search_summarygroup:search_summary:table\">\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:shared_object\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:lalwrapper_cvs_tag\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:lal_cvs_tag\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:comment\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:in_start_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:in_start_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:in_end_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:in_end_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:out_start_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:out_start_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:out_end_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:out_end_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:nevents\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:nnodes\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Stream Name=\"search_summarygroup:search_summary:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SEARCH_SUMMARY_ROW \
"         \"process:process_id:0\",\"standalone\",\" \",\"%s\",\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%d,%d"

#define PRINT_LIGOLW_XML_SEARCH_SUMMVARS(fp) (  \
fputs( "   <Table Name=\"search_summvarsgroup:search_summvars:table\">\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summvarsgroup:search_summvars:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summvarsgroup:search_summvars:name\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summvarsgroup:search_summvars:string\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summvarsgroup:search_summvars:value\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Stream Name=\"search_summvarsgroup:search_summvars:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SEARCH_SUMMVARS_ROW \
"         \"process:process_id:0\",\"%s\",\"%s\",%22.16e"

#define PRINT_LIGOLW_XML_SNGL_BURST(fp) ( \
fputs( "   <Table Name=\"sngl_burstgroup:sngl_burst:table\">\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_burstgroup:sngl_burst:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_burstgroup:sngl_burst:ifo\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_burstgroup:sngl_burst:search\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_burstgroup:sngl_burst:channel\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_burstgroup:sngl_burst:start_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_burstgroup:sngl_burst:start_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_burstgroup:sngl_burst:peak_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_burstgroup:sngl_burst:peak_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_burstgroup:sngl_burst:duration\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_burstgroup:sngl_burst:central_freq\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_burstgroup:sngl_burst:bandwidth\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_burstgroup:sngl_burst:amplitude\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_burstgroup:sngl_burst:snr\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_burstgroup:sngl_burst:confidence\" Type=\"real_4\"/>\n", fp ) == EOF || \
fputs( "      <Stream Name=\"sngl_burstgroup:sngl_burst:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SNGL_BURST_ROW \
"         \"process:process_id:0\",\"%s\",\"%s\",\"%s\",%d,%d,%d,%d,%e,%e,%e,%e,%e,%e"

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
fputs( "      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:sigmasq\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Stream Name=\"sngl_inspiralgroup:sngl_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SNGL_INSPIRAL_ROW \
"         \"process:process_id:0\",\"%s\",\"%s\",\"%s\",%d,%d,%22.16e,%d,%d,%22.16e,%22.16e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%22.16e"

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
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:ifo1_eff_distance\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:ifo2_eff_distance\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:eff_distance\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:coa_phase\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:mass1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:mass2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:mchirp\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:eta\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:tau0\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:tau2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:tau3\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:tau4\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:tau5\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:ttotal\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:ifo1_snr\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:ifo2_snr\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:snr\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:chisq\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:sigmasq\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:ligo_axis_ra\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:ligo_axis_dec\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:ligo_angle\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"multi_inspiralgroup:multi_inspiral:ligo_angle_sig\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Stream Name=\"multi_inspiralgroup:multi_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define MULTI_INSPIRAL_ROW \
"         \"process:process_id:0\",\"%s\",\"%s\",%d,%d,%22.16e,%d,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%e,%e,%e,%e,%e"

#define PRINT_LIGOLW_XML_SIM_INSPIRAL(fp) ( \
fputs( "   <Table Name=\"sim_inspiralgroup:sim_inspiral:table\">\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:process_id\" Type=\"ilwd:char\"/>\n", fp ) &&  \
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
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:spin1y\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:spin1z\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:spin2x\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:spin2y\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:spin2z\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:theta0\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:phi0\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:f_lower\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:f_final\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:mchirp\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:eff_dist_h\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:eff_dist_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:eff_dist_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:eff_dist_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:eff_dist_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inspiralgroup:sim_inspiral:simulation_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Stream Name=\"sim_inspiralgroup:sim_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SIM_INSPIRAL_ROW \
"         \"process:process_id:0\",\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%22.16e,\"%s\",%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,\"sim_inspiral:simulation_id:0\""

#define PRINT_LIGOLW_XML_SIM_BURST(fp) ( \
fputs( "   <Table Name=\"sim_burstgroup:sim_burst:table\">\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_burstgroup:sim_burst:process_id\" Type=\"ilwd:char\"/>\n", fp ) &&  \
fputs( "      <Column Name=\"sim_burstgroup:sim_burst:waveform\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_burstgroup:sim_burst:geocent_peak_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_burstgroup:sim_burst:geocent_peak_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_burstgroup:sim_burst:h_peak_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_burstgroup:sim_burst:h_peak_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_burstgroup:sim_burst:l_peak_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_burstgroup:sim_burst:l_peak_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_burstgroup:sim_burst:peak_time_gmst\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_burstgroup:sim_burst:dtminus\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_burstgroup:sim_burst:dtplus\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_burstgroup:sim_burst:longitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_burstgroup:sim_burst:latitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_burstgroup:sim_burst:coordinates\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_burstgroup:sim_burst:polarization\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_burstgroup:sim_burst:hrss\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_burstgroup:sim_burst:hpeak\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_burstgroup:sim_burst:freq\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_burstgroup:sim_burst:tau\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_burstgroup:sim_burst:zm_number\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_burstgroup:sim_burst:simulation_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Stream Name=\"sim_burstgroup:sim_burst:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SIM_BURST_ROW \
"         \"process:process_id:0\",\"%s\",%d,%d,%d,%d,%d,%d,%22.16e,%e,%e,%e,%e,\"%s\",%e,%e,%e,%e,%e,%d,\"sim_burst:simulation_id:0\""

#define PRINT_LIGOLW_XML_SUMM_VALUE(fp) ( \
fputs( "   <Table Name=\"summ_valuegroup:summ_value:table\">\n" , fp ) == EOF || \
fputs( "      <Column Name=\"summ_valuegroup:summ_value:program\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"summ_valuegroup:summ_value:process_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"summ_valuegroup:summ_value:start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"summ_valuegroup:summ_value:start_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"summ_valuegroup:summ_value:end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"summ_valuegroup:summ_value:end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"summ_valuegroup:summ_value:ifo\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"summ_valuegroup:summ_value:name\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"summ_valuegroup:summ_value:value\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"summ_valuegroup:summ_value:comment\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Stream Name=\"summ_valuegroup:summ_value:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SUMM_VALUE_ROW \
"         \"%s\",\"process:process_id:0\",%d,%d,%d,%d,\"%s\",\"%s\",%e,\"%s\""

#define PRINT_LIGOLW_XML_SIM_INST_PARAMS(fp) ( \
fputs( "   <Table Name=\"sim_inst_paramsgroup:sim_inst_params:table\">\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inst_paramsgroup:sim_inst_params:simulation_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inst_paramsgroup:sim_inst_params:name\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inst_paramsgroup:sim_inst_params:comment\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inst_paramsgroup:sim_inst_params:value\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Stream Name=\"sim_inst_paramsgroup:sim_inst_params:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SIM_INST_PARAMS_ROW \
"         \"sim_inst:simulation_id:0\",\"%s\",\"%s\",%22.16e"

#define PRINT_LIGOLW_XML_STOCHASTIC(fp) ( \
fputs( "   <Table Name=\"stochasticgroup:stochastic:table\">\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:process_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:ifo_one\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:ifo_two\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:channel_one\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:channel_two\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:start_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:duration\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:duration_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:f_min\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:f_max\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:cc_stat\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:cc_sigma\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Stream Name=\"stochasticgroup:stochastic:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define STOCHASTIC_ROW \
"         \"process:process_id:0\",\"%s\",\"%s\",\"%s\",\"%s\",%d,%d,%d,%d,%.2f,%.2f,%22.16e,%22.16e"

#define PRINT_LIGOLW_XML_STOCH_SUMM(fp) ( \
fputs( "   <Table Name=\"stochsummgroup:stochsumm:table\">\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:process_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:ifo_one\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:ifo_two\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:channel_one\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:channel_two\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:start_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:f_min\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:f_max\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:y_opt\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:error\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Stream Name=\"stochsummgroup:stochsumm:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define STOCH_SUMM_ROW \
"         \"process:process_id:0\",\"%s\",\"%s\",\"%s\",\"%s\",%d,%d,%d,%d,%.2f,%.2f,%22.16e,%22.16e"

#define PRINT_LIGOLW_XML_EXT_TRIGGERS(fp) ( \
fputs( " <Table Name=\"external_trigger:table\">\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:process_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:det_alts\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:det_band\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:det_fluence\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:det_fluence_int\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:det_name\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:det_peak\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:det_peak_int\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:det_snr\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:email_time\" Type=\"int_4s\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:event_dec\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:event_dec_err\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:event_epoch\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:event_err_type\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:event_ra\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:event_ra_err\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:start_time\" Type=\"int_4s\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:start_time_ns\" Type=\"int_4s\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:event_type\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:event_z\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:event_z_err\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:notice_comments\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:notice_id\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:notice_sequence\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:notice_time\" Type=\"int_4s\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:notice_type\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:notice_url\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:obs_fov_dec\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:obs_fov_dec_width\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:obs_fov_ra\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:obs_fov_ra_width\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:obs_loc_ele\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:obs_loc_lat\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:obs_loc_long\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Stream Name=\"external_trigger:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )


#define EXT_TRIGGERS_ROW \
"         \"process:process_id:0\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\",\"%s\",%d, %e, %e, \"%s\",\"%s\", %e , %e, %d, %d, \"%s\",%e, %e, \"%s\",\"%s\",\"%s\",%d, \"%s\",\"%s\",%e, %e,%e,%e, %e, %e, %e"


#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _LIGOLWXMLHEADERS_H */
