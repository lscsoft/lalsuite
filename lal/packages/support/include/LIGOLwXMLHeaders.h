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
with changes in the LIGO database table definitions.
\vfill{\footnotesize\input{LIGOLwXMLHeadersHV}
</lalLaTeX>
#endif

#ifndef _LIGOLWXMLHEADERS_H
#define _LIGOLWXMLHEADERS_H

NRCSID( LIFOLWXMLHEADERSH, "$Id$" );

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

#define LIGOLW_XML_HEADER \
"<?xml version='1.0' encoding='utf-8' ?>\n" \
"<!DOCTYPE LIGO_LW [\n" \
"<!ELEMENT LIGO_LW ((LIGO_LW|Comment|Param|Table|Array|Stream|IGWDFrame|AdcData|AdcInterval|Time|Detector)*)>\n" \
"<!ATTLIST LIGO_LW\n" \
"         Name CDATA #IMPLIED\n" \
"         Type CDATA #IMPLIED>\n" \
"\n" \
"<!ELEMENT Comment (#PCDATA)>\n" \
"\n" \
"<!ELEMENT Param (#PCDATA|Comment)*>\n" \
"<!ATTLIST Param \n" \
"         Name CDATA #IMPLIED\n" \
"         Type CDATA #IMPLIED\n" \
"         Start CDATA #IMPLIED\n" \
"         Scale CDATA #IMPLIED\n" \
"         Unit CDATA #IMPLIED\n" \
"         DataUnit CDATA #IMPLIED>\n" \
"\n" \
"<!ELEMENT Table (Comment?,Column*,Stream?)>\n" \
"<!ATTLIST Table \n" \
"         Name CDATA #IMPLIED\n" \
"         Type CDATA #IMPLIED>\n" \
"\n" \
"<!ELEMENT Column EMPTY>\n" \
"<!ATTLIST Column\n" \
"         Name CDATA #IMPLIED\n" \
"         Type CDATA #IMPLIED\n" \
"         Unit CDATA #IMPLIED>\n" \
"\n" \
"<!ELEMENT Array (Dim*,Stream?)>\n" \
"<!ATTLIST Array \n" \
"         Name CDATA #IMPLIED\n" \
"         Type CDATA #IMPLIED\n" \
"         Unit CDATA #IMPLIED>\n" \
"\n" \
"<!ELEMENT Dim (#PCDATA)>\n" \
"<!ATTLIST Dim \n" \
"         Name  CDATA #IMPLIED\n" \
"         Unit CDATA #IMPLIED\n" \
"         Start CDATA #IMPLIED\n" \
"         Scale CDATA #IMPLIED>\n" \
"\n" \
"<!ELEMENT Stream (#PCDATA)>\n" \
"<!ATTLIST Stream \n" \
"         Name      CDATA #IMPLIED\n" \
"         Type      (Remote|Local) \"Local\"\n" \
"         Delimiter CDATA \",\"\n" \
"         Encoding  CDATA #IMPLIED\n" \
"         Content   CDATA #IMPLIED>\n" \
"\n" \
"<!ELEMENT IGWDFrame ((Comment|Param|Time|Detector|AdcData|LIGO_LW|Stream?|Array|IGWDFrame)*)>\n" \
"<!ATTLIST IGWDFrame \n" \
"         Name CDATA #IMPLIED>\n" \
"\n" \
"<!ELEMENT Detector ((Comment|Param|LIGO_LW)*)>\n" \
"<!ATTLIST Detector \n" \
"         Name CDATA #IMPLIED>\n" \
"\n" \
"<!ELEMENT AdcData ((AdcData|Comment|Param|Time|LIGO_LW|Array)*)>\n" \
"<!ATTLIST AdcData \n" \
"         Name CDATA #IMPLIED>\n" \
"\n" \
"<!ELEMENT AdcInterval ((AdcData|Comment|Time)*)>\n" \
"<!ATTLIST AdcInterval \n" \
"         Name CDATA #IMPLIED\n" \
"         StartTime CDATA #IMPLIED\n" \
"         DeltaT CDATA #IMPLIED>\n" \
"\n" \
"<!ELEMENT Time (#PCDATA)>\n" \
"<!ATTLIST Time \n" \
"         Name CDATA #IMPLIED\n" \
"         Type (GPS|Unix|ISO-8601) \"ISO-8601\">\n" \
"]>\n" \
"\n" \
"<LIGO_LW>\n" \
"   <Comment>metadata</Comment>\n"


#define LIGOLW_XML_FOOTER \
"</LIGO_LW>"

#define LIGOLW_XML_TABLE_FOOTER \
"\n" \
"      </Stream>\n" \
"   </Table>\n"

#define LIGOLW_XML_PROCESS \
"   <Table Name=\"processgroup:process:table\">\n" \
"      <Column Name=\"processgroup:process:program\" Type=\"lstring\"/>\n" \
"      <Column Name=\"processgroup:process:version\" Type=\"lstring\"/>\n" \
"      <Column Name=\"processgroup:process:cvs_repository\" Type=\"lstring\"/>\n" \
"      <Column Name=\"processgroup:process:cvs_entry_time\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"processgroup:process:comment\" Type=\"lstring\"/>\n" \
"      <Column Name=\"processgroup:process:is_online\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"processgroup:process:node\" Type=\"lstring\"/>\n" \
"      <Column Name=\"processgroup:process:username\" Type=\"lstring\"/>\n" \
"      <Column Name=\"processgroup:process:unix_procid\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"processgroup:process:start_time\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"processgroup:process:end_time\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"processgroup:process:jobid\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"processgroup:process:domain\" Type=\"lstring\"/>\n" \
"      <Column Name=\"processgroup:process:ifos\" Type=\"ilwd:char\"/>\n" \
"      <Stream Name=\"processgroup:process:table\" Type=\"Local\" Delimiter=\",\">\n"

#define PROCESS_ROW \
"         \"%s\",\"%s\",\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",%d,%d,%d,%d,\"%s\",\"%s\",\"process:process_id:0\""

#define LIGOLW_XML_PROCESS_PARAMS \
"   <Table Name=\"process_paramsgroup:process_params:table\">\n" \
"      <Column Name=\"process_paramsgroup:process_params:program\" Type=\"lstring\"/>\n" \
"      <Column Name=\"process_paramsgroup:process_params:process_id\" Type=\"ilwd:char\"/>\n" \
"      <Column Name=\"process_paramsgroup:process_params:param\" Type=\"lstring\"/>\n" \
"      <Column Name=\"process_paramsgroup:process_params:type\" Type=\"lstring\"/>\n" \
"      <Column Name=\"process_paramsgroup:process_params:value\" Type=\"lstring\"/>\n" \
"      <Stream Name=\"process_paramsgroup:process_params:table\" Type=\"Local\" Delimiter=\",\">\n"

#define PROCESS_PARAMS_ROW \
"         \"%s\",\"process:process_id:0\",\"%s\",\"%s\",\"%s\""

#define LIGOLW_XML_SEARCH_SUMMARY \
"   <Table Name=\"search_summarygroup:search_summary:table\">\n" \
"      <Column Name=\"search_summarygroup:search_summary:process_id\" Type=\"ilwd:char\"/>\n" \
"      <Column Name=\"search_summarygroup:search_summary:shared_object\" Type=\"lstring\"/>\n" \
"      <Column Name=\"search_summarygroup:search_summary:lalwrapper_cvs_tag\" Type=\"lstring\"/>\n" \
"      <Column Name=\"search_summarygroup:search_summary:lal_cvs_tag\" Type=\"lstring\"/>\n" \
"      <Column Name=\"search_summarygroup:search_summary:comment\" Type=\"lstring\"/>\n" \
"      <Column Name=\"search_summarygroup:search_summary:in_start_time\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"search_summarygroup:search_summary:in_start_time_ns\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"search_summarygroup:search_summary:in_end_time\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"search_summarygroup:search_summary:in_end_time_ns\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"search_summarygroup:search_summary:out_start_time\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"search_summarygroup:search_summary:out_start_time_ns\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"search_summarygroup:search_summary:out_end_time\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"search_summarygroup:search_summary:out_end_time_ns\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"search_summarygroup:search_summary:nevents\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"search_summarygroup:search_summary:nnodes\" Type=\"int_4s\"/>\n" \
"      <Stream Name=\"search_summarygroup:search_summary:table\" Type=\"Local\" Delimiter=\",\">\n"

#define SEARCH_SUMMARY_ROW \
"         \"process:process_id:0\",\"shared_object\",\"lal_cvs_tag\",\"lalwrapper_cvs_tag\",\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%d,%d"

#define LIGOLW_XML_SEARCH_SUMMVARS \
"   <Table Name=\"search_summvarsgroup:search_summvars:table\">\n" \
"      <Column Name=\"search_summvarsgroup:search_summvars:process_id\" Type=\"ilwd:char\"/>\n" \
"      <Column Name=\"search_summvarsgroup:search_summvars:name\" Type=\"lstring\"/>\n" \
"      <Column Name=\"search_summvarsgroup:search_summvars:string\" Type=\"lstring\"/>\n" \
"      <Column Name=\"search_summvarsgroup:search_summvars:value\" Type=\"real_8\"/>\n" \
"      <Stream Name=\"search_summvarsgroup:search_summvars:table\" Type=\"Local\" Delimiter=\",\">\n"

#define SEARCH_SUMMVARS_ROW \
"         \"process:process_id:0\",\"%s\",\"%s\",%22.16e"

#define LIGOLW_XML_SNGL_BURST \
"   <Table Name=\"sngl_burstgroup:sngl_burst:table\">\n" \
"      <Column Name=\"sngl_burstgroup:sngl_burst:process_id\" Type=\"ilwd:char\"/>\n" \
"      <Column Name=\"sngl_burstgroup:sngl_burst:ifo\" Type=\"lstring\"/>\n" \
"      <Column Name=\"sngl_burstgroup:sngl_burst:search\" Type=\"lstring\"/>\n" \
"      <Column Name=\"sngl_burstgroup:sngl_burst:channel\" Type=\"lstring\"/>\n" \
"      <Column Name=\"sngl_burstgroup:sngl_burst:start_time\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"sngl_burstgroup:sngl_burst:start_time_ns\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"sngl_burstgroup:sngl_burst:duration\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_burstgroup:sngl_burst:central_freq\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_burstgroup:sngl_burst:bandwidth\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_burstgroup:sngl_burst:amplitude\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_burstgroup:sngl_burst:snr\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_burstgroup:sngl_burst:confidence\" Type=\"real_4\"/>\n" \
"      <Stream Name=\"sngl_burstgroup:sngl_burst:table\" Type=\"Local\" Delimiter=\",\">\n"

#define SNGL_BURST_ROW \
"         \"process:process_id:0\",\"%s\",\"%s\",\"%s\",%d,%d,%e,%e,%e,%e,%e,%e"

#define LIGOLW_XML_SNGL_INSPIRAL \
"   <Table Name=\"sngl_inspiralgroup:sngl_inspiral:table\">\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:process_id\" Type=\"ilwd:char\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:ifo\" Type=\"lstring\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:search\" Type=\"lstring\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:channel\" Type=\"lstring\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:end_time\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:end_time_ns\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:impulse_time\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:impulse_time_ns\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:template_duration\" Type=\"real_8\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:event_duration\" Type=\"real_8\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:amplitude\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:eff_distance\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:coa_phase\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:mass1\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:mass2\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:mchirp\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:eta\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:tau0\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:tau2\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:tau3\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:tau4\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:tau5\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:ttotal\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:snr\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:chisq\" Type=\"real_4\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:chisq_dof\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"sngl_inspiralgroup:sngl_inspiral:sigmasq\" Type=\"real_8\"/>\n" \
"      <Stream Name=\"sngl_inspiralgroup:sngl_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n"

#define SNGL_INSPIRAL_ROW \
"         \"process:process_id:0\",\"%s\",\"%s\",\"%s\",%d,%d,%d,%d,%22.16e,%22.16e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%22.16e"

#define LIGOLW_XML_SUMM_VALUE \
"   <Table Name=\"summ_valuegroup:summ_value:table\">\n" \
"      <Column Name=\"summ_valuegroup:summ_value:program\" Type=\"lstring\"/>\n" \
"      <Column Name=\"summ_valuegroup:summ_value:process_id\" Type=\"ilwd:char\"/>\n" \
"      <Column Name=\"summ_valuegroup:summ_value:start_time\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"summ_valuegroup:summ_value:start_time_ns\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"summ_valuegroup:summ_value:end_time\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"summ_valuegroup:summ_value:end_time_ns\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"summ_valuegroup:summ_value:ifo\" Type=\"lstring\"/>\n" \
"      <Column Name=\"summ_valuegroup:summ_value:name\" Type=\"lstring\"/>\n" \
"      <Column Name=\"summ_valuegroup:summ_value:value\" Type=\"real_4\"/>\n" \
"      <Column Name=\"summ_valuegroup:summ_value:comment\" Type=\"lstring\"/>\n" \
"      <Stream Name=\"summ_valuegroup:summ_value:table\" Type=\"Local\" Delimiter=\",\">\n"

#define SUMM_VALUE_ROW \
"         \"%s\",\"process:process_id:0\",%d,%d,%d,%d,\"%s\",\"%s\",%e,\"%s\""

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _LIGOLWXMLHEADERS_H */
