/*----------------------------------------------------------------------- 
 * 
 * File Name: girdmatchxml.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#define LIGO_LW_HEADER \
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


#define LIGO_LW_FOOTER \
"</LIGO_LW>"

#define TABLE_FOOTER \
"      </Stream>\n" \
"   </Table>\n"

#define PROCESS_HEADER \
"   <Table Name=\"processgroup:process:table\">\n" \
"      <Column Name=\"processgroup:process:program\" Type=\"lstring\"/>\n" \
"      <Column Name=\"processgroup:process:version\" Type=\"lstring\"/>\n" \
"      <Column Name=\"processgroup:process:cvs_repository\" Type=\"lstring\"/>\n" \
"      <Column Name=\"processgroup:process:cvs_entry_time\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"processgroup:process:is_online\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"processgroup:process:node\" Type=\"lstring\"/>\n" \
"      <Column Name=\"processgroup:process:username\" Type=\"lstring\"/>\n" \
"      <Column Name=\"processgroup:process:unix_procid\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"processgroup:process:start_time\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"processgroup:process:end_time\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"processgroup:process:jobid\" Type=\"int_4s\"/>\n" \
"      <Column Name=\"processgroup:process:domain\" Type=\"lstring\"/>\n" \

#define PROCESS_COMMENT_HEADER \
"      <Column Name=\"processgroup:process:comment\" Type=\"lstring\"/>\n"

#define PROCESS_STREAM_HEADER \
"      <Column Name=\"processgroup:process:process_id\" Type=\"ilwd:char\"/>\n" \
"      <Stream Name=\"processgroup:process:table\" Type=\"Local\" Delimiter=\",\">\n"

#define PROCESS_PARAMS_HEADER \
"   <Table Name=\"process_paramsgroup:process_params:table\">\n" \
"      <Column Name=\"process_paramsgroup:process_params:program\" Type=\"lstring\"/>\n" \
"      <Column Name=\"process_paramsgroup:process_params:process_id\" Type=\"ilwd:char\"/>\n" \
"      <Column Name=\"process_paramsgroup:process_params:param\" Type=\"lstring\"/>\n" \
"      <Column Name=\"process_paramsgroup:process_params:type\" Type=\"lstring\"/>\n" \
"      <Column Name=\"process_paramsgroup:process_params:value\" Type=\"lstring\"/>\n" \
"      <Stream Name=\"process_paramsgroup:process_params:table\" Type=\"Local\" Delimiter=\",\">\n"

#define PPARAMS \
"         \"lalapps_inspiral\",\"process:process_id:0\","

