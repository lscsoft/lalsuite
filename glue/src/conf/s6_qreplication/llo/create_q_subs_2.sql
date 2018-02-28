--Beginning of script 2--   DatabaseDB2LUOW (SEG6_LLO) [WARNING***Please do not alter this line]--
-- CONNECT TO SEG6_LLO USER XXXX using XXXX#
CREATE BUFFERPOOL BP8K2 SIZE 125 PAGESIZE 8192#
CREATE TABLESPACE QPASN2 PAGESIZE 8192 MANAGED BY SYSTEM USING 
('QPASN2_TSC') BUFFERPOOL BP8K2#
CREATE TABLE ASN.IBMQREP_DELTOMB
(
 TARGET_OWNER VARCHAR(30) NOT NULL,
 TARGET_NAME VARCHAR(128) NOT NULL,
 VERSION_TIME TIMESTAMP NOT NULL,
 VERSION_NODE SMALLINT NOT NULL,
 KEY_HASH INTEGER NOT NULL,
 PACKED_KEY VARCHAR(4096) FOR BIT DATA NOT NULL
)
 IN QPASN2#
ALTER TABLE ASN.IBMQREP_DELTOMB
 VOLATILE CARDINALITY#
CREATE INDEX ASN.IX1DELTOMB ON ASN.IBMQREP_DELTOMB
(
 VERSION_TIME ASC,
 TARGET_NAME ASC,
 TARGET_OWNER ASC,
 KEY_HASH ASC
)#
ALTER TABLE LDBD.PROCESS
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
CREATE TRIGGER LDBD.APROCESSESS1 NO CASCADE BEFORE INSERT ON
 LDBD.PROCESS
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 1*4; 
END#
CREATE TRIGGER LDBD.BPROCESSESS1 NO CASCADE BEFORE UPDATE ON
 LDBD.PROCESS
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 1*4; 
END#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('PROCESS0001', 'LDBD', 'PROCESS', 'ASN.QM2_TO_QM1.DATAQ', 'P', 'N',
 'N', 'N', 'I', 'I', '000001', 1, 2, 'NNNN', 'N', 'SEG6_LHO',
 'SEG6_LHO', 'LDBD', 'PROCESS', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0001', 'CREATOR_DB', 2, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0001', 'PROGRAM', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0001', 'VERSION', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0001', 'CVS_REPOSITORY', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0001', 'CVS_ENTRY_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0001', 'COMMENT', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0001', 'IS_ONLINE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0001', 'NODE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0001', 'USERNAME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0001', 'UNIX_PROCID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0001', 'START_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0001', 'END_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0001', 'JOBID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0001', 'DOMAIN', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0001', 'PROCESS_ID', 1, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0001', 'PARAM_SET', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0001', 'IFOS', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0001', 'INSERTION_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0001', 'ibmqrepVERTIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0001', 'ibmqrepVERNODE', 0, 'NNNNNNNNNN')#
ALTER TABLE LDBD.PROCESS_PARAMS
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
CREATE TRIGGER LDBD.APROCESS_PARAMAMS1 NO CASCADE BEFORE INSERT ON
 LDBD.PROCESS_PARAMS
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 1*4; 
END#
CREATE TRIGGER LDBD.BPROCESS_PARAMAMS1 NO CASCADE BEFORE UPDATE ON
 LDBD.PROCESS_PARAMS
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 1*4; 
END#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('PROCESS_PARAMS0001', 'LDBD', 'PROCESS_PARAMS',
 'ASN.QM2_TO_QM1.DATAQ', 'P', 'N', 'N', 'N', 'I', 'I', '000001', 1, 2,
 'NNNN', 'N', 'SEG6_LHO', 'SEG6_LHO', 'LDBD', 'PROCESS_PARAMS', 1,
 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0001', 'CREATOR_DB', 2, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0001', 'PROGRAM', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0001', 'PROCESS_ID', 1, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0001', 'PARAM', 3, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0001', 'TYPE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0001', 'VALUE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0001', 'INSERTION_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0001', 'ibmqrepVERTIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0001', 'ibmqrepVERNODE', 0, 'NNNNNNNNNN')#
ALTER TABLE LDBD.SEGMENT
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
CREATE TRIGGER LDBD.ASEGMENTENT1 NO CASCADE BEFORE INSERT ON
 LDBD.SEGMENT
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 1*4; 
END#
CREATE TRIGGER LDBD.BSEGMENTENT1 NO CASCADE BEFORE UPDATE ON
 LDBD.SEGMENT
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 1*4; 
END#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SEGMENT0001', 'LDBD', 'SEGMENT', 'ASN.QM2_TO_QM1.DATAQ', 'P', 'N',
 'N', 'N', 'I', 'I', '000001', 1, 2, 'NNNN', 'N', 'SEG6_LHO',
 'SEG6_LHO', 'LDBD', 'SEGMENT', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0001', 'CREATOR_DB', 2, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0001', 'PROCESS_ID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0001', 'SEGMENT_ID', 1, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0001', 'START_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0001', 'END_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0001', 'INSERTION_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0001', 'SEGMENT_DEF_ID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0001', 'SEGMENT_DEF_CDB', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0001', 'ibmqrepVERTIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0001', 'ibmqrepVERNODE', 0, 'NNNNNNNNNN')#
ALTER TABLE LDBD.SEGMENT_DEFINER
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
CREATE TRIGGER LDBD.ASEGMENT_DEFINNER1 NO CASCADE BEFORE INSERT ON
 LDBD.SEGMENT_DEFINER
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 1*4; 
END#
CREATE TRIGGER LDBD.BSEGMENT_DEFINNER1 NO CASCADE BEFORE UPDATE ON
 LDBD.SEGMENT_DEFINER
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 1*4; 
END#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SEGMENT_DEFINER0001', 'LDBD', 'SEGMENT_DEFINER',
 'ASN.QM2_TO_QM1.DATAQ', 'P', 'N', 'N', 'N', 'I', 'I', '000001', 1, 2,
 'NNNN', 'N', 'SEG6_LHO', 'SEG6_LHO', 'LDBD', 'SEGMENT_DEFINER', 1,
 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0001', 'CREATOR_DB', 2, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0001', 'PROCESS_ID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0001', 'SEGMENT_DEF_ID', 1, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0001', 'IFOS', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0001', 'NAME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0001', 'VERSION', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0001', 'COMMENT', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0001', 'INSERTION_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0001', 'ibmqrepVERTIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0001', 'ibmqrepVERNODE', 0, 'NNNNNNNNNN')#
ALTER TABLE LDBD.SEGMENT_SUMMARY
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
CREATE TRIGGER LDBD.ASEGMENT_SUMMAARY1 NO CASCADE BEFORE INSERT ON
 LDBD.SEGMENT_SUMMARY
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 1*4; 
END#
CREATE TRIGGER LDBD.BSEGMENT_SUMMAARY1 NO CASCADE BEFORE UPDATE ON
 LDBD.SEGMENT_SUMMARY
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 1*4; 
END#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SEGMENT_SUMMARY0001', 'LDBD', 'SEGMENT_SUMMARY',
 'ASN.QM2_TO_QM1.DATAQ', 'P', 'N', 'N', 'N', 'I', 'I', '000001', 1, 2,
 'NNNN', 'N', 'SEG6_LHO', 'SEG6_LHO', 'LDBD', 'SEGMENT_SUMMARY', 1,
 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0001', 'SEGMENT_SUM_ID', 1, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0001', 'CREATOR_DB', 2, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0001', 'START_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0001', 'END_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0001', 'COMMENT', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0001', 'SEGMENT_DEF_ID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0001', 'SEGMENT_DEF_CDB', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0001', 'PROCESS_ID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0001', 'ibmqrepVERTIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0001', 'ibmqrepVERNODE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, modelq, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase, source_type)
 VALUES
 ('PROCESS0002', 'ASN.QM1_TO_QM2.DATAQ', 'LDBD', 'PROCESS', 'LDBD',
 'PROCESS', 'IBMQREP.SPILL.MODELQ', 'SEG6_LHO', 'SEG6_LHO', 1, 'I',
 'P', 'V', 'F', 'S', '000001', 2, 1, 0, 'I', 'D')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0002', 'ASN.QM1_TO_QM2.DATAQ', 'CREATOR_DB', 'CREATOR_DB',
 'Y', 0, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0002', 'ASN.QM1_TO_QM2.DATAQ', 'PROGRAM', 'PROGRAM', 'N', 1,
 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0002', 'ASN.QM1_TO_QM2.DATAQ', 'VERSION', 'VERSION', 'N', 2,
 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0002', 'ASN.QM1_TO_QM2.DATAQ', 'CVS_REPOSITORY',
 'CVS_REPOSITORY', 'N', 3, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0002', 'ASN.QM1_TO_QM2.DATAQ', 'CVS_ENTRY_TIME',
 'CVS_ENTRY_TIME', 'N', 4, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0002', 'ASN.QM1_TO_QM2.DATAQ', 'COMMENT', 'COMMENT', 'N', 5,
 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0002', 'ASN.QM1_TO_QM2.DATAQ', 'IS_ONLINE', 'IS_ONLINE', 'N'
, 6, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0002', 'ASN.QM1_TO_QM2.DATAQ', 'NODE', 'NODE', 'N', 7, 'R',
 null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0002', 'ASN.QM1_TO_QM2.DATAQ', 'USERNAME', 'USERNAME', 'N',
 8, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0002', 'ASN.QM1_TO_QM2.DATAQ', 'UNIX_PROCID', 'UNIX_PROCID',
 'N', 9, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0002', 'ASN.QM1_TO_QM2.DATAQ', 'START_TIME', 'START_TIME',
 'N', 10, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0002', 'ASN.QM1_TO_QM2.DATAQ', 'END_TIME', 'END_TIME', 'N',
 11, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0002', 'ASN.QM1_TO_QM2.DATAQ', 'JOBID', 'JOBID', 'N', 12,
 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0002', 'ASN.QM1_TO_QM2.DATAQ', 'DOMAIN', 'DOMAIN', 'N', 13,
 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0002', 'ASN.QM1_TO_QM2.DATAQ', 'PROCESS_ID', 'PROCESS_ID',
 'Y', 14, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0002', 'ASN.QM1_TO_QM2.DATAQ', 'PARAM_SET', 'PARAM_SET', 'N'
, 15, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0002', 'ASN.QM1_TO_QM2.DATAQ', 'IFOS', 'IFOS', 'N', 16, 'R',
 null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0002', 'ASN.QM1_TO_QM2.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 17, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0002', 'ASN.QM1_TO_QM2.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0002', 'ASN.QM1_TO_QM2.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, modelq, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase, source_type)
 VALUES
 ('PROCESS_PARAMS0002', 'ASN.QM1_TO_QM2.DATAQ', 'LDBD',
 'PROCESS_PARAMS', 'LDBD', 'PROCESS_PARAMS', 'IBMQREP.SPILL.MODELQ',
 'SEG6_LHO', 'SEG6_LHO', 1, 'I', 'P', 'V', 'F', 'S', '000001', 2, 1, 0
, 'I', 'D')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0002', 'ASN.QM1_TO_QM2.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0002', 'ASN.QM1_TO_QM2.DATAQ', 'PROGRAM', 'PROGRAM',
 'N', 1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0002', 'ASN.QM1_TO_QM2.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'Y', 2, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0002', 'ASN.QM1_TO_QM2.DATAQ', 'PARAM', 'PARAM', 'Y',
 3, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0002', 'ASN.QM1_TO_QM2.DATAQ', 'TYPE', 'TYPE', 'N', 4
, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0002', 'ASN.QM1_TO_QM2.DATAQ', 'VALUE', 'VALUE', 'N',
 5, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0002', 'ASN.QM1_TO_QM2.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 6, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0002', 'ASN.QM1_TO_QM2.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0002', 'ASN.QM1_TO_QM2.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, modelq, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase, source_type)
 VALUES
 ('SEGMENT0002', 'ASN.QM1_TO_QM2.DATAQ', 'LDBD', 'SEGMENT', 'LDBD',
 'SEGMENT', 'IBMQREP.SPILL.MODELQ', 'SEG6_LHO', 'SEG6_LHO', 1, 'I',
 'P', 'V', 'F', 'S', '000001', 2, 1, 0, 'I', 'D')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0002', 'ASN.QM1_TO_QM2.DATAQ', 'CREATOR_DB', 'CREATOR_DB',
 'Y', 0, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0002', 'ASN.QM1_TO_QM2.DATAQ', 'PROCESS_ID', 'PROCESS_ID',
 'N', 1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0002', 'ASN.QM1_TO_QM2.DATAQ', 'SEGMENT_ID', 'SEGMENT_ID',
 'Y', 2, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0002', 'ASN.QM1_TO_QM2.DATAQ', 'START_TIME', 'START_TIME',
 'N', 3, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0002', 'ASN.QM1_TO_QM2.DATAQ', 'END_TIME', 'END_TIME', 'N',
 4, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0002', 'ASN.QM1_TO_QM2.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 5, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0002', 'ASN.QM1_TO_QM2.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'N', 6, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0002', 'ASN.QM1_TO_QM2.DATAQ', 'SEGMENT_DEF_CDB',
 'SEGMENT_DEF_CDB', 'N', 7, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0002', 'ASN.QM1_TO_QM2.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0002', 'ASN.QM1_TO_QM2.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, modelq, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase, source_type)
 VALUES
 ('SEGMENT_DEFINER0002', 'ASN.QM1_TO_QM2.DATAQ', 'LDBD',
 'SEGMENT_DEFINER', 'LDBD', 'SEGMENT_DEFINER', 'IBMQREP.SPILL.MODELQ',
 'SEG6_LHO', 'SEG6_LHO', 1, 'I', 'P', 'V', 'F', 'S', '000001', 2, 1, 0
, 'I', 'D')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0002', 'ASN.QM1_TO_QM2.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0002', 'ASN.QM1_TO_QM2.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'N', 1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0002', 'ASN.QM1_TO_QM2.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'Y', 2, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0002', 'ASN.QM1_TO_QM2.DATAQ', 'IFOS', 'IFOS', 'N',
 3, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0002', 'ASN.QM1_TO_QM2.DATAQ', 'NAME', 'NAME', 'N',
 4, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0002', 'ASN.QM1_TO_QM2.DATAQ', 'VERSION', 'VERSION',
 'N', 5, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0002', 'ASN.QM1_TO_QM2.DATAQ', 'COMMENT', 'COMMENT',
 'N', 6, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0002', 'ASN.QM1_TO_QM2.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 7, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0002', 'ASN.QM1_TO_QM2.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0002', 'ASN.QM1_TO_QM2.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, modelq, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase, source_type)
 VALUES
 ('SEGMENT_SUMMARY0002', 'ASN.QM1_TO_QM2.DATAQ', 'LDBD',
 'SEGMENT_SUMMARY', 'LDBD', 'SEGMENT_SUMMARY', 'IBMQREP.SPILL.MODELQ',
 'SEG6_LHO', 'SEG6_LHO', 1, 'I', 'P', 'V', 'F', 'S', '000001', 2, 1, 0
, 'I', 'D')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0002', 'ASN.QM1_TO_QM2.DATAQ', 'SEGMENT_SUM_ID',
 'SEGMENT_SUM_ID', 'Y', 0, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0002', 'ASN.QM1_TO_QM2.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0002', 'ASN.QM1_TO_QM2.DATAQ', 'START_TIME',
 'START_TIME', 'N', 2, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0002', 'ASN.QM1_TO_QM2.DATAQ', 'END_TIME',
 'END_TIME', 'N', 3, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0002', 'ASN.QM1_TO_QM2.DATAQ', 'COMMENT', 'COMMENT',
 'N', 4, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0002', 'ASN.QM1_TO_QM2.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'N', 5, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0002', 'ASN.QM1_TO_QM2.DATAQ', 'SEGMENT_DEF_CDB',
 'SEGMENT_DEF_CDB', 'N', 6, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0002', 'ASN.QM1_TO_QM2.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'N', 7, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0002', 'ASN.QM1_TO_QM2.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0002', 'ASN.QM1_TO_QM2.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, modelq, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase, source_type)
 VALUES
 ('PROCESS0005', 'ASN.QM3_TO_QM2.DATAQ', 'LDBD', 'PROCESS', 'LDBD',
 'PROCESS', 'IBMQREP.SPILL.MODELQ', 'SEG6_CIT', 'SEG6_CIT', 1, 'I',
 'P', 'V', 'F', 'S', '000001', 3, 1, 0, 'I', 'D')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0005', 'ASN.QM3_TO_QM2.DATAQ', 'CREATOR_DB', 'CREATOR_DB',
 'Y', 0, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0005', 'ASN.QM3_TO_QM2.DATAQ', 'PROGRAM', 'PROGRAM', 'N', 1,
 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0005', 'ASN.QM3_TO_QM2.DATAQ', 'VERSION', 'VERSION', 'N', 2,
 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0005', 'ASN.QM3_TO_QM2.DATAQ', 'CVS_REPOSITORY',
 'CVS_REPOSITORY', 'N', 3, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0005', 'ASN.QM3_TO_QM2.DATAQ', 'CVS_ENTRY_TIME',
 'CVS_ENTRY_TIME', 'N', 4, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0005', 'ASN.QM3_TO_QM2.DATAQ', 'COMMENT', 'COMMENT', 'N', 5,
 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0005', 'ASN.QM3_TO_QM2.DATAQ', 'IS_ONLINE', 'IS_ONLINE', 'N'
, 6, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0005', 'ASN.QM3_TO_QM2.DATAQ', 'NODE', 'NODE', 'N', 7, 'R',
 null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0005', 'ASN.QM3_TO_QM2.DATAQ', 'USERNAME', 'USERNAME', 'N',
 8, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0005', 'ASN.QM3_TO_QM2.DATAQ', 'UNIX_PROCID', 'UNIX_PROCID',
 'N', 9, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0005', 'ASN.QM3_TO_QM2.DATAQ', 'START_TIME', 'START_TIME',
 'N', 10, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0005', 'ASN.QM3_TO_QM2.DATAQ', 'END_TIME', 'END_TIME', 'N',
 11, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0005', 'ASN.QM3_TO_QM2.DATAQ', 'JOBID', 'JOBID', 'N', 12,
 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0005', 'ASN.QM3_TO_QM2.DATAQ', 'DOMAIN', 'DOMAIN', 'N', 13,
 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0005', 'ASN.QM3_TO_QM2.DATAQ', 'PROCESS_ID', 'PROCESS_ID',
 'Y', 14, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0005', 'ASN.QM3_TO_QM2.DATAQ', 'PARAM_SET', 'PARAM_SET', 'N'
, 15, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0005', 'ASN.QM3_TO_QM2.DATAQ', 'IFOS', 'IFOS', 'N', 16, 'R',
 null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0005', 'ASN.QM3_TO_QM2.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 17, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0005', 'ASN.QM3_TO_QM2.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0005', 'ASN.QM3_TO_QM2.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, modelq, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase, source_type)
 VALUES
 ('PROCESS_PARAMS0005', 'ASN.QM3_TO_QM2.DATAQ', 'LDBD',
 'PROCESS_PARAMS', 'LDBD', 'PROCESS_PARAMS', 'IBMQREP.SPILL.MODELQ',
 'SEG6_CIT', 'SEG6_CIT', 1, 'I', 'P', 'V', 'F', 'S', '000001', 3, 1, 0
, 'I', 'D')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0005', 'ASN.QM3_TO_QM2.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0005', 'ASN.QM3_TO_QM2.DATAQ', 'PROGRAM', 'PROGRAM',
 'N', 1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0005', 'ASN.QM3_TO_QM2.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'Y', 2, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0005', 'ASN.QM3_TO_QM2.DATAQ', 'PARAM', 'PARAM', 'Y',
 3, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0005', 'ASN.QM3_TO_QM2.DATAQ', 'TYPE', 'TYPE', 'N', 4
, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0005', 'ASN.QM3_TO_QM2.DATAQ', 'VALUE', 'VALUE', 'N',
 5, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0005', 'ASN.QM3_TO_QM2.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 6, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0005', 'ASN.QM3_TO_QM2.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0005', 'ASN.QM3_TO_QM2.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, modelq, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase, source_type)
 VALUES
 ('SEGMENT0005', 'ASN.QM3_TO_QM2.DATAQ', 'LDBD', 'SEGMENT', 'LDBD',
 'SEGMENT', 'IBMQREP.SPILL.MODELQ', 'SEG6_CIT', 'SEG6_CIT', 1, 'I',
 'P', 'V', 'F', 'S', '000001', 3, 1, 0, 'I', 'D')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0005', 'ASN.QM3_TO_QM2.DATAQ', 'CREATOR_DB', 'CREATOR_DB',
 'Y', 0, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0005', 'ASN.QM3_TO_QM2.DATAQ', 'PROCESS_ID', 'PROCESS_ID',
 'N', 1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0005', 'ASN.QM3_TO_QM2.DATAQ', 'SEGMENT_ID', 'SEGMENT_ID',
 'Y', 2, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0005', 'ASN.QM3_TO_QM2.DATAQ', 'START_TIME', 'START_TIME',
 'N', 3, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0005', 'ASN.QM3_TO_QM2.DATAQ', 'END_TIME', 'END_TIME', 'N',
 4, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0005', 'ASN.QM3_TO_QM2.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 5, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0005', 'ASN.QM3_TO_QM2.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'N', 6, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0005', 'ASN.QM3_TO_QM2.DATAQ', 'SEGMENT_DEF_CDB',
 'SEGMENT_DEF_CDB', 'N', 7, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0005', 'ASN.QM3_TO_QM2.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0005', 'ASN.QM3_TO_QM2.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, modelq, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase, source_type)
 VALUES
 ('SEGMENT_DEFINER0005', 'ASN.QM3_TO_QM2.DATAQ', 'LDBD',
 'SEGMENT_DEFINER', 'LDBD', 'SEGMENT_DEFINER', 'IBMQREP.SPILL.MODELQ',
 'SEG6_CIT', 'SEG6_CIT', 1, 'I', 'P', 'V', 'F', 'S', '000001', 3, 1, 0
, 'I', 'D')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0005', 'ASN.QM3_TO_QM2.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0005', 'ASN.QM3_TO_QM2.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'N', 1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0005', 'ASN.QM3_TO_QM2.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'Y', 2, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0005', 'ASN.QM3_TO_QM2.DATAQ', 'IFOS', 'IFOS', 'N',
 3, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0005', 'ASN.QM3_TO_QM2.DATAQ', 'NAME', 'NAME', 'N',
 4, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0005', 'ASN.QM3_TO_QM2.DATAQ', 'VERSION', 'VERSION',
 'N', 5, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0005', 'ASN.QM3_TO_QM2.DATAQ', 'COMMENT', 'COMMENT',
 'N', 6, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0005', 'ASN.QM3_TO_QM2.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 7, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0005', 'ASN.QM3_TO_QM2.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0005', 'ASN.QM3_TO_QM2.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, modelq, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase, source_type)
 VALUES
 ('SEGMENT_SUMMARY0005', 'ASN.QM3_TO_QM2.DATAQ', 'LDBD',
 'SEGMENT_SUMMARY', 'LDBD', 'SEGMENT_SUMMARY', 'IBMQREP.SPILL.MODELQ',
 'SEG6_CIT', 'SEG6_CIT', 1, 'I', 'P', 'V', 'F', 'S', '000001', 3, 1, 0
, 'I', 'D')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0005', 'ASN.QM3_TO_QM2.DATAQ', 'SEGMENT_SUM_ID',
 'SEGMENT_SUM_ID', 'Y', 0, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0005', 'ASN.QM3_TO_QM2.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0005', 'ASN.QM3_TO_QM2.DATAQ', 'START_TIME',
 'START_TIME', 'N', 2, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0005', 'ASN.QM3_TO_QM2.DATAQ', 'END_TIME',
 'END_TIME', 'N', 3, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0005', 'ASN.QM3_TO_QM2.DATAQ', 'COMMENT', 'COMMENT',
 'N', 4, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0005', 'ASN.QM3_TO_QM2.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'N', 5, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0005', 'ASN.QM3_TO_QM2.DATAQ', 'SEGMENT_DEF_CDB',
 'SEGMENT_DEF_CDB', 'N', 6, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0005', 'ASN.QM3_TO_QM2.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'N', 7, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0005', 'ASN.QM3_TO_QM2.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0005', 'ASN.QM3_TO_QM2.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('PROCESS0006', 'LDBD', 'PROCESS', 'ASN.QM2_TO_QM3.DATAQ', 'P', 'N',
 'N', 'N', 'I', 'I', '000001', 1, 3, 'NNNN', 'N', 'SEG6_CIT',
 'SEG6_CIT', 'LDBD', 'PROCESS', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0006', 'CREATOR_DB', 2, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0006', 'PROGRAM', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0006', 'VERSION', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0006', 'CVS_REPOSITORY', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0006', 'CVS_ENTRY_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0006', 'COMMENT', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0006', 'IS_ONLINE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0006', 'NODE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0006', 'USERNAME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0006', 'UNIX_PROCID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0006', 'START_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0006', 'END_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0006', 'JOBID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0006', 'DOMAIN', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0006', 'PROCESS_ID', 1, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0006', 'PARAM_SET', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0006', 'IFOS', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0006', 'INSERTION_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0006', 'ibmqrepVERTIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0006', 'ibmqrepVERNODE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('PROCESS_PARAMS0006', 'LDBD', 'PROCESS_PARAMS',
 'ASN.QM2_TO_QM3.DATAQ', 'P', 'N', 'N', 'N', 'I', 'I', '000001', 1, 3,
 'NNNN', 'N', 'SEG6_CIT', 'SEG6_CIT', 'LDBD', 'PROCESS_PARAMS', 1,
 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0006', 'CREATOR_DB', 2, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0006', 'PROGRAM', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0006', 'PROCESS_ID', 1, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0006', 'PARAM', 3, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0006', 'TYPE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0006', 'VALUE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0006', 'INSERTION_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0006', 'ibmqrepVERTIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0006', 'ibmqrepVERNODE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SEGMENT0006', 'LDBD', 'SEGMENT', 'ASN.QM2_TO_QM3.DATAQ', 'P', 'N',
 'N', 'N', 'I', 'I', '000001', 1, 3, 'NNNN', 'N', 'SEG6_CIT',
 'SEG6_CIT', 'LDBD', 'SEGMENT', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0006', 'CREATOR_DB', 2, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0006', 'PROCESS_ID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0006', 'SEGMENT_ID', 1, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0006', 'START_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0006', 'END_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0006', 'INSERTION_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0006', 'SEGMENT_DEF_ID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0006', 'SEGMENT_DEF_CDB', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0006', 'ibmqrepVERTIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0006', 'ibmqrepVERNODE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SEGMENT_DEFINER0006', 'LDBD', 'SEGMENT_DEFINER',
 'ASN.QM2_TO_QM3.DATAQ', 'P', 'N', 'N', 'N', 'I', 'I', '000001', 1, 3,
 'NNNN', 'N', 'SEG6_CIT', 'SEG6_CIT', 'LDBD', 'SEGMENT_DEFINER', 1,
 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0006', 'CREATOR_DB', 2, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0006', 'PROCESS_ID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0006', 'SEGMENT_DEF_ID', 1, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0006', 'IFOS', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0006', 'NAME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0006', 'VERSION', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0006', 'COMMENT', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0006', 'INSERTION_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0006', 'ibmqrepVERTIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0006', 'ibmqrepVERNODE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SEGMENT_SUMMARY0006', 'LDBD', 'SEGMENT_SUMMARY',
 'ASN.QM2_TO_QM3.DATAQ', 'P', 'N', 'N', 'N', 'I', 'I', '000001', 1, 3,
 'NNNN', 'N', 'SEG6_CIT', 'SEG6_CIT', 'LDBD', 'SEGMENT_SUMMARY', 1,
 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0006', 'SEGMENT_SUM_ID', 1, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0006', 'CREATOR_DB', 2, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0006', 'START_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0006', 'END_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0006', 'COMMENT', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0006', 'SEGMENT_DEF_ID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0006', 'SEGMENT_DEF_CDB', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0006', 'PROCESS_ID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0006', 'ibmqrepVERTIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0006', 'ibmqrepVERNODE', 0, 'NNNNNNNNNN')#
-- COMMIT#