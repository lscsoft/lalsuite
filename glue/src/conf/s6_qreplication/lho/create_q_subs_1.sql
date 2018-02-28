--Beginning of script 1--   DatabaseDB2LUOW (SEG6_LHO) [WARNING***Please do not alter this line]--
-- CONNECT TO SEG6_LHO USER XXXX using XXXX#
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
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, modelq, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase, source_type)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD', 'PROCESS', 'LDBD',
 'PROCESS', 'IBMQREP.SPILL.MODELQ', 'SEG6_LLO', 'SEG6_LLO', 1, 'I',
 'P', 'V', 'F', 'S', '000001', 1, 2, 0, 'I', 'D')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB', 'CREATOR_DB',
 'Y', 0, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROGRAM', 'PROGRAM', 'N', 1,
 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'VERSION', 'VERSION', 'N', 2,
 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'CVS_REPOSITORY',
 'CVS_REPOSITORY', 'N', 3, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'CVS_ENTRY_TIME',
 'CVS_ENTRY_TIME', 'N', 4, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'COMMENT', 'COMMENT', 'N', 5,
 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'IS_ONLINE', 'IS_ONLINE', 'N'
, 6, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'NODE', 'NODE', 'N', 7, 'R',
 null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'USERNAME', 'USERNAME', 'N',
 8, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'UNIX_PROCID', 'UNIX_PROCID',
 'N', 9, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_TIME', 'START_TIME',
 'N', 10, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'END_TIME', 'END_TIME', 'N',
 11, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'JOBID', 'JOBID', 'N', 12,
 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'DOMAIN', 'DOMAIN', 'N', 13,
 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROCESS_ID', 'PROCESS_ID',
 'Y', 14, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'PARAM_SET', 'PARAM_SET', 'N'
, 15, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'IFOS', 'IFOS', 'N', 16, 'R',
 null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 17, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', -1, 'R', null, null)#
CREATE TRIGGER LDBD.APROCESSESS NO CASCADE BEFORE INSERT ON
 LDBD.PROCESS
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
END#
CREATE TRIGGER LDBD.BPROCESSESS NO CASCADE BEFORE UPDATE ON
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
new."ibmqrepVERNODE" = 2*4; 
END#
ALTER TABLE LDBD.PROCESS_PARAMS
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, modelq, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase, source_type)
 VALUES
 ('PROCESS_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD',
 'PROCESS_PARAMS', 'LDBD', 'PROCESS_PARAMS', 'IBMQREP.SPILL.MODELQ',
 'SEG6_LLO', 'SEG6_LLO', 1, 'I', 'P', 'V', 'F', 'S', '000001', 1, 2, 0
, 'I', 'D')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROGRAM', 'PROGRAM',
 'N', 1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'Y', 2, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'PARAM', 'PARAM', 'Y',
 3, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'TYPE', 'TYPE', 'N', 4
, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'VALUE', 'VALUE', 'N',
 5, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 6, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', -1, 'R', null, null)#
CREATE TRIGGER LDBD.APROCESS_PARAMAMS NO CASCADE BEFORE INSERT ON
 LDBD.PROCESS_PARAMS
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
END#
CREATE TRIGGER LDBD.BPROCESS_PARAMAMS NO CASCADE BEFORE UPDATE ON
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
new."ibmqrepVERNODE" = 2*4; 
END#
ALTER TABLE LDBD.SEGMENT
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, modelq, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase, source_type)
 VALUES
 ('SEGMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD', 'SEGMENT', 'LDBD',
 'SEGMENT', 'IBMQREP.SPILL.MODELQ', 'SEG6_LLO', 'SEG6_LLO', 1, 'I',
 'P', 'V', 'F', 'S', '000001', 1, 2, 0, 'I', 'D')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB', 'CREATOR_DB',
 'Y', 0, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROCESS_ID', 'PROCESS_ID',
 'N', 1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'SEGMENT_ID', 'SEGMENT_ID',
 'Y', 2, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_TIME', 'START_TIME',
 'N', 3, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'END_TIME', 'END_TIME', 'N',
 4, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 5, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'N', 6, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'SEGMENT_DEF_CDB',
 'SEGMENT_DEF_CDB', 'N', 7, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', -1, 'R', null, null)#
CREATE TRIGGER LDBD.ASEGMENTENT NO CASCADE BEFORE INSERT ON
 LDBD.SEGMENT
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
END#
CREATE TRIGGER LDBD.BSEGMENTENT NO CASCADE BEFORE UPDATE ON
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
new."ibmqrepVERNODE" = 2*4; 
END#
ALTER TABLE LDBD.SEGMENT_DEFINER
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, modelq, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase, source_type)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD',
 'SEGMENT_DEFINER', 'LDBD', 'SEGMENT_DEFINER', 'IBMQREP.SPILL.MODELQ',
 'SEG6_LLO', 'SEG6_LLO', 1, 'I', 'P', 'V', 'F', 'S', '000001', 1, 2, 0
, 'I', 'D')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'N', 1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'Y', 2, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'IFOS', 'IFOS', 'N',
 3, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'NAME', 'NAME', 'N',
 4, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'VERSION', 'VERSION',
 'N', 5, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'COMMENT', 'COMMENT',
 'N', 6, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 7, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', -1, 'R', null, null)#
CREATE TRIGGER LDBD.ASEGMENT_DEFINNER1 NO CASCADE BEFORE INSERT ON
 LDBD.SEGMENT_DEFINER
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
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
new."ibmqrepVERNODE" = 2*4; 
END#
ALTER TABLE LDBD.SEGMENT_SUMMARY
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, modelq, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase, source_type)
 VALUES
 ('SEGMENT_SUMMARY0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD',
 'SEGMENT_SUMMARY', 'LDBD', 'SEGMENT_SUMMARY', 'IBMQREP.SPILL.MODELQ',
 'SEG6_LLO', 'SEG6_LLO', 1, 'I', 'P', 'V', 'F', 'S', '000001', 1, 2, 0
, 'I', 'D')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0001', 'ASN.QM2_TO_QM1.DATAQ', 'SEGMENT_SUM_ID',
 'SEGMENT_SUM_ID', 'Y', 0, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_TIME',
 'START_TIME', 'N', 2, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0001', 'ASN.QM2_TO_QM1.DATAQ', 'END_TIME',
 'END_TIME', 'N', 3, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0001', 'ASN.QM2_TO_QM1.DATAQ', 'COMMENT', 'COMMENT',
 'N', 4, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0001', 'ASN.QM2_TO_QM1.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'N', 5, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0001', 'ASN.QM2_TO_QM1.DATAQ', 'SEGMENT_DEF_CDB',
 'SEGMENT_DEF_CDB', 'N', 6, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'N', 7, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', -1, 'R', null, null)#
CREATE TRIGGER LDBD.ASEGMENT_SUMMAARY NO CASCADE BEFORE INSERT ON
 LDBD.SEGMENT_SUMMARY
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
END#
CREATE TRIGGER LDBD.BSEGMENT_SUMMAARY NO CASCADE BEFORE UPDATE ON
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
new."ibmqrepVERNODE" = 2*4; 
END#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('PROCESS0002', 'LDBD', 'PROCESS', 'ASN.QM1_TO_QM2.DATAQ', 'P', 'N',
 'N', 'N', 'I', 'I', '000001', 2, 1, 'NNNN', 'N', 'SEG6_LLO',
 'SEG6_LLO', 'LDBD', 'PROCESS', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0002', 'CREATOR_DB', 2, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0002', 'PROGRAM', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0002', 'VERSION', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0002', 'CVS_REPOSITORY', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0002', 'CVS_ENTRY_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0002', 'COMMENT', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0002', 'IS_ONLINE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0002', 'NODE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0002', 'USERNAME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0002', 'UNIX_PROCID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0002', 'START_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0002', 'END_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0002', 'JOBID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0002', 'DOMAIN', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0002', 'PROCESS_ID', 1, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0002', 'PARAM_SET', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0002', 'IFOS', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0002', 'INSERTION_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0002', 'ibmqrepVERTIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0002', 'ibmqrepVERNODE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('PROCESS_PARAMS0002', 'LDBD', 'PROCESS_PARAMS',
 'ASN.QM1_TO_QM2.DATAQ', 'P', 'N', 'N', 'N', 'I', 'I', '000001', 2, 1,
 'NNNN', 'N', 'SEG6_LLO', 'SEG6_LLO', 'LDBD', 'PROCESS_PARAMS', 1,
 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0002', 'CREATOR_DB', 2, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0002', 'PROGRAM', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0002', 'PROCESS_ID', 1, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0002', 'PARAM', 3, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0002', 'TYPE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0002', 'VALUE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0002', 'INSERTION_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0002', 'ibmqrepVERTIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0002', 'ibmqrepVERNODE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SEGMENT0002', 'LDBD', 'SEGMENT', 'ASN.QM1_TO_QM2.DATAQ', 'P', 'N',
 'N', 'N', 'I', 'I', '000001', 2, 1, 'NNNN', 'N', 'SEG6_LLO',
 'SEG6_LLO', 'LDBD', 'SEGMENT', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0002', 'CREATOR_DB', 2, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0002', 'PROCESS_ID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0002', 'SEGMENT_ID', 1, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0002', 'START_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0002', 'END_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0002', 'INSERTION_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0002', 'SEGMENT_DEF_ID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0002', 'SEGMENT_DEF_CDB', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0002', 'ibmqrepVERTIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0002', 'ibmqrepVERNODE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SEGMENT_DEFINER0002', 'LDBD', 'SEGMENT_DEFINER',
 'ASN.QM1_TO_QM2.DATAQ', 'P', 'N', 'N', 'N', 'I', 'I', '000001', 2, 1,
 'NNNN', 'N', 'SEG6_LLO', 'SEG6_LLO', 'LDBD', 'SEGMENT_DEFINER', 1,
 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0002', 'CREATOR_DB', 2, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0002', 'PROCESS_ID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0002', 'SEGMENT_DEF_ID', 1, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0002', 'IFOS', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0002', 'NAME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0002', 'VERSION', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0002', 'COMMENT', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0002', 'INSERTION_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0002', 'ibmqrepVERTIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0002', 'ibmqrepVERNODE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SEGMENT_SUMMARY0002', 'LDBD', 'SEGMENT_SUMMARY',
 'ASN.QM1_TO_QM2.DATAQ', 'P', 'N', 'N', 'N', 'I', 'I', '000001', 2, 1,
 'NNNN', 'N', 'SEG6_LLO', 'SEG6_LLO', 'LDBD', 'SEGMENT_SUMMARY', 1,
 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0002', 'SEGMENT_SUM_ID', 1, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0002', 'CREATOR_DB', 2, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0002', 'START_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0002', 'END_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0002', 'COMMENT', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0002', 'SEGMENT_DEF_ID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0002', 'SEGMENT_DEF_CDB', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0002', 'PROCESS_ID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0002', 'ibmqrepVERTIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0002', 'ibmqrepVERNODE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, modelq, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase, source_type)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD', 'PROCESS', 'LDBD',
 'PROCESS', 'IBMQREP.SPILL.MODELQ', 'SEG6_CIT', 'SEG6_CIT', 1, 'I',
 'P', 'V', 'F', 'S', '000001', 3, 2, 0, 'I', 'D')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB', 'CREATOR_DB',
 'Y', 0, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROGRAM', 'PROGRAM', 'N', 1,
 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'VERSION', 'VERSION', 'N', 2,
 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'CVS_REPOSITORY',
 'CVS_REPOSITORY', 'N', 3, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'CVS_ENTRY_TIME',
 'CVS_ENTRY_TIME', 'N', 4, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'COMMENT', 'COMMENT', 'N', 5,
 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'IS_ONLINE', 'IS_ONLINE', 'N'
, 6, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'NODE', 'NODE', 'N', 7, 'R',
 null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'USERNAME', 'USERNAME', 'N',
 8, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'UNIX_PROCID', 'UNIX_PROCID',
 'N', 9, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_TIME', 'START_TIME',
 'N', 10, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'END_TIME', 'END_TIME', 'N',
 11, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'JOBID', 'JOBID', 'N', 12,
 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'DOMAIN', 'DOMAIN', 'N', 13,
 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROCESS_ID', 'PROCESS_ID',
 'Y', 14, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'PARAM_SET', 'PARAM_SET', 'N'
, 15, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'IFOS', 'IFOS', 'N', 16, 'R',
 null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 17, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, modelq, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase, source_type)
 VALUES
 ('PROCESS_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD',
 'PROCESS_PARAMS', 'LDBD', 'PROCESS_PARAMS', 'IBMQREP.SPILL.MODELQ',
 'SEG6_CIT', 'SEG6_CIT', 1, 'I', 'P', 'V', 'F', 'S', '000001', 3, 2, 0
, 'I', 'D')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROGRAM', 'PROGRAM',
 'N', 1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'Y', 2, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'PARAM', 'PARAM', 'Y',
 3, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'TYPE', 'TYPE', 'N', 4
, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'VALUE', 'VALUE', 'N',
 5, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 6, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('PROCESS_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, modelq, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase, source_type)
 VALUES
 ('SEGMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD', 'SEGMENT', 'LDBD',
 'SEGMENT', 'IBMQREP.SPILL.MODELQ', 'SEG6_CIT', 'SEG6_CIT', 1, 'I',
 'P', 'V', 'F', 'S', '000001', 3, 2, 0, 'I', 'D')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB', 'CREATOR_DB',
 'Y', 0, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROCESS_ID', 'PROCESS_ID',
 'N', 1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'SEGMENT_ID', 'SEGMENT_ID',
 'Y', 2, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_TIME', 'START_TIME',
 'N', 3, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'END_TIME', 'END_TIME', 'N',
 4, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 5, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'N', 6, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'SEGMENT_DEF_CDB',
 'SEGMENT_DEF_CDB', 'N', 7, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, modelq, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase, source_type)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD',
 'SEGMENT_DEFINER', 'LDBD', 'SEGMENT_DEFINER', 'IBMQREP.SPILL.MODELQ',
 'SEG6_CIT', 'SEG6_CIT', 1, 'I', 'P', 'V', 'F', 'S', '000001', 3, 2, 0
, 'I', 'D')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'N', 1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'Y', 2, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'IFOS', 'IFOS', 'N',
 3, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'NAME', 'NAME', 'N',
 4, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'VERSION', 'VERSION',
 'N', 5, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'COMMENT', 'COMMENT',
 'N', 6, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 7, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, modelq, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase, source_type)
 VALUES
 ('SEGMENT_SUMMARY0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD',
 'SEGMENT_SUMMARY', 'LDBD', 'SEGMENT_SUMMARY', 'IBMQREP.SPILL.MODELQ',
 'SEG6_CIT', 'SEG6_CIT', 1, 'I', 'P', 'V', 'F', 'S', '000001', 3, 2, 0
, 'I', 'D')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0003', 'ASN.QM3_TO_QM1.DATAQ', 'SEGMENT_SUM_ID',
 'SEGMENT_SUM_ID', 'Y', 0, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_TIME',
 'START_TIME', 'N', 2, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0003', 'ASN.QM3_TO_QM1.DATAQ', 'END_TIME',
 'END_TIME', 'N', 3, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0003', 'ASN.QM3_TO_QM1.DATAQ', 'COMMENT', 'COMMENT',
 'N', 4, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0003', 'ASN.QM3_TO_QM1.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'N', 5, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0003', 'ASN.QM3_TO_QM1.DATAQ', 'SEGMENT_DEF_CDB',
 'SEGMENT_DEF_CDB', 'N', 6, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'N', 7, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo
, mapping_type, src_col_map, bef_targ_colname)
 VALUES
 ('SEGMENT_SUMMARY0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', -1, 'R', null, null)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('PROCESS0004', 'LDBD', 'PROCESS', 'ASN.QM1_TO_QM3.DATAQ', 'P', 'N',
 'N', 'N', 'I', 'I', '000001', 2, 3, 'NNNN', 'N', 'SEG6_CIT',
 'SEG6_CIT', 'LDBD', 'PROCESS', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0004', 'CREATOR_DB', 2, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0004', 'PROGRAM', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0004', 'VERSION', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0004', 'CVS_REPOSITORY', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0004', 'CVS_ENTRY_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0004', 'COMMENT', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0004', 'IS_ONLINE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0004', 'NODE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0004', 'USERNAME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0004', 'UNIX_PROCID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0004', 'START_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0004', 'END_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0004', 'JOBID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0004', 'DOMAIN', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0004', 'PROCESS_ID', 1, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0004', 'PARAM_SET', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0004', 'IFOS', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0004', 'INSERTION_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0004', 'ibmqrepVERTIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS0004', 'ibmqrepVERNODE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('PROCESS_PARAMS0004', 'LDBD', 'PROCESS_PARAMS',
 'ASN.QM1_TO_QM3.DATAQ', 'P', 'N', 'N', 'N', 'I', 'I', '000001', 2, 3,
 'NNNN', 'N', 'SEG6_CIT', 'SEG6_CIT', 'LDBD', 'PROCESS_PARAMS', 1,
 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0004', 'CREATOR_DB', 2, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0004', 'PROGRAM', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0004', 'PROCESS_ID', 1, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0004', 'PARAM', 3, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0004', 'TYPE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0004', 'VALUE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0004', 'INSERTION_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0004', 'ibmqrepVERTIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('PROCESS_PARAMS0004', 'ibmqrepVERNODE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SEGMENT0004', 'LDBD', 'SEGMENT', 'ASN.QM1_TO_QM3.DATAQ', 'P', 'N',
 'N', 'N', 'I', 'I', '000001', 2, 3, 'NNNN', 'N', 'SEG6_CIT',
 'SEG6_CIT', 'LDBD', 'SEGMENT', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0004', 'CREATOR_DB', 2, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0004', 'PROCESS_ID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0004', 'SEGMENT_ID', 1, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0004', 'START_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0004', 'END_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0004', 'INSERTION_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0004', 'SEGMENT_DEF_ID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0004', 'SEGMENT_DEF_CDB', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0004', 'ibmqrepVERTIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT0004', 'ibmqrepVERNODE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SEGMENT_DEFINER0004', 'LDBD', 'SEGMENT_DEFINER',
 'ASN.QM1_TO_QM3.DATAQ', 'P', 'N', 'N', 'N', 'I', 'I', '000001', 2, 3,
 'NNNN', 'N', 'SEG6_CIT', 'SEG6_CIT', 'LDBD', 'SEGMENT_DEFINER', 1,
 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0004', 'CREATOR_DB', 2, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0004', 'PROCESS_ID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0004', 'SEGMENT_DEF_ID', 1, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0004', 'IFOS', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0004', 'NAME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0004', 'VERSION', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0004', 'COMMENT', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0004', 'INSERTION_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0004', 'ibmqrepVERTIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_DEFINER0004', 'ibmqrepVERNODE', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SEGMENT_SUMMARY0004', 'LDBD', 'SEGMENT_SUMMARY',
 'ASN.QM1_TO_QM3.DATAQ', 'P', 'N', 'N', 'N', 'I', 'I', '000001', 2, 3,
 'NNNN', 'N', 'SEG6_CIT', 'SEG6_CIT', 'LDBD', 'SEGMENT_SUMMARY', 1,
 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0004', 'SEGMENT_SUM_ID', 1, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0004', 'CREATOR_DB', 2, 'YNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0004', 'START_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0004', 'END_TIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0004', 'COMMENT', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0004', 'SEGMENT_DEF_ID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0004', 'SEGMENT_DEF_CDB', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0004', 'PROCESS_ID', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0004', 'ibmqrepVERTIME', 0, 'NNNNNNNNNN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key, col_options_flag)
 VALUES
 ('SEGMENT_SUMMARY0004', 'ibmqrepVERNODE', 0, 'NNNNNNNNNN')#
-- COMMIT#