--Beginning of script 1--   DatabaseDB2LUOW (SEG6_LHO) [WARNING***Please do not alter this line]--
-- CONNECT TO SEG6_LHO USER XXXX using XXXX;
INSERT INTO ASN.IBMQREP_SENDQUEUES
 (pubqmapname, sendq, message_format, msg_content_type, state,
 error_action, heartbeat_interval, max_message_size, description,
 apply_alias, apply_schema, recvq, apply_server, sendraw_iferror)
 VALUES
 ('SEG6_LHO_ASN_TO_SEG6_CIT_ASN', 'ASN.QM1_TO_QM3.DATAQ', 'C', 'T',
 'A', 'S', 60, 64, '', 'SEG6_CIT', 'ASN', 'ASN.QM1_TO_QM3.DATAQ',
 'SEG6_CIT', 'N');
-- COMMIT;
