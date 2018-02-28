-- This file is part of the Grid LSC User Environment (GLUE)
-- 
-- GLUE is free software: you can redistribute it and/or modify it under the
-- terms of the GNU General Public License as published by the Free Software
-- Foundation, either version 3 of the License, or (at your option) any later
-- version.
-- 
-- This program is distributed in the hope that it will be useful, but WITHOUT
-- ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
-- FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
-- details.
-- 
-- You should have received a copy of the GNU General Public License along
-- with this program.  If not, see <http://www.gnu.org/licenses/>.

--Beginning of script 1--   DatabaseDB2LUOW (SEG_LHO) [WARNING***Please do not alter this line]--
-- CONNECT TO SEG_LHO USER XXXX using XXXX;
INSERT INTO ASN.IBMQREP_SIGNAL
 (signal_type, signal_subtype, signal_input_in)
 VALUES
 ('CMD', 'CAPSTART', 'PROCESS0002');
INSERT INTO ASN.IBMQREP_SIGNAL
 (signal_type, signal_subtype, signal_input_in)
 VALUES
 ('CMD', 'CAPSTART', 'PROCESS_PARAMS0002');
INSERT INTO ASN.IBMQREP_SIGNAL
 (signal_type, signal_subtype, signal_input_in)
 VALUES
 ('CMD', 'CAPSTART', 'SEGMENT_DEFINER0002');
INSERT INTO ASN.IBMQREP_SIGNAL
 (signal_type, signal_subtype, signal_input_in)
 VALUES
 ('CMD', 'CAPSTART', 'SEGMENT0002');
INSERT INTO ASN.IBMQREP_SIGNAL
 (signal_type, signal_subtype, signal_input_in)
 VALUES
 ('CMD', 'CAPSTART', 'SEGMENT_SUMMARY0002');
-- COMMIT;
