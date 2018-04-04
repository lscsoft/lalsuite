#!/usr/bin/python

import sys
from glue.segments import segment, segmentlist
from glue.ligolw import utils as ligolw_utils
from glue.segmentdb import segmentdb_utils


def vn_ldas_gap_check(vn_xml_file,vn_flag,ldas_gap_xml_file):
    # read in summary intervals of vn_flag
    vndoc = ligolw_utils.load_url(vn_xml_file, gz = (vn_xml_file or "stdin").endswith(".gz"))
    vn_sums=segmentdb_utils.find_segments(vndoc,"%s" % vn_flag,use_segment_table=False)

    # read in segment and summary intervals of LDAS allowed gaps
    ldasdoc = ligolw_utils.load_url(ldas_gap_xml_file, gz = (ldas_gap_xml_file or "stdin").endswith(".gz"))
    ldas_gaps=segmentdb_utils.find_segments(ldasdoc,'H1:RESULT:1')
    ldas_sums=segmentdb_utils.find_segments(ldasdoc,'H1:DCH-MISSING_LDAS_C02_L2:1',use_segment_table=False)
    if len(ldas_sums) == 0:
       message = "no LDAS info was found within the specified time in the database!"
       print message
       return

    # if valid Vn file and ldas files both found, do further testings:
    if len(vn_sums) > 1 :
       vn_sum_gaps=segmentlist([segment(vn_sums[0][0],vn_sums[-1][1])]) - vn_sums
       if vn_sum_gaps == ldas_gaps:
          message = ""
       else:
          message = "test failed. Found unexpected gaps in LDAS frames: %s \n" % vn_sum_gaps - ldas_gaps
          print message 
          return
    else: 
       if (len(vn_sums & ldas_gaps) == 0):
          message = ""
       else:
          message = "test failed. Found overlap between %s summary intervals and LDAS allowed gaps: \n" % (vn_flag, vn_sums - ldas_gaps)
          print message
          return
    print message


if __name__ == "__main__":
   vn_ldas_gap_check()
   sys.exit(0)
