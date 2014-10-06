# from C. Pankow

import sys, os
import glob
import re
from pylal.dq.dqTriggerUtils import fromkwfile
from glue.ligolw.utils import write_filename, write_fileobj
from glue.ligolw.ligolw import LIGO_LW, Document
from glue.ligolw import lsctables
from glue.ligolw import ilwd

if len(glob.glob("%s/*.trg" % sys.argv[1])) == 0:
	print "No .trg files found. Nothing to convert."
	sys.exit(0)

for f in glob.glob("%s/*.trg" % sys.argv[1]):
	#print f
	xmldoc = Document()
	xmldoc.appendChild( LIGO_LW() )
	sbt = lsctables.New(lsctables.SnglBurstTable,
            ["ifo", "peak_time", "peak_time_ns", "event_id", "process_id",
			"start_time", "start_time_ns", "confidence", "chisq", "chisq_dof",
			"amplitude", 
            "duration",  "search", "central_freq", "channel", "snr",
            "bandwidth"])
	#H1_TCS-ITMY_PD_ISS_OUT_AC_1_1024.xml
	fspl = os.path.basename(f).split("_")
	ifo = fspl[0]
	channel = "_".join( fspl[1:-2] )
	sbt += fromkwfile( f, ifo=ifo, channel=channel, columns = ["duration", "start_time", "peak_time", "central_freq", "bandwidth", "snr", "confidence"] )
	for i, sb in enumerate(sbt): 
		sb.search = "KleineWelle"
		sb.process_id = ilwd.ilwdchar("process:process_id:0")
		sb.event_id = ilwd.ilwdchar("sngl_burst:event_id:%d"% i )
		#sb.confidence = 0
		sb.chisq_dof = 0
		sb.chisq = 0
		sb.amplitude = 0
	xmldoc.childNodes[0].appendChild( sbt )
	write_filename( xmldoc, re.sub( "trg", "xml", f ) )
	#write_fileobj( xmldoc, sys.stdout )

