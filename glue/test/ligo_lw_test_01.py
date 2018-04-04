import matplotlib
matplotlib.use("Agg")
from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import numpy
import sys

from glue.ligolw import ligolw
from glue.ligolw import array as ligolw_array
from glue.ligolw import param as ligolw_param
from glue.ligolw import utils as ligolw_utils

class ContentHandler(ligolw.LIGOLWContentHandler):
	pass
ligolw_array.use_in(ContentHandler)
ligolw_param.use_in(ContentHandler)

xmldoc = ligolw_utils.load_filename("ligo_lw_test_01.xml", contenthandler = ContentHandler, verbose = True)
ligolw_utils.write_filename(xmldoc, "/dev/null")

t, = xmldoc.getElementsByTagName(ligolw.Time.tagName)
print >>sys.stderr, "%s: %s" % (t.Name, t.pcdata)

for n, a in enumerate(xmldoc.getElementsByTagName(ligolw.Array.tagName)):
	print >>sys.stderr, "found %s array '%s'" % ("x".join(map(str, a.array.shape)), a.Name)
	fig = figure.Figure()
	FigureCanvas(fig)
	axes = fig.gca()
	axes.loglog()
	axes.grid(True)
	for i in range(1, a.array.shape[0]):
		axes.plot(numpy.fabs(a.array[0]), numpy.fabs(a.array[i]))
	axes.set_title(a.Name)
	print >>sys.stderr, "saving as 'ligo_lw_test_01_%d.png' ..." % n
	fig.savefig("ligo_lw_test_01_%d.png" % n)
	print >>sys.stderr, "done."

	# try turning it back into XML
	ligolw_array.from_array(a.Name, a.array)
