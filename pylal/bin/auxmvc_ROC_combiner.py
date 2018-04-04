#!/usr/bin/python

from optparse import *
import matplotlib
matplotlib.use('Agg')
import pylab
import pickle
import numpy

parser=OptionParser(usage="""
Supply several pickle files with false alarm percentages v. true alarm percentages, this code will plot them all on the same image
Usage: 
auxmvc_ROC_combiner.py --tag S6_all_data --output-dir auxmvc_results_plots file1:label1 file2:label2 file3:label3 ...
""", version="Kari Hodge")
parser.add_option("","--tag", help="filenames will be combined_ROC_tag.png")
parser.add_option("","--title", help="make an informative title here")
parser.add_option("","--output-dir", help="directory where output files will be written to")

(opts,args)=parser.parse_args()

picklefiles=[]
labels=[]
for arg in args:
	picklefiles.append(arg.split(':')[0])
	labels.append(arg.split(':')[1])

fig = pylab.figure(1)
ax = pylab.axes()
matplotlib.rcParams.update({'legend.fontsize': 12})
ax.set_color_cycle(['b','c','g','y',(0.9,0.5,0.2),'r','m',(0.5,0.,0.8),'k'])

min_eff = 1.0
min_fap = 1.0

for file in picklefiles:
	data = pickle.load(open(file))
	pylab.loglog(data[0],data[1])
        min_eff = min(min_eff, min(data[1]))
        min_fap = min(min_fap, min(data[0]))
pylab.loglog(numpy.arange(min_fap,1.0, 0.01), numpy.arange(min(data[0]),1.0, 0.01), linestyle="--", color="k",label='random')
pylab.xlabel('False Alarm Probability')
pylab.ylabel('Efficiency')
pylab.xlim([min_fap,1])
pylab.ylim([min_eff,1])
pylab.hold(True)
pylab.legend(labels,loc=4)
pylab.title(opts.title)
#pylab.text(0.005,1.0,typelabel,horizontalalignment='center')
pylab.savefig(opts.output_dir+'/loglog_combined_ROC_'+opts.tag+'.png')

fig = pylab.figure(2)
ax = pylab.axes()
matplotlib.rcParams.update({'legend.fontsize': 12})
ax.set_color_cycle(['b','c','g','y',(0.9,0.5,0.2),'r','m',(0.5,0.,0.8),'k'])

for file in picklefiles:
        data = pickle.load(open(file))
        pylab.semilogx(data[0],data[1])
pylab.semilogx(numpy.arange(min_fap,1.0, 0.01), numpy.arange(min(data[0]),1.0, 0.01), linestyle="--", color="k",label='random')
pylab.xlabel('False Alarm Probability')
pylab.ylabel('Efficiency')
pylab.xlim([min_fap,1])
pylab.ylim([min_eff,1])
pylab.hold(True)
pylab.legend(labels,loc=2)
pylab.title(opts.title)
#pylab.text(0.005,1.0,typelabel,horizontalalignment='center')
pylab.savefig(opts.output_dir+'/combined_ROC_'+opts.tag+'.png')


