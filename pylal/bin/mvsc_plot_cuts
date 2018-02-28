#!/usr/bin/python
"""Plots the cuts made on any two dimensions by SprBaggerDecisionTreeApp."""

__author__ = 'Tristan Miller <tmiller@caltech.edu>'
__prog__ = 'mvsc_plot_cuts'

##############################################################################
# import modules
import math, re, sys
from optparse import OptionParser
from mvsc_plots import patread
from matplotlib import pyplot

##############################################################################
# parse arguments

parser = OptionParser()
parser.add_option("-f","--filepath",help="The filepath of the test results")
parser.add_option("-z","--zeropath",default=False, \
                  help="The filepath of zero-lag test results (optional)")
parser.add_option("-t","--treepath",help="The filepath of the decision trees")
parser.add_option("-x",help="Name of the x dimension")
parser.add_option("-y",help="Name of the y dimension")
parser.add_option("-p",default=False,action="store_true",\
                  help="Print out the choices of dimensions and quit")
parser.add_option("","--xlin",action="store_true",default=False, \
                  help="Make linear in x dimension (log by default)")
parser.add_option("","--ylin",action="store_true",default=False, \
                  help="Make linear in y dimension (log by default)")
parser.add_option("-c","--cuts",action="store_true",default=False, \
                  help="Create a plot with all 1D cuts shown")
parser.add_option("-r","--curve",action="store_true",default=False, \
    help="Create a plot with 2d curve cuts shown.") 

opts,args = parser.parse_args()

##############################################################################

def find_cuts(filename,dim1_header,dim2_header):
    """Generates a list of all the cuts made on the two given dimensions."""
    
    try:
        f = open(filename)
    except IOError:
        print '***Error!*** Trouble opening file', filename
        return

    p = re.compile(r'Id: \S+ Score: \S+ Dim: (\S+) Cut: (\S+)')
    p2 = re.compile(r'Dimensions:')
    p3 = re.compile(r'\s+(\S+)\s+(\S+)')

    dim1cuts = []
    dim2cuts = []
    cutcols = {}

    while True:
        if p2.match(f.readline()):
            break

    while True:
        n = f.readline()
        m = p3.match(n)
        if m:
            cutcols[m.group(2)] = m.group(1)
        else:
            break

    f.seek(0)
    if cutcols.has_key(dim1_header):
        cutcol1 = cutcols[dim1_header]
    else:
        cutcol1 = -2
    if cutcols.has_key(dim2_header):
        cutcol2 = cutcols[dim2_header]
    else:
        cutcol2 = -2
        
    while True:
        n = f.readline()
        m = p.match(n)
        if m:
            if m.group(1) == cutcol1:
                dim1cuts.append(float(m.group(2)))
            elif m.group(1) == cutcol2:
                dim2cuts.append(float(m.group(2)))
        elif p2.match(n):
            break
        elif not n:
            print '***Error!*** Unexpected format in',filename
            return

    f.close()
        
    return dim1cuts, dim2cuts

##############################################################################

def separate_triggers(data,cols,dim1_header,dim2_header):
    """Separates triggers into timeslides and injections, only keeping the
    given dimensions."""
    
    if not (cols.has_key(dim1_header) & cols.has_key(dim2_header)):
        print 'Invalid dimensions given to plot-cuts option'
        
        header = patread( filename,headeronly=True )
        print 'Dimensions available:'
        for i in range(len(header)):
            print header[i]
            
        sys.exit()
        
    dim1 = cols[dim1_header]
    dim2 = cols[dim2_header]

    inj = [[],[]]
    ts = [[],[]]

    #separate into timeslides and injections
    for i in range(len(data[0])):
        if data[1][i] == 0:
            ts[0].append( data[dim1][i] )
            ts[1].append( data[dim2][i] )
        else:
            inj[0].append(data[dim1][i])
            inj[1].append(data[dim2][i])

    return inj,ts

##############################################################################
# The main function: plot_cuts

def plot_cuts(inj,ts,dim1cuts,dim2cuts,dim1_header,dim2_header,dim1log, \
              dim2log, zerodata=None):
    """Attempts to plot the decision tree cuts in two dimensions.

    Filename is the path to decision tree file.  If not given, will not plot
    any cuts."""

    #Plot injections and timeslides
    pyplot.figure()
    pyplot.plot(ts[0],ts[1],'xk',label='Timeslides')
    pyplot.plot(inj[0],inj[1],'+r',mec='r',label='Injections')

    if zerodata:
        zero = [[],[]]
        zerodim2 = []
        for i in range(len(zerodata[0])):
            zero[0].append( zerodata[dim1][i])
            zero[1].append( zerodata[dim2][i])

        pyplot.plot(zero[0],zero[1],'.g',mec='g', \
            label='Zero lag')
    
    #pyplot.legend(loc='lower right')
    pyplot.xlabel(dim1_header)
    pyplot.ylabel(dim2_header)
    
    if dim1log & dim2log:
        pyplot.loglog()
    elif dim1log:
        pyplot.semilogx()
    elif dim2log:
        pyplot.semilogy()
        
    xmin,xmax = pyplot.xlim()
    ymin,ymax = pyplot.ylim()
    
    #Plot cuts
    for i in range(len(dim1cuts)):
        pyplot.plot([dim1cuts[i],dim1cuts[i]],[ymin,ymax],'b',alpha=0.2)
    for i in range(len(dim2cuts)):
        pyplot.plot([xmin,xmax],[dim2cuts[i],dim2cuts[i]],'b',alpha=0.2)

    pyplot.xlim(xmin,xmax)
    pyplot.ylim(ymin,ymax)
    pyplot.title('Decision tree cuts on "'+dim1_header+'" and "' \
                 +dim2_header+'" dimensions' )

##############################################################################

def curve_cut(inj,ts,dim1cuts,dim2cuts,dim1_header,dim2_header,dim1log,\
              dim2log,pos_slope,zerodata=None):
    """As opposed to plot_cuts, which draws multiple rectangular lines, this
    function attempts to draw a curve which represents what cut the trees are
    making in two dimensions.

    This plot is entirely experimental at this point."""

    numcurves = 10

    n1 = len(dim1cuts)
    n2 = len(dim2cuts)
    if n1 == 0:
        print "There are no cuts on '" + dim1_header + "'.  Aborting..."
        return
    if n2 == 0:
        print "There are no cuts on '" + dim2_header + "'.  Aborting..."
        return
    
    dim1cuts.sort()
    dim2cuts.sort()
    if not pos_slope:
        dim1cuts.reverse()
    
    pyplot.figure()
    pyplot.plot(ts[0],ts[1],'xk',label='Timeslides')
    pyplot.plot(inj[0],inj[1],'+r',mec='r',label='Injections')

    if zerodata:
        zero = [[],[]]
        zerodim2 = []
        for i in range(len(zerodata[0])):
            zero[0].append( zerodata[dim1][i])
            zero[1].append( zerodata[dim2][i])

        pyplot.plot(zero[0],zero[1],'.g',mec='g', \
            label='Zero lag')
    
    #pyplot.legend(loc='lower right')
    pyplot.xlabel(dim1_header)
    pyplot.ylabel(dim2_header)
    
    if dim1log & dim2log:
        pyplot.loglog()
    elif dim1log:
        pyplot.semilogx()
    elif dim2log:
        pyplot.semilogy()
    
    for i in range(1,numcurves+1):
        k = int( i*(n1 + n2)/float(numcurves+1) ) - n1

        curvecut = [[],[]]
        
        if k < 0:
            curvecut[0] = dim1cuts[-k:]
            curvecut[1] = dim2cuts[:len(curvecut[0])]
            curvecut[0] = curvecut[0][:len(curvecut[1])]
        else:
            curvecut[1] = dim2cuts[k:]
            curvecut[0] = dim1cuts[:len(curvecut[1])]
            curvecut[1] = curvecut[1][:len(curvecut[0])]

        pyplot.plot(curvecut[0],curvecut[1],'b')
            
    pyplot.title('2D representation of cuts on "'+dim1_header+'" and "' \
                 +dim2_header+'" dimensions' )
    
##############################################################################
# execute plot_cuts

if not opts.filepath:
    print 'Filepath option (-f) required'
    sys.exit()

if opts.p:
    header = patread( opts.filepath,headeronly=True )
    print 'Dimensions available:'
    for i in range(len(header)):
        print header[i]
    sys.exit()

if not opts.treepath:
    print 'Treepath option (-t) required'
    sys.exit()
elif not opts.x:
    print 'X dimension option (-x) required'
    sys.exit()
elif not opts.y:
    print 'Y dimension option (-y) required'
    sys.exit()

data,cols = patread( opts.filepath )

if opts.zeropath:
    zerodata,temp1 = patread( opts.zeropath )
else:
    zerodata = None

inj,ts = separate_triggers(data,cols,opts.x,opts.y)
dim1cuts,dim2cuts = find_cuts( opts.treepath, opts.x, opts.y )

if opts.cuts:
    plot_cuts(inj,ts,dim1cuts,dim2cuts,opts.x,opts.y,not opts.xlin, \
              not opts.ylin, zerodata)
    
    ymin,ymax = pyplot.ylim()
    pyplot.ylim(10.**(-8.),ymax)

if opts.curve:
    curve_cut(inj,ts,dim1cuts,dim2cuts,opts.x,opts.y,not opts.xlin, \
              not opts.ylin,True,zerodata)
    
    ymin,ymax = pyplot.ylim()
    pyplot.ylim(10.**(-8.),ymax)
    
    curve_cut(inj,ts,dim1cuts,dim2cuts,opts.x,opts.y,not opts.xlin, \
              not opts.ylin,False,zerodata)
    
    ymin,ymax = pyplot.ylim()
    pyplot.ylim(10.**(-8.),ymax)

pyplot.show()
