#!/usr/bin/env python
import sys, getopt
import matplotlib.cm
from matplotlib.patches     import Patch
from matplotlib.axes        import Axes
from matplotlib.collections import PolyCollection
from matplotlib.colors      import normalize, Colormap
from optparse import * 
from pylab    import *
import numpy
from numpy import power
from pylal import rate

#####################################################################
# use tex labels
params =  {'text.usetex': True }
rcParams.update(params)


def square_axis(axes=None):
  """
  Expands the x- and y-limits on the given axes so that they are equal
  (defaults to the current axes).
  """
  if axes is None:
    axes = gca()
  tmpv = axes.axis()
  xmax = max([tmpv[1], tmpv[3]])
  xmin = min([tmpv[0], tmpv[2]])
  axes.axis([xmin, xmax, xmin, xmax])

def simpleplot(*args):
  if len(args)==3:
      mytable, col1, col2 = args
  else:
      raise TypeError, 'Illegal arguments to simpleplot; see help(rectfill)'
   
  if mytable is not None: assert(isinstance(mytable, metaDataTable))
  tmpvar1 = mytable.mkarray(col1)
  tmpvar2 = mytable.mkarray(col2)
  plot(tmpvar1,tmpvar2,'rx')
  xlabel(col1.replace("_"," "), size='x-large')
  ylabel(col2.replace("_"," "), size='x-large')
  title(col1.replace("_"," ") + ' v ' + col2.replace("_"," "))
  grid(True)

#######################################################################
# function to read in a column from a given table 
def readcol(table, col_name, ifo=None ):
  """
  function to read in a column from a given table.  If the column is ifo
  dependent (eg, end_time or eff_dist) then the ifo is used to select the
  appropriate value.

  this function can also read in two values not stored in the sngl_inspiral
  table:
    - snr_chi: read in value of snr/chi
    - s3_snr_chi_stat: read in value of snr^4 / ( chisq ( snr^2 + 250 ))

  @param table: metaDataTable
  @param col_name: name of column to read in 
  @param ifo: name of ifo (used to extract the appropriate end_time/eff_dist)
  """
  
  col_names = [col_name]
  if ifo:
    col_names.append(col_name + '_' + ifo[0].lower())
    col_names.append(ifo[0].lower() + '_' + col_name)

  if 'dist' in col_name:
    col_names.append(col_name + 'ance')

  col_data = []
  
  if len(table):
    for c_name in col_names:
      try: 
        col_data = table.get_column(c_name)
        if 'time' in c_name:
          col_data = col_data + 1e-9 * table.get_column(c_name + '_ns')
      except: x = 0
  return col_data

#######################################################################
# function to read in a column from two tables 
def readcolfrom2tables(table1, table2, col_name ):
  """
  function to read in a column from a two tables.  If the column is ifo
  dependent (eg, end_time or eff_dist) then the ifo of one table is used to 
  select the appropriate value in the other table.  The number of events in
  the two tables must be equal.
  
  @param table1: metaDataTable
  @param table2: metaDataTable
  @param col_name: name of column to read in 
  """
  
  if len(table1) != len(table2):
    raise ValueError, "number of events in table1 and table2 must be equal"
 
  if len(table1):
    if ("ifo" in table1.validcolumns.keys()):
      ifo = table1[0].ifo
    elif ("ifo" in table2.validcolumns.keys()):
      ifo = table2[0].ifo
    else:
      ifo = None
  else:
    ifo = None

  col1 = readcol(table1, col_name, ifo)
  col2 = readcol(table2, col_name, ifo)

  cols = []
  cols.append(col1)
  cols.append(col2)
  cols.append(ifo)

  return cols


#######################################################################
# function to read in a column from two tables 
def timeindays(col_data ):
  """
  function to re-express the time in days after the start of the run
  known runs:
    - ligo_virgo: [700000000, 700086400]
    - S2:         [729273613, 734367613]
    - S3:         [751658413, 757699213]
    - S4:         [793130413, 795679213]
    - S5:         [815119213, 875232014]
    - E13:        [924606015, 924865215]
    - E14:        [928875615, 929134815]
    - S6:         [930960015, ?????????]
  @param col_data: array containing times in GPS seconds
  """
  lvtimes = [700000000, 700086400]
  v1times = [811132263, 811143059]
  s2times = [729273613, 734367613]
  s3times = [751658413, 757699213]
  s4times = [793130413, 795679213]
  s5times = [815119213, 875232014]
  ninjatimes = [900000000, 905000000]
  e13times = [924606015, 924865215]
  e14times = [928875615, 929134815]
  s6time = 930960015

  if len(col_data) == 0: return col_data

  if col_data[0] > s6time:
    start = s6time
  elif col_data[0] > s2times[0] and col_data[0] < s2times[1]:
    start = s2times[0]
  elif col_data[0] > s3times[0] and col_data[0] < s3times[1]:
    start = s3times[0]
  elif col_data[0] > s4times[0] and col_data[0] < s4times[1]:
    start = s4times[0]
  elif col_data[0] > s5times[0] and col_data[0] < s5times[1]:
      start = s5times[0]
  elif col_data[0] > lvtimes[0] and col_data[0] < lvtimes[1]:
    start = lvtimes[0]
  elif col_data[0] > v1times[0] and col_data[0] < v1times[1]:
    start = v1times[0]
  elif col_data[0] > e13times[0] and col_data[0] < e13times[1]:
    start = e13times[0]
  elif col_data[0] > e14times[0] and col_data[0] < e14times[1]:
    start = e14times[0]
  elif col_data[0] > ninjatimes[0] and col_data[0] < ninjatimes[1]:
    start = ninjatimes[0]
  else:
    raise ValueError, "events not from a known science run"

  col_data = (col_data - start)/(60 * 60 * 24.0)

  return col_data

#############################################################################
# make steps so that fill will work properly
def makesteps(left_edges, y1, y2=None):
  """
  Return ordinates and coordinates for the vertices of a polygon with
  "stairsteps" with upper heights y1 and lower heights y2 (zero if y2 is
  None).  The right edge of the polygon is extrapolated from the left_edges,
  assuming uniform spacing.
  """
  left_edges = numpy.asanyarray(left_edges)
  n = len(left_edges)
  if len(y1) != n or (y2 is not None and len(y2) != n):
    raise ValueError, "left_edges, y1, and y2 (if specified) must have "\
      "matching lengths"
  if n < 2:
    raise ValueError, "don't know how to form stairs with less than 2 points"

  if y2 is None: y2 = y1
  y1 = numpy.asanyarray(y1)
  y2 = numpy.asanyarray(y2)
  
  right_edge = left_edges[-1] + (left_edges[-1] - left_edges[-2])
  
  # fill x
  x_new = numpy.empty(4*n, dtype=left_edges.dtype)
  x_new[:2*n-1:2] = left_edges
  x_new[1:2*n-1:2] = left_edges[1:]
  x_new[2*n-1] = right_edge
  x_new[2*n:] = x_new[2*n-1::-1]  # make it a palindrome
  
  # fill y
  y_new = numpy.empty(4*n, dtype=min(y1.dtype, y2.dtype))
  y_new[:2*n:2] = y1
  y_new[1:2*n:2] = y1
  y_new[2*n::2] = y2[::-1]
  y_new[2*n+1::2] = y2[::-1]
  
  return x_new, y_new



#######################################################################
# function to plot the col1 vs col2 from the table
def plot_a_v_b(table, col_name_a, col_name_b, plot_type = 'linear', 
  plot_sym = 'kx', plot_label = None, output_name = None, ifo = None,
  x_min = None, x_max = None, y_min = None, y_max = None):
  """
  function to plot the values of col_name_a vs col_name_b from the table

  @param table: metaDataTable
  @param col_name_a: name of first column (x-axis)
  @param col_name_b: name of second column (y-axis)
  @param plot_type: One of 'linear' (default) ,'logx','logy','loglog','seconds'
  @param plot_sym : Symbol for plot, default 'kx'
  @param plot_label: Name for the plot in the legend
  @param output_name: If given, save the plot, with output_name as prefix
  @param ifo: Name of ifo
  @param x_min: Minimum value of x axis
  @param x_max: Maximum value of x axis
  @param y_min: Minimum value of y axis
  @param y_max: Maximum value of y axis
  """
  if not ifo and len(table):
    if ("ifo" in table.validcolumns.keys()):
      ifo = table[0].ifo
  
  col_a = readcol(table, col_name_a, ifo )
  col_b = readcol(table, col_name_b, ifo )

  if plot_type != 'seconds':
    if 'time' in col_name_a:
      col_a = timeindays(col_a)
    if 'time' in col_name_b:
      col_b = timeindays(col_b)
   
  if (plot_type == 'linear') or (plot_type == 'seconds'):
    plot(col_a, col_b, plot_sym, markersize=12,markeredgewidth=1,\
        markerfacecolor='None', label = plot_label)
  elif plot_type == 'logx':
    semilogx(col_a, col_b, plot_sym, markersize=12,markeredgewidth=1,\
        markerfacecolor='None', label = plot_label)
  elif plot_type == 'logy':
    semilogy(col_a, col_b, plot_sym, markersize=12,markeredgewidth=1,\
        markerfacecolor='None', label = plot_label)
  elif plot_type == 'loglog':
    loglog(col_a, col_b, plot_sym, markersize=12,markeredgewidth=1,\
        markerfacecolor='None', label = plot_label)
    xlim(0.95 * min(col_a), 1.05 * max(col_a))
    ylim(0.95 * min(col_b), 1.05 * max(col_b))

  xlabel(col_name_a.replace("_"," "),size='x-large')
  ylabel(col_name_b.replace("_"," "),size='x-large')

  xticks(fontsize='x-large')
  yticks(fontsize='x-large')

  if ifo:
    title(ifo + ' ' + col_name_b.replace("_"," ") + ' vs ' + \
        col_name_a.replace("_"," "), size='x-large')
  else:
    title(col_name_b.replace("_"," ") + ' vs ' + col_name_a.replace("_"," "), \
        size='x-large')

  grid(True)

  if plot_label:
    legend()

  if x_min:
    xlim(xmin=x_min)
  if x_max:
    xlim(xmax=x_max)
    
  if y_min:
    ylim(ymin=y_min)
  if y_max:
    ylim(ymax=y_max)

  if output_name:
    if ifo:
      output_name += '_' + ifo
    output_name += '_' + col_name_a + '_vs_' + col_name_b + '.png'
    savefig(output_name)

#################################################################
# function to plot the difference between values of 'col_name' in
# two tables, table1 and table2
def plotdiff(table1, table2, col_name, plot_type = 'linear', plot_sym = 'kx',\
    plot_name = ''):
  """
  function to plot the difference between the value of col_name stored in 2
  tables (of equal length).  
  
  @param table1: metaDataTable
  @param table2: metaDataTable
  @param col_name: name of column to plot
  @param plot_type: either 'linear' (default)or 'log' plot on x-axis
  @param plot_sym:  the symbol to use when plotting, default = 'kx'
  @param plot_name: name of the plot (for the legend)
  """

  [tmpvar1, tmpvar2, ifo ] = readcolfrom2tables(table1, table2, col_name)

  tmp_diff = tmpvar2 - tmpvar1

  if 'time' in col_name:
    tmpvar1 = timeindays(tmpvar1)

  if plot_type == 'linear':
    plot(tmpvar1, tmp_diff, plot_sym, markersize=12,markerfacecolor='None',
      markeredgewidth=1, label = plot_name)
  elif plot_type == 'log':
    semilogx(tmpvar1, tmp_diff, plot_sym, markersize=12,markerfacecolor='None',
      markeredgewidth=1, label = plot_name)
    
#################################################################
# function to label above plot
def labeldiff(col_name, units = None, axis = [0,0,0,0], leg = None, 
  title_text = None, output_name = None):
  """
  function to label the output of plotdiff
  
  @param col_name: name of column to plot
  @param units: the units of the column
  @param axis: axis limits [xmin,xmax,ymin,ymax].  If both min and max of x or
               y is zero then that axis is not set.
  @param leg: add legend to plot
  @param title_text: text to add at start of title, if given then title is
                     title_text + col_name + Accuracy
  @param output_name: used in naming output file
  """
  
  if units:
    xlabel(col_name.replace("_"," ") + ' (' + units +')', size='x-large')
    ylabel(col_name.replace("_"," ") + ' difference (' + units +')',\
        size='x-large')
  else:
    xlabel(col_name.replace("_"," "), size='x-large')
    ylabel(col_name.replace("_"," ") + ' difference', size='x-large')

  xticks(fontsize='x-large')
  yticks(fontsize='x-large')
  
  if axis[0] or axis[1]:
    xlim(axis[0], axis[1])

  if axis[2] or axis[3]:
    ylim(axis[2], axis[3])

  if leg:
    legend()
 
  grid(True)

  if title_text:
    title(title_text + ' ' + col_name.replace("_"," ") + '  Accuracy', \
        size='x-large', weight='bold')
  # else:
  #   title(col_name + ' Accuracy', size='x-large',weight='bold')
  
  if output_name:
    output_name += '_' + col_name + '_accuracy.png'
    savefig(output_name)


############################################################################
# function to plot the fractional difference between values of 'col_name' in
# two tables, table1 and table2
def plotfracdiff(table1, table2, col_name, plot_type = 'linear', 
    plot_sym = 'kx', plot_name = ''):
  """
  function to plot the fractional difference between the value of 
  col_name stored in 2 tables (of equal length).  
  
  @param table1: metaDataTable
  @param table2: metaDataTable
  @param col_name: name of column to plot
  @param plot_type: either 'linear' (default) or 'log' plot on x-axis
  @param plot_sym: the symbol to use when plotting, default = 'kx'
  @param plot_name: name of the plot (for the legend)  
  """

  [tmpvar1, tmpvar2, ifo ] = readcolfrom2tables(table1, table2, col_name)

  frac_diff = (tmpvar2 - tmpvar1)/tmpvar1

  if 'time' in col_name:
    tmpvar1 = timeindays(tmpvar1)

  if plot_type == 'linear':
    plot(tmpvar1, frac_diff,plot_sym,markersize=12,markerfacecolor='None',\
        markeredgewidth=1, label = plot_name)
  elif plot_type == 'log':
    semilogx(tmpvar1, frac_diff,plot_sym,markersize=12,markerfacecolor='None',\
        markeredgewidth=1, label = plot_name)


#################################################################
# function to label above plot
def labelfracdiff(col_name, units = None, axis = [0,0,0,0], leg = None, 
  title_text = None, output_name = None):
  """
  function to label the output of plotfracdiff
  
  @param col_name: name of column to plot
  @param units: the units of the column
  @param axis: axis limits [xmin,xmax,ymin,ymax].  If both min and max of x or
               y is zero then that axis is not set.
  @param leg: add legend to plot
  @param title_text: text to add at start of title
  @param output_name: used in naming output file
  """

  if units:
    xlabel(col_name.replace("_"," ") + ' (' + units +')', size='x-large')
  else:
    xlabel(col_name.replace("_"," "), size='x-large')
  
  ylabel(col_name.replace("_"," ") + ' fractional difference', size='x-large')

  xticks(fontsize='x-large')
  yticks(fontsize='x-large')

  if axis[0] or axis[1]:
    xlim(axis[0], axis[1])

  if axis[2] or axis[3]:
    ylim(axis[2], axis[3])

  if leg:
    legend()
 
  grid(True)

  if title_text:
    title(title_text + ' ' + col_name.replace("_"," ") + '  Accuracy', \
        size='x-large', weight='bold')
  else:
    title(col_name.replace("_"," ") + ' Accuracy', \
        size='x-large',weight='bold')
  
  if output_name:
    output_name += '_' + col_name + '_frac_accuracy.png'
    savefig(output_name)


############################################################################
# function to plot the fractional difference between values of 'col_name_a' in
# two tables, table1 and table2 against the values of 'col_name_b' in table1
def plotdiffa_vs_b(table1, table2, col_name_a, col_name_b, \
    plot_type = 'linear', plot_sym = 'kx'):
  """
  function to plot the difference if col_name_a in two tables against the
  value of col_name_b in table1.  
  
  @param table1: metaDataTable
  @param table2: metaDataTable
  @param col_name_a: name of column to plot difference of on y-axis
  @param col_name_b: name of column to plot on x-axis
  @param plot_type: either 'linear' or 'log' plot on x-axis
  @param plot_sym: the symbol to use when plotting
  """
  
  [tmpvar1, tmpvar2, ifo ] = readcolfrom2tables(table1, table2, col_name_a)

  diff_a = (tmpvar2 - tmpvar1)
  col_b = readcol(table1, col_name_b, ifo ) 

  if 'time' in col_name_b:
    col_b = timeindays(col_b)

  if plot_type == 'linear':
    plot(col_b, diff_a,plot_sym,markersize=12,markerfacecolor='None',
      markeredgewidth=1)
  elif plot_type == 'log':
    semilogx(col_b, diff_a,plot_sym,markersize=12,markerfacecolor='None',
      markeredgewidth=1)
  
  
#################################################################
# function to label above plot
def labeldiffa_vs_b(col_name_a, col_name_b, units_a=None, units_b=None,
  axis=[0,0,0,0], leg = None, title_text = None, output_name = None):
  """
  function to label the output of plotdiffa_vs_b
  
  @param col_name_a: name of column to plot
  @param col_name_b: name of column to plot
  @param units_a: the units of column a
  @param units_b: the units of column b
  @param axis: axis limits [xmin,xmax,ymin,ymax].  If both min and max of x or
               y is zero then that axis is not set.
  @param leg: legend to add to plot
  @param title_text: text to add at start of title
  @param output_name: used in naming output file
  """
    
  if units_b:
    xlabel(col_name_b.replace("_"," ") + ' (' + units_b +')', size='x-large')
  else:
    xlabel(col_name_b.replace("_"," "), size='x-large')
  
  if units_a:
    ylabel(col_name_a.replace("_"," ") + ' difference (' + units_a + ')', \
        size='x-large')
  else:
    ylabel(col_name_a.replace("_"," ") + ' difference', size='x-large')

  xticks(fontsize='x-large')
  yticks(fontsize='x-large')

  if axis[0] or axis[1]:
    xlim(axis[0], axis[1])

  if axis[2] or axis[3]:
    ylim(axis[2], axis[3])

  if leg:
    legend(leg)
 
  grid(True)

  if title_text:
    title(title_text + ' ' + col_name_a.replace("_"," ") + ' Accuracy vs ' + \
        col_name_b.replace("_"," "), \
        size='x-large', weight='bold')
  else:
    title(col_name_a.replace("_"," ") + ' Accuracy vs ' + \
        col_name_b.replace("_"," "),
        size='x-large',weight='bold')

  if output_name:
    output_name += '_' + col_name_a + '_vs_' + col_name_b + '_accuracy.png'
    savefig(output_name)

############################################################################
# function to plot the fractional difference between values of 'col_name_a' in
# two tables, table1 and table2 against the values of 'col_name_b' in table1
def plotfracdiffa_vs_b(table1, table2, col_name_a, col_name_b, \
    plot_type = 'linear', plot_sym = 'kx'):
  """
  function to plot the difference if col_name_a in two tables against the
  value of col_name_b in table1.

  @param table1: metaDataTable
  @param table2: metaDataTable
  @param col_name_a: name of column to plot difference of on y-axis
  @param col_name_b: name of column to plot on x-axis
  @param plot_type: either 'linear' or 'log' plot on x-axis
  @param plot_sym: the symbol to use when plotting
  """

  [tmpvar1, tmpvar2, ifo ] = readcolfrom2tables(table1, table2, col_name_a)

  diff_a = 2.*(tmpvar2 - tmpvar1)/(tmpvar2 + tmpvar1)
  col_b = readcol(table1, col_name_b, ifo )

  if 'time' in col_name_b:
    col_b = timeindays(col_b)

  if plot_type == 'linear':
    plot(col_b, diff_a,plot_sym,markersize=12,markerfacecolor='None',
      markeredgewidth=1)
  elif plot_type == 'log':
    semilogx(col_b, diff_a,plot_sym,markersize=12,markerfacecolor='None',
      markeredgewidth=1)


#################################################################
# function to label above plot
def labelfracdiffa_vs_b(col_name_a, col_name_b, units_a=None, units_b=None,
  axis=[0,0,0,0], leg = None, title_text = None, output_name = None):
  """
  function to label the output of plotfracdiffa_vs_b

  @param col_name_a: name of column to plot
  @param col_name_b: name of column to plot
  @param units_a: the units of column a
  @param units_b: the units of column b
  @param axis: axis limits [xmin,xmax,ymin,ymax].  If both min and max of x or
               y is zero then that axis is not set.
  @param leg: legend to add to plot
  @param title_text: text to add at start of title
  @param output_name: used in naming output file
  """

  if units_b:
    xlabel(col_name_b.replace("_"," ") + ' (' + units_b +')', size='x-large')
  else:
    xlabel(col_name_b.replace("_"," "), size='x-large')

  ylabel(col_name_a.replace("_"," ") + ' fractional difference', \
      size='x-large')

  xticks(fontsize='x-large')
  yticks(fontsize='x-large')

  if axis[0] or axis[1]:
    xlim(axis[0], axis[1])

  if axis[2] or axis[3]:
    ylim(axis[2], axis[3])

  if leg:
    legend(leg)

  grid(True)

  if title_text:
    title(title_text + ' ' + col_name_a.replace("_"," ") + \
        ' Fractional Accuracy vs ' + \
        col_name_b.replace("_"," "), \
        size='x-large', weight='bold')
  else:
    title(col_name_a.replace("_"," ") + ' Fractional Accuracy vs ' + \
        col_name_b.replace("_"," "),
        size='x-large',weight='bold')

  if output_name:
    output_name += '_' + col_name_a + '_vs_' + col_name_b + '_frac_accuracy.png'
    savefig(output_name)
 
###################################################
# function to plot the value of 'col_name' in table1 vs its
# value in table2 
def plotval(table1, table2, col_name, plot_type, plot_sym):
  """
  function to plot the value of col_name in table1 vs its value in table2  
  
  @param table1: metaDataTable
  @param table2: metaDataTable
  @param col_name: name of column to plot
  @param plot_type: either 'linear' or 'log' plot on x-axis
  @param plot_sym: the symbol to use when plotting
  """
  
  [tmpvar1, tmpvar2, ifo ] = readcolfrom2tables(table1, table2, col_name)

  if plot_type == 'linear':
    plot(tmpvar1, tmpvar2,plot_sym,markersize=12,markerfacecolor='None', \
      markeredgewidth=1)
  elif plot_type == 'log':
    loglog(tmpvar1, tmpvar2,plot_sym,markersize=12,markerfacecolor='None', \
      markeredgewidth=1)

   

#################################################################
# function to label above plot
def labelval(col_name, units = None, axis = [0,0,0,0], xlab = None, \
  ylab = None, leg = None, title_text = None, output_name = None):
  """
  function to label the output of plotval
  
  @param col_name: name of column to plot
  @param units: the units of the column
  @param axis: axis limits [xmin,xmax,ymin,ymax].  If both min and max of x or
               y is zero then that axis is not set.
  @param xlab: label for x-axis
  @param ylab: label for y-axis
  @param leg: legend to add to plot
  @param title_text: text to add at start of title
  @param output_name: used in naming output file
  """

  if xlab:
    xlab += ' ' + col_name
  else:
    xlab = col_name

  if ylab:
    ylab += ' ' + col_name
  else:
    ylab = col_name
    
  if units:
    xlab +=' (' + units +')'
    ylab += ' (' + units +')'
    
  xlabel(xlab.replace("_"," "), size='x-large')
  ylabel(ylab.replace("_"," "), size='x-large')
 
  if axis[0] or axis[1]:
    xlim(axis[0], axis[1])

  if axis[2] or axis[3]:
    ylim(axis[2], axis[3])

  if leg:
    legend(leg)
 
  grid(True)

  if title_text:
    title(title_text + ' ' + col_name.replace("_"," "), size='x-large', \
        weight='bold')
  else:
    title(col_name.replace("_"," "),size='x-large',weight='bold')

  if output_name:
    output_name += '_' + col_name + '_plot.png'
    savefig(output_name)
  

###########################################################
# function to plot the value of 'col_name' in ifo1  vs its
# value in ifo2 
def plotcoincval(coinctable, col_name, ifo1, ifo2, plot_sym, plot_type,\
  additional_ifo=None):
  
  if additional_ifo:
    ifo_coinc = coinctable.coinctype([ifo1,ifo2,additional_ifo])
  else:
    ifo_coinc = coinctable.coinctype([ifo1,ifo2])

  ifo1_val = readcol(ifo_coinc.getsngls(ifo1),col_name,ifo1)
  ifo2_val = readcol(ifo_coinc.getsngls(ifo2),col_name,ifo2)
  
  if plot_type == 'linear':
    plot(ifo1_val, ifo2_val,plot_sym,markersize=12,markeredgewidth=1,\
        markerfacecolor='None')
  elif plot_type == 'log':
    loglog(ifo1_val, ifo2_val,plot_sym,markersize=12,markeredgewidth=1,\
        markerfacecolor='None')


###########################################################
# function to plot the value of 'col_name' in hanford  vs its
# value in ifo.  
def plotcoinchanford(coinctable, col_name, ifo, \
  hanford_method, plot_sym,plot_type, additional_ifo=None):

  if additional_ifo:
    ifo_coinc = coinctable.coinctype([ifo,additional_ifo,'H1','H2'])
  else:
    ifo_coinc = coinctable.coinctype([ifo,'H1','H2'])

  h1_val = readcol(ifo_coinc.getsngls('H1'),col_name,'H1')
  h2_val = readcol(ifo_coinc.getsngls('H2'),col_name,'H2')
  ifo_val = readcol(ifo_coinc.getsngls(ifo),col_name,ifo)
  
  if hanford_method == 'sum':
    h_val = h1_val + h2_val
  if hanford_method == 'mean':
    h_val = (h1_val + h2_val)/2
  if hanford_method == 'sqrtsqr':
    h_val = []
    for idx in range(len(h1_val)):
      h_val.append(sqrt(h1_val[idx] * h1_val[idx] + h2_val[idx] * h2_val[idx]))
   
  if plot_type == 'linear':
    plot(h_val, ifo_val, plot_sym,markersize=12,markerfacecolor='None',\
        markeredgewidth=1)
  elif plot_type == 'log':
    loglog(h_val, ifo_val,plot_sym,markersize=12,markerfacecolor='None',\
        markeredgewidth=1)

    
######################################################################
# helper function to get the axes right when plotting logarithmic data
def setAxisForLog( ax , limit):
  """
  Helper function to get the axes right when plotting logarithmic data
  """

  # get the range
  ticks = range( int(limit[0]+0.5), int(limit[1])+1)

  # set the ticks and the ticklabels
  ax.set_ticks(ticks)
  ax.set_ticklabels(power(10.0,ticks))

######################################################################
# function to histogram the difference between values of 'col_name' in
# two tables, table1 and table2
def histcol(table1, col_name,nbins = None, width = None, output_name = None, xlimit=[0,0], plot_type='normal', color='blue'):
  """

  """ 
  if ("ifo" in table1.validcolumns.keys()):
    ifo = table1[0].ifo
  else:
    ifo = None

  data = readcol(table1, col_name, ifo )
  
  if len(data) > 0:
    if not nbins:
      nbins = 10
    
    bins = []
    if width:
      for i in range(-nbins,nbins + 1):
        bins.append(width * i/nbins)
    
    # creates the histogram and take plot_type into account  
    if plot_type == 'loglog' or plot_type=='logx':
      data = log10(data)


    if len(bins) != 0:
      ydata, xdata = numpy.histogram(data, bins)
    else:
      ydata, xdata = numpy.histogram(data, nbins)
    xdata = xdata[:-1]

    width = xdata[1] - xdata[0]

    if plot_type == 'loglog' or plot_type=='logy':
      indexPositive = find(ydata>0)
      ydata = log10( ydata[indexPositive] )
      xdata = xdata[indexPositive]

      clf()
      # ydata may be zero. if this is the case bar fails. 
      # let cheat and increment ydata by 1, 
      # which we will take into account in the tick_labels 
      bar( xdata, 1 + ydata, width ,color=color)
    else:
      clf()
      bar( xdata, ydata, width ,color=color)
  
    ax=axes()
    xlabel(col_name.replace("_"," "), size='x-large')
    ylabel('Number', size='x-large')

    # now let us set the ticks and ylabels.
    # First to be human readable (not in log10), 
    # let us come back to power of 10 instead of log10 values
    # on x, which is easy....:
    if plot_type=='logx' or plot_type=='loglog':
      locs, labels = xticks()
      l = len(locs)
      lim1 = floor(log10(power(10, locs[0])))
      lim2 = ceil(log10(power(10, locs[l-1])))
      this = range(int(lim1), int(lim2), 1)
      ticks_labels = power(10.0, [x for x in this])
      ticks = range(int(lim1), int(lim2), 1)
      this = ax.get_xaxis()
      this.set_ticks(ticks)
      this.set_ticklabels([ str(x) for x in ticks_labels])

    if plot_type=='logy' or plot_type=='loglog':
      # and on y 
      locs, labels = yticks()
      l = len(locs)
      lim1 = floor(log10(power(10, locs[0])))
      # note the +1 here (in anticipation to a shift later on
      lim2 = ceil(log10(power(10, locs[l-1]+1)))
      this = range(int(lim1), int(lim2), 1)
      # note the -1 here to compensate with the bar(xdata, 1+ydata) call
      ticks_labels = power(10.0, [x-1 for x in this])
      # finally, if a tick is 0.1, it should be 0
      if ticks_labels[0]==0.1: ticks_labels[0]  = 0

      ticks = range(int(lim1), int(lim2), 1)
      this = ax.get_yaxis()
      this.set_ticks(ticks)
      this.set_ticklabels([ str(x) for x in ticks_labels])

  #set the x axis limits taking into account if we use log10(data) 
  if xlimit[0] and  xlimit[1]:
    xlim(xlimit)
  if plot_type == 'loglog'  or plot_type=='logx':
    xlimit[0] = log10(xlimit[0])
    if xlimit[1] < xlimit[0]:
      xlimit[1] = max(data)
      xlim(xlimit)
  else:
    if xlimit[1] < xlimit[0]:
      xlimit[1] = max(data)
      xlim(xlimit)

  
  if ifo:
    title(ifo + ' ' + col_name.replace("_"," ") + ' histogram', size='x-large')
  else:
    title(col_name.replace("_"," ") + ' histogram', size='x-large')

  grid(True)

  if output_name:
    if ifo:
      output_name += '_' + ifo
    output_name += '_' + col_name + '_histogram.png'
    savefig(output_name)


######################################################################
# function to histogram the difference between values of 'col_name' in
# two tables, table1 and table2
def cumhistcol(table1, col_name, plot_type = 'logy', normalization=None, \
    output_name = None, ifo=None, xlimit = [0,0], color='k'):
 
  if not ifo and ("ifo" in table1.validcolumns.keys()):
    ifo = table1[0].ifo

  data = readcol(table1, col_name, ifo )

  data_len = len(data)
  data_sort = sort(data)
  data_range = arange(len(data))

  y_data = data_len - data_range
  if normalization:
    y_data = y_data/float(normalization)

  symbol='-'+color  
  if plot_type == 'logy':
    semilogy(data_sort, y_data,symbol,linewidth=1)
  elif plot_type == 'logx':
    semilogx(data_sort, y_data,symbol,linewidth=1)
  elif plot_type == 'loglog':
    loglog(data_sort, y_data,symbol,linewidth=1)
  else:
    plot(data_sort, y_data,symbol,linewidth=1)

  # set the x-axis limits
  if xlimit[0] and  xlimit[1]:
    xlim(xlimit)
  if xlimit[1] < xlimit[0]:
    xlimit[1] = max(data)
    xlim(xlimit)

  xlabel(col_name.replace("_"," "), size='x-large')
  
  if normalization:
    ylabel('Probability', size='x-large')
  else:  
    ylabel('Cumulative Number', size='x-large')

  if ifo:
    title_string = ifo + ' ' + col_name.replace("_"," ")
  else:
    title_string = col_name.replace("_"," ")
  if normalization:
    title_string += ' normalized'
  title_string += ' cumulative histogram'
  title(title_string, size='x-large')

  grid(True)
  if output_name:
    if ifo:
      output_name += '_' + ifo
    if normalization:
      output_name += '_' + col_name + '_norm_hist.png'
    else:
      output_name += '_' + col_name + '_cum_hist.png'
    savefig(output_name)


######################################################################
# function to make a cumulative histogram of the coinc inspiral statistic
def cumhiststat(trigs=None, slide_trigs=None,ifolist = None, min_val = None, \
  max_val = None, nbins = 20, stat=None,scalebkg=False):
  """
  function to plot a cumulative histogram of the snr of coincident events
  in the zero lag and times slides
  
  @param trigs: coincInspiralTable
  @param slide_trigs: dictionary of time slide triggers
  @param ifolist: list of ifos
  @param min_val: minimum of snr to be plotted
  @param max_val: maximum of snr to be plotted
  @param nbins: number of bins to use in histogram
  @param stat: the statistic being used
  @param scalebkg: Use this option if plotting playground zero lag against
  full data time slides (it will rescale the time slides).
  """
  internal_min = numpy.inf
  internal_max = -numpy.inf
  
  # filter down to ifolist
  slide_trig_list = [s["coinc_trigs"] for s in slide_trigs]
  if ifolist:
    if trigs:
      trigs = trigs.coinctype(ifolist)
    slide_trig_list = [s.coinctype(ifolist) for s in slide_trig_list]
  
  # read in zero lag triggers
  if trigs:
    snr = trigs.getstat()
    if len(snr) > 0:
      internal_max = max(internal_max, snr.max())
      internal_min = min(internal_min, snr.min())
  
  # read in slide triggers
  if slide_trig_list:
    slide_snr_list = []
    for this_slide in slide_trig_list:
      slide_snr = this_slide.getstat()
      slide_snr_list.append(slide_snr)
      
      if len(slide_snr) > 0:
        internal_max = max(internal_max, slide_snr.max())
        internal_min = min(internal_min, slide_snr.min())

  # set up the bin boundaries
  if not max_val:
    max_val = internal_max
  if not min_val:
    min_val = internal_min
  
  if min_val >= max_val:
    # CHECKME: what should we do without any trigs or slide_trigs?
    # This is the old behavior.
    min_val = 5.
    max_val = 10.
  if min_val == max_val:
    # NB: this is numpy.histogram's default behavior for equal max and min
    min_val -= 0.5
    max_val += 0.5
    
  bins = numpy.linspace(min_val, max_val, nbins + 1, endpoint=True)

  # hist of the zero lag:
  if trigs:
    zero_dist, _ = numpy.histogram(snr, bins)
    cum_dist_zero = zero_dist[::-1].cumsum()[::-1]

  # hist of the slides:
  if slide_trig_list:
    cum_dist_slide = []
    for slide_snr in slide_snr_list:
      num_slide, _ = numpy.histogram(slide_snr, bins)
      cum_slide = num_slide[::-1].cumsum()[::-1]
      cum_dist_slide.append(cum_slide)
    cum_dist_slide = numpy.array(cum_dist_slide)
    slide_mean = cum_dist_slide.mean(axis=0)
    slide_std = cum_dist_slide.std(axis=0)
    if scalebkg:
      slide_mean *= 600./6370.
      slide_std *= sqrt(600./6370.)

  ds = (bins[1] - bins[0]) / 2
  if "bitten_l" in stat:
    lefts = bins[:-1]
    centers = bins[:-1] + ds
  else:
    lefts = bins[:-1]**2
    centers = (bins[:-1] + ds)**2

  figure()
  # plot zero lag
  if trigs:
    semilogy(centers,cum_dist_zero+0.0001,'r^',markerfacecolor="b",\
        markersize=12)

  # plot time slides
  if slide_trig_list and len(slide_snr_list):
    slide_min = []
    for i in range( len(slide_mean) ):
      slide_min.append( max(slide_mean[i] - slide_std[i], 0.0001) )
      slide_mean[i] = max(slide_mean[i], 0.0001)
    semilogy(centers,asarray(slide_mean), 'r+', markersize=12)
    tmpx, tmpy = makesteps(lefts, slide_min, slide_mean + slide_std)
    fill(tmpx, tmpy, facecolor='y', alpha=0.3)

  if stat == 'coherent_snr': xlab = 'Coherent SNR$^{2}$'
  elif stat: xlab = 'combined ' + stat.replace('_',' ')
  else: xlab = 'Combined Statistic'
  xlabel(xlab, size='x-large')
  ylabel('Number of events', size='x-large')
  title_text = 'Cum. hist. of num events vs ' + xlab
  if ifolist:
    title_text += ' for ' 
    for ifo in ifolist:
      title_text += ifo + ' '
 
  title(title_text, size='x-large')

######################################################################
# function to make a histogram of the coinc inspiral statistic
def histstat(trigs=None, slide_trigs=None,ifolist = None, min_val = None, \
  max_val = None, nbins = 20, stat=None):
  """
  function to plot a histogram of the snr of coincident events
  in the zero lag and times slides

  @param trigs: coincInspiralTable
  @param slide_trigs: dictionary of time slide triggers
  @param ifolist: list of ifos
  @param min_val: minimum of snr to be plotted
  @param max_val: maximum of snr to be plotted
  @param nbins: number of bins to use in histogram
  @param stat: the statistic being used
  """
  internal_min = numpy.inf
  internal_max = -numpy.inf
  
  # filter down to ifolist
  slide_trig_list = [s["coinc_trigs"] for s in slide_trigs]
  if ifolist:
    if trigs:
      trigs = trigs.coinctype(ifolist)
    slide_trig_list = [s.coinctype(ifolist) for s in slide_trig_list]
  
  # read in zero lag triggers
  if trigs:
    snr = trigs.getstat()
    if len(snr) > 0:
      internal_max = max(internal_max, snr.max())
      internal_min = min(internal_min, snr.min())
  
  # read in slide triggers
  if slide_trig_list:
    slide_snr_list = []
    for this_slide in slide_trig_list:
      slide_snr = this_slide.getstat()
      slide_snr_list.append(slide_snr)
      
      if len(slide_snr) > 0:
        internal_max = max(internal_max, slide_snr.max())
        internal_min = min(internal_min, slide_snr.min())

  # set up the bin boundaries
  if not max_val:
    max_val = internal_max
  if not min_val:
    min_val = internal_min
  
  if min_val >= max_val:
    # CHECKME: what should we do without any trigs or slide_trigs?
    # This is the old behavior.
    min_val = 5.
    max_val = 10.
  if min_val == max_val:
    # NB: this is numpy.histogram's default behavior for equal max and min
    min_val -= 0.5
    max_val += 0.5

  bins = numpy.linspace(min_val, max_val, nbins + 1, endpoint=True)

  # hist of the zero lag:
  if trigs:
    hist_zero, _ = numpy.histogram(snr, bins)

  # hist of the slides:
  if slide_trigs:
    slide_dist = []
    hist_slide = []
    for slide_snr in slide_snr_list:
      num_slide, _ = numpy.histogram(slide_snr, bins)
      hist_slide.append(num_slide)
    hist_slide = numpy.array(hist_slide)
    slide_mean = hist_slide.mean(axis=0)
    slide_std = hist_slide.std(axis=0)

  ds = (bins[1] - bins[0]) / 2
  if "bitten_l" in stat:
    lefts = bins[:-1]
    centers = bins[:-1] + ds
  else:
    lefts = bins[:-1]**2
    centers = (bins[:-1] + ds)**2

  figure()
  # plot zero lag
  if trigs and len(trigs):
    semilogy(centers, hist_zero + 0.0001, 'r^', markerfacecolor="b",\
        markersize=12)

  # plot time slides
  if slide_trigs and len(slide_snr_list):
    slide_min = []
    for i in range( len(slide_mean) ):
      slide_min.append( max(slide_mean[i] - slide_std[i], 0.0001) )
      slide_mean[i] = max(slide_mean[i], 0.0001)
    semilogy(centers, asarray(slide_mean), 'r+', markersize=12)
    tmpx, tmpy = makesteps(lefts, slide_min, slide_mean + slide_std)
    fill(tmpx, tmpy, facecolor='y', alpha=0.3)

  if stat == 'coherent_snr': xlab = 'Coherent SNR$^{2}$'
  else: xlab = 'Combined Statistic'
  xlabel(xlab, size='x-large')
  ylabel('Number of events', size='x-large')
  title_text = 'Histogram of Number of events vs ' + xlab
  if ifolist:
    title_text += ' for '
    for ifo in ifolist:
      title_text += ifo + ' '

  title(title_text, size='x-large')

######################################################################
# function to determine the efficiency as a function of distance
def efficiencyplot(found, missed, col_name, ifo=None, plot_type = 'linear', \
    nbins = 40, output_name = None, plotsym = 'k-', plot_name = '', \
    title_string = '', errors = False):
  """
  function to plot the difference if col_name_a in two tables against the
  value of col_name_b in table1.  
  
  @param found:  metaDataTable containing found injections
  @param missed: metaDataTable containing missed injections
  @param col_name: name of column used to plot efficiency
  @param ifo: name of ifo (default = None), 
              used in extracting information (e.g. which eff_dist)
  @param plot_type: either 'linear' or 'log' plot on x-axis
  @param plotsym:  the symbol to use when plotting, default = 'k-'
  @param plot_name: name of the plot (for the legend)
  @param title_string: extra info for title
  @param errors: plot errorbars on the efficiencies (using binomial errors)
                 default = False
  """
  
  if not ifo and ("ifo" in found.validcolumns.keys()):
    ifo = found[0].ifo

  foundVal = readcol(found,col_name, ifo)
  missedVal = readcol(missed,col_name, ifo)

  if len(foundVal) or len(missedVal):
    # we have found or missed injections so we can generate the plot
    if plot_type == 'log':
      foundVal = log10(foundVal)
      missedVal = log10(missedVal)

    if len(foundVal):
      bins = numpy.linspace(min(foundVal), max(foundVal), nbins + 1,
                            endpoint=True)
      step = bins[1] - bins[0]
      plotbins = bins[0:-1] + step/2.
      if step == 0:
        bins = array([foundVal[0]/2.0, foundVal[0], foundVal[0] * 3.0/2.0])
        plotbins = bins[0:-1] + foundVal[0]/4.0
    else:
      bins = numpy.linspace(min(missedVal), max(missedVal), nbins + 1,
                            endpoint=True)
      step = bins[1] - bins[0]
      plotbins = bins[0:-1] + step/2.
      if step == 0:
        bins = array([missedVal[0]/2.0, missedVal[0], missedVal[0] * 3.0/2.0])
        plotbins = bins[0:-1] + missedVal[0]/4.0
    num_found, _ = numpy.histogram(foundVal, bins)
    num_missed, _ = numpy.histogram(missedVal, bins)

    fig_num = gcf().number
    figure(fig_num)
    num_found = array(num_found,'d')
    eff = num_found / (num_found + num_missed)
    total_num = num_found + num_missed
    yerr_common = total_num * (2 * num_found + 1)
    yerr_denom = 2*total_num*(total_num + 1)
    yerr_vary = 4 * total_num * num_found * (total_num - num_found) +\
                total_num**2
    yerr_vary = yerr_vary**0.5
    yerr_low = (yerr_common - yerr_vary)/yerr_denom
    yerr_high = (yerr_common + yerr_vary)/yerr_denom
   
    error = [(num_found/total_num) - yerr_low,yerr_high - (num_found/total_num)]

    if plot_type == 'log':
      plotbins = 10**plotbins
      if errors:
        errorbar(plotbins, eff, yerr=error, markersize=12, markerfacecolor='None',\
            markeredgewidth=1, linewidth=2, label = plot_name, \
            fmt = plotsym)
        semilogx()
      else:
        semilogx(plotbins, eff, plotsym, markersize=12, markerfacecolor='None',\
            markeredgewidth=1, linewidth=2, label = plot_name)
              
    else:
      if errors:
        errorbar(plotbins, eff, yerr=error, markersize=12, markerfacecolor='None',\
            markeredgewidth=1, linewidth=2, label = plot_name, \
            fmt = plotsym)
      else:
        plot(plotbins, eff, plotsym, markersize=12, markerfacecolor='None',\
            markeredgewidth=1, linewidth=2, label = plot_name)

    xlabel(col_name.replace("_"," "), size='x-large')
    ylabel('Efficiency', size='x-large')
    ylim(0,1.1)
  else:
    # no found or missed injections
    figtext(0,0,'No found or missed injections',fontsize=32)

  if ifo:
    title_string += ' ' + ifo  
  
  title_string += ' ' + col_name.replace("_"," ")
  title_string += ' efficiency plot'
  title(title_string, size='x-large')

  grid(True)

  if output_name:
    if ifo:
      output_name += '_' + ifo
    output_name += '_' + col_name + '_eff.png'
    savefig(output_name)


######################################################################
# function to histogram the difference between values of 'col_name' in
# two tables, table1 and table2
def histdiff(table1, table2, col_name, plot_type, hist_num, 
  total_hists=1,  nbins=20, hist_width=[0,0], hist_norm=None):
  """
  function to plot a histogram of the difference of the value of col_name
  between table1 and table2  
  
  @param table1: metaDataTable
  @param table2: metaDataTable
  @param col_name: name of column to plot
  @param plot_type: either 'hist' or 'frac_hist' 
  @param hist_num: the number of the histogram (used for color + placement)
  @param total_hists: total number of histograms (used to place bars)
  @param nbins: number of bins to plot in histogram (default = 20)
  @param hist_width:the minimum/maximum difference to be shown (2 element list)
  @param hist_norm: normalization of the histogram (total nevents)
  """

  histcolors = ['b','r','g','k']

  [tmpvar1, tmpvar2, ifo ] = readcolfrom2tables(table1, table2, col_name)
  tmp_diff = tmpvar2 - tmpvar1

  if (plot_type == 'frac_hist'):
    tmp_diff /= tmpvar1

  if hist_width[0] and hist_width[1]:
    bins = numpy.linspace(hist_width[0], hist_width[1], nbins + 1,
                          endpoint=True)
    height, _ = numpy.histogram(tmp_diff, bins=bins)
  else:
    height, bins = numpy.histogram(tmp_diff, bins=nbins)
  bins = bins[:-1]

  fig_num = gcf().number
  figure(fig_num)
  
  width = (bins[1] - bins[0]) / total_hists
  
  if hist_norm:
    height = height / float(hist_norm)
  left = []
  for val in bins:
    val = val + (width * hist_num)/2
    left.append(val)
 
  bar(left,height,width,color=histcolors[hist_num])

  # figtext(0.13,0.8 - 0.1* hist_num," mean = %6.3e" % mean(tmp_diff))
  # figtext(0.13,0.75 - 0.1 * hist_num,'sigma = %6.3e' % std(tmp_diff))
  
######################################################################
# function to histogram the difference between values of 'col_name' in
# two tables, table1 and table2
def labelhistdiff(col_name, plot_type, units, leg = None, title_text = None, 
  output_name = None):
  """
  function to label the output of histdiff

  @param col_name: name of column to plot
  @param plot_type: either 'hist' or 'frac_hist' 
  @param units: the units of the column
  @param leg: legend to add to plot
  @param title_text: text to add at start of title
  @param output_name: used in naming output file
  """


  label = col_name 
  if (plot_type == 'frac_hist'):
    label += ' fractional'
  label += ' difference'
  if units and not (plot_type == 'frac_hist'):
    label += ' (' + units +')'
  xlabel(label.replace("_"," "), size='x-large')

  ylabel('Number', size='x-large')
  
  if title_text:
    title(title_text + ' ' + col_name.replace("_"," ") + '  Histogram', \
        size='x-large')
  else:
    title(col_name.replace("_"," ") + ' Histogram', size='x-large')
  
  grid(True)

  if output_name:
    if (plot_type == 'frac_hist'):  
      output_name += '_' + col_name + '_frac_histogram.png'
    else:
      output_name += '_' + col_name + '_histogram.png'
    savefig(output_name)


###################################################
# function to plot the difference between values of 'col_name' in
# two tables, table1 and table2
def histdiffdiff(ifo1_trig, ifo2_trig, inj, col_name, sym, units=None, 
  nbins = None, width = None, output_name = None):
 
  histcolors = ['b','r','k']

  [tmpvar1, injvar1, ifo1 ] = readcolfrom2tables(ifo1_trig, inj, col_name)
  [tmpvar2, injvar2, ifo2 ] = readcolfrom2tables(ifo2_trig, inj, col_name)
  
  diff1 = injvar1 - tmpvar1
  diff2 = injvar2 - tmpvar2
  
  diffdiff = diff1 - diff2
  
  if not nbins:
    nbins = 10
  
  bins = []
  if width and width[0]!=width[1]:
    for i in range(-nbins,nbins):
      bins.append(width * i/nbins)

  if bins:
    out = hist(diffdiff,bins)
  else:
    out = hist(diffdiff,nbins)

  width = out[1][1] - out[1][0]
  bar(out[1],out[0],width,color=histcolors[sym])

  diffdiff = asarray(diffdiff)
  figtext(0.13,0.8 - 0.1* sym," mean = %6.3e" % mean(diffdiff))
  figtext(0.13,0.75 - 0.1 * sym,'sigma = %6.3e' % std(diffdiff))
 
  label = col_name 
  label += ' difference'
  if units:
    label += ' (' + units +')'
  xlabel(label.replace("_"," "), size='x-large')

  ylabel('Number', size='x-large')
  
  title(ifo1 + ' - ' + ifo2 + ' ' + col_name.replace("_"," ") + \
      ' difference histogram', size='x-large')
  
  grid(True)

  if output_name:
    output_name += '_' + ifo1 + '_' + ifo2 
    output_name += '_' + col_name + '_histogram.png'
    savefig(output_name)
  

######################################################################
# function to histogram the difference between values of 'col_name' in
# two tables, table1 and table2
def histslides(slide_trigs, zerolag_trigs = None, ifolist = None, scalebkg = None):
  """
  function to make a histogram of the number of triggers per time slide
  
  @param slide_trigs: dictionary of time slide triggers
  @param zerolag_trigs: coincInspiralTable
  @param ifolist: list of ifos
  """
  nevents = []
  slides = []
  for slide in slide_trigs:
    if ifolist:
      nevents.append( len(slide["coinc_trigs"].coinctype(ifolist)) )
    else:  
      nevents.append( len(slide["coinc_trigs"]) )
      
    slides.append(slide["slide_num"])

  mean_events = mean(nevents)
  std_events = std(nevents)
  if scalebkg:
    for i in range(len(nevents)):
      nevents[i] = nevents[i] * ( 600.0 / 6370.0 )
    mean_events = mean(nevents)
    std_events = std(nevents) * ( 6370.0 / 600.0 ) ** ( 0.5 )

  hist(nevents, align='center')

  figtext(0.13,0.8, " mean = %6.3e" % mean_events)
  figtext(0.13,0.75,"sigma = %6.3e" % std_events)
  if zerolag_trigs:
    hold(True)
    if ifolist:
      nfgevents = len(zerolag_trigs.coinctype(ifolist))
    else:
      nfgevents = len(zerolag_trigs)
    figtext(0.13,0.70,"zero lag = %6.3e" % nfgevents )
    axvline(nfgevents,color='r',linewidth=2)

  v = axis()
  
  if scalebkg and mean_events>0:
    x_gauss = []
    y_gauss = []
    min_events = min(nevents)
    max_events = max(nevents)
    bin_size = ( max_events - min_events ) / 10.0
    normalization = bin_size * len(nevents) / \
      ( std_events * ( 2.0 * pi ) ** ( 0.5 ) )
    for i in range(101):
      x_gauss.append( mean_events + 3.0 * ( ( i / 50.0 ) - 1.0 ) * std_events )
      y_gauss.append( exp( -( ( x_gauss[i] - mean_events ) ** ( 2.0 ) ) / \
        ( 2.0 * ( std_events ** ( 2.0 ) ) ) ) * normalization )
    plot(x_gauss,y_gauss,'b',linewidth=2)
    axis(v)
  
  xlabel('Number of triggers',size='x-large')
  title_text = 'Histogram of number coincident '
  if ifolist:
    for ifo in ifolist:
      title_text += ifo + ' '
  title_text += 'triggers per time slide'
  title(title_text, size='x-large')
  
  
######################################################################
# function to histogram the difference between values of 'col_name' in
# two tables, table1 and table2
def plotslides(slide_trigs, zerolag_trigs = None, ifolist = None, \
    scalebkg = None):
  """
  function to make a histogram of the number of triggers per time slide
  
  @param slide_trigs: dictionary of time slide triggers
  @param zerolag_trigs: coincInspiralTable
  @param ifolist: list of ifos
  @param scalebkg: if True then scale background by 600/6370
  """
  nevents = []
  slides = []
  for slide in slide_trigs:
    if ifolist:
      nevents.append( len(slide["coinc_trigs"].coinctype(ifolist)) )
    else:  
      nevents.append(len(slide["coinc_trigs"]))
    slides.append(slide["slide_num"])

  mean_events = mean(nevents)
  std_events = sqrt(mean_events)

  if scalebkg:
    for i in range(len(nevents)):
      nevents[i] = nevents[i] * ( 600.0 / 6370.0 )
    mean_events = mean(nevents)
    std_events = std(nevents) * ( 6370.0 / 600.0 ) ** ( 0.5 )

  bar(slides, nevents, 0.8, 0, color='b', align="center") 
  axhline(mean_events,color='k',linewidth=2)
  axhline(mean_events + std_events,color='k',linestyle='--',linewidth=2)
  axhline(max(mean_events - std_events,0),color='k',linestyle='--',linewidth=2)
 
  if zerolag_trigs:
    hold(True)
    if ifolist:
      nfgevents = len(zerolag_trigs.coinctype(ifolist))
    else:
      nfgevents = len(zerolag_trigs)
    bar([0], [nfgevents], 0.8, 0, color='r', align="center") 
 
  xlim(min(slides)-0.5,max(slides)+0.5) 
  xlabel('Number of time slide',size='x-large')
  ylabel('Number of triggers in slide',size='x-large')
  title_text = 'Plot of number coincident '
  if ifolist:
    for ifo in ifolist:
      title_text += ifo + ' '
  title_text += 'triggers per time slide'
  title(title_text, size='x-large')

   
######################################################################
def tfplot(*args, **kwargs):
  """
  tfplot(x, y, s=20, c='b', marker='o', cmap=None, norm=None,
    vmin=None, vmax=None, alpha=1.0)

  Supported function signatures:

  TFPLOT(x, y)  : make a scatter plot of x vs y

  TFPLOT(x, y, s)  : make a scatter plot of x vs y with size in area
  given by s

  TFPLOT(x, y, s, c) : make a scatter plot of x vs y with size in area
  given by s and colors given by c

  TFPLOT(x, y, s, c, **kwargs) : control colormapping and scaling
  with keyword args; see below

  Make a scatter plot of x versus y.  s is a size in points^2 a scalar
  or an array of the same length as x or y.  c is a color and can be a
  """
  shading = kwargs.get('shading', 'faceted')
  cmap = kwargs.get('cmap', cm.get_cmap())
  norm = kwargs.get('norm', normalize())
  alpha = kwargs.get('alpha', 1.0)
  vmin = kwargs.get('vmin', None)
  vmax = kwargs.get('vmax', None)  
  a = kwargs.get('axes', gca())

  if len(args)==5:
      X, dX, Y, dY, C = args
  else:
      raise TypeError, 'Illegal arguments to rectfill; see help(rectfill)'
  
  Nx, = X.shape
  verts = [ ( (X[i,] , Y[i,]) , (X[i,]+dX[i,] , Y[i,]),
              (X[i,]+dX[i,] , Y[i,]+dY[i,]), 
              (X[i,] , Y[i,]+dY[i,]) )
            for i in range(Nx-1) ] 
  C = array([C[i,] for i in range(Nx-1)])
              
  if shading == 'faceted': edgecolors =  (0,0,0,1), 
  else: edgecolors = 'None'
  
  collection = PolyCollection(
          verts,
          edgecolors   = edgecolors,
          antialiaseds = (0,),
          linewidths   = (0.25,),
          )
  collection.set_alpha(alpha)
  collection.set_array(C)
  if norm is not None: assert(isinstance(norm, normalize))
  if cmap is not None: assert(isinstance(cmap, Colormap))
  collection.set_cmap(cmap)
  collection.set_norm(norm)
  if norm is not None: collection.set_clim(vmin, vmax)
  minx = amin(X)
  maxx = amax(X)
  miny = amin(Y)
  maxy = amax(Y)
  corners = (minx, miny), (maxx, maxy)      
  a.update_datalim( corners )
  a.autoscale_view()
  # add the collection last
  a.add_collection(collection)
  return collection

######################################################################
def rescale_axis(limNew, limOld):    
  """
  Rescaling axes, i.e. converting one axis to the other.
  This is useful to make a contour plot with the correct
  axis labels.
  @param limOld: Two-component vector containing the new limits
  @param limNew: Two component vector containing the original limits
  
  Example:
  
  If you plot a 20x20 matrix as a contour plot the axis limits
  go from 0 to 20. But the x-axis really spans values from e.g.
  3 to 10. Then, to get the tick positions and the tick labels
  to correct for this, you do:
  
  tickPos, tickLabels = rescale_axis([0, 20], [3, 10])
  tickPos = [1.4799985170355541, 5.1799948096244393, \
             8.8799911022133244, 12.57998739480221, \
             16.279983687391095, 19.97997997997998]
  tickLimits = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
  
  tickPos contains the position at which you have to put the labels
  in the vector tickLimits    
  """

  # get the difference of the new labels
  Delta=limNew[1]-limNew[0]

  # calculate a first order step width for the tick labels        
  a = round( log10(3*Delta/20.0) )
  step = pow(10, a)

  # fine-tuning of the step-width on the new axis,
  # so there are neither too many tick nor too few
  nc = Delta/step
  if nc<4:
    step/=2
  if nc>15:
    step*=2

  # get the starting and end-points of the axes
  # adding one to cover the entire senseful range
  startTick = step*int(limNew[0]/step+1.0)
  endTick = step*int(limNew[1]/step+1.0)

  # calculate the linear transformation parameters
  slope = Delta/(limOld[1]-limOld[0])
  intercept = limNew[0]-slope*limOld[0]
        
  # compute the axis positions and labels
  tickPos = []
  tickLabel = []
  for t in arange( startTick, endTick, step):
    tickLabel.append(t)
    tickPos.append((t-intercept)/slope)
    
  # returning the results
  return tickPos, tickLabel, intercept, slope


###############################
def create_contour( datax, datay ):

    # create the array
    xbins = rate.LinearBins(min(datax), max(datax), 20)
    ybins = rate.LinearBins(min(datay), max(datay), 20)    
    h = rate.BinnedArray(rate.NDBins((ybins,xbins)))

    # fill this array
    for x, y in zip(datax, datay):
        h[y,x]+=1

    mat = h.array
    defX = xbins.centres()
    defY = ybins.centres()

    # get the dimensions
    n=len(numpy.unique(defX))
    m=len(numpy.unique(defY))

    # get the boundaries
    minX=min(defX)
    maxX=max(defX)  
    minY=min(defY)
    maxY=max(defY)

    # create the plot
    contourf(mat,40, alpha=0.3)

    # set the axes
    ax = gca()
    [xtick, xlabel, interX, slopeX]=rescale_axis([minX, maxX], [0,n])
    ax.set_xticks(xtick)
    ax.set_xticklabels(xlabel)
    [ytick, ylabel, interY, slopeY]=rescale_axis([minY, maxY],[0, m])
    ax.set_yticks(ytick)
    ax.set_yticklabels(ylabel)

   

######################################################################
def test_rescale_axis():
  """
  Unit test(?) for rescale_axis
  """
  
  # define input values
  input = []
  input.append([[0,20],[0,20]])
  input.append([[0,200],[0,20]])
  input.append([[-45.5, 300], [0, 20]])
  
  # define output values
  output = []
  output.append([[2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0],\
                 [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]])
  output.append([[2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0],\
                 [20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0]])
  output.append([[2.6338639652677283, 5.5282199710564406, 8.4225759768451525,\
                  11.316931982633864, 14.211287988422578, 17.105643994211288, 20.0],\
                 [0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0]])
  
  # accuracy of this test
  eps = 1.0e-5
  
  # do the test
  test=[None,None]
  for inp, out in zip(input, output):

    # test the unit
    test[0], test[1] = rescale_axis(inp[0], inp[1])

    for k in range(2):
      if len(test[k])!=len(out[k]):
        print "Length of output vectors inequal"
        return 0
      for x, y in zip( out[k], test[k]):
        if abs(x-y)>eps:
          print "Test variables inequal. Expected %f got %f " %\
                (y,x)
          return 0

  # module succesfull tested
  return 1
      
######################################################################
def main():
  # define usage and command line options and arguments - parse
  usage = "usage: %prog ..."
  parser = OptionParser( usage )

  opts_snglInsp = OptionGroup( parser, "Single Inspiral plotting functions",\
	"Example ..." )
  opts_snglInsp.add_option( "-a", "--snglInsp_snrVtime",\
	action="store_true", default=False,\
	help="plot snr vs time from a single inspiral xml" )
  opts_snglInsp.add_option( "-b", "--snglInsp_snrVchisq",\
        action="store_true", default=False,\
        help="plot snr vs chisq from single inspiral xml")
  opts_snglInsp.add_option( "-c", "--snglInsp_histHistc_snr",\
        action="store_true", default=False,\
        help="plot snr histograms from single inspiral xml" )
  opts_snglInsp.add_option( "-d", "--snglInsp_summary",\
        action="store_true", default=False,\
        help="plot summary info from single inspiral xml" )

  parser.add_option_group( opts_snglInsp )

  parser.add_option( "-p", "--show_plot",\
        action="store_true", default=False,\
        help="display plot" )
  # change this so that by default the fig is always saved using 
  # the name convention already implemented. Now instead of --save-fig
  # you have --save-off and --save-as flag to override
  # the standard name. Also add the date+time to the standard name OR
  # check for the existence of the standard name + 001, 002, 003, ...
  # Remove where possible the defns of dest in favour of the long option name
  parser.add_option( "-s", "--save_fig",\
        action="store_true", default=False,\
        help="save figure in .png and .ps format" )
  parser.add_option( "-t", "--temporary-test",\
        action="store_true", default=False,\
        help="only for developers to test this program" )

  ( options, xml_files ) = parser.parse_args()
  
  # check xml_files have been given as required arguments 
  if not xml_files:
    print >> sys.stderr, "No trigger file specified"
    sys.exit(1)

  # read data files and call plotting function desired
  if   options.snglInsp_snrVtime:
    trigs = snglInspiral(xml_files[0])
    trigs.plot_snr_v_time()
  elif options.snglInsp_snrVchisq:
    trigs = snglInspiral(xml_files[0])
    trigs.plot_snr_v_chisq()
  elif options.snglInsp_histHistc_snr:
    trigs = snglInspiral(xml_files[0])
    trigs.histHistc_snr()
  elif options.snglInsp_summary:
    trigs = snglInspiral(xml_files[0])
    trigs.summary()
  else:
    print >> sys.stderr, "No plot option specified"
    sys.exit(1)
  
  # save and show plot if desired
  if options.save_fig:
    png_file = xml_file[:-3] + plot_type + ".png"
    ps_file  = xml_file[:-3] + plot_type + ".ps"
    savefig(png_file)
    savefig(ps_file)
    print "Saved plot to file %s" % (png_file)
    print "Saved plot to file %s" % (ps_file)
  if options.show_plot:
    show()

# execute main if this module is explicitly invoked by the user
if __name__=="__main__":
        main()
