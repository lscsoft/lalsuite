#!/usr/bin/env python

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

__prog__ = 'lalapps_cbc_plotrates'
description = 'Creates plots of histogram counts and cumulative rates vs. a given statistic. Statistics can be either single or coincident, and can be queried or computed on the fly from either a sngl table or a coinc table. NOTE: This is a modified version of that in lalsuite. This also plots extended background triggers, which must be stored in text files. It also writes out text files of the data that was plotted in the cumulative histogram.'
usage = '%s [options] file1.sqlite file2.sqlite ...'

from optparse import OptionParser
try:
    import sqlite3
except ImportError:
    # pre 2.5.x
    from pysqlite2 import dbapi2 as sqlite3
import sys
import os
import re
import math
import numpy
from scipy import special
import bisect
import matplotlib
matplotlib.use('Agg')
import pylab
pylab.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
pylab.rc('text', usetex=True)

from glue import lal
from glue import segments
from glue.ligolw import dbtables
from glue.ligolw import lsctables
from glue import git_version

from pylal import ligolw_sqlutils as sqlutils
from pylal import CoincInspiralUtils
from pylal import InspiralUtils
from pylal import printutils
from pylal.ligolw_dataUtils import get_row_stat

# =============================================================================
#
#                                   Set Options
#
# =============================================================================


def parse_command_line(argv = None):
    """
    Parse the command line, return options and check for consistency among the
    options.
    """
    parser = OptionParser( version = git_version.verbose_msg, usage = usage, description = description )

    # following are related to file input and output naming
    parser.add_option( "-o", "--output-path", action = "store", type = "string", default = './',
        help =
            "Output path to save figures to."
            )
    parser.add_option( "-u", "--user-tag", action = "store", type = "string", default = '',
        help =
            "User-tag for plot file names"
            )
    parser.add_option( "-t", "--tmp-space", action = "store", type = "string", default = None,
        metavar = "PATH",
        help =
            "Location of local disk on which to do work. " +
            "This is used to enhance performance in a networked " +
            "environment."
            )
    parser.add_option( "-c", "--cache-file", action = "append", default = [],
        help =
            "Read files from the given cache file. File must be lal formated. To specify multiple, give the argument multiple times."
            )
    parser.add_option( "-v", "--verbose", action = "store_true", default = False,
        help =
            "Print progress information"
           )
    parser.add_option( "", "--debug", action = "store_true", default = False,
        help =
            "Print SQLite queries used to obtain data."
           )
    # following are specific to this program
    parser.add_option("-s", "--statistic", action = "store", type="string",
        default=None, metavar="STATTYPE(:label)",
        help =
            "Stat to plot. Can be any column or any function of columns in the given table. Syntax for functions is python; functions available are anything from python's math module. To specify a label for the stat when plotting, append a colon and the label."
        )
    parser.add_option("-d", "--data-table", action = "store", type = "string", default = None, 
        help =
            "Table to get triggers from. Can be any table with either a coinc_event_id or an event_id. If the table has a coinc_event_id, it will be assumed to be a coincident table and so must have an 'ifos' column. If the table has an event_id, it will be assumed to be a single event table, and so must have an 'ifo' column. Additionally, if a combining-method or function is specified, single events will be groupped together via their event_ids in the coinc_event_map table. Note: Plotting using single event tables can cause the program to take siginificantly longer to execute than coincident event tables. If a coincident statistic has already been computed and exists in the coinc_event table, always use that instead."
        )
    parser.add_option("-C", "--combining-method", action = "store", type="string", metavar = 'sum | quad_sum | mean | median | min | max', default=None,
        help =
            "How to combine the given statistic for each ifo. This is only needed if trying to plot coincident statistics from a single event table. Supported methods are sum, quad_sum (the square root of the sum of the squares), mean, median, min, or max values."
        )
    parser.add_option("-F", "--combining-function", action = "store", type="string", default=None,
        help =
            "A more explicit alternative to combining-method, this allows you to specify exactly how to combine values from each ifo. The expression should be a function of the ifos to be combined; if an ifo isn't present in a given coinc, it will be given a 0 value. For instance, if --statistic is set to snr, and the combining-method is set to sqrt(H1**2+L1**2+(.5*V1)**2), then the sum of the squares of H1, L1, and V1 will be plotted, with a weight of .5 given to V1. If doubles are present, they will also be plotted. Only ifos that are specified will be included in coincs. In other words, if only H1 and L1 are listed in the equation only H1 and L1 will be included, even if H2 was listed as being coincident with them in the database. If combining-method is also specified, combining-function takes precedence."
        )
#    parser.add_option("-S", "--reference-snr", action = "store", type = "float",
#        default = 8.,
#        help =
#            "Reference snr to plot. Default is 8.0."
#        )
    parser.add_option("-f", "--foreground-datatype", action = "store", type = "string", default = "slide",
        help =
            "Plot the given foreground type. Options are all_data, exclude_play, playground, or simulation. Default is to plot no foreground. Background expected rates are determined by the livetime of the specified datatype. If not foreground datatype specified, expected rates are based on the average background time."
       )
    parser.add_option( "-p", "--param-name", metavar = "PARAMETER[:label]",
        action = "store", default = None,
        help =
            "Specifying this and param-ranges will only select triggers that fall within the given range for each plot. " +
            "This will apply to all plots. Any column(s) in the simulation or recovery table may be used. " +
            "To specify a parameter in the simulation table, prefix with 'injected_'; " +
            "for recovery table, 'recovered_'. As with variables, math operations are permitted between columns, e.g.," +
            "injected_mass1+injected_mass2. Any function in the python math module may be used; syntax is python. The parameter name " +
            "will be in the title of each plot; to give a label for the name, add a colon and the label. If no label is specified, the " +
            "the parameter name as given will be used. " +
            "WARNING: if you use any recovered parameters, missed injections will not be included in plots " +
            "as recovered parameters are undefined for them. For example, if you chose to bin by recovered_mchirp " +
            "missed injections will not be plotted on any plot. If you chose to bin by injected_mchirp, however, missed injections " +
            "will be plotted. "
        )
    parser.add_option( "-q", "--param-ranges", action = "store", default = None,
        metavar = " [ LOW1, HIGH1 ); ( LOW2, HIGH2]; !VAL3; etc.",
        help = 
            "Requires --param-name. Specify the parameter ranges " +
            "to select triggers in. A '(' or ')' implies an open " +
            "boundary, a '[' or ']' a closed boundary. To specify " +
            "multiple ranges, separate each range by a ';'. To " +
            "specify a single value, just type that value with no " +
            "parentheses or brackets. To specify not equal to a single " +
            "value, put a '!' before the value. If " +
            "multiple ranges are specified, a separate plot for each range will be generated."
        )
    parser.add_option("", "--param-combining-method", default = "mean",
        help =
            "For sngl-ifo tables, how to combine the --param-name. Either a combining-function or method may be used. Default is mean."
        )
    parser.add_option( "", "--exclude-coincs", action = "store", type = "string", default = None,
        metavar = " [COINC_INSTRUMENTS1 + COINC_INSTRUMENTS2 in INSTRUMENTS_ON1];"
            "[ALL in INSTRUMENTS_ON2]; etc.",
        help = 
            "Exclude coincident types in specified detector times, " +
            "e.g., '[H2,L1 in H1,H2,L1]'. Some rules: " +
                "* Coinc-types and detector time must be separated by " +
                "an ' in '. When specifying a coinc_type or detector " +
                "time, detectors and/or ifos must be separated by " +
                "commas, e.g. 'H1,L1' not 'H1L1'. " +
                "* To specify multiple coinc-types in one type of time, " +
                "separate each coinc-type by a '+', e.g., " +
                "'[H1,H2 + H2,L1 in H1,H2,L1]'. " +
                "* To exclude all coincs in a specified detector time " +
                "or specific coinc-type in all times, use 'ALL'. E.g., " +
                "to exclude all H1,H2 triggers, use '[H1,H2 in ALL]' " +
                "or to exclude all H2,L1 time use '[ALL in H2,L1]'. " + 
                "* To specify multiple exclusions, separate each " +
                "bracket by a ';'. " +
                "* Order of the instruments nor case of the letters " +
                "matter. So if your pinky is broken and you're " +
                "dyslexic you can type '[h2,h1 in all]' without a " +
                "problem."
            )
    parser.add_option( "", "--include-only-coincs", action = "store", type = "string", default = None,
        metavar = " [COINC_INSTRUMENTS1 + COINC_INSTRUMENTS2 in INSTRUMENTS_ON1];" +
            "[ALL in INSTRUMENTS_ON2]; etc.",
        help =
            "Opposite of --exclude-coincs: only plot the specified coinc types. WARNING: If you type 'ALL' for on-instruments, e.g., '[H1,L1 in ALL]', and plot-by-instrument-time is not on, livetime from *all* of the instrument times will be included in the background, even if it is from an instrument time that could not produce H1,L1 triggers. If this isn't desired, either be explicit about what instrument times you really want -- i.e., don't use 'ALL' -- or use --exclude-coincs to exclude the instrument times you don't want. In the example, adding --exclude-coincs [ALL in L1,V1];[ALL in H1,V1] will prevent any H1,V1 or L1,V1 time from being added."
            )
    parser.add_option( "", "--plot-special-time", action = "append", type = "float", default = [],
        help =
            "Exclude all single triggers that are within N seconds of the given end-time from the background, where N is set by --plot-special-window. All foreground triggers that fall within the range will be plotted as stars. Note: this can only be done when reading from sngl-ifo tables."
            )
    parser.add_option( "", "--plot-special-window", action = "store", type = "float", default = 8.0,
        help =
            "Specify the time window, in seconds, to use around the plot-special-time. Default is 8.0 seconds."
            )
    parser.add_option( "", "--add-background-file", action = "append", default = [], metavar = " FILE,BACKGROUND_LIVETIME(:LABEL)",
        help =
            "Add a point to the background. Any point added will be plotted, but will not be added to the background distribution obtained from input files for cumulative rate computation nor for extrapolation. If the point added is louder than the loudest background trigger, this point will be used for extending the probability contours (if plotted). Points should be specified by first providing the x-value of the point, then the total background livetime used to acquire that point, in the same units as --time-units. To specify multiple points, give the argument multiple times."
            )
    parser.add_option( "", "--add-special-background-file", action = "append", default = [], metavar = " FILE,BACKGROUND_LIVETIME(:LABEL)",
        help =
            "Add a point to the background. Any point added will be plotted, but will not be added to the background distribution obtained from input files for cumulative rate computation nor for extrapolation. If the point added is louder than the loudest background trigger, this point will be used for extending the probability contours (if plotted). Points should be specified by first providing the x-value of the point, then the total background livetime used to acquire that point, in the same units as --time-units. To specify multiple points, give the argument multiple times."
            )
    parser.add_option( "", "--single-table", action = "store", default = None,
        help =
            "If using plot-special-time with a coinc table, a --single-table must be specified to look in for the un-slid single event end-times."
            )
    parser.add_option( "-I", "--plot-by-ifos", action = "store_true", default = False,
        help =
            "Create a separate plot for each coincidence type (if plotting a combined statistic) or single ifo (if plotting single ifo statistics). Default is to plot all ifo(s) together."
            )
    parser.add_option( "-T", "--plot-by-instrument-time", action = "store_true", default = False,
        help =
            "Create a separate plot for each instrument time. Default is to plot all instrument times together."
            )
    parser.add_option( "-M", "--plot-by-multiplicity", action = "store_true", default = False,
        help =
            "Create a separate plot based on the number of coincident instruments. (For single ifo plotting, this is a no-op.) In other words, doubles will be plotted with doubles, triples with triples, etc. Default is to plot all coincident ifos together. If this and plot-by-ifos specified, plot-by-ifos takes precedence."
            )
    parser.add_option( "-U", "--time-units", type = "string", default = "yr",
        help =
            "What unit of time to plot the rates in. Can be either 'yr', 'days', 'hr', 'min', or 's'. Default is yr."
            )
    parser.add_option( "-R", "--use-ifo-colors", action = "store_true", default = False,
        help =
            "Color foreground points by ifo color. Default is to plot no color (all triggers have a blue outline)."
            )
    parser.add_option( "-H", "--plot-poisson-pdf", action = "store_true", default = False,
        help =
            "Plot the Poisson probability density of the background."
            )
    parser.add_option( "-K", "--plot-significance-bands", action = "store_true", default = False,
        help =
            "Plot the Poisson probability bands in terms of the number of sigma. The values of the n sigma are the points in the poisson PDF that equal n*erf(y)."
            )
    parser.add_option( "-S", "--max-sigma", type = 'int', default = 5,
        help =
            "Maximum number of significance bands to plot. Default=5"
            )
    parser.add_option( "-E", "--extrapolate", default = None, metavar = "Gaussian|Power|erfc",
        help =
            "Extrapolate the background using a the given model out to the largest foreground statistic. If this is beyond the measured background, the probability density beyond the loudest background point will be calculated using this extrapolation. Options for models are 'Gaussian', 'Power', or 'erfc'. If 'Gaussian' or 'Power' a least-squares fit is done to both the non-cumulative and the cumulative histograms: for 'Gaussian' y = A*exp(-B*x**(2./x-power)) is fit; if 'Power', y = A*x**B is fit. If 'erfc', a Gaussian is fitted to the non-cumulative histogram to get sigmasq and the scalling factor A. These values are then used in the complitmentary error function y = A*\int{dx*exp(-x^2/2*sigmasq)} to extrapolate the cumulative distribution."
            )
    parser.add_option( "-P", "--x-power", type = "float", default = 1.,
        help =
            "Power of the combined statistic to use when fitting a gaussian to the data and binning for the non-cumulative plot.  The natural log of the Gaussian is taken so that the actual fitted equation is log(y) = Beta[1] + Beta[2]*x**(2/x-power). For example, if plotting linearly in your statistic, x-power should be 1. If plotting against the statistic squared, then x-power would be 2, thus the fit would be log(y) = Beta[1] + Beta[2]*x. Default=1."
            )
    parser.add_option( "", "--min-x-fit", type = 'float', default = None,
        help =
            "Background stat at which to begin fitting to for extrapolation. Default is the first background point."
            )
    parser.add_option( "", "--lin-x", action = "store_true", default = False,
        help =
            "Plot x-axis on a linear scale. Default is to plot on a log scale."
            )
    parser.add_option( "", "--lin-y", action = "store_true", default = False,
        help =
            "Plot y-axis on a linear scale. Default is to plot on a log scale."
            )
    parser.add_option('-a', '--xmin', action = 'store', type = 'float', default = None,
        help =
            'Set a minimum value for the x-axis. If lin-x not set, must be > 0. This will apply to all plots generated.'
        )
    parser.add_option('-b', '--xmax', action = 'store', type = 'float', default = None,
        help =
            'Set a maximum value for the x-axis. If lin-x not set, must be > 0. This will apply to all plots generated.'
        )
    parser.add_option('-A', '--ymin', action = 'store', type = 'float', default = None,
        help =
            'Set a minimum value for the y-axis. If lin-y not set, must be > 0. This will apply to all plots generated.'
        )
    parser.add_option('-B', '--ymax', action = 'store', type = 'float', default = None,
        help =
            'Set a maximum value for the y-axis. If lin-y not set, must be > 0. This will apply to all plots generated.'
        )
    parser.add_option( '', '--nbins', type = 'int', default = 100,
        help =
            'The number of bins to use for the non-cumulative histogram.'
        )
    parser.add_option( '', '--title', action = 'store_true', default = False,
        help =
            'Add a title to plots.'
        )
    parser.add_option( '', '--dpi', action = 'store', type = 'float', default = 150.,
        help =
            'Set dpi for plots; default is 150.'
        )

    if argv is None:
        (options, args) = parser.parse_args()
    else:
        (options, args) = parser.parse_args(argv)

    if options.extrapolate is not None and options.extrapolate.upper() == "GAUSSIAN" and options.x_power is None:
        raise ValueError, "If extrapolating a Gaussian, must specify a fit-power."

    return options, args, sys.argv[1:]


# =============================================================================
#
#                       Function Definitions
#
# =============================================================================

def fix_zarray(zarray,sigmas):
    for n,s in enumerate(sorted(sigmas)):
        for ii in range(zarray.shape[0]):
            if zarray[ii,n] >= s:
                zarray[ii,n] = s
        #for jj in range(zarray.shape[1]):
        #    if zarray[n,jj] >= s:
        #        zarray[n,jj] = s
        #for ii in range(zarray.shape[0]):
        #    if zarray[ii,-(n+1)] >= s:
        #        zarray[ii,-(n+1)] = s


def createDataHolder( tableName, columns = None ):
    """
    Creates a DataHolder object. If tableName is the same as a table in lsctables, the DataHolder class will be an instance of that table's RowType.
    """
    if tableName in lsctables.TableByName:
        base = lsctables.TableByName[ tableName ].RowType
    else:
        base = object
    # define the class
    class DataHolder( base ):
        def __init__(self, columns = None):
            # override the __slots__ class (if it exists) with columns, if they are specified
            if columns is not None:
                self.__slots__ = columns

        def store( self, data ):
            """
            Takes a list of data and assigns the values to the object's variables.

            @data: a list of tuples in which the first element is the column name and the second is the value to assign.
            """
            for col, val in data:
                setattr( self, col, val )

        def get_value(self, arg):
            """
            Returns the result of some operation on the elements in self. 
            'arg' can be the name of any defined function in self's base class,
            a slot in self, or a function of either or both. 
            See ligolw.dataUtils.get_row_stat for more info.
            """
            return get_row_stat( self, arg )

    return DataHolder

def combineRowStats( function, rows ):
    """
    Performs the desired function on the list of single statistics. Note: this can only combine one statistic from each row.
    @function: can be either a known pre-set (see below) or an arbitrary function. If an arbitrary function, it must be in terms of the ifo names.
    @rows: a dictionary of statistics keyed by the ifos
    """
    # check if the function is a known pre-sets
    if function == 'sum':
        return sum(rows.values())
    if function == 'quad_sum':
        return math.sqrt(sum([x**2. for x in rows.values()]))
    if function == 'min':
        return min(rows.values())
    if function == 'max':
        return max(rows.values())
    if function == 'mean':
        return numpy.mean(numpy.array(rows.values()))
    if function == 'median':
        return numpy.median(numpy.array(rows.values()))
    if function == 'alpha_min':
        return rows[min(rows.keys())]
    if function == 'sorted_keys':
        return ','.join(sorted(rows.keys()))
    if function == 'sorted_values':
        return ';'.join(sorted(map( str, rows.values() )))
    if function == 'echo':
        return rows

    # otherwise, evaluate the function explicitly
    safe_dict = dict([ [name,val] for name,val in rows.items() + math.__dict__.items() if not name.startswith('__') ])

    try:
        return eval( function, {"__builtins__":None}, safe_dict )
    except NameError:
        # this can happen if an ifo that's specified in the combining function is not in the coincident ifos; in this case, just return None
        return None

class OffsetVector(dict):
    weak_equality = False
    def __init__(self, offset_dict):
        for ifo in offset_dict:
            self[ifo] = offset_dict[ifo]

    def __eq__(self, other):
        """
        The default equality test is to consider two vectors to be equal only if all ifos are the same and all offsets are the same. If one vector is a subset of the other vector, they will not be considered equal. However, if the class attribute weak_equality is set to True, only offsets of the ifos that are both in self and other will be checked. For example:
        >>> a = OffsetVector({'H1': 0, 'L1': 5})
        >>> b = OffsetVector({'H1': 0, 'L1': 5, 'V1': 10})
        >>> a == b
        False
        >>> OffsetVector.weak_equality = True
        >>> a == b
        True
        """
        if type(other) != type(self):
            return False
        if OffsetVector.weak_equality:
            return all( self[ifo] == other[ifo] for ifo in set(self.keys()) & set(other.keys()) )
        else:
            return self.__hash__() == other.__hash__()

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        if OffsetVector.weak_equality:
            return 1
        else:
            return hash(tuple(sorted(self.items())))


class Category:
    """
    Class to store category information.
    """
    default_match_criteria = ['offset_vector', 'datatype', 'veto_cat', 'on_instruments', 'ifos', 'param_group']

    def __init__(self, offset_vector = {}, datatype = None, veto_cat = None, on_instruments = frozenset([u'ALL']), ifos = frozenset([u'ALL']), param_group = None):
        self.offset_vector = OffsetVector(offset_vector)
        self.datatype = unicode(datatype)
        self.veto_cat = unicode(veto_cat)
        self.on_instruments = frozenset(on_instruments)
        self.ifos = frozenset(ifos)
        self.param_group = param_group
        self.livetime = 0

    def add_livetime(self, time):
        self.livetime += time

    def get_livetime(self, time_units = 'yr'):
        return sqlutils.convert_duration( self.livetime, time_units )

    def selective_eq(self, other, check_me):
        """
        Only checks the values listed in check_me to figure out whether or not self is equal to other.
        """
        if type(other) != type(self):
            return False
        return all(getattr(self,x) == getattr(other,x) for x in check_me)

    def __eq__(self, other):
        """
        For default equality check, uses class attribute default_match_criteria to check what parameters should be considered.
        """
        b = type(self) == type(other) and self.__hash__() == other.__hash__()
        if b and OffsetVector.weak_equality and 'offset_vector' in Category.default_match_criteria:
                b = self.offset_vector == other.offset_vector
        return b

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(tuple(getattr(self,x) for x in Category.default_match_criteria))


class Data( dict ):
    """
    Class to store statistics and livetime for plotting.
    """
    class DataElement:
        """
        Sub-class to store individual data elements.

        @categories: a list of instances of the Category class defining which categories this data element belongs to
        @data: an instance of the DataHolder class listing statistics and methods associated with this element
        """
        def __init__(self, thisid, categories, data):
            self._id = thisid
            self.categories = categories
            self.data = data
            self.cumrates = {}

        def update(self, _id = None, categories = [], data = None,  addToExistingCat = True):
            # update id
            if _id is not None:
                self._id = _id
            # update data
            if data is not None:
                self.data = data
            # update categories
            if addToExistingCat:
                self.categories.extend(categories)
            else:
                self.categories = categories

    def __init__(self):
        """
        A list of all the data elements is kept as an index.
        """
        self.data_index = {}

    def add_data(self, _id, categories, data):
        """
        Adds a new DataElement to self.

        @id: some unique value to identify the data element
        @categories: a list of categories that this data falls in. If one or more of these categories are equal (equality determined by the default Category match_criteria) to a category already in all_categories, the category is set to that category. This results in distinct categories only being saved once in memory, with all DataElements that share that category pointing to the same memory address.
        """
        d = self.DataElement( _id, categories, data )
        self.data_index[_id] = d
        for c in categories:
            self.setdefault(c, [])
            self[c].append( d )

    def update(self, _id, categories = [], data = None, addToExistingCat = True, errOnMissing = True):
        """
        Updates all DataElements in self that have the given id. If no DataElement is found with the given id and errOnMissing is False, adds a new entry.
        """
        if _id not in self.data_index:
            if errOnMissing:
                raise ValueError, "An element with id %s could not be found." % str(_id)
            else:
                self.add_data( _id, categories, data )
        else:
            self.data_index[_id].update( categories = categories, data = data, addToExistingCat = addToExistingCat)
            self.refresh_categories( [self.data_index[_id]] )

    def refresh_categories( self, data_list ):
        """
        Updates self with categories that are in data_list that are not in self.
        """
        for d in data_list:
            self.data_index[d._id] = d
            for c in d.categories:
                self.setdefault(c, [])
                if d not in self[c]:
                    self[c].append(d)
    
    def add_livetime(self, livetime, category, match_criteria = []):
        """
        Adds livetime to all categories in self that match the given criteria.
        """
        if match_criteria == []:
            match_criteria = Category.default_match_criteria
        for cat in [cat for cat in self if cat.selective_eq(category, match_criteria)]:
            cat.livetime += livetime

    def get_livetime(self, category, match_criteria = [], time_units = 'yr'):
        """
        Returns the sum of all the livetimes of categories that match the given category via the given match_criteria.
        """
        if match_criteria == []:
            match_criteria = Category.default_match_criteria
        return sqlutils.convert_duration(sum([cat.livetime for cat in self if cat.selective_eq(category, match_criteria)]), time_units)

    def create_background(self, match_criteria = []):
        """
        Creates background categories out of the slide categories and adds this to all slide elements' categories lists. Default action is to create a background for each veto-category, on_instruments, ifos, and param_group. However, this can be overridden with the match_criteria argument.
        """
        if match_criteria == []:
            match_criteria = ['veto_cat', 'on_instruments', 'ifos', 'param_group']
        for vals in set([ tuple(getattr(c, x) for x in match_criteria) for c in self if c.datatype == 'slide' ]):
            # create the background category
            bkg_cat = Category( offset_vector = {}, datatype = 'background' )
            [setattr(bkg_cat, x, y) for x, y in zip(match_criteria, vals)]
            bkg_cat.livetime = sum([c.livetime for c in self if c.datatype == 'slide' and bkg_cat.selective_eq(c, match_criteria)  ])
            # add this background category to each matching slide's categories
            self[bkg_cat] = list(set([x for c in self if c.datatype == 'slide' and bkg_cat.selective_eq(c, match_criteria) for x in self[c]]))

    def compute_cumrates(self, stat, foreground_datatype, rank_by = 'max', group_by = [], num_slides = 100.):
        """
        Computes the cumulative rates for all the distinct groups that exist in self. Distinct groups are determined by group_by.
        """
        if group_by == []:
            group_by = ['datatype', 'veto_cat', 'on_instruments', 'ifos', 'param_group']
        distinct_groups = set([ tuple(getattr(c,x) for x in group_by) for c in self])
        for group in distinct_groups:
            this_group = Category()
            [setattr(this_group, x, y) for (x,y) in zip( group_by, group )]
            this_group.livetime = self.get_livetime( this_group, group_by, time_units = 's' )
            # get the list of all stats that fall in this category
            this_data = []
            for c in self:
                if c.selective_eq(this_group, group_by):
                    this_data.extend( self[c] )
            this_data = sorted(set(this_data), key = lambda x: getattr(x.data, stat), reverse = rank_by == 'min')
            d = [getattr(x.data, stat) for x in this_data]
            #   we need to figure out the number of trials that were done in this category; we do this by taking the ratio
            #   of the this category's livetime to the foreground datatype associated with this category's livetime
            # temporarily set this_group's datatype to foreground_datatype in order to get the right livetime
            orig_dt = this_group.datatype
            this_group.datatype = foreground_datatype
            fg_livetime = self.get_livetime( this_group, match_criteria = 'datatype' not in group_by and group_by+['datatype'] or group_by, time_units = 's' )
            nTrials =  float(this_group.livetime) / fg_livetime
            # set datatype back to correct
            this_group.datatype = orig_dt
            # compute the cum-rates
            these_cumrates = [ (len(d) - bisect.bisect_left(d, x))/nTrials for x in d ]
            # assign to data
            for data_elem, rate in zip( this_data, these_cumrates ):
                data_elem.cumrates[this_group] = rate
       
    def get_cumrates(self, group, stat, rank_by ='max'):
        """
        Returns a sorted list (by stat) of stats, cumrates, and ids for the given group.
        """
        return sorted([(getattr(d.data, stat), d.cumrates[group], d._id) for d in self.data_index.values() if group in d.cumrates], reverse = rank_by == 'min')

    def get_data(self, _id = None, category = None, category_match_criteria = []):
        """
        Returns a list of DataElements that matches a given id, a given category, or both. If category_match_criteria is specified, will get data that matches the specified elements in category. Otherwise, will use Category.default_match_criteria for comparing category to the stored categories.
        """
        if category_match_criteria == []:
            category_match_criteria = Category.default_match_criteria
        return set([x for c in self for x in self[c] if (category is None or c.selective_eq(category, category_match_criteria)) and (_id is None or _id == x._id)])

    def get_categories(self, category, match_criteria = []):
        """
        Returns a list of categories in self that match the given category via the match_criteria.
        """
        if match_criteria == []:
            match_criteria = Category.default_match_criteria
        return [x for x in self if x.selective_eq(category, match_criteria)] 

    def collapse(self, args):
        """
        Cycles over the DataElements in self, keeping only the given args.
        
        @args: A list of tuples. In each tuple, the first element is the name to give the new collapsed value and the second element is the argument to carry out (either a name or a function) on the uncollapsed row to get the collapsed value.
        """
        cols = [arg[0] for arg in args]
        fns = [arg[1] for arg in args]
        collapsedRow = createDataHolder( 'collapsedRow' )
        for n,origElement in enumerate(self.data_index.values()):
            d = collapsedRow( cols )
            d.store([(col, origElement.data.get_value(fn)) for col, fn in zip(cols, fns)]) 
            origElement.data = d
       
def combineData(dataObj, match_column, args, param_grouping_function, verbose = False):
    """
    Cycles over the DataElements in dataObj, combining any DataElements with the same match_column value via the given args and returns a new Data object in which the element's ids are the values of the match_column. Note: the categories of the DataElements in the new Data object are just the concatenation of the older objects individual categories. These might need to be updated depending on the paramters of the newer category.

    @dataObj: the instace of Data to carry the combination on
    @match_column: name of column in the DataElements to use to match rows to combine; e.g., 'coinc_event_id'
    @args: a list of tuples. In each tuple the first element is the name to give the new combined value, the second element is the column in each row to identify that row by, the third is the column or function of columns in each row to combine, and the final element is the way to combine them, which can be either a predefined method or a function in terms of values of the first element. For example, if you wanted the average chirp mass and the sum of the squares of new_snr over H1 and L1, the args should look like:
        args = [ (combined_newsnr_sq, ifo, get_new_snr, H1**2.+L1**2.), (combined_mchirp, ifo, mchirp, mean) ]
    """
    cols = [arg[0] for arg in args]
    colkeys = [arg[1] for arg in args]
    sngl_stats = [arg[2] for arg in args]
    cmb_fncs = [arg[3] for arg in args]
    newData = Data()
    combinedRow = createDataHolder( 'combinedRow' )
    # get the unique match values
    match_vals = {}
    for d in dataObj.data_index.values():
        this_id = d.data.get_value(match_column)
        match_vals.setdefault(this_id, [])
        match_vals[this_id].append(d)
    ii = 0
    for idcol, combine_data in match_vals.items():
        ii += 1
        if verbose:
            if ii != len(match_vals):
                print "%i/%i (%.2f%%)\r" % (ii, len(match_vals), 100*float(ii)/len(match_vals)),
            else:
                print '' 
        newRow = combinedRow( cols )
        stats = [ dict([ [x.data.get_value(colkey), x.data.get_value(snglstat)] for x in combine_data ]) for colkey, snglstat in zip(colkeys, sngl_stats) ]
        newRow.store( [( col, combineRowStats( fnc, stat_dict )) for col, fnc, stat_dict in zip(cols, cmb_fncs, stats)] )
        orig_cats = [y for x in combine_data for y in x.categories] 
        ifos_param = 'ifos' in dir(newRow) and 'ifos' or 'ifo'
        new_cats = [Category( c.offset_vector, c.datatype, c.veto_cat, c.on_instruments, getattr(newRow, ifos_param), param_grouping_function(newRow.param) ) for c in orig_cats]
        newData.add_data( id(newRow), new_cats, newRow )

    return newData

def join_experiment_tables_to_sngl_table( table ):
    return """ 
    JOIN
        experiment, experiment_summary, experiment_map, coinc_event_map 
    ON ( 
        experiment.experiment_id == experiment_summary.experiment_id
        AND experiment_summary.experiment_summ_id == experiment_map.experiment_summ_id
        AND experiment_map.coinc_event_id == coinc_event_map.coinc_event_id
        AND coinc_event_map.event_id == %s.event_id )""" % table

def coinc_in_filter( on_instruments, ifos, filter, sngl_ifos = False ):
    if on_instruments is None or on_instruments == 'ALL':
        on_instruments = frozenset(['ALL'])
    if not (isinstance(on_instruments,set) or isinstance(on_instruments, frozenset)):
        on_instruments = lsctables.instrument_set_from_ifos(on_instruments)
    # following means to only check the on_instruments
    if ifos is None or ifos == 'ALL':
        return frozenset(on_instruments) in filter or frozenset(['ALL']) in filter

    if not (isinstance(ifos,set) or isinstance(ifos, frozenset)):
        if sngl_ifos:
            ifos = set([ifos])
        else:
            ifos = lsctables.instrument_set_from_ifos(ifos)
    # if given on_instruments is 'ALL', means to only check if ifos are in any filter
    if on_instruments == frozenset(['ALL']):
        if sngl_ifos:
            return set(['ALL']) in filter.values() or any(ifos.issubset(x for x in filter.values()))
        else:
            return set(['ALL']) in filter.values() or ifos in filter.values()
    # do full test
    if frozenset(['ALL']) in filter: # on-instruments don't matter, just check ifos
        if sngl_ifos:
            return set(['ALL']) in filter[frozenset(['ALL'])] or ifos.issubset(filter[frozenset(['ALL'])])
        else:
            return set(['ALL']) in filter[frozenset(['ALL'])] or ifos in filter[frozenset(['ALL'])]
    elif frozenset(on_instruments) in filter:
        if sngl_ifos:
            return set(['ALL']) in filter[frozenset(on_instruments)] or ifos.issubset(filter[frozenset(on_instruments)])
        else:
            return set(['ALL']) in filter[frozenset(on_instruments)] or ifos in filter[frozenset(on_instruments)]
    # if get to here, fails all tests
    else:
        return False
        
def poissonPDFlog10( k,lmbda ):
    return k*numpy.log10( lmbda ) - (lmbda + special.gammaln(k+1.))/numpy.log(10.)

def ColorFormatter(y, pos):
    return "$10^{%.1f}$" % y

#
#   Helper functions to compute stats on the fly when querying the database
#

class checkInstruments:
    def __init__(self, filter_name, filter, is_single = False):
        self.filter_name = filter_name
        self.filter = filter
        self.is_single = is_single
    def apply_test(self, on_instruments, ifos):
        return coinc_in_filter( on_instruments, ifos, self.filter, sngl_ifos = self.is_single )
    def create_apply_test(self, sqlite_connection):
        sqlite_connection.create_function(self.filter_name, 2, self.apply_test)


def createCombineRowsMethod( tableName, columns, functionList ):
    """
    Creates a CombineRows class that can be used in a sqlite database to combine rows on the fly. Takes in a sngl_function, which is the function used to combine columns within a single row, and a combining_function, which is the function used to combine the results of the sngl_functions across rows.

    @tableName: the name of the table that will be reading from. If it is a table in lsctables.py, all methods and columns from that table will be inherited.
    @columns: the list of columns that will be storing data to. This list must be in the same order that will be reading data in from the database with.
    @functionList: a list of tuples. The first item should be the combining function to use, in terms of the ifos to combine, and the second item should be the sngl function to use, in terms of columns or methods of the sngl_row.
    """

    sngl_row = createDataHolder(tableName, columns)

    class CombineRows:
        def __init__(self):
            """
            Initializes variables needed for the step process.
            """
            self.this_coinc = dict([ [x, {}] for x in functionList ])

        def step(self, *args):
            """
            Populates self.this_coinc
            """
            this_row = sngl_row(columns)
            this_row.store(zip(columns,args))
            for combine_func, sngl_function in functionList:
                self.this_coinc[(combine_func, sngl_function)][this_row.ifo] = this_row.get_value(sngl_function)

        def finalize(self):
            """
            Once all the singles for the coinc have been gathered, applies the desired combining function(s) to them and returns the result. Results are returned as a comma seperated string.
            """
            return ','.join([str(combineRowStats( cfunc, self.this_coinc[(cfunc, sfunc)] )) for cfunc, sfunc in functionList])

    return CombineRows


class getRowStatsFromDatabase:
    """
    Helper class to get row stat directly from the database.
    """
    def __init__(self, tableName, columns, functionList):
        """
        @tableName: see createCombineRowsMethod.
        @columns: see createCombineRowsMethod.
        @functionList: a list of functions to retrieve from the columns in tableName.
        """
        self.dataRow = createDataHolder(tableName, columns)
        self.functionList = functionList
        self.columns = columns
    def get_dbrow_stat(self, *args ):
        """
        Gets the stat.

        @args: data must be in the same order as columns.
        """
        this_row = self.dataRow(self.columns)
        this_row.store(zip(self.columns,args))
        return ','.join([ str(this_row.get_value(func)) for func in self.functionList ])


# =============================================================================
#
#                                     Main
#
# =============================================================================

#
#       Generic Initilization
#
def main(argv = None):

    opts, filenames, input_args = parse_command_line(argv)
    for this_cache in opts.cache_file:
        this_cache = lal.Cache().fromfile(file(this_cache))
        filenames.extend([x.path() for x in this_cache])
    sqlite3.enable_callback_tracebacks(True)

    statistic, stat_label = len(opts.statistic.split(':')) == 2 and opts.statistic.split(':') or [opts.statistic, opts.statistic.replace('_', ' ').title()]
    data_table = sqlutils.validate_option(opts.data_table)

    # get any added background points
    add_background = []
    for addme in opts.add_background_file:
        # compute combined_newsnr
        f = open(addme, 'r')
        bkg_time = f.readline()
        stats = sorted([ float(line) for line in f ])
        f.close()
        lbl = ''
        if len(addme.split(':')) > 1:
            lbl = addme.split(':')[1]
        add_background.append((stats, float(bkg_time) / 3.15569e7, lbl))
    add_background_sp = []
    for addme in opts.add_special_background_file:
        # compute combined_newsnr
        f = open(addme, 'r')
        bkg_time = f.readline()
        stats = sorted([ float(line) for line in f ])
        f.close()
        lbl = ''
        if len(addme.split(':')) > 1:
            lbl = addme.split(':')[1]
        add_background_sp.append((stats, float(bkg_time)/ 3.15569e7, lbl))

    combining_function = None
    if opts.combining_method is not None:
        combining_function = opts.combining_method
    if opts.combining_function is not None:
        combining_function = opts.combining_function

    # get foreground datatype
    fg_datatype = sqlutils.validate_option( opts.foreground_datatype, lower = True )
    if fg_datatype not in lsctables.ExperimentSummaryTable.datatypes:
        raise ValueError, "unrecognized foreground-datatype %s" % opts.foreground_datatype
    # for the sqlqueries:
    if fg_datatype != 'slide':
        add_foreground = ' OR experiment_summary.datatype == "%s"' % fg_datatype
    else:
        add_foreground = ''

    # Get param and param-ranges if specified
    if opts.param_name:
        param_name, param_label = len(opts.param_name.split(':')) == 2 and opts.param_name.split(':') or [opts.param_name, opts.param_name.replace('_', ' ').title()]
        param_parser = sqlutils.parse_param_ranges( data_table, param_name, opts.param_ranges, verbose = opts.verbose )

    # Get plot-special ranges if desired
    if opts.plot_special_time != []:
        plot_special_times = segments.segmentlist([ segments.segment( t - opts.plot_special_window, t + opts.plot_special_window ) for t in opts.plot_special_time ])
       
    # Get coincs to include/exclude
    #FIXME: parse_coinc_options should use frozenset(['ALL]), not 'ALL'
    if opts.exclude_coincs:
        exclude_filter = sqlutils.parse_coinc_options( opts.exclude_coincs ).coinc_types
        if 'ALL' in exclude_filter:
            exclude_filter[frozenset(['ALL'])] = exclude_filter['ALL']
            del exclude_filter['ALL']
        for on_inst, coinc_ifos in exclude_filter.items():
            exclude_filter[on_inst] = [ ifos == 'ALL' and set(['ALL']) or ifos for ifos in coinc_ifos ]
    if opts.include_only_coincs:
        include_filter = sqlutils.parse_coinc_options( opts.include_only_coincs ).coinc_types
        if 'ALL' in include_filter:
            include_filter[frozenset(['ALL'])] = include_filter['ALL']
            del include_filter['ALL']
        for on_inst, coinc_ifos in include_filter.items():
            include_filter[on_inst] = [ ifos == 'ALL' and set(['ALL']) or ifos for ifos in coinc_ifos ]

    # since we don't care about offset vectors, remove it from Category's default_match_criteria
    Category.default_match_criteria = [x for x in Category.default_match_criteria if x != 'offset_vector']

    # create things needed to store data
    data = Data()
    rowType = createDataHolder( 'data_holder', ['end_time', stat_label] )
    last_table_type = ''
    gps_start_time = numpy.inf
    gps_end_time = -numpy.inf
    analyzed_instruments = set()
    num_slides = None

    if opts.verbose:
        print >> sys.stdout, "Analyzing file:"
    for filename in filenames:
        if opts.verbose:
            print >> sys.stdout, "\t%s" % filename
        working_filename = dbtables.get_connection_filename( filename, tmp_path = opts.tmp_space, verbose = False )
        connection = sqlite3.connect( working_filename )
        if opts.tmp_space:
            dbtables.set_temp_store_directory(connection, opts.tmp_space, verbose = False)

        # figure out the data table type
        data_cols = sqlutils.get_column_names_from_table(connection, data_table)
        if 'coinc_event_id' in data_cols:
            table_type = 'coinc'
        if 'event_id' in data_cols:
            table_type = 'sngl'
        
        # sanity checks
        if table_type == 'coinc' and 'ifos' not in data_cols or table_type == 'sngl' and 'ifo' not in data_cols:
            raise ValueError, 'Could not find %s column in the data table (%s)' % (table_type == 'sngl' and 'ifo' or 'ifos', data_table)
        if last_table_type != '' and last_table_type != table_type:
            raise ValueError, '%s table in file %s does not match the types in the other files' %( data_table, filename )
        if table_type == 'coinc' and opts.plot_special_time != [] and opts.single_table is None:
            raise ValueError, "Must specify a --single-table if using --plot-special-time with a coinc-ifos table."
        if table_type == 'sngl' and combining_function is None:
            raise ValueError, "Must provide a combining function or method if querying sngl-ifo tables."

        #
        #   create/update the background categories in data
        #
        if opts.verbose:
            print "\tgetting background categories..."

        # if excluding instrument times or coincs, create the test for them
        in_desired_times = ''
        if opts.include_only_coincs:
            db_include_filter = checkInstruments( 'in_include_filter', include_filter, is_single = False )
            db_include_filter.create_apply_test(connection)
            in_desired_times = ' AND in_include_filter(experiment.instruments, "ALL")'
        if opts.exclude_coincs:
            db_exclude_filter = checkInstruments( 'in_exclude_filter', exclude_filter, is_single = False )
            db_exclude_filter.create_apply_test(connection)
            in_desired_times = ' '.join([ in_desired_times, 'AND NOT in_exclude_filter(experiment.instruments, "ALL")'])

        # establish how many independent backgrounds we'll need
        tsids = set()
        match_criteria = ['datatype', 'veto_cat', 'param_group']
        if opts.plot_by_instrument_time:
            match_criteria.append( 'on_instruments' )
        if opts.param_name:
            n_param_groups = len(param_parser.param_ranges)
        else:
            n_param_groups = 1
        sqlquery = """
            SELECT
                experiment.instruments,
                experiment.gps_start_time,
                experiment.gps_end_time,
                experiment_summary.time_slide_id,
                experiment_summary.veto_def_name,
                experiment_summary.duration,
                experiment_summary.datatype
            FROM
                experiment_summary
            JOIN
                experiment
            ON
                experiment.experiment_id == experiment_summary.experiment_id
            WHERE
                (experiment_summary.datatype == 'slide'""" + add_foreground + ')' + in_desired_times
        if opts.debug:
            print >> sys.stderr, "Sqlite query used to get categories:"
            print >> sys.stderr, sqlquery
        for on_instruments, this_gps_start, this_gps_end, tsid, veto_cat, duration, datatype in connection.cursor().execute(sqlquery):
            on_instruments = frozenset(lsctables.instrument_set_from_ifos(on_instruments))
            analyzed_instruments.update( on_instruments )
            gps_start_time = min(gps_start_time, this_gps_start)
            gps_end_time = max(gps_end_time, this_gps_end)
            if datatype == "slide":
                datatypes = ["slide", "background"]
                if opts.plot_special_time != []:
                    datatypes.append("background_special")
                # store the slide tsid's to keep track of the number of slides
                tsids.add(tsid)
            else:
                datatypes = [datatype]
            if opts.plot_by_ifos:
                 ifo_list = [frozenset([x for x in ifos]) for ifos in CoincInspiralUtils.get_ifo_combos(list(on_instruments))]
            else:
                ifo_list = [frozenset(['ALL'])]
            # if there is more than one param-group or we are plotting by ifos, add the needed backgrounds
            for datatype in datatypes:
                for param_group in range(n_param_groups):
                    for ifos in ifo_list:
                        this_cat = Category( offset_vector = {}, datatype = datatype, veto_cat = veto_cat, ifos = ifos, param_group = param_group)
                        if opts.plot_by_instrument_time:
                            this_cat.on_instruments = on_instruments
                        # add the background category
                        data.setdefault( this_cat, [] )
                        # add the livetime
                        # if doing plot special, subtract the duration of the removed time from the background livetime.
                        # Note that this isn't exactly correct: the range could overlap a slid veto, in which case only some of the time should be subtracted.
                        # However, we assume this is a small affect.
                        if datatype == "background_special" and opts.plot_special_time != [] and plot_special_times.intersects_segment( segments.segment(this_gps_start, this_gps_end) ):
                            data.add_livetime( duration - len(on_instruments)*sum([abs(seg) for seg in plot_special_times]), this_cat, match_criteria )
                        else:
                            data.add_livetime( duration, this_cat, match_criteria )
        
        # add and check the slide count
        # FIXME: We don't have to have the same number of slides in each database in order to do this. However, for that to happen,
        # weighted averages have to be computed.
        if num_slides is None and len(tsids) != 0:
            num_slides = len(tsids)
        elif len(tsids) != num_slides and len(tsids) != 0:
            raise ValueError, "Database %s has a different number of slides (%i) than the other databases (%i)." %( filename, len(tsids), num_slides )

        # set the slide's livetime to the average background time; we do this because all slides have been thrown into the same bucket
        if len(tsids) != 0.:
            for slide_cat in data.get_categories( Category(datatype = 'slide'), match_criteria = ['datatype'] ):
                slide_cat.livetime = float(slide_cat.livetime) / len(tsids)

        #
        #   retrieve the desired statistic
        #
        if opts.verbose:
           print >> sys.stdout, "\tgetting data..." 
       
        quick_query = False
        if table_type == 'sngl':
            # create a function list of functions we'll need to combine
            functionList = []
            if opts.param_name:
                functionList.append((opts.param_combining_method, param_name))
            functionList.extend([('alpha_min', 'end_time+end_time_ns*1e-9'), ('mean', 'mtotal'), (combining_function, statistic)])
            if opts.plot_special_time != []:
                functionList.append( ('sorted_values', 'end_time+end_time_ns*1e-9') )
            get_combined_data = createCombineRowsMethod( data_table, data_cols, functionList )
            connection.create_aggregate('get_combined_data', len(data_cols), get_combined_data)
            # create an aggregate function to get the ifos
            get_ifos = createCombineRowsMethod( 'ifo_holder', ['ifo'], [('sorted_keys', 'ifo')] )
            connection.create_aggregate('get_ifos', 1, get_ifos)
             
            #  If including or excluding coincs, need to create a function to combine ifos as well as discriminate single ifos.
            #  By not allowing any single ifos by that are not in the desired coincident ifos, we're able to pick out doubles from triples (if desired).
            #  We do this by taking the coinc filters and breaking up the coincident ifos that are in them, then pass them through coinc_in_filter
            add_sngl_filter = ''
            add_coinc_filter = ''
            if opts.include_only_coincs:
                include_sngl_filter = dict([ [on_inst, set([coinc == 'ALL' and 'ALL' or ifo for coinc in coincs for ifo in coinc])] for on_inst,coincs in include_filter.items() ])
                db_include_sngl_filter = checkInstruments( 'in_include_sngl_filter', include_sngl_filter, is_single = True )
                db_include_sngl_filter.create_apply_test(connection)
                add_sngl_filter = ''.join([ 'AND in_include_sngl_filter(experiment.instruments, ', data_table, '.ifo', ') '])
                add_coinc_filter = ''.join([ 'HAVING in_include_filter(experiment.instruments, ifos)' ])
            if opts.exclude_coincs:
                exclude_sngl_filter = dict([ [on_inst, set([coinc == 'ALL' and 'ALL' or ifo for coinc in coincs for ifo in coinc])] for on_inst,coincs in exclude_filter.items() ])
                db_exclude_sngl_filter = checkInstruments( 'in_exclude_sngl_filter', exclude_sngl_filter, is_single = True )
                db_exclude_sngl_filter.create_apply_test(connection)
                add_sngl_filter = ''.join([ add_sngl_filter, 'AND NOT in_exclude_sngl_filter(experiment.instruments, ', data_table, '.ifo', ') '])
                add_coinc_filter = ''.join([ add_coinc_filter, add_coinc_filter == '' and 'HAVING ' or 'AND ', 'NOT in_exclude_filter(experiment.instruments, ifos)'])
            sqlquery = ''.join([ """
                SELECT
                    experiment.instruments,
                    experiment_summary.datatype,
                    experiment_summary.veto_def_name,
                    coinc_event_map.coinc_event_id,
                    get_ifos(""", data_table, """.ifo) AS ifos,
                    get_combined_data(""", ','.join(['.'.join([data_table,col]) for col in data_cols]), """)
                FROM
                    """, data_table, """
                """, join_experiment_tables_to_sngl_table( data_table ), """
                WHERE
                    (experiment_summary.datatype == 'slide'""", add_foreground, """)
                    """, add_sngl_filter, """
                GROUP BY
                    coinc_event_map.coinc_event_id
                """, add_coinc_filter ])

        else:
            # if doing plot special, create a table of the coinc events to ignore in the database
            add_time_filter = ''
            if opts.plot_special_time != []:
                sngl_table = sqlutils.validate_option( opts.single_table )
                # create a function in the database to check whether a single time falls in the plot sepcial windows or not
                def check_sngl_time( end_time, end_time_ns ):
                    return end_time+end_time_ns*1e-9 in plot_special_times
                connection.create_function('check_sngl_time', 2, check_sngl_time)
                # populate the table
                sqlquery = ''.join([ """
                        SELECT DISTINCT
                            coinc_event_map.coinc_event_id AS ceid
                        FROM
                            coinc_event_map
                        JOIN
                            """, sngl_table, """ AS sngl_table
                        ON (
                            coinc_event_map.table_name == '""", sngl_table, """' AND
                            sngl_table.event_id == coinc_event_map.event_id )
                        JOIN
                            experiment_summary, experiment_map
                        ON (
                            experiment_summary.experiment_summ_id == experiment_map.experiment_summ_id AND
                            experiment_map.coinc_event_id == coinc_event_map.coinc_event_id )
                        WHERE
                            experiment_summary.datatype == "slide" AND
                            check_sngl_time( sngl_table.end_time, sngl_table.end_time_ns )
                        """])
                if opts.debug:
                    print >> sys.stderr, "SQLite query used to get conics falling in plot special time:"
                    print >> sys.stderr, sqlquery
                remove_coincs = [ceid[0] for ceid in connection.cursor().execute(sqlquery)]

                #add_time_filter = '\nAND %s.coinc_event_id NOT IN (SELECT ceid FROM skip_coincs)' % data_table
                             
            # speed-up: if param_name and statistic are in the data_cols list, just query them directly
            if (opts.param_name is not None and param_name in data_cols) and statistic in data_cols:
                quick_query = True
                get_stats = ''.join([ data_table, ".", param_name, ", ", data_table, ".end_time+", data_table, ".end_time_ns*1e-9, ", data_table, ".mass, ", data_table, ".", statistic ])
            elif opts.param_name is None and statistic in data_cols:
                quick_query = true
                get_stats = ''.join([ data_table, ".end_time+", data_table, ".end_time_ns*1e-9, ", data_table, ".mass, ", data_table, ".", statistic ])
            else:
                # create method to get the stats and parameter
                functionList = []
                if opts.param_name:
                    functionList.append(param_name)
                functionList.extend(['end_time+end_time_ns*1e-9', 'mass', statistic])
                getStatDB = getRowStatsFromDatabase( data_table, data_cols, functionList )
                connection.create_function( 'get_stat', len(data_cols),  getStatDB.get_dbrow_stat )
                get_stats = ''.join([ 'get_stat(', ','.join([ '.'.join([data_table,col]) for col in data_cols ]), ')' ])

            add_coinc_filter = ''
            if opts.include_only_coincs:
                add_coinc_filter = ''.join([ '\nAND in_include_filter(experiment.instruments, ', data_table, '.ifos)' ])
            if opts.exclude_coincs:
                add_coinc_filter = ''.join([ add_coinc_filter, '\nAND NOT in_exclude_filter(experiment.instruments, ', data_table, '.ifos)'])
            sqlquery = ''.join([ """
                SELECT
                    experiment.instruments,
                    experiment_summary.datatype,
                    experiment_summary.veto_def_name,
                    """, data_table, """.coinc_event_id,
                    """, data_table, """.ifos,
                    """, get_stats, """
                FROM
                    """, data_table, """ 
                """,
                    sqlutils.join_experiment_tables_to_coinc_table( data_table ),
                """
                WHERE
                    (experiment_summary.datatype == 'slide'""", add_foreground, ')', add_coinc_filter, add_time_filter ])
        if opts.debug:
            print >> sys.stderr, "Sqlite query used to get data:"
            print >> sys.stderr, sqlquery
        results = connection.cursor().execute(sqlquery).fetchall()
        for result in results:
            if opts.plot_by_instrument_time:
                on_instruments = frozenset(lsctables.instrument_set_from_ifos( result[0] ))
            else:
                on_instruments = frozenset(["ALL"])
            datatype = result[1]
            # if these are slid triggers, check if any fall in the plot-special window (for sngl tables only)
            if datatype == "slide" and opts.plot_special_time != [] and table_type == 'sngl':
                sngl_end_times = [float(t) for t in result[5].split(',')[-1].split(';')]
                if any(t in plot_special_times for t in sngl_end_times):
                    continue
            veto_cat = result[2]
            ceid = result[3]
            if opts.plot_by_ifos:
                ifos = frozenset(lsctables.instrument_set_from_ifos(result[4]))
            else:
                ifos = frozenset(["ALL"])
            if opts.param_name is not None:
                if quick_query:
                    param_group, end_time, mtotal, stat = [result[5], result[6], result[7], result[8]]
                else:
                    param_group, end_time, mtotal, stat = [float(x) for x in result[5].split(',') if ';' not in x]
                param_group = param_parser.group_by_param_range( param_group )
                # skip this if it doesn't fall in any desired param group
                if param_group is None:
                    continue
            else:
                param_group = 0.
                if quick_query:
                   end_time, mtotal, stat = [result[5], result[6], result[7]] 
                else:
                   end_time, mtotal, stat = [float(x) for x in result[5].split(',') if ';' not in x]
            # Memory improvement: throw out any data less than x-fit and xmin
            if opts.xmin is not None:
                if stat < opts.xmin and opts.min_x_fit is not None and stat < opts.min_x_fit:
                    continue
                elif stat < opts.xmin and opts.min_x_fit is None:
                    continue
            # passed any tests, store data
            this_row = rowType()
            this_row.store([ ('end_time', end_time), ('mtotal', mtotal), ('stat', stat), ('ifos', result[4]) ])
            # construct categories
            categories = [ Category( {}, datatype, veto_cat, on_instruments, ifos, param_group ) ]
            if datatype == "slide":
                categories.append(Category( {}, 'background', veto_cat, on_instruments, ifos, param_group ))
                if opts.plot_special_time != [] and ceid not in remove_coincs:
                    categories.append(Category( {}, 'background_special', veto_cat, on_instruments, ifos, param_group ))
                    
            # add to data
            thisid = id(this_row) 
            data.add_data( thisid, categories, this_row )

        # go on to the next database
        connection.close()
        dbtables.discard_connection_filename( filename, working_filename, verbose = False)

    for cat in data:
        print cat.datatype, len(data[cat])

    #
    #   Sort and calculate cumrate
    #
    if opts.verbose:
        print "Computing cumulative rates..."

    # figure out how to populate the background categories
    group_by = ['veto_cat', 'param_group']
    if opts.plot_by_instrument_time:
        group_by.append( 'on_instruments' )
    if opts.plot_by_ifos:
        group_by.append( 'ifos' )

    data.compute_cumrates( 'stat', fg_datatype, rank_by = 'max', group_by = group_by+['datatype'])

    #
    #   Plot/Compute statistics
    #
    if opts.verbose:
        print >> sys.stdout, "Plotting..."
    opts.gps_start_time = gps_start_time
    opts.gps_end_time = gps_end_time
    opts.ifo_times = ''.join(sorted(analyzed_instruments)) 
    opts.ifo_tag = ''
    opts.enable_output = True
    InspiralUtilsOpts = InspiralUtils.initialise( opts, __prog__, git_version.verbose_msg )
    fnameList = []
    tagList = []
    comments = ''

    for bkg_category in data.get_categories( Category(datatype = 'background'), ['datatype'] ):
        if opts.verbose:
            print >> sys.stdout, "\tgroup: %s, %s, %s, %s" %(bkg_category.veto_cat, str(bkg_category.param_group), ''.join(bkg_category.on_instruments), ''.join(bkg_category.ifos))
        # set a match category to pick out data with
        match_category = Category()
        match_criteria = group_by+['datatype']
        [setattr(match_category, x, getattr(bkg_category, x)) for x in match_criteria]

        # get foreground and  background livetime
        T_bkg = data.get_livetime(bkg_category, time_units = opts.time_units)
        match_category.datatype = fg_datatype
        T_fg = data.get_livetime( match_category, time_units = opts.time_units )

        #
        #
        #   Non-cumulative plot
        #
        #
        plot_non_cum = True
        if plot_non_cum:
            if opts.verbose:
                print >> sys.stdout, "\tcreating histogram..."
            pylab.figure()
            pylab.hold(True)
            xmin = numpy.inf
            xmax = -numpy.inf
            ymin = numpy.inf
            ymax = -numpy.inf

            #
            # create histogram
            #
            # get foreground data
            if fg_datatype != 'slide':
                fg_data = data.get_data( category = match_category, category_match_criteria = match_criteria )
                fg_stats = numpy.array([x.data.stat for x in fg_data])
                plot_fg = fg_stats != numpy.array([])
            else:
                plot_fg = False
            # get background data
            bkg_data = data[bkg_category]
            bkg_stats = numpy.array([x.data.stat for x in bkg_data])
            # FIXME: create an empty plot
            if bkg_stats == numpy.array([]):
                continue

            # get min/max for bins values
            if plot_fg:
                min_val = min( bkg_stats.min(), fg_stats.min() )
                max_val = max( bkg_stats.max(), fg_stats.max() )
            else:
                min_val = min( bkg_stats.min(), bkg_stats.min() )
                max_val = max( bkg_stats.max(), bkg_stats.max() )

            # create the bins
            nbins = opts.nbins
            if opts.lin_x:
                #bins = numpy.linspace(min_val, max_val, nbins+1, endpoint = True)
                bins = numpy.linspace(numpy.power(min_val, 1./opts.x_power), numpy.power(max_val, 1./opts.x_power), nbins+1, endpoint = True)
            else:
                bins = numpy.logspace( numpy.log10(min_val), numpy.log10(max_val), nbins+1, endpoint = True )
            bins = numpy.power(bins, opts.x_power)
            # bump up the last bin by a little bit to make sure there are no numerical errors
            bins[-1] += .0001
            ds = numpy.array([bins[i+1] - bins[i] for i in range(len(bins)-1)])/2.
            xvals = bins[:-1] + ds

            #
            # plot background
            #
            bkg_count, _ = numpy.histogram(bkg_stats, bins)
            # divide the background count by the ratio of the background time to the foreground time
            bkg_count = bkg_count * (T_fg / T_bkg )
            remove_me = pylab.find(bkg_count == 0.)
            bkg_count = numpy.array([bkg_count[i] for i in range(len(bkg_count)) if i not in remove_me])
            bkg_hist_stats = numpy.array([xvals[i] for i in range(len(xvals)) if i not in remove_me])
            bkg_ds = numpy.array([ds[i] for i in range(len(xvals)) if i not in remove_me])

            # reset ymin/ymax
            if bkg_count != []:
                ymin = min(min(bkg_count),ymin)
                ymax = max(max(bkg_count),ymax)

            #
            # compute and plot background error
            #
            bkg_error = numpy.sqrt( bkg_count / (T_bkg/T_fg) )
            pylab.errorbar( bkg_hist_stats, bkg_count, xerr = bkg_ds, yerr = bkg_error, color = 'k', ecolor = 'gray', marker = 'o', markersize = 4, linestyle = 'None', zorder = 2, label = 'Background' )

            #
            # plot foreground
            #
            loudest_events = []
            if plot_fg:

                fg_count, _ = numpy.histogram(fg_stats, bins)
                remove_me = pylab.find(fg_count == 0.)
                fg_count = numpy.array([fg_count[i] for i in range(len(fg_count)) if i not in remove_me])
                fg_hist_stats = numpy.array([xvals[i] for i in range(len(xvals)) if i not in remove_me])

                # reset ymin/ymax
                if fg_count != []:
                    ymin = min(min(fg_count),ymin)
                    ymax = max(max(fg_count),ymax)

                fclr = 'b'
                lbl = fg_datatype.replace( '_', '-').title()
                if fg_datatype == 'simulation':
                    edgclr = 'r'
                elif opts.plot_poisson_pdf:
                    edgclr = 'white'
                else:
                    edgclr = 'none'
                pylab.scatter( fg_hist_stats, fg_count, marker = 'o', edgecolor = edgclr, facecolor = fclr, label = lbl, s = 20, linewidth = .5, zorder = 4 )

                #
                # plot foreground triggers that fall in plot-special window as stars
                #
                if opts.plot_special_time != []:
                    plot_sp_vals = [x.data.stat for x in fg_data if x.data.end_time in plot_special_times]
                    spx = []
                    spy = []
                    if plot_sp_vals != []:
                        # in which bin does each plot_sp_val fall?
                        for x in plot_sp_vals:
                            if x > fg_hist_stats.max():
                                spx.append(fg_hist_stats.max())
                                spy.append(fg_count[-1])
                            else:
                                idxs = set(pylab.find( x >= fg_hist_stats )) & set(pylab.find( x <= fg_hist_stats))
                                spx.extend([ fg_hist_stats[i] for i in idxs ])
                                spy.extend([ fg_count[i] for i in idxs ])
                        pylab.scatter( spx, spy, c = 'yellow', marker = (5,1,0), s = 50, linewidth = .5, label = '_no_legend', zorder = 5 )


            #
            # plot extrapolation
            #
            if plot_fg and opts.extrapolate is not None:
                if opts.verbose:
                    print >> sys.stdout, "\textrapolating background..."
                ### For Gaussian fit:
                    # We're going to fit a Gaussian with 0 mean to the background. Since this is given by y = A*exp(-x^2./(2*sigma^2.))/sqrt(2*pi*sigma^2.), we take the natural log, and thus fit
                    # y = beta[0] + beta[2]*x^2 where beta[0] = ln(A) - ln(2*pi*sigma^2.)/2., beta[2] = -1/(2.*sigma^2.) and y = ln(bkg_count)
                    # To do this, we'll do a weighted linear least squares fit using Beta_hat = (transpose(X).W.X)^(-1).transpose(X).W.Y. where W is the nxn matrix of weights (which are just taken
                    # as the error in each background point), X is the nx2 matrix of [1 x[i]^2.], Y is the nx1 matrix of ln(bkg_cound[i]), and Beta_hat is the 2x1 matrix of fitted betas.
                    # That assumes that the plotted x is linear in whatever the desired statistic is. However, if the plotted x is actually the desired statistic squared, then we should
                    # fit y = beta[0] + beta[2]*x. Generalizing, we actually fit y = beta[0] + beta[2]*x^((2./opts.x_power)), since only the user can know what function of their desired statistic they've
                    # punched into the command line.
                ### For power-law fit:
                    # We're going to fit bkg_count/T_fg = A*statistic^B by doing a least squares fit to y = beta[0] + beta[1]*x where y = ln(bkg_count/T_fg), x = ln(statistic), beta[0] = ln(A) 
                    # and beta[1] = B.
                if opts.min_x_fit is None:
                    min_x_fit = bkg_hist_stats.min()
                else:
                    min_x_fit = opts.min_x_fit
                start_i = pylab.find( bkg_hist_stats >= min_x_fit )[0]
                n = len(bkg_count) - start_i
                X = numpy.matrix(numpy.zeros((n,2), dtype = float))
                Y = numpy.matrix(numpy.zeros((n,1), dtype = float))
                W = numpy.matrix(numpy.zeros((n,n), dtype = float))
                # populate the input matrices
                for i in range(n):
                    if bkg_count[start_i + i] == 0:
                        continue
                    Y[i,0] = numpy.log( bkg_count[start_i + i] )
                    W[i,i] = 1./numpy.power( bkg_error[start_i + i]/bkg_count[start_i + i], 2. )
                    X[i,0] = 1.
                    if opts.extrapolate.upper() == 'POWER':
                        X[i,1] = numpy.log( bkg_hist_stats[start_i + i] )
                    else:
                        X[i,1] = bkg_hist_stats[start_i + i]**(2./opts.x_power)

                # calculate the Beta_hats
                Beta_hat = (X.T * W * X).I * X.T * W * Y

                # Compute the Chisq of the fit
                noncum_chisq = ( W * numpy.power((Y - X * Beta_hat), 2.) ).sum()
                # the number of degrees of freedom = Ndof = number of points - number of fitting params - 1
                noncum_Ndof = n - 2 - 1

                # now that we have the fit parameters, plot the fitted values 
                n_fit = 50
                if opts.lin_x:
                    x_fit = numpy.linspace( bkg_hist_stats[start_i], bkg_hist_stats.max(), num = n_fit )
                else:
                    x_fit = numpy.logspace( numpy.log10(bg_hist_stats[start_i]), numpy.log10( bkg_hist_stats.max() ), num = n_fit )
                Xfit = numpy.matrix(numpy.zeros((n_fit,2), dtype = float))
                for i,x in enumerate(x_fit):
                    Xfit[i,0] = 1.
                    if opts.extrapolate.upper() == 'POWER':
                        Xfit[i,1] = numpy.log(x)
                    else:
                        Xfit[i,1] = x**(2./opts.x_power)
                y_fit = numpy.exp( numpy.array(Xfit * Beta_hat) )

                # plot the extrapolated line
                pylab.plot( x_fit, y_fit, 'g-', linewidth = 2, zorder = 5, label = 'Fitted Bkg.' )

                # reset ymin/ymax
                if y_fit.any():
                    ymin = min(y_fit.min(),ymin)
                    ymax = max(y_fit.max(),ymax)

            #
            # finalize plot settings
            #
            pylab.grid()
            pylab.legend(loc = 'lower left', scatterpoints = 1, numpoints = 1)

            pylab.xlabel( stat_label )
            pylab.ylabel( "Count per Foreground Experiment" )
            if opts.title:
                if opts.param_name is not None:
                    t = 'Histogram %s %s' % (param_label.title(), param_parser.param_range_by_group(bkg_category.param_group))
                else:
                    t = 'Histogram'
                pylab.title(t)

            # set axes limits and scale
            xmin = min(xvals)
            xmax = max(xvals)
            if not opts.lin_x:
                pylab.gca().semilogx()
                xmax = xmax*10**(0.3)
            else:
                xmax = xmax+.1*xmax
            if not opts.lin_y:
                pylab.gca().semilogy()
                ymin, ymax = ymin*10**(-0.5), ymax*10**(0.5)
            else:
                ymin, ymax = ymin - .1*ymin, ymax + .1*ymax
            # overrule with input options
            if opts.xmin is not None:
                xmin = opts.xmin
            if opts.xmax is not None:
                xmax = opts.xmax
            if opts.ymin is not None:
                ymin = opts.ymin
            if opts.ymax is not None:
                ymax = opts.ymax
            

            pylab.xlim( xmin, xmax )
            pylab.ylim( ymin, ymax )
            
            # store plot info
            plot_description = 'F%i' % ( len(fnameList) )
            name = InspiralUtils.set_figure_tag( plot_description, datatype_plotted = fg_datatype.upper(), open_box = 'all_data' in fg_datatype or 'exclude_play' in fg_datatype)
            fname = InspiralUtils.set_figure_name(InspiralUtilsOpts, name)
            #fname = re.sub('.png', '.pdf', fname)
            fname_thumb = InspiralUtils.savefig_pylal( filename=fname, dpi=opts.dpi  )
            fnameList.append(fname)
            tagList.append(name)


        #
        #
        #   Cumulative plot
        #
        #
        if opts.verbose:
            print >> sys.stdout, "\tcreating cumulative plot..."
        pylab.figure()
        pylab.hold(True)
        xmin = numpy.inf
        xmax = -numpy.inf
        ymin = numpy.inf
        ymax = -numpy.inf

        #
        #   XXX: For paper, create text files of data. One file is created
        #   for each data type: regular background (black dots), extended background
        #   (black crosses), foreground (blue triangles), regular background with
        #   the loud event removed (gray dots), and extended background with
        #   the loud event removed (gray crosses).
        #

        # the basic header used in all of the text files
        txthead = """#
    # Data used to generate the %s
    # in the cumulative rate plot (Figure 3) in LIGO Document P1100034.
    #
    # For questions regarding the use of this material, please contact the LSC
    # spokesperson Gaby Gonzalez (gonzalez@lsu.edu) or corresponding authors
    # Collin Capano (collin.capano@ligo.org) and Stephen Privitera (sprivite@caltech.edu)
    #
    """

        #
        # plot special background
        #
        if opts.plot_special_time != []:
            # get background data
            match_category.datatype = 'background_special'
            plotvals = data.get_cumrates( match_category, 'stat', rank_by = 'max' )
            bkgsp_stats = numpy.array([x[0] for x in plotvals])
            bkgsp_cumnum = numpy.array([y[1] for y in plotvals])

            # rest xmmin/xmax
            if bkgsp_stats != []:
                xmin = min(min(bkgsp_stats),xmin)
                xmax = max(max(bkgsp_stats),xmax)
            # reset ymin/ymax
            if bkgsp_cumnum != []:
                ymin = min(min(bkgsp_cumnum / T_fg),ymin)
                ymax = max(max(bkgsp_cumnum / T_fg),ymax)

            #
            # compute and plot background error
            #
            bkgsp_error = numpy.sqrt( bkgsp_cumnum / (T_fg * T_bkg) )
            plot_errs = numpy.array([ bkgsp_error,bkgsp_error ])
            for n, (y, err) in enumerate(zip(bkgsp_cumnum.tolist(),bkgsp_error.tolist())):
                if y/T_fg - err <= 0.0:
                    plot_errs[0,n] = y/T_fg - 1e-15
            pylab.errorbar( bkgsp_stats, bkgsp_cumnum / T_fg, yerr = plot_errs, color = '#777777', ecolor = '#777777', marker = 'o', markersize = 4, markeredgecolor = '#777777', linestyle = 'None', zorder = -1, label = '_no_legend')
            # print data to text file
            f = open('PlotRatesData-HundredSlideBkg-NoEvent.txt', 'w')
            print >> f, txthead % 'background estimate from 100 slides with the candidate event removed (gray dots)'
            print >> f, '#Combined NewSNR\tCumulative Rate (yr^-1)\ty-error (+/- yr^-1)'
            for stat, cumrate, err in zip(bkgsp_stats.tolist(), (bkgsp_cumnum / T_fg).tolist(), bkgsp_error.tolist()): 
                print >> f, '%f\t%f\t%f' % (stat, cumrate, err)
            f.close()
            
            # store the ids of the loudest background events
            loudest_bkgsp_events = []
            ii = 0
            last_stat = None
            for x in plotvals[::-1]:
                if ii > 10:
                    break
                if x[0] != last_stat:
                    ii += 1
                loudest_bkgsp_events.append(x[2])

            #
            # plot added background triggers
            #
            added_bkgsp = []
            for newsnrs, bkg_time, lbl in add_background_sp:
                nTrials =  bkg_time / T_fg
                # compute the cum-rates
                these_cumnums = numpy.array([ (len(newsnrs) - bisect.bisect_left(newsnrs, x))/nTrials for x in newsnrs ])
                this_error = numpy.sqrt( these_cumnums / (T_fg * bkg_time) )
                plot_errs = numpy.array([this_error, this_error])
                for n, (y, err) in enumerate(zip(these_cumnums.tolist(),this_error.tolist())):
                    if y/T_fg - err <= 0.0:
                        plot_errs[0,n] = y/T_fg - 1e-15
                T_add_bkgsp = bkg_time
                pylab.errorbar( newsnrs, these_cumnums / T_fg, yerr = plot_errs, color = '#777777', ecolor = '#777777', marker = 'x', markersize = 6, mew = 1, linestyle = 'None', zorder = -1, label = '_no_legend' )

                # print data to text file
                f = open('PlotRatesData-ExtendedBkg-NoEvent.txt', 'w')
                print >> f, txthead % 'extended background estimate with the candidate event removed (gray crosses)'
                print >> f, '#Combined NewSNR\tCumulative Rate (yr^-1)\ty-error (+/- yr^-1)'
                for stat, cumrate, err in zip(newsnrs, (these_cumnums / T_fg).tolist(), this_error.tolist()): 
                    print >> f, '%f\t%f\t%f' % (stat, cumrate, err)
                f.close()

                #if lbl != '':
                #    if opts.lin_y:
                #        ytext = yval+.05*(ymax-ymin)
                #    else:
                #        ytext = yval * 10**.25
                #    pylab.annotate( lbl, (xval,yval), xytext=(xval,ytext) )
                # reset xmin/xmax
                xmin = min(min(newsnrs), xmin)
                xmax = max(max(newsnrs), xmax)
                ymin = min(these_cumnums.min(), ymin)
                ymax = max(these_cumnums.max(), ymax)
                # save x, y points for poisson pdf calc.
                added_bkgsp = [x for x in zip(newsnrs, these_cumnums) if x[0] > bkgsp_stats.max()] #.append((xval, T_fg / pt_livetime))

        #
        # plot background
        #

        # get background data
        match_category.datatype = 'background'
        plotvals = data.get_cumrates( match_category, 'stat', rank_by = 'max' )
        bkg_stats = numpy.array([x[0] for x in plotvals])
        bkg_cumnum = numpy.array([y[1] for y in plotvals])

        # rest xmmin/xmax
        if bkg_stats != []:
            xmin = min(min(bkg_stats),xmin)
            xmax = max(max(bkg_stats),xmax)
        # reset ymin/ymax
        if bkg_cumnum != []:
            ymin = min(min(bkg_cumnum / T_fg),ymin)
            ymax = max(max(bkg_cumnum / T_fg),ymax)

        #
        # compute and plot background error
        #
        bkg_error = numpy.sqrt( bkg_cumnum / (T_fg * T_bkg) )
        plot_errs = numpy.array([ bkg_error,bkg_error ])
        for n, (y, err) in enumerate(zip(bkg_cumnum.tolist(),bkg_error.tolist())):
            if y/T_fg - err <= 0.0:
                plot_errs[0,n] = y/T_fg - 1e-15
        pylab.errorbar( bkg_stats, bkg_cumnum / T_fg, yerr = plot_errs, color = 'k', ecolor = 'k', marker = 'o', markersize = 4, markeredgecolor = 'k', linestyle = 'None', zorder = 2, label = 'Background' )

        # print data to text file
        f = open('PlotRatesData-HundredSlideBkg.txt', 'w')
        print >> f, txthead % 'background estimate from 100 slides (black dots)'
        print >> f, '#Combined NewSNR\tCumulative Rate (yr^-1)\ty-error (+/- yr^-1)'
        for stat, cumrate, err in zip(bkg_stats.tolist(), (bkg_cumnum / T_fg).tolist(), bkg_error.tolist()): 
            print >> f, '%f\t%f\t%f' % (stat, cumrate, err)
        f.close()
            
        # store the ids of the loudest background events
        loudest_bkg_events = []
        ii = 0
        last_stat = None
        for x in plotvals[::-1]:
            if ii > 10:
                break
            if x[0] != last_stat:
                ii += 1
            loudest_bkg_events.append(x[2])

        #
        # plot added background triggers
        #
        added_bkg = []
        for newsnrs, bkg_time, lbl in add_background:
            nTrials =  bkg_time / T_fg
            # compute the cum-rates
            these_cumnums = numpy.array([ (len(newsnrs) - bisect.bisect_left(newsnrs, x))/nTrials for x in newsnrs ])
            this_error = numpy.sqrt( these_cumnums / (T_fg * bkg_time) )
            plot_errs = numpy.array([this_error, this_error])
            for n, (y, err) in enumerate(zip(these_cumnums.tolist(),this_error.tolist())):
                if y/T_fg - err <= 0.0:
                    plot_errs[0,n] = y/T_fg - 1e-15
            T_add_bkg = bkg_time
            pylab.errorbar( newsnrs, these_cumnums / T_fg, yerr = plot_errs, color = 'k', ecolor = 'k', marker = 'x', markersize = 6, mew = 1, linestyle = 'None', zorder = 2, label = lbl )

            # print data to text file
            f = open('PlotRatesData-ExtendedBkg.txt', 'w')
            print >> f, txthead % 'extended background estimate (black crosses)'
            print >> f, '#Combined NewSNR\tCumulative Rate (yr^-1)\ty-error (+/- yr^-1)'
            for stat, cumrate, err in zip(newsnrs, (these_cumnums / T_fg).tolist(), this_error.tolist()): 
                print >> f, '%f\t%f\t%f' % (stat, cumrate, err)
            f.close()

            #if lbl != '':
            #    if opts.lin_y:
            #        ytext = yval+.05*(ymax-ymin)
            #    else:
            #        ytext = yval * 10**.25
            #    pylab.annotate( lbl, (xval,yval), xytext=(xval,ytext) )
            # reset xmin/xmax
            xmin = min(min(newsnrs), xmin)
            xmax = max(max(newsnrs), xmax)
            ymin = min(these_cumnums.min(), ymin)
            ymax = max(these_cumnums.max(), ymax)
            # save x, y points for poisson pdf calc.
            added_bkg = [x for x in zip(newsnrs, these_cumnums) if x[0] > bkg_stats.max()] #.append((xval, T_fg / pt_livetime))


        #
        # plot foreground
        #
        loudest_events = []
        if plot_fg:
            match_category.datatype = fg_datatype
            plotvals = data.get_cumrates( match_category, 'stat', rank_by = 'max' ) 
            fg_stats = numpy.array([x[0] for x in plotvals])
            fg_cumrate = numpy.array([y[1] for y in plotvals]) / T_fg 
            # FIXME: add ability to color by ifo colors?
            fclr = 'b'
            if fg_datatype == 'simulation':
                edgclr = 'r'
                edgwth = .5
            elif opts.plot_poisson_pdf:
                edgclr = 'white'
                edgwth = .5
            else:
                edgclr = 'b'
                edgwth = 0.5
            lbl = 'Foreground' #fg_datatype.replace( '_', '-').title()
            pylab.plot( fg_stats, fg_cumrate, marker = '^', markersize = 6.5, markerfacecolor = fclr, markeredgecolor = edgclr, markeredgewidth = edgwth, linestyle = 'None', label = lbl, zorder = 4 )

            # print data to text file
            f = open('PlotRatesData-Foreground.txt', 'w')
            print >> f, txthead % 'foreground coincident events (blue triangles)'
            print >> f, '#Combined NewSNR\tCumulative Rate (yr^-1)'
            for stat, cumrate in zip(fg_stats.tolist(), fg_cumrate.tolist()): 
                print >> f, '%f\t%f' % (stat, cumrate)
            f.close()
            
            # rest xmmin/xmax
            if fg_stats.any():
                xmin = min(min(fg_stats),xmin)
                xmax = max(max(fg_stats),xmax)
            # reset ymin/ymax
            if fg_cumrate.any():
                ymin = min(min(fg_cumrate),ymin)
                ymax = max(max(fg_cumrate),ymax)

            # store the ids of the loudest events
            for x in plotvals[::-1]:
                if x[0] != max(fg_stats):
                    break
                loudest_events.append(x[2])

            #
            # plot foreground triggers that fall in plot-special window as stars
            #
            #if opts.plot_special_time != []:
            #    plot_sp_vals = [x for x in plotvals if data.data_index[x[2]].data.end_time in plot_special_times]
            #    if plot_sp_vals != []:
            #        pylab.scatter( [x[0] for x in plot_sp_vals], [y[1] / T_fg for y in plot_sp_vals], c = 'yellow', marker = (5,1,0), s = 80, linewidth = .5, label = 'Candidate', zorder = 5 )

        #
        # plot extrapolation
        #
        x_ext = numpy.array([])
        y_ext = numpy.array([])
        if plot_fg and opts.extrapolate is not None:
            if opts.verbose:
                print >> sys.stdout, "\textrapolating background..."
            if opts.extrapolate.upper() != 'ERFC':
                # figure out how many points to use for the fit; the default is to start at the first background point that drops below the foreground; this can be overridden with options
                min_x_fit = opts.min_x_fit is None and bkg_hist_stats.min() or opts.min_x_fit
                ### For Gaussian fit:
                    # We're going to fit a Gaussian with 0 mean to the background. Since this is given by y = A*exp(-x^2./(2*sigma^2.))/sqrt(2*pi*sigma^2.), we take the natural log, and thus fit
                    # y = beta[0] + beta[2]*x^2 where beta[0] = ln(A) - ln(2*pi*sigma^2.)/2., beta[2] = -1/(2.*sigma^2.) and y = ln(bkg_cumnum/T_fg)
                    # To do this, we'll do a weighted linear least squares fit using Beta_hat = (transpose(X).W.X)^(-1).transpose(X).W.Y. where W is the nxn matrix of weights (which are just taken
                    # as the error in each background point), X is the nx2 matrix of [1 x[i]^2.], Y is the nx1 matrix of ln(bkg_cumnum[i]/T_fg), and Beta_hat is the 2x1 matrix of fitted betas.
                    # That assumes that the plotted x is linear in whatever the desired statistic is. However, if the plotted x is actually the desired statistic squared, then we should
                    # fit y = beta[0] + beta[2]*x. Generalizing, we actually fit y = beta[0] + beta[2]*x^((2./opts.x_power)), since only the user can know what function of their desired statistic they've
                    # punched into the command line.
                ### For power-law fit:
                    # We're going to fit bkg_cumnum/T_fg = A*statistic^B by doing a least squares fit to y = beta[0] + beta[1]*x where y = ln(bkg_cumnum/T_fg), x = ln(statistic), beta[0] = ln(A) 
                    # and beta[1] = B.
                start_i = pylab.find( bkg_stats >= min_x_fit )[0]
                n = len(bkg_cumnum) - start_i
                X = numpy.matrix(numpy.zeros((n,2), dtype = float))
                Y = numpy.matrix(numpy.zeros((n,1), dtype = float))
                W = numpy.matrix(numpy.zeros((n,n), dtype = float))
                # populate the input matrices
                for i in range(n):
                    Y[i,0] = numpy.log( bkg_cumnum[start_i + i] ) - numpy.log( T_fg )
                    W[i,i] = 1./numpy.power( bkg_error[start_i + i]/bkg_cumnum[start_i + i], 2. )
                    X[i,0] = 1.
                    if opts.extrapolate.upper() == 'POWER':
                        X[i,1] = numpy.log( bkg_stats[start_i + i] )
                    else:
                        X[i,1] = bkg_stats[start_i + i]**(2./opts.x_power)

                # calculate the Beta_hats
                Beta_hat = (X.T * W * X).I * X.T * W * Y

                # Compute the Chisq of the fit
                cum_chisq = ( W * numpy.power((Y - X * Beta_hat), 2.) ).sum()
                # the number of degrees of freedom = Ndof = number of points - number of fitting params - 1
                cum_Ndof = n - 2 - 1

                # now that we have the fit parameters, extrapolate out to the maximum foreground value
                n_ext = 50
                if opts.lin_x:
                    x_ext = numpy.linspace( bkg_stats[start_i], max(fg_stats.max(), bkg_stats.max()), num = n_ext )
                else:
                    x_ext = numpy.logspace( numpy.log10(bkg_stats[start_i]), numpy.log10( max(fg_stats.max(), bkg_stats.max()) ), num = n_ext )
                Xext = numpy.matrix(numpy.zeros((n_ext,2), dtype = float))
                for i,x in enumerate(x_ext):
                    Xext[i,0] = 1.
                    if opts.extrapolate.upper() == 'POWER':
                        Xext[i,1] = numpy.log(x)
                    else:
                        Xext[i,1] = x**(2./opts.x_power)
                y_ext = numpy.exp( numpy.array(Xext * Beta_hat) )


            else:
                # for error function, we use the results of the gaussian fit to the non-cumulative histogram
                n_ext = 50
                sigmasq = -1./(2.*Beta_hat[1,0])

                amp = numpy.exp(Beta_hat[0,0]) * numpy.sqrt(2*numpy.pi*sigmasq) / ( numpy.power(bins[1], 1./opts.x_power) - numpy.power(bins[0], 1./opts.x_power) ) #numpy.power(ds*2, 1./opts.x_power)[0] # (ds* 2.)
                if opts.lin_x:
                    x_ext = numpy.linspace( min_x_fit, fg_stats.max(), num = n_ext )
                else:
                    x_ext = numpy.logspace( numpy.log10(min_x_fit), numpy.log10( fg_stats.max() ), num = n_ext )

                y_ext = amp * special.erfc( numpy.power(x_ext, 1./opts.x_power)/numpy.sqrt(2*sigmasq) )  / (2. * T_fg)
            
            # plot the extrapolated line
            pylab.plot( x_ext, y_ext, 'g--', linewidth = 2, zorder = 5, label = 'Extrapolated Bkg.' )

            # rest xmmin/xmax
            if x_ext.any():
                xmin = min(min(x_ext),xmin)
                xmax = max(max(x_ext),xmax)
            # reset ymin/ymax
            if y_ext.any():
                ymin = min(y_ext.min(),ymin)
                ymax = max(y_ext.max(),ymax)


        #
        # plot poisson pdf
        #
        if opts.plot_poisson_pdf or opts.plot_significance_bands:
            if opts.verbose:
                print "\tcomputing poisson pdf..."
            xrange = 101
            yrange = 2*xrange
            # add added_bkg points to bkg stats
            bkg_stats = numpy.array( bkg_stats.tolist() + [x[0] for x in added_bkg] )
            bkg_cumnum = numpy.array( bkg_cumnum.tolist() + [x[1] for x in added_bkg] )

            # if extrapolating, use the extrpolation beyond the largeest measured background point
            if opts.extrapolate is not None and plot_fg and x_ext.any():
                xcoords = numpy.linspace(bkg_stats.min(), max(x_ext.max(), bkg_stats.max()), xrange)
            else:
                xcoords = numpy.linspace(bkg_stats.min(), bkg_stats.max(), xrange)
            ycoords = numpy.logspace(0, numpy.log10(ymax+0.5), yrange)

            xarray = numpy.zeros((xrange,yrange), dtype=float)
            yarray = numpy.zeros((xrange,yrange), dtype=float)
            zarray = numpy.zeros((xrange,yrange), dtype=float)

            for xi,x in enumerate(xcoords):
                if x <= bkg_stats.max():
                    # lmbda is the cumnum at the closest snr to (the right of) of x in the background
                    lmbda = bkg_cumnum[pylab.find( bkg_stats >= x )[0]]
                # if x greater than bkg_stats.max, must mean we have extrapolated and we are in the region of extrapolation
                elif opts.extrapolate is not None and opts.extrapolate.upper() == "POWER":
                    lmbda = numpy.exp( (Beta_hat[0,0] + Beta_hat[1,0]*numpy.log(x)) ) * T_fg
                elif opts.extrapolate is not None and opts.extrapolate.upper() == "GAUSSIAN":
                    lmbda = numpy.exp( (Beta_hat[0,0] + Beta_hat[1,0]*x**(2./opts.x_power)) ) * T_fg
                elif opts.extrapolate is not None and opts.extrapolate.upper() == "ERFC":
                    lmbda = amp * special.erfc( numpy.power(x,1./opts.x_power)/numpy.sqrt(2*sigmasq) ) / 2. 
                for yi,y in enumerate(ycoords):
                    xarray[xi,yi] += x
                    yarray[xi,yi] += y
                    p = poissonPDFlog10( y, lmbda )
                    zarray[xi, yi] += p
            
            # replace -infs with smallest non-inf value of z - 1
            minz = min([z for z in zarray.flatten() if z != -numpy.inf])
            for xi in range(len(xcoords)):
                for yi in range(len(ycoords)):
                    if zarray[xi,yi] == -numpy.inf:
                        zarray[xi,yi] = minz - 1

            # rescale the y-coordinates by T_fg
            yarray = yarray / T_fg

            # plot it
            if opts.plot_poisson_pdf and 'NaN' not in zarray:
                hxplt = pylab.hexbin(xarray.flatten(), yarray.flatten(), C = zarray.flatten(), gridsize=xrange-1, xscale='linear', yscale='log', edgecolors = 'none', vmin=-10.)
                cb = pylab.colorbar(hxplt, format = pylab.FuncFormatter(ColorFormatter))
                cb.ax.set_ylabel( 'Probability Density' )

            #
            # plot contours
            #
            if opts.plot_significance_bands and opts.max_sigma > 0.:
                if opts.verbose:
                    print "\tplotting contours..."
                sigmas = sorted([numpy.log10(special.erfc((n+1)/numpy.sqrt(2))/2.) for n in range(opts.max_sigma)])
                if opts.plot_poisson_pdf:
                    contplt = pylab.contour(xarray, yarray, zarray, sigmas, colors = 'k', linestyles='dashdot', label='N-\sigma', zorder = 1)
                    cb.add_lines(contplt)
                else:
                    fix_zarray(zarray,sigmas)
                    pylab.contourf(xarray, yarray, zarray, sigmas+[0.0], cmap = pylab.cm.gray_r, alpha = 0.7) 
                    #for n in range(len(sigmas)-1):
                    #    pylab.contourf(xarray, yarray, zarray, [sorted(sigmas)[n], sorted(sigmas)[n+1]], cmap = pylab.cm.gray_r, label = 'N-\sigma', alpha = .2*n+.1)

            # plot a line to show where the extrapolated data begins to be used
            if plot_fg and opts.extrapolate is not None:
                pylab.plot( [bkg_stats.max(), bkg_stats.max()], [fg_cumrate.min(), ymax*10**0.5], 'k-', linewidth = 2, zorder = 3 )


        #
        # finalize plot settings
        #
        pylab.grid()
        pylab.legend(loc = 'lower left', scatterpoints = 1, numpoints = 3)
        pylab.xlabel( stat_label, fontsize=16 )
        # reset xlabel ticks
        for xl in pylab.gca().get_xticklabels():
            xl.set_size(16)
        pylab.ylabel( "Cumulative Rate (%s$^{-1}$)" % opts.time_units, fontsize=16 )
        # reset ylabel ticks
        for yl in pylab.gca().get_yticklabels():
            yl.set_size(16)
        if opts.title:
            if opts.param_name is not None:
                t = 'Cumulative Histogram %s %s' % (param_label.title(), param_parser.param_range_by_group(bkg_category.param_group))#, match_category.veto_cat.replace('_',' ').title())
            else:
                t = 'Cumulative Histogram'
            t = 'Cumulative Histogram of H1L1 triggers with $3.48\mathrm{M_\odot} \le \mathcal{M} < 7.40\mathrm{M_\odot}$'
            pylab.title(t)

        # set axes limits and scale
        if plot_fg:
            xmin = min(fg_stats.min(), bkg_stats.min())
        else:
            xmin = bkg_stats.min()
        if not opts.lin_x:
            pylab.gca().semilogx()
            xmax = xmax*10**(0.3)
        else:
            xmax = xmax+.1*xmax
        if not opts.lin_y:
            pylab.gca().semilogy()
            ymin, ymax = ymin*10**(-0.5), ymax*10**(0.5)
        else:
            ymin, ymax = ymin - .1*ymin, ymax + .1*ymax
        # overrule with input options
        if opts.xmin is not None:
            xmin = opts.xmin
        if opts.xmax is not None:
            xmax = opts.xmax
        if opts.ymin is not None:
            ymin = opts.ymin
        if opts.ymax is not None:
            ymax = opts.ymax

        pylab.xlim( xmin, xmax )
        pylab.ylim( ymin, ymax )
        
        # store plot info
        plot_description = 'cumulative_F%i' % ( len(fnameList) )
        name = InspiralUtils.set_figure_tag( plot_description, datatype_plotted = fg_datatype.upper(), open_box = 'all_data' in fg_datatype or 'exclude_play' in fg_datatype)
        fname = InspiralUtils.set_figure_name(InspiralUtilsOpts, name)
        #fname = re.sub('.png', '.pdf', fname)
        fname_thumb = InspiralUtils.savefig_pylal( filename=fname, dpi=opts.dpi )
        fnameList.append(fname)
        tagList.append(name)

        # get loudest event info for this plot group and add to comments
        comments += 'On instruments: <b>%s</b>   Coinc. ifos: <b>%s</b>   Veto Category: <b>%s</b>' %(''.join(sorted(match_category.on_instruments)), ''.join(sorted(match_category.ifos)), match_category.veto_cat)
        if opts.param_name is not None:
            comments += '   Param-group: <b>%s</b>' %(param_parser.param_range_by_group(match_category.param_group))
        comments += '<br />\n'
        comments += "Background livetime: %.2f %s<br />\n" % (T_bkg, opts.time_units)
        if opts.add_background_file:
            comments+= "Added background livetime: %.2f %s<br />\n" % (T_add_bkg, opts.time_units)
        comments += "<b>Loudest background event%s:</b><br />\n" % (len(loudest_bkg_events) > 1 and 's' or '')
        comments += '<table cellpadding="5", border="1">\n'
        comments += '<tr><th>End Time</th><th>UTC</th><th>ifos</th><th>Total Mass</th><th>%s</th></tr>\n' % stat_label.replace(r'$', '')
        for le_id in loudest_bkg_events:
            le = data.data_index[le_id]
            comments += '<tr><td>%.2f</td><td>%s</td><td>%s</td><td>%.2f</td><td>%.2f</td></tr>' %( le.data.get_value('end_time'), printutils.format_end_time_in_utc(int(le.data.get_value('end_time'))), le.data.get_value('ifos'), le.data.get_value('mtotal'), le.data.get_value('stat') ) 
        comments += '</table>\n'
        if plot_fg:
            comments += "<br />%s livetime: %.3f %s<br />\n" % (fg_datatype.replace('_', '-').title(), T_fg, opts.time_units)
            comments += "<b>Loudest %s event%s:</b><br />\n" % (fg_datatype, len(loudest_events) > 1 and 's' or '')
            comments += '<table cellpadding="5", border="1">\n'
            comments += '<tr><th>End Time</th><th>UTC</th><th>ifos</th><th>Total Mass</th><th>%s</th><th>Measured FAR</th>%s<th>Probability given background</tr>\n' % (stat_label.replace(r'$', ''), opts.extrapolate is not None and '<th>Estimated FAR</th>' or '')
            for le_id in loudest_events:
                le = data.data_index[le_id]
                if opts.add_background_file and le.data.stat >= added_bkg[0][0]:
                    T_far = T_add_bkg
                else:
                    T_far = T_bkg
                far = (len(bkg_stats)-bisect.bisect_left(bkg_stats, le.data.stat))/T_far
                if far == 0:
                    far = '< 1/%.2f %s<sup>-1</sup>' % ( T_far, opts.time_units )
                else:
                    far = '1/%.2e %s<sup>-1</sup>' % (1./far, opts.time_units)
                comments += '<tr><td>%.2f</td><td>%s</td><td>%s</td><td>%.2f</td><td>%.2f</td><td>%s</td>' %( le.data.get_value('end_time'), printutils.format_end_time_in_utc(int(le.data.get_value('end_time'))), le.data.get_value('ifos'), le.data.get_value('mtotal'), le.data.get_value('stat'), far ) 
                if opts.extrapolate is not None:
                    if opts.extrapolate.upper() == "POWER":
                        far_est = numpy.exp( Beta_hat[0,0] + Beta_hat[1,0]*numpy.log(le.data.stat) )
                    elif opts.extrapolate.upper() == "GAUSSIAN":
                        far_est = numpy.exp( Beta_hat[0,0] + Beta_hat[1,0]*le.data.stat**(2./opts.x_power) )
                    else:
                        far_est = amp * special.erfc( numpy.power(le.data.stat, 1./opts.x_power)/numpy.sqrt(2*sigmasq) ) / (2. * T_fg)
                    comments += '<td>1/%.2e %s<sup>-1</sup></td>' % ( 1./far_est, opts.time_units )
                x = le.data.stat
                if x <= bkg_stats.max():
                    lmbda = bkg_cumnum[pylab.find( bkg_stats >= x )[0]]
                    comments += "<td>%.2e</td>" % numpy.power(10, poissonPDFlog10( 1, lmbda ))
                else:
                    lmbda = bkg_cumnum[-1]
                    comments += "<td>< %.2e</td>" % numpy.power(10, poissonPDFlog10( 1, lmbda ))
                comments += '</tr>'
            comments += '</table>\n'
            if opts.extrapolate is not None:
                if opts.extrapolate.upper() == 'POWER':
                    comments += "<br />Fitted <i>y = Ax<sup>B</sup></i><br />Results:<br />"
                    comments += "A = %.2e  B = %.2f<br />" % (numpy.exp(Beta_hat[0,0]), Beta_hat[1,0])
                elif opts.extrapolate.upper() == 'GAUSSIAN':
                    comments += "<br />Fitted <i>y = Aexp(-x<sup>2</sup>/2&#963;<sup>2</sup>)/sqrt(2&#960;&#963;<sup>2</sup>)</i><br />Results:<br />"
                    sigmasq = -1./(2.*Beta_hat[1,0])
                    amp = numpy.exp(Beta_hat[0,0] + numpy.log(2*numpy.pi*sigmasq)/2.)
                    comments += "&#963;<sup>2</sup> = %.2f  A = %.2f<br />" %( sigmasq, amp )
                    comments += "Extrapolated from:<br />&#946;<sub>0</sub> = %.2f  &#946;<sub>1</sub> = %.2f<br />" %( Beta_hat[0,0], Beta_hat[1,0] )
                else:
                    comments += "<br />Fitted <i>y = Aexp(-x<sup>2</sup>/2&#963;<sup>2</sup>)/sqrt(2&#960;&#963;<sup>2</sup>)</i> to non-cumulative plot.<br />Results:<br />"
                    comments += "&#963;<sup>2</sup> = %.2f  A = %.2f<br />" %( sigmasq, amp )
                    comments += "<br />Goodness of fit:<br />"
                    comments += "&chi;<sup>2</sup> = %.2f  Reduced &chi;<sup>2</sup> = %.2f" %( noncum_chisq, noncum_chisq / noncum_Ndof )

        comments += '<hr />\n'
                
    if opts.verbose:
        print >> sys.stdout, "Writing html file and cache."

    # create html of closed box plots
    comment = comments
    plothtml = InspiralUtils.write_html_output( InspiralUtilsOpts, input_args, fnameList, tagList, comment = comments, add_box_flag = True )
    InspiralUtils.write_cache_output( InspiralUtilsOpts, plothtml, fnameList )

    
if __name__ == "__main__":
    main()

