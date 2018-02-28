#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

'''
A collection of utilities to assist in storing and sorting data queried from a database
or xml document.
'''

import sys, re, math
import time, datetime

from glue.ligolw.utils import print_tables
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue import git_version

from pylal.xlal.date import XLALGPSToUTC
try:
    from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
except ImportError:
    # s6 code
    from pylal.xlal.date import LIGOTimeGPS
from pylal import tools
from pylal import ligolw_sqlutils as sqlutils


__author__ = "Collin Capano <collin.capano@ligo.org>"
__version__ = git_version.id


# =============================================================================
#
#                           Utilities
#
# =============================================================================

#
#
#   Tools for manipulating statistics
#
#

def get_row_stat(row, arg):
    """
    Method to evaluate the desired operation on columns from a row in a table.
    The desired operation can be either a pre-defined function (if it exists in
    the row's namespace) or a function of the elements in the row's name space.
    Syntax for the arg is python.  The name space available to eval is limited
    to the columns in the row and functions in the math module.
    """
    # speed up: if the arg exists in the row, just return it
    try:
        val = getattr(row, arg)
        try:
            # we'll try returning val as a function
            return val()
        except TypeError:
            # that failed, so val must not be a function; just return it
            return val
    except AttributeError:
        # arg isn't a simple argument of row, so we'll have to eval it
        try:
            row_dict = row.__dict__
        except AttributeError:
            row_dict = dict([ [name, getattr(row,name)] for name in dir(row) ])
        safe_dict = {}
        safe_dict.update(row_dict)
        safe_dict.update(math.__dict__)
        return eval(arg, {"__builtins__": None}, safe_dict)

def get_needed_columns(column_list, function):
    """
    Returns a list of columns from the given column list that are needed by
    the given function string. This can be used to reduce the number of columns
    that are passed to get_row_stat.

    Parameters
    ----------
    column_list: list
        A list of strings given the possible columns to pull out.
    function: str
        A string specifying the match criteria. Can either be:
            * "endTime" or "startTime": in this case, all columns in
              column_list with "end_time", "start_time", or "ifo" in their name
              will be retrieved.
            * "eThinca": in this case, all columns with mass1, mass2, mchirp,
              eta, tau[0-9], time, or [Gg]amma[0-9] in their name will be
              retrieved.
            * a python string that is a function of one or more of the columns
              in column_list.

    Returns
    -------
    needed_columns: str
        The subset of columns needed by the function.
    """
    if function == 'eThinca':
        regex = re.compile(
            'ifo|mass1|mass2|mchirp|eta|time|tau[0-9]|[Gg]amma[0-9]')
        needed_cols = [col for col in column_list \
            if regex.search(col) is not None]
    elif function == 'endTime' or function == 'startTime':
        # if endTime, we'll need all columns with "end_time" in it; this
        # will get both an end_time column from a sngl_inspiral table and
        # a {site}_end_time from a sim_inspiral table; ditto start time
        regex = re.compile('ifo|end_time|start_time')
        needed_cols = [col for col in column_list \
            if regex.search(col) is not None]
    else:
        needed_cols = [col for col in column_list \
            if re.search('(%s)' % col, function) is not None]
    return needed_cols


def createDataRowClass( classTableName, baseClass = None, columns = [] ):
    """
    Creates a DataRow class. If classTableName is the same as a table in
    lsctables, and baseClass is not specified, the DataRow class will inherit
    from that table's RowType. If baseClass is specified, the DataRow class
    that is created will inherit from that base class. If baseClass is not
    specified and classTableName is not in lsctables, then DataRow will just
    inherit from object. Regardless of inheritance, the DataRow will have a
    __slots__ attribute. Any columns that are specified in columns will be
    added to the __slots__ class. The DataRow will also have a tableName
    attribute. This will be whatever classTableName is, regardless of what
    baseClass is set to. Think of the DataRow class as a more arbitrary version
    of an lsctable row.
    
    @classTableName: a string specifying the DataRow's tableName. If baseClass
    is not specified, and tableName is the same as a table in lsctables, the
    DataRow will inherit from that table's RowType. Example: 'sngl_inspiral'
    
    @baseClass: specify what class the DataRow should inherit from. Example:
    lsctables.SnglInspiral

    @columns: a list of strings specifying columns to add to the DataRow's
    __slots__ attribute. All columns in __slots__ will also be an attribute of
    the DataRow class. Only columns not in the base class's __slots__ attribute
    will be added.
    
    Note that this function returns a class, not an instance of a class.
    """
    # determine the base class
    if baseClass is not None:
        base = baseClass
    elif classTableName in lsctables.TableByName:
        base = lsctables.TableByName[ classTableName ].RowType
    else:
        base = object

    # define the class
    class DataRow( base ):

        tableName = classTableName
        if '__slots__' not in dir( base ):
            __slots__ = columns
        else:
            __slots__ = base.__slots__ + [c for c in set(columns)-set(base.__slots__)]
        
        def __init__(self):
            """
            If all slots are not populated, we will get an AttributeError when
            using get_value. To avoid this, initialize all slots as None.
            """
            for column in self.__slots__:
                setattr(self, column, None)

        def store(self, dataTuple):
            """
            Takes a list of tuples of (column_name, data) and assigns the
            values to the object's variables. The column_name must be in 
            self.__slots__.
            @dataTuple: a list of tuples in which the first element is the
             column name and the second is the value to assign.
            """
            for col, val in dataTuple:
                setattr( self, col, val )

        def get_value(self, arg):
            """
            Returns the result of some operation on the elements in self.
            @arg: can be the name of any defined function in self's base class,
             a slot in self, or a function of either or both. See get_row_stat
             for more info.

            Example:
            >>> from glue.ligolw import lsctables
            >>> SnglInspRow = lsctables.createDataRowClass('sngl_inspiral')
            >>> test = SnglInspRow()
            >>> test.store([('snr', 6.), ('chisq', 32.), ('chisq_dof', 16.)])
            >>> test.get_value('snr**2.')
            36.0
            >>> test.get_value('get_new_snr')
            5.8993671171391338
            >>> test.get_value('log(get_new_snr())')
            1.7748450768765174
            """
            return get_row_stat( self, arg )

    return DataRow


def combineRowStats( function, rows ):
    """
    Performs the desired function on the list of single statistics. Note: this
    can only combine one statistic from each row.

    @function: can be either a known pre-set (see below) or an arbitrary
    function. If an arbitrary function, it must be in terms of the ifo names.
    
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

    # otherwise, evaulate the function explicitly
    safe_dict = dict([ [name,val] for name,val in rows.items() + math.__dict__.items() if not name.startswith('__') ])

    try:
        return eval( function, {"__builtins__":None}, safe_dict )
    except NameError:
        # this can happen if an ifo that's specified in the combining function is not in the coincident ifos; in this case, just return None
        return None


def createCombineRowsMethod( tableName, columns, functionList ):
    """
    Creates a CombineRows class that can be used in a sqlite database to
    combine rows on the fly. Takes in a sngl_function, which is the function
    used to combine columns within a single row, and a combining_function,
    which is the function used to combine the results of the sngl_functions
    across rows.

    @tableName: the name of the table that will be reading from. If it is a table in lsctables.py, all methods and columns from that table will be inherited.
    @columns: the list of columns that will be storing data to. This list must be in the same order that will be reading data in from the database with.
    @functionList: a list of tuples. The first item should be the combining function to use, in terms of the ifos to combine, and the second item should be the sngl function to use, in terms of columns or methods of the sngl_row.
    """

    sngl_row = createDataRowClass(tableName, columns)

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


#
#
#   Utilities for storing data
#
#
class dbDataRow:
    """
    A class to assist in loading data from and performing functions on tables
    in a SQLite database.
    """
    def __init__(self, connection, tableName, baseClass = None ):
        self.connection = connection
        self.columns = sqlutils.get_column_names_from_table( self.connection, tableName )
        self.rowClass = createDataRowClass( self.tableName, baseClass, self.columns )
        self._function = None

    def set_function( self, function ):
        self._function = function

    def store( self, *rowData):
        thisRow = self.rowClass()
        thisRow.store( zip( self.columns, rowData ) )
        return thisRow

    def get_db_value( self, *rowData ):
        thisRow = self.store( rowData )
        return thisRow.get_value( self._function )

    def create_db_func( self, function, funcName ):
        self.set_function( function )
        self.connection.create_function( funcName, len(self.columns), self.get_db_value ) 


class CompareDataRows:
    """
    A class that can perform various types of comparison tests between
    arbitrary DataRow classes. The class has the following attributes:
    *classA: A DataRow class. Note: this is a class, not an instance
      of a class.
    *classB: A DataRow class. This can be the same type of class as
      classA, or different. Like classA, this is a class, not an instance
      of that class.
    *matchCriteriaA: What column, or function of columns, to get from
      classA when doing a comparison between an instance of classA and
      an instance of classB.
    *matchCriteriaB: What column, or function of columns, to get from
      classB when doing a comparison between an instance of classA and
      an instance of classB.
    *_diffFunc: What function to perform to differentiate classA from classB.
      This function should be one of the functions below; it takes data to
      populate an instance of classA and an instance of classB, and returns a
      numerical value >= 0 representing the difference between these instances of
      classA and classB. This value can then be compared to the window size to
      determine if A and B are the same or not.
    *window: How large of a window to use to consider an instance of
     classA equal to an instance of classB.

    Example:
     >>> classA = createDataRowClass( 'sngl_inspiral' )
     >>> classB = createDataRowClass( 'sngl_inspiral' )
     >>> compF = CompareDataRows( classA, classB )
     >>> compF.set_diffFunc( compF.diffRowARowB )
     >>> compF.set_window( 0.1 )
     >>> compF.set_matchCriteriaA('mass1/mass2')
     >>> compF.set_matchCriteriaB = ('mass1/mass2')
     >>> dataA = [('mass1', '10.0'), ('mass2', '5.0')]
     >>> dataB = [('mass1', '10.1'), ('mass2', '5.0')]
     >>> compF.compare( dataA, dataB )
     True
     >>> compF.set_window(0)
     >>> compF.compare( dataA, dataB )
     False
    """
    def __init__(self, RowClassA = None, RowClassB = None):
        self.classA = RowClassA
        self.classB = RowClassB
        self._matchCriteriaA = None
        self._matchCriteriaB = None
        self._neededColumnsA = None
        self._neededColumnsB = None
        self.diffFunc = None
        self.window = None

    def set_classA(self, DataRowClass):
        self.classA = DataRowClass

    def set_classB(self, DataRowClass):
        self.classB = DataRowClass

    def set_matchCriteriaA(self, match_criteria):
        """
        Sets the match criteria for classA. Also sets the columns needed for
        the given match criteria.
        """
        self._matchCriteriaA = match_criteria
        # set the needed columns for the given match criteria
        self.set_neededColumnsA()

    @property
    def matchCriteriaA(self):
        return self._matchCriteriaA

    def set_matchCriteriaB(self, match_criteria):
        """
        Sets the match criteria for classB. Also sets the columns needed for
        the given match criteria.
        """
        self._matchCriteriaB = match_criteria
        # set the needed columns for the given match criteria
        self.set_neededColumnsB()

    @property
    def matchCriteriaB(self):
        return self._matchCriteriaB

    def get_needed_columnsAB(self, AorB):
        """
        Retrieves which columns in the desired class is needed for the match
        criteria.

        Parameters
        ----------
        AorB: str
            Either 'A' or 'B'; which class to get the columns for.

        Returns
        -------
        needed_cols: list
            The list of needed columns; see get_needed_columns for
            details.
        """
        return get_needed_columns(
                getattr(self, 'class%s' % AorB).__slots__,
                getattr(self, 'matchCriteria%s' %AorB))

    def set_neededColumnsA(self):
        self._neededColumnsA = self.get_needed_columnsAB('A')

    @property
    def neededColumnsA(self):
        return self._neededColumnsA

    def set_neededColumnsB(self):
        self._neededColumnsB = self.get_needed_columnsAB('B')

    @property
    def neededColumnsB(self):
        return self._neededColumnsB

    def set_diffFunc( self, function ):
        self.diffFunc = function

    def set_window(self, window):
        self.window = window
    #
    #   Functions
    #
    def _diff( self, a, b ):
        """
        Returns the absolute value of the difference between a and b.

        Parameters
        ----------
        a: float or integer
        b: float or integer

        Returns
        -------
        difference: float or integer
            The abs difference between a and b.
        """
        return abs(a - b)

    def compare(self, a, b):
        """
        Runs self.diffFunc on a and b and checks that that is <= self.window.

        Parameters
        ----------
        a: instance of classA row
            The data passed to the first argument of self.diffFunc.
        b: instance of classB row
            The data passed to the second argument of self.diffFunc.

        Returns
        -------
        comparison: bool
            True if self.diffFunc(a, b) is <= self.window; False otherwise.
        """
        return self.diffFunc(a, b) <= self.window

    def dbWrapper(self, *args):
        """
        A database wrapper for the compare functions.

        Parameters
        ----------
        args: list
            A list of values. The first len(self.neededColumnsA) is assumed to
            be the data for classA, in the order that neededColumnsA is in.
            The rest of the values are assumed to be the data for classB, in
            the order that neededColumnsB is in.

        Returns
        -------
        comparison: bool
            The result of self.compare, where the first argument passed is
            the data from classA and the second is data from classB.
        """
        dataA = [args[i] for i in range(len(self.neededColumnsA))]
        dataB = [args[i] for i in range(len(self.neededColumnsA), len(args))]
        dataA = zip(self.neededColumnsA, dataA)
        dataB = zip(self.neededColumnsB, dataB)
        return self.compare(dataA, dataB)
        
    def create_dbCompF(self, connection, diffFunc, compFuncName, window):
        """
        Creates a function in the given connection to a database that allows
        the given diffFunc to be performed on classA and classB on the fly.
        The matchCriteria and the neededColumns for each class must be already
        set (this should happen simultaneously by using set_matchCriteria(A|B).

        Parameters
        ----------
        connection: sqlite3.connection
            A connection to SQLite database.
        diffFunc: function
            The function to use to do comparisons; must be one of the
            functions defined in this class.
        compFuncName: str
            What to call the call function in the database; must be unique.
        window: float
            The size of the window to use when determining whether or not
            classA and classB are the same.
        """
        if self._matchCriteriaA is None:
            raise ValueError("matchCriteriaA not set! " +\
                "Run self.set_matchCriteriaA with appropriate arguments.")
        if self._neededColumnsA is None:
            raise ValueError("neededColumnsA not set! " +\
                "Run self.set_matchCriteriaA to set the needed columns and " +\
                "the match criteria.")
        if self._matchCriteriaB is None:
            raise ValueError("matchCriteriaB not set! " +\
                "Run self.set_matchCriteriaB with appropriate arguments.")
        if self._neededColumnsB is None:
            raise ValueError("neededColumnsB not set! " +\
                "Run self.set_matchCriteriaB to set the needed columns and " +\
                "the match criteria.")
        self.set_diffFunc(diffFunc)
        self.set_window(window)
        connection.create_function(compFuncName,
            len(self.neededColumnsA)+len(self.neededColumnsB), self.dbWrapper)

    def diffRowARowB(self, dataA, dataB):
        """

        Runs self.diff on self.classA and self.classB using self.matchCriteriA
        and self.matchCriteriaB. A or B can be any DataRow class; the only
        requirement is that their match criteria (set by
        self.matchCriteria(A|B)) be a function of their slots. Special match
        criteria are 'startTime' and 'endTime'. In this case,
        (start|end)_time+1e-9*(start|end)_time_ns will calculated.

        Parameters
        ----------
        dataA: list
            A list of tuples with data to populate this instance of classA.
            The first value of each tuple is the column name, the second the
            value, e.g., ('ifo', 'H1').
        dataB: list
            A list of data tuples to populate this instance of classB.

        Returns
        -------
        diff: float
            The return of self._diff(a,b), where a(b) is the matchCritieraA(B)
            function run on dataA(B).
        """
        # store the data
        rowA = self.classA()
        rowA.store(dataA)
        rowB = self.classB()
        rowB.store(dataB)
        # analyze the sngl functions
        if self.matchCriteriaA == 'startTime':
            a = rowA.start_time + 1e-9*rowA.start_time_ns
        elif self.matchCriteriaA == 'endTime':
            a = rowA.end_time + 1e-9*rowA.end_time_ns
        else:
            a = rowA.get_value( self.matchCriteriaA )
        if self.matchCriteriaB == 'startTime':
            b = rowB.start_time + 1e-9*rowB.start_time_ns
        elif self.matchCriteriaB == 'endTime':
            b = rowB.end_time + 1e-9*rowB.end_time_ns
        else:
            b = rowB.get_value( self.matchCriteriaB )
        return self._diff(a, b)

    def diffSimSngl( self, simData, snglData ):
        """
        Same as diffRowARowB, except that classA is assumed to be some sort of
        simulation table (e.g., sim_inspiral) and classB is assumed to be some
        sort of single-IFO table (e.g., sngl_inspiral). This assumption only
        matters if 'startTime' or 'endTime' are the match criteria for classA.
        In that case, the observatory that recorded the event in classB is
        retrieved from classB.ifo. This is then used to pick out the
        appropriate end|start time to use from classA. For example, if H1 is
        the ifo in the snglData, then
        h_(end|start)_time+1e-9*h_(end|start)_time_ns will be retrieved from
        the simData.
        @simData: a list of tuples with data to populate this instance of
         classA. If self.matchCriteriaA is 'endTime' or 'startTime', classA is
         assumed to be a row in a simulation table, and must have
         {site}_(start|end)_time(_ns) columns.
        @snglData: a list of tuples with data to populate this instance of
         classB. If self.matchCriteriaB is 'endTime' or 'startTime', classB is
         assumed to be a rown in a single-IFO table, and must have an ifo
         column.
        """
        # store the data
        simRow = self.classA()
        simRow.store(simData)
        snglRow = self.classB()
        snglRow.store(snglData)
        # analyze the sim function
        if self.matchCriteriaA == 'startTime':
            site = snglRow.ifo.lower()[0]
            a = getattr( simRow, '%s_start_time' % site ) + 1e-9*getattr( simRow, '%s_start_time_ns' % site )
        elif self.matchCriteriaA == 'endTime':
            site = snglRow.ifo.lower()[0]
            a = getattr( simRow, '%s_end_time' % site ) + 1e-9*getattr( simRow, '%s_end_time_ns' % site )
        else:
            a = simRow.get_value( self.matchCriteriaA )
        # analyze the sngl function
        if self.matchCriteriaB == 'startTime':
            b = snglRow.start_time + 1e-9*snglRow.start_time_ns
        elif self.matchCriteriaB == 'endTime':
            b = snglRow.end_time + 1e-9*snglRow.end_time_ns
        else:
            b = snglRow.get_value( self.matchCriteriaB )
        return self._diff(a, b)

    def eThincaSim( self, simData, snglData):
        """
        Computes the eThinca distance between an instance of self.classA and an
        instance of self.classB. This assumes that classA inherited from the
        SimInspiral class and classB inherited from the SnglInspiral class.
        @simData: List of data tuples (column_name, value) with which to
         populate this instance of self.classA.
        @snglData: List of data tuples (column_name, value) with which to
         populate this instance of self.classB.
        """
        simRow = self.classA()
        simRow.store(simData)
        snglRow = self.classB()
        snglRow.store(snglData)
        # lal expects the event_id and simulation_id to be integers
        simRow.simulation_id = 0
        snglRow.event_id = 0
        return tools.XLALEThincaParameterForInjection( simRow, snglRow )

    def eThincaSngl( self, snglDataA, snglDataB ):
        """
        Computes the eThinca distance between an instance of self.classA and an
        instance of self.classB. This assumes that both classA and classB
        inherited from the SnglInspiral class.
        @snglDataA: List of data tuples (column_name, value) with which to
         populate this instance of self.classA.
        @snglDataB: List of data tuples (column_name, value) with which to
         populate this instance of self.classB.
        """
        snglRowA = self.classA()
        snglRowA.store(snglDataA)
        snglRowB = self.classB()
        snglRowB.store(snglDataB)
        # lal expects the event_ids to be integers
        snglRowA.event_id = 0
        snglRowB.event_id = 0
        try:
            ethincaVal = tools.XLALCalculateEThincaParameter( snglRowA, snglRowB )
        except ValueError:
            # not coincident, just return inf
            ethincaVal = float('inf')
        return ethincaVal 


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

    def __init__(self, offset_vector = {}, datatype = None, veto_cat = None, on_instruments = frozenset(['ALL']), ifos = frozenset(['ALL']), param_group = None):
        self.offset_vector = OffsetVector(offset_vector)
        self.datatype = datatype
        self.veto_cat = veto_cat
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
        @data: an instance of the DataRow class listing statistics and methods associated with this element
        """
        def __init__(self, thisid, data):
            self._id = thisid
            self.data = data
            self.cumrates = {}

        def update(self, _id = None, data = None):
            # update id
            if _id is not None:
                self._id = _id
            # update data
            if data is not None:
                self.data = data

    def __init__(self):
        """
        A list of all the data elements is kept as an index.
        """
        self.data_index = {}

    def add_data(self, _id, categories, data):
        """
        Adds a new DataElement to self.

        @_id: some unique value to identify the data element
        @categories: a list of categories that this data falls in. If one or more of these categories are equal (equality determined by the default Category match_criteria) to a category already in all_categories, the category is set to that category. This results in distinct categories only being saved once in memory, with all DataElements that share that category pointing to the same memory address.
        """
        d = self.DataElement( _id, data )
        self.data_index[d._id] = d
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
            self.data_index[_id].update( data = data)
            self.refresh_categories( [self.data_index[_id]] )

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
        collapsedRow = createDataRowClass( 'collapsedRow' )
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
    combinedRow = createDataRowClass( 'combinedRow' )
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



