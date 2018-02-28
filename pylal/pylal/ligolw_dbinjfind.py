import sys

from glue.ligolw import lsctables
from glue.ligolw import ilwd

from pylal import ligolw_sqlutils as sqlutils
from pylal import ligolw_dataUtils as dataUtils

def make_rec_sngls_table( connection, recovery_table ):
    """
    Makes a temporary table containing events from the given recovery table
    that could potentially be injections --- i.e., events from the "simulation"
    datatype --- and the process id of the injection jobs that created them.
    This allows for quick matching between injections and single events later
    on.
    """
    sqlquery = ''.join(['''
        CREATE TEMP TABLE rec_sngls AS
            SELECT
                experiment_summary.sim_proc_id AS sim_proc_id,
            ''', recovery_table, '''.*
        FROM
            ''', recovery_table, '''
        ''', sqlutils.join_experiment_tables_to_sngl_table( recovery_table ), '''
        WHERE
            experiment_summary.datatype == "simulation"''' ])
    connection.cursor().execute(sqlquery)


def dbinjfind( connection, simulation_table, recovery_table, match_criteria, rough_match = None, rejection_criteria = [], rough_rejection = None, verbose = False ):

    # validate simulation_table and recovery_table
    simulation_table = sqlutils.validate_option( simulation_table )
    recovery_table = sqlutils.validate_option( recovery_table )

    # create DataRow classes to store data for each table
    simColumns = sqlutils.get_column_names_from_table( connection, simulation_table )
    recColumns = sqlutils.get_column_names_from_table( connection, recovery_table )
    
    SimDataRow = dataUtils.createDataRowClass( simulation_table, columns = simColumns )
    RecDataRow = dataUtils.createDataRowClass( recovery_table, columns = recColumns )

    # create a temporary table to store the eligible foreground events that can be matched
    if verbose:
        print >> sys.stdout, "Getting eligible events..."
    make_rec_sngls_table( connection, recovery_table )

    # if using rough match, create an index for it
    rough_match_test = ''
    if rough_match is not None:
        simRough, recRough, winRough = rough_match
        simRough = sqlutils.validate_option( simRough )
        recRough = sqlutils.validate_option( recRough )
        sqlquery = "CREATE INDEX rs_rmtch_idx ON rec_sngls (%s)" % recRough
        connection.cursor().execute( sqlquery )
        rough_match_test = "rec_sngls.%s >= sim.%s - %f AND rec_sngls.%s <= sim.%s + %f AND\n" %( recRough, simRough, winRough, recRough, simRough, winRough ) 


    #
    # Remove triggers that match all_data triggers
    #
    if rejection_criteria != []:
        if verbose:
            print >> sys.stdout, "Applying rejection test to eligible events..."
        # If comparing to all_data for rejection, create a temp table of all data events
        # This rejection test only uses the single-ifo triggers from coincident events
        sqlquery = ''.join(['''
            CREATE TEMP TABLE all_data_sngls AS
                SELECT
                    ''', recovery_table, '''.*
                FROM
                ''', recovery_table, '''
                ''', sqlutils.join_experiment_tables_to_sngl_table( recovery_table ), '''
                WHERE
                    experiment_summary.datatype == "all_data"''' ])
        connection.cursor().execute(sqlquery)

        rough_test = ''
        if rough_rejection is not None:
            rejRough, rejRoughWin = rough_rejection
            rejRough = sqlutils.validate_option( rejRough )
            sqlquery = "CREATE INDEX ads_rmtch_idx ON all_data_sngls (%s)" % rejRough
            connection.cursor().execute( sqlquery )
            rough_test = "all_data_sngls.%s >= rec_sngls.%s - %f AND all_data_sngls.%s <= rec_sngls.%s + %f AND\n" % ( rejRough, rejRough, rejRoughWin, rejRough, rejRough, rejRoughWin )

        # cycle over the rejection criteria, creating a function in the database for each
        rejection_tests = []
        for n,(thisFunc, window) in enumerate(rejection_criteria):
            compF = dataUtils.CompareDataRows(RecDataRow, RecDataRow)
            funcName = 'matches_all_data%i' % n
            # Note: setting the match criteria also sets the needed columns
            compF.set_matchCriteriaA(thisFunc)
            compF.set_matchCriteriaB(thisFunc)
            # need different diff function if using eThinca
            if thisFunc == 'eThinca':
                diffFunc = compF.eThincaSngl
            else:
                diffFunc = compF.diffRowARowB
            compF.create_dbCompF(connection, diffFunc, funcName, window)
            simSnglCols = ','.join(['rec_sngls.%s' %(col) for col in compF.neededColumnsA])
            allSnglCols = ','.join(['all_data_sngls.%s' %(col) for col in compF.neededColumnsB])
            rejection_tests.append( '%s(%s, %s)' %(funcName, simSnglCols, allSnglCols) ) 

        # now remove triggers
        sqlquery = ''.join([ '''
            DELETE FROM
                rec_sngls
            WHERE EXISTS (
                SELECT
                    *
                FROM
                    all_data_sngls
                WHERE
                    ''', rough_test, '\nAND '.join( rejection_tests ), ')' ])
        connection.cursor().execute(sqlquery)
        connection.commit()

    #
    #   Determine Sim-Sngl matches
    #
    
    if verbose:
        print >> sys.stdout, "Applying match criteria to find sim-sngl maps..."
    # cycle over the match criteria, creating a function in the database for each
    match_tests = []
    for n,(simFunc, snglFunc, window) in enumerate(match_criteria):
        compF = dataUtils.CompareDataRows(SimDataRow, RecDataRow)
        # set the name of the compare function to use in the database
        funcName = 'are_match%i' % n
        compF.set_matchCriteriaA(simFunc)
        compF.set_matchCriteriaB(snglFunc)
        # need different diff function if using eThinca
        if simFunc == 'eThinca':
            diffFunc = compF.eThincaSim
        else:
            diffFunc = compF.diffSimSngl
        compF.create_dbCompF(connection, diffFunc, funcName, window)
        simCols = ','.join(['sim.%s'%(col) for col in compF.neededColumnsA])
        snglCols = ','.join(['rec_sngls.%s'%(col) for col in compF.neededColumnsB])
        match_tests.append( '%s(%s, %s)' %(funcName, simCols, snglCols) )

    # determine matches
    sqlquery = ''.join(["""
    CREATE TEMP TABLE found_inj AS
    SELECT
        sim.simulation_id AS sim_id,
        rec_sngls.event_id AS event_id
    FROM
        """, simulation_table, """ AS sim, rec_sngls
    WHERE
        sim.process_id == rec_sngls.sim_proc_id AND
        """, rough_match_test, '\n\t\tAND'.join( match_tests ) ])
    connection.cursor().execute(sqlquery)
    connection.commit()

def strlst_is_subset(stringA, stringB):
    return set(stringA.split(',')).issubset(set(stringB.split(',')))

def write_coincidences(connection, map_label, search, process_id, verbose = False):
    """
    Writes coincidences to coinc_event_map table.
    """
    # for all the maps, see if there is another coincidence
    if verbose:
        print >> sys.stdout, "Getting mapped sngls belonging to a coincident event..."
    connection.create_aggregate("ag_cat", 1, sqlutils.aggregate_concatenate)
    connection.create_function("issubset", 2, strlst_is_subset)

    sqlquery = '''
    CREATE INDEX finj_simid_idx ON found_inj (sim_id);
    CREATE INDEX finj_eid_idx ON found_inj (event_id);
    CREATE TEMP TABLE coinc_inj AS
        SELECT
            found_inj.sim_id AS sid,
            found_inj.event_id AS evid,
            coinc_event_map.coinc_event_id AS ceid
        FROM
            found_inj
        JOIN
            coinc_event_map
        ON (
            coinc_event_map.event_id == evid )
        WHERE issubset(
                (
                SELECT ag_cat(c.event_id)
                FROM coinc_event_map AS c
                WHERE c.coinc_event_id == ceid
                GROUP BY c.coinc_event_id
                ORDER BY c.event_id ASC),
                (
                SELECT ag_cat(b.event_id)
                FROM found_inj AS b
                WHERE b.sim_id == sid
                GROUP BY b.sim_id
                ORDER BY b.event_id ASC)
                );
    CREATE INDEX cij_eid_idx ON coinc_inj (evid);
    '''
    connection.cursor().executescript(sqlquery)
    # get the sim_coincs
    sqlquery = "SELECT DISTINCT sid, ceid FROM coinc_inj"
    sim_coincs = [(ilwd.ilwdchar(sim_id), ilwd.ilwdchar(ceid)) for ceid, sim_id in connection.cursor().execute( sqlquery ).fetchall()]

    # get the sim_sngls
    sqlquery = "SELECT sim_id, event_id FROM found_inj WHERE event_id NOT IN (SELECT DISTINCT evid FROM coinc_inj)"
    sim_sngls = [(ilwd.ilwdchar(sim_id), ilwd.ilwdchar(eid)) for sim_id, eid in connection.cursor().execute( sqlquery ).fetchall()]

    # create a new coinc_def id for this map label, if it already doesn't exist
    coinc_def_id = sqlutils.write_newstyle_coinc_def_entry( connection, map_label, search=search )

    # get the time_slide id
    # XXX: NOTE: We are assuming that all simulation entries have the same time_slide id
    sqlquery = 'SELECT DISTINCT time_slide_id FROM experiment_summary WHERE datatype LIKE "simulation%"'
    time_slide_id = connection.cursor().execute(sqlquery).fetchall()
    if len(time_slide_id) > 1:
        raise ValueError, "more than one time_slide_id found for the simulation datatype"
    elif len(time_slide_id) == 0:
        raise ValueError, "no time_slide_id found for the simulation datatype"
    time_slide_id = ilwd.ilwdchar(time_slide_id.pop()[0])
    
    # write the number of new entries needed for the sim_coincs to the coinc_event table
    if verbose:
        print >> sys.stdout, "Adding injection maps to coinc_event table..."
    new_ceids = sqlutils.add_coinc_event_entries( connection, process_id, coinc_def_id, time_slide_id, len(sim_coincs) ) 

    # add these new entries to coinc_event_map table
    if verbose:
        print >> sys.stdout, "Adding injection-coinc_event maps to coinc_event_map table..."
    sqlquery = 'INSERT INTO coinc_event_map (coinc_event_id, table_name, event_id) VALUES (?,?,?)'
    connection.cursor().executemany( sqlquery, [(str(ceid), sim_id.table_name, str(sim_id)) for ceid, (sim_id, _) in zip(new_ceids, sim_coincs)] )
    connection.cursor().executemany( sqlquery, [(str(ceid), coinc_ceid.table_name, str(coinc_ceid)) for ceid, (_, coinc_ceid) in zip(new_ceids, sim_coincs)] )

    # ditto for the sim-sngls
    if verbose:
        print >> sys.stdout, "Adding injection-sngl maps to coinc_event_map table..."
    new_ceids = sqlutils.add_coinc_event_entries( connection, process_id, coinc_def_id, time_slide_id, len(sim_sngls) ) 
    connection.cursor().executemany( sqlquery, [(str(ceid), sim_id.table_name, str(sim_id)) for ceid, (sim_id, _) in zip(new_ceids, sim_sngls)] )
    connection.cursor().executemany( sqlquery, [(str(ceid), eid.table_name, str(eid)) for ceid, (_, eid) in zip(new_ceids, sim_sngls)] )

    # update the number of events in the coinc_event table
    if verbose:
        print >> sys.stdout, "Updating coinc_event nevents column..."
    sqlutils.update_coinctab_nevents( connection )
