import os

from glue.segmentdb import segmentdb_utils
from glue.segmentdb import query_engine

from glue.ligolw import dbtables
from glue.ligolw.utils import ligolw_sqlite

from glue.segments import segment, segmentlist


def test_basic_seg_summary(engine):
    res = segmentdb_utils.query_segments( engine, 'segment_summary', [ ('H1','DMT-TESTSEG_1',1,924900000,924900016,0,0) ] )

    if res != [ segmentlist( [segment(924900000, 924900016)] ) ]:
        return False

    res = segmentdb_utils.query_segments( engine, 'segment_summary', [ ('H1','DMT-TESTSEG_1',2,924900000,924900016,0,0) ] )
    
    if res != [ segmentlist( [segment(924900008, 924900010)] ) ]:
        return False

    return True


def test_expand_versions(engine):
    res = segmentdb_utils.expand_version_number( engine, ('H1','DMT-TESTSEG_1','*',924900000,924900016,0,0) )

    if len(res) != 3:
        return False

    values = [('H1', 'DMT-TESTSEG_1', 2, 924900008, 924900010, 0, 0), 
              ('H1', 'DMT-TESTSEG_1', 1, 924900000, 924900008, 0, 0), 
              ('H1', 'DMT-TESTSEG_1', 1, 924900010, 924900016, 0, 0)]

    for v in values:
        found = False
        for r in res:
            found = found or v == r

        if not found:
            return False

    return True


def test_optimized_query(engine):
    res = segmentdb_utils.query_segments( engine, 'segment_summary', [ ('H1','DMT-TESTSEG_2',1,924900000,924900016,0,0),
                                                                      ('H1','DMT-TESTSEG_3',1,924900000,924900016,0,0) ] )

    if res[0] != segmentlist([segment(924900000, 924900010)]): 
        return False

    if res[1] != segmentlist([segment(924900008, 924900016)]):
        return False

    return True


if __name__ == '__main__':
    db_name    = 'seg_test_db.sqlite'
    target     = dbtables.get_connection_filename(db_name, None, True, False)
    connection = ligolw_sqlite.setup(target)

    engine     = query_engine.SqliteQueryEngine(connection)

    ligolw_sqlite.insert(connection, ['test_segdb_utils.xml'])

    print "Testing basic segment summary...",
    print test_basic_seg_summary(engine) and "suceeded." or "FAILED."

    print "Testing expanding version numbers...",
    print test_expand_versions(engine) and "suceeded." or "FAILED."

    print "Testing optimized segment query...",
    print test_optimized_query(engine) and "suceeded." or "FAILED."


    connection.close()
    os.remove(db_name)

