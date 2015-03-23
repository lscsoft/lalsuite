import nose
from ligo.gracedb.rest import GraceDb, HTTPError
from nose.tools import assert_true, assert_equal
import sys, os
import json
import voeventparse
import StringIO

# XXX Might be better to use option parsing for these things
# No logging unless user specifies --log-file
TEST_SERVICE = "https://moe.phys.uwm.edu/branson/api/"
DATAFILE = "cbc-lm.xml"
LOGFILE = "/tmp/voevent_test.log"

service = os.environ.get('TEST_SERVICE', TEST_SERVICE)
testdatadir = os.path.join(os.path.dirname(__file__), "data")
testdatadir = os.environ.get('TEST_DATA_DIR', testdatadir)
datafile = os.path.join(testdatadir, DATAFILE)

# Module level global variables
g = GraceDb(service)

# Variables to which we will need access at the module level
event = None
graceid = None
update_voevent = None
retraction_voevent = None
preliminary_voevent = None
preliminary_voevent_text = None

f = open(LOGFILE, "w")

# Utility for getting out a dictionary of ivorns and citation types
def get_citations_dict(v):
    citations_dict = {}
    for e in v.Citations.iterchildren():
        f.write("Got tag, value: %s, %s" % (e.tag, e.text))
        if e.tag == 'EventIVORN':
            ivorn = e.text
            citation_type = e.attrib['cite']
            citations_dict[ivorn] = citation_type
    return citations_dict

def setup_module():
    global graceid
    r = g.createEvent("Test", "gstlal", datafile, "LowMass")
    event = r.json()
    graceid = event['graceid']
    f.write("created event %s\n" % graceid)
    # Upload fake skymap file to use later
    # XXX May want some more error handling.
    r = g.writeLog(graceid, "Fake skymap file.", filename = "fake_skymap.txt",
        filecontents = "Fake skymap.", tagname = "sky_loc")
    r = g.writeLog(graceid, "Fake skymap image file.", filename = "fake_skymap_image.txt",
        filecontents = "Fake skymap image.", tagname = "sky_loc")
    f.write("successfully wrote log\n")

def teardown_module():
    global f
#    pass
    f.close()

def test_create_preliminary_voevent():
    global preliminary_voevent
    global preliminary_voevent_text
    f.write("inside test prelim, graceid  = %s\n" % graceid)
    try:
        r = g.createVOEvent(graceid, "Preliminary")
        rdict = r.json()
        assert_true('voevent_type' in rdict.keys())
        f.write('got text = %s\n' % rdict['text'])
        preliminary_voevent_text = rdict['text']
        preliminary_voevent = voeventparse.load(StringIO.StringIO(rdict['text']))
    except HTTPError, e:
        outfile = open('tmp.html', 'w')
        outfile.write(str(e))
        outfile.close()

def test_retrieve_voevent():
    r = g.voevents(graceid)
    voevent_list = r.json()['voevents']
    voevent_list = [v['text'] for v in voevent_list]
    assert_true(len(voevent_list) == 1 and preliminary_voevent_text in voevent_list)

def test_create_update_voevent():
    global update_voevent
    r = g.createVOEvent(graceid, "Update", skymap_filename = "fake_skymap.txt",
        skymap_type = "FAKE", skymap_image_filename = "fake_skymap_image.txt")
    rdict = r.json()
    f.write("got update text = %s\n"  % rdict['text'])
    assert_true('voevent_type' in rdict.keys())
    update_voevent = voeventparse.load(StringIO.StringIO(rdict['text']))

def test_ivorns_unique():
    preliminary_ivorn = preliminary_voevent.attrib['ivorn']
    f.write("preliminary ivorn = %s\n" % preliminary_ivorn)
    update_ivorn = update_voevent.attrib['ivorn']
    f.write("update ivorn = %s\n" % update_ivorn)
    assert_true(update_ivorn != preliminary_ivorn)

def test_citation_section():
    update_citations = get_citations_dict(update_voevent)
    preliminary_ivorn = preliminary_voevent.attrib['ivorn']
    assert_equal(update_citations[preliminary_ivorn], 'supersedes')

def test_create_retraction_voevent():
    global retraction_voevent
    r = g.createVOEvent(graceid, "Retraction")
    rdict = r.json()
    assert_true('voevent_type' in rdict.keys())
    f.write("got retraction text = %s" % rdict['text'])
    retraction_voevent = voeventparse.load(StringIO.StringIO(rdict['text']))

def test_retraction_citations():
    # Parse retraction voevent and check for correct citations
    retraction_citations = get_citations_dict(retraction_voevent)
    preliminary_ivorn = preliminary_voevent.attrib['ivorn']
    update_ivorn = update_voevent.attrib['ivorn']
    cond = retraction_citations[preliminary_ivorn] == 'retraction' 
    cond = cond and retraction_citations[update_ivorn] == 'retraction'
    assert_true(cond)

nose.runmodule()

