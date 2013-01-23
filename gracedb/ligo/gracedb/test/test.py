
import unittest
import random
import os

from ligo.gracedb.rest import GraceDb

# Test the GraceDb REST API class.
#
#  To run:
#
#     python $PATH_TO_GRACEDB_LIB/test/test.py
# 
#  Environment Variables:
#
#     TEST_SERVICE
#       defaults to https://moe.phys.uwm.edu/gracedb/api/
#       live site would be https://gracedb.ligo.org/api/
#
#     TEST_DATA_DIR
#       defaults to $PATH_TO_GRACEDB_LIB/test/data/
#
#       Files expected:
#
#          burst-cwb.txt
#          cbc-lm2.xml
#          cbc-lm.xml
#          cbc-mbta.gwf
#          upload2.data
#          upload.data
#          upload.data.gz
#
#     X509_USER_PROXY
#
#     X509_USER_CERT
#     X509_USER_KEY


TEST_SERVICE = "https://moe.phys.uwm.edu/gracedb/api/"

class TestGracedb(unittest.TestCase):
    """
        This is as much a test of the REST API on the server
        as it is of the REST client.

        Requires a valid x509 cert in a place that the GraceDb
        class can find it and that the cert DN is accepted by
        the test service.

        This is not a complete set of unit tests,
        but it is a decent start.
    """

    def test_root(self):
        """Did root resource retrieval succeed?"""
        self.assertTrue("CBC" in gracedb.groups)
        pass

    def test_get_events(self):
        """Get the events resource.  Make sure first event looks like an event."""
        events = gracedb.events()
        for event in events:
            self.assertTrue('graceid' in event)
            break

    def test_create_log(self):
        """Create an event log message"""
        message = "Message is {0}".format(random.random())
        resp = gracedb.writeLog(eventId, message)
        self.assertEqual(resp.status, 201)
        new_log_uri = resp.getheader('Location')
        new_log = resp.json()
        self.assertEqual(new_log_uri, new_log['self'])
        check_new_log = gracedb.get(new_log_uri).json()
        self.assertEqual(check_new_log['comment'], message)

    def test_get_log(self):
        """Retrieve event log"""
        logs = gracedb.logs(eventId).json()
        self.assertTrue('numRows' in logs)
        pass

    def test_upload_file(self):
        """Upload and re-upload a file"""

        uploadFile = os.path.join(testdatadir, "upload.data")
        r = gracedb.writeFile(eventId, uploadFile)
        self.assertEqual(r.status, 201) # CREATED
        r_content = r.json()
        link = r_content['permalink']

        self.assertEqual(
                open(uploadFile, 'r').read(),
                gracedb.get(gracedb.files(eventId).json()['upload.data']).read()
                )

        self.assertEqual(
                open(uploadFile, 'r').read(),
                gracedb.get(link).read()
                )

        # Re-upload slightly different file.
        uploadFile2 = os.path.join(testdatadir, "upload2.data")
        r = gracedb.writeFile(
                eventId,
                filename="upload.data",
                filecontents=open(uploadFile2, 'r'))
        self.assertEqual(r.status, 201) # CREATED
        r_content = r.json()
        link2 = r_content['permalink']

        self.assertEqual(
                open(uploadFile2, 'r').read(),
                gracedb.get(gracedb.files(eventId).json()['upload.data']).read()
                )

        self.assertEqual(
                open(uploadFile2, 'r').read(),
                gracedb.get(link2).read()
                )

        self.assertNotEqual(link, link2)


    def test_files(self):
        """Get file info"""
        r = gracedb.files(eventId)
        event = r.json()
        self.assertEqual(r.status, 200)
        self.assertTrue(isinstance(event, dict))

    def test_label_event(self):
        """Label an event"""
        r = gracedb.writeLabel(eventId, "DQV")
        self.assertEqual(r.status, 201) # CREATED
        r = gracedb.labels(eventId, "DQV")
        self.assertEqual(r.status, 200)
        label = r.json()
        self.assertEqual("DQV", label['name'])

    def test_slot_event(self):
        """Create a slot"""
        r = gracedb.createSlot(eventId, "newslot", "event.log")
        self.assertEqual(r.status, 201) # CREATED
        r = gracedb.slot(eventId, "newslot")
        self.assertEqual(r.status, 200)
        slotname = r.json()['value']
        self.assertTrue(slotname.endswith("event.log"))

    def test_create_cwb(self):
        """Create a CWB event"""
        """burst-cwb.txt"""
        eventFile = os.path.join(testdatadir, "burst-cwb.txt")
        r = gracedb.createEvent("Test", "CWB", eventFile)
        self.assertEqual(r.status, 201) # CREATED
        cwb_event = r.json()
        self.assertEqual(cwb_event['group'], "Test")
        self.assertEqual(cwb_event['analysisType'], "CWB")
        self.assertEqual(cwb_event['gpstime'], 1042312876)

    def test_create_lowmass(self):
        """Create a Low Mass event"""
        """cbc-lm.xml"""
        # This is done with the initially created event.
        pass

    def test_create_mbta(self):
        """Create an MBTA event"""
        """cbc-mbta.gwf"""
        eventFile = os.path.join(testdatadir, "cbc-mbta.gwf")
        mbta_event = gracedb.createEvent(
                "Test", "MBTA", eventFile).json()
        self.assertEqual(mbta_event['group'], "Test")
        self.assertEqual(mbta_event['analysisType'], "MBTAOnline")
        self.assertEqual(mbta_event['gpstime'], 1011992635)
        self.assertEqual(mbta_event['far'], 0.000245980441198379)

    def test_replace_event(self):
        graceid = eventId

        old_event = gracedb.event(graceid).json()
        self.assertEqual(old_event['group'], "Test")
        self.assertEqual(old_event['analysisType'], "LowMass")
        self.assertEqual(old_event['gpstime'], 971609248)

        replacementFile = os.path.join(testdatadir, "cbc-lm2.xml")

        response = gracedb.replaceEvent(graceid, replacementFile)
        self.assertEqual(response.status, 202)

        new_event = gracedb.event(graceid).json()
        self.assertEqual(new_event['group'], "Test")
        self.assertEqual(new_event['analysisType'], "LowMass")
        self.assertEqual(new_event['gpstime'], 971609249)

    def test_upload_binary(self):
        """
        Test workaround for Python bug
        http://bugs.python.org/issue11898
        Raises exception if workaround fails.
        """
        uploadFile = os.path.join(testdatadir, "upload.data.gz")
        r = gracedb.writeFile(eventId, uploadFile)
        self.assertEqual(r.status, 201) # CREATED

    def test_unicode_param(self):
        """
        Test workaround for Python bug
        http://bugs.python.org/issue11898
        Raises exception if workaround fails.
        """
        uploadFile = os.path.join(testdatadir, "upload.data.gz")
        r = gracedb.writeFile(eventId, uploadFile)
        self.assertEqual(r.status, 201) # CREATED

    def test_logger(self):
        import logging
        import ligo.gracedb.rest
        import ligo.gracedb.logging
     
        logging.basicConfig()
        log = logging.getLogger('testing')
        log.propagate = False   # Don't write to console

        #gracedb = ligo.gracedb.rest.GraceDb()
        graceid = eventId
     
        log.addHandler(ligo.gracedb.logging.GraceDbLogHandler(gracedb, graceid))

        message = "Message is {0}".format(random.random())
        log.warn(message)

        event_logs = gracedb.logs(graceid).read()
        self.assertTrue(message in event_logs)

if __name__ == "__main__":

    global gracedb, testdatadir, createdEvent, eventId

#   Hacky "global" fixture.
#   Do not want to create a millions events
#   in order to test out event-updating features,
#   which is what a normal test case setUp() would do.

    testdatadir = os.path.join(os.path.dirname(__file__), "data")

    service = os.environ.get('TEST_SERVICE', TEST_SERVICE)
    testdatadir = os.environ.get('TEST_DATA_DIR', testdatadir)

    gracedb = GraceDb(service)
    print "Using service", service

    eventFile = os.path.join(testdatadir, "cbc-lm.xml")
    createdEvent = gracedb.createEvent(
            "Test", "LowMass", eventFile).json()
    eventId = createdEvent["graceid"]

    unittest.main()
