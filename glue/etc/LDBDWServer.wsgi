import DB2
from glue import LDBDWServer

def application(environ, start_response):
    try:
        server = LDBDWServer.Server()
        return server(environ, start_response)
    except Exception, e:
        print >> environ['wsgi.errors'], "%s" % e
        start_response("500 Internal Server Error", [('Content-Type', 'text/plain')])
        msg = "500 Internal Server Error\n\n%s" % e
        return [ msg ]

