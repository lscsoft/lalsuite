"""
This module is intended to make it easier to build web clients
written in Python using the urllib2 module that can
interoperate with the @LIGO.ORG infrastructure.
"""

import re
import urllib2
import exceptions

class LIGOSAMLClientException(exceptions.Exception):
    """
    """
    pass

try:
    import kerberos
except ImportError as e:
    msg = """\n
The glue.auth.saml module requires the python-kerberos Python module to
be installed. On both Debian and RedHat-based systems like CentOS and
Scientific Linux the name of the package to install is 'python-kerberos'.
"""
    raise LIGOSAMLClientException(msg)

class HTTPNegotiateAuthHandler(urllib2.BaseHandler):
    """
    This class uses an existing Kerberos ticket to authenticate
    via HTTP Negotiate Authentication. An instance of this class
    can be passed into the build_opener function from the urllib2
    module.

    Modified from source found at

    http://selenic.com/pipermail/mercurial/2008-June/019776.html
    """

    rx = re.compile('(?:.*,)*\s*Negotiate\s*([^,]*),?', re.I)
    handler_order = 480  # before Digest auth

    def __init__(self, service_principal):
        """
        service_principal is the Kerberos principal of the
        host against which the client authenticates. It 
        should usually be the string 'HTTP@login.ligo.org'.
        """
        self.retried = 0
        self.context = None
        self.service_principal = service_principal

    def negotiate_value(self, headers):
        authreq = headers.get('www-authenticate', None)

        if authreq:
            mo = HTTPNegotiateAuthHandler.rx.search(authreq)
            if mo:
                return mo.group(1)

        return None

    def generate_request_header(self, req, headers):
        neg_value = self.negotiate_value(headers)
        if neg_value is None:
            self.retried = 0
            return None

        if self.retried > 5:
            raise urllib2.HTTPError(req.get_full_url(), 401, "negotiate auth failed", headers, None)

        self.retried += 1

        result, self.context = kerberos.authGSSClientInit(self.service_principal)

        if result < 1:
            return None

        result = kerberos.authGSSClientStep(self.context, neg_value)

        if result < 0:
            return None

        response = kerberos.authGSSClientResponse(self.context)
        
        return "Negotiate %s" % response

    def authenticate_server(self, headers):
        neg_value = self.negotiate_value(headers)
        if neg_value is None:
            return None

        if kerberos.authGSSClientStep(self.context, neg_value) < 1:
            pass

    def clean_context(self):
        if self.context is not None:
            kerberos.authGSSClientClean(self.context)

    def http_error_401(self, req, fp, code, msg, headers):
        try:
            neg_hdr = self.generate_request_header(req, headers)

            if neg_hdr is None:
                return None

            req.add_unredirected_header('Authorization', neg_hdr)
            resp = self.parent.open(req)

            self.authenticate_server(resp.info())

            return resp
        
        finally:
            self.clean_context()

