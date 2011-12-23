#
# (C) Copyright 2003-2004 Jacek Konieczny <jajcus@jajcus.net>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License Version
# 2.1 as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#

"""PubSub XMPP stanza handling

Normative reference:
   ????
"""

__revision__="$Id$"
__docformat__="restructuredtext en"

import libxml2

## Needed because openfire has a minor bug associated with results
## being ignored in the handshake stage. There is a one line fix
## in this overloaded class.
from ligo.lvalert.lvstanzaprocessor import LVAlertStanzaProcessor


from pyxmpp.stanza import Stanza
from pyxmpp.utils import to_utf8,from_utf8
from pyxmpp.xmlextra import common_ns
from pyxmpp.objects import StanzaPayloadObject
from pyxmpp.iq import Iq
from pyxmpp.xmlextra import get_node_ns_uri
from pyxmpp.exceptions import ClientError, FatalClientError

DISCO_INFO_NS = 'http://jabber.org/protocol/disco#info'
DISCO_ITEMS_NS = 'http://jabber.org/protocol/disco#items'
JABBER_X_DATA_NS = 'jabber:x:data'

PUBSUB_NS = 'http://jabber.org/protocol/pubsub'
PUBSUB_OWNER_NS = 'http://jabber.org/protocol/pubsub#owner'
PUBSUB_NODE_CONFIG_NS = 'http://jabber.org/protocol/pubsub#node_config'

class PubSub(Iq):
    xml_element_name = "pubsub"
    xml_element_namespace = PUBSUB_NS

    def new_pubsub(self,ns_uri=PUBSUB_NS,name="pubsub"):
        ps = Stanza(name)
        ns = ps.xmlnode.newNs(ns_uri,None)
        ps.xmlnode.setNs(ns)
        return ps

    def get_pubsub(self):
	c = self.xmlnode.children
        while c:
            try:
                if c.ns():
                    return c
            except libxml2.treeError:
                pass
            c = c.next
        return None

    def get_pubsub_ns(self):
        """Get a namespace of the stanza payload.

        :return: XML namespace URI of the payload or None if there is no
            payload.
        :returntype: `str`"""
        q=self.get_query()
        if q:
            return get_node_ns_uri(q)
        else:
            return None

    def get_affiliations(self):
        ps = self.new_pubsub()
        ps.add_new_content(None,"affiliations")
        self.add_content(ps.xmlnode)

    def get_nodes(self,ns_uri="http://jabber.org/protocol/disco#items",name="query"):
        self.set_new_content(ns_uri,name)
        self.set_type("get")

    def get_nodes_result(self,cresult):
        c = cresult.xmlnode.children
        c = c.children
        print "List of nodes"
        while c:
            print "Node: %s, %s" % (c.prop("node"), c.prop("name")) 
            c = c.next
        #print "End list of nodes"

    def create_node(self,node_name):
        ps = self.new_pubsub()
        c = ps.add_new_content(None,"create")
        c.newProp("node",node_name)
        ps.add_new_content(None,"configure")
        self.add_content(ps.xmlnode)
        self.set_type("set")

    def generic_result(self,cresult):
        print "Successfully completed operation"

    def create_timeout(self,ctimeout):
        """Process session request time out.  
        :raise FatalClientError:"""
        raise FatalClientError("Timeout while tryin to create node")

    def create_error(self,ctimeout):
        """Process session request time out.  
        :raise FatalClientError:"""
        raise FatalClientError("Error while trying to create node")

    def get_node_config(self,node_name):
        ps = self.new_pubsub(ns_uri=PUBSUB_OWNER_NS)
        c = ps.add_new_content(None,"configure")
        c.newProp("node",node_name)
        self.add_content(ps.xmlnode)
        self.set_type("get")

    def delete_node(self,node_name):
        ps = self.new_pubsub(ns_uri=PUBSUB_OWNER_NS)
        c = ps.add_new_content(None,"delete")
        c.newProp("node",node_name)
        self.add_content(ps.xmlnode)
        self.set_type("set")

    def publish(self,content,node_name):
        ps = self.new_pubsub(ns_uri=PUBSUB_NS)
        c = ps.add_new_content(None,"publish")
        c.newProp("node",node_name)
        c = c.newChild(None,"item",None)
        c = c.newTextChild(None,"entry",content)
        self.add_content(ps.xmlnode)
        self.set_type("set")

    def subscriptions(self,jid):
        ps = self.new_pubsub(ns_uri=PUBSUB_NS)
        c = ps.add_new_content(None,"subscriptions")
        self.add_content(ps.xmlnode)
        self.set_type("set")

    def subscriptions_result(self,cresult):
        c = cresult.xmlnode.children.children.children
        while c:
            print "Node: %s [ %s subid=%s]" % ( c.prop("node"),\
                c.prop("subscription"), c.prop("subid") )
            c = c.next

    def publisher(self,jid,node_name,affiliation):
        ps = self.new_pubsub(ns_uri=PUBSUB_OWNER_NS)
        c = ps.add_new_content(None,"affiliations")
        c.newProp("node",node_name)
        c = c.newChild(None,"affiliation",None)
        c.newProp("jid",jid.bare().as_utf8())
        c.newProp("affiliation",affiliation)
        self.add_content(ps.xmlnode)
        self.set_type("set")

    def affiliations(self,jid,node_name):
        ps = self.new_pubsub(ns_uri=PUBSUB_OWNER_NS)
        c = ps.add_new_content(None,"affiliations")
        c.newProp("node",node_name)
        self.add_content(ps.xmlnode)
        self.set_type("get")

    def affiliations_result(self,cresult):
        c = cresult.xmlnode.children.children.children
        while c:
            print "Affiliation: %s [ %s ]" % ( c.prop("jid"),\
                c.prop("affiliation"))
            c = c.next

    def subscribe(self,jid,node_name):
        ps = self.new_pubsub(ns_uri=PUBSUB_NS)
        c = ps.add_new_content(None,"subscribe")
        c.newProp("node",node_name)
        c.newProp("jid",jid.bare().as_utf8())
        self.add_content(ps.xmlnode)
        self.set_type("set")

    def unsubscribe(self,jid,node_name,node_subid=None):
        ps = self.new_pubsub(ns_uri=PUBSUB_NS)
        c = ps.add_new_content(None,"unsubscribe")
        c.newProp("node",node_name)
        if node_subid:
            c.newProp("subid",node_subid)
        c.newProp("jid",jid.bare().as_utf8())
        self.add_content(ps.xmlnode)
        self.set_type("set")
        self.set_from(jid)

    def subscribe_result(self,cresult):
        node_name = cresult.get_query().children.prop("node")
        print "Successfully subscribed to node %s" % node_name

    def subscribe_timeout(self,ctimeout):
        """Process session request time out.  
        :raise FatalClientError:"""
        raise FatalClientError("Timeout while tryin to subscribe to node")

    def subscribe_error(self,ctimeout):
        """Process session request time out.  
        :raise FatalClientError:"""
        raise FatalClientError("Error while trying to subscribe to node")

    def generic_timeout(self):
        """Process session request time out.  
        :raise FatalClientError:"""
        raise FatalClientError("Timeout while tryin to perform action")

    def generic_error(self,iq):
        """Process session request failure.
        :Parameters:
            - `iq`: IQ error stanza received as result of the session request.
        :Types:
            - `iq`: `pyxmpp.Iq`

        :raise FatalClientError:"""
        err=iq.get_error()
        msg=err.get_message()
        raise FatalClientError("Failed to establish a session: "+msg)

class PubSubMessage(Stanza):
    """Wraper object for <message /> stanzas."""
    stanza_type="message"
    def __init__(self, xmlnode = None, from_jid = None, to_jid = None,\
                stanza_type = None, stanza_id = None, error = None,\
                error_cond = None, stream = None):
        """Initialize a `PubSubMessage` object..........
            - `error_cond`: `unicode`"""

        self.xmlnode=None
        if isinstance(xmlnode,PubSubMessage):
            pass
        elif isinstance(xmlnode,Stanza):
            raise TypeError, "Couldn't make PubSubMessage from other Stanza"
        elif isinstance(xmlnode,libxml2.xmlNode):
            pass
        elif xmlnode is not None:
            raise TypeError, "Couldn't make PubSubMessage from %r" % (type(xmlnode),)

        if xmlnode is None:
            xmlnode="message"

        Stanza.__init__(self, xmlnode, from_jid = from_jid, to_jid = to_jid,\
                    stanza_type = stanza_type, stanza_id = stanza_id, \
                    error = error, error_cond = error_cond, stream = stream)

    def get_event(self):
        """Get the message event

        :return: the message event or `None` if there is no event
        :returntype: `unicode`"""
        n=self.xpath_eval("ns:event")
        if n:
            return from_utf8(n[0].getContent())
        else:
            return None

# vi: sts=4 et sw=4 
