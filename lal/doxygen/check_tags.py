# check_tags.py - check Doxygen tags

__author__ = 'Karl Wette <karl.wette@ligo.org>'
__copyright__ = 'Copyright (C) 2015 Karl Wette'

import sys
from xml.etree.cElementTree import ElementTree

# print error message and exit
def fail(msg):
    sys.stderr.write('%s: %s\n' % (sys.argv[0], msg))
    sys.exit(1)

# get input arguments
tagfile, = sys.argv[1:]

# parse tag file
tree = ElementTree()
try:
    tree.parse(tagfile)
except:
    fail("could not parse XML input from '%s'" % tagfile)

# find duplicate documentation anchors
anchors = dict()
for elem in tree.iter('docanchor'):
    if not elem.text in anchors:
        anchors[elem.text] = []
    anchors[elem.text].append(elem)
dup_anchors = dict((anchor, elems) for anchor, elems in anchors.items() if len(elems) > 1)

# remove duplicate documentation anchors, preferring namespace and group anchors to other kinds
parent_map = dict((elem, parent_elem) for parent_elem in tree.iter() for elem in parent_elem)
for anchor in dup_anchors:
    elem_rank = dict()
    for elem in dup_anchors[anchor]:
        parent_elem = elem
        while parent_elem in parent_map:
            parent_elem = parent_map[parent_elem]
            if parent_elem.tag == "compound":
                parent_elem_kind = parent_elem.get("kind", "")
                if parent_elem_kind == "namespace":
                    elem_rank[elem] = 0
                elif parent_elem_kind == "group":
                    elem_rank[elem] = 1
                else:
                    elem_rank[elem] = 10
    dup_anchors[anchor].sort(key=lambda elem: elem_rank[elem])
    dup_anchors[anchor].pop(0)
    for elem in dup_anchors[anchor]:
        parent_map[elem].remove(elem)

# write tag file
tree.write(tagfile, encoding="UTF-8", xml_declaration=True)
