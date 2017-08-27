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

# warn on duplicate documentation anchors
anchors = dict()
for elem in tree.iter('docanchor'):
    if elem.text in anchors:
        anchors[elem.text] = anchors[elem.text] + 1
    else:
        anchors[elem.text] = 1
for anchor in anchors:
    if anchors[anchor] > 1:
        print('%s: warning: duplicate anchor %s' % (tagfile, anchor))
