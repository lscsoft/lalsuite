from xml.etree import ElementTree as ET
import copy
from math import exp,log
import random
import sys
from pylal import bayespputils
from pylal.bayespputils import vo_nest2pos

xmlns='http://www.ivoa.net/xml/VOTable/v1.1'

USAGE = """%prog [options] <results.xml>
Produce posterior samples from nested samples
stored in results.xml as a VOTable.

WARNING: Will add a new table to results.xml and overwrite
unless output file specified. All existing tables are preserved.
"""

try:
    register_namespace=ET.register_namespace
except AttributeError:
    def register_namespace(prefix,uri):
        ET._namespace_map[uri]=prefix
register_namespace('vot',xmlns)

if __name__ == '__main__':
    from optparse import OptionParser
    parser=OptionParser(USAGE)
    parser.add_option("-o","--out",help="Output to FILE, will not over-write results.xml",metavar="FILE")
    parser.add_option("-N","--nlive",type="int",help="Over-ride the number of live points in results.xml",metavar="NUM")
    (opts,args)=parser.parse_args()
    
    if len(args)==1:
      infile=args[0]
    else:
      print USAGE
      sys.exit(1)
    
    tree = ET.ElementTree()
    tree.parse(infile)
    nsresource = [node for node in tree.findall('{%s}RESOURCE'%(xmlns)) if node.get('name')=='Nested sampling run'][0]
    postable = vo_nest2pos(nsresource,opts.nlive)
    nsresource.append(postable)
    if opts.out:
      tree.write(opts.out)
    else:
      tree.write(infile)
