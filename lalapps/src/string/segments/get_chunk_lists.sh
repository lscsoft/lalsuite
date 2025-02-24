#!/usr/bin/env python

import igwn_segments as segments
from igwn_segments import utils as segmentsUtils

S5h1_seg=segmentsUtils.fromsegwizard(open("../S5-H1segments-cat1.txt")).coalesce()
S5h2_seg=segmentsUtils.fromsegwizard(open("../S5-H2segments-cat1.txt")).coalesce()
S5l1_seg=segmentsUtils.fromsegwizard(open("../S5-L1segments-cat1.txt")).coalesce()
S5v1_seg=segmentsUtils.fromsegwizard(open("../S5-V1segments-cat1.txt")).coalesce()
S6h1_seg=segmentsUtils.fromsegwizard(open("../S6-H1segments-cat1.txt")).coalesce()
S6l1_seg=segmentsUtils.fromsegwizard(open("../S6-L1segments-cat1.txt")).coalesce()
S6v1_seg=segmentsUtils.fromsegwizard(open("../S6-V1segments-cat1.txt")).coalesce()

for i, chunk in enumerate(segmentsUtils.fromsegwizard(open("chunks.txt")), 1):
	print "Building segment lists for chunk %s ..." % i
	# make segment into a segmentlist
	chunk = segments.segmentlist([chunk])
	if chunk[0][1] < 920000000:
		segmentsUtils.tosegwizard(open('S5-H1segments-cat1_c%d.txt' % i,'w'), chunk & S5h1_seg)
		segmentsUtils.tosegwizard(open('S5-H2segments-cat1_c%d.txt' % i,'w'), chunk & S5h2_seg)
                segmentsUtils.tosegwizard(open('S5-L1segments-cat1_c%d.txt' % i,'w'), chunk & S5l1_seg)
 		segmentsUtils.tosegwizard(open('S5-V1segments-cat1_c%d.txt' % i,'w'), chunk & S5v1_seg)

		h1_lt = abs(chunk & S5h1_seg)
		h2_lt = abs(chunk & S5h2_seg)
		l1_lt = abs(chunk & S5l1_seg)
		v1_lt = abs(chunk & S5v1_seg)
	else:
		segmentsUtils.tosegwizard(open('S6-H1segments-cat1_c%d.txt' % i,'w'), chunk & S6h1_seg)
		segmentsUtils.tosegwizard(open('S6-L1segments-cat1_c%d.txt' % i,'w'), chunk & S6l1_seg)
		segmentsUtils.tosegwizard(open('S6-V1segments-cat1_c%d.txt' % i,'w'), chunk & S6v1_seg)

		h1_lt = abs(chunk & S6h1_seg)
		h2_lt = 0
		l1_lt = abs(chunk & S6l1_seg)
		v1_lt = abs(chunk & S6v1_seg)

        print "H1 livetime= %d sec" % h1_lt
	print "H2 livetime= %d sec" % h2_lt
        print "L1 livetime= %d sec" % l1_lt
        print "V1 livetime= %d sec" % v1_lt
