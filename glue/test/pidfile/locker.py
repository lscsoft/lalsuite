#!/usr/bin/env python

from glue import pidfile
import time
import sys

pidfile.get_lock("lock.pid")

if len(sys.argv) > 1 and sys.argv[1].isdigit():
    secs=int(sys.argv[1])
else:
    secs=0
time.sleep(secs)

if pidfile.confirm_lock("lock.pid"):
    print "confirmed!"
