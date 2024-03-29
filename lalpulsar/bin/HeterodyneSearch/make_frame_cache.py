"""
Matt Pitkin - 13/11/06
code to take scan a directory of frame files and output it as a frame cache file between given
times
"""

from __future__ import print_function

# import modules
import sys
import os
import getopt

# program usage


def usage():
    msg = """\
Usage: make_frame_cache [options]

  -h, --help              display this message
  -d, --frame-dir         directory of frame file
  -s, --gps-start-time    start time of frames to output
  -e, --gps-end-time      end time fo frames to output
  -o, --output-file       file to output frame cache to
"""
    print(msg, file=sys.stderr)


# parse command line
shortop = "hd:s:e:o:"
longop = ["help", "frame-dir=", "gps-start-time=", "gps-end-time=", "output-file="]

try:
    opts, args = getopt.getopt(sys.argv[1:], shortop, longop)
except getopt.GetoptError:
    usage()
    sys.exit(1)


# process options
for o, a in opts:
    if o in ("-h", "--help"):
        usage()
        sys.exit(1)
    elif o in ("-d", "--frame-dir"):
        frame_dir = a
    elif o in ("-s", "--gps-start-time"):
        start = int(a)
    elif o in ("-e", "--gps-end-time"):
        end = int(a)
    elif o in ("-o", "--output-file"):
        output = a
    else:
        print("Unknown option: {}".format(o), file=sys.stderr)
        usage()
        sys.exit(1)

# get all files in frame dir
try:
    files = os.listdir(frame_dir)
    print(files[0], file=sys.stderr)
except Exception as e:
    print("Problem listing directory {}".format(frame_dir), file=sys.stderr)
    sys.exit(1)

# open output file
try:
    f = open(output, "w")
except Exception as e:
    print("Can't open file {}".format(output), file=sys.stderr)
    sys.exit(1)

files.sort()

# scan through directory and output frame files with given time range
i = 0
j = 0
while i < len(files):
    # check if frame file
    if ".gwf" in files[i]:
        # get frame channel, time and duration
        ifo = files[i][0]
        frinfo = files[i].split("-")
        channel = frinfo[1]  # channel should be first field
        time = int(frinfo[2])  # time should be the second field

        # find the - before the duration (last value)
        index1 = files[i].rfind("-")
        # find the . before gwf
        index2 = files[i].rfind(".")

        duration = int(files[i][index1 + 1 : index2])

        # check if file is within time range and output
        if time + duration > start and time <= end:
            cache = (
                ifo
                + " "
                + channel
                + " "
                + str(time)
                + " "
                + str(duration)
                + " "
                + "file://localhost"
                + frame_dir
                + "/"
                + files[i]
                + "\n"
            )
            f.write(cache)
            j += 1

    i += 1

f.close()

# if no frame files were found then say so
if j == 0:
    print("No frames files between {} and {}.".format(start, end), file=sys.stderr)
    os.remove(output)

sys.exit(0)
