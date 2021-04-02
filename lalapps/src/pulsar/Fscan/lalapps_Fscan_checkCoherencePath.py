import glob
import time
import sys

# Check if spec*.txt files exist under the directory given on the command line.
# If these files exist, then the SFT jobs have finished.
# Try again after 2 hr

checkPath = sys.argv[1] + '/spec*.txt'

fileList = glob.glob(checkPath)

if len(fileList) < 1:
   # Try again in 2 hr
   time.sleep(7200)
   fileList = glob.glob(checkPath)

if len(fileList) < 1:
   print('Timout: SFTs jobs under' + sys.argv[1] + 'did not finish.\n')
   exit(1)
else:
   exit(0)

