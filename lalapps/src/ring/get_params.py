# Script to read in the names of the injection directories
# from the ini file and write them to a file that is readable
# by Matlab.

import ConfigParser

cp = ConfigParser.ConfigParser()
cp.readfp(open('ihope_ring.ini'))
injdetails=cp.items('injections')

fp=file('runparams.txt','w')
for (dir,seed) in injdetails:
  injType=cp.get(dir,"injection-type")
  print(injType, seed,dir)
  fp.write(injType +"\t" + dir +"\t" +seed+"\n")
fp.close()


