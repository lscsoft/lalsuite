import os

for f in os.listdir('.'):
  if f.split('.')[-1] == 'gz':
    newf = f.replace('.gz','')
    os.rename(f,newf)
    print "renaming " + f + " to " + newf
  
