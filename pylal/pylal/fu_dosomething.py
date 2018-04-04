from pylal.fu_utils import *

##############################################################################
# Examples of follow up functions
##############################################################################
def dosomething(trig):
  # all functions should first initialize their HTML container
  # in the following way this sets up the appropriate output 
  # files and directories
  container = HTMLcontainer(trig,__prog__)
  container.txt = "your favorite module output"
  linkFile = open(container.locallink,'w')
  writeIULHeader(linkFile)
  linkFile.write("<b>Put your cool stuff here<b>")
  linkFile.close()
  # Your function body should go here and it should interface with
  # the container class to determine where output goes 
  # at the end make sure and return the container
  return container
