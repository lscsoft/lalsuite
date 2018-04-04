#!/usr/bin/python
import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
import pylal.bayespputils as bppu
from optparse import OptionParser
from math import pi
matplotlib.rcdefaults()


parser=OptionParser()
parser.add_option("-i", "--input",dest="i")
parser.add_option("-o", "--output",dest="o")
parser.add_option("-j","--injFile",dest="j")
parser.add_option("-e","--event",type="int",dest="e")
parser.add_option("-u","--url",action="store_true",dest="u",default=False)
parser.add_option("-n","--name",dest="n")
parser.add_option("-p","--password",dest="p")
parser.add_option("--samples_per_bin",type="int",dest="spb", default=32)
parser.add_option("--plot",action="store_true",dest="x",default=False)
(options, args) = parser.parse_args()

if (options.i == None or options.o == None or options.j == None or options.e == None):
    parser.error("not enough arguments")
input=options.i
output=options.o
injFile=options.j
event=options.e
urlUsed=options.u
userName = options.n
userPassword = options.p
makePlot = options.x
samplesPerBin = options.spb

###############                                                                                                                            
def open_url_wget(url,folder,un=userName,pw=userPassword,eventNum=0, args=[]):
    import subprocess
    import urlparse
    name = folder+"/posteriors/posterior_"+str(eventNum)+".dat"
    if un is not None and pw is not None:
        args+=["-O",name,"--user",un,"--password",pw,"--no-check-certificate"]
    retcode=subprocess.call(['wget']+[url]+args)

    return retcode

def plot_kdtree(tiles,injCoords = None):
    myfig = plt.figure()  
    ax = plt.axes()
    ax.set_xlim((0,2*pi))
    ax.set_ylim((-pi/2,pi/2))
    plt.xlabel('RA (rad)')
    plt.ylabel('dec (rad)')

    for tile in tiles:
        patches = []

        art = mpatches.Rectangle((tile[0],tile[2]),tile[1] - tile[0],tile[3] - tile[2],fill = True, facecolor=matplotlib.cm.jet(1.-tile[4]),linewidth=0)#,edgecolor = 'k') 
        ax.add_patch(art)
    if injCoords is not None:
        plt.scatter(injCoords[0],injCoords[1],marker=(5,1),s=250,c='k')
    myfig.savefig(output+'/plot'+str(event))
###############load posterior#########                                                                                                     

data=bppu.PEOutputParser('common')
if(urlUsed):
    open_url_wget(input,output,eventNum=event)
    inputFileObj = open(output+'/posteriors/posterior_'+str(event)+'.dat')
else:
    inputFileObj = open(input)
dataObj0 = data.parse(inputFileObj)

                                                                                                                          
from pylal import SimInspiralUtils
injections = SimInspiralUtils.ReadSimInspiralFromFiles([injFile])
if(len(injections)<event):
    print "Error: You asked for event %d, but %s contains only %d injections" %(event,injfile,len(injections))
    sys.exit(1)
else:
    injection = injections[event]

posterior = bppu.Posterior(dataObj0,injection)
outFile = open(output+'/kdresult' + str(event), 'w')
outFile.write('label injection_cl injection_area\n')
confidenceLevels = [0.1,0.2,0.3,0.4,0.5,0.6,0.68,0.7,0.8,0.9]

######################################

def mapping(ra,dec):
    return (ra,dec)

def isEdge(bounds):
    if bounds[0][0] == 0. or bounds[1][0]==2*pi or bounds[0][1] == -pi/2. or bounds[1][1]==pi/2:
        return True
    return False


#set up evrything for running kd algorithm
if 'ra' and 'dec' not in posterior.names:
    if 'rightascension' and 'declination' in posterior.names:
        new_param_names = ['ra','dec']
        post_names = ['rightascension','declination']
        posterior.append_multiD_mapping(new_param_names, mapping , post_names)
    else:
        print 'param not found'
print 'out test --------'
injCoordinates=[posterior['ra'].injval,posterior['dec'].injval]

#create kd node list
nodeList, areas, injInfo = bppu.kdtree_bin2Step(posterior,['ra','dec'],confidenceLevels,samples_per_bin = samplesPerBin,skyCoords=True,injCoords = injCoordinates)
print 'details'
print injInfo
print areas

rollingCL = 0
edgeCL = 100.
totalSamples = 0.
for node in nodeList:
    totalSamples+=node[2]

for node in nodeList:
    rollingCL += node[2]
    if isEdge(node[0]):
        if float(rollingCL)/totalSamples < edgeCL:
            edgeCL = float(rollingCL)/totalSamples
            break

for cl in confidenceLevels:
    temp_area = areas[cl]*(180/pi)**2.
    if cl > edgeCL:
        outFile.write('kd_areaCL' +str(cl) + ' ' + str(temp_area) + ' EDGE\n')
    else:
        outFile.write('kd_areaCL' +str(cl) + ' ' + str(temp_area) + '\n')

outFile.write('kd_injCL ' + str(injInfo[3])+' \n')                                                                                                      
temp_area = injInfo[4]*(180/pi)**2.
outFile.write('kd_injCLarea ' + str(temp_area) + '\n')

if makePlot:
    tiles = []
    temp = 0
    total = 0
    for node in nodeList:
        total +=  node[2]
    for node in nodeList:
        temp += node[2]
        tiles.append([node[0][0][0],node[0][1][0],node[0][0][1],node[0][1][1],float(temp)/total])
    plot_kdtree(tiles,injCoordinates)
