import pylal.bayespputils as bppu
from optparse import OptionParser
import math
from math import pi
parser=OptionParser()
parser.add_option("-i", "--input",dest="i")
parser.add_option("-o", "--output",dest="o")
parser.add_option("-j","--injFile",dest="j")
parser.add_option("-e","--event",type="int",dest="e")
parser.add_option("-c","--code",dest="c")

(options, args) = parser.parse_args()

if (options.i == None or options.o == None or options.j == None or options.e == None):
    parser.error("not enough arguments")
input=options.i
output=options.o
injFile=options.j
event=options.e
code=options.c

file = open(input,'r')

from pylal import SimInspiralUtils
injections = SimInspiralUtils.ReadSimInspiralFromFiles([injFile])

if(len(injections)<event):
    print "Error: You asked for event %d, but %s contains only %d injections" %(event,injfile,len(injections))
    sys.exit(1)
else:
    injection = injections[event]


def grabParam(posterior,injVal,param_name):
    counter = 0
    total = 0
    countSamples = 0
    
    column = None
    for line in posterior:
        words = line.split()
        if counter ==0:
            counter = 1
            for word in range(len(words)):
                if words[word] == param_name:
                    column = word
        elif column == None:
            file.seek(0)
            return None
        else:
            if float(words[column]) < injVal:
                countSamples += 1
            total += 1
    file.seek(0)
    return float(countSamples)/total

def grabMass(posterior,injMc,injEta):
    counter = 0
    total = 0
    countSamples1 = 0
    countSamples2 = 0

    columnM1 = None
    columnM2 = None
    for line in posterior:
        words = line.split()
        if counter == 0:
            counter = 1
            for word in range(len(words)):
                if words[word] == 'm1':
                    columnM1 = word
                if words[word] == 'm2':
                    columnM2 = word
        elif columnM1 == None or columnM2 ==None:
            print 'failed to find m1 and m2'
            file.seek(0)
            return None,None
        else:
            m1 = float(words[columnM1])
            m2 = float(words[columnM2])
            mchirp,eta = bppu.ms2mc(m1,m2) 
            if mchirp < injMc:
                countSamples1 += 1
            if eta < injEta:
                countSamples2 += 1
            total += 1
    file.seek(0)
    return float(countSamples1)/total,float(countSamples2)/total

def grabMassq(posterior,injMc,injEta):
    counter = 0
    total = 0
    countSamples1 = 0
    countSamples2 = 0

    columnM1 = None
    columnM2 = None
    for line in posterior:
        words = line.split()
        if counter == 0:
            counter = 1
            for word in range(len(words)):
                if words[word] == 'mc' or words[word] == 'chirpmass' or words[word] == 'mchirp':
                    columnM1 = word
                if words[word] == 'q' or words[word] == 'asym_massratio':
                    columnM2 = word
        elif columnM1 == None or columnM2 ==None:
            print 'failed'
            file.seek(0)
            return None,None
        else:
            mc = float(words[columnM1])
            q = float(words[columnM2])
            eta = bppu.q2eta(mc,q)
            if mc < injMc:
                countSamples1 += 1
            if eta < injEta:
                countSamples2 += 1
            total += 1
    file.seek(0)
    return float(countSamples1)/total,float(countSamples2)/total

if code == 'm':
    param_names = {'mchirp':'mchirp','eta':'eta','dec':'dec','ra':'ra','polarisation':'psi','dist':'dist','time':'time','phase':'phi0','iota':'iota'}
elif code == 'mct' or code == 'mcf':
    param_names = {'mchirp':'mc','eta':'eta','dec':'dec','ra':'ra','polarisation':'psi','dist':'dist','time':'time','phase':'phi_orb','iota':'iota'}
elif code == 'i':
    param_names = {'mchirp':'chirpmass','eta':'massratio','dec':'declination','ra':'rightascension','polarisation':'polarisation','dist':'distance','time':'time','phase':'phase','iota':'inclination','logprior':'logPrior','logL':'logL'}

if code != None:
    if code == 'm':
        mchirpVal,etaVal = grabMass(file,injection.mchirp,injection.eta)
    elif code == 'mct' or code == 'mcf' or code == 'i':
        mchirpVal,etaVal = grabMassq(file,injection.mchirp,injection.eta)
    else:
        mchirpVal = grabParam(file,injection.mchirp,param_names['mchirp'])
        etaVal = grabParam(file,injection.eta,param_names['eta'])
    decVal = grabParam(file,injection.latitude,param_names['dec'])
    raVal = grabParam(file,injection.longitude,param_names['ra'])
    polarisation = grabParam(file,math.fmod(injection.polarization,pi),param_names['polarisation'])
    timeval = grabParam(file,injection.geocent_end_time+1e-9*injection.geocent_end_time_ns,param_names['time'])
    distval = grabParam(file,injection.distance,param_names['dist'])
    iotaVal = grabParam(file,injection.inclination,param_names['iota'])
    phaseVal = grabParam(file,injection.coa_phase,param_names['phase'])

outputFile = open(output + '1d_' + code + str(event),'w')
outputFile.write('mchirp '+str(mchirpVal)+' \n')
outputFile.write('eta '+str(etaVal)+' \n')
outputFile.write('ra '+str(raVal)+' \n')
outputFile.write('dec '+str(decVal)+' \n')
outputFile.write('pol '+str(polarisation)+' \n')
outputFile.write('time '+str(timeval)+' \n')
outputFile.write('dist '+str(distval)+' \n')
outputFile.write('iota '+str(iotaVal)+' \n')
outputFile.write('phase '+str(phaseVal)+' \n')
