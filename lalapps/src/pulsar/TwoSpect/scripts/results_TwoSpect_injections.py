#!/usr/bin/python

import argparse, sys, re

parser = argparse.ArgumentParser(description='Script that collects the output files into one file')
parser.add_argument('--dir', required=True, type=str, help='Base path to output directories')
parser.add_argument('--jobs', required=True, type=int, help='Number of jobs (also the number of directories)')
parser.add_argument('--IFO', required=True, type=str, help='Interferometers tested (specify multiple as CSV)')
parser.add_argument('--ULoutput', type=str, help='Filename for UL file')
parser.add_argument('--candidateOutput', type=str, help='Filename for candidates file')
args = parser.parse_args()

if args.ULoutput: uloutput = open('{}/{}'.format(args.dir, args.ULoutput), 'w')
if args.candidateOutput: candidates = open('{}/{}'.format(args.dir, args.candidateOutput), 'w')

ifoval = 0
IFOList = args.IFO.split(',')
for ii in range(0, len(IFOList)):
    if IFOList[ii]=='H1': ifoval = ifoval | (1 << 0)
    elif IFOList[ii]=='L1': ifoval = ifoval | (1 << 1)
    elif IFOList[ii]=='V1': ifoval = ifoval | (1 << 2)
    elif IFOList[ii]=='H2' or IFOList[ii]=='H2r': ifoval = ifoval | (1 << 3)
    else: sys.exit('Invalid IFO value!\n')

for ii in range(0, args.jobs):
    errfile = open('{}/err/err.{}'.format(args.dir, ii),'r')
    err = errfile.readlines()
    if not re.match('^Program finished', err[len(err)-1]):
        errfile.close()
        continue
    errfile.close()
    
    injectionfile = open('{}/{}/injections.dat'.format(args.dir, ii),'r')
    injections = injectionfile.readlines()
    numberinjections = len(injections)
    if numberinjections<10:
        injectionfile.close()
        continue

    reversedinjections = injections[::-1]
    numberinjections = 0
    reversedrecentinjections = []
    for injection in reversedinjections:
        reversedrecentinjections.append(injection)
        numberinjections += 1
        if numberinjections==10: break
    if numberinjections<10: continue

    recentinjections = reversedrecentinjections[::-1]

    if len(recentinjections)!=len(set(recentinjections)): continue

    for jj in range(0, 10):
        recentinjections[jj] = recentinjections[jj].rstrip()

        if args.ULoutput:
            ulfile = open('{}/{}/uls_{}.dat'.format(args.dir, ii, jj),'r')
            uls = ulfile.readlines()
            uloutput.write(uls)
            ulfile.close

        if args.candidateOutput:
            logfile = open('{}/{}/logfile_{}.txt'.format(args.dir, ii, jj),'r')
            log = logfile.readlines()
            log = log[::-1]
            kk = 1
            foundoutlier = 0
            while (re.match('^fsig = (\d+.\d+), period = (\d+.\d+), df = (\d+.\d+), RA = (\d+.\d+), DEC = (-?\d+.\d+), R = (\d+.\d+), h0 = (\d+.\d+e-\d+), Prob = (-\d+.\d+), TF norm = (\d+.\d+e\+\d+)', log[kk])):
                candidate = re.match('^fsig = (\d+.\d+), period = (\d+.\d+), df = (\d+.\d+), RA = (\d+.\d+), DEC = (-?\d+.\d+), R = (\d+.\d+), h0 = (\d+.\d+e-\d+), Prob = (-\d+.\d+), TF norm = (\d+.\d+e\+\d+)', log[kk])
                candidates.write('{} {} {} {} {} {} {} {} {} {} {} {}\n'.format(jj, ii, recentinjections[jj], candidate.group(1), candidate.group(2), candidate.group(3), candidate.group(4), candidate.group(5), candidate.group(6), candidate.group(7), candidate.group(8), ifoval))
                kk += 1
                foundoutlier = 1

            if foundoutlier==0 and re.match('^system lalapps_TwoSpect failed', log[0]): candidates.write('{} {} {} {} {} {} {} {} {} {} {} {}\n'.format(jj, ii, recentinjections[jj], -1, -1, -1, -1, -1, -1, -1, ifoval))
            elif foundoutlier==0: candidates.write('{} {} {} {} {} {} {} {} {} {} {} {}\n'.format(jj, ii, recentinjections[jj], 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', ifoval))

            logfile.close()

    injectionfile.close()

if args.ULoutput: uloutput.close()
if args.candidateOutput: candidates.close()
        
