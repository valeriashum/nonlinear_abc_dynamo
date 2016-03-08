#! /usr/bin/python

'''
run_stats.py
This script is intended to allow quick extraction of statistics about
the flow from the timevar file produced by SNOOPY.
'''

import sys
import argparse

# Start of statistical function definitions
def mean(list):
    return sum(list) / len(list)

def variance(list):
    squaredList = [x*x for x in list]
    return sum(squaredList)/len(list) - (sum(list)/len(list))**2

def covariance(list1, list2):
    if len(list1) != len(list2):
        return 0
    else:
        return sum([x*y for x,y in zip(list1, list2)]) / len(list1) -  sum(list1)*sum(list2)/(len(list1)**2)
# End of statistical function definitions

parser = argparse.ArgumentParser(description='Extract statistics from the timevar file produced by SNOOPY.')

parser.add_argument('timevarFilename', metavar='timevar', type=str, help='timevar file output by SNOOPY')
parser.add_argument('--start', metavar='start time', type=float, help='time from which to start taking statistics')
parser.add_argument('--end', metavar='end time', type=float, help='time up to which we should take statistics')

args = parser.parse_args()

if args.start != None:
    t1 = args.start
else:
    t1 = 0.0

# Open the filename given

try:
    timevar = open(args.timevarFilename, 'r')
except IOError:
    print 'Unable to open ' + sys.argv[1] + ' for reading. Quitting...'
    sys.exit(2)

headerLine = timevar.readline().split()

if 't' in headerLine:
    tColumn = headerLine.index('t')
    print 'Found t in column ', tColumn
else:
    tColumn = -1

if 'intCRxy' in headerLine:
    intCRxyColumn = headerLine.index('intCRxy')
    print 'Found intCRxy in column ', intCRxyColumn
else:
    intCRxyColumn = -1

if not(intCRxyColumn < 0 and tColumn < 0):
    print '(Note that columns are zero-indexed)\n'

# Run through the file and process each line
if tColumn >=0: tList = []
if intCRxyColumn >=0: intCRxyList = []

for line in timevar:
    line = line.strip().split()
    if tColumn >=0: tList.append(float(line[tColumn]))
    if intCRxyColumn >=0: intCRxyList.append(float(line[intCRxyColumn]))

# Delete lines before the start time
if args.start != None and tColumn >= 0:
    for t in tList:
        if t >= args.start:
		    # Do all other columns first so the index of t doesn't change!
            if intCRxyColumn >= 0: intCRxyList = intCRxyList[tList.index(t):]
            tList = tList[tList.index(t):]
            break

# Delete lines after the end time
if args.end != None and tColumn >= 0:
    for t in tList:
        if t > args.end:
		    # Do all other columns first so the index of t doesn't change!
            if intCRxyColumn >= 0: intCRxyList = intCRxyList[:tList.index(t)]
            tList = tList[:tList.index(t)]
            break

# If possible, do a linear regression fit for the intCRxy time series
if tColumn >=0 and intCRxyColumn >=0:
    tVar = variance(tList)
    tMean = mean(tList)
    intCRxyMean = mean(intCRxyList)
    intCRxyVar = variance(intCRxyList)
    tintCRxyCov = covariance(tList, intCRxyList)
    w1 = tintCRxyCov / tVar
    w2 = intCRxyMean - tMean * tintCRxyCov/tVar
    # Pearson Correlation Coefficient
    PCC = tintCRxyCov / (tVar*intCRxyVar)**0.5
    # Signal to noise ratio measure
    rmsResiduals = (sum([(w1*t + w2 - intCRxy)**2 for t, intCRxy in zip(tList, intCRxyList)]) / len(tList))**0.5
    print 'Least squares fit gives\nintCRxy = ', w1, ' * t + ', w2
    print 'Pearson\'s correlation coefficient is ', PCC
    print 'Noise to signal ratio is ', rmsResiduals / (w1*(tList[-1]-tList[0])), '\n'
