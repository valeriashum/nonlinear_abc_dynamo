#! /usr/bin/python2.7

import argparse
import sys

def cov(xList, yList):
    if len(xList) != len(yList):
        raise ValueError('Lists must be the same length.')
    if len(xList) == 0:
        raise ValueError('Lists must be non-empty.')
    n = len(xList)
    xyBar = sum([x*y for x,y in zip(xList, yList)])/n
    xBaryBar = sum(xList)*sum(yList)/n**2
    return (xyBar - xBaryBar)

def var(xList):
    return cov(xList, xList)

def main():
    parser = argparse.ArgumentParser("Get statistics for a single timevar file")
    parser.add_argument('filename', metavar = 'file name', help = 'Name of the timevar file')
    parser.add_argument('--start', metavar = 'start time', type = float, help = 'We discard times earlier than this')
    parser.add_argument('--C', action = 'store_true', help = 'Output the gradient of the lease squares fit to the intCRxy curve.')
    parser.add_argument('--S', action = 'store_true', help = 'Output the gradient of the lease squares fit to the intSRxy curve.')

    args = parser.parse_args()

    filename = args.filename

    if args.start is not None: tStart = args.start
    else: tStart = 0.0

    try:
        timevar = open(filename, "r")
    except IOError:
        print "Unable to open " + filename + " for reading."
        sys.exit()

    line = timevar.readline().split()
    try:
        tIndex = line.index('t')
        CRxyIndex = line.index('intCRxy')
        SRxyIndex = line.index('intSRxy')
    except ValueError:
        print "Unable to locate all the necessary columns in " + filename
        sys.exit()
     
    tList = []
    CRxyList = []
    SRxyList = []
    for line in timevar:
        line = line.split()
        tList.append(float(line[tIndex]))
        CRxyList.append(float(line[CRxyIndex]))
        SRxyList.append(float(line[SRxyIndex]))

    CRxyList = [x for t,x in zip(tList, CRxyList) if t >= tStart]
    SRxyList = [x for t,x in zip(tList, SRxyList) if t >= tStart]
    tList    = [t for t in tList if t >= tStart]

    Cbeta1 = cov(CRxyList, tList) / var(tList)
    Sbeta1 = cov(SRxyList, tList) / var(tList)

    if args.C:
        print Cbeta1

    if args.S:
        print Sbeta1

if __name__ == '__main__':
    main()
