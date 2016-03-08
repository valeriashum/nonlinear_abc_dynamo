#! /usr/bin/python2.7

## Program to calculate the gradient of a least squares fit to the
## integrals of sin(wt) R_xy and cos(wt) R_xy from the timevar files
## This program is capable of doing this for files from multiple
## runs, and averaging the results

import argparse
from math import sqrt

def main():
    parser = argparse.ArgumentParser("Perform statistic on timevar files.")
    parser.add_argument('firstIndex', metavar = 'first index', type=int, help = 'N s.t that first file we read is ../timevarN.dat')
    parser.add_argument('lastIndex', metavar = 'last index', type=int, help = 'N s.t that last file we read is ../timevarN.dat')
    parser.add_argument('--start', metavar = 'start time', type=float, help = 'time from which to start taking the fit')
    parser.add_argument('--end', metavar = 'end time', type=float, help = 'time from up to we take the fit')
    parser.add_argument('--verbose', dest = 'verbose', action='store_true', help = 'Print human readable output for each file.')

    args = parser.parse_args()

    if args.start != None: tStart = args.start
    else: tStart = 0.0

    if args.end != None: tEnd = args.end
    else: tEnd = None

    if args.verbose: verbose = True
    else: verbose = False

    ## Lists to hold the gradients of the best-fit lines
    Sbeta1List = []
    Cbeta1List = []
    KEMeanList = []

    for fileIndex in range(args.firstIndex, args.lastIndex+1):
        filename = "../timevar" + str(fileIndex) + ".dat"
        ## Open the file
        try:
            timevar = open(filename, 'r')
        except IOError:
            print "Unable to open " + filename + " for reading. Aborting reading files."
            break
        ## Find out which columns the data is in
        try:
            headerLine = timevar.readline().split()
            tIndex = headerLine.index('t')
            SIndex = headerLine.index('intSRxy')
            CIndex = headerLine.index('intCRxy')
            KEIndex = headerLine.index('ev')
        except ValueError:
            print "File " + filename + " does not appear to contain the required data. Skipping."
        ## Read in all of the data points
        dataList = [(float(line.split()[tIndex]), float(line.split()[SIndex]),
                     float(line.split()[CIndex]), float(line.split()[KEIndex])) for line in timevar]
        timevar.close()
        ## Throw away any before our start time
        dataList = [x for x in dataList if x[0] >= tStart]
        ## Throw away any after our end time
        if tEnd is not None: dataList = [x for x in dataList if x[0] <= tEnd]
        ## Define new lists for convenience
        tList = [x[0] for x in dataList]
        SList = [x[1] for x in dataList]
        CList = [x[2] for x in dataList]
        KEList = [x[3] for x in dataList]
        ## Compute the slope of the regression line for integral[cos(wt) * Rxy] and integral[cos(wt) * Rxy]
        Cbeta1 = beta1(tList, CList)
        Cbeta0 = beta0(tList, CList)
        Cbeta1List.append(Cbeta1)
        if verbose:
            print "Best fit line for integral[Rxy * cos(omega * t) dt]:"
            print str(Cbeta0) + " + " + str(Cbeta1) + " * t"
        Sbeta1 = beta1(tList, SList)
        Sbeta0 = beta0(tList, SList)
        Sbeta1List.append(Sbeta1)
        if verbose:
            print "Best fit line for integral[Rxy * sin(omega * t) dt]:"
            print str(Sbeta0) + " + " + str(Sbeta1) + " * t"
        KE = mean(KEList)
        KEMeanList.append(KE)
        if verbose:
            print "Mean kinetic energy:"
            print str(KE)

    ## Compute the mean and standard deviation of these regession line slopes
    print "From fitting intCRxy:"
    Cbeta1bar = mean(Cbeta1List)
    print "beta1 = " + str(Cbeta1bar)
    if len(Cbeta1List) > 1:
        Cbeta1error = sqrt(variance(Cbeta1List)/(len(Cbeta1List)-1))
        print "Standard error: " + str(Cbeta1error)
    print "From fitting intSRxy:"
    Sbeta1bar = mean(Sbeta1List)
    print "beta1 = " + str(Sbeta1bar)
    if len(Sbeta1List) > 1:
        Sbeta1error = sqrt(variance(Sbeta1List)/(len(Sbeta1List)-1))
        print "Standard error: " + str(Sbeta1error)
    print "Kinetic energy:"
    KEMean = mean(KEMeanList)
    print "mean = " + str(KEMean)
    if len(KEMeanList) > 1:
        KEMeanError = sqrt(variance(KEMeanList)/(len(KEMeanList)-1))
        print "Standard error: " + str(KEMeanError)

    ## Print in a format convenient for copying into my gradlist file
    if len(Cbeta1List) > 1 and len(Sbeta1List) > 1 and len(KEMeanList) > 1:
        print "Easy copy-and-paste line:"
        print str(Cbeta1bar) + "\t" + str(Cbeta1error) + "\t" + str(Sbeta1bar) + "\t" + str(Sbeta1error) + "\t" + str(KEMean) + "\t" + str(KEMeanError)
        


## This function isn't being used at present
## Since Snoopy attempts to do output as close as possible to
## a constant rate (e.g. every 0.1 time units) the smallest dt
## is never particularly small. Equally, it never gets very large.
## Therefore I don't think much error will be introduced by just using
## the raw points with no interpolation.
def smallestdt(dataList):
    ## Scan through the list and find the closest pair of times
    diffList = [a[0] - b[0] for a,b in zip(dataList[1:], dataList[0:-1])]
    #print diffList
    return reduce(min, diffList)

def mean(x):
    return (sum(x)/len(x))

## Compute the sample variance of x
## i.e. with the 1/(n-1) factor
def variance(x):
    n = len(x)
    xbar = mean(x)
    return ((1/float(n-1)) * sum([(xi - xbar)**2 for xi in x]))

## Compute the sample covariance of x and y
## i.e. with the 1/(n-1) factor
def covariance(x, y):
    if len(x) == len(y):
        n = len(x)
    ## FIXME - raise an exception if lengths differ
    return ((sum([a*b for a,b in zip(x,y)]) - n*mean(x)*mean(y))/(n-1))

## Compute the slope of the least squares line for the data points (x,y)
## This assumes the existence of a beta0 that is also free to vary
def beta1(x, y):
    return (covariance(x,y)/variance(x))

## Compute the intercept of the least squares line for the data points (x,y)
## This assumes the existence of a beta1 that is also free to vary
def beta0(x, y):
    return (mean(y) - beta1(x, y)*mean(x))

if __name__ == "__main__":
    main()
