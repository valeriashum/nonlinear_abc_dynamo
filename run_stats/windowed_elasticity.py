#! /usr/bin/python2.7

import argparse

windowWidth = 10.0
a0=0.02
gradFrac = 2.0/3.0  # Limit used to choose when we accept the gradient as 'close' to the expected one

outfilename = "gradientClassifiedIsotropy.dat"
highHistFilename = "highHist.dat"
lowHistFilename = "lowHist.dat"
windowSemiWidth = windowWidth / 2.0

## Histogram stuff
h = 0.001 # bin size
## Domain limits
xiL = -1.0/6.0
xiR = 1.0/3.0
etaB = 0.0
etaT = 1.0/3.0
## Number of bins in the xi direction
Nxi = (int)((xiR - xiL)/h)
## Number of bins in the eta direction
Neta = (int)((etaT - etaB)/h)
## Total number of bins
N = Nxi * Neta

def bin(xi, eta):
    nxi  = (int)((xi  - xiL )  / h)
    neta = (int)((eta - etaB) / h)
    return nxi + neta*Nxi

def etaLocation(binNum):
    neta = (int)(binNum / Nxi)
    return neta*h + etaB

def xiLocation(binNum):
    neta = (int)(binNum / Nxi)
    nxi = binNum - neta*Nxi
    return nxi*h + xiL

def main():
    parser = argparse.ArgumentParser("Look at the timevar file, run a sliding window along it looking for regions where the gradient is 'close' to the predicted value, and then output the Rotta isotropy tensors based on this classification")
    parser.add_argument('firstIndex', metavar = 'first index', type=int, help = 'N s.t that first file we read is ../timevarN.dat')
    parser.add_argument('lastIndex', metavar = 'last index', type=int, help = 'N s.t that last file we read is ../timevarN.dat')
    parser.add_argument('KE', metavar = 'kinetic energy', type=float, help = 'the kinetic energy of the unsheared flow')
    parser.add_argument('--start', metavar = 'start time', type=float, help = 'time from which to start taking the fit')
    parser.add_argument('--end', metavar = 'end time', type=float, help = 'time from up to we take the fit')
    parser.add_argument('--hist', action = 'store_true', help = 'produce 2D histogram files')

    args = parser.parse_args()

    if args.start != None: tStart = args.start
    else: tStart = 0.0

    if args.end != None: tEnd = args.end
    else: tEnd = None

    outputHistograms = args.hist
    if outputHistograms:
        HHigh = [0 for x in range(Nxi) for y in range(Neta)]
        HLow  = [0 for x in range(Nxi) for y in range(Neta)]

    isotropicGradient = -(a0/2.0)*(4.0/15.0)*args.KE

    outfile = open(outfilename, "w")
    ## Write the header
    outfile.write("#t\tintSRxy\tRottaXi\tRottaEta\thighElast\n")

    ## So we can print out some stats at the end
    gradHighCount = 0
    gradLowCount = 0

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
            etaIndex = headerLine.index('RottaEtaSquared')
            xiIndex = headerLine.index('RottaXiCubed')
        except ValueError:
            print "File " + filename + " does not appear to contain the required data. Skipping."
        ## Read in all of the data points
        dataList = [map(float, [line.split()[i] for i in [tIndex, SIndex, etaIndex, xiIndex]]) for line in timevar]
        timevar.close()
        ## Throw away any before our start time
        dataList = [x for x in dataList if x[0] >= tStart]
        ## Throw away any after our end time
        if tEnd is not None: dataList = [x for x in dataList if x[0] <= tEnd]
        ## Set up the windowing pointers
        winStart = 0
        winCentre = winStart + 1
        while dataList[winCentre][0] < dataList[winStart][0] + windowSemiWidth:
            winCentre += 1
        winEnd = winCentre + 1
        while winEnd < len(dataList) and dataList[winEnd][0] < dataList[winCentre][0] + windowSemiWidth:
            winEnd += 1
        ## Loop over the points in the dataList, moving the window
        while winEnd < len(dataList):
            ## Process the data
            windowedGradient = (dataList[winEnd][1] - dataList[winStart][1])/(dataList[winEnd][0] - dataList[winStart][0])
            t = dataList[winCentre][0]
            intSRxy = dataList[winCentre][1]
            RottaEta = dataList[winCentre][2]**0.5
            RottaXi = (dataList[winCentre][3]/abs(dataList[winCentre][3]))*abs(dataList[winCentre][3])**(1.0/3.0)
            if windowedGradient / isotropicGradient >= gradFrac:
                highElast = 1
                gradHighCount += 1
                if outputHistograms:
                    HHigh[bin(RottaXi, RottaEta)] += 1
            else:
                highElast = 0
                gradLowCount += 1
                if outputHistograms:
                    HLow[bin(RottaXi, RottaEta)] += 1
            ## Write the datapoint to the file
            outfile.write("{}\t{}\t{}\t{}\t{}\n".format(t, intSRxy, RottaEta, RottaXi, highElast))
            ## Advance the window
            winCentre += 1
            while dataList[winStart][0] < dataList[winCentre][0] - windowSemiWidth:
                winStart += 1
            winEnd = winCentre + 1
            while winEnd < len(dataList) and dataList[winEnd][0] < dataList[winCentre][0] + windowSemiWidth:
                winEnd += 1
        timevar.close()

    outfile.close()

    if outputHistograms:
        ## Normalise
        CHigh = 0.0
        for x in HHigh: CHigh += x
        HHigh = [float(x)/CHigh for x in HHigh]
        CLow = 0.0
        for x in HLow: CLow += x
        HLow = [float(x)/CLow for x in HLow]
        ## Write the output to files
        outfileHigh = open(highHistFilename, "w")
        outfileHigh.write("#xi\teta\tdensity\n")
        for x in range(len(HHigh)):
            eta = etaLocation(x)
            xi  = xiLocation(x)
            outfileHigh.write("{}\t{}\t{}\n".format(xi, eta, HHigh[x]))
        outfileHigh.close()
        outfileLow = open(lowHistFilename, "w")
        outfileLow.write("#xi\teta\tdensity\n")
        for x in range(len(HLow)):
            eta = etaLocation(x)
            xi  = xiLocation(x)
            outfileLow.write("{}\t{}\t{}\n".format(xi, eta, HLow[x]))
        outfileLow.close()

    gradCountFrac = float(gradHighCount)/float(gradHighCount + gradLowCount)
    print "Gradient was accepted as 'close' to the predicted elasticity for a fraction {} of the data".format(gradCountFrac)

if __name__ == "__main__":
    main()
